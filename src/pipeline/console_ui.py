"""
Terminal UI abstraction for the batch pipeline.

Rich prints phase banners on stdout with the same Unicode rule style as
:mod:`src.reporting.galaxy_log`. Long-running loops use a **fixed-width**
``Processing ━━━…`` line on stderr (no embedded phase title). Pass ``--logs`` for
per-galaxy provenance from :mod:`src.reporting.galaxy_log`.

``PlainPipelineConsole`` is retained for tests that inject ``console=`` into
:func:`src.pipeline.runner.run_pipeline`.
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
import time
from abc import ABC, abstractmethod
from collections.abc import Iterable, Iterator, Sequence
from contextlib import AbstractContextManager, contextmanager, nullcontext
from typing import TypeVar

from src.reporting.galaxy_log import format_rule_banner

T = TypeVar("T")

# Session + phase titles (stdout); progress bar on stderr does not repeat these strings.
PIPELINE_SESSION_RUN_CENTER = "E-G model pipeline"
PIPELINE_PHASE_1_DISCOVERY_TITLE = "Phase 1 — Discovery"
PIPELINE_PHASE_3_ANALYSIS_TITLE = "Phase 3 — Analysis"
PIPELINE_PHASE_4_GRAPHS_TITLE = "Phase 4 — Graphs"
PIPELINE_PHASE_5_VIZ_TITLE = "Phase 5 — 3D"

# Back-compat aliases (older modules / tests may still import these names).
PIPELINE_PHASE_3_POPULATION_GRAPHS_TITLE = PIPELINE_PHASE_4_GRAPHS_TITLE
PIPELINE_PHASE_4_CONSISTENCY_ANALYSIS_TITLE = PIPELINE_PHASE_3_ANALYSIS_TITLE
PIPELINE_PHASE_3_PUBLICATION_TITLE = PIPELINE_PHASE_4_GRAPHS_TITLE
PIPELINE_PHASE_4_AUDIT_TITLE = PIPELINE_PHASE_3_ANALYSIS_TITLE


def pipeline_phase_2_validation_title(n_scenarios: int) -> str:
    """Short phase-2 title; ``n_scenarios`` is accepted for call-site parity (counts live in details)."""
    _ = n_scenarios
    return "Phase 2 — Validation"


_INDENT_L1 = " "
_INDENT_L2 = "    "

# stderr progress line aligns with phase body lines (same indent as phase detail text).
_PROGRESS_BODY_INDENT = _INDENT_L2

# stderr progress bar (heavy fill, light empty, optional tick at partial cell).
_BAR_FILL = "\u2501"
_BAR_EMPTY = "\u2504"
_BAR_TICK = "\u2578"
# Erase full row (CSI EL) then redraw; avoids ghost lines when stdout interleaves or the bar wraps.
_PROGRESS_LINE_CLEAR = "\r\x1b[2K"

# Session explainer (stdout): keep each logical sentence on one physical line where possible.
_SESSION_OVERVIEW_LINES = (
    "Ingests SPARC rotation curves, NED metadata, and Planck Stokes maps.",
    "Runs the physics engine and writes audited artefacts under results/.",
    "This target set drives step_discovery, step_validation, galaxy_analysis, and population graphs in one invocation.",
)


def _stderr_terminal_columns() -> int:
    try:
        return max(48, shutil.get_terminal_size((88, 24)).columns)
    except OSError:
        return 88


def pipeline_terminal_columns() -> int:
    """TTY width for pipeline rules and stderr progress (uncapped; matches widest bar)."""
    return _stderr_terminal_columns()


def pipeline_session_start_block(*, rule_width: int | None = None) -> str:
    """Top-level session banner plus overview lines (same indent as phase details)."""
    w = rule_width if rule_width is not None else pipeline_terminal_columns()
    banner = "\n" + format_rule_banner(PIPELINE_SESSION_RUN_CENTER, "session", w) + "\n"
    body = "\n".join(f"{_INDENT_L2}{line}" for line in _SESSION_OVERVIEW_LINES)
    return f"{banner}{body}"


def format_processing_progress_line(
    n: int,
    total: int,
    *,
    elapsed_s: float,
    terminal_cols: int | None = None,
    line_indent: str | None = None,
) -> str:
    """
    Single-line stderr progress: indented ``Processing`` + Unicode bar + ``pct (n/total)``.

    ``line_indent`` defaults to the same width as phase detail lines so stderr aligns with
    stdout phase blocks. Total width matches :func:`pipeline_terminal_columns` so rules and the
    bar share the same right edge as indented phase banners (one-space rule prefix).
    """
    cols = terminal_cols if terminal_cols is not None else _stderr_terminal_columns()
    indent = _PROGRESS_BODY_INDENT if line_indent is None else line_indent
    # One column shy of full width avoids terminal soft-wrap that splits one logical line into two rows.
    max_len = max(48, cols - 1)
    total_eff = max(int(total), 1)
    n_clamped = min(max(int(n), 0), total_eff)
    pct = 100.0 * n_clamped / total_eff

    prefix = indent + "Processing "
    rhs_core = f" {pct:3.0f}% ({n_clamped}/{total_eff})"
    eta_s = 0.0
    if n_clamped > 0 and n_clamped < total_eff:
        eta_s = max(0.0, (elapsed_s / n_clamped) * (total_eff - n_clamped))
    tail = f" [{int(elapsed_s)}s<{int(eta_s)}s]"

    # Reserve space for tail only if the full line still fits.
    trial_rhs = rhs_core + tail
    bar_budget = max_len - len(prefix) - len(trial_rhs)
    if bar_budget < 8:
        trial_rhs = rhs_core
        bar_budget = max_len - len(prefix) - len(trial_rhs)
    bar_len = max(8, bar_budget)

    frac = bar_len * n_clamped / total_eff
    filled = int(math.floor(frac))
    partial = frac - filled
    if partial > 1e-9 and filled < bar_len:
        mid = _BAR_FILL * filled + _BAR_TICK + _BAR_EMPTY * (bar_len - filled - 1)
    else:
        mid = _BAR_FILL * filled + _BAR_EMPTY * (bar_len - filled)

    line = prefix + mid + trial_rhs
    if len(line) < max_len:
        line = line + " " * (max_len - len(line))
    return line[:max_len]


def _banner_width_for_indent(indent: str, explicit: int | None) -> int:
    term = explicit if explicit is not None else pipeline_terminal_columns()
    return max(48, term - len(indent))


def format_pipeline_summary_section(
    heading: str,
    bracket: str,
    detail_lines: Sequence[str],
    *,
    rule_width: int | None = None,
    use_detail_arrow: bool = False,
) -> str:
    """Indented rule banner plus detail lines; optional ``→`` prefix for closing-tree / docs only."""
    wd = _banner_width_for_indent(_INDENT_L1, rule_width)
    banner_line = _INDENT_L1 + format_rule_banner(heading, bracket, wd)
    if use_detail_arrow:
        bodies = [f"{_INDENT_L2}\u2192 {t}" for t in detail_lines]
    else:
        bodies = [f"{_INDENT_L2}{t}" for t in detail_lines]
    return "\n".join([banner_line, *bodies])


def _split_detail_parts(detail: str) -> list[str]:
    parts = [p.strip() for p in detail.split(" · ") if p.strip()]
    if parts:
        return parts
    stripped = detail.strip()
    return [stripped] if stripped else []


def _phase_start_lines(title: str, *, detail: str, rule_width: int | None = None) -> str:
    """Phase open: rule banner only; stderr progress is the first runtime feedback."""
    _ = detail
    wd = _banner_width_for_indent(_INDENT_L1, rule_width)
    banner_line = _INDENT_L1 + format_rule_banner(title, "phase", wd)
    return f"\n{banner_line}"


def _phase_end_lines(title: str, *, detail: str, rule_width: int | None = None) -> str:
    """Completion bullet lines only (no leading newline: ``phase_start`` already ends on a fresh stdout row)."""
    _ = title
    _ = rule_width
    if detail:
        parts = _split_detail_parts(detail)
        if not parts:
            body_lines = [f"{_INDENT_L2}Phase finished"]
        else:
            body_lines = [f"{_INDENT_L2}{p}" for p in parts]
    else:
        body_lines = [f"{_INDENT_L2}Phase finished"]
    return "\n".join(body_lines)


def _iter_with_fixed_processing_bar(
    iterable: Iterable[T],
    *,
    total: int | None,
) -> Iterator[T]:
    """Drive stderr progress with CSI line erase between redraws; newline when done."""
    if total is None:
        seq = list(iterable)
        total_i = len(seq)
    else:
        seq = iterable
        total_i = int(total)
    total_i = max(total_i, 1)
    t0 = time.perf_counter()
    n = 0
    sys.stdout.flush()
    for item in seq:
        n += 1
        elapsed = time.perf_counter() - t0
        line = format_processing_progress_line(n, total_i, elapsed_s=elapsed)
        sys.stderr.write(_PROGRESS_LINE_CLEAR + line)
        sys.stderr.flush()
        yield item
    sys.stderr.write("\n")
    sys.stderr.flush()


def build_pipeline_console(_args: argparse.Namespace) -> PipelineConsole:
    return RichPipelineConsole()


class PipelineConsole(ABC):
    """Console strategy for :func:`src.pipeline.runner.run_pipeline`."""

    @abstractmethod
    def emit(self, message: str) -> None:
        """General informational line."""

    @abstractmethod
    def phase_start(self, title: str, *, detail: str = "") -> None:
        """Macro-phase entry: indented rule banner only (``detail`` reserved for API stability)."""

    @abstractmethod
    def phase_end(self, title: str, *, detail: str = "") -> None:
        """Macro-phase exit: overview lines on stdout (no second banner)."""

    @abstractmethod
    def iter_progress(
        self,
        iterable: Iterable[T],
        *,
        desc: str,
        total: int | None = None,
        leave: bool = False,
        enabled: bool,
    ) -> Iterator[T]:
        """Wrap ``iterable`` with a fixed stderr progress line when ``enabled`` (``desc`` ignored)."""

    @abstractmethod
    def print_final_summary(self, lines: list[str]) -> None:
        """Closing block."""

    @abstractmethod
    def planck_activity(self, message: str) -> AbstractContextManager[None]:
        """Context manager around slow Planck / HEALPix setup (Rich status when available)."""


class PlainPipelineConsole(PipelineConsole):
    def emit(self, message: str) -> None:
        print(message, flush=True)

    def phase_start(self, title: str, *, detail: str = "") -> None:
        print(_phase_start_lines(title, detail=detail), flush=True)

    def phase_end(self, title: str, *, detail: str = "") -> None:
        print(_phase_end_lines(title, detail=detail), flush=True)

    def iter_progress(
        self,
        iterable: Iterable[T],
        *,
        desc: str,
        total: int | None = None,
        leave: bool = False,
        enabled: bool,
    ) -> Iterator[T]:
        if not enabled:
            yield from iterable
            return
        _ = desc
        _ = leave
        yield from _iter_with_fixed_processing_bar(iterable, total=total)

    def print_final_summary(self, lines: list[str]) -> None:
        print("\n".join(lines), flush=True)

    def planck_activity(self, message: str) -> AbstractContextManager[None]:
        return nullcontext()


class RichPipelineConsole(PipelineConsole):
    def __init__(self) -> None:
        from rich.console import Console

        self._console = Console(stderr=False)

    def emit(self, message: str) -> None:
        self._console.print(message)

    def phase_start(self, title: str, *, detail: str = "") -> None:
        self._console.print(_phase_start_lines(title, detail=detail))

    def phase_end(self, title: str, *, detail: str = "") -> None:
        self._console.print(_phase_end_lines(title, detail=detail))

    def iter_progress(
        self,
        iterable: Iterable[T],
        *,
        desc: str,
        total: int | None = None,
        leave: bool = False,
        enabled: bool,
    ) -> Iterator[T]:
        if not enabled:
            yield from iterable
            return
        _ = desc
        _ = leave
        yield from _iter_with_fixed_processing_bar(iterable, total=total)

    def print_final_summary(self, lines: list[str]) -> None:
        # Full-width Rich print avoids narrow Panel reflow on long paths.
        self._console.print("\n".join(lines))

    @contextmanager
    def planck_activity(self, message: str) -> Iterator[None]:
        from rich.status import Status

        with self._console.status(f"[magenta]{message}[/]", spinner="dots"):
            yield
