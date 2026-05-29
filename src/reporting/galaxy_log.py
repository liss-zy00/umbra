"""
Per-galaxy provenance log block.

Each galaxy processed by the pipeline (validation or discovery) emits a compact,
aligned block that cross-references every data source consulted before the
physics engine runs. Centralising the format here keeps discovery vs validation
styling consistent and makes the layout trivially testable.

Output convention
-----------------
Unicode light horizontal rule (U+2500) for the title line (no ``[phase]`` suffix).
Each provider row uses ``│`` only after the source tag; fields inside the row use
narrow no-break space + middle dot (U+202F U+00B7 U+202F) between ``label: value`` pairs.
When emitted from the pipeline with ``--logs``, every line is prefixed with four spaces
so the block nests under the macro-phase banner (same indent as phase body lines).

Example inner layout (each stored line also has a four-space phase indent prefix)::

    ── DDO154 ───────────────────────────────────────────
     SPARC  │ dist: 4.04 Mpc · shells: 12 · file: DDO154_rotmod.dat
     MORPH  │ type: Im · T: n/d · incl: 64.0° · vflat: 47.0 km/s
       EDD  │ pgc: 12345 · tier: 2 · src: edd_cf4_DM · D: 4.04 Mpc · σ_D: 0.050 Mpc
       NED  │ status: … · ra: … · dec: … · z: … · obj: …
    PLANCK  │ Q: … · U: … · scale: … kpc/px

When the same galaxy already printed a full discovery block in the current run,
validation logging may print a shortened stub per galaxy during the grid. Phase 2
(``--logs``) prints a scenario legend (two intro lines plus one ``→ Scenario k`` line per
``Δx`` / ``f_ion`` pair), then stderr galaxy progress, then a deep-indented name grid when the
loop finishes. Phase 3/4 print short stdout gloss blocks when ``--logs`` is enabled.
Discovery still emits the five-line provenance blocks per galaxy when ``--logs`` is set.

Missing values are rendered as ``n/d`` rather than ``None`` or ``nan`` so
downstream readers (humans grepping the terminal, log scrapers) get a stable
sentinel instead of a Python ``repr``.
"""

from __future__ import annotations

import math
import shutil
from typing import Any, Literal, Sequence

import numpy as np

_RULE_CHAR = "\u2500"  # BOX DRAWINGS LIGHT HORIZONTAL
_MIN_RULE_WIDTH = 88
_MAX_RULE_WIDTH = 100
_SOURCE_TAG_WIDTH = 6
# Same as pipeline phase body indent so per-galaxy ``--logs`` blocks read as nested content.
_PROVENANCE_BLOCK_INDENT = "    "
_VALIDATION_LOG_DEEP_INDENT = "        "  # Phase 2 ``--logs`` legend + galaxy table (nested under banner).
# Between label:value pairs inside a provider row (after the tag ``│``).
_FIELD_SEP = "\u202f\u00b7\u202f"

_MISSING: str = "n/d"

PipelineLogPhase = Literal["discovery", "validation"]

VALIDATION_PROVENANCE_STUB_LINE = (
    "→ Provenance verified in discovery. Running the physical scenario grid..."
)


def _terminal_rule_width() -> int:
    try:
        cols = shutil.get_terminal_size((_MIN_RULE_WIDTH, 24)).columns
    except OSError:
        return _MIN_RULE_WIDTH
    return max(_MIN_RULE_WIDTH, min(_MAX_RULE_WIDTH, cols))


def format_rule_banner(center_text: str, bracket: str, width: int | None = None) -> str:
    """``── {center_text} ───… [{bracket}]`` padded to terminal width (or ``width`` when given)."""
    target = width if width is not None else _terminal_rule_width()
    suffix = f" [{bracket}]"
    prefix = f"{_RULE_CHAR}{_RULE_CHAR} {center_text} "
    min_total = len(prefix) + 2 + len(suffix)
    total = max(target, min_total)
    filler_len = total - len(prefix) - len(suffix)
    filler_len = max(2, filler_len)
    return prefix + (_RULE_CHAR * filler_len) + suffix


def format_rule_fill(width: int | None = None) -> str:
    """Solid horizontal rule (same width as :func:`format_rule_banner` defaults)."""
    w = width if width is not None else _terminal_rule_width()
    return _RULE_CHAR * w


def format_rule_banner_plain(center_text: str, width: int | None = None) -> str:
    """Light rule with ``center_text`` only (no trailing `` [bracket]``); width matches ``width``."""
    target = width if width is not None else _terminal_rule_width()
    prefix = f"{_RULE_CHAR}{_RULE_CHAR} {center_text} "
    min_total = len(prefix) + 2
    total = max(target, min_total)
    filler_len = total - len(prefix)
    filler_len = max(2, filler_len)
    return prefix + (_RULE_CHAR * filler_len)


def _tag_cell(source: str) -> str:
    return source.rjust(_SOURCE_TAG_WIDTH)


def _fmt_float(value: Any, fmt: str) -> str:
    """Format ``value`` with ``fmt`` if it is a finite number, else ``n/d``."""
    if value is None:
        return _MISSING
    try:
        f = float(value)
    except (TypeError, ValueError):
        return _MISSING
    if not np.isfinite(f):
        return _MISSING
    return format(f, fmt)


def _fmt_int(value: Any) -> str:
    if value is None:
        return _MISSING
    try:
        f = float(value)
    except (TypeError, ValueError):
        return _MISSING
    if not np.isfinite(f):
        return _MISSING
    return str(int(f))


def _fmt_str(value: Any) -> str:
    if value is None:
        return _MISSING
    text = str(value).strip()
    return text if text else _MISSING


def _pipe_row(source: str, segments: list[tuple[str, str]]) -> str:
    body = _FIELD_SEP.join(f"{label}: {value}" for label, value in segments)
    return f"{_tag_cell(source)}  │ {body}"


def _edd_row(edd: dict[str, Any] | None) -> str:
    """Extragalactic Distance Database join / tier / adopted luminosity distance (see ``edd_provider``)."""
    e = edd if edd is not None else {}
    if not e.get("edd_tables_loaded", False):
        segs: list[tuple[str, str]] = [("status", "catalogues not loaded")]
        if e.get("pgc") is not None:
            segs.insert(0, ("pgc", _fmt_int(e.get("pgc"))))
        return _pipe_row("EDD", segs)
    if e.get("pgc") is None:
        return _pipe_row(
            "EDD",
            [
                ("pgc", _MISSING),
                ("status", "no PGC map entry"),
            ],
        )
    pgc_s = _fmt_int(e.get("pgc"))
    if e.get("applied"):
        return _pipe_row(
            "EDD",
            [
                ("pgc", pgc_s),
                ("tier", _fmt_int(e.get("tier"))),
                ("src", _fmt_str(e.get("source_table"))),
                ("D", f"{_fmt_float(e.get('distance_mpc'), '.2f')} Mpc"),
                ("σ_D", f"{_fmt_float(e.get('distance_uncertainty_mpc'), '.3f')} Mpc"),
            ],
        )
    return _pipe_row(
        "EDD",
        [
            ("pgc", pgc_s),
            ("status", "no catalogue distance row"),
            ("D", f"{_fmt_float(e.get('distance_mpc'), '.2f')} Mpc"),
        ],
    )


def format_galaxy_provenance_block(
    galaxy_name: str,
    *,
    sparc: dict[str, Any] | None = None,
    morph: dict[str, Any] | None = None,
    edd: dict[str, Any] | None = None,
    ned: dict[str, Any] | None = None,
    planck: dict[str, Any] | None = None,
    pipeline_phase: PipelineLogPhase = "discovery",
    rule_width: int | None = None,
) -> str:
    """
    Build the multi-line provenance block as a single string.

    ``pipeline_phase`` is accepted for API stability; the rule line no longer prints a
    phase bracket. The five provider rows (SPARC, MORPH, EDD, NED, PLANCK) are identical
    for both phases (callers choose whether to print them during validation via
    :func:`emit_validation_provenance_stub`; the cohort banner is emitted from the runner
    when ``--logs`` is enabled).
    """
    _ = pipeline_phase
    tw = rule_width if rule_width is not None else _terminal_rule_width()
    inner_w = max(48, tw - len(_PROVENANCE_BLOCK_INDENT))
    header = format_rule_banner_plain(galaxy_name, inner_w)

    sparc = sparc or {}
    morph = morph or {}
    edd = edd if edd is not None else {}
    ned = ned or {}
    planck = planck or {}

    scale_val = _fmt_float(planck.get("pixel_scale_kpc"), ".2f")
    scale_cell = f"{scale_val} kpc/px" if scale_val != _MISSING else _MISSING

    lines = [
        header,
        _pipe_row(
            "SPARC",
            [
                ("dist", f"{_fmt_float(sparc.get('distance_mpc'), '.2f')} Mpc"),
                ("shells", _fmt_int(sparc.get("n_shells"))),
                ("file", _fmt_str(sparc.get("source_file"))),
            ],
        ),
        _pipe_row(
            "MORPH",
            [
                ("type", _fmt_str(morph.get("morphology_class"))),
                ("T", _fmt_str(morph.get("hubble_type"))),
                ("incl", f"{_fmt_float(morph.get('inclination_deg'), '.1f')}°"),
                ("vflat", f"{_fmt_float(morph.get('vflat_kms'), '.1f')} km/s"),
            ],
        ),
        _edd_row(edd),
        _pipe_row(
            "NED",
            [
                ("status", _fmt_str(ned.get("status"))),
                ("ra", f"{_fmt_float(ned.get('right_ascension_degrees'), '.3f')}°"),
                ("dec", f"{_fmt_float(ned.get('declination_degrees'), '.3f')}°"),
                ("z", _fmt_float(ned.get("redshift"), ".5f")),
                ("obj", _fmt_str(ned.get("object_type"))),
            ],
        ),
        _pipe_row(
            "PLANCK",
            [
                ("Q", _fmt_float(planck.get("q_center"), ".2e")),
                ("U", _fmt_float(planck.get("u_center"), ".2e")),
                ("scale", scale_cell),
            ],
        ),
    ]
    return "\n".join(_PROVENANCE_BLOCK_INDENT + ln for ln in lines)


def _scenario_matrix_line_label(delta_x_kpc: float, f_ion: float) -> str:
    """One matrix cell as explicit Δx and f_ion strings (values come from scenario pairs)."""
    return f"Δx={float(delta_x_kpc):.3f} kpc f_ion={float(f_ion):.3f}"


def _scenario_legend_compact_knobs(delta_x_kpc: float, f_ion: float) -> str:
    """Δx / f_ion segment for scenario lines (same separators as provider rows)."""
    dx = float(delta_x_kpc)
    fi = float(f_ion)
    return f"Δx={dx:.3f} kpc{_FIELD_SEP}f_ion={fi:.3f}"


def format_validation_phase_scenario_legend(
    scenario_pairs: Sequence[tuple[float, float]],
    *,
    n_galaxies: int,
    line_indent: str | None = None,
) -> str:
    """
    One line per grid scenario, printed once at Phase 2 start when ``--logs`` tracks execution.

    Order matches :func:`src.pipeline.per_galaxy.scenario_grid_runner.run_scenario_grid` evaluation.
    """
    ind = _PROVENANCE_BLOCK_INDENT if line_indent is None else line_indent
    pairs = [(float(g), float(f_)) for g, f_ in scenario_pairs]
    n = len(pairs)
    if n == 0:
        return f"{ind}(no grid scenario pairs)\n\n"
    intro_a = f"{ind}GRADIENT_SCALES_KPC × IONIZATION_FRACTIONS (nested);"
    intro_b = f"{ind}per file: scenarios 1→{n}, then next galaxy."
    _ = n_galaxies
    lines = [
        f"{ind}→ Scenario {i + 1}{_FIELD_SEP}{_scenario_legend_compact_knobs(g, f_)}"
        for i, (g, f_) in enumerate(pairs)
    ]
    return intro_a + "\n" + intro_b + "\n\n" + "\n".join(lines) + "\n\n"


def format_phase3_graphs_verbose_block(
    *,
    n_targets: int,
    graphs_rel: str,
    line_indent: str | None = None,
) -> str:
    """Stdout helper for ``--logs``: what Phase 3 (graphs) reads and writes before matplotlib work."""
    ind = _PROVENANCE_BLOCK_INDENT if line_indent is None else line_indent
    body = [
        f"{ind}galaxy_analysis → per-galaxy PNGs + {graphs_rel}/population cohort plots ({n_targets} names).",
        f"{ind}Stderr Processing line = matplotlib export queue.",
    ]
    return "\n".join(body) + "\n"


def format_phase4_analysis_verbose_block(
    *,
    audit_rel: str,
    n_summary_scan: int,
    line_indent: str | None = None,
) -> str:
    """Stdout helper for ``--logs``: read-only audit pass before Markdown is written."""
    ind = _PROVENANCE_BLOCK_INDENT if line_indent is None else line_indent
    body = [
        f"{ind}read {n_summary_scan} *_summary.json (evidence|reliable|not_reliable) → bias + covariate stats.",
        f"{ind}Write {audit_rel}.",
    ]
    return "\n".join(body) + "\n"


def format_scenario_matrix_one_per_line(
    scenario_pairs: Sequence[tuple[float, float]],
    *,
    line_prefix: str = "  → ",
) -> list[str]:
    """Emit one scenario combination per output line (``--logs`` validation batch)."""
    if not scenario_pairs:
        return [f"{line_prefix}(no scenario pairs)"]
    return [f"{line_prefix}{_scenario_matrix_line_label(float(g), float(f_))}" for g, f_ in scenario_pairs]


def format_validation_grid_batch_block(
    *,
    n_cohort_galaxies: int,
    n_discovery_reliable: int | None,
    n_scenario_combinations: int,
    scenario_pairs: Sequence[tuple[float, float]],
    rule_width: int | None = None,
) -> str:
    """One banner plus compact cohort / matrix / total lines for the validation grid (``--logs``)."""
    title = "Validation (multi-scenario matrix)"
    banner = format_rule_banner(title, "validation", rule_width)
    total_eval = n_cohort_galaxies * n_scenario_combinations
    if n_discovery_reliable is not None:
        rel_seg = f"{n_discovery_reliable} discovery-reliable in census"
    else:
        rel_seg = "discovery-reliable count n/a"

    matrix_lines = format_scenario_matrix_one_per_line(scenario_pairs)
    lines: list[str] = [
        banner,
        f"→ Cohort (pre-filter passed): {n_cohort_galaxies} galaxies ({rel_seg}).",
        f"→ Scenario matrix: {n_scenario_combinations} combinations/galaxy (Δx_kpc × f_ion).",
        *matrix_lines,
        f"→ Total evaluations (upper bound): {total_eval} (cohort × combinations).",
    ]
    return "\n".join(lines)


def emit_validation_grid_batch_header(
    *,
    n_cohort_galaxies: int,
    n_discovery_reliable: int | None,
    n_scenario_combinations: int,
    scenario_pairs: Sequence[tuple[float, float]],
    enabled: bool = True,
    rule_width: int | None = None,
) -> None:
    """Print the validation cohort banner once before the first grid galaxy starts."""
    if not enabled or n_cohort_galaxies <= 0:
        return
    block = format_validation_grid_batch_block(
        n_cohort_galaxies=n_cohort_galaxies,
        n_discovery_reliable=n_discovery_reliable,
        n_scenario_combinations=n_scenario_combinations,
        scenario_pairs=scenario_pairs,
        rule_width=rule_width,
    )
    print(f"\n{block}", flush=True)


def format_completed_galaxy_grid(
    names: Sequence[str],
    *,
    rule_width: int | None = None,
    target_rows: int = 7,
    gap: int = 2,
    line_indent: str | None = None,
) -> str:
    """
    Lay out galaxy names in fixed-width columns, preferring about ``target_rows`` rows.

    Column count starts at ``ceil(n / target_rows)`` and decreases until rows fit
    ``rule_width`` minus ``line_indent`` (defaults to :data:`_VALIDATION_LOG_DEEP_INDENT`, i.e. two
    levels under the phase banner). Output ends with a blank line (caller should use
    ``print(..., end="")`` when printing this string as-is).
    """
    ind = _VALIDATION_LOG_DEEP_INDENT if line_indent is None else line_indent
    names_list = [str(n).strip() for n in names if str(n).strip()]
    if not names_list:
        return f"{ind}(no galaxies in completion list)\n\n"
    w = rule_width if rule_width is not None else _terminal_rule_width()
    usable = max(24, w - len(ind) - 1)
    max_len = max(len(n) for n in names_list)
    sep = " " * gap
    n = len(names_list)
    target_rows_eff = max(1, int(target_rows))
    chosen_cols = max(1, math.ceil(n / target_rows_eff))
    while chosen_cols > 1:
        row_width = chosen_cols * max_len + (chosen_cols - 1) * len(sep)
        if row_width <= usable:
            break
        chosen_cols -= 1
    n_rows = math.ceil(n / chosen_cols)
    body_lines: list[str] = []
    for i in range(0, n, chosen_cols):
        row = names_list[i : i + chosen_cols]
        padded = [name.ljust(max_len) for name in row]
        body_lines.append(ind + sep.join(padded).rstrip())
    head = f"{ind}Galaxies processed ({n}); rows={n_rows}; columns={chosen_cols}"
    return head + "\n" + "\n".join(body_lines) + "\n\n"


def emit_validation_grid_progress_rewrite(
    completed_names: Sequence[str],
    *,
    total: int,
    enabled: bool = True,
    rule_width: int | None = None,
) -> None:
    """Legacy hook: inline ``\\r`` counter was removed from the ``--logs`` UX (no-op)."""
    _ = completed_names
    _ = total
    _ = enabled
    _ = rule_width
    return


def emit_validation_grid_progress_finalize(
    *,
    enabled: bool = True,
    completed_names: Sequence[str] | None = None,
    rule_width: int | None = None,
) -> None:
    """Print the full galaxy grid (indented, trailing blank line); no ``\\r`` clear step."""
    if not enabled:
        return
    if completed_names is not None and len(completed_names) > 0:
        text = format_completed_galaxy_grid(completed_names, rule_width=rule_width)
        print(text, end="", flush=True)


def format_validation_provenance_stub(galaxy_name: str, *, rule_width: int | None = None) -> str:
    """Two-line validation stub (plain rule header + stub line), indented like full provenance."""
    tw = rule_width if rule_width is not None else _terminal_rule_width()
    inner_w = max(48, tw - len(_PROVENANCE_BLOCK_INDENT))
    head = format_rule_banner_plain(galaxy_name, inner_w)
    return (
        _PROVENANCE_BLOCK_INDENT
        + head
        + "\n"
        + _PROVENANCE_BLOCK_INDENT
        + VALIDATION_PROVENANCE_STUB_LINE
    )


def emit_galaxy_provenance_block(
    galaxy_name: str,
    *,
    sparc: dict[str, Any] | None = None,
    morph: dict[str, Any] | None = None,
    edd: dict[str, Any] | None = None,
    ned: dict[str, Any] | None = None,
    planck: dict[str, Any] | None = None,
    pipeline_phase: PipelineLogPhase = "discovery",
    enabled: bool = True,
    rule_width: int | None = None,
) -> None:
    """Print the provenance block when ``enabled`` is True (default); no-op otherwise."""
    if not enabled:
        return
    block = format_galaxy_provenance_block(
        galaxy_name,
        sparc=sparc,
        morph=morph,
        edd=edd,
        ned=ned,
        planck=planck,
        pipeline_phase=pipeline_phase,
        rule_width=rule_width,
    )
    print(f"\n{block}", flush=True)


def emit_validation_provenance_stub(
    galaxy_name: str,
    *,
    enabled: bool = True,
    rule_width: int | None = None,
) -> None:
    """
    Print the shortened validation banner when discovery already emitted the full
    provenance block for this galaxy in the same process.

    Assumes provider inputs are unchanged since the discovery pass; use the full
    :func:`emit_galaxy_provenance_block` path when that contract does not hold.
    """
    if not enabled:
        return
    print("\n" + format_validation_provenance_stub(galaxy_name, rule_width=rule_width), flush=True)
