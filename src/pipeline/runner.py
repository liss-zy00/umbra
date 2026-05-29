"""
Pipeline orchestration entry point.

Wires the data providers (:mod:`src.providers`), the per-galaxy processing chain
(:mod:`src.pipeline.per_galaxy`), the cross-galaxy aggregation layer
(:mod:`src.pipeline.aggregation`) and the reporting writers (:mod:`src.reporting`).

**Unified flow** (``python src/main.py``):

1. **Phase 1 — Discovery** — one engine pass per selected galaxy at the physical Planck
   pixel scale with ``f_ion = 1.0``; writes census CSVs under ``results/step_discovery/``.
2. **Automatic pre-filter** — galaxies whose Discovery row carries a blocking audit
   (see :func:`src.pipeline.per_galaxy.discovery_processor.discovery_audits_clean`) are
   excluded from the grid; all others enter **Phase 2 — Validation**.
3. **Phase 3 — Analysis** — population consistency markdown under ``results/step_validation/``
   (runs after validation artefacts exist and **before** publication figures).
4. **Phase 4 — Graphs** — :mod:`src.reporting.publication_bundle` (per-galaxy and cohort figures).
5. **Phase 5 — 3D** — supreme-cohort PyVista PNG + interactive HTML under ``results/graphs/viz_3d/``.

Phases 1–3 are :func:`src.pipeline.workflows.run_physics_workflow`; phases 4–5 are
:func:`src.pipeline.workflows.run_visualization_workflow`. Use ``--compute-only`` or ``--viz-only``
to run a single branch.

``results/galaxy_analysis/{validation_scientific_evidence|validation_reliable|validation_not_reliable}/<NAME>/``
is written during **Phase 2 (validation)**.

Galaxy selection is controlled solely by ``--galaxies`` (``all``, ``default``, or explicit names),
except in ``--viz-only`` mode where cohorts are read from existing CSVs on disk.

Each scientific decision still lives in a SOLID-compliant module behind a focused
class (``MagneticFieldCalibrator``, ``GasDensityEstimator``, ``BaryonicGapCalculator``,
``MagneticPressureGradientEstimator``, ``PlasmaDiagnostics``, ``ReliabilityMaskBuilder``)
in :mod:`src.physics`.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any, Iterable

logger = logging.getLogger(__name__)

from src.pipeline.constants import (
    DEFAULT_CONTROL_GALAXIES,
    RESULTS_DIRNAME,
    RESULTS_GRAPHS_DIRNAME,
    RESULTS_GRAPHS_POPULATION_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
    RESULTS_STEP_DISCOVERY_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    SPARC_DIR,
)
from src.pipeline.bootstrap import build_pipeline_dependencies
from src.pipeline.console_ui import (
    PipelineConsole,
    build_pipeline_console,
    format_pipeline_summary_section,
    pipeline_session_start_block,
    pipeline_terminal_columns,
)
from src.pipeline.workflows import (
    PhysicsDiscoveryOnly,
    assert_viz_only_prerequisites,
    run_physics_workflow,
    run_visualization_workflow,
)
from src.pipeline.aggregation.readme_builder import build_readme_manual
from src.providers import SparcMorphologyProvider
from src.reporting.galaxy_log import format_rule_banner


_RESERVED_GALAXY_KEYWORDS: frozenset[str] = frozenset({"all", "default"})


def _resolve_target_set(
    galaxies_tokens: Iterable[str],
    sparc_dir: Path,
    morphology_provider: SparcMorphologyProvider,
) -> list[Path]:
    """
    Resolve the ``--galaxies`` token list into a concrete list of ``*_rotmod.dat`` paths.

    Tokens are case-insensitive. The single tokens ``all`` and ``default`` are reserved
    presets; any other value is treated as an explicit galaxy name. Mixing presets with
    explicit names raises a CLI error so the operator does not silently get an unexpected
    galaxy set.
    """
    tokens = [str(t).strip() for t in galaxies_tokens if str(t).strip()]
    if not tokens:
        tokens = ["all"]

    lowered = [t.lower() for t in tokens]
    presets_present = [t for t in lowered if t in _RESERVED_GALAXY_KEYWORDS]
    names_present = [tokens[i] for i, lt in enumerate(lowered) if lt not in _RESERVED_GALAXY_KEYWORDS]

    if presets_present and names_present:
        raise SystemExit(
            "error: --galaxies cannot mix presets ('all'/'default') with explicit galaxy names; "
            f"got {tokens!r}"
        )
    if len(presets_present) > 1:
        raise SystemExit(
            f"error: --galaxies accepts at most one preset; got {presets_present!r}"
        )

    all_sparc = sorted(sparc_dir.glob("*_rotmod.dat"))
    if not all_sparc:
        raise FileNotFoundError(f"No SPARC galaxy files found under {sparc_dir}.")

    if presets_present == ["all"]:
        return all_sparc
    if presets_present == ["default"]:
        target_keys = {morphology_provider.normalize_galaxy_name(g) for g in DEFAULT_CONTROL_GALAXIES}
    else:
        target_keys = {morphology_provider.normalize_galaxy_name(g) for g in tokens}

    selected = [
        p
        for p in all_sparc
        if morphology_provider.normalize_galaxy_name(p.stem.replace("_rotmod", "")) in target_keys
    ]
    if not selected:
        raise FileNotFoundError(
            f"No SPARC galaxy files matched --galaxies {tokens!r} under {sparc_dir}."
        )
    return selected


def _short_pipeline_closing_lines(
    *,
    discovery_only: bool,
    processed_root: Path,
    repo_root: Path,
    n_sparc_targets: int,
    n_disc_reliable: int | None = None,
    n_census_rows: int | None = None,
    n_skipped_grid: int | None = None,
    n_grid_galaxies: int | None = None,
    n_direct_evidence: int | None = None,
    n_validation_reliable: int | None = None,
    rule_width: int | None = None,
) -> list[str]:
    """Terminal closing block aligned with :mod:`src.reporting.galaxy_log` rule aesthetics."""
    try:
        results_rel = str(processed_root.relative_to(repo_root))
    except ValueError:
        results_rel = str(processed_root)

    title = "Discovery phase complete" if discovery_only else "Pipeline Finished Successfully"
    term_w = rule_width if rule_width is not None else pipeline_terminal_columns()
    closing_banner = format_rule_banner(title, "results", term_w)

    folder_rows: list[tuple[str, str]] = [
        (RESULTS_GALAXY_ANALYSIS_DIRNAME, "tables, JSON and publication PNGs"),
        (f"{RESULTS_GRAPHS_DIRNAME}/{RESULTS_GRAPHS_POPULATION_DIRNAME}", "population diagnostic figures"),
        ("grid_results", "scenario knob combinations"),
        (RESULTS_STEP_DISCOVERY_DIRNAME, "census (all / reliable / not_reliable)"),
        (RESULTS_STEP_VALIDATION_DIRNAME, "evidence, anchor sweep and audit md"),
    ]
    key_width = max(len(k) for k, _ in folder_rows)
    output_detail_lines = [
        f"{results_rel}/",
        *[f"{key.ljust(key_width)}  |  {desc}" for key, desc in folder_rows],
    ]

    if discovery_only:
        assert n_disc_reliable is not None and n_census_rows is not None and n_skipped_grid is not None
        left1 = f"{n_sparc_targets} SPARC targets"
        right1 = f"{n_census_rows} discovery census rows"
        left2 = f"{n_disc_reliable} passing discovery gates"
        right2 = f"{n_skipped_grid} deferred for grid audits"
        left_w = max(len(left1), len(left2))
        stat_line_1 = f"{left1.ljust(left_w)}  |  {right1}"
        stat_line_2 = f"{left2.ljust(left_w)}  |  {right2}"
    else:
        assert n_grid_galaxies is not None
        assert n_direct_evidence is not None and n_validation_reliable is not None
        de = n_direct_evidence
        de_show = f"{de:02d}" if de < 100 else str(de)
        left1 = f"{n_sparc_targets} SPARC targets"
        right1 = f"{n_grid_galaxies} galaxies in validation grid"
        left2 = f"{n_validation_reliable} reliable validation rows"
        right2 = f"{de_show} direct-evidence galaxies"
        left_w = max(len(left1), len(left2))
        stat_line_1 = f"{left1.ljust(left_w)}  |  {right1}"
        stat_line_2 = f"{left2.ljust(left_w)}  |  {right2}"

    doc_detail_lines = [
        "docs/results_interpretation/Data_Interpretation_Guide_PT.md",
        "docs/results_interpretation/Data_Interpretation_Guide_EN.md",
    ]

    outputs_block = format_pipeline_summary_section(
        "Results — Outputs",
        "tree",
        output_detail_lines,
        rule_width=rule_width,
        use_detail_arrow=True,
    )
    stats_block = format_pipeline_summary_section(
        "Results — Summary",
        "stats",
        [stat_line_1, stat_line_2],
        rule_width=rule_width,
        use_detail_arrow=False,
    )
    docs_block = format_pipeline_summary_section(
        "Results — Documentation",
        "docs",
        doc_detail_lines,
        rule_width=rule_width,
        use_detail_arrow=True,
    )

    lines: list[str] = [
        "",
        outputs_block,
        "",
        stats_block,
        "",
        docs_block,
        "",
        closing_banner,
    ]
    return lines


def run_pipeline(
    args: argparse.Namespace | None = None,
    *,
    project_root: Path | None = None,
    console: PipelineConsole | None = None,
) -> None:
    if args is None:
        from src.pipeline.cli import parse_args

        args = parse_args()

    root = project_root if project_root is not None else Path(__file__).resolve().parent.parent.parent
    processed_root = root / RESULTS_DIRNAME

    ui = console if console is not None else build_pipeline_console(args)
    emit_logs = bool(getattr(args, "logs", False))
    tqdm_on = not emit_logs

    readme_path = root / "README.md"

    if bool(getattr(args, "viz_only", False)):
        ui.emit(pipeline_session_start_block())
        plot_names = assert_viz_only_prerequisites(root)
        viz_note = run_visualization_workflow(
            root,
            processed_root,
            args,
            ui,
            plot_names,
            tqdm_on=tqdm_on,
        )
        readme_path.write_text(build_readme_manual(), encoding="utf-8")
        tail = [viz_note] if viz_note else []
        ui.print_final_summary(
            [
                "",
                "Visualization-only run: publication figures and supreme-cohort 3D/HTML refreshed from disk.",
                "",
                *tail,
            ]
        )
        return

    with ui.planck_activity("Initializing Planck HEALPix maps"):
        deps = build_pipeline_dependencies(
            root,
            verbose_planck=emit_logs,
            emit_verbose_galaxy_logs=emit_logs,
        )

    sparc_files = _resolve_target_set(args.galaxies, root / SPARC_DIR, deps.morphology_provider)

    ui.emit(pipeline_session_start_block())

    if deps.icrs_supplement and emit_logs:
        ui.emit(
            f"Loaded ICRS supplement: {len(deps.icrs_supplement)} galaxy position(s) available for NED cone-search seeding.",
        )

    phys_out = run_physics_workflow(
        root,
        processed_root,
        args,
        ui,
        deps,
        sparc_files,
        tqdm_on=tqdm_on,
    )

    if isinstance(phys_out, PhysicsDiscoveryOnly):
        summary_lines = _short_pipeline_closing_lines(
            discovery_only=True,
            processed_root=processed_root,
            repo_root=root,
            n_sparc_targets=len(phys_out.sparc_files),
            n_disc_reliable=phys_out.n_disc_reliable,
            n_census_rows=len(phys_out.discovery_rows),
            n_skipped_grid=phys_out.n_skipped_grid,
        )
        ui.print_final_summary(summary_lines)
        return

    viz_note = ""
    if not bool(getattr(args, "compute_only", False)):
        viz_note = run_visualization_workflow(
            root,
            processed_root,
            args,
            ui,
            phys_out.plot_names,
            tqdm_on=tqdm_on,
        )

    readme_path.write_text(build_readme_manual(), encoding="utf-8")

    summary_lines = _short_pipeline_closing_lines(
        discovery_only=False,
        processed_root=processed_root,
        repo_root=root,
        n_sparc_targets=len(phys_out.sparc_files),
        n_grid_galaxies=len(phys_out.phase2_files),
        n_direct_evidence=len(phys_out.direct_evidence_rows),
        n_validation_reliable=phys_out.n_val_reliable,
    )
    if viz_note:
        summary_lines = list(summary_lines) + ["", viz_note]
    ui.print_final_summary(summary_lines)


def run_pipeline_with_root(root: Path, args: argparse.Namespace) -> None:
    """Test helper: ``chdir(root)`` so relative ``sources/`` paths resolve, then run the pipeline."""
    import os

    prev = Path.cwd()
    try:
        os.chdir(root)
        run_pipeline(args, project_root=root.resolve())
    finally:
        os.chdir(prev)
