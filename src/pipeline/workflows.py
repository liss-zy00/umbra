"""
Physics vs visualization workflow split for :func:`src.pipeline.runner.run_pipeline`.

Physics covers Discovery, Validation, and the population **Analysis** report (Phase 3).
Visualization covers publication figures (Phase 4) and supreme-cohort 3D + HTML (Phase 5).
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Sequence

import numpy as np

from src.physics import GRADIENT_SCALES_KPC, IONIZATION_FRACTIONS
from src.pipeline.constants import (
    RESULTS_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
    RESULTS_GRAPHS_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_RELIABLE_FILENAME,
)
from src.pipeline.console_ui import (
    PIPELINE_PHASE_3_ANALYSIS_TITLE,
    PIPELINE_PHASE_4_GRAPHS_TITLE,
    PIPELINE_PHASE_5_VIZ_TITLE,
    PipelineConsole,
    pipeline_phase_2_validation_title,
)
from src.pipeline.per_galaxy.discovery_processor import discovery_audits_clean, process_galaxy_discovery
from src.pipeline.per_galaxy.galaxy_processor import (
    MutablePipelineState,
    PipelineDependencies,
    process_sparc_file,
)
from src.pipeline.aggregation.global_aggregator import build_direct_evidence_rows
from src.reporting.galaxy_log import (
    emit_validation_grid_progress_finalize,
    format_phase3_graphs_verbose_block,
    format_phase4_analysis_verbose_block,
    format_validation_phase_scenario_legend,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PhysicsDiscoveryOnly:
    """Returned when no galaxies pass the Discovery pre-filter for the validation grid."""

    sparc_files: list[Path]
    discovery_rows: list[dict[str, Any]]
    n_disc_reliable: int
    n_skipped_grid: int


@dataclass(frozen=True)
class PhysicsWorkflowResult:
    """Successful physics path including validation grid and population Analysis report."""

    sparc_files: list[Path]
    phase2_files: list[Path]
    plot_names: list[str]
    discovery_rows: list[dict[str, Any]]
    n_disc_reliable: int
    n_skipped_grid: int
    state: MutablePipelineState
    direct_evidence_rows: list[dict[str, Any]]
    n_val_reliable: int


def _run_discovery_block(
    deps: PipelineDependencies,
    sparc_files: list[Path],
    ui: PipelineConsole,
    *,
    tqdm_on: bool,
    persister: Any,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], Any, int]:
    from src.pipeline.console_ui import PIPELINE_PHASE_1_DISCOVERY_TITLE
    from src.pipeline.constants import RESULTS_STEP_DISCOVERY_DIRNAME

    ui.phase_start(PIPELINE_PHASE_1_DISCOVERY_TITLE)
    discovery_rows: list[dict[str, Any]] = []
    for file_path in ui.iter_progress(
        sparc_files,
        desc="",
        total=len(sparc_files),
        leave=False,
        enabled=tqdm_on,
    ):
        row = process_galaxy_discovery(deps, file_path)
        if row is not None:
            discovery_rows.append(row)
    *_, elite_band_info = persister.save_discovery_evidence(discovery_rows)
    discovery_lookup: dict[str, dict[str, Any]] = {
        str(r["galaxy_name"]): r for r in discovery_rows if r.get("galaxy_name") is not None
    }
    reliable_disc = [
        r for r in discovery_rows if bool(r.get("discovery_is_reliable", r.get("is_reliable")))
    ]
    n_disc_reliable = len(reliable_disc)
    ui.phase_end(
        PIPELINE_PHASE_1_DISCOVERY_TITLE,
        detail=(
            f"{len(sparc_files)} targets processed\u202f\u00b7\u202f{n_disc_reliable} passed minimum reliability gates"
            f" · Census data saved to {RESULTS_DIRNAME}/{RESULTS_STEP_DISCOVERY_DIRNAME}/"
        ),
    )
    return discovery_rows, discovery_lookup, elite_band_info, n_disc_reliable


def run_physics_workflow(
    root: Path,
    processed_root: Path,
    args: argparse.Namespace,
    ui: PipelineConsole,
    deps: PipelineDependencies,
    sparc_files: list[Path],
    *,
    tqdm_on: bool,
) -> PhysicsWorkflowResult | PhysicsDiscoveryOnly:
    """
    Phases 1–3: Discovery, Validation, population Analysis (consistency markdown).

    Publication figures and 3D exports are excluded (see :func:`run_visualization_workflow`).
    """
    persister = deps.persister
    discovery_rows, discovery_lookup, elite_band_info, n_disc_reliable = _run_discovery_block(
        deps, sparc_files, ui, tqdm_on=tqdm_on, persister=persister
    )

    phase2_files: list[Path] = []
    for p in sparc_files:
        gname = p.stem.replace("_rotmod", "")
        drow = discovery_lookup.get(gname, {"audit_flags": ""})
        if discovery_audits_clean(drow):
            phase2_files.append(p)
    n_skipped_grid = len(sparc_files) - len(phase2_files)

    if not phase2_files:
        return PhysicsDiscoveryOnly(
            sparc_files=sparc_files,
            discovery_rows=discovery_rows,
            n_disc_reliable=n_disc_reliable,
            n_skipped_grid=n_skipped_grid,
        )

    scenario_rows: dict[tuple[float, float], list[dict[str, Any]]] = {
        (s, i): [] for s in GRADIENT_SCALES_KPC for i in IONIZATION_FRACTIONS
    }
    scenario_success_counts = np.zeros((len(GRADIENT_SCALES_KPC), len(IONIZATION_FRACTIONS)), dtype=int)
    state = MutablePipelineState(
        scenario_rows=scenario_rows,
        scenario_success_counts=scenario_success_counts,
    )
    n_scenario_matrix = len(GRADIENT_SCALES_KPC) * len(IONIZATION_FRACTIONS)
    deps.validation_grid_total_galaxies = len(phase2_files)
    deps.validation_discovery_reliable_count = n_disc_reliable
    deps.validation_completed_names.clear()

    phase2_title = pipeline_phase_2_validation_title(n_scenario_matrix)
    ui.phase_start(phase2_title)
    if deps.emit_verbose_galaxy_logs:
        scenario_pairs = [(float(s), float(i)) for s in GRADIENT_SCALES_KPC for i in IONIZATION_FRACTIONS]
        print(
            format_validation_phase_scenario_legend(
                scenario_pairs,
                n_galaxies=len(phase2_files),
            ),
            end="",
            flush=True,
        )
    for file_path in ui.iter_progress(
        phase2_files,
        desc="",
        total=len(phase2_files),
        leave=False,
        enabled=tqdm_on,
    ):
        process_sparc_file(deps, state, file_path)

    if deps.emit_verbose_galaxy_logs:
        emit_validation_grid_progress_finalize(
            enabled=True,
            completed_names=list(deps.validation_completed_names),
        )

    for (scale, ion), rows in scenario_rows.items():
        scenario_name = f"scenario_scale_{scale:.3f}_ionization_{ion:.3f}"
        persister.save_grid_scenario_report(scenario_name, rows)

    persister.save_full_catalog_inventory(state.master_inventory_rows)
    persister.save_validation_reliability_partitions(state.master_inventory_rows)

    direct_evidence_rows = build_direct_evidence_rows(
        state.supreme_rows,
        deps.morphology_lookup,
        deps.morphology_provider,
        discovery_lookup=discovery_lookup,
        elite_band_info=elite_band_info,
    )
    persister.save_scientific_evidence_report(direct_evidence_rows)

    n_val_reliable = sum(
        1
        for row in state.master_inventory_rows
        if bool(row.get("validation_is_reliable", row.get("is_reliable")))
    )
    n_eval = n_scenario_matrix * len(phase2_files)
    ui.phase_end(
        phase2_title,
        detail=(
            f"Matrix: {n_scenario_matrix} scenarios per galaxy\u202f\u00b7\u202f{n_eval:,} evaluations"
            f" · {n_val_reliable} galaxies verified as reliable\u202f\u00b7\u202f{len(direct_evidence_rows)} supreme matches identified"
            f" · Physical bundles synced to {RESULTS_DIRNAME}/{RESULTS_GALAXY_ANALYSIS_DIRNAME}/"
        ),
    )

    from src.physics.diagnostics.anchor_sensitivity import build_anchor_sensitivity_rows

    base_ratios: dict[str, float] = {}
    for row in state.supreme_rows:
        gname = str(row.get("galaxy_name") or "")
        if not gname:
            continue
        try:
            base_ratios[gname] = float(row.get("theoretical_gap_coverage_ratio"))
        except (TypeError, ValueError):
            base_ratios[gname] = float("nan")
    anchor_rows = build_anchor_sensitivity_rows(base_ratios)
    persister.save_anchor_sensitivity_sweep(anchor_rows)

    # Phase 3 — Analysis (population consistency; runs before publication so summaries stay in reliable tier paths).
    from src.pipeline.aggregation.audit import count_galaxy_summary_json_paths, write_consistency_audit_report

    ui.phase_start(PIPELINE_PHASE_3_ANALYSIS_TITLE)
    if deps.emit_verbose_galaxy_logs:
        ga_root = processed_root / RESULTS_GALAXY_ANALYSIS_DIRNAME
        n_sum = count_galaxy_summary_json_paths(ga_root)
        print(
            format_phase4_analysis_verbose_block(
                audit_rel=f"{RESULTS_DIRNAME}/{RESULTS_STEP_VALIDATION_DIRNAME}/validation_consistency_audit.md",
                n_summary_scan=n_sum,
            ),
            end="",
            flush=True,
        )
    write_consistency_audit_report(
        processed_root=processed_root,
        supreme_rows=state.supreme_rows,
        master_inventory_rows=state.master_inventory_rows,
    )
    ui.phase_end(
        PIPELINE_PHASE_3_ANALYSIS_TITLE,
        detail=(
            "Aggregating supreme rows and metadata synthesis"
            f" · Consistency analysis delivered: {RESULTS_DIRNAME}/{RESULTS_STEP_VALIDATION_DIRNAME}/validation_consistency_audit.md"
        ),
    )

    plot_names = sorted({p.stem.replace("_rotmod", "") for p in phase2_files})
    return PhysicsWorkflowResult(
        sparc_files=sparc_files,
        phase2_files=phase2_files,
        plot_names=plot_names,
        discovery_rows=discovery_rows,
        n_disc_reliable=n_disc_reliable,
        n_skipped_grid=n_skipped_grid,
        state=state,
        direct_evidence_rows=direct_evidence_rows,
        n_val_reliable=n_val_reliable,
    )


def _publication_progress(ui: PipelineConsole, names: Sequence[str], *, tqdm_on: bool) -> Iterator[str]:
    seq = list(names)
    return ui.iter_progress(
        seq,
        desc="",
        total=len(seq),
        leave=False,
        enabled=tqdm_on,
    )


def run_visualization_workflow(
    root: Path,
    processed_root: Path,
    args: argparse.Namespace,
    ui: PipelineConsole,
    plot_names: Sequence[str],
    *,
    tqdm_on: bool,
) -> str:
    """
    Phases 4–5: publication figures then supreme-cohort 3D PNG + interactive HTML.

    Returns a short note for the final summary (may be empty).
    """
    _ = processed_root  # reserved for future path checks

    def _pub_progress(names: Sequence[str]) -> Iterator[str]:
        return _publication_progress(ui, names, tqdm_on=tqdm_on)

    ui.phase_start(PIPELINE_PHASE_4_GRAPHS_TITLE)
    if getattr(args, "logs", False):
        print(
            format_phase3_graphs_verbose_block(
                n_targets=len(plot_names),
                graphs_rel=f"{RESULTS_DIRNAME}/{RESULTS_GRAPHS_DIRNAME}",
            ),
            end="",
            flush=True,
        )
    try:
        from src.reporting.publication_bundle import run_publication_figures_for_run

        run_publication_figures_for_run(
            root,
            plot_names,
            write_companion_md=False,
            per_galaxy_progress=_pub_progress,
        )
    except Exception as exc:
        logger.warning("Publication figures skipped: %s", exc)
    ui.phase_end(
        PIPELINE_PHASE_4_GRAPHS_TITLE,
        detail=(
            "Exporting scientific figures & population diagnostics"
            f" · Plots stored in {RESULTS_DIRNAME}/{RESULTS_GRAPHS_DIRNAME}/"
        ),
    )

    ui.phase_start(PIPELINE_PHASE_5_VIZ_TITLE)
    from src.visualization.batch import run_supreme_museum_batch

    note = run_supreme_museum_batch(project_root=root, ui=ui)
    ui.phase_end(
        PIPELINE_PHASE_5_VIZ_TITLE,
        detail="Supreme-cohort 3D PNG + interactive HTML under results/graphs/viz_3d/ (when PyVista is available).",
    )
    return note


def assert_viz_only_prerequisites(project_root: Path) -> list[str]:
    """
    Ensure validation artefacts exist for a standalone visualization re-run.

    Returns the galaxy name list to pass to the publication bundle (validation-reliable cohort).
    """
    rel = project_root / RESULTS_DIRNAME / RESULTS_STEP_VALIDATION_DIRNAME / VALIDATION_RELIABLE_FILENAME
    if not rel.is_file():
        raise SystemExit(
            "error: --viz-only requires an existing validation cohort table at "
            f"{rel.as_posix()} (run the physics workflow first)."
        )
    from src.reporting.publication_bundle import read_cohort_galaxy_names

    return read_cohort_galaxy_names(project_root, "validation_reliable")
