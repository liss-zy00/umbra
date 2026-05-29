"""Post-pipeline publication figures (per-galaxy PNGs beside CSV/JSON; cohort panels under ``graphs/``)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable, Iterator, Literal, Sequence

import numpy as np
import pandas as pd

from src.pipeline.constants import (
    DISCOVERY_RELIABLE_FILENAME,
    RESULTS_DIRNAME,
    RESULTS_GRAPHS_DIRNAME,
    RESULTS_GRAPHS_POPULATION_DIRNAME,
    RESULTS_STEP_DISCOVERY_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_RELIABLE_FILENAME,
    VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME,
)
from src.reporting.galaxy_from_disk import (
    analysis_csv_path,
    build_galaxy_from_galaxy_logs,
    resolve_per_galaxy_analysis_dir,
)
from src.reporting.plot_style import apply_academic_mpl_style
from src.reporting.plotting.plotter import GalaxyPlotter
from src.reporting.population_figures import (
    save_discovery_rho_req_vs_wim_ratio,
    save_evidence_coverage_vs_pearson,
    save_evidence_rar_vs_btfr_offsets,
    save_pearson_histogram_from_evidence,
    save_resolution_diagnostic_scatter,
    save_validation_mean_gap_vs_hubble,
)

logger = logging.getLogger(__name__)

CohortKey = Literal["discovery_reliable", "validation_reliable", "scientific_evidence"]

_COHORT_SPEC: dict[CohortKey, tuple[str, str]] = {
    "discovery_reliable": (RESULTS_STEP_DISCOVERY_DIRNAME, DISCOVERY_RELIABLE_FILENAME),
    "validation_reliable": (RESULTS_STEP_VALIDATION_DIRNAME, VALIDATION_RELIABLE_FILENAME),
    "scientific_evidence": (RESULTS_STEP_VALIDATION_DIRNAME, VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME),
}


def cohort_csv_path(project_root: Path, cohort: CohortKey) -> Path:
    subdir, fname = _COHORT_SPEC[cohort]
    return project_root / RESULTS_DIRNAME / subdir / fname


def read_cohort_galaxy_names(project_root: Path, cohort: CohortKey) -> list[str]:
    path = cohort_csv_path(project_root, cohort)
    if not path.is_file():
        raise FileNotFoundError(f"Cohort CSV not found: {path}")
    df = pd.read_csv(path)
    if "galaxy_name" not in df.columns:
        raise ValueError(f"Column galaxy_name missing in {path}")
    names = [str(x).strip() for x in df["galaxy_name"].tolist() if str(x).strip()]
    return sorted(set(names))


def format_summary_metadata_footer(summary: dict[str, Any] | None) -> str | None:
    """One-line provenance string from ``*_summary.json`` for embedding in figure captions."""
    if not summary:
        return None
    best = summary.get("best_theoretical_scenario") or {}
    unc = summary.get("uncertainty_propagation") or {}
    parts: list[str] = []
    sn = best.get("scenario_name")
    if sn:
        parts.append(f"Best scenario: {sn}")
    r = best.get("theoretical_gap_coverage_ratio")
    if r is not None:
        try:
            rf = float(r)
            if np.isfinite(rf):
                parts.append(f"Gap coverage={rf:.3f}")
        except (TypeError, ValueError):
            pass
    p16, p84 = unc.get("theoretical_gap_coverage_ratio_p16"), unc.get("theoretical_gap_coverage_ratio_p84")
    if p16 is not None and p84 is not None:
        try:
            lo, hi = float(p16), float(p84)
            if np.isfinite(lo) and np.isfinite(hi):
                parts.append(f"MC p16–p84: {lo:.3f}–{hi:.3f}")
        except (TypeError, ValueError):
            pass
    p = best.get("pearson_correlation")
    if p is not None:
        try:
            pf = float(p)
            if np.isfinite(pf):
                parts.append(f"Pearson r={pf:.3f}")
        except (TypeError, ValueError):
            pass
    return " · ".join(parts) if parts else None


def compute_cohort_force_y_limits(project_root: Path, galaxy_names: list[str]) -> tuple[float, float] | None:
    """Shared log-axis limits for |gap| and |a_B| across the cohort (2–98 percentile of log10)."""
    vals: list[float] = []
    for name in galaxy_names:
        p = analysis_csv_path(project_root, name)
        if not p.is_file():
            continue
        try:
            df = pd.read_csv(p, usecols=["acceleration_gap", "magnetic_pressure_acceleration"])
        except ValueError:
            df = pd.read_csv(p)
            if "acceleration_gap" not in df.columns or "magnetic_pressure_acceleration" not in df.columns:
                continue
        ga = np.abs(df["acceleration_gap"].to_numpy(dtype=float))
        mb = np.abs(df["magnetic_pressure_acceleration"].to_numpy(dtype=float))
        for arr in (ga, mb):
            v = arr[np.isfinite(arr) & (arr > 0.0)]
            vals.extend(v.astype(float).tolist())
    if len(vals) < 3:
        return None
    a = np.log10(np.asarray(vals, dtype=float))
    lo, hi = float(np.nanquantile(a, 0.02)), float(np.nanquantile(a, 0.98))
    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        return None
    return (10.0**lo, 10.0**hi)


def compute_cohort_rar_log_limits(project_root: Path, galaxy_names: list[str]) -> tuple[tuple[float, float], tuple[float, float]] | None:
    """Shared log limits for RAR scatter (g_bar on x, max(g_obs, a_B) on y)."""
    xs: list[float] = []
    ys: list[float] = []
    for name in galaxy_names:
        p = analysis_csv_path(project_root, name)
        if not p.is_file():
            continue
        df = pd.read_csv(p)
        need = {"radius_kpc", "observed_velocity_kms", "baryonic_velocity_kms"}
        if not need.issubset(df.columns):
            continue
        radius_m = df["radius_kpc"].to_numpy(dtype=float) * 3.0857e19
        v_obs_ms = df["observed_velocity_kms"].to_numpy(dtype=float) * 1000.0
        v_bar_ms = df["baryonic_velocity_kms"].to_numpy(dtype=float) * 1000.0
        with np.errstate(divide="ignore", invalid="ignore"):
            g_obs = np.where(radius_m > 0, v_obs_ms**2 / radius_m, np.nan)
            g_bar = np.where(radius_m > 0, v_bar_ms**2 / radius_m, np.nan)
        ok = np.isfinite(g_obs) & np.isfinite(g_bar) & (g_obs > 0) & (g_bar > 0)
        xs.extend(g_bar[ok].tolist())
        ys.extend(g_obs[ok].tolist())
        if "magnetic_pressure_acceleration" in df.columns:
            a_b = np.abs(df["magnetic_pressure_acceleration"].to_numpy(dtype=float))
            m2 = ok & np.isfinite(a_b) & (a_b > 0)
            ys.extend(a_b[m2].tolist())
            xs.extend(g_bar[m2].tolist())
    if len(xs) < 3:
        return None
    lx = np.log10(np.asarray(xs, dtype=float))
    ly = np.log10(np.asarray(ys, dtype=float))
    xlo, xhi = float(np.nanquantile(lx, 0.02)), float(np.nanquantile(lx, 0.98))
    ylo, yhi = float(np.nanquantile(ly, 0.02)), float(np.nanquantile(ly, 0.98))
    if not all(np.isfinite([xlo, xhi, ylo, yhi])) or xlo >= xhi or ylo >= yhi:
        return None
    pad = 0.12
    return ((10.0 ** (xlo - pad), 10.0 ** (xhi + pad)), (10.0 ** (ylo - pad), 10.0 ** (yhi + pad)))


def run_publication_figures_for_run(
    project_root: Path,
    galaxy_names: Sequence[str],
    *,
    population_dirname: str = RESULTS_GRAPHS_POPULATION_DIRNAME,
    write_companion_md: bool = False,
    per_galaxy_progress: Callable[[Sequence[str]], Iterator[str]] | None = None,
) -> list[Path]:
    """
    Build per-galaxy PNGs next to each galaxy's analysis CSV/JSON plus cohort panels.

    Per-galaxy PNGs are written only for galaxies listed in ``step_validation/validation_reliable.csv``
    (intersection with the passed ``galaxy_names``); the not-reliable tier receives no PNG pass.
    After population panels, galaxies listed in ``validation_scientific_evidence_report.csv`` are
    moved from ``galaxy_analysis/validation_reliable/`` to ``galaxy_analysis/validation_scientific_evidence/``.

    Population PNGs are written under ``results/graphs/<population_dirname>/`` (default
    ``population`` for a full pipeline run).
    """
    apply_academic_mpl_style()
    names = sorted({str(n).strip() for n in galaxy_names if str(n).strip()})
    graphs_root = project_root / RESULTS_DIRNAME / RESULTS_GRAPHS_DIRNAME
    pop_dir = graphs_root / population_dirname
    pop_dir.mkdir(parents=True, exist_ok=True)

    val_rel_path = project_root / RESULTS_DIRNAME / RESULTS_STEP_VALIDATION_DIRNAME / VALIDATION_RELIABLE_FILENAME
    validation_reliable_set: set[str] = set()
    if val_rel_path.is_file():
        try:
            validation_reliable_set = set(read_cohort_galaxy_names(project_root, "validation_reliable"))
        except (OSError, ValueError, FileNotFoundError) as exc:
            logger.warning("Could not load validation reliable cohort for per-galaxy PNG filter: %s", exc)
    per_galaxy_names = [n for n in names if n in validation_reliable_set] if validation_reliable_set else list(names)

    y_force = compute_cohort_force_y_limits(project_root, per_galaxy_names)
    rar_limits = compute_cohort_rar_log_limits(project_root, per_galaxy_names)
    x_rar = rar_limits[0] if rar_limits else None
    y_rar = rar_limits[1] if rar_limits else None

    plotter = GalaxyPlotter()
    written: list[Path] = []

    name_iter: Iterator[str] = (
        per_galaxy_progress(per_galaxy_names) if per_galaxy_progress is not None else iter(per_galaxy_names)
    )
    for name in name_iter:
        try:
            galaxy, summary = build_galaxy_from_galaxy_logs(project_root, name)
            gdir = resolve_per_galaxy_analysis_dir(project_root, name)
        except (OSError, ValueError, KeyError, FileNotFoundError) as exc:
            logger.warning("Skip %s: %s", name, exc)
            continue
        best = (summary or {}).get("best_theoretical_scenario") or {}
        sk = best.get("scenario_gradient_scale_kpc")
        ion = best.get("scenario_ionization_fraction")
        sk_f = float(sk) if sk is not None and np.isfinite(float(sk)) else None
        ion_f = float(ion) if ion is not None and np.isfinite(float(ion)) else None
        k_raw = best.get("resolution_correction_factor")
        k_plot_f = float(k_raw) if k_raw is not None and np.isfinite(float(k_raw)) else None
        written.append(
            plotter.save_force_plot(
                galaxy,
                gdir,
                scenario_gradient_scale_kpc=sk_f,
                scenario_ionization_fraction=ion_f,
                resolution_correction_factor_k=k_plot_f,
                write_companion=write_companion_md,
                y_accel_log_limits=y_force,
            )
        )
        written.append(
            plotter.save_rar_plot(
                galaxy,
                gdir / f"{name}_rar_plot.png",
                write_companion=write_companion_md,
                x_log_limits=x_rar,
                y_log_limits=y_rar,
            )
        )
        written.append(
            plotter.save_acceleration_vs_radius_sample(
                galaxy,
                gdir / f"{name}_signed_accel_vs_radius.png",
                write_companion=write_companion_md,
            )
        )

    disc_path = project_root / RESULTS_DIRNAME / RESULTS_STEP_DISCOVERY_DIRNAME / DISCOVERY_RELIABLE_FILENAME
    if disc_path.is_file():
        ddf = pd.read_csv(disc_path)
        written.append(
            save_discovery_rho_req_vs_wim_ratio(
                ddf, pop_dir / "discovery_rho_req_vs_wim_ratio.png", write_companion=write_companion_md
            )
        )
        written.append(
            save_resolution_diagnostic_scatter(
                ddf, pop_dir / "discovery_beam_vs_resolution_factor.png", write_companion=write_companion_md
            )
        )

    ev_path = project_root / RESULTS_DIRNAME / RESULTS_STEP_VALIDATION_DIRNAME / VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME
    if ev_path.is_file():
        edf = pd.read_csv(ev_path)
        written.append(
            save_evidence_coverage_vs_pearson(
                edf, pop_dir / "evidence_coverage_vs_pearson.png", write_companion=write_companion_md
            )
        )
        written.append(
            save_evidence_rar_vs_btfr_offsets(
                edf, pop_dir / "evidence_rar_vs_btfr_offsets.png", write_companion=write_companion_md
            )
        )
        written.append(
            save_pearson_histogram_from_evidence(
                edf, pop_dir / "evidence_pearson_histogram.png", write_companion=write_companion_md
            )
        )

    val_path = project_root / RESULTS_DIRNAME / RESULTS_STEP_VALIDATION_DIRNAME / VALIDATION_RELIABLE_FILENAME
    if val_path.is_file():
        vdf = pd.read_csv(val_path)
        written.append(
            save_validation_mean_gap_vs_hubble(
                vdf, pop_dir / "validation_mean_gap_vs_hubble.png", write_companion=write_companion_md
            )
        )

    from src.reporting.galaxy_analysis_layout import relocate_scientific_evidence_galaxy_dirs

    relocate_scientific_evidence_galaxy_dirs(project_root)

    return written


def run_publication_bundle(
    project_root: Path,
    cohort: CohortKey,
    *,
    write_companion_md: bool = False,
) -> list[Path]:
    """
    Generate per-galaxy PNGs plus population diagnostics for one cohort table.

    Per-galaxy PNGs are still filtered to ``validation_reliable.csv`` when that file exists.
    Cohort panels use ``results/graphs/population_<cohort>/`` so alternate cohort bundles do
    not overwrite the default ``population/`` folder from a full run.
    """
    names = read_cohort_galaxy_names(project_root, cohort)
    return run_publication_figures_for_run(
        project_root,
        names,
        population_dirname=f"{RESULTS_GRAPHS_POPULATION_DIRNAME}_{cohort}",
        write_companion_md=write_companion_md,
    )
