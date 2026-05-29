"""Post-run aggregates for population-level direct-evidence synthesis.

``build_direct_evidence_rows`` merges Discovery and validation metrics. Publication-quality
figures are written by :mod:`src.reporting.publication_bundle` at the end of
:func:`src.pipeline.runner.run_pipeline`.

The data-dictionary JSON output was retired: the **conceptual** meaning of each column lives
in ``EG_MODEL_THEORY.md`` (timeless theory) while the **operational** column listing is
regenerated into ``README.md`` per run.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from src.pipeline.constants import PHOTOMETRIC_LSB_LABEL
from src.pipeline.per_galaxy.morphology import include_in_scientific_evidence_report, is_photometric_lsb
from src.physics import resolution_correction_factor
from src.providers.sparc import SparcMorphologyProvider

from src.config.thresholds import (
    RELIABILITY_CUTS,
    wim_density_match_score,
)


def _validation_status_supreme(row: dict[str, Any]) -> str:
    """
    ``SUPREME`` labels the high-provenance validation cohort: EDD-backed distance, finite PA,
    no blocking ``fail_*`` audits on the merged row, no BTFR outlier flag, beam not larger than
    the rotation-curve diameter proxy, and ``calibrated_convergence_ratio`` within one dex of unity.

    ``DATA_LIMITED`` marks systems where the Planck beam footprint exceeds the disk diameter
    inferred from ``r_max_kpc`` — the gradient reconstruction cannot claim high-fidelity structure
    inside an unresolved object.
    """
    beam = row.get("beam_size_kpc")
    r_max = row.get("r_max_kpc")
    try:
        beam_f = float(beam)
        r_max_f = float(r_max)
    except (TypeError, ValueError):
        beam_f = float("nan")
        r_max_f = float("nan")
    diam = float(RELIABILITY_CUTS.GALAXY_DIAMETER_FROM_RMAX_MULTIPLIER) * r_max_f
    if (
        np.isfinite(beam_f)
        and np.isfinite(r_max_f)
        and r_max_f > 0.0
        and np.isfinite(diam)
        and diam > 0.0
        and beam_f > diam
    ):
        return "DATA_LIMITED"

    if not bool(row.get("prov_high_precision_dist")):
        return "STANDARD"
    if bool(row.get("warn_no_pa")):
        return "STANDARD"
    if bool(row.get("warn_btfr_outlier")):
        return "STANDARD"
    flags = row.get("audit_flags") or []
    if any(str(c).startswith("fail_") for c in flags):
        return "STANDARD"
    pa = row.get("position_angle_deg")
    try:
        pa_f = float(pa)
    except (TypeError, ValueError):
        return "STANDARD"
    if not np.isfinite(pa_f):
        return "STANDARD"

    cal = row.get("calibrated_convergence_ratio")
    try:
        cal_f = float(cal)
    except (TypeError, ValueError):
        cal_f = float("nan")
    lo = float(RELIABILITY_CUTS.SUPREME_CALIBRATED_CONVERGENCE_RATIO_MIN)
    hi = float(RELIABILITY_CUTS.SUPREME_CALIBRATED_CONVERGENCE_RATIO_MAX)
    if not np.isfinite(cal_f) or cal_f <= 0.0 or cal_f < lo or cal_f > hi:
        return "STANDARD"
    return "SUPREME"


def _direct_evidence_sort_key(row: dict[str, Any]) -> tuple[int, float]:
    """Prefer ordering by WIM-calibrated convergence (physical anchor); fall back to gap distance."""
    cal = row.get("calibrated_convergence_ratio")
    try:
        cf = float(cal)
    except (TypeError, ValueError):
        cf = float("nan")
    gap_dist = row.get("ratio_distance_to_unity")
    try:
        gd = float(gap_dist)
    except (TypeError, ValueError):
        gd = float("nan")
    if np.isfinite(cf) and cf > 0.0:
        return (
            0,
            abs(
                np.log10(
                    max(cf, float(RELIABILITY_CUTS.BEST_SCENARIO_LOG_FLOOR)),
                )
            ),
        )
    return (1, gd if np.isfinite(gd) else 1.0e9)


def build_direct_evidence_rows(
    supreme_rows: list[dict[str, Any]],
    morphology_lookup: dict[str, dict[str, Any]],
    morphology_provider: SparcMorphologyProvider,
    *,
    discovery_lookup: dict[str, dict[str, Any]] | None = None,
    elite_band_info: "DiscoveryEliteBandInfo | None" = None,
) -> list[dict[str, Any]]:
    """
    Build the merged "direct-evidence" rows that feed ``validation_scientific_evidence_report.csv``.

    Each row consists of the best-scenario grid output augmented with the matching
    Discovery columns (``required_density_for_match``, beam-smearing diagnostics, etc.)
    plus the boolean flags ``discovery_is_reliable`` and ``discovery_is_elite``. The
    ``discovery_is_elite`` flag follows **only** the **dynamic band** defined by
    :class:`DiscoveryEliteBandInfo`: a galaxy is flagged when its
    ``required_density_for_match`` lies within ``[band_min_kg_m3, band_max_kg_m3]`` —
    a multiplicative window around the population median computed over the reliable
    Discovery subset. When the band is not valid (e.g. no reliable Discovery rows),
    ``discovery_is_elite`` is ``False`` and ``discovery_elite_band_source`` is ``none``.
    """
    from src.reporting.persistence.exporter import DiscoveryEliteBandInfo

    rows: list[dict[str, Any]] = []
    for row in supreme_rows:
        if not include_in_scientific_evidence_report(row, morphology_lookup, morphology_provider):
            continue
        ratio = row.get("theoretical_gap_coverage_ratio")
        if ratio is None or not np.isfinite(float(ratio)):
            continue
        ratio_f = float(ratio)
        if (
            RELIABILITY_CUTS.DIRECT_EVIDENCE_RATIO_MIN
            <= ratio_f
            <= RELIABILITY_CUTS.DIRECT_EVIDENCE_RATIO_MAX
        ):
            merged = {**row, "ratio_distance_to_unity": abs(np.log10(ratio_f))}
            name = str(row.get("galaxy_name", ""))
            drow = discovery_lookup.get(name) if discovery_lookup else None
            if drow:
                merged["required_density_for_match"] = drow.get("required_density_for_match")
                merged["convergence_ratio_to_wim"] = drow.get("convergence_ratio_to_wim")
                merged["delta_r_m"] = drow.get("delta_r_m")
                merged["physical_pixel_scale_kpc"] = drow.get("physical_pixel_scale_kpc")
                merged["distance_mpc"] = drow.get("distance_mpc")
                # Beam-smearing diagnostic: ``D_beam = D · θ_beam`` and the dilution
                # factor ``k = D_beam / L_grad`` are calculated in
                # :mod:`src.physics.diagnostics.resolution`. Two flavours are exposed:
                #   • ``resolution_correction_factor_discovery`` — k against the
                #     **data-derived** physical pixel separation (Δx_phys, Discovery).
                #   • ``resolution_correction_factor_grid`` — k against the **grid**
                #     scenario knob ``scenario_gradient_scale_kpc`` picked for this row.
                # The legacy ``resolution_correction_factor`` column keeps the grid value
                # so consumers reading "how much is the magnetic gradient diluted at the
                # ~100 pc scale assumed by the engine" see the meaningful number directly.
                merged["beam_size_kpc"] = drow.get("beam_size_kpc")
                merged["resolution_correction_factor_discovery"] = drow.get(
                    "resolution_correction_factor"
                )
                # Recompute against the best-scenario grid scale exposed on the supreme row.
                grid_scale_kpc = row.get("scenario_gradient_scale_kpc")
                k_grid = resolution_correction_factor(
                    distance_mpc=drow.get("distance_mpc"),
                    gradient_scale_kpc=grid_scale_kpc,
                )
                merged["resolution_correction_factor_grid"] = (
                    float(k_grid) if np.isfinite(k_grid) else float("nan")
                )
                merged["resolution_correction_factor"] = merged[
                    "resolution_correction_factor_grid"
                ]
                merged["discovery_is_reliable"] = bool(
                    drow.get("discovery_is_reliable", drow.get("is_reliable"))
                )
                # Dynamic cohort band on ``required_density_for_match`` only (no alternate rule).
                rho_value = drow.get("required_density_for_match")
                try:
                    rho_f = float(rho_value)
                except (TypeError, ValueError):
                    rho_f = float("nan")
                if (
                    elite_band_info is not None
                    and elite_band_info.is_valid
                    and np.isfinite(elite_band_info.band_min_kg_m3)
                    and np.isfinite(elite_band_info.band_max_kg_m3)
                ):
                    merged["discovery_is_elite"] = bool(
                        drow.get("discovery_is_reliable", drow.get("is_reliable"))
                        and np.isfinite(rho_f)
                        and rho_f > 0.0
                        and elite_band_info.band_min_kg_m3 <= rho_f <= elite_band_info.band_max_kg_m3
                    )
                    merged["discovery_elite_band_source"] = "dynamic_rho_band"
                else:
                    merged["discovery_is_elite"] = False
                    merged["discovery_elite_band_source"] = "none"
            else:
                merged["required_density_for_match"] = float("nan")
                merged["convergence_ratio_to_wim"] = float("nan")
                merged["delta_r_m"] = float("nan")
                merged["physical_pixel_scale_kpc"] = float("nan")
                merged["distance_mpc"] = float("nan")
                merged["beam_size_kpc"] = float("nan")
                merged["resolution_correction_factor"] = float("nan")
                merged["resolution_correction_factor_discovery"] = float("nan")
                merged["resolution_correction_factor_grid"] = float("nan")
                merged["discovery_is_reliable"] = False
                merged["discovery_is_elite"] = False
                merged["discovery_elite_band_source"] = "none"
            # WIM convergence MC band: ``convergence_ratio_to_wim`` is proportional to
            # ``required_density_for_match`` which inherits the same multiplicative scaling
            # in (D, i, V_flat) as ``theoretical_gap_coverage_ratio``. The relative MC band
            # therefore carries through one-to-one — we apply the same fractional widths to
            # the WIM convergence ratio when both numbers are finite.
            p16 = row.get("theoretical_gap_coverage_ratio_p16")
            p84 = row.get("theoretical_gap_coverage_ratio_p84")
            conv = merged.get("convergence_ratio_to_wim")
            try:
                p16_f = float(p16)
                p84_f = float(p84)
                conv_f = float(conv)
                if np.isfinite(p16_f) and np.isfinite(p84_f) and np.isfinite(conv_f) and ratio_f > 0:
                    merged["convergence_ratio_to_wim_p16"] = conv_f * (p16_f / ratio_f)
                    merged["convergence_ratio_to_wim_p84"] = conv_f * (p84_f / ratio_f)
                else:
                    merged["convergence_ratio_to_wim_p16"] = float("nan")
                    merged["convergence_ratio_to_wim_p84"] = float("nan")
            except (TypeError, ValueError):
                merged["convergence_ratio_to_wim_p16"] = float("nan")
                merged["convergence_ratio_to_wim_p84"] = float("nan")
            merged["wim_density_match_score"] = float(
                wim_density_match_score(
                    float(merged.get("required_density_for_match", float("nan"))),
                    float(merged.get("resolution_correction_factor", float("nan"))),
                )
            )
            merged["calibrated_convergence_ratio"] = row.get("calibrated_convergence_ratio")
            merged["validation_status"] = _validation_status_supreme(merged)
            rows.append(merged)
    rows.sort(key=_direct_evidence_sort_key)
    return rows


def _morphology_comparison_stats(
    supreme_rows: list[dict[str, Any]],
    morphology_lookup: dict[str, dict[str, Any]],
    morphology_provider: SparcMorphologyProvider,
) -> dict[str, dict[str, float]]:
    """
    Aggregate ``theoretical_gap_coverage_ratio × 100`` per morphology bucket.

    Returns a mapping ``{label: {"mean_pct": float, "count": int}}`` covering every morphology
    class that appears in the reliable supreme rows, plus a single ``Photometric LSB`` pseudo-bucket
    that gathers galaxies cut by ``SBdisk < PHOTOMETRIC_LSB_SB_THRESHOLD`` regardless of Hubble bin.
    Buckets with zero contributing galaxies are still surfaced (count=0) so the plot can decide
    whether to skip them; the ``mean_pct`` is then 0.0.
    """
    accumulator: dict[str, list[float]] = {}

    def _push(label: str, ratio: float) -> None:
        accumulator.setdefault(str(label), []).append(float(ratio) * 100.0)

    for row in supreme_rows or []:
        ratio = row.get("theoretical_gap_coverage_ratio")
        if ratio is None:
            continue
        try:
            ratio_f = float(ratio)
        except (TypeError, ValueError):
            continue
        if not np.isfinite(ratio_f):
            continue
        morphology_label = str(row.get("morphology_class") or "Unclassified")
        _push(morphology_label, ratio_f)
        if is_photometric_lsb(row, morphology_lookup, morphology_provider):
            _push(PHOTOMETRIC_LSB_LABEL, ratio_f)

    return {
        label: {
            "mean_pct": float(np.mean(values)) if values else 0.0,
            "count": int(len(values)),
        }
        for label, values in accumulator.items()
    }


def _build_median_dm_equivalence_grid(
    scenario_rows: dict[tuple[float, float], list[dict[str, Any]]] | None,
    gradient_scales: tuple[float, ...],
    ion_fractions: tuple[float, ...],
) -> np.ndarray:
    """
    For every (Δx, f_ion) cell, compute the median ``theoretical_gap_coverage_ratio``
    across the *reliable* galaxies in that scenario, expressed in percent.

    Cells with no reliable galaxies become NaN (rendered as a masked square in the heatmap).
    """
    grid = np.full((len(gradient_scales), len(ion_fractions)), np.nan, dtype=float)
    if scenario_rows is None:
        return grid
    for r_idx, scale in enumerate(gradient_scales):
        for c_idx, ion in enumerate(ion_fractions):
            cell_rows = scenario_rows.get((float(scale), float(ion)), [])
            ratios = []
            for row in cell_rows:
                if not bool(row.get("is_reliable", False)):
                    continue
                ratio = row.get("theoretical_gap_coverage_ratio")
                if ratio is None:
                    continue
                try:
                    val = float(ratio)
                except (TypeError, ValueError):
                    continue
                if not np.isfinite(val):
                    continue
                ratios.append(val * 100.0)
            if ratios:
                grid[r_idx, c_idx] = float(np.median(ratios))
    return grid


