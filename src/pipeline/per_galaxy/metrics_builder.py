"""Assemble per-scenario metric dictionaries used in CSV / JSON exports.

Audit-flag contract
-------------------
Every scenario row carries a list ``audit_flags`` enumerating the **exact** reasons it failed
the reliability checks. ``is_reliable`` is the derived signal ``len(audit_flags) == 0``. Codes:

    * ``fail_lat``                 — galactic latitude inside the foreground band.
    * ``fail_reliable_fraction``   — too few shells survived the per-point ``reliable_measurement``
                                     gate (driven by the magnetic-to-centripetal MAD test).
    * ``fail_density_clamped``     — too many radial shells fell back to the ρ_gas floor.
    * ``fail_stellar``             — too many shells dominated by V_disk/V_bulge (gas-mass tracer
                                     becomes invalid; magnetic pressure not applied there).
    * ``fail_skipped_physics``     — engine never produced ``a_B`` (missing Stokes / Δx / f_ion).
    * ``fail_ned_resolution``      — NED + SIMBAD could not resolve sky coordinates for the SPARC label.
    * ``warn_b_band``              — informational: median |B| outside the 0.1–100 µG band.
    * ``warn_no_pa`` / ``warn_pure_pressure`` / ``warn_btfr_outlier`` — engine / population diagnostics (non-blocking).

    * ``warn_unknown_morphology`` — SPARC morphology bucket missing or ``Unclassified`` (tension toy may be ambiguous).
    * ``warn_stokes_neighbor_imputed_center`` — only central Stokes amplitude was finite; neighbour mean copied from centre.

The audit is **per scenario**. The galaxy-level audit (in ``galaxy_processor``) inherits the
best-scenario flags and adds galaxy-only codes (``fail_inclination``, ``fail_calibration``).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from src.physics.diagnostics.statistics import edge_radial_mean, pearson_correlation_coefficient

from src.config.thresholds import (
    AUDIT_CODES,
    RELIABILITY_CUTS,
    calibrated_convergence_ratio_from_rho_and_k,
)


# Module-level alias preserved as the public API so callers / tests do not rebind to the
# centralised dataclass — the canonical token strings live in ``src.config.thresholds``.
AUDIT_FLAG_FAIL_LATITUDE = AUDIT_CODES.FAIL_LATITUDE
AUDIT_FLAG_FAIL_RELIABLE_FRACTION = AUDIT_CODES.FAIL_RELIABLE_FRACTION
AUDIT_FLAG_FAIL_DENSITY_CLAMPED = AUDIT_CODES.FAIL_DENSITY_CLAMPED
AUDIT_FLAG_FAIL_STELLAR = AUDIT_CODES.FAIL_STELLAR
AUDIT_FLAG_FAIL_SKIPPED_PHYSICS = AUDIT_CODES.FAIL_SKIPPED_PHYSICS
AUDIT_FLAG_FAIL_INCLINATION = AUDIT_CODES.FAIL_INCLINATION
AUDIT_FLAG_FAIL_CALIBRATION = AUDIT_CODES.FAIL_CALIBRATION
AUDIT_FLAG_FAIL_NED_RESOLUTION = AUDIT_CODES.FAIL_NED_RESOLUTION
AUDIT_FLAG_WARN_B_BAND = AUDIT_CODES.WARN_MAGNETIC_FIELD_BAND
AUDIT_FLAG_WARN_NO_POSITION_ANGLE = AUDIT_CODES.WARN_NO_POSITION_ANGLE
AUDIT_FLAG_WARN_PURE_PRESSURE_GEOMETRY = AUDIT_CODES.WARN_PURE_PRESSURE_GEOMETRY
AUDIT_FLAG_WARN_BTFR_OUTLIER = AUDIT_CODES.WARN_BTFR_OUTLIER
AUDIT_FLAG_WARN_UNKNOWN_MORPHOLOGY = AUDIT_CODES.WARN_UNKNOWN_MORPHOLOGY
AUDIT_FLAG_WARN_STOKES_NEIGHBOR_IMPUTED = AUDIT_CODES.WARN_STOKES_NEIGHBOR_IMPUTED_FROM_CENTER

ALL_AUDIT_FLAG_CODES = (
    AUDIT_FLAG_FAIL_LATITUDE,
    AUDIT_FLAG_FAIL_RELIABLE_FRACTION,
    AUDIT_FLAG_FAIL_DENSITY_CLAMPED,
    AUDIT_FLAG_FAIL_STELLAR,
    AUDIT_FLAG_FAIL_SKIPPED_PHYSICS,
    AUDIT_FLAG_FAIL_INCLINATION,
    AUDIT_FLAG_FAIL_CALIBRATION,
    AUDIT_FLAG_FAIL_NED_RESOLUTION,
    AUDIT_FLAG_WARN_B_BAND,
    AUDIT_FLAG_WARN_NO_POSITION_ANGLE,
    AUDIT_FLAG_WARN_PURE_PRESSURE_GEOMETRY,
    AUDIT_FLAG_WARN_BTFR_OUTLIER,
    AUDIT_FLAG_WARN_UNKNOWN_MORPHOLOGY,
    AUDIT_FLAG_WARN_STOKES_NEIGHBOR_IMPUTED,
)


def _is_blocking_flag(code: str) -> bool:
    """``warn_*`` codes are informational-only; everything else blocks reliability."""
    return not str(code).startswith("warn_")


def blocking_flags(audit_flags: list[str] | tuple[str, ...]) -> list[str]:
    """
    Return the blocking subset of ``audit_flags``.

    Centralised so the convention "everything not prefixed by ``warn_`` blocks reliability"
    has a single, testable implementation. Callers (``metrics_for_successful_physics``,
    ``galaxy_processor``) reuse this helper instead of re-implementing the prefix check.
    """
    return [code for code in audit_flags if _is_blocking_flag(code)]


def aggregate_galaxy_audit_flags(
    *,
    scenario_audit_flags: list[str] | tuple[str, ...] | None,
    inclination_quality_passed: bool,
    calibration_status: str,
    suboptimal_calibration_status: str,
    ned_unresolved: bool = False,
) -> tuple[list[str], list[str], bool]:
    """
    Compose the galaxy-level audit-flag list from the best scenario plus galaxy-only gates.

    The function inherits the scenario-level codes verbatim and **appends** the galaxy-only
    failures that scenario rows cannot see:

    * :data:`AUDIT_FLAG_FAIL_INCLINATION` when ``inclination_quality_passed`` is ``False``
      (face-on systems destabilise the :math:`V_{\\mathrm{obs}}/\\sin i` propagation).
    * :data:`AUDIT_FLAG_FAIL_CALIBRATION` when ``calibration_status`` equals the configured
      ``suboptimal_calibration_status`` (the Stokes amplitude reference fell back to the
      emergency value and any inferred :math:`B` is therefore non-measured).
    * :data:`AUDIT_FLAG_FAIL_NED_RESOLUTION` when ``ned_unresolved`` is ``True`` (no NED/SIMBAD
      coordinates for the catalogue label).

    Parameters
    ----------
    scenario_audit_flags
        Best-scenario ``audit_flags`` list (may be ``None`` or empty).
    inclination_quality_passed
        ``False`` if the galaxy inclination is below the configured floor.
    calibration_status
        The Stokes calibration status emitted for this galaxy.
    suboptimal_calibration_status
        The string token that signals a fallback (typically
        :data:`src.config.thresholds.CALIBRATION_STATUS_SUBOPTIMAL`).
    ned_unresolved
        ``True`` when NED and SIMBAD both failed to resolve sky coordinates for this SPARC name.

    Returns
    -------
    tuple[list[str], list[str], bool]
        ``(galaxy_audit_flags, blocking_subset, is_reliable)``.
    """
    galaxy_flags: list[str] = list(scenario_audit_flags or [])
    if not inclination_quality_passed and AUDIT_FLAG_FAIL_INCLINATION not in galaxy_flags:
        galaxy_flags.append(AUDIT_FLAG_FAIL_INCLINATION)
    if calibration_status == suboptimal_calibration_status and AUDIT_FLAG_FAIL_CALIBRATION not in galaxy_flags:
        galaxy_flags.append(AUDIT_FLAG_FAIL_CALIBRATION)
    if ned_unresolved and AUDIT_FLAG_FAIL_NED_RESOLUTION not in galaxy_flags:
        galaxy_flags.append(AUDIT_FLAG_FAIL_NED_RESOLUTION)
    blocking = blocking_flags(galaxy_flags)
    return galaxy_flags, blocking, len(blocking) == 0


def metrics_for_successful_physics(
    *,
    galaxy_name: str,
    scenario_name: str,
    gradient_scale_kpc: float,
    ionization_fraction: float,
    scenario_df: Any,
    mean_acceleration_gap_edge: float | None,
    galactic_latitude_degrees: float | None,
    morphology_class: str,
    engine_audit_flags: list[str] | tuple[str, ...] | None = None,
    required_density_for_match_median: float | None = None,
    resolution_correction_factor_k: float | None = None,
) -> dict[str, Any]:
    reliable_measurement = scenario_df["reliable_measurement"].to_numpy(dtype=bool)
    reliable_gap = scenario_df["acceleration_gap"].to_numpy(dtype=float)[reliable_measurement]
    reliable_mag = scenario_df["magnetic_pressure_acceleration"].to_numpy(dtype=float)[reliable_measurement]
    pearson = pearson_correlation_coefficient(reliable_gap, reliable_mag)
    mean_mag_edge = edge_radial_mean(scenario_df["magnetic_pressure_acceleration"])
    theoretical_gap_coverage_ratio = (
        abs(mean_mag_edge) / abs(mean_acceleration_gap_edge)
        if (mean_acceleration_gap_edge not in (None, 0) and mean_mag_edge is not None)
        else np.nan
    )
    dark_matter_equivalent_percentage = (
        theoretical_gap_coverage_ratio * 100.0 if np.isfinite(theoretical_gap_coverage_ratio) else np.nan
    )

    # Full-Lorentz parallel coverage ratio (pressure gradient + toy magnetic-tension term).
    # Reported alongside the headline ``theoretical_gap_coverage_ratio`` so readers can see how
    # much the toy spiral geometry inflates / depresses the coverage.
    if "accel_full_lorentz_m_s2" in scenario_df.columns:
        mean_full_lorentz_edge = edge_radial_mean(scenario_df["accel_full_lorentz_m_s2"])
        theoretical_gap_coverage_ratio_full_lorentz = (
            abs(mean_full_lorentz_edge) / abs(mean_acceleration_gap_edge)
            if (mean_acceleration_gap_edge not in (None, 0) and mean_full_lorentz_edge is not None)
            else np.nan
        )
    else:
        theoretical_gap_coverage_ratio_full_lorentz = np.nan
    # Scenarios that skipped the physics path emit an all-NaN column; ``np.nanmedian``
    # warns on an empty slice in that case, so we silence the warning at the call site
    # (the ``None`` fallback already captures the empty-slice semantics).
    with np.errstate(invalid="ignore"):
        if "accel_pressure_only_m_s2" in scenario_df.columns:
            arr = scenario_df["accel_pressure_only_m_s2"].to_numpy(dtype=float)
            if np.any(np.isfinite(arr)):
                accel_pressure_only_median = float(np.nanmedian(arr))
                if not np.isfinite(accel_pressure_only_median):
                    accel_pressure_only_median = None
            else:
                accel_pressure_only_median = None
        else:
            accel_pressure_only_median = None
        if "accel_full_lorentz_m_s2" in scenario_df.columns:
            arr = scenario_df["accel_full_lorentz_m_s2"].to_numpy(dtype=float)
            if np.any(np.isfinite(arr)):
                accel_full_lorentz_median = float(np.nanmedian(arr))
                if not np.isfinite(accel_full_lorentz_median):
                    accel_full_lorentz_median = None
            else:
                accel_full_lorentz_median = None
        else:
            accel_full_lorentz_median = None

    # Density-clamp audit: count radial shells whose ρ_gas was synthesized from the emergency floor.
    # If too many points were clamped, the scenario forfeits reliability regardless of other gates.
    if "density_floor_clamped" in scenario_df.columns:
        density_clamped_arr = scenario_df["density_floor_clamped"].to_numpy(dtype=bool)
    else:
        density_clamped_arr = np.zeros(len(scenario_df), dtype=bool)
    n_points_density_clamped = int(np.count_nonzero(density_clamped_arr))
    n_points_total = int(density_clamped_arr.size)
    density_clamped_fraction = (
        float(n_points_density_clamped) / float(n_points_total) if n_points_total > 0 else 0.0
    )
    density_audit_passed = (
        density_clamped_fraction <= RELIABILITY_CUTS.MAX_DENSITY_CLAMPED_FRACTION_FOR_RELIABILITY
    )

    # Stellar-dominance audit: if too many shells are dominated by V_disk + V_bulge, the
    # ``M_enc ≈ V_gas² r / G`` density estimator is invalid and the scenario can not be trusted.
    if "stellar_dominated" in scenario_df.columns:
        stellar_dominated_arr = scenario_df["stellar_dominated"].to_numpy(dtype=bool)
    else:
        stellar_dominated_arr = np.zeros(len(scenario_df), dtype=bool)
    n_points_stellar_dominated = int(np.count_nonzero(stellar_dominated_arr))
    stellar_dominated_fraction = (
        float(n_points_stellar_dominated) / float(n_points_total) if n_points_total > 0 else 0.0
    )
    stellar_dominance_audit_passed = (
        stellar_dominated_fraction
        <= RELIABILITY_CUTS.MAX_STELLAR_DOMINANCE_FRACTION_FOR_RELIABILITY
    )

    # Per-radius RAR residuals (if the engine populated them) — surface the median |residual|
    # in dex for the scenario row. Diagnostic only; does NOT participate in is_reliable.
    if "rar_residual_log10" in scenario_df.columns:
        rar_arr = scenario_df["rar_residual_log10"].to_numpy(dtype=float)
        rar_arr = rar_arr[np.isfinite(rar_arr)]
        rar_median_abs_residual_dex = (
            float(np.median(np.abs(rar_arr))) if rar_arr.size > 0 else None
        )
    else:
        rar_median_abs_residual_dex = None

    # Sanity diagnostics from the scenario dataframe (median μG, median plasma β).
    # ``magnetic_field_microgauss`` is the per-shell B field expressed in microgauss
    # (1 μG = 1e-10 T) — the unit galactic-disk literature uses. The flag
    # ``magnetic_field_within_plausibility_band`` reports whether the median falls
    # inside the observational band [0.1, 100] μG configured in RELIABILITY_CUTS.
    magnetic_field_microgauss_median = (
        float(np.nanmedian(scenario_df["magnetic_field_microgauss"].to_numpy(dtype=float)))
        if "magnetic_field_microgauss" in scenario_df.columns
        else None
    )
    if magnetic_field_microgauss_median is not None and not np.isfinite(magnetic_field_microgauss_median):
        magnetic_field_microgauss_median = None
    magnetic_field_within_band = (
        magnetic_field_microgauss_median is not None
        and RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MIN_MICROGAUSS
        <= magnetic_field_microgauss_median
        <= RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MAX_MICROGAUSS
    )
    # Plasma β = 2·μ₀·P_thermal / B². β < 1 ⇒ magnetic pressure dominates the dynamics
    # (the central MHD-pressure-gradient hypothesis is locally consistent). β ≫ 1 ⇒
    # thermal pressure dwarfs the magnetic stress and ``a_B`` should be interpreted as
    # an upper bound only. We surface the median, the flag, and the fraction of shells
    # below unity so reports can quantify how much of the radial range is in the
    # magnetically-dominated regime.
    if "plasma_beta_thermal" in scenario_df.columns:
        beta_arr = scenario_df["plasma_beta_thermal"].to_numpy(dtype=float)
        beta_finite = beta_arr[np.isfinite(beta_arr)]
        plasma_beta_median = float(np.median(beta_finite)) if beta_finite.size > 0 else None
        plasma_beta_below_unity_fraction = (
            float(np.mean(beta_finite < 1.0)) if beta_finite.size > 0 else None
        )
    else:
        plasma_beta_median = None
        plasma_beta_below_unity_fraction = None
    magnetic_dominance_flag = (
        plasma_beta_median is not None and plasma_beta_median < 1.0
    )

    # Build the structured audit-flag list — every code corresponds to a SPECIFIC, traceable
    # reason this scenario was rejected. Empty list ⇒ reliable. Galaxy-level codes (inclination,
    # calibration) are appended later in ``galaxy_processor`` since they are not visible here.
    reliable_fraction = (
        float(np.mean(reliable_measurement)) if reliable_measurement.size else 0.0
    )

    audit_flags: list[str] = []
    if (
        galactic_latitude_degrees is None
        or abs(galactic_latitude_degrees) < RELIABILITY_CUTS.MIN_GALACTIC_LATITUDE_DEG
    ):
        audit_flags.append(AUDIT_FLAG_FAIL_LATITUDE)
    if reliable_fraction < RELIABILITY_CUTS.MIN_RELIABLE_MEASUREMENT_FRACTION:
        audit_flags.append(AUDIT_FLAG_FAIL_RELIABLE_FRACTION)
    if not density_audit_passed:
        audit_flags.append(AUDIT_FLAG_FAIL_DENSITY_CLAMPED)
    if not stellar_dominance_audit_passed:
        audit_flags.append(AUDIT_FLAG_FAIL_STELLAR)
    if (
        magnetic_field_microgauss_median is not None
        and not magnetic_field_within_band
    ):
        audit_flags.append(AUDIT_FLAG_WARN_B_BAND)
    # Forward warning-only engine flags (``warn_no_pa`` / ``warn_pure_pressure``) so the scenario
    # row's ``audit_flags`` advertises the missing supplementary information without blocking
    # reliability — both codes are ``warn_*`` per the engine's contract.
    for engine_code in engine_audit_flags or ():
        if engine_code and engine_code not in audit_flags:
            audit_flags.append(str(engine_code))

    blocking = blocking_flags(audit_flags)
    is_reliable = len(blocking) == 0
    calibrated_convergence_ratio = float("nan")
    if (
        required_density_for_match_median is not None
        and resolution_correction_factor_k is not None
        and np.isfinite(float(required_density_for_match_median))
        and np.isfinite(float(resolution_correction_factor_k))
    ):
        calibrated_convergence_ratio = float(
            calibrated_convergence_ratio_from_rho_and_k(
                float(required_density_for_match_median),
                float(resolution_correction_factor_k),
            )
        )
    explained_fraction = (
        float(np.mean(np.abs(reliable_mag) >= RELIABILITY_CUTS.EXPLAINED_GAP_HALF_FACTOR * np.abs(reliable_gap)))
        if reliable_gap.size > 0
        else 0.0
    )
    return {
        "galaxy_name": galaxy_name,
        "scenario_name": scenario_name,
        "gradient_scale_kpc": gradient_scale_kpc,
        "ionization_fraction": ionization_fraction,
        "mean_acceleration_gap_edge": mean_acceleration_gap_edge,
        "mean_magnetic_pressure_acceleration_edge": mean_mag_edge,
        "theoretical_gap_coverage_ratio": theoretical_gap_coverage_ratio,
        "theoretical_gap_coverage_ratio_full_lorentz": theoretical_gap_coverage_ratio_full_lorentz,
        "accel_pressure_only_m_s2": accel_pressure_only_median,
        "accel_full_lorentz_m_s2": accel_full_lorentz_median,
        "dark_matter_equivalent_percentage": dark_matter_equivalent_percentage,
        "pearson_correlation": pearson,
        "is_reliable": is_reliable,
        "audit_flags": audit_flags,
        "reliability_failure_codes": list(blocking),
        "reliable_measurement_fraction": reliable_fraction,
        "galactic_latitude_degrees": galactic_latitude_degrees,
        "morphology_class": morphology_class,
        "explained_fraction": explained_fraction,
        "n_points_density_clamped": n_points_density_clamped,
        "n_points_total": n_points_total,
        "density_clamped_fraction": density_clamped_fraction,
        "density_audit_passed": density_audit_passed,
        "n_points_stellar_dominated": n_points_stellar_dominated,
        "stellar_dominated_fraction": stellar_dominated_fraction,
        "stellar_dominance_audit_passed": stellar_dominance_audit_passed,
        "rar_median_abs_residual_dex": rar_median_abs_residual_dex,
        "magnetic_field_microgauss_median": magnetic_field_microgauss_median,
        "magnetic_field_within_plausibility_band": bool(magnetic_field_within_band),
        "plasma_beta_thermal_median": plasma_beta_median,
        "plasma_beta_below_unity_fraction": plasma_beta_below_unity_fraction,
        "magnetic_dominance_flag": bool(magnetic_dominance_flag),
        "required_density_for_match_median": (
            float(required_density_for_match_median)
            if required_density_for_match_median is not None
            and np.isfinite(float(required_density_for_match_median))
            else None
        ),
        "resolution_correction_factor": (
            float(resolution_correction_factor_k)
            if resolution_correction_factor_k is not None
            and np.isfinite(float(resolution_correction_factor_k))
            else float("nan")
        ),
        "calibrated_convergence_ratio": calibrated_convergence_ratio,
        "_increment_success": bool(
            is_reliable
            and np.isfinite(theoretical_gap_coverage_ratio)
            and RELIABILITY_CUTS.SCENARIO_SUCCESS_RATIO_MIN
            <= theoretical_gap_coverage_ratio
            <= RELIABILITY_CUTS.SCENARIO_SUCCESS_RATIO_MAX
        ),
    }


def metrics_for_skipped_physics(
    *,
    galaxy_name: str,
    scenario_name: str,
    gradient_scale_kpc: float,
    ionization_fraction: float,
    mean_acceleration_gap_edge: float | None,
    galactic_latitude_degrees: float | None,
    morphology_class: str,
) -> dict[str, Any]:
    audit_flags: list[str] = [AUDIT_FLAG_FAIL_SKIPPED_PHYSICS]
    if (
        galactic_latitude_degrees is None
        or abs(galactic_latitude_degrees) < RELIABILITY_CUTS.MIN_GALACTIC_LATITUDE_DEG
    ):
        audit_flags.append(AUDIT_FLAG_FAIL_LATITUDE)
    return {
        "galaxy_name": galaxy_name,
        "scenario_name": scenario_name,
        "gradient_scale_kpc": gradient_scale_kpc,
        "ionization_fraction": ionization_fraction,
        "mean_acceleration_gap_edge": mean_acceleration_gap_edge,
        "mean_magnetic_pressure_acceleration_edge": None,
        "theoretical_gap_coverage_ratio": np.nan,
        "theoretical_gap_coverage_ratio_full_lorentz": np.nan,
        "accel_pressure_only_m_s2": None,
        "accel_full_lorentz_m_s2": None,
        "dark_matter_equivalent_percentage": np.nan,
        "pearson_correlation": None,
        "is_reliable": False,
        "audit_flags": audit_flags,
        "reliability_failure_codes": [c for c in audit_flags if _is_blocking_flag(c)],
        "reliable_measurement_fraction": 0.0,
        "galactic_latitude_degrees": galactic_latitude_degrees,
        "morphology_class": morphology_class,
        "explained_fraction": 0.0,
        "n_points_density_clamped": 0,
        "n_points_total": 0,
        "density_clamped_fraction": 0.0,
        "density_audit_passed": True,
        "n_points_stellar_dominated": 0,
        "stellar_dominated_fraction": 0.0,
        "stellar_dominance_audit_passed": True,
        "rar_median_abs_residual_dex": None,
        "magnetic_field_microgauss_median": None,
        "magnetic_field_within_plausibility_band": False,
        "plasma_beta_thermal_median": None,
        "plasma_beta_below_unity_fraction": None,
        "magnetic_dominance_flag": False,
        "required_density_for_match_median": None,
        "resolution_correction_factor": float("nan"),
        "calibrated_convergence_ratio": float("nan"),
        "_increment_success": False,
    }
