"""
Centralized physical constants, astrophysical defaults, Planck calibration policy, and reliability cuts.

Values are grouped for clarity; ``UNIVERSAL_CONSTANTS`` pulls :math:`\\mu_0` and :math:`G` from
``astropy.constants`` (CODATA) as plain high-precision ``float`` in SI units.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from astropy import constants as const
from astropy import units as astropy_u


@dataclass(frozen=True, slots=True)
class UniversalConstants:
    """SI fundamental / CODATA constants used in kinematic and magnetic formulas.

    Sourced from :mod:`astropy.constants` (CODATA) at import time and frozen so the
    pipeline never accidentally rebinds them.
    """

    # Vacuum magnetic permeability (T·m/A ≡ N/A²).
    MU0_SI: float
    # Newtonian constant of gravitation (m³·kg⁻¹·s⁻²).
    G_SI: float
    # Microgauss → Tesla conversion. Exact by definition: ``1 μG = 1·10⁻¹⁰ T``.
    # Centralised so every diagnostic / reporting layer reads the same value.
    TESLA_PER_MICROGAUSS_SI: float


def _build_universal_constants() -> UniversalConstants:
    try:
        mu0_q = const.mu0
        mu0_si = float(mu0_q.to_value(astropy_u.T * astropy_u.m / astropy_u.A))
    except AttributeError:
        # Astropy without ``mu0`` — SI exact definition 4π×10⁻⁷ H/m.
        mu0_si = float(np.float64(4.0 * np.pi * 1e-7))
    g_q = const.G
    g_si = float(g_q.to_value(astropy_u.m**3 / (astropy_u.kg * astropy_u.s**2)))
    return UniversalConstants(
        MU0_SI=mu0_si,
        G_SI=g_si,
        TESLA_PER_MICROGAUSS_SI=float(np.float64(1.0e-10)),
    )


UNIVERSAL_CONSTANTS: Final[UniversalConstants] = _build_universal_constants()


@dataclass(frozen=True, slots=True)
class AstrophysicalDefaults:
    """Fixed numerical defaults for disk geometry, grid search, and CMB anchoring scale."""

    # CMB energy density (J/m³) — dimensionful anchor for Stokes → B (see EG_MODEL_THEORY.md §3.3); not a field model.
    CMB_ENERGY_DENSITY_J_M3: float
    # Floor on reconstructed gas density (kg/m³) when kinematic inversion yields non-positive shells.
    MIN_GAS_MASS_DENSITY_KG_M3: float
    # Cylindrical shell thickness h(R_max) = clip(f × R_max, h_min, h_max) in kpc.
    SHELL_THICKNESS_FRACTION_OF_RMAX: float
    SHELL_THICKNESS_CLIP_MIN_KPC: float
    SHELL_THICKNESS_CLIP_MAX_KPC: float
    # Scenario grid (gradient scale in kpc × effective ionized fraction).
    GRADIENT_SCALES_KPC: tuple[float, ...]
    IONIZATION_FRACTIONS: tuple[float, ...]
    # Radial Acceleration Relation (McGaugh+2016) acceleration scale a0 in m/s².
    RAR_A0_M_S2: float
    # Plasma sanity assumptions for the β diagnostic: representative warm-ionized-medium
    # temperature (Kelvin) and proton mass for n = ρ/m_p.
    PLASMA_THERMAL_TEMPERATURE_K: float
    PROTON_MASS_KG: float
    BOLTZMANN_J_PER_K: float
    # Linear-propagation factor in δa_c / a_c = K · |δV / V|. The factor is 2 for the
    # second-power dependence ``a_c = V²/r``; centralised so unit tests can assert the
    # algebraic origin, not a hard-coded ``2.0`` literal scattered around the kinematics
    # module.
    INCLINATION_SYSTEMATIC_PROPAGATION_FACTOR: float
    # Toy logarithmic-spiral pitch angle (deg) used by ``TensionTermEstimator`` to derive the
    # field-line curvature radius ``R_c = r / sin²(p)``. 12° matches the median pitch reported
    # for SAB/SB galaxies in Yuan & Kuo 1996 and Diaz-Garcia+2019.
    SPIRAL_PITCH_DEG: float
    # Turbulent-equipartition magnetic energy-density anchor (J/m³). Computed as
    # ``0.5 · ρ_gas_typical · σ_v²`` with ``σ_v = 10 km/s`` (warm-neutral-medium dispersion)
    # and ``ρ_gas_typical = 1e-24 kg/m³`` (mean disk gas density). This is an alternative
    # to ``CMB_ENERGY_DENSITY_J_M3`` for the anchor-sensitivity sweep — both anchors are
    # dimensionful only (no physical claim of origin).
    TURBULENT_EQUIPARTITION_J_M3: float
    # McGaugh+2012 Baryonic Tully–Fisher relation ridge:
    # ``log10(M_bar / M_sun) = SLOPE · log10(V_flat / (km/s)) + INTERCEPT``.
    BTFR_RIDGE_SLOPE: float
    BTFR_RIDGE_INTERCEPT: float
    # Galaxies sitting farther than this distance (dex) from the BTFR ridge are flagged
    # with ``warn_btfr_outlier``. The threshold is informational; it does not block.
    BTFR_OUTLIER_THRESHOLD_DEX: float
    # Spitzer 3.6 µm mass-to-light ratio (solar units) used to convert ``L36_1e9solLum``
    # to ``M_star``. 0.5 is the canonical value adopted by SPARC (Lelli+2016, Schombert+2014).
    STELLAR_MASS_TO_LIGHT_3P6UM: float
    # Helium correction factor applied to the catalogued HI mass to obtain the cold-gas
    # contribution to the BTFR baryonic mass. ``M_gas = 1.33 · M_HI`` (helium fraction 0.25).
    HELIUM_CORRECTION_FACTOR_HI_TO_GAS: float
    # Reference warm-ionised-medium density (kg/m³) used by Discovery convergence and the
    # scientific-evidence ``wim_density_match_score`` (ρ_eff relative to this anchor).
    WIM_REFERENCE_DENSITY_KG_M3: float


ASTROPHYSICAL_DEFAULTS: Final[AstrophysicalDefaults] = AstrophysicalDefaults(
    CMB_ENERGY_DENSITY_J_M3=float(np.float64(4.17e-14)),
    MIN_GAS_MASS_DENSITY_KG_M3=float(np.float64(1e-25)),
    SHELL_THICKNESS_FRACTION_OF_RMAX=float(np.float64(0.05)),
    SHELL_THICKNESS_CLIP_MIN_KPC=float(np.float64(0.2)),
    SHELL_THICKNESS_CLIP_MAX_KPC=float(np.float64(2.0)),
    GRADIENT_SCALES_KPC=(float(np.float64(0.05)), float(np.float64(0.1)), float(np.float64(0.5))),
    IONIZATION_FRACTIONS=(float(np.float64(0.001)), float(np.float64(0.01)), float(np.float64(0.1))),
    RAR_A0_M_S2=float(np.float64(1.2e-10)),
    PLASMA_THERMAL_TEMPERATURE_K=float(np.float64(1.0e4)),
    PROTON_MASS_KG=float(const.m_p.to_value(astropy_u.kg)),
    BOLTZMANN_J_PER_K=float(const.k_B.to_value(astropy_u.J / astropy_u.K)),
    INCLINATION_SYSTEMATIC_PROPAGATION_FACTOR=float(np.float64(2.0)),
    SPIRAL_PITCH_DEG=float(np.float64(12.0)),
    # 0.5 × 1e-24 kg/m³ × (10 000 m/s)² = 5e-17 J/m³.
    TURBULENT_EQUIPARTITION_J_M3=float(np.float64(5.0e-17)),
    BTFR_RIDGE_SLOPE=float(np.float64(4.0)),
    BTFR_RIDGE_INTERCEPT=float(np.float64(1.99)),
    BTFR_OUTLIER_THRESHOLD_DEX=float(np.float64(0.5)),
    STELLAR_MASS_TO_LIGHT_3P6UM=float(np.float64(0.5)),
    HELIUM_CORRECTION_FACTOR_HI_TO_GAS=float(np.float64(1.33)),
    WIM_REFERENCE_DENSITY_KG_M3=float(np.float64(1.0e-27)),
)

# Minor-axis sky offsets are stretched by 1/cos(i) onto the disk plane; floor avoids LSQ blow-up as i → 90°.
DISK_PLANE_COS_INCLINATION_NUMERICAL_FLOOR: Final[float] = float(np.float64(1.0e-3))


@dataclass(frozen=True, slots=True)
class UncertaintyConstants:
    """Monte-Carlo propagation defaults for the gap-coverage uncertainty band."""

    # Number of MC samples drawn per galaxy when propagating (D, i, V_flat) systematics.
    # 300 keeps the runtime per galaxy well under 10 ms with vectorised numpy ops while
    # leaving enough samples to estimate the p16/p84 quantiles with ~3% precision.
    MC_DEFAULT_SAMPLES: int
    # RNG seed for the Monte-Carlo step. Determinism is mandatory so reports are byte-stable
    # across re-runs on the same input data.
    MC_DEFAULT_SEED: int


UNCERTAINTY: Final[UncertaintyConstants] = UncertaintyConstants(
    MC_DEFAULT_SAMPLES=300,
    MC_DEFAULT_SEED=42,
)


@dataclass(frozen=True, slots=True)
class PlanckCalibration:
    """
    Planck map calibration: global polarization amplitude reference from the map itself.

    ``stokes_amp_reference`` is set to the **median** of :math:`\\sqrt{Q^2+U^2}` over a strided subsample
    of all HEALPix rows (bounded by ``STOKES_MEDIAN_MAX_SAMPLES`` for RAM). If that is non-finite or
    non-positive, ``STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY`` is used and the calibration status is
    flagged as ``SUBOPTIMAL_FALLBACK`` downstream (see ``PlanckProvider`` and ``galaxy_processor``).
    """

    # Emergency-only constant: used **only** when the dynamic median cannot be computed.
    # Reaching this value triggers a SUBOPTIMAL_FALLBACK flag on the calibration metadata.
    STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY: float
    STOKES_MEDIAN_MAX_SAMPLES: int
    # Minimum positive median to accept (guards against all-zero maps).
    STOKES_MEDIAN_MIN_POSITIVE: float
    # Planck GNILC thermal-dust map effective beam FWHM (arcmin). Values reported by the
    # Planck Collaboration (2016) for the unires GNILC product at 353 GHz / 545 GHz lie
    # in the 5'-band; the constant feeds the beam-smearing diagnostic that compares the
    # physical projection of the beam (``D_beam = distance · θ_beam``) with the assumed
    # disk-plane gradient scale used by the engine. Reading pointer:
    # :mod:`src.physics.diagnostics.resolution`.
    BEAM_FWHM_ARCMIN: float


PLANCK_CALIBRATION: Final[PlanckCalibration] = PlanckCalibration(
    STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY=float(np.float64(1e-4)),
    STOKES_MEDIAN_MAX_SAMPLES=5_000_000,
    STOKES_MEDIAN_MIN_POSITIVE=float(np.float64(1e-30)),
    BEAM_FWHM_ARCMIN=float(np.float64(5.0)),
)


# Provenance / calibration status enum (string constants, not Enum, to keep JSON serialization trivial).
CALIBRATION_REFERENCE_METHOD_DYNAMIC: Final[str] = "dynamic_map_median"
CALIBRATION_REFERENCE_METHOD_FALLBACK: Final[str] = "default_fallback_warning"
CALIBRATION_STATUS_OPTIMAL: Final[str] = "OPTIMAL"
CALIBRATION_STATUS_SUBOPTIMAL: Final[str] = "SUBOPTIMAL_FALLBACK"


@dataclass(frozen=True, slots=True)
class AuditCodes:
    """
    String tokens emitted in ``audit_flags`` / ``reliability_failure_codes`` JSON arrays.

    Centralising these lets tests and downstream consumers reference the same tokens
    without grep-able literals scattered across the metrics builder. They are *strings*
    (not :class:`enum.Enum`) so the per-galaxy JSON stays trivially serialisable.
    """

    FAIL_LATITUDE: str
    FAIL_RELIABLE_FRACTION: str
    FAIL_DENSITY_CLAMPED: str
    FAIL_STELLAR: str
    FAIL_SKIPPED_PHYSICS: str
    FAIL_INCLINATION: str
    FAIL_CALIBRATION: str
    FAIL_NED_RESOLUTION: str
    WARN_MAGNETIC_FIELD_BAND: str
    # Phase-1 refinements: warning-level codes that do NOT block reliability.
    WARN_NO_POSITION_ANGLE: str
    # Irregular / BCD hosts: pressure-only magnetic geometry (toy spiral tension disabled).
    WARN_PURE_PRESSURE_GEOMETRY: str
    WARN_BTFR_OUTLIER: str
    # Morphology bucket missing from the SPARC morphology JSON (treated as unknown, not Irregular).
    WARN_UNKNOWN_MORPHOLOGY: str
    # Only the central Stokes amplitude was finite; neighbour mean was imputed equal to centre.
    WARN_STOKES_NEIGHBOR_IMPUTED_FROM_CENTER: str


AUDIT_CODES: Final[AuditCodes] = AuditCodes(
    FAIL_LATITUDE="fail_lat",
    FAIL_RELIABLE_FRACTION="fail_reliable_fraction",
    FAIL_DENSITY_CLAMPED="fail_density_clamped",
    FAIL_STELLAR="fail_stellar",
    FAIL_SKIPPED_PHYSICS="fail_skipped_physics",
    FAIL_INCLINATION="fail_inclination",
    FAIL_CALIBRATION="fail_calibration",
    FAIL_NED_RESOLUTION="fail_ned_resolution",
    WARN_MAGNETIC_FIELD_BAND="warn_b_band",
    WARN_NO_POSITION_ANGLE="warn_no_pa",
    WARN_PURE_PRESSURE_GEOMETRY="warn_pure_pressure",
    WARN_BTFR_OUTLIER="warn_btfr_outlier",
    WARN_UNKNOWN_MORPHOLOGY="warn_unknown_morphology",
    WARN_STOKES_NEIGHBOR_IMPUTED_FROM_CENTER="warn_stokes_neighbor_imputed_center",
)

# One-line English prose for CSV ``note`` / ``reliability_note`` when a blocking ``fail_*`` gate fires.
AUDIT_BLOCKING_NOTE_EN: Final[dict[str, str]] = {
    "fail_lat": "Galactic latitude inside the Milky-Way foreground band; polarized dust can bias the Stokes-based field.",
    "fail_reliable_fraction": "Too few shells retained reliable magnetometry after robust magnetic-pressure outlier rejection.",
    "fail_density_clamped": "Gas density hit the numerical floor on too many shells, so the inferred ion density is model-dominated.",
    "fail_stellar": "Stellar gravity dominates baryonic acceleration on too many shells for a gas-only density estimate.",
    "fail_skipped_physics": "Magnetic pressure gradient could not be evaluated (missing physics inputs or numerical failure).",
    "fail_inclination": "Inclination below the pipeline floor; V_obs/sin(i) deprojection is unstable.",
    "fail_calibration": "Stokes Q/U calibration used the suboptimal or fallback reference level.",
    "fail_ned_resolution": "NED/SIMBAD did not resolve coordinates for this label; sky position is unknown.",
}

# One-line English prose for CSV ``note`` / ``reliability_note`` when only non-blocking ``warn_*``
# audit codes are present (these document provenance and interpretation limits).
AUDIT_WARN_NOTE_EN: Final[dict[str, str]] = {
    "warn_b_band": "Median inferred |B| lies outside the nominal 0.1–100 microgauss plausibility band; treat |B| as order-of-magnitude only.",
    "warn_no_pa": "No curated position angle in the SPARC supplement; magnetic gradient uses the axisymmetric sky-plane estimate with cos(i) correction instead of spiral-arm geometry.",
    "warn_pure_pressure": "Irregular or BCD morphology: the toy logarithmic-spiral magnetic-tension term is disabled (pressure channel only).",
    "warn_btfr_outlier": "Baryonic Tully–Fisher residual exceeds 0.5 dex from the McGaugh+2012 ridge; kinematic mass scaling may be atypical for this system.",
    "warn_unknown_morphology": "Morphology bucket missing or Unclassified; spiral-tension assumptions may not apply.",
    "warn_stokes_neighbor_imputed_center": "Only the central Stokes amplitude was finite on the neighbour ring; neighbour mean was copied from the centre (noisier gradient estimate).",
}


@dataclass(frozen=True, slots=True)
class ReliabilityCuts:
    """Quality gates, evidence windows, and narrative highlight thresholds (decision policy)."""

    # Legacy hard ceiling on |a_B|/|a_c| used only when too few shells remain to build a robust
    # per-galaxy log-ratio distribution (see ``MAG_PRESSURE_RATIO_ROBUST_STATS_MIN_SHELLS``).
    MAX_MAGNETIC_PRESSURE_TO_CENTRIPETAL_ACCEL_FACTOR: float
    # Robust outlier gate on R = |a_B|/|a_c|: keep shells with log10(R) <= median + Nσ·σ_robust,
    # where σ_robust is MAD×1.4826 on log10(R) over gas-traced shells (not density-clamped, not
    # stellar-dominated). The outlier rejection runs **per shell** — a single noisy radius never
    # discards the whole galaxy; failed shells are simply dropped from ``reliable_measurement``.
    # Bumped from 2σ → 3σ (Median Absolute Deviation) so only severe magnetic spikes are dropped,
    # leaving the galaxy-level audit to ``MIN_RELIABLE_MEASUREMENT_FRACTION`` instead of the gate.
    MAG_PRESSURE_RATIO_OUTLIER_SIGMA: float
    MAG_PRESSURE_RATIO_ROBUST_STATS_MIN_SHELLS: int
    # Absolute linear ceiling on R even after the robust gate (numerical safety / pathological B).
    MAG_PRESSURE_RATIO_ABSOLUTE_CEILING: float
    # Galactic latitude cut (deg) for foreground / MW dust caution.
    MIN_GALACTIC_LATITUDE_DEG: float
    MIN_RELIABLE_MEASUREMENT_FRACTION: float
    EXPLAINED_GAP_HALF_FACTOR: float
    SCENARIO_SUCCESS_RATIO_MIN: float
    SCENARIO_SUCCESS_RATIO_MAX: float
    BEST_SCENARIO_LOG_FLOOR: float
    PHYSICAL_INTERPRETATION_STRONG_MIN: float
    PHYSICAL_INTERPRETATION_STRONG_MAX: float
    PHYSICAL_INTERPRETATION_SPIKE_MIN: float
    DIRECT_EVIDENCE_RATIO_MIN: float
    DIRECT_EVIDENCE_RATIO_MAX: float
    PHOTOMETRIC_LSB_SB_THRESHOLD: float
    NARRATIVE_HIGHLIGHT_COVERAGE_MIN: float
    NARRATIVE_HIGHLIGHT_COVERAGE_MAX: float
    NARRATIVE_HIGHLIGHT_PEARSON_MIN: float
    # If more than this fraction of radial shells had their gas density clamped to MIN_GAS_MASS_DENSITY_KG_M3,
    # the scenario is forced is_reliable=False (the inferred ρ_ion floor is dominating, not measured gas).
    MAX_DENSITY_CLAMPED_FRACTION_FOR_RELIABILITY: float
    # If more than this fraction of radial shells are dominated by stellar contribution
    # (V_disk² + V_bulge² > stellar_fraction · V_bar²), the gas-density-from-V_gas approximation
    # is unreliable and the scenario forfeits is_reliable.
    MAX_STELLAR_DOMINANCE_FRACTION_FOR_RELIABILITY: float
    # Threshold (per shell) used to mark a row as ``stellar_dominated``: dimensionless ratio of
    # stellar V² to total baryonic V². 0.85 means stars supply ≥85% of the squared-velocity budget
    # (relaxed vs 0.70 so intermediate spirals are not over-flagged while still guarding shells
    # where the V_gas-only mass estimator is a poor tracer).
    STELLAR_DOMINATED_PER_SHELL_THRESHOLD: float
    # Plausibility band for the inferred Galactic B field, in microgauss (μG = 1e-10 T).
    # Galactic-disk literature typically reports 1–10 μG; we accept 0.1–100 μG as the
    # ``not implausible'' band before flagging the run.
    MAGNETIC_FIELD_PLAUSIBILITY_MIN_MICROGAUSS: float
    MAGNETIC_FIELD_PLAUSIBILITY_MAX_MICROGAUSS: float
    # Tighter "literature-realistic" Galactic-disk band used by the reporting layer to
    # shade the *expected* range on the magnetic-field-vs-scale plot. This band is a
    # **display annotation only** — the reliability gate above is what controls audit flags.
    MAGNETIC_FIELD_REALISTIC_BAND_MIN_MICROGAUSS: float
    MAGNETIC_FIELD_REALISTIC_BAND_MAX_MICROGAUSS: float
    # Pearson-correlation interpretation thresholds applied uniformly across plots and
    # the narrative. Values below ``MODERATE_MIN`` are reported as "noise"; between
    # ``MODERATE_MIN`` and ``STRONG_MIN`` as "moderate"; above ``STRONG_MIN`` as
    # "strong".
    PEARSON_INTERPRETATION_STRONG_MIN: float
    PEARSON_INTERPRETATION_MODERATE_MIN: float
    # Plasma β = 2 μ₀ P_thermal / B². β ≪ 1 ⇒ magnetic pressure dominates dynamics; β ≫ 1 ⇒ thermal
    # pressure swamps magnetic and the MHD-pressure-driven acceleration hypothesis is questionable.
    MAX_PLASMA_BETA_FOR_MHD_DOMINANCE: float
    # Inclination quality cut: very face-on systems (i < this) inflate any V_obs / sin(i) propagation.
    MIN_INCLINATION_DEG_FOR_RELIABILITY: float
    # Floor on the beam-smearing factor ``k`` when forming ρ_eff = ρ_req / max(k, k_floor) for the
    # WIM density match score (avoids division by zero when the grid-scale diagnostic is missing).
    RESOLUTION_CORRECTION_FACTOR_MIN_FOR_WIM_SCORE: float
    # ``SUPREME`` validation cohort: require the beam-corrected WIM convergence ratio
    # ``(ρ_req / k) / ρ_WIM`` to lie within one order of magnitude of unity (policy gate).
    SUPREME_CALIBRATED_CONVERGENCE_RATIO_MIN: float
    SUPREME_CALIBRATED_CONVERGENCE_RATIO_MAX: float
    # Outer rotation-curve radius proxy: treat the disk diameter as ``multiplier × r_max_kpc``
    # when comparing to the Planck beam footprint (``DATA_LIMITED`` when the beam exceeds it).
    GALAXY_DIAMETER_FROM_RMAX_MULTIPLIER: float
    # Per-shell ``ρ_req`` can hit zeros outside the under-performing mask; when rescaling
    # ``|a_B|`` to the WIM-anchored "honest" diagnostic curve, clip ``ρ_req`` relative to
    # ``ρ_WIM`` so the plot stays finite without hiding physically large requirements.
    HONEST_ACCEL_RHO_REQ_RELATIVE_FLOOR: float


RELIABILITY_CUTS: Final[ReliabilityCuts] = ReliabilityCuts(
    MAX_MAGNETIC_PRESSURE_TO_CENTRIPETAL_ACCEL_FACTOR=float(np.float64(200.0)),
    MAG_PRESSURE_RATIO_OUTLIER_SIGMA=float(np.float64(3.0)),
    MAG_PRESSURE_RATIO_ROBUST_STATS_MIN_SHELLS=int(5),
    MAG_PRESSURE_RATIO_ABSOLUTE_CEILING=float(np.float64(500.0)),
    MIN_GALACTIC_LATITUDE_DEG=float(np.float64(20.0)),
    MIN_RELIABLE_MEASUREMENT_FRACTION=float(np.float64(0.7)),
    EXPLAINED_GAP_HALF_FACTOR=float(np.float64(0.5)),
    SCENARIO_SUCCESS_RATIO_MIN=float(np.float64(0.1)),
    SCENARIO_SUCCESS_RATIO_MAX=float(np.float64(10.0)),
    BEST_SCENARIO_LOG_FLOOR=float(np.float64(1e-12)),
    PHYSICAL_INTERPRETATION_STRONG_MIN=float(np.float64(0.5)),
    PHYSICAL_INTERPRETATION_STRONG_MAX=float(np.float64(1.5)),
    PHYSICAL_INTERPRETATION_SPIKE_MIN=float(np.float64(10.0)),
    DIRECT_EVIDENCE_RATIO_MIN=float(np.float64(0.2)),
    DIRECT_EVIDENCE_RATIO_MAX=float(np.float64(5.0)),
    PHOTOMETRIC_LSB_SB_THRESHOLD=float(np.float64(100.0)),
    NARRATIVE_HIGHLIGHT_COVERAGE_MIN=float(np.float64(0.85)),
    NARRATIVE_HIGHLIGHT_COVERAGE_MAX=float(np.float64(1.25)),
    NARRATIVE_HIGHLIGHT_PEARSON_MIN=float(np.float64(0.85)),
    MAX_DENSITY_CLAMPED_FRACTION_FOR_RELIABILITY=float(np.float64(0.5)),
    MAX_STELLAR_DOMINANCE_FRACTION_FOR_RELIABILITY=float(np.float64(0.6)),
    STELLAR_DOMINATED_PER_SHELL_THRESHOLD=float(np.float64(0.85)),
    MAGNETIC_FIELD_PLAUSIBILITY_MIN_MICROGAUSS=float(np.float64(0.1)),
    MAGNETIC_FIELD_PLAUSIBILITY_MAX_MICROGAUSS=float(np.float64(100.0)),
    MAGNETIC_FIELD_REALISTIC_BAND_MIN_MICROGAUSS=float(np.float64(1.0)),
    MAGNETIC_FIELD_REALISTIC_BAND_MAX_MICROGAUSS=float(np.float64(10.0)),
    PEARSON_INTERPRETATION_STRONG_MIN=float(np.float64(0.7)),
    PEARSON_INTERPRETATION_MODERATE_MIN=float(np.float64(0.4)),
    MAX_PLASMA_BETA_FOR_MHD_DOMINANCE=float(np.float64(1.0)),
    MIN_INCLINATION_DEG_FOR_RELIABILITY=float(np.float64(30.0)),
    RESOLUTION_CORRECTION_FACTOR_MIN_FOR_WIM_SCORE=float(np.float64(1.0e-12)),
    SUPREME_CALIBRATED_CONVERGENCE_RATIO_MIN=float(np.float64(0.1)),
    SUPREME_CALIBRATED_CONVERGENCE_RATIO_MAX=float(np.float64(10.0)),
    GALAXY_DIAMETER_FROM_RMAX_MULTIPLIER=float(np.float64(2.0)),
    HONEST_ACCEL_RHO_REQ_RELATIVE_FLOOR=float(np.float64(1.0e-30)),
)


def luminosity_distance_mpc_from_distance_modulus(distance_modulus_mag: float) -> float:
    """
    Convert a luminosity distance modulus μ (mag) to D (Mpc) using the conventional definition
    ``μ = 5 log10(D_Mpc) + 25`` for extragalactic work in empty space (no ΛCDM lookup).
    """
    mu = float(distance_modulus_mag)
    return float(np.power(10.0, (mu - 25.0) / 5.0))


def luminosity_distance_uncertainty_mpc_from_modulus(
    distance_mpc: float, sigma_distance_modulus_mag: float
) -> float:
    """Propagate a symmetric Gaussian uncertainty on μ to σ_D via ``dD/dμ = D ln(10) / 5``."""
    d = float(distance_mpc)
    sigma_mu = float(sigma_distance_modulus_mag)
    return float(abs(d * np.log(10.0) / 5.0 * sigma_mu))


def effective_required_density_after_beam_kg_m3(
    rho_req_kg_m3: float,
    resolution_correction_factor_k: float,
) -> float:
    """
    Beam-corrected required density ``ρ_eff = ρ_req / max(k, k_floor)`` (kg/m³).

    The Planck beam dilutes gradients by ``k ≈ D_beam / L_grad`` so the reconstructed
    ``ρ_req`` is biased high by the same factor; dividing removes that instrumental tilt
    before comparing to the WIM anchor.
    """
    rho = float(rho_req_kg_m3)
    k = float(resolution_correction_factor_k)
    k_floor = float(RELIABILITY_CUTS.RESOLUTION_CORRECTION_FACTOR_MIN_FOR_WIM_SCORE)
    if not np.isfinite(rho) or rho <= 0.0 or not np.isfinite(k):
        return float("nan")
    denom = max(k, k_floor)
    rho_eff = rho / denom
    return float(rho_eff) if np.isfinite(rho_eff) and rho_eff > 0.0 else float("nan")


def calibrated_convergence_ratio_from_rho_and_k(
    rho_req_kg_m3: float,
    resolution_correction_factor_k: float,
) -> float:
    """
    ``(ρ_req / max(k, k_floor)) / ρ_WIM`` — unity means the match density matches the WIM anchor
    after Planck beam dilution is discounted (grid-search / validation scenario scalar).
    """
    rho_eff = effective_required_density_after_beam_kg_m3(rho_req_kg_m3, resolution_correction_factor_k)
    anchor = float(ASTROPHYSICAL_DEFAULTS.WIM_REFERENCE_DENSITY_KG_M3)
    if not np.isfinite(rho_eff) or rho_eff <= 0.0 or anchor <= 0.0:
        return float("nan")
    return float(rho_eff / anchor)


@dataclass(frozen=True, slots=True)
class VisualizationDisplayDefaults:
    """Purely presentational knobs for ``src.visualization`` (no physics policy)."""

    # Maximum galaxies when resolving capped auto-lists in legacy tooling/tests.
    MAX_AUTO_GALAXIES: int
    # Annulus count for pseudo-3D WIM density rings mapped from ``effective_ionized_mass_density``.
    WIM_RING_BIN_COUNT: int
    # Structured-grid resolution for illustrative B-field streamlines (in-plane).
    STREAMLINE_GRID_SIZE: int
    # Log a warning when an exported Museum HTML exceeds this byte size (browser payload guardrail).
    MUSEUM_HTML_WARN_BYTES: int


VIZ_DISPLAY: Final[VisualizationDisplayDefaults] = VisualizationDisplayDefaults(
    MAX_AUTO_GALAXIES=12,
    WIM_RING_BIN_COUNT=8,
    STREAMLINE_GRID_SIZE=32,
    MUSEUM_HTML_WARN_BYTES=15_000_000,
)


def wim_density_match_score(rho_req_kg_m3: float, resolution_correction_factor: float) -> float:
    """
    Map ``ρ_req`` and the grid beam-smearing factor ``k`` to ``[0, 1]`` via
    ``1 / (1 + |log10( ρ_eff / ρ_anchor )|)`` with ``ρ_eff = ρ_req / max(k, k_floor)``.
    """
    anchor = float(ASTROPHYSICAL_DEFAULTS.WIM_REFERENCE_DENSITY_KG_M3)
    rho_eff = effective_required_density_after_beam_kg_m3(rho_req_kg_m3, resolution_correction_factor)
    if not np.isfinite(rho_eff) or rho_eff <= 0.0 or anchor <= 0.0:
        return float("nan")
    return float(1.0 / (1.0 + abs(np.log10(rho_eff / anchor))))
