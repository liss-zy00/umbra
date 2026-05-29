"""
Auxiliary plasma- and reliability diagnostics for the magnetic-pressure pipeline.

Each class in this module has a **single responsibility** and is composed (not subclassed)
into :class:`src.physics.engine.PhysicsEngine`. Splitting these out of the main engine keeps
the orchestrator small and lets the diagnostics be unit-tested in isolation.

Physical content (in academic notation, SI throughout)
------------------------------------------------------
* **Plasma β**:
      β_thermal = (2 μ₀ · n · k_B · T) / B²,
  where ``n = ρ_gas / m_p`` is the proton number density. β ≪ 1 marks the regime in which
  the ideal-MHD pressure-gradient ansatz (``a_B = -∇P_B / ρ_ion``) is locally consistent;
  β ≫ 1 marks shells where thermal pressure dominates and the MHD ansatz is *out of
  regime* — those shells should not be over-interpreted.

* **Alfvén velocity**:
      v_A = B / √(μ₀ · ρ_ion),
  the speed of small-amplitude transverse magnetic perturbations in an ideal-MHD fluid;
  reported per shell so the report can flag galaxies whose disks would require Alfvén
  speeds incompatible with observed turbulence levels.

* **Magnetic-field plausibility band**:
      B ∈ [MAGNETIC_FIELD_PLAUSIBILITY_MIN_MICROGAUSS, MAX_MICROGAUSS]
  A configurable band centred on the literature-consistent disk B-field range. Galaxies
  whose inferred B falls outside the band raise a ``warn_b_band`` audit flag — the
  CMB-energy-density anchor is an order-of-magnitude bridge, not a measurement.

* **Reliability mask** (per-shell, *never* galaxy-level):
  Robust outlier rejection on ``log10(|a_B|/|a_c|)`` using median + Nσ × MAD with a hard
  ceiling on the linear ratio. Galaxy-level reliability is the responsibility of
  :mod:`src.pipeline.per_galaxy.metrics_builder` (``audit_flags``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS, RELIABILITY_CUTS, UNIVERSAL_CONSTANTS

__all__ = [
    "TESLA_PER_MICROGAUSS",
    "PlasmaDiagnostics",
    "MagneticFieldPlausibilityCheck",
    "ReliabilityMaskBuilder",
]

#: One microgauss expressed in Tesla. Single source of truth lives in
#: ``UNIVERSAL_CONSTANTS.TESLA_PER_MICROGAUSS_SI``; this re-export keeps callers free of
#: the dataclass coupling.
TESLA_PER_MICROGAUSS: float = UNIVERSAL_CONSTANTS.TESLA_PER_MICROGAUSS_SI


@dataclass(frozen=True)
class PlasmaDiagnostics:
    """
    Compute thermal plasma-β and Alfvén speed for an MHD disk shell.

    Parameters of the diagnostics are defaulted to the pipeline-wide constants but kept as
    explicit fields so unit tests can override them without monkey-patching globals.

    Attributes
    ----------
    mu0 : float
        Vacuum magnetic permeability in SI units (H/m).
    proton_mass_kg : float
        Proton rest mass (kg). Sets ``n = ρ / m_p``.
    boltzmann_j_per_k : float
        Boltzmann constant (J/K).
    plasma_temperature_k : float
        Thermal temperature (K) used for the warm-ionized-medium proxy. The current pipeline
        uses a single, configurable temperature; per-galaxy thermometry is future work
        (see ``docs/theory/EG_MODEL_THEORY_EN.md`` Section 5).
    """

    mu0: float = UNIVERSAL_CONSTANTS.MU0_SI
    proton_mass_kg: float = ASTROPHYSICAL_DEFAULTS.PROTON_MASS_KG
    boltzmann_j_per_k: float = ASTROPHYSICAL_DEFAULTS.BOLTZMANN_J_PER_K
    plasma_temperature_k: float = ASTROPHYSICAL_DEFAULTS.PLASMA_THERMAL_TEMPERATURE_K

    def thermal_pressure(self, gas_mass_density_kg_m3: np.ndarray) -> np.ndarray:
        """Return ``P_thermal = n_p · k_B · T`` in pascals, per shell."""
        rho = np.asarray(gas_mass_density_kg_m3, dtype=float)
        return (rho / self.proton_mass_kg) * self.boltzmann_j_per_k * self.plasma_temperature_k

    def plasma_beta_thermal(
        self,
        gas_mass_density_kg_m3: np.ndarray,
        magnetic_field_tesla: float,
    ) -> np.ndarray:
        """Return per-shell ``β_thermal = 2 μ₀ P_thermal / B²``.

        When ``B² ≤ 0`` the result is NaN — β is meaningless in that limit.
        """
        thermal_pressure = self.thermal_pressure(gas_mass_density_kg_m3)
        b_sq = float(magnetic_field_tesla) ** 2
        if b_sq > 0:
            return 2.0 * self.mu0 * thermal_pressure / b_sq
        return np.full_like(thermal_pressure, np.nan, dtype=float)

    def alfven_velocity_ms(
        self,
        magnetic_field_tesla: float,
        effective_ionized_mass_density_kg_m3: np.ndarray,
    ) -> np.ndarray:
        """Return per-shell ``v_A = |B| / √(μ₀ · ρ_ion)`` in m/s.

        ``ρ_ion ≤ 0`` shells return NaN. The pipeline reports the median of finite values.
        """
        denom = self.mu0 * np.asarray(effective_ionized_mass_density_kg_m3, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(
                denom > 0,
                np.abs(float(magnetic_field_tesla)) / np.sqrt(denom),
                np.nan,
            )


@dataclass(frozen=True)
class MagneticFieldPlausibilityCheck:
    """
    Boolean gate: is the inferred μG field inside the literature-consistent disk band?

    The band is defined by ``RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MIN/MAX_MICROGAUSS``.
    Galaxies whose inferred field falls outside the band are not automatically rejected —
    instead the metrics builder records a ``warn_b_band`` audit flag so the assumption
    surfaces in the per-galaxy summary and the population-level audit report.
    """

    min_microgauss: float = RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MIN_MICROGAUSS
    max_microgauss: float = RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MAX_MICROGAUSS

    def is_within_band(self, magnetic_field_microgauss: float) -> bool:
        """Return ``True`` when ``B`` (μG) sits inside ``[min_microgauss, max_microgauss]``."""
        b_value = float(magnetic_field_microgauss)
        if not np.isfinite(b_value):
            return False
        return self.min_microgauss <= b_value <= self.max_microgauss


@dataclass(frozen=True)
class ReliabilityMaskBuilder:
    """
    Build a per-shell reliability mask for the magnetic-pressure / centripetal ratio.

    The estimator is **per-shell, never galaxy-level**. Shells flagged as outliers are
    dropped from ``reliable_measurement``; the rest of the galaxy keeps contributing.
    Galaxy-level reliability (audit flags) lives in :mod:`src.pipeline.per_galaxy.metrics_builder`.

    Algorithm
    ---------
    1. Compute ``ratio = |a_B| / |a_c|`` and ``log_r = log10(ratio)`` on shells where both
       accelerations are finite and ``a_c > 0``.
    2. Restrict to shells that are *gas-traced* (not ``density_clamped`` and not
       ``stellar_dominated``); call the count ``n_stat``.
    3. If ``n_stat ≥ MAG_PRESSURE_RATIO_ROBUST_STATS_MIN_SHELLS``:
           threshold = median(log_r) + Nσ · σ_robust,
           σ_robust = max(MAD × 1.4826, 0.05),
           Nσ = MAG_PRESSURE_RATIO_OUTLIER_SIGMA  (currently 3).
       A linear ``MAG_PRESSURE_RATIO_ABSOLUTE_CEILING`` cap is applied as a numerical
       safety net.
    4. Otherwise (sparse curves) use the fixed factor
       ``MAX_MAGNETIC_PRESSURE_TO_CENTRIPETAL_ACCEL_FACTOR`` on every finite shell.
    """

    eps_centripetal: float = 1e-30
    log_floor: float = 1e-300
    sigma_min: float = 0.05
    default_sigma: float = 0.1
    mad_to_sigma: float = 1.4826

    def build_mask(
        self,
        magnetic_pressure_acceleration: np.ndarray,
        centripetal_acceleration: np.ndarray,
        density_clamped_mask: np.ndarray,
        stellar_dominated_mask: np.ndarray,
    ) -> np.ndarray:
        """Return the boolean ``within`` mask described in the class docstring."""
        cent = np.abs(np.asarray(centripetal_acceleration, dtype=float))
        mag = np.abs(np.asarray(magnetic_pressure_acceleration, dtype=float))
        n = cent.size
        finite_base = np.isfinite(cent) & np.isfinite(mag) & (cent > self.eps_centripetal)
        ratio = np.full(n, np.nan, dtype=float)
        ratio[finite_base] = mag[finite_base] / cent[finite_base]
        with np.errstate(divide="ignore", invalid="ignore"):
            log_r = np.full(n, np.nan, dtype=float)
            log_r[finite_base] = np.log10(np.maximum(ratio[finite_base], self.log_floor))

        dcm = np.asarray(density_clamped_mask, dtype=bool)
        sdm = np.asarray(stellar_dominated_mask, dtype=bool)
        stat_mask = finite_base & ~dcm & ~sdm
        n_stat = int(np.count_nonzero(stat_mask))
        min_s = int(RELIABILITY_CUTS.MAG_PRESSURE_RATIO_ROBUST_STATS_MIN_SHELLS)

        if n_stat >= min_s:
            sample = log_r[stat_mask]
            med = float(np.median(sample))
            mad = float(np.median(np.abs(sample - med)))
            if mad > 1e-12:
                sigma_robust = self.mad_to_sigma * mad
            elif sample.size > 1:
                sigma_robust = float(np.std(sample))
            else:
                sigma_robust = self.default_sigma
            sigma_robust = max(sigma_robust, self.sigma_min)
            thresh_log = med + float(RELIABILITY_CUTS.MAG_PRESSURE_RATIO_OUTLIER_SIGMA) * sigma_robust
            ceiling = float(RELIABILITY_CUTS.MAG_PRESSURE_RATIO_ABSOLUTE_CEILING)
            return finite_base & (log_r <= thresh_log) & (ratio <= ceiling)
        cap = float(RELIABILITY_CUTS.MAX_MAGNETIC_PRESSURE_TO_CENTRIPETAL_ACCEL_FACTOR)
        return finite_base & (ratio <= cap)
