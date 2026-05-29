"""Unit tests for the SOLID diagnostic collaborators decomposed out of PhysicsEngine."""

from __future__ import annotations

import numpy as np
import pytest

from src.config.thresholds import RELIABILITY_CUTS, UNIVERSAL_CONSTANTS
from src.physics.diagnostics.plasma import (
    TESLA_PER_MICROGAUSS,
    MagneticFieldPlausibilityCheck,
    PlasmaDiagnostics,
    ReliabilityMaskBuilder,
)


def test_tesla_per_microgauss_centralised_in_universal_constants() -> None:
    """The conversion must come from ``UNIVERSAL_CONSTANTS`` — no shadow constant in diagnostics."""
    assert TESLA_PER_MICROGAUSS == UNIVERSAL_CONSTANTS.TESLA_PER_MICROGAUSS_SI
    assert TESLA_PER_MICROGAUSS == pytest.approx(1.0e-10)


# ---------------------------- PlasmaDiagnostics ----------------------------


def test_plasma_thermal_pressure_proportional_to_density() -> None:
    """``P_thermal = (ρ/m_p) k_B T`` is linear in ρ — assert the proportionality directly."""
    diag = PlasmaDiagnostics()
    rho = np.array([1.0e-22, 2.0e-22, 0.0], dtype=float)
    p = diag.thermal_pressure(rho)
    assert np.allclose(p[1] / p[0], 2.0)
    assert p[2] == 0.0


def test_plasma_beta_thermal_inverse_square_in_b() -> None:
    """β = 2 μ₀ P / B² → halving B quadruples β at fixed ρ."""
    diag = PlasmaDiagnostics()
    rho = np.array([1.0e-22], dtype=float)
    beta_b = diag.plasma_beta_thermal(rho, magnetic_field_tesla=2.0e-10)
    beta_b_half = diag.plasma_beta_thermal(rho, magnetic_field_tesla=1.0e-10)
    assert np.isclose(beta_b_half[0] / beta_b[0], 4.0)


def test_plasma_beta_thermal_returns_nan_for_zero_field() -> None:
    diag = PlasmaDiagnostics()
    rho = np.array([1.0e-22], dtype=float)
    out = diag.plasma_beta_thermal(rho, magnetic_field_tesla=0.0)
    assert np.isnan(out[0])


def test_alfven_velocity_increases_with_b_decreases_with_density() -> None:
    """v_A = B / sqrt(μ₀ ρ_ion). Doubling B doubles v_A; quadrupling ρ halves v_A."""
    diag = PlasmaDiagnostics()
    rho = np.array([1.0e-22], dtype=float)
    v1 = diag.alfven_velocity_ms(magnetic_field_tesla=1.0e-10, effective_ionized_mass_density_kg_m3=rho)
    v2 = diag.alfven_velocity_ms(magnetic_field_tesla=2.0e-10, effective_ionized_mass_density_kg_m3=rho)
    assert np.isclose(v2[0] / v1[0], 2.0)
    rho4 = np.array([4.0e-22], dtype=float)
    v_rho4 = diag.alfven_velocity_ms(magnetic_field_tesla=1.0e-10, effective_ionized_mass_density_kg_m3=rho4)
    assert np.isclose(v_rho4[0] / v1[0], 0.5)


def test_alfven_velocity_returns_nan_for_zero_density() -> None:
    diag = PlasmaDiagnostics()
    rho = np.array([0.0, -1.0], dtype=float)
    out = diag.alfven_velocity_ms(magnetic_field_tesla=1.0e-10, effective_ionized_mass_density_kg_m3=rho)
    assert np.all(np.isnan(out))


# ---------------------- MagneticFieldPlausibilityCheck ----------------------


def test_plausibility_band_default_uses_central_thresholds() -> None:
    chk = MagneticFieldPlausibilityCheck()
    assert chk.is_within_band(RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MIN_MICROGAUSS) is True
    assert chk.is_within_band(RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MAX_MICROGAUSS) is True
    assert chk.is_within_band(RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MIN_MICROGAUSS / 10.0) is False
    assert chk.is_within_band(RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MAX_MICROGAUSS * 10.0) is False


def test_plausibility_band_rejects_nan() -> None:
    chk = MagneticFieldPlausibilityCheck()
    assert chk.is_within_band(float("nan")) is False


def test_plausibility_band_with_custom_overrides() -> None:
    """Demonstration of the dependency-injection contract: override is honored without globals."""
    chk = MagneticFieldPlausibilityCheck(min_microgauss=5.0, max_microgauss=7.0)
    assert chk.is_within_band(6.0) is True
    assert chk.is_within_band(4.999) is False


# ----------------------------- ReliabilityMaskBuilder -----------------------------


def test_reliability_mask_drops_only_outlier_shell() -> None:
    """One outrageously high |a_B|/|a_c| shell must be dropped — the others survive."""
    builder = ReliabilityMaskBuilder()
    # 9 shells with |a_B|/|a_c| ≈ 1, plus one extreme spike at 1e6.
    base = np.array([1.0e-10] * 9 + [1.0e-4], dtype=float)
    cent = np.full_like(base, 1.0e-10)
    mask = builder.build_mask(
        magnetic_pressure_acceleration=base,
        centripetal_acceleration=cent,
        density_clamped_mask=np.zeros_like(base, dtype=bool),
        stellar_dominated_mask=np.zeros_like(base, dtype=bool),
    )
    assert mask[-1] == False  # outlier dropped
    assert int(np.count_nonzero(mask)) == 9


def test_reliability_mask_uses_legacy_cap_when_too_few_gas_shells() -> None:
    """If most shells are stellar-dominated/clamped the robust statistic cannot be built and the
    legacy linear factor applies. Below the cap they should all survive."""
    builder = ReliabilityMaskBuilder()
    a_b = np.array([1.0e-10] * 4, dtype=float)
    a_c = np.array([1.0e-10] * 4, dtype=float)
    # Only one shell is gas-traced (well below MIN_SHELLS=5) — the legacy linear cap kicks in.
    dcm = np.array([True, True, True, False], dtype=bool)
    sdm = np.zeros_like(dcm)
    mask = builder.build_mask(a_b, a_c, dcm, sdm)
    # All four shells (including clamped) are still finite and below the linear ceiling, so the
    # builder returns ``finite & ratio<=cap`` (the higher-level engine OR-combines this with the
    # exclusion masks separately).
    assert int(np.count_nonzero(mask)) == 4


def test_reliability_mask_drops_non_finite_shells() -> None:
    builder = ReliabilityMaskBuilder()
    a_b = np.array([1.0e-10, np.nan, np.inf, 1.0e-10], dtype=float)
    a_c = np.array([1.0e-10, 1.0e-10, 1.0e-10, 0.0], dtype=float)
    mask = builder.build_mask(
        a_b, a_c, np.zeros_like(a_b, dtype=bool), np.zeros_like(a_b, dtype=bool)
    )
    assert mask[0] == True
    assert mask[1] == False
    assert mask[2] == False
    assert mask[3] == False  # a_c == 0 fails the eps gate
