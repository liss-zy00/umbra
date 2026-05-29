"""
Magnetic-tension toy term (logarithmic-spiral curvature).

Tests cover three facets of :mod:`src.physics.diagnostics.tension`:

* the closed-form curvature radius ``R_c = r / sin² p``;
* the radial acceleration magnitude on a spiral morphology against the analytical
  identity ``|a_tension,r| = cos p · B² / (μ₀ ρ_ion R_c)``;
* the morphology gate — irregular / BCD morphologies force the tension array to zero
  (the caller is expected to raise ``warn_pure_pressure`` in that case).
"""

from __future__ import annotations

import numpy as np

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS, UNIVERSAL_CONSTANTS
from src.physics.diagnostics.tension import (
    curvature_radius_m,
    supports_logarithmic_spiral,
    tension_acceleration_radial,
)


def test_curvature_radius_matches_log_spiral_identity():
    pitch_deg = 12.0
    radius_m = np.array([1.0e19, 5.0e19, 1.0e20], dtype=float)
    expected = radius_m / (np.sin(np.deg2rad(pitch_deg)) ** 2)
    out = curvature_radius_m(radius_m, pitch_deg=pitch_deg)
    assert np.allclose(out, expected, rtol=1.0e-12)


def test_curvature_radius_rejects_invalid_pitch():
    radius_m = np.array([1.0e19, 5.0e19], dtype=float)
    for bad in (0.0, 90.0, -1.0, 95.0):
        out = curvature_radius_m(radius_m, pitch_deg=bad)
        assert np.all(np.isnan(out))


def test_tension_acceleration_matches_closed_form_on_spiral():
    """For a spiral host, ``|a_tension,r| = cos p · B² / (μ₀ ρ R_c)`` shell-by-shell."""
    radius_m = np.array([1.0e19, 2.0e19, 5.0e19], dtype=float)
    b_tesla = 5.0e-10
    rho = np.array([1.0e-22, 5.0e-23, 1.0e-23], dtype=float)
    pitch_deg = ASTROPHYSICAL_DEFAULTS.SPIRAL_PITCH_DEG
    mu0 = UNIVERSAL_CONSTANTS.MU0_SI

    rc = radius_m / (np.sin(np.deg2rad(pitch_deg)) ** 2)
    cos_p = float(np.cos(np.deg2rad(pitch_deg)))
    expected = -cos_p * (b_tesla ** 2) / (mu0 * rho * rc)

    out = tension_acceleration_radial(
        radius_m=radius_m,
        magnetic_field_tesla=b_tesla,
        effective_ionized_mass_density=rho,
        morphology_class="Late-type Spirals (Sc, Scd)",
        pitch_deg=pitch_deg,
    )
    assert np.allclose(out, expected, rtol=1.0e-12)
    # The radial component must be negative for an outward-pointing pressure gradient,
    # i.e. tension *opposes* the centrifugal direction.
    assert np.all(out < 0.0)


def test_tension_collapses_to_zero_on_irregular_morphology():
    radius_m = np.array([1.0e19, 5.0e19], dtype=float)
    rho = np.array([1.0e-22, 5.0e-23], dtype=float)
    out = tension_acceleration_radial(
        radius_m=radius_m,
        magnetic_field_tesla=1.0e-10,
        effective_ionized_mass_density=rho,
        morphology_class="Dwarf Irregulars (Im)",
    )
    assert np.allclose(out, 0.0)


def test_morphology_gate_excludes_bcd_and_im():
    assert supports_logarithmic_spiral("Late-type Spirals (Sb, Sbc)") is True
    assert supports_logarithmic_spiral("Dwarf Irregulars (Im)") is False
    assert supports_logarithmic_spiral("BCD") is False
    assert supports_logarithmic_spiral(None) is False
    assert supports_logarithmic_spiral("") is False


def test_tension_returns_nan_where_density_is_nonpositive():
    """``ρ_ion ≤ 0`` shells must propagate NaN, not zero (caller flags as unreliable)."""
    radius_m = np.array([1.0e19, 2.0e19, 5.0e19], dtype=float)
    rho = np.array([1.0e-22, 0.0, -1.0e-23], dtype=float)
    out = tension_acceleration_radial(
        radius_m=radius_m,
        magnetic_field_tesla=5.0e-10,
        effective_ionized_mass_density=rho,
        morphology_class="Sa",
    )
    assert np.isfinite(out[0])
    assert np.isnan(out[1])
    assert np.isnan(out[2])
