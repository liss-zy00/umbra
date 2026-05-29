"""
Unit tests for :mod:`src.physics.diagnostics.resolution`.

Covers:

* closed-form value of ``D_beam = D · θ_beam`` (5'  → 1.454e-3 rad) against a hand-computed
  reference at ``D = 4 Mpc`` (≈ 5.8 kpc);
* dimensionless dilution factor ``k = D_beam / L_grad`` for a representative
  ``L_grad = 100 pc``;
* propagation of ``NaN`` / non-positive inputs through the public API;
* vectorised population helper.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.config.thresholds import PLANCK_CALIBRATION
from src.physics import (
    beam_fwhm_radians,
    beam_physical_size_kpc,
    resolution_correction_factor,
    vectorised_resolution_correction,
)


def test_beam_fwhm_radians_matches_planck_calibration():
    """The exported FWHM in radians must equal ``π · arcmin / (180 · 60)``."""
    expected = PLANCK_CALIBRATION.BEAM_FWHM_ARCMIN * math.pi / (180.0 * 60.0)
    assert beam_fwhm_radians() == pytest.approx(expected, rel=1e-12)


def test_beam_physical_size_kpc_closed_form_at_4Mpc():
    """``D_beam = D · θ_beam`` at D=4 Mpc, θ=5'  ≈ 5.818 kpc (4·1.4544e-3·1000)."""
    out = beam_physical_size_kpc(distance_mpc=4.0)
    expected = 4.0 * (5.0 * math.pi / (180.0 * 60.0)) * 1000.0
    assert math.isclose(out, expected, rel_tol=1e-12)
    # Sanity: the 5'  FWHM scales to single-digit kpc for typical SPARC distances.
    assert 5.5 <= out <= 6.1


def test_beam_physical_size_kpc_handles_missing_or_invalid_distance():
    for bad in (None, float("nan"), -1.0, 0.0, "not-a-number"):
        result = beam_physical_size_kpc(distance_mpc=bad)  # type: ignore[arg-type]
        assert math.isnan(result), f"expected NaN for distance={bad!r}, got {result}"


def test_resolution_correction_factor_grows_linearly_with_distance():
    """At fixed ``L_grad`` the dilution factor must scale linearly with ``D``."""
    k1 = resolution_correction_factor(distance_mpc=4.0, gradient_scale_kpc=0.1)
    k2 = resolution_correction_factor(distance_mpc=8.0, gradient_scale_kpc=0.1)
    assert math.isfinite(k1) and math.isfinite(k2)
    assert k2 == pytest.approx(2.0 * k1, rel=1e-10)


def test_resolution_correction_factor_inverse_in_gradient_scale():
    """At fixed ``D`` doubling ``L_grad`` must halve ``k``."""
    k_fine = resolution_correction_factor(distance_mpc=4.0, gradient_scale_kpc=0.05)
    k_coarse = resolution_correction_factor(distance_mpc=4.0, gradient_scale_kpc=0.1)
    assert math.isfinite(k_fine) and math.isfinite(k_coarse)
    assert k_coarse == pytest.approx(0.5 * k_fine, rel=1e-10)


def test_resolution_correction_factor_typical_galaxy_is_severely_diluted():
    """For NGC-3741-like distance (~3 Mpc) and a 100 pc gradient scale, k must be >> 1."""
    k = resolution_correction_factor(distance_mpc=3.0, gradient_scale_kpc=0.1)
    assert k > 40.0  # 3 · 5'  / 100 pc ≈ 43.6 — well above 1.


@pytest.mark.parametrize(
    "distance_mpc,gradient_scale_kpc",
    [
        (None, 0.1),
        (4.0, None),
        (float("nan"), 0.1),
        (4.0, float("nan")),
        (-1.0, 0.1),
        (4.0, 0.0),
    ],
)
def test_resolution_correction_factor_returns_nan_on_bad_inputs(distance_mpc, gradient_scale_kpc):
    out = resolution_correction_factor(
        distance_mpc=distance_mpc, gradient_scale_kpc=gradient_scale_kpc
    )
    assert math.isnan(out)


def test_vectorised_resolution_correction_matches_scalar_path():
    distances = [3.0, 4.0, float("nan"), 0.0, 10.0]
    scales = [0.1, 0.1, 0.1, 0.1, 0.1]
    vec = vectorised_resolution_correction(distance_mpc=distances, gradient_scale_kpc=scales)
    assert vec.shape == (5,)
    for i, d in enumerate(distances):
        expected = resolution_correction_factor(
            distance_mpc=d, gradient_scale_kpc=scales[i]
        )
        if math.isnan(expected):
            assert math.isnan(vec[i])
        else:
            assert vec[i] == pytest.approx(expected, rel=1e-12)


def test_resolution_correction_factor_with_non_default_beam_size():
    """Doubling the beam FWHM must double the dilution factor."""
    k_default = resolution_correction_factor(distance_mpc=4.0, gradient_scale_kpc=0.1)
    k_wider = resolution_correction_factor(
        distance_mpc=4.0,
        gradient_scale_kpc=0.1,
        beam_fwhm_arcmin=2.0 * PLANCK_CALIBRATION.BEAM_FWHM_ARCMIN,
    )
    assert k_wider == pytest.approx(2.0 * k_default, rel=1e-12)
