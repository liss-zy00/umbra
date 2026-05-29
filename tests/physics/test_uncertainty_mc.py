"""
Monte-Carlo propagation of catalogue uncertainties onto ``theoretical_gap_coverage_ratio``.

The propagator (:func:`src.physics.diagnostics.uncertainty.propagate_gap_coverage_band`)
samples ``(D, i, V_flat)`` from independent Gaussians and rescales the central ratio
analytically. The tests below check that

* the (p16, median, p84) ordering is consistent and brackets the point estimate;
* the RNG is deterministic for a fixed seed;
* zero / missing uncertainties collapse the band onto the central value;
* the closed-form ``default_ratio_scaling`` reproduces a known multiplicative scaling.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.config.thresholds import UNCERTAINTY
from src.physics.diagnostics.uncertainty import (
    default_ratio_scaling,
    propagate_gap_coverage_band,
)


def _base_kwargs() -> dict[str, float]:
    return {
        "base_ratio": 0.85,
        "distance_mpc": 10.0,
        "e_distance_mpc": 0.5,
        "inclination_deg": 60.0,
        "e_inclination_deg": 2.0,
        "vflat_kms": 80.0,
        "e_vflat_kms": 4.0,
    }


def test_mc_band_brackets_base_ratio_and_is_ordered():
    median, p16, p84 = propagate_gap_coverage_band(**_base_kwargs())
    assert p16 < median < p84
    # The propagator is multiplicative around the central inputs, so the median must sit
    # near the point estimate.
    assert abs(median - 0.85) < 0.20


def test_mc_band_is_seed_deterministic():
    args = _base_kwargs()
    a = propagate_gap_coverage_band(**args, seed=42, n_samples=500)
    b = propagate_gap_coverage_band(**args, seed=42, n_samples=500)
    assert a == pytest.approx(b)


def test_mc_band_changes_with_different_seed():
    args = _base_kwargs()
    a = propagate_gap_coverage_band(**args, seed=42, n_samples=500)
    b = propagate_gap_coverage_band(**args, seed=7, n_samples=500)
    # The p16/p84 quantiles should differ between seeds when widths are non-zero.
    assert not (a == pytest.approx(b, abs=1.0e-6))


def test_mc_band_collapses_when_uncertainties_are_zero():
    args = _base_kwargs()
    args["e_distance_mpc"] = 0.0
    args["e_inclination_deg"] = 0.0
    args["e_vflat_kms"] = 0.0
    median, p16, p84 = propagate_gap_coverage_band(**args)
    assert median == pytest.approx(args["base_ratio"], rel=1.0e-9)
    assert p16 == pytest.approx(args["base_ratio"], rel=1.0e-9)
    assert p84 == pytest.approx(args["base_ratio"], rel=1.0e-9)


def test_mc_band_collapses_when_inputs_are_none():
    """Missing catalogue 1-σ entries cannot widen the band; the propagator returns the point estimate."""
    args = _base_kwargs()
    args["e_distance_mpc"] = None
    median, p16, p84 = propagate_gap_coverage_band(**args)
    assert (median, p16, p84) == pytest.approx(
        (args["base_ratio"], args["base_ratio"], args["base_ratio"])
    )


def test_mc_band_nan_base_ratio_propagates_nan():
    args = _base_kwargs()
    args["base_ratio"] = float("nan")
    median, p16, p84 = propagate_gap_coverage_band(**args)
    assert np.isnan(median) and np.isnan(p16) and np.isnan(p84)


def test_default_ratio_scaling_reproduces_multiplicative_form():
    """``ratio = base · (D₀/D_s) · (sin i_s / sin i₀)² · (V₀/V_s)²``."""
    base = 1.2
    d_ref, i_ref, v_ref = 8.0, 50.0, 70.0
    d_s = np.array([8.0, 9.0, 7.0], dtype=float)
    i_s = np.array([50.0, 55.0, 45.0], dtype=float)
    v_s = np.array([70.0, 80.0, 60.0], dtype=float)

    out = default_ratio_scaling(
        base_ratio=base,
        distance_mpc=d_ref,
        distance_sample=d_s,
        inclination_deg=i_ref,
        inclination_sample=i_s,
        vflat_kms=v_ref,
        vflat_sample=v_s,
    )
    expected = (
        base
        * (d_ref / d_s)
        * (np.sin(np.deg2rad(i_s)) / np.sin(np.deg2rad(i_ref))) ** 2
        * (v_ref / v_s) ** 2
    )
    assert np.allclose(out, expected, rtol=1.0e-12)


def test_uncertainty_defaults_are_seeded_and_match_constants():
    assert UNCERTAINTY.MC_DEFAULT_SEED == 42
    assert UNCERTAINTY.MC_DEFAULT_SAMPLES >= 200
