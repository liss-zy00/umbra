"""Unit tests for the McGaugh+2016 RAR interpolating function."""

from __future__ import annotations

import numpy as np

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS
from src.physics.diagnostics.rar import rar_predicted_observed_acceleration, rar_residual_log10


def test_rar_high_acceleration_limit_recovers_baryonic():
    """For ``g_bar >> a0`` the RAR collapses to ``g_obs ≈ g_bar``."""
    a0 = ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2
    g_bar = 1e3 * a0  # deep in the high-acceleration regime
    g_obs = float(rar_predicted_observed_acceleration(g_bar))
    assert np.isclose(g_obs, g_bar, rtol=1e-3)


def test_rar_low_acceleration_limit_recovers_sqrt_a0_g_bar():
    """For ``g_bar << a0`` the RAR collapses to ``g_obs ≈ sqrt(g_bar · a0)`` (deep MOND-like)."""
    a0 = ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2
    g_bar = 1e-4 * a0
    g_obs = float(rar_predicted_observed_acceleration(g_bar))
    assert np.isclose(g_obs, np.sqrt(g_bar * a0), rtol=5e-2)
    # And the prediction is strictly above g_bar in this regime (DM-like uplift).
    assert g_obs > g_bar


def test_rar_handles_nonpositive_and_nan_inputs_without_crash():
    g_bar = np.array([np.nan, -1.0, 0.0, 1e-12], dtype=float)
    g_obs = rar_predicted_observed_acceleration(g_bar)
    assert g_obs.shape == g_bar.shape
    assert np.isnan(g_obs[0])
    assert np.isnan(g_obs[1])
    assert np.isnan(g_obs[2])
    assert np.isfinite(g_obs[3])


def test_rar_residual_zero_when_curve_matches_prediction():
    a0 = ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2
    g_bar = np.array([1e-2, 1.0, 100.0]) * a0
    g_obs = rar_predicted_observed_acceleration(g_bar)
    residual = rar_residual_log10(g_obs, g_bar)
    assert np.allclose(residual, 0.0, atol=1e-12)
