"""
Baryonic Tully–Fisher Relation (McGaugh+2012) per-galaxy check.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS
from src.physics.diagnostics.btfr import (
    MCGAUGH_2012_RIDGE,
    baryonic_mass_solar,
    btfr_residual_dex,
    is_btfr_outlier,
)


def test_mcgaugh_ridge_constants_match_catalogue():
    """The ridge baked into the module must reproduce the canonical McGaugh+2012 values."""
    assert MCGAUGH_2012_RIDGE.slope == pytest.approx(4.0)
    assert MCGAUGH_2012_RIDGE.intercept == pytest.approx(1.99)


def test_residual_zero_on_canonical_ridge_galaxy():
    """A galaxy whose ``M_bar`` sits exactly on the ridge yields ``residual = 0``."""
    vflat = 100.0
    log_m_pred = MCGAUGH_2012_RIDGE.slope * math.log10(vflat) + MCGAUGH_2012_RIDGE.intercept
    m_bar_target = 10.0 ** log_m_pred
    # We force a galaxy whose baryonic mass equals the ridge prediction (HI leg only).
    m_hi_1e9 = m_bar_target / (ASTROPHYSICAL_DEFAULTS.HELIUM_CORRECTION_FACTOR_HI_TO_GAS * 1.0e9)
    residual = btfr_residual_dex(
        vflat_kms=vflat,
        hi_mass_1e9_solar=m_hi_1e9,
        luminosity_3p6um_1e9_solar=None,
    )
    assert residual == pytest.approx(0.0, abs=1.0e-9)
    assert is_btfr_outlier(residual) is False


def test_outlier_flag_flips_at_configured_threshold():
    """Residuals just above ``BTFR_OUTLIER_THRESHOLD_DEX`` must flip the outlier flag."""
    threshold = ASTROPHYSICAL_DEFAULTS.BTFR_OUTLIER_THRESHOLD_DEX
    assert is_btfr_outlier(threshold + 0.01) is True
    assert is_btfr_outlier(-(threshold + 0.01)) is True
    assert is_btfr_outlier(threshold - 0.01) is False
    assert is_btfr_outlier(0.0) is False
    assert is_btfr_outlier(None) is False
    assert is_btfr_outlier(float("nan")) is False


def test_baryonic_mass_combines_hi_and_stellar_legs():
    m = baryonic_mass_solar(hi_mass_1e9_solar=1.0, luminosity_3p6um_1e9_solar=2.0)
    helium = ASTROPHYSICAL_DEFAULTS.HELIUM_CORRECTION_FACTOR_HI_TO_GAS
    mass_to_light = ASTROPHYSICAL_DEFAULTS.STELLAR_MASS_TO_LIGHT_3P6UM
    expected = helium * 1.0 * 1.0e9 + mass_to_light * 2.0 * 1.0e9
    assert m == pytest.approx(expected, rel=1.0e-12)


def test_baryonic_mass_handles_missing_legs():
    m_hi_only = baryonic_mass_solar(hi_mass_1e9_solar=2.0, luminosity_3p6um_1e9_solar=None)
    m_star_only = baryonic_mass_solar(hi_mass_1e9_solar=None, luminosity_3p6um_1e9_solar=2.0)
    assert m_hi_only > 0.0
    assert m_star_only > 0.0
    assert np.isnan(baryonic_mass_solar(hi_mass_1e9_solar=None, luminosity_3p6um_1e9_solar=None))


def test_residual_returns_nan_on_missing_inputs():
    assert np.isnan(
        btfr_residual_dex(vflat_kms=None, hi_mass_1e9_solar=1.0, luminosity_3p6um_1e9_solar=1.0)
    )
    assert np.isnan(
        btfr_residual_dex(vflat_kms=-10.0, hi_mass_1e9_solar=1.0, luminosity_3p6um_1e9_solar=1.0)
    )
    assert np.isnan(
        btfr_residual_dex(vflat_kms=100.0, hi_mass_1e9_solar=None, luminosity_3p6um_1e9_solar=None)
    )
