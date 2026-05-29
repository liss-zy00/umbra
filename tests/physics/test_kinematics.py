"""Unit tests for the BaryonicGapCalculator (extracted from PhysicsEngine)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from astropy import units as u

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS
from src.physics.kinematics import BaryonicGapCalculator


def _base_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "radius_kpc": [1.0, 2.0, 3.0],
            "observed_velocity_kms": [100.0, 110.0, 120.0],
            "baryonic_velocity_kms": [60.0, 65.0, 70.0],
        }
    )


def test_compute_returns_si_columns_with_correct_conversion() -> None:
    calc = BaryonicGapCalculator()
    out = calc.compute(_base_df())
    expected_radius_m = np.array([1.0, 2.0, 3.0]) * u.kpc.to(u.m)
    assert np.allclose(out["radius_m"], expected_radius_m)
    assert np.allclose(out["observed_velocity_ms"], np.array([100.0, 110.0, 120.0]) * 1000.0)
    assert np.allclose(out["baryonic_velocity_ms"], np.array([60.0, 65.0, 70.0]) * 1000.0)


def test_acceleration_gap_equals_centripetal_minus_baryonic() -> None:
    calc = BaryonicGapCalculator()
    out = calc.compute(_base_df())
    cent = out["observed_velocity_ms"] ** 2 / out["radius_m"]
    bar = out["baryonic_velocity_ms"] ** 2 / out["radius_m"]
    assert np.allclose(out["centripetal_acceleration"], cent)
    assert np.allclose(out["baryonic_acceleration"], bar)
    assert np.allclose(out["acceleration_gap"], cent - bar)


def test_compute_does_not_mutate_input() -> None:
    """The estimator returns a copy — the input dataframe must remain unchanged."""
    df = _base_df()
    snapshot = df.copy()
    BaryonicGapCalculator().compute(df)
    pd.testing.assert_frame_equal(df, snapshot)


def test_missing_columns_raise_value_error() -> None:
    calc = BaryonicGapCalculator()
    df = pd.DataFrame({"radius_kpc": [1.0]})
    with pytest.raises(ValueError, match="Missing columns"):
        calc.compute(df)


def test_invalid_radius_raises() -> None:
    calc = BaryonicGapCalculator()
    df = pd.DataFrame(
        {
            "radius_kpc": [0.0],
            "observed_velocity_kms": [100.0],
            "baryonic_velocity_kms": [50.0],
        }
    )
    with pytest.raises(ValueError, match="R must be finite"):
        calc.compute(df)


def test_inclination_systematic_propagated_with_central_factor() -> None:
    """δa_c / a_c = K · |δV/V| with K from ASTROPHYSICAL_DEFAULTS."""
    calc = BaryonicGapCalculator()
    df = _base_df()
    df["observed_velocity_systematic_kms"] = [10.0, 11.0, 12.0]
    out = calc.compute(df)
    rel = np.array([10.0 / 100.0, 11.0 / 110.0, 12.0 / 120.0])
    expected = (
        ASTROPHYSICAL_DEFAULTS.INCLINATION_SYSTEMATIC_PROPAGATION_FACTOR
        * np.abs(out["centripetal_acceleration"].to_numpy(dtype=float))
        * rel
    )
    assert np.allclose(out["centripetal_acceleration_systematic"], expected)


def test_inclination_systematic_zero_when_column_missing() -> None:
    calc = BaryonicGapCalculator()
    out = calc.compute(_base_df())
    assert np.all(out["centripetal_acceleration_systematic"] == 0.0)
