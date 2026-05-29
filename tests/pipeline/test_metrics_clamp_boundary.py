"""
Boundary tests for the density-clamp and stellar-dominance reliability gates inside
``metrics_for_successful_physics``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.config.thresholds import RELIABILITY_CUTS
from src.pipeline.per_galaxy.metrics_builder import metrics_for_successful_physics


def _scenario_df(
    *,
    n: int = 10,
    n_clamped: int = 0,
    n_stellar_dominated: int = 0,
    reliable_per_row: bool = True,
) -> pd.DataFrame:
    """Build a minimal scenario dataframe with the columns the metrics builder reads."""
    clamped = np.zeros(n, dtype=bool)
    clamped[:n_clamped] = True
    stellar = np.zeros(n, dtype=bool)
    stellar[:n_stellar_dominated] = True
    reliable = np.full(n, reliable_per_row, dtype=bool)
    # Make reliable points consistent: a row that was clamped or stellar-dominated cannot still be
    # ``reliable`` from the engine's standpoint.
    reliable[clamped] = False
    reliable[stellar] = False
    return pd.DataFrame(
        {
            "reliable_measurement": reliable,
            "acceleration_gap": np.linspace(1e-12, 1e-11, n),
            "magnetic_pressure_acceleration": np.linspace(1e-12, 1e-11, n),
            "density_floor_clamped": clamped,
            "stellar_dominated": stellar,
            "rar_residual_log10": np.zeros(n, dtype=float),
            "magnetic_field_microgauss": np.full(n, 1.0, dtype=float),
            "plasma_beta_thermal": np.full(n, 1e-3, dtype=float),
        }
    )


def test_metrics_passes_when_clamped_fraction_below_threshold():
    # Use n=20 so we can clamp a small fraction and still keep enough reliable rows to pass
    # MIN_RELIABLE_MEASUREMENT_FRACTION (which is 0.7).
    n = 20
    n_clamped = 2  # 10% — well below MAX_DENSITY_CLAMPED_FRACTION_FOR_RELIABILITY (0.5).
    df = _scenario_df(n=n, n_clamped=n_clamped)
    out = metrics_for_successful_physics(
        galaxy_name="X",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-11,
        galactic_latitude_degrees=45.0,
        morphology_class="dummy",
    )
    assert out["density_audit_passed"] is True
    assert out["is_reliable"] is True


def test_metrics_fails_when_clamped_fraction_above_threshold():
    n = 10
    # Just above threshold ⇒ scenario forfeits is_reliable regardless of other gates.
    n_clamped = int(np.ceil(RELIABILITY_CUTS.MAX_DENSITY_CLAMPED_FRACTION_FOR_RELIABILITY * n)) + 1
    df = _scenario_df(n=n, n_clamped=min(n_clamped, n))
    out = metrics_for_successful_physics(
        galaxy_name="X",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-11,
        galactic_latitude_degrees=45.0,
        morphology_class="dummy",
    )
    assert out["density_audit_passed"] is False
    assert out["is_reliable"] is False


def test_metrics_fails_when_stellar_dominance_fraction_above_threshold():
    n = 10
    n_stellar = int(np.ceil(RELIABILITY_CUTS.MAX_STELLAR_DOMINANCE_FRACTION_FOR_RELIABILITY * n)) + 1
    df = _scenario_df(n=n, n_stellar_dominated=min(n_stellar, n))
    out = metrics_for_successful_physics(
        galaxy_name="X",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-11,
        galactic_latitude_degrees=45.0,
        morphology_class="dummy",
    )
    assert out["stellar_dominance_audit_passed"] is False
    assert out["is_reliable"] is False
