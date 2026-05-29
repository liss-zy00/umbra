"""
Diagnostics added to ``PhysicsEngine`` in the high-relevance refinement pass:
μG plausibility, plasma β, Alfvén velocity, RAR columns, physical-Δx variant, and
inclination-systematic propagation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.config.thresholds import RELIABILITY_CUTS
from src.models.galaxy import Galaxy
from src.physics import PhysicsEngine, TESLA_PER_MICROGAUSS


def _base_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "radius_kpc": [1.0, 2.0, 3.0],
            "observed_velocity_kms": [80.0, 90.0, 95.0],
            "baryonic_velocity_kms": [50.0, 55.0, 60.0],
            "gas_velocity_kms": [40.0, 45.0, 48.0],
            "disk_velocity_kms": [30.0, 33.0, 36.0],
            "bulge_velocity_kms": [0.0, 0.0, 0.0],
        }
    )


def test_force_gap_round_trip_zero_when_v_obs_equals_v_bar():
    """If V_obs == V_bar at every radius, the gap must collapse to zero (no DM signal)."""
    df = _base_df()
    df["observed_velocity_kms"] = df["baryonic_velocity_kms"].copy()
    g = Galaxy(name="RoundTrip", rotation_dataframe=df)
    PhysicsEngine(g).calculate_force_gap()
    gap = g.rotation_dataframe["acceleration_gap"].to_numpy(dtype=float)
    assert np.allclose(gap, 0.0)


def test_rar_columns_populated_after_force_gap():
    df = _base_df()
    g = Galaxy(name="RAR", rotation_dataframe=df)
    PhysicsEngine(g).calculate_force_gap()
    out = g.rotation_dataframe
    assert "rar_predicted_acceleration" in out.columns
    assert "rar_residual_log10" in out.columns
    assert np.all(np.isfinite(out["rar_predicted_acceleration"]))


def test_inclination_systematic_propagates_into_centripetal_band():
    """When an inclination systematic is on V_obs, the centripetal acceleration band must scale."""
    df = _base_df()
    # 5% systematic on V_obs ⇒ 10% on V_obs² ⇒ 10% systematic on a_c.
    df["observed_velocity_systematic_kms"] = 0.05 * df["observed_velocity_kms"]
    g = Galaxy(name="Tilt", rotation_dataframe=df)
    PhysicsEngine(g).calculate_force_gap()
    out = g.rotation_dataframe
    rel = (
        out["centripetal_acceleration_systematic"].to_numpy(dtype=float)
        / out["centripetal_acceleration"].to_numpy(dtype=float)
    )
    assert np.allclose(rel, 0.10, rtol=1e-6)


def test_lorentz_force_emits_microgauss_beta_and_alfven():
    df = _base_df()
    g = Galaxy(name="Sanity", rotation_dataframe=df)
    eng = PhysicsEngine(g)
    eng.calculate_force_gap()
    summary = eng.calculate_magnetic_pressure_acceleration(
        stokes_amp_center=2e-4,
        stokes_amp_neighbors_mean=1e-4,
        scenario_gradient_scale_kpc=0.1,
        scenario_ionization_fraction=0.01,
        stokes_amp_reference=1e-4,
    )
    assert summary["status"] == "ok"
    assert summary["magnetic_field_microgauss"] is not None
    # Per-row μG must equal Tesla / 1e-10.
    out = g.rotation_dataframe
    assert np.isclose(
        out["magnetic_field_microgauss"].iloc[0],
        out["magnetic_field_tesla"].iloc[0] / TESLA_PER_MICROGAUSS,
    )
    assert "plasma_beta_thermal" in out.columns
    assert "alfven_velocity_ms" in out.columns
    # When the field is in the plausibility band, the summary flag must be True.
    if (
        RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MIN_MICROGAUSS
        <= summary["magnetic_field_microgauss"]
        <= RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MAX_MICROGAUSS
    ):
        assert summary["magnetic_field_within_plausibility_band"] is True


def test_lorentz_force_physical_dx_column_is_finite_when_separation_provided():
    df = _base_df()
    g = Galaxy(name="PhysDx", rotation_dataframe=df)
    eng = PhysicsEngine(g)
    eng.calculate_force_gap()
    summary = eng.calculate_magnetic_pressure_acceleration(
        stokes_amp_center=2e-4,
        stokes_amp_neighbors_mean=1e-4,
        scenario_gradient_scale_kpc=0.1,
        scenario_ionization_fraction=0.01,
        stokes_amp_reference=1e-4,
        physical_neighbor_separation_kpc=5.0,
    )
    out = g.rotation_dataframe
    assert "magnetic_pressure_acceleration_physical_dx" in out.columns
    physical_arr = out["magnetic_pressure_acceleration_physical_dx"].to_numpy(dtype=float)
    assert np.all(np.isfinite(physical_arr))
    # The two accelerations must differ by exactly the ratio of the two Δx values.
    standard_arr = out["magnetic_pressure_acceleration"].to_numpy(dtype=float)
    ratio = physical_arr / standard_arr
    assert np.allclose(ratio, 0.1 / 5.0, rtol=1e-12)
    assert summary["physical_gradient_scale_kpc"] == 5.0


def test_engine_summary_exposes_new_refinement_keys():
    """Phase-1 refinements must expose the new diagnostic keys in the summary dict."""
    df = _base_df()
    g = Galaxy(name="Keys", rotation_dataframe=df, distance_mpc=10.0)
    eng = PhysicsEngine(g)
    eng.calculate_force_gap()
    summary = eng.calculate_magnetic_pressure_acceleration(
        stokes_amp_center=2e-4,
        stokes_amp_neighbors_mean=1e-4,
        scenario_gradient_scale_kpc=0.1,
        scenario_ionization_fraction=0.01,
        stokes_amp_reference=1e-4,
        physical_neighbor_separation_kpc=5.0,
        morphology_class="Dwarf Irregulars (Im)",
    )
    for key in (
        "accel_pressure_only_m_s2",
        "accel_full_lorentz_m_s2",
        "tension_term_used",
        "deprojected_gradient_used",
        "audit_flags",
    ):
        assert key in summary, f"Missing key in engine summary: {key}"
    out = g.rotation_dataframe
    for col in ("accel_pressure_only_m_s2", "accel_full_lorentz_m_s2", "tension_acceleration_m_s2"):
        assert col in out.columns
    # Irregular morphology: tension is zero and the engine flags ``warn_pure_pressure``.
    assert "warn_pure_pressure" in summary["audit_flags"]
    assert summary["tension_term_used"] is False


def test_stellar_dominated_shells_are_not_reliable():
    """Shells flagged stellar-dominated must drop out of ``reliable_measurement``."""
    # V_disk dominates so stellar fraction is ~1.
    df = _base_df()
    df["disk_velocity_kms"] = [200.0, 200.0, 200.0]
    df["bulge_velocity_kms"] = [50.0, 50.0, 50.0]
    df["gas_velocity_kms"] = [10.0, 10.0, 10.0]
    df["baryonic_velocity_kms"] = np.sqrt(
        df["gas_velocity_kms"] ** 2 + df["disk_velocity_kms"] ** 2 + df["bulge_velocity_kms"] ** 2
    )
    g = Galaxy(name="StellarDom", rotation_dataframe=df)
    eng = PhysicsEngine(g)
    eng.calculate_force_gap()
    eng.calculate_magnetic_pressure_acceleration(2e-4, 1e-4, 0.1, 0.01, stokes_amp_reference=1e-4)
    out = g.rotation_dataframe
    assert np.all(out["stellar_dominated"])
    assert not np.any(out["reliable_measurement"])
