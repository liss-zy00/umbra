import numpy as np
import pandas as pd
from astropy import constants as const
from astropy import units as u

from src.models.galaxy import Galaxy
from src.physics.engine import MU0, ENERGY_DENSITY_CMB, PhysicsEngine


def test_calculate_force_gap():
    df = pd.DataFrame(
        {
            "radius_kpc": [1.0],
            "observed_velocity_kms": [100.0],
            "baryonic_velocity_kms": [50.0],
            "gas_velocity_kms": [10.0],
        }
    )
    galaxy = Galaxy(name="TestGal", rotation_dataframe=df)
    engine = PhysicsEngine(galaxy)
    engine.calculate_force_gap()

    res_df = engine.galaxy.rotation_dataframe
    r_m = 1.0 * u.kpc.to(u.m)
    v_obs_ms = 100000.0
    v_bar_ms = 50000.0

    expected_centripetal = (v_obs_ms**2) / r_m
    expected_gap = expected_centripetal - (v_bar_ms**2) / r_m

    assert np.isclose(res_df["centripetal_acceleration"].iloc[0], expected_centripetal)
    assert np.isclose(res_df["acceleration_gap"].iloc[0], expected_gap)


def test_calculate_magnetic_pressure_acceleration_values():
    df = pd.DataFrame(
        {
            "radius_kpc": [1.0],
            "observed_velocity_kms": [100.0],
            "baryonic_velocity_kms": [50.0],
            "gas_velocity_kms": [10.0],
        }
    )
    galaxy = Galaxy(name="TestGal", rotation_dataframe=df)
    engine = PhysicsEngine(galaxy)
    engine.calculate_force_gap()

    stokes_center = 2e-4
    stokes_neighbors = 1e-4
    scale_kpc = 0.1
    ion_frac = 0.01

    result = engine.calculate_magnetic_pressure_acceleration(
        stokes_amp_center=stokes_center,
        stokes_amp_neighbors_mean=stokes_neighbors,
        scenario_gradient_scale_kpc=scale_kpc,
        scenario_ionization_fraction=ion_frac,
        stokes_amp_reference=1e-4,
    )

    res_df = engine.galaxy.rotation_dataframe

    b_cmb = np.sqrt(2.0 * MU0 * ENERGY_DENSITY_CMB)
    expected_b_center = b_cmb * (stokes_center / 1e-4)
    assert np.isclose(res_df["magnetic_field_tesla"].iloc[0], expected_b_center)

    expected_energy_center = (expected_b_center**2) / (2.0 * MU0)
    assert np.isclose(res_df["magnetic_energy_density_center"].iloc[0], expected_energy_center)

    expected_b_neighbors = b_cmb * (stokes_neighbors / 1e-4)
    expected_energy_neigh = (expected_b_neighbors**2) / (2.0 * MU0)
    delta_u = expected_energy_neigh - expected_energy_center
    delta_x = scale_kpc * u.kpc.to(u.m)
    expected_grad = delta_u / delta_x
    assert np.isclose(res_df["magnetic_pressure_gradient"].iloc[0], expected_grad)

    rho_ion = res_df["effective_ionized_mass_density"].iloc[0]
    expected_accel = -expected_grad / rho_ion
    assert np.isclose(res_df["magnetic_pressure_acceleration"].iloc[0], expected_accel)
    assert np.isfinite(expected_accel)
    assert result["status"] == "ok"
    assert result["magnetic_pressure_acceleration"] is not None

    # Delta x audit: scenario_gradient_scale_m must match astropy's u.kpc.to(u.m) value.
    assert np.isclose(
        res_df["scenario_gradient_scale_m"].iloc[0],
        scale_kpc * u.kpc.to(u.m),
    )
    assert np.isclose(result["scenario_gradient_scale_m"], scale_kpc * u.kpc.to(u.m))
    # No physical Delta x supplied -> summary value is None and the column is NaN.
    assert result["physical_gradient_scale_m"] is None
    assert np.isnan(res_df["physical_gradient_scale_m"].iloc[0])


def test_physical_gradient_scale_m_audit_column():
    """When a physical Delta x is supplied, the SI-broadcast column matches astropy's conversion."""
    df = pd.DataFrame(
        {
            "radius_kpc": [1.0],
            "observed_velocity_kms": [100.0],
            "baryonic_velocity_kms": [50.0],
            "gas_velocity_kms": [10.0],
        }
    )
    galaxy = Galaxy(name="TestGal", rotation_dataframe=df)
    engine = PhysicsEngine(galaxy)
    engine.calculate_force_gap()
    physical_kpc = 0.42
    result = engine.calculate_magnetic_pressure_acceleration(
        stokes_amp_center=2e-4,
        stokes_amp_neighbors_mean=1e-4,
        scenario_gradient_scale_kpc=0.1,
        scenario_ionization_fraction=0.01,
        stokes_amp_reference=1e-4,
        physical_neighbor_separation_kpc=physical_kpc,
    )
    assert np.isclose(
        engine.galaxy.rotation_dataframe["physical_gradient_scale_m"].iloc[0],
        physical_kpc * u.kpc.to(u.m),
    )
    assert np.isclose(result["physical_gradient_scale_m"], physical_kpc * u.kpc.to(u.m))


def test_magnetic_pressure_skipped_on_missing_data():
    df = pd.DataFrame(
        {
            "radius_kpc": [1.0],
            "observed_velocity_kms": [100.0],
            "baryonic_velocity_kms": [50.0],
            "gas_velocity_kms": [10.0],
        }
    )
    galaxy = Galaxy(name="TestGal", rotation_dataframe=df)
    engine = PhysicsEngine(galaxy)
    engine.calculate_force_gap()

    result = engine.calculate_magnetic_pressure_acceleration(None, None, 0.1, 0.01)

    assert "Skipped" in result["status"]
    assert not bool(engine.galaxy.rotation_dataframe["reliable_measurement"].iloc[0])
    assert np.isnan(engine.galaxy.rotation_dataframe["magnetic_pressure_acceleration"].iloc[0])


def test_gas_mass_enclosed_monotonic_for_increasing_rotation_curve():
    """Rough check: faster rotation at larger radius ⇒ non-decreasing enclosed gas mass."""
    radii_kpc = np.linspace(1.0, 5.0, 8)
    v_gas = 20.0 + 3.0 * radii_kpc
    df = pd.DataFrame(
        {
            "radius_kpc": radii_kpc,
            "observed_velocity_kms": v_gas + 50.0,
            "baryonic_velocity_kms": np.sqrt(v_gas**2 + 10.0**2),
            "gas_velocity_kms": v_gas,
        }
    )
    g = Galaxy(name="Mono", rotation_dataframe=df)
    engine = PhysicsEngine(g)
    engine.calculate_force_gap()
    engine.calculate_magnetic_pressure_acceleration(1e-4, 1e-4, 0.1, 0.01)

    r_m = df["radius_kpc"].to_numpy(dtype=float) * u.kpc.to(u.m)
    v_ms = v_gas * 1000.0
    g_si = const.G.to_value(u.m**3 / (u.kg * u.s**2))
    m_enc = (v_ms**2 * r_m) / g_si
    assert np.all(np.diff(m_enc) >= -1e-30)
