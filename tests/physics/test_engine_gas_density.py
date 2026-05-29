import numpy as np
from astropy import units as u

from src.physics.engine import GasDensityEstimator, MIN_GAS_MASS_DENSITY


def test_shell_density_positive_and_finite():
    est = GasDensityEstimator()
    radius_kpc = np.array([1.0, 2.0, 3.0], dtype=float)
    radius_m = radius_kpc * u.kpc.to(u.m)
    v_gas_kms = np.array([30.0, 35.0, 40.0], dtype=float)
    (
        rho,
        thickness_kpc,
        clamped_mask,
        stellar_fraction,
        stellar_mask,
        residual_kg,
    ) = est.shell_gas_mass_density(radius_m, v_gas_kms, radius_kpc)
    assert np.all(np.isfinite(rho))
    assert np.all(rho > 0)
    assert thickness_kpc >= 0.2
    assert clamped_mask.shape == rho.shape
    assert clamped_mask.dtype == bool
    # Without stellar inputs, fraction is identically zero and nothing is flagged stellar-dominated.
    assert np.all(stellar_fraction == 0.0)
    assert not np.any(stellar_mask)
    # No baryonic decomposition supplied → residual cross-check is identically zero by design.
    assert residual_kg.shape == rho.shape
    assert np.allclose(residual_kg, 0.0)


def test_zero_gas_velocity_clamps_all_shells_to_floor():
    """V_gas = 0 implies enclosed mass = 0 → every shell density falls back to the floor.

    The provider must report this honestly via ``clamped_mask`` so downstream code can mark
    the scenario as unreliable instead of treating the floor as a real measurement.
    """
    est = GasDensityEstimator()
    radius_kpc = np.array([1.0, 2.0, 3.0], dtype=float)
    radius_m = radius_kpc * u.kpc.to(u.m)
    v_gas_kms = np.zeros_like(radius_kpc)
    rho, _thickness, clamped_mask, *_rest = est.shell_gas_mass_density(
        radius_m, v_gas_kms, radius_kpc
    )
    assert np.all(clamped_mask)
    assert np.allclose(rho, MIN_GAS_MASS_DENSITY)


def test_stellar_dominance_flag_picks_up_disk_dominated_shells():
    """When V_disk² ≫ V_gas² we expect the stellar fraction near 1 and the mask to fire."""
    est = GasDensityEstimator()
    radius_kpc = np.array([1.0, 2.0, 3.0], dtype=float)
    radius_m = radius_kpc * u.kpc.to(u.m)
    v_gas_kms = np.array([10.0, 10.0, 10.0], dtype=float)
    v_disk_kms = np.array([100.0, 100.0, 100.0], dtype=float)
    _rho, _thickness, _clamp, stellar_fraction, stellar_mask, _residual = est.shell_gas_mass_density(
        radius_m, v_gas_kms, radius_kpc, disk_velocity_kms=v_disk_kms
    )
    assert np.all(stellar_fraction > 0.95)
    assert np.all(stellar_mask)


def test_baryonic_residual_cross_check_is_consistent():
    """When V_bar = sqrt(V_gas² + V_disk² + V_bulge²) the two enclosed-mass paths must agree."""
    est = GasDensityEstimator()
    radius_kpc = np.array([1.0, 2.0, 3.0], dtype=float)
    radius_m = radius_kpc * u.kpc.to(u.m)
    v_gas_kms = np.array([30.0, 40.0, 45.0], dtype=float)
    v_disk_kms = np.array([10.0, 20.0, 25.0], dtype=float)
    v_bulge_kms = np.array([5.0, 10.0, 15.0], dtype=float)
    v_bar_kms = np.sqrt(v_gas_kms**2 + v_disk_kms**2 + v_bulge_kms**2)
    _rho, _thickness, _clamp, _stellar_fraction, _stellar_mask, residual_kg = (
        est.shell_gas_mass_density(
            radius_m,
            v_gas_kms,
            radius_kpc,
            disk_velocity_kms=v_disk_kms,
            bulge_velocity_kms=v_bulge_kms,
            baryonic_velocity_kms=v_bar_kms,
        )
    )
    # The residual must be ~0 (catalog-consistent decomposition).
    assert np.allclose(residual_kg, 0.0, atol=1e-6 * np.abs(residual_kg).max() + 1e30)


def test_observed_velocity_strict_stellar_gate_fires():
    """When V_stellar²/V_obs² is high the stricter gate must mark the shell as dominated."""
    est = GasDensityEstimator()
    radius_kpc = np.array([1.0, 2.0, 3.0], dtype=float)
    radius_m = radius_kpc * u.kpc.to(u.m)
    v_gas_kms = np.array([20.0, 22.0, 25.0], dtype=float)
    v_disk_kms = np.array([90.0, 95.0, 100.0], dtype=float)
    v_bulge_kms = np.array([20.0, 20.0, 20.0], dtype=float)
    v_obs_kms = np.array([95.0, 100.0, 105.0], dtype=float)
    _rho, _thickness, _clamp, _frac, stellar_mask, _residual = est.shell_gas_mass_density(
        radius_m,
        v_gas_kms,
        radius_kpc,
        disk_velocity_kms=v_disk_kms,
        bulge_velocity_kms=v_bulge_kms,
        observed_velocity_kms=v_obs_kms,
    )
    assert np.all(stellar_mask)
