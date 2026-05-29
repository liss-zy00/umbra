import numpy as np

from src.physics.engine import MagneticPressureGradientEstimator


def test_pressure_gradient_changes_sign_when_neighbor_energy_flips():
    est = MagneticPressureGradientEstimator()
    dx = 100.0
    g_high = est.pressure_gradient(est.delta_energy_density(1.0, 3.0), dx)
    g_low = est.pressure_gradient(est.delta_energy_density(3.0, 1.0), dx)
    assert g_high > 0
    assert g_low < 0


def test_mhd_acceleration_opposes_positive_gradient():
    rho = np.array([1e-24, 2e-24], dtype=float)
    grad = 1e-30
    a = MagneticPressureGradientEstimator.mhd_acceleration(grad, rho)
    assert np.allclose(a, -grad / rho)
