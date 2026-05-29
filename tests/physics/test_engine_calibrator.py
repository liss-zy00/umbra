import numpy as np

from src.physics.engine import ENERGY_DENSITY_CMB, MU0, MagneticFieldCalibrator


def test_field_scales_linearly_with_stokes_amplitude():
    cal = MagneticFieldCalibrator()
    ref = 1e-4
    b1 = cal.field_tesla(1e-4, ref)
    b2 = cal.field_tesla(2e-4, ref)
    assert np.isclose(b2, 2.0 * b1)


def test_b_scale_matches_sqrt_two_mu_zero_u_cmb():
    cal = MagneticFieldCalibrator()
    expected = float(np.sqrt(2.0 * MU0 * ENERGY_DENSITY_CMB))
    assert np.isclose(cal.b_scale_tesla(), expected)
