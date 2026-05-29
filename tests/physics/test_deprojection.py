"""
Disk-frame deprojected gradient estimator
=========================================

The engine replaces the scalar ``Δu_B / Δx`` finite-difference on the physical-Δx path
by a planar least-squares fit ``u_B,i = u_B,0 + g_x · Δx + g_y · Δy`` over the Planck
neighbour ring (see :func:`src.physics.magnetic_pressure_gradient.disk_plane_pressure_gradient`).

The tests exercise three regimes:

* **Closed-form correctness** — a synthetic ring with a pure ``+x`` gradient must
  reproduce the slope exactly, both in the sky frame (PA=0) and after a 90° rotation
  (PA=90° projects ``+x_sky`` onto the disk minor axis, which is then deprojected by
  ``1/cos i``).
* **Axisymmetric fallback** — when ``position_angle_deg=None``, the function still
  returns a finite magnitude (no crash) and applies the ``cos i`` rescaling.
* **Degenerate input** — fewer than two finite neighbours collapses the fit to NaN
  without raising.
"""

from __future__ import annotations

import numpy as np

from src.physics.magnetic_pressure_gradient import disk_plane_pressure_gradient


def _ring(n: int = 8, radius_m: float = 1.0e19) -> tuple[np.ndarray, np.ndarray]:
    """Symmetric ring of ``n`` points at fixed angular separation in (dx, dy) metres."""
    angles = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    dx = radius_m * np.cos(angles)
    dy = radius_m * np.sin(angles)
    return dx, dy


def test_planar_gradient_recovers_pure_x_slope_without_pa():
    """A field ``u_B = u_0 + g · Δx`` is recovered exactly by the LSQ fit (axisymmetric)."""
    dx, dy = _ring()
    g_true = 1.5e-23
    u_b_center = 1.0e-13
    u_b_ring = u_b_center + g_true * dx
    grad_mag, used_pa = disk_plane_pressure_gradient(
        u_b_per_neighbour=u_b_ring,
        u_b_center=u_b_center,
        dx_east_m=dx,
        dy_north_m=dy,
        position_angle_deg=None,
        inclination_deg=None,
    )
    assert used_pa is False
    assert np.isclose(grad_mag, abs(g_true), rtol=1.0e-9)


def test_planar_gradient_pa_zero_aligns_with_sky_north():
    """PA=0 means major-axis on sky-north; the fit must recover the ``+y`` slope cleanly."""
    dx, dy = _ring()
    g_true = 2.5e-23
    u_b_center = 5.0e-14
    u_b_ring = u_b_center + g_true * dy
    grad_mag, used_pa = disk_plane_pressure_gradient(
        u_b_per_neighbour=u_b_ring,
        u_b_center=u_b_center,
        dx_east_m=dx,
        dy_north_m=dy,
        position_angle_deg=0.0,
        inclination_deg=30.0,
    )
    assert used_pa is True
    # PA = 0: major axis = +y_sky, minor axis = -x_sky. The +y component sits entirely on
    # the disk major axis (no inclination stretch), so the fitted magnitude equals g_true.
    assert np.isclose(grad_mag, abs(g_true), rtol=1.0e-9)


def test_planar_gradient_pa_ninety_stretches_minor_axis():
    """PA=90 puts the gradient on the minor axis; the magnitude is rescaled by ``1/cos i``."""
    dx, dy = _ring()
    g_true = 2.0e-23
    incl_deg = 60.0
    u_b_center = 4.0e-14
    u_b_ring = u_b_center + g_true * dy
    grad_mag, used_pa = disk_plane_pressure_gradient(
        u_b_per_neighbour=u_b_ring,
        u_b_center=u_b_center,
        dx_east_m=dx,
        dy_north_m=dy,
        position_angle_deg=90.0,
        inclination_deg=incl_deg,
    )
    assert used_pa is True
    expected = abs(g_true) / float(np.cos(np.deg2rad(incl_deg)))
    assert np.isclose(grad_mag, expected, rtol=1.0e-6)


def test_planar_gradient_degenerate_input_returns_nan():
    """Fewer than two informative neighbours collapses the fit to NaN."""
    dx = np.array([1.0e19, np.nan])
    dy = np.array([0.0, np.nan])
    u_b = np.array([1.0e-13, np.nan])
    grad_mag, used_pa = disk_plane_pressure_gradient(
        u_b_per_neighbour=u_b,
        u_b_center=1.0e-13,
        dx_east_m=dx,
        dy_north_m=dy,
        position_angle_deg=None,
        inclination_deg=None,
    )
    assert not np.isfinite(grad_mag)
    assert used_pa is False


def test_planar_gradient_axisymmetric_fallback_applies_cos_inclination():
    """No PA, but inclination provided: the magnitude must be rescaled by ``cos i``."""
    dx, dy = _ring()
    g_true = 3.0e-23
    incl_deg = 45.0
    u_b_center = 2.0e-14
    u_b_ring = u_b_center + g_true * dx
    grad_mag, used_pa = disk_plane_pressure_gradient(
        u_b_per_neighbour=u_b_ring,
        u_b_center=u_b_center,
        dx_east_m=dx,
        dy_north_m=dy,
        position_angle_deg=None,
        inclination_deg=incl_deg,
    )
    assert used_pa is False
    expected = abs(g_true) * float(np.cos(np.deg2rad(incl_deg)))
    assert np.isclose(grad_mag, expected, rtol=1.0e-9)
