"""
Magnetic energy-density gradient helpers and scalar finite-difference estimator.

Planar LSQ on HEALPix neighbour rings estimates ``||∇u_B||`` in the disk frame when PA and
inclination are known; otherwise the sky-plane fit is rescaled by ``cos i`` (axisymmetric
approximation). Scalar ``Δu_B/Δx`` is retained as the exploratory scenario-scale path.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u

from src.config.thresholds import DISK_PLANE_COS_INCLINATION_NUMERICAL_FLOOR

__all__ = [
    "MagneticPressureGradientEstimator",
    "disk_plane_pressure_gradient",
    "sky_offsets_metres",
]


def sky_offsets_metres(
    *,
    center_theta_rad: float,
    center_phi_rad: float,
    neighbour_theta_rad: np.ndarray,
    neighbour_phi_rad: np.ndarray,
    distance_mpc: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Small-angle tangent-plane offsets of HEALPix neighbours from the central pixel.

    Returns ``(dx_east_m, dy_north_m)`` in SI metres. The east direction is taken as
    ``+φ`` (RA-like) and north as decreasing colatitude (``-θ``). The conversion uses the
    galaxy luminosity distance ``D_galaxy`` so the tangent-plane displacement matches the
    physical Planck-pixel projection on the disk.
    """
    d_phi = np.asarray(neighbour_phi_rad, dtype=float) - float(center_phi_rad)
    # Wrap Δφ to [-π, +π] so the pixel ring around the wrap-around longitude does not
    # introduce a 2π jump in the linear LSQ system.
    d_phi = (d_phi + np.pi) % (2.0 * np.pi) - np.pi
    dx_east_rad = np.sin(float(center_theta_rad)) * d_phi
    dy_north_rad = float(center_theta_rad) - np.asarray(neighbour_theta_rad, dtype=float)

    one_mpc_in_m = float((1.0 * u.Mpc).to_value(u.m))
    scale_m = float(distance_mpc) * one_mpc_in_m
    return dx_east_rad * scale_m, dy_north_rad * scale_m


def disk_plane_pressure_gradient(
    *,
    u_b_per_neighbour: np.ndarray,
    u_b_center: float,
    dx_east_m: np.ndarray,
    dy_north_m: np.ndarray,
    position_angle_deg: float | None,
    inclination_deg: float | None,
) -> tuple[float, bool]:
    """
    Planar least-squares estimate of ``||∇u_B||`` in the galaxy disk frame.

    When both ``position_angle_deg`` and ``inclination_deg`` are finite, the neighbour
    offsets are rotated so the disk major axis aligns with the local ``+x_disk`` axis and
    the minor-axis component is stretched by ``1/cos i`` to compensate for the projection
    onto the sky. A planar least-squares fit
    ``u_B,i = u_B,0 + g_x · Δx_disk + g_y · Δy_disk`` then returns the gradient
    magnitude ``√(g_x² + g_y²)`` in J/m⁴.

    When PA is missing the fit runs in the raw sky frame and the magnitude is rescaled by
    ``cos i`` (linear approximation of the disk-plane gradient under axisymmetric assumption);
    callers are expected to flag ``warn_no_pa`` in that case.

    The ``1/cos i`` stretch uses ``max(cos i, DISK_PLANE_COS_INCLINATION_NUMERICAL_FLOOR)`` so
    the LSQ system stays well-conditioned as ``i → 90°`` (geometric singularity of edge-on
    projection, not a physical preference).

    Returns
    -------
    tuple[float, bool]
        ``(gradient_magnitude_j_m4, used_pa)``. ``gradient_magnitude_j_m4`` is ``nan`` when
        the LSQ system is degenerate (fewer than two informative neighbours).
    """
    floor = float(DISK_PLANE_COS_INCLINATION_NUMERICAL_FLOOR)
    u_neighbours = np.asarray(u_b_per_neighbour, dtype=float)
    dx = np.asarray(dx_east_m, dtype=float)
    dy = np.asarray(dy_north_m, dtype=float)
    finite = np.isfinite(u_neighbours) & np.isfinite(dx) & np.isfinite(dy)
    if int(np.count_nonzero(finite)) < 2 or not np.isfinite(u_b_center):
        return float("nan"), False
    dx = dx[finite]
    dy = dy[finite]
    u_neighbours = u_neighbours[finite]

    # Anchor the LSQ system with the central pixel (Δx=Δy=0, u_B=u_B,center). It pins the
    # intercept to a physically meaningful value and improves conditioning at small N.
    dx_full = np.concatenate(([0.0], dx))
    dy_full = np.concatenate(([0.0], dy))
    u_full = np.concatenate(([float(u_b_center)], u_neighbours))

    use_pa = (
        position_angle_deg is not None
        and inclination_deg is not None
        and np.isfinite(float(position_angle_deg))
        and np.isfinite(float(inclination_deg))
        and 0.0 < float(inclination_deg) < 90.0
    )

    if use_pa:
        pa_rad = float(np.deg2rad(float(position_angle_deg)))
        cos_pa, sin_pa = np.cos(pa_rad), np.sin(pa_rad)
        # PA is measured east of north on the sky. The major-axis unit vector in
        # (east, north) coordinates is (sin PA, cos PA); the minor axis is the orthogonal
        # (-cos PA, sin PA). Project onto disk axes, then deproject the minor leg.
        dx_major = dx_full * sin_pa + dy_full * cos_pa
        dy_minor_sky = -dx_full * cos_pa + dy_full * sin_pa
        cos_i = float(np.cos(np.deg2rad(float(inclination_deg))))
        cos_i = cos_i if cos_i > floor else floor
        dy_minor_disk = dy_minor_sky / cos_i
        coords = np.column_stack([np.ones_like(dx_major), dx_major, dy_minor_disk])
    else:
        coords = np.column_stack([np.ones_like(dx_full), dx_full, dy_full])

    try:
        coeffs, *_ = np.linalg.lstsq(coords, u_full, rcond=None)
    except np.linalg.LinAlgError:
        return float("nan"), False
    if coeffs.size < 3 or not np.all(np.isfinite(coeffs)):
        return float("nan"), False
    g_x, g_y = float(coeffs[1]), float(coeffs[2])
    grad_mag = float(np.hypot(g_x, g_y))

    if not use_pa:
        cos_i = (
            float(np.cos(np.deg2rad(float(inclination_deg))))
            if inclination_deg is not None and np.isfinite(float(inclination_deg))
            and 0.0 < float(inclination_deg) < 90.0
            else 1.0
        )
        grad_mag = grad_mag * cos_i

    return grad_mag, use_pa


class MagneticPressureGradientEstimator:
    """
    Finite-difference estimate of ``∇u_B`` along one axis using sky-center vs sky-neighbor ``u_B``.

    ``pressure_gradient_j_m4`` has SI units J/m⁴ (same as Pa/m). MHD acceleration uses ``a = -(∇P_B)/ρ_ion``.
    """

    @staticmethod
    def delta_energy_density(u_center: float, u_neighbors: float) -> float:
        """``Δu_B = u_B(neighbors) - u_B(center)`` in J/m³."""
        return float(u_neighbors - u_center)

    @staticmethod
    def pressure_gradient(delta_u_j_m3: float, gradient_scale_m: float) -> float:
        """Scalar approximation ``∇u_B ≈ Δu_B / Δx``."""
        return float(delta_u_j_m3 / gradient_scale_m)

    @staticmethod
    def mhd_acceleration(pressure_gradient_j_m4: float, effective_ionized_density_kg_m3: np.ndarray) -> np.ndarray:
        """``a_B = -(∇u_B)/ρ_ion`` element-wise."""
        return -pressure_gradient_j_m4 / np.asarray(effective_ionized_density_kg_m3, dtype=float)
