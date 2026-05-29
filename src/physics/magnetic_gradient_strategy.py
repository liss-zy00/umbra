"""
Strategies for the physical-scale magnetic energy-density gradient used in ``a_B = -(∇u_B)/ρ``.

The scenario column always uses scalar ``Δu_B / Δx_scenario``. The physical-Δx column can use
either the same finite-difference on the true pixel separation or a disk-plane LSQ fit when a
full HEALPix neighbour ring is available.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from src.config.thresholds import AUDIT_CODES

from .magnetic_field_calibration import MagneticFieldCalibrator
from .magnetic_pressure_gradient import (
    MagneticPressureGradientEstimator,
    disk_plane_pressure_gradient,
    sky_offsets_metres,
)


@dataclass(frozen=True)
class PhysicalGradientInputs:
    """Inputs shared by physical-gradient strategies (SI / radians / dimensionless Stokes)."""

    delta_magnetic_energy_density_j_m3: float
    physical_gradient_scale_m: float
    magnetic_energy_density_center_j_m3: float
    stokes_amp_reference: float
    neighbour_ring: list[dict[str, float]] | None
    center_theta_rad: float | None
    center_phi_rad: float | None
    distance_mpc: float
    position_angle_deg: float | None
    inclination_deg: float | None
    field_calibrator: MagneticFieldCalibrator
    mu0_si: float


@dataclass(frozen=True)
class PhysicalGradientOutcome:
    """Scalar ``∇u_B`` (J/m⁴) for the physical-Δx MHD path and provenance."""

    physical_pressure_gradient_j_m4: float
    deprojected_gradient_used: bool
    audit_flags: tuple[str, ...]


class MagneticPhysicalGradientStrategy(Protocol):
    def compute(self, inp: PhysicalGradientInputs) -> PhysicalGradientOutcome:
        ...


class ScalarDeltaOverPhysicalDxStrategy:
    """Finite-difference ``∇u_B ≈ Δu_B / Δx_phys`` using mean neighbour vs centre Stokes."""

    def __init__(self, pressure_gradient_estimator: MagneticPressureGradientEstimator) -> None:
        self._pressure = pressure_gradient_estimator

    def compute(self, inp: PhysicalGradientInputs) -> PhysicalGradientOutcome:
        g = self._pressure.pressure_gradient(
            inp.delta_magnetic_energy_density_j_m3,
            inp.physical_gradient_scale_m,
        )
        return PhysicalGradientOutcome(
            physical_pressure_gradient_j_m4=float(g),
            deprojected_gradient_used=False,
            audit_flags=(),
        )


class DiskPlaneRingLeastSquaresStrategy:
    """
    Planar LSQ gradient of ``u_B`` on the sky ring, rotated/deprojected into the disk frame.

    Falls back to :class:`ScalarDeltaOverPhysicalDxStrategy` when the LSQ system is degenerate.
    """

    def __init__(
        self,
        pressure_gradient_estimator: MagneticPressureGradientEstimator,
        *,
        scalar_fallback: ScalarDeltaOverPhysicalDxStrategy,
    ) -> None:
        self._pressure = pressure_gradient_estimator
        self._scalar_fallback = scalar_fallback

    def compute(self, inp: PhysicalGradientInputs) -> PhysicalGradientOutcome:
        ring = inp.neighbour_ring
        if ring is None or len(ring) < 2:
            return self._scalar_fallback.compute(inp)
        theta_n = np.array([float(p.get("theta_rad", np.nan)) for p in ring], dtype=float)
        phi_n = np.array([float(p.get("phi_rad", np.nan)) for p in ring], dtype=float)
        q_n = np.array([float(p.get("q", np.nan)) for p in ring], dtype=float)
        u_n = np.array([float(p.get("u", np.nan)) for p in ring], dtype=float)
        amp_n = np.hypot(q_n, u_n)
        b_per_neighbour_t = np.array(
            [inp.field_calibrator.field_tesla(float(a), inp.stokes_amp_reference) for a in amp_n],
            dtype=float,
        )
        u_b_per_neighbour = (b_per_neighbour_t**2) / (2.0 * inp.mu0_si)
        if inp.center_theta_rad is None or inp.center_phi_rad is None:
            return self._scalar_fallback.compute(inp)
        dx_east_m, dy_north_m = sky_offsets_metres(
            center_theta_rad=float(inp.center_theta_rad),
            center_phi_rad=float(inp.center_phi_rad),
            neighbour_theta_rad=theta_n,
            neighbour_phi_rad=phi_n,
            distance_mpc=float(inp.distance_mpc),
        )
        gradient_magnitude, used_pa = disk_plane_pressure_gradient(
            u_b_per_neighbour=u_b_per_neighbour,
            u_b_center=float(inp.magnetic_energy_density_center_j_m3),
            dx_east_m=dx_east_m,
            dy_north_m=dy_north_m,
            position_angle_deg=inp.position_angle_deg,
            inclination_deg=inp.inclination_deg,
        )
        if not np.isfinite(gradient_magnitude):
            return self._scalar_fallback.compute(inp)
        du = inp.delta_magnetic_energy_density_j_m3
        signed = (
            float(np.sign(du) * gradient_magnitude)
            if du != 0.0
            else float(gradient_magnitude)
        )
        flags: list[str] = []
        if not used_pa:
            flags.append(AUDIT_CODES.WARN_NO_POSITION_ANGLE)
        return PhysicalGradientOutcome(
            physical_pressure_gradient_j_m4=signed,
            deprojected_gradient_used=True,
            audit_flags=tuple(flags),
        )


def resolve_physical_gradient_strategy(
    *,
    neighbour_ring: list[dict[str, float]] | None,
    center_theta_rad: float | None,
    center_phi_rad: float | None,
    distance_mpc: float,
    pressure_gradient_estimator: MagneticPressureGradientEstimator,
) -> MagneticPhysicalGradientStrategy:
    """
    Pick LSQ ring strategy when geometry is complete; otherwise scalar Δu_B/Δx_phys.
    """
    scalar = ScalarDeltaOverPhysicalDxStrategy(pressure_gradient_estimator)
    ring_ok = (
        neighbour_ring is not None
        and len(neighbour_ring) >= 2
        and center_theta_rad is not None
        and center_phi_rad is not None
        and np.isfinite(distance_mpc)
        and distance_mpc > 0
    )
    if ring_ok:
        return DiskPlaneRingLeastSquaresStrategy(pressure_gradient_estimator, scalar_fallback=scalar)
    return scalar
