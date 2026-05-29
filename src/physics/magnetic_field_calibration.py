"""
Stokes-amplitude to magnetic-field calibration (dimensionless Q,U → Tesla).

The scale ``B ∝ √(2 μ₀ u_ref)`` uses ``u_ref`` as a dimensionful reference energy density
(CMB-scale by default); that is a units anchor, not a causal claim about field origin.
"""

from __future__ import annotations

import numpy as np

from .defaults import DEFAULT_STOKES_AMPLITUDE_REFERENCE, ENERGY_DENSITY_CMB, MU0

__all__ = ["MagneticFieldCalibrator"]


class MagneticFieldCalibrator:
    """
    Single responsibility: convert dimensionless Stokes amplitudes to magnetic field strength (Tesla).

    Uses ``B_scale = √(2 μ₀ u_ref)`` times ``S/S_ref``, where ``u_ref`` defaults to ``ENERGY_DENSITY_CMB``.
    """

    def __init__(self, mu0: float = MU0, reference_energy_density_j_m3: float = ENERGY_DENSITY_CMB) -> None:
        self.mu0 = mu0
        self.reference_energy_density_j_m3 = reference_energy_density_j_m3

    def b_scale_tesla(self) -> float:
        """Constant multiplier √(2 μ₀ u_ref) linking Stokes ratio to field strength."""
        return float(np.sqrt(2.0 * self.mu0 * self.reference_energy_density_j_m3))

    def field_tesla(self, stokes_amp: float, stokes_reference: float) -> float:
        """Return ``B`` (T) for amplitude ``S`` given reference ``S_ref``."""
        ref = (
            float(stokes_reference)
            if np.isfinite(stokes_reference) and stokes_reference > 0
            else DEFAULT_STOKES_AMPLITUDE_REFERENCE
        )
        return self.b_scale_tesla() * (float(stokes_amp) / ref)
