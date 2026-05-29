"""
Magnetic-tension term ``(B · ∇) B / μ₀`` under a toy logarithmic-spiral field geometry.

Background
----------
Under ideal MHD the volume Lorentz force on a magnetised fluid is

.. math::

    \\vec f_{\\mathrm{Lorentz}} \\;=\\; -\\nabla\\!\\left(\\frac{B^2}{2\\mu_0}\\right)
      \\;+\\;\\frac{(\\vec B\\cdot\\nabla)\\vec B}{\\mu_0}.

The pipeline already evaluates the first term (the magnetic-*pressure* gradient) inside
:class:`src.physics.engine.PhysicsEngine`. This module implements the second term — the
magnetic *tension* — under a deliberately schematic assumption: that the disk field traces
a **logarithmic spiral** with pitch angle ``p``, in which case the field-line curvature
radius at galactocentric radius ``r`` collapses to the closed-form expression

.. math::

    R_c(r) \\;=\\; \\frac{r}{\\sin^2 p}.

(Identity for a logarithmic spiral ``r(\\theta) = a\\,e^{b\\theta}`` with ``b = \\tan p``;
see Binney & Tremaine 2008, Ch. 6.) The radial component of the tension term — the only
component that competes with gravity in the disk-plane equilibrium — is then

.. math::

    a_{\\mathrm{tension},r}(r) \\;=\\; -\\cos p \\cdot \\frac{B^2(r)}{\\mu_0\\,\\rho_{\\mathrm{ion}}(r)\\,R_c(r)}.

The sign is negative because the tension pulls the field line back towards the curvature
centre, i.e. inwards in the disk plane.

Caveats — read first
--------------------
This is a *toy* geometry meant to bound the magnitude of the tension contribution, not
a realistic field reconstruction. Irregular and BCD galaxies do not host a coherent
logarithmic spiral; for those morphologies :func:`tension_acceleration_radial` returns a
zero array and the caller is expected to surface ``warn_pure_pressure`` on the audit-flag
list (see :data:`src.config.thresholds.AUDIT_CODES.WARN_PURE_PRESSURE_GEOMETRY`).
"""

from __future__ import annotations

from typing import Iterable

import numpy as np

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS, UNIVERSAL_CONSTANTS


# Hubble-stage tokens (or substrings) that mark morphologies without a coherent spiral
# field geometry. The tension toy collapses to zero for these galaxies; see notes above.
_NON_SPIRAL_MORPHOLOGY_TOKENS: tuple[str, ...] = (
    "irregular",
    "irr",
    "im",
    "bcd",
)


def supports_logarithmic_spiral(morphology_class: str | None) -> bool:
    """
    Return ``True`` when ``morphology_class`` is consistent with a coherent disk spiral.

    The check is intentionally permissive: any morphology that does not match an irregular
    or BCD token is treated as spiral-like. This matches the SPARC catalogue convention
    where late-type spirals (Sa–Sd) all carry a measurable pitch angle in the literature.
    """
    if not morphology_class:
        return False
    label = str(morphology_class).strip().lower()
    if not label:
        return False
    return not any(tok in label for tok in _NON_SPIRAL_MORPHOLOGY_TOKENS)


def curvature_radius_m(radius_m: np.ndarray | Iterable[float], pitch_deg: float) -> np.ndarray:
    """
    Logarithmic-spiral field-line curvature radius ``R_c(r) = r / sin²(p)``.

    Parameters
    ----------
    radius_m:
        Per-shell galactocentric radius in metres (SI).
    pitch_deg:
        Spiral pitch angle in degrees. Must be in ``(0, 90)``; a zero pitch corresponds to
        a purely radial field (no curvature) and the curvature radius diverges.

    Returns
    -------
    np.ndarray
        Curvature radius in metres, same shape as ``radius_m``. NaN where ``radius_m`` is
        non-positive / non-finite or where ``pitch_deg`` is outside the valid open interval.
    """
    r = np.asarray(radius_m, dtype=float)
    if not (0.0 < float(pitch_deg) < 90.0):
        return np.full_like(r, np.nan, dtype=float)
    sin_p = float(np.sin(np.deg2rad(float(pitch_deg))))
    if sin_p <= 0.0:
        return np.full_like(r, np.nan, dtype=float)
    out = np.full_like(r, np.nan, dtype=float)
    valid = np.isfinite(r) & (r > 0.0)
    out[valid] = r[valid] / (sin_p**2)
    return out


def tension_acceleration_radial(
    *,
    radius_m: np.ndarray | Iterable[float],
    magnetic_field_tesla: float | np.ndarray | Iterable[float],
    effective_ionized_mass_density: np.ndarray | Iterable[float],
    morphology_class: str | None = None,
    pitch_deg: float = ASTROPHYSICAL_DEFAULTS.SPIRAL_PITCH_DEG,
    mu0: float = UNIVERSAL_CONSTANTS.MU0_SI,
) -> np.ndarray:
    """
    Radial component of the magnetic-tension acceleration under the toy spiral geometry.

    Parameters
    ----------
    radius_m:
        Per-shell galactocentric radius in metres.
    magnetic_field_tesla:
        Per-shell field magnitude in Tesla. Scalar broadcasts to the radius shape.
    effective_ionized_mass_density:
        Per-shell ionised density ``ρ_ion`` in kg/m³ (the same denominator used by the
        magnetic-pressure-gradient acceleration in :mod:`src.physics.engine`).
    morphology_class:
        SPARC morphology bucket label (e.g. ``"Late-type Spirals (Sc, Scd)"``,
        ``"Dwarf Irregulars (Im)"``). When the label fails the
        :func:`supports_logarithmic_spiral` test the function returns an all-zero array.
    pitch_deg:
        Spiral pitch angle in degrees. Defaults to
        :attr:`src.config.thresholds.AstrophysicalDefaults.SPIRAL_PITCH_DEG`.
    mu0:
        Vacuum permeability (T·m/A). Defaults to
        :attr:`src.config.thresholds.UniversalConstants.MU0_SI`.

    Returns
    -------
    np.ndarray
        Per-shell ``a_tension,r`` in m/s² (negative ⇒ inward), shape matches ``radius_m``.
        Returns zeros for irregular / BCD morphologies and NaN at shells where the inputs
        are non-finite or ``ρ_ion`` is non-positive.
    """
    r = np.asarray(radius_m, dtype=float)
    if not supports_logarithmic_spiral(morphology_class):
        return np.zeros_like(r, dtype=float)

    b = np.broadcast_to(np.asarray(magnetic_field_tesla, dtype=float), r.shape).astype(float, copy=True)
    rho = np.asarray(effective_ionized_mass_density, dtype=float)
    if rho.shape != r.shape:
        rho = np.broadcast_to(rho, r.shape).astype(float, copy=True)

    rc = curvature_radius_m(r, pitch_deg=float(pitch_deg))
    cos_p = float(np.cos(np.deg2rad(float(pitch_deg))))

    out = np.full_like(r, np.nan, dtype=float)
    valid = (
        np.isfinite(r) & (r > 0.0)
        & np.isfinite(b)
        & np.isfinite(rho) & (rho > 0.0)
        & np.isfinite(rc) & (rc > 0.0)
    )
    if not np.any(valid):
        return out
    out[valid] = -cos_p * (b[valid] ** 2) / (mu0 * rho[valid] * rc[valid])
    return out


__all__ = [
    "supports_logarithmic_spiral",
    "curvature_radius_m",
    "tension_acceleration_radial",
]
