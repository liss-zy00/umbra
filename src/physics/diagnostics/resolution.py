"""
Instrumental-resolution (beam-smearing) diagnostic for the Planck-anchored gradient.

Background
----------
The magnetic-pressure gradient ``∇u_B`` is reconstructed from the central Planck pixel
and its angular neighbours. Each pixel integrates the polarised emission over an effective
beam of full-width-half-maximum (FWHM) ``θ_beam ≈ 5'`` (Planck GNILC R3.00). When the
galaxy sits at distance ``D``, the physical footprint of the beam on the disk plane is

.. math::

    D_{\\rm beam}(D) \\;=\\; D \\cdot \\theta_{\\rm beam}.

For ``D = 4`` Mpc and ``θ_beam = 5'`` this gives ``D_beam ≈ 5.8`` kpc — orders of
magnitude larger than the disk-scale-height-driven gradient the model assumes (a few
hundred parsecs). The pixel therefore *averages* magnetic-field structures with
characteristic scale ``L_int ≪ D_beam`` and the reconstructed gradient is **diluted**
relative to the underlying field. A point-mass MHD argument gives

.. math::

    \\frac{\\nabla u_B^{\\rm true}}{\\nabla u_B^{\\rm beam}}
        \\;\\approx\\;
        \\frac{D_{\\rm beam}}{L_{\\rm grad}},

where ``L_grad`` is the disk-plane gradient scale (the engine knob
``scenario_gradient_scale_kpc``, ~50–500 pc). The same ratio enters the inverse of the
*required-density* estimate; conversely

.. math::

    k(D)
        \\;\\equiv\\;
        \\frac{D_{\\rm beam}}{L_{\\rm grad}}
        \\;=\\;
        \\rho_{\\rm req}^{\\rm true} \\bigl/ \\rho_{\\rm req}^{\\rm reconstructed},

so ``required_density_for_match`` reported by the engine is overestimated by the same
factor whenever beam smearing dominates the reconstruction. A clean falsification test is
therefore *expected*: if ``ρ_req`` correlates positively with ``D`` (or equivalently with
``D_beam / L_grad``), the residual is instrumental — the underlying field is just being
averaged out by a too-wide beam — and not a failure of the E-G hypothesis.

This module exposes pure functions consumed by the engine, the persistence layer and
optional offline diagnostics (see ``src.support.beam_smearing_plot``). They are designed to be cheap (vectorised numpy ops on
per-galaxy scalars) and to never raise on missing inputs — non-finite arguments propagate
``NaN`` so the downstream report can flag the row without crashing.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np

from src.config.thresholds import PLANCK_CALIBRATION


# 1 arcminute expressed in radians (``π / (180 · 60)``); used to convert the catalogued
# Planck beam FWHM (arcminutes) into a small angle without round-trips through astropy.
_ARCMIN_TO_RAD: float = math.pi / (180.0 * 60.0)


def beam_fwhm_radians(beam_fwhm_arcmin: float = PLANCK_CALIBRATION.BEAM_FWHM_ARCMIN) -> float:
    """Return the Planck effective beam FWHM in radians."""
    return float(beam_fwhm_arcmin) * _ARCMIN_TO_RAD


def beam_physical_size_kpc(
    *,
    distance_mpc: float | None,
    beam_fwhm_arcmin: float = PLANCK_CALIBRATION.BEAM_FWHM_ARCMIN,
) -> float:
    """
    Physical projection of the Planck beam on the disk plane (kpc).

    ``D_beam = D · θ_beam`` with ``D`` in Mpc and ``θ_beam`` in radians. The factor 1000
    converts Mpc·rad into kpc. Returns ``NaN`` for missing / non-finite / non-positive
    distance — downstream consumers treat ``NaN`` as "diagnostic not evaluable".
    """
    if distance_mpc is None:
        return float("nan")
    try:
        d = float(distance_mpc)
    except (TypeError, ValueError):
        return float("nan")
    if not math.isfinite(d) or d <= 0.0:
        return float("nan")
    theta = beam_fwhm_radians(beam_fwhm_arcmin=beam_fwhm_arcmin)
    return d * theta * 1000.0


def resolution_correction_factor(
    *,
    distance_mpc: float | None,
    gradient_scale_kpc: float | None,
    beam_fwhm_arcmin: float = PLANCK_CALIBRATION.BEAM_FWHM_ARCMIN,
) -> float:
    """
    Beam-smearing dilution factor ``k = D_beam / L_grad`` (dimensionless).

    The factor estimates **how much the reconstructed magnetic gradient is diluted by
    the finite Planck beam**:

    * ``k ≈ 1`` — the beam matches the assumed gradient scale; the reconstruction is
      not contaminated by smearing.
    * ``k > 1`` — the beam is *wider* than the gradient scale; the reconstructed
      ``∇u_B`` is averaged out by a factor ~``k`` and the engine *over-estimates*
      ``required_density_for_match`` by the same factor.
    * ``k < 1`` — the beam is *narrower* than the gradient scale; the reconstruction
      is locally over-resolved (rare under SPARC distances).

    Parameters
    ----------
    distance_mpc:
        Galaxy luminosity distance (Mpc).
    gradient_scale_kpc:
        The disk-plane gradient scale assumed by the engine (kpc). For the
        Discovery / headline path this is the physical pixel projection
        ``D · θ_neighbour``; for the scenario grid it is the
        ``scenario_gradient_scale_kpc`` knob.
    beam_fwhm_arcmin:
        Planck effective beam FWHM (arcmin); defaults to
        :attr:`src.config.thresholds.PlanckCalibration.BEAM_FWHM_ARCMIN`.

    Returns
    -------
    float
        Dilution factor (dimensionless). ``NaN`` when either input is missing /
        non-finite / non-positive.
    """
    if gradient_scale_kpc is None:
        return float("nan")
    try:
        l_grad = float(gradient_scale_kpc)
    except (TypeError, ValueError):
        return float("nan")
    if not math.isfinite(l_grad) or l_grad <= 0.0:
        return float("nan")
    d_beam_kpc = beam_physical_size_kpc(
        distance_mpc=distance_mpc, beam_fwhm_arcmin=beam_fwhm_arcmin
    )
    if not math.isfinite(d_beam_kpc) or d_beam_kpc <= 0.0:
        return float("nan")
    return d_beam_kpc / l_grad


def vectorised_resolution_correction(
    *,
    distance_mpc: Iterable[float | None],
    gradient_scale_kpc: Iterable[float | None],
    beam_fwhm_arcmin: float = PLANCK_CALIBRATION.BEAM_FWHM_ARCMIN,
) -> np.ndarray:
    """
    Vectorised counterpart of :func:`resolution_correction_factor` for population tables.

    Both arrays are coerced to float; missing / non-positive entries map to ``NaN``.
    """
    d_arr = np.asarray(list(distance_mpc), dtype=float)
    l_arr = np.asarray(list(gradient_scale_kpc), dtype=float)
    valid = np.isfinite(d_arr) & np.isfinite(l_arr) & (d_arr > 0.0) & (l_arr > 0.0)
    theta = beam_fwhm_radians(beam_fwhm_arcmin=beam_fwhm_arcmin)
    out = np.full_like(d_arr, np.nan, dtype=float)
    out[valid] = (d_arr[valid] * theta * 1000.0) / l_arr[valid]
    return out


__all__ = [
    "beam_fwhm_radians",
    "beam_physical_size_kpc",
    "resolution_correction_factor",
    "vectorised_resolution_correction",
]
