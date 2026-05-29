"""
Monte-Carlo propagation of input systematics onto the gap-coverage ratio.

The headline column ``theoretical_gap_coverage_ratio`` is computed from three input
quantities that each carry a catalogued 1-σ uncertainty: distance ``D`` (Mpc),
inclination ``i`` (deg), and asymptotic flat rotation velocity ``V_flat`` (km/s).
Re-running the full engine ``N`` times per galaxy would be wasteful — the per-shell
kinematic rescaling is analytic, so we sample the inputs independently and propagate
them through a closed-form transformation:

* The observed centripetal acceleration scales like ``V_obs² / r``. SPARC's ``V_obs``
  is already deprojected (``V_obs = V_LOS / sin i``); resampling ``i_sample`` therefore
  rescales ``V_obs → V_obs · sin(i_true) / sin(i_sample)``.
* Distance enters linearly in the conversion ``r_kpc → r_m``; resampling ``D_sample``
  rescales ``r → r · (D_sample / D)``.
* The flat-rotation amplitude ``V_flat`` is a normalisation knob on the outer disk —
  it multiplies the kinematic gap coverage roughly as ``(V_flat_sample / V_flat)²``.

The combined scaling reads

.. math::

    R_{\\mathrm{ratio}}(D_s, i_s, V_s)
      \\;=\\;
      R_{\\mathrm{ratio},0}
      \\cdot \\left(\\frac{D_s}{D_0}\\right)^{-1}
      \\cdot \\left(\\frac{\\sin i_s}{\\sin i_0}\\right)^{2}
      \\cdot \\left(\\frac{V_0}{V_s}\\right)^{2},

which is a multiplicative model in the inputs (the engine's denominator and Δx scale
linearly with ``D``, the gap scales like ``V_obs² / r``). The caller may override the
default closure to plug in a different scaling — e.g. zero distance dependence, or a
non-quadratic ``V_flat`` exponent — through ``base_ratio_callable``.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from src.config.thresholds import UNCERTAINTY


def _safe_sin(angle_deg: float | np.ndarray, *, floor: float = 1.0e-3) -> np.ndarray:
    """``sin(angle_deg)`` clamped to ``floor`` to avoid 1/sin(i) blow-ups in the band."""
    sin_val = np.sin(np.deg2rad(np.asarray(angle_deg, dtype=float)))
    return np.where(np.abs(sin_val) < floor, np.copysign(floor, sin_val), sin_val)


def default_ratio_scaling(
    *,
    base_ratio: float,
    distance_mpc: float,
    distance_sample: np.ndarray,
    inclination_deg: float,
    inclination_sample: np.ndarray,
    vflat_kms: float,
    vflat_sample: np.ndarray,
) -> np.ndarray:
    """
    Closed-form multiplicative rescaling for ``theoretical_gap_coverage_ratio``.

    Returns ``base_ratio * (D_s/D)^-1 * (sin i_s / sin i)^2 * (V/V_s)^2``.

    Each input array is expected to have shape ``(n_samples,)``. The reference values
    ``distance_mpc``, ``inclination_deg`` and ``vflat_kms`` are point estimates.
    """
    if not np.isfinite(base_ratio):
        return np.full(distance_sample.shape, np.nan, dtype=float)

    d_ref = float(distance_mpc) if np.isfinite(distance_mpc) and distance_mpc > 0 else np.nan
    v_ref = float(vflat_kms) if np.isfinite(vflat_kms) and vflat_kms > 0 else np.nan
    sin_i_ref = float(_safe_sin(inclination_deg))
    sin_i_s = _safe_sin(inclination_sample)

    out = np.full(distance_sample.shape, np.nan, dtype=float)
    valid_idx = (
        np.isfinite(distance_sample) & (distance_sample > 0)
        & np.isfinite(inclination_sample)
        & np.isfinite(vflat_sample) & (vflat_sample > 0)
        & np.isfinite(d_ref) & np.isfinite(v_ref) & (sin_i_ref != 0.0)
    )
    if not np.any(valid_idx):
        return out

    d_factor = d_ref / distance_sample[valid_idx]
    i_factor = (sin_i_s[valid_idx] / sin_i_ref) ** 2
    v_factor = (v_ref / vflat_sample[valid_idx]) ** 2
    out[valid_idx] = float(base_ratio) * d_factor * i_factor * v_factor
    return out


def propagate_gap_coverage_band(
    *,
    base_ratio: float,
    distance_mpc: float | None,
    e_distance_mpc: float | None,
    inclination_deg: float | None,
    e_inclination_deg: float | None,
    vflat_kms: float | None,
    e_vflat_kms: float | None,
    base_ratio_callable: Callable[..., np.ndarray] = default_ratio_scaling,
    n_samples: int = UNCERTAINTY.MC_DEFAULT_SAMPLES,
    seed: int = UNCERTAINTY.MC_DEFAULT_SEED,
) -> tuple[float, float, float]:
    """
    Monte-Carlo p16 / median / p84 band on ``theoretical_gap_coverage_ratio``.

    Samples ``(D, i, V_flat)`` independently from Gaussian priors centred on the point
    estimates with widths given by the catalogue 1-σ uncertainties, propagates each draw
    through ``base_ratio_callable``, and returns the empirical (median, p16, p84) tuple.

    Degenerate inputs collapse the band to the point estimate:

    * Any missing / non-finite uncertainty input → returns ``(base_ratio, base_ratio,
      base_ratio)`` (no propagation possible — the band cannot be widened).
    * Zero uncertainty on all three inputs → returns the same triple (band collapses to
      the point estimate by construction).

    Parameters
    ----------
    base_ratio:
        Engine-computed ``theoretical_gap_coverage_ratio`` at the central inputs.
    distance_mpc, e_distance_mpc:
        Galaxy luminosity distance and its 1-σ catalogue uncertainty (Mpc).
    inclination_deg, e_inclination_deg:
        Disk inclination and its 1-σ catalogue uncertainty (deg).
    vflat_kms, e_vflat_kms:
        Asymptotic flat rotation velocity and its catalogue uncertainty (km/s).
    base_ratio_callable:
        Function applied to each MC draw to recompute the ratio. Defaults to
        :func:`default_ratio_scaling`. Custom callables receive the keyword arguments
        ``base_ratio``, ``distance_mpc``, ``distance_sample``, ``inclination_deg``,
        ``inclination_sample``, ``vflat_kms``, ``vflat_sample`` and must return an
        ``np.ndarray`` of shape ``(n_samples,)``.
    n_samples:
        Number of Monte-Carlo draws (defaults to
        :attr:`src.config.thresholds.UncertaintyConstants.MC_DEFAULT_SAMPLES`).
    seed:
        RNG seed for ``np.random.default_rng``. Determinism is mandatory.

    Returns
    -------
    tuple[float, float, float]
        ``(median, p16, p84)`` of the propagated ratio. All three entries equal
        ``base_ratio`` when propagation cannot be performed.
    """
    if not np.isfinite(base_ratio):
        return float("nan"), float("nan"), float("nan")

    def _ok(val: float | None) -> bool:
        return val is not None and np.isfinite(float(val))

    if not (
        _ok(distance_mpc) and _ok(inclination_deg) and _ok(vflat_kms)
        and _ok(e_distance_mpc) and _ok(e_inclination_deg) and _ok(e_vflat_kms)
    ):
        return float(base_ratio), float(base_ratio), float(base_ratio)

    n = int(max(n_samples, 1))
    rng = np.random.default_rng(int(seed))

    d_samples = rng.normal(loc=float(distance_mpc), scale=max(float(e_distance_mpc), 0.0), size=n)
    i_samples = rng.normal(loc=float(inclination_deg), scale=max(float(e_inclination_deg), 0.0), size=n)
    v_samples = rng.normal(loc=float(vflat_kms), scale=max(float(e_vflat_kms), 0.0), size=n)
    # Constrain physical priors: distances and velocities must remain positive; inclinations
    # within (0, 90]. We reflect at the boundaries with a small floor to keep sample sizes.
    d_samples = np.clip(d_samples, a_min=1.0e-3 * float(distance_mpc), a_max=None)
    v_samples = np.clip(v_samples, a_min=1.0e-3 * float(vflat_kms), a_max=None)
    i_samples = np.clip(i_samples, a_min=1.0, a_max=89.9)

    ratios = base_ratio_callable(
        base_ratio=float(base_ratio),
        distance_mpc=float(distance_mpc),
        distance_sample=d_samples,
        inclination_deg=float(inclination_deg),
        inclination_sample=i_samples,
        vflat_kms=float(vflat_kms),
        vflat_sample=v_samples,
    )
    finite = ratios[np.isfinite(ratios)]
    if finite.size == 0:
        return float(base_ratio), float(base_ratio), float(base_ratio)
    median = float(np.median(finite))
    p16, p84 = (float(q) for q in np.quantile(finite, [0.16, 0.84]))
    return median, p16, p84


__all__ = [
    "default_ratio_scaling",
    "propagate_gap_coverage_band",
]
