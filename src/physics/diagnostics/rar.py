"""
Radial Acceleration Relation (McGaugh, Lelli & Schombert 2016).

Empirical fit linking the **observed** centripetal acceleration ``g_obs = V_obs² / r``
to the **purely baryonic** acceleration ``g_bar = V_bar² / r``:

.. math::
    g_{\\mathrm{obs}}(r) \\;=\\; \\frac{g_{\\mathrm{bar}}(r)}{1 - \\exp\\!\\bigl(-\\sqrt{g_{\\mathrm{bar}}(r)/a_0}\\bigr)},

with ``a_0 ≈ 1.2 × 10⁻¹⁰ m/s²`` (taken from ``ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2``).

Why we compute it
-----------------
This is the standard empirical relation that any rotation-curve-residual model is expected to
**reproduce** if it is to compete with the dark-matter explanation on the SPARC dataset. The
pipeline does NOT use this for fitting the magnetic-pressure model — it adds it as a per-radius
column so reports can show, side by side, the **gap that DM is invoked to close** (already in
``acceleration_gap``) and the **RAR-predicted** centripetal acceleration that any successful
alternative model has to recover.

The function is numerically safe at small ``g_bar`` (where the ``1 - exp(-x)`` denominator goes
to zero like ``x`` and the prediction asymptotes to ``√(g_bar · a₀)``).
"""

from __future__ import annotations

import numpy as np

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS


def rar_predicted_observed_acceleration(
    baryonic_acceleration_si: np.ndarray | float,
    a0_m_s2: float = ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2,
) -> np.ndarray:
    """
    McGaugh+2016 RAR interpolating function ``g_obs(g_bar)``.

    Parameters
    ----------
    baryonic_acceleration_si:
        Per-radius baryonic acceleration ``g_bar = V_bar² / r`` in SI (m/s²). Scalar or array.
    a0_m_s2:
        Acceleration scale ``a₀`` in SI; defaults to the value in
        ``ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2``.

    Returns
    -------
    np.ndarray
        Predicted observed centripetal acceleration ``g_obs(r)`` in SI (m/s²) with the same
        shape as the input. Non-finite inputs propagate as ``NaN``.
    """
    g_bar = np.asarray(baryonic_acceleration_si, dtype=float)
    out = np.full_like(g_bar, np.nan, dtype=float)
    valid = np.isfinite(g_bar) & (g_bar > 0) & (a0_m_s2 > 0)
    if not np.any(valid):
        return out
    x = np.sqrt(g_bar[valid] / a0_m_s2)
    # ``1 - exp(-x)`` is well-behaved for x>0; np.expm1 keeps small-x precision.
    denom = -np.expm1(-x)  # = 1 - exp(-x)
    # When denom would underflow to zero (g_bar -> 0), the limit is sqrt(g_bar * a0).
    safe = denom > 0
    out_valid = np.where(
        safe,
        g_bar[valid] / np.where(safe, denom, 1.0),
        np.sqrt(g_bar[valid] * a0_m_s2),
    )
    out[valid] = out_valid
    return out


def rar_residual_log10(
    observed_acceleration_si: np.ndarray | float,
    baryonic_acceleration_si: np.ndarray | float,
    a0_m_s2: float = ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2,
) -> np.ndarray:
    """
    Per-radius residual of the rotation curve against the RAR, in dex.

    Returns ``log10(g_obs_measured / g_obs_RAR_predicted)``: zero means perfect agreement with
    the RAR, positive means the galaxy needs *more* support than RAR predicts, negative means
    *less*. Non-finite or non-positive inputs return ``NaN``.
    """
    g_obs = np.asarray(observed_acceleration_si, dtype=float)
    g_pred = rar_predicted_observed_acceleration(baryonic_acceleration_si, a0_m_s2=a0_m_s2)
    out = np.full_like(g_obs, np.nan, dtype=float)
    ok = np.isfinite(g_obs) & np.isfinite(g_pred) & (g_obs > 0) & (g_pred > 0)
    out[ok] = np.log10(g_obs[ok] / g_pred[ok])
    return out


def rar_median_abs_residual_dex(
    observed_acceleration_si: np.ndarray | float,
    baryonic_acceleration_si: np.ndarray | float,
    a0_m_s2: float = ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2,
) -> float | None:
    """
    Single-number RAR consistency statistic: ``nanmedian(|rar_residual_log10|)`` in dex.

    The median absolute residual is the cleanest scalar summary of how far a galaxy sits
    from the McGaugh RAR; it is robust to a couple of outlier shells (unlike the mean
    absolute residual) and is the standard quantity reported in the literature
    (Lelli+2017, eq. 4). Returns ``None`` when no shells produce a finite residual.

    Parameters
    ----------
    observed_acceleration_si, baryonic_acceleration_si:
        Per-shell ``g_obs`` and ``g_bar`` in SI (m/s²).
    a0_m_s2:
        Acceleration scale ``a₀`` passed through to :func:`rar_residual_log10`.

    Returns
    -------
    float | None
        Median absolute residual in dex, or ``None`` if no finite shell residual exists.
    """
    residuals = rar_residual_log10(
        observed_acceleration_si, baryonic_acceleration_si, a0_m_s2=a0_m_s2
    )
    finite = np.asarray(residuals, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None
    return float(np.median(np.abs(finite)))
