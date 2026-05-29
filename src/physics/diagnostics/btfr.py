"""
Baryonic Tully–Fisher Relation (McGaugh 2012).

The BTFR connects the asymptotic flat rotation velocity ``V_flat`` to the *baryonic*
mass ``M_bar = M_gas + M_star`` of a galaxy through

.. math::

    \\log_{10}\\!\\bigl(M_{\\mathrm{bar}}/M_\\odot\\bigr)
      \\;=\\; A\\,\\log_{10}\\!\\bigl(V_{\\mathrm{flat}}/(\\mathrm{km/s})\\bigr) + B,

with ``A = 4.0`` and ``B = 1.99`` in the McGaugh+2012 calibration (SPARC + extended HI
sample). This module exposes the ridge as a pure function plus a per-galaxy residual
``btfr_residual_dex`` and a warning flag for galaxies sitting more than
:attr:`src.config.thresholds.AstrophysicalDefaults.BTFR_OUTLIER_THRESHOLD_DEX` from the
ridge.

The baryonic-mass reconstruction follows the SPARC convention:

* ``M_gas = 1.33 · M_HI`` (helium correction; Lelli+2016).
* ``M_star = Υ_{3.6} · L_{3.6}`` with ``Υ_{3.6} = 0.5 M_⊙/L_⊙`` (Schombert+2014,
  adopted in :attr:`src.config.thresholds.AstrophysicalDefaults.STELLAR_MASS_TO_LIGHT_3P6UM`).

When either ``M_HI`` or ``L_{3.6}`` is missing the corresponding term contributes zero and
the residual carries the partial baryonic mass — flagged but not silently set to NaN so
the audit can still see the dwarf-galaxy regime where HI dominates and the stellar leg is
absent.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS


@dataclass(frozen=True, slots=True)
class BTFRRidgeFit:
    """Linear ridge ``log10(M_bar) = slope · log10(V_flat) + intercept``."""

    slope: float
    intercept: float

    def log_mass_at(self, vflat_kms: float | np.ndarray) -> np.ndarray:
        v = np.asarray(vflat_kms, dtype=float)
        out = np.full_like(v, np.nan, dtype=float)
        valid = np.isfinite(v) & (v > 0.0)
        out[valid] = self.slope * np.log10(v[valid]) + self.intercept
        return out


# Canonical McGaugh+2012 ridge baked into the configuration constants.
MCGAUGH_2012_RIDGE: BTFRRidgeFit = BTFRRidgeFit(
    slope=float(ASTROPHYSICAL_DEFAULTS.BTFR_RIDGE_SLOPE),
    intercept=float(ASTROPHYSICAL_DEFAULTS.BTFR_RIDGE_INTERCEPT),
)


def baryonic_mass_solar(
    *,
    hi_mass_1e9_solar: float | None,
    luminosity_3p6um_1e9_solar: float | None,
    helium_correction: float = ASTROPHYSICAL_DEFAULTS.HELIUM_CORRECTION_FACTOR_HI_TO_GAS,
    mass_to_light_3p6um: float = ASTROPHYSICAL_DEFAULTS.STELLAR_MASS_TO_LIGHT_3P6UM,
) -> float:
    """
    SPARC baryonic mass ``M_bar = 1.33 · M_HI + Υ_{3.6} · L_{3.6}`` in solar masses.

    Inputs are in units of ``10⁹ M_⊙`` (HI mass) and ``10⁹ L_⊙`` (3.6 µm luminosity), to
    match the SPARC catalogue columns ``MHI_1e9solMass`` and ``L36_1e9solLum``. Either
    leg may be missing — the function returns ``NaN`` when *both* legs are missing or
    non-finite, and otherwise returns the available partial sum.
    """
    legs: list[float] = []
    if hi_mass_1e9_solar is not None and np.isfinite(float(hi_mass_1e9_solar)) and float(hi_mass_1e9_solar) > 0.0:
        legs.append(float(helium_correction) * float(hi_mass_1e9_solar) * 1.0e9)
    if (
        luminosity_3p6um_1e9_solar is not None
        and np.isfinite(float(luminosity_3p6um_1e9_solar))
        and float(luminosity_3p6um_1e9_solar) > 0.0
    ):
        legs.append(float(mass_to_light_3p6um) * float(luminosity_3p6um_1e9_solar) * 1.0e9)
    if not legs:
        return float("nan")
    return float(sum(legs))


def btfr_residual_dex(
    *,
    vflat_kms: float | None,
    hi_mass_1e9_solar: float | None,
    luminosity_3p6um_1e9_solar: float | None,
    ridge: BTFRRidgeFit = MCGAUGH_2012_RIDGE,
) -> float:
    """
    Signed distance (dex) of the galaxy from the BTFR ridge.

    Positive residual ⇒ ``M_bar`` heavier than the ridge prediction at the given
    ``V_flat`` (galaxy lies *above* the BTFR locus); negative ⇒ lighter than predicted.
    Returns ``NaN`` when either ``V_flat`` or the baryonic-mass leg cannot be evaluated.
    """
    m_bar = baryonic_mass_solar(
        hi_mass_1e9_solar=hi_mass_1e9_solar,
        luminosity_3p6um_1e9_solar=luminosity_3p6um_1e9_solar,
    )
    if not np.isfinite(m_bar) or m_bar <= 0.0:
        return float("nan")
    if vflat_kms is None or not np.isfinite(float(vflat_kms)) or float(vflat_kms) <= 0.0:
        return float("nan")
    log_obs = float(np.log10(m_bar))
    log_pred = float(ridge.slope * np.log10(float(vflat_kms)) + ridge.intercept)
    return log_obs - log_pred


def is_btfr_outlier(
    residual_dex: float | None,
    threshold_dex: float = ASTROPHYSICAL_DEFAULTS.BTFR_OUTLIER_THRESHOLD_DEX,
) -> bool:
    """
    Return ``True`` when ``|residual_dex|`` exceeds the configured BTFR outlier threshold.

    Missing / non-finite residuals are *not* flagged (the warning is reserved for cases
    where the BTFR can be evaluated and disagrees with the ridge).
    """
    if residual_dex is None or not np.isfinite(float(residual_dex)):
        return False
    return abs(float(residual_dex)) > float(threshold_dex)


__all__ = [
    "BTFRRidgeFit",
    "MCGAUGH_2012_RIDGE",
    "baryonic_mass_solar",
    "btfr_residual_dex",
    "is_btfr_outlier",
]
