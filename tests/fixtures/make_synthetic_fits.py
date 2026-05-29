"""
Create minimal HEALPix-compatible FITS tables for PlanckProvider tests.

The production Planck maps use FITS extension 1 with columns ``Q_STOKES`` and ``U_STOKES``
indexed by HEALPix pixel (one row per pixel). This helper writes the same layout at low NSIDE.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits


def write_polarization_fits(path: Path, *, nside: int = 16, ordering: str = "RING") -> Path:
    """
    Write a synthetic polarization map: Q and U constant over the sphere (uniform polarization).

    Parameters
    ----------
    path:
        Output ``.fits`` path (parent directories are created).
    nside:
        HEALPix resolution parameter (must be a power of two).
    ordering:
        ``RING`` or ``NESTED`` — must match what ``PlanckProvider`` expects for ``ang2pix``.
    """
    npix = 12 * nside * nside
    q = np.full(npix, 1e-5, dtype=np.float64)
    u = np.full(npix, 2e-5, dtype=np.float64)
    # One table row per HEALPix pixel (matches Planck ``Q_STOKES`` / ``U_STOKES`` columns).
    col_q = fits.Column(name="Q_STOKES", format="D", array=q)
    col_u = fits.Column(name="U_STOKES", format="D", array=u)
    hdu1 = fits.BinTableHDU.from_columns([col_q, col_u])
    hdu1.header["NSIDE"] = nside
    hdu1.header["ORDERING"] = ordering
    hdu0 = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu0, hdu1])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    hdul.writeto(path, overwrite=True)
    return path
