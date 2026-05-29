"""Planck provider behaviour when HEALPix bindings are unavailable."""

from pathlib import Path

import pytest

from src.providers.planck import PlanckProvider


def test_initialize_raises_when_healpy_missing_with_fits(monkeypatch, synthetic_planck_fits: Path):
    """A FITS map requires healpy for HEALPix indexing; missing healpy must fail fast."""
    monkeypatch.setattr(PlanckProvider, "_ensure_healpy", lambda self: False)
    prov = PlanckProvider([synthetic_planck_fits], verbose=False)
    with pytest.raises(RuntimeError, match="healpy"):
        prov.initialize()


def test_extract_returns_missing_inputs_when_healpy_false_without_initialize(monkeypatch, synthetic_planck_fits: Path):
    """Direct extract without a successful initialize still short-circuits when healpy is absent."""
    monkeypatch.setattr(PlanckProvider, "_ensure_healpy", lambda self: False)
    prov = PlanckProvider([synthetic_planck_fits], verbose=False)
    prov.fits_path = synthetic_planck_fits
    prov.nside = 16
    prov.healpix_nested = False
    out = prov.extract_stokes_qu_mean(180.0, 45.0)
    assert out.get("status") == "missing inputs"
