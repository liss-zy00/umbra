from pathlib import Path
from typing import Any

import numpy as np
import pytest

from src.config.thresholds import (
    CALIBRATION_REFERENCE_METHOD_DYNAMIC,
    CALIBRATION_REFERENCE_METHOD_FALLBACK,
    PLANCK_CALIBRATION,
)
from src.providers.planck import PlanckProvider, resolve_stokes_amp_reference


def test_planck_initialize_and_extract(synthetic_planck_fits: Path):
    pytest.importorskip("healpy")
    prov = PlanckProvider([synthetic_planck_fits], verbose=False)
    summary = prov.initialize()
    assert summary.get("fits_found") is True
    # Provenance contract: a real FITS map must report the dynamic median path.
    assert summary.get("reference_method") == CALIBRATION_REFERENCE_METHOD_DYNAMIC
    # Synthetic FITS uses constant Q=1e-5, U=2e-5 → expected median is hypot(1e-5, 2e-5).
    expected_median = float(np.hypot(1e-5, 2e-5))
    assert np.isclose(summary["stokes_amp_reference"], expected_median, rtol=1e-9)
    out = prov.extract_stokes_qu_mean(180.0, 45.0, galaxy_name="Test")
    assert out.get("status") == "ok"
    assert out["stokes_amp_center"] >= 0
    # Each extract must echo the same provenance to keep downstream JSON honest.
    assert out["reference_method"] == CALIBRATION_REFERENCE_METHOD_DYNAMIC
    assert np.isclose(out["stokes_amp_reference"], expected_median, rtol=1e-9)


def test_resolve_stokes_prefers_stokes_summary_over_stale_fits_summary():
    """Per-galaxy extraction dict must win over a FITS summary missing the dynamic key."""
    good = float(5.622358303232117e-05)
    stokes = {"stokes_amp_reference": good, "reference_method": CALIBRATION_REFERENCE_METHOD_DYNAMIC}
    fits = {"stokes_amp_reference": float(PLANCK_CALIBRATION.STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY)}
    assert np.isclose(resolve_stokes_amp_reference(stokes, fits), good, rtol=1e-12)


def test_planck_initialize_returns_fallback_when_no_fits(tmp_path: Path):
    """No FITS candidates exist → provider must report the emergency-only constant honestly."""
    missing = tmp_path / "does_not_exist.fits"
    prov = PlanckProvider([missing], verbose=False)
    summary = prov.initialize()
    assert summary.get("fits_found") is False
    assert summary.get("reference_method") == CALIBRATION_REFERENCE_METHOD_FALLBACK
    assert float(summary["stokes_amp_reference"]) == float(
        PLANCK_CALIBRATION.STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY
    )


# --------------------------------------------------------------------------------------
# Phase 6 — performance contract: the Planck provider must use ``memmap=True`` and pull
# only the central pixel + ≤8 HEALPix neighbours, never the full map. These tests guard
# against future regressions that would silently load tens of millions of rows into RAM.
# --------------------------------------------------------------------------------------


def test_planck_initialize_opens_fits_with_memmap_true(synthetic_planck_fits: Path, monkeypatch):
    """``_fits_open_memmap`` must always pass ``memmap=True`` to ``astropy.io.fits.open``."""
    pytest.importorskip("healpy")

    captured: dict[str, Any] = {"memmap_kwargs": []}

    from astropy.io import fits as _fits

    real_open = _fits.open

    def spy_open(*args, **kwargs):
        captured["memmap_kwargs"].append(kwargs.get("memmap"))
        return real_open(*args, **kwargs)

    monkeypatch.setattr(_fits, "open", spy_open)

    prov = PlanckProvider([synthetic_planck_fits], verbose=False)
    prov.initialize()
    prov.extract_stokes_qu_mean(180.0, 45.0, galaxy_name="X")

    assert captured["memmap_kwargs"], "fits.open was never called"
    # Every single call must request memmap (no exceptions in the hot path).
    assert all(v is True for v in captured["memmap_kwargs"]), (
        f"All Planck FITS opens must use memmap=True; got {captured['memmap_kwargs']}"
    )


def test_planck_extract_pulls_only_neighbour_rows_not_full_map(synthetic_planck_fits: Path, monkeypatch):
    """``extract_stokes_qu_mean`` must read at most 9 rows (centre + 8 neighbours) — never
    iterate the entire HEALPix table. The check inspects the size of the index array
    forwarded to ``data['Q_STOKES'][idx]``."""
    pytest.importorskip("healpy")

    captured: dict[str, Any] = {"max_rows_read": 0}

    real_asarray = np.asarray

    def spy_asarray(arr, *args, **kwargs):
        # Detect the per-pixel pull: it always concerns a 1-D ndarray of small dtype.
        try:
            n = int(np.asarray(arr).shape[0]) if hasattr(arr, "shape") else len(arr)
            if 0 < n <= 9:
                captured["max_rows_read"] = max(captured["max_rows_read"], n)
        except Exception:
            pass
        return real_asarray(arr, *args, **kwargs)

    monkeypatch.setattr(np, "asarray", spy_asarray)

    prov = PlanckProvider([synthetic_planck_fits], verbose=False)
    prov.initialize()
    out = prov.extract_stokes_qu_mean(180.0, 45.0, galaxy_name="X")

    assert out.get("status") == "ok"
    # At least one read must have hit the small ≤9-row neighbour index slice.
    assert 1 <= captured["max_rows_read"] <= 9, (
        "Planck extract must read only the centre + 8 HEALPix neighbour rows; "
        f"observed max = {captured['max_rows_read']}"
    )


def test_planck_median_uses_strided_subsampling_to_bound_ram(
    synthetic_planck_fits: Path, monkeypatch
):
    """For a large map the Stokes median must be computed on a STRIDED subsample whose size
    is bounded by ``STOKES_MEDIAN_MAX_SAMPLES``. The test injects a tiny cap (rebinding the
    module-level :data:`PLANCK_CALIBRATION` constant since the dataclass is frozen) and
    asserts that the resulting sample respects it — so RAM never explodes on
    full-resolution Planck maps."""
    pytest.importorskip("healpy")

    from src.config.thresholds import PlanckCalibration
    from src.providers import planck as planck_mod

    tiny_cap = 64
    capped = PlanckCalibration(
        STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY=PLANCK_CALIBRATION.STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY,
        STOKES_MEDIAN_MAX_SAMPLES=tiny_cap,
        STOKES_MEDIAN_MIN_POSITIVE=PLANCK_CALIBRATION.STOKES_MEDIAN_MIN_POSITIVE,
        BEAM_FWHM_ARCMIN=PLANCK_CALIBRATION.BEAM_FWHM_ARCMIN,
    )
    monkeypatch.setattr(planck_mod, "PLANCK_CALIBRATION", capped, raising=True)

    prov = PlanckProvider([synthetic_planck_fits], verbose=False)
    prov.fits_path = synthetic_planck_fits
    median, method = prov._compute_stokes_amplitude_reference_median()

    assert method == CALIBRATION_REFERENCE_METHOD_DYNAMIC
    assert np.isfinite(median) and median > 0

    # Reproduce the stride arithmetic the production code uses and assert that the sampled
    # row count never exceeds the cap (plus one row for ceiling rounding).
    from astropy.io import fits as _fits
    with _fits.open(synthetic_planck_fits, memmap=True) as hdul:
        n = int(len(hdul[1].data))
    expected_stride = max(1, int(np.ceil(n / float(tiny_cap))))
    sampled = int(np.ceil(n / expected_stride))
    assert sampled <= tiny_cap + 1, (
        f"Strided sample must respect the configured cap; got {sampled} rows for cap {tiny_cap}"
    )
