"""
Negative-control e2e test: a synthetic galaxy with a pure NFW dark-matter rotation curve must
*not* be matched by the magnetic-pressure model when the input Planck map is uniform (no real
gradient). This is the falsifiability check for the MHD-pressure hypothesis.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

from src.pipeline.cli import parse_args
from src.pipeline.runner import run_pipeline_with_root
from src.reporting.galaxy_from_disk import summary_json_path
from tests.fixtures.synthetic_nfw import write_nfw_rotmod


@pytest.fixture
def nfw_sources(
    tmp_path: Path,
    sample_morphology_json: Path,
    synthetic_planck_fits: Path,
):
    """Mirror the per-provenance ``sources/`` layout under ``tmp_path`` for the NFW control."""
    root = tmp_path / "sources"
    sparc_root = root / "sparc"
    rot_dir = sparc_root / "rotation_curves"
    rot_dir.mkdir(parents=True)
    write_nfw_rotmod(rot_dir / "SYN_NFW_rotmod.dat", galaxy_name="SYN_NFW", distance_mpc=10.0)
    # Reuse the small morphology JSON; the galaxy is NOT in it, so it will fall through to
    # ``Unclassified`` morphology (intentional: we don't want morphology to bias this test).
    shutil.copy(sample_morphology_json, sparc_root / "SPARC_Lelli2016c.json")
    dust = root / "planck"
    dust.mkdir(parents=True)
    shutil.copy(
        synthetic_planck_fits,
        dust / "COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits",
    )
    (root / "ned").mkdir(parents=True, exist_ok=True)
    return tmp_path


def test_synthetic_nfw_is_not_matched_by_uniform_polarization(
    monkeypatch: pytest.MonkeyPatch, nfw_sources: Path, mock_ned_success
):
    """Uniform Q/U map → ΔP_B ≈ 0 → a_B ≈ 0 → cannot reproduce the NFW gap.

    This is the *falsifiability* check: when the polarization map carries no spatial structure,
    the engine must NOT report the synthetic NFW halo as a successful magnetic-pressure match.
    Specifically, the gap-coverage ratio must remain well below the success band, OR the
    scenario must fail other reliability gates.
    """
    pytest.importorskip("healpy")
    monkeypatch.setattr(sys, "argv", ["pytest", "--galaxies", "SYN_NFW"])
    args = parse_args()
    run_pipeline_with_root(nfw_sources, args)

    summary_path = summary_json_path(nfw_sources, "SYN_NFW")
    assert summary_path.exists(), "Pipeline must produce a per-galaxy summary even for the NFW control."
    summary = json.loads(summary_path.read_text(encoding="utf-8"))

    # Best scenario: ratio must be way off the success band (the engine cannot match a real
    # DM-equivalent residual using a uniform polarization map).
    best = summary["best_theoretical_scenario"]
    ratio = best.get("theoretical_gap_coverage_ratio")
    # Ratio is either NaN (skipped) or ≪ 1 (massively below the matching band).
    if ratio is not None:
        try:
            ratio_f = float(ratio)
        except (TypeError, ValueError):
            ratio_f = float("nan")
        # Require: the synthetic NFW must NOT pass the strong-match band.
        assert not (0.5 <= ratio_f <= 1.5), (
            f"Synthetic NFW reached suspicious match ratio {ratio_f}; the uniform polarization "
            "map should not reproduce a true dark-matter rotation curve."
        )
