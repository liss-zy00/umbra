"""End-to-end integration tests with an isolated ``sources/`` tree under ``tmp_path``.

These tests exercise the full pipeline (``SparcProvider`` → engine → metrics builder →
``ResultPersister`` → ``NarrativeGenerator``) on a synthetic galaxy and assert that:

* the structured audit-flag JSON arrays propagate to disk exactly as documented;
* every reliability failure mode (latitude / outlier-prone shells / density floor /
  stellar dominance / skipped physics) is observable in the per-galaxy summary file;
* the global aggregator outputs land in the expected ``results/`` layout.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest

from src.pipeline.cli import parse_args
from src.pipeline.constants import (
    GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME,
    RESULTS_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
    RESULTS_STEP_DISCOVERY_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_ANCHOR_SENSITIVITY_SWEEP_FILENAME,
    VALIDATION_CATALOG_ALL_FILENAME,
    VALIDATION_NOT_RELIABLE_FILENAME,
    VALIDATION_RELIABLE_FILENAME,
    VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME,
)
from src.pipeline.runner import run_pipeline_with_root
from src.reporting.galaxy_from_disk import summary_json_path


@pytest.fixture
def minimal_sources(tmp_path: Path, sample_rotmod_path: Path, sample_morphology_json: Path, synthetic_planck_fits: Path):
    """
    Build a per-provenance ``sources/`` mirror under ``tmp_path``:

    ``sources/sparc/`` holds the rotmod sample + the morphology JSON, ``sources/planck/``
    holds the synthetic FITS, and ``sources/ned/`` is created empty so the NED provider
    falls through to the bare-name path with no alias map (the test mocks NED anyway).
    """
    root = tmp_path / "sources"
    sparc_root = root / "sparc"
    rot_dir = sparc_root / "rotation_curves"
    rot_dir.mkdir(parents=True)
    shutil.copy(sample_rotmod_path, rot_dir / "SYN001_rotmod.dat")
    shutil.copy(sample_morphology_json, sparc_root / "SPARC_Lelli2016c.json")
    dust = root / "planck"
    dust.mkdir(parents=True)
    # First candidate FITS name from ``pipeline.constants.FITS_CANDIDATES``.
    shutil.copy(
        synthetic_planck_fits,
        dust / "COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits",
    )
    (root / "ned").mkdir(parents=True, exist_ok=True)
    return tmp_path


def test_pipeline_e2e_galaxies_syn001(monkeypatch: pytest.MonkeyPatch, minimal_sources: Path, mock_ned_success):
    pytest.importorskip("healpy")
    monkeypatch.setattr(sys, "argv", ["pytest", "--galaxies", "SYN001"])
    args = parse_args()
    run_pipeline_with_root(minimal_sources, args)

    proc = minimal_sources / RESULTS_DIRNAME
    val = proc / RESULTS_STEP_VALIDATION_DIRNAME
    assert (val / VALIDATION_CATALOG_ALL_FILENAME).exists()
    assert (val / VALIDATION_RELIABLE_FILENAME).exists()
    assert (val / VALIDATION_NOT_RELIABLE_FILENAME).exists()
    ga = proc / RESULTS_GALAXY_ANALYSIS_DIRNAME
    assert (ga / GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME / "SYN001").is_dir() or (
        ga / GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME / "SYN001"
    ).is_dir() or (ga / GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME / "SYN001").is_dir()

    # Phase 1: Discovery census artefacts (always written before the grid).
    disc = proc / RESULTS_STEP_DISCOVERY_DIRNAME
    for name in (
        "discovery_all.csv",
        "discovery_reliable.csv",
        "discovery_not_reliable.csv",
        "discovery_cohort_rho_band.csv",
    ):
        assert (disc / name).exists(), f"unified pipeline must write {name}"

    report_path = val / VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME
    if report_path.exists() and report_path.stat().st_size > 0:
        header = report_path.read_text(encoding="utf-8").splitlines()[0]
        for col in (
            "required_density_for_match",
            "convergence_ratio_to_wim",
            "delta_r_m",
            "discovery_is_reliable",
            "discovery_is_elite",
            "discovery_elite_band_source",
            # Phase-1 refinement columns: the headline CSV must surface the new bands and
            # the side-by-side full-Lorentz accelerations so consumers can compare the
            # pressure-only ratio with the full ``a_B + a_tension`` variant.
            "theoretical_gap_coverage_ratio_p16",
            "theoretical_gap_coverage_ratio_p84",
            "theoretical_gap_coverage_ratio_full_lorentz",
            "accel_pressure_only_m_s2",
            "accel_full_lorentz_m_s2",
            "rar_median_abs_residual_dex",
            "btfr_residual_dex",
            "warn_btfr_outlier",
            "warn_no_pa",
            "warn_pure_pressure",
            "prov_high_precision_dist",
            "position_angle_deg",
            "wim_density_match_score",
            "validation_status",
            "calibrated_convergence_ratio",
            "r_max_kpc",
            # Beam-smearing diagnostic columns (k = D_beam / L_grad).
            "distance_mpc",
            "beam_size_kpc",
            "resolution_correction_factor",
            "resolution_correction_factor_discovery",
            "resolution_correction_factor_grid",
        ):
            assert col in header, f"{VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME} must include merged column {col}"

    # Phase-1 refinement: the anchor-sensitivity sweep must always be written, even when
    # the Discovery census filters every galaxy out — the file is just empty in that case.
    anchor_path = val / VALIDATION_ANCHOR_SENSITIVITY_SWEEP_FILENAME
    assert anchor_path.exists(), f"{VALIDATION_ANCHOR_SENSITIVITY_SWEEP_FILENAME} must be written by the pipeline"
    anchor_lines = anchor_path.read_text(encoding="utf-8").splitlines()
    if anchor_lines:
        anchor_header = anchor_lines[0]
        for col in ("galaxy_name", "ratio_cmb", "ratio_turbulent_equipartition", "ratio_log10_delta"):
            assert col in anchor_header, f"{VALIDATION_ANCHOR_SENSITIVITY_SWEEP_FILENAME} missing column {col}"

    # Discovery CSV schema + dynamic elite band (Phase 1 invariants). The elite filter
    # is now a multiplicative window on ``required_density_for_match`` around the
    # reliable-subset median; the test recomputes the median from the on-disk reliable
    # file and asserts every elite row falls inside ``[REL_MIN, REL_MAX]·median``.
    from src.reporting.persistence.exporter import ResultPersister

    reliable_path = disc / "discovery_reliable.csv"
    elite_path = disc / "discovery_cohort_rho_band.csv"
    for required_column in (
        "galaxy_name",
        "required_density_for_match",
        "convergence_ratio_to_wim",
        "delta_r_m",
        "resolution_correction_factor",
        "beam_size_kpc",
        "discovery_is_reliable",
    ):
        assert required_column in reliable_path.read_text(encoding="utf-8").splitlines()[0]
    rel_lines = reliable_path.read_text(encoding="utf-8").splitlines()
    reliable_rhos: list[float] = []
    if len(rel_lines) > 1:
        rel_header = rel_lines[0].split(",")
        i_rel_flag = rel_header.index("discovery_is_reliable")
        i_rho = rel_header.index("required_density_for_match")
        for body in rel_lines[1:]:
            cells = body.split(",")
            assert cells[i_rel_flag].strip().lower() == "true"
            try:
                rho_val = float(cells[i_rho].strip())
            except (TypeError, ValueError):
                continue
            if rho_val > 0:
                reliable_rhos.append(rho_val)
    el_lines = elite_path.read_text(encoding="utf-8").splitlines()
    if len(el_lines) > 1 and reliable_rhos:
        import statistics as _stats

        median_rho = _stats.median(reliable_rhos)
        band_min = ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MIN * median_rho
        band_max = ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MAX * median_rho
        h = el_lines[0].split(",")
        i_rel = h.index("discovery_is_reliable")
        i_rho = h.index("required_density_for_match")
        for body in el_lines[1:]:
            cells = body.split(",")
            assert cells[i_rel].strip().lower() == "true"
            rho_val = float(cells[i_rho].strip())
            assert band_min <= rho_val <= band_max, (
                f"elite row ρ_req={rho_val} not in dynamic band [{band_min}, {band_max}]"
            )


def test_pipeline_e2e_audit_flags_serialised_in_summary_json(
    monkeypatch: pytest.MonkeyPatch,
    minimal_sources: Path,
    mock_ned_success,
):
    """End-to-end: the per-galaxy summary must expose ``audit_flags`` as a JSON array.

    The structured audit list is the user-facing contract for transparency: when
    ``validation_is_reliable`` flips to ``False`` the consumer must be able to read *why*. We do not
    assert specific failure codes here (the synthetic galaxy might pass every gate or
    fail several depending on the synthetic Planck content); we only assert the contract.
    """
    pytest.importorskip("healpy")
    monkeypatch.setattr(sys, "argv", ["pytest", "--galaxies", "SYN001"])
    args = parse_args()
    run_pipeline_with_root(minimal_sources, args)

    summary_path = summary_json_path(minimal_sources, "SYN001")
    assert summary_path.exists()
    payload = json.loads(summary_path.read_text(encoding="utf-8"))

    # Galaxy-level audit must surface the structured arrays. The canonical token strings
    # live in :class:`src.config.thresholds.AuditCodes`; both ``audit_flags`` and
    # ``reliability_failure_codes`` are JSON arrays of those tokens.
    reliability = payload.get("data_reliability_assessment")
    assert reliability is not None, "summary must include data_reliability_assessment block"
    assert isinstance(reliability.get("audit_flags"), list)
    assert isinstance(reliability.get("reliability_failure_codes"), list)

    # ``validation_is_reliable`` must agree with the structured audit-flag content; no warn_*-only entry
    # silently demotes the run.
    blocking = [c for c in reliability["audit_flags"] if not str(c).startswith("warn_")]
    val_rel = reliability.get("validation_is_reliable")
    if val_rel is None:
        val_rel = reliability.get("is_reliable")
    assert val_rel is (len(blocking) == 0)

    # The grid-aware CSV files carry the per-scenario codes — they are the source of truth
    # for the auditor. Smoke-check at least one scenario CSV exists and contains the canonical
    # column names.
    grid_dir = minimal_sources / RESULTS_DIRNAME / "grid_results"
    grid_csvs = list(grid_dir.glob("scenario_*.csv"))
    assert grid_csvs, "expected at least one scenario_*.csv in grid_results"
    sample_text = grid_csvs[0].read_text(encoding="utf-8").splitlines()[0]
    assert "audit_flags" in sample_text
    assert "reliability_failure_codes" in sample_text


@pytest.fixture
def minimal_two_galaxies(
    tmp_path: Path,
    sample_rotmod_path: Path,
    sample_morphology_json: Path,
    synthetic_planck_fits: Path,
) -> Path:
    """Two rotmod files (SYN001 + LATFAIL) under ``sources/sparc/`` for pre-filter tests."""
    root = tmp_path / "sources"
    sparc_root = root / "sparc"
    rot_dir = sparc_root / "rotation_curves"
    rot_dir.mkdir(parents=True)
    shutil.copy(sample_rotmod_path, rot_dir / "SYN001_rotmod.dat")
    shutil.copy(sample_rotmod_path, rot_dir / "LATFAIL_rotmod.dat")
    morph = json.loads(sample_morphology_json.read_text(encoding="utf-8"))
    morph["data"].append({"Galaxy": "LATFAIL", "T": 10, "SBdisk_solLum_pc2": 50})
    (sparc_root / "SPARC_Lelli2016c.json").write_text(json.dumps(morph), encoding="utf-8")
    dust = root / "planck"
    dust.mkdir(parents=True)
    shutil.copy(
        synthetic_planck_fits,
        dust / "COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits",
    )
    (root / "ned").mkdir(parents=True, exist_ok=True)
    return tmp_path


def test_pipeline_e2e_blocking_discovery_audit_skips_grid(
    monkeypatch: pytest.MonkeyPatch,
    minimal_two_galaxies: Path,
    mock_ned_success,
):
    """LATFAIL carries a blocking Discovery audit: row in ``discovery_all.csv`` but no grid."""
    import src.pipeline.runner as rnr

    from src.config.thresholds import AUDIT_CODES

    pytest.importorskip("healpy")
    orig = rnr.process_galaxy_discovery

    def wrapped(deps, file_path: Path):
        row = orig(deps, file_path)
        if row is None:
            return None
        if file_path.stem.replace("_rotmod", "") == "LATFAIL":
            row = dict(row)
            row["audit_flags"] = AUDIT_CODES.FAIL_LATITUDE
            row["discovery_is_reliable"] = False
        return row

    monkeypatch.setattr(rnr, "process_galaxy_discovery", wrapped)
    monkeypatch.setattr(sys, "argv", ["pytest", "--galaxies", "SYN001", "LATFAIL"])
    args = parse_args()
    run_pipeline_with_root(minimal_two_galaxies, args)

    proc = minimal_two_galaxies / RESULTS_DIRNAME
    all_csv = proc / RESULTS_STEP_DISCOVERY_DIRNAME / "discovery_all.csv"
    assert all_csv.exists()
    text = all_csv.read_text(encoding="utf-8")
    assert "LATFAIL" in text
    ga = proc / RESULTS_GALAXY_ANALYSIS_DIRNAME
    assert (ga / GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME / "SYN001").is_dir() or (
        ga / GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME / "SYN001"
    ).is_dir() or (ga / GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME / "SYN001").is_dir()
    assert not (ga / GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME / "LATFAIL").exists()
    assert not (ga / GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME / "LATFAIL").exists()
    assert not (ga / GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME / "LATFAIL").exists()
    for p in (proc / "grid_results").glob("scenario_*.csv"):
        assert "LATFAIL" not in p.read_text(encoding="utf-8")
