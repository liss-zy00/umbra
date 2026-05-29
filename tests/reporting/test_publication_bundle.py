"""Smoke tests for post-hoc publication plotting from on-disk ``results/`` layout."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from src.pipeline.constants import (
    DISCOVERY_RELIABLE_FILENAME,
    GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
    RESULTS_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
    RESULTS_GRAPHS_DIRNAME,
    RESULTS_GRAPHS_POPULATION_DIRNAME,
    RESULTS_STEP_DISCOVERY_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_RELIABLE_FILENAME,
    VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME,
)
from src.reporting.galaxy_from_disk import build_galaxy_from_galaxy_logs
from src.reporting.publication_bundle import (
    compute_cohort_force_y_limits,
    format_summary_metadata_footer,
    read_cohort_galaxy_names,
    run_publication_bundle,
)


def _write_minimal_analysis(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "radius_kpc": 1.0,
            "observed_velocity_kms": 50.0,
            "baryonic_velocity_kms": 40.0,
            "acceleration_gap": 1e-11,
            "magnetic_pressure_acceleration": 8e-12,
            "required_density_for_match": 2e-27,
        },
        {
            "radius_kpc": 2.0,
            "observed_velocity_kms": 55.0,
            "baryonic_velocity_kms": 42.0,
            "acceleration_gap": 9e-12,
            "magnetic_pressure_acceleration": 7e-12,
            "required_density_for_match": 2.5e-27,
        },
    ]
    pd.DataFrame(rows).to_csv(path, index=False)


def _write_minimal_summary(path: Path) -> None:
    payload = {
        "galaxy_identity": {"galaxy_name": "SYNP", "distance_local_mpc": 10.0},
        "astrophysical_metadata": {"distance_mpc": 10.0},
        "best_theoretical_scenario": {
            "scenario_name": "scale_0.100_ion_0.010",
            "scenario_gradient_scale_kpc": 0.1,
            "scenario_ionization_fraction": 0.01,
            "theoretical_gap_coverage_ratio": 0.5,
            "pearson_correlation": 0.25,
            "resolution_correction_factor": 2.0,
        },
        "uncertainty_propagation": {
            "theoretical_gap_coverage_ratio_p16": 0.45,
            "theoretical_gap_coverage_ratio_p84": 0.55,
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_build_galaxy_from_galaxy_logs_roundtrip(tmp_path: Path) -> None:
    gdir = tmp_path / RESULTS_DIRNAME / RESULTS_GALAXY_ANALYSIS_DIRNAME / GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME / "SYNP"
    _write_minimal_analysis(gdir / "SYNP_analysis.csv")
    _write_minimal_summary(gdir / "SYNP_summary.json")
    gal, summ = build_galaxy_from_galaxy_logs(tmp_path, "SYNP")
    assert gal.name == "SYNP"
    assert gal.distance_mpc == 10.0
    assert len(gal.rotation_dataframe) == 2
    assert summ is not None
    footer = format_summary_metadata_footer(summ)
    assert footer is not None
    assert "scale_0.100" in footer


def test_run_publication_bundle_writes_pngs(tmp_path: Path) -> None:
    root = tmp_path
    disc_root = root / RESULTS_DIRNAME / RESULTS_STEP_DISCOVERY_DIRNAME
    disc_root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"galaxy_name": "SYNP", "required_density_for_match": 1e-25, "convergence_ratio_to_wim": 10.0, "galactic_latitude_deg": 30.0, "beam_size_kpc": 5.0, "resolution_correction_factor": 2.0}]).to_csv(
        disc_root / DISCOVERY_RELIABLE_FILENAME, index=False
    )
    val_root = root / RESULTS_DIRNAME / RESULTS_STEP_VALIDATION_DIRNAME
    val_root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"galaxy_name": "SYNP", "mean_acceleration_gap_edge": 1e-11, "hubble_type": 10.0}]).to_csv(val_root / VALIDATION_RELIABLE_FILENAME, index=False)
    pd.DataFrame(
        [
            {
                "galaxy_name": "OTHER",
                "theoretical_gap_coverage_ratio": 0.9,
                "pearson_correlation": 0.5,
                "rar_median_abs_residual_dex": 0.1,
                "btfr_residual_dex": -0.05,
                "validation_status": "STANDARD",
                "morphology_class": "Dwarf Irregulars (Im)",
            }
        ]
    ).to_csv(val_root / VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME, index=False)

    gdir = root / RESULTS_DIRNAME / RESULTS_GALAXY_ANALYSIS_DIRNAME / GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME / "SYNP"
    _write_minimal_analysis(gdir / "SYNP_analysis.csv")
    _write_minimal_summary(gdir / "SYNP_summary.json")

    names = read_cohort_galaxy_names(root, "discovery_reliable")
    assert names == ["SYNP"]
    yl = compute_cohort_force_y_limits(root, names)
    assert yl is not None

    paths = run_publication_bundle(root, "discovery_reliable", write_companion_md=False)
    assert any(p.name == "SYNP_acceleration_comparison.png" for p in paths)
    out = root / RESULTS_DIRNAME / RESULTS_GRAPHS_DIRNAME / f"{RESULTS_GRAPHS_POPULATION_DIRNAME}_discovery_reliable"
    assert (gdir / "SYNP_acceleration_comparison.png").is_file()
    assert (out / "discovery_rho_req_vs_wim_ratio.png").is_file()
