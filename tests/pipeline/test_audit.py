"""Smoke + behavior tests for the consistency-audit module."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from src.pipeline.aggregation.audit import (
    build_consistency_audit_markdown,
    write_consistency_audit_report,
)
from src.pipeline.constants import (
    GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
)


def _write_summary(galaxy_analysis_root: Path, name: str, payload: dict) -> None:
    rel = bool((payload.get("data_reliability_assessment") or {}).get("validation_is_reliable", False))
    tier = GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME if rel else GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME
    galaxy_dir = galaxy_analysis_root / tier / name
    galaxy_dir.mkdir(parents=True, exist_ok=True)
    (galaxy_dir / f"{name}_summary.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=True), encoding="utf-8"
    )


def _make_payload(
    *,
    galaxy: str,
    morphology: str,
    validation_is_reliable: bool,
    ratio: float,
    galactic_lat: float,
    inclination: float,
    b_microgauss: float,
    plasma_beta: float,
    rar_residual: float,
    distance: float = 5.0,
) -> dict:
    return {
        "galaxy_identity": {
            "galaxy_name": galaxy,
            "morphology_class": morphology,
            "hubble_t_type": 9,
            "distance_local_mpc": distance,
        },
        "data_reliability_assessment": {
            "validation_is_reliable": validation_is_reliable,
            "galactic_latitude_degrees": galactic_lat,
        },
        "calibration_metadata": {"calibration_status": "OPTIMAL"},
        "best_theoretical_scenario": {
            "theoretical_gap_coverage_ratio": ratio,
            "pearson_correlation": 0.5,
        },
        "plasma_diagnostics": {
            "magnetic_field_microgauss": b_microgauss,
            "magnetic_field_within_plausibility_band": 0.1 <= b_microgauss <= 100.0,
            "plasma_beta_thermal_median": plasma_beta,
        },
        "rar_diagnostics": {"median_absolute_residual_dex": rar_residual},
        "density_audit": {
            "density_clamped_fraction": 0.05,
            "stellar_dominated_fraction": 0.1,
        },
        "kinematic_metadata": {"inclination_deg": inclination},
        "gradient_scale_audit": {"physical_neighbor_separation_kpc": 0.3},
    }


def test_audit_handles_empty_run(tmp_path: Path):
    processed = tmp_path / "results"
    (processed / RESULTS_GALAXY_ANALYSIS_DIRNAME).mkdir(parents=True)
    out = write_consistency_audit_report(processed_root=processed)
    assert out.exists()
    assert out.parent.name == "step_validation"
    assert out.name == "validation_consistency_audit.md"
    text = out.read_text(encoding="utf-8")
    assert "No reliable galaxies" in text
    assert "Total summaries loaded: **0**" in text


def test_audit_aggregates_morphology_and_correlations(tmp_path: Path):
    processed = tmp_path / "results"
    galaxy_root = processed / RESULTS_GALAXY_ANALYSIS_DIRNAME
    galaxy_root.mkdir(parents=True)

    # 3 reliable Im, 2 reliable Sm/Sd, 5 unreliable.
    for i in range(3):
        _write_summary(
            galaxy_root,
            f"IM{i:02d}",
            _make_payload(
                galaxy=f"IM{i:02d}",
                morphology="Dwarf Irregulars (Im)",
                validation_is_reliable=True,
                ratio=0.5 + 0.2 * i,
                galactic_lat=20.0 + 5.0 * i,
                inclination=55.0 + 5.0 * i,
                b_microgauss=1.0 + 0.5 * i,
                plasma_beta=10.0 - 2.0 * i,
                rar_residual=0.1 + 0.05 * i,
            ),
        )
    for i in range(2):
        _write_summary(
            galaxy_root,
            f"SD{i:02d}",
            _make_payload(
                galaxy=f"SD{i:02d}",
                morphology="Sm/Sd Spirals",
                validation_is_reliable=True,
                ratio=1.5 + 0.3 * i,
                galactic_lat=60.0 + 5.0 * i,
                inclination=65.0 + 5.0 * i,
                b_microgauss=3.0 + 0.5 * i,
                plasma_beta=2.0 + 0.5 * i,
                rar_residual=0.05 + 0.02 * i,
            ),
        )
    for i in range(5):
        _write_summary(
            galaxy_root,
            f"UNR{i:02d}",
            _make_payload(
                galaxy=f"UNR{i:02d}",
                morphology="Late-type Spirals (Sb, Sbc, Sc)",
                validation_is_reliable=False,
                ratio=0.001,
                galactic_lat=10.0,
                inclination=20.0,
                b_microgauss=0.5,
                plasma_beta=100.0,
                rar_residual=0.5,
            ),
        )

    out = write_consistency_audit_report(processed_root=processed)
    assert out.parent.name == "step_validation"
    text = out.read_text(encoding="utf-8")
    assert "Reliable galaxies: **5**" in text
    # Morphology table contains both groups.
    assert "Dwarf Irregulars (Im)" in text
    assert "Sm/Sd Spirals" in text
    # Foreground bias section + correlation table headers.
    assert "Foreground bias check" in text
    assert "Pearson" in text and "Spearman" in text


def test_audit_markdown_handles_no_reliable(tmp_path: Path):
    df = pd.DataFrame(
        [
            {
                "galaxy_name": "X",
                "morphology_class": "Dwarf",
                "hubble_t_type": 9,
                "distance_mpc": 3.0,
                "galactic_latitude_degrees": 10.0,
                "abs_galactic_latitude": 10.0,
                "validation_is_reliable": False,
                "calibration_status": "OPTIMAL",
                "ratio": 0.001,
                "pearson": 0.0,
                "magnetic_field_microgauss": 0.1,
                "magnetic_field_within_band": True,
                "plasma_beta_median": 100.0,
                "rar_residual_dex": 0.5,
                "density_clamped_fraction": 0.5,
                "stellar_dominated_fraction": 0.5,
                "inclination_deg": 30.0,
                "physical_neighbor_separation_kpc": 0.4,
            }
        ]
    )
    md = build_consistency_audit_markdown(df)
    assert "No reliable galaxies" in md
