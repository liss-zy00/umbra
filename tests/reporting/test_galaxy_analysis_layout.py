"""Tests for ``galaxy_analysis`` tier relocation."""

from __future__ import annotations

import pandas as pd

from src.pipeline.constants import (
    GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME,
    RESULTS_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME,
)
from src.reporting.galaxy_analysis_layout import relocate_scientific_evidence_galaxy_dirs


def test_relocate_moves_reliable_folder_to_scientific_evidence(tmp_path) -> None:
    root = tmp_path
    ga = root / RESULTS_DIRNAME / RESULTS_GALAXY_ANALYSIS_DIRNAME
    rel = ga / GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME / "GALX1"
    rel.mkdir(parents=True)
    (rel / "GALX1_analysis.csv").write_text("radius_kpc\n1.0\n", encoding="utf-8")

    val = root / RESULTS_DIRNAME / RESULTS_STEP_VALIDATION_DIRNAME
    val.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"galaxy_name": "GALX1", "theoretical_gap_coverage_ratio": 1.0}]).to_csv(
        val / VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME, index=False
    )

    n = relocate_scientific_evidence_galaxy_dirs(root)
    assert n == 1
    dst = ga / GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME / "GALX1"
    assert dst.is_dir()
    assert (dst / "GALX1_analysis.csv").is_file()
    assert not rel.exists()
