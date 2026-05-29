import pandas as pd

from src.pipeline.constants import VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME
from src.reporting.persistence.raw_data import DataPersister


def test_save_scientific_evidence_report_columns(tmp_path):
    persister = DataPersister(tmp_path)
    rows = [
        {
            "galaxy_name": "G1",
            "right_ascension_degrees": 1.0,
            "declination_degrees": 2.0,
            "best_scenario_name": "s",
            "scenario_gradient_scale_kpc": 0.1,
            "scenario_ionization_fraction": 0.01,
            "theoretical_gap_coverage_ratio": 1.0,
            "pearson_correlation": 0.5,
            "morphology_class": "Im",
            "ratio_distance_to_unity": 0.0,
        }
    ]
    path = persister.save_scientific_evidence_report(rows)
    assert path.parent.name == "step_validation"
    assert path.name == VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME
    df = pd.read_csv(path)
    assert list(df.columns) == [
        "galaxy_name",
        "right_ascension_degrees",
        "declination_degrees",
        "best_scenario_name",
        "scenario_gradient_scale_kpc",
        "scenario_ionization_fraction",
        "theoretical_gap_coverage_ratio",
        "pearson_correlation",
        "morphology_class",
        "ratio_distance_to_unity",
    ]


def test_save_text_file(tmp_path):
    p = DataPersister(tmp_path)
    out = p.save_text_file("note.md", "# hi")
    assert out.read_text(encoding="utf-8") == "# hi"
