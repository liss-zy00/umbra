from src.pipeline.constants import RESULTS_STEP_DISCOVERY_DIRNAME, RESULTS_STEP_VALIDATION_DIRNAME
from src.reporting.persistence.exporter import ResultPersister


def test_result_persister_grid_and_inventory(tmp_path):
    rp = ResultPersister(tmp_path)
    grid_path = rp.save_grid_scenario_report("scenario_x", [{"galaxy_name": "A", "k": 1}])
    assert grid_path.exists()
    inv_path = rp.save_full_catalog_inventory([])
    assert inv_path.exists()


def test_result_persister_does_not_create_validation_subdirs_on_init(tmp_path):
    """Discovery CSV writers only need ``processed_root``; subfolders are created lazily."""
    proc = tmp_path / "results"
    ResultPersister(proc)
    assert proc.is_dir()
    assert not (proc / "galaxy_analysis").exists()
    assert not (proc / "grid_results").exists()
    assert not (proc / "reports").exists()
    assert not (proc / "graphs").exists()
    assert not (proc / RESULTS_STEP_DISCOVERY_DIRNAME).exists()
    assert not (proc / RESULTS_STEP_VALIDATION_DIRNAME).exists()


def test_save_galaxy_summary_json_is_strict_and_nulls_non_finite_floats(tmp_path):
    import json

    import pandas as pd

    from src.models.galaxy import Galaxy

    proc = tmp_path / "results"
    rp = ResultPersister(proc)
    df = pd.DataFrame({"radius_kpc": [1.0]})
    g = Galaxy(name="T1", distance_mpc=1.0, rotation_dataframe=df)
    summary = {
        "data_reliability_assessment": {"validation_is_reliable": False},
        "theoretical_gap_coverage_ratio": float("nan"),
        "nested": {"x": 1.0},
    }
    _, _, json_path = rp.save_galaxy(g, summary)
    text = json_path.read_text(encoding="utf-8")
    assert "NaN" not in text
    payload = json.loads(text)
    assert payload["theoretical_gap_coverage_ratio"] is None
    assert payload["nested"]["x"] == 1.0

