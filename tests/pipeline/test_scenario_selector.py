import numpy as np

from src.pipeline.per_galaxy.scenario_selector import select_best_scenario


def test_select_best_prefers_calibrated_convergence_near_one():
    metrics = [
        {"is_reliable": True, "theoretical_gap_coverage_ratio": 1.0, "calibrated_convergence_ratio": 50.0},
        {"is_reliable": True, "theoretical_gap_coverage_ratio": 10.0, "calibrated_convergence_ratio": 1.05},
        {"is_reliable": False, "theoretical_gap_coverage_ratio": 1.0, "calibrated_convergence_ratio": 1.0},
    ]
    best = select_best_scenario(metrics)
    assert np.isclose(best["calibrated_convergence_ratio"], 1.05)


def test_select_best_falls_back_to_gap_ratio_when_calibrated_missing():
    metrics = [
        {"is_reliable": True, "theoretical_gap_coverage_ratio": 10.0, "calibrated_convergence_ratio": np.nan},
        {"is_reliable": True, "theoretical_gap_coverage_ratio": 1.05, "calibrated_convergence_ratio": np.nan},
    ]
    best = select_best_scenario(metrics)
    assert np.isclose(best["theoretical_gap_coverage_ratio"], 1.05)


def test_select_fallback_when_none_reliable():
    metrics = [{"is_reliable": False, "theoretical_gap_coverage_ratio": np.nan, "calibrated_convergence_ratio": np.nan}]
    assert select_best_scenario(metrics) is metrics[0]
