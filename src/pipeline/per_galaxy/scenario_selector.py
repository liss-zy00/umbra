"""Pick the preferred (scale, ionization) scenario from per-galaxy metrics."""

from __future__ import annotations

from typing import Any

import numpy as np

from src.config.thresholds import RELIABILITY_CUTS


def _log_distance_to_unity(value: float) -> float:
    return abs(
        np.log10(
            max(
                float(value),
                RELIABILITY_CUTS.BEST_SCENARIO_LOG_FLOOR,
            )
        )
    )


def select_best_scenario(all_scenario_metrics: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Among metrics marked ``is_reliable`` with finite positive ``calibrated_convergence_ratio``,
    minimize ``|log10(max(ratio, floor))|`` so the beam-corrected required density tracks the WIM
    anchor (Planck dilution discounted via ``k``).

    Falls back to reliable scenarios with finite ``theoretical_gap_coverage_ratio`` (legacy gap
    amplitude screen), then to the first grid cell.
    """
    valid_calibrated = [
        m
        for m in all_scenario_metrics
        if m["is_reliable"]
        and np.isfinite(float(m.get("calibrated_convergence_ratio", float("nan"))))
        and float(m["calibrated_convergence_ratio"]) > 0.0
    ]
    if valid_calibrated:
        return min(valid_calibrated, key=lambda m: _log_distance_to_unity(float(m["calibrated_convergence_ratio"])))

    valid_gap = [
        m
        for m in all_scenario_metrics
        if m["is_reliable"] and np.isfinite(m["theoretical_gap_coverage_ratio"])
    ]
    if valid_gap:
        return min(
            valid_gap,
            key=lambda m: _log_distance_to_unity(float(m["theoretical_gap_coverage_ratio"])),
        )
    return all_scenario_metrics[0]
