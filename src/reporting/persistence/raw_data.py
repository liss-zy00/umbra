from __future__ import annotations

"""
CSV / JSON persistence layer for the per-galaxy and per-grid scenario records.

The :class:`DataPersister` is a thin façade over :mod:`pandas` and :mod:`json` that exposes
the directory contract used by the rest of the pipeline:

    results/
    ├── galaxy_analysis/  ─ ``{validation_reliable|validation_not_reliable}/<NAME>/<NAME>_analysis.csv``
    │                       and ``<NAME>_summary.json`` (structured summary + audit inputs for
    │                       :mod:`src.pipeline.aggregation.audit`)
    ├── grid_results/     ─ one CSV per scenario row of the magnetic-field grid
    ├── step_validation/  ─ population-level CSV / Markdown synthesis (see :class:`ResultPersister`)
    └── graphs/           ─ figures produced by :mod:`src.reporting.plotting`

The class is intentionally **decoupled from any plotting or narrative concern** (Single
Responsibility Principle): figure writers live in :mod:`src.reporting.plotting`, narrative
generation in :mod:`src.reporting.narrative`, and CSV/JSON synthesis in
:mod:`src.reporting.persistence.exporter`. The shared persistence layer here covers raw
per-galaxy artifacts only — high-level synthesis tables flow through the dedicated
:class:`src.reporting.persistence.exporter.ResultPersister`.
"""

import json
from pathlib import Path
from typing import Any

import pandas as pd

from src.pipeline.constants import (
    GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME,
)
from src.models.galaxy import Galaxy
from src.reporting.persistence.exporter import _sanitize_for_json


class DataPersister:
    """
    Single-responsibility writer for the per-galaxy raw-data artifacts.

    The persister never inspects the scientific content; it serialises whatever the
    pipeline gives it. Plotting, narrative composition and population-level synthesis
    each live in their own SRP-compliant module under :mod:`src.reporting`.
    """

    def __init__(self, processed_root: Path) -> None:
        self.processed_root = processed_root
        self.individual_root = self.processed_root / RESULTS_GALAXY_ANALYSIS_DIRNAME
        self.grid_root = self.processed_root / "grid_results"
        self.validation_root = self.processed_root / RESULTS_STEP_VALIDATION_DIRNAME
        self.graphs_root = self.processed_root / "graphs"

        self.processed_root.mkdir(parents=True, exist_ok=True)

    def save_galaxy(self, galaxy: Galaxy, summary: dict[str, Any]) -> tuple[Path, Path, Path]:
        """Saves the specific rotation curve and metadata for a single galaxy."""
        self.individual_root.mkdir(parents=True, exist_ok=True)
        dra = summary.get("data_reliability_assessment")
        reliable = bool((dra or {}).get("validation_is_reliable", False))
        tier = GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME if reliable else GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME
        galaxy_dir = self.individual_root / tier / galaxy.name
        galaxy_dir.mkdir(parents=True, exist_ok=True)

        csv_path = galaxy_dir / f"{galaxy.name}_analysis.csv"
        json_path = galaxy_dir / f"{galaxy.name}_summary.json"

        galaxy.rotation_dataframe.to_csv(csv_path, index=False)
        json_path.write_text(
            json.dumps(
                _sanitize_for_json(summary),
                indent=2,
                ensure_ascii=True,
                allow_nan=False,
            ),
            encoding="utf-8",
        )

        return galaxy_dir, csv_path, json_path

    def _save_report(self, filename: str, rows: list[dict[str, Any]], columns: list[str]) -> Path:
        """Internal helper to save population-level CSV tables under ``step_validation/``."""
        self.validation_root.mkdir(parents=True, exist_ok=True)
        output_path = self.validation_root / filename
        if not rows:
            pd.DataFrame(columns=columns).to_csv(output_path, index=False)
        else:
            pd.DataFrame(rows).reindex(columns=columns).to_csv(output_path, index=False)
        return output_path

    def save_dataframe_report(self, filename: str, rows: list[dict[str, Any]], columns: list[str]) -> Path:
        """Public wrapper for typed CSV exports (inventory / synthesis tables)."""
        return self._save_report(filename, rows, columns)

    def save_scientific_evidence_report(self, rows: list[dict[str, Any]]) -> Path:
        """Saves the final filtered list of galaxies that successfully match the theoretical model."""
        columns = [
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
        return self._save_report(VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME, rows, columns)

    def save_text_file(self, filename: str, content: str) -> Path:
        """Write Markdown or plain text under ``step_validation/`` (thin helper used by tests)."""
        self.validation_root.mkdir(parents=True, exist_ok=True)
        output_path = self.validation_root / filename
        output_path.write_text(content, encoding="utf-8")
        return output_path
