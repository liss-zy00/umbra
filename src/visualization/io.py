"""Load per-galaxy pipeline artefacts for visualization."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from src.reporting.galaxy_from_disk import analysis_csv_path, load_summary_json, summary_json_path


def load_analysis_dataframe(project_root: Path, galaxy_name: str) -> pd.DataFrame:
    """Return the per-shell analysis table for ``galaxy_name``."""
    path = analysis_csv_path(project_root, galaxy_name)
    return pd.read_csv(path)


def load_optional_summary(project_root: Path, galaxy_name: str) -> dict[str, Any] | None:
    """Best-effort JSON summary next to the analysis CSV (may be absent in tests)."""
    return load_summary_json(summary_json_path(project_root, galaxy_name))
