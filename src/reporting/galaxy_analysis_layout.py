"""Filesystem layout helpers for ``results/galaxy_analysis/`` tier directories."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pandas as pd

from src.pipeline.constants import (
    GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME,
    RESULTS_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME,
)

logger = logging.getLogger(__name__)


def relocate_scientific_evidence_galaxy_dirs(project_root: Path) -> int:
    """
    Move ``galaxy_analysis/validation_reliable/<NAME>/`` into
    ``galaxy_analysis/validation_scientific_evidence/<NAME>/`` for every ``galaxy_name`` listed in
    ``validation_scientific_evidence_report.csv``.

    Keeps ``validation_reliable/`` focused on reliable galaxies that are *not* in the merged
    direct-evidence cohort; the scientific-evidence tier holds the same CSV/JSON/PNG bundle for
    the highlighted systems. Idempotent when a destination already exists (replaced).
    """
    evidence_csv = (
        project_root
        / RESULTS_DIRNAME
        / RESULTS_STEP_VALIDATION_DIRNAME
        / VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME
    )
    if not evidence_csv.is_file():
        return 0
    try:
        df = pd.read_csv(evidence_csv)
    except (OSError, ValueError) as exc:
        logger.warning("Could not read scientific evidence report %s: %s", evidence_csv, exc)
        return 0
    if "galaxy_name" not in df.columns:
        return 0
    names = [str(x).strip() for x in df["galaxy_name"].tolist() if str(x).strip()]
    if not names:
        return 0
    root = project_root / RESULTS_DIRNAME / RESULTS_GALAXY_ANALYSIS_DIRNAME
    src_tier = root / GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME
    dst_tier = root / GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME
    dst_tier.mkdir(parents=True, exist_ok=True)
    moved = 0
    for name in names:
        src = src_tier / name
        if not src.is_dir():
            continue
        dst = dst_tier / name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.move(str(src), str(dst))
        moved += 1
    if moved:
        logger.info(
            "Relocated %d galaxy artefact folder(s) to %s/",
            moved,
            GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME,
        )
    return moved
