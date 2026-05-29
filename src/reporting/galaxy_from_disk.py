"""Rehydrate a minimal :class:`~src.models.galaxy.Galaxy` from on-disk pipeline artefacts."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from src.models.galaxy import Galaxy
from src.pipeline.constants import (
    GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME,
    RESULTS_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
)

logger = logging.getLogger(__name__)


def galaxy_analysis_root(project_root: Path) -> Path:
    return project_root / RESULTS_DIRNAME / RESULTS_GALAXY_ANALYSIS_DIRNAME


def resolve_per_galaxy_analysis_dir(project_root: Path, galaxy_name: str) -> Path:
    """
    Return the directory that holds ``{galaxy_name}_analysis.csv``.

    Checks the scientific-evidence tier first (post-run relocation), then validation-reliable,
    then the not-reliable grid cohort.
    """
    root = galaxy_analysis_root(project_root)
    csv_name = f"{galaxy_name}_analysis.csv"
    for tier in (
        GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME,
        GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
        GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME,
    ):
        candidate_dir = root / tier / galaxy_name
        if (candidate_dir / csv_name).is_file():
            return candidate_dir
    raise FileNotFoundError(
        f"Missing {csv_name} under {root}/{{{GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME},"
        f"{GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME},{GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME}}}/{galaxy_name}/"
    )


def analysis_csv_path(project_root: Path, galaxy_name: str) -> Path:
    return resolve_per_galaxy_analysis_dir(project_root, galaxy_name) / f"{galaxy_name}_analysis.csv"


def summary_json_path(project_root: Path, galaxy_name: str) -> Path:
    return resolve_per_galaxy_analysis_dir(project_root, galaxy_name) / f"{galaxy_name}_summary.json"


def load_summary_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Could not read summary JSON %s: %s", path, exc)
        return None


def build_galaxy_from_galaxy_logs(project_root: Path, galaxy_name: str) -> tuple[Galaxy, dict[str, Any] | None]:
    """
    Load ``*_analysis.csv`` and optional ``*_summary.json`` under
    ``results/galaxy_analysis/{tier}/{name}/``.
    """
    csv_path = analysis_csv_path(project_root, galaxy_name)
    if not csv_path.is_file():
        raise FileNotFoundError(f"Missing analysis CSV: {csv_path}")

    df = pd.read_csv(csv_path)
    json_path = summary_json_path(project_root, galaxy_name)
    summary = load_summary_json(json_path)

    distance_mpc: float | None = None
    ra: float | None = None
    dec: float | None = None
    redshift: float | None = None
    inclination_deg: float | None = None
    inclination_uncertainty_deg: float | None = None
    metadata: dict[str, Any] = {}
    vflat_kms: float | None = None
    vflat_uncertainty_kms: float | None = None
    position_angle_deg: float | None = None
    position_angle_uncertainty_deg: float | None = None
    distance_uncertainty_mpc: float | None = None

    if summary:
        ident = summary.get("galaxy_identity") or {}
        astro = summary.get("astrophysical_metadata") or {}
        if ident.get("morphology_class") is not None:
            metadata["morphology_class"] = ident.get("morphology_class")
        if ident.get("hubble_t_type") is not None:
            metadata["hubble_t_type"] = ident.get("hubble_t_type")
        dm = ident.get("distance_local_mpc")
        if dm is not None:
            try:
                distance_mpc = float(dm)
            except (TypeError, ValueError):
                distance_mpc = None
        if distance_mpc is None and astro.get("distance_mpc") is not None:
            try:
                distance_mpc = float(astro["distance_mpc"])
            except (TypeError, ValueError):
                pass
        if astro.get("right_ascension_degrees") is not None:
            try:
                ra = float(astro["right_ascension_degrees"])
            except (TypeError, ValueError):
                ra = None
        if astro.get("declination_degrees") is not None:
            try:
                dec = float(astro["declination_degrees"])
            except (TypeError, ValueError):
                dec = None
        if astro.get("redshift") is not None:
            try:
                redshift = float(astro["redshift"])
            except (TypeError, ValueError):
                redshift = None
        kin = summary.get("kinematic_metadata") or {}
        if kin.get("inclination_deg") is not None:
            try:
                inclination_deg = float(kin["inclination_deg"])
            except (TypeError, ValueError):
                inclination_deg = None
        if kin.get("inclination_uncertainty_deg") is not None:
            try:
                inclination_uncertainty_deg = float(kin["inclination_uncertainty_deg"])
            except (TypeError, ValueError):
                inclination_uncertainty_deg = None
        if kin.get("vflat_kms") is not None:
            try:
                vflat_kms = float(kin["vflat_kms"])
            except (TypeError, ValueError):
                vflat_kms = None
        if kin.get("vflat_uncertainty_kms") is not None:
            try:
                vflat_uncertainty_kms = float(kin["vflat_uncertainty_kms"])
            except (TypeError, ValueError):
                vflat_uncertainty_kms = None
        if kin.get("position_angle_deg") is not None:
            try:
                position_angle_deg = float(kin["position_angle_deg"])
            except (TypeError, ValueError):
                position_angle_deg = None
        if kin.get("position_angle_uncertainty_deg") is not None:
            try:
                position_angle_uncertainty_deg = float(kin["position_angle_uncertainty_deg"])
            except (TypeError, ValueError):
                position_angle_uncertainty_deg = None
        dist_unc = kin.get("distance_uncertainty_mpc")
        if dist_unc is not None:
            try:
                distance_uncertainty_mpc = float(dist_unc)
            except (TypeError, ValueError):
                distance_uncertainty_mpc = None

    if distance_mpc is None or not (distance_mpc > 0.0):
        distance_mpc = 1.0

    return (
        Galaxy(
            name=galaxy_name,
            distance_mpc=float(distance_mpc),
            distance_uncertainty_mpc=distance_uncertainty_mpc,
            rotation_dataframe=df,
            right_ascension_degrees=ra,
            declination_degrees=dec,
            redshift=redshift,
            inclination_deg=inclination_deg,
            inclination_uncertainty_deg=inclination_uncertainty_deg,
            position_angle_deg=position_angle_deg,
            position_angle_uncertainty_deg=position_angle_uncertainty_deg,
            vflat_kms=vflat_kms,
            vflat_uncertainty_kms=vflat_uncertainty_kms,
            metadata=metadata,
        ),
        summary,
    )
