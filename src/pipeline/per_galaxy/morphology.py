"""Morphology catalog lookups and scientific-evidence inclusion rules."""

from __future__ import annotations

from typing import Any

import numpy as np

from src.providers.sparc import SparcMorphologyProvider

from src.config.thresholds import RELIABILITY_CUTS


def catalog_row_for_galaxy(
    galaxy_name: str,
    morphology_lookup: dict[str, dict[str, Any]],
    morphology_provider: SparcMorphologyProvider,
) -> dict[str, Any]:
    key = morphology_provider.normalize_galaxy_name(galaxy_name)
    return morphology_lookup.get(key) or morphology_lookup.get(galaxy_name.replace(" ", "").upper(), {})


def is_photometric_lsb(
    row: dict[str, Any],
    morphology_lookup: dict[str, dict[str, Any]],
    morphology_provider: SparcMorphologyProvider,
) -> bool:
    entry = catalog_row_for_galaxy(str(row.get("galaxy_name", "")), morphology_lookup, morphology_provider)
    sb = entry.get("disk_surface_brightness_solar_per_pc2")
    if sb is None:
        return False
    sbf = float(sb)
    return bool(np.isfinite(sbf) and sbf < RELIABILITY_CUTS.PHOTOMETRIC_LSB_SB_THRESHOLD)


def include_in_scientific_evidence_report(
    row: dict[str, Any],
    morphology_lookup: dict[str, dict[str, Any]],
    morphology_provider: SparcMorphologyProvider,
) -> bool:
    if row.get("morphology_class") == "Dwarf Irregulars (Im)":
        return True
    return is_photometric_lsb(row, morphology_lookup, morphology_provider)
