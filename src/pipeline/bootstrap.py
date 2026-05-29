"""
Composition root for pipeline services (Factory-style wiring of providers and reporters).

Keeps :mod:`src.pipeline.runner` focused on phase orchestration and CLI concerns.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np

from src.pipeline.constants import (
    EDD_CF4_INDIVIDUAL_CSV,
    EDD_CMD_TRGB_CSV,
    EDD_COSMICFLOWS3_CSV,
    FITS_CANDIDATES,
    MORPHOLOGY_DISK_PA_OVERLAY_JSON,
    NED_ICRS_CONE_SEEDS_JSON,
    RESULTS_DIRNAME,
    SPARC_CATALOG_JSON,
    SPARC_NAME_TO_PGC_JSON,
)
from src.pipeline.per_galaxy.galaxy_processor import PipelineDependencies
from src.providers import NedProvider, PlanckProvider, SparcMorphologyProvider, SparcProvider
from src.providers.edd_provider import EddProvider, build_sparc_name_to_pgc_lookup

logger = logging.getLogger(__name__)


def load_icrs_supplement(path: Path) -> dict[str, tuple[float, float]]:
    """
    Read ``NED_ICRS_CONE_SEEDS_JSON`` (see :mod:`src.pipeline.constants`) and return a flat
    ``{galaxy_name -> (ra_deg, dec_deg)}`` map.
    """
    if not path.is_file():
        return {}
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Could not read ICRS supplement %s: %s", path, exc)
        return {}

    positions = raw.get("positions") if isinstance(raw, dict) else None
    if not isinstance(positions, dict):
        return {}

    out: dict[str, tuple[float, float]] = {}
    for name, value in positions.items():
        if str(name).startswith("_") or not isinstance(value, dict):
            continue
        try:
            ra = float(value["ra_deg"])
            dec = float(value["dec_deg"])
        except (KeyError, TypeError, ValueError):
            continue
        if np.isfinite(ra) and np.isfinite(dec):
            out[str(name).strip()] = (ra, dec)
    return out


def build_pipeline_dependencies(
    root: Path,
    *,
    verbose_planck: bool = False,
    emit_verbose_galaxy_logs: bool = False,
) -> PipelineDependencies:
    """
    Construct :class:`PipelineDependencies` with initialized Planck maps and morphology tables.

    ``deps.fits_summary`` holds the dict from :meth:`PlanckProvider.initialize` for downstream
    provenance on global artefacts.

    ``verbose_planck`` controls Planck map load chatter; ``emit_verbose_galaxy_logs`` gates the
    per-galaxy provenance terminal block (Unicode rule line + SPARC/MORPH/EDD/NED/PLANCK rows; see ``--logs``).
    """
    from src.reporting.persistence.exporter import ResultPersister

    processed_root = root / RESULTS_DIRNAME

    sparc = SparcProvider()
    morphology_provider = SparcMorphologyProvider(
        root / SPARC_CATALOG_JSON,
        pa_supplement_path=root / MORPHOLOGY_DISK_PA_OVERLAY_JSON,
    )
    morphology_lookup = morphology_provider.load_morphology_lookup()
    ned = NedProvider()
    planck = PlanckProvider([root / rel for rel in FITS_CANDIDATES], verbose=verbose_planck)
    persister = ResultPersister(processed_root)

    fits_summary = planck.initialize()

    icrs_supplement = load_icrs_supplement(root / NED_ICRS_CONE_SEEDS_JSON)

    edd: EddProvider | None = None
    try:
        if (root / EDD_CMD_TRGB_CSV).is_file():
            edd = EddProvider.from_csv_paths(
                root / EDD_CMD_TRGB_CSV,
                root / EDD_CF4_INDIVIDUAL_CSV,
                root / EDD_COSMICFLOWS3_CSV,
            )
    except (OSError, ValueError) as exc:
        logger.warning("EDD tables could not be loaded: %s", exc)

    try:
        pgc_by_sparc_name = build_sparc_name_to_pgc_lookup(
            root / SPARC_NAME_TO_PGC_JSON,
            root / EDD_COSMICFLOWS3_CSV,
            morphology_provider,
        )
    except (OSError, ValueError) as exc:
        logger.warning("PGC crosswalk build failed: %s", exc)
        pgc_by_sparc_name = {}

    deps = PipelineDependencies(
        sparc=sparc,
        morphology_provider=morphology_provider,
        morphology_lookup=morphology_lookup,
        ned=ned,
        planck=planck,
        persister=persister,
        fits_summary=fits_summary,
        icrs_supplement=icrs_supplement,
        edd=edd,
        pgc_by_sparc_name=pgc_by_sparc_name,
        emit_verbose_galaxy_logs=emit_verbose_galaxy_logs,
    )
    return deps
