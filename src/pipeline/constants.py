"""Filesystem locations and static galaxy lists for the SPARC + Planck pipeline."""

from __future__ import annotations

from pathlib import Path

# Single source of truth for pipeline artefact root (``<project_root>/<RESULTS_DIRNAME>/``).
RESULTS_DIRNAME = "results"

# Phase folders under ``results/`` use a common ``step_`` prefix so they sort adjacently.
RESULTS_STEP_DISCOVERY_DIRNAME = "step_discovery"
RESULTS_STEP_VALIDATION_DIRNAME = "step_validation"
# Discovery cohort CSV filenames under ``results/step_discovery/``.
DISCOVERY_RELIABLE_FILENAME = "discovery_reliable.csv"
# Publication / global figures root under ``results/``.
RESULTS_GRAPHS_DIRNAME = "graphs"
# Cohort panels from :mod:`src.reporting.publication_bundle` (flat under ``results/graphs/``).
RESULTS_GRAPHS_POPULATION_DIRNAME = "population"
# Per-galaxy Phase-3 artefacts: ``results/galaxy_analysis/{tier}/{name}/`` (CSV, JSON, optional PNGs).
RESULTS_GALAXY_ANALYSIS_DIRNAME = "galaxy_analysis"
GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME = "validation_reliable"
# Mirrors ``validation_not_reliable.csv`` (``validation_is_reliable == False``); CSV+JSON only in normal runs.
GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME = "validation_not_reliable"
# Populated after publication figures: galaxies listed in ``validation_scientific_evidence_report.csv``.
GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME = "validation_scientific_evidence"

# Validation-phase population artefacts (all under ``results/step_validation/``).
VALIDATION_CATALOG_ALL_FILENAME = "validation_all.csv"
# Same inventory schema as ``validation_all.csv``, split on ``validation_is_reliable`` for convenience.
VALIDATION_RELIABLE_FILENAME = "validation_reliable.csv"
VALIDATION_NOT_RELIABLE_FILENAME = "validation_not_reliable.csv"
# One row per galaxy in the **supreme** cohort (validation-reliable + inclination gate): best grid
# row merged with Discovery headline columns (rho_req, beam diagnostics, discovery_is_elite, MC bands, …).
VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME = "validation_scientific_evidence_report.csv"
VALIDATION_ANCHOR_SENSITIVITY_SWEEP_FILENAME = "validation_anchor_sensitivity_sweep.csv"
# Markdown population audit (re-reads ``galaxy_analysis/*/*/*_summary.json`` after the run).
VALIDATION_CONSISTENCY_AUDIT_FILENAME = "validation_consistency_audit.md"

# Lelli+2016c SPARC ``*_rotmod.dat`` collection (one file per galaxy): Newtonian mass-model
# rotation curves (radius, velocities, surface densities). Header carries the SPARC distance.
SPARC_DIR = Path("sources") / "sparc" / "rotation_curves"
SPARC_CATALOG_JSON = Path("sources") / "sparc" / "SPARC_Lelli2016c.json"
# Disk PA overlay (HyperLEDA, NED, …) keyed by SPARC galaxy name — not part of the Lelli+2016c drop.
MORPHOLOGY_DISK_PA_OVERLAY_JSON = Path("sources") / "morphology" / "disk_position_angle_overlay.json"
# Optional ICRS seeds for NED cone search when name resolution fails (consumed via NedProvider).
NED_ICRS_CONE_SEEDS_JSON = Path("sources") / "ned" / "ned_icrs_cone_seeds.json"

# Extragalactic Distance Database (EDD) — HyperLEDA-formatted CSV drops
# (see ``sources/extragalactic_distances/README.md``).
EXTRAGALACTIC_DISTANCES_DIR = Path("sources") / "extragalactic_distances"
EDD_CMD_TRGB_CSV = EXTRAGALACTIC_DISTANCES_DIR / "edd_cmd_trgb.csv"
EDD_CF4_INDIVIDUAL_CSV = EXTRAGALACTIC_DISTANCES_DIR / "edd_cf4_individual_distances.csv"
EDD_COSMICFLOWS3_CSV = EXTRAGALACTIC_DISTANCES_DIR / "edd_cosmicflows3_distances.csv"
SPARC_NAME_TO_PGC_JSON = EXTRAGALACTIC_DISTANCES_DIR / "sparc_name_to_pgc.json"

# Display label for the photometric-LSB pseudo-class used in the morphology comparison plot.
# This is *not* a SPARC Hubble bin; it tags any galaxy whose disk surface brightness sits below
# ``RELIABILITY_CUTS.PHOTOMETRIC_LSB_SB_THRESHOLD`` (regardless of its T-type), so the figure can
# show the historical "LSB galaxies vs the rest" contrast alongside the standard Hubble bins.
PHOTOMETRIC_LSB_LABEL = "Photometric LSB (SBdisk < 100)"

FITS_CANDIDATES = [
    Path("sources") / "planck" / "COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits",
    Path("sources") / "planck" / "COM_CompMap_QU-thermaldust-commander_2048_R3.00_full.fits",
]

DEFAULT_CONTROL_GALAXIES = [
    "NGC5005",
    "NGC2841",
    "NGC3198",
    "DDO154",
    "UGC4483",
    "IC2574",
    "NGC6946",
    "NGC5907",
    "NGC3741",
    "F563-V2",
]

# E2E regression audit pair (project policy). Every name must appear in ``DEFAULT_CONTROL_GALAXIES``.
E2E_AUDIT_CONTROL_GALAXIES: tuple[str, ...] = ("DDO154", "NGC3741")

for _e2e_name in E2E_AUDIT_CONTROL_GALAXIES:
    if _e2e_name not in DEFAULT_CONTROL_GALAXIES:
        raise ValueError(
            f"E2E_AUDIT_CONTROL_GALAXIES entry {_e2e_name!r} must be listed in DEFAULT_CONTROL_GALAXIES"
        )
del _e2e_name
