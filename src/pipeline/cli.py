"""Command-line interface for batch galaxy processing.

The pipeline runs in a single unified flow (see :func:`src.pipeline.runner.run_pipeline`):

1. **Discovery census** — one engine pass per target galaxy at the physical Planck pixel scale
   with ``f_ion = 1.0``; writes ``discovery_all.csv``, ``discovery_reliable.csv``, ``discovery_not_reliable.csv``,
   and ``discovery_cohort_rho_band.csv`` under ``results/step_discovery/``.
2. **Automatic pre-filter** — galaxies with blocking Discovery audits (latitude, skipped physics,
   NED resolution) are excluded from the grid; all others enter the 9-scenario search.
3. **Validation grid** — full per-galaxy artefacts and population reports for the filtered set.

By default, ``python src/main.py`` runs **physics** (Phases 1–3: Discovery, Validation, population
Analysis report) then **visualization** (Phases 4–5: publication figures and supreme-cohort 3D + HTML).

Use ``--compute-only`` to run physics only, or ``--viz-only`` to regenerate figures and 3D from existing
``results/step_validation/`` tables (no SPARC / Planck pass).

Galaxy selection:

- ``--galaxies``: preset ``all`` (full SPARC catalogue, ~175), ``default`` (10 curated control
  galaxies), or one or more explicit galaxy names (e.g. ``NGC3741 DDO154``).
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="E-G galaxy rotation-curve analysis pipeline (SPARC + Planck).",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--compute-only",
        action="store_true",
        help="Run physics only (Discovery, Validation, Phase 3 population Analysis). Skip graphs and 3D.",
    )
    mode.add_argument(
        "--viz-only",
        action="store_true",
        help=(
            "Regenerate publication figures and supreme-cohort 3D + HTML from existing "
            "results/step_validation/*.csv (no SPARC / Planck physics pass)."
        ),
    )
    parser.add_argument(
        "--galaxies",
        nargs="+",
        default=["all"],
        metavar="ALL|DEFAULT|NAME",
        help=(
            "Which galaxies to process. Use the literal 'all' (full SPARC catalogue, ~175), "
            "'default' (10 curated control galaxies), or one or more explicit galaxy names "
            "(e.g. 'NGC3741 DDO154'). Presets cannot be mixed with names. Defaults to 'all'. "
            "Ignored when --viz-only is set (cohorts are read from disk)."
        ),
    )
    parser.add_argument(
        "--logs",
        action="store_true",
        help=(
            "Print per-galaxy provenance (SPARC / MORPH / NED / Planck, Unicode rule header), ICRS supplement "
            "notice, and verbose Planck map initialization chatter. "
            "Disables tqdm progress bars so those lines are not mixed with bar redraws."
        ),
    )
    return parser.parse_args()
