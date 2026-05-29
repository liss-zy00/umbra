from __future__ import annotations

"""
Domain model: a single SPARC galaxy as a typed mutable record.

Carries provider-supplied catalog metadata (catalogue name, distance, ICRS coordinates,
inclination, asymptotic ``V_flat``) plus the per-shell rotation table as a pandas
``DataFrame``. The :class:`Galaxy` is intentionally agnostic of physics or persistence —
estimators in :mod:`src.physics` and writers in :mod:`src.reporting` operate on the
:class:`Galaxy` instance through composition, never through inheritance.

Single source of truth for column conventions
---------------------------------------------
The ``rotation_dataframe`` holds SPARC kinematic inputs in km/s and kpc. After the engine
runs, additional SI-unit columns (``radius_m``, ``observed_velocity_ms``,
``centripetal_acceleration``, ``acceleration_gap``, ``magnetic_pressure_acceleration`` …)
are appended. The conceptual / mathematical meaning of each derived quantity is documented
in :doc:`EG_MODEL_THEORY` (§3, *Mathematical Implementation*). Operational column listings
live in the runtime ``README.md`` so this domain model stays implementation-agnostic.
"""

from dataclasses import dataclass, field
from typing import Any

import pandas as pd


@dataclass
class Galaxy:
    """
    Mutable kinematic record for one SPARC object.

    Attributes
    ----------
    name:
        Catalog identifier matching the ``*_rotmod.dat`` stem (without suffix).
    distance_mpc:
        Luminosity-distance anchor from the SPARC header (NED may override metadata context only).
    right_ascension_degrees / declination_degrees:
        ICRS coordinates when NED resolves the target (required for Planck sampling).
    redshift:
        Optional heliocentric redshift from NED.
    inclination_deg / inclination_uncertainty_deg:
        Disk inclination from the Lelli+2016 catalog and its 1-σ error. SPARC rotmod V_obs is
        already deprojected, so the inclination error becomes the dominant kinematic systematic
        propagated downstream by ``SparcProvider.attach_inclination_uncertainty``.
    vflat_kms / vflat_uncertainty_kms:
        Catalogued asymptotic flat rotation velocity and its error (used by RAR/BTFR diagnostics).
    distance_mpc / distance_uncertainty_mpc:
        SPARC luminosity distance and its 1-σ catalogued error ``e_D_Mpc``. The uncertainty is
        consumed by the Monte-Carlo band on ``theoretical_gap_coverage_ratio`` (see
        :mod:`src.physics.diagnostics.uncertainty`).
    position_angle_deg / position_angle_uncertainty_deg:
        Disk major-axis position angle (degrees, east-of-north) and its 1-σ error. SPARC does
        not publish a PA; values are loaded from ``sources/morphology/disk_position_angle_overlay.json``
        (HyperLEDA / dynamical-fit curation). When PA is missing the engine falls back to the
        axisymmetric gradient and emits ``warn_no_pa``.
    rotation_dataframe:
        Per-radius SPARC columns plus derived SI quantities after processing.
    metadata:
        Arbitrary auxiliary payload (debug hooks, provider blobs).
    """

    name: str
    distance_mpc: float | None = None
    distance_uncertainty_mpc: float | None = None
    right_ascension_degrees: float | None = None
    declination_degrees: float | None = None
    redshift: float | None = None
    inclination_deg: float | None = None
    inclination_uncertainty_deg: float | None = None
    position_angle_deg: float | None = None
    position_angle_uncertainty_deg: float | None = None
    vflat_kms: float | None = None
    vflat_uncertainty_kms: float | None = None
    rotation_dataframe: pd.DataFrame = field(default_factory=pd.DataFrame)
    metadata: dict[str, Any] = field(default_factory=dict)
