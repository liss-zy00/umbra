"""
Shared pytest fixtures: synthetic galaxies, SPARC paths, Planck FITS, and patched NED.

All file-producing fixtures write under ``tmp_path`` only — never under ``sources/`` or ``results/``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from src.models.galaxy import Galaxy


@pytest.fixture
def fixtures_dir() -> Path:
    return Path(__file__).resolve().parent / "fixtures"


@pytest.fixture
def sample_rotmod_path(fixtures_dir: Path) -> Path:
    return fixtures_dir / "SYN001_rotmod.dat"


@pytest.fixture
def sample_morphology_json(fixtures_dir: Path) -> Path:
    return fixtures_dir / "sample_morphology.json"


@pytest.fixture
def synthetic_planck_fits(tmp_path: Path) -> Path:
    from tests.fixtures.make_synthetic_fits import write_polarization_fits

    return write_polarization_fits(tmp_path / "planck_synthetic.fits", nside=16, ordering="RING")


@pytest.fixture
def minimal_rotation_dataframe() -> pd.DataFrame:
    """Single-radius kinematic row with all velocities needed by ``PhysicsEngine``."""
    return pd.DataFrame(
        {
            "radius_kpc": [1.0],
            "observed_velocity_kms": [100.0],
            "baryonic_velocity_kms": [50.0],
            "gas_velocity_kms": [10.0],
            "disk_velocity_kms": [40.0],
            "bulge_velocity_kms": [0.0],
        }
    )


@pytest.fixture
def minimal_galaxy(minimal_rotation_dataframe: pd.DataFrame) -> Galaxy:
    return Galaxy(name="TestGal", distance_mpc=10.0, rotation_dataframe=minimal_rotation_dataframe)


@pytest.fixture
def mock_ned_success(monkeypatch: pytest.MonkeyPatch) -> None:
    """Patch ``NedProvider.fetch_metadata`` to avoid HTTP calls."""

    def fake_fetch(
        self: Any,
        galaxy_name: str,
        *,
        icrs_deg: tuple[float, float] | None = None,
    ) -> dict[str, Any]:
        _ = galaxy_name, icrs_deg
        return {
            "right_ascension_degrees": 180.0,
            "declination_degrees": 45.0,
            "redshift": 0.001,
            "distance_mpc": 10.0,
            "distance_source": "mock",
            "object_type": "G",
            "raw_status": "ok",
        }

    from src.providers import ned as ned_mod

    monkeypatch.setattr(ned_mod.NedProvider, "fetch_metadata", fake_fetch)
