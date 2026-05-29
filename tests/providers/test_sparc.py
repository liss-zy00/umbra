import json
from pathlib import Path

import numpy as np

from src.providers.sparc import (
    SparcMorphologyProvider,
    SparcProvider,
    _load_pa_supplement,
)


def test_load_rotation_curve_from_fixture(sample_rotmod_path):
    sparc = SparcProvider()
    df, d_mpc = sparc.load_rotation_curve(sample_rotmod_path)
    assert d_mpc == 10.0
    assert "baryonic_velocity_kms" in df.columns
    v_gas = df["gas_velocity_kms"].to_numpy(dtype=float)
    v_disk = df["disk_velocity_kms"].to_numpy(dtype=float)
    v_bulge = df["bulge_velocity_kms"].to_numpy(dtype=float)
    expected_bar = np.sqrt(v_gas**2 + v_disk**2 + v_bulge**2)
    assert np.allclose(df["baryonic_velocity_kms"].to_numpy(dtype=float), expected_bar)


def test_pa_supplement_loader_returns_normalised_keys(tmp_path: Path):
    """Supplement keys must be normalised so they join the morphology lookup directly."""
    payload = {
        "description": "synthetic",
        "positions": {
            "PA 001": {
                "position_angle_deg": 41.0,
                "position_angle_uncertainty_deg": 5.0,
                "source": "HyperLEDA",
            },
            "PA 002": {
                "position_angle_deg": 230.0,
                "position_angle_uncertainty_deg": 5.0,
                "source": "HyperLEDA",
            },
        },
    }
    path = tmp_path / "disk_position_angle_overlay.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    out = _load_pa_supplement(path)
    assert "PA1" in out
    assert "PA2" in out
    assert out["PA1"]["position_angle_deg"] == 41.0
    assert out["PA1"]["source"] == "HyperLEDA"
    assert out["PA2"]["position_angle_uncertainty_deg"] == 5.0


def test_pa_supplement_missing_file_returns_empty():
    """Absent supplement file is not an error; engine falls back to axisymmetric."""
    out = _load_pa_supplement(Path("/no/such/path/disk_position_angle_overlay.json"))
    assert out == {}


def test_morphology_provider_injects_pa_from_supplement(tmp_path: Path, sample_morphology_json: Path):
    """When the supplement sits next to the morphology JSON, PA must surface in the lookup."""
    # Copy the sample morphology JSON to a tmp dir and add a PA supplement next to it.
    target = tmp_path / "sample_morphology.json"
    target.write_text(Path(sample_morphology_json).read_text(encoding="utf-8"), encoding="utf-8")
    (tmp_path / "disk_position_angle_overlay.json").write_text(
        json.dumps(
            {
                "positions": {
                    "SYN001": {
                        "position_angle_deg": 73.5,
                        "position_angle_uncertainty_deg": 4.0,
                        "source": "HyperLEDA",
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    prov = SparcMorphologyProvider(target)
    lookup = prov.load_morphology_lookup()
    key = prov.normalize_galaxy_name("SYN001")
    assert key in lookup
    assert lookup[key]["position_angle_deg"] == 73.5
    assert lookup[key]["position_angle_uncertainty_deg"] == 4.0


def test_morphology_provider_handles_absent_supplement(sample_morphology_json: Path):
    """No supplement file → ``position_angle_deg`` defaults to ``None``."""
    prov = SparcMorphologyProvider(sample_morphology_json)
    lookup = prov.load_morphology_lookup()
    key = prov.normalize_galaxy_name("SYN001")
    assert key in lookup
    assert lookup[key].get("position_angle_deg") is None
