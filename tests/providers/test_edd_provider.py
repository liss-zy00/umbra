from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from src.config.thresholds import luminosity_distance_mpc_from_distance_modulus
from src.providers.edd_provider import EddProvider, build_edd_provenance_log_dict, build_sparc_name_to_pgc_lookup
from src.providers.sparc import SparcMorphologyProvider


def test_edd_hierarchy_prefers_trgb_over_cf4() -> None:
    trgb = pd.DataFrame(
        [
            {
                "pgc": 1,
                "name_cmd": "x",
                "m_trgb": 0.0,
                "dm_tip": 30.0,
                "edm_lo": 0.1,
                "edm_hi": 0.1,
            }
        ]
    )
    cf4 = pd.DataFrame(
        [
            {
                "pgc": 1,
                "pgc1": 1,
                "vcmb": 0,
                "DM": 40.0,
                "eDM": 0.1,
                "DMtrgb": np.nan,
                "eDMt": np.nan,
                "DMcep": np.nan,
                "eDMcep": np.nan,
            }
        ]
    )
    cf3 = pd.DataFrame([{"pgc": 1, "Dist": 100.0, "DM": 35.0, "eDM": 0.2, "Ty": 0.0, "Name": "X"}])
    prov = EddProvider(trgb, cf4, cf3)
    hit = prov.lookup_distance_for_pgc(1)
    assert hit is not None
    assert hit.tier == 1
    assert hit.source_table == "edd_cmd_trgb"
    assert hit.distance_mpc == pytest.approx(luminosity_distance_mpc_from_distance_modulus(30.0), rel=1e-9)


def test_edd_cf4_dmcep_used_when_no_trgb() -> None:
    trgb = pd.DataFrame(
        [
            {
                "pgc": 999,
                "name_cmd": "x",
                "m_trgb": 0.0,
                "dm_tip": 40.0,
                "edm_lo": 0.1,
                "edm_hi": 0.1,
            }
        ]
    )
    cf4 = pd.DataFrame(
        [
            {
                "pgc": 7,
                "pgc1": 7,
                "vcmb": 0,
                "DM": 40.0,
                "eDM": 0.5,
                "DMtrgb": np.nan,
                "eDMt": np.nan,
                "DMcep": 29.0,
                "eDMcep": 0.05,
            }
        ]
    )
    cf3 = pd.DataFrame([{"pgc": 7, "Dist": 2.0, "DM": 31.0, "eDM": 0.1, "Ty": 0.0, "Name": "G7"}])
    prov = EddProvider(trgb, cf4, cf3)
    hit = prov.lookup_distance_for_pgc(7)
    assert hit is not None
    assert hit.tier == 2
    assert "DMcep" in hit.source_table


def test_build_edd_provenance_log_dict_reads_pgc_and_metadata() -> None:
    from types import SimpleNamespace

    galaxy = SimpleNamespace(
        metadata={
            "prov_high_precision_dist": True,
            "edd_distance_tier": 2,
            "edd_distance_source_table": "edd_cf4_DM",
        },
        distance_mpc=10.25,
        distance_uncertainty_mpc=0.4,
    )
    morph_key = "ngc0001"
    deps = SimpleNamespace(edd=object(), pgc_by_sparc_name={morph_key: 42})
    d = build_edd_provenance_log_dict(deps, morphology_key=morph_key, galaxy_name="NGC 1", galaxy=galaxy)
    assert d["edd_tables_loaded"] is True
    assert d["pgc"] == 42
    assert d["applied"] is True
    assert d["tier"] == 2
    assert d["source_table"] == "edd_cf4_DM"
    assert d["distance_mpc"] == 10.25
    assert d["distance_uncertainty_mpc"] == 0.4


def test_pgc_lookup_json_overrides_cosmicflows(tmp_path, sample_morphology_json) -> None:
    cf3 = tmp_path / "cf3.csv"
    cf3.write_text(
        "\n".join(
            [
                "HYPERLEDA,CF3",
                "pgc,Dist,DM,eDM,Ty,Name",
                "%7s,%6.2f,%5.2f,%4.2f,%4.1f,%11s",
                "---,Mpc,mag,,,",
                '"note"',
                "99,1.0,30.0,0.1,5.0,UGC00099",
            ]
        ),
        encoding="utf-8",
    )
    override = tmp_path / "pgc.json"
    override.write_text(
        json.dumps({"positions": {"UGC99": {"pgc": 42}}}),
        encoding="utf-8",
    )
    morph = SparcMorphologyProvider(sample_morphology_json)
    m = build_sparc_name_to_pgc_lookup(override, cf3, morph)
    assert m[morph.normalize_galaxy_name("UGC99")] == 42
