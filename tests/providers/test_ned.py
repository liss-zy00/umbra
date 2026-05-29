import json

import pytest

from src.providers.ned import NedProvider


# ---------------------------------------------------------------------------
# Test doubles
# ---------------------------------------------------------------------------


class _Table:
    """Minimal duck-typed stand-in for ``astropy.table.Table`` rows used by NedProvider."""

    def __init__(self, columns: dict[str, list]) -> None:
        self._columns = columns
        self.colnames = list(columns.keys())

    def __len__(self):
        if not self._columns:
            return 0
        return len(next(iter(self._columns.values())))

    def __getitem__(self, key):
        return self._columns[key]


def _empty_table() -> _Table:
    return _Table({})


def _stub_query_region_always_empty(monkeypatch) -> None:
    """Neutralise the cone-search path in unit tests that do not exercise it."""
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    monkeypatch.setattr(
        ipac_ned.Ned,
        "query_region",
        lambda *args, **kwargs: _empty_table(),
    )


def _stub_simbad_always_empty(monkeypatch) -> None:
    simbad_mod = pytest.importorskip("astroquery.simbad")
    monkeypatch.setattr(simbad_mod.Simbad, "query_object", lambda name: _empty_table())


# ---------------------------------------------------------------------------
# Basic contract
# ---------------------------------------------------------------------------


def test_ned_graceful_when_empty_name():
    ned = NedProvider()
    meta = ned.fetch_metadata("   ")
    assert meta["right_ascension_degrees"] is None
    assert meta["raw_status"] == "empty_name"


def test_ned_mocked_success(monkeypatch):
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    monkeypatch.setattr(
        ipac_ned.Ned,
        "query_object",
        lambda name: _Table({"RA": [10.0], "DEC": [20.0], "Redshift": [0.01]}),
    )
    _stub_query_region_always_empty(monkeypatch)
    _stub_simbad_always_empty(monkeypatch)

    ned = NedProvider()
    meta = ned.fetch_metadata("dummy")
    assert meta["raw_status"] == "ok"
    assert meta["right_ascension_degrees"] == 10.0
    assert meta["declination_degrees"] == 20.0


def test_ned_handles_empty_table(monkeypatch):
    """Empty NED response must fall through to SIMBAD; unresolved when both return empty."""
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    monkeypatch.setattr(ipac_ned.Ned, "query_object", lambda name: _empty_table())
    _stub_query_region_always_empty(monkeypatch)
    _stub_simbad_always_empty(monkeypatch)

    ned = NedProvider()
    meta = ned.fetch_metadata("ghost-galaxy")
    assert meta["raw_status"] == "unresolved"
    assert meta["right_ascension_degrees"] is None


def test_ned_handles_object_type_and_invalid_redshift(monkeypatch):
    """Cover the object-type fallback column and the redshift parser graceful-degrade path."""
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    monkeypatch.setattr(
        ipac_ned.Ned,
        "query_object",
        lambda name: _Table({"RA": [3.5], "DEC": [4.5], "z": ["--"], "Type": ["G"]}),
    )
    _stub_query_region_always_empty(monkeypatch)
    _stub_simbad_always_empty(monkeypatch)

    ned = NedProvider()
    meta = ned.fetch_metadata("M-something")
    assert meta["raw_status"] == "ok"
    assert meta["redshift"] is None
    assert meta["object_type"] == "G"


def test_ned_swallows_unexpected_exceptions(monkeypatch):
    """Generic exceptions per attempt are debug-logged and resolution proceeds to ``unresolved``."""
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")

    def explode(name):
        raise RuntimeError("boom")

    monkeypatch.setattr(ipac_ned.Ned, "query_object", explode)
    _stub_query_region_always_empty(monkeypatch)
    _stub_simbad_always_empty(monkeypatch)

    ned = NedProvider()
    meta = ned.fetch_metadata("X")
    assert meta["raw_status"] == "unresolved"
    assert meta["right_ascension_degrees"] is None


# ---------------------------------------------------------------------------
# Alias map (exact and legacy flat schema)
# ---------------------------------------------------------------------------


def test_ned_uses_legacy_flat_alias_map(tmp_path, monkeypatch):
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    alias_path = tmp_path / "aliases.json"
    alias_path.write_text(json.dumps({"GX1": "NED_CANONICAL"}), encoding="utf-8")

    captured: list[str] = []

    def capture_query(name):
        captured.append(str(name))
        return _Table({"RA": [1.0], "DEC": [2.0]})

    monkeypatch.setattr(ipac_ned.Ned, "query_object", capture_query)
    _stub_query_region_always_empty(monkeypatch)
    _stub_simbad_always_empty(monkeypatch)

    ned = NedProvider(alias_map_path=alias_path)
    meta = ned.fetch_metadata("GX1")
    assert captured == ["NED_CANONICAL"]
    assert meta["raw_status"] == "ok_via_alias_map"
    assert meta["right_ascension_degrees"] == 1.0


def test_ned_uses_sectioned_alias_map(tmp_path, monkeypatch):
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    alias_path = tmp_path / "aliases.json"
    alias_path.write_text(
        json.dumps({"alias_map": {"GX2": "CANONICAL_2"}}),
        encoding="utf-8",
    )

    captured: list[str] = []

    def capture_query(name):
        captured.append(str(name))
        return _Table({"RA": [5.0], "DEC": [6.0]})

    monkeypatch.setattr(ipac_ned.Ned, "query_object", capture_query)
    _stub_query_region_always_empty(monkeypatch)
    _stub_simbad_always_empty(monkeypatch)

    ned = NedProvider(alias_map_path=alias_path)
    meta = ned.fetch_metadata("GX2")
    assert captured == ["CANONICAL_2"]
    assert meta["raw_status"] == "ok_via_alias_map"


# ---------------------------------------------------------------------------
# Pattern-rule + prefix permutation
# ---------------------------------------------------------------------------


def test_ned_pattern_rule_resolves_when_bare_name_fails(tmp_path, monkeypatch):
    """``F561-1`` matches the LSB pattern; the templated query is tried before any prefix."""
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    alias_path = tmp_path / "aliases.json"
    alias_path.write_text(
        json.dumps(
            {
                "alias_map": {},
                "pattern_rules": [
                    {"regex": r"^F\d{3}-.+$", "ned_query_template": "LSB {name}"}
                ],
                "ned_prefix_attempts": ["LSB ", "ESO "],
            }
        ),
        encoding="utf-8",
    )

    captured: list[str] = []

    def fake_query(name):
        captured.append(str(name))
        if name == "LSB F561-1":
            return _Table({"RA": [123.4], "DEC": [-12.3]})
        return _empty_table()

    monkeypatch.setattr(ipac_ned.Ned, "query_object", fake_query)
    _stub_query_region_always_empty(monkeypatch)
    _stub_simbad_always_empty(monkeypatch)

    ned = NedProvider(alias_map_path=alias_path)
    meta = ned.fetch_metadata("F561-1")
    assert captured[0] == "LSB F561-1"
    assert meta["raw_status"] == "ok_via_pattern_rule"
    assert meta["right_ascension_degrees"] == 123.4


def test_ned_prefix_attempt_resolves_when_pattern_does_not_match(tmp_path, monkeypatch):
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    alias_path = tmp_path / "aliases.json"
    alias_path.write_text(
        json.dumps(
            {
                "alias_map": {},
                "pattern_rules": [],
                "ned_prefix_attempts": ["LSB ", "ESO "],
            }
        ),
        encoding="utf-8",
    )

    def fake_query(name):
        if name == "ESO MyDwarf":
            return _Table({"RA": [42.0], "DEC": [-17.0]})
        return _empty_table()

    monkeypatch.setattr(ipac_ned.Ned, "query_object", fake_query)
    _stub_query_region_always_empty(monkeypatch)
    _stub_simbad_always_empty(monkeypatch)

    ned = NedProvider(alias_map_path=alias_path)
    meta = ned.fetch_metadata("MyDwarf")
    assert meta["raw_status"] == "ok_via_prefix_attempt"
    assert meta["right_ascension_degrees"] == 42.0


# ---------------------------------------------------------------------------
# SIMBAD fallback
# ---------------------------------------------------------------------------


def test_ned_simbad_fallback_when_every_ned_candidate_empty(tmp_path, monkeypatch):
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    simbad_mod = pytest.importorskip("astroquery.simbad")

    alias_path = tmp_path / "aliases.json"
    alias_path.write_text(
        json.dumps({"alias_map": {}, "pattern_rules": [], "ned_prefix_attempts": ["LSB "]}),
        encoding="utf-8",
    )

    monkeypatch.setattr(ipac_ned.Ned, "query_object", lambda name: _empty_table())

    def simbad_q(name):
        if name == "LSB MyGal":
            return _Table({"ra": [100.0], "dec": [10.0], "main_id": ["LSB X"]})
        return _empty_table()

    monkeypatch.setattr(simbad_mod.Simbad, "query_object", simbad_q)

    # SIMBAD hit will trigger cone search. Return empty cone so we exit via
    # ``ok_via_simbad`` and not ``ok_via_ned_cone``.
    monkeypatch.setattr(ipac_ned.Ned, "query_region", lambda *a, **k: _empty_table())

    ned = NedProvider(alias_map_path=alias_path)
    meta = ned.fetch_metadata("MyGal")
    assert meta["raw_status"] == "ok_via_simbad"
    assert meta["right_ascension_degrees"] == 100.0


# ---------------------------------------------------------------------------
# NED cone search
# ---------------------------------------------------------------------------


def test_ned_cone_search_seeded_by_simbad(tmp_path, monkeypatch):
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    simbad_mod = pytest.importorskip("astroquery.simbad")

    alias_path = tmp_path / "aliases.json"
    alias_path.write_text(
        json.dumps({"alias_map": {}, "pattern_rules": [], "ned_prefix_attempts": []}),
        encoding="utf-8",
    )

    monkeypatch.setattr(ipac_ned.Ned, "query_object", lambda name: _empty_table())
    monkeypatch.setattr(
        simbad_mod.Simbad,
        "query_object",
        lambda name: _Table({"ra": [101.0], "dec": [-15.0], "main_id": ["LSB Y"]}),
    )

    captured_regions: list = []

    def fake_region(coord, radius):
        captured_regions.append((coord, radius))
        return _Table({"RA": [101.001], "DEC": [-15.002], "Redshift": [0.005]})

    monkeypatch.setattr(ipac_ned.Ned, "query_region", fake_region)

    ned = NedProvider(alias_map_path=alias_path)
    meta = ned.fetch_metadata("HiddenGal")
    assert meta["raw_status"] == "ok_via_ned_cone"
    assert meta["right_ascension_degrees"] == 101.001
    assert meta["redshift"] == 0.005
    assert len(captured_regions) == 1


def test_ned_cone_search_seeded_by_icrs_supplement(tmp_path, monkeypatch):
    """When SIMBAD also returns nothing, an explicit ``icrs_deg`` seed must still trigger the cone."""
    ipac_ned = pytest.importorskip("astroquery.ipac.ned")
    simbad_mod = pytest.importorskip("astroquery.simbad")

    alias_path = tmp_path / "aliases.json"
    alias_path.write_text(
        json.dumps({"alias_map": {}, "pattern_rules": [], "ned_prefix_attempts": []}),
        encoding="utf-8",
    )

    monkeypatch.setattr(ipac_ned.Ned, "query_object", lambda name: _empty_table())
    monkeypatch.setattr(simbad_mod.Simbad, "query_object", lambda name: _empty_table())
    monkeypatch.setattr(
        ipac_ned.Ned,
        "query_region",
        lambda *a, **k: _Table({"RA": [50.5], "DEC": [-2.5]}),
    )

    ned = NedProvider(alias_map_path=alias_path)
    meta = ned.fetch_metadata("OrphanGal", icrs_deg=(50.5, -2.5))
    assert meta["raw_status"] == "ok_via_ned_cone"
    assert meta["right_ascension_degrees"] == 50.5


# ---------------------------------------------------------------------------
# Inventory column wiring
# ---------------------------------------------------------------------------


def test_save_full_catalog_inventory_contains_ned_status_column(tmp_path):
    """The CSV inventory must expose ``ned_status`` so downstream consumers can audit resolution outcomes."""
    import pandas as pd

    from src.pipeline.constants import DEFAULT_CONTROL_GALAXIES
    from src.reporting.persistence.exporter import ResultPersister

    g_ok = DEFAULT_CONTROL_GALAXIES[0]
    persister = ResultPersister(tmp_path)
    output_path = persister.save_full_catalog_inventory(
        [
            {
                "galaxy_name": g_ok,
                "right_ascension_degrees": 1.0,
                "declination_degrees": 2.0,
                "galactic_latitude_degrees": 60.0,
                "validation_is_reliable": True,
                "ned_status": "OK_NED_SPARC_LABEL_DIRECT",
                "mean_acceleration_gap_edge": -1e-11,
                "morphology_class": "Dwarf Irregulars (Im)",
                "hubble_type": 10.0,
                "audit_flags": [],
                "reliability_failure_codes": [],
                "note": "Reliable for synthesis (no blocking audit gates; audit_flags empty).",
            },
            {
                "galaxy_name": "Phantom",
                "right_ascension_degrees": None,
                "declination_degrees": None,
                "galactic_latitude_degrees": None,
                "validation_is_reliable": False,
                "ned_status": "NOT_FOUND_NED",
                "mean_acceleration_gap_edge": None,
                "morphology_class": "Unclassified",
                "hubble_type": None,
                "audit_flags": ["fail_ned_resolution"],
                "reliability_failure_codes": ["fail_ned_resolution"],
                "note": "NED unresolved.",
            },
        ]
    )
    df = pd.read_csv(output_path)
    assert "ned_status" in df.columns
    assert "audit_flags" in df.columns
    assert "note" in df.columns
    assert df.loc[df["galaxy_name"] == g_ok, "ned_status"].iloc[0] == "OK_NED_SPARC_LABEL_DIRECT"
    assert df.loc[df["galaxy_name"] == "Phantom", "ned_status"].iloc[0] == "NOT_FOUND_NED"
    assert df.loc[df["galaxy_name"] == "Phantom", "audit_flags"].iloc[0] == "fail_ned_resolution"
    assert "Reliable" in str(df.loc[df["galaxy_name"] == g_ok, "note"].iloc[0])
