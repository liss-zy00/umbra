import pandas as pd

from src.models.galaxy import Galaxy
from src.reporting.persistence.exporter import ResultPersister


def test_save_evidence_and_validation_partitions(tmp_path):
    rp = ResultPersister(tmp_path)
    ev = rp.save_scientific_evidence_report([])
    assert ev.exists()
    inv = [
        {
            "galaxy_name": "G_OK",
            "right_ascension_degrees": 1.0,
            "declination_degrees": 2.0,
            "galactic_latitude_degrees": 40.0,
            "validation_is_reliable": True,
            "ned_status": "OK",
            "mean_acceleration_gap_edge": 1.0e-11,
            "morphology_class": "Im",
            "hubble_type": 10,
            "audit_flags": [],
            "reliability_failure_codes": [],
        },
        {
            "galaxy_name": "G_BAD",
            "right_ascension_degrees": 3.0,
            "declination_degrees": 4.0,
            "galactic_latitude_degrees": 41.0,
            "validation_is_reliable": False,
            "ned_status": "OK",
            "mean_acceleration_gap_edge": 1.0e-11,
            "morphology_class": "Im",
            "hubble_type": 10,
            "audit_flags": ["fail_lat"],
            "reliability_failure_codes": ["fail_lat"],
        },
    ]
    rel_path, bad_path = rp.save_validation_reliability_partitions(inv)
    assert rel_path.exists() and bad_path.exists()
    assert "G_OK" in rel_path.read_text(encoding="utf-8")
    assert "G_BAD" in bad_path.read_text(encoding="utf-8")


def test_save_galaxy_writes_csv_and_summary(tmp_path):
    rp = ResultPersister(tmp_path)
    df = pd.DataFrame(
        {
            "radius_kpc": [1.0, 2.0],
            "observed_velocity_kms": [50.0, 60.0],
        }
    )
    galaxy = Galaxy(name="TST", rotation_dataframe=df)
    summary = {"validation_is_reliable": True, "audit_flags": []}
    galaxy_dir, csv_path, json_path = rp.save_galaxy(galaxy, summary)
    assert galaxy_dir.is_dir()
    assert csv_path.exists()
    assert json_path.exists()


def test_save_full_catalog_partitions_and_evidence_report(tmp_path):
    rp = ResultPersister(tmp_path)
    rows = [
        {
            "galaxy_name": "G1",
            "right_ascension_degrees": 10.0,
            "declination_degrees": -5.0,
            "galactic_latitude_degrees": 40.0,
            "validation_is_reliable": True,
            "ned_status": "OK",
            "mean_acceleration_gap_edge": 1.0e-11,
            "morphology_class": "Im",
            "hubble_type": 10,
            "audit_flags": [],
            "reliability_failure_codes": [],
        }
    ]
    catalog_path = rp.save_full_catalog_inventory(rows)
    text = catalog_path.read_text(encoding="utf-8")
    assert "G1" in text and "morphology_class" in text and "validation_is_reliable" in text

    rel, nrel = rp.save_validation_reliability_partitions(rows)
    assert "G1" in rel.read_text(encoding="utf-8")
    ev_path = rp.save_scientific_evidence_report([])
    assert ev_path.exists()
