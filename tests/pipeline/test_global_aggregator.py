from pathlib import Path

import numpy as np
import pytest

from src.config.thresholds import wim_density_match_score
from src.pipeline.aggregation.global_aggregator import build_direct_evidence_rows
from src.providers.sparc import SparcMorphologyProvider


def test_direct_evidence_merges_discovery_columns(sample_morphology_json):
    from src.reporting.persistence.exporter import (
        DiscoveryEliteBandInfo,
        ResultPersister,
    )

    prov = SparcMorphologyProvider(sample_morphology_json)
    lookup = prov.load_morphology_lookup()
    supreme = [
        {
            "galaxy_name": "SYN001",
            "theoretical_gap_coverage_ratio": 1.0,
            "morphology_class": "Dwarf Irregulars (Im)",
            # Best-scenario grid scale is one of the engine knobs; the aggregator
            # recomputes ``resolution_correction_factor_grid`` against this number.
            "scenario_gradient_scale_kpc": 0.1,
        },
    ]
    disc = {
        "SYN001": {
            "required_density_for_match": 1.1e-27,
            "convergence_ratio_to_wim": 1.1,
            "delta_r_m": 1.0e20,
            "physical_pixel_scale_kpc": 0.05,
            "resolution_correction_factor": 70.0,
            "beam_size_kpc": 3.5,
            "distance_mpc": 4.0,
            "is_reliable": True,
        }
    }
    # Dynamic band only: ``elite_band_info`` must be supplied for ``discovery_is_elite``.
    band_info = DiscoveryEliteBandInfo(
        median_required_density_kg_m3=1.0e-27,
        band_min_kg_m3=0.7e-27,
        band_max_kg_m3=1.3e-27,
        n_reliable=1,
        n_elite=1,
        relative_min=ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MIN,
        relative_max=ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MAX,
    )
    rows = build_direct_evidence_rows(
        supreme,
        lookup,
        prov,
        discovery_lookup=disc,
        elite_band_info=band_info,
    )
    assert len(rows) == 1
    r0 = rows[0]
    assert r0["required_density_for_match"] == 1.1e-27
    assert r0["convergence_ratio_to_wim"] == 1.1
    assert r0["delta_r_m"] == 1.0e20
    assert r0["physical_pixel_scale_kpc"] == 0.05
    assert r0["beam_size_kpc"] == 3.5
    # ``resolution_correction_factor_discovery`` is the value forwarded from the
    # Discovery row (computed against the physical pixel separation).
    assert r0["resolution_correction_factor_discovery"] == 70.0
    # ``resolution_correction_factor_grid`` is recomputed against the best grid scale
    # (0.1 kpc here); for D=4 Mpc / θ=5'  / L_grad=0.1 kpc this should be ≈ 58.18.
    assert r0["resolution_correction_factor_grid"] == pytest.approx(58.178, rel=1e-3)
    assert r0["resolution_correction_factor"] == r0["resolution_correction_factor_grid"]
    assert r0["distance_mpc"] == 4.0
    assert r0["discovery_is_reliable"] is True
    # ρ_req = 1.1e-27 lies inside ``[0.7e-27, 1.3e-27]`` → elite.
    assert r0["discovery_is_elite"] is True
    assert r0["discovery_elite_band_source"] == "dynamic_rho_band"
    assert r0["validation_status"] == "STANDARD"
    assert np.isfinite(float(r0["wim_density_match_score"]))


def test_direct_evidence_supreme_validation_status(sample_morphology_json):
    from src.reporting.persistence.exporter import (
        DiscoveryEliteBandInfo,
        ResultPersister,
    )

    prov = SparcMorphologyProvider(sample_morphology_json)
    lookup = prov.load_morphology_lookup()
    supreme = [
        {
            "galaxy_name": "SYN001",
            "theoretical_gap_coverage_ratio": 1.0,
            "morphology_class": "Dwarf Irregulars (Im)",
            "scenario_gradient_scale_kpc": 0.1,
            "prov_high_precision_dist": True,
            "position_angle_deg": 45.0,
            "warn_no_pa": False,
            "warn_btfr_outlier": False,
            "audit_flags": [],
            "calibrated_convergence_ratio": 1.0,
            "r_max_kpc": 50.0,
        },
    ]
    disc = {
        "SYN001": {
            "required_density_for_match": 1.0e-27,
            "convergence_ratio_to_wim": 1.0,
            "delta_r_m": 1.0e20,
            "physical_pixel_scale_kpc": 0.05,
            "resolution_correction_factor": 1.0,
            "beam_size_kpc": 3.5,
            "distance_mpc": 4.0,
            "is_reliable": True,
        }
    }
    band_info = DiscoveryEliteBandInfo(
        median_required_density_kg_m3=1.0e-27,
        band_min_kg_m3=0.7e-27,
        band_max_kg_m3=1.3e-27,
        n_reliable=1,
        n_elite=1,
        relative_min=ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MIN,
        relative_max=ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MAX,
    )
    rows = build_direct_evidence_rows(
        supreme,
        lookup,
        prov,
        discovery_lookup=disc,
        elite_band_info=band_info,
    )

    k_grid = float(rows[0]["resolution_correction_factor"])
    assert rows[0]["validation_status"] == "SUPREME"
    assert float(rows[0]["wim_density_match_score"]) == pytest.approx(
        wim_density_match_score(1.0e-27, k_grid), rel=1e-9
    )


def test_direct_evidence_data_limited_when_beam_exceeds_diameter(sample_morphology_json):
    """Planck beam larger than 2×R_max ⇒ ``DATA_LIMITED`` (no ``SUPREME`` label)."""
    from src.reporting.persistence.exporter import (
        DiscoveryEliteBandInfo,
        ResultPersister,
    )

    prov = SparcMorphologyProvider(sample_morphology_json)
    lookup = prov.load_morphology_lookup()
    supreme = [
        {
            "galaxy_name": "SYN001",
            "theoretical_gap_coverage_ratio": 1.0,
            "morphology_class": "Dwarf Irregulars (Im)",
            "scenario_gradient_scale_kpc": 0.1,
            "prov_high_precision_dist": True,
            "position_angle_deg": 45.0,
            "warn_no_pa": False,
            "warn_btfr_outlier": False,
            "audit_flags": [],
            "calibrated_convergence_ratio": 1.0,
            "r_max_kpc": 1.0,
        },
    ]
    disc = {
        "SYN001": {
            "required_density_for_match": 1.0e-27,
            "convergence_ratio_to_wim": 1.0,
            "delta_r_m": 1.0e20,
            "physical_pixel_scale_kpc": 0.05,
            "resolution_correction_factor": 1.0,
            "beam_size_kpc": 5.0,
            "distance_mpc": 4.0,
            "is_reliable": True,
        }
    }
    band_info = DiscoveryEliteBandInfo(
        median_required_density_kg_m3=1.0e-27,
        band_min_kg_m3=0.7e-27,
        band_max_kg_m3=1.3e-27,
        n_reliable=1,
        n_elite=1,
        relative_min=ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MIN,
        relative_max=ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MAX,
    )
    rows = build_direct_evidence_rows(
        supreme,
        lookup,
        prov,
        discovery_lookup=disc,
        elite_band_info=band_info,
    )
    assert rows[0]["validation_status"] == "DATA_LIMITED"


def test_direct_evidence_standard_when_calibrated_out_of_band(sample_morphology_json):
    from src.reporting.persistence.exporter import (
        DiscoveryEliteBandInfo,
        ResultPersister,
    )

    prov = SparcMorphologyProvider(sample_morphology_json)
    lookup = prov.load_morphology_lookup()
    supreme = [
        {
            "galaxy_name": "SYN001",
            "theoretical_gap_coverage_ratio": 1.0,
            "morphology_class": "Dwarf Irregulars (Im)",
            "scenario_gradient_scale_kpc": 0.1,
            "prov_high_precision_dist": True,
            "position_angle_deg": 45.0,
            "warn_no_pa": False,
            "warn_btfr_outlier": False,
            "audit_flags": [],
            "calibrated_convergence_ratio": 100.0,
            "r_max_kpc": 50.0,
        },
    ]
    disc = {
        "SYN001": {
            "required_density_for_match": 1.0e-27,
            "convergence_ratio_to_wim": 1.0,
            "delta_r_m": 1.0e20,
            "physical_pixel_scale_kpc": 0.05,
            "resolution_correction_factor": 1.0,
            "beam_size_kpc": 3.5,
            "distance_mpc": 4.0,
            "is_reliable": True,
        }
    }
    band_info = DiscoveryEliteBandInfo(
        median_required_density_kg_m3=1.0e-27,
        band_min_kg_m3=0.7e-27,
        band_max_kg_m3=1.3e-27,
        n_reliable=1,
        n_elite=1,
        relative_min=ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MIN,
        relative_max=ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MAX,
    )
    rows = build_direct_evidence_rows(
        supreme,
        lookup,
        prov,
        discovery_lookup=disc,
        elite_band_info=band_info,
    )
    assert rows[0]["validation_status"] == "STANDARD"


def test_direct_evidence_no_elite_when_dynamic_band_missing(sample_morphology_json):
    """Without ``elite_band_info`` there is no cohort band: ``discovery_is_elite`` is always False."""
    prov = SparcMorphologyProvider(sample_morphology_json)
    lookup = prov.load_morphology_lookup()
    supreme = [
        {
            "galaxy_name": "SYN001",
            "theoretical_gap_coverage_ratio": 1.0,
            "morphology_class": "Dwarf Irregulars (Im)",
        },
    ]
    disc = {
        "SYN001": {
            "required_density_for_match": 1.1e-27,
            "convergence_ratio_to_wim": 1.1,
            "is_reliable": True,
        }
    }
    rows = build_direct_evidence_rows(supreme, lookup, prov, discovery_lookup=disc)
    assert rows[0]["discovery_is_elite"] is False
    assert rows[0]["discovery_elite_band_source"] == "none"


def test_direct_evidence_filters(sample_morphology_json):
    prov = SparcMorphologyProvider(sample_morphology_json)
    lookup = prov.load_morphology_lookup()
    supreme = [
        {
            "galaxy_name": "SYN001",
            "theoretical_gap_coverage_ratio": 1.0,
            "morphology_class": "Dwarf Irregulars (Im)",
        },
        {
            "galaxy_name": "SYN002",
            "theoretical_gap_coverage_ratio": 100.0,
            "morphology_class": "Early-type Spirals (Sa, Sab)",
        },
    ]
    rows = build_direct_evidence_rows(supreme, lookup, prov)
    assert len(rows) >= 1
    assert rows[0]["galaxy_name"] == "SYN001"


def test_median_dm_equivalence_grid_basic():
    """Median grid: only reliable rows contribute and the value is in percent."""
    from src.pipeline.aggregation.global_aggregator import _build_median_dm_equivalence_grid

    scenario_rows = {
        (0.05, 0.001): [
            {"is_reliable": True, "theoretical_gap_coverage_ratio": 0.5},
            {"is_reliable": True, "theoretical_gap_coverage_ratio": 1.5},
            {"is_reliable": False, "theoretical_gap_coverage_ratio": 5.0},
        ],
        (0.05, 0.01): [],
        (0.10, 0.001): [
            {"is_reliable": True, "theoretical_gap_coverage_ratio": 1.0},
        ],
        (0.10, 0.01): [{"is_reliable": False, "theoretical_gap_coverage_ratio": 0.9}],
    }
    grid = _build_median_dm_equivalence_grid(scenario_rows, (0.05, 0.10), (0.001, 0.01))
    assert grid.shape == (2, 2)
    assert abs(grid[0, 0] - 100.0) < 1e-6  # median(0.5, 1.5) * 100
    assert np.isnan(grid[0, 1])
    assert abs(grid[1, 0] - 100.0) < 1e-6
    assert np.isnan(grid[1, 1])  # only reliable rows count
