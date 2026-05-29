import numpy as np
import pandas as pd

from src.config.thresholds import AUDIT_CODES, CALIBRATION_STATUS_OPTIMAL, CALIBRATION_STATUS_SUBOPTIMAL
from src.pipeline.per_galaxy.metrics_builder import (
    ALL_AUDIT_FLAG_CODES,
    AUDIT_FLAG_FAIL_CALIBRATION,
    AUDIT_FLAG_FAIL_DENSITY_CLAMPED,
    AUDIT_FLAG_FAIL_INCLINATION,
    AUDIT_FLAG_FAIL_LATITUDE,
    AUDIT_FLAG_FAIL_RELIABLE_FRACTION,
    AUDIT_FLAG_FAIL_SKIPPED_PHYSICS,
    AUDIT_FLAG_FAIL_STELLAR,
    AUDIT_FLAG_FAIL_NED_RESOLUTION,
    AUDIT_FLAG_WARN_B_BAND,
    AUDIT_FLAG_WARN_BTFR_OUTLIER,
    AUDIT_FLAG_WARN_NO_POSITION_ANGLE,
    AUDIT_FLAG_WARN_PURE_PRESSURE_GEOMETRY,
    AUDIT_FLAG_WARN_STOKES_NEIGHBOR_IMPUTED,
    AUDIT_FLAG_WARN_UNKNOWN_MORPHOLOGY,
    aggregate_galaxy_audit_flags,
    blocking_flags,
    metrics_for_skipped_physics,
    metrics_for_successful_physics,
)


def test_successful_metrics_structure():
    df = pd.DataFrame(
        {
            "acceleration_gap": np.array([1e-10, 2e-10], dtype=float),
            "magnetic_pressure_acceleration": np.array([5e-11, 1e-10], dtype=float),
            "reliable_measurement": np.array([True, True]),
        }
    )
    m = metrics_for_successful_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-10,
        galactic_latitude_degrees=30.0,
        morphology_class="Im",
    )
    assert "_increment_success" in m
    m.pop("_increment_success")
    assert m["galaxy_name"] == "G"
    assert m["is_reliable"] is True
    # Reliable scenarios carry an EMPTY ``audit_flags`` list — that's the source of truth.
    assert m["audit_flags"] == []
    assert m["reliability_failure_codes"] == []
    assert m["reliable_measurement_fraction"] == 1.0


def test_skipped_metrics():
    m = metrics_for_skipped_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        mean_acceleration_gap_edge=None,
        galactic_latitude_degrees=None,
        morphology_class="Im",
    )
    assert m.pop("_increment_success") is False
    assert np.isnan(m["theoretical_gap_coverage_ratio"])
    # Skipped scenarios always carry the structured ``fail_skipped_physics`` code.
    assert AUDIT_FLAG_FAIL_SKIPPED_PHYSICS in m["audit_flags"]
    # Latitude is missing too → both codes should be present.
    assert AUDIT_FLAG_FAIL_LATITUDE in m["audit_flags"]


def test_density_clamp_audit_forces_unreliable_when_floor_dominates():
    """If most shells were clamped to the gas-density floor, the scenario must be flagged unreliable.

    All other reliability gates pass here (latitude, reliable_measurement fraction), so the only
    source of failure must be the density-clamp audit.
    """
    df = pd.DataFrame(
        {
            "acceleration_gap": np.array([1e-10, 2e-10, 3e-10, 4e-10], dtype=float),
            "magnetic_pressure_acceleration": np.array([5e-11, 1e-10, 1.5e-10, 2e-10], dtype=float),
            "reliable_measurement": np.array([True, True, True, True]),
            "density_floor_clamped": np.array([True, True, True, False]),
        }
    )
    m = metrics_for_successful_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-10,
        galactic_latitude_degrees=30.0,
        morphology_class="Im",
    )
    assert m["n_points_density_clamped"] == 3
    assert m["n_points_total"] == 4
    assert m["density_clamped_fraction"] == 0.75
    assert m["density_audit_passed"] is False
    assert m["is_reliable"] is False
    # And the EXACT failure mode is recorded in the structured flag list.
    assert AUDIT_FLAG_FAIL_DENSITY_CLAMPED in m["audit_flags"]
    assert AUDIT_FLAG_FAIL_DENSITY_CLAMPED in m["reliability_failure_codes"]


def test_plasma_beta_below_unity_fraction_and_dominance_flag():
    """When β<1 on the majority of shells the dominance flag flips to True."""
    df = pd.DataFrame(
        {
            "acceleration_gap": np.array([1e-10, 2e-10, 3e-10, 4e-10], dtype=float),
            "magnetic_pressure_acceleration": np.array([5e-11, 1e-10, 1.5e-10, 2e-10], dtype=float),
            "reliable_measurement": np.array([True, True, True, True]),
            "magnetic_field_microgauss": np.array([1.0, 1.5, 2.0, 2.5], dtype=float),
            "plasma_beta_thermal": np.array([0.1, 0.5, 0.8, 1.0e-3], dtype=float),
        }
    )
    m = metrics_for_successful_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-10,
        galactic_latitude_degrees=30.0,
        morphology_class="Im",
    )
    assert m["plasma_beta_thermal_median"] is not None
    assert m["plasma_beta_below_unity_fraction"] == 1.0
    assert m["magnetic_dominance_flag"] is True
    assert m["magnetic_field_microgauss_median"] is not None
    assert m["magnetic_field_within_plausibility_band"] is True


def test_audit_flags_records_low_latitude_explicitly():
    """A scenario at |b| < 20° must still pass density / reliable_measurement gates and ONLY tag fail_lat."""
    df = pd.DataFrame(
        {
            "acceleration_gap": np.array([1e-10, 2e-10], dtype=float),
            "magnetic_pressure_acceleration": np.array([5e-11, 1e-10], dtype=float),
            "reliable_measurement": np.array([True, True]),
        }
    )
    m = metrics_for_successful_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-10,
        galactic_latitude_degrees=5.0,
        morphology_class="Im",
    )
    assert m["audit_flags"] == [AUDIT_FLAG_FAIL_LATITUDE]
    assert m["is_reliable"] is False


def test_audit_flags_records_low_reliable_fraction_explicitly():
    """A scenario where MOST shells fail the per-shell gate should tag fail_reliable_fraction."""
    df = pd.DataFrame(
        {
            "acceleration_gap": np.array([1e-10, 2e-10, 3e-10, 4e-10], dtype=float),
            "magnetic_pressure_acceleration": np.array([5e-11, 1e-10, 1.5e-10, 2e-10], dtype=float),
            "reliable_measurement": np.array([True, False, False, False]),
        }
    )
    m = metrics_for_successful_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-10,
        galactic_latitude_degrees=45.0,
        morphology_class="Im",
    )
    assert AUDIT_FLAG_FAIL_RELIABLE_FRACTION in m["audit_flags"]
    assert m["reliable_measurement_fraction"] == 0.25
    assert m["is_reliable"] is False


def test_audit_flags_stellar_dominance_independently_recorded():
    df = pd.DataFrame(
        {
            "acceleration_gap": np.array([1e-10, 2e-10, 3e-10, 4e-10], dtype=float),
            "magnetic_pressure_acceleration": np.array([5e-11, 1e-10, 1.5e-10, 2e-10], dtype=float),
            "reliable_measurement": np.array([True, True, True, True]),
            "stellar_dominated": np.array([True, True, True, False]),
        }
    )
    m = metrics_for_successful_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-10,
        galactic_latitude_degrees=45.0,
        morphology_class="Im",
    )
    assert AUDIT_FLAG_FAIL_STELLAR in m["audit_flags"]
    assert m["stellar_dominance_audit_passed"] is False


def test_plasma_beta_dominance_flag_is_false_when_thermal_swamps():
    df = pd.DataFrame(
        {
            "acceleration_gap": np.array([1e-10, 2e-10], dtype=float),
            "magnetic_pressure_acceleration": np.array([5e-11, 1e-10], dtype=float),
            "reliable_measurement": np.array([True, True]),
            "magnetic_field_microgauss": np.array([0.5, 0.5], dtype=float),
            "plasma_beta_thermal": np.array([50.0, 200.0], dtype=float),
        }
    )
    m = metrics_for_successful_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-10,
        galactic_latitude_degrees=30.0,
        morphology_class="Im",
    )
    assert m["plasma_beta_below_unity_fraction"] == 0.0
    assert m["magnetic_dominance_flag"] is False


# --------------------------------------------------------------------------------------
# Phase 5b — explicit per-failure-mode assertions on the structured audit-flag contract.
# Every gate must produce a SPECIFIC, traceable token; multiple gates may coexist; warn
# tokens must never block reliability. The galaxy-level aggregator append rules
# (inclination, calibration) are also exercised against the centralised AUDIT_CODES.
# --------------------------------------------------------------------------------------


def test_audit_flag_aliases_match_centralised_audit_codes():
    """Aliases in ``metrics_builder`` are 1:1 with :data:`AUDIT_CODES` — single source of truth."""
    assert AUDIT_FLAG_FAIL_LATITUDE == AUDIT_CODES.FAIL_LATITUDE
    assert AUDIT_FLAG_FAIL_RELIABLE_FRACTION == AUDIT_CODES.FAIL_RELIABLE_FRACTION
    assert AUDIT_FLAG_FAIL_DENSITY_CLAMPED == AUDIT_CODES.FAIL_DENSITY_CLAMPED
    assert AUDIT_FLAG_FAIL_STELLAR == AUDIT_CODES.FAIL_STELLAR
    assert AUDIT_FLAG_FAIL_SKIPPED_PHYSICS == AUDIT_CODES.FAIL_SKIPPED_PHYSICS
    assert AUDIT_FLAG_FAIL_INCLINATION == AUDIT_CODES.FAIL_INCLINATION
    assert AUDIT_FLAG_FAIL_CALIBRATION == AUDIT_CODES.FAIL_CALIBRATION
    assert AUDIT_FLAG_FAIL_NED_RESOLUTION == AUDIT_CODES.FAIL_NED_RESOLUTION
    assert AUDIT_FLAG_WARN_B_BAND == AUDIT_CODES.WARN_MAGNETIC_FIELD_BAND
    assert AUDIT_FLAG_WARN_NO_POSITION_ANGLE == AUDIT_CODES.WARN_NO_POSITION_ANGLE
    assert AUDIT_FLAG_WARN_PURE_PRESSURE_GEOMETRY == AUDIT_CODES.WARN_PURE_PRESSURE_GEOMETRY
    assert AUDIT_FLAG_WARN_BTFR_OUTLIER == AUDIT_CODES.WARN_BTFR_OUTLIER
    assert AUDIT_FLAG_WARN_UNKNOWN_MORPHOLOGY == AUDIT_CODES.WARN_UNKNOWN_MORPHOLOGY
    assert AUDIT_FLAG_WARN_STOKES_NEIGHBOR_IMPUTED == AUDIT_CODES.WARN_STOKES_NEIGHBOR_IMPUTED_FROM_CENTER
    # ``ALL_AUDIT_FLAG_CODES`` enumerates the full canonical set for downstream consumers.
    assert set(ALL_AUDIT_FLAG_CODES) == {
        AUDIT_CODES.FAIL_LATITUDE,
        AUDIT_CODES.FAIL_RELIABLE_FRACTION,
        AUDIT_CODES.FAIL_DENSITY_CLAMPED,
        AUDIT_CODES.FAIL_STELLAR,
        AUDIT_CODES.FAIL_SKIPPED_PHYSICS,
        AUDIT_CODES.FAIL_INCLINATION,
        AUDIT_CODES.FAIL_CALIBRATION,
        AUDIT_CODES.FAIL_NED_RESOLUTION,
        AUDIT_CODES.WARN_MAGNETIC_FIELD_BAND,
        AUDIT_CODES.WARN_NO_POSITION_ANGLE,
        AUDIT_CODES.WARN_PURE_PRESSURE_GEOMETRY,
        AUDIT_CODES.WARN_BTFR_OUTLIER,
        AUDIT_CODES.WARN_UNKNOWN_MORPHOLOGY,
        AUDIT_CODES.WARN_STOKES_NEIGHBOR_IMPUTED_FROM_CENTER,
    }


def test_blocking_flags_excludes_warn_pure_pressure():
    """``warn_pure_pressure`` is informational and must not block reliability."""
    flags = [AUDIT_CODES.FAIL_LATITUDE, AUDIT_CODES.WARN_PURE_PRESSURE_GEOMETRY]
    blocked = blocking_flags(flags)
    assert AUDIT_CODES.WARN_PURE_PRESSURE_GEOMETRY not in blocked
    assert AUDIT_CODES.FAIL_LATITUDE in blocked


def test_blocking_flags_excludes_warn_prefix():
    """``warn_*`` codes are informational and must never appear in ``reliability_failure_codes``."""
    flags = [
        AUDIT_FLAG_FAIL_LATITUDE,
        AUDIT_FLAG_WARN_B_BAND,
        AUDIT_FLAG_FAIL_DENSITY_CLAMPED,
    ]
    blocked = blocking_flags(flags)
    assert AUDIT_FLAG_WARN_B_BAND not in blocked
    assert AUDIT_FLAG_FAIL_LATITUDE in blocked
    assert AUDIT_FLAG_FAIL_DENSITY_CLAMPED in blocked


def test_warn_b_band_emitted_outside_plausibility_band_but_does_not_block_reliability():
    """When |B| sits outside the 0.1–100 µG band, the ``warn_b_band`` warning must surface
    but ``is_reliable`` must remain ``True`` because the band is a *display* annotation only."""
    df = pd.DataFrame(
        {
            "acceleration_gap": np.array([1e-10, 2e-10], dtype=float),
            "magnetic_pressure_acceleration": np.array([5e-11, 1e-10], dtype=float),
            "reliable_measurement": np.array([True, True]),
            "magnetic_field_microgauss": np.array([0.001, 0.002], dtype=float),
        }
    )
    m = metrics_for_successful_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-10,
        galactic_latitude_degrees=30.0,
        morphology_class="Im",
    )
    assert AUDIT_FLAG_WARN_B_BAND in m["audit_flags"]
    assert m["magnetic_field_within_plausibility_band"] is False
    # Warn-only code: must NOT block reliability and must NOT appear in ``reliability_failure_codes``.
    assert AUDIT_FLAG_WARN_B_BAND not in m["reliability_failure_codes"]
    assert m["is_reliable"] is True


def test_multiple_failure_modes_coexist_in_audit_flags():
    """Latitude + reliable-fraction + density-clamp + stellar all fail simultaneously."""
    df = pd.DataFrame(
        {
            "acceleration_gap": np.array([1e-10, 2e-10, 3e-10, 4e-10], dtype=float),
            "magnetic_pressure_acceleration": np.array([5e-11, 1e-10, 1.5e-10, 2e-10], dtype=float),
            "reliable_measurement": np.array([True, False, False, False]),
            "density_floor_clamped": np.array([True, True, True, False]),
            "stellar_dominated": np.array([True, True, True, False]),
        }
    )
    m = metrics_for_successful_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        scenario_df=df,
        mean_acceleration_gap_edge=1e-10,
        galactic_latitude_degrees=5.0,
        morphology_class="Im",
    )
    assert AUDIT_FLAG_FAIL_LATITUDE in m["audit_flags"]
    assert AUDIT_FLAG_FAIL_RELIABLE_FRACTION in m["audit_flags"]
    assert AUDIT_FLAG_FAIL_DENSITY_CLAMPED in m["audit_flags"]
    assert AUDIT_FLAG_FAIL_STELLAR in m["audit_flags"]
    # All four blocking codes appear in the failure-codes subset.
    for code in (
        AUDIT_FLAG_FAIL_LATITUDE,
        AUDIT_FLAG_FAIL_RELIABLE_FRACTION,
        AUDIT_FLAG_FAIL_DENSITY_CLAMPED,
        AUDIT_FLAG_FAIL_STELLAR,
    ):
        assert code in m["reliability_failure_codes"]
    assert m["is_reliable"] is False
    # ``_increment_success`` must be False whenever any blocking flag is set.
    assert m["_increment_success"] is False


def test_skipped_physics_propagates_blocking_codes_to_failure_codes():
    """``metrics_for_skipped_physics`` must propagate every blocking code into
    ``reliability_failure_codes`` (the column tools use to filter failures)."""
    m = metrics_for_skipped_physics(
        galaxy_name="G",
        scenario_name="s",
        gradient_scale_kpc=0.1,
        ionization_fraction=0.01,
        mean_acceleration_gap_edge=None,
        galactic_latitude_degrees=2.0,
        morphology_class="Im",
    )
    assert AUDIT_FLAG_FAIL_SKIPPED_PHYSICS in m["reliability_failure_codes"]
    assert AUDIT_FLAG_FAIL_LATITUDE in m["reliability_failure_codes"]
    assert m["is_reliable"] is False


def test_aggregate_galaxy_audit_flags_inherits_scenario_codes():
    """Galaxy-level aggregator inherits scenario codes verbatim (no double-append)."""
    flags, blocking, ok = aggregate_galaxy_audit_flags(
        scenario_audit_flags=[AUDIT_FLAG_FAIL_LATITUDE],
        inclination_quality_passed=True,
        calibration_status=CALIBRATION_STATUS_OPTIMAL,
        suboptimal_calibration_status=CALIBRATION_STATUS_SUBOPTIMAL,
    )
    assert flags == [AUDIT_FLAG_FAIL_LATITUDE]
    assert blocking == [AUDIT_FLAG_FAIL_LATITUDE]
    assert ok is False


def test_aggregate_galaxy_audit_flags_appends_inclination_failure():
    flags, blocking, ok = aggregate_galaxy_audit_flags(
        scenario_audit_flags=[],
        inclination_quality_passed=False,
        calibration_status=CALIBRATION_STATUS_OPTIMAL,
        suboptimal_calibration_status=CALIBRATION_STATUS_SUBOPTIMAL,
    )
    assert flags == [AUDIT_FLAG_FAIL_INCLINATION]
    assert AUDIT_FLAG_FAIL_INCLINATION in blocking
    assert ok is False


def test_aggregate_galaxy_audit_flags_appends_calibration_failure():
    flags, blocking, ok = aggregate_galaxy_audit_flags(
        scenario_audit_flags=[AUDIT_FLAG_WARN_B_BAND],
        inclination_quality_passed=True,
        calibration_status=CALIBRATION_STATUS_SUBOPTIMAL,
        suboptimal_calibration_status=CALIBRATION_STATUS_SUBOPTIMAL,
    )
    assert AUDIT_FLAG_WARN_B_BAND in flags
    assert AUDIT_FLAG_FAIL_CALIBRATION in flags
    # ``warn_b_band`` does NOT block; ``fail_calibration`` does.
    assert AUDIT_FLAG_WARN_B_BAND not in blocking
    assert AUDIT_FLAG_FAIL_CALIBRATION in blocking
    assert ok is False


def test_aggregate_galaxy_audit_flags_no_double_append_when_already_present():
    """If the scenario list already contains ``fail_inclination``, the aggregator must not duplicate it."""
    flags, _, _ = aggregate_galaxy_audit_flags(
        scenario_audit_flags=[AUDIT_FLAG_FAIL_INCLINATION],
        inclination_quality_passed=False,
        calibration_status=CALIBRATION_STATUS_OPTIMAL,
        suboptimal_calibration_status=CALIBRATION_STATUS_SUBOPTIMAL,
    )
    assert flags.count(AUDIT_FLAG_FAIL_INCLINATION) == 1


def test_aggregate_galaxy_audit_flags_reliable_when_all_gates_pass():
    flags, blocking, ok = aggregate_galaxy_audit_flags(
        scenario_audit_flags=[],
        inclination_quality_passed=True,
        calibration_status=CALIBRATION_STATUS_OPTIMAL,
        suboptimal_calibration_status=CALIBRATION_STATUS_SUBOPTIMAL,
    )
    assert flags == []
    assert blocking == []
    assert ok is True


def test_aggregate_galaxy_audit_flags_appends_ned_unresolved():
    flags, blocking, ok = aggregate_galaxy_audit_flags(
        scenario_audit_flags=[],
        inclination_quality_passed=True,
        calibration_status=CALIBRATION_STATUS_OPTIMAL,
        suboptimal_calibration_status=CALIBRATION_STATUS_SUBOPTIMAL,
        ned_unresolved=True,
    )
    assert flags == [AUDIT_FLAG_FAIL_NED_RESOLUTION]
    assert blocking == [AUDIT_FLAG_FAIL_NED_RESOLUTION]
    assert ok is False
