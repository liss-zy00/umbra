"""Per-galaxy orchestration: SPARC → NED → Planck → physics grid → persistence under ``galaxy_analysis/``."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from src.models.galaxy import Galaxy
from src.physics import PhysicsEngine, edge_radial_mean
from src.config.thresholds import (
    ASTROPHYSICAL_DEFAULTS,
    AUDIT_BLOCKING_NOTE_EN,
    AUDIT_WARN_NOTE_EN,
    CALIBRATION_REFERENCE_METHOD_DYNAMIC,
    CALIBRATION_REFERENCE_METHOD_FALLBACK,
    CALIBRATION_STATUS_OPTIMAL,
    CALIBRATION_STATUS_SUBOPTIMAL,
    PLANCK_CALIBRATION,
    RELIABILITY_CUTS,
    UNCERTAINTY,
)
from src.pipeline.constants import SPARC_CATALOG_JSON
from src.pipeline.per_galaxy.metrics_builder import (
    AUDIT_FLAG_FAIL_CALIBRATION,
    AUDIT_FLAG_FAIL_INCLINATION,
    aggregate_galaxy_audit_flags,
)
from src.pipeline.per_galaxy.scenario_grid_runner import run_scenario_grid
from src.pipeline.per_galaxy.scenario_selector import select_best_scenario
from src.providers.planck import PlanckProvider, resolve_stokes_amp_reference
from src.providers.sparc import SparcMorphologyProvider, SparcProvider
from src.providers.ned import NedProvider, raw_status_to_ned_status
from src.providers.edd_provider import apply_edd_distance_over_rotmod, build_edd_provenance_log_dict
from src.reporting.galaxy_log import (
    emit_galaxy_provenance_block,
    emit_validation_provenance_stub,
)

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from src.reporting.persistence.exporter import ResultPersister


@dataclass
class PipelineDependencies:
    """Injectable services for ``process_sparc_file`` (tests swap mocks here)."""

    sparc: SparcProvider
    morphology_provider: SparcMorphologyProvider
    morphology_lookup: dict[str, dict[str, Any]]
    ned: NedProvider
    planck: PlanckProvider
    persister: ResultPersister
    fits_summary: dict[str, Any]
    # Optional SPARC ``galaxy_name -> (ra_deg, dec_deg)`` ICRS seeds. Loaded once from
    # ``NED_ICRS_CONE_SEEDS_JSON`` by the runner; empty when the file is
    # missing or the ``positions`` section is empty. When a key matches, the cone
    # search inside ``NedProvider`` runs even if every name-based lookup fails.
    icrs_supplement: dict[str, tuple[float, float]] = field(default_factory=dict)
    # Optional EDD-backed luminosity distances (``None`` when CSV drops are absent).
    edd: Any | None = None
    # ``SparcMorphologyProvider.normalize_galaxy_name`` → HyperLEDA ``pgc`` for EDD joins.
    pgc_by_sparc_name: dict[str, int] = field(default_factory=dict)
    # When False (default interactive), per-galaxy provenance terminal blocks are suppressed.
    emit_verbose_galaxy_logs: bool = False
    # Galaxies that already printed the full five-line provenance block in Discovery (--logs);
    # validation reuses a one-line stub instead of repeating identical provider rows.
    provenance_full_log_galaxy_names: set[str] = field(default_factory=set)
    # Set by the runner before the validation grid loop; one cohort banner plus a single CR-based progress line.
    validation_grid_total_galaxies: int | None = None
    validation_discovery_reliable_count: int | None = None
    validation_completed_names: list[str] = field(default_factory=list)


@dataclass
class MutablePipelineState:
    """Accumulator mutated across galaxies (inventory rows, scenario grid, audit roll-up)."""

    master_inventory_rows: list[dict[str, Any]] = field(default_factory=list)
    supreme_rows: list[dict[str, Any]] = field(default_factory=list)
    scenario_rows: dict[tuple[float, float], list[dict[str, Any]]] = field(default_factory=dict)
    scenario_success_counts: Any | None = None
    # Audit summary — every processed galaxy contributes one entry with its name, latitude,
    # morphology and ``audit_flags`` (empty list = reliable). Drives the failure-mode breakdown
    # in the consistency audit report.
    galaxy_audit_summaries: list[dict[str, Any]] = field(default_factory=list)


def process_sparc_file(deps: PipelineDependencies, state: MutablePipelineState, file_path: Path) -> None:
    """Run the full analysis chain for one ``*_rotmod.dat`` path."""
    galaxy_name = file_path.stem.replace("_rotmod", "")

    rotation_dataframe, distance_mpc = deps.sparc.load_rotation_curve(file_path)
    if not np.any(np.isfinite(rotation_dataframe["gas_velocity_kms"].to_numpy(dtype=float))):
        msg = f"[{galaxy_name}] Missing valid gas velocity data - skipped."
        if deps.emit_verbose_galaxy_logs:
            print(msg, flush=True)
        else:
            logger.warning("%s", msg)
        return

    # Resolve morphology now (we need ``Inc_deg`` to propagate the kinematic systematic onto the
    # rotation dataframe before any acceleration math runs).
    morphology_key = deps.morphology_provider.normalize_galaxy_name(galaxy_name)
    morphology_entry = deps.morphology_lookup.get(morphology_key) or deps.morphology_lookup.get(
        galaxy_name.replace(" ", "").upper(), {}
    )
    if not morphology_entry:
        morphology_entry = {}
    morphology_class = str(morphology_entry.get("morphology_class", "Unclassified"))
    inclination_deg = morphology_entry.get("inclination_deg")
    inclination_uncertainty_deg = morphology_entry.get("inclination_uncertainty_deg")
    vflat_kms = morphology_entry.get("vflat_kms")
    vflat_uncertainty_kms = morphology_entry.get("vflat_uncertainty_kms")
    # PA + distance uncertainty supplements (Phase-1 refinement): PA comes from the curated
    # HyperLEDA supplement merged by ``SparcMorphologyProvider``; distance error from the
    # SPARC ``e_D_Mpc`` column. Missing values gracefully degrade — the engine falls back to
    # the axisymmetric gradient and emits ``warn_no_pa``; MC propagation reports the point
    # estimate without a band when uncertainties are absent.
    position_angle_deg = morphology_entry.get("position_angle_deg")
    position_angle_uncertainty_deg = morphology_entry.get("position_angle_uncertainty_deg")
    distance_uncertainty_mpc = morphology_entry.get("distance_uncertainty_mpc")
    hi_mass_1e9_solar = morphology_entry.get("hi_mass_1e9_solar")
    luminosity_3p6um_1e9_solar = morphology_entry.get("luminosity_3p6um_1e9_solar")

    rotation_dataframe = SparcProvider.attach_inclination_uncertainty(
        rotation_dataframe,
        inclination_deg=inclination_deg if isinstance(inclination_deg, (int, float)) else None,
        inclination_uncertainty_deg=(
            inclination_uncertainty_deg if isinstance(inclination_uncertainty_deg, (int, float)) else None
        ),
    )

    galaxy = Galaxy(
        name=galaxy_name,
        distance_mpc=distance_mpc,
        distance_uncertainty_mpc=(
            float(distance_uncertainty_mpc)
            if isinstance(distance_uncertainty_mpc, (int, float)) and float(distance_uncertainty_mpc) > 0
            else None
        ),
        rotation_dataframe=rotation_dataframe,
        inclination_deg=float(inclination_deg) if isinstance(inclination_deg, (int, float)) else None,
        inclination_uncertainty_deg=(
            float(inclination_uncertainty_deg)
            if isinstance(inclination_uncertainty_deg, (int, float))
            else None
        ),
        position_angle_deg=(
            float(position_angle_deg)
            if isinstance(position_angle_deg, (int, float)) and np.isfinite(float(position_angle_deg))
            else None
        ),
        position_angle_uncertainty_deg=(
            float(position_angle_uncertainty_deg)
            if isinstance(position_angle_uncertainty_deg, (int, float))
            and np.isfinite(float(position_angle_uncertainty_deg))
            else None
        ),
        vflat_kms=float(vflat_kms) if isinstance(vflat_kms, (int, float)) and float(vflat_kms) > 0 else None,
        vflat_uncertainty_kms=(
            float(vflat_uncertainty_kms)
            if isinstance(vflat_uncertainty_kms, (int, float)) and float(vflat_uncertainty_kms) > 0
            else None
        ),
    )
    icrs_seed = deps.icrs_supplement.get(galaxy_name)
    ned_data = deps.ned.fetch_metadata(galaxy_name, icrs_deg=icrs_seed)
    if ned_data.get("distance_mpc") is None:
        ned_data["distance_mpc"] = distance_mpc
        ned_data["distance_source"] = "SPARC header fallback"
    galaxy.right_ascension_degrees = ned_data.get("right_ascension_degrees")
    galaxy.declination_degrees = ned_data.get("declination_degrees")
    galaxy.redshift = ned_data.get("redshift")
    prov_high_precision_dist = apply_edd_distance_over_rotmod(
        deps,
        morphology_key=morphology_key,
        galaxy_name=galaxy_name,
        rotmod_distance_mpc=float(distance_mpc),
        galaxy=galaxy,
    )
    if not prov_high_precision_dist:
        galaxy.distance_mpc = distance_mpc

    ned_unresolved = ned_data.get("raw_status") == "unresolved"
    ned_status_label = raw_status_to_ned_status(ned_data.get("raw_status"))

    stokes_summary: dict[str, Any] = {"q_mean": None, "u_mean": None, "n_pixels_used": 0, "status": "not_attempted"}
    if galaxy.right_ascension_degrees is not None and galaxy.declination_degrees is not None:
        stokes_summary = deps.planck.extract_stokes_qu_mean(
            galaxy.right_ascension_degrees,
            galaxy.declination_degrees,
            galaxy_name=galaxy_name,
        )

    # Physical Δx: angular separation × galaxy distance. Reported in kpc and forwarded to the
    # engine so it can compute a parallel a_B based on the *real* sky-projected distance instead
    # of the scenario-grid Δx knob.
    mean_neighbor_separation_rad = stokes_summary.get("mean_neighbor_separation_rad")
    physical_neighbor_separation_kpc: float | None = None
    if (
        isinstance(mean_neighbor_separation_rad, (int, float))
        and np.isfinite(float(mean_neighbor_separation_rad))
        and float(mean_neighbor_separation_rad) > 0
        and galaxy.distance_mpc is not None
        and np.isfinite(float(galaxy.distance_mpc))
        and float(galaxy.distance_mpc) > 0
    ):
        # 1 Mpc = 1000 kpc; θ in rad, distance in Mpc → physical separation in Mpc → kpc.
        physical_neighbor_separation_kpc = float(
            float(galaxy.distance_mpc) * float(mean_neighbor_separation_rad) * 1000.0
        )

    if deps.emit_verbose_galaxy_logs and galaxy_name in deps.provenance_full_log_galaxy_names:
        if not isinstance(deps, PipelineDependencies):
            emit_validation_provenance_stub(galaxy_name, enabled=True)
        # Real deps: Discovery already printed the full block — no duplicate provenance at grid start.
    elif deps.emit_verbose_galaxy_logs:
        emit_galaxy_provenance_block(
            galaxy_name,
            enabled=True,
            sparc={
                "source_file": file_path.name,
                "distance_mpc": distance_mpc,
                "n_shells": int(len(rotation_dataframe)),
            },
            morph={
                "morphology_class": morphology_class,
                "hubble_type": morphology_entry.get("hubble_type"),
                "inclination_deg": galaxy.inclination_deg,
                "vflat_kms": galaxy.vflat_kms,
            },
            edd=build_edd_provenance_log_dict(
                deps,
                morphology_key=morphology_key,
                galaxy_name=galaxy_name,
                galaxy=galaxy,
            ),
            ned={
                "status": ned_status_label,
                "right_ascension_degrees": galaxy.right_ascension_degrees,
                "declination_degrees": galaxy.declination_degrees,
                "redshift": galaxy.redshift,
                "object_type": ned_data.get("object_type"),
            },
            planck={
                "q_center": stokes_summary.get("q_center"),
                "u_center": stokes_summary.get("u_center"),
                "pixel_scale_kpc": physical_neighbor_separation_kpc,
            },
            pipeline_phase="validation",
        )

    base_engine = PhysicsEngine(galaxy)
    base_engine.calculate_force_gap()
    base_dataframe = galaxy.rotation_dataframe.copy()
    mean_acceleration_gap_edge_global = edge_radial_mean(base_dataframe["acceleration_gap"])

    assert state.scenario_success_counts is not None
    all_scenario_metrics = run_scenario_grid(
        galaxy_name=galaxy_name,
        base_dataframe=base_dataframe,
        galaxy_coords=galaxy,
        stokes_summary=stokes_summary,
        fits_summary=deps.fits_summary,
        morphology_class=morphology_class,
        scenario_rows=state.scenario_rows,
        scenario_success_counts=state.scenario_success_counts,
        physical_neighbor_separation_kpc=physical_neighbor_separation_kpc,
    )

    best = select_best_scenario(all_scenario_metrics)
    best_scale = float(best["gradient_scale_kpc"])
    best_ion = float(best["ionization_fraction"])
    best_name = str(best["scenario_name"])
    best_ratio = best.get("theoretical_gap_coverage_ratio")
    best_ratio_value = float(best_ratio) if best_ratio is not None and np.isfinite(float(best_ratio)) else np.nan

    best_galaxy = Galaxy(
        name=galaxy.name,
        distance_mpc=galaxy.distance_mpc,
        distance_uncertainty_mpc=galaxy.distance_uncertainty_mpc,
        right_ascension_degrees=galaxy.right_ascension_degrees,
        declination_degrees=galaxy.declination_degrees,
        redshift=galaxy.redshift,
        inclination_deg=galaxy.inclination_deg,
        inclination_uncertainty_deg=galaxy.inclination_uncertainty_deg,
        position_angle_deg=galaxy.position_angle_deg,
        position_angle_uncertainty_deg=galaxy.position_angle_uncertainty_deg,
        vflat_kms=galaxy.vflat_kms,
        vflat_uncertainty_kms=galaxy.vflat_uncertainty_kms,
        rotation_dataframe=base_dataframe.copy(),
    )
    best_engine = PhysicsEngine(best_galaxy)
    # Calibration provenance: ``S_ref`` must match the live Planck extraction dict first, then the
    # FITS-level summary from ``initialize()``. Never silently substitute the emergency constant
    # when ``stokes_summary`` already carries the dynamic map median from ``extract_stokes_qu_mean``.
    stokes_reference_used = resolve_stokes_amp_reference(stokes_summary, deps.fits_summary)
    reference_method = str(
        stokes_summary.get("reference_method")
        or deps.fits_summary.get("reference_method", CALIBRATION_REFERENCE_METHOD_FALLBACK)
    )
    calibration_status = (
        CALIBRATION_STATUS_OPTIMAL
        if reference_method == CALIBRATION_REFERENCE_METHOD_DYNAMIC
        else CALIBRATION_STATUS_SUBOPTIMAL
    )

    best_lorentz_summary = best_engine.calculate_magnetic_pressure_acceleration(
        stokes_summary.get("stokes_amp_center"),
        stokes_summary.get("stokes_amp_neighbors_mean"),
        best_scale,
        best_ion,
        stokes_amp_reference=stokes_reference_used,
        physical_neighbor_separation_kpc=physical_neighbor_separation_kpc,
        neighbour_ring=stokes_summary.get("neighbour_ring"),
        center_theta_rad=stokes_summary.get("center_theta_rad"),
        center_phi_rad=stokes_summary.get("center_phi_rad"),
        position_angle_deg=best_galaxy.position_angle_deg,
        position_angle_uncertainty_deg=best_galaxy.position_angle_uncertainty_deg,
        morphology_class=morphology_class,
    )
    best_df = best_galaxy.rotation_dataframe
    if "magnetic_pressure_acceleration" not in best_df.columns:
        best_df["magnetic_pressure_acceleration"] = np.nan
    if "reliable_measurement" not in best_df.columns:
        best_df["reliable_measurement"] = False
    best_galaxy.rotation_dataframe = best_df

    galactic_latitude_degrees = None
    if galaxy.right_ascension_degrees is not None and galaxy.declination_degrees is not None:
        galactic_latitude_degrees = float(
            SkyCoord(
                ra=galaxy.right_ascension_degrees * u.deg,
                dec=galaxy.declination_degrees * u.deg,
                frame="icrs",
            ).galactic.b.deg
        )

    # Density-clamp + stellar-dominance audit for the best scenario, surfaced into the per-galaxy
    # summary so downstream readers can see *why* a galaxy is (or is not) reliable.
    n_points_density_clamped = int(best_lorentz_summary.get("n_points_density_clamped", 0) or 0)
    n_points_total = int(best_lorentz_summary.get("n_points_total", 0) or 0)
    density_clamped_fraction = float(best_lorentz_summary.get("density_clamped_fraction", 0.0) or 0.0)
    n_points_stellar_dominated = int(best_lorentz_summary.get("n_points_stellar_dominated", 0) or 0)
    stellar_dominated_fraction = float(best_lorentz_summary.get("stellar_dominated_fraction", 0.0) or 0.0)

    # Inclination-quality flag: face-on systems (i below the configured floor) make the
    # systematic V_obs / sin(i) propagation unstable.
    inclination_quality_passed = (
        galaxy.inclination_deg is None
        or float(galaxy.inclination_deg) >= float(RELIABILITY_CUTS.MIN_INCLINATION_DEG_FOR_RELIABILITY)
    )

    # RAR per-galaxy summary statistic: median |residual| in dex on the reliable shells. This is
    # the cleanest one-number indicator of how far this galaxy sits from the McGaugh fit.
    rar_residual_arr = (
        best_galaxy.rotation_dataframe.get("rar_residual_log10").to_numpy(dtype=float)
        if "rar_residual_log10" in best_galaxy.rotation_dataframe.columns
        else np.array([], dtype=float)
    )
    rar_residual_arr = rar_residual_arr[np.isfinite(rar_residual_arr)]
    rar_median_abs_residual_dex = (
        float(np.median(np.abs(rar_residual_arr))) if rar_residual_arr.size > 0 else None
    )

    # BTFR per-galaxy distance-to-ridge (McGaugh+2012). Uses M_HI from the SPARC catalogue and
    # the 3.6 µm luminosity with the canonical Υ_{3.6} = 0.5 M_⊙/L_⊙. Outliers (|residual| above
    # the configured threshold) get an informational ``warn_btfr_outlier`` audit flag.
    from src.physics.diagnostics.btfr import btfr_residual_dex, is_btfr_outlier

    btfr_residual_value = btfr_residual_dex(
        vflat_kms=galaxy.vflat_kms,
        hi_mass_1e9_solar=hi_mass_1e9_solar,
        luminosity_3p6um_1e9_solar=luminosity_3p6um_1e9_solar,
    )
    btfr_is_outlier = is_btfr_outlier(btfr_residual_value)

    # MC band on the headline coverage ratio. Propagates (D, i, V_flat) systematics through
    # the analytic ``base_ratio_callable`` so the engine does not need to re-run N times.
    from src.physics.diagnostics.uncertainty import propagate_gap_coverage_band

    base_ratio = best.get("theoretical_gap_coverage_ratio")
    base_ratio_f = float(base_ratio) if base_ratio is not None and np.isfinite(float(base_ratio)) else float("nan")
    mc_median, mc_p16, mc_p84 = propagate_gap_coverage_band(
        base_ratio=base_ratio_f,
        distance_mpc=galaxy.distance_mpc,
        e_distance_mpc=galaxy.distance_uncertainty_mpc,
        inclination_deg=galaxy.inclination_deg,
        e_inclination_deg=galaxy.inclination_uncertainty_deg,
        vflat_kms=galaxy.vflat_kms,
        e_vflat_kms=galaxy.vflat_uncertainty_kms,
    )
    # WIM convergence band: same MC reshaping rescales ``convergence_ratio_to_wim`` since it is
    # proportional to ``rho_req`` which inherits the same linear scaling in (D, i, V_flat).
    # We forward the same multiplicative factor as on ``theoretical_gap_coverage_ratio``.
    convergence_p16: float | None = None
    convergence_p84: float | None = None
    # ``theoretical_gap_coverage_ratio_full_lorentz`` is supplied by the metrics_builder via
    # the best-scenario row (computed from the per-shell ``accel_full_lorentz_m_s2`` column).
    theoretical_gap_coverage_ratio_full_lorentz = best.get("theoretical_gap_coverage_ratio_full_lorentz")
    accel_pressure_only_m_s2 = best.get("accel_pressure_only_m_s2")
    accel_full_lorentz_m_s2 = best.get("accel_full_lorentz_m_s2")

    # Galaxy-level audit flags — composition delegated to ``aggregate_galaxy_audit_flags``
    # (single source of truth for the warn_*-vs-fail_* convention and the inclination /
    # calibration appends). Keeping the call here means the orchestration logic stays in
    # ``galaxy_processor`` while the *policy* lives in ``metrics_builder``.
    galaxy_audit_flags, galaxy_blocking_flags, galaxy_is_reliable = aggregate_galaxy_audit_flags(
        scenario_audit_flags=best.get("audit_flags"),
        inclination_quality_passed=inclination_quality_passed,
        calibration_status=calibration_status,
        suboptimal_calibration_status=CALIBRATION_STATUS_SUBOPTIMAL,
        ned_unresolved=ned_unresolved,
    )

    # Append the per-galaxy BTFR outlier flag (warning-only — does NOT block reliability;
    # caller treats it as informational, mirroring ``warn_b_band``).
    from src.config.thresholds import AUDIT_CODES as _AC

    if btfr_is_outlier and _AC.WARN_BTFR_OUTLIER not in galaxy_audit_flags:
        galaxy_audit_flags.append(_AC.WARN_BTFR_OUTLIER)

    reliability_note = "Reliable for synthesis."
    if galactic_latitude_degrees is not None and abs(galactic_latitude_degrees) < RELIABILITY_CUTS.MIN_GALACTIC_LATITUDE_DEG:
        reliability_note = "High probability of foreground contamination from Milky Way dust"
    elif not inclination_quality_passed:
        reliability_note = (
            f"Inclination {galaxy.inclination_deg:.1f}° below the "
            f"{RELIABILITY_CUTS.MIN_INCLINATION_DEG_FOR_RELIABILITY:.0f}° floor — V_obs / sin(i) "
            "deprojection is unstable for this galaxy."
        )
    elif not galaxy_is_reliable:
        codes = (
            list(galaxy_blocking_flags)
            if galaxy_blocking_flags
            else [c for c in galaxy_audit_flags if str(c).startswith("fail_")]
        )
        if codes:
            reliability_note = " ".join(
                AUDIT_BLOCKING_NOTE_EN.get(str(code), f"Blocking gate triggered ({code}).") for code in codes
            )
        else:
            reliability_note = "Validation marked not reliable; inspect audit_flags for detail."

    if galaxy_is_reliable:
        info_codes = [str(c) for c in galaxy_audit_flags if str(c).startswith("warn_")]
        if info_codes:
            warn_text = "; ".join(
                AUDIT_WARN_NOTE_EN.get(code, f"Informational audit tag ({code}).") for code in info_codes
            )
            reliability_note = (
                "Reliable for synthesis (no blocking audit gates). "
                "These informational flags do not invalidate the headline metrics: "
                f"{warn_text}"
            )
        else:
            reliability_note = "Reliable for synthesis (no blocking audit gates; audit_flags empty)."

    if (
        np.isfinite(best_ratio_value)
        and RELIABILITY_CUTS.PHYSICAL_INTERPRETATION_STRONG_MIN
        <= best_ratio_value
        <= RELIABILITY_CUTS.PHYSICAL_INTERPRETATION_STRONG_MAX
    ):
        physical_interpretation = "Strong electromagnetic dominance candidate"
    elif np.isfinite(best_ratio_value) and best_ratio_value > RELIABILITY_CUTS.PHYSICAL_INTERPRETATION_SPIKE_MIN:
        physical_interpretation = "Potential foreground noise or localized plasma spike"
    else:
        physical_interpretation = "Intermediate regime requiring additional observational constraints"

    compact_scenarios = [
        {
            "scenario_name": s.get("scenario_name"),
            "theoretical_gap_coverage_ratio": s.get("theoretical_gap_coverage_ratio"),
            "pearson_correlation": s.get("pearson_correlation"),
        }
        for s in all_scenario_metrics
    ]

    galaxy_summary = {
        "galaxy_identity": {
            "galaxy_name": galaxy_name,
            "morphology_class": morphology_class,
            "hubble_t_type": morphology_entry.get("hubble_type"),
            "distance_local_mpc": galaxy.distance_mpc,
            "sparc_source_file": file_path.name,
            "sparc_catalog_file": SPARC_CATALOG_JSON.name,
        },
        "astrophysical_metadata": {
            "metadata_provider": "NASA/IPAC Extragalactic Database (NED)",
            "right_ascension_degrees": galaxy.right_ascension_degrees,
            "declination_degrees": galaxy.declination_degrees,
            "redshift": galaxy.redshift,
            "distance_mpc": ned_data.get("distance_mpc"),
            "distance_source": ned_data.get("distance_source"),
            "name_resolution_raw_status": ned_data.get("raw_status"),
        },
        "planck_instrumental_source": {
            "instrument": "Planck polarization maps",
            "fits_filename": deps.fits_summary.get("fits_filename"),
            "nside": deps.fits_summary.get("nside"),
            "unit": deps.fits_summary.get("unit"),
            "stokes_amplitude_reference": deps.fits_summary.get("stokes_amp_reference"),
            "healpy_ready": bool(deps.fits_summary.get("healpy_ready"))
            or bool(stokes_summary.get("healpy_ready")),
            "stokes_local_measurement": stokes_summary,
        },
        "calibration_metadata": {
            # If it is not exact, the report must say so.
            "stokes_reference_used": stokes_reference_used,
            "calibration_status": calibration_status,
            "reference_method": reference_method,
            # Audit constant only — not applied when calibration_status is OPTIMAL. Present so
            # reviewers can compare against stokes_reference_used without hunting the codebase.
            "stokes_reference_emergency_value": float(
                PLANCK_CALIBRATION.STOKES_AMPLITUDE_REFERENCE_EMERGENCY_ONLY
            ),
        },
        "density_audit": {
            "n_points_density_clamped": n_points_density_clamped,
            "n_points_total": n_points_total,
            "density_clamped_fraction": density_clamped_fraction,
            "max_clamped_fraction_for_reliability": float(
                RELIABILITY_CUTS.MAX_DENSITY_CLAMPED_FRACTION_FOR_RELIABILITY
            ),
            "n_points_stellar_dominated": n_points_stellar_dominated,
            "stellar_dominated_fraction": stellar_dominated_fraction,
            "max_stellar_dominance_fraction_for_reliability": float(
                RELIABILITY_CUTS.MAX_STELLAR_DOMINANCE_FRACTION_FOR_RELIABILITY
            ),
        },
        "kinematic_metadata": {
            "inclination_deg": galaxy.inclination_deg,
            "inclination_uncertainty_deg": galaxy.inclination_uncertainty_deg,
            "min_inclination_deg_for_reliability": float(
                RELIABILITY_CUTS.MIN_INCLINATION_DEG_FOR_RELIABILITY
            ),
            "inclination_quality_passed": bool(inclination_quality_passed),
            "position_angle_deg": galaxy.position_angle_deg,
            "position_angle_uncertainty_deg": galaxy.position_angle_uncertainty_deg,
            "position_angle_source": morphology_entry.get("position_angle_source"),
            "vflat_kms": galaxy.vflat_kms,
            "vflat_uncertainty_kms": galaxy.vflat_uncertainty_kms,
            "distance_uncertainty_mpc": galaxy.distance_uncertainty_mpc,
            "hi_mass_1e9_solar": morphology_entry.get("hi_mass_1e9_solar"),
            "luminosity_3p6um_1e9_solar": morphology_entry.get("luminosity_3p6um_1e9_solar"),
        },
        "rar_diagnostics": {
            "a0_m_s2": float(ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2),
            "median_absolute_residual_dex": rar_median_abs_residual_dex,
        },
        "btfr_diagnostics": {
            "ridge_slope": float(ASTROPHYSICAL_DEFAULTS.BTFR_RIDGE_SLOPE),
            "ridge_intercept": float(ASTROPHYSICAL_DEFAULTS.BTFR_RIDGE_INTERCEPT),
            "residual_dex": btfr_residual_value if np.isfinite(btfr_residual_value) else None,
            "is_outlier": bool(btfr_is_outlier),
            "outlier_threshold_dex": float(ASTROPHYSICAL_DEFAULTS.BTFR_OUTLIER_THRESHOLD_DEX),
        },
        "uncertainty_propagation": {
            "theoretical_gap_coverage_ratio_median": (mc_median if np.isfinite(mc_median) else None),
            "theoretical_gap_coverage_ratio_p16": (mc_p16 if np.isfinite(mc_p16) else None),
            "theoretical_gap_coverage_ratio_p84": (mc_p84 if np.isfinite(mc_p84) else None),
            "mc_samples": int(UNCERTAINTY.MC_DEFAULT_SAMPLES),
            "mc_seed": int(UNCERTAINTY.MC_DEFAULT_SEED),
        },
        "plasma_diagnostics": {
            "magnetic_field_microgauss": best_lorentz_summary.get("magnetic_field_microgauss"),
            "magnetic_field_within_plausibility_band": best_lorentz_summary.get(
                "magnetic_field_within_plausibility_band"
            ),
            "magnetic_field_plausibility_min_microgauss": float(
                RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MIN_MICROGAUSS
            ),
            "magnetic_field_plausibility_max_microgauss": float(
                RELIABILITY_CUTS.MAGNETIC_FIELD_PLAUSIBILITY_MAX_MICROGAUSS
            ),
            "plasma_beta_thermal_median": best_lorentz_summary.get("plasma_beta_thermal_median"),
            "plasma_beta_above_unity_fraction": best_lorentz_summary.get("plasma_beta_above_unity_fraction"),
            "plasma_beta_below_unity_fraction": best_lorentz_summary.get("plasma_beta_below_unity_fraction"),
            "magnetic_dominance_flag": bool(best_lorentz_summary.get("magnetic_dominance_flag", False)),
            "max_plasma_beta_for_mhd_dominance": float(
                RELIABILITY_CUTS.MAX_PLASMA_BETA_FOR_MHD_DOMINANCE
            ),
            "alfven_velocity_median_kms": best_lorentz_summary.get("alfven_velocity_median_kms"),
            "thermal_temperature_assumed_K": float(ASTROPHYSICAL_DEFAULTS.PLASMA_THERMAL_TEMPERATURE_K),
        },
        "gradient_scale_audit": {
            "scenario_gradient_scale_kpc": best_scale,
            "physical_neighbor_separation_kpc": physical_neighbor_separation_kpc,
            "scenario_to_physical_ratio": (
                float(best_scale / physical_neighbor_separation_kpc)
                if physical_neighbor_separation_kpc and physical_neighbor_separation_kpc > 0
                else None
            ),
            "magnetic_pressure_acceleration_physical_dx_median": best_lorentz_summary.get(
                "magnetic_pressure_acceleration_physical_dx"
            ),
        },
        "best_theoretical_scenario": {
            "scenario_name": best_name,
            "scenario_gradient_scale_kpc": best_scale,
            "scenario_ionization_fraction": best_ion,
            "theoretical_gap_coverage_ratio": best.get("theoretical_gap_coverage_ratio"),
            "theoretical_gap_coverage_ratio_full_lorentz": theoretical_gap_coverage_ratio_full_lorentz,
            "accel_pressure_only_m_s2": accel_pressure_only_m_s2,
            "accel_full_lorentz_m_s2": accel_full_lorentz_m_s2,
            "pearson_correlation": best.get("pearson_correlation"),
            "physical_interpretation": physical_interpretation,
            "required_density_for_match_median": best.get("required_density_for_match_median"),
            "resolution_correction_factor": best.get("resolution_correction_factor"),
            "calibrated_convergence_ratio": best.get("calibrated_convergence_ratio"),
            "magnetic_pressure_engine_summary": best_lorentz_summary,
        },
        "data_reliability_assessment": {
            # ``validation_is_reliable`` is the derived signal; the structured ``audit_flags`` list is the
            # source of truth — every code corresponds to a SPECIFIC, traceable failure reason
            # so consumers can answer "why was this galaxy excluded?" without re-running gates.
            "validation_is_reliable": bool(galaxy_is_reliable),
            "audit_flags": list(galaxy_audit_flags),
            "reliability_failure_codes": list(galaxy_blocking_flags),
            "galactic_latitude_degrees": galactic_latitude_degrees,
            "reliability_note": reliability_note,
        },
        "all_scenarios": compact_scenarios,
    }

    # Phase 3 only: logs bundle the **best** grid scenario's per-shell table + JSON (Discovery never picks a knob).
    galaxy_dir, csv_path, json_path = deps.persister.save_galaxy(best_galaxy, galaxy_summary)
    _ = (csv_path, json_path)  # paths reported in the final summary block, not per-galaxy

    state.master_inventory_rows.append(
        {
            "galaxy_name": galaxy_name,
            "right_ascension_degrees": galaxy.right_ascension_degrees,
            "declination_degrees": galaxy.declination_degrees,
            "galactic_latitude_degrees": galactic_latitude_degrees,
            "validation_is_reliable": bool(galaxy_is_reliable),
            "ned_status": ned_status_label,
            "audit_flags": list(galaxy_audit_flags),
            "reliability_failure_codes": list(galaxy_blocking_flags),
            "mean_acceleration_gap_edge": mean_acceleration_gap_edge_global,
            "morphology_class": morphology_class,
            "hubble_type": morphology_entry.get("hubble_type"),
            "note": reliability_note,
        }
    )

    # Audit roll-up — every processed galaxy contributes one entry, regardless of pass/fail,
    # so the consistency report can summarise failure-mode frequencies (latitude, stellar,
    # density, inclination, calibration) across the full processed sample.
    state.galaxy_audit_summaries.append(
        {
            "galaxy_name": galaxy_name,
            "morphology_class": morphology_class,
            "galactic_latitude_degrees": galactic_latitude_degrees,
            "inclination_deg": galaxy.inclination_deg,
            "validation_is_reliable": bool(galaxy_is_reliable),
            "audit_flags": list(galaxy_audit_flags),
            "reliability_failure_codes": list(galaxy_blocking_flags),
            "calibration_status": calibration_status,
        }
    )

    if deps.emit_verbose_galaxy_logs and isinstance(deps, PipelineDependencies):
        total_g = deps.validation_grid_total_galaxies
        if total_g is not None and total_g > 0:
            deps.validation_completed_names.append(galaxy_name)

    if galaxy_is_reliable and inclination_quality_passed:
        state.supreme_rows.append(
            {
                "galaxy_name": galaxy_name,
                "right_ascension_degrees": galaxy.right_ascension_degrees,
                "declination_degrees": galaxy.declination_degrees,
                "best_scenario_name": best_name,
                "scenario_gradient_scale_kpc": best_scale,
                "scenario_ionization_fraction": best_ion,
                "theoretical_gap_coverage_ratio": best["theoretical_gap_coverage_ratio"],
                "theoretical_gap_coverage_ratio_full_lorentz": theoretical_gap_coverage_ratio_full_lorentz,
                "theoretical_gap_coverage_ratio_p16": (mc_p16 if np.isfinite(mc_p16) else None),
                "theoretical_gap_coverage_ratio_p84": (mc_p84 if np.isfinite(mc_p84) else None),
                "accel_pressure_only_m_s2": accel_pressure_only_m_s2,
                "accel_full_lorentz_m_s2": accel_full_lorentz_m_s2,
                "btfr_residual_dex": btfr_residual_value if np.isfinite(btfr_residual_value) else None,
                "warn_btfr_outlier": bool(btfr_is_outlier),
                "warn_no_pa": galaxy.position_angle_deg is None,
                "warn_pure_pressure": morphology_class is not None
                and morphology_class.lower().startswith(("dwarf irreg", "blue compact")),
                "position_angle_deg": galaxy.position_angle_deg,
                "prov_high_precision_dist": bool(galaxy.metadata.get("prov_high_precision_dist")),
                "pearson_correlation": best["pearson_correlation"],
                "morphology_class": morphology_class,
                "calibration_status": calibration_status,
                "reference_method": reference_method,
                "stokes_reference_used": stokes_reference_used,
                "audit_flags": list(galaxy_audit_flags),
                "reliability_failure_codes": list(galaxy_blocking_flags),
                "n_points_density_clamped": n_points_density_clamped,
                "density_clamped_fraction": density_clamped_fraction,
                "n_points_stellar_dominated": n_points_stellar_dominated,
                "stellar_dominated_fraction": stellar_dominated_fraction,
                "inclination_deg": galaxy.inclination_deg,
                "inclination_quality_passed": inclination_quality_passed,
                "magnetic_field_microgauss": best_lorentz_summary.get("magnetic_field_microgauss"),
                "B_median_microGauss": best_lorentz_summary.get("magnetic_field_microgauss"),
                "magnetic_field_tesla_median": (
                    float(np.nanmedian(best_df["magnetic_field_tesla"].to_numpy(dtype=float)))
                    if "magnetic_field_tesla" in best_df.columns
                    and np.any(np.isfinite(best_df["magnetic_field_tesla"].to_numpy(dtype=float)))
                    else None
                ),
                "magnetic_field_within_plausibility_band": best_lorentz_summary.get(
                    "magnetic_field_within_plausibility_band"
                ),
                "plasma_beta_thermal_median": best_lorentz_summary.get("plasma_beta_thermal_median"),
                "plasma_beta_below_unity_fraction": best_lorentz_summary.get("plasma_beta_below_unity_fraction"),
                "magnetic_dominance_flag": bool(best_lorentz_summary.get("magnetic_dominance_flag", False)),
                "alfven_velocity_median_kms": best_lorentz_summary.get("alfven_velocity_median_kms"),
                "rar_median_abs_residual_dex": rar_median_abs_residual_dex,
                "physical_neighbor_separation_kpc": physical_neighbor_separation_kpc,
                "r_max_kpc": (
                    float(np.nanmax(best_df["radius_kpc"].to_numpy(dtype=float)))
                    if "radius_kpc" in best_df.columns
                    and np.any(np.isfinite(best_df["radius_kpc"].to_numpy(dtype=float)))
                    else None
                ),
                "required_density_for_match_median": best.get("required_density_for_match_median"),
                "resolution_correction_factor": best.get("resolution_correction_factor"),
                "calibrated_convergence_ratio": best.get("calibrated_convergence_ratio"),
            }
        )
