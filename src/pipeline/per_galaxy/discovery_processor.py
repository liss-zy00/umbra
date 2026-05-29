"""
Lean per-galaxy processor for **Phase 1 (Discovery census)** of the unified pipeline.

The discovery pass does not exercise the 9-scenario grid, the per-galaxy plots,
the audit report, or the population-level aggregation. It answers a single question:
*given the magnetic pressure gradient measured at the physical Planck pixel scale,
what plasma density ``rho_req`` would be required to close the acceleration gap?*

The expected empirical result is that ``median(rho_req) ~ 1e-27 kg/m^3``, i.e. the
density of the Warm Ionised Medium. Convergence to that value across the SPARC
catalogue is the falsifiable prediction documented in
:doc:`../../docs/results_interpretation/Scientific_Evidence_Narrative_EN`.

For each galaxy the helper:

1. Loads the SPARC rotation curve and the SPARC morphology row.
2. Resolves NED coordinates and inclination (required for Planck sampling).
3. Pulls the Stokes amplitudes at the galaxy's sky position from the Planck map.
4. Computes the **physical** projected pixel separation ``Δx_phys = D_galaxy · θ_neighbour``.
5. Runs the engine **once** with ``scenario_gradient_scale_kpc = Δx_phys`` (so the
   primary ``magnetic_pressure_gradient`` column carries the physical Δx) and
   ``scenario_ionization_fraction = 1.0`` (``rho_req`` is independent of f_ion;
   any positive value works).
6. Aggregates ``required_density_for_match`` over the shells where the engine could
   compute a finite, underperforming value, and returns a single dictionary.

A row is returned for **every** processed galaxy (including those that fail the
audit) so the discovery CSV doubles as a coverage inventory. Audit flags surface
the exclusion reason (``fail_lat``, ``fail_skipped_physics``, etc.).

The runner persists the resulting rows as four CSVs under ``results/step_discovery/``:

- ``discovery_all.csv`` — every processed galaxy.
- ``discovery_reliable.csv`` — only rows with ``discovery_is_reliable == True``.
- ``discovery_not_reliable.csv`` — complement of ``discovery_reliable`` (same schema; ``audit_flags`` explain why).
- ``discovery_cohort_rho_band.csv`` — internally consistent subset: reliable rows whose
  ``required_density_for_match`` lies inside ``[REL_MIN · median, REL_MAX · median]`` over the
  reliable subset, with ``REL_MIN`` / ``REL_MAX`` from
  :attr:`src.reporting.persistence.exporter.ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MIN`
  and ``DISCOVERY_ELITE_BAND_RELATIVE_MAX``.

All four files share the same schema, carrying the headline column
``required_density_for_match`` (kg/m³, the median across the shells where the
engine produced a finite, underperforming ``ρ_req``), the SI ``delta_r_m`` for
Δx_phys, and the convergence ratio against the Warm Ionised Medium reference
density.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS, AUDIT_CODES, RELIABILITY_CUTS
from src.models.galaxy import Galaxy
from src.physics import PhysicsEngine, beam_physical_size_kpc, resolution_correction_factor
from src.pipeline.per_galaxy.galaxy_processor import PipelineDependencies
from src.providers.edd_provider import apply_edd_distance_over_rotmod, build_edd_provenance_log_dict
from src.providers.ned import raw_status_to_ned_status
from src.providers.planck import resolve_stokes_amp_reference
from src.providers.sparc import SparcProvider
from src.reporting.galaxy_log import emit_galaxy_provenance_block

logger = logging.getLogger(__name__)

# Discovery audit codes that block entry into Phase 2 (grid search). Galaxies with only
# non-blocking issues (e.g. ``is_reliable == False`` because every shell is "overperforming"
# with |a_B| > |gap|) still pass: they are strong candidates and must receive the full grid.
DISCOVERY_BLOCKING_AUDIT_CODES: frozenset[str] = frozenset(
    {
        AUDIT_CODES.FAIL_LATITUDE,
        AUDIT_CODES.FAIL_SKIPPED_PHYSICS,
        AUDIT_CODES.FAIL_NED_RESOLUTION,
    }
)


def discovery_audits_clean(row: dict[str, Any]) -> bool:
    """
    Return ``True`` iff the Discovery row carries **no** blocking audit code.

    Galaxies where Discovery could not invert ``rho_req`` because ``|a_B| > |gap|`` on every
    shell (``is_reliable == False`` with an empty ``audit_flags``) intentionally pass: they
    still enter the validation grid.
    """
    flags = {f.strip() for f in str(row.get("audit_flags", "")).split(",") if f.strip()}
    return not (DISCOVERY_BLOCKING_AUDIT_CODES & flags)


def _galactic_latitude_deg(ra_deg: float | None, dec_deg: float | None) -> float | None:
    if ra_deg is None or dec_deg is None:
        return None
    return float(
        SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs").galactic.b.deg
    )


def _safe_float(value: Any) -> float:
    """Cast to ``float`` returning ``nan`` for ``None``/non-finite inputs."""
    if value is None:
        return float("nan")
    try:
        f = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return f if np.isfinite(f) else float("nan")


def _nan_row(
    galaxy_name: str,
    sparc_filename: str,
    galactic_latitude_deg: float | None,
    distance_mpc: float | None,
    audit_flags: list[str],
    note: str,
    ned_status: str = "NOT_FOUND_NED",
) -> dict[str, Any]:
    """Coverage-preserving row when the engine cannot be exercised."""
    return {
        "galaxy_name": galaxy_name,
        "sparc_source_file": sparc_filename,
        "n_shells_used": 0,
        # ``required_density_for_match`` is the headline discovery column. The
        # ``median_*`` alias is retained for backward compatibility with earlier
        # tooling that consumed the verbose name.
        "required_density_for_match": float("nan"),
        "median_required_density_for_match_kg_m3": float("nan"),
        "median_magnetic_pressure_acceleration_physical_dx_m_s2": float("nan"),
        "median_acceleration_gap_m_s2": float("nan"),
        "physical_pixel_scale_kpc": float("nan"),
        "delta_r_m": float("nan"),
        # Beam-smearing dilution factor ``k = D_beam / L_grad``. NaN here means the
        # diagnostic was not evaluable (no distance / no pixel scale resolved).
        "resolution_correction_factor": float("nan"),
        "beam_size_kpc": float("nan"),
        "convergence_ratio_to_wim": float("nan"),
        "wim_reference_density_kg_m3": float(ASTROPHYSICAL_DEFAULTS.WIM_REFERENCE_DENSITY_KG_M3),
        "galactic_latitude_deg": _safe_float(galactic_latitude_deg),
        "distance_mpc": _safe_float(distance_mpc),
        "audit_flags": ",".join(audit_flags) if audit_flags else "",
        "discovery_is_reliable": False,
        "ned_status": ned_status,
        "note": note,
    }


def process_galaxy_discovery(
    deps: PipelineDependencies,
    file_path: Path,
) -> dict[str, Any] | None:
    """
    Build one discovery row for the galaxy at ``file_path``. The row is appended
    to ``discovery_all.csv`` by the runner, additionally to
    ``discovery_reliable.csv`` when ``discovery_is_reliable`` is ``True``, to
    ``discovery_not_reliable.csv`` when ``False``, and
    further to ``discovery_cohort_rho_band.csv`` when the row is reliable and
    ``required_density_for_match`` lies inside the dynamic median-relative band
    defined on :class:`src.reporting.persistence.exporter.ResultPersister`.

    Returns ``None`` only when the SPARC file lacks a usable gas-velocity column
    (the validation branch silently skips those too); every other galaxy yields a
    row, even if the engine could not produce a finite ``rho_req`` — the audit
    flags then explain why.
    """
    galaxy_name = file_path.stem.replace("_rotmod", "")

    rotation_dataframe, distance_mpc = deps.sparc.load_rotation_curve(file_path)
    if not np.any(np.isfinite(rotation_dataframe["gas_velocity_kms"].to_numpy(dtype=float))):
        msg = f"[discovery][{galaxy_name}] missing valid gas velocity — skipped."
        if deps.emit_verbose_galaxy_logs:
            print(msg, flush=True)
        else:
            logger.warning("%s", msg)
        return None

    morphology_key = deps.morphology_provider.normalize_galaxy_name(galaxy_name)
    morphology_entry = deps.morphology_lookup.get(morphology_key) or deps.morphology_lookup.get(
        galaxy_name.replace(" ", "").upper(), {}
    ) or {}
    inclination_deg = morphology_entry.get("inclination_deg")
    inclination_uncertainty_deg = morphology_entry.get("inclination_uncertainty_deg")
    position_angle_deg = morphology_entry.get("position_angle_deg")
    position_angle_uncertainty_deg = morphology_entry.get("position_angle_uncertainty_deg")
    distance_uncertainty_mpc = morphology_entry.get("distance_uncertainty_mpc")
    morphology_class = morphology_entry.get("morphology_class")

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
    )

    icrs_seed = deps.icrs_supplement.get(galaxy_name)
    ned_data = deps.ned.fetch_metadata(galaxy_name, icrs_deg=icrs_seed)
    if ned_data.get("distance_mpc") is None:
        ned_data["distance_mpc"] = distance_mpc
    galaxy.right_ascension_degrees = ned_data.get("right_ascension_degrees")
    galaxy.declination_degrees = ned_data.get("declination_degrees")
    galaxy.redshift = ned_data.get("redshift")
    _ = apply_edd_distance_over_rotmod(
        deps,
        morphology_key=morphology_key,
        galaxy_name=galaxy_name,
        rotmod_distance_mpc=float(distance_mpc),
        galaxy=galaxy,
    )
    if galaxy.distance_mpc is not None and np.isfinite(float(galaxy.distance_mpc)):
        distance_mpc = float(galaxy.distance_mpc)

    ned_status_label = raw_status_to_ned_status(ned_data.get("raw_status"))

    galactic_latitude_deg = _galactic_latitude_deg(
        galaxy.right_ascension_degrees, galaxy.declination_degrees
    )

    audit_flags: list[str] = []
    if galactic_latitude_deg is not None and abs(galactic_latitude_deg) < float(
        RELIABILITY_CUTS.MIN_GALACTIC_LATITUDE_DEG
    ):
        audit_flags.append(AUDIT_CODES.FAIL_LATITUDE)

    if galaxy.right_ascension_degrees is None or galaxy.declination_degrees is None:
        extra = [AUDIT_CODES.FAIL_SKIPPED_PHYSICS]
        if ned_data.get("raw_status") == "unresolved":
            extra.append(AUDIT_CODES.FAIL_NED_RESOLUTION)
        return _nan_row(
            galaxy_name,
            file_path.name,
            galactic_latitude_deg,
            distance_mpc,
            audit_flags + extra,
            note=(
                "NED and SIMBAD did not resolve sky coordinates for this SPARC label; Planck sampling is impossible."
                if ned_data.get("raw_status") == "unresolved"
                else "NED did not resolve sky coordinates; Planck sampling is impossible."
            ),
            ned_status=ned_status_label,
        )

    stokes_summary = deps.planck.extract_stokes_qu_mean(
        galaxy.right_ascension_degrees,
        galaxy.declination_degrees,
        galaxy_name=galaxy_name,
    )
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
        physical_neighbor_separation_kpc = float(
            float(galaxy.distance_mpc) * float(mean_neighbor_separation_rad) * 1000.0
        )

    emit_galaxy_provenance_block(
        galaxy_name,
        enabled=deps.emit_verbose_galaxy_logs,
        sparc={
            "source_file": file_path.name,
            "distance_mpc": distance_mpc,
            "n_shells": int(len(rotation_dataframe)),
        },
        morph={
            "morphology_class": morphology_entry.get("morphology_class"),
            "hubble_type": morphology_entry.get("hubble_type"),
            "inclination_deg": galaxy.inclination_deg,
            "vflat_kms": morphology_entry.get("vflat_kms"),
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
        pipeline_phase="discovery",
    )
    if deps.emit_verbose_galaxy_logs:
        deps.provenance_full_log_galaxy_names.add(galaxy_name)

    if physical_neighbor_separation_kpc is None:
        return _nan_row(
            galaxy_name,
            file_path.name,
            galactic_latitude_deg,
            distance_mpc,
            audit_flags + [AUDIT_CODES.FAIL_SKIPPED_PHYSICS],
            note="Planck pixel angular separation could not be resolved at this sky position.",
            ned_status=ned_status_label,
        )

    engine = PhysicsEngine(galaxy)
    engine.calculate_force_gap()
    stokes_reference_used = resolve_stokes_amp_reference(stokes_summary, deps.fits_summary)

    lorentz_summary = engine.calculate_magnetic_pressure_acceleration(
        stokes_summary.get("stokes_amp_center"),
        stokes_summary.get("stokes_amp_neighbors_mean"),
        # Anchor the scenario knob to the physical Planck pixel separation. The engine then
        # computes ``magnetic_pressure_gradient`` (and the dependent ``required_density_for_match``
        # column) at the same physical scale used for ``magnetic_pressure_acceleration_physical_dx``.
        physical_neighbor_separation_kpc,
        # f_ion = 1.0 (WIM is fully ionised). ``required_density_for_match`` is
        # |∇u_B| / |gap_abs| and does not depend on f_ion; any positive value is acceptable.
        1.0,
        stokes_amp_reference=stokes_reference_used,
        physical_neighbor_separation_kpc=physical_neighbor_separation_kpc,
        neighbour_ring=stokes_summary.get("neighbour_ring"),
        center_theta_rad=stokes_summary.get("center_theta_rad"),
        center_phi_rad=stokes_summary.get("center_phi_rad"),
        position_angle_deg=galaxy.position_angle_deg,
        position_angle_uncertainty_deg=galaxy.position_angle_uncertainty_deg,
        morphology_class=morphology_class,
    )

    if lorentz_summary.get("status") != "ok":
        return _nan_row(
            galaxy_name,
            file_path.name,
            galactic_latitude_deg,
            distance_mpc,
            audit_flags + [AUDIT_CODES.FAIL_SKIPPED_PHYSICS],
            note=str(lorentz_summary.get("status") or "engine reported no valid magnetic pressure shells"),
            ned_status=ned_status_label,
        )

    df = galaxy.rotation_dataframe
    rho_req_per_shell = (
        df["required_density_for_match"].to_numpy(dtype=float)
        if "required_density_for_match" in df.columns
        else np.array([], dtype=float)
    )
    rho_req_finite = rho_req_per_shell[np.isfinite(rho_req_per_shell) & (rho_req_per_shell > 0)]
    median_rho_req = float(np.median(rho_req_finite)) if rho_req_finite.size > 0 else float("nan")

    a_phys_per_shell = (
        df["magnetic_pressure_acceleration_physical_dx"].to_numpy(dtype=float)
        if "magnetic_pressure_acceleration_physical_dx" in df.columns
        else np.array([], dtype=float)
    )
    a_phys_finite = a_phys_per_shell[np.isfinite(a_phys_per_shell)]
    median_a_phys = float(np.median(np.abs(a_phys_finite))) if a_phys_finite.size > 0 else float("nan")

    gap_per_shell = (
        df["acceleration_gap"].to_numpy(dtype=float)
        if "acceleration_gap" in df.columns
        else np.array([], dtype=float)
    )
    gap_finite = gap_per_shell[np.isfinite(gap_per_shell)]
    median_gap = float(np.median(np.abs(gap_finite))) if gap_finite.size > 0 else float("nan")

    convergence_ratio = (
        median_rho_req / float(ASTROPHYSICAL_DEFAULTS.WIM_REFERENCE_DENSITY_KG_M3)
        if np.isfinite(median_rho_req)
        else float("nan")
    )

    is_reliable = (
        not audit_flags
        and np.isfinite(median_rho_req)
        and rho_req_finite.size > 0
    )

    # SI counterpart of ``physical_pixel_scale_kpc``: same number expressed in metres so
    # downstream analyses do not have to re-do the unit conversion to consume Δx_phys.
    delta_r_m: float = float("nan")
    if (
        physical_neighbor_separation_kpc is not None
        and np.isfinite(float(physical_neighbor_separation_kpc))
        and float(physical_neighbor_separation_kpc) > 0
    ):
        delta_r_m = float(
            (float(physical_neighbor_separation_kpc) * u.kpc).to_value(u.m)
        )

    # Beam-smearing diagnostic: the Planck FWHM projects to ``D_beam = D · θ_beam`` on
    # the disk plane; comparing that footprint against the gradient scale used by the
    # engine (the physical pixel separation in the Discovery path) yields the dilution
    # factor ``k = D_beam / L_grad``. ``k >> 1`` ⇒ the reconstructed ``∇u_B`` is
    # averaged out by a too-wide beam and ``ρ_req`` is over-estimated by the same factor.
    beam_size_kpc_value = beam_physical_size_kpc(distance_mpc=distance_mpc)
    correction_factor = resolution_correction_factor(
        distance_mpc=distance_mpc,
        gradient_scale_kpc=physical_neighbor_separation_kpc,
    )

    return {
        "galaxy_name": galaxy_name,
        "sparc_source_file": file_path.name,
        "n_shells_used": int(rho_req_finite.size),
        # ``required_density_for_match`` is the headline discovery column (kg/m^3).
        # ``median_required_density_for_match_kg_m3`` duplicates that value for legacy CSV
        # consumers only — new integrations should read ``required_density_for_match``.
        "required_density_for_match": median_rho_req,
        "median_required_density_for_match_kg_m3": median_rho_req,
        "median_magnetic_pressure_acceleration_physical_dx_m_s2": median_a_phys,
        "median_acceleration_gap_m_s2": median_gap,
        "physical_pixel_scale_kpc": _safe_float(physical_neighbor_separation_kpc),
        "delta_r_m": delta_r_m,
        "resolution_correction_factor": _safe_float(correction_factor),
        "beam_size_kpc": _safe_float(beam_size_kpc_value),
        "convergence_ratio_to_wim": convergence_ratio,
        "wim_reference_density_kg_m3": float(ASTROPHYSICAL_DEFAULTS.WIM_REFERENCE_DENSITY_KG_M3),
        "galactic_latitude_deg": _safe_float(galactic_latitude_deg),
        "distance_mpc": _safe_float(distance_mpc),
        "audit_flags": ",".join(audit_flags) if audit_flags else "",
        "discovery_is_reliable": bool(is_reliable),
        "ned_status": ned_status_label,
        "note": "" if is_reliable else (
            "Latitude cut applied (Milky Way foreground)" if AUDIT_CODES.FAIL_LATITUDE in audit_flags else ""
        ),
    }
