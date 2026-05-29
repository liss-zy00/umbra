"""Execute the (gradient scale × ionization fraction) grid for one galaxy."""

from __future__ import annotations

from typing import Any

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

from src.models.galaxy import Galaxy
from src.physics import GRADIENT_SCALES_KPC, IONIZATION_FRACTIONS, PhysicsEngine, edge_radial_mean, resolution_correction_factor
from src.providers.planck import resolve_stokes_amp_reference
from src.pipeline.per_galaxy.metrics_builder import metrics_for_skipped_physics, metrics_for_successful_physics


def galactic_latitude_degrees(ra_deg: float | None, dec_deg: float | None) -> float | None:
    if ra_deg is None or dec_deg is None:
        return None
    return float(
        SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs").galactic.b.deg
    )


def run_scenario_grid(
    *,
    galaxy_name: str,
    base_dataframe: Any,
    galaxy_coords: Galaxy,
    stokes_summary: dict[str, Any],
    fits_summary: dict[str, Any],
    morphology_class: str,
    scenario_rows: dict[tuple[float, float], list[dict[str, Any]]],
    scenario_success_counts: np.ndarray,
    physical_neighbor_separation_kpc: float | None = None,
) -> list[dict[str, Any]]:
    """Return ``all_scenario_metrics`` and mutate ``scenario_rows`` / ``scenario_success_counts``."""
    lat_deg = galactic_latitude_degrees(
        galaxy_coords.right_ascension_degrees,
        galaxy_coords.declination_degrees,
    )
    mean_acceleration_gap_edge = edge_radial_mean(base_dataframe["acceleration_gap"])

    all_scenario_metrics: list[dict[str, Any]] = []
    stokes_ref = resolve_stokes_amp_reference(stokes_summary, fits_summary)
    neighbour_ring = stokes_summary.get("neighbour_ring")
    center_theta_rad = stokes_summary.get("center_theta_rad")
    center_phi_rad = stokes_summary.get("center_phi_rad")

    for scale_idx, gradient_scale_kpc in enumerate(GRADIENT_SCALES_KPC):
        for ion_idx, ionization_fraction in enumerate(IONIZATION_FRACTIONS):
            scenario_name = f"scale_{gradient_scale_kpc:.3f}_ion_{ionization_fraction:.3f}"
            scenario_galaxy = Galaxy(
                name=galaxy_coords.name,
                distance_mpc=galaxy_coords.distance_mpc,
                distance_uncertainty_mpc=galaxy_coords.distance_uncertainty_mpc,
                right_ascension_degrees=galaxy_coords.right_ascension_degrees,
                declination_degrees=galaxy_coords.declination_degrees,
                redshift=galaxy_coords.redshift,
                inclination_deg=galaxy_coords.inclination_deg,
                inclination_uncertainty_deg=galaxy_coords.inclination_uncertainty_deg,
                position_angle_deg=galaxy_coords.position_angle_deg,
                position_angle_uncertainty_deg=galaxy_coords.position_angle_uncertainty_deg,
                vflat_kms=galaxy_coords.vflat_kms,
                vflat_uncertainty_kms=galaxy_coords.vflat_uncertainty_kms,
                rotation_dataframe=base_dataframe.copy(),
            )
            scenario_engine = PhysicsEngine(scenario_galaxy)
            lorentz_summary = scenario_engine.calculate_magnetic_pressure_acceleration(
                stokes_summary.get("stokes_amp_center"),
                stokes_summary.get("stokes_amp_neighbors_mean"),
                gradient_scale_kpc,
                ionization_fraction,
                stokes_amp_reference=stokes_ref,
                physical_neighbor_separation_kpc=physical_neighbor_separation_kpc,
                neighbour_ring=neighbour_ring,
                center_theta_rad=center_theta_rad,
                center_phi_rad=center_phi_rad,
                position_angle_deg=galaxy_coords.position_angle_deg,
                position_angle_uncertainty_deg=galaxy_coords.position_angle_uncertainty_deg,
                morphology_class=morphology_class,
            )
            scenario_df = scenario_galaxy.rotation_dataframe

            if lorentz_summary.get("status") == "ok":
                k_cell = float(
                    resolution_correction_factor(
                        distance_mpc=galaxy_coords.distance_mpc,
                        gradient_scale_kpc=gradient_scale_kpc,
                    )
                )
                rho_med = lorentz_summary.get("required_density_for_match")
                rho_med_f: float | None
                if rho_med is None:
                    rho_med_f = None
                else:
                    rho_med_f = float(rho_med) if np.isfinite(float(rho_med)) else None
                scenario_metrics = metrics_for_successful_physics(
                    galaxy_name=galaxy_name,
                    scenario_name=scenario_name,
                    gradient_scale_kpc=gradient_scale_kpc,
                    ionization_fraction=ionization_fraction,
                    scenario_df=scenario_df,
                    mean_acceleration_gap_edge=mean_acceleration_gap_edge,
                    galactic_latitude_degrees=lat_deg,
                    morphology_class=morphology_class,
                    engine_audit_flags=lorentz_summary.get("audit_flags") or [],
                    required_density_for_match_median=rho_med_f,
                    resolution_correction_factor_k=k_cell if np.isfinite(k_cell) else None,
                )
            else:
                scenario_metrics = metrics_for_skipped_physics(
                    galaxy_name=galaxy_name,
                    scenario_name=scenario_name,
                    gradient_scale_kpc=gradient_scale_kpc,
                    ionization_fraction=ionization_fraction,
                    mean_acceleration_gap_edge=mean_acceleration_gap_edge,
                    galactic_latitude_degrees=lat_deg,
                    morphology_class=morphology_class,
                )

            increment = bool(scenario_metrics.pop("_increment_success", False))
            if increment:
                scenario_success_counts[scale_idx, ion_idx] += 1

            scenario_rows[(gradient_scale_kpc, ionization_fraction)].append(scenario_metrics)
            all_scenario_metrics.append(scenario_metrics)

    return all_scenario_metrics
