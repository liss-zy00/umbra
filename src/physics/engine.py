from __future__ import annotations

"""
Magnetic **pressure-gradient** engine for galactic rotation curves.

This module hosts :class:`PhysicsEngine`, the orchestrator for per-galaxy MHD diagnostics.
Collaborators live in focused modules (``gas_density``, ``magnetic_field_calibration``,
``magnetic_pressure_gradient``, ``magnetic_gradient_strategy``).

Important naming note
---------------------
The public estimator is named ``calculate_magnetic_pressure_acceleration`` to match the term
it actually evaluates: the volume Lorentz force on a magnetized fluid reduces, under ideal-MHD
assumptions (no magnetic tension, isotropic pressure), to the magnetic-pressure-gradient
divided by the ionized mass density.

Calibration caveat
------------------
``MagneticFieldCalibrator`` anchors Stokes amplitudes to Tesla via ``√(2 μ₀ u_CMB)``. This uses
the CMB energy density only as a **dimensionful scale** — it is not a claim that Galactic fields
originate from the CMB.

Gradient-scale caveat
---------------------
Sky neighbors from HEALPix probe ~arcminute scales on the map, while ``scenario_gradient_scale_kpc``
sets **Δx** in the disk (50–500 pc). Mixing these scales is a deliberate exploratory shortcut;
results should be interpreted with strong caution (flagged in ``unreliable_reason`` when noise dominates).
"""

from typing import Any, Iterator

import numpy as np
from astropy import units as u

from src.config.thresholds import AUDIT_CODES, RELIABILITY_CUTS

from ..models.galaxy import Galaxy
from .defaults import (
    DEFAULT_STOKES_AMPLITUDE_REFERENCE,
    ENERGY_DENSITY_CMB,
    GRADIENT_SCALES_KPC,
    IONIZATION_FRACTIONS,
    MIN_GAS_MASS_DENSITY,
    MU0,
    RAR_A0_M_S2,
    SHELL_THICKNESS_CLIP_MAX_KPC,
    SHELL_THICKNESS_CLIP_MIN_KPC,
    SHELL_THICKNESS_FRACTION_OF_RMAX,
)
from .diagnostics.plasma import (
    TESLA_PER_MICROGAUSS,
    MagneticFieldPlausibilityCheck,
    PlasmaDiagnostics,
    ReliabilityMaskBuilder,
)
from .diagnostics.tension import (
    supports_logarithmic_spiral,
    tension_acceleration_radial,
)
from .gas_density import GasDensityEstimator
from .kinematics import BaryonicGapCalculator
from .magnetic_field_calibration import MagneticFieldCalibrator
from .magnetic_gradient_strategy import (
    PhysicalGradientInputs,
    resolve_physical_gradient_strategy,
)
from .magnetic_pressure_gradient import (
    MagneticPressureGradientEstimator,
)


def _finite_stokes_amp(value: Any) -> float | None:
    """Coerce Planck / JSON Stokes amplitude to a finite float, or ``None``."""
    if value is None:
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f if np.isfinite(f) else None


def iter_magnetic_scenarios() -> Iterator[tuple[float, float, str]]:
    """Yield ``(gradient_scale_kpc, ionization_fraction, scenario_name)`` for the full grid."""
    for g_kpc in GRADIENT_SCALES_KPC:
        for ion in IONIZATION_FRACTIONS:
            name = f"scale_{g_kpc:.3f}_ion_{ion:.3f}"
            yield g_kpc, ion, name


class PhysicsEngine:
    """
    SOLID composer for the per-galaxy MHD pipeline.

    The engine delegates each physical concern to a single-responsibility collaborator,
    injected through the constructor (open/closed + dependency-inversion):

    * :class:`BaryonicGapCalculator` — V_obs, V_bar → ``acceleration_gap`` (kinematics).
    * :class:`MagneticFieldCalibrator` — Stokes amplitudes → magnetic field strength.
    * :class:`GasDensityEstimator` — V_gas (and stellar components) → ``ρ_gas`` with
      stellar-dominance gate and density-floor audit.
    * :class:`MagneticPressureGradientEstimator` — finite-difference ``∇u_B`` and
      ``a = -∇P_B / ρ_ion``.
    * :class:`PlasmaDiagnostics` — plasma-β and Alfvén-velocity diagnostics.
    * :class:`MagneticFieldPlausibilityCheck` — μG plausibility-band gate.
    * :class:`ReliabilityMaskBuilder` — per-shell ``reliable_measurement`` mask.

    Public workflow:
        1. :py:meth:`calculate_force_gap` — SI accelerations and ``acceleration_gap``.
        2. :py:meth:`calculate_magnetic_pressure_acceleration` — Planck-anchored ``a_B``
           together with the full sanity-diagnostic block.
    """

    def __init__(
        self,
        galaxy: Galaxy,
        *,
        gap_calculator: BaryonicGapCalculator | None = None,
        field_calibrator: MagneticFieldCalibrator | None = None,
        gas_density: GasDensityEstimator | None = None,
        pressure_gradient: MagneticPressureGradientEstimator | None = None,
        plasma_diagnostics: PlasmaDiagnostics | None = None,
        plausibility_check: MagneticFieldPlausibilityCheck | None = None,
        reliability_mask_builder: ReliabilityMaskBuilder | None = None,
    ) -> None:
        self.galaxy = galaxy
        self._gap_calculator = gap_calculator or BaryonicGapCalculator()
        self._field_calibrator = field_calibrator or MagneticFieldCalibrator()
        self._gas_density = gas_density or GasDensityEstimator()
        self._pressure_gradient_estimator = pressure_gradient or MagneticPressureGradientEstimator()
        self._plasma = plasma_diagnostics or PlasmaDiagnostics()
        self._plausibility = plausibility_check or MagneticFieldPlausibilityCheck()
        self._reliability = reliability_mask_builder or ReliabilityMaskBuilder()

    def calculate_force_gap(self) -> None:
        """
        Populate the SI acceleration columns and the kinematic ``acceleration_gap``.

        Delegates to :class:`BaryonicGapCalculator`. ``V_bar`` must already be present
        (SPARC provider sets ``sqrt(V_gas² + V_disk² + V_bulge²)``). When
        ``observed_velocity_systematic_kms`` is provided by the SPARC provider, a
        ``centripetal_acceleration_systematic`` column carries the linear inclination
        propagation ``δa_c / a_c = 2 |δV/V|``.
        """
        self.galaxy.rotation_dataframe = self._gap_calculator.compute(
            self.galaxy.rotation_dataframe
        )

    def calculate_magnetic_pressure_acceleration(
        self,
        stokes_amp_center: float | None,
        stokes_amp_neighbors_mean: float | None,
        scenario_gradient_scale_kpc: float,
        scenario_ionization_fraction: float,
        stokes_amp_reference: float = DEFAULT_STOKES_AMPLITUDE_REFERENCE,
        physical_neighbor_separation_kpc: float | None = None,
        *,
        neighbour_ring: list[dict[str, float]] | None = None,
        center_theta_rad: float | None = None,
        center_phi_rad: float | None = None,
        position_angle_deg: float | None = None,
        position_angle_uncertainty_deg: float | None = None,
        morphology_class: str | None = None,
    ) -> dict[str, Any]:
        """
        Estimate the MHD magnetic-pressure acceleration ``a_B`` from Planck Stokes amplitudes
        and the per-scenario grid parameters.

        Under ideal MHD the volume Lorentz force on a magnetized fluid reduces to the negative
        gradient of the magnetic pressure ``P_B = B²/(2 μ₀)`` divided by the ionized mass density
        ``ρ_ion`` (the rest of the Lorentz term — the magnetic-tension contribution
        ``(B · ∇) B / μ₀ ρ_ion`` — is dropped in the ideal-MHD pressure-gradient branch; see
        ``docs/theory/EG_MODEL_THEORY_EN.md`` Section 5 (limitations) and Section 2--3 for the
        pressure ansatz versus tension.

            a_B(r) = − ∇ u_B(r) / ρ_ion(r),     u_B = B² / (2 μ₀).

        ``u_B`` is reconstructed from Planck polarization Stokes amplitudes ``Q,U`` calibrated to
        Tesla via ``√(2 μ₀ u_ref)`` (see ``MagneticFieldCalibrator``); ``ρ_ion`` is the gas mass
        density returned by ``GasDensityEstimator`` multiplied by the ionization fraction knob.

        Parameters
        ----------
        stokes_amp_center, stokes_amp_neighbors_mean:
            Dimensionless Stokes amplitudes at the central HEALPix pixel and (mean of) its
            neighbors.
        scenario_gradient_scale_kpc:
            Disk-scale Δx (kpc) used by the **primary** ``magnetic_pressure_acceleration`` column.
            This is the scenario-grid exploratory Δx knob.
        scenario_ionization_fraction:
            Multiplicative factor applied to the gas mass density to obtain the inertial ionized
            density used in ``a = -∇P_B / ρ_ion``.
        stokes_amp_reference:
            ``S_ref`` used by the field calibrator. Defaults to the emergency-only constant; the
            real value should always be passed by the pipeline (see provenance contract).
        physical_neighbor_separation_kpc:
            Optional **physical** projected distance (kpc) corresponding to the HEALPix neighbor
            separation, computed as ``D_galaxy · θ_neighbor`` in the pipeline. When provided, the
            engine adds a parallel ``magnetic_pressure_acceleration_physical_dx`` column that uses
            this real sky-to-disk Δx instead of the arbitrary scenario knob — letting reports
            quote both numbers and judge which Δx matters.
        neighbour_ring:
            Optional ``[{"q", "u", "theta_rad", "phi_rad", "separation_rad"}, ...]`` list
            produced by :meth:`src.providers.planck.PlanckProvider.extract_stokes_qu_mean`.
            When provided together with the galaxy distance, the engine replaces the scalar
            ``Δu_B/Δx`` gradient on the physical-Δx path with a **disk-frame planar
            least-squares** gradient (see :func:`disk_plane_pressure_gradient`).
        center_theta_rad, center_phi_rad:
            HEALPix angles of the central pixel (radians). Required when ``neighbour_ring``
            is supplied — used to convert ring offsets into a tangent-plane (Δx_east, Δy_north).
        position_angle_deg, position_angle_uncertainty_deg:
            Disk major-axis position angle (east of north, deg) and its 1-σ catalogue error.
            When PA is missing, the engine flags ``warn_no_pa`` and falls back to the
            axisymmetric sky-plane gradient magnitude.
        morphology_class:
            SPARC morphology bucket label. Required to gate the magnetic-tension toy: any
            morphology consistent with a coherent disk spiral (i.e. not Irregular/BCD) gets
            the toy logarithmic-spiral curvature term added on top of the pressure
            acceleration; irregulars get ``warn_pure_pressure`` and the tension is zero.
        """
        _ = position_angle_uncertainty_deg  # Reserved for Monte-Carlo propagation upstream.
        df = self.galaxy.rotation_dataframe.copy()
        if "radius_m" not in df.columns:
            raise ValueError(
                "Run calculate_force_gap before calculate_magnetic_pressure_acceleration."
            )

        stokes_neighbor_imputed = False
        s_center = _finite_stokes_amp(stokes_amp_center)
        s_neighbors = _finite_stokes_amp(stokes_amp_neighbors_mean)
        if s_center is not None and s_neighbors is None:
            s_neighbors = s_center
            stokes_neighbor_imputed = True

        stokes_ok = (
            s_center is not None
            and s_neighbors is not None
            and np.isfinite(scenario_gradient_scale_kpc)
            and np.isfinite(scenario_ionization_fraction)
            and scenario_gradient_scale_kpc > 0
            and scenario_ionization_fraction > 0
        )

        if not stokes_ok:
            out = df.copy()
            out["stokes_amp_center"] = np.nan
            out["stokes_amp_neighbors_mean"] = np.nan
            out["magnetic_pressure_acceleration"] = np.nan
            out["magnetic_pressure_acceleration_physical_dx"] = np.nan
            out["accel_pressure_only_m_s2"] = np.nan
            out["accel_full_lorentz_m_s2"] = np.nan
            out["tension_acceleration_m_s2"] = np.nan
            out["magnetic_field_microgauss"] = np.nan
            out["plasma_beta_thermal"] = np.nan
            out["alfven_velocity_ms"] = np.nan
            out["stellar_velocity_fraction"] = out.get("stellar_velocity_fraction", 0.0)
            out["stellar_dominated"] = False
            out["reliable_measurement"] = False
            out["density_floor_clamped"] = False
            out["n_points_density_clamped"] = 0
            out["gas_mass_residual_consistency_kg"] = np.nan
            out["scenario_gradient_scale_m"] = np.nan
            out["physical_gradient_scale_m"] = np.nan
            self.galaxy.rotation_dataframe = out
            return {
                "status": "Skipped - Missing Magnetic Data",
                "magnetic_field_tesla": None,
                "magnetic_field_microgauss": None,
                "magnetic_field_within_plausibility_band": None,
                "magnetic_pressure_acceleration": None,
                "magnetic_pressure_acceleration_physical_dx": None,
                "accel_pressure_only_m_s2": None,
                "accel_full_lorentz_m_s2": None,
                "tension_term_used": False,
                "physical_gradient_scale_kpc": physical_neighbor_separation_kpc,
                "plasma_beta_thermal_median": None,
                "plasma_beta_above_unity_fraction": None,
                "plasma_beta_below_unity_fraction": None,
                "magnetic_dominance_flag": False,
                "alfven_velocity_median_kms": None,
                "stellar_dominated_fraction": 0.0,
                "n_points_stellar_dominated": 0,
                "n_points_density_clamped": 0,
                "n_points_total": int(len(out)),
                "density_clamped_fraction": 0.0,
                "audit_flags": [],
                "deprojected_gradient_used": False,
            }

        stokes_ref = (
            float(stokes_amp_reference)
            if np.isfinite(stokes_amp_reference) and stokes_amp_reference > 0
            else DEFAULT_STOKES_AMPLITUDE_REFERENCE
        )
        magnetic_field_center_t = self._field_calibrator.field_tesla(float(s_center), stokes_ref)
        magnetic_field_neighbors_t = self._field_calibrator.field_tesla(float(s_neighbors), stokes_ref)
        magnetic_field_tesla = magnetic_field_center_t

        radius_m = df["radius_m"].to_numpy(dtype=float)
        v_disk_kms = df["disk_velocity_kms"].to_numpy(dtype=float) if "disk_velocity_kms" in df.columns else None
        v_bulge_kms = df["bulge_velocity_kms"].to_numpy(dtype=float) if "bulge_velocity_kms" in df.columns else None
        v_bar_kms = (
            df["baryonic_velocity_kms"].to_numpy(dtype=float)
            if "baryonic_velocity_kms" in df.columns
            else None
        )
        v_obs_kms = (
            df["observed_velocity_kms"].to_numpy(dtype=float)
            if "observed_velocity_kms" in df.columns
            else None
        )
        (
            gas_mass_density,
            thickness_kpc,
            density_clamped_mask,
            stellar_fraction,
            stellar_dominated_mask,
            gas_mass_residual_consistency_kg,
        ) = self._gas_density.shell_gas_mass_density(
            radius_m,
            df["gas_velocity_kms"].to_numpy(dtype=float),
            df["radius_kpc"].to_numpy(dtype=float),
            disk_velocity_kms=v_disk_kms,
            bulge_velocity_kms=v_bulge_kms,
            baryonic_velocity_kms=v_bar_kms,
            observed_velocity_kms=v_obs_kms,
        )
        n_points_density_clamped = int(np.count_nonzero(density_clamped_mask))
        n_points_total = int(density_clamped_mask.size)
        density_clamped_fraction = (
            float(n_points_density_clamped) / float(n_points_total) if n_points_total > 0 else 0.0
        )
        n_points_stellar_dominated = int(np.count_nonzero(stellar_dominated_mask))
        stellar_dominated_fraction = (
            float(n_points_stellar_dominated) / float(n_points_total) if n_points_total > 0 else 0.0
        )

        magnetic_energy_density = (magnetic_field_tesla**2) / (2.0 * MU0)
        magnetic_energy_density_center = (magnetic_field_center_t**2) / (2.0 * MU0)
        magnetic_energy_density_neighbors = (magnetic_field_neighbors_t**2) / (2.0 * MU0)
        delta_magnetic_energy_density = self._pressure_gradient_estimator.delta_energy_density(
            magnetic_energy_density_center,
            magnetic_energy_density_neighbors,
        )

        scenario_gradient_scale_m = float((scenario_gradient_scale_kpc * u.kpc).to_value(u.m))
        magnetic_pressure_gradient = self._pressure_gradient_estimator.pressure_gradient(
            delta_magnetic_energy_density,
            scenario_gradient_scale_m,
        )

        effective_ionized_mass_density = gas_mass_density * scenario_ionization_fraction
        magnetic_pressure_acceleration = self._pressure_gradient_estimator.mhd_acceleration(
            magnetic_pressure_gradient,
            effective_ionized_mass_density,
        )

        engine_audit_flags: list[str] = []
        if stokes_neighbor_imputed:
            engine_audit_flags.append(AUDIT_CODES.WARN_STOKES_NEIGHBOR_IMPUTED_FROM_CENTER)
        morph_label = "" if morphology_class is None else str(morphology_class).strip().lower()
        if morph_label == "" or morph_label == "unclassified":
            engine_audit_flags.append(AUDIT_CODES.WARN_UNKNOWN_MORPHOLOGY)

        deprojected_gradient_used = False
        if (
            physical_neighbor_separation_kpc is not None
            and np.isfinite(physical_neighbor_separation_kpc)
            and physical_neighbor_separation_kpc > 0
        ):
            physical_gradient_scale_m = float(
                (float(physical_neighbor_separation_kpc) * u.kpc).to_value(u.m)
            )
            distance_mpc_val = float(self.galaxy.distance_mpc) if self.galaxy.distance_mpc else float("nan")
            strategy = resolve_physical_gradient_strategy(
                neighbour_ring=neighbour_ring,
                center_theta_rad=center_theta_rad,
                center_phi_rad=center_phi_rad,
                distance_mpc=distance_mpc_val,
                pressure_gradient_estimator=self._pressure_gradient_estimator,
            )
            phys_in = PhysicalGradientInputs(
                delta_magnetic_energy_density_j_m3=delta_magnetic_energy_density,
                physical_gradient_scale_m=physical_gradient_scale_m,
                magnetic_energy_density_center_j_m3=float(magnetic_energy_density_center),
                stokes_amp_reference=stokes_ref,
                neighbour_ring=neighbour_ring,
                center_theta_rad=center_theta_rad,
                center_phi_rad=center_phi_rad,
                distance_mpc=distance_mpc_val,
                position_angle_deg=position_angle_deg,
                inclination_deg=self.galaxy.inclination_deg,
                field_calibrator=self._field_calibrator,
                mu0_si=MU0,
            )
            out_grad = strategy.compute(phys_in)
            physical_pressure_gradient = out_grad.physical_pressure_gradient_j_m4
            deprojected_gradient_used = out_grad.deprojected_gradient_used
            engine_audit_flags.extend(list(out_grad.audit_flags))
            magnetic_pressure_acceleration_physical_dx = self._pressure_gradient_estimator.mhd_acceleration(
                physical_pressure_gradient,
                effective_ionized_mass_density,
            )
        else:
            physical_gradient_scale_m = float("nan")
            physical_pressure_gradient = float("nan")
            magnetic_pressure_acceleration_physical_dx = np.full_like(
                magnetic_pressure_acceleration, np.nan, dtype=float
            )

        radius_m_arr = df["radius_m"].to_numpy(dtype=float)
        if morphology_class is not None and not supports_logarithmic_spiral(morphology_class):
            tension_acceleration = np.zeros_like(radius_m_arr, dtype=float)
            engine_audit_flags.append(AUDIT_CODES.WARN_PURE_PRESSURE_GEOMETRY)
        else:
            tension_acceleration = tension_acceleration_radial(
                radius_m=radius_m_arr,
                magnetic_field_tesla=magnetic_field_tesla,
                effective_ionized_mass_density=effective_ionized_mass_density,
                morphology_class=morphology_class,
            )

        accel_pressure_only = np.asarray(magnetic_pressure_acceleration, dtype=float).copy()
        accel_full_lorentz = accel_pressure_only + np.nan_to_num(tension_acceleration, nan=0.0)

        magnetic_field_microgauss = magnetic_field_tesla / TESLA_PER_MICROGAUSS
        magnetic_field_in_band = self._plausibility.is_within_band(magnetic_field_microgauss)

        thermal_pressure = self._plasma.thermal_pressure(gas_mass_density)
        plasma_beta_thermal = self._plasma.plasma_beta_thermal(
            gas_mass_density, magnetic_field_tesla
        )
        alfven_velocity_ms = self._plasma.alfven_velocity_ms(
            magnetic_field_tesla, effective_ionized_mass_density
        )

        df["stokes_amp_center"] = float(s_center)
        df["stokes_amp_neighbors_mean"] = float(s_neighbors)
        df["magnetic_field_tesla"] = magnetic_field_tesla
        df["magnetic_field_neighbors_tesla"] = magnetic_field_neighbors_t
        df["magnetic_field_microgauss"] = magnetic_field_microgauss
        df["gas_mass_density"] = gas_mass_density
        df["effective_ionized_mass_density"] = effective_ionized_mass_density
        df["thermal_pressure"] = thermal_pressure
        df["plasma_beta_thermal"] = plasma_beta_thermal
        df["alfven_velocity_ms"] = alfven_velocity_ms
        df["magnetic_energy_density"] = magnetic_energy_density
        df["magnetic_energy_density_center"] = magnetic_energy_density_center
        df["magnetic_energy_density_neighbors"] = magnetic_energy_density_neighbors
        df["delta_magnetic_energy_density"] = delta_magnetic_energy_density
        df["magnetic_pressure_gradient"] = magnetic_pressure_gradient
        df["physical_pressure_gradient"] = physical_pressure_gradient
        df["scenario_gradient_scale_kpc"] = scenario_gradient_scale_kpc
        df["scenario_gradient_scale_m"] = float(scenario_gradient_scale_m)
        df["physical_gradient_scale_kpc"] = (
            float(physical_neighbor_separation_kpc) if physical_neighbor_separation_kpc is not None else np.nan
        )
        df["physical_gradient_scale_m"] = float(physical_gradient_scale_m)
        df["scenario_ionization_fraction"] = scenario_ionization_fraction
        df["magnetic_pressure_acceleration"] = magnetic_pressure_acceleration
        df["magnetic_pressure_acceleration_physical_dx"] = magnetic_pressure_acceleration_physical_dx
        df["accel_pressure_only_m_s2"] = accel_pressure_only
        df["accel_full_lorentz_m_s2"] = accel_full_lorentz
        df["tension_acceleration_m_s2"] = tension_acceleration
        df["density_floor_clamped"] = density_clamped_mask.astype(bool)
        df["n_points_density_clamped"] = n_points_density_clamped
        df["stellar_velocity_fraction"] = stellar_fraction
        df["stellar_dominated"] = stellar_dominated_mask.astype(bool)
        df["gas_mass_residual_consistency_kg"] = gas_mass_residual_consistency_kg

        within_acceleration_envelope = self._reliability.build_mask(
            df["magnetic_pressure_acceleration"].to_numpy(dtype=float),
            df["centripetal_acceleration"].to_numpy(dtype=float),
            density_clamped_mask,
            stellar_dominated_mask,
        )
        df["reliable_measurement"] = (
            within_acceleration_envelope & ~density_clamped_mask & ~stellar_dominated_mask
        )

        gap_abs = np.abs(df["acceleration_gap"].to_numpy(dtype=float))
        use_physical = (
            physical_neighbor_separation_kpc is not None
            and np.isfinite(physical_gradient_scale_m)
            and physical_gradient_scale_m > 0
        )
        if use_physical:
            req_gradient = np.abs(physical_pressure_gradient)
            req_acceleration_abs = np.abs(
                df["magnetic_pressure_acceleration_physical_dx"].to_numpy(dtype=float)
            )
        else:
            req_gradient = np.abs(magnetic_pressure_gradient)
            req_acceleration_abs = np.abs(
                df["magnetic_pressure_acceleration"].to_numpy(dtype=float)
            )
        required_density_for_match_point = np.where(
            gap_abs > 0,
            req_gradient / gap_abs,
            np.nan,
        )
        underperforming = req_acceleration_abs < gap_abs
        req_valid = required_density_for_match_point[
            np.isfinite(required_density_for_match_point) & underperforming & (required_density_for_match_point > 0)
        ]
        required_density_for_match = float(np.nanmedian(req_valid)) if req_valid.size > 0 else None

        df["required_density_for_match"] = required_density_for_match_point
        self.galaxy.rotation_dataframe = df

        unreliable_reason = (
            "Sky-map neighbor scale mixed with disk gradient scale (Δx); interpret as exploratory. "
            "Also screen for foreground / noise if |a_B| ≫ |a_centripetal|."
        )

        beta_finite = plasma_beta_thermal[np.isfinite(plasma_beta_thermal)]
        if beta_finite.size > 0:
            plasma_beta_above_unity_fraction = float(
                np.count_nonzero(beta_finite > RELIABILITY_CUTS.MAX_PLASMA_BETA_FOR_MHD_DOMINANCE)
            ) / float(beta_finite.size)
            plasma_beta_below_unity_fraction = float(
                np.count_nonzero(beta_finite < 1.0)
            ) / float(beta_finite.size)
            plasma_beta_median_value = float(np.median(beta_finite))
        else:
            plasma_beta_above_unity_fraction = None
            plasma_beta_below_unity_fraction = None
            plasma_beta_median_value = None
        magnetic_dominance_flag = (
            plasma_beta_median_value is not None and plasma_beta_median_value < 1.0
        )

        return {
            "status": "ok",
            "magnetic_field_tesla": magnetic_field_tesla,
            "magnetic_field_neighbors_tesla": magnetic_field_neighbors_t,
            "magnetic_field_microgauss": float(magnetic_field_microgauss),
            "magnetic_field_within_plausibility_band": bool(magnetic_field_in_band),
            "magnetic_energy_density": magnetic_energy_density,
            "delta_magnetic_energy_density": delta_magnetic_energy_density,
            "magnetic_pressure_gradient": magnetic_pressure_gradient,
            "scenario_gradient_scale_kpc": scenario_gradient_scale_kpc,
            "scenario_gradient_scale_m": float(scenario_gradient_scale_m),
            "physical_gradient_scale_kpc": (
                float(physical_neighbor_separation_kpc)
                if physical_neighbor_separation_kpc is not None
                else None
            ),
            "physical_gradient_scale_m": (
                float(physical_gradient_scale_m)
                if np.isfinite(physical_gradient_scale_m)
                else None
            ),
            "scenario_ionization_fraction": scenario_ionization_fraction,
            "stokes_amplitude_local": float(s_center),
            "stokes_amplitude_neighbors": float(s_neighbors),
            "stokes_amplitude_reference": stokes_ref,
            "required_density_for_match": required_density_for_match,
            "disk_thickness_kpc": thickness_kpc,
            "unreliable_points_count": int(np.count_nonzero(~df["reliable_measurement"])),
            "unreliable_reason": unreliable_reason,
            "magnetic_pressure_acceleration": float(np.nanmedian(magnetic_pressure_acceleration)),
            "magnetic_pressure_acceleration_physical_dx": (
                float(np.nanmedian(magnetic_pressure_acceleration_physical_dx))
                if np.any(np.isfinite(magnetic_pressure_acceleration_physical_dx))
                else None
            ),
            "accel_pressure_only_m_s2": (
                float(np.nanmedian(accel_pressure_only))
                if np.any(np.isfinite(accel_pressure_only))
                else None
            ),
            "accel_full_lorentz_m_s2": (
                float(np.nanmedian(accel_full_lorentz))
                if np.any(np.isfinite(accel_full_lorentz))
                else None
            ),
            "tension_acceleration_median_m_s2": (
                float(np.nanmedian(tension_acceleration))
                if np.any(np.isfinite(tension_acceleration))
                else None
            ),
            "tension_term_used": bool(
                np.any(np.isfinite(tension_acceleration))
                and np.any(np.asarray(tension_acceleration) != 0.0)
            ),
            "deprojected_gradient_used": bool(deprojected_gradient_used),
            "audit_flags": list(engine_audit_flags),
            "plasma_beta_thermal_median": plasma_beta_median_value,
            "plasma_beta_above_unity_fraction": plasma_beta_above_unity_fraction,
            "plasma_beta_below_unity_fraction": plasma_beta_below_unity_fraction,
            "magnetic_dominance_flag": bool(magnetic_dominance_flag),
            "alfven_velocity_median_kms": (
                float(np.nanmedian(alfven_velocity_ms) / 1000.0)
                if np.any(np.isfinite(alfven_velocity_ms))
                else None
            ),
            "stellar_dominated_fraction": stellar_dominated_fraction,
            "n_points_stellar_dominated": n_points_stellar_dominated,
            "n_points_density_clamped": n_points_density_clamped,
            "n_points_total": n_points_total,
            "density_clamped_fraction": density_clamped_fraction,
        }
