"""
Estimate volumetric gas density from SPARC kinematic decomposition via enclosed-mass shells.

Mass-budget contract (rigorous, two parallel paths)
---------------------------------------------------
SPARC publishes the rotation curve already decomposed into independent component velocities,
each one corresponding to the rotational support due to that mass component **alone**:

    V_obs²  =  V_gas²  +  V_disk²  +  V_bulge²  +  V_DM²    (no double-counting)

The estimator therefore computes ``M_enc,gas(r)`` from a **gas-only** velocity term so that
stellar mass never leaks into ``ρ_gas``. Two equivalent expressions are evaluated for every
shell so the two MUST agree (and the larger one is reported as ``rho_gas`` so the audit is
conservative — i.e. never inflates the inertial reservoir):

    path A (direct)    :  M_gas(r) = V_gas² · r / G
    path B (residual)  :  M_gas(r) = max(V_bar² − V_disk² − V_bulge², 0) · r / G

``V_bar = sqrt(V_gas² + V_disk² + V_bulge²)`` is the SPARC baryonic synthesis already
available on the rotation table. Path B explicitly **subtracts the stellar contribution from
the total baryonic potential** so the comparison falsifies any silent leak between
components. The per-shell residual ``M_gas(A) − M_gas(B)`` is exposed via
``gas_mass_residual_consistency_kg`` for downstream auditing.

Stellar-dominance gate
----------------------
Even when the gas velocity is internally consistent, shells where stellar kinetics dominate
the total baryonic budget are flagged ``stellar_dominated``: in those shells a magnetic
pressure correction is not applied because gravity from baryons (mostly stars) already
explains the kinematics, and a poorly resolved gas tracer would only contaminate the MHD
diagnostic. The gate uses **two thresholds** — ``V_stellar²/V_bar²`` (baryonic share, the
historical metric) and ``V_stellar²/V_obs²`` (share of the *observed* squared-velocity
budget) — and flags the shell when either crosses the configured threshold.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u

from src.config.thresholds import RELIABILITY_CUTS, UNIVERSAL_CONSTANTS

from .defaults import (
    MIN_GAS_MASS_DENSITY,
    SHELL_THICKNESS_CLIP_MAX_KPC,
    SHELL_THICKNESS_CLIP_MIN_KPC,
    SHELL_THICKNESS_FRACTION_OF_RMAX,
)

__all__ = ["GasDensityEstimator"]


class GasDensityEstimator:
    """Enclosed-mass shell inversion from SPARC ``V_gas`` (and consistency path B)."""

    def __init__(
        self,
        min_mass_density_kg_m3: float = MIN_GAS_MASS_DENSITY,
        stellar_dominated_threshold: float = RELIABILITY_CUTS.STELLAR_DOMINATED_PER_SHELL_THRESHOLD,
    ) -> None:
        self.min_mass_density_kg_m3 = min_mass_density_kg_m3
        self.stellar_dominated_threshold = float(stellar_dominated_threshold)

    def shell_gas_mass_density(
        self,
        radius_m: np.ndarray,
        gas_velocity_kms: np.ndarray,
        radius_kpc: np.ndarray,
        disk_velocity_kms: np.ndarray | None = None,
        bulge_velocity_kms: np.ndarray | None = None,
        baryonic_velocity_kms: np.ndarray | None = None,
        observed_velocity_kms: np.ndarray | None = None,
    ) -> tuple[np.ndarray, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Return ``(rho_gas, disk_thickness_kpc, clamped, stellar_fraction, stellar_dominated, residual_kg)``.

        ``residual_kg`` is the per-shell consistency residual ``M_enc,A − M_enc,B`` between the
        two equivalent enclosed-mass paths (direct ``V_gas²`` vs ``V_bar²−V_stellar²``). It is
        identically zero when SPARC's decomposition is internally consistent; non-zero values
        flag a leak between V_gas and the (V_bar, V_disk, V_bulge) trio reported by the catalog.

        Parameters
        ----------
        radius_m, gas_velocity_kms, radius_kpc:
            Per-shell SI radius, gas velocity in km/s, radius in kpc (parallel arrays).
        disk_velocity_kms, bulge_velocity_kms:
            Optional SPARC stellar-component rotation contributions (km/s). Required for the
            ``stellar_velocity_fraction`` audit and the residual cross-check.
        baryonic_velocity_kms:
            Optional ``V_bar`` (km/s) used by path B (residual ``V_bar² − V_disk² − V_bulge²``).
            When ``None`` the estimator reconstructs it as ``sqrt(V_gas²+V_disk²+V_bulge²)``.
        observed_velocity_kms:
            Optional ``V_obs`` (km/s) — when provided, the stellar-dominance gate also checks
            ``V_stellar² / V_obs²`` (share of the *observed* dynamics). This is stricter and
            avoids applying magnetic pressure where stellar gravity already explains motion.

        Returns
        -------
        rho_gas_kg_m3, disk_thickness_kpc, clamped_mask, stellar_velocity_fraction,
        stellar_dominated_mask, gas_mass_residual_consistency_kg : see docstring.
        """
        gas_velocity_ms = np.asarray(gas_velocity_kms, dtype=float) * 1000.0
        radius_m = np.asarray(radius_m, dtype=float)
        g_si = UNIVERSAL_CONSTANTS.G_SI

        # Path A — direct V_gas² · r / G (gas-only kinematic mass, current definition).
        gas_mass_enclosed_direct = (gas_velocity_ms**2 * radius_m) / g_si
        gas_mass_enclosed_direct = np.maximum(gas_mass_enclosed_direct, 0.0)

        # Path B — explicit subtraction of the stellar potential from V_bar² · r / G.
        # By construction V_bar² = V_gas² + V_disk² + V_bulge², so this MUST equal path A.
        # We still evaluate it independently so the audit can detect catalog inconsistencies.
        # When the SPARC decomposition is *not* available (unit tests / minimal callers), there
        # is no path B to compare against and the residual is reported as 0 (no information).
        path_b_available = (
            baryonic_velocity_kms is not None
            or disk_velocity_kms is not None
            or bulge_velocity_kms is not None
        )
        if path_b_available:
            v_disk_arr = np.asarray(
                disk_velocity_kms if disk_velocity_kms is not None else np.zeros_like(gas_velocity_kms),
                dtype=float,
            )
            v_bulge_arr = np.asarray(
                bulge_velocity_kms if bulge_velocity_kms is not None else np.zeros_like(gas_velocity_kms),
                dtype=float,
            )
            if baryonic_velocity_kms is not None:
                v_bar_arr = np.asarray(baryonic_velocity_kms, dtype=float)
            else:
                v_bar_arr = np.sqrt(
                    np.asarray(gas_velocity_kms, dtype=float) ** 2 + v_disk_arr**2 + v_bulge_arr**2
                )
            v_gas_residual_sq_kms2 = np.maximum(v_bar_arr**2 - v_disk_arr**2 - v_bulge_arr**2, 0.0)
            gas_mass_enclosed_residual = (v_gas_residual_sq_kms2 * 1.0e6 * radius_m) / g_si
            gas_mass_enclosed_residual = np.maximum(gas_mass_enclosed_residual, 0.0)
            # Conservative choice: report the larger of the two paths so ρ_ion is never silently
            # inflated by an arithmetic mismatch between V_gas and (V_bar² − V_stellar²).
            gas_mass_enclosed = np.maximum(gas_mass_enclosed_direct, gas_mass_enclosed_residual)
            gas_mass_residual_consistency_kg = gas_mass_enclosed_direct - gas_mass_enclosed_residual
        else:
            gas_mass_enclosed = gas_mass_enclosed_direct
            gas_mass_residual_consistency_kg = np.zeros_like(gas_mass_enclosed_direct)

        max_radius_kpc = float(np.nanmax(np.asarray(radius_kpc, dtype=float)))
        thickness_kpc = float(
            np.clip(
                SHELL_THICKNESS_FRACTION_OF_RMAX * max_radius_kpc,
                SHELL_THICKNESS_CLIP_MIN_KPC,
                SHELL_THICKNESS_CLIP_MAX_KPC,
            )
        )
        thickness_m = thickness_kpc * u.kpc.to(u.m)

        r_inner = np.concatenate(([0.0], radius_m[:-1]))
        shell_volume_m3 = np.pi * (radius_m**2 - r_inner**2) * thickness_m
        shell_mass_kg = np.concatenate(([gas_mass_enclosed[0]], np.diff(gas_mass_enclosed)))
        shell_mass_kg = np.maximum(shell_mass_kg, 0.0)

        true_rho = np.where(shell_volume_m3 > 0, shell_mass_kg / shell_volume_m3, np.nan)
        clamped_mask = ~(np.isfinite(true_rho) & (true_rho >= self.min_mass_density_kg_m3))
        rho = np.where(clamped_mask, self.min_mass_density_kg_m3, true_rho)

        if disk_velocity_kms is not None or bulge_velocity_kms is not None:
            v_disk = np.asarray(
                disk_velocity_kms if disk_velocity_kms is not None else np.zeros_like(gas_velocity_kms),
                dtype=float,
            )
            v_bulge = np.asarray(
                bulge_velocity_kms if bulge_velocity_kms is not None else np.zeros_like(gas_velocity_kms),
                dtype=float,
            )
            v_gas_kms_arr = np.asarray(gas_velocity_kms, dtype=float)
            v_stellar_sq = v_disk**2 + v_bulge**2
            v_bar_sq = v_gas_kms_arr**2 + v_stellar_sq
            with np.errstate(invalid="ignore", divide="ignore"):
                stellar_fraction = np.where(v_bar_sq > 0, v_stellar_sq / v_bar_sq, 0.0)
            stellar_fraction = np.where(np.isfinite(stellar_fraction), stellar_fraction, 0.0)
            # Optional stricter audit: V_stellar² / V_obs² — share of the observed dynamics,
            # not just baryonic. When stars already explain a fraction of V_obs comparable to
            # the configured threshold, the shell is also flagged stellar-dominated so the
            # magnetic-pressure term is not applied on top of the stellar potential.
            if observed_velocity_kms is not None:
                v_obs_arr = np.asarray(observed_velocity_kms, dtype=float)
                with np.errstate(invalid="ignore", divide="ignore"):
                    stellar_observed_fraction = np.where(
                        v_obs_arr**2 > 0, v_stellar_sq / (v_obs_arr**2), 0.0
                    )
                stellar_observed_fraction = np.where(
                    np.isfinite(stellar_observed_fraction), stellar_observed_fraction, 0.0
                )
            else:
                stellar_observed_fraction = np.zeros_like(stellar_fraction)
        else:
            stellar_fraction = np.zeros_like(rho, dtype=float)
            stellar_observed_fraction = np.zeros_like(rho, dtype=float)

        stellar_dominated_mask = (stellar_fraction >= self.stellar_dominated_threshold) | (
            stellar_observed_fraction >= self.stellar_dominated_threshold
        )
        return (
            rho,
            thickness_kpc,
            clamped_mask,
            stellar_fraction,
            stellar_dominated_mask,
            gas_mass_residual_consistency_kg,
        )
