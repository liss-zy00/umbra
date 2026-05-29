"""
Rotation-curve kinematics: centripetal vs baryonic acceleration and the gap.

Encapsulates the **first half** of :class:`src.physics.engine.PhysicsEngine` — the part that
turns SPARC's published velocities into SI accelerations and computes the kinematic
acceleration deficit conventionally attributed to dark matter.

Physical content
----------------
For a circular orbit of radius ``r`` the centripetal acceleration is
``a_c(r) = V_obs²(r) / r``. SPARC publishes the rotation curve already decomposed into
non-overlapping component velocities, so the *Newtonian baryonic prediction* is

    a_bar(r) = V_bar²(r) / r,    V_bar = √(V_gas² + V_disk² + V_bulge²).

The **acceleration gap** is the residual ``a_gap(r) = a_c(r) − a_bar(r)``. In the
standard ΛCDM picture the gap is fitted with a dark-matter halo. The E-G hypothesis
asks whether the gap can instead be supplied by an MHD magnetic-pressure term
``a_B(r) = − ∇ u_B(r) / ρ_ion(r)``; the ``PhysicsEngine`` evaluates the latter and the
metrics layer compares them.

Inclination uncertainty (linear propagation)
--------------------------------------------
SPARC distributes ``V_obs = V_LOS / sin i`` already deprojected. The inclination has its
own catalog uncertainty ``δi``; the SparcProvider writes a per-shell systematic
``δV_obs = V_obs · |cos i / sin i| · δi``. This module propagates it onto the centripetal
acceleration via ``δa_c / a_c = 2 |δV/V|``, exposing
``centripetal_acceleration_systematic`` for the report layer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from astropy import units as u

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS

from .diagnostics.rar import rar_predicted_observed_acceleration, rar_residual_log10

__all__ = ["BaryonicGapCalculator"]


@dataclass(frozen=True)
class BaryonicGapCalculator:
    """
    Convert SPARC rotation-curve velocities (km/s) to SI accelerations and compute the gap.

    The class is stateless aside from configuration knobs (currently none with safe defaults
    other than the inclination-systematic factor). Inputs are read from the provided
    ``rotation_dataframe`` columns; outputs are returned as a *new* dataframe so the
    estimator never mutates its input — this preserves SRP and lets the engine decide when
    to commit results back to the ``Galaxy`` model.
    """

    inclination_systematic_factor: float = (
        ASTROPHYSICAL_DEFAULTS.INCLINATION_SYSTEMATIC_PROPAGATION_FACTOR
    )

    REQUIRED_COLUMNS: tuple[str, ...] = (
        "radius_kpc",
        "observed_velocity_kms",
        "baryonic_velocity_kms",
    )

    def _ensure_columns(self, df: pd.DataFrame) -> None:
        missing = set(self.REQUIRED_COLUMNS) - set(df.columns)
        if missing:
            raise ValueError(f"Missing columns for force gap: {sorted(missing)}")

    @staticmethod
    def _validate_si_inputs(radius_m: np.ndarray, velocity_ms: np.ndarray) -> None:
        if np.any(~np.isfinite(radius_m)) or np.any(radius_m <= 0):
            raise ValueError("R must be finite and in meters.")
        if np.any(~np.isfinite(velocity_ms)):
            raise ValueError("V must be finite and in m/s.")

    def compute(self, rotation_dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Return a copy of the rotation dataframe enriched with SI acceleration columns:

        * ``radius_m``, ``observed_velocity_ms``, ``baryonic_velocity_ms``,
        * ``centripetal_acceleration``, ``baryonic_acceleration``, ``acceleration_gap``,
        * ``rar_predicted_acceleration`` and ``rar_residual_log10`` (McGaugh+2016 fit),
        * ``centripetal_acceleration_systematic`` (linear propagation of the SPARC
          inclination uncertainty when ``observed_velocity_systematic_kms`` is present).

        Raises
        ------
        ValueError
            If required SPARC columns are missing or if the SI conversions yield
            non-finite / non-positive radii or non-finite velocities.
        """
        self._ensure_columns(rotation_dataframe)
        df = rotation_dataframe.copy()

        radius_m = df["radius_kpc"].to_numpy(dtype=float) * u.kpc.to(u.m)
        observed_velocity_ms = df["observed_velocity_kms"].to_numpy(dtype=float) * 1000.0
        baryonic_velocity_ms = df["baryonic_velocity_kms"].to_numpy(dtype=float) * 1000.0

        self._validate_si_inputs(radius_m, observed_velocity_ms)
        self._validate_si_inputs(radius_m, baryonic_velocity_ms)

        df["radius_m"] = radius_m
        df["observed_velocity_ms"] = observed_velocity_ms
        df["baryonic_velocity_ms"] = baryonic_velocity_ms
        df["centripetal_acceleration"] = (observed_velocity_ms**2) / radius_m
        df["baryonic_acceleration"] = (baryonic_velocity_ms**2) / radius_m
        df["acceleration_gap"] = df["centripetal_acceleration"] - df["baryonic_acceleration"]

        # RAR predicted centripetal (g_obs as a function of g_bar, McGaugh+2016) and the
        # per-radius log-residual log10(g_obs_measured / g_obs_RAR). These columns let the
        # report compare each rotation curve to the empirical RAR.
        baryonic_acc = df["baryonic_acceleration"].to_numpy(dtype=float)
        observed_acc = df["centripetal_acceleration"].to_numpy(dtype=float)
        df["rar_predicted_acceleration"] = rar_predicted_observed_acceleration(baryonic_acc)
        df["rar_residual_log10"] = rar_residual_log10(observed_acc, baryonic_acc)

        # Linear inclination-systematic propagation: δa_c / a_c = 2 · |δV / V|.
        # When the upstream provider could not propagate the inclination error the column
        # is zero — *not* missing — so downstream code can broadcast it as an additive band.
        if "observed_velocity_systematic_kms" in df.columns:
            v_sys_ms = df["observed_velocity_systematic_kms"].to_numpy(dtype=float) * 1000.0
            with np.errstate(invalid="ignore", divide="ignore"):
                rel = np.where(
                    observed_velocity_ms != 0.0,
                    np.abs(v_sys_ms / observed_velocity_ms),
                    0.0,
                )
            df["centripetal_acceleration_systematic"] = (
                self.inclination_systematic_factor * np.abs(df["centripetal_acceleration"]) * rel
            )
        else:
            df["centripetal_acceleration_systematic"] = 0.0

        return df
