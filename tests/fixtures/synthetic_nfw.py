"""
Synthetic NFW rotation-curve generator used as a *control* for the magnetic-pressure model.

A canonical NFW halo with concentration ``c`` and virial velocity ``V_200`` produces a rotation
curve that rises, peaks, and asymptotes — the classic dark-matter signature. We simulate a galaxy
where the **observed** velocity comes entirely from this NFW halo (plus a tiny baryonic disk so
``V_bar`` is non-degenerate), then run it through the pipeline.

Expected pipeline behavior for this fixture:
    * ``acceleration_gap`` is large at outer radii (the canonical DM signature).
    * With **uniform** Stokes amplitudes (no real gradient), ``magnetic_pressure_acceleration``
      is essentially zero — the MHD-pressure model cannot fake an NFW residual when the input
      polarization map has no spatial structure. ``theoretical_gap_coverage_ratio`` should
      therefore stay well below the success band (``SCENARIO_SUCCESS_RATIO_MIN``).

This is the *negative control* / falsifiability check that the rest of the test suite needs.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def nfw_velocity_kms(
    radius_kpc: np.ndarray,
    v200_kms: float = 100.0,
    concentration: float = 10.0,
    r200_kpc: float = 100.0,
) -> np.ndarray:
    """
    Circular velocity of an NFW halo at ``radius_kpc``.

    Uses the standard expression :math:`V^2(r) = V_{200}^2 \\, \\frac{1}{x} \\, \\frac{\\ln(1+cx) - cx/(1+cx)}{\\ln(1+c) - c/(1+c)}`,
    where ``x = r / R_200`` and ``c`` is the concentration.
    """
    x = np.asarray(radius_kpc, dtype=float) / float(r200_kpc)
    c = float(concentration)
    g_c = np.log(1.0 + c) - c / (1.0 + c)
    cx = c * x
    g_cx = np.log(1.0 + cx) - cx / (1.0 + cx)
    v_sq_over_v200_sq = g_cx / (g_c * np.where(x > 0, x, 1.0))
    return float(v200_kms) * np.sqrt(np.where(x > 0, v_sq_over_v200_sq, 0.0))


def write_nfw_rotmod(
    path: Path,
    *,
    galaxy_name: str = "SYN_NFW",
    distance_mpc: float = 10.0,
    radii_kpc: tuple[float, ...] | None = None,
    v200_kms: float = 100.0,
    concentration: float = 10.0,
    r200_kpc: float = 100.0,
    v_disk_kms: float = 5.0,
    v_gas_kms: float = 5.0,
) -> Path:
    """
    Write a SPARC-format ``*_rotmod.dat`` file whose ``V_obs`` is purely the NFW circular velocity
    (plus a tiny baryonic disk + gas component so V_bar is non-zero).

    The file follows the canonical SPARC layout:
        ``# Distance = D Mpc`` header line, then per-radius rows with columns
        ``radius observed_velocity err_v gas disk bulge``.
    """
    if radii_kpc is None:
        radii_kpc = tuple(float(x) for x in np.linspace(1.0, 20.0, 20))
    radius_arr = np.asarray(radii_kpc, dtype=float)
    v_nfw = nfw_velocity_kms(
        radius_arr, v200_kms=v200_kms, concentration=concentration, r200_kpc=r200_kpc
    )
    # Total V_obs reflects "pure DM": V_obs² ≈ V_NFW² + V_baryonic² (disk + gas).
    v_disk = np.full_like(radius_arr, float(v_disk_kms))
    v_gas = np.full_like(radius_arr, float(v_gas_kms))
    v_bulge = np.zeros_like(radius_arr)
    v_obs = np.sqrt(v_nfw**2 + v_disk**2 + v_gas**2 + v_bulge**2)
    err_v = np.full_like(radius_arr, 1.0)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"# Distance = {float(distance_mpc):.2f} Mpc\n")
        handle.write(f"# Synthetic NFW-only fixture for galaxy {galaxy_name}\n")
        handle.write("# radius_kpc V_obs err V_gas V_disk V_bulge\n")
        for i in range(radius_arr.size):
            handle.write(
                f"{radius_arr[i]:.6f} {v_obs[i]:.6f} {err_v[i]:.6f} "
                f"{v_gas[i]:.6f} {v_disk[i]:.6f} {v_bulge[i]:.6f}\n"
            )
    return path
