"""
Sensitivity of ``theoretical_gap_coverage_ratio`` to the magnetic-field calibration anchor.

The pipeline calibrates Stokes amplitudes to Tesla through
:class:`src.physics.engine.MagneticFieldCalibrator`, which uses the relation

.. math::

    B_{\\mathrm{scale}} \\;=\\; \\sqrt{2\\,\\mu_0\\,u_{\\mathrm{ref}}},

where ``u_ref`` is a dimensional anchor only — the physical claim is *not* that the
Galactic field originates from the CMB; the constant simply fixes the units. This module
swaps the anchor between two physically motivated choices and reports the per-galaxy
coverage ratio under each, so reports can quantify how much of the headline number is
intrinsic vs. anchor-dependent.

Two anchors are supported today:

* ``"cmb"`` — CMB energy density ``u_CMB ≈ 4.17 × 10⁻¹⁴ J/m³``, the current default
  (:attr:`src.config.thresholds.AstrophysicalDefaults.CMB_ENERGY_DENSITY_J_M3`).
* ``"turbulent_equipartition"`` — warm-neutral-medium turbulent pressure floor
  ``0.5 · ρ_gas · σ_v²`` with ``σ_v = 10 km/s`` and ``ρ_gas ≈ 10⁻²⁴ kg/m³``
  (:attr:`src.config.thresholds.AstrophysicalDefaults.TURBULENT_EQUIPARTITION_J_M3`).

A larger anchor inflates ``B_scale`` ⇒ inflates ``u_B = B²/(2μ₀)`` ⇒ inflates the
magnetic-pressure gradient ⇒ inflates the coverage ratio. Reporting the log10 delta
makes the multiplicative scaling obvious.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS


ANCHORS: Mapping[str, float] = {
    "cmb": float(ASTROPHYSICAL_DEFAULTS.CMB_ENERGY_DENSITY_J_M3),
    "turbulent_equipartition": float(ASTROPHYSICAL_DEFAULTS.TURBULENT_EQUIPARTITION_J_M3),
}


def anchor_scale_factor(anchor: str) -> float:
    """
    Return the multiplicative factor ``√(u_anchor / u_CMB)`` applied to ``B_scale``.

    Because the calibration is linear in ``B_scale`` and the magnetic-pressure gradient
    scales as ``B² ∝ B_scale²``, the per-galaxy ``theoretical_gap_coverage_ratio`` under
    a non-CMB anchor is the CMB value times ``(u_anchor / u_CMB)``. The helper exposes
    the multiplier so callers can apply it directly to the engine output without
    re-instantiating the calibrator.
    """
    u_ref = ANCHORS.get(str(anchor))
    if u_ref is None or not np.isfinite(u_ref) or u_ref <= 0.0:
        raise ValueError(f"Unknown anchor key: {anchor!r}")
    return float(u_ref / ANCHORS["cmb"])


def sweep_coverage_ratio(
    *,
    galaxy_name: str,
    base_ratio_cmb: float,
) -> dict[str, Any]:
    """
    Recompute ``theoretical_gap_coverage_ratio`` under each anchor in :data:`ANCHORS`.

    Parameters
    ----------
    galaxy_name:
        Name of the galaxy — passed through verbatim so the writer can identify the row.
    base_ratio_cmb:
        The engine-computed ratio under the default CMB anchor. Non-finite values
        propagate to NaN columns; the row is still emitted so the CSV's coverage
        inventory is preserved.

    Returns
    -------
    dict[str, Any]
        Row dict with keys ``galaxy_name``, ``ratio_cmb``, ``ratio_turbulent_equipartition``,
        ``ratio_log10_delta`` (= ``log10(ratio_turbulent_equipartition / ratio_cmb)``).
    """
    base = float(base_ratio_cmb) if base_ratio_cmb is not None else float("nan")

    row: dict[str, Any] = {"galaxy_name": str(galaxy_name)}
    row["ratio_cmb"] = base
    finite = np.isfinite(base) and base > 0.0

    for key in ANCHORS:
        if key == "cmb":
            continue
        if not finite:
            row[f"ratio_{key}"] = float("nan")
            continue
        try:
            factor = anchor_scale_factor(key)
        except ValueError:
            row[f"ratio_{key}"] = float("nan")
            continue
        row[f"ratio_{key}"] = float(base * factor)

    other = row.get("ratio_turbulent_equipartition", float("nan"))
    if finite and np.isfinite(other) and other > 0.0:
        row["ratio_log10_delta"] = float(np.log10(other / base))
    else:
        row["ratio_log10_delta"] = float("nan")
    return row


def build_anchor_sensitivity_rows(
    galaxy_to_ratio: Mapping[str, float],
) -> list[dict[str, Any]]:
    """
    Compose the per-galaxy anchor-sensitivity rows from a ``{name -> ratio_cmb}`` map.

    The output is stable-ordered by galaxy name so the CSV is byte-deterministic across
    re-runs on the same input set. NaN ratios are preserved so the inventory matches the
    number of galaxies fed to the pipeline.
    """
    return [
        sweep_coverage_ratio(galaxy_name=str(name), base_ratio_cmb=float(ratio) if ratio is not None else float("nan"))
        for name, ratio in sorted(galaxy_to_ratio.items(), key=lambda kv: str(kv[0]))
    ]


__all__ = [
    "ANCHORS",
    "anchor_scale_factor",
    "sweep_coverage_ratio",
    "build_anchor_sensitivity_rows",
]
