"""
Beam-smearing falsifiability plot: ``convergence_ratio_to_wim`` vs ``distance_mpc``.

Reads a Discovery-phase CSV (e.g. ``results/step_discovery/discovery_reliable.csv``),
fits a log–log relation, and writes a PNG. Invoked only via import (no standalone CLI).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import pandas as pd

from src.config.thresholds import PLANCK_CALIBRATION
from src.physics import beam_physical_size_kpc


def _load_reliable_subset(input_path: Path) -> pd.DataFrame:
    if not input_path.exists():
        raise FileNotFoundError(
            f"Required Discovery CSV not found: {input_path}. Run the pipeline first."
        )
    df = pd.read_csv(input_path)
    for col in ("distance_mpc", "convergence_ratio_to_wim", "required_density_for_match"):
        if col not in df.columns:
            raise KeyError(f"{input_path} is missing the required column '{col}'.")
    df = df.copy()
    df["distance_mpc"] = pd.to_numeric(df["distance_mpc"], errors="coerce")
    df["convergence_ratio_to_wim"] = pd.to_numeric(
        df["convergence_ratio_to_wim"], errors="coerce"
    )
    df["required_density_for_match"] = pd.to_numeric(
        df["required_density_for_match"], errors="coerce"
    )
    finite_mask = (
        np.isfinite(df["distance_mpc"])
        & (df["distance_mpc"] > 0)
        & np.isfinite(df["convergence_ratio_to_wim"])
        & (df["convergence_ratio_to_wim"] > 0)
    )
    return df.loc[finite_mask].reset_index(drop=True)


def _powerlaw_fit_loglog(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """
    Return ``(slope, intercept_at_x=1, pearson_r)`` of the OLS fit
    ``log10(y) = slope · log10(x) + intercept``. ``NaN`` is returned when the sample
    size is below two valid points.
    """
    log_x = np.log10(x)
    log_y = np.log10(y)
    valid = np.isfinite(log_x) & np.isfinite(log_y)
    if int(np.count_nonzero(valid)) < 2:
        return float("nan"), float("nan"), float("nan")
    log_x = log_x[valid]
    log_y = log_y[valid]
    slope, intercept = np.polyfit(log_x, log_y, deg=1)
    pearson = float(np.corrcoef(log_x, log_y)[0, 1])
    return float(slope), float(intercept), pearson


def _annotate(ax: plt.Axes, slope: float, pearson: float, beam_arcmin: float) -> None:
    annotation = (
        f"Power-law fit log–log:\n"
        f"  slope = {slope:.3f}\n"
        f"  Pearson r = {pearson:.3f}\n"
        f"Beam-smearing prediction:\n"
        f"  slope_expected ≈ +1.0 (ρ_req ∝ D, since D_beam ∝ D)\n"
        f"Planck beam FWHM: {beam_arcmin:.1f}'"
    )
    ax.text(
        0.02,
        0.98,
        annotation,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.85, edgecolor="0.6"),
    )


def plot_beam_smearing(input_path: Path, output_path: Path) -> Path:
    """
    Build the beam-smearing diagnostic figure and write ``output_path``.

    Returns the resolved output path.
    """
    df = _load_reliable_subset(input_path)
    if df.empty:
        raise RuntimeError(
            f"No reliable galaxies with finite (distance, convergence_ratio_to_wim) in {input_path}."
        )

    distance_mpc = df["distance_mpc"].to_numpy(dtype=float)
    conv_ratio = df["convergence_ratio_to_wim"].to_numpy(dtype=float)

    slope, _intercept, pearson = _powerlaw_fit_loglog(distance_mpc, conv_ratio)

    beam_arcmin = float(PLANCK_CALIBRATION.BEAM_FWHM_ARCMIN)
    distance_grid_mpc = np.linspace(distance_mpc.min(), distance_mpc.max(), 64)
    beam_grid_kpc = np.array(
        [beam_physical_size_kpc(distance_mpc=float(d)) for d in distance_grid_mpc], dtype=float
    )

    fig, (ax_main, ax_beam) = plt.subplots(
        2,
        1,
        figsize=(9.0, 7.5),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    ax_main.scatter(distance_mpc, conv_ratio, s=28, color="#1f77b4", alpha=0.75, edgecolor="0.2")
    ax_main.axhline(1.0, color="red", linestyle="--", linewidth=1.0, label="WIM anchor (ratio = 1)")
    if np.isfinite(slope):
        fit_x = np.linspace(distance_mpc.min(), distance_mpc.max(), 64)
        fit_log_y = slope * np.log10(fit_x) + _intercept
        ax_main.plot(
            fit_x,
            np.power(10.0, fit_log_y),
            color="#444444",
            linewidth=1.2,
            label=f"OLS log–log fit (slope {slope:+.2f})",
        )
    ax_main.set_xscale("log")
    ax_main.set_yscale("log")
    ax_main.set_ylabel(r"$\mathrm{convergence\_ratio\_to\_wim}$  (= $\rho_\mathrm{req} / 10^{-27}\,\mathrm{kg\,m^{-3}}$)")
    ax_main.set_title("Beam-smearing falsifiability test: ρ_req vs distance")
    ax_main.grid(True, which="both", alpha=0.25)
    _annotate(ax_main, slope, pearson, beam_arcmin)
    ax_main.legend(loc="lower right", fontsize=8)

    ax_beam.plot(distance_grid_mpc, beam_grid_kpc, color="#d62728", linewidth=1.5)
    ax_beam.set_xlabel("Distance (Mpc)")
    ax_beam.set_ylabel(r"$D_\mathrm{beam}$ (kpc)")
    ax_beam.set_title(f"Planck beam footprint on disk plane ({beam_arcmin:.1f}'  FWHM)")
    ax_beam.grid(True, which="both", alpha=0.25)
    ax_beam.set_yscale("log")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return output_path
