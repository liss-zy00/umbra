"""Population-level figures built from aggregated ``results/step_*`` CSV tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.reporting.plotting.companion import write_plot_companion_markdown


def save_discovery_rho_req_vs_wim_ratio(
    discovery_df: pd.DataFrame,
    output_path: Path,
    *,
    write_companion: bool = False,
) -> Path:
    """Scatter log10(required density) vs log10(WIM convergence ratio); colour by |galactic latitude|."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    need = {"required_density_for_match", "convergence_ratio_to_wim"}
    if not need.issubset(set(discovery_df.columns)):
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        ax.text(0.5, 0.5, "Missing columns for discovery scatter", ha="center", va="center", transform=ax.transAxes)
    else:
        rho = discovery_df["required_density_for_match"].to_numpy(dtype=float)
        ratio = discovery_df["convergence_ratio_to_wim"].to_numpy(dtype=float)
        ok = np.isfinite(rho) & np.isfinite(ratio) & (rho > 0.0) & (ratio > 0.0)
        if not np.any(ok):
            fig, ax = plt.subplots(figsize=(7.0, 5.0))
            ax.text(0.5, 0.5, "No finite positive rows", ha="center", va="center", transform=ax.transAxes)
        else:
            lx = np.log10(rho[ok])
            ly = np.log10(ratio[ok])
            fig, ax = plt.subplots(figsize=(7.0, 5.0))
            if "galactic_latitude_deg" in discovery_df.columns:
                cvals = np.abs(discovery_df["galactic_latitude_deg"].to_numpy(dtype=float)[ok])
                sc = ax.scatter(lx, ly, c=cvals, cmap="viridis", s=36, alpha=0.85, edgecolors="none")
                cbar = fig.colorbar(sc, ax=ax, shrink=0.82)
                cbar.set_label(r"$|b_{\mathrm{Gal}}|$ (deg)")
            else:
                ax.scatter(lx, ly, s=36, alpha=0.85, color="#1f77b4", edgecolors="none")
            ax.set_xlabel(r"$\log_{10}(\rho_{\mathrm{req}} / \mathrm{kg\,m^{-3}})$")
            ax.set_ylabel(r"$\log_{10}(\rho_{\mathrm{req}} / \rho_{\mathrm{WIM}})$")
            ax.set_title("Discovery census — required density vs WIM anchor ratio")
            ax.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    if write_companion:
        write_plot_companion_markdown(
            output_path,
            title="Discovery — rho_req vs WIM ratio",
            purpose="Population view of the Discovery-phase closure density and its ratio to the weakly interacting matter anchor.",
            data_shown="Rows from ``discovery_*.csv`` with finite positive ``required_density_for_match`` and ``convergence_ratio_to_wim``.",
            axes_and_units="Log10 SI density on x; log10 dimensionless WIM ratio on y; optional colour by Galactic latitude (deg).",
            how_to_read="Each point is one galaxy at the Discovery pass; extreme ratios flag objects whose inferred closure density departs strongly from the reference anchor.",
        )
    return output_path


def save_resolution_diagnostic_scatter(
    discovery_df: pd.DataFrame,
    output_path: Path,
    *,
    write_companion: bool = False,
) -> Path:
    """Beam FWHM scale (kpc) vs resolution correction factor (Discovery definition)."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    beam = discovery_df.get("beam_size_kpc")
    kfac = discovery_df.get("resolution_correction_factor")
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    if beam is None or kfac is None:
        ax.text(0.5, 0.5, "Missing beam_size_kpc or resolution_correction_factor", ha="center", va="center", transform=ax.transAxes)
    else:
        x = beam.to_numpy(dtype=float)
        y = kfac.to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if not np.any(ok):
            ax.text(0.5, 0.5, "No finite points", ha="center", va="center", transform=ax.transAxes)
        else:
            ax.scatter(x[ok], y[ok], s=38, alpha=0.8, color="#6d4c41", edgecolors="none")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("Planck beam FWHM scale  (kpc)")
            ax.set_ylabel("Resolution correction factor  k = D_beam / L_grad")
            ax.set_title("Discovery — beam scale vs resolution correction")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    if write_companion:
        write_plot_companion_markdown(
            output_path,
            title="Discovery — resolution diagnostic",
            purpose="Visualises how the beam smoothing scale compares to the adopted gradient scale proxy.",
            data_shown="``beam_size_kpc`` vs ``resolution_correction_factor`` from Discovery tables.",
            axes_and_units="Both axes logarithmic when data allow.",
            how_to_read="Large k means the beam is wide compared to the gradient length used in the proxy.",
        )
    return output_path


def save_evidence_coverage_vs_pearson(
    evidence_df: pd.DataFrame,
    output_path: Path,
    *,
    write_companion: bool = False,
) -> Path:
    """Headline coverage ratio vs Pearson r (validation scientific-evidence cohort)."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cov = evidence_df.get("theoretical_gap_coverage_ratio")
    pear = evidence_df.get("pearson_correlation")
    morph = evidence_df.get("morphology_class")
    status = evidence_df.get("validation_status")
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    if cov is None or pear is None:
        ax.text(0.5, 0.5, "Missing coverage or Pearson column", ha="center", va="center", transform=ax.transAxes)
    else:
        x = cov.to_numpy(dtype=float)
        y = pear.to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        if not np.any(ok):
            ax.text(0.5, 0.5, "No finite rows", ha="center", va="center", transform=ax.transAxes)
        else:
            sub = evidence_df.loc[ok]
            xs = sub["theoretical_gap_coverage_ratio"].to_numpy(dtype=float)
            ys = sub["pearson_correlation"].to_numpy(dtype=float)
            if status is not None and "validation_status" in sub.columns:
                for lab, grp in sub.groupby("validation_status"):
                    ax.scatter(
                        grp["theoretical_gap_coverage_ratio"].to_numpy(dtype=float),
                        grp["pearson_correlation"].to_numpy(dtype=float),
                        s=44,
                        alpha=0.88,
                        label=str(lab),
                    )
                ax.legend(title="validation_status", fontsize=8, title_fontsize=9)
            elif morph is not None and "morphology_class" in sub.columns:
                cmap = plt.get_cmap("tab20")
                for i, (lab, grp) in enumerate(sub.groupby("morphology_class")):
                    ax.scatter(
                        grp["theoretical_gap_coverage_ratio"].to_numpy(dtype=float),
                        grp["pearson_correlation"].to_numpy(dtype=float),
                        s=40,
                        alpha=0.85,
                        color=cmap(i % 20),
                        label=str(lab)[:32],
                    )
                ax.legend(title="morphology", fontsize=7, title_fontsize=8, loc="best")
            else:
                ax.scatter(xs, ys, s=44, alpha=0.88, color="#08519c", edgecolors="none")
            ax.axvline(1.0, color="#888888", linestyle="--", linewidth=1.0, alpha=0.7)
            ax.set_xlabel("Theoretical gap coverage ratio  (modulus vs baryonic deficit)")
            ax.set_ylabel("Pearson r  (residual vs magnetic drive, reliable shells)")
            ax.set_title("Validation supreme cohort — coverage vs correlation")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    if write_companion:
        write_plot_companion_markdown(
            output_path,
            title="Evidence — coverage vs Pearson",
            purpose="Synthesis scatter for galaxies in the scientific-evidence table.",
            data_shown="``theoretical_gap_coverage_ratio`` vs ``pearson_correlation``.",
            axes_and_units="Linear axes; vertical reference at unity coverage.",
            how_to_read="Points near unity on x with high positive r indicate tight amplitude-and-shape agreement between the residual and the magnetic proxy on retained shells.",
        )
    return output_path


def save_evidence_rar_vs_btfr_offsets(
    evidence_df: pd.DataFrame,
    output_path: Path,
    *,
    write_companion: bool = False,
) -> Path:
    """RAR median absolute residual (dex) vs BTFR residual (dex) when both present."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rar = evidence_df.get("rar_median_abs_residual_dex")
    btfr = evidence_df.get("btfr_residual_dex")
    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    if rar is None:
        ax.text(0.5, 0.5, "Missing rar_median_abs_residual_dex", ha="center", va="center", transform=ax.transAxes)
    else:
        xr = rar.to_numpy(dtype=float)
        if btfr is None:
            ax.text(0.5, 0.5, "Missing btfr_residual_dex", ha="center", va="center", transform=ax.transAxes)
        else:
            yb = btfr.to_numpy(dtype=float)
            ok = np.isfinite(xr) & np.isfinite(yb)
            if not np.any(ok):
                ax.text(0.5, 0.5, "No paired finite RAR/BTFR offsets", ha="center", va="center", transform=ax.transAxes)
            else:
                ax.scatter(xr[ok], yb[ok], s=44, alpha=0.88, color="#2ca02c", edgecolors="none")
                ax.set_xlabel("RAR median absolute residual  (dex)")
                ax.set_ylabel("BTFR residual  (dex)")
                ax.set_title("Scientific-evidence cohort — RAR vs BTFR residuals")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    if write_companion:
        write_plot_companion_markdown(
            output_path,
            title="RAR vs BTFR offsets",
            purpose="Cross-check kinematic scaling relations for the supreme merged cohort.",
            data_shown="``rar_median_abs_residual_dex`` vs ``btfr_residual_dex``.",
            axes_and_units="Both axes in dex (log10 space offsets).",
            how_to_read="Large |dex| values mark galaxies whose scaling-relation residuals dominate the interpretation.",
        )
    return output_path


def save_pearson_histogram_from_evidence(
    evidence_df: pd.DataFrame,
    output_path: Path,
    *,
    write_companion: bool = False,
) -> Path:
    pear = evidence_df.get("pearson_correlation")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    if pear is None:
        ax.text(0.5, 0.5, "Missing pearson_correlation", ha="center", va="center", transform=ax.transAxes)
    else:
        vals = pear.to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size > 0:
            ax.hist(vals, bins=min(20, max(5, vals.size)), alpha=0.78, color="#9467bd", edgecolor="white")
        ax.set_xlabel("Pearson r  (reliable shells)")
        ax.set_ylabel("Galaxy count")
        ax.set_title("Scientific-evidence cohort — Pearson r distribution")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    if write_companion:
        write_plot_companion_markdown(
            output_path,
            title="Pearson distribution",
            purpose="Histogram of best-grid Pearson correlations for the evidence table.",
            data_shown="``pearson_correlation`` column.",
            axes_and_units="Linear Pearson on x; counts on y.",
            how_to_read="Broad distribution indicates heterogeneous agreement between residual and magnetic drive shapes.",
        )
    return output_path


def save_validation_mean_gap_vs_hubble(
    validation_df: pd.DataFrame,
    output_path: Path,
    *,
    write_companion: bool = False,
) -> Path:
    """Mean outer-disk acceleration gap vs Hubble stage for the validation-reliable inventory."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    need = {"mean_acceleration_gap_edge", "hubble_type"}
    if not need.issubset(set(validation_df.columns)):
        ax.text(0.5, 0.5, "Missing mean_acceleration_gap_edge or hubble_type", ha="center", va="center", transform=ax.transAxes)
    else:
        x = validation_df["hubble_type"].to_numpy(dtype=float)
        y = validation_df["mean_acceleration_gap_edge"].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y) & (y > 0.0)
        if not np.any(ok):
            ax.text(0.5, 0.5, "No finite rows", ha="center", va="center", transform=ax.transAxes)
        else:
            ax.scatter(x[ok], y[ok], s=42, alpha=0.85, c="#17becf", edgecolors="none")
            ax.set_yscale("log")
            ax.set_xlabel("Hubble stage  T")
            ax.set_ylabel("Mean outer-disk acceleration gap  (m/s²)")
            ax.set_title("Validation reliable cohort — kinematic gap vs Hubble type")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    if write_companion:
        write_plot_companion_markdown(
            output_path,
            title="Validation gap vs Hubble T",
            purpose="Inventory-level view of the mean acceleration gap across reliable validation galaxies.",
            data_shown="``validation_reliable.csv`` (or equivalent schema).",
            axes_and_units="Hubble T linear; gap magnitude logarithmic (SI).",
            how_to_read="Each point is one galaxy; spread in y reflects heterogeneous outer-disk residuals after baryonic closure.",
        )
    return output_path
