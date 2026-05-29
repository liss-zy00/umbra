from __future__ import annotations

"""
Matplotlib-backed figures for per-galaxy diagnostics and cohort-level panels (Pearson histogram,
plasma β, morphology bars, etc.) used by tests and optional tooling. The main-run publication
PNG bundle is assembled in :mod:`src.reporting.publication_bundle` and :mod:`src.reporting.population_figures`.

Design conventions enforced here (data-visualization style guide):

* **Human nomenclature.** All axes / titles / legends spell out the physical quantity
  in English. Internal pipeline column names like ``acceleration_gap`` or
  ``magnetic_pressure_acceleration`` are *only* used in code; the figures expose
  the reader-facing names below (``LABEL_*``).
* **Reader-facing captions.** Cohort-level figures may carry a one-sentence footnote
  (``_apply_caption``). Per-galaxy publication PNGs omit the bottom figure caption to keep
  panels visually clean; optional Markdown companions document provenance when enabled.
* **Robust to sample size.** Marker alphas, log-binning and KDE / rug overlays are
  tuned so the same figure stays legible from 10 reliable galaxies up to several
  hundred shells.
* **Non-interactive backend ``Agg``** to avoid GUI dependencies in headless CI / Conda.
"""

from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402 — after backend selection
import numpy as np

from src.config.thresholds import ASTROPHYSICAL_DEFAULTS, RELIABILITY_CUTS
from src.models.galaxy import Galaxy
from .companion import write_plot_companion_markdown

# --- Reader-facing axis / legend labels (kept in one place so the language is consistent) ---
LABEL_GAP_ACCELERATION = "Gravitational Acceleration Deficiency  (m/s²)"
LABEL_MAGNETIC_DRIVE = "Magnetic Pressure Drive  (m/s²)"
LABEL_DM_EQUIVALENCE_PCT = "Dark Matter Equivalence  (%)"
LABEL_PLASMA_BETA = "Plasma Magnetization Index  (β)"
LABEL_MAGNETIC_FIELD_UG = "Inferred Magnetic Field Intensity  (µG)"
LABEL_RADIUS_KPC = "Galactocentric Radius  R  (kpc)"
LABEL_GALAXY_MAX_RADIUS_KPC = "Galaxy Outer Radius  R_max  (kpc)"
LABEL_PEARSON = "Pearson Correlation  r  (reliable shells)"
LABEL_BARYONIC_ACCEL = "Baryonic Centripetal Acceleration  g_bar  (m/s²)"
LABEL_OBSERVED_ACCEL = "Observed Centripetal Acceleration  g_obs  (m/s²)"

# Empirical Galactic disk B-field range used in nearby-galaxy literature (1–10 µG).
# Centralised in :mod:`src.config.thresholds`; aliased here for legibility / locality.
GALACTIC_REALISTIC_B_MIN_UG = RELIABILITY_CUTS.MAGNETIC_FIELD_REALISTIC_BAND_MIN_MICROGAUSS
GALACTIC_REALISTIC_B_MAX_UG = RELIABILITY_CUTS.MAGNETIC_FIELD_REALISTIC_BAND_MAX_MICROGAUSS

# Pearson interpretation bands (Cohen-style heuristic for noisy astronomical samples).
_PEARSON_STRONG_MIN = RELIABILITY_CUTS.PEARSON_INTERPRETATION_STRONG_MIN
_PEARSON_MODERATE_MIN = RELIABILITY_CUTS.PEARSON_INTERPRETATION_MODERATE_MIN

_FOOTNOTE_KW: dict[str, Any] = {
    "ha": "center",
    "va": "bottom",
    "fontsize": 8.5,
    "color": "#404040",
    "wrap": True,
}


def _annotation_color_for_colormap_cell(
    cmap: Any,
    norm: Any,
    value: float,
    *,
    luminance_threshold: float = 0.58,
) -> str:
    """
    Pick black or white annotation text so it stays readable on diverging / sequential fills.

    Uses relative luminance of the mapped RGBA face color (Rec. 709 weights).
    """
    from matplotlib import colors as mcolors

    rgba = mcolors.to_rgba(cmap(norm(float(value))))
    r, g, b = float(rgba[0]), float(rgba[1]), float(rgba[2])
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "#1a1a1a" if lum > luminance_threshold else "#ffffff"


def _apply_caption(fig: plt.Figure, caption: str) -> None:
    """Attach a single-line italic caption at the bottom of the figure."""
    fig.text(0.5, 0.005, caption, style="italic", **_FOOTNOTE_KW)
    fig.subplots_adjust(bottom=max(0.16, fig.subplotpars.bottom))


class GalaxyPlotter:
    """Figure writers only — never mutate ``Galaxy`` kinematic tables."""

    # ------------------------------------------------------------------ per-galaxy plots

    def save_force_plot(
        self,
        galaxy: Galaxy,
        galaxy_dir: Path,
        *,
        scenario_gradient_scale_kpc: float | None = None,
        scenario_ionization_fraction: float | None = None,
        scale_kpc: float | None = None,
        ion_frac: float | None = None,
        resolution_correction_factor_k: float | None = None,
        write_companion: bool = True,
        metadata_footer: str | None = None,
        y_accel_log_limits: tuple[float, float] | None = None,
    ) -> Path:
        """
        Per-galaxy radial comparison of the gravitational gap and the magnetic-pressure
        drive (log–log).

        When ``centripetal_acceleration_systematic`` exists in the dataframe, error bars
        propagate the inclination-systematic envelope onto |gap|.

        When ``resolution_correction_factor_k`` is finite and ``required_density_for_match`` is
        present, a third series overlays the WIM-anchored diagnostic: ``|a_B|`` rescaled so a shell
        whose closure density equals ``ρ_WIM · k`` would match the headline pressure drive — making
        the gap–versus–``ρ_req`` tension visible on the same radial axis as the gap curve.
        """
        sk = scenario_gradient_scale_kpc if scenario_gradient_scale_kpc is not None else scale_kpc
        ion = scenario_ionization_fraction if scenario_ionization_fraction is not None else ion_frac

        df = galaxy.rotation_dataframe
        fig, ax = plt.subplots(figsize=(8.5, 5.0))

        radius = df["radius_kpc"].to_numpy(dtype=float)
        gap_abs = df["acceleration_gap"].abs().to_numpy(dtype=float)
        mag_abs = df["magnetic_pressure_acceleration"].abs().to_numpy(dtype=float)

        gap_err = None
        if "centripetal_acceleration_systematic" in df.columns:
            sys_arr = np.abs(df["centripetal_acceleration_systematic"].to_numpy(dtype=float))
            sys_arr = np.where(np.isfinite(sys_arr), sys_arr, 0.0)
            if np.any(sys_arr > 0):
                lower = np.minimum(sys_arr, gap_abs * 0.99)
                gap_err = np.vstack([lower, sys_arr])

        gap_label = "Gravitational Acceleration Deficiency  ± inclination σ" if gap_err is not None else "Gravitational Acceleration Deficiency"
        if gap_err is not None:
            ax.errorbar(
                radius,
                gap_abs,
                yerr=gap_err,
                fmt="o-",
                markersize=3.8,
                linewidth=1.6,
                elinewidth=1.0,
                capsize=2.5,
                ecolor="#3182bd",
                color="#08519c",
                label=gap_label,
            )
        else:
            ax.plot(radius, gap_abs, "o-", color="#08519c", markersize=3.8, linewidth=1.6, label=gap_label)

        ax.plot(radius, mag_abs, "s-", color="#d95f0e", markersize=3.8, linewidth=1.6, label="Magnetic Pressure Drive")

        k_plot = resolution_correction_factor_k
        if (
            k_plot is not None
            and np.isfinite(float(k_plot))
            and float(k_plot) > 0.0
            and "required_density_for_match" in df.columns
        ):
            rho_anchor = float(ASTROPHYSICAL_DEFAULTS.WIM_REFERENCE_DENSITY_KG_M3)
            k_safe = max(float(k_plot), float(RELIABILITY_CUTS.RESOLUTION_CORRECTION_FACTOR_MIN_FOR_WIM_SCORE))
            rho_req = df["required_density_for_match"].to_numpy(dtype=float)
            floor = rho_anchor * float(RELIABILITY_CUTS.HONEST_ACCEL_RHO_REQ_RELATIVE_FLOOR)
            rho_eff = np.maximum(rho_req, floor)
            scale = (rho_anchor * k_safe) / rho_eff
            mag_honest = mag_abs * np.where(np.isfinite(scale) & np.isfinite(rho_req) & (rho_req > 0.0), scale, np.nan)
            ax.plot(
                radius,
                mag_honest,
                "--",
                color="#31a354",
                linewidth=1.6,
                markersize=0,
                label="Magnetic drive at WIM density × k (diagnostic)",
            )

        ax.set_xlabel(LABEL_RADIUS_KPC)
        ax.set_ylabel("Acceleration  (m/s²)")
        ax.set_yscale("log")
        if y_accel_log_limits is not None:
            ylo, yhi = (float(y_accel_log_limits[0]), float(y_accel_log_limits[1]))
            if np.isfinite(ylo) and np.isfinite(yhi) and ylo > 0.0 and yhi > ylo:
                ax.set_ylim(ylo, yhi)
        if sk is not None and ion is not None:
            scale_pc = int(round(float(sk) * 1000.0))
            ion_pct = float(ion) * 100.0
            title = (
                f"{galaxy.name} — Gravitational Deficit vs Magnetic Drive\n"
                f"(Δx scenario: {scale_pc} pc · ionized fraction: {ion_pct:.3g}%)"
            )
        else:
            title = f"{galaxy.name} — Gravitational Deficit vs Magnetic Drive"
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize=9)

        path = galaxy_dir / f"{galaxy.name}_acceleration_comparison.png"
        fig.tight_layout()
        fig.savefig(path, dpi=160)
        plt.close(fig)
        if write_companion:
            write_plot_companion_markdown(
                path,
                title=f"{galaxy.name} — Gravitational deficit vs magnetic drive",
                purpose=(
                    "Compares the magnitude of the gravitational acceleration deficit (kinematic excess over baryons) "
                    "to the magnetic pressure drive along the best-scenario rotation curve."
                ),
                data_shown=(
                    f"- Per-radius shells from ``Galaxy.rotation_dataframe`` for **{galaxy.name}**.\n"
                    "- **Gravitational Acceleration Deficiency:** ``|acceleration_gap|`` (SI m/s²).\n"
                    "- **Magnetic Pressure Drive:** ``|magnetic_pressure_acceleration|``.\n"
                    "- **WIM × k diagnostic (dashed):** ``|a_B| · (ρ_WIM · k / max(ρ_req, floor))`` when ``ρ_req`` and ``k`` are available.\n"
                    "- Optional error bars: ``centripetal_acceleration_systematic`` propagated from inclination uncertainty."
                ),
                axes_and_units=(
                    f"- X: {LABEL_RADIUS_KPC} (linear).\n"
                    "- Y: acceleration magnitude (m/s²), **logarithmic** scale."
                ),
                how_to_read=(
                    "Curves that track each other indicate the MHD pressure term can match the gap in amplitude across radius. "
                    "A magnetic curve far below the gap leaves an unexplained residual; far above suggests over-compensation or calibration issues. "
                    "The dashed green series shows the same magnetic drive rescaled to a WIM closure density multiplied by the beam dilution factor "
                    "so you can judge whether an observationally anchored density can follow the gap curve."
                ),
            )
        return path

    def save_rar_plot(
        self,
        galaxy: Galaxy,
        output_path: Path,
        *,
        a0_m_s2: float | None = None,
        write_companion: bool = True,
        x_log_limits: tuple[float, float] | None = None,
        y_log_limits: tuple[float, float] | None = None,
    ) -> Path:
        """
        Per-galaxy Radial Acceleration Relation (RAR) plot.

        X: ``g_bar`` (baryonic). Y: ``g_obs`` measured plus the McGaugh+2016 empirical fit
        and a 1:1 reference; magnetic-pressure |a_B| is overlaid as a second series.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        a0 = float(a0_m_s2) if a0_m_s2 is not None else float(ASTROPHYSICAL_DEFAULTS.RAR_A0_M_S2)
        df = galaxy.rotation_dataframe

        radius_m = df["radius_kpc"].to_numpy(dtype=float) * 3.0857e19  # kpc → m
        v_obs_ms = df["observed_velocity_kms"].to_numpy(dtype=float) * 1000.0
        v_bar_ms = df["baryonic_velocity_kms"].to_numpy(dtype=float) * 1000.0
        with np.errstate(divide="ignore", invalid="ignore"):
            g_obs = np.where(radius_m > 0, v_obs_ms**2 / radius_m, np.nan)
            g_bar = np.where(radius_m > 0, v_bar_ms**2 / radius_m, np.nan)

        a_b = (
            np.abs(df["magnetic_pressure_acceleration"].to_numpy(dtype=float))
            if "magnetic_pressure_acceleration" in df.columns
            else np.full_like(g_obs, np.nan)
        )

        ok = np.isfinite(g_obs) & np.isfinite(g_bar) & (g_obs > 0) & (g_bar > 0)

        fig, ax = plt.subplots(figsize=(7.2, 6.2))
        if not np.any(ok):
            ax.text(0.5, 0.5, "No reliable RAR points", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
        else:
            x_grid = np.logspace(np.log10(np.nanmin(g_bar[ok])) - 0.3, np.log10(np.nanmax(g_bar[ok])) + 0.3, 200)
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                y_fit = x_grid / (1.0 - np.exp(-np.sqrt(x_grid / a0)))
            ax.plot(x_grid, y_fit, color="#525252", linestyle="--", linewidth=1.6, label="RAR (McGaugh+2016)")
            ax.plot(x_grid, x_grid, color="#bdbdbd", linestyle=":", linewidth=1.0, label="1:1 (Newtonian)")
            ax.scatter(g_bar[ok], g_obs[ok], color="#08519c", s=26, alpha=0.85, label="Measured  g_obs")
            mask_ab = ok & np.isfinite(a_b) & (a_b > 0)
            if np.any(mask_ab):
                ax.scatter(g_bar[mask_ab], a_b[mask_ab], color="#d95f0e", marker="^", s=26, alpha=0.85, label="Magnetic Pressure Drive")

            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.axvline(a0, color="#999999", linewidth=0.8, alpha=0.45)
            ax.set_xlabel(LABEL_BARYONIC_ACCEL)
            ax.set_ylabel(LABEL_OBSERVED_ACCEL + "   /   Magnetic Drive")
            ax.set_title(f"{galaxy.name} — Radial Acceleration Relation")
            ax.legend(loc="lower right", fontsize=9)
            ax.grid(True, which="both", alpha=0.25)
            if x_log_limits is not None:
                xlo, xhi = (float(x_log_limits[0]), float(x_log_limits[1]))
                if np.isfinite(xlo) and np.isfinite(xhi) and xlo > 0.0 and xhi > xlo:
                    ax.set_xlim(xlo, xhi)
            if y_log_limits is not None:
                ylo, yhi = (float(y_log_limits[0]), float(y_log_limits[1]))
                if np.isfinite(ylo) and np.isfinite(yhi) and ylo > 0.0 and yhi > ylo:
                    ax.set_ylim(ylo, yhi)

        fig.tight_layout()
        fig.savefig(output_path, dpi=160)
        plt.close(fig)
        if write_companion:
            write_plot_companion_markdown(
                output_path,
                title=f"{galaxy.name} — Radial Acceleration Relation (RAR)",
                purpose=(
                    "Places each galaxy on the empirical Radial Acceleration Relation (McGaugh+2016 reference) "
                    "and overlays the inferred magnetic pressure drive so the reader can judge whether MHD terms "
                    "could shift baryonic points toward the observed acceleration."
                ),
                data_shown=(
                    f"- **{galaxy.name}** rotation-curve shells: ``g_obs = V_obs²/r`` and ``g_bar = V_bar²/r`` (SI m/s²).\n"
                    "- Reference curve: McGaugh+2016 RAR with ``a0`` from pipeline astrophysical defaults.\n"
                    "- Orange triangles: ``|magnetic_pressure_acceleration|`` at each radius (same x as ``g_bar``)."
                ),
                axes_and_units=(
                    f"- X: {LABEL_BARYONIC_ACCEL} (log scale).\n"
                    f"- Y: {LABEL_OBSERVED_ACCEL} — both ``g_obs`` and magnetic drive share the vertical log axis in the legend."
                ),
                how_to_read=(
                    "Points near the gray dashed RAR curve behave like typical disk galaxies in the SPARC literature. "
                    "If orange triangles fall near that curve while blue points sit above the baryonic 1:1 line, "
                    "the magnetic term is a candidate partial substitute for a dark-matter halo in that object."
                ),
            )
        return output_path

    def save_acceleration_vs_radius_sample(
        self, galaxy: Galaxy, output_path: Path, *, write_companion: bool = True
    ) -> Path:
        """Signed acceleration components vs radius (diagnostic; preserves engine sign)."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df = galaxy.rotation_dataframe
        fig, ax = plt.subplots(figsize=(8.5, 5.0))
        ax.plot(df["radius_kpc"], df["acceleration_gap"], label="Gravitational Acceleration Deficiency (signed)", linewidth=1.7)
        ax.plot(df["radius_kpc"], df["magnetic_pressure_acceleration"], label="Magnetic Pressure Drive (signed)", linewidth=1.7)
        ax.set_xlabel(LABEL_RADIUS_KPC)
        ax.set_ylabel("Acceleration  (m/s²)")
        ax.set_title(f"{galaxy.name} — Signed Accelerations vs Radius")
        ax.axhline(0.0, color="k", linewidth=0.6, alpha=0.4)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(output_path, dpi=160)
        plt.close(fig)
        if write_companion:
            write_plot_companion_markdown(
                output_path,
                title=f"{galaxy.name} — Signed accelerations vs radius",
                purpose="Diagnostic plot preserving the engine sign convention for gap and magnetic drive along radius.",
                data_shown=(
                    f"- **{galaxy.name}**: signed ``acceleration_gap`` and signed ``magnetic_pressure_acceleration`` vs ``radius_kpc``."
                ),
                axes_and_units=(
                    f"- X: {LABEL_RADIUS_KPC} (linear).\n"
                    "- Y: acceleration (m/s²), **linear** scale; zero line marks sign flip."
                ),
                how_to_read=(
                    "Same-sign overlap at a radius means the magnetic pressure gradient pushes in the same direction as the "
                    "kinematic deficit at that shell; opposite signs flag a physical inconsistency or a bad scenario choice."
                ),
            )
        return output_path

    # ------------------------------------------------------------------ population-level plots

    def save_plasma_beta_vs_coverage(
        self,
        supreme_rows: list[dict[str, Any]],
        output_path: Path,
    ) -> Path:
        """
        Median plasma β vs Dark Matter Equivalence (%) for the reliable synthesis rows.
        Vertical line at β = 1 separates magnetically-dominated (β < 1) from thermally-dominated (β > 1).
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        betas: list[float] = []
        equiv_pct: list[float] = []
        for row in supreme_rows:
            b = row.get("plasma_beta_thermal_median")
            r = row.get("theoretical_gap_coverage_ratio")
            if (
                isinstance(b, (int, float))
                and np.isfinite(float(b))
                and float(b) > 0
                and isinstance(r, (int, float))
                and np.isfinite(float(r))
                and float(r) > 0
            ):
                betas.append(float(b))
                equiv_pct.append(float(r) * 100.0)

        fig, ax = plt.subplots(figsize=(7.5, 5.6))
        arr_b = np.asarray(betas, dtype=float)
        arr_y = np.asarray(equiv_pct, dtype=float)
        if arr_b.size == 0:
            ax.text(0.5, 0.5, "No synthesis rows with β and Dark Matter Equivalence", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
        else:
            ax.scatter(arr_b, arr_y, c="#31a354", s=44, alpha=0.85, edgecolors="white", linewidths=0.6)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.axvline(1.0, color="#737373", linestyle="--", linewidth=1.2, label="β = 1  (magnetic ↔ thermal balance)")
            ax.axhline(100.0, color="#525252", linestyle=":", linewidth=1.0, alpha=0.6, label="100 % equivalence")
            ax.set_xlabel(LABEL_PLASMA_BETA)
            ax.set_ylabel(LABEL_DM_EQUIVALENCE_PCT)
            ax.set_title("Plasma Magnetization vs Dark Matter Equivalence  (reliable sample)")
            ax.legend(loc="best", fontsize=9, framealpha=0.85)
            ax.grid(True, which="both", alpha=0.25)
        fig.tight_layout()
        _apply_caption(
            fig,
            "Galaxies on the left of β = 1 are in the regime where the MHD-pressure assumption holds; "
            "a clustering at high y on the left would be the strongest empirical support for the magnetic-explanation hypothesis.",
        )
        fig.savefig(output_path, dpi=160)
        plt.close(fig)
        n_gal = int(arr_b.size)
        write_plot_companion_markdown(
            output_path,
            title="Plasma magnetization vs Dark Matter Equivalence",
            purpose="Tests whether magnetically dominated disks (low plasma β) systematically achieve better gap coverage.",
            data_shown=(
                f"- **{n_gal}** reliable synthesis galaxies with finite ``plasma_beta_thermal_median`` and ``theoretical_gap_coverage_ratio``.\n"
                "- Y-axis values are ``theoretical_gap_coverage_ratio × 100`` (Dark Matter Equivalence %)."
            ),
            axes_and_units=(
                f"- X: {LABEL_PLASMA_BETA} (log). Reference line at β = 1.\n"
                f"- Y: {LABEL_DM_EQUIVALENCE_PCT} (log). Dotted line at 100 %."
            ),
            how_to_read=(
                "If magnetic dominance were the sole driver of success, galaxies would cluster at low β and high equivalence. "
                "A diffuse cloud falsifies that simple story and motivates mixed (β, calibration, morphology) explanations."
            ),
        )
        return output_path

    def save_magnetic_field_vs_galaxy_scale(
        self,
        supreme_rows: list[dict[str, Any]],
        output_path: Path,
    ) -> Path:
        """
        Inferred Magnetic Field Intensity (µG) versus the galaxy outer radius (kpc).

        Shaded band marks the 1–10 µG ``Realistic Galactic Range`` reported in nearby-galaxy
        polarization surveys. Galaxies sitting outside this band over a wide range of R_max
        are an empirical red flag for the calibration anchor / Δx choice.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        radii: list[float] = []
        b_field: list[float] = []
        names: list[str] = []
        for row in supreme_rows:
            r_max = row.get("r_max_kpc")
            b = row.get("magnetic_field_microgauss")
            if (
                isinstance(r_max, (int, float))
                and np.isfinite(float(r_max))
                and float(r_max) > 0
                and isinstance(b, (int, float))
                and np.isfinite(float(b))
                and float(b) > 0
            ):
                radii.append(float(r_max))
                b_field.append(float(b))
                names.append(str(row.get("galaxy_name", "")))

        fig, ax = plt.subplots(figsize=(8.0, 5.6))
        rx = np.asarray(radii, dtype=float)
        by = np.asarray(b_field, dtype=float)
        if rx.size == 0:
            ax.text(0.5, 0.5, "No reliable galaxies with R_max and B field", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
        else:
            ax.axhspan(GALACTIC_REALISTIC_B_MIN_UG, GALACTIC_REALISTIC_B_MAX_UG, color="#a1d99b", alpha=0.30,
                       label="Realistic Galactic Range  1–10 µG")
            ax.scatter(rx, by, c="#08519c", s=42, alpha=0.80, edgecolors="white", linewidths=0.6,
                       label="Reliable galaxies")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(LABEL_GALAXY_MAX_RADIUS_KPC)
            ax.set_ylabel(LABEL_MAGNETIC_FIELD_UG)
            ax.set_title("Inferred Magnetic Field vs Galaxy Outer Radius")
            ax.legend(loc="best", fontsize=9, framealpha=0.85)
            ax.grid(True, which="both", alpha=0.25)
        fig.tight_layout()
        _apply_caption(
            fig,
            "Points inside the green band reproduce the µG-scale fields measured in nearby spirals; "
            "galaxies straying systematically outside reveal calibration / Δx artefacts rather than a real population trend.",
        )
        fig.savefig(output_path, dpi=160)
        plt.close(fig)
        n_pts = int(rx.size)
        write_plot_companion_markdown(
            output_path,
            title="Inferred magnetic field vs galaxy outer radius",
            purpose="Sanity check that inferred B fields stay in the µG range expected for disk galaxies as a function of spatial scale.",
            data_shown=(
                f"- **{n_pts}** reliable galaxies with ``r_max_kpc`` from the rotation curve and ``magnetic_field_microgauss`` from the best scenario.\n"
                "- Shaded band: 1–10 µG (literature reference range for nearby disks; illustrative, not a hard cut)."
            ),
            axes_and_units=(
                f"- X: {LABEL_GALAXY_MAX_RADIUS_KPC} (log).\n"
                f"- Y: {LABEL_MAGNETIC_FIELD_UG} (log)."
            ),
            how_to_read=(
                "Points inside the band are broadly consistent with observed disk fields; a systematic locus far above or below "
                "the band across many objects flags calibration or Δx issues rather than random scatter."
            ),
        )
        return output_path

    def save_morphology_comparison(
        self,
        stats_by_type: dict[str, dict[str, float]],
        output_path: Path,
    ) -> Path:
        """
        Horizontal bar chart of mean Dark Matter Equivalence (%) for every morphology bucket
        present in the run, plus the photometric-LSB pseudo-bucket.

        Parameters
        ----------
        stats_by_type:
            Mapping ``{label: {"mean_pct": float, "count": int}}``. Empty buckets (count=0)
            are skipped so the figure does not advertise classes the run never produced.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Drop empty buckets — keeps the bar chart compact and avoids implying we have data
        # for morphology classes that are absent from the current run.
        populated = {
            str(label): {
                "mean_pct": float(payload.get("mean_pct", 0.0) or 0.0),
                "count": int(payload.get("count", 0) or 0),
            }
            for label, payload in stats_by_type.items()
            if int(payload.get("count", 0) or 0) > 0
        }
        labels = sorted(
            populated.keys(),
            key=lambda lab: (-populated[lab]["count"], lab),
        )

        fig, ax = plt.subplots(figsize=(8.6, max(4.6, 0.55 * max(len(labels), 1) + 2.2)))
        if not labels:
            ax.text(
                0.5,
                0.5,
                "No reliable galaxies — morphology comparison skipped.",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=11,
            )
            ax.set_axis_off()
            fig.tight_layout()
            _apply_caption(fig, "Awaiting reliable data from at least one morphology bucket.")
            fig.savefig(output_path, dpi=160)
            plt.close(fig)
            write_plot_companion_markdown(
                output_path,
                title="Morphology comparison — mean Dark Matter Equivalence",
                purpose=(
                    "Population-level breakdown of mean Dark Matter Equivalence (%) by SPARC "
                    "morphology class plus a photometric-LSB pseudo-class."
                ),
                data_shown=(
                    "- Empty state: no reliable galaxy contributed a finite "
                    "``theoretical_gap_coverage_ratio`` for any morphology class in this run."
                ),
                axes_and_units=(
                    f"- Y: morphology class labels.\n"
                    f"- X: {LABEL_DM_EQUIVALENCE_PCT} (linear)."
                ),
                how_to_read=(
                    "Empty figure: the audit gates rejected every galaxy — see the "
                    "``results/step_validation/validation_consistency_audit.md`` failure-mode breakdown."
                ),
            )
            return output_path

        means = np.asarray([populated[lab]["mean_pct"] for lab in labels], dtype=float)
        counts = np.asarray([populated[lab]["count"] for lab in labels], dtype=int)
        y_pos = np.arange(len(labels), dtype=float)

        # Tab10 colour cycle — distinct hue per class, deterministic regardless of label order.
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i % 10) for i in range(len(labels))]
        bars = ax.barh(y_pos, means, color=colors, edgecolor="white", linewidth=0.6)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([f"{lab}  (N={counts[i]})" for i, lab in enumerate(labels)], fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel(LABEL_DM_EQUIVALENCE_PCT)
        ax.set_title("Morphology Comparison — mean Dark Matter Equivalence per class")
        ax.grid(True, axis="x", alpha=0.25)

        # Reference line at 100 % so readers can spot classes that are over- or under-compensating.
        ax.axvline(100.0, color="#252525", linestyle="--", linewidth=1.0, alpha=0.7, label="Perfect Equivalence (100 %)")

        # Annotate each bar with the mean value (placed inside if the bar is long enough,
        # otherwise just outside, so short bars stay legible).
        if means.size:
            x_max = float(np.nanmax(means)) if np.any(np.isfinite(means)) else 0.0
            for rect, value in zip(bars, means, strict=True):
                w = float(rect.get_width())
                cx = rect.get_y() + rect.get_height() / 2.0
                inside = w > max(0.18 * x_max, 5.0)
                ax.text(
                    w + (-0.01 * x_max if inside else 0.01 * x_max),
                    cx,
                    f"{value:.1f} %",
                    ha="right" if inside else "left",
                    va="center",
                    fontsize=8.5,
                    color="white" if inside else "#1a1a1a",
                )
        ax.legend(loc="lower right", fontsize=8, framealpha=0.85)

        fig.tight_layout()
        # Reserve extra room at the bottom so the wrapped caption never collides with the x-axis label.
        fig.subplots_adjust(bottom=max(0.28, fig.subplotpars.bottom))
        _apply_caption(
            fig,
            "Each bar is the mean Dark Matter Equivalence (%) for reliable galaxies in that morphology class. "
            "The dashed line marks 100 % (perfect coverage); the photometric-LSB row is a pseudo-class.",
        )
        fig.savefig(output_path, dpi=160)
        plt.close(fig)

        lbl_txt = "; ".join(
            f"{lab}: {populated[lab]['mean_pct']:.2f} % (N={populated[lab]['count']})" for lab in labels
        )
        write_plot_companion_markdown(
            output_path,
            title="Morphology comparison — mean Dark Matter Equivalence per class",
            purpose=(
                "Population-level breakdown of mean Dark Matter Equivalence (%) by SPARC morphology class, "
                "plus a photometric-LSB pseudo-class. Earlier versions of this figure compared only LSB and "
                "early-type spirals; the current figure preserves that historical contrast (LSB galaxies "
                "demand the largest dark-matter halos in the standard picture, so they are the strongest "
                "single test for the magnetic hypothesis) while showing every class actually present in the run."
            ),
            data_shown=(
                f"- One bar per non-empty morphology bucket (count shown next to the label).\n"
                f"- Source: ``theoretical_gap_coverage_ratio × 100`` averaged across reliable galaxies "
                f"that share the bucket label.\n"
                f"- Sample stats: {lbl_txt}."
            ),
            axes_and_units=(
                f"- Y: morphology class labels (sorted by sample size).\n"
                f"- X: {LABEL_DM_EQUIVALENCE_PCT} (linear).\n"
                "- Dashed vertical: 100 % perfect equivalence reference."
            ),
            how_to_read=(
                "Bars near 100 % indicate classes where the magnetic-pressure drive supplies, on average, "
                "the missing dark-matter equivalent. Bars far above 100 % suggest over-compensation "
                "(possible foreground contamination); bars near 0 % suggest the model under-explains that "
                "morphology bucket. A flat profile across all classes is consistent with a morphology-blind "
                "explanation; strong asymmetries flag the populations where the hypothesis works (or fails)."
            ),
        )
        return output_path

    def save_pearson_correlation_distribution(self, pearson_values: list[float], output_path: Path) -> Path:
        """
        Histogram of Pearson correlation between |gap| and |magnetic drive| across reliable
        galaxies, with shaded interpretation bands (Strong, Moderate, Noise).
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        arr = np.asarray(pearson_values, dtype=float)
        arr = arr[np.isfinite(arr)]

        fig, ax = plt.subplots(figsize=(8.0, 4.6))

        # Plot interpretation bands first so bars sit on top.
        ax.axvspan(_PEARSON_STRONG_MIN, 1.0, color="#a1d99b", alpha=0.40, label="Strong  (|r| ≥ 0.7)")
        ax.axvspan(_PEARSON_MODERATE_MIN, _PEARSON_STRONG_MIN, color="#fee391", alpha=0.50, label="Moderate  (0.4 ≤ |r| < 0.7)")
        ax.axvspan(-_PEARSON_MODERATE_MIN, _PEARSON_MODERATE_MIN, color="#d9d9d9", alpha=0.55, label="Noise  (|r| < 0.4)")
        ax.axvspan(-_PEARSON_STRONG_MIN, -_PEARSON_MODERATE_MIN, color="#fee391", alpha=0.40)
        ax.axvspan(-1.0, -_PEARSON_STRONG_MIN, color="#a1d99b", alpha=0.30)

        if arr.size:
            ax.hist(arr, bins=min(30, max(5, arr.size)), color="#54278f", edgecolor="white", alpha=0.95, zorder=3)
            ax.axvline(float(np.median(arr)), color="#252525", linestyle="--", linewidth=1.2, label=f"median r = {float(np.median(arr)):.2f}")

        ax.set_xlim(-1.05, 1.05)
        ax.set_xlabel(LABEL_PEARSON)
        ax.set_ylabel("Number of reliable galaxies")
        ax.set_title("Population Distribution — Pearson Correlation between Magnetic Drive and Gap")
        ax.legend(loc="upper left", fontsize=8, framealpha=0.85)
        ax.grid(True, axis="y", alpha=0.20)
        fig.tight_layout()
        _apply_caption(
            fig,
            "Higher r means the magnetic pressure drive tracks the gravitational deficit shape across radii, "
            "not just its amplitude — the fundamental falsifiability test of the MHD-replaces-DM hypothesis.",
        )
        fig.savefig(output_path, dpi=160)
        plt.close(fig)
        write_plot_companion_markdown(
            output_path,
            title="Pearson correlation — gap vs magnetic drive",
            purpose="Distribution of how tightly |gap| tracks |magnetic drive| in radius across each reliable galaxy.",
            data_shown=(
                f"- **{int(arr.size)}** Pearson coefficients (one per reliable galaxy best scenario).\n"
                "- Correlation computed on shells that pass ``reliable_measurement`` in the pipeline."
            ),
            axes_and_units=(
                f"- X: {LABEL_PEARSON} (linear, −1 to +1).\n"
                "- Y: histogram count.\n"
                "- Background bands: Strong (|r|≥0.7), Moderate (0.4–0.7), Noise (<0.4)."
            ),
            how_to_read=(
                "Positive r toward +1 means the magnetic drive rises and falls with the gap in the same sense across radius — "
                "the shape test of the hypothesis. Values near zero mean no linear relationship in log or linear space as encoded by the metric builder."
            ),
        )
        return output_path

    # ------------------------------------------------------------------ magnetic field distribution (standalone)

    def save_magnetic_field_histogram(
        self,
        microgauss_values: list[float],
        output_path: Path,
        *,
        plausibility_band: tuple[float, float] | None = None,
    ) -> Path:
        """
        Log-binned histogram of the per-galaxy median **Inferred Magnetic Field Intensity (µG)**
        for all reliable galaxies. Optional shaded band marks the pipeline plausibility range
        (often aligned with 1–10 µG). Rug marks at the bottom show individual galaxies.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        arr = np.asarray(microgauss_values, dtype=float)
        arr = arr[np.isfinite(arr) & (arr > 0)]

        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        if arr.size == 0:
            ax.text(0.5, 0.5, "No reliable galaxies with B field", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
        else:
            log_arr = np.log10(arr)
            n_bins = max(5, min(25, int(np.sqrt(arr.size) * 2)))
            counts, edges, _patches = ax.hist(log_arr, bins=n_bins, color="#31a354", edgecolor="white", alpha=0.85)
            ymax = float(np.max(counts)) if counts.size else 1.0
            if plausibility_band is not None:
                b_min, b_max = float(plausibility_band[0]), float(plausibility_band[1])
                if b_min > 0 and b_max > b_min:
                    ax.axvspan(np.log10(b_min), np.log10(b_max), color="#74c476", alpha=0.18, label="Realistic Galactic Range  1–10 µG")
            median_val = float(np.median(arr))
            ax.axvline(np.log10(median_val), color="#252525", linestyle="--", linewidth=1.2, label=f"median = {median_val:.2f} µG")
            # Rug plot: one tick per galaxy at the histogram baseline (does not occlude other figures).
            for lv in log_arr:
                ax.plot([lv, lv], [0.0, 0.02 * ymax], color="#1b5e20", linewidth=0.9, clip_on=True)
            ax.set_xlabel("log₁₀(Inferred Magnetic Field Intensity  /  µG)")
            ax.set_ylabel("Number of reliable galaxies")
            ax.set_title("Inferred Magnetic Field Intensity — Reliable Sample Distribution")
            ax.legend(loc="upper right", fontsize=9)
            ax.grid(True, which="major", alpha=0.25)
        fig.tight_layout()
        _apply_caption(
            fig,
            "The distribution summarizes the median µG field each reliable galaxy needs for its best scenario; "
            "a median outside the shaded band flags systematic calibration or Δx stress.",
        )
        fig.savefig(output_path, dpi=160)
        plt.close(fig)
        n_gal = int(arr.size)
        band_note = ""
        if plausibility_band is not None:
            band_note = f"\n- Shaded band: [{plausibility_band[0]}, {plausibility_band[1]}] µG (pipeline plausibility / reference)."
        write_plot_companion_markdown(
            output_path,
            title="Inferred magnetic field — reliable sample distribution",
            purpose="Shows the spread of median inferred B (µG) across galaxies that passed reliability and inclination gates.",
            data_shown=(
                f"- **{n_gal}** per-galaxy median values (µG), one count per reliable synthesis row.{band_note}\n"
                "- X-axis: log₁₀(B / µG); green rug: individual galaxies."
            ),
            axes_and_units=(
                f"- X: log₁₀({LABEL_MAGNETIC_FIELD_UG}).\n"
                "- Y: galaxy count per bin."
            ),
            how_to_read=(
                "A narrow peak inside the realistic band means inferred fields are astrophysically modest; "
                "a long tail to high µG suggests many objects require extreme fields under the current MHD parametrization."
            ),
        )
        return output_path

