"""
Off-screen PyVista renders for pipeline analysis tables.

Single-panel toy-disk composition from tabular profiles, not an N-body or self-consistent
Poisson+MHD run. Optional layers (gap glyphs, WIM annuli, streamlines) default to **auto**
from CSV columns and morphology unless overridden explicitly (tests).
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Any

import numpy as np

from src.config.thresholds import VIZ_DISPLAY
from src.physics.diagnostics.tension import supports_logarithmic_spiral
from src.visualization.batch import write_viz_manifest
from src.visualization.io import load_analysis_dataframe, load_optional_summary
from src.visualization.scene import (
    glyphs_from_layer,
    illustrative_b_streamlines,
    layer_from_radial_profile,
    reference_disk,
    wim_density_annulus_meshes,
)

logger = logging.getLogger(__name__)

_CAPTION_FOOTER = (
    "Graphic composition from pipeline CSV profiles — not N-body / not self-consistent "
    "Poisson+MHD. Pedagogical toy disk geometry."
)

_CAPTION_STREAMLINES = (
    "Gold streamlines: illustrative in-plane azimuthal field scaled by "
    "magnetic_field_microgauss(r) with optional PA offset — not a full 3D Planck B map."
)


def _scenario_caption(summary: dict[str, Any] | None) -> str:
    if not summary:
        return "Scenario: (summary JSON not found — see analysis CSV columns)."
    scen = summary.get("scenario") or summary.get("best_scenario") or {}
    if isinstance(scen, dict):
        parts = []
        for key in ("scenario_gradient_scale_kpc", "scenario_ionization_fraction"):
            if scen.get(key) is not None:
                parts.append(f"{key}={scen.get(key)}")
        if parts:
            return "Scenario knobs: " + "; ".join(parts)
    return "Scenario: (see per-galaxy summary JSON)."


def _sanitize_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "galaxy"


def _column_si(df: Any, name: str) -> np.ndarray:
    if name not in df.columns:
        return np.zeros(len(df), dtype=float)
    return np.nan_to_num(df[name].to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)


def _position_angle_deg(summary: dict[str, Any] | None) -> float | None:
    if not summary:
        return None
    kin = summary.get("kinematic_metadata") or {}
    v = kin.get("position_angle_deg")
    if v is None:
        return None
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def _morphology_class_from_summary(summary: dict[str, Any] | None) -> str | None:
    if not summary:
        return None
    ident = summary.get("galaxy_identity") or {}
    v = ident.get("morphology_class")
    if v is None:
        return None
    s = str(v).strip()
    return s or None


def resolve_auto_viz_layers(
    df: Any,
    summary: dict[str, Any] | None,
    *,
    morphology_class_override: str | None = None,
    gap_glyphs: bool | None = None,
    wim_rings: bool | None = None,
    streamlines: bool | None = None,
) -> tuple[bool, bool, bool]:
    """
    When a flag is None, derive it from CSV columns and spiral-like morphology policy.

    ``morphology_class_override`` is used for streamlines when summary JSON lacks identity
    (e.g. supreme cohort row merged from validation tables).
    """
    has_gap_col = "acceleration_gap" in df.columns
    rho_ion = _column_si(df, "effective_ionized_mass_density")
    wim_data_ok = (
        "effective_ionized_mass_density" in df.columns
        and bool(np.any(np.isfinite(rho_ion)))
        and bool(np.any(rho_ion > 0.0))
    )
    morph = morphology_class_override if morphology_class_override else _morphology_class_from_summary(summary)

    use_gap = has_gap_col if gap_glyphs is None else bool(gap_glyphs)
    use_wim = wim_data_ok if wim_rings is None else bool(wim_rings)
    use_sl = supports_logarithmic_spiral(morph) if streamlines is None else bool(streamlines)
    return use_gap, use_wim, use_sl


def _caption_extras(
    *,
    has_gap_column: bool,
    gap_glyphs: bool,
    wim_rings: bool,
    streamlines: bool,
) -> str:
    parts: list[str] = []
    if has_gap_column:
        parts.append("acceleration_gap column present in CSV (kinematic deficit proxy).")
        if gap_glyphs:
            parts.append("Green arrows: acceleration_gap (toy radial display scale).")
    if wim_rings:
        parts.append("Annuli: log10 median effective_ionized_mass_density per radial bin (pseudo-3D).")
    if streamlines:
        parts.append(_CAPTION_STREAMLINES)
    return " ".join(parts) if parts else ""


def _add_wim_annuli(
    pl: Any,
    r_kpc: np.ndarray,
    rho: np.ndarray,
    *,
    n_bins: int,
    z_plot: float,
) -> None:
    for mesh in wim_density_annulus_meshes(r_kpc, rho, n_bins=n_bins, z_plot=z_plot):
        pl.add_mesh(
            mesh,
            scalars="log10_rho_ion_kg_m3",
            cmap="plasma",
            opacity=0.34,
            show_scalar_bar=False,
        )


def _add_streamlines(
    pl: Any,
    r_max: float,
    r_kpc: np.ndarray,
    b_microgauss: np.ndarray,
    pa_deg: float | None,
    *,
    grid_n: int,
    z_plot: float,
) -> None:
    sl = illustrative_b_streamlines(
        r_max,
        r_kpc,
        b_microgauss,
        pa_deg,
        grid_n=grid_n,
        z_plot=z_plot,
    )
    if sl is not None and sl.n_points > 0:
        try:
            tubed = sl.tube(radius=max(r_max * 0.012, 1.0e-6))
            pl.add_mesh(tubed, smooth_shading=True, color="gold")
        except (AttributeError, ValueError):
            pl.add_mesh(sl, line_width=2.0, color="gold")


def render_galaxy_png(
    project_root: Path,
    galaxy_name: str,
    output_path: Path,
    *,
    window_size: tuple[int, int] = (1280, 720),
    gap_glyphs: bool | None = None,
    wim_rings: bool | None = None,
    streamlines: bool | None = None,
    morphology_class_override: str | None = None,
    wim_bins: int | None = None,
    streamline_grid: int | None = None,
    export_html: bool = True,
) -> None:
    import pyvista as pv

    df = load_analysis_dataframe(project_root, galaxy_name)
    summary = load_optional_summary(project_root, galaxy_name)

    r_kpc = df["radius_kpc"].to_numpy(dtype=float)
    g_bar = df["baryonic_acceleration"].to_numpy(dtype=float)
    g_obs = df["centripetal_acceleration"].to_numpy(dtype=float)
    has_gap = "acceleration_gap" in df.columns
    g_gap = _column_si(df, "acceleration_gap") if has_gap else np.zeros_like(g_bar)
    rho_ion = _column_si(df, "effective_ionized_mass_density")
    b_mu = _column_si(df, "magnetic_field_microgauss")
    pa_deg = _position_angle_deg(summary)

    gap_glyphs, wim_rings, streamlines = resolve_auto_viz_layers(
        df,
        summary,
        morphology_class_override=morphology_class_override,
        gap_glyphs=gap_glyphs,
        wim_rings=wim_rings,
        streamlines=streamlines,
    )

    r_max = float(np.nanmax(r_kpc)) if r_kpc.size else 1.0
    n_wim = int(wim_bins) if wim_bins is not None else int(VIZ_DISPLAY.WIM_RING_BIN_COUNT)
    grid_n = int(streamline_grid) if streamline_grid is not None else int(VIZ_DISPLAY.STREAMLINE_GRID_SIZE)

    extras = _caption_extras(
        has_gap_column=has_gap,
        gap_glyphs=gap_glyphs,
        wim_rings=wim_rings,
        streamlines=streamlines,
    )
    extra_block = ("\n" + extras) if extras else ""

    pl = pv.Plotter(off_screen=True, window_size=window_size)
    if wim_rings:
        _add_wim_annuli(pl, r_kpc, rho_ion, n_bins=n_wim, z_plot=0.02)
    disk = reference_disk(r_max)
    pl.add_mesh(disk, color="#c8c8c8", opacity=0.35, smooth_shading=True)
    layer_bar = layer_from_radial_profile(r_kpc, g_bar, r_max_kpc=r_max)
    layer_obs = layer_from_radial_profile(r_kpc, g_obs, r_max_kpc=r_max, z_offset=0.06)
    pl.add_mesh(glyphs_from_layer(layer_bar), color="#1f77b4")
    pl.add_mesh(glyphs_from_layer(layer_obs), color="#d62728")
    if gap_glyphs and has_gap:
        layer_gap = layer_from_radial_profile(r_kpc, g_gap, r_max_kpc=r_max, z_offset=0.11)
        pl.add_mesh(glyphs_from_layer(layer_gap), color="#2ca02c")
    if streamlines:
        _add_streamlines(pl, r_max, r_kpc, b_mu, pa_deg, grid_n=grid_n, z_plot=0.14)
    pl.add_text(
        f"{galaxy_name} — baryonic (blue) vs observed centripetal (red) accelerations{extra_block}\n"
        f"{_CAPTION_FOOTER}\n{_scenario_caption(summary)}",
        position="upper_left",
        font_size=9,
        color="white",
        viewport=True,
    )
    pl.camera_position = "xy"
    pl.enable_parallel_projection()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pl.screenshot(str(output_path))
    if export_html:
        html_path = output_path.with_suffix(".html")
        try:
            pl.export_html(str(html_path))
            try:
                sz = html_path.stat().st_size
                if sz > int(VIZ_DISPLAY.MUSEUM_HTML_WARN_BYTES):
                    logger.warning(
                        "Museum HTML for %s is large (%s bytes); browser load may be slow.",
                        galaxy_name,
                        sz,
                    )
            except OSError:
                pass
        except (OSError, AttributeError, ValueError, ImportError) as exc:
            logger.warning("HTML export failed for %s: %s", galaxy_name, exc)
    pl.close()


def run_cli(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Supreme-cohort 3D toy-disk renders from pipeline analysis CSV (optional PyVista).",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path.cwd(),
        help="Repository root (expects results/galaxy_analysis/...).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/graphs/viz_3d"),
        help="Directory for PNG and HTML exports.",
    )
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="Write PNG only (skip interactive Museum HTML export).",
    )
    args = parser.parse_args(argv)

    root = args.project_root.resolve()
    out_dir = args.output_dir
    if not out_dir.is_absolute():
        out_dir = (root / out_dir).resolve()

    from src.reporting.publication_bundle import cohort_csv_path, read_cohort_galaxy_names

    ev_path = cohort_csv_path(root, "scientific_evidence")
    if not ev_path.is_file():
        raise SystemExit(f"Missing supreme cohort table: {ev_path}")

    import pandas as pd

    ev_df = pd.read_csv(ev_path)
    morph_by_name: dict[str, str | None] = {}
    if "galaxy_name" in ev_df.columns and "morphology_class" in ev_df.columns:
        for _, row in ev_df.iterrows():
            gn = str(row.get("galaxy_name") or "").strip()
            if not gn:
                continue
            mc = row.get("morphology_class")
            morph_by_name[gn] = str(mc).strip() if mc is not None and str(mc).strip() else None

    names = read_cohort_galaxy_names(root, "scientific_evidence")
    written: list[str] = []
    for name in names:
        safe = _sanitize_filename(name)
        target = out_dir / f"{safe}_viz_3d.png"
        render_galaxy_png(
            root,
            name,
            target,
            morphology_class_override=morph_by_name.get(name),
            export_html=not bool(args.no_html),
        )
        print(f"Wrote {target}", flush=True)
        written.append(str(target.resolve()))
        if not args.no_html:
            hp = target.with_suffix(".html")
            if hp.is_file():
                written.append(str(hp.resolve()))
    if written:
        mp = write_viz_manifest(out_dir, written)
        print(f"Manifest {mp}", flush=True)


def main() -> None:
    run_cli()


if __name__ == "__main__":
    main()
