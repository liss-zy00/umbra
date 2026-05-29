"""
Construct PyVista meshes for toy disk + acceleration glyphs.

Geometry is a communication aid: radii use ``radius_kpc`` as display coordinates
(kpc treated as numeric length in the scene). Accelerations remain SI magnitudes
for relative scaling of arrow glyphs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ArrowLayer:
    """Point cloud with one vector per row for ``PolyData.glyph``."""

    points_xyz: np.ndarray
    vectors: np.ndarray


def reference_disk(r_max_kpc: float, *, n_angular: int = 96) -> "pv.PolyData":
    import pyvista as pv

    return pv.Disc(inner=0.0, outer=float(r_max_kpc), r_res=2, c_res=n_angular, normal=(0.0, 0.0, 1.0))


def _subsample_indices(n_rows: int, max_points: int) -> np.ndarray:
    cap = max(2, int(max_points))
    if n_rows <= cap:
        return np.arange(n_rows, dtype=int)
    return np.unique(np.linspace(0, n_rows - 1, num=cap, dtype=int))


def layer_from_radial_profile(
    radius_kpc: np.ndarray,
    acceleration_m_s2: np.ndarray,
    *,
    max_arrows: int = 40,
    z_offset: float = 0.02,
    r_max_kpc: float | None = None,
) -> ArrowLayer:
    """
    Place arrows along +x at ``(r_kpc, 0, z_offset)`` with vector ``(g, 0, 0)``.

    Vector magnitudes are rescaled into *display kpc* using ``0.35 * r_max_kpc``
    times ``g / max(|g|)`` so glyphs stay within the disk for typical SPARC ranges.
    """
    r = np.asarray(radius_kpc, dtype=float)
    g = np.asarray(acceleration_m_s2, dtype=float)
    g = np.nan_to_num(g, nan=0.0, posinf=0.0, neginf=0.0)
    r_cap = float(np.nanmax(r)) if r_max_kpc is None else float(r_max_kpc)
    gabs = float(np.nanmax(np.abs(g))) if g.size else 0.0
    denom = gabs if gabs > 0.0 else 1.0
    disp_scale = 0.35 * max(r_cap, 1e-6)
    idx = _subsample_indices(int(r.size), max_arrows)
    r_s = r[idx]
    g_s = (g[idx] / denom) * disp_scale
    pts = np.column_stack([r_s, np.zeros_like(r_s), np.full_like(r_s, z_offset)])
    vec = np.column_stack([g_s, np.zeros_like(g_s), np.zeros_like(g_s)])
    return ArrowLayer(points_xyz=pts, vectors=vec)


def glyphs_from_layer(layer: ArrowLayer, *, factor: float = 1.0) -> "pv.PolyData":
    import pyvista as pv

    pdata = pv.PolyData(layer.points_xyz)
    pdata["vectors"] = layer.vectors
    return pdata.glyph(orient="vectors", scale="vectors", factor=float(factor), geom=pv.Arrow())


def scalar_disk_colors(
    disk: "pv.PolyData",
    radius_kpc: np.ndarray,
    scalar: np.ndarray,
    *,
    name: str,
) -> "pv.PolyData":
    """Attach cell/point data sampled radially for disk coloring (Panel B context)."""
    r = np.asarray(radius_kpc, dtype=float)
    s = np.nan_to_num(np.asarray(scalar, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    # Disk cell centers approximate mid-radius per ring — use point data on disc points.
    disk = disk.copy(deep=True)
    x = disk.points[:, 0]
    y = disk.points[:, 1]
    rad = np.hypot(x, y)
    vals = np.interp(rad, r, s, left=s[0] if s.size else 0.0, right=s[-1] if s.size else 0.0)
    disk.point_data[name] = vals
    return disk


def wim_density_annulus_meshes(
    radii_kpc: np.ndarray,
    rho_kg_m3: np.ndarray,
    *,
    n_bins: int,
    z_plot: float = 0.03,
) -> list["pv.PolyData"]:
    """
    Concentric annular discs coloured by shell-median ``effective_ionized_mass_density``.

    This is a lightweight pseudo-3D cue (not a self-consistent volume render).
    """
    import pyvista as pv

    r = np.asarray(radii_kpc, dtype=float)
    rho = np.asarray(rho_kg_m3, dtype=float)
    meshes: list[pv.PolyData] = []
    if r.size == 0 or n_bins < 2:
        return meshes
    r_max = float(np.nanmax(r))
    if not np.isfinite(r_max) or r_max <= 0.0:
        return meshes
    edges = np.linspace(0.0, r_max, int(n_bins) + 1)
    for i in range(int(n_bins)):
        r0, r1 = float(edges[i]), float(edges[i + 1])
        if r1 <= r0:
            continue
        inner = r0 if r0 > 0.0 else max(r1 * 1.0e-6, r_max * 1.0e-9)
        if i < int(n_bins) - 1:
            mask = (r >= r0) & (r < r1)
        else:
            mask = (r >= r0) & (r <= r1)
        if not np.any(mask):
            vals_bin = float("nan")
        else:
            vals_bin = float(np.nanmedian(rho[mask]))
        disc = pv.Disc(inner=inner, outer=r1, r_res=2, c_res=72, normal=(0.0, 0.0, 1.0))
        disc = disc.translate((0.0, 0.0, float(z_plot)), inplace=False)
        disc = disc.copy(deep=True)
        logv = float(np.log10(max(vals_bin, 1.0e-40))) if np.isfinite(vals_bin) else 0.0
        disc.point_data["log10_rho_ion_kg_m3"] = np.full(disc.n_points, logv, dtype=float)
        meshes.append(disc)
    return meshes


def illustrative_b_streamlines(
    r_max_kpc: float,
    r_kpc: np.ndarray,
    b_microgauss: np.ndarray,
    pa_deg: float | None,
    *,
    grid_n: int,
    z_plot: float,
    max_time: float | None = None,
) -> "pv.PolyData | None":
    """
    In-plane toy azimuthal field scaled by ``magnetic_field_microgauss(r)`` with optional PA offset.

    Not a Planck vector map — pedagogical streamlines only (see caption strings in ``render.py``).
    """
    import pyvista as pv

    r_line = np.asarray(r_kpc, dtype=float)
    b_line = np.nan_to_num(np.asarray(b_microgauss, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    if r_line.size == 0:
        return None
    r_cap = float(r_max_kpc) if np.isfinite(float(r_max_kpc)) else float(np.nanmax(r_line))
    if not np.isfinite(r_cap) or r_cap <= 0.0:
        return None
    half = r_cap * 1.05
    gn = max(8, int(grid_n))
    xs = np.linspace(-half, half, gn)
    ys = np.linspace(-half, half, gn)
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    zz = np.zeros_like(xx, dtype=float) + float(z_plot)
    grid = pv.StructuredGrid(xx, yy, zz)
    rr = np.hypot(xx, yy)
    phi = np.arctan2(yy, xx)
    b_flat = np.interp(
        rr.ravel(),
        r_line,
        b_line,
        left=float(b_line[0]) if b_line.size else 0.0,
        right=float(b_line[-1]) if b_line.size else 0.0,
    )
    bmu = b_flat.reshape(rr.shape)
    bmu = np.nan_to_num(bmu, nan=0.0, posinf=0.0, neginf=0.0)
    pa = np.deg2rad(float(pa_deg)) if pa_deg is not None else 0.0
    cph, sph = np.cos(phi + pa), np.sin(phi + pa)
    vx = -bmu * sph
    vy = bmu * cph
    vz = np.zeros_like(vx)
    grid.point_data["vectors"] = np.column_stack([vx.ravel(), vy.ravel(), vz.ravel()])
    mtime = float(max_time) if max_time is not None else max(0.45 * r_cap, 0.5)
    try:
        try:
            lines = grid.streamlines(
                vectors="vectors",
                n_points=120,
                max_length=float(mtime),
            )
        except TypeError:
            lines = grid.streamlines(
                vectors="vectors",
                n_points=120,
                max_time=float(mtime),
            )
    except (AttributeError, ValueError):
        return None
    if lines is None or lines.n_points == 0:
        return None
    return lines
