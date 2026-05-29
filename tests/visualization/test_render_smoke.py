from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from src.visualization.io import load_analysis_dataframe


def test_load_analysis_dataframe_missing_raises(tmp_path: Path):
    (tmp_path / "results").mkdir()
    with pytest.raises(FileNotFoundError):
        load_analysis_dataframe(tmp_path, "MissingGalaxy")


def test_render_galaxy_png_smoke(tmp_path: Path):
    pv = pytest.importorskip("pyvista")
    _ = pv  # silence lint

    galaxy = "TestGal"
    tier = tmp_path / "results" / "galaxy_analysis" / "validation_reliable" / galaxy
    tier.mkdir(parents=True)
    rows = {
        "radius_kpc": [0.5, 1.0, 2.0],
        "baryonic_acceleration": [1e-11, 8e-12, 5e-12],
        "centripetal_acceleration": [2e-11, 1.5e-11, 9e-12],
        "magnetic_pressure_acceleration": [1e-12, 2e-12, 3e-12],
        "acceleration_gap": [1e-11, 7e-12, 4e-12],
        "effective_ionized_mass_density": [1e-24, 2e-24, 3e-24],
        "magnetic_field_microgauss": [2.0, 2.5, 3.0],
    }
    pd.DataFrame(rows).to_csv(tier / f"{galaxy}_analysis.csv", index=False)

    from src.visualization.render import render_galaxy_png

    out = tmp_path / "out" / "t.png"
    render_galaxy_png(tmp_path, galaxy, out, export_html=True)
    assert out.is_file()
    html = out.with_suffix(".html")
    if html.is_file():
        assert html.stat().st_size > 0

    out2 = tmp_path / "out" / "manual.png"
    render_galaxy_png(
        tmp_path,
        galaxy,
        out2,
        gap_glyphs=False,
        wim_rings=False,
        streamlines=False,
        export_html=False,
    )
    assert out2.is_file()
