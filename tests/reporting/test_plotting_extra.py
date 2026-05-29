import pandas as pd

from src.models.galaxy import Galaxy
from src.reporting.plotting.plotter import GalaxyPlotter


def test_pearson_empty(tmp_path):
    path = tmp_path / "p.png"
    GalaxyPlotter().save_pearson_correlation_distribution([], path)
    assert path.exists()


def test_acceleration_sample(tmp_path):
    df = pd.DataFrame(
        {
            "radius_kpc": [1.0, 2.0],
            "acceleration_gap": [1e-11, -1e-11],
            "magnetic_pressure_acceleration": [2e-11, -2e-11],
        }
    )
    gal = Galaxy(name="X", rotation_dataframe=df)
    path = tmp_path / "a.png"
    GalaxyPlotter().save_acceleration_vs_radius_sample(gal, path)
    assert path.exists()
