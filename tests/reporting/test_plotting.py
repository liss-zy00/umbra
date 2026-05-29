import pandas as pd

from src.models.galaxy import Galaxy
from src.reporting.plotting.plotter import GalaxyPlotter


def test_force_plot_writes_png(tmp_path):
    df = pd.DataFrame(
        {
            "radius_kpc": [1.0, 2.0],
            "acceleration_gap": [1e-11, 5e-12],
            "magnetic_pressure_acceleration": [1e-11, 4e-12],
        }
    )
    gal = Galaxy(name="P", rotation_dataframe=df)
    path = GalaxyPlotter().save_force_plot(gal, tmp_path, scenario_gradient_scale_kpc=0.1, scenario_ionization_fraction=0.01)
    assert path.suffix == ".png"
    assert path.exists()
