import numpy as np
import pandas as pd

from src.models.galaxy import Galaxy
from src.physics import GRADIENT_SCALES_KPC, IONIZATION_FRACTIONS, PhysicsEngine
from src.pipeline.per_galaxy.scenario_grid_runner import run_scenario_grid


def test_scenario_grid_runs_without_crash():
    df = pd.DataFrame(
        {
            "radius_kpc": [1.0],
            "observed_velocity_kms": [80.0],
            "baryonic_velocity_kms": [50.0],
            "gas_velocity_kms": [40.0],
        }
    )
    gal = Galaxy(name="Mini", distance_mpc=10.0, rotation_dataframe=df)
    PhysicsEngine(gal).calculate_force_gap()
    base_df = gal.rotation_dataframe.copy()

    scenario_rows = {(s, i): [] for s in GRADIENT_SCALES_KPC for i in IONIZATION_FRACTIONS}
    counts = np.zeros((len(GRADIENT_SCALES_KPC), len(IONIZATION_FRACTIONS)), dtype=int)

    coords = Galaxy(name="Mini", distance_mpc=10.0, right_ascension_degrees=180.0, declination_degrees=45.0)
    metrics = run_scenario_grid(
        galaxy_name="Mini",
        base_dataframe=base_df,
        galaxy_coords=coords,
        stokes_summary={"stokes_amp_center": 1e-4, "stokes_amp_neighbors_mean": 1e-4, "status": "ok"},
        fits_summary={"stokes_amp_reference": 1e-4},
        morphology_class="Im",
        scenario_rows=scenario_rows,
        scenario_success_counts=counts,
    )
    assert len(metrics) == len(GRADIENT_SCALES_KPC) * len(IONIZATION_FRACTIONS)
    assert sum(len(v) for v in scenario_rows.values()) == len(metrics)
