from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd

from src.models.galaxy import Galaxy
from src.pipeline.per_galaxy.galaxy_processor import MutablePipelineState, PipelineDependencies, process_sparc_file
from src.physics import GRADIENT_SCALES_KPC, IONIZATION_FRACTIONS


def test_process_sparc_file_skips_without_gas(tmp_path):
    rot = tmp_path / "EMPTY_rotmod.dat"
    rot.write_text("# Distance = 1 Mpc\n1.0 10 1 nan nan nan\n", encoding="utf-8")

    deps = MagicMock()
    deps.sparc.load_rotation_curve.return_value = (
        pd.DataFrame(
            {
                "radius_kpc": [1.0],
                "observed_velocity_kms": [10.0],
                "gas_velocity_kms": [np.nan],
                "disk_velocity_kms": [0.0],
                "bulge_velocity_kms": [0.0],
                "baryonic_velocity_kms": [0.0],
            }
        ),
        1.0,
    )

    scenario_rows = {(s, i): [] for s in GRADIENT_SCALES_KPC for i in IONIZATION_FRACTIONS}
    state = MutablePipelineState(
        scenario_rows=scenario_rows,
        scenario_success_counts=np.zeros((len(GRADIENT_SCALES_KPC), len(IONIZATION_FRACTIONS)), dtype=int),
    )

    process_sparc_file(deps, state, rot)
    deps.ned.fetch_metadata.assert_not_called()
