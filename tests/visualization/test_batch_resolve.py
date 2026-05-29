from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.pipeline.constants import (
    RESULTS_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_RELIABLE_FILENAME,
)
from src.visualization.batch import expand_cli_galaxy_names


def test_expand_cli_galaxy_names_auto_reads_reliable(tmp_path: Path):
    vdir = tmp_path / RESULTS_DIRNAME / RESULTS_STEP_VALIDATION_DIRNAME
    vdir.mkdir(parents=True)
    df = pd.DataFrame({"galaxy_name": ["A", "B", "C"], "validation_is_reliable": [True, True, True]})
    df.to_csv(vdir / VALIDATION_RELIABLE_FILENAME, index=False)

    names = expand_cli_galaxy_names(tmp_path, ["auto"], tier="reliable", max_n=2)
    assert names == ["A", "B"]


def test_expand_cli_galaxy_names_explicit_passthrough(tmp_path: Path):
    names = expand_cli_galaxy_names(tmp_path, ["X", " Y "], tier="reliable", max_n=99)
    assert names == ["X", "Y"]
