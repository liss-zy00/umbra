import sys

import pytest

from src.pipeline.cli import parse_args
from src.pipeline.constants import DEFAULT_CONTROL_GALAXIES
from src.pipeline.runner import _resolve_target_set, run_pipeline
from src.providers.sparc import SparcMorphologyProvider


def test_run_pipeline_raises_when_no_rotmod_files(tmp_path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(sys, "argv", ["prog"])
    with pytest.raises(FileNotFoundError, match="No SPARC galaxy files"):
        run_pipeline(parse_args(), project_root=tmp_path)


def test_resolve_target_set_rejects_mixed_presets_and_names(tmp_path):
    sparc_dir = tmp_path / "rotation_curves"
    sparc_dir.mkdir(parents=True)
    g0 = DEFAULT_CONTROL_GALAXIES[0]
    (sparc_dir / f"{g0}_rotmod.dat").write_text("# stub\n")
    morphology = SparcMorphologyProvider(tmp_path / "missing.json")

    with pytest.raises(SystemExit, match="cannot mix presets"):
        _resolve_target_set(["all", g0], sparc_dir, morphology)


def test_resolve_target_set_rejects_two_presets(tmp_path):
    sparc_dir = tmp_path / "rotation_curves"
    sparc_dir.mkdir(parents=True)
    g0 = DEFAULT_CONTROL_GALAXIES[0]
    (sparc_dir / f"{g0}_rotmod.dat").write_text("# stub\n")
    morphology = SparcMorphologyProvider(tmp_path / "missing.json")

    with pytest.raises(SystemExit, match="at most one preset"):
        _resolve_target_set(["all", "default"], sparc_dir, morphology)


def test_resolve_target_set_explicit_name_filters_set(tmp_path):
    sparc_dir = tmp_path / "rotation_curves"
    sparc_dir.mkdir(parents=True)
    g0, g1 = DEFAULT_CONTROL_GALAXIES[0], DEFAULT_CONTROL_GALAXIES[1]
    (sparc_dir / f"{g0}_rotmod.dat").write_text("# stub\n")
    (sparc_dir / f"{g1}_rotmod.dat").write_text("# stub\n")
    morphology = SparcMorphologyProvider(tmp_path / "missing.json")

    selected = _resolve_target_set([g0], sparc_dir, morphology)
    assert [p.name for p in selected] == [f"{g0}_rotmod.dat"]


def test_resolve_target_set_preset_all_returns_full_set(tmp_path):
    sparc_dir = tmp_path / "rotation_curves"
    sparc_dir.mkdir(parents=True)
    g0, g1 = DEFAULT_CONTROL_GALAXIES[0], DEFAULT_CONTROL_GALAXIES[1]
    (sparc_dir / f"{g0}_rotmod.dat").write_text("# stub\n")
    (sparc_dir / f"{g1}_rotmod.dat").write_text("# stub\n")
    morphology = SparcMorphologyProvider(tmp_path / "missing.json")

    selected = _resolve_target_set(["all"], sparc_dir, morphology)
    assert len(selected) == 2
