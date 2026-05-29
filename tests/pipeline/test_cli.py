import sys

import pytest

from src.pipeline.cli import parse_args
from src.pipeline.constants import DEFAULT_CONTROL_GALAXIES


def test_parse_args_default_is_all(monkeypatch: pytest.MonkeyPatch):
    """Running without any flag must default to the full SPARC catalogue."""
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.setattr(sys, "argv", ["prog"])
    ns = parse_args()
    assert ns.galaxies == ["all"]


def test_parse_args_galaxies_default_preset(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.setattr(sys, "argv", ["prog", "--galaxies", "default"])
    ns = parse_args()
    assert ns.galaxies == ["default"]


def test_parse_args_galaxies_explicit_list(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("CI", raising=False)
    g0, g1 = DEFAULT_CONTROL_GALAXIES[0], DEFAULT_CONTROL_GALAXIES[1]
    monkeypatch.setattr(sys, "argv", ["prog", "--galaxies", g0, g1])
    ns = parse_args()
    assert ns.galaxies == [g0, g1]


def test_parse_args_rejects_unknown_flag(monkeypatch: pytest.MonkeyPatch):
    """Argparse must fail fast on unsupported options (no silent ignores)."""
    monkeypatch.setattr(sys, "argv", ["prog", "--not-a-supported-cli-option"])
    with pytest.raises(SystemExit):
        parse_args()


def test_parse_args_logs_flag(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.setattr(sys, "argv", ["prog", "--galaxies", "default", "--logs"])
    ns = parse_args()
    assert ns.logs is True


def test_parse_args_logs_default_false(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.setattr(sys, "argv", ["prog", "--galaxies", "default"])
    ns = parse_args()
    assert ns.logs is False


def test_parse_args_ci_does_not_mutate_flags(monkeypatch: pytest.MonkeyPatch):
    """CI must not override Rich console or logging flags (no hidden ``plain`` mode)."""
    monkeypatch.setenv("CI", "true")
    monkeypatch.setattr(sys, "argv", ["prog", "--galaxies", "default"])
    ns = parse_args()
    assert ns.logs is False


def test_parse_args_compute_only_and_viz_only_mutex(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.setattr(sys, "argv", ["prog", "--compute-only", "--viz-only"])
    with pytest.raises(SystemExit):
        parse_args()


def test_parse_args_compute_only_sets_flag(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.setattr(sys, "argv", ["prog", "--galaxies", "default", "--compute-only"])
    ns = parse_args()
    assert ns.compute_only is True
    assert ns.viz_only is False


def test_parse_args_viz_only_sets_flag(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.setattr(sys, "argv", ["prog", "--viz-only"])
    ns = parse_args()
    assert ns.viz_only is True
    assert ns.compute_only is False
