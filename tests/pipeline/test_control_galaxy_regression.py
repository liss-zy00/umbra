"""Regression guard-rails for gold-standard SPARC control galaxies (project audit policy)."""

from __future__ import annotations

from src.pipeline.constants import DEFAULT_CONTROL_GALAXIES, E2E_AUDIT_CONTROL_GALAXIES


def test_default_control_preset_includes_e2e_audit_pair():
    """E2E audit policy galaxies must remain in the CLI ``default`` preset list."""
    names = {str(g).strip().upper() for g in DEFAULT_CONTROL_GALAXIES}
    for g in E2E_AUDIT_CONTROL_GALAXIES:
        assert g.upper() in names
