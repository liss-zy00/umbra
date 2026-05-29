"""Tests for the per-galaxy multi-source provenance block."""

from __future__ import annotations

import math

from src.pipeline.constants import DEFAULT_CONTROL_GALAXIES
from src.reporting.galaxy_log import (
    VALIDATION_PROVENANCE_STUB_LINE,
    emit_galaxy_provenance_block,
    emit_validation_grid_batch_header,
    emit_validation_grid_progress_finalize,
    emit_validation_grid_progress_rewrite,
    emit_validation_provenance_stub,
    format_completed_galaxy_grid,
    format_galaxy_provenance_block,
    format_phase3_graphs_verbose_block,
    format_phase4_analysis_verbose_block,
    format_validation_grid_batch_block,
    format_validation_phase_scenario_legend,
    format_validation_provenance_stub,
)


def test_format_full_block_validation_phase_fixed_width():
    g = DEFAULT_CONTROL_GALAXIES[1]
    width = 88
    block = format_galaxy_provenance_block(
        g,
        sparc={"source_file": f"{g}_rotmod.dat", "distance_mpc": 14.1, "n_shells": 21},
        morph={
            "morphology_class": "Late-type Spirals (Sc, Scd)",
            "hubble_type": 4,
            "inclination_deg": 68.0,
            "vflat_kms": 283.0,
        },
        ned={
            "status": "OK",
            "right_ascension_degrees": 140.5104,
            "declination_degrees": 50.9761,
            "redshift": 0.00211,
            "object_type": "G",
        },
        planck={"q_center": 2.78e-05, "u_center": 7.86e-05, "pixel_scale_kpc": 14.7},
        pipeline_phase="validation",
        rule_width=width,
        edd={
            "edd_tables_loaded": True,
            "pgc": 7004,
            "applied": True,
            "tier": 2,
            "source_table": "edd_cf4_DM",
            "distance_mpc": 14.1,
            "distance_uncertainty_mpc": 0.12,
        },
    )

    lines = block.splitlines()
    assert lines[0].startswith("    \u2500\u2500 ")
    assert f" {g} " in lines[0]
    assert " [validation]" not in lines[0]
    assert len(lines[0]) == width
    assert " SPARC  │ " in lines[1]
    assert "\u202f\u00b7\u202f" in lines[1]
    assert "dist: 14.10 Mpc" in lines[1]
    assert " EDD  │ " in lines[3] or "   EDD  │ " in lines[3]
    assert "PLANCK" in lines[5]
    assert "scale: 14.70 kpc/px" in lines[5]


def test_format_block_discovery_header_tag():
    width = 88
    block = format_galaxy_provenance_block("F561-1", pipeline_phase="discovery", rule_width=width)
    first = block.splitlines()[0]
    assert " [discovery]" not in first
    assert len(first) == width


def test_format_validation_stub_two_lines():
    width = 88
    stub = format_validation_provenance_stub("DDO154", rule_width=width)
    lines = stub.splitlines()
    assert len(lines) == 2
    assert " [validation]" not in lines[0]
    assert lines[1] == "    " + VALIDATION_PROVENANCE_STUB_LINE


def test_format_block_missing_values_render_as_n_slash_d():
    block = format_galaxy_provenance_block(
        "Phantom",
        sparc={"source_file": None, "distance_mpc": None, "n_shells": None},
        morph={
            "morphology_class": None,
            "hubble_type": None,
            "inclination_deg": float("nan"),
            "vflat_kms": None,
        },
        ned={
            "status": None,
            "right_ascension_degrees": None,
            "declination_degrees": None,
            "redshift": None,
            "object_type": None,
        },
        planck={"q_center": math.nan, "u_center": None, "pixel_scale_kpc": None},
        rule_width=80,
    )
    assert "None" not in block
    assert "nan" not in block.lower()
    assert "n/d" in block


def test_pipe_rows_have_provider_pipe():
    block = format_galaxy_provenance_block("X", rule_width=72, edd={"edd_tables_loaded": False})
    data_lines = block.splitlines()[1:]
    assert len(data_lines) == 5
    for line in data_lines:
        assert "  │ " in line
        if "   EDD  │ " in line or " EDD  │ " in line:
            assert "catalogues not loaded" in line or "\u202f\u00b7\u202f" in line
        else:
            assert "\u202f\u00b7\u202f" in line


def test_emit_skips_when_disabled(capsys):
    emit_galaxy_provenance_block(
        "Silent",
        sparc={"source_file": "x.dat", "distance_mpc": 1.0, "n_shells": 1},
        enabled=False,
    )
    assert capsys.readouterr().out == ""


def test_emit_writes_discovery_block(capsys):
    emit_galaxy_provenance_block(
        "TestGal",
        sparc={"source_file": "TestGal_rotmod.dat", "distance_mpc": 5.0, "n_shells": 7},
        ned={"status": "OK_NED_VIA_PATTERN_RULE", "right_ascension_degrees": 10.0, "declination_degrees": 20.0},
        pipeline_phase="discovery",
        rule_width=88,
    )
    captured = capsys.readouterr().out
    assert "[discovery]" not in captured
    assert "EDD" in captured
    assert " SPARC  │ " in captured
    assert "status: OK_NED_VIA_PATTERN_RULE" in captured


def test_emit_validation_stub(capsys):
    emit_validation_provenance_stub("UGCA281", enabled=True, rule_width=88)
    out = capsys.readouterr().out
    assert "UGCA281" in out and "\u2500\u2500 UGCA281" in out
    assert "    " + VALIDATION_PROVENANCE_STUB_LINE in out


def test_format_validation_grid_batch_block_counts():
    from src.physics import GRADIENT_SCALES_KPC, IONIZATION_FRACTIONS

    n_sc = len(GRADIENT_SCALES_KPC) * len(IONIZATION_FRACTIONS)
    pairs = [(float(s), float(i)) for s in GRADIENT_SCALES_KPC for i in IONIZATION_FRACTIONS]
    block = format_validation_grid_batch_block(
        n_cohort_galaxies=3,
        n_discovery_reliable=158,
        n_scenario_combinations=n_sc,
        scenario_pairs=pairs,
        rule_width=88,
    )
    lines = block.splitlines()
    assert "Validation (multi-scenario matrix)" in lines[0]
    assert lines[0].endswith(" [validation]")
    assert "3 galaxies" in lines[1] and "158 discovery-reliable" in lines[1]
    assert f"{n_sc} combinations/galaxy" in lines[2]
    assert lines[-1].startswith("→ Total evaluations")
    assert str(3 * n_sc) in lines[-1]
    scenario_body = lines[3:-1]
    assert len(scenario_body) == n_sc
    for sl in scenario_body:
        assert sl.startswith("  → ")
        assert "Δx=" in sl and "kpc" in sl and "f_ion=" in sl


def test_emit_validation_grid_batch_header_and_progress(capsys):
    from src.physics import GRADIENT_SCALES_KPC, IONIZATION_FRACTIONS

    pairs = [(float(s), float(i)) for s in GRADIENT_SCALES_KPC for i in IONIZATION_FRACTIONS]
    n_sc = len(pairs)
    emit_validation_grid_batch_header(
        n_cohort_galaxies=2,
        n_discovery_reliable=10,
        n_scenario_combinations=n_sc,
        scenario_pairs=pairs,
        enabled=True,
        rule_width=88,
    )
    emit_validation_grid_progress_rewrite(["A"], total=2, enabled=True)
    emit_validation_grid_progress_rewrite(["A", "B"], total=2, enabled=True)
    emit_validation_grid_progress_finalize(enabled=True, completed_names=["A", "B"], rule_width=88)
    out = capsys.readouterr().out
    assert "Validation (multi-scenario matrix)" in out
    assert "completed (2/2)" not in out
    assert "Galaxies processed (2)" in out
    assert "rows=" in out and "columns=" in out
    assert "A" in out and "B" in out


def test_format_completed_galaxy_grid_all_names_visible():
    names = ["DDO154", "NGC3741", "UGCA281", "X", "Y", "Z"]
    text = format_completed_galaxy_grid(names, rule_width=120, target_rows=7)
    for n in names:
        assert n in text
    assert "        Galaxies processed" in text
    assert "rows=" in text and "columns=" in text
    assert text.endswith("\n\n")


def test_format_completed_galaxy_grid_seven_rows_when_width_allows():
    names = [f"G{i:03d}" for i in range(158)]
    text = format_completed_galaxy_grid(names, rule_width=320, target_rows=7)
    assert "rows=7" in text
    assert "columns=23" in text


def test_format_validation_phase_scenario_legend_matches_matrix_size():
    from src.physics import GRADIENT_SCALES_KPC, IONIZATION_FRACTIONS

    pairs = [(float(s), float(i)) for s in GRADIENT_SCALES_KPC for i in IONIZATION_FRACTIONS]
    text = format_validation_phase_scenario_legend(pairs, n_galaxies=158)
    n_sc = len(pairs)
    assert "GRADIENT_SCALES_KPC × IONIZATION_FRACTIONS (nested);" in text
    assert "per file: scenarios 1→" in text
    assert text.count("Scenario ") == n_sc
    assert "    →" in text
    assert "→ Scenario 1" in text
    assert f"Scenario {n_sc}" in text
    assert "Scenario 1/" not in text
    assert "GRADIENT_SCALES_KPC" in text and "IONIZATION_FRACTIONS" in text
    assert "next galaxy" in text


def test_format_phase3_graphs_verbose_block_mentions_outputs():
    t = format_phase3_graphs_verbose_block(n_targets=12, graphs_rel="results/graphs")
    assert "galaxy_analysis" in t
    assert "12 names" in t
    assert "matplotlib" in t
    assert "Phase 3" not in t


def test_format_phase4_analysis_verbose_block_mentions_audit():
    t = format_phase4_analysis_verbose_block(audit_rel="results/step_validation/x.md", n_summary_scan=40)
    assert "40 *_summary.json" in t
    assert "x.md" in t
    assert "Phase 4" not in t
