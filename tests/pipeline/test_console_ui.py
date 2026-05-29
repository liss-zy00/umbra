import argparse
from pathlib import Path

from src.pipeline.console_ui import (
    PIPELINE_PHASE_1_DISCOVERY_TITLE,
    PlainPipelineConsole,
    RichPipelineConsole,
    build_pipeline_console,
    format_pipeline_summary_section,
    format_processing_progress_line,
    pipeline_session_start_block,
)


def test_build_console_always_rich():
    for ns in (
        argparse.Namespace(logs=False, galaxies=["default"]),
        argparse.Namespace(logs=True, galaxies=["all"]),
    ):
        c = build_pipeline_console(ns)
        assert isinstance(c, RichPipelineConsole)


def test_plain_console_class_exists_for_injection():
    """Tests may inject ``PlainPipelineConsole()`` into ``run_pipeline(..., console=...)``."""
    assert isinstance(PlainPipelineConsole(), PlainPipelineConsole)


def test_plain_phase_start_end_use_rule_banners(capsys):
    c = PlainPipelineConsole()
    c.phase_start(PIPELINE_PHASE_1_DISCOVERY_TITLE)
    c.phase_end(
        PIPELINE_PHASE_1_DISCOVERY_TITLE,
        detail="First overview sentence · Second fact · Third metric",
    )
    out = capsys.readouterr().out
    assert "[phase]" in out
    assert "Starting" not in out
    assert out.count("── Phase 1 — Discovery") == 1
    assert "First overview sentence" in out
    assert "Second fact" in out


def test_format_processing_progress_line_shape():
    line = format_processing_progress_line(51, 158, elapsed_s=3.2, terminal_cols=100)
    assert line.startswith("    Processing ")
    assert len(line) <= 99
    assert "%" in line
    assert "(51/158)" in line


def test_format_pipeline_summary_section_uses_indented_rule_banner():
    block_plain = format_pipeline_summary_section(
        "Results — Summary",
        "stats",
        ["175 SPARC  |  158 grid", "28 reliable  |  04 evidence"],
        rule_width=88,
    )
    assert block_plain.startswith(" ── Results — Summary")
    assert "[stats]" in block_plain
    assert "\u2192" not in block_plain
    assert "    175 SPARC" in block_plain

    block_arrow = format_pipeline_summary_section(
        "Results — Outputs",
        "tree",
        ["results/", "galaxy_analysis  |  tables"],
        rule_width=88,
        use_detail_arrow=True,
    )
    assert block_arrow.startswith(" ── Results — Outputs")
    assert "[tree]" in block_arrow
    assert "    \u2192 results/" in block_arrow
    assert "    \u2192 galaxy_analysis  |  tables" in block_arrow


def test_closing_summary_puts_finish_banner_after_sections():
    from src.pipeline.runner import _short_pipeline_closing_lines

    text = "\n".join(
        _short_pipeline_closing_lines(
            discovery_only=False,
            processed_root=Path("results"),
            repo_root=Path("."),
            n_sparc_targets=1,
            n_grid_galaxies=1,
            n_direct_evidence=0,
            n_validation_reliable=0,
            rule_width=88,
        )
    )
    assert text.index("Results — Outputs") < text.index("Pipeline Finished Successfully")
    assert text.index("Results — Summary") < text.index("Pipeline Finished Successfully")
    assert text.index("Results — Documentation") < text.index("Pipeline Finished Successfully")
    tail_lines = [ln for ln in text.splitlines() if ln.strip()]
    assert tail_lines[-1].startswith("\u2500\u2500 Pipeline Finished Successfully")


def test_pipeline_session_block_contains_session_bracket():
    block = pipeline_session_start_block(rule_width=88)
    assert "[session]" in block
    assert "E-G model pipeline" in block
    assert "Ingests SPARC rotation curves" in block
    assert "Runs the physics engine" in block
    assert "step_discovery" in block
