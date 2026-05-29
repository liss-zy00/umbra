"""Sidecar Markdown for each figure: same basename as the PNG, ``.md`` extension.

Keeps long-form documentation out of the figure canvas while staying co-located with
the artifact for sanity checks and publication workflows.
"""

from __future__ import annotations

from pathlib import Path


def write_plot_companion_markdown(
    figure_png_path: Path,
    *,
    title: str,
    purpose: str,
    data_shown: str,
    axes_and_units: str,
    how_to_read: str,
) -> Path:
    """
    Write ``<stem>.md`` next to ``<stem>.png`` describing purpose, data, axes, and interpretation.

    The Markdown links to the PNG in the same directory so viewers in GitHub / VS Code resolve
    the pair without path gymnastics.
    """
    figure_png_path = Path(figure_png_path)
    figure_png_path.parent.mkdir(parents=True, exist_ok=True)
    md_path = figure_png_path.with_suffix(".md")
    png_name = figure_png_path.name
    lines = [
        f"# {title}",
        "",
        f"**Figure:** [{png_name}](./{png_name})",
        "",
        "## Purpose",
        "",
        purpose.strip(),
        "",
        "## Data shown",
        "",
        data_shown.strip(),
        "",
        "## Axes and units",
        "",
        axes_and_units.strip(),
        "",
        "## How to read this figure",
        "",
        how_to_read.strip(),
        "",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")
    return md_path
