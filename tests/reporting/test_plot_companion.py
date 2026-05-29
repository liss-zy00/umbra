"""Tests for Markdown sidecars next to PNG figures."""
from __future__ import annotations

from pathlib import Path

from src.reporting.plotting.companion import write_plot_companion_markdown


def test_write_plot_companion_creates_md_next_to_png(tmp_path: Path):
    png = tmp_path / "figures" / "demo_plot.png"
    png.parent.mkdir(parents=True)
    png.write_bytes(b"\x89PNG\r\n\x1a\n")
    md = write_plot_companion_markdown(
        png,
        title="Demo plot",
        purpose="Unit test fixture.",
        data_shown="- Synthetic bytes only.",
        axes_and_units="- Not applicable.",
        how_to_read="N/A.",
    )
    assert md == png.with_suffix(".md")
    assert md.exists()
    text = md.read_text(encoding="utf-8")
    assert "Demo plot" in text
    assert "[demo_plot.png](./demo_plot.png)" in text
