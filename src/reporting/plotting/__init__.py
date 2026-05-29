"""
Matplotlib output for the pipeline.

* :mod:`plotter` — :class:`GalaxyPlotter`: per-galaxy figures (acceleration comparison,
  RAR plot, sample) and global summaries (scatter, scenario heatmap, morphology
  comparison, B-field histogram, plausibility scatter, plasma-β vs coverage, Pearson
  distribution).
* :mod:`companion` — :func:`write_plot_companion_markdown`: writes a sidecar ``.md``
  description file next to every figure so reviewers always see purpose, axes, and
  reading guide alongside the PNG.
"""

from __future__ import annotations

from .companion import write_plot_companion_markdown
from .plotter import GalaxyPlotter

__all__ = ["GalaxyPlotter", "write_plot_companion_markdown"]
