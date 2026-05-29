"""Matplotlib style presets for publication-quality figures."""

from __future__ import annotations

from typing import Any


def apply_academic_mpl_style() -> None:
    """
    Serif-forward defaults that read well in PDF/thesis without requiring a LaTeX install.

    Uses DeJaVu Serif (bundled with Matplotlib) so CI and headless runs stay reproducible.
    Call once at the start of a batch plotting job.
    """
    import matplotlib.pyplot as plt

    rc: dict[str, Any] = {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "Nimbus Roman", "Computer Modern Roman", "serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "legend.fontsize": 9,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "figure.titlesize": 13,
        "figure.dpi": 100,
        "savefig.dpi": 160,
        "axes.grid": True,
        "grid.alpha": 0.25,
    }
    plt.rcParams.update(rc)
