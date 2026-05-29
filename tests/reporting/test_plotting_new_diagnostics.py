"""Tests for the new plotting helpers: error bars, RAR plot, and B-field histogram."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from src.models.galaxy import Galaxy
from src.reporting.plotting.plotter import GalaxyPlotter


def _make_galaxy_dataframe(*, with_systematic: bool, with_rar_columns: bool) -> pd.DataFrame:
    radii = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    df = pd.DataFrame(
        {
            "radius_kpc": radii,
            "observed_velocity_kms": np.array([60.0, 80.0, 95.0, 100.0, 105.0]),
            "baryonic_velocity_kms": np.array([55.0, 70.0, 80.0, 82.0, 83.0]),
            "gas_velocity_kms": np.array([30.0, 35.0, 38.0, 40.0, 41.0]),
            "disk_velocity_kms": np.array([45.0, 60.0, 70.0, 71.0, 72.0]),
            "bulge_velocity_kms": np.zeros_like(radii),
            "acceleration_gap": np.array([1e-11, 9e-12, 8e-12, 7e-12, 6e-12]),
            "magnetic_pressure_acceleration": np.array([1.2e-11, 1e-11, 9e-12, 7.5e-12, 6.5e-12]),
            "reliable_measurement": np.array([True, True, True, True, True]),
        }
    )
    if with_systematic:
        df["centripetal_acceleration_systematic"] = df["acceleration_gap"].abs() * 0.2
    if with_rar_columns:
        df["rar_predicted_acceleration"] = df["magnetic_pressure_acceleration"] * 0.95
        df["rar_residual_log10"] = np.array([0.0, -0.05, 0.05, -0.1, 0.02])
    return df


def _make_galaxy(df: pd.DataFrame) -> Galaxy:
    return Galaxy(name="TEST_GAL", distance_mpc=10.0, rotation_dataframe=df)


def test_force_plot_with_inclination_error_bars(tmp_path: Path):
    df = _make_galaxy_dataframe(with_systematic=True, with_rar_columns=False)
    galaxy = _make_galaxy(df)
    plotter = GalaxyPlotter()
    out = plotter.save_force_plot(
        galaxy,
        tmp_path,
        scenario_gradient_scale_kpc=0.1,
        scenario_ionization_fraction=0.01,
    )
    assert out.exists()
    assert out.stat().st_size > 1024


def test_force_plot_without_systematic_falls_back(tmp_path: Path):
    """Plot still renders when the engine never populated the systematic column."""
    df = _make_galaxy_dataframe(with_systematic=False, with_rar_columns=False)
    galaxy = _make_galaxy(df)
    out = GalaxyPlotter().save_force_plot(galaxy, tmp_path)
    assert out.exists()


def test_rar_plot_renders(tmp_path: Path):
    df = _make_galaxy_dataframe(with_systematic=True, with_rar_columns=True)
    galaxy = _make_galaxy(df)
    out = GalaxyPlotter().save_rar_plot(galaxy, tmp_path / "rar.png")
    assert out.exists()
    assert out.stat().st_size > 1024


def test_rar_plot_handles_empty_dataframe(tmp_path: Path):
    """Galaxy with no finite g_obs/g_bar still produces a valid (empty-state) figure."""
    df = pd.DataFrame(
        {
            "radius_kpc": [0.0, 0.0],
            "observed_velocity_kms": [0.0, 0.0],
            "baryonic_velocity_kms": [0.0, 0.0],
            "gas_velocity_kms": [0.0, 0.0],
            "disk_velocity_kms": [0.0, 0.0],
            "bulge_velocity_kms": [0.0, 0.0],
            "acceleration_gap": [np.nan, np.nan],
            "magnetic_pressure_acceleration": [np.nan, np.nan],
            "reliable_measurement": [False, False],
        }
    )
    out = GalaxyPlotter().save_rar_plot(_make_galaxy(df), tmp_path / "empty_rar.png")
    assert out.exists()


def test_magnetic_field_histogram_renders(tmp_path: Path):
    """Histogram with a band shade renders for typical μG values."""
    out = GalaxyPlotter().save_magnetic_field_histogram(
        [0.5, 1.2, 2.0, 3.5, 5.0, 8.0, 12.0, 25.0, 0.8, 1.5],
        tmp_path / "b_hist.png",
        plausibility_band=(0.1, 100.0),
    )
    assert out.exists()
    assert out.stat().st_size > 1024


def test_magnetic_field_histogram_handles_empty(tmp_path: Path):
    """Empty input still produces a valid empty-state PNG (no crash)."""
    out = GalaxyPlotter().save_magnetic_field_histogram([], tmp_path / "empty_hist.png")
    assert out.exists()


def test_plasma_beta_vs_coverage_plot(tmp_path: Path):
    rows = [
        {"plasma_beta_thermal_median": 0.5, "theoretical_gap_coverage_ratio": 1.2},
        {"plasma_beta_thermal_median": 10.0, "theoretical_gap_coverage_ratio": 0.4},
    ]
    out = GalaxyPlotter().save_plasma_beta_vs_coverage(rows, tmp_path / "pb.png")
    assert out.exists()


def test_magnetic_field_vs_galaxy_scale_plot(tmp_path: Path):
    """The Realistic Galactic Range band must render with a small reliable sample."""
    rows = [
        {"galaxy_name": "G1", "r_max_kpc": 5.0, "magnetic_field_microgauss": 1.5},
        {"galaxy_name": "G2", "r_max_kpc": 12.0, "magnetic_field_microgauss": 4.0},
        {"galaxy_name": "G3", "r_max_kpc": 25.0, "magnetic_field_microgauss": 0.6},
        {"galaxy_name": "G4", "r_max_kpc": 0.0, "magnetic_field_microgauss": 5.0},  # filtered
        {"galaxy_name": "G5", "r_max_kpc": 8.0, "magnetic_field_microgauss": None},  # filtered
    ]
    path = tmp_path / "b_vs_scale.png"
    GalaxyPlotter().save_magnetic_field_vs_galaxy_scale(rows, path)
    assert path.exists()
    assert path.stat().st_size > 1024


def test_magnetic_field_vs_galaxy_scale_handles_empty(tmp_path: Path):
    """Empty input still produces a valid empty-state PNG (no crash)."""
    path = tmp_path / "b_vs_scale_empty.png"
    GalaxyPlotter().save_magnetic_field_vs_galaxy_scale([], path)
    assert path.exists()


def test_pearson_distribution_renders_bands(tmp_path: Path):
    """Pearson histogram renders even without samples (just the bands and labels)."""
    path = tmp_path / "pearson.png"
    GalaxyPlotter().save_pearson_correlation_distribution([0.85, 0.6, 0.2, -0.3, 0.95], path)
    assert path.exists()
    assert path.stat().st_size > 1024


def test_morphology_comparison_renders_all_buckets(tmp_path: Path):
    """The generalized morphology bar chart must render every populated bucket and skip empties."""
    stats = {
        "Dwarf Irregulars (Im)": {"mean_pct": 95.4, "count": 4},
        "Late-type Spirals (Sc, Scd)": {"mean_pct": 110.2, "count": 6},
        "Photometric LSB (SBdisk < 100)": {"mean_pct": 88.0, "count": 3},
        "Lenticular/S0": {"mean_pct": 0.0, "count": 0},
    }
    path = tmp_path / "morph_bars.png"
    GalaxyPlotter().save_morphology_comparison(stats, path)
    assert path.exists()
    assert path.stat().st_size > 1024
    md_text = path.with_suffix(".md").read_text(encoding="utf-8")
    # Every bucket with count>0 should be cited in the companion markdown.
    assert "Dwarf Irregulars" in md_text
    assert "Photometric LSB" in md_text
    assert "Late-type Spirals" in md_text
    # Empty bucket must NOT mislead readers into thinking it has data.
    assert "Lenticular/S0" not in md_text


def test_morphology_comparison_handles_no_data(tmp_path: Path):
    """All-empty buckets still produce a valid empty-state PNG and companion file."""
    stats = {"Lenticular/S0": {"mean_pct": 0.0, "count": 0}}
    path = tmp_path / "morph_empty.png"
    GalaxyPlotter().save_morphology_comparison(stats, path)
    assert path.exists()
    assert path.with_suffix(".md").exists()


