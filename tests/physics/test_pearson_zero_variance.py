"""Regression test: ``pearson_correlation_coefficient`` returns ``None`` instead of crashing
on degenerate (constant) inputs."""

from __future__ import annotations

import numpy as np

from src.physics.diagnostics.statistics import pearson_correlation_coefficient


def test_pearson_returns_none_for_constant_array_a():
    a = np.array([1.0, 1.0, 1.0])
    b = np.array([1.0, 2.0, 3.0])
    assert pearson_correlation_coefficient(a, b) is None


def test_pearson_returns_none_for_constant_array_b():
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([5.0, 5.0, 5.0])
    assert pearson_correlation_coefficient(a, b) is None


def test_pearson_returns_none_for_too_few_valid_points():
    a = np.array([1.0, np.nan, np.nan])
    b = np.array([2.0, 3.0, 4.0])
    assert pearson_correlation_coefficient(a, b) is None


def test_pearson_returns_one_for_perfectly_aligned_arrays():
    a = np.linspace(0.0, 1.0, 10)
    b = 2.0 * a + 3.0
    assert np.isclose(pearson_correlation_coefficient(a, b), 1.0, atol=1e-12)
