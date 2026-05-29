"""
Post-engine analytics and numerical helpers.

Concerns living here:

* :mod:`plasma` — plasma diagnostics (Plasma-β, Alfvén speed), magnetic-field plausibility
  band check, and the per-shell reliability mask builder.
* :mod:`rar` — Radial Acceleration Relation (McGaugh+2016) prediction and per-shell
  log-residuals plus the per-galaxy median |residual| statistic.
* :mod:`btfr` — Baryonic Tully–Fisher Relation (McGaugh+2012) ridge and per-galaxy
  distance-to-ridge residual.
* :mod:`tension` — Toy logarithmic-spiral magnetic-tension acceleration ``(B·∇)B/μ₀``.
* :mod:`uncertainty` — Monte-Carlo propagation of the ``(D, i, V_flat)`` band onto the
  gap-coverage ratio.
* :mod:`anchor_sensitivity` — CMB vs turbulent-equipartition anchor sweep on the
  Stokes → Tesla calibration constant.
* :mod:`resolution` — instrumental beam-smearing diagnostic ``k = D_beam / L_grad`` that
  quantifies how much the finite Planck beam dilutes the reconstructed magnetic gradient.
* :mod:`statistics` — small numerical helpers (Pearson, edge-radial mean) shared between
  the engine, the metrics builder, and the reporting layer.

The names are re-exported from the parent :mod:`src.physics` package; importers should
prefer that public surface over deep paths into this subpackage.
"""

from __future__ import annotations

from .anchor_sensitivity import (
    ANCHORS,
    anchor_scale_factor,
    build_anchor_sensitivity_rows,
    sweep_coverage_ratio,
)
from .btfr import (
    BTFRRidgeFit,
    MCGAUGH_2012_RIDGE,
    baryonic_mass_solar,
    btfr_residual_dex,
    is_btfr_outlier,
)
from .plasma import (
    TESLA_PER_MICROGAUSS,
    MagneticFieldPlausibilityCheck,
    PlasmaDiagnostics,
    ReliabilityMaskBuilder,
)
from .rar import (
    rar_median_abs_residual_dex,
    rar_predicted_observed_acceleration,
    rar_residual_log10,
)
from .resolution import (
    beam_fwhm_radians,
    beam_physical_size_kpc,
    resolution_correction_factor,
    vectorised_resolution_correction,
)
from .statistics import edge_radial_mean, pearson_correlation_coefficient
from .tension import (
    curvature_radius_m,
    supports_logarithmic_spiral,
    tension_acceleration_radial,
)
from .uncertainty import default_ratio_scaling, propagate_gap_coverage_band

__all__ = [
    "ANCHORS",
    "BTFRRidgeFit",
    "MCGAUGH_2012_RIDGE",
    "TESLA_PER_MICROGAUSS",
    "MagneticFieldPlausibilityCheck",
    "PlasmaDiagnostics",
    "ReliabilityMaskBuilder",
    "anchor_scale_factor",
    "baryonic_mass_solar",
    "beam_fwhm_radians",
    "beam_physical_size_kpc",
    "btfr_residual_dex",
    "build_anchor_sensitivity_rows",
    "curvature_radius_m",
    "default_ratio_scaling",
    "edge_radial_mean",
    "is_btfr_outlier",
    "pearson_correlation_coefficient",
    "propagate_gap_coverage_band",
    "rar_median_abs_residual_dex",
    "rar_predicted_observed_acceleration",
    "rar_residual_log10",
    "resolution_correction_factor",
    "supports_logarithmic_spiral",
    "sweep_coverage_ratio",
    "tension_acceleration_radial",
    "vectorised_resolution_correction",
]
