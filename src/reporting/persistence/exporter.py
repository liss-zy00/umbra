from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from src.models.galaxy import Galaxy
from src.pipeline.constants import (
    GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
    RESULTS_STEP_DISCOVERY_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_ANCHOR_SENSITIVITY_SWEEP_FILENAME,
    VALIDATION_CATALOG_ALL_FILENAME,
    VALIDATION_NOT_RELIABLE_FILENAME,
    VALIDATION_RELIABLE_FILENAME,
    VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME,
)


@dataclass(frozen=True, slots=True)
class DiscoveryEliteBandInfo:
    """
    Summary of the dynamic ρ-req band computed from the Discovery census.

    The cohort (``discovery_cohort_rho_band.csv``) lists reliable galaxies whose
    ``required_density_for_match`` lies within a multiplicative window around the
    *population* median of the reliable subset, so internal consistency is judged relative
    to the catalogue scale rather than an absolute WIM anchor.
    """

    median_required_density_kg_m3: float
    band_min_kg_m3: float
    band_max_kg_m3: float
    n_reliable: int
    n_elite: int
    relative_min: float
    relative_max: float

    @property
    def is_valid(self) -> bool:
        return bool(
            np.isfinite(self.median_required_density_kg_m3)
            and self.median_required_density_kg_m3 > 0.0
        )


def _inventory_row_for_csv(row: dict[str, Any]) -> dict[str, Any]:
    """Flatten list-valued audit fields to comma-separated strings for CSV export."""
    out = dict(row)
    af = out.get("audit_flags")
    if isinstance(af, (list, tuple)):
        out["audit_flags"] = ",".join(str(x) for x in af)
    elif af is None:
        out["audit_flags"] = ""
    else:
        out["audit_flags"] = str(af)
    rfc = out.get("reliability_failure_codes")
    if isinstance(rfc, (list, tuple)):
        out["reliability_failure_codes"] = ",".join(str(x) for x in rfc)
    elif rfc is None:
        out["reliability_failure_codes"] = ""
    else:
        out["reliability_failure_codes"] = str(rfc)
    note = out.get("note")
    if note is None:
        note = out.get("reliability_note", "")
    out["note"] = str(note) if note is not None else ""
    return out


def _sanitize_for_json(obj: Any) -> Any:
    """Replace non-finite floats with ``None`` so ``json.dumps(..., allow_nan=False)`` is valid JSON."""
    if isinstance(obj, float):
        return obj if np.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_for_json(v) for v in obj]
    return obj


class ResultPersister:
    """CSV / JSON persistence under ``results/`` (step_discovery, step_validation, galaxy_analysis, grid_results)."""

    def __init__(self, processed_root: Path) -> None:
        self.processed_root = processed_root
        self.discovery_root = self.processed_root / RESULTS_STEP_DISCOVERY_DIRNAME
        self.validation_root = self.processed_root / RESULTS_STEP_VALIDATION_DIRNAME
        self.individual_root = self.processed_root / RESULTS_GALAXY_ANALYSIS_DIRNAME
        self.grid_root = self.processed_root / "grid_results"
        self.processed_root.mkdir(parents=True, exist_ok=True)

    def save_galaxy(self, galaxy: Galaxy, summary: dict[str, Any]) -> tuple[Path, Path, Path]:
        """
        Persist the **best-scenario** per-shell table plus the structured galaxy summary.

        Output directory is ``galaxy_analysis/{validation_reliable|validation_not_reliable}/<galaxy.name>/``,
        chosen from ``summary["data_reliability_assessment"]["validation_is_reliable"]``.

        The CSV and JSON are **not** redundant:

        * ``*_analysis.csv`` — one row per radial shell after the winning grid scenario; wide numeric
          kinematic / field / density columns (``Galaxy.rotation_dataframe``).
        * ``*_summary.json`` — nested audit and provenance payload (reliability, calibration, NED/Planck
          context, MC bands on the headline ratio, optional ``neighbour_ring`` geometry). This is what
          :mod:`src.pipeline.aggregation.audit` re-reads to build ``validation_consistency_audit.md``;
          those fields are not all present or not all recoverable from the shell CSV alone.

        Removing the JSON would break the population audit unless the inventory schema is extended
        and the audit loader is rewritten to consume that alternative.
        """
        self.individual_root.mkdir(parents=True, exist_ok=True)
        dra = summary.get("data_reliability_assessment")
        reliable = bool((dra or {}).get("validation_is_reliable", False))
        tier = GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME if reliable else GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME
        galaxy_dir = self.individual_root / tier / galaxy.name
        galaxy_dir.mkdir(parents=True, exist_ok=True)
        csv_path = galaxy_dir / f"{galaxy.name}_analysis.csv"
        json_path = galaxy_dir / f"{galaxy.name}_summary.json"
        galaxy.rotation_dataframe.to_csv(csv_path, index=False)
        json_path.write_text(
            json.dumps(
                _sanitize_for_json(summary),
                indent=2,
                ensure_ascii=True,
                allow_nan=False,
            ),
            encoding="utf-8",
        )
        return galaxy_dir, csv_path, json_path

    def save_full_catalog_inventory(self, rows: list[dict[str, Any]]) -> Path:
        import pandas as pd

        self.validation_root.mkdir(parents=True, exist_ok=True)
        output_path = self.validation_root / VALIDATION_CATALOG_ALL_FILENAME
        columns = [
            "galaxy_name",
            "right_ascension_degrees",
            "declination_degrees",
            "galactic_latitude_degrees",
            "validation_is_reliable",
            "ned_status",
            "mean_acceleration_gap_edge",
            "morphology_class",
            "hubble_type",
            "audit_flags",
            "reliability_failure_codes",
            "note",
        ]
        if not rows:
            empty = pd.DataFrame(columns=columns)
            empty.to_csv(output_path, index=False)
            return output_path
        normalized = [_inventory_row_for_csv(r) for r in rows]
        frame = pd.DataFrame(normalized).reindex(columns=columns)
        frame.to_csv(output_path, index=False)
        return output_path

    def save_validation_reliability_partitions(self, rows: list[dict[str, Any]]) -> tuple[Path, Path]:
        """
        Split the validation inventory on ``validation_is_reliable`` (same columns as
        :meth:`save_full_catalog_inventory`) for operators who prefer two CSVs over filtering.
        """
        import pandas as pd

        self.validation_root.mkdir(parents=True, exist_ok=True)
        reliable_path = self.validation_root / VALIDATION_RELIABLE_FILENAME
        not_reliable_path = self.validation_root / VALIDATION_NOT_RELIABLE_FILENAME
        columns = [
            "galaxy_name",
            "right_ascension_degrees",
            "declination_degrees",
            "galactic_latitude_degrees",
            "validation_is_reliable",
            "ned_status",
            "mean_acceleration_gap_edge",
            "morphology_class",
            "hubble_type",
            "audit_flags",
            "reliability_failure_codes",
            "note",
        ]
        if not rows:
            empty = pd.DataFrame(columns=columns)
            empty.to_csv(reliable_path, index=False)
            empty.to_csv(not_reliable_path, index=False)
            return reliable_path, not_reliable_path
        normalized = [_inventory_row_for_csv(r) for r in rows]
        frame = pd.DataFrame(normalized).reindex(columns=columns)
        ok = frame["validation_is_reliable"].map(lambda x: bool(x))
        frame.loc[ok].to_csv(reliable_path, index=False)
        frame.loc[~ok].to_csv(not_reliable_path, index=False)
        return reliable_path, not_reliable_path

    def save_scientific_evidence_report(self, rows: list[dict[str, Any]]) -> Path:
        """
        Wide CSV: one row per **direct-evidence** galaxy.

        Rows are produced by :func:`src.pipeline.aggregation.global_aggregator.build_direct_evidence_rows`
        from *supreme* inventory rows (``validation_is_reliable`` plus the inclination-quality gate),
        then filtered by morphology (Im or photometric LSB) and by the coverage-ratio window
        ``DIRECT_EVIDENCE_RATIO_MIN/MAX``. Each row merges the best grid scenario with Discovery
        headline columns and uncertainty bands.
        """
        import pandas as pd

        self.validation_root.mkdir(parents=True, exist_ok=True)
        output_path = self.validation_root / VALIDATION_SCIENTIFIC_EVIDENCE_REPORT_FILENAME
        columns = [
            "galaxy_name",
            "right_ascension_degrees",
            "declination_degrees",
            "best_scenario_name",
            "scenario_gradient_scale_kpc",
            "scenario_ionization_fraction",
            "theoretical_gap_coverage_ratio",
            "theoretical_gap_coverage_ratio_p16",
            "theoretical_gap_coverage_ratio_p84",
            "theoretical_gap_coverage_ratio_full_lorentz",
            "accel_pressure_only_m_s2",
            "accel_full_lorentz_m_s2",
            "pearson_correlation",
            "required_density_for_match",
            "convergence_ratio_to_wim",
            "convergence_ratio_to_wim_p16",
            "convergence_ratio_to_wim_p84",
            "calibrated_convergence_ratio",
            "rar_median_abs_residual_dex",
            "btfr_residual_dex",
            "warn_btfr_outlier",
            "warn_no_pa",
            "warn_pure_pressure",
            "prov_high_precision_dist",
            "position_angle_deg",
            "wim_density_match_score",
            "validation_status",
            "delta_r_m",
            "physical_pixel_scale_kpc",
            # Beam-smearing diagnostic: ``k = D_beam / L_grad`` quantifies how much the
            # finite Planck beam (5'  FWHM, GNILC R3.00) dilutes the reconstructed gradient
            # at the grid-search scale ``scenario_gradient_scale_kpc``. ``k >> 1`` ⇒ the
            # reconstructed ``ρ_req`` is over-estimated by the same factor. The "_grid"
            # column matches the engine's grid knob (the ~50–500 pc scale assumed by the
            # MHD ansatz); the "_discovery" column reports k against the data-derived
            # physical pixel separation. See :mod:`src.physics.diagnostics.resolution`.
            "distance_mpc",
            "beam_size_kpc",
            "resolution_correction_factor",
            "resolution_correction_factor_discovery",
            "resolution_correction_factor_grid",
            "discovery_is_reliable",
            "discovery_is_elite",
            "discovery_elite_band_source",
            "morphology_class",
            "ratio_distance_to_unity",
            "r_max_kpc",
        ]
        if not rows:
            pd.DataFrame(columns=columns).to_csv(output_path, index=False)
            return output_path
        pd.DataFrame(rows).reindex(columns=columns).to_csv(output_path, index=False)
        return output_path

    def save_anchor_sensitivity_sweep(self, rows: list[dict[str, Any]]) -> Path:
        """
        Persist the CMB-vs-turbulent-equipartition anchor-sensitivity sweep.

        One row per galaxy, columns ``galaxy_name``, ``ratio_cmb``,
        ``ratio_turbulent_equipartition``, ``ratio_log10_delta`` (the log10 ratio between
        the two anchors). Lives at ``results/step_validation/validation_anchor_sensitivity_sweep.csv``
        alongside ``validation_scientific_evidence_report.csv`` — see
        :mod:`src.physics.diagnostics.anchor_sensitivity`.
        """
        import pandas as pd

        self.validation_root.mkdir(parents=True, exist_ok=True)
        output_path = self.validation_root / VALIDATION_ANCHOR_SENSITIVITY_SWEEP_FILENAME
        columns = [
            "galaxy_name",
            "ratio_cmb",
            "ratio_turbulent_equipartition",
            "ratio_log10_delta",
        ]
        if not rows:
            pd.DataFrame(columns=columns).to_csv(output_path, index=False)
            return output_path
        pd.DataFrame(rows).reindex(columns=columns).to_csv(output_path, index=False)
        return output_path

    def save_thesis_summary(self, rows: list[dict[str, Any]]) -> Path:
        import pandas as pd

        out = self.processed_root / "EG_Final_Thesis_Summary.csv"
        pd.DataFrame(rows).to_csv(out, index=False)
        return out

    def save_grid_scenario_report(self, scenario_name: str, rows: list[dict[str, Any]]) -> Path:
        import pandas as pd

        self.grid_root.mkdir(parents=True, exist_ok=True)
        output_path = self.grid_root / f"{scenario_name}.csv"
        pd.DataFrame(rows).to_csv(output_path, index=False)
        return output_path

    # Dynamic ρ-req band: galaxies are flagged ``discovery_is_elite`` when their per-galaxy
    # ``required_density_for_match`` lies within ``[REL_MIN · median, REL_MAX · median]``,
    # where ``median`` is computed over the **reliable** subset only. The 0.7–1.3 window
    # selects galaxies that are internally consistent with the population — i.e. those
    # whose density agrees with the population scale within ±30 % — and therefore probes
    # the model's coherence independently of where the absolute density actually lands.
    DISCOVERY_ELITE_BAND_RELATIVE_MIN: float = 0.7
    DISCOVERY_ELITE_BAND_RELATIVE_MAX: float = 1.3

    def save_discovery_evidence(
        self, rows: list[dict[str, Any]]
    ) -> tuple[Path, Path, Path, Path, DiscoveryEliteBandInfo]:
        """
        Persist **Phase 1 (Discovery census)** artefacts for the unified pipeline.

        Writes four CSVs under ``results/step_discovery/``:

        - ``discovery_all.csv`` — every processed galaxy (coverage inventory).
        - ``discovery_reliable.csv`` — only rows with ``discovery_is_reliable == True``,
          i.e. the subset whose median ``required_density_for_match`` is finite and not
          flagged by the latitude / Planck-skip audits.
        - ``discovery_not_reliable.csv`` — complement of ``discovery_reliable`` (same schema;
          ``audit_flags`` lists the blocking reason).
        - ``discovery_cohort_rho_band.csv`` — the **internally consistent** subset:
          reliable rows whose ``required_density_for_match`` lies inside
          ``[REL_MIN · median, REL_MAX · median]``, where ``median`` is the population
          median of ``required_density_for_match`` over the reliable subset and
          ``REL_MIN`` / ``REL_MAX`` are :attr:`DISCOVERY_ELITE_BAND_RELATIVE_MIN` and
          :attr:`DISCOVERY_ELITE_BAND_RELATIVE_MAX`.

        Returns the four persisted paths plus a :class:`DiscoveryEliteBandInfo` block
        describing the dynamic band actually applied (median + numerical limits).
        """
        import pandas as pd

        self.processed_root.mkdir(parents=True, exist_ok=True)
        self.discovery_root.mkdir(parents=True, exist_ok=True)
        frame = pd.DataFrame(rows)
        all_path = self.discovery_root / "discovery_all.csv"
        passed_path = self.discovery_root / "discovery_reliable.csv"
        rejected_path = self.discovery_root / "discovery_not_reliable.csv"
        cohort_path = self.discovery_root / "discovery_cohort_rho_band.csv"
        frame.to_csv(all_path, index=False)

        if "discovery_is_reliable" in frame.columns:
            reliable_mask = frame["discovery_is_reliable"].astype(bool)
        else:
            reliable_mask = pd.Series(False, index=frame.index)
        passed_frame = frame.loc[reliable_mask]
        rejected_frame = frame.loc[~reliable_mask]
        passed_frame.to_csv(passed_path, index=False)
        rejected_frame.to_csv(rejected_path, index=False)

        # Dynamic band: median of ``required_density_for_match`` over the passed rows,
        # then ``[REL_MIN · median, REL_MAX · median]``. NaN / non-positive values are
        # dropped from the median so a single bad row never contaminates the band.
        median_value = float("nan")
        if "required_density_for_match" in passed_frame.columns:
            rho_series = pd.to_numeric(
                passed_frame["required_density_for_match"], errors="coerce"
            )
            rho_finite = rho_series[np.isfinite(rho_series) & (rho_series > 0.0)]
            if not rho_finite.empty:
                median_value = float(rho_finite.median())

        if np.isfinite(median_value) and median_value > 0.0:
            band_min = float(self.DISCOVERY_ELITE_BAND_RELATIVE_MIN) * median_value
            band_max = float(self.DISCOVERY_ELITE_BAND_RELATIVE_MAX) * median_value
            rho_for_mask = pd.to_numeric(
                passed_frame["required_density_for_match"], errors="coerce"
            )
            elite_mask = (
                rho_for_mask.notna()
                & (rho_for_mask > 0.0)
                & (rho_for_mask >= band_min)
                & (rho_for_mask <= band_max)
            )
        else:
            band_min = float("nan")
            band_max = float("nan")
            elite_mask = pd.Series(False, index=passed_frame.index)

        elite_frame = passed_frame.loc[elite_mask]
        elite_frame.to_csv(cohort_path, index=False)

        band_info = DiscoveryEliteBandInfo(
            median_required_density_kg_m3=median_value,
            band_min_kg_m3=band_min,
            band_max_kg_m3=band_max,
            n_reliable=int(passed_frame.shape[0]),
            n_elite=int(elite_frame.shape[0]),
            relative_min=float(self.DISCOVERY_ELITE_BAND_RELATIVE_MIN),
            relative_max=float(self.DISCOVERY_ELITE_BAND_RELATIVE_MAX),
        )
        return all_path, passed_path, rejected_path, cohort_path, band_info
