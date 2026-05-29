"""
Population-level consistency audit of pipeline outputs.

This module re-reads ``results/galaxy_analysis/*/*/*_summary.json`` after the run
finishes and writes ``results/step_validation/validation_consistency_audit.md`` with:

* sample-wide reliability counts and calibration provenance,
* ``theoretical_gap_coverage_ratio`` distribution by morphology,
* foreground-bias check (mean/median ratio for ``|b| < 30°`` vs ``|b| >= 30°``),
* coverage by inclination quartile and ``|galactic latitude|`` quartile,
* Pearson + Spearman correlations between coverage ratio and the major covariates
  (latitude, inclination, μG, plasma β, density-clamp/stellar fraction,
  RAR residual, physical Δx, distance), and
* basic μG-distribution stats on the reliable subset.

The audit is intentionally read-only: it never mutates per-galaxy outputs and only
writes the markdown report (plus an optional CSV merge for downstream notebooks).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from src.pipeline.constants import (
    GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
    GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME,
    RESULTS_GALAXY_ANALYSIS_DIRNAME,
    RESULTS_STEP_VALIDATION_DIRNAME,
    VALIDATION_CONSISTENCY_AUDIT_FILENAME,
)

# Columns that participate in the bias-correlation table. Ordered so the most
# physically meaningful covariates land at the top of the printed table.
_CANDIDATE_COVARIATES: tuple[str, ...] = (
    "galactic_latitude_degrees",
    "abs_galactic_latitude",
    "inclination_deg",
    "magnetic_field_microgauss",
    "plasma_beta_median",
    "density_clamped_fraction",
    "stellar_dominated_fraction",
    "rar_residual_dex",
    "physical_neighbor_separation_kpc",
    "distance_mpc",
)


def _safe_get(payload: dict[str, Any] | None, *keys: str, default: Any = None) -> Any:
    """Walk a nested mapping returning ``default`` whenever a key is missing."""
    current: Any = payload
    for key in keys:
        if not isinstance(current, dict):
            return default
        current = current.get(key)
        if current is None:
            return default
    return current


def _iter_galaxy_summary_paths(galaxy_analysis_root: Path) -> list[Path]:
    paths: list[Path] = []
    for tier in (
        GALAXY_ANALYSIS_TIER_SCIENTIFIC_EVIDENCE_DIRNAME,
        GALAXY_ANALYSIS_TIER_RELIABLE_DIRNAME,
        GALAXY_ANALYSIS_TIER_NOT_RELIABLE_DIRNAME,
    ):
        tier_dir = galaxy_analysis_root / tier
        if not tier_dir.is_dir():
            continue
        paths.extend(sorted(tier_dir.glob("*/*_summary.json")))
    return sorted(paths)


def count_galaxy_summary_json_paths(galaxy_analysis_root: Path) -> int:
    """Count ``*_summary.json`` files that the consistency audit will scan (tiers under ``galaxy_analysis``)."""
    return len(_iter_galaxy_summary_paths(galaxy_analysis_root))


def _load_summaries(galaxy_analysis_root: Path) -> pd.DataFrame:
    """Read every ``*_summary.json`` under both validation tiers and flatten audit columns."""
    rows: list[dict[str, Any]] = []
    for path in _iter_galaxy_summary_paths(galaxy_analysis_root):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        rows.append(
            {
                "galaxy_name": _safe_get(payload, "galaxy_identity", "galaxy_name"),
                "morphology_class": _safe_get(payload, "galaxy_identity", "morphology_class"),
                "hubble_t_type": _safe_get(payload, "galaxy_identity", "hubble_t_type"),
                "distance_mpc": _safe_get(payload, "galaxy_identity", "distance_local_mpc"),
                "galactic_latitude_degrees": _safe_get(
                    payload, "data_reliability_assessment", "galactic_latitude_degrees"
                ),
                "validation_is_reliable": bool(
                    _safe_get(
                        payload,
                        "data_reliability_assessment",
                        "validation_is_reliable",
                        default=_safe_get(
                            payload,
                            "data_reliability_assessment",
                            "is_reliable",
                            default=False,
                        ),
                    )
                ),
                "audit_flags": list(
                    _safe_get(payload, "data_reliability_assessment", "audit_flags", default=[]) or []
                ),
                "reliability_failure_codes": list(
                    _safe_get(
                        payload, "data_reliability_assessment", "reliability_failure_codes", default=[]
                    )
                    or []
                ),
                "calibration_status": _safe_get(payload, "calibration_metadata", "calibration_status"),
                "ratio": _safe_get(
                    payload, "best_theoretical_scenario", "theoretical_gap_coverage_ratio"
                ),
                "pearson": _safe_get(
                    payload, "best_theoretical_scenario", "pearson_correlation"
                ),
                "magnetic_field_microgauss": _safe_get(
                    payload, "plasma_diagnostics", "magnetic_field_microgauss"
                ),
                "magnetic_field_within_band": _safe_get(
                    payload, "plasma_diagnostics", "magnetic_field_within_plausibility_band"
                ),
                "plasma_beta_median": _safe_get(
                    payload, "plasma_diagnostics", "plasma_beta_thermal_median"
                ),
                "rar_residual_dex": _safe_get(
                    payload, "rar_diagnostics", "median_absolute_residual_dex"
                ),
                "btfr_residual_dex": _safe_get(payload, "btfr_diagnostics", "residual_dex"),
                "btfr_is_outlier": bool(
                    _safe_get(payload, "btfr_diagnostics", "is_outlier", default=False) or False
                ),
                "ratio_p16": _safe_get(
                    payload, "uncertainty_propagation", "theoretical_gap_coverage_ratio_p16"
                ),
                "ratio_p84": _safe_get(
                    payload, "uncertainty_propagation", "theoretical_gap_coverage_ratio_p84"
                ),
                "density_clamped_fraction": _safe_get(
                    payload, "density_audit", "density_clamped_fraction"
                ),
                "stellar_dominated_fraction": _safe_get(
                    payload, "density_audit", "stellar_dominated_fraction"
                ),
                "inclination_deg": _safe_get(
                    payload, "kinematic_metadata", "inclination_deg"
                ),
                "physical_neighbor_separation_kpc": _safe_get(
                    payload, "gradient_scale_audit", "physical_neighbor_separation_kpc"
                ),
            }
        )
    df = pd.DataFrame(rows)
    if "galactic_latitude_degrees" in df.columns:
        df["abs_galactic_latitude"] = df["galactic_latitude_degrees"].abs()
    else:
        df["abs_galactic_latitude"] = np.nan
    return df


def _format_distribution_table(grouped: pd.DataFrame) -> str:
    """Render a small grouped describe() into a simple markdown table."""
    lines: list[str] = []
    lines.append("| group | count | mean | median | std | min | max |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for label, row in grouped.iterrows():
        try:
            count = int(row.get("count", 0))
        except (TypeError, ValueError):
            count = 0
        mean = float(row.get("mean", float("nan")))
        median = float(row.get("median", float("nan")))
        std = float(row.get("std", float("nan")))
        mn = float(row.get("min", float("nan")))
        mx = float(row.get("max", float("nan")))
        lines.append(
            f"| {label} | {count} | {mean:.3f} | {median:.3f} | {std:.3f} | {mn:.3f} | {mx:.3f} |"
        )
    return "\n".join(lines)


def _format_quartile_table(label: str, df: pd.DataFrame, ratio_col: str = "ratio") -> str:
    """Render a quartile breakdown of ``ratio_col`` using ``label`` as the column name."""
    sub = df.dropna(subset=[label, ratio_col])
    if sub.shape[0] < 4:
        return f"_(insufficient {label} data — N={sub.shape[0]})_"
    quartiles = pd.qcut(sub[label], q=4, duplicates="drop")
    grouped = (
        sub.assign(_q=quartiles)
        .groupby("_q", observed=False)[ratio_col]
        .agg(["count", "mean", "median"])
    )
    lines: list[str] = []
    lines.append(f"| {label} quartile | count | mean ratio | median ratio |")
    lines.append("|---|---:|---:|---:|")
    for q, row in grouped.iterrows():
        lines.append(
            f"| {q} | {int(row['count'])} | {row['mean']:.3f} | {row['median']:.3f} |"
        )
    return "\n".join(lines)


def _format_correlation_table(df: pd.DataFrame, ratio_col: str = "ratio") -> str:
    """Pearson + Spearman correlation of the ratio against each candidate covariate."""
    import warnings
    from scipy.stats import ConstantInputWarning  # local import — avoids hard dep on import
    lines: list[str] = ["| column | Pearson | Spearman | N |", "|---|---:|---:|---:|"]
    for col in _CANDIDATE_COVARIATES:
        if col not in df.columns:
            continue
        sub = df[[ratio_col, col]].dropna()
        if sub.shape[0] < 5:
            lines.append(f"| {col} | – | – | {sub.shape[0]} |")
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ConstantInputWarning)
                warnings.simplefilter("ignore", RuntimeWarning)
                r_p = float(sub[ratio_col].corr(sub[col], method="pearson"))
                r_s = float(sub[ratio_col].corr(sub[col], method="spearman"))
        except Exception:  # pragma: no cover — defensive against pathological data
            r_p, r_s = float("nan"), float("nan")
        lines.append(f"| {col} | {r_p:.3f} | {r_s:.3f} | {sub.shape[0]} |")
    return "\n".join(lines)


def build_consistency_audit_markdown(
    df: pd.DataFrame,
    *,
    supreme_rows: list[dict[str, Any]] | None = None,
) -> str:
    """Produce the full markdown text for ``validation_consistency_audit.md``."""
    n_total = int(df.shape[0])
    if n_total > 0 and "validation_is_reliable" in df.columns:
        reliable = df[df["validation_is_reliable"].astype(bool)].copy()
    elif n_total > 0 and "is_reliable" in df.columns:
        reliable = df[df["is_reliable"].astype(bool)].copy()
    else:
        reliable = df.iloc[0:0].copy()
    n_reliable = int(reliable.shape[0])

    calibration_counts = (
        df["calibration_status"].value_counts(dropna=False).to_dict() if "calibration_status" in df.columns else {}
    )
    band_counts = (
        df["magnetic_field_within_band"].value_counts(dropna=False).to_dict()
        if "magnetic_field_within_band" in df.columns
        else {}
    )

    parts: list[str] = []
    parts.append("# Consistency Audit Report")
    parts.append("")
    parts.append(
        "Read-only audit produced after every pipeline run. It re-loads each "
        "`results/galaxy_analysis/{validation_scientific_evidence|validation_reliable|validation_not_reliable}/<galaxy>/<galaxy>_summary.json` and looks for "
        "background biases in the gap-coverage signal: foreground (galactic latitude), "
        "morphology selection, inclination deprojection, calibration provenance, and "
        "the magnetic-pressure self-fulfillment risk."
    )
    parts.append("")
    parts.append("## Sample overview")
    parts.append("")
    parts.append(f"- Total summaries loaded: **{n_total}**")
    parts.append(f"- Reliable galaxies: **{n_reliable}** ({(n_reliable / n_total * 100.0) if n_total else 0.0:.1f}%)")
    if calibration_counts:
        parts.append(
            "- Calibration provenance: "
            + ", ".join(f"`{k}` = {v}" for k, v in calibration_counts.items())
        )
    if band_counts:
        parts.append(
            "- μG plausibility band: "
            + ", ".join(f"`{k}` = {v}" for k, v in band_counts.items())
        )

    # Failure-mode breakdown — a frequency table of every audit code emitted across the run.
    # This is the "if it fails, tell me WHY it failed" answer at the population level.
    parts.append("")
    parts.append("## Failure-mode breakdown (audit_flags frequency)")
    parts.append("")
    parts.append("Every processed galaxy contributes its full ``audit_flags`` list. A galaxy can")
    parts.append("appear under multiple codes when several gates fail simultaneously.")
    parts.append("")
    if "audit_flags" in df.columns and n_total > 0:
        flag_counts: dict[str, int] = {}
        for flags in df["audit_flags"]:
            if not flags:
                continue
            for code in flags:
                flag_counts[str(code)] = flag_counts.get(str(code), 0) + 1
        if flag_counts:
            parts.append("| audit code | galaxies tagged | share of sample |")
            parts.append("|---|---:|---:|")
            for code in sorted(flag_counts.keys(), key=lambda c: (-flag_counts[c], c)):
                count = int(flag_counts[code])
                share = count / n_total * 100.0 if n_total else 0.0
                parts.append(f"| `{code}` | {count} | {share:.1f}% |")
        else:
            parts.append("_No audit codes recorded — every processed galaxy is reliable._")
    else:
        parts.append("_audit_flags column not present in summaries (older pipeline run?)._")

    parts.append("")
    if reliable.empty:
        parts.append("> No reliable galaxies in this run — bias audit skipped.")
        return "\n".join(parts) + "\n"

    parts.append("## 1. Coverage ratio by morphology (reliable subset)")
    parts.append("")
    grp = (
        reliable.groupby("morphology_class")["ratio"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .sort_values("count", ascending=False)
    )
    parts.append(_format_distribution_table(grp))

    parts.append("")
    parts.append("## 2. Foreground bias check — galactic latitude split")
    parts.append("")
    low = reliable[reliable["abs_galactic_latitude"] < 30.0]["ratio"].dropna()
    high = reliable[reliable["abs_galactic_latitude"] >= 30.0]["ratio"].dropna()
    parts.append("| `|b|` group | N | mean ratio | median ratio |")
    parts.append("|---|---:|---:|---:|")
    parts.append(
        f"| < 30° (high foreground risk) | {int(low.size)} | "
        f"{(low.mean() if low.size else float('nan')):.3f} | "
        f"{(low.median() if low.size else float('nan')):.3f} |"
    )
    parts.append(
        f"| ≥ 30° (clean) | {int(high.size)} | "
        f"{(high.mean() if high.size else float('nan')):.3f} | "
        f"{(high.median() if high.size else float('nan')):.3f} |"
    )
    if low.size and high.size and high.mean() != 0:
        parts.append("")
        parts.append(
            f"_Foreground inflation factor (mean ratio low/high)_: **{low.mean() / high.mean():.2f}× —** values >> 1 indicate the "
            "magnetic-pressure signal is contaminated by Milky Way dust polarization in low-latitude windows."
        )

    parts.append("")
    parts.append("## 3. Coverage by `|galactic latitude|` quartile (reliable subset)")
    parts.append("")
    parts.append(_format_quartile_table("abs_galactic_latitude", reliable))

    parts.append("")
    parts.append("## 4. Coverage by inclination quartile (reliable subset)")
    parts.append("")
    parts.append(_format_quartile_table("inclination_deg", reliable))

    parts.append("")
    parts.append("## 5. Correlations vs theoretical_gap_coverage_ratio (reliable subset)")
    parts.append("")
    parts.append(_format_correlation_table(reliable))
    parts.append("")
    parts.append(
        "_Interpretation cues._ Strong positive Spearman with `magnetic_field_microgauss` "
        "is structurally trivial (a_B grows with B), but a strong correlation with "
        "`abs_galactic_latitude` would point to remaining foreground contamination, and "
        "a strong correlation with `inclination_deg` to deprojection systematics."
    )

    parts.append("")
    parts.append("## 6. Correlations vs theoretical_gap_coverage_ratio (entire sample)")
    parts.append("")
    parts.append(_format_correlation_table(df))

    # Phase-1 refinement: per-population RAR and BTFR summary.
    parts.append("")
    parts.append("## 7a. RAR consistency on the reliable sample")
    parts.append("")
    rar_arr = reliable["rar_residual_dex"].dropna().to_numpy(dtype=float) if "rar_residual_dex" in reliable.columns else np.array([], dtype=float)
    if rar_arr.size > 0:
        med = float(np.median(np.abs(rar_arr)))
        mad = float(np.median(np.abs(np.abs(rar_arr) - med))) * 1.4826
        parts.append(
            f"- N = {int(rar_arr.size)}  median |residual| = {med:.3f} dex  "
            f"MAD-robust σ = {mad:.3f} dex"
        )
    else:
        parts.append("- (no reliable galaxy carried a finite RAR residual)")

    parts.append("")
    parts.append("## 7b. BTFR distance-to-ridge on the reliable sample")
    parts.append("")
    btfr_arr = (
        reliable["btfr_residual_dex"].dropna().to_numpy(dtype=float)
        if "btfr_residual_dex" in reliable.columns
        else np.array([], dtype=float)
    )
    if btfr_arr.size > 0:
        outliers = int(np.count_nonzero(np.abs(btfr_arr) > 0.5))
        outlier_pct = outliers / float(btfr_arr.size) * 100.0
        parts.append(
            f"- N = {int(btfr_arr.size)}  median |residual| = {float(np.median(np.abs(btfr_arr))):.3f} dex"
        )
        parts.append(
            f"- |residual| > 0.5 dex: **{outliers} galaxies ({outlier_pct:.1f}%)**"
        )
    else:
        parts.append("- (no reliable galaxy carried a finite BTFR residual — HI mass / 3.6 µm leg missing)")

    parts.append("")
    parts.append("## 8. Magnetic-field distribution on the reliable sample")
    parts.append("")
    b_arr = reliable["magnetic_field_microgauss"].dropna().to_numpy(dtype=float)
    if b_arr.size > 0:
        parts.append(
            f"- N = {int(b_arr.size)}  median = {float(np.median(b_arr)):.3f} μG  "
            f"mean = {float(np.mean(b_arr)):.3f} μG  std = {float(np.std(b_arr)):.3f} μG"
        )
        parts.append(
            f"- p5 = {float(np.quantile(b_arr, 0.05)):.3f} μG  "
            f"p95 = {float(np.quantile(b_arr, 0.95)):.3f} μG"
        )
        parts.append("")
        parts.append(
            "See `graphs/magnetic_field_distribution.png` and `graphs/magnetic_field_distribution.md` for the "
            "histogram + rug, and `graphs/magnetic_field_vs_galaxy_scale.png` / `.md` for the B-field vs outer-radius view."
        )

    if supreme_rows:
        parts.append("")
        parts.append("## 9. Supreme-row consistency cross-check")
        parts.append("")
        parts.append(
            f"- Supreme rows (passed inclination + reliability): **{len(supreme_rows)}**"
        )
        parts.append(
            f"- Reliable summaries on disk: **{n_reliable}**  "
            f"(mismatch ⇒ stale or interrupted run)"
        )

    parts.append("")
    parts.append(
        "> *This report is regenerated automatically every run. To inspect raw inputs "
        "see the per-galaxy `*_summary.json` files under `results/galaxy_analysis/`.*"
    )
    return "\n".join(parts) + "\n"


def write_consistency_audit_report(
    *,
    processed_root: Path,
    supreme_rows: list[dict[str, Any]] | None = None,
    master_inventory_rows: list[dict[str, Any]] | None = None,
) -> Path:
    """Compute the audit and persist ``results/step_validation/validation_consistency_audit.md``."""
    galaxy_analysis_root = processed_root / RESULTS_GALAXY_ANALYSIS_DIRNAME
    df = _load_summaries(galaxy_analysis_root)
    text = build_consistency_audit_markdown(df, supreme_rows=supreme_rows)
    out_dir = processed_root / RESULTS_STEP_VALIDATION_DIRNAME
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / VALIDATION_CONSISTENCY_AUDIT_FILENAME
    out.write_text(text, encoding="utf-8")
    return out
