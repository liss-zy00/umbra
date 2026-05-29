from __future__ import annotations

"""
Markdown narrative generator for human-readable evidence summaries.

Journal-style prose for ``Scientific_Evidence_Narrative.md``; thresholds for
``notable result`` callouts use thresholds from ``src.config.thresholds``.

The narrative is intentionally **galaxy-agnostic**: every reliable galaxy gets the same
treatment. There are no per-object spotlight sections so individual systems are not
privileged in the report.
"""

from typing import Any

import numpy as np

from src.config.thresholds import (
    CALIBRATION_STATUS_SUBOPTIMAL,
    RELIABILITY_CUTS,
)


_FALLBACK_CALIBRATION_NOTE = (
    "Nota: Resultados baseados em calibração padrão por falha na leitura da mediana local."
)


def _is_notable_result(row: dict[str, Any]) -> bool:
    """Editorial highlight: strong edge coverage near unity and high Pearson on reliable points."""
    ratio = row.get("theoretical_gap_coverage_ratio")
    pearson = row.get("pearson_correlation")
    if ratio is None or pearson is None:
        return False
    if not np.isfinite(float(ratio)) or not np.isfinite(float(pearson)):
        return False
    r = float(ratio)
    p = float(pearson)
    return (
        RELIABILITY_CUTS.NARRATIVE_HIGHLIGHT_COVERAGE_MIN
        <= r
        <= RELIABILITY_CUTS.NARRATIVE_HIGHLIGHT_COVERAGE_MAX
        and p >= RELIABILITY_CUTS.NARRATIVE_HIGHLIGHT_PEARSON_MIN
    )


class NarrativeGenerator:
    """Builds Scientific Evidence Narrative Markdown from per-galaxy synthesis rows."""

    @staticmethod
    def generate_evidence_markdown(reliable_rows: list[dict[str, Any]]) -> str:
        """
        Compose a short, journal-style narrative: Abstract-style opening followed by an
        equally-weighted paragraph per reliable galaxy. No per-object spotlight section.
        """
        lines: list[str] = [
            "# Evidence for a radial magnetic pressure contribution to the outer-disk acceleration budget",
            "",
            "## Abstract",
            "",
            "We summarize rotation-curve residuals from the SPARC sample after strict baryonic "
            "closure, and compare them—at large galactocentric radii—to an ideal-MHD estimate of "
            "**magnetic pressure**–driven acceleration, \\(a_B \\approx -\\nabla P_B/\\rho_{\\mathrm{ion}}\\), "
            "constructed from Planck Stokes amplitudes and a two-parameter grid in effective gradient length "
            "and ionized gas fraction. Each entry below reports the **best-grid** scenario under internal "
            "reliability cuts (Galactic latitude and per-radius consistency flags). Correlation coefficients "
            "are Pearson \\(r\\) between the residual acceleration and \\(a_B\\) along radii retained as reliable.",
            "",
            "## Results",
            "",
        ]

        notable: list[str] = []
        for row in reliable_rows:
            name = str(row.get("galaxy_name", ""))
            morph = str(row.get("morphology_class", "Unclassified"))
            ratio = row.get("theoretical_gap_coverage_ratio")
            scale_kpc = row.get("scenario_gradient_scale_kpc")
            ion = row.get("scenario_ionization_fraction")
            scenario = str(row.get("best_scenario_name", ""))
            pearson = row.get("pearson_correlation")

            ratio_pct = float(ratio) * 100.0 if ratio is not None and np.isfinite(float(ratio)) else float("nan")
            scale_pc = (
                int(round(float(scale_kpc) * 1000.0))
                if scale_kpc is not None and np.isfinite(float(scale_kpc))
                else None
            )
            ion_pct = float(ion) * 100.0 if ion is not None and np.isfinite(float(ion)) else float("nan")

            ratio_s = f"{ratio_pct:.1f}%" if np.isfinite(ratio_pct) else "n/a"
            scale_s = f"{scale_pc} pc (nominal Δx = {float(scale_kpc):.3f} kpc)" if scale_pc is not None else "n/a"
            ion_s = f"{ion_pct:.4g}% (effective fraction {float(ion):.6f})" if ion is not None and np.isfinite(float(ion)) else "n/a"
            pearson_s = f"{float(pearson):.3f}" if pearson is not None and np.isfinite(float(pearson)) else "n/a"

            # Optional sanity diagnostics propagated from the engine; rendered only when present
            # and finite, so older summaries (or skipped scenarios) do not leak the literal "None".
            mu_g = row.get("magnetic_field_microgauss")
            mu_g_band = row.get("magnetic_field_within_plausibility_band")
            beta = row.get("plasma_beta_thermal_median")
            alfven = row.get("alfven_velocity_median_kms")
            rar_dex = row.get("rar_median_abs_residual_dex")
            phys_dx = row.get("physical_neighbor_separation_kpc")

            def _fmt(v: Any, spec: str) -> str:
                try:
                    fv = float(v)
                except (TypeError, ValueError):
                    return "n/a"
                if not np.isfinite(fv):
                    return "n/a"
                return format(fv, spec)

            para = (
                f"**{name}** ({morph}). Under scenario `{scenario}`, the outer-disk mean modulus of the "
                f"magnetic-pressure acceleration reaches **{ratio_s}** of the mean modulus of the baryonic "
                f"acceleration deficit at the same radii. Grid parameters: radial step proxy **{scale_s}**, "
                f"ionized mass fraction **{ion_s}**. The Pearson correlation between the residual and "
                f"\\(a_B\\) along reliable radii is **{pearson_s}**, indicating how closely the radial shapes align."
            )
            sanity_bits: list[str] = []
            mu_g_str = _fmt(mu_g, ".2f")
            if mu_g_str != "n/a":
                band_tag = "" if mu_g_band is None else (" (within μG plausibility band)" if mu_g_band else " (outside μG plausibility band)")
                sanity_bits.append(f"inferred field ≈ **{mu_g_str} μG**{band_tag}")
            beta_str = _fmt(beta, ".2e")
            if beta_str != "n/a":
                sanity_bits.append(f"median plasma β = **{beta_str}**")
            alfven_str = _fmt(alfven, ".2f")
            if alfven_str != "n/a":
                sanity_bits.append(f"v_A ≈ **{alfven_str} km/s**")
            rar_str = _fmt(rar_dex, ".2f")
            if rar_str != "n/a":
                sanity_bits.append(f"RAR offset ≈ **{rar_str} dex**")
            phys_dx_str = _fmt(phys_dx, ".2f")
            if phys_dx_str != "n/a":
                sanity_bits.append(f"physical Δx (sky→disk) ≈ **{phys_dx_str} kpc**")
            if sanity_bits:
                para += " Sanity diagnostics: " + "; ".join(sanity_bits) + "."
            lines.append(para)
            # If the Stokes calibration fell back to the emergency-only constant, the report must say so.
            if str(row.get("calibration_status", "")) == CALIBRATION_STATUS_SUBOPTIMAL:
                lines.append("")
                lines.append(f"> **{_FALLBACK_CALIBRATION_NOTE}**")
            lines.append("")

            if _is_notable_result(row):
                notable.append(name)

        if notable:
            lines.extend(
                [
                    "### Notable systems",
                    "",
                    "The following objects exhibit **particularly tight** agreement between outer-disk "
                    f"shape and amplitude of the magnetic-pressure proxy and the kinematic residual "
                    f"(coverage ratio near unity and Pearson \\(r \\ge {RELIABILITY_CUTS.NARRATIVE_HIGHLIGHT_PEARSON_MIN:.2f}\\)): "
                    f"**{', '.join(notable)}**. These cases merit follow-up with higher-resolution polarization "
                    "and three-dimensional MHD modeling beyond the exploratory parameterization used here.",
                    "",
                ]
            )

        lines.append("---")
        lines.append("")
        lines.append(
            "*Interpretation:* Evidence here supports **consistency checks** between a magnetic pressure "
            "proxy and SPARC residuals under stated assumptions; it does not, by itself, establish uniqueness "
            "against other non-baryonic or modified-gravity explanations. See `EG_MODEL_THEORY.md` for the "
            "conceptual hypothesis, the MHD bridge, and the inherent refinement parameters."
        )

        return "\n".join(lines)
