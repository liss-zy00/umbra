# Consistency Audit Report

Read-only audit produced after every pipeline run. It re-loads each `results/galaxy_analysis/{validation_scientific_evidence|validation_reliable|validation_not_reliable}/<galaxy>/<galaxy>_summary.json` and looks for background biases in the gap-coverage signal: foreground (galactic latitude), morphology selection, inclination deprojection, calibration provenance, and the magnetic-pressure self-fulfillment risk.

## Sample overview

- Total summaries loaded: **158**
- Reliable galaxies: **28** (17.7%)
- Calibration provenance: `OPTIMAL` = 158
- μG plausibility band: `True` = 157, `False` = 1

## Failure-mode breakdown (audit_flags frequency)

Every processed galaxy contributes its full ``audit_flags`` list. A galaxy can
appear under multiple codes when several gates fail simultaneously.

| audit code | galaxies tagged | share of sample |
|---|---:|---:|
| `warn_no_pa` | 130 | 82.3% |
| `fail_reliable_fraction` | 128 | 81.0% |
| `fail_stellar` | 79 | 50.0% |
| `warn_pure_pressure` | 40 | 25.3% |
| `warn_btfr_outlier` | 27 | 17.1% |
| `fail_inclination` | 12 | 7.6% |
| `fail_density_clamped` | 1 | 0.6% |
| `warn_b_band` | 1 | 0.6% |

## 1. Coverage ratio by morphology (reliable subset)

| group | count | mean | median | std | min | max |
|---|---:|---:|---:|---:|---:|---:|
| Dwarf Irregulars (Im) | 17 | 105.454 | 0.257 | 263.413 | 0.000 | 975.278 |
| Sm/Sd Spirals | 8 | 12.921 | 0.442 | 27.876 | 0.001 | 79.459 |
| Blue Compact Dwarfs (BCD) | 3 | 0.330 | 0.176 | 0.405 | 0.025 | 0.790 |

## 2. Foreground bias check — galactic latitude split

| `|b|` group | N | mean ratio | median ratio |
|---|---:|---:|---:|
| < 30° (high foreground risk) | 4 | 3.587 | 0.068 |
| ≥ 30° (clean) | 24 | 78.447 | 0.261 |

_Foreground inflation factor (mean ratio low/high)_: **0.05× —** values >> 1 indicate the magnetic-pressure signal is contaminated by Milky Way dust polarization in low-latitude windows.

## 3. Coverage by `|galactic latitude|` quartile (reliable subset)

| abs_galactic_latitude quartile | count | mean ratio | median ratio |
|---|---:|---:|---:|
| (20.895999999999997, 36.492] | 7 | 141.679 | 0.094 |
| (36.492, 48.461] | 7 | 26.496 | 0.018 |
| (48.461, 70.082] | 7 | 0.518 | 0.265 |
| (70.082, 89.41] | 7 | 102.317 | 21.406 |

## 4. Coverage by inclination quartile (reliable subset)

| inclination_deg quartile | count | mean ratio | median ratio |
|---|---:|---:|---:|
| (29.999, 53.75] | 7 | 139.597 | 0.025 |
| (53.75, 60.5] | 7 | 2.558 | 0.726 |
| (60.5, 67.25] | 7 | 11.388 | 0.153 |
| (67.25, 90.0] | 7 | 117.467 | 0.790 |

## 5. Correlations vs theoretical_gap_coverage_ratio (reliable subset)

| column | Pearson | Spearman | N |
|---|---:|---:|---:|
| galactic_latitude_degrees | -0.189 | -0.317 | 28 |
| abs_galactic_latitude | -0.039 | 0.091 | 28 |
| inclination_deg | -0.132 | 0.341 | 28 |
| magnetic_field_microgauss | -0.137 | 0.435 | 28 |
| plasma_beta_median | -0.073 | -0.396 | 28 |
| density_clamped_fraction | 0.524 | 0.747 | 28 |
| stellar_dominated_fraction | -0.254 | -0.421 | 28 |
| rar_residual_dex | -0.110 | 0.038 | 28 |
| physical_neighbor_separation_kpc | -0.162 | 0.009 | 28 |
| distance_mpc | -0.158 | 0.032 | 28 |

_Interpretation cues._ Strong positive Spearman with `magnetic_field_microgauss` is structurally trivial (a_B grows with B), but a strong correlation with `abs_galactic_latitude` would point to remaining foreground contamination, and a strong correlation with `inclination_deg` to deprojection systematics.

## 6. Correlations vs theoretical_gap_coverage_ratio (entire sample)

| column | Pearson | Spearman | N |
|---|---:|---:|---:|
| galactic_latitude_degrees | -0.263 | -0.267 | 158 |
| abs_galactic_latitude | -0.136 | -0.028 | 158 |
| inclination_deg | -0.019 | 0.178 | 158 |
| magnetic_field_microgauss | 0.221 | 0.318 | 158 |
| plasma_beta_median | -0.040 | -0.535 | 158 |
| density_clamped_fraction | 0.237 | 0.494 | 158 |
| stellar_dominated_fraction | -0.012 | 0.243 | 158 |
| rar_residual_dex | -0.048 | -0.030 | 158 |
| physical_neighbor_separation_kpc | 0.161 | 0.194 | 158 |
| distance_mpc | 0.168 | 0.201 | 158 |

## 7a. RAR consistency on the reliable sample

- N = 28  median |residual| = 0.162 dex  MAD-robust σ = 0.109 dex

## 7b. BTFR distance-to-ridge on the reliable sample

- N = 16  median |residual| = 0.317 dex
- |residual| > 0.5 dex: **4 galaxies (25.0%)**

## 8. Magnetic-field distribution on the reliable sample

- N = 28  median = 1.531 μG  mean = 5.484 μG  std = 9.303 μG
- p5 = 0.235 μG  p95 = 26.955 μG

See `graphs/magnetic_field_distribution.png` and `graphs/magnetic_field_distribution.md` for the histogram + rug, and `graphs/magnetic_field_vs_galaxy_scale.png` / `.md` for the B-field vs outer-radius view.

## 9. Supreme-row consistency cross-check

- Supreme rows (passed inclination + reliability): **28**
- Reliable summaries on disk: **28**  (mismatch ⇒ stale or interrupted run)

> *This report is regenerated automatically every run. To inspect raw inputs see the per-galaxy `*_summary.json` files under `results/galaxy_analysis/`.*
