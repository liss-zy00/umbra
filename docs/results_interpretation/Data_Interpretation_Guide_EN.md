# Data interpretation guide (E-G pipeline)

## Purpose

This document complements [`Scientific_Evidence_Narrative_EN.md`](Scientific_Evidence_Narrative_EN.md) with **short question-and-answer blocks** for the main headline columns in Discovery outputs and in `validation_scientific_evidence_report.csv`. It does **not** duplicate the full semantics of `audit_flags`; for that, see §5 of the narrative and the methodology §4.4 in [`Draft_Paper_EG_Model_EN.md`](../academic/Draft_Paper_EG_Model_EN.md).

**Cohorts.** `validation_scientific_evidence_report.csv` only contains *supreme rows* (reliable + inclination quality) that pass the Im / photometric-LSB morphology filter **and** the coverage-ratio window **0.2–5.0** (`build_direct_evidence_rows`). For the full `validation_is_reliable` cohort, use `validation_reliable.csv` + `validation_consistency_audit.md` (narrative §1).

---

## `required_density_for_match`

- **Question.** What ionised-plasma density would the magnetic-pressure channel need so that the gradient closes the kinematic gap shell by shell?
- **Physical meaning.** Shell-wise ratio of magnetic-energy-gradient scale to gap magnitude, expressed as a density \(\rho_{\mathrm{req}}\) that would make the ideal-MHD pressure-gradient ansatz self-consistent at that radius.
- **Where to find it.** `discovery_*.csv` columns; merged into `validation_scientific_evidence_report.csv` only for rows accepted by `build_direct_evidence_rows` (Im or photometric LSB + coverage ratio inside the configured window).
- **What a healthy value looks like.** Values **near** \(\sim 10^{-27}\,\mathrm{kg\,m^{-3}}\) *after* crossing beam diagnostics, `resolution_correction_factor_*`, and the anchor sweep support the diffuse-ionised (WIM) coupling story. **Population medians far above** that anchor are the usual empirical outcome in documented SPARC+Planck runs — read as tension between ansatz and data, not as a silent disproof of the formalism (see [`Scientific_Evidence_Narrative_EN.md`](Scientific_Evidence_Narrative_EN.md) §5.3).

---

## `convergence_ratio_to_wim`

- **Question.** How close is the inferred requirement to the empirical WIM reference density?
- **Physical meaning.** Ratio \(\rho_{\mathrm{req}} / \rho_{\mathrm{WIM}}\) (or equivalent headline statistic from the Discovery row).
- **Where to find it.** Discovery CSVs; merged column in `validation_scientific_evidence_report.csv`.
- **What a healthy value looks like.** Unity is a perfect match to the anchor; the cohort list in `discovery_cohort_rho_band.csv` is defined primarily from a dynamic band on `required_density_for_match` (see `discovery_is_elite` above), not from this ratio alone.

---

## `theoretical_gap_coverage_ratio`

- **Question.** How much of the gravitational deficit did the magnetic-pressure proxy close at the outer disk, at the best grid scenario?
- **Physical meaning.** Modulus ratio \(|\langle a_B\rangle_{\mathrm{edge}}| / |\langle \mathrm{gap}\rangle_{\mathrm{edge}}|\) after shell masking and scenario selection.
- **Where to find it.** Per-scenario metrics and the best row in `validation_scientific_evidence_report.csv`.
- **What a healthy value looks like.** Order-of-magnitude consistency with unity (e.g. within \([0.1, 10]\) in log space) keeps the ansatz falsifiable but not immediately ruled out at that scenario; values far outside that band indicate mismatch at the chosen \((\Delta x, f_{\mathrm{ion}})\).

---

## `pearson_correlation`

- **Question.** Do the **shapes** of the radial profiles agree, not only the edge magnitudes?
- **Physical meaning.** Pearson correlation between `magnetic_pressure_acceleration(r)` and `acceleration_gap(r)` over reliable shells at the scenario under test.
- **Where to find it.** Scenario tables and `validation_scientific_evidence_report.csv` for the best scenario.
- **What a healthy value looks like.** Positive correlation indicates profiles rise and fall together; near-zero or strongly negative values undermine a simple shared-radial driver interpretation even if a ratio at one radius looks acceptable.

---

## `delta_r_m` and `physical_pixel_scale_kpc`

- **Question.** At what physical step is the pressure-gradient estimate taken?
- **Physical meaning.** `delta_r_m` is the neighbour separation in SI used in the finite-difference estimate of \(\nabla u_B\); `physical_pixel_scale_kpc` reports the same scale in kpc for quick comparison to disk scale heights and scenario gradient knobs.
- **Where to find it.** Discovery rows; both appear in the merged `validation_scientific_evidence_report.csv` when Discovery data exist for that galaxy.
- **What a healthy value looks like.** In Discovery, \(\Delta r\) should match the physical Planck pixel scale when `Δx_phys` is available; large discrepancies between physical pixel scale and the scenario-only \(\Delta x\) are expected to shift `ρ_req` and should be read together with `audit_flags`.

---

## `magnetic_field_microgauss`

- **Question.** Is the inferred field strength in a literature plausibility band for the target class?
- **Physical meaning.** Order-of-magnitude field implied by the Stokes-based calibration path feeding \(u_B\).
- **Where to find it.** Per-galaxy summaries and Discovery evidence tables as exposed by the pipeline exporters.
- **What a healthy value looks like.** Values wildly inconsistent with typical \(\mu\mathrm{G}\) disk fields should trigger scrutiny together with calibration audits (e.g. `warn_b_band` / calibration fallbacks documented in the narrative).

---

## `plasma_beta_thermal_median`

- **Question.** Is the magnetic-pressure ansatz locally appropriate (thermal pressure vs magnetic pressure)?
- **Physical meaning.** Median plasma \(\beta\) of thermal gas pressure to magnetic pressure over contributing shells — a local MHD validity probe, not a substitute for full tension and anisotropic transport.
- **Where to find it.** Galaxy-level metrics in pipeline outputs where computed.
- **What a healthy value looks like.** Very high \(\beta\) everywhere suggests magnetic pressure is dynamically negligible at those radii; very low \(\beta\) everywhere may indicate an over-strong \(B\) or density floor effects — read next to `audit_flags`.

---

## `discovery_is_reliable`, `validation_is_reliable`, and `audit_flags`

- **Question.** Which reliability flag applies to which pipeline phase?
- **Physical meaning.** `discovery_is_reliable` is the Phase 1 (Discovery census) row flag: minimal gates (finite \(\rho_{\mathrm{req}}\), latitude / Planck / NED audits). `validation_is_reliable` is the galaxy-level outcome after the 9-scenario grid and the full shell / stellar / inclination / calibration gates (`validation_all.csv` and per-galaxy `*_summary.json` under `data_reliability_assessment`). Per-scenario grid CSV rows still expose a scenario-local `is_reliable` boolean derived from the same `audit_flags` vocabulary at that \((\Delta x, f_{\mathrm{ion}})\) cell. `audit_flags` remains the structured source of truth in every case.
- **Where to find it.** `results/step_discovery/discovery_all.csv` lists every processed galaxy; `discovery_reliable.csv` / `discovery_not_reliable.csv` split on `discovery_is_reliable`. After the grid, `results/step_validation/validation_all.csv` lists every galaxy that entered the grid with `validation_is_reliable` on each row (same schema as before, including comma-separated `audit_flags`, `reliability_failure_codes`, and human-readable `note`; filter in tooling rather than separate pass/fail CSVs). Grid tables: `results/grid_results/scenario_*.csv`.
- **What a healthy value looks like.** A galaxy may be Discovery-passed (enters the grid) yet validation-failed after stricter masking; treat population headline statistics separately for Discovery vs Validation cohorts.

---

## `discovery_is_elite` (report column)

- **Question.** How is the short “internal consistency” list defined after Discovery?
- **Physical meaning.** When the dynamic band is available, `discovery_is_elite` is true for Discovery-reliable rows whose `required_density_for_match` lies in `[0.7, 1.3]` times the **median** `required_density_for_match` over the reliable subset (see `ResultPersister.DISCOVERY_ELITE_BAND_*`). When that band cannot be formed, `discovery_is_elite` is false and `discovery_elite_band_source` is `none`.
- **Where to find it.** `results/step_discovery/discovery_cohort_rho_band.csv`; mirrored as `discovery_is_elite` on `validation_scientific_evidence_report.csv` for galaxies that pass the direct-evidence dynamic ρ band.
- **What a healthy value looks like.** Elite rows highlight coherent \(\rho_{\mathrm{req}}\) scales within the Discovery-passed population; they are not a substitute for `validation_is_reliable`.

---

## `beam_size_kpc` and resolution dilution (`resolution_correction_factor_*`)

- **Question.** How does finite beam size relate to the magnetic gradient step?
- **Physical meaning.** `beam_size_kpc` approximates the physical FWHM scale on the sky converted to kpc at the adopted distance. `resolution_correction_factor_discovery` stores \(k \approx D_{\mathrm{beam}} / L_{\mathrm{grad}}\) with \(L_{\mathrm{grad}}\) from the **physical** Planck pixel separation (Discovery). `resolution_correction_factor_grid` recomputes the same ratio using the **best-scenario** `scenario_gradient_scale_kpc` knob. The legacy `resolution_correction_factor` column copies the grid value for backward compatibility.
- **Where to find it.** Discovery CSVs and merged columns in `validation_scientific_evidence_report.csv`.
- **What a healthy value looks like.** Large \(k\) warns that the smoothed Stokes map may dilute gradients relative to the assumed disk step — read `convergence_ratio_to_wim` and `required_density_for_match` together with these columns, not in isolation.

---

## `theoretical_gap_coverage_ratio_p16` / `_p84`

- **Question.** What does the catalogue uncertainty on distance, inclination and `V_flat` do to the headline ratio?
- **Physical meaning.** Monte-Carlo band on `theoretical_gap_coverage_ratio` (`N = 300` samples, seed `42`). The propagator rescales the central ratio analytically as `R = R_0 \cdot (D_0/D_s) \cdot (\sin i_s / \sin i_0)^2 \cdot (V_0/V_s)^2`. The same multiplicative rescaling is applied to `convergence_ratio_to_wim` to populate `_p16` / `_p84` for the WIM convergence column.
- **Where to find it.** `validation_scientific_evidence_report.csv`.
- **What a healthy value looks like.** A symmetric band centred on the point estimate (`median ≈ point`) with width controlled by the dominant catalogue uncertainty (typically inclination). A wildly asymmetric band signals that one input is sitting near a non-linearity (e.g. very low `sin i`).

---

## `accel_pressure_only_m_s2`, `accel_full_lorentz_m_s2`, `theoretical_gap_coverage_ratio_full_lorentz`

- **Question.** Does the magnetic-tension term change the conclusion?
- **Physical meaning.** Side-by-side per-galaxy medians of (i) the pressure-only acceleration `a_B = -\nabla u_B / \rho_\mathrm{ion}` and (ii) the same plus the toy tension `a_\mathrm{tension,r}(r) = -\cos p \cdot B^2 / (\mu_0 \rho_\mathrm{ion} R_c(r))` under a logarithmic-spiral geometry with pitch `p = 12°` and `R_c = r / \sin^2 p`. The parallel coverage ratio is reported as `theoretical_gap_coverage_ratio_full_lorentz`.
- **Where to find it.** `validation_scientific_evidence_report.csv`.
- **What a healthy value looks like.** For spiral hosts the two ratios should differ by at most order unity — the tension term is a known second-order effect. For irregular / BCD hosts the engine forces tension to zero (`warn_pure_pressure`) and the two columns coincide.

---

## `rar_median_abs_residual_dex`

- **Question.** Is the kinematics self-consistent with the Radial Acceleration Relation?
- **Physical meaning.** Median of `|log10(g_obs / g_pred)|` over reliable shells against the McGaugh+2016 RAR. A value `\lesssim 0.1` dex sits inside the published intrinsic scatter; a value `\gtrsim 0.3` dex flags a galaxy whose kinematics disagree with the empirical RAR locus.
- **Where to find it.** `validation_scientific_evidence_report.csv`; population median + MAD reported in `results/step_validation/validation_consistency_audit.md` (§7a).
- **What a healthy value looks like.** Most SPARC galaxies have `\le 0.13` dex.

---

## `btfr_residual_dex` and `warn_btfr_outlier`

- **Question.** Does the galaxy sit on the Baryonic Tully–Fisher ridge?
- **Physical meaning.** Signed distance (dex) of the galaxy from the McGaugh+2012 ridge `log10(M_bar / M_sun) = 4.0 \cdot log10(V_flat) + 1.99`, with `M_bar = 1.33 \cdot M_HI + \Upsilon_{3.6} \cdot L_{3.6}` and `\Upsilon_{3.6} = 0.5 M_sun/L_sun`. `warn_btfr_outlier` is true when `|residual| > 0.5` dex.
- **Where to find it.** `validation_scientific_evidence_report.csv`; population outlier fraction in `results/step_validation/validation_consistency_audit.md` (§7b).
- **What a healthy value looks like.** `|residual| \lesssim 0.2` dex.

---

## `warn_no_pa` and `warn_pure_pressure`

- **Question.** Which assumptions did the engine fall back on for this galaxy?
- **Physical meaning.** `warn_no_pa` means the curated supplement (`sources/morphology/disk_position_angle_overlay.json`) does not list a position angle for the galaxy; the engine recomputes `∇u_B` in the axisymmetric (sky-plane) frame with a `cos i` correction. `warn_pure_pressure` means the morphology is Irregular / BCD: the toy magnetic-tension term is forced to zero so `accel_full_lorentz_m_s2 = accel_pressure_only_m_s2`.
- **Where to find it.** `validation_scientific_evidence_report.csv`.
- **What a healthy value looks like.** Both flags are non-blocking — they only carry the information that the auxiliary refinement could not be exercised.

---

## `validation_anchor_sensitivity_sweep.csv`

- **Question.** How much does the headline ratio depend on the calibration anchor `u_ref`?
- **Physical meaning.** Two ratios per galaxy: `ratio_cmb` (the default CMB-anchored ratio) and `ratio_turbulent_equipartition` (the same ratio with the calibrator re-anchored to a turbulent equipartition pressure floor `u_turb = 0.5 \cdot \rho_gas \cdot \sigma_v^2`). `ratio_log10_delta = log10(ratio_turbulent / ratio_cmb)` captures the dex shift between the two anchors.
- **Where to find it.** `results/step_validation/validation_anchor_sensitivity_sweep.csv`.
- **What a healthy value looks like.** `|ratio_log10_delta| \lesssim 0.5` dex — the two anchors should agree within an order of magnitude.

---

## Figure gallery (publication bundle)

- **When it runs.** After **Phase 2 (Validation)** completes for the run (grid + per-galaxy artefacts), the runner calls :mod:`src.reporting.publication_bundle` (**Phase 3 — Graphs**). The only CLI entry point remains ``python src/main.py --galaxies …``; there is no separate figure subcommand.
- **CLI flags.** Galaxy selection uses ``--galaxies`` (see `src/pipeline/cli.py`). Optional ``--logs`` prints verbose per-galaxy provenance and **disables** tqdm progress bars so logs stay readable.

- **Galaxy list.** Per-galaxy PNGs are built for every object that entered the validation grid in that run (the same set as the ``*_rotmod.dat`` files that passed the Discovery pre-filter). Shared population panels read the freshly written ``results/step_*/*.csv`` tables on disk.

- **Output layout.** Per-galaxy PNGs (`*_acceleration_comparison.png`, `*_rar_plot.png`, `*_signed_accel_vs_radius.png`) are emitted only for galaxies listed in `step_validation/validation_reliable.csv`, next to each galaxy's CSV/JSON under `results/galaxy_analysis/{validation_scientific_evidence|validation_reliable}/<NAME>/` after a full run (the direct-evidence cohort is moved under `validation_scientific_evidence/`). Shared population panels sit in `results/graphs/population/`.

- **Markdown policy.** By default **no** `.md` sidecars are written next to these PNGs: axes, titles, and caption text on the figure carry the physics. Long-form interpretation stays in this folder (`Scientific_Evidence_Narrative_*.md`, this guide).

- **Typography and scales.** The job applies serif-forward Matplotlib defaults (DejaVu Serif, no external LaTeX requirement). Log-axis limits for the gap–magnetic comparison and for the RAR scatter use quantiles **across the galaxies processed in that run** so multi-panel thesis figures share comparable dynamic range.

- **Tests / alternate cohort panels.** The helper ``run_publication_bundle(..., cohort=…)`` writes population PNGs under ``results/graphs/population_<cohort>/`` when invoked from unit tests that load a specific cohort CSV.

---

## Where to look first

- **`results/step_validation/validation_consistency_audit.md`** — population summary and `audit_flags` histograms (regenerated every run).
- **`discovery_cohort_rho_band.csv`** — internal-consistency subset after the dynamic \(\rho_{\mathrm{req}}\) band.
- **`discovery_reliable.csv`** — Discovery `discovery_is_reliable` cohort for median \(\rho_{\mathrm{req}}\) and ratio statistics.
- **`discovery_all.csv` / `discovery_not_reliable.csv`** — full census and explicit complement.
- **`validation_scientific_evidence_report.csv`** — cross-phase view: Discovery headline fields **plus** best-scenario gap coverage and correlation.
- **`results/graphs/population/`** — cohort figures from the automatic publication step at the end of a full run (see Figure gallery above).
- **Per-galaxy artefacts** — `results/galaxy_analysis/{validation_scientific_evidence|validation_reliable|validation_not_reliable}/<NAME>/` holds CSV and JSON; publication PNGs follow the reliable cohort (and `validation_scientific_evidence/` after the run's relocation step).

---

## Related reading

- [`Scientific_Evidence_Narrative_EN.md`](Scientific_Evidence_Narrative_EN.md) — stable reading guide, plot semantics, and `audit_flags` discussion.
- [`Scientific_Evidence_Narrative_PT.md`](Scientific_Evidence_Narrative_PT.md) — Portuguese mirror.
