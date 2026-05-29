# How to Read the E-G Pipeline Results

> **Purpose.** This document is the **stable reading guide** for the
> artefacts produced by the pipeline (the per-galaxy CSV summary, the
> publication figure bundle under `results/graphs/population/` and co-located per-galaxy PNGs under `results/galaxy_analysis/`, and the per-galaxy `audit_flags` list). The module
> [`src/reporting/narrative.py`](../../src/reporting/narrative.py) can generate Markdown prose from reliable rows, but **the main CLI (`python src/main.py`) does not write** that file under `results/` by default. **This guide**
> explains *what each column means*, *how to read each plot*, and *what
> conclusions are licensed by each combination of values*. It deliberately
> avoids freezing per-galaxy numbers from any single execution.
> **Portuguese version:** [`Scientific_Evidence_Narrative_PT.md`](Scientific_Evidence_Narrative_PT.md).
> **Column Q&A (headline fields):** [`Data_Interpretation_Guide_EN.md`](Data_Interpretation_Guide_EN.md) / [`Data_Interpretation_Guide_PT.md`](Data_Interpretation_Guide_PT.md).

---

## 1. Which files to open first (reading order)

The pipeline writes many artefacts. For **statistics on the full validation cohort**
(`validation_is_reliable`), read in this order:

1. **`results/step_validation/validation_consistency_audit.md`** — `audit_flags` histograms, cohort RAR/BTFR summaries, morphology splits, and bias probes (regenerated every run).
2. **`results/step_validation/validation_reliable.csv`** — one row per galaxy with `validation_is_reliable == True` (best scenario already folded into the inventory row).
3. **`results/step_validation/validation_all.csv`** — every galaxy that entered the grid (reliable or not).

**`results/step_validation/validation_scientific_evidence_report.csv` does *not* replace step 2.** It is a **direct-evidence slice** built in [`build_direct_evidence_rows`](../../src/pipeline/aggregation/global_aggregator.py) from *supreme rows* (galaxies that are `validation_is_reliable` **and** pass the processor inclination-quality gate). A galaxy is written only if it additionally satisfies:

- morphology **Im (Dwarf Irregulars)** **or** photometric LSB (threshold `RELIABILITY_CUTS.PHOTOMETRIC_LSB_SB_THRESHOLD`; see [`include_in_scientific_evidence_report`](../../src/pipeline/per_galaxy/morphology.py));
- a finite `theoretical_gap_coverage_ratio` inside **[0.2, 5.0]** (`RELIABILITY_CUTS.DIRECT_EVIDENCE_RATIO_MIN/MAX` in [`thresholds.py`](../../src/config/thresholds.py)).

Expect **far fewer rows** than `validation_reliable.csv`. Use this CSV for Im/LSB synthesis and population panels that read it; reviewers who need the **entire** reliable cohort should start from `validation_reliable.csv` + the audit.

- **`results/graphs/population/`** — cohort diagnostic PNGs emitted at the end of a successful validation run (built from the on-disk `step_*` tables and `galaxy_analysis` metadata).

## 2. Columns of `validation_scientific_evidence_report.csv`

Each CSV row is produced from a *supreme row* (§1), enriched in [`build_direct_evidence_rows`](../../src/pipeline/aggregation/global_aggregator.py), and written by [`save_scientific_evidence_report`](../../src/reporting/persistence/exporter.py). The full column schema (including `validation_status`, `wim_density_match_score`, beam diagnostics) lives in the exporter — the subsections below focus on the highest-physics-weight fields.

The columns described here are the ones that carry physical meaning for
the test. They are produced by
[`src/pipeline/per_galaxy/metrics_builder.py`](../../src/pipeline/per_galaxy/metrics_builder.py)
and persisted by
[`src/reporting/persistence/exporter.py`](../../src/reporting/persistence/exporter.py).

The table also carries **Discovery headline fields** merged from Phase 1
(`required_density_for_match`, `convergence_ratio_to_wim`, `delta_r_m`,
`physical_pixel_scale_kpc`, `discovery_is_reliable`, `discovery_is_elite`);
see [`Data_Interpretation_Guide_EN.md`](Data_Interpretation_Guide_EN.md) for
compact Q&A on those columns.

### 2.1 Identification

- `galaxy_name` — SPARC catalogue name.
- `morphology_class` — Hubble-stage class group from the Lelli+2016 JSON
  lookup. Read this together with the coverage ratio: dwarfs and LSBs are
  the sample where the magnetic-pressure hypothesis is *expected to matter
  most*, because the gravitational baryonic contribution there is small.

### 2.2 Best-scenario parameters

- `scenario_gradient_scale_kpc` — the value of Δx (in kpc) used to convert
  the HEALPix neighbour difference `Δu_B` into a gradient. Three values
  are explored: 0.050, 0.100 and 0.500.
- `scenario_ionization_fraction` — the WIM ionisation fraction `f_ion`
  multiplying `ρ_gas` to obtain the responsive density `ρ_ion`. Three
  values: 0.001, 0.010, 0.100.

A scenario where the best fit prefers small `scenario_gradient_scale_kpc`
(e.g. 50 pc) and small `scenario_ionization_fraction` (e.g. 0.1 %) is
physically consistent with a **filamentary magnetic geometry coupled to
the tenuous WIM**. A best fit at large `scenario_gradient_scale_kpc` (e.g.
500 pc) and large `f_ion` (e.g. 10 %) is consistent with a **smoothly
distributed disk field acting on a denser ionised layer**. Both are
permissible in principle, but they point at different physical regimes
that subsequent refinements (anchor sensitivity, magnetic tension) must
disambiguate.

### 2.3 The two headline numbers

- `theoretical_gap_coverage_ratio` —
  `|⟨a_B⟩_edge| / |⟨gap⟩_edge|`. The natural quantity to look at is its
  base-10 logarithm. Values inside `[10^{-1}, 10^{+1}]` (i.e. within one
  order of magnitude of unity) mean the magnetic-pressure proxy is of the
  right order of magnitude to close the kinematic deficit at the outer
  disk; values outside that band falsify the ansatz at the scenario it
  was computed under.
- `pearson_correlation` — Pearson correlation coefficient between the
  radial profiles of `magnetic_pressure_acceleration` and
  `acceleration_gap` over reliable shells. A positive value indicates that
  the two profiles *increase together*: not just that the magnitudes agree,
  but that the **shapes** agree. A negative value, or a value statistically
  consistent with zero, weakens the case even if the modulus ratio is
  close to one.

### 2.4 Physical sanity columns

- `magnetic_field_microgauss` — median inferred `|B|` for the galaxy, in
  µG. Values in the band `[0.1, 100]` are consistent with the
  literature for spiral and dwarf disk fields. Values outside that band
  raise an informational `warn_b_band` audit flag (not blocking).
- `plasma_beta_thermal_median` — median plasma-β over reliable shells.
  Values `< 1` mean magnetic pressure dominates thermal pressure (the
  regime where the ideal-MHD pressure-gradient ansatz is locally
  consistent); `> 1` means thermal pressure dominates. A high `β` is *not*
  by itself a failure — it just signals that the diagnostic is locally
  out of regime and should not be taken as confirmatory for that galaxy.

### 2.5 Audit columns

- `audit_flags` — string list carried from the supreme-row snapshot; empty means **no blocking `fail_*` codes** on the validation inventory (`warn_*` remain informational). Exact codes: §4.
- `validation_status` — label from [`_validation_status_supreme`](../../src/pipeline/aggregation/global_aggregator.py): `SUPREME` (high-provenance subset), `DATA_LIMITED` (Planck beam larger than the disk-diameter proxy), or `STANDARD`. Do **not** confuse with the inventory column `validation_is_reliable` (every row here already passed reliability + inclination before the morphology/ratio filters).

On per-scenario grid CSVs under `results/grid_results/`, the boolean **`is_reliable`** still means “no `fail_*` in that \((\Delta x, f_{\mathrm{ion}})\) cell”.

### 2.6 Phase-1 refinement columns

- `theoretical_gap_coverage_ratio_p16`, `theoretical_gap_coverage_ratio_p84`
  — 16th and 84th percentiles of the Monte-Carlo band on the headline
  coverage ratio. The MC propagates the catalogue 1-σ uncertainties on
  distance `D`, inclination `i` and asymptotic flat rotation velocity
  `V_flat` through the closed-form multiplicative rescaling
  `R(D_s, i_s, V_s) = R_0\,(D_0/D_s)\,(\sin i_s/\sin i_0)^2\,(V_0/V_s)^2`,
  drawing `N = 300` samples with seed `42` (deterministic).
- `convergence_ratio_to_wim_p16`, `convergence_ratio_to_wim_p84` — analogous
  band on the WIM convergence ratio (Discovery census column scaled by the
  same multiplicative kinematic factor).
- `theoretical_gap_coverage_ratio_full_lorentz` — coverage ratio when the
  toy magnetic-tension term (under a logarithmic-spiral geometry, pitch
  angle 12°) is added to the pressure-only acceleration. The headline
  `theoretical_gap_coverage_ratio` keeps the pressure-only convention; this
  column reports the parallel "full Lorentz" budget so reports can compare
  the two side by side.
- `accel_pressure_only_m_s2`, `accel_full_lorentz_m_s2` — per-galaxy median
  of the pressure-only and full-Lorentz acceleration profiles, in SI.
- `rar_median_abs_residual_dex` — median `|log10(g_obs / g_pred)|` over
  reliable shells against the McGaugh+2016 Radial Acceleration Relation
  (RAR). The population statistic over the reliable sample is surfaced in
  `results/step_validation/validation_consistency_audit.md` (Section 7a).
- `btfr_residual_dex` — signed distance (dex) of the galaxy from the
  McGaugh+2012 Baryonic Tully–Fisher ridge. `warn_btfr_outlier` flips when
  `|residual| > 0.5` dex. The population outlier fraction appears in
  `results/step_validation/validation_consistency_audit.md` (Section 7b).
- `warn_no_pa`, `warn_pure_pressure`, `warn_btfr_outlier` — Phase-1 informational audit tags (non-blocking; see §4).

The companion file
`results/step_validation/validation_anchor_sensitivity_sweep.csv` lists, per galaxy,
`ratio_cmb` (the default CMB-anchored ratio) and
`ratio_turbulent_equipartition` (the ratio recomputed against a turbulent
equipartition pressure floor), together with `ratio_log10_delta`. Use this
file to audit the headline ratio against the choice of dimensional anchor.

## 3. Scenario grid vs headline row

During **Phase 2 (Validation)** the engine evaluates the full discrete \((\Delta x, f_{\mathrm{ion}})\) grid for each galaxy that passed the Discovery pre-filter. For the **entire** reliable cohort, the winning cell lives in each `*_summary.json` / inventory row (`validation_reliable.csv`). A row in **`validation_scientific_evidence_report.csv`** exists only for the Im/LSB + coverage-ratio subset (§1) and already collapses that search: `best_scenario_name`, `scenario_gradient_scale_kpc`, and `scenario_ionization_fraction` identify the winning cell, while `theoretical_gap_coverage_ratio` and `pearson_correlation` describe how well that cell closes the gap in modulus and in shape.

The pipeline **no longer** writes a standalone 3×3 PNG heatmap of the grid. For a population view of gap coverage, Pearson \(r\), Discovery–Validation cross-checks, and resolution diagnostics, open **`results/graphs/population/`** after a full run (see [`Data_Interpretation_Guide_EN.md`](Data_Interpretation_Guide_EN.md), Figure gallery).

### How to read the headline scenario fields

- **Coverage near unity** on the best row means the magnetic-pressure proxy matched the kinematic deficit in order of magnitude for that \((\Delta x, f_{\mathrm{ion}})\) choice — necessary but not sufficient for a positive interpretation.
- **Coverage \(\ll 1\)** means the proxy is systematically too small at the chosen knobs.
- **Coverage \(\gg 1\)** means the proxy overshoots — often a sign that \(\Delta x\) is too aggressive or \(f_{\mathrm{ion}}\) is too low for the assumed ISM model.
- **Pearson \(r\)** near zero with good coverage would indicate a modulus match without a credible radial shape match; treat that combination as physically suspect even if the ratio looks acceptable.

The discrete grid alone does **not** confirm the hypothesis: it only maps where the modulus test is consistent. The cohort panels under `graphs/population/` complement the CSV by showing how those headline fields distribute across the processed sample.

## 4. Reading `audit_flags`

Every grid galaxy carries `audit_flags` as the structured gate vocabulary. Rows in
**`validation_scientific_evidence_report.csv`** are supreme-row snapshots (reliable +
inclination quality) **after** the morphology and coverage-ratio filters from §1; they
may still include informational `warn_*` tags. The canonical code list is in
[`src/pipeline/per_galaxy/metrics_builder.py`](../../src/pipeline/per_galaxy/metrics_builder.py)
and the matching section of [`README.md`](../../README.md) ("Reliability audit").

| Code                       | Meaning                                                                                              | Blocking? |
| -------------------------- | ---------------------------------------------------------------------------------------------------- | --------- |
| `fail_lat`                 | Galactic latitude inside the foreground-contamination band.                                          | yes       |
| `fail_reliable_fraction`   | Too few radial shells survived the per-shell MAD outlier gate on `magnetic_pressure_acceleration`.   | yes       |
| `fail_density_clamped`     | Too many shells fell back to the density floor in `GasDensityEstimator`.                             | yes       |
| `fail_stellar`             | Too many shells dominated by `V_disk`/`V_bulge`; gas-mass inversion is not a meaningful tracer.       | yes       |
| `fail_skipped_physics`     | The engine could not produce `a_B` (missing Stokes, Δx, or `f_ion`).                                  | yes       |
| `fail_inclination`         | Catalogued inclination below the configured floor; `V_obs / sin i` propagation destabilises.          | yes       |
| `fail_calibration`         | The Stokes reference fell back to the emergency value; the inferred `B` is non-measured.             | yes       |
| `warn_b_band`              | Median \|B\| outside the 0.1–100 µG plausibility band; informational only.                            | no        |
| `warn_no_pa`               | No curated disk position angle in `sources/morphology/disk_position_angle_overlay.json`; the engine fell back to the axisymmetric gradient. | no        |
| `warn_pure_pressure`           | Irregular / BCD morphology; the magnetic-tension toy collapses to zero so `accel_full_lorentz_m_s2 = accel_pressure_only_m_s2`. | no        |
| `warn_btfr_outlier`        | Galaxy sits more than 0.5 dex from the McGaugh+2012 BTFR ridge.                                       | no        |

On **`validation_all.csv` / `validation_reliable.csv`**, `validation_is_reliable` is `True`
iff no `fail_*` code is present (`warn_*` never block). The JSON summary under each
`galaxy_analysis/.../<NAME>/` folder mirrors the same gate list; **`results/step_validation/validation_consistency_audit.md`**
aggregates how often each code appears across the full processed sample so
systematic biases in the rejected cohort stay visible without mixing them
into the headline reliable statistics.

## 5. Data Cleaning and Noise

Scientific reading rests on three operational exclusion rules plus the physical
interpretation of `required_density_for_match` in **Phase 1 (Discovery)**. The first three
subsections cover foreground, noisy shells, and the Discovery census; the population median
of \(\rho_{\mathrm{req}}\) must **not** be read in isolation without the beam diagnostics and
anchor sweep (§5.3).

### 5.1 Galactic-latitude cut ($|b| < 20^\circ$)

The pipeline emits the audit code `fail_lat` and rejects the galaxy whenever
the galactic latitude satisfies $|b| < 20^\circ$. The canonical threshold lives
in
[`RELIABILITY_CUTS.MIN_GALACTIC_LATITUDE_DEG = 20.0`](../../src/config/thresholds.py).
The motivation is instrumental, not cosmological: the equatorial band of the
Milky Way concentrates polarised dust whose thermal emission saturates the
Stokes signal in the Planck band used as input (353 GHz). At low latitudes,
$\sqrt{Q^2 + U^2}$ measured at the central pixel of the target galaxy is
dominated by the galactic foreground rather than by the magnetic field of the
external object. The consequence is direct on the calibrator: the relation

\[
B(\hat{n}) \;=\; \sqrt{2\, \mu_0\, u_{\mathrm{CMB}}}\,\frac{S(\hat{n})}{S_{\mathrm{ref}}}
\]

stops measuring the target's field and starts measuring the Milky Way disk
field along the line of sight. Retaining such galaxies would produce a high
systematic bias in `B`, `u_B`, and — by construction — in `a_B`, mimicking a
positive E-G signal for a purely geographical reason. The cut is therefore
part of the experimental design.

### 5.2 "Noisy" galaxy: operational criterion

A galaxy is classified as noisy through two complementary mechanisms, both
built on robust statistics over the radial-shell ensemble.

- **Per-shell physical outliers.** The `ReliabilityMaskBuilder` gate computes
  $\log_{10}(|a_B|/|a_c|)$ at each shell and rejects those falling outside
  $\mathrm{median} \pm N\,\sigma_{\mathrm{MAD}}$, with $N = 3.0$ and
  $\sigma_{\mathrm{MAD}} = 1.4826 \times \mathrm{MAD}$. The parameter $N$ is
  [`MAG_PRESSURE_RATIO_OUTLIER_SIGMA = 3.0`](../../src/config/thresholds.py).
  Shells outside this window are dropped *individually*: one anomalous shell
  never invalidates a whole galaxy.
- **Galaxy-level population filter.** If more than
  $1 - \texttt{MIN\_RELIABLE\_MEASUREMENT\_FRACTION} = 0.30$ of the shells fail
  the per-shell gate, the galaxy receives `fail_reliable_fraction`, is **not**
  `validation_is_reliable`, and never enters `state.supreme_rows` — so it cannot
  appear in `validation_scientific_evidence_report.csv` or `validation_reliable.csv`.
  The row remains in `validation_all.csv` for audit purposes.

The practical consequence is that a rotation curve contaminated by
beam-smearing, by unmodelled warps, or by local density spikes does not
disappear silently: it surfaces in the galaxy's audit flags and still
contributes to the "ghost layer" of the global scatter, keeping the
transparency about what was and was not rejected.

### 5.3 Column `required_density_for_match` and joint reading with WIM / resolution

During **Phase 1 (Discovery census)** the engine fixes $f_{\mathrm{ion}} = 1.0$ and evaluates
the magnetic-pressure gradient at the **physical Planck pixel step**
$\Delta r_{\mathrm{phys}} = D_{\mathrm{gal}}\cdot\theta_{\mathrm{HEALPix}}$,
rather than the exploratory `scenario_gradient_scale_kpc` knob used in validation. Each
galaxy contributes rows to the four Discovery artefacts below.

- [`results/step_discovery/discovery_all.csv`](../../results/step_discovery/discovery_all.csv)
  — every processed galaxy, with the exclusion reason preserved in
  `audit_flags`.
- [`results/step_discovery/discovery_reliable.csv`](../../results/step_discovery/discovery_reliable.csv)
  — the `discovery_is_reliable == True` subset; use this for population medians of
  \(\rho_{\mathrm{req}}\) over the clean Discovery cohort.
- [`results/step_discovery/discovery_not_reliable.csv`](../../results/step_discovery/discovery_not_reliable.csv)
  — the complement of `discovery_reliable` (same schema).
- [`results/step_discovery/discovery_cohort_rho_band.csv`](../../results/step_discovery/discovery_cohort_rho_band.csv)
  — the *internal-consistency* short list: Discovery-reliable rows whose
  `required_density_for_match` lies in the multiplicative band
  `[0.7, 1.3]` around the **median** `required_density_for_match` over the reliable subset
  (see `ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MIN/MAX` in
  [`exporter.py`](../../src/reporting/persistence/exporter.py)).
  When that dynamic band cannot be formed, there is no additional cohort: all dependent
  rows receive `discovery_is_elite = False` and `discovery_elite_band_source = none`.

The headline column `required_density_for_match` answers: what plasma density
\(\rho_{\mathrm{req}}\) would be required for the **reconstructed** magnetic-pressure
gradient to close the kinematic gap shell by shell?

\[
\rho_{\mathrm{req}}(r) \;=\;
\frac{\bigl|\,\nabla u_B(r)\,\bigr|}{\bigl|\,\Delta a(r)\,\bigr|}.
\]

The column `delta_r_m` stores $\Delta r_{\mathrm{phys}}$ in SI metres. The legacy column
`median_required_density_for_match_kg_m3` duplicates `required_density_for_match` for
backward compatibility.

`convergence_ratio_to_wim` compares \(\rho_{\mathrm{req}}\) to the fixed WIM reference
density \(10^{-27}\,\mathrm{kg\,m^{-3}}\) (order-of-magnitude diffuse ionised ISM; see
H\(\alpha\) literature). **In full SPARC+Planck runs documented after Phase 1**, the
population median of `required_density_for_match` over `discovery_is_reliable` is typically
**many orders of magnitude above** \(10^{-27}\,\mathrm{kg\,m^{-3}}\), so `convergence_ratio_to_wim ≫ 1`
is the usual empirical state, not an automatic “lock” onto the WIM anchor. That does **not**
invalidate the formalism; it forces joint reading with:

- `beam_size_kpc`, `resolution_correction_factor_discovery`, `resolution_correction_factor_grid`
  (finite-beam dilution of the gradient at the assumed reconstruction step);
- `validation_anchor_sensitivity_sweep.csv` (dex shift between the CMB dimensional anchor and a
  turbulent equipartition floor in the Stokes→Tesla bridge);
- the discussion of finite-beam dilution and calibration leverage in Section 5 of
  [`EG_MODEL_THEORY_EN.md`](../theory/EG_MODEL_THEORY_EN.md) and in the data interpretation guides under
  [`docs/results_interpretation/`](../results_interpretation/).

Physical reading: when \(\rho_{\mathrm{req}}\) **approaches** \(10^{-27}\,\mathrm{kg\,m^{-3}}\)
*after* these diagnostics, the WIM-coupling narrative remains plausible; when it stays high,
the number is **information about the ansatz + data**, not silent confirmation. Cold HI
(\(\rho_{\mathrm{HI}} \sim 10^{-21}\,\mathrm{kg\,m^{-3}}\)) still marks the naive “dense gas”
upper reference — the MHD channel in this formulation is not pretending to close the gap with
that component.

At the end of a run the pipeline prints the Discovery-reliable median of `required_density_for_match`,
the `N_passed` count, and the dynamic-band `N_elite` tally as a quick order-of-magnitude check.

### 5.4 Discovery as the empirical foundation of the validation grid

The project treats Discovery and the grid as sequential phases of a single argument,
not as isolated runs:

1. **Discovery census** answers *"what density does the magnetic signal demand?"* with no
   grid knobs — Δx is the physical Planck pixel and $f_{\mathrm{ion}}=1$. Absolute
   \(\rho_{\mathrm{req}}\) values should be confronted with the WIM **as a reference**, not as
   a guaranteed population-median target (§5.3).
2. **Grid search (Phase 2 — Validation)** asks how gap closure changes when effective gradient scale and
   ionised fraction in \(\rho_{\mathrm{ion}}\) are allowed to vary. The discrete
   `scenario_gradient_scale_kpc ∈ {0.050, 0.100, 0.500}` kpc and
   `scenario_ionization_fraction ∈ {0.001, 0.010, 0.100}` values are an explicit prior on
   small-scale magnetic structure and ionisation; they do not replace Discovery.

In practice `discovery_reliable.csv` sets the order-of-magnitude context before the grid, and
`discovery_cohort_rho_band.csv` lists galaxies whose \(\rho_{\mathrm{req}}\) is internally
consistent with the Discovery-reliable population median.

**Do not confuse with validation cohorts:** filter `step_validation/validation_all.csv` on
`validation_is_reliable` for post-grid headline statistics (typically a strictly smaller universe
than Discovery-reliable). Merged validation tables still expose `discovery_is_reliable` /
`discovery_is_elite` with **Discovery** semantics.

Changing the definition of \(\rho_{\mathrm{req}}\) in §4.2 of
[`Draft_Paper_EG_Model_EN.md`](../academic/Draft_Paper_EG_Model_EN.md), or
the elite-band constants (`DISCOVERY_ELITE_BAND_*`), requires updating this
paragraph and its PT mirror in the same commit.

## 6. What a coverage ratio close to unity does *not* prove

Even if `theoretical_gap_coverage_ratio ≈ 1` for many reliable galaxies
and the Pearson correlation is positive, the result is not yet a
displacement of the dark-matter interpretation. The following limitations
must be cleared before the claim is statistically robust (see Section 5 of
[`EG_MODEL_THEORY_EN.md`](../theory/EG_MODEL_THEORY_EN.md) and Section 6 of this narrative); several already have **partial mitigations** in the current code — the open question is whether they satisfy an external referee:

1. **Disk vs sky gradient geometry.** The engine already estimates `||∇u_B||` with a disk-plane LSQ fit when PA and inclination exist (`disk_plane_pressure_gradient`); when PA is missing it falls back to an axisymmetric sky-plane path (`warn_no_pa`). Mixing HEALPix neighbour geometry, `scenario_gradient_scale_kpc`, and disk physics nevertheless remains a reading limitation — see [`PhysicsEngine`](../../src/physics/engine.py) caveats.
2. **Uncertainty on the headline ratio.** Catalogue uncertainties on \(D\), \(i\), and \(V_{\mathrm{flat}}\) are already propagated into `theoretical_gap_coverage_ratio_p16` / `_p84` via closed-form Monte Carlo. This is **not** a full Bayesian posterior and does not propagate correlated Stokes noise or local \(B\) calibration errors.
3. **Stokes→Tesla anchor.** `validation_anchor_sensitivity_sweep.csv` quantifies the dex shift between the default CMB dimensional anchor and a turbulent equipartition floor; it does **not** implement observational synchrotron-based anchors per galaxy.
4. **Magnetic tension.** A logarithmic-spiral toy tension branch exists (`theoretical_gap_coverage_ratio_full_lorentz`) but is disabled for irregular/BCD hosts (`warn_pure_pressure`). Predictive tension still needs sampled 3D geometry beyond the toy.
5. **Synthetic controls.** NFW-style controls bound only specific false-positive classes; they do not exclude all non-MHD mechanisms or instrumental systematics.

A near-unity coverage ratio is therefore a **necessary** but not
**sufficient** outcome: it licenses further work; it does not by itself
overturn ΛCDM.

## 7. What a failed coverage ratio *does* prove

If the headline numbers come back orders of magnitude away from unity
and the Pearson correlation is consistent with zero across reliable
galaxies, the conclusion is **definite for that scenario**: the
ideal-MHD pressure-gradient ansatz with the configured anchor and the
configured Δx/`f_ion` is not the mechanism behind the rotation-curve
anomaly. A negative result of this kind is, in scientific terms, just
as valuable as a positive one — it removes a candidate from the
hypothesis space without requiring an external experiment.

## 8. Where to look next

After `validation_reliable.csv` and the audit, use
`validation_scientific_evidence_report.csv` only when the scientific argument is explicitly
Im/LSB **and** the 0.2–5.0 coverage-ratio window matters.

- **`results/galaxy_analysis/{validation_scientific_evidence|validation_reliable|validation_not_reliable}/<NAME>/`** — radial acceleration comparison, RAR, and signed-acceleration panels for each galaxy processed in that run (written next to the CSV/JSON tables on disk; the `validation_scientific_evidence` subfolder lists galaxies named in the direct-evidence CSV after the publication relocation step).
- **`results/step_validation/validation_consistency_audit.md`** — population-level audit
  including the frequency of every `audit_flags` code across the run, so
  that systematic biases in the rejected sample are visible.
- **Optional Markdown narrative** — `NarrativeGenerator` in `src/reporting/narrative.py` can render template prose from reliable rows, but the **CLI does not write** it under `results/` automatically; for durable reporting use this guide + CSVs + the audit.
