# E-G Rotation-Curve Pipeline — Software Guide

This repository implements a batch pipeline that fuses **SPARC** kinematics, **NED** coordinates, **Planck** polarization (HEALPix, memmap), and an **MHD magnetic-pressure** test of the **acceleration gap** per galaxy. Scientific motivation and equations are in [`docs/theory/EG_MODEL_THEORY_EN.md`](docs/theory/EG_MODEL_THEORY_EN.md) (PT: [`EG_MODEL_THEORY_PT.md`](docs/theory/EG_MODEL_THEORY_PT.md)); limitations of the ansatz are summarized in Section 5 of those theory files. For how to read pipeline outputs, see [`docs/results_interpretation/`](docs/results_interpretation/) (narratives and data guides). Theory files are not overwritten by runs.

## Pipeline logic (data flow)

End-to-end flow from disk to reports:

```text
SPARC *_rotmod.dat
    └─> SparcProvider.load_rotation_curve()
            → DataFrame (R, V_obs, V_gas, V_disk, V_bulge, V_bar)
            → Galaxy.rotation_dataframe

SPARC_Lelli2016c.json
    └─> SparcMorphologyProvider.load_morphology_lookup()
            → morphology_class, SBdisk (for LSB proxy)

NED (HTTP, optional)
    └─> NedProvider.fetch_metadata(name)
            → RA/Dec, z, distance hints → Galaxy fields

sources/morphology/disk_position_angle_overlay.json (optional)
    └─> merged into morphology lookup for disk-plane gradient / warn_no_pa

sources/extragalactic_distances/*.csv (optional)
    └─> EddProvider (bootstrap) → distance cross-checks when files exist

Planck polarization FITS (memmap)
    └─> PlanckProvider.initialize()  # headers, NSIDE, healpy readiness
    └─> PlanckProvider.extract_stokes_qu_mean(RA, Dec)
            → Stokes amplitudes (center + neighbor ring) → Discovery census + scenario grid

PhysicsEngine (src/physics/engine.py)
    └─> calculate_force_gap()   # SI accelerations, acceleration_gap
    └─> calculate_magnetic_pressure_acceleration(...)  # MHD a_B = -∇u_B/ρ_ion column

Pipeline orchestration (src/pipeline/)
    └─> scenario_grid_runner, scenario_selector, metrics_builder
    └─> galaxy_processor.process_sparc_file()

Reporting
    └─> ResultPersister / DataPersister → CSV, JSON, Markdown
    └─> publication_bundle → per-galaxy PNGs under `galaxy_analysis/` + cohort PNGs under `graphs/population/`
```

The entry point is [`src/main.py`](src/main.py): it adds the repo root to `sys.path` and calls `run_pipeline(parse_args())` from [`src/pipeline/runner.py`](src/pipeline/runner.py).

## How to run (unified pipeline)

Run from the **repository root** (paths such as `sources/` are relative to CWD). External scientific inputs live in [`sources/`](sources/), organised by provenance — `sources/sparc/` (Lelli+2016c rotation curves and morphology JSON), `sources/morphology/` (curated disk position angle supplement for deprojection), `sources/planck/` (Planck Legacy Archive thermal-dust polarisation maps), `sources/ned/` (NED/SIMBAD name-resolution config), and `sources/extragalactic_distances/` (EDD-style distance tables consumed by `edd_provider` when those CSVs are present). Each subdirectory has its own `README.md` with the canonical citation and the consumer module.

Every `python src/main.py` invocation runs **both** phases in order: (1) **Discovery census** on the full `--galaxies` target set (`f_ion = 1.0`, physical Planck pixel scale, four `step_discovery/*.csv` files), then (2) **automatic pre-filter** (galaxies with blocking Discovery audits `fail_lat`, `fail_skipped_physics`, or `fail_ned_resolution` skip the grid), then (3) **9-scenario grid search** and full validation artefacts only for galaxies that passed the pre-filter. The console ends with one summary block. Column semantics: [`docs/results_interpretation/Data_Interpretation_Guide_EN.md`](docs/results_interpretation/Data_Interpretation_Guide_EN.md) and [`Data_Interpretation_Guide_PT.md`](docs/results_interpretation/Data_Interpretation_Guide_PT.md).

```bash
conda activate dark-matter   # or: conda run -n dark-matter --no-capture-output python ...

# ``healpy`` must be importable in the active environment: when a Planck FITS is found,
# ``PlanckProvider.initialize()`` raises ``RuntimeError`` if healpy is missing (avoids silent NaN runs).

# ------------------------------------------------------------------
# UNIFIED PIPELINE — Discovery census + grid (pre-filtered) + reports
# ------------------------------------------------------------------
python src/main.py                                              # ALL galaxies (~175)
python src/main.py --galaxies default                           # 10 control galaxies
python src/main.py --galaxies NGC3741                           # ONE galaxy (debug)
python src/main.py --galaxies NGC3741 DDO154 NGC2403            # explicit list
python src/main.py --logs --galaxies NGC3741                   # verbose per-galaxy logs (disables tqdm redraws)
```

- **`--galaxies`**: which galaxies to process. Use the literal `all` (full SPARC folder scan), `default` (subset matched against `DEFAULT_CONTROL_GALAXIES` after name normalization), or one or more explicit galaxy names (e.g. `NGC3741 DDO154`). Presets cannot be mixed with names. Defaults to `all`.
- **`--logs`**: print per-galaxy provenance (SPARC / morphology / NED / Planck, Unicode rule banner) and verbose provider chatter; **disables tqdm** progress bars so log lines are not interleaved with bar redraws.

### Optional 3D toy-disk renders (PyVista)

Pedagogical PNGs from per-galaxy ``*_analysis.csv`` (graphic composition — not N-body / not self-consistent Poisson+MHD). Use either the stand-alone module or optional flags on the unified pipeline:

```bash
python -m src.visualization --galaxies DDO154 --project-root .
python -m src.visualization --galaxies auto --max-n 8 --tier reliable --viz-comparison --wim-rings --gap-glyphs
python src/main.py --galaxies DDO154 --viz-galaxies DDO154 --viz-wim-rings --viz-gap-glyphs
```

Flag naming: the stand-alone module uses ``--gap-glyphs``, ``--wim-rings``, ``--streamlines``, ``--manifest`` (no ``viz_`` prefix). The unified ``python src/main.py`` mirror flags are prefixed with ``--viz-`` (e.g. ``--viz-gap-glyphs``).

- **`python -m src.visualization`**: ``--galaxies auto`` reads ``validation_reliable.csv`` or, with ``--tier evidence``, ``validation_scientific_evidence_report.csv``; ``--max-n`` caps the list (default ``12`` from ``VIZ_DISPLAY`` in ``src/config/thresholds.py``). ``--gap-glyphs``, ``--wim-rings``, ``--streamlines``, ``--manifest`` mirror the in-process renderer.
- **`--viz-galaxies` … on `python src/main.py`**: runs **after** README regeneration when PyVista is importable; otherwise prints a skip notice. Pair with ``--viz-comparison``, ``--viz-wim-rings``, ``--viz-gap-glyphs``, ``--viz-streamlines``, ``--viz-tier``, ``--viz-max``, ``--viz-manifest``, ``--viz-output-dir`` as needed.

## SOLID-oriented layout

| Layer | Responsibility | Key modules |
|-------|----------------|-------------|
| **Providers** | I/O to external or packaged data only; no science policy | `src/providers/sparc.py`, `planck.py`, `ned.py`, `edd_provider.py` |
| **Physics** | SI kinematics, MHD acceleration, statistics | `src/physics/engine.py`, `statistics.py` |
| **Pipeline** | Orchestration, thresholds, scenario grid, CSV shapes | `src/pipeline/*.py` |
| **Reporting** | Persistence and plots; no change to physics numerics | `src/reporting/persistence.py`, `exporter.py`, `plotting.py`, `narrative.py` |

- **Single Responsibility**: e.g. `PlanckProvider` never writes `results/`; `DataPersister` never queries NED.
- **Open/Closed**: new providers can be swapped via `PipelineDependencies` in tests.
- **Dependency inversion**: the pipeline depends on provider interfaces, not on FITS layout details (those stay inside `PlanckProvider`).

## Physical basis (summary)

At each radius, `centripetal_acceleration` from `V_obs` exceeds `baryonic_acceleration` from `V_bar`; their difference is **`acceleration_gap`**. The engine estimates an MHD-style **magnetic pressure** acceleration `magnetic_pressure_acceleration` (not single-particle `q v × B`). See [`docs/theory/EG_MODEL_THEORY_EN.md`](docs/theory/EG_MODEL_THEORY_EN.md) §2 (the MHD bridge) and §3 (mathematical implementation; PT mirror: [`EG_MODEL_THEORY_PT.md`](docs/theory/EG_MODEL_THEORY_PT.md)).

Core formulas in code:

- `u_B = B² / (2 μ₀)`
- `Δu_B = u_B(neighbors) − u_B(center)` on the sky map
- `magnetic_pressure_gradient ≈ Δu_B / Δx` with `Δx` from the scenario grid (kpc → m)
- `magnetic_pressure_acceleration = −(magnetic_pressure_gradient) / ρ_ion_eff`

**Scale-mixing caveat:** HEALPix neighbors are separated by ~arcminutes on the sky; `scenario_gradient_scale_kpc` sets a disk spatial step. This is documented as exploratory.

### CMB anchoring (dimensionful scale)

Stokes amplitudes map to Tesla via `√(2 μ₀ u_ref)` with `u_ref = ENERGY_DENSITY_CMB` in `src/physics/engine.py`. This is a **calibration anchor**, not a claim that disk fields originate in the CMB.

## Morphology and LSB proxy

Hubble stage `T` from the JSON maps to `morphology_class`. Photometric LSB uses `SBdisk_solLum_pc2 < 100` (solar luminosities per pc²) for evidence filtering.

## Grid search

Nine scenarios: `scenario_gradient_scale_kpc` ∈ {0.05, 0.1, 0.5} kpc, `scenario_ionization_fraction` ∈ {0.001, 0.01, 0.1}. Among **reliable** scenarios, the best minimizes `|log10(theoretical_gap_coverage_ratio)|` toward unity.

## Outputs (under `results/`)

- `step_validation/` — population tables after the grid: `validation_all.csv` (all grid entrants), `validation_reliable.csv` / `validation_not_reliable.csv` (same schema, partitioned on `validation_is_reliable`), `validation_scientific_evidence_report.csv` (**supreme** cohort: best grid **plus** merged Discovery columns; headline ratio is always CMB-anchored), `validation_anchor_sensitivity_sweep.csv` (CMB vs turbulent-equipartition sensitivity), `validation_consistency_audit.md` (bias audit; re-reads `galaxy_analysis` JSON).
- `grid_results/` — one CSV per (scale, ion) scenario (galaxies that passed the Discovery pre-filter only)
- `step_discovery/` — `discovery_all.csv` / `discovery_reliable.csv` / `discovery_not_reliable.csv` / `discovery_cohort_rho_band.csv` — four per-galaxy tables from **Phase 1** (always written). All carry the headline column `required_density_for_match` (kg/m³, median over underperforming shells), the SI `delta_r_m`, the deprecated duplicate column `median_required_density_for_match_kg_m3` (same value as `required_density_for_match`), and `convergence_ratio_to_wim`. `discovery_reliable` is the `discovery_is_reliable == True` subset; `discovery_not_reliable` is its complement; `discovery_cohort_rho_band` applies the dynamic ±30% band around the reliable-subset median of `required_density_for_match` (internal-consistency list; see `ResultPersister.DISCOVERY_ELITE_BAND_*`). The `step_` prefixes keep phase folders adjacent when sorted.
- `graphs/population/` — cohort diagnostic PNGs (Discovery + validation CSVs), written automatically at the end of a full pipeline run.
- `galaxy_analysis/{`validation_scientific_evidence`|`validation_reliable`|`validation_not_reliable`}/<NAME>/` — per-galaxy CSV, JSON summary, and (for reliable / evidence tiers) publication PNGs (force, RAR, signed accelerations) after the run. **Phase 2 — Validation** writes per-galaxy CSV and JSON under this tree once the best grid scenario is known; **Phase 3 — Graphs** runs `publication_bundle` for cohort PNGs under `results/graphs/population/` and per-galaxy publication figures where the pipeline tier allows PNG export. Galaxies with `validation_is_reliable == True` start under `validation_reliable/`; those that appear in `validation_scientific_evidence_report.csv` are then moved to `validation_scientific_evidence/` (same artefacts). All other grid entrants use `validation_not_reliable/` (CSV+JSON only in normal runs — no per-galaxy PNG pass).
- column semantics: this README (Outputs section) is the operational dictionary; the physics behind each column is in [`docs/theory/EG_MODEL_THEORY_EN.md`](docs/theory/EG_MODEL_THEORY_EN.md) §3 (PT: [`EG_MODEL_THEORY_PT.md`](docs/theory/EG_MODEL_THEORY_PT.md))

## Reliability audit (`audit_flags`)

Every per-galaxy summary and every scenario row carries a structured `audit_flags` list. Codes prefixed with `fail_` block ``validation_is_reliable`` (and, on scenario rows, mark ``is_reliable`` false); codes prefixed with `warn_` are **non-blocking** provenance diagnostics and may appear even when the galaxy is validation-reliable. Blocking codes spell out **exactly which gates failed**: `fail_lat` (foreground latitude), `fail_reliable_fraction` (too few shells survived the per-shell outlier gate), `fail_density_clamped` (gas-density floor dominates), `fail_stellar` (stellar gravity dominates), `fail_skipped_physics` (engine could not produce a_B), `fail_inclination` (galaxy inclination too low), `fail_calibration` (Stokes reference fell back to the emergency value). Non-blocking warnings: `warn_b_band` (inferred `|B|` outside the 0.1–100 µG plausibility band), `warn_no_pa` (no curated position angle — engine fell back to the axisymmetric gradient), `warn_pure_pressure` (irregular / BCD morphology — magnetic-tension toy collapses to zero), `warn_btfr_outlier` (galaxy lies more than 0.5 dex from the McGaugh+2012 BTFR ridge). Per-scenario CSV rows expose a derived `is_reliable` boolean; the galaxy-level inventory uses `validation_is_reliable` (see `validation_all.csv` in `step_validation/`); Discovery census rows use `discovery_is_reliable`. In all cases `audit_flags` is the source of truth.

Visualization for the latest run: per-galaxy PNGs sit under ``galaxy_analysis/`` next to each galaxy's CSV/JSON; cohort panels live under ``results/graphs/population/`` (built from the step_* CSVs).

## Tests

Pytest uses `tmp_path` for all generated files. Run:

```bash
pytest tests/ -q --cov=src --cov-fail-under=75
```

Planck-heavy tests require **healpy** in the environment.