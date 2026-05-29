# `sources/morphology/` — disk geometry overlay (multi-archive, SPARC-keyed)

**Purpose.** SPARC does not publish disk major-axis position angles. The pipeline still needs a PA
to project the Planck polarisation gradient into the disk plane for the spiral-arm channel.

**Why one JSON instead of `sources/hyperleda/` vs `sources/ned/` splits.** Each galaxy row carries
its own `source` string (HyperLEDA, NED, …). Splitting into two provider folders would duplicate
keys and require merge precedence rules in code; a single overlay file keeps one join key
(SPARC name) and a clear per-row provenance field.

| File | Role |
| --- | --- |
| `disk_position_angle_overlay.json` | Curated `position_angle_deg` (east of north), optional uncertainty, plus optional audit fields (e.g. NED kinematic rows). Consumed by [`SparcMorphologyProvider`](../../src/providers/sparc.py) using [`MORPHOLOGY_DISK_PA_OVERLAY_JSON`](../../src/pipeline/constants.py). |

## Citation policy (literature-only PA)

**Only** add or change `position_angle_deg` when a **verifiable primary reference** exists (published bibcode / DOI / stable catalogue table or figure at object level). Each row must carry a `source` string that a reviewer can trace (table/column, bibcode, or stable URL). **Do not** fill gaps with uncited guesses (for example legacy PDFs that do not list PA for that object). Galaxies without a citable PA keep the axisymmetric gradient path and the engine’s `warn_no_pa` provenance flag.

Offline refresh of NED-hosted rows: [`src/support/ned_pa_supplement.py`](../../src/support/ned_pa_supplement.py) (`refresh_disk_pa_overlay_from_validation_csv`, import-only). Any future manual rows written into the JSON must satisfy the same citation rule as above.

## Consumers

- [`src/providers/sparc.py`](../../src/providers/sparc.py) — merges overlay PA into each morphology lookup entry.
