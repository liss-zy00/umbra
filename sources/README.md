# `sources/` — external scientific inputs

This directory holds every external data product the pipeline consumes. It is
organised **by scientific provenance** rather than by file format: each
subdirectory groups together every artefact that originates from the same
public archive, so the lineage of every byte under `results/` is
traceable in one step.

```
sources/
  sparc/      SPARC rotation curves and Lelli+2016c morphology catalogue.
  morphology/ Curated disk PA (and related notes) keyed by SPARC name — multi-archive overlay.
  planck/     Planck Legacy Archive thermal-dust polarisation maps.
  ned/        Static configuration for NED / SIMBAD name resolution.
```

Each subdirectory carries its own `README.md` with the canonical citation, a
link to the public archive, the file inventory, and the Python module that
consumes it. None of these files are auto-generated — they are versioned
inputs and should only be replaced when the upstream archive is updated.

## Conventions

- The pipeline addresses every path relative to the repository root (see
  [`src/pipeline/constants.py`](../src/pipeline/constants.py)). Run scripts
  from the root so the relative paths resolve.
- Output artefacts go to `results/`, never under `sources/`.
- The `docs/personal/` directory (gitignored) is the only place for private
  notes about these inputs; do not commit annotations into the `sources/`
  tree itself.
