# `sources/sparc/` — SPARC rotation curves and morphology

**Origin.** SPARC ("Spitzer Photometry and Accurate Rotation Curves"),
Lelli, McGaugh & Schombert 2016, *AJ* **152**, 157. Public release archived at
the SPARC project page (<https://astroweb.case.edu/SPARC/>).

The SPARC sample provides 175 nearby disk galaxies with Spitzer 3.6 µm
photometric mass models and high-quality HI/Hα rotation curves. The kinematic
columns (`V_gas`, `V_disk`, `V_bulge`, `V_obs`) are the input to the
baryonic gap calculation that drives the E-G model test.

## Files

| File | Content |
| --- | --- |
| `rotation_curves/` | One `*_rotmod.dat` per galaxy (Lelli+2016c SPARC Newtonian mass-model tables: radius, velocities, surface densities). Header carries the SPARC distance. |
| `SPARC_Lelli2016c.json` | Parsed morphology/mass catalogue (Hubble type, inclination, V_flat, HI mass, photometric LSB flag) used by [`SparcMorphologyProvider`](../../src/providers/sparc.py). |
| `SPARC_Lelli2016c.mrt` | Raw machine-readable table from the SPARC release. Source of truth — the JSON is derived from this file. |

Disk position angles are **not** part of the SPARC archive; the curated overlay lives under
[`sources/morphology/`](../morphology/README.md) (see ``MORPHOLOGY_DISK_PA_OVERLAY_JSON`` in
[`constants.py`](../../src/pipeline/constants.py)).

## Consumers

- [`src/providers/sparc.py`](../../src/providers/sparc.py) — `SparcProvider` parses the rotmod files; `SparcMorphologyProvider` parses the JSON.
- [`src/pipeline/constants.py`](../../src/pipeline/constants.py) — pins `SPARC_DIR`, `SPARC_CATALOG_JSON`, `MORPHOLOGY_DISK_PA_OVERLAY_JSON`, and `NED_ICRS_CONE_SEEDS_JSON`.

## Citation

> Lelli F., McGaugh S. S., Schombert J. M., *SPARC: Mass Models for 175
> Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves*,
> AJ 152, 157 (2016).
