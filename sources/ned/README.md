# `sources/ned/` — NED / SIMBAD name resolution config

**Origin.** Static configuration for the
[NASA/IPAC Extragalactic Database](https://ned.ipac.caltech.edu/) (NED)
and [SIMBAD](https://simbad.cds.unistra.fr/simbad/) astronomical databases,
which the pipeline queries at runtime through `astroquery` to obtain
equatorial coordinates and redshift for each SPARC galaxy.

NED's bare-name interpreter rejects some SPARC labels (notably the Schombert
LSB F-series like `F561-1`). The two JSON files in this directory teach the
provider how to resolve those labels deterministically without needing a
network probe per attempt.

## Files

| File | Content |
| --- | --- |
| `ned_alias_map.json` | Three-section config used by [`NedProvider`](../../src/providers/ned.py): `alias_map` (exact SPARC → NED canonical name), `pattern_rules` (regex → templated query, e.g. `^F\d{3}-.+$` → `LSB {name}`), `ned_prefix_attempts` (ordered list of prefixes to try when the bare name fails). The full schema is documented in the file's `_comment` field. |
| `ned_icrs_cone_seeds.json` | Optional `galaxy_name -> {ra_deg, dec_deg}` ICRS seeds for the NED cone search. Empty by default; populate only when both NED and SIMBAD fail to resolve a SPARC label by name. Path constant: `NED_ICRS_CONE_SEEDS_JSON`. |

## Consumers

- [`src/providers/ned.py`](../../src/providers/ned.py) — loads `ned_alias_map.json` at construction time; the runner passes seeds from `ned_icrs_cone_seeds.json` as optional `icrs_deg` when resolving names.

## Citation

> Helou G., Madore B. F., Schmitz M., Bicay M. D., Wu X., Bennett J., *The
> NASA/IPAC Extragalactic Database*, in *Databases & On-line Data in
> Astronomy*, eds. Albrecht & Egret (1991).
>
> Wenger M. et al., *The SIMBAD astronomical database*, A&AS 143, 9 (2000).
