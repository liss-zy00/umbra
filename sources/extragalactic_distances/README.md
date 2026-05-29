# Extragalactic Distance Database (EDD) artefacts

This folder holds static CSV extracts from the [Extragalactic Distance Database](https://edd.ifa.hawaii.edu/) (University of Hawaii / Institute for Astronomy), used only for **luminosity-distance hierarchy** joins in the pipeline (no physics beyond distance modulus → Mpc conversion documented in `src.config.thresholds`).

## Citation

EDD is described in Tully, R. B. et al., 2009, *AJ*, **138**, 323 (NASA ADS bibcode **2009AJ....138..323T**). Please cite that reference when publishing results that rely on these tables.

## Files

| File | Role in pipeline |
| --- | --- |
| `edd_cmd_trgb.csv` | Tier **1**: TRGB/CMD distance moduli (`DM_tip`, asymmetric `eDM_lo` / `eDM_hi` → conservative symmetric σ in modulus space). |
| `edd_cf4_individual_distances.csv` | Tier **2**: CF4 individual moduli — prefer `DMcep`/`eDMcep`, else `DMtrgb`/`eDMt`, else aggregate `DM`/`eDM`. |
| `edd_cosmicflows3_distances.csv` | Tier **3**: Cosmicflows-3 — use `DM`/`eDM` when finite, else `Dist` (Mpc) with σ propagated from `eDM` when a modulus is recoverable. |
| `sparc_name_to_pgc.json` | Optional **overrides** for `normalize_galaxy_name` → integer `pgc`. If `positions` is empty, `pgc` is inferred from the Cosmicflows-3 `Name` column at bootstrap (first match wins). |

## Hierarchy

For each galaxy the resolver loads `pgc` from the crosswalk, then applies **TRGB → CF4 → Cosmicflows-3**; if no tier yields a finite distance, the pipeline keeps the **SPARC rotmod header** distance (tier 4 fallback). When tiers 1–3 supply the adopted distance, `prov_high_precision_dist` is set on the galaxy metadata and surfaced in validation exports.

## Asymmetric TRGB errors

When both `eDM_lo` and `eDM_hi` are finite, the modulus uncertainty is taken as `max(eDM_lo, eDM_hi)` so the propagated luminosity-distance band is conservative.
