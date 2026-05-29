# Modelo E-G — Refinamentos V1 (já implementado no código)

> **Uso.** Notas pessoais em `docs/personal/refinements/done/`. Não pões links para este ficheiro em README, papers ou guias públicos.  
> **Autoridade física.** Se algo aqui discordar da teoria base, prevalece [`EG_MODEL_THEORY_PT.md`](../../theory/EG_MODEL_THEORY_PT.md) (ou EN equivalente).

Isto é o inventário **do que o motor e o pipeline já fazem hoje**, em **seis** blocos claros.

| # | O quê | Onde ver (código / CSV) |
|---|--------|-------------------------|
| **1** | Gradiente de pressão magnética **no plano do disco**: ajuste planar (LSQ) ao anel Planck; se faltar *position angle*, `warn_no_pa` e fallback axissimétrico com `cos i`. | `PhysicsEngine`, gradiente magnético; `sources/morphology/disk_position_angle_overlay.json` |
| **2** | **Incerteza** em distância, inclinação e `V_flat`: bandas `p16` / `p84` na razão de cobertura e no `convergence_ratio_to_wim` (Monte Carlo fechado, barato). | `validation_scientific_evidence_report.csv`; `UncertaintyConstants` em `src/config/thresholds.py` |
| **3** | **Sensibilidade da âncora** Stokes→Tesla: CMB vs *floor* de equipartição turbulenta (dex entre escalas). | `results/step_validation/validation_anchor_sensitivity_sweep.csv` |
| **4** | **Tensão magnética** (*toy* espiral log) além do canal só-pressão: acelerações paralelas, `theoretical_gap_coverage_ratio_full_lorentz`; Im/BCD com `warn_pure_pressure` (só pressão). | `src/physics/diagnostics/tension.py` |
| **5** | **RAR** como cheque lateral: mediana de \|resíduo\| em dex por galáxia; agregação no *audit* populacional. | Coluna `rar_median_abs_residual_dex`; `validation_consistency_audit.md` §7a |
| **6** | **BTFR** + **resolução / feixe**: *ridge* McGaugh+2012, `btfr_residual_dex`, `warn_btfr_outlier`; módulo `resolution.py` (`beam_size_kpc`, `resolution_correction_factor*`), coorte *elite* dinâmica em Discovery. | `src/physics/diagnostics/btfr.py`, `resolution.py`; `discovery_*.csv` |

## Notas rápidas (decisões que importam)

- PA extra vem do overlay curado; o resto dispara `warn_no_pa` (não bloqueia sozinho).
- MC usa `N` e *seed* fixos em `thresholds.py` — reprodutibilidade vem dos catálogos, não de ruído em runtime.
- **Pitch** espiral do *toy* de tensão está em constantes de configuração (valor conservador).
- BTFR usa correção de hélio e massa estelar como no código (`btfr.py`).
- O *sweep* de âncora turbulenta não polui o run principal: vive no CSV dedicado.

## Onde estão os resultados

- `results/step_discovery/` — censo Discovery (todos / fiáveis / não fiáveis / banda elite).  
- `results/step_validation/` — grelha, evidência, *sweep*, `validation_consistency_audit.md`.  
- `results/galaxy_analysis/` — por galáxia (CSV, JSON, figuras).  
- O `README.md` na raiz do repo é **gerado** no fim do pipeline (`readme_builder.py`).
