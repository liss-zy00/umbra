# Modelo E-G — Refinamentos V2 (melhorias propostas: física, dados, estatística)

> **Uso.** Notas pessoais em `docs/personal/refinements/to_do/`. Não cites este ficheiro em documentação pública. Teoria base: [`EG_MODEL_THEORY_PT.md`](../../theory/EG_MODEL_THEORY_PT.md).  
> **Conteúdo.** Primeiro um *snapshot* do **último** `python src/main.py --galaxies all` (ambiente conda `dark-matter`, com `healpy`). Depois uma lista de **melhorias** que fariam sentido a seguir — independentemente de qualquer outro ficheiro de roadmap na mesma pasta.

---

## 1. Leitura do último run completo

**Comando de referência:** `conda run --no-capture-output -n dark-matter python src/main.py --galaxies all`

### 1.1 Resumo que a consola imprime

- **175** alvos SPARC processados.  
- **158** passam os *gates* mínimos do Discovery.  
- **158** entram na grelha de validação (9 cenários por galáxia → **1 422** avaliações de matriz).  
- **28** galáxias com `validation_is_reliable == True` após a grelha.  
- **4** entradas no extracto de “evidência directa” (`validation_scientific_evidence_report.csv` — filtros Im/LSB + janela na razão de cobertura; **não** é a mesma coisa que as 28 fiáveis).

**Agregados detalhados:** `results/step_validation/validation_consistency_audit.md` (gerado de novo em cada run).

### 1.2 Números úteis do *audit* (sub-amostra fiável, N = 28)

- RAR: mediana de \(\|\mathrm{residual}\|\) ≈ **0,162** dex (MAD σ ≈ 0,109 dex).  
- BTFR: N = 16 com residual finito; **4** galáxias acima de 0,5 dex ⇒ **25 %** *outliers* (mediana \|residual\| ≈ 0,32 dex).  
- *Foreground* (latitude): factor médio baixa / alta latitude ≈ **0,05×** (ver nota no *audit*).  
- Pearson `galactic_latitude_degrees` *vs.* `theoretical_gap_coverage_ratio`: **−0,19** (N = 28); em **N = 158** resumos carregados: **−0,26**.

### 1.3 *Flags* mais frequentes (N = 158 resumos)

| Código | Galáxias | % |
|--------|----------|---:|
| `warn_no_pa` | 130 | 82,3 |
| `fail_reliable_fraction` | 128 | 81,0 |
| `fail_stellar` | 79 | 50,0 |
| `warn_pure_pressure` | 40 | 25,3 |
| `warn_btfr_outlier` | 27 | 17,1 |

### 1.4 Booleans que não se devem confundir

- **`discovery_is_reliable`** — censo Discovery com \(f_{\mathrm{ion}}=1\) e escala física Planck; alimenta o pré-filtro da grelha.  
- **`validation_is_reliable`** — pós-grelha; é o conjunto das **28** acima.  
- **`is_reliable` em `scenario_*.csv`** — booleano **por** cenário \((\Delta x, f_{\mathrm{ion}})\).

### 1.5 Discovery fiável — medianas de feixe e densidade pedida

Valores lidos de `results/step_discovery/discovery_reliable.csv` **após** o mesmo run:

| Métrica | Valor |
|---------|--------|
| Mediana `required_density_for_match` | **2,46 × 10⁻²⁵ kg m⁻³** |
| Mediana `beam_size_kpc` | **22,4 kpc** |
| Mediana `resolution_correction_factor` (*k* = \(D_{\mathrm{beam}}/L_{\mathrm{grad}}\) no Discovery) | **2,29** (min ≈ 1,90 · max ≈ 2,40) |
| Linhas na coorte *elite* dinâmica (`discovery_cohort_rho_band.csv`) | **13** |

### 1.6 *Beam smearing* — teste log–log (OLS)

Script auxiliar: `src.support.beam_smearing_plot.plot_beam_smearing` → `results/graphs/beam_smearing_vs_distance.png`.

| Subconjunto | N | declive log₁₀–log₁₀ | Pearson r |
|-------------|---|----------------------|-------------|
| `discovery_reliable.csv` | **158** | **−0,68** | **−0,23** |
| `validation_scientific_evidence_report.csv` | **4** | **+1,37** | **0,82** |

O segundo subconjunto é **pequeno** (filtros do extracto de evidência); interpretar com cautela.

**Leitura física (resumo).** O feixe Planck (~5′) dilui o gradiente reconstruído; o factor mediano **~2,3** explica parte do salto em relação ao *floor* WIM **1 × 10⁻²⁷ kg m⁻³**, mas **não** explica sozinho **todas** as ordens de magnitude até **~10⁻²⁵ kg m⁻³**. Entram calibração da energia de referência (ver *sweep* CMB *vs.* equipartição, `ratio_log10_delta` ≈ **−2,92 dex** nos controlos do `validation_anchor_sensitivity_sweep.csv`), geometria/tensão, e densidades ISM reais.

**Próximo passo forte:** polarização rádio resolvida (JVLA / SKA) num punhado de discos — §3.1.

### 1.7 *Smoke* local

Para regressão rápida: `--galaxies default` (10 alvos). As contagens **não** substituem o run completo para estatística populacional.

---

## 2. Melhorias: MHD “2ª ordem” no diagnóstico

O motor já endureceu o gradiente de pressão magnética. Ainda fora do escopo:

1. **Viscosidade de Spitzer** — estimar tempos `τ_visc(r)` por casca para marcar onde o MHD ideal é ingénuo.  
2. **Condução térmica anisotrópica** — coluna de tempo paralelo *vs.* `plasma_beta_thermal`.  
3. **Pressão de raios cósmicos** — terceira pressão além de térmica + magnética; *proxy* de SFR; coluna `plasma_beta_with_cr`.

Cada uma pode entrar como *diagnostic* novo sem rebentar a *spine* `provider → engine → exporter`.

---

## 3. Melhorias: dados externos

1. **JVLA / MeerKAT / SKA-mid** — mapas polarizados resolvidos (~100 pc); *provider* alternativo ao Planck com a mesma interface de amostragem.  
2. **Gaia DR4** — distâncias e cinemática mais apertadas para o volume local; regenerar overlay de PA/inclinação.  
3. **Síncrotron 1,4 GHz** (WSRT-SINGS, VLASS, …) — amarrar amplitude de \(B\) ao disco com fluxo observado, reduzindo dependência de uma âncora global.

---

## 4. Melhorias: estatística, grelha e comunicação

- **Multiplicidade da grelha 3×3** — documentar critério (*BIC*, *hold-out* radial, ou par fixo *a priori*); publicar sempre métricas num **cenário físico default** antes da grelha exploratória.  
- **Gates de forma** — ex.: exigir correlação de Pearson > 0 para classe “evidência forte” quando fizer sentido.  
- **Estratificar por `warn_pure_pressure`** — tabelas separadas Im/BCD (*só pressão*) *vs.* espirais com ramo de tensão activo, para não misturar ordens de grandeza.  
- **RAR *vs.* narrativa MHD** — deixar explícito que RAR pequena é cheque lateral, não prova do canal `a_B`.  
- **BTFR** — com N baixo, intervalos binomiais na fração de *outliers*.  
- **Pacote *paper-ready*** — publicar `validation_anchor_sensitivity_sweep.csv` junto das tabelas principais e uma figura âncora CMB *vs.* equipartição por controlo (DDO154, NGC3741, …).

---

## 5. O que já existe *vs.* o que este ficheiro pede

| Área | Já no `src/` / `results/` | Próximo passo (este doc) |
|------|---------------------------|--------------------------|
| Diluição de feixe | `resolution.py`, colunas `beam_size_kpc`, `resolution_correction_factor*` | JVLA *vs.* Planck (§3.1) |
| OLS *beam* | Texto §1.6 + PNG regenerável | Re-gerar figura após cada run completo |
| MHD 2ª ordem | — | §2 |
| Novos *providers* | Planck + SPARC + NED + … | §3 |
| Rigor / grelha | *Audit* + `thresholds.py` | §4 |

---

*Roadmaps entregues: `refinements/done/`. Pendências UX + *viz* avançado: [`Model_Refinements_V4.md`](../Model_Refinements_V4.md).*
