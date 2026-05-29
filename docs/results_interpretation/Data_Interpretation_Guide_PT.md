# Guia de interpretação de dados (pipeline E-G)

## Objectivo

Este documento complementa [`Scientific_Evidence_Narrative_PT.md`](Scientific_Evidence_Narrative_PT.md) com **blocos breves em formato pergunta–resposta** para as colunas principais dos outputs Discovery e de `validation_scientific_evidence_report.csv`. **Não** duplica a semântica completa de `audit_flags`; para isso, ver §5 da narrativa e a metodologia §4.4 em [`Draft_Paper_EG_Model_PT.md`](../academic/Draft_Paper_EG_Model_PT.md).

**Coortes.** `validation_scientific_evidence_report.csv` contém apenas *supreme rows* (fiável + inclinação) que passam o filtro morfológico Im/LSB fotométrico **e** a janela de razão de cobertura 0,2–5,0 (`build_direct_evidence_rows`). Para toda a coorte `validation_is_reliable`, use `validation_reliable.csv` + `validation_consistency_audit.md` (ver §1 da narrativa).

---

## `required_density_for_match`

- **Pergunta.** Que densidade de plasma ionizado o canal de pressão magnética precisaria para fechar a lacuna cinemática casca a casca?
- **Significado físico.** Razão casca-a-casca entre a escala do gradiente de energia magnética e a magnitude do *gap*, expressa como uma densidade \(\rho_{\mathrm{req}}\) que tornaria o ansatz de gradiente de pressão magnética MHD ideal auto-consistente nesse raio.
- **Onde encontrar.** Colunas em `discovery_*.csv`; fundidas em `validation_scientific_evidence_report.csv` apenas para linhas incluídas por `build_direct_evidence_rows` (Im ou LSB fotométrico + razão de cobertura na janela configurada).
- **O que é um valor saudável.** Valores **perto** de \(\sim 10^{-27}\,\mathrm{kg\,m^{-3}}\) *depois* de cruzar feixe, `resolution_correction_factor_*` e o *sweep* de âncora alinham com a narrativa de acoplamento ao meio ionizado difuso (WIM). **Medianas populacionais muito acima** dessa âncora são o estado empírico habitual em varreduras SPARC+Planck documentadas — leia como tensão entre modelo e dados, não como falha silenciosa do formalismo (ver [`Scientific_Evidence_Narrative_PT.md`](Scientific_Evidence_Narrative_PT.md) §5.3).

---

## `convergence_ratio_to_wim`

- **Pergunta.** Quão perto está o requisito inferido da densidade de referência empírica do WIM?
- **Significado físico.** Razão \(\rho_{\mathrm{req}} / \rho_{\mathrm{WIM}}\) (ou estatística equivalente na linha Discovery).
- **Onde encontrar.** CSVs Discovery; coluna fundida em `validation_scientific_evidence_report.csv`.
- **O que é um valor saudável.** A unidade é correspondência perfeita à âncora; a lista *elite* em `discovery_cohort_rho_band.csv` define-se sobretudo por uma banda dinâmica em `required_density_for_match` (ver `discovery_is_elite` abaixo), não apenas por esta razão.

---

## `theoretical_gap_coverage_ratio`

- **Pergunta.** Quanto do déficit gravitacional o proxy de pressão magnética fechou no disco externo, no melhor cenário do grid?
- **Significado físico.** Razão de módulos \(|\langle a_B\rangle_{\mathrm{edge}}| / |\langle \mathrm{gap}\rangle_{\mathrm{edge}}|\) após máscaras de casca e selecção de cenário.
- **Onde encontrar.** Métricas por cenário e a melhor linha em `validation_scientific_evidence_report.csv`.
- **O que é um valor saudável.** Consistência em ordem de grandeza com a unidade (por exemplo dentro de \([0.1, 10]\) em espaço logarítmico) mantém o ansatz falsificável mas não imediatamente excluído nesse cenário; valores muito fora dessa banda indicam desvio ao par \((\Delta x, f_{\mathrm{ion}})\) escolhido.

---

## `pearson_correlation`

- **Pergunta.** Concordam as **formas** dos perfis radiais, não só as magnitudes na borda?
- **Significado físico.** Correlação de Pearson entre `magnetic_pressure_acceleration(r)` e `acceleration_gap(r)` sobre cascas confiáveis no cenário em teste.
- **Onde encontrar.** Tabelas de cenário e `validation_scientific_evidence_report.csv` para o melhor cenário.
- **O que é um valor saudável.** Correlação positiva indica que os perfis sobem e descem em conjunto; valores perto de zero ou fortemente negativos enfraquecem uma interpretação simples de condutor radial partilhado mesmo que uma razão num único raio pareça aceitável.

---

## `delta_r_m` e `physical_pixel_scale_kpc`

- **Pergunta.** Com que passo físico é estimado o gradiente de pressão?
- **Significado físico.** `delta_r_m` é a separação entre vizinhos em SI usada na diferença finita de \(\nabla u_B\); `physical_pixel_scale_kpc` reporta a mesma escala em kpc para comparação rápida com escalas do disco e com os botões de gradiente do cenário.
- **Onde encontrar.** Linhas Discovery; ambas aparecem em `validation_scientific_evidence_report.csv` quando existem dados Discovery para essa galáxia.
- **O que é um valor saudável.** No Discovery, \(\Delta r\) deve coincidir com a escala física do pixel Planck quando `Δx_phys` está disponível; discrepâncias grandes entre a escala física do pixel e o \(\Delta x\) só de cenário deslocam `ρ_req` e devem ser lidas em conjunto com `audit_flags`.

---

## `magnetic_field_microgauss`

- **Pergunta.** O campo inferido está numa banda de plausibilidade da literatura para a classe alvo?
- **Significado físico.** Ordem de grandeza do campo implícita pelo caminho de calibração Stokes que alimenta \(u_B\).
- **Onde encontrar.** Resumos por galáxia e tabelas Discovery conforme exposto pelos exportadores do pipeline.
- **O que é um valor saudável.** Valores em \(\mu\mathrm{G}\) totalmente fora do típico para discos exigem escrutínio conjunto com auditorias de calibração (por exemplo `warn_b_band` / *fallbacks* documentados na narrativa).

---

## `plasma_beta_thermal_median`

- **Pergunta.** O ansatz de pressão magnética é localmente apropriado (pressão térmica vs magnética)?
- **Significado físico.** Mediana do plasma \(\beta\) de pressão gasosa térmica sobre pressão magnética nas cascas contribuintes — sonda MHD local, não substituto de tensão completa e transporte anisotrópico.
- **Onde encontrar.** Métricas ao nível da galáxia nos outputs do pipeline quando calculado.
- **O que é um valor saudável.** \(\beta\) muito elevado em todo o lado sugere pressão magnética dinamicamente negligenciável nesses raios; \(\beta\) muito baixo em todo o lado pode indicar \(B\) excessivo ou efeitos de *floor* de densidade — ler junto com `audit_flags`.

---

## `discovery_is_reliable`, `validation_is_reliable` e `audit_flags`

- **Pergunta.** Qual *flag* de confiabilidade corresponde a qual fase do pipeline?
- **Significado físico.** `discovery_is_reliable` é o *flag* da Fase 1 (censura Discovery): *gates* mínimos (\(\rho_{\mathrm{req}}\) finita, auditorias de latitude / Planck / NED). `validation_is_reliable` é o resultado ao nível da galáxia após o grid 9-cenário e os *gates* completos de casca / estrelar / inclinação / calibração (`validation_all.csv` e `*_summary.json` em `data_reliability_assessment`). As linhas do grid por cenário continuam a expor um booleano `is_reliable` local derivado do mesmo vocabulário `audit_flags` na célula \((\Delta x, f_{\mathrm{ion}})\). `audit_flags` permanece a fonte de verdade estruturada em todos os casos.
- **Onde encontrar.** `results/step_discovery/discovery_all.csv` lista todas as galáxias processadas; `discovery_reliable.csv` / `discovery_not_reliable.csv` particionam por `discovery_is_reliable`. Após o grid, `results/step_validation/validation_all.csv` lista todas as galáxias que entraram no grid com `validation_is_reliable` em cada linha (mesmo *schema* de antes, incluindo `audit_flags`, `reliability_failure_codes` e `note` em texto separados por vírgulas onde aplicável; filtre em ferramentas em vez de CSVs separados de passou/não passou). Tabelas de grid: `results/grid_results/scenario_*.csv`.
- **O que é um valor saudável.** Uma galáxia pode estar Discovery-*passed* (entra no grid) e falhar na validação após máscaras mais estritas; tratar estatísticas populacionais em coortes Discovery vs Validação separadamente.

---

## `discovery_is_elite` (coluna do reporte)

- **Pergunta.** Como se define a lista curta de “consistência interna” após o Discovery?
- **Significado físico.** Quando a banda dinâmica está disponível, `discovery_is_elite` é verdadeiro para linhas Discovery-*reliable* cujo `required_density_for_match` cai em `[0.7, 1.3]` vezes a **mediana** de `required_density_for_match` no subconjunto *reliable* (ver `ResultPersister.DISCOVERY_ELITE_BAND_*`). Quando essa banda não pode ser formada, `discovery_is_elite` é falso e `discovery_elite_band_source` é `none`.
- **Onde encontrar.** `results/step_discovery/discovery_cohort_rho_band.csv`; espelhado como `discovery_is_elite` em `validation_scientific_evidence_report.csv` para galáxias que passam a banda dinâmica em \(\rho_{\mathrm{req}}\).
- **O que é um valor saudável.** Linhas *elite* destacam escalas \(\rho_{\mathrm{req}}\) coerentes dentro da população Discovery-*passed*; não substituem `validation_is_reliable`.

---

## `beam_size_kpc` e diluição por resolução (`resolution_correction_factor_*`)

- **Pergunta.** Como o feixe finito se relaciona com o passo do gradiente magnético?
- **Significado físico.** `beam_size_kpc` aproxima a escala física FWHM no céu convertida para kpc à distância adoptada. `resolution_correction_factor_discovery` guarda \(k \approx D_{\mathrm{beam}} / L_{\mathrm{grad}}\) com \(L_{\mathrm{grad}}\) da separação física do pixel Planck (Discovery). `resolution_correction_factor_grid` recalcula a mesma razão usando o botão `scenario_gradient_scale_kpc` do **melhor cenário**. A coluna legada `resolution_correction_factor` copia o valor de *grid* por compatibilidade retroativa.
- **Onde encontrar.** CSVs Discovery e colunas fundidas em `validation_scientific_evidence_report.csv`.
- **O que é um valor saudável.** \(k\) elevado avisa que o mapa Stokes suavizado pode diluir gradientes relativamente ao passo assumido do disco — ler `convergence_ratio_to_wim` e `required_density_for_match` em conjunto com estas colunas, não isoladamente.

---

## `theoretical_gap_coverage_ratio_p16` / `_p84`

- **Pergunta.** O que a incerteza catalogada em distância, inclinação e `V_flat` provoca na razão *headline*?
- **Significado físico.** Banda Monte Carlo sobre `theoretical_gap_coverage_ratio` (`N = 300` amostras, semente `42`). O propagador reescala a razão central analiticamente como `R = R_0 \cdot (D_0/D_s) \cdot (\sin i_s / \sin i_0)^2 \cdot (V_0/V_s)^2`. O mesmo reescalonamento multiplicativo é aplicado a `convergence_ratio_to_wim` para preencher `_p16` / `_p84` da coluna de convergência ao WIM.
- **Onde encontrar.** `validation_scientific_evidence_report.csv`.
- **O que é um valor saudável.** Banda simétrica centrada no ponto (`mediana ≈ ponto`) com largura controlada pela incerteza catalogada dominante (tipicamente a inclinação). Banda fortemente assimétrica sinaliza que uma entrada está perto de uma não linearidade (ex.: `sin i` muito baixo).

---

## `accel_pressure_only_m_s2`, `accel_full_lorentz_m_s2`, `theoretical_gap_coverage_ratio_full_lorentz`

- **Pergunta.** O termo de tensão magnética muda a conclusão?
- **Significado físico.** Medianas por galáxia lado a lado de (i) a aceleração só-de-pressão `a_B = -\nabla u_B / \rho_\mathrm{ion}` e (ii) a mesma somada do *toy* de tensão `a_\mathrm{tension,r}(r) = -\cos p \cdot B^2 / (\mu_0 \rho_\mathrm{ion} R_c(r))` sob geometria de espiral logarítmica com passo `p = 12°` e `R_c = r / \sin^2 p`. A razão de cobertura paralela é reportada como `theoretical_gap_coverage_ratio_full_lorentz`.
- **Onde encontrar.** `validation_scientific_evidence_report.csv`.
- **O que é um valor saudável.** Para hospedeiros espirais as duas razões devem diferir por no máximo ordem de unidade — o termo de tensão é um efeito de segunda ordem conhecido. Para hospedeiros irregulares / BCD o motor força a tensão a zero (`warn_pure_pressure`) e as duas colunas coincidem.

---

## `rar_median_abs_residual_dex`

- **Pergunta.** A cinemática é autoconsistente com a *Radial Acceleration Relation*?
- **Significado físico.** Mediana de `|log10(g_obs / g_pred)|` sobre cascas confiáveis contra a RAR de McGaugh+2016. Um valor `\lesssim 0,1` dex está dentro da dispersão intrínseca publicada; `\gtrsim 0,3` dex sinaliza uma galáxia cuja cinemática discorda do locus empírico da RAR.
- **Onde encontrar.** `validation_scientific_evidence_report.csv`; mediana populacional + MAD em `results/step_validation/validation_consistency_audit.md` (§7a).
- **O que é um valor saudável.** A maioria das galáxias SPARC tem `\le 0,13` dex.

---

## `btfr_residual_dex` e `warn_btfr_outlier`

- **Pergunta.** A galáxia está sobre o *ridge* Baryonic Tully–Fisher?
- **Significado físico.** Distância sinalizada (dex) da galáxia ao *ridge* de McGaugh+2012 `log10(M_bar / M_sun) = 4,0 \cdot log10(V_flat) + 1,99`, com `M_bar = 1,33 \cdot M_HI + \Upsilon_{3,6} \cdot L_{3,6}` e `\Upsilon_{3,6} = 0,5 M_sun/L_sun`. `warn_btfr_outlier` é verdadeira quando `|residual| > 0,5` dex.
- **Onde encontrar.** `validation_scientific_evidence_report.csv`; fração populacional de *outliers* em `results/step_validation/validation_consistency_audit.md` (§7b).
- **O que é um valor saudável.** `|residual| \lesssim 0,2` dex.

---

## `warn_no_pa` e `warn_pure_pressure`

- **Pergunta.** Em que pressupostos o motor recaiu para esta galáxia?
- **Significado físico.** `warn_no_pa` significa que o suplemento curado (`sources/morphology/disk_position_angle_overlay.json`) não lista um ângulo de posição para a galáxia; o motor recalcula `∇u_B` no referencial axissimétrico (plano do céu) com a correção `cos i`. `warn_pure_pressure` significa que a morfologia é Irregular / BCD: o *toy* de tensão magnética é forçado a zero, portanto `accel_full_lorentz_m_s2 = accel_pressure_only_m_s2`.
- **Onde encontrar.** `validation_scientific_evidence_report.csv`.
- **O que é um valor saudável.** Ambos os *flags* não bloqueiam — só carregam a informação de que o refinamento auxiliar não foi exercitado.

---

## `validation_anchor_sensitivity_sweep.csv`

- **Pergunta.** Quanto a razão *headline* depende da âncora de calibração `u_ref`?
- **Significado físico.** Duas razões por galáxia: `ratio_cmb` (a razão ancorada à CMB padrão) e `ratio_turbulent_equipartition` (a mesma razão com o calibrador reancorado num piso de pressão por equipartição turbulenta `u_turb = 0,5 \cdot \rho_gas \cdot \sigma_v^2`). `ratio_log10_delta = log10(ratio_turbulent / ratio_cmb)` captura a diferença em dex entre as duas âncoras.
- **Onde encontrar.** `results/step_validation/validation_anchor_sensitivity_sweep.csv`.
- **O que é um valor saudável.** `|ratio_log10_delta| \lesssim 0,5` dex — as duas âncoras devem concordar dentro de uma ordem de magnitude.

---

## Galeria de figuras (*publication bundle*)

- **Quando corre.** Depois de a **Fase 2 (Validação)** completar o *run* (grelha + artefactos por galáxia), o *runner* chama automaticamente :mod:`src.reporting.publication_bundle` (**Fase 3 — Gráficos**). O único ponto de entrada CLI continua a ser ``python src/main.py --galaxies …``; não existe subcomando separado para figuras.
- **Flags da CLI.** A selecção de galáxias usa ``--galaxies`` (ver `src/pipeline/cli.py`). O opcional ``--logs`` imprime proveniência verbose por galáxia e **desactiva** as barras tqdm para manter os logs legíveis.

- **Lista de galáxias.** Os PNGs por objecto são gerados para todas as galáxias que entraram na grelha de validação nesse *run* (o mesmo conjunto de ficheiros ``*_rotmod.dat`` que passaram o pré-filtro Discovery). Os painéis de população leem as tabelas ``results/step_*/*.csv`` recém-escritas em disco.

- **Layout de saída.** Os PNGs por galáxia (`*_acceleration_comparison.png`, `*_rar_plot.png`, `*_signed_accel_vs_radius.png`) são gerados apenas para galáxias listadas em `step_validation/validation_reliable.csv` e ficam ao lado do CSV/JSON em `results/galaxy_analysis/{validation_scientific_evidence|validation_reliable}/<NOME>/` após o *run* (a coorte de evidência científica é movida para `validation_scientific_evidence/`). Os painéis de população ficam em `results/graphs/population/`.

- **Política Markdown.** Por defeito **não** se escrevem `.md` ao lado destes PNG: eixos, títulos e *caption* no próprio gráfico transportam a física. A interpretação longa permanece nesta pasta (`Scientific_Evidence_Narrative_*.md`, este guia).

- **Tipografia e escalas.** O *job* aplica predefinições Matplotlib *serif* (DejaVu Serif, sem LaTeX externo). Os limites log usam quantis **nas galáxias processadas nesse run**.

- **Testes / pasta alternativa de coorte.** ``run_publication_bundle(..., cohort=…)`` grava os PNGs de população em ``results/graphs/population_<cohort>/`` nos testes unitários com um CSV de coorte específico.

---

## Onde olhar primeiro

- **`results/step_validation/validation_consistency_audit.md`** — resumo populacional e frequências de `audit_flags` (regenerado a cada *run*).
- **`discovery_cohort_rho_band.csv`** — subconjunto de consistência interna após a banda dinâmica em \(\rho_{\mathrm{req}}\).
- **`discovery_reliable.csv`** — coorte Discovery `discovery_is_reliable` para medianas de \(\rho_{\mathrm{req}}\) e estatísticas de razão.
- **`discovery_all.csv` / `discovery_not_reliable.csv`** — censo completo e complemento explícito.
- **`validation_scientific_evidence_report.csv`** — vista cruzada: campos de cabeçalho Discovery **mais** melhor cenário de cobertura do *gap* e correlação.
- **`results/graphs/population/`** — figuras de coorte do passo automático de publicação no fim de um *run* completo (ver Galeria acima).
- **Artefactos por galáxia** — `results/galaxy_analysis/{validation_scientific_evidence|validation_reliable|validation_not_reliable}/<NOME>/` contém CSV e JSON; PNGs de publicação acompanham as galáxias fiáveis (e a subpasta `validation_scientific_evidence` após reorganização do *run*).

---

## Leitura relacionada

- [`Scientific_Evidence_Narrative_PT.md`](Scientific_Evidence_Narrative_PT.md) — guia de leitura estável, semântica de gráficos e discussão de `audit_flags`.
- [`Scientific_Evidence_Narrative_EN.md`](Scientific_Evidence_Narrative_EN.md) — espelho em inglês.
