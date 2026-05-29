# Como Ler os Resultados do Pipeline E-G

> **Propósito.** Este documento é o **guia de leitura estável** para os
> artefatos produzidos pelo pipeline (o CSV de resumo por galáxia, o pacote
> de figuras em `results/graphs/population/` e PNGs por galáxia em `results/galaxy_analysis/`, e a lista `audit_flags` por galáxia). O módulo
> [`src/reporting/narrative.py`](../../src/reporting/narrative.py) pode gerar Markdown de resumo a partir de linhas fiáveis, mas **a CLI principal (`python src/main.py`) não grava** esse ficheiro em `results/` por defeito. **Este guia**
> explica *o que cada coluna significa*, *como ler cada gráfico* e *quais
> conclusões são autorizadas por cada combinação de valores*. Deliberadamente
> não congela números por galáxia de qualquer execução particular.
> **Versão em inglês:** [`Scientific_Evidence_Narrative_EN.md`](Scientific_Evidence_Narrative_EN.md).
> **Perguntas e respostas por coluna (campos principais):** [`Data_Interpretation_Guide_PT.md`](Data_Interpretation_Guide_PT.md) / [`Data_Interpretation_Guide_EN.md`](Data_Interpretation_Guide_EN.md).

---

## 1. Que ficheiros abrir primeiro (ordem de leitura)

O pipeline escreve muitos artefatos. Para **estatística da coorte de validação completa**
(`validation_is_reliable`), a ordem recomendada é:

1. **`results/step_validation/validation_consistency_audit.md`** — frequências de `audit_flags`, RAR/BTFR agregados, tabelas por morfologia e *checks* de viés (regenerado a cada *run*).
2. **`results/step_validation/validation_reliable.csv`** — uma linha por galáxia com `validation_is_reliable == True` (melhor cenário já colapsado no inventário).
3. **`results/step_validation/validation_all.csv`** — todas as galáxias que entraram na grelha (fiáveis ou não), para auditar rejeitadas.

**`results/step_validation/validation_scientific_evidence_report.csv`** não substitui o passo 2: é uma **sub-amostra de “evidência directa”** construída em [`build_direct_evidence_rows`](../../src/pipeline/aggregation/global_aggregator.py) a partir das *supreme rows* (galáxias com `validation_is_reliable` **e** *gate* de qualidade de inclinação no processador). Só entram linhas que cumulativamente:

- têm morfologia **Im (Dwarf Irregulars)** **ou** LSB fotométrico (limiar em `RELIABILITY_CUTS.PHOTOMETRIC_LSB_SB_THRESHOLD`, ver [`include_in_scientific_evidence_report`](../../src/pipeline/per_galaxy/morphology.py));
- têm `theoretical_gap_coverage_ratio` finita na janela **0,2–5,0** (`RELIABILITY_CUTS.DIRECT_EVIDENCE_RATIO_MIN/MAX` em [`thresholds.py`](../../src/config/thresholds.py)).

Por isso o CSV de evidência tem **muito menos linhas** que `validation_reliable.csv`. Serve para síntese Im/LSB e figuras de população que leem esse ficheiro; um revisor que queira a coorte fiável **inteira** deve começar por `validation_reliable.csv` + *audit*.

- **`results/graphs/population/`** — PNGs de diagnóstico de coorte ao fim de um *run* de validação bem-sucedido (tabelas `step_*` + metadados em `galaxy_analysis`).

## 2. Colunas de `validation_scientific_evidence_report.csv`

Cada linha deste CSV é produzida a partir de uma *supreme row* (ver §1), enriquecida em
[`build_direct_evidence_rows`](../../src/pipeline/aggregation/global_aggregator.py) e escrita
por [`save_scientific_evidence_report`](../../src/reporting/persistence/exporter.py). O *schema*
completo (incl. `validation_status`, `wim_density_match_score`, diagnósticos de feixe) está no
código do exportador — a lista abaixo cobre os campos com maior carga física.

As colunas descritas aqui são as que carregam significado físico para o
teste. São produzidas por
[`src/pipeline/per_galaxy/metrics_builder.py`](../../src/pipeline/per_galaxy/metrics_builder.py)
e persistidas por
[`src/reporting/persistence/exporter.py`](../../src/reporting/persistence/exporter.py).

A tabela inclui também **campos de cabeçalho do Discovery** fundidos a partir da Fase 1
(`required_density_for_match`, `convergence_ratio_to_wim`, `delta_r_m`,
`physical_pixel_scale_kpc`, `discovery_is_reliable`, `discovery_is_elite`);
ver [`Data_Interpretation_Guide_PT.md`](Data_Interpretation_Guide_PT.md) para
perguntas e respostas compactas sobre essas colunas.

### 2.1 Identificação

- `galaxy_name` — nome no catálogo SPARC.
- `morphology_class` — grupo de classe Hubble-stage a partir do *lookup*
  JSON do Lelli+2016. Leia junto com a razão de cobertura: anãs e LSBs
  são a amostra onde a hipótese de pressão magnética *deve importar mais*,
  porque a contribuição bariônica gravitacional ali é pequena.

### 2.2 Parâmetros do melhor cenário

- `scenario_gradient_scale_kpc` — o valor de Δx (em kpc) usado para
  converter a diferença entre vizinhos HEALPix `Δu_B` em um gradiente.
  Três valores são explorados: 0,050, 0,100 e 0,500.
- `scenario_ionization_fraction` — a fração de ionização do WIM `f_ion`
  multiplicando `ρ_gas` para obter a densidade responsiva `ρ_ion`. Três
  valores: 0,001, 0,010 e 0,100.

Um cenário em que o melhor ajuste prefere `scenario_gradient_scale_kpc`
pequeno (e.g. 50 pc) e `scenario_ionization_fraction` pequena (e.g. 0,1 %)
é fisicamente consistente com uma **geometria magnética filamentar acoplada
ao WIM rarefeito**. Um melhor ajuste em `scenario_gradient_scale_kpc`
grande (e.g. 500 pc) e `f_ion` grande (e.g. 10 %) é consistente com um
**campo de disco distribuído suavemente atuando sobre uma camada ionizada
mais densa**. Ambos são, em princípio, permitidos, mas apontam para
regimes físicos diferentes que refinamentos subsequentes (sensibilidade
de âncora, tensão magnética) precisam desambiguar.

### 2.3 Os dois números principais

- `theoretical_gap_coverage_ratio` —
  `|⟨a_B⟩_edge| / |⟨gap⟩_edge|`. A grandeza natural a ser observada é seu
  logaritmo na base 10. Valores dentro de `[10^{-1}, 10^{+1}]` (isto é,
  a uma ordem de magnitude da unidade) significam que o *proxy* de pressão
  magnética está na ordem de magnitude correta para fechar o déficit
  cinemático no disco externo; valores fora dessa banda falsificam o
  *ansatz* no cenário sob o qual foi calculado.
- `pearson_correlation` — coeficiente de correlação de Pearson entre os
  perfis radiais de `magnetic_pressure_acceleration` e
  `acceleration_gap` sobre as cascas confiáveis. Um valor positivo indica
  que os dois perfis *crescem juntos*: não apenas que as magnitudes
  concordam, mas que as **formas** concordam. Um valor negativo, ou um
  valor estatisticamente consistente com zero, enfraquece o caso mesmo
  que a razão de módulos esteja próxima de um.

### 2.4 Colunas de sanidade física

- `magnetic_field_microgauss` — mediana de `|B|` inferida para a galáxia,
  em µG. Valores na banda `[0,1; 100]` µG são consistentes com a literatura
  para campos de disco em espirais e anãs. Valores fora dessa banda
  acionam um *audit flag* informacional `warn_b_band` (não bloqueante).
- `plasma_beta_thermal_median` — mediana de plasma-β sobre as cascas
  confiáveis. Valores `< 1` significam que a pressão magnética domina a
  pressão térmica (o regime em que o *ansatz* de gradiente de pressão em
  MHD ideal é localmente consistente); `> 1` significa que a pressão
  térmica domina. Um `β` alto **não** é, por si só, uma falha — apenas
  sinaliza que o diagnóstico está localmente fora de regime e não deve
  ser tomado como confirmatório para aquela galáxia.

### 2.5 Colunas de auditoria

- `audit_flags` — lista de strings na linha *supreme* correspondente; lista vazia significa *sem códigos bloqueantes* nos *gates* de validação. Os códigos exatos são explicados na §4 abaixo.
- `validation_status` — rótulo derivado em [`_validation_status_supreme`](../../src/pipeline/aggregation/global_aggregator.py): `SUPREME` (subconjunto de alta proveniência), `DATA_LIMITED` (feixe Planck maior que o diâmetro do disco *proxy*), ou `STANDARD`. **Não** confundir com a coluna `validation_is_reliable` do inventário (aqui todas as linhas já passaram fiabilidade + inclinação antes do filtro morfológico/janela de razão).

Nas tabelas `scenario_*.csv` do *grid*, o booleano **`is_reliable`** continua a significar “sem `fail_*` naquela célula \((\Delta x, f_{\mathrm{ion}})\)”.

### 2.6 Colunas do refinamento Fase 1

- `theoretical_gap_coverage_ratio_p16`, `theoretical_gap_coverage_ratio_p84`
  — percentis 16 e 84 da banda Monte Carlo da razão *headline* de cobertura.
  O MC propaga as incertezas 1-σ catalogadas em distância `D`, inclinação
  `i` e velocidade plana assintótica `V_flat` através do reescalonamento
  multiplicativo de forma fechada
  `R(D_s, i_s, V_s) = R_0\,(D_0/D_s)\,(\sin i_s/\sin i_0)^2\,(V_0/V_s)^2`,
  com `N = 300` amostras e semente `42` (determinístico).
- `convergence_ratio_to_wim_p16`, `convergence_ratio_to_wim_p84` — banda
  análoga sobre a razão de convergência ao WIM (a coluna *headline* do
  Discovery escalada pelo mesmo fator cinemático multiplicativo).
- `theoretical_gap_coverage_ratio_full_lorentz` — razão de cobertura quando
  o termo *toy* de tensão magnética (geometria de espiral logarítmica,
  ângulo de passo 12°) é somado à aceleração só-de-pressão. A coluna
  `theoretical_gap_coverage_ratio` mantém a convenção só-de-pressão; esta
  coluna reporta o orçamento "Lorentz total" paralelo para que os relatórios
  possam comparar os dois lado a lado.
- `accel_pressure_only_m_s2`, `accel_full_lorentz_m_s2` — mediana por
  galáxia dos perfis de aceleração só-de-pressão e Lorentz-total (SI).
- `rar_median_abs_residual_dex` — mediana de
  `|log10(g_obs / g_pred)|` sobre as cascas confiáveis contra a *Radial
  Acceleration Relation* (RAR) de McGaugh+2016. A estatística populacional
  sobre a amostra confiável aparece em `results/step_validation/validation_consistency_audit.md` (secção 7a).
- `btfr_residual_dex` — distância sinalizada (dex) da galáxia ao *ridge*
  baryonic Tully–Fisher de McGaugh+2012. `warn_btfr_outlier` é levantada
  quando `|residual| > 0,5` dex. A fração de *outliers* populacional
  aparece em `results/step_validation/validation_consistency_audit.md` (secção 7b).
- `warn_no_pa`, `warn_pure_pressure`, `warn_btfr_outlier` — etiquetas de auditoria informativas da Fase 1 (não bloqueantes; ver secção 4).

O ficheiro companheiro
`results/step_validation/validation_anchor_sensitivity_sweep.csv` lista, por galáxia,
`ratio_cmb` (a razão ancorada à CMB padrão) e
`ratio_turbulent_equipartition` (a razão recalculada contra um piso de
pressão por equipartição turbulenta), com `ratio_log10_delta`. Use este
ficheiro para auditar a razão *headline* contra a escolha da âncora
dimensional.

## 3. Grelha de cenários vs linha *headline*

Durante a **Fase 2 (Validação)** o motor avalia a grelha discreta completa em \((\Delta x, f_{\mathrm{ion}})\) para cada galáxia que passou o pré-filtro Discovery. Para **toda** a coorte fiável, o melhor cenário está no `*_summary.json` / inventário (`validation_reliable.csv`). A linha em **`validation_scientific_evidence_report.csv`** só existe para o subconjunto Im/LSB + janela de razão (§1) e resume essa procura: `best_scenario_name`, `scenario_gradient_scale_kpc` e `scenario_ionization_fraction` identificam a célula vencedora, enquanto `theoretical_gap_coverage_ratio` e `pearson_correlation` descrevem o fecho da lacuna em módulo e em forma.

O pipeline **já não** grava um PNG *heatmap* 3×3 autónomo dessa grelha. Para vistas de população (cobertura, Pearson, cruzamentos Discovery–Validação, resolução), abra **`results/graphs/population/`** após um *run* completo (ver [`Data_Interpretation_Guide_PT.md`](Data_Interpretation_Guide_PT.md), Galeria de figuras).

### Como ler os campos do melhor cenário

- **Cobertura perto da unidade** na linha vencedora significa que o *proxy* de pressão magnética casou, em ordem de grandeza, com o déficit cinemático para aquele par \((\Delta x, f_{\mathrm{ion}})\) — condição necessária, não suficiente.
- **Cobertura \(\ll 1\)** indica *sub-estima* sistemática do *proxy*.
- **Cobertura \(\gg 1\)** indica *super-estima* — frequentemente \(\Delta x\) demasiado agressivo ou \(f_{\mathrm{ion}}\) demasiado baixo para o modelo de MEI assumido.
- **Pearson \(r\)** perto de zero com boa cobertura módulo sugere casamento de amplitude sem casamento de forma radial; trate essa combinação como fisicamente suspeita.

A grelha sozinha **não** confirma a hipótese: apenas mapeia onde o teste de módulo é consistente. Os painéis em `graphs/population/` complementam o CSV mostrando como esses campos se distribuem na amostra processada.

## 4. Lendo `audit_flags`

Cada galáxia no *grid* carrega `audit_flags` como fonte de verdade dos *gates*.
Na **`validation_scientific_evidence_report.csv`**, as linhas são *snapshots* das
*supreme rows* (já fiáveis + inclinação) **depois** dos filtros morfológico e de janela
de razão (§1); podem ainda conter códigos `warn_*` informativos. A referência canónica
dos códigos é
[`src/pipeline/per_galaxy/metrics_builder.py`](../../src/pipeline/per_galaxy/metrics_builder.py)
e a secção correspondente em [`README.md`](../../README.md)
("Reliability audit").

| Código                     | Significado                                                                                       | Bloqueante? |
| -------------------------- | -------------------------------------------------------------------------------------------------- | ----------- |
| `fail_lat`                 | Latitude galáctica dentro da banda de contaminação por *foreground*.                                | sim         |
| `fail_reliable_fraction`   | Poucas cascas radiais sobreviveram ao *gate* de outliers MAD em `magnetic_pressure_acceleration`.  | sim         |
| `fail_density_clamped`     | Cascas demais caíram no piso de densidade em `GasDensityEstimator`.                                | sim         |
| `fail_stellar`             | Cascas demais dominadas por `V_disk`/`V_bulge`; a inversão de massa de gás perde significado.       | sim         |
| `fail_skipped_physics`     | O motor não conseguiu produzir `a_B` (Stokes, Δx ou `f_ion` faltando).                              | sim         |
| `fail_inclination`         | Inclinação catalogada abaixo do piso configurado; a propagação `V_obs / sin i` desestabiliza.       | sim         |
| `fail_calibration`         | A referência de Stokes caiu no valor de emergência; o `B` inferido não é medido.                    | sim         |
| `warn_b_band`              | Mediana de \|B\| fora da banda de plausibilidade 0,1–100 µG; apenas informacional.                  | não         |
| `warn_no_pa`               | Sem ângulo de posição de disco curado em `sources/morphology/disk_position_angle_overlay.json`; o motor recaiu sobre o gradiente axissimétrico. | não         |
| `warn_pure_pressure`           | Morfologia irregular / BCD; o *toy* de tensão magnética colapsa para zero, `accel_full_lorentz_m_s2 = accel_pressure_only_m_s2`.  | não         |
| `warn_btfr_outlier`        | Galáxia desviada mais de 0,5 dex do *ridge* BTFR de McGaugh+2012.                                                 | não         |

No **inventário** `validation_all.csv` / `validation_reliable.csv`, `validation_is_reliable`
é `True` se, e somente se, nenhum código `fail_*` estiver presente (`warn_*` não bloqueiam).
O resumo JSON em cada `galaxy_analysis/.../<NOME>/` espelha a mesma lista; **`results/step_validation/validation_consistency_audit.md`**
agrega a frequência de cada código em toda a amostra processada, mantendo visíveis os vieses sistemáticos da coorte rejeitada sem os misturar nas estatísticas *headline* das galáxias confiáveis.

## 5. Limpeza de Dados e Ruído

A leitura científica do pipeline depende de três regras de exclusão operacionais e da
**interpretação física** da coluna `required_density_for_match` na **Fase 1 (Discovery)**.
Os três primeiros subsecções cobrem *foreground*, cascas ruidosas e o censo Discovery;
a mediana populacional de \(\rho_{\mathrm{req}}\) **não** deve ser lida como “prova”
isolada sem `resolution_correction_factor_*` e sem o *sweep* de âncora (ver §5.3).

### 5.1 Corte por latitude galáctica ($|b| < 20^\circ$)

O pipeline emite o código `fail_lat` e rejeita a galáxia sempre que a
latitude galáctica satisfaz $|b| < 20^\circ$. O limiar canônico vive em
[`RELIABILITY_CUTS.MIN_GALACTIC_LATITUDE_DEG = 20.0`](../../src/config/thresholds.py).
A motivação é instrumental, não cosmológica: a banda equatorial da Via
Láctea concentra a poeira polarizada cuja emissão térmica satura o sinal
de Stokes na banda do Planck usada como entrada (353 GHz). Em latitudes
baixas, $\sqrt{Q^2 + U^2}$ medido no pixel central da galáxia-alvo é
dominado pelo *foreground* galáctico, não pelo campo magnético da
galáxia externa. A consequência é direta sobre o calibrador: a relação

\[
B(\hat{n}) \;=\; \sqrt{2\, \mu_0\, u_{\mathrm{CMB}}}\,\frac{S(\hat{n})}{S_{\mathrm{ref}}}
\]

deixa de medir o campo da galáxia e passa a medir o campo do disco da
Via Láctea ao longo da linha de visada. Reter essas galáxias produziria
um viés alto sistemático em `B`, `u_B` e — por construção — em `a_B`,
imitando um sinal positivo do modelo E-G por uma razão puramente
geográfica. O corte é, portanto, parte da definição do experimento.

### 5.2 Galáxia "ruidosa": critério operacional

Uma galáxia é classificada como ruidosa por dois mecanismos
complementares, ambos baseados em estatística robusta sobre o conjunto
de cascas radiais.

- **Outliers físicos por casca.** O *gate* `ReliabilityMaskBuilder`
  computa $\log_{10}(|a_B|/|a_c|)$ em cada casca e rejeita aquelas que
  caem fora de $\mathrm{mediana} \pm N\,\sigma_{\mathrm{MAD}}$, com
  $N = 3.0$ e $\sigma_{\mathrm{MAD}} = 1.4826 \times \mathrm{MAD}$. O
  parâmetro $N$ é
  [`MAG_PRESSURE_RATIO_OUTLIER_SIGMA = 3.0`](../../src/config/thresholds.py).
  Cascas fora dessa janela são *individualmente* descartadas: uma única
  casca anômala não invalida a galáxia inteira.
- **Filtro populacional por galáxia.** Se mais de
  $1 - \texttt{MIN\_RELIABLE\_MEASUREMENT\_FRACTION} = 0{,}30$ das
  cascas falharem no *gate* acima, a galáxia recebe
  `fail_reliable_fraction`, deixa de ser `validation_is_reliable` e **não**
  entra em `state.supreme_rows` — portanto **não** pode aparecer no CSV de evidência
  directa nem em `validation_reliable.csv`. A linha permanece em
  `validation_all.csv` para fins de auditoria.

A consequência prática é que uma curva de rotação contaminada por
*beam-smearing*, por *warps* não modelados ou por spikes locais de
densidade não desaparece silenciosamente: ela é exposta nos *flags* da
galáxia e ainda contribui para a "camada fantasma" do scatter global,
preservando a transparência sobre o que foi e o que não foi rejeitado.

### 5.3 Coluna `required_density_for_match` e leitura conjunta com WIM / resolução

Na **Fase 1 (censo Discovery)**, o motor usa $f_{\mathrm{ion}} = 1{,}0$ e o gradiente de
pressão magnética avaliado **na escala física do pixel Planck**,
$\Delta r_{\mathrm{phys}} = D_{\mathrm{gal}}\cdot\theta_{\mathrm{HEALPix}}$,
em vez do *knob* `scenario_gradient_scale_kpc` da grelha de validação. Cada galáxia gera
linhas nos quatro artefatos Discovery listados abaixo.

- [`results/step_discovery/discovery_all.csv`](../../results/step_discovery/discovery_all.csv)
  — todas as galáxias processadas, com motivo de exclusão preservado
  via `audit_flags`.
- [`results/step_discovery/discovery_reliable.csv`](../../results/step_discovery/discovery_reliable.csv)
  — subconjunto `discovery_is_reliable == True`; é o ficheiro a consultar para medianas
  de \(\rho_{\mathrm{req}}\) sobre a coorte Discovery limpa.
- [`results/step_discovery/discovery_not_reliable.csv`](../../results/step_discovery/discovery_not_reliable.csv)
  — complemento de `discovery_reliable` (mesmo *schema*).
- [`results/step_discovery/discovery_cohort_rho_band.csv`](../../results/step_discovery/discovery_cohort_rho_band.csv)
  — lista curta de *consistência interna*: linhas Discovery-*reliable* cuja
  `required_density_for_match` cai na banda multiplicativa
  $[0.7, 1.3]$ em torno da **mediana** de `required_density_for_match` no subconjunto *reliable*
  (ver `ResultPersister.DISCOVERY_ELITE_BAND_RELATIVE_MIN/MAX` em
  [`exporter.py`](../../src/reporting/persistence/exporter.py)).
  Quando essa banda dinâmica não pode ser formada, não há coorte adicional: todas as
  linhas relevantes recebem `discovery_is_elite = False` e `discovery_elite_band_source = none`.

A coluna principal, `required_density_for_match`, responde à pergunta
direta: qual densidade $\rho_{\mathrm{req}}$ o plasma precisaria ter para
que o gradiente de pressão magnética **reconstruído** fechasse o *gap* casca a casca?

\[
\rho_{\mathrm{req}}(r) \;=\;
\frac{\bigl|\,\nabla u_B(r)\,\bigr|}{\bigl|\,\Delta a(r)\,\bigr|}.
\]

A coluna `delta_r_m` regista $\Delta r_{\mathrm{phys}}$ em metros (SI). A coluna legacy
`median_required_density_for_match_kg_m3` duplica `required_density_for_match` por
compatibilidade com análises antigas.

A coluna `convergence_ratio_to_wim` compara \(\rho_{\mathrm{req}}\) à densidade de referência
**fixa** do WIM \(10^{-27}\,\mathrm{kg\,m^{-3}}\) (ordem de grandeza clássica do meio ionizado
difuso; ver literatura de emissão $\mathrm{H}\alpha$). **Em corridas SPARC+Planck completas
documentadas após a Fase 1**, a mediana populacional de `required_density_for_match` sobre
`discovery_is_reliable` situa-se tipicamente **muitas ordens de magnitude acima** de
\(10^{-27}\,\mathrm{kg\,m^{-3}}\) — logo `convergence_ratio_to_wim` \(\gg 1\) é o estado
empírico habitual, não a “coincidência” automática com o WIM. Isso **não** invalida o
formalismo; obriga a leitura conjunta com:

- `beam_size_kpc`, `resolution_correction_factor_discovery`, `resolution_correction_factor_grid`
  (diluição finita do gradiente por feixe e por escala de reconstrução);
- `validation_anchor_sensitivity_sweep.csv` (mudança de dex entre âncora CMB e piso de
  equipartição turbulenta no calibrador Stokes→Tesla);
- `validation_anchor_sensitivity_sweep.csv` (desvio em dex entre a âncora dimensional CMB e um
  patamar de equipartição turbulenta na ponte Stokes→Tesla);
- a Secção 5 de [`EG_MODEL_THEORY_PT.md`](../theory/EG_MODEL_THEORY_PT.md) e os guias em
  [`docs/results_interpretation/`](../results_interpretation/) sobre colunas de resolução e leitura conjunta com o *sweep* de âncoras.

Interpretação física útil: quando \(\rho_{\mathrm{req}}\) **se aproxima** de \(10^{-27}\,\mathrm{kg\,m^{-3}}\)
*depois* destes diagnósticos, o argumento de acoplamento ao plasma ionizado difuso permanece
plausível; quando permanece elevada, o resultado é **informação sobre o ansatz + dados**,
não uma confirmação silenciosa do WIM. O HI frio do disco (\(\rho_{\mathrm{HI}} \sim 10^{-21}\,\mathrm{kg\,m^{-3}}\))
continua a marcar o limite superior “ingénuo” — o canal MHD do modelo não pretende substituir
esse gás denso como responsivo ao gradiente na mesma formulação.

Ao terminar a execução, o pipeline imprime no console a mediana de `required_density_for_match`
sobre o subconjunto Discovery-*reliable*, a contagem `N_passed` e `N_elite` da banda dinâmica,
como verificação rápida da ordem de grandeza.

### 5.4 Discovery como base empírica do grid de validação

O fluxo do projeto trata o censo Discovery e o grid como etapas sequenciais de um
mesmo argumento, não como execuções isoladas:

1. **Censo Discovery** responde *"qual densidade o sinal magnético pede?"* sem *knobs* de
   grelha — Δx é o pixel Planck físico e $f_{\mathrm{ion}}=1$. Os valores absolutos de
   \(\rho_{\mathrm{req}}\) devem ser confrontados com o WIM **como referência**, não como
   destino garantido da mediana populacional (§5.3).
2. **Grelha (Fase 2 — Validação)** pergunta *"como varia o fecho do *gap* quando se permitem escalas
   efectivas de gradiente e frações ionizadas diferentes no denominador \(\rho_{\mathrm{ion}}\)?"*.
   As escalas `scenario_gradient_scale_kpc ∈ {0,050; 0,100; 0,500}` kpc e as frações
   `scenario_ionization_fraction ∈ {0,001; 0,010; 0,100}` são um *prior* explícito sobre
   estrutura magnética de pequena escala e grau de ionização; não substituem o Discovery.

Em termos práticos: `discovery_reliable.csv` ancora a ordem de grandeza do problema antes da
grelha; `discovery_cohort_rho_band.csv` lista exemplos **internamente consistentes** com a
mediana de \(\rho_{\mathrm{req}}\) na população Discovery-*reliable*.

**Não confundir com coortes de validação:** filtre `step_validation/validation_all.csv` por
`validation_is_reliable` para estatísticas *headline* pós-*gates* do *grid* (universo tipicamente
mais pequeno que Discovery-*reliable*). As colunas `discovery_is_reliable` / `discovery_is_elite`
fundidas nas tabelas de validação conservam sempre a semântica **Discovery**.

Mudar a definição de \(\rho_{\mathrm{req}}\) na §4.2 do
[`Draft_Paper_EG_Model_PT.md`](../academic/Draft_Paper_EG_Model_PT.md)
ou as constantes da banda *elite* (`DISCOVERY_ELITE_BAND_*`) exige actualizar
simultaneamente este parágrafo e o espelho EN.

## 6. O que uma razão de cobertura próxima da unidade *não* prova

Mesmo que `theoretical_gap_coverage_ratio ≈ 1` para muitas galáxias
confiáveis e a correlação de Pearson seja positiva, o resultado ainda não
desloca a interpretação de matéria escura. As limitações a seguir
devem ser sanadas antes de uma reivindicação estatisticamente robusta (Secção 5 de
[`EG_MODEL_THEORY_PT.md`](../theory/EG_MODEL_THEORY_PT.md) e esta Secção 6); várias já têm **mitigação parcial** no código actual — o ponto em aberto é o grau de robustez perante um revisor externo:

1. **Gradiente no disco vs. céu.** O motor já estima `||∇u_B||` com LSQ no plano do disco
   quando PA e inclinação existem (`disk_plane_pressure_gradient`); quando o PA falta,
   recai-se no modo axissimétrico (`warn_no_pa`). Ainda assim, a mistura entre escala HEALPix,
   `scenario_gradient_scale_kpc` e física do disco permanece uma limitação de leitura — ver
   comentários em [`PhysicsEngine`](../../src/physics/engine.py).
2. **Incerteza na razão *headline*.** Já existe propagação Monte Carlo fechada em \((D, i, V_{\mathrm{flat}})\)
   para `theoretical_gap_coverage_ratio_p16` / `_p84`. Isto **não** substitui um posterior completo
   nem propaga incertezas correlacionadas no mapa Stokes ou na calibração local de \(B\).
3. **Âncora Stokes→Tesla.** O *sweep* `validation_anchor_sensitivity_sweep.csv` quantifica o
   deslocamento em dex entre âncora CMB e piso de equipartição turbulenta; **não** varre
   âncoras observacionais locais (síncrotron) nem calibradores externos pixel a pixel.
4. **Tensão magnética.** O ramo *toy* de espiral logarítmica está implementado
   (`theoretical_gap_coverage_ratio_full_lorentz`); é desactivado para Im/BCD (`warn_pure_pressure`).
   Falta geometria 3D amostrada para tensão *predictiva* fora do *toy*.
5. **Controles sintéticos.** Controles NFW e similares limitam apenas classes específicas de
   falso positivo; não excluem mecanismos não-MHD nem todos os artefactos instrumentais.

Uma razão de cobertura próxima da unidade é, portanto, um desfecho
**necessário** mas não **suficiente**: ela autoriza trabalho adicional;
não, por si só, derruba ΛCDM.

## 7. O que uma razão de cobertura falha *prova*

Se os números principais voltam ordens de magnitude longe da unidade e
a correlação de Pearson é consistente com zero ao longo das galáxias
confiáveis, a conclusão é **definitiva para aquele cenário**: o *ansatz*
de gradiente de pressão em MHD ideal, com a âncora configurada e os Δx
e `f_ion` configurados, não é o mecanismo por trás da anomalia da curva
de rotação. Um resultado negativo desse tipo é, em termos científicos,
tão valioso quanto um positivo — remove um candidato do espaço de
hipóteses sem exigir um experimento externo.

## 8. Para onde olhar em seguida

Depois de `validation_reliable.csv` e do *audit*, use
`validation_scientific_evidence_report.csv` apenas quando o seu argumento for explicitamente
Im/LSB **e** a janela de razão 0,2–5,0 for relevante.

- **`results/galaxy_analysis/{validation_scientific_evidence|validation_reliable|validation_not_reliable}/<NOME>/`** — comparação radial de acelerações, RAR e acelerações assinadas para cada galáxia processada naquele *run* (ao lado do CSV/JSON em disco; o subdiretório `validation_scientific_evidence` reúne a coorte listada no CSV de evidência directa após o passo de publicação).
- **`results/step_validation/validation_consistency_audit.md`** — auditoria
  populacional incluindo a frequência de cada código `audit_flags` na
  execução, para que vieses sistemáticos na amostra rejeitada fiquem
  visíveis.
- **Narrativa Markdown opcional** — `NarrativeGenerator` em `src/reporting/narrative.py` existe para texto *template* a partir de linhas fiáveis, mas **não** é escrito automaticamente para `results/` pela CLI; para relatórios permanentes use este guia + CSVs + *audit*.
