# Um Teste Magnetohidrodinâmico Falsificável para Anomalias de Curvas de Rotação Galáctica: o Modelo E-G

> **Status: rascunho de trabalho.** Modelo estrutural para uma submissão
> futura a um periódico revisado por pares (e.g. MNRAS ou A&A). Os resultados
> numéricos que populam §Metodologia e qualquer seção de Resultados futura
> são produzidos pelo pipeline open-source deste repositório; eles não são
> embutidos aqui para que o rascunho permaneça desacoplado de qualquer
> execução individual. A versão em inglês está em
> [`Draft_Paper_EG_Model_EN.md`](Draft_Paper_EG_Model_EN.md).

---

## Resumo (Abstract)

Os perfis de rotação planos de galáxias espirais e anãs são convencionalmente
interpretados como evidência de um halo de matéria escura não-bariônica e
sem colisões, que fornece a aceleração centrípeta faltante. Este trabalho
testa uma contribuição alternativa para o orçamento de força radial: um
gradiente de pressão magnética magnetohidrodinâmico (MHD) atuando sobre o
componente ionizado morno do meio interestelar (WIM). Para cada galáxia do
catálogo SPARC, o pipeline (i) calcula o déficit cinemático de aceleração
`Δa(r) = V_obs²(r)/r − V_bar²(r)/r`; (ii) inverte o componente de rotação
só-de-gás em uma densidade volumétrica de gás ionizado via decomposição em
cascas cilíndricas; (iii) calibra a amplitude de Stokes de polarização
linear do Planck na posição celeste da galáxia em uma intensidade de campo
magnético via uma âncora dimensional baseada na densidade de energia da CMB;
e (iv) avalia a aceleração de pressão magnética
`a_B(r) = −∇(B²/2μ₀)/ρ_ion(r)` casca a casca. Reporta-se a razão de módulos
no disco externo `|a_B|/|Δa|` e a correlação de Pearson entre os perfis
radiais de `a_B` e `Δa` sobre uma grade de cenários que abrange três escalas
de gradiente (50, 100, 500 pc) e três frações de ionização (0,1, 1, 10 %).
O pipeline é exercitado contra controles negativos NFW sintéticos e reporta
uma auditoria explícita de confiabilidade por galáxia e por cenário. O
resultado é enquadrado como um teste de falsificação do *ansatz* de pressão
MHD, não como uma refutação da matéria escura enquanto fenômeno.

**Palavras-chave:** galáxias: cinemática e dinâmica — galáxias: campos
magnéticos — magnetohidrodinâmica (MHD) — meio interestelar: campos
magnéticos — métodos: análise de dados.

---

## 1. Introdução

Curvas de rotação de galáxias espirais e anãs exibem uma anomalia empírica
robusta: em raios onde a distribuição de massa visível (bariônica) prevê um
decaimento Kepleriano `V ∝ r^{−1/2}`, medidas espectroscópicas mostram, de
forma consistente, perfis de velocidade planos (e.g. Rubin & Ford 1970;
Bosma 1981). Dentro do arcabouço cosmológico padrão (ΛCDM), essa anomalia
é fechada por um halo de matéria escura sem colisões e esférico cuja
gravidade fornece a aceleração centrípeta faltante.

Campos magnéticos galácticos com intensidades características de um a dez
microgauss estão observacionalmente estabelecidos em galáxias espirais e
anãs por emissão síncrotron difusa, rotação de Faraday de fontes polarizadas
de fundo e os mapas de polarização linear do satélite Planck (Planck
Collaboration 2016). A contribuição dinâmica desses campos para o problema
das curvas de rotação é, no entanto, tipicamente tratada como uma correção
de segunda ordem a uma descrição apenas gravitacional.

O arcabouço E-G (Eletrodynamic–Gravity) aqui apresentado **não** é uma
proposta de gravidade modificada. É um teste quantitativo de uma única
afirmação falsificável: a de que o gradiente de pressão magnética
`−∇(B²/2μ₀)`, atuando sobre o componente ionizado morno do meio interestelar,
contribui com uma aceleração radial da ordem de
`Δa(r) = V_obs²(r)/r − V_bar²(r)/r` nos raios onde halos de matéria escura
são convencionalmente invocados. O teste produz estimativas de ordem de
magnitude para `a_B(r)` a partir de *inputs* observacionais independentes
(a cinemática SPARC de um lado e os mapas de polarização do Planck do outro)
e reporta o desfecho com o mesmo rigor para confirmação e falsificação.

O artigo está organizado da seguinte forma. A Seção 2 enuncia a hipótese
como uma decomposição em dois regimes do orçamento de força radial. A
Seção 3 deriva o *ansatz* de gradiente de pressão magnética a partir da
força de Lorentz volumétrica e discute o acoplamento com o WIM. A Seção 4
descreve os dados, o esquema de calibração e o pipeline de estimadores por
casca. A Seção 5 discute o protocolo de falsificação e a auditoria explícita
de confiabilidade. A Seção 6 enumera as limitações conhecidas e os
refinamentos planejados.

---

## 2. Hipótese — Decomposição em Dois Regimes

A aceleração radial total no raio galactocêntrico `r` é decomposta como

```text
a_total(r) = a_grav(r) + a_B(r),
```

onde `a_grav(r) = V_bar²(r)/r` é a contribuição Newtoniana bariônica da
decomposição de massa SPARC `V_bar = √(V_gas² + V_disk² + V_bulge²)` e
`a_B(r)` é a contribuição de pressão MHD definida na Seção 3.

Dois regimes são distinguidos pelo peso relativo dos dois termos:

- **Regime A (compacto/bariônico).** Nos poucos kpc internos de um disco
  galáctico e em sistemas estelares, `a_grav` domina: `|a_grav| ≫ |a_B|`.
  A gravidade Newtoniana é recuperada como limite local e o arcabouço E-G
  se reduz trivialmente a ela.
- **Regime B (estendido/magnetizado).** No disco externo, onde a lacuna
  cinemática `Δa(r) = V_obs²(r)/r − V_bar²(r)/r` se abre, a questão é se
  `a_B(r)` está na ordem de `Δa(r)`. Uma resposta positiva identificaria
  o gradiente de pressão MHD como um contribuinte viável para a aceleração
  faltante; uma resposta negativa (em particular `|a_B| ≪ |Δa|` por ordens
  de magnitude) falsificaria este *ansatz* específico.

A hipótese **não** afirma que a totalidade da anomalia da curva de rotação
seja de origem eletromagnética; testa-se a magnitude de uma contribuição
específica que a análise convencional rotineiramente negligencia.

---

## 3. Mecanismo Proposto — Gradiente de Pressão MHD sobre o WIM

### 3.1 Da força de Lorentz microscópica à pressão magnética

Uma única partícula carregada de carga `q` e velocidade `v` em um campo
eletromagnético `(E, B)` sofre a força de Lorentz `F = q(E + v × B)`. Para
um plasma interestelar, as escalas de tempo relevantes são longas comparadas
ao período de gira e a grandeza dinâmica de interesse é a aceleração de
massa de um elemento de fluido. Somando sobre o conteúdo de carga de uma
parcela de fluido e assumindo quase-neutralidade, a força de Lorentz
volumétrica por unidade de volume se reduz, no limite MHD não-relativístico,
a `f = J × B = (∇ × B) × B / μ₀`. Uma identidade vetorial padrão decompõe
isso em um termo de gradiente de pressão `−∇(B²/2μ₀)` e um termo de tensão
magnética `(B · ∇)B/μ₀` (Jackson 1999).

### 3.2 O *ansatz* de gradiente de pressão em MHD ideal

Este artigo retém apenas o componente de gradiente de pressão e trata o
plasma como ideal (alta condutividade, sem resistividade, sem tensor de
pressão anisotrópico, sem contribuição de raios cósmicos, sem pressão
turbulenta). A aceleração que o componente ionizado do meio interestelar
sofre por causa desse gradiente é

```text
a_B(x) ≈ − ∇u_B(x) / ρ_ion(x),     com u_B(x) = B(x)² / (2 μ₀).
```

A justificativa para descartar o termo de tensão é discutida como um eixo
de refinamento na §6.

### 3.3 O WIM como componente responsivo

A densidade de massa no denominador **não** é a densidade total de gás. É
a densidade do componente que responde eletromagneticamente na escala de
tempo relevante, modelada como

```text
ρ_ion(r) = f_ion · ρ_gas(r),
```

com `f_ion ∈ (0, 1]` uma fração de ionização livre. O componente
observacionalmente relevante é o Meio Ionizado Morno (WIM): hidrogênio
ionizado a `T ~ 10⁴ K` permeando discos e halos galácticos com densidades
de massa características ordens de magnitude abaixo do gás neutro HI
detectado por surveys de 21 cm (Reynolds 1990; Haffner et al. 2009).

### 3.4 Geometria de escala

A aceleração gravitacional de uma massa bariônica concentrada centralmente
decai como `1/r²`. O gradiente de pressão magnética em uma geometria
filamentar ou de disco decai substancialmente mais devagar (`B ∝ 1/r` ao
longo de um filamento reto idealizado conduzindo corrente, com o gradiente
de pressão correspondente não-trivial em distâncias inter-filamentares).
A competição entre as duas contribuições depende, portanto, explicitamente
da geometria radial, e o raio de cruzamento é, em princípio, uma grandeza
derivada do modelo, não um parâmetro livre.

---

## 4. Metodologia

> **Nota código-como-verdade.** Os nomes de classe referenciados abaixo são
> os estimadores Single-Responsibility implementados em `src/physics/`. O
> `.cursorrules` do repositório impõe que qualquer mudança nesses estimadores
> dispare uma sugestão de sincronização contra esta seção, mantendo
> metodologia e código coerentes.

### 4.1 Dados

- **Cinemática SPARC.** O catálogo Spitzer Photometry and Accurate Rotation
  Curves (Lelli, McGaugh & Schombert 2016) fornece decomposições por raio
  de `V_obs`, `V_gas`, `V_disk` e `V_bulge` para 175 galáxias espirais e
  anãs, além de inclinação, distância e morfologia.
- **Polarização do Planck.** Mapas de Stokes `Q` e `U` da emissão de poeira
  galáctica do Planck Legacy Archive (Planck Collaboration 2016), acessados
  via a pixelização HEALPix (Górski et al. 2005). O pipeline usa leituras
  memory-mapped da carga útil FITS para escalar até o céu inteiro.
- **NED.** O NASA/IPAC Extragalactic Database é consultado para coordenadas
  celestes `(RA, Dec)`, redshift e pistas de distância não presentes no
  SPARC, com cache em disco para manter a execução reproduzível.

### 4.2 Estimadores

O pipeline compõe quatro estimadores Single-Responsibility, orquestrados
pelo `PhysicsEngine`:

- `BaryonicGapCalculator` avalia `a_obs(r) = V_obs²/r`,
  `a_bar(r) = V_bar²/r` com `V_bar² = V_gas² + V_disk² + V_bulge²`, e o
  `acceleration_gap = a_obs − a_bar` por casca. A incerteza de inclinação
  é propagada linearmente em uma banda de ±1σ sobre `a_obs`.
- `GasDensityEstimator` inverte o componente de rotação só-de-gás em uma
  densidade volumétrica de gás ionizado. Dois caminhos equivalentes
  `M_gas^A = V_gas²·r/G` e
  `M_gas^B = max(V_bar² − V_disk² − V_bulge², 0)·r/G` são avaliados e
  reconciliados com `max(A,B)` para impedir super-inflação; a densidade
  volumétrica decorre de uma decomposição em cascas cilíndricas com volume
  de casca `ΔV_i = π(r_i² − r_{i−1}²)·h`. Um piso de densidade configurável
  e um *gate* de dominância estelar suprimem cascas onde a inversão é
  desprovida de significado.
- `MagneticFieldCalibrator` converte a amplitude de Stokes
  `S = √(Q² + U²)` na posição celeste da galáxia em uma intensidade de
  campo magnético em Tesla `B = √(2 μ₀ u_ref) · S / S_ref`, onde `u_ref`
  é a densidade de energia da CMB (uma âncora dimensional, não uma
  afirmação física sobre a origem do campo de disco) e `S_ref` é a mediana
  de `S` sobre o mapa de céu carregado (calibração auto-referencial). A
  escolha da âncora é um eixo explícito de sensibilidade.
- `MagneticPressureGradientEstimator` estima `∇u_B` como uma diferença
  finita escalar entre o pixel HEALPix central e a média de seus vizinhos
  angulares, dividida por um passo físico `Δx`. Dois valores paralelos de
  `Δx` são rastreados: um botão de cenário
  `scenario_gradient_scale_kpc ∈ {0.05, 0.1, 0.5}` e um físico
  `Δx_phys = D_galaxy · θ_neighbour` derivado da distância catalogada e do
  espaçamento angular dos vizinhos HEALPix.
- **Gradiente deprojetado no referencial do disco (refinamento Fase 1).**
  No caminho Δx-físico, a diferença escalar `Δu_B/Δx` é substituída por um
  gradiente direccional obtido a partir de um ajuste planar por mínimos
  quadrados `u_{B,i} = u_{B,0} + g_x\,\Delta x_i + g_y\,\Delta y_i` sobre o
  anel de vizinhos HEALPix. Quando o ângulo de posição do disco (PA) e a
  inclinação `i` estão ambos disponíveis, os deslocamentos no plano do céu
  são rotacionados de modo que o eixo maior alinhe com `+x_{disk}` e o
  componente do eixo menor é reescalado por `1/\cos i`. A magnitude do
  gradiente reportada é a norma no referencial do disco
  `\sqrt{g_x^2 + g_y^2}`. Quando o PA não está disponível, o motor recai
  sobre a magnitude axissimétrica (com a correção linear `\cos i`) e levanta
  `warn_no_pa`.
- **Termo *toy* de tensão magnética** (`TensionTermEstimator`,
  `src/physics/diagnostics/tension.py`). A componente radial de
  `(\mathbf B\cdot\nabla)\mathbf B/\mu_0` é avaliada sob uma geometria de
  linha de campo em espiral logarítmica com ângulo de passo `p`, para a
  qual o raio de curvatura reduz-se a `R_c(r) = r/\sin^2 p`. A aceleração
  reportada é então
  `a_{\mathrm{tension},r}(r) = -\cos p\cdot B^2(r)/[\mu_0\,\rho_{\mathrm{ion}}(r)\,R_c(r)]`.
  Morfologias sem espiral coerente (`Im`, `BCD`) forçam o termo a zero e
  levantam `warn_pure_pressure`. O motor escreve as colunas por casca lado a
  lado `accel_pressure_only_m_s2` (alias da aceleração de pressão no Δx do
  cenário) e
  `accel_full_lorentz_m_s2 = accel_pressure_only + a_{\mathrm{tension},r}`,
  junto a uma `theoretical_gap_coverage_ratio_full_lorentz` paralela para
  permitir auditar o impacto do termo de tensão *no mesmo Δx da razão
  *headline*.
- **Propagação Monte Carlo de incertezas**
  (`src/physics/diagnostics/uncertainty.py`). As incertezas 1-σ catalogadas
  em distância `D`, inclinação `i` e velocidade plana assintótica
  `V_{\rm flat}` são propagadas para `theoretical_gap_coverage_ratio` e
  `convergence_ratio_to_wim` através de uma forma fechada multiplicativa
  `R(D_s, i_s, V_s) = R_0\,(D_0/D_s)\,(\sin i_s/\sin i_0)^2\,(V_0/V_s)^2`,
  com `N=300` amostras e semente `42`. O resultado é exposto como
  `theoretical_gap_coverage_ratio_p16`, `theoretical_gap_coverage_ratio_p84`,
  `convergence_ratio_to_wim_p16`, `convergence_ratio_to_wim_p84` no CSV
  *headline*.
- **Sweep de sensibilidade da âncora**
  (`src/physics/diagnostics/anchor_sensitivity.py`). A âncora de calibração
  MHD `u_{\rm ref}` é variada entre a densidade de energia da CMB (padrão)
  e um piso de pressão por equipartição turbulenta
  `u_{\rm turb} = 0.5\,\rho_{\rm gas}\,\sigma_v^2`. As duas razões resultantes
  por galáxia são persistidas em uma tabela dedicada
  `results/step_validation/validation_anchor_sensitivity_sweep.csv` para auditar a
  razão *headline* contra a escolha da âncora.

O motor também reporta uma `required_density_for_match` por casca, isto é,
a densidade de plasma ionizado `ρ_req` que fecharia a lacuna cinemática
dado o gradiente de pressão magnética observado:

\[
\rho_{\mathrm{req}}(r) \;=\;
\frac{\bigl|\,\nabla u_B(r)\,\bigr|}{\bigl|\,\Delta a(r)\,\bigr|}.
\]

Quando `Δx_phys` está disponível (e é finito), o `∇u_B` no numerador é
tomado na *separação física* dos vizinhos HEALPix, e a máscara de cascas
*underperforming* que seleciona os pontos contribuintes à mediana usa
`|a_B(Δx_phys)|`. Quando `Δx_phys` está ausente, o motor recai sobre o
gradiente exploratório do cenário e sobre o `a_B` correspondente. Essa
ramificação é necessária para que a **Fase 1 (censo Discovery)** (que fixa
`Δx = Δx_phys` e `f_ion = 1`) ancore a verificação de convergência ao WIM
no pixel físico real do Planck, enquanto a **Fase 2 (Validação — grelha de cenários)** continua
explorando o *prior* (Δx, f_ion).

A tabela populacional `validation_scientific_evidence_report.csv` (escrita após o grid) funde, para cada
galáxia que passa na banda de razão de evidência directa, os campos de cabeçalho do Discovery
(`required_density_for_match`, `convergence_ratio_to_wim`, `delta_r_m`, `physical_pixel_scale_kpc`)
juntamente com os booleanos `discovery_is_reliable` e `discovery_is_elite`, lado a lado com a
melhor combinação de cenário `theoretical_gap_coverage_ratio` e `pearson_correlation`.

### 4.3 Grade de cenários

Por galáxia, o motor avalia o produto cartesiano de
`scenario_gradient_scale_kpc ∈ {0.05, 0.1, 0.5}` (kpc) e
`scenario_ionization_fraction ∈ {0.001, 0.01, 0.1}`, totalizando nove
cenários. Para cada cenário, a `magnetic_pressure_acceleration` por casca
é comparada ao `acceleration_gap` por casca via a razão de módulos no
disco externo
`theoretical_gap_coverage_ratio = |⟨a_B⟩_edge| / |⟨gap⟩_edge|` e a
correlação de Pearson dos dois perfis radiais sobre as cascas confiáveis.

### 4.4 Auditoria de confiabilidade

Cada galáxia carrega uma lista explícita `audit_flags` que documenta os
*gates* em que falhou: `fail_lat` (latitude de *foreground*),
`fail_reliable_fraction` (poucas cascas sobreviveram ao *gate* de outliers),
`fail_density_clamped` (densidade de gás *clamped* no piso na maioria das
cascas), `fail_stellar` (gravidade estelar domina a cinemática),
`fail_skipped_physics` (o motor não conseguiu produzir `a_B`),
`fail_inclination` (inclinação catalogada abaixo de um limiar),
`fail_calibration` (a referência de Stokes caiu em um valor de emergência).
*Warnings* não bloqueantes incluem `warn_b_band` (`B` inferido fora da
banda de plausibilidade da literatura), `warn_no_pa` (sem ângulo de posição
de disco curado — o motor recaiu sobre o gradiente axissimétrico),
`warn_pure_pressure` (morfologia irregular / BCD — o *toy* de tensão magnética
colapsa para zero) e `warn_btfr_outlier` (galáxia desviada do *ridge*
baryonic Tully–Fisher de McGaugh+2012 por mais que o limiar configurado em
dex). O booleano derivado `is_reliable` em cada **linha de cenário** é verdadeiro se, e somente se,
nenhum código `fail_*` estiver presente; códigos `warn_*` são informativos
e não rebaixam o `is_reliable` daquele cenário. O inventário ao nível da galáxia usa
`validation_is_reliable` em `validation_all.csv` (e em
`data_reliability_assessment` dentro de cada `*_summary.json`). As tabelas do censo Discovery (Fase 1) usam
`discovery_is_reliable` (`results/step_discovery/discovery_reliable.csv` / `discovery_not_reliable.csv`). A auditoria de consistência também surfaca
dois diagnósticos populacionais sobre a amostra confiável: o resíduo
absoluto mediano contra a *Radial Acceleration Relation*
(`rar_median_abs_residual_dex`) e a fração de *outliers* BTFR
(`warn_btfr_outlier`).

### 4.5 Controles negativos sintéticos

O pipeline também é exercitado em curvas de rotação sintéticas geradas a
partir de halos de matéria escura Navarro–Frenk–White (NFW) (Navarro, Frenk
& White 1996). Por construção, essas curvas não têm gradiente magnético
real; o pipeline deve, portanto, *falhar* em encontrar uma correspondência
de pressão magnética para elas. Esse único teste limita uma classe específica
de falsos positivos.

---

## 5. Protocolo de Falsificação

O modelo é considerado **falsificado no cenário configurado** se ocorrer
qualquer um dos seguintes:

1. `theoretical_gap_coverage_ratio` difere da unidade em mais de uma ordem
   de magnitude na população de galáxias confiáveis; ou
2. a correlação de Pearson entre os perfis radiais de `a_B` e `Δa` é
   consistente com zero sobre os raios confiáveis.

O modelo é considerado **consistente no cenário configurado** se ambas as
grandezas ficarem na banda descrita na §4 simultaneamente sobre uma fração
não-trivial do subconjunto confiável. Consistência em um cenário não impede
refinamento; em particular, os itens da §6 devem ser endereçados antes que
uma razão de cobertura próxima da unidade possa ser considerada
*estatisticamente* distinguível do acaso.

O protocolo se aplica igualmente à confirmação e à falsificação, e o pipeline
reporta ambos os desfechos com a mesma trilha de auditoria.

---

## 6. Limitações e Refinamentos Planejados

Uma enumeração completa está na Secção 5 de [`EG_MODEL_THEORY_PT.md`](../theory/EG_MODEL_THEORY_PT.md)
e nos guias em [`docs/results_interpretation/`](../results_interpretation/). Os itens principais são:

1. **Deprojeção para o plano do disco.** O gradiente escalar `Δu_B / Δx` é
   projetado no céu; o componente dinamicamente relevante é o radial no
   plano do disco. Uma amostragem em anel de múltiplos vizinhos com
   derivada direcional deprojetada é o próximo passo.
2. **Propagação Bayesiana de incertezas.** Incertezas de inclinação,
   distância e `V_flat` devem ser propagadas através da razão de cobertura
   por Monte-Carlo antes que qualquer valor próximo da unidade seja
   interpretado como estatisticamente significativo.
3. **Cross-check de massa de gás.** A integral
   `gas_mass_density × shell_volume` deve ser validada contra a massa de
   HI SPARC `MHI_1e9_solar`.
4. **Sensibilidade à âncora.** A âncora baseada na densidade de energia da
   CMB é dimensional; âncoras alternativas (equipartição turbulenta,
   equipartição de raios cósmicos, B fixo de disco galáctico) devem ser
   varridas e reportadas.
5. **Tensão magnética.** O *ansatz* só-pressão descarta o termo de tensão
   `(B · ∇)B/μ₀`; uma geometria-brinquedo de linhas de campo em espiral
   logarítmica permitiria ao motor reportar colunas `pressure-only` e
   `pressure+tension` lado a lado.
6. **Pressão de raios cósmicos.** Para anãs formadoras de estrelas, a
   pressão de raios cósmicos pode competir com a térmica e deslocar o
   diagnóstico plasma-β.
7. **Auto-consistência RAR/BTFR.** Galáxias para as quais o modelo de
   pressão magnética é declarado consistente devem também satisfazer a
   Radial Acceleration Relation empírica (McGaugh, Lelli & Schombert 2016)
   e a Baryonic Tully–Fisher Relation (McGaugh et al. 2000).
8. **Aposentadoria do botão de cenário.** Uma vez implementadas a
   deprojeção em disco e a propagação Bayesiana, a grade legada
   `scenario_gradient_scale_kpc` se torna uma varredura de sensibilidade,
   não o resultado primário.

Esses itens não bloqueiam a afirmação de falsificabilidade da §5, mas são
necessários antes que qualquer alegação de razão de cobertura próxima da
unidade possa ser promovida a um resultado cosmológico positivo.

---

## 7. Itens a Revisitar — Predições Adjacentes

Três linhas de investigação ficam na fronteira da análise atual e merecem um
acompanhamento dedicado. Nenhuma delas é necessária para o teste de
falsificabilidade da §5, mas um resultado positivo daquele teste tornaria cada
uma delas consideravelmente mais interessante.

### 7.1 Plasma Lensing

Caso halos estendidos contenham um reservatório substancial de plasma
ionizado morno, fótons atravessando-os experimentam uma deflexão refrativa
dependente da frequência, ausente da lente puramente gravitacional. No limite
de plasma frio, o índice de refração é

\[
n(\omega) \;=\; \sqrt{1 - \omega_p^{\,2}/\omega^{2}}, \qquad
\omega_p \;=\; \sqrt{\dfrac{n_e\, e^{2}}{\varepsilon_0\, m_e}},
\]

de modo que o ângulo de deflexão escala como `Δψ ∝ ∇⊥ n_e / ν²` e é, portanto,
**cromático**, em contraste com a deflexão gravitacional acromática da
Relatividade Geral. Para um halo com `n_e ∼ 10⁻³ cm⁻³` e comprimento de
gradiente de 10 kpc, isso prevê deslocamentos mensuráveis dependentes de
frequência em fontes rádio de fundo nas bandas GHz–MHz, acessíveis com
interferometria VLBI multi-frequência e arrays de timing de pulsares. Um
resultado nulo nas linhas de visada mais limpas (afastadas do plano galáctico e
de aglomerados interpostos) apertaria os limites sobre a densidade eletrônica
integrada na coluna ao longo das linhas de visada relevantes e, portanto, sobre
o orçamento de massa do WIM invocado na §3.3.

Referências: Pen & King (2012); Cordes et al. (2017); Vedantham (2019).

### 7.2 Primordial Magnetogenesis

O *ansatz* MHD da §3 toma a existência de um campo magnético coerente da ordem
de `B ∼ 1`–`10 μG` como uma *condição de contorno* amostrada pelos mapas de
polarização do Planck. Sua **origem** fica em aberto pela análise atual. Três
mecanismos consagrados são tipicamente invocados para semear campos magnéticos
cosmológicos:

- A **bateria de Biermann**, que gera campo magnético sempre que gradientes de
  densidade e de temperatura eletrônicas não são paralelos,
  \[
  \dfrac{\partial \mathbf{B}}{\partial t} \;=\; -\,\dfrac{c\, k_B}{e\, n_e}\,
  \nabla T_e \times \nabla n_e.
  \]
- Mecanismos **inflacionários** que quebram a invariância conforme do
  eletromagnetismo durante a época inflacionária.
- Mecanismos de **transição de fase** (eletrofraca ou QCD) que produzem um campo
  primordial estocástico via colisões de paredes de bolha.

As intensidades de semente previstas por esses canais cobrem muitas ordens de
grandeza; a amplificação subsequente pelo dínamo galáctico (§7.3) os eleva até
a banda μG amostrada pelo pipeline atual. Estabelecer qual canal de semente é
consistente com o espectro e o comprimento de coerência do campo reconstruído
pelo calibrador do modelo E-G fecharia uma das maiores lacunas conceituais em
aberto.

Referências: Subramanian (2016); Durrer & Neronov (2013); Widrow (2002).

### 7.3 Self-Excited Dynamo

O uso de amplitudes de Stokes do Planck como sonda de **estado estacionário**
do campo magnético do disco só se justifica se o campo for sustentado contra
difusão ôhmica em escalas de tempo muito menores que a idade do disco. A
equação relevante é a equação de indução magnetohidrodinâmica,

\[
\dfrac{\partial \mathbf{B}}{\partial t} \;=\; \nabla \times (\mathbf{v} \times \mathbf{B})
\;+\; \eta\, \nabla^{2} \mathbf{B},
\]

cujo primeiro termo advecta e estica o campo com o fluxo de velocidade `v`
enquanto o segundo termo o amortece com difusividade magnética `η`. Em um
disco com rotação diferencial, isso admite um dínamo `αΩ` auto-excitado: a
rotação diferencial `Ω` gera a componente toroidal, enquanto a turbulência
helicoidal de pequena escala (o efeito-`α`) regenera a poloidal. Valores
típicos de disco galáctico `τ_amp ∼ 10⁸`–`10⁹ yr ≪ τ_ohm` certificam que o
campo na escala μG está em equilíbrio dinâmico com a curva de rotação. Dentro
do modelo E-G, esse é o gancho de auto-consistência que permite que o mesmo
perfil `V_obs(r)` tanto *entre* no gap de aceleração (lado esquerdo) quanto
*sustente* o campo que o fecha (lado direito); um teste quantitativo desse
acoplamento — por exemplo, correlacionando `pearson_correlation` com proxies
do cisalhamento `Ω(r)` — é um próximo passo natural.

Referências: Beck (2015); Brandenburg & Subramanian (2005);
Ruzmaikin, Shukurov & Sokoloff (1988).

---

## 8. Referências

- Alfvén, H., 1970, *Cosmology in the Plasma Universe* (Palestra Nobel).
- Beck, R., 2015, *A&A Rev.*, 24, 4 (campos magnéticos galácticos).
- Bosma, A., 1981, *AJ*, 86, 1825.
- Brandenburg, A., Subramanian, K., 2005, *Phys. Rep.*, 417, 1 (campos
  magnéticos astrofísicos e o dínamo αΩ).
- Cordes, J. M., Wasserman, I., Hessels, J. W. T., et al., 2017, *ApJ*, 842, 35
  (plasma lensing).
- Durrer, R., Neronov, A., 2013, *A&A Rev.*, 21, 62 (campos magnéticos
  cosmológicos).
- Górski, K. M., Hivon, E., Banday, A. J., et al., 2005, *ApJ*, 622, 759.
- Haffner, L. M., Dettmar, R.-J., Beckman, J. E., et al., 2009, *Rev. Mod.
  Phys.*, 81, 969.
- Jackson, J. D., 1999, *Classical Electrodynamics*, 3rd edn (Wiley).
- Lelli, F., McGaugh, S. S., Schombert, J. M., 2016, *AJ*, 152, 157
  (catálogo SPARC).
- McGaugh, S. S., Schombert, J. M., Bothun, G. D., de Blok, W. J. G., 2000,
  *ApJ*, 533, L99.
- McGaugh, S. S., Lelli, F., Schombert, J. M., 2016, *Phys. Rev. Lett.*,
  117, 201101 (Radial Acceleration Relation).
- Navarro, J. F., Frenk, C. S., White, S. D. M., 1996, *ApJ*, 462, 563.
- Pen, U.-L., King, L., 2012, *MNRAS*, 421, L132 (lente gravitacional de
  plasma).
- Peratt, A. L., 1986, *IEEE Trans. Plasma Sci.*, 14, 763.
- Planck Collaboration, 2016, *A&A*, 594, A1; A&A 596, A105.
- Reynolds, R. J., 1990, *ApJ*, 349, L17.
- Rubin, V. C., Ford, W. K., 1970, *ApJ*, 159, 379.
- Ruzmaikin, A. A., Shukurov, A. M., Sokoloff, D. D., 1988, *Magnetic Fields
  of Galaxies* (Kluwer).
- Subramanian, K., 2016, *Rep. Prog. Phys.*, 79, 076901 (origem dos campos
  magnéticos cósmicos).
- Vedantham, H. K., 2019, *MNRAS*, 488, 4757 (refração cromática por meios
  ionizados).
- Widrow, L. M., 2002, *Rev. Mod. Phys.*, 74, 775 (origem dos campos
  magnéticos galácticos e extragalácticos).
