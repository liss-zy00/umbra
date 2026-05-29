# Avaliação de ponta a ponta — Modelo E-G, código e resultados

> **Tipo de documento.** Nota pessoal em `docs/personal/`. Sintetiza a proposta teórica
> em `docs/theory/`, a arquitectura em `src/`, e a leitura dos artefactos em `results/`
> e `docs/results_interpretation/`, à luz do estado do repositório nesta sessão (maio de 2026).
> Não substitui o manifesto atemporal em `docs/theory/` nem números congelados noutros
> relatórios; quando reexecutar o pipeline, cruze sempre com `validation_consistency_audit.md`
> regenerado.

---

## 1. Resumo executivo

O projeto implementa um **teste falsificável**, não uma “prova” de matéria escura alternativa:
pergunta se um **gradiente de pressão magnética** (MHD ideal, ancorado a polarização Planck
e cinemática SPARC) consegue, em ordem de grandeza e tendência radial, acompanhar o
**déficit cinemático** \(a_{\mathrm{obs}}^2/r - a_{\mathrm{bar}}^2/r\) onde a curva de
rotação exige aceleração além do bariónico newtoniano. A documentação canónica deixa explícito
que um resultado “sim” **não refuta** matéria escura fria; um “não” refutaria **este** ansatz
de pressão com os dados e a escala de resolução usados.

Sobre os **melhores dados** que o pipeline consegue hoje: a sub-amostra com
`validation_is_reliable` (28 de 158 entradas no grid, de 175 rotmods processados) é
**útil para exploração e auditoria interna**, mas **não é**, por si só, evidência definitiva
de física nova. Motivos centrais: grelha exploratória \((\Delta x, f_{\mathrm{ion}})\),
âncora dimensional CMB para Stokes→Tesla (não física do disco), resolução Planck vs disco,
e fração elevada de avisos (`warn_no_pa`, `warn_pure_pressure`, BTFR com 25% de outliers
|residual| > 0.5 dex onde aplicável).

**Coerência interna:** a tradução matemática básica (Lorentz → pressão + tensão opcional)
em `docs/theory/EG_MODEL_THEORY_EN.md` alinha-se com o que `PhysicsEngine` calcula; o código
é modular (provider → motor → diagnósticos → exportação). Os **maiores riscos** não são
“erros de sinal” triviais, mas **decisões de modelagem** documentadas (mistura de escala
céu–disco no gradiente, calibrador global, tensão em espiral logarítmica como *toy*).

**Impacto no meio da física atual:** em termos realistas, isto posiciona-se como
**astrofísica galáctica / MHD observacional** com ligação ao puzzle das curvas de rotação —
não como revisão de gravidade fundamental ou solução única ao problema da matéria escura.

**Publicação:** um artigo bem feito poderia contribuir como **método + conjunto de
diagnósticos + limitações honestas**; sem dados de maior resolução e sem reduzir graus de
liberdade exploratórios, o caminho típico seria **periferia de alto impacto** ou
**journals especializados**, com revisão exigente sobre calibração e degenerescências.

---

## 2. O que a teoria propõe (leitura de `docs/theory/`)

Ficheiros centrais: `EG_MODEL_THEORY_EN.md` / `EG_MODEL_THEORY_PT.md`.

- **Hipótese de estados estendidos:** em plasmas magnetizados e estruturas filamentares,
  stresses magnéticos podem decair mais lentamente que campos gravitacionais tipo ponto;
  o “teste E-G” pergunta se o **canal MHD** (via \(u_B = B^2/(2\mu_0)\) e \(\nabla u_B\))
  fecha **parte** do *gap* cinemático às distâncias relevantes.
- **Tradução para o código:** o motor usa a forma **pressão isotrópica** do volume de
  Lorentz em MHD ideal; a tensão \((\mathbf{B}\cdot\nabla)\mathbf{B}\) aparece como
  extensão *toy* (espiral logarítmica) onde a morfologia o permite — ver comentários em
  `src/physics/engine.py`.
- **Postura científica:** falhas de *match* (razões longe de 1) são resultados válidos;
  incertezas de catálogo propagam-se por Monte Carlo nas colunas `p16`/`p84` onde
  implementado.

Isto está alinhado com a filosofia do próprio repositório (falsificabilidade, sem *ajuste*
opaco para forçar *match*).

---

## 3. Da fonte de dados ao resultado (arquitectura)

### 3.1 Ingestão

| Etapa | Módulo / pasta | Função |
|--------|----------------|--------|
| Cinemática | `src/providers/sparc.py` | Curvas SPARC, \(V_{\mathrm{bar}}\) a partir de gás+disco+bulge, distância do cabeçalho, incerteza de inclinação para coluna de erro sistemático em \(V\). |
| Morfologia | `SparcMorphologyProvider` + JSON Lelli+2016 | Classes de Hubble agregadas, LSB fotométrico para regras de relatório. |
| Céu / polarização | `src/providers/planck.py` (e suplementos NED) | Stokes, anel HEALPix, provisão de PA quando existir. |
| Constantes | `src/config/thresholds.py`, `src/physics/defaults.py` | Reduz *hardcoding* físico no corpo das funções. |

### 3.2 Motor físico

| Componente | Ficheiro | Papel |
|-------------|----------|--------|
| *Gap* | `src/physics/kinematics.py` | \(a_c = V_{\mathrm{obs}}^2/r\), \(a_{\mathrm{bar}} = V_{\mathrm{bar}}^2/r\), `acceleration_gap`; RAR McGaugh+2016 como diagnóstico paralelo. |
| Campo | `src/physics/magnetic_field_calibration.py` | \(B \propto \sqrt{2\mu_0 u_{\mathrm{ref}}}\,(S/S_{\mathrm{ref}})\) com \(u_{\mathrm{ref}}\) tipo CMB como **âncora dimensional**. |
| Densidade | `src/physics/gas_density.py` | \(\rho_{\mathrm{gas}}\) a partir das componentes SPARC; \(f_{\mathrm{ion}}\) multiplica para \(\rho_{\mathrm{ion}}\). |
| Gradiente | `src/physics/magnetic_pressure_gradient.py` | LSQ no plano do disco quando PA e \(i\) existem; fallback axisimétrico com `warn_no_pa`. |
| Orquestração | `src/physics/engine.py` | Compõe dependências, cenários, colunas paralelas (`magnetic_pressure_acceleration_physical_dx`, Lorentz total onde aplicável). |

### 3.3 Pipeline e resultados

| Etapa | Saída típica |
|--------|----------------|
| Discovery | `results/step_discovery/discovery_*.csv` — coorte, `required_density_for_match`, factores de resolução, `discovery_is_reliable` / elite. |
| Validação (grelha 3×3) | `results/grid_results/scenario_*.csv`, `validation_all.csv`. |
| Fiabilidade | `validation_reliable.csv` / `validation_not_reliable.csv`, `validation_is_reliable`. |
| Auditoria populacional | `results/step_validation/validation_consistency_audit.md`. |
| “Evidência directa” (subconjunto) | `validation_scientific_evidence_report.csv` — *supreme rows* (fiável + *gate* de inclinação), filtro Im **ou** LSB fotométrico, e `theoretical_gap_coverage_ratio` na janela **0,2–5,0** (`build_direct_evidence_rows`). |
| Por galáxia | `results/galaxy_analysis/.../*_analysis.csv`, figuras, `*_summary.json`. |

**Nota importante para leitura dos CSV:** o ficheiro `validation_scientific_evidence_report.csv`
é **sub-amostrado** (morfologia + janela de razão; ver `docs/results_interpretation/`). Para estatística da
coorte fiável completa use `validation_reliable.csv` + `validation_all.csv` e o *audit*.

---

## 4. O que os resultados dizem (evidências extraíveis)

Os números abaixo reflectem o **último** `validation_consistency_audit.md` presente no
repositório; se alterar dados ou código, regenere o pipeline.

### 4.1 Tamanhos de amostra

- **175** rotmods processados (`discovery_all.csv`).
- **158** com `discovery_is_reliable` (17 bloqueadas antes do grid).
- **158** entradas no grid de validação; **28** com `validation_is_reliable` (~17.7%).

Isto mostra que a **maioria** das galáxias falha *gates* de fiabilidade (fração de cascas
fiáveis, domínio estelar, inclinação, etc. — ver frequências de `audit_flags` no *audit*).

### 4.2 Métricas populacionais na sub-amostra fiável (N = 28)

- **RAR:** mediana de \(\|\mathrm{residual}\|\) ≈ **0.162 dex**; MAD ≈ **0.109 dex** — compatível
  com a literatura RAR a nível grosso, mas **não** isola o mecanismo magnético (o RAR já
  está construído sobre \(g_{\mathrm{obs}}\)).
- **BTFR:** com residual finito (N = 16), **25%** têm \(|\mathrm{residual}| > 0.5\) dex —
  sinal de que **parte** da coorte fiável ainda é atípica na escala bariónica clássica.
- **`theoretical_gap_coverage_ratio`:** medianas por grupo morfológico no *audit* variam
  entre ~0.18 e ~0.79 nos grupos com poucos elementos; **Im** mostra média enorme (~105)
  com desvio padrão ~263 por **outliers** — ou seja, a cauda domina a média; a mediana
  (~0.257) diz que, para muitos Im, o canal de **pressão sozinha** fica **abaixo** de
  fechar o *gap* por um factor típico ~4× (ordem 0.26 em \(|\langle a_B\rangle|/|\langle\mathrm{gap}\rangle|\)).
- **Correlação Pearson** entre `magnetic_pressure_acceleration` e `acceleration_gap` é
  **informativa por galáxia**; valores negativos (ex.: NGC3741, UGC05764 no extracto do
  CSV de evidência) **enfraquecem** a narrativa de “forma acompanha déficit” mesmo quando
  a razão de magnitudes é moderada.

### 4.3 Discovery vs WIM

- `required_density_for_match` na coorte *passed* situa-se muito acima do *anchor* WIM
  típico — a leitura honesta é que **beam smearing** mediano (~factor 2–3) **não** explica sozinho
  **três ordens de magnitude**; restam calibração, refinamentos físicos (deprojecção, tensão onde activa), ou densidades ionizadas reais
  mais altas que o *floor* WIM (ver Secção 5 de `EG_MODEL_THEORY_PT.md` e guias de interpretação).

### 4.4 Sinais de sistema e auditoria

- **`warn_no_pa`:** ~82% das galáxias no *audit* carregam este aviso — gradiente no plano
  do céu com correcção cos\(i\) em vez de geometria de braço espiral completa.
- **`fail_stellar` / `fail_reliable_fraction`:** muito frequentes — o teste **exige** regiões
  onde o *gap* e o proxy magnético são simultaneamente bem definidos; discos dominados por
  estrelas ou com poucas cascas úteis removem poder estatístico.
- **Latitude galáctica:** correlações fracas a moderadas conforme o subconjunto; o *audit*
  separa \(|b| < 30^\circ\) vs \(\geq 30^\circ\) com N pequeno em baixa latitude — interpretar
  com cautela (medianas próximas, médias dominadas por outliers).

### 4.5 Controles de ouro (DDO154, NGC3741)

Ambos aparecem como **Im** fiáveis com `warn_pure_pressure` — a documentação consolidada
já nota: para estudar **tensão** de espiral de forma limpa, precisa de espirais onde o *toy*
de tensão não está desactivado por morfologia. Continuam úteis como **regressão** do
pipeline e para bandas MC estreitas quando os erros de catálogo são pequenos.

---

## 5. Os “melhores dados” são bons?

**Sim, com o significado correcto:** são o melhor que o desenho actual permite — coorte
SPARC de referência, Planck GNILC, propagação de erros onde implementada, *flags* de
proveniência e relatório de auditoria automático. Isso é **nível de trabalho sério de
pré-publicação** em termos de engenharia e transparência.

**Não** no sentido de “fechamento definitivo do *gap* com polarização all-sky”: ainda são
**ruins** para uma reivindicação forte porque (i) a grelha \((\Delta x, f_{\mathrm{ion}})\)
funciona como **knobs** que podem absorver erro sistemático; (ii) a escala física do
gradiente no disco não coincide com a separação nativa do mapa; (iii) a calibragem
Stokes→\(B\) não está amarrada a síncrotron local; (iv) tensão e pressão competem de
forma dependente da geometria assumida.

Em uma frase: **bons para delimitar o ansatz e para falsificação interna; insuficientes,
sem dados de alta resolução e calibração local, para reescrever o consenso cosmológico.**

---

## 6. Lógica de cálculo e física — coerência e pontos fracos

### 6.1 O que é coerente

- Uso consistente de \(a = V^2/r\) e decomposição SPARC para \(V_{\mathrm{bar}}\).
- Passagem de Lorentz em MHD ideal para \(-\nabla P_B\) documentada e implementada.
- Separação SOLID (motor vs providers vs exportação) e ausência de “ajuste manual” por
  galáxia nos resultados agregados.

### 6.2 Onde o revisor mais forte vai bater

1. **Interpretação do `theoretical_gap_coverage_ratio`:** razão de médias (ou agregados)
   no bordo do disco pode **mascarar** cancelamentos radiais; o `pearson_correlation` tenta
   corrigir isso parcialmente — deve estar sempre **junto** na narrativa.
2. **Âncora CMB:** dimensionalmente limpa, fisicamente **não** justificada como amplitude
   do campo de disco; o *sweep* `validation_anchor_sensitivity_sweep.csv` deve ser
   citação obrigatória em qualquer texto público.
3. **Grelha de cenários:** nove combinações por galáxia implicam **multiplicidade de testes**;
   sem correcção explícita ou *hold-out*, o *p*-valor informal de “encontrei um cenário
   que encaixa” é optimista.
4. **Tensão logarítmica:** útil como sensibilidade; não confundir com reconstrução 3D do
   campo vectorial real.

### 6.3 Possíveis erros lógicos a evitar na comunicação

- Não dizer “o modelo prevê a RAR” só porque os resíduos RAR são pequenos — a RAR já está
  construída sobre observáveis que entram no *gap*.
- Não misturar **28 fiáveis** com **4 linhas** no CSV de “evidência científica” sem explicar
  o filtro Im/LSB.

---

## 7. Impacto no “meio” da física actual

- **Cosmologia e matéria escura (CMB, lensing, cluster):** este trabalho **não** os
  confronta directamente; mesmo que o MHD fechasse o *gap* em algumas galáxias, **não**
  substitui halos em escala de aglomerados ou BAO.
- **Galáxias e baryons:** o impacto potencial é **reclassificar parte do “missing”
  cinemático** como **sistema acoplado plasma–campo–feedback** em discos; isso intersecta
  literatura sobre **fields, cosmic rays, turbulence**, mas ainda como **correcção
  sub-dominante ou regional**, não como revolução única.
- **MOND / emergência fenomenológica:** a RAR entra como **diagnóstico lateral**; o projecto
  não implementa uma modificação newtoniana explícita — qualquer paralelo é narrativo,
  não estrutural no código.

---

## 8. Se publicasse os estudos — cenário realista

**Caminho positivo plausível:** submissão para *Monthly Notices*, *A&A*, *ApJ* ou
*journals* de plasma astrofísico, com ênfase em **método reprodutível**, **dados públicos**
(SPARC, Planck), e **limitações** (resolução, âncora, multiplicidade de cenários). Aceite
provável: **revisão dura** pedindo JVLA/SKA, síncrotron, ou redução de graus de liberdade.

**Caminho negativo plausível:** rejeição por “modelo demasiado flexível” ou “calibração
arbitrária”, ou transferência para repositório / *data note* sem peer review de alto
  factor.

**Risco reputacional:** afirmações fortes sem dados de alta resolução são **mais** danosas
  do que um abstract conservador com *clear falsifiability*.

**Institucional / carreira:** depende do contexto (afiliação, coautores estabelecidos,
  divulgação prévia em arXiv com moderação). O código aberto e o *audit* automático são
  **pontos a favor** de integridade percebida.

---

## 9. Conclusão

O repositório cumpre bem o papel de **laboratório computacional disciplinado** para testar
um ansatz MHD específico contra curvas de rotação públicas. Os **resultados actuais**
mostram **compatibilidade RAR modesta** na sub-amostra fiável, **tensão com BTFR** em uma
fração material, e **cobertura fraca do *gap* pelo canal de pressão sozinha** para muitos
Im — tudo isso são **sinais científicos legítimos**, não falhas de implementação.

O passo seguinte honesto é **dados e geometria**: polarização rádio resolvida, melhor PA/distância (Gaia),
calibração por síncrotron, e relatórios que **sempre** acoplem âncora de energia e escala
de gradiente à conclusão.

---

## Referências internas rápidas

| Tópico | Onde |
|--------|------|
| Teoria | `docs/theory/EG_MODEL_THEORY_EN.md` |
| Leitura de colunas | `docs/results_interpretation/Data_Interpretation_Guide_PT.md` (sincronizado com o código em maio de 2026) |
| Narrativa estável | `docs/results_interpretation/Scientific_Evidence_Narrative_PT.md` (idem) |
| Rascunho académico | `docs/academic/Draft_Paper_EG_Model_PT.md` |
| Motor | `src/physics/engine.py` |
| Auditoria numérica | `results/step_validation/validation_consistency_audit.md` |
