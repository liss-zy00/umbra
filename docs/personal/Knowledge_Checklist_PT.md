# Knowledge Checklist — Modelo E-G (uso pessoal)

> **Aviso.** Este arquivo é privado (`docs/personal/` está no `.gitignore`).
> É o guia de estudo da Leticia para sustentar o modelo E-G em qualquer
> debate técnico — formal ou informal. Não tem o rigor de documento
> público; pode (e deve) ser editado livremente conforme a compreensão
> evolui.

---

## A. Conceitos a dominar (em ordem crescente de profundidade)

Cada conceito vem com (i) um teste de "explique em uma frase" e (ii) o
ponteiro para onde a teoria formal está descrita.

1. **Curva de rotação galáctica.** Plot `V(r)` para uma galáxia espiral
   ou anã. O fato observacional é que `V(r)` fica plana para `r` grande,
   contrariando a previsão Kepleriana `V ∝ r^{-1/2}`. Teste: explique
   em uma frase por que isso é um problema para a gravidade Newtoniana
   só com matéria visível.
   *Referência:* [`../theory/EG_MODEL_THEORY_PT.md`](../theory/EG_MODEL_THEORY_PT.md) §1.

2. **Aceleração centrípeta vs aceleração bariônica.**
   `a_obs = V_obs²/r`, `a_bar = V_bar²/r`, e o **gap** `Δa = a_obs − a_bar`.
   Teste: dada uma linha do SPARC com `V_obs`, `V_gas`, `V_disk`, `V_bulge`
   e `r`, calcule `Δa` mentalmente em ordem de magnitude.
   *Referência:* §1, equação (1.1).

3. **Matéria escura como variável de fechamento.** O halo de matéria
   escura é, no ΛCDM, *a solução* para o gap. Teste: nomeie três tipos
   de candidato (WIMP, áxion, MACHO) e diga por que nenhum foi detectado
   diretamente.

4. **Plasma e meio interestelar (ISM).** Quatro fases do ISM: molecular,
   neutra atômica (HI), ionizada morna (WIM) e ionizada quente (HIM).
   Teste: nomeie temperatura e densidade típicas de cada fase.

5. **WIM (Warm Ionised Medium).** O componente que importa para o E-G:
   `T ~ 10⁴ K`, `ρ ~ 10⁻²⁷ kg/m³` no halo. Teste: explique por que ele
   é invisível para a linha de 21 cm mas relevante para a força de Lorentz.

6. **Força de Lorentz microscópica vs volumétrica.**
   `F = q(E + v×B)` é microscópica; `f = J×B = (∇×B)×B/μ₀` é volumétrica.
   Teste: derive a segunda forma a partir da primeira somando sobre as
   cargas de uma parcela de fluido quase-neutra.
   *Referência:* §2.

7. **Decomposição da força de Lorentz: pressão + tensão.**
   `J×B/μ₀ = −∇(B²/2μ₀) + (B·∇)B/μ₀`. Teste: explique fisicamente o que
   cada termo faz. Por que o E-G atual descarta a tensão?
   *Referência:* §2.3 e §5.4 das limitações.

8. **Densidade de energia magnética como pressão.**
   `u_B = B²/(2μ₀)` tem unidade de pressão (J/m³ = Pa). Teste: para
   `B = 1 µG`, qual é `u_B` em SI? Compare com a pressão térmica do WIM.

9. **Aceleração de pressão magnética.**
   `a_B ≈ −∇u_B / ρ_ion`. Teste: explique por que o denominador é
   `ρ_ion` e não `ρ_total` ou `ρ_HI`.

10. **MHD ideal e congelamento do fluxo.** O resultado de Alfvén:
    `E + v×B = 0` em fluido condutor perfeito; linhas de campo congeladas
    no plasma. Teste: dê uma condição (escala, condutividade) sob a qual
    MHD ideal falha em galáxias.

11. **Plasma-β.** `β = 2μ₀ n k_B T / B² = P_thermal / P_magnetic`. Teste:
    `β < 1` versus `β > 1` — qual é o regime onde o *ansatz* E-G faz
    sentido? Por quê?

12. **Velocidade de Alfvén.** `v_A = B / √(μ₀ ρ_ion)`. Teste: para
    `B = 1 µG` e `ρ_ion = 10⁻²⁷ kg/m³`, calcule `v_A` em km/s. (Resposta:
    centenas a milhares de km/s, dependendo da densidade.)

13. **Parâmetros de Stokes Q e U.** Polarização linear: `S = √(Q² + U²)`
    é a amplitude, `ψ = atan2(U, Q)/2` é o ângulo. Teste: por que o
    Planck mede `Q/U` em vez de medir `B` diretamente?

14. **A ponte dimensional Stokes → Tesla.**
    `B = √(2μ₀ u_ref) · S/S_ref`, com `u_ref = u_CMB` como **âncora
    dimensional** (não afirmação física). Teste: o que muda em `a_B` se
    `u_ref` for trocada por equipartição turbulenta?

15. **Pixelização HEALPix.** `NSIDE`, `NPIX = 12·NSIDE²`. Teste: dado
    `NSIDE = 2048`, qual é a resolução angular típica em arcminutos?
    Como isso vira `Δx_phys` para uma galáxia a 10 Mpc?

16. **Gradiente como diferença finita escalar.**
    `Δu_B = u_B^vizinhos − u_B^centro`; `∇u_B · ê ≈ Δu_B / Δx`. Teste:
    por que essa diferença finita é uma cota superior do gradiente radial
    *no plano do disco*?

17. **`scenario_gradient_scale_kpc` e `scenario_ionization_fraction`.**
    Os dois botões do grid. Teste: para o cenário (0,1 kpc, 0,01), qual
    é a `ρ_ion` esperada em uma anã rica em gás?

18. **`theoretical_gap_coverage_ratio`.** A métrica de comparação. Teste:
    valores `~ 1`, `~ 0,01` e `~ 100` autorizam que conclusão?

19. **Correlação de Pearson.** Testa a forma, não a magnitude. Teste:
    dê um exemplo onde a razão é 1 mas a correlação é negativa — o que
    significa fisicamente?

20. **`audit_flags`.** Lista de códigos bloqueantes. Teste: nomeie os
    sete códigos bloqueantes e o único informacional `warn_*`.
    *Referência:* [`../results_interpretation/Scientific_Evidence_Narrative_PT.md`](../results_interpretation/Scientific_Evidence_Narrative_PT.md) §4.

21. **Controles negativos NFW.** Curvas sintéticas que o pipeline
    **não** deve confirmar. Teste: explique por que isso fecha apenas
    *uma* classe de falsos positivos.

22. **RAR e BTFR.** Relação radial de aceleração (McGaugh+ 2016) e
    relação Tully-Fisher bariônica (McGaugh+ 2000). Teste: por que um
    modelo que pretende substituir matéria escura precisa satisfazer
    ambas simultaneamente?

---

## B. Equações-chave (memorizar)

| # | Equação | Significado |
| - | ------- | ----------- |
| (1) | `a_obs(r) = V_obs²(r) / r` | aceleração centrípeta observada |
| (2) | `V_bar² = V_gas² + V_disk² + V_bulge²` | composição bariônica SPARC |
| (3) | `gap(r) = a_obs(r) − a_bar(r)` | lacuna que matéria escura preenche |
| (4) | `F = q(E + v × B)` | Lorentz microscópica |
| (5) | `f = J × B = (∇×B) × B / μ₀` | Lorentz volumétrica (MHD) |
| (6) | `J × B / μ₀ = −∇(B²/2μ₀) + (B·∇)B/μ₀` | decomposição pressão + tensão |
| (7) | `u_B = B² / (2μ₀)` | densidade de energia magnética |
| (8) | `a_B ≈ −∇u_B / ρ_ion` | equação operacional central |
| (9) | `ρ_ion = f_ion · ρ_gas` | acoplamento WIM |
| (10) | `B = √(2μ₀ u_ref) · S / S_ref` | ponte dimensional Stokes → Tesla |
| (11) | `β = 2μ₀ n k_B T / B²` | plasma-β (termal vs magnético) |
| (12) | `v_A = B / √(μ₀ ρ_ion)` | velocidade de Alfvén |
| (13) | `θ_coverage = \|⟨a_B⟩_edge\| / \|⟨gap⟩_edge\|` | métrica principal |

---

## C. Antíteses e respostas (perguntas de revisor)

### C.1 "Lentes gravitacionais provam que a matéria escura existe."

**Resposta defendível:** "O E-G não nega a existência de massa
gravitacional adicional onde quer que ela seja necessária. O que estou
afirmando é que a **dinâmica de rotação em discos galácticos** pode ser
explicada — em parte ou totalmente — sem invocar essa massa, via pressão
MHD. Lentes em aglomerados (como o Bullet Cluster) e dinâmica de
rotação são problemas físicos *diferentes*; o E-G é silencioso sobre o
primeiro por design (veja §8 do `Simplified_Theory_PT.md`)."

Plano B (mais especulativo): mencionar refração de plasma como mecanismo
secundário de desvio de luz, sem fazer disso o argumento principal.

### C.2 "ΛCDM explica a CMB perfeitamente."

**Resposta defendível:** "Concordo. O modelo E-G é restrito ao problema
das curvas de rotação. A CMB envolve a física do universo primordial, que
está fora do escopo desta proposta. Magnetogênese primordial (campos B
nos primeiros segundos) é um tópico em aberto na literatura, mas eu não
estou alegando que o E-G explica a CMB."

### C.3 "Campos B galácticos são fracos demais (microGauss) para produzir aceleração da ordem de `10^{-11} m/s²`."

**Resposta defendível:** "A magnitude da aceleração depende da razão
entre o gradiente de pressão magnética e a densidade do plasma responsivo.
Em `ρ_ion ~ 10^{-27} kg/m³` (WIM), um campo de microGauss produz
aceleração na ordem correta. O pipeline calcula isso shell por shell e
reporta o resultado dentro do `audit_flags`."

### C.4 "A resolução do Planck (arcminutos) é grande demais para amostrar gradientes em parsecs."

**Resposta defendível:** "Esse é um problema real, listado como
limitação §5.1 do `EG_MODEL_THEORY_PT.md`. É exatamente por isso que o
pipeline reporta **duas** colunas paralelas: uma usando `Δx` físico
(distância × ângulo HEALPix) e outra usando `scenario_gradient_scale_kpc`
como botão exploratório. A comparação entre as duas quantifica o
desencontro de escala."

### C.5 "A âncora CMB é arbitrária."

**Resposta defendível:** "Concordo que é fenomenológica. Ela serve para
dar **dimensão** correta a `B`; a magnitude absoluta depende dessa escolha.
A análise de sensibilidade de âncora (`validation_anchor_sensitivity_sweep.csv`, ligada ao
item 4 da Fase 1 no código) quantifica o deslocamento entre âncoras; relatar explicitamente
esse *sweep* em qualquer narrativa pública." "Você está descartando o termo de tensão magnética. Isso pode mudar tudo."

**Resposta defendível:** "Sim, está descartado neste *ansatz*
(`EG_MODEL_THEORY_PT.md` §2.3 e §5.4). O argumento é que o termo de
pressão é o **maior em ordem de magnitude** quando o campo é
predominantemente turbulento; em campos altamente ordenados a tensão
pode competir. O roadmap (item 5) adiciona uma correção de tensão *toy*
em geometria espiral logarítmica para limitar a ordem de magnitude do
termo faltante."

### C.7 "Por que o WIM e não o gás HI denso?"

**Resposta defendível:** "Porque a Força de Lorentz, na forma volumétrica,
atua sobre partículas carregadas. O gás HI é majoritariamente neutro e
responde ao magnetismo apenas indiretamente, via colisões com o plasma.
A `f_ion ∈ (0, 1]` no pipeline encapsula essa fração ionizada, e o grid
search testa três valores plausíveis."

### C.8 "Você está superajustando ao varrer o grid."

**Resposta defendível:** "Por isso o resultado primário é o **heatmap
populacional**, não a galáxia individual. Se um único cenário concentrasse
todas as galáxias confiáveis em torno de `theoretical_gap_coverage_ratio ≈ 1`,
isso seria uma escala física emergente. Se a distribuição estiver
pulverizada por todos os nove cenários, é ajuste manual disfarçado e
o resultado é não-conclusivo."

---

## D. Perguntas-armadilha de revisor (mais sutis)

1. **"O `theoretical_gap_coverage_ratio` é avaliado no `edge` do disco.
   Como `edge` é definido operacionalmente, e ele é robusto a galáxias
   truncadas?"** → resposta envolve falar de `edge_radial_mean` em
   [`src/physics/diagnostics/statistics.py`](../../src/physics/diagnostics/statistics.py)
   e da fração das cascas externas usada.

2. **"A inversão de cascas (3.3)–(3.4) assume simetria esférica para gás
   que está em um disco. Quanto isso polui `ρ_ion`?"** → resposta: é
   exatamente a limitação §5.2; o cross-check `MHI_1e9_solar` (item 3 do
   roadmap) é o que vai responder em produção.

3. **"Por que `max(M_gas^A, M_gas^B)` em vez de soma ou média?"** →
   resposta: porque `M_gas^B` subtrai o potencial estelar; se as estrelas
   já explicam a dinâmica naquela casca, `M_gas^B → 0` e o `max` evita
   inflar o reservatório inercial.

4. **"A normalização `S_ref` como mediana do mapa torna a calibração
   auto-referencial. Por que isso é uma vantagem, e não um problema?"** →
   resposta: porque elimina dependência da escala absoluta do mapa
   (sensível ao processamento Commander vs GNILC), preservando apenas o
   contraste.

5. **"Se `audit_flags` é vazio mas a correlação de Pearson é negativa,
   o que isso significa?"** → resposta: `is_reliable` cobre apenas
   critérios de qualidade, não a hipótese física. Uma correlação negativa
   com razão de cobertura próxima de 1 sugere que `a_B` e `gap` têm
   magnitudes iguais por coincidência mas formas radiais diferentes —
   evidência *contra* o modelo nessa galáxia.

6. **"Por que a grade tem só nove cenários e não, digamos, vinte e sete?"**
   → resposta: porque os botões legados serão *aposentados* (item 9 do
   roadmap) quando deprojeção em disco + propagação Bayesiana estiverem
   prontas. A grade atual existe como exploração, não como resultado
   científico final.

---

## E. Cronograma sugerido de aprofundamento

| Mês | Tópico | Recurso |
| --- | ------ | ------- |
| 1 | Equações de Maxwell, vetorial, MHD ideal | Jackson cap. 6–8; "Plasma Physics" de Chen |
| 2 | Física de plasma interestelar, WIM, Faraday rotation | Draine "Physics of the Interstellar and Intergalactic Medium" |
| 3 | Curvas de rotação, SPARC, RAR, BTFR | Lelli+ 2016; McGaugh+ 2016 |
| 4 | Polarização de poeira, HEALPix, Planck | Planck Coll. 2016 *A&A* 594 A1; healpy docs |
| 5 | Aprofundar §5 das limitações + guias de resultados | `../theory/EG_MODEL_THEORY_PT.md` e `../results_interpretation/` |

---

## F. O que ler **antes** de qualquer apresentação pública

1. `../theory/EG_MODEL_THEORY_PT.md` (foundational, atemporal).
2. `../educational/Simplified_Theory_PT.md` (para alinhar o nível de
   complexidade ao público).
3. `../results_interpretation/Scientific_Evidence_Narrative_PT.md`
   (para responder qualquer pergunta sobre os números do CSV).
4. Este arquivo (Knowledge_Checklist_PT) — pelo menos a §C (antíteses).

Se a apresentação for acadêmica formal, abrir também `../academic/Draft_Paper_EG_Model_PT.md`.
