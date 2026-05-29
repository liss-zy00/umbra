# Evolution Log — Relatório Interno da Refatoração

> **Aviso.** Este arquivo é privado (`docs/personal/` está no `.gitignore`).
> Ele **não é** um documento científico; é o relatório operacional da
> refatoração de documentação do projeto E-G: o que foi **movido**, o que
> foi **removido**, o que foi **refutado/substituído** durante a evolução
> da hipótese.
>
> Edite livremente. Em caso de conflito com a documentação canônica em
> `docs/theory/`, a documentação canônica prevalece.

---

## 1. Movido (refatoração estrutural)

| Origem | Destino | Motivo |
| ------ | ------- | ------ |
| `EG_MODEL_THEORY.md` (raiz) | `docs/theory/EG_MODEL_THEORY_EN.md` | Centralizar toda documentação técnica em `docs/`. |
| `EG_MODEL_THEORY.md` (raiz) | `docs/theory/EG_MODEL_THEORY_PT.md` (novo) | Versão bilíngue para uso interno e apresentações. |
| `MODEL_REFINEMENTS.md` (raiz) | `docs/personal/` (notas privadas de roadmap; não citadas na documentação pública do repositório) | Separar *backlog* operacional da teoria canónica em `docs/theory/`. |
| `MODEL_REFINEMENTS.md` (raiz) | `docs/personal/MODEL_REFINEMENTS_PT.md` (legado) | Tradução pessoal antiga; conteúdo absorvido pelo V1 onde aplicável. |
| `Condensed_theory.txt` (raiz, ~1357 linhas) | `docs/personal/Condensed_Theory.md` (Markdown; pasta `docs/personal/` no `.gitignore`) | Arquivo de arquivo pessoal; conteúdo técnico redistribuído — ver secção abaixo. |

Conteúdo redistribuído a partir de `Condensed_theory.txt`:

- Núcleo teórico em inglês (linhas 462–510 e 752–1248) → `docs/theory/EG_MODEL_THEORY_EN.md`.
- Resumos didáticos em PT (linhas 1–60, 263–294, 1250–1295) → base do
  `docs/educational/Simplified_Theory_PT.md`.
- "Arquitetura Híbrida do Cosmo" + abstract em EN (linhas 263–294 e
  462–510) → `docs/academic/Draft_Paper_EG_Model_*.md`.
- Glossário, pré-requisitos e antíteses (linhas 123–164, 707–750,
  1296–1323) → `docs/personal/Knowledge_Checklist_PT.md`.
- Interpretação de resultados (linhas 365–460) → `docs/results_interpretation/Scientific_Evidence_Narrative_*.md`.

---

## 2. Removido da exposição pública (sem perda de informação)

Cada item abaixo foi extraído do *raw* `Condensed_theory.txt` mas
**não** está presente em nenhum documento de `docs/theory/`,
`docs/educational/`, `docs/academic/` ou `docs/results_interpretation/`.
A razão de remoção é registrada para que a decisão possa ser revertida
se algum desses pontos voltar a fazer sentido.

| Item | Onde estava no raw | Por que foi removido |
| ---- | ------------------- | -------------------- |
| Inspiração religiosa/mítica (Saturno bíblico, Caibalion, petróglifos) | linhas 49, 259, 715–717 | Não falsificável; cria dependência argumentativa em ideias não-científicas que enfraquecem a exposição técnica. |
| Citação "Instabilidades de Peratt como descargas testemunhadas pela humanidade pré-histórica" | linha 49 | Especulação arqueológica fora do escopo do modelo. |
| "Universo elétrico" como rótulo da hipótese | linhas 85, 259, 583 | Carrega bagagem pseudocientífica; o termo correto e neutro é "MHD". |
| Frases motivacionais sobre o feito pessoal da Leticia (~30+ ocorrências) | linhas 374–425, 510–534, 1338–1354, etc. | Pertinentes para a autora, irrelevantes para o leitor científico. Mantidas apenas neste arquivo. |
| Estimativas isoladas de cobertura por galáxia (e.g. "DDO154 49%", "ESO444 97.52%", "match de 132.47%") | linhas 372, 443–447, 454 | Valores de execução não devem viver em documentação foundational/atemporal. Permanecem nos artefatos de runtime sob `results/`. |
| Análise sociológica/política da publicação ("Mondianos", "burburinho na comunidade") | linhas 1325–1354 | Especulação sobre recepção, não conteúdo técnico. |

---

## 3. Refutado / substituído durante a evolução da hipótese

Cada ideia abaixo foi explicitamente abandonada em favor de uma versão
mais defensável. Esta seção é o **histórico intelectual** do projeto:
serve para que apresentações futuras saibam quais armadilhas evitar e
por quê.

### 3.1 "Finite State Machine" literal → regimes dinâmicos co-atuantes

**Ideia inicial.** O Universo opera como uma FSM (Máquina de Estados
Finitos), com um *gatilho* de densidade de corrente crítica `J_critical`
que comuta a dominância entre gravidade e Lorentz.

**Por que foi substituída.** A formulação FSM sugere transição brusca
("interruptor"), o que não tem base na MHD. A natureza descreve uma
**superposição contínua**: ambas as forças coexistem sempre e seu peso
relativo muda gradualmente com a geometria (decaimento `1/r²` vs `1/r`).

**Forma atual.** *Regime A* (compacto/bariônico) e *Regime B*
(estendido/magnetizado) como **decomposição contínua** do orçamento de
força radial, sem gatilho discreto. Ver
`docs/theory/EG_MODEL_THEORY_PT.md` §1 e
`docs/academic/Draft_Paper_EG_Model_PT.md` §2.

---

### 3.2 Cálculo direto de `J = (∇×B)/μ₀` pixel-a-pixel → gradiente de pressão `−∇(B²/2μ₀)`

**Ideia inicial.** Tomar os mapas Stokes do Planck, reconstruir `B` em
3D e calcular o rotacional `∇×B` pixel-a-pixel para obter `J`; depois
aplicar `F = J×B` como Força de Lorentz volumétrica.

**Por que foi substituída.** A resolução angular do Planck é de
arcminutos. O rotacional discreto numa grade dessa resolução *dilui* o
sinal de corrente em filamentos de parsecs por ordens de magnitude. O
resultado bate na medida do ruído, não da física.

**Forma atual.** Usar a identidade vetorial
`J×B/μ₀ = −∇(B²/2μ₀) + (B·∇)B/μ₀` e calcular apenas o gradiente escalar
de `u_B`, que sobrevive à resolução do Planck. Ver §2.3 da teoria.

---

### 3.3 Força de Lorentz microscópica `F = q(v × B)` sobre estrelas → forma volumétrica MHD `J × B`

**Ideia inicial.** Tratar cada estrela como uma partícula carregada
sentindo `F = q(v×B)` diretamente.

**Por que foi substituída.** Estrelas são **eletricamente neutras** em
escalas macroscópicas. Lorentz microscópica não as afeta. O que sofre
Lorentz é o **plasma** que permeia o disco; estrelas são arrastadas por
acoplamento gravitacional/colisional com o plasma.

**Forma atual.** Usar a forma volumétrica `f = J × B` para o fluido
ionizado, decompor em pressão + tensão e tratar `ρ_ion` como o
denominador apropriado em `a_B ≈ −∇u_B / ρ_ion`. Ver §2 da teoria.

---

### 3.4 Empurrar gás HI denso (`ρ_HI ~ 10⁻²² kg/m³`) → acoplar com WIM (`ρ_WIM ~ 10⁻²⁷ kg/m³`)

**Ideia inicial.** Calcular `a_B` usando a densidade do gás neutro HI
(que o radiotelescópio mede via linha de 21 cm).

**Por que foi substituída.** O gás HI é **neutro**: a Força de Lorentz
não atua sobre ele diretamente. Além disso, sua densidade é alta demais
para que o `B` em microGauss gere a aceleração de ordem `10⁻¹¹ m/s²`
exigida pelo gap.

**Forma atual.** Acoplar Lorentz ao **WIM** (Warm Ionised Medium), 5
ordens de magnitude menos denso. A densidade efetiva no denominador é
`ρ_ion = f_ion · ρ_gas`, com `f_ion` varrendo `{0.001, 0.01, 0.1}` no
grid. Ver §3.2 da teoria e a discussão do WIM em
`docs/educational/Simplified_Theory_PT.md` §2.

---

### 3.5 Ionização fixa em 1% → grid search agnóstico sobre `scenario_ionization_fraction`

**Ideia inicial.** Assumir `f_ion = 0.01` (1%) como valor "educado da
literatura" e calcular `a_B` apenas com esse valor.

**Por que foi substituída.** Risco de **viés de confirmação**: ajustar
um único parâmetro até que o resultado dê 100% de match seria justamente
o tipo de prática que invalidaria a publicação.

**Forma atual.** Grid search sobre nove cenários (3 escalas × 3 ionizações)
com cada cenário reportado, não apenas o melhor. O heatmap de sucesso
mostra qual **região** do espaço de parâmetros é consistente, não um
ponto único. Ver `docs/results_interpretation/Scientific_Evidence_Narrative_PT.md` §3.

---

### 3.6 Expectativa de "match 100%" como critério de sucesso → leitura agnóstica via `theoretical_gap_coverage_ratio`

**Ideia inicial.** Considerar o modelo confirmado se e somente se
`a_B ≈ Δa` exato (cobertura 100%).

**Por que foi substituída.** "100% exato" é numericamente implausível
em dados astronômicos ruidosos. Além disso, esperar uma única resposta
binária (substitui-ou-não a matéria escura) é cientificamente ingênuo.

**Forma atual.** A métrica é o logaritmo da razão de módulos
`theoretical_gap_coverage_ratio`. Valores em `[10^{-1}, 10^{+1}]` (uma
ordem de magnitude da unidade) já contam como "ordem de magnitude
correta" — uma evidência de contribuição **significativa**, não de
substituição completa. Confirmar substituição total exige propagação
Bayesiana de incertezas (item 2 do roadmap). Ver §5–6 do
`Scientific_Evidence_Narrative_PT.md`.

---

### 3.7 Aceitação acrítica do plasma como "neutro líquido total" → quase-neutralidade com correntes residuais

**Ideia inicial.** Aceitar que o plasma interestelar é estritamente
neutro (cargas se cancelam ⇒ Lorentz = 0).

**Por que foi substituída.** "Quase-neutro" significa `n_e ≈ n_i` em
média, mas correntes elétricas `J` podem existir e existem (Faraday
rotation, emissão síncrotron, mapas de polarização). Confundir os dois
conceitos é o erro pelo qual a astrofísica convencional descartou
prematuramente o eletromagnetismo galáctico.

**Forma atual.** Quasi-neutralidade entra na derivação para anular o
termo elétrico `qE` em (2.1), preservando o termo magnético `qv×B`.
A versão volumétrica `J×B` decorre disso sem contradição. Ver §2.2 da
teoria.

---

### 3.8 Calibração com fallback fixo → calibração dinâmica auto-referencial + auditoria de proveniência

**Ideia inicial.** Usar uma constante hardcoded para `S_ref` se o
cálculo dinâmico falhasse.

**Por que foi substituída.** O *fallback* mascarava falhas reais na
leitura do FITS do Planck e produzia `B` infundado para várias galáxias
sem que o relatório alertasse.

**Forma atual.** `S_ref` é a mediana do mapa de céu no `MagneticFieldCalibrator`.
Quando a calibração dinâmica falha, o pipeline marca a galáxia com
`fail_calibration` em `audit_flags`, e a coluna `is_reliable` se torna
`False`. A proveniência (`OPTIMAL` vs `SUBOPTIMAL`) é registrada no
`results/step_validation/validation_consistency_audit.md`. Ver §4 do `Scientific_Evidence_Narrative_PT.md`.

---

## 4. Renomeações terminológicas

| Antes | Agora | Por quê |
| ----- | ----- | ------- |
| "Estado A / Estado B" (Máquina de Estados) | "Regime A / Regime B" (decomposição contínua) | Linguagem mais próxima da MHD acadêmica. |
| "Universo Elétrico" | "Modelo E-G (Electrodynamic–Gravity)" / "MHD galáctica" | Distanciamento do EU pop-científico. |
| "Densidade de Corrente Crítica `J_critical`" | (removido; ver §3.1) | Sem suporte em equações de Maxwell. |
| "Direct_Evidence_Report.csv" (nome usado em apresentações iniciais) | `validation_scientific_evidence_report.csv` (nome actual sob `results/step_validation/`) | Alinhar nomes do código com a documentação. |

---

## 5. Itens a re-visitar em uma futura revisão

(Não foram refutados — apenas adiados.)

- Quantificação de **lentes gravitacionais por refração de plasma**:
  argumento bibliográfico esboçado no Knowledge_Checklist §C.1 mas não
  desenvolvido em nenhum documento canônico. Vale escrever um *appendix*
  do Draft Paper se a literatura `Plasma Lensing` (Pen & King 2012;
  Bannister et al. 2016) suportar.
- **Magnetogênese primordial** como contraponto à pergunta sobre CMB:
  mencionado no Knowledge_Checklist §C.2 mas sem citação. Subbarao et al.
  (1996) e Durrer & Neronov (2013) são os pontos de entrada na literatura.
- **Dínamo galáctico auto-excitado** como mecanismo de manutenção da
  estrutura magnética: mencionado em uma frase no raw (linha 579) e
  perdido na refatoração. Adicionar uma subseção em
  `EG_MODEL_THEORY_PT.md` §4 se houver demanda.

---

## 6. Como manter este arquivo atualizado

- Sempre que uma ideia for **explicitamente refutada/substituída** em
  conversa ou em código, adicionar uma entrada em §3 com data implícita
  via *commit history* (este arquivo é gitignored, então o histórico
  fica apenas no editor local).
- Sempre que um conteúdo for **removido** de um documento público, mover
  um resumo para §2 com a justificativa.
- Sempre que um conteúdo for **movido** entre arquivos públicos,
  registrar em §1.
- Em caso de divergência entre este arquivo e qualquer documento de
  `docs/theory/`, `docs/educational/`, `docs/academic/` ou
  `docs/results_interpretation/`, o documento público vence — este aqui
  é apenas memória operacional pessoal.
