# O Modelo E-G — Exposição Técnica em Linguagem Acessível

> **Público-alvo.** Leitores com cultura científica geral, mas sem formação
> em física de plasmas ou dinâmica galáctica. O texto é técnico no vocabulário,
> mas linear na exposição: cada seção introduz apenas o que a próxima precisa.
> Metáforas são deliberadamente evitadas; quando há analogias, elas são físicas.
> **Versão em inglês:** [`Simplified_Theory_EN.md`](Simplified_Theory_EN.md).
> **Versão formal (com a matemática completa):** [`../theory/EG_MODEL_THEORY_PT.md`](../theory/EG_MODEL_THEORY_PT.md).

---

## 1. O problema observacional

Galáxias espirais e anãs giram de um modo que não corresponde à gravidade
de sua matéria visível. Se apenas as estrelas, a poeira e o gás detectáveis
por telescópios estivessem presentes, a gravidade Newtoniana preveria que
estrelas nas regiões externas de uma galáxia desacelerassem com a distância,
da mesma forma que os planetas externos do Sistema Solar orbitam mais lentamente
que os internos. O que se mede, na verdade, é o oposto: a velocidade de
rotação permanece aproximadamente **plana** até os maiores raios em que se
consegue observar de forma confiável.

A formulação numérica da discrepância usa duas acelerações no mesmo raio `r`:

- a aceleração centrípeta **observada**, calculada a partir da velocidade
  de rotação: `a_obs = V_obs² / r`;
- a aceleração **bariônica**, calculada apenas a partir do perfil de massa
  visível: `a_bar = V_bar² / r`, com `V_bar² = V_gas² + V_disk² + V_bulge²`.

No disco externo, `a_obs > a_bar`. A diferença

```text
gap(r) = a_obs(r) − a_bar(r)
```

é a **lacuna de aceleração** (acceleration gap). A interpretação padrão da
cosmologia atribui essa lacuna a um halo de *matéria escura* não-bariônica:
uma substância que possui gravidade mas não interage com a luz, nunca foi
detectada diretamente e é invocada especificamente para preencher esse tipo
de aceleração faltante.

## 2. Do que o espaço é realmente feito

A interpretação da lacuna depende das suposições feitas sobre o meio que
preenche uma galáxia. Dois fatos da astrofísica estabelecida são relevantes:

1. **A maior parte do Universo visível é plasma.** Plasma é um gás no qual
   uma fração substancial dos átomos perdeu um ou mais elétrons e, portanto,
   carrega carga elétrica. Os meios interestelar e intergaláctico são
   predominantemente hidrogênio ionizado com traços de elementos mais
   pesados; conduzem eletricidade e respondem a campos magnéticos.
2. **Galáxias possuem campos magnéticos mensuráveis.** Vários canais
   observacionais independentes — rotação de Faraday de fontes polarizadas
   de fundo, emissão síncrotron difusa e os mapas de polarização linear do
   satélite Planck — convergem para a mesma conclusão: galáxias espirais e
   anãs são permeadas por campos magnéticos coerentes com intensidades da
   ordem de um a dez microgauss.

O componente de plasma mais relevante para esta discussão é o **Meio Ionizado
Morno (WIM, *Warm Ionised Medium*)**: hidrogênio ionizado pela radiação
ultravioleta de estrelas massivas, em temperatura da ordem de dez mil Kelvin
e densidade muito baixa (da ordem de `10⁻²⁷ kg/m³` nos halos galácticos).
Ele coexiste com o gás neutro (`HI`, da ordem de `10⁻²² kg/m³`) que
radiotelescópios detectam por meio da linha de 21 cm, mas é o WIM — não o
gás neutro — que responde diretamente a forças magnéticas.

## 3. Energia magnética como uma pressão

Um campo magnético de intensidade `B` armazena energia no espaço. A
quantidade de energia por unidade de volume é

```text
u_B(x) = B(x)² / (2 μ₀)
```

com `μ₀` sendo a permeabilidade do vácuo. Essa grandeza se comporta, nas
equações da magnetohidrodinâmica (MHD), exatamente como uma **pressão** —
tem as mesmas unidades de uma pressão (joules por metro cúbico) e, quando
sua distribuição espacial é não-uniforme, empurra fluido magnetizado de
regiões de `u_B` mais alto para regiões de `u_B` mais baixo. A grandeza
diretora não é, portanto, `u_B` em si, mas seu **gradiente** `∇u_B`.

A aceleração que uma parcela de fluido ionizado sofre por causa desse
gradiente é

```text
a_B(x) ≈ − ∇u_B(x) / ρ_ion(x)
```

onde `ρ_ion` é a densidade de massa do componente responsivo (ionizado).
O sinal negativo significa que a aceleração aponta de `u_B` alto para `u_B`
baixo. Isto **não** é a força de Lorentz microscópica dos livros-texto sobre
uma única partícula carregada; é a densidade de força macroscópica obtida
quando a mesma física subjacente é mediada sobre uma parcela de fluido. A
derivação completa está em `../theory/EG_MODEL_THEORY_PT.md`, §2.

## 4. Por que a escala importa

A gravidade de uma massa concentrada centralmente decai como `1/r²`. A
pressão magnética não decai da mesma forma: em uma estrutura organizada em
filamentos ou em geometria de disco, o gradiente de energia magnética
relevante decai muito mais lentamente que `1/r²`. O decaimento exato depende
da geometria — um filamento reto idealizado conduzindo corrente tem
`B ∝ 1/r`, de forma que seu gradiente de densidade de energia decai como
`1/r³` ao longo do eixo do filamento, mas muito mais lentamente
transversalmente ao lençol filamentar.

A consequência prática é que, em grandes raios galactocêntricos, onde a
gravidade da matéria visível ficou fraca, um gradiente de pressão magnética
que começou menor no disco interno pode acabar comparável em tamanho ou
maior. O raio de transição depende da geometria real do campo e da densidade
real do componente ionizado.

## 5. A hipótese E-G

O modelo E-G propõe que a dinâmica galáctica é descrita por dois regimes
co-existentes, não por duas teorias de gravidade concorrentes:

- **Regime A** domina onde quer que o potencial gravitacional bariônico
  seja forte: interiores estelares, sistemas planetários, os poucos kpc
  internos de um disco galáctico. A gravidade Newtoniana é exata dentro do
  erro observacional e a lacuna é pequena ou nula.
- **Regime B** é relevante onde a lacuna já se abriu: as partes externas
  dos discos galácticos e os halos difusos. Ali, o gradiente de pressão
  magnética definido acima age sobre o componente de gás ionizado (o WIM)
  e pode, em princípio, contribuir com uma aceleração de magnitude adequada
  para reduzir — ou substituir totalmente — a contribuição atribuída à
  matéria escura.

A hipótese é **falsificável**. Se o pipeline que calcula `gap` e `a_B` a
partir de dados reais encontrar `|a_B|` muitas ordens de magnitude menor
que `|gap|`, o modelo falha. Se `|a_B|` ficar na mesma faixa de `|gap|` e
os dois perfis compartilharem a mesma forma radial nos raios confiáveis,
o modelo passou em um teste observacional independente — embora insuficiente,
por si só, para deslocar a interpretação de matéria escura.

## 6. O que o pipeline efetivamente faz

O repositório implementa esse teste sobre o catálogo público SPARC (175
galáxias espirais e anãs com curvas de rotação bem medidas e decomposição
de massa) e sobre os mapas de polarização do Planck (os parâmetros de
Stokes de polarização linear `Q` e `U` da emissão de poeira galáctica).

Os passos, em ordem, são:

1. **Calcular a lacuna.** Para cada galáxia do catálogo, avaliar
   `a_obs(r)`, `a_bar(r)` e `gap(r)` a partir da tabela de rotação SPARC.
2. **Estimar a densidade responsiva.** Inverter a curva de velocidade só-de-gás
   em um perfil de densidade volumétrica de gás `ρ_gas(r)`. Multiplicar por
   uma fração de ionização livre `f_ion ∈ (0, 1]` para obter
   `ρ_ion(r) = f_ion · ρ_gas(r)`.
3. **Calibrar o campo magnético.** Converter a amplitude de Stokes do Planck
   `S(n̂) = √(Q² + U²)` na posição celeste da galáxia em uma intensidade
   de campo `B(n̂)` usando uma âncora dimensional (a densidade de energia
   da CMB, normalizada pela mediana do mapa; ver
   `../theory/EG_MODEL_THEORY_PT.md`, §3.3).
4. **Estimar o gradiente.** Comparar `u_B` no pixel central que aponta para
   a galáxia com a média de `u_B` dos pixels vizinhos HEALPix; dividir pela
   separação física para obter `∇u_B`.
5. **Calcular `a_B`.** Combinar o gradiente com `ρ_ion` para obter `a_B(r)`
   por casca radial.
6. **Comparar.** Reportar a razão de módulos
   `|⟨a_B⟩_edge| / |⟨gap⟩_edge|` e a correlação de Pearson entre os perfis
   radiais de `a_B` e `gap` sobre os raios que passam pelos *gates* de
   confiabilidade.

O pipeline faz isso para uma grade de nove cenários — três escalas de
gradiente (50 pc, 100 pc, 500 pc) cruzadas com três frações de ionização
(0,1 %, 1 %, 10 %) — e reporta o resultado de cada cenário, não apenas o
mais favorável.

## 7. Lendo o resultado honestamente

Uma razão de cobertura próxima de um significa que o *proxy* de pressão
magnética está na ordem de magnitude correta para fechar o déficit cinemático
nos raios onde a matéria escura é convencionalmente invocada. Uma razão
ordens de magnitude menor ou maior falsifica o *ansatz* no cenário
configurado.

Duas precauções estão embutidas no relatório:

- **Audit flags.** Cada galáxia é anotada com uma lista de condições sob as
  quais o resultado não deveria ser confiado: baixa latitude galáctica
  (contaminação por *foreground*), poucas cascas radiais confiáveis, piso
  de densidade dominando a inversão, gravidade estelar já explicando a
  dinâmica, inclinação baixa demais, calibração caindo em um valor de
  emergência.
- **Controles negativos sintéticos.** O pipeline é exercitado contra curvas
  de rotação sintéticas geradas a partir de um halo de matéria escura
  Navarro–Frenk–White (uma curva que, por construção, não tem gradiente
  magnético real). O pipeline **não** deve encontrar correspondência para
  esses controles.

A lista completa de limitações conhecidas — discrepância plano-do-céu vs
plano-do-disco, inversão indireta de densidade de gás, calibração
fenomenológica, termo de tensão magnética descartado, incertezas de
distância e inclinação — vive em [`EG_MODEL_THEORY_PT.md`](../theory/EG_MODEL_THEORY_PT.md), §5, e nos
guias de [`docs/results_interpretation/`](../results_interpretation/) que acompanham colunas CSV e *flags* de auditoria.

## 8. O que o projeto não está afirmando

O modelo E-G não afirma:

- que a gravidade Newtoniana esteja errada;
- que o modelo cosmológico padrão (ΛCDM) esteja errado como um todo;
- que observações de lentes gravitacionais (como o Bullet Cluster) sejam
  explicadas por efeitos magnéticos;
- que as flutuações da Radiação Cósmica de Fundo sejam explicadas sem um
  componente de matéria escura fria.

O modelo é restrito, por *design*, ao **problema das curvas de rotação
galáctica** — a anomalia empírica específica resumida na Seção 1. Sua
falsificação ou confirmação nesse único domínio é a única pergunta que o
pipeline responde.
