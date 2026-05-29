# O Modelo E-G — Teoria Conceitual e Matemática

> **Status: documento foundational.** Este arquivo descreve a *proposição teórica*
> por trás do pipeline e a tradução matemática dos princípios físicos para o código.
> Ele é intencionalmente **atemporal**: não contém resultados de execução, valores
> por galáxia, nem caminhos de artefatos de runtime. Números concretos de qualquer
> execução específica vivem em `results/` e no [`README.md`](../../README.md)
> operacional. As limitações do ansatz estão resumidas na Secção 5 deste documento;
> para ler saídas do pipeline, use `docs/results_interpretation/`. A versão em inglês desta exposição está em
> [`EG_MODEL_THEORY_EN.md`](EG_MODEL_THEORY_EN.md).

---

## 1. Introdução — A Hipótese dos Estados Transicionais

As curvas de rotação galácticas são a anomalia empírica mais limpa da astrofísica
moderna. Em raios galactocêntricos onde a distribuição de massa visível (bariônica)
prevê um decaimento Kepleriano \(V \propto r^{-1/2}\), as medições espectroscópicas
mostram, de forma consistente, perfis de velocidade **planos**. A interpretação
padrão fecha essa lacuna com um **halo de matéria escura** esférico e
não-bariônico cuja gravidade fornece a aceleração centrípeta faltante

\[
\Delta a(r) \;=\; \underbrace{\frac{V_{\mathrm{obs}}^{2}(r)}{r}}_{a_{\mathrm{obs}}}
              \;-\; \underbrace{\frac{V_{\mathrm{bar}}^{2}(r)}{r}}_{a_{\mathrm{bar}}}
              \;>\; 0
\quad\text{(em }r\text{ grande)},
\qquad
V_{\mathrm{bar}} \;=\; \sqrt{V_{\mathrm{gas}}^{2}+V_{\mathrm{disk}}^{2}+V_{\mathrm{bulge}}^{2}} .
\tag{1.1}
\]

A **hipótese E-G** não contesta a gravidade Newtoniana em seu domínio de validade.
Ela adota uma postura complementar: a dinâmica galáctica pode ser
**dependente do estado em função da escala**, com a gravidade e o eletromagnetismo
atuando como forças fundamentais co-atuantes e mutuamente consistentes, cujo peso
relativo muda com a geometria do sistema.

* **Estado A — regime local / compacto.** Em corpos ligados gravitacionalmente e
  não magnetizados, o potencial Newtoniano domina. O campo de aceleração de uma
  massa pontual decai como \(1/r^{2}\); centros galácticos, sistemas estelares e
  órbitas planetárias são bem descritos neste limite.
* **Estado B — regime estendido / filamentar.** Em estruturas de plasma
  magnetizado (folhas de corrente, campos magnéticos de discos galácticos, a
  rede filamentar da cosmic web observada em surveys de estrutura em larga escala),
  a *estrutura espacial do campo magnético* permite tensões que decaem mais
  lentamente que \(1/r^{2}\). Para um filamento idealizado conduzindo corrente,
  a magnitude do campo escala como \(B \propto 1/r\) e o gradiente de pressão
  magnética correspondente não se anula nos raios onde os halos de matéria
  escura são convencionalmente invocados.

O teste E-G faz a pergunta **falsificável**: *pode um gradiente de pressão de
fluido magnetizado, ancorado em medições reais de polarização, fornecer uma
aceleração da ordem de magnitude correta e com a tendência radial correta nos
raios onde a lacuna (1.1) se abre?* Um "sim" não refuta a matéria escura; um
"não" — ou seja, um campo magnético inferido dimensionalmente implausível ou
cujo gradiente espacial é fraco demais por ordens de magnitude — falsificaria
diretamente este *ansatz* específico de pressão MHD. O pipeline é construído
de modo que ambos os desfechos sejam reportados com o mesmo rigor científico
(Seção 5).

---

## 2. A Ponte MHD — Da Força de Lorentz à Pressão Magnética

A ponte conceitual do modelo E-G é a passagem da força de Lorentz **microscópica**
sobre partículas carregadas para a tensão magnética **macroscópica** sobre um
elemento de fluido. Essa tradução é a física central que o pipeline implementa;
sem ela, a proposta seria metafísica.

### 2.1 Origem microscópica

Uma única partícula carregada de carga \(q\) e velocidade \(\mathbf{v}\) em um
campo eletromagnético \((\mathbf{E},\mathbf{B})\) sofre

\[
\mathbf{F}_{\mathrm{Lorentz}} \;=\; q\,\big(\mathbf{E} + \mathbf{v}\times\mathbf{B}\big).
\tag{2.1}
\]

Esta é a força de Lorentz dos livros-texto. Para um plasma interestelar, porém,
(2.1) não é o nível correto de descrição: as escalas de tempo relevantes são
longas comparadas ao período de gira das partículas, espécies de elétrons e
íons estão fortemente acopladas por colisões e feedback eletrostático, e a
grandeza dinâmica de interesse não é a trajetória de uma partícula, mas a
aceleração de massa de um elemento de fluido.

### 2.2 Forma de meio contínuo: a força de Lorentz volumétrica

Somando (2.1) sobre o conteúdo de carga de uma parcela de fluido e assumindo
quase-neutralidade (\(\rho_{e} \to 0\), de forma que o termo elétrico se anula
no referencial comóvel), a **força de Lorentz volumétrica por unidade de volume**
se reduz a

\[
\mathbf{f}_{\mathrm{Lorentz}} \;=\; \mathbf{J}\times\mathbf{B}
\;=\; \frac{1}{\mu_{0}}\big(\nabla\times\mathbf{B}\big)\times\mathbf{B},
\tag{2.2}
\]

onde a lei de Ampère (com corrente de deslocamento desprezada, o **limite MHD
não-relativístico**) foi usada para escrever a densidade de corrente em termos
do campo magnético. Uma identidade vetorial padrão decompõe (2.2) em duas
contribuições fisicamente distintas:

\[
\frac{1}{\mu_{0}}(\nabla\times\mathbf{B})\times\mathbf{B}
\;=\;
-\,\nabla\!\left(\frac{B^{2}}{2\mu_{0}}\right)
\;+\;
\frac{1}{\mu_{0}}(\mathbf{B}\cdot\nabla)\mathbf{B} .
\tag{2.3}
\]

O primeiro termo é uma **pressão** — age como o gradiente da densidade escalar
de **energia magnética**

\[
u_{B}(\mathbf{x}) \;=\; \frac{B(\mathbf{x})^{2}}{2\mu_{0}} \;=\; P_{B}(\mathbf{x}),
\tag{2.4}
\]

e aponta de regiões de maior energia magnética para regiões de menor energia
magnética. O segundo termo é a **tensão magnética**, anisotrópica na direção
de \(\mathbf{B}\); age como uma tensão restauradora ao longo de linhas de campo
curvas, fazendo o papel de uma componente "elástica".

### 2.3 O *ansatz* de gradiente de pressão em MHD ideal

O pipeline adota a **aproximação de gradiente de pressão em MHD ideal**: ao
longo da direção em que (2.3) é amostrada, o termo de tensão é descartado.
Dividindo a força de gradiente de pressão remanescente pela densidade de massa
inercial do fluido magnetizado, obtém-se a **aceleração de pressão magnética**

\[
\boxed{\;
\mathbf{a}_{B}(\mathbf{x})
\;\approx\;
-\,\frac{\nabla u_{B}(\mathbf{x})}{\rho_{\mathrm{ion}}(\mathbf{x})}
\;=\;
-\,\frac{\nabla\big(B^{2}/2\mu_{0}\big)}{\rho_{\mathrm{ion}}(\mathbf{x})} .
\;}
\tag{2.5}
\]

A equação (2.5) é a **equação operacional central** do teste E-G. O denominador
é a densidade de massa *responsiva* — a componente de gás ionizado que se
acopla eletromagneticamente ao campo. Sua modelagem é o assunto da Seção 3.

A convenção de sinal é física: uma *queda* de energia magnética no sentido
externo (\(\Delta u_{B} < 0\)) produz uma aceleração radial **positiva** para
fora, na leitura unidimensional simplificada de (2.5). Daí a analogia conceitual
com um "vento" de pressão magnética empurrando o plasma magnetizado para fora,
onde o campo está enfraquecendo.

A justificativa para descartar o termo de tensão é, ela mesma, uma afirmação
falsificável (ver Seção 5.4); é o *ansatz* MHD de primeira ordem mais simples
que permite a um pipeline puramente de dados produzir uma estimativa de ordem
de magnitude de \(\mathbf{a}_{B}\) sem exigir uma reconstrução tridimensional
global da topologia do campo magnético galáctico.

---

## 3. Implementação Matemática — Da Teoria ao Código

O pipeline é a tradução numérica literal de (1.1) e (2.5) em uma sequência de
estimadores, cada um com uma única responsabilidade física. Esta seção declara
o que cada estimador computa e *por quê*; o código-fonte em `src/physics/`
corresponde a esta seção linha por linha.

### 3.1 O déficit cinemático de aceleração

O primeiro estimador (`BaryonicGapCalculator`) lê a tabela de rotação SPARC e
produz acelerações em SI:

\[
a_{\mathrm{obs}}(r) = \frac{V_{\mathrm{obs}}^{2}(r)}{r},\qquad
a_{\mathrm{bar}}(r) = \frac{V_{\mathrm{bar}}^{2}(r)}{r},\qquad
\mathrm{gap}(r) = a_{\mathrm{obs}}(r) - a_{\mathrm{bar}}(r).
\tag{3.1}
\]

A síntese bariônica compõe explicitamente as componentes SPARC

\[
V_{\mathrm{bar}}(r) \;=\; \sqrt{V_{\mathrm{gas}}^{2}(r) + V_{\mathrm{disk}}^{2}(r) + V_{\mathrm{bulge}}^{2}(r)} ,
\tag{3.2}
\]

que é o análogo rotacional da velocidade bariônica usada na literatura
RAR/BTFR. A incerteza de inclinação se propaga linearmente em uma banda
simétrica de ±1σ sobre \(a_{\mathrm{obs}}\).

### 3.2 O estimador de densidade de gás

O denominador de (2.5) requer uma densidade de massa **ionizada** volumétrica.
O estimador (`GasDensityEstimator`) usa a inversão de cascas esféricas da
velocidade SPARC só-de-gás, avaliada por **dois caminhos equivalentes** que
devem concordar por construção:

\[
M_{\mathrm{gas}}^{(A)}(r)
\;=\;
\frac{V_{\mathrm{gas}}^{2}(r)\,r}{G},
\qquad
M_{\mathrm{gas}}^{(B)}(r)
\;=\;
\frac{\max\!\big(V_{\mathrm{bar}}^{2}-V_{\mathrm{disk}}^{2}-V_{\mathrm{bulge}}^{2},\,0\big)\,r}{G}.
\tag{3.3}
\]

Pela ortogonalidade SPARC (3.2), os dois caminhos são matematicamente
idênticos; o estimador avalia ambos para que qualquer vazamento entre
componentes do catálogo seja auditável como um resíduo por casca não-nulo.
O caminho B subtrai explicitamente o **potencial estelar** do potencial
bariônico total — uma casca onde estrelas já explicam a dinâmica resulta em
um resíduo nulo, e a escolha conservadora \(\max(A,B)\) impede que o estimador
infle o reservatório inercial.

A densidade volumétrica decorre de uma decomposição em cascas cilíndricas.
Com espessura de casca \(h(R)\) configurada a partir da extensão do catálogo
e volume de casca \(\Delta V_{i} = \pi (r_{i}^{2}-r_{i-1}^{2})\,h\), a massa
e a densidade da casca são

\[
\Delta M_{i} \;=\; M_{\mathrm{gas}}(r_{i}) - M_{\mathrm{gas}}(r_{i-1}),\qquad
\rho_{\mathrm{gas}}(r_{i}) \;=\; \frac{\Delta M_{i}}{\Delta V_{i}} .
\tag{3.4}
\]

A densidade ionizada efetiva usada em (2.5) é

\[
\rho_{\mathrm{ion}}(r) \;=\; f_{\mathrm{ion}} \cdot \rho_{\mathrm{gas}}(r) ,
\tag{3.5}
\]

onde \(f_{\mathrm{ion}} \in (0, 1]\) é um botão de cenário que explora a fração
desconhecida de ionização do meio ionizado morno que responde ao gradiente
magnético na escala de tempo relevante.

Um **piso de densidade** impede que cascas onde a inversão cinemática produz
densidades não-finitas ou não-positivas corrompam (2.5); um **gate de
dominância estelar** remove cascas onde a gravidade estelar fornece a maior
parte do orçamento de velocidade ao quadrado, já que a inversão só-de-gás
deixa de ser um traçador significativo de densidade ali.

### 3.3 O calibrador Planck — Amplitudes de Stokes para Tesla

O satélite Planck mede a polarização linear via os parâmetros de Stokes
\(Q(\hat n)\) e \(U(\hat n)\) no céu. Sua amplitude

\[
S(\hat n) \;=\; \sqrt{Q(\hat n)^{2} + U(\hat n)^{2}}
\tag{3.6}
\]

é dimensional em unidades de mapa, *não* em Tesla no meio interestelar. O
calibrador (`MagneticFieldCalibrator`) converte (3.6) em uma intensidade de
campo magnético via uma **ponte dimensional**:

\[
B(\hat n) \;=\; \sqrt{2\,\mu_{0}\,u_{\mathrm{ref}}}\;\cdot\;\frac{S(\hat n)}{S_{\mathrm{ref}}} ,
\tag{3.7}
\]

onde \(u_{\mathrm{ref}}\) é uma densidade de energia de referência fixa (em SI
J/m³) e \(S_{\mathrm{ref}}\) é uma escala de calibração derivada do próprio mapa.

Dois pontos interpretativos importantes:

1. \(u_{\mathrm{ref}}\) é a **densidade de energia da CMB**
   \(u_{\mathrm{CMB}} = a_{\mathrm{rad}} T_{\mathrm{CMB}}^{4}\) usada como
   *âncora dimensional*. Isto **não** é uma afirmação de que o campo magnético
   galáctico é igual ao campo da CMB. É uma escolha fenomenológica: a CMB
   fornece uma escala cosmologicamente uniforme de densidade de energia com
   as unidades corretas, a amplitude de polarização é normalizada em relação
   a si mesma e o \(B\) resultante carrega as dimensões apropriadas de Tesla
   para uso em (2.4) e (2.5).
2. \(S_{\mathrm{ref}}\) é a **mediana** de (3.6) sobre o mapa de céu carregado,
   o que torna a calibração **auto-referencial** ao conteúdo do mapa e
   independente da escala absoluta do mapa.

A sensibilidade à escolha de \(u_{\mathrm{ref}}\) é um eixo explícito de
refinamento (Seção 5.3).

### 3.4 O gradiente de pressão magnética

O gradiente \(\nabla u_{B}\) em (2.5) é estimado por uma **diferença finita
escalar** entre o pixel HEALPix central apontando para a galáxia e a média
de seus vizinhos angulares:

\[
\Delta u_{B} \;=\; u_{B}^{(\mathrm{neighbours})} - u_{B}^{(\mathrm{centre})},
\qquad
\frac{\nabla u_{B}\,\cdot\,\hat e}{1} \;\approx\; \frac{\Delta u_{B}}{\Delta x} ,
\tag{3.8}
\]

onde \(\Delta x\) é a separação física entre o centro e seus vizinhos, e
\(\hat e\) é o vetor unitário (desconhecido) que liga essas duas posições no
referencial do disco. O pipeline rastreia **duas** escolhas paralelas de
\(\Delta x\): um botão exploratório de cenário (50–500 pc, variado pela grid)
e um \(\Delta x_{\mathrm{phys}} = D_{\mathrm{galaxy}}\,\theta_{\mathrm{neighbour}}\)
**físico**, calculado a partir da distância catalogada e do espaçamento
angular dos vizinhos HEALPix. Os dois valores produzem duas colunas paralelas
de aceleração; a leitura honesta dos dados é a coluna física, e a comparação
entre elas quantifica a discrepância geométrica de ordem de magnitude.

### 3.5 A aceleração MHD integrada

Combinando (3.4)–(3.8), o pipeline avalia (2.5) casca a casca:

\[
a_{B}(r) \;=\; -\,\frac{\Delta u_{B}}{\Delta x \cdot \rho_{\mathrm{ion}}(r)} .
\tag{3.9}
\]

Diagnósticos de sanidade por casca derivados de (3.9) incluem:

* a razão **plasma-β** de pressão térmica e magnética,
  \(\beta = 2\mu_{0} n k_{B} T / B^{2}\), que marca o regime onde a hipótese
  de gradiente de pressão em MHD ideal é localmente consistente (\(\beta < 1\))
  versus o regime onde a pressão térmica suplanta a pressão magnética
  (\(\beta \gg 1\), onde o *ansatz* está *fora de regime*);
* a **velocidade de Alfvén** \(v_{A} = B / \sqrt{\mu_{0}\,\rho_{\mathrm{ion}}}\),
  velocidade na qual perturbações magnéticas transversais se propagam;
* uma **banda de plausibilidade** sobre o \(B\) inferido, de forma que valores
  muito fora de campos de disco consistentes com a literatura apareçam como
  avisos de auditoria.

A métrica de comparação do teste é a razão de módulos no disco externo

\[
\theta_{\mathrm{coverage}}
\;=\;
\frac{|\langle a_{B}\rangle_{\mathrm{edge}}|}{|\langle\mathrm{gap}\rangle_{\mathrm{edge}}|}
\tag{3.10}
\]

avaliada nos raios onde halos de matéria escura são convencionalmente
invocados. Valores próximos da unidade indicam que o *proxy* de pressão
magnética está na *ordem de magnitude correta* para fechar o déficit cinemático
naqueles raios; valores muito distantes da unidade falsificam diretamente o
*ansatz* no cenário configurado. Uma correlação de Pearson complementar, ao
longo de raios confiáveis, verifica independentemente se a **forma** dos dois
perfis concorda, e não apenas sua amplitude integrada.

---

## 4. Premissas e Evidências de Apoio

O *ansatz* E-G se apoia em cinco pilares empíricos da astrofísica estabelecida.
Nenhum deles é exclusivo deste trabalho; o pipeline apenas os combina em um
único teste falsificável.

### 4.1 Campos magnéticos galácticos existem

Emissão síncrotron difusa, rotação de Faraday de fontes polarizadas de fundo
e mapas diretos de polarização do Planck e de experimentos CMB anteriores
confirmam todos campos magnéticos coerentes permeando galáxias espirais e
anãs com intensidades típicas \(\mathcal O(1\text{–}10\,\mu\mathrm{G})\). A
*existência* de um \(B(\mathbf{x})\) galáctico não é, portanto, uma premissa,
mas um *input* observacional; o que (2.5) quantifica é sua *contribuição*
dinâmica.

### 4.2 Magnetohidrodinâmica ideal se aplica em escalas galácticas

Para condições interestelares típicas (alta condutividade, baixa resistividade,
grandes escalas de comprimento), a aproximação **MHD ideal**
\(\mathbf{E} + \mathbf{v}\times\mathbf{B}=0\) é válida com alta precisão. A
descoberta de Hannes Alfvén de que linhas de campo magnético são *congeladas*
em um fluido perfeitamente condutor é a licença para escrever (2.2) como uma
divergência de tensões magnéticas e reter apenas a componente de pressão nas
maiores escalas. Correções resistivas, tensores de pressão anisotrópicos e
pressão de raios cósmicos são reais, mas de segunda ordem; tratá-los é o
assunto do programa de refinamento da Seção 5.

### 4.3 O Universo é filamentar nas escalas relevantes

Surveys de estrutura em larga escala (SDSS, eBOSS, DESI) e simulações N-body
revelam consistentemente a **cosmic web**: uma rede de filamentos, lâminas e
voids na qual matéria escura e bariônica estão arranjadas. Em um filamento
reto idealizado que conduz corrente, a magnitude do campo magnético escala
como \(B \propto 1/r\), de modo que o gradiente de pressão magnética
correspondente é não-trivial em distâncias comparáveis à separação
inter-filamentar. A hipótese E-G propõe que as partes *externas* dos discos
galácticos vivem nessa geometria, em vez da geometria de halo esférico
assumida pelo ΛCDM — uma alternativa geométrica, não de partículas exóticas.

### 4.4 Os parâmetros de Stokes carregam informação de campo magnético

Os mapas de polarização \((Q,U)\) do satélite Planck medem a polarização
linear da emissão térmica de poeira. Para emissão opticamente fina com
polarização intrínseca aleatória alinhada por um campo ordenado, o ângulo de
polarização rastreia a orientação projetada do campo magnético e a *fração*
de polarização depende monotonicamente da intensidade do campo. A equação
(3.7) é a leitura **dimensional** mínima dessa informação: não assume modelo
específico de grãos de poeira, apenas que a amplitude de polarização é
monotônica em \(B\) e uniforme sobre o suporte de calibração.

### 4.5 Gravidade Newtoniana é o limite local

A lacuna cinemática (1.1) é construída assumindo gravidade Newtoniana a partir
da decomposição de massa SPARC. A hipótese E-G **não** modifica isso. No regime
local — interiores estelares, sistemas planetários, os poucos kpc internos
de uma galáxia onde o componente bariônico domina — Newton é exato dentro
da precisão observacional e a lacuna é pequena ou nula. O termo de pressão
MHD só se torna relevante nos raios onde o déficit já se abriu.

---

## 5. Parâmetros de Refinamento (Limitações Inerentes)

A lista abaixo cataloga as **fragilidades conhecidas** da formulação atual.
Cada item é um eixo ao longo do qual uma versão refinada do modelo poderia
ser falsificada com mais precisão; nenhum deles invalida a estrutura
conceitual subjacente das Seções 1–4. Corrigir ou mitigar cada limitação exige alargar o
pipeline (*providers*, diagnósticos e relatórios estatísticos) mantendo este núcleo teórico.

### 5.1 Discrepância geométrica entre plano do céu e plano do disco

O gradiente em (3.8) é amostrado ao longo do eixo **projetado no céu**
definido pelos vizinhos HEALPix, enquanto o gradiente *físico* que contribui
para a aceleração radial em (2.5) é a componente radial no **plano do disco**.
O escalar projetado \(\Delta u_{B}/\Delta x\) é, portanto, um limite superior
sobre a componente radial do verdadeiro gradiente; resolver a diferença
direcional exigiria reconstrução multi-vizinhos ou uma deprojeção do mapa de
polarização para o referencial do disco.

### 5.2 Densidade de gás indireta a partir da inversão de cascas esféricas

A equação (3.4) inverte uma velocidade de rotação em densidade volumétrica
sob a suposição implícita de simetria esférica da massa contida. Um disco
de HI viola explicitamente essa suposição: o gás está concentrado em um disco
fino, não em uma esfera. As verificações cruzadas dos caminhos A e B (3.3) e
o *gate* de dominância estelar impedem o estimador de inflar massa sob
contaminação estelar, mas não transformam a densidade de gás residual em uma
densidade real do ISM. Uma validação cruzada contra perfis de massa de HI
catalogados é o próximo passo natural.

### 5.3 Status fenomenológico da âncora de Stokes

A equação (3.7) fornece as **dimensões** corretas para \(B\); não fornece
uma magnitude de primeiros princípios. Uma escolha diferente de
\(u_{\mathrm{ref}}\) — digamos, equipartição com uma pressão turbulenta fixa
— reescalaria \(B\) e, portanto, \(a_{B}\) por um fator constante. Uma
varredura de sensibilidade sobre escolhas de âncora fisicamente motivadas é
necessária antes que a razão *absoluta* de cobertura (3.10) possa ser
interpretada com confiança.

### 5.4 *Ansatz* MHD ideal: tensão e pressão não-térmica ausentes

O pipeline retém apenas a componente de gradiente de pressão magnética de
(2.3), descartando o termo de tensão magnética \(({\bf B}\cdot\nabla){\bf B}/\mu_{0}\)
e tratando o plasma como ideal (sem resistividade, sem tensor de pressão
anisotrópico, sem contribuição de raios cósmicos, sem pressão turbulenta).
Onde a topologia do campo magnético é altamente ordenada, a tensão pode
competir com a pressão; raios cósmicos em anãs ativamente formadoras de
estrelas podem carregar uma pressão não-térmica comparável à térmica. O
diagnóstico plasma-β expõe o regime em que o *ansatz* está localmente fora
de regime, mas não o corrige.

### 5.5 Incertezas de inclinação e distância como bandas simétricas

A incerteza de inclinação é propagada em uma banda simétrica de ±1σ sobre
\(a_{\mathrm{obs}}\) — mas a razão de cobertura (3.10) em si é calculada a
partir do valor central, não das bordas da banda. Distâncias SPARC carregam
erros não-triviais que afetam linearmente o \(\Delta x_{\mathrm{phys}}\)
físico em (3.8). Ambos devem ser propagados através de (3.10) usando um
Monte-Carlo Bayesiano antes que qualquer valor de cobertura próximo da unidade
possa ser considerado *estatisticamente* distinguível do nulo de não-efeito.

### 5.6 Escopo de falsificabilidade dos controles sintéticos

O modelo é exercitado contra controles negativos NFW sintéticos — curvas de
rotação de matéria escura sem gradiente magnético real — para garantir que
o pipeline **não consiga** combiná-las por construção. Esse único teste limita
uma classe específica de falsos positivos. Ele não exclui mecanismos
alternativos não-MHD (inércia modificada, estrutura bariônica não-esférica,
artefatos de *beam-smearing* em \(V_{\mathrm{obs}}\) reais) que poderiam
reproduzir a forma da lacuna observada sem invocar pressão magnética.
Campanhas de falsificabilidade mais fortes — varreduras de sensibilidade de
âncora, *cross-checks* contra RAR/BTFR empíricos, correção de tensão
magnética — estão implementadas ou acompanhadas no código e nos artefactos de
validação em `results/` (ver README do repositório e `docs/results_interpretation/`).

---

## 6. Mapa de Leitura

* **Programador entrando no codebase.** A Seção 3 explica *o que
  `calculate_magnetic_pressure_acceleration` está computando* — ela avalia
  (3.9) usando `MagneticFieldCalibrator` (Seção 3.3),
  `GasDensityEstimator` (Seção 3.2), e `MagneticPressureGradientEstimator`
  (Seção 3.4), compostos pelo orquestrador `PhysicsEngine`. A lacuna
  cinemática está na Seção 3.1 (`BaryonicGapCalculator`).
* **Físico avaliando a proposição.** As Seções 1–2 estabelecem a hipótese, a
  ponte MHD de (2.1)–(2.2) até o *ansatz* de gradiente de pressão (2.5) e
  os pilares empíricos (Seção 4). A Seção 5 lista as fragilidades honestas;
  o *software* e os guias de resultados descrevem o que já está exercitado em código
  versus o que permanece em aberto.
* **Revisor perguntando o que o teste conclui para um dataset específico.**
  Essa pergunta é respondida pelos artefatos de runtime sob `results/`
  e pelo resumo operacional em [`README.md`](../../README.md). Esses arquivos
  *não* fazem parte deste documento por design — a teoria deve permanecer
  reproduzível entre execuções e desacoplada de qualquer execução particular.
