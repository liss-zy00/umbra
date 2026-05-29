# Modelo E-G — Refinamentos V3 (UX terminal + 3D opcional) — **entregue**

> **Arquivo.** Este ficheiro em `docs/personal/refinements/done/` regista o que **já está** no código relativamente ao antigo roadmap V3: fases nomeadas, *rich* onde aplicável, `tqdm`, `--logs`, README gerado, **indicadores de carga** (objectivo 2) e **entrada 3D opcional** via `python -m src.visualization` (objectivo 3), conforme a proposta foi refinada na implementação.  
> **Implementação actual (Maio 2026):** glifos `acceleration_gap`, anéis WIM (`effective_ionized_mass_density`), *streamlines* ilustrativas, `--galaxies auto`, manifest JSON, e `--viz-galaxies` no `main.py` — ver `src/visualization/render.py`, `scene.py`, `batch.py` e `VIZ_DISPLAY` em `src/config/thresholds.py`.  
> **Pendências** (UX terminal + *viz* avançado): [`../Model_Refinements_V4.md`](../Model_Refinements_V4.md). Teoria: [`../../theory/EG_MODEL_THEORY_PT.md`](../../theory/EG_MODEL_THEORY_PT.md).

---

## Estado entregue (vs `src/`)

| Entrega | Estado | Onde / notas |
|---------|--------|--------------|
| Fases nomeadas (Discovery → Validação → Gráficos → *Audit*) | **Feito** | `console_ui.py` + `runner.py` (`PIPELINE_PHASE_*`, `pipeline_phase_2_validation_title`) |
| Banners estilo regra + texto *rich* | **Parcial** | `console_ui.py` (`Rich` *Status* / *Console* onde aplicável) |
| `tqdm` por galáxia / lote | **Feito** | `runner.py` (`iter_progress`, desligado com `--logs`) |
| `--logs` (verbose, sem barras) | **Feito** | `cli.py` + `runner.py`; documentado no README **gerado** (`readme_builder.py`) |

---

## Objectivos 2 e 3 (este arquivo)

2. **Indicadores de carga** (I/O Planck, etc.) sem esconder erros — ver §3.  
3. **PyVista opcional** — entrada **`python -m src.visualization`** (`src/visualization/`), *bundle* curado; ver §4.

---

## Barras de progresso (§2.2)

- **Por galáxia** dentro de cada fase (ex.: loop SPARC com `tqdm`, opcionalmente
  silenciável com `--quiet` para CI).
- **Por cenário** na Fase 2 — Validação (9 células), se útil para *debug* de uma galáxia
  isolada (`--galaxies UGC…`).

## Modos de saída (§2.3)

- **`--pretty` ou `RICH_DISABLE=0` por defeito local** — definir política clara:
  em CI (GitHub Actions) manter saída minimalista; em terminal interactivo,
  activar *rich* por defeito ou via *flag*.
- **Níveis de *log***: `INFO` resumido para painel; detalhe técnico completo
  continua em ficheiro (`*.log` opcional) ou `DEBUG` com `--verbose`.

---

## 3. Carga estética e “lúdico” informativo

### 3.1 Indicadores de carregamento

- **I/O pesado:** leitura Planck HEALPix, parse SPARC, escrita em massa de CSV —
  mostrar *spinner* ou *pulse* com **etiqueta textual** (“Planck: Q/U ring …”),
  não apenas “Loading…”.
- **Transições de fase:** *banner* curto (1–3 linhas) com ícone/estilo consistente
  e **duração** opcional da fase anterior.

### 3.2 Movimento sem perda de rigor

- Animações **subtis** (actualização de contadores, barra suave); evitar *flash*
  contínuo que canse ou esconda *tracebacks*.
- Em caso de **excepção**, congelar o painel *pretty* e imprimir *stack* completo
  em texto normal (nunca esconder o erro atrás de UI).

### 3.3 Resumo final

Ao concluir: tabela *rich* com as mesmas métricas que hoje aparecem no bloco
`=== E-G pipeline complete ===`, mais **hiperligações relativas** aos CSV
gerados (opcional, se o terminal suportar *OSC 8*).

---

## 4. Visualização 3D (PyVista)

### 4.1 Âmbito

- **Opcional (empacotamento nomeado):** grupo `viz` em `pyproject.toml` / `pip install .[viz]` — **ainda em aberto**; hoje `pyvista` está listado em `requirements.txt`. Ver [`Model_Refinements_V4.md`](../Model_Refinements_V4.md).
- **Primeira entrega (feita):** módulo `src/visualization/` + `python -m src.visualization` + *hooks* `--viz-galaxies` no `main.py`; lê `results/galaxy_analysis/.../*_analysis.csv`, `*_summary.json`, e colunas de cenário vencedor via resumo JSON.

### 4.2 Cena 3D — camadas prioritárias (ordenadas)

1. **Disco / volume de referência** — cilindro ou disco achatado (raio máximo
   `r_max_kpc`, escala de espessura *toy*) para contexto.
2. **Plasma WIM (volume)** — *isosurface* ou volume semi-transparente (spec original); **implementação actual:** anéis pseudo-3D com mediana de `effective_ionized_mass_density` por bin radial (`--wim-rings`). *Isosurface* VTK completa → V4.
3. **Linhas de campo / *streamlines*** — seed regular no plano do disco;
   direcção do campo a partir de \(\mathbf{B}\) se disponível por casca, ou
   campo sintético alinhado ao gradiente de \(u_B\) com **etiqueta explícita**
   “ilustrativo” se o vector 3D completo não existir nos dados.
4. **Vectores de força** — setas radiais amostradas: \(|\mathbf{g}_{\mathrm{bar}}|\),
   \(|\mathbf{g}_{\mathrm{obs}}|\), magnitude do *gap*, e
   \(|a_B|\) (pressão ou Lorentz total conforme colunas); comparação visual
   “gravidade vs MHD” com escala comum e legenda SI.

### 4.3 Selecção de galáxias e modo comparação

- *Default:* intersecção das galáxias com PNG de evidência já gerados, ou
  primeiras *N* de `validation_reliable.csv`.
- *Override:* `--viz-galaxies DDO154,NGC3741,…` (lista curada; galáxias de controlo em `docs/theory/` e testes E2E).
- **`--viz-comparison` (opcional):** quando activo, o *render* exporta **dois
  *frames* lado a lado** (mesma câmara, mesma escala de cores / mesma barra de
  magnitudes) para cada galáxia:
  - **Painel A — “só bariões” (*toy* cinemático):** campo de aceleração radial
    derivado apenas de \(V_{\mathrm{bar}}\) (ex.: \(g_{\mathrm{bar}} = V_{\mathrm{bar}}^2/r\)
    a partir das colunas do `*_analysis.csv`). Visualmente: o “suporte” que a
    massa **visível** prediz, sem canal MHD.
  - **Painel B — “bariões + canal MHD do pipeline”:** o mesmo \(g_{\mathrm{bar}}\)
    **mais** a magnitude (e, se implementado, a componente radial assinada) de
    \(a_B\) / Lorentz conforme o melhor cenário — ou, de forma equivalente para
    narrativa, \(g_{\mathrm{obs}}\) *vs.* \(g_{\mathrm{bar}}\) com sobreposição de
    \(|a_B|\), desde que a legenda diga explicitamente o que está somado. O
    objectivo é **contraste visual imediato** no disco: “onde o *gap* existe no
    painel A, o que o modelo coloca no painel B”.

**Nota científica (obrigatória no *caption*):** isto **não** é uma simulação
N-body nem uma solução auto-consistente de Poisson + MHD; é uma **composição
gráfica** a partir dos perfis tabulares já produzidos pelo motor. O “argumento de
*venda*” é **pedagógico** — ajuda audiências a ver o *gap* e o termo magnético na
mesma geometria — e deve conviver com o *disclaimer* da §4.5.

### 4.4 Saída

- PNG / MP4 offscreen (`pyvista.Plotter(off_screen=True)`); resolução configurável.
- Documentar no *caption* do *render*: cenário \((\Delta x, f_{\mathrm{ion}})\),
  âncora, e `audit_flags` relevantes.

### 4.5 Riscos a documentar no *viewer*

- Planck dá **escalar** Stokes-derived \(B\) e gradiente de \(u_B\); **não** há
  mapa vectorial 3D completo do disco — qualquer *streamline* “bonito” pode ser
  **hipótese geométrica**, não observação directa.
- O *viewer* não substitui o motor; é **comunicação e exploração**.
- **`--viz-comparison`:** evitar linguagem “gravidade pura *vs.* gravidade real”
  ambígua; preferir rótulos do tipo **“aceleração prevista só com \(V_{\mathrm{bar}}\)”**
  *vs.* **“perfil observado / perfil com termo \(a_B\) do pipeline”**, para não
  confundir com teorias de modificação da gravidade.

### 4.6 Futuro — portfólio / apresentação pública (HTML + Three.js)

**Fora do escopo imediato** do MVP PyVista (§§4.1–4.5), mas a manter em mente: se o
projecto for usado como **portfólio** ou **demo pública**, considerar uma exportação
leve para **HTML interactivo** (por exemplo **Three.js** no browser). O *ecosystem*
PyVista/VTK inclui caminhos para **exportar cenas para visualização web** (p.ex.
`pyvista.plotting.itkwidgets`, *trame*, ou pipelines que geram `*.html` / glTF com
dependências adicionais — a escolha exacta fica para a fase de *spike*).

**Benefício:** um visitante abre um único ficheiro ou URL e roda as galáxias
**Supreme** (ou a sub-amostra de evidência) em 3D **sem instalar** Python, VTK nem
clone do repositório — útil para revisores leigos, *outreach* e *landing page* do
modelo. **Requisitos:** mesmos *disclaimers* pedagógicos do §4.5; não substituir o
motor nem implicar validação científica extra só por ter *spin* no browser.

---

## Apêndice — Plano PyVista secção 4 (ex-`to_do/Viz_PyVista_Section4_Plan_PT.md`)

Este bloco integra o plano de auditoria / fases que vivia em `docs/personal/refinements/to_do/Viz_PyVista_Section4_Plan_PT.md` (ficheiro removido após integração). **Regra de código** (do plano): nenhuma constante física nova no corpo das funções de *viz*; colunas do `*_analysis.csv` e `VIZ_DISPLAY` em `src/config/thresholds.py` para limites puramente de *display*.

### A.1 Auditoria inicial (histórico — pré-extensão Maio 2026)

| Item do spec (§4) | Estado na altura | Onde |
|---------------------|------------------|------|
| Módulo `src/visualization/` + `python -m src.visualization` | Feito | `__main__.py` → `render.run_cli` |
| Leitura `*_analysis.csv` + `*_summary.json` | Feito | `io.py` + `galaxy_from_disk.analysis_csv_path` |
| Disco + glifos bariões *vs.* observado | Feito | `scene.py` |
| Painel A/B + `--viz-comparison` + *captions* | Feito | `render.py` |
| PNG off-screen | Feito | `Plotter(off_screen=True).screenshot` |
| Teste de fumo | Feito | `tests/visualization/test_render_smoke.py` |

### A.2 Lacunas originais do plano — resultado após código

| § plano | Pedido original | Resultado (Maio 2026) |
|---------|-------------------|------------------------|
| 4.2.2 | WIM volume / *isosurface* \(\rho_{\mathrm{ion}}\) | **Parcial:** anéis 2.5D com mediana por bin (`wim_density_annulus_meshes`); não há `UniformGrid` + *marching cubes*. |
| 4.2.3 | *Streamlines* \(\mathbf{B}\) *toy* | **Feito** (opcional `--streamlines`); *caption* dedicada em `render.py`. |
| 4.2.4 | *gap* explícito nos vectores | **Feito** (opcional `--gap-glyphs` / `--viz-gap-glyphs`). |
| 4.3 | *Default* lista galáxias (CSV coorte) | **Feito:** `--galaxies auto`, `--tier`, `--max-n`. Não: intersecção automática com PNG de evidência. |
| 4.4 | MP4 | **Não** → [`Model_Refinements_V4.md`](../Model_Refinements_V4.md). |
| 4.1 | `dark-matter[viz]` opcional | **Não** → V4 / `requirements.txt` actual. |
| 4.6 | HTML / Three.js | **Não** → V4. |
| `main.py` | `--viz-galaxies` pós-pipeline | **Feito:** `cli.py`, `runner.py`, `batch.run_visualization_after_pipeline`. |

### A.3 Fases do plano (A–E) — estado

| Fase | Conteúdo | Estado |
|------|-----------|--------|
| **A** | `acceleration_gap` + anéis WIM + testes | **Feito** |
| **B** | Campo ilustrativo + *streamlines* + *caption* | **Feito** |
| **C** | `--galaxies auto`, manifest, README | **Feito** (manifest sem hash *git*) |
| **D** | *Packaging* + `main.py` | **Feito** (*packaging* `[viz]` ficou fora; ver V4) |
| **E** | MP4, resolução CLI, HTML | **Não** → V4 |

### A.4 Riscos e ordem (do plano)

- *Streamlines* exigem *disclaimer* → texto em `render.py` (`_CAPTION_STREAMLINES` + rodapé).  
- Volumes VTK pesados → estratégia adoptada: anéis antes de *isosurface* completa.  
- Ordem de trabalho seguida: A–D concluídos; E e *packaging* opcional ficam no V4.

### A.5 Constantes de *display*

`VIZ_DISPLAY` (`MAX_AUTO_GALAXIES`, `WIM_RING_BIN_COUNT`, `STREAMLINE_GRID_SIZE`) em `src/config/thresholds.py`.

