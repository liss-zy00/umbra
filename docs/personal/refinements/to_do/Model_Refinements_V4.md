# Modelo E-G — Refinamentos V4 (pendente: *viz* avançado)

> **Uso.** Notas pessoais em `docs/personal/refinements/`. Não cites em docs públicos.  
> **Âmbito.** Tudo o que **ainda não** está fechado depois da integração PyVista **A–D** (ver apêndice em [`done/Model_Refinements_V3.md`](done/Model_Refinements_V3.md)). Teoria: [`../theory/EG_MODEL_THEORY_PT.md`](../theory/EG_MODEL_THEORY_PT.md).

---

## 1. Visualização 3D — além do que o código faz hoje

Referência de implementação actual: `src/visualization/render.py`, `scene.py`, `batch.py`; *hooks* `src/pipeline/cli.py` + `runner.py`.

| Ideia (ex-plan §4 / plano PyVista) | Estado | Onde ir |
|------|--------|---------|
| WIM como *isosurface* / volume VTK completo (§4.2.2 “volume”) | **Não** | Hoje: anéis pseudo-3D (`wim_density_annulus_meshes`) — V4 se quiseres *marching cubes* |
| *Default* de lista = intersecção com PNG de evidência já gerados | **Não** | Hoje: `--galaxies auto` lê CSV *reliable* / *evidence* |
| MP4 / vídeo off-screen | **Não** | Fase E do plano integrado |
| Resolução de janela / export configurável por CLI | **Parcial** | `window_size` fixo em `render_galaxy_png` |
| `pip install .[viz]` / `[project.optional-dependencies]` | **Não** | `pyvista` em `requirements.txt` |
| Manifest JSON com checksum *git* opcional | **Não** | Hoje: só lista de paths + *timestamp* UTC |
| Export HTML / Three.js (§4.6) | **Não** | *Spike* futuro |

---

## 2. Onde ir a seguir

- Entregue (motor + V3 base + *viz* A–D): [`done/Model_Refinements_V1.md`](done/Model_Refinements_V1.md), [`done/Model_Refinements_V3.md`](done/Model_Refinements_V3.md).  
- Física / dados: [`to_do/Model_Refinements_V2.md`](to_do/Model_Refinements_V2.md).  
- README operador (gerado): `python src/main.py` após alterações em `readme_builder.py`.
