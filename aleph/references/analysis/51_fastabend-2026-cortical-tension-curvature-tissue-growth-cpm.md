---
id: 51_fastabend-2026-cortical-tension-curvature-tissue-growth-cpm
paper_n: 51
title: "Cortical tension links curvature to tissue growth in the cellular Potts model"
authors: "Fastabend et al."
year: "2026"
venue: "Physical Review E 113:024403"
doi: "10.1103/z152-x4l1"
paper_type: theory/simulation (cellular Potts model)
ffn_relevance: High (bridge)
ffn_themes: [cortical-tension, spheroid-scale-context, validation-oracle, numerics/methods]
entities: [cortical-tension, cellular-potts-model, perimeter-tension, curvature, tissue-growth, epithelium]
methods: [cellular-potts-model, young-laplace, interfacial-mechanics]
measurables: [perimeter-tension, surface-tension, curvature, growth-rate, radius]
keywords: [cortical-tension, CPM, perimeter-term, young-laplace, R-equals-lambda-over-sigma, growth-proportional-to-curvature, 1-over-R, A/A0-law, bridge-relation]
tags: ["#cortical-tension", "#bridge-relation", "#spheroid-mechanics", "#validation-oracle", "#young-laplace", "#curvature", "#cellular-potts"]
has_transferable_params: true
membrane_vs_cortical: CORTICAL — uses actomyosin cortical tension as the CPM perimeter term (not membrane tension)
---

# [51] Cortical Tension Links Curvature to Tissue Growth in the Cellular Potts Model

**Tags:** #cortical-tension #bridge-relation #spheroid-mechanics #validation-oracle #young-laplace #curvature #cellular-potts

> ⭐ **The published mechanism for the Layer-2 `A/A₀ = a + b/R + c/R²` law.** A paper-model ⇒ enters
> ffn_cellsim only as an acceptance oracle + bridge-relation, never as a runtime mechanism.

| Field | Value |
|---|---|
| Authors | Fastabend et al. |
| Year / Venue | 2026 / Physical Review E 113:024403 |
| DOI / ID | 10.1103/z152-x4l1 (`z152-x4l1.pdf`) |
| Type | theory / cellular Potts model |
| Pages | 13 |
| ffn_cellsim relevance | High (bridge) — supplies the published derivation linking single-cell cortical tension → curvature → growth, i.e. the source of the L2 spheroid `A/A₀(R)` law |

## 1. Summary
Embeds cortical tension as the CPM perimeter-energy term `λ_P(P−P₀)²` and shows that the resulting
interfacial mechanics make **tissue growth rate proportional to local interfacial curvature** (∝ 1/R),
with cortical-tension stretch ΔP serving as the proliferation signal and growth occurring only on
concave surfaces. A 2D Young-Laplace analog `R = λ/σ` sets the equilibrium radius from the line/surface
tension balance.

## 2. Key anchors (value + units + provenance)
| Relation | Form | Role |
|---|---|---|
| Cortical tension = CPM perimeter term | `λ_P (P − P₀)²` | how single-cell γ enters the tissue model |
| 2D Young-Laplace analog | **`R = λ/σ`** (radius set by line/surface-tension balance) | → 1/R dependence |
| Growth rate ∝ local interfacial curvature | growth ∝ 1/R; ΔP (cortical-tension stretch) is the proliferation signal; concave-only | **the published mechanism for `A/A₀ = a + b/R + c/R²`** |

(Quantitative parameters are CPM internal units — bridge-relation forms transfer, absolute a.u. do not.)

## 3. Relevance to ffn_cellsim
- **Bridge-relation oracle for L2 D2.** Our Layer-2 line already reproduces the endpoint
  (`A/A₀ = a + b/R + c/R²`, L2.5 G3 PASS, r²=0.98). Fastabend supplies the published derivation of *why*
  cortical tension produces that 1/R + 1/R² curvature dependence — closing the loop between the in-band
  `g_rigid` single-cell cortical tension (0.57 mN/m) and the emergent spheroid surface-tension law.
- **Connects to the in-band channel.** The bridge is NOT blocked by the active-γ floor: it runs on the
  already-in-band `g_rigid` structural cortical tension. This is the runnable-now Track-2.
- Pairs with Roffay (Young-Laplace `ΔP = σ(1/R+1/R')`) and Okuda (`Γ = cortical − adhesion`) as the
  L2-surface-tension oracle set.

## 4. Classification
Acceptance oracle + bridge-relation (single-cell cortical tension → curvature → growth → `A/A₀(R)`).
Paper-model ⇒ **never a runtime mechanism.** CORTICAL tension (actomyosin), not membrane tension.
