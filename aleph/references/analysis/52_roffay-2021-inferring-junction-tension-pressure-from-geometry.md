---
id: 52_roffay-2021-inferring-junction-tension-pressure-from-geometry
paper_n: 52
title: "Inferring cell junction tension and pressure from cell geometry"
authors: "Roffay C., Chan C.J., Guirao B., Hiiragi T., Graner F."
year: "2021"
venue: "Development 148:dev192773"
doi: "10.1242/dev.192773"
paper_type: theory/methods (geometric force inference) + experimental validation
ffn_relevance: High (bridge)
ffn_themes: [cortical-tension, spheroid-scale-context, validation-oracle, numerics/methods]
entities: [junction-tension, cortical-tension, hydrostatic-pressure, mouse-embryo, blastomere, curvature]
methods: [vertex-force-balance, young-laplace, micropipette-aspiration, geometric-inference]
measurables: [junction-tension, pressure, curvature, surface-tension-ratio]
keywords: [young-laplace, vertex-force-balance, tension-inference, outer-vs-inner-tension-ratio, dP-equals-t-K, 3D-curvature, bridge-relation, surface-tension-observable]
tags: ["#cortical-tension", "#bridge-relation", "#young-laplace", "#validation-oracle", "#tension-inference", "#spheroid-mechanics"]
has_transferable_params: true
membrane_vs_cortical: CORTICAL/junctional (actomyosin cortex sets junction tension) — not membrane tension
---

# [52] Inferring Cell Junction Tension and Pressure from Cell Geometry

**Tags:** #cortical-tension #bridge-relation #young-laplace #validation-oracle #tension-inference #spheroid-mechanics

> ⭐ **The Young-Laplace pressure↔tension↔curvature law for a real 3D cap + the surface-tension
> inference observable for L2.** Paper-model / measurement protocol ⇒ acceptance oracle + observable,
> never a runtime mechanism. ⚠️ Mouse-embryo system — the tension RATIO transfers, absolute a.u. do not.

| Field | Value |
|---|---|
| Authors | Roffay, Chan, Guirao, Hiiragi, Graner |
| Year / Venue | 2021 / Development 148:dev192773 |
| DOI / ID | 10.1242/dev.192773 (`DownloadCombinedArticleAndSupplmentPdf.pdf`, article+supplement combined) |
| Type | geometric force-inference theory/methods + micropipette validation |
| Pages | 19 (article + supplement) |
| ffn_cellsim relevance | High (bridge) — supplies the exact pressure↔tension↔curvature law for a 3D cap + a vertex-tension-inference observable to MEASURE the emergent spheroid surface tension as an L2 validation gate |

## 1. Summary
Develops and pipette-validates a method to infer cell-cell junction tension and hydrostatic pressure
from cell geometry alone, using vertex force balance and Young-Laplace. Applied to mouse early embryos,
it finds the outer (cell-medium) surface is a higher-tension boundary than interior cell-cell junctions —
the empirical signature that an aggregate's SURFACE tension is cortex-derived.

## 2. Key anchors (value + units + provenance)
| Relation / value | Form | Role |
|---|---|---|
| Vertex force balance | `Σ tᵢ = 0` | measurement protocol for tension inference |
| Young-Laplace (3D) | **`ΔP = t·K`, with `K = 1/R + 1/R'`** | exact pressure↔tension↔curvature law for a real 3D cap |
| Outer-vs-inner tension ratio | **outer (cell-medium) ≈ 1.6–2× interior cell-cell tension** (pipette-validated) | the empirical signature: aggregate surface = higher-tension boundary ⇒ surface tension is cortex-derived |

(⚠️ Mouse-embryo, arbitrary-unit tensions: the **~1.6× ratio transfers**, the absolute a.u. do NOT.)

## 3. Relevance to ffn_cellsim
- **Young-Laplace 3D oracle for L2 D2.** `ΔP = σ(1/R + 1/R')` is the pressure↔tension↔curvature law for a
  3D spheroid cap (supports the L2.6 "MCF7 = 3D cap not monolayer" geometry).
- **Surface-tension observable.** The vertex-tension-inference protocol lets us *measure* the emergent
  surface tension from simulated geometry — a candidate L2 surface-tension validation gate, in the same
  shape-observable family as the existing L2 hull/core-area metrics (read from GSD trajectories).
- **Bridge corroboration.** The outer>inner ~1.6× ratio is independent evidence that aggregate surface
  tension is cortex-derived — supporting the PI thesis that single-cell γ is the root of spheroid surface
  tension.
- Pairs with Fastabend (growth∝curvature, `R=λ/σ`) and Okuda (`Γ = cortical − adhesion`).

## 4. Classification
Acceptance oracle (Young-Laplace 3D) + measurement-protocol observable (vertex-tension inference) +
constant-anchor (outer/inner ratio, transfers as ratio only). Paper-model/protocol ⇒
**never a runtime mechanism.** CORTICAL/junctional tension, not membrane tension.
