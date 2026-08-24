---
id: 50_chugh-2017-actin-cortex-architecture-cell-surface-tension
paper_n: 50
title: "Actin cortex architecture regulates cell surface tension"
authors: "Chugh P., Clark A.G., Smith M.B., Cassani D.A.D., Dierkes K., Ragab A., Roux P.P., Charras G., Salbreux G., Paluch E.K."
year: "2017"
venue: "Nature Cell Biology 19(6):689-697"
doi: "10.1038/ncb3525"
paper_type: experimental + simulation
ffn_relevance: High
ffn_themes: [cortex, cortical-tension, parameter-source, validation-oracle, motor/myosin]
entities: [actin, myosin-ii, alpha-actinin, fascin, arp2-3, formin, cortical-actin, MYH9, MYH10, TPM4, hela]
methods: [atomic-force-microscopy, parallel-plate-confinement, dual-color-confocal-linescan, 3d-actomyosin-network-simulation, sirna-screen, pharmacology]
measurables: [cortical-tension, cortex-thickness, filament-length, network-connectivity, surface-tension]
keywords: [cortical-tension, cortex-architecture, filament-length-optimum, stall-force, connectivity, T0-230pNum, magnitude-oracle, KU-3.5]
tags: ["#cortical-tension", "#cortex-structure", "#parameter-source", "#validation-oracle", "#magnitude-oracle", "#KU-3.5", "#myosin-ii"]
has_transferable_params: true
membrane_vs_cortical: CORTICAL (actomyosin) — this IS the cortical-tension quantity KU-3.5 anchors
---

# [50] Actin Cortex Architecture Regulates Cell Surface Tension

**Tags:** #cortical-tension #cortex-structure #parameter-source #validation-oracle #magnitude-oracle #KU-3.5 #myosin-ii

> ⭐ **The KU-3.5 cortical-tension magnitude oracle.** This measures ACTOMYOSIN CORTICAL tension —
> the quantity the KU-3.5 band [0.35, 0.65] mN/m is about. Distinct from the membrane-tension
> papers (De Belly, Keren, Lüchtefeld) in this batch, which measure a different physical quantity.

| Field | Value |
|---|---|
| Authors | Chugh, Clark, Smith, Cassani, Dierkes, Ragab, Roux, Charras, Salbreux, Paluch |
| Year / Venue | 2017 / Nature Cell Biology 19(6):689-697 |
| DOI / ID | 10.1038/ncb3525 (manuscript PDF `emss-72183.pdf`) |
| Type | experimental (AFM + confinement + screens) + 3D actomyosin network simulation |
| Pages | 26 (author-manuscript) |
| ffn_cellsim relevance | High — validates the KU-3.5 band, supplies T₀ normalization + filament-length architecture optimum, and corroborates "tension = stall-force × connectivity, not myosin count" |

## 1. Summary
Combines AFM/confinement tension measurement, dual-colour cortex-thickness linescans, an siRNA/drug
architecture screen, and a 3D mechanistic actin+myosin+crosslinker simulation to show that cell
**surface (cortical) tension is set by cortex ARCHITECTURE** — filament length and network connectivity
— not by cortex thickness or raw myosin amount. Tension is non-monotonic in filament length, peaking
at an intermediate length. The simulation, anchored to myosin stall force, reproduces the magnitude and
the two-condition law for net tension.

## 2. Key anchors (value + units + provenance)
| Quantity | Value | Method / provenance | Maps to |
|---|---|---|---|
| Cortical tension, interphase HeLa | ~few hundred pN/µm | AFM constant-height (Fischer-Friedrich Eq.1), Fig 1c | band sanity: 1 mN/m = 1000 pN/µm ⇒ KU [0.35,0.65] = **350–650 pN/µm** |
| Cortical tension, mitotic HeLa | ~1000–1500 pN/µm | AFM, Fig 1c/3d | ≈1.0–1.5 mN/m, ~2–4× over band (rounded mitotic) |
| **Model normalization T₀** | **230 pN/µm**; coherent peak T/T₀≈1.6 ⇒ **~0.37 mN/m** | 3D sim anchored to myosin stall, Fig 4 legend/4c-d | **the magnitude oracle** — a correctly-built mechanistic cortex peaks at the band FLOOR (0.35–0.37) |
| Filament-length optimum for max γ | **~400–500 nm** (non-monotonic; ≈0 below ~250 nm) | sim, Fig 4d | cortex segment-length / mesh architecture target |
| Cortex thickness | interphase ~300–450 nm, mitotic ~200 nm | dual-colour linescan, Fig 1c | thickness ↑ ⇏ tension ↑ — thickness is NOT the knob |
| Two-condition rule | net tension needs (1) dense+connected enough AND (2) tensile-stress asymmetry | sim + theory, p.6 | the architectural law |
| Myosin NOT an architecture regulator | MYH9/10, blebbistatin, TPM4 perturbations leave thickness unchanged; tension drops on length perturbation with myosin unchanged | siRNA/drug screen, Fig 2a/3e-f | tension scale = stall-force × connectivity, not myosin count alone |

## 3. Relevance to ffn_cellsim
- **Validates the KU-3.5 band.** A coherent intermediate-length mechanistic cortex peaks at ~0.37 mN/m =
  the band floor; T₀=230 pN/µm. The band is right; the active-γ deficit measured in STAGE-1/2 is real.
- **Architecture constant.** Filament-length optimum ~400–500 nm is a direct cortex segment-length /
  mesh target for H.1.
- **Mechanism corroboration.** "Tension set by stall force × connectivity, not myosin count" independently
  matches STAGE-2's "force generation / soft-coupling ceiling" conclusion (connectivity z≈3 already OK,
  coherence OK ⇒ remaining factor is stall-force EXPRESSION / force delivery). The head-tension diagnostic
  (44% stall, ~3% bound, s_grip≈0) is the concrete mechanism by which stall force fails to express.
- **NOT a runtime fix.** Raising expressed force is core contractile physics, PI-gated; a k-stiffening
  sweep is hard-rule-forbidden without a PI nod. Enters only as oracle + architecture anchor.

## 4. Classification
Acceptance oracle (magnitude) + literature constant (filament-length architecture) + mechanism corroboration.
**Never a runtime mechanism.** CORTICAL (actomyosin) tension — the genuine KU-3.5 anchor.
