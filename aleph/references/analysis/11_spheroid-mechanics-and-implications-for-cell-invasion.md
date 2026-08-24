---
id: 11_spheroid-mechanics-and-implications-for-cell-invasion
paper_n: 11
title: "Spheroid mechanics and implications for cell invasion"
authors: "Ruben C. Boot et al."
year: "2021"
venue: "Advances in Physics: X"
doi: "10.1080/23746149.2021.1978316"
paper_type: review
ffn_relevance: Medium
ffn_themes: [spheroid-scale-context, cortex, FA/clutch, ECM/collagen, nucleus, cytoplasm, validation-oracle, parameter-source]
entities: [actin, myosin-ii, cadherin, e-cadherin, n-cadherin, p-cadherin, actomyosin-cortex, vimentin, fibronectin, collagen-i, matrigel, integrin, nucleus, cytoplasm, mcf7, mcf-10a, mda-mb-231, mda-mb-436, spheroid]
methods: [atomic-force-microscopy, micropipette-aspiration, tissue-surface-tensiometry, cavitation-rheology, optical-tweezers, traction-force-microscopy, microtweezers, hydrogel-mechanosensors, spheroid-fusion]
measurables: [surface-tension, cortical-tension, elastic-modulus, shear-modulus, viscosity, traction-stress, invasion-speed, migration-persistence, oxygen-tension, cell-shape-index]
keywords: [spheroid mechanics, tissue surface tension, cell sorting, differential adhesion hypothesis, cortical tension, cell invasion, jamming, EMT, viscoelasticity, breast cancer, collagen density, mechanobiology]
tags: ["#spheroid-mechanics", "#validation-oracle", "#parameter-source", "#breast-cancer", "#cell-sorting", "#cortical-tension", "#surface-tension", "#ECM-collagen", "#review"]
has_transferable_params: true
---

# [11] Spheroid mechanics and implications for cell invasion

**Tags:** #spheroid-mechanics #validation-oracle #parameter-source #breast-cancer #cell-sorting #cortical-tension #surface-tension #ECM-collagen #review

| Field | Value |
|---|---|
| Authors | Ruben C. Boot, Gijsje H. Koenderink, Pouyan E. Boukany |
| Year / Venue | 2021 / Advances in Physics: X (6:1, 1978316) |
| DOI / ID | 10.1080/23746149.2021.1978316 |
| Type | review |
| Pages | 36 |
| ffn_cellsim relevance | Medium — supracellular/spheroid-scale context, terminology, and measurement-oracle map for the multicell cohesion-overlay target; no fine-grained mechanism |

## 1. Summary
This is a soft-matter/mechanobiology review that surveys (i) the experimental toolbox for measuring the mechanics of 3D multicellular spheroids "from without" (tissue surface tensiometry, micropipette aspiration, AFM, microtweezers, spheroid fusion) and "from within" (cavitation rheology, hydrogel mechanosensors, optical-tweezer microrheology), and (ii) how the resulting mechanical parameters — effective tissue surface tension, viscoelastic moduli, and cortical tension — govern cell sorting and cell invasion into the ECM. The conceptual spine is the progression of cell-sorting models: the Differential Adhesion Hypothesis (DAH, adhesion-only), the Differential Interfacial Tension Hypothesis (DITH, adhesion balanced against actomyosin cortical tension), and the High heterotypic Interfacial Tension (HIT, ephrin-Eph repulsion) model, plus the jamming/unjamming (cell-shape-index) framework. The invasion half links these to breast-cancer spheroid assays: sorting architecture, EMT, collagen density, interstitial flow, and gap-junction-driven peripheral cell swelling/softening all modulate whether cells stay jammed (non-invasive) or unjam (invasive).

## 2. Problem & motivation
Tissue-scale mechanical response emerges from single-cell cytoskeletal/membrane properties plus cell-cell adhesion across length scales, but no single technique probes all scales. Spheroids are the dominant in-vitro 3D model that reproduces physiological structure (necrotic core / quiescent / proliferative shells with O2, pH, nutrient gradients). The review's stated gap: a concise survey of the experimental tools that quantify spheroid mechanics and how those mechanics drive cell detachment and invasion (cancer metastasis) was missing. Understanding the adhesion-vs-cortical-tension balance and its link to invasion can reveal biomarkers/therapeutic targets.

## 3. Methods / model
Pure literature review — no new data or simulation. It is organized as a methods/concept survey:
- **Probing from without:** tissue surface tensiometry (parallel-plate, Young-Laplace fit for apparent surface tension), micropipette aspiration (viscoelastic-drop tongue-flow), AFM (cantilever indentation -> apparent/pseudo-elastic modulus; dynamic sinusoidal mode for viscoelasticity at 0.5–4 Hz), microtweezers (dual force-sensing cantilevers -> Young's modulus), spheroid fusion (coalescence time -> bulk fluidity + surface tension).
- **Probing from within:** cavitation rheology (needle-induced bubble pressure-growth -> elastic modulus + cortical tension via binding-energy comparison), hydrogel mechanosensors (PAA microbeads, PNiPAAM thermoresponsive beads -> local pressure/elasticity via volume strain, Mooney-Rivlin beyond linear regime), optical tweezers (drag endocytosed bead -> cytoplasmic stiffness; insensitive to actin cortex).
- **Theory class for sorting:** liquid-drop analogy + free-energy minimization (DAH/DITH), interfacial-tension balance (DITH ratio of adhesion tension to cortical tension), repulsion (HIT), and vertex-/SPV-type jamming with a dimensionless cell shape index.
- **Length/time scales:** subcellular (optical tweezers, nm/fN) through cellular (cavitation, beads) to whole-spheroid (TST, MPA). Time scales: short-time elastic / long-time viscous; fusion over days.

## 4. Key results (quantitative)
This is a review, so numbers are cited from the surveyed literature, not measured here:
- **Diffusion/necrosis limit:** inner cells beyond ~200 µm from the spheroid edge undergo apoptosis -> necrotic core; in-vivo cells lie ~100–200 µm from nearest capillary (Sec. 2.1, p.5).
- **AFM viscoelastic regime:** dynamic AFM operated at functionally relevant frequencies 0.5–4 Hz; AFM indentation limited to depths < 10 µm (surface-only) (Sec. 2.2.2, p.12).
- **Micropipette aspiration geometry:** pipette radius ~3–4× smaller than cell diameter for single-cell aspiration (p.13).
- **Surface tension is mechanosensitive (MPA):** aspirated-spheroid surface tension increases with applied force (Guevorkian et al.); TST instead found surface tension independent of applied force (p.13–14) — an explicit oracle discrepancy.
- **Surface tension ∝ cadherin:** spheroid surface tension is a direct, linear function of cadherin expression level (DAH regime) (Sec. 2.2.1, p.8).
- **Optical-tweezer resolution:** positioning ~1 nm, force sensitivity ~50 fN (Sec. 2.3.2, p.16).
- **Breast-cancer microrheology:** spheroids seeded in alginate/Matrigel hydrogel with shear modulus ~300 Pa mimic in-vivo breast carcinoma; peripheral/invasive cells are softer, larger, longer, more dynamic than core cells (Han et al., Sec. 3.2, p.25).
- **Collagen density gates invasion mode:** low collagen ~1 mg/ml -> single-cell, unjammed, fluid-like invasive periphery; high collagen ~4 mg/ml -> collective migration / non-invasive solid-like (Fig. 7c, p.24–25).
- **Cell lines used in surveyed sorting/invasion work:** MCF-10A (non-tumorigenic), MDA-MB-231 and MDA-MB-436 (metastatic), spanning E-/N-/P-cadherin EMT shift; heterogeneous 1:1 MDA-MB-231:MCF-10A spheroids in collagen, with MCF-10A enclosing the malignant core by ~4 days due to higher proliferation, suppressing invasion (Fig. 7a, p.23).

## 5. Parameters & constants of interest
| Quantity | Value + units | Source in paper |
|---|---|---|
| Oxygen/nutrient diffusion limit -> apoptosis | ~200 µm from spheroid edge (in vivo ~100–200 µm from capillary) | Sec. 2.1, p.5 |
| AFM dynamic-mode frequency (viscoelastic) | 0.5–4 Hz | Sec. 2.2.2, p.12 |
| AFM max indentation depth | < 10 µm (surface technique) | Sec. 2.2.2, p.12 |
| MPA pipette radius vs cell | ~3–4× smaller than cell diameter | Sec. 2.3, p.13 |
| Optical-tweezer positioning resolution | ~1 nm | Sec. 2.3.2, p.16 |
| Optical-tweezer force sensitivity | ~50 fN | Sec. 2.3.2, p.16 |
| In-vivo breast-carcinoma-mimic hydrogel shear modulus | ~300 Pa (alginate/Matrigel) | Sec. 3.2, p.25 |
| Collagen-I density (single-cell invasion) | ~1 mg/ml | Fig. 7c, p.24–25 |
| Collagen-I density (collective/non-invasive) | ~4 mg/ml | Fig. 7c, p.24–25 |
| Microtweezer force range | <100 nN to ~1 mN | Sec. 2.2.2, p.12 |
| Surface tension vs cadherin | linear (qualitative law) | Sec. 2.2.1, p.8 |

Note: surface tension and cortical tension are given as concepts with SI units (N/m) in the definitions box (p.22), but the review does not tabulate specific numeric tension/modulus values for named cell lines; those live in the cited primary papers ([20], [70], [124], [135], [155], [165] etc.), which are the actual quantitative oracles to chase.

## 6. Relevance to ffn_cellsim
Classification: primarily **(d) spheroid/multicellular-scale context** with secondary **(b) parameter source / (a) validation-oracle pointer** value. It is NOT a mechanism reference for a fine-grained particle simulator and contains no equations, force laws, or rate constants that map to a runtime bond/motor/clutch.

Where it is useful:
- **Multicell cohesion overlay target (PI MCF7-spheroid program):** this is the conceptual scaffold for the platform's multi-cell validation tier. It frames tissue surface tension and cortical tension as the emergent observables that a mechanistic cohesion model (catch-bond cadherins KU-4.2 + actomyosin cortex H.1) must reproduce, and it names the exact measurement modalities (TST, MPA, cavitation, optical-tweezer microrheology) whose outputs become acceptance oracles. The c-coefficient (multi-cell cohesion) term of the PI's A/A0 = a + b/R + c/R² fit lives in this regime.
- **DAH -> DITH -> HIT progression as the right-altitude oracle hierarchy:** ffn_cellsim's emergent cell-sorting/cohesion should ideally recover the DITH balance (surface tension = ratio of adhesion tension to cortical tension) — exactly the mechanism (cadherin adhesion + myosin-II cortical contraction) that H.1 cortex + cadherin units build from particles. This review is the citation that justifies why both adhesion AND cortical tension are first-class, not adhesion alone (DAH breakdown when surface cells elongate).
- **EXTEND-track motivation:** explicitly flags the nucleus as an intracellular mechano-gauge (nuclear-membrane tension triggers actomyosin contractility) and cytoplasmic stiffness (optical-tweezer microrheology) — direct context for H.8/9/10 membrane/nucleus/cytoplasm additions, and for why cytoplasm/nucleus stiffness modulate invasion.
- **Cell-line + ECM grounding:** MCF-10A / MDA-MB-231 / MDA-MB-436 and collagen-I density gating (1 vs 4 mg/ml) match the project's breast-cancer overlay vocabulary; the ~300 Pa hydrogel and ~200 µm diffusion limit are usable anchor numbers for any spheroid-scale sanity context.

Where it is off-topic: the entire review treats cells as effective viscoelastic/liquid-drop objects (continuum + vertex/SPV altitude). It offers no filament-, motor-, or clutch-level constitutive content, so it cannot seed any H-unit runtime mechanism. Its role is to define WHAT the mesoscale must emerge into and WHICH experiments measure it — a target-and-oracle map, not a mechanism source.

## 7. Limitations & caveats
- Review, not primary data: every number is second-hand; quantitative oracles must be pulled from the cited primaries, not this paper.
- Scale gap: liquid-drop / surface-tension / jamming framing is two-to-three coarse-graining levels above a fine-grained HOOMD particle/bond model; no bridge equation is provided to descend from N/m tissue tension to per-cadherin/per-motor forces.
- The DAH/DITH liquid analogy is repeatedly flagged as an oversimplification (cancer aggregates don't round up, rough surfaces, active motility, jamming, EMT all break the equilibrium-liquid assumption) — so even the "oracle" laws are regime-limited.
- Several techniques (microtweezers, hydrogel mechanosensors) are described as novel/under-validated; mechanosensors lack a time-dependent component and may perturb the tissue (foreign-body response, pH/temperature sensitivity).
- In-vitro vs in-vivo sorting can disagree (boundary polarization, Xenopus long-timescale cadherin reorganization), cautioning against over-fitting any single in-vitro oracle.

## 8. Key figures / tables
- **Table 1 (p.4):** technique compendium — analyzed parameter (elastic modulus / surface tension / viscosity / interfacial tension / cytoplasmic stiffness) per method, from-without vs from-within; the master oracle-modality map.
- **Figure 6 (p.19):** DAH vs DITH vs HIT tension schematics + N-cadherin-ratio sorting confocal + actin-depolymerization (latrunculin A / cytochalasin D) surface-tension drop — the adhesion-vs-cortical-tension coaction evidence.
- **Figure 7 (p.24):** spheroid invasion — heterogeneous MDA-MB-231/MCF-10A architecture gating invasion; collagen vs Matrigel EMT-marker switch; collagen-density (1 vs 4 mg/ml) single-vs-collective migration phase separation.
- **Definitions box (p.22):** SI-unit definitions of cortical tension, stress, strain, elasticity (Young/shear/bulk), viscosity, viscoelasticity, tissue fluidity, tissue surface tension — useful canonical terminology anchor.

## 9. Notable quotes / citable claims
- "the tissue surface tension actually depends on a balance of adhesion, cortical tension and cortical elasticity ... the overall surface tension of a multicellular aggregate is determined by the ratio of adhesion tension to cortical tension, indicating a crossover from adhesion-dominated to cortical tension-dominated behavior" (DITH, p.18).
- "Spheroids have a surface tension that is a direct and linear function of their cadherin expression level" (DAH regime, p.8).
- "Cortical tension: The apparent surface tension of a cell, presumed to be dominated by myosin motor-driven contraction of the actin cortex and the interaction of the actin cortex with the membrane. SI unit: N/m." (definitions box, p.22).
- "Low collagen densities result in a locally unjammed invasive periphery resembling a fluid-like phase, while high collagen densities ensure non-invasive solid-like behavior" (p.25).
