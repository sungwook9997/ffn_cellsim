---
id: 14_three-dimensional-in-vitro-tumor-spheroid-models-for-evaluat
paper_n: 14
title: "Three-Dimensional In Vitro Tumor Spheroid Models for Evaluation of Anticancer Therapy: Recent Updates"
authors: "Nayak, P. et al."
year: "2023"
venue: "Cancers (MDPI)"
doi: "10.3390/cancers15194846"
paper_type: review
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, tangential]
entities: [mcf7, mda-mb-231, e-cadherin, vimentin, integrin, fibroblast, ecm, collagen-i, spheroid, organoid, tumor-microenvironment]
methods: [hanging-drop, liquid-overlay, microfluidics, magnetic-levitation, spinner-flask, confocal-live-imaging, flow-cytometry]
measurables: [proliferation-rate, cell-viability, drug-ic50, spheroid-size, oxygen-tension]
keywords: [3d-spheroid, tumor-microenvironment, drug-screening, nuclear-medicine, nanocarrier, photodynamic-therapy, organoid, anticancer-therapy, necrotic-core, cell-cell-adhesion]
tags: ["#spheroid-scale-context", "#tumor-spheroid", "#drug-screening", "#review", "#breast-cancer", "#tangential"]
has_transferable_params: false
---

# [14] Three-Dimensional In Vitro Tumor Spheroid Models for Evaluation of Anticancer Therapy: Recent Updates

**Tags:** #spheroid-scale-context #tumor-spheroid #drug-screening #review #breast-cancer #tangential

| Field | Value |
|---|---|
| Authors | Nayak, P.; Bentivoglio, V.; Varani, M.; Signore, A. |
| Year / Venue | 2023 / Cancers (MDPI), 15, 4846 |
| DOI / ID | 10.3390/cancers15194846 |
| Type | review (biomedical/pharmacology narrative review) |
| Pages | 23 |
| ffn_cellsim relevance | Low — pure oncology/drug-screening review; no mechanics, no equations, no transferable parameters |

## 1. Summary
This is a narrative biomedical review of 3D in vitro tumor-spheroid models as a middle ground between 2D monolayer culture and animal models for anticancer-therapy evaluation. It surveys (i) why 3D spheroids better replicate solid-tumor features (cell-cell adhesion, necrotic/proliferating/quiescent zones, drug-resistance, gene-expression patterns); (ii) spheroid fabrication and characterization techniques (hanging-drop, liquid-overlay, spinner flask, magnetic levitation/bio-printing, microfluidics, with confocal/multi-photon/UV-Vis/flow-cytometry readouts); (iii) spheroid use across cancer types (prostate, liver, breast, pancreatic, thyroid, lung, ovarian); and (iv) theragnostic applications — nuclear-medicine therapy, stem-cell therapy, photodynamic therapy, immunotherapy — and nanocarrier drug delivery (dendrimers, quantum dots, carbon nanotubes, liposomes, polymeric micelles, silver nanoparticles, nanogels, nanodiamonds). The conclusions are qualitative: 3D spheroids reduce animal use and improve drug-screening predictivity. No quantitative biophysical/mechanical data and no governing equations are presented.

## 2. Problem & motivation
2D cell cultures and animal models poorly recapitulate solid-tumor architecture, microenvironment, hypoxia, and drug-resistance, causing ineffective drugs to advance into expensive in vivo testing. The review motivates 3D spheroids/organoids as more physiologically faithful, animal-sparing in vitro screening platforms for anticancer drug discovery and theragnostics.

## 3. Methods / model
Not a model or experimental paper — a literature review. The "methods" it catalogues are spheroid-fabrication and -characterization techniques used by the cited primary studies:
- **Fabrication**: spontaneous aggregation, matrix embedding, spinner flasks, ultra-low-attachment plates, micro-patterned plates, magnetic levitation, magnetic 3D printing, hanging-drop, matrix-on-top, matrix encapsulation, liquid overlay, microfluidics (Table 3, advantages/disadvantages).
- **Characterization/readout**: optical and electron microscopy, flow cytometry, Western blotting, UV/Vis spectroscopy, fluorescence spectroscopy, multi-photon microscopy, confocal laser microscopy (Figure 2).
- **Biological readouts** in cited studies: proliferation, viability, IC50, EMT markers (E-cadherin, vimentin), HIF-1/hypoxia, metabolic dependencies. Cell lines mentioned across the cancer-type tour include MDA-MB-231 (breast), 22Rv1/PC-3/C4-2B/DU145 (prostate), A549 (lung), PDAC lines, thyroid lines, and OC lines (A2780, OVCAR3, etc.). No length/time scales, no numerics, no solver — there is no model.

## 4. Key results (quantitative)
The review is overwhelmingly qualitative. The only numeric values are clinical/epidemiological or drug-screening, not biophysical:
- Lung cancer: ~1.6 million deaths/year globally, ~10% five-year survival, >80% NSCLC (Section 4.6).
- Erlotinib-loaded nanoemulsion IC50 was ~2.8x lower than free erlotinib in A549 NSCLC (Section 4.6).
- Ovarian cancer ~40% five-year survival; >80% detected at severe stages (Section 4.7).
- Qualitative spheroid-formation criterion: "cell-to-cell adhesion must be higher than cells to substrates" (Section 3) — relevant conceptually but no numbers.
- No spheroid diameters, no stiffness/modulus values, no traction, no diffusion/oxygen-gradient numbers anywhere; "scaffold stiffness" (prostate section) is mentioned only qualitatively.

## 5. Parameters & constants of interest
None (no transferable constants). The paper reports no mechanical moduli, no spheroid sizes with units, no adhesion energies, no traction stresses, no diffusion coefficients, no oxygen tensions — nothing usable as an oracle or validation anchor for a fine-grained particle simulator. All numbers present are epidemiological survival statistics or a single drug-IC50 fold-change.

## 6. Relevance to ffn_cellsim
**(f) tangential — Low.** This paper is off-topic for a mechanistic single-cell cytoskeleton+ECM simulator. It contains:
- **No cortex / FA-clutch / motor-myosin / membrane / nucleus / cytoplasm mechanism content** — none of the H.1-H.10 mechanistic units are touched.
- **No transferable parameters** for any oracle or validation gate.
- **No numerics/methods** relevant to HOOMD, BAOAB integration, or particle/bond dynamics.

The only thin connection is **spheroid-scale context**: ffn_cellsim's experimental overlay target is the PI's MCF7-spheroid-on-pV4D4/collagen-I platform, and this review (a) names MCF7-adjacent breast lines (MDA-MB-231) and (b) states the qualitative aggregation rule that cell-cell adhesion must exceed cell-substrate adhesion to form spheroids — a useful framing sentence when motivating the multi-cell cohesion overlay (the "c" term of A/A0 = a + b/R + c/R^2). But this is background narrative, not a quantitative anchor. The multicellular/spheroid scale here sits far above the single-cell, filament-resolved regime ffn_cellsim simulates; the review offers no cell-cell catch-bond cadherin parameters, no spheroid mechanics, and no ECM mechanics that could feed KU-4.2 cadherin or the FA/ECM units. Keep as a citable framing reference for "why spheroids / 3D matters" in REPORT motivation sections; do NOT mine it for physics.

## 7. Limitations & caveats
- Pure narrative review: no original data, no meta-analysis, no quantitative synthesis.
- Scale gap: spheroid/organoid/tissue-engineering scale, orders of magnitude coarser than a single-cell filament-resolved simulator; provides no bridging quantities.
- Drug-delivery / nuclear-medicine / nanocarrier emphasis is entirely outside the mechanobiology scope of ffn_cellsim.
- Figures are conceptual schematics (organoid-vs-spheroid comparison, characterization-technique overview), not data plots.

## 8. Key figures / tables
- **Figure 1** — Schematic comparison of organoid vs spheroid models (conceptual, no data).
- **Figure 2** — Overview of characterization techniques for 3D tumor spheroids (microscopy/spectroscopy/flow-cytometry; conceptual).
- **Table 3** — 3D spheroid fabrication techniques with advantages/disadvantages (hanging-drop, spinner flask, magnetic levitation/bio-printing, liquid overlay, microfluidics, etc.) — methods catalogue, qualitative.

## 9. Notable quotes / citable claims
- "The primary condition for forming cellular spheroids is that cell-to-cell adhesion must be higher than cells to substrates." (Section 3) — useful framing for multi-cell cohesion overlay.
- "Three-dimensional spheroids can correctly replicate some features of solid tumors (such as the secretion of soluble mediators, drug resistance mechanisms, gene expression patterns and physiological responses) better than 2D cell cultures or animal models." (Abstract)
- "cellular spheroids have distinct cell phenotypes that match the structure of actual tumors, such as necrotic, proliferating, and non-proliferating cells." (Section 3) — the canonical spheroid zonation, but stated without dimensions or gradients.
- "epithelial ... formed small round spheres ... whereas those with high vimentin and low E-cadherin expression levels (mesenchymal) formed large grape-like spheres" (Section 4.4, PDAC) — qualitative cadherin/EMT-to-morphology link.
