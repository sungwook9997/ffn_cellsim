---
id: 21_3d-tumor-spheroid-and-organoid-to-model-tumor-microenvironme
paper_n: 21
title: "3D Tumor Spheroid and Organoid to Model Tumor Microenvironment for Cancer Immunotherapy"
authors: "Zhu, Y.; Kang, E. et al."
year: "2022"
venue: "Organoids (MDPI)"
doi: "10.3390/organoids1020012"
paper_type: review
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, tangential]
entities: [mcf-7, pc3, hct116, e-cadherin, collagen, laminin, fibronectin, matrigel, tumor-associated-macrophage, cancer-associated-fibroblast, car-t-cell, ecm]
methods: [3d-cell-culture, spheroid-culture, organoid-culture, microfluidics, 3d-bioprinting, drug-screening, immunohistochemistry]
measurables: [oxygen-tension, ph-gradient, drug-penetration-depth, lactate-production, hypoxia]
keywords: [tumor-microenvironment, 3d-spheroid, organoid, patient-derived-organoid, hypoxia, immunotherapy, car-t, tumor-associated-macrophage, drug-screening, ecm, necrotic-core, cancer-stem-cell]
tags: ["#tumor-microenvironment", "#3d-spheroid", "#organoid", "#immunotherapy", "#spheroid-scale-context", "#review", "#biology-not-mechanics"]
has_transferable_params: true
---

# [21] 3D Tumor Spheroid and Organoid to Model Tumor Microenvironment for Cancer Immunotherapy

**Tags:** #tumor-microenvironment #3d-spheroid #organoid #immunotherapy #spheroid-scale-context #review #biology-not-mechanics

| Field | Value |
|---|---|
| Authors | Zhu, Y.; Kang, E.; Wilson, M.; Basso, T.; Chen, E.; Yu, Y.; Li, Y.-R. (UCLA) |
| Year / Venue | 2022 / Organoids (MDPI), vol. 1, 149–167 |
| DOI / ID | 10.3390/organoids1020012 |
| Type | review (biology/immunology) |
| Pages | 19 |
| ffn_cellsim relevance | Low — biology/immunology review of 3D culture methods; provides spheroid-scale TME context (hypoxia/pH gradients, ECM, cell–cell adhesion) but no mechanistic cytoskeletal/clutch/motor physics |

## 1. Summary
This is a narrative biology review surveying in-vitro 3D tumor models — cell-line-derived spheroids, patient-derived spheroids, and patient-derived organoids (PDOs) — and how they recapitulate the tumor microenvironment (TME) for cancer drug screening and immunotherapy. It catalogs (i) the cellular/acellular components of the TME (tumor cells, tumor-associated macrophages/TAMs, myeloid-derived suppressor cells/MDSCs, cancer-associated fibroblasts/CAFs, endothelial cells, ECM); (ii) culture methods (ultra-low-attachment plates, hanging drop, spinner flask, liquid overlay, microfluidics, 3D bioprinting, Matrigel/Geltrex embedding); and (iii) emergent spheroid features that mimic solid tumors — radial hypoxia/pH/necrosis gradients, the Warburg metabolic shift, increased ECM protein expression, E-cadherin–driven aggregation, cancer-stem-cell enrichment, and barriers to drug/immune-cell penetration. It concludes spheroids/organoids are becoming the gold-standard in-vitro TME model for immunotherapy (CAR-T, checkpoint-inhibitor) screening, while noting reproducibility, murine-ECM artifact, and imaging-depth limitations.

## 2. Problem & motivation
2D monolayer cultures fail to reproduce cell–cell and cell–ECM interactions, drug-resistance, and the spatial chemical gradients of real solid tumors, while animal/xenograft models are costly, low-throughput, and lack functional human immune systems. The review motivates 3D spheroid/organoid models as a middle ground that captures TME complexity (hypoxia, acidosis, immunosuppression, ECM barriers) for more predictive preclinical drug and immunotherapy screening.

## 3. Methods / model
No computational model. This is a literature review. The "methods" it discusses are wet-lab 3D-culture techniques:
- **Spheroid formation**: ultra-low-attachment 96-well plates, magnetic levitation, hanging drop, spinner flask, liquid overlay (agarose), agitation, microfluidic microwells with concentration-gradient generators.
- **Scaffolds**: ultra-porous cellulose (PC3 prostate), fibrous scaffolds (MCF-7 breast), alginate encapsulation (NCI-H157 NSCLC), collagen-coated Teflon membranes (multicellular layers).
- **Organoids**: single-cell suspensions / minced tissue embedded in ECM substitutes (Matrigel, Geltrex) per the Sato et al. protocol; air–liquid-interface variants retain native TILs/CAFs.
- **Readouts**: drug-penetration assays, pH/oxygen-gradient mapping, CSC marker assays, IF staining (DAPI/PD-L1/7-AAD), coculture immune-killing assays (CAR-T, checkpoint inhibitors).
Length scale = multicellular spheroid (hundreds of µm); the only quantitative physical scales are diffusion/penetration distances and chemical gradient ranges (below).

## 4. Key results (quantitative)
Numbers are cited from the reviewed primary literature, not original measurements:
- Doxorubicin penetrates only **~40–50 µm** from blood vessels into tumor tissue (p.151, refs 26,27).
- Chemotherapy penetration through tumor tissue can drop to as low as **~20%** due to the multicellular layer (MCL) containing laminin + collagen (p.151, ref 23).
- TME hypoxia threshold: oxygen tension **< 5–10 mmHg** drives TAM/Treg/MDSC accumulation (p.151).
- Intratumoral / spheroid-core pH can fall to **~5.6** (acidic), contributing to drug resistance (p.153).
- Only **~5%** of anticancer drugs reach the clinical stage (p.151–152, ref 42).
- Spheroid internal architecture: outer proliferative layer / middle senescent layer / necrotic core (p.155).
- Spheroids show **increased ECM protein expression** (fibronectin, laminin, collagens) vs. 2D — glioma U-118 MG and thyroid HTh-7 lines (p.154, ref 95).
- Confocal microscopy has limited penetration depth, preventing imaging of large spheroids (p.161).

## 5. Parameters & constants of interest
The only transferable values are coarse spheroid-scale TME boundary conditions — useful as soft context for an EXTEND/multicellular overlay, NOT as mechanistic single-cell constants.

| Quantity | Value + units | Source in paper |
|---|---|---|
| Drug (doxorubicin) penetration depth from vessel | 40–50 µm | p.151 (refs 26,27) |
| Multicellular-layer drug penetration fraction | as low as ~20% | p.151 (ref 23) |
| Hypoxia threshold (TME) | < 5–10 mmHg O₂ | p.151 |
| Spheroid/tumor core pH (acidic floor) | ~5.6 | p.153 (refs 65,72,73) |
| Spheroid ECM proteins enriched vs 2D | fibronectin, laminin, collagen (qualitative) | p.154 (ref 95) |
| Spheroid radial structure | proliferative shell / senescent mid / necrotic core | p.155 |

No elastic moduli, traction stresses, bond rates, motor velocities, filament parameters, or any constant usable as a fine-grained-simulator oracle.

## 6. Relevance to ffn_cellsim
**Category: (d) spheroid/multicellular-scale context, mostly tangential.** This is a biology/immunology review with no mechanics, no cytoskeleton, no continuum or particle model, and no transferable mechanistic constants. It does NOT inform the H.1→H.7 mechanistic single-cell chain (cortex, FA/clutch, motor/myosin, ECM cross-links, integrator/numerics) in any way.

Its only marginal value is as **background context for the multicellular / spheroid-scale overlay** that sits above the single-cell simulator:
- Reinforces that the project's MCF7-spheroid-on-collagen-I validation target lives in a TME with radial O₂/pH/necrosis gradients and ECM (collagen/laminin/fibronectin) enrichment — relevant when reasoning about why a multicellular cohesion overlay (the "c" coefficient in A/A₀ = a + b/R + c/R²) behaves differently from single-cell traction (the "b" coefficient).
- Confirms **E-cadherin–driven cell–cell aggregation** as the dominant spheroid-cohesion adhesion (p.154, ref 91) — consistent with the project's catch-bond cadherin (KU-4.2) mechanism, though this paper gives no kinetic numbers.
- Provides order-of-magnitude diffusion/penetration scales (40–50 µm, hypoxia <5–10 mmHg) that could seed a future reaction-diffusion oxygen/nutrient field IF the project ever adds spheroid-scale biochemistry — far outside the current Phase 1 mechanistic scope.

For a fine-grained particle simulator of a single cell's cytoskeleton + adhesions, this paper is off-topic: it operates at the tissue/organ scale, discusses immunotherapy pharmacology, and lists no force, stiffness, rate, or length constant at the molecular/filament scale. Keep it filed as multicellular-context background only. has_transferable_params is set true ONLY because of the coarse TME boundary-condition numbers (penetration depth, hypoxia threshold, core pH), which are context anchors, not mechanistic oracles.

## 7. Limitations & caveats
- Pure narrative review — no original data, no model, no equations; all numbers are secondhand from cited primaries.
- Scale gap is total: spheroid/organoid (hundreds of µm, many cells) vs. ffn_cellsim's single-cell ×40-mesoscopic filament resolution. Nothing here constrains a particle/bond model.
- TME numbers are coarse ranges (O₂ <5–10 mmHg, pH ~5.6, 40–50 µm) without spatial-resolution or per-cell detail.
- Heavy immunology/pharmacology focus (CAR-T, checkpoint inhibitors, TAM/MDSC polarization) is orthogonal to mechanobiology.

## 8. Key figures / tables
- **Figure 1 (p.150–151)**: schematic of the immunosuppressive TME — tumor cells with TAMs, MDSCs, CAFs, TILs (panel a); hypoxia as a key TME driver of immunosuppression and proliferation (panel b). Context for what surrounds a spheroid.
- **Figure 2 (p.158–159)**: patient-derived spheroid/organoid formation workflow (a); IF staining of clear-cell renal carcinoma organoid showing retained TILs/TAM (b); chordoma organoid + nivolumab enhanced TIL killing (c). Method/immunotherapy illustration.
- **Table 1 (p.160)**: advantages/disadvantages of cell-line spheroids vs patient-derived spheroids vs organoids (cost, ECM definition, CSC enrichment, reproducibility). The single most compact summary of the review.

## 9. Notable quotes / citable claims
- "doxorubicin … is unable to penetrate more than 40–50 µm from blood vessels to sufficiently target tumor cells" (p.151).
- "hypoxic conditions of lower than 5–10 mmHg in the TME is responsible for TAM, Treg, and MDSC accumulation" (p.151).
- "Acidic conditions inside the TME, which may subside to as low as a pH of 5.6, contribute to drug resistance" (p.153).
- "3D spheroids show increased expression of ECM proteins, including fibronectin, laminin, and collagens" (p.154, on glioma U-118 MG / thyroid HTh-7 lines).
- "E-cadherin … would increase its expression in the multicellular spheroids of ovarian cancer cells … induces spheroid formation, maintenance, and drug resistance" (p.154).
