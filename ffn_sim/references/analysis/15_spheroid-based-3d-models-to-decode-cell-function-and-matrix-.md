---
id: 15_spheroid-based-3d-models-to-decode-cell-function-and-matrix-
paper_n: 15
title: "Spheroid-Based 3D Models to Decode Cell Function and Matrix Effectors in Breast Cancer"
authors: "Mangani, S. et al."
year: "2025"
venue: "Cancers (MDPI), 17, 3512"
doi: "10.3390/cancers17213512"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, ECM/collagen, tangential]
entities: [mcf7, mda-mb-231, e-cadherin, vimentin, slug, jam-a, syndecan-1, syndecan-4, mmp2, mmp7, mmp9, mt1-mmp, egfr, igf1r, estrogen-receptor-alpha, estrogen-receptor-beta, integrin]
methods: [scanning-electron-microscopy, rt-qpcr, immunofluorescence, confocal-live-imaging, wound-healing-assay, spheroid-dissemination-assay, kaplan-meier-analysis, string-network]
measurables: [migration-speed, wound-closure-rate, spheroid-core-area, dissemination-area, gene-fold-change, overall-survival]
keywords: [breast-cancer, spheroid, mcf7, mda-mb-231, emt, matrix-metalloproteinase, syndecan, e-cadherin, tumor-microenvironment, matrix-free-3d-culture, collective-migration]
tags: ["#spheroid-biology", "#breast-cancer", "#emt", "#ecm-effectors", "#mmp", "#syndecan", "#multicellular-context", "#tangential"]
has_transferable_params: false
---

# [15] Spheroid-Based 3D Models to Decode Cell Function and Matrix Effectors in Breast Cancer

**Tags:** #spheroid-biology #breast-cancer #emt #ecm-effectors #mmp #syndecan #multicellular-context #tangential

| Field | Value |
|---|---|
| Authors | Mangani, S. et al. (Koutsakis, Koletsis, Piperigkou, Franchi, Götte, Karamanos) |
| Year / Venue | 2025 / Cancers (MDPI), 17, 3512 |
| DOI / ID | 10.3390/cancers17213512 |
| Type | experimental (molecular/cell biology) |
| Pages | 23 |
| ffn_cellsim relevance | Low — molecular-expression/phenotype study of the same MCF-7/MDA-MB-231 spheroid system the project validates against, but carries no mechanistic/biophysical parameters |

## 1. Summary
The authors develop matrix-free 3D spheroids from two breast cancer cell lines of distinct ER status and metastatic potential (MCF-7: Luminal A, ERα+, low-metastatic, epithelial; MDA-MB-231: TNBC, ERα-/ERβ+, high-metastatic, mesenchymal-like) by self-aggregation in ultra-low-adhesion U-bottom 96-well plates (72 h). They characterize the spheroids morphologically (phase-contrast + SEM), measure gene expression of EMT markers, receptors (ERs, EGFR, IGF1R), syndecans (SDC1/4) and MMPs (MMP2/7/9/14) by RT-qPCR comparing 2D vs 3D, run immunofluorescence/confocal for E-cadherin protein, and assess function via spheroid dissemination and wound-healing migration assays. Key biology: MDA-MB-231 undergoes a partial mesenchymal→epithelial shift in 3D (E-cadherin/CDH1 up 3.8-fold, vimentin down to 0.1-fold), all MMPs are upregulated in 3D, and spheroid-derived cells migrate faster than 2D monolayers. Bioinformatic tools (STRING, Kaplan–Meier, Human Protein Atlas) tie these matrix regulators to patient survival.

## 2. Problem & motivation
2D monolayer cultures fail to reproduce the cell–cell and cell–matrix interactions of solid tumors (uniform nutrient/O2 access, no avascular gradients, no 3D architecture). The paper asks whether a simple, scaffold-free spheroid platform recapitulates tumor-microenvironment features and whether ECM effectors (MMPs, syndecans) and receptors are differentially regulated in 3D vs 2D in two breast cancer subtypes — to provide a more predictive preclinical model.

## 3. Methods / model
Purely experimental, no physical/computational model.
- **Cell lines:** MCF-7 (HTB-22) and MDA-MB-231 (HTB-26), DMEM + 10% FBS, 37 °C, 5% CO2.
- **Spheroids:** 5,000–15,000 cells/well in ultra-low-adhesion U-bottom 96-well plates, self-aggregation, 72 h, then 16–20 h serum starvation. Scaffold/matrix-free.
- **Morphology:** phase-contrast microscopy; SEM (Karnovsky fixation, OsO4 post-fix, critical-point drying, Pd/Au sputter, Philips 515 SEM, secondary-electron mode).
- **Expression:** RT-qPCR, 2^(-ΔΔCt) normalized to ACTB (primer Table 1: ESR1/2, EGFR, IGF1R, CDH1, F11R, VIM, SNAI2, SDC1/4, MMP2/7/9/14). E-cadherin immunofluorescence (Alexa-488) + DAPI; confocal (Leica TCS SP8) on spheroids.
- **Function:** spheroid dissemination on flat-bottom plates (core area + dissemination area via ImageJ, 0/24/48 h); wound-healing/scratch assay with cytarabine (10 µM) to suppress proliferation, quantifying wound area at 0/24/48 h.
- **Bioinformatics:** Kaplan–Meier plotter (OS, Luminal A n=141, TNBC n=201), STRING PPI network, Human Protein Atlas.
- **Scales:** spheroids form over 72 h; intercellular gaps 5–10 µm (SEM); assay readouts in hours/days. No force, modulus, or rheology measurement.

## 4. Key results (quantitative)
- MCF-7 spheroids: dense, compact, intercellular gaps 5–10 µm wide (Fig 1C,D, SEM). MDA-MB-231 spheroids: loosely organized outer layer with detaching/migrating surface cells, compact core, mesenchymal→globular epithelial-like shift (Fig 1F,H).
- EMT shift (MDA-MB-231 3D vs 2D, RT-qPCR, Fig 2A): E-cadherin/CDH1 +3.8-fold (p≤0.001); JAM-A/F11R +1.1-fold (p≤0.0001); SNAI2/SLUG down to 0.3-fold (p≤0.001); VIM down to 0.1-fold (p≤0.05). 3D E-cadherin still below epithelial 2D MCF-7 levels (Fig 2B).
- Receptors (3D vs 2D, Fig 3A): ESR1 down to 0.7-fold in MCF-7 (p≤0.001); ESR2 +4.5-fold in MDA-MB-231 (p≤0.01); IGF1R up in both (MCF-7 0.5-fold change, MDA 2-fold; p≤0.01); EGFR +9-fold in MCF-7 (p≤0.001) but down to 0.3-fold in MDA-MB-231 (p≤0.01).
- Syndecans (Fig 3B): SDC1 down to 0.2-fold in MDA-MB-231 (p≤0.05); SDC4 down to 0.25-fold in MCF-7 (p≤0.01) but up 0.5-fold change in MDA-MB-231 (p≤0.05).
- MMPs all up in 3D (Fig 3C): MMP2 +0.4-fold MCF-7 / +10-fold MDA (p≤0.01); MMP9 ~+1-fold both (p≤0.01); MMP7 +11-fold MCF-7 (p≤0.001) / +1.6-fold MDA (p≤0.01); MMP14/MT1-MMP +27-fold MCF-7 / +1.5-fold MDA (p≤0.01).
- Migration (wound-healing, Fig 6): 2D MCF-7 ~5%/10% closure at 24/48 h vs spheroid-derived MCF-7 ~30%/40%; 2D MDA ~30%/45% vs spheroid-derived MDA ~35%/55%. Spheroid-derived cells migrate faster than 2D for both lines.
- Dissemination (Fig 5): core area shrinks and dissemination area expands 24→48 h; MDA-MB-231 disseminates faster than MCF-7.
- Survival (Kaplan–Meier, Fig 4): low ESR1 → lower OS in Luminal A (p≤0.05); high ESR2 → poorer OS in TNBC (p≤0.01); high SDC4 → poorer OS in TNBC (p≤0.001); high MMP2/MMP9 → poorer OS in TNBC (p≤0.01/p≤0.001); MMP14 trend in both subtypes (p=0.065 Luminal A, p=0.053 TNBC).
- THPA reference values: IGF1R high in MCF-7 (nTPM 52.8); EGFR predominant in MDA-MB-231 (nTPM 61.6).

## 5. Parameters & constants of interest
None transferable to a mechanistic single-cell cytoskeleton+ECM simulator. The numbers here are gene-expression fold-changes, qualitative morphology, and assay percentages — no elastic moduli, traction stresses, bond rates, motor velocities, or geometric/mechanical constants usable as an oracle.

| Quantity | Value+units | Source-in-paper |
|---|---|---|
| Intercellular gap width (MCF-7 spheroid, SEM) | 5–10 µm | Fig 1C,D / §3.1 |
| Spheroid seeding density | 5,000–15,000 cells/well | §2.2 |
| Spheroid maturation time | 72 h | §2.2 |
| Wound-closure % (MCF-7 spheroid-derived, 24/48 h) | ~30% / ~40% | Fig 6C / §3.6 |
| Wound-closure % (MDA-MB-231 spheroid-derived, 24/48 h) | ~35% / ~55% | Fig 6D / §3.6 |

These are at best loose multicellular-scale context anchors (gap width, time scales), not physical constants. Marked here only because they are the few dimensional numbers in the paper.

## 6. Relevance to ffn_cellsim
**Category: (d) spheroid/multicellular-scale context, and largely (f) tangential.** This is a molecular/phenotypic cell-biology study, not a biophysical or modeling one. It contains no quantity a fine-grained particle/bond single-cell simulator can ingest — no cortical tension, traction stress, adhesion bond rate, motor force-velocity, ECM stiffness, or cytoskeletal parameter.

Where it touches the project:
- **Same experimental system as the validation target.** ffn_cellsim's overlay target is MCF-7-spheroid-on-pV4D4/collagen-I (single-cell traction, multi-cell cohesion). This paper studies the same two cell lines (MCF-7, MDA-MB-231) as spheroids, so it is useful *context* for what those lines do phenotypically in 3D — e.g. MDA-MB-231's mesenchymal→epithelial shift, E-cadherin upregulation in 3D, faster spheroid-derived migration. This informs expectations for cohesion (cadherin) and dissemination at the multi-cell scale, which maps loosely to the project's later multi-cell/cohesion overlay and the catch-bond cadherin (KU-4.2) mechanism — but only qualitatively.
- **ECM/collagen theme is nominal.** Despite "matrix effectors" in the title, the platform is *matrix-free* (scaffold-free self-aggregation); ECM appears only as gene-expression of MMPs/syndecans, not as a mechanical/structural element. No collagen mechanics, no fiber data — so it does not inform the project's explicit ECM cross-link/collagen-I particle model.
- **Not a validation oracle, not a parameter source, not a mechanism reference, not numerics.** It cannot anchor a sanity gate or supply a constant.

Plainly: off-topic for the runtime mechanics of a mechanistic single-cell simulator. Keep as background reading for the cancer-cell-type/EMT/cadherin/MMP biology of the MCF-7 & MDA-MB-231 lines the project overlays against, and as a citable anchor that E-cadherin-mediated cohesion and migratory phenotype differ sharply between these two lines in 3D. Relevance: Low.

## 7. Limitations & caveats
- Scaffold/matrix-free spheroids: no defined ECM, no controlled stiffness — opposite of a mechanistic ECM model; "matrix effectors" measured as mRNA only.
- Read-out is gene expression (RT-qPCR fold-change), with limited protein validation (E-cadherin only); no mechanical, force, or rheological measurement at any scale.
- Cancer-cell-only spheroids (no stromal/fibroblast compartment), which the authors note may itself alter SDC1 expression — so even the biology is not a full TME.
- Fold-changes reported relative to 2D with n=3 biological replicates; several "trends" not statistically significant (e.g. MMP14 survival).
- Scale gap vs ffn_cellsim is total: this is a multicellular tissue-biology assay, not single-cell particle/bond physics.

## 8. Key figures / tables
- **Figure 1** — Phase-contrast + SEM morphology of MCF-7 (compact, 5–10 µm gaps) vs MDA-MB-231 (loose outer layer, detaching cells) spheroids; the clearest physical/morphological content.
- **Figure 3 (A–C + D)** — RT-qPCR fold-changes for ERs/RTKs, MMPs, syndecans (2D vs 3D) plus the STRING PPI network linking receptors and matrix effectors.
- **Figure 6** — Wound-healing migration: spheroid-derived cells migrate faster than 2D for both lines (the only quasi-kinetic functional readout).
- **Table 1** — qPCR primer panel (the gene list studied: EMT markers, ERs, EGFR/IGF1R, SDC1/4, MMP2/7/9/14).

## 9. Notable quotes / citable claims
- "a statistically significant increase in the epithelial markers E-cadherin/CDH1 (3.8-fold change, p ≤ 0.001) ... was observed in MDA-MB-231 spheroids as compared to their levels in 2D cultures" (§3.2).
- "MCF-7 spheroids showed compacted cells which were individually indistinguishable, and developed intercellular gaps (5–10 µm wide), potentially to facilitate nutrient/O2 diffusion" (§3.1, Fig 1).
- "spheroid-derived cells demonstrated significantly faster migration compared to their 2D counterparts of both cell lines" (§3.6, Fig 6).
- "during tumor spheroid formation, cells initially aggregate through loose integrin-ECM interactions, while over time, they establish tighter cell–cell contacts, primarily mediated by increased E-cadherin expression and reduced N-cadherin levels" (§4, citing ref [37]) — relevant qualitative picture for cadherin-mediated cohesion.
