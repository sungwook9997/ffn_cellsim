---
id: 07_in-vitro-3d-spheroid-model-preserves-tumor-microenvironment-
paper_n: 7
title: "In Vitro 3D Spheroid Model Preserves Tumor Microenvironment of Hot and Cold Breast Cancer Subtypes"
authors: "Hemavathi Dhandapani et al."
year: "2023"
venue: "Advanced Healthcare Materials"
doi: "10.1002/adhm.202300164"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, tangential]
entities: [mcf7, mda-mb-231, e-cadherin, n-cadherin, mmp-9, vimentin, hif-1a, fibroblast, cancer-associated-fibroblast, macrophage, monocyte, t-cell, ecm]
methods: [confocal-live-imaging, flow-cytometry, qrt-pcr, scanning-electron-microscopy, elisa, liquid-overlay-spheroid, mtt-assay]
measurables: [spheroid-radius, ic50, gene-expression-fold-change, cell-cycle-fraction, immune-infiltration]
keywords: [breast-cancer, tumor-microenvironment, 3d-spheroid, liquid-overlay, EMT, cancer-stemness, immunosuppression, cancer-associated-fibroblast, tumor-associated-macrophage, MCF-7, MDA-MB-231, doxorubicin]
tags: ["#spheroid-tumor-biology", "#breast-cancer", "#tumor-microenvironment", "#immuno-oncology", "#tangential", "#mcf7", "#mda-mb-231", "#cell-line-context"]
has_transferable_params: false
---

# [7] In Vitro 3D Spheroid Model Preserves Tumor Microenvironment of Hot and Cold Breast Cancer Subtypes

**Tags:** #spheroid-tumor-biology #breast-cancer #tumor-microenvironment #immuno-oncology #tangential #mcf7 #mda-mb-231 #cell-line-context

| Field | Value |
|---|---|
| Authors | Hemavathi Dhandapani, Armaan Siddiqui, Shivam Karadkar, Prakriti Tayalia (IIT Bombay) |
| Year / Venue | 2023 / Advanced Healthcare Materials, 12, 2300164 |
| DOI / ID | 10.1002/adhm.202300164 |
| Type | experimental (cancer cell biology / immuno-oncology) |
| Pages | 13 |
| ffn_cellsim relevance | Low — pure tumor-immunology/gene-expression biology; uses the same MCF-7 cell line as the project's validation overlay but contributes no mechanics, no physical constants, no cytoskeleton/ECM model |

## 1. Summary
The authors build a 3D multicellular spheroid model of breast tumor microenvironment (TME) via the liquid-overlay (agarose-coated well) technique, comparing a "hot" triple-negative subtype (MDA-MB-231) against a "cold" luminal-A subtype (MCF-7). Mono-, co- (cancer + human dermal fibroblasts, HDF), and tri-culture (+ THP-1 monocytes or PBMCs) spheroids are characterized for hypoxia (HIF-1α), epithelial-to-mesenchymal transition (E-/N-cadherin, MMP-9, vimentin), cancer stemness (OCT3/4, NANOG, ABCB6), chemoresistance to doxorubicin (DOX), immunomodulator expression (PDL-1, IDO-1, TGF-β, MHC I/II), CAF conversion (CXCL12, FAP, FSP-1, α-SMA), and immune-cell infiltration/polarization (M2 macrophages CD68/CD206, Tregs FoxP3). MDA spheroids show a more mesenchymal, stem, immunosuppressive phenotype with earlier immune infiltration; both subtypes show suppressive TME features that are partly reversed by the IDO-1 inhibitor 1-methyl-tryptophan (1-MT), most strongly in MCF-7 tricultures. The conclusion is that the spheroid model recapitulates in-vivo subtype-specific TME for drug/immunomodulator screening.

## 2. Problem & motivation
Standard 2D culture and animal models fail to recapitulate the human TME (cell–cell contact, ECM crosstalk, immune cells) that drives therapy resistance and immune evasion in breast cancer. The question: can a tractable in-vitro 3D spheroid platform preserve the distinct microenvironments of "hot" (TNBC) versus "cold" (luminal) subtypes well enough to screen immunomodulatory and chemo-therapeutics? Motivation is translational/immuno-oncology drug screening, not biophysics.

## 3. Methods / model
Purely experimental; no computational/physical model. Cell lines: MCF-7 (luminal), MDA-MB-231 (TNBC), HDF (stromal fibroblasts), THP-1 (monocytes), and donor PBMCs. Spheroids formed by liquid-overlay on 1.6% ultrapure agarose in 96-well plates, cells centrifuged for contact, grown ≤11 days. Assays/readouts: bright-field/ImageJ size tracking; calcein/PI live–dead confocal (LSM 780); SEM (Phenom ProX) for surface/ECM morphology; qRT-PCR (ΔΔCt vs GAPDH) for EMT/stemness/CAF/immunomodulator transcripts; flow cytometry (BD FACSAria) for protein markers and cell-cycle (PI); MTT IC50 of DOX in 2D vs 3D; ELISA for IL-10, CCL2, kynurenine. Scales are tissue/cell-population (microns–millimeters, days); no sub-cellular mechanics, no forces, no time-resolved dynamics relevant to a particle simulator.

## 4. Key results (quantitative)
- Spheroid seeding densities: MCF-7 2000–10000 cells/well; MDA-MB-231 500–4000 cells/well (Experimental Section, p.11). Coculture ratio HDF:cancer 1:1 (MCF 3000:3000; MDA 1500:1500); triculture 1:1:1 (Section 4).
- MCF spheroids compact until day 9, grow to day 11 only for 8000/10000-cell seeds; MDA spheroids grow to day 9 then disintegrate by day 11 (Fig S2). Size increases with seeding density.
- Viability: healthy to day 3, necrotic-core/PI-positive death by day 6 (Fig 1a).
- HIF-1α (hypoxia) absent until day 3, marked by day 6 in MCF; gradual day-wise rise in MDA (Fig 1b).
- EMT: E-cadherin up in MCF, significantly down in MDA by day 6; N-cadherin and MMP-9 up in MDA, reversed in MCF (Fig 1c–e); vimentin up in MDA (Fig S4).
- Cell cycle: 3D spheroids have more quiescent G1/G0 and less G2/M than 2D (Fig 1f–i).
- DOX IC50 significantly higher in 3D than 2D (Fig S4d,e); 2D dose range 97.6–25000 ×10⁻⁹ M (Section 4). DOX penetration increases over ~4 h (Fig S5).
- 1-MT dose used: 0.5 ×10⁻⁹ M; downregulates CD68/CD206 in MCF tricultures (Fig 6d,e).
- cDNA from 780 ng RNA; calcein 1 µg mL⁻¹, PI 1 mg mL⁻¹; cell-cycle PI 20 µg mL⁻¹ + RNase 10 µg mL⁻¹ (Section 4) — all assay reagent concentrations, not physical constants.
- Immune: significant THP-1 infiltration into MDA tricultures by 48 h (earlier than MCF) (Fig 6b,c); M2 (CD68⁺/CD206⁺) and Treg (CD3⁺CD4⁺FoxP3⁺) elevated in MDA; strong CD206–IL-10 correlation, none with kynurenine (Fig 8c).

## 5. Parameters & constants of interest
None (no transferable physical constants). All quantitative values are gene-expression fold-changes, flow-cytometry marker percentages, drug IC50/dose concentrations, seeding densities, and reagent concentrations — none are mechanical, kinetic, or material parameters usable as an oracle or anchor for a fine-grained cytoskeleton+ECM particle simulator. There are no elastic moduli, traction stresses, bond rates, viscosities, diffusion coefficients, or geometric mechanical descriptors.

## 6. Relevance to ffn_cellsim  <-- MOST IMPORTANT
This is **(f) tangential**, with a faint **(d) spheroid/multicellular-scale context** flavor. The paper sits squarely in immuno-oncology cancer biology: it measures transcriptional/protein phenotype, immune infiltration, and drug response in 3D spheroids. ffn_cellsim is a mechanistic single-cell mechanobiology simulator (explicit filaments, motors, clutches, cross-links, Bell-Evans/Hill/Stam-Hocky physics) — none of which appears here. There is **no governing physics, no force/length/time scale, no parameter** that maps to any H-unit (cortex H.1, FA/clutch, motor/myosin, ECM/collagen, membrane/nucleus/cytoplasm EXTEND, numerics).

The only thread to the project: it uses **MCF-7 and MDA-MB-231**, the same breast-cancer lines as the project's PI experimental validation overlay (MCF-7-spheroid-on-pV4D4/collagen-I), and it characterizes spheroid-scale phenomena (compaction vs disintegration, EMT, hypoxic core, ECM deposition seen by SEM). That makes it weak background/context for the *multi-cell cohesion* end of the validation roadmap — e.g., qualitative expectations that MCF-7 forms cohesive E-cadherin-rich compact spheroids while MDA-MB-231 is mesenchymal/disintegrating. But it is **not** a validation oracle (no quantitative mechanical target), **not** a parameter source, **not** a mechanism reference, and **not** numerics. It should be filed as cell-line phenotype context only; do not mine it for simulator constants.

## 7. Limitations & caveats
- Bulk-population readouts (qRT-PCR, flow cytometry, ELISA): no single-cell spatial mechanics, no force/displacement, no time-resolved cytoskeletal dynamics.
- Spheroid "size" is a 2D-projected ImageJ area, not a calibrated mechanical/geometric measurement; no modulus or traction data.
- Phenotype is descriptive/correlative (e.g., CD206–IL-10 correlation); no causal mechanistic model.
- Cell-line and donor-PBMC variability; n = 3–6 typical; in-vitro model, not validated against matched patient tissue here.
- Scale gap vs ffn_cellsim is total: this is tissue/transcriptomics scale, the simulator is nanometre–micron particle/bond scale.

## 8. Key figures / tables
- Figure 1 — spheroid viability (calcein/PI), HIF-1α hypoxia, EMT markers (E-/N-cadherin, MMP-9), and 2D-vs-3D cell-cycle distribution for both subtypes.
- Figure 3 — immunosuppressive immunomodulators (PDL-1, IDO-1, TGF-β) in 2D vs 3D ± DOX, with regulation schematic (Fig 3e).
- Figure 5 — CAF-marker activation (CXCL12, FSP-1, FAP) and spatial HDF distribution when fibroblasts are cocultured (HDF central in MCF, dispersed in MDA).
- Figure 6 — confocal infiltration of CMTPX-labeled THP-1 monocytes and their CD68/CD206 M2 differentiation in tricultures ± 1-MT.

## 9. Notable quotes / citable claims
- "an in vitro 3D spheroid model is designed using liquid overlay method to simulate hot (MDA-MB-231) and cold (MCF-7) breast tumor microenvironment (TME)" (Abstract, p.1).
- "MCF spheroids showed initial compaction until day 9 ... whereas MDA spheroids gradually increased in size until day 9 and eventually disintegrated by day 11" (Section 2.1, p.2).
- "E-cadherin was significantly upregulated ... in formation of MCF spheroids via tight cell–cell adhesions ... on the contrary, expression of E-cadherin was significantly downregulated in MDA spheroids by day 6 indicating their mesenchymal and metastatic phenotype" (Section 2.1, p.2).
- "The concentration of drug (IC50) required to cause cell death in 3D spheroids was found to be significantly higher as compared to that in 2D culture" (Section 2.2, p.3).
