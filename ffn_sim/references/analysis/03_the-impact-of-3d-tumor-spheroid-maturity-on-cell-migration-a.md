---
id: 03_the-impact-of-3d-tumor-spheroid-maturity-on-cell-migration-a
paper_n: 3
title: "The impact of 3D tumor spheroid maturity on cell migration and invasion dynamics"
authors: "Lingke Feng et al."
year: "2025"
venue: "Biochemical Engineering Journal"
doi: "10.1016/j.bej.2024.109567"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, ECM/collagen, tangential]
entities: [du-145, mcf7, emt-6, e-cadherin, vimentin, mmp9, collagen-i, matrigel, integrin, focal-adhesion, actin, reactive-oxygen-species, tumor-spheroid]
methods: [confocal-live-imaging, rna-seq, qpcr, time-lapse-imaging, side-view-imaging, fluorescence-probe]
measurables: [invasion-speed, migration-distance, spheroid-radius, contact-angle, spheroid-thickness, gene-expression, ros-level, cell-viability]
keywords: [tumor-spheroid, spheroid-maturity, cell-migration, cell-invasion, reactive-oxygen-species, doxorubicin, du-145, breast-cancer, matrigel-invasion, necrotic-core]
tags: ["#spheroid-mechanics", "#breast-cancer", "#cell-migration", "#tumor-spheroid", "#ros", "#drug-resistance", "#tangential"]
has_transferable_params: false
---

# [3] The impact of 3D tumor spheroid maturity on cell migration and invasion dynamics

**Tags:** #spheroid-mechanics #breast-cancer #cell-migration #tumor-spheroid #ros #drug-resistance #tangential

| Field | Value |
|---|---|
| Authors | Lingke Feng et al. (corresponding: Yan Li, Ling Yu) |
| Year / Venue | 2025 / Biochemical Engineering Journal 213 (2025) 109567 |
| DOI / ID | 10.1016/j.bej.2024.109567 |
| Type | experimental |
| Pages | 10 |
| ffn_cellsim relevance | Low — purely phenomenological multicellular spheroid-scale assay; no fine-grained single-cell mechanism or transferable mechanical constant |

## 1. Summary
The authors culture three cancer cell lines — human prostate DU 145, human breast MCF-7 (non-metastatic), and murine breast EMT-6 (metastatic) — as 3D spheroids in agarose micro-wells for varying durations (3, 5/7, 11 days) and ask how spheroid "maturity"/aging affects subsequent migration (onto flat plastic) and invasion (into a Matrigel coating). The core finding is that older/more-mature spheroids migrate and invade more: e.g. DU 145 migration distance rises from ~1290 µm (3-day) to ~1496 µm (11-day) over 48 h, and invasion distance from ~1076 µm to ~1382 µm. RNA-seq + qPCR show that 10/11-day spheroids up-regulate cell-adhesion, focal-adhesion, actin-cytoskeleton, and Rap1-signalling genes plus CDH1/Vimentin/MMP-9, while a single-cell optical-fiber probe and confocal imaging show that mature spheroids have higher intra-spheroid reactive oxygen species (ROS), peaking at the necrotic core. Mature spheroids are also more resistant to doxorubicin (DOX). The take-home message is a methodological caution: spheroid culture time is a hidden confounder that must be standardized before spheroid-based migration/invasion or drug-screening assays.

## 2. Problem & motivation
2D scratch/trans-well migration assays don't reproduce the 3D architecture, hypoxic/nutrient gradients, necrotic core, or ECM interactions of solid tumors. 3D spheroids do, but the "starting condition" of a spheroid (especially culture time) is rarely controlled, and spheroids are heterogeneous (proliferative outer / quiescent middle / necrotic core zones that grow with age). The authors hypothesize that spheroid maturity mirrors physiological solid-tumor conditions (hypoxia, functional heterogeneity) and therefore systematically modulates migration and invasion. The practical motivation is reproducibility of spheroid-based functional and drug-screening assays.

## 3. Methods / model
Pure wet-lab experiment, no computational/physical model. Cell lines: DU 145, MCF-7, EMT-6 (Procell). Spheroids formed by seeding 1×10^4 cells/well into 2% agarose micro-wells (U-shaped, 2 mm height × 2 mm diameter, DLP-3D-printed mold), cultured 3/7/11 days (EMT-6 at 3/5/7 days). Assays: (i) migration — spheroid pipetted onto a flat 96-well bottom, cells escape and spread; (ii) invasion — spheroid placed on a 30 µL Matrigel (1:10 in DMEM) coating. Readouts: top-view time-lapse (ImageJ area/perimeter/distance, 0.1 µm accuracy) and a "homemade side-view device" measuring vertical thickness and the spheroid-substrate contact angle every 12 h. Molecular: RNA-seq (BGISEQ 500, DESeq2) on 2D vs 3-day vs 10-day; qPCR of CDH1, VIM, MMP-9 (GAPDH reference, 2^-ΔΔCt). ROS: DCFH-DA → DCF fluorescence read by a fiber-optic single-cell multimode analyzer penetrating the spheroid in 50 µm steps (ex 488 nm / em 525 nm) plus confocal Z-stacks. Drug: 20 µM DOX, 48 h, CCK-8 viability at 450 nm. n=6 spheroids (migration/invasion), n=3 (vertical/ROS/qPCR); Student's t-test, p<0.05. Length scales are tissue/multicellular (~150–300 µm spheroid thickness; ~1000–1500 µm migration fields); time scale is days (culture) and hours (assay).

## 4. Key results (quantitative)
- DU 145 migration distance over 48 h: 1289.67±39.15 µm (3-day), 1391.34±73.59 µm (7-day), 1495.50±99.83 µm (11-day) (Fig 1C).
- DU 145 migrated cell-occupied area after 24 h: (9.29±0.29)×10^5 µm² (3-day) → (11.99±1.07)×10^5 µm² (11-day) (Fig 1E). Migration-area perimeter 3872→4314 µm (Fig 1F).
- EMT-6 (metastatic) migration distance: 833.45±68.88 µm (3-day), 929.38±78.33 µm (5-day), 1175.11±65.75 µm (7-day) (Fig 2B).
- MCF-7 (non-metastatic) migration distance only 202.70±52.78 to 323.25±38.16 µm — much smaller; 7-day MCF-7 spheroids partly decompose (Fig 2).
- DU 145 invasion distance into Matrigel: 1075.83±40.04 µm (3-day) → 1382.13±106.99 µm (11-day); area (6.86±0.38)×10^5 → (11.36±0.88)×10^5 µm² (Fig 3C,D). Migration generally exceeds invasion.
- Vertical / contact angle: DU 145 spheroid-substrate contact angle after 48 h migration 51.6±4.64° (3-day), 66.45±2.39° (7-day), 88.83±4.28° (11-day) (Fig 4C). Initial spheroid thickness 145.47±3.02 µm (3-day), 297.90±8.60 µm (7-day), 301.73±18.70 µm (11-day); thickness reduction rate 0.43 (3-day) → 0.26 (11-day) (Fig 4D). Older spheroids spread more horizontally but thin less vertically.
- RNA-seq: 254 DEGs common to 2D-vs-3-day and 2D-vs-10-day; 68 up / 202 down in 10-day vs 3-day; KEGG enrichment in cell cycle, focal adhesion, actin cytoskeleton, Rap1 signalling; MAPK pathway down-regulated (Fig 5). qPCR: CDH1, Vimentin, MMP-9 all higher at 11 days (Fig 5E).
- Intra-spheroid ROS (photon counts, 7-day): surface (6.82±1.21)×10^4 → 50 µm (11.85±2.17)×10^4 → center (19.35±3.10)×10^4 → beyond radius (11.54±4.74)×10^4 (Fig 6A). 3-day edge/center (0.86/2.18)×10^4 vs 11-day edge/center (13.22/21.99)×10^4 — ROS rises with age and toward the core.
- DOX (20 µM, 48 h): 11-day DU 145 still migrated 537.31±34.27 µm vs 3-day 333.95±31.54 µm; invasion 524.20±34.31 (11-day) vs 311.93±33.62 µm (3-day); older spheroids more DOX-resistant by CCK-8 (Fig 7).

## 5. Parameters & constants of interest

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| DU 145 migration distance, 11-day, 48 h | 1495.50±99.83 µm | Fig 1C / p.5 |
| DU 145 invasion distance into Matrigel, 11-day | 1382.13±106.99 µm | Fig 3C / p.6 |
| Spheroid initial thickness (3/7/11-day) | 145 / 298 / 302 µm | Fig 4 / p.6 |
| Spheroid-substrate contact angle, 11-day | 88.83±4.28° | Fig 4C / p.6 |
| Effective whole-spheroid migration "speed" (inferred) | ~1300–1500 µm / 48 h ≈ 27–31 µm/h (collective front, NOT single-cell crawl) | inferred from Fig 1C/3C |

Note: these are collective multicellular front-spreading and tissue-geometry numbers, not single-cell mechanical constants. They are NOT directly usable as oracles for a fine-grained single-cell simulator. **Effectively: None (no transferable single-cell/cytoskeletal constants).**

## 6. Relevance to ffn_cellsim
This paper is **(d) spheroid/multicellular-scale context** with a faint thread of **(b) potential overlay relevance** to the project's MCF-7-spheroid validation target — but for the current mechanistic single-cell roadmap it is essentially **(f) tangential**. Why:

- It is phenomenological tissue-scale biology: spheroid front-spreading, contact angle, RNA-seq/qPCR gene panels, ROS gradients, DOX resistance. None of this exposes a per-filament, per-motor, per-clutch, or per-cross-link mechanism that `ffn_cellsim` models. There are no force-velocity curves, no bond off-rates, no moduli, no traction stresses.
- The MCF-7 line appears (relevant to the project's PI-experimental overlay target, which is MCF-7 spheroids on pV4D4/collagen-I), but here MCF-7 is the non-metastatic control with tiny migration (~200–320 µm) and partial spheroid decomposition. There is no traction-force or cohesion measurement that maps onto the project's single-cell-traction / multi-cell-cohesion overlay layers. The substrate is Matrigel/plastic, not the project's pV4D4/collagen-I.
- It does, however, supply **qualitative multicellular context** the project's spheroid-scale narrative can cite: (i) spheroids have proliferative-outer / quiescent-middle / necrotic-core zones whose collagen accumulation and ROS rise with maturity; (ii) collective invasion involves MMP-9-mediated ECM (Matrigel/collagen) degradation; (iii) focal-adhesion and actin-cytoskeleton pathways are transcriptionally up-regulated in motile mature spheroids — consistent (loosely) with the project's FA/clutch and cortex emphasis, but only as motivation, not as a parameter or oracle.
- H-unit mapping: touches ECM/collagen (Matrigel invasion, MMP-9, collagen accumulation in core) and spheroid-scale-context only at the narrative level. Nothing for cortex (H.1), FA/clutch (H.4), motor/myosin, or membrane/nucleus/cytoplasm EXTEND in a quantitative sense.

Bottom line: keep it as background/citation for the spheroid-scale and ECM-degradation framing and as a note that spheroid culture-age confounds any future multi-cell validation comparison; do not treat it as a validation oracle or parameter source. Off-topic for the fine-grained single-cell engine.

## 7. Limitations & caveats
- No mechanics: no stiffness, traction, force, or velocity-force data; "migration speed" is a collective front displacement, not single-cell motility.
- Scale gap: everything is tissue-scale (10²–10³ µm, days), orders of magnitude above the project's particle/bond resolution and the ×40 mesoscale filament coarse-graining; nothing bridges to per-filament dynamics.
- Correlational, not causal: ROS↑ and motility↑ are correlated; no perturbation establishes ROS as the driver. Gene-expression changes are descriptive.
- Confounds: spheroid age co-varies necrotic-core size, collagen accumulation, hypoxia, proliferation, and ROS — the assay cannot isolate which drives motility.
- Mixed cell lines/time points (DU 145 3/7/11 d; EMT-6 3/5/7 d) complicate cross-line comparison; n=3–6, single lab, t-test only.
- Invasion substrate is Matrigel (basement-membrane-like), not the fibrillar collagen-I the project's ECM module targets.

## 8. Key figures / tables
- Fig 1 (p.4–5): DU 145 horizontal migration vs culture time — distance/area/perimeter all rise with age (the central quantitative result).
- Fig 3 (p.5–6): Matrigel invasion distance/area for DU 145, EMT-6, MCF-7 vs culture time — invasion < migration, both increase with age.
- Fig 4 (p.6): Side-view vertical dynamics — contact angle 51.6°→88.83° and thickness-reduction rate 0.43→0.26 with age; H&E showing heterogeneous core; the most "geometry/mechanics-adjacent" figure.
- Fig 5 (p.7–8) & Fig 6 (p.8): RNA-seq/qPCR (focal-adhesion, actin, CDH1/VIM/MMP-9 up) and intra-spheroid ROS depth map (core-peaked, age-amplified).

## 9. Notable quotes / citable claims
- "11-day-old DU145 spheroids exhibited the most robust and fastest migration ability, followed by the 7-day-old spheroids, with the 3-day-old spheroids showing the slowest and weakest migration capacity." (p.6)
- "up-regulated genes were mainly involved in processes like cell cycle progression, cell adhesion, actin cytoskeleton dynamics, tight junction formation... and the Rap1 signaling pathway" / KEGG "enriched in pathways such as cell cycle, focal adhesion, and actin cytoskeleton regulation." (p.7–8)
- "the highest ROS content was detected from both older spheroids and cells migrated from the older spheroids, suggesting a correlation between the elevated ROS levels in older spheroids (11-day-old) and the increased motility of migrating DU 145 cells." (Conclusion, p.9)
- "it is crucial to systematically analyze spheroids' growth conditions, particularly the cultivation duration, before applying them to drug screening assays." (p.9) — the methodological caution that confounds any multicellular validation comparison.
