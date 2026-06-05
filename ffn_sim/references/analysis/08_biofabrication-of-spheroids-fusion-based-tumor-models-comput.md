---
id: 08_biofabrication-of-spheroids-fusion-based-tumor-models-comput
paper_n: 8
title: "Biofabrication of spheroids fusion-based tumor models: computational simulation of glucose effects"
authors: "David J Bustamante et al."
year: "2021"
venue: "Biofabrication"
doi: "10.1088/1758-5090/abe025"
paper_type: agent-based
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, tangential]
entities: [cadherin, integrin, fibronectin, spheroid, tumor-cell, mcf7, glucose]
methods: [agent-based, cellular-potts, monte-carlo, reaction-diffusion]
measurables: [spheroid-radius, proliferation-rate, invasion-speed, oxygen-tension, cohesion]
keywords: [cellular-potts, compucell3d, glazier-graner-hogeweg, spheroid-fusion, glucose-metabolism, tumor-microenvironment, scaffold-free-bioprinting, kenzan, diabetes, cell-sorting, center-of-mass]
tags: ["#spheroid-fusion", "#cellular-potts", "#agent-based", "#tumor-model", "#glucose-metabolism", "#tangential", "#multicellular-scale"]
has_transferable_params: false
---

# [8] Biofabrication of spheroids fusion-based tumor models: computational simulation of glucose effects

**Tags:** #spheroid-fusion #cellular-potts #agent-based #tumor-model #glucose-metabolism #tangential #multicellular-scale

| Field | Value |
|---|---|
| Authors | David J Bustamante et al. (E J Basile, B M Hildreth, N W Browning, S A Jensen, L Moldovan, H I Petrache, N I Moldovan) |
| Year / Venue | 2021 / Biofabrication 13 035010 |
| DOI / ID | 10.1088/1758-5090/abe025 |
| Type | agent-based (Cellular Potts / GGH) |
| Pages | 15 |
| ffn_cellsim relevance | Low — multicellular-tissue-scale Potts model with lumped cells; no subcellular cytoskeleton/clutch/motor mechanism |

## 1. Summary
The authors use the open-source CompuCell3D (CC3D) platform — a Cellular Potts / Glazier–Graner–Hogeweg (GGH) agent-based modeling framework — to simulate the **fusion of cell spheroids** (tumor and "normal" cells) as a function of **glucose availability**. They extend a published single-tumor-spheroid CC3D model (Swat et al 2015) by reverse-engineering "normal" cells (muting proliferation/stemness) and modeling binary doublet fusions and a 3×3 microenvironment array (a tumor nodule surrounded by 8 normal spheroids). They introduce a shape-independent fusion metric, the **Distance between Centers of Mass (DCM)**. Heuristic findings: at hypoglycemic glucose, tumor spheroids abruptly collapse (amplified by contact with normal spheroids); at hyperglycemic glucose, cancerous cells intermix more, show strong anti-phase proliferating/quiescent oscillations, structural re-fragmentation, and escape/invasion of the surrounding normal "box." The work is explicitly framed as generating testable hypotheses about cancer–nutrition–diabetes links and informing scaffold-free bioprinting ("Kenzan" method).

## 2. Problem & motivation
In vitro tumor spheroid models are limited by nutrient (glucose) diffusion into deeper layers when non-vascularized, and scaffold-free biofabrication (cells-only bioprinting) creates larger constructs that further expose cells to nutrient deprivation. The question: how does **glucose availability** shape spheroid fusion dynamics, intermixing, structural stability, and tumor-cell invasion in single- and multi-spheroid constructs? Motivation is biofabrication/tissue-engineering design plus the medical link between diabetes/diet and cancer progression.

## 3. Methods / model
- **Model class**: Cellular Potts / GGH agent-based model on CompuCell3D. Cells are flat multi-pixel domains (effectively 2D disks despite the "3D" name); dynamics via Monte Carlo (Metropolis) energy minimization driven by adhesion (cadherin–cadherin, integrin–fibronectin) plus volume constraints.
- **Cell types**: proliferative, quiescent, necrotic, proliferative stem, quiescent stem, plus a "medium" type representing stroma/ECM. Stem types were disabled; "normal" cells were created by muting cancer proliferation/stemness parameters.
- **Glucose field**: steady-state reaction–diffusion field with diffusion constant, decay constant, secretion constant, and cell uptake rate. Glucose modeled as a "production/generation rate"; diffusion constant varied 0.1× to 10× of physiologic. Each MCS, a cell's center glucose is checked vs a threshold; below-threshold accumulates a "damage point," and on enough damage the cell turns irreversibly necrotic.
- **Time/length scale**: time is in Monte Carlo Steps (MCS); calibration ≈ 20 s/MCS (taken from Swat et al), so 5000 MCS ≈ 24 h, ~4000 MCS ≈ 15 h. Cells ~3 pixels wide; spheroids "muted" at 15 pixels radius; doublet centers placed 32 pixels apart. Runs up to 30 000 MCS.
- **Metrics**: DCM (Eq 1, Euclidean distance between spheroid centers of mass); Feret-style doublet length/width (max/min contiguous non-medium pixel runs via a custom Python applet on AVI frames); aspect ratio; cell-type subpopulation counts; dispersal fraction of tumor cells outside a dynamically sized circle.
- **Controls/scans**: triplicate runs; cadherin parameter scan for normal cells (values 9, 12, 15) vs tumor cadherin = 8; t-test for normal-vs-tumor significance.

## 4. Key results (quantitative)
- Normal-spheroid doublet length decreased quasi-linearly after contact, ~**15% reduction** over the considered interval (Fig 4A, lower curve); width increased correspondingly (Fig 4B).
- Normal vs tumor length/width curves statistically different (**p < 0.05**) except first two time points (Fig 4 caption); width relative increase ~**2× larger** than length (p. 8 / Fig S4B).
- Glucose scanned over **0.1×, 0.2×, 0.5×, 1×, 5×, 10×** of normal; 5× and 10× roughly mapped to type II and type I diabetes respectively (p. 12). Diffusion constant varied 0.1×–10× (p. 3).
- "Metabolic shock" from fusion itself drops nutrient supply as early as **50 MCS**; a cooperative DCM collapse appears at ~**400 MCS** in low glucose (p. 5–6, Fig S2).
- Tumoral fusions stabilize at ~**4000 MCS** (≈15 h) before a second disintegration phase (p. 11), in agreement with empirical data of Susienka et al.
- Fusion progress insensitive to glucose in early phase (up to ~**500 MCS**), then proliferation/turnover dominates (p. 12).
- Microenvironment model: at basal/low glucose the central tumor disappears; at hyperglycemia it stabilizes and tumor cells escape the normal "box" even when normal cadherin (9) > tumor cadherin (8) (Figs 5–7).
- Initial spheroid cell counts varied **76–81** cells despite identical initial diameter/position (p. 12, randomness).

## 5. Parameters & constants of interest
The "parameters" are dimensionless CC3D/GGH model knobs (adhesion energies in arbitrary units, MCS time), not transferable physical constants for a fine-grained particle simulator. Listed for completeness only:

| Quantity | Value+units | Source-in-paper |
|---|---|---|
| Cadherin relative expression, cancer cells | 8 (a.u.) | §2.2, from Swat et al |
| Cadherin relative expression, normal ("muted") cells | 16 (a.u.); scan 9/12/15 | §2.2 |
| Fibronectin expression density | 16 (a.u.) | §2.2 |
| Cadherin–cadherin binding constant | 2 (a.u., energy) | §2.2 |
| Cadherin-vs-integrin binding strength ratio | changed 0.2 → 0.02 | §2.2 |
| Cadherin mutation chance per division | 10% | §2.2 |
| Damage points before necrosis (muted cells) | raised 102 → 1020 | §2.3 |
| Time calibration | ≈ 20 s / MCS | §4.1, from Swat et al [14] |
| Glucose scan range | 0.1× – 10× physiologic | §2.4 |
| Spheroid muting radius | 15 pixels | §2.5 |
| Doublet center spacing | 32 pixels | §2.5 |

These are not physically grounded for ffn_cellsim (no SI units; adhesion is a Potts contact energy, not a Bell-Evans/catch-bond off-rate). **Not usable as a mechanistic oracle or parameter source.**

## 6. Relevance to ffn_cellsim  <-- MOST IMPORTANT
**Classification: (d) spheroid/multicellular-scale context, leaning (f) tangential.** This is a coarse, lumped-cell tissue-scale model — the exact opposite of ffn_cellsim's fine-grained mechanistic philosophy. Key mismatches:
- **Wrong mechanism class.** Cells are single Potts agents whose behavior emerges from Metropolis energy minimization over adhesion contact energies. ffn_cellsim's hard rule explicitly rejects Metropolis as a detailed-balance proxy in favor of Bell-Evans force-dependent off-rates, and rejects lumped cell agents in favor of explicit filaments/motors/clutches. There is **no cytoskeleton, no motor, no clutch, no ECM cross-link** at the particle level here.
- **Wrong scale.** This is multicellular (spheroid fusion, 76–81 cells/spheroid, tissue "boxes"), 2D, with cells as flat pixel domains. ffn_cellsim is single-cell, 3D, sub-filament resolution. The ×40 mesoscale coarse-graining ffn_cellsim sanctions is still far finer than a Potts cell.
- **No transferable physics.** Adhesion is in arbitrary Potts energy units; time is MCS. Nothing here calibrates a traction stress, modulus, off-rate, motor force-velocity, or branching angle.

**Where it is marginally useful**: (i) as *spheroid-scale context* for the eventual multi-cell cohesion overlay (PI's MCF7-spheroid platform) — it documents how cadherin-differential cell sorting and glucose metabolism shape fusion and tumor-cell escape, which is conceptual background for the multi-cell cohesion (the "c" term) layer, not the single-cell traction ("b") layer. (ii) The **DCM metric** (Eq 1) and Feret length/width fusion tracking are simple, reusable *measurement protocols* if ffn_cellsim ever simulates aggregate fusion. (iii) It is a clean example of the "paper-model-as-runtime" approach ffn_cellsim deliberately inverts — useful as a contrast reference. It is **off-topic as a mechanistic single-cell oracle or parameter source** and should be tagged tangential. Cited cadherin = cell–cell adhesion strength is conceptually adjacent to ffn_cellsim's catch-bond cadherins (KU-4.2) but provides no numbers usable for that mechanism.

## 7. Limitations & caveats
- Authors' own limitations (§4.4): ultra-simplified cellular structure/behavior intrinsic to the GGH class; **2D** (despite "CC3D" name — cells are pixel disks, not voxel spheres); model not spatially/temporally re-calibrated (inherited from Swat et al); cell–cell interaction strengths and their large-range variations are "arbitrary."
- Heuristic, hypothesis-generating only — no new experimental validation in this paper (comparisons are qualitative against prior literature, e.g. Susienka 2016).
- Time calibration is method-dependent and approximate (MCS↔seconds mapping is a stated assumption).
- For a fine-grained particle simulator: cells have no internal mechanics, so traction, cortical tension, retrograde flow, motor stall, etc. are entirely absent — the scale gap is fundamental, not bridgeable by parameter import.

## 8. Key figures / tables
- **Figure 1** — Fusion of normal/heterogeneous/tumor spheroid pairs at normal glucose: stages, glucose fields, subpopulation kinetics, and DCM trajectories. The core qualitative result figure.
- **Figure 4** — Doublet length and width vs time (≈24 h interval) for normal (green) vs tumor (red); means ±SD, n=3, p<0.05. Quantitative fusion-geometry readout.
- **Figures 5–7** — 3×3 microenvironment "box": tumor disappearance vs stabilization vs escape as a function of glucose and cadherin (normal 9/12/15 vs tumor 8); tumor-cell dispersion fraction outside a reference circle.
- **Eq 1** — DCM definition (Euclidean distance between two spheroids' centers of mass) — the paper's introduced metric.

## 9. Notable quotes / citable claims
- "the spheroids are formed and maintained through the general (bio)physical principle of the system's energy minimization, by maximal engagement of intercellular adhesion molecules (mainly cadherins)" (Introduction, p. 2).
- "Although the program is called CC3D, for computational reasons in practice it considers as 'cells' flat domains composed of pixels (rather than voxels), which by energy minimization generate 2D disks (rather than true 3D 'spheres')" (§2.5, p. 3).
- "Among the limitations of this model, as of those from all the GGH class, is the ultra-simplified cellular structure and behavior, yet sufficient for capturing essential features of cellular dynamics" (§4.4, p. 12).
- "the post-fusion tumor destabilization … increases the surface-to-volume … ratio, and thus the chance for single cells to leave the structure, migrate, and form metastases" (§4.3, p. 12).
