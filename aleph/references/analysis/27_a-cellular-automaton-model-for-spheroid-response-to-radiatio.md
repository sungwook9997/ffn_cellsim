---
id: 27_a-cellular-automaton-model-for-spheroid-response-to-radiatio
paper_n: 27
title: "A cellular automaton model for spheroid response to radiation and hyperthermia treatments"
authors: "Brüningk, S. C. et al."
year: "2019"
venue: "Scientific Reports"
doi: "10.1038/s41598-019-54117-x"
paper_type: agent-based
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, tangential]
entities: [hct116, spheroid, necrotic-core, oxygen, extracellular-matrix]
methods: [agent-based, cellular-automaton, reaction-diffusion, clonogenic-assay, histology, deep-learning-segmentation]
measurables: [spheroid-radius, oxygen-tension, proliferation-rate, clonogenic-survival-fraction, doubling-time, necrotic-core-radius]
keywords: [cellular-automaton, tumour-spheroid, radiotherapy, hyperthermia, oxygen-diffusion, hypoxia, clonogenic-survival, linear-quadratic-model, mitotic-catastrophe, systems-oncology]
tags: ["#spheroid-scale-context", "#tumour-spheroid", "#cellular-automaton", "#radiation-hyperthermia", "#oxygen-diffusion", "#tangential"]
has_transferable_params: false
---

# [27] A cellular automaton model for spheroid response to radiation and hyperthermia treatments

**Tags:** #spheroid-scale-context #tumour-spheroid #cellular-automaton #radiation-hyperthermia #oxygen-diffusion #tangential

| Field | Value |
|---|---|
| Authors | Brüningk, S. C., Ziegenhein P., Rivens I., Oelfke U., ter Haar G. |
| Year / Venue | 2019 / Scientific Reports 9:17674 |
| DOI / ID | 10.1038/s41598-019-54117-x |
| Type | agent-based (cellular automaton, voxel-per-cell) |
| Pages | 12 |
| ffn_cellsim relevance | Low — cell-as-voxel oncology growth model; no cytoskeleton/clutch/ECM mechanics at all |

## 1. Summary
The authors build a 3D cellular automaton (CA) "systems oncology" model in which each tumour cell occupies a single voxel on a 200×200×200 lattice (12 µm voxels) and progresses through a G1/S/G2/M/G0 cell cycle by an individualized timer. They couple this to an iteratively solved oxygen reaction-diffusion field to produce a hypoxic/necrotic core and a growth plateau, then add radiation response (linear-quadratic clonogenic survival, proliferation-dependent mitotic-catastrophe death) and hyperthermia response (thermal-dose CEM43 survival, proliferation-independent fast death with reoxygenation). The model is calibrated and validated against HCT116 colorectal-carcinoma spheroid growth curves for 0–10 Gy radiation, 0–240 CEM43 heat, and combinations, reproducing experimental spheroid-diameter dynamics with R² > 0.85 across all cases. The central biological finding: radiation and hyperthermia give comparable clonogenic survival yet very different spheroid growth dynamics because radiation kills only at division (proliferation-dependent) while heat kills quickly and independently of proliferation, allowing reoxygenation and repopulation.

## 2. Problem & motivation
Combined radiotherapy + hyperthermia is promising for radio-resistant (e.g. hypoxic) tumours, but standard biologically-equivalent-dose (BEQD) weighting based on 2D clonogenic survival cannot explain why isoeffective radiation vs heat doses produce different 3D spheroid growth dynamics. The authors argue the difference lies in the cell-death *mechanism* (proliferation-dependent mitotic catastrophe for RT vs fast proliferation-independent death for HT) and that a 3D agent-based model with explicit microenvironment (oxygen) is needed to capture treatment-induced growth-delay/shrinkage dynamics that lumped BEQD misses.

## 3. Methods / model
- **Model class**: discrete 3D cellular automaton, one cell = one cubic voxel (12 µm edge), grid 200³ voxels. Cells progress through G1/S/G2/M with an individualized timer; division searches Moore/von-Neumann neighbourhoods (up to 2nd order), favouring central sites for compactness; no free space → reversible quiescence (G0). No cell migration.
- **Microenvironment**: oxygen partial pressure pO₂ solved iteratively from a reaction-diffusion PDE, ∂pO₂/∂t = D_O₂·Δ(pO₂) − Φ(x), with constant per-voxel consumption Φ where occupied; Dirichlet boundary pO₂ = 100 mmHg outside the spheroid; equilibrated each cell-cycle step (≤0.01% max change). Cells below 11 mmHg labelled hypoxic; die with probability p_hypoxiaDeath; necrotic cells cleared with p_clearNecrotic; cells shuffled inward to maintain a dense sphere.
- **Radiation**: clonogenic surviving fraction S_RT from the linear-quadratic model S = exp[−γ·(α·d_OER + β·d_OER²)], with cell-cycle weighting γ and an oxygen-enhancement-ratio (OER 1–3) dose modifier d_OER. "Dying" cells die by mitotic catastrophe with p_mitoticDeath only when they next attempt division (two-stage: <0.5 then >0.5 after delay t_delayRT).
- **Hyperthermia**: thermal-dose (CEM43) survival S_HT from the authors' prior "AlphaR" model (step-wise function of t43, plateau S_HT,plateau). Dying cells die after a normally-distributed delay (mean t_delayHTtodeath) then are immediately cleared, enabling reoxygenation/repopulation; G0 cells 3× more heat-resistant than cycling cells.
- **Combination**: combined survival S_RTHT·S_HT with thermal-dose-dependent α,β for synergy.
- **Validation/numerics**: spheroid diameter = mean of 100 largest cell-to-centre distances (MATLAB 2017a); fits by grid search maximizing R² with a "one-parameter-at-a-time" strategy. Compared to the analytical Grimes et al. oxygenation model (Eqs 6–7).
- **Experimental side**: HCT116 spheroids (seed 300 cells, ULA U-bottom plates, 4-day pre-growth); 220 kVp X-rays 0–10 Gy; 47 °C heating; propidium-iodide live/dead; Celigo/Incucyte imaging; deep-learning segmentation; pimonidazole/Ki-67/DAPI histology.

## 4. Key results (quantitative)
- Untreated growth best fit R² = 0.99 with p_hypoxiaDeath = 0.01, p_clearNecrotic = 0.001, cell doubling time 28 h (p.7).
- Fitted oxygen diffusion coefficient D_O₂ = 3.8×10⁻⁹ m²/s; consumption rate 22.1 mmHg/s (from Grimes et al.) (p.6).
- Simulated vs analytical oxygen distributions agree to <15% (max deviation in outer cell layers); spheroid sizes modelled 284, 500, 664, 792 µm at days 4/11/18/25 (Fig 3, p.6).
- Radiation: 10 Gy calibration R² = 0.98 with p_mitoticDeath<0.5 = 0.44 (first 3.5 days = t_delayRT) then 0.59; validation R²(5 Gy) = 0.95, R²(2 Gy) = 0.98 (Fig 6 top, p.7).
- Hyperthermia: fits S_HT,plateau = 0.0005, mean death delay t_delayHTtodeath = 96 h, mean cell-cycle arrest 30 h; validation R²(240 CEM43) = 0.95, R²(80) = 0.94, R²(40) = 0.98; S(80 CEM43) ≈ 0.04 (Fig 6 middle, p.7–8).
- Combination (2 Gy + 40/80/160 CEM43): R² = 0.96 / 0.88 / 0.86 (Fig 6 bottom, p.8).
- Oxygen-related geometry (Grimes HCT116): diffusion limit r_l = 233 µm, necrotic-core radius r_n = 155 µm; hypoxia threshold 11 mmHg, medium pO₂ ≈ 100 mmHg (p.5–6).

## 5. Parameters & constants of interest
None transferable to a fine-grained single-cell cytoskeleton+ECM simulator. The constants here are tumour-population / oncology parameters at the wrong scale:

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| HCT116 doubling time | 28 h | p.7 (fit) |
| O₂ diffusion coeff D_O₂ | 3.8×10⁻⁹ m²/s | p.6 (fit) |
| O₂ consumption rate | 22.1 mmHg/s | p.6 (Grimes et al.) |
| Hypoxia threshold | 11 mmHg | p.2–3 |
| Medium pO₂ (boundary) | ~100 mmHg | p.2, Eq 6 |
| HCT116 voxel/cell size | 12 µm | p.2 |
| O₂ diffusion limit r_l / necrotic r_n | 233 µm / 155 µm | p.6 (Grimes 2016) |

These describe avascular-spheroid oxygen transport and clonogenic survival, not particle/bond mechanics — they are not oracles for ffn_cellsim.

## 6. Relevance to ffn_cellsim  <-- MOST IMPORTANT
**Category: (d) spheroid/multicellular-scale context, but mostly (f) tangential.** This is a population-scale oncology growth model where the entire cell is a single voxel with no internal mechanics — the exact opposite end of the modelling spectrum from ffn_cellsim, which resolves every filament, motor head, clutch, and cross-link as explicit particles/bonds. There is **no cytoskeleton, no cortex, no focal adhesion/clutch, no myosin, no ECM mechanics, and no membrane/nucleus mechanics** in this model. Mapped to the H-unit chain (H.1 cortex → H.2 → H.3 → H.5 → H.7) and the EXTEND track (membrane/nucleus/cytoplasm), it informs **none** of them mechanistically.

Its only thin connection: it is a **spheroid-scale context** reference. ffn_cellsim's validation overlay target is MCF7-spheroid-on-pV4D4/collagen-I (traction at single-cell, cohesion at multi-cell). This paper is a different cell line (HCT116 colorectal) and a different physical question (radiation/heat treatment response, oxygen-limited growth), so it is not a validation oracle and not a parameter source for the platform. It is useful only as background on (i) how the field models the multicellular regime where many single-cell models must eventually tile up, and (ii) what NOT to do for a mechanistic simulator — the authors explicitly note the CA's inability to model "loosening of the cellular structure", "mechanical compression", varying cell density, or ECM breakdown precisely because cells are dimensionless voxels (Discussion, p.9). That limitation is exactly the regime ffn_cellsim exists to resolve, so the paper is a clean illustration of the scale gap the project's fine-grained approach is designed to close, but it contributes no mechanism, parameter, numerics, or oracle. Recommend tagging Low / tangential with a spheroid-scale-context note; do not cite as a physics anchor.

## 7. Limitations & caveats
- Cell = one fixed voxel ⇒ cannot represent cell deformation, density change, mechanical compression, or ECM breakdown (authors flag this; suggest a future lattice-gas model). This is the central scale gap vs a particle simulator.
- "Oxygenation" is a lumped proxy that the authors admit may stand in for glucose, waste, pH and other diffusive factors — explicitly a deliberate parameter-minimizing simplification.
- All clonogenic-survival mechanisms (LQ model, CEM43/AlphaR, OER) are phenomenological closed forms fitted to data, not derived from cell mechanics — in ffn_cellsim terms these are population-level oracles at best, and even then for a treatment-response question outside the platform's scope.
- Parameters are fitted one-at-a-time to a specific artificial system; the authors state absolute values "should not be used for direct interpretation of biological meaning."
- No cell migration; spherical compactness imposed by inward shuffling rather than emerging from mechanics.
- Companion experimental dataset "to be published separately"; framework not yet public at time of writing.

## 8. Key figures / tables
- **Fig 1 (p.3–4)**: schematic + probability-driven response cascades for RT (proliferation-dependent, death at division) vs HT (random delay, immediate clearance) — the conceptual core distinguishing the two death mechanisms.
- **Fig 3 (p.6)**: iteratively simulated vs analytical (Grimes) oxygen partial-pressure maps for 4 spheroid sizes; difference map shows <15% agreement.
- **Fig 5 (p.7)**: untreated growth curves (sim vs experiment) + live/dead cross-sections vs PI time-lapse microscopy — necrotic core agreement.
- **Fig 6 (p.7–8)**: the validation centrepiece — sim vs experimental spheroid-diameter growth curves for RT, HT, and RTHT combinations (all R² > 0.85).

## 9. Notable quotes / citable claims
- "Despite comparable clonogenic survival, spheroid growth differed significantly following radiation or hyperthermia." (Abstract) — the headline mechanism-vs-survival decoupling.
- "Heat-induced cell death was implemented as a fast, proliferation-independent process, allowing reoxygenation and repopulation, whereas radiation was modelled as proliferation-dependent mitotic catastrophe." (Abstract).
- "Due to the fixed lattice size, it was not easy to implement treatment-induced variations in cell size, a loosening of the cellular structure (i.e. cell density variation), or mechanical compression of cells within the necrotic spheroid core." (Discussion, p.9) — the explicit mechanics scale-gap.
- "absolute parameter values refer entirely to an artificial, simplified biological system and thus should not be used for direct interpretation of biological meaning." (Methods, p.5).
