---
id: 25_simulation-based-inference-of-cell-migration-dynamics-in-com
paper_n: 25
title: "Simulation-based inference of cell migration dynamics in complex spatial environments"
authors: "Jonas Arruda et al."
year: "2026"
venue: "npj Systems Biology and Applications"
doi: "10.1038/s41540-026-00648-9"
paper_type: agent-based
ffn_relevance: Low
ffn_themes: [numerics/methods, spheroid-scale-context, tangential]
entities: [dendritic-cell, ccl19, ccr7, nucleus]
methods: [cellular-potts, agent-based, monte-carlo, bayesian-inference, neural-posterior-estimation, approximate-bayesian-computation, confocal-live-imaging, single-particle-tracking]
measurables: [migration-speed, migration-persistence, turning-angle, normalized-displacement, cell-area, waiting-time-rate]
keywords: [cellular-potts-model, simulation-based-inference, neural-posterior-estimation, approximate-bayesian-computation, dendritic-cell-migration, chemotaxis, pillar-forest, summary-statistics, normalizing-flow, parameter-inference, persistent-random-walk]
tags: ["#cellular-potts", "#simulation-based-inference", "#neural-posterior-estimation", "#cell-migration", "#chemotaxis", "#parameter-inference", "#dendritic-cell", "#methods"]
has_transferable_params: false
---

# [25] Simulation-based inference of cell migration dynamics in complex spatial environments

**Tags:** #cellular-potts #simulation-based-inference #neural-posterior-estimation #cell-migration #chemotaxis #parameter-inference #dendritic-cell #methods

| Field | Value |
|---|---|
| Authors | Jonas Arruda et al. (Alamoudi, Mueller, Vaisband, Molkenbur, Merrin, Kiermaier, Hasenauer) |
| Year / Venue | 2026 / npj Systems Biology and Applications (2026)12:20 |
| DOI / ID | 10.1038/s41540-026-00648-9 |
| Type | agent-based (Cellular Potts model) + machine-learning simulation-based inference (methods) |
| Pages | 15 |
| ffn_cellsim relevance | Low — a lattice CPM + Bayesian-inference methods paper at the cell-shape/centroid scale; no cytoskeletal/clutch/motor mechanism, no transferable molecular constants |

## 1. Summary
The authors record dendritic-cell (DC) migration through a microfabricated PDMS "pillar forest" (10 µm gaps) under a CCL19 chemokine gradient, tracking nuclear positions every 30 s by wide-field microscopy. They build a 2D extended Cellular Potts model (CPM, Graner-Glazier base) with four free parameters — chemokine attraction strength (m_dir), intrinsic random-motion strength (m_rand), waiting-time rate (λ) of directional changes, and target cell area (a) — and calibrate it against the trajectory data using simulation-based inference (SBI). They benchmark four likelihood-free strategies: ABC-SMC with hand-crafted summary statistics, ABC with posterior-mean neural summaries (ABC-PM), ABC with inference-tailored neural summaries (ABC-NPE), and full neural posterior estimation (NPE) with jointly trained summary networks. NPE achieves the lowest normalized error (NRMSE 0.20) with ~10× fewer simulations (~32,000 vs 66,000–432,000 for ABC) and amortized ~1-s inference, while hand-crafted summaries (MSD-like, turning-angle) fail to capture richer spatiotemporal structure. A simulation study concludes that directional persistence — not just chemokine guidance — is the key driver of cell dispersion past obstacles, and that constant (vs transient) chemokine can paradoxically trap cells in the pillar forest.

## 2. Problem & motivation
How do physical geometry (obstacles, confinement) and chemical gradients jointly shape immune-cell migration, and how can mechanistic migration models be calibrated from experimental trajectories when no explicit likelihood exists? CPM and similar models have no tractable likelihood, so classical Bayesian methods (MCMC) are unavailable; ABC depends on hand-crafted summary statistics that may discard information. The paper's central methodological contribution is showing that neural, inference-tailored summary statistics outperform hand-crafted ones for calibrating spatially-structured cell-migration models.

## 3. Methods / model
- **Model class**: 2D extended Cellular Potts model (lattice-based, Metropolis copy-attempt dynamics), built in MorpheusML (model M6342) and driven via the FitMultiCell pipeline. Spatial resolution 1.31 µm/node; cells are connected lattice domains; pillar sites reject copy attempts (hard obstacles).
- **Energy function** (Eq. 1–4): classical adhesion term J·(1−δ) + volume constraint λ_v·(v_c − a)² , plus two extensions — chemotaxis ΔH_chemotaxis = −m_dir·(f_xi − f_xj) (Eq. 2, biases moves up the chemoattractant concentration) and persistence ΔH_persistence = −m_rand·v_c·⟨(cos α, sin α)ᵀ, s_c⟩ (Eq. 3, run-and-tumble-like persistent random walk; direction α resampled after waiting time t_α ~ Exp(λ)).
- **Chemokine field**: time-independent Gaussian, peak 7×10⁶ µm⁻² at the chemokine hole, σ = 550 µm.
- **Inference**: likelihood-free. Hand-crafted summaries = normalized displacement, velocity, turning angle, angle-degree distributions, compared with Wasserstein-1 distance in pyABC SMC (1000 particles, stop at acceptance rate 0.01 or 15 generations). NPE uses BayesFlow neural spline normalizing flows (6–8 layers) with a summary net = 1D-conv + GRU(32 units) + temporal attention pooling → fixed 8-dim population summary; trained by minimizing KL divergence; validated by simulation-based calibration (SBC).
- **Experiment**: bone-marrow-derived DCs (C57BL/6J mice), matured with LPS, exposed to a stable CCL19 gradient (7 ng, 1 g/mL in hole; CCR7 receptor) in a PDMS pillar chip (~5 µm height, 10 µm gaps); nuclei stained NucBlue, tracked at 30-s intervals over 2 h with TrackMate/Fiji. 143–147 cell trajectories.
- **Compute**: one CPM simulation 10–30 s; runs on a CPU cluster (6×32-core Intel Xeon Sapphire Rapids) + one Nvidia A40 GPU for NN training.

## 4. Key results (quantitative)
- NPE: NRMSE 0.20 using only 32,100 simulations, ~4 h training, ~1 s inference (Table 1).
- ABC (hand-crafted): NRMSE 0.23, 66,000–413,000 simulations, 2.1–19.2 h inference, 12–16 generations (Table 1).
- ABC-NPE: NRMSE 0.23, 120,000–428,000 simulations (Table 1). ABC-PM: NRMSE 0.48 as feature extractor (systematic bias in m_rand) though 0.18 when directly predicting parameters (Table 1, text p.4).
- ABC requires ~10× more simulations than NPE for comparable accuracy (Discussion, p.8).
- Experimental data: 143 trajectories, 18–120 frames/cell, mean duration 33.2 min (p.6).
- Inferred median target cell area a = 41.73 µm² (≈ 7.29 µm effective circular diameter), below the 10 µm pillar gap and the 10–15 µm mature-DC size — interpreted as confinement response (p.7, Fig 5C; ref 48 / BNID 113239).
- Inferred median waiting-time rate λ = 2.1 ms⁻¹ → mean waiting time 7.9 min (>>30-s imaging interval), limiting λ identifiability (p.7).
- Simulation study: persistence drives dispersion even without a gradient; constant chemokine can trap cells in the pillar forest; transient/pulsed chemokine (active first 50–75% of time) yields higher arrival fraction at the source (Fig 6).
- Inference time for full ABC on experimental data: 200,000–300,000 simulations, 27–50 h, converged after 13 generations (p.5–6).

## 5. Parameters & constants of interest

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Pillar gap (confinement spacing) | 10 µm | p.2 / Methods, Fig 1A |
| Imaging interval | 30 s | Methods, p.9 |
| CPM lattice resolution | 1.31 µm / node | Model description, p.10 |
| Chemokine field (Gaussian) | peak 7×10⁶ µm⁻², σ = 550 µm | Model description, p.10 |
| Inferred DC target area a | 41.73 µm² (≈ 7.29 µm diameter) | Fig 5C, p.7 |
| Mature DC reference size | 10–15 µm | ref 48 / BNID 113239, p.7 |
| Inferred waiting-time rate λ | 2.1 ms⁻¹ (≈ 7.9 min mean wait) | p.7 |
| Cell migration speed scaling | √T random-walk displacement | Eq. 6, p.11 |

These are phenomenological CPM/cell-shape parameters and assay geometry, not molecular/mechanistic constants (no stiffness, no bond rate, no motor force). Not usable as ffn_cellsim oracles. The DC size (10–15 µm) and pillar gap (10 µm) are generic confinement-scale references only.

## 6. Relevance to ffn_cellsim  — MOST IMPORTANT
This paper is **(f) largely tangential** to ffn_cellsim's mechanistic single-cell mandate, with a thin **(e) numerics/methods** thread.

- **Wrong mechanism class for the project's hard rule.** ffn_cellsim is explicitly fine-grained: every filament, motor head, clutch, and cross-link is an explicit particle/bond, and lumped/paper-model wrappers are forbidden at runtime. A Cellular Potts model is the canonical *lumped* mechanism — cells are lattice-domain energy functionals with phenomenological adhesion (J), volume (λ_v), chemotaxis (m_dir), and persistence (m_rand) terms. CPM is precisely the abstracted, mesoscale approach the project's architectural principle rejects. There is no cortex, no FA/clutch, no myosin, no ECM cross-link physics here.
- **Cell type and assay are off-target.** Dendritic-cell amoeboid chemotaxis through pillar forests is not the MCF7/MDA breast-cancer spheroid-on-pV4D4/collagen-I validation platform; no traction, no spreading, no cohesion readout overlaps.
- **No transferable molecular constants.** All inferred quantities are CPM-scaling parameters tied to one lattice resolution (the authors explicitly state the parameters are "optimized only for this resolution and cell size"). They cannot anchor a HOOMD particle/bond oracle. `has_transferable_params: false`.
- **The one genuinely useful angle is methodological (low priority).** The simulation-based inference machinery — likelihood-free calibration of an expensive stochastic simulator via neural posterior estimation (BayesFlow normalizing flows + learned summary networks), benchmarked against ABC-SMC — is directly relevant to a *future* problem ffn_cellsim will face: calibrating an expensive, likelihood-free HOOMD cell simulator against experimental trajectories/traction without hand-crafting summary statistics. If the project ever needs to infer mechanistic parameters (e.g. motor density, clutch off-rate) from PI's experimental overlays, the NPE-with-learned-summaries pattern (and the ~10×-fewer-simulations result, amortized inference, SBC calibration check) is a credible template. This maps to a hypothetical numerics/methods + parameter-inference track, not to any current H-unit.
- **Spheroid/multicellular-scale context (weak).** The "geometry modulates migration; persistence drives dispersion through obstacles" finding is generic confinement physics; it could loosely inform expectations for the spheroid-scale-context EXTEND track but contributes no quantitative anchor.

Bottom line: keep as a methods reference for the SBI/parameter-inference toolchain (BayesFlow, pyABC, FitMultiCell, SBC); do not treat as a mechanism, parameter, or validation-oracle source for the cytoskeleton+ECM core.

## 7. Limitations & caveats
- CPM is a coarse lattice energy model: it captures cell shape/area and emergent motility but **not** sub-cellular mechanics — no explicit cytoskeleton, no force generation, no traction, no clutch. The "target area" is a soft compactness constraint that cannot represent 3D/vertical deformation (authors note BMDCs deform out of plane, which the 2D model misses).
- Environment is homogeneous and static: no chemokine diffusion/decay over time, no cell-driven ECM remodeling or self-generated gradients (authors flag this as future work).
- Parameters are resolution-locked (1.31 µm/node) and not portable to other grids or to particle-based simulators.
- λ is poorly identifiable because the 30-s imaging interval is far below the ~7.9-min mean waiting time.
- Inference validated mostly on synthetic data; only 143 experimental trajectories, with high cell-to-cell variability and tracking limited to the pillar-forest window (cannot observe arrival at the source).

## 8. Key figures / tables
- **Table 1** (p.3): the headline benchmark — simulations, training/inference time, generations, and NRMSE for ABC / ABC-PM / ABC-NPE / NPE. NPE: 32,100 sims, ~1 s inference, NRMSE 0.20.
- **Fig 1** (p.2): experimental + modeling schematic — DC entry hole → CCL19 gradient → pillar forest; the four CPM parameters (m_dir, m_rand, a, λ).
- **Fig 2** (p.3): inference pipeline; summary-network architecture (conv + GRU + attention) and joint NPE training objective; LASSO latent-space view.
- **Fig 6** (p.7): simulation study — trajectories + turning-angle distributions across conditions (with/without chemokine, transient pulse, no persistence, no pillars); shows persistence and pulsed chemokine effects on arrival fraction.

## 9. Notable quotes / citable claims
- "classical summary statistics, such as mean squared displacement and turning angle distributions, can resolve key mechanistic features but fail to extract richer spatiotemporal patterns, limiting accurate parameter inference." (Abstract, p.1)
- "ABC requires 10 times more simulations to attain comparable accuracy, thereby highlighting the efficacy of NPE." (Discussion, p.8)
- "the inferred cell migration parameters are optimized only for this resolution and cell size" (Model description, p.10) — explicit non-portability of the CPM parameters.
- "it is the combination of chemokine signaling and persistence that can lead to cell entrapment ... temporally modulated guidance cues can facilitate more efficient exploration in the presence of obstacles." (Results, p.7–8)
