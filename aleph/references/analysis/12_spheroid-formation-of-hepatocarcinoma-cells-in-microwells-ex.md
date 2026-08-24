---
id: 12_spheroid-formation-of-hepatocarcinoma-cells-in-microwells-ex
paper_n: 12
title: "Spheroid Formation of Hepatocarcinoma Cells in Microwells: Experiments and Monte Carlo Simulations"
authors: "Wang Y et al."
year: "2016"
venue: "PLoS ONE"
doi: "10.1371/journal.pone.0161915"
paper_type: agent-based
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, cortex, tangential]
entities: [huh-7.5, hepatocarcinoma, f-actin, ki-67, spheroid, silicone-elastomer]
methods: [monte-carlo, lattice-gas, confocal-live-imaging, immunofluorescence, microwell-culture]
measurables: [spheroid-radius, proliferation-rate, occupation-fraction, cell-cell-adhesion-energy]
keywords: [spheroid-formation, monte-carlo, lattice-model, cell-cell-adhesion, cell-surface-adhesion, gravity, proliferation, diffusion, huh-7.5, microwell, cortical-actin, ki-67]
tags: ["#spheroid-mechanics", "#monte-carlo", "#lattice-model", "#hepatocarcinoma", "#multicellular-scale", "#cell-aggregation", "#tangential"]
has_transferable_params: false
---

# [12] Spheroid Formation of Hepatocarcinoma Cells in Microwells: Experiments and Monte Carlo Simulations

**Tags:** #spheroid-mechanics #monte-carlo #lattice-model #hepatocarcinoma #multicellular-scale #cell-aggregation #tangential

| Field | Value |
|---|---|
| Authors | Wang Y, Kim MH, Tabaei SR, Park JH, Na K, Chung S, Zhdanov VP, Cho NJ (et al.) |
| Year / Venue | 2016 / PLoS ONE |
| DOI / ID | 10.1371/journal.pone.0161915 |
| Type | experimental + agent-based (lattice Monte Carlo) |
| Pages | 13 |
| ffn_cellsim relevance | Low — multicellular spheroid pattern-formation at the population (lattice) scale; no single-cell cytoskeletal mechanism, no transferable physical constants |

## 1. Summary
The authors culture Huh-7.5 human hepatocarcinoma cells in a microfabricated low-adhesion silicone-elastomer microwell platform (SpheroFilm; 500 μm well diameter, 1000 μm depth) and observe spheroid formation over 10 days at three seeding densities (0.05/0.1/0.2 × 10⁶ cells/ml), tracked by light microscopy plus confocal immunofluorescence of F-actin (phalloidin), the proliferation marker Ki-67, and DAPI nuclei. Experimentally they find cells aggregate toward the well center/bottom, proliferation (Ki-67) is fastest at day 1 then suppressed by day 4 and confined to the spheroid periphery, and F-actin shifts from stress-fiber (monolayer) to cortical (spheroid) distribution. They complement this with a 3D lattice (cubic, 60×60×80, hemisphere-bottom + cylinder-top) Monte Carlo model in which each site is empty or holds one cell, cells diffuse (nn jumps) and divide, with dimensionless cell–cell (cc = −1.2) and cell–surface (cs) attractions, a division/diffusion balance parameter p_div, and an optional gravity bias (g = −0.03). Varying these reproduces the qualitative regimes: gravity drives accumulation at the well bottom matching experiment, strong cell–cell + weak cell–surface adhesion favors compact central spheroids, and growth kinetics are exponential-then-linear.

## 2. Problem & motivation
3D cell culture in scaffolds matters for tissue engineering, disease modeling, and drug screening, and cellular spheroids are the most common aggregate shape. The paper asks: in a specific low-adhesion microwell, how do cell diffusion, division, cell–cell vs cell–surface adhesion, and gravity together shape the formation, location, and size of Huh-7.5 spheroids? It pairs direct imaging with a coarse-grained statistical-physics model to give an integrated mechanistic-at-the-population-level picture.

## 3. Methods / model
- **Experimental:** Huh-7.5 hepatocarcinoma line in DMEM + 10% FBS; seeded at 0.05/0.1/0.2 × 10⁶ cells/ml into SpheroFilm silicone-elastomer microwells (oxygen-permeable, shear-free triangular-wall design; 500 μm inner diameter hemispheres, 1000 μm total depth, 361 wells in 19×19 array). Readouts: light microscopy over 10 days; confocal immunofluorescence (LSM 710) of Ki-67 (Alexa 488), F-actin (Alexa 555 phalloidin), nuclei (DAPI), plus DIC.
- **Model class:** 3D lattice-gas / kinetic Monte Carlo (agent-based at the cell level — one cell per lattice site). NOT a cytoskeletal or particle-mechanics model; cells are dimensionless lattice occupants.
- **Mechanisms:** cell diffusion via nn jumps, cell division into vacant nn sites (no-flux boundaries), dimensionless (kBT-normalized) cell–cell attraction cc < 0 and cell–surface attraction cs ≤ 0. "Initial-state" dynamics: jump rate reduced by exp(n·cc) for n neighbors and exp(m·cs) for m wall contacts. Division/diffusion balance via p_div (rates ∝ p_div and 1−p_div). Gravity modeled by reducing upward-jump probability by exp(g), g < 0.
- **Length/time scales:** lattice spacing ≈ one cell; well represented by 60×60×80 sites (hemisphere at z<30, cylinder at 30≤z≤80). Time in Monte Carlo steps (MCS); one MCS ≈ N trials; runs to occupation ϑ = 0.4; started from 5 random seed cells.
- **Numerics:** Gillespie-like time increment Δt = |ln ρ|/N per trial; random-site selection; aggregation threshold known to occur at |cc| > 0.89 (ref. lattice-gas literature).

## 4. Key results (quantitative)
- Microwell geometry: 500 μm inner diameter, 1000 μm total depth (depth = 2× diameter), 361 (19×19) wells; observed spheroids up to ~500 μm scale (p.2; Fig 1).
- Cell–cell aggregation threshold |cc| > 0.89; chosen cc = −1.2 (dimensionless, /kBT) (p.6).
- Gravity bias g = −0.03 reproduces bottom-of-well accumulation matching experiment (Fig 9; p.7).
- p_div regimes (cs = 0): p_div = 1 (no diffusion) → few aggregates around seeds (Fig 5); p_div = 10⁻² (slow diffusion) → many small aggregates (Fig 6); p_div = 10⁻⁴ (rapid diffusion) → fewer, larger aggregates (Fig 7) (p.7).
- Cell–surface adhesion cs = −1.2 with p_div = 10⁻⁴ → cells line the microwell walls (Fig 8; p.7).
- Growth kinetics: exponential initially, then approximately linear; diffusion accelerates kinetics (Fig 10; p.9).
- Proliferation timeline (experiment): Ki-67 highest at day 1, markedly reduced by day 4, roughly constant days 4–10; after day 4 strong Ki-67 only at spheroid periphery (white arrows, Fig 3); cells outside wells (flat surface) keep high Ki-67 (Fig 4) (p.4–5).
- Nutrient/O₂ limitation expected only above ~200–400 μm aggregate scale (ref. 36); observed aggregates ≤ this, so growth not nutrient-terminated (p.6).
- Spheroid size literature anchors cited: up to ~1 cm (colon carcinoma, rat hepatocytes), ~100 μm (mammary epithelial), ~200 μm (hepatocytes) (p.1–2).

## 5. Parameters & constants of interest
None (no transferable physical constants for a fine-grained single-cell simulator). The model parameters are dimensionless lattice-gas quantities (cc, cs in units of kBT; p_div; g) with no mapping to forces, energies in SI, motor rates, bond lifetimes, or moduli. The only dimensional values are the device geometry and the qualitative spheroid size scales.

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Microwell inner diameter | 500 μm | Methods, p.2 / Fig 1 |
| Microwell total depth | 1000 μm | Methods, p.2 / Fig 1 |
| Cell–cell interaction (dimensionless) | cc = −1.2 (/kBT) | p.6 |
| Aggregation threshold | |cc| > 0.89 (/kBT) | p.6 (ref. 37) |
| Cell–surface interaction (dimensionless) | cs = −1.2 (or 0) (/kBT) | p.7 |
| Gravity bias energy | g = −0.03 (/kBT) | p.7 |
| Division/diffusion balance | p_div = 1, 10⁻², 10⁻⁴ | p.7 |
| Run-termination occupancy | ϑ = 0.4 | p.6 |
| Nutrient/O₂ limit length scale | ~200–400 μm | p.6 (ref. 36) |

## 6. Relevance to ffn_cellsim  — MOST IMPORTANT
This is primarily **(d) spheroid/multicellular-scale context** and otherwise **(f) tangential** to a fine-grained mechanistic single-cell simulator. It is **not** a validation oracle, parameter source, mechanism reference, or numerics reference for ffn_cellsim:

- **No single-cell cytoskeletal mechanism.** Cells are structureless lattice occupants; there is no actin filament, motor, clutch, adhesion bond, or membrane in the model. This directly conflicts with ffn_cellsim's architectural principle (every filament/motor/clutch/cross-link is an explicit particle/bond). The paper's "adhesion" is a dimensionless kBT-normalized Ising-like attraction, the antithesis of the mechanistic Bell-Evans / catch-bond clutches used here.
- **No transferable parameters.** cc, cs, g, p_div are abstract lattice energies/rates with no SI mapping; they cannot seed Hill motor, Bell-Evans off-rate, LJ EV, or cortex parameters.
- **Numerics are off-model.** Kinetic Monte Carlo on a cubic lattice with initial-state jump dynamics is a coarse-grained population method, unrelated to the Leimkuhler-Matthews BAOAB Langevin integrator used in ffn_cellsim.
- **Where it touches the project, weakly:**
  - *Spheroid-scale context for the PI validation overlay.* ffn_cellsim's experimental overlay target is MCF7 spheroids on pV4D4/collagen-I, with cohesion read at the multi-cell scale. This paper is an example of the spheroid-formation phenomenology (peripheral-only proliferation, central core suppression, role of cell–cell vs cell–surface adhesion balance) the platform's multicellular EXTEND layer would eventually need to reproduce — but as background phenomenology, not a quantitative anchor. Different cell line (Huh-7.5 hepatocarcinoma, not MCF7), different substrate.
  - *Cortex (H.1) qualitative note.* The paper cites/observes the monolayer→spheroid F-actin transition (stress fibers → cortical actin distributed at the cell outline, ref. Chang & Hughes-Fulford). This is consistent with H.1's cortical-actin focus but is a qualitative immunofluorescence observation, not a measurement ffn_cellsim can calibrate against.

Net: keep as multicellular spheroid-formation background. Do not mine it for parameters, mechanisms, or oracles. It sits well above the single-cell fine-grained scale ffn_cellsim operates at.

## 7. Limitations & caveats
- **Scale gap.** The model is a population-level lattice gas (one cell = one site), orders of magnitude coarser than ffn_cellsim's per-filament particles. No mechanics, no forces, no geometry within a cell.
- **Self-declared omissions.** Cell–cell communication ignored (no data); nutrient/oxygen limitation ignored (justified only for <200–400 μm aggregates); aggregate sedimentation/diffusion ignored.
- **Dimensionless, non-calibrated parameters.** cc, cs, g, p_div are swept qualitatively "over a wide range" to illustrate regimes, not fit to data — so even the model-internal numbers are illustrative, not measured.
- **Gravity term is ad hoc.** An exp(g) upward-jump penalty is a phenomenological stand-in, not a force balance.
- **Experimental readouts are imaging/IF only** (Ki-67, F-actin, DAPI, light microscopy); no mechanical measurements (no traction, modulus, or adhesion force).

## 8. Key figures / tables
- **Fig 1** — SpheroFilm scaffold scheme: 19×19 = 361 microwells, 500 μm diameter / 1000 μm depth geometry; seeding and spheroid-formation cartoon. (Geometry source.)
- **Fig 3** — Ki-67 / F-actin / DAPI confocal of Huh-7.5 spheroids over days 1/4/7/10: shows proliferation moving to the periphery and the F-actin (cortical) signal evolving.
- **Figs 5–9** — MC cross-sections through the well at occupancy W = 0.2/0.4 across regimes (no diffusion; slow/rapid diffusion; cell–surface adhesion; gravity), showing where aggregates localize.
- **Fig 10** — MC growth kinetics (curves 1–5 for Figs 5–9): exponential→linear population growth, diffusion speeds it up.

## 9. Notable quotes / citable claims
- "In this case, the aggregation is known to occur at |cc| > 0.89. In our calculations, we set cc = −1.2." (p.6) — the model's cell–cell adhesion threshold and chosen value.
- "the cell spheroids would form on the basis of strong cell-cell interaction and weak cell-surface adhesion, and growth of spheroids mainly depends on the cell proliferation on their periphery." (Conclusion, p.10) — the headline mechanism claim.
- "F-actin stress fibers were found in cells of monolayer culture, while cells in spheroids featured cortical actin that clearly distributed at the outline of the cells." (p.5, citing Chang & Hughes-Fulford) — the monolayer→spheroid cortical-actin transition.
- "Such limitations [nutrient/oxygen] are expected to be significant roughly on the aggregate length scale larger than 200–400 μm." (p.6) — the scale below which growth is not diffusion-limited.
