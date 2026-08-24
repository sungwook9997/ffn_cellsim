---
id: 23_development-of-a-scoring-function-for-comparing-simulated-an
paper_n: 23
title: "Development of a scoring function for comparing simulated and experimental tumor spheroids"
authors: "Herold J, Behle E, Rosenbauer J, Ferruzzi J, Schug A"
year: "2023"
venue: "PLOS Computational Biology"
doi: "10.1371/journal.pcbi.1010471"
paper_type: methods-software
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, ECM/collagen, validation-oracle, numerics/methods, tangential]
entities: [mda-mb-231, collagen-i, matrigel, dapi, ecm, spheroid]
methods: [cellular-potts, agent-based, monte-carlo, multiphoton-microscopy, voronoi-tessellation, marching-cubes, wasserstein-distance, bayesian-inference]
measurables: [spheroid-radius, spheroid-diameter, cell-density-distribution, voronoi-cell-volume, spheroid-surface-area, gaslike-fraction, deviation-score]
keywords: [tumor-spheroid, scoring-function, deviation-score, point-cloud, cellular-potts-model, cells-in-silico, nastjapy, collagen-density, invasion-phenotype, wasserstein-distance, mda-mb-231, sim-vs-experiment-comparison]
tags: ["#spheroid-mechanics", "#validation-oracle", "#breast-cancer", "#cellular-potts", "#agent-based", "#ecm-collagen", "#methods-software"]
has_transferable_params: true
---

# [23] Development of a scoring function for comparing simulated and experimental tumor spheroids

**Tags:** #spheroid-mechanics #validation-oracle #breast-cancer #cellular-potts #agent-based #ecm-collagen #methods-software

| Field | Value |
|---|---|
| Authors | Herold J, Behle E, Rosenbauer J, Ferruzzi J, Schug A |
| Year / Venue | 2023 / PLOS Computational Biology |
| DOI / ID | 10.1371/journal.pcbi.1010471 |
| Type | methods-software (with coarse-grained agent-based/Cellular-Potts simulation + experimental overlay) |
| Pages | 26 |
| ffn_cellsim relevance | Low — multicellular spheroid-scale comparison toolbox; lumped CPM cells, not a fine-grained cytoskeleton; useful only as far-field context + a possible sim-vs-experiment scoring idea |

## 1. Summary
The paper presents a **method/software** ("Nastjapy") for quantitatively comparing 3D tumor spheroids — regardless of whether they come from simulation or experiment — by reducing each spheroid to a **point cloud of cell centers** and extracting five spatial features (cell density distribution, gaslike-cell distribution, Voronoi cell-volume distribution, spheroid surface area, spheroid surface deformation). Each feature gets a distance metric (MSE, Euclidean, or 1-Wasserstein/earth-mover's distance), which are standardized and linearly combined (with optimized weights λf) into a single scalar **overall deviation score Di,j**. They generate four simulated phenotypes ("spherical", "deformed", "spherical with far gaslikes", "disordered") with their in-house **Cells in Silico (CiS)** framework (a Cellular Potts Model + agent layer), use them to calibrate the weights and validate the score (rotation/translation invariance + monotony under noise/deformation/scaling transformations), then apply the score to experimental MDA-MB-231 breast-cancer spheroids invading collagen-I at 1–4 mg/ml (data from Kang et al. 2021). The score reproduces the known experimental invasion transition between 2 and 3 mg/ml collagen at day 3, and can also score simulated-vs-experimental pairs.

## 2. Problem & motivation
Both experimentalists and modelers produce 3D spheroid data, but there is no systematic way to **quantitatively compare** two spheroids of arbitrary origin. Such a scalar distance is needed (a) for experimentalists to quantify how a single variable (e.g. collagen density) changes invasion, and (b) for theorists to use as an **objective function** when fitting simulation parameters to experimental data ("closing the loop between modeling and experiments"). The high dimensionality of even a single spheroid snapshot (e.g. R^3000 for 1000 cells) means any scalar metric is necessarily lossy and must be use-case-adapted — hence the weight-optimization scheme.

## 3. Methods / model
- **Comparison method (the paper's core)**: spheroid → point cloud of cell-center coordinates at one time. Five features extracted:
  1. *Cell density distribution* — fraction of cells in concentric spherical shells (central local density). Metric: 1-Wasserstein.
  2. *Gaslike cell distribution* — detached "gaslike" cells defined via Eq 1: distance from center > Dcrit AND min nearest-neighbor distance > dcrit (Dcrit = 125 μm, dcrit = 19 μm from Kang et al. data). Feature is a 2-tuple (fraction detached, normalized mean distance). Metric: Euclidean.
  3. *Voronoi cell-volume distribution* — Voronoi tessellation of cell centers; histogram of cell volumes (proxy for confinement/deformation). Metric: 1-Wasserstein.
  4. *Spheroid surface area* — voxelize non-gaslike points → marching-cubes triangulation → triangle-mesh area. Metric: MSE.
  5. *Surface deformation* — per-vertex scalar product of normal vs origin vector (=1 for perfect sphere); histogram. Metric: 1-Wasserstein.
- **Score**: standardize each metric distance (Eq 4–5), then Di,j = Σ λf · d*i,j,f (Eq 6). Weights λf found by maximizing inter-phenotype minus intra-phenotype distance (Eq 7) subject to Σλf² = 1 (Eq 8), solved with SciPy SLSQP. Framed as "inverse clustering" / contrastive-learning-like with a linear model.
- **Simulation model (Cells in Silico, CiS)**: NOT fine-grained. A **3D Cellular Potts Model** (Graner-Glazier 1992) at the microscale (Hamiltonian Eq 14 with cell-volume, cell-surface, and cell-cell-adhesion terms), nutrient/signal mesoscale (unused here), and an **agent-based macroscale** for division, motility-direction update, mutation. Motility added as a directional potential (Eq 15, persistent random walk via Wiener process, persistence p∈[0,1]). ECM = rigid, non-displaceable overlapping fibers that cells adhere to and degrade (lattice-point removal at 50% probability per degradation event). Time: 250,000 Monte-Carlo steps per run, 1 MC step ≈ 1 s real time. System (800 μm)³, ~2000 cells, initial spheroid diameter 200 μm, 5 replicates per phenotype.
- **Experimental data (from Kang et al., re-used)**: triple-negative MDA-MB-231 breast cancer cells, ~1000 cells coalesced into 300–400 μm spheroids (2.5% Matrigel, low-attachment 96-well, 48 h), embedded in rat-tail collagen-I at 1/2/3/4 mg/ml, cultured 1 h / 24 h / 48 h / 72 h (days 0–3), fixed, optically cleared, DAPI-stained, imaged on a Bruker Ultima multiphoton microscope (16× water objective, 0.8 NA, 3 mm WD). 3D nuclei positions extracted as point clouds (n=3/group with exceptions).

## 4. Key results (quantitative)
- Four simulated phenotypes separated by **3 parameters** (Table 1): cell-ECM adhesion 50 vs 450, ECM degradation period (disabled vs 5000 MC steps), motility magnitude 0 vs 100. (Other params fixed: ECM density 1, radially aligned, cell-cell adhesion 50, persistence 0, division enabled.)
- Fitted feature weights λf (Table 3): cell density 0.41, gaslike 0.50, Voronoi 0.43, surface area 0.34, surface deformation 0.52.
- "Spherical" phenotype: cell density drops sharply at radius ~150 μm; "spherical with far gaslikes" has non-zero density above 175 μm (Fig 2a).
- Voronoi volumes: sharp peak ~4000 μm³ + smaller peak ~2000 μm³ for ordered phenotypes; "disordered" spread over wide range (Fig 2c).
- Overall deviation score scales **approximately linearly** with transformation strength (noise up to level 200; deformation amplitude up to 120) — an emergent, not imposed, property (Fig 5, §2.5).
- Experimental: statistically significant separation between collagen concentrations only appears at **day 3** (Fig 7c); deviation lowest between 3 and 4 mg/ml, reproducing the **invasion transition between 2 and 3 mg/ml** reported by Kang et al. (single-cell "gas-like" invasion at 1–2 mg/ml vs collective "liquid-like" invasion at 3–4 mg/ml).
- Sim-vs-experiment: lowest deviation between 4 mg/ml collagen and "deformed" phenotype; highest between 2 mg/ml and "spherical" (Fig 8).
- Cell number ~doubles after 3 days; division probability set so cell number doubles after 250,000 MC steps; each cell allowed to divide once; division requires volume ≥ 90% of target.

## 5. Parameters & constants of interest
| Quantity | Value + units | Source in paper |
|---|---|---|
| Spheroid initial diameter (sim) | 200 μm | §4.2 |
| Spheroid diameter (experimental, MDA-MB-231) | 300–400 μm | §4.3 |
| Cells per spheroid | ~1000 (exp) / ~2000 (sim) | §4.2, §4.3 |
| Simulation box | (800 μm)³ | §4.2 |
| Collagen-I concentrations | 1, 2, 3, 4 mg/ml | §2.7, §4.3 |
| Invasion phase transition (collagen) | between 2 and 3 mg/ml | §2.7 (cf. Kang et al.) |
| Voronoi cell-volume mode (ordered cells) | ~4000 μm³ (2nd peak ~2000 μm³) | Fig 2c, §2.2 |
| Cell-density bulk edge ("spherical") | ~150 μm radius | Fig 2a |
| Gaslike thresholds (from Kang et al. data) | Dcrit = 125 μm, dcrit = 19 μm | §2.2 |
| MC-step ↔ real time | 1 MC step ≈ 1 s | Fig 1 caption, §4.2 |
| Time to phenotypic differentiation | ~3 days (cell number doubles) | §4.2, §4.3 |
| Spheroid formation time | ~48 h | §4.3 |

Note: the CPM coupling factors (λV, λS, adhesion matrix A) and "adhesion 50/450", "motility 100" are **dimensionless lattice-Hamiltonian** numbers, not transferable SI mechanism constants. The transferable values above are geometric/concentration anchors of the MDA-MB-231-in-collagen system.

## 6. Relevance to ffn_cellsim
This is primarily **(d) spheroid/multicellular-scale context** plus a thin **(b/e) methods touch**, and it is **off-target as a mechanism source** for a fine-grained single-cell simulator. Concretely:

- **Mechanism mismatch (why Low)**: The simulation engine is a **Cellular Potts Model + agent layer** — exactly the *lumped, coarse-grained* class ffn_cellsim's architectural principle rejects as the runtime mechanism. Cells are voxel-aggregates governed by an energy Hamiltonian; there are no explicit filaments, no motor heads, no clutch bonds, no Bell-Evans off-rates, no Hill force-velocity. ECM here is **rigid non-displaceable fibers with stochastic lattice degradation**, not explicit cross-linked collagen bonds. So none of the CiS mechanisms or constants port into the H.1→H.7 mechanistic stack.
- **Where it does touch the project**:
  - *Spheroid-scale validation context*: The platform's validation overlay target (memory: PI-exp = MCF7-spheroid-on-pV4D4/collagen-I; single-cell traction → multi-cell cohesion) is the **same experimental genre** — breast-cancer spheroids in collagen-I. This paper re-uses **Kang et al. (ref 13)** MDA-MB-231 + collagen-I (1–4 mg/ml) data, which is a sibling dataset to the project's intended overlay. The **2→3 mg/ml collagen invasion transition** and **300–400 μm / ~1000-cell** geometry are reusable far-field anchors when ffn_cellsim eventually scales to multicellular cohesion.
  - *Validation-oracle / scoring idea*: The **overall deviation score** (standardize feature distances → Wasserstein/Euclidean/MSE → SLSQP-weighted linear combination) is a reusable **sim-vs-experiment scoring recipe** for the project's "validation overlay" step if/when ffn_cellsim produces multicellular point clouds. It is a comparison/acceptance methodology, consistent with the project's "closed-form models = acceptance oracles, not runtime" philosophy — here the *score* (not a mechanism) is what could be borrowed.
  - *Numerics note (tangential)*: Wasserstein/earth-mover's distance, Voronoi tessellation, marching-cubes surface extraction are standard point-cloud analysis tools that could feed an ffn_cellsim multicell-output analyzer.
- **Bottom line**: Keep as context for the multicellular/EXTEND far horizon and as a candidate sim-vs-experiment scoring function; do **not** treat any CiS parameter or the CPM/ECM mechanism as a single-cell oracle or parameter source. It is several abstraction levels above ffn_cellsim's fine-grained particle/bond runtime.

## 7. Limitations & caveats
- **Static snapshots only** — features capture one time point; no cell-trajectory dynamics, velocity correlation, or persistence (authors flag this).
- **Single cell type** — no distinction between cell populations.
- **ECM is rigid** — fiber alignment is fixed and cannot be displaced by cells (only degraded), which the authors call "overly simplified"; real collagen remodeling/strain-stiffening is absent.
- **"Disordered" phenotype is an edge case** — surface-extraction features break when there is no solid bulk core; authors question whether "disordered" is biologically realistic.
- **Coarse-graining gap vs ffn_cellsim** — CPM voxel cells + dimensionless Hamiltonian weights are ~3 scale-levels above explicit filament/motor/clutch particles; 1 MC step ≈ 1 s is a loose mapping, not a physical integrator.
- Score is **use-case-dependent by construction** (the weights are calibrated to the four CiS phenotypes), so the absolute Di,j values are not universal.

## 8. Key figures / tables
- **Fig 1 + Table 1/2** — the four simulated phenotypes and the parameter space / phenotype-defining parameter differences (cell-ECM adhesion, ECM degradation period, motility).
- **Fig 2** — the three cell-based features (cell density distribution, gaslike distribution, Voronoi volume histogram) across phenotypes; gives the ~4000 μm³ Voronoi mode and ~150 μm bulk edge.
- **Fig 7** — experimental MDA-MB-231 spheroids in 1–4 mg/ml collagen across days 1–3; deviation score reproduces the 2→3 mg/ml invasion transition at day 3 (most project-relevant figure).
- **Fig 8 + Table 3** — sim-vs-experiment deviation-score matrix and the fitted feature weights λf.

## 9. Notable quotes / citable claims
- "there is currently no strategy for systematically comparing 3D structural data between two spheroids, which can be obtained in vitro using stacked multiphoton microscopy images" (§1).
- "One of the possible applications of the deviation score is to use it as a reliable objective function for fitting simulated spheroids to experimental data" (§3 Discussion).
- "by varying the collagen concentration between 1 and 4 mg/ml, one can tune the fiber density and overall mechanical properties of the collagen network surrounding each tumor spheroid ... single cell invasion in 1–2 mg/ml collagen (gas-like phase) and collective invasion in 3–4 mg/ml collagen (liquid-like invasion)" (§4.3, citing Kang et al.).
- "A single MC step corresponds to roughly 1 s of real time in the context of this study" (Fig 1 caption / §4.2).
