---
id: 19_multiscale-modeling-of-spheroid-tumors-effect-of-nutrient-av
paper_n: 19
title: "Multiscale Modeling of Spheroid Tumors: Effect of Nutrient Availability on Tumor Evolution"
authors: "Jakob Rosenbauer et al. (Berghoff, Glazier, Schug)"
year: "2023"
venue: "Journal of Physical Chemistry B (J. Phys. Chem. B, 127, 3607-3615)"
doi: "10.1021/acs.jpcb.2c08114"
paper_type: agent-based
ffn_relevance: Low
ffn_themes: [spheroid-scale-context, tangential]
entities: [cancer-cell, healthy-cell, tumor-spheroid, oxygen, glucose, mcf7]
methods: [agent-based, cellular-potts, monte-carlo, reaction-diffusion]
measurables: [proliferation-rate, oxygen-tension, spheroid-radius, migration-persistence, invasion-speed]
keywords: [cellular-potts-model, tumor-evolution, cell-cell-adhesion, cell-motility, nutrient-fluctuation, spheroid, fitness-landscape, metropolis, CellsInSilico, low-adhesion-selection]
tags: ["#spheroid-mechanics", "#cellular-potts", "#tumor-evolution", "#multicellular-scale", "#tangential"]
has_transferable_params: false
---

# [19] Multiscale Modeling of Spheroid Tumors: Effect of Nutrient Availability on Tumor Evolution

**Tags:** #spheroid-mechanics #cellular-potts #tumor-evolution #multicellular-scale #tangential

| Field | Value |
|---|---|
| Authors | Jakob Rosenbauer, Marco Berghoff, James A. Glazier, Alexander Schug |
| Year / Venue | 2023 / J. Phys. Chem. B, 127, 3607-3615 (Onuchic Festschrift special issue) |
| DOI / ID | 10.1021/acs.jpcb.2c08114 |
| Type | agent-based (cellular Potts model / lattice-based multicellular) |
| Pages | 9 |
| ffn_cellsim relevance | Low — multicellular CPM tumor-evolution study; cells are lumped Potts "spins", no subcellular cytoskeleton/clutch mechanism that a fine-grained single-cell simulator could borrow |

## 1. Summary
The authors run large-scale 3D cellular Potts model (CPM) simulations of an avascular spheroid tumor embedded in healthy host tissue, using the parallelized CellsInSilico/NAStJA framework on a supercomputer. Each cancer cell has two heritable phenotype parameters — cell-cell adhesion and motility — that mutate by small increments at a 4% probability per division. They track the centroid and spread of the tumor's phenotype distribution (a 12x12 adhesion x motility matrix) over time under three nutrient regimes: (i) competition for space only, (ii) static radial nutrient gradient coupling division/death rates, and (iii) a spatiotemporally oscillating nutrient sink. The robust finding across all regimes is that purely mechanical + spatial constraints drive selection toward low-adhesion phenotypes (which mechanically sort to the nutrient-rich tumor surface), favoring invasion; no consistent selection along the motility axis appears. Nutrient-dependent division/death accelerates evolution, and periodic nutrient fluctuations produce a resonance: a distinct peak in evolutionary speed at fluctuation period T = 100-200 kMCS (~15-30 cell generations, i.e. ~15-30 days assuming a 24 h cycle).

## 2. Problem & motivation
Tumor progression can be framed as Darwinian evolution: heterogeneous cell phenotypes compete for limited resources, and the "fittest" subpopulation dominates. Fitness, however, is not intrinsic — it is confounded by microenvironmental position (nutrient access, mechanical pressure). The authors ask whether purely geometric/mechanical constraints (adhesion + motility, plus a position-dependent nutrient field) are sufficient to drive evolution toward malignant traits (low adhesion = detachment/invasion; high motility = metastasis), and how *unstable* (fluctuating) nutrient supply — a hallmark of real avascular tumors where solid stress collapses vessels — shifts the evolutionary trajectory and speed.

## 3. Methods / model
- **Model class**: Cellular Potts Model (CPM, Graner-Glazier 1992), lattice-based. One cell = a connected set of lattice voxels labeled by a unique integer spin. Dynamics = Metropolis-accepted spin flips minimizing a Hamiltonian.
- **Hamiltonian terms**: Volume, Surface, Adhesion (proportional to shared contact area between cell types — explicitly *not* focal-adhesion-quantized), Random motility (per-cell preferential direction reassigned every 100 MCS), and a Central potential pulling cells toward the box center (and the nutrient sink).
- **Nutrient field**: radially linear decay within a central sphere; the negative source can be static or moved on a circle (amplitude A=50, period T) for the dynamic case. Nutrient is a growth-limiting proxy for oxygen/glucose.
- **Division/death**: cells divide if age > 2 kMCS and V > 0.9 V0 (450 of 500 um^3); two nutrient couplings tested — (1) constant rates, (2) linear: division rate 0 -> ~0.005 increasing with nutrient, death rate 0.001 -> 0 decreasing with nutrient.
- **Mutation**: at each division a daughter has 4% chance to increment/decrement one of {adhesion, motility} by one discrete step.
- **Scales**: grid 200x200x200 voxels, 1 voxel ~ 1 um, periodic BC; cell volume V0 = 500 um^3; saturation ~16000 cells; runs to ~580 kMCS; N = 20-21 replicates; mean cell lifetime ~7 kMCS.
- **Framework / numerics**: CellsInSilico (Berghoff et al. 2020) on the NAStJA backend, parallelized for the JUWELS supercomputer; open source at gitlab.com/nastja/nastja.

## 4. Key results (quantitative)
- Across all three nutrient regimes the phenotype-distribution centroid moves directionally toward **low adhesion**; no significant trend on the motility axis (Fig 3, Fig 4, p.4-7).
- Phenotype space is a **12x12 matrix** of adhesion x motility combinations (144 phenotypes; Fig 2, p.4).
- Adhesion and motility each scanned over discrete values [0,10,20,...,110], range [0...2*TMC] with Metropolis temperature TMC = 55 (Table 1, p.3).
- Linear nutrient dependence of division/death **accelerates** evolution vs. the constant case (Fig 3 bottom, p.5); a threshold-based dependence instead *slows* it (Fig S2).
- Fluctuating nutrients produce a non-monotone evolution-speed vs. period curve: speed is depressed for fast fluctuations (T = 100 MCS - 50 kMCS), then **peaks at T = 100-200 kMCS** (max at ~200 kMCS), then relaxes to the constant-case level (Fig 5b, p.7). Significance by t-test vs. constant case, p < 0.04.
- The resonant period corresponds to **15-30 cell generations** (lifetime ~7 kMCS), ~15-30 days at a 24 h cell cycle (p.7).
- Distribution spread (phenotypes with >10 cells) increases for fast fluctuations but decreases above T = 5 kMCS, indicating more directed evolution (Fig 5a, p.6).

## 5. Parameters & constants of interest
None (no transferable physical constants for a fine-grained single-cell simulator). All "parameters" are dimensionless CPM Hamiltonian couplings (Metropolis temperature, adhesion/surface/volume coupling constants) and abstract division/death rates per Monte Carlo sweep — they are not SI-grounded mechanical or kinetic constants. The only physically anchored numbers are geometric: voxel ~1 um, cell volume V0 = 500 um^3, grid 200^3 um, saturation ~16000 cells (Table 1 / Results, p.3-4) — useful only as loose spheroid-scale context, not as mechanism oracles.

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Voxel size | ~1 um | Results, p.4 |
| Cell target volume V0 | 500 um^3 | Table 1 / Results, p.3-4 |
| Division volume threshold | 450 um^3 (0.9 V0) | Results, p.4 |
| Simulation box | 200x200x200 um | Table 1, p.3 |
| Cell saturation count | ~16000 cells | Results, p.4 |
| Mean cell lifetime | ~7 kMCS | p.7 / Fig S1 |
| Mutation probability | 4% per division | Methods/Results, p.2,4 |
| Metropolis temperature TMC | 55 (dimensionless) | Table 1, p.3 |
| Resonant fluctuation period | T = 100-200 kMCS (~15-30 generations) | Fig 5b / p.7 |

(All rate/coupling values are CPM-internal, not SI; do not import as oracle constants.)

## 6. Relevance to ffn_cellsim
This is **(d) spheroid/multicellular-scale context only**, and largely **(f) tangential** for a mechanistic single-cell simulator. The paper is explicitly the *opposite* modeling philosophy from ffn_cellsim: cells are coarse lumped Potts "spins" on a lattice with energy-functional adhesion, and the authors state outright that adhesion "is not limited or quantized by focal adhesion" — i.e. exactly the lumped abstraction that ffn_cellsim's architectural principle forbids at runtime (no Bell-Evans clutches, no actin/myosin, no cytoskeleton). There is nothing here to use as a mechanism reference for cortex (H.1), FA/clutch, motor/myosin, or ECM.

Where it has marginal value:
- **Spheroid-scale narrative context** for the project's multi-cell/cohesion overlay ambition (the MCF7-spheroid validation track): it documents the qualitative biology that low-adhesion phenotypes mechanically sort to the spheroid rim and gain a proliferation/invasion advantage, and that proliferation concentrates at the rim while death concentrates in the hypoxic core. That rim-proliferation / core-death radial structure is a useful coarse expectation when ffn_cellsim eventually composes single cells into a spheroid — but only as a high-level sanity narrative, not a quantitative oracle.
- **Not a validation oracle**: its outputs (evolution speed, phenotype-space centroid drift, resonance period in kMCS) are in Monte-Carlo-sweep time and abstract phenotype units with no SI mapping, so they cannot anchor any ffn_cellsim gate.
- **Not a parameter source**: no force constants, bond lifetimes, moduli, or traction values.
- **CPM as a method** is noted in CLAUDE.md only as the kind of paper-model wrapper ffn_cellsim deliberately rejects, so it is a useful contrast example but not a numerics reference.

Net: keep for spheroid-scale-context tagging and the rim-proliferation/core-death biology; do not mine for mechanism or parameters. Relevance = Low.

## 7. Limitations & caveats
- **Lumped cell representation**: no subcellular cytoskeleton, no explicit adhesion molecules, no mechanics beyond the CPM Hamiltonian — a fundamental scale and fidelity gap vs. a particle/bond simulator.
- **Nutrient is a pure geometric proxy** (radial linear field), not a diffusion-reaction solve tied to consumption; the authors acknowledge real vasculature/angiogenesis dynamics are not modeled and that moving localized sources would need larger domains.
- **Time/length units are simulation-internal** (MCS, voxels) with only loose physical mapping; rates are not kinetically calibrated.
- **Motility result is null/inconclusive** — the predicted fitness advantage for motile cells under fluctuation could not be confirmed; only true-random-walk motility was implemented (no persistent random walk).
- Periodic boundary conditions + a central confining potential are artificial constructs to keep the tumor centered, not physiological.

## 8. Key figures / tables
- **Figure 1 (p.3)**: central slice of the 3D spheroid colored by cancer phenotype clone (wedge-shaped clonal domains, rim proliferation / core death) and by nutrient availability — the core spatial-structure result.
- **Figure 3 (p.5)**: phenotype-centroid trajectories drifting toward low adhesion, plus evolution-speed time series comparing no-dependency vs. linear nutrient dependency.
- **Figure 5 (p.7)**: macroscopic readouts vs. fluctuation period T — (a) distribution width, (b) evolution speed showing the resonant peak at T = 100-200 kMCS with t-test p-values.
- **Table 1 (p.3)**: full CPM parameter set (adhesion/motility ranges, TMC=55, box 200^3, V0=500, coupling constants).

## 9. Notable quotes / citable claims
- "Regardless of nutrient availability, we find a fitness advantage of low-adhesion cells, which are favorable for tumor invasion." (Abstract, p.1)
- "The strength is not limited or quantized by focal adhesion but is only determined by the adhesion parameter between the cell types and the shared area." (Methods, Adhesion, p.3) — explicit statement of the lumped-adhesion abstraction ffn_cellsim rejects.
- "There is a critical time scale for fluctuation of nutrient availability that provides a distinct peak in maximal evolution speed, which we find to be between T = 100 kMCS and T = 200 kMCS. This time range is between 15 and 30 cell generations." (Conclusions, p.7)
- "Mechanical properties alone increase proliferation at the edge of the tumor and cell death in the center." (Conclusions, p.6)
