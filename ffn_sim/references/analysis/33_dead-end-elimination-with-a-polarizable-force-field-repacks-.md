---
id: 33_dead-end-elimination-with-a-polarizable-force-field-repacks-
paper_n: 33
title: "Dead-End Elimination with a Polarizable Force Field Repacks PCNA Structures"
authors: "LuCore et al."
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.06.062"
paper_type: methods-software
ffn_relevance: Low
ffn_themes: [numerics/methods, tangential]
entities: [pcna, dna, side-chain-rotamer, amoeba-force-field, ubiquitin, sumo]
methods: [dead-end-elimination, x-ray-crystallography, polarizable-force-field, molecular-mechanics, particle-mesh-ewald]
measurables: [potential-energy, r-free, molprobity-score, rotamer-energy, polarization-energy]
keywords: [dead-end-elimination, AMOEBA, polarizable-force-field, side-chain-repacking, PCNA, x-ray-refinement, many-body-energy, rotamer-library, Goldstein-criteria, MolProbity, protein-design]
tags: ["#structural-biology", "#crystallography-refinement", "#polarizable-force-field", "#many-body-energy", "#protein-side-chain-packing", "#off-topic-for-cellsim"]
has_transferable_params: false
---

# [33] Dead-End Elimination with a Polarizable Force Field Repacks PCNA Structures

**Tags:** #structural-biology #crystallography-refinement #polarizable-force-field #many-body-energy #protein-side-chain-packing #off-topic-for-cellsim

| Field | Value |
|---|---|
| Authors | LuCore et al. (S.D. LuCore, J.M. Litman, K.T. Powers, S. Gao, A.M. Lynn, W.T.A. Tollefson, T.D. Fenn, M.T. Washington, M.J. Schnieders) |
| Year / Venue | 2015 / Biophysical Journal 109(4) 816–826 |
| DOI / ID | 10.1016/j.bpj.2015.06.062 |
| Type | methods-software (computational structural biology) |
| Pages | 11 |
| ffn_cellsim relevance | Low — atomistic crystallographic side-chain repacking / polarizable force-field theory; no cytoskeleton, cell, or ECM mechanics |

## 1. Summary
The authors extend the dead-end elimination (DEE) algorithm — a deterministic combinatorial method for finding the global minimum energy conformation (GMEC) of protein side chains — from pairwise-additive energy functions to **many-body energy functions** (i.e., the polarizable AMOEBA force field, where polarization energy is not pairwise decomposable). They re-derive the classical DEE (Desmet) and the more-stringent Goldstein elimination criteria to include three-body and higher-order energy terms, then apply the method, truncated at three-body interactions, to seven mid-to-low-resolution (2.5–3.8 Å, mean 3.0 Å) PCNA (proliferating cell nuclear antigen) x-ray crystal structures using a hybrid x-ray + force-field target function in the open-source FORCE FIELD X (FFX) software. Many-body AMOEBA DEE produced higher-quality models than PDB_REDO, local minimization, or pairwise OPLS-AA/L DEE: it lowered R_free and pushed the mean MOLPROBITY score to atomic-resolution quality (1.25 Å, 100th percentile). The refined models gave new structural/mechanistic insight into how the E113G and G178S separation-of-function mutations destabilize the PCNA trimer.

## 2. Problem & motivation
Most PDB structures are at mid-to-low resolution where the diffraction data alone underdetermine side-chain conformations, so refinement leans on force-field "prior chemical knowledge." A brute-force search over discrete rotamers is combinatorially intractable; DEE solves it deterministically but, until this work, only for **pairwise-additive** force fields. Real molecular driving forces — the hydrophobic effect and electronic polarization — are inherently many-body (the strength of an interaction between two residues depends on their mutual environment), so pairwise force fields systematically misrepresent them. The paper asks: can DEE be made compatible with polarizable / many-body potentials, and does doing so measurably improve crystal-structure quality?

## 3. Methods / model
- **Model class**: atomistic protein side-chain packing as a discrete combinatorial optimization (rotamers) over a continuous polarizable molecular-mechanics energy surface, coupled to experimental x-ray density.
- **Energy expansion**: total energy decomposed into environment + self-energy E_self(r_i) + two-body E2(r_i,r_j) + three-body E3(r_i,r_j,r_k) + ... (Eqs. 1, 7–9). For pairwise force fields the ≥3-body terms vanish; for AMOEBA the three-body term captures mutual polarization (Fig. 2).
- **Elimination criteria**: classical DEE (Eq. 4–6) and Goldstein (Eq. 12–13) rotamer/rotamer-pair criteria are re-derived to arbitrary n-body order (Eqs. 10–11, 14–15); derivations in Supporting Material.
- **Force fields**: polarizable AMOEBA (atomic multipoles + induced dipoles) vs pairwise OPLS-AA/L baseline. Electrostatics via particle-mesh Ewald.
- **Target function**: hybrid maximum-likelihood E_tot = E_chem + w_A·E_x-ray (Eq. 16); real-space density maps via Read/Cowtan formalism, Catmull-Rom spline (t=0.25). AMOEBA: energy and density weighted equally; OPLS-AA/L: w_A = 2.
- **Workflow**: minimize → divide unit cell into 4 Å boxes with 3 Å overlap → per-box DEE using Richardson/Ponder-Richards rotamer library (augmented with the initial coordinates as an extra rotamer) → AMOEBA truncated after trimer (three-body) interactions → re-minimize → iterate. Self-energy pruning: rotamers >30 kcal/mol above the residue's self-energy minimum pruned before two/three-body calculation.
- **Software**: FORCE FIELD X (FFX), http://ffx.biochem.uiowa.edu. Run on the Iowa NEON cluster.
- **System**: PCNA — a ring-shaped homotrimer DNA-replication processivity clamp; seven structures including wild-type, ubiquitin/SUMO-modified, and mismatch-repair / TLS-deficient mutants.

## 4. Key results (quantitative)
- Input PCNA data: resolutions 2.5–3.8 Å, mean 2.96 Å; deposited mean MOLPROBITY ≈ 2.86 Å (68th percentile), mean R/R_free = 25.5/28.7 (Table 1).
- Many-body AMOEBA DEE: mean MOLPROBITY score improved to **1.25 Å (100th percentile)**; mean R_free reduced to 26.70 (Table 2). Improvement in R_free was 3.0 for AMOEBA DEE vs 2.5 for pairwise OPLS-AA/L DEE; near-zero for PDB_REDO (p.821).
- MOLPROBITY score improvement: ~1.0 from local minimization, 1.24 from pairwise OPLS-AA/L, **1.61 from many-body AMOEBA** rotamer optimization (p.821).
- Poor-rotamer percentage cut from a deposited mean of 8.4% to **1.7%** (AMOEBA DEE) vs 3.5% (OPLS-AA/L DEE) (Table 2; p.823). PDB_REDO mean poor-rotamer 6.6%.
- Relative to OPLS-AA/L pairwise DEE, AMOEBA DEE adds: R_free lower by 0.5, MOLPROBITY lower by 0.37, 1.8% fewer poor rotamers (p.822).
- AMOEBA / OPLS-AA/L DEE models favored over locally minimized baseline by **>200 kcal/mol** per structure in force-field energy (p.821).
- **Many-body energy magnitude (Table 3)**: truncating the AMOEBA expansion at pairwise neglects ~1 kcal/mol/residue; truncating at three-body reduces neglected energy by an order of magnitude to <0.1 kcal/mol/residue. Largest single three-body energy >10 kcal/mol; example three-body energy 1.55 kcal/mol (Fig. 2). A 7 Å induced-dipole vector = 1 Debye (Fig. 2 scale).
- Mechanistic finding: wild-type PCNA trimer more stable than monomer by **1667 kcal/mol** of AMOEBA energy; drops to **1424 kcal/mol** in the E113G mutant — rationalizing reduced trimer stability and loss of translesion-synthesis support (p.823).

## 5. Parameters & constants of interest
None (no transferable constants for a mechanistic single-cell cytoskeleton/ECM simulator). The energy figures here (kcal/mol/residue, R_free, MOLPROBITY, Debye dipoles) are atomistic crystallographic-refinement quantities specific to protein side-chain packing in PCNA and do not map onto coarse-grained filament/motor/clutch/ECM physics.

| Quantity | Value+units | Source-in-paper |
|---|---|---|
| Pairwise-truncation neglected energy | ~1 kcal/mol/residue | Table 3, p.822 |
| Three-body-truncation neglected energy | <0.1 kcal/mol/residue | Table 3, p.822 |
| Self-energy prune threshold | 30 kcal/mol above per-residue minimum | Methods, p.820 |
| PCNA trimer-vs-monomer stabilization (WT) | 1667 kcal/mol (AMOEBA) | p.823 |
| Induced-dipole display scale | 7 Å = 1 Debye | Fig. 2 |

(These are listed for completeness, not as oracle anchors — none feed ffn_cellsim physics.)

## 6. Relevance to ffn_cellsim  <-- MOST IMPORTANT
**(f) Tangential / off-topic.** This is an atomistic structural-biology methods paper: deterministic combinatorial side-chain repacking (DEE/Goldstein) generalized to a polarizable atomistic force field (AMOEBA), validated against x-ray crystallography of PCNA, a DNA-replication clamp. It contains **no** cytoskeletal filament, motor, adhesion-clutch, ECM, membrane, nucleus-as-mechanical-body, or cell-scale mechanics. None of the ffn_cellsim H-units (H.1 cortex, H.2, H.3, H.5, H.7; EXTEND H.8/9/10 membrane/nucleus/cytoplasm) or the MCF7-spheroid validation overlay are touched.

Two extremely weak, non-actionable adjacencies, recorded only for completeness:
1. **Numerics/methods (very loose)**: it is a particle-based, force-field-driven, energy-minimization method (atomic multipoles, induced dipoles, particle-mesh Ewald). ffn_cellsim's mechanistic stance ("every interaction explicit, no lumped proxies") philosophically rhymes with the paper's argument that many-body effects shouldn't be approximated as pairwise — but ffn_cellsim's coarse-grained ×40 filament particles, LJ excluded volume, and BAOAB Langevin dynamics share no algorithmic machinery with DEE rotamer enumeration or polarizable crystallographic refinement. There is nothing to port.
2. **Tooling**: FFX/AMOEBA are atomistic protein-MD/refinement tools, not relevant to a HOOMD-blue coarse-grained mechanobiology engine.

Bottom line: **Low relevance.** Useful only as a negative example / boundary marker in the RAG corpus (atomistic protein structure refinement ≠ mesoscale cell mechanics). It is not a validation oracle, parameter source, mechanism reference, or scale-context paper for this project.

## 7. Limitations & caveats
- **Scale mismatch**: operates at single-amino-acid / atomic resolution inside a single protein; ffn_cellsim's smallest explicit unit is a coarse-grained filament segment — ~9+ orders of magnitude apart in scope.
- The method itself truncates the many-body expansion at three-body (four-body+ "infeasible computational cost" for PCNA-sized systems), so even within its domain it is an approximation; individual neglected ≥4-body energies can exceed 10 kcal/mol.
- Backbone held essentially fixed (only local minimization); no large conformational change, no dynamics, no kinetics — purely a static GMEC packing problem.
- Model bias from phase-dependent target functions acknowledged but not removed.
- No mechanical readouts (no forces, moduli, traction, velocities) of any kind relevant to cell mechanobiology.

## 8. Key figures / tables
- **Table 2** (p.821–822): per-structure R/R_free, force-field ΔE, and MOLPROBITY metrics for PDB_REDO vs OPLS-AA/L(±DEE) vs AMOEBA(±DEE) — the central result; AMOEBA DEE wins on most metrics.
- **Table 3** (p.822): neglected higher-order energy when truncating at two- vs three-body — quantifies the ~1 vs <0.1 kcal/mol/residue many-body contribution; the core justification for the method.
- **Fig. 2** (p.818): schematic of self/two-body/three-body energy terms; induced-dipole vectors showing AMOEBA polarization response (3-body energy 1.55 kcal/mol example).
- **Fig. 3 / Fig. 4** (p.822–824): PCNA β-strand electrostatic networks and F_o–F_c density maps showing AMOEBA recovers additional (bifurcated / CαH···O) hydrogen bonds vs PDB_REDO.

## 9. Notable quotes / citable claims
- "Beginning from dead-end elimination, we derive the first algorithm, to our knowledge, capable of deterministic global repacking of side chains compatible with many-body energy functions." (Abstract, p.816)
- "important molecular driving forces such as the hydrophobic effect and electronic polarization, which are fundamentally many-body in nature, have been implicitly approximated or neglected entirely." (Intro, p.816)
- "truncation of the many-body expansion at pairwise interactions neglects ~1 kcal/mol/residue of interaction energy under the AMOEBA polarizable force field." (Results, p.822 / Table 3)
- "many-body DEE using AMOEBA reduced the percentage of poor rotamers to 1.7% while simultaneously improving overall MOLPROBITY score and lowering both Rfree and AMOEBA potential energy." (Conclusions, p.823)
