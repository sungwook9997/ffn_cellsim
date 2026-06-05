---
id: 36_thermodynamics-of-micelle-formation-and-membrane-fusion-modu
paper_n: 36
title: "Thermodynamics of Micelle Formation and Membrane Fusion Modulate Antimicrobial Lipopeptide Activity"
authors: "Lin & Grossfield"
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.07.011"
paper_type: particle-MD
ffn_relevance: Low
ffn_themes: [membrane, numerics/methods, tangential]
entities: [lipopeptide, c16-kggk, popc, pope, popg, micelle, lipid-bilayer, martini-cg]
methods: [molecular-dynamics, coarse-grained-martini, umbrella-sampling, weighted-histogram-analysis, hamiltonian-replica-exchange, steered-md, string-method, free-energy-pmf]
measurables: [binding-free-energy, free-energy-barrier, potential-of-mean-force, micelle-size-distribution, membrane-selectivity]
keywords: [antimicrobial-lipopeptide, amlp, micelle, membrane-binding, martini, umbrella-sampling, pmf, wham, bacterial-vs-mammalian-membrane, kinetic-selectivity, free-energy-barrier]
tags: ["#membrane-biophysics", "#coarse-grained-md", "#martini", "#free-energy-pmf", "#umbrella-sampling", "#antimicrobial-peptide", "#tangential"]
has_transferable_params: false
---

# [36] Thermodynamics of Micelle Formation and Membrane Fusion Modulate Antimicrobial Lipopeptide Activity

**Tags:** #membrane-biophysics #coarse-grained-md #martini #free-energy-pmf #umbrella-sampling #antimicrobial-peptide #tangential

| Field | Value |
|---|---|
| Authors | Lin & Grossfield |
| Year / Venue | 2015 / Biophysical Journal 109(4):750–759 |
| DOI / ID | 10.1016/j.bpj.2015.07.011 |
| Type | particle-MD (coarse-grained MARTINI, free-energy/umbrella-sampling) |
| Pages | 10 |
| ffn_cellsim relevance | Low — membrane-biophysics drug-design study; no cytoskeleton/ECM/single-cell-mechanics content |

## 1. Summary
The authors use coarse-grained (MARTINI 2.2P) molecular dynamics plus umbrella-sampling free-energy methods to study how the antimicrobial lipopeptide (AMLP) C16-KGGK forms micelles in water and how a 48-mer micelle binds to two model bilayers: an anionic bacterial-like membrane (POPE:POPG 2:1) and a neutral mammalian-like membrane (pure POPC). They introduce a novel reaction coordinate based on the number of hydrophobic contacts (C16–C16 and C16–lipid-tail) rather than the usual center-of-mass distance, run umbrella sampling along it, and reconstruct potentials of mean force (PMFs) with a custom high-precision WHAM. Key finding: micelle binding to membranes is thermodynamically favorable but, unlike the barrierless monomer case, presents a large free-energy barrier whose height depends strongly on membrane composition (~1.3 kcal/mol into POPE:POPG vs ~79 kcal/mol into POPC). This implies AMLP selectivity for bacterial membranes may be as much kinetic as thermodynamic, and that oligomeric (micellar) state in solution is a tunable design variable for antibiotic optimization.

## 2. Problem & motivation
Antimicrobial peptides/lipopeptides preferentially disrupt microbial membranes and resist evolved resistance, but their clinical use is limited by cost, size, protease susceptibility, and the poorly understood link between solution oligomerization and biological activity/selectivity. The question: (1) what is the equilibrium oligomer-size distribution of C16-KGGK in solution, and (2) does oligomerization (micellization) alter binding affinity and mechanism to bacterial vs mammalian membranes? Understanding this could rationalize peptide/lipopeptide drug optimization.

## 3. Methods / model
- Model class: coarse-grained MD with the MARTINI 2.2P force field (each CG bead ≈ 4 heavy atoms) and polarizable MARTINI water; GROMACS 4.6.3.
- Systems: bacterial-like bilayer = POPE:POPG 2:1; mammalian-like = pure POPC; 480 lipids (240/leaflet). C16-KGGK = palmitoyl chain + KGGK tetrapeptide (two lysines, two glycines), treated as random coil. Micelle binding systems = 51,300 CG particles (single 48-mer micelle placed ~60 Å from bilayer COM, 10:1 lipid/peptide); micelle-formation-in-water systems = 12,730 particles. ~100 mM NaCl.
- Reaction coordinate: novel "hydrophobic contact number" — a smooth sigmoidal contact function S_ij(r_ij) = 1/(1+(r_ij/r0)^n) summed over pairs (Eqs. 1–2), with harmonic umbrella restraints U = (k/2)(C_AB − C0_AB)^2 (Eq. 3). Neighbor list with cutoff Rcut updated every 5 steps.
- Sampling/analysis: umbrella sampling + WHAM (custom optimized C++ implementation "GWHAM" using Polak-Ribière conjugate gradient, Brent line search, MPFR multiple-precision library for numerical stability); Hamiltonian replica exchange (Gibbs sampling, attempted every 500 steps) for micelle formation; steered MD (SMD) to seed windows; string method for minimum free energy paths (MFEPs) on 2D PMFs.
- MD parameters: 20-fs time step, NPT (Nosé-Hoover thermostat 300 K, Parrinello-Rahman barostat 1 bar), Coulomb cutoff 12 Å, vdW switch 9 Å / cutoff 12 Å.
- Scale: nanometer lengths, hundreds-of-ns to microsecond aggregate sampling; total simulation time ~348 µs (micelle formation, 696 windows × 500 ns), ~1626 µs (POPE:POPG binding), ~1434 µs (POPC binding).

## 4. Key results (quantitative)
- Micelle-formation PMF (Fig. 1A) has three minima; the global minimum is NOT the 48-mer but a coexistence of a 17-mer + 31-mer; the 48-mer is metastable by ~9.0 kcal/mol; barrier to fuse 17+31 → 48 is ~22 kcal/mol. Solution is polydisperse (oligomers ~10–38 lipopeptides).
- ΔΔG_binding (bacterial − mammalian) ≈ −246.7 − (−19.3) ≈ −227.4 kcal/mol for the 48-mer micelle (Fig. 2); i.e. binding to anionic bacterial membrane is far more favorable.
- Per-lipopeptide binding: mammalian (POPC) ≈ −0.40 kcal/mol (less than k_BT); bacterial (POPE:POPG) ≈ −5.14 kcal/mol/molecule. (For comparison, isolated monomer binding to anionic membrane in prior work = −14.5 kcal/mol — micelle stability reduces the per-molecule affinity.)
- Membrane-entry free-energy barriers (Fig. 3, along MFEP): POPE:POPG ≈ 1.3 kcal/mol; POPC ≈ 79 kcal/mol. Barrier difference 77.7 kcal/mol ⇒ ~10^56 difference in binding rate (inferred from Arrhenius/Boltzmann scaling).
- Mechanism: into bacterial membrane the micelle flattens and inserts cooperatively, stabilized by Lys–PG-phosphate electrostatics (also drives POPG lateral demixing/clustering even before contact, Fig. S10); into mammalian membrane lipopeptides transfer one-at-a-time until the residual micelle is too small to shield acyl chains, then insert simultaneously (origin of the large POPC barrier).

## 5. Parameters & constants of interest
None transferable to ffn_cellsim. The reported numbers are CG-MARTINI free energies and barriers for an antimicrobial lipopeptide binding lipid bilayers — they characterize peptide-drug/membrane thermodynamics, not cytoskeletal filaments, motors, adhesion clutches, ECM cross-links, or cell-scale mechanics. The lipid species (POPC/POPE/POPG) and MARTINI mapping are membrane-chemistry details with no bearing on the platform's H.1–H.7 mechanistic units. (No transferable constants for a fine-grained single-cell cytoskeleton+ECM simulator.)

| Quantity | Value+units | Source-in-paper |
|---|---|---|
| (none usable as ffn_cellsim oracle/anchor) | — | — |

## 6. Relevance to ffn_cellsim  <-- MOST IMPORTANT
**Category: (f) tangential**, with a very thin (e) numerics/methods footnote.

This paper is off-topic for ffn_cellsim's mechanistic single-cell scope. ffn_cellsim models cytoskeletal filaments (actin/AFINES), motors (Hill/Stam-Hocky myosin), adhesion clutches (Bell-Evans/integrin), cadherin catch-bonds, Arp2/3 branching, and ECM cross-links as explicit HOOMD particles/bonds, validated against MCF7-spheroid traction/cohesion. This study is about a small antimicrobial lipopeptide drug binding model lipid bilayers — a membrane drug-design/biophysics problem. None of its entities (C16-KGGK, POPC/POPE/POPG micelles), measurables (binding/insertion free energies), or biological questions (bacterial-vs-mammalian membrane selectivity) map onto the cortex (H.1), FA/clutch (H.4), motor/myosin, ECM/collagen, or even the EXTEND membrane-track (H.8). The EXTEND membrane unit, if/when built, would concern the cell's own plasma membrane as a mechanical envelope (bending, area-compression, cortex coupling), not lipopeptide partitioning thermodynamics — so even the "membrane" theme tag is only nominal here, not a usable parameter source.

The only marginally transferable content is methodological: this group's free-energy machinery — umbrella sampling, a custom numerically-stable WHAM (GWHAM, MPFR multiple-precision), Hamiltonian replica exchange via Gibbs sampling, steered MD seeding, and the string method for minimum-free-energy paths on 2D PMFs — is a competent template if ffn_cellsim ever needs PMF/free-energy reconstruction along a collective variable (e.g. clutch engagement or filament-bundle coordinate). But ffn_cellsim's integrator is the Leimkuhler-Matthews BAOAB plugin on HOOMD with explicit mechanistic kinetics (Bell-Evans off-rates), not equilibrium free-energy sampling, so even this is a distant analogy rather than a directly reusable recipe. The "hydrophobic contact number" reaction coordinate is clever but specific to amphiphile aggregation. Net: keep as a Low-relevance, tangential reference; no oracle, no parameters, no mechanism for the platform.

## 7. Limitations & caveats
- Coarse-grained MARTINI: no atomic detail, no chirality (acknowledged; argued irrelevant for a 4-mer), smoothed PES makes kinetics ~4× faster (authors decline to rescale since focus is thermodynamics).
- Strong finite-size effects: the simulation cell is too small for true equilibrium oligomer distributions; adding a lipopeptide to the micelle artificially raises the entropic penalty by depleting the bath. They apply an analytical correction (~−24 kcal/mol favoring the larger aggregate) but note large uncertainty; cannot resolve mesoscopic fibrils seen experimentally.
- Single fixed 48-mer / 10:1 lipid:peptide chosen as the "representative" micelle; the membrane PMFs are conditional on this choice. Intermediate (20–30-mer) micelles could shift mammalian-membrane kinetics, unaddressed quantitatively.
- WHAM over ~600–800 windows with ~200–300 kcal/mol dynamic range required custom high-precision solvers — convergence is delicate.
- Scale gap vs ffn_cellsim is total: nm-scale lipid/peptide thermodynamics vs µm-scale cell mechanics; no overlapping observables.

## 8. Key figures / tables
- Fig. 1 (A,B): PMF of 48 C16-KGGK aggregating vs C16–C16 contact number; joint probability of cluster size vs contacts — shows three minima, global minimum = 17-mer + 31-mer coexistence (polydispersity).
- Fig. 2 (A,B): 2D PMFs (color maps) of micelle binding to POPE:POPG (A) vs POPC (B) as functions of C16–C16 and C16–lipid contacts, with the MFEP (black line) and snapshots of surface-bound (1) / transition (2) / inserted (3) states.
- Fig. 3: 1D PMFs along the MFEP — directly contrasts the ~1.3 kcal/mol bacterial-membrane barrier with the ~79 kcal/mol mammalian-membrane barrier (the paper's headline kinetic-selectivity result).
- Table S1 (Supporting): RC definitions, Eqs. 1–3 parameters, window counts (supplemental).

## 9. Notable quotes / citable claims
- "the binding of these AMLP micelles to membranes is thermodynamically favorable, but in contrast to the monomeric case, there are significant free energy barriers" (Abstract, p.750).
- "The barrier to entering a POPE:POPG bilayer is relatively small (1.3 kcal/mol)… in contrast to the barrier to enter a zwitterionic POPC bilayer (79 kcal/mol). The difference in barrier height is 77.7 kcal/mol, suggesting a difference in binding rates of 10^56." (p.756).
- "the most likely oligomerization state for 48 C16-KGGK molecules is the formation of a 17- and a 31-mer… more favorable than the 48-mer micelle… by ≈9.0 kcal/mol" (p.755).
- "controlling the oligomeric state in solution will vary the mechanism of binding and thus the binding kinetics in ways not readily predictable by considering the monomer alone." (Conclusions, p.757).
