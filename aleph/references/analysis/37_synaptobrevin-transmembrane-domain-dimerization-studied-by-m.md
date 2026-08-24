---
id: 37_synaptobrevin-transmembrane-domain-dimerization-studied-by-m
paper_n: 37
title: "Synaptobrevin Transmembrane Domain Dimerization Studied by Multiscale Molecular Dynamics Simulations"
authors: "Han, J. et al."
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.06.049"
paper_type: particle-MD
ffn_relevance: Low
ffn_themes: [tangential, membrane, numerics/methods]
entities: [synaptobrevin-2, snare-complex, syntaxin-1, snap-25, popc-bilayer, transmembrane-helix, lipid-bilayer]
methods: [molecular-dynamics, coarse-grained, martini-force-field, atomistic-md, potential-of-mean-force, umbrella-sampling]
measurables: [bond-lifetime, free-energy, helix-crossing-angle, association-rate, rmsd, interaction-energy]
keywords: [synaptobrevin, snare, transmembrane-domain, dimerization, membrane-fusion, martini, coarse-grained-md, pmf, popc, helix-helix-association]
tags: ["#membrane-biophysics", "#snare", "#molecular-dynamics", "#coarse-grained", "#transmembrane-helix", "#tangential"]
has_transferable_params: false
---

# [37] Synaptobrevin Transmembrane Domain Dimerization Studied by Multiscale Molecular Dynamics Simulations

**Tags:** #membrane-biophysics #snare #molecular-dynamics #coarse-grained #transmembrane-helix #tangential

| Field | Value |
|---|---|
| Authors | Han, J., Pluhackova, K., Wassenaar, T.A., Böckmann, R.A. |
| Year / Venue | 2015 / Biophysical Journal 109(4):760–771 |
| DOI / ID | 10.1016/j.bpj.2015.06.049 |
| Type | particle-MD (coarse-grained + atomistic molecular dynamics) |
| Pages | 12 |
| ffn_cellsim relevance | Low — molecular-scale membrane-protein structural biophysics; no cytoskeleton/ECM/adhesion content |

## 1. Summary
The authors use a multiscale MD pipeline (MARTINI coarse-grained for unbiased association, then atomistic CHARMM36 and AMBER/SLIPIDS for validation) to study how the transmembrane (TM) domain of synaptobrevin-2 (sybII), one of the SNARE fusion proteins, homodimerizes inside a POPC lipid bilayer. Using the DAFT high-throughput docking-assay protocol, they ran 500 unbiased 5-µs CG simulations per sequence (wild-type plus W89A/W90A, poly-Val, poly-Ile, poly-Leu mutants). The wild-type helix forms a stable right-handed dimer (crossing angle ~30° in CG, ~50° in atomistic for the dominant I-c1 interface) whose interface (L99, C103, L107, I111, plus M95/M96) matches a prior mutagenesis prediction; two additional binding interfaces (I-c2, III) are found and proposed to enable the experimentally observed higher-order oligomers (trimer/tetramer). Mutations tune dimerization: poly-Leu lowers it (correlating with low fusogenicity), poly-Val raises it via a promiscuous, nonspecific interface, and the interfacial-tryptophan double mutant behaves like wild-type. The work links TM-sequence specificity to peptide aggregation and SNARE-mediated membrane fusion.

## 2. Problem & motivation
Synaptic-vesicle fusion is driven by the SNARE complex (SNAP-25, syntaxin-1, synaptobrevin-2). Beyond the soluble four-helix bundle, the single-pass TM domains of SNAREs were shown experimentally to dimerize/oligomerize in a sequence-specific way and to be needed for the hemifusion-to-fusion transition, yet the structural basis and functional role of sybII TMD self-association were poorly understood and partly contradictory across experiments (weak vs. substantial dimerization). The paper aims to resolve the sybII TM dimer geometry, find alternative interfaces that could explain higher-order oligomers, and connect sequence to dimerization propensity/kinetics.

## 3. Methods / model
- **Model class:** explicit-particle molecular dynamics of a peptide pair in a solvated lipid bilayer; no continuum or agent-based component.
- **Coarse-grained stage:** MARTINI force field; two idealized sybII TM helices (residues 85–116) placed ~4.5 nm apart with random in-plane rotation (unbiased start) in a POPC bilayer (~110 lipids, 1280 CG waters, ~2900 CG beads total). Setup via PYMOL → MARTINIZE → DAFT checkerboard layout → INSANE membrane builder → MARTINATE minimization/equilibration. 500 production runs × 5 µs each per sequence, NpT, 310 K (Bussi v-rescale thermostat), 1 bar (Berendsen barostat, semi-isotropic), 20 fs CG time step, GROMACS 4.5.2. One 5-µs CG run ≈ 17 h on 8 cores.
- **Atomistic stage:** selected CG dimers reverse-mapped with BACKWARD; 200-ns all-atom runs (>31,000 atoms) with CHARMM36 and with AMBER99sb-ILDN + SLIPIDS for force-field cross-check; PME electrostatics, Parrinello-Rahman barostat, Nosé-Hoover thermostat, LINCS constraints, 2 fs step. One 200-ns AA run ≈ 84 h on 96 cores.
- **Analysis/readouts:** PMF(-like) profiles vs. inter-helical COM distance (population-density method, 0.05 nm bins); helix crossing-angle distributions; packing-density maps of helix B around helix A; 2D orientational distributions (binding position β, binding phase φ, relative binding position ω); per-residue van der Waals + Coulomb interaction-energy decomposition; backbone RMSD; first-order association rate constants from a 0–2 µs linear fit.

## 4. Key results (quantitative)
- Dimers formed in nearly all runs: e.g. 491/500 compact dimers for WT, 490/500 WWAA, 500/500 PolyV, 434/497 PolyL (Table 1).
- Relative association rate constants (WT = 100%): WWAA 75.4%, PolyV 190.7%, PolyL 55.8% (Table 1, 0–2 µs first-order fit).
- RH-dimer crossing angle ~30° in CG (Fig. 4 A), consistent with a prior simulated-annealing prediction (ref 19); WWAA shifts the RH distribution +8° to a larger angle.
- Atomistic crossing angles for WT interfaces: I-c1 ~50°, I-c2 ~30°, region III ~19° (CHARMM36; Fig. 11 caption).
- Initial helix-helix contact typically at 100–400 ns; no dissociation observed over ~10 µs total per WT (text, Results).
- PMF bound-state minima: WT ~0.75 and ~1.0 nm COM (two interfaces at 500 ns; the 1.0 nm minimum vanishes after 3 µs); WWAA minimum ~0.6 nm (closer packing); PolyV/PolyL minima ~1.0 nm (looser packing) (Figs. 5–6).
- WT dimer is right-handed at ~98% of final ensemble; PolyV LH-dimer population reaches ~39% (well-maintained), WWAA shows no LH dimer at 5 µs (Results).
- Interface populations at 5 µs (Table 2): WT I-c1 48.6% / III 26.2% / I-c2 21.2%; PolyV I-c1 41.1% / III 11.0% / I-c2 2.7% / extra "E" 40.6%; PolyL I-c1 73.9% / III 6.6% / E 8.1%.
- Dimer backbone RMSD 2–4 Å (stable) for most configs; PolyI most flexible, up to ~6 Å (Fig. 9), in both force fields.
- Dominant WT interface residues by interaction energy: L99, C103, L107, I111 (+ M95, M96) (Fig. 10) — matches mutagenesis predictions (refs 17,19).

## 5. Parameters & constants of interest
None transferable to ffn_cellsim. The quantitative outputs are molecular-scale membrane-protein observables (helix crossing angles, sub-nm COM PMF minima, ns–µs association times, kBT-scale interface energies) for a synaptic-fusion protein, with no cytoskeletal filament, motor, adhesion-clutch, ECM, or cell-mechanics constant. Listed below for completeness only; not usable as an ffn oracle.

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| WT RH-dimer crossing angle | ~30° (CG); ~50° I-c1 (AA) | Fig. 4 A; Fig. 11 caption |
| Inter-helix COM PMF minimum (WT) | ~0.75 nm | Fig. 5, text |
| Initial association time | 100–400 ns | Results, text |
| Relative dimerization rate (PolyV/PolyL) | 190.7% / 55.8% of WT | Table 1 |
| Dimer backbone RMSD | 2–4 Å (up to ~6 Å PolyI) | Fig. 9 |
| CG production timescale | 5 µs/run, 500 runs, POPC ~110 lipids | Methods |

## 6. Relevance to ffn_cellsim
This is **(f) tangential** for a mechanistic single-cell cytoskeleton+ECM simulator. The biology is synaptic-vesicle SNARE fusion at the molecular/membrane-protein scale — synaptobrevin TM-helix homodimerization in a POPC bilayer — which is entirely disjoint from ffn_cellsim's mechanisms (AFINES actin filaments, Hill/Stam-Hocky myosin motors, Bell-Evans adhesion clutches, catch-bond cadherins, Arp2/3 branching, collagen ECM). It maps to no H-unit on the Phase-1 chain (H.1 cortex → H.2 → H.3 → H.5 → H.7), and even the EXTEND membrane track (H.8) concerns whole-cell plasma-membrane mechanics (tension, bending, area), not specific transmembrane-protein dimer interfaces.

The only thin, secondary connection is **(e) numerics/methods**: the paper is a clean exemplar of the MARTINI coarse-grained MD + DAFT high-throughput unbiased-association workflow with reverse-mapping (BACKWARD) to atomistic validation, and of extracting PMF/free-energy profiles from population densities of many short unbiased runs rather than umbrella sampling. ffn_cellsim is HOOMD-based (not GROMACS/MARTINI) and operates at the mesoscale (~1000 effective filaments/cell), so this is conceptual cross-reference at most — the population-density-PMF idea and the "many unbiased short runs to build a free-energy/association profile" ensemble strategy loosely echo how one might validate binding-rate statistics, but nothing here ports as code or parameter. Not a validation oracle, not a parameter source, not a mechanism reference for this project.

Recommendation: index as Low / tangential; do not cite in any H-unit brief or validation gate.

## 7. Limitations & caveats
- **Scale mismatch:** molecular MD of a two-peptide system (~2900 CG beads); no organelle-, cell-, or tissue-scale mechanics — irrelevant to a single-cell mechanobiology engine.
- The CG "PMF" is strictly PMF-*like*: no dissociation events were observed on the 5-µs timescale, so the bound-state ensemble is not fully equilibrated and absolute binding free energies are not extracted.
- MARTINI's elastic-network/secondary-structure restraints fix helicity, so it cannot probe unfolding or sequence effects on backbone conformation; PolyL vs. PolyI are nearly indistinguishable in MARTINI (only slight bond-length differences), limiting mutant resolution.
- Results are for a specific protein (sybII) in a single-component POPC bilayer — no cytoskeleton, cortex, adhesion, or ECM coupling that ffn_cellsim cares about.

## 8. Key figures / tables
- **Table 1** — per-sequence monomer/C-dimer/compact-dimer counts (out of 500) and relative association rate constants; the headline kinetics comparison.
- **Fig. 4 / Fig. 6** — helix crossing-angle distribution (LH vs RH) and the 1D PMF/kinetics comparison across WT and mutants (free-energy + rate readout).
- **Fig. 10 / Fig. 11** — per-residue interaction-energy decomposition and the I-c1/I-c2/III dimer interface geometries (structural conclusion: L99/C103/L107/I111 + M95/M96).
- **Fig. 12** — modeled sybII homotrimer and tetramer built from the alternative dimer interfaces (the oligomerization argument).

## 9. Notable quotes / citable claims
- "The wild-type helix is shown to form a stable, right-handed dimer with the most populated helix-helix interface, including key residues predicted in a previous mutagenesis study." (Abstract)
- "In almost all simulations (464 out of 500 runs) for wt sybII, the dimer formation was initiated by C-terminal contacts between the helices, forming a V-shaped dimer." (Results, TM helix association)
- "Dissociation events leading to separation of the helices were not observed in the ensemble over the entire 5 µs of simulation time, i.e., no dissociation was seen for a total simulation time of ~10 ms [µs]." (Results)
- "This study provides molecular insight into the role of TM sequence specificity for peptide aggregation in membranes." (Abstract)
