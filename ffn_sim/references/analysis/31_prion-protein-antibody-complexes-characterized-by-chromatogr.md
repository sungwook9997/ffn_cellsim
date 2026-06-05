---
id: 31_prion-protein-antibody-complexes-characterized-by-chromatogr
paper_n: 31
title: "Prion Protein-Antibody Complexes Characterized by Chromatography-Coupled Small-Angle X-Ray Scattering"
authors: "Carter, Kim, Schneidman-Duhovny et al."
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.06.065"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [tangential]
entities: [prion-protein, fab-antibody, recombinant-mouse-prp]
methods: [small-angle-x-ray-scattering, size-exclusion-chromatography, integrative-modeling, monte-carlo-conformational-sampling]
measurables: [radius-of-gyration, maximum-dimension, molecular-weight, porod-volume]
keywords: [prion, prp, fplc-saxs, fab-fragment, radius-of-gyration, multi-state-modeling, foxs, multifoxs, integrative-structural-biology, protein-misfolding]
tags: ["#structural-biology", "#saxs", "#protein-misfolding", "#tangential", "#not-cell-mechanics"]
has_transferable_params: false
---

# [31] Prion Protein-Antibody Complexes Characterized by Chromatography-Coupled Small-Angle X-Ray Scattering

**Tags:** #structural-biology #saxs #protein-misfolding #tangential #not-cell-mechanics

| Field | Value |
|---|---|
| Authors | Carter, Kim, Schneidman-Duhovny et al. (Prusiner & Sali labs) |
| Year / Venue | 2015 / Biophysical Journal 109(4):793-805 |
| DOI / ID | 10.1016/j.bpj.2015.06.065 |
| Type | experimental (solution structural biology + integrative computational modeling) |
| Pages | 13 |
| ffn_cellsim relevance | Low — single-protein/antibody SAXS structure determination; no cell mechanics, cytoskeleton, ECM, or adhesion physics |

## 1. Summary
The authors use a fast-protein-liquid-chromatography-coupled small-angle X-ray scattering (FPLC-SAXS) pipeline at SSRL Beamline BL4-2 to obtain monodisperse solution-scattering data for recombinant mouse prion protein (full-length recPrP(23-230) and N-terminally truncated recPrP(89-230)) and for their complexes with two Fab antibody fragments (HuM-P/Fab-P and HuM-R1/Fab-R1) that recognize N- and C-terminal epitopes respectively. In-line SEC immediately before the X-ray beam suppresses the aggregation/non-monodispersity that defeated a conventional autosampler setup. SAXS profiles are interpreted by integrative atomic modeling (RRT conformational sampling of disordered regions, rigid-body antibody-antigen docking via PatchDock, FoXS/MultiFoXS single- and multi-state fitting). They find PrP populates ~60-75% compact "closed" and ~25-40% extended "open" conformational states, and propose two structural mechanisms by which the Fabs may block the PrP^C->PrP^Sc misfolding conversion (Fab-P clamps the conversion-prone 95-105 region; Fab-R1 additionally contacts the a2-b2 loop, residues 165-175).

## 2. Problem & motivation
Prion diseases arise from templated misfolding of cellular prion protein (PrP^C, mostly a-helical) into a b-sheet-rich, aggregation-prone, protease-resistant isoform (PrP^Sc). Anti-PrP antibodies clear PrP^Sc from infected cell cultures, but the structural mechanism is unknown. The flexible N-terminal domain of PrP and PrP^Sc itself resist crystallography/NMR, and amyloidogenic proteins aggregate at the high concentrations SAXS needs. The goal is to obtain artifact-free solution structural envelopes of PrP and PrP-Fab complexes to rationalize antibody-mediated inhibition of conversion.

## 3. Methods / model
- **Samples**: recombinant mouse recPrP(23-230), recPrP(89-230); Fab fragments HuM-P, HuM-R1; and the four recPrP-Fab complexes. Buffer 20 mM sodium acetate pH 5.1 (running) / Tris 50 mM NaCl reported in Table 1, 150 mM NaCl, 283 K.
- **Assay/readout**: FPLC-SAXS — Superdex PC3.2/200 SEC column (flow 0.05-0.08 mL/min) directly upstream of the X-ray capillary (1.5 mm quartz); Mar225 / Rayonix MX225-HE CCD; 11 keV beam, 1.7 m sample-detector distance, 1 s exposures, ~700 images per run; q range 0.01-0.5 A^-1; AgBe calibrant.
- **Data processing**: SasTool buffer subtraction, ATSAS/PRIMUS merging of 10-50 consecutive low-variance profiles near the UV peak, GNOM pair-distribution, SAXS MOW molecular-weight estimation, DATPOROD Porod volumes.
- **Computational modeling (IMP)**: disordered regions sampled with rapidly-exploring random trees (RRT) independently over phi/psi, 10,000 (validated against 50,000) models; Fab elbow angles (130-180 deg, 28 templates) modeled with MODELLER; recPrP-Fab complexes docked with PatchDock antibody-antigen protocol (>400,000 models), scored by SOAP-PP statistical potential; SAXS fit by FoXS (optimizing hydration-layer density c2 and excluded-volume c1) and chi-score; multi-state ensembles by MultiFoXS branch-and-bound enumeration (N=1-5 states).
- **Length/time scale**: nanometer single-protein/complex structure in dilute solution; no dynamics, no cellular scale.

## 4. Key results (quantitative)
- recPrP(89-230): Rg = 21.8 +/- 0.9 A, Dmax = 76.5 A, MW ~15.4-16.2 kDa, Porod volume 20,548 A^3 (Table 2).
- recPrP(23-230): Rg = 30.7 +/- 1.2 A, Dmax = 104.4 A, MW ~22.9-26.0 kDa, Porod volume 39,832 A^3 (Table 2).
- Fab-P: Rg = 26.9 +/- 0.1 A, Dmax = 84.8 A, ~47-48 kDa; Fab-R1: Rg = 26.3 +/- 1.1 A, Dmax = 96.8 A (Table 2, Fig 2 C,D).
- recPrP(89-230)-Fab-P: Rg = 36.2 +/- 0.6 A, Dmax 121.1 A; recPrP(23-230)-Fab-P: Rg = 38.9 +/- 0.6 A, Dmax 136.1 A (Table 2).
- recPrP-Fab-R1 complexes more compact: Rg ~34.0 A for both (Fig 2 G,H).
- For comparison, the NMR structure (PDB 2L39) of recPrP(120-231) gives Rg = 15.0 A (p.797) — much smaller, reflecting the disordered N-terminus missing from NMR.
- Conformational populations: first Rg peak (compact closed) ~60-75% of population; second peak (extended open) ~25-40% (Figs 3 E,F and 5 E,F).
- Best-fit chi scores: single-state recPrP RRT models 1.33-3.05 (89-230) and 1.26-5.70 (23-230); two-state models 1.19-1.22 and 1.01-1.20; Fab-P best model chi=1.30 (elbow 177 deg) vs crystal 2HH0 chi=1.79 (elbow 128 deg); Fab-R1 chi=1.65; complexes chi 1.13-1.51.
- 75% of recPrP-Fab-R1 interface residues match the known 221-230 epitope; remainder from a2-b2 loop (165-175); cluster Ca RMSD 2.1 A (Fig 6 E).

## 5. Parameters & constants of interest
None (no transferable constants). All quantities are single-protein solution-structure descriptors (Rg, Dmax, MW, Porod volume, SAXS chi-scores, antibody elbow angles) for prion protein and immunoglobulin Fab fragments. There are no cytoskeletal, motor, adhesion, ECM, membrane-mechanics, or cell-scale mechanical parameters. The FoXS hydration-layer (c2) and excluded-volume (c1) fit parameters are SAXS-profile fitting knobs, not physics inputs usable by a particle simulator.

## 6. Relevance to ffn_cellsim
**(f) Tangential / off-topic.** This is a structural-biology paper about determining the solution shape of a single misfolding-prone protein (prion PrP) and its antibody complexes by SAXS, plus the integrative-modeling toolchain (IMP/FoXS/MultiFoXS/PatchDock) used to fit those envelopes. ffn_cellsim is a fine-grained mechanistic simulator of single-cell mechanobiology — actin/myosin/clutch/ECM particles and bonds, traction, cortex tension, spreading. None of the H-units (H.1 cortex, H.2, H.3, H.5, H.7; EXTEND membrane/nucleus/cytoplasm) touch prion biology, antibody binding, amyloid conversion, or atomic-resolution SAXS structure determination.

There is no usable mechanism, oracle, parameter, or numerics overlap:
- Not a validation oracle (no mechanical closed-form like Bell-Evans/Hill/Stam-Hocky).
- Not a parameter source (no force constants, rates, moduli, lengths/times at the cell scale).
- Not a mechanism reference (no cytoskeleton, adhesion, or ECM physics).
- Not spheroid/multicellular context.
- The only faint methodological adjacency is the abstract notion of "multi-state ensembles fit to experimental data" and "excluded-volume in a scattering calc" — but FoXS's c1 excluded-volume is a SAXS form-factor correction, conceptually unrelated to ffn_cellsim's LJ excluded-volume in dynamics. The RRT/Monte-Carlo conformational sampling is a static structure-search algorithm, not a dynamics integrator like BAOAB.

Keep at Low and tag tangential. The paper landed in the cellpress_bundle by topic-keyword overlap ("protein", "antibody", "membrane-bound PrP^C"), not by mechanistic relevance.

## 7. Limitations & caveats
- SAXS gives only low-resolution (nanometer) shape/size, not atomic structure; models are under-determined and the authors lean on Occam's razor plus multi-state ensembles to avoid overfitting.
- Recombinant, non-glycosylated, GPI-anchorless PrP in dilute buffer — not the native membrane-tethered, glycosylated protein; PrP^Sc itself was never measured.
- Proposed inhibition mechanisms are explicitly "speculation" from static envelopes; no kinetics, no force, no dynamics.
- Entirely outside the length/time/physics regime of a particle-based cell-mechanics simulator: single protein vs ~1000-filament cell; equilibrium solution scattering vs Langevin dynamics with active stress.

## 8. Key figures / tables
- **Table 2** — full SAXS structural summary: Rg (real), Guinier points, q*Rg, Dmax, MW, MW estimation, Porod volume for all 8 samples. The one data-dense table.
- **Fig 1** — autosampler (aggregated, blue) vs FPLC-SAXS (linear Guinier, red) profiles, motivating the in-line SEC pipeline.
- **Fig 2** — per-sample SAXS profiles + Kratky/Porod-Debye flexibility plots + FPLC Rg/I0 traces for all 8 samples.
- **Fig 6** — recPrP-Fab-R1 complex model showing the a2-b2 loop (residues 165-175) contact, the proposed conversion-influencing site.

## 9. Notable quotes / citable claims
- "Misfolding of the mostly a-helical cellular prion protein (PrP^C) into a b-sheet-rich disease-causing isoform (PrP^Sc) is the key molecular event in the formation of PrP^Sc aggregates." (Abstract)
- "In-line measurements by fast protein liquid chromatography coupled with SAXS minimized data artifacts caused by a non-monodispersed sample." (Abstract)
- "For both recPrP constructs, the first Rg peak (~60-70% of the population) corresponded to compact closed conformations, whereas the second peak (~30-40%) corresponded to extended open conformations." (p.797-798)
- "Fab-P binds to the N-terminal region (residues 95-105) that is known to undergo structural conversion from PrP^C to PrP^Sc ... acting as a structural clamp on the epitope." (Discussion, p.802)
