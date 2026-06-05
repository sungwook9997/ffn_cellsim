---
id: 32_protein-conformational-changes-are-detected-and-resolved-sit
paper_n: 32
title: "Protein Conformational Changes Are Detected and Resolved Site Specifically by Second-Harmonic Generation"
authors: "Moree, B. et al."
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.07.016"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [tangential, membrane, validation-oracle]
entities: [calmodulin, maltose-binding-protein, dihydrofolate-reductase, lipid-bilayer, supported-lipid-bilayer, shg-dye]
methods: [second-harmonic-generation, mass-spectrometry, x-ray-crystallography, supported-lipid-bilayer, fluorescence-microscopy]
measurables: [shg-intensity-change, conformational-change-magnitude, ec50, binding-affinity, dye-tilt-angle]
keywords: [second-harmonic-generation, protein-conformational-change, supported-lipid-bilayer, calmodulin, maltose-binding-protein, dihydrofolate-reductase, ligand-binding, site-specific-labeling, his-tag-tethering, drug-screening]
tags: ["#structural-biology", "#biophysical-methods", "#shg", "#supported-lipid-bilayer", "#conformational-change", "#tangential"]
has_transferable_params: false
---

# [32] Protein Conformational Changes Are Detected and Resolved Site Specifically by Second-Harmonic Generation

**Tags:** #structural-biology #biophysical-methods #shg #supported-lipid-bilayer #conformational-change #tangential

| Field | Value |
|---|---|
| Authors | Moree, B. et al. (Connell, Mortensen, Liu, Benkovic, Salafsky) |
| Year / Venue | 2015 / Biophysical Journal 109(4):806–815 |
| DOI / ID | 10.1016/j.bpj.2015.07.016 |
| Type | experimental (biophysical methods / structural biology) |
| Pages | 10 |
| ffn_cellsim relevance | Low — in vitro single-protein conformational spectroscopy; no cytoskeleton / cell-mechanics / ECM content |

## 1. Summary
The authors present a surface-based optical method for real-time detection of ligand-induced protein conformational changes in solution. Purified proteins are labeled with a second-harmonic-generation (SHG)-active dye (amine-reactive on lysines or thiol-reactive on cysteines), tethered via an N-terminal poly-histidine tag to a Ni-NTA supported lipid bilayer (SLB) on glass, and irradiated by a pulsed Ti:Sapphire fundamental beam through total internal reflection. Because SHG intensity depends sensitively on the net average tilt angle of the dye relative to the surface normal, a conformational change that reorients the labeled side chain produces a measurable change in SHG signal. They validate the technique on three well-characterized proteins — calmodulin (CaM), maltose-binding protein (MBP), and E. coli dihydrofolate reductase (DHFR), including an engineered single-cysteine M20C DHFR mutant — showing that different ligands produce distinct signal magnitudes and directions, correlate the observed motions with published x-ray crystal structures and mass-spec-identified labeling sites, and recover a published nanomolar EC50 for an inhibitor.

## 2. Problem & motivation
Proteins are dynamic and change conformation upon ligand binding; these changes underlie allostery, signal transduction, enzyme catalysis, and motor-protein movement. Existing tools (NMR, x-ray crystallography, AFM) for probing conformational landscapes in real time and physiological conditions require high protein amounts, soluble/non-aggregating samples, and substantial structural a-priori knowledge. SHG is intrinsically surface-selective and orientation-sensitive but had not been realized as an easily accessible, biocompatible platform. The goal is a broadly applicable, low-protein-requirement, real-time assay for ligand-induced conformational change suitable for drug screening / structure-activity-relationship (SAR) work.

## 3. Methods / model
- Class: experimental nonlinear-optics biophysical assay (not a computational model).
- Mechanism / readout: %ΔSH = (I_tmax − I_t0)/I_t0 (Eq. 1); SHG intensity ∝ net average orientation of the dye relative to surface normal (z axis).
- Platform: Ni-NTA supported lipid bilayer (SLB) formed from small unilamellar vesicles incubated with 1 mM NiCl2; His-tagged protein captured on the bilayer; 16-well silicone gasket (~14 µL/well, 384-well-plate geometry).
- Labeling: SHG1-SE (amine/NHS-ester, on lysines) and SHG2-maleimide (thiol, on cysteines); conjugation pH tunes which lysine dominates (e.g., MBP K15 vs K88).
- Instrument: Mira 900 Ti:Sapphire oscillator pumped by Millenia V DPSS; p-polarized fundamental into a Dove prism at 69° (critical angle for TIR), ~100 µm spot; photon-counting PMT, 1 s integration.
- Orthogonal validation: liquid-chromatography tandem MS to identify labeled residues; overlay onto published PDB crystal structures (e.g., CaM 1CFD/1CLL/1CDL; MBP 1ANF/1JW4; DHFR 1RB3/1RX1).
- Length/time scales: single ~10–50 kDa proteins; signal changes resolved over seconds-to-minutes after ligand injection (mass-transport-limited in wells).

## 4. Key results (quantitative)
- CaM + CBP (calmodulin-binding peptide) with Ca²⁺: +18.1% ± 0.6% SHG change vs +0.1% ± 0.4% for buffer (N = 9) (p.808, Fig. 2A).
- CaM + CBP without Ca²⁺: 3.1% ± 1.6% (essentially no change) (Fig. 2B); CaCl2 injection alone: −12.6% ± 1.4% (Fig. 2C).
- MBP (labeled pH 8.3) + 1 mM maltose: −33.8% ± 1.0%; buffer +0.12% ± 0.61%; lactose (non-binding stereoisomer) +0.44% ± 0.79% (Fig. 3A).
- MBP (labeled pH 7.5, K15-dominant) + maltose: +27.4% ± 0.61% — sign reversal vs pH 8.3 prep (Fig. 3B), attributed to switching dominant label from K88 to K15.
- DHFR amine-labeled + 1 mM MTX (methotrexate): −60.2% ± 1.26% vs −2.27% ± 0.48% buffer (N = 8) (p.810, Fig. 4A); +100 µM TMP: −33.9% ± 5.31%; MTX after TMP: −3.48% ± 1.31% (competition / blocked, Fig. 4B).
- DHFR cysteine-labeled (M20C) + MTX: −81.03% ± 2.4% vs −0.6% ± 0.4% buffer; TMP alone: −88.05% ± 1.5%; MTX after TMP: −1.9% ± 1.0%; natural substrate DHF: −90.2% ± 0.4% (Fig. 5A–C).
- TMP dose-response EC50: logEC50 = 8.05 ± 0.07 → EC50 = 8.9 nM, "in excellent agreement with published values" (Fig. S3) — demonstrates nanomolar-affinity sensitivity.
- DHFR hydride-transfer rates (kinetic QC of constructs): WT His-tag 205 ± 20 s⁻¹; M20C 180 ± 10 s⁻¹; non-His-tag WT 220 s⁻¹ at pH 7, 25 °C (Materials & Methods, p.807).
- Labeling degrees (dye/protein): CaM 1.0; MBP 1.3; amine-DHFR 1.1; cysteine-M20C-DHFR 1.0.

## 5. Parameters & constants of interest
None transferable to ffn_cellsim. The measured quantities are SHG signal percentages, single-protein binding affinities (EC50 8.9 nM for TMP–DHFR), and DHFR hydride-transfer rate constants — all properties of isolated soluble enzymes/messengers (CaM, MBP, DHFR) with no role in the cytoskeleton, cortex, focal adhesions, motors, or ECM that ffn_cellsim simulates.

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| TMP–DHFR EC50 | 8.9 nM (logEC50 8.05 ± 0.07) | Fig. S3 / p.811 |
| DHFR hydride-transfer rate (WT) | 205 ± 20 s⁻¹ (pH 7, 25 °C) | Materials & Methods p.807 |
| TIR critical angle (this instrument) | 69° | SHG instrumentation p.808 |
| Bilayer Ni-loading | 1 mM NiCl2, 30 min | Sample prep p.808 |

These are recorded for completeness only; none is an oracle or validation anchor for a single-cell mechanics simulator.

## 6. Relevance to ffn_cellsim  <-- MOST IMPORTANT
This paper is (f) tangential. It is an in vitro, single-molecule structural-biology / nonlinear-optics methods paper about detecting ligand-induced conformational changes of three soluble model proteins (CaM, MBP, DHFR) tethered to a supported lipid bilayer. None of its mechanisms, proteins, length/time scales, or measurables map onto ffn_cellsim's domain:
- No cytoskeletal filaments (actin/microtubules), motors (myosin-II/Hill), adhesion clutches (Bell-Evans/integrin/talin/vinculin), cadherins, Arp2/3 branching, cortex, or ECM/collagen — the H.1–H.7 chain and the membrane/nucleus/cytoplasm EXTEND track find no anchor here.
- The "membrane" element is a synthetic Ni-NTA SLB used purely as an immobilization/orientation substrate, not a cell plasma membrane with mechanics — so it does not even inform the H.8 membrane EXTEND in any mechanistic way.
- The SHG technique is a static-protein conformational reporter; it is not a force probe (no traction, no bond lifetime, no retrograde flow), so it offers nothing as a validation oracle for the MCF7-spheroid traction/cohesion overlays.
- No transferable physical constants, no force-velocity or force-dependent off-rate data, no continuum or particle model to cross-check the BAOAB integrator or LJ excluded-volume implementation.
The only thread connecting it to the project's vocabulary is that DHFR/CaM are referenced in passing as motor-/allostery-relevant exemplars in the intro (citing kinesin/myosin force generation, ref. 3) — but the paper itself studies none of that. Conclusion: keep as a tagged "tangential" entry; do not use as a parameter source, oracle, or mechanism reference. Useful only as an example of SHG-on-SLB as an experimental modality, which is outside ffn_cellsim's simulation scope.

## 7. Limitations & caveats
- Heterogeneous labeling: amine chemistry labels multiple lysines (MBP labeled at ≥7 sites; degree 1.0–1.3), so the reported signal is a population-weighted average dominated by the most-labeled residue — site assignment is inferential (MS + crystal-structure overlay), not direct.
- Signal magnitude depends on label position and protein orientation on the bilayer; sign can flip with conjugation pH (MBP), so the readout is qualitative/relative, not an absolute conformational metric.
- Mass-transport-limited kinetics in static wells (the authors note two CaM experiments would not reach the same endpoint in 2 min), so reported kinetics are apparatus-limited.
- Entirely in vitro on isolated, His-tagged, dye-modified proteins on a synthetic bilayer — no cellular, mechanical, or multi-protein-complex context. Scale gap vs a fine-grained particle simulator is total: no forces, no network, no coarse-grained filament correspondence.

## 8. Key figures / tables
- Fig. 1: schematic of the TIR-SHG experiment — His-tagged labeled protein on Ni-NTA SLB; SHG magnitude set by net dye tilt vs surface normal (z).
- Fig. 2 (A–E): CaM time courses + summary + crystal-structure overlay (apo 1CFD, Ca-bound 1CLL, Ca/CBP 1CDL) with MS-identified lysines K13/K94/K115.
- Fig. 3 (A–D): MBP maltose vs lactose specificity, pH-dependent sign reversal, and apo/maltose-bound crystal overlay (1JW4/1ANF) with K15/K88.
- Fig. 5 (A–D): cysteine-M20C-labeled DHFR with MTX/TMP/DHF (largest changes, up to −90.2%), M20 reorientation in DHFR holoenzyme structures (1RX1/1RB3).

## 9. Notable quotes / citable claims
- "the magnitude of the SHG signal is proportional to the net, average orientation of the dye label relative to the surface normal" (p.808).
- "the addition of 1 mM of maltose resulted in a rapid decrease of 33.8% ± 1.0%, whereas ... lactose resulted in a negligible change" — specificity control (p.809, Fig. 3A).
- "we determined its EC50 for binding to DHFR as logEC50 of 8.05 ± 0.07 (EC50 = 8.9 nM) ... in excellent agreement with published values ... SHG has the sensitivity to detect binders with even nanomolar affinity" (p.811).
- "no a priori structural knowledge or mutagenesis is required ... there are no restrictions on the size or type of protein that can be studied by SHG" — scope claim (Discussion, p.813).
