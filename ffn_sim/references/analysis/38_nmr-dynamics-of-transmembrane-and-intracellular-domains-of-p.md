---
id: 38_nmr-dynamics-of-transmembrane-and-intracellular-domains-of-p
paper_n: 38
title: "NMR Dynamics of Transmembrane and Intracellular Domains of p75NTR in Lipid-Protein Nanodiscs"
authors: "Mineev et al."
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.07.009"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [tangential, membrane]
entities: [p75ntr, neurotrophin-receptor, death-domain, chopper-domain, transmembrane-domain, lipid-bilayer, nanodisc, dpc-micelle, dmpc, dmpg, membrane-scaffold-protein, ngf]
methods: [solution-nmr, nmr-relaxation, cell-free-expression, size-exclusion-chromatography, translational-diffusion-nmr, talos-prediction]
measurables: [rotational-correlation-time, transverse-relaxation-time, longitudinal-relaxation-time, heteronuclear-noe, hydrodynamic-radius, cross-correlated-relaxation-rate]
keywords: [p75ntr, neurotrophin-receptor, lipid-protein-nanodisc, solution-nmr, intrinsically-disordered, death-domain, chopper-domain, transmembrane-helix, membrane-mimetic, receptor-activation, snail-tong-model, msp]
tags: ["#structural-biology", "#solution-nmr", "#membrane-protein", "#nanodisc", "#intrinsically-disordered", "#neurotrophin-receptor", "#tangential"]
has_transferable_params: false
---

# [38] NMR Dynamics of Transmembrane and Intracellular Domains of p75NTR in Lipid-Protein Nanodiscs

**Tags:** #structural-biology #solution-nmr #membrane-protein #nanodisc #intrinsically-disordered #neurotrophin-receptor #tangential

| Field | Value |
|---|---|
| Authors | Mineev, Goncharuk, Kuzmichev, Vilar, Arseniev |
| Year / Venue | 2015 / Biophysical Journal 109(4):772–782 |
| DOI / ID | 10.1016/j.bpj.2015.07.009 |
| Type | experimental (solution NMR structural biology) |
| Pages | 11 |
| ffn_cellsim relevance | Low — single-receptor atomistic NMR dynamics; no cytoskeleton/ECM/adhesion mechanics |

## 1. Summary
The authors used solution NMR to characterize the spatial structure and backbone dynamics of the transmembrane domain (TMD), juxtamembrane "chopper" + linker region, and intracellular death domain (DD) of the p75 neurotrophin receptor (p75NTR, a type I single-span membrane protein in the TNFR family). Constructs (p75-TMDCD residues 245–308; p75-ΔECD residues 245–425; p75ICD 288–425) were reconstituted into detergent micelles (DPC), bicelles (DMPC/DHPC, DMPC/CHAPS), and lipid-protein nanodiscs (LPNs) of three membrane-scaffold-protein (MSP) sizes. Key findings: the chopper/linker region is intrinsically disordered and highly flexible (subnanosecond/picosecond motions) in every membrane mimetic; the folded DD's motions are completely uncoupled from the TM helix because of this long flexible linker; none of the intracellular domains showed stable membrane binding or self-association/dimerization even when several copies were crowded into one nanodisc under reducing conditions. The results are argued against the rigid "snail-tong" activation model (which requires DD–DD interaction coupled to TMD conformation) and in favor of alternatives (clustering/compartmentalization or proteolytic processing).

## 2. Problem & motivation
p75NTR is a drug target (Alzheimer's, neurodegeneration) but full-length structural data are sparse and contradictory. The mechanism connecting ligand binding at the extracellular domain to intracellular signaling (the "snail-tong" model, ligand-induced dimerization, etc.) is unresolved. The work asks: in a near-native lipid environment, are the intracellular domains rigid and coupled to the TMD (as the snail-tong model needs), and do they self-associate? It also serves as a methods demonstration that LPNs are a superior membrane mimetic for multidomain membrane proteins with large soluble domains.

## 3. Methods / model
- Class: experimental atomistic/molecular structural biology (no continuum or particle-mechanics model).
- Proteins from rat wild-type p75: expressed by continuous-exchange cell-free system (¹⁵N / ¹⁵N-¹³C labeling) or in E. coli BL21(DE3); purified on Ni-Sepharose, thrombin tag cleavage.
- Membrane mimetics: DPC micelles, DMPC/DHPC and DMPC/CHAPS bicelles, and LPNs assembled with MSP1 (≈124 kDa disc), MSP1D1ΔH5 (≈95 kDa), MSP1D1ΔH4H5H6 (≈52 kDa). Lipids DMPC/DMPG 4:1, POPC.
- Readouts: solution NMR on 600/800 MHz Bruker spectrometers, 30 °C, pH 6.5–7.2. Chemical-shift assignment via triple-resonance + NOESY 3D with nonuniform sampling and BEST-TROSY. Relaxation parameters (T1, T2, heteronuclear ¹H-¹⁵N NOE, cross-correlated rate η_xy) via pseudo-3D HSQC. LPN size from translational diffusion (stimulated-echo bipolar-gradient NMR). Secondary structure / mobility from chemical shifts via TALOS+. PDB reference 1NGR for the DD.
- Length/time scales: atomic (per-residue); dynamics on ps–ns (fast backbone motions) and µs–ms (slow conformational transitions). Particle sizes: nanodiscs ~4–12 nm hydrodynamic radius; extended ICD up to ~20 nm.

## 4. Key results (quantitative)
- Transmembrane helix spans L253–W276; chopper domain N277–I308 is disordered (Fig 1A). TM chemical shifts differ between micelle/bicelle, but chopper signals are invariant across all three mimetics → chopper does not interact with membrane (Figs 1C, S4).
- η_xy of chopper residues ≈ 0 within error (except N277–Q281 in micelles, N277–K283 in nanodiscs near the membrane surface) → disordered, subnanosecond motions (Fig 1B).
- DD rotational correlation time τ_C = 6.3 ± 1 ns in p75-ΔECD, equal within error to the isolated domain (6.5 ± 1 ns); NOE ≈ 0.9 → folded, stable DD uncoupled from TM helix (Figs 2, 3).
- Linker-chopper region: τ_C ≈ 2 ns, negative NOE, T2 ≈ 0.4 s, η_xy ≈ 2 s⁻¹ → extremely high mobility on a picosecond timescale (p.776).
- LPN hydrodynamic radii: empty/compact DMPG MSP1 LPNs diffuse ~4 nm; MSP1D1ΔH4H5H6 particles ~5.5 nm radius; crowded discs do not exceed ~10–12 nm (Figs 2C, 4B).
- Extended juxtamembrane region may reach ~20 nm, exceeding nanodisc diameter (p.776).
- Crowding 6 p75-ΔECD copies per LPN → effective DD concentration >1 mM, yet no stable DD dimerization; DD τ_C increased only ~1.5 ns in 6:1 vs 1:1 MSP1 discs (Fig 4). Seven DD residues show accelerated transverse relaxation (slow µs–ms motions): N-terminal site S339, L342, A394; C-terminal site E412, V413, I406, E349 (Fig 4D).
- DD structure identical to PDB 1NGR; BMRB depositions 25647, 25648, 25646.

## 5. Parameters & constants of interest
None transferable to a mechanistic cytoskeleton/ECM/adhesion simulator. The reported quantities (per-residue NMR relaxation times, rotational correlation times, nanodisc hydrodynamic radii) are properties of an isolated signaling-receptor's intracellular domains in a membrane mimetic, not of the actin cortex, motors, clutches, or ECM crosslinks that ffn_cellsim models.

| Quantity | Value + units | Source in paper |
|---|---|---|
| DD rotational correlation time τ_C | 6.3 ± 1 ns | p.776 / Fig 2D, 3 |
| Linker-chopper τ_C | ~2 ns | p.776 |
| Nanodisc hydrodynamic radius | ~4 nm (compact) to ~10–12 nm (crowded) | Figs 2C, 4B |
| Extended ICD length | up to ~20 nm | p.776 |
| (Not usable by ffn_cellsim — receptor biophysics, not cell mechanics) | — | — |

## 6. Relevance to ffn_cellsim
Category (f) tangential. This is an atomic-resolution membrane-protein structural biology / receptor-signaling paper. ffn_cellsim is a coarse-grained mechanistic particle simulator of the actin cytoskeleton, myosin motors, adhesion clutches, and ECM crosslinks at the ~µm/cell scale; the only sanctioned coarse-graining is the ×40 mesoscopic filament scale. p75NTR is a neurotrophin receptor with no role in the cortex/FA/motor/ECM machinery, and the paper's observables (per-residue NMR relaxation, intrinsic disorder of a juxtamembrane linker, nanodisc dynamics) operate ~3 orders of magnitude below the simulator's smallest particle and concern signaling-domain conformational dynamics, not mechanics. It is not a validation oracle, parameter source, mechanism reference for any H-unit (cortex H.1, FA/clutch H.4, motor H.3, ECM, membrane/nucleus/cytoplasm EXTEND H.8/9/10), nor a numerics reference. The single weak thematic adjacency is "membrane" — but this is membrane-mimetic NMR methodology (nanodiscs/MSPs), not the lipid-bilayer mechanics an EXTEND membrane unit would need (no bending modulus, tension, or area-elasticity values are reported). Recommendation: keep as a low-priority tangential reference; do not wire into any unit or oracle.

## 7. Limitations & caveats
- In vitro, no ligand (NGF) or NGF-binding domain, no cellular factors, no post-translational modifications — the authors explicitly note the snail-tong model "should not be completely dismissed" because of these absences.
- Could not assemble the C257 S–S-linked dimer of p75-ΔECD (the physiologically activating covalent dimer), so the key activated state was not directly observed.
- TMD and 5–10 adjacent juxtamembrane residues were not observed in LPNs (line broadening), so TMD conformation itself is uncharacterized here.
- Scale gap vs ffn_cellsim is fundamental: atomic per-residue dynamics in a nanodisc vs µm-scale particle/bond cell mechanics; no mechanical moduli or forces are reported.

## 8. Key figures / tables
- Fig 1 — chopper domain disorder in p75-TMDCD: secondary chemical shifts (A), τ_C from cross-correlated relaxation (B), invariant glycine HSQC region across three mimetics (C).
- Fig 3 — full per-residue NMR relaxation panel for p75-ΔECD in MSP1 LPNs (T1, T2, η_xy, NOE, α-helix propensity vs residue) showing uncoupled folded DD and disordered linker.
- Fig 4 — crowding series: DD signal intensity, LPN hydrodynamic radius, and τ_C vs copies-per-disc; slow-motion DD residues mapped onto PDB 1NGR (D).
- Fig 5 — four candidate p75NTR activation mechanisms (snail-tong, ligand-induced dimerization, compartment migration, proteolytic processing).

## 9. Notable quotes / citable claims
- "Our data reveal a high level of flexibility and disorder in the juxtamembrane chopper domain of p75NTR, which results in the motions of the receptor death domain being uncoupled from the motions of the transmembrane helix." (Abstract)
- "none of the intracellular domains of p75NTR demonstrated a propensity to interact with the membrane or to self-associate under the experimental conditions." (Abstract)
- "the snail-tong model fits the in vivo data well, but contradicts the structural information; therefore, we have to consider other possible mechanisms of p75NTR activation." (Discussion, p.781)
- "we presume that such information on the TMD can be obtained from samples with selectively protonated methyl groups or segmental isotope labeling." (Discussion, p.779)
