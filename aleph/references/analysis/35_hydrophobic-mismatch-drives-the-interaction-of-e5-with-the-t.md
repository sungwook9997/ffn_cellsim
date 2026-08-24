---
id: 35_hydrophobic-mismatch-drives-the-interaction-of-e5-with-the-t
paper_n: 35
title: "Hydrophobic Mismatch Drives the Interaction of E5 with the Transmembrane Segment of PDGF Receptor"
authors: "Windisch, D. et al."
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.07.022"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [membrane, tangential]
entities: [e5-protein, pdgf-receptor, transmembrane-helix, lipid-bilayer, phospholipid, bovine-papillomavirus]
methods: [solid-state-nmr, synchrotron-radiation-circular-dichroism, oriented-circular-dichroism, pisa-wheel-analysis, lineshape-simulation]
measurables: [helix-tilt-angle, bilayer-hydrophobic-thickness, order-parameter, aggregated-fraction, secondary-structure-content]
keywords: [hydrophobic-mismatch, transmembrane-helix, e5-oncoprotein, pdgfr, lipid-bilayer-thickness, helix-tilt, golgi-membrane, receptor-tyrosine-kinase, dimerization]
tags: ["#membrane-biophysics", "#structural-biology", "#transmembrane-helix", "#hydrophobic-mismatch", "#tangential", "#not-cytoskeleton"]
has_transferable_params: false
---

# [35] Hydrophobic Mismatch Drives the Interaction of E5 with the Transmembrane Segment of PDGF Receptor

**Tags:** #membrane-biophysics #structural-biology #transmembrane-helix #hydrophobic-mismatch #tangential #not-cytoskeleton

| Field | Value |
|---|---|
| Authors | Windisch, D. et al. (Ziegler, Grage, Bürck, Zeitler, Gor'kov, Ulrich) |
| Year / Venue | 2015 / Biophysical Journal 109(4) 737–749 |
| DOI / ID | 10.1016/j.bpj.2015.07.022 |
| Type | experimental (solid-state NMR + CD structural biophysics) |
| Pages | 13 |
| ffn_cellsim relevance | Low — single-protein transmembrane-helix structural biophysics; no cytoskeleton, FA, motor, ECM, or cell-scale mechanics |

## 1. Summary
The authors study how the 44-amino-acid bovine papillomavirus E5 oncoprotein (and its cysteine-free truncated variant ΔE5) inserts into lipid bilayers of varying thickness, alone and co-reconstituted with the transmembrane domain (TMD) of the PDGF receptor β. Using synchrotron-radiation CD (SRCD), oriented CD (OCD), and solid-state 1D/2D ¹⁵N-NMR (PISA-wheel analysis), they find E5 forms an unusually long (~26-residue) α-helix whose tilt angle and stability depend critically on membrane thickness. In very thick DNPC bilayers (38.2 Å) the helix is nearly upright (~12° tilt) with little aggregation; as bilayers thin (DErPC→DEiPC→DOPC) the helix tilts modestly (14°→17°→21°), then distorts/kinks, then aggregates severely. E5 cannot tilt enough to fully relieve the mismatch. Co-reconstituting E5 with the PDGFR-TMD — which has the same hydrophobic length and behaves identically vs. thickness — rescues E5 in thin (DOPC) membranes, driving it back into an upright, well-aligned transmembrane state. The conclusion: intrinsic hydrophobic mismatch in thin Golgi/ER membranes is the long-range driving force that brings E5 and PDGFR together into a closely packed, mutually aligned helix bundle, explaining ligand-independent receptor activation in the secretory pathway.

## 2. Problem & motivation
PDGFR is a receptor tyrosine kinase activated normally by ligand-induced dimerization, but the E5 viral oncoprotein activates it in a ligand-independent manner purely through transmembrane helix-helix interactions, leading to oncogenic transformation. The recognition mechanism is puzzling because almost any hydrophobic residue in the E5 TMD can be substituted without losing function — so sequence is not the primary code. The authors test the hypothesis that the *physical* property of hydrophobic length/mismatch (rather than a specific sequence motif) is what drives E5–PDGFR recognition and explains why activation happens in the thin membranes of the ER/Golgi.

## 3. Methods / model
Pure experimental structural biophysics — no continuum or particle model.
- **Constructs**: truncated ΔE5 (first 34 aa, cysteine-free, native first 34 residues) and PDGFR-TMD, expressed in E. coli as Trp-ΔLE fusions, uniformly ¹⁵N-labeled in M9 + (¹⁵NH₄)₂SO₄.
- **Membranes**: reconstituted in SUVs / macroscopically aligned bilayers at protein/lipid 1:50 (mol/mol); lipids span a hydrophobic-thickness ladder — DNPC (di-C24:1, 38.2 Å), DErPC (di-C22:1, 34.4 Å), DEiPC (di-C20:1, 30.6 Å), DOPC (di-C18:1, 26.8 Å), POPC, DMPC (di-C14:0, 25.6 Å); thicknesses per Marsh.
- **SRCD** (UV-CD12 beamline, ANKA): secondary-structure content via DichroWeb/CONTIN.
- **OCD**: tilt diagnostic from the orientation-sensitive 208 nm band (positive = upright TM, negative = surface/tilted).
- **Solid-state ¹⁵N-NMR**: 1D lineshape deconvolution into well-oriented / misaligned / aggregated fractions (using ³¹P-NMR of the lipid to fix alignment quality); 2D SAMMY separated-local-field spectra fitted to simulated PISA wheels to read tilt angle τ and azimuthal rotation ρ. ¹⁵N-CSA tensor = (53.6, 77.1, 218.5) ppm; ideal α-helix φ=−61°, ψ=−45°; max ¹H–¹⁵N dipolar coupling 8.8 kHz.
- **Mixed 15N/14N experiments**: ¹⁵N-PDGFR-TMD + ¹⁴N-ΔE5 (and reverse) to isolate each partner's orientation in the heterocomplex.
Length scale: single transmembrane α-helix (~26–28 residues, ~39 Å) in a ~25–38 Å bilayer. No time dynamics, no mechanics, no cell.

## 4. Key results (quantitative)
- ΔE5 helix tilt angle (1D-NMR, vs membrane normal): **12°±3° (DNPC) → 14°±2° (DErPC) → 17°±1° (DEiPC) → 21°±5° (DOPC)** (Table 2).
- Aggregated fraction f_agg of ΔE5 rises with thinning: **13%±7% (DNPC), 10%±3% (DErPC), 35%±5% (DEiPC), 77%±8% (DOPC)** (Table 2).
- Order parameter S_mol ≈ 0.95–0.97 (close to 1.0) in all lipids — no helix wobble (Table 2).
- 2D PISA-wheel tilt: 0–10°(+2°) in DNPC, 12°±2° DErPC, 15°±5° DEiPC (kink/bend); DOPC not determinable (powder) (Table 2, Fig 4).
- PDGFR-TMD tilt mirrors E5: **3° (DNPC), 10° (DErPC), 14° (DEiPC), 22° (DOPC)** but with far less aggregation (Table 2, Fig 3).
- Secondary structure of ΔE5: total α-helix **91% in DNPC, 89% in DErPC**; ideal-helix (αR) **73% ≈ 26 aa in DNPC** vs 82% helix in detergent micelles (Table 1, p.744).
- E5 TMD = 28 residues (24 hydrophobic + 4 anchoring), ideal 26-aa helix ≈ **39 Å**, matched to DNPC's 38.2 Å (p.744–745).
- Geometric expectation vs observed tilt (E5 tilts far less than needed): expected 28° (DErPC, obs 14°), 38° (DEiPC, obs 17°), 47° (DOPC, obs 21°); residual unrelieved mismatch ≈ **3.4 Å (DErPC), 6.8 Å (DEiPC), 9.6 Å (DOPC)** (p.745).
- Heterocomplex in thin DOPC: PDGFR rescues E5 to an upright **7°±1° tilt** (vs aggregated alone); PDGFR itself unchanged (12°±2°) (Table 3, Fig 5D).

## 5. Parameters & constants of interest
| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Bilayer hydrophobic thickness ladder | 25.6 Å (DMPC), 26.8 Å (DOPC), 30.6 Å (DEiPC), 34.4 Å (DErPC), 38.2 Å (DNPC) | p.740 (per Marsh) |
| E5 TM helix tilt vs normal | 12°→21° as bilayer thins (DNPC→DOPC) | Table 2 |
| E5 ideal-helix length | ~26 aa ≈ 39 Å | p.744 |
| Order parameter S_mol | 0.95–0.97 | Table 2 |
| ¹⁵N-CSA tensor (E5 powder) | 53.6 / 77.1 / 218.5 ppm | p.739 |
| Max ¹H–¹⁵N dipolar coupling | 8.8 kHz | p.739 |

These are membrane-protein structural parameters (helix tilt vs bilayer thickness, CSA tensor, dipolar coupling). **None are transferable to ffn_cellsim's cytoskeleton/FA/ECM mechanics.**

## 6. Relevance to ffn_cellsim
**(f) Tangential / off-topic.** This is a single-protein membrane-biophysics paper: it characterizes how one viral transmembrane α-helix and a receptor TMD respond to lipid-bilayer hydrophobic thickness via NMR/CD. It contains nothing about actin/myosin, focal adhesions/clutch, ECM/collagen, cell cortex, traction, spreading, or any cell-scale mechanics — the entire content of the ffn_cellsim mechanistic chain (H.1 cortex → H.2 → H.3 → H.5 → H.7, FA/clutch, motor/myosin, ECM, spheroid overlays).

The only conceivable touch-point is the **membrane EXTEND track (H.8)**: ffn_cellsim plans to add a membrane compartment additively. If that membrane were ever modeled at the level of bilayer thickness and embedded transmembrane proteins (it is not, per current plans — the membrane is a mechanical/elastic shell, not a molecular bilayer with explicit lipids), then hydrophobic-mismatch tilt data could in principle be a reference. But that is a different physics regime (sub-helix molecular structure) than the coarse-grained cell-mechanics membrane envisioned, and a coarse-grained particle simulator of the cytoskeleton gains no oracle, parameter, mechanism, or numeric method from this work. Recommendation: **archive as Low; do not wire into any H-unit.** Not a validation oracle, not a parameter source, not a mechanism reference for the cytoskeleton.

## 7. Limitations & caveats
- Model lipids (di-C20:1–C24:1) are non-physiological synthetic phospholipids chosen for accessible phase-transition temperatures; thickness ladder is an in vitro proxy for ER/Golgi, not real organelle membranes.
- Uses the truncated, cysteine-free ΔE5, not full-length E5 (justified by prior equivalence but still a model peptide); native dimer also stabilized by disulfides absent here.
- DOPC/POPC/DMPC data are degraded by aggregation, absorption flattening, and light scattering — secondary-structure deconvolution reliable only for DNPC/DErPC.
- Purely static/structural: no dynamics, kinetics, forces, or cell-level readouts. Scale gap to a fine-grained particle cell simulator is total — this is sub-nanometer helix orientation, not micron-scale cytoskeletal mechanics.

## 8. Key figures / tables
- **Table 2** — tilt angles and aggregated fractions of ΔE5 and PDGFR-TMD across the four lipids (the quantitative core).
- **Fig 6** — schematic model: E5, PDGFR-TMD, and the E5+PDGFR complex tilting/aggregating vs bilayer thickness; complex rescues E5 in thin (Golgi-like) membranes.
- **Fig 4** — 2D SAMMY/PISA-wheel spectra of ΔE5: coherent wheel (DNPC) → kink (DEiPC) → powder/aggregate (DOPC).
- **Table 1** — ΔE5 secondary-structure content (91% helix DNPC) by SRCD/CONTIN.

## 9. Notable quotes / citable claims
- "Hydrophobic mismatch between lipid membranes and integral proteins can be an important regulator of protein function." (Introduction, p.737)
- "The intrinsic hydrophobic mismatch of these two interaction partners drives them together. They seem to recognize each other by forming a closely packed bundle of mutually aligned transmembrane helices." (Abstract, p.737)
- "E5 can tilt to a certain extent, but then the protein aggregates because it is unable to tilt any more." (Discussion, p.745)
- "In thin bilayers, E5 binds to the tilted and monomeric TMD of PDGFR to counteract the hydrophobic mismatch... The need to avoid hydrophobic mismatch may thus be the underlying driving force for the long-range recognition and binding of E5 to the receptor." (Discussion, p.746)
