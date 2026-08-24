---
id: 41_anionic-lipids-modulate-the-activity-of-the-aquaglyceroporin
paper_n: 41
title: "Anionic Lipids Modulate the Activity of the Aquaglyceroporin GlpF"
authors: "Klein, Hellmann & Schneider"
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.06.063"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [membrane, tangential]
entities: [glpf, aquaglyceroporin, aquaporin, lipid-bilayer, phosphatidylethanolamine, phosphatidylglycerol, cardiolipin, phosphatidylcholine, phosphatidylserine, e-coli, liposome]
methods: [stopped-flow-light-scattering, proteoliposome-reconstitution, semi-native-sds-page, dynamic-light-scattering]
measurables: [channel-permeability, weighted-rate-constant, tetramer-fraction, bilayer-thickness, spontaneous-curvature, acyl-chain-length]
keywords: [anionic-lipids, aquaglyceroporin, glpf, lateral-pressure-profile, negative-surface-charge, lipid-protein-interaction, transmembrane-protein, ribitol-flux, e-coli-membrane, tetramer-stability]
tags: ["#membrane-biophysics", "#lipid-protein-interaction", "#aquaporin", "#transmembrane-channel", "#structural-biology", "#tangential"]
has_transferable_params: false
---

# [41] Anionic Lipids Modulate the Activity of the Aquaglyceroporin GlpF

**Tags:** #membrane-biophysics #lipid-protein-interaction #aquaporin #transmembrane-channel #structural-biology #tangential

| Field | Value |
|---|---|
| Authors | Klein, Hellmann & Schneider |
| Year / Venue | 2015 / Biophysical Journal 109(4):722–731 |
| DOI / ID | 10.1016/j.bpj.2015.06.063 |
| Type | experimental (in vitro reconstitution biophysics) |
| Pages | 10 |
| ffn_cellsim relevance | Low — single-protein lipid-channel biophysics; no cytoskeleton, adhesion, motor, ECM, or cell-scale mechanics |

## 1. Summary
The authors reconstituted the *E. coli* aquaglyceroporin GlpF (a polytopic, tetrameric transmembrane water/glycerol channel) into liposomes of systematically varied lipid composition and measured channel activity via stopped-flow light scattering of ribitol flux. They find that (i) GlpF requires a minimal acyl-chain length (inactive in diC14:1/diC16:1-PC, active and increasing from diC18:1 to diC22:1-PC); (ii) GlpF activity is **not** affected by changes in the membrane lateral pressure profile (increasing the curvature-stress lipid PE up to high mol% leaves per-protein activity unchanged); and (iii) GlpF activity is **strongly impaired by negatively charged (anionic) lipid headgroups** (PG, cardiolipin/CL, and PS) regardless of headgroup chemistry, i.e., the effect tracks negative surface-charge density. Notably GlpF is poorly active in lipid bilayers that most closely mimic the native *E. coli* membrane (EPL extract and a PE/PG/CL ternary mix), and most active in non-physiological PC bilayers. The proposed mechanism is that electrostatic interaction between basic protein residues and anionic lipid headgroups induces subtle conformational changes (possibly affecting the conserved selectivity-filter arginine or the half-helix dipoles) that inactivate the channel without destabilizing the tetramer.

## 2. Problem & motivation
How does the lipid bilayer environment regulate the activity of an embedded polytopic transmembrane protein? Prior MD simulations predicted GlpF should be *more* active in PE bilayers (via lateral-pressure / channel-radius effects). The authors test, in vitro, which membrane variable — acyl-chain length (hydrophobic thickness), lateral pressure (non-bilayer PE curvature stress), or headgroup charge — actually controls GlpF function, and why GlpF appears suboptimally tuned to its own native membrane.

## 3. Methods / model
- **Class:** in vitro single-protein biophysics; not a model/simulation paper (cites MD work but does not perform it).
- **System:** purified *E. coli* GlpF reconstituted into liposomes at lipid/protein molar ratio 400:1 by detergent (n-octyl-β-D-glucopyranoside) dialysis.
- **Lipids:** Δ9-cis monounsaturated PC of chain lengths diC14:1 → diC22:1; diC18:1-PE, -PG, -PS, tetraC18:1-cardiolipin; LysoPC; *E. coli* polar lipid (EPL) extract (67% PE / 23.2% PG / 9.8% CL); PE/PG/CL 70/20/10 ternary mix.
- **Activity readout:** SX20 stopped-flow; proteoliposomes mixed with hypertonic 600 mM ribitol; 90° light scattering at 600 nm tracks shrink/re-swell as ribitol enters via GlpF. Decay fit to double-exponential; a weighted rate constant kw = (A1/(A1+A2))·k1 + (A2/(A1+A2))·k2 (Eq. 1) is normalized to incorporated-protein amount (kw,rel.) then to diC18:1-PC (kw,norm.).
- **Oligomeric state / incorporation:** semi-native SDS-PAGE (preserves native tetramer), Coomassie + ImageJ quantification of tetramer/dimer/monomer fractions.
- **Liposome sizing:** dynamic light scattering, mean radius 37 ± 11 nm.
- **Length/time scales:** molecular/nanoscale (single channel, ~Å channel radius, ~nm bilayer, ~30 nm vesicles); millisecond-to-second stopped-flow kinetics. No cell-, tissue-, or cytoskeletal-scale content.

## 4. Key results (quantitative)
- Hydrophobic bilayer core thickness spans 29.6 → 45.5 Å (dP-P) from diC14:1-PC to diC22:1-PC (Table S1; p.3).
- GlpF inactive in diC14:1-PC / diC16:1-PC; activity increases monotonically from diC18:1-PC to diC22:1-PC (Fig. 1C; p.4).
- Tetramer fraction: 0.78 (diC18:1-PC), 0.77 (diC20:1-PC), 0.65 (diC22:1-PC); 0.77 (EPL), 0.80 (PE/PG/CL) (Figs. S1/S2; p.3).
- PE up to 50 mol% keeps tetramer ~0.80; 60 mol% → 0.62, 70 mol% → 0.54 (Figs. S1B/S2B; p.4–5). Per-incorporated-protein activity unaffected by PE → lateral pressure does not modulate activity (Fig. 2C; p.5).
- Anionic-lipid incorporation peaks bell-shaped: max GlpF incorporation at 40 mol% PG and 20 mol% CL; PC/PG increased incorporation ×1.3, PC/CL ×3.2 vs pure PC (Fig. 3A,D; p.5).
- Highest per-protein activity at low anionic content (≈10 mol% PG, ≈2 mol% CL); activity decreased already at ~20 mol% PG/PS or ~5 mol% CL (Figs. 3,4; p.5–6, Discussion p.7).
- Cited spontaneous monolayer curvatures: J0m(DOPE) = −0.399, J0m(POPE) = −0.316 nm⁻¹ (from ref. 47; p.7) — used to argue both lipids raise acyl-chain lateral pressure.
- Native *E. coli* inner-membrane lipid fractions quoted: PE ~70–75%, PG ~15–20%, CL ~5–10% (p.1).

## 5. Parameters & constants of interest
| Quantity | Value + units | Source-in-paper |
|---|---|---|
| Liposome mean radius | 37 ± 11 nm | p.3 (DLS) |
| Bilayer hydrophobic thickness range | 29.6–45.5 Å (dP-P) | Table S1, p.3 |
| DOPE / POPE spontaneous curvature | −0.399 / −0.316 nm⁻¹ | cited ref. 47, p.7 |
| *E. coli* inner-membrane lipid mix | PE 70–75% / PG 15–20% / CL 5–10% | p.1 |
| Lipid:protein reconstitution ratio | 400:1 (mol) | p.2 |
| GlpF extinction coeff. (280 nm) | 37,930 M⁻¹cm⁻¹ | p.2 |

These are membrane-lipid / bacterial-channel constants. **None are transferable to ffn_cellsim** as cytoskeleton, adhesion, motor, ECM, or cell-mechanics oracles. (The bilayer-thickness and spontaneous-curvature numbers concern molecular lipid packing, not the continuum/coarse membrane mechanics an EXTEND-track membrane module would use.)

## 6. Relevance to ffn_cellsim
**Category (f): tangential / off-topic.** This is single-protein channel biophysics of a *bacterial* aquaglyceroporin in artificial liposomes. ffn_cellsim is a fine-grained but **cytoskeleton-centric** single-cell mechanobiology simulator (actin filaments, myosin minifilaments, adhesion clutches, ECM cross-links, cortex/cell-spreading mechanics validated against mammalian MCF7 spheroid traction/cohesion). Nothing here maps to:
- the Phase-1 chain (H.1 cortex → H.2 → H.3 → H.5 → H.7),
- FA/clutch, motor/myosin, ECM/collagen mechanisms, or
- the validation overlays (mammalian spheroid traction/cohesion).

The single thread of contact is the **EXTEND-track membrane module (H.8)** — but even there this paper is not useful: it treats the membrane as a host for a transmembrane *transport* protein and is about water/glycerol permeability and lipid-headgroup electrostatics, not the mechanical (bending modulus, tension, area-elasticity) or cortex-coupling properties an ffn_cellsim plasma-membrane module would need. The water/glycerol channel function, anionic-lipid charge regulation, GlpF tetramer stability, and *E. coli* lipidomics are all irrelevant to a mammalian cytoskeleton+ECM particle simulator. It is **not** a validation oracle, parameter source, mechanism reference, multicellular-scale context, or numerics reference for this project.

Keep at Low; cataloged for completeness so the RAG index records that this Cell Press bundle item was reviewed and deemed out of scope.

## 7. Limitations & caveats
- In vitro reconstitution in synthetic liposomes; per-protein activity normalization depends on SDS-PAGE incorporation quantification (semi-native, indirect).
- Proposed inactivation mechanism (electrostatic perturbation of selectivity-filter arginine / half-helix dipoles) is explicitly speculative; "will likely require high-resolution structures in the presence of negatively charged lipids" (p.8).
- Bacterial channel; no mammalian-cell, cytoskeletal, or tissue-scale relevance.
- Even as a membrane-mechanics reference, the paper measures transport/charge effects, not bilayer mechanical moduli; scale and physics gap vs. a particle-based cell membrane.

## 8. Key figures / tables
- **Fig. 1 (p.3–4):** GlpF incorporation + activity vs acyl-chain length (diC14:1→diC22:1-PC) plus EPL and PE/PG/CL — establishes the minimal-chain-length requirement and low activity in native-mimic membranes.
- **Fig. 2 (p.5):** GlpF incorporation/activity vs increasing PE — per-protein activity flat ⇒ lateral pressure does not modulate activity.
- **Fig. 3 (p.6):** GlpF vs increasing PG and cardiolipin — bell-shaped incorporation, sharp activity loss with anionic charge.
- **Fig. 5 (p.8):** GlpF crystal structure (PDB 1FX8) marking glycerol-interacting (green) vs lipid-headgroup-interacting (red) residues.

## 9. Notable quotes / citable claims
- "GlpF activity was severely affected by negatively charged lipids regardless of the exact chemical nature of the lipid headgroup, whereas GlpF was not sensitive to changes in the lateral membrane pressure." (Abstract, p.1)
- "the negative charge itself, and not the PG and CL headgroup chemistry, had caused inactivation of GlpF." (p.6)
- "the GlpF activity appears to be resistant to changes in the lateral membrane pressure ... implies an enormous structural stability." (p.7)
- "our work illustrates a mechanism by which a cell could control the activity of an α-helical membrane protein merely by varying the negative charge density around the protein." (Conclusions, p.8)
