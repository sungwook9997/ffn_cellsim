---
id: 43_through-thick-and-thin-interfilament-communication-in-nbsp-m
paper_n: 43
title: "Through Thick and Thin—Interfilament Communication in Muscle"
authors: "Woodhead & Craig"
year: "2015"
venue: "Biophysical Journal (New and Notable)"
doi: "10.1016/j.bpj.2015.07.019"
paper_type: review
ffn_relevance: Low
ffn_themes: [motor/myosin, tangential]
entities: [myosin-ii, actin, troponin, tropomyosin, mybp-c, regulatory-light-chain, thick-filament, thin-filament, sarcomere]
methods: [electron-microscopy, fluorescence-polarization, x-ray-diffraction, three-dimensional-reconstruction]
measurables: [lever-arm-orientation, atpase-rate, head-conformation-population]
keywords: [interacting-heads-motif, super-relaxed-state, sarcomere, myosin-regulatory-light-chain, thick-filament-activation, MyBP-C, sentinel-head, striated-muscle, calcium-activation]
tags: ["#muscle-contraction", "#myosin-structure", "#sarcomere", "#thick-filament-regulation", "#tangential", "#structural-biology"]
has_transferable_params: false
---

# [43] Through Thick and Thin—Interfilament Communication in Muscle

**Tags:** #muscle-contraction #myosin-structure #sarcomere #thick-filament-regulation #tangential #structural-biology

| Field | Value |
|---|---|
| Authors | John L. Woodhead & Roger Craig |
| Year / Venue | 2015 / Biophysical Journal (New and Notable commentary) |
| DOI / ID | 10.1016/j.bpj.2015.07.019 |
| Type | review (a "New & Notable" dispatch commenting on Fusi et al. 2015) |
| Pages | 3 |
| ffn_cellsim relevance | Low — sarcomeric striated-muscle myosin regulation, not non-muscle cytoskeletal mechanics |

## 1. Summary
This is a 3-page "New and Notable" commentary in *Biophysical Journal* that contextualizes a research article by Fusi, Huang & Irving (2015, same issue) on how thin-filament (actin) activation by calcium is communicated back to thick (myosin) filaments in striated muscle. The puzzle: in vertebrate muscle, myosin thick filaments are intrinsically calcium-insensitive (isolated thick filaments stay relaxed even in calcium), yet they switch on rapidly during a twitch. The relaxed state is characterized by the "interacting heads motif" (IHM), where two myosin heads pair asymmetrically — one "blocked" head occludes the other's actin-binding, producing a super-relaxed, low-ATPase ordered array. Fusi et al. used bifunctional fluorescent labels on the myosin regulatory light chain (RLC) to read lever-arm orientation by polarization in intact skinned rabbit skeletal fibers, finding three relaxed-state RLC orientations (RX1/RX2/RX3); RX1 (~30% of heads) is rigor-like / perpendicular to the filament axis, suggesting a population of "sentinel" heads that probe the thin filament for activation and feed the signal back to disorder the IHM array, cooperatively switching the thick filament on. The commentary frames this alongside MyBP-C bridging and EM data.

## 2. Problem & motivation
The "missing link" in muscle regulation: thin filaments carry the calcium switch (via troponin/tropomyosin), but thick filaments must also be activated, and vertebrate myosin is not directly calcium-sensitive. How is thin-filament activation transmitted to the thick filament to recruit myosin heads into the cross-bridge cycle? Proposed answers: (i) MyBP-C bridges in the C-zone, and (ii) "sentinel" myosin heads that sample actin and propagate a structural signal along helical IHM tracks. The motivation matters for muscle physiology / cardiomyopathy, not for non-muscle cell mechanics.

## 3. Methods / model
This is a literature commentary; no new data. It reviews methods from the papers it discusses:
- **Fluorescence polarization** of bifunctional probes bound to the myosin RLC, in intact, skinned rabbit skeletal muscle fibers at near-physiological temperature — reads lever-arm (RLC) orientation, distinguishes RX1/RX2/RX3 relaxed conformations.
- **Electron microscopy (EM) + 3D reconstruction** of isolated cardiac/invertebrate thick filaments — source of the IHM atomic model and C-zone/D-zone disorder maps.
- **X-ray diffraction** — cited for filament periodicities and earliest post-calcium structural change (tropomyosin shift).
Length scales are molecular/structural (Ångström-level filament lattice geometry). No simulation, no governing equations, no numerics.

## 4. Key results (quantitative)
- Axial spacing of IHM crowns ~145 Å; helical repeat close to 430 Å (p.665, refs 1,7).
- Three relaxed-state RLC lever-arm orientations identified: RX1, RX2, RX3 (Fusi et al., reported on p.665).
- RX2 ≈ free-head and RX3 ≈ blocked-head orientations of the IHM atomic model (lever arms nearly parallel to filament axis) (p.665).
- RX1 ≈ rigor-like, RLC roughly perpendicular to filament axis; estimated **~30% of myosin heads** carry this rigor-like orientation under relaxing, near-physiological conditions (p.665–666).
- "Every third crown" of heads in the cardiac C-zone is more disordered than the IHM crowns (EM 3D reconstruction, ref 5) — if all those had the rigor-like RLC orientation it could account for the ~30% (p.666).

## 5. Parameters & constants of interest
| Quantity | Value + units | Source in paper |
|---|---|---|
| IHM crown axial spacing | ~145 Å (14.5 nm) | p.665 |
| Thick-filament helical repeat | ~430 Å (43 nm) | p.665 |
| Rigor-like (RX1) head population, relaxed | ~30% of heads | p.665–666 |

These are sarcomeric thick-filament lattice geometries and a striated-muscle conformational-population estimate. They are **not transferable** to a non-muscle cytoskeletal/ECM particle simulator: ffn_cellsim's myosin-II is the Stam-Hocky non-muscle bipolar minifilament with Hill force-velocity, not a striated sarcomere lattice with IHM/super-relaxed regulation. has_transferable_params: false.

## 6. Relevance to ffn_cellsim  <-- MOST IMPORTANT
**Tangential / Low.** This paper is structural-biology muscle physiology: regulation of striated-muscle thick filaments via the interacting-heads motif, MyBP-C, and RLC lever-arm orientation. It maps loosely to the **motor/myosin** theme only in the most general sense (it is about myosin-II head mechanics), but the specific physics — calcium-troponin-tropomyosin thin-filament regulation, the super-relaxed/IHM ordered helical array, sarcomeric C-zone/D-zone organization — is **absent from and irrelevant to** ffn_cellsim's model:
- ffn_cellsim uses the **Stam-Hocky bipolar non-muscle myosin-II minifilament** with **Hill (1938) force-velocity** and **Bell-Evans** detachment as the runtime mechanism. There is no sarcomere, no thin-filament calcium gate, no IHM/super-relaxed regulation, no MyBP-C in the model.
- The paper offers **no force constants, no rate constants, no force-velocity data, no off-rates** — nothing that could seed a motor oracle or a binding parameter. The only numbers are sarcomere lattice spacings and a population fraction, none of which inform a non-muscle cortex/clutch simulator.
- It is **not** a validation oracle (no closed-form law), **not** a parameter source, **not** a numerics/methods reference, and **not** spheroid/multicellular context.

Classification: **(f) tangential.** It is a muscle-myosin regulation commentary; useful only as background on how heads can be ordered/auto-inhibited in muscle, which does not apply to the non-muscle cortical actomyosin ffn_cellsim simulates. Safe to file as low-priority context for the motor/myosin theme.

## 7. Limitations & caveats
- A commentary, not primary data — all quantitative claims are secondhand (Fusi et al. 2015 and cited EM/X-ray work).
- The fluorescence-polarization technique reports only RLC (lever-arm) orientation, **not** motor-domain conformation, and gives **no spatial location** of the differently-oriented heads (noted explicitly on p.666). So the "sentinel head" hypothesis is structurally suggestive, not localized.
- Cardiac vs skeletal differences flagged: relaxed skeletal thick/thin filaments do not show the actin-binding interactions seen in cardiac, so cross-system comparison warrants caution (p.666).
- Entirely scale-mismatched to a fine-grained non-muscle cell simulator: striated-muscle sarcomere lattice, no relevance to cortex/FA/ECM mechanics.

## 8. Key figures / tables
- **Figure 1** (only figure): schematic of two communication routes from thin-filament activation to thick filament — via MyBP-C bridges and via "sentinel" myosin heads — labeling actin (A), troponin (Tn), tropomyosin (Tm), myosin motor domain (MD), RLC, lever arm (LA), IHM, sentinel head (SH), MyBP-C. The conceptual core of the paper.
- No tables.

## 9. Notable quotes / citable claims
- "A characteristic head conformation, known as the interacting heads motif (IHM), is associated with the relaxed state of myosin in all muscle types studied so far" (p.665).
- "myosin in vertebrate muscles is not directly sensitive to calcium ... So how is calcium activation of the thin filaments communicated to the thick filaments in the myofibril? This missing link has been a puzzle for muscle researchers over several decades." (p.665).
- "it is estimated that ~30% of the myosin heads have this rigorlike RLC orientation" under near-physiological relaxing conditions (p.665).
- "these heads act as sensors that continually probe the thin filaments to detect whether actin binding sites have become available" — the sentinel-head hypothesis (p.665).
