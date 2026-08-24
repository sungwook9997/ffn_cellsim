---
id: 34_the-conformation-of-myosin-heads-in-relaxed-skeletal-muscle-
paper_n: 34
title: "The Conformation of Myosin Heads in Relaxed Skeletal Muscle: Implications for Myosin-Based Regulation"
authors: "Fusi, Huang & Irving"
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.06.038"
paper_type: experimental
ffn_relevance: Low
ffn_themes: [motor/myosin, tangential]
entities: [myosin-ii, regulatory-light-chain, thick-filament, thin-filament, actin, troponin, tropomyosin, myosin-binding-protein-c, blebbistatin, rabbit-psoas-muscle]
methods: [fluorescence-polarization, maximum-entropy-analysis, x-ray-diffraction, electron-microscopy, single-muscle-fiber-mechanics]
measurables: [order-parameter, transition-temperature, rlc-orientation-angle, isometric-force, fraction-ihm-heads]
keywords: [interacting-heads-motif, super-relaxed-state, myosin-regulatory-light-chain, thick-filament-regulation, lever-arm-orientation, blebbistatin, lattice-spacing, skeletal-muscle, bifunctional-rhodamine-probe, OFF-state]
tags: ["#motor-myosin", "#structural-biophysics", "#muscle-regulation", "#interacting-heads-motif", "#tangential"]
has_transferable_params: false
---

# [34] The Conformation of Myosin Heads in Relaxed Skeletal Muscle: Implications for Myosin-Based Regulation

**Tags:** #motor-myosin #structural-biophysics #muscle-regulation #interacting-heads-motif #tangential

| Field | Value |
|---|---|
| Authors | Fusi, Huang & Irving |
| Year / Venue | 2015 / Biophysical Journal 109(4):783–792 |
| DOI / ID | 10.1016/j.bpj.2015.06.038 |
| Type | experimental (structural biophysics of muscle) |
| Pages | 10 |
| ffn_cellsim relevance | Low — sarcomeric thick-filament myosin regulation; structural, not a transferable mechanism for AFINES non-muscle myosin-II minifilaments |

## 1. Summary
The authors measured the *in situ* orientation of the regulatory light chain (RLC) C-terminal lobe of myosin heads in relaxed, demembranated rabbit psoas (skeletal) muscle fibers using fluorescence polarization from bifunctional sulforhodamine (BSR) probes at four sites on the RLC. They report a temperature-dependent conformational equilibrium among three preferred RLC orientations (RX1, RX2, RX3). At low temperature (2.5 °C) the dominant state RX1 has the myosin lever arm roughly perpendicular to the filament axis (rigor-like); on warming, two roughly-parallel orientations (RX2, RX3) emerge that match the free and blocked heads of the interacting-heads motif (IHM) seen by EM in invertebrate thick filaments. The fraction of parallel (IHM-like) heads rises sigmoidally with temperature (transition ~19 °C), reaching ~70% near physiological temperature, with ~30% remaining in the perpendicular RX1 state. Lattice compression (5% dextran T500) and the myosin inhibitor blebbistatin (25 µM) each lower the transition temperature ~5 °C (additively, ~10 °C together) and stabilize the parallel/IHM states without raising their high-temperature ceiling. On Ca²⁺ activation all three ordered orientations vanish, replaced by a broad intermediate distribution. The work connects the IHM, the super-relaxed (SRX) low-ATP-turnover state, and thick-filament "OFF" regulation in vertebrate skeletal muscle.

## 2. Problem & motivation
Striated-muscle contraction is classically gated by the thin filament (Ca²⁺ → troponin/tropomyosin uncovering actin sites). But the thick filament also has a regulatory "OFF" state — the IHM, where each myosin's two heads fold back against the backbone, inhibiting ATPase and actin binding. The IHM is well characterized in smooth/invertebrate muscle and cardiac C-zones, but for *vertebrate skeletal* muscle there are no high-resolution thick-filament structures, only indirect x-ray evidence. The question: in relaxed skeletal fibers, what fraction of heads are in IHM-like conformations, what controls the equilibrium, and how does activation disrupt it — i.e., is there a thick-filament-based regulatory switch parallel to the thin-filament switch?

## 3. Methods / model
- **System:** demembranated (skinned) ~6 mm rabbit psoas muscle fiber segments, mounted between a force transducer and loudspeaker motor; sarcomere length 2.42 ± 0.05 µm, cross-section 4542 ± 1186 µm² (n=12); Peltier temperature control to ±0.1 °C.
- **Probe/label:** four double-cysteine chicken skeletal RLC mutants (E-helix D95C/V103C; G-helix E131C/A138C; H-helix K151C/T158C; FG K134C/T122C), cross-linked with bifunctional sulforhodamine (BSR) to give 1:1 conjugates; the E-helix dipole is ~parallel to the myosin lever arm. Exchanged into fibers by a mild protocol replacing ~26% of native RLC (39 ± 10 µM BSR-RLC; native RLC ~150 µM). Force recovery >95% confirms function preserved.
- **Readout:** polarized fluorescence intensities → 2nd- and 4th-rank order parameters ⟨P2⟩, ⟨P4⟩ of the dipole distribution vs. filament axis; ⟨P2di⟩ characterizes fast probe wobble. Temperature swept 2.5–33.0 °C in 3 °C steps in relaxing solution.
- **Orientation model:** an "EG frame" defined from the nucleotide-free chicken skeletal myosin S1 crystal structure (Rayment 1993, PDB-equivalent), angles β (E-helix vs. filament axis) and γ (twist about E-helix). A maximum-entropy (ME) algorithm computes the smoothest (β,γ) distribution consistent with the four probes' ⟨P2⟩,⟨P4⟩; plotted on a sphere (β=latitude, γ=longitude).
- **Perturbations:** 5% dextran T500 (osmotic lattice compression; restores intact-cell thick-thin spacing), 25 µM blebbistatin (myosin inhibitor), ionic strength 150 vs. 190 mM; comparison states: Ca²⁺ activation (temperature-jump protocol, pCa 4.7) and rigor.
- **Validation:** IHM reference is the 3DTP tarantula EM reconstruction (Alamo 2008); confocal immunofluorescence confirms uniform A-band exchange.

## 4. Key results (quantitative)
- ⟨P2⟩ vs. temperature is sigmoidal for each probe, **transition temperature ~19 °C** in standard relaxing solution (Fig. 2).
- **E-helix** ⟨P2⟩ increases by ~0.4 units from low to high T (lever arm tilts toward parallel-to-filament); FG also becomes more parallel; G and H become more perpendicular (Fig. 2).
- **5% dextran T500** lowers transition T by ~5 °C; **25 µM blebbistatin** by ~5 °C; both together ~10 °C; effects additive (Fig. 2).
- Ionic strength: transition T ~4 °C lower at 150 mM than 190 mM (Fig. S5).
- Three ME orientations: **RX1** equatorial (β=90°, γ=90°, ~20° dispersion); **RX2** near pole (β=165°, γ=70°); **RX3** (β=130°, γ=0°) (Fig. 3). IHM reference: free head β=158°,γ=60°; blocked head β=131°,γ=0° — RX2≈free, RX3≈blocked.
- **Fraction in RX2+RX3** (FRX(2+3)): ~35% at 2.5 °C → max ~70% at 25 °C, sigmoidal (Fig. 4). ~30% remain in perpendicular RX1.
- This ~70% parallel fraction matches the literature **SRX-state occupancy ~75%** of heads at physiological T in rabbit fibers (refs 13,14).
- Active isometric force ~46 ± 4% higher at 24.8 °C than 11.0 °C; max isometric force 224 ± 23 kPa at 11 °C (Methods).
- On activation: RX1/RX2/RX3 disappear → broad peak β≈110–115°, γ≈50–55°, ~30–53° dispersion (Fig. 3I, S7).
- Rigor: single broad peak (β=90°, γ=80°), near RX1; E-helix ⟨P2⟩ drops ~0.2 (11 °C) to ~0.5 (24.8 °C) (Fig. 2).
- Lattice geometry: dextran shrinks thick-thin center-to-center spacing 28→23 nm (surface-to-surface 17→12 nm; thin radius 3.5 nm, thick 7.5 nm).

## 5. Parameters & constants of interest
These are sarcomere/structural-biology numbers, not parameters a non-muscle AFINES myosin-II minifilament model would ingest, but listed for completeness:

| Quantity | Value + units | Source-in-paper |
|---|---|---|
| RLC conformational transition temperature | ~19 °C (relaxing) | Fig. 2 |
| IHM-like (parallel RLC) head fraction, physiological T | ~70% | Fig. 4 / Discussion |
| Perpendicular (RX1, non-IHM) head fraction | ~30% | Discussion |
| Native RLC concentration in skinned rabbit psoas | 150 µM | Methods (ref 24) |
| Exchanged BSR-RLC concentration | 39 ± 10 µM (~26% exchange) | Methods |
| Max isometric force (skinned fiber, 11 °C) | 224 ± 23 kPa | Methods |
| Sarcomere length / cross-section | 2.42 µm / 4542 µm² | Methods |
| Thick–thin filament center-center spacing (relaxed→dextran) | 28 → 23 nm | Discussion (ref 31) |
| Thin / thick filament radius | 3.5 / 7.5 nm | Discussion |
| Blebbistatin dose used | 25 µM | Methods |

(None of these are transferable as fine-grained simulator inputs; see §6.)

## 6. Relevance to ffn_cellsim
**Tangential / Low.** This is a structural-biophysics study of the *regulatory OFF state of sarcomeric (skeletal) myosin-II thick filaments* — the interacting-heads motif and its temperature/lattice/blebbistatin-dependent equilibrium. ffn_cellsim's myosin is **non-muscle myosin-II** modeled as Stam-Hocky bipolar minifilaments with Hill (1938) force-velocity heads and Bell-Evans off-rates on AFINES actin in a cortex/stress-fiber context — not the highly ordered, MyBP-C-gated, IHM-regulated *sarcomeric* assembly studied here. Mapping to ffn themes:
- (c) **mechanism reference — weak/indirect.** The paper documents an autoinhibited folded-back ("OFF") head conformation and a graded equilibrium between an ATPase-silent IHM/SRX state and an active-competent state. If a future EXTEND track ever wanted to model myosin-II *activation kinetics* or an inactive/sequestered minifilament pool (non-muscle myosin-II also forms 10S autoinhibited monomers), this is background context for the *concept* of a folded OFF state — but the specific geometry (β,γ angles, RLC lever-arm orientations) and the sarcomeric/MyBP-C/titin machinery are muscle-specific and not represented in a cortical AFINES model.
- (b) **parameter source — no.** No force-velocity curve, no duty ratio, no step size, no stiffness, no off-rate, no bond energy — i.e., none of the motor constants ffn_cellsim's Hill/Bell-Evans/Stam-Hocky oracles need. The numbers here (order parameters, transition temperatures, RLC angles, kPa forces per fiber cross-section) describe regulation geometry, not single-motor mechanics.
- (a) **validation oracle — no.** There is no closed-form motor law to use as an acceptance oracle.
- (d/e) Not spheroid-scale, not numerics/methods relevant.

Bottom line: useful only as *biological-context literature* on the existence of a regulated, autoinhibited myosin-II OFF state; off-topic for the current Phase 1 cortex/clutch/motor build, which uses constitutively active minifilaments. File under motor/myosin background, flagged tangential. Do NOT mine it for runtime constants.

## 7. Limitations & caveats
- **Skeletal sarcomeric system**, not the non-muscle/cortical context ffn_cellsim simulates — different myosin isoform, different filament architecture (43-nm repeat, C-zone, MyBP-C, titin), thin-filament Ca²⁺ regulation present.
- ME analysis cannot resolve the *relative* probabilities of two closely spaced orientations precisely; RX2 vs RX3 apparent inequality is interpretive.
- Probe reports only the RLC C-lobe (b,γ); it is insensitive to catalytic-domain conformation, radial/azimuthal head position — so "IHM-like" ≠ proven identical to the EM IHM.
- Demembranated fibers, unphosphorylated RLC, ~26% exchanged labeled RLC; in vivo phosphorylation/lattice differences modulate the equilibrium and are only partly mimicked (dextran).
- All measurements are static-equilibrium orientation distributions; no single-molecule kinetics (rates, forces) are extracted — exactly the quantities a fine-grained particle simulator would need.

## 8. Key figures / tables
- **Fig. 2** — ⟨P2⟩ vs. temperature (2.5–33 °C) for E/G/H/FG probes, ± dextran ± blebbistatin, with active and rigor reference points; defines the ~19 °C sigmoidal transition and the −5/−5/−10 °C shifts.
- **Fig. 3** — ME spherical contour maps of RLC (β,γ) distributions across temperature and conditions, overlaid with IHM free/blocked-head positions (3DTP); shows RX1→RX2/RX3 shift and loss of order on activation.
- **Fig. 4** — FRX(2+3) (parallel/IHM-like fraction) vs. temperature; sigmoidal rise to ~70% at 25 °C; the central quantitative result.
- **Table 1** — physiological solution compositions (relaxing/preactivating/activating/rigor); experimental conditions reference.

## 9. Notable quotes / citable claims
- "the majority of the myosin heads in relaxed muscle fibers at near-physiological temperature have their RLC regions in conformations similar to those in the IHM" (Discussion, p.788).
- "~70% of the myosin heads have parallel RLC (RX2 and RX3) orientations, similar to those in the IHM/SRX state, with the myosin heads folded back onto the thick filament surface in a conformation that inhibits the myosin ATPase and the interaction with the overlapping thin filament" (p.789).
- "The remaining 30% of the heads have an RLC orientation (RX1) perpendicular to the filament axis and inconsistent with the IHM state." (p.789).
- "blebbistatin … favors the closed conformation and stabilizes the IHM … the helical order of the thick filament, and the SRX state of myosin. Our results show that blebbistatin also stabilizes the parallel orientations of RLC in relaxed muscle fibers by lowering the transition temperature" (p.789).
