---
id: 49_mechanistic-investigation-of-the-arrhythmogenic-role-of-oxid
paper_n: 49
title: "Mechanistic Investigation of the Arrhythmogenic Role of Oxidized CaMKII in the Heart"
authors: "Foteinou, Greenstein & Winslow"
year: "2015"
venue: "Biophysical Journal"
doi: "10.1016/j.bpj.2015.06.064"
paper_type: continuum-model
ffn_relevance: Low
ffn_themes: [tangential]
entities: [camkii, l-type-ca-channel, ryanodine-receptor, sodium-channel, ncx, calmodulin, phospholamban, cardiomyocyte, ros, h2o2]
methods: [stochastic-simulation, ode-electrophysiology-model, markov-channel-gating, local-control-ca-release]
measurables: [action-potential-duration, ion-current, calcium-concentration, reversal-potential, ca-leak-flux, pacing-cycle-length]
keywords: [camkii, oxidative-stress, early-afterdepolarizations, cardiac-electrophysiology, excitation-contraction-coupling, ncx-reversal, late-sodium-current, l-type-calcium-channel, reactive-oxygen-species]
tags: ["#cardiac-electrophysiology", "#camkii", "#arrhythmia", "#off-topic", "#computational-cardiology", "#tangential"]
has_transferable_params: false
---

# [49] Mechanistic Investigation of the Arrhythmogenic Role of Oxidized CaMKII in the Heart

**Tags:** #cardiac-electrophysiology #camkii #arrhythmia #off-topic #computational-cardiology #tangential

| Field | Value |
|---|---|
| Authors | Foteinou, Greenstein & Winslow |
| Year / Venue | 2015 / Biophysical Journal 109(4):838–849 |
| DOI / ID | 10.1016/j.bpj.2015.06.064 |
| Type | continuum-model (stochastic cellular electrophysiology / ECC) |
| Pages | 12 |
| ffn_cellsim relevance | Low — cardiac ion-channel / Ca-signaling electrophysiology, no cytoskeleton/ECM/adhesion mechanics |

## 1. Summary
The authors build a stochastic computational model of cardiac CaMKII activation that includes BOTH the classical phosphorylation activation pathway AND a newly-identified oxidation (methionine-oxidation by ROS/H2O2) pathway, parameterized from CaMKIIδ (cardiac isoform) data. This CaMKII model is embedded into their prior local-control canine ventricular myocyte model, which stochastically simulates ~12,500 Ca2+ release units and CaMKII-mediated phosphorylation of L-type Ca2+ channels (LCCs), ryanodine receptors (RyRs), and Na+ channels. The model reproduces the experimentally observed slow-pacing-rate dependence of H2O2-induced early afterdepolarizations (EADs): EADs appear on every action potential (AP) at PCL = 6 s but are absent at PCL = 1 s under 200 µM H2O2. The central mechanistic finding is that EADs require SYNERGISTIC activation of both late Na+ current (INaL) and L-type Ca2+ current (ICaL) — neither alone suffices — and that the Na+/Ca2+ exchanger (NCX) is a key amplifier: EAD emergence correlates with a shift in the timing of NCX current reversal (reverse→forward mode) from the repolarization phase toward the AP plateau. Blocking NCX abolishes EADs.

## 2. Problem & motivation
Oxidative stress (high ROS) drives cardiac disease and arrhythmias, and ROS can directly oxidize CaMKII to render it persistently active. Because CaMKII phosphorylates many excitation-contraction-coupling (ECC) proteins, its integrative effect on whole-myocyte electrophysiology under oxidative stress is hard to dissect experimentally. The goal is a quantitative model that links cellular ROS → oxidative CaMKII activation → ion-channel modulation → arrhythmogenic AP behavior (EADs), to identify which ionic changes most influence ROS-related arrhythmia propensity.

## 3. Methods / model
- **Model class:** stochastic cellular cardiac electrophysiology / ECC model (Markov-state channel gating + ODE membrane dynamics), NOT a spatial particle/mechanics model.
- **CaMKII module:** 12-subunit holoenzyme; each subunit transitions among states I (inactive), B (Ca2+/CaM-bound), P (autophosphorylated), OxB/OxA (oxidized), A (autonomous), OxP (oxidized+phosphorylated). Autophosphorylation restricted to adjacent subunits (holoenzyme geometry constraint). Built on Hashambhoy/Chiba frameworks; kox fit to Erickson et al. H2O2 dose-response; CaM affinity (KD) fit to Gaertner et al.
- **Whole-cell module:** local-control canine ventricular myocyte (Greenstein-Winslow), one 12-subunit CaMKII holoenzyme tethered per LCC; CaMKII phosphorylation shifts LCC gating mode 1→mode 2, increases RyR sensitivity to [Ca2+]dyad, and increases INaL (via the a8 background→bursting transition rate in the Grandi Na+ channel model).
- **Scales:** single cardiomyocyte; ~12,500 Ca2+ release units; AP timescales (ms–s); pacing-cycle lengths 1–6 s; 60 consecutive APs per simulation, statistics over final 50 APs.
- **Analysis:** APD90 ("APD"), box plots, Kruskal-Wallis ANOVA for significance; in-silico interventions = NCX block, RyR block, NKA inhibition, selective single-target phosphorylation.

## 4. Key results (quantitative)
- Rate dependence: no EADs at PCL = 1 s; EADs on every AP at PCL = 6 s under 200 µM H2O2 (Fig 2 A,D).
- Max fraction of open LCCs in mode-2 gating: ~0.35% (control) → ~7% with H2O2 (Fig 4 A).
- Diastolic RyR Ca2+ leak (JRyR): 0.6 µM/s (control) → 8 µM/s with H2O2 (Fig 4 B); elsewhere stated 9.2 → 64.3 µM/s (control → 200 µM H2O2), ~15-fold less than experiment (Limitations).
- INaL augmentation at long PCL ≈ 0.45% of peak INa at PCL = 6 s (text, p.842).
- NCX block: diastolic [Ca2+]i 67 nM → 208 nM; systolic [Ca2+]i 0.85 µM → 1.08 µM; [Na+]i drops ~2 mM; APD ~2-fold reduced (Fig 6).
- At time of NCX reversal: INaL 0.37 → 0.53 pA/pF after NCX block (p<0.001); ICaL decreases after NCX block (Fig 7).
- EAD-associated NCX reversal potential lies in a narrow window ~ −10 to 0 mV ("EAD voltage range"), median ~ −2 ± 2 mV; stable APs (NCX block) shift to ~ +9 ± 2 mV (Fig 8, Fig 9).
- NKA inhibition (~50% reduction) → overall ~6 mM gain in [Na+]i (matches Wagner et al.; Fig S17).
- H2O2 levels: ~35 µM normal human blood, up to ~100× in oxidative stress; 200 µM H2O2 used as pathophysiological proxy (p.840).

## 5. Parameters & constants of interest
None (no transferable constants). All quantities are cardiac ion-channel currents, calcium concentrations, action-potential voltages, NCX reversal potentials, and CaMKII kinetic rates — none map to the cytoskeletal/adhesion/ECM mechanics of ffn_cellsim. The CaMKII catch-/autophosphorylation kinetics are biochemical-signaling rates, not the mechanical bond/motor/clutch kinetics (Bell-Evans, Hill, Stam-Hocky) the simulator uses.

## 6. Relevance to ffn_cellsim
**(f) Tangential / off-topic.** This is a cardiac-myocyte electrophysiology and Ca2+-signaling modeling paper. ffn_cellsim is a fine-grained mechanistic simulator of single-cell cytoskeletal *mechanics* — actin filaments, myosin-II motors, adhesion clutches, ECM cross-links, cortex, with HOOMD particle/bond dynamics. There is essentially no overlap:
- No cytoskeleton, no actin/myosin, no focal adhesions/integrins, no ECM/collagen, no membrane mechanics, no cell spreading, no spheroid context.
- The "mechanistic" framing here refers to ionic/biochemical mechanisms of arrhythmia, not the mechanical force-generating mechanisms ffn_cellsim resolves.
- CaMKII appears, but as a signaling kinase acting on ion channels — not as a mechanical/cytoskeletal regulator.
- One very loose methodological echo: it uses a stochastic state-transition framework constrained by holoenzyme geometry and fit to single-molecule/biochemical rate data, and it embeds local stochastic units (Ca2+ release units) inside a whole-cell model — conceptually analogous to ffn_cellsim's "explicit stochastic sub-units inside a single cell" philosophy. But this is at most a distant design-pattern parallel, not a usable oracle, parameter source, or mechanism reference.

Bottom line: not a validation oracle, not a parameter source, not a mechanism reference for any H-unit (cortex / FA-clutch / motor-myosin / ECM-collagen / membrane-nucleus-cytoplasm EXTEND / numerics). Safe to file as tangential; retain only so the RAG index can explicitly exclude it from cytoskeleton/mechanics queries.

## 7. Limitations & caveats
- Scale/domain gap: this is a point/compartmental electrophysiology model (ODE + Markov gating + local-control Ca2+ units), no spatial mechanics, no particle dynamics — orthogonal to a fine-grained mechanobiology particle simulator.
- The model itself notes it cannot reproduce DADs (lacks spatial [Ca2+] gradient / Ca2+-wave machinery) and underestimates SR Ca2+ leak ~15-fold vs experiment.
- Na+ channel CaMKII effect reduced to a single rate (a8) modifying INaL amplitude; other CaMKII-INa kinetic effects omitted.
- Species mixing: canine myocyte model parameterized with some rabbit/CaMKIIδ data.

## 8. Key figures / tables
- Fig 1: CaMKII stochastic state diagram (phosphorylation + oxidation pathways) with CaM-affinity and H2O2 dose-response fits.
- Fig 2: rate dependence of H2O2-induced EADs (none at PCL 1 s, every-beat at 6 s); INaL-driven APD rate adaptation.
- Fig 5: comparative single-target analysis — INaL-only and LCC-only insufficient; both together produce EADs (synergy).
- Fig 9: EAD emergence correlated with NCX reversal-potential shift into the −10 to 0 mV "EAD voltage range."

## 9. Notable quotes / citable claims
- "selective activation of INaL ... although it prolongs action potential duration, is not by itself sufficient to produce EADs ... However, EADs emerge upon simultaneous activation of both LCCs and Na+ channels." (Abstract / Results, p.838, 842)
- "the emergence of H2O2-induced EADs was correlated with a shift in the timing of NCX current reversal toward the plateau phase earlier in the action potential." (Abstract, p.838)
- "Blocking the INCX ablates the formation of EADs, confirming the experimental findings of Zhao et al." (Results, p.843)
- "the arrhythmogenic effects of ROS are clearly multifactorial." (Discussion/Limitations, p.848)
