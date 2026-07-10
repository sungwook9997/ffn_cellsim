# SourceEvidence / KnowledgeClaim / Parameter / ValidationGate candidates — Kim + Miyazaki corpus (2026-07-10)

**Scope.** 77 unique papers (Taeyoon Kim lab, Makito Miyazaki lab, + close-neighbour active-matter /
division / durotaxis / reconstitution) from `ffn_sim/references/2026_07_10`. All are ingested into the
TAG content layer (`paper_chunks`, BM25; corpus 9376→12810 chunks, 278 refs). This manifest stages the
*contract-graph* rows for **PI ratification into Notion** (the SoT). **New gates/contracts are
PI-authored — nothing here is auto-created.** Full mechanistic detail + per-engine directions:
`ffn_sim/docs/v2_audit/ENGINE_ROADMAP_KIM_MIYAZAKI_2026-07-10.md`. Per-paper dossiers + 3-lens analyses:
session scratchpad `dossiers/`, `lenses/`, `synth/digest_{dcm,ff,tag}.md`.

**Rules honoured.** Experimental values = overlay/cross-check, never fit targets. Bands are contracts
(not loosened). No magic numbers. Simulation-origin constants are flagged, not registered. Before citing
a source in a deliverable, confirm its `source_audit.verdict` is OK.

---

## Part A — SourceEvidence rows (bibliographic; 77 papers)
Appended below by `build_se_table.py` as a table: `pNN | citation_key | lab | type | year | journal |
DOI | title`. These are the least-gated rows (pure bibliography); recommend batch-adding to
SourceEvidence, then linking to the KnowledgeClaims in Part B. DOIs should be run through
`verify_sources.py` before deliverable citation.

---

## Part B — KnowledgeClaim candidates (new KB-x.y, PI to assign ids)

Grouped by mechanism; each lists the anchoring paper(s). These extend the existing KB (162→ rows) with
the mechanistic content the corpus adds. `→ engine` marks the primary consumer.

**Actin turnover & fragmentation (→ FF):**
- KC-cand-1: Cortical stress relaxation is dominated by *filament turnover* (treadmilling; link dies when
  host segment depolymerizes), with η_eff ∝ 1/turnover; force-dependent crosslinker unbinding is
  secondary below ~10 Pa. [p03 BPJ2014, p04 BPJ2018, p24 Cytoskeleton2019, p60 ncomms10323]
- KC-cand-2: A second, orthogonal relaxation channel is buckling-curvature-gated stochastic F-actin
  severing k_sev=k0·exp(λθ); both channels are required together for physiological cortical pulses.
  [p42 ACSMacroLett2016, p35 SoftMatter2017, p26 Cytoskeleton2024, p48 PLOSCB2025, p27 Cytoskeleton2025]
- KC-cand-3: A single dimensionless control Π = turnover-rate / motor-walk-rate collapses cortex
  stress-vs-turnover to a master curve and bounds sustained prestress from above. [p60, p04]

**Crosslinker mechanics (→ FF):**
- KC-cand-4: Load-bearing crosslinkers are *catch-slip* (α-actinin-4 lifetime peaks ~4 pN); catch bonds
  "dissociate on demand", homogenize per-bond force, and make weaker networks stronger + more
  deformable — requires crosslinker mobility (unbind→hop≤5µm→rebind) and N≳10. [p74 NatMater2022, p10
  ActaBiomater2025]
- KC-cand-5: Cross-link *unbinding* (not Ig-domain unfolding) governs cortex rheology at physiological
  strain rate (~0.1 s⁻¹); unfolding dominates only above ~1 s⁻¹. Filamin-A unbind 70±23 pN, Ig-unfold
  57±19 pN, sawtooth 28 nm. [p02 BPJ2011, p67 CMBE2009]
- KC-cand-6: Sustained cortical tension is set by the crosslinker off-rate ratio (clutch), not motor
  count; filament buckling is REQUIRED for contractile tension (100×κ_b → tension→0). Max-tension ∝
  R_M·R_ACP, sustainability ∝ R_ACP/R_M. [p63 BMMB2015, p69 CompPartMech2015, p44 PLOSCB2017, p15
  biorxiv2025]

**Myosin motor (→ FF):**
- KC-cand-7: Non-muscle myosin II is better modelled as an Erdmann-Schwarz parallel-cluster catch-bond
  motor (walk & unbind both fall under load) on an explicit bipolar minifilament (bare zone + arm count
  + spacing); network force scales as √N_M (isotropic slab) vs N_M (bundle). [p41 eLife2025, p61
  ncomms12615, p05 BPJ2019, p36 SoftMatter2020, p15 biorxiv2025]

**Formin / nucleation (→ FF):**
- KC-cand-8: Formin mDia1 barbed-end elongation is a Bell-type force sensor v=k_on·C_A(1−p_c)−k_off·p_c,
  p_c=p_c⁰·exp(−fd/kT), d=5.5 nm — tension *accelerates* elongation; discrete ~2.7 nm stepping.
  [p57 BPJ2017, p62 NanoLett2018]
- KC-cand-9: A finite shared G-actin pool makes filament length *emergent* (self-limiting), replacing an
  imposed length distribution; polymerization growth-pressure can buckle stalled barbed ends (motor-free
  templating). [p21 AdvFunctMater2019, p53 CellRep2022]

**Cell division / cortex (→ DCM):**
- KC-cand-10: Cytokinetic-ring constriction under enclosed-volume conservation produces large axial
  elongation from small volume loss (~5% loss → ~60% elongation); the spindle is not required for
  elongation; division is neighbour-transmitted. [p46 JCB2021, p22 AdvSci2021, p19 biorxiv2026]
- KC-cand-11: In cell-sized confinement (R ≲ ℓ_p) an equatorial contractile ring self-organizes;
  constriction rate dP/dt=−k_c(P−P*), k_c=0.39 µm/min per µm, P*≈mean filament length; effective myosin
  concentration ∝ (duty ratio)². [p59 NatCellBiol2015, p28 review]
- KC-cand-12: Cortical tension is an emergent actomyosin-density field with in-plane flow (σ_act=(ζμ)₀ρ/
  (ρ₀+ρ); reaction-advection-diffusion transport); the contractile ring emerges as a cortical wave.
  [p30 PhysRevResearch2023, p61 ncomms12615, p60 ncomms10323]
- KC-cand-13: Membrane–cortex coupling is a layer of breakable Bell linkers; blebs select between
  Detachment and Rupture via dimensionless coupling R_C and connectivity R_X; detachment tension oracle
  γ_C^{D*}∝R_C. [p14 biorxiv2025, p55 BPJ2015]

**Aggregate compaction rate / junctions (→ DCM):**
- KC-cand-14: Loose→compact aggregate compaction rate is set by a percolation-gated maturation kernel
  τ_p=2^N·τ (exponential-in-size), NOT by cadherin single-bond kinetics (which max at ~6 s and do not
  gate 24–48 h compaction); junction/cortex plasticity = Bell off-rate + rebind → permanent deformation.
  [p70 ncomms2020, p47 PLOSCB2022, p07 MatrixBiol2019, p03 BPJ2014]  **Resolves the open item in
  project-dcm-cadherin-cluster-redesign ("(b) T1 escalation needed").**

**ECM / collagen (→ FF + DCM):**
- KC-cand-15: Collagen ECM viscoelasticity/plasticity emerges from transient Bell-slip crosslink
  turnover + fibril rearrangement; stress-relaxation rate is set by connectivity/heterogeneity at
  invariant stiffness; long-range ~1/r stress transmission is anisotropic (radial n=1.2 vs tangential
  n=2.2) and is a *fibrous-network* property. [p39 SoftMatter2021, p07 MatrixBiol2019, p75 Nature2024,
  p13 biorxiv2025, p22 AdvSci2021, p19 biorxiv2026]

**Microtubule/kinesin (→ FF future compartment):**
- KC-cand-16: A first-class MT+kinesin+MAP compartment: three-state plus-end dynamic instability +
  angle-dependent zippering/catastrophe + stress→plus-end catch-bond feedback; extensile kinesin bundles
  buckle under confinement. [p68 BMCPlantBiol2023, p77 PNAS2017, p08 CurrBiol2025, p40 ACSNano2022]

---

## Part C — Parameter candidates (register with provenance; flag supersession)

| name | value | unit | condition / source | note |
|---|---|---|---|---|
| formin mDia1 k_on (ATP/ADP) | 21 / 3.6 | µM⁻¹s⁻¹ | 26°C, single-mol [p57] | FF formin |
| formin mDia1 k_off (ATP/ADP) | 1.3 / 3.9 | s⁻¹ | [p57] | FF formin |
| formin p_c⁰ (ATP/ADP) | 0.29 / 0.61 | — | closed-state prob [p57] | FF formin |
| formin working distance d | 5.5 | nm | FH2 half [p57] | Bell length |
| formin barbed-end speed | 1.1–1.3 | µm/s (~450 mon/s) | [p53] | vs FF v0=0.623 (reconcile) |
| formin step size | 2.7 | nm | [p62] | discrete stepping |
| α-actinin-4 catch peak force | ~4 | pN | single-mol [p74] | FF catch crosslinker |
| filamin-A unbinding force | 70 ± 23 | pN | OTFS [p67] | FF crosslinker F_b anchor |
| filamin-A Ig-unfold force | 57 ± 19 | pN | [p67] | secondary channel |
| filamin-A k0_off | 0.115 | s⁻¹ | Ferrer 2008 [p02,p03,p53,p69] | **CONFLICT: distinct isoform vs KB α-actinin 0.066** |
| NMII stall force | 5.0–5.7 | pN/head | PCM [p05,p14,p26,p69] | **CONFLICT vs PARAM-F_stall_motor 2 pN "per motor"** |
| NMII unloaded velocity | 140 | nm/s | PCM [p05,p14] | **vs PARAM-v_unloaded 100** |
| NMII minifilament length | 0.64 ± 0.30 | µm | TIRF n=550 [p14] | FF minifilament (SMM 1.32 µm [p73]) |
| NMII heads per minifilament | 64 (8×8) / 32 (8×4) | — | [p26,p41,p53] | FF discretization |
| myosin duty ratio | 0.04–0.05 (II) / 0.5–0.7 (V) | — | [p59] | duty² law |
| F-actin severing threshold | 300–500 | pN | Tsuda 1996 [p14,p26] | FF severing gate |
| actin persistence length | 9–18 | µm | [p14,p59] | distinct from collagen ℓ_p 17 |
| cortex tension plateau (reconstituted) | ~0.2 (200 pN/µm) | mN/m | ∝R_C [p14] | corroborates interphase γ~0.15 |
| ring constriction rate | 0.39 | µm/min per µm | [p59] | division |
| critical contractable perimeter P* | ~4.5 (≈filament length) | µm | [p59] | division |
| division volume loss → elongation | ~5% → ~60% | — | MDCK n=22 [p46] | division gate |
| ring / spindle force | ~100 / ~100 | pN/side | [p22] | division |
| collagen G'(conc) | 10 / 100 / 400 | Pa @ 1/3/5 mg/mL | [p22] | extends VG-U1-G0 |
| ECM crosslink k0_u, x_u | 3e-6, 1e-10 | s⁻¹, m (F_b≈41 pN) | [p07,p39] | ECM (not cortex) |
| percolation maturation kernel | τ_p = 2^N·τ | — | [p70] | DCM compaction-rate |
| cell-sized confinement threshold | R ≈ ℓ_p (~15–20) | µm | [p59,p28,p70] | confinement scenario |
| MT persistence length | 1–2 | mm | [p77,p68] | FF MT compartment |
| MT bundle elongation rate | 3.58 ± 1.61 | µm/min | [p77] | FF MT compartment |

Per-segment cylindrical drag `ζ=3πμr_c(3+2r₀/r_c)/5` (p03,p13,p14,p19,p22,p44) — register as a
*derived form* (not a scalar), candidate fix for the FF crawl-speed grid-drag problem.

---

## Part D — ValidationGate candidates (PI-authored; new lanes)

| proposed vg_id | observable | band / value | engine | source |
|---|---|---|---|---|
| VG-DCM-division-elongation | axial elongation vs enclosed-volume loss | ~60% elong @ ~5% vol loss; ≥40% loss to abolish | DCM | p46 |
| VG-DCM-ring-constriction | ring constriction rate vs perimeter | dP/dt=−0.39 µm/min·(P−P*)/µm, P*≈filament length | DCM | p59 |
| VG-DCM-division-force | ring/spindle force + elongation split | ~100 pN/side; 20% (ring) / 40% (ring+spindle) | DCM | p22 |
| VG-FF-formin-fv | mDia1 elongation vs load | Bell law d=5.5 nm; ADP 3.4→10.9 subs/s at 0.13→4 pN | FF / oracle | p57 |
| VG-FF-crosslinker-catch | α-actinin-4 lifetime vs force | catch peak ~4 pN; K255E monotone slip | FF | p74 |
| VG-FF-network-rupture | rupture stress & strain vs bond kinetics | catch 24.5 Pa/221% vs slip 6.5–8.1 Pa/63–129% | FF | p74 |
| VG-FF-contraction-percolation | εmax vs myosin density | percolation ρc≈0.56 µm⁻²; telescopic law | FF | p61 |
| VG-FF-severing-threshold | tensile F-actin severing | F_sev 300–500 pN | FF | p14,p26 |
| VG-U1-ECM-viscoplasticity | creep-recovery permanent strain @ invariant modulus | 80%→10% as covalent crosslink↑ | FF/DCM ECM | p07,p75 |
| VG-U1-ECM-stress-relaxation | τ½ vs fibril/bundle length, angle, connectivity | faster relaxation at shorter/less-connected, stiffness invariant | FF/DCM ECM | p75,p39 |
| VG-U1-force-propagation-aniso | radial vs tangential decay exponent | n=1.2 (radial) vs 2.2 (tangential) — TIGHTENS existing 1/r gate | ECM | p13 |
| VG-DCM-confined-motility | confined migration speed vs height | non-monotonic V(h), optimum ~60 µm; V_drop closed form | DCM / oracle | p76 |
| VG-DCM-bleb-detachment | bleb onset tension vs coupling | γ_C^{D*}∝R_C (Bell + Young-Laplace) | DCM / oracle | p14,p55 |

**Oracle candidates** (analytic cross-checks): formin Bell elongation law (p57); confined-motility V_drop
(p76); bleb detachment-tension γ_C^{D*} (p14); duty² effective-concentration + constant-volume
constriction invariants (p59); ⟨ω⟩∝1/R confined rotational flow (p77).

---

## Part E — Conflicts / decisions to surface to PI (do NOT silently reconcile)
1. Crosslinker k0=0.115 s⁻¹ (filamin-A) vs KB 0.066 s⁻¹ (α-actinin) — same Ferrer 2008 paper, different
   isoform → register as **distinct species**, do not overwrite the α-actinin anchor.
2. Myosin stall 5–5.7 pN/head vs `PARAM-F_stall_motor` default 2 pN "per motor" (per-head vs per-motor);
   NMII v0 140 nm/s vs `PARAM-v_unloaded` 100 → PI to decide per-head vs per-motor convention.
3. `PARAM-F_stall_motor` default 2 pN and `PARAM-v_unloaded` 100 nm/s are currently **uncalibrated
   (`calibrated=No`)** — this corpus supplies sourced values; recommend re-anchoring under PI sign-off.

## Part F — SI-fetch backlog (fetch before registering these force constants)
- p74 Mulla — Suppl. Table 1 (1D Bell k_catch0/k_slip0/f_catch/f_slip) + Table 2 (3D network κ's).
- p75 Fan — Suppl. Table 5 (fibril/crosslinker/bundler κ, Δt, k_off0, x_b); also github.com/ktyman2/liverCancer.
- p22 Nam — Table S1 (division-in-collagen sim params).
- p26 Matsuda — Table S1 (bipolar minifilament geometry).
- p05 Kim — SI numeric tables (motor mobility ζ_M, PCM rates).


### Part A table (auto-generated)

| pNN | lab | type | DOI | title (truncated) |
|---|---|---|---|---|
| p01 | Miyazaki | research | 10.1063/1.3516587 | RESEARCH ARTICLE /  FEBRUARY 28 2011 |
| p02 | Kim | research | 10.1016/j.bpj.2011.08.033 | Dynamic Role of Cross-Linking Proteins in Actin Rheology |
| p03 | Kim | research | 10.1016/j.bpj.2013.12.031 | Determinants of Fluidlike Behavior and Effective Viscosity in Cr |
| p04 | Kim | research | 10.1016/j.bpj.2018.10.008 | Balance between Force Generation and Relaxation |
| p05 | Kim | research | 10.1016/j.bpj.2019.04.018 | Mobility of Molecular Motors Regulates Contractile |
| p06 | Kim | review | 10.1016/j.ydbio.2018.12.004 | Contents lists available at ScienceDirect |
| p07 | Kim | research | 10.1016/j.matbio.2019.05.006 | Covalent cross-linking of basement |
| p08 | Miyazaki | research | 10.1016/j.cub.2025.11.073 | Chimeras of kinesin-6 and kinesin-14 reveal head- |
| p09 | Kim | research | 10.1016/j.actbio.2025.05.069 | Reciprocal folding dynamics in cellular networks at the stroma-b |
| p10 | Kim | research | 10.1016/j.actbio.2025.06.004 | Fine-tuning of material properties by catch bonds |
| p11 | Miyazaki | research | 10.1063/1.3574396 | RESEARCH ARTICLE /  APRIL 07 2011 |
| p12 | Miyazaki | review | 10.2142/biophysico.bppb-v19.0031 | Honmachi, Sakyo-ku, Kyoto 606-8501, Japan. ORCID iD: https://orc |
| p13 | Kim | research | 10.1101/2025.01.31.635980 | Swirling motion of breast cancer cells radially aligns collagen  |
| p14 | JOINT | research | 10.1101/2025.05.18.654456 | Reconstitution of actomyosin networks in cell-sized |
| p15 | Kim | research | 10.1101/2025.07.09.663940 | Dissecting Molecular Origins of the Mechano-Adaptive |
| p16 | Kim | research | 10.1101/2025.09.25.678704 | Computational Modeling of Motor-Driven |
| p17 | Kim | research | 10.64898/2026.02.08.704688 | Kinetic Control of Out-Of-Equilibrium |
| p18 | Kim | research | 10.64898/2026.05.19.726229 | Durotactic Migration Driven by Anisotropic Matrix |
| p19 | Kim | research | 10.64898/2026.06.16.732593 | Mechanical Checkpoint for Cell Division in |
| p20 | Miyazaki | experimental | 10.1119/5.0059810 | PAPERS /  SEPTEMBER 01 2022 |
| p21 | Kim | research | 10.1002/adfm.201905243 | Dr. V. Yadav, Dr. A. P. Tabatabai, Prof. M. P. Murrell |
| p22 | Kim | research | 10.1002/advs.202000403 | Cellular Pushing Forces during Mitosis Drive Mitotic |
| p23 | Miyazaki | method | 10.21769/bioprotoc.5656 | Bio-protocol 16(8): e5656. DOI: 10.21769/BioProtoc.5656 |
| p24 | Kim | research | 10.1002/cm.21582 | R E S E A R C H A R T I C L E |
| p25 | Kim | research | 10.1002/cm.21808 | R E S E A R C H A R T I C L E |
| p26 | Kim | research | 10.1002/cm.21848 | R E S E A R C H A R T I C L E |
| p27 | Kim | research | 10.1002/cm.70014 | Cytoskeleton, 2026; 83:407–426 |
| p28 | Miyazaki | review | 10.1080/15421406.2017.1289445 | Molecular Crystals and Liquid Crystals |
| p29 | Kim | research | 10.1103/physrevlett.133.158402 | Mechanics and Morphology of Proliferating Cell Collectives with  |
| p30 | Miyazaki | research | 10.1103/physrevresearch.5.013208 | PHYSICAL REVIEW RESEARCH 5, 013208 (2023) |
| p31 | Kim | research | 10.1103/physrevresearch.7.013312 | PHYSICAL REVIEW RESEARCH 7, 013312 (2025) |
| p32 | Kim | review | 10.1115/1.4046863 | Weldon School of Biomedical Engineering, |
| p33 | Miyazaki | experimental | 10.1039/c3ra41112e | Accurate polarity control and parallel alignment of actin |
| p34 | Kim | review | 10.1039/c5ib00043b | Integr. Biol., 2015, 7, 1093--1108 / 1093 |
| p35 | Kim | research | 10.1039/c6sm02703b | Soft Matter, 2017, 13, 3213--3220 / 3213 |
| p36 | Kim | research | 10.1039/c9sm02082a | 1548 / Soft Matter, 2020, 16, 1548--1559 |
| p37 | Kim | research | 10.1091/mbc.e23-09-0364 | Molecular Biology of the Cell • 35:ar47, 1–12, April 1, 2024 |
| p38 | Miyazaki | research | 10.1021/acs.nanolett.3c02742 | Controlling Physical and Biochemical Parameters of Actin |
| p39 | Kim | research | 10.1039/d0sm01911a | 10274 /  Soft Matter, 2021, 17, 10274–10285 |
| p40 | Kim | research | 10.1021/acsnano.2c05593 | Durability of Aligned Microtubules Dependent |
| p41 | Kim | research | 10.7554/elife.105236 | Ding, Chou et al. eLife 2025;14:RP105236. DOI: https://doi.org/1 |
| p42 | Kim | research | 10.1021/acsmacrolett.6b00232 | F‑Actin Fragmentation Induces Distinct Mechanisms of Stress |
| p43 | Kim | research | 10.1371/journal.pcbi.1000439 | Computational Analysis of Viscoelastic Properties of |
| p44 | Kim | research | 10.1371/journal.pcbi.1005277 | Morphological Transformation and Force |
| p45 | Kim | research | 10.1167/iovs.66.2.65 | Characterization, Enrichment, and Computational |
| p46 | Kim | research | 10.1083/jcb.202011106 | The nature of cell division forces in epithelial |
| p47 | Kim | research | 10.1371/journal.pcbi.1010257 | Role of actin filaments and cis binding in |
| p48 | Kim | research | 10.1371/journal.pcbi.1013572 | PLOS Computational Biology / https://doi.org/10.1371/journal.pcb |
| p49 | Miyazaki | research | 10.1021/acs.biochem.5c00495 | Kinetic Scheme of Myosin Phosphorylation by ZIP Kinase |
| p50 | : | ? | 10.1091/mbc.e22-06-0243 | Molecular Biology of the Cell • 34:ar67, 1–10, June 1, 2023 |
| p51 | Kim | research | 10.1021/acsbiomaterials.8b01365 | Mechanical Model for Durotactic Cell Migration |
| p52 | Miyazaki | research | 10.1021/acs.nanolett.1c02955 | Microscopic Temperature Control Reveals Cooperative Regulation |
| p53 | Kim | research | 10.1016/j.celrep.2022.110868 | Rapid assembly of a polar network architecture |
| p54 | Miyazaki | experimental | 10.1016/j.bpj.2025.12.013 | on the number and size of liposomes formed |
| p55 | Miyazaki | experimental | 10.1016/j.bpj.2015.06.016 | Directional Bleb Formation in Spherical Cells under Temperature  |
| p56 | Miyazaki | method | 10.1016/j.bpj.2014.05.039 | Quantitative Analysis of the Lamellarity of Giant Liposomes Prep |
| p57 | Miyazaki | experimental | 10.1016/j.bpj.2017.06.012 | Biphasic Effect of Proﬁlin Impacts the Formin mDia1 |
| p58 | Miyazaki | research | 10.1021/acs.nanolett.5c02558 | Myosin-Driven Advection and Actin Reorganization Control the |
| p59 | Miyazaki | research | 10.1038/ncb3142 | Cell-sized spherical conﬁnement induces the |
| p60 | Kim | research | 10.1038/ncomms10323 | Interplay of active processes modulates tension |
| p61 | Kim | research | 10.1038/ncomms12615 | Disordered actomyosin networks are sufﬁcient to |
| p62 | Miyazaki | research | 10.1021/acs.nanolett.8b03277 | Processive Nanostepping of Formin mDia1 Loosely Coupled with |
| p63 | Kim | research | 10.1007/s10237-014-0608-2 | Biomech Model Mechanobiol (2015) 14:345–355 |
| p64 | Kim | research | 10.1007/s10237-015-0660-6 | Biomech Model Mechanobiol (2015) 14:1143–1155 |
| p65 | Kim | research | 10.1007/s11340-007-9091-3 | Computational Analysis of a Cross-linked |
| p66 | Kim | research | 10.1007/s11538-019-00585-1 | Bulletin of Mathematical Biology (2019) 81:3301–3321 |
| p67 | Kim | research | 10.1007/s12195-009-0048-8 | Cytoskeletal Deformation at High Strains and the Role of Cross-l |
| p68 | Kim | research | 10.1186/s12870-023-04252-5 | sharing, adaptation, distribution and reproduction in any medium |
| p69 | Kim | research | 10.1007/s40571-015-0052-9 | Comp. Part. Mech. (2015) 2:317–327 |
| p70 | Miyazaki | research | 10.1038/s41467-020-16677-9 | Tug-of-war between actomyosin-driven |
| p71 | Kim | research | 10.1038/s41467-022-34715-6 | F-actin architecture determines constraints |
| p72 | Miyazaki | research | 10.1038/s41467-025-62076-3 | Protein design of two-component tubular |
| p73 | Miyazaki | research | 10.1038/s41467-025-62653-6 | Optogenetic actin network assembly on lipid |
| p74 | Kim | research | 10.1038/s41563-022-01288-0 | 1Living Matter Department, AMOLF, Amsterdam, The Netherlands. 2I |
| p75 | Kim | research | 10.1038/s41586-023-06991-9 | Nature  /  Vol 626  /  15 February 2024  /  635 |
| p76 | Miyazaki | research | 10.1073/pnas | BIOPHYSICS AND COMPUTATIONAL BIOLOGY |
| p77 | Miyazaki | research | 10.1073/pnas.1616001114 | Spatial confinement of active microtubule networks |

_77 papers. citation_keys in `paper_refs` (BM25). Labs: Kim=51, Miyazaki=24, Other=0._
