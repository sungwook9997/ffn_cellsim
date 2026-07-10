# Engine roadmap from the Kim Taeyoon + Miyazaki corpus (77 papers) — 2026-07-10

**Mandate (PI):** absorb *every page* of *every* Kim Taeyoon– and Makito Miyazaki–lab paper (plus
close-neighbour active-matter / division / durotaxis / reconstitution work) placed in
`ffn_sim/references/2026_07_10`, and turn them into concrete improvement directions for our engines
(DCM, FF), a judgment on whether a new engine is warranted, and the validation data TAG still needs —
all at *our* single-cell fine-grained scale.

**How this was produced.** 80 PDFs → 77 unique (3 byte-dups dropped). Each paper was read page-by-page
by a dedicated subagent that wrote a complete faithful **dossier** (equations, force laws, parameter
tables, quantitative results, per-page coverage log) under `scratchpad/dossiers/pNN.md`; then **three
independent lens subagents** (DCM engine / FF engine / TAG-validation-gap) mined each dossier, writing
`scratchpad/lenses/pNN.{dcm,ff,tag}.md`. 308 agents total (the 15 tail lenses that hit the account
session-limit were reconstructed by the Lead from their full dossiers). Structured summaries per paper
live in `scratchpad/synth/digest_{dcm,ff,tag}.md`. All 77 PDFs were also ingested into the TAG content
layer (`paper_chunks`, BM25). Nothing was skipped.

**One-line verdict.** The corpus is, to a first approximation, **the published source code of the FF
engine**: Taeyoon Kim's agent-based Brownian-dynamics actomyosin/collagen models are the *exact physics
family* FF implements (cylindrical filament segments + hinged crosslinkers + bipolar myosin + Bell-Evans
kinetics). It hands FF the mechanisms it is currently missing — **turnover, buckling-gated severing,
catch-slip crosslinkers, PCM myosin, force-dependent formin, principled per-segment drag** — most of
which directly attack our three standing open problems (γ-floor, aggregate-compaction rate/T1,
crawl-speed grid-drag). Miyazaki's in-vitro reconstitution + Kim's surface-mesh division paper hand
**DCM** its single biggest missing module (cytokinetic division) plus emergent-γ, bleb, and a
mechanistic compaction-rate clock. **No new engine is required** — the corpus fills the existing two;
the only *new compartment* worth staging is a first-class **microtubule/kinesin module** inside FF.

Paper key: `pNN` ↔ `scratchpad/worklist.json`. Scores below are the lens relevance/priority (0–5).

---

## 1. FF engine roadmap (Filament-FEM, Cytosim physics) — the corpus's centre of mass

38 of 77 papers score ≥3.5 for FF. They converge on a small number of high-leverage mechanisms. Ordered
by impact × readiness.

### T1 — Filament turnover + buckling-gated severing = the missing stress-relaxation channel  ★★★ (top priority)
FF filaments are currently unbreakable, non-treadmilling beams. Every Kim rheology/contraction paper
says this is *the* gap. Two orthogonal, both-required channels:
- **Treadmilling turnover** (barbed-polymerize = pointed-depolymerize), releasing *all* bound
  crosslinkers/motors on a removed segment — p03 (BPJ 2014, 10.1016/j.bpj.2013.12.031; single global
  rate `n`, Eq.5 `k+_B=6n…`, η_eff ∝ 1/turnover), p04 (BPJ 2018), p24 (Jung, Cytoskeleton 2019,
  10.1002/cm.21582; turnover-rescue threshold k_t,A≈20 s⁻¹), p60 (ncomms10323 2016; nondimensional
  control Π=k_turnover/v_motor collapses stress to a master curve — the single cortex-stability knob).
- **Buckling-curvature-gated stochastic severing** `k_sev = k0_sev·exp(λ_sev·θ)` — p42 (ACS Macro Lett
  2016, 10.1021/acsmacrolett.6b00232), p35 (Soft Matter 2017), p26 (Matsuda, Cytoskeleton 2024,
  10.1002/cm.21848; tension-severing F_frag≈500 pN), p48 (PLOS CB 2025, angle-severing unjams
  myosin-compacted bundles), p27 (Cytoskeleton 2025, full turnover kit). Severed ends feed formin/Arp2/3
  regrowth → closed turnover balance.
**Why it matters here:** this is the direct mechanistic route into the **γ-floor** — Kim shows sustained
cortical tension is *turnover*- and *crosslinker-clutch*-limited, not motor-count-limited (see T3, T4).

### T2 — Catch-slip crosslinker off-rate + crosslinker mobility (unbind→hop→rebind)  ★★★
FF uses a slip-only Bell-Evans α-actinin off-rate. Upgrade to the two-pathway catch/slip form and add
mobility:
- p10 (Acta Biomater 2025, 10.1016/j.actbio.2025.06.004) — catch-slip off-rate `k=k_s0·e^{+F/λ_s} +
  k_c0·e^{−F/λ_c}`; crosslinker **turnover** (unbind→hop ≤5 µm→rebind) is a capability *absent in
  AFINES/MEDYAN/Cytosim* and is *required* for the catch advantage; force-dependent CL redistribution
  amplifies myosin contraction ~3–4×. Store primitives {k_s0,k_c0,λ_s,λ_c}; slip-only recovered at
  k_c0=0 (back-compat).
- p74 (Mulla/Kim, Nat Mater 2022, 10.1038/s41563-022-01288-0) — single-molecule α-actinin-4 catch peak
  ~4 pN; "dissociation on demand" homogenizes per-bond force and postpones rupture; catch needs
  mobility + N≳10. Validation: rupture stress 24.5 vs 6.5–8.1 Pa, strain 221 vs 63–129%.
- Single-molecule crosslinker constants to anchor F_b: p67 (CMBE 2009; filamin-A unbind 70±23 pN,
  Ig-unfold 57±19 pN, 28 nm sawtooth, angle-rupture >45°), p02 (Kim-Kamm, BPJ 2011; WLC arm law,
  unbinding — not unfolding — governs rheology at physiological strain rate; filamin k0=0.115 s⁻¹).
  **Conflict to flag PI:** k0=0.115 s⁻¹ is *filamin-A* (Ferrer 2008); FF's α-actinin anchor is 0.066 s⁻¹
  from the *same* paper, different isoform — register filamin as a **distinct species**, do not overwrite.

### T3 — PCM (parallel-cluster-model) catch-bond myosin, finer than Hill  ★★★
Replace/augment the imposed Hill force-velocity with the Erdmann-Schwarz per-head Markov PCM where both
walking and unbinding rates fall under load (a genuine catch-bond motor), on an explicit bipolar
minifilament geometry:
- p41 (Ding/Chou, eLife 2025, 10.7554/eLife.105236) — explicit bare-zone L_bz + arm-count N_a +
  arm-spacing; **bundle-vs-network force scaling F_tot ∝ N_M (bundle) vs √N_M (isotropic slab)** — a
  direct lever on the cortical-tension γ-floor; emergent ACP-partitioned contractile units (a free
  self-test FF should already pass).
- p15 (biorxiv 2025, single stress fiber; per-head v0=56 nm/s, F_s=24 pN, 30 pN slip cutoff; buckling =
  compression rate > unloaded walking speed), p05 (BPJ 2019; motor **mobility** drag ζ_M as explicit
  DOF, one-vs-two-arm connectivity, single ATP knob k₂₀ with stall∝1/k₂₀), p36 (Soft Matter 2020),
  p61 (ncomms12615 2016; hinged 3-seg×42 nm backbone, 4 arms, 7 nm lattice — pins FF's Stam-Hocky
  discretization). Keep FF's Stam-Hocky+Hill as the sanctioned cross-check oracle; PCM is the runtime
  refinement. Anchors: NMII v0≈140 nm/s, stall 5.7 pN/arm (p14, p69).

### T4 — Crosslinker-as-clutch: sustained tension is set by off-rate ratio, not motor force  ★★★ (γ-floor)
p63 (BMMB 2015, 10.1007/s10237-014-0608-2), p69 (Comp Part Mech 2015), p44 (PLOS CB 2017,
10.1371/journal.pcbi.1005277): a transient crosslinker clutch (its Bell-slip off-rate), not the motor
count, decides whether motor tension is *sustained* or dissipated; **filament buckling is REQUIRED for
contractile tension** (100×κ_b → tension collapses to ~0), validating FF's bending-beam as necessary.
Max-tension ∝ R_M·R_ACP, sustainability ∝ R_ACP/R_M — the mechanistic sustainability law for the open
γ/contraction problem. p61 adds the falsification guardrail: buckling is a *rate* modulator, not the
mechanism (ℓ_p 20 µm→2.6 mm slows but preserves contraction).

### T5 — Force-dependent formin elongation + discrete nanostepping  ★★
FF's barbed-end growth is a force-blind constant-velocity ratchet. Upgrade:
- p57 (Miyazaki, BPJ 2017, 10.1016/j.bpj.2017.06.012) — mDia1 is a **Bell-type force sensor**,
  `v = k_on·C_A·(1−p_c) − k_off·p_c`, `p_c = p_c⁰·e^{−fd/kT}`, d=5.5 nm; tension *accelerates* (mDia1 ≠
  Bni1p); nucleotide-resolved kon/koff/p_c⁰; profilin biphasic force switch (−1.8→+40.2 subs/s at ~4 pN).
  **The closed form is an analytic-oracle candidate.**
- p62 (Miyazaki, Nano Lett 2018) — discrete ~2.7 nm Poisson barbed-end stepping (replace lumped drift).
- p21 (Yadav/Murrell, AdvFunctMater 2019, 10.1002/adfm.201905243) — finite shared G-actin pool →
  **emergent** filament-length distribution (kills the imposed exponential, kills a magic number);
  polymerization growth-pressure buckling (motor-free templating). p53 (Cell Rep 2022): barbed-end pool
  depletion, emergent length law L=v/(1/τ_actin+1/τ_formin)≈6 µm, formin v 1.1–1.3 µm/s.

### T6 — Principled per-segment cylindrical drag = candidate fix for the crawl-speed grid-drag problem  ★★
The open item from [[project-ff-crawl-mechanism]] (native crawl 500× slow, Σγ ∝ Nc). Kim uses a
geometry-derived per-segment Stokes drag `ζ = 3πμ·r_c·(3+2r₀/r_c)/5` (p03 Eq.4, p13, p22 Eq.4, p14 S3
Underhill-Doyle, p19, p44) at physiological μ. This is a *physical* reframing of FF's count-scaling
`--bulk-drag`; adopting it (with the anisotropic ζ⊥/ζ∥≈1.64 of p43) is the principled route to a stable,
grid-invariant drag instead of the destabilizing bulk-drag hack.

### T7 — Crosslinker angular mechanics + torsion/helical registry (emergent bundle-vs-mesh)  ★
FF lumps the crosslinker into one scalar stiffness. Kim resolves three modes (arm stretch + two bending
angles): p43 (PLOS CB 2009, 10.1371/journal.pcbi.1000439; θ1,eq = π bundle vs 0.417π orthogonal mesh),
p64 (BMMB 2015), p65 (Exp Mech 2009; filament **torsion** κ_t + helical binding-site registry make the
cross-link angle emergent). This is how FF should build filopodial/stress-fiber *bundles* mechanistically
(cf. p75's two-angle bundler for collagen).

### T8 — ECM/collagen side of FF (shares the exact same BD physics)  ★★
Kim's collagen models are FF with fibrils↔filaments, bundlers↔fascin: p13 (biorxiv 2025,
10.1101/2025.01.31.635980; negative normal stress, annular swirl→radial alignment, full oracle param
set), p39 (Soft Matter 2021; transient Bell ECM crosslinks → emergent viscoelasticity/plasticity),
p07 (Matrix Biol 2019; covalent = k0→0 limit, creep-recovery plasticity), p19/p22 (division-in-collagen),
p75 (Nature 2024; two-angle bundler, stiffness-invariant relaxation). Directly extends the existing
DCM⊗ECM 1/r work into FF's fine-grained collagen substrate.

### T9 — New FF compartment: microtubule + kinesin (extensile)  ★ (stage, don't rush)
p68 (BMC Plant Biol 2023; MT three-state dynamic instability rate set + angle-dependent
zippering/catastrophe + stress→plus-end Bell-catch feedback), p77 (Suzuki/Miyazaki, PNAS 2017; extensile
kinesin bundles, MT ℓ_p 1–2 mm, elongation 3.58 µm/min, buckling under confinement). FF already exposes
`--microtubules`; this is the mechanistic content to make it first-class (see §3).

---

## 2. DCM engine roadmap (Deformable Cell Model, SimuCell3D physics)

Fewer but very high-leverage papers. The headline is **cell division**, which is DCM's single biggest
missing module — and the corpus contains a *published Taeyoon-Kim triangulated-surface-mesh model that is
essentially DCM* plus Miyazaki's in-vitro ring reconstitution.

### D1 — Cytokinetic division module = DCM's #1 gap, now drop-in  ★★★ (top priority)
- p46 (Kim, JCB 2021, 10.1083/jcb.202011106) — **a T. Kim triangulated-surface-mesh DCM**: prescribed
  inward constriction of a 1-µm equatorial node band + **enclosed-volume conservation** yields ~60%
  axial elongation from only ~5% volume loss (≥40% loss to abolish). No ring line-tension magic number;
  node-area-proportional drag ζ_i=λA_i (grid-invariant); membrane-reservoir soft area law (κ_a
  10⁻⁴…10⁻⁶ N/m, +22% area); **spindle not required** for elongation; division is neighbour-transmitted
  (build it in the N≈400 aggregate). This is a near drop-in DCM module, equations + validation oracle.
- p59 (Miyazaki, Nat Cell Biol 2015, 10.1038/ncb3142) — *emergent* equatorial furrow (energy-min great
  circle, R_ring/R_droplet=0.96), constriction law `dP/dt = −k_c(P−P*)`, k_c=0.39 µm/min per µm,
  **P*≈mean filament length 5.4 µm**, constant-volume constriction, duty-ratio² effective-tension,
  motor-escalation life cycle to complete R_ring→0. Confinement criterion R_cell ≲ ℓ_p.
- p22 (Nam, AdvSci 2021, 10.1002/advs.202000403) — ring + spindle forces (100 pN/side), elongation 20%
  (ring-only) / 40% (ring+spindle), volume-conserved, divide against real collagen G'=10/100/400 Pa.
- p19 (biorxiv 2026, 10.64898/2026.06.16.732593) — growth-via-V0 ramp (×2 over G1) then CRC/ISE force
  fields, stop at 58.74%, AR targets 1.41 (2D)/~2 (3D).
- p30 (Phys Rev Research 2023) & p28 (review) — ring as an *emergent cortical wave* (flow→periodic-ring),
  ring closure `dD/dt ∝ D₀`.

### D2 — Cortex tension as an emergent actomyosin-density field (not a lumped scalar γ)  ★★★ (γ-floor, DCM side)
Replace DCM's lumped γ with a per-face density field + in-plane cortical flow:
- p30 (Phys Rev Research 5, 013208 2023) — `σ_act = (ζμ)₀·ρ/(ρ₀+ρ)` (Eq.3), surface-Stokes flow
  `∇²v+λ∇(∇·v)+∇f=v` (Eq.6), reaction-advection-diffusion density transport `∂ρ/∂t+∇·(ρv)=D∇²ρ+k_p−k_d ρ`
  (Eq.4); Péclet Pe=(ζμ)₀/(Dγ) as a grid-invariant contractility diagnostic (flow 38 / wave 115).
- p61 (ncomms12615) — active-stress buildup `σ_a=σ₀(1−e^{−t/τ_a})` (exponential γ ramp, not step-jump);
  per-face localized activated cap → heterogeneous cortex (cap/ring). p60 (ncomms10323) — turnover-gated
  dynamic γ ∝ n_bound·f_motor with sustainability diagnostic. This makes DCM's γ *emergent and dynamic*,
  and the density field is exactly what carries the contractile ring and blebs.

### D3 — Membrane–cortex coupling as breakable Bell linkers → mechanistic bleb module  ★★
DCM currently co-locates membrane == cortex. Split them with explicit breakable linkers:
- p14 (biorxiv 2025, 10.1101/2025.05.18.654456) — Bell-kinetics membrane-cortex linkers (Eq.5),
  two-mechanism bleb selection (Detachment vs Rupture) set by dimensionless R_C (coupling) and R_X
  (connectivity); **closed-form detachment-tension oracle γ_C^{D*} ≃ 0.35·(k_BT·ρ_H·R/r₀)·ln(ρ_A/K_m) ∝
  R_C**; portable surface-mesh potential set (per-triangle area, enclosed-volume, dihedral bending,
  steepened node-face contact); Young-Laplace split ΔP=2(γ_m+γ_c)/R (== DCM γ=ΔP·R/2, a gate); cortex
  network tension plateau ~0.2 mN/m ∝ R_C (bottom-up interphase-γ anchor).
- p55 (BPJ 2015, 10.1016/j.bpj.2015.06.016) — directional bleb from a locally broken γ field:
  per-face heterogeneous γ + cortex-damage scalar + Marangoni density advection; γ(θ)=γ_base[1+A cosθ].

### D4 — Aggregate-compaction RATE / cadherin-T1 clock = the mechanism for our open problem  ★★★
This directly resolves the open item I left in [[project-dcm-compaction-turgor-blocked]] /
[[project-dcm-cadherin-cluster-redesign]] (aggregate-σ compaction is a mechanical *end-state* only; the
real 24–48 h *rate* needs T1 rearrangement, not aggregate-σ pull):
- p70 (Miyazaki, ncomms 2020, 10.1038/s41467-020-16677-9) — **percolation-gated maturation kernel
  `τ_p = 2^N·τ`** (exponential-in-size slow timescale) as the mechanistic RATE law for loose→compact;
  occupancy `p = 1−[1−p_site^N]^{T/τ}` as a fine-grained replacement for any binary junction latch
  (the forbidden `cad_mult` switch); size-regime boundary R_c/L = log₂(T/τ) (derivable, no tuning). Also
  gives the passive nucleus a real emergent *positioning* role (cortex-vs-bulk tug-of-war).
- p07 (Matrix Biol 2019) & p03 (BPJ 2014) — junction/cortex **plasticity = Bell off-rate + rebind →
  permanent (not elastic) deformation**; turnover-limited plastic-flow rate `ε̇=k⁻·Δε_avg`; biphasic
  yield-stress creep (σ_th≈1 Pa). The bond-turnover driver the compaction-rate problem needs.
- p47 (Kim, PLOS CB 2022, 10.1371/journal.pcbi.1010257) — cadherin AJ cluster lifetime τ∝[cad]²; and an
  **independent confirmation that cadherin junction kinetics (~6 s max) do NOT gate 24–48 h compaction**
  — corroborates my S4 finding. Escalate compaction by *density/recruitment*, not by lowering off-rate.

### D5 — Nucleus + cytoplasm two-regime constitutive upgrades  ★
p32 (Kim review, J Biomech Eng 2020, 10.1115/1.4046863) — nucleus two-regime (chromatin soft <3 µm,
lamin-A load-bearing above, +2× lamina area reserve); poroelastic cytoplasm (exp<0.05 s → power-law,
D_p=E·ξ²/η); superelastic cortex plateau; power-law (soft-glassy) G'(ω),G''(ω) instead of single-time
Maxwell. p28 (review) — three-timescale active viscoelastic inclusion (elastic<1 s / viscous~10 s /
active~min), anomalous Poisson ν≈2. These feed the H.8/H.9/H.10 membrane/nucleus/cytoplasm gates.

### D6 — Confined single-cell motility (adhesion-free friction law)  ★
p76 (Sakamoto/Miyazaki, PNAS 2022) — cortex-contraction → interfacial-friction propulsion; `σ_act=α·v_act`;
closed-form V_drop with a geometric trade-off (friction ∝ h⁻¹, drag ∝ h⁻²) → optimum height ~60 µm. A
DCM confined-motility scenario + analytic cross-check (α≈1.7×10² Pa·s·m⁻¹).

---

## 3. Should a new engine be added?  → **No new engine. One new FF compartment + two cross-engine couplings.**

The corpus does not point to a third engine; it points to *completing the two we have*. Specifically:

- **New engine: not warranted.** Every mechanism in these 77 papers lands inside DCM (surface-mesh
  continuum cell) or FF (fine-grained filaments). Kim's own models *are* the FF paradigm; Miyazaki's
  reconstitutions are FF/DCM validation targets. Introducing a third runtime would duplicate FF.
- **New FF compartment — microtubule/kinesin (stage next after the actin turnover work).** p68 + p77 (+
  p08 kinesin chimeras, p40 aligned-MT durability, p28 spindle viscoelasticity) give a self-contained
  MT+kinesin+MAP module: three-state dynamic instability, angle-dependent zippering/catastrophe,
  extensile bundle sliding, stress→plus-end catch-bond feedback. FF already exposes `--microtubules`;
  this makes it first-class and unlocks spindle forces for the DCM division module (D1) and confined
  rotational-flow physics. **Recommend: build after FF-T1/T2/T3 land, PI-gated.**
- **Two cross-engine couplings worth formalizing** (not new engines): (i) **FF-cortex → DCM-γ**: use FF
  (T3/T4) to *derive* the emergent density-field γ that DCM (D2) consumes, replacing the lumped set-point
  and closing the γ-floor from both sides; (ii) **DCM-cell → FF-collagen ECM** (T8): the same Kim BD
  physics on both sides of the cell membrane, extending the DCM⊗ECM 1/r result to fine-grained collagen.

---

## 4. TAG validation-gap — what data the KB still needs (candidates, PI-authored)

New gates/contracts are PI-authored; these are **staged candidates** (full manifest:
`ffn_sim/outputs/tag_kb/SE_REGISTRATION_CANDIDATES_2026-07-10_kim-miyazaki.md`). All experimental values
are overlay/cross-check, never fit targets (oracle-is-crosscheck rule).

**Whole new validation lanes (no existing gate):**
1. **Cell division / cytokinesis** (p46,p22,p59,p28) — ~5% volume loss → ~60% axial elongation (MDCK
   n=22); ring constriction 0.39 µm/min per µm perimeter; F_ring≈F_spindle≈100 pN/side; elongation 20%
   (ring)/40% (ring+spindle); P*≈filament length 5.4 µm; R*≈ℓ_p. → **VG-DCM-division-*** family.
2. **Formin force-velocity** (p57,p62,p53) — Bell-gated elongation law (analytic-oracle candidate),
   mDia1 kon 21/3.6, koff 1.3/3.9 s⁻¹, p_c⁰ 0.29/0.61 (ATP/ADP), v 1.1–1.3 µm/s, 2.7 nm step. →
   **VG-FF-formin-fv** + **oracle candidate**.
3. **Crosslinker single-molecule + network rupture** (p67,p74,p02) — filamin unbind 70±23 pN / unfold
   57±19 pN; α-actinin-4 catch peak ~4 pN; network rupture 24.5 vs 6.5–8.1 Pa / 221 vs 63–129%. →
   **VG-FF-crosslinker-catch**, **VG-FF-network-rupture**.
4. **Cooperative/percolation contraction** (p61) — ρc≈0.56 µm⁻², telescopic law, εmax transition
   (n_H=11 is a FIT, overlay only). → **VG-FF-contraction-percolation**.
5. **Matrix viscoplasticity, stiffness-decoupled** (p07,p75,p39) — creep-recovery permanent strain
   80%→10% at invariant modulus; τ½ relaxation. → **VG-U1-ECM-viscoplasticity / stress-relaxation**.
6. **Confined actomyosin ring reconstitution** (p59) — duty² law, constant-volume invariant, R*≈Lp.
7. **Severing thresholds** (p14 F_sev=300 pN; p26 F_frag=500 pN, Tsuda 1996).
8. **Confined single-cell motility** (p76) — V_drop closed form (oracle candidate), α, optimum h.

**Confirms / tightens existing gates:** collagen G'=10/100/400 Pa at 1/3/5 mg/mL (p22) extends
VG-U1-G0; rBM ~100 Pa (p07) confirms VG-U1-G0 upper band; radial n=1.2 vs tangential n=2.2 (p13)
confirms + anisotropizes VG-U1-force-propagation ~1/r; ~0.2 mN/m cortex tension (p14) corroborates
interphase γ~0.15; catch bond (p74) extends VG-U2-catch-slip-bond from integrin (~30 pN) to crosslinker
(~4 pN).

**Conflicts to surface to PI (do NOT silently reconcile):**
- Crosslinker k0=**0.115 s⁻¹** (filamin-A) vs KB **0.066 s⁻¹** (α-actinin) — same Ferrer 2008 paper,
  different isoform. Keep as distinct species (p02,p03,p53,p69).
- Myosin stall **5–5.7 pN per head** (p05,p14,p26,p69) vs `PARAM-F_stall_motor` default **2 pN "per
  motor"** — per-head vs per-motor; NMII v0 **140 nm/s** (p05,p14) vs `PARAM-v_unloaded` default 100.
- α-actinin K_d/F-actin 1.4 µM, rupture 1.4–80 pN (p14) — anchor F_b before registering.
- Several sim-origin constants sit behind **SI tables not in the downloaded PDFs** (Mulla Suppl. Table 1;
  Fan Suppl. Table 5 / github.com/ktyman2/liverCancer; Nam Table S1; Matsuda Table S1) — **fetch SI
  before KB-registering any force constant** from p74, p75, p22, p26.

---

## 5. How this resolves three standing project open problems

- **γ-floor** ([[project-gamma-floor-layered-resolution]]): Kim's answer is that sustained cortical
  tension is **turnover- and crosslinker-clutch-limited, and requires filament buckling** — not a
  motor-magnitude deficit. Concrete levers: T1 (turnover+severing), T4 (clutch sustainability law
  ∝R_ACP/R_M), T3 (bundle vs √N_M network force scaling, p41), D2 (emergent density-field γ). This is a
  mechanism, not a tuning knob — consistent with the no-magic-number rule.
- **Aggregate-compaction rate / cadherin T1** (my [[project-dcm-cadherin-cluster-redesign]] left this
  open as "(b) T1 rearrangement escalation needed"): p70 gives the exact RATE law (percolation
  maturation τ_p=2^N·τ), p07/p03 give the bond-turnover plasticity driver, and p47 independently
  confirms cadherin kinetics don't gate 24–48 h compaction. → D4 is now a mechanistic build, not an
  open question.
- **Crawl-speed grid-drag** ([[project-ff-crawl-mechanism]], native crawl 500× slow, `--bulk-drag`
  destabilizes): T6 replaces the count-scaling drag with Kim's geometry-derived per-segment cylindrical
  Stokes drag at physiological μ (p03/p13/p14/p22) — a physical, grid-invariant reframing.

---

## 6. Recommended sequencing (PI to ratify)

1. **FF-T1** turnover + buckling-severing (unblocks γ-floor + stress-relaxation) — highest ROI.
2. **FF-T2/T3/T4** catch-slip crosslinker + PCM myosin + clutch-sustainability (γ-floor mechanism).
3. **DCM-D1** cytokinetic division from p46 (drop-in) — DCM's biggest capability gap.
4. **DCM-D2 + D4** emergent-density γ and percolation compaction-rate clock.
5. **FF-T6** cylindrical drag (crawl-speed) and **FF-T5** force-dependent formin.
6. **DCM-D3/D5**, **FF-T7/T8**; then stage **FF-T9 microtubule compartment** (PI-gated).

Each step must land at **native + full-compartment** scale and pass its (PI-authored) gate before it is
reported — no coarse/stripped validation, no gate-loosening, no fitting.
