# Mechanism Audit & Literature-Strengthening Report — ffn_cellsim

**Date:** 2026-05-30
**Scope:** Consolidation of a 33-agent audit (8 subsystem audits + 20 per-mechanism literature dossiers + 5 cross-cutting critiques: completeness, wiring_register, lumped_violations, param_verify, fa_design).
**Trigger:** H.4 focal adhesions (`ffn_sim/bridge/`) are implemented + unit-tested but NOT wired into `ffn_sim/cell/cell.py`, so the integrated cortex is a floating shell with no substrate anchor. KU-3.5 cortical tension fails ~3,000× under target and KU-5.1 dendritic density fails ~10,000× under. This audit found every other FA-class gap and strengthened every mechanism against current literature.

## Executive summary

The fine-grained mechanistic core of `ffn_cellsim` is structurally sound: nearly every individual mechanism is the correct field-standard mechanistic choice (Mikado/Broedersz networks, Bell-Evans slip, Hill f-v, Stam-Hocky minifilament, Leimkuhler-Matthews BAOAB + constrained g-BAOAB/SHAKE/Fixman, angle-harmonic Arp2/3, Brownian-ratchet elongation, force-dependent capping, WCA excluded volume), and the published closed-forms (Storm-MacKintosh, Hill, Pereverzev, Chan-Odde) are correctly held as acceptance oracles in most places. The dominant failure is **architectural, not parametric**: the entire H.4 focal-adhesion bridge and the entire H.1 ECM Mikado substrate are fully built + tested **but have zero call sites in `Cell.build()`** — verified directly: zero `ecm.*`/`bridge.*` imports in `cell.py`, `with_fa` (line 584) is a dead no-op, `tag_ranges()` hard-zeros the FA slot (`offsets['fa_integrin']=(cur,cur)`, line 902). The cortex therefore has no mechanical ground to pull against (→ KU-3.5) and lamellipodial barbed ends carry no load (→ KU-5.1). The completeness critic additionally surfaced a **class of biophysical mechanisms with zero runtime implementation**: actin turnover (cofilin/severing/treadmilling), intracellular pressure / enclosed-volume (the actual KU-3.1 mechanism), explicit plasma membrane, talin/vinculin mechanosensing (dead config), substrate compliance spring, and force-dependent debranching. The lumped/param critics found a small number of real runtime violations: the myosin Hill-stepping **binned-rest-length ratchet** makes emergent tension grid-dependent; the ERM tether is a static Hookean pin where the field uses a Bell-Evans slip clutch; and **filamin (the 70% majority crosslinker) is modeled as slip when it is a documented catch bond** (with a Sanity-Gate test locking the wrong sign). A handful of real parameter/provenance corrections were confirmed (myosin v0 ~5× too fast; capping δ_cap ~9× too small; α-actinin k_off0 ~15× too fast; capping miscited to a nonexistent "Funk 2022"; branch angle 72° vs 70°). The FA-integration design (§6) is the restart deliverable: a layered `substrate→integrin→talin→vinculin→clutch→cortex/myosin↔membrane→lamellipodium` chain, every link explicit, all closed-forms oracle-only.

**Headline counts.** Wiring gaps: **6 critical, 7 high, 3 medium, 2 low** (FA + ECM + EV are the floating-shell triad). Completeness gaps (zero-runtime mechanisms): **2 critical, 4 high, 3 medium/low**. Real runtime architectural/lumped violations: **2 critical, 3 high, 4 medium**. Real parameter corrections needing action: **~8**. Distinct paywalled papers flagged for gbook/KAIST fetch: **~30**.

---

## §1 Wiring-gap register

Every implemented-but-unwired or specified-but-unimplemented force-pathway, ranked critical → low. "Wired" = reachable from `Cell.build()` in the integrated whole-cell simulation. Verified against the live tree.

| # | Mechanism | KU | Impl. | Wired | Sev. | What whole-cell behaviour/gate it silently breaks | Evidence (file:line) |
|---|-----------|----|-------|-------|------|----------------------------------------------------|----------------------|
| 1 | H.4 FA bridge subsystem (`build_h4_state`, `build_h4_simulation`, `IntegrinBondUpdater`, `build_minifilaments_for_fas`, `FAGrowthMonitor`) | KU-2.4/2.5/2.17/2.18, KU-3.5 | yes | **no** | **critical** | Floating-shell root cause: cortex has nothing to pull against → KU-3.5 γ≈1.4e-4 vs [0.35,0.65] mN/m (~3,000×); barbed ends carry no load → capping stays maximal → KU-5.1 ~0.012 vs ~100/μm² (~10,000×); whole rigidity-sensing axis inert | `bridge/{fa,integrin_bonds,motor,fa_growth}.py` (fully tested); `cell.py:584` `with_fa` dead no-op, `:902` `offsets['fa_integrin']=(cur,cur)`, `:866` `fa_integrin=0`, `:638` `fa` slot never populated; zero `bridge.*` import in `cell.py` |
| 2 | H.1 ECM / Mikado substrate (`mikado`, `cross_links`, `shear_protocol`, `equilibrate`) | KU-1.1–1.30 | yes | **no** | **critical** | Other half of floating shell: FA integrins have no ECM ligand to anchor to; no E_substrate input to the clutch → no rigidity threshold; compounds KU-3.5/KU-5.1 | `ecm/mikado.py:393-431` etc. (≤1e-6 energy oracle, tested); zero `ecm.*` import in `cell.py`; no `with_ecm` flag, no builder path |
| 3 | Substrate anchor at z=0 (fixed mechanical ground) | KU-2.4 / Chan-Odde | **no** | **no** | **critical** | No fixed reference frame → myosin stress integrates to ~zero net tension (the direct KU-3.5 deficit) | no immobile/position-restrained ligand layer anywhere; `fa.py` ligand is a single fixed point under each FA (no group exclusion / pinning in `cell.py`) |
| 4 | Molecular-clutch loading loop (cortex actin → integrin → talin → ligand → anchor, resisting myosin retrograde flow) | KU-2.x / Chan-Odde | parts exist, **uncoupled** | **no** | **critical** | The tension-building engine: parts exist (myosin wired, integrin updater implemented) but no clutch bond couples cortex actin to substrate → myosin only locally rearranges a floating shell | `cortex/myosin.py` wired; `bridge/integrin_bonds.py` implemented; no clutch bond between cortex actin and substrate in `cell.py` |
| 5 | Explicit membrane object → lamellipodial barbed-end load (Brownian-ratchet brake) | KU-5.1/5.2/5.6 | **no** | **no** | **critical** | Barbed ends see F=0 so `exp(-Fδ/kT)` capping/elongation laws are inert → KU-5.1 collapse; lumped into a CFL-magic-number WAVE plane pin | no Helfrich/membrane object; `lamellipodium.py` `k_wave_pin=1e-4 N/m` (no lit anchor); ERM is cortex-radial only |
| 6 | FA per-FA minifilaments (`build_minifilaments_for_fas`) + bond/WCA param registration for bridge types | D5, D7 | yes | **no** | **critical** | `with_fa=True` would crash at step 0 (no integrin_ligand/motor bond params, no WCA cross-pairs registered) | `bridge/motor.py:215-374`; only caller is `fa.py` standalone; `cell.py:301-336` registers cortex/xlink/myosin bonds only |
| 7 | Talin multi-domain WLC mechanosensor (rod unfolding 5–25 pN, R3 ~5 pN, exposes VBS) | KU-2.6 | **dead config** | **no** | high | Open catch/reinforcement loop; no rigidity-sensing threshold → FA chain under-builds tension on stiff substrate even once wired | `phase1_h4.yaml` `talin{enabled:false,...}` plumbed into `ResolvedH4` (fa.py:157-160,326-329); NO updater consumes it |
| 8 | Vinculin force-gated reinforcement (discrete recruitment + vinculin-actin directional catch) | KU-2.7 | **dead config + lumped** | **no** | high | Adhesions can't mature/reinforce under load; reinforcement is a forbidden scalar proxy `k_int_eff=k_int_bare·(1+α·N_vin)` | `phase1_h4.yaml` `vinculin{enabled:false,...}`, `integrin.alpha_vin=0.05`; no updater; scalar k_eff is the lumped violation |
| 9 | FAGrowthMonitor → growth feedback (KU-2.17 Hill threshold ~50 pN → N_engaged≈10) | KU-2.17 | observer only | **no** | high | FA size/number non-adaptive; observes n_engaged but never feeds back into growth | `bridge/fa_growth.py:83-206`; instantiated only in `test_ku217_fa_growth.py`; no Cell attribute/monitor |
| 10 | Substrate compliance spring (E_substrate series spring) | Chan-Odde / Bangasser | **no** | **no** | high | Rigid ligand → biphasic traction-vs-stiffness optimum (~1 kPa) and rigidity sensing cannot emerge; no E sweep driver | `fa.py` ligand fixed at z=0; no `k_sub`/`E_substrate` runtime element |
| 11 | Equilibration prelude (`equilibrate_no_shear` soft-start + BAOAB drain) not called in Cell.build | M2-rest (PI 2026-05-20) | yes | **no** | high | If ECM/FA wired without it, construction-time xlink overlap → step-0 force explosion → BAOAB divergence | `ecm/equilibrate.py:179-309`; used only in `test_h1_shear.py`/`test_ku130.py` |
| 12 | Cortex excluded-volume factory orphaned + WCA shift-mode bug | D7 / KU-1.2 | yes | partial/**no** | high | `build_excluded_volume` (cortex.py:875) never called by `build_cortex_full_simulation`; HOOMD LJ defaults to mode='none' → truncated-not-shifted WCA injects +ε force discontinuity into BAOAB; if KU-3.5 measured with EV off, second compounding cause | `cortex.py:875` orphaned + duplicated at `cell.py:111`; LJ shift mode unset; σ passed as bare 2·R_bead not ×40 bundle diameter |
| 13 | ERM modeled as static Hookean pin, not Bell-Evans slip clutch | KU-3.18 | wired but **under-resolved** | partial | high | Wired but lumped: membrane-cortex tension not transmitted as a regenerating force-sensitive clutch population (~100 linkers/μm²); weakens membrane-load leg feeding KU-5.x and KU-3.1 confinement | `cortex/erm.py` fixed-k radial spring, no `k_off`/`exp()`; `k_ERM=1e-4` softened 1000× (but 1e-4 is the correct per-linker value, §3) |
| 14 | Variable-length cortex filaments (1–5 μm) | H.3 brief | yes | **no** (config parsed, ignored) | medium | Config `variable_length.enabled=false` never checked; blocks Chugh-2017 non-monotonic tension-vs-length validation of KU-3.5 | `cortex.py:1176-1231`; `phase1_h3.yaml:128-140`; Cell.build never routes there |
| 15 | Global cross-subsystem CFL / dt reconciliation | KU-1.26 | partial | **no** | medium | ECM (stiffer bonds) under cortex dt may violate ECM CFL; SHAKE convergence not asserted vs dt; no unified `dt=min(τ)` gate | `mikado.py`/`cortex.py` compute dt_cfl independently; `baoab.py` asserts dt>0 but τ_min delegated to caller (audit 5, sev high) |
| 16 | Shear-strain (Lees-Edwards) protocol unwired to whole-cell | KU-1.30 #2 | yes | **no** | medium | No whole-cell oscillatory-shear rheology path (validation capability, not a cell mechanism) | `ecm/shear_protocol.py:102-150`; `test_h1_shear.py` only |
| 17 | Static crosslinker bonds (legacy stage-1) | — | yes | **no** (default off) | low | Legacy; dynamic D2 kinetics is the live deliverable | `cortex.py:700-755`; `phase1_h3.yaml crosslinkers.enabled=false` |
| 18 | Cadherin junction catch bond (cell-cell) | KU-4.2 | **no** (empty dir) | **no** | low | Single-cell scope: no second cell; correctly deferred to Phase 2 | `junction/` only `__init__.py`; sole artifact is read-only v1 oracle `validation/oracles/junction/types.py` |

**Subsumption note.** The ECM audit additionally lists ~7 *individual* H.1 mechanisms (Mikado topology, xl bonding, ECM WCA, bond/angle harmonics, ECM BAOAB, ECM CFL, ECM validation gates) each as `wired=false, severity=high`. These are not independent gaps — they are all symptoms of row 2 (H.1 not wired into Cell.build) and become live the moment the `with_ecm` builder path exists. Counted once.

**Wiring-gap tally:** critical 6 · high 7 · medium 3 · low 2.

---

## §2 Architectural-principle violations (lumped/proxy as runtime) + magic numbers / placeholder bands

Architecturally clean in the large (BAOAB, constrained g-BAOAB+SHAKE+Fixman, Mikado, cortex topology, Stam-Hocky topology, dynamic Bell-Evans crosslinkers, lamellipodium branch/elongate/cap are all genuine particle/bond dynamics; Storm-MacKintosh/Hill/Pereverzev correctly held as oracles in most places). The seven runtime violations + magic numbers found:

| Location | Violation | Severity | Fix |
|----------|-----------|----------|-----|
| `cortex/myosin.py:720-738, 965-1002` Hill-stepping kernel | **Grid-dependent runtime construct.** Hill stepping does NOT translocate the head; it ratchets a **binned bond rest-length** (r0) toward bin 0 via a fractional accumulator (`KU-3.5 단계-4 fix`). True per-tick displacement v·batch_dt ≈ 70 pm ≪ bin_width ≈ 33 nm (= max_bind_dist 327 nm / n_bins 10), so emergent contractile tension depends on the arbitrary grid (n_bins, bin_width), not physical displacement — directly feeds the failing KU-3.5 gate. | **critical** | Replace bin-ratchet with physical sub-bin processive stepping (translate the bound anchor along the actin segment by v·batch_dt). Add an n_bins-sweep grid-invariance Sanity Gate (γ_total invariant within tolerance). Keep `hill_velocity_clamped` oracle-only. Needs PI (frozen-integrator-adjacent + mechanism change). |
| `cortex/erm.py` ERMHarmonic | **Lumped elastic proxy where field standard is a stochastic clutch.** Static force-independent Hookean radial pin `U=½k(|r|−R)²`, no unbinding/turnover/density — but the membrane-cortex tether is a Bell-Evans SLIP clutch population (Korkmazhan-Dunn 2022; Alert 2015), and CLAUDE.md's own architecture table names Bell-Evans as the right choice over a static proxy. | high | Convert to a stochastic Bell-Evans SLIP clutch population (reuse `IntegrinBondUpdater`): per-linker k=1e-4 N/m, force-dependent off-rate, gated rebinding, density ~100/μm². Keep stiff element on constrained-BAOAB. |
| `cortex/crosslinkers.py:456-458, 592-598` + `test_crosslinkers.py:238-247` | **Wrong mechanism class for the majority crosslinker, with a Sanity-Gate test locking the wrong sign.** Filamin (70% of the cortex pool, KU-3.19) is modeled as pure Bell **slip** (`exp(+F·x_β/kT)`), and a Sanity-Gate test asserts `dk_off/dF>0` using filamin params — but filamin is a documented **catch** bond (Gieseke/Rief 2013, Ehrlicher 2011, Rognoni 2012). Under load a pure-slip cortex *sheds* crosslinks where the real cortex *stiffens* → plausible contributor to the KU-3.5 collapse. | high | Flip filamin to the two-pathway Pereverzev catch-slip form already implemented for integrins. The slip-sign Sanity-Gate test must be revised to a catch-slip contract — **gate-contract change, surface to PI, do not edit inline.** Reassess α-actinin isoform (ACTN4 cortex is also catch). |
| `cortex/myosin.py` head-actin bond (separate from above) | NMII head-actin modeled Bell-Evans **slip-only**; literature shows NMII has biphasic catch-slip character. Slip-only under-sustains tension under load. (Lower priority than the binning issue.) | medium | Introduce a load-dependent duty-ratio / catch term as oracle first, then runtime. |
| `cell.py` `build_cortex_full_simulation` inter-subsystem LJ disabled | Myosin×actin, xlink×actin steric exclusion **globally disabled** to allow Bell-Evans binding without WCA blockade. Deliberate + documented, but permits non-physical overlap at higher density. | medium | Binding-aware potential (attractive well at bond range + WCA core) or per-bonded-pair exclusion lists rather than global disable. |
| `phase1_h5.yaml:41-43` `angle_branch_k` Arp2/3 branch-angle stiffness | **Magic number in the live force path.** `k_angle=1e-19 J/rad²` explicitly tagged "literature TBD"; no acceptance band on the daughter-angle distribution, so unconstrained by any gate. | medium | Derive `k_angle=kT/Var(θ)` from cryo-ET 68±9° (Fäßler 2020) → ~1.7e-19 J/rad²; Magic-Number Block + branch-angle-distribution oracle. Also fix runtime t0 72°→70° (§3). |
| `phase1_h5.yaml:30-32` `k_wave_pin` WAVE membrane pin | **CFL-driven magic number, no literature anchor** (Bieling/Funk don't specify it); live runtime confinement force; brief intent was to NOT parameterize it. | medium | Replace with the explicit membrane object (§6 #5), OR document explicitly as a non-mechanistic CFL-bounded confinement device with a swept-insensitivity Magic-Number Block. Surface to PI. |
| `phase1_h1.yaml:91-115` ECM acceptance bands rebanded for D4 | **Authorized gate-softening.** G_0 [15,200]→[1,50] Pa and |β|∈[1.5,2.5] to absorb D4 (21-bead) rigidity-percolation prefactor uncertainty. PI-authorized 2026-05-21 + first-principles-derived (estimate ~5.6 Pa, band ~10× wider) — flagged for transparency. | low | Keep PI-ratified; optionally re-derive the D4 prefactor to tighten. |

**Confirmed NOT violations (oracle-clean / derivable):** Hill f-v as a per-head closed-form (ensemble emerges) — but note the *binned-r0 stepping* surrogate IS the violation above. Myosin `head_actin_max_bind_dist=327 nm` is geometrically derived (√((ℓ₀/2)²+capture_perp²)), not magic. WCA cutoff 2^(1/6)σ, LJ ε=0.5 kT (D7), H.2 equipartition target 0.9898 kT, all H.1/H.3 acceptance bands — first-principles. BAOAB/constrained-BAOAB integrators are genuine mechanistic. `shake_tol=1e-10` and `INT32_GUARD=1e8` are empirically justified numerics (no lit citation, low severity).

---

## §3 Parameter-verification table

Disagreements the `param_verify` critic judged **REAL and actionable** (recomputed analytic landmarks; not taken on the lit agents' word):

| Param | Project value | Literature value | Source | Verdict | Recommended action |
|-------|---------------|------------------|--------|---------|--------------------|
| Filamin mechanism class | slip (Bell, +x_β) | **CATCH** (off-rate falls then rises with F) | Gieseke/Rief 2013; Ehrlicher 2011; Rognoni 2012 | **real — wrong sign, majority crosslinker** | Flip to two-pathway catch-slip; revise the slip Sanity-Gate test → catch contract (PI gate-contract) |
| Capping δ_cap | 0.3 nm | 2.7 nm (= δ_elong; shared barbed-end gap) | Li/Bieling eLife 2022 | **real — ~9× too small, internally inconsistent** | Set δ_cap = δ_elong = 2.7 nm; fix the "0.3 pN" units mislabel (length, not force). Single most useful capping fix for KU-5.1 |
| Capping citation | "Funk 2022 Nature" | Li/Bieling eLife 2022 e73145 (no such Funk paper) | Li 2022; Funk 2021 Nat Commun | **real — miscitation** | Re-attribute capping force-law to Li/Bieling 2022; keep Funk 2021 for branching/NPF only; fix every "Funk 2022" string |
| Myosin v0 | 1.0 μm/s | NMIIA ~0.1–0.3 μm/s (~5× too fast) | Kovács 2003; Wang 2003; Norstrom 2010 | **real — ~5× over, source-mismatch** | Set v0 ≈ 0.2 μm/s for NMIIA + correct citation, OR document as a deliberate mesoscale effective value with derivation (PI) |
| α-actinin k_off0 | 1.0 /s | ~0.066/s single-molecule (Ferrer 2008); ~0.66/s bulk (Wachsstock) | Ferrer 2008; Wachsstock 1994 | **real — ~15× too fast, miscited** | Re-anchor to single-molecule (Ferrer 2008 ~0.066/s) or justify bulk; correct citation; re-derive with the filamin catch fix |
| Arp2/3 branch angle t0 | 72° (5π/12) | 70° (68±9° cryo-ET) | Mullins 1998; Fäßler 2020; Ding 2022 | real (low sev) — stale 2° offset | Set t0 = 70° (1.2217 rad); derive k_angle from the 68±9° distribution |
| Integrin catch peak F* | 6.99 pN (analytic from Pereverzev) | ~10–30 pN (α5β1-FN) | Kong 2009; Roca-Cusachs 2009; Chen 2017 | real — but **self-flagged + PI-deferred to Phase 2** (gate written against analytic value, honest) | Migrate Pereverzev two-pathway → Kong 2009 two-state (F*~30 pN) when wiring FA |
| Talin unfold force / vinculin | absent (config disabled) | R3 ~5 pN, ladder 5–25 pN; vinculin = discrete catch, not scalar | Yao 2016; del Rio 2009; Tapia-Rojo 2020; Huang 2017 | real (architectural — needs the §6 implementation) | Pin per-domain forces when implementing talin/vinculin (§6 #3,#4) |

**Real-but-architectural (fix by wiring/adding mechanism, NOT by re-parameterizing):**

| Param | Project (measured) | Literature | Source | Verdict |
|-------|--------------------|-----------|--------|---------|
| Cortical tension γ | 1.4e-4 mN/m | 0.1–1.0 mN/m (~0.5) | Tinevez 2009; Salbreux 2012; Chugh 2017 | inputs correct; floating-shell + no turnover/membrane produce wrong emergent γ |
| Dendritic density | 0.012 /μm² | ~100 /μm² | Abraham 1999; Bieling 2016 | caused by missing barbed-end load + no turnover, not a parameter error |
| Substrate stiffness E | absent (fixed ligand) | 0.1–100 kPa (optimum ~1 kPa) | Chan-Odde 2008; Bangasser 2013 | needs explicit compliance spring (§5, §6 #6) |
| Intracellular pressure ΔP | absent | ~40 Pa interphase → ~400 Pa metaphase | Stewart 2011; Fischer-Friedrich 2014 | needs enclosed-volume term (§5) |

**Benign / confirmed-OK (no action beyond citing):**
- **k_ERM 1e-4 vs brief 0.1 N/m** — softened 1000× for CFL (PI-ratified), but **1e-4 N/m is actually the correct per-linker MCA value** (Alert 2015); only a provenance/comment fix (the brief's 0.1 N/m was the mis-scaled areal number).
- η_water 1e-3 vs ~0.69e-3 Pa·s at 310 K — conservative (higher drag → smaller dt; affects kinetics, not equilibrium).
- Collagen fibril radius 50 nm vs 50–250 nm (Ushiki 2002; Lindström 2010) — within range (low end).
- Actin L_p 17 μm vs 17.7 μm at RT (Gittes 1993; Ott 1993) — within 9–20 μm spread; document the 310 K vs RT kBT scaling.
- κ_B 7e-26 vs ~7.3e-26 N·m² — within uncertainty.
- NMII k_off0 10/s, x_β 0.6 nm (Veigel 2002) — within range (verify page/table).
- Filamin k_off0 0.1/s (Furuike 2001 — but Furuike measures Ig *unfolding*, wrong observable; re-cite).
- F_stall 0.5 pN/head — within optical-trap NMII range; a/F_s=0.5 lacks an NMII anchor (muscle ~0.25) — lower-severity, document.
- Myosin backbone 700 nm, density 100/cell (3/μm², Salbreux 2012), bond/angle k = μ/ℓ₀, κ/ℓ₀, WCA params — confirmed.

---

## §4 Per-mechanism literature strengthening (20 dossiers)

| # | Mechanism | Best mechanistic choice? | Field-standard assessment | Better alternative | Top strengthening recs |
|---|-----------|--------------------------|----------------------------|--------------------|------------------------|
| 0 | Mikado ECM (Broedersz/Storm) | No (refine within class) | The CLASS (discrete network, emergent stiffening) is field-standard & correct | Move to a 3D off-lattice fiber network (z_c=6 in 3D vs 4 in 2D, so the sub-isostatic margin that sets G_0/onset-strain is wrong in 2D) — Licup 2015, Jansen 2018; add tension/compression asymmetry (fiber buckling) | Decide+document the ECM target biopolymer (params are internally inconsistent: actin L_p=17 μm + collagen r=50 nm/E=1.1 MPa + entropic 3/2 oracle); promote stiffening to a K∝σ^{3/2} oracle (Storm 2005; Licup 2015; Sharma 2016); plan 3D off-lattice network for Phase 2. Flagged param: strain-stiffening band |β|≈1.5 vs theoretical 3/2 (agree=No per agent, though numerically the same — re-anchor band to 3/2). |
| 1 | Actin L_p / bending modulus | **Yes** | Unambiguous WLC standard (AFINES/Cytosim/MEDYAN all use it) | None (optional twist/torsion later) | Document L_p 17 μm @310 K vs 17.7 μm @RT (Gittes 1993); keep κ=L_p·kBT internally consistent |
| 2 | BAOAB Leimkuhler-Matthews | **Yes** | Current best-practice overdamped Langevin; faithful LM-limit | None at scale (g-BAOAB already present for rigid path) | Make equipartition + FDT diffusion checks **enforced** CI gates; cite L-M 2013 |
| 3 | Constrained BD (SHAKE + Fixman) | **Yes** | Legitimate fine-grained rigid-bond approach; correct stiff-bond choice | r-RESPA / implicit integrator for stiff-not-rigid (optional) | Add equipartition gate on constrained trimer to validate Fixman sign +1 (currently deferred); cite Fixman 1974 |
| 4 | Bell-Evans slip off-rate | No (class caveat) | Correct slip primitive for α-actinin / NMII / capping / elongation | Where the real bond is CATCH (filamin, integrin, vinculin, cadherin) slip is the WRONG class | Reclassify each bond; switch filamin → catch-slip; add Dudko-Hummer-Szabo high-force oracle |
| 5 | Integrin catch (Pereverzev two-pathway) | No | Catch is real & essential; two-pathway no longer best | **Kong 2009 / Chakrabarti-Hummer two-STATE allosteric** (bent↔extended) → right F*~30 pN | Migrate Pereverzev→two-state; pin F*~30 pN; add Roca-Cusachs 2009 / Chen 2017 |
| 6 | Molecular clutch (Chan-Odde) | No (Chan-Odde is lumped) | Right target; Chan-Odde/Elosegui-Artola are mean-field → ORACLES only | Explicit per-clutch: integrin two-state catch + talin multi-domain WLC + discrete vinculin + substrate compliance spring | Implement explicit chain; add E_substrate spring → biphasic optimum ~1 kPa emerges; Chan-Odde oracle-only |
| 7 | FA growth / talin-vinculin reinforcement | No (dead config) | Correct current picture; growth laws (Riveline/Balaban) = oracles | Full Tapia-Rojo 2020 talin folding kinetics | 13-domain serial WLC, R3 first ~5 pN, +60–80 nm/unfold exposing VBS (Yao 2016); discrete force-gated vinculin (Huang 2017 catch), not k_eff scalar |
| 8 | Hill f-v for myosin | No (per-head runtime) | Hill is an empirical *muscle-level* relation; per-head runtime use is a phenomenological closed-form — OK as oracle, not best mechanistic | Discrete strain-dependent cross-bridge cycle (Huxley 1957 / Duke 1999 / Erdmann-Schwarz 2013) so f-v emerges | **Demote Hill from per-head runtime to ensemble oracle**; per-head attach/powerstroke/detach cycle; fix v0 (~0.2 μm/s); anchor a/F_s to NMII or document as effective (muscle ~0.25, current 0.5 unanchored) |
| 9 | Stam-Hocky minifilament | No (refine within class) | Field-standard fine-grained NMII (AFINES/MEDYAN-class) | Topology is right; refine stepping + stoichiometry | **Replace binned-r0 ratchet with physical translocation (grid-invariance)**; reconcile heads/side 10→~15 (Billington 2013, ~30/minifilament); fix v0 |
| 10 | α-actinin + filamin crosslinkers | No (within-class fix) | Explicit divalent force-dependent xl IS the standard (right side of the line) | Within-class: filamin → catch; optional WLC spacer | Switch filamin off-rate slip→catch (Ehrlicher 2011; Rognoni 2012; Gieseke 2013); re-anchor α-actinin k_off0 to Ferrer 2008 |
| 11 | ERM membrane-cortex tether | No (static pin ≠ standard) | Current static Hookean pin is NOT field-standard | **Explicit stochastic Bell-Evans SLIP clutch population** (k~1e-4 N/m, ρ~100/μm², gated rebind) | Convert pin → clutch (Korkmazhan-Dunn 2022; Alert 2015); re-anchor k_ERM as the per-linker value |
| 12 | Cortical tension generation | **Yes** | Emergent tension from explicit motors+network IS the standard | None (additions, not replacement) | **Add actin turnover** (Chugh 2017 — the dominant missing lever); add FA anchor; cite Salbreux 2012 |
| 13 | Mitotic rounding / blebbistatin | No (missing pressure) | Right intent; missing the central enclosed-pressure term | Add osmotic/enclosed-volume term so ΔP=2γ/R is mechanistic | Enclosed-volume term (Stewart 2011); blebbistatin = graded per-head duty-ratio reduction (Ramanathan 2015), not binary; note metaphase γ~1.6 mN/m > KU-3.5 ceiling (gate question) |
| 14 | Arp2/3 dendritic branching | No (within-class fix) | Correct three-part explicit treatment; γ-Phase-1 CP-NPF current | Add force-dependent debranching | **Add debranching** (Pandit 2020 junction Bell-Evans, 0–2 pN→~100×; also bounds n_total); fix t0→70°; derive k_angle from 68±9° (Fäßler 2020) |
| 15 | Branched network f-v + WAVE/WRC nucleation | No (missing membrane) | Right fine-grained intent; γ-Phase-1 current | Explicit deformable membrane + WRC trans-autoinhibition/Rac1-PIP3 gating | Add membrane object (barbed-end load → Bieling 2016 f-v/load-adaptation); replace plane pin with WRC activation gate |
| 16 | Force-dependent capping (CP/CapZ) | **Yes** | Correct force-dependent-rate runtime class | None (refinements only) | Set δ_cap=δ_elong=2.7 nm (Li/Bieling 2022); re-attribute from "Funk 2022"; make k_cap0 = k_on·[CP] derivable (~0.63/s @100 nM, Schafer 1996) |
| 17 | Barbed-end elongation (Brownian ratchet) | No (missing membrane load) | k_on0[G]exp(−Fδ/kT)−k_off is field-standard; δ=2.7 nm sound | None (refinements) | Provide physical membrane load (so F is real); dynamic monomer reservoir [G]; cite Demoulin 2014, Footer 2007 |
| 18 | Excluded volume (WCA) | **Yes** | Universal field-standard repulsion primitive | None | **Set LJ shift mode so U(2^(1/6)σ)=0** (true WCA) + force-continuity unit test; derive σ from ×40 bundle diameter, not bare 2·R_bead (Weeks-Chandler-Andersen 1971) |
| 19 | Integrated force-transmission chain | No (UNWIRED) | Exactly the right architecture; matches field consensus | None — the architecture IS the standard; fix is wiring + per-link fidelity | Wire the full chain into Cell.build (§6); keep Chan-Odde/Bell/Hill/Riveline/Balaban as oracles only |

---

## §5 Completeness gaps (entirely missing / under-modeled biophysical mechanisms)

From the completeness critic — a CLASS of processes with **zero runtime hits** (regex-confirmed across `cell.py`, `lamellipodium.py`, `cortex/*`, `bridge/*`, `ecm/*`):

1. **(critical) Actin turnover — severing / depolymerization / treadmilling (cofilin/ADF).** Zero hits for `cofilin|treadmill|sever|depolymeri`. The cortex is a frozen, permanently-connected network (xlinks turn over, filaments never disassemble); the lamellipodium grows + recycles but has no depolymerization channel. The active-gel literature (Chugh 2017, Saha 2016, Fritzsche 2013) identifies turnover (half-life ~5–30 s) as the **single dominant determinant of steady-state cortical tension** — without it KU-3.5 cannot be reproduced *even after the FA anchor is wired* (a turnover-free network jams or never builds steady stress), and KU-5.1 has no disassembly to balance growth. → Add a force/age-aware disassembly Updater (cofilin severing + pointed-end depoly) to cortex AND lamellipodium; Magic-Number-Block the turnover time from primary data; co-requisite for the tension/density gates.

2. **(critical) No intracellular pressure / enclosed-volume / Laplace balance — KU-3.1 is non-mechanistic.** No pressure/osmotic/hydrostatic/Laplace/enclosed-volume term anywhere (the only "pressure" hits are the unrelated lamellipodium abortive threshold). KU-3.1 mitotic rounding is physically ΔP=2γ/R (Stewart 2011) — the H.3 brief's own KU-3.5 criterion reads "Laplace law from internal pressure" — yet rounding is currently held only by the ERM radial pin + an aspect-ratio band. → Add an enclosed-volume/osmotic-pressure term on the closed shell so ΔP=2γ/R emerges and aspect≤1.2 is demoted to an oracle (Fischer-Friedrich 2014: metaphase ~400 Pa, ~1.6 mN/m). **Surface to PI:** metaphase ~1.6 mN/m exceeds the KU-3.5 ceiling (1.0) — KU-3.5 may be interphase-biased (a gate-contract question).

3. **(high) Talin rod unfolding + vinculin reinforcement is dead config** (KU-2.6/2.7). `talin{enabled:false}` / `vinculin{enabled:false}` plumbed into `ResolvedH4` fields with NO updater; the only reinforcement is the forbidden scalar `k_int_eff=k_int_bare·(1+α·N_vin)`. This is the Elosegui-Artola 2016 mechanism that lets tension BUILD on a stiff/anchored substrate — without it the FA chain stays catch-bond-only (2008-era) and under-builds tension even once wired. → §6 #3,#4.

4. **(high) No explicit substrate compliance spring (E_substrate)** — ligand is a fixed point at z=0; no series spring encoding E, so the Chan-Odde biphasic rigidity optimum (~1 kPa) and rigidity sensing cannot emerge, and there is no E-sweep driver. → Anchor each ligand through a harmonic spring whose k encodes E via a Hertz/Bangasser-Odde scaling (Magic-Number Block) + a 0.1–100 kPa sweep.

5. **(high) Force-dependent debranching absent** — zero `debranch` hits; branches only abort pre-nucleation, never rupture under load, so n_total is unbounded (the H5/γ design doc itself flags this). Pandit 2020: 0–2 pN junction force accelerates debranching >2 orders; ADP-Pi "young" vs ADP "old" junctions differ ~20×. → Model the Arp2/3-mother junction as a Bell-Evans slip bond whose off-rate rises with junction force (Pandit 2020) — adds the physics AND bounds n_total. Revise the "bond count grows monotonically" Sanity-Gate invariant via PI.

6. **(medium) No explicit membrane object** — no Helfrich/bending/area/tension; membrane lumped into the WAVE plane pin (lamellipodial load) + the ERM radial pin (cortex coupling). → Minimal deformable surface or tension/area reservoir the barbed ends push against (so KU-5.2 f-v emerges as a membrane reaction) + cortex coupled via ERM. κ_m ~20 kT, tension ~10–300 pN/μm.

7. **(medium) Nucleotide state (ATP/ADP-Pi/ADP) unmodeled** — no filament-age tag, so the young-vs-old debranching force ladder and cofilin age-selectivity can't be expressed (prerequisite for turnover/debranching to be quantitatively correct). The BAOAB tag-space extension hook (commit e860f31) is the natural integration point. Record as a scoped deferral.

8. **(low) No distinct stress-fiber contractile bundle** — the FA→interior load path runs through the cortex shell only; for single-cell-on-substrate, stress-fiber bundles anchored at FAs would give traction a direct interior path (Phase-1.5 fidelity item).

9. **(low, infrastructure) Syncthing working-tree corruption observed during the audit** — `.sync-conflict-*` artifacts + runtime files transiently read 0 bytes on disk while committed git blobs were intact. A production run launched against a mid-sync tree can silently import a truncated/0-byte module. → Add a pre-run integrity guard (assert each imported runtime module's on-disk size matches its committed git-blob size; pause launches while Syncthing reports active sync); clean up the `.sync-conflict-*` files in `outputs/h3/figs/` and `production/`.

---

## §6 FA-integration design synthesis (the restart deliverable)

The `fa_design` critic's full architecture for fixing KU-3.5 (cortical tension) and KU-5.1 (dendritic density). **Principle:** every link is an explicit HOOMD bond with mechanistic kinetics — no lumped Chan-Odde clutch at runtime; Chan-Odde's biphasic traction-vs-stiffness curve is the acceptance oracle only. The root cause is architectural: the cortex is a FLOATING SHELL with no substrate anchor.

**Full chain (build/wiring order behind `with_ecm` + `with_fa` in `build_cortex_full_simulation`):**

```
ECM substrate / fixed ligand layer at z=0  (immobile, position-restrained = mechanical ground)
        │  ← integrin-ligand CATCH bond (migrate Pereverzev→Kong 2009 two-state; F*~30 pN)
        ▼
  talin rod (13-domain serial WLC + per-domain two-state unfold/refold; R3 ~5 pN, +60–80 nm/unfold → exposes VBS)
        │  ← vinculin (discrete force-gated recruitment to exposed VBS; vinculin-actin DIRECTIONAL catch bond)
        ▼
  molecular clutch (cortex/stress-fiber actin → engaged integrin; load from myosin retrograde flow)
        ▼
  cortex + lamellipodium actin  ←ERM clutch→  membrane (Brownian-ratchet barbed-end load; tension ~10–300 pN/μm)
```

**The 11 design findings (ordered as the build sequence):**

1. **(critical) Substrate anchor at z=0 — the missing mechanical ground.** A layer of FIXED (zero-mobility / position-restrained) ligand particles at z=0 (HOOMD: exclude from the integration group, or pin with a stiff harmonic to initial position). Every FA nucleates from these. Without a fixed reference, myosin stress integrates to ~zero (the ~3,000× deficit). Anchor density per FA from KU-2.4 (50 integrins/FA). First step of the `with_fa` path.

2. **(critical) Integrin-ligand catch-bond clutch engaging cortex actin to the anchor.** Use the existing `IntegrinBondUpdater` (Pereverzev) PLUS a clutch bond from cortex/stress-fiber actin to the engaged integrin. Catch lifetime peaks at F*~30 pN (KU-2.5) → load-and-fail transmits cortical stress to the substrate. Migrate Pereverzev (current F*~7 pN, self-flagged) → Kong 2009 two-state (~30 pN). Chan-Odde = oracle.

3. **(high) Talin multi-domain WLC mechanosensor.** Explicit serial chain of two-state (folded/unfolded) WLC bonds between integrin and clutch/actin; per-domain Bell-Evans unfold k_u(F)=k_u0·exp(F·dx/kT) + refold k_f(F); R3 first ~5 pN; +60–80 nm contour per unfold (Yao 2016 / Tapia-Rojo 2020). Activate the existing (dead) talin config block. Creates the rigidity threshold; each unfold exposes a VBS. WLC/Tapia-Rojo = oracle.

4. **(high) Vinculin force-gated reinforcement (discrete, not a scalar multiplier).** Discrete vinculin particles bind force-exposed talin VBSs (gated by talin unfold state) and form parallel reinforcing bonds; model vinculin-actin as a **directional catch bond** (Huang 2017, Owen 2022). Activate the vinculin config block. **Replace** the forbidden `k_int_eff=k_int_bare·(1+α·N_vin)` scalar. FA reinforcement then emerges from the explicit population.

5. **(high) Molecular-clutch loading loop — the tension-building engine.** Myosin (Stam-Hocky, Hill, already wired) pulls cortex/stress-fiber actin into retrograde flow; engaged clutches (integrin-talin) resist → tension builds in the network = cortical tension. The parts exist but are uncoupled (no clutch bond between cortex actin and substrate). Closing the loop is integration of findings 1–3 + existing myosin, not new physics. Verify force balance: myosin contraction resisted by engaged clutches ⇒ network tension = γ. Chan-Odde = oracle.

6. **(high) Substrate compliance spring (E_substrate) for rigidity sensing.** Anchor each ligand to its fixed point THROUGH a harmonic spring whose k encodes E (Hertz/Bangasser-Odde derivable scaling, Magic-Number Block). Add a sweep driver over E~0.1–100 kPa. Validate traction-vs-E biphasic (optimum ~1 kPa) and retrograde flow monotonically decreasing with E, vs Chan-Odde / Bangasser-Odde oracles. Couple to talin so reinforcement-vs-E shows the threshold.

7. **(high → critical for KU-5.1) Explicit membrane object for lamellipodial barbed-end load.** Barbed ends carry no load → `exp(-Fδ/kT)` laws inert → capping maximal → density collapses. Introduce a deformable bounding surface or membrane-tension/area reservoir (Helfrich κ~20 kT, tension ~10–300 pN/μm) the barbed ends push against, so KU-5.2 f-v emerges as a membrane reaction. Couple cortex to it via ERM (finding 8). This is the KU-5.1 fix once FA anchoring provides cortex tension.

8. **(medium) ERM as a stochastic Bell-Evans slip clutch population (not a static pin).** Convert ERMHarmonic to an explicit force-sensitive SLIP clutch population (reuse `IntegrinBondUpdater`): per-linker k=1e-4 N/m (Alert 2015 per-linker value), force-dependent slip off-rate (Korkmazhan-Dunn 2022), density ~100/μm², gated rebinding. Keep the stiff element on constrained-BAOAB. SLIP (not integrin CATCH) kinetics.

9. **(medium) FA growth emergent from clutch population dynamics.** FA size/maturation emerges from the engaged-clutch population + reinforced state (KU-2.17 Hill threshold ~50 pN → N_engaged~10). `FAGrowthMonitor` becomes a feedback driver, not just an observer. Riveline 2001 / Balaban (constant stress ~5.5 nN/μm²) = oracles. Emerges once 1–5 are wired.

10. **(high) Integration & wiring into Cell.build with equilibration prelude.** Wire stepwise behind `with_fa`/`with_ecm`: substrate anchor → integrin/ligand particles + catch bonds → talin WLC → vinculin → clutch bonds to cortex actin → (membrane) → ERM clutch. Register all new bond/WCA params, extend `gamma_map`, replace `offsets['fa_integrin']=(cur,cur)` with the real tag range, and run `equilibrate_no_shear` as a mandatory prelude before the BAOAB production loop. Re-run KU-3.5/KU-5.1 only after the full chain is closed and equilibrated.

11. **(medium) Parameter sourcing + oracle separation.** Per-link Magic-Number Block from primary literature; every closed-form an oracle: integrin catch F*~30 pN (Kong 2009), talin 5–25 pN + contour 60–80 nm (Yao 2016/Tapia-Rojo 2020), vinculin-actin directional catch (Huang 2017), E~0.1–100 kPa (Chan-Odde sweep), FA constant stress ~5.5 nN/μm² (Balaban), Chan-Odde traction-stiffness biphasic optimum ~1 kPa as oracle.

**Co-requisites surfaced by the completeness critic** (must accompany the FA wiring, not follow it, for KU-3.5/KU-5.1 to be mechanistically reachable): **actin turnover** (§5 #1), **enclosed-volume pressure** (§5 #2), and the **membrane object** (#7 above). A turnover-free, pressure-free, membrane-free cell will still miss the gates even with the FA anchor wired.

---

## §7 gbook / KAIST institutional full-text fetch list

Deduplicated paywalled / institutional-access papers flagged by the lit agents (citation — DOI/PMID where known):

**FA / mechanotransduction (highest value for the restart):**
1. Kong et al. 2009, *J Cell Biol* — α5β1-FN catch bond, lifetime-vs-force, two-state — `10.1083/jcb.200810002` (PMID 19564406)
2. Chakrabarti, Hummer, Thirumalai 2017 — two-state allosteric catch-bond theory
3. del Rio et al. 2009, *Science* 323:638 — stretching single talin activates vinculin binding — `10.1126/science.1162912`
4. Yao et al. 2016, *Nat Commun* 7:11966 — talin rod-domain unfolding hierarchy — `10.1038/ncomms11966` (PMID 27384396)
5. Tapia-Rojo et al. 2020 — talin folding "tuning fork" of mechanotransduction (Sci Adv / PNAS)
6. Huang et al. 2017, *Science* 357:703 — vinculin directionally asymmetric (catch) bond with F-actin — `10.1126/science.aan2556`
7. Elosegui-Artola et al. 2016, *Nat Cell Biol* 18:540 — molecular-clutch force transmission vs rigidity (talin unfolding) — `10.1038/ncb3336`
8. Chan & Odde 2008, *Science* 322:1687 — clutch traction vs stiffness (ORACLE) — `10.1126/science.1163595`
9. Bangasser-Odde 2013 / Bangasser et al. 2017, *Nat Commun* 8:15313 — motor-clutch max force / optimal stiffness — `10.1038/ncomms15313`
10. Roca-Cusachs et al. 2009, *PNAS* — integrin force / α5β1
11. Owen et al. 2022 — vinculin-actin catch bond / reinforcement; Hu et al. 2024 (ForceChrono) — per-clutch dF/dt loading rate

**Cortex / membrane / tension / rounding:**
12. Chugh et al. 2017, *Nat Cell Biol* 19:689 — cortex thickness/turnover sets surface tension — `10.1038/ncb3525` (PMID 28530659)
13. Salbreux, Charras, Paluch 2012, *Trends Cell Biol* 22:536 — actin cortex mechanics — `10.1016/j.tcb.2012.07.001`
14. Tinevez et al. 2009, *PNAS* 106:18581 — cortical tension — `10.1073/pnas.0903353106`
15. Stewart et al. 2011, *Nature* 469:226 — hydrostatic pressure + cortex drive mitotic rounding — `10.1038/nature09642` (PMID 21196934)
16. Fischer-Friedrich et al. 2014, *Sci Rep* 4:6213 — surface tension + internal pressure of mitotic cells — `10.1038/srep06213` (PMID 25169063)
17. Saha et al. 2016 — cortex turnover/flow; Fritzsche et al. 2013/2014 — ERM/ezrin turnover dynamics
18. Korkmazhan & Dunn 2022, *Curr Biol* — ERM force-sensitive (slip) membrane-cortex clutch
19. Alert & Casademunt 2015, *PRL* 115:098101 — bleb nucleation / per-linker membrane-cortex stiffness & density (~100/μm²) — `10.1103/PhysRevLett.115.098101`
20. Ramanathan et al. 2015 — graded blebbistatin / cortex solid-vs-fluid response

**Cytoskeleton / network / myosin / branching:**
21. Broedersz & MacKintosh 2014, *Rev Mod Phys* 86:995 — semiflexible polymer networks (canonical) — `10.1103/RevModPhys.86.995`
22. Licup et al. 2015, *PNAS* — stress controls collagen mechanics (K∝σ^{3/2} 3D oracle) — `10.1073/pnas.1504258112`
23. Sharma et al. 2016, *Nat Phys* — strain-controlled rigidity transition / critical strain — `10.1038/nphys3628`
24. Jansen et al. 2018, *Biophys J* — architecture→G_0 / onset-strain — `10.1016/j.bpj.2018.07.020`; Lindström 2010 / Stein 2011 — 3D crosslinked-network provenance
25. Stam et al. 2015, *PNAS* 112:E4061 — filament rigidity/connectivity tune active-network deformation — `10.1073/pnas.1509663112`
26. Billington et al. 2013, *J Biol Chem* 288:33398 — NMII minifilament stoichiometry (~15 heads/side) — `10.1074/jbc.M113.499848`
27. Kovács et al. 2003, *J Biol Chem* — human cytoplasmic NMIIA kinetics / v0 — `10.1074/jbc.M204935200`
28. Erdmann & Schwarz 2012/2013, *PRL / Biophys J* — stochastic ensembles of non-processive motors (Hill emergence)
29. Ferrer et al. 2008, *PNAS* — single-molecule α-actinin/filamin off-rates + x_β — `10.1073/pnas.0710258105`
30. Ehrlicher et al. 2011, *Nature* — filamin force-dependent binding (catch) / FilGAP — `10.1038/nature10316`
31. Rognoni et al. 2012, *PNAS* — single-molecule force sensing of filamin (catch) — `10.1073/pnas.1117368109`; Gieseke/Schwaiger/Rief 2013 four-bead OT
32. Pandit, Pollard, De La Cruz 2020, *PNAS* — force-dependent debranching of actin networks — `10.1073/pnas.1911708117` (PMID 32461373)
33. Fäßler et al. 2020 / Ding et al. 2022 / Chou-Pollard 2022 — cryo-ET Arp2/3 branch junction, 70–71° (68±9°)
34. Bieling et al. 2016, *Cell* 164:115 — force feedback in self-assembling branched actin networks — `10.1016/j.cell.2015.11.057`
35. Li/Bieling/Mullins/Fletcher 2022, *eLife* 11:e73145 — force-dependent insertional capping, shared 2.7 nm gap — `10.7554/eLife.73145`
36. Demoulin et al. 2014, *PNAS* — actin ratcheting against a load; Footer et al. 2007, *PNAS* — single-filament stall force

**Cadherin (deferred, multi-cell scope):**
37. Rakshit 2012 / Manibog et al. 2014, *Nat Commun* 5:3941 — E-cadherin X-dimer catch — `10.1038/ncomms4941`

*(Open-access / already-resolvable anchors, no fetch needed: Bell 1978, Evans-Ritchie 1997, Dudko-Hummer-Szabo 2006, Hill 1938, Gittes 1993, Ott 1993, Mullins 1998, Schafer 1996, Pollard 1986, Pollard-Borisy 2003, Mogilner-Oster 2003, Leimkuhler-Matthews 2013/2016, Fixman 1974, Wachsstock 1994, Furuike 2001, Veigel 2002, Weeks-Chandler-Andersen 1971, Funk 2021 Nat Commun, Riveline 2001, Balaban 2001.)*

---

## §8 Prioritized action plan

Tags: **[doc]** documentation-only · **[code-free]** analysis/design with no behaviour change · **[needs-PI: frozen-integrator | param-change | gate-contract]** requires PI sign-off.

**Tier 0 — process / safety (do first, cheap):**
1. **[code-free]** Add a pre-run integrity guard to every production driver (assert on-disk module size matches committed git-blob; pause while Syncthing is mid-sync); clean up the `.sync-conflict-*` artifacts (§5 #9).

**Tier 1 — unblock the whole-cell gates (the restart, §6):**
2. **[code]** `with_ecm` builder path: Mikado substrate / fixed ligand layer at z=0 as position-restrained mechanical ground (§1 #2,#3; §6 #1).
3. **[code]** `with_fa` builder path: integrin catch bonds → clutch bonds → cortex/stress-fiber actin; register integrin_ligand + motor bond params + WCA cross-pairs; replace the zero-width `fa_integrin` tag slot; per-FA minifilaments; attach `IntegrinBondUpdater` + `FAGrowthMonitor` (§1 #1,#4,#6; §6 #2,#5).
4. **[code]** Explicit plasma-membrane object (Helfrich κ_m~20 kT, tension ~10–300 pN/μm) coupled to cortex via ERM and to lamellipodial barbed ends as Brownian-ratchet load — the KU-5.1/KU-5.6 fix (§5 #6; §6 #7).
5. **[code][needs-PI: frozen-integrator]** Call `equilibrate_no_shear` in Cell.build after ECM+FA instantiation; add a single global `dt=min(τ)` reconciliation + assert SHAKE convergence vs dt (§1 #11,#15; §6 #10).
6. **[code]** Add actin turnover (cofilin severing + treadmilling, ~5–30 s, Chugh 2017) to cortex AND lamellipodium — co-requisite for KU-3.5/KU-5.1 (§5 #1).
7. **[code]** Add enclosed-volume / osmotic-pressure term so ΔP=2γ/R emerges (Stewart 2011) — KU-3.1 (§5 #2).

**Tier 2 — complete the mechanotransduction loop:**
8. **[code]** Talin 13-domain WLC unfolding (Yao 2016, R3 ~5 pN) + discrete force-gated vinculin + vinculin-actin directional catch (Huang 2017); activate the dead talin/vinculin config; remove the scalar k_eff (§5 #3; §6 #3,#4).
9. **[code]** Substrate compliance spring (E_substrate) + 0.1–100 kPa sweep driver → biphasic Chan-Odde optimum (§5 #4; §6 #6).
10. **[code][needs-PI: param-change]** Migrate integrin (and later cadherin) catch Pereverzev→Kong 2009 two-state; pin F*~30 pN, lifetime, k_off0 (§3; §4 #5,#15; §6 #2,#11).
11. **[code]** Convert ERM static pin → Bell-Evans slip clutch population (Korkmazhan-Dunn 2022; Alert 2015 per-linker k) (§2; §4 #11; §6 #8).

**Tier 3 — mechanism-fidelity & magic-number corrections:**
12. **[code][needs-PI: gate-contract]** Replace the myosin **binned-r0 Hill ratchet** with physical sub-bin processive stepping; add an n_bins grid-invariance Sanity Gate (§2; §4 #9).
13. **[code][needs-PI: gate-contract]** Flip filamin slip→catch (two-pathway); revise the slip-sign Sanity-Gate test to a catch contract (do NOT edit inline) (§2; §3; §4 #4,#10).
14. **[code]** Add force-dependent Arp2/3 debranching (Pandit 2020 junction Bell-Evans) — bounds n_total; revise the "monotonic bond count" invariant via PI (§5 #5; §4 #14).
15. **[code][needs-PI: param-change]** Parameter corrections: capping δ_cap 0.3→2.7 nm + fix "0.3 pN" units mislabel (Li/Bieling 2022); myosin v0 1→~0.2 μm/s; α-actinin k_off0 1.0→~0.066/s (Ferrer 2008); branch t0 72°→70°; derive k_angle from 68±9° (§3; §4 #16,#9,#10,#14).
16. **[code]** Set LJ shift mode so U(2^(1/6)σ)=0 (true WCA) + force-continuity unit test; wire/consolidate the orphaned `build_excluded_volume`; derive σ from ×40 bundle diameter (§1 #12; §4 #18).
17. **[code]** Reconcile myosin heads/side 10→~15 (Billington 2013) with documented ×40 meso scaling (§3; §4 #9).
18. **[code]** Add per-bead nucleotide-state tag (ATP/ADP-Pi/ADP) via the BAOAB tag-space hook — prerequisite for age-aware turnover/debranching (§5 #7).

**Tier 4 — documentation / governance:**
19. **[doc][needs-PI: ratification]** Document the Hill fractional-step accumulator + the KU-3.5 Option C placement fix as ratified mechanism refinements; confirm no mean-velocity bias (§2).
20. **[doc]** Re-attribute every "Funk 2022" citation → Li/Bieling eLife 2022 (capping) / Furuike 2001→correct observable for filamin / Kovács for v0 (§3).
21. **[code-free][needs-PI: gate-contract]** Re-derive the D4 rigidity-percolation prefactor to justify/tighten the rebanded G_0 [1,50] Pa (§2). Surface that metaphase γ~1.6 mN/m exceeds the KU-3.5 ceiling (§5 #2).
22. **[doc]** Stand up a single gated-production entry point that runs STATIC + EMPIRICAL + RUNTIME tiers (incl. equipartition/FDT/grid-invariance) before any production run is accepted (lumped critic, sev low but enabling) (§2).

---

### Auditor notes on agent quality / contradictions

- **Subsumption:** the ECM audit's ~7 per-mechanism `high` rows are all symptoms of the single "H.1 not wired" gap (§1 row 2) — counted once.
- **Line-number drift, verified against the live tree:** the FA audit's `fa` slot is at `cell.py:638` and `with_fa` at `:584` (confirmed); the live tree has *more* dormant flags than any single audit named — but verification shows **only `with_fa` exists as a dormant flag**; there is **no `with_ecm`/`with_membrane`/`with_volume`/`with_junction` flag yet** (an earlier draft of this report mistakenly listed those — corrected here). The dormant artifacts that DO exist: `with_fa` (line 584, dead no-op), `fa` slot (638, never populated), `fa_integrin=0` (866), `offsets['fa_integrin']=(cur,cur)` (902). The `bridge/` files are `fa.py`, `fa_growth.py`, `integrin_bonds.py`, `motor.py` (no `fa_builder.py`/`clutch.py`/`talin.py` — talin/vinculin live only as disabled config consumed by `fa.py::ResolvedH4`). `junction/` is empty but for `__init__.py`.
- **No material contradictions** between audits, lits, and critics. The param_verify critic adversarially re-derived analytic landmarks (e.g. Pereverzev F*=6.992 pN recomputed) and its benign-vs-real triage is internally consistent with the lit dossiers' `agree` fields. Two lit/critic agreements worth flagging as *non-obvious wins*: (a) the k_ERM 1000× "softening" accidentally lands on the correct per-linker value (Alert 2015) — provenance fix, not a physics error; (b) the Pereverzev F*~7 pN mismatch is honestly self-flagged in-config and PI-deferred, with the gate written against the analytic value.
