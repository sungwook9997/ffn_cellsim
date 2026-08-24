---
archived_on: 2026-07-28
superseded_by: aleph/docs/v2_audit/AC_EXECUTION_PLAN_2026-07-25.md
reason: >
  Written BEFORE the 2026-07-25 PI reframe, i.e. for a different objective — forward prediction
  and magnitude matching, rather than inferring per-cell-type parameters with the gate on
  mechanical connectedness. Archived, not deleted: its measurements and reasoning stand as a
  record of what was true then. Nothing in it may be quoted as current state; STATE.md is that.
  Selected mechanically: pre-reframe AND cited by no live file (code, STATE.md, CLAUDE.md,
  cell_engine/, gate_contracts/, tests, Makefile). Citations from run outputs and from other
  pre-reframe documents were not treated as protective.
---

# Compartment Build Plan — FF single-cell full-fidelity program (2026-07-07)

**Origin:** PI 2026-07-07 "전부 다 구현 계획 세우고 일단 전부 구현하자" — after the compartment audit
(36-agent workflow: 6 present / 16 partial / 10 absent). Two workflows produced this: `compartment-audit`
(gap-finding, adversarial rg-verified) → `compartment-buildplan` (20 per-compartment designs → this master
plan). Correctly-absent multicellular units (desmosome/hemidesmosome — MCF7 lacks α6β4 KB-PIV-5; tight
junctions — DCM territory) are excluded by design.

**Execution rule:** phase-by-phase, each unit closes on its analytic/lit VERIFICATION GATE (a number, not a
vibe). Register-first units do NOT get a coded constant until SourceEvidence→KnowledgeClaim lands (PI-gated).
BLOCKED items (§3) surface to PI before build — never silently defaulted.

---

# MASTER BUILD PLAN — FF Single-Cell Compartment Program (20 units)

Scope: 20 compartment designs for the going-forward FF engine (`aleph/laws/`, Warp, µm·pN·s). This plan sequences them by dependency and ROI, isolates the shared solver/infra threads that gate multiple units, and surfaces every KB/PI block that must clear before a unit can carry a number. Nothing here invents a constant; register-first units stay gated.

---

## 1. Dependency graph (prose)

### Live substrate this program builds on (already on disk, not re-built)
FF cortex/SF network (`network_warp.py`, `cortex_assembly.py`, native N=70,686), FA clutch (`fa_clutch_warp.py`, `fa_anchor.py`), Hand-KMC (`hand_kmc.py`), weave/architecture (`weave.py`, `architecture_spec.py`), nucleus **bead-cloud** (`nucleus_shell_kernel`, `_seed_nucleus_cloud`), **independent** plasma-membrane sheet (`membrane_surface.py`, commit feac927 — exists as a compartment but is NOT yet the production surface; the driver still runs a lumped cortex-hull ΔP proxy), Mikado ECM (`ecm_mikado.py`), implicit GPU solver (`implicit_ff.py`).

### Five cross-cutting infra threads (each unblocks several compartments — build once, reuse)
- **Thread A — implicit anchor DOFs (Winkler/Hookean movable anchor in `ff_implicit_step_gpu`).** Gated by the substrate unit; the rigid-limit byte-identity gate is the guardrail. Prereq for: **Substrate (1)** and conceptually shared with **Fibrillar-ECM (14)** (both make the clutch reaction land on a real, movable DOF instead of an E=∞ pin).
- **Thread B — implicit nonlinear per-bond tangent** (K assembled from `k_tan(λ)`, not constant-k only). Prereq for: **IF cage (5)**, **Nuclear envelope (6)** strain-stiffening lamina.
- **Thread C — implicit multi-node-set / off-diagonal coupling block** (a second node-set enters K with cross-block rank-1 terms). Prereq for: **Nuclear lamina (7)** (cortex↔lamina LINC), **LINC (9)** (nucleus-bead DOFs become columns in K).
- **Thread D — promote the independent membrane sheet (Template-2 mesh) into the production driver** (real membrane node block + area kernels + containment + ERM index map, replacing the cortex-hull ΔP proxy). Prereq for: **ERM (10)**, **Membrane reservoir/caveolae (11)**, **Spectrin (19)**, **Septin (18)** curvature target, and the *local* per-face tension form of **Piezo (15)**.
- **Thread E — per-node mobility/drag array** (replace uniform `dt_mu` with a µ-array). Prereq for: **Nucleolus (20)** viscous clock, **Poroelastic interior (8)** drained/undrained beads.

### Compartment-to-compartment edges
- **Substrate-E (1) → FA-maturation (3) → Fibrillar-ECM (14)** is the PI experimental spine: a compliant substrate gives talin a real resting load; matured clutches carry the traction that remodels the ECM. (3) and (14) can each run on a synthetic/fixed anchor for their own gates, so the hard edge is conceptual ROI ordering, not a build block.
- **Myosin-Hill (2)** is independent (cortex only) but is a *loading source* for (1)/(3)/(14) durotaxis/maturation gates — if run static, those units drive anchors kinematically at v_u=120 nm/s instead.
- **Nuclear stack ordering: NE-surface reuse (6) → Lamina (7) → Poroelastic interior (8); IF cage (5) and LINC (9) attach once a nucleus shell exists; Nucleolus (20) is a nested leaf of the nucleus.** (6) reformulates the bead-cloud into a real surface DOF reusing the membrane kernels; (7) refines the shell with the lamin-A/B split; (8) fills the interior it pushes against; (5) cages it and LINC-anchors to the cloud; (9) couples nucleus↔cortex; (20) seeds droplets inside it.
- **Membrane stack: Thread-D → ERM (10) + Reservoir/caveolae (11) + Spectrin (19)** all ride the promoted sheet. **Piezo (15)** reporter reads membrane tension (compartment already exists) but its *local* form wants Thread-D too.
- **Microtubule aster (4)**: static stiff-aster/buckling/L_p gates are **standalone** (reuse live bending+constraint+implicit). Containment-strut gate needs the membrane (live). Nucleus-**centering** gate needs an **FF nucleus that does not exist yet** → deferred behind the nuclear stack.
- **Stress-fiber array (13)**: hard-blocked on a **spread/adherent ventral morphology** (the FF cortex is a sphere). Interim analytic flattened-cap unblocks the gate; final number needs true spread shape.
- **KB-blocked minors — Fascin (16), Glycocalyx (17), Septin (18), Nucleolus (20)** — depend only on their host (membrane/cortex/nucleus, all live) plus their own literature registration; they are last purely because the KB rows do not exist.

---

## 2. Phased sequence

Effort key: S≈0.5, M≈1, L≈2 units. "Register" = SourceEvidence→KnowledgeClaim in Notion, PI-gated, CrossRef-verified, verdict=OK before any constant is coded.

### Phase 1 — PI substrate axis + hard-rule fixes + a pure-port quick win  (~5.5 units)
Highest ROI: the PI's stiffness/ligand experimental axis, the two architectural-rule violations (rigid-pin clutch, constant myosin), and a standalone port that yields analytic gates immediately.

| Unit | Effort | Register-first? | Lands |
|---|---|---|---|
| **(2) Myosin-II Hill FV + Stam-Hocky minifilament** | L | Yes — NMII hyperbolic Hill law (Hill 1938 / Kovács 2003 a/F_s≈0.5); reconcile per-head stall 0.5 vs 2 pN *in the claim* | Emergent γ from explicit duty×Hill minifilaments, replacing the swept constant `f_myo` (hard-rule fix) |
| **(1) Compliant substrate + ligand coat + Thread-A** | L | Partial — ρ_L Bare/Pre/Lam4 areal density (magnitude only); confirm KB-PIV-5/6 verdict OK | Movable Winkler substrate anchor; ligand presets (LAMININ_A6B1, COLLAGEN_A2B1); durotaxis re-seeding |
| **(4) Microtubule aster + MTOC (static MVP)** | M | No for MVP (KAPPA_MT live); DI/centering deferred | Radial aster, Euler-buckling + L_p gates on the live bending/constraint/implicit path |

**Gates that close Phase 1:** (2) Tier-1 Hill-fit recovers a/F_s=0.50±0.02 and f_s=n·F_stall; Tier-2 per-minifil 17–100 pN; Tier-3 γ_active vs KB-3.5 ±30% (may honestly REFUTE → report the deviation factor, do not tune). (1) rigid-limit byte-identity max|Δ|<1e-14, Bangasser-Odde biphasic within KB 2× band + saturating plateau ≳5 kPa, durotaxis sign net-positive across ≥8 seeds, ligand ordering Bare≪Lam4≲Pre. (4) F_crit=π²κ/L²=7.90 pN within 2%, L_p≈4.7 mm in 1–8 mm band, GPU↔CPU <1e-12.

### Phase 2 — substrate downstream, reciprocity, cheap wins  (~3.5 units)
Everything that hangs off Phase-1 substrate/clutch work or is anchored wiring.

| Unit | Effort | Register-first? | Lands |
|---|---|---|---|
| **(14) Fibrillar adhesions remodeling ECM** | M | No (KB-2.5 anchored) — but reconcile α5β1-FN vs Mikado α2β1 ligand mismatch; carry F*≈7 vs 30 pN flag | Reciprocal clutch (Newton-3rd onto ECM node), hash-grid re-bind, emergent ECM remodeling |
| **(3) FA maturation: talin + vinculin + Hill growth** | M | Yes — vinculin k_diss off-rate; reconcile Han2021/TapiaRojo2019 DOIs | Per-clutch k_int_eff array, emergent maturation on the resting clutch load |
| **(15) Piezo1 reporter** | S | Reporter no (KB-3.10 ready); **feedback stub default-OFF** until density/conductance/Ca-gain registered | Tension-gated P_open field, area-weighted open fraction |
| **(12) Cytosol viscoelastic (G′-memory + Darcy)** | L (start) | Core no (KB-3.B3.1/.2); **crowding piece stubbed off** until Ellis/Delarue registered | SLS memory tether via diag_extra + implicit Darcy pore-pressure |

**Gates:** (14) Newton-3rd Σ-force <1e-8 rel, 1/r mid-field decay n∈[0.5,1.5], catch-slip peak F*, compliance monotonicity. (3) talin Bell parity <1e-10, Hill fixed point crosses at F=5.0±0.2 pN, vinculin reinforcement 1.1–2× + τ_ratio>1. (15) six analytic logistic gates incl. P_open(rest)=0.0347 (the pN/nm↔pN/µm unit-trap guard). (12) SLS stress-relaxation G′∈[33,79] Pa & τ∈[0.5,3.5] s, grid-invariance <5%, Darcy Dp∈[40,60] µm²/s, k_mem=0 bit-identity.

### Phase 3 — nuclear stack (register-heavy, Threads B & C & E)  (~7 units)
Ordered NE→lamina→interior; IF cage and LINC attach in parallel once the shell exists.

| Unit | Effort | Register-first? | Lands |
|---|---|---|---|
| **(6) Nuclear envelope surface DOF + Thread-B** | L | Yes — K_A,NE, σ_NE0, κ_NE, ρ_NPC, LINC force; **ε_rupt/σ_rupt may be irreducible → PI proxy** | Bead-cloud → real surface shell (area+bending+turgor), emergent rupture mask |
| **(7) Nuclear lamina (lamin A/B split) + Thread-C** | L | Yes — **Moore2016 R_nuc=0.77 + Kim2020 E_nuc≈399 Pa**, PI ratify in-situ vs isolated | Strain-stiffening bilinear shell, cortex↔lamina LINC off-diagonal |
| **(8) Chromatin/nucleoplasm poroelastic + Thread-E** | L | Yes — Dp_nuc, Skempton B, ξ_nuc, permeability (cytoplasm Dp FORBIDDEN) | Darcy pore-pressure interior, drained/undrained relaxation |
| **(5) Intermediate-filament cage + Thread-B** | L | Yes — Kreplak/Block **digitized (λ,f) curve**; HALT if only qualitative | Tabulated nonlinear backbone (zero free params), rupture, nucleus protection |
| **(9) LINC nesprin-SUN coupling + Thread-C** | M | Yes — **entire LINC lit absent** (Déjardin 8 pN, spectrin WLC Lp proxy PI-gate) | WLC nucleus↔cortex bonds; k_linc unset until registered |
| **(20) Nucleolus / LLPS inclusions** | M | Yes — σ_no, η_no (Caragine 2018 / Feric 2016) | Young-Laplace droplets, per-node viscous mobility |
| **(4b) MT nucleus-centering** | S | DI/centering KB rows (Gittes/Walker/Dmitrieff) | MTOC anchors to the now-existing FF nucleus |

**Gates:** (6) Young-Laplace ΔP·R/2 <1%, bulk E_nuc∈1–10 kPa, rupture EMERGES at <10% cross-section & zero at ≥50%, Helfrich κ <5%. (7) round-trip E_eff within 10% of 399 Pa + grid-invariant, strain-stiffening knee at 3µm, φ_A^0.5 scaling. (8) Fickian σ²=6Dp·t within 10%, Hu-2010 poroelastic master-curve collapse, τ_p in AFM band. (5) small-strain tangent 942 pN/µm ±1%, K(1.5)/K(1) = digitized ratio, force-free at t=0, NE-strain drops cage-ON. (9) WLC host↔device <1e-10, momentum closure <1e-12, resting nesprin tension 8 pN (band 6–10). (20) τ=η_no·R/σ_no within 0.5–2×, grid-invariant.

### Phase 4 — membrane-associated stack (Thread-D)  (~3.5 units)
Promote the independent sheet into production first (the bulk of the cost + a parity-regression risk — keep membrane-off bit-identical), then the linkers/reservoir/skeleton.

| Unit | Effort | Register-first? | Lands |
|---|---|---|---|
| **Thread-D wiring** (prereq sub-task) | M | No | Real membrane node block in `simulate_compressed_shell_on_device`, ERM index map, containment |
| **(10) ERM catch/slip Hand-KMC** | L | Yes — ERM k0, x_β, areal density (Fritzsche FRAP / ezrin spectroscopy); catch-vs-slip open | Force-free bound tethers, Bell/Pereverzev KMC, emergent bleb |
| **(11) Membrane reservoir + caveolae** | L | Yes — caveolar area fraction, σ_flat, kinetics; **delete magic RESERVOIR_STRAIN=0.60** | A0_eff(t) reservoir, caveolae flattening entity |
| **(19) Spectrin-ankyrin skeleton** | M | Yes — μ0 (RBC proxy, PI-gated provenance tag) | 2D WLC triangulated shear network |

**Gates:** Thread-D Laplace closure + membrane-off bit-identity. (10) single-bond Bell within 5%, discrete→continuum γ_MCA=10 pN/µm in band. (11) reservoir plateau length = seeded inventory <1% + post-plateau slope=K_A ±2%, caveolae buffer attenuation, tether f_t 5–40 pN. (19) virtual-shear μ0 matches Dao-2006 closed form <5% + 6–9 µN/m RBC band, force-free build, whole-cell shift <2% (the "minor" falsifiable claim).

### Phase 5 — spread-morphology-dependent + KB-blocked fidelity minors  (~4.5 units)
Lowest ROI / most gated. Build once their prerequisite geometry or literature lands.

| Unit | Effort | Register-first? | Lands |
|---|---|---|---|
| **(13) Stress-fiber array (ventral/dorsal/transverse)** | M | Register only if arc taxonomy promoted to gate; reconcile F*≈7 vs 30 pN | 3-population SF array on a spread ventral footprint; the H.7 net-traction→0 sign fix. **Blocked on spread morphology** (interim flattened-cap unblocks the gate) |
| **(16) Filopodial fascin + tip clutch** | M | Yes — fascin kinetics (Courson/Aratyn) + **fascin single-molecule link_k genuine missing datum**; new sub-id (KB-3.8 is occupied) | Fascin Bell-slip crosslink, tip integrin clutch; N²-tight Euler-buckling discriminator |
| **(17) Glycocalyx brush** | M | Yes — L_g + **grafting spacing s** (may be absent for MCF7 → HALT, do not back-solve) | Alexander-de Gennes steric coat, emergent FA-maturation priming |
| **(18) Septin filaments** | L | Yes — κ_sep, density, k_off0, c*, σ_c; CrossRef-verify authorship | Curvature-sensing attach + cortical-rigidity crosslink; confined-migration payoff needs a channel geometry (absent) |

**Gates:** (13) net |Σtraction|→0 with CI narrowing as N_FA↑, inward-radial>0 per seed, per-cell 10–100 nN, force closure T≈5–6 nN/fiber <1%. (16) F_c within 15% of N²-tight Euler (~7.7 pN) and ~N× stiffer than the loose stand-in, fascin off-rate host↔GPU parity. (17) AdG P(D) analytic ≤1e-12, standoff 10–100 nm, catch-lifetime rises brush-ON. (18) curvature-preference peak at R*∈[1,3] µm with ≥2× enrichment, κ_eff beam-theory increase.

---

## 3. BLOCKED — needs PI decision or KB registration before build

Surface all of these; none may be silently defaulted (physiological-baseline + no-magic-number rules).

**Hard KB blocks (entire literature absent — build cannot start the numeric path):**
- **(9) LINC** — no nesprin/SUN/LINC claim exists at all. Need Déjardin 2020 (8 pN resting), Autore/Rief spectrin WLC L_p (proxy from erythroid spectrin — **PI must ratify the proxy**), Crisp gap/density. **k_linc stays unset.**
- **(17) Glycocalyx** — no claim. Need Paszek/Shurer L_g and especially **grafting spacing s** (sets the AdG prefactor); s may not exist for MCF7 → **HALT rather than back-solve s to a standoff.**
- **(18) Septin** — no claim, no code. Need κ_sep/L_p, cortical density/stoichiometry, k_off0, and curvature preference c*/σ_c. σ_c must be *derived* from the amphipathic-helix model, not fitted. Confined-migration gate also needs a channel-confinement geometry that does not exist.
- **(16) Fascin** — kinetics absent; **fascin single-molecule stiffness link_k has no AFM anchor** (unlike α-actinin/filamin) → derive from bundle bending or HALT; do **not** inherit filamin 8.2e5. Note KB-3.8 is already "filopodium structure" — register a NEW sub-id.

**Register-first (anchors partial; specific data missing):**
- **(5) IF cage** — the nonlinear law must be the **digitized Kreplak/Block (λ,f) curve**. If the papers give only qualitative curves with no extractable points, **HALT** — a non-magic nonlinear bond cannot be built.
- **(6) Nuclear envelope** — K_A,NE, σ_NE0, κ_NE, ρ_NPC, LINC force; the **rupture threshold ε_rupt/σ_rupt may be genuinely irreducible to a single number** (Denais/Raab report it as curvature+lamin-context) → PI-approved documented proxy, like the α6β1 precedent.
- **(7) Nuclear lamina** — **R_nuc/E_nuc audit dispute**: Moore2016 R_nuc/R_cell≈0.77 and Kim2020 in-situ E_nuc≈399 Pa are flagged by audit#18/#19 but **not yet registered**; the in-situ-vs-isolated (399 Pa vs 1–5 kPa) choice is a **measurement-protocol PI decision** (physiological-baseline favors in-situ). LINC rupture force likely absent.
- **(8) Poroelastic interior** — Dp_nuc, Skempton B/ν_u, ξ_nuc, permeability all absent; **borrowing cytoplasm Dp (KB-3.B3.2) is a forbidden magic number** (source says not nucleus-specific).
- **(2) Myosin-Hill** — the hyperbolic Hill FV law is not yet a KnowledgeClaim; **duty ratio has no MCF7 value** (lit ~0.1–0.3); **per-head stall 0.5 pN (Kovács single-head) vs 2 pN (KB-3.18) conflict** must be reconciled *in the claim by lit*, not by which passes 0.1–1 mN/m.
- **(11) Membrane reservoir/caveolae** — caveolar area fraction, σ_flat, flattening kinetics absent (Raucher-Sheetz itself unregistered). **Delete the magic RESERVOIR_STRAIN=0.60**; if inventory numbers can't be sourced for MCF7, surface to PI.
- **(10) ERM** — k0, x_β, areal density absent; catch-vs-slip itself open. Also fix the existing bug conflating tube-pull f_t with single-bond rupture.
- **(20) Nucleolus** — σ_no, η_no (Caragine 2018 / Feric 2016) must be extracted from PDFs + DOI-verified; no MCF7 value likely → HALT rather than default.
- **(19) Spectrin** — no nucleated-cell μ0; **RBC 6.3 µN/m is a PI-gated provenance-tagged proxy** (bounded because it moves whole-cell mechanics <2%).
- **(1) Substrate** — ρ_L Bare/Pre/Lam4 areal densities absent → **only the ordering gate (4), not a magnitude match, is claimable**; confirm KB-PIV-5/6 verdict=OK first.
- **(12) Cytosol** — the **crowding piece is hard-blocked** (colloid-osmotic φ, Π(φ), Ellis/Delarue) → ship default-OFF; G′/Darcy core proceeds.
- **(3) FA maturation** — vinculin **k_diss off-rate not in KB**; Han2021/TapiaRojo2019 are NO_DOI_FOUND → reconcile DOI before a cited deliverable.
- **(4) MT** — DI rates (Walker, keep default-OFF, in-vivo differs) + Dmitrieff centering anchor must be registered before those gates; **FF nucleus absent** → centering is proxy-only until Phase 3.
- **(15) Piezo** — **feedback stub ships default-OFF, gain=0**; density/conductance/Ca-gain must be registered before any nonzero coupling.

**Engineering (not KB) blocks:**
- **(13) Stress-fibers** — needs a **spread/adherent ventral morphology**; the FF cortex is a sphere. Interim analytic flattened-cap + FA ring unblocks the gate; the final traction number needs the true spread shape.
- **Threads A–E** are prerequisite solver/driver changes, each with a parity guardrail (rigid-limit byte-identity, membrane-off bit-identity) protecting the validated cortex path.

**Recurring flags to carry (not blocks, but never retune to a gate):** α5β1 F*≈7 pN vs Kong-2009 ~30 pN (surface if traction floors ~4× low — same force-magnitude-deficit family as the γ-floor/SF-force HALTs); R_FA=1.5 µm NeedsSource/InternalConsistency; NMII areal density Nie-2015 HeLa 0.625/µm² (not MCF7); Chan-Odde saturating-vs-peaked carried PROVISIONAL, diagnostic-only.

---

## 4. Definition of done (whole program)

The FF single-cell is a fully mechanistic, GPU-native (Warp, µm·pN·s) cell in which every compartment — compliant ligand-coated substrate, Hill/Stam-Hocky myosin minifilaments, maturing talin-vinculin focal adhesions remodeling a deformable ECM, microtubule aster, intermediate-filament cage, a real nuclear-envelope surface with strain-stiffening lamina, poroelastic chromatin interior, LINC coupling, nucleolar LLPS inclusions, ERM-tethered plasma membrane with an emergent reservoir/caveolae buffer, spectrin skeleton, viscoelastic-poroelastic cytosol, Piezo tension reporter, stress-fiber array, glycocalyx, and septins — is an **explicit particle/bond/kernel system initialized at its physiological operating point**, with **zero lumped proxies and zero constants tuned to pass a gate**. Every quantitative constant is KB-derivable with verdict=OK provenance (or the unit is explicitly gated OFF and surfaced to PI as a missing datum); every unit passes its analytic ground-truth gate (parity, force-balance, buckling, Laplace, relaxation, or grid-invariance) as the primary always-on test, with literature oracles as cross-checks only; each landed milestone emits an interactive HTML morphology viewer plus a REPORT.md Figures section; and all shared solver threads (anchor DOFs, nonlinear tangent, multi-node-set coupling, membrane-in-production, per-node mobility) preserve byte/bit-identity to the pre-existing validated cortex baseline. Where a headline number still refutes (e.g. the ~530× γ-floor persisting under emergent Hill myosin), the deliverable is the **honest mechanistic measurement of the deviation factor**, reported not tuned away.

---

## Appendix — per-compartment design summary (approach · effort · gate · KB)

- **FF Compliant Elastic Substrate + Ligand-Density Coat (E_sub + Bare/Pre/Lam4) — Bangasser-Odde molecular clutch on a real, deformable substrate** — `build-new-KB-anchored` · effort L · deps: FF cortex/SF network (ff/fiber_network.py + ff/cortex_assembly.py) must exist to
  - gate: FOUR concrete gates, all numbers not vibes: (1) RIGID-LIMIT off-path invariance — E_sub→∞ ⇒ k_sub→∞ ⇒ anchor immovable ⇒ force field byte-identical (max|Δ|<1e-1
  - KB: STIFFNESS AXIS FULLY ANCHORED — port-ready, no registration needed: KB-1.5 (E_sub biological range 0.1–100 kPa, PAA 0.1–50 kPa, Phase-1 default 5 kPa, ν=0.45), 
- **Myosin II Hill force-velocity + explicit Stam-Hocky multi-head bipolar minifilament (FF cortex)** — `register-KB-then-build` · effort L · deps: HARD PREREQ: FF cortex actin network (ff/network_warp.py CrosslinkedCortex + gam
  - gate: Tier-1 (ANALYTIC ground truth = the gate number): isolate ONE minifilament between two actin anchors, sweep imposed load → measured force-velocity curve must FI
  - KB: PARTIAL — structural anchors exist, the force-velocity law must be registered first. HAVE: KB-3.18/KB-3.1 (minifilament ~30 myosin, ~300nm bare zone, per-head ~
- **Focal-adhesion maturation: talin unfolding + vinculin reinforcement + Hill FA growth (FF)** — `fine-grain-existing` · effort M · deps: SAME compartment: ff/fa_clutch_warp.py integrin clutch (EXISTS — hard prereq, is
  - gate: THREE analytic ground-truths (numbers, not vibes): (1) TALIN Bell parity — kernel k_unfold(5pN) must equal host closed-form k_u0·exp(F·dx/kBT); with k_u0=3e-4 s
  - KB: ALL anchors present with numbers. KB-2.6 talin: k_unfold(F)=k_u0·exp(F·dx/kBT); one-state=R3 F_th~5pN, k_u0~3e-4 s^-1, dx~0.5nm (13-domain R1-R13 heterogeneity 
- **Microtubule aster + MTOC (FF, Warp) — centrosomal radial array, compression strut, dynamic instability** — `port-from-archive` · effort M · deps: STANDALONE for the static stiff-aster + Euler-buckling + L_p gates (reuses only 
  - gate: Three analytic ground-truths (numbers, not vibes): (1) EULER BUCKLING (primary, NF2007 p9 = the bending kernel's own stated anchor): a clamped aster arm under a
  - KB: CORE mechanics anchor ALREADY LIVE: ff/units.py KAPPA_MT=20.0 pN·µm² grounded to Nédélec-Foethke 2007 NJP 9:427 p9, cross-checks Gittes 1993 EI=2.2e-23 N·m²=22 
- **Intermediate-filament (keratin/vimentin) perinuclear cage — nonlinear strain-stiffening, FF/Warp** — `register-KB-then-build` · effort L · deps: NUCLEUS compartment — EXISTS in FF (ResolvedNucleus / _seed_nucleus_cloud / nucl
  - gate: (1) SMALL-STRAIN TANGENT: df/dx|_{λ→1} = k_bb = 942 pN/µm (analytic: E_if=6e6 Pa · A_if=π(5nm)²=7.854e-17 m² / l_seg=0.5µm = 9.42e-4 N/m → 942 pN/µm in FF units
  - KB: MUST register SourceEvidence→KnowledgeClaim FIRST — no IF claim exists in knowledge_claim (verified: keratin/vimentin/Kreplak/Block/Mucke/Lichtenstern all absen
- **Nuclear envelope (double bilayer + lamina + NPC) as a distinct FF surface DOF, with emergent rupture** — `register-KB-then-build` · effort L · deps: Cortex (H.3) must land first — LINC tethers attach NE nodes to cortex nodes and 
  - gate: FOUR analytic/lit gates (numbers, not vibes): (1) LAPLACE ground truth (always-on, per 'oracle is cross-check' memo): resting NE with chromatin turgor ΔP_nuc ba
  - KB: MUST register SourceEvidence→KnowledgeClaim first — envelope-surface constitutive constants are genuinely absent (BM25 corpus probe returns only review text, no
- **Nuclear lamina (lamin A/C + B meshwork, GPU-native FF shell)** — `register-KB-then-build` · effort L · deps: (1) KB registration MUST land first — Moore2016 R_nuc/R_cell≈0.77 + Kim2020 in-s
  - gate: PRIMARY analytic round-trip (a number): build lamina from E_nuc=399 Pa, apply a known uniform compressive strain ε (δ=ε·R_nuc), sum restoring force, recover eff
  - KB: PARTIAL — physics-law anchors exist: KB-3.B2.1 (E_nuc band 1-10kPa in-situ~5k/isolated~8k, cyto 0.5-1k, N:C 1.4-5x in-situ, "do NOT hard-code 10x"; Caille2002),
- **Chromatin / nucleoplasm poroelastic-viscoelastic nuclear interior (FF)** — `register-KB-then-build` · effort L · deps: HARD PREREQ: KB registration of Dp_nuc + undrained/drained ratio (SE->KC, PI-gat
  - gate: Two grid-invariant analytic ground-truths. GATE-1 (diffusion kernel alone): seed a Gaussian pore-pressure blob, measure variance growth — must obey sigma^2(t)=s
  - KB: HAVE (elastic/viscous skeleton): KB-3.B2.1 E_nuc 1-10 kPa (in-situ ~5, isolated ~8) = DRAINED chromatin modulus anchor; KB-3.B2.2 nuclear viscosity~[lamin-A]^3,
- **LINC complex + perinuclear coupling (nesprin-SUN cortex↔nucleus springs), FF/Warp** — `register-KB-then-build` · effort M · deps: Nucleus compartment must be ON — SATISFIED (already in FF: resolve_nucleus + _se
  - gate: Three concrete numbers, all pass-before-commit: (1) DEVICE-VS-HOST WLC parity — linc_wlc_force on cuda:0 reproduces the host Marko-Siggia value to relative erro
  - KB: HARD BLOCK — the entire LINC/nesprin literature is ABSENT from the KB. source_evidence returns 0 rows for dejardin/arsenovic/autore/lombardi/crisp/nesprin/linc;
- **ERM membrane-cortex linkers — production catch/slip Hand-KMC** — `register-KB-then-build` · effort L · deps: HARD PREREQ: independent plasma-membrane SHEET promoted into the production loop
  - gate: TWO analytic numbers. GATE-1 (single-bond Bell): under constant pull F, KMC survival exponential with rate p0·exp(F·x_β/kBT) — measured off-rate within 5% of an
  - KB: PARTIAL — register first. HAVE: k_ERM stiffness 0.1 N/m = 1e5 pN/µm (KB-3.1, KB-3.18, KB-3.23); membrane tube-tether force f_t 5-40 pN + γ_MCA 1e-6–1e-4 J/m² + 
- **Membrane reservoir engaged in production + caveolae flattening entity (FF)** — `register-KB-then-build` · effort L · deps: CORTEX must land first (present) — ERM tethers + containment couple membrane nod
  - gate: Four numeric gates. (1) LAPLACE closure: closed membrane sphere R, uniform σ → summed nodal tension force = 0 to <1e-6·(σR); inferred inward ΔP = 2σ/R to <1%. (
  - KB: HALT-gated. KB has membrane anchors KB-3.B1.1 (γ_mem 3-40 pN/µm bilayer band), KB-3.B1.2 (κ_m 0.0828 pN·µm), KB-3.B1.3 (K_A 2.35e5 pN/µm = 0.235 N/m Rawicz 2000
- **FF cytosol: viscoelastic G′-memory (SLS/Zener) + poroelastic Darcy pore-pressure + crowding** — `build-new-KB-anchored` · effort L · deps: NONE hard-blocking — the cytosol is a background medium touching every FF node, 
  - gate: FOUR concrete numeric gates. (1) STEP-STRAIN STRESS RELAXATION vs analytic SLS: apply shear strain γ₀, record σ(t); SLS analytic ground truth is single-exponent
  - KB: CORE anchored, ONE sub-piece must register-first. G′ storage: KB-3.B3.1 (Dessard 2024, MRS Burgers fit) MCF7 G=32.9±6.0 Pa τ=3.23s / MDA G=38.6 τ=0.57s / MCF10A
- **Stress-fiber array (ventral/dorsal/transverse arcs) assembled into the spread cell** — `wire-existing-pieces` · effort M · deps: HARD: a SPREAD/adherent cell ventral footprint (flattened basal disk with a FA s
  - gate: SIGN RESOLUTION (the H.7 fix): the −127±156 non-convergence was measuring the NET traction VECTOR — which for a radially-symmetric spread cell MUST → 0 by force
  - KB: Anchors AVAILABLE — no new registration needed for the core: KB-2.12 (per-cell traction 10-100 nN, per-FA 1-10 nN typ 5, per-clutch 5-20 pN, stress 100Pa-10kPa)
- **Fibrillar adhesions pulling/remodeling the deformable ECM (reciprocal FA traction)** — `wire-existing-pieces` · effort M · deps: ecm_mikado.py (Mikado collagen ECM) — DONE. Combined cortex+ECM single-array ass
  - gate: FOUR gates, analytic/lit numbers not vibes: (1) Newton's-3rd-law force balance — Σ clutch force on cortex nodes = −Σ clutch force on ECM nodes to machine precis
  - KB: Anchors AVAILABLE — no new registration needed for the core: KB-2.5 (Kong2009, α5β1-FN catch-slip: k_slip=0.5/s Fs=30pN, k_catch=0.4/s Fc=7pN — already encoded 
- **Piezo1 mechanosensitive channel (tension-gated reporter) + Ca→contractility feedback stub** — `register-KB-then-build` · effort S · deps: H.8 plasma membrane (membrane_surface.py, commit feac927) MUST be present — it i
  - gate: Analytic ground-truth gates on the KB-3.10 logistic (FF units, γ_half=5000, γ_s=1500 pN/µm): (1) at σ=γ_half exactly P_open=0.5 (machine-precision). (2) PHYSIOL
  - KB: KB-3.10 READY and sufficient for the reporter: P_open=1/(1+exp(-(γ-γ_half)/γ_s)), γ_half~5 pN/nm=5000 pN/µm, γ_s~1-2 pN/nm=1000-2000 pN/µm; core logistic anchor
- **Filopodial fascin crosslinker (replace filamin stand-in) + tip FA clutch** — `register-KB-then-build` · effort M · deps: HARD BLOCKER: fascin SE→KC KB registration must land + PI-ratify FIRST. Then dep
  - gate: Euler buckling N²-scaling discriminator (analytic, KB-3.8: L_c=π√(κ_bundle/F); tight κ=N²·κ_single, loose κ=N·κ_single; κ_single=ℓ_p·kBT=17µm×4.1e-21J=7e-26 N·m
  - KB: MUST register SourceEvidence→KnowledgeClaim FIRST for fascin kinetics — gating dependency, HALT until PI ratifies. KB-3.8 is NOT fascin: it is "Filopodium struc
- **Glycocalyx pericellular brush (Alexander-de Gennes steric coat + emergent FA-maturation priming)** — `register-KB-then-build` · effort M · deps: Membrane compartment (aleph/laws/membrane_surface.py — PRESENT, commit feac927):
  - gate: TWO gates. (1) ANALYTIC (always-on, primary): the kernel P(D) must reproduce the closed-form Alexander-de Gennes law P(D)=(kT/s³)[(L_g/D)^{9/4}−(D/L_g)^{3/4}] t
  - KB: HALT-gated: no glycocalyx KnowledgeClaim exists (verified — KC search for glycocalyx/brush/pericellular/Paszek all empty; Paszek only in paper_chunks reference 
- **Septin filaments — cortical-rigidity reinforcement + micron-curvature-sensing membrane compartment (FF)** — `register-KB-then-build` · effort L · deps: Membrane compartment (ff/membrane_surface.py, commit feac927) — REQUIRED: the cu
  - gate: PRIMARY (analytic, no migration machinery needed) — CURVATURE-PREFERENCE PEAK: run septin_curvature_attach on fixed test spheres of R = {0.5, 1, 2, 3, 5, 10} µm
  - KB: EMPTY — confirmed by direct DuckDB query: no knowledge_claim row matches septin / diffusion barrier / cortical rigidity / constricted migration / curvature sens
- **Spectrin-ankyrin sub-membrane skeleton (2D triangulated shear network under the plasma membrane)** — `register-KB-then-build` · effort M · deps: membrane_surface.py plasma-membrane mesh must be the ACTIVE surface (commit feac
  - gate: Virtual simple-shear on a flat periodic triangulated patch built with the SAME kernels: apply small pure-shear γ_shear=0.02, measure induced 2D shear stress via
  - KB: ABSENT — must register SourceEvidence→KnowledgeClaim FIRST (HALT-gate before any μ0 constant enters code). No shear-modulus KC exists; KB-3.18 only lists "spect
- **Nucleolus / sub-nuclear bodies as viscous LLPS inclusions in the FF nucleus** — `register-KB-then-build` · effort M · deps: Nucleus compartment (resolve_nucleus + nucleus_shell_kernel + whole-cell driver)
  - gate: Viscous-capillary relaxation time (analytic ground truth, Caragine 2018 inverse-capillary-velocity): seed an ellipsoidal nucleolus (or two touching droplets rad
  - KB: ABSENT — must register SourceEvidence→KnowledgeClaim FIRST (this is the gating step; nothing to anchor a build to today). DuckDB scan of knowledge_claim + paper
