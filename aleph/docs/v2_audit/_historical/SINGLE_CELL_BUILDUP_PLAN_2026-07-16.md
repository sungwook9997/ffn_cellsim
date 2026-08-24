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

# Single-Cell Dynamic FEM+CFD — Incremental Build-Up Plan (2026-07-16)

**Professor's bottom-up method × PI's 3-thread**, grounded in the 183-gap inventory
(`wf_072fc794`) and the now-registered KB drafts (79 SE + 73 KU, `KB_UPLOAD_MANIFEST_2026-07-15.md`).
Supersedes nothing; it is the execution ladder that sits on top of `HIERARCHICAL_MECHANICS_VALIDATION_PLAN_2026-07-14.md`
(S1–S3 done) + `DYNAMIC_FEM_CFD_UPGRADE_ROADMAP_2026-07-15.md` (the two upgrade axes) +
`CELL_MECHANICS_FRAMEWORK_2026-07-15.md` (the 10-compartment spec).

---

## 0. Method (do not deviate)

1. **Smallest first, one physics per stage, validate-before-next.** Never add two compartments
   at once. The professor's rule; sharpened by CLAUDE.md.
2. **KB-first gate.** A stage may not start until its source model is a **registered KU** (the
   79 SE + 73 KU drafts are in; each stage cites its KU below). Magic-number params HALT to PI.
3. **`validate` ≠ experimental comparison** (PI clarification). Each stage's validation is two
   build-fidelity checks, NOT "does it match a measured MCF7 cell":
   - **(i) KB-fidelity** — implemented exactly as the TAG-sourced model says (right law, right
     param, right mechanism). Cross-checked against the KU + an analytic oracle.
   - **(ii) phenomenon-completeness** — captures every feature the reference physics exhibits;
     a checklist of "did we miss anything like the AFM time-force saturation?" (S1 Finding 7 was
     exactly this class of catch). Nothing in the reference behaviour is silently dropped.
4. **3-thread per stage** (PI 2026-07-15 — this is NOT mere mechanical part-assembly):
   - **(a) Mechanics gate** — analytic oracle + the two fidelity checks above.
   - **(b) Role-by-ablation** — run WITH and WITHOUT the new compartment at native scale; quantify
     what the cell *gains/loses*. This shows the compartment's **role for the cell**, not just its
     stiffness. (S3 already did this: remove nucleus → −12% force.)
   - **(c) Dynamic-capability arc** — show how the new piece moves the engine toward expressing
     **dynamic behaviour** (local stretch/shrink → … → protrusion/migration). The local-deform
     capability (S2) is the seed; each stage adds what a real protrusion needs.
5. **Native + full, always.** Every authoritative run at `--cortex-fil 38000 --from-resting`
   NF≈70686 with all already-built compartments ON at physiological values. Coarse = dev smoke
   only, labelled non-authoritative (CLAUDE.md HARD).
6. **Terminology restraint** (professor). Before the active machinery exists (< S9), an outward
   deformation is **"local outward deformation under load"**, never "protrusion". The word
   "protrusion" is earned only at S9.
7. **Double-count guards** (framework §3) checked at every stage that could trip one — especially
   **retire the 6πηR node drag the moment the Biot fluid lands (S7), same commit** (η double-count).

Every stage lands: a native run + the two-gate result + the ablation delta + a stress/strain
interactive HTML + MP4 (per the viz rules) + a `REPORT.md` stanza + a plan Finding.

---

## 1. Current position (validated base)

| Stage | Added | Status | Key result |
|---|---|---|---|
| **S1** | bare cortex shell | ✅ | drained cortex = clean LINEAR-SHELL (r²=0.992, E≈3.6 kPa); loading (not layer) sets modulus |
| **S2** | localized load (2A/2B) | ✅ | sign-reversal PASS; compression stiffens / tension ≈linear; localized boundary layer (ratio 2695). **= the dynamic-capability seed: the engine CAN localize a stretch/shrink** |
| **S3** | nucleus (incompressible) | ✅ | soft inclusion, +12% at ε=0.45; volume-conserved bulge; role-by-ablation already demonstrated |

---

## 2. The build-up ladder (single cell)

Each stage adds ONE compartment/physics on the validated base, KB-fidelity + phenomenon-complete
gated, with role-by-ablation and a dynamic-capability contribution. Ordered smallest/most-foundational
first and by dependency.

### S4 — Plasma membrane as a real 2-D shell (currently lumped 2γ/R)
- **Add:** wire `membrane_surface.py` (its own mesh) into the native driver — Helfrich bending κ_m,
  area-elasticity K_A (Rawicz), ERM membrane–cortex breakable tethers. (Machinery exists, un-wired —
  smallest real addition.) **KU:** KB-3.B1.2 (κ_m), KB-3.B1.3 (K_A), KB-3.B1.4 (ERM) *verified*.
- **(a) mechanics:** tether force/radius `f_t=2π√(2κ_m T_m)` + area-change oracle. *phenomenon:* area
  stays ~constant (very stiff in area) while bending is soft → shape change **unfolds reservoir**, does
  not stretch the bilayer; lysis at ~2–3% areal strain.
- **(b) role:** ablate membrane → cell loses the area reservoir + osmotic barrier; with it, a squeezed
  cell accommodates by unwrinkling, not by cortex over-stretch. **Role = area accommodation + barrier.**
- **(c) dynamic-capability:** the membrane is the surface a protrusion must push — wiring it is the
  prerequisite for any bleb/protrusion to have something to inflate.
- **Guard:** HAZARD-B1 (membrane tension + cortex tension double-add) — membrane sheet REPLACES the
  lumped 2γ/R, not adds to it.

### S5 — Pre-stressed osmotic baseline + LINC (the physiological baseline + nucleus load-path)
- **Add:** resting osmotic pressure ON at the physiological van't Hoff Π₀ (not zero); LINC nucleus↔cortex
  coupling so the nucleus is load-coupled, **not floating** (fixes the S3 "nucleus pokes through" artifact).
  **KU:** KB-DRAFT-3.B-04, -09 (LINC); volume/osmotic KUs; framework §7.
- **(a) mechanics:** the cell rests at a **pre-stressed equilibrium** — cortex pre-tensioned by turgor
  (γ=ΔP·R/2), not force-free. *phenomenon:* baseline is pre-stressed (not a relaxed bag); LINC transmits
  cortex indentation load into the nucleus (measurable nuclear strain under cortex load).
- **(b) role:** ablate turgor → floppy bag, cortex un-tensioned (the 2026-06-04 failure); ablate LINC →
  nucleus decouples / pokes through the cortex (the S3 artifact). **Role: osmotic = the pre-stress
  generator; LINC = the load path to the nucleus.**
- **(c) dynamic-capability:** internal pressure is what *drives* blebs and pressure-based protrusion —
  no pre-stress, no protrusion.
- **Guard:** HAZARD-B2 (turgor's Young–Laplace partner re-added as passive tension) + LINC/plate
  double-feed (land LINC + plate dim=Nc together).

### S6 — Microtubules + intermediate filaments (the internal load frame)
- **Add:** MT stiff polar beams with Euler buckling `F_crit=π²·EI/L²` (EI=k_BT·L_p) + MTOC; IF
  strain-stiffening cage (nucleus↔membrane). **KU:** KB-DRAFT-3-01/-02 (κ_MT, buckling); KB-DRAFT-3.B-06..13 (IF).
- **(a) mechanics:** MT buckling oracle (`F_crit`); IF nonlinear strain-stiffening. *phenomenon:* MT
  buckle at the predicted force under compression; IF stiffen at large strain (not before).
- **(b) role:** ablate MT → cell loses compression-bearing struts + nucleus centering; ablate IF →
  local overstrain propagates to global failure (IF distributes load). **Role: MT = compression frame +
  positioning; IF = integrity/load-distribution.**
- **(c) dynamic-capability:** MT position the nucleus + set polarity axis; IF keep the cell intact
  through the large deformations a protrusion/migration imposes.

### S7 — Cytoplasm CFD: the Biot pore-pressure field + FSI (**FEM → FEM+CFD milestone**)
- **Add:** spatially-resolved Biot pore-pressure field p(x,t) on a device grid (Darcy pore-flow, D≈40–60
  µm²/s), α-Biot coupling, immersed-boundary two-way FSI. **Retire the 6πηR node drag in the same commit.**
  **KU:** poroelastic (materialized) + KB-DRAFT-3.B-18/-19/-20/-21/-27 (Darcy field, RPY mobility,
  Brinkman, two-way coupling, spatial pressure).
- **(a) mechanics:** poroelastic indentation relaxation τ_p ~ L²/D emerges (~0.5–1 s cell-scale).
  *phenomenon (the PI's "AFM saturation" concern, now spatially resolved):* the indentation force
  relaxes on the **poroelastic τ_p**, not only the 0-D efflux τ_osm — and the drainage front is spatial,
  not uniform. Nothing of the measured poroelastic response is dropped.
- **(b) role:** ablate the fluid field → no cytoplasmic streaming, no bleb inflation, no spatial pressure
  equilibration (only a 0-D volume clock); with it → pressure redistributes through the cell and a bleb
  can inflate. **Role: the fluid is the medium that transmits pressure + inflates protrusions.**
- **(c) dynamic-capability:** the pore-pressure field is what physically **inflates a bleb / pressure-driven
  protrusion** — the first genuine fluid back-reaction.
- **Guard:** HAZARD-B4 — the 6πηR particle drag and the Biot fluid are the SAME dissipation; retire the
  drag the instant the fluid lands (this is the single biggest correctness trap).

### S8 — Dynamic fiber remodeling: KMC pool+mask (**static → living network**)
- **Add:** node/bond **pool + active mask**, a KMC event loop (polymerize / depolymerize / Arp2/3
  nucleate / cofilin sever / re-crosslink), G-actin monomer reservoir, explicit fiber-id, on the device
  spatial grid from S7. **KU:** KB-DRAFT-3-25/-26/-27/-28 (KMC scheduler, topology change, severing,
  capping/formin).
- **(a) mechanics:** turnover → emergent Maxwell stress relaxation (τ from the turnover rate); retrograde
  flow emerges from stress balance. *phenomenon:* the network **fluidizes past the turnover time (~30 s)** —
  elastic short-time, fluid long-time (the cortex's known crossover), not a permanently elastic solid.
- **(b) role:** ablate remodeling → the fixed-topology network cannot reconfigure (this is exactly the
  crawl-disp=0 wall we hit); with it → the network flows and rebuilds. **Role: remodeling = the ability
  to change shape irreversibly / move mass.**
- **(c) dynamic-capability:** **this unblocks migration** — a static network can deform elastically but
  cannot crawl; remodeling is what lets a protrusion persist and the rear retract.
- **Guard:** the pool+mask + device grid are shared infra with S7's IBM — build the grid once (S7), reuse.

### S9 — Active machinery wired + mechanosensing→contractility (**protrusion EMERGES**)
- **Add:** close the two dead feedback gains (mechanosensing→active stress, polarizer→traction); upgrade
  myosin toward the head-cycle where fidelity demands it. **KU:** KB-DRAFT-3-13..17 (myosin head-cycle,
  duty ratio, minifilament assembly), KB-DRAFT-3.B-14/-15/-17.
- **(a) mechanics:** contractile tension now **responds to load** (mechanosensing closes the loop).
  *phenomenon:* contractility rises with load / a polarized cell breaks symmetry.
- **(b) role:** ablate active → passive cell, no directed motion; with active + polarization + S8
  remodeling + S7 fluid + S4 membrane → **a directed protrusion forms and the cell can crawl.**
- **(c) dynamic-capability:** **the arc completes.** From S9 an outward, active-machinery-driven
  deformation IS a protrusion — the terminology-restraint gate lifts. Migration is now expressible.

---

## 3. The dynamic-capability arc (why the local stretch/shrink matters — PI 2026-07-15)

The professor's plan validates *mechanics*; the PI adds that we must **show the engine can express
dynamic movement (protrusion)**, and the S2 local stretch/shrink is the first step of that proof:

```
S2 local deform (CAN localize)  →  S4 membrane (a surface to push)  →  S5 pre-stress (a driver: pressure)
→  S7 fluid field (inflates)     →  S8 remodeling (persists / moves mass)  →  S9 active+polarity (DIRECTS it)
                                                                              = protrusion / migration EMERGES
```

Each stage's thread-(c) is one piece the eventual protrusion requires. We never *claim* protrusion early
(terminology restraint); we show the **capability accreting**, and demonstrate it as real motion only at S9.

---

## 4. Ordering rationale, dependencies, guards

- **Smallest/most-foundational first:** S4 (wire existing membrane) < S5 (turn on physiological baseline)
  < S6 (add MT/IF struts) < S7 (build the fluid field — the big one) < S8 (dynamic topology on S7's grid)
  < S9 (close feedback loops). Each is one physics.
- **Dependencies:** S5 LINC needs S3 nucleus (done). S7 fluid needs the device grid → also serves S8's
  partner search + IBM (build once). S8 needs S7's grid. S9 needs S8 (a static network can't crawl) + S4
  (a membrane to protrude) + S5 (a pressure driver).
- **KB-first:** every stage's KU is registered (drafts, PI to ratify). No stage starts on an unverified
  premise; magic numbers (η_s, k_LINC, MT DI rates, per-head F_stall, IF n_fil, RESERVOIR_STRAIN) HALT
  to PI (they are in the 40-item PI-sign-off list).
- **Double-count guards** at the trigger stages: B1 (S4 membrane replaces lumped tension), B2 (S5 turgor),
  **B4 (S7 retire 6πηR — same commit)**, B3 (S9 active-stress must not co-fire with explicit myosin).
- **Multicell (Part II)** — adhesion (clutch exists) → ECM (Mikado exists) → cadherin junction (FF, from
  the DCM catch-bond) → spheroid — runs *after* the single cell is complete, per the professor's S5–S8.

---

## 5. Per-stage closeout template (each Sx)

1. Native `--from-resting` run WITH and WITHOUT the new compartment (role-by-ablation).
2. Gate (a): analytic-oracle match (KB-fidelity) + phenomenon-completeness checklist signed.
3. Delta (b): quantified role (what the cell gains/loses).
4. Thread (c): the dynamic-capability contribution, demonstrated (not claimed as protrusion pre-S9).
5. Viz: stress/strain interactive HTML + MP4; `REPORT.md` stanza; plan Finding; figure via `mech_hier_vis.py`.
6. KU ratification check: the stage's KU is VERIFIED (PI) before the result is authoritative.
7. Only then unlock Sx+1.

---

### Change log
- 2026-07-16: created. Professor's method × PI 3-thread, on the S1–S3 base + the FEM+CFD/dynamic axes,
  gated on the 79 SE + 73 KU KB drafts. Awaits PI review of the ladder + priority (start S4, or jump the
  device-grid/S7 fluid, or the S8 dynamic axis first) + KU ratification before build begins.
