# AC Engine integration plan — component/connector architecture as canonical (2026-07-22)

**PI decision (2026-07-22):** `aleph/engine/` (the component/connector "component cell engine") is the
**canonical target architecture**. All new work — and the governance docs (CLAUDE.md, STRUCTURE.md, memory)
— align to it. This document is the integration plan; the final section is its adversarial audit.

Lead: Claude (this session). Branch `codex/ff-ac-codex`. Grounded in `aleph/engine/contracts.py` (read in full),
the existing `docs/v2_audit/cell_engine/` design set, the R0/R1 checkpoint Notion page, the ECM audit, and
this session's resting-baseline diagnosis.

> **Scope note.** The engine ARCHITECTURE is already documented by Codex in `docs/v2_audit/cell_engine/`
> (12 files: `CELL_ENGINE_ARCHITECTURE.md` = the architecture contract, `ROLLING_ROADMAP.md` = R0–R15+/R16
> with the 8-state evidence ladder, `COMPONENT_BUILD_MATRIX.md`, + 9 per-component plans). This plan does NOT
> re-specify the architecture — it **defers to those as the architecture SoT** and covers only what is not yet
> settled: (i) aligning the project constitution (CLAUDE.md/STRUCTURE.md/memory) to the ratified direction,
> (ii) how the OPEN resting baseline reconciles with the engine, (iii) sequencing/coordination + the PI-gated
> items. Where this plan and the `cell_engine/` docs disagree, the `cell_engine/` docs win.

**Evidence ladder (from `ROLLING_ROADMAP.md`, adopt this vocabulary):**
`CONTRACTED → SEAMED → KERNEL_BOUND → CUDA_UNIT → CONNECTED → NATIVE → OPTIMISED → PRODUCTION` (skipping a
state is forbidden). Today the engine is **predominantly SEAMED**; only `load_path.py` (the SF–FA–ECM series
joint) is KERNEL_BOUND. "Structural tests pass" ≠ production (`COMPONENT_BUILD_MATRIX.md:16`).

---

## 1. What ac/engine IS (the design, from `contracts.py`)

A declarative **ownership + connectivity** graph — *"never treats co-location in an array as a mechanical
connection"* (contracts.py:4). Two first-class abstractions:

**Components (13)** — state-owning bodies, each a `ComponentRole`:
- SURFACE_BODY: `membrane` (fluid Helfrich surface FEM), `cortex` (explicit crosslinked F-actin; gated surface condensation)
- FLUID_VOLUME: `cytosol` (Biot/Darcy finite-volume)
- CORE_BODY: `nucleus` (native lamina/chromatin FEM; gated modal reduction)
- ACTIVE_LOAD_PATH: `sf_arc`, `focal_adhesion` (clutch KMC, no geometry), `lamellipodium`, `filopodium`, `nmii` (Stam-Hocky backbone + individual heads)
- STRUCTURAL_RIG: `microtubule` (dynamic rod graph), `intermediate_filament` (nonlinear cable graph)
- ENVIRONMENT: `ecm` (fiber/crosslink world), `world_boundary` (far-field frame; static)

**Connectors (~32)** — bidirectional inter-component force paths, each a `ConnectorFamily` (ERM, FA_CLUTCH,
ACTIN_ANCHOR, TRANSIENT_ACTIN, LINC, PLECTIN, SPECTRAPLAKIN, MOTOR, IMMERSED_TRANSFER, FLUID_BOUNDARY,
FIBER_CROSSLINK, ENVIRONMENT_BOUNDARY, CONTACT). Each carries: `kinetics`+`commit_on_accept` (kinetic joints
commit ONLY on an accepted physical step), enforced `bidirectional`+`adjoint_transfer_required` (Newton's
3rd law), a `chemistry_card` (the molecular law, e.g. `nmii_head_actin_hill_bell`, `alpha2beta1_collagen`),
and optional `mechanical_group` (semantic edges resolved by one composite joint — e.g. `fa_actin_anchor` +
`integrin_collagen_clutch` → the single `alpha2beta1_collagen_series` runtime), plus `generation_required`/
`remap_on_accept`/`blocks_sleep_refine` for topology-changing joints.

**Why this is the right architecture** (it operationalizes existing HARD rules):
- It is the structural encoding of the PI cell-mechanics framework (cell = load-transfer; 10 compartments;
  pre-stressed; the six double-count traps). Components = load-bearing compartments; connectors = the ONLY
  mechanical load paths. `force_path(a,b)` / `mechanical_group()` make load transfer explicit and auditable.
- `co-location ≠ connection` is the structural form of the **anti-lumping / no-double-count** rule: sharing a
  device array is not a mechanical connection — only an explicit `ConnectorContract` transfers load.
- `bidirectional` + `adjoint_transfer_required` structurally enforce Newton's 3rd law (no one-way loads).
- `kinetics ⇒ commit_on_accept` structurally enforces the one-outer-physical-clock / accepted-step transaction.

## 2. Current state (grounded, not aspirational)

| Layer | What it is | Status |
|---|---|---|
| `ac/engine/` | contracts (13 comp / 32 conn) + runtime/world/dispatch (7 candidate phases, device-resident accepted-step transaction, one `rng_seed`/step, exact-once canonical-facade dispatch across 8 facades) + per-component seams | **~173 engine tests (172 CPU-green + 1 CUDA-gated on load_path)**. **SEAMED**, not production: only `load_path.py` is KERNEL_BOUND (real segment-joint force/commit/rollback kernels); the other 9 rigs are protocol/spy seams that prove ownership/transaction coverage and are DESIGNED to launch existing `ac/*`+`ff/*` kernels via injected delegates but are not yet bound. The rigs even reject their own landed backends as non-production (MT = STATIC bead adapter, IF = monolithic-Hookean reference, NMII landed MyosinForce = reference/transition only). |
| `ac/cell/` | build_cell / driver / compartments — the **running physics** | cortex·myosin·nucleus·membrane·fluid·steric·pressure compose + build at native 70,686. **Resting-baseline convergence OPEN** (see 2026-07-22c diagnosis). |
| `ff/` | audited Warp kernels (Cytosim mechanics, ECM library, network) | reusable; the kernels the engine seams must bind to. |
| `ac/solid`, `ac/weave`, `ac/motor`, `ac/fluid` | isolated compartment physics (MT/IF/clutch/steric, SF/lamellipodium, NMII, Biot) | built + tested + `accumulate(pos,f)`-ready; not all wired into build_cell. |

Branches (unmerged): `codex/ecm-gpu-topology` (+3: CUDA Mikado topology + transactional remodeling),
`ac/sf-wiring` (SF), `ac/mt-if-linc-vertical`, `ac/sf-fa-nmii-vertical`.

**ECM audit (Codex):** FF ECM material cards + Warp constitutive kernels reusable; quarantine all host
NumPy/cKDTree builders, host relaxation/readback, CPU integrators, host FA attachment; build a full GPU
device schema + topology/contact pipeline. 🔴 **Collagen modulus gate is self-inconsistent** — `aleph/laws/ecm_library.py`
30–100 Pa vs `ff_ecm_validate.py` band 5–100 Pa (its note even cites the 30–100 KB-1.30 anchor). No production
collagen constitutive law/gate is closed. **PI/Notion SoT normalization required** (gate-contract change).

## 3. The reconciliation model — engine OWNS, kernels EXECUTE

ac/engine does **not** reimplement physics. Each component/connector seam **binds injected delegates that
launch the existing `ac/cell` + `ff/` Warp kernels**, while the engine adds: single-owner authoritative CUDA
state, the accepted-step transaction + deterministic RNG epoch, and the adjoint bidirectional force scatter.
`ac/cell/build_cell` stays the authoritative runtime until a given vertical slice is bound and gated.

**Component → existing kernel binding (target):**

| Engine component | Existing physics to bind | Status |
|---|---|---|
| membrane | `aleph/components/incumbent/compartments.py` Helfrich + `membrane_pressure.py` | ac/cell live; engine seam unbound |
| cortex | `aleph/laws/cortex_assembly.py` + `aleph/laws/network_warp.py` link/bending | ac/cell live; seam unbound |
| cytosol | `ac/fluid/` Biot/Darcy | ac/cell live; seam unbound |
| nucleus | `aleph/components/incumbent/compartments.py` lamina/LINC | ac/cell live; seam unbound |
| nmii | `ac/motor/` Stam-Hocky heads (Hill/Bell) | isolated; seam unbound |
| microtubule / IF | `ac/solid/{microtubule,intermediate_filament}.py` | isolated + `codex/*-vertical`; seam unbound |
| focal_adhesion / ecm | `aleph/components/solid/adhesion_clutch.py` + `ff/ecm_*` (+ topology branch) | isolated; **gated on collagen SoT** |
| sf_arc / lamellipodium / filopodium | `ac/weave/{stress_fiber,lamellipodium,...}` | isolated / `ac/sf-wiring`; seam unbound |

**Connector → chemistry:** each `chemistry_card` (ezrin, α2β1-collagen, nmii-Hill-Bell, nesprin-LINC,
plectin, spectraplakin, brownian-ratchet) is a KB-anchored molecular law; the connector seam binds the
existing per-bond Warp kinetics (Bell/Hill/catch-slip) and enforces adjoint scatter + commit-on-accept.

## 4. Integration sequencing — DEFER to `ROLLING_ROADMAP.md`; slices map to its R-waves

The authoritative sequence is `ROLLING_ROADMAP.md` (R0–R15+). The slices below are the SAME sequence named for
this plan's discussion — they must not diverge from the roadmap. Each slice = compose its components, bind its
connectors to native kernels (SEAMED→KERNEL_BOUND→CUDA_UNIT→CONNECTED→NATIVE), and pass a native accepted-step
gate on the ONE A5000 (the gate chain is strictly serial — `AC_PARALLEL_SESSIONS_2026-07-16.md`) before the
next. Mapping: Slice 1 = R2 (surface/fluid binding); Slice 2 = R3+R4 (ECM GPU port + SF-FA-ECM/NMII); Slice 3
= R5 (MT/IF/LINC); Slice 4 folds into R4; Slice 5 = R6; terminal = R8 whole-cell native baseline.

1. **Slice 1 — membrane ⊕ ERM ⊕ cortex ⊕ cytosol (the resting baseline).** This is both the engine's first
   real physics binding AND the resting-baseline closure. **Under the engine, the resting pre-stress is
   established through the explicit ERM/cortex connectors AT BUILD** (physiological baseline: cortex hoop
   pre-tension, ERM pre-holding turgor), NOT by relaxing a force-free mesh with a static solver — which
   resolves the architectural drift this session surfaced (§6).
2. **Slice 2 — sf_arc ⊕ focal_adhesion ⊕ ecm** (the `alpha2beta1_collagen_series` mechanical group). Gated on
   the collagen-gate SoT normalization + the ECM GPU topology (codex/ecm-gpu-topology).
3. **Slice 3 — microtubule / intermediate_filament ⊕ LINC ⊕ nucleus** (tensegrity: MT compression, IF cage).
4. **Slice 4 — nmii MOTOR connectors** onto cortex/sf/lamellipodium/filopodium (active tension emerges).
5. **Slice 5 — lamellipodium / filopodium CONTACT + nascent FA** (protrusion, brownian ratchet).
6. **Terminal — full native population run + emergence detectors + bleb validation.**

## 5. Governance doc updates (the "모두 다 수정" deliverable)

- **CLAUDE.md**: (a) Stack §Engine — name `ac/engine` the canonical composition layer over `ac/cell`+`ff/`
  kernels, and point to `docs/v2_audit/cell_engine/` (esp. `CELL_ENGINE_ARCHITECTURE.md` + `ROLLING_ROADMAP.md`)
  as the architecture SoT. (b) Add an **Architectural-principle** entry: component/connector — components own
  state, connectors are the ONLY mechanical connections (`co-location ≠ connection`), bidirectional+adjoint,
  kinetic⇒commit-on-accept, one device accepted-step transaction. (c) Extend the I0-A contract to reference
  `reference_cell_architecture()` (13 components / 32 connectors) + the 8-state evidence ladder as the
  composition contract and the "structural tests ≠ production" bar. (d) Repo-layout: add `ac/engine/` as the
  orchestration layer and demote `ac/cell` build path to "incumbent runtime / compatibility adapter" per
  `CELL_ENGINE_ARCHITECTURE.md:161-175` (no existing runtime is deleted mid-migration).
- **STRUCTURE.md**: add the `ac/engine/` layer and a component→kernel map; mark `ac/cell/build_cell` as the
  incumbent runtime being progressively bound under the engine.
- **Memory**: a `project` memory recording ac/engine as the ratified canonical architecture + the seam→kernel
  binding model, linking [[reference-cell-mechanics-framework]] and [[project-active-cell-c-redesign]].

## 6. How the resting baseline relates to ac/engine — it does NOT close "by construction" (corrected)

⚠️ **Corrected after adversarial audit.** An earlier draft of this section claimed the resting baseline is
"balanced by construction" under the engine and deprioritized the solver work. **That is refuted by this
session's own `RESTING_BASELINE_DIAGNOSIS_2026-07-22c.md`:** pre-stressing at build (ERM preload + cortex
pre-tension = exactly "construction") does NOT balance t0 — it relocates the residual into the cortex reaction
(~2.22 pN) and the coupled rest-length preload "stuck at 2.22 / diverges; option (1) reduces to option (2)".
`evaluate_preload_capacity`/NG-3 are a pre-tension **capacity/existence** check, NOT the strict projected-force
**convergence** gate (~0.21 pN) that is actually failing — citing them as closure was a gate substitution.

**Correct position:** the resting baseline is an OPEN physics problem in `ac/cell`, independent of the engine
re-architecture. Moving the arrays under a component owner does not change the ill-conditioned `K δx = −f`.
It must be closed FIRST in ac/cell (next probe = the FD-vs-operator Jacobian check from the diagnosis), then
wrapped. See §9.

## 7. PI-gated items — DECIDED 2026-07-22 (see `AC_DECISION_CARDS_2026-07-22.md`)
1. **Collagen modulus** → (c) concentration-resolved conditional gate (Card 1); P2 executes.
2. **MCF7 ERM density** → production HOLD; capacity-gate-only; f_rupt≠single-ERM correction (Card 2); P2.
3. **MCF7 γ_cortex** → re-register Hosseini 2020 direct MCF-7 value + resolve the γ_mem double-count (Card 3); P2.
4. **Sequencing** → run only P0/P1/P2 now; P3–P7 deferred (Card 4).
5. **ac/cell retirement** → gate-based 3-stage strangler; feature-freeze now (Card 5).

## 8. Adversarial audit outcome (3 independent auditors + Lead verification)

Three adversarial auditors red-teamed this plan (architecture/physics, sequencing/feasibility, governance/
fact-check). Verdicts: **FUNDAMENTALLY-FLAWED / FUNDAMENTALLY-FLAWED / SOUND-WITH-FIXES**. The fact-check
auditor confirmed every hard number (13 comp / 32 conn, ~173 tests, the collagen conflict, the 2.22 pN, all
branch/file claims) — the defects are governance/sequencing, not facts. Convergent findings, after Lead
verification (I did NOT accept adversarial claims at face value):

| # | Finding | Verdict after verification |
|---|---|---|
| A1 | The 13-component set makes cortex/sf_arc/lamellipodium/filopodium **separate owned bodies**, which the running `aleph/components/weave/woven_cell.py` (ONE unified network, `region_id` diagnostic-only) and the PI emergence rule treat as label-blind roles. Two auditors called this "inverts the no-double-count HARD rule." | **REFRAMED, then RESOLVED.** The "inverts a rule" framing is overstated: `LOAD_PATH_PLAN.md:12,105,233` *deliberately* forbids one giant WovenCell and shared-node welds, and the engine ships anti-double-count guards (`BorrowedSegmentActorView` non-owning views, dispatch dedup, `mechanical_group`, "double-force guard mandatory"). It is a genuine **cortex/SF ownership disagreement** (ac/weave emergent-partition vs ac/engine separate-component), not a bug. **PI resolved 2026-07-22: "SF는 별도로 한다" → separate component is ratified.** |
| A2 | §6 "resting baseline closes by construction" is refuted by this session's own 22c diagnosis; it relocates the residual to the cortex (2.22 pN); NG-3 is a capacity/existence check, not the failing projected-force gate — a gate substitution. | **CONFIRMED.** §6 corrected above. The resting baseline is OPEN in ac/cell; the engine does not fix it. |
| A3 | Architecture-ahead-of-physics: canonicalizing + rewriting the constitution now, while the ONE gate for all dynamics is OPEN and the engine has **never run native** (173 tests all mock the device; the single real force-kernel test is CUDA-skipped and has never executed). | **CONFIRMED.** Governance scope reduced (§10): record the ratified DIRECTION + principle + evidence ladder; do NOT claim PRODUCTION or deprecate ac/cell. |
| A4 | "Engine binds existing kernels (not reimplement)" is loose: no engine module imports `ff`/`ac.cell`/`ac.solid`/`ac.motor`/`ac.weave`; only `load_path.py` has real kernels. | **REFRAMED.** Per the per-component plans the binding is a **spectrum**: thin-adapter reuse (surface/membrane), kernel-reuse-with-relayout (NMII), and full SoA **port/re-host** (ECM, `ROLLING_ROADMAP` R3). "Non-invasive wrap" was too clean; each component's binding cost is real and stated in its plan. |
| A5 | Shared-branch merge risk already materialized (`implicit_mechanics.py` reverted mid-session by Codex/Syncthing); Codex owns unmerged `codex/ecm-gpu-topology`, `ac/*-vertical`. | **CONFIRMED.** §11 adds a fencing/coordination protocol. |
| A6 | Authority overreach: "canonical" is not on-disk; several open PI calls (deprecation, deprioritize-solver, MCF7 gaps) were written as decided. | **CONFIRMED.** Corrected — this plan records only what the PI ratified (direction + SF-separate) and lists the rest as open. |

## 9. Revised recommendation (post-audit + PI SF decision) — physics-first, scoped governance

1. **Close the resting baseline in `ac/cell` FIRST — it is the critical path and it is a physics problem, not
   an architecture one.** Next concrete step is the diagnosis's prescribed **FD-vs-operator Jacobian probe**
   (does `PKP` deliver the measured stiffness into the Newton step for the turgor/ERM mode?), then land a
   conditioned-implicit or targeted per-node coupled statics step, and pass the strict projected-force native
   gate (~0.21 pN) at subdiv ≥6 on the A5000. Only then is Slice 1 real. (Lead-owned, serial, single-GPU.)
2. **SF (and the protrusions) are separate components — RATIFIED (PI 2026-07-22).** Build SF/arc as its own
   filament population owned by `sf_arc`, coupled to cortex via the transient connector (never shared nodes),
   using `BorrowedSegmentActorView` for the joint. This SUPERSEDES the ac/weave "SF = emergent label-blind
   cortex partition" implementation for the engine. The no-double-count invariant is preserved by construction:
   **each physical filament belongs to exactly ONE component; components couple only through explicit
   connectors.** Update the CLAUDE.md rule accordingly (§10).
3. **Scoped governance now** (§10): record ac/engine as the ratified TARGET architecture + the component/
   connector principle + the 8-state evidence ladder + the reframed no-double-count rule. Do NOT mark ac/engine
   production or ac/cell deprecated (evidence ladder: it is SEAMED, never run native).
4. **Unblock ECM** (parallel): PI/KB collagen-gate SoT normalization (30–100 vs 5–100 Pa) + land the
   `codex/ecm-gpu-topology` SoA port.
5. **Parallelize the seam→kernel authoring** (§11): each component's binding is a disjoint module authored
   CPU-green to CUDA_UNIT in its own worktree; NATIVE/CONNECTED gates serialize behind Slice 1 on the one GPU.

## 10. Governance edits (scoped — apply now)

- **CLAUDE.md**: (a) §Stack Engine — name `ac/engine` the ratified canonical composition layer, SoT =
  `docs/v2_audit/cell_engine/`; keep `ac/cell` as the running incumbent (NOT "deprecated"). (b) New
  architectural-principle row: component/connector — components own state, connectors are the ONLY mechanical
  connections (`co-location ≠ connection`), bidirectional+adjoint, kinetic⇒commit-on-accept, one device
  accepted-step transaction. (c) **Reframe the unified-actin / no-double-count rule** (PI 2026-07-22): the
  cortex/SF/arc/lamellipodium/filopodium are now **separate components, each owning a DISJOINT filament
  population**; the invariant becomes "each physical filament belongs to exactly one component; components
  couple only via explicit connectors, never shared nodes." (d) Add the 8-state evidence ladder + "structural
  tests ≠ production." (e) Repo-layout: add `ac/engine/`. Do NOT edit the 70,686 native-population or
  Warp-only I0-A clauses (still binding).
- **STRUCTURE.md**: add the `ac/engine/` layer + a one-line component/connector map + the `cell_engine/` doc set.
- **Memory**: supersede [[project-ac-sf-wiring]] (SF-emergent) with SF-separate-component; add ac/engine-canonical.

## 11. Parallel launch pack — copy-paste boot prompts

**Coordination (BINDING — every session reads first):** (1) Runtime = Warp-CUDA only (I0-A); dev Mac authors
kernel SOURCE + CPU/NumPy oracles only; native gates run on the ONE gbook A5000 and **serialize** (each slice's
native cell contains the prior's live state — adding sessions does NOT compress the native chain; it only keeps
every module CPU-green + structurally tested before its GPU slot opens). (2) **Your worktree, your disjoint
`ac/engine/<module>.py` + its `tests/ac/engine/test_<module>.py` ONLY.** Do NOT edit `ff/`, `ac/cell/`, or
another session's files; write any shared-file change as an `INTEGRATION.md` patch-note for the Lead to apply
in spine order. (3) Target evidence state = **CUDA_UNIT** (KERNEL_BOUND + a CUDA unit gate authored; the
NATIVE/CONNECTED gate is Lead-owned + serial). (4) Bind to the EXISTING physics named in your per-component
plan under `docs/v2_audit/cell_engine/`; reuse laws/kernels, do not invent constitutive law; no magic numbers.
(5) gbook: `ssh gbook`, `PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim`, python `~/miniconda3/envs/ffn_sim/bin/python`.

> **PI 2026-07-22 — run ONLY P0 / P1 / P2 now.** Lead runs **P0 (critical path)** solo on the A5000; **P1**
> (SF component) is an isolated worktree; **P2** is expanded to the single **three-card evidence/SoT lane**
> (collagen + ERM density + Hosseini γ — see `AC_DECISION_CARDS_2026-07-22.md`). **P3–P7 below are DEFERRED**
> to their R-waves after P0's native GO — do NOT open them as sessions now (kept for reference/sequencing).

### P0 — Lead, critical path (resting baseline; NOT parallel — Lead owns the GPU)
```
You are the Lead. Branch codex/ff-ac-codex. Close the OPEN resting-baseline gate — it blocks ALL dynamics and
the whole ac/engine slice chain. Read docs/v2_audit/RESTING_BASELINE_DIAGNOSIS_2026-07-22c.md. Do the prescribed
FD-vs-operator probe FIRST: at subdiv 6 + overlap_free + radial ERM pairing + ERM preload, instrument the
analytic_implicit CG at the Newton-step level — dump ‖P F‖ vs ‖F‖ (is the residual in the projection
nullspace?), dump proposed dx on membrane/loaded-cortex nodes vs the 0.13 nm expected, and finite-difference
the true force to validate the _central_action crosslink/ERM tangent assembly in PKP (implicit_mechanics.py).
That decides conditioned-implicit vs targeted-per-node preload. Land the fix, drive the strict projected-force
native gate (~0.21 pN) under band at subdiv≥6 on the gbook A5000, then stop and report. No gate loosening.
```

### P1 — SF as a SEPARATE component (RATIFIED 2026-07-22)  [worktree: ac/engine-sf]
```
Worktree session. Own aleph/engine/stress_fiber.py + tests/ac/engine/test_stress_fiber.py ONLY. Read
docs/v2_audit/cell_engine/LOAD_PATH_PLAN.md + COMPONENT_BUILD_MATRIX.md. PI ratified SF as a SEPARATE
state-owning component (NOT an emergent label-blind partition of the cortex WovenCell). Build sf_arc as its
own DISJOINT F-actin population (ventral/dorsal bundles anchored at FAs), coupled to cortex via the transient
crosslink connector using BorrowedSegmentActorView (non-owning) — NEVER shared nodes / welds (LOAD_PATH_PLAN
lines 105,233). Invariant to prove in a test: each filament belongs to exactly one component (no cortex/SF node
aliasing). Bind the SF segment mechanics to the existing ac/weave/stress_fiber physics as a delegate; CPU-green
oracle + structural tests to CUDA_UNIT. Do NOT edit ac/cell or ff. INTEGRATION.md for any shared change.
```

### P2 — THREE-CARD evidence/SoT lane (collagen + ERM density + Hosseini γ)  [no worktree; KB/doc]
```
Own the three parameter/gate decision cards in docs/v2_audit/AC_DECISION_CARDS_2026-07-22.md. Do NOT pick
values to pass a gate; surface each to PI as a gate-contract change. Deliver three decision-cards/memos:

(1) COLLAGEN = concentration-resolved gate (Card 1). Split into VG-ECM-{G0(c), c-scaling n≈2.0-2.1,
stiffening, normal-stress}, each with T/frequency/strain-amplitude/3D-bulk. Demote 5-100 Pa to an exploration
envelope; move 30-100 Pa to 3 mg/mL; extract the 1.5 mg/mL band directly from the source figure. Fix the code
conflicts: ecm_library.py has band (30,100) AT ref_conc 1.5 mg/mL AND conc_exponent=0.5 while the scaling
target is n≈2.0-2.1; the KB-1.30 note is stale (use KB-1.32 / KB-1.V.2.1); re-check the Yang-Kaufman DOI.
Historical 11-15 Pa results = re-evaluation-pending until the new contract is PI-approved.

(2) ERM DENSITY = production HOLD (Card 2). rho_ERM*F_single >= residual pressure is a NECESSARY capacity gate
ONLY; do not adopt the back-computed min density as physiological. Keep 1/node + zebrafish 600/um^2 diagnostic-
only. Correct the record: f_rupt (11.4 pN) is a CONTINUUM tube-extraction scale, NOT single-ERM; single
ezrin-F-actin (Braunger 2014) = ~4.6 pN/nm stiffness (=k_erm), ~50 pN unbinding, k_off~1.3/s. Approve density +
bound/active fraction + Bell(k_on,k_off,F0,capture) as one bundle. Proxy priority: MCF7 proteomics x localization
> other human epithelial areal density > other-species embryo (diagnostic-only).

(3) HOSSEINI gamma (Card 3). Digitize Hosseini 2020 (MCF-7 suspended interphase AFM, n=27; ~0.27 mN/m, IQR
0.18-0.40) into a dedicated KnowledgeClaim; Hosseini 2021 (~0.41) = upper cross-check. It is currently only
figure-locked in KB-6.1.3 with no numeric gamma. Reconcile against the existing 0.35-0.65 KU-3.5 target.
RESOLVE the double-count: is Hosseini gamma total-effective (dP=2*gamma_eff/R) or cortex-only
(dP=2*(gamma_cortex+gamma_mem)/R)? First reading at R=7.5um -> dP~48-107 Pa (central ~72). Keep 40 Pa as the
HeLa diagnostic proxy until this lands. Use outputs/tag_kb/tag_query.py + the reference PDFs.
```

### (DEFERRED — do NOT open now; after P0 native GO, per R-waves) P3–P7

### P3 — Surface Body binding (membrane ⊕ cortex)  [worktree: ac/engine-surface]
```
Worktree session. Own aleph/engine/surface_body.py + tests/ac/engine/test_surface_body.py ONLY. Read
docs/v2_audit/cell_engine/SURFACE_BODY_PLAN.md. Bind the Surface Body seam's injected delegates to the EXISTING
Warp kernels: aleph/components/incumbent/compartments.py Helfrich/area + membrane_pressure.py traction + aleph/laws/cortex_assembly.py /
aleph/laws/network_warp.py link/bending (reuse "directly or through a thin adapter" per the plan §6). Keep membrane +
cortex as distinct SURFACE_BODY owners joined by the ERM connector; the facade never becomes a 3rd state owner.
CPU-green oracle + structural tests to CUDA_UNIT. Do NOT edit ac/cell or ff (patch-notes only).
```

### P4 — Fluid Core binding (cytosol ⊕ nucleus)  [worktree: ac/engine-fluid]
```
Worktree session. Own aleph/engine/fluid_core.py + tests/ac/engine/test_fluid_core.py ONLY. Read
docs/v2_audit/cell_engine/FLUID_CORE_PLAN.md. Inject the EXISTING ac/fluid Biot/Darcy field solver + domain +
IBM transfers and the ac/cell nucleus lamina/LINC kernels behind the Fluid Volume + Core Body delegates; close
mass / no-flux / adjoint-work in ONE candidate_iteration. Moving membrane/nucleus boundaries via the
FLUID_BOUNDARY connectors. CPU-green oracle + structural tests to CUDA_UNIT. Do NOT edit ac/fluid or ac/cell.
```

### P5 — ECM World GPU port (coordinate with Codex codex/ecm-gpu-topology)  [worktree: ac/engine-ecm]
```
Worktree session, GATED on P2 (collagen SoT). Own aleph/engine/ecm_world.py + tests/ac/engine/test_ecm_world.py.
Read docs/v2_audit/cell_engine/ECM_WORLD_PLAN.md + ROLLING_ROADMAP R3. PORT aleph/laws/ecm_library.py collagen-I
material/topology + reusable axial/bending/link primitives into ECM-owned CUDA SoA state (device schema,
biological IDs + generation counters, active/free-list accounting, GPU topology construction, segment
broadphase, reaction-preserving far-field boundary ledger). Quarantine ALL host NumPy/cKDTree builders + host
relaxation/readback. The α2β1-collagen clutch + collagen crosslink are graph connectors. COORDINATE with Codex
(codex/ecm-gpu-topology, +3 commits: CUDA Mikado topology + transactional remodeling) — do not duplicate; the
Lead merges. Do NOT close the collagen constitutive gate until P2 lands.
```

### P6 — NMII actuator binding (head-resolved)  [worktree: ac/engine-nmii]
```
Worktree session. Own aleph/engine/nmii_actuator.py + tests/ac/engine/test_nmii_actuator.py ONLY. Read
docs/v2_audit/cell_engine/NMII_ACTUATOR_PLAN.md. Bind the explicit Stam-Hocky backbone + individual heads
(Hill FV, per-head Bell kinetics) from ac/motor behind the nmii component + MOTOR connectors to sf_arc / cortex
/ lamellipodium / filopodium material-point ports; reject aggregate/linear runtime substitutes (they are
diagnostics only, I0-A). Reuse the ac/motor kernels, re-own the layout per the plan. CPU-green + structural
tests to CUDA_UNIT. Do NOT edit ac/motor.
```

### P7 — MT + IF rigs: promote from static/monolithic to dynamic  [worktree: ac/engine-mtif]
```
Worktree session. Own ac/engine/{microtubule_rig,intermediate_filament_rig}.py + their two test files ONLY.
Read MICROTUBULE_PLAN.md + INTERMEDIATE_FILAMENT_PLAN.md. The landed ac/solid/microtubule is a STATIC bead-chain
adapter and ac/solid/intermediate_filament is a monolithic-Hookean reference — the plans reject both as
production backends. Author the dynamic MT rod graph (dynamic instability, MTOC, plus-end) and the nonlinear IF
cable graph with turnover, behind the STRUCTURAL_RIG owners + their LINC/plectin/spectraplakin connectors.
Also fold the unmerged ac/mt-if-linc-vertical (+3) work. CPU-green + structural tests to CUDA_UNIT.
```

## 12. Per-stage visualization gate + P-VIS lane (PI 2026-07-22)

Visualization is a **ladder checkpoint, not only a PRODUCTION artifact** (CLAUDE.md). Every component reaching
`CUDA_UNIT`/`CONNECTED` must (a) render **itself + its connectors in ISOLATION** and (b) refresh the
**CUMULATIVE composition** render. Tool: `scripts/ac_cell_assembled_viz.py` — now carries per-compartment
**isolation scenes** (cortex-only / NMII-only / membrane-only / nucleus-only / load-paths-only) alongside the
composed scenes, data-driven (a compartment absent from the dump is skipped; a new component = one added scene
once its nodes are dumped), full-native resolution (no downsample), one self-contained WebGL HTML, browser-
verified via `scripts/browser_check.py`. Dump = gbook `dump_state.py`; render + browser-check = dev Mac (no CUDA).

**P-VIS runs as a 4th lane** (P0/P1/P2 are already running). It owns the viz tooling + the per-stage renders;
it does NOT touch physics.

### P-VIS — per-stage / per-compartment visualization lane  [worktree: ac-viz; no physics]
```
Worktree/viz session. Own scripts/ac_cell_assembled_viz.py, scripts/ac_cell_dynamics_viz.py,
scripts/dump_state.py, scripts/browser_check.py ONLY (+ new viz helpers). Read CLAUDE.md "Per-stage
visualization gate" rule. Do NOT edit any physics (ac/cell, ac/engine components, ff, solver, scheduler).

Mission: make the ac/engine staged build VISUALLY legible per component.
1. ac_cell_assembled_viz.py already has per-compartment ISOLATION scenes (cortex/NMII/membrane/nucleus/load-
   paths) + composed scenes; keep it data-driven so a NEW component (SF, MT, IF, ECM, lamellipodium, filopodium)
   becomes one added isolation scene the moment its nodes appear in the dump. Extend dump_state.py so each
   component's nodes/faces/connector-endpoints are dumped with a stable per-component key + a global unique
   filament/actor ID (so SF-separate etc. never double-draws a cortex node).
2. Add a CONNECTOR/load-path scene family: draw each connector family (ERM, FA-clutch, LINC, MOTOR, plectin,
   spectraplakin, transient-actin, immersed-transfer) as its own coloured line set between the two components,
   so the force PATH is visible (co-location != connection — the viz must show the explicit joint, not proximity).
3. Add a CUMULATIVE composition view that grows as slices bind, and a per-component |F|/load colour scale.
4. Every render: full native resolution (NO downsample), ONE self-contained interactive WebGL HTML, then
   browser_check.py to prove it renders (a grep is NOT verification; a visual screenshot catches wrong-sign /
   outlier / interpenetration bugs — e.g. the steric-hotspot scene).
5. Wire it so P0 (resting cell + residual heatmap) and P1 (SF isolation + its FA/cortex connectors) each emit
   their stage render at CUDA_UNIT/GO. Publish to the gh-pages gallery per the existing flagship workflow.

Dump on gbook (PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim); render + browser-check on the dev Mac
(no CUDA). Big HTML is fine (full-res). Commit end: Co-Authored-By: Claude Opus 4.8 (1M context)
<noreply@anthropic.com>. 한글 답변.
```


## 13. Overnight autonomous progress (2026-07-23) — vs the ROLLING_ROADMAP

Lead orchestration (P0 GPU + parallel CPU authoring + adversarial verify + integrate, looped). State now:
- **R0/R1 CLOSED**: the 13-component / 32-connector composed world is REGISTERED + dispatched (`composition.py` builder, exact-once facades), global `GlobalCellLedger` with the Newton-3rd force-balance gate, and `dump_composed_world_state` with global-unique filament IDs (no double-draw). Commits 43cee611 / 234963f5 / 94156518.
- **R2–R5 kernel-binding, STRUCTURAL**: 8 components authored one increment (surface/fluid/ecm/nmii/MT/IF/protrusion/linc); the two-array adjoint NMII crossbridge is bound; **all component CUDA gates PASS on the A5000 — `tests/ac/engine` = 256 passed, 0 failures** (load_path + fluid_core device bugs were test-oracle / guard bugs, fixed; kernels were correct). Honest states preserved (adversarial verify corrected 3 overstated KERNEL_BOUND claims).
- **Resting baseline (the R2 prerequisite)**: NOT closed by any solver — 3 preconditioner families + passive pre-tension + preload all REFUTED at native. Root cause characterized as MODEL COMPLETENESS: the missing pieces are a resting cortical-tension SOURCE + a force-transmitting ERM. The **candidate-① mechanism (resting bound-myosin) is BUILT + DEMONSTRATED** (`9bb91e0a`: membrane residual 40.74→2.17 pN, minimum exactly at the Laplace γ_cortex=ΔP·R/2−γ_mem). Opt-in, default-OFF; needs ONLY the PI-GAP `resting_bound_myosin_fraction + per-head force` to activate → close the gate → unblock R2 native binding.
- **Parameters**: P2 cards applied (collagen conc-gate, Hosseini γ KB-6.1.7, ERM HOLD; kb-check green). Viz: per-compartment on/off + composed-world render.

**The single blocker to resume R2→R8 is the resting-gate PI-GAP** (bound-myosin fraction/force). Everything else is done + green. See `_historical/OVERNIGHT_ORCHESTRATION_STATE_2026-07-23.md` §MORNING.
