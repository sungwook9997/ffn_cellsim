# COMPLETENESS CRITIQUE — what is missing from the 10-track plan

Ranked by damage if left out.

---

## RANK 1 — There is no coupled-solve track, and the plan's answer (per-component implicit solvers) *is* the block-diagonal failure in a new costume

**The audit claim is not just true, it is understated.** The composed native step does not "let each component relax on its own" — it does not integrate at all.

`ac/engine/composed_native.py:347-358` — `solve_candidate()` is: zero every owner's force array → `self.world.pipeline.run_candidate()` → `_load_compute()`. No position update, no iteration, no residual, no convergence latch. Its own docstring calls it "One composed candidate-force **pass**."

`ac/engine/transaction.py:129` — step 3 is `solve()`, a caller-supplied `Callable[[], None]`. `CellTransaction` owns no solver by design.

`scripts/ac_gate_b_composed_native.py:232` — `accepted_d = wp.array(np.ones(1, np.int32))`, and no `ledger`/`tol_sq_d` is passed. The composed cell's "accepted physical step" unconditionally accepts a state in which nothing moved.

`docs/v2_audit/cell_engine/NATIVE_COMPOSITION_SCOPE_2026-07-25.md` (table, §(a)) — the sf_arc driver "runs its **own** overdamped relax loop (`axpy_kernel`, `:164-240`) … Transaction: **none** — a bare relax loop."

**Connector Jacobians are required, not optional — and they already exist, but only in the merged-array form the engine forbids.** `ac/cell/implicit_mechanics.py:1718` `_stiffness` launches genuine cross-component tangent blocks: `add_erm_stiffness_kernel` (membrane↔cortex, `:1775`), `add_linc_stiffness_kernel` (nucleus↔actin, `:1770`), `add_segment_crossbridge_stiffness_kernel` (myosin↔actin, `:1764`). These work *only* because every component shares one array — `ac/cell/assemble.py:754`, "compose the global node array `[actin | myosin | nucleus | membrane]`". **De-aliasing into component-private arrays — which every track demands — destroys the only coupled solver in the tree.** No track notices this.

The one correct N-array pattern exists and is 2-way hard-coded: `ac/engine/sf_implicit.py:452-491` `SFImplicitCG.step` gathers two owners into contiguous scratch (`gather_two_kernel`, `:143`), solves the coupled operator, scatters back (`add_scatter_two_kernel`, `:173`).

What the plan proposes instead: T3 `mt_implicit.py` + `if_implicit.py`, T6 `NucleusImplicitCG`, existing `sf_implicit.py` — **N separate per-component SPD solves with connector forces exchanged between them.** That is component-level Jacobi/Gauss-Seidel; its convergence factor is set by connector-stiffness / component-stiffness, and the mechanically important connectors are precisely the stiff ones (k_erm ≈ 4600 pN/µm; α-actinin 4.6e5 pN/µm dominating λ_max ≈ 3.7e6 per `sf_implicit.py:11-13`). It will not converge for exactly the couplings that matter, and a non-converged split produces the *same* signature as the aliasing bug (39.8 nm vs 0.0 nm) — so **T2-7 / T3's transmission matrix / T6-N8 cannot distinguish "connector is dead" from "coupled iteration didn't converge."** Force-only coupling is not merely weak here; it makes every connectedness gate in six tracks unfalsifiable.

**Missing track (Track 0):** N-owner coupled implicit solve. Generalise `gather_two_kernel` to an owner index-map; require every force-real connector to supply `add_stiffness(vector, out)` + `add_preconditioner(diag)` alongside `accumulate` (make it a `ConnectorContract` runtime-API requirement so a connector cannot be registered force-real without its Jacobian); wire the predicate to `ledger.balance_ok_d` with a *declared* tolerance (`transaction.py:139-146` supports it; nothing uses it). **Negative control that doubles as the acceptance gate:** solve one 2-component system three ways — (a) merged-array monolithic, (b) N-owner coupled, (c) per-component alternating. Require (a) ≡ (b) to CG tolerance **and (c) to demonstrably FAIL**. If (c) passes, the gate is vacuous.

---

## RANK 2 — PI worked example (ii) is not covered at any link. There is no protrusion track.

Trace the chain lamellipodium → forward motion:

| Link | Owner |
|---|---|
| Protrusive force generation (lamellipodium) | **NO TRACK** |
| `lamellipodium_membrane_contact` (Brownian ratchet: polymerization → force) | **NO TRACK** — T1 lists it "explicitly OUT of scope" |
| Actin polymerization / barbed-end ratchet / treadmilling | **NO TRACK** — zero mentions across all ten |
| `lamellipodium_nascent_fa` cell side | **NO TRACK** — T10 ships ECM side only, T8 "STAGED after the protrusion track" |
| `nmii_lamellipodium_motor` (rear contraction = other half of the dipole) | **NO TRACK** — T4: "remain DECLARED-ONLY" |
| Plastic ECM remodelling | **NO TRACK** — T10's own words: "has no representation at all: `remodel_runtime` does resolution refine/coarsen, not plastic remodelling" |
| Forward motion gate | T10-S10, whose own deps concede it can "only run arms (ii)/(iii) as the null structure" |

Two tracks (T8, T10) declare a hard dependency on "the protrusion track." **It does not exist.** Yet the owners are already in-tree: `ac/engine/protrusion.py:385` `LamellipodiumStateOwner`, `:648` `FilopodiumStateOwner`, `:303` `LamellipodiumBranchAngleMechanics`, `:543` `FilopodiumFascinBundleMechanics`; plus `ff/polymerization_warp.py`, `ff/motility_warp.py`, `ff/polarization_activegel.py`, `ac/weave/lamellipodium.py:201` (debranching/treadmilling cycle). Arp2/3 angle-harmonic branching — a HARD-RULE worked example — is already launched by the incumbent (`ac/cell/driver.py:152-156`, `branch_angle_kernel` with `ARP23_THETA0_RAD`) and no track validates it.

Consequence: **the plan cannot exercise its own flagship negative control.** "Traction was isotropic radial grip that cancels laterally; F_net was one quarter of one clutch" can only be *disproven* by a run with a protrusion at the front and a motor at the rear. The plan has neither. T10's traction first-moment tensor will be measured on a cell that has no polarised channel — it will read isotropic, correctly, and prove nothing.

Also unowned: cortex actin/crosslink turnover (`ff/motility_warp.xl_turnover_kernel` via `ac/weave/crosslink_kmc.py:20`; ECM crosslink turnover *is* owned by T10-S3).

---

## RANK 3 — Connector coverage: 7 of 32 are never made force-real by any track

Union across all `connectors_involved`, against the 32 in `ac/engine/contracts.py:253-448`:

**Zero coverage (explicitly excluded by the tracks that touch their endpoints):**
1. `lamellipodium_membrane_contact` (`contracts.py:344`) — T1 "explicitly OUT of scope"
2. `filopodium_membrane_tip` (`:364`) — same
3. `nmii_lamellipodium_motor` (`:412`) — T4 "remain DECLARED-ONLY"
4. `nmii_filopodium_motor` (`:420`) — same

**Split so that no track owns the cell side (functionally holes):**
5. `lamellipodium_nascent_fa` (`:358`) — T10: ECM/coat side only; T8: staged after nonexistent track
6. `filopodium_nascent_fa` (`:391`) — same
7. `membrane_ecm_contact` (`:446`) — T10 supplies ECM side and says "the membrane-side quadrature depends on the membrane track"; T1 defers past T1.9; T8 calls it "adjacent, not closed by this track". **Nobody builds the membrane-side quadrature** — and the PI baseline is an adherent cell resting on collagen.

**Ambiguous:** `dorsal_arc_crosslink` (`:428`) — T3 says "must either give it the full API or reclassify"; T4 excludes it. An unresolved either/or is not ownership.

**Correctly flagged-but-unowned PI decisions:** no cell↔`world_boundary` edge exists in the 32 (only `ecm_far_field_anchor`, `:439`). T10 flags the coverslip; T2 flags the probe as a possible 33rd. Neither track owns getting the decision made, and an adherent baseline needs the cell to rest on *something*.

---

## RANK 4 — Cortex-local arrays and the `_accumulate_all` omit-mask: the two most-cited dependencies, owned by nobody

- **Omit-mask / driver seam:** declared a HARD blocker by T1, T3, T6, T7. `ac/cell/driver.py:141` `_accumulate_all(cell, pos, f)` has no `omit` parameter; the only "omit" in the file is build-time CLI (`:1288-1312 --no-myosin/--no-nucleus/--no-membrane`), which cannot hand one component to the engine while keeping the rest. Four tracks blocked, zero owners.
- **Cortex-private arrays:** T1 lands *membrane*-local only; T2 explicitly declines ("Track 2 does not require it"); T3 calls it a "HARD BLOCKER" assigned to "cortex/membrane track"; T6 says "whichever track lands owner-private arrays… should land the same mechanism." `INTEGRATION_surface_body.md:22-39` lists five un-done shared-file steps assigned to "the Lead," not to any track.
- **Why this is worse than it looks:** 7 of the 32 connectors terminate on the cortex (`membrane_erm_cortex`, `surface_porous_transfer`, `sf_cortex_transient`, `mt_cortex_capture`, `lamellipodium_cortex_seam`, `filopodium_cortex_root`, `nmii_cortex_motor`). While the cortex is a bind-target *port* over `cell.pos_d`/`cell.f_d`, all seven adjoint claims are structurally unverifiable — and T3 itself says its gates "will keep reporting the aliasing until it lands." **No track promises to make the largest component a component.**
- `ac/engine/sf_implicit.py` is untracked in git today. T4-M5 requires it committed. No track owns commits.

---

## RANK 5 — Instrumentation must be its own track, landed first; otherwise the plan silently stalls *and* diverges

~118 missing observables across ten tracks, with heavy duplication of the *same* quantity under different names:

| Observable | Tracks specifying it independently |
|---|---|
| Per-connector two-sided transmitted resultant + adjoint-work closure | T1, T3, T5, T6, T7, T8, T10 (**seven**) |
| Per-component rest-vs-load displacement table / transmission matrix | T2-7, T3 `transmission_matrix.py`, T6-N8, T7-V8 (**four separate implementations**) |
| Net vs gross (Σf vs Σ\|f\|) non-vacuity denominator | T3, T8, T10 |
| Device enclosed volume | T1, T7 (both note the reduce kernel already exists — `ac/nucleus/envelope.py:177-183`) |
| Machine-measured evidence rung + BLOCKED stamp + JSON artifact | T1, T3 ("audit critical-path item 0"), T4, T8, T10, T11 |
| Measured noise floor / null harness | T1, T6 |

If each track builds its own, cross-track comparison — which is the entire content of "compartments deform each other" — becomes impossible, and four incompatible transmission matrices will disagree with no arbiter. `ac/engine/forces_manifest.py` already exists (49 KB) and T3 reports **no driver calls it**, so the SettlingForce-recurrence gate credited as closed still does not run.

**Missing track (Track 0b):** one `ac/engine/observe/` package (transmission matrix, per-connector work ledger, Σf/Σ\|f\|, noise-floor harness, rigid-mode projector, device volume/area reducers), one artifact schema with machine-measured rung + provenance + BLOCKED stamps, one negative-control harness. Land before T1-T10 physics. Everything else in the plan reduces to a specialisation of these.

---

## RANK 6 — Order and dependency cycles

**Genuine cycles:**
- **T8 ⇄ T10.** T8 needs "ECM component at CUDA_UNIT with a quasi-2D collagen-coated dish config" (= T10 S1-S3); T10 needs "focal_adhesion component owner MUST land first (or concurrently)" (= T8). Break: T10 S1-S5 → T8 A0-A2 → merged A3/S6.
- **File collision inside that cycle:** both tracks create `ac/engine/fa_clutch_connector.py` with different internals (T8 wraps `LoadPathJointRuntime`; T10 builds on `load_path.segment_pair_reference` + `_segment_joint_force_kernel`). One file, two designs — the shared-tree divergence failure mode.

**Blocked on a deliverable no track produces:**
- **T5-S4 — "THE DECISION GATE," and T5 says "nothing beyond S4 should be committed before the PI sees the S4 number."** Its dep requires a *converged cortex resting baseline*. T2 explicitly does not deliver one (its linear gates "do not require a converged resting state"; T2-6 capped at ~200 nm). GATE A is open on a PI modelling decision. So the plan's single most consequential decision gate is unreachable.
- **T11-S1/S2/S5** wait on "the PI-gated smoothness-preserving overlap fix" (`docs/v2_audit/CORTEX_KINK_AUDIT_2026-07-25.md`, 3.42 → 0.13). Unowned.
- **T2's load instrument** needs an additive hook in frozen `ac/cell/driver.py:141-167`. If PI refuses and the omit-mask isn't landed, T2 cannot apply a load at all — and worked example (i) becomes undeliverable.

**Ordering contradiction:** T2-9, T7-V9 and T10-S11 each need the ~40× default-off wins *before* their sweeps; T11 states "everything from S1 onward depends on the assembly tracks landing first." Fix: split T11 into **T11a** (S0 profiler + FIX1/FIX2 bit-identity A/B on the current frozen config, land now) and **T11b** (rest, post-assembly).

**Purpose inversion:** PI purpose #5 (infer molecular parameters by sweeping against macroscopic data) lives entirely in **T11-S8 — the last step of the last track**, which self-admits it "cannot bootstrap itself." Truncate the plan anywhere and the project's stated purpose is the first casualty.

---

## RANK 7 — PI worked example (i): covered, with two unowned gates in the path

**Delivered end to end by T2:** T2-0 absolute connectivity → T2-1 exact Green's-function reciprocity vs dense NumPy oracle → T2-2 native reciprocity + transmission at 494,802 nodes → T2-3 decay field, functional-form discrimination, mesh independence → T2-5 cut-plane transmitted force with per-family decomposition (the sums-to-zero catcher). This is the strongest part of the plan and correctly needs no magnitude.

**Two blockers, both unowned:** (a) the frozen-driver load hook above; (b) the fine-mesh half — T2-6 is capped at ~200 nm pending multigrid, and the plateau that motivated multigrid was an overlap_free WCA kink whose PI-gated fix no track owns. Verdict: **(i) ≈ 80% covered, reachable once those two land. (ii) not covered at any link.**

---

## RANK 8 — What no track owns (residual list)

- `focal_adhesion` as a registered **component state owner** (`contracts.py:219`, `owns_geometry=False`). T8 builds a connector + a site population; T10 says the owner "MUST land first." Neither says "I build the component."
- **The 13/32 composition gate itself.** Every composed build today passes `require_complete=False` (`composed_native.py`, "bring-up composed world … not the whole 13/32 cell"). Nobody schedules `require_complete=True`, so "everything is assembled" — T11's own precondition — has no owner and no gate.
- **`CellState` axis.** `ac/engine/cell_state.py:22`: "ratified instance ships with EMPTY registries." PI purpose #5 is inference per cell type **and per cell STATE**; no track populates or validates the state axis, and `ff/cell_type.py` (KB-grounded presets) is in the reuse list unwired.
- **Identifiability.** Nothing checks that the observable set the plan produces can even *identify* the molecular parameters. A local Fisher-rank / profile-likelihood check on each track's observable vector is cheap and is the actual gate for PI purpose #5.
- The 4 pending CHEMICAL_FLUX connectors (32→36, PI-gated): `ac/engine/monomer_flux.py` exists; T7 uses it only as a *pattern*.
- `ff/piezo.py` (mechanotransduction): one passing mention in T7, no owner.
- `scripts/ac_composed_world_dump.py` `TypeError` (missing `propose_events` on `membrane_erm_cortex`): T3 correctly flags it as blocking every cumulative render — i.e. the per-stage viz gate is currently unsatisfiable for **all ten tracks**, and only T3 lists the fix.