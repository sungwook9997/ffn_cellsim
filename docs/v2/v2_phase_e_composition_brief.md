# V2 Closed-Loop ECM Gate — Phase E Composition — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-05 KST
**Status**: design-discussion brief, **opening position only —
NOT a lock, NOT a code commitment**.
**Author**: implementation-work Claude, drafted under (B-trigger)
precedent established by HB#3 / HB#4 / Phase D no-op / B1 / HB#1+#2
/ HB#5 (all locked + Sanity Gated + implemented + Codex PASS in
this session).
**Source**: locked phased plan
`docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md` §1 Phase E
+ Codex `id=1522` HB#5 final PASS routing directive ("move the
next unit to Phase E composition planning").
**Unblock condition** (now satisfied): all 5 Phase E upstream
hard blockers locked + impl PASS:
- HB#1+#2 (constitutive direction + saturation merged):
  `f05abef` (Codex `id=1495` PASS).
- HB#3 (FA→ECM scattering): pre-session lock + commits.
- HB#4 (ECM→FA bias target neutral): pre-session lock + commits.
- HB#5 (Lyapunov-like metric): `7b2a186` (Codex `id=1522` PASS).
- Phase D no-op orchestrator: `60f726f` + `b1b2b92`.

**Hard contract**: this brief is the **opening position** for
the adversarial design round on the Phase E composition cycle.
It does NOT lock the composition shape, does NOT commit the
implementation, and does NOT authorize closed-loop ECM gate
satisfaction wording. The locked phased plan §1 closed-loop
ECM gate items 1–6 require explicit per-item satisfaction
evidence; this brief frames how Phase E composition tests them.

---

## 0. Why this brief exists + Phase E v1 vs v2 scope boundary

Phase E composition is the **first user** of all 5 Phase E
upstream hard blockers + Phase D no-op orchestrator. It is the
unit that demonstrates the closed-loop ECM gate is functional
under the locked-orientation-only first active law (HB#1+#2)
+ neutral ECM→FA bias (HB#4 Phase D default).

**Phase E scope boundary** (per HB#5 Y4 + HB#4 Phase D-default
locked plan §3):

- **Phase E v1 (this brief's scope)**: combines Phase D wrapper
  + HB#1+#2 active orientation update + HB#4 still-neutral bias
  + HB#5 metric. Provides **ECM-orientation half** evidence for
  the closed-loop ECM gate (Items 1, 2, 3, 4 per the gate
  matrix). Item 5 (sensitivity sweep) is a separate sweep
  harness unit. Item 6 (Rule 10 unit-chain) is already
  proven inline in each HB lock + Sanity Gate.
- **Phase E v2 (separate later cycle)**: replaces neutral HB#4
  with an active variant (`compute_ecm_to_fa_bias_active`
  function, separate Sanity Gate + lock). Demonstrates **full
  FA→ECM→FA closed loop**. Per HB#5 Y4 explicit boundary, this
  is OUT of HB#5 scope and OUT of Phase E v1 scope.

The Phase E v1 vs v2 boundary is the **central design decision**
of this brief: should v1 land first as a stand-alone milestone,
or should v1 be delayed until HB#4-active is also designed
(merged as v1+v2)?

---

## 1. What is already locked (do not redebate here)

Per locked phased plan + completed unit locks:

- **Phase D no-op orchestrator** (`docs/v2/v2_phase_d_no_op_scaffolding_locked.md`,
  `60f726f` + `b1b2b92`): `step_phase_d_no_op(adhesions, ecm)`
  returns `PhaseDNoOpStepResult` with `fa_to_ecm` (scatter +
  identity ECM) + `ecm_to_fa` (neutral bias) + `updated_ecm`
  (object identity in Phase D). Deterministic call order:
  HB#3 scatter first → HB#4 bias only if HB#3 succeeds.
- **HB#1+#2 active orientation update**
  (`docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`,
  `f05abef`): `step_ecm_orientation_response(ecm,
  traction_density_xy, dt_s, *, traction_ref_nN_per_um2,
  k_orient_per_s)` returns `ECMOrientationResponseResult` with
  updated ECM + 13-field diagnostics. Two-sided validation
  contract (`ecm.validate()` at entry +
  `updated_ecm.validate()` at return) per Codex `id=1486`.
  5 owned `failure_kind`s via `ECMConstitutiveResponseError`.
- **HB#3 FA→ECM scatter** (sister-precedent locked + impl):
  `scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)` returns
  `(nx, ny, 2)` traction density.
- **HB#4 ECM→FA bias** (sister-precedent locked + impl):
  `compute_ecm_to_fa_bias_neutral(adhesions, ecm)` returns
  `ECMToFABiasResult` with all-1.0 multipliers (Phase D default).
- **HB#5 Lyapunov-like metric**
  (`docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`,
  `7b2a186`): `compute_ecm_orientation_lyapunov_metric(ecm,
  traction_density_xy)` returns
  `ECMOrientationLyapunovMetricResult` with `v_active_um2` +
  10-field diagnostics. Pure read-only, active-cells-only V.
- **§3 effective-stiffness guard**: NO `effective_stiffness()`
  helper consuming `fiber_density` or `orientation_tensor` in
  Phase E first law.
- **Closed-loop ECM gate Item 6 (Rule 10)** is already proven
  inline in each HB lock + Sanity Gate.

This brief is the next-level-down design choice **inside** that
already-locked frame.

---

## 2. The four design questions

### Q1 — Module / function shape

Phase D lock §1 instructs separate module
`acs/v2/dynamics/closed_loop_phase_e.py` with its own lock.
Three candidate function shapes:

- **(a) Single composition function**:
  `step_closed_loop_phase_e(adhesions, ecm, dt_s, *,
  traction_ref_nN_per_um2, k_orient_per_s) ->
  PhaseEStepResult` — combines all in one call. Caller is
  responsible for FA-side updates (6.3a/6.3b) before/after the
  Phase E ECM-side step.
- **(b) Pluggable orchestrator**:
  `step_closed_loop(adhesions, ecm, dt_s, *, ecm_response_law,
  ecm_to_fa_bias_law, ...)` — function takes the response and
  bias laws as parameters. v1 calls with neutral bias; v2 with
  active. **Rejected by HB#1+#2 lock §1 forbidden list**:
  optional Phase E params with default values are silent-
  activation risk. Same risk applies here.
- **(c) Two functions**: `step_closed_loop_phase_e_v1` (with
  neutral bias hard-wired) + future `step_closed_loop_phase_e_v2`
  (with active bias hard-wired). Mirrors HB#1+#2 / HB#4 Phase D
  / Phase E function naming separation precedent.

### Q2 — Composition return shape

Once the function is called, what does it return?

- **(a) Single composite dataclass** `PhaseEStepResult` with
  fields:
  - `traction_density_xy: np.ndarray (nx, ny, 2)` (HB#3 output)
  - `orientation_response: ECMOrientationResponseResult`
    (HB#1+#2 output: updated_ecm + 13-field diagnostics)
  - `ecm_to_fa_bias: ECMToFABiasResult` (HB#4 neutral output)
  - `lyapunov_metric: ECMOrientationLyapunovMetricResult`
    (HB#5 output: v_active_um2 + 10-field diagnostics)
  - `updated_ecm: ECMSubstrateState` (= `orientation_response.updated_ecm`)
- **(b) Tuple of sub-results**: `(traction_density_xy,
  orientation_response, ecm_to_fa_bias, lyapunov_metric,
  updated_ecm)` — caller unpacks.
- **(c) Minimal**: just `(updated_ecm, lyapunov_metric_v_active)`.
  Loss of intermediate state.

Sister-pattern: Phase D `PhaseDNoOpStepResult` chose (a)
composite dataclass with explicit fields (no diagnostics dict).
HB#1+#2 chose (a). HB#5 chose (a). Strong sister-pattern for (a).

### Q3 — Composition order + dependency map

What is the canonical order of sub-call within
`step_closed_loop_phase_e`?

Logical order (per closed-loop ECM gate semantics):

1. HB#3 scatter: `traction_density_xy =
   scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)` —
   reads FA state + current ECM (geometry only).
2. HB#1+#2 active update:
   `orientation_response = step_ecm_orientation_response(ecm,
   traction_density_xy, dt_s, ...)` — reads current ECM +
   scatter, produces updated ECM with new orientation.
3. HB#4 bias readout:
   `ecm_to_fa_bias = compute_ecm_to_fa_bias_neutral(adhesions,
   updated_ecm)` — reads updated ECM (post-HB#1+#2). Per Phase D
   no-op, multipliers are all 1.0; readout is for
   diagnostic/sweep purposes.
4. HB#5 metric:
   `lyapunov_metric = compute_ecm_orientation_lyapunov_metric(
   updated_ecm, traction_density_xy)` — reads updated ECM +
   scatter, returns scalar V_active + diagnostics.

**Open question**: HB#4 reads `updated_ecm` (post-HB#1+#2) or
`ecm` (pre-HB#1+#2)? Per Phase D no-op semantics, neutral
multipliers are state-independent — same answer either way.
But for forward-compatibility with HB#4-active in Phase E v2,
the bias readout should match what the next FA dynamics step
would see (post-HB#1+#2 ECM).

**Validation order**: Each sub-call already validates its own
inputs; the composition wrapper does NOT add new validation.
HB#1+#2 calls `ecm.validate()` at entry; HB#5 calls it again at
its entry (after HB#1+#2 produces updated_ecm). Two
`ecm.validate()` calls per Phase E step is acceptable
(microsecond-scale check for typical grid sizes); if
performance becomes an issue (Phase F sweep), batched
validation is a separate optimization.

### Q4 — Closed-loop ECM gate items 1-6 testing scope

Per locked phased plan §1: Items 1-6 must be satisfied for the
gate. Phase E composition tests **which** items?

| Item | Description | Phase E v1 responsibility |
|---|---|---|
| 1 | Response monotonicity under prescribed traction | **Yes** — test under fixed traction |
| 2 | Saturation behavior under repeated traction | **Yes** — test under sustained traction (convex saturation) |
| 3 | No response when traction = 0 | **Yes** — test (HB#1+#2 `w=0`, HB#5 V=0) |
| 4 | Bounded feedback in single-cell loop | **Yes** — test (HB#5 boundedness channel) |
| 5 | Sensitivity sweep over grid spacing + dt_ecm | **No** — separate sweep harness unit |
| 6 | Rule 10 unit-chain proof | **No** — already proven inline in each HB lock |

So Phase E v1 tests Items 1-4 (4 of 6); Items 5 + 6 are
out-of-scope for this brief.

**Open question**: should Phase E v1's test catalog include
explicit `test_phase_e_satisfies_item_X` tests for each of
Items 1-4, or are the items implicit in the existing HB tests?

My lean: explicit `test_phase_e_satisfies_item_X` tests, one
per gate item, that asserts the GATE-LEVEL property under
composition (not just per-HB). This makes the gate satisfaction
auditable post-Phase E commit.

---

## 3. Recommended opening positions (NOT a lock)

| Q | Recommendation | Rationale |
|---|---|---|
| Q1 module / function | (c) `step_closed_loop_phase_e_v1` hard-wires neutral bias; future v2 separate function | Mirrors HB#1+#2 + HB#4 Phase D / Phase E function naming separation precedent; rejects (b) silent-activation risk |
| Q2 return shape | (a) `PhaseEStepResult` composite dataclass with 5 fields (traction_density_xy + orientation_response + ecm_to_fa_bias + lyapunov_metric + updated_ecm) | Sister-pattern with Phase D / HB#1+#2 / HB#5 (all chose composite dataclass) |
| Q3 composition order | HB#3 → HB#1+#2 → HB#4 (post-HB#1+#2 ECM) → HB#5 (post-HB#1+#2 ECM) | Forward-compat with HB#4-active in v2; matches what next FA dynamics step would see |
| Q4 testing scope | Phase E v1 tests Items 1-4 explicitly; Items 5 + 6 deferred (Item 5 separate sweep, Item 6 already proven) | Explicit test names = auditable; separate sweep harness avoids scope creep into Phase E v1 |

**Caveats / unresolved-by-this-brief**:

- Phase E v1 with HB#4-neutral does NOT close the FA loop (per
  HB#5 Y4 explicit boundary). The "closed-loop ECM gate" wording
  for v1 must be careful: it satisfies the **ECM-side** Items
  1-4 only, not the full FA→ECM→FA loop. Round 1 should lock
  the satisfaction-wording boundary.
- Item 4 (bounded feedback) under HB#4-neutral only tests
  ECM→ECM (orientation update under prescribed traction) — NOT
  ECM→FA→ECM. The "single-cell loop" wording in Item 4 is
  potentially misleading without HB#4-active. Round 1 should
  decide whether to gate-pass Item 4 under v1 or defer to v2.
- The composition function does NOT update FA state. Caller is
  responsible for 6.3a/6.3b FA dynamics before/after Phase E.
  This is the "ECM-side step" framing — explicit in the
  function name + docstring.

---

## 4. Allowed Phase E composition candidates vs Forbidden shortcuts

### Allowed (Phase E v1 candidates)

- `step_closed_loop_phase_e_v1(adhesions, ecm, dt_s, *,
  traction_ref_nN_per_um2, k_orient_per_s) -> PhaseEStepResult`:
  pure orchestrator, calls HB#3 scatter + HB#1+#2 update + HB#4
  neutral bias + HB#5 metric in locked order, returns composite
  dataclass.
- `PhaseEStepResult` typed dataclass (frozen+slots, 5 fields).
- Module path `acs/v2/dynamics/closed_loop_phase_e.py` per Phase
  D lock §1 instruction.
- Module docstring forbidden-list guard per B1 `id=1458`
  precedent.

### Forbidden (per Hard Rules + sister patterns + Phase E v1
scope boundary)

- **HB#4-active variant**
  (`compute_ecm_to_fa_bias_active`) — DEFERRED to Phase E v2;
  Phase E v1 uses neutral only. Optional `bias_law` parameter
  is rejected (silent-activation risk per HB#1+#2 lock).
- **Mutation of `adhesions`** — Phase E v1 does not update FA
  state; caller owns 6.3a/6.3b dynamics.
- **Mutation of `ecm` in place** — HB#1+#2 returns a new
  `updated_ecm`; Phase E composition uses that (no in-place
  mutation).
- **`effective_stiffness()` shortcut** — phased plan §3 guard.
- **Empirical tuning constants at composition layer** —
  Magic-Number Block; Phase E v1 introduces zero new constants
  (all HB constants inherited).
- **PI-data fitting** — Hard Rule 1.
- **"Full FA→ECM→FA closed loop" wording without HB#4-active**
  — HB#5 Y4 explicit boundary; v1 satisfies ECM-side Items 1-4
  only.
- **Item 5 (sweep) + Item 6 (Rule 10 inline) tests at Phase E
  layer** — Item 5 is separate sweep harness unit; Item 6
  already proven inline in each HB lock.
- **Generic[...] parameterization or preemptive `# type: ignore`**
  (B1 sister-precedent).
- **Pluggable strategy / dict diagnostics** (Y5 sister-pattern).

---

## 5. Reference biology / control table

Phase E composition is mostly orchestration; the biology is
embedded in HB#1+#2 (literature-pinned in HB#1+#2 lock §1:
Munevar 2001 + Tan 2003 for `TRACTION_REF`, Hall 2016 +
Notbohm 2015 for `TAU_ALIGN_RANGE_S`). Phase E composition adds
no new literature parameters.

The closed-loop ECM gate Items 1-6 are the discipline source
(per `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`)
rather than direct biology references.

---

## 6. What is *out of scope* for this brief

- HB#4-active variant
  (`compute_ecm_to_fa_bias_active`) — Phase E v2 unit.
- FA dynamics integration at Phase E layer — caller's
  responsibility (6.3a/6.3b precedent).
- `fiber_density` Lyapunov component — HB#1+#2 didn't update
  `fiber_density`; future Phase F constitutive law that does
  is a separate unit.
- `stiffness_kpa` Lyapunov component — phased plan §3 guard.
- Cell-substrate work rate (phased plan §5 candidate 3).
- Item 5 sensitivity sweep — separate sweep harness unit.
- Item 6 Rule 10 unit-chain — already proven inline in each
  HB lock.
- PI-data parameter fitting — Hard Rule 1.
- Full closed-loop FA→ECM→FA stability evidence — HB#5 Y4
  explicit boundary; needs HB#4-active.

---

## 7. Sanity Gate items the future Phase E v1 lock+code will
need

Per CLAUDE.md Sanity Gate Protocol applied to a Phase E
composition orchestrator:

- **§1 Dimensional**: Phase E composition introduces no new
  unit chain; inherits HB#1+#2 + HB#3 + HB#4 + HB#5 chains.
  Rule 10 verifies via composition at the dataclass-field
  level (each field's units locked at the source HB).
- **§2 Boundary**: empty FA list → all-zero scatter (HB#3) →
  HB#1+#2 identity update → HB#4 (0, 3) multipliers → HB#5
  V_active = 0 → composite dataclass with all-zero/identity
  fields. dt=0 → HB#1+#2 identity update; HB#5 reads unchanged
  ECM.
- **§3 Conservation**: HB#1+#2 has Y2 distance-to-target
  monotone per fixed target; HB#5 boundedness channel. Phase
  E composition under fixed traction direction → V_after ≤
  V_before per step (composition of HB#1+#2 update + HB#5
  read). Under varying traction direction → V_active ≤ V_max
  bound per step.
- **§4 Numerical**: float64 throughout (inherited); no new
  tolerance constants; per-call work = HB#3 + HB#1+#2 + HB#4
  + HB#5 sequential.
- **§5 Sign**: composition introduces no sign; inherited from
  sub-calls.
- **§6 Measurement-protocol**: Phase E composite dataclass
  exposes all 4 sub-results — no scalarization at composition
  layer. Hard Rule 11 boundary preserved (HB#5 active-cells-only
  V; HB#1+#2 update-domain matching).

**Closed-loop ECM gate items satisfaction tests** (Q4):

- `test_phase_e_satisfies_item_1_response_monotonicity`:
  prescribed-traction step yields measurable orientation
  response (V_active changes monotonically toward V_active=0).
- `test_phase_e_satisfies_item_2_saturation`: sustained traction
  → orientation reaches `T_target` asymptotically (HB#1+#2
  convex-weight saturation).
- `test_phase_e_satisfies_item_3_no_response_when_traction_zero`:
  zero traction → HB#1+#2 `w=0` everywhere → identity update;
  HB#5 V_active = 0; ECM bytewise unchanged.
- `test_phase_e_satisfies_item_4_bounded_feedback_single_cell_loop`:
  varying-traction sequence → V_active ≤ V_max bound at every
  step. **Wording caution per §3 caveat**: this is ECM-side
  bounded feedback, not full FA→ECM→FA loop (HB#5 Y4).

**Magic-Number Block**: zero new tunables at composition
layer. Inherits HB#1+#2 + HB#5 constants.

---

## 8. What this brief is *not*

- Not a lock. The composition shape + Item 4 satisfaction
  wording are decided by the design-discussion adversarial
  round.
- Not a Sanity Gate doc.
- Not a closed-loop ECM gate v2 satisfaction. v1 satisfies
  Items 1-4 (ECM-side); v2 satisfies the full closed loop
  (separate cycle).
- Not authoritative for biological parameters (none new).
- Not authoritative for the Item 5 sweep harness — that is a
  separate unit.

---

## 9. Process expectation (mirrors HB#3/#4/Phase D/B1/HB#1+#2/HB#5
precedents)

1. **Cross-room dispatch** to design-discussion. Status:
   `decision-needed`.
2. **Adversarial round**. HB#3 4 / HB#4 4 / Phase D 2 / B1 4 /
   HB#1+#2 4 / HB#5 3 rounds. Phase E composition has narrow
   scope (orchestrator only; all HBs locked) but item-1-4
   satisfaction wording + v1/v2 boundary clarity is
   substantive. **2-3 rounds plausible**.
3. **Lock artifact** delivered as
   `docs/v2/v2_phase_e_composition_locked.md`.
4. **Cross-room dispatch back to implementation-work**.
5. **Sanity Gate doc** drafted by impl-work, Codex review.
6. **Code commit** — new module
   `acs/v2/dynamics/closed_loop_phase_e.py` implementing
   `step_closed_loop_phase_e_v1` + `PhaseEStepResult` dataclass
   + module-docstring forbidden-list guard; new test catalog
   including item-1-4 satisfaction tests; exports through both
   `__init__.py` (2 new symbols if Q2=a + Q1=c: function +
   dataclass).
7. **Codex review** of code commit.
8. **Phase E v1 unit complete** = closed-loop ECM gate Items
   1-4 satisfied (ECM-side); Item 5 (sweep) + HB#4-active +
   Phase E v2 are next.

implementation-work stays idle/review-capable while this design
round runs.

---

## 10. References

- Locked phased plan:
  `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md` §1 Phase
  E + §3 effective-stiffness guard + §6 Sanity Gate matrix
  Phase E.
- Forward roadmap: `docs/v2/v2_phase1_forward_roadmap.md` (commit
  `da40eef`).
- Phase D lock + impl (closed_loop_phase_e.py module-naming
  precedent):
  - `docs/v2/v2_phase_d_no_op_scaffolding_locked.md`
  - `acs/v2/dynamics/closed_loop_phase_d.py`
- HB#1+#2 lock + impl (active orientation update):
  - `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`
  - `acs/v2/dynamics/ecm_constitutive_response.py`
- HB#3 lock + impl (FA→ECM scatter):
  - `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `acs/v2/dynamics/fa_to_ecm_scattering.py`
- HB#4 lock + impl (ECM→FA bias neutral):
  - `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `acs/v2/dynamics/ecm_to_fa_bias.py`
- HB#5 lock + impl (Lyapunov-like metric):
  - `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`
  - `acs/v2/dynamics/ecm_lyapunov_metric.py`
- B1 audit fix (typed FocalAdhesionDynamicsResult):
  - `docs/v2/v2_focal_adhesion_dynamics_result_typed_locked.md`
  - `acs/v2/dynamics/focal_adhesion.py`
- Memory rules:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, per Codex `id=1428`).
  - `rule10_unit_derivation_in_docs.md`.
  - `hard_rule_11_wording_boundary_meta_test.md`.
  - `feedback_aggressive_design_debate.md`.
  - `cadence_promise_must_send_even_when_idle.md`.

---

## 11. Cross-room dispatch instruction (impl-work →
design-discussion)

The MCP message that accompanies this brief:
- `to`: `codex`
- `room`: `design-discussion`
- `topic`: `v2-phase-e-composition`
- `status`: `decision-needed`
- `body`: brief 4-question summary + Phase E v1 vs v2 boundary
  + recommended opening positions
- `refs`: implementation-work `id=1522` (HB#5 PASS routing
  directive) + locked phased plan §1 + 5 HB locks + Phase D
  lock as composition prerequisites

design-discussion round produces a lock artifact that supersedes
this brief; this file remains as historical opening-position
context.
