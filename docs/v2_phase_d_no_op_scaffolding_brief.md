# Phase D No-Op Scaffolding Integrator — Design-Discussion Brief (opening, NOT a lock)

**Date**: 2026-05-04 KST
**Status**: design-discussion brief, **opening position only — NOT a
lock, NOT a code commitment**.
**Author**: implementation-work Claude, drafted under (B-trigger)
precedent established by Hard Blocker #3 / #4 (both locked + Sanity
Gated + implemented + Codex PASS in this session). Codex impl
`id=1389` ack with five brief emphases (see §0).
**Source**: closed-loop ECM gate phased plan
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §1 Phase D +
§2 Phase D entry contract, plus the now-completed Hard Blocker #3
(`docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`) and
Hard Blocker #4
(`docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`) lock
artifacts.

**Hard contract**: this brief is the **opening position** for the
adversarial design round on Phase D no-op scaffolding. It does not
lock the integrator interface, does not commit the implementation,
and does not authorize any Phase E activation. The locked phased
plan §1 Phase E remains blocked until all 5 Hard Blockers + the
effective_stiffness law decision are resolved.

---

## 0. Why this brief exists + Codex 5 emphases + Step 6

The Phase D entry blockers (HB#3 + HB#4) are now both locked +
Sanity-Gated + implemented + Codex-PASSed. The phased plan §1
Phase D scope is "schema/API hooks for FA→ECM and ECM→FA paths,
default OFF or identity, no behavior change" plus an integrator
wrapper "cell step → 6.3b FA step → step_fa_to_ecm_response →
step_ecm_to_fa_bias → next step". The named functions
`step_fa_to_ecm_response` and `step_ecm_to_fa_bias` are API
commitments that have not yet been written; this brief opens
their design lock.

**Per Codex impl `id=1389`, five brief emphases enforced
throughout this brief**:

1. Phase D default identity ECM update, no active response law.
2. Explicit return dataclass (or named tuple) exposing
   `updated_ecm`, `traction_density_xy`, and the
   `ECMToFABiasResult` — no hidden measurement modality at the
   integrator layer.
3. Phase E hook **structurally absent** in Phase D, OR explicitly
   no-op. Avoid optional `mapping_law` / `response_law` params
   that silently activate behavior.
4. Integrator must not mutate input ECM in place unless the
   adversarial round explicitly chooses that semantics.
5. HB#3 / HB#4 validation ordering as sister-patterns to preserve
   (lock-specific failure_kind first → schema `fa.validate()` →
   out-of-grid raise → deposit/sample).

**Step 6 sister-gate-mirror check (baked into pre-commit
checklist for this brief; NOT yet promoted to the shared memory
rule per Codex `id=1389`)**: when an integrator copies an
interface or validation pattern from HB#3 / HB#4, verify
line-by-line that the sister-module pattern is mirrored —
particularly the validation order
(`len(position) check → unpack → non_finite_fa_position →
fa.validate() → out-of-grid → deposit/sample`).

---

## 1. What is already locked (do not redebate here)

Per `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`:

- §1 Phase D scope: schema/API hooks, identity / neutral default,
  no behavior change.
- §2 Phase D entry contract:
  - `step_fa_to_ecm_response` named API: `(adhesions, ecm_state,
    scattering_params)` → `traction_density_nN_per_um2[grid_shape]`
    (per-cell vector) per locked plan §2 (mirroring the HB#3
    output type).
  - `step_ecm_to_fa_bias` named API: `(adhesions, ecm_state,
    bias_params)` → per-FA bias multipliers (matching the
    `ECMToFABiasResult` shape from HB#4).
  - Default no-op: zero traction density / 1.0 multipliers.
- §3 effective_stiffness side-door guard: any helper consuming
  ECM geometry fields is ECM→FA bias by another name; the Phase
  D integrator must not introduce such a side-door beyond
  HB#4's locked sampler.

Per
`docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`:
- `scatter_fa_traction_to_ecm_bilinear(adhesions, ecm) ->
  np.ndarray (nx, ny, 2) nN/μm²` is the FA→ECM scatter
  primitive. Pure function, no ECM mutation, no scalarization.

Per
`docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`:
- `compute_ecm_to_fa_bias_neutral(adhesions, ecm) ->
  ECMToFABiasResult` is the ECM→FA bias primitive (Phase D
  no-op default). `sample_ecm_at_fa_positions(adhesions, ecm) ->
  ECMSampledAtFAs` is the diagnostic sampler.
- `RATE_NAMES = ("k_maturity_per_s", "k_bind_per_s",
  "k_unbind_per_s")` fixes column ordering.
- Phase D / Phase E function naming separation
  (`compute_ecm_to_fa_bias_neutral` vs the future
  `compute_ecm_to_fa_bias_active`) is the silent-activation
  guard.

This brief is the next-level-down design choice **inside** that
already-locked frame.

---

## 2. The four design questions

### Q1 — `step_fa_to_ecm_response` return type

The locked phased plan §2 says the function returns
`traction_density_nN_per_um2[grid_shape]`. But the function is
the **whole** Phase D step: scatter + (Phase D no-op identity
response) → ECM. Two candidate return types:

- **(a) Locked plan literal**: return only the
  `traction_density_xy` array `(nx, ny, 2)`. Caller is
  responsible for any subsequent ECM update (Phase E response
  law). Phase D no-op = caller does nothing.
- **(b) Composite dataclass**: return a `PhaseDStepResult`
  dataclass with `traction_density_xy: np.ndarray (nx, ny, 2)`
  + `updated_ecm: ECMSubstrateState` (= input `ecm` in Phase D,
  identity) + `bias_result: ECMToFABiasResult` (one call to
  `compute_ecm_to_fa_bias_neutral`). All Phase D-relevant
  quantities exposed in one return.
- **(c) Hybrid**: return dataclass with only the scatter
  output + the bias result (omit `updated_ecm` since it's
  trivially identical to input in Phase D).

Codex emphasis #2 prefers (b) for its explicit measurement
modality. The phased plan §2's literal return type is
reachable from (b) via the dataclass field; nothing structural
breaks.

### Q2 — `step_ecm_to_fa_bias` alias / wrapper / different function

`compute_ecm_to_fa_bias_neutral` is already the Phase D no-op
default per the HB#4 lock. The phased plan §2 names
`step_ecm_to_fa_bias` as the Phase D API surface. Three
candidates:

- **(a) Alias**: `step_ecm_to_fa_bias = compute_ecm_to_fa_bias_neutral`
  module-level rebinding. Caller uses the `step_*` name; HB#4
  ownership preserved.
- **(b) Thin wrapper**: `def step_ecm_to_fa_bias(adhesions, ecm,
  bias_params=None): return compute_ecm_to_fa_bias_neutral(adhesions, ecm)`
  — adds a `bias_params` slot for future Phase E. **Codex
  emphasis #3 challenges this** because optional `bias_params`
  with default `None` silently routes to neutral, but a Phase E
  caller could pass an active law without changing function name.
- **(c) Different function**: keep `step_ecm_to_fa_bias` as a
  separate named API that calls `compute_ecm_to_fa_bias_neutral`
  internally with no extra params. Phase E would add a
  `step_ecm_to_fa_bias_active` separate function. Mirrors HB#4
  Phase D / Phase E function naming separation precedent.

Recommended: (c) — preserves the structural Phase D / Phase E
guard from HB#4. (a) is acceptable but hides the future
expansion path.

### Q3 — Phase E hook structural shape

Codex emphasis #3: "Phase E hook must be structurally absent or
explicitly no-op in Phase D; avoid optional `mapping_law`
params."

Three candidates:

- **(a) Structurally absent**: Phase D has zero Phase-E-related
  parameters. Phase E lands as a separate function
  (`step_fa_to_ecm_response_active`) with its own Sanity Gate +
  lock + commit. Mirrors HB#4 `_neutral` / `_active` separation.
- **(b) Explicit Phase D no-op marker**: a sentinel object
  (e.g., `PhaseDNoOpResponse()`) is the only allowed
  `response_law` value in Phase D; Phase E commit changes the
  type signature to accept other concrete laws. Verbose.
- **(c) Optional kwarg with default `None`** (rejected per
  Codex emphasis): `def step_fa_to_ecm_response(adhesions, ecm,
  response_law=None)` — silent activation risk if a future
  caller passes a non-None value before Phase E gates land.

Recommended: (a) structurally absent. Matches HB#4 precedent.

### Q4 — Integrator wrapper sequencing

The phased plan §1 Phase D integrator wrapper:
`cell step → 6.3b FA step → step_fa_to_ecm_response →
step_ecm_to_fa_bias → next step`.

In-place vs return-new ECM:

- **(a) In-place**: the integrator mutates the ECM state
  field-by-field. Codex emphasis #4 challenges this because
  Phase D should not mutate input ECM unless the lock
  explicitly chooses that semantics.
- **(b) Return-new**: the integrator returns a new
  `ECMSubstrateState` (in Phase D, byte-identical to input).
  Phase E response law then mutates the returned state if
  needed. Pure-function semantics preserved.

Recommended: (b). Matches the HB#3 scatter (pure return) and
HB#4 sampler (pure read) sister patterns.

---

## 3. Recommended opening positions (NOT a lock)

| Q | Recommendation | Rationale |
|---|---|---|
| Q1 return type | (b) composite `PhaseDStepResult` dataclass with `traction_density_xy` + `updated_ecm` + `bias_result` | Codex emphasis #2: explicit measurement modality, no hidden state |
| Q2 step_ecm_to_fa_bias | (c) different function, calls `compute_ecm_to_fa_bias_neutral` internally; Phase E adds `step_ecm_to_fa_bias_active` separately | Codex emphasis #3: HB#4 Phase D / Phase E function naming separation precedent |
| Q3 Phase E hook | (a) structurally absent in Phase D | Codex emphasis #3: silent-activation guard; mirrors HB#4 `_neutral` / `_active` |
| Q4 integrator sequencing | (b) return-new ECM, no in-place mutation | Codex emphasis #4: HB#3 / HB#4 sister-pattern (both pure functions) |

**Caveats / unresolved-by-this-brief**:

- The Phase D no-op integrator's `updated_ecm` is trivially
  identical to input. The dataclass field is a future-proofing
  for Phase E response laws that update ECM fields. The round
  may argue (c) hybrid (omit `updated_ecm` until Phase E adds
  it) — both are defensible.
- `bias_params` argument in Q2: the recommendation has the Phase
  D function take only `(adhesions, ecm)`. If a future caller
  wants to pass any Phase D-tuneable knob (e.g., a multiplier
  scale, but Phase D has no such knob by design), the function
  name change is forced. Acceptable for Phase D's "no
  behavior change" scope.

---

## 4. Phase D no-op interface contract (restated, locked elsewhere)

Per `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §2 Phase
D entry contract (NOT changeable in this round):

- `step_fa_to_ecm_response(adhesions, ecm, scattering_params)
  -> {traction_density_xy, ...}`
  - Default no-op: returns zero traction density (or composite
    dataclass per Q1 with zero scatter + identity ecm + neutral
    bias).
- `step_ecm_to_fa_bias(adhesions, ecm, bias_params)
  -> ECMToFABiasResult`
  - Default no-op: returns 1.0 multipliers neutral.

The round defines:
- Concrete `step_*` function signatures matching the recommended
  return types (Q1, Q2).
- `PhaseDStepResult` dataclass shape (Q1).
- Phase E hook treatment (Q3).
- Integrator wrapper return semantics (Q4).
- Empty FA list → zero/neutral result (mirrors HB#3 / HB#4
  empty-list precedent).
- Out-of-grid FA position raise (mirrors HB#3 / HB#4).
- Pure-function constraint (no ECM mutation, no FA mutation, no
  RNG, no state-label transitions, no traction_scale_nN
  scaling, no scalarization at integrator layer).

---

## 5. Sanity Gate items Phase D will need (preview, NOT commitment)

Per Phase D row of locked phased plan §6 Sanity Gate matrix
applied to the integrator:

- §1 Dimensional: integrator inherits HB#3 / HB#4 unit chains;
  no new chain at the integrator layer.
- §2 Boundary: empty FA list returns zero/neutral result; FA
  position validation order matches HB#3 / HB#4 (sister-pattern
  preservation per Codex emphasis #5 + Step 6 sister-gate-mirror).
- §3 Conservation: scatter conservation inherited from HB#3;
  bias multiplier structural invariant inherited from HB#4;
  Phase D integrator returns new ECM = input ECM (identity).
- §4 Numerical: float64 throughout; per-FA work is the sum of
  HB#3 + HB#4 work; no dt-rate gate (single step).
- §5 Sign: scatter sign from HB#3; multiplier sign from HB#4.
- §6 Measurement-protocol: Phase D return shape exposes scatter
  + bias + (identity) updated_ecm without scalarization. A
  runtime meta-test
  `test_phase_d_step_does_not_satisfy_phase_e_response_law`
  enforces the Phase D / Phase E boundary at runtime.

Magic-Number Block: zero new tunables at integrator layer.
Inherits from HB#3 / HB#4.

---

## 6. What this brief is *not*

- Not a lock. The integrator interface is decided by the
  design-discussion adversarial round.
- Not a Sanity Gate document. The Sanity Gate is written by the
  future Phase D code commit, after the lock artifact lands.
- Not a Phase E activation. Phase D no-op default is identity /
  neutral; activating non-identity / non-neutral values is Phase
  E (separate function, separate Sanity Gate, separate lock).
- Not a closed-loop ECM gate satisfaction. Phase D no-op
  integrator does not satisfy any of Items 1-6 of the
  closed-loop ECM gate (those are Phase E's responsibility).
- Not authoritative for biological parameters. Any Phase E
  parameter (e.g., constitutive law coefficient) requires
  literature-first extraction per Magic-Number Block.

---

## 7. Process expectation (mirrors HB#3 / HB#4 precedents)

1. **Cross-room dispatch** to design-discussion with this brief
   + the 4 design questions + the recommended opening positions.
   Status: `decision-needed`.
2. **Adversarial round** in design-discussion (Codex/Claude per
   memory `feedback_aggressive_design_debate.md`). HB#3 took 4
   rounds; HB#4 took 4 rounds; Phase D scaffolding is a smaller
   surface (composition + hooks) so 2-3 rounds plausible.
3. **Unresolved disagreements list** before lock — particularly
   the Q1 return-type composite-vs-literal and Q4 in-place-vs-new.
4. **Lock artifact** delivered to implementation-work as
   `docs/v2_phase_d_no_op_scaffolding_locked.md`, parallel to
   `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md` /
   `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`.
5. **Cross-room dispatch back to implementation-work** with the
   locked integrator contract. Sanity Gate doc + code/tests
   follow per HB#3 / HB#4 sequencing precedent.
6. **Phase D no-op scaffolding committed** = closed-loop ECM
   gate Phase D milestone reached. Phase E unlock is the next
   design-discussion target (constitutive law + saturation +
   Lyapunov metric, locked phased plan §2 Hard Blockers #1 / #2
   / #5).

implementation-work stays idle/review-capable while this design
round runs. Per Codex impl `id=1389` ack, this brief is the
impl-work side's "opening dispatch" matching the HB#3 / HB#4
(B-trigger) pattern.

---

## 8. References

- Locked phased plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- Hard Blocker #3 (FA→ECM scattering, locked + implemented):
  `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md` +
  `docs/v2_fa_to_ecm_scattering_sanity_gate.md` +
  `acs/v2/dynamics/fa_to_ecm_scattering.py` +
  `tests/test_v2_fa_to_ecm_scattering.py`
- Hard Blocker #4 (ECM→FA bias, locked + implemented):
  `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` +
  `docs/v2_ecm_to_fa_bias_sanity_gate.md` +
  `acs/v2/dynamics/ecm_to_fa_bias.py` +
  `tests/test_v2_ecm_to_fa_bias.py`
- Sister-brief precedents (mirror this brief's structure):
  `docs/v2_hard_blocker_3_fa_to_ecm_scattering_brief.md` +
  `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_brief.md`
- 6.3a / 6.3b dynamics (FA step the integrator wraps):
  `acs/v2/dynamics/focal_adhesion.py` +
  `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py`
- ECM substrate schema: `acs/v2/ecm_substrate.py`
- FA schema: `acs/v2/focal_adhesion.py`
- Memory rules:
  - `rule10_unit_derivation_in_docs.md`
  - `hard_rule_11_wording_boundary_meta_test.md`
  - `design_note_pre_commit_batch.md` (5-step pre-commit batch +
    proposed Step 6 sister-gate-mirror baked into this brief's
    pre-commit checklist per Codex `id=1389`)
  - `feedback_aggressive_design_debate.md`

---

## 9. Cross-room dispatch instruction (impl-work → design-discussion)

The MCP message that accompanies this brief:
- `to`: `codex` (design-discussion-pane Codex receives + Claude
  pane reads)
- `room`: `design-discussion`
- `topic`: `v2-layer-2-phase-d-no-op-scaffolding`
- `status`: `decision-needed`
- `body`: brief 4-question summary + opening recommendations +
  reference to this file
- `refs`: implementation-work
  `id=1387, 1388, 1389` + locked phased plan reference + HB#3 /
  HB#4 lock artifacts as precedent

design-discussion round produces a lock artifact that supersedes
this brief; this file remains as historical opening-position
context.
