# V2 Phase 1 Phase D No-Op Scaffolding — Closed-Loop Integrator (LOCKED)

**Date**: 2026-05-04 KST
**Status**: design-discussion lock complete after 2 rounds.
Round 1 = `id=1394` (Codex proposal) + `id=1395` (Claude
counter-challenges), Round 2 = `id=1396` (Codex ACCEPT all 3 +
C2 amendment) + `id=1397` (Claude ACCEPT amendment + lock readiness).
Round 3 ACK pending from Codex after this lock commit.
**Authors**: Claude (impl-work) + Codex (impl-work / design-discussion),
both rooms agreeing.
**Source**: `docs/v2/v2_phase_d_no_op_scaffolding_brief.md` (opening
brief, commit `eafbdb9`) + `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
§1 Phase D + §2 Hard Blocker #3 / #4 entry contracts.
**Supersedes**: the opening brief (which remains as historical
opening-position context).

This lock pins the Phase D no-op closed-loop scaffolding integrator
interface so impl-work can proceed to Sanity Gate + code without
any unresolved design choice. **Phase E activation is BLOCKED** on
all 5 Hard Blockers + the effective_stiffness law decision; this
lock does NOT authorize any active response law, any active bias
law, or any Phase E behavior.

---

## 0. Scope

Phase D no-op scaffolding = the composition layer that wires HB#3
(FA→ECM scatter primitive) and HB#4 (ECM→FA bias primitive) into
a single integrator step entry point. Phase D semantics:

- **FA→ECM leg**: scatter per-FA traction onto the ECM grid as
  vector traction density (`(nx, ny, 2)` `nN/μm²`); ECM state is
  **not updated** — `updated_ecm is ecm` (object identity).
- **ECM→FA leg**: read ECM fields at FA positions and return
  neutral (all-1.0) rate multipliers per HB#4
  `compute_ecm_to_fa_bias_neutral`.
- **Composition wrapper**: call FA→ECM first, then ECM→FA only if
  FA→ECM succeeds. Returns a typed composite dataclass with both
  leg results plus the Phase-D-identity ECM.

The locked phased plan §1 Phase D scope ("schema/API hooks for
FA→ECM and ECM→FA paths, default OFF or identity, no behavior
change") is satisfied by these three pure functions.

**Out of scope** for this lock:

- Active constitutive response law (Hard Blocker #1).
- Active saturation form (Hard Blocker #2).
- Lyapunov-like metric (Hard Blocker #5).
- effective_stiffness helper (locked phased plan §3 guard).
- Any Phase E function or parameter.

---

## 1. Final Scope LOCK

### Module

`acs/v2/dynamics/closed_loop_phase_d.py` (new file). Phase E adds
a separate module (e.g., `closed_loop_phase_e.py`) with its own
lock. **Do not** add Phase E variants to this module.

### Dataclasses

```python
@dataclass(frozen=True, slots=True)
class FAToECMResponseResult:
    traction_density_xy: np.ndarray        # (nx, ny, 2), float64, nN/μm²
    updated_ecm: ECMSubstrateState         # Phase D: same object as input ecm

@dataclass(frozen=True, slots=True)
class PhaseDNoOpStepResult:
    fa_to_ecm: FAToECMResponseResult
    ecm_to_fa: ECMToFABiasResult           # imported from HB#4 module
    updated_ecm: ECMSubstrateState         # Phase D: same object as input ecm
```

`PhaseDNoOpStepResult.updated_ecm is fa_to_ecm.updated_ecm is ecm`
(triple-identity invariant in Phase D — three references to the
same object).

### Functions

```python
def step_fa_to_ecm_response(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
) -> FAToECMResponseResult:
    """Phase D FA→ECM leg: scatter traction + identity ECM response.

    Calls scatter_fa_traction_to_ecm_bilinear(adhesions, ecm) for the
    scatter array. Returns FAToECMResponseResult with the scatter
    array + the input ecm as updated_ecm (object identity).

    Failure: any FA validation failure from HB#3 propagates with the
    HB#3 failure_kind ("fa_position_outside_ecm_grid",
    "non_finite_fa_position", "non_finite_fa_traction") and
    FAToECMScatteringError exception class.
    """


def step_ecm_to_fa_bias(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
) -> ECMToFABiasResult:
    """Phase D ECM→FA leg: neutral bias wrapper.

    Calls compute_ecm_to_fa_bias_neutral(adhesions, ecm) directly.

    Failure: any FA validation failure from HB#4 propagates with the
    HB#4 failure_kind ("fa_bias_position_outside_ecm_grid",
    "non_finite_fa_position") and FAToECMBiasError exception class.
    """


def step_phase_d_no_op(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
) -> PhaseDNoOpStepResult:
    """Phase D no-op closed-loop step: composition wrapper.

    Calls step_fa_to_ecm_response first, then step_ecm_to_fa_bias
    only if FA→ECM succeeds. Returns PhaseDNoOpStepResult with both
    leg results plus updated_ecm (= ecm in Phase D).

    Failure: surfaces the FIRST leg's failure (HB#3 first, then
    HB#4 only if HB#3 succeeded). Wrapper introduces no new
    failure_kinds.
    """
```

All three functions are **pure** — no ECM mutation, no FA
mutation, no RNG, no state-label transitions, no traction scaling,
no scalarization, no Phase E hook.

### Error classes

The wrapper introduces **no new error class**. Failures propagate
from HB#3 (`FAToECMScatteringError`) and HB#4 (`FAToECMBiasError`)
with their locked failure_kinds.

### Forbidden in Phase D

- Loose `diagnostics_dict: dict[str, ...]` on `FAToECMResponseResult`
  (typed-contract regression — rejected in round 1).
- Optional Phase E params (`response_law`, `mapping_law`,
  `bias_params`) — silent activation risk (Codex emphasis #3).
- Active response law of any kind.
- Active bias law of any kind.
- In-place ECM mutation.
- Deep-copy ECM that pretends to be an update.
- New failure_kinds at the wrapper layer (delegate fully to HB#3 /
  HB#4 primitives).
- Scalarization of `traction_density_xy` at any layer (Hard Rule
  11 — Phase E response law owns the reduction).
- effective_stiffness helper or any geometry-field consuming
  shortcut (locked phased plan §3 side-door guard).

---

## 2. Coordinate Convention

Inherits from HB#3 + HB#4 locks (cell-centered grid, inclusive
physical footprint, boundary half-cell truncated stencil with
renormalized weights). The wrapper introduces no new coordinate
convention.

---

## 3. Call Ordering (deterministic)

Inside `step_phase_d_no_op`:

1. `fa_to_ecm = step_fa_to_ecm_response(adhesions, ecm)`
   - If this raises any HB#3 error, propagate immediately. Do NOT
     call `step_ecm_to_fa_bias`.
2. `ecm_to_fa = step_ecm_to_fa_bias(adhesions, ecm)`
   - If this raises any HB#4 error, propagate (but at this point
     HB#3 has already validated the same adhesions, so HB#4 is
     called on validated input).
3. Return `PhaseDNoOpStepResult(fa_to_ecm=fa_to_ecm,
   ecm_to_fa=ecm_to_fa, updated_ecm=ecm)`.

**Why this order**: HB#3 and HB#4 use distinct failure_kinds for
the same physical condition (`fa_position_outside_ecm_grid` vs
`fa_bias_position_outside_ecm_grid`). Without deterministic order,
the wrapper's external behavior is non-reproducible. FA→ECM first
because the scatter is the primary action of the closed-loop
step; ECM→FA bias is the rate-modulation read for the next step.

---

## 4. Empty-FA-List Invariants (LOCKED)

For `adhesions = ()` (or `[]`) input:

- `fa_to_ecm.traction_density_xy.shape == (nx, ny, 2)`, dtype
  `np.float64`.
- `np.all(fa_to_ecm.traction_density_xy == 0.0)` (exact zero, not
  tolerance-bound).
- `fa_to_ecm.updated_ecm is ecm` (object identity).
- `ecm_to_fa.multipliers_per_fa.shape == (0, 3)`, dtype
  `np.float64`.
- `ecm_to_fa.rate_names == RATE_NAMES` per HB#4 lock
  (`("k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s")`).
- `ecm_to_fa.sampled_diagnostics` is a non-`None` empty
  `ECMSampledAtFAs` (preserving HB#4 sampler semantics — never
  `None`).
- `result.updated_ecm is ecm` (object identity).

These are runtime-meta-test targets; see §6 test catalog.

---

## 5. Sanity Gate 6 Items (preview — Sanity Gate doc owns)

The Sanity Gate doc `docs/v2/v2_phase_d_no_op_scaffolding_sanity_gate.md`
(future, written by impl-work) must address:

- **§1 Dimensional**: integrator inherits HB#3 / HB#4 unit chains;
  no new chain at the wrapper layer.
- **§2 Boundary**: empty FA list returns the locked invariants in
  §4; FA position validation order matches HB#3 / HB#4 (sister-pattern
  preservation per Codex emphasis #5 + Step 6 sister-gate-mirror).
- **§3 Conservation**: scatter conservation inherited from HB#3;
  bias multiplier structural invariant inherited from HB#4; wrapper
  preserves both component-wise.
- **§4 Numerical**: float64 throughout; per-FA work = HB#3 work +
  HB#4 work; no dt-rate gate (single step).
- **§5 Sign**: scatter sign from HB#3; multiplier sign from HB#4;
  wrapper preserves both.
- **§6 Measurement-protocol**: wrapper exposes scatter array + bias
  multipliers + (identity) updated_ecm without scalarization.
  Runtime meta-test
  `test_phase_d_step_does_not_satisfy_phase_e_response_law`
  enforces the Phase D / Phase E boundary.

**Magic-Number Block**: zero new tunables at wrapper layer. Inherits
fully from HB#3 / HB#4.

---

## 6. Test Catalog (preview ~12 tests)

Owned by the future code commit. Each test maps to a locked
invariant.

1. `test_phase_d_step_returns_typed_dataclass` — return type
   `PhaseDNoOpStepResult` with `frozen=True, slots=True`.
2. `test_phase_d_step_object_identity_updated_ecm` — `result.updated_ecm
   is ecm and result.fa_to_ecm.updated_ecm is ecm`.
3. `test_phase_d_step_no_input_mutation` — input `ecm` and FA states
   byte-identical before and after the call.
4. `test_phase_d_step_empty_fa_list` — empty-FA-list invariants per
   §4 locked list.
5. `test_phase_d_step_call_order_deterministic` — out-of-grid FA
   surfaces as HB#3 `fa_position_outside_ecm_grid` (NOT HB#4
   `fa_bias_position_outside_ecm_grid`); when HB#3 raises, HB#4 is
   not called (mockable).
6. `test_phase_d_step_no_phase_e_kwargs` — `step_phase_d_no_op`,
   `step_fa_to_ecm_response`, `step_ecm_to_fa_bias` all have
   exactly two positional parameters (`adhesions`, `ecm`); zero
   optional kwargs (introspect via `inspect.signature`).
7. `test_phase_d_step_legs_match_primitives` — `step_fa_to_ecm_response`
   output equals `scatter_fa_traction_to_ecm_bilinear` output;
   `step_ecm_to_fa_bias` output equals `compute_ecm_to_fa_bias_neutral`
   output (no mutation, identical traction array, identical
   multipliers).
8. `test_phase_d_step_propagates_hb3_failure_kinds` — non-finite FA
   position raises `FAToECMScatteringError` with
   `failure_kind="non_finite_fa_position"`.
9. `test_phase_d_step_propagates_hb4_failure_kinds_when_bias_leg_raises`
   — after FA→ECM succeeds, a mocked/stubbed HB#4 failure (for
   example monkeypatching `step_ecm_to_fa_bias` to raise
   `FAToECMBiasError`) is propagated unchanged. Natural
   geometry-boundary inputs should have matching HB#3/HB#4
   pass/fail behavior; a geometry case that passes HB#3 but fails
   HB#4 is a sister-gate-mirror bug, not a valid propagation
   scenario.
10. `test_phase_d_step_scatter_conservation_inherited` — sum over
    grid * cell area equals sum over FAs (HB#3 invariant
    propagated through wrapper).
11. `test_phase_d_step_does_not_satisfy_phase_e_response_law`
    (Hard Rule 11 runtime meta-test) — wrapper does NOT update ECM
    fields based on traction density (assert ECM fields unchanged
    even with non-trivial scatter input).
12. `test_phase_d_step_no_diagnostics_dict_field` — assert
    `FAToECMResponseResult.__dataclass_fields__.keys() ==
    {"traction_density_xy", "updated_ecm"}` exactly.

---

## 7. Implementation Files

- `acs/v2/dynamics/closed_loop_phase_d.py` — module home for the
  three Phase D functions + the two new dataclasses.
- `tests/v2/test_v2_closed_loop_phase_d.py` — test file for the ~12
  tests above.
- `docs/v2/v2_phase_d_no_op_scaffolding_sanity_gate.md` — Sanity Gate
  doc (future, impl-work writes after Codex ACK on this lock).

Imports the wrapper module needs:

- HB#3: `from acs.v2.dynamics.fa_to_ecm_scattering import
  scatter_fa_traction_to_ecm_bilinear, FAToECMScatteringError`
- HB#4: `from acs.v2.dynamics.ecm_to_fa_bias import
  compute_ecm_to_fa_bias_neutral, ECMToFABiasResult,
  FAToECMBiasError`
- Schemas: `from acs.v2.ecm_substrate import ECMSubstrateState`,
  `from acs.v2.focal_adhesion import FocalAdhesionState`

The wrapper module does **not** import any Phase E module
(structural absence guard).

---

## 8. Sequencing & Cross-Room Dispatch

### Phase D code unblock — after this lock + Codex ACK

1. This lock doc commits to `docs/v2/v2_phase_d_no_op_scaffolding_locked.md`.
2. Codex ACK in design-discussion (one round 4 message confirming
   the lock).
3. Cross-room dispatch back to implementation-work: "Phase D no-op
   scaffolding lock complete; impl-work begins Sanity Gate doc."
4. impl-work Sanity Gate doc draft:
   `docs/v2/v2_phase_d_no_op_scaffolding_sanity_gate.md`.
5. Codex review of Sanity Gate doc.
6. impl-work code commit:
   `acs/v2/dynamics/closed_loop_phase_d.py` +
   `tests/v2/test_v2_closed_loop_phase_d.py`.
7. Codex review of code commit.
8. Phase D no-op scaffolding committed = closed-loop ECM gate
   Phase D milestone reached.

### Phase E unblock — much later

Phase E remains BLOCKED on all 5 Hard Blockers (#1, #2, #3 ✓, #4 ✓,
#5) plus the effective_stiffness law decision. Phase E lock will be
a separate design-discussion round, separate lock artifact, separate
module (`closed_loop_phase_e.py`).

### Cadence rule

implementation-work cadence: ~10 min check during active design
rounds, ~25 min during impl review wait, ~30 min during PI watchdog
windows. Per memory `cadence_promise_must_send_even_when_idle`,
heartbeat at promised cadence even when no movement.

---

## 9. References

- Phased plan: `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
  §1 Phase D + §2 Hard Blocker #3/#4 entry contracts
- Hard Blocker #3 lock + impl + tests:
  - `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2/v2_fa_to_ecm_scattering_sanity_gate.md`
  - `acs/v2/dynamics/fa_to_ecm_scattering.py`
  - `tests/v2/test_v2_fa_to_ecm_scattering.py`
- Hard Blocker #4 lock + impl + tests:
  - `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2/v2_ecm_to_fa_bias_sanity_gate.md`
  - `acs/v2/dynamics/ecm_to_fa_bias.py`
  - `tests/v2/test_v2_ecm_to_fa_bias.py`
- Opening brief (superseded by this lock, kept for history):
  `docs/v2/v2_phase_d_no_op_scaffolding_brief.md`
- Design-discussion message thread:
  `id=1391` (impl-work dispatch) →
  `id=1394` (Codex round 1 proposal) →
  `id=1395` (Claude round 1 counter-challenges) →
  `id=1396` (Codex round 2 ACCEPT all 3 + C2 amendment) →
  `id=1397` (Claude round 2 ACCEPT amendment + lock readiness)
- Memory rules informing this lock:
  - `rule10_unit_derivation_in_docs.md`
  - `hard_rule_11_wording_boundary_meta_test.md`
  - `design_note_pre_commit_batch.md` (5-step + Step 6 sister-gate-mirror
    baked into Phase D brief / lock pre-commit)
  - `feedback_aggressive_design_debate.md`
  - `cadence_promise_must_send_even_when_idle.md`
  - `no_standby_after_completion.md`

---

## 10. Confirmation log

| Round | Sender | Message id | Action |
|---|---|---|---|
| 1 (Codex open) | Codex (design-discussion) | `id=1394` | Proposed 3-function split + object-identity + deterministic order + `diagnostics_dict` |
| 1 (Claude reply) | Claude (design-discussion) | `id=1395` | ACCEPT 3-function split + object-identity + ordering; counter-challenge `diagnostics_dict` (drop), empty-FA invariants, file location |
| 2 (Codex reply) | Codex (design-discussion) | `id=1396` | ACCEPT all 3 counter-challenges; amend C2 sampled_diagnostics empty (not None) |
| 2 (Claude close) | Claude (design-discussion) | `id=1397` | ACCEPT amendment; lock-ready statement |
| Lock commit | Claude (impl-work) | (this commit) | Lock artifact `docs/v2/v2_phase_d_no_op_scaffolding_locked.md` |
| 3 (Codex ACK pending) | Codex (design-discussion) | TBD | Confirm lock holds |

This lock is binding once Codex's round-3 ACK lands. impl-work
begins Sanity Gate doc only after that confirmation.
