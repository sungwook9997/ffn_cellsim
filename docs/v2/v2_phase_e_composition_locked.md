# V2 Closed-Loop ECM Gate — Phase E Composition v1 (ECM-Side Evidence Only) — Locked Design

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (4-round adversarial lock,
PI id=809 aggressive debate posture applied, PI id=1008/1057/1161
autonomy + visible deliverable focus)
**Source unit**: design-discussion `topic=v2-phase-e-composition`,
MCP id 1526–1532 (rounds 1–4)
**Brief**: `docs/v2/v2_phase_e_composition_brief.md` (commit `7529d3e`)
**Parent locked plan**: `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
(Phase E entry — final upstream composition step)
**Upstream sister locks (all 5 hard blockers + Phase D)**:
- `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
- `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
- `docs/v2/v2_phase_d_no_op_scaffolding_locked.md`
- `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`
- `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`
- `docs/v2/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1)
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for the Phase E v1 Sanity Gate doc + code entry.

---

## 0. Scope

This unit locks the **Phase E v1 composition** — the first integrated
closed-loop ECM gate composition combining the 5 upstream hard blockers
+ Phase D infrastructure into a single per-step function.

### v1 vs v2 boundary (the central design decision)

**v1 (this lock's scope)**:
- HB#4 hard-wired to **neutral** (`compute_ecm_to_fa_bias_neutral`); FA
  rates are NOT modulated by ECM state in v1
- ECM **does** evolve under FA traction (HB#1+#2 active orientation
  update)
- Provides **ECM-side composition evidence** for parent plan's closed-
  loop gate Items 1–4
- **Does NOT close or satisfy the full FA→ECM→FA loop** — full
  closed-loop satisfaction requires HB#4-active design (separate v2
  cycle) AND Item 5 sweep harness AND an explicit PI-approved decision
  if ECM-side-only is enough
- effective_stiffness law: **NOT used** (no helper imported, no symbol
  referenced; deferred to v2 ECM→FA mechanosensing — Y8/Y13)

**v2 (separate later cycle, NOT this unit)**:
- HB#4-active variant (`compute_ecm_to_fa_bias_active` — does not yet
  exist) replaces neutral readout
- Full FA→ECM→FA closure
- Possibly consumes effective_stiffness law (separate decision)
- Will get its own lock + function `step_closed_loop_phase_e_v2`

### Wording discipline (Y1 — central anchor)

This lock follows the parent plan's "precursor not satisfaction"
discipline. Per Codex `id=1527` C1:

- ✅ "Phase E v1 **demonstrates ECM-side composition evidence** for
  Items 1–4 under prescribed/scattered FA traction and HB#4-neutral
  readout"
- ✅ "Phase E v1 **does not close or satisfy** the full closed-loop ECM
  gate"
- ✅ "Full gate satisfaction remains pending Item 5 sweep harness AND
  Phase E v2/HB#4-active OR an explicit PI-approved decision that
  ECM-side-only is enough"

The brief's `test_phase_e_satisfies_item_X` naming is **superseded** by
`test_phase_e_v1_provides_item_X_*_evidence` (Y2) — test names are the
contract enforcement layer; satisfaction-naming would silently overclaim.

---

## 1. Final Lock Summary

### Function

```python
"""V2 Phase E closed-loop ECM gate composition (v1: ECM-side evidence only).

Composes HB#3 scatter → HB#1+#2 orientation update → HB#4 neutral bias
readout → HB#5 Lyapunov-like metric. v1 hard-wires HB#4 neutral; v2
with active HB#4 is a separate cycle.
"""

import numpy as np
from dataclasses import dataclass

from acs.v2.focal_adhesion import FocalAdhesionState                 # Y15: corrected import path
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.dynamics.fa_to_ecm_scattering import scatter_fa_traction_to_ecm_bilinear
from acs.v2.dynamics.ecm_constitutive_response import (
    K_ORIENT_PER_S,                    # imported only, NOT redefined (Y14)
    TRACTION_REF_NN_PER_UM2,           # imported only, NOT redefined (Y14)
    ECMOrientationResponseResult,
    step_ecm_orientation_response,
)
from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMToFABiasResult,
    compute_ecm_to_fa_bias_neutral,
)
from acs.v2.dynamics.ecm_lyapunov_metric import (
    ECMOrientationLyapunovMetricResult,
    compute_ecm_orientation_lyapunov_metric,
)


@dataclass(frozen=True, slots=True)
class PhaseEStepResult:
    """Outputs of one Phase E v1 closed-loop ECM gate composition step.

    Identity invariant (Y4): ``updated_ecm is orientation_response.updated_ecm``
    (no copy, no mutation in this composition wrapper).
    """
    traction_density_xy: np.ndarray                       # HB#3 output (nx, ny, 2) nN/μm²
    orientation_response: ECMOrientationResponseResult    # HB#1+#2 output
    ecm_to_fa_bias: ECMToFABiasResult                     # HB#4 output (neutral in v1)
    lyapunov_metric: ECMOrientationLyapunovMetricResult   # HB#5 output (ECM-side evidence)
    updated_ecm: ECMSubstrateState                        # = orientation_response.updated_ecm


def step_closed_loop_phase_e_v1(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
    dt_s: float,
    *,
    traction_ref_nN_per_um2: float = TRACTION_REF_NN_PER_UM2,
    k_orient_per_s: float = K_ORIENT_PER_S,
) -> PhaseEStepResult:
    """Phase E v1 ECM-side composition step.

    Composes:
      1. HB#3 scatter (FA traction → ECM grid; empty FA → all-zero per
         source line 170)
      2. HB#1+#2 orientation update (instantaneous traction stimulus)
      3. HB#4 neutral bias readout (post-update ECM)
      4. HB#5 Lyapunov-like metric (post-update ECM, original scatter)

    No wrapper failure kinds; sub-call errors propagate as-is.
    Deterministic stop order: HB#3 → HB#1+#2 → HB#4 → HB#5.
    """
    traction_density_xy = scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)
    orientation_response = step_ecm_orientation_response(
        ecm, traction_density_xy, dt_s,
        traction_ref_nN_per_um2=traction_ref_nN_per_um2,
        k_orient_per_s=k_orient_per_s,
    )
    ecm_to_fa_bias = compute_ecm_to_fa_bias_neutral(adhesions, orientation_response.updated_ecm)
    lyapunov_metric = compute_ecm_orientation_lyapunov_metric(
        orientation_response.updated_ecm, traction_density_xy
    )
    return PhaseEStepResult(
        traction_density_xy=traction_density_xy,
        orientation_response=orientation_response,
        ecm_to_fa_bias=ecm_to_fa_bias,
        lyapunov_metric=lyapunov_metric,
        updated_ecm=orientation_response.updated_ecm,  # object identity (Y4)
    )
```

### Composition order (Y5 — locked)

1. **HB#3 scatter**: `scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)`
   — FA traction → ECM grid traction density. Empty FA list returns
   all-zero per source contract `acs/v2/dynamics/fa_to_ecm_scattering.py:170`.
2. **HB#1+#2 orientation update**:
   `step_ecm_orientation_response(ecm, traction_density_xy, dt_s, ...)`
   — instantaneous traction stimulus, exact-exponential convex update,
   returns fresh ECM with no aliasing.
3. **HB#4 neutral bias readout**:
   `compute_ecm_to_fa_bias_neutral(adhesions, orientation_response.updated_ecm)`
   — uses **post-update** ECM. Future-proof for v2 (in v1, HB#4 neutral
   is state-independent so order is moot, but locked for forward-compat).
4. **HB#5 Lyapunov-like metric**:
   `compute_ecm_orientation_lyapunov_metric(orientation_response.updated_ecm, traction_density_xy)`
   — uses **post-update** ECM (so V is measured AT the new state) and
   the **original scatter output** (NOT recomputed; scatter is
   deterministic, but the call-order test enforces this).

### Error propagation (Y6 — locked)

- Composition wrapper introduces **NO new failure kinds** (Y17 source
  guard enforces no `PhaseE*Error` class, no `failure_kind =`
  assignment in module).
- Sub-call errors propagate **as-is** (HB#3 / HB#1+#2 / HB#4 / HB#5
  errors retain their original types and messages).
- Deterministic stop order:
  1. HB#3 failure → stop before HB#1+#2 / HB#4 / HB#5
  2. HB#1+#2 failure → stop before HB#4 / HB#5
  3. HB#4 failure → stop before HB#5
  4. HB#5 failure → only after successful upstream

### Identity invariant (Y4 — locked)

- `result.updated_ecm is result.orientation_response.updated_ecm` —
  **object identity**, not just bytewise equality
- Composition wrapper does NOT copy or mutate sub-results
- HB#1+#2 no-aliasing (per its lock Y12) guarantees `updated_ecm` is
  fresh from input; Phase E does NOT add another copy layer

### Empty-FA chain (Y7+Y11 — locked behavior)

For `adhesions = ()`:
- HB#3 returns all-zero traction (source `fa_to_ecm_scattering.py:170`)
- HB#1+#2 with zero traction → `w = 0` everywhere → identity update;
  but **`result.updated_ecm` is still a fresh ECMSubstrateState** with
  all 5 array fields equal-but-not-identical to input (no shared
  memory, per HB#1+#2 Y12)
- HB#4 neutral with empty adhesions → `multipliers_per_fa.shape == (0, 3)`
- HB#5 with zero traction → `active.any() == False` → `V_active = 0.0`,
  passive_orientation_magnitude reflects input orientation magnitude
- `result.updated_ecm is ecm` is **False** (fresh wrapper from HB#1+#2)
- `result.updated_ecm.orientation_tensor == ecm.orientation_tensor`
  bytewise (no in-place edit)

### Exports + constant guard (Y14+Y16 — locked)

- `acs/v2/dynamics/__init__.py` adds 2 new exports + alphabetical
  reorder of `__all__`:
  - `PhaseEStepResult`
  - `step_closed_loop_phase_e_v1`
- `acs/v2/__init__.py` mirrors the 2 additions
- `closed_loop_phase_e.py` does NOT define new uppercase numeric
  constants beyond imported HB constants
- `TRACTION_REF_NN_PER_UM2` and `K_ORIENT_PER_S` imported via
  `from acs.v2.dynamics.ecm_constitutive_response import ...`; tests
  enforce `ce.TRACTION_REF_NN_PER_UM2 is cr.TRACTION_REF_NN_PER_UM2`
  (object identity = re-import, not redefinition)
- HB#4 result attribute: `multipliers_per_fa` (NDArray shape `(N_FA, 3)`),
  NOT `.multipliers` (Y16 — verified at
  `acs/v2/dynamics/ecm_to_fa_bias.py:138`)

### Forbidden in Phase E v1

- "Items 1-4 satisfied" / "full closed-loop" / "full FA→ECM→FA" /
  "gate satisfied" wording in module docstring, function docstring, or
  return-class docstrings (Y1 + Y2 forward guards via meta-test 14)
- Test names containing `satisfies_item_X` (Y2 — use
  `provides_item_X_*_evidence` instead)
- v2 placeholder, enum, hook, strategy, or optional bias law parameter
  (Y3 — silent-activation risk; v2 gets a separate function lock when
  HB#4-active exists)
- In-place ECM mutation, adhesion mutation, or copying of sub-results
  (Y4 — composition wrapper is a thin orchestrator)
- Wrapper failure kinds: any `PhaseE*Error` class definition or
  `failure_kind = ...` assignment in `closed_loop_phase_e.py` (Y17 —
  sub-call errors propagate verbatim)
- Reading or referencing `effective_stiffness` anywhere in
  `closed_loop_phase_e.py` (Y8 + Y13 — string + AST guard via meta-test
  15; the lock doc itself can explain, but the .py file cannot)
- Local redefinition of `TRACTION_REF_NN_PER_UM2` or `K_ORIENT_PER_S`
  (Y14 — must be imported from `ecm_constitutive_response`; meta-test
  16 enforces object identity)
- Object identity claim `result.updated_ecm is ecm` (Y11 — false;
  HB#1+#2 always returns a fresh wrapper)
- Use of `result.ecm_to_fa_bias.multipliers` (Y16 — wrong field name;
  use `multipliers_per_fa`)
- Reordering composition (Y5 — locked HB#3 → HB#1+#2 → HB#4 → HB#5;
  call-order spy test 4 enforces)
- Recomputing scatter inside HB#5 call (Y5 — must pass the original
  `traction_density_xy` from step 1, not a recomputed one)

---

## 2. Reasoned-acceptance trace (Y1–Y17)

The lock converged after 4 rounds. Each Y is a Claude concession with
the round in which it was accepted, preserving the adversarial-debate
audit trail (PI id=809 posture).

**Y1 (round 2, accepted Codex C1 = ECM-side evidence wording, NOT "Items 1-4 satisfied")**:
- Claude opening lean: tests named `test_phase_e_satisfies_item_X`,
  framing v1 as full gate satisfaction.
- Codex catch: with HB#4 neutral and no FA dynamics update, v1 cannot
  truthfully satisfy "bounded feedback in single-cell loop" (parent
  plan Item 4). It can provide ECM-side composition evidence only.
  Same discipline as parent plan's "precursor not satisfaction"
  language for earlier phases.
- Resolution: lock §0 + §1 use "ECM-side composition evidence" wording
  consistently; tests renamed to `provides_*_evidence` (Y2);
  meta-test 14 forward-guards against overclaim wording in source.

**Y2 (round 2, accepted Codex C2 = test rename `satisfies_*` →
`provides_*_evidence`)**:
- Test names are the contract enforcement layer.
- Resolution: 4 evidence tests named with `provides_item_X_*_evidence`;
  meta-test 14 asserts module/function docstrings do NOT contain
  satisfaction-claim phrases.

**Y3 (round 2, accepted Codex C3 = (c) v1-only, no v2 placeholder)**:
- Codex: do NOT implement v2 placeholder, enum, hook, strategy, or
  optional bias law. Future v2 gets a separate lock and function.
- Resolution: function name `step_closed_loop_phase_e_v1`; future
  `_v2` is a separate lock cycle; HB#1+#2 sister-pattern (Y9 there).

**Y4 (round 2, accepted Codex C4 = identity invariant + no mutation)**:
- Codex required: `result.updated_ecm is
  result.orientation_response.updated_ecm` object identity invariant;
  composition wrapper does NOT copy or mutate sub-results.
- Resolution: dataclass docstring + test 2 enforce.

**Y5 (round 2, accepted Codex C5 = composition order + call-order spy
test)**:
- Order: HB#3 → HB#1+#2 → HB#4(post-update) → HB#5(post-update,
  original scatter).
- Call-order test via monkeypatch with spies — verifies HB#5 receives
  `orientation_response.updated_ecm` and the **original** scatter
  output (not recomputed).
- Resolution: §1 composition order block + test 4.

**Y6 (round 2, accepted Codex C6 = error-propagation contract +
deterministic stop order)**:
- No wrapper failure kinds; sub-call errors propagate as-is.
- Stop-before-next-leg deterministic order.
- Resolution: §1 error propagation block + tests 6–9.

**Y7 (round 2, accepted Codex C7 = empty FA list contract verified at HB#3 source)**:
- Verified at `acs/v2/dynamics/fa_to_ecm_scattering.py:170`:
  ```python
  out = np.zeros((nx, ny, 2), dtype=np.float64)
  if not adhesions:
      return out
  ```
- Resolution: empty-FA chain test 5 cites the source line; Phase E
  does NOT invent empty-handling.

**Y8 (round 2, accepted Codex C8 = effective_stiffness "not used" +
meta-test enforcement)**:
- Codex: state explicit deferral + add source/meta test blocking
  `effective_stiffness` string/name in `closed_loop_phase_e.py`.
- Resolution: §0 deferral wording + meta-test 15 (string + AST check);
  Y13 strict tightening.

**Y9 (round 3, accepted Codex A1 = multi-step boundedness wording —
ECM-side ONLY, no monotone-decrease claim)**:
- Multi-step rotating-traction test: assert `V_active ≤ V_bound` at
  every step + neutral multipliers preserved across all steps.
- Do NOT assert monotone decrease (would be wrong claim under
  varying-target regime).
- Resolution: test 13 framed as "ECM-side boundedness evidence" only;
  no convergence wording.

**Y10 (round 3, accepted Codex A2 = call-order spy test scaffold
correction)**:
- Claude round 2 pseudocode had two bugs:
  1. `calls[0][2]['return']` won't work (fakes don't write returns
     into kwargs)
  2. `calls[2][1][0] is mock_response.updated_ecm` was off-by-one (HB#4
     signature is `(adhesions, ecm)`, so `args[1]` is the ECM)
- Resolution: corrected scaffold pattern with explicit `scatter_out` /
  `updated_ecm` test-fixture variables and explicit signature-aware
  index assertions.

**Y11 (round 3, accepted Codex A3 = empty-FA wording fresh-but-equal,
NOT object identity)**:
- HB#1+#2 always returns fresh ECMSubstrateState (per its Y12
  no-aliasing); `result.updated_ecm is ecm` is **False** even on
  identity step.
- Resolution: §1 empty-FA section uses "fresh ECMSubstrateState with
  all 5 array fields equal-but-not-identical to input"; test 5
  asserts `not np.shares_memory(...)` and `result.updated_ecm is not ecm`.

**Y12 (round 3, accepted Codex A4 = monkeypatch IMPORTED symbol in
`closed_loop_phase_e`, not source module)**:
- Standard mock pattern: monkeypatch where the symbol is **used**, not
  where it's defined.
- Resolution: tests 4 + 6–9 use `import acs.v2.dynamics.closed_loop_phase_e as ce` +
  `monkeypatch.setattr(ce, "scatter_fa_traction_to_ecm_bilinear", ...)` etc.

**Y13 (round 3, accepted Codex A5 = `effective_stiffness` term
completely excluded from module body/docstring)**:
- Strict: `closed_loop_phase_e.py` source must NOT contain the string
  `effective_stiffness` anywhere — function body, docstring, comment,
  module docstring.
- The lock doc explains; the .py file does not need to.
- Resolution: meta-test 15 does both AST walk (for `Attribute.attr` /
  `Name.id`) AND raw `"effective_stiffness" not in src` string check.

**Y14 (round 3, accepted Codex A6 = lean exports + no constant
redefinition guard)**:
- Exports: only `PhaseEStepResult` + `step_closed_loop_phase_e_v1`
  through both `acs.v2.dynamics` and `acs.v2`.
- No new local constants; HB#1+#2 constants imported only.
- Resolution: meta-test 16 asserts object identity
  (`ce.TRACTION_REF_NN_PER_UM2 is cr.TRACTION_REF_NN_PER_UM2`) and
  scans for surprise uppercase numeric constants.

**Y15 (round 4, accepted Codex B1 = FocalAdhesionState import path
correction)**:
- Claude round 3 pseudocode: `from acs.v2.adhesion_state import FocalAdhesionState`
- Verified: `acs/v2/adhesion_state.py` does NOT exist; class lives at
  `acs/v2/focal_adhesion.py:51`.
- Resolution: corrected to `from acs.v2.focal_adhesion import FocalAdhesionState`.

**Y16 (round 4, accepted Codex B2 = HB#4 result field is
`multipliers_per_fa`, NOT `multipliers`)**:
- Claude round 3 test pseudocode used `result.ecm_to_fa_bias.multipliers`
  which does not exist.
- Verified: `ECMToFABiasResult.multipliers_per_fa` at
  `acs/v2/dynamics/ecm_to_fa_bias.py:138` (NDArray shape `(N_FA, 3)`).
- Resolution: corrected all tests to use
  `result.ecm_to_fa_bias.multipliers_per_fa`; empty-FA test asserts
  `.shape == (0, 3)`.

**Y17 (round 4, accepted Codex B3 = wrapper-failure-kind test as
source/meta guard, NOT runtime test)**:
- Claude round 3 had `test_phase_e_v1_no_wrapper_failure_kinds` as a
  runtime test, redundant with tests 6–9.
- Codex: simpler form is source/meta guard — assert no
  `PhaseE*Error` class definition + no `failure_kind =` assignment in
  module source.
- Resolution: test 3 reworded to
  `test_phase_e_v1_module_does_not_define_wrapper_failure_kinds`;
  AST walk + string check; runtime behavior covered by tests 6–9.

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

1. **Units** (composition-level invariants):
   - Inputs: `adhesions: tuple[FocalAdhesionState, ...]` (validated by
     HB#3); `ecm: ECMSubstrateState` (validated by HB#1+#2 entry);
     `dt_s: float [s]` (validated by HB#1+#2)
   - Outputs: each sub-result carries its own units per its own lock
     (HB#3 nN/μm² traction; HB#1+#2 ECM with orientation
     dimensionless; HB#4 dimensionless multipliers; HB#5 V_active in
     μm²)
   - Composition wrapper performs NO unit transformation
2. **Boundary**:
   - Empty FA list: chain produces zero-traction → identity orientation
     → empty bias `(0, 3)` → V_active=0 (Y7+Y11)
   - Single-cell ECM, single FA: each sub-call handles per its own
     boundary
   - dt_s == 0: HB#1+#2 valid no-op → identity update → V_active = 0.5
     · cell_area · Σ ‖T_old − n⊗n‖²_F (current state); HB#4 neutral
     unchanged
   - dt_s < 0 / nonfinite / boolean: HB#1+#2 raises (propagates)
3. **Conservation**: composition wrapper conserves nothing of its own
   (it has no state); each sub-call conserves per its lock. The
   identity invariant `result.updated_ecm is
   result.orientation_response.updated_ecm` is the composition's
   structural invariant.
4. **Numerical**: composition wrapper introduces no numerical
   approximation. Each sub-call uses its own numerical regime
   (HB#1+#2 expm1; HB#5 einsum). Float widths consistent across the
   chain (np.float64 enforced by upstream).
5. **Sign**: each sub-call enforces its own sign convention. The
   wrapper has no sign decisions of its own.
6. **Measurement-protocol** (Hard Rule 11, the central anchor for this
   composition):
   - HB#5 measurement domain (active cells only) matches HB#1+#2
     update domain (active cells only) — sister-pattern guard already
     enforced at HB#5 lock Y1
   - Phase E v1 wording (this lock §0 + §1 + §2) explicitly limits
     evidence to **ECM-side**; full FA→ECM→FA loop satisfaction is
     out of scope
   - Test names (`provides_*_evidence`, NOT `satisfies_*`) are the
     contract enforcement layer (Y2); meta-test 14 enforces that the
     module source does not silently overclaim via docstring wording
   - Items 1–4 evidence framing maps to:
     - Item 1 (response monotonicity) → ECM-side: monotone in
       traction strength (single-cell test)
     - Item 2 (saturation) → ECM-side: bounded by
       `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 4.5` per
       active cell
     - Item 3 (zero-stimulus zero-response) → ECM-side: zero traction
       → identity orientation update + V_active = 0
     - Item 4 (bounded feedback in loop) → ECM-side:
       multi-step rotating-traction stays in `[0,
       v_active_max_bound_um2]` and HB#4 multipliers stay neutral

---

## 4. Test Catalog (17 tests)

### Composition invariants (5)

1. `test_phase_e_v1_returns_phase_e_step_result_dataclass` — return is
   `PhaseEStepResult` instance with all 5 expected fields populated
2. `test_phase_e_v1_identity_invariant_updated_ecm_is_orientation_response_updated_ecm`
   (Y4) — `result.updated_ecm is result.orientation_response.updated_ecm`
3. `test_phase_e_v1_module_does_not_define_wrapper_failure_kinds`
   (Y17) — source/meta guard: no `PhaseE*Error` class, no
   `failure_kind =` assignment in module body (AST walk + string
   check)
4. `test_phase_e_v1_call_order_via_monkeypatch` (Y10 corrected
   pattern) — verifies HB#3 → HB#1+#2 → HB#4 → HB#5 order; HB#4
   receives `(adhesions, post-update ECM)`; HB#5 receives `(post-update
   ECM, original scatter)`
5. `test_phase_e_v1_empty_adhesions_chain` (Y11) — empty FA list →
   all-zero scatter (cite `fa_to_ecm_scattering.py:170`); orientation
   bytewise unchanged but fresh array (no aliasing); HB#4 result shape
   `(0, 3)`; V_active = 0; identity invariant holds

### Error propagation (4)

6. `test_phase_e_v1_propagates_hb3_error_stops_before_orient` —
   monkeypatch HB#3 to raise; verify HB#1+#2/HB#4/HB#5 not called
7. `test_phase_e_v1_propagates_hb12_error_stops_before_bias` —
   monkeypatch HB#1+#2 to raise; verify HB#4/HB#5 not called
8. `test_phase_e_v1_propagates_hb4_error_stops_before_lyap` —
   monkeypatch HB#4 to raise; verify HB#5 not called
9. `test_phase_e_v1_propagates_hb5_error_at_end` — monkeypatch HB#5 to
   raise; verify HB#3/HB#1+#2/HB#4 all called first

### ECM-side evidence (Items 1-4 EVIDENCE, NOT satisfaction) (4)

10. `test_phase_e_v1_provides_item_1_response_monotonicity_evidence` —
    apply with two traction magnitudes; assert larger traction →
    larger `convex_weight` and faster orientation alignment
11. `test_phase_e_v1_provides_item_2_saturation_evidence` — apply with
    extreme traction; assert `v_active_um2 ≤ v_active_max_bound_um2`
    (saturation built-in to HB#1+#2 convex weight ∈ [0,1])
12. `test_phase_e_v1_provides_item_3_zero_traction_no_response_evidence` —
    apply with zero traction; assert orientation unchanged + V_active
    = 0
13. `test_phase_e_v1_provides_item_4_ecm_side_boundedness_evidence`
    (Y9) — apply 50 steps with rotating traction direction (θ uniform
    [0, 2π]); assert `v_active_um2 ≤ v_active_max_bound_um2` at every
    step AND `multipliers_per_fa` all-1.0 (neutral) at every step;
    docstring states this is ECM-side evidence, NOT closed-loop
    convergence

### Meta + guard tests (3)

14. `test_phase_e_v1_module_doc_does_not_overclaim_full_loop` (Y2
    forward guard) — module docstring + function docstring + return
    class docstring do NOT contain `"full closed-loop"`,
    `"full FA→ECM→FA"`, `"gate satisfied"`, `"Items 1-4 satisfied"`
    strings
15. `test_phase_e_v1_module_does_not_reference_effective_stiffness`
    (Y13) — both string-search (`"effective_stiffness" not in src`)
    AND AST walk (no `Attribute.attr` / `Name.id` reference)
16. `test_phase_e_v1_imports_hb12_constants_no_redefinition` (Y14) —
    `ce.TRACTION_REF_NN_PER_UM2 is cr.TRACTION_REF_NN_PER_UM2`
    (object identity); `ce.K_ORIENT_PER_S is cr.K_ORIENT_PER_S`; no
    surprise uppercase numeric constants

### Exports (1)

17. `test_phase_e_v1_exports_through_both_init` — `PhaseEStepResult`
    + `step_closed_loop_phase_e_v1` importable from BOTH
    `acs.v2.dynamics` AND `acs.v2`; no other new exports

---

## 5. Files

- `docs/v2/v2_phase_e_composition_brief.md` (existing, opening brief,
  commit `7529d3e`)
- `docs/v2/v2_phase_e_composition_locked.md` (this file, source of truth)
- `acs/v2/dynamics/closed_loop_phase_e.py` (NEW):
  - `PhaseEStepResult` dataclass (frozen, slots, 5 fields with
    identity invariant)
  - `step_closed_loop_phase_e_v1` function
  - 5 imports from `ecm_constitutive_response` (2 constants + 1 result
    type + 1 function), 2 imports from `ecm_to_fa_bias`, 2 imports from
    `ecm_lyapunov_metric`, 1 import from `fa_to_ecm_scattering`, 1
    import from `focal_adhesion`, 1 import from `ecm_substrate`
  - NO new constants, NO new failure kinds, NO `effective_stiffness`
    references
- `acs/v2/dynamics/__init__.py` (2 new exports +
  alphabetical reorder of `__all__`)
- `acs/v2/__init__.py` (mirror exports)
- `tests/v2/test_v2_closed_loop_phase_e.py` (NEW, 17 tests per §4)

---

## 6. References

- Brief (opening position):
  `docs/v2/v2_phase_e_composition_brief.md` (commit `7529d3e`)
- Parent locked plan:
  `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- Sister locks (5 hard blockers + Phase D + B1):
  - `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2/v2_phase_d_no_op_scaffolding_locked.md`
  - `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`
  - `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`
  - `docs/v2/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1)
- Verified at source for B1+B2 catches:
  - `acs/v2/focal_adhesion.py:51` (FocalAdhesionState import path)
  - `acs/v2/dynamics/ecm_to_fa_bias.py:138` (ECMToFABiasResult.multipliers_per_fa)
  - `acs/v2/dynamics/fa_to_ecm_scattering.py:170` (HB#3 empty-FA contract)
- Sister implementation precedents:
  - `acs/v2/dynamics/closed_loop_phase_d.py` (Phase D no-op orchestrator)
  - `acs/v2/dynamics/ecm_constitutive_response.py` (HB#1+#2 active law)
- v2 measurement registry (separate concern):
  `acs/v2/metrics.py` (file, not directory)
- Adversarial debate posture: memory
  `feedback_aggressive_design_debate.md`, PI id=809
- Hard Rule 10 inline derivation: memory
  `rule10_unit_derivation_in_docs.md`, Codex id=1234
- Hard Rule 11 wording-boundary meta-test: memory
  `hard_rule_11_wording_boundary_meta_test.md`
- Step 6 pre-commit batch: memory `design_note_pre_commit_batch.md`,
  Codex id=1428

---

## 7. Cross-room dispatch

This file is the design-team input to implementation-work for:

1. impl Claude writes Phase E v1 Sanity Gate doc
   (`docs/v2/v2_phase_e_composition_sanity_gate.md`) from this lock — 6
   sections per §3, with Hard Rule 11 measurement-protocol section as
   the central anchor + the v1-vs-v2 boundary (§0) reproduced verbatim.
2. impl Codex review (6 focus per Codex `id=1531` final SEAL):
   - **Wording boundary**: ECM-side composition evidence ONLY; tests
     are `provides_*_evidence`, NOT `satisfies_*`; meta-test 14
     enforces no satisfaction-claim wording in module
   - **Composition order + identity invariant**: HB#3 → HB#1+#2 → HB#4
     (post-update ECM) → HB#5 (post-update ECM, original scatter);
     `result.updated_ecm is result.orientation_response.updated_ecm`
   - **No wrapper failure kinds**: meta-test 3 source guard (no
     `PhaseE*Error` class, no `failure_kind =` assignment)
   - **API surface verifications**: `from acs.v2.focal_adhesion import
     FocalAdhesionState` (Y15); `result.ecm_to_fa_bias.multipliers_per_fa`
     (Y16); empty-FA `(0, 3)` shape
   - **effective_stiffness exclusion**: meta-test 15 (string + AST
     check)
   - **No constant redefinition**: meta-test 16
     (`ce.TRACTION_REF_NN_PER_UM2 is cr.TRACTION_REF_NN_PER_UM2`)
3. On Sanity Gate PASS: Phase E v1 code commit
   (`acs/v2/dynamics/closed_loop_phase_e.py` new + 2 export updates +
   1 new test file).
4. After commit: design-discussion idle on Phase E v1; **all 5 Phase E
   upstream hard blockers + Phase D no-op + Phase E v1 composition
   sealed**. Next blockers (separate cycles):
   - Phase E v2 (HB#4-active variant) design — requires
     `compute_ecm_to_fa_bias_active` design + lock first
   - Item 5 sweep harness (separate unit)
   - effective_stiffness law decision (separate unit, only if v2
     consumes it)

Rounds 1–4 of the design lock are MCP id 1526–1532.
