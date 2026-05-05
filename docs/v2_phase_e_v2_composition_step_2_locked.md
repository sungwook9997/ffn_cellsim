# V2 Closed-Loop ECM Gate — Phase E v2 Composition Step 2 (Active Bias in the Loop) — Locked Design (B-tier compressed)

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (B-tier compressed cycle
per PI risk-tier policy `id=1653`: 1 design round + 1 review pass +
seal, NOT 3-round full A-tier debate)
**Source unit**: design-discussion `topic=v2-phase-e-v2-composition`,
MCP id 1646–1664 (intent / round 1 / single review / seal ack)
**Brief**: `docs/v2_phase_e_v2_composition_brief.md` (commit `f64dccf`)
**Parent locked plan**: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
**Upstream sister locks**:
- `docs/v2_phase_e_composition_locked.md` (Phase E v1; structural
  twin with HB#4-neutral)
- `docs/v2_hard_blocker_4_active_locked.md` (HB#4-active step 1; the
  active bias law swapped in here)
- `docs/v2_hard_blocker_5_lyapunov_metric_locked.md`
- `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`
- `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
- `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` (HB#4 v1
  neutral; co-resident in `ecm_to_fa_bias.py`)
**Tier**: **B-tier** — pure composition wrapper over locked upstream
contracts. PI policy `id=1653` mandates compressed review (NOT
adversarial multi-round debate). impl-work pane: lock → targeted
tests + one Codex review = milestone (skip Sanity Gate doc).
**PI ratify status**: full delegation per PI id=939/1008.

---

## 0. Scope + Wording Boundary (Q4 = (c) hybrid neutral)

This unit locks **Phase E v2 step 2** — the composition wrapper that
swaps Phase E v1's HB#4-neutral readout for HB#4-active in the
FA→ECM→FA composition.

### Wording boundary (Codex C4 verbatim)

> "Phase E v2 step 2 achieves per-step **composition-structure closure**
> by placing HB#4-active non-neutral ECM→FA bias inside the FA→ECM→FA
> loop. It does **not** by itself establish Items 1-4 multi-step
> satisfaction over time; that requires a sweep harness or
> simulation-engine trajectory."

This preserves the Phase E v1 Y1+Y2 "evidence not satisfaction"
discipline. Test names use `provides_*_evidence` / `anti_collapse_*`,
NOT `satisfies_item_X` (forbidden — would silently overclaim).

### What this unit IS

- Public function `step_closed_loop_phase_e_v2(adhesions, ecm, dt_s, *, k_active, ...)`
- Reuses `PhaseEStepResult` (Q1 = (a); HB#4-active Y15 sister-pattern)
- Lives in `closed_loop_phase_e.py` alongside v1 (Q2 = (a); sister
  with `ecm_to_fa_bias.py` carrying both neutral + active variants)
- Composition order: HB#3 scatter → HB#1+#2 orientation update →
  HB#4-**active** bias readout → HB#5 Lyapunov metric (v1 sister with
  active bias swap)
- Validation order: `ecm.validate() → _validate_k_active(k_active) →
  scatter` (Codex C3; sister with HB#4-active Y12)
- No wrapper failure kinds (Y17 inherited from v1)
- Identity invariant: `result.updated_ecm is
  result.orientation_response.updated_ecm` (Y4 inherited from v1)
- Anti-collapse meta-test enforces FA differentiation under
  anisotropic ECM (Codex C5 strengthened)

### What this unit is NOT

- NOT Items 1-4 satisfaction (multi-step concern; this is per-step
  composition closure only — Q4 = (c) hybrid wording)
- NOT a new result dataclass (Q1 = (a) reuse)
- NOT a separate module (Q2 = (a) modify)
- NOT a public validator surface for `_validate_k_active` (Codex C3:
  same-package private import via internal alias only)

---

## 1. Final Lock Summary

### Function

```python
"""V2 Phase E closed-loop ECM gate composition (v2: active bias step 2).

Structural swap of v1 neutral HB#4 readout for HB#4-active orientation-driven
bias. Composition closure structure satisfied; multi-step Items 1-4
satisfaction remains separate concern (sweep harness / simulation engine).
"""

# Existing closed_loop_phase_e.py — MODIFY to add v2 below v1

from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMToFABiasResult,
    _validate_k_active as _validate_hb4_active_k_active,  # internal alias (Codex C3)
    compute_ecm_to_fa_bias_active,
    compute_ecm_to_fa_bias_neutral,                       # already imported by v1
)


def step_closed_loop_phase_e_v2(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
    dt_s: float,
    *,
    k_active: float,                                      # REQUIRED, no default (Q3 = (b))
    traction_ref_nN_per_um2: float = TRACTION_REF_NN_PER_UM2,
    k_orient_per_s: float = K_ORIENT_PER_S,
) -> PhaseEStepResult:
    """Phase E v2 ECM-side composition step (active HB#4 in the loop).

    Composes:
      1. HB#3 scatter
      2. HB#1+#2 orientation update
      3. HB#4-**active** bias readout (deviatoric Rayleigh, non-trivial multipliers)
      4. HB#5 Lyapunov-like metric

    Validation order (Codex C3, sister with HB#4-active Y12):
      a. ecm.validate()
      b. _validate_k_active(k_active)  — cheap parameter failure first
      c. sub-call sequence (HB#3 → HB#1+#2 → HB#4-active → HB#5)

    No wrapper failure kinds (Y17 inheritance from v1).
    Identity invariant: result.updated_ecm is result.orientation_response.updated_ecm (Y4).
    """
    # Validation order (Codex C3)
    ecm.validate()
    _validate_hb4_active_k_active(k_active)

    # Composition (sister with v1, with HB#4-active swap at step 3)
    traction_density_xy = scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)
    orientation_response = step_ecm_orientation_response(
        ecm, traction_density_xy, dt_s,
        traction_ref_nN_per_um2=traction_ref_nN_per_um2,
        k_orient_per_s=k_orient_per_s,
    )
    # HB#4-**active** swap (vs v1 neutral)
    ecm_to_fa_bias = compute_ecm_to_fa_bias_active(
        adhesions, orientation_response.updated_ecm, k_active=k_active,
    )
    lyapunov_metric = compute_ecm_orientation_lyapunov_metric(
        orientation_response.updated_ecm, traction_density_xy
    )
    return PhaseEStepResult(
        traction_density_xy=traction_density_xy,
        orientation_response=orientation_response,
        ecm_to_fa_bias=ecm_to_fa_bias,
        lyapunov_metric=lyapunov_metric,
        updated_ecm=orientation_response.updated_ecm,  # identity invariant Y4
    )
```

### Forbidden in Phase E v2 step 2

- `test_phase_e_v2_satisfies_item_X` test names (Q4 = (c) hybrid;
  silent overclaim — use `provides_*_evidence` / `anti_collapse_*`
  instead)
- "Phase E v2 satisfaction achieved" / "full closed-loop satisfied" /
  "Items 1-4 satisfied" wording in module docstring, function
  docstring, or test docstrings (meta-test 13 forward guard)
- New `PhaseEV2StepResult` dataclass (Q1 = (a) reuse; HB#4-active Y15
  sister)
- New `closed_loop_phase_e_v2.py` module (Q2 = (a) modify)
- Public export of `_validate_k_active` (Codex C3; internal alias
  only)
- Default value for `k_active` in v2 wrapper signature (Q3 = (b)
  required-no-default sister with HB#4-active Y3)
- Wrapper failure kinds: `PhaseE*Error` class definition or
  `failure_kind = ...` assignment in `closed_loop_phase_e.py` v2
  function body (Y17 inherited; sub-call errors propagate verbatim)
- Anti-collapse test that only checks non-neutrality but not FA
  differentiation (Codex C5 strengthened)
- Recomputing scatter inside HB#5 call (call-order spy test enforces
  HB#5 receives original `traction_density_xy` object — Codex C6)
- Validation order other than `ecm.validate() →
  _validate_k_active(k_active) → scatter` (Codex C3)

---

## 2. Test Catalog (14 tests)

### Composition invariants (4, sister with v1)

1. `test_phase_e_v2_returns_phase_e_step_result_dataclass` — return is
   `PhaseEStepResult` instance with all 5 expected fields populated
2. `test_phase_e_v2_identity_invariant_updated_ecm_is_orientation_response_updated_ecm`
   (Y4 inherited)
3. `test_phase_e_v2_module_does_not_define_wrapper_failure_kinds`
   (Y17 inherited source/meta guard)
4. `test_phase_e_v2_call_order_via_monkeypatch` — verify HB#3 → HB#1+#2
   → HB#4-**active** → HB#5; HB#4-active receives `(adhesions,
   post-update ECM, k_active=k_active)`; HB#5 receives `(post-update
   ECM, original scatter object)` — also enforces Codex C6 traction
   object-identity preservation

### Validation (3, Codex C3)

5. `test_phase_e_v2_invalid_k_active_raises_before_scatter` — valid
   ECM + invalid k_active; assert `_validate_hb4_active_k_active`
   raises BEFORE `scatter_fa_traction_to_ecm_bilinear` is called
   (monkeypatch scatter to track invocation; should NOT be called)
6. `test_phase_e_v2_propagates_hb12_error_stops_before_bias` —
   monkeypatch HB#1+#2 to raise; verify HB#4-active and HB#5 not
   called
7. `test_phase_e_v2_propagates_hb5_error_at_end` — monkeypatch HB#5
   to raise; verify HB#3/HB#1+#2/HB#4-active all called first

### HB#4-active swap (3)

8. `test_phase_e_v2_invokes_active_not_neutral` — monkeypatch both
   `compute_ecm_to_fa_bias_active` and
   `compute_ecm_to_fa_bias_neutral`; assert active is called, neutral
   is NOT
9. `test_phase_e_v2_diagnostics_mechanism_is_deviatoric_rayleigh_orientation`
   — `result.ecm_to_fa_bias.diagnostics_dict["mechanism"] ==
   "deviatoric_rayleigh_orientation"`
10. `test_phase_e_v2_diagnostics_records_k_active` —
    `result.ecm_to_fa_bias.diagnostics_dict["k_active"] == passed value`
    (sweep traceback)

### Anti-collapse (1, Codex C5 strengthened)

11. `test_phase_e_v2_anti_collapse_under_anisotropic_ecm_with_two_perpendicular_fas`
    — THE central new invariant. Concrete spec:
    - Anisotropic ECM `T = [[1.0, 0.0], [0.0, -1.0]]` (already
      traceless → Q = T) at FA position
    - Two FAs at SAME cell-centered position (e.g., domain center)
    - Traction directions `(1.0, 0.0)` and `(0.0, 1.0)`
    - `dt_s = 0.0` (HB#1+#2 exact no-op; isolates v2 active readout)
    - `k_active = 1.0`
    - Math derivation:
      - FA 0: `n = (1, 0)`, `n^T Q n = 1·1·1 + 0 + 0 + 0·(-1)·0 = 1` → multiplier = `exp(+1)`
      - FA 1: `n = (0, 1)`, `n^T Q n = 0 + 0 + 0 + 1·(-1)·1 = -1` → multiplier = `exp(-1)`
    - Assertions:
      - `result.ecm_to_fa_bias.multipliers_per_fa[0, :]` all equal
        `exp(+1)` (3 rates, all same per HB#4-active Y4 scaffolding)
      - `result.ecm_to_fa_bias.multipliers_per_fa[1, :]` all equal
        `exp(-1)`
      - `not np.allclose(result.ecm_to_fa_bias.multipliers_per_fa, 1.0)`
        — neither row collapses to neutral
      - `not np.allclose(result.ecm_to_fa_bias.multipliers_per_fa[0],
        result.ecm_to_fa_bias.multipliers_per_fa[1])` — rows
        differentiated (the FA differentiation guard)

### Composition closure evidence (2, Q4 = (c) hybrid)

12. `test_phase_e_v2_provides_composition_closure_structure_evidence`
    — assert `result.ecm_to_fa_bias.diagnostics_dict["mechanism"]` is
    active (NOT neutral) under non-isotropic ECM; document this as
    **per-step composition closure evidence**, NOT Items 1-4
    satisfaction
13. `test_phase_e_v2_module_doc_distinguishes_composition_closure_from_items_1_4_satisfaction`
    — Q4 (c) forward guard. AST/string check that v2 function +
    module docstring contain the explicit boundary phrase or both of:
    - `"composition-structure closure"` (or close paraphrase)
    - `"does not by itself establish Items 1-4 multi-step satisfaction"`
      (or close paraphrase)
    AND do NOT contain `"Items 1-4 satisfied"`, `"full closed-loop
    satisfied"`, `"satisfies_item"`

### Exports (1)

14. `test_phase_e_v2_exports_through_both_init` —
    `step_closed_loop_phase_e_v2` importable from BOTH
    `acs.v2.dynamics` AND `acs.v2`; no other new exports (no new
    types since Q1 = (a) reuse)

---

## 3. Files

- `docs/v2_phase_e_v2_composition_brief.md` (existing intent note,
  commit `f64dccf`)
- `docs/v2_phase_e_v2_composition_step_2_locked.md` (this file,
  source of truth)
- `acs/v2/dynamics/closed_loop_phase_e.py` (MODIFY):
  - Add internal alias import `from acs.v2.dynamics.ecm_to_fa_bias
    import _validate_k_active as _validate_hb4_active_k_active`
  - Add import `compute_ecm_to_fa_bias_active`
  - Add `step_closed_loop_phase_e_v2` function
  - **NO** new dataclasses; **NO** new failure kinds; **NO**
    `effective_stiffness` references; **NO** local redefinition of HB
    constants
- `acs/v2/dynamics/__init__.py` + `acs/v2/__init__.py` — 1 new export
  (`step_closed_loop_phase_e_v2`)
- `tests/test_v2_closed_loop_phase_e.py` (EXTEND existing v1 test
  file with v2 tests) OR `tests/test_v2_closed_loop_phase_e_v2.py`
  (NEW — impl decision per implementation review)

---

## 4. References

- Brief (intent note):
  `docs/v2_phase_e_v2_composition_brief.md` (commit `f64dccf`)
- PI risk-tier policy: MCP `id=1653` (B-tier compressed protocol for
  composition wrappers)
- Parent locked plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- Sister locks:
  - `docs/v2_phase_e_composition_locked.md` (Phase E v1; structural
    twin)
  - `docs/v2_hard_blocker_4_active_locked.md` (HB#4-active step 1;
    Y3 no-default k_active, Y4 all-3-rates same scaffolding, Y8
    function-scoped guard, Y12 validation order, Y15 result reuse)
  - `docs/v2_hard_blocker_5_lyapunov_metric_locked.md`
  - `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`
  - `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` (HB#4 v1
    neutral; co-resident with active in `ecm_to_fa_bias.py`)
- Verified at source for Codex C5 anti-collapse spec:
  - `acs/v2/dynamics/ecm_to_fa_bias.py:128` (ECMToFABiasResult schema
    with `multipliers_per_fa: np.ndarray`)
  - `acs/v2/dynamics/ecm_constitutive_response.py` (HB#1+#2;
    `dt_s=0.0` is valid no-op per Y15 of HB#1+#2 lock)
- Adversarial debate posture: memory
  `feedback_aggressive_design_debate.md`, PI id=809 (now scoped to
  A-tier only per `id=1653`)
- Hard Rule 11 wording-boundary meta-test: memory
  `hard_rule_11_wording_boundary_meta_test.md` (Q4 = (c) hybrid is
  the natural sister of Phase E v1 Y1+Y2 + HB#4-active Y4)
- Step 6 pre-commit batch: memory `design_note_pre_commit_batch.md`,
  Codex id=1428

---

## 5. Cross-room dispatch (B-tier compressed)

This file is the design-team input to implementation-work for:

1. impl Claude implements `step_closed_loop_phase_e_v2` per §1
   pseudocode + §2 14-test catalog directly. **Skip Sanity Gate doc**
   per PI B-tier policy `id=1653`.
2. impl Codex performs **single review pass** on the resulting
   commit covering:
   - Composition order + identity invariant + no wrapper failure kinds
     (sister with v1 review focus)
   - Anti-collapse test 11 with concrete `T = [[1, 0], [0, -1]]` +
     two perpendicular FAs + dt_s=0 + math `exp(±1)`
   - Validation order `ecm.validate() →
     _validate_hb4_active_k_active(k_active) → scatter` (test 5
     enforces; scatter NOT called when k_active invalid)
   - Q4 (c) hybrid wording in docstrings (test 13 enforces)
   - 1 new export only (`step_closed_loop_phase_e_v2`)
3. On PASS: Phase E v2 step 2 milestone reached.
4. After milestone: design-discussion idle on Phase E v2 step 2.
   Next design entries (per PI north star `id=1653`: v2 production-like
   pilot runner + visual evidence):
   - v2 production-like pilot runner (likely B-tier compressed)
   - Visual evidence dashboard / artifact viewer (likely B-tier
     compressed, possibly C-tier docs)
   - effective_stiffness law / Phase E v3 (A-tier; only if consumer
     needs it)

Compressed cycle MCP id 1646–1664 (intent → round 1 → single review
pass → seal ack).
