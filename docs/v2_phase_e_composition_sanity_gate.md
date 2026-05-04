# V2 Closed-Loop ECM Gate — Phase E Composition v1 Sanity Gate

**Date**: 2026-05-05 KST
**Status**: pre-execution Sanity Gate for the Phase E v1
composition locked design (commit `c8b9550`,
`docs/v2_phase_e_composition_locked.md`). Required by CLAUDE.md
"Sanity Gate Protocol" before the first execution of any new
physics / numerics module. Phase E v1 is the first integrated
closed-loop ECM gate composition combining all 5 upstream hard
blockers + Phase D infrastructure into a single per-step
function; Step 6 sister-gate-mirror layered against HB#3 / HB#4
/ Phase D / B1 / HB#1+#2 / HB#5 (the typed-dataclass +
frozen+slots + validation-before-act + AST-meta-test +
wording-boundary sister precedents).
**Source lock**: `docs/v2_phase_e_composition_locked.md` (commit
`c8b9550`).
**Source brief** (superseded by lock):
`docs/v2_phase_e_composition_brief.md` (commit `7529d3e`).
**Parent locked plan**:
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`.
**Target module**: `acs/v2/dynamics/closed_loop_phase_e.py`
(new; not yet committed).
**Test catalog**: `tests/test_v2_closed_loop_phase_e.py` (new;
not yet committed).

This Sanity Gate is the impl-work side's gate before code lands.
Phase E v1 is **ECM-side composition evidence ONLY, NOT full
closed-loop ECM gate satisfaction** (locked §0 wording boundary
Y1 + Y2). Full closed-loop satisfaction remains pending Phase E
v2 (HB#4-active variant) + Item 5 sweep harness OR an explicit
PI-approved decision that ECM-side-only is enough.

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

The new module `acs/v2/dynamics/closed_loop_phase_e.py`
containing per locked §1:

- `PhaseEStepResult` dataclass (frozen+slots, 5 fields):
  `traction_density_xy`, `orientation_response`,
  `ecm_to_fa_bias`, `lyapunov_metric`, `updated_ecm`. Identity
  invariant (Y4): `updated_ecm is
  orientation_response.updated_ecm`.
- `step_closed_loop_phase_e_v1(adhesions, ecm, dt_s, *,
  traction_ref_nN_per_um2, k_orient_per_s) -> PhaseEStepResult`
  — single composition function; HB#4 hard-wired neutral
  (`compute_ecm_to_fa_bias_neutral`); no v2 placeholder/hooks
  (Y3); no wrapper failure kinds.
- 11 imports from sister modules (HB#3 1 + HB#1+#2 4 = 2
  constants + 1 result type + 1 function + HB#4 2 + HB#5 2 +
  ECM substrate 1 + focal adhesion 1; HB#1+#2 constants
  re-imported, NOT redefined per Y14).
- Exports through `acs/v2/dynamics/__init__.py` and
  `acs/v2/__init__.py` (only 2 new symbols per Y14:
  `PhaseEStepResult` + `step_closed_loop_phase_e_v1`).

### Explicitly out of scope (will FAIL gate if introduced)

Per locked §1 Forbidden:

- "Items 1-4 satisfied" / "full closed-loop" / "full
  FA→ECM→FA" / "gate satisfied" wording in module docstring,
  function docstring, or return-class docstrings (Y1 + Y2
  forward guards via meta-test 14).
- Test names containing `satisfies_item_X` (Y2 — must use
  `provides_item_X_*_evidence` instead; locked §4 catalog uses
  the latter exclusively).
- v2 placeholder, enum, hook, strategy, or optional bias law
  parameter (Y3 — silent-activation risk; v2 gets a separate
  lock when HB#4-active exists).
- In-place ECM mutation, adhesion mutation, or copying of
  sub-results (Y4 — composition wrapper is a thin orchestrator).
- Wrapper failure kinds: any `PhaseE*Error` class definition
  or `failure_kind = ...` assignment in `closed_loop_phase_e.py`
  (Y17 — sub-call errors propagate verbatim).
- Reading or referencing `effective_stiffness` anywhere in
  `closed_loop_phase_e.py` (Y8 + Y13 — string + AST guard via
  meta-test 15; the lock doc itself can explain, but the .py
  file cannot).
- Local redefinition of `TRACTION_REF_NN_PER_UM2` or
  `K_ORIENT_PER_S` (Y14 — must be imported from
  `ecm_constitutive_response`; meta-test 16 enforces object
  identity).
- Object identity claim `result.updated_ecm is ecm` (Y11 —
  false; HB#1+#2 always returns a fresh wrapper).
- Use of `result.ecm_to_fa_bias.multipliers` (Y16 — wrong
  field name; use `multipliers_per_fa`).
- Reordering composition (Y5 — locked HB#3 → HB#1+#2 → HB#4
  → HB#5; call-order spy test 4 enforces).
- Recomputing scatter inside HB#5 call (Y5 — must pass the
  original `traction_density_xy` from step 1, not a recomputed
  one).

If any of the above appear in code, the gate FAILs and the unit
halts to surface the scope creep to PI per CLAUDE.md "Sanity
Gate Failure handling".

---

## 1. Dimensional analysis (Hard Rule 10, inline)

Phase E v1 composition introduces **no new unit chain**; all
sub-calls preserve their own locked unit chains:

- HB#3 scatter:
  `traction_density_xy [nN/μm²]` from FA traction (validated by
  HB#3 unit chain at
  `docs/v2_fa_to_ecm_scattering_sanity_gate.md` §1).
- HB#1+#2 active update: `dt_s [s]`, `traction_ref_nN_per_um2
  [nN/μm²]`, `k_orient_per_s [1/s]`; per-cell `K · S · dt`
  dimensionless argument to `expm1`; `orientation_tensor`
  components dimensionless `[-1, 1]` (validated by HB#1+#2
  unit chain at
  `docs/v2_hard_blocker_1_2_constitutive_response_sanity_gate.md`
  §1).
- HB#4 neutral bias:
  `multipliers_per_fa` dimensionless 1.0 (no Rule 10 chain at
  this layer; sampled diagnostics inherit ECM unit chain).
- HB#5 metric: `cell_area [μm²]`, `||T - n⊗n||²_F`
  dimensionless, `V_active_um2 [μm²]` (validated by HB#5 unit
  chain at `docs/v2_hard_blocker_5_lyapunov_metric_sanity_gate.md`
  §1).

The composition wrapper performs **no unit transformation**
between sub-calls; outputs are passed through as typed
dataclasses.

**Magic-Number Block 3-test verdict** (zero new tunables at
composition layer per Y14):

1. **Derivable** ✓: zero new constants. Both
   `TRACTION_REF_NN_PER_UM2` and `K_ORIENT_PER_S` are imported
   from `ecm_constitutive_response`, NOT redefined; meta-test
   16 enforces `ce.TRACTION_REF_NN_PER_UM2 is
   cr.TRACTION_REF_NN_PER_UM2` (object identity = re-import,
   not redefinition).
2. **Grid-invariant** ✓ by inheritance.
3. **Not fitted** ✓ by inheritance.

### Status

**PASS**. No new dimensional chain at composition layer. 5
sub-call unit chains preserved by the composition pseudocode
(no transformation between calls). Magic-Number Block compliant
by inheritance.

---

## 2. Boundary cases

Per locked §3 item 2, plus the validation contracts inherited
from each sub-call:

### 2.1 Empty FA list (`adhesions = ()`)

Per Y7 + Y11 locked behavior + verified at HB#3 source line 170:

- HB#3 returns all-zero traction
  (`np.zeros((nx, ny, 2), dtype=np.float64)`).
- HB#1+#2 with zero traction → `w = 0` everywhere → identity
  update; **`result.updated_ecm` is still a fresh
  ECMSubstrateState** with all 5 array fields equal-but-not-
  identical to input (no shared memory, per HB#1+#2 Y12).
- HB#4 neutral with empty adhesions →
  `multipliers_per_fa.shape == (0, 3)` (verified at HB#4
  source line 138).
- HB#5 with zero traction → `active.any() == False` →
  `V_active = 0.0`, `passive_orientation_magnitude_um2`
  reflects input orientation magnitude.
- **`result.updated_ecm is ecm` is False** (Y11 — fresh wrapper
  from HB#1+#2 even on identity step).
- `result.updated_ecm.orientation_tensor == ecm.orientation_tensor`
  bytewise (no in-place edit).

Tested by test 5 (`test_phase_e_v1_empty_adhesions_chain`)
which cites `fa_to_ecm_scattering.py:170` + asserts `not
np.shares_memory(...)` and `result.updated_ecm is not ecm`.

### 2.2 `dt_s == 0` (HB#1+#2 Y15 valid no-op)

HB#1+#2 lock §1 + Sanity Gate §2.2 already established `dt_s
== 0` is valid no-op (`w = 0` everywhere → identity update).
Phase E v1 inherits: orientation unchanged; HB#5 V_active
reflects current `‖T_old − n⊗n‖²_F` integrated over active
cells; HB#4 neutral unchanged.

### 2.3 `dt_s < 0` / nonfinite / `bool` (HB#1+#2 propagation)

HB#1+#2 raises
`ECMConstitutiveResponseError(failure_kind="dt_invalid")` per
locked HB#1+#2 §1. Phase E v1 propagates this **as-is**
(no wrapper failure kind per Y17).

### 2.4 Single-cell ECM, single FA

Each sub-call handles per its own boundary; composition has no
special-case branching. Test fixtures with `grid_shape = (1,
1)` + 1 FA verify the chain end-to-end without compositional
gotchas.

### 2.5 Composition stop-before-next-leg (Y6 + tests 6-9)

When any sub-call raises, Phase E v1 stops **before** the next
leg:

- HB#3 raises → HB#1+#2 / HB#4 / HB#5 NOT called (test 6).
- HB#1+#2 raises → HB#4 / HB#5 NOT called (test 7).
- HB#4 raises → HB#5 NOT called (test 8).
- HB#5 raises → all upstream legs already completed (test 9).

Each test uses `monkeypatch.setattr(ce, "<sub_call>", _raiser)`
on the imported symbol in `acs.v2.dynamics.closed_loop_phase_e`
(Y12 — monkeypatch where used, not where defined).

### Status

**PASS**. 5 boundary classes (empty FA / dt=0 / dt invalid /
single-cell / stop-before-next-leg) all locked + tested.
Composition wrapper introduces no new boundary handling.

---

## 3. Conservation invariants

### 3.1 Identity invariant (Y4 — locked composition structural
invariant)

`result.updated_ecm is result.orientation_response.updated_ecm`
— **object identity**, NOT just bytewise equality. The
composition wrapper does NOT copy or mutate sub-results; the
HB#1+#2 lock Y12 no-aliasing guarantees `updated_ecm` is fresh
from input, and Phase E v1 passes that reference through to
the result dataclass.

Tested by test 2
(`test_phase_e_v1_identity_invariant_updated_ecm_is_orientation_response_updated_ecm`).

### 3.2 No wrapper conservation (composition delegates fully)

The composition wrapper has NO state, NO arithmetic, NO
internal accumulators. Each sub-call's conservation invariants
propagate through unchanged:

- HB#3 component-wise scatter conservation (per HB#3 lock).
- HB#1+#2 schema invariant preservation (`|T_ij| ≤ 1`,
  symmetry; `updated_ecm.validate()` runtime guard).
- HB#4 multiplier neutral structural invariant
  (all-1.0 `multipliers_per_fa`).
- HB#5 V_active ≥ 0 + bounded by
  `4.5 · active_area_um2`.
- ECM unchanged-fields bytewise equality (HB#1+#2 Y12 — 4
  fields outside `orientation_tensor` are bytewise-equal but
  fresh-copied).

### 3.3 Stop-before-next-leg propagation (Y6)

The composition is **stop-on-error**: any raise in sub-call N
prevents sub-calls N+1 .. final. This is tested explicitly by
tests 6-9 (per §2.5).

### 3.4 No FA mutation, no input ECM mutation

Phase E v1 does NOT update FA state (`adhesions` argument is
read-only). Phase E v1 does NOT mutate input `ecm` (HB#1+#2
returns a fresh `updated_ecm`; the wrapper passes it through).
Test 5 asserts no shared memory between input and result for
all 5 ECM array fields.

### Status

**PASS**. Composition wrapper conserves nothing of its own
(stateless); sub-call conservation invariants propagate
unchanged; identity invariant `updated_ecm is
orientation_response.updated_ecm` is the wrapper's structural
invariant; stop-before-next-leg ensures no partial-update
contamination.

---

## 4. Numerical sanity

### 4.1 Float precision

`np.float64` enforced throughout the chain (inherited from each
sub-call). Composition wrapper introduces no precision changes.

### 4.2 No new tolerance constants (Y14)

The composition wrapper introduces zero new numerical
tolerance / round-off constants. Both
`TRACTION_REF_NN_PER_UM2` and `K_ORIENT_PER_S` are
**imported** from `ecm_constitutive_response`, NOT redefined.
Meta-test 16 enforces object identity:
`ce.TRACTION_REF_NN_PER_UM2 is cr.TRACTION_REF_NN_PER_UM2`
(re-import, not redefinition).

### 4.3 Per-call work

Per-call work = sum of HB#3 + HB#1+#2 + HB#4 + HB#5 work +
O(1) dataclass instantiation. Linear in grid size + FA count.

### 4.4 Stability

No new time integration at composition layer. HB#1+#2 owns the
stability bound (unconditionally stable in `dt` via convex-
weight construction per HB#1+#2 §4.5). Composition wrapper
does not introduce a CFL-equivalent.

### Status

**PASS**. Float64 throughout; no new tolerance; no new time
integration; no new approximation.

---

## 5. Sign / sense check

The composition wrapper performs no arithmetic. Each sub-call
enforces its own sign convention:

- HB#3 scatter sign preserves traction-direction.
- HB#1+#2 convex weight `w ∈ [0, 1]` always; orientation
  `|T_ij| ≤ 1` preserved.
- HB#4 multipliers all `1.0` (positive).
- HB#5 `V_active ≥ 0` always; passive diagnostic `≥ 0`.

The wrapper has no sign decisions of its own.

### Status

**PASS**. No composition-introduced sign. Sub-call sign
conventions inherited.

---

## 6. Measurement-protocol consistency (Hard Rule 11) — central anchor

The locked §6 calls this the central anchor for Phase E v1
because the wording boundary is what prevents silent overclaim
of closed-loop ECM gate satisfaction.

### 6.1 ECM-side wording boundary (Y1, central)

Per locked §0:

- ✅ "Phase E v1 demonstrates ECM-side composition evidence
  for Items 1-4 under prescribed/scattered FA traction and
  HB#4-neutral readout"
- ✅ "Phase E v1 does not close or satisfy the full closed-
  loop ECM gate"
- ❌ FORBIDDEN: "satisfies Items 1-4" / "full closed-loop
  satisfaction" / "full FA→ECM→FA loop satisfied"

The brief's `test_phase_e_satisfies_item_X` naming is
SUPERSEDED by `test_phase_e_v1_provides_item_X_*_evidence`
(Y2). Test names are the contract enforcement layer; meta-test
14 forward-guards against satisfaction-claim wording in module
docstring, function docstring, or return-class docstring.

### 6.2 Domain-matching across sub-calls (HB#5 Y1 inheritance)

HB#5 measurement domain (active cells only) matches HB#1+#2
update domain (active cells only) per HB#5 lock Y1. Phase E v1
preserves this by the locked composition order:

- HB#3 produces `traction_density_xy`.
- HB#1+#2 updates ECM in active cells only (`w=0` elsewhere).
- HB#5 receives the **post-HB#1+#2 ECM** + the **original
  `traction_density_xy`**, so its active-cells mask is
  computed against the same scatter that drove the update.

This is why locked §1 Composition Order item 4 says HB#5 takes
the **original scatter output** (NOT recomputed); call-order
spy test 4 enforces this.

### 6.3 Items 1-4 evidence framing (Y9)

| Closed-loop ECM gate Item | Phase E v1 evidence | Evidence test |
|---|---|---|
| Item 1 (response monotonicity) | ECM-side: monotone in traction strength (single-cell test) | test 10 `provides_item_1_response_monotonicity_evidence` |
| Item 2 (saturation) | ECM-side: bounded by `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 4.5` per active cell; HB#1+#2 convex-weight saturation | test 11 `provides_item_2_saturation_evidence` |
| Item 3 (zero-stimulus zero-response) | ECM-side: zero traction → identity orientation update + V_active = 0 | test 12 `provides_item_3_zero_traction_no_response_evidence` |
| Item 4 (bounded feedback in single-cell loop) | ECM-side: multi-step rotating-traction stays in `[0, v_active_max_bound_um2]` AND HB#4 multipliers stay neutral | test 13 `provides_item_4_ecm_side_boundedness_evidence` |

**Test 13 wording note** (Y9): docstring states this is
ECM-side evidence, NOT closed-loop convergence; no
monotone-decrease claim under varying-target regime.

### 6.4 Source-level meta-tests (Y2 + Y13 + Y14 + Y17)

Four meta-tests forward-guard against future regressions:

- Test 14 (Y2 forward guard) — module docstring + function
  docstring + return class docstring do NOT contain `"full
  closed-loop"`, `"full FA→ECM→FA"`, `"gate satisfied"`,
  `"Items 1-4 satisfied"` strings.
- Test 15 (Y13) — both string-search (`"effective_stiffness"
  not in src`) AND AST walk (no `Attribute.attr` / `Name.id`
  reference); the lock doc explains effective_stiffness in
  prose, but the .py file cannot.
- Test 16 (Y14) — `ce.TRACTION_REF_NN_PER_UM2 is
  cr.TRACTION_REF_NN_PER_UM2` (object identity); `ce.K_ORIENT_PER_S
  is cr.K_ORIENT_PER_S`; no surprise uppercase numeric
  constants.
- Test 3 (Y17) — source/meta guard: AST walk forbids local
  `PhaseE*Error` class definitions; string check forbids
  `failure_kind =` assignments. Docstring explanatory phrases
  (e.g., "No wrapper failure kinds") are allowed per Codex
  round-4 clarification.

### 6.5 Step 6 sister-gate-mirror application

Per memory `design_note_pre_commit_batch.md` 6-step, Step 6
sister-gate-mirror layered against HB#3 / HB#4 / Phase D / B1
/ HB#1+#2 / HB#5 (the typed-dataclass + frozen+slots +
validation-before-act + AST-meta-test sister precedents):

- **Code (validation ordering)**: composition wrapper does NO
  validation of its own; each sub-call validates per its own
  lock (HB#1+#2 has `ecm.validate()` at entry; HB#5 has
  `ecm.validate()` at entry). Two `ecm.validate()` calls per
  Phase E v1 step (HB#1+#2 entry + HB#5 entry on post-update
  ECM).
- **Design (typed dataclass + frozen+slots)**: `PhaseEStepResult`
  follows Phase D / HB#1+#2 / HB#5 / B1 `@dataclass(frozen=True,
  slots=True)` precedent.
- **API surface (`__init__.py` exports)**: only 2 new symbols
  per Y14, mirroring the lean-export discipline. Test 17
  enforces.
- **Failure-kind discipline**: NO error class introduced (Y17
  source guard); sub-call errors propagate verbatim. Different
  from HB#1+#2's `ECMConstitutiveResponseError` because the
  composition wrapper is a thin orchestrator with no
  domain-specific failure modes — same rationale as HB#5's
  no-error-class decision.
- **Wording-boundary meta-test**: Phase E v1 inherits the
  `hard_rule_11_wording_boundary_meta_test` precedent
  (Phase B test 4 + Phase C test 5 pattern + HB#5 Y4 boundary).
  Meta-test 14 + 15 are the wording-boundary meta-tests for
  Phase E v1.

### Status

**PASS**. Six layers of measurement-protocol guard:
- Locked §0 ECM-side wording boundary (text-level).
- Test names `provides_*_evidence` (NOT `satisfies_*`)
  enforce the boundary.
- Meta-test 14 forward-guards against satisfaction-claim
  wording in source.
- Meta-test 15 forward-guards against `effective_stiffness`
  references in source.
- Meta-test 16 forward-guards against constant redefinition.
- Meta-test 3 forward-guards against wrapper failure kinds.
- Step 6 sister-gate-mirror applied at all 4 layers (code /
  design / API surface / failure-kind) + wording-boundary.

---

## 7. Magic-Number Block check

**Zero new tunables introduced at the composition layer.**

The module imports:
- `TRACTION_REF_NN_PER_UM2` from `ecm_constitutive_response`
  (HB#1+#2 lock).
- `K_ORIENT_PER_S` from `ecm_constitutive_response` (HB#1+#2
  lock).

Both are used as default values in `step_closed_loop_phase_e_v1`'s
keyword-only parameters. Meta-test 16 enforces object-identity
re-import.

`MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0` and
`MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 4.5` are HB#5
constants used implicitly (HB#5 returns
`v_active_max_bound_um2 = 4.5 · active_area_um2` via its own
diagnostics field); Phase E v1 does not reference these
directly.

Magic-Number Block 3-test verdict (per locked §1 + Y14):

1. **Derivable** ✓ by inheritance — both imported constants
   are derivable per HB#1+#2 lock §1 derivation chain.
2. **Grid-invariant** ✓ by inheritance.
3. **Not fitted** ✓ by inheritance.

### Status

**PASS**. Zero new tunables; both used constants are
re-imported (object-identity verified by meta-test 16), not
redefined. Meta-test 16 also scans for surprise uppercase
numeric constants and asserts none exist beyond the imports.

---

## 8. Test catalog (17 tests per locked §4)

Owned by `tests/test_v2_closed_loop_phase_e.py` (not yet
committed). Each test maps to a locked invariant in
`docs/v2_phase_e_composition_locked.md` §4.

### Composition invariants (5)

1. `test_phase_e_v1_returns_phase_e_step_result_dataclass` —
   return is `PhaseEStepResult` instance with all 5 expected
   fields populated.
2. `test_phase_e_v1_identity_invariant_updated_ecm_is_orientation_response_updated_ecm`
   (Y4) — `result.updated_ecm is
   result.orientation_response.updated_ecm`.
3. `test_phase_e_v1_module_does_not_define_wrapper_failure_kinds`
   (Y17) — source/meta guard: no `PhaseE*Error` class, no
   `failure_kind =` assignment in module body (AST walk +
   string check).
4. `test_phase_e_v1_call_order_via_monkeypatch` (Y10
   corrected pattern + Y12 monkeypatch where used) — verifies
   HB#3 → HB#1+#2 → HB#4 → HB#5 order; HB#4 receives
   `(adhesions, post-update ECM)`; HB#5 receives `(post-update
   ECM, original scatter)`.
5. `test_phase_e_v1_empty_adhesions_chain` (Y7 + Y11) —
   empty FA list → all-zero scatter (cite
   `fa_to_ecm_scattering.py:170`); orientation bytewise
   unchanged but fresh array (no aliasing); HB#4 result shape
   `(0, 3)`; V_active = 0; identity invariant holds; `result.updated_ecm
   is not ecm`.

### Error propagation (4)

6. `test_phase_e_v1_propagates_hb3_error_stops_before_orient`
   — monkeypatch HB#3 to raise; verify HB#1+#2 / HB#4 / HB#5
   not called.
7. `test_phase_e_v1_propagates_hb12_error_stops_before_bias`
   — monkeypatch HB#1+#2 to raise; verify HB#4 / HB#5 not
   called.
8. `test_phase_e_v1_propagates_hb4_error_stops_before_lyap`
   — monkeypatch HB#4 to raise; verify HB#5 not called.
9. `test_phase_e_v1_propagates_hb5_error_at_end` —
   monkeypatch HB#5 to raise; verify HB#3 / HB#1+#2 / HB#4
   all called first.

### ECM-side evidence (Items 1-4 EVIDENCE, NOT satisfaction) (4)

10. `test_phase_e_v1_provides_item_1_response_monotonicity_evidence`
    — apply with two traction magnitudes; assert larger
    traction → larger `convex_weight` and faster orientation
    alignment.
11. `test_phase_e_v1_provides_item_2_saturation_evidence` —
    apply with extreme traction; assert `v_active_um2 ≤
    v_active_max_bound_um2` (saturation built-in to HB#1+#2
    convex weight ∈ [0, 1]).
12. `test_phase_e_v1_provides_item_3_zero_traction_no_response_evidence`
    — apply with zero traction; assert orientation unchanged
    + V_active = 0.
13. `test_phase_e_v1_provides_item_4_ecm_side_boundedness_evidence`
    (Y9) — apply 50 steps with rotating traction direction
    (θ uniform `[0, 2π]`); assert `v_active_um2 ≤
    v_active_max_bound_um2` at every step AND
    `multipliers_per_fa` all-1.0 (neutral) at every step;
    docstring states this is ECM-side evidence, NOT
    closed-loop convergence.

### Meta + guard tests (3)

14. `test_phase_e_v1_module_doc_does_not_overclaim_full_loop`
    (Y2 forward guard) — module docstring + function docstring
    + return class docstring do NOT contain `"full closed-loop"`,
    `"full FA→ECM→FA"`, `"gate satisfied"`, `"Items 1-4
    satisfied"` strings.
15. `test_phase_e_v1_module_does_not_reference_effective_stiffness`
    (Y13) — both string-search (`"effective_stiffness" not in
    src`) AND AST walk (no `Attribute.attr` / `Name.id`
    reference).
16. `test_phase_e_v1_imports_hb12_constants_no_redefinition`
    (Y14) — `ce.TRACTION_REF_NN_PER_UM2 is
    cr.TRACTION_REF_NN_PER_UM2` (object identity); `ce.K_ORIENT_PER_S
    is cr.K_ORIENT_PER_S`; no surprise uppercase numeric
    constants.

### Exports (1)

17. `test_phase_e_v1_exports_through_both_init` —
    `PhaseEStepResult` + `step_closed_loop_phase_e_v1`
    importable from BOTH `acs.v2.dynamics` AND `acs.v2`; no
    other new exports.

---

## 9. Gate verdict

| Item | Status | Notes |
|---|---|---|
| §1 Dimensional | **PASS** | No new unit chain; 5 sub-call chains preserved by the composition pseudocode |
| §2 Boundary | **PASS** | 5 boundary classes locked + tested (empty FA / dt=0 / dt invalid / single-cell / stop-before-next-leg) |
| §3 Conservation | **PASS** | Composition stateless; identity invariant `updated_ecm is orientation_response.updated_ecm` is the wrapper's structural invariant; stop-before-next-leg ensures no partial-update contamination |
| §4 Numerical | **PASS** | float64 throughout; no new tolerance; no new time integration |
| §5 Sign | **PASS** | No composition-introduced sign |
| §6 Measurement-protocol | **PASS** | 6-layer guard (locked §0 wording + test names + 4 meta-tests + Step 6 sister-gate-mirror); ECM-side evidence framing for Items 1-4 |
| §7 Magic-Number Block | **PASS** | Zero new tunables; both used constants re-imported (object-identity verified) |

**Overall**: gate PASS. impl-work is clear to commit
`acs/v2/dynamics/closed_loop_phase_e.py` + 2 export updates +
`tests/test_v2_closed_loop_phase_e.py` after Codex review of
this Sanity Gate doc, mirroring HB#3 / HB#4 / Phase D / B1 /
HB#1+#2 / HB#5 Sanity Gate review precedent.

### Outstanding before code lands

- Codex review of this Sanity Gate doc per locked §7 6 review
  focus areas:
  1. Wording boundary: ECM-side composition evidence ONLY;
     tests are `provides_*_evidence`, NOT `satisfies_*`;
     meta-test 14 enforces no satisfaction-claim wording in
     module
  2. Composition order + identity invariant: HB#3 → HB#1+#2 →
     HB#4 (post-update ECM) → HB#5 (post-update ECM, original
     scatter); `result.updated_ecm is
     result.orientation_response.updated_ecm`
  3. No wrapper failure kinds: meta-test 3 source guard (no
     `PhaseE*Error` class, no `failure_kind =` assignment)
  4. API surface verifications: `from acs.v2.focal_adhesion
     import FocalAdhesionState` (Y15);
     `result.ecm_to_fa_bias.multipliers_per_fa` (Y16); empty-FA
     `(0, 3)` shape
  5. effective_stiffness exclusion: meta-test 15 (string +
     AST check)
  6. No constant redefinition: meta-test 16
     (`ce.TRACTION_REF_NN_PER_UM2 is cr.TRACTION_REF_NN_PER_UM2`)
- **Code commit module-docstring intent guard** (Codex
  `id=1546` blocker resolution + Codex `id=1549` example-text
  refinement + B1 `id=1458` precedent divergence note): the
  module docstring must describe the Phase E v1 wording
  boundary in **paraphrased** form. Example safe wording per
  Codex `id=1549`: "this composition provides ECM-side
  evidence only; complete bidirectional closure remains future
  work requiring Phase E v2, the Item 5 sweep harness, and PI
  approval." This wording avoids ALL four literal forbidden
  substrings — `"Items 1-4 satisfied"`, `"full closed-loop"`,
  `"full FA→ECM→FA"`, `"gate satisfied"` — which would trip
  meta-test 14 if reproduced in module source. The literal
  forbidden-string list is owned by meta-test 14 + lock doc +
  this Sanity Gate doc + test docstrings/assertions where
  they are explicitly the targets of the string check, NOT
  the module docstring itself. **Divergence note from B1
  `id=1458` precedent**: B1's forbidden list could be
  reproduced verbatim because no string-matching meta-test
  enforced absence of the forbidden phrases in source. Phase
  E v1 has meta-test 14 specifically scanning for those
  phrases, so verbatim reproduction is contradictory;
  paraphrased intent guard with "complete bidirectional
  closure" / "ECM-side evidence only" wording is the
  resolution.
- All 17 tests must pass at first commit; no `TODO test_X`
  placeholders.
- `pytest tests/test_v2_closed_loop_phase_e.py` + combined
  sister-gate regression (HB#3 + HB#4 + Phase D + 6.3a + 6.3b
  + HB#1+#2 + HB#5 + Phase E v1) must PASS at commit time.

---

## 10. References

- Phase E v1 lock:
  `docs/v2_phase_e_composition_locked.md` (commit `c8b9550`).
- Phase E v1 brief (superseded by lock):
  `docs/v2_phase_e_composition_brief.md` (commit `7529d3e`).
- Parent locked plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` (Phase
  E entry — final upstream composition step).
- Sister-gate Sanity Gate precedents:
  - `docs/v2_fa_to_ecm_scattering_sanity_gate.md` (HB#3)
  - `docs/v2_ecm_to_fa_bias_sanity_gate.md` (HB#4)
  - `docs/v2_phase_d_no_op_scaffolding_sanity_gate.md` (Phase D)
  - `docs/v2_focal_adhesion_dynamics_result_typed_sanity_gate.md`
    (B1)
  - `docs/v2_hard_blocker_1_2_constitutive_response_sanity_gate.md`
    (HB#1+#2)
  - `docs/v2_hard_blocker_5_lyapunov_metric_sanity_gate.md`
    (HB#5)
- Verified at source (Y15+Y16 verification + Y7 empty-FA
  contract):
  - `acs/v2/focal_adhesion.py` (FocalAdhesionState — locked
    §1 corrected import path).
  - `acs/v2/dynamics/ecm_to_fa_bias.py`
    (`ECMToFABiasResult.multipliers_per_fa` — locked §1
    corrected field name).
  - `acs/v2/dynamics/fa_to_ecm_scattering.py` (HB#3 empty-FA
    contract).
- Sister implementation precedents:
  - `acs/v2/dynamics/closed_loop_phase_d.py` (Phase D no-op
    orchestrator).
  - `acs/v2/dynamics/ecm_constitutive_response.py` (HB#1+#2
    active law).
  - `acs/v2/dynamics/ecm_lyapunov_metric.py` (HB#5 metric).
- Memory rules informing this gate:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, per Codex `id=1428`).
  - `rule10_unit_derivation_in_docs.md`.
  - `hard_rule_11_wording_boundary_meta_test.md`.
  - `feedback_aggressive_design_debate.md`.
- Design-discussion thread: `id=1526` → `1532` (4-round
  adversarial lock); `id=1533` Codex SEAL ack; `id=1535`
  terminal ack.
- impl-work review thread (this Sanity Gate cycle): `id=1542`
  (Codex state-sync) → impl-work Sanity Gate doc commit (this
  file) → Codex impl-work review → code commit.
