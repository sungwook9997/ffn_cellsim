# V2 Phase 1 — Phase D No-Op Scaffolding Sanity Gate

**Date**: 2026-05-05 KST
**Status**: pre-execution Sanity Gate for the Phase D no-op
closed-loop scaffolding integrator. Required by CLAUDE.md "Sanity
Gate Protocol" before the first execution of any new physics /
numerics module. Authored by impl-work Claude after the Phase D
lock (commit `18b5430` / `7786b20`) and Codex round-3 ACK
(`id=1411` design-discussion / `id=1412` implementation-work).
**Source lock**: `docs/v2_phase_d_no_op_scaffolding_locked.md`
(commit `7786b20`).
**Target module**: `acs/v2/dynamics/closed_loop_phase_d.py` (not
yet committed).
**Test catalog**: `tests/test_v2_closed_loop_phase_d.py` (not yet
committed).

This Sanity Gate is the impl-work side's gate before code lands.
Phase E activation is BLOCKED on all 5 Hard Blockers + the
effective_stiffness law decision; nothing in this gate authorizes
any active response or active bias law.

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

The composition wrapper module
`acs/v2/dynamics/closed_loop_phase_d.py` containing exactly:

- `FAToECMResponseResult` dataclass (`traction_density_xy: np.ndarray (nx, ny, 2) float64 nN/μm²` + `updated_ecm: ECMSubstrateState` Phase-D-identity).
- `PhaseDNoOpStepResult` dataclass (`fa_to_ecm: FAToECMResponseResult` + `ecm_to_fa: ECMToFABiasResult` + `updated_ecm: ECMSubstrateState`).
- `step_fa_to_ecm_response(adhesions, ecm) -> FAToECMResponseResult` — FA→ECM leg only (delegates to HB#3 `scatter_fa_traction_to_ecm_bilinear`).
- `step_ecm_to_fa_bias(adhesions, ecm) -> ECMToFABiasResult` — ECM→FA leg only (delegates to HB#4 `compute_ecm_to_fa_bias_neutral`).
- `step_phase_d_no_op(adhesions, ecm) -> PhaseDNoOpStepResult` — composition wrapper (calls FA→ECM first, then ECM→FA only if FA→ECM succeeds).

All three functions are **pure**: no ECM mutation, no FA mutation,
no RNG, no state-label transitions, no scalarization, no Phase E
hook. Object identity `result.updated_ecm is ecm` everywhere.

### Explicitly out of scope (will FAIL gate if introduced)

- Active constitutive response law (Hard Blocker #1).
- Active saturation form (Hard Blocker #2).
- Lyapunov-like metric (Hard Blocker #5).
- effective_stiffness helper or any geometry-field consuming
  shortcut (locked phased plan §3 side-door guard).
- Diagnostics dict on `FAToECMResponseResult` (rejected in lock
  round 1 C1 — typed-contract regression).
- Optional Phase E params (`response_law`, `mapping_law`,
  `bias_params`) — silent activation risk.
- In-place ECM mutation.
- Deep-copy ECM that pretends to be an update.
- New failure_kinds at the wrapper layer (delegate fully to HB#3
  / HB#4 primitives).
- Scalarization of `traction_density_xy` at any layer.
- Any Phase E function or parameter.
- Use of HB#3 `scatter_fa_traction_to_ecm_bilinear` for anything
  other than the FA→ECM leg call site.
- Use of HB#4 `compute_ecm_to_fa_bias_neutral` for anything other
  than the ECM→FA leg call site (sampler-only HB#4
  `sample_ecm_at_fa_positions` is not consumed at the Phase D
  layer — diagnostics propagate via `ECMToFABiasResult.sampled_diagnostics`).

If any of these are introduced, the gate FAILs and the unit halts
to surface the scope creep to PI per CLAUDE.md "Sanity Gate
Failure handling".

---

## 1. Dimensional analysis

The Phase D wrapper introduces **no new unit chain**. All units
delegate to the locked sister-gate analyses:

- HB#3 `scatter_fa_traction_to_ecm_bilinear`:
  `[nN]·[dimensionless]/[μm²] = [nN/μm²]` per
  `docs/v2_fa_to_ecm_scattering_sanity_gate.md` §1.
- HB#4 `compute_ecm_to_fa_bias_neutral`:
  `[dimensionless]` multiplier × `[1/s]` rate = `[1/s]` rate, with
  identity multiplier preserving the rate scale, per
  `docs/v2_ecm_to_fa_bias_sanity_gate.md` §1.

The wrapper's only "operation" is composition: it reads the leg
outputs and packs them into `PhaseDNoOpStepResult` without arithmetic.

### Status

**PASS**. No new dimensional chain. Composition wrapper preserves
the unit chains of both sister gates by construction. No CFL /
stability bound applies (no time integration at the wrapper
layer; single step).

---

## 2. Boundary cases

The wrapper's boundary contract has three locked cases:

### 2.1 Empty FA list (`adhesions = ()` or `[]`)

Per locked plan §4, the empty-FA-list invariants are:

- `fa_to_ecm.traction_density_xy.shape == (nx, ny, 2)`, dtype
  `np.float64`.
- `np.all(fa_to_ecm.traction_density_xy == 0.0)` (exact zero —
  composition of HB#3 empty-list zero scatter and HB#4
  empty-list neutral bias).
- `fa_to_ecm.updated_ecm is ecm` (object identity).
- `ecm_to_fa.multipliers_per_fa.shape == (0, 3)`, dtype
  `np.float64`.
- `ecm_to_fa.rate_names == RATE_NAMES` per HB#4 lock.
- `ecm_to_fa.sampled_diagnostics` is non-`None` empty
  `ECMSampledAtFAs` (HB#4 lock semantics — never `None`,
  verified `acs/v2/dynamics/ecm_to_fa_bias.py` §`if n_fa == 0`
  branch).
- `result.updated_ecm is ecm` (object identity).

Tests 4 + 11 cover these.

### 2.2 Single-FA interior

Wrapper produces:
- `fa_to_ecm.traction_density_xy` = HB#3 bilinear scatter (4-cell
  stencil, weights summing to 1.0 in interior).
- `ecm_to_fa.multipliers_per_fa` = `np.ones((1, 3), dtype=np.float64)`
  per HB#4 neutral.
- `ecm_to_fa.sampled_diagnostics.fa_ids == (fa.adhesion_id,)`.
- `result.updated_ecm is ecm`.

Test 7 covers leg-output equivalence with primitives. Tests 2 + 3
cover identity invariants.

### 2.3 FA at out-of-grid position

The wrapper's deterministic call ordering (locked §3) means:
- HB#3 raises `FAToECMScatteringError(failure_kind="fa_position_outside_ecm_grid")` first; HB#4 is **not** called.
- The same physical out-of-grid position would also fail HB#4 if reached, with `FAToECMBiasError(failure_kind="fa_bias_position_outside_ecm_grid")`. Per the locked HB#3+HB#4 sister-pattern (both consume the inclusive footprint + `_BOUNDARY_TOL_RELATIVE = 1e-12 * max(nx*dx, ny*dy)` tolerance), any geometry-natural case that passes HB#3 but fails HB#4 is a sister-gate-mirror bug, not a propagation scenario (locked §6 test 9 amendment per `7786b20`).

Test 5 covers deterministic ordering. Test 9 covers HB#4 propagation via mocked/stubbed HB#4 failure (not natural geometry).

### Status

**PASS**. All three boundary cases are explicitly locked with
runtime invariants. The validation order
(`len(position) check → unpack → non_finite_fa_position → fa.validate() → out-of-grid raise`)
delegates fully to HB#3 + HB#4 primitives — wrapper introduces no
new validation. Sister-gate-mirror preservation verified at the
brief pre-commit batch (commit `eafbdb9` Step 6) and lock
pre-commit batch (commit `18b5430` Step 6) — both passed.

---

## 3. Conservation invariants

The wrapper preserves all sister-gate conservation laws by
construction (composition without arithmetic):

### 3.1 Pure-function constraint

The locked §1 forbidden list explicitly forbids:
- ECM mutation (in-place or deep-copy substitute).
- FA mutation.
- RNG.
- State-label transitions.
- Scalarization.
- Phase E hook.

Test 3 (`test_phase_d_step_no_input_mutation`) asserts both `ecm`
and FA states are byte-identical before and after the wrapper
call. Test 2 asserts the triple-identity invariant
`result.updated_ecm is fa_to_ecm.updated_ecm is ecm`.

### 3.2 HB#3 scatter conservation (component-wise vector)

HB#3 lock guarantees `Σ_grid density · cell_area == Σ_FA traction`
component-wise to float64 round-off. The wrapper passes `result.fa_to_ecm.traction_density_xy`
through unchanged from `scatter_fa_traction_to_ecm_bilinear`'s return.

Test 7 (`test_phase_d_step_legs_match_primitives`) asserts the
wrapper's leg output equals the primitive output. Test 10
(`test_phase_d_step_scatter_conservation_inherited`) asserts the
HB#3 conservation directly on the wrapper's output.

### 3.3 HB#4 multiplier structural invariant

HB#4 lock guarantees `multipliers_per_fa.shape == (n_fa, 3)`,
dtype `float64`, and `np.all(multipliers_per_fa == 1.0)` in the
neutral default (Phase D). The wrapper passes `result.ecm_to_fa`
through unchanged from `compute_ecm_to_fa_bias_neutral`'s return.

Test 7 covers the leg-primitive match. Test 4 covers empty-FA invariants.

### 3.4 Composition determinism

The wrapper is deterministic by construction:
- Inputs: pure functions of `(adhesions, ecm)` only.
- No global state.
- No RNG.
- HB#3 and HB#4 primitives are themselves deterministic.

Test 2 + Test 3 cover identity / no-mutation. The wrapper's
determinism is the conjunction of HB#3 + HB#4 determinism and
zero-side-effect composition.

### Status

**PASS**. Wrapper introduces zero new conservation laws. All
sister-gate invariants are preserved by construction. Test 7 is
the central sister-gate-mirror conservation test.

---

## 4. Numerical sanity

### 4.1 Float precision

`np.float64` throughout. Wrapper's `traction_density_xy` is the
output of HB#3 (`np.zeros((nx, ny, 2), dtype=np.float64)`).
Wrapper's `multipliers_per_fa` is the output of HB#4
(`np.ones((n_fa, 3), dtype=np.float64)`). No precision conversion
at the wrapper layer.

### 4.2 No new tolerance constants

The wrapper introduces no new numerical tolerance. It inherits
HB#3's `_BOUNDARY_TOL_RELATIVE = 1e-12 * max(nx*dx, ny*dy)` and
HB#4's identical relative tolerance for the inclusive-footprint
check. **Magic-Number Block: zero new tunables** — see §7.

### 4.3 Per-call work

Wrapper per-call work = HB#3 per-call work + HB#4 per-call work +
O(1) dataclass instantiation overhead. Per-FA work is `O(8 + 12)
= O(20)` (HB#3 4-corner stencil × 2 components = 8; HB#4 4-corner
stencil × 3 fields = 12). No new asymptotic complexity.

### 4.4 Stability

No time integration at the wrapper layer; single step. CFL /
explicit-stepping bounds do not apply. Sister gates own those
checks at their respective primitive layers (and the primitives
themselves are single-step, so no CFL there either).

### Status

**PASS**. Float64 throughout; no new tolerance; per-call work is
the additive sum of sister-gate work; no time integration.

---

## 5. Sign / sense check

The wrapper performs no arithmetic on outputs, so sign / sense is
inherited:

- **`traction_density_xy`**: HB#3 lock §5 — output preserves
  input traction sign component-wise (bilinear weights are
  non-negative; renormalized weights inherit; sign of `f_x` /
  `f_y` propagates to output cells).
- **`multipliers_per_fa`**: HB#4 lock §5 — neutral multiplier is
  exactly `1.0` (positive, sign-preserving on rates). Phase D
  semantics is identity on rates.
- **`updated_ecm is ecm`**: not a "sign" but a semantic
  invariant — Phase D performs no ECM update; `is` identity is
  the cleanest no-op proof.

### Status

**PASS**. No wrapper-introduced sign. Sister gates own the
sign-preservation.

---

## 6. Measurement-protocol consistency (Hard Rule 11)

This is the most critical Sanity Gate item for the wrapper, since
the wrapper's primary risk is silent Phase E activation under the
guise of "Phase D no-op".

### 6.1 No scalarization at wrapper layer

The wrapper exposes:
- `traction_density_xy: np.ndarray (nx, ny, 2)` (vector field).
- `multipliers_per_fa: np.ndarray (n_fa, 3)` (per-FA per-rate
  multiplier).
- `updated_ecm: ECMSubstrateState` (full state; identity in
  Phase D).

No magnitude / norm reduction at the wrapper. Phase E response
law (separate function, separate lock, separate Sanity Gate) owns
the reduction choice per Hard Rule 11. This mirrors HB#3 lock's
"no scalarization at scatter time" boundary.

### 6.2 Geometry consistency with HB#3 + HB#4

The wrapper composes two primitives that read the same ECM grid
with the same coordinate convention (cell-centered grid,
inclusive footprint, boundary half-cell truncated stencil with
renormalized weights). No new geometry at the wrapper layer.

A natural geometry-boundary FA position that passes HB#3 but
fails HB#4 (or vice versa) would be a sister-gate-mirror bug —
the locked §6 test 9 amendment (commit `7786b20`) explicitly
calls this out and routes test 9 through a mocked/stubbed HB#4
failure scenario instead.

### 6.3 Runtime meta-test (Hard Rule 11 boundary)

`test_phase_d_step_does_not_satisfy_phase_e_response_law` is the
explicit Hard Rule 11 boundary meta-test. It asserts that the
wrapper does NOT update ECM fields based on traction density —
even with a non-trivial scatter input, the ECM state is
unchanged after the wrapper call (object identity + field-level
byte-identity).

This test fails the moment any future commit:
- adds an active response law to Phase D (silent activation);
- adds a deep-copy + field-update path to the wrapper;
- introduces an effective_stiffness shortcut;
- adds optional Phase E params with non-default values.

### 6.4 Structural meta-tests (no Phase E hook)

Two structural meta-tests enforce the lock's "Phase E hook
structurally absent" guarantee at runtime:

- `test_phase_d_step_no_phase_e_kwargs`: introspects
  `inspect.signature(step_phase_d_no_op)`,
  `inspect.signature(step_fa_to_ecm_response)`, and
  `inspect.signature(step_ecm_to_fa_bias)`; asserts each has
  exactly two positional parameters (`adhesions`, `ecm`) and
  zero optional kwargs. Any future PR adding `response_law=None`
  / `mapping_law=None` / `bias_params=None` fails this.
- `test_phase_d_step_no_diagnostics_dict_field`: introspects
  `FAToECMResponseResult.__dataclass_fields__.keys() ==
  {"traction_density_xy", "updated_ecm"}` exactly. Any future
  PR adding a `diagnostics_dict` field fails this.

### Status

**PASS**. Three layers of measurement-protocol guard:
- Locked §1 forbidden list (text-level).
- §6.3 runtime meta-test (behavior-level).
- §6.4 structural meta-tests (signature / dataclass-shape level).

This is the strongest Hard Rule 11 boundary in the v2 closed-loop
ECM gate so far.

---

## 7. Magic-Number Block check

**Zero new tunables introduced at the wrapper layer.**

The wrapper inherits:
- HB#3 `_BOUNDARY_TOL_RELATIVE = 1e-12` (relative float64
  tie-break tolerance, derivable from float64 precision,
  grid-invariant in relative terms, NOT chosen to fit any test
  target). Documented in `docs/v2_fa_to_ecm_scattering_sanity_gate.md`
  §7.
- HB#4 identical `_BOUNDARY_TOL_RELATIVE = 1e-12` per HB#4 lock
  (sister-gate-mirror with HB#3).

Magic-Number Block 3-test verdict for the wrapper:
1. **Derivable**: PASS by inheritance — wrapper has no constants
   to derive.
2. **Grid-invariant**: PASS by inheritance — wrapper inherits the
   sister-gate relative tolerance.
3. **Fitting**: PASS — wrapper introduces zero numbers, so
   nothing was chosen to fit any target.

### Status

**PASS**. Zero wrapper-introduced magic numbers. Inherits sister
gates' three-test compliance.

---

## 8. Test catalog (~12)

Owned by `tests/test_v2_closed_loop_phase_d.py` (not yet
committed). Each test maps to a locked invariant in
`docs/v2_phase_d_no_op_scaffolding_locked.md` (commit `7786b20`)
§6.

1. `test_phase_d_step_returns_typed_dataclass` — return type is
   `PhaseDNoOpStepResult`, frozen + slotted dataclass with the
   three locked fields (`fa_to_ecm`, `ecm_to_fa`, `updated_ecm`).
2. `test_phase_d_step_object_identity_updated_ecm` —
   `result.updated_ecm is ecm and result.fa_to_ecm.updated_ecm
   is ecm` (triple-identity invariant per locked §1).
3. `test_phase_d_step_no_input_mutation` — input `ecm` field
   arrays and FA states byte-identical before and after the
   wrapper call (forbidden mutation guard).
4. `test_phase_d_step_empty_fa_list` — empty-FA-list invariants
   per locked §4 (zero scatter + (0, 3) multipliers + non-`None`
   empty `ECMSampledAtFAs` + identity ECM).
5. `test_phase_d_step_call_order_deterministic` — out-of-grid
   FA position surfaces as HB#3
   `FAToECMScatteringError(failure_kind="fa_position_outside_ecm_grid")`
   (NOT HB#4); when HB#3 raises, HB#4 is not called (assertable
   via monkeypatch counter on `compute_ecm_to_fa_bias_neutral`).
6. `test_phase_d_step_no_phase_e_kwargs` — structural meta-test
   per §6.4: each of the three functions has exactly two
   positional parameters and zero optional kwargs.
7. `test_phase_d_step_legs_match_primitives` — wrapper's leg
   outputs equal the underlying primitive outputs (no mutation,
   identical traction array, identical multipliers).
8. `test_phase_d_step_propagates_hb3_failure_kinds` — non-finite
   FA position raises `FAToECMScatteringError` with
   `failure_kind="non_finite_fa_position"`.
9. `test_phase_d_step_propagates_hb4_failure_kinds_when_bias_leg_raises`
   — after FA→ECM succeeds, a mocked/stubbed HB#4 failure
   (monkeypatching `step_ecm_to_fa_bias` to raise
   `FAToECMBiasError`) is propagated unchanged. Per locked §6
   test 9 amendment (commit `7786b20`): natural HB#3/HB#4
   geometry pass/fail mismatch is a sister-gate bug, not a
   propagation scenario.
10. `test_phase_d_step_scatter_conservation_inherited` — HB#3
    component-wise sum-over-grid·cell-area equals sum-over-FAs
    holds on the wrapper's output (sister-gate conservation
    inheritance).
11. `test_phase_d_step_does_not_satisfy_phase_e_response_law` —
    Hard Rule 11 runtime meta-test per §6.3: ECM fields
    unchanged even with non-trivial scatter input.
12. `test_phase_d_step_no_diagnostics_dict_field` — structural
    meta-test per §6.4:
    `FAToECMResponseResult.__dataclass_fields__.keys() ==
    {"traction_density_xy", "updated_ecm"}` exactly.

If a 13th test surfaces during code commit (e.g., a regression
caught by Codex review), it is added with explicit lock
reference; the count is not capped at 12.

---

## 9. Gate verdict

| Item | Status | Notes |
|---|---|---|
| §1 Dimensional | **PASS** | Wrapper has no new unit chain; delegates to HB#3 + HB#4 |
| §2 Boundary | **PASS** | Empty-FA + single-FA + out-of-grid all locked; sister-gate-mirror validation order |
| §3 Conservation | **PASS** | Pure composition; no new conservation; HB#3 + HB#4 invariants preserved |
| §4 Numerical | **PASS** | float64 throughout; no new tolerance; additive per-call work |
| §5 Sign | **PASS** | No wrapper arithmetic; sister-gate signs inherited |
| §6 Measurement-protocol | **PASS** | Three-layer Hard Rule 11 guard (text + behavior + structural meta-tests) |
| §7 Magic-Number Block | **PASS** | Zero new tunables |

**Overall**: gate PASS. impl-work is clear to commit
`acs/v2/dynamics/closed_loop_phase_d.py` +
`tests/test_v2_closed_loop_phase_d.py` after Codex review of
this Sanity Gate doc.

### Outstanding before code lands

- Codex review of this Sanity Gate doc (mirrors HB#3 / HB#4
  Sanity Gate review precedent).
- Code commit must reproduce the locked §1 forbidden list at the
  module docstring level (text-level guard layered on top of
  the runtime meta-tests).
- Test catalog ~12 tests must all be present at first commit; no
  "TODO test_X" placeholders.

---

## 10. References

- Phase D lock: `docs/v2_phase_d_no_op_scaffolding_locked.md`
  (commit `18b5430` initial + `7786b20` Y11 correction).
- Phase D opening brief (superseded by lock):
  `docs/v2_phase_d_no_op_scaffolding_brief.md` (commit `eafbdb9`).
- Hard Blocker #3 (FA→ECM scattering) lock + Sanity Gate +
  impl + tests:
  - `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2_fa_to_ecm_scattering_sanity_gate.md`
  - `acs/v2/dynamics/fa_to_ecm_scattering.py`
  - `tests/test_v2_fa_to_ecm_scattering.py`
- Hard Blocker #4 (ECM→FA bias) lock + Sanity Gate + impl +
  tests:
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2_ecm_to_fa_bias_sanity_gate.md`
  - `acs/v2/dynamics/ecm_to_fa_bias.py`
  - `tests/test_v2_ecm_to_fa_bias.py`
- Phased plan: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`.
- Schemas:
  - `acs/v2/ecm_substrate.py` (`ECMSubstrateState`).
  - `acs/v2/focal_adhesion.py` (`FocalAdhesionState`).
- Memory rules informing this gate:
  - `rule10_unit_derivation_in_docs.md`
  - `hard_rule_11_wording_boundary_meta_test.md`
  - `design_note_pre_commit_batch.md` (5-step + Step 6
    sister-gate-mirror baked into Phase D brief / lock /
    Sanity-Gate pre-commit; promotion-to-memory-rule candidacy
    pending post-Sanity-Gate review per Codex `id=1389`).
- Design-discussion thread: `id=1394 → 1395 → 1396 → 1397 →
  1411 → 1412`.
