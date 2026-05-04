# V2 Closed-Loop ECM Gate — HB#1+#2 ECM Orientation Constitutive Response Sanity Gate

**Date**: 2026-05-05 KST
**Status**: pre-execution Sanity Gate for the HB#1+#2 merged
locked design (commit `7b1d3cf`,
`docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`).
Required by CLAUDE.md "Sanity Gate Protocol" before the first
execution of any new physics / numerics module. This is the
**first Phase E active constitutive law** Sanity Gate — no prior
Phase E sister gate exists, so Step 6 sister-gate-mirror applies
to HB#3 / HB#4 / Phase D / B1 (which are pre-Phase-E units) at
the discipline-mirror layer (typed dataclasses, frozen+slots,
no aliasing, validate-before-return, deterministic ordering)
rather than as a direct algebra sister.
**Source lock**: `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`
(commit `7b1d3cf`).
**Source brief** (superseded by lock):
`docs/v2_hard_blocker_1_constitutive_direction_brief.md` (commit
`4950302`).
**Parent locked plan**:
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` (Phase E
entry, HB#1+#2 merged).
**Target module**: `acs/v2/dynamics/ecm_constitutive_response.py`
(new; not yet committed).
**Test catalog**: `tests/test_v2_ecm_constitutive_response.py`
(new; not yet committed).

This Sanity Gate is the impl-work side's gate before code lands.
HB#1+#2 is a **shape-only orientation update** under
instantaneous traction stimulus: the law is exact-exponential
convex (`T_new = (1−w)·T_old + w·n⊗n` with `w = -expm1(-K·S·dt)`),
so saturation is built into the convex weight. Phase E
activation remains BLOCKED on HB#5 (Lyapunov metric) +
effective_stiffness law decision; this gate does NOT authorize
Phase E composition with the rest of the closed-loop ECM gate.

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

The new module `acs/v2/dynamics/ecm_constitutive_response.py`
containing per locked §1:

- 3 module constants (literature-pinned, Magic-Number Block
  compliant):
  - `TRACTION_REF_NN_PER_UM2: Final[float] = 1.0` (Munevar 2001
    + Tan 2003; `1 nN/μm² = 1 kPa` inline derivation).
  - `TAU_ALIGN_RANGE_S: Final[tuple[float, float]] = (3600.0, 7200.0)`
    (Hall 2016 + Notbohm 2015).
  - `K_ORIENT_PER_S: Final[float] = 1.964e-4` (geometric midpoint
    of `TAU_ALIGN_RANGE_S`; derivable, `1 / sqrt(3600 * 7200)`).
- 2 typed dataclasses (frozen, slots):
  - `ECMOrientationResponseDiagnostics` (13 fields covering
    cell-count, input-parameter traceback, Rule 10 bridge,
    convex-weight stats, distance-to-target stats, schema
    sanity).
  - `ECMOrientationResponseResult` (`updated_ecm: ECMSubstrateState`
    + `diagnostics: ECMOrientationResponseDiagnostics`).
- 1 pure function `step_ecm_orientation_response(ecm,
  traction_density_xy, dt_s, *, traction_ref_nN_per_um2,
  k_orient_per_s) -> ECMOrientationResponseResult`.
- Locked error class **`ECMConstitutiveResponseError(ValueError)`**
  with machine-readable `failure_kind` attribute (Codex
  `id=1488` Step 6 sister-gate-mirror with HB#3
  `FAToECMScatteringError` + HB#4 `FAToECMBiasError`). Owned
  failure kinds: `dt_invalid`, `traction_ref_invalid`,
  `k_orient_invalid`, `traction_shape_mismatch`,
  `non_finite_traction`. Inherited `ecm.validate()` failures
  remain plain schema `ValueError` (test 20).
- Exports through `acs/v2/dynamics/__init__.py` and
  `acs/v2/__init__.py` (**7 new symbols** per Codex `id=1488`
  6→7 expansion to include the locked error class: error + 3
  constants + 2 dataclasses + 1 function).

### Explicitly out of scope (will FAIL gate if introduced)

Per locked §1 Forbidden:

- Update of `stiffness_kpa`, `fiber_density`, `ligand_density`,
  or `accumulated_traction_nNs_per_um2` (orientation-only first
  law, Y1).
- Consumption of `accumulated_traction_nNs_per_um2` as driver
  (instantaneous traction only, Y6).
- ε / epsilon as physical or numerical regularizer in direction
  computation (branch-defined zero-traction only, Y4).
- Hard cap or clamp on `orientation_tensor` components
  (saturation built into convex weight, Y10).
- PSD / eigenvalue / trace claims beyond the schema's
  componentwise `|T_ij| ≤ 1` invariant (Y5).
- Aliasing of unchanged arrays into the returned ECM (Y12 —
  all 5 array fields must be fresh copies).
- `Generic[...]` parameterization or preemptive `# type: ignore`
  on the result/diagnostics dataclasses (B1 sister-precedent).
- Silent zero-K no-op (`k_orient_per_s == 0` raises; Y11 — no
  silent path).
- Mechanosensing shortcut: any future downstream consumer that
  uses `orientation_tensor` without explicit decision on
  whether/how to weight by `fiber_density` is disallowed by
  the **latent-orientation contract** (Y7).

If any of the above appear in code, the gate FAILs and the unit
halts to surface the scope creep to PI per CLAUDE.md "Sanity
Gate Failure handling".

---

## 1. Dimensional analysis (Hard Rule 10, inline)

The full unit chain is locked at `7b1d3cf` §1. Reproduced inline
here per the `rule10_unit_derivation_in_docs` memory rule (every
Rule 10 unit comparison must be inline-derivable in the cited
markdown):

- **Constant `TRACTION_REF_NN_PER_UM2` literature equivalence**:
  `1.0 nN/μm² = (1 × 10⁻⁹ N) / (1 × 10⁻¹² m²) = 1 × 10³ Pa = 1.0 kPa` ✓
  — matches collagen substrate FA traction order from Munevar
  2001 + Tan 2003.
- **Stimulus normalization**:
  `traction_norm [nN/μm²] / TRACTION_REF_NN_PER_UM2 [nN/μm²] =
  S [-]` dimensionless ✓.
- **Exponential argument**:
  `K_ORIENT_PER_S [1/s] · S [-] · dt_s [s] = [-]` dimensionless ✓
  (input to `np.expm1`).
- **Convex weight**:
  `w = -np.expm1(-K · S · dt) ∈ [0, 1]` dimensionless ✓.
- **Convex update**:
  `(1 − w) · T_old + w · T_target → orientation_tensor` —
  dimensionless components in `[-1, 1]` (schema invariant;
  closed under componentwise convex combination).
- **`k_orient_per_s` literature derivation**:
  `K_ORIENT_PER_S = 1 / √(τ_low · τ_high) = 1 / √(3600 × 7200 s²)
  ≈ 1.964 × 10⁻⁴ /s` (geometric midpoint, Y9).

**Magic-Number Block 3-test verdict** (zero new tunables beyond
the locked 3 constants):

1. **Derivable** ✓:
   - `TRACTION_REF_NN_PER_UM2 = 1.0` from collagen substrate
     literature (Munevar 2001 + Tan 2003);
   - `TAU_ALIGN_RANGE_S = (3600, 7200)` from collagen alignment
     literature (Hall 2016 + Notbohm 2015);
   - `K_ORIENT_PER_S = 1.964e-4` from geometric midpoint of the
     above range.
2. **Grid-invariant** ✓: no `dx`, `dy`, `dt`, or `grid_n`
   dependence in any constant.
3. **Not fitted** ✓: chosen before any A/A₀ comparison or test
   target gating.

### Status

**PASS**. Full unit chain inline-derived per Rule 10. All 3
constants Magic-Number Block compliant. No new tolerance,
threshold, or per-call numerical regularizer.

---

## 2. Boundary cases

Per locked §3 item 2, plus the 9 validation/failure contracts
locked at Y11 + Y15:

### 2.1 Empty / zero traction (`adhesions = 0` upstream, traction array all zeros)

`traction_norm == 0` everywhere → `nonzero` mask all False →
`n_xy = 0` (preserved zero from `np.zeros_like`) →
`T_target = 0` → `S = 0` → `w = -expm1(0) = 0` everywhere →
`orientation_new = orientation_old` (componentwise convex
combination with `w=0`) → ECM byte-identical except all 5 arrays
fresh-copied per Y12.

Diagnostics: `n_nonzero_cells = 0`, all `*_nonzero` fields = 0,
`max_convex_weight = 0`, `max_orientation_delta_frobenius = 0`.

Tested by test 1 `test_orientation_response_zero_traction_is_identity`.

### 2.2 `dt_s == 0` (Y15 valid no-op)

`w = -expm1(-K · S · 0) = -expm1(0) = 0` everywhere → identity
update + zero-weight diagnostics. **No special branch needed**;
the convex update naturally collapses to identity.

Tested by test 2 `test_orientation_response_dt_zero_is_identity_no_op`.

### 2.3 `dt_s < 0` / bool / nonfinite (rejected per Y11)

- `dt_s < 0`: rejected as nonsensical (orientation cannot
  rewind under causal stimulus).
- `dt_s` of type `bool`: rejected per the open-loop precedent
  `acs/v2/dynamics/ecm_open_loop.py` `_validate_dt`.
- `dt_s` non-finite (`nan`, `±inf`): rejected.

Tested by tests 3-5 (`test_orientation_response_dt_negative_raises`,
`test_orientation_response_dt_bool_raises`,
`test_orientation_response_dt_nonfinite_raises`).

### 2.4 `traction_ref_nN_per_um2` and `k_orient_per_s` strict-positive

Both must be strictly `> 0` and finite. `traction_ref ≤ 0`
violates the dimensional-stimulus interpretation. `k_orient == 0`
would silently produce an identity update for all `dt > 0`,
which is a forbidden silent zero-K no-op (Y11) — caller must use
an explicit disabled-state mechanism if needed.

Tested by tests 6 + 7
(`test_orientation_response_traction_ref_nonpositive_raises`,
`test_orientation_response_k_orient_nonpositive_raises`).

### 2.5 `traction_density_xy` shape mismatch / non-finite

Shape must be `(*ecm.grid_shape, 2)`; mismatch raises. Non-finite
values (`nan`, `±inf`) raise.

Tested by tests 8 + 9.

### 2.6 `traction_ref → ∞` / `K_ORIENT → ∞` / `dt → ∞` (asymptotic limits)

- `traction_ref → ∞` ⇒ `S → 0` ⇒ `w → 0` ⇒ identity update
  (no NaN).
- `K_ORIENT → ∞` ⇒ same as `dt → ∞` ⇒ `w → 1` everywhere
  `S > 0` ⇒ `T → T_target` immediately (asymptotic-bound is
  the saturation, no NaN).
- `dt → ∞` ⇒ same. `w` is `1 − exp(-∞) = 1` cleanly.

These limits are not run as tests (production callers won't hit
them), but the algebra handles them by the convex-weight
construction without special-casing.

### Status

**PASS**. Six boundary classes locked + 7 of 9 covered by
explicit tests (the asymptotic limits in §2.6 are
algebra-handled, not test-required). No NaN paths.

---

## 3. Conservation invariants

### 3.1 ECM schema invariant preservation (locked §3 item 3)

The convex update `T_new = (1−w)·T_old + w·n⊗n` preserves the
ECM schema invariant under the contract **"valid input ECM →
valid output ECM"** (Codex `id=1486` silent-heal-path closure).
Invalid input ECM is **rejected** at function entry by
`ecm.validate()` per Step 6 sister-gate-mirror with HB#3 / HB#4
(both call `ecm.validate()` at function entry as local dynamics
precedent), NOT "healed" by the convex update toward
`T_target` under high traction × dt.

Under valid input, the proof:

- Symmetry: `n⊗n` is symmetric by construction; `T_old` is
  symmetric (input was validated). Convex combination of
  symmetric matrices is symmetric. ✓
- Componentwise bound: schema enforces `|T_ij| ≤ 1`. `n⊗n` has
  components in `[-1, 1]` (since `|n_i| ≤ 1` after
  normalization). `T_old` components in `[-1, 1]` per validated
  input. Convex combination preserves the box `[-1, 1]^4`. ✓
- `updated_ecm.validate()` is called before return as the
  output-side guard, catching any float-rounding drift (e.g.,
  pushing a component to `1 + 1e-17`).

**Two-sided validation contract** (per Codex `id=1486`):

| Validation | When | Purpose |
|---|---|---|
| `ecm.validate()` | function entry | reject invalid input — closes silent-heal path |
| `updated_ecm.validate()` | before return | catch float-rounding drift in output |

Tested by test 12
(`test_orientation_response_componentwise_bound_preserved`) +
test 13 (`test_orientation_response_symmetry_preserved`) + test
20 (`test_orientation_response_invalid_input_ecm_rejected_before_update`,
new per Codex `id=1486`).

### 3.2 Unchanged-fields bytewise equality (Y12, locked §3 item 3)

`stiffness_kpa`, `ligand_density`, `fiber_density`, and
`accumulated_traction_nNs_per_um2` are NOT updated by this law.
The wrapper copies them via `np.array(..., copy=True)` so the
returned ECM has fresh arrays (Y12 no-aliasing) but their values
are bytewise identical to input.

Tested by test 11
(`test_orientation_response_unchanged_fields_bytewise_equal`).

### 3.3 No-aliasing contract (Y12)

All 5 array fields of the returned `ECMSubstrateState` must be
fresh copies (no shared memory with input). This prevents
contamination via downstream mutation of the returned ECM
affecting the input.

Tested by test 10
(`test_orientation_response_no_aliasing_all_arrays`) using
`not np.shares_memory(...)` for all 5 fields.

### 3.4 Distance-to-target monotonicity (Y2)

For a fixed `T_target` (constant traction direction), per-step
`‖T_new − T_target‖_F ≤ ‖T_old − T_target‖_F` because the
convex update with `w ∈ [0, 1]` and a fixed second argument is
contractive toward `T_target` (and equal-to-target iff `w = 0`
or `T_old == T_target`).

This is **per-step under fixed target** (NOT global, NOT under
varying `T_target`). The wording is locked at Y11.

Tested by test 14
(`test_orientation_response_distance_to_target_monotone_per_step`).

### 3.5 Convergence to target under constant traction (Y2 corollary)

Many steps with constant traction → `w` accumulates → `T → T_target`.
Numerically: `‖T_N − T_target‖_F → 0` as N grows.

Tested by test 15
(`test_orientation_response_uniform_traction_converges_to_target`).

### Status

**PASS**. Convex algebra preserves schema invariant by
construction; `updated_ecm.validate()` is the runtime guard;
unchanged-fields bytewise equality + no-aliasing tested
explicitly; distance-to-target monotonicity is the contractive
property of the convex map under fixed target.

---

## 4. Numerical sanity

### 4.1 `expm1` numerical stability (Y10)

`-np.expm1(-x)` is mathematically identical to `1 - np.exp(-x)`
but numerically stable for small `x`. The locked default
parameters give `K · S · dt ≈ 2e-4 · 1 · 60 = 1.2e-2` for
typical Phase 1 timestep — well within the regime where
`exp(-x) ≈ 1 - x + x²/2 - …` and direct `1 - exp(-x)` suffers
catastrophic cancellation.

`expm1(-x)` is the canonical `numpy` primitive for this; no
separate implementation.

### 4.2 Float precision

`np.float64` everywhere. Validators enforce dtype consistency on
`traction_density_xy` (must be array-like coercible to
`float64`). All array fields in `ECMSubstrateState` are already
`float64` per schema.

### 4.3 No new tolerance constants

The only new constants are the 3 locked literature parameters
(`TRACTION_REF_NN_PER_UM2`, `TAU_ALIGN_RANGE_S`,
`K_ORIENT_PER_S`). Zero numerical-tolerance / round-off
constants introduced (the lock explicitly drops the `ε` from
the brief — Y4).

### 4.4 Per-call work

Per-call work = 5 array copies (O(grid_size)) + 3 reductions
(traction_norm, mask, `T_target = einsum`) + 1 convex update
+ diagnostics aggregations. Linear in `grid_size`; no
asymptotic surprise.

### 4.5 Stability (CFL-equivalent)

The convex combination is **unconditionally stable** in `dt`
because `w ∈ [0, 1]` regardless of `K · S · dt`. The bound is
the convex-weight construction itself — no separate CFL
condition.

This is a key reason the lock chose exact-exponential convex
over a linear-in-dt update like `T_new = T_old + dt · K · S ·
(T_target − T_old)`: the linear form has CFL `K · S · dt < 1`
for stability, while the exact-exponential is stable for any
`dt`. (The linear form would also need a separate hard cap to
prevent overshoot, which the lock §1 forbids.)

### Status

**PASS**. `expm1` chosen for stability near zero argument;
float64 throughout; no new tolerance; unconditionally stable
in `dt` via convex-weight construction.

---

## 5. Sign / sense check

### 5.1 Convex weight `w ≥ 0`

`w = -np.expm1(-K · S · dt)` with `K, S, dt ≥ 0` (validated)
implies `-K · S · dt ≤ 0`, so `np.expm1(...) ≤ 0`, so
`w = -np.expm1(...) ≥ 0`. Combined with `w ≤ 1` (since
`np.expm1(-x) ≥ -1` for `x ≥ 0`), we have `w ∈ [0, 1]`.

### 5.2 `T_target = n⊗n` PSD-ish (Y5 weakened)

`n⊗n` is rank-1 with eigenvalues `(|n|², 0)`, so PSD with trace
`|n|²`. Schema does NOT enforce PSD on `orientation_tensor`,
only componentwise bound + symmetry — but `n⊗n` happens to be
PSD by construction. The Sanity Gate makes no PSD claim per
Y5; the eigenvalue/trace property is a side effect of the
construction, not a contract.

### 5.3 Distance-to-target monotone decreasing per step (Y2)

`‖T_new − T_target‖_F = (1−w)·‖T_old − T_target‖_F` (convex
combination algebra, with fixed target). For `w ∈ [0, 1]`,
`(1−w) ∈ [0, 1]`, so `‖T_new − T_target‖_F ≤ ‖T_old − T_target‖_F`.

This is the per-step Lyapunov-like property. It is **NOT** the
HB#5 closed-loop Lyapunov metric (which requires a Phase E
composition with FA dynamics), but it is the algebraic
monotonicity that HB#5 will inherit at the constitutive layer.

### Status

**PASS**. Convex weight non-negative + bounded by construction;
distance-to-target monotone decreasing per step under fixed
target. No sign ambiguity introduced.

---

## 6. Measurement-protocol consistency (Hard Rule 11)

The measurement protocol locked at §3 item 6: orientation tensor
is **per-cell pointwise**; the law is strictly local; the
algebraic proof (convex combination preserves bound + symmetry +
distance-to-target) holds at every cell independently.

### 6.1 No off-proof points

Unlike the Stage 1a v13 episode (`_build_curvature` band-averaged
measurement included off-peak cells), this law has no spatial
coupling. Every cell's `T_new` depends only on its own `T_old`,
its own `n_stim`, the global scalar `K · S · dt`, and the cell's
own `S`. No averaging window, no kernel, no neighbor cells.

### 6.2 Diagnostics fields capture Rule 10 bridge (Y13)

`max_traction_norm_nN_per_um2` + `max_dimensionless_stimulus`
are explicitly captured to make the Rule 10 unit chain
inspectable at every step. Sweep diagnostics can verify the
bridge `[nN/μm²] / [nN/μm²] = [-]` is preserved across runs.

### 6.3 Diagnostics fields capture input parameters (Y13)

`traction_ref_nN_per_um2` + `k_orient_per_s` are recorded in
the diagnostics so a sweep over either parameter can be
reconstructed from the output stream alone (no need to consult
the caller's config).

### 6.4 Schema-sanity diagnostics

`symmetry_residual_max` and `componentwise_bound_residual_max`
are computed at every step. If either exceeds float-rounding
tolerance, the corresponding `updated_ecm.validate()` would
raise — but the diagnostics surface the residual magnitude
even when validation passes, allowing a sweep to detect
near-violations.

### 6.5 Step 6 sister-gate-mirror application

Per memory `design_note_pre_commit_batch.md` 6-step, Step 6
sister-gate-mirror layered against HB#3 / HB#4 / Phase D / B1
(the only typed-dataclass pre-Phase-E precedents):

- **Code (validation ordering)**: HB#1+#2 has different inputs
  (no `position_um_xy`, no FA list); the mirror is at the
  schema-level: **`ecm.validate()` first** (per HB#3/HB#4 local
  dynamics precedent + Codex `id=1486` silent-heal-path closure)
  → `_validate_finite_*` → `_validate_dt_*` →
  `_validate_positive_*` → `_validate_*_shape` before any
  computation, mirroring HB#3/#4's "validate-before-act"
  discipline.
- **Design (typed dataclass + frozen+slots)**: HB#1+#2 follows
  B1 / HB#4 / Phase D `@dataclass(frozen=True, slots=True)`
  precedent for `ECMOrientationResponseDiagnostics` +
  `ECMOrientationResponseResult`.
- **API surface (`__init__.py` exports)**: **7 new symbols**
  at both `acs/v2/__init__.py` and `acs/v2/dynamics/__init__.py`
  with `__all__` inclusion (per Codex `id=1488` 6→7 expansion
  to include the locked `ECMConstitutiveResponseError`). Test 18
  (`test_orientation_response_exports_through_both_init`)
  enforces.
- **Failure-kind discipline** (locked per Codex `id=1488`):
  HB#1+#2 introduces `ECMConstitutiveResponseError(ValueError)`
  with machine-readable `failure_kind` attribute, mirroring
  HB#3 `FAToECMScatteringError` + HB#4 `FAToECMBiasError`
  pattern. Owned failure kinds are exhaustively enumerated:
  `dt_invalid` (covers negative / bool / nonfinite),
  `traction_ref_invalid` (nonpositive / nonfinite),
  `k_orient_invalid` (nonpositive / nonfinite — strict
  no-zero-no-op), `traction_shape_mismatch`,
  `non_finite_traction`. Inherited `ecm.validate()` failures
  remain plain schema `ValueError` and are NOT wrapped in
  `ECMConstitutiveResponseError` (Codex `id=1488` —
  schema-level errors stay schema-level, not masked as
  constitutive-law-level errors). Tests 3-9 assert the typed
  error + specific failure_kind; test 20 asserts plain
  `ValueError` for the inherited schema rejection path.

### 6.6 Y12-guard meta-test (extended per Codex `id=1486`)

Test 19
(`test_orientation_response_uses_current_schema_not_stale_naming`)
is a static check via `inspect.getsource()` regex or AST walk
that the implementation:

- does NOT reference `dx_um` (stale field name);
- does NOT pass `grid_shape` to the `ECMSubstrateState`
  constructor (it is a derived property, not a field);
- calls **both** `ecm.validate()` at function entry (Codex
  `id=1486` silent-heal-path closure) AND
  `updated_ecm.validate()` before return (Y12 output guard).

This locks Codex's Y12 + `id=1486` catches into runtime
evidence. The two-sided validation contract is enforced
statically.

### Status

**PASS**. Measurement protocol is strictly local (no off-proof
points); Rule 10 bridge fields captured in typed diagnostics;
input parameters recorded for sweep traceback; schema-sanity
residuals surfaced; Step 6 sister-gate-mirror applied at all
four layers; Y12-guard locked into a runtime meta-test.

---

## 7. Test catalog (~20 tests per locked §4 + Codex `id=1486` test 20)

Owned by `tests/test_v2_ecm_constitutive_response.py` (not yet
committed). Each test maps to a locked invariant in
`docs/v2_hard_blocker_1_2_constitutive_direction_locked.md` §4.

1. `test_orientation_response_zero_traction_is_identity` (§2.1)
2. `test_orientation_response_dt_zero_is_identity_no_op` (§2.2;
   Y15)
3. `test_orientation_response_dt_negative_raises` (§2.3) —
   asserts `ECMConstitutiveResponseError` with
   `failure_kind="dt_invalid"`
4. `test_orientation_response_dt_bool_raises` (§2.3) —
   asserts `ECMConstitutiveResponseError` with
   `failure_kind="dt_invalid"`
5. `test_orientation_response_dt_nonfinite_raises` (§2.3) —
   asserts `ECMConstitutiveResponseError` with
   `failure_kind="dt_invalid"`
6. `test_orientation_response_traction_ref_nonpositive_raises`
   (§2.4) — asserts `ECMConstitutiveResponseError` with
   `failure_kind="traction_ref_invalid"`
7. `test_orientation_response_k_orient_nonpositive_raises`
   (§2.4; strict, Y11) — asserts
   `ECMConstitutiveResponseError` with
   `failure_kind="k_orient_invalid"`
8. `test_orientation_response_traction_shape_mismatch_raises`
   (§2.5) — asserts `ECMConstitutiveResponseError` with
   `failure_kind="traction_shape_mismatch"`
9. `test_orientation_response_traction_nonfinite_raises`
   (§2.5) — asserts `ECMConstitutiveResponseError` with
   `failure_kind="non_finite_traction"`
10. `test_orientation_response_no_aliasing_all_arrays` (§3.3;
    Y12) — `not np.shares_memory(...)` for all 5 array fields
11. `test_orientation_response_unchanged_fields_bytewise_equal`
    (§3.2)
12. `test_orientation_response_componentwise_bound_preserved`
    (§3.1)
13. `test_orientation_response_symmetry_preserved` (§3.1)
14. `test_orientation_response_distance_to_target_monotone_per_step`
    (§3.4; Y2)
15. `test_orientation_response_uniform_traction_converges_to_target`
    (§3.5)
16. `test_orientation_response_diagnostics_typed_dataclass` (B1
    sister meta-test: `isinstance` + `__getitem__` raises)
17. `test_orientation_response_diagnostics_records_input_params`
    (§6.3; Y13)
18. **`test_orientation_response_exports_through_both_init`**
    (§6.5; sister-gate-mirror API surface — **7 new symbols**
    per Codex `id=1488` 6→7 expansion at both `acs.v2` and
    `acs.v2.dynamics`: `ECMConstitutiveResponseError`,
    `ECMOrientationResponseDiagnostics`,
    `ECMOrientationResponseResult`,
    `step_ecm_orientation_response`,
    `TRACTION_REF_NN_PER_UM2`, `K_ORIENT_PER_S`,
    `TAU_ALIGN_RANGE_S`)
19. **`test_orientation_response_uses_current_schema_not_stale_naming`**
    (§6.6; Y12-guard static check: no `dx_um`, no
    `grid_shape` constructor arg; **both** `ecm.validate()` at
    entry AND `updated_ecm.validate()` before return called per
    Codex `id=1486` two-sided validation contract)
20. **`test_orientation_response_invalid_input_ecm_rejected_before_update`**
    (§3.1; Codex `id=1486`) — passing an `ECMSubstrateState`
    whose `orientation_tensor` violates `validate()` (e.g., a
    component `> 1.0` that the convex update could "heal"
    toward `T_target` under high traction × dt) must raise
    plain schema `ValueError` from `ecm.validate()` at function
    entry, BEFORE any algebra runs. Asserts plain `ValueError`
    (NOT `ECMConstitutiveResponseError`) per Codex `id=1488`:
    schema-level errors stay schema-level, not wrapped.
    Closes the silent-heal path.

If a regression case surfaces during code commit (e.g., Codex
review catches an edge case in the einsum or expm1 call), an
extra test is added with explicit lock reference; the count is
not capped at 20.

---

## 8. Gate verdict

| Item | Status | Notes |
|---|---|---|
| §1 Dimensional | **PASS** | Full Hard Rule 10 unit chain inline-derived; 3 constants Magic-Number Block compliant |
| §2 Boundary | **PASS** | Empty traction / dt=0 / dt<0/bool/nonfinite / strict-positive params / shape mismatch / asymptotic limits all locked |
| §3 Conservation | **PASS** | Schema invariant preserved by convex algebra; unchanged-fields bytewise equal; no-aliasing; distance-to-target monotone per step |
| §4 Numerical | **PASS** | `expm1` for stability near 0; float64; no new tolerance; unconditionally stable in `dt` |
| §5 Sign | **PASS** | Convex weight `∈ [0, 1]`; distance-to-target monotone decreasing per step under fixed target |
| §6 Measurement-protocol | **PASS** | Strictly local algebra (no off-proof points); Rule 10 bridge in diagnostics; Step 6 sister-gate-mirror across 4 layers; Y12-guard meta-test |
| §7 Test catalog | 20 tests planned per locked §4 + Y12-guard + sister-pattern meta-tests + Codex `id=1486` invalid-input rejection test |

**Overall**: gate PASS. impl-work is clear to commit
`acs/v2/dynamics/ecm_constitutive_response.py` + 6 export
updates + `tests/test_v2_ecm_constitutive_response.py` after
Codex review of this Sanity Gate doc, mirroring HB#3 / HB#4 /
Phase D / B1 Sanity Gate review precedent.

### Outstanding before code lands

- Codex review of this Sanity Gate doc per locked §7 6 review
  focus areas:
  1. exact-exponential convex update form + `expm1` numerical
     stability
  2. no aliasing of any of the 5 ECM array fields;
     `validate()` called
  3. typed diagnostics records all input parameters + Rule 10
     bridge fields
  4. `dt_s == 0` valid no-op, others rejected; `k_orient` and
     `traction_ref` strict positive
  5. Y12-guard: `spacing_um` (not `dx_um`), no `grid_shape`
     arg, no mechanosensing shortcut
  6. **7 new symbols** exported through both `acs.v2.dynamics`
     and `acs.v2` (per Codex `id=1488` 6→7 expansion to include
     the locked `ECMConstitutiveResponseError`)
- Code commit must reproduce the locked §1 forbidden list at
  module-docstring level (text-level guard layered on top of
  runtime tests, per B1 `id=1458` precedent).
- All 20 tests must pass at first commit; no `TODO test_X`
  placeholders.
- `pytest tests/test_v2_ecm_constitutive_response.py` +
  combined sister-gate regression
  (`tests/test_v2_fa_to_ecm_scattering.py`,
  `tests/test_v2_ecm_to_fa_bias.py`,
  `tests/test_v2_closed_loop_phase_d.py`,
  `tests/test_v2_focal_adhesion_dynamics.py`,
  `tests/test_v2_protrusion_coupled_focal_adhesion.py`) must
  PASS at commit time.

---

## 9. References

- HB#1+#2 lock:
  `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`
  (commit `7b1d3cf`).
- HB#1 brief (superseded by lock):
  `docs/v2_hard_blocker_1_constitutive_direction_brief.md`
  (commit `4950302`).
- Parent locked plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`.
- Sister-gate Sanity Gate precedents:
  - `docs/v2_fa_to_ecm_scattering_sanity_gate.md` (HB#3).
  - `docs/v2_ecm_to_fa_bias_sanity_gate.md` (HB#4).
  - `docs/v2_phase_d_no_op_scaffolding_sanity_gate.md` (Phase D).
  - `docs/v2_focal_adhesion_dynamics_result_typed_sanity_gate.md` (B1).
- ECM substrate schema (verified for Y12):
  `acs/v2/ecm_substrate.py`.
- ECM open-loop precedent for array-copy-on-return:
  `acs/v2/dynamics/ecm_open_loop.py`.
- Literature (`TRACTION_REF`):
  - Munevar S, Wang Y, Dembo M (2001). "Traction force
    microscopy of migrating normal and H-ras transformed 3T3
    fibroblasts." Biophys J 80(4):1744-1757.
  - Tan JL, Tien J, Pirone DM, Gray DS, Bhadriraju K, Chen CS
    (2003). "Cells lying on a bed of microneedles: an
    approach to isolate mechanical force." PNAS 100(4):1484-1489.
- Literature (`TAU_ALIGN`):
  - Hall MS, Alisafaei F, Ban E, Feng X, Hui CY, Shenoy VB,
    Wu M (2016). "Fibrous nonlinear elasticity enables
    positive mechanical feedback between cells and ECMs." PNAS
    113(49):14043-14048.
  - Notbohm J, Lesman A, Rosakis P, Tirrell DA, Ravichandran G
    (2015). "Microbuckling of fibrin provides a mechanism for
    cell mechanosensing." J R Soc Interface 12(108):20150320.
- Memory rules informing this gate:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, Codex `id=1428`).
  - `rule10_unit_derivation_in_docs.md`.
  - `hard_rule_11_wording_boundary_meta_test.md`.
  - `feedback_aggressive_design_debate.md`.
- Design-discussion thread:
  `id=1466` → `1476` (4-round adversarial lock + seal ack).
- impl-work review thread (this Sanity Gate cycle): `id=1480`
  (design-discussion → impl-work dispatch) → impl-work
  Sanity Gate doc commit (this file) → Codex impl-work review
  → code commit.
