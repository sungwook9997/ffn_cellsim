# V2 Closed-Loop ECM Gate — HB#5 Lyapunov-Like Metric Sanity Gate

**Date**: 2026-05-05 KST
**Status**: pre-execution Sanity Gate for the HB#5 locked design
(commit `340c353`,
`docs/v2_hard_blocker_5_lyapunov_metric_locked.md`). Required by
CLAUDE.md "Sanity Gate Protocol" before the first execution of
any new physics / numerics module. Final Phase E upstream hard
blocker; Step 6 sister-gate-mirror applies against HB#3 / HB#4 /
Phase D / B1 / HB#1+#2 (the typed-dataclass + frozen+slots +
validation-before-act + AST-meta-test sister precedents).
**Source lock**: `docs/v2_hard_blocker_5_lyapunov_metric_locked.md`
(commit `340c353`).
**Source brief** (superseded by lock):
`docs/v2_hard_blocker_5_lyapunov_metric_brief.md` (commit `84b2bcd`).
**Parent locked plan**:
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` (Phase E
entry, HB#5 — final upstream blocker).
**Target module**: `acs/v2/dynamics/ecm_lyapunov_metric.py`
(new; not yet committed).
**Test catalog**: `tests/test_v2_ecm_lyapunov_metric.py` (new;
not yet committed).

This Sanity Gate is the impl-work side's gate before code lands.
HB#5 is a **pure read-only metric**: no ECM mutation, no FA
mutation, no copies needed. The Lyapunov-like contraction
property is a **test-time** assertion (per locked §3 item 3),
NOT an in-function invariant. Phase E activation remains BLOCKED
on the closed-loop ECM gate composition cycle (separate unit).

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

The new module `acs/v2/dynamics/ecm_lyapunov_metric.py`
containing per locked §1:

- 2 module constants (Magic-Number Block compliant per locked §1
  + Y8 derivation chain):
  - `MAX_SQ_FROBENIUS_DIFF_PER_CELL: Final[float] = 9.0`
    (analytical adversarial bound at θ=π/4 with adversarial T;
    schema-tight per §6 derivation).
  - `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL: Final[float] =
    0.5 * MAX_SQ_FROBENIUS_DIFF_PER_CELL` (= 4.5; algebraic
    derivation visible in code).
- 2 typed dataclasses (frozen, slots):
  - `ECMOrientationLyapunovMetricDiagnostics` (10 fields:
    `n_active_cells`, `n_total_cells`, `cell_area_um2`,
    `active_area_um2`, `grid_area_um2`,
    `v_active_max_bound_um2`, `v_active_per_cell_max_um2`,
    `v_active_per_cell_mean_um2`,
    `passive_orientation_magnitude_um2`,
    `max_traction_norm_nN_per_um2`).
  - `ECMOrientationLyapunovMetricResult` (`v_active_um2: float`
    + `diagnostics: ECMOrientationLyapunovMetricDiagnostics`).
- 1 pure function
  `compute_ecm_orientation_lyapunov_metric(ecm,
  traction_density_xy) -> ECMOrientationLyapunovMetricResult`.
- Exports through `acs/v2/dynamics/__init__.py` and
  `acs/v2/__init__.py` (5 new symbols total: 2 constants + 2
  dataclasses + 1 function).

### Explicitly out of scope (will FAIL gate if introduced)

Per locked §1 Forbidden:

- Inclusion of zero-traction cells with `T_target = 0`
  placeholder in the `V_active` sum (Y1 — orientation-magnitude-
  penalty smuggling via placeholder; metric domain MUST match
  HB#1+#2 update domain).
- Reading or referencing `accumulated_traction_nNs_per_um2`,
  `stiffness_kpa`, `fiber_density`, or `ligand_density` (Y7 —
  metric is for orientation-only; AST-walk meta-test enforces
  this).
- Mutation of any ECM array field (function is pure read-only;
  no copies needed because no return-value ECM).
- Pluggable strategy / composite tuple return / dict diagnostics
  (Y5 — single typed result + typed diagnostics, B1 sister-pattern).
- ε / epsilon as physical or numerical regularizer in the active
  mask (Y7 — branch-defined `traction_norm > 0.0` only).
- Hard cap or clamp on `v_active_um2` (the bound is a diagnostic
  field for testing, not an in-line clamp).
- "Positive-definite over full grid" claim (Y2 — actual property
  is positive-semidefinite relative to active-input measurement
  protocol).
- "Cumulative bound = Lyapunov decrease" wording (Y3 —
  boundedness and per-step contraction are separate evidence
  channels and must be tested separately).
- "Full closed-loop FA→ECM→FA stability satisfied by HB#5" claim
  (Y4 — HB#5 is ECM-response-half evidence only; full loop
  requires HB#4-active design, which is currently neutral per
  Phase D no-op).

If any of the above appear in code, the gate FAILs and the unit
halts to surface the scope creep to PI per CLAUDE.md "Sanity
Gate Failure handling".

---

## 1. Dimensional analysis (Hard Rule 10, inline)

The full unit chain is locked at `340c353` §1. Reproduced inline
here per the `rule10_unit_derivation_in_docs` memory rule (every
Rule 10 unit comparison must be inline-derivable in the cited
markdown):

- **Inputs**:
  - `traction_density_xy`: `nN/μm²` (HB#3 output contract).
  - `traction_norm = ||traction_density_xy||₂` (per cell):
    `nN/μm²` componentwise vector norm.
  - `cell_area_um2 = ecm.spacing_um² [μm²]`.
- **Direction (active cells only)**:
  - `n_xy = traction / traction_norm` (unit vector,
    dimensionless) where `traction_norm > 0`.
  - `T_target = n_xy ⊗ n_xy` (rank-1 projection,
    dimensionless components in `[-1, 1]`).
- **Per-cell distance²**:
  - `||T - T_target||²_F = Σ_ij (T_ij - T_target_ij)²`:
    dimensionless (T components dimensionless per schema).
- **V_active (active cells only)**:
  - `V_active_um2 = 0.5 · cell_area_um2 · Σ_active
    ||T - T_target||²_F`
  - Units: `[μm²] · [-] = [μm²]` ✓ (the energy is reported in
    `μm²`, the natural area unit for a per-cell dimensionless
    integrand on the ECM grid).
- **Bound**:
  - `MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0` (dimensionless).
  - `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 0.5 · 9.0 =
    4.5` (dimensionless).
  - `v_active_max_bound_um2 = 4.5 · active_area_um2` (units
    `[-] · [μm²] = [μm²]` ✓).

**Magic-Number Block 3-test verdict** (zero new tunables beyond
the 2 locked algebraic constants):

1. **Derivable** ✓:
   - `MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0` derived analytically
     from schema invariant `|T_ij| ≤ 1` + symmetry; worst-case
     `T = ((-1, -1), (-1, -1))` against unit-vector target at
     `θ = π/4` (see §6 below). Pure mathematical bound, no
     literature parameter.
   - `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 4.5` is
     algebraically `0.5 · 9.0` (the energy-form factor); the
     `0.5` comes from the standard quadratic-form definition
     `V = 0.5 · ||...||²`.
2. **Grid-invariant** ✓: no `dx` / `dy` / `dt` / `grid_n`
   dependence in either constant. The bound `v_active_max_bound`
   is `4.5 · active_area_um2` — proportional to actual measured
   area, which scales with grid + spacing changes.
3. **Not fitted** ✓: chosen analytically from schema invariants
   before any A/A₀ comparison or test target gating.

### Status

**PASS**. Full unit chain inline-derived per Rule 10. 2 constants
Magic-Number Block compliant. No new tolerance, threshold, or
per-call numerical regularizer. The 9 → 4.5 algebraic chain
makes derivation traceability visible in the code (Y8 lock).

---

## 2. Boundary cases

Per locked §3 item 2, plus the validation contracts locked at
Y7 + Y10:

### 2.1 All-zero traction (`active.any() == False`)

Per the lock pseudocode:
- `traction_norm` all 0 → `active` mask all False.
- `n_xy` left as `np.zeros_like` zeros (no division).
- `T_target = einsum(zeros, zeros) = zeros`.
- `sq_diff_per_cell = ||T - 0||²_F = ||T||²_F` per cell.
- The `if active.any():` branch is False → `v_active = 0.0`,
  `v_per_cell_max = 0.0`, `v_per_cell_mean = 0.0`.
- The `if passive.any():` branch is True (all cells passive) →
  `passive_orientation_magnitude = 0.5 · cell_area · ||T||²_F`
  (full-grid sum).

Tested by test 1 (`test_lyapunov_metric_zero_traction_returns_zero`)
+ test 2
(`test_lyapunov_metric_no_active_cells_no_passive_penalty_in_v`).

### 2.2 All-active traction (`passive.any() == False`)

`passive_orientation_magnitude = 0.0` (no passive cells to
contribute). `V_active` computed normally over all cells. Tested
implicitly via test 6
(`test_lyapunov_metric_varying_target_boundedness`) which uses
all-active traction.

### 2.3 At-target orientation (`T == n⊗n` per active cell)

`sq_diff_per_cell = 0` per active cell → `v_active = 0.0` exactly.
Per-cell max = 0, per-cell mean = 0.

Tested by test 3 (`test_lyapunov_metric_at_target_returns_zero`).

### 2.4 Adversarial schema-saturation

`T = ((-1, -1), (-1, -1))` everywhere (componentwise schema-
extreme), traction along `(1/√2, 1/√2)` (`θ = π/4` unit
direction):
- Per-cell `||T - T_target||²_F = 9.0` exactly (verified by §6
  derivation).
- `v_active_um2 = 0.5 · cell_area · n_active · 9.0 = 4.5 ·
  active_area_um2` (saturates the locked bound).

Tested by test 4
(`test_lyapunov_metric_v_active_within_bound_4_5_times_active_area`)
with **dual assertion** per Y10:
- per-cell raw squared Frobenius diff = 9.0
- `v_active_um2 ≈ 4.5 · active_area_um2`

### 2.5 Single-cell ECM (`grid_shape = (1, 1)`)

All counts + bounds well-defined (no special-case branching in
the pseudocode). Sanity check covered implicitly by general
test cases; no separate test.

### 2.6 Invalid input ECM

`ecm.validate()` runs FIRST at function entry per HB#1+#2
sister-precedent (`id=1486` silent-heal-path closure). Invalid
input ECM raises plain `ValueError` from schema; NOT wrapped in
any HB#5-specific error type (no error type defined; no
wrapping).

Tested by test 10 (`test_lyapunov_metric_invalid_ecm_raises_at_validate`).

### 2.7 Non-finite / shape-mismatch traction

- `traction_density_xy` containing `nan` / `±inf`: raises (test
  8 `test_lyapunov_metric_traction_nonfinite_raises`).
- `traction_density_xy.shape != (*ecm.grid_shape, 2)`: raises
  (test 9 `test_lyapunov_metric_traction_shape_mismatch_raises`).

The lock pseudocode shows `_validate_finite_traction` +
`_validate_traction_shape` calls; impl-work decides whether to
import these from `ecm_constitutive_response.py` or duplicate
minimally (per locked §5). **My lean**: duplicate minimally with
plain `ValueError` raises (no HB#5-specific error type per the
lock — no error type defined). This keeps the modules
independent and avoids cross-dependency on validators that may
evolve in the constitutive-response module.

### Status

**PASS**. All 7 boundary classes (zero/all-active/at-target/
adversarial-saturation/single-cell/invalid-ECM/non-finite-or-
shape) explicitly locked + 6 of 7 covered by tests (single-cell
covered implicitly). No NaN paths.

---

## 3. Conservation invariants

### 3.1 Pure read-only

The function performs **no mutation** of `ecm` or any of its
fields. No copies are needed because the function does NOT
return an ECM (only a typed result `ECMOrientationLyapunovMetricResult`
with `v_active_um2` + diagnostics). This is a direct corollary
of "pure read-only metric" per locked §1 + §3 item 3.

Tested by test 7 (`test_lyapunov_metric_pure_no_mutation_of_ecm`)
which checks all 5 ECM array fields bytewise unchanged after the
metric call.

### 3.2 Fixed-target per-step contraction (Y3 evidence channel A)

Per locked §6: the contraction property is **per-step under
fixed target**:

```
T_new − T_target = (1−w)(T_old − T_target)
‖T_new − T_target‖_F ≤ ‖T_old − T_target‖_F  for w ∈ [0, 1]
```

This is a **test-time** assertion (per locked §3 item 3), NOT
an in-function invariant. The metric itself only computes the
current-state distance; the contraction is observable by
applying one HB#1+#2 update and computing the metric before +
after.

Tested by test 5
(`test_lyapunov_metric_fixed_target_contraction_under_hb12_step`).

### 3.3 Varying-target boundedness (Y3 evidence channel B)

Per locked §1 + §3 item 5: `0 ≤ V_active ≤
v_active_max_bound_um2` at every measured step under any
varying-target sequence. The bound is exact (`4.5 ·
active_area_um2`) because per-cell `||T - T_target||²_F ≤ 9` by
the §6 derivation.

Tested by test 6
(`test_lyapunov_metric_varying_target_boundedness`) which
applies many steps with rotating traction direction and asserts
the bound at every measured step.

### 3.4 Boundedness + contraction NOT conflated (Y3 hard boundary)

The lock §1 forbidden list explicitly rejects "cumulative bound
= Lyapunov decrease" wording. The two evidence channels are
**separate** (test 5 vs test 6). No claim that V is monotone
under varying target; no claim that the bound implies
contraction.

### Status

**PASS**. Pure read-only — nothing to conserve in the function
itself. Two evidence channels locked + tested separately + the
"don't conflate" forbidden item enforces the discipline at the
lock-text level. Test 7 (no-mutation) is the runtime guard.

---

## 4. Numerical sanity

### 4.1 Float precision

`np.float64` enforced everywhere (matches HB#1+#2 + ECM schema).
Validators coerce input `traction_density_xy` to `float64` if
needed.

### 4.2 No new tolerance constants

The 2 locked constants (9.0 + 4.5) are dimensionless algebraic
bounds, NOT tolerances. Zero numerical-tolerance / round-off
constants introduced (the lock explicitly drops ε from the
brief — Y7).

### 4.3 Per-call work

Per-call work = 1 `linalg.norm` (O(grid_size)) + 1 `einsum`
outer product (O(grid_size · 2 · 2 = 4·grid_size)) + 1 `T -
T_target` subtraction + 1 `**2` square + 1 `sum(axis=(-2,-1))`
+ active-mask reductions. Linear in grid_size; no asymptotic
surprise.

### 4.4 `einsum` numerical stability

`np.einsum("...i,...j->...ij", n_xy, n_xy)` is exact algebra (2x2
outer product per cell); no accumulation order issues.

### 4.5 Frobenius² stability

`np.sum(diff**2, axis=(-2, -1))` is componentwise positive; sum
order is stable for typical grid sizes (no Kahan summation
required at this scale per locked §3 item 4).

### 4.6 `cell_area_um2` exact float square

`cell_area_um2 = float(ecm.spacing_um) ** 2` — exact float
squaring, no rounding artifact.

### Status

**PASS**. Float64 throughout; no new tolerance; per-call work is
linear in grid_size; einsum + Frobenius² numerically stable for
the 2x2-tensor + grid-size regime.

---

## 5. Sign / sense check

Per locked §3 item 5:

- `V_active ≥ 0` always (sum of squares times nonneg `0.5 ·
  cell_area`).
- Per-cell components in `sq_diff_per_cell ≥ 0` (sum of
  squares).
- Bound `v_active_max_bound_um2 ≥ 0` always (`4.5 · area ≥ 0`).
- Passive diagnostic `passive_orientation_magnitude_um2 ≥ 0`
  (sum of squares times nonneg factor).
- All 10 diagnostics fields are non-negative quantities by
  construction (counts ≥ 0, areas ≥ 0, energies ≥ 0, max ≥
  mean ≥ 0).

### Status

**PASS**. No sign ambiguity; all output fields non-negative by
construction; sign-of-difference (`T - T_target`) absorbed by
the squaring.

---

## 6. Measurement-protocol consistency (Hard Rule 11) — central anchor

The locked §6 calls this "the central anchor" because the v13
Stage 1a `_build_curvature` Laplacian off-peak episode is the
canonical anti-pattern: an analytical proof valid in an active
region was extended past the active region by the band-averaged
measurement, yielding a 5× worse result.

### 6.1 Domain matching (Y1)

The HB#1+#2 update domain is **active cells only**: zero-traction
cells have `w = 0`, target placeholder algebraically unused.
HB#5's `V_active` MUST also be restricted to active cells.

Including zero-traction cells with `T_target = 0` placeholder in
the V sum would smuggle an orientation-magnitude penalty into
the Lyapunov energy — exactly the v13 Stage 1a anti-pattern.

The lock's `if active.any():` branch enforces this: the V sum
runs only over `sq_diff_per_cell[active]`. Passive-cell
orientation magnitude is reported as a **separate diagnostic
field** `passive_orientation_magnitude_um2`, NOT part of `V`.

Tested by test 12
(`test_lyapunov_metric_passive_diagnostic_separate_from_v_active`)
which sets T to large magnitude in passive cells and verifies
`V_active` is unchanged while the passive diagnostic reflects
the magnitude.

### 6.2 Adversarial bound derivation (the §1 constant 9.0) — inline

For unit `n = (cos θ, sin θ)`:

```
||T - T_target||²_F
  = (T_xx − cos²θ)² + 2·(T_xy − cos θ · sin θ)² + (T_yy − sin²θ)²
```

Schema invariant: `|T_ij| ≤ 1` componentwise; T symmetric. No
PSD assumption. Adversarial T (componentwise extreme against
target):

- `T_xx = -1`, `T_yy = -1`, `T_xy = -sign(sin 2θ)`.

Substitute at `θ = π/4` (`cos²θ = sin²θ = 0.5`,
`cos θ · sin θ = 0.5`, `sign(sin 2θ) = +1`):

```
||T - T_target||²_F
  = (-1 − 0.5)² + 2·(-1 − 0.5)² + (-1 − 0.5)²
  = (1.5)² · 4
  = 9.0
```

**Verification via Frobenius identity** (per locked §6):

```
||T - n⊗n||²_F = ||T||²_F + ||n⊗n||²_F − 2 · trace(T · n⊗n)
              = ||T||²_F + 1 − 2 · n^T T n
```

For the adversarial setup: `||T||²_F = 4` (4 entries each at
`±1`); `||n⊗n||²_F = (cos²θ)² + 2·(cos θ sin θ)² + (sin²θ)² =
0.25 + 0.5 + 0.25 = 1`; `n^T T n` at `θ = π/4` with `T =
((-1, -1), (-1, -1))`: `(1/√2 −1/√2)·(−1, −1)·(1/√2, 1/√2) =
−2`. So:

```
||T - n⊗n||²_F = 4 + 1 − 2·(−2) = 4 + 1 + 4 = 9.0 ✓
```

Therefore `MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0` is
**schema-tight** and **derivable from first principles**.

### 6.3 Two-constant traceability (Y8)

The lock uses 2 constants (`9.0` + `4.5 = 0.5 · 9.0`) instead of
just `4.5` to make the `9 → 4.5` algebraic derivation visible in
the code itself. Test 4 dual-asserts both: per-cell raw squared
Frobenius diff = 9.0 AND `v_active_um2 ≈ 4.5 · active_area_um2`.
This guards against the 9-vs-4.5 naming mixup reappearing later
(Y8 was triggered by Claude round-2 incorrectly naming
`MAX_SQ_FROBENIUS_DIFF_PER_CELL = 4.5`).

### 6.4 Two-channel evidence separation (Y3)

Per §3.4: contraction (test 5) and boundedness (test 6) are
**separate** evidence channels. The lock §1 forbidden list
explicitly rejects "cumulative bound = Lyapunov decrease"
wording.

### 6.5 Closed-loop coupling depth (Y4)

HB#5 is **ECM-response-half evidence only**. The locked §0 +
§1 forbidden list explicitly preclude "full FA→ECM→FA stability
satisfied by HB#5" claims. Phase E composition can use HB#5 as
one input to a bounded-feedback evaluation, but full closed-loop
stability requires HB#4-active design (currently HB#4 = Phase D
no-op neutral multipliers).

### 6.6 AST-walk meta-test scope (Y11)

Test 14 is an AST walk on
`inspect.getsource(compute_ecm_orientation_lyapunov_metric)`
wrapped with `textwrap.dedent(...)` for robustness against
decorators / formatting changes.

The walk asserts that none of `accumulated_traction_nNs_per_um2`,
`stiffness_kpa`, `fiber_density`, `ligand_density` appear as
`ast.Attribute.attr` or `ast.Name.id` in the function body
(NOT module-level docstrings or comments — Y11 explicitly
scoped the walk to the function body to avoid false failures
from forbidden-field names appearing in docstrings).

The lock's pseudocode body itself does not embed forbidden
field names. The lock-document (this Sanity Gate too) does
mention them in §0 forbidden lists, but those are markdown
text, not Python source — they don't affect the AST walk.

### 6.7 Sister-gate-mirror application (Step 6 of 6-step batch)

Per memory `design_note_pre_commit_batch.md` 6-step, Step 6
sister-gate-mirror layered against HB#3 / HB#4 / Phase D / B1 /
HB#1+#2:

- **Code (validation ordering)**: HB#5 mirrors HB#1+#2's
  `ecm.validate()` first (id=1486 silent-heal-path closure
  precedent) → finite traction → traction shape. Same
  validate-before-act discipline.
- **Design (typed dataclass + frozen+slots)**: HB#5 follows
  HB#1+#2 / B1 / Phase D `@dataclass(frozen=True, slots=True)`
  precedent for `ECMOrientationLyapunovMetricDiagnostics` +
  `ECMOrientationLyapunovMetricResult`.
- **API surface (`__init__.py` exports)**: 5 new symbols at
  both `acs/v2/__init__.py` and
  `acs/v2/dynamics/__init__.py` with `__all__` inclusion. Test
  13 enforces.
- **Failure-kind discipline**: HB#5 has NO HB#5-specific error
  class (per lock §5 + the "no new error class" implicit
  decision; locked validators raise plain `ValueError`).
  Different from HB#1+#2 which has `ECMConstitutiveResponseError`
  with 5 owned `failure_kind`s. Rationale per locked §5: HB#5
  is a pure read-only metric with no domain-specific failure
  modes; the failure modes (non-finite traction, shape
  mismatch, invalid ECM) are all schema/input issues that are
  best surfaced as plain `ValueError`. Codex `id=1488`-style
  failure_kind discipline does NOT apply because the metric
  does not introduce any constitutive-law-level error.

### Status

**PASS**. Domain matching (active cells only); adversarial bound
9.0 derived analytically + verified via Frobenius identity;
two-constant traceability (9 → 4.5) visible in code; two
evidence channels separated (Y3); closed-loop coupling depth
restricted to ECM half (Y4); AST-walk meta-test scoped to
function body (Y11); Step 6 sister-gate-mirror applied at all 4
discipline layers.

---

## 7. Magic-Number Block check

**Zero new tunables introduced beyond the 2 locked algebraic
constants.**

The module declares:
- `MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0` (analytical adversarial
  bound; §1 + §6 derived).
- `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 0.5 *
  MAX_SQ_FROBENIUS_DIFF_PER_CELL` (= 4.5; algebraic chain).

Magic-Number Block 3-test verdict (per locked §1 +
`design_note_pre_commit_batch.md`):

1. **Derivable** ✓: `9.0` from schema invariant + adversarial
   maximization (§6); `4.5` from `0.5 · 9.0` (energy-form
   factor).
2. **Grid-invariant** ✓: no `dx` / `dy` / `dt` / `grid_n`
   dependence. The bound `v_active_max_bound = 4.5 · area`
   scales with measured area, not grid resolution.
3. **Not fitted** ✓: chosen analytically before any A/A₀
   comparison or test target gating.

### Status

**PASS**. Zero non-derived tunables. Both constants
algebra-traceable; the `9 → 4.5` chain is visible in the code
itself.

---

## 8. Test catalog (~15 tests per locked §4 after 14/15 merge)

Owned by `tests/test_v2_ecm_lyapunov_metric.py` (not yet
committed). Each test maps to a locked invariant in
`docs/v2_hard_blocker_5_lyapunov_metric_locked.md` §4.

1. `test_lyapunov_metric_zero_traction_returns_zero` (§2.1) —
   `V_active = 0.0`; passive diagnostic = full-grid orientation
   magnitude.
2. `test_lyapunov_metric_no_active_cells_no_passive_penalty_in_v`
   (§2.1; Y1 forward guard) — V_active does not include
   passive cells; passive diagnostic does.
3. `test_lyapunov_metric_at_target_returns_zero` (§2.3) —
   `T == n⊗n` per active cell → `V_active = 0.0` exactly.
4. `test_lyapunov_metric_v_active_within_bound_4_5_times_active_area`
   (§2.4; Y8 + Y10 dual assertion) — adversarial setup; per-cell
   raw squared Frobenius diff = 9.0 AND `v_active_um2 ≈
   4.5 · active_area_um2`.
5. `test_lyapunov_metric_fixed_target_contraction_under_hb12_step`
   (§3.2; Y3 channel A) — apply one
   `step_ecm_orientation_response` with same traction; assert
   `V_active_after ≤ V_active_before`.
6. `test_lyapunov_metric_varying_target_boundedness` (§3.3; Y3
   channel B) — many steps with rotating traction; assert
   `0 ≤ V_active ≤ V_active_max_bound_um2` at every step.
7. `test_lyapunov_metric_pure_no_mutation_of_ecm` (§3.1) — all 5
   ECM array fields bytewise unchanged after metric call.
8. `test_lyapunov_metric_traction_nonfinite_raises` (§2.7).
9. `test_lyapunov_metric_traction_shape_mismatch_raises` (§2.7).
10. `test_lyapunov_metric_invalid_ecm_raises_at_validate` (§2.6;
    HB#1+#2 `id=1486` sister-pattern) — schema-invalid input ECM
    raises plain `ValueError` from `ecm.validate()` BEFORE any
    algebra.
11. `test_lyapunov_metric_diagnostics_typed_dataclass` (B1
    sister-pattern meta-test) — `isinstance(...)` + `with
    pytest.raises(TypeError):
    result.diagnostics["v_active_max_bound_um2"]`.
12. `test_lyapunov_metric_passive_diagnostic_separate_from_v_active`
    (§6.1; Y1 forward guard) — set T to large in passive cells;
    `V_active` unchanged; `passive_orientation_magnitude_um2`
    reflects the magnitude.
13. `test_lyapunov_metric_exports_through_both_init` (§6.7
    sister-gate-mirror API surface — 5 new symbols at both
    `acs.v2` and `acs.v2.dynamics`).
14. `test_lyapunov_metric_does_not_reference_other_ecm_fields`
    (§6.6; Y11 + Codex `id=1506` review focus 5) — AST walk on
    `textwrap.dedent(inspect.getsource(...))` of the function
    body; assert none of `accumulated_traction_nNs_per_um2`,
    `stiffness_kpa`, `fiber_density`, `ligand_density` appear
    as `ast.Attribute.attr` or `ast.Name.id`.
15. (Test 14 absorbed all 4 forbidden names; tests 14/15 merged
    per locked §4.)
16. `test_lyapunov_metric_active_area_consistent_with_n_active`
    — `active_area_um2 == n_active_cells · cell_area_um2`
    (sanity).
17. `test_lyapunov_metric_max_bound_consistent_with_active_area`
    (§6.3; Y8 derivation chain check) —
    `v_active_max_bound_um2 == 4.5 · active_area_um2`.

(15 active tests; tests 14/15 merged per locked §4.)

If a regression case surfaces during code commit (e.g., Codex
review catches an edge case in einsum or active-mask handling),
an extra test is added with explicit lock reference; the count
is not capped at 15.

---

## 9. Gate verdict

| Item | Status | Notes |
|---|---|---|
| §1 Dimensional | **PASS** | Full Hard Rule 10 unit chain inline-derived; 2 constants Magic-Number Block compliant; 9 → 4.5 algebraic chain visible |
| §2 Boundary | **PASS** | 7 boundary classes locked + 6 of 7 test-covered (single-cell implicit) |
| §3 Conservation | **PASS** | Pure read-only; two evidence channels (contraction + boundedness) tested separately per Y3 |
| §4 Numerical | **PASS** | float64; einsum + Frobenius² stable; no new tolerance |
| §5 Sign | **PASS** | All output fields non-negative by construction |
| §6 Measurement-protocol | **PASS** | Domain matching (active cells only) is the central anchor; adversarial bound 9.0 derived analytically; AST-walk meta-test scoped to function body; Step 6 sister-gate-mirror applied at all 4 layers |
| §7 Magic-Number Block | **PASS** | Zero non-derived tunables; both constants algebra-traceable |

**Overall**: gate PASS. impl-work is clear to commit
`acs/v2/dynamics/ecm_lyapunov_metric.py` + 5 export updates +
`tests/test_v2_ecm_lyapunov_metric.py` after Codex review of
this Sanity Gate doc, mirroring HB#3 / HB#4 / Phase D / B1 /
HB#1+#2 Sanity Gate review precedent.

### Outstanding before code lands

- Codex review of this Sanity Gate doc per locked §7 6 review
  focus areas:
  1. active-cells-only V (Hard Rule 11) — passive cells in
     separate diagnostic, NOT in V
  2. two-constant form `MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0` +
     `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 4.5`
     (derivation traceability)
  3. branch-defined active mask (no epsilon)
  4. `ecm.validate()` at entry; pure read-only (no mutation, no
     copies)
  5. AST-walk meta-test scoped to function body via
     `textwrap.dedent(inspect.getsource(...))`
  6. 5 new symbols exported through both `acs.v2.dynamics` and
     `acs.v2`
- Code commit must reproduce the locked §1 forbidden list at
  module-docstring level (text-level guard layered on top of
  runtime tests, per B1 `id=1458` precedent).
- All 15 tests must pass at first commit; no `TODO test_X`
  placeholders.
- `pytest tests/test_v2_ecm_lyapunov_metric.py` + combined
  sister-gate regression (HB#3 / HB#4 / Phase D / 6.3a/6.3b /
  HB#1+#2 + HB#5) must PASS at commit time.

---

## 10. References

- HB#5 lock:
  `docs/v2_hard_blocker_5_lyapunov_metric_locked.md` (commit
  `340c353`).
- HB#5 brief (superseded by lock):
  `docs/v2_hard_blocker_5_lyapunov_metric_brief.md` (commit
  `84b2bcd`).
- Parent locked plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` (Phase E
  entry, HB#5 — final upstream blocker).
- Sister-gate Sanity Gate precedents:
  - `docs/v2_fa_to_ecm_scattering_sanity_gate.md` (HB#3)
  - `docs/v2_ecm_to_fa_bias_sanity_gate.md` (HB#4)
  - `docs/v2_phase_d_no_op_scaffolding_sanity_gate.md` (Phase D)
  - `docs/v2_focal_adhesion_dynamics_result_typed_sanity_gate.md`
    (B1)
  - `docs/v2_hard_blocker_1_2_constitutive_response_sanity_gate.md`
    (HB#1+#2)
- ECM substrate schema: `acs/v2/ecm_substrate.py`
  (`spacing_um`, componentwise `|T_ij| ≤ 1`, no PSD).
- HB#1+#2 implementation (sister precedent for validation
  pattern + sister module under `acs/v2/dynamics/`):
  `acs/v2/dynamics/ecm_constitutive_response.py`.
- v2 measurement registry (separate concern, NOT this unit's
  home): `acs/v2/metrics.py` (file, not directory).
- CLAUDE.md Sanity Gate §6 (measurement-protocol consistency,
  v13 Stage 1a `_build_curvature` off-peak anti-pattern).
- Memory rules informing this gate:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, per Codex `id=1428`).
  - `rule10_unit_derivation_in_docs.md`.
  - `hard_rule_11_wording_boundary_meta_test.md`.
  - `feedback_aggressive_design_debate.md`.
- Design-discussion thread: `id=1498` → `1506` (3-round
  adversarial lock + seal ack); `id=1509` reinforcement.
- impl-work review thread (this Sanity Gate cycle): `id=1512`
  (design-discussion → impl-work dispatch) → impl-work Sanity
  Gate doc commit (this file) → Codex impl-work review → code
  commit.
