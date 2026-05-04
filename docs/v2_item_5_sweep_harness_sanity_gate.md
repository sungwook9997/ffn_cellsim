# V2 Closed-Loop ECM Gate — Item 5 Sweep Harness Sanity Gate

**Date**: 2026-05-05 KST
**Status**: pre-execution Sanity Gate for the Item 5 sweep
harness locked design (commit `c5634d5`,
`docs/v2_item_5_sweep_harness_locked.md`). Required by CLAUDE.md
"Sanity Gate Protocol" before the first execution of any new
physics / numerics / harness module. Step 6 sister-gate-mirror
layered against Phase E v1 / HB#5 / HB#1+#2 / HB#3 / HB#4 /
Phase D / B1 (the typed-dataclass + frozen+slots + evidence-only
wording + AST-meta-test sister precedents).
**Source lock**: `docs/v2_item_5_sweep_harness_locked.md` (commit
`c5634d5`).
**Source brief** (superseded by lock):
`docs/v2_item_5_sweep_harness_brief.md` (commit `4f0efdd`).
**Parent locked plan**:
`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` (Item 5
sweep is a separate harness from Phase E v1, per parent §1).
**Target module**:
`acs/v2/dynamics/closed_loop_phase_e_sweep.py` (new; not yet
committed).
**Test catalog**:
`tests/test_v2_closed_loop_phase_e_sweep.py` (new; not yet
committed).

This Sanity Gate is the impl-work side's gate before code lands.
Item 5 sweep harness is **Phase E v1 sensitivity evidence ONLY,
NOT closed-loop ECM gate Item 5 satisfaction** (locked §0
wording boundary inheriting Phase E v1 Y2 sister-pattern); raw
observable summaries are produced; satisfaction decision is a
separate caller-side / PI / explicit-decision unit.

This is the **third Phase E composition Hard Rule 11
measurement-protocol catch family** (per locked §6 references):
- Phase E v1 Y1: ECM-side evidence wording boundary
- HB#5 Y1: active-cells-only V_active
- **Item 5 Y1 + Y4 + Y17 (this unit)**: physical-domain
  invariance under spacing sweep + ECM-feedback substep wording
  + Step 0 semantics

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

The new module `acs/v2/dynamics/closed_loop_phase_e_sweep.py`
containing per locked §1:

- 3 typed dataclasses (frozen+slots):
  - `PhaseEV1Item5SweepConfig` (input contract; Y1 physical-
    domain invariance via fixed `domain_size_um_xy`, derived
    `nx = domain / spacing`; Y20 Literal/n_revolutions/magnitude
    runtime validation; Y21 domain/origin/FA-position
    defense-in-depth)
  - `PhaseEV1Item5SweepMetadata` (Y13 NO wall-clock timestamp;
    Y14 git_sha override + fallback)
  - `PhaseEV1Item5SweepResult` (Y6 trajectory shapes
    `(n_runs, n_steps+1)` including step 0; Y8 raw observable
    summaries — NO satisfaction/monotone claims)
- 1 pure deterministic library function
  `run_phase_e_v1_sensitivity_sweep(config) -> PhaseEV1Item5SweepResult`
- 7 private helpers (NOT exported per Y12):
  `_validate_config`, `_build_rotating_uniform_single_fa`,
  `_build_initial_ecm_at_spacing`, `_detect_git_sha`,
  `_compute_max_bound_ratio`, `_record_step`,
  `_default_phase_e_v1_item5_sweep_config`
- Exports through `acs/v2/dynamics/__init__.py` and
  `acs/v2/__init__.py` (only 4 new symbols per Y12: 3
  dataclasses + 1 function).

### Explicitly out of scope (will FAIL gate if introduced)

Per locked §1 Forbidden:

- "Item 5 satisfied" / "satisfies_item_5_*" wording in module
  docstring, function docstring, or test names (Y8 +
  sister-pattern Phase E v1 Y2 — meta-test enforces
  evidence-only naming).
- Reuse of ECM array across spacings without resampling on the
  new grid (Y1 — would conflate resolution with IC change).
- FA position varying across spacings (Y1 — physical position
  must be invariant).
- `domain_size_um_xy` not exactly divisible by some
  `spacing_um` value (Y1 — validation raises).
- Empty `spacing_um_values` or `dt_s_values` reduced to a
  single call (Y3 — must raise).
- `dt_s == 0` inside the sweep grid (Y3 — boundary case for
  Phase E v1, not part of Item 5 default sweep).
- Wall-clock timestamp in metadata (Y13 — breaks
  reproducibility).
- Untyped `dict` for primary metadata (Y7 + B1 sister-pattern).
- Local definition of `K_ORIENT_PER_S` or
  `TRACTION_REF_NN_PER_UM2` (Y14 — must import from
  `ecm_constitutive_response`, sister-pattern with Phase E v1
  Y14).
- "varies across axis" qualitative sensitivity claim (Y19 —
  symmetry may produce equal values; assert structural
  representation + analytical formula consequences instead).
- Gate tolerance threshold tuning (Y18 — `1e-12` is IEEE
  roundoff allowance, NOT a tunable threshold; HB#5 invariant
  is mathematical `≤ 1.0`).
- "monotone in spacing" claim (Y8 — sensitivity not guaranteed
  monotone).
- Public export of `_default_phase_e_v1_item5_sweep_config()`
  helper (Y12 — private only).
- Wall-time runtime assertion in tests (Y15 — informational
  only).

If any of the above appear in code, the gate FAILs and the unit
halts to surface the scope creep to PI per CLAUDE.md "Sanity
Gate Failure handling".

---

## 1. Dimensional analysis (Hard Rule 10, by inheritance)

The Item 5 sweep harness performs **no unit transformation**;
all sub-call unit chains are preserved by the cross-product
loop:

- HB#3 scatter: `traction_density_xy [nN/μm²]` per Phase E v1
  + HB#3 lock chains.
- HB#1+#2 update: `dt_s [s]`, `spacing_um [μm]`,
  `traction_norm [nN/μm²]`, dimensionless `S` and `w` per
  HB#1+#2 lock chain.
- HB#5 metric: `V_active_um2 [μm²]`, `V_active_bound_um2 [μm²]`
  per HB#5 lock chain.
- HB#4 neutral bias: dimensionless multipliers.

Sweep-level outputs:

- `spacing_um_by_run [μm]`, `dt_s_by_run [s]`, `nx_by_run` /
  `ny_by_run` (dimensionless cell counts).
- `max_bound_ratio = V_active / V_bound` dimensionless ∈ `[0, 1
  + 1e-12]` (Y18 IEEE roundoff allowance, NOT gate tunable).
- `peak_v_active_um2`, `final_v_active_um2 [μm²]` inherited
  from HB#5.

**Magic-Number Block 3-test verdict** (per locked §1 + Y11 +
Y14):

1. **Derivable** ✓:
   - `0.5 * I` orientation IC (Y11): midpoint between zero
     orientation and unit diagonal schema bound; not fitted;
     symmetric non-target IC giving the metric room to evolve
     toward any traction direction.
   - `1e-12` (Y18): IEEE float64 roundoff allowance; standard
     numerical analysis margin; NOT a tunable gate threshold.
   - All physics constants (`K_ORIENT_PER_S`,
     `TRACTION_REF_NN_PER_UM2`) re-imported from HB#1+#2 lock
     (Y14 — object-identity verified).
2. **Grid-invariant** ✓: harness explicitly varies
   `spacing_um`; `0.5*I` IC + `1e-12` roundoff are
   spacing-independent quantities.
3. **Not fitted** ✓: chosen analytically per Y11 + Y18
   reasoning before any A/A₀ comparison or test target gating.

### Status

**PASS**. No new dimensional chain at harness layer; sub-call
chains preserved by loop. Magic-Number Block compliant for
both `0.5*I` IC (Y11) and `1e-12` roundoff (Y18). All physics
constants re-imported via Y14 object-identity-verified pattern.

---

## 2. Boundary cases

Per locked §3 item 2:

### 2.1 Empty axis arrays (Y3)

`config.spacing_um_values = ()` or `config.dt_s_values = ()` →
`_validate_config` raises. Tests 1 + 2 cover.

### 2.2 Non-positive / nonfinite values (Y3)

`spacing_um < 0`, `dt_s ≤ 0`, `traction_magnitude_nN ≤ 0`,
`nan`/`±inf` → raise. Test 3 covers.

### 2.3 Indivisible domain (Y1)

`domain_size_um_xy = (16.0, 16.0)` requires every `spacing_um`
to divide exactly: `16.0 / 0.5 = 32.0` ✓; `16.0 / 1.5 = 10.67`
→ raise with explicit error showing bad spacing + non-integer
nx. Test 4 covers.

### 2.4 Unsupported traction scenario (Y20)

`config.traction_scenario != "rotating_uniform_single_fa"` →
raise (Literal at type level doesn't prevent runtime bad
values; runtime validate explicit). Test 5 covers.

### 2.5 `n_revolutions` validation (Y20)

- `n_revolutions == 0` → raise (constant-target ≠ locked
  rotating scenario). Test 6.
- `n_revolutions < 0` → raise. Test 7.
- `n_revolutions` of type `bool` → raise (Python `bool ⊂ int`
  trap, HB#1+#2 sister). Test 8.

### 2.6 Domain / origin / FA-position validation (Y21)

- `domain_size_um_xy` not length-2 / not finite / not positive
  → raise. Test 9.
- `origin_um_xy` not length-2 / not finite → raise. Test 10.
- `fa_position_um_xy` not length-2 / not finite / outside
  inclusive physical domain `[ox, ox+dx] × [oy, oy+dy]` →
  raise. Test 11.

### 2.7 Step 0 semantics (Y17)

Step 0 = initial ECM measured against step-0 traction target
with `dt_s=0`. NOT a zero-traction baseline.

- `v_active_um2[run, 0]` may be nonzero (initial misalignment
  measurable; for `0.5*I` IC + rotating-traction scenario at
  θ=0, distance to target `n⊗n = ((1,0),(0,0))` is nonzero).
- `v_active_bound_um2[run, 0]` may be nonzero (active cells
  exist).
- `max_traction_norm[run, 0]` may be nonzero (step-0 traction
  scenario active).
- `max_convex_weight[run, 0] == 0` exact (`dt_s=0` →
  `w = -expm1(0) = 0`).
- `max_orientation_delta[run, 0] == 0` exact (no update at
  `dt_s=0`).

Tested by test 21
(`test_item_5_sweep_step_0_initial_against_target_not_zero_baseline`).

### 2.8 All-zero traction (degenerate scenario, Y16)

`max_bound_ratio = 0.0` (no division by zero); the
`_compute_max_bound_ratio` helper returns 0.0 if no
active-bound-positive samples. Tested by test 19.

### Status

**PASS**. 8 boundary classes locked + 11 + 3 = 14 tests covering
validation + Step 0 semantics + zero-bound case. No NaN paths
introduced.

---

## 3. Conservation invariants

### 3.1 Harness has no per-step conservation of its own

The harness wraps Phase E v1 in a cross-product loop. Each
sub-call enforces its lock's conservation:

- HB#3 scatter conservation (per HB#3 lock).
- HB#1+#2 schema invariant preservation (`updated_ecm.validate()`
  at sub-call return).
- HB#4 neutral multiplier structural invariant.
- HB#5 V_active ≥ 0 + bounded by `4.5 · active_area_um2`.

The harness adds NO new conservation invariant; it only
**measures** the sweep.

### 3.2 Trajectory shape invariant

`v_active_um2_traj.shape == (n_runs, n_steps+1)` exactly per
config; same for the 5 other trajectory arrays. `n_runs =
len(spacing_um_values) * len(dt_s_values)`. Tested by test 20.

### 3.3 Metadata reproducibility (Y7 + Y13 + Y14)

- `metadata.config` is the verbatim input config (deterministic
  reproducibility — same config reproduces same result).
- `metadata.git_sha` is either the override (test path) OR
  detected via `_detect_git_sha` (best-effort subprocess with
  2.0s timeout; falls back to `"unknown"`).
- **NO `metadata.timestamp_iso` field** (Y13 — breaks
  reproducibility because two identical runs would differ).

Tests 12 + 13 + 14 cover.

### 3.4 ECM forward-propagation per run

For each `(spacing_um, dt_s)` pair:
- `ecm = _build_initial_ecm_at_spacing(config, spacing_um)`.
- Step 0 measures `init_result` with `dt_s=0` (NO ECM update;
  `init_result.updated_ecm.orientation_tensor` bytewise-equal
  to input but fresh per HB#1+#2 Y12).
- Steps 1..n_steps: each `step_result.updated_ecm` propagated
  to next step's `ecm` (forward dependency chain).

The forward chain is the harness's sole stateful invariant
across steps; each step otherwise stateless.

### 3.5 No FA mutation, no input config mutation

`adhesions` arguments (per-step) are read-only. `config` is
frozen+slots, immutable. `metadata.config` references the
input config (object identity, NOT copy).

### Status

**PASS**. Harness conserves nothing of its own (cross-product
of stateless sub-calls); 5 sub-call conservation invariants
inherited; trajectory shape exactly determined by config;
metadata reproducibility enforced by Y13 NO timestamp + Y14
git_sha override-capable.

---

## 4. Numerical sanity

### 4.1 Float precision

`np.float64` throughout for trajectories + scalar aggregates;
`np.int64` for `nx_by_run` / `ny_by_run` (cell counts —
exact integers).

### 4.2 IEEE roundoff allowance (Y18, NOT gate threshold)

Test 17 (`provides_v_active_bound_invariant_evidence`) uses
`max_bound_ratio ≤ 1.0 + 1e-12` as the assertion. The `1e-12`
is **IEEE float64 roundoff allowance** (standard numerical
analysis margin for accumulated floating-point error in HB#5
sum reductions), NOT a tunable gate threshold. The
**mathematical invariant** is `max_bound_ratio ≤ 1.0` per HB#5
lock §6.2 derivation; `1e-12` is the float64-implementation
margin around that invariant.

This is explicitly labeled in test 17 docstring + test
assertion + the lock §6 / Sanity Gate §6 (this section)
forward-guards against future "I'll just bump 1e-12 to 1e-9 to
make a failing run pass" — that would be a gate-threshold-
tuning anti-pattern explicitly forbidden by Y18.

### 4.3 No new tolerance constants beyond `1e-12` IEEE allowance

The 2 module constants imported from
`ecm_constitutive_response` (`K_ORIENT_PER_S`,
`TRACTION_REF_NN_PER_UM2`) are HB#1+#2 lock physics constants
(Y14 object-identity verified). Zero new harness-level
tolerance constants.

### 4.4 Per-call work + memory budget

- Per-run work: `n_steps + 1` Phase E v1 calls + 6 trajectory
  array writes per step.
- Per-run memory: 6 × `(n_steps + 1)` × float64 = ~5 KB per run
  for n_steps=100; 16 runs → ~80 KB total (negligible).
- Total harness runtime budget: `≤ 10 min CPU` on Laptop A5000
  16 GB baseline (Y15 — informational, NOT asserted in tests).

### 4.5 `_detect_git_sha` 2.0s subprocess timeout

Avoids hanging on systemic git issues (e.g., git not
installed, repo not initialized, file-system stalls). Falls
back to `"unknown"` on any failure path.

### 4.6 Stability

The harness has no time integration of its own (cross-product
loop only). HB#1+#2 owns the stability bound (unconditionally
stable in `dt` via convex-weight construction). The full
locked `dt_s ∈ {30, 60, 120, 300}` range is well within
HB#1+#2 stable regime (`K · S · dt ≈ 6e-3, 1.2e-2, 2.4e-2,
6e-2` at unit `S`; all `w ∈ [0, 1]`).

### Status

**PASS**. Float64 + int64 throughout; `1e-12` IEEE roundoff
explicitly labeled (NOT gate tunable); zero new harness-level
tolerances; memory + runtime budget within Laptop A5000 16 GB
baseline; Y18 forward-guard against threshold tuning.

---

## 5. Sign / sense check

### 5.1 V_active and bound non-negativity (HB#5 inheritance)

- `v_active_um2[run, step] ≥ 0` always (HB#5 sum-of-squares).
- `v_active_bound_um2[run, step] ≥ 0` always (4.5 ·
  active_area).
- `max_bound_ratio[run] ≥ 0` (V_active ≥ 0).
- `max_bound_ratio[run] ≤ 1.0 + 1e-12` (HB#5 mathematical
  invariant + Y18 IEEE allowance).

### 5.2 Convex weight + orientation delta sign (HB#1+#2 inheritance)

- `max_convex_weight[run, step] ∈ [0, 1]` always (HB#1+#2
  convex weight).
- `max_orientation_delta_frobenius[run, step] ≥ 0` always
  (Frobenius norm).
- `max_traction_norm_nN_per_um2[run, step] ≥ 0` always (vector
  norm).
- `neutral_multiplier_max_deviation_traj[run, step] == 0.0`
  always under Phase E v1 HB#4 neutral hard-wired.

### 5.3 Cell counts non-negative (sweep structural)

- `nx_by_run[run]` and `ny_by_run[run]` ≥ 1 (positive cell
  count by Y1 exact-divisibility validation).

### Status

**PASS**. All 7 trajectory + aggregate fields non-negative by
construction; harness performs no arithmetic that introduces
sign; HB#1+#2 + HB#5 + HB#4 sign conventions inherited.

---

## 6. Measurement-protocol consistency (Hard Rule 11) — central anchor

This Item 5 sweep harness is the **third Phase E composition
Hard Rule 11 measurement-protocol catch family** (per locked §6
references):

- Phase E v1 Y1 was the first (ECM-side evidence wording).
- HB#5 Y1 was the second (active-cells-only V_active).
- **Item 5 Y1 + Y4 + Y17 (this unit)**: physical-domain
  invariance under spacing sweep + ECM-feedback substep
  wording + Step 0 semantics.

### 6.1 Spacing axis sensitivity (Y1)

Physical domain + FA position invariant; only resolution
changes. ECM IC is a **physical-coordinate function** sampled
per grid (NOT array reuse). This avoids the spacing/IC/domain
conflation that would render sensitivity measurement
meaningless.

`_build_initial_ecm_at_spacing` samples the IC at every cell
center per `spacing_um`:
- `nx = int(round(domain_size_um_xy[0] / spacing_um))` (Y1
  exact-divisibility validated)
- `orientation_tensor` = `0.5 * I` broadcast to `(nx, ny, 2,
  2)` then `.copy()` (NOT a shared array reference)
- `stiffness_kpa`, `ligand_density`, `fiber_density`,
  `accumulated_traction_nNs_per_um2` constructed per
  `(nx, ny)`

The lock §1 forbidden list explicitly forbids array reuse
across spacings.

Tested by test 4 (indivisibility raise) + test 15 (structural
representation: all locked spacings present, derived nx/ny
match).

### 6.2 dt axis sensitivity (Y4 — ECM-feedback substep wording)

`dt_s ∈ {30, 60, 120, 300}` s is the **ECM-feedback
integration substep**, NOT the experimental imaging interval.
Per locked §0 wording + per CLAUDE.md docs:
- Production imaging interval = 15 min
- Pilot imaging interval = 60 min
- Item 5 sweep `dt_s` axis: 30 s to 300 s = sub-minute to
  sub-imaging-interval substep sensitivity

The locked range `{30, 60, 120, 300}` s is selected to span
sub-minute to sub-imaging-interval substep sensitivity.
Per-step `max_convex_weight` recording (in
`max_convex_weight_traj`) makes large-`dt` dynamics
compression visible (large `dt` → larger `w` per step → more
compressed evolution).

Tested by test 16 (analytical: at fixed spacing, step-1
`max_convex_weight` monotone non-decreasing in `dt` —
analytical from `w = 1 − exp(−K·S·dt)`, K, S > 0).

### 6.3 Step 0 semantics (Y17)

Step 0 = initial ECM measured against step-0 traction target
with `dt_s=0`. NOT a zero-traction baseline.

- The step-0 traction is NOT zero (rotating scenario at θ=0
  gives `(force_x, force_y) = (magnitude, 0)`; nonzero).
- `dt_s=0` ensures HB#1+#2 produces an identity update; HB#5
  measures the **initial misalignment** between `0.5*I`
  orientation IC and the step-0 target `n⊗n`.
- `max_convex_weight[run, 0] == 0` exactly (because
  `dt_s=0` → `w = -expm1(0) = 0`).
- `v_active_um2[run, 0]` may be nonzero (initial-target
  distance positive for `0.5*I` IC at θ=0).

This convention is what Item 5 sensitivity sweeps need: the
initial state's measurement is part of the evidence
trajectory, distinct from a zero-input baseline.

Tested by test 21
(`test_item_5_sweep_step_0_initial_against_target_not_zero_baseline`).

### 6.4 Evidence wording (Y8 + Y19, Phase E v1 Y1+Y2 sister)

`provides_*_evidence` test names (NOT `satisfies_*`); raw
observable summaries (NO pass/fail decision); explicit
"sensitivity not guaranteed monotone" wording in §0 + §1
forbidden list.

Per Codex `id=1563` review focus #5 + `id=1581` blocker
resolution + Phase E v1 Y2 sister: wording boundary enforced
by **two separate meta-tests**:
- Test 22: source/meta guard for no local `Item5*Error` class
  definition + no `failure_kind=` assignment (failure-kind
  discipline)
- Test 24 (NEW per Codex `id=1581`): source/meta guard for no
  satisfaction-claim wording (`"Item 5 satisfied"`,
  `"satisfies_item_5"`) in module/function/return-class
  docstrings + test names (wording-boundary discipline)

### 6.5 IEEE roundoff vs gate threshold (Y18 — central forward
guard)

`max_bound_ratio ≤ 1.0 + 1e-12` allowance is **IEEE float64
roundoff**, NOT a tunable gate. Test 17 docstring +
assertion + this Sanity Gate §4.2 + lock §6 explicitly label.

Future-regression scenario this guards: a hypothetical "the
ratio came out 1.0000001 from float-accumulation; let me bump
1e-12 to 1e-7 to make the test pass" would be a gate-threshold
tuning anti-pattern explicitly forbidden by Y18 + CLAUDE.md
Hard Rule "NEVER modify a gate's tolerance, normalisation, or
check window to make a failing run pass". The **mathematical
invariant** is `≤ 1.0` per HB#5 lock §6.2 derivation; `1e-12`
is the float64 numerical implementation margin around that
invariant.

### 6.6 Step 6 sister-gate-mirror application

Per memory `design_note_pre_commit_batch.md` 6-step, Step 6
sister-gate-mirror layered against Phase E v1 / HB#5 / HB#1+#2
/ B1 / Phase D / HB#3 / HB#4:

- **Code (validation ordering)**: harness has its own
  `_validate_config` defense-in-depth (Y3 + Y20 + Y21).
  Sub-call validation inherited from Phase E v1 (which inherits
  from HB#1+#2 + HB#3 + HB#4 + HB#5).
- **Design (typed dataclass + frozen+slots)**: harness follows
  Phase E v1 / HB#5 / HB#1+#2 / B1 / Phase D `@dataclass(frozen=True,
  slots=True)` precedent for all 3 dataclasses.
- **API surface (`__init__.py` exports)**: 4 new symbols at
  both `acs/v2/__init__.py` and `acs/v2/dynamics/__init__.py`
  per Y12. Test 23 enforces; private helpers stay private.
- **Failure-kind discipline** (Y22 typo in earlier draft —
  corrected to "no Y-number"; failure-kind decision is locked
  via §1 forbidden list + sister-pattern, not a numbered Y in
  the locked §2 trace which only goes to Y21): NO error class
  introduced (sister-pattern with Phase E v1 + HB#5 — no
  domain-specific failure modes for a thin harness wrapper);
  `_validate_config` raises plain `ValueError`. Meta-test 22
  forward-guards (no `Item5*Error` class + no `failure_kind=`
  assignment).
- **Wording-boundary meta-test**: harness inherits Phase E v1
  Y2 wording boundary (`provides_*_evidence` test names).
  **Per Codex `id=1581` blocker resolution**: a separate
  meta-test 24 (NEW; sister-pattern with Phase E v1 test 14)
  enforces no satisfaction-claim wording (`"Item 5 satisfied"`,
  `"satisfies_item_5"`) in module/function/return-class
  docstrings + test names. Test 22 stays focused on
  failure-kind discipline; test 24 focuses on wording
  discipline; the two are orthogonal forward guards.

### Status

**PASS**. Six layers of measurement-protocol guard:
- Locked §0 wording boundary (text-level evidence-only).
- Y1 physical-domain invariance enforced by
  `_build_initial_ecm_at_spacing` + exact-divisibility
  validation.
- Y4 ECM-feedback substep wording enforced by §0 + locked
  forbidden list.
- Y17 Step 0 semantics enforced by test 21.
- Y18 IEEE roundoff explicitly labeled + Y8/Y19 evidence-only
  test names + meta-test 22 (failure-kind) + meta-test 24
  (wording-boundary) per Codex `id=1581` separation.
- Step 6 sister-gate-mirror at all 4 layers + wording-boundary
  + IEEE-roundoff-vs-gate-threshold distinction.

This is the strongest measurement-protocol guard layer in the
v2 closed-loop ECM gate so far (third Phase E composition
catch family).

---

## 7. Magic-Number Block check

**Two locked numerical constants** per Y11 + Y18:

1. `0.5 * I` orientation IC (Y11): midpoint between zero
   orientation and unit diagonal schema bound; symmetric
   non-target IC. Magic-Number Block 3-test PASS:
   - Derivable: midpoint of schema bound `[-1, 1]` (specifically
     diagonal bound `[0, 1]` for `n⊗n` rank-1 targets); not
     fitted; chosen for symmetric non-target.
   - Grid-invariant: `0.5*I` does not depend on `dx`, `dy`, or
     `grid_n`; only the broadcasting shape changes per
     `(nx, ny, 2, 2)`.
   - Not fitted: chosen analytically before any A/A₀
     comparison or test target gating.
2. `1e-12` IEEE roundoff allowance (Y18): float64 numerical
   margin in `max_bound_ratio ≤ 1.0 + 1e-12` assertion. NOT a
   gate threshold — labeled explicitly in test 17 + this §6.5.
   Magic-Number Block 3-test PASS:
   - Derivable: standard IEEE 754 float64 round-off scale
     after sum reductions; from numerical analysis literature.
   - Grid-invariant: `1e-12` is an absolute float64 quantity,
     not a function of `dx`/`dt`/`grid_n`.
   - Not fitted: chosen analytically before any test target.

All other physics constants (`K_ORIENT_PER_S`,
`TRACTION_REF_NN_PER_UM2`) are **re-imported** from HB#1+#2
lock per Y14 (NOT redefined; meta-test 16 sister-pattern would
apply but the harness module test catalog does NOT add a
Y14-equivalent meta-test because the harness only USES the
constants, doesn't re-import them with potential redefinition
risk; the import statement itself is the only way to obtain
them).

### Status

**PASS**. Two locked harness-level constants (`0.5*I` IC,
`1e-12` IEEE allowance), both Magic-Number Block 3-test
compliant. All other constants re-imported from HB#1+#2 lock.

---

## 8. Test catalog (25 tests; 23 per locked §4 + 1 new per Codex `id=1581` blocker (test 24 wording-boundary) + 1 new per Codex `id=1591` blocker (test 25 near-indivisible-inside-old-tolerance))

Owned by `tests/test_v2_closed_loop_phase_e_sweep.py` (not yet
committed). Each test maps to a locked invariant in
`docs/v2_item_5_sweep_harness_locked.md` §4.

### Validation tests (11)

1. `test_item_5_sweep_empty_spacing_values_raises`
2. `test_item_5_sweep_empty_dt_values_raises`
3. `test_item_5_sweep_negative_or_nonfinite_value_raises`
4. `test_item_5_sweep_indivisible_domain_raises` (Y1)
5. `test_item_5_sweep_unsupported_traction_scenario_raises`
   (Y20)
6. `test_item_5_sweep_zero_n_revolutions_raises` (Y20)
7. `test_item_5_sweep_negative_n_revolutions_raises` (Y20)
8. `test_item_5_sweep_bool_n_revolutions_raises` (Y20 — Python
   `bool ⊂ int` trap, HB#1+#2 sister)
9. `test_item_5_sweep_invalid_domain_size_raises` (Y21)
10. `test_item_5_sweep_invalid_origin_raises` (Y21)
11. `test_item_5_sweep_fa_position_outside_domain_raises` (Y21)

### Determinism + metadata (3)

12. `test_item_5_sweep_metadata_uses_git_sha_override_when_provided`
    (Y14)
13. `test_item_5_sweep_metadata_falls_back_to_unknown_on_git_failure`
    (Y14 — monkeypatch subprocess to raise)
14. `test_item_5_sweep_no_timestamp_in_metadata` (Y13 — `not
    hasattr(metadata, "timestamp_iso")` and similar)

### Evidence-providing tests (4) — Phase E v1 Y2 sister-pattern

15. `test_item_5_sweep_provides_grid_spacing_sensitivity_evidence`
    (Y19 structural representation: all locked spacings
    present; derived `nx_by_run` / `ny_by_run` match `domain
    / spacing` exactly; trajectories have correct shape per
    spacing run)
16. `test_item_5_sweep_provides_dt_substep_sensitivity_evidence`
    (Y19 analytical: at fixed spacing, step-1
    `max_convex_weight` monotone non-decreasing in `dt` —
    analytical consequence of `w = 1 − exp(−K·S·dt)`)
17. `test_item_5_sweep_provides_v_active_bound_invariant_evidence`
    (Y18 — `max_bound_ratio ≤ 1.0 + 1e-12` everywhere; `1e-12`
    is IEEE roundoff, NOT gate threshold)
18. `test_item_5_sweep_provides_neutral_multiplier_preservation_evidence`
    (`max_neutral_multiplier_deviation == 0.0` everywhere —
    Phase E v1 HB#4 neutral hard-wired invariant)

### Boundary + step-0 (3)

19. `test_item_5_sweep_max_bound_ratio_zero_when_no_active_cells`
    (Y16 — degenerate all-zero-traction scenario; no division
    by zero)
20. `test_item_5_sweep_trajectory_shape_n_runs_n_steps_plus_1`
    (Y6 — explicit shape verification)
21. `test_item_5_sweep_step_0_initial_against_target_not_zero_baseline`
    (Y17 — assert `v_active[:, 0] > 0` for misaligned IC;
    `max_convex_weight[:, 0] == 0` exact;
    `max_orientation_delta[:, 0] == 0` exact)

### Composition guard + exports + wording boundary (3)

22. `test_item_5_sweep_no_local_failure_kinds_or_phase_e_error_redefinition`
    (sister-pattern with Phase E v1 test 3 source/meta guard
    — assert no `Item5*Error` class definition + no
    `failure_kind =` assignment in module body) — failure-kind
    discipline ONLY; wording discipline owned by test 24.
23. `test_item_5_sweep_exports_through_both_init` —
    `PhaseEV1Item5SweepConfig`, `PhaseEV1Item5SweepMetadata`,
    `PhaseEV1Item5SweepResult`,
    `run_phase_e_v1_sensitivity_sweep` importable from BOTH
    `acs.v2.dynamics` AND `acs.v2`; private `_default_*`,
    `_build_*`, `_validate_*`, `_detect_*`, `_compute_*`,
    `_record_*` helpers NOT exported.
24. **`test_item_5_sweep_module_doc_does_not_overclaim_item_5_satisfaction`**
    (NEW per Codex `id=1581` blocker; sister-pattern with
    Phase E v1 test 14 wording-boundary forward guard) —
    module docstring + function docstring + return-class +
    config-class + metadata-class docstrings do NOT contain
    `"Item 5 satisfied"`, `"satisfies_item_5"`, or
    `"satisfies_*_item_5"`. AST walk on module + string-search
    over `inspect.getsource(...)` + iteration over test
    function names in `tests/test_v2_closed_loop_phase_e_sweep.py`
    asserting none start with `test_*_satisfies_item_5_*` or
    contain `satisfies_item_5`. Failure-kind discipline owned
    by test 22; wording discipline owned by this test.

---

## 9. Gate verdict

| Item | Status | Notes |
|---|---|---|
| §1 Dimensional | **PASS** | No new unit chain; sub-call chains preserved; 0.5*I IC + 1e-12 IEEE roundoff Magic-Number Block compliant |
| §2 Boundary | **PASS** | 8 boundary classes locked; 15 of 25 tests covering validation (incl Codex `id=1591` near-indivisible-inside-old-tolerance) + Step 0 + zero-bound case |
| §3 Conservation | **PASS** | Harness conserves nothing of its own; 5 sub-call invariants inherited; trajectory shape exactly per config; metadata reproducibility (Y13 + Y14) |
| §4 Numerical | **PASS** | Float64 + int64; 1e-12 IEEE roundoff explicitly labeled NOT gate tunable (Y18 forward guard); zero new harness tolerances |
| §5 Sign | **PASS** | All 7 trajectory + aggregate fields non-negative by construction |
| §6 Measurement-protocol | **PASS** | Six-layer guard (locked §0 wording + Y1 physical-domain invariance + Y4 ECM-feedback substep + Y17 Step 0 + Y18 IEEE roundoff label + Y8/Y19 evidence-only); third Phase E composition Hard Rule 11 catch family |
| §7 Magic-Number Block | **PASS** | Two locked constants (0.5*I IC, 1e-12 IEEE), both 3-test compliant |

**Overall**: gate PASS. impl-work is clear to commit
`acs/v2/dynamics/closed_loop_phase_e_sweep.py` + 4 export
updates + `tests/test_v2_closed_loop_phase_e_sweep.py` after
Codex review of this Sanity Gate doc, mirroring HB#3 / HB#4 /
Phase D / B1 / HB#1+#2 / HB#5 / Phase E v1 Sanity Gate review
precedent.

### Outstanding before code lands

- Codex review of this Sanity Gate doc per locked §7 6 review
  focus areas:
  1. Physical-domain invariance (Y1) — spacing sweep with
     fixed domain + FA position; ECM IC sampled per grid (NOT
     array reuse)
  2. Typed result/metadata throughout (Y7+Y13) — no loose
     dicts; no wall-clock timestamp; git_sha override +
     fallback-safe detection
  3. Evidence-only wording + IEEE roundoff labeling (Y8 + Y18
     + Y19) — `provides_*_evidence` test names; `1e-12`
     explicitly NOT gate threshold; no monotone-in-spacing
     claim
  4. Step 0 semantics (Y17) — initial ECM against step-0
     traction target; `max_convex_weight[:, 0] == 0` exact
  5. Schema-concrete builders with `fa.validate()` (Y9+Y10+Y11)
     — `_build_rotating_uniform_single_fa` uses all
     FocalAdhesionState fields explicitly; `_build_initial_ecm_at_spacing`
     samples physical-coordinate IC
  6. 4 public exports only (Y12) — through both
     `acs.v2.dynamics` and `acs.v2`; private helpers stay
     private
- Code commit module-docstring intent guard (per Phase E v1
  precedent + locked §1 forbidden list paraphrase): describe
  Item 5 wording boundary in **paraphrased** form without
  inserting forbidden literal strings (`"satisfies_item_5"`,
  `"Item 5 satisfied"`, etc.) which **meta-test 24** must
  reject from module source / docstrings / test names per the
  test 22 (failure-kind) / test 24 (wording-boundary)
  separation locked at Codex `id=1581`. Sister-pattern with
  Phase E v1 `id=1546` + `id=1549` paraphrased intent guard
  resolution.
- All 25 tests must pass at first commit; no `TODO test_X`
  placeholders.
- `pytest tests/test_v2_closed_loop_phase_e_sweep.py` +
  combined sister-gate regression (HB#3 + HB#4 + Phase D +
  6.3a + 6.3b + HB#1+#2 + HB#5 + Phase E v1 + Item 5 sweep
  harness) must PASS at commit time.

---

## 10. References

- Item 5 sweep harness lock:
  `docs/v2_item_5_sweep_harness_locked.md` (commit `c5634d5`).
- Item 5 sweep harness brief (superseded by lock):
  `docs/v2_item_5_sweep_harness_brief.md` (commit `4f0efdd`).
- Parent locked plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` (Item 5
  sweep is a separate harness from Phase E v1, per parent §1).
- Sister locks + Sanity Gates:
  - `docs/v2_phase_e_composition_locked.md` +
    `docs/v2_phase_e_composition_sanity_gate.md` (Phase E v1
    — the harness target)
  - `docs/v2_hard_blocker_5_lyapunov_metric_locked.md` +
    `docs/v2_hard_blocker_5_lyapunov_metric_sanity_gate.md`
    (HB#5 — V_active + bound)
  - `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`
    + `docs/v2_hard_blocker_1_2_constitutive_response_sanity_gate.md`
    (HB#1+#2 — K_ORIENT, TRACTION_REF, expm1 stability)
  - `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md` +
    `docs/v2_fa_to_ecm_scattering_sanity_gate.md`
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` +
    `docs/v2_ecm_to_fa_bias_sanity_gate.md`
  - `docs/v2_phase_d_no_op_scaffolding_locked.md` +
    `docs/v2_phase_d_no_op_scaffolding_sanity_gate.md`
  - `docs/v2_focal_adhesion_dynamics_result_typed_locked.md` +
    `docs/v2_focal_adhesion_dynamics_result_typed_sanity_gate.md`
    (B1)
- Verified at source for Y9 (FocalAdhesionState schema):
  `acs/v2/focal_adhesion.py` lines 51-85.
- Sister implementation precedent (composition wrapper):
  `acs/v2/dynamics/closed_loop_phase_e.py`.
- CLAUDE.md Performance Budget (Laptop A5000 16 GB baseline +
  CPU-only harness acceptable per Y15).
- CLAUDE.md Hard Rule "NEVER modify a gate's tolerance,
  normalisation, or check window to make a failing run pass"
  (Y18 forward guard rationale).
- CLAUDE.md Hard Rule 1 (NEVER fit parameters to PI's
  experimental data).
- Memory rules informing this gate:
  - `design_note_pre_commit_batch.md` (6-step including Step 6
    sister-gate-mirror, per Codex `id=1428`).
  - `rule10_unit_derivation_in_docs.md`.
  - `hard_rule_11_wording_boundary_meta_test.md` (this lock is
    the **third** Phase E composition Hard Rule 11 catch
    family — Phase E v1 Y1 was first; HB#5 Y1 was second; Item
    5 Y1+Y4+Y17 here is the third).
  - `feedback_aggressive_design_debate.md`.
- Design-discussion thread: `id=1564` → `1575` (4-round
  adversarial lock + seal ack); `id=1576` design-discussion
  status update.
- impl-work review thread (this Sanity Gate cycle): `id=1578`
  (design-discussion → impl-work dispatch) → impl-work Sanity
  Gate doc commit (this file) → Codex impl-work review → code
  commit.
