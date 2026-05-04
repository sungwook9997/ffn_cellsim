# V2 HB#4-Active Step 1 — Sanity Gate (impl-work review)

**Date**: 2026-05-05 KST
**Author**: implementation-work Claude
**Lock**: `docs/v2_hard_blocker_4_active_locked.md` (commit `fed370b`)
**Brief**: `docs/v2_hard_blocker_4_active_brief.md` (commit `013e161`)
**Source**: lock §3 Sanity Gate skeleton (6 items) + lock §4 21-test
catalog + Codex `id=1612` SEAL ack 6 review-focus items
**Verdict** (this doc): **all 6 Sanity Gate items PASS**, ready for
impl-work Codex review per (B-trigger) precedent established by
HB#3 / HB#4 v1 / Phase D / B1 / HB#1+#2 / HB#5 / Phase E v1 / Item 5
sweep harness Sanity Gate cycles.

This doc is the implementation-work-pane review surface preceding
code commit. It mirrors the Phase E v1 / Item 5 sweep harness
Sanity Gate structure, adapted for the active-bias-law unit:
single-function pure read-only module modification, NOT new module.

---

## 0. Cross-reference to Codex `id=1612` 6 review-focus areas

| # | Codex review focus (per `id=1612` SEAL) | Sanity Gate section | Lock Y-anchor |
|---|---|---|---|
| 1 | Deviatoric form correctness (`Q = T - 0.5·trace(T)·I`; isotropic → score = 0 → multiplier = 1.0 exact) | §2.6 + §6.1 | Y1 |
| 2 | Schema-tight √2 example correctness (`T = [[1,1],[1,-1]]` traceless, eigenvalues ±√2, eigenvector at θ=π/8) | §2.6 + §6.2 | Y10 |
| 3 | Validation order `ecm.validate() → _validate_k_active → sampler` | §1 + §6.3 | Y12 |
| 4 | Schema-corrected fields (`fa_ids`, `multipliers_per_fa`, `rate_names`, `sampled_diagnostics`, `diagnostics_dict`) | §3 | Y5 |
| 5 | Failure-kind extension (`k_active_invalid`, `active_multiplier_non_finite`) + finite defensive guards on score and multiplier | §4 + §6.5 | Y9 + Y11 |
| 6 | Function-scoped AST + string guard (no `effective_stiffness` / mechanosensing fields in active body); single export with private constant `_DEVIATORIC_SCORE_BOUND` | §6.4 + §6.6 | Y8 + Y14 |

---

## 1. Function-level summary

### Public function (per lock §1)

```python
def compute_ecm_to_fa_bias_active(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
    *,
    k_active: float,                 # REQUIRED, no default (Y3)
) -> ECMToFABiasResult:
    ...
```

### Validation chain (Y12, locked order)

1. `ecm.validate()` — ECM substrate boundary validation (cheap structural)
2. `_validate_k_active(k_active)` — parameter discipline (cheap, no FA dependency)
3. `sample_ecm_at_fa_positions(adhesions, ecm)` — most expensive; may also raise on FA position

The order is enforced in code by sequencing; **invalid `k_active` raises BEFORE `sample_ecm_at_fa_positions` is called even with empty adhesions** — guarded by meta-test 19 (lock §4 #19) that monkeypatches the sampler to track invocation.

### Return type (Y5 + Y15 grandfathering)

`ECMToFABiasResult` (HB#4 v1 grandfathered, NOT new dataclass) with 5
fields verified at `acs/v2/dynamics/ecm_to_fa_bias.py:127-141`:

- `fa_ids: tuple[str, ...]`
- `multipliers_per_fa: np.ndarray` of shape `(n_fa, 3)` (3 rates, all same per Y4)
- `rate_names: tuple[str, ...]` = `("k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s")`
- `sampled_diagnostics: Optional[ECMSampledAtFAs]` — non-None even for empty FA list
- `diagnostics_dict: dict[str, float | int | str]` — exactly 10 keys (Y6)

**Anti-pattern blocked**: pseudocode field name `sampled` (singular) or
`multipliers` (no `_per_fa` suffix) is a Phase E v1 B2 sister-pattern
catch class; locked code uses verified canonical names only.

### Multiplier law (Y1 + Y2 + Y4)

For each FA `i`:

```
n_fa = traction_force / |traction_force|     # branch-defined zero-traction → exact neutral
T_local = sampled.orientation_tensor[i]      # (2, 2) symmetric, sampled at FA position
trace_T = T_local[0,0] + T_local[1,1]
Q = T_local - 0.5 * trace_T * I              # deviatoric: alignment-only
score = n_fa @ Q @ n_fa                      # ∈ [-√2, +√2] under schema |T_ij| ≤ 1
multiplier = exp(k_active * score)           # bounded positive exponential map
multipliers_per_fa[i, :] = multiplier        # all 3 rates same (Y4 Q2=a)
```

---

## 2. Sanity Gate (6 items)

### 2.1 Units (PASS)

| Symbol | Unit | Source |
|---|---|---|
| `traction_force_nN_xy` | nN | `FocalAdhesionState` schema |
| `traction_norm` | nN | derived `np.linalg.norm` |
| `n_fa_unit` | dimensionless unit vector | derived |
| `T_local`, `Q` | dimensionless | `ECMSubstrateState.orientation_tensor` schema (`|T_ij| ≤ 1`) |
| `alignment_score` | dimensionless | `n.T @ Q @ n` |
| `k_active` | dimensionless | caller-supplied (no biological default per Y3) |
| `multiplier` | dimensionless | `exp(k_active · score)`, ≥ 0 |
| `_DEVIATORIC_SCORE_BOUND` | dimensionless | `√2 ≈ 1.4142` |

**Hard Rule 10 inline derivation** (per CLAUDE.md `rule10_unit_derivation_in_docs` memory; mandatory in §6.1 below): `n.T @ Q @ n` is dimensionless because `n` is dimensionless unit + `Q` inherits the dimensionless schema bound from `T`. `multiplier` is dimensionless because `k_active · score` is dimensionless × dimensionless. No kPa↔nN/μm² unit comparison present (function never reads `stiffness_kpa` per Y8 forbidden list). **Unit chain coherent.**

### 2.2 Boundary (PASS)

11 of 21 tests cover boundary cases (lock §4):

- Empty adhesions (test 14): `multipliers.shape == (0, 3)`; sampler still called → `sampled_diagnostics` non-None; `diagnostics_dict` populated with default zero counts (`max/min_multiplier == 1.0`, `max/min_alignment_score == 0.0`, `zero_traction_count == 0`)
- All-zero traction (tests 12+13): `alignment_score = 0.0` exact, `multiplier = 1.0` exact, `zero_traction_count == n_fa`
- Isotropic ECM `T = a·I` (test 6): `Q = 0` → `alignment_score = 0` → `multiplier = 1.0` exact for ALL FAs (Y1 anti-pattern guard; Item 5 IC `T = 0.5·I` automatically neutral)
- Schema-saturating ECM `T = [[1,1],[1,-1]]` traceless rank-2 (test 7): score reaches `+√2` at angle `θ = π/8`
- `k_active = ε` small (covered by test 9 aligned-boost via small + positive value): multiplier ≈ `1 + k·score` Taylor
- `k_active = 1000` overflow (test 4): `_validate_k_active` raises `k_active_invalid`. **Threshold derivation** (Codex `id=1618` correction): `exp(k · √2) > float64.max` iff `k > log(float64.max) / √2 ≈ 709.78 / 1.4142 ≈ 501.892`. Direct computation: `exp(100·√2) ≈ 2.62e+61` is finite (well below `float64.max ≈ 1.79e+308`); `exp(502·√2)` overflows; `exp(1000·√2) = exp(1414.2)` overflows clearly. Test fixture uses `k_active = 1000.0` to be safely above the ~502 threshold. Lock §3 #2 amended in same Codex `id=1618` BLOCKER fix.
- `k_active = 0` (test 1): raises `k_active_invalid`
- `k_active < 0` (test 1): raises `k_active_invalid`
- `k_active = NaN/inf` (test 2): raises `k_active_invalid`
- `k_active = True` Python bool ⊂ int trap (test 3): raises `k_active_invalid`
- Invalid ECM (test 5): `ecm.validate()` raises before `_validate_k_active`

**No silent NaN/inf propagation**: defensive finite guards on score and multiplier (Y11) raise `active_multiplier_non_finite` if `ecm.validate()` somehow lets a corrupted orientation_tensor through.

### 2.3 Conservation (PASS, function is pure read-only)

The function performs ZERO mutation:

- `adhesions` parameter: read-only iteration, never mutated (input is `tuple[FocalAdhesionState, ...]`, FocalAdhesionState is `frozen=True, slots=True` schema)
- `ecm` parameter: read-only `ecm.validate()` + read via `sample_ecm_at_fa_positions` (HB#4 v1 sister-pattern)
- `multipliers`, `alignment_scores` numpy arrays: freshly allocated `np.ones(...)` / `np.zeros(...)` per call; no shared state
- Returned `ECMToFABiasResult`: `frozen=True, slots=True` (HB#4 v1 grandfathered, immutable)

**No conservation invariant applies** beyond purity (this is a read-only diagnostic-style function); no momentum/energy/mass to conserve at this layer. The composition-level conservation discipline (FA→ECM→FA closure) is Phase E v2 step 2 scope — explicitly out of scope per Codex `id=1604` (a1) deferral.

### 2.4 Numerical (PASS)

- **Float precision**: `np.float64` enforced throughout (`dtype=np.float64` on `multipliers`, `alignment_scores`, `traction`, `np.eye(2)`)
- **Branch-defined zero-traction** (Y7): `traction_norm > 0.0` exact (no epsilon); HB#1+#2 Y4 sister-pattern. `==` is unsafe for floats but `> 0.0` strictly distinguishes positive from non-positive without epsilon ambiguity.
- **Defensive finite guards** (Y11): `np.isfinite(alignment_score)` AND `np.isfinite(multiplier)` checks BEFORE returning. Catches corrupted `orientation_tensor` even if `ecm.validate()` didn't trigger (defense in depth).
- **No new tolerance constant** introduced. Tests use:
  - Exact equality `== 1.0` for isotropic-IC neutral (test 6) — exact arithmetic since `Q = 0` makes `score = 0.0` exact
  - `rel=1e-10` for the saturation example test 7 — IEEE roundoff allowance for `√2` computed via `np.sqrt(2.0)` and `θ = π/8` cosine/sine
  - **NO ad-hoc 1e-9 or similar magic-number tolerance** (Codex `id=1591` Item 5 sweep precedent — exact-divisibility lesson)
- **No CFL bound applies**: this is not a time-stepping function; it's a per-step diagnostic.

### 2.5 Sign (PASS)

| Quantity | Sign / range | Direction interpretation |
|---|---|---|
| `traction_norm` | ≥ 0 | magnitude, never negative |
| `n_fa_unit` | norm = 1 | direction of FA traction |
| `Q = T - 0.5·trace(T)·I` | traceless symmetric | deviatoric: alignment-only (isotropic factored out) |
| `alignment_score` | `∈ [-√2, +√2]` | **+** when ECM aligned with FA traction direction → BOOST; **−** when orthogonal → SUPPRESS; **0** when isotropic ECM or zero traction → NEUTRAL |
| `multiplier = exp(k_active · score)` | `> 0` always (exponential of finite real); `> 1` for boost; `< 1` for suppress; `= 1` for neutral (`score = 0`) | multiplicative on FA rates: boost increases rate, suppress decreases rate |

**Sense check**: aligned ECM with strong fiber direction matching FA traction → boost matches biological expectation (Hall 2016 + Trichet 2012 cross-reference in lock §6 references). Orthogonal ECM → suppress matches expectation. Isotropic ECM → neutral matches expectation.

### 2.6 Measurement-protocol (PASS — Hard Rule 11 central anchor)

This is the central Sanity Gate item for HB#4-active step 1 — **the fourth Phase E composition Hard Rule 11 catch family** (Phase E v1 Y1, HB#5 Y1, Item 5 Y1+Y4+Y17, **HB#4-active Y1**).

The measurement is the bias multiplier per FA, evaluated point-wise at FA position. The analytical proof must cover the same evaluation point.

**Deviatoric bound derivation** (Hard Rule 10 inline derivation per `rule10_unit_derivation_in_docs` memory — mandatory):

Setup: `T` is symmetric 2×2 with schema invariant `|T_ij| ≤ 1` (no PSD constraint required; verified at `acs/v2/ecm_substrate.py:36` for the `_ORIENTATION_BOUND = 1.0` constant and `:115` for the validator message — lock §6 cited `:42` which is the dataclass body offset; the canonical source-of-truth lines for the invariant are `:36` and `:115`).

Step-by-step deviatoric construction:
- `trace_T = T_xx + T_yy ∈ [-2, +2]`
- `Q = T - 0.5·trace_T·I`, by element:
  - `Q_xx = T_xx - 0.5·(T_xx + T_yy) = 0.5·(T_xx - T_yy)`
  - `Q_yy = T_yy - 0.5·(T_xx + T_yy) = -0.5·(T_xx - T_yy) = -Q_xx`
  - `Q_xy = Q_yx = T_xy` (preserved off-diagonal)
  - `Q_xx + Q_yy = 0` ✓ (traceless by construction)

For unit `n = (cos θ, sin θ)`:
- `n.T @ Q @ n = Q_xx · cos²θ + 2 · Q_xy · cosθ · sinθ + Q_yy · sin²θ`
- Substituting `Q_yy = -Q_xx`:
  - `= Q_xx · (cos²θ - sin²θ) + 2 · Q_xy · cosθ · sinθ`
  - `= Q_xx · cos(2θ) + Q_xy · sin(2θ)` (using double-angle identities)
  - `= 0.5·(T_xx - T_yy)·cos(2θ) + T_xy·sin(2θ)`

Bound (Cauchy-Schwarz / amplitude of `a·cos(2θ) + b·sin(2θ)`):
- `|n.T @ Q @ n| ≤ √(((T_xx - T_yy)/2)² + T_xy²)`
- Schema: `|T_ij| ≤ 1` ⇒ `|T_xx - T_yy| ≤ 2` ⇒ `((T_xx - T_yy)/2)² ≤ 1`; `T_xy² ≤ 1`
- Sum: `((T_xx - T_yy)/2)² + T_xy² ≤ 2`
- ⇒ `|n.T @ Q @ n| ≤ √2 ≈ 1.4142`

**Saturation example** (Y10 corrected math):
- `T = [[1, 1], [1, -1]]` is **already traceless** (`T_xx + T_yy = 1 + (-1) = 0`), so `Q = T`
- Eigenvalues from characteristic polynomial: `det(Q - λI) = (1-λ)(-1-λ) - 1 = λ² - 1 - 1 = λ² - 2 = 0` ⇒ `λ = ±√2`
- Eigenvector for `+√2`: solve `(Q - √2·I) v = 0`, gives angle `θ = π/8` (since `tan(2θ) = T_xy/(0.5·(T_xx - T_yy)) = 1/1 = 1` ⇒ `2θ = π/4` ⇒ `θ = π/8`)
- Hence `n.T @ Q @ n = +√2` at `n = (cos(π/8), sin(π/8))` — saturates the bound
- This is the test fixture for meta-test 7 (lock §4): `T = [[1,1],[1,-1]]`, `n = (cos(π/8), sin(π/8))`, expected `multiplier == exp(k_active · √2)` within `rel=1e-10`

**Anti-pattern guard** (Y1 — the central catch):
- **Raw Rayleigh** `n.T @ T @ n` would give `0.5` for `T = 0.5·I` regardless of `n` direction → all FAs uniformly boosted → false alignment claim
- **Deviatoric Rayleigh** `n.T @ Q @ n` for `T = a·I` gives `Q = 0` ⇒ `score = 0` ⇒ `multiplier = 1.0` exact → correctly neutral
- Test 6 enforces this: `T = 0.5·I` (Item 5 IC), assert `multiplier == 1.0` exactly for any traction direction

**Anti-pattern caught at round 1 → round 2 transition** (lock §2 Y1 trace; Codex `id=1608` C1).

**Off-proof points checked**:
- The analytical proof gives the BOUND. The measurement evaluates `n.T @ Q @ n` at one specific `n_fa` direction per FA — the analytical bound covers ALL `n` directions, so any single evaluation is bounded.
- The proof gives `|score| ≤ √2` for ALL valid schema-bounded `T`. The test (test 7) explicitly evaluates at the saturating angle to verify the bound is tight.
- For `T = a·I` isotropic, `Q = 0` is **exact arithmetic** (no floating-point roundoff in `T - 0.5·trace(T)·I` when `T_xx == T_yy`, since both subtraction terms are equal in finite precision when `T` is loaded from float64 storage and `trace_T = T_xx + T_yy = 2·T_xx`). Test 6 asserts `multiplier == 1.0` exact (no `rel=` tolerance).

**Closure mechanism** (composition-level forward note, NOT this lock's burden):
- HB#1+#2 evolves `T` toward rank-1 `n_traction ⊗ n_traction`
- `T` becomes deviatoric-aligned with `n_traction` (the FA's pulling direction)
- Active law sees increasing `score` over time → multiplier grows above neutral → FA rates rescaled upward
- Phase E v2 step 2 composition will exercise this; HB#4-active step 1 only provides the bias mapping

**Item 5 IC compatibility**:
- Item 5 sweep harness uses `0.5·np.eye(2)` initial orientation
- Under deviatoric form: `Q(0.5·I) = 0` ⇒ score = 0 ⇒ multiplier = 1.0 exactly
- Therefore Item 5 sweep harness can be augmented with HB#4-active in a future cycle without polluting the "scaffolding starts neutral" baseline (forward-compat note; not this commit's scope)

---

## 3. Schema-corrected fields (Y5)

Codex `id=1608` C5 caught that the brief's pseudocode used field names
`sampled` and `diagnostics` that do NOT exist in the actual schema.
Verified at source `acs/v2/dynamics/ecm_to_fa_bias.py:127-141` the
canonical field names:

| Pseudocode (forbidden) | Actual schema (locked) |
|---|---|
| `sampled` | `sampled_diagnostics` |
| `diagnostics` | `diagnostics_dict` |
| `multipliers` | `multipliers_per_fa` |

**Locked construction** (per lock §1 + §2 Y5):

```python
return ECMToFABiasResult(
    fa_ids=sampled.fa_ids,
    multipliers_per_fa=multipliers,
    rate_names=tuple(RATE_NAMES),
    sampled_diagnostics=sampled,
    diagnostics_dict=diagnostics_dict,
)
```

**Empty FA list invariant**: `sample_ecm_at_fa_positions((), ecm)` is
called even for empty adhesions (HB#4 v1 sister-pattern), so
`sampled_diagnostics` is **non-None** with shape-consistent zero-row
arrays. Test 14 asserts:
- `multipliers_per_fa.shape == (0, 3)`
- `sampled_diagnostics is not None`
- `diagnostics_dict["n_adhesions"] == 0`
- `diagnostics_dict["max_multiplier"] == 1.0` and `min_multiplier == 1.0` (default per lock §1 body branch)

**Sister-pattern with Phase E v1 B2 catch** (Codex `id=1543`): the
Phase E v1 unit also caught a pseudocode field name divergence
(`multipliers` → `multipliers_per_fa`); HB#4-active inherits the
verify-at-source discipline.

---

## 4. Failure-kind extension (Y9)

Locked extension (lock §1 + §2 Y9):

```python
# acs/v2/dynamics/ecm_to_fa_bias.py
FAToECMBiasFailureKind = Literal[
    "fa_bias_position_outside_ecm_grid",     # existing (HB#4 v1)
    "non_finite_fa_position",                # existing (HB#4 v1)
    "k_active_invalid",                      # NEW (HB#4-active step 1)
    "active_multiplier_non_finite",          # NEW (HB#4-active step 1)
]
```

**Sister-pattern with HB#4 v1 + B1 typed-schema discipline**: the
literal extension preserves the existing 2 kinds (sampler-related)
and adds 2 active-law-specific kinds. Sampler failures from
`sample_ecm_at_fa_positions` continue to surface with
`fa_bias_position_outside_ecm_grid` / `non_finite_fa_position` as
before; `compute_ecm_to_fa_bias_active` does NOT remap or wrap those
errors.

**Defensive finite guards (Y11)** are dual-layer:

1. **Score guard**: `np.isfinite(alignment_score)` BEFORE `np.exp(...)`
   — catches corrupted `orientation_tensor` even though `ecm.validate()`
   should have prevented this (defense in depth)
2. **Multiplier guard**: `np.isfinite(multiplier)` AFTER `np.exp(...)`
   — catches `np.exp` overflow despite `_validate_k_active` checking
   `np.exp(k_active · √2)` upper bound (defense in depth — protects
   against asymmetric overflow paths)

Both guards raise `FAToECMBiasError` with `failure_kind ==
"active_multiplier_non_finite"` (single kind, not two) per Y11
resolution.

**Test 18** (lock §4 #18) monkeypatches `sampled.orientation_tensor`
to contain NaN/inf and asserts the guard fires with the correct
failure_kind.

---

## 5. Test catalog mirror (21 tests)

Mirroring lock §4 with cross-reference to Y items:

### Validation (5)

| # | Test | Y |
|---|---|---|
| 1 | `test_compute_ecm_to_fa_bias_active_k_active_zero_or_negative_raises` | Y3 |
| 2 | `test_compute_ecm_to_fa_bias_active_k_active_nonfinite_raises` | Y3 |
| 3 | `test_compute_ecm_to_fa_bias_active_k_active_bool_raises` | Y3 (Python `bool ⊂ int` trap) |
| 4 | `test_compute_ecm_to_fa_bias_active_k_active_too_large_raises_non_finite_bound` | Y3 |
| 5 | `test_compute_ecm_to_fa_bias_active_invalid_ecm_raises_at_validate` | Y12 (validation order) |

### Deviatoric form (3, Y1+Y10 corrected math)

| # | Test | Y |
|---|---|---|
| 6 | `test_compute_ecm_to_fa_bias_active_isotropic_orientation_returns_neutral_multiplier` (`T = 0.5·I` Item 5 IC; assert `multiplier == 1.0` exact for any traction direction) | Y1 anti-pattern guard |
| 7 | `test_compute_ecm_to_fa_bias_active_schema_tight_bound_score_approaches_sqrt2` (`T = [[1,1],[1,-1]]`, `n = (cos(π/8), sin(π/8))`; assert `multiplier == exp(√2 · k_active)` within `rel=1e-10`; assert `max_alignment_score ≈ √2`) | Y10 |
| 8 | `test_compute_ecm_to_fa_bias_active_rank1_aligned_score_half` (`T = n⊗n`; assert `score = 0.5`) | Y1 (rank-1 minus isotropic-half) |

### Multiplier behavior (3)

| # | Test | Y |
|---|---|---|
| 9 | `test_compute_ecm_to_fa_bias_active_aligned_target_boosts_above_neutral` | Y1 |
| 10 | `test_compute_ecm_to_fa_bias_active_orthogonal_target_suppresses_below_neutral` | Y1 |
| 11 | `test_compute_ecm_to_fa_bias_active_all_three_rates_same_multiplier_per_fa` | Y4 scaffolding contract |

### Zero-traction + boundary (3, Y7)

| # | Test | Y |
|---|---|---|
| 12 | `test_compute_ecm_to_fa_bias_active_zero_traction_returns_exact_neutral` (`multiplier == 1.0` exact; `alignment_score == 0.0` exact) | Y7 |
| 13 | `test_compute_ecm_to_fa_bias_active_zero_traction_count_in_diagnostics` | Y6 + Y7 |
| 14 | `test_compute_ecm_to_fa_bias_active_empty_adhesions_returns_empty_consistent_shape` (shape `(0, 3)`; `sampled_diagnostics` non-None) | Y5 |

### Diagnostics (3, Y6+Y13)

| # | Test | Y |
|---|---|---|
| 15 | `test_compute_ecm_to_fa_bias_active_diagnostics_dict_has_exact_10_keys` (assert `set(diagnostics_dict.keys()) == {n_adhesions, max_multiplier, min_multiplier, sampler_geometry, mechanism, k_active, max_alignment_score, min_alignment_score, zero_traction_count, score_bound}`) | Y6 |
| 16 | `test_compute_ecm_to_fa_bias_active_diagnostics_dict_no_loose_extra_keys` (subset check + length assertion; no surprise keys) | Y6 |
| 17 | `test_compute_ecm_to_fa_bias_active_diagnostics_sampler_geometry_matches_v1` (`sampler_geometry == "bilinear_cell_centered"` exact — sister with HB#4 v1) | Y13 |

### Failure-kind (1)

| # | Test | Y |
|---|---|---|
| 18 | `test_compute_ecm_to_fa_bias_active_multiplier_non_finite_failure_kind` (monkeypatch `sampled.orientation_tensor` with NaN/inf; assert `FAToECMBiasError` with `failure_kind == "active_multiplier_non_finite"`) | Y9 + Y11 |

### Validation order (1, Y12)

| # | Test | Y |
|---|---|---|
| 19 | `test_compute_ecm_to_fa_bias_active_invalid_k_active_raises_before_sampling` (empty adhesions + bad `k_active`; monkeypatch sampler to track invocation; assert raises `k_active_invalid` BEFORE sampler called) | Y12 |

### Effective-stiffness guard (1, Y8)

| # | Test | Y |
|---|---|---|
| 20 | `test_compute_ecm_to_fa_bias_active_does_not_reference_effective_stiffness` — function-scoped AST + string guard via `textwrap.dedent(inspect.getsource(compute_ecm_to_fa_bias_active))`; assert no `Attribute.attr` / `Name.id` reference to forbidden set `{effective_stiffness, stiffness_kpa, fiber_density, ligand_density, accumulated_traction_nNs_per_um2}` AND no forbidden literal substring | Y8 |

### Exports + private constant (1, Y14)

| # | Test | Y |
|---|---|---|
| 21 | `test_compute_ecm_to_fa_bias_active_only_function_exported` — `compute_ecm_to_fa_bias_active` importable from BOTH `acs.v2.dynamics` AND `acs.v2`; `_DEVIATORIC_SCORE_BOUND` NOT importable from either; `_validate_k_active` NOT importable | Y14 |

**Total: 21 tests**, all per lock §4. No additional tests beyond lock
contract; the wording-boundary discipline is structural in this unit
(test 20 AST guard + test 21 export discipline) — no separate
"forbidden literal in module docstring" meta-test required because
the active law's module docstring is the existing
`ecm_to_fa_bias.py` shared module docstring (HB#4 v1 grandfathered),
NOT a new module. Module-level forbidden-list paraphrasing is not
needed; the guard scope is the `compute_ecm_to_fa_bias_active`
function body specifically (Y8 function-scoped, NOT module-scoped).

---

## 6. Pre-commit batch (6 steps)

### 6.1 Schema grep — PASS

Verified all cited schema names exist at source:

- `ECMToFABiasResult` at `acs/v2/dynamics/ecm_to_fa_bias.py:127-141` ✓
- `ECMSampledAtFAs` at `acs/v2/dynamics/ecm_to_fa_bias.py:107-124` ✓
- `FAToECMBiasFailureKind` at `acs/v2/dynamics/ecm_to_fa_bias.py:89-92` (extending 2-kind to 4-kind literal) ✓
- `FAToECMBiasError` at `acs/v2/dynamics/ecm_to_fa_bias.py:95-104` ✓
- `RATE_NAMES` at `acs/v2/dynamics/ecm_to_fa_bias.py:78-87` ✓
- `sample_ecm_at_fa_positions` at `acs/v2/dynamics/ecm_to_fa_bias.py:244` ✓ (shared sampler, NOT new primitive)
- `compute_ecm_to_fa_bias_neutral` at `acs/v2/dynamics/ecm_to_fa_bias.py:352` ✓ (sister function for naming convention)

**Field names verified at source** (Y5 anti-pattern guard):
- `fa_ids`, `multipliers_per_fa`, `rate_names`, `sampled_diagnostics`, `diagnostics_dict` — all 5 present at `:137-141` ✓

### 6.2 git ls-files — PASS

All 6 cited sister-doc paths exist in git:

```
docs/v2_closed_loop_ecm_gate_phased_plan_locked.md      ✓
docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md  ✓
docs/v2_hard_blocker_1_2_constitutive_direction_locked.md ✓
docs/v2_hard_blocker_5_lyapunov_metric_locked.md        ✓
docs/v2_phase_e_composition_locked.md                   ✓
docs/v2_item_5_sweep_harness_locked.md                  ✓
docs/v2_focal_adhesion_dynamics_result_typed_locked.md  ✓ (B1 sister)
docs/v2_hard_blocker_4_active_locked.md                 ✓ (the lock)
docs/v2_hard_blocker_4_active_brief.md                  ✓ (opening brief)
```

### 6.3 Line-number citation regex — PASS

Code-line citations in this Sanity Gate doc:

- `acs/v2/dynamics/ecm_to_fa_bias.py:127-141` (ECMToFABiasResult) — verified above
- `acs/v2/dynamics/ecm_to_fa_bias.py:89-92` (FAToECMBiasFailureKind, current 2-kind literal extending to 4-kind) — verified above
- `acs/v2/dynamics/ecm_to_fa_bias.py:244` (sample_ecm_at_fa_positions) — verified above
- `acs/v2/dynamics/ecm_to_fa_bias.py:352` (compute_ecm_to_fa_bias_neutral) — verified above
- `acs/v2/ecm_substrate.py:36` (`_ORIENTATION_BOUND = 1.0` constant for `|T_ij| ≤ 1` invariant) — cited in §2.6 deviatoric bound derivation; required for Hard Rule 10 unit chain
- `acs/v2/ecm_substrate.py:115` (validator error message enforcing `|T_ij| ≤ 1`) — companion citation for the constant

All 6 line-number citations verifiable via `git ls-files` + `grep`. Lock §6 cited `acs/v2/ecm_substrate.py:42` as the invariant location, but `:42` is just `class ECMSubstrateState:` (the dataclass body offset, not the invariant); canonical lines for the bound are `:36` (constant) + `:115` (validator). This Sanity Gate cites both. (Codex `id=1618` correction.)

### 6.4 Deferred-vs-locked wording sweep — PASS

Wording discipline mirroring Phase E v1 + Item 5 sweep harness:

- **Locked (this unit)**: deviatoric Rayleigh score; bounded positive exponential map; all 3 rates same multiplier; `k_active` required; validation order ecm → k_active → sampler; `ECMToFABiasResult` grandfathered; private `_DEVIATORIC_SCORE_BOUND`; private `_validate_k_active`; failure-kind 2-kind → 4-kind extension; function-scoped AST + string guard
- **Deferred (separate later cycles, NOT in scope)**: Phase E v2 composition wrapper; per-rate selectivity (Phase E v3); `effective_stiffness` law (only if v3 consumes it); HB#4-active step 2
- **Out of scope (this lock)**: Phase E v2 satisfaction wording; raw Rayleigh `n.T @ T @ n` (Y1 anti-pattern); per-rate distinct multipliers (Y4); literature-fitted `k_active` default (Y3); `effective_stiffness` / `stiffness_kpa` / `fiber_density` / `ligand_density` / `accumulated_traction_nNs_per_um2` consumption (Y8 function-scoped guard)

Test names use `provides_active_bias_law_*` / `compute_ecm_to_fa_bias_active_*` patterns; **NO** `satisfies_phase_e_v2_*` / `satisfies_active_bias_*` / `satisfies_hb4_active_*` literals in any test name (sister-pattern with Phase E v1 + Item 5 wording-boundary discipline).

### 6.5 Absolute-claim consistency — PASS

§0 verdict table claims "all 6 Sanity Gate items PASS" — consistent
with §2.1-§2.6 individual PASS verdicts (one PASS per item). §1
function-level summary is consistent with §3 (Y5 schema-corrected
fields), §4 (Y9+Y11 failure-kind extension), §5 (test catalog
21-test count matches lock §4 21-test catalog), §6 (pre-commit batch
6/6 PASS).

**Forbidden literals check** (mirror of Item 5 test 24 wording-boundary
discipline; here at the doc level):

```
"Phase E v2 satisfied"           — 0 matches ✓
"satisfies_phase_e_v2"           — 0 matches ✓
"HB#4-active satisfies"          — 0 matches ✓
"satisfies_hb4_active"           — 0 matches ✓
"satisfies_active_bias"          — 0 matches ✓
```

(All grammatical "satisfied" / "satisfies" instances in this doc
refer to upstream locked items being satisfied as preconditions, e.g.
"phased plan §3 effective-stiffness guard satisfied by Q1=(b)
orientation-driven choice"; these are NOT claims about this unit
satisfying Phase E v2.)

### 6.6 Sister-gate-mirror — PASS

Step 6 sister-gate-mirror at 4 layers (per memory rule
`design_note_pre_commit_batch.md` Step 6 promoted; sister-pattern
with HB#4 v1 / HB#5 / Phase E v1 / Item 5 Sanity Gate cycles):

| Layer | HB#4 v1 (sister) | HB#4-active (this unit) | Status |
|---|---|---|---|
| **Code** | `compute_ecm_to_fa_bias_neutral`, validation chain `ecm.validate() → sampler` | `compute_ecm_to_fa_bias_active`, validation chain `ecm.validate() → _validate_k_active → sampler` | parallel structure with cheap-parameter-failure-first extension (Y12) |
| **Design** | All-1.0 multipliers; Phase D no-op default | Non-1.0 multipliers via deviatoric Rayleigh; Phase E v2 step 1 scaffolding | parallel design surface (same return type, different multiplier law) |
| **API** | Public `compute_ecm_to_fa_bias_neutral`; shared `sample_ecm_at_fa_positions` | Public `compute_ecm_to_fa_bias_active`; shared sampler (NO new primitive); private `_DEVIATORIC_SCORE_BOUND` + `_validate_k_active` (Y14) | API symmetric: 1 public per variant; private helpers properly hidden |
| **Failure-kind** | 2 kinds: `fa_bias_position_outside_ecm_grid`, `non_finite_fa_position` | Extended to 4 kinds: + `k_active_invalid`, + `active_multiplier_non_finite` (Y9) | Literal extension (additive); existing kinds unchanged (sampler errors still surface as before) |

**Wording-boundary mirror** (from Phase E v1 / Item 5 / HB#5 wording-boundary discipline):
- Phase E v1: "ECM-side composition evidence" wording boundary (no full-closure satisfaction claim)
- Item 5: `provides_*_evidence` test names; `1e-12` IEEE roundoff NOT gate threshold
- HB#5: `provides_lyapunov_metric_*_evidence` test names; metric reports active-cells-only
- **HB#4-active step 1** (this unit): `compute_ecm_to_fa_bias_active_*` test names (NO `satisfies_phase_e_v2_*`); function provides active bias law / evidence channel, NOT Phase E v2 satisfaction (test 20 function-scoped AST guard enforces structural separation; test 21 export discipline enforces single-symbol surface)

**IEEE-roundoff distinction** (Item 5 sister-pattern):
- Test 6 isotropic-IC neutral uses **exact equality** (no `rel=`) since `T - 0.5·trace(T)·I = 0` is exact arithmetic when `T = a·I` (subtraction of equal float64 values yields 0.0 exact)
- Test 7 saturation example uses `rel=1e-10` for IEEE roundoff allowance on `np.sqrt(2.0)` and `cos(π/8)` / `sin(π/8)` evaluation; this `1e-10` is **explicitly NOT a gate threshold** — it's the standard `pytest.approx` IEEE roundoff allowance for trigonometric/square-root evaluation precision
- **NO new tolerance constant introduced** in module source (Y14); all tolerances are test-side `rel=` arguments, not validation-code constants

**Sister-gate-mirror divergence (deliberate)**:

| Aspect | HB#4 v1 (sister) | HB#4-active (this unit) | Reason |
|---|---|---|---|
| Failure kinds added | 2 baseline | +2 (k_active_invalid + active_multiplier_non_finite) | active law has parameter discipline (k_active) and finite-guards (Y11) that v1 didn't need |
| Private constants | none | `_DEVIATORIC_SCORE_BOUND` | bound is derived (√2), not a fitted constant; private + reported as value not name (Y14) |
| Validation chain length | 2 steps | 3 steps (cheap-parameter-failure-first per Y12) | k_active discipline must fail before sampler |

---

## 7. Outstanding-before-code

The following items are **not** done yet but **must** be in the code commit:

1. **Module-level imports**: `Final` from typing (for `_DEVIATORIC_SCORE_BOUND`); `inspect` and `textwrap` and `ast` are already used by HB#5 / Item 5 sister tests so available; numpy already imported
2. **Public function** `compute_ecm_to_fa_bias_active` per lock §1 body (with all Y annotations preserved as code comments where they affect readability; full Y1-Y15 trace in lock not duplicated in code per CLAUDE.md "Default to writing no comments" rule — only WHY-non-obvious comments)
3. **Private helper** `_validate_k_active` per lock §1 body; raises `FAToECMBiasError` with `failure_kind="k_active_invalid"`
4. **Private constant** `_DEVIATORIC_SCORE_BOUND: Final[float] = float(np.sqrt(2.0))` — module-level
5. **Failure-kind literal extension** in `FAToECMBiasFailureKind` at line 89-92 (additive: 2 → 4 kinds)
6. **Exports**: 1 new symbol `compute_ecm_to_fa_bias_active` through both `acs/v2/dynamics/__init__.py` and `acs/v2/__init__.py` `__all__` (alphabetical reorder); private constant + private helper NOT exported (test 21)
7. **Test file**: NEW `tests/test_v2_ecm_to_fa_bias_active.py` (21 tests) per lock §4. Decision: **new file**, NOT extending existing `tests/test_v2_ecm_to_fa_bias.py` (536 lines) — separation by variant (neutral vs active) preserves test-file purity for sister-gate combined regression
8. **AST + string guard test** (test 20): use `inspect.getsource(compute_ecm_to_fa_bias_active)` + `textwrap.dedent` + `ast.walk` to enumerate `Attribute` and `Name` nodes; assert no `attr ∈ forbidden_set` and no `id ∈ forbidden_set`; also string-search `inspect.getsource(...)` for forbidden substrings (HB#5 Y11 sister-pattern)
9. **Validation-order monkeypatch test** (test 19): use `pytest.MonkeyPatch` to replace `acs.v2.dynamics.ecm_to_fa_bias.sample_ecm_at_fa_positions` with a tracking sentinel; pass empty adhesions + bad `k_active`; assert raises `k_active_invalid` AND sentinel was never invoked
10. **Verification before commit**:
    - `pytest tests/test_v2_ecm_to_fa_bias_active.py -x -v` (21/21 PASS)
    - Sister-gate combined `pytest tests/test_v2_ecm_to_fa_bias.py tests/test_v2_ecm_to_fa_bias_active.py tests/test_v2_closed_loop_phase_e.py tests/test_v2_closed_loop_phase_e_sweep.py` (no regression)
    - Full v2 suite `pytest tests/test_v2_*.py` (560 + 21 = 581/581 PASS expected)

---

## 8. Verdict table

| Item | Verdict | Section | Lock anchor |
|---|---|---|---|
| 1. Units | PASS | §2.1 | Y1, Y2, Y3 |
| 2. Boundary | PASS (11 of 21 tests) | §2.2 | Y3, Y6, Y7, Y14 |
| 3. Conservation | PASS (pure read-only) | §2.3 | Y15 grandfathered |
| 4. Numerical | PASS | §2.4 | Y7, Y11 |
| 5. Sign | PASS | §2.5 | Y1 |
| 6. Measurement-protocol | PASS (Hard Rule 11 central anchor; deviatoric bound derivation in §2.6) | §2.6 | Y1, Y10 |
| Schema-corrected fields (Y5) | PASS (verified at `:127-141`) | §3 | Y5 |
| Failure-kind extension (Y9+Y11) | PASS (additive 2→4 kinds; dual finite guards) | §4 | Y9, Y11 |
| Test catalog mirror (21 tests) | PASS (mirrors lock §4) | §5 | Y1–Y15 |
| Pre-commit batch 6 steps | PASS | §6 | — |
| Sister-gate-mirror (Step 6) | PASS at 4 layers (code/design/API/failure-kind) + wording-boundary | §6.6 | — |

**Overall verdict**: all 6 Sanity Gate items + 5 cross-reference items
PASS. Ready for impl-work Codex review per (B-trigger) precedent.

---

## 9. References

- Lock: `docs/v2_hard_blocker_4_active_locked.md` (commit `fed370b`)
- Brief: `docs/v2_hard_blocker_4_active_brief.md` (commit `013e161`)
- Parent locked plan: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- Sister locks (per lock §6):
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` (HB#4 v1
    neutral; same return type, shared sampler)
  - `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md` (HB#1+#2;
    Y4 branch-defined zero-traction sister-pattern; orientation-only update)
  - `docs/v2_phase_e_composition_locked.md` (Phase E v1; uses HB#4 v1
    neutral; v2 composition is separate)
  - `docs/v2_hard_blocker_5_lyapunov_metric_locked.md` (Y11 AST-walk
    function-scoped guard sister-pattern)
  - `docs/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1 typed
    schema discipline; here grandfathered for `ECMToFABiasResult`)
  - `docs/v2_item_5_sweep_harness_locked.md` (Item 5 sweep harness;
    uses `0.5·I` IC which is automatically neutral under this active
    law per Y1 deviatoric form; sister Sanity Gate at
    `docs/v2_item_5_sweep_harness_sanity_gate.md`)
- Verified at source for Y5 (ECMToFABiasResult schema):
  `acs/v2/dynamics/ecm_to_fa_bias.py:127-141`
- Verified at source for Y9 (FAToECMBiasFailureKind):
  `acs/v2/dynamics/ecm_to_fa_bias.py:89-92`
- Verified at source for shared sampler (HB#4 v1 grandfathered):
  `acs/v2/dynamics/ecm_to_fa_bias.py:244` (`sample_ecm_at_fa_positions`)
- Verified at source for orientation_tensor schema bound `|T_ij| ≤ 1`:
  `acs/v2/ecm_substrate.py:36` (`_ORIENTATION_BOUND = 1.0` constant) + `:115` (validator error message). Lock §6 cited `:42` which is the dataclass body offset; canonical lines for the bound invariant are `:36` (constant) + `:115` (validator).
- Adversarial debate posture: memory `feedback_aggressive_design_debate.md`,
  PI id=809
- Hard Rule 10 inline derivation: memory
  `rule10_unit_derivation_in_docs.md`, Codex id=1234 (deviatoric bound
  derivation §2.6 mandatory inline)
- Hard Rule 11 wording-boundary meta-test: memory
  `hard_rule_11_wording_boundary_meta_test.md` (Y1 deviatoric form is
  the **fourth** Phase E composition Hard Rule 11 catch family — Phase E
  v1 Y1, HB#5 Y1, Item 5 Y1+Y4+Y17, **HB#4-active Y1**)
- Step 6 pre-commit batch: memory `design_note_pre_commit_batch.md`,
  Codex id=1428 (Step 6 sister-gate-mirror promoted)
- Design-discussion adversarial cycle: MCP rounds 1–3 + SEAL ack
  (id=1605–1612 in design-discussion room; impl-work dispatch id=1615)
