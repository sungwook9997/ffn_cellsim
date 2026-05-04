# V2 Closed-Loop ECM Gate — HB#4-Active Step 1 (Deviatoric-Rayleigh Orientation Bias) — Locked Design

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (3-round adversarial lock,
PI id=809 aggressive debate posture applied, PI id=1008/1057/1161
autonomy + visible deliverable focus)
**Source unit**: design-discussion `topic=v2-hard-blocker-4-active`,
MCP id 1605–1612 (rounds 1–3 + seal ack)
**Brief**: `docs/v2_hard_blocker_4_active_brief.md` (commit `013e161`)
**Parent locked plan**: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
**Upstream sister locks**:
- `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` (HB#4 v1
  neutral; same `ECMToFABiasResult` return type + shared sampler)
- `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md` (HB#1+#2
  orientation-only update; Y4 branch-defined zero-traction sister-pattern)
- `docs/v2_phase_e_composition_locked.md` (Phase E v1 composition
  using HB#4 v1 neutral; v2 composition is a separate later cycle)
- `docs/v2_hard_blocker_5_lyapunov_metric_locked.md` (HB#5 sister:
  Y11 AST-walk source-guard pattern + Y9/Y13 effective_stiffness
  exclusion sister-pattern)
- `docs/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1 typed
  schema discipline; here grandfathered for `ECMToFABiasResult`)
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for the HB#4-active step 1 Sanity Gate doc + code entry.

---

## 0. Scope

This unit locks **HB#4-active step 1** — the first orientation-driven
non-neutral active bias law for `compute_ecm_to_fa_bias_active`. It
provides the active-bias **scaffolding** for Phase E v2, but is **NOT
itself Phase E v2 satisfaction**.

### v1 (neutral) → active step 1 → step 2 (composition) split

- **HB#4 v1 (locked + impl)**: `compute_ecm_to_fa_bias_neutral` returns
  all-1.0 multipliers; consumed by Phase E v1 (also locked + impl).
- **HB#4-active step 1 (this lock)**:
  `compute_ecm_to_fa_bias_active(adhesions, ecm, *, k_active)` returns
  non-1.0 multipliers via deviatoric-Rayleigh orientation bias. Same
  `ECMToFABiasResult` return type as v1 (forward-compat). Pure
  read-only.
- **HB#4-active step 2 (separate later cycle)**:
  `step_closed_loop_phase_e_v2` composition wrapper that swaps v1
  neutral for active step 1 in the Phase E composition. Out of scope
  here.

### Wording discipline (Y4 — central anchor)

This lock follows Phase E v1 + Item 5 sweep harness sister-pattern
"evidence not satisfaction" discipline:

- Test names use `provides_active_bias_law_*` / `provides_*_evidence`
  pattern, **NOT** `satisfies_*`
- "All-three-rates same multiplier" is honestly framed as **active-bias
  scaffolding** — closed-loop machinery demonstration, NOT per-rate
  selectivity
- Per-rate selective biology (assembly vs disassembly differential
  favoring mature stable FAs) is **Phase E v3+ scope** with its own
  per-rate literature anchors

### What this unit IS

- Pure read-only function `compute_ecm_to_fa_bias_active(adhesions, ecm, *, k_active)`
- Returns `ECMToFABiasResult` (HB#4 v1 grandfathered, NO new result type)
- Uses shared `sample_ecm_at_fa_positions` sampler (NO new sampling
  primitive)
- Bias law: deviatoric Rayleigh quotient → bounded positive exponential
  multiplier
- All 3 rates (`k_maturity_per_s`, `k_bind_per_s`, `k_unbind_per_s`)
  receive same multiplier per FA (scaffolding, not per-rate selective)

### What this unit is NOT

- **NOT a Phase E v2 satisfaction step**: composition wrapper is the
  separate step 2 cycle
- **NOT consumer of `effective_stiffness`, `stiffness_kpa`,
  `fiber_density`, `ligand_density`, or `accumulated_traction_nNs_per_um2`**:
  function-scoped AST + string guard enforces this (Y8); ONLY consumes
  `sampled.orientation_tensor`
- **NOT per-rate selective**: all 3 rates receive same multiplier
  (Y4 honest scaffolding wording)
- **NOT a literature-fitted parameter unit**: `k_active` is required
  (no default); tests use values as mathematical fixtures, NOT
  biological defaults (Y3)

---

## 1. Final Lock Summary

### Constants

```python
_DEVIATORIC_SCORE_BOUND: Final[float] = float(np.sqrt(2.0))   # private
# = √2 ≈ 1.4142, schema-tight bound on |n^T Q n| where Q is deviatoric of T
# Derivation: see §3 measurement-protocol section
_RATE_COUNT: Final[int] = len(RATE_NAMES)  # = 3
```

`_DEVIATORIC_SCORE_BOUND` is **private** — NOT exported (Y14). Its
value is reported in `diagnostics_dict["score_bound"]` as a value, NOT
the constant name.

### Function signature + body

```python
def compute_ecm_to_fa_bias_active(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
    *,
    k_active: float,                 # REQUIRED, no default (Y3)
) -> ECMToFABiasResult:
    """HB#4-active step 1: orientation-driven Rayleigh-quotient bias on
    deviatoric tensor.

    For each FA i:
      n_fa = traction_force / |traction_force|  (branch-defined zero-traction → exact neutral)
      T_local = ECM orientation tensor sampled at FA position (shared sampler)
      Q = T_local - 0.5 · trace(T_local) · I  (deviatoric: alignment-only, isotropic factored out)
      score = n_fa.T @ Q @ n_fa  (dimensionless; |score| ≤ √2 schema bound)
      multiplier = exp(K_active · score)  (bounded positive exponential map; all 3 rates same)

    Range: multiplier ∈ [exp(-√2·K_active), exp(+√2·K_active)]; ≡ 1.0 at score=0
    (i.e., isotropic T = a·I → Q=0 → score=0 → neutral; Item 5 IC `T = 0.5·I`
    is automatically neutral under this law).
    """
    # Validation order Y12: ecm → k_active → sampler (cheap parameter failure first)
    ecm.validate()
    _validate_k_active(k_active)
    sampled = sample_ecm_at_fa_positions(adhesions, ecm)

    n_fa = len(adhesions)
    multipliers = np.ones((n_fa, _RATE_COUNT), dtype=np.float64)
    alignment_scores = np.zeros(n_fa, dtype=np.float64)
    zero_traction_count = 0

    for i, fa in enumerate(adhesions):
        traction = np.array(fa.traction_force_nN_xy, dtype=np.float64)
        traction_norm = float(np.linalg.norm(traction))
        if traction_norm > 0.0:                           # branch-defined, no epsilon (HB#1+#2 Y4 sister)
            n_fa_unit = traction / traction_norm
            T_local = sampled.orientation_tensor[i]       # (2, 2) sampled at FA position
            trace_T = float(T_local[0, 0] + T_local[1, 1])
            Q = T_local - 0.5 * trace_T * np.eye(2, dtype=np.float64)
            alignment_score = float(n_fa_unit @ Q @ n_fa_unit)
        else:
            alignment_score = 0.0                          # exact (Y7)
            zero_traction_count += 1

        # Defensive finite guards (Y11): catches corrupted inputs (defense in depth)
        if not np.isfinite(alignment_score):
            raise FAToECMBiasError(
                "active_multiplier_non_finite",
                f"FA {fa.adhesion_id}: alignment_score {alignment_score} non-finite "
                f"(suggests corrupted orientation_tensor; ecm.validate() should have caught)",
            )
        multiplier = float(np.exp(k_active * alignment_score))
        if not np.isfinite(multiplier):
            raise FAToECMBiasError(
                "active_multiplier_non_finite",
                f"FA {fa.adhesion_id}: multiplier {multiplier} non-finite "
                f"(k_active={k_active}, score={alignment_score})",
            )

        multipliers[i, :] = multiplier                     # all 3 rates same (Y4 Q2=a)
        alignment_scores[i] = alignment_score

    # Diagnostics dict (Y6 exact 10 keys; Y13 sampler_geometry sister-consistent)
    diagnostics_dict = {
        "n_adhesions": int(n_fa),
        "max_multiplier": float(multipliers.max()) if n_fa > 0 else 1.0,
        "min_multiplier": float(multipliers.min()) if n_fa > 0 else 1.0,
        "sampler_geometry": "bilinear_cell_centered",       # Y13: HB#4 v1 sister
        "mechanism": "deviatoric_rayleigh_orientation",
        "k_active": float(k_active),
        "max_alignment_score": float(alignment_scores.max()) if n_fa > 0 else 0.0,
        "min_alignment_score": float(alignment_scores.min()) if n_fa > 0 else 0.0,
        "zero_traction_count": int(zero_traction_count),
        "score_bound": float(_DEVIATORIC_SCORE_BOUND),      # Y14: value, not name
    }

    # Result construction (Y5 schema-corrected: fa_ids, multipliers_per_fa, rate_names,
    # sampled_diagnostics, diagnostics_dict — verified at ecm_to_fa_bias.py:128)
    return ECMToFABiasResult(
        fa_ids=sampled.fa_ids,
        multipliers_per_fa=multipliers,
        rate_names=tuple(RATE_NAMES),
        sampled_diagnostics=sampled,
        diagnostics_dict=diagnostics_dict,
    )
```

### Validation helper

```python
def _validate_k_active(k_active: float) -> None:
    """k_active must be float (not bool — Python bool ⊂ int trap),
    finite, > 0, and produce a finite derived multiplier bound."""
    if isinstance(k_active, bool):
        raise FAToECMBiasError(
            "k_active_invalid",
            f"k_active must be float (not bool — Python bool ⊂ int trap), "
            f"got {type(k_active).__name__}",
        )
    if not (np.isfinite(k_active) and k_active > 0.0):
        raise FAToECMBiasError(
            "k_active_invalid",
            f"k_active must be finite positive, got {k_active}",
        )
    derived_bound = np.exp(k_active * _DEVIATORIC_SCORE_BOUND)
    if not np.isfinite(derived_bound):
        raise FAToECMBiasError(
            "k_active_invalid",
            f"k_active={k_active} produces non-finite multiplier bound "
            f"exp({k_active}·√2) = {derived_bound}; reduce k_active",
        )
```

### Failure-kind literal extension (Y9)

```python
# acs/v2/dynamics/ecm_to_fa_bias.py
FAToECMBiasFailureKind = Literal[
    "fa_bias_position_outside_ecm_grid",
    "non_finite_fa_position",
    "k_active_invalid",                    # NEW (Y9, HB#4-active step 1)
    "active_multiplier_non_finite",        # NEW (Y9, HB#4-active step 1)
]
```

### Forbidden in HB#4-active step 1

- Raw Rayleigh quotient `n^T T n` (Y1 — conflates isotropic with
  alignment; Item 5 IC `T=0.5·I` would falsely boost all FAs)
- Per-rate distinct multipliers (Y4 — scaffolding only; per-rate
  selectivity is Phase E v3+ scope)
- "Asymptotic-to-bound saturation" wording (Y2 — this is bounded
  exponential map over a schema-bounded score, not asymptotic
  saturation as t→∞ in HB#1+#2 sense)
- Reading or referencing `effective_stiffness`, `stiffness_kpa`,
  `fiber_density`, `ligand_density`, or `accumulated_traction_nNs_per_um2`
  in the active function body or docstring (Y8 — function-scoped AST +
  string guard via meta-test 18; ONLY `sampled.orientation_tensor` is
  consumed)
- Default value for `k_active` (Y3 — caller must supply explicit value
  with external provenance note before execution)
- Mutable input handling (function is pure; uses HB#4 v1 grandfathered
  result type which already enforces frozen+slots immutability)
- Use of epsilon in zero-traction branch (Y7 — exact `traction_norm > 0.0`
  branch; HB#1+#2 Y4 sister-pattern)
- Phase E v2 composition wrapper, hooks, strategy, or any silent
  composition path inside this function (Y3 — composition is a
  separate cycle)
- Adding new result dataclass (Y15 — `ECMToFABiasResult` v1
  grandfathered, active-specific metadata in `diagnostics_dict`)
- Loose extra keys in `diagnostics_dict` beyond the locked 10 keys
  (Y6 — tests assert exact key set)
- Public export of `_DEVIATORIC_SCORE_BOUND` or `_validate_k_active`
  (Y14 — private only; only `compute_ecm_to_fa_bias_active` is
  exported)
- Sampler geometry string other than `"bilinear_cell_centered"` (Y13 —
  HB#4 v1 sister-consistency)
- Validation order other than `ecm.validate() → _validate_k_active()
  → sample_ecm_at_fa_positions()` (Y12)

---

## 2. Reasoned-acceptance trace (Y1–Y15)

The lock converged after 3 rounds. Each Y is a Claude concession with
the round in which it was accepted, preserving the adversarial-debate
audit trail (PI id=809 posture).

**Y1 (round 2, accepted Codex C1 = deviatoric Rayleigh score)**:
- Claude opening lean (round 1): raw Rayleigh quotient
  `bias = n_fa^T T_local n_fa`.
- Codex catch: raw form conflates isotropic orientation magnitude with
  directional alignment. Item 5 IC `T = 0.5·I` gives `n^T T n = 0.5`
  for every traction direction → all FAs uniformly boosted, NO
  alignment-with-traction discrimination. Anti-pattern.
- Resolution: deviatoric form `Q = T - 0.5·trace(T)·I`,
  `alignment_score = n_fa^T Q n_fa`. Properties: isotropic `T = a·I` →
  `Q = 0` → score = 0 → multiplier = 1.0 (neutral); rank-1 aligned →
  `+0.5`; rank-1 orthogonal → `-0.5`; schema-valid `|T_ij| ≤ 1` →
  `|score| ≤ √2`.

**Y2 (round 2, accepted Codex C2 = "bounded positive exponential map" wording)**:
- Brief used "asymptotic-to-bound exponential" sister-claim with
  HB#1+#2 Y10 `1 - exp(-x)` saturation form.
- Codex catch: `exp(k · score)` over bounded score is **bounded
  exponential map over a schema-bounded input**, NOT saturation in the
  HB#1+#2 sense (which approaches a limit as t→∞).
- Resolution: §1 + §3 use the precise wording. Bounds derive from
  schema: `multiplier ∈ [exp(-√2·k_active), exp(+√2·k_active)]`.

**Y3 (round 2, accepted Codex C3 = no-default k_active + finite-bound
check + bool guard)**:
- Codex agreed with K-ii (no default) and added robustness:
  - Validate `k_active > 0`, finite, not bool (Python `bool ⊂ int` trap)
  - Validate derived bound finite: `np.isfinite(np.exp(k_active · √2))`;
    raise `k_active_invalid` if not
  - Tests use values as mathematical fixtures, NOT biological defaults
- Resolution: `_validate_k_active` helper + 4 validation tests (5 with
  too-large k_active overflow case).

**Y4 (round 2, accepted Codex C4 = all-three-same as scaffolding ONLY)**:
- Claude proposal already framed as "minimum viable active-bias
  scaffolding"; Codex confirmed and asked for explicit lock wording.
- Resolution: §0 + §1 wording: "active-bias scaffolding;
  per-rate selectivity is Phase E v3+ scope".

**Y5 (round 2, accepted Codex C5 = corrected ECMToFABiasResult fields)**:
- Claude round 1 pseudocode used `sampled` and `diagnostics` — neither
  exists in the actual schema.
- Verified at `acs/v2/dynamics/ecm_to_fa_bias.py:128`: actual fields
  are `fa_ids`, `multipliers_per_fa`, `rate_names`,
  `sampled_diagnostics: Optional[ECMSampledAtFAs]`,
  `diagnostics_dict: dict[str, float | int | str]`.
- Resolution: corrected result construction. Empty FA list still calls
  `sample_ecm_at_fa_positions((), ecm)` so `sampled_diagnostics` is
  non-None and shape-consistent (HB#4 v1 sister).
- Sister-pattern with Phase E v1 B2 catch (`multipliers_per_fa` not
  `multipliers`).

**Y6 (round 2, accepted Codex C6 = exact 10-key diagnostics_dict)**:
- Codex specified exact key set; tests assert no loose extra keys.
- Resolution: 10 keys: `n_adhesions`, `max_multiplier`,
  `min_multiplier`, `sampler_geometry`, `mechanism`, `k_active`,
  `max_alignment_score`, `min_alignment_score`, `zero_traction_count`,
  `score_bound`. Test 13 asserts exact set.

**Y7 (round 2, accepted Codex C7 = zero-traction exact-neutral with count)**:
- For `traction_norm == 0`: `alignment_score = 0.0` exact, `multiplier
  = 1.0` exact. No direction sampling. Count in `zero_traction_count`
  diagnostic.

**Y8 (round 2, accepted Codex C8 = function-scoped guard, NOT
module-wide)**:
- `ecm_to_fa_bias.py` already contains policy text mentioning
  `effective_stiffness` (HB#4 v1 docstrings/comments).
- Module-wide string check would false-positive.
- Resolution: function-scoped AST + string guard via
  `textwrap.dedent(inspect.getsource(compute_ecm_to_fa_bias_active))`
  (HB#5 Y11 sister-pattern). Forbidden set:
  `effective_stiffness`, `stiffness_kpa`, `fiber_density`,
  `ligand_density`, `accumulated_traction_nNs_per_um2`. ONLY
  `sampled.orientation_tensor` may be consumed.

**Y9 (round 2, accepted Codex C9 = failure-kind literal extension)**:
- Current `FAToECMBiasFailureKind` covers only sampler geometry
  failures.
- Resolution: extend with `k_active_invalid` and
  `active_multiplier_non_finite`; sampler failures unchanged.

**Y10 (round 3, accepted Codex A1 = schema-tight √2 test matrix
math fix)**:
- Claude round 2 test matrix `T = [[1, -1], [-1, 1]]` had **trace 2**,
  NOT 0. So `Q = T - I = [[0, -1], [-1, 0]]` has eigenvalues `±1`,
  NOT `±√2`. My math was wrong.
- Codex correction: `T = [[1, 1], [1, -1]]` is **already traceless**
  (trace 0), so `Q = T`. Eigenvalues are `±√2` (from `λ² - 2 = 0`).
  Eigenvector for `+√2` at angle `π/8` (since `tan(2θ) = b/a = 1/1 =
  1` → `2θ = π/4` → `θ = π/8`).
- Resolution: corrected test matrix + bound derivation example in §3.

**Y11 (round 3, accepted Codex A2 = finite-score defensive guard)**:
- Defensive check `np.isfinite(alignment_score)` BEFORE `np.exp(...)`
  catches corrupted orientation_tensor (defense in depth, even though
  `ecm.validate()` should prevent this).
- Resolution: dual finite guards on score and multiplier; both raise
  `active_multiplier_non_finite` failure_kind (no new kind needed).

**Y12 (round 3, accepted Codex A3 = validation order ecm → k_active →
sampler)**:
- Cheap parameter failure first (`_validate_k_active` independent of
  FA list); ECM boundary validation at public entry; sampler last
  (most expensive + may also raise on FA position).
- Resolution: explicit order in §1; test 17 asserts invalid `k_active`
  raises before sampler is called (e.g., empty adhesions + bad
  k_active still raises).

**Y13 (round 3, accepted Codex A4 = sampler_geometry sister-consistent
string)**:
- Claude round 2 used `"bilinear"`; HB#4 v1 uses
  `"bilinear_cell_centered"`.
- Resolution: `"bilinear_cell_centered"` — exact v1 string;
  test 14 asserts sister consistency.

**Y14 (round 3, accepted Codex A5 = single export + private constant)**:
- Public export: `compute_ecm_to_fa_bias_active` only (1 new symbol)
  through both `acs.v2.dynamics` and `acs.v2`.
- `_DEVIATORIC_SCORE_BOUND` and `_validate_k_active` private (NOT
  exported).
- `score_bound` value in `diagnostics_dict` is the value, NOT the
  constant name reference.
- Test 19 asserts `_DEVIATORIC_SCORE_BOUND` NOT importable from public
  surfaces.

**Y15 (round 3, accepted Codex A6 = ECMToFABiasResult reuse, no new
result dataclass)**:
- Active-specific metadata lives in the `diagnostics_dict` field
  because HB#4 v1 grandfathered that surface.
- Forward-compat: future `step_closed_loop_phase_e_v2` can swap
  `compute_ecm_to_fa_bias_neutral` for `compute_ecm_to_fa_bias_active`
  (both return same `ECMToFABiasResult` type).

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

1. **Units**:
   - `traction_force_nN_xy`: nN
   - `traction_norm`: nN
   - `n_fa_unit`: dimensionless unit vector
   - `T_local`, `Q`, `alignment_score`: dimensionless
   - `k_active`: dimensionless (caller-supplied)
   - `multiplier`: dimensionless (≥ 0, multiplicative on FA rates)
   - `_DEVIATORIC_SCORE_BOUND`: dimensionless = √2 ≈ 1.4142
2. **Boundary**:
   - Empty adhesions: `multipliers.shape == (0, 3)`; sampler still
     called → `sampled_diagnostics` non-None; `diagnostics_dict`
     populated with default zero counts; `max/min_multiplier` default
     to 1.0
   - All-zero traction: `alignment_score = 0.0` exact, `multiplier =
     1.0` exact, `zero_traction_count == n_fa`
   - Isotropic ECM `T = a·I` (any `a`): `Q = 0` → `alignment_score = 0`
     → `multiplier = 1.0` exactly for ALL FAs (Y1 anti-pattern guard;
     Item 5 IC automatic neutral)
   - Schema-saturating ECM `T = [[1, 1], [1, -1]]` traceless rank-2:
     score reaches `±√2` at angle `π/8`
   - `k_active = ε` small: multiplier ≈ 1.0 + k·score (Taylor)
   - `k_active = 1000` (overflow): `_validate_k_active` raises
     `k_active_invalid` because `exp(1000·√2) = exp(1414.2)` overflows
     float64. Threshold derivation: `exp(k·√2) > float64.max` iff
     `k > log(float64.max) / √2 ≈ 709.78 / 1.4142 ≈ 501.892`. Use `1000`
     to be safely above threshold. (Codex `id=1618` correction: prior
     `k_active = 100` claim was false — `exp(100·√2) ≈ 2.62e61` is
     finite, well below `float64.max ≈ 1.79e308`.)
3. **Conservation**: function is pure read-only — no state mutation.
   The HB#4 v1 grandfathered result `ECMToFABiasResult` is frozen+slots
   (immutable). `multipliers_per_fa` is freshly allocated.
4. **Numerical**:
   - `np.float64` enforced throughout
   - Branch `traction_norm > 0.0` exact (no epsilon)
   - Defensive finite guards on `alignment_score` and `multiplier`
     (Y11) catch corrupted inputs without propagating NaN/inf
5. **Sign**:
   - `alignment_score ∈ [-√2, +√2]`; positive when ECM aligned with FA
     traction direction (boost), negative when orthogonal (suppress)
   - `multiplier > 0` always (exponential of any finite real)
   - `multiplier > 1` for `score > 0` (boost); `multiplier < 1` for
     `score < 0` (suppress); `multiplier = 1` for `score = 0` (neutral)
6. **Measurement-protocol** (Hard Rule 11, the central anchor):
   - **Deviatoric score bound derivation** (in lock body, mandatory
     per CLAUDE.md `rule10_unit_derivation_in_docs` memory):
     - `T` is symmetric 2×2 with schema invariant `|T_ij| ≤ 1` (no PSD
       constraint per `acs/v2/ecm_substrate.py:42`)
     - `trace_T = T_xx + T_yy ∈ [-2, +2]`
     - `Q = T - 0.5·trace(T)·I` is **traceless symmetric**:
       - `Q_xx = T_xx - 0.5·(T_xx + T_yy) = 0.5·(T_xx - T_yy)`
       - `Q_yy = T_yy - 0.5·(T_xx + T_yy) = -0.5·(T_xx - T_yy) = -Q_xx`
       - `Q_xy = T_xy`
     - For unit `n = (cos θ, sin θ)`:
       - `n^T Q n = 0.5·(T_xx - T_yy)·cos(2θ) + T_xy·sin(2θ)`
       - `|n^T Q n| ≤ √(((T_xx - T_yy)/2)² + T_xy²)` (Cauchy-Schwarz)
     - Schema bounds: `|T_ij| ≤ 1` → `|T_xx - T_yy| ≤ 2` →
       `((T_xx - T_yy)/2)² ≤ 1`; `T_xy² ≤ 1`. Sum: `≤ 2` → `|n^T Q n|
       ≤ √2`.
     - **Saturation example** (Y10 corrected): `T = [[1, 1], [1, -1]]`
       (already traceless), `Q = T`, eigenvalues `±√2`. At eigenvector
       angle `θ = π/8` (since `tan(2θ) = 1/1`), `n^T Q n = +√2`
       (saturates the bound).
   - **Anti-pattern guard (Y1)**: raw Rayleigh `n^T T n` would give
     `0.5` for `T = 0.5·I` regardless of `n`, falsely boosting all FAs.
     Deviatoric form correctly produces `score = 0` for any isotropic
     `T`. Test 5 enforces this.
   - **Closure mechanism**: HB#1+#2 evolves `T` toward rank-1 `n_traction
     ⊗ n_traction`; `T` becomes deviatoric-aligned with `n_traction`;
     active law sees increasing `score` over time → multiplier grows
     above neutral → FA rates rescaled.
   - **Item 5 IC compatibility**: Item 5 sweep harness uses
     `0.5·np.eye(2)` as initial orientation; under deviatoric form,
     this is automatically neutral (no spurious initial boost). Item 5
     sweep harness can therefore use HB#4-active without polluting the
     "scaffolding starts neutral" baseline.

---

## 4. Test Catalog (21 tests)

### Validation (5)

1. `test_compute_ecm_to_fa_bias_active_k_active_zero_or_negative_raises`
2. `test_compute_ecm_to_fa_bias_active_k_active_nonfinite_raises`
3. `test_compute_ecm_to_fa_bias_active_k_active_bool_raises` (Python
   `bool ⊂ int` trap)
4. `test_compute_ecm_to_fa_bias_active_k_active_too_large_raises_non_finite_bound`
5. `test_compute_ecm_to_fa_bias_active_invalid_ecm_raises_at_validate`

### Deviatoric form (3, Y1+Y10 corrected math)

6. `test_compute_ecm_to_fa_bias_active_isotropic_orientation_returns_neutral_multiplier`
   — `T = 0.5·I` (Item 5 IC); assert `multiplier == 1.0` exactly for
   any traction direction (Y1 anti-pattern guard)
7. `test_compute_ecm_to_fa_bias_active_schema_tight_bound_score_approaches_sqrt2`
   — `T = [[1, 1], [1, -1]]` (traceless), `n = (cos(π/8), sin(π/8))`;
   assert `multiplier == exp(√2 · k_active)` within `rel=1e-10`;
   assert `diagnostics_dict["max_alignment_score"] ≈ √2` (Y10 corrected)
8. `test_compute_ecm_to_fa_bias_active_rank1_aligned_score_half`
   — `T = n⊗n` (rank-1 aligned); assert `score = 0.5` (rank-1 minus
   isotropic-half)

### Multiplier behavior (3)

9. `test_compute_ecm_to_fa_bias_active_aligned_target_boosts_above_neutral`
10. `test_compute_ecm_to_fa_bias_active_orthogonal_target_suppresses_below_neutral`
11. `test_compute_ecm_to_fa_bias_active_all_three_rates_same_multiplier_per_fa`
    (Y4 scaffolding contract — assert `multipliers_per_fa[i, 0] ==
    multipliers_per_fa[i, 1] == multipliers_per_fa[i, 2]` for every i)

### Zero-traction + boundary (3, Y7)

12. `test_compute_ecm_to_fa_bias_active_zero_traction_returns_exact_neutral`
    (multiplier == 1.0 exact; alignment_score == 0.0 exact)
13. `test_compute_ecm_to_fa_bias_active_zero_traction_count_in_diagnostics`
14. `test_compute_ecm_to_fa_bias_active_empty_adhesions_returns_empty_consistent_shape`
    (`multipliers_per_fa.shape == (0, 3)`; `sampled_diagnostics` non-None)

### Diagnostics (3, Y6+Y13)

15. `test_compute_ecm_to_fa_bias_active_diagnostics_dict_has_exact_10_keys`
    (assert `set(result.diagnostics_dict.keys()) == {...10 keys...}`)
16. `test_compute_ecm_to_fa_bias_active_diagnostics_dict_no_loose_extra_keys`
    (assert no surprise keys via subset check + length)
17. `test_compute_ecm_to_fa_bias_active_diagnostics_sampler_geometry_matches_v1`
    (Y13: `result.diagnostics_dict["sampler_geometry"] ==
    "bilinear_cell_centered"` exactly — sister with HB#4 v1)

### Failure-kind (1) [tests 1–4 already cover `k_active_invalid`; this
tests `active_multiplier_non_finite` separately]

18. `test_compute_ecm_to_fa_bias_active_multiplier_non_finite_failure_kind`
    (monkeypatch `sampled.orientation_tensor` to contain NaN/inf;
    assert `FAToECMBiasError` with `failure_kind ==
    "active_multiplier_non_finite"`)

### Validation order (1, Y12)

19. `test_compute_ecm_to_fa_bias_active_invalid_k_active_raises_before_sampling`
    (pass empty adhesions + bad `k_active`; assert raises
    `k_active_invalid` before `sample_ecm_at_fa_positions` is called —
    monkeypatch sampler to track invocation)

### Effective-stiffness guard (1, Y8)

20. `test_compute_ecm_to_fa_bias_active_does_not_reference_effective_stiffness`
    — function-scoped AST + string guard via
    `textwrap.dedent(inspect.getsource(compute_ecm_to_fa_bias_active))`;
    assert no `Attribute.attr` / `Name.id` reference to forbidden set
    `{effective_stiffness, stiffness_kpa, fiber_density, ligand_density,
    accumulated_traction_nNs_per_um2}` AND no `"effective_stiffness"`
    in function source string (HB#5 Y11 sister)

### Exports + private constant (1, Y14)

21. `test_compute_ecm_to_fa_bias_active_only_function_exported`
    — `compute_ecm_to_fa_bias_active` importable from BOTH
    `acs.v2.dynamics` AND `acs.v2`; `_DEVIATORIC_SCORE_BOUND` NOT
    importable from either; `_validate_k_active` NOT importable

---

## 5. Files

- `docs/v2_hard_blocker_4_active_brief.md` (existing, opening brief,
  commit `013e161`)
- `docs/v2_hard_blocker_4_active_locked.md` (this file, source of truth)
- `acs/v2/dynamics/ecm_to_fa_bias.py` (MODIFY):
  - Add `_DEVIATORIC_SCORE_BOUND: Final[float] = float(np.sqrt(2.0))`
    private constant
  - Add `_validate_k_active(k_active: float) -> None` private helper
  - Add `compute_ecm_to_fa_bias_active(adhesions, ecm, *, k_active)
    -> ECMToFABiasResult` public function
  - Extend `FAToECMBiasFailureKind` Literal: add `"k_active_invalid"`
    and `"active_multiplier_non_finite"`
- `acs/v2/dynamics/__init__.py` (1 new export +
  alphabetical reorder of `__all__`)
- `acs/v2/__init__.py` (mirror export)
- `tests/test_v2_ecm_to_fa_bias_active.py` (NEW, 21 tests per §4) OR
  add to existing `tests/test_v2_ecm_to_fa_bias.py` (impl decision per
  Sanity Gate)

---

## 6. References

- Brief (opening position):
  `docs/v2_hard_blocker_4_active_brief.md` (commit `013e161`)
- Parent locked plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- Sister locks:
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` (HB#4 v1
    neutral; same return type, shared sampler)
  - `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md` (HB#1+#2;
    Y4 branch-defined zero-traction sister)
  - `docs/v2_phase_e_composition_locked.md` (Phase E v1; uses HB#4 v1
    neutral; v2 composition is separate)
  - `docs/v2_hard_blocker_5_lyapunov_metric_locked.md` (Y11 AST-walk
    function-scoped guard sister-pattern)
  - `docs/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1 typed
    schema discipline; here grandfathered for `ECMToFABiasResult`)
  - `docs/v2_item_5_sweep_harness_locked.md` (Item 5 sweep harness;
    uses `0.5·I` IC which is automatically neutral under this active
    law per Y1 deviatoric form)
- Verified at source for Y5 (ECMToFABiasResult schema):
  `acs/v2/dynamics/ecm_to_fa_bias.py:128` (5 fields: `fa_ids`,
  `multipliers_per_fa`, `rate_names`, `sampled_diagnostics`,
  `diagnostics_dict`)
- Verified at source for Y9 (FAToECMBiasFailureKind):
  `acs/v2/dynamics/ecm_to_fa_bias.py:89` (extending 2-kind literal to
  4-kind)
- Sister implementation precedent (HB#4 v1 grandfathered shared
  sampler): `acs/v2/dynamics/ecm_to_fa_bias.py:244`
  (`sample_ecm_at_fa_positions`)
- Adversarial debate posture: memory
  `feedback_aggressive_design_debate.md`, PI id=809
- Hard Rule 10 inline derivation: memory
  `rule10_unit_derivation_in_docs.md`, Codex id=1234 (deviatoric bound
  derivation §3 mandatory)
- Hard Rule 11 wording-boundary meta-test: memory
  `hard_rule_11_wording_boundary_meta_test.md` (Y1 deviatoric form is
  the fourth Phase E composition Hard Rule 11 catch family — Phase E
  v1 Y1, HB#5 Y1, Item 5 Y1+Y4+Y17, HB#4-active Y1)
- Step 6 pre-commit batch: memory `design_note_pre_commit_batch.md`,
  Codex id=1428

---

## 7. Cross-room dispatch

This file is the design-team input to implementation-work for:

1. impl Claude writes HB#4-active step 1 Sanity Gate doc
   (`docs/v2_hard_blocker_4_active_sanity_gate.md`) from this lock —
   6 sections per §3, with deviatoric bound derivation as central
   anchor + isotropic-IC anti-pattern guard prominent.
2. impl Codex review (6 focus per Codex `id=1612` final SEAL):
   - **Deviatoric form correctness**: `Q = T - 0.5·trace(T)·I`; isotropic
     `T = a·I` → `score = 0` → multiplier = 1.0 exactly (Y1)
   - **Schema-tight √2 example correctness**: `T = [[1, 1], [1, -1]]`
     traceless, eigenvalues ±√2, eigenvector at θ = π/8 (Y10 corrected
     math)
   - **Validation order**: `ecm.validate() → _validate_k_active → sampler`
     (Y12); invalid k_active raises before sampler called
   - **Schema-corrected fields**: `fa_ids`, `multipliers_per_fa`,
     `rate_names`, `sampled_diagnostics`, `diagnostics_dict` (Y5)
   - **Failure-kind extension**: `k_active_invalid`,
     `active_multiplier_non_finite` (Y9); finite defensive guards on
     score and multiplier (Y11)
   - **Function-scoped AST + string guard**: no
     `effective_stiffness` / raw mechanosensing field references in
     `compute_ecm_to_fa_bias_active` body (Y8)
   - **Single export**: only `compute_ecm_to_fa_bias_active` exported
     through both `__init__.py`s (Y14); `_DEVIATORIC_SCORE_BOUND`
     private
3. On Sanity Gate PASS: HB#4-active step 1 code commit
   (`acs/v2/dynamics/ecm_to_fa_bias.py` MODIFY: add private constant +
   private helper + public function + 2 new failure kinds; 1 new
   export through both `__init__.py`s; new test file or extend
   existing).
4. After commit: design-discussion idle on HB#4-active step 1; **all
   upstream Phase E v2 step 1 design sealed**. Next design entries
   (separate cycles):
   - Phase E v2 composition step 2 (`step_closed_loop_phase_e_v2`
     wrapper that swaps v1 neutral for active step 1) — separate
     adversarial-debate cycle, requires HB#4-active step 1 impl PASS
     first
   - Per-rate selectivity (Phase E v3) — separate cycle requiring
     per-rate literature anchors
   - effective_stiffness law (separate unit, only if v3 consumes it)

Rounds 1–3 of the design lock are MCP id 1605–1612.
