# V2 Closed-Loop ECM Gate — Hard Blocker #1+#2 Constitutive Direction (Orientation Response) — Locked Design

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (4-round adversarial lock,
PI id=809 aggressive debate posture applied, PI id=1008/1057/1161
autonomy + visible deliverable focus)
**Source unit**: design-discussion `topic=v2-hard-blocker-1-constitutive-direction`,
MCP id 1466–1476 (rounds 1–4 + seal ack)
**Brief**: `docs/v2_hard_blocker_1_constitutive_direction_brief.md` (commit `4950302`)
**Parent locked plan**: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
(Phase E entry, Hard Blockers #1 and #2 merged)
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for the HB#1+#2 Sanity Gate doc + code entry.

---

## 0. Scope

This unit locks the **first active constitutive law** for closed-loop
ECM evolution under FA traction stimulus. Hard Blocker #1 (constitutive
direction) and Hard Blocker #2 (saturation form) are **merged** because
the chosen exact-exponential convex update has saturation built-in
(convex weight `w ∈ [0,1]` is itself the bound).

**What is updated**: `orientation_tensor` only. `stiffness_kpa`,
`fiber_density`, `ligand_density`, and `accumulated_traction_nNs_per_um2`
are NOT modified by this law. This is a deliberate **orientation-only
first active law** decision (Y1, see §2).

**What drives it**: instantaneous `traction_density_xy` from HB#3
(`scatter_fa_traction_to_ecm_bilinear`). `accumulated_traction_nNs_per_um2`
is NOT consumed (Y6).

**Why merged HB#1+#2**: the convex combination
`T_new = (1−w)·T_old + w·(n⊗n)` with `w = 1 − exp(−k·S·dt)` is itself
the saturation; no separate hard cap is needed (Q4 / Y10).

**Why now**: Phase E activation requires the constitutive direction
locked before any ECM-mutation can be permitted. Per the locked phased
plan §1, Phase D no-op scaffolding is committed; Phase E is BLOCKED on
HB#1, HB#2, HB#5 + the effective_stiffness law decision. This unit
unblocks HB#1+#2.

---

## 1. Final Lock Summary

### Constants (literature-pinned)

```python
TRACTION_REF_NN_PER_UM2: Final[float] = 1.0   # = 1.0 kPa
TAU_ALIGN_RANGE_S:       Final[tuple[float, float]] = (3600.0, 7200.0)
K_ORIENT_PER_S:          Final[float] = float(
    1.0 / np.sqrt(TAU_ALIGN_RANGE_S[0] * TAU_ALIGN_RANGE_S[1])
)  # = 1.964e-4 1/s ≈ 2.0e-4, geometric midpoint of declared range
```

**Hard Rule 10 unit chain** (inline, mandatory per
`rule10_unit_derivation_in_docs` memory):

- `1.0 nN/μm² = (1×10⁻⁹ N) / (1×10⁻¹² m²) = 1×10³ Pa = 1.0 kPa` ✓
  (matches collagen substrate FA traction order from Munevar/Wang/Dembo
  2001 + Tan/Tien/Pirone/Chen 2003)
- `traction_norm [nN/μm²] / TRACTION_REF_NN_PER_UM2 [nN/μm²] = S [-]`
  dimensionless ✓
- `K_ORIENT_PER_S [1/s] · S [-] · dt [s] = [-]` dimensionless argument
  to `exp(·)` ✓

**Magic-Number Block compliance**:
- `TRACTION_REF_NN_PER_UM2 = 1.0`: derivable (literature ✓), grid-invariant
  (no dx/dt dependence ✓), not fitted (chosen before any A/A₀
  comparison ✓)
- `K_ORIENT_PER_S = 1.964e-4`: derivable (geometric midpoint of declared
  literature range `[3600, 7200]` s ✓), grid-invariant ✓, not fitted ✓

### Typed result + diagnostics

```python
@dataclass(frozen=True, slots=True)
class ECMOrientationResponseDiagnostics:
    """Typed diagnostics for one HB#1+#2 ECM orientation response step."""
    n_nonzero_cells: int
    n_total_cells: int
    traction_ref_nN_per_um2: float                          # caller-passed value
    k_orient_per_s: float                                   # caller-passed value
    max_traction_norm_nN_per_um2: float                     # Rule 10 bridge
    max_dimensionless_stimulus: float                       # Rule 10 bridge
    max_convex_weight: float                                # saturation indicator
    mean_convex_weight_nonzero: float
    max_distance_to_target_frobenius_nonzero: float
    mean_distance_to_target_frobenius_nonzero: float
    max_orientation_delta_frobenius: float
    symmetry_residual_max: float                            # schema sanity
    componentwise_bound_residual_max: float                 # schema sanity


@dataclass(frozen=True, slots=True)
class ECMOrientationResponseResult:
    """Outputs of one HB#1+#2 ECM orientation response step."""
    updated_ecm: ECMSubstrateState
    diagnostics: ECMOrientationResponseDiagnostics
```

### Function

```python
def step_ecm_orientation_response(
    ecm: ECMSubstrateState,
    traction_density_xy: np.ndarray,                # (nx, ny, 2) nN/μm², HB#3 output
    dt_s: float,
    *,
    traction_ref_nN_per_um2: float = TRACTION_REF_NN_PER_UM2,
    k_orient_per_s: float = K_ORIENT_PER_S,
) -> ECMOrientationResponseResult:
    # --- Validate (ecm.validate() FIRST per HB#3/HB#4 sister precedent + Codex
    #     id=1486 silent-heal-path closure: invalid input ECM is rejected,
    #     not healed) ---
    ecm.validate()                                                        # symmetry + componentwise + finite, on input
    _validate_finite_traction(traction_density_xy)                       # nan/inf raise
    _validate_dt_nonneg_finite_nonbool(dt_s)                             # dt>=0, finite, not bool
    _validate_positive_finite(traction_ref_nN_per_um2, "traction_ref")   # strict > 0
    _validate_positive_finite(k_orient_per_s, "k_orient")                # strict > 0, no zero no-op
    _validate_traction_shape(traction_density_xy, ecm.grid_shape)        # (*ecm.grid_shape, 2)

    # --- Stimulus + direction (branch-defined, NO epsilon) ---
    traction_norm = np.linalg.norm(traction_density_xy, axis=-1)
    nonzero = traction_norm > 0.0                                         # exact componentwise
    S = traction_norm / traction_ref_nN_per_um2                           # dimensionless

    n_xy = np.zeros_like(traction_density_xy, dtype=np.float64)
    n_xy[nonzero] = traction_density_xy[nonzero] / traction_norm[nonzero, None]
    T_target = np.einsum("...i,...j->...ij", n_xy, n_xy)

    # --- Convex weight (exact exponential, expm1 numerically stable near 0) ---
    w = -np.expm1(-k_orient_per_s * S * dt_s)                             # in [0,1]

    # --- Orientation update (componentwise convex combination) ---
    orientation_old = ecm.orientation_tensor
    orientation_new = (
        (1.0 - w[..., None, None]) * orientation_old
        + w[..., None, None] * T_target
    )

    # --- Build new ECM (NO aliasing — copy all arrays per open-loop precedent) ---
    updated_ecm = ECMSubstrateState(
        origin_um_xy=tuple(ecm.origin_um_xy),
        spacing_um=float(ecm.spacing_um),
        stiffness_kpa=np.array(ecm.stiffness_kpa, dtype=np.float64, copy=True),
        ligand_density=np.array(ecm.ligand_density, dtype=np.float64, copy=True),
        fiber_density=np.array(ecm.fiber_density, dtype=np.float64, copy=True),
        orientation_tensor=np.array(orientation_new, dtype=np.float64, copy=True),
        accumulated_traction_nNs_per_um2=np.array(
            ecm.accumulated_traction_nNs_per_um2, dtype=np.float64, copy=True
        ),
        source=ecm.source,
    )
    updated_ecm.validate()

    # --- Diagnostics (typed, sister-pattern of B1) ---
    diff = orientation_new - T_target
    delta = orientation_new - orientation_old
    if nonzero.any():
        mean_w_nz = float(w[nonzero].mean())
        diff_nz_norms = np.linalg.norm(diff[nonzero], axis=(-2, -1))
        max_dist = float(diff_nz_norms.max())
        mean_dist = float(diff_nz_norms.mean())
    else:
        mean_w_nz = 0.0
        max_dist = 0.0
        mean_dist = 0.0
    diagnostics = ECMOrientationResponseDiagnostics(
        n_nonzero_cells=int(nonzero.sum()),
        n_total_cells=int(np.prod(ecm.grid_shape)),
        traction_ref_nN_per_um2=float(traction_ref_nN_per_um2),
        k_orient_per_s=float(k_orient_per_s),
        max_traction_norm_nN_per_um2=float(traction_norm.max()),
        max_dimensionless_stimulus=float(S.max()),
        max_convex_weight=float(w.max()),
        mean_convex_weight_nonzero=mean_w_nz,
        max_distance_to_target_frobenius_nonzero=max_dist,
        mean_distance_to_target_frobenius_nonzero=mean_dist,
        max_orientation_delta_frobenius=float(np.linalg.norm(delta, axis=(-2, -1)).max()),
        symmetry_residual_max=float(
            np.abs(orientation_new - np.swapaxes(orientation_new, -1, -2)).max()
        ),
        componentwise_bound_residual_max=float(max(0.0, np.abs(orientation_new).max() - 1.0)),
    )

    return ECMOrientationResponseResult(
        updated_ecm=updated_ecm,
        diagnostics=diagnostics,
    )
```

### Forbidden in HB#1+#2

- Update of `stiffness_kpa`, `fiber_density`, `ligand_density`, or
  `accumulated_traction_nNs_per_um2` (orientation-only first law)
- Consumption of `accumulated_traction_nNs_per_um2` as driver
  (instantaneous traction only; future memory-driven law is a separate
  unit and must redo the unit chain because units are nN·s/μm²)
- ε / epsilon as physical or numerical regularizer in direction
  computation (branch-defined zero-traction only)
- Hard cap or clamp on `orientation_tensor` components (saturation is
  built into the convex weight)
- PSD / eigenvalue / trace claims beyond the schema's componentwise
  `|T_ij| ≤ 1` invariant (schema strengthening is a separate unit)
- Aliasing of unchanged arrays into the returned ECMSubstrateState
  (all 5 array fields must be fresh copies)
- Generic[...] parameterization or `# type: ignore` preemptive on the
  result/diagnostics dataclasses (B1 sister-precedent)
- Silent zero-K no-op (`k_orient_per_s == 0` raises; disabled state must
  use a separate entry-point or explicit `enabled=False` flag in a
  future wrapper)
- Mechanosensing shortcut: any future downstream consumer that uses
  `orientation_tensor` without explicit decision on whether/how to
  weight by `fiber_density` is **disallowed by the latent-orientation
  contract** (Y7, see §2)

---

## 2. Reasoned-acceptance trace (Y1–Y16)

The lock converged after 4 rounds. Each Y is a Claude concession with
the round in which it was accepted, preserving the adversarial-debate
audit trail (PI id=809 posture).

**Y1 (round 2, accepted Codex C1 = orientation-only first law)**:
- Claude opening lean: both `fiber_density` and `orientation_tensor`
  update in HB#1.
- Codex challenge: Hall 2016 / Notbohm-type evidence supports
  force-driven fiber alignment, NOT local density compaction or
  deposition. Pointwise `fiber_density += ...` creates collagen with no
  mass/deposition bookkeeping — biologically incoherent without an
  explicit semantic choice (compaction-flux conservative / deposition
  source / image-density proxy).
- Resolution: HB#1 first active law updates `orientation_tensor` only;
  `fiber_density` is unchanged. Future `fiber_density` law is a
  separate unit and must declare its semantic choice.

**Y2 (round 2, accepted Codex C2 = monotone-toward-target wording)**:
- Brief said "monotone-positive" for orientation, which is incorrect
  componentwise.
- Codex catch: tensor components can decrease (e.g., isotropic →
  rank-1 update). The defensible invariant is "distance-to-target
  decreasing" (Frobenius or any norm).
- Resolution: convex update with `w ∈ [0,1]` guarantees
  `‖T_new − T_target‖ ≤ ‖T_old − T_target‖` per step (fixed target).

**Y3 (round 2, accepted Codex C3 = dimensionless stimulus + k_orient
[1/s])**:
- Brief proposed `α [μm²/(nN·s)]` derived "from a timescale", which is
  a Hard Rule 10 violation (a timescale gives `k [1/s]`, not `α`).
- Codex fix: define `S = |traction| / τ_ref` dimensionless with `τ_ref
  [nN/μm²]` from literature; then `dT/dt = k_orient · S · (T_target −
  T)` with `k_orient [1/s]` from alignment timescale.
- Resolution: adopted dimensionless-stimulus formulation. α never
  appears.

**Y4 (round 2, accepted Codex C4 = branch-defined zero-traction, no
epsilon)**:
- Brief used `n_stim = traction / max(|traction|, ε)`.
- Codex catch: ε is a magic number (Magic-Number Block test 1 fails:
  no derivation). Constitutive law must not include numerical
  regularizers.
- Resolution: branch logic — `if |traction| == 0: rate = 0 exactly;
  else: n = traction / |traction|`. ε only used in dtype finite
  validation, not in physics.

**Y5 (round 2, accepted Codex C5 = schema invariant only)**:
- Brief claimed `eigenvalues in [0,1]` for orientation target, which
  is stronger than the actual schema invariant.
- Codex catch: current `ECMSubstrateState.validate()` enforces
  symmetry + componentwise `|T_ij| ≤ 1`, NOT PSD / eigenvalue / trace.
  Convex update is component-safe under the actual invariant.
- Resolution: lock uses schema's componentwise bound; PSD strengthening
  is a separate unit.

**Y6 (round 2, accepted Codex C6 = instantaneous traction; accumulated
stays diagnostic)**:
- Brief was ambiguous on whether `accumulated_traction_nNs_per_um2` or
  instantaneous `traction_density_xy` drives the law.
- Codex requested explicit pin.
- Resolution: HB#1 first law consumes instantaneous traction only;
  accumulated remains HB#4 sampled diagnostic, not consumed by HB#1.
  Future memory-driven law is a separate unit.

**Y7 (round 3, accepted Codex A1 = (5b) ungated + latent-orientation
caveat)**:
- Open question Q5: should orientation update gate on `fiber_density >
  0`?
- Codex chose (5b) ungated with reinforcing caveat: "orientation_tensor
  is a latent/local collagen-network orientation state even where
  current density is zero or unobservable. Any future downstream
  physical effect that consumes orientation must decide whether/how to
  weight by fiber_density; HB#1 does not grant an orientation-only
  mechanosensing shortcut."
- Resolution: ungated update + the latent-orientation caveat is
  pinned as a forbidden-list item in §1 (no mechanosensing shortcut).

**Y8 (round 3, accepted Codex A2 = rename TAU_REF → TRACTION_REF)**:
- Claude naming `TAU_REF_NN_PER_UM2` was confusing (`tau` reads like
  time constant, conflicts with `τ_align`).
- Resolution: `TRACTION_REF_NN_PER_UM2`. Dimension chain explicit:
  `traction_norm [nN/μm²] / TRACTION_REF_NN_PER_UM2 [nN/μm²] = S [-]`.

**Y9 (round 3, accepted Codex A3 = K_ORIENT geometric midpoint
derivation)**:
- Claude proposed `K_ORIENT_PER_S = 2.0e-4` as exact value, which would
  appear arbitrary.
- Codex fix: derive default as geometric midpoint of declared range
  `tau_align_range_s = [3600, 7200]` → `K_default = 1 /
  sqrt(3600·7200) = 1.964e-4 ≈ 2.0e-4`. Rationale: rate constants
  spanning ~2× factors have principled log-scale midpoint; geometric
  midpoint is proportion-symmetric, arithmetic is high-end-biased.
- Resolution: declared range is the literature pin; default is
  geometric midpoint, derivable. Magic-Number Block test 1 passes.

**Y10 (round 3, accepted Codex A4 = expm1 numerical stability)**:
- Claude pseudocode used `1.0 - np.exp(-x)`.
- Codex fix: `-np.expm1(-x)` is mathematically identical but
  numerically stable near `x=0`. ECM evolution often runs at small
  `K·S·dt` (e.g., K=2e-4, S=1, dt=60 → x=0.012), where catastrophic
  cancellation in `exp(-x) - 1` matters.
- Resolution: lock uses `-np.expm1(-K·S·dt)`.

**Y11 (round 3, accepted Codex A5 = 5 validation/failure contracts)**:
- Brief had loose validation specification.
- Codex specified: nonfinite traction reject; negative/bool/nonfinite
  `dt_s` reject; nonpositive/nonfinite `traction_ref` reject;
  `k_orient` strict `> 0` (no zero no-op silent path); shape mismatch
  reject; no-mutation contract; component bound proof cites schema
  invariant only; distance-to-target wording is per-step (fixed target),
  NOT global.
- Resolution: all 5 contracts pinned in §1 + §3 + §4.

**Y12 (round 4, accepted Codex B1 = schema correction + no-aliasing)**:
- Claude round 3 pseudocode had wrong field names (`dx_um` vs actual
  `spacing_um`) and aliased unchanged arrays into the new ECM.
- Codex catch: verified against `acs/v2/ecm_substrate.py:42` —
  fields are `origin_um_xy`, `spacing_um`, `stiffness_kpa`,
  `ligand_density`, `fiber_density`, `orientation_tensor`,
  `accumulated_traction_nNs_per_um2`, `source`; `grid_shape` is a
  property, not a constructor arg. Numpy arrays in frozen dataclasses
  remain mutable — sharing arrays risks contamination via downstream
  mutation. Open-loop precedent (`ecm_open_loop.py`) copies all arrays
  on return.
- Resolution: lock pseudocode uses correct schema, copies all 5 array
  fields with `np.array(..., copy=True)`, calls `updated_ecm.validate()`
  before returning. Tests assert bytewise equality + `not
  np.shares_memory(...)` for all 5 array fields.

**Y13 (round 4, accepted Codex B2 = 13-field typed diagnostics)**:
- Claude proposed 8 fields; Codex revised to 13.
- Added: `traction_ref_nN_per_um2` + `k_orient_per_s` (sweep
  reproducibility traceback), `max_traction_norm_nN_per_um2` +
  `max_dimensionless_stimulus` (Rule 10 bridge to weight),
  `_nonzero` suffix on means/distances (avoids zero-traction placeholder
  confusion).
- Resolution: 13 fields. Field count justified by HB#1+#2 being the
  most complex constitutive law (vs B1 4-field /
  `FocalAdhesionDynamicsDiagnostics`, Phase D 6-field).

**Y14 (round 4, accepted Codex B3 = ECMOrientationResponseResult)**:
- Codex named the result dataclass shape; Claude accepted.
- Phase D `FAToECMResponseResult` sister-pattern.

**Y15 (round 4, accepted Codex B4 = dt_s == 0 valid no-op)**:
- Claude proposed rejecting zero `dt_s`.
- Codex correction: `dt_s == 0` should be valid → `w = -expm1(0) = 0`
  → identity update + diagnostics weights all zero. No special branch
  needed. `dt_s < 0` rejected (nonsensical), boolean rejected
  (open-loop precedent), nonfinite rejected.
- Resolution: `_validate_dt_nonneg_finite_nonbool`.

**Y16 (round 4, accepted Codex B5 = lock structure + module path)**:
- Codex confirmed lock filename
  `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md` and
  module path `acs/v2/dynamics/ecm_constitutive_response.py`.
- Resolution: §0–§7 layout adopted (HB#3/#4/Phase D/B1 sister-pattern).

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

1. **Units**:
   - `traction_density_xy`: nN/μm² (HB#3 output)
   - `traction_ref_nN_per_um2`: nN/μm² (= kPa via inline derivation)
   - `S = traction_norm / traction_ref`: dimensionless
   - `k_orient_per_s`: 1/s (geometric midpoint of literature alignment
     timescale range)
   - `dt_s`: s
   - `k · S · dt`: dimensionless (argument to exp)
   - `w`: dimensionless ∈ [0, 1]
   - `T_target = n⊗n`: dimensionless, components ∈ [-1, 1]
   - `orientation_new`: dimensionless, components ∈ [-1, 1] by
     componentwise bound proof
2. **Boundary**:
   - `n_adhesions == 0` (caller did not run HB#3): traction array all
     zero, all `w == 0`, identity update, diagnostics zero
   - `dt_s == 0`: identity update (Y15)
   - `dt_s → ∞`: `w → 1` everywhere `S > 0`, `T → T_target`
     immediately (asymptotic-bound is the saturation, no NaN)
   - `K_ORIENT → ∞`: same as `dt → ∞`
   - `traction_ref → ∞`: `S → 0`, `w → 0`, identity update (no NaN)
   - All-zero traction: branch handles cleanly (no division)
   - Boolean `dt_s`: rejected
3. **Conservation**: orientation update is a convex combination —
   conserves nothing physical (orientation has no conservation law),
   but **preserves** the schema invariant under the contract "valid
   input ECM → valid output ECM" (componentwise bound, symmetry).
   Invalid input ECM is **rejected** by `ecm.validate()` at
   function entry per Codex `id=1486`, not "healed" toward
   `T_target` by the convex update. `stiffness_kpa`,
   `fiber_density`, `ligand_density`,
   `accumulated_traction_nNs_per_um2` all bytewise unchanged.
4. **Numerical**:
   - `expm1` chosen for stability near zero argument (Y10)
   - `np.float64` enforced on all arrays (no f32/f64 mixing)
   - Convex combination is unconditionally stable in `dt` (no CFL)
5. **Sign**: convex weight `w ≥ 0` always (since `K, S, dt ≥ 0` and
   `expm1` of nonpositive arg is in [-1, 0]). Distance to target
   monotonically decreasing per step under fixed target (Y2).
6. **Measurement-protocol**: orientation tensor is a per-cell
   dimensionless tensor on the ECM grid; the law is per-cell pointwise
   (no spatial coupling, no averaging window). The off-proof concern
   from the v13 Stage 1a episode does not apply because the law is
   strictly local and the "proof" is the convex-combination algebra
   that holds at every cell independently.

---

## 4. Test Catalog (~20 tests; count not capped)

1. `test_orientation_response_zero_traction_is_identity`
2. `test_orientation_response_dt_zero_is_identity_no_op`
3. `test_orientation_response_dt_negative_raises` —
   `ECMConstitutiveResponseError`, `failure_kind="dt_invalid"`
4. `test_orientation_response_dt_bool_raises` —
   `ECMConstitutiveResponseError`, `failure_kind="dt_invalid"`
5. `test_orientation_response_dt_nonfinite_raises` —
   `ECMConstitutiveResponseError`, `failure_kind="dt_invalid"`
6. `test_orientation_response_traction_ref_nonpositive_raises`
   — `ECMConstitutiveResponseError`,
   `failure_kind="traction_ref_invalid"`
7. `test_orientation_response_k_orient_nonpositive_raises`
   (strict, no zero no-op) —
   `ECMConstitutiveResponseError`,
   `failure_kind="k_orient_invalid"`
8. `test_orientation_response_traction_shape_mismatch_raises` —
   `ECMConstitutiveResponseError`,
   `failure_kind="traction_shape_mismatch"`
9. `test_orientation_response_traction_nonfinite_raises` —
   `ECMConstitutiveResponseError`,
   `failure_kind="non_finite_traction"`
10. `test_orientation_response_no_aliasing_all_arrays` —
    `not np.shares_memory(...)` for all 5 array fields
11. `test_orientation_response_unchanged_fields_bytewise_equal` —
    `stiffness_kpa`, `ligand_density`, `fiber_density`,
    `accumulated_traction_nNs_per_um2` bytewise equal
12. `test_orientation_response_componentwise_bound_preserved` —
    `|T_new_ij| ≤ 1` everywhere, `updated_ecm.validate()` PASS
13. `test_orientation_response_symmetry_preserved` —
    `T_new == T_new.T` everywhere
14. `test_orientation_response_distance_to_target_monotone_per_step` —
    fixed traction direction, `‖T_new − T_target‖_F ≤ ‖T_old −
    T_target‖_F`
15. `test_orientation_response_uniform_traction_converges_to_target` —
    apply many steps with constant traction → `T → n⊗n`
16. `test_orientation_response_diagnostics_typed_dataclass` —
    `isinstance(result.diagnostics, ECMOrientationResponseDiagnostics)`,
    `with pytest.raises(TypeError): result.diagnostics["max_convex_weight"]`
    (B1 sister-pattern meta-test)
17. `test_orientation_response_diagnostics_records_input_params` —
    `result.diagnostics.traction_ref_nN_per_um2 == passed value`,
    `result.diagnostics.k_orient_per_s == passed value`
18. **`test_orientation_response_exports_through_both_init`** —
    7 symbols importable from both `acs.v2.dynamics` and
    `acs.v2`: `ECMConstitutiveResponseError`,
    `ECMOrientationResponseDiagnostics`,
    `ECMOrientationResponseResult`, `step_ecm_orientation_response`,
    `TRACTION_REF_NN_PER_UM2`, `K_ORIENT_PER_S`,
    `TAU_ALIGN_RANGE_S` (per Codex `id=1488` 6→7 expansion to
    include the locked error class)
19. **`test_orientation_response_uses_current_schema_not_stale_naming`** —
    Y12-guard meta-test: implementation must not reference `dx_um` or
    pass `grid_shape` to the constructor; must call **both**
    `ecm.validate()` at function entry and `updated_ecm.validate()`
    before returning. (Static check via `inspect.getsource()` regex
    or AST walk.) Updated per Codex `id=1486` to require both
    validations, not only the output one.
20. **`test_orientation_response_invalid_input_ecm_rejected_before_update`**
    (Codex `id=1486`) — passing an `ECMSubstrateState` whose
    `orientation_tensor` violates `validate()` (e.g., a component
    `> 1.0` that the convex update could "heal" toward `T_target`
    under high traction × dt) must raise plain schema
    `ValueError` from `ecm.validate()` at function entry, BEFORE
    any algebra runs. Closes the silent-heal path: invalid
    input ECM is rejected, not healed. Inherited `ValueError`
    is NOT wrapped in `ECMConstitutiveResponseError` (per Codex
    `id=1488`: schema-level errors stay schema-level, not masked
    as constitutive-law-level errors).

---

## 5. Files

- `docs/v2_hard_blocker_1_constitutive_direction_brief.md` (existing,
  opening brief, commit `4950302`)
- `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md` (this
  file, source of truth)
- `acs/v2/dynamics/ecm_constitutive_response.py` (new):
  - `TRACTION_REF_NN_PER_UM2`, `TAU_ALIGN_RANGE_S`, `K_ORIENT_PER_S`
    module constants
  - `ECMOrientationResponseDiagnostics` dataclass (frozen, slots, 13
    fields)
  - `ECMOrientationResponseResult` dataclass (frozen, slots)
  - `step_ecm_orientation_response` function
  - **`ECMConstitutiveResponseError(ValueError)`** with
    machine-readable `failure_kind` attribute (locked per Codex
    `id=1488` Step 6 sister-gate-mirror with HB#3
    `FAToECMScatteringError` + HB#4 `FAToECMBiasError`). Owned
    failure kinds (raised by validators inside
    `step_ecm_orientation_response`): `dt_invalid` (covers
    negative / bool / nonfinite), `traction_ref_invalid`
    (nonpositive / nonfinite), `k_orient_invalid` (nonpositive /
    nonfinite — strict no-zero-no-op), `traction_shape_mismatch`,
    `non_finite_traction`. Inherited `ecm.validate()` failures
    remain plain schema `ValueError` (per Codex `id=1488` to
    avoid masking schema-level errors as constitutive-law-level
    errors).
- `acs/v2/dynamics/__init__.py` (**7 new exports** per Codex
  `id=1488` API-surface lock + alphabetical reorder of
  `__all__`): error + 3 constants + 2 dataclasses + 1 function
- `acs/v2/__init__.py` (mirror 7 exports)
- `tests/test_v2_ecm_constitutive_response.py` (new, ~19 tests per §4)

---

## 6. References

- Brief (opening position):
  `docs/v2_hard_blocker_1_constitutive_direction_brief.md` (commit
  `4950302`)
- Parent locked plan:
  `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- Sister locks:
  - `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2_phase_d_no_op_scaffolding_locked.md`
  - `docs/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1)
- ECM substrate schema (verified for Y12):
  `acs/v2/ecm_substrate.py:42`
- ECM open-loop precedent for array-copy on return:
  `acs/v2/dynamics/ecm_open_loop.py`
- Literature (TRACTION_REF):
  - Munevar S, Wang Y, Dembo M (2001). "Traction force microscopy of
    migrating normal and H-ras transformed 3T3 fibroblasts." Biophys J
    80(4):1744–1757.
  - Tan JL, Tien J, Pirone DM, Gray DS, Bhadriraju K, Chen CS (2003).
    "Cells lying on a bed of microneedles: an approach to isolate
    mechanical force." PNAS 100(4):1484–1489.
- Literature (TAU_ALIGN):
  - Hall MS, Alisafaei F, Ban E, Feng X, Hui CY, Shenoy VB, Wu M
    (2016). "Fibrous nonlinear elasticity enables positive mechanical
    feedback between cells and ECMs." PNAS 113(49):14043–14048.
  - Notbohm J, Lesman A, Rosakis P, Tirrell DA, Ravichandran G (2015).
    "Microbuckling of fibrin provides a mechanism for cell mechanosensing."
    J R Soc Interface 12(108):20150320.
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

1. impl Claude writes HB#1+#2 Sanity Gate doc
   (`docs/v2_hard_blocker_1_2_constitutive_response_sanity_gate.md`)
   from this lock — 6 sections per §3.
2. impl Codex review (5 focus per Codex `id=1476` + Y12-guard):
   - exact-exponential convex update form + `expm1` numerical stability
   - no aliasing of any of the 5 ECM array fields; `validate()` called
   - typed diagnostics records all input parameters + Rule 10 bridge
     fields
   - `dt_s == 0` valid no-op, `dt_s < 0` / bool / nonfinite rejected;
     `k_orient` and `traction_ref` strict positive
   - Y12-guard: `spacing_um` (not `dx_um`), no `grid_shape` arg, no
     mechanosensing shortcut
   - all 6 new symbols exported through both `acs.v2.dynamics` and
     `acs.v2`
3. On Sanity Gate PASS: HB#1+#2 code commit
   (`acs/v2/dynamics/ecm_constitutive_response.py` new, exports,
   tests).
4. After commit: design-discussion idle on HB#1+#2; next Phase E
   blockers are HB#5 (Lyapunov metric) and the effective_stiffness law
   decision. Phase E activation requires both blockers locked AND
   HB#1+#2 implementation merged.

Rounds 1–4 of the design lock (+ ack) are MCP id 1466–1476.
