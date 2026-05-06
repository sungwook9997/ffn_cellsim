# V2 Closed-Loop ECM Gate — Hard Blocker #5 Lyapunov-Like Metric (ECM-Orientation Stability Evidence) — Locked Design

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (3-round adversarial lock,
PI id=809 aggressive debate posture applied, PI id=1008/1057/1161
autonomy + visible deliverable focus)
**Source unit**: design-discussion `topic=v2-hard-blocker-5-lyapunov-metric`,
MCP id 1498–1506 (rounds 1–3 + seal ack)
**Brief**: `docs/v2/v2_hard_blocker_5_lyapunov_metric_brief.md` (commit `84b2bcd`)
**Parent locked plan**: `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
(Phase E entry, Hard Blocker #5 — final upstream blocker)
**Upstream sister locks**:
- `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md` (HB#3, traction scattering geometry)
- `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` (HB#4, neutral multiplier sampler)
- `docs/v2/v2_phase_d_no_op_scaffolding_locked.md` (Phase D no-op orchestrator)
- `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md` (HB#1+#2, orientation-only convex update)
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for the HB#5 Sanity Gate doc + code entry.

---

## 0. Scope

This unit locks the **Lyapunov-like metric** that provides the
ECM-orientation stability evidence channel for the closed-loop ECM
gate. It is the **final Phase E upstream hard blocker**; HB#1–#4 are
locked + implemented at this commit.

**What HB#5 is**:
- Pure read-only function on an `ECMSubstrateState` snapshot + a
  `traction_density_xy` field
- Returns a scalar `v_active_um2` (Lyapunov-like energy over **active
  cells only**) + a typed diagnostics record
- Provides two evidence channels: (a) **fixed-target per-step
  contraction** under one HB#1+#2 update with the same traction; (b)
  **varying-target boundedness** `0 ≤ V_active ≤ V_active_max_bound`
  at every measured step

**What HB#5 is NOT**:
- NOT a state-update function (no ECM mutation)
- NOT a full FA→ECM→FA closed-loop stability proof — HB#5 alone
  provides evidence for the **ECM-orientation half** only. Phase E
  composition can use it for bounded-feedback evaluation; full closed
  loop requires HB#4-active design (currently HB#4 locked at neutral
  multipliers — Phase D no-op).
- NOT a measurement of `accumulated_traction_nNs_per_um2` evolution,
  `stiffness_kpa` evolution, `fiber_density` evolution, or
  `ligand_density` evolution — these are NOT updated by HB#1+#2 and
  therefore have no Lyapunov interpretation under this metric.

**Why now**: Phase E activation requires the stability evidence channel
locked before any ECM-mutation can be permitted. HB#1+#2 final PASS
`id=1495` on `f05abef` closed the last upstream blocker. This unit
unblocks HB#5; together with HB#1+#2 + HB#3 + HB#4 + Phase D no-op,
the Phase E composition cycle can begin (a separate unit).

**Hard Rule 11 anchor (measurement-protocol consistency)**: the metric
domain must match the update domain. HB#1+#2 update zero-traction
cells with `w = 0` (target placeholder algebraically unused). HB#5
therefore restricts `V_active` to `active = (traction_norm > 0)`
cells; passive-cell orientation magnitude is reported as a **separate
diagnostic field**, NOT part of `V_active` (Y1, see §2). This is the
direct anti-pattern of the v13 Stage 1a `_build_curvature` Laplacian
off-peak episode (CLAUDE.md Sanity Gate §6).

---

## 1. Final Lock Summary

### Constants (derivation chain explicit)

```python
# Adversarial bound for ||T - n⊗n||²_F under schema invariant
# (T symmetric, |T_ij| ≤ 1, no PSD) is 9.0 — see §3 derivation.
MAX_SQ_FROBENIUS_DIFF_PER_CELL: Final[float] = 9.0

# Energy factor: V = 0.5 · cell_area · Σ ||T - n⊗n||²_F  ⇒ per-cell max = 4.5 · cell_area.
MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL: Final[float] = (
    0.5 * MAX_SQ_FROBENIUS_DIFF_PER_CELL  # = 4.5
)
```

**Hard Rule 10 unit chain** (inline, mandatory per
`rule10_unit_derivation_in_docs` memory):

- `traction_density_xy [nN/μm²]` (from HB#3)
- `cell_area_um2 = ecm.spacing_um² [μm² ]`
- `||T - n⊗n||²_F` is dimensionless (T components dimensionless per
  schema)
- `V_active_um2 = 0.5 · cell_area · Σ_active (dimensionless)` →
  `[μm²]` ✓ (the energy is reported in `μm²`, the natural area unit
  for a per-cell dimensionless integrand on the ECM grid)

**Magic-Number Block compliance**:
- `MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0`: derivable (analytical worst
  case under schema invariant ✓), grid-invariant (no dx/dt dependence
  ✓), not fitted (purely mathematical bound ✓)
- `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 4.5`: derived
  algebraically from `0.5 · 9.0` ✓

### Typed result + diagnostics

```python
@dataclass(frozen=True, slots=True)
class ECMOrientationLyapunovMetricDiagnostics:
    """Typed diagnostics for one HB#5 Lyapunov-like metric snapshot."""
    n_active_cells: int                                         # cells with traction_norm > 0
    n_total_cells: int                                          # = nx * ny
    cell_area_um2: float                                        # = ecm.spacing_um²
    active_area_um2: float                                      # = n_active_cells · cell_area
    grid_area_um2: float                                        # = n_total_cells · cell_area
    v_active_max_bound_um2: float                               # = 4.5 · active_area_um2
    v_active_per_cell_max_um2: float                            # max over active cells
    v_active_per_cell_mean_um2: float                           # mean over active cells
    passive_orientation_magnitude_um2: float                    # diagnostic-only, NOT part of V
    max_traction_norm_nN_per_um2: float                         # Rule 10 bridge


@dataclass(frozen=True, slots=True)
class ECMOrientationLyapunovMetricResult:
    """Outputs of one HB#5 metric snapshot."""
    v_active_um2: float                                         # the Lyapunov-like energy
    diagnostics: ECMOrientationLyapunovMetricDiagnostics
```

### Function

```python
def compute_ecm_orientation_lyapunov_metric(
    ecm: ECMSubstrateState,
    traction_density_xy: np.ndarray,                # (nx, ny, 2) nN/μm², HB#3 output
) -> ECMOrientationLyapunovMetricResult:
    # --- Validate (HB#1+#2 sister-pattern; ecm.validate() FIRST per id=1486 silent-heal closure) ---
    ecm.validate()
    _validate_finite_traction(traction_density_xy)
    _validate_traction_shape(traction_density_xy, ecm.grid_shape)

    # --- Active mask (branch-defined, NO epsilon — sister-pattern with HB#1+#2 Y4) ---
    traction_norm = np.linalg.norm(traction_density_xy, axis=-1)
    active = traction_norm > 0.0  # exact componentwise

    # --- Direction + per-cell distance² (active cells only) ---
    cell_area_um2 = float(ecm.spacing_um) ** 2
    n_xy = np.zeros_like(traction_density_xy, dtype=np.float64)
    n_xy[active] = traction_density_xy[active] / traction_norm[active, None]
    T_target = np.einsum("...i,...j->...ij", n_xy, n_xy)
    diff = ecm.orientation_tensor - T_target
    sq_diff_per_cell = np.sum(diff ** 2, axis=(-2, -1))  # (nx, ny), Frobenius²

    # --- V_active (active cells only — Hard Rule 11 measurement-protocol consistency) ---
    if active.any():
        v_active = 0.5 * cell_area_um2 * float(sq_diff_per_cell[active].sum())
        per_cell_active = sq_diff_per_cell[active]
        v_per_cell_max = 0.5 * cell_area_um2 * float(per_cell_active.max())
        v_per_cell_mean = 0.5 * cell_area_um2 * float(per_cell_active.mean())
    else:
        v_active = 0.0
        v_per_cell_max = 0.0
        v_per_cell_mean = 0.0

    # --- Passive diagnostic (NOT part of V — diagnostic-only) ---
    passive = ~active
    if passive.any():
        passive_mag_sq = float(np.sum(ecm.orientation_tensor[passive] ** 2))
        passive_orientation_magnitude = 0.5 * cell_area_um2 * passive_mag_sq
    else:
        passive_orientation_magnitude = 0.0

    # --- Bound + diagnostics ---
    n_active = int(active.sum())
    n_total = int(np.prod(ecm.grid_shape))
    active_area = n_active * cell_area_um2
    grid_area = n_total * cell_area_um2
    v_active_max_bound = MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL * active_area

    diagnostics = ECMOrientationLyapunovMetricDiagnostics(
        n_active_cells=n_active,
        n_total_cells=n_total,
        cell_area_um2=cell_area_um2,
        active_area_um2=active_area,
        grid_area_um2=grid_area,
        v_active_max_bound_um2=v_active_max_bound,
        v_active_per_cell_max_um2=v_per_cell_max,
        v_active_per_cell_mean_um2=v_per_cell_mean,
        passive_orientation_magnitude_um2=passive_orientation_magnitude,
        max_traction_norm_nN_per_um2=float(traction_norm.max()),
    )

    return ECMOrientationLyapunovMetricResult(
        v_active_um2=v_active,
        diagnostics=diagnostics,
    )
```

### Forbidden in HB#5

- Inclusion of zero-traction cells with `T_target = 0` placeholder
  in the V sum (Y1 — orientation-magnitude-penalty smuggling via
  placeholder; measurement domain must match update domain)
- Reading or referencing `accumulated_traction_nNs_per_um2`,
  `stiffness_kpa`, `fiber_density`, or `ligand_density` (Y7 — the
  metric is for orientation only; AST-walk meta-test enforces this,
  see §4)
- Mutation of any ECM array field (function is pure read-only; no
  copies needed because no return-value ECM)
- Pluggable strategy / composite tuple return / dict diagnostics
  (Y5 — single typed result + typed diagnostics, B1 sister-pattern)
- ε / epsilon as physical or numerical regularizer in the active mask
  (Y7 — branch-defined `traction_norm > 0.0` only)
- Hard cap or clamp on `v_active_um2` (the bound is a diagnostic
  field for testing, not an in-line clamp)
- "Positive-definite over full grid" claim (Y2 — actual property is
  positive-semidefinite relative to the active-input measurement
  protocol, with equilibrium set including all-passive states)
- "Cumulative bound = Lyapunov decrease" wording (Y3 — boundedness and
  per-step contraction are **separate evidence channels** and must be
  tested separately)
- "Full closed-loop FA→ECM→FA stability satisfied by HB#5" claim
  (Y4 — HB#5 is ECM-response-half evidence only; full loop requires
  HB#4-active design)

---

## 2. Reasoned-acceptance trace (Y1–Y11)

The lock converged after 3 rounds. Each Y is a Claude concession with
the round in which it was accepted, preserving the adversarial-debate
audit trail (PI id=809 posture).

**Y1 (round 2, accepted Codex C1 = active-cells-only V; passive
diagnostic-only)**:
- Claude opening lean: `V = 0.5·Σ‖T-T_target‖²_F·cell_area` over the
  full grid, with `T_target = 0` for zero-traction cells.
- Codex catch: HB#1+#2 update has `w = 0` for zero-traction cells, so
  the target placeholder is algebraically unused; including those
  cells in V with `T_target = 0` retroactively assigns physical meaning
  to a placeholder, smuggling an orientation-magnitude penalty into
  the Lyapunov energy. This is the v13 Stage 1a `_build_curvature`
  off-peak anti-pattern (analytical proof valid in active region;
  measurement extends past active region).
- Resolution: `V_active = 0.5 · cell_area · Σ_active ||T - n⊗n||²_F`
  over `active = (traction_norm > 0)` only. Passive orientation
  magnitude reported as a separate diagnostic field, NOT part of V.

**Y2 (round 2, accepted Codex C2 = positive semidefinite, not
positive definite)**:
- Brief claimed "V is positive-definite over full grid orientation
  states".
- Codex catch: under the active-cells-only restriction, the
  equilibrium set includes all passive cells (zero traction → zero
  update) AND active cells exactly at target. The honest control
  property is **nonnegative + zero on the equilibrium set**
  (positive-semidefinite relative to the active-input measurement
  protocol).
- Resolution: lock §0 + §3 use the semidefinite wording.

**Y3 (round 2, accepted Codex C3 = contraction vs boundedness as
separate evidence channels)**:
- Brief conflated "cumulative bound" with "Lyapunov-like decrease".
- Codex correction: `V ≤ V_max` is **boundedness evidence**, not
  Lyapunov decrease. Lyapunov-like decrease is per-step under fixed
  target. Lock both channels separately:
  - **fixed-target contraction**: `V_active_after ≤ V_active_before`
    under one HB#1+#2 update with the same traction field
  - **varying-target boundedness**: `0 ≤ V_active ≤ V_active_max_bound`
    at every measured step
- Resolution: lock §1 (forbidden list) + §4 (separate tests 5 and 6)
  preserve the Lyapunov-like part without overclaiming closed-loop
  monotonicity.

**Y4 (round 2, accepted Codex C4 = HB#5 is ECM-response-half evidence,
NOT full closed-loop satisfaction)**:
- Brief was loose on the closed-loop coupling depth.
- Codex correction: HB#5 alone provides evidence for the **ECM-response
  half**. Phase E composition can use it for bounded-feedback
  evaluation, but HB#5 does NOT prove full FA→ECM→FA stability while
  HB#4-active is absent.
- Resolution: lock §0 hard boundary + §1 forbidden-list item.

**Y5 (round 2, accepted Codex C5 = single typed result + 10-field
typed diagnostics; no pluggable strategy)**:
- Codex specified the dataclass shape; Claude accepted.
- Diagnostics rationale:
  - `n_active_cells`/`n_total_cells`/`cell_area_um2`/`active_area_um2`/
    `grid_area_um2`: counting + area accounting for V interpretation
  - `v_active_max_bound_um2`: schema-tight bound (= 4.5 · active_area)
    for varying-target boundedness test
  - `v_active_per_cell_max_um2`/`v_active_per_cell_mean_um2`: per-cell
    summaries for failure diagnosis
  - `passive_orientation_magnitude_um2`: explicit reporting that
    passive-cell orientation magnitude is NOT in V (Y1 forward guard)
  - `max_traction_norm_nN_per_um2`: Rule 10 bridge to stimulus scale

**Y6 reserved** (no Y6 in this unit's trace; Codex's challenge
numbering C6 maps to round-2 push-back, not a Y; see §1 derivation
correction).

**Y7 (round 2, accepted Codex C7 = pure-function validation +
no-other-field-reference meta-test)**:
- Codex specified: `ecm.validate()` at entry, finite + shape traction
  validation, branch exact (no epsilon), no copies needed (read-only),
  AST-walk meta-test that the function body does NOT reference
  `accumulated_traction_nNs_per_um2`, `stiffness_kpa`, `fiber_density`,
  or `ligand_density`.
- Resolution: §1 validation block + §4 tests 14/15.

**Y8 (round 3, accepted Codex C8 = two-constant form for derivation
traceability)**:
- Claude round 2 pseudocode named `MAX_SQ_FROBENIUS_DIFF_PER_CELL =
  4.5`, which was wrong (the squared Frobenius diff max is 9.0; 4.5
  is the post-`0.5·||...||²` energy factor).
- Codex correction: two constants — `MAX_SQ_FROBENIUS_DIFF_PER_CELL =
  9.0` and `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 0.5 *
  MAX_SQ_FROBENIUS_DIFF_PER_CELL` (= 4.5) — makes the derivation chain
  9 → 4.5 explicit in code itself.
- Resolution: §1 constants + §4 test 4 dual-assertion (raw Frobenius
  9 + energy 4.5).

**Y9 (round 3, accepted Codex Q5 answer = (5a) `acs/v2/dynamics/`
sister-pattern)**:
- Open question: pure read-only metric belongs in `dynamics/` (sister
  with `ecm_constitutive_response.py`), or a new `diagnostics/`
  directory, or a new `metrics/` directory?
- Resolution: (5a) `acs/v2/dynamics/ecm_lyapunov_metric.py`. The
  function consumes HB#3 traction + HB#1+#2 orientation semantics, so
  it belongs beside `ecm_constitutive_response.py`. Function name
  `compute_ecm_orientation_lyapunov_metric` is enough to signal
  read-only behavior (no `step_*` prefix).
- Note: Codex initially mentioned an "existing `acs/v2/metrics/`
  package" which Claude verified does not exist as a directory; Codex
  corrected — `acs/v2/metrics.py` exists as a **file** (pure
  measurement-registry module for v2 output metrics, NOT a Phase E
  dynamics-gate evidence home). Either way, the (5a) decision stands
  on the sister-pattern argument.

**Y10 (round 3, accepted Codex test 4 refinement = dual 9 + 4.5
assertion)**:
- Test 4 should assert both `per_cell_max ≈ 4.5 · cell_area` and the
  per-cell raw squared Frobenius diff = 9.0 for the adversarial cell.
  This guards against the 9-vs-4.5 naming mixup reappearing later.
- Resolution: §4 test 4 split into two assertions.

**Y11 (round 3, accepted Codex tests 14/15 refinement = AST-walk on
function body, not module-level regex)**:
- Naive regex over `inspect.getsource(module)` would catch references
  in module docstrings or comments (e.g., a forbidden-field name
  appearing in a docstring "this function does NOT consume X"), causing
  false failures.
- Codex specified: AST walk on
  `inspect.getsource(compute_ecm_orientation_lyapunov_metric)` (function
  body only), wrapped with `textwrap.dedent(...)` for robustness
  against decorators/formatting changes.
- Resolution: §4 tests 14/15 use AST `ast.Attribute` + `ast.Name` walk
  rather than regex; lock-document forbidden-field-name wording is also
  scoped (this lock §1 docstrings/comments do not embed forbidden field
  names where they would propagate into the source file).

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

1. **Units**:
   - `traction_density_xy`: nN/μm² (HB#3 output)
   - `traction_norm`: nN/μm² (componentwise vector norm)
   - `cell_area_um2`: μm² (= `ecm.spacing_um²`)
   - `||T - n⊗n||²_F`: dimensionless (T components dimensionless per
     schema, n unit vector)
   - `V_active_um2`: μm² (= `0.5 · cell_area · Σ_active dimensionless`)
   - `MAX_SQ_FROBENIUS_DIFF_PER_CELL`: dimensionless (= 9.0)
   - `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL`: dimensionless (= 4.5)
   - `v_active_max_bound_um2`: μm² (= 4.5 · active_area_um2)
2. **Boundary**:
   - All-zero traction: `active.any() == False` → `v_active = 0.0`,
     `v_per_cell_max = 0.0`, `v_per_cell_mean = 0.0`,
     `passive_orientation_magnitude = 0.5 · cell_area · ||T||²_F`
     (diagnostic-only)
   - All-active traction: `passive.any() == False` →
     `passive_orientation_magnitude = 0.0`
   - Single-cell ECM (`grid_shape = (1, 1)`): all bounds + counts well
     defined
   - At-target orientation: `T == n⊗n` everywhere active →
     `v_active = 0.0`
   - Adversarial schema-saturating: `T = ((-1, -1), (-1, -1))`
     everywhere, `traction = (1/√2, 1/√2) · const > 0` →
     per-cell sq diff = 9.0, `v_active = 4.5 · active_area`
     (saturation of bound)
3. **Conservation**: function is pure read-only — nothing to conserve
   in the function itself. The Lyapunov-like contraction property is a
   **test-time** assertion (V_after ≤ V_before under fixed target +
   one HB#1+#2 step), NOT an in-function invariant.
4. **Numerical**:
   - `np.float64` enforced everywhere (matches HB#1+#2 + ECM schema)
   - `np.einsum("...i,...j->...ij", n_xy, n_xy)` is exact algebra; no
     accumulation order issues for the small per-cell tensor outer
     product
   - Frobenius² computed via `np.sum(diff**2, axis=(-2, -1))` is
     componentwise positive; sum order is stable for typical grid
     sizes (no Kahan summation required at this scale)
   - `cell_area_um2 = float(ecm.spacing_um) ** 2` — exact float
     squaring, no rounding artifact
5. **Sign**: `V_active ≥ 0` always (sum of squares times nonneg
   factors). Per-cell components in `sq_diff_per_cell ≥ 0`. Bound
   `v_active_max_bound_um2 ≥ 0` always (4.5 nonneg, area nonneg).
   Passive diagnostic ≥ 0.
6. **Measurement-protocol** (Hard Rule 11, the central anchor):
   - The HB#1+#2 update domain is **active cells only** (zero-traction
     cells have `w = 0`, target placeholder algebraically unused)
   - Therefore the Lyapunov metric MUST also be restricted to active
     cells (Y1)
   - Passive-cell orientation magnitude is reported as a **separate
     diagnostic** so the field is visible without polluting V
   - Adversarial bound derivation (= the §1 constant 9.0):
     - `||T - T_target||²_F = (T_xx - cos²θ)² + 2(T_xy - cosθ·sinθ)² +
       (T_yy - sin²θ)²` for unit `n = (cos θ, sin θ)`
     - Adversarial T (componentwise extreme against target): `T_xx =
       T_yy = -1`, `T_xy = -sign(sin 2θ)`
     - Max sum at `θ = π/4`: `(3/2)² · 4 = 9` (verified via Frobenius
       identity: `||T||²_F + 1 − 2 · n^T T n = 4 + 1 + 4 = 9`)
     - Therefore `MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0` is
       schema-tight
   - The contraction property is **per-step under fixed target**:
     `T_new − T_target = (1−w)(T_old − T_target)` ⇒ `‖T_new − T_target‖
     ≤ ‖T_old − T_target‖` in any norm (Frobenius default)
   - Global monotonicity over time is NOT claimed (traction direction
     can change between steps); varying-target boundedness is the
     companion evidence channel

---

## 4. Test Catalog (~16 tests)

1. `test_lyapunov_metric_zero_traction_returns_zero` — all-zero
   traction → `V_active = 0.0`; `passive_orientation_magnitude_um2 =
   0.5 · grid_area · ||T||²_F`.
2. `test_lyapunov_metric_no_active_cells_no_passive_penalty_in_v` —
   confirm V_active does NOT include passive cells (set T to nonzero
   in passive cells; verify V_active still 0; passive diagnostic
   reflects the magnitude).
3. `test_lyapunov_metric_at_target_returns_zero` — orientation tensor
   equals `n⊗n` for the prescribed traction direction → `V_active =
   0.0` exactly.
4. `test_lyapunov_metric_v_active_within_bound_4_5_times_active_area` —
   adversarial setup: T = `((-1,-1),(-1,-1))` everywhere, traction
   along `(1/√2, 1/√2)` (unit vector at θ = π/4). **Dual assertion**:
   - per-cell raw squared Frobenius diff = 9.0 (matches
     `MAX_SQ_FROBENIUS_DIFF_PER_CELL`)
   - `v_active_um2 ≈ 4.5 · active_area_um2` (matches
     `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL · active_area`)
5. `test_lyapunov_metric_fixed_target_contraction_under_hb12_step` —
   apply one `step_ecm_orientation_response` with the same traction
   field; assert `V_active_after ≤ V_active_before` (per-step
   Lyapunov-like contraction under fixed target).
6. `test_lyapunov_metric_varying_target_boundedness` — apply many steps
   with rotating traction direction; assert `0 ≤ V_active ≤
   V_active_max_bound_um2` at every measured step.
7. `test_lyapunov_metric_pure_no_mutation_of_ecm` — ecm
   orientation/stiffness/fiber/ligand/accumulated bytewise unchanged
   after metric call (no copy needed since read-only).
8. `test_lyapunov_metric_traction_nonfinite_raises`.
9. `test_lyapunov_metric_traction_shape_mismatch_raises`.
10. `test_lyapunov_metric_invalid_ecm_raises_at_validate` —
    schema-invalid input ECM raises (silent-heal path closure, HB#1+#2
    `id=1486` sister-pattern).
11. `test_lyapunov_metric_diagnostics_typed_dataclass` —
    `isinstance(result.diagnostics,
    ECMOrientationLyapunovMetricDiagnostics)`,
    `with pytest.raises(TypeError):
    result.diagnostics["v_active_max_bound_um2"]`
    (B1 sister-pattern meta-test).
12. `test_lyapunov_metric_passive_diagnostic_separate_from_v_active` —
    set T to large magnitude in passive cells; verify `V_active`
    unchanged; `passive_orientation_magnitude_um2` reflects the
    magnitude (Y1 forward guard).
13. `test_lyapunov_metric_exports_through_both_init` — symbols
    `ECMOrientationLyapunovMetricDiagnostics`,
    `ECMOrientationLyapunovMetricResult`,
    `compute_ecm_orientation_lyapunov_metric`,
    `MAX_SQ_FROBENIUS_DIFF_PER_CELL`,
    `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL` importable from BOTH
    `acs.v2.dynamics` AND `acs.v2`.
14. `test_lyapunov_metric_does_not_reference_other_ecm_fields` —
    AST-walk meta-test on `inspect.getsource(
    compute_ecm_orientation_lyapunov_metric)` wrapped with
    `textwrap.dedent(...)`; assert none of `accumulated_traction_nNs_per_um2`,
    `stiffness_kpa`, `fiber_density`, `ligand_density` appear as
    `ast.Attribute.attr` or `ast.Name.id` in the function body.
    Implementation pattern:
    ```python
    import ast, inspect, textwrap
    src = textwrap.dedent(inspect.getsource(compute_ecm_orientation_lyapunov_metric))
    tree = ast.parse(src)
    referenced = (
        {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)} |
        {n.id   for n in ast.walk(tree) if isinstance(n, ast.Name)}
    )
    forbidden = {
        "accumulated_traction_nNs_per_um2",
        "stiffness_kpa",
        "fiber_density",
        "ligand_density",
    }
    assert not (forbidden & referenced)
    ```
15. (merged into 14 — single AST walk covers all 4 forbidden names).
16. `test_lyapunov_metric_active_area_consistent_with_n_active` —
    `result.diagnostics.active_area_um2 ==
    result.diagnostics.n_active_cells *
    result.diagnostics.cell_area_um2` (sanity).
17. `test_lyapunov_metric_max_bound_consistent_with_active_area` —
    `result.diagnostics.v_active_max_bound_um2 == 4.5 *
    result.diagnostics.active_area_um2` (Y8 derivation chain check).

(15 active tests; tests 14/15 merged.)

---

## 5. Files

- `docs/v2/v2_hard_blocker_5_lyapunov_metric_brief.md` (existing,
  opening brief, commit `84b2bcd`)
- `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md` (this file,
  source of truth)
- `acs/v2/dynamics/ecm_lyapunov_metric.py` (new):
  - `MAX_SQ_FROBENIUS_DIFF_PER_CELL`,
    `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL` module constants
  - `ECMOrientationLyapunovMetricDiagnostics` dataclass (frozen,
    slots, 10 fields)
  - `ECMOrientationLyapunovMetricResult` dataclass (frozen, slots)
  - `compute_ecm_orientation_lyapunov_metric` function
  - Validation helpers can be reused from the HB#1+#2 module if
    available (`_validate_finite_traction`, `_validate_traction_shape`)
    or duplicated minimally to keep modules independent — impl-work
    decision per Sanity Gate
- `acs/v2/dynamics/__init__.py` (5 new exports + alphabetical reorder
  of `__all__`)
- `acs/v2/__init__.py` (mirror exports)
- `tests/v2/test_v2_ecm_lyapunov_metric.py` (new, ~15 tests per §4)

---

## 6. References

- Brief (opening position):
  `docs/v2/v2_hard_blocker_5_lyapunov_metric_brief.md` (commit `84b2bcd`)
- Parent locked plan:
  `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
- Sister locks:
  - `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2/v2_phase_d_no_op_scaffolding_locked.md`
  - `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`
  - `docs/v2/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1)
- ECM substrate schema: `acs/v2/ecm_substrate.py:42`
  (`spacing_um`, componentwise `|T_ij| ≤ 1`, no PSD)
- HB#1+#2 implementation (sister precedent for validation pattern):
  `acs/v2/dynamics/ecm_constitutive_response.py`
- v2 measurement registry (separate concern, NOT this unit's home):
  `acs/v2/metrics.py` (file, not directory)
- CLAUDE.md Sanity Gate §6 (measurement-protocol consistency, v13
  Stage 1a `_build_curvature` off-peak anti-pattern)
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

1. impl Claude writes HB#5 Sanity Gate doc
   (`docs/v2/v2_hard_blocker_5_lyapunov_metric_sanity_gate.md`) from this
   lock — 6 sections per §3, full bound derivation 9 → 4.5 inline,
   Hard Rule 11 measurement-protocol section as the central anchor.
2. impl Codex review (6 focus per Codex `id=1506` final SEAL):
   - active-cells-only V (Hard Rule 11) — passive cells in separate
     diagnostic, NOT in V
   - two-constant form `MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0` +
     `MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL = 4.5` (derivation
     traceability)
   - branch-defined active mask (no epsilon)
   - `ecm.validate()` at entry; pure read-only (no mutation, no
     copies)
   - AST-walk meta-test scoped to function body via
     `textwrap.dedent(inspect.getsource(...))`
   - 5 new symbols exported through both `acs.v2.dynamics` and
     `acs.v2`
3. On Sanity Gate PASS: HB#5 code commit
   (`acs/v2/dynamics/ecm_lyapunov_metric.py` new + 2 export updates +
   1 new test file).
4. After commit: design-discussion idle on HB#5; **all 5 Phase E
   upstream hard blockers (HB#1+#2+#3+#4+#5) locked + impl**. Phase E
   composition becomes the next blocker (separate cycle: combining
   Phase D no-op wrapper + HB#1+#2 active law + HB#5 metric into a
   closed-loop ECM gate Phase E module). The effective_stiffness law
   decision is also still pending — to be picked up as a separate unit
   if/when consumers need it.

Rounds 1–3 of the design lock (+ ack) are MCP id 1498–1506.
