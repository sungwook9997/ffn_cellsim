# V2 Phase 1 Hard Blocker #4 — ECM→FA Bias Target (LOCKED)

**Date**: 2026-05-04 KST
**Authors**: Claude + Codex design-discussion (4-round adversarial lock,
PI id=809 aggressive debate posture, PI id=1008/1057 autonomy)
**Source unit**: design-discussion `topic=v2-layer-2-hard-blocker-4-ecm-to-fa-bias-target`,
MCP id 1365-1370
**Parent unit**: closed-loop ECM gate phased plan
(`docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`), Hard Blocker #4.
**Sister lock**: Hard Blocker #3 FA→ECM scattering geometry
(`docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`).
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for Phase D no-op scaffolding (with Hard Blocker #3 also locked).

---

## 0. Scope

This document locks the ECM→FA bias target interface for closed-loop
ECM Phase D entry. Scope is **interface lock only** — Phase D no-op
scaffolding can use this bias function with default neutral
(all 1.0 multipliers) behavior. Active mapping (banded / continuous /
anisotropic) is Phase E and locked separately.

The bias function is **interface + diagnostic sampler in #4**. It
returns per-FA, per-rate multipliers (all 1.0 in Phase D). Sampler
returns local ECM state at each FA position for audit/test only.
Active mapping from sampled values to non-neutral multipliers is Phase E.

**Phase D entry status after #3 + #4 both locked**: UNBLOCKED.

---

## 1. Final Scope LOCK

### Module
`acs/v2/dynamics/ecm_to_fa_bias.py` (new)

### Module constants
```python
RATE_NAMES: Final[tuple[Literal["k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s"], ...]] = (
    "k_maturity_per_s",
    "k_bind_per_s",
    "k_unbind_per_s",
)
```

`multipliers_per_fa[:, k]` corresponds to `RATE_NAMES[k]`. Caller does
not supply or reorder rate names. Result `rate_names` field equals
`RATE_NAMES` constant.

### Dataclasses

```python
@dataclass(frozen=True, slots=True)
class ECMSampledAtFAs:
    """Bilinear-interpolated ECM state at FA positions. Diagnostic-only
    in Phase D; Phase E active mapping consumes these for non-neutral
    multipliers."""
    fa_ids: tuple[str, ...]                            # input order preserved
    stiffness_kpa: np.ndarray                          # (N_FA,) kPa
    fiber_density: np.ndarray                          # (N_FA,) [0,1]
    ligand_density: np.ndarray                         # (N_FA,) [0,1]
    orientation_tensor: np.ndarray                     # (N_FA, 2, 2) traceless symmetric
    accumulated_traction_nNs_per_um2: np.ndarray       # (N_FA,) scalar history (matches ECM schema)


@dataclass(frozen=True, slots=True)
class ECMToFABiasResult:
    fa_ids: tuple[str, ...]                            # input order preserved
    multipliers_per_fa: np.ndarray                     # (N_FA, 3) for RATE_NAMES order
    rate_names: tuple[str, ...]                        # == RATE_NAMES module constant
    sampled_diagnostics: Optional[ECMSampledAtFAs] = None  # always present in Phase D
    diagnostics_dict: dict[str, float | int | str] = field(default_factory=dict)
```

### Functions

```python
def sample_ecm_at_fa_positions(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
) -> ECMSampledAtFAs:
    """Bilinear interpolation read of ECM fields at FA positions.

    Convention matches scatter_fa_traction_to_ecm_bilinear (#3):
    - cell-centered grid: index (i,j) = origin + ((i+0.5)*dx, (j+0.5)*dy)
    - inclusive physical footprint
    - boundary half-cell uses truncated stencil with renormalized weights
    - out-of-grid raises FAToECMBiasError(failure_kind="fa_bias_position_outside_ecm_grid")

    Pure read function. No ECM mutation. Empty FA list returns empty arrays.
    """


def compute_ecm_to_fa_bias_neutral(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
) -> ECMToFABiasResult:
    """Phase D no-op default. Returns 1.0 multipliers for all FAs and all
    RATE_NAMES rates. Includes sampled diagnostics from
    sample_ecm_at_fa_positions. Sampled values do NOT alter multipliers
    in Phase D (effective-stiffness side-door guard).

    Pure function. No FA/ECM mutation. No active mapping.

    Phase E active mapping (banded / continuous / anisotropic) is a
    SEPARATE function `compute_ecm_to_fa_bias_active(...)` with its own
    Sanity Gate and lock — silent activation guard.
    """
```

### Error class

```python
class FAToECMBiasError(ValueError):
    def __init__(self, message: str, *, failure_kind: FailureKind):
        super().__init__(message)
        self.failure_kind = failure_kind


FailureKind = Literal[
    "fa_bias_position_outside_ecm_grid",  # distinct from #3 scatter context
    "non_finite_fa_position",
]
```

### Forbidden in Phase D
- Active non-1.0 multipliers (banded, continuous, anisotropic mapping)
- `traction_scale_nN` mutation (no 6.3a magnitude semantics change)
- FA nucleation
- RNG / stochastic
- State label transitions
- FA / ECM state mutation
- Effective_stiffness helper using geometry fields (parent §3 guard)
- Single-function with optional `mapping_law` param (silent activation
  risk; Phase E uses different function name `compute_ecm_to_fa_bias_active`)
- Vector accumulated_traction sampling (schema is scalar; vector field
  would require ECM schema change + separate lock)
- `empty_adhesions` raise (return empty result instead, consistent with
  #3 zero-FA scatter and 6.3a empty list)
- Caller-supplied rate_names reordering (use `RATE_NAMES` constant)
- Arrays in `diagnostics_dict` (use typed fields)

---

## 2. Coordinate Convention

Identical to #3 scatter:
- Cell centers: `x_i = origin_x + (i + 0.5) * dx`, `y_j = origin_y + (j + 0.5) * dy`
- Inclusive physical footprint: `[origin_x, origin_x + nx*dx] × [origin_y, origin_y + ny*dy]`
- `_FA_POSITION_BOUNDARY_TOL_UM` reused from #3 (same numerical tie-break)
- Boundary half-cell regions use truncated stencil with renormalized weights

This consistency prevents grid hysteresis: #3 writes with bilinear,
#4 reads with bilinear → no asymmetry in Item 5 sensitivity sweep
(closed-loop side, Phase E).

---

## 3. Bilinear Read Stencil

For FA at `(x_fa, y_fa)`:

### Interior (FA strictly inside center-bracket)
```
i = floor((x_fa - origin_x) / dx - 0.5)
j = floor((y_fa - origin_y) / dy - 0.5)

fx = (x_fa - x_i) / dx     # in [0, 1]
fy = (y_fa - y_j) / dy

w_00 = (1 - fx) * (1 - fy)   for (i, j)
w_10 = fx * (1 - fy)         for (i+1, j)
w_01 = (1 - fx) * fy         for (i, j+1)
w_11 = fx * fy               for (i+1, j+1)

sampled_value = sum(w * field[k, l]) for (k, l) in 4 stencil cells
```

For tensor field `orientation_tensor[nx, ny, 2, 2]`, bilinear-interpolate
each component independently → result `(2, 2)` per FA.

### Boundary half-cell (FA in `[origin, origin + 0.5*dx)` etc.)
1. Compute raw bilinear weights
2. Mask invalid centers (i < 0, i >= nx, j < 0, j >= ny)
3. Renormalize valid weights to sum 1.0
4. Read from valid centers only

### Outside footprint
`raise FAToECMBiasError(failure_kind="fa_bias_position_outside_ecm_grid", ...)`

No silent clamp. No nearest-cell fallback.

---

## 4. Diagnostics Dict (locked keys)

```python
diagnostics_dict = {
    "n_adhesions": int,
    "max_multiplier": float,        # = 1.0 in Phase D neutral
    "min_multiplier": float,        # = 1.0 in Phase D neutral
    "sampler_geometry": str,        # "bilinear_cell_centered"
}
```

JSON-friendly, frame_dump round-trip compatible. No arrays.

---

## 5. Sanity Gate 6 Items

1. **Units**: stiffness kPa, density dimensionless [0,1], orientation
   tensor dimensionless [-1,1] component-wise, multiplier dimensionless
   [≥0], accumulated_traction nN·s/μm². No Rule 10 cross-unit
   comparison (read-only sampler).
2. **Boundary**: out-of-grid raise (`fa_bias_position_outside_ecm_grid`);
   boundary half-cell renormalized read; **empty FA list → empty result,
   no raise**; finite FA position guards.
3. **Conservation**: pure read function (no state mutation possible).
   Sampler bilinear preserves total interpolation weight = 1.0 within
   float64 epsilon (boundary half-cell renormalized).
4. **Numerical**: float64 throughout, bilinear weights ≥ 0,
   renormalization preserves sum=1, no auto-shrink.
5. **Sign**: bilinear weights non-negative; sampled values preserve sign
   of underlying field (orientation can be ±); multipliers default 1.0
   (neutral, sign-positive). Multipliers ≥ 0 per RATE_NAMES.
6. **Measurement**: pure function, no side effects, no scalarization,
   neutral default forbids active behavior. Sampled diagnostics for
   audit only (Hard Rule 11 + parent §3 effective_stiffness guard).

---

## 6. Test Catalog (~17)

- `test_neutral_multipliers_all_ones_for_arbitrary_fa_list`
- `test_neutral_multipliers_match_RATE_NAMES_column_order`
- `test_sampler_at_cell_center_returns_field_value_exactly`
- `test_sampler_at_grid_crossing_bilinear_average`
- `test_sampler_boundary_half_cell_renormalized`
- `test_sampler_outside_footprint_raises_fa_bias_position_outside_ecm_grid`
- `test_sampler_input_order_preserved`
- `test_sampler_accumulated_traction_is_scalar_shape_n_fa`
- `test_sampler_orientation_tensor_shape_n_fa_2_2`
- `test_neutral_bias_no_ecm_mutation`
- `test_neutral_bias_no_fa_mutation`
- `test_diagnostics_dict_only_serializable_plain_values`
- `test_diagnostics_dict_keys_are_locked_set` (n_adhesions, max/min_multiplier, sampler_geometry)
- `test_empty_fa_list_returns_empty_result_not_raise`
- `test_sampler_consistency_with_scatter_geometry` — same FA position
  scattered by #3 then sampled by #4 = same convention (cross-blocker
  round-trip)
- `test_returned_shapes_match_locked_signature`
- `test_no_active_law_invoked_in_phase_d`
- `test_non_finite_fa_position_raises`

---

## 7. Phase D Integrator (separate future unit)

Phase D no-op integrator wraps #3 + #4 with identity ECM update:

```python
def step_fa_to_ecm_response(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
) -> tuple[ECMSubstrateState, ECMToFABiasResult]:
    """Phase D no-op integrator. Computes scattered traction (#3) and
    neutral bias multipliers (#4). Default Phase D returns ECM unchanged
    (identity). Phase E active law plugs in here.

    Returns (updated_ecm, bias_result). Phase D: updated_ecm == ecm.
    """
    traction_density_xy = scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)  # #3
    bias_result = compute_ecm_to_fa_bias_neutral(adhesions, ecm)              # #4
    # Phase E response law plugs in here (no-op identity for Phase D)
    return ecm, bias_result
```

Phase D no-op scaffolding (`step_fa_to_ecm_response`) is a separate
unit, dispatched after both #3 and #4 are committed.

---

## 8. Implementation Files

- `acs/v2/dynamics/ecm_to_fa_bias.py` (new)
  - `RATE_NAMES` module constant
  - `ECMSampledAtFAs` dataclass
  - `ECMToFABiasResult` dataclass
  - `sample_ecm_at_fa_positions(adhesions, ecm) -> ECMSampledAtFAs`
  - `compute_ecm_to_fa_bias_neutral(adhesions, ecm) -> ECMToFABiasResult`
  - `FAToECMBiasError(ValueError)` with typed `failure_kind`
- `tests/v2/test_v2_ecm_to_fa_bias.py` (new, ~17 tests)
- `acs/v2/dynamics/__init__.py` (export)
- Optional: `acs/v2/__init__.py` (top-level export)

---

## 9. Sequencing & Cross-Room Dispatch

### Phase D entry — UNBLOCKED after #3 and #4 both committed

#3 (`scatter_fa_traction_to_ecm_bilinear`) lock complete (commit pending in impl-work).
#4 (this lock) design complete (impl pending).

After both committed: Phase D no-op `step_fa_to_ecm_response` dispatchable.

### Cross-room dispatch checklist
1. Write this lock artifact (done)
2. Cross-room dispatch to implementation-work room
3. impl Claude writes Sanity Gate doc, then code+tests
4. impl Codex review (5 focus per §6: RATE_NAMES + column order /
   empty FA empty result / diagnostics_dict serializable /
   accumulated_traction scalar / pure functions only)
5. PASS → commit + cross-room backflow
6. After #4 committed AND #3 also committed: Phase D no-op
   scaffolding dispatch

### Cadence rule
impl-work in-progress reports per `docs/implementation_workflow.md` §8
(promised cadence even when no change).

### Process improvement (sequencing rule applied)
This dispatch is design-first lock → cross-room dispatch (per PI
id=1088/1089 sequencing fix). impl entry only after cross-room dispatch.

Rounds 1-4 of the Hard Blocker #4 design lock are MCP id 1365-1370.
