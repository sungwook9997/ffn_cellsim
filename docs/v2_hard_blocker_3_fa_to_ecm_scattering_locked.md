# V2 Phase 1 Hard Blocker #3 — FA→ECM Scattering Geometry (LOCKED)

**Date**: 2026-05-04 KST
**Authors**: Claude + Codex design-discussion (4-round adversarial lock,
PI id=809 aggressive debate posture, PI id=1008/1057 autonomy)
**Source unit**: design-discussion `topic=v2-layer-2-hard-blocker-3-fa-to-ecm-scattering`,
MCP id 1331-1342
**Parent unit**: closed-loop ECM gate phased plan
(`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`), Hard Blocker #3.
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for Phase D no-op scaffolding (with Hard Blocker #4 also locked).

---

## 0. Scope

This document locks the FA→ECM scattering geometry interface for
closed-loop ECM Phase D entry. Scope is **interface lock only** — Phase
D no-op scaffolding can use this scattering function with default
identity behavior. Active constitutive response (Phase E) consumes the
output but is locked separately under Hard Blocker #1.

The scatter function is **interface-only for #3**. It produces vector
traction density. It does NOT update ECM state, does NOT call the
accumulator, does NOT decide response modality.

---

## 1. Final Scope LOCK

### Function signature

```python
def scatter_fa_traction_to_ecm_bilinear(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
) -> np.ndarray:  # shape (nx, ny, 2), nN/μm²
    """Bilinear single-point scatter of per-FA cell-on-substrate traction
    onto ECM grid as a vector traction-density field.

    Pure function; does not mutate ecm or call accumulate_prescribed_traction.

    Coordinate convention: cell-centered. Grid index (i,j) corresponds
    to center at origin + ((i+0.5)*dx, (j+0.5)*dy).

    Out-of-grid policy: raise FAToECMScatteringError(
        failure_kind="fa_position_outside_ecm_grid", ...
    ) if any FA position is outside physical footprint
    [origin_x, origin_x + nx*dx] x [origin_y, origin_y + ny*dy].

    Boundary half-cell uses truncated stencil with renormalized weights
    (preserves component-wise vector conservation at boundaries).

    Single-point scatter; FA footprint convolution deferred to Phase E.
    """
```

### Inputs
- `adhesions: tuple[FocalAdhesionState, ...]` — each provides
  `position_um_xy` (μm) + `traction_force_nN_xy` (nN, cell-on-substrate
  per 6.3a sign convention)
- `ecm: ECMSubstrateState` — provides `origin_um_xy`, `spacing_um`
  (`dx == dy`), `grid_shape (nx, ny)`

### Output
- `np.ndarray` shape `(nx, ny, 2)`, dtype `float64`, units `nN/μm²`
  vector traction density field
- Component 0 = x-component traction density, component 1 = y-component

### Pure function constraints
- No `ecm` mutation (input ECM unchanged after call)
- No `accumulate_prescribed_traction` call
- No response field update (`stiffness_kpa`, `fiber_density`,
  `orientation_tensor` untouched)
- No scalarization (output is vector field, NOT magnitude/norm)
- No side effects beyond returning the array

---

## 2. Coordinate Convention (LOCKED)

### Cell centers
- `x_i = origin_x + (i + 0.5) * dx`, for `i = 0, 1, ..., nx-1`
- `y_j = origin_y + (j + 0.5) * dy`, for `j = 0, 1, ..., ny-1`
- `dx == dy` per existing ECM schema (uniform spacing)

### Physical footprint (inclusive boundary)
- `[origin_x, origin_x + nx * dx]` × `[origin_y, origin_y + ny * dy]`
- FAs at the boundary inclusive (boundary positions accepted)
- FAs outside this footprint: raise (see §4)

### Boundary half-cell regions
- x: `[origin_x, origin_x + 0.5*dx)` and `(origin_x + (nx-0.5)*dx, origin_x + nx*dx]`
- y: same pattern
- These regions use truncated stencil with renormalized weights (see §3)

### Numerical tie-break tolerance
- `_FA_POSITION_BOUNDARY_TOL_UM = 1e-12 * max(nx*dx, ny*dy)` for
  floating-point boundary inclusion check
- Documented as numerical tie-break, NOT physics; passes
  Magic-Number Block (derivable from float64 precision, grid-invariant
  in relative terms, not fitted)

---

## 3. Bilinear Stencil + Boundary Renormalization

### Interior bilinear (FA strictly inside center-bracket region)

For FA at `(x_fa, y_fa)`:

1. Find bracketing center indices:
   `i = floor((x_fa - origin_x) / dx - 0.5)`,
   `j = floor((y_fa - origin_y) / dy - 0.5)`
   Bracket: centers `(i, j), (i+1, j), (i, j+1), (i+1, j+1)`
2. Compute fractional offsets:
   `fx = (x_fa - x_i) / dx`, `fy = (y_fa - y_j) / dy` (in [0, 1])
3. Bilinear weights:
   - `w_00 = (1 - fx) * (1 - fy)` for `(i, j)`
   - `w_10 = fx * (1 - fy)` for `(i+1, j)`
   - `w_01 = (1 - fx) * fy` for `(i, j+1)`
   - `w_11 = fx * fy` for `(i+1, j+1)`
4. Deposit to grid:
   `traction_density_xy[k, l, :] += w * traction_force_nN_xy / cell_area_um2`
   for each of the 4 (k, l) ∈ stencil with corresponding weight

`cell_area_um2 = dx * dy`

### Boundary half-cell renormalization (FA in boundary half-cell)

For FA in boundary region (e.g., `x_fa < origin_x + 0.5*dx`), some of the
bracketing centers are outside the grid (`i < 0` or `i+1 >= nx`).

Algorithm:
1. Compute raw bilinear weights as in interior case (some weights
   reference centers with `i < 0` or `i >= nx`)
2. Mask weights for valid in-grid centers only:
   `valid_centers = {(k, l) : 0 <= k < nx and 0 <= l < ny}`
   `w_valid_kl = w_kl if (k, l) in valid_centers else 0`
3. Renormalize: `w_normalized_kl = w_valid_kl / sum(w_valid)`
4. Deposit `w_normalized_kl * traction_force_nN_xy / cell_area_um2`
   to valid centers

This preserves total weight = 1.0 (vector conservation), distributing
the FA's full traction across fewer cells. Component-wise conservation
holds within float64 epsilon.

### Edge cases
- **Exact cell center**: `fx = 0` or `fy = 0` → weight 1.0 to that center,
  others 0 (degenerate bilinear)
- **Exact midway between centers**: `fx = 0.5` → 0.5/0.5 split along x
- **Exact 4-center crossing**: `fx = 0.5, fy = 0.5` → 0.25 each
- **Exact physical corner** (e.g., `(origin_x, origin_y)`): boundary
  half-cell truncation gives 100% to nearest corner center
- **Exact physical edge midpoint** (e.g., `(origin_x, origin_y + ny*dy/2)`):
  truncated stencil along edge

---

## 4. Out-of-Grid Policy — Raise (No Clamp)

### Validation
For each FA:
1. Check `position_um_xy` is finite (NaN, ±inf rejected)
2. Check `origin_x - tol <= x_fa <= origin_x + nx*dx + tol`
3. Check `origin_y - tol <= y_fa <= origin_y + ny*dy + tol`
   where `tol = _FA_POSITION_BOUNDARY_TOL_UM`
4. If violated: `raise FAToECMScatteringError(
   failure_kind="fa_position_outside_ecm_grid",
   message=f"FA {adhesion.adhesion_id} at {position_um_xy} is outside "
           f"ECM physical footprint [{origin_x}, {origin_x + nx*dx}] × "
           f"[{origin_y}, {origin_y + ny*dy}]"
   )`

### No clamp
- Silent clamping to nearest cell hides setup errors and breaks
  conservation
- Caller responsibility to ensure FA positions are within ECM domain
  (e.g., via cell-cluster initialization)

---

## 5. Conservation (Numerical Lock)

### Component-wise conservation
For all FA list inputs:
```
sum over grid (i, j): traction_density_xy[i, j, :] * cell_area_um2
== sum over FA: traction_force_nN_xy
```

### Numerical tolerance
- Float64 accumulation, expected error `O(N_FA * eps_float64 * |max FA traction|)`
- Test tolerance: `1e-9 * max(|component sums|)` (loose enough for
  N_FA up to ~10^4, tight enough to catch real bugs)

---

## 6. Sanity Gate 6 Items

1. **Units**: nN traction (FA), `nN/μm²` density (output), μm position,
   μm² cell area. Single Rule 10 chain `[nN] / [μm²] = [nN/μm²]`.
2. **Boundary**: out-of-grid raise (specific failure_kind); boundary
   half-cell renormalized stencil; empty FA list → zero field.
   Numerical tie-break `_FA_POSITION_BOUNDARY_TOL_UM`.
3. **Conservation**: component-wise sum equals FA traction sum (numerical
   tolerance 1e-9 · |max component sum|). Pure function — no state
   mutation possible.
4. **Numerical**: float64 throughout, no auto-shrink, no clipping;
   bilinear weights `w >= 0` always; renormalization preserves
   sum(weights) = 1.0 within float64 epsilon.
5. **Sign**: bilinear weights non-negative (no negative weights);
   deposit preserves input traction sign component-wise (positive Fx
   → positive density at downstream center).
6. **Measurement**: Pure function, no side effects (Hard Rule 11
   protection). Vector output preserves response modality decision for
   downstream Phase E law.

---

## 7. Test Catalog (~16)

- `test_single_fa_at_center_deposits_full_traction_to_one_cell`
- `test_single_fa_at_grid_crossing_splits_25_each`
- `test_single_fa_midway_along_x_splits_50_50`
- `test_single_fa_at_physical_corner_deposits_100_to_corner_center`
- `test_single_fa_on_physical_edge_midpoint_uses_truncated_stencil`
- `test_two_fas_distinct_positions_sum_correctly`
- `test_zero_fa_returns_zero_field`
- `test_conservation_component_wise_with_random_fas`
- `test_out_of_grid_raises_fa_position_outside_ecm_grid`
- `test_inside_boundary_tol_accepted`
- `test_outside_boundary_tol_raises`
- `test_boundary_half_cell_renormalized_weights_sum_to_one`
- `test_pure_function_no_ecm_mutation`
- `test_pure_function_no_accumulator_call_via_mock`
- `test_negative_traction_input_preserved_in_sign_x_and_y`
- `test_returned_shape_is_nx_ny_2_float64`
- `test_non_finite_fa_position_raises`

---

## 8. Phase D Integrator (separate future unit)

The scatter function is the input layer. Phase D no-op integrator is a
separate wrapper that calls scatter and applies the response law (Phase
E). Phase D default behavior is identity (no ECM state change).

```python
def step_fa_to_ecm_response(
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
    scattering_params,  # interface lock #3
    response_params,    # interface lock #1/#2 (Phase E)
) -> ECMSubstrateState:
    """Phase D no-op default integrator. Calls scatter, accumulates per
    response law (Phase E unit decides). Phase D default = identity (no
    state change)."""
    traction_density_xy = scatter_fa_traction_to_ecm_bilinear(adhesions, ecm)
    # Phase E response law plugs in here (no-op identity for Phase D)
    return ecm  # default no-op identity
```

Phase D no-op scaffolding cannot land until Hard Blocker #4 (ECM→FA bias
target) is also locked.

---

## 9. Implementation Files

- `acs/v2/dynamics/fa_to_ecm_scattering.py` (new)
  - `scatter_fa_traction_to_ecm_bilinear(adhesions, ecm) -> np.ndarray`
  - `FAToECMScatteringError(ValueError)` with `failure_kind: Literal["fa_position_outside_ecm_grid", "non_finite_fa_position", "non_finite_fa_traction"]`
  - Module constant `_FA_POSITION_BOUNDARY_TOL_UM = 1e-12 * <derived>`
- `tests/test_v2_fa_to_ecm_scattering.py` (new, ~16 tests)
- `acs/v2/dynamics/__init__.py` (export)
- Optional: `acs/v2/__init__.py` (top-level export)

---

## 10. Sequencing & Cross-Room Dispatch

### #4 immediately after #3 lock
Hard Blocker #4 (ECM→FA bias target) design round can start immediately
after this #3 lock. Both #3 and #4 must be locked before Phase D no-op
scaffolding can land.

### Cross-room dispatch checklist
1. Write this lock artifact (done)
2. Cross-room dispatch to implementation-work room with brief
3. impl Claude writes Sanity Gate doc, then code+tests
4. impl Codex review (5 focus per §6: vector output / coordinate
   convention / boundary renormalization / out-of-grid raise / pure
   function)
5. PASS → commit + cross-room backflow
6. After #3 + #4 both committed: Phase D no-op scaffolding dispatch

### Cadence rule
impl-work in-progress reports per `docs/implementation_workflow.md` §8
(promised cadence even when no change).

### Process improvement (sequencing rule applied)
This dispatch is design-first lock → cross-room dispatch (per PI
id=1088/1089 sequencing fix). impl entry only after cross-room dispatch.

Rounds 1-4 of the Hard Blocker #3 design lock are MCP id 1331-1342.
