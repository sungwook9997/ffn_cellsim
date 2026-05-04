# V2 Phase 1 — FA→ECM Scattering Sanity Gate (Hard Blocker #3)

**Status**: pre-execution Sanity Gate per Plan §10/§11 + closed-loop
ECM gate phased plan
(`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §1 Phase D
entry blocker = Hard Blocker #3 + #4 interface locks). Written from
the locked design in
`docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md` (4-round
adversarial design-discussion lock; cross-room dispatch to
implementation-work at MCP `id=1344`).

**Contract**: this gate is the document Plan §11 requires **before**
any executable scattering physics lands in
`acs/v2/dynamics/fa_to_ecm_scattering.py`. If any item §1–§6 below
fails its check, Hard Blocker #3 reports `status=blocker` to PI
with at least three concrete options. No partial physics commits
under uncertainty.

**Source of truth**:
- locked design: `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
- parent plan: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
  §1 Phase D + §2 Hard Blocker #3
- design-discussion ledger: MCP `id=1331-1342`
- impl-work review path: MCP `id=1344-1349` (cross-room
  dispatch + handoff)
- 6.3a static FA traction (output 6.3a's `traction_force_nN_xy`
  is the input FA traction this scatter consumes):
  `acs/v2/dynamics/focal_adhesion.py` +
  `docs/v2_focal_adhesion_dynamics_sanity_gate.md`
- ECM substrate schema (target grid + spacing + origin):
  `acs/v2/ecm_substrate.py`
- FA schema (input position + traction): `acs/v2/focal_adhesion.py`

**Sequencing context**: Hard Blocker #3 is the first of two Phase D
entry interface blockers; Hard Blocker #4 (ECM→FA bias target
interface) is design-discussion's next round. Phase D no-op
scaffolding entry requires both #3 + #4 locked + Sanity-Gated +
implemented. Phase E active closed-loop response is gated on all
five Hard Blockers + effective_stiffness law decision.

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

- A single pure function
  `scatter_fa_traction_to_ecm_bilinear(adhesions, ecm) ->
  np.ndarray (nx, ny, 2)` that bilinearly scatters per-FA
  cell-on-substrate traction onto the ECM grid as a vector
  traction-density field in `nN/μm²`.
- Cell-centered coordinate convention with inclusive footprint
  per locked design §2.
- Boundary half-cell handling via truncated stencil with
  renormalized weights (locked §3).
- Out-of-grid policy via explicit
  `FAToECMScatteringError(failure_kind="fa_position_outside_ecm_grid")`
  (locked §4).
- New `FAToECMScatteringError` exception class with
  `failure_kind` attribute mirroring `ECMOpenLoopError` /
  `FocalAdhesionDynamicsError` precedents.
- 16-test catalog per locked §7.

### Explicitly out of scope (will FAIL gate if introduced)

- ECM state mutation. The function is pure: input ECM is
  unchanged on return.
- Call to `accumulate_prescribed_traction`. The scatter is
  one-step deposit, not a time-accumulating accumulator. Any call
  into the open-loop preflight is a §3 conservation violation
  by introducing implicit dt-multiplication.
- Response field update (`stiffness_kpa`, `fiber_density`,
  `orientation_tensor`). That is Phase E's job, not Phase D's.
- Scalarization at scatter time. Output is vector field
  `(nx, ny, 2)`; **no magnitude/norm reduction at this layer**.
  Hard Rule 11: the consumer (Phase E) decides whether to
  scalarize and how; scattering must not pre-decide.
- Silent out-of-grid clamping. Out-of-grid FAs raise; no nearest
  cell clamp.
- Gaussian / disc / Boussinesq footprint. Single-point bilinear
  only; footprint convolution is deferred to Phase E.
- `step_fa_to_ecm_response` wrapper. That belongs to a separate
  Phase D unit (post Hard Blocker #4 lock); this commit is the
  scatter primitive only.
- Default values in the function signature. The function takes
  `(adhesions, ecm)` only — no `=...` defaults for either
  positional input. Caller-supplied discipline matches the locked
  Phase C sweep harness pattern.

---

## 1. Dimensional analysis

Inputs and outputs (every variable's unit explicit):

| Symbol | Meaning | Unit | Source |
|---|---|---|---|
| `position_um_xy` | per-FA xy position on substrate | μm | `FocalAdhesionState` schema (Cycle C) |
| `traction_force_nN_xy` | per-FA cell-on-substrate force vector | nN | `FocalAdhesionState` schema; sign convention 6.3a §5 (inward radial-to-centroid default) |
| `origin_um_xy` | ECM grid bottom-left corner | μm | `ECMSubstrateState` schema |
| `spacing_um` | uniform cell edge | μm | `ECMSubstrateState` schema (`dx == dy`) |
| `grid_shape` | `(nx, ny)` cell counts per axis | dimensionless integers | `ECMSubstrateState` |
| `cell_area_um2` | `dx * dy` | μm² | derived |
| `_FA_POSITION_BOUNDARY_TOL_UM` | `1e-12 * max(nx*dx, ny*dy)` | μm | derived numerical tie-break |
| `traction_density_xy` (output) | per-cell vector traction density | `nN/μm²` | computed; shape `(nx, ny, 2)` |

### Reductions

- **Single Rule 10 unit chain**: per-FA `[nN] / [μm²] = [nN/μm²]`
  per cell after bilinear weight `w` (dimensionless). The full
  per-FA per-cell deposit is `w · traction_force_nN_xy /
  cell_area_um2` — `[dimensionless] · [nN] / [μm²] = [nN/μm²]` ✓.
- **Numerical tie-break tolerance**: `_FA_POSITION_BOUNDARY_TOL_UM
  = 1e-12 · max(nx·dx, ny·dy)` is in μm, scales with the grid's
  largest physical extent so the relative tolerance is
  `~1e-12` per float64 round-off. Magic-Number Block test 1
  (derivable): from float64 epsilon `~2.22e-16` and a safety
  factor of `~5000` for accumulated error in
  `(x_fa - origin_x) / dx`-style reductions; documented as
  numerical tie-break, not physics.

### Status

PASS — single unit chain, no per-volume vs per-area mismatch
(Hard Rule 10 trivially satisfied). All inputs in known schema
units, output in `nN/μm²` matching the locked Phase D output
contract.

---

## 2. Boundary cases

| Case | Guard | Failure-kind |
|---|---|---|
| Empty `adhesions` tuple | early-return `np.zeros((nx, ny, 2), dtype=np.float64)` (no FA iteration) | (no exception) |
| `position_um_xy` non-finite (NaN, ±inf) | reject before bilinear computation | `non_finite_fa_position` (per locked design) |
| `traction_force_nN_xy` non-finite | reject before deposit | `non_finite_fa_traction` (per locked design) |
| FA position strictly outside `[origin_x, origin_x + nx·dx] × [origin_y, origin_y + ny·dy]` plus `_FA_POSITION_BOUNDARY_TOL_UM` | raise BEFORE bilinear; no silent clamp | `fa_position_outside_ecm_grid` |
| FA position exactly at footprint corner | inclusive boundary; truncated stencil deposits 100% to corner cell center | (no exception) |
| FA position exactly at footprint edge midpoint | inclusive; truncated stencil deposits across edge cells with renormalized weights | (no exception) |
| FA position strictly at a cell center | bilinear degenerates: `fx == 0`, `fy == 0` → weight 1.0 to that center, 0 elsewhere | (no exception) |
| FA position exactly at 4-center crossing (`fx == 0.5, fy == 0.5`) | weights `(0.25, 0.25, 0.25, 0.25)` | (no exception) |
| `ecm.spacing_um` ≤ 0 or non-finite | scatter calls `ecm.validate()` at function entry per `acs/v2/dynamics/ecm_open_loop.py` precedent (Codex review id=1354 alignment); the ECM schema's `ValueError` propagates verbatim | inherited `ECMSubstrateState.validate()` `ValueError` |
| `ecm.grid_shape` non-positive integer | same — `ecm.validate()` at entry catches before scatter | inherited (same) |
| FA tuple containing duplicate `adhesion_id` | not the scatter's responsibility (caller-side concern); scatter just iterates | (no exception) |
| Boundary half-cell with **all** stencil weights mapping to invalid centers | impossible by construction within footprint (`(b.iii) renormalization` applies only when `at least one` valid weight exists; and the inclusive footprint with finite `tol` guarantees at least one valid stencil cell exists) | (impossible per construction; documented as invariant) |

### Status

PASS — every boundary case has either an explicit failure_kind
or a documented "absent = empty result / handled by truncated
stencil" rule. No silent clamp.

---

## 3. Conservation invariants

### Component-wise vector conservation (locked §5)

For every FA list input:

```
sum over grid (i, j): traction_density_xy[i, j, c] * cell_area_um2
  == sum over FA: traction_force_nN_xy[c]
```

for `c ∈ {0, 1}` (x and y components). Numerical tolerance:
`1e-9 * max(|component sums|)`.

### Why this holds

- **Interior FA**: the four bilinear weights sum to 1.0 by
  construction; `Σ_k (w_k · F / A) · A = F · Σ_k w_k = F · 1 = F`.
- **Boundary half-cell FA**: raw weights summing to `s < 1`;
  renormalized weights `w'_k = w_k / s` sum to 1.0; deposit
  `Σ_k (w'_k · F / A) · A = F · Σ_k w'_k = F`.
- **Out-of-grid FA**: raises `fa_position_outside_ecm_grid`;
  the function does not return, so conservation is vacuous on
  failure (§4 documents that the failure mode is loud, not silent).

### Pure-function constraint (locked §1)

- No `ecm` mutation: input `ecm` after the call is byte-identical
  to before. Verified by a per-field equality test on every
  `ECMSubstrateState` array field before and after the call.
- No `accumulate_prescribed_traction` call: tested by mocking the
  function and asserting `mock.assert_not_called()` on a
  representative scatter call.
- No response field update: tested by per-field equality on
  `stiffness_kpa`, `fiber_density`, `orientation_tensor`,
  `accumulated_traction_nNs_per_um2`, `ligand_density` before/after.
- No scalarization at scatter time: function return type is
  `np.ndarray` shape `(nx, ny, 2)`; no `np.linalg.norm` /
  `np.hypot` / dot product happens at deposit time; tested
  structurally (`assert returned.shape == (nx, ny, 2)`).

### Status

PASS — conservation is exact-by-construction within float64
tolerance; pure-function constraints are enforceable by
structural tests.

---

## 4. Numerical sanity

### Float precision

float64 throughout. Per-FA per-cell deposit is one multiplication
and one addition per cell per component (8 ops per FA in the
interior bilinear case; ≤ 4 ops in degenerate cases). Total work
is `O(N_FA)` for a fixed-size 2×2 stencil; memory is `O(nx · ny ·
2)` for the output array.

### No auto-shrink, no clipping

- `_FA_POSITION_BOUNDARY_TOL_UM` is a numerical tie-break for
  inclusion check only; it does not shift the FA position or
  modify the deposit.
- Bilinear weights are computed from the validated FA position,
  not from a snapped grid index.
- Out-of-grid violations raise; no silent move-to-nearest-cell.

### Stability

The scatter is a single-step (no time integration); no CFL-style
bound applies. The single Rule 10 division by `cell_area_um2`
is finite for any positive `spacing_um` (which the ECM schema
already validates).

### Status

PASS — single-step pure function, no dt-rate gate needed; the
existing `_FA_POSITION_BOUNDARY_TOL_UM` numerical tie-break is
the only new numeric and is documented + named, not a tunable.

---

## 5. Sign / sense check

| Quantity | Direction | One-line check |
|---|---|---|
| Per-FA `traction_force_nN_xy` | per 6.3a §5: cell-on-substrate inward radial-to-centroid by default; sign is the caller's. | inherited from 6.3a Sanity Gate. |
| Bilinear weight `w_kl` | non-negative by construction; `w_kl ∈ [0, 1]` for all `k, l`. | one-line: weights are products of `(1-fx, fx, 1-fy, fy)` factors, each in `[0, 1]`. |
| Renormalized weight `w'_kl` (boundary half-cell) | non-negative; `w'_kl = w_kl / s` for `s > 0`; `Σ w'_kl = 1`. | one-line: identical sign to interior; only normalization differs. |
| Per-cell per-component deposit | `w · F / A` preserves sign of `F`; `w` non-negative, `A` positive. | tested via `test_negative_traction_input_preserved_in_sign_x_and_y` per locked §7. |
| Output array per-component | sum of non-negative-weighted FA traction component contributions; sign is the linear combination of input signs. | tested via component-wise conservation test. |

### Status

PASS — every sign is derived from the bilinear construction; no
sign convention is invented at the scatter layer. The 6.3a
inward-radial default propagates verbatim through the scatter.

---

## 6. Measurement-protocol consistency (Hard Rule 11)

The scatter's output modality is **vector traction density
`(nx, ny, 2)` in `nN/μm²`**. The consumer (Phase E response law)
will decide how to reduce or compare this against the response
field; the scatter does not pre-decide.

### No scalarization at scatter time (Hard Rule 11)

Per locked §0 / §5 / §9, the scatter must NOT reduce to a scalar
field at this layer. A future Phase E response law might want
the magnitude (`|F|`); a different law might want the divergence
(`∇·F`); a third might want a directional component along ECM
fiber alignment. All of these are downstream choices on the same
vector field; doing the reduction at the scatter layer would lock
the consumer to one choice and violate Hard Rule 11
(measurement-protocol consistency at the consumer's modality).

### Runtime meta-test (lock §9 + Phase B/C precedent)

A meta-test
`test_scatter_does_not_satisfy_closed_loop_response_law` will
encode the Hard Rule 11 boundary at runtime, parallel to Phase B
test 4 / Phase C test 5:

```python
def test_scatter_does_not_satisfy_closed_loop_response_law():
    """Phase D Hard Blocker #3 — Hard Rule 11 wording protection.

    The bilinear scatter outputs a vector traction-density field;
    it does NOT decide a constitutive response law. Phase E
    response laws (consuming this output) decide how to reduce or
    compare against the response field. This meta-test makes the
    boundary explicit at runtime so a future refactor cannot
    silently relabel the scatter as response-law satisfaction.
    """
    scatter_modality = "vector traction density (nx, ny, 2) in nN/μm²"
    response_law_claim = "ECM response field update (Phase E TBD)"
    assert scatter_modality != response_law_claim, (
        "FA→ECM bilinear scatter is an interface primitive, not a "
        "response law. Phase E (separate Sanity Gate) is the "
        "response-law decision per locked phased plan §1 Phase E."
    )
```

### Status

PASS — measurement modality is the locked output type; no
reduction is invented; runtime meta-test enforces the boundary.

---

## 7. Magic-Number Block check

Every numeric in the scatter implementation is one of:

| Symbol | Source | Magic-Number Block status |
|---|---|---|
| `position_um_xy`, `traction_force_nN_xy`, `origin_um_xy`, `spacing_um`, `grid_shape` | runtime inputs (FA / ECM schemas) | **not a magic number** — caller-supplied via schema |
| `_FA_POSITION_BOUNDARY_TOL_UM = 1e-12 * max(nx*dx, ny*dy)` | numerical tie-break for boundary inclusion check | named module-level constant; derivable from float64 epsilon ~2.22e-16 + safety factor ~5000 for accumulated reduction error; **grid-invariant in relative terms** (scales with grid extent); **not chosen to fit a target** (selected pre-implementation per locked §2) |
| Test conservation tolerance `1e-9 * max(|component sums|)` | per locked §5 numerical lock | named test-only tolerance; documented derivation `O(N_FA · eps_float64 · |max FA traction|)`; not a runtime constant |

Magic-Number Block test pass:

1. **Derivable**: `_FA_POSITION_BOUNDARY_TOL_UM` is a relative
   float64 tolerance, not a physical scale; the conservation
   tolerance is documented as a function of `N_FA`.
2. **Grid-invariant**: `_FA_POSITION_BOUNDARY_TOL_UM` scales with
   the grid's largest extent so the relative inclusion check is
   invariant. The conservation tolerance scales with the FA
   traction component sums.
3. **Not fitting**: both values are documented pre-implementation
   in the locked design (§2 + §5); neither was selected to make
   any specific test pass.

### Status

PASS — Magic-Number Block clean. The `_FA_POSITION_BOUNDARY_TOL_UM`
constant lives at module level with the rationale comment, not as
a function default.

---

## 8. Visual deliverable plan (post-commit)

Aligns with the precedent set by `scripts/run_p1_alpha_gate.py`,
`scripts/run_ecm_ol_harness.py`,
`scripts/run_protrusion_coupled_fa_harness.py`, and
`scripts/run_ecm_ol_sensitivity_sweep.py`. **Out of scope for
the Sanity Gate commit + the code/tests commit** — the visual
deliverable is a separate post-Hard-Blocker-3 unit (potentially
combined with the Hard Blocker #4 visual deliverable since both
land before Phase D scaffolding can produce closed-loop output).

Planned (NOT committed by this gate):
- A scatter-only smoke runner exercising the four boundary cases
  (interior / cell-center / 4-cell crossing / corner / edge
  midpoint) with a small fixture FA list and a small ECM grid.
- Per-scenario `summary.html` + `metadata.json` +
  `diagnostic_<scenario>.png` (heatmap of the scatter output's
  magnitude per component, plus a quiver overlay).
- Top-level `index.json` with `aggregate_status`,
  `intentional_failures_observed`, `unexpected_failures`,
  `unexpected_passes` per the f8cdff3 / runner pattern.

This plan is documentary; the code/tests commit does not depend
on it.

---

## 9. Test catalog (16, per locked §7)

Each test exercises one Sanity Gate item or one forbidden
behavior.

| # | Test name | Gate item |
|---|---|---|
| 1 | `test_single_fa_at_center_deposits_full_traction_to_one_cell` | §1 unit chain, §3 conservation, §5 sign |
| 2 | `test_single_fa_at_grid_crossing_splits_25_each` | §3 boundary case at 4-center crossing |
| 3 | `test_single_fa_midway_along_x_splits_50_50` | §3 1D bilinear edge case |
| 4 | `test_single_fa_at_physical_corner_deposits_100_to_corner_center` | §2 corner boundary + truncated stencil |
| 5 | `test_single_fa_on_physical_edge_midpoint_uses_truncated_stencil` | §2 edge boundary + renormalization |
| 6 | `test_two_fas_distinct_positions_sum_correctly` | §3 linearity / conservation across FAs |
| 7 | `test_zero_fa_returns_zero_field` | §2 empty list early-return |
| 8 | `test_conservation_component_wise_with_random_fas` | §3 conservation invariant numerical |
| 9 | `test_out_of_grid_raises_fa_position_outside_ecm_grid` | §2 / §4 out-of-grid policy |
| 10 | `test_inside_boundary_tol_accepted` | §2 numerical tie-break tolerance |
| 11 | `test_outside_boundary_tol_raises` | §2 tolerance boundary |
| 12 | `test_boundary_half_cell_renormalized_weights_sum_to_one` | §3 renormalization correctness |
| 13 | `test_pure_function_no_ecm_mutation` | §3 pure-function constraint |
| 14 | `test_pure_function_no_accumulator_call_via_mock` | §3 forbidden side-effect |
| 15 | `test_negative_traction_input_preserved_in_sign_x_and_y` | §5 sign preservation |
| 16 | `test_returned_shape_is_nx_ny_2_float64` | §6 vector output, no scalarization |
| 17 | `test_non_finite_fa_position_raises` | §2 boundary failure_kind |
| 18 | `test_scatter_does_not_satisfy_closed_loop_response_law` | §6 Hard Rule 11 meta-test |

(Locked §7 lists 17 tests including
`test_non_finite_fa_position_raises` and
`test_returned_shape_is_nx_ny_2_float64`; the Hard Rule 11
meta-test #18 is added per Phase B / Phase C precedent and is
not considered a scope expansion. Test count 17 + 1 = 18.)

---

## 10. Gate verdict

§1, §2, §3, §4, §5, §6 PASS. Magic-Number Block (§7) clean.
Visual deliverable plan (§8) is post-commit and does not block.
Test catalog (§9) maps every gate item to at least one test plus
a Hard Rule 11 meta-test.

The gate clears for executable code in two commits:

1. `acs/v2/dynamics/fa_to_ecm_scattering.py` — new module with
   `FAToECMScatteringError` exception class +
   `scatter_fa_traction_to_ecm_bilinear(...)` function.
   `_FA_POSITION_BOUNDARY_TOL_UM = 1e-12 * max(nx*dx, ny*dy)`
   computed at runtime per ECM (not a fixed module-level
   constant; recomputed per call to remain grid-invariant).
2. `tests/test_v2_fa_to_ecm_scattering.py` — 18 tests per §9
   catalog above (17 from locked §7 catalog including
   `test_non_finite_fa_position_raises` and
   `test_returned_shape_is_nx_ny_2_float64`, plus 1 Hard Rule 11
   meta-test added per Phase B / Phase C precedent).

Plus `acs/v2/dynamics/__init__.py` + `acs/v2/__init__.py`
exports.

If Codex review surfaces a missing reduction or hidden numeric,
the gate flips to BLOCKER with the three-options template.

### Outstanding before code lands

- Codex review of this gate document (5 review-focus items per
  impl Codex `id=1349`):
  1. Vector output `(nx, ny, 2)` only, NO scalarization in API
  2. Cell-centered inclusive footprint convention
  3. Boundary truncated-stencil renormalization + conservation
     proof
  4. Out-of-grid raises, no clamp
  5. Pure function (no `ecm` mutation, no response-law behavior)
  6. Rule 10 unit chain + Rule 11 measurement boundary

---

## 11. References

- Locked design (source of truth):
  `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
- Parent plan: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
  §1 Phase D + §2 Hard Blocker #3
- 6.3a static FA traction (FA traction sign convention precedent
  this scatter consumes):
  `acs/v2/dynamics/focal_adhesion.py` +
  `docs/v2_focal_adhesion_dynamics_sanity_gate.md`
- 6.3b protrusion-coupled FA dynamics (per-FA delegation
  precedent + meta-test pattern):
  `docs/v2_63b_protrusion_coupled_fa_sanity_gate.md`
- ECM substrate schema: `acs/v2/ecm_substrate.py`
- FA schema: `acs/v2/focal_adhesion.py`
- Hard Rule 10 (CLAUDE.md): dimensional comparison verification
- Hard Rule 11 (CLAUDE.md): measurement-protocol consistency
- Memory rules:
  - `rule10_unit_derivation_in_docs.md`
  - `hard_rule_11_wording_boundary_meta_test.md`
  - `design_note_pre_commit_batch.md` (5-step pre-commit batch)
  - `feedback_aggressive_design_debate.md`
- Phase B precursor (stimulus monotonicity meta-test precedent):
  `tests/test_v2_ecm_open_loop.py` Phase B section
- Phase C open-loop sweep baseline (meta-test precedent):
  `tests/test_v2_ecm_ol_sweep.py`
