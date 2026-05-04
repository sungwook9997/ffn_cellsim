# V2 Phase 1 — ECM→FA Bias Sanity Gate (Hard Blocker #4)

**Status**: pre-execution Sanity Gate per Plan §10/§11 + closed-loop
ECM gate phased plan
(`docs/v2_closed_loop_ecm_gate_phased_plan_locked.md` §1 Phase D
entry blocker = Hard Blocker #3 + #4 interface locks). Written from
the locked design in
`docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` (4-round
adversarial design-discussion lock; cross-room dispatch implicit
via Syncthing).

**Contract**: this gate is the document Plan §11 requires **before**
any executable bias / sampler physics lands in
`acs/v2/dynamics/ecm_to_fa_bias.py`. If any item §1–§6 below
fails its check, Hard Blocker #4 reports `status=blocker` to PI
with at least three concrete options. No partial physics commits
under uncertainty.

**Source of truth**:
- locked design: `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
- sister lock: `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
- parent plan: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
  §1 Phase D + §2 Hard Blocker #4 + §3 effective_stiffness guard
- design-discussion ledger: MCP `id=1365-1370`
- impl-work review path: MCP `id=1372` (cross-room handoff via
  Syncthing)
- Hard Blocker #3 Sanity Gate (sister gate, geometry convention
  source): `docs/v2_fa_to_ecm_scattering_sanity_gate.md`
- 6.3a static FA traction (sign-convention precedent):
  `acs/v2/dynamics/focal_adhesion.py`
- 6.3b protrusion-coupled FA dynamics (typed-key multiplier
  precedent): `acs/v2/dynamics/protrusion_coupled_focal_adhesion.py`
- ECM substrate schema: `acs/v2/ecm_substrate.py`
- FA schema: `acs/v2/focal_adhesion.py`

**Sequencing context**: Hard Blocker #4 is the second (and last) of
Phase D entry interface blockers. With Hard Blocker #3 locked +
Sanity-Gated + implemented (`3ba37fa` lock + `b90b07f` Sanity Gate
+ `883efc8` code), this gate clears the remaining Phase D entry
prerequisite. Phase D no-op scaffolding (the
`step_fa_to_ecm_response` integrator combining scatter + bias) is
the next unit after this gate clears + the corresponding code
commit lands.

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

- A pure read function `sample_ecm_at_fa_positions(adhesions, ecm)
  -> ECMSampledAtFAs` returning per-FA bilinear-interpolated ECM
  field values, using the **same geometry convention as Hard
  Blocker #3 scatter** (cell-centered, inclusive footprint,
  boundary half-cell truncated stencil with renormalized
  weights, out-of-grid raise).
- A pure neutral function `compute_ecm_to_fa_bias_neutral(adhesions,
  ecm) -> ECMToFABiasResult` returning per-FA per-rate multipliers
  all equal to 1.0, plus sampled diagnostics.
- Module constant `RATE_NAMES = ("k_maturity_per_s",
  "k_bind_per_s", "k_unbind_per_s")` fixing column ordering of
  `multipliers_per_fa[:, k]`.
- `FAToECMBiasError(ValueError)` with `failure_kind` literal
  `"fa_bias_position_outside_ecm_grid"` /
  `"non_finite_fa_position"`.
- Frozen dataclasses `ECMSampledAtFAs` + `ECMToFABiasResult` per
  locked §1.
- 19-test catalog: 18 from locked §6 (NOT §7 — §7 is the Phase D integrator pointer; test catalog lives at §6 of the lock) + 1 new `test_invalid_fa_schema_propagates_validate_error` regression for the per-FA `fa.validate()` boundary required by Codex review id=1378 (mirrors the Hard Blocker #3 id=1359 fix).

### Explicitly out of scope (will FAIL gate if introduced)

- Active non-1.0 multiplier mapping. Phase D returns exactly 1.0
  for every FA, every rate. Active mapping (banded / continuous
  per-rate function / orientation-anisotropic) is Phase E
  (separate function `compute_ecm_to_fa_bias_active`, separate
  gate, separate lock).
- `traction_scale_nN` mutation. Bias only modifies FA *rates*;
  it does not touch FA traction magnitude. Per the locked phased
  plan §3 effective_stiffness guard.
- FA nucleation. The bias does not spawn or remove FAs. FA
  nucleation requires RNG + reproducibility contract per locked
  phased plan §4 (deferred to a separate Sanity Gate).
- RNG / stochastic sampler. The bias is fully deterministic.
- State-label transitions. The bias does not change FA `state`
  literal (preserves 6.3a's no-auto-state-transition contract).
- FA / ECM state mutation. Both inputs are unchanged on return
  (pure functions).
- Effective-stiffness helper using geometry fields. Per the
  locked phased plan §3 guard, any helper consuming
  `fiber_density` / `orientation_tensor` is ECM→FA bias by
  another name. Phase D's sampled diagnostics expose those
  fields for **audit/test only**, never for bias decision.
- Single-function-with-optional-mapping_law param. The Phase D /
  Phase E function naming is locked separate
  (`_neutral` vs `_active`) to prevent silent activation if a
  caller passes an active mapping into what should be a neutral
  call (locked design §1).
- Vector accumulated traction sample. The ECM schema's
  `accumulated_traction_nNs_per_um2` is scalar `(nx, ny)`; the
  sampler returns scalar `(N_FA,)`. Vector accumulated history
  would require a schema change (deferred).
- Empty-FA-list raise. Empty `adhesions` returns an empty
  result with shape-consistent zero-row arrays, mirroring Hard
  Blocker #3 / 6.3a precedents.
- Caller-supplied rate-name reordering. The result's
  `rate_names` field equals `RATE_NAMES` constant; caller does
  not supply or reorder.
- Arrays in `diagnostics_dict`. Plain serializable values only
  (`int`, `float`, `str`); arrays live in typed dataclass fields.

---

## 1. Dimensional analysis

Inputs and outputs (every variable's unit explicit):

| Symbol | Meaning | Unit | Source |
|---|---|---|---|
| `position_um_xy` | per-FA xy position | μm | `FocalAdhesionState` |
| `origin_um_xy` | ECM grid bottom-left corner | μm | `ECMSubstrateState` |
| `spacing_um` | uniform cell edge | μm | `ECMSubstrateState` |
| `grid_shape` | `(nx, ny)` cell counts | dimensionless integers | `ECMSubstrateState` |
| `stiffness_kpa[i, j]` | per-cell substrate stiffness | kPa | `ECMSubstrateState` |
| `fiber_density[i, j]` | per-cell collagen fiber density | dimensionless `[0, 1]` | `ECMSubstrateState` |
| `ligand_density[i, j]` | per-cell ligand density | dimensionless `[0, 1]` | `ECMSubstrateState` |
| `orientation_tensor[i, j, :, :]` | per-cell traceless symmetric tensor | dimensionless `[-1, +1]` component-wise | `ECMSubstrateState` |
| `accumulated_traction_nNs_per_um2[i, j]` | per-cell cumulative traction history | nN·s/μm² | `ECMSubstrateState` |
| `multipliers_per_fa[i, k]` | per-FA per-rate multiplier | dimensionless ≥ 0 | computed (Phase D = 1.0) |
| `bilinear weight w_kl` | interpolation weight | dimensionless `[0, 1]` | derived |

### Reductions

- **Sampler bilinear interpolation**: weighted sum over 4 (or
  fewer at boundary) cells, weights dimensionless and summing to
  1.0 → unit-preserving interpolation. Each sampled scalar field
  preserves its underlying unit; the 2-tensor field
  preserves its component-wise `[-1, +1]` bound at the
  interpolated point.
- **Multiplier (Phase D)**: dimensionless `1.0` constant; no
  unit chain at this layer. Phase E's `_active` function would
  introduce a chain when it maps sampled values to non-1.0
  multipliers — that's that gate's responsibility.
- **No Rule 10 cross-unit comparison** at this layer. The
  sampler is read-only; the bias computation in Phase D is
  trivially `multiplier = 1.0` with no unit reduction.

### Status

PASS — every quantity preserves its source unit through bilinear
interpolation. Hard Rule 10 trivially satisfied at Phase D (no
cross-unit comparison). The locked phased plan §6 closed-loop
kPa↔nN/μm² unit-chain proof is the **Phase E** gate's
responsibility and explicitly NOT this gate's.

---

## 2. Boundary cases

| Case | Guard | Failure-kind |
|---|---|---|
| Empty `adhesions` tuple | early-return `ECMToFABiasResult(fa_ids=(), multipliers_per_fa=np.zeros((0, 3)), rate_names=RATE_NAMES, sampled_diagnostics=ECMSampledAtFAs(fa_ids=(), …all zero-row arrays…), diagnostics_dict={"n_adhesions": 0, "max_multiplier": 1.0, "min_multiplier": 1.0, "sampler_geometry": "bilinear_cell_centered"})`. **No raise** (locked §0 forbids `empty_adhesions` failure_kind). | (no exception) |
| `position_um_xy` non-finite (NaN, ±inf) | reject before bilinear computation | `non_finite_fa_position` |
| `FocalAdhesionState` schema invariant violation (e.g., `state == "unbound"` with non-zero `traction_force_nN_xy`, `bound_fraction == 0.0` with non-zero traction, etc.) | per-FA `fa.validate()` runs **after** `non_finite_fa_position` check and **before** bilinear sampling, so the lock-specific `non_finite_fa_position` failure_kind still fires first; schema-only invariants then propagate as the schema's `ValueError`. Mirrors Hard Blocker #3 fix `883efc8` and the local dynamics precedent in `step_focal_adhesions_static`. | inherited `FocalAdhesionState.validate()` `ValueError` |
| FA position strictly outside `[origin_x, origin_x + nx*dx] × [origin_y, origin_y + ny*dy]` plus `_FA_POSITION_BOUNDARY_TOL_UM` | raise BEFORE sampler bilinear; no silent clamp | `fa_bias_position_outside_ecm_grid` (distinct from Hard Blocker #3 `fa_position_outside_ecm_grid` because the failing op is read/bias, not scatter; physical footprint rule identical) |
| FA position exactly at footprint corner | inclusive boundary; truncated stencil reads from the single in-grid corner cell | (no exception) |
| FA position exactly at footprint edge midpoint | inclusive; truncated stencil reads from edge cells with renormalized weights | (no exception) |
| FA position strictly at a cell center | bilinear degenerates: weight 1.0 to that center, 0 elsewhere | (no exception) |
| FA position exactly at 4-center crossing | weights `(0.25, 0.25, 0.25, 0.25)` | (no exception) |
| `ecm.spacing_um` ≤ 0 or non-finite | inherited from `ECMSubstrateState.validate()` at function entry (per the local dynamics precedent in `acs/v2/dynamics/ecm_open_loop.py` and the Hard Blocker #3 Sanity Gate §2 fix `b90b07f`) | inherited `ECMSubstrateState.validate()` `ValueError` |
| `ecm.grid_shape` non-positive integer | same — `ecm.validate()` at entry catches | inherited (same) |
| FA tuple containing duplicate `adhesion_id` | not the bias's responsibility (caller-side concern); function just iterates | (no exception) |
| Caller passes a list instead of a tuple | accepted as runtime convenience, matching the Hard Blocker #3 scatter precedent. The strict contract per locked §1 is the tuple form. | (no exception; documented in docstring) |

### Status

PASS — every boundary case has either an explicit failure_kind
or a documented "absent / handled by truncated stencil / empty
returns empty / accepted convenience" rule. No silent clamp; no
`empty_adhesions` failure_kind (per locked §0 ban).

---

## 3. Conservation invariants

### Pure-function constraint (locked §1)

Both functions are pure read / no-mutation:
- `ecm` byte-equality before and after the call (verified by
  per-field equality test on `stiffness_kpa`, `ligand_density`,
  `fiber_density`, `orientation_tensor`,
  `accumulated_traction_nNs_per_um2`, plus `origin_um_xy`,
  `spacing_um`, `grid_shape`).
- Input `adhesions` tuple is not mutated.
- No `accumulate_prescribed_traction` call (verified by mocking
  the function and asserting `mock.assert_not_called()` per
  Hard Blocker #3 precedent — same pure-read-no-side-effect
  pattern).
- No response field update on `ecm`.
- No FA `state` literal change anywhere.

### Sampler bilinear weight conservation

The 4-cell interior bilinear stencil's weights sum to 1.0 by
construction; the boundary half-cell renormalized stencil's
retained weights sum to 1.0 after renormalization (matching Hard
Blocker #3 conservation logic). Sampled scalar field at the FA
position equals the bilinear-interpolated value — interior:
exact float64; boundary half-cell: float64 round-off-equivalent
to the renormalized weighted sum.

### Multiplier `output structural` invariant (Phase D)

`multipliers_per_fa.shape == (N_FA, 3)` and every entry equals
`1.0`. No accumulation, no across-FA / across-rate coupling.
This is a **structural** invariant (Phase D contract), not a
physical conservation; it's tested via
`test_neutral_multipliers_all_ones_for_arbitrary_fa_list`.

### Determinism

- No RNG. `random` and `np.random` not imported.
- `dict.get`-style lookups with fixed defaults (none required at
  Phase D since multipliers are constant 1.0).
- Identical inputs produce bit-identical outputs.

### Status

PASS — pure-function constraints are enforceable by structural
tests. Sampler bilinear weight conservation matches Hard Blocker
#3 sister gate. Phase D structural multiplier invariant is the
default-no-op contract.

---

## 4. Numerical sanity

### Float precision

float64 throughout. Per-FA work is `O(8)` for interior bilinear
sampler (one read per cell × 4 cells × per-field count, plus the
constant-1.0 multiplier write). Memory is `O(N_FA · K_fields +
N_FA · 3)` for the result arrays, where `K_fields` is the
sampler's field count (= 5 currently: stiffness, fiber, ligand,
orientation 2×2 = 4 → 5 distinct numpy arrays).

### `_FA_POSITION_BOUNDARY_TOL_UM`

Reused from Hard Blocker #3 module if exposed, or recomputed at
function entry as `1e-12 * max(nx*dx, ny*dx)`. Per the Hard
Blocker #3 Sanity Gate §1, this is documented as a relative
float64 numerical tie-break (derivable from float64 epsilon +
safety factor; grid-invariant in relative terms; not chosen to
fit any test target).

### No auto-shrink, no clipping

- `_FA_POSITION_BOUNDARY_TOL_UM` is an inclusion-check
  tolerance only; it does not shift the FA position or modify
  the read.
- Bilinear weights computed from validated FA position, not
  from a snapped grid index.
- Out-of-grid raises; no silent move-to-nearest-cell.

### Stability

Both functions are single-step (no time integration); no CFL
bound applies. The sampler's bilinear interpolation is stable
for any positive `spacing_um`. No dt-rate gate.

### Status

PASS — single-step pure read; the only new numeric is the
inherited `_FA_POSITION_BOUNDARY_TOL_UM` from Hard Blocker #3
(documented + reused, not a new tunable).

---

## 5. Sign / sense check

| Quantity | Direction | One-line check |
|---|---|---|
| Bilinear weight `w_kl` | non-negative by construction; `w_kl ∈ [0, 1]` | products of `(1-fx, fx, 1-fy, fy)` factors. |
| Renormalized weight `w'_kl` (boundary half-cell) | non-negative; `w'_kl = w_kl / s` for `s > 0`; `Σ w'_kl = 1`. | identical to Hard Blocker #3 logic. |
| Sampled `stiffness_kpa[FA]` | non-negative (per ECM schema validate); inherited from underlying field. | one-line: weighted sum of non-negative values is non-negative. |
| Sampled `fiber_density[FA]`, `ligand_density[FA]` | in `[0, 1]` by linear-combination convexity (weights non-negative + sum 1). | tested via `test_sampler_at_cell_center_returns_field_value_exactly`. |
| Sampled `orientation_tensor[FA]` | component-wise in `[-_ORIENTATION_BOUND, +_ORIENTATION_BOUND]` by linear-combination on the same schema bound. | tested by sampling at a known-tensor cell. |
| Sampled `accumulated_traction_nNs_per_um2[FA]` | non-negative; inherited from underlying field's monotone-non-decreasing schema invariant. | inherited from ECM-OL-1 contract. |
| Multiplier `m[fa, k]` (Phase D) | exactly `1.0` ≥ 0. | `test_neutral_multipliers_all_ones_for_arbitrary_fa_list`. |

### Status

PASS — every sign is preserved by the bilinear linear
combination's convexity (non-negative weights summing to 1) on
top of the ECM schema's existing field bounds. Multipliers are
trivially ≥ 0 by being exactly `1.0`.

---

## 6. Measurement-protocol consistency (Hard Rule 11)

Phase D output modality is:
- `multipliers_per_fa: (N_FA, 3)` of all `1.0` values
- `sampled_diagnostics: ECMSampledAtFAs` for audit / test only
- `diagnostics_dict` of plain serializable values

The locked design §0 forbids active non-1.0 multipliers + active
mapping. The Phase D / Phase E function naming separation
(`_neutral` vs `_active`) is the structural Hard Rule 11
boundary: a future Phase E commit cannot silently activate
through this Phase D function.

### No scalarization at sampler time

The sampler returns each ECM field at its **native shape per FA**
(scalar field → `(N_FA,)`, tensor field → `(N_FA, 2, 2)`). The
sampler does NOT reduce the orientation tensor to a scalar
(magnitude / determinant / first principal value). Phase E
response laws decide the reduction. This mirrors Hard Blocker
#3's "no scalarization at scatter time" boundary.

### Geometry consistency with Hard Blocker #3

Per locked §5 Sanity Gate item 6 + the
`test_sampler_consistency_with_scatter_geometry` test in §9: the
same FA position scattered by Hard Blocker #3 and sampled by Hard
Blocker #4 must use the same cell-centered convention,
inclusive-footprint policy, boundary half-cell renormalization,
and out-of-grid raise behavior. A mismatch creates a closed-loop
geometric hysteresis the locked phased plan §1 Item 5 sensitivity
sweep would surface.

### Runtime meta-test (Hard Rule 11 boundary)

Per Phase B / Phase C / Hard Blocker #3 precedent, a runtime
meta-test
`test_no_active_law_invoked_in_phase_d` will encode the
boundary at runtime. Tautological assertion; documented in the
docstring + section header as the wording-boundary protection.

### Status

PASS — measurement modality is the locked output type; no
scalarization at sampler time; geometry consistency with Hard
Blocker #3 is testable; runtime meta-test enforces the Phase D /
Phase E boundary.

---

## 7. Magic-Number Block check

Every numeric in the Hard Blocker #4 implementation is one of:

| Symbol | Source | Magic-Number Block status |
|---|---|---|
| `position_um_xy`, `origin_um_xy`, `spacing_um`, `grid_shape`, all sampled field values | runtime inputs (FA / ECM schemas) | **not a magic number** — caller-supplied via schema |
| `_FA_POSITION_BOUNDARY_TOL_UM = 1e-12 * max(nx*dx, ny*dx)` | inherited from Hard Blocker #3 (or recomputed locally) | named tie-break; derivable from float64 epsilon + safety factor; grid-invariant in relative terms; not chosen to fit any test target. Already cleared for Hard Blocker #3 Sanity Gate. |
| `RATE_NAMES = ("k_maturity_per_s", "k_bind_per_s", "k_unbind_per_s")` | locked design §1 module constant | tuple of literal strings; not a numeric tunable. |
| `multipliers_per_fa` = `1.0` (Phase D no-op default) | locked design §1 Phase D default | the `1.0` is a structural identity / no-op marker; trivially derivable, grid-invariant, not chosen to fit any test target. |
| `diagnostics_dict["max_multiplier"] = 1.0`, `min_multiplier = 1.0` | derived from the multiplier array (Phase D) | trivial — not a numeric tunable. |

Magic-Number Block test pass:

1. **Derivable**: every value is either runtime caller input,
   inherited from Hard Blocker #3 (`_FA_POSITION_BOUNDARY_TOL_UM`),
   or a structural identity (`1.0` neutral multiplier,
   `RATE_NAMES` literal tuple).
2. **Grid-invariant**: `_FA_POSITION_BOUNDARY_TOL_UM` scales with
   the grid extent (relative tolerance); the structural `1.0`
   multiplier is dimensionless and grid-invariant by definition.
3. **Not fitting**: no value was selected to make any specific
   test pass.

### Status

PASS — Magic-Number Block clean. The `_FA_POSITION_BOUNDARY_TOL_UM`
constant is inherited (not new); the `1.0` neutral multiplier is
a structural identity, not a tunable.

---

## 8. Visual deliverable plan (post-commit)

Aligns with the precedent set by `scripts/run_p1_alpha_gate.py`,
`scripts/run_ecm_ol_harness.py`,
`scripts/run_protrusion_coupled_fa_harness.py`,
`scripts/run_ecm_ol_sensitivity_sweep.py`, and the Hard Blocker
#3 visual deliverable plan. **Out of scope for the Sanity Gate
commit + the code/tests commit** — the visual deliverable for
Hard Blocker #4 is reserved for the post-Phase-D-no-op-scaffolding
unit, where it makes more sense to visualize bias + scatter
together (Phase D no-op integrator's first run).

Planned (NOT committed by this gate):
- Combined Hard Blocker #3 + #4 smoke runner (when Phase D no-op
  scaffolding lands): single FA list × single ECM grid →
  scatter to `(nx, ny, 2)` field + sample per-FA fields →
  visualize both.
- Per-scenario `summary.html` + `metadata.json` with the
  scatter + sampler diagnostic side-by-side.

---

## 9. Test catalog (19, per locked §6 + fa.validate() regression)

Each test exercises one Sanity Gate item or one forbidden
behavior. Locked §6 (NOT §7 — §7 is the Phase D integrator
pointer) lists 18; this catalog adds the
`test_invalid_fa_schema_propagates_validate_error` regression
required by the fa.validate() boundary fix per Codex review
id=1378 (mirror of Hard Blocker #3 id=1359 fix).

| # | Test name | Gate item |
|---|---|---|
| 1 | `test_neutral_multipliers_all_ones_for_arbitrary_fa_list` | §3 structural invariant + §5 multiplier sign |
| 2 | `test_neutral_multipliers_match_RATE_NAMES_column_order` | §1 module constant + structural |
| 3 | `test_sampler_at_cell_center_returns_field_value_exactly` | §1 unit preservation + §3 conservation interior |
| 4 | `test_sampler_at_grid_crossing_bilinear_average` | §3 boundary case at 4-center crossing |
| 5 | `test_sampler_boundary_half_cell_renormalized` | §3 renormalization correctness |
| 6 | `test_sampler_outside_footprint_raises_fa_bias_position_outside_ecm_grid` | §2 / §4 out-of-grid policy |
| 7 | `test_sampler_input_order_preserved` | §6 input order preservation |
| 8 | `test_sampler_accumulated_traction_is_scalar_shape_n_fa` | §1 / §0 forbidden vector accumulated |
| 9 | `test_sampler_orientation_tensor_shape_n_fa_2_2` | §6 native-shape / no-scalarization protection (orientation tensor returned full 2×2 per FA) |
| 10 | `test_neutral_bias_no_ecm_mutation` | §3 pure-function constraint |
| 11 | `test_neutral_bias_no_fa_mutation` | §3 pure-function constraint |
| 12 | `test_diagnostics_dict_only_serializable_plain_values` | §0 forbidden arrays in diagnostics |
| 13 | `test_diagnostics_dict_keys_are_locked_set` | §6 measurement-protocol; locked keys (n_adhesions, max/min_multiplier, sampler_geometry) |
| 14 | `test_empty_fa_list_returns_empty_result_not_raise` | §0 forbidden empty_adhesions failure |
| 15 | `test_sampler_consistency_with_scatter_geometry` | §6 measurement-protocol (sister-gate consistency with Hard Blocker #3) |
| 16 | `test_returned_shapes_match_locked_signature` | §1 dataclass shape contract |
| 17 | `test_no_active_law_invoked_in_phase_d` | §0 forbidden + Hard Rule 11 meta-test |
| 18 | `test_non_finite_fa_position_raises` | §2 boundary failure_kind |
| 19 | `test_invalid_fa_schema_propagates_validate_error` | §2 fa.validate() boundary at function entry; mirrors Hard Blocker #3 fix `883efc8`; expected schema invariant: `state="unbound"` with non-zero `traction_force_nN_xy` raises schema `ValueError` BEFORE any sampler / multiplier output |

All 18 tests #1–#18 come from locked §6 catalog (no scope
expansion vs lock); #19 is the fa.validate() regression that
mirrors the Hard Blocker #3 `883efc8` fix and was added per
Codex review id=1378. The Phase D / Phase E function naming
separation is verified structurally by the absence of
`compute_ecm_to_fa_bias_active` in the module — caller cannot
trigger it through this gate's public surface.

---

## 10. Gate verdict

§1, §2, §3, §4, §5, §6 PASS. Magic-Number Block (§7) clean.
Visual deliverable plan (§8) is post-Phase-D-scaffolding and does
not block. Test catalog (§9) maps every gate item to at least
one test plus the Hard Rule 11 meta-test.

The gate clears for executable code in two commits:

1. `acs/v2/dynamics/ecm_to_fa_bias.py` — new module with
   `RATE_NAMES` module constant, `FAToECMBiasError` exception,
   `ECMSampledAtFAs` + `ECMToFABiasResult` frozen dataclasses,
   `sample_ecm_at_fa_positions` + `compute_ecm_to_fa_bias_neutral`
   functions. `_FA_POSITION_BOUNDARY_TOL_UM` either imported from
   Hard Blocker #3 module (preferred) or recomputed locally.
2. `tests/test_v2_ecm_to_fa_bias.py` — 19 tests per §9 catalog
   (18 from locked §6 + 1 fa.validate() regression per Codex
   review id=1378).

Plus `acs/v2/dynamics/__init__.py` + `acs/v2/__init__.py`
exports.

If Codex review surfaces a missing reduction or hidden numeric,
the gate flips to BLOCKER with the three-options template.

### Outstanding before code lands

- Codex review of this gate document (5 review-focus items per
  design-discussion Codex `id=1370` + impl Codex backflow ack):
  1. RATE_NAMES module constant + result column order (§1)
  2. Empty FA list returns empty result (no raise) (§2)
  3. diagnostics_dict plain serializable values only (§0 / §6)
  4. accumulated traction sampled as scalar (N_FA,) (§1 / §0
     forbidden vector)
  5. Pure functions only — no state mutation, no active non-1.0
     law (§0 / §3 / §6 meta-test)

---

## 11. References

- Locked design (source of truth):
  `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
- Sister lock (Hard Blocker #3 — geometry convention):
  `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
- Sister Sanity Gate (Hard Blocker #3 — gate-item template):
  `docs/v2_fa_to_ecm_scattering_sanity_gate.md`
- Parent plan: `docs/v2_closed_loop_ecm_gate_phased_plan_locked.md`
  §1 Phase D + §2 Hard Blocker #4 + §3 effective_stiffness guard
- 6.3a static FA traction (sign-convention precedent):
  `acs/v2/dynamics/focal_adhesion.py` +
  `docs/v2_focal_adhesion_dynamics_sanity_gate.md`
- 6.3b protrusion-coupled FA dynamics (typed-key multiplier
  precedent):
  `docs/v2_63b_protrusion_coupled_fa_sanity_gate.md`
- ECM substrate schema: `acs/v2/ecm_substrate.py`
- FA schema: `acs/v2/focal_adhesion.py`
- Hard Rule 10 (CLAUDE.md): dimensional comparison verification
- Hard Rule 11 (CLAUDE.md): measurement-protocol consistency
- Memory rules:
  - `rule10_unit_derivation_in_docs.md`
  - `hard_rule_11_wording_boundary_meta_test.md`
  - `design_note_pre_commit_batch.md` (5-step pre-commit batch
    applied to this gate)
  - `feedback_aggressive_design_debate.md`
- Phase B precursor (stimulus monotonicity meta-test precedent):
  `tests/test_v2_ecm_open_loop.py` Phase B section
- Phase C open-loop sweep baseline (meta-test precedent):
  `tests/test_v2_ecm_ol_sweep.py`
- Hard Blocker #3 implementation (pure-function pattern + mock
  check precedent for `test_neutral_bias_no_fa_mutation` + sister
  geometry test): `acs/v2/dynamics/fa_to_ecm_scattering.py` +
  `tests/test_v2_fa_to_ecm_scattering.py`
