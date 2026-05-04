# V2 Phase 1 — Focal Adhesion Dynamics 6.3a Sanity Gate

**Status**: pre-execution Sanity Gate per Plan §10/§11 + forward
roadmap "Separated dynamics modes". Written from the
design-discussion lock in MCP id=1070→1077 (Claude+Codex 5-round
adversarial lock, two implementation guards locked in id=1074).

**Contract**: this gate is the document Plan §11 requires **before**
any executable focal-adhesion physics lands in
`acs/v2/dynamics/focal_adhesion.py`. If any item §1–§6 below fails
its check, 6.3a reports `status=blocker` to PI with at least three
concrete options. No partial physics commits under uncertainty.

**Source of truth**: design-discussion `topic=v2-layer-2-separated-dynamics`,
MCP id=1070, 1072, 1073, 1074, 1076, 1077.

**Sequencing context**: this is the active design-locked next unit
after the P1 alpha completion (`f28c109`). The ECM-OL-1
(`29a360f`) and ECM-OL-2 (`2c76c0c`) sidecar commits landed in
parallel from a separate impl-work decision branch and are
considered "sidecar landed units" per the reconciliation in
design-discussion id=1086/1088. The active design-locked unit is
this one (6.3a).

---

## 0. Scope (locked, no scope creep allowed)

### In scope (this unit only)

- Static-boundary deterministic FA traction algebra.
- Explicit-rate deterministic maturation / binding only when
  caller-supplied rates are present.
- Per-FA cell-on-substrate traction vector computation.
- Per-FA substrate-on-cell reaction vector (= negation, action-reaction).
- Radial/tangential decomposition relative to a supplied frozen
  centroid.
- Optional validation-only `max_traction_nN` upper bound (raises;
  never clamps).

### Explicitly out of scope (will FAIL gate if introduced)

- Any active-contour `step()` call. The cell contour/centroid is a
  frozen geometric reference for this unit.
- ECM substrate update. ECM-OL-1's `accumulate_prescribed_traction`
  may be used downstream by a different unit; this unit produces
  the per-FA cell-on-substrate traction only.
- Stochastic nucleation hazards / lifetime samplers.
- Protrusion-coupled FA formation.
- ECM-biased FA formation.
- Automatic state-label transitions
  (e.g. nascent→mature without an explicit caller-supplied rule).
- Magic clipping. Any default-without-supplied-value cap is a Magic
  Number Block test-3 fail.

---

## 1. Dimensional analysis

Inputs and outputs (every variable's unit explicit):

| Symbol | Meaning | Unit | Source |
|---|---|---|---|
| `position_um_xy` | FA center on substrate plane | μm | `FocalAdhesionState` (Cycle C schema) |
| `centroid_um_xy` | frozen cell centroid | μm | caller |
| `maturity` | adhesion maturity | dimensionless [0, 1] | `FocalAdhesionState` |
| `bound_fraction` | clutch bound fraction | dimensionless [0, 1] | `FocalAdhesionState` |
| `traction_scale_nN` | per-FA traction magnitude scale | nN | `FocalAdhesionDynamicsParameters` (caller-supplied, no default) |
| `dt_fa_s` | integration step | s | `FocalAdhesionDynamicsParameters` (caller-supplied, no default) |
| `k_maturity_per_s` | maturity-relaxation rate | 1/s | optional, no default |
| `k_bind_per_s` | bound_fraction toward 1 | 1/s | optional, no default |
| `k_unbind_per_s` | bound_fraction toward 0 | 1/s | optional, no default |
| `max_traction_nN` | per-FA validation upper bound | nN | optional, no default; raise on exceed (no clamp) |
| `traction_axis_xy` (caller-supplied) | direction override | dimensionless 2-vector | optional |
| `traction_force_nN_xy` | per-FA cell-on-substrate force | nN | computed |
| `substrate_reaction_nN_xy` | per-FA substrate-on-cell force | nN | `-traction_force_nN_xy` |

### Reductions

- **Maturity update**: `maturity_new = maturity + dt_fa_s · k_maturity_per_s · (1 - maturity)` reduces as `[s] · [1/s] · [dimensionless] = [dimensionless]` ✓ on the increment.
- **Bound_fraction update** (optional, mirror form): identical
  dimensional pattern ✓.
- **Traction magnitude (default radial)**:
  `traction_force_nN_xy = traction_scale_nN · maturity · bound_fraction · unit(centroid - position)`
  reduces as `[nN] · [dimensionless] · [dimensionless] · [dimensionless] = [nN]` ✓.
- **Substrate reaction**: `[nN]` from negation. ✓
- **Radial/tangential decomposition**: dimensionless unit vectors
  multiply per-FA force `[nN]` to give per-FA scalar components in `[nN]`. ✓

### Status

PASS — no kPa conversion is performed in this unit (FA traction is
already in `nN`); the closed-loop ECM gate's Rule 10 unit-chain
proof for kPa↔nN/μm² lives downstream and is *not* this gate's
responsibility. Per-FA force is the only stress-like quantity.

---

## 2. Boundary cases

| Case | Guard | Failure-kind |
|---|---|---|
| `state == "unbound"` | force traction to (0, 0) and skip rate updates that would change `bound_fraction` away from 0 | (no exception; logical guard) |
| `state == "released"` | force traction to (0, 0); permit no further bound_fraction increase | (no exception; logical guard) |
| `bound_fraction == 0.0` | traction must be (0, 0) regardless of `state` (no traction without bound) | per-FA assertion in update |
| Non-finite `position_um_xy` / `centroid_um_xy` | reject before compute | `non_finite_position` / `non_finite_centroid` |
| Empty / non-string `adhesion_id` or `cell_id` | already enforced by schema | (inherited) |
| `centroid_um_xy` length != 2 | reject | `centroid_shape` |
| `traction_axis_xy` length != 2 (when supplied) | reject | `axis_shape` |
| `traction_axis_xy` zero magnitude (when supplied) | reject | `axis_zero_magnitude` |
| Position equals centroid (any FA) | radial decomposition (§6) requires a defined radial unit; reject regardless of `traction_axis_xy` | `radial_axis_undefined` |
| `dt_fa_s < 0` or non-finite | reject | `dt_invalid` |
| `dt_fa_s` is bool or non-int/float | reject | `dt_invalid` |
| `traction_scale_nN < 0` or non-finite | reject | `traction_scale_invalid` |
| `traction_scale_nN` is bool | reject | `traction_scale_invalid` |
| Any supplied rate (`k_maturity_per_s`, `k_bind_per_s`, `k_unbind_per_s`) negative, non-finite, bool, or non-numeric | reject **before** `dt·max(rate)` evaluated, so a negative rate cannot flip §5 sign convention | `rate_invalid` |
| `max_traction_nN < 0`, non-finite, bool, or non-numeric (when supplied) | reject; `max_traction_nN == 0` is allowed and is interpreted as "only zero-magnitude traction passes" | `max_traction_invalid` |
| `max_traction_nN` supplied with computed traction `||F|| > max_traction_nN` | raise (do **not** clamp) | `max_traction_exceeded` |
| `dt_fa_s · max(rate) > 0.5` (after rate validation) | raise (no auto-shrink) | `dt_rate_violation` |

### Status

PASS — every boundary case has either a logical guard (no
exception, just sets traction to zero) or an explicit
failure-kind exception. No silent clamping or auto-shrink.

---

## 3. Conservation invariants

### Per-FA action-reaction (exact)

For every FA `i`:
`traction_force_nN_xy[i] + substrate_reaction_nN_xy[i] = (0, 0)` exactly (float64 round-off zero, since the substrate reaction is computed as the negation of the traction).

### Aggregate Newton 3

`sum_i traction_force_nN_xy[i] + sum_i substrate_reaction_nN_xy[i] = (0, 0)`. Trivially true given per-FA cancellation, but enforced as a separate test so a future change that, e.g., applies traction
asymmetrically would surface immediately.

### Provenance / immutability

- The input adhesions tuple is never mutated; the result returns a
  new tuple of `FocalAdhesionState` instances (frozen dataclasses).
- The supplied `centroid_um_xy` and `params` are not mutated.

### Status

PASS — both per-FA and aggregate conservation invariants are
testable by exact equality (float64).

---

## 4. Numerical sanity

### Stability bound for explicit-rate updates

The maturity / bound_fraction explicit-Euler update of the form
`x_new = x + dt · k · (1 - x)` (or analogous toward zero) is
linearly stable for `dt · k < 2` and energy-monotone (toward
equilibrium) for `dt · k ≤ 1`. The runtime contract
**`dt_fa_s · max_rate ≤ 0.5`** is a named pre-run safety margin
(factor of 2 below the energy-monotone bound), reusing the
identical rationale documented in
`docs/v2_p1_active_contour_sanity_gate.md` §2 for the active
contour. The 0.5 value is a named module-level constant
(`_DT_RATE_SAFETY_MARGIN = 0.5`), not chosen to make any test
pass.

If `max_rate == 0` (i.e. all rate updates absent), the runtime
contract is vacuously satisfied (`dt · 0 = 0 ≤ 0.5`).

### Float precision

float64 throughout. The traction direction is computed via
`(centroid - position) / ||centroid - position||` with explicit
finite check on the magnitude before division.

### Status

PASS — stability bound is documented, named, and reused from the
active-contour Sanity Gate's already-reviewed safety margin
rationale.

---

## 5. Sign / sense check

| Term | Direction | One-line check |
|---|---|---|
| Cell-on-substrate traction (default radial) | inward toward `centroid_um_xy` from each FA position. Sign: `traction_force_nN_xy[i] = +scale · maturity · bound_fraction · (centroid - position)/||centroid - position||`. | a single FA at (1, 0) with centroid (0, 0), `scale=1`, `maturity=1`, `bound=1` → `traction = (-1, 0)` (inward). |
| Substrate-on-cell reaction | outward from each FA position (negation of traction). | for the same FA → `reaction = (+1, 0)` (outward). |
| Maturity update under positive `k_maturity_per_s` | monotonic increase toward 1. `maturity_new ≥ maturity` for `0 ≤ maturity ≤ 1` and `dt · k ≤ 1`. | with `maturity = 0.4`, `dt · k = 0.5` → `maturity_new = 0.4 + 0.5 · 0.6 = 0.7 > 0.4`. |
| Bound-fraction update | toward 1 under `k_bind`, toward 0 under `k_unbind`. The state table below decides which rate (if any) applies per FA per step. Negative rates rejected at validation (§2 `rate_invalid`). | tested per state. |
| `state == "slipping"` | `bound_fraction` cannot increase (only decrease via `k_unbind` toward 0). | guarded. |

### State-specific rate policy (locked)

The per-FA per-step rate update follows this table. Rates that
"can apply" are still optional — they only update the field if the
caller actually supplied that rate parameter; absent rate means no
change. Rates that "do not apply" are silently ignored even when
supplied, so the caller can pass a single shared `params` for a
heterogeneous adhesion list without runtime errors.

| State | `k_maturity_per_s` | `k_bind_per_s` | `k_unbind_per_s` | Traction non-zero? |
|---|---|---|---|---|
| `unbound` | not applied (maturity stays) | not applied (no rate-driven binding without an explicit nucleation event, which is out of scope) | not applied | no — `bound_fraction == 0` is enforced and traction forced to (0, 0) |
| `nascent` | applied if supplied (maturity ↑ toward 1) | applied if supplied (bound_fraction ↑ toward 1) | not applied | yes if `bound_fraction > 0`; magnitude follows §5 formula |
| `mature` | applied if supplied (maturity ↑ toward 1) | applied if supplied (bound_fraction ↑ toward 1) | not applied | yes if `bound_fraction > 0` |
| `slipping` | not applied (maturity stays; structural compromise) | not applied (Guard 2: slipping never increases bound) | applied if supplied (bound_fraction ↓ toward 0) | yes if `bound_fraction > 0`; will trend to zero under `k_unbind` |
| `released` | not applied | not applied | not applied | no — terminal state, traction forced to (0, 0) |

Automatic state-label transitions (e.g., `nascent → mature` when
`maturity` crosses a threshold; `mature → slipping` under stress)
are explicitly **out of scope** for 6.3a. The state label only
changes if the caller assigns a new state on the input adhesion;
the dynamics step does not modify it.

### Status

PASS — every sign is derivable from the locked formula and
testable by single-FA direct check.

---

## 6. Measurement-protocol consistency (Hard Rule 11)

This unit produces:

- `traction_force_nN_xy[i]` as a per-FA xy vector in `nN`,
  conventionally **cell-on-substrate**.
- `substrate_reaction_nN_xy[i]` as the negation, conventionally
  **substrate-on-cell**.
- `radial_components[i]` (scalar) and `tangential_components[i]`
  (scalar) of the traction with respect to the supplied centroid,
  with the convention `radial = traction · n_hat_radial_outward`
  (so a default radial-inward traction has `radial < 0`).

### Round-trip contract

`reconstructed = radial_components[i] · n_hat_radial_outward + tangential_components[i] · n_hat_tangential_ccw`
must equal `traction_force_nN_xy[i]` to within float64 round-off.

### Schema docstring patch

`acs/v2/focal_adhesion.py` `FocalAdhesionState.traction_force_nN_xy`
docstring is updated in the same commit to read:

> ``traction_force_nN_xy``: per-FA cell-on-substrate traction
> vector in nN. The substrate-on-cell reaction is the negation;
> Newton 3 holds per FA. Convention locked in
> `docs/v2_focal_adhesion_dynamics_sanity_gate.md` §6.
>
> State conventions (locked in the same gate, §5 state table):
>
> - ``state == "unbound"``: ``bound_fraction`` is forced to 0;
>   traction is forced to (0, 0). The existing schema validator
>   (``no_traction_when_unbound``) already enforces this on
>   construction.
> - ``state == "released"``: terminal absorbing state; the 6.3a
>   dynamics step forces traction to (0, 0) regardless of stored
>   ``bound_fraction``. Schema permits ``bound_fraction > 0`` to
>   record a "what was bound right before release" history value;
>   downstream consumers must read traction, not bound_fraction,
>   for the active force.
> - ``state == "slipping"``: traction follows §5 magnitude
>   formula but ``bound_fraction`` cannot increase (decreases only
>   under explicit ``k_unbind_per_s``).
> - ``state in {"nascent", "mature"}``: traction follows §5
>   formula; both maturity and bound_fraction can increase under
>   explicit caller-supplied rates.

This patch is contract-defining; consumers downstream (frame_dump,
viz, future ECM-OL coupling) read the same field with the same
sign.

### Status

PASS — every per-FA force has a documented direction and a
round-trip-tested decomposition.

---

## 7. Magic-Number Block

| Symbol | Source | Block status |
|---|---|---|
| `traction_scale_nN` | runtime config (no project default) | not a magic number — caller-supplied, fail to validate if absent. |
| `dt_fa_s` | runtime config (no project default) | not a magic number — caller-supplied, fail to validate if absent. |
| `k_maturity_per_s`, `k_bind_per_s`, `k_unbind_per_s` | optional runtime config (no project default) | not a magic number — absent means "no rate update" (explicit zero). |
| `max_traction_nN` | optional validation-only upper bound | not a magic number — never clamps; raises on exceed. |
| `traction_axis_xy` (when supplied) | runtime config | not a magic number — caller direction. |
| `_DT_RATE_SAFETY_MARGIN = 0.5` | reused from active contour | named module-level constant; rationale identical to the active-contour Sanity Gate (linear-stability `<2` / energy-monotone `≤1` / factor-of-2 nonlinear safety `≤0.5`); pre-run, not chosen to make tests pass. |

All three Magic-Number Block tests pass on every numeric:
1. **Derivable**: every value is either runtime config or the
   reused 0.5 stability margin.
2. **Grid-invariant**: trivially — no grid in this unit.
3. **Not fitting**: every value was chosen *before* any test was
   run; none was tuned to pass a gate.

### Status

PASS — Magic-Number Block clean. The lock pre-removed the natural
gate-tuning temptations by forbidding default rates and clipping.

---

## 8. Gate verdict

§1 dimensional, §2 boundary, §3 conservation, §4 numerical, §5
sign, §6 measurement-protocol — all PASS. §7 Magic-Number Block
PASS.

The gate clears for executable code in two small commits or one
combined commit:

- `acs/v2/dynamics/focal_adhesion.py` (the static FA dynamics
  module).
- `acs/v2/focal_adhesion.py` (schema docstring patch for the
  cell-on-substrate sign convention).
- `tests/test_v2_focal_adhesion_dynamics.py` (the Sanity Gate's
  testable items).

If Codex review surfaces a missing reduction, hidden tunable, or
sign-convention drift, the gate flips to BLOCKER with the
three-options template (defer 6.3a, narrow scope further, switch
to a different next unit).

### Outstanding before code lands

- Codex review of this gate (the FA-side spec match for the
  design lock id=1076).
- A single docs commit landing this gate file alongside the code
  is acceptable. The locked design discussion (already on
  `docs/v2_phase1_forward_roadmap.md` and the design-discussion
  ledger) is the spec source; this file is the impl-side gate.

PI awareness: this unit is mid-cycle; PI will see commit hashes +
fast-suite test counts as the cadence cycles complete.
