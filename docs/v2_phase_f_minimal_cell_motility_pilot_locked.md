# V2 Phase F Minimal Cell Motility Pilot — Locked Design (A-tier-compressed per PI ASAP)

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (A-tier-compressed per PI
`id=1903` ASAP directive — 2 design rounds + Codex single review pass +
seal, NOT 3-round full A-tier debate)
**Source unit**: design-discussion
`topic=v2-phase-f-minimal-cell-motility-pilot`, MCP id 1914–1938
(round 1 / Codex review C1-C6 / round 2 Z1-Z6 corrections / seal ack)
**Parent PI directive**: PI `id=1894` ("세포 움직이는 것까지 최대한
빠르게 보고싶다") + PI `id=1903` (α 채택)
**Upstream sister locks**:
- `docs/v2_p1_derivation_locked.md` (P1 active contour; Phase F-local
  helper reuses P1 internal force functions WITHOUT modifying locked P1
  step signature — Codex Z2/Z5)
- `docs/v2_phase_e_v2_composition_step_2_locked.md` (Phase E v2; called
  per Phase F step for ECM-side diagnostics in the loop)
- `docs/v2_hard_blocker_4_active_locked.md` (HB#4-active; multipliers
  **diagnostic-only here** — Codex C4 critical guard against magic-number
  rate-to-force conversion)
- `docs/v2_hard_blocker_5_lyapunov_metric_locked.md` (HB#5; v_active
  diagnostic carried via Phase E v2 result)
- `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md` (HB#1+#2)
- `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md` (HB#3)
- `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` (HB#4 v1
  neutral, sister with active)
**Tier**: **A-tier-compressed**. New physics (FA-vertex Newton 3rd
coupling) but PI ASAP directive justifies compressed cycle (2 rounds +
single review). Future Phase F proper for ECM→cell feedback is a
separate A-tier full-debate cycle.
**PI ratify status**: full delegation per PI `id=1903`. impl-work uses
this for the Phase F minimal pilot Sanity Gate doc + code entry.

---

## 0. Scope + Wording Boundary (Z6 honest-scope anchor)

### What this pilot IS

- Minimal Phase F cell motility pilot demonstrating **FA-traction-driven
  contour vertex motion** in v2 architecture
- Reuses P1 active contour mechanics (cortex + area forces) +
  Phase E v2 wrapper (ECM-side diagnostics in the loop)
- New physics: single layer = FA traction → Newton's 3rd law reaction
  on attached cell vertex → vertex displacement via P1-style overdamped
  Euler step
- Output: visible cell-boundary motion in v2 mp4 (PI ASAP deliverable)

### What this pilot is NOT (Codex Z6 wording verbatim)

> "This pilot shows FA-traction-driven contour motion with Phase E v2
> ECM-side diagnostics running in the loop. It **does not yet implement
> ECM→cell motility feedback**, because HB#4 rate multipliers are
> **diagnostic-only here**."

- ❌ NOT ECM→cell motility feedback (HB#4 rate multipliers strictly
  diagnostic-only; using them as traction-force scalers would be a
  hidden new constitutive law / magic-number — Codex C4/Y4 critical
  catch)
- ❌ NOT FA dynamics (assembly/disassembly not modeled; FA position
  follows attached vertex but FA count/state remain fixed)
- ❌ NOT durotaxis or alignment-based motility — those require explicit
  ECM stiffness/orientation → vertex coupling laws (Phase F proper,
  separate A-tier full-debate cycle)
- ❌ NOT validated against PI experimental data (synthetic pilot only;
  Hard Rule 1)
- ❌ NOT a `satisfies_*` claim — test names use `provides_*_evidence` /
  scoped-fixture naming (Phase E v1 Y1+Y2 sister discipline)

---

## 1. Final Lock Summary

### Function signature

```python
"""V2 Phase F minimal cell motility pilot — first cell-moving v2 sample.

Couples P1 active contour + Phase E v2 wrapper via FA-vertex external
force (Newton's 3rd law reaction). HB#4 rate multipliers are
diagnostic-only; ECM→cell motility feedback is NOT implemented here
(future Phase F proper, separate A-tier cycle).
"""

# acs/v2/dynamics/phase_f_minimal_motility.py (NEW module)

import numpy as np
from dataclasses import dataclass, replace

from acs.v2.active_contour import (
    _DT_RATE_SAFETY_MARGIN,
    ActiveContourState,
    compute_rate_max,
)
from acs.v2.dynamics.active_contour import (
    ActiveContourStepError,
    compute_area_forces,
    compute_cortex_forces,
)
from acs.v2.dynamics.closed_loop_phase_e import (
    PhaseEStepResult,
    step_closed_loop_phase_e_v2,
)
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState


@dataclass(frozen=True, slots=True)
class PhaseFStepResult:
    """One Phase F minimal motility step result.

    HB#4 multipliers in `phase_e_v2.ecm_to_fa_bias.diagnostics_dict`
    are diagnostic-only; no traction-force scaling (Y4/Z5 guard).
    """
    contour: ActiveContourState                              # post-step P1 state
    adhesions: tuple[FocalAdhesionState, ...]                # FA positions follow attached vertices
    phase_e_v2: PhaseEStepResult                             # ECM-side diagnostics (in-the-loop)
    updated_ecm: ECMSubstrateState                           # = phase_e_v2.updated_ecm (identity invariant)


def step_phase_f_minimal_motility(
    contour_state: ActiveContourState,
    adhesions: tuple[FocalAdhesionState, ...],
    ecm: ECMSubstrateState,
    fa_to_vertex_index: dict[str, int],
    *,
    k_active: float,
) -> PhaseFStepResult:
    """One synchronous Phase F minimal motility step.

    Uses contour_state.params.dt_cell_s as the synchronous dt for both
    P1 active contour and Phase E v2 wrapper (no subcycling — future
    Phase F proper). HB#4 rate multipliers are diagnostic-only; FA
    traction drives cell vertex motion through Newton's 3rd law
    reaction directly (no multiplier scaling).
    """
    # Y1 (Z1 corrected fields): validate fa_to_vertex_index fail-closed
    n_vertices = contour_state.vertices_xy_um.shape[0]
    _validate_fa_to_vertex_index(fa_to_vertex_index, adhesions, n_vertices)

    # Synchronous small dt (Y3): use P1's dt_cell_s, NOT a separate dt argument
    dt_s = float(contour_state.params.dt_cell_s)

    # 1. Phase E v2 wrapper (ECM diagnostics in the loop, multipliers diagnostic-only)
    pe_v2_result = step_closed_loop_phase_e_v2(
        adhesions, ecm, dt_s, k_active=k_active,
    )

    # 2. Build per-vertex external force (Y2 sign + Y4 NO multiplier scaling + Y10 wrapper seam)
    external_forces = np.zeros((n_vertices, 2), dtype=np.float64)
    for fa in adhesions:
        v_idx = fa_to_vertex_index[fa.adhesion_id]
        # Newton's 3rd: cell vertex receives -F where F = cell-on-substrate traction
        # NO HB#4 multiplier scaling (diagnostic-only — Y4 critical guard)
        external_forces[v_idx] += -np.asarray(fa.traction_force_nN_xy, dtype=np.float64)

    # 3. Phase F-local helper — reuse P1 internal forces + add validated external (Y5)
    contour_new = _step_active_contour_with_external_force(
        contour_state, external_forces_per_vertex_nN=external_forces,
    )

    # 4. Update FA positions to follow attached vertices (Y6 corrected field name)
    adhesions_new = []
    for fa in adhesions:
        v_idx = fa_to_vertex_index[fa.adhesion_id]
        new_position = tuple(contour_new.vertices_xy_um[v_idx].tolist())
        adhesions_new.append(replace(fa, position_um_xy=new_position))

    return PhaseFStepResult(
        contour=contour_new,
        adhesions=tuple(adhesions_new),
        phase_e_v2=pe_v2_result,
        updated_ecm=pe_v2_result.updated_ecm,
    )


def _step_active_contour_with_external_force(
    contour_state: ActiveContourState,
    external_forces_per_vertex_nN: np.ndarray,
) -> ActiveContourState:
    """Phase F-local helper (Y5). Same overdamped Euler update as P1 step,
    plus external force per vertex. Reuses locked compute_cortex_forces /
    compute_area_forces / compute_rate_max without modifying P1 step
    signature.

    Failure semantics sister-consistent with P1:
    raises ActiveContourStepError("dt_violation") on dt_rate_product
    safety margin breach (Y7).
    """
    # Y7 (Z2 corrected): P1 sister-consistent failure
    rate_info = compute_rate_max(contour_state)
    if rate_info["dt_rate_product"] > _DT_RATE_SAFETY_MARGIN:
        raise ActiveContourStepError(
            "dt_violation",
            f"dt_cell_s={contour_state.params.dt_cell_s}·rate_max product "
            f"{rate_info['dt_rate_product']:.4f} > safety margin "
            f"{_DT_RATE_SAFETY_MARGIN}",
        )

    # Internal forces (P1 reused, no modification)
    internal = compute_cortex_forces(contour_state) + compute_area_forces(contour_state)

    # Validate external force shape
    if external_forces_per_vertex_nN.shape != internal.shape:
        raise ValueError(
            f"external_forces_per_vertex_nN shape {external_forces_per_vertex_nN.shape} "
            f"!= internal force shape {internal.shape}"
        )

    # Total force
    forces = internal + external_forces_per_vertex_nN

    # Same overdamped Euler update as P1 step
    zeta = rate_info["zeta_nN_s_per_um"]
    new_vertices = (
        contour_state.vertices_xy_um
        + contour_state.params.dt_cell_s * forces / zeta[:, None]
    )

    return replace(contour_state, vertices_xy_um=new_vertices)


def _validate_fa_to_vertex_index(
    fa_to_vertex_index: dict[str, int],
    adhesions: tuple[FocalAdhesionState, ...],
    n_vertices: int,
) -> None:
    """Y1 fail-closed validation of FA-to-vertex mapping."""
    fa_ids = {fa.adhesion_id for fa in adhesions}
    mapping_ids = set(fa_to_vertex_index.keys())
    unknown = mapping_ids - fa_ids
    if unknown:
        raise ValueError(
            f"fa_to_vertex_index references unknown FA ids: {unknown}"
        )
    missing = fa_ids - mapping_ids
    if missing:
        raise ValueError(
            f"fa_to_vertex_index missing FA ids: {missing}"
        )
    for fa_id, v_idx in fa_to_vertex_index.items():
        if isinstance(v_idx, bool):
            raise ValueError(
                f"vertex index for FA {fa_id} must be int (not bool — "
                f"Python bool ⊂ int trap), got {type(v_idx).__name__}"
            )
        if not isinstance(v_idx, int):
            raise ValueError(
                f"vertex index for FA {fa_id} must be int, "
                f"got {type(v_idx).__name__}"
            )
        if not (0 <= v_idx < n_vertices):
            raise ValueError(
                f"FA {fa_id} vertex index {v_idx} outside [0, {n_vertices})"
            )
```

### Forbidden in Phase F minimal pilot

- HB#4 rate multipliers used as traction-force scalers in any form
  (Codex C4/Y4 — critical guard against hidden new constitutive law /
  magic-number; multipliers strictly **diagnostic-only**, monkeypatch
  test 8 enforces at wrapper seam)
- Sign-flipping FA traction force to make movie look better (Codex
  Y2 — schema convention `traction_force_nN_xy` is cell-on-substrate;
  Newton 3rd reaction on cell is `-traction`. For visible +x motion,
  set fixture traction to `(-1, 0)` and document)
- `dt_s` argument to Phase F step (Codex Y3 — uses
  `contour_state.params.dt_cell_s` synchronously; no subcycling in
  minimal pilot)
- Patching P1 active_contour `step` signature (Codex Y5 — Phase F-local
  helper reuses internal force functions without modifying P1 lock)
- New `dt_violation` failure semantics (Codex Y7 — reuses
  `ActiveContourStepError("dt_violation", ...)` from P1)
- Field name `vertices_um_xy` (Codex Z1 — actual schema is
  `vertices_xy_um`)
- `dt_cell_s = 0` (Codex Z3 — `ActiveContourParameters.validate()`
  requires strictly positive; for no-motion test use
  `lambda_c=0, k_a=0, traction=(0,0)` instead)
- ECM→cell direct coupling (durotaxis, alignment-based motility, etc.)
  — separate Phase F proper A-tier full-debate cycle
- FA dynamics (assembly/disassembly) — separate cycle (currently 6.3a
  static / 6.3b protrusion-coupled both treat FA as input)
- "satisfies_*" or "Items 1-4 satisfied" wording (Phase E v1 Y1+Y2
  sister discipline; use `provides_*_evidence` / scoped-fixture naming)
- PI experimental data fitting (Hard Rule 1)

---

## 2. Reasoned-acceptance trace (Y1–Y11)

The lock converged after 2 design rounds + 1 Codex final review pass +
seal ack (A-tier-compressed per PI ASAP).

**Y1 (round 2, accepted Codex C1)**: Q1=(b) attached vertex with
fail-closed validation (`fa_to_vertex_index` fail on unknown FA id /
out-of-range vertex / bool-int trap).

**Y2 (round 2, accepted Codex C2)**: Newton 3rd sign + corrected
fixture direction. `traction_force_nN_xy` is cell-on-substrate per
`acs/v2/focal_adhesion.py`; Newton reaction on cell = `-traction`. For
visible +x cell motion, fixture FA traction = `(-1.0, 0.0)`.
Sign-flip-for-cosmetic-effect explicitly forbidden.

**Y3 (round 2, accepted Codex C3)**: Synchronous small `dt_cell_s` from
P1 contour state, NOT separate `dt_s` argument. P1 enforces `dt_cell ·
rate_max ≤ 0.5`; Phase E v2 sub-call uses same small dt (no
subcycling). Pilot config: `dt_cell_s=1e-3 s, n_steps=10000,
frame_interval=100` for visible cumulative motion.

**Y4 (round 2, accepted Codex C4 — CRITICAL)**: HB#4 rate multipliers
strictly diagnostic-only. Mean-of-multipliers-as-traction-scaler would
be a hidden new constitutive law (rate-to-force conversion factor =
magic number; Magic-Number Block violations across all 3 tests). Lock
forbids any multiplier-based force scaling. ECM→cell feedback deferred
to Phase F proper (separate A-tier cycle).

**Y5 (round 2, accepted Codex C5)**: Phase F-local helper
`_step_active_contour_with_external_force` reuses
`compute_cortex_forces`, `compute_area_forces`, `compute_rate_max`
from P1 without patching P1 `step` signature. Phase F module owns the
new physics; P1 lock untouched.

**Y6 (round 3, accepted Codex Z1)**: `vertices_xy_um` field name
(NOT `vertices_um_xy`). Verified at `acs/v2/active_contour.py:241`.
Pseudocode + tests + FA-position follow logic all updated.

**Y7 (round 3, accepted Codex Z2)**: `compute_rate_max` returns
`dt_rate_product` + `margin_ratio` + `zeta_nN_s_per_um` (NOT
`dt_violation`). P1 sister-consistent failure: reuse
`ActiveContourStepError("dt_violation", ...)` and
`_DT_RATE_SAFETY_MARGIN` from P1.

**Y8 (round 3, accepted Codex Z3)**: `dt_cell_s = 0` invalid per
`ActiveContourParameters.validate()` (strictly positive). Replace
`zero_dt_no_movement` test with `lambda_c=0 + k_a=0 + finite positive
dt_cell_s + traction (0,0) → no motion`.

**Y9 (round 3, accepted Codex Z4)**: ECM evolution test split into:
(a) invariant test (Phase E v2 always invoked, identity invariant
holds) — always asserted; (b) fixture-specific evolution test — only
asserted for fixtures with guaranteed-nonzero HB#1+#2 update. Don't
promise broad ECM evolution invariant.

**Y10 (round 3, accepted Codex Z5)**: Y4 guard test monkeypatches at
**wrapper seam** (`step_closed_loop_phase_e_v2` in Phase F module
namespace), NOT deep internals (`compute_ecm_to_fa_bias_active`).
Robust against import/cache patterns. Two runs with identical ECM but
neutral vs huge multipliers must produce IDENTICAL vertex displacement.

**Y11 (round 3, accepted Codex Z6)**: Honest-scope wording verbatim
in module docstring + function docstring + return-class docstring:
"This pilot shows FA-traction-driven contour motion with Phase E v2
ECM-side diagnostics running in the loop. It does not yet implement
ECM→cell motility feedback, because HB#4 rate multipliers are
diagnostic-only here." Forward guard against future audit/PI
misinterpretation.

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

1. **Units**:
   - `vertices_xy_um`: μm (per P1 schema)
   - `traction_force_nN_xy`: nN (cell-on-substrate per FA schema)
   - `external_forces_per_vertex_nN`: nN (Newton reaction on vertex)
   - `dt_cell_s`: s (per P1 schema, strictly positive)
   - `zeta_nN_s_per_um`: nN·s/μm (per P1)
   - vertex displacement = `dt · forces / zeta` → s · nN / (nN·s/μm) = μm ✓
2. **Boundary**:
   - Empty FA list: `external_forces` all-zero → P1 step alone (cortex
     + area dynamics)
   - Zero traction (FA exists but `traction_force_nN_xy=(0,0)`): no
     external force; P1 dynamics only
   - Single FA single vertex: simplest minimal pilot setup
   - All vertices have FAs: each vertex receives one Newton reaction;
     summation handles correctly
   - dt safety margin breach: `ActiveContourStepError("dt_violation")`
     raised before vertex update (P1 sister-consistent)
   - Invalid `fa_to_vertex_index`: fail-closed `ValueError` (Y1)
3. **Conservation**: P1 active contour conserves nothing (overdamped);
   external force adds momentum input from FA-substrate interaction
   (Newton's 3rd law). Phase E v2 ECM-side conservation per its lock.
4. **Numerical**:
   - `np.float64` throughout (matches P1 + Phase E v2)
   - Synchronous small `dt_cell_s` ~1e-3 s; many steps + frame_interval
     for visible motion
   - No subcycling in minimal pilot (Y3); future Phase F proper
5. **Sign**:
   - Newton's 3rd law: `cell_vertex_external_force = -fa.traction_force`
   - `traction_force=(1,0)` → `vertex displacement -x` (cell pulls
     substrate +x, substrate pulls cell -x)
   - For visible +x movie: fixture sets `traction=(-1,0)`, with
     explicit docstring (Y2)
6. **Measurement-protocol** (Hard Rule 11, the central anchor):
   - **Y4 magic-number guard**: HB#4 rate multipliers
     (`k_maturity_per_s`, `k_bind_per_s`, `k_unbind_per_s`) are NOT
     traction-force multipliers. Any mean/sum/product reinterpretation
     would be a new constitutive law requiring its own derivation +
     literature anchor. Lock test 8 (wrapper seam monkeypatch) enforces
     "huge multipliers vs neutral → identical vertex displacement".
   - **Y2 sign convention**: per `acs/v2/focal_adhesion.py` schema,
     `traction_force_nN_xy` is cell-on-substrate. Reaction direction
     fixed by Newton's 3rd law, not by movie aesthetics.
   - **Y11 honest-scope wording**: module/function docstrings
     explicitly state "ECM→cell motility feedback NOT implemented here".
     Future audit / PI review cannot misinterpret minimal pilot mp4 as
     full ECM-driven motility.
   - **Phase E v2 v1+v2 sister-pattern**: tests use
     `provides_*_evidence` / scoped-fixture naming; meta-test enforces
     no `satisfies_*_item` test names (Phase E v1 Y2 forward guard).

---

## 4. Test Catalog (12 tests)

### Composition + validation (4)

1. `test_phase_f_step_returns_phase_f_step_result_dataclass` — return
   is `PhaseFStepResult` with all 4 fields populated; identity
   invariant `result.updated_ecm is result.phase_e_v2.updated_ecm`
2. `test_phase_f_step_unknown_fa_id_raises` (Y1) —
   `fa_to_vertex_index` references FA id not in adhesions → `ValueError`
3. `test_phase_f_step_vertex_index_out_of_range_raises` (Y1) — vertex
   index outside `[0, n_vertices)` → `ValueError`
4. `test_phase_f_step_dt_violation_raises_active_contour_step_error`
   (Y7) — `dt_cell_s · rate_max > _DT_RATE_SAFETY_MARGIN` → reuses
   P1 `ActiveContourStepError("dt_violation", ...)` (NOT a
   wrapper-specific failure kind)

### Sign + force coupling (3)

5. `test_phase_f_step_traction_one_zero_moves_vertex_negative_x`
   (Y2 + isolated `lambda_c=0, k_a=0`) — fixture FA traction `(1, 0)`
   with zero internal force; assert vertex moves to **−x** (Newton
   3rd reaction on cell-on-substrate convention)
6. `test_phase_f_step_zero_traction_and_zero_internal_force_no_movement`
   (Y8 — replaces invalid `dt_cell_s=0` test) —
   `lambda_c=0, k_a=0, dt_cell_s=1e-3, traction=(0,0)` → vertex
   positions unchanged within float64 tolerance
7. `test_phase_f_step_internal_force_alone_relaxes_per_p1` — sanity:
   external force zero + nonzero `lambda_c, k_a` → vertex motion
   matches P1 active_contour `step` directly (helper is sister-consistent)

### HB#4 multipliers diagnostic-only guard (1, Y4 + Y10 wrapper seam)

8. `test_phase_f_step_multipliers_diagnostic_only_not_used_for_force_scaling`
   (Y10 wrapper seam) — monkeypatch `pf.step_closed_loop_phase_e_v2`
   to return two results with identical ECM but neutral vs huge
   multipliers (1e6); assert vertex displacement IDENTICAL within
   `atol=1e-15`. Forward guard against silent magic-number injection.

### FA position + ECM evolution (2)

9. `test_phase_f_step_fa_position_follows_attached_vertex` — after
   step, `result.adhesions[i].position_um_xy ==
   tuple(result.contour.vertices_xy_um[v_idx])` for the attached
   vertex
10. `test_phase_f_step_phase_e_v2_invoked_and_ecm_carried` (Y9
    invariant) — `result.phase_e_v2 is not None`,
    `result.updated_ecm is result.phase_e_v2.updated_ecm` (identity
    invariant inherited from Phase E v2 Y4)

### Fixture-specific evolution (1, Y9 — Codex bookkeeping ack)

11. `test_phase_f_step_ecm_orientation_evolves_under_nonzero_traction_fixture`
    (Y9) — fixture: nonzero FA traction, isotropic `0.5*I` IC,
    `k_active=1.5, dt_cell_s=1e-3`. Assert `np.abs(result.updated_ecm.orientation_tensor
    - ecm.orientation_tensor).max() > 0.0` (at least one cell evolved
    via HB#1+#2). Locked as required test per Codex `id=1938`
    bookkeeping (12 tests total, not optional).

### Exports (1)

12. `test_phase_f_step_exports_through_both_init` —
    `step_phase_f_minimal_motility` + `PhaseFStepResult` importable
    from BOTH `acs.v2.dynamics` AND `acs.v2`; private helpers
    (`_step_active_contour_with_external_force`,
    `_validate_fa_to_vertex_index`) NOT exported.

---

## 5. Files

- `docs/v2_phase_f_minimal_cell_motility_pilot_locked.md` (this file)
- `acs/v2/dynamics/phase_f_minimal_motility.py` (NEW):
  - `PhaseFStepResult` dataclass (frozen, slots, 4 fields)
  - `step_phase_f_minimal_motility` function
  - `_step_active_contour_with_external_force` private helper (Y5)
  - `_validate_fa_to_vertex_index` private helper (Y1)
  - Imports: P1 internal force functions + Phase E v2 wrapper +
    schemas
  - **NO** `dt_s` argument to public function (Y3 synchronous)
  - **NO** HB#4 multiplier reads in force computation (Y4/Z5)
  - **NO** P1 step signature modification (Y5)
- `acs/v2/dynamics/__init__.py` + `acs/v2/__init__.py` — 2 new exports
  (`PhaseFStepResult`, `step_phase_f_minimal_motility`)
- `tests/test_v2_phase_f_minimal_motility.py` (NEW, 12 tests per §4)
- (Future, separate cycle) `acs/v2/phase_f_minimal_motility_pilot_runner.py`
  — sister with `phase_e_v2_pilot_runner.py`. Either same impl cycle
  (B-tier compressed bundled with this lock's impl) or follow-up cycle
  (impl-work decision per dispatch).

---

## 6. References

- Brief / round 1 / round 2 / seal ack: design-discussion MCP id
  1914–1938
- PI directives:
  - `id=1894` "세포 움직이는 것까지 최대한 빠르게 보고싶다"
  - `id=1903` "알파로 진행 부탁" (option α from id=1899)
  - `id=1653` PI risk-tier policy (B-tier compressed default; A-tier
    full debate exception, here compressed by ASAP exception)
- Sister locks (all 5 HB + Phase D + B1 + Phase E v1+v2 +
  HB#4-active):
  - `docs/v2_p1_derivation_locked.md` (P1 active contour mechanics
    — Y5 reuse without modification)
  - `docs/v2_phase_e_v2_composition_step_2_locked.md` (Phase E v2
    wrapper — Y4 multipliers diagnostic-only)
  - `docs/v2_hard_blocker_4_active_locked.md` (HB#4-active step 1 —
    Y4 magic-number guard)
  - `docs/v2_hard_blocker_5_lyapunov_metric_locked.md`
  - `docs/v2_hard_blocker_1_2_constitutive_direction_locked.md`
  - `docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
- Verified at source for Y6/Y7/Y8 (P1 schema):
  - `acs/v2/active_contour.py:241` (`vertices_xy_um` field name)
  - `acs/v2/dynamics/active_contour.py:53`
    (`ActiveContourStepError` class)
  - `acs/v2/dynamics/active_contour.py:136` (P1 step uses
    `_DT_RATE_SAFETY_MARGIN`, `dt_rate_product` key)
- v1 reference (cell motility precedent):
  - `runs/20260503T1833Z_p1_alpha/test_4_coupled_ellipse/summary.html`
    (P1 active contour ellipse → circular convergence dynamics —
    validated; Phase F minimal pilot extends this with FA-vertex
    Newton 3rd coupling)
- Hard Rule 1 (no PI data fitting): CLAUDE.md
- Hard Rule 11 (sim-experiment measurement matching): CLAUDE.md +
  memory `hard_rule_11_wording_boundary_meta_test.md`
- Magic-Number Block: CLAUDE.md (Y4 critical guard against
  rate-to-force conversion magic number)
- Phase E v1 Y1+Y2 evidence-not-satisfaction discipline: lock
  `docs/v2_phase_e_composition_locked.md` (Y11 sister)

---

## 7. Cross-room dispatch

This file is the design-team input to implementation-work for:

1. impl Claude implements `step_phase_f_minimal_motility` per §1
   pseudocode + §4 12-test catalog directly. **Skip Sanity Gate doc**
   per A-tier-compressed mode (PI ASAP).
2. impl Codex performs **single review pass** on the resulting
   commit covering:
   - Y4 critical guard: HB#4 multipliers strictly diagnostic-only;
     test 8 (wrapper seam monkeypatch) verifies neutral-vs-huge
     multipliers produce identical vertex displacement
   - Y2 sign convention: fixture traction `(1, 0)` produces vertex
     displacement to **−x** (Newton 3rd reaction)
   - Y3 synchronous small dt: `dt_s = contour_state.params.dt_cell_s`,
     no separate `dt_s` argument
   - Y5 P1 untouched: Phase F-local helper reuses P1 internal force
     functions; P1 `step` signature unchanged
   - Y6/Y7/Y8 schema correctness: `vertices_xy_um` field,
     `ActiveContourStepError("dt_violation")` failure semantics,
     `dt_cell_s = 0` rejected at P1 schema layer
   - Y11 honest-scope wording: module/function/return-class docstrings
     contain explicit "does not yet implement ECM→cell motility
     feedback" phrase
3. On PASS: Phase F minimal pilot milestone reached.
4. After milestone: design-discussion idle on Phase F minimal pilot.
   Next forward units (per PI ASAP follow-on):
   - **Pilot runner** (B-tier compressed bundled or separate):
     `acs/v2/phase_f_minimal_motility_pilot_runner.py` sister with
     `phase_e_v2_pilot_runner.py`; CLI; persistent run with HDF5
     frames + dashboard + replay + mp4 (reuse Cycles 13-15
     visualization platform); PI mp4 with cell visibly moving
   - **Phase F proper** (separate A-tier full-debate cycle): ECM→cell
     motility feedback; constitutive law for HB#4 multipliers (or
     ECM stiffness/orientation directly) → vertex force; literature
     anchors required; explicit Magic-Number Block compliance

Compressed cycle MCP id 1914–1938 (round 1 / Codex review C1-C6 /
round 2 Z1-Z6 / seal ack).
