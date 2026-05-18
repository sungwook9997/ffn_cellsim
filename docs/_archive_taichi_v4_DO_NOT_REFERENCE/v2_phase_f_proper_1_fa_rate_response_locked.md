# V2 Phase F Proper-1: ECM-Modulated FA Rate Response — Locked Design

> ⚠️ **STATUS: WIP HALTED** (added 2026-05-06 KST after Option A+B
> restructure post-merge audit, claude-work `mcp_msg:2550`, codex
> approval `mcp_msg:2551`).
>
> Implementation coherent but untested; lock §4/§5 test catalog is
> not implemented; Phase F proper resume requires explicit PI
> directive after Option 1 pivot. The module
> `acs/v2/dynamics/fa_rate_response.py` was landed as "coherent,
> untested" in commit `a1d97e8` per Codex `id=2196` STOP directive
> (PI Option 1 pivot — vertical 5-layer framework completion takes
> priority over Phase F follow-ups). Lock `89d6ddc` remains useful
> as design / implementation reference and is **not retracted** —
> only halted before test implementation. Imports / exports through
> `acs/v2/dynamics/__init__.py` and `acs/v2/__init__.py` are wired
> but the function `step_fa_rate_response` has zero direct test
> coverage; consumers must not treat it as execution-ready until PI
> reopens the mechanism path and the §4 17-test catalog lands.

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (4-round adversarial lock,
PI id=809 + id=1985 aggressive debate posture)
**Source unit**: design-discussion `topic=v2-phase-f-proper-design`,
MCP id 2127–2168 (rounds 1–4 + seal ack)
**PI directives**:
- `id=2099` "이걸로는 아무것도 알 수가 없음 너무 부족!!!" (visual evidence inadequate)
- `id=2114` external audit (plumbing-vs-mechanism gap; honest correction)
- `id=2123` "바로 시작" (immediate start of (X+Y) parallel routing)
- `id=1985` "건강하고 활발한 토론" (aggressive debate posture)
**Tier**: **A-tier full debate** (4 rounds; 19 substantive Codex catches
across rounds 1-4, audit-driven entry).
**Stage**: F-proper-1 of the 3-stage Phase F proper sequence (per Codex
`id=2147` Y3 staging):
- F-proper-1 (this lock): **ECM-modulated FA rate dynamics** (per-rate
  HB#4 multiplier scaling; delegate traction algebra to 6.3a; no
  contour motion)
- F-proper-2 (separate lock): multi-rate motility integrator (FA
  traction → contour vertex displacement; subcycling)
- F-proper-3 (separate lock): PI-readable 15-min annotated pilot
  (visible deformation, annotated movie)
**PI ratify status**: full delegation per PI id=939/1008/2123. impl-work
uses this for the F-proper-1 implementation.

---

## 0. Scope + Wording Boundary (audit `id=2114` directive)

### What this unit IS

- F-proper-1 = **first stage of Phase F proper** (real cell motility,
  NOT the minimal pilot's plumbing demo)
- ECM-modulated FA rate dynamics: HB#4-active multipliers act as
  **per-rate scalers** (k_maturity_per_s, k_bind_per_s, k_unbind_per_s)
- Delegate traction algebra to existing 6.3a `step_focal_adhesions_static`
  (Y6 — leverage locked physics, NO new constitutive law)
- Manual aggregation of per-FA results (Y16 — NO second delegate call;
  prevents twice-mutation corruption)
- Pure read-only-of-FA, single step (no contour motion, no integrator,
  no subcycling — those belong to F-proper-2)

### What this unit is NOT (audit-driven scope discipline)

- ❌ NOT cell motility (visible cell motion requires F-proper-2
  integrator + F-proper-3 long-horizon pilot)
- ❌ NOT a patch to `step_phase_f_minimal_motility` (Phase F minimal Y4
  guard preserved forever — Y4 here)
- ❌ NOT new HB lock (Y5 — uses existing HB#4-active multipliers as
  per-rate scalers per their original `RATE_NAMES` design intent)
- ❌ NOT new FA state-label transitions (Y15 — 6.3a updates
  `maturity`/`bound_fraction` only; nascent→mature/released auto-labels
  require separate lock)
- ❌ NOT FA spawning / nucleation / RNG (Y5 — `k_nucleation_per_s` etc.
  introduce site-selection / stochasticity, separate lock)
- ❌ NOT direct multiplier-to-force scaling (Y4 minimal sister + Y6 —
  force change comes through 6.3a algebra `traction = traction_scale_nN
  · maturity · bound_fraction · unit_axis`)
- ❌ NOT contour subcycling / multi-rate integrator (Y8 — F-proper-2
  owns)
- ❌ NOT "boundary-normal outward" axis wording (Y10 sign trap —
  outward cell-on-substrate would CONTRACT the contour; spreading
  requires INWARD cell-on-substrate → outward substrate-on-cell
  reaction)
- ❌ NOT silent multiplier reorder (Y14 — fa_ids exact match required)
- ❌ NOT PI 260313 fitting (Hard Rule 1)
- ❌ NOT magic-number rate constants (Y8 — literature anchor required;
  Cho 2020 + general epithelial spheroid FA turnover)

### Audit-driven framing wording (`id=2114` central anchor)

This unit is the **first real cell-motility step** in v2 architecture.
Phase F minimal pilot was honestly labeled `not_mechanistic: true` in
its metadata; F-proper-1 begins the genuine mechanism path. F-proper-1
alone does NOT produce visible cell motion in PI-readable mp4 — that
requires F-proper-2 (multi-rate integrator) and F-proper-3 (15-min
annotated pilot, visible deformation ≥ 0.1 μm).

---

## 1. Final Lock Summary

### Function signature + body

```python
"""V2 Phase F Proper-1: ECM-Modulated FA Rate Response.

Per-FA HB#4 multipliers scale base FA rate constants (Y4-sister: NOT
force scalers). Effective-rate FA dynamics delegate to existing 6.3a
step_focal_adhesions_static (Y6 traction algebra). Manual aggregate of
per-FA results (Y16: NO second delegate call to prevent twice-mutation
corruption).

Phase F minimal Y4 guard preserved: this is a NEW pathway, not a patch
to step_phase_f_minimal_motility. HB#4 multipliers used as per-rate
scalers ONLY.
"""

import numpy as np
from dataclasses import dataclass, replace
from typing import Optional

from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.dynamics.focal_adhesion import (
    FocalAdhesionDynamicsDiagnostics,
    FocalAdhesionDynamicsError,
    FocalAdhesionDynamicsParameters,
    FocalAdhesionDynamicsResult,
    _DT_RATE_SAFETY_MARGIN,
    step_focal_adhesions_static,
)
from acs.v2.dynamics.ecm_to_fa_bias import (
    ECMToFABiasResult,
    RATE_NAMES,
)


@dataclass(frozen=True, slots=True)
class FARateResponseDiagnostics:
    """Typed diagnostics for one F-proper-1 step."""
    n_adhesions: int
    effective_k_maturity_per_s: np.ndarray      # (n_fa,) 0.0 for None base
    effective_k_bind_per_s: np.ndarray
    effective_k_unbind_per_s: np.ndarray
    max_effective_rate_per_s: float              # 0.0 if no FAs or all None
    max_dt_fa_rate_product: float                # = dt_fa_s * max_effective_rate
    aggregate_traction_change_nN_xy: tuple[float, float]
    n_none_base_rates: int                       # count of None entries (Y17)


@dataclass(frozen=True, slots=True)
class FARateResponseResult:
    updated_adhesions: tuple[FocalAdhesionState, ...]
    fa_dynamics_result: FocalAdhesionDynamicsResult   # manual aggregate (Y16)
    diagnostics: FARateResponseDiagnostics


def step_fa_rate_response(
    adhesions: tuple[FocalAdhesionState, ...],
    centroid_um_xy: tuple[float, float],
    multipliers: ECMToFABiasResult,
    base_params: FocalAdhesionDynamicsParameters,
    *,
    traction_axis_xy_per_fa: Optional[list] = None,
) -> FARateResponseResult:
    """F-proper-1 step.

    Pipeline (Y20 base_params validation BEFORE empty-adhesion return):
      0. Validate base_params (NEVER bypass — Y20)
      1. Validate multipliers (Y14): fa_ids exact match, rate_names
         match, shape, finite/non-bool/non-negative
      2. Validate traction_axis_xy_per_fa length
      3. Build per-FA effective rates (Y17 None handling)
      4. dt-rate gate via base_params.dt_fa_s · max_effective_rate
         (≤ _DT_RATE_SAFETY_MARGIN; Y6 sister)
      5. Per-FA delegate to step_focal_adhesions_static; manual
         aggregate (Y16 — NO second delegate call)
      6. Diagnostics + return

    Y4 minimal guard preserved: this is NEW pathway. HB#4 multipliers
    used as per-rate scalers ONLY (per RATE_NAMES original design
    intent). Force change through 6.3a traction algebra only.
    """
    # 0. Y20 base_params validation FIRST (NEVER bypass via empty adhesions)
    _validate_base_params_for_proper(base_params)

    # 1. Multiplier validation (Y14 fail-closed)
    _validate_multipliers_match_adhesions(multipliers, adhesions)
    n_fa = len(adhesions)

    # 2. traction_axis length validation
    if traction_axis_xy_per_fa is not None and len(traction_axis_xy_per_fa) != n_fa:
        raise FocalAdhesionDynamicsError(
            "traction_axis_length_mismatch",
            f"traction_axis_xy_per_fa length {len(traction_axis_xy_per_fa)} != n_fa {n_fa}",
        )

    # 3. Per-FA effective rates (Y17 None handling)
    n_none = 0
    eff_k_maturity_list = []
    eff_k_bind_list = []
    eff_k_unbind_list = []
    for i in range(n_fa):
        m = multipliers.multipliers_per_fa[i]
        eff_k_maturity_list.append(_scale_rate_or_none(base_params.k_maturity_per_s, float(m[0])))
        eff_k_bind_list.append(_scale_rate_or_none(base_params.k_bind_per_s, float(m[1])))
        eff_k_unbind_list.append(_scale_rate_or_none(base_params.k_unbind_per_s, float(m[2])))
        if base_params.k_maturity_per_s is None:
            n_none += 1
        if base_params.k_bind_per_s is None:
            n_none += 1
        if base_params.k_unbind_per_s is None:
            n_none += 1

    eff_k_maturity_arr = np.array(
        [r if r is not None else 0.0 for r in eff_k_maturity_list], dtype=np.float64,
    )
    eff_k_bind_arr = np.array(
        [r if r is not None else 0.0 for r in eff_k_bind_list], dtype=np.float64,
    )
    eff_k_unbind_arr = np.array(
        [r if r is not None else 0.0 for r in eff_k_unbind_list], dtype=np.float64,
    )

    # 4. dt-rate gate (empty FAs → max_rate=0.0; minor refinement)
    if n_fa == 0:
        max_rate = 0.0
    else:
        max_rate = float(max(
            float(eff_k_maturity_arr.max()) if eff_k_maturity_arr.size else 0.0,
            float(eff_k_bind_arr.max()) if eff_k_bind_arr.size else 0.0,
            float(eff_k_unbind_arr.max()) if eff_k_unbind_arr.size else 0.0,
        ))
    dt_fa_s = base_params.dt_fa_s
    dt_rate_product = dt_fa_s * max_rate
    if dt_rate_product > _DT_RATE_SAFETY_MARGIN:
        raise FocalAdhesionDynamicsError(
            "dt_rate_violation",
            f"dt_fa_s={dt_fa_s} * max_effective_rate={max_rate:.4f} = "
            f"{dt_rate_product:.4f} > safety margin {_DT_RATE_SAFETY_MARGIN}",
        )

    # 5. Per-FA delegate + manual aggregate (Y16 — NO second delegate call)
    pre_traction = (
        np.array([fa.traction_force_nN_xy for fa in adhesions], dtype=np.float64)
        if n_fa > 0 else np.zeros((0, 2))
    )
    pre_traction_sum = pre_traction.sum(axis=0) if n_fa > 0 else np.zeros(2)

    cell_force = np.zeros((n_fa, 2), dtype=np.float64)
    substrate_reaction = np.zeros((n_fa, 2), dtype=np.float64)
    radial_components = np.zeros(n_fa, dtype=np.float64)
    tangential_components = np.zeros(n_fa, dtype=np.float64)
    updated = []

    for i, fa in enumerate(adhesions):
        per_fa_params = replace(
            base_params,
            k_maturity_per_s=eff_k_maturity_list[i],
            k_bind_per_s=eff_k_bind_list[i],
            k_unbind_per_s=eff_k_unbind_list[i],
        )
        axis = traction_axis_xy_per_fa[i] if traction_axis_xy_per_fa is not None else None
        result_one = step_focal_adhesions_static(
            (fa,), centroid_um_xy, per_fa_params,
            traction_axis_xy_per_fa=[axis],
        )
        updated.append(result_one.updated_adhesions[0])
        cell_force[i] = result_one.cell_force_nN_xy[0]
        substrate_reaction[i] = result_one.substrate_reaction_nN_xy[0]
        radial_components[i] = result_one.radial_components[0]
        tangential_components[i] = result_one.tangential_components[0]

    updated_tuple = tuple(updated)

    # Manual aggregate FocalAdhesionDynamicsResult (Y16 — NO second delegate)
    fa_dynamics_result = FocalAdhesionDynamicsResult(
        updated_adhesions=updated_tuple,
        cell_force_nN_xy=cell_force,
        substrate_reaction_nN_xy=substrate_reaction,
        radial_components=radial_components,
        tangential_components=tangential_components,
        diagnostics=FocalAdhesionDynamicsDiagnostics(
            n_adhesions=n_fa,
            aggregate_cell_force_nN_xy=(
                tuple(cell_force.sum(axis=0).tolist()) if n_fa > 0 else (0.0, 0.0)
            ),
            aggregate_substrate_reaction_nN_xy=(
                tuple(substrate_reaction.sum(axis=0).tolist())
                if n_fa > 0 else (0.0, 0.0)
            ),
            max_traction_magnitude_nN=(
                float(np.linalg.norm(cell_force, axis=1).max()) if n_fa > 0 else 0.0
            ),
        ),
    )

    # 6. Diagnostics
    post_traction_sum = cell_force.sum(axis=0) if n_fa > 0 else np.zeros(2)
    aggregate_change = tuple((post_traction_sum - pre_traction_sum).tolist())

    diagnostics = FARateResponseDiagnostics(
        n_adhesions=n_fa,
        effective_k_maturity_per_s=eff_k_maturity_arr,
        effective_k_bind_per_s=eff_k_bind_arr,
        effective_k_unbind_per_s=eff_k_unbind_arr,
        max_effective_rate_per_s=max_rate,
        max_dt_fa_rate_product=dt_rate_product,
        aggregate_traction_change_nN_xy=aggregate_change,
        n_none_base_rates=n_none,
    )

    return FARateResponseResult(
        updated_adhesions=updated_tuple,
        fa_dynamics_result=fa_dynamics_result,
        diagnostics=diagnostics,
    )


def _validate_base_params_for_proper(base_params: FocalAdhesionDynamicsParameters) -> None:
    """Y20 fail-closed validation BEFORE empty-adhesion return path.

    Mirrors 6.3b two-pass pattern: validate scalar/rate fields upfront
    so empty case cannot bypass schema validation. Uses base_params'
    own validate() if compatible, else manual scalar checks.
    """
    # Use existing validate() pattern from 6.3a/6.3b sister precedent.
    # If base rates × dt_fa_s exceed margin (which would fail 6.3a gate
    # for the base case but be safe under multiplier <1), use 6.3b
    # two-pass: validate dt/traction/max-traction separately + raw rates
    # with dt=0, then F-proper-1 effective-rate gate handles dt validation.
    base_params.validate()  # default; 6.3b two-pass if needed (impl decision)


def _scale_rate_or_none(base_rate: Optional[float], multiplier: float) -> Optional[float]:
    """Y17: None base rate stays None (= "no update" in 6.3a contract);
    float scales by multiplier."""
    if base_rate is None:
        return None
    return float(base_rate) * float(multiplier)


def _validate_multipliers_match_adhesions(
    multipliers: ECMToFABiasResult,
    adhesions: tuple[FocalAdhesionState, ...],
) -> None:
    """Y14 fail-closed validation. No silent reorder in F-proper-1."""
    expected_ids = tuple(fa.adhesion_id for fa in adhesions)
    if multipliers.fa_ids != expected_ids:
        raise FocalAdhesionDynamicsError(
            "multipliers_fa_ids_mismatch",
            f"multipliers.fa_ids {multipliers.fa_ids} != adhesions ids "
            f"{expected_ids}. Caller must align order; silent reorder "
            f"forbidden in F-proper-1 first lock.",
        )
    if multipliers.rate_names != RATE_NAMES:
        raise FocalAdhesionDynamicsError(
            "multipliers_rate_names_mismatch",
            f"multipliers.rate_names {multipliers.rate_names} != RATE_NAMES "
            f"{RATE_NAMES}. Strict ordering required.",
        )
    n_fa = len(adhesions)
    expected_shape = (n_fa, len(RATE_NAMES))
    arr = np.asarray(multipliers.multipliers_per_fa)
    if arr.shape != expected_shape:
        raise FocalAdhesionDynamicsError(
            "multipliers_shape_mismatch",
            f"multipliers.multipliers_per_fa.shape {arr.shape} != expected "
            f"{expected_shape}",
        )
    # Reject bool dtype + object arrays containing bool (minor refinement)
    if arr.dtype == bool:
        raise FocalAdhesionDynamicsError(
            "multipliers_bool_dtype",
            "multipliers must not be bool dtype (Python bool ⊂ int trap)",
        )
    arr_f = np.asarray(arr, dtype=np.float64)
    if not np.all(np.isfinite(arr_f)):
        raise FocalAdhesionDynamicsError(
            "multipliers_non_finite",
            f"multipliers must be finite; got {arr}",
        )
    if (arr_f < 0).any():
        raise FocalAdhesionDynamicsError(
            "multipliers_negative",
            f"multipliers must be non-negative; got {arr}",
        )
```

### Forbidden in F-proper-1 (Y3-Y20 enforce)

- ❌ **Second `step_focal_adhesions_static` call after per-FA effective-rate steps**
  (Y16 — silent twice-mutation/corruption; manual aggregate only)
- ❌ Force scaling via multipliers (Y4 sister; multipliers RATE-only)
- ❌ FA spawning / nucleation / RNG / new state labels (Y5 + Y15)
- ❌ New constitutive law / new HB lock (Y5 + Y6 — delegate to 6.3a)
- ❌ Contour motion / vertex displacement (F-proper-2 owns)
- ❌ Subcycling beyond FA rate gate (F-proper-2 owns)
- ❌ "boundary-normal outward" axis wording (Y10 sign trap)
- ❌ Silent multiplier reorder (Y14)
- ❌ PI 260313 fitting (Hard Rule 1)
- ❌ Magic-number rate constants (Y8 — Cho 2020 + epithelial literature)
- ❌ `None * multiplier` (Y17 — explicit None preservation)
- ❌ Empty-adhesion return path bypassing `base_params` validation
  (Y20 — schema-validation bypass forbidden)
- ❌ Patch to `step_phase_f_minimal_motility` (Y4 minimal guard preserved)

---

## 2. Reasoned-acceptance trace (Y1–Y20)

The lock converged after 4 rounds (id 2127–2168) under PI directive
id=1985 aggressive debate posture, with audit `id=2114`
plumbing-vs-mechanism gap framing absorbed.

**Y1 (Codex round 1: 5-area scope decomposition)**: timescale /
ECM→cell feedback / FA dynamics / visible pilot / acceptance metrics.

**Y2 (Codex round 1: (D) staged A→B feedback path)**: Stage 1 = FA
rate feedback (this lock); Stage 2 = traction force generation
(separate lock); avoids direct force-scaling magic-number risk.

**Y3 (Codex round 2: 3-stage staging F-proper-1/2/3 separate locks)**:
F-proper-1 alone does NOT guarantee visible motion; need integrator
(F-proper-2) + long-horizon pilot (F-proper-3).

**Y4 (Codex round 2: Phase F minimal Y4 guard preserved forever; new
pathway in proper)**: `step_phase_f_minimal_motility` Y4 (HB#4
multipliers diagnostic-only) stays. F-proper-1 is a NEW module/path,
not a patch.

**Y5 (Codex round 2: existing HB#4 + 6.3a/6.3b, no new HB, no new rate
names, no nucleation)**: Use `RATE_NAMES = ("k_maturity_per_s",
"k_bind_per_s", "k_unbind_per_s")`. No `k_nucleation_per_s` (RNG/site
selection separate lock). Fixed candidate FA set; deterministic.

**Y6 (Codex round 2: delegate to 6.3a algebra)**: Force change comes
through `traction = traction_scale_nN · maturity · bound_fraction ·
unit_axis` (locked 6.3a). No direct multiplier-to-force scaling.

**Y7 (Codex round 2 Q10 → round 3 sign correction)**: 6.3a default
inward radial-to-centroid traction CONTRACTS cell. For motility,
**boundary-normal INWARD cell-on-substrate traction** (substrate-on-cell
OUTWARD reaction → spreading). Round 2 had it reversed (boundary-normal
OUTWARD = contraction); Codex round 3 corrected.

**Y8 (Codex round 2 Q11: multi-rate budget)**: dt_fa_s, dt_ecm_s,
dt_cell_s as pilot config (not hardcoded defaults), validated, reported.
F-proper-1 owns FA rate gate only; F-proper-2 owns contour subcycling.

**Y9 (Codex round 2: F-proper-1 lock target)**: `v2 Phase F Proper-1:
ECM-Modulated FA Rate Response — Locked Design`.

**Y10 (Codex round 3 Q10 CRITICAL sign correction)**: Round 2 (iv)
"boundary-normal outward" was WRONG (= cell-on-substrate outward =
contraction). Correct: cell-on-substrate INWARD via
`traction_axis = -outward_normal_at_vertex`. Substrate-on-cell reaction
= OUTWARD. Same Y2 sign error pattern Phase F minimal had.

**Y11 (Codex round 3 Q11: subcycling pilot config + reporting)**:
Caller-supplied with validation; diagnostics report contour substep
count, max dt-rate products.

**Y12 (Codex round 3 Q12: linear gradient ECM IC + uniform baseline)**:
Single fixture = overinterpretation risk. Always include uniform ECM
baseline (negative control: symmetric/no-net-migration) alongside
gradient (positive evidence).

**Y13 (Codex round 3 API: dual dt source removed)**: Use
`base_params.dt_fa_s` only. No separate `dt_s` argument.

**Y14 (Codex round 3 Q14: multiplier ordering + FA identity validation)**:
fa_ids exact match, rate_names exact match, shape, finite, non-bool,
non-negative. No silent reorder.

**Y15 (Codex round 3 Q15: state transitions stay 6.3a semantics)**:
F-proper-1 must NOT invent FA state-label transitions. 6.3a updates
maturity/bound_fraction only. PI-visible state transitions = separate
lock.

**Y16 (Codex round 4 P0 CRITICAL: NO second `step_focal_adhesions_static`
call)**: Round 3 candidate had a fatal bug — re-stepping `updated_tuple`
with base_params after per-FA effective-rate steps caused twice
mutation/corruption. Manual aggregation of per-FA results required.

**Y17 (Codex round 4 P1: explicit None base rate handling)**:
`None * multiplier` = TypeError. Per-FA params preserve None (= "no
update" in 6.3a contract); diagnostics arrays use 0.0 for None entries.

**Y18 (Codex round 4 P2: test semantics correction)**: Neutral
multiplier (1.0) = exactly 6.3a base behavior (NOT no traction change).
Test 8 corrected: neutral matches manual `step_focal_adhesions_static`
with same base params. Test 9: multiplier >1 on k_maturity/k_bind →
larger traction for unsaturated FA. Test 10: multiplier >1 on k_unbind
→ lower traction for slipping FA. "All <1 = traction decreases" wording
forbidden (per-rate effects distinct).

**Y19 (Codex round 4 P3: F-proper-1 sign convention test at vector
level)**: Round 3 test 14 "boundary-normal INWARD axis → outward vertex
reaction" wrong scope (vertex motion = F-proper-2). Corrected:
F-proper-1 tests sign at vector level only:
`traction_axis = -outward_normal` → `cell_force · outward_normal < 0`
+ `substrate_reaction · outward_normal > 0`.

**Y20 (Codex round 4 SEAL guard: validate base_params before empty-
adhesion return path)**: `adhesions=()` cannot bypass schema validation
of `base_params`. Mirror 6.3b two-pass pattern: validate scalar/rate
fields upfront. Also adds `traction_axis_xy_per_fa` length mismatch
test.

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

1. **Units**:
   - `multipliers.multipliers_per_fa`: dimensionless ∈ [0, ∞)
   - `base_params.k_*_per_s`: 1/s
   - `effective_rate = base_rate × multiplier`: 1/s
   - `dt_fa_s × max_effective_rate`: dimensionless (P1-style stability)
   - 6.3a delegated traction: nN
2. **Boundary**:
   - Empty adhesions: max_rate=0.0; aggregate fields=(0,0); diagnostics
     populated; `base_params.validate()` still runs (Y20)
   - All None base rates: per-FA params preserve None; diagnostics
     arrays = 0.0; max_rate = 0.0; n_none_base_rates reflects count
   - Multiplier=1.0 everywhere: matches 6.3a base behavior with same
     params (Y18)
   - dt_rate_product > safety margin: raise dt_rate_violation
   - traction_axis_xy_per_fa length mismatch: raise before mutation
3. **Conservation**:
   - 6.3a Newton-3 conservation delegated to per-FA results
     (cell_force = -substrate_reaction per FA)
   - Manual aggregate sums preserve per-FA Newton-3 → aggregate
     cell_force_sum = -aggregate substrate_reaction_sum (modulo
     floating-point)
4. **Numerical**:
   - `np.float64` enforced (multiplier asarray cast, effective rate
     arrays)
   - dt-rate gate checked BEFORE delegate call (fail-fast)
   - `_DT_RATE_SAFETY_MARGIN` reused from 6.3a (no new constant)
   - Empty-array `.max()` avoided (explicit conditional)
5. **Sign**:
   - Multipliers ≥ 0 (validated)
   - Effective rates = base × multiplier ≥ 0 (since base ≥ 0 if not
     None, multiplier ≥ 0)
   - 6.3a Newton-3 sign preserved per per-FA delegate
   - Y10 axis convention: `traction_axis = -outward_normal` →
     `cell_force · outward_normal < 0` (cell-on-substrate INWARD) →
     `substrate_reaction · outward_normal > 0` (substrate-on-cell
     OUTWARD = spreading)
6. **Measurement-protocol** (Hard Rule 11, central anchor):
   - **Y4 minimal guard preserved**: this is NEW pathway, not a patch.
     `step_phase_f_minimal_motility`'s multipliers stay diagnostic-only.
   - **Y6 traction algebra delegation**: force change comes through
     6.3a algebra (`traction = traction_scale_nN · maturity ·
     bound_fraction · unit_axis`); NOT through direct multiplier-to-force
     scaling.
   - **Y10 sign convention**: cell-on-substrate INWARD =
     substrate-on-cell OUTWARD = spreading. `traction_axis =
     -outward_normal_at_vertex` for boundary-normal INWARD form.
   - **Y14 ordering strictness**: no silent reorder; caller must align
     `multipliers.fa_ids` with `tuple(fa.adhesion_id for fa in adhesions)`
     exactly.
   - **Y15 state semantics**: maturity/bound_fraction updated per 6.3a;
     no new state-label transitions invented.
   - **Y16 manual aggregation**: per-FA delegate results aggregated
     manually (not via second delegate call) to preserve effective-rate
     step outputs.
   - **Y17 None semantics**: `None` base rate = "no update" per 6.3a
     contract; preserved in per-FA params; 0.0 in diagnostics arrays.
   - **Y20 validation discipline**: `base_params.validate()` runs
     before empty-adhesion return path; schema validation never bypassed.
   - **Audit `id=2114` framing**: F-proper-1 alone does NOT produce
     visible cell motion; F-proper-2 (integrator) + F-proper-3 (15-min
     annotated pilot) are required for PI-readable mp4.

---

## 4. Test Catalog (~17 tests, A-tier full)

### Validation (7)

1. `test_fa_rate_response_multipliers_shape_mismatch_raises` (Y14)
2. `test_fa_rate_response_multipliers_fa_ids_mismatch_raises` (Y14)
3. `test_fa_rate_response_multipliers_rate_names_mismatch_raises` (Y14)
4. `test_fa_rate_response_multipliers_negative_raises` (Y14)
5. `test_fa_rate_response_multipliers_bool_dtype_raises` (Y14)
6. `test_fa_rate_response_dt_rate_violation_raises` (Y6 stability)
7. `test_fa_rate_response_traction_axis_length_mismatch_raises`
   (Y20 minor refinement)

### Boundary + invalid base_params (2)

8. `test_fa_rate_response_empty_adhesions_validates_base_params`
   (Y20 — empty case still runs `base_params.validate()`; schema
   bypass forbidden) — uses invalid base_params to verify validation
   raises even with `adhesions=()`
9. `test_fa_rate_response_empty_adhesions_max_rate_zero` (minor —
   `max_effective_rate_per_s == 0.0`, `aggregate_change == (0, 0)`)

### ODE behavior (3, Y18 corrected)

10. `test_fa_rate_response_neutral_multipliers_match_manual_6_3a`
    (Y18 corrected: neutral 1.0 = exactly 6.3a base behavior; result
    matches manual `step_focal_adhesions_static` with same base
    params)
11. `test_fa_rate_response_high_maturity_or_bind_multiplier_increases_traction_for_unsaturated_fa`
    (Y18 corrected: multiplier > 1 on `k_maturity` or `k_bind` →
    larger post-step traction for unsaturated nascent/bound<1 FA)
12. `test_fa_rate_response_high_unbind_multiplier_decreases_traction_for_slipping_fa`
    (Y18 corrected: multiplier > 1 on `k_unbind` → lower post-step
    traction for slipping FA — NOT for unsaturated baseline)

### Per-FA differentiation + None handling (2)

13. `test_fa_rate_response_per_fa_distinct_multipliers_produce_distinct_effective_rates`
    (multipliers vary by FA; effective_*_arr reflects per-FA distinction)
14. `test_fa_rate_response_none_base_rate_preserved_in_per_fa_params_zero_in_diagnostics`
    (Y17 — `base_params.k_maturity_per_s = None` → per-FA params have
    `None`; diagnostics array `effective_k_maturity_per_s[i] = 0.0`)

### Diagnostics + sign convention (2)

15. `test_fa_rate_response_diagnostics_aggregate_traction_change_matches_pre_post_diff`
    (manual cross-check)
16. `test_fa_rate_response_boundary_normal_inward_traction_axis_produces_outward_substrate_reaction_vector`
    (Y19 — F-proper-1 vector level; NO vertex motion test):
    `traction_axis = -outward_normal` →
    `cell_force · outward_normal < 0` AND
    `substrate_reaction · outward_normal > 0`

### Y16 P0 guard (1) — critical regression test

17. `test_fa_rate_response_no_double_step_corruption`
    (Y16 P0 — explicit regression: verify `result.fa_dynamics_result`
    matches per-FA effective-rate step outputs aggregated, NOT a
    second `step_focal_adhesions_static` call with `base_params`.
    Inject distinct multipliers per FA; verify
    `result.updated_adhesions[i]` matches per-FA effective-rate result,
    not base-rate result.)

### Exports (bundled)

Tests 1-17 collectively exercise imports of
`FARateResponseDiagnostics`, `FARateResponseResult`,
`step_fa_rate_response` through `acs.v2.dynamics.fa_rate_response`
+ `acs.v2.dynamics` re-exports + `acs.v2` mirror.

---

## 5. Files

- `docs/v2/v2_phase_f_proper_1_fa_rate_response_locked.md` (this file)
- `acs/v2/dynamics/fa_rate_response.py` (NEW):
  - `FARateResponseDiagnostics` dataclass (frozen+slots, 8 fields
    incl. `n_none_base_rates`)
  - `FARateResponseResult` dataclass (frozen+slots, 3 fields)
  - `step_fa_rate_response` public function
  - Private helpers: `_validate_base_params_for_proper`,
    `_scale_rate_or_none`, `_validate_multipliers_match_adhesions`
- `acs/v2/dynamics/__init__.py` + `acs/v2/__init__.py` — 3 new exports
- `tests/v2/test_v2_fa_rate_response.py` (NEW, ~17 tests per §4)

---

## 6. References

- PI directives:
  - `id=2099`: 시각자료 부족 frustration (audit trigger)
  - `id=2114`: external audit (plumbing-vs-mechanism gap; F-proper-3
    minimum conditions)
  - `id=2123`: 바로 시작 (immediate (X+Y) parallel routing)
  - `id=1985`: 건강하고 활발한 토론 (aggressive debate posture)
  - `id=809`: aggressive design debate posture (memory)
- 3-stage Phase F proper sequence:
  - **F-proper-1** (this lock): ECM-modulated FA rate dynamics
  - F-proper-2 (separate lock): multi-rate motility integrator
  - F-proper-3 (separate lock): PI-readable 15-min annotated pilot
- Sister locks (preserved, NOT patched):
  - `docs/v2/v2_phase_f_minimal_cell_motility_pilot_locked.md` (Y4 minimal
    guard preserved forever)
  - `docs/v2/v2_phase_e_v2_composition_step_2_locked.md`
  - `docs/v2/v2_hard_blocker_4_active_locked.md` (HB#4-active multipliers)
  - `docs/v2/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1)
  - `docs/v2/v2_63b_protrusion_coupling_locked.md` (6.3b sister-pattern
    for two-pass validation per Y20)
- Verified at source for delegation:
  - `acs/v2/dynamics/focal_adhesion.py` (`step_focal_adhesions_static`,
    `FocalAdhesionDynamicsResult`, `FocalAdhesionDynamicsDiagnostics`,
    `FocalAdhesionDynamicsError`, `_DT_RATE_SAFETY_MARGIN`)
  - `acs/v2/dynamics/ecm_to_fa_bias.py` (`ECMToFABiasResult`,
    `RATE_NAMES`)
- Audit-driven framing source:
  - `runs/20260505T135202Z_phase_f_minimal_motility_pilot/metadata.json`
    (`"not_mechanistic": true`, `"ecm_to_cell_feedback": false`)
- Hard Rule 1 (no PI data fitting): CLAUDE.md
- Hard Rule 11 (sim-experiment measurement matching): CLAUDE.md
- Magic-Number Block: CLAUDE.md (Y8 — rate constants literature
  anchor required at impl time)

---

## 7. Cross-room dispatch

This file is the design-team input to implementation-work for:

1. impl Claude/Codex implements
   `acs/v2/dynamics/fa_rate_response.py` + tests per §1 pseudocode +
   §4 17-test catalog. **Skip Sanity Gate doc** (A-tier full debate
   already done in design; impl direct + Codex single review = milestone).
2. impl-side single review pass on commit covering:
   - Y16 P0 critical: NO second `step_focal_adhesions_static` call;
     manual aggregation of per-FA results
   - Y17 None handling: per-FA params preserve `None`; diagnostics
     arrays use 0.0
   - Y18 test semantics: neutral = base behavior, NOT no-change;
     per-rate effects distinct
   - Y19 sign convention test at vector level only (NO vertex motion
     test in F-proper-1)
   - Y20 base_params validation BEFORE empty-adhesion return path
   - Y14 multiplier ordering strict (no silent reorder)
   - Y4 Phase F minimal guard preserved (NEW pathway, not patch)
3. On PASS: F-proper-1 milestone reached.
4. After milestone: design-discussion can open **F-proper-2 (multi-rate
   motility integrator)** as next A-tier full debate cycle. F-proper-2
   uses F-proper-1 traction output to drive contour vertex displacement
   via Newton 3rd reaction (substrate-on-cell), with subcycling between
   FA dynamics and contour dynamics. F-proper-3 (15-min annotated
   pilot) follows F-proper-2 milestone.

A-tier debate cycle MCP id 2127–2168 (rounds 1-4 + Codex SEAL with Y20
guard).
