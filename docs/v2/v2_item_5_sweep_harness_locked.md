# V2 Closed-Loop ECM Gate — Item 5 Sweep Harness (Phase E v1 Sensitivity Evidence) — Locked Design

**Date**: 2026-05-05 KST
**Authors**: Claude + Codex design-discussion (4-round adversarial lock,
PI id=809 aggressive debate posture applied, PI id=1008/1057/1161
autonomy + visible deliverable focus)
**Source unit**: design-discussion `topic=v2-item-5-sweep-harness`,
MCP id 1564–1575 (rounds 1–4 + seal ack)
**Brief**: `docs/v2/v2_item_5_sweep_harness_brief.md` (commit `4f0efdd`)
**Parent locked plan**: `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md`
**Upstream sister locks**:
- `docs/v2/v2_phase_e_composition_locked.md` (Phase E v1, the harness target)
- `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md` (HB#5, V_active definition)
- `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md` (HB#1+#2, K_ORIENT/TRACTION_REF constants)
- `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md` (HB#3, scatter contract)
- `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md` (HB#4 neutral readout)
- `docs/v2/v2_phase_d_no_op_scaffolding_locked.md`
**PI ratify status**: full delegation per PI id=939/1008. impl-work uses
this for the Item 5 Sanity Gate doc + code entry.

---

## 0. Scope

This unit locks the **Item 5 sensitivity sweep harness** for Phase E v1.
It runs Phase E v1 across a (`spacing_um`, `dt_s`) cross-product
under a **fixed physical domain + invariant FA scenario**, and
returns raw observable summaries.

### What this harness IS

- A pure deterministic library function `run_phase_e_v1_sensitivity_sweep`
  that wraps `step_closed_loop_phase_e_v1` over a sweep
- Input: a single typed `PhaseEV1Item5SweepConfig` dataclass
- Output: a typed `PhaseEV1Item5SweepResult` with full `(n_runs,
  n_steps+1)` trajectories + raw aggregate summaries + typed
  reproducibility metadata

### What this harness is NOT

- **NOT a satisfaction or pass/fail decision for Item 5**: the harness
  produces raw evidence; whether observed sensitivity satisfies any
  specific Item 5 gate criterion is a **separate caller-side decision**
  (PI / explicit-decision unit). Tests are named
  `provides_*_evidence`, NOT `satisfies_*` (Y2 sister-pattern with
  Phase E v1 lock).
- **NOT a Phase E v2 (HB#4-active) harness**: HB#4 stays neutral via
  Phase E v1; the test catalog enforces `max_neutral_multiplier_deviation
  == 0.0` everywhere.
- **NOT a fitting/tuning utility**: per Codex `id=1563` review focus,
  no PI-data fitting; no gate tolerance edits; literature-derived
  parameters only.

### Wording discipline (Y4 — central anchor)

`dt_s` axis is the **ECM-feedback integration substep**, NOT the
experimental imaging interval (which is 15 min production / 60 min pilot
per v1 docs). The chosen range `{30, 60, 120, 300}` s is selected to
span sub-minute to sub-imaging-interval substep sensitivity. This is
the second Hard Rule 11 measurement-protocol catch in Phase E composition
(Phase E v1 Y1 was the first; Item 5 Y1 here is the third — physical
domain invariance under spacing sweep).

---

## 1. Final Lock Summary

### Typed config + result + metadata

```python
import sys
from dataclasses import dataclass
from typing import Literal

import numpy as np

from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.dynamics.closed_loop_phase_e import (
    PhaseEStepResult,
    step_closed_loop_phase_e_v1,
)
from acs.v2.dynamics.ecm_constitutive_response import (
    K_ORIENT_PER_S,
    TRACTION_REF_NN_PER_UM2,
)


@dataclass(frozen=True, slots=True)
class PhaseEV1Item5SweepConfig:
    """Typed config for one Item 5 sweep harness invocation.

    Y1 contract: physical domain + FA position invariant across spacing
    sweep; only resolution changes. Domain must be exactly divisible by
    every sweep spacing (validation raises otherwise).
    """
    spacing_um_values: tuple[float, ...]                   # NON-empty, finite positive
    dt_s_values: tuple[float, ...]                          # NON-empty, finite positive
    n_steps: int                                             # positive
    domain_size_um_xy: tuple[float, float]                  # finite positive
    origin_um_xy: tuple[float, float]                       # finite
    traction_scenario: Literal["rotating_uniform_single_fa"]
    traction_magnitude_nN: float                            # finite positive
    fa_position_um_xy: tuple[float, float]                  # finite, inside domain
    n_revolutions: int = 1                                  # positive int (no bool)
    git_sha_override: str | None = None                     # tests pass fixed value


@dataclass(frozen=True, slots=True)
class PhaseEV1Item5SweepMetadata:
    """Typed run metadata for reproducibility (Y13 — NO wall-clock timestamp)."""
    config: PhaseEV1Item5SweepConfig
    git_sha: str                                            # detected or override (Y14)
    numpy_version: str
    python_version: str
    k_orient_per_s: float                                   # imported HB#1+#2 constant
    traction_ref_nN_per_um2: float                          # imported HB#1+#2 constant
    # NO timestamp_iso — caller can wrap with own timestamp if file-write concerns


@dataclass(frozen=True, slots=True)
class PhaseEV1Item5SweepResult:
    """Typed sweep result. Trajectory shapes (n_runs, n_steps+1) include step 0 initial."""
    # Per-run metadata (length n_runs = len(spacing_um_values) * len(dt_s_values))
    spacing_um_by_run: np.ndarray                          # (n_runs,) float64
    dt_s_by_run: np.ndarray                                # (n_runs,) float64
    nx_by_run: np.ndarray                                  # (n_runs,) int64
    ny_by_run: np.ndarray                                  # (n_runs,) int64
    # Per-run per-step trajectories (n_steps+1 includes step 0)
    v_active_um2_traj: np.ndarray                          # (n_runs, n_steps+1)
    v_active_bound_um2_traj: np.ndarray                    # (n_runs, n_steps+1)
    max_convex_weight_traj: np.ndarray                     # (n_runs, n_steps+1)
    max_orientation_delta_frobenius_traj: np.ndarray       # (n_runs, n_steps+1)
    max_traction_norm_nN_per_um2_traj: np.ndarray          # (n_runs, n_steps+1)
    neutral_multiplier_max_deviation_traj: np.ndarray      # (n_runs, n_steps+1)
    # Per-run aggregates (raw, NOT pass/fail)
    max_bound_ratio: np.ndarray                            # (n_runs,)
    max_neutral_multiplier_deviation: np.ndarray           # (n_runs,)
    peak_v_active_um2: np.ndarray                          # (n_runs,)
    final_v_active_um2: np.ndarray                         # (n_runs,)
    metadata: PhaseEV1Item5SweepMetadata
```

### Harness function

```python
def run_phase_e_v1_sensitivity_sweep(
    config: PhaseEV1Item5SweepConfig,
) -> PhaseEV1Item5SweepResult:
    """Run cross-product sweep over (spacing_um × dt_s).

    For each (spacing_um, dt_s) pair:
      1. Build initial ECM at this spacing via deterministic
         physical-coordinate IC function (Y1)
      2. Capture step 0 V_active by calling Phase E v1 with dt_s=0
         (initial state measured against step-0 traction target — Y17)
      3. Loop n_steps with rotating traction direction, propagating
         updated_ecm forward
      4. Record per-step trajectories + aggregate summaries
    """
    _validate_config(config)  # Y3 + Y20 + Y21 + Y1 divisibility

    n_runs = len(config.spacing_um_values) * len(config.dt_s_values)
    n_steps_plus_1 = config.n_steps + 1

    # Allocate
    spacing_by = np.empty(n_runs, dtype=np.float64)
    dt_by = np.empty(n_runs, dtype=np.float64)
    nx_by = np.empty(n_runs, dtype=np.int64)
    ny_by = np.empty(n_runs, dtype=np.int64)
    v_active_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    v_bound_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    weight_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    delta_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    traction_norm_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    multiplier_dev_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)

    run_idx = 0
    for spacing_um in config.spacing_um_values:
        for dt_s in config.dt_s_values:
            ecm = _build_initial_ecm_at_spacing(config, spacing_um)
            spacing_by[run_idx] = spacing_um
            dt_by[run_idx] = dt_s
            nx_by[run_idx] = ecm.grid_shape[0]
            ny_by[run_idx] = ecm.grid_shape[1]

            # Step 0 = initial ECM measured against step-0 traction (Y17)
            init_adhesions = _build_rotating_uniform_single_fa(config, 0, dt_s)
            init_result = step_closed_loop_phase_e_v1(init_adhesions, ecm, dt_s=0.0)
            _record_step(
                init_result, run_idx, 0, v_active_traj, v_bound_traj,
                weight_traj, delta_traj, traction_norm_traj, multiplier_dev_traj,
            )
            # Note: step 0 max_convex_weight + max_orientation_delta == 0 exactly
            # (because dt_s=0 → w = -expm1(0) = 0); v_active may be nonzero.

            # Steps 1..n_steps
            for step in range(1, config.n_steps + 1):
                adhesions = _build_rotating_uniform_single_fa(config, step, dt_s)
                step_result = step_closed_loop_phase_e_v1(adhesions, ecm, dt_s=dt_s)
                _record_step(
                    step_result, run_idx, step, v_active_traj, v_bound_traj,
                    weight_traj, delta_traj, traction_norm_traj, multiplier_dev_traj,
                )
                ecm = step_result.updated_ecm  # propagate forward

            run_idx += 1

    # Aggregates (Y8, Y16 zero-bound case)
    max_bound_ratio = np.array(
        [_compute_max_bound_ratio(v_active_traj[i], v_bound_traj[i]) for i in range(n_runs)]
    )

    metadata = PhaseEV1Item5SweepMetadata(
        config=config,
        git_sha=_detect_git_sha(config.git_sha_override),
        numpy_version=np.__version__,
        python_version=sys.version,
        k_orient_per_s=K_ORIENT_PER_S,
        traction_ref_nN_per_um2=TRACTION_REF_NN_PER_UM2,
    )

    return PhaseEV1Item5SweepResult(
        spacing_um_by_run=spacing_by,
        dt_s_by_run=dt_by,
        nx_by_run=nx_by,
        ny_by_run=ny_by,
        v_active_um2_traj=v_active_traj,
        v_active_bound_um2_traj=v_bound_traj,
        max_convex_weight_traj=weight_traj,
        max_orientation_delta_frobenius_traj=delta_traj,
        max_traction_norm_nN_per_um2_traj=traction_norm_traj,
        neutral_multiplier_max_deviation_traj=multiplier_dev_traj,
        max_bound_ratio=max_bound_ratio,
        max_neutral_multiplier_deviation=multiplier_dev_traj.max(axis=1),
        peak_v_active_um2=v_active_traj.max(axis=1),
        final_v_active_um2=v_active_traj[:, -1],
        metadata=metadata,
    )
```

### Private builders (locked spec, NOT exported)

```python
def _build_rotating_uniform_single_fa(
    config: PhaseEV1Item5SweepConfig, step_index: int, dt_s: float,
) -> tuple[FocalAdhesionState, ...]:
    """Single-FA scenario with deterministic rotating traction direction (Y9, Y10).

    age_s = step_index * dt_s preserves physical-time metadata without
    affecting HB#3 scatter (scatter is age-independent). Returns (fa,) only.
    """
    theta = 2.0 * np.pi * config.n_revolutions * step_index / config.n_steps
    force_x = config.traction_magnitude_nN * np.cos(theta)
    force_y = config.traction_magnitude_nN * np.sin(theta)
    fa = FocalAdhesionState(
        adhesion_id="item5-fa-0",
        cell_id="item5-cell-0",
        position_um_xy=tuple(config.fa_position_um_xy),
        age_s=float(step_index) * dt_s,
        maturity=1.0,
        bound_fraction=1.0,
        state="mature",
        traction_force_nN_xy=(float(force_x), float(force_y)),
        linked_protrusion_id=None,
        source="simulated",
    )
    fa.validate()
    return (fa,)


def _build_initial_ecm_at_spacing(
    config: PhaseEV1Item5SweepConfig, spacing_um: float,
) -> ECMSubstrateState:
    """Sample physical-coordinate IC function onto specified grid (Y1, Y11).

    - Cell-centered grid: x_i = origin_x + (i + 0.5) * spacing_um
    - Stiffness: uniform 1.0 kPa across domain
    - Ligand density: uniform 0.5
    - Fiber density: uniform 0.3
    - Orientation: 0.5 * I (latent isotropic IC; Y11 Magic-Number Block-compliant)
    - Accumulated traction: zero
    """
    nx = int(round(config.domain_size_um_xy[0] / spacing_um))
    ny = int(round(config.domain_size_um_xy[1] / spacing_um))
    return ECMSubstrateState(
        origin_um_xy=tuple(config.origin_um_xy),
        spacing_um=float(spacing_um),
        stiffness_kpa=np.ones((nx, ny), dtype=np.float64),
        ligand_density=np.full((nx, ny), 0.5, dtype=np.float64),
        fiber_density=np.full((nx, ny), 0.3, dtype=np.float64),
        orientation_tensor=np.broadcast_to(
            0.5 * np.eye(2, dtype=np.float64), (nx, ny, 2, 2)
        ).copy(),
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
    )


def _detect_git_sha(override: str | None) -> str:
    """Best-effort git SHA detection (Y14).

    Tests pass override for determinism. Caller-supplied None →
    best-effort subprocess invocation; any failure → "unknown".
    """
    if override is not None:
        return override
    try:
        import subprocess
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True, timeout=2.0,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return "unknown"


def _compute_max_bound_ratio(v_active_traj: np.ndarray, v_bound_traj: np.ndarray) -> float:
    """Max V_active / V_bound over active-bound-positive samples (Y16).

    Returns 0.0 if no active-bound samples (no division by zero).
    Mathematical invariant: ratio <= 1.0 by HB#5 bound (4.5 · active_area).
    """
    mask = v_bound_traj > 0.0
    if not mask.any():
        return 0.0
    return float((v_active_traj[mask] / v_bound_traj[mask]).max())
```

### Recommended config (locked grid table)

`domain_size_um_xy = (16.0, 16.0)` + `spacing_um ∈ {0.5, 1.0, 2.0, 4.0}`:

| spacing_um | nx | ny | n_cells | notes                              |
|-----------:|---:|---:|--------:|------------------------------------|
| 0.5        | 32 | 32 | 1024    | finest resolution case             |
| 1.0        | 16 | 16 | 256     | mid-resolution case                |
| 2.0        | 8  | 8  | 64      | coarse case                        |
| 4.0        | 4  | 4  | 16      | extrapolation limit case (intentional small grid for resolution-sensitivity) |

Total **4 × 4 = 16 runs**, design budget **≤ 10 min CPU on Laptop A5000
16 GB baseline** (CPU-only harness; VRAM N/A). Actual runtime
informational, NOT asserted in tests (Y15).

### Forbidden in Item 5 sweep harness

- "Item 5 satisfied" / "satisfies_item_5_*" wording in module
  docstring, function docstring, or test names (Y8 + sister-pattern
  Phase E v1 Y2 — meta-test enforces evidence-only naming)
- Reuse of ECM array across spacings without resampling on the new
  grid (Y1 — would conflate resolution with IC change)
- FA position varying across spacings (Y1 — physical position must be
  invariant)
- `domain_size_um_xy` not exactly divisible by some `spacing_um` value
  (Y1 — validation raises)
- Empty `spacing_um_values` or `dt_s_values` reduced to a single call
  (Y3 — must raise)
- `dt_s == 0` inside the sweep grid (Y3 — boundary case for Phase E
  v1, not part of Item 5 default sweep)
- Wall-clock timestamp in metadata (Y13 — breaks reproducibility)
- Untyped `dict` for primary metadata (Y7 + B1 sister-pattern)
- Local definition of `K_ORIENT_PER_S` or `TRACTION_REF_NN_PER_UM2`
  (must import from `ecm_constitutive_response`, sister-pattern with
  Phase E v1 Y14)
- "varies across axis" qualitative sensitivity claim (Y19 — symmetry
  may produce equal values; assert structural representation +
  analytical formula consequences instead)
- Gate tolerance threshold tuning (Y18 — `1e-12` is IEEE roundoff
  allowance, NOT a tunable threshold; HB#5 invariant is mathematical
  `≤ 1.0`)
- "monotone in spacing" claim (Y8 — sensitivity not guaranteed
  monotone)
- Public export of `_default_phase_e_v1_item5_sweep_config()` helper
  (Y12 — private only)
- Wall-time runtime assertion in tests (Y15 — informational only)

---

## 2. Reasoned-acceptance trace (Y1–Y21)

The lock converged after 4 rounds. Each Y is a Claude concession with
the round in which it was accepted, preserving the adversarial-debate
audit trail (PI id=809 posture).

**Y1 (round 2, accepted Codex C1 = physical-domain invariance contract)**:
- Claude opening lean: vary `spacing_um` while `nx, ny` stay fixed
  (would change physical domain) AND reuse ECM array across spacings
  (would conflate IC with spacing).
- Codex catch: Item 5 measures **grid-spacing sensitivity** —
  resolution must vary with all else held physically invariant. This
  is the Hard Rule 11 measurement-protocol consistency family caught
  also at Phase E v1 Y1 + HB#5 Y1.
- Resolution: fixed `domain_size_um_xy`; derive `nx = domain / spacing`
  with exact-divisibility validation; FA position invariant in physical
  μm; ECM IC sampled per spacing via deterministic physical-coordinate
  function.

**Y2 (round 2, accepted Codex C2 = typed `PhaseEV1Item5SweepConfig`)**:
- Codex: function should take typed scenario/config, not raw
  `ecm_init, fa_init` plus grids.
- Resolution: dataclass with all parameters + harness owns deterministic
  state construction. Custom `ecm_init/fa_init` argument rejected as
  conflict with Y1.

**Y3 (round 2, accepted Codex C3 = empty axis arrays raise, NOT silent
reduce)**:
- Brief said empty axis arrays reduce to a single Phase E v1 call —
  surprising hides errors.
- Resolution: `_validate_config` raises on empty/non-positive/nonfinite
  values.

**Y4 (round 2, accepted Codex C4 = ECM-feedback substep wording, NOT
imaging interval)**:
- Brief implied `dt_s` aligned with PI imaging interval (15 min/60 min).
- Codex: `dt_s` here is the integration substep, not observation
  interval. Frame range `{30, 60, 120, 300}` s as substep sensitivity.
- Resolution: §0 wording + per-step `max_convex_weight` recording so
  large-dt dynamics compression is visible.

**Y5 (round 2, accepted Codex C5 = grid table from divisibility)**:
- `domain_size_um_xy = (16.0, 16.0)` + `spacing ∈ {0.5, 1.0, 2.0, 4.0}`
  → grids `{32, 16, 8, 4}`. 4×4 is intentional limit case for
  extrapolation behavior.

**Y6 (round 2, accepted Codex C6 = typed result with `(n_runs,
n_steps+1)` trajectories)**:
- Step 0 (initial) + per-step values; per-run metadata arrays of length
  `n_runs`.

**Y7 (round 2, accepted Codex C7 = typed metadata, no loose dict)**:
- B1 sister-pattern: typed dataclass throughout; `to_json_dict()`
  helper if file-write needed (separate concern).

**Y8 (round 2, accepted Codex C8 = raw evidence summaries, NO
satisfaction/monotone claims)**:
- `max_bound_ratio`, `max_neutral_multiplier_deviation`,
  `peak_v_active_um2`, `final_v_active_um2`, `max_convex_weight` per
  run — raw observables, no pass/fail decision.

**Y9 (round 3, accepted Codex A1 = concrete `FocalAdhesionState`
builder, schema-explicit)**:
- Schema verified at `acs/v2/focal_adhesion.py:51-85`: `adhesion_id`,
  `cell_id`, `position_um_xy`, `age_s`, `maturity`, `bound_fraction`,
  `state ∈ Literal[...]`, `traction_force_nN_xy`,
  `linked_protrusion_id`, `source`. Builder uses `state="mature"`,
  `maturity=1.0`, `bound_fraction=1.0`, `source="simulated"`.

**Y10 (round 3, accepted Codex A2 = builder return `(fa,)` only)**:
- Round 2 pseudocode `return ((fa,), ...)` invalid (`...` ellipsis was
  placeholder). Phase E v1 derives traction internally via HB#3;
  scenario builder returns adhesions tuple only.

**Y11 (round 3, accepted Codex A3 = `0.5*I` IC reasoning Magic-Number
Block-compliant)**:
- Midpoint between zero orientation and unit diagonal schema bound;
  not fitted; chosen for symmetric non-target IC giving the metric room
  to evolve toward any traction direction. Test 1 / Test 2 / Test 3 of
  Magic-Number Block all pass.

**Y12 (round 3, accepted Codex A4 = default helper PRIVATE only)**:
- Public exports: config + result + metadata + harness function (4
  total). `_default_*` helper is private; tests construct config
  explicitly.

**Y13 (round 3, accepted Codex A5 = no wall-clock timestamp in
metadata)**:
- Two identical runs differ if `datetime.now()` recorded; breaks
  reproducibility. Caller can wrap with own timestamp if file-write
  concerns.

**Y14 (round 3, accepted Codex A6 = `git_sha_override` injectable +
fallback-safe `_detect_git_sha`)**:
- Tests pass fixed override for determinism; runtime detection
  best-effort, returns `"unknown"` on any failure.

**Y15 (round 3, accepted Codex A7 = drop exact-second runtime
estimates)**:
- Round 2 estimates "~30s, ~10s, ~3s, ~1s" un-benchmarked.
- Resolution: lock keeps cell counts + ≤ 10 min design budget; actual
  runtime measured informationally, NOT asserted in tests.

**Y16 (round 3, accepted Codex A8 = `max_bound_ratio` zero-active-bound
case)**:
- `_compute_max_bound_ratio` uses active-bound-positive mask; returns
  0.0 if no active samples (no division by zero).

**Y17 (round 4, accepted Codex B1 = step 0 semantics explicit)**:
- Step 0 = initial ECM measured against step-0 traction target with
  `dt_s=0`. NOT a zero-traction baseline. `v_active_um2[run, 0]`,
  `v_active_bound_um2[run, 0]`, `max_traction_norm[run, 0]` may be
  nonzero (reflect step-0 scenario); `max_convex_weight[run, 0]` and
  `max_orientation_delta[run, 0]` are exactly 0 (because `dt_s=0` →
  `w = -expm1(0) = 0`).

**Y18 (round 4, accepted Codex B2 = IEEE roundoff `1e-12`, NOT gate
threshold)**:
- Test 10 (`v_active_bound invariant evidence`) uses `1e-12` allowance,
  explicitly labeled as IEEE roundoff. HB#5 mathematical invariant
  remains `≤ 1.0`; the margin is for float64 accumulation, NOT a
  tunable gate.

**Y19 (round 4, accepted Codex B3 = evidence tests as structural
representation + analytical-formula consequence, NOT qualitative
"varies")**:
- Test 8 (spacing): assert all spacings present + `nx/ny` derived
  correctly + trajectory shapes; NOT "values vary across spacings".
- Test 9 (dt): assert step-1 `max_convex_weight` monotone non-decreasing
  in dt at fixed spacing — analytical from `w = 1 − exp(−K·S·dt)`,
  K, S > 0.

**Y20 (round 4, accepted Codex B4 = runtime Literal validation +
`n_revolutions` positive int)**:
- `traction_scenario` Literal at type level doesn't prevent runtime
  bad values; validate `== "rotating_uniform_single_fa"` at runtime.
- `n_revolutions` must be `int` (not `bool` per Python `bool ⊂ int`
  trap, HB#1+#2 sister) and `> 0`.
- `traction_magnitude_nN` finite positive.

**Y21 (round 4, accepted Codex B5 = domain/origin/FA-position
defense-in-depth validation)**:
- `domain_size_um_xy` length-2 + finite + positive
- `origin_um_xy` length-2 + finite
- `fa_position_um_xy` length-2 + finite + inside inclusive physical
  domain `[ox, ox+dx] × [oy, oy+dy]`

---

## 3. Sanity Gate (impl writes in module docstring or sibling sanity doc)

1. **Units**:
   - All inputs carry units per their dataclass annotations
     (`spacing_um` μm, `dt_s` s, `traction_magnitude_nN` nN,
     `domain_size_um_xy` μm, `fa_position_um_xy` μm)
   - Outputs: trajectories carry HB#5/HB#1+#2 units (V_active μm²,
     traction nN/μm², dimensionless multipliers/weights)
   - `max_bound_ratio` dimensionless (V_active/V_bound)
   - `git_sha` string identifier; numpy/python_version strings
2. **Boundary**:
   - Empty axis arrays: `_validate_config` raises
   - Non-positive / nonfinite values: raise
   - `n_revolutions == 0`: raise (constant-target ≠ locked rotating
     scenario)
   - `n_revolutions` boolean: raise (Python `bool ⊂ int` trap)
   - Indivisible domain: raise with explicit error showing the bad
     spacing and resulting non-integer nx
   - FA position outside physical domain: raise
   - Step 0 with `dt_s=0`: identity update; `max_convex_weight == 0`
     exact (Y17)
   - All-zero traction (degenerate scenario): `max_bound_ratio == 0.0`,
     no division by zero (Y16)
3. **Conservation**: harness wraps Phase E v1 with no per-step
   conservation invariant of its own. Each sub-call enforces its lock's
   conservation. The harness's structural invariants:
   - Trajectory shape `(n_runs, n_steps+1)` matches config
   - `metadata.config` is the verbatim input (deterministic
     reproducibility)
4. **Numerical**:
   - `np.float64` throughout; `np.int64` for `nx_by_run`/`ny_by_run`
   - `1e-12` IEEE roundoff allowance for `max_bound_ratio` invariant
     (Y18 — explicitly NOT a gate threshold)
   - `_detect_git_sha` 2.0s subprocess timeout to avoid hanging on
     systemic git issues
5. **Sign**: per-cell measurements ≥ 0 by HB#5 construction;
   `max_bound_ratio ≥ 0` by construction (V_active ≥ 0); aggregate
   summaries inherit sign from underlying observables.
6. **Measurement-protocol** (Hard Rule 11, the central anchor):
   - **Spacing axis sensitivity (Y1)**: physical domain + FA position
     invariant; only resolution changes. ECM IC is a physical-coordinate
     function sampled per grid (NOT array reuse). This avoids the
     spacing/IC/domain conflation that would render sensitivity
     measurement meaningless.
   - **dt axis sensitivity (Y4)**: dt is the ECM-feedback integration
     substep, NOT the experimental imaging interval. Frame chosen range
     accordingly in §0 wording + lock §6 references.
   - **Step 0 semantics (Y17)**: initial ECM against step-0 traction
     target, NOT a zero-traction baseline. Step 0 makes initial
     misalignment visible.
   - **Evidence wording (Y8 + Y19 + sister-pattern Phase E v1 Y1+Y2)**:
     `provides_*_evidence` test names; raw observable summaries; no
     "satisfies", no "monotone-in-spacing" claim.
   - **IEEE roundoff vs gate threshold (Y18)**: `1e-12` margin labeled
     as float64 roundoff allowance, not a tunable gate. HB#5 invariant
     remains mathematical `≤ 1.0`.

---

## 4. Test Catalog (23 tests)

### Validation tests (11)

1. `test_item_5_sweep_empty_spacing_values_raises`
2. `test_item_5_sweep_empty_dt_values_raises`
3. `test_item_5_sweep_negative_or_nonfinite_value_raises`
4. `test_item_5_sweep_indivisible_domain_raises` (Y1)
5. `test_item_5_sweep_unsupported_traction_scenario_raises` (Y20)
6. `test_item_5_sweep_zero_n_revolutions_raises` (Y20)
7. `test_item_5_sweep_negative_n_revolutions_raises` (Y20)
8. `test_item_5_sweep_bool_n_revolutions_raises` (Y20 — Python
   `bool ⊂ int` trap, HB#1+#2 sister)
9. `test_item_5_sweep_invalid_domain_size_raises` (Y21)
10. `test_item_5_sweep_invalid_origin_raises` (Y21)
11. `test_item_5_sweep_fa_position_outside_domain_raises` (Y21)

### Determinism + metadata (3)

12. `test_item_5_sweep_metadata_uses_git_sha_override_when_provided` (Y14)
13. `test_item_5_sweep_metadata_falls_back_to_unknown_on_git_failure`
    (Y14 — monkeypatch subprocess to raise)
14. `test_item_5_sweep_no_timestamp_in_metadata` (Y13 — `not
    hasattr(metadata, "timestamp_iso")` and similar)

### Evidence-providing tests (4) — Phase E v1 Y2 sister-pattern

15. `test_item_5_sweep_provides_grid_spacing_sensitivity_evidence`
    (Y19 structural representation: all locked spacings present;
    derived `nx_by_run`/`ny_by_run` match `domain / spacing` exactly;
    trajectories have correct shape per spacing run)
16. `test_item_5_sweep_provides_dt_substep_sensitivity_evidence`
    (Y19 analytical: at fixed spacing, step-1 `max_convex_weight`
    monotone non-decreasing in dt — analytical consequence of
    `w = 1 − exp(−K·S·dt)`)
17. `test_item_5_sweep_provides_v_active_bound_invariant_evidence`
    (Y18 — `max_bound_ratio ≤ 1.0 + 1e-12` everywhere; `1e-12` is
    IEEE roundoff, NOT gate threshold)
18. `test_item_5_sweep_provides_neutral_multiplier_preservation_evidence`
    (`max_neutral_multiplier_deviation == 0.0` everywhere — Phase E v1
    HB#4 neutral hard-wired invariant)

### Boundary + step-0 (3)

19. `test_item_5_sweep_max_bound_ratio_zero_when_no_active_cells` (Y16
    — degenerate all-zero-traction scenario; no division by zero)
20. `test_item_5_sweep_trajectory_shape_n_runs_n_steps_plus_1` (Y6 —
    explicit shape verification)
21. `test_item_5_sweep_step_0_initial_against_target_not_zero_baseline`
    (Y17 — assert `v_active[:, 0] > 0` for misaligned IC;
    `max_convex_weight[:, 0] == 0` exact;
    `max_orientation_delta[:, 0] == 0` exact)

### Composition guard + exports (2)

22. `test_item_5_sweep_no_local_failure_kinds_or_phase_e_error_redefinition`
    (sister-pattern with Phase E v1 test 3 source/meta guard — assert
    no `Item5*Error` class definition + no `failure_kind =` assignment
    in module body)
23. `test_item_5_sweep_exports_through_both_init` — `PhaseEV1Item5SweepConfig`,
    `PhaseEV1Item5SweepMetadata`, `PhaseEV1Item5SweepResult`,
    `run_phase_e_v1_sensitivity_sweep` importable from BOTH
    `acs.v2.dynamics` AND `acs.v2`; private `_default_*`,
    `_build_*`, `_validate_*`, `_detect_*`, `_compute_*` helpers NOT
    exported.

---

## 5. Files

- `docs/v2/v2_item_5_sweep_harness_brief.md` (existing, opening brief,
  commit `4f0efdd`)
- `docs/v2/v2_item_5_sweep_harness_locked.md` (this file, source of truth)
- `acs/v2/dynamics/closed_loop_phase_e_sweep.py` (NEW):
  - `PhaseEV1Item5SweepConfig` dataclass
  - `PhaseEV1Item5SweepMetadata` dataclass
  - `PhaseEV1Item5SweepResult` dataclass
  - `run_phase_e_v1_sensitivity_sweep` function
  - Private helpers: `_validate_config`, `_build_rotating_uniform_single_fa`,
    `_build_initial_ecm_at_spacing`, `_detect_git_sha`,
    `_compute_max_bound_ratio`, `_record_step`,
    `_default_phase_e_v1_item5_sweep_config`
  - Imports: HB#1+#2 constants (NOT redefined), Phase E v1 result +
    function, FocalAdhesionState, ECMSubstrateState
  - **NO** local constants, **NO** new failure kinds, **NO**
    `effective_stiffness` references (sister-pattern with Phase E v1)
- `acs/v2/dynamics/__init__.py` + `acs/v2/__init__.py` — 4 new exports
- `tests/v2/test_v2_closed_loop_phase_e_sweep.py` (NEW, **23 tests** per §4)

---

## 6. References

- Brief (opening position):
  `docs/v2/v2_item_5_sweep_harness_brief.md` (commit `4f0efdd`)
- Parent locked plan:
  `docs/v2/v2_closed_loop_ecm_gate_phased_plan_locked.md` (Item 5 sweep
  is a separate harness from Phase E v1, per parent §1)
- Sister locks:
  - `docs/v2/v2_phase_e_composition_locked.md` (Phase E v1, the harness
    target)
  - `docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`
  - `docs/v2/v2_hard_blocker_1_2_constitutive_direction_locked.md`
  - `docs/v2/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`
  - `docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`
  - `docs/v2/v2_phase_d_no_op_scaffolding_locked.md`
  - `docs/v2/v2_focal_adhesion_dynamics_result_typed_locked.md` (B1)
- Verified at source for Y9 (FocalAdhesionState schema):
  `acs/v2/focal_adhesion.py:51-85`
- Sister implementation precedent (composition wrapper):
  `acs/v2/dynamics/closed_loop_phase_e.py`
- Adversarial debate posture: memory
  `feedback_aggressive_design_debate.md`, PI id=809
- Hard Rule 10 inline derivation: memory
  `rule10_unit_derivation_in_docs.md`, Codex id=1234
- Hard Rule 11 wording-boundary meta-test: memory
  `hard_rule_11_wording_boundary_meta_test.md` (this lock is the
  **third** Phase E composition Hard Rule 11 catch — Phase E v1 Y1
  was first; HB#5 Y1 was second; Item 5 Y1+Y4+Y17 here is the third
  family of measurement-protocol corrections)
- Magic-Number Block compliance for `0.5*I` IC: see Y11 reasoning in
  §2; CLAUDE.md Magic-Number Block section
- Step 6 pre-commit batch: memory `design_note_pre_commit_batch.md`,
  Codex id=1428

---

## 7. Cross-room dispatch

This file is the design-team input to implementation-work for:

1. impl Claude writes Item 5 sweep harness Sanity Gate doc
   (`docs/v2/v2_item_5_sweep_harness_sanity_gate.md`) from this lock — 6
   sections per §3, with Hard Rule 11 measurement-protocol section as
   the central anchor (the third Phase E composition catch).
2. impl Codex review (6 focus per Codex `id=1575` final SEAL):
   - **Physical-domain invariance**: spacing sweep with fixed domain +
     FA position; ECM IC sampled per grid (NOT array reuse)
   - **Typed result/metadata throughout**: no loose dicts; no
     wall-clock timestamp; git_sha override + fallback-safe detection
   - **Evidence-only wording**: `provides_*_evidence` test names;
     `1e-12` IEEE roundoff explicitly labeled NOT gate threshold;
     no "monotone-in-spacing" claim
   - **Step 0 semantics**: initial ECM against step-0 traction target;
     `max_convex_weight[:, 0] == 0` exact (Y17)
   - **Schema-concrete builders**: `_build_rotating_uniform_single_fa`
     uses all FocalAdhesionState fields explicitly + `fa.validate()`;
     `_build_initial_ecm_at_spacing` samples physical-coordinate IC
   - **4 public exports only**: `PhaseEV1Item5SweepConfig`,
     `PhaseEV1Item5SweepMetadata`, `PhaseEV1Item5SweepResult`,
     `run_phase_e_v1_sensitivity_sweep` through both
     `acs.v2.dynamics` and `acs.v2`; private helpers stay private
3. On Sanity Gate PASS: Item 5 sweep harness code commit
   (`acs/v2/dynamics/closed_loop_phase_e_sweep.py` new + 2 export
   updates + 1 new test file).
4. After commit: design-discussion idle on Item 5 sweep harness;
   **Phase E composition + Item 5 evidence channel both sealed**. Next
   design entries (separate cycles):
   - Phase E v2 (HB#4-active variant) design — requires
     `compute_ecm_to_fa_bias_active` design + lock first
   - effective_stiffness law decision (separate unit, only if v2
     consumes it)

Rounds 1–4 of the design lock are MCP id 1564–1575.
