"""V2 Phase E v1 Item 5 sweep harness — ECM-side sensitivity evidence.

Implements the locked Item 5 sweep harness design from
``docs/v2/v2_item_5_sweep_harness_locked.md`` (commit ``c5634d5``)
and the impl-work Sanity Gate
``docs/v2/v2_item_5_sweep_harness_sanity_gate.md`` (commits
``338be1b`` + ``19825e9`` + ``774d963``; Codex Sanity Gate PASS
at ``id=1587``).

Module-docstring intent guard (paraphrased per Phase E v1
``id=1546``/``id=1549`` + Codex ``id=1581`` test 24 wording-
boundary precedent): this harness produces ECM-side sensitivity
evidence for the closed-loop ECM gate Item 5 by sweeping
Phase E v1 across grid spacing and ECM-feedback substep dt.
Whether the observed evidence supports the closed-loop ECM
gate Item 5 conclusion is a separate caller-side / PI /
explicit-decision unit. The harness reports raw observables
without any pass/fail decision built in.

This is the **third Phase E composition Hard Rule 11
measurement-protocol catch family**:
- Phase E v1 Y1 was the first (ECM-side evidence wording)
- HB#5 Y1 was the second (active-cells-only V_active)
- Item 5 Y1 + Y4 + Y17 (this module): physical-domain
  invariance under spacing sweep + ECM-feedback substep
  wording (NOT imaging interval) + Step 0 semantics
  (initial-against-target, NOT zero baseline)

The module deliberately does **not**:

- claim closed-loop ECM gate Item 5 has been resolved by this
  harness anywhere in module/function/return-class docstrings
  or test names (Y8 + Phase E v1 Y2 sister-pattern; meta-test
  24 enforces by string-matching forbidden literal phrases —
  this harness produces evidence only, separate caller-side /
  PI / explicit-decision unit determines whether Item 5 is
  resolved);
- reuse ECM array across spacings (Y1 — would conflate
  resolution with IC change; physical-coordinate IC sampled
  per grid via ``_build_initial_ecm_at_spacing``);
- vary FA position across spacings (Y1 — physical position
  invariant; ``_validate_config`` raises if outside domain);
- accept ``domain_size_um_xy`` not exactly divisible by some
  ``spacing_um`` value (Y1 — ``_validate_config`` raises);
- reduce empty axis arrays to a single Phase E v1 call (Y3
  — must raise);
- include ``dt_s == 0`` inside the sweep grid (Y3 — boundary
  case for Phase E v1, not part of Item 5 default sweep;
  ``_validate_config`` raises);
- record wall-clock timestamp in metadata (Y13 — breaks
  reproducibility; caller can wrap with own timestamp);
- use untyped dict for primary metadata (Y7 + B1 sister-
  pattern; all dataclasses are frozen+slots);
- redefine ``K_ORIENT_PER_S`` or ``TRACTION_REF_NN_PER_UM2``
  locally (Y14 — re-imported from
  ``ecm_constitutive_response``; sister-pattern with
  Phase E v1 Y14);
- claim "varies across axis" qualitative sensitivity (Y19 —
  symmetry may produce equal values; structural representation
  + analytical formula consequences instead);
- tune the ``1e-12`` IEEE roundoff allowance (Y18 — float64
  numerical margin around the mathematical ``≤ 1.0`` invariant,
  NOT a tunable gate threshold);
- claim "monotone in spacing" (Y8 — sensitivity not guaranteed
  monotone);
- expose ``_default_phase_e_v1_item5_sweep_config`` helper
  publicly (Y12 — private only);
- assert wall-time runtime in tests (Y15 — informational only).

Validation contract: ``_validate_config`` runs defense-in-depth
checks per Y3 (empty arrays, non-positive/nonfinite values) +
Y20 (Literal scenario, n_revolutions, magnitude) + Y21
(domain/origin/FA-position) + Y1 (exact divisibility). All
validation failures raise plain ``ValueError`` (sister-pattern
with Phase E v1 + HB#5 — no domain-specific failure modes for
a thin harness).

Sanity Gate scope (acs/v2/dynamics/closed_loop_phase_e_sweep.py):

- §1 dimensional: harness performs no unit transformation;
  inherits Phase E v1 + HB#5 unit chains.
- §2 boundary: 8 boundary classes covered by 16 of 26 tests
  (validation 11 + Step 0 + zero-bound + invalid-cases +
  divisibility-fix regressions per Codex id=1591/id=1594).
- §3 conservation: harness has no per-step conservation of
  its own; trajectory shape exactly per config; metadata
  reproducibility (Y13 + Y14).
- §4 numerical: float64 + int64; ``1e-12`` IEEE roundoff
  allowance explicitly NOT gate tunable; zero new harness
  tolerances.
- §5 sign: all 7 trajectory + aggregate fields non-negative
  by construction.
- §6 measurement-protocol (Hard Rule 11, central anchor):
  six-layer guard — locked §0 wording + Y1 physical-domain
  invariance + Y4 ECM-feedback substep wording + Y17 Step 0
  semantics + Y18 IEEE-roundoff-vs-gate-threshold separation
  + Y8/Y19 evidence-only test names; meta-test 22 (failure-
  kind discipline) and meta-test 24 (wording-boundary
  discipline) per Codex ``id=1581`` test separation.

Magic-Number Block: 2 locked harness-level constants
(``0.5 * I`` orientation IC per Y11; ``1e-12`` IEEE roundoff
allowance per Y18); both 3-test compliant.
"""

from __future__ import annotations

import math
import subprocess
import sys
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

import numpy as np

from acs.v2.dynamics.closed_loop_phase_e import (
    PhaseEStepResult,
    step_closed_loop_phase_e_v1,
    step_closed_loop_phase_e_v2,
)
from acs.v2.dynamics.ecm_constitutive_response import (
    K_ORIENT_PER_S,
    TRACTION_REF_NN_PER_UM2,
)
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState


_GIT_SUBPROCESS_TIMEOUT_S: float = 2.0


@dataclass(frozen=True, slots=True)
class PhaseEV1Item5SweepConfig:
    """Typed config for one Item 5 sweep harness invocation.

    Y1 contract: physical domain + FA position invariant across
    spacing sweep; only resolution changes. Domain must be
    exactly divisible by every sweep spacing
    (``_validate_config`` raises otherwise).

    Attributes:
        spacing_um_values: NON-empty tuple of finite positive
            ``spacing_um`` values to sweep.
        dt_s_values: NON-empty tuple of finite positive
            ``dt_s`` values to sweep (ECM-feedback integration
            substep, NOT imaging interval per Y4).
        n_steps: positive integer number of HB#1+#2 update
            steps per run; total trajectory length is
            ``n_steps + 1`` to include step 0 initial.
        domain_size_um_xy: physical domain size in μm; finite
            positive 2-vector.
        origin_um_xy: physical domain origin in μm; finite
            2-vector.
        traction_scenario: Literal-restricted scenario name;
            currently only ``"rotating_uniform_single_fa"``
            (Y20 runtime validation).
        traction_magnitude_nN: finite positive traction
            magnitude in nN.
        fa_position_um_xy: FA physical position in μm; finite
            2-vector inside the inclusive physical domain
            (Y21 validation).
        n_revolutions: positive int (NOT bool per Y20 Python
            ``bool ⊂ int`` trap); number of full revolutions
            of traction direction across n_steps.
        git_sha_override: optional fixed git SHA for tests
            (Y14 deterministic reproducibility); None →
            best-effort runtime detection.
    """

    spacing_um_values: tuple[float, ...]
    dt_s_values: tuple[float, ...]
    n_steps: int
    domain_size_um_xy: tuple[float, float]
    origin_um_xy: tuple[float, float]
    traction_scenario: Literal["rotating_uniform_single_fa"]
    traction_magnitude_nN: float
    fa_position_um_xy: tuple[float, float]
    n_revolutions: int = 1
    git_sha_override: str | None = None


@dataclass(frozen=True, slots=True)
class PhaseEV1Item5SweepMetadata:
    """Typed run metadata for reproducibility (Y13 — NO wall-clock timestamp)."""

    config: PhaseEV1Item5SweepConfig
    git_sha: str
    numpy_version: str
    python_version: str
    k_orient_per_s: float
    traction_ref_nN_per_um2: float


@dataclass(frozen=True, slots=True)
class PhaseEV1Item5SweepResult:
    """Typed sweep result. Trajectory shapes (n_runs, n_steps+1) include step 0."""

    spacing_um_by_run: np.ndarray
    dt_s_by_run: np.ndarray
    nx_by_run: np.ndarray
    ny_by_run: np.ndarray
    v_active_um2_traj: np.ndarray
    v_active_bound_um2_traj: np.ndarray
    max_convex_weight_traj: np.ndarray
    max_orientation_delta_frobenius_traj: np.ndarray
    max_traction_norm_nN_per_um2_traj: np.ndarray
    neutral_multiplier_max_deviation_traj: np.ndarray
    max_bound_ratio: np.ndarray
    max_neutral_multiplier_deviation: np.ndarray
    peak_v_active_um2: np.ndarray
    final_v_active_um2: np.ndarray
    metadata: PhaseEV1Item5SweepMetadata


@dataclass(frozen=True, slots=True)
class PhaseEV2Item5SweepConfig:
    """Typed config for the Phase E v2 Item 5 sweep harness variant.

    This B-tier sister harness reuses the Phase E v1 sweep geometry and
    validation contract, but swaps the per-step composition call to
    :func:`acs.v2.dynamics.closed_loop_phase_e.step_closed_loop_phase_e_v2`.
    It records active HB#4 multiplier summaries as evidence only; it
    does not claim Item 5 is resolved.
    """

    spacing_um_values: tuple[float, ...]
    dt_s_values: tuple[float, ...]
    n_steps: int
    domain_size_um_xy: tuple[float, float]
    origin_um_xy: tuple[float, float]
    traction_scenario: Literal["rotating_uniform_single_fa"]
    traction_magnitude_nN: float
    fa_position_um_xy: tuple[float, float]
    k_active: float
    n_revolutions: int = 1
    git_sha_override: str | None = None


@dataclass(frozen=True, slots=True)
class PhaseEV2Item5SweepMetadata:
    """Typed v2 sweep metadata for reproducibility; no wall-clock timestamp."""

    config: PhaseEV2Item5SweepConfig
    git_sha: str
    numpy_version: str
    python_version: str
    k_orient_per_s: float
    traction_ref_nN_per_um2: float


@dataclass(frozen=True, slots=True)
class PhaseEV2Item5SweepResult:
    """Typed v2 sweep result; raw evidence only, no pass/fail decision."""

    spacing_um_by_run: np.ndarray
    dt_s_by_run: np.ndarray
    nx_by_run: np.ndarray
    ny_by_run: np.ndarray
    v_active_um2_traj: np.ndarray
    v_active_bound_um2_traj: np.ndarray
    max_convex_weight_traj: np.ndarray
    max_orientation_delta_frobenius_traj: np.ndarray
    max_traction_norm_nN_per_um2_traj: np.ndarray
    active_multiplier_min_traj: np.ndarray
    active_multiplier_mean_traj: np.ndarray
    active_multiplier_max_traj: np.ndarray
    active_multiplier_max_deviation_traj: np.ndarray
    max_bound_ratio: np.ndarray
    peak_active_multiplier_max_deviation: np.ndarray
    peak_v_active_um2: np.ndarray
    final_v_active_um2: np.ndarray
    metadata: PhaseEV2Item5SweepMetadata


def _validate_config(config: PhaseEV1Item5SweepConfig) -> None:
    """Defense-in-depth config validation (Y3 + Y20 + Y21 + Y1)."""

    if len(config.spacing_um_values) == 0:
        raise ValueError("spacing_um_values must be non-empty")
    if len(config.dt_s_values) == 0:
        raise ValueError("dt_s_values must be non-empty")

    for spacing in config.spacing_um_values:
        if isinstance(spacing, bool) or not isinstance(spacing, (int, float)):
            raise ValueError(
                f"spacing_um_values must contain finite positive floats; "
                f"got {type(spacing).__name__} {spacing!r}"
            )
        if not math.isfinite(float(spacing)) or float(spacing) <= 0.0:
            raise ValueError(
                f"spacing_um_values must contain finite positive floats; "
                f"got {spacing!r}"
            )

    for dt_s in config.dt_s_values:
        if isinstance(dt_s, bool) or not isinstance(dt_s, (int, float)):
            raise ValueError(
                f"dt_s_values must contain finite positive floats; "
                f"got {type(dt_s).__name__} {dt_s!r}"
            )
        if not math.isfinite(float(dt_s)) or float(dt_s) <= 0.0:
            raise ValueError(
                f"dt_s_values must contain finite positive floats; "
                f"got {dt_s!r}"
            )

    if (
        isinstance(config.n_steps, bool)
        or not isinstance(config.n_steps, int)
        or config.n_steps <= 0
    ):
        raise ValueError(
            f"n_steps must be a positive int (not bool); got "
            f"{type(config.n_steps).__name__} {config.n_steps!r}"
        )

    if config.traction_scenario != "rotating_uniform_single_fa":
        raise ValueError(
            f"traction_scenario must be 'rotating_uniform_single_fa'; "
            f"got {config.traction_scenario!r}"
        )

    if (
        isinstance(config.n_revolutions, bool)
        or not isinstance(config.n_revolutions, int)
        or config.n_revolutions <= 0
    ):
        raise ValueError(
            f"n_revolutions must be a positive int (not bool); got "
            f"{type(config.n_revolutions).__name__} {config.n_revolutions!r}"
        )

    if (
        isinstance(config.traction_magnitude_nN, bool)
        or not isinstance(config.traction_magnitude_nN, (int, float))
        or not math.isfinite(float(config.traction_magnitude_nN))
        or float(config.traction_magnitude_nN) <= 0.0
    ):
        raise ValueError(
            f"traction_magnitude_nN must be a finite positive float; got "
            f"{config.traction_magnitude_nN!r}"
        )

    if (
        len(config.domain_size_um_xy) != 2
        or any(
            isinstance(v, bool)
            or not isinstance(v, (int, float))
            or not math.isfinite(float(v))
            or float(v) <= 0.0
            for v in config.domain_size_um_xy
        )
    ):
        raise ValueError(
            f"domain_size_um_xy must be a length-2 finite positive tuple; "
            f"got {config.domain_size_um_xy!r}"
        )

    if (
        len(config.origin_um_xy) != 2
        or any(
            isinstance(v, bool)
            or not isinstance(v, (int, float))
            or not math.isfinite(float(v))
            for v in config.origin_um_xy
        )
    ):
        raise ValueError(
            f"origin_um_xy must be a length-2 finite tuple; got "
            f"{config.origin_um_xy!r}"
        )

    if (
        len(config.fa_position_um_xy) != 2
        or any(
            isinstance(v, bool)
            or not isinstance(v, (int, float))
            or not math.isfinite(float(v))
            for v in config.fa_position_um_xy
        )
    ):
        raise ValueError(
            f"fa_position_um_xy must be a length-2 finite tuple; got "
            f"{config.fa_position_um_xy!r}"
        )
    ox = float(config.origin_um_xy[0])
    oy = float(config.origin_um_xy[1])
    dx = float(config.domain_size_um_xy[0])
    dy = float(config.domain_size_um_xy[1])
    fx = float(config.fa_position_um_xy[0])
    fy = float(config.fa_position_um_xy[1])
    if not (ox <= fx <= ox + dx) or not (oy <= fy <= oy + dy):
        raise ValueError(
            f"fa_position_um_xy {config.fa_position_um_xy!r} must lie inside "
            f"physical domain [{ox}, {ox + dx}] x [{oy}, {oy + dy}]"
        )

    for spacing in config.spacing_um_values:
        spacing_f = float(spacing)
        # Y1 exact divisibility per Codex id=1591 + id=1594 BLOCKER fix:
        # use exact decimal representation via Fraction(str(value)).
        # NOT float64 binary-quotient semantics (`x / y; .is_integer()`)
        # — that rejects mathematically exact decimal values like
        # 0.3 / 0.1 = 2.9999999999999996 (because 0.1 is not exactly
        # representable in float64). For a YAML/user-facing float
        # config, "exact divisibility" means exact decimal divisibility
        # (user-intended), not binary divisibility.
        # Fraction(str(value)) parses via Python's shortest round-
        # trippable repr → exact decimal interpretation:
        #   - Fraction(str(0.3)) = Fraction(3, 10)
        #   - Fraction(str(0.1)) = Fraction(1, 10)
        #   - Fraction(3, 10) / Fraction(1, 10) = Fraction(3, 1)
        #     → denominator == 1 → exact integer
        # Catches near-indivisible-inside-prior-tolerance:
        #   - Fraction(str(1.0 + 1e-11)) ≈ Fraction(100000000001, 1e11)
        #     → Fraction(16) / Fraction(...) has non-1 denominator
        #     → raises (Codex id=1591 test 25 still PASSes).
        spacing_frac = Fraction(str(spacing_f))
        dx_frac = Fraction(str(dx))
        dy_frac = Fraction(str(dy))
        if spacing_frac == 0:
            raise ValueError(
                f"spacing_um {spacing_f!r} cannot be zero (caught earlier "
                f"by finiteness/positivity check; this guard is defensive)"
            )
        quotient_x = dx_frac / spacing_frac
        quotient_y = dy_frac / spacing_frac
        if quotient_x.denominator != 1 or quotient_y.denominator != 1:
            raise ValueError(
                f"domain_size_um_xy {config.domain_size_um_xy!r} must be "
                f"exactly divisible by every spacing_um; got spacing="
                f"{spacing_f} with quotient_x={quotient_x}, "
                f"quotient_y={quotient_y} (must be exact positive integers)"
            )
        nx = int(quotient_x)
        ny = int(quotient_y)
        if nx < 1 or ny < 1:
            raise ValueError(
                f"domain_size_um_xy {config.domain_size_um_xy!r} produces "
                f"non-positive grid dimensions at spacing={spacing_f}: "
                f"nx={nx}, ny={ny} (must be ≥ 1)"
            )


def _as_phase_e_v1_config(
    config: PhaseEV2Item5SweepConfig,
) -> PhaseEV1Item5SweepConfig:
    """Reuse the locked v1 sweep geometry validation for the v2 harness."""

    return PhaseEV1Item5SweepConfig(
        spacing_um_values=config.spacing_um_values,
        dt_s_values=config.dt_s_values,
        n_steps=config.n_steps,
        domain_size_um_xy=config.domain_size_um_xy,
        origin_um_xy=config.origin_um_xy,
        traction_scenario=config.traction_scenario,
        traction_magnitude_nN=config.traction_magnitude_nN,
        fa_position_um_xy=config.fa_position_um_xy,
        n_revolutions=config.n_revolutions,
        git_sha_override=config.git_sha_override,
    )


def _validate_v2_config(config: PhaseEV2Item5SweepConfig) -> None:
    """Validate v2 sweep config while preserving v1 sister contracts."""

    _validate_config(_as_phase_e_v1_config(config))


def _build_rotating_uniform_single_fa(
    config: PhaseEV1Item5SweepConfig | PhaseEV2Item5SweepConfig,
    step_index: int,
    dt_s: float,
) -> tuple[FocalAdhesionState, ...]:
    """Single-FA scenario with deterministic rotating traction direction (Y9 + Y10)."""

    theta = 2.0 * math.pi * config.n_revolutions * step_index / config.n_steps
    force_x = config.traction_magnitude_nN * math.cos(theta)
    force_y = config.traction_magnitude_nN * math.sin(theta)
    fa = FocalAdhesionState(
        adhesion_id="item5-fa-0",
        cell_id="item5-cell-0",
        position_um_xy=(
            float(config.fa_position_um_xy[0]),
            float(config.fa_position_um_xy[1]),
        ),
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
    config: PhaseEV1Item5SweepConfig | PhaseEV2Item5SweepConfig,
    spacing_um: float,
) -> ECMSubstrateState:
    """Sample physical-coordinate IC function onto specified grid (Y1 + Y11)."""

    nx = int(round(config.domain_size_um_xy[0] / spacing_um))
    ny = int(round(config.domain_size_um_xy[1] / spacing_um))
    return ECMSubstrateState(
        origin_um_xy=(
            float(config.origin_um_xy[0]),
            float(config.origin_um_xy[1]),
        ),
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
    """Best-effort git SHA detection (Y14)."""

    if override is not None:
        return override
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=_GIT_SUBPROCESS_TIMEOUT_S,
        )
        return result.stdout.strip()
    except (
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
        FileNotFoundError,
        OSError,
    ):
        return "unknown"


def _compute_max_bound_ratio(
    v_active_traj: np.ndarray, v_bound_traj: np.ndarray,
) -> float:
    """Max V_active / V_bound over active-bound-positive samples (Y16)."""

    mask = v_bound_traj > 0.0
    if not mask.any():
        return 0.0
    return float((v_active_traj[mask] / v_bound_traj[mask]).max())


def _record_step(
    step_result: PhaseEStepResult,
    run_idx: int,
    step_idx: int,
    v_active_traj: np.ndarray,
    v_bound_traj: np.ndarray,
    weight_traj: np.ndarray,
    delta_traj: np.ndarray,
    traction_norm_traj: np.ndarray,
    multiplier_dev_traj: np.ndarray,
) -> None:
    """Write one step's diagnostics into the trajectory arrays."""

    diag_lyap = step_result.lyapunov_metric.diagnostics
    diag_orient = step_result.orientation_response.diagnostics
    multipliers = step_result.ecm_to_fa_bias.multipliers_per_fa

    v_active_traj[run_idx, step_idx] = step_result.lyapunov_metric.v_active_um2
    v_bound_traj[run_idx, step_idx] = diag_lyap.v_active_max_bound_um2
    weight_traj[run_idx, step_idx] = diag_orient.max_convex_weight
    delta_traj[run_idx, step_idx] = diag_orient.max_orientation_delta_frobenius
    traction_norm_traj[run_idx, step_idx] = diag_lyap.max_traction_norm_nN_per_um2
    if multipliers.size > 0:
        multiplier_dev_traj[run_idx, step_idx] = float(
            np.abs(multipliers - 1.0).max()
        )
    else:
        multiplier_dev_traj[run_idx, step_idx] = 0.0


def _record_step_v2(
    step_result: PhaseEStepResult,
    run_idx: int,
    step_idx: int,
    v_active_traj: np.ndarray,
    v_bound_traj: np.ndarray,
    weight_traj: np.ndarray,
    delta_traj: np.ndarray,
    traction_norm_traj: np.ndarray,
    multiplier_min_traj: np.ndarray,
    multiplier_mean_traj: np.ndarray,
    multiplier_max_traj: np.ndarray,
    multiplier_dev_traj: np.ndarray,
) -> None:
    """Write one v2 step's diagnostics, including active multiplier summaries."""

    diag_lyap = step_result.lyapunov_metric.diagnostics
    diag_orient = step_result.orientation_response.diagnostics
    multipliers = step_result.ecm_to_fa_bias.multipliers_per_fa

    v_active_traj[run_idx, step_idx] = step_result.lyapunov_metric.v_active_um2
    v_bound_traj[run_idx, step_idx] = diag_lyap.v_active_max_bound_um2
    weight_traj[run_idx, step_idx] = diag_orient.max_convex_weight
    delta_traj[run_idx, step_idx] = diag_orient.max_orientation_delta_frobenius
    traction_norm_traj[run_idx, step_idx] = diag_lyap.max_traction_norm_nN_per_um2
    if multipliers.size > 0:
        multiplier_min_traj[run_idx, step_idx] = float(multipliers.min())
        multiplier_mean_traj[run_idx, step_idx] = float(multipliers.mean())
        multiplier_max_traj[run_idx, step_idx] = float(multipliers.max())
        multiplier_dev_traj[run_idx, step_idx] = float(
            np.abs(multipliers - 1.0).max()
        )
    else:
        multiplier_min_traj[run_idx, step_idx] = 1.0
        multiplier_mean_traj[run_idx, step_idx] = 1.0
        multiplier_max_traj[run_idx, step_idx] = 1.0
        multiplier_dev_traj[run_idx, step_idx] = 0.0


def _default_phase_e_v1_item5_sweep_config() -> PhaseEV1Item5SweepConfig:
    """Default config matching the locked grid table (private; Y12)."""

    return PhaseEV1Item5SweepConfig(
        spacing_um_values=(0.5, 1.0, 2.0, 4.0),
        dt_s_values=(30.0, 60.0, 120.0, 300.0),
        n_steps=100,
        domain_size_um_xy=(16.0, 16.0),
        origin_um_xy=(0.0, 0.0),
        traction_scenario="rotating_uniform_single_fa",
        traction_magnitude_nN=2.0,
        fa_position_um_xy=(8.0, 8.0),
        n_revolutions=1,
        git_sha_override=None,
    )


def run_phase_e_v1_sensitivity_sweep(
    config: PhaseEV1Item5SweepConfig,
) -> PhaseEV1Item5SweepResult:
    """Run cross-product sweep over (spacing_um × dt_s).

    For each ``(spacing_um, dt_s)`` pair:

    1. Build initial ECM at this spacing via deterministic
       physical-coordinate IC function (Y1).
    2. Capture step 0 V_active by calling Phase E v1 with
       ``dt_s=0`` (initial state measured against step-0
       traction target — Y17).
    3. Loop ``n_steps`` with rotating traction direction,
       propagating ``updated_ecm`` forward.
    4. Record per-step trajectories + aggregate summaries.

    Args:
        config: typed :class:`PhaseEV1Item5SweepConfig`. Validated
            via ``_validate_config`` (Y3 + Y20 + Y21 + Y1).

    Returns:
        :class:`PhaseEV1Item5SweepResult` with full trajectories
        + per-run aggregates + typed metadata. NO pass/fail
        decision; raw observable evidence only.
    """

    _validate_config(config)

    n_runs = len(config.spacing_um_values) * len(config.dt_s_values)
    n_steps_plus_1 = config.n_steps + 1

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
            ecm = _build_initial_ecm_at_spacing(config, float(spacing_um))
            spacing_by[run_idx] = float(spacing_um)
            dt_by[run_idx] = float(dt_s)
            nx_by[run_idx] = ecm.grid_shape[0]
            ny_by[run_idx] = ecm.grid_shape[1]

            init_adhesions = _build_rotating_uniform_single_fa(
                config, 0, float(dt_s)
            )
            init_result = step_closed_loop_phase_e_v1(
                init_adhesions, ecm, dt_s=0.0
            )
            _record_step(
                init_result,
                run_idx,
                0,
                v_active_traj,
                v_bound_traj,
                weight_traj,
                delta_traj,
                traction_norm_traj,
                multiplier_dev_traj,
            )

            for step in range(1, config.n_steps + 1):
                adhesions = _build_rotating_uniform_single_fa(
                    config, step, float(dt_s)
                )
                step_result = step_closed_loop_phase_e_v1(
                    adhesions, ecm, dt_s=float(dt_s)
                )
                _record_step(
                    step_result,
                    run_idx,
                    step,
                    v_active_traj,
                    v_bound_traj,
                    weight_traj,
                    delta_traj,
                    traction_norm_traj,
                    multiplier_dev_traj,
                )
                ecm = step_result.updated_ecm

            run_idx += 1

    max_bound_ratio = np.array(
        [
            _compute_max_bound_ratio(v_active_traj[i], v_bound_traj[i])
            for i in range(n_runs)
        ],
        dtype=np.float64,
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


def run_phase_e_v2_sensitivity_sweep(
    config: PhaseEV2Item5SweepConfig,
) -> PhaseEV2Item5SweepResult:
    """Run the Phase E v2 Item 5 sweep harness variant.

    This function is a B-tier sister of
    :func:`run_phase_e_v1_sensitivity_sweep`: it preserves the same
    physical-domain spacing/dt sweep and records raw evidence only, but
    calls :func:`step_closed_loop_phase_e_v2` with the configured
    ``k_active``. It does not make a pass/fail decision and does not
    claim Item 5 is resolved.
    """

    _validate_v2_config(config)

    n_runs = len(config.spacing_um_values) * len(config.dt_s_values)
    n_steps_plus_1 = config.n_steps + 1

    spacing_by = np.empty(n_runs, dtype=np.float64)
    dt_by = np.empty(n_runs, dtype=np.float64)
    nx_by = np.empty(n_runs, dtype=np.int64)
    ny_by = np.empty(n_runs, dtype=np.int64)
    v_active_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    v_bound_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    weight_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    delta_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    traction_norm_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)
    multiplier_min_traj = np.ones((n_runs, n_steps_plus_1), dtype=np.float64)
    multiplier_mean_traj = np.ones((n_runs, n_steps_plus_1), dtype=np.float64)
    multiplier_max_traj = np.ones((n_runs, n_steps_plus_1), dtype=np.float64)
    multiplier_dev_traj = np.zeros((n_runs, n_steps_plus_1), dtype=np.float64)

    run_idx = 0
    for spacing_um in config.spacing_um_values:
        for dt_s in config.dt_s_values:
            ecm = _build_initial_ecm_at_spacing(config, float(spacing_um))
            spacing_by[run_idx] = float(spacing_um)
            dt_by[run_idx] = float(dt_s)
            nx_by[run_idx] = ecm.grid_shape[0]
            ny_by[run_idx] = ecm.grid_shape[1]

            init_adhesions = _build_rotating_uniform_single_fa(
                config, 0, float(dt_s)
            )
            init_result = step_closed_loop_phase_e_v2(
                init_adhesions,
                ecm,
                dt_s=0.0,
                k_active=float(config.k_active),
            )
            _record_step_v2(
                init_result,
                run_idx,
                0,
                v_active_traj,
                v_bound_traj,
                weight_traj,
                delta_traj,
                traction_norm_traj,
                multiplier_min_traj,
                multiplier_mean_traj,
                multiplier_max_traj,
                multiplier_dev_traj,
            )

            for step in range(1, config.n_steps + 1):
                adhesions = _build_rotating_uniform_single_fa(
                    config, step, float(dt_s)
                )
                step_result = step_closed_loop_phase_e_v2(
                    adhesions,
                    ecm,
                    dt_s=float(dt_s),
                    k_active=float(config.k_active),
                )
                _record_step_v2(
                    step_result,
                    run_idx,
                    step,
                    v_active_traj,
                    v_bound_traj,
                    weight_traj,
                    delta_traj,
                    traction_norm_traj,
                    multiplier_min_traj,
                    multiplier_mean_traj,
                    multiplier_max_traj,
                    multiplier_dev_traj,
                )
                ecm = step_result.updated_ecm

            run_idx += 1

    max_bound_ratio = np.array(
        [
            _compute_max_bound_ratio(v_active_traj[i], v_bound_traj[i])
            for i in range(n_runs)
        ],
        dtype=np.float64,
    )

    metadata = PhaseEV2Item5SweepMetadata(
        config=config,
        git_sha=_detect_git_sha(config.git_sha_override),
        numpy_version=np.__version__,
        python_version=sys.version,
        k_orient_per_s=K_ORIENT_PER_S,
        traction_ref_nN_per_um2=TRACTION_REF_NN_PER_UM2,
    )

    return PhaseEV2Item5SweepResult(
        spacing_um_by_run=spacing_by,
        dt_s_by_run=dt_by,
        nx_by_run=nx_by,
        ny_by_run=ny_by,
        v_active_um2_traj=v_active_traj,
        v_active_bound_um2_traj=v_bound_traj,
        max_convex_weight_traj=weight_traj,
        max_orientation_delta_frobenius_traj=delta_traj,
        max_traction_norm_nN_per_um2_traj=traction_norm_traj,
        active_multiplier_min_traj=multiplier_min_traj,
        active_multiplier_mean_traj=multiplier_mean_traj,
        active_multiplier_max_traj=multiplier_max_traj,
        active_multiplier_max_deviation_traj=multiplier_dev_traj,
        max_bound_ratio=max_bound_ratio,
        peak_active_multiplier_max_deviation=multiplier_dev_traj.max(axis=1),
        peak_v_active_um2=v_active_traj.max(axis=1),
        final_v_active_um2=v_active_traj[:, -1],
        metadata=metadata,
    )
