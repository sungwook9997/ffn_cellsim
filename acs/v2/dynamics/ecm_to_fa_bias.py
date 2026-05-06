"""ECM→FA bias target (Hard Blocker #4) — Phase D no-op default.

Implements the locked Hard Blocker #4 design from
``docs/v2/v2_hard_blocker_4_ecm_to_fa_bias_target_locked.md`` and the
pre-execution Sanity Gate
``docs/v2/v2_ecm_to_fa_bias_sanity_gate.md``: a pure read function
that bilinearly samples ECM fields at FA positions plus a pure
neutral function that returns all-1.0 per-FA per-rate multipliers
(Phase D no-op default).

The module deliberately does **not**:

- mutate the input ECM or FA states (pure functions);
- compute non-1.0 multipliers via :func:`compute_ecm_to_fa_bias_neutral`
  (non-1.0 multipliers are available **only** through the explicit
  :func:`compute_ecm_to_fa_bias_active` entry point — locked at
  ``docs/v2/v2_hard_blocker_4_active_locked.md``, Sanity Gated at
  ``docs/v2/v2_hard_blocker_4_active_sanity_gate.md``; function-naming
  separation preserves the silent-activation guard);
- mutate ``traction_scale_nN`` (locked phased plan §3
  effective-stiffness guard);
- spawn or remove FAs (FA nucleation deferred — locked phased plan
  §4);
- use any RNG (deterministic);
- alter FA ``state`` literal (preserves 6.3a no-auto-state contract);
- scalarize the orientation tensor (returned full ``(N_FA, 2, 2)``
  per Hard Rule 11 measurement-protocol consistency);
- carry vector accumulated traction memory (ECM schema is scalar
  ``(nx, ny)`` per ``accumulated_traction_nNs_per_um2``);
- raise on empty FA list (returns shape-consistent empty result).

The HB#4-active variant ``compute_ecm_to_fa_bias_active`` is now also
defined in this module (Phase E v2 step 1; locked at
``docs/v2/v2_hard_blocker_4_active_locked.md``, Sanity Gated at
``docs/v2/v2_hard_blocker_4_active_sanity_gate.md``). Function-naming
separation preserves the silent-activation guard: a caller cannot
trigger active behavior through ``compute_ecm_to_fa_bias_neutral``;
the two variants are distinct entry points with distinct validation
chains and distinct failure kinds.

Sanity Gate scope (acs/v2/dynamics/ecm_to_fa_bias.py):

- §1 dimensional: bilinear linear-combination preserves underlying
  field units; multiplier dimensionless 1.0; no Rule 10 cross-unit
  comparison at this layer.
- §2 boundary: every failure_kind from the Sanity Gate §2 is
  raised before any sampler/multiplier output;
  ``ecm.validate()`` at function entry (local dynamics precedent);
  ``fa.validate()`` per FA after the lock-specific
  ``non_finite_fa_position`` check (mirrors Hard Blocker #3 fix
  ``883efc8`` per Codex review id=1378).
- §3 conservation: pure read / no-mutation; sampler bilinear
  weight conservation matches Hard Blocker #3 (interior trivial,
  boundary half-cell renormalized); Phase D structural multiplier
  invariant 1.0 everywhere.
- §4 numerical: float64; per-FA work O(8) for sampler bilinear;
  no dt-rate gate; ``_FA_POSITION_BOUNDARY_TOL_UM`` recomputed per
  call (matches Hard Blocker #3 grid-invariant pattern).
- §5 sign: bilinear weights non-negative + sum 1.0 (linear
  combination convexity preserves underlying field bounds);
  multipliers exactly 1.0.
- §6 measurement-protocol: locked output modality; no
  scalarization at sampler time; orientation tensor returned full
  ``(N_FA, 2, 2)``; runtime meta-test
  ``test_no_active_law_invoked_in_phase_d`` enforces Hard Rule 11
  boundary (Phase B / Phase C / Hard Blocker #3 precedent).

Magic-Number Block: ``_BOUNDARY_TOL_RELATIVE = 1e-12`` (Hard Blocker
#3 inheritance pattern) and ``_DEVIATORIC_SCORE_BOUND = sqrt(2.0)``
(derived schema bound on ``|n.T @ Q @ n|`` under ECM
``|T_ij| <= 1`` invariant; HB#4-active step 1 lock Y10) are the
two module-level numeric constants. Both are derivable from
upstream invariants — neither is fitted to a target, neither
requires re-tuning under grid changes (Magic-Number Block tests 1
and 2 pass; not chosen to fit any specific datapoint, test 3
pass). ``RATE_NAMES`` is a literal-string tuple, not a numeric
tunable; the all-1.0 neutral multiplier is a structural identity.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final, Literal, Optional

import numpy as np

from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState

_BOUNDARY_TOL_RELATIVE: Final[float] = 1e-12

RATE_NAMES: Final[
    tuple[
        Literal["k_maturity_per_s"],
        Literal["k_bind_per_s"],
        Literal["k_unbind_per_s"],
    ]
] = (
    "k_maturity_per_s",
    "k_bind_per_s",
    "k_unbind_per_s",
)

FAToECMBiasFailureKind = Literal[
    "fa_bias_position_outside_ecm_grid",
    "non_finite_fa_position",
    "k_active_invalid",
    "active_multiplier_non_finite",
]


class FAToECMBiasError(ValueError):
    """Validation error raised by
    :func:`sample_ecm_at_fa_positions` /
    :func:`compute_ecm_to_fa_bias_neutral` /
    :func:`compute_ecm_to_fa_bias_active` when an FA's position
    or active-bias parameter violates the locked bias contract.
    Carries a machine-readable :attr:`failure_kind`."""

    def __init__(self, failure_kind: str, message: str) -> None:
        self.failure_kind = failure_kind
        super().__init__(message)


@dataclass(frozen=True, slots=True)
class ECMSampledAtFAs:
    """Bilinear-interpolated ECM field values at FA positions.

    Diagnostic-only in Phase D; Phase E active mapping (separate
    function ``compute_ecm_to_fa_bias_active`` with its own
    Sanity Gate) consumes these for non-neutral multipliers.

    Attributes preserve the input FA order. Empty FA list yields
    zero-row arrays with shape ``(0, ...)``.
    """

    fa_ids: tuple[str, ...]
    stiffness_kpa: np.ndarray
    fiber_density: np.ndarray
    ligand_density: np.ndarray
    orientation_tensor: np.ndarray
    accumulated_traction_nNs_per_um2: np.ndarray


@dataclass(frozen=True, slots=True)
class ECMToFABiasResult:
    """Per-FA per-rate multipliers + diagnostics for Phase D no-op.

    ``multipliers_per_fa[i, k]`` is the multiplier for FA ``i``
    on rate ``RATE_NAMES[k]``. In Phase D the default
    implementation returns all 1.0 (neutral). Phase E active
    mapping (separate function) may return non-1.0 multipliers.
    """

    fa_ids: tuple[str, ...]
    multipliers_per_fa: np.ndarray
    rate_names: tuple[str, ...]
    sampled_diagnostics: Optional[ECMSampledAtFAs] = None
    diagnostics_dict: dict[str, float | int | str] = field(default_factory=dict)


def _boundary_tol_um(ecm: ECMSubstrateState) -> float:
    """Recompute the numerical tie-break tolerance per ECM call so
    the inclusion check stays grid-invariant in relative terms.
    Matches the Hard Blocker #3 ``_boundary_tol_um`` formula
    exactly so the two modules share identical inclusion behavior."""

    nx, ny = ecm.grid_shape
    dx = float(ecm.spacing_um)
    return _BOUNDARY_TOL_RELATIVE * max(nx * dx, ny * dx)


def _validate_fa_position_and_schema(
    fa: FocalAdhesionState, x_min: float, x_max: float, y_min: float, y_max: float
) -> tuple[float, float]:
    """Run the locked validation order: structural len-check on
    ``position_um_xy`` first (defer to schema validator for the
    canonical error message — Codex review id=1385), then
    lock-specific ``non_finite_fa_position`` failure_kind, then
    per-FA ``fa.validate()`` for other schema invariants
    (Codex review id=1378), then out-of-grid inclusion. Returns
    the validated ``(x_fa, y_fa)`` tuple."""

    # Defensive structural len check before unpacking. If
    # `fa.position_um_xy` has length != 2 (e.g., a 1-tuple), the
    # subsequent unpack would raise Python's "not enough values
    # to unpack" ValueError instead of the schema's canonical
    # "position_um_xy must contain exactly 2 finite values"
    # message. Defer to fa.validate() in that case so the schema
    # error surfaces with its locked wording.
    if len(fa.position_um_xy) != 2:
        fa.validate()
        # Unreachable: fa.validate() must raise on len != 2.
        raise FAToECMBiasError(
            "non_finite_fa_position",
            f"FA {fa.adhesion_id!r} position_um_xy length != 2 "
            f"and FocalAdhesionState.validate() did not raise; "
            f"this should be unreachable",
        )

    x_fa, y_fa = fa.position_um_xy
    x_fa = float(x_fa)
    y_fa = float(y_fa)
    if not (math.isfinite(x_fa) and math.isfinite(y_fa)):
        raise FAToECMBiasError(
            "non_finite_fa_position",
            f"FA {fa.adhesion_id!r} position_um_xy={fa.position_um_xy!r} "
            f"is not finite",
        )
    # Per-FA schema validation per the Hard Blocker #3 fix
    # (`883efc8`) and the local dynamics precedent. Schema-only
    # invariants such as "no traction without attachment"
    # (`state == "unbound"` with non-zero traction) propagate as
    # the schema's ``ValueError`` after the locked failure_kind
    # has already had its chance to fire.
    fa.validate()
    if not (x_min <= x_fa <= x_max) or not (y_min <= y_fa <= y_max):
        raise FAToECMBiasError(
            "fa_bias_position_outside_ecm_grid",
            f"FA {fa.adhesion_id!r} at ({x_fa!r}, {y_fa!r}) is outside "
            f"ECM physical footprint",
        )
    return x_fa, y_fa


def _bilinear_weights_and_indices(
    x_fa: float,
    y_fa: float,
    origin_x: float,
    origin_y: float,
    dx: float,
    nx: int,
    ny: int,
) -> tuple[list[tuple[int, int, float]], float]:
    """Compute the 4-point bilinear stencil indices and weights for
    an FA at ``(x_fa, y_fa)`` on a cell-centered grid; filter to
    in-grid centers; return the in-grid `(k, l, weight)` list and
    the raw weight sum (before renormalization). Matches the Hard
    Blocker #3 scatter geometry exactly so a sampled point and a
    scattered point at the same FA position use identical
    interpolation weights."""

    u = (x_fa - origin_x) / dx - 0.5
    v = (y_fa - origin_y) / dx - 0.5
    i_lo = int(math.floor(u))
    j_lo = int(math.floor(v))
    fx = u - i_lo
    fy = v - j_lo
    candidates = (
        (i_lo, j_lo, (1.0 - fx) * (1.0 - fy)),
        (i_lo + 1, j_lo, fx * (1.0 - fy)),
        (i_lo, j_lo + 1, (1.0 - fx) * fy),
        (i_lo + 1, j_lo + 1, fx * fy),
    )
    valid = [
        (k, l, w) for (k, l, w) in candidates if 0 <= k < nx and 0 <= l < ny
    ]
    weight_sum = sum(w for (_k, _l, w) in valid)
    return valid, weight_sum


def sample_ecm_at_fa_positions(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
) -> ECMSampledAtFAs:
    """Bilinear interpolation read of ECM fields at FA positions.

    Convention matches
    :func:`acs.v2.dynamics.fa_to_ecm_scattering.scatter_fa_traction_to_ecm_bilinear`:
    cell-centered grid (``x_i = origin_x + (i + 0.5) * dx``),
    inclusive physical footprint, boundary half-cell uses
    truncated stencil with renormalized weights, out-of-grid
    raises ``FAToECMBiasError(failure_kind=
    "fa_bias_position_outside_ecm_grid")`` (distinct from Hard
    Blocker #3's ``fa_position_outside_ecm_grid`` because the
    failing operation is read/bias, not scatter; physical
    footprint rule identical).

    Pure read function: ``ecm`` is not mutated, no FA state is
    mutated, no scalarization is performed (orientation tensor
    returned full ``(N_FA, 2, 2)``).

    Empty FA list returns an :class:`ECMSampledAtFAs` with all
    arrays of shape ``(0, ...)`` and ``fa_ids = ()``.

    Args:
        adhesions: tuple or list of :class:`FocalAdhesionState`.
            The locked design accepts a tuple; list is accepted
            as a runtime convenience.
        ecm: :class:`ECMSubstrateState`. Validated at entry.

    Returns:
        :class:`ECMSampledAtFAs` with input order preserved.
    """

    ecm.validate()

    nx, ny = ecm.grid_shape
    dx = float(ecm.spacing_um)
    origin_x, origin_y = ecm.origin_um_xy
    origin_x = float(origin_x)
    origin_y = float(origin_y)
    tol = _boundary_tol_um(ecm)
    x_min = origin_x - tol
    x_max = origin_x + nx * dx + tol
    y_min = origin_y - tol
    y_max = origin_y + ny * dx + tol

    n_fa = len(adhesions)
    fa_ids_list: list[str] = []
    stiffness = np.zeros((n_fa,), dtype=np.float64)
    fiber = np.zeros((n_fa,), dtype=np.float64)
    ligand = np.zeros((n_fa,), dtype=np.float64)
    orientation = np.zeros((n_fa, 2, 2), dtype=np.float64)
    accumulated = np.zeros((n_fa,), dtype=np.float64)

    if n_fa == 0:
        return ECMSampledAtFAs(
            fa_ids=(),
            stiffness_kpa=stiffness,
            fiber_density=fiber,
            ligand_density=ligand,
            orientation_tensor=orientation,
            accumulated_traction_nNs_per_um2=accumulated,
        )

    stiffness_kpa = np.asarray(ecm.stiffness_kpa, dtype=np.float64)
    fiber_density = np.asarray(ecm.fiber_density, dtype=np.float64)
    ligand_density = np.asarray(ecm.ligand_density, dtype=np.float64)
    orientation_tensor = np.asarray(ecm.orientation_tensor, dtype=np.float64)
    accumulated_traction = np.asarray(
        ecm.accumulated_traction_nNs_per_um2, dtype=np.float64
    )

    for i, fa in enumerate(adhesions):
        x_fa, y_fa = _validate_fa_position_and_schema(
            fa, x_min, x_max, y_min, y_max
        )
        valid, weight_sum = _bilinear_weights_and_indices(
            x_fa, y_fa, origin_x, origin_y, dx, nx, ny
        )
        if not valid:
            # Defensive — should be unreachable inside the
            # inclusive footprint plus tolerance.
            raise FAToECMBiasError(
                "fa_bias_position_outside_ecm_grid",
                f"FA {fa.adhesion_id!r} bilinear stencil has no in-grid "
                f"centers; this should be unreachable within the "
                f"inclusive footprint plus tolerance",
            )
        for k, l, w in valid:
            w_norm = w / weight_sum
            stiffness[i] += w_norm * stiffness_kpa[k, l]
            fiber[i] += w_norm * fiber_density[k, l]
            ligand[i] += w_norm * ligand_density[k, l]
            orientation[i] += w_norm * orientation_tensor[k, l]
            accumulated[i] += w_norm * accumulated_traction[k, l]
        fa_ids_list.append(fa.adhesion_id)

    return ECMSampledAtFAs(
        fa_ids=tuple(fa_ids_list),
        stiffness_kpa=stiffness,
        fiber_density=fiber,
        ligand_density=ligand,
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=accumulated,
    )


def compute_ecm_to_fa_bias_neutral(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
) -> ECMToFABiasResult:
    """Phase D no-op default. Returns 1.0 multipliers for all FAs
    and all rate names in :data:`RATE_NAMES` order.

    Includes sampled diagnostics from
    :func:`sample_ecm_at_fa_positions` for audit/test only.
    Sampled values do **not** alter multipliers in Phase D
    (effective-stiffness side-door guard, locked phased plan §3).

    Phase E active mapping (banded / continuous / anisotropic) is
    a SEPARATE function ``compute_ecm_to_fa_bias_active(...)``
    with its own Sanity Gate and lock. The structural separation
    is the silent-activation guard: a caller cannot trigger
    active behavior through this Phase D function.

    Args:
        adhesions: tuple or list of :class:`FocalAdhesionState`.
        ecm: :class:`ECMSubstrateState`. Validated at entry.

    Returns:
        :class:`ECMToFABiasResult` with all-1.0 multipliers
        (shape ``(N_FA, 3)``) and sampled diagnostics.
    """

    sampled = sample_ecm_at_fa_positions(adhesions, ecm)
    n_fa = len(sampled.fa_ids)
    multipliers = np.ones((n_fa, len(RATE_NAMES)), dtype=np.float64)
    if n_fa == 0:
        diagnostics_dict: dict[str, float | int | str] = {
            "n_adhesions": 0,
            "max_multiplier": 1.0,
            "min_multiplier": 1.0,
            "sampler_geometry": "bilinear_cell_centered",
        }
    else:
        diagnostics_dict = {
            "n_adhesions": n_fa,
            "max_multiplier": float(multipliers.max()),
            "min_multiplier": float(multipliers.min()),
            "sampler_geometry": "bilinear_cell_centered",
        }
    return ECMToFABiasResult(
        fa_ids=sampled.fa_ids,
        multipliers_per_fa=multipliers,
        rate_names=tuple(RATE_NAMES),
        sampled_diagnostics=sampled,
        diagnostics_dict=diagnostics_dict,
    )


_DEVIATORIC_SCORE_BOUND: Final[float] = float(np.sqrt(2.0))


def _validate_k_active(k_active: float) -> None:
    if isinstance(k_active, bool):
        raise FAToECMBiasError(
            "k_active_invalid",
            f"k_active must be float (not bool — Python bool subset int trap), "
            f"got {type(k_active).__name__}",
        )
    if not isinstance(k_active, (int, float, np.floating, np.integer)):
        raise FAToECMBiasError(
            "k_active_invalid",
            f"k_active must be a numeric scalar (int / float / numpy scalar), "
            f"got {type(k_active).__name__}",
        )
    if not (np.isfinite(k_active) and k_active > 0.0):
        raise FAToECMBiasError(
            "k_active_invalid",
            f"k_active must be finite positive, got {k_active!r}",
        )
    derived_bound = float(np.exp(k_active * _DEVIATORIC_SCORE_BOUND))
    if not np.isfinite(derived_bound):
        raise FAToECMBiasError(
            "k_active_invalid",
            f"k_active={k_active!r} produces non-finite multiplier bound "
            f"exp({k_active!r}*sqrt(2)) = {derived_bound!r}; reduce k_active",
        )


def compute_ecm_to_fa_bias_active(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
    *,
    k_active: float,
) -> ECMToFABiasResult:
    """HB#4-active step 1: orientation-driven Rayleigh-quotient bias on
    the deviatoric ECM orientation tensor.

    For each FA i:
        n_fa = traction / |traction|  (branch-defined zero-traction → exact neutral)
        T_local = ECM orientation tensor sampled at FA position
        Q = T_local - 0.5 * trace(T_local) * I  (deviatoric: alignment-only)
        score = n_fa.T @ Q @ n_fa  (dimensionless; |score| <= sqrt(2) under |T_ij| <= 1)
        multiplier = exp(k_active * score)  (bounded positive exponential map)

    All three FA rates (k_maturity_per_s, k_bind_per_s, k_unbind_per_s)
    receive the SAME multiplier per FA — active-bias scaffolding only;
    per-rate selectivity is Phase E v3+ scope.

    Range: multiplier in [exp(-sqrt(2)*k_active), exp(+sqrt(2)*k_active)];
    multiplier == 1.0 exactly at score == 0 (isotropic T = a*I gives Q = 0
    so any traction direction yields neutral; Item 5 IC T = 0.5*I is
    therefore automatically neutral under this law).

    Validation order (ecm -> k_active -> sampler) lets cheap parameter
    failures fire before the expensive sampler is called even with empty
    adhesions.

    Reads only ``orientation_tensor`` from the sampled ECM; the other
    sampled ECM fields (mechanosensing scalars and accumulated-traction
    memory) are excluded by the function-scoped AST + string guard
    in the test catalog (HB#5 Y11 sister-pattern).

    Args:
        adhesions: tuple or list of :class:`FocalAdhesionState`.
        ecm: :class:`ECMSubstrateState`. Validated at entry.
        k_active: dimensionless coupling strength. REQUIRED, no
            default. Must be finite, positive, and produce a finite
            ``exp(k_active*sqrt(2))`` upper bound.

    Returns:
        :class:`ECMToFABiasResult` with deviatoric-Rayleigh
        multipliers (shape ``(N_FA, 3)``, all three rates same per
        FA) and 10-key diagnostics dict including the alignment
        score range and the score bound value.

    Raises:
        FAToECMBiasError: with ``failure_kind == "k_active_invalid"``
            if ``k_active`` is invalid; ``failure_kind ==
            "active_multiplier_non_finite"`` if a per-FA score or
            multiplier becomes non-finite (defense in depth);
            sampler failure_kinds propagate from
            :func:`sample_ecm_at_fa_positions`.
    """

    ecm.validate()
    _validate_k_active(k_active)
    sampled = sample_ecm_at_fa_positions(adhesions, ecm)

    n_fa = len(sampled.fa_ids)
    n_rates = len(RATE_NAMES)
    multipliers = np.ones((n_fa, n_rates), dtype=np.float64)
    alignment_scores = np.zeros(n_fa, dtype=np.float64)
    zero_traction_count = 0
    identity_2x2 = np.eye(2, dtype=np.float64)

    for i, fa in enumerate(adhesions):
        traction = np.asarray(fa.traction_force_nN_xy, dtype=np.float64)
        traction_norm = float(np.linalg.norm(traction))
        if traction_norm > 0.0:
            n_fa_unit = traction / traction_norm
            t_local = sampled.orientation_tensor[i]
            trace_t = float(t_local[0, 0] + t_local[1, 1])
            q_dev = t_local - 0.5 * trace_t * identity_2x2
            alignment_score = float(n_fa_unit @ q_dev @ n_fa_unit)
        else:
            alignment_score = 0.0
            zero_traction_count += 1

        if not np.isfinite(alignment_score):
            raise FAToECMBiasError(
                "active_multiplier_non_finite",
                f"FA {fa.adhesion_id!r}: alignment_score {alignment_score!r} "
                f"non-finite (suggests corrupted orientation_tensor; "
                f"ecm.validate() should have caught)",
            )
        multiplier = float(np.exp(k_active * alignment_score))
        if not np.isfinite(multiplier):
            raise FAToECMBiasError(
                "active_multiplier_non_finite",
                f"FA {fa.adhesion_id!r}: multiplier {multiplier!r} non-finite "
                f"(k_active={k_active!r}, score={alignment_score!r})",
            )

        multipliers[i, :] = multiplier
        alignment_scores[i] = alignment_score

    if n_fa > 0:
        max_multiplier = float(multipliers.max())
        min_multiplier = float(multipliers.min())
        max_score = float(alignment_scores.max())
        min_score = float(alignment_scores.min())
    else:
        max_multiplier = 1.0
        min_multiplier = 1.0
        max_score = 0.0
        min_score = 0.0

    diagnostics_dict: dict[str, float | int | str] = {
        "n_adhesions": int(n_fa),
        "max_multiplier": max_multiplier,
        "min_multiplier": min_multiplier,
        "sampler_geometry": "bilinear_cell_centered",
        "mechanism": "deviatoric_rayleigh_orientation",
        "k_active": float(k_active),
        "max_alignment_score": max_score,
        "min_alignment_score": min_score,
        "zero_traction_count": int(zero_traction_count),
        "score_bound": float(_DEVIATORIC_SCORE_BOUND),
    }

    return ECMToFABiasResult(
        fa_ids=sampled.fa_ids,
        multipliers_per_fa=multipliers,
        rate_names=tuple(RATE_NAMES),
        sampled_diagnostics=sampled,
        diagnostics_dict=diagnostics_dict,
    )
