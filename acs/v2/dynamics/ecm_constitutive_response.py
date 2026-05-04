"""ECM orientation constitutive response (HB#1+#2 merged lock).

Implements the locked HB#1+#2 design from
``docs/v2_hard_blocker_1_2_constitutive_direction_locked.md``
(commits ``7b1d3cf`` initial + ``c18a5dc`` ``ecm.validate()`` at
entry + ``42d34e7`` error class + 7 exports + ``7154342``
stale-count cleanup) and the impl-work Sanity Gate
``docs/v2_hard_blocker_1_2_constitutive_response_sanity_gate.md``
(commit ``186fb75`` + same amendments; Codex Sanity Gate PASS at
``id=1492``).

This is the **first Phase E active constitutive law**: under
instantaneous FA traction stimulus, the ECM ``orientation_tensor``
field updates toward a per-cell rank-1 target ``n⊗n`` via
exact-exponential convex combination. Saturation is built into
the convex weight ``w ∈ [0, 1]``; no separate hard cap is needed.

The module deliberately does **not**:

- update ``stiffness_kpa``, ``fiber_density``, ``ligand_density``,
  or ``accumulated_traction_nNs_per_um2`` (orientation-only first
  law per locked Y1);
- consume ``accumulated_traction_nNs_per_um2`` as driver
  (instantaneous traction only per locked Y6);
- use any ε / epsilon as physical or numerical regularizer
  (branch-defined zero-traction only per locked Y4);
- apply hard caps or clamps on ``orientation_tensor`` components
  (saturation built into convex weight per locked Y10);
- claim PSD / eigenvalue / trace properties beyond the schema's
  componentwise ``|T_ij| ≤ 1`` invariant (locked Y5);
- alias unchanged arrays into the returned ECM (all 5 array
  fields fresh-copied per locked Y12);
- introduce ``Generic[...]`` parameterization or preemptive
  ``# type: ignore`` (B1 sister-precedent);
- silently no-op on ``k_orient_per_s == 0`` (raises strict per
  locked Y11);
- expose any mechanosensing shortcut (latent-orientation contract
  per locked Y7: future downstream consumers using
  ``orientation_tensor`` MUST decide whether/how to weight by
  ``fiber_density``).

Two-sided validation contract (Codex ``id=1486`` silent-heal-path
closure):

- ``ecm.validate()`` at function entry — invalid input ECM is
  REJECTED, not "healed" by the convex update toward
  ``T_target`` under high traction × dt. Mirrors HB#3 / HB#4
  local dynamics precedent.
- ``updated_ecm.validate()`` before return — catches any
  float-rounding drift in output.

Failure-kind discipline (Codex ``id=1488`` API-surface lock):
``ECMConstitutiveResponseError(ValueError)`` with 5 owned
``failure_kind`` values: ``dt_invalid``,
``traction_ref_invalid``, ``k_orient_invalid``,
``traction_shape_mismatch``, ``non_finite_traction``. Inherited
``ecm.validate()`` failures stay plain schema ``ValueError`` —
schema-level errors stay schema-level, not wrapped.

Sanity Gate scope (acs/v2/dynamics/ecm_constitutive_response.py):

- §1 dimensional: ``traction_norm [nN/μm²] / TRACTION_REF
  [nN/μm²] = S [-]`` dimensionless; ``K [1/s] · S [-] · dt [s]
  = [-]`` dimensionless argument to ``np.expm1``. Full Hard Rule
  10 unit chain inline-derived in the locked §1 +
  Sanity-Gate §1.
- §2 boundary: empty traction / ``dt_s == 0`` / ``dt_s < 0`` or
  ``bool`` or nonfinite / strict-positive params /
  shape mismatch / asymptotic limits — all locked + tested.
- §3 conservation: convex algebra preserves schema invariant
  under valid input; invalid input rejected at
  ``ecm.validate()``; unchanged-fields bytewise-equal +
  no-aliasing per locked Y12; per-step distance-to-target
  monotone decreasing under fixed target per locked Y2.
- §4 numerical: ``-np.expm1(-K·S·dt)`` for stability near zero
  argument per locked Y10; ``np.float64`` throughout; no new
  tolerance; unconditionally stable in ``dt`` via convex-weight
  construction.
- §5 sign: convex weight ``∈ [0, 1]`` by construction;
  distance-to-target monotone decreasing per step under fixed
  target.
- §6 measurement-protocol: strictly local per-cell algebra (no
  off-proof points like Stage 1a v13); Rule 10 bridge fields in
  typed diagnostics (Y13); two-sided validation contract
  enforced via test 19 static check + test 20 invalid-input
  rejection (Codex ``id=1486`` + ``id=1488``).

Magic-Number Block: 3 module constants — all derivable from
literature (Munevar 2001 + Tan 2003 for ``TRACTION_REF``; Hall
2016 + Notbohm 2015 for ``TAU_ALIGN_RANGE_S``;
``K_ORIENT_PER_S`` is the geometric midpoint of the locked
range per locked Y9), grid-invariant in relative terms (no
``dx`` / ``dt`` / ``grid_n`` dependence), and chosen before
any A/A₀ comparison or test target gating.

Phase E activation remains BLOCKED on HB#5 (Lyapunov metric) +
the effective_stiffness law decision. This module landing does
NOT authorize Phase E composition with the rest of the
closed-loop ECM gate; it provides the constitutive-law primitive
that a future Phase E composition will call.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np

from acs.v2.ecm_substrate import ECMSubstrateState

TRACTION_REF_NN_PER_UM2: Final[float] = 1.0
TAU_ALIGN_RANGE_S: Final[tuple[float, float]] = (3600.0, 7200.0)
K_ORIENT_PER_S: Final[float] = float(
    1.0 / math.sqrt(TAU_ALIGN_RANGE_S[0] * TAU_ALIGN_RANGE_S[1])
)


class ECMConstitutiveResponseError(ValueError):
    """Validation error raised by
    :func:`step_ecm_orientation_response` when its inputs violate
    the locked HB#1+#2 contract. Carries a machine-readable
    :attr:`failure_kind` attribute (HB#3 / HB#4 sister-pattern).

    Owned ``failure_kind`` values:

    - ``"dt_invalid"`` — ``dt_s`` is negative, of type ``bool``,
      or non-finite.
    - ``"traction_ref_invalid"`` — ``traction_ref_nN_per_um2`` is
      non-positive or non-finite.
    - ``"k_orient_invalid"`` — ``k_orient_per_s`` is non-positive
      or non-finite (strict ``> 0``; no zero no-op silent path).
    - ``"traction_shape_mismatch"`` — ``traction_density_xy`` does
      not have shape ``(*ecm.grid_shape, 2)``.
    - ``"non_finite_traction"`` — ``traction_density_xy`` contains
      ``nan`` or ``±inf``.

    Inherited ``ecm.validate()`` failures from the input ECM stay
    plain schema ``ValueError`` and are NOT wrapped in this
    error type (Codex ``id=1488`` — schema-level errors stay
    schema-level).
    """

    def __init__(self, failure_kind: str, message: str) -> None:
        self.failure_kind = failure_kind
        super().__init__(message)


@dataclass(frozen=True, slots=True)
class ECMOrientationResponseDiagnostics:
    """Typed diagnostics for one HB#1+#2 ECM orientation response step.

    13 fields covering cell-count, input-parameter traceback for
    sweep reproducibility, Hard Rule 10 bridge fields,
    convex-weight statistics, distance-to-target statistics
    (``_nonzero`` suffixed to avoid zero-traction placeholder
    confusion), and schema-sanity residuals.

    Attributes:
        n_nonzero_cells: count of grid cells with non-zero
            traction stimulus this step (count, dimensionless).
        n_total_cells: total grid cell count (count,
            dimensionless).
        traction_ref_nN_per_um2: caller-passed ``traction_ref``
            value, recorded for sweep traceback (``nN/μm²``).
        k_orient_per_s: caller-passed ``k_orient`` value
            (``1/s``).
        max_traction_norm_nN_per_um2: max L2 norm of per-cell
            traction across the grid (``nN/μm²``); Rule 10
            bridge to weight.
        max_dimensionless_stimulus: max of `S` per cell
            (dimensionless); Rule 10 bridge.
        max_convex_weight: max of `w` per cell (dimensionless,
            ``∈ [0, 1]``); saturation indicator.
        mean_convex_weight_nonzero: mean of `w` over nonzero
            cells (dimensionless); ``0.0`` when no nonzero
            cells.
        max_distance_to_target_frobenius_nonzero: max Frobenius
            norm of `T_new - T_target` over nonzero cells
            (dimensionless); ``0.0`` when no nonzero cells.
        mean_distance_to_target_frobenius_nonzero: mean
            Frobenius norm of `T_new - T_target` over nonzero
            cells; ``0.0`` when no nonzero cells.
        max_orientation_delta_frobenius: max Frobenius norm of
            `T_new - T_old` per cell across the whole grid.
        symmetry_residual_max: max of
            ``|T_new - T_new.T|`` componentwise; schema sanity
            residual (always ≤ float-rounding tolerance under
            valid input).
        componentwise_bound_residual_max: max of
            ``max(0, |T_new_ij| - 1)`` componentwise; schema
            sanity residual (always ≤ float-rounding tolerance
            under valid input).
    """

    n_nonzero_cells: int
    n_total_cells: int
    traction_ref_nN_per_um2: float
    k_orient_per_s: float
    max_traction_norm_nN_per_um2: float
    max_dimensionless_stimulus: float
    max_convex_weight: float
    mean_convex_weight_nonzero: float
    max_distance_to_target_frobenius_nonzero: float
    mean_distance_to_target_frobenius_nonzero: float
    max_orientation_delta_frobenius: float
    symmetry_residual_max: float
    componentwise_bound_residual_max: float


@dataclass(frozen=True, slots=True)
class ECMOrientationResponseResult:
    """Outputs of one HB#1+#2 ECM orientation response step."""

    updated_ecm: ECMSubstrateState
    diagnostics: ECMOrientationResponseDiagnostics


def _validate_dt_nonneg_finite_nonbool(dt_s: float) -> float:
    if isinstance(dt_s, bool):
        raise ECMConstitutiveResponseError(
            "dt_invalid",
            f"dt_s must be a numeric float, not bool; got {dt_s!r}",
        )
    if not isinstance(dt_s, (int, float)):
        raise ECMConstitutiveResponseError(
            "dt_invalid",
            f"dt_s must be a numeric float; got {type(dt_s).__name__}",
        )
    value = float(dt_s)
    if not math.isfinite(value):
        raise ECMConstitutiveResponseError(
            "dt_invalid", f"dt_s must be finite; got {value!r}"
        )
    if value < 0.0:
        raise ECMConstitutiveResponseError(
            "dt_invalid", f"dt_s must be >= 0.0; got {value!r}"
        )
    return value


def _validate_positive_finite(value: float, name: str) -> float:
    if isinstance(value, bool):
        raise ECMConstitutiveResponseError(
            f"{name}_invalid",
            f"{name} must be a numeric float, not bool; got {value!r}",
        )
    if not isinstance(value, (int, float)):
        raise ECMConstitutiveResponseError(
            f"{name}_invalid",
            f"{name} must be a numeric float; got {type(value).__name__}",
        )
    v = float(value)
    if not math.isfinite(v):
        raise ECMConstitutiveResponseError(
            f"{name}_invalid", f"{name} must be finite; got {v!r}"
        )
    if v <= 0.0:
        raise ECMConstitutiveResponseError(
            f"{name}_invalid",
            f"{name} must be strictly positive; got {v!r} (no zero no-op)",
        )
    return v


def _validate_finite_traction(traction_density_xy: np.ndarray) -> np.ndarray:
    arr = np.asarray(traction_density_xy, dtype=np.float64)
    if not np.isfinite(arr).all():
        raise ECMConstitutiveResponseError(
            "non_finite_traction",
            "traction_density_xy contains non-finite values "
            "(nan / ±inf)",
        )
    return arr


def _validate_traction_shape(
    traction_density_xy: np.ndarray, grid_shape: tuple[int, int]
) -> None:
    expected = (*grid_shape, 2)
    if traction_density_xy.shape != expected:
        raise ECMConstitutiveResponseError(
            "traction_shape_mismatch",
            f"traction_density_xy must have shape {expected!r}; got "
            f"{traction_density_xy.shape!r}",
        )


def step_ecm_orientation_response(
    ecm: ECMSubstrateState,
    traction_density_xy: np.ndarray,
    dt_s: float,
    *,
    traction_ref_nN_per_um2: float = TRACTION_REF_NN_PER_UM2,
    k_orient_per_s: float = K_ORIENT_PER_S,
) -> ECMOrientationResponseResult:
    """One step of the HB#1+#2 ECM orientation constitutive response.

    Computes a per-cell convex update of ``ecm.orientation_tensor``
    toward the rank-1 target ``n_stim ⊗ n_stim`` driven by
    instantaneous traction stimulus ``traction_density_xy`` (HB#3
    output). The update is exact-exponential
    (``-np.expm1(-K · S · dt)``); saturation is built into the
    convex weight ``w ∈ [0, 1]``.

    Two-sided validation contract per Codex ``id=1486``:
    ``ecm.validate()`` rejects invalid input ECM at entry;
    ``updated_ecm.validate()`` catches float-rounding drift in
    output.

    Failure modes (raised as
    :class:`ECMConstitutiveResponseError` with the listed
    ``failure_kind``):

    - ``dt_invalid``: ``dt_s`` is negative / ``bool`` / non-finite.
    - ``traction_ref_invalid``: ``traction_ref_nN_per_um2`` is
      non-positive / non-finite / ``bool``.
    - ``k_orient_invalid``: ``k_orient_per_s`` is non-positive /
      non-finite / ``bool``.
    - ``traction_shape_mismatch``: ``traction_density_xy`` shape
      does not match ``(*ecm.grid_shape, 2)``.
    - ``non_finite_traction``: ``traction_density_xy`` contains
      ``nan`` or ``±inf``.

    Inherited ``ecm.validate()`` failures (e.g., invalid
    ``orientation_tensor``) remain plain schema ``ValueError``
    and are NOT wrapped.

    Args:
        ecm: input :class:`ECMSubstrateState`. Validated at entry.
        traction_density_xy: ``(*ecm.grid_shape, 2)`` float64
            array in ``nN/μm²``, the per-cell vector traction
            density (HB#3 output).
        dt_s: timestep in seconds. ``0.0`` is valid no-op
            (``w == 0`` everywhere).
        traction_ref_nN_per_um2: literature-pinned reference
            traction scale (default ``TRACTION_REF_NN_PER_UM2``,
            ``1.0 nN/μm² = 1.0 kPa`` from Munevar 2001 + Tan
            2003).
        k_orient_per_s: alignment rate constant (default
            ``K_ORIENT_PER_S`` ≈ ``1.964e-4 /s``, geometric
            midpoint of ``TAU_ALIGN_RANGE_S`` from Hall 2016 +
            Notbohm 2015).

    Returns:
        :class:`ECMOrientationResponseResult` with the updated
        ECM (5 fresh-copied array fields, no aliasing per locked
        Y12) and the typed 13-field diagnostics.
    """

    ecm.validate()
    arr = _validate_finite_traction(traction_density_xy)
    dt_value = _validate_dt_nonneg_finite_nonbool(dt_s)
    traction_ref_value = _validate_positive_finite(
        traction_ref_nN_per_um2, "traction_ref"
    )
    k_orient_value = _validate_positive_finite(k_orient_per_s, "k_orient")
    _validate_traction_shape(arr, ecm.grid_shape)

    traction_norm = np.linalg.norm(arr, axis=-1)
    nonzero = traction_norm > 0.0
    S = traction_norm / traction_ref_value

    n_xy = np.zeros_like(arr, dtype=np.float64)
    n_xy[nonzero] = arr[nonzero] / traction_norm[nonzero, None]
    T_target = np.einsum("...i,...j->...ij", n_xy, n_xy)

    w = -np.expm1(-k_orient_value * S * dt_value)

    orientation_old = np.asarray(ecm.orientation_tensor, dtype=np.float64)
    orientation_new = (
        (1.0 - w[..., None, None]) * orientation_old
        + w[..., None, None] * T_target
    )

    updated_ecm = ECMSubstrateState(
        origin_um_xy=tuple(ecm.origin_um_xy),
        spacing_um=float(ecm.spacing_um),
        stiffness_kpa=np.array(ecm.stiffness_kpa, dtype=np.float64, copy=True),
        ligand_density=np.array(
            ecm.ligand_density, dtype=np.float64, copy=True
        ),
        fiber_density=np.array(ecm.fiber_density, dtype=np.float64, copy=True),
        orientation_tensor=np.array(
            orientation_new, dtype=np.float64, copy=True
        ),
        accumulated_traction_nNs_per_um2=np.array(
            ecm.accumulated_traction_nNs_per_um2,
            dtype=np.float64,
            copy=True,
        ),
        source=ecm.source,
    )
    updated_ecm.validate()

    diff = orientation_new - T_target
    delta = orientation_new - orientation_old
    if nonzero.any():
        mean_w_nz = float(w[nonzero].mean())
        diff_nz_norms = np.linalg.norm(diff[nonzero], axis=(-2, -1))
        max_dist = float(diff_nz_norms.max())
        mean_dist = float(diff_nz_norms.mean())
    else:
        mean_w_nz = 0.0
        max_dist = 0.0
        mean_dist = 0.0
    diagnostics = ECMOrientationResponseDiagnostics(
        n_nonzero_cells=int(nonzero.sum()),
        n_total_cells=int(np.prod(ecm.grid_shape)),
        traction_ref_nN_per_um2=traction_ref_value,
        k_orient_per_s=k_orient_value,
        max_traction_norm_nN_per_um2=float(traction_norm.max()),
        max_dimensionless_stimulus=float(S.max()),
        max_convex_weight=float(w.max()),
        mean_convex_weight_nonzero=mean_w_nz,
        max_distance_to_target_frobenius_nonzero=max_dist,
        mean_distance_to_target_frobenius_nonzero=mean_dist,
        max_orientation_delta_frobenius=float(
            np.linalg.norm(delta, axis=(-2, -1)).max()
        ),
        symmetry_residual_max=float(
            np.abs(orientation_new - np.swapaxes(orientation_new, -1, -2)).max()
        ),
        componentwise_bound_residual_max=float(
            max(0.0, np.abs(orientation_new).max() - 1.0)
        ),
    )

    return ECMOrientationResponseResult(
        updated_ecm=updated_ecm,
        diagnostics=diagnostics,
    )
