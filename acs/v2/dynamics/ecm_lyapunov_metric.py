"""Pure read-only ECM orientation Lyapunov-like metric (HB#5 lock).

Implements the locked HB#5 design from
``docs/v2/v2_hard_blocker_5_lyapunov_metric_locked.md`` (commit
``340c353``) and the impl-work Sanity Gate
``docs/v2/v2_hard_blocker_5_lyapunov_metric_sanity_gate.md`` (commit
``1c9bb0a``; Codex Sanity Gate PASS at ``id=1518``).

This is the **final Phase E upstream hard blocker**: a pure
read-only metric that provides the ECM-orientation stability
evidence channel for the closed-loop ECM gate. Together with
HB#3 (FA→ECM scattering), HB#4 (ECM→FA bias neutral), HB#1+#2
(active orientation update), and Phase D (no-op orchestrator),
it completes the Phase E upstream blocker set; the next blocker
is the closed-loop ECM gate composition cycle (separate unit).

The module deliberately does **not**:

- include zero-traction cells with ``T_target = 0`` placeholder
  in the V sum (Y1 — orientation-magnitude-penalty smuggling
  via placeholder; Hard Rule 11 measurement-protocol
  consistency: metric domain MUST match HB#1+#2 update domain
  which is active cells only);
- read or reference ``accumulated_traction_nNs_per_um2``,
  ``stiffness_kpa``, ``fiber_density``, or ``ligand_density``
  (Y7 — metric is for orientation-only; AST-walk meta-test
  enforces this);
- mutate any ECM array field (function is pure read-only; no
  copies needed because no return-value ECM);
- offer a pluggable strategy / composite tuple return / dict
  diagnostics (Y5 — single typed result + typed diagnostics,
  B1 sister-pattern);
- use any ε / epsilon as physical or numerical regularizer in
  the active mask (Y7 — branch-defined ``traction_norm > 0.0``
  only);
- apply any hard cap or clamp on ``v_active_um2`` (the bound
  is a diagnostic field for testing, not an in-line clamp);
- claim "positive-definite over full grid" (Y2 — actual
  property is positive-semidefinite relative to the
  active-input measurement protocol with equilibrium set
  including all-passive states);
- conflate "cumulative bound" with "Lyapunov-like decrease"
  (Y3 — boundedness and per-step contraction are separate
  evidence channels and must be tested separately);
- claim "full closed-loop FA→ECM→FA stability satisfied by
  HB#5" (Y4 — HB#5 is ECM-response-half evidence only; full
  loop requires HB#4-active design which is currently HB#4
  Phase D no-op neutral multipliers);
- introduce an error class with ``failure_kind`` attribute
  (HB#5 has NO domain-specific failure modes; all input
  validation failures raise plain ``ValueError`` per locked
  §5 — different from HB#1+#2's
  ``ECMConstitutiveResponseError`` because the metric is a
  pure read-only computation, not a constitutive law).

Validation contract per HB#1+#2 ``id=1486`` sister-precedent:
``ecm.validate()`` at function entry rejects invalid input ECM
BEFORE any algebra runs (silent-heal-path closure analog).
Plain schema ``ValueError`` from ``ecm.validate()`` is NOT
wrapped.

Sanity Gate scope (acs/v2/dynamics/ecm_lyapunov_metric.py):

- §1 dimensional: ``traction_density_xy [nN/μm²]``, ``cell_area
  [μm²]``, ``||T - n⊗n||²_F`` dimensionless ⇒ ``V_active_um2
  [μm²]``. Full Hard Rule 10 unit chain inline-derived in the
  locked §1 + Sanity Gate §1.
- §2 boundary: 7 boundary classes (zero/all-active/at-target/
  adversarial-saturation/single-cell/invalid-ECM/non-finite-or-
  shape) all locked + 6 of 7 test-covered.
- §3 conservation: pure read-only; two evidence channels
  (fixed-target per-step contraction Y3-A + varying-target
  boundedness Y3-B) tested separately; locked forbidden item
  rejects conflation.
- §4 numerical: ``np.float64`` throughout; ``np.einsum`` outer
  product + ``np.sum(diff**2)`` Frobenius² numerically stable
  for the 2x2-tensor + grid-size regime; no new tolerance
  (Y7 ε-drop).
- §5 sign: ``V_active ≥ 0`` always (sum of squares times
  nonneg ``0.5 · cell_area``); 10 diagnostics fields all
  non-negative by construction.
- §6 measurement-protocol (Hard Rule 11, central anchor):
  domain matching active cells only (Y1); adversarial bound
  9.0 derived analytically + verified via Frobenius identity
  ``||T - n⊗n||²_F = ||T||²_F + 1 − 2 n^T T n``; two-constant
  9 → 4.5 algebraic chain visible in code (Y8); two-channel
  evidence separation (Y3); closed-loop coupling depth
  restricted to ECM half (Y4); AST-walk meta-test scoped to
  function body via ``textwrap.dedent(inspect.getsource(...))``
  (Y11).

Magic-Number Block: 2 module constants — both algebra-traceable.
``MAX_SQ_FROBENIUS_DIFF_PER_CELL = 9.0`` is the analytical
adversarial bound (T = ((-1, -1), (-1, -1)) against unit-vector
target at θ = π/4); ``MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL =
0.5 · MAX_SQ_FROBENIUS_DIFF_PER_CELL = 4.5`` is the energy-form
factor. Both grid-invariant in absolute terms (no dx/dt/grid_n
dependence); not chosen to fit any test target.

Phase E activation remains BLOCKED on the closed-loop ECM gate
composition cycle (separate unit combining Phase D no-op
wrapper + HB#1+#2 active law + HB#5 metric). This module
landing does NOT authorize Phase E composition itself; it
provides the Lyapunov-like evidence-channel primitive that
future Phase E composition will call.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np

from acs.v2.ecm_substrate import ECMSubstrateState

MAX_SQ_FROBENIUS_DIFF_PER_CELL: Final[float] = 9.0
MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL: Final[float] = (
    0.5 * MAX_SQ_FROBENIUS_DIFF_PER_CELL
)


@dataclass(frozen=True, slots=True)
class ECMOrientationLyapunovMetricDiagnostics:
    """Typed diagnostics for one HB#5 Lyapunov-like metric snapshot.

    10 fields covering counting + area accounting (for
    `V_active` interpretation), the schema-tight bound
    `v_active_max_bound_um2 = 4.5 · active_area_um2`, per-cell
    summaries for failure diagnosis, the explicit
    passive-cell orientation-magnitude diagnostic (Y1 forward
    guard — passive-cell magnitude is reported here, NOT in
    `v_active_um2`), and the Rule 10 stimulus-scale bridge.

    Attributes:
        n_active_cells: count of grid cells with non-zero
            traction stimulus this snapshot (count,
            dimensionless).
        n_total_cells: total grid cell count = `nx · ny`
            (count, dimensionless).
        cell_area_um2: per-cell area = `ecm.spacing_um²` (μm²).
        active_area_um2: `n_active_cells · cell_area_um2`
            (μm²).
        grid_area_um2: `n_total_cells · cell_area_um2` (μm²).
        v_active_max_bound_um2: schema-tight bound
            `4.5 · active_area_um2` (μm²) from the locked
            adversarial-bound derivation (per-cell
            `||T - n⊗n||²_F ≤ 9` ⇒ `V_per_cell ≤ 4.5`).
        v_active_per_cell_max_um2: max V contribution over
            active cells (μm²).
        v_active_per_cell_mean_um2: mean V contribution over
            active cells (μm²); 0.0 if no active cells.
        passive_orientation_magnitude_um2: passive-cell
            orientation magnitude as a diagnostic-only field
            (μm²); NOT part of `v_active_um2`. Y1 forward
            guard.
        max_traction_norm_nN_per_um2: max L2 norm of
            traction across the grid (nN/μm²); Rule 10 bridge
            to stimulus scale.
    """

    n_active_cells: int
    n_total_cells: int
    cell_area_um2: float
    active_area_um2: float
    grid_area_um2: float
    v_active_max_bound_um2: float
    v_active_per_cell_max_um2: float
    v_active_per_cell_mean_um2: float
    passive_orientation_magnitude_um2: float
    max_traction_norm_nN_per_um2: float


@dataclass(frozen=True, slots=True)
class ECMOrientationLyapunovMetricResult:
    """Outputs of one HB#5 Lyapunov-like metric snapshot."""

    v_active_um2: float
    diagnostics: ECMOrientationLyapunovMetricDiagnostics


def _validate_finite_traction(traction_density_xy: np.ndarray) -> np.ndarray:
    arr = np.asarray(traction_density_xy, dtype=np.float64)
    if not np.isfinite(arr).all():
        raise ValueError(
            "traction_density_xy contains non-finite values (nan / ±inf)"
        )
    return arr


def _validate_traction_shape(
    traction_density_xy: np.ndarray, grid_shape: tuple[int, int]
) -> None:
    expected = (*grid_shape, 2)
    if traction_density_xy.shape != expected:
        raise ValueError(
            f"traction_density_xy must have shape {expected!r}; got "
            f"{traction_density_xy.shape!r}"
        )


def compute_ecm_orientation_lyapunov_metric(
    ecm: ECMSubstrateState,
    traction_density_xy: np.ndarray,
) -> ECMOrientationLyapunovMetricResult:
    """Compute the HB#5 Lyapunov-like metric on an ECM snapshot.

    Pure read-only function. Returns a scalar Lyapunov-like
    energy ``v_active_um2`` over **active cells only** (cells
    with `traction_norm > 0`) plus a typed 10-field diagnostics
    record. Domain matching with HB#1+#2 update is enforced (Y1)
    so that zero-traction cells (where HB#1+#2 has `w = 0` and
    target placeholder is algebraically unused) do NOT
    contribute to V via a `T_target = 0` smuggling path.

    Two evidence channels (Y3, NOT to be conflated):

    - **fixed-target per-step contraction**: under one
      :func:`acs.v2.dynamics.ecm_constitutive_response.step_ecm_orientation_response`
      with the same traction field, `V_active_after ≤
      V_active_before` (test-time assertion, NOT in-function
      invariant). This is the Lyapunov-like decrease channel.
    - **varying-target boundedness**: under any varying-target
      sequence, `0 ≤ V_active ≤ v_active_max_bound_um2 = 4.5 ·
      active_area_um2` at every measured step. This is the
      boundedness channel.

    Validation contract (HB#1+#2 ``id=1486`` sister-precedent):
    ``ecm.validate()`` at function entry rejects invalid input
    ECM BEFORE any algebra; plain schema ``ValueError`` from
    ``ecm.validate()`` is NOT wrapped. Traction validation
    raises plain ``ValueError`` (no HB#5-specific error class
    per locked §5).

    Failure modes:

    - non-finite traction (`nan` / `±inf`): plain ``ValueError``.
    - traction shape mismatch with `(*ecm.grid_shape, 2)`:
      plain ``ValueError``.
    - inherited `ecm.validate()` failures: plain schema
      ``ValueError``.

    Args:
        ecm: input :class:`ECMSubstrateState`. Validated at
            entry. NOT mutated.
        traction_density_xy: ``(*ecm.grid_shape, 2)`` float64
            array in ``nN/μm²`` (HB#3 output contract).

    Returns:
        :class:`ECMOrientationLyapunovMetricResult` with the
        scalar `v_active_um2` (Lyapunov-like energy over active
        cells, μm²) and 10-field typed diagnostics.
    """

    ecm.validate()
    arr = _validate_finite_traction(traction_density_xy)
    _validate_traction_shape(arr, ecm.grid_shape)

    traction_norm = np.linalg.norm(arr, axis=-1)
    active = traction_norm > 0.0

    cell_area_um2 = float(ecm.spacing_um) ** 2
    n_xy = np.zeros_like(arr, dtype=np.float64)
    n_xy[active] = arr[active] / traction_norm[active, None]
    T_target = np.einsum("...i,...j->...ij", n_xy, n_xy)
    diff = np.asarray(ecm.orientation_tensor, dtype=np.float64) - T_target
    sq_diff_per_cell = np.sum(diff ** 2, axis=(-2, -1))

    if active.any():
        v_active = 0.5 * cell_area_um2 * float(sq_diff_per_cell[active].sum())
        per_cell_active = sq_diff_per_cell[active]
        v_per_cell_max = 0.5 * cell_area_um2 * float(per_cell_active.max())
        v_per_cell_mean = 0.5 * cell_area_um2 * float(per_cell_active.mean())
    else:
        v_active = 0.0
        v_per_cell_max = 0.0
        v_per_cell_mean = 0.0

    passive = ~active
    if passive.any():
        passive_mag_sq = float(
            np.sum(
                np.asarray(ecm.orientation_tensor, dtype=np.float64)[passive]
                ** 2
            )
        )
        passive_orientation_magnitude = 0.5 * cell_area_um2 * passive_mag_sq
    else:
        passive_orientation_magnitude = 0.0

    n_active = int(active.sum())
    n_total = int(np.prod(ecm.grid_shape))
    active_area = n_active * cell_area_um2
    grid_area = n_total * cell_area_um2
    v_active_max_bound = (
        MAX_LYAPUNOV_ENERGY_FACTOR_PER_ACTIVE_CELL * active_area
    )

    diagnostics = ECMOrientationLyapunovMetricDiagnostics(
        n_active_cells=n_active,
        n_total_cells=n_total,
        cell_area_um2=cell_area_um2,
        active_area_um2=active_area,
        grid_area_um2=grid_area,
        v_active_max_bound_um2=v_active_max_bound,
        v_active_per_cell_max_um2=v_per_cell_max,
        v_active_per_cell_mean_um2=v_per_cell_mean,
        passive_orientation_magnitude_um2=passive_orientation_magnitude,
        max_traction_norm_nN_per_um2=float(traction_norm.max()),
    )

    return ECMOrientationLyapunovMetricResult(
        v_active_um2=v_active,
        diagnostics=diagnostics,
    )
