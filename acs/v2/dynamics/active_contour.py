"""Active contour dynamics: cortex + area forces and overdamped Euler step.

Implements the locked P1 alpha force law from
``docs/v2_p1_derivation_locked.md`` §1, §2, §3 and the runtime contract
from §4. The schema layer (``acs.v2.active_contour``) owns the
parameter contract and the analytic Gershgorin rate bound; this
module only computes per-vertex forces and advances the state by one
overdamped Euler step.

Sanity Gate scope (dynamics/active_contour.py):

- §1 dimensional analysis: ``compute_cortex_forces`` returns
  per-vertex force in nN (verified via the unit chain
  ``λ_c [nN] · (unit tangent change [dimensionless]) → nN``);
  ``compute_area_forces`` returns per-vertex force in nN
  (``p_A [nN/μm] · ∇_i A [μm] → nN``); ``step`` divides by
  ``ζ_i [nN·s/μm]`` and multiplies by ``dt_cell [s]`` to give a
  vertex displacement in μm.
- §2 timestep stability: ``step`` re-runs ``compute_rate_max`` from
  the schema module and aborts when ``dt_rate_product`` exceeds the
  named ``_DT_RATE_SAFETY_MARGIN``. No auto-shrink.
- §3 boundary cases: every step exports through
  ``to_measurement_boundary`` so a self-intersection / pinched / too-short
  edge introduced by the step raises ``MeasurementBoundaryError``
  rather than silently corrupting the state.
- §5 sign / sense: cortex force per the locked tangent-difference
  form is contractile (regular polygon shrinks under cortex-only).
  Area force per the gradient form pushes outward when ``A < A_0``
  and inward when ``A > A_0``. Both are unit-tested in
  ``tests/test_v2_active_contour_dynamics.py``.
- §6 measurement-protocol consistency: Hard Rule 11. The post-step
  validation routes through ``MeasurementBoundary.from_array`` so
  the canonical world_um_y_up CCW polygon is enforced after every
  integration step.

Magic-Number Block: this module declares no tunable numeric. The
runtime safety margin lives in ``acs.v2.active_contour`` as
``_DT_RATE_SAFETY_MARGIN`` and is reused here, not redefined.
"""

from __future__ import annotations

import numpy as np

from acs.v2.active_contour import (
    _DT_RATE_SAFETY_MARGIN,
    ActiveContourParametersError,
    ActiveContourState,
    compute_rate_max,
)


class ActiveContourStepError(RuntimeError):
    """Raised when a step violates the locked runtime contract.

    Carries a ``failure_kind`` so the caller / test harness can switch
    on the specific contract violation.
    """

    def __init__(self, failure_kind: str, message: str) -> None:
        self.failure_kind = failure_kind
        super().__init__(message)


def compute_cortex_forces(state: ActiveContourState) -> np.ndarray:
    """Per-vertex cortex force ``F_i^c = λ_c · (t_i − t_{i−1})`` in nN.

    Locked form from ``docs/v2_p1_derivation_locked.md`` §1: cortex
    energy ``E_c = λ_c · P`` with perimeter ``P``; the variational
    derivative gives a vertex force equal to the change in unit
    tangents at that vertex. The form never computes a curvature
    first (gate §1 anti-pattern note).

    Returns an ``(N, 2)`` float64 array. When ``λ_c == 0`` the
    returned array is exactly zero everywhere.
    """

    lambda_c = state.params.lambda_c_resolved_nN
    verts = state.vertices_xy_um
    n = state.params.n_vertices
    edge_vec = np.roll(verts, -1, axis=0) - verts  # r_{i+1} - r_i
    edge_len = np.linalg.norm(edge_vec, axis=1)
    if (edge_len == 0.0).any():
        raise ActiveContourParametersError(
            "edge_length_zero",
            "active contour state has zero-length edges; cannot compute cortex force",
        )
    t_i = edge_vec / edge_len[:, None]  # unit tangent on edge i→i+1
    t_im1 = np.roll(t_i, 1, axis=0)  # unit tangent on edge i-1→i
    return lambda_c * (t_i - t_im1)


def compute_area_forces(state: ActiveContourState) -> np.ndarray:
    """Per-vertex area-restoring force ``F_i^A = -p_A · ∇_i A`` in nN.

    Locked form from ``docs/v2_p1_derivation_locked.md`` §2 with
    pressure-like scalar ``p_A = K_A · (A − A_0) / A_0`` and area
    gradient ``∇_i A = ½ (y_{i+1} − y_{i−1}, x_{i−1} − x_{i+1})``.

    Returns an ``(N, 2)`` float64 array. When ``K_A == 0`` or
    ``A == A_0`` the returned array is exactly zero everywhere.
    """

    k_a = state.params.k_a_nN_per_um
    target_area = state.params.target_area_um2
    if k_a == 0.0:
        return np.zeros_like(state.vertices_xy_um)
    area = state.area_um2()
    p_a = k_a * (area - target_area) / target_area
    grad = np.column_stack(
        (
            0.5 * (np.roll(state.vertices_xy_um[:, 1], -1)
                   - np.roll(state.vertices_xy_um[:, 1], 1)),
            0.5 * (np.roll(state.vertices_xy_um[:, 0], 1)
                   - np.roll(state.vertices_xy_um[:, 0], -1)),
        )
    )
    return -p_a * grad


def energy_cortex(state: ActiveContourState) -> float:
    """``E_c = λ_c · P`` in nN·μm."""

    return float(state.params.lambda_c_resolved_nN * state.perimeter_um())


def energy_area(state: ActiveContourState) -> float:
    """``E_A = ½ K_A · (A − A_0)² / A_0`` in nN·μm."""

    k_a = state.params.k_a_nN_per_um
    target_area = state.params.target_area_um2
    delta = state.area_um2() - target_area
    return float(0.5 * k_a * delta * delta / target_area)


def step(state: ActiveContourState) -> ActiveContourState:
    """Advance the state by one overdamped Euler integration step.

    Raises :class:`ActiveContourStepError` (``failure_kind="dt_violation"``)
    when ``dt_cell · rate_max`` exceeds ``_DT_RATE_SAFETY_MARGIN``;
    the locked SOP is "do NOT auto-adjust" — the caller must reduce
    ``dt_cell`` or switch integrator.

    Raises :class:`MeasurementBoundaryError` (via
    :meth:`ActiveContourState.to_measurement_boundary`) when the post-step
    geometry violates the canonical contract (self-intersection,
    pinched polygon, sub-``min_edge_um`` edge, …). The state passed
    in is never mutated; the returned state is a fresh
    :class:`ActiveContourState`.
    """

    rate_info = compute_rate_max(state)
    if rate_info["dt_rate_product"] > _DT_RATE_SAFETY_MARGIN:
        raise ActiveContourStepError(
            "dt_violation",
            f"dt_cell · rate_max = {rate_info['dt_rate_product']:.6e} exceeds "
            f"safety margin {_DT_RATE_SAFETY_MARGIN!r}; reduce dt_cell or "
            f"switch to an implicit integrator (auto-shrink is not allowed)",
        )

    forces = compute_cortex_forces(state) + compute_area_forces(state)
    zeta = rate_info["zeta_nN_s_per_um"]  # nN·s/μm, per vertex
    velocity = forces / zeta[:, None]  # μm/s
    new_vertices = state.vertices_xy_um + state.params.dt_cell_s * velocity

    new_state = ActiveContourState(
        cell_id=state.cell_id,
        vertices_xy_um=new_vertices,
        params=state.params,
        step_count=state.step_count + 1,
    )
    # Post-step measurement-boundary validation (Hard Rule 11): a step
    # that produces a self-intersection or sub-floor edge surfaces
    # immediately rather than silently corrupting state.
    new_state.to_measurement_boundary()
    return new_state
