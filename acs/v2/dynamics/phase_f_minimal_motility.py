"""V2 Phase F minimal cell motility pilot — first cell-moving v2 sample.

Couples P1 active contour + Phase E v2 wrapper via FA-vertex external
force (Newton's 3rd law reaction). HB#4 rate multipliers are
diagnostic-only; ECM->cell motility feedback is NOT implemented here
(future Phase F proper, separate A-tier cycle).

Locked at ``docs/v2_phase_f_minimal_cell_motility_pilot_locked.md``
(commit ``7f8c171``); A-tier-compressed cycle per PI ASAP directive
(``id=1894``/``id=1903``). 2 design rounds + Codex single review pass
+ seal; no Sanity Gate doc per compressed protocol.

This pilot shows FA-traction-driven contour motion with Phase E v2
ECM-side diagnostics running in the loop. It does not yet implement
ECM->cell motility feedback, because HB#4 rate multipliers are
diagnostic-only here.

The locked design **forbids** (Y4/Z5 critical guard):

- Using HB#4 rate multipliers as traction-force scalers in any form
  (mean / sum / product / lookup). Multipliers are strictly carried
  through ``result.phase_e_v2.ecm_to_fa_bias.diagnostics_dict`` for
  observation only.
- Sign-flipping ``traction_force_nN_xy`` to make the movie look
  better. The schema convention is cell-on-substrate; Newton's 3rd
  law reaction on the cell vertex is ``-traction_force_nN_xy``. For
  visible +x cell motion, set the fixture FA traction to ``(-1, 0)``.
- A separate ``dt_s`` argument. The synchronous timestep comes from
  ``contour_state.params.dt_cell_s`` so P1 and Phase E v2 evolve on
  the same dt without subcycling (deferred to Phase F proper).
- Modifying P1 ``acs.v2.dynamics.active_contour.step`` signature.
  This module reuses P1's internal force functions
  (``compute_cortex_forces``, ``compute_area_forces``,
  ``compute_rate_max``) without touching the locked P1 step contract.

Failure semantics are sister-consistent with P1: dt safety-margin
breaches raise :class:`acs.v2.dynamics.active_contour.ActiveContourStepError`
with ``failure_kind="dt_violation"`` (Y7).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

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

    HB#4 multipliers in
    ``phase_e_v2.ecm_to_fa_bias.diagnostics_dict`` are
    diagnostic-only; no traction-force scaling is performed (Y4/Z5
    guard). This pilot shows FA-traction-driven contour motion with
    Phase E v2 ECM-side diagnostics running in the loop. It does
    not yet implement ECM->cell motility feedback, because HB#4
    rate multipliers are diagnostic-only here.

    Attributes:
        contour: post-step :class:`ActiveContourState` (P1 vertices
            updated by overdamped Euler with combined internal +
            external forces).
        adhesions: tuple of :class:`FocalAdhesionState` with
            positions snapped to attached vertices' new positions
            (no FA dynamics; count/state preserved).
        phase_e_v2: :class:`PhaseEStepResult` from the locked
            ``step_closed_loop_phase_e_v2`` call. The
            ``ecm_to_fa_bias.diagnostics_dict`` carries HB#4
            multiplier summaries for observation only.
        updated_ecm: equal to ``phase_e_v2.updated_ecm`` by object
            identity (Phase E v2 Y4 invariant inherited).
    """

    contour: ActiveContourState
    adhesions: tuple[FocalAdhesionState, ...]
    phase_e_v2: PhaseEStepResult
    updated_ecm: ECMSubstrateState


def _validate_fa_to_vertex_index(
    fa_to_vertex_index: dict[str, int],
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    n_vertices: int,
) -> None:
    """Y1 fail-closed validation of FA-to-vertex mapping."""

    fa_ids = {fa.adhesion_id for fa in adhesions}
    mapping_ids = set(fa_to_vertex_index.keys())
    unknown = mapping_ids - fa_ids
    if unknown:
        raise ValueError(
            f"fa_to_vertex_index references unknown FA ids: {sorted(unknown)}"
        )
    missing = fa_ids - mapping_ids
    if missing:
        raise ValueError(
            f"fa_to_vertex_index missing FA ids: {sorted(missing)}"
        )
    for fa_id, v_idx in fa_to_vertex_index.items():
        if isinstance(v_idx, bool):
            raise ValueError(
                f"vertex index for FA {fa_id!r} must be int (not bool — "
                f"Python bool subset int trap), got {type(v_idx).__name__}"
            )
        if not isinstance(v_idx, int):
            raise ValueError(
                f"vertex index for FA {fa_id!r} must be int, "
                f"got {type(v_idx).__name__}"
            )
        if not (0 <= v_idx < n_vertices):
            raise ValueError(
                f"FA {fa_id!r} vertex index {v_idx} outside [0, {n_vertices})"
            )


def _step_active_contour_with_external_force(
    contour_state: ActiveContourState,
    external_forces_per_vertex_nN: np.ndarray,
) -> ActiveContourState:
    """Phase F-local helper (Y5). Same overdamped Euler update as P1
    step, plus external force per vertex. Reuses locked
    :func:`compute_cortex_forces` / :func:`compute_area_forces` /
    :func:`compute_rate_max` without modifying P1 step signature.

    Failure semantics sister-consistent with P1 (Y7): raises
    :class:`ActiveContourStepError` with ``failure_kind="dt_violation"``
    on dt-rate-product safety margin breach (same threshold
    ``_DT_RATE_SAFETY_MARGIN`` as P1).
    """

    # First gate (P1 sister-consistent baseline): internal-only dt safety
    # margin via compute_rate_max. This preserves the P1 step's contract
    # so a Phase F caller's contour_state alone cannot violate it
    # silently.
    rate_info = compute_rate_max(contour_state)
    if rate_info["dt_rate_product"] > _DT_RATE_SAFETY_MARGIN:
        raise ActiveContourStepError(
            "dt_violation",
            f"dt_cell_s={contour_state.params.dt_cell_s}*rate_max product "
            f"{rate_info['dt_rate_product']:.4f} > safety margin "
            f"{_DT_RATE_SAFETY_MARGIN} (internal-only)",
        )

    internal = compute_cortex_forces(contour_state) + compute_area_forces(
        contour_state
    )

    if external_forces_per_vertex_nN.shape != internal.shape:
        raise ValueError(
            f"external_forces_per_vertex_nN shape "
            f"{external_forces_per_vertex_nN.shape!r} != internal force "
            f"shape {internal.shape!r}"
        )

    forces = internal + external_forces_per_vertex_nN
    zeta = rate_info["zeta_nN_s_per_um"]
    dt_cell_s = float(contour_state.params.dt_cell_s)

    # Second gate (Phase F new physics): combined internal+external dt
    # safety margin. Required because the new external-force layer is
    # outside P1's compute_rate_max scope; without this gate a large FA
    # traction would escape the dt_violation contract and surface as a
    # downstream MeasurementBoundary geometry error (Codex id=1986 P0
    # BLOCKER fix). Per-vertex rate magnitude = |forces_i| / zeta_i.
    per_vertex_rate = np.linalg.norm(forces, axis=-1) / zeta
    combined_rate_max = float(per_vertex_rate.max())
    combined_dt_rate_product = dt_cell_s * combined_rate_max
    if combined_dt_rate_product > _DT_RATE_SAFETY_MARGIN:
        raise ActiveContourStepError(
            "dt_violation",
            f"dt_cell_s={dt_cell_s}*combined_rate_max product "
            f"{combined_dt_rate_product:.4f} > safety margin "
            f"{_DT_RATE_SAFETY_MARGIN} (internal+external)",
        )

    new_vertices = (
        contour_state.vertices_xy_um + dt_cell_s * forces / zeta[:, None]
    )

    return replace(contour_state, vertices_xy_um=new_vertices)


def step_phase_f_minimal_motility(
    contour_state: ActiveContourState,
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
    fa_to_vertex_index: dict[str, int],
    *,
    k_active: float,
) -> PhaseFStepResult:
    """One synchronous Phase F minimal motility step.

    Uses ``contour_state.params.dt_cell_s`` as the synchronous dt for
    both P1 active contour and the Phase E v2 wrapper (no subcycling
    — that is deferred to a future Phase F proper A-tier cycle). HB#4
    rate multipliers are diagnostic-only; FA traction drives cell
    vertex motion through Newton's 3rd law reaction directly (no
    multiplier scaling).

    This pilot shows FA-traction-driven contour motion with Phase E v2
    ECM-side diagnostics running in the loop. It does not yet
    implement ECM->cell motility feedback, because HB#4 rate
    multipliers are diagnostic-only here.

    Args:
        contour_state: P1 :class:`ActiveContourState` with finite
            positive ``dt_cell_s``.
        adhesions: tuple or list of :class:`FocalAdhesionState`. Empty
            list valid (no external force; pure P1 dynamics).
        ecm: input :class:`ECMSubstrateState`.
        fa_to_vertex_index: dict mapping each FA's
            ``adhesion_id`` to a vertex index in
            ``contour_state.vertices_xy_um``. Validated fail-closed
            (Y1) — every FA must have an entry, and every entry must
            be a non-bool int in ``[0, n_vertices)``.
        k_active: HB#4-active coupling strength forwarded to the
            Phase E v2 wrapper. Multipliers from this readout are
            **diagnostic-only** in Phase F minimal pilot (Y4 guard).

    Returns:
        :class:`PhaseFStepResult` with post-step contour + adhesions
        snapped to new vertex positions + Phase E v2 result + ECM
        identity-equal to ``phase_e_v2.updated_ecm``.

    Raises:
        ValueError: if ``fa_to_vertex_index`` is invalid (Y1) or if
            external force shape mismatches internal force shape.
        ActiveContourStepError: with ``failure_kind="dt_violation"``
            on dt-rate safety-margin breach (Y7 P1 sister-consistent).
        Errors from sub-calls (HB#3 / HB#1+#2 / HB#4-active / HB#5)
            propagate verbatim (no wrapper-level failure kind).
    """

    n_vertices = int(contour_state.vertices_xy_um.shape[0])
    _validate_fa_to_vertex_index(fa_to_vertex_index, adhesions, n_vertices)

    dt_s = float(contour_state.params.dt_cell_s)

    pe_v2_result = step_closed_loop_phase_e_v2(
        adhesions,
        ecm,
        dt_s,
        k_active=k_active,
    )

    external_forces = np.zeros((n_vertices, 2), dtype=np.float64)
    for fa in adhesions:
        v_idx = fa_to_vertex_index[fa.adhesion_id]
        # Newton's 3rd law: cell vertex receives -F where F is the
        # cell-on-substrate traction (FA schema convention). Y4 guard:
        # NO HB#4 multiplier scaling here (diagnostic-only).
        external_forces[v_idx] += -np.asarray(
            fa.traction_force_nN_xy, dtype=np.float64
        )

    contour_new = _step_active_contour_with_external_force(
        contour_state,
        external_forces_per_vertex_nN=external_forces,
    )

    adhesions_new: list[FocalAdhesionState] = []
    for fa in adhesions:
        v_idx = fa_to_vertex_index[fa.adhesion_id]
        new_position = (
            float(contour_new.vertices_xy_um[v_idx, 0]),
            float(contour_new.vertices_xy_um[v_idx, 1]),
        )
        adhesions_new.append(replace(fa, position_um_xy=new_position))

    return PhaseFStepResult(
        contour=contour_new,
        adhesions=tuple(adhesions_new),
        phase_e_v2=pe_v2_result,
        updated_ecm=pe_v2_result.updated_ecm,
    )
