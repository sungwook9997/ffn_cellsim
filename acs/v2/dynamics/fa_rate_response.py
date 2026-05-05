"""V2 Phase F Proper-1: ECM-Modulated FA Rate Response.

Per-FA HB#4 multipliers scale base FA rate constants (Y4-sister: NOT
force scalers). Effective-rate FA dynamics delegate to existing 6.3a
``step_focal_adhesions_static`` (Y6 traction algebra). Manual
aggregate of per-FA results (Y16: NO second delegate call to prevent
twice-mutation corruption).

Phase F minimal Y4 guard preserved: this is a NEW pathway, not a
patch to ``step_phase_f_minimal_motility``. HB#4 multipliers used as
per-rate scalers ONLY.

Locked at ``docs/v2_phase_f_proper_1_fa_rate_response_locked.md``
(commit ``89d6ddc``, Y1-Y20 reasoned-acceptance trace, A-tier full
debate output).

This module is **F-proper-1, Stage 1 of Phase F proper 3-stage**:
- F-proper-1 (this lock): rate dynamics, NO contour motion → cell
  motility 가시 NOT 발생
- F-proper-2 (separate lock): multi-rate motility integrator
- F-proper-3 (separate lock): PI-readable 15-min annotated pilot

PI mp4 (cell motility 가시)는 F-proper-3까지 가야 함. F-proper-1
milestone만으로 PI에게 "cell motility 보임" wording NOT.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np

from acs.v2.active_contour import _DT_RATE_SAFETY_MARGIN
from acs.v2.dynamics.ecm_to_fa_bias import ECMToFABiasResult, RATE_NAMES
from acs.v2.dynamics.focal_adhesion import (
    FocalAdhesionDynamicsDiagnostics,
    FocalAdhesionDynamicsError,
    FocalAdhesionDynamicsParameters,
    FocalAdhesionDynamicsResult,
    step_focal_adhesions_static,
)
from acs.v2.focal_adhesion import FocalAdhesionState


@dataclass(frozen=True, slots=True)
class FARateResponseDiagnostics:
    """Typed diagnostics for one F-proper-1 step."""

    n_adhesions: int
    effective_k_maturity_per_s: np.ndarray
    effective_k_bind_per_s: np.ndarray
    effective_k_unbind_per_s: np.ndarray
    max_effective_rate_per_s: float
    max_dt_fa_rate_product: float
    aggregate_traction_change_nN_xy: tuple[float, float]
    n_none_base_rates: int


@dataclass(frozen=True, slots=True)
class FARateResponseResult:
    """One F-proper-1 step result.

    ``fa_dynamics_result`` is a manually-aggregated
    :class:`FocalAdhesionDynamicsResult` (Y16: NO second delegate call
    to ``step_focal_adhesions_static`` after the per-FA effective-rate
    steps; that would silently twice-mutate / corrupt the per-FA
    result).
    """

    updated_adhesions: tuple[FocalAdhesionState, ...]
    fa_dynamics_result: FocalAdhesionDynamicsResult
    diagnostics: FARateResponseDiagnostics


def _validate_base_params_for_proper(
    base_params: FocalAdhesionDynamicsParameters,
) -> None:
    """Y20 fail-closed validation BEFORE empty-adhesion return path.

    Schema-validation cannot be bypassed via empty adhesions; calling
    ``base_params.validate()`` at function entry guarantees the locked
    contract holds even when the per-FA loop is a no-op.
    """
    base_params.validate()


def _scale_rate_or_none(
    base_rate: Optional[float], multiplier: float
) -> Optional[float]:
    """Y17: ``None`` base rate stays ``None`` (= "no update" in 6.3a
    contract); float scales by ``multiplier``."""
    if base_rate is None:
        return None
    return float(base_rate) * float(multiplier)


def _validate_multipliers_match_adhesions(
    multipliers: ECMToFABiasResult,
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
) -> None:
    """Y14 fail-closed validation; no silent reorder in F-proper-1."""

    expected_ids = tuple(fa.adhesion_id for fa in adhesions)
    if multipliers.fa_ids != expected_ids:
        raise FocalAdhesionDynamicsError(
            "multipliers_fa_ids_mismatch",
            f"multipliers.fa_ids {multipliers.fa_ids!r} != adhesions ids "
            f"{expected_ids!r}. Caller must align order; silent reorder "
            f"forbidden in F-proper-1 first lock.",
        )
    if tuple(multipliers.rate_names) != tuple(RATE_NAMES):
        raise FocalAdhesionDynamicsError(
            "multipliers_rate_names_mismatch",
            f"multipliers.rate_names {multipliers.rate_names!r} != "
            f"RATE_NAMES {RATE_NAMES!r}. Strict ordering required.",
        )
    n_fa = len(adhesions)
    expected_shape = (n_fa, len(RATE_NAMES))
    arr = np.asarray(multipliers.multipliers_per_fa)
    if arr.shape != expected_shape:
        raise FocalAdhesionDynamicsError(
            "multipliers_shape_mismatch",
            f"multipliers.multipliers_per_fa.shape {arr.shape!r} != "
            f"expected {expected_shape!r}",
        )
    if arr.dtype == bool:
        raise FocalAdhesionDynamicsError(
            "multipliers_bool_dtype",
            "multipliers must not be bool dtype (Python bool subset int trap)",
        )
    arr_f = np.asarray(arr, dtype=np.float64)
    if not np.all(np.isfinite(arr_f)):
        raise FocalAdhesionDynamicsError(
            "multipliers_non_finite",
            f"multipliers must be finite; got {arr!r}",
        )
    if (arr_f < 0).any():
        raise FocalAdhesionDynamicsError(
            "multipliers_negative",
            f"multipliers must be non-negative; got {arr!r}",
        )


def step_fa_rate_response(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    centroid_um_xy: tuple[float, float],
    multipliers: ECMToFABiasResult,
    base_params: FocalAdhesionDynamicsParameters,
    *,
    traction_axis_xy_per_fa: Optional[list] = None,
) -> FARateResponseResult:
    """F-proper-1 step.

    Pipeline (Y20 base_params validation BEFORE empty-adhesion return):

    0. Validate base_params (NEVER bypass — Y20)
    1. Validate multipliers (Y14): fa_ids exact match, rate_names match,
       shape, finite/non-bool/non-negative
    2. Validate ``traction_axis_xy_per_fa`` length
    3. Build per-FA effective rates (Y17 None handling)
    4. dt-rate gate via ``base_params.dt_fa_s * max_effective_rate``
       (≤ ``_DT_RATE_SAFETY_MARGIN``; Y6 sister)
    5. Per-FA delegate to :func:`step_focal_adhesions_static`; manual
       aggregate (Y16 — NO second delegate call)
    6. Diagnostics + return

    Y4 minimal guard preserved: this is a NEW pathway. HB#4 multipliers
    are used as per-rate scalers ONLY (per ``RATE_NAMES`` original
    design intent). Force change through 6.3a traction algebra only.
    """

    _validate_base_params_for_proper(base_params)

    _validate_multipliers_match_adhesions(multipliers, adhesions)
    n_fa = len(adhesions)

    if traction_axis_xy_per_fa is not None and len(traction_axis_xy_per_fa) != n_fa:
        raise FocalAdhesionDynamicsError(
            "traction_axis_length_mismatch",
            f"traction_axis_xy_per_fa length "
            f"{len(traction_axis_xy_per_fa)} != n_fa {n_fa}",
        )

    n_none = 0
    eff_k_maturity_list: list[Optional[float]] = []
    eff_k_bind_list: list[Optional[float]] = []
    eff_k_unbind_list: list[Optional[float]] = []
    for i in range(n_fa):
        m = multipliers.multipliers_per_fa[i]
        eff_k_maturity_list.append(
            _scale_rate_or_none(base_params.k_maturity_per_s, float(m[0]))
        )
        eff_k_bind_list.append(
            _scale_rate_or_none(base_params.k_bind_per_s, float(m[1]))
        )
        eff_k_unbind_list.append(
            _scale_rate_or_none(base_params.k_unbind_per_s, float(m[2]))
        )
        if base_params.k_maturity_per_s is None:
            n_none += 1
        if base_params.k_bind_per_s is None:
            n_none += 1
        if base_params.k_unbind_per_s is None:
            n_none += 1

    eff_k_maturity_arr = np.array(
        [r if r is not None else 0.0 for r in eff_k_maturity_list],
        dtype=np.float64,
    )
    eff_k_bind_arr = np.array(
        [r if r is not None else 0.0 for r in eff_k_bind_list],
        dtype=np.float64,
    )
    eff_k_unbind_arr = np.array(
        [r if r is not None else 0.0 for r in eff_k_unbind_list],
        dtype=np.float64,
    )

    if n_fa == 0:
        max_rate = 0.0
    else:
        candidates: list[float] = []
        if eff_k_maturity_arr.size:
            candidates.append(float(eff_k_maturity_arr.max()))
        if eff_k_bind_arr.size:
            candidates.append(float(eff_k_bind_arr.max()))
        if eff_k_unbind_arr.size:
            candidates.append(float(eff_k_unbind_arr.max()))
        max_rate = float(max(candidates)) if candidates else 0.0

    dt_fa_s = float(base_params.dt_fa_s)
    dt_rate_product = dt_fa_s * max_rate
    if dt_rate_product > _DT_RATE_SAFETY_MARGIN:
        raise FocalAdhesionDynamicsError(
            "dt_rate_violation",
            f"dt_fa_s={dt_fa_s!r} * max_effective_rate={max_rate:.4f} = "
            f"{dt_rate_product:.4f} > safety margin {_DT_RATE_SAFETY_MARGIN}",
        )

    if n_fa > 0:
        pre_traction = np.array(
            [fa.traction_force_nN_xy for fa in adhesions], dtype=np.float64
        )
        pre_traction_sum = pre_traction.sum(axis=0)
    else:
        pre_traction = np.zeros((0, 2), dtype=np.float64)
        pre_traction_sum = np.zeros(2, dtype=np.float64)

    cell_force = np.zeros((n_fa, 2), dtype=np.float64)
    substrate_reaction = np.zeros((n_fa, 2), dtype=np.float64)
    radial_components = np.zeros(n_fa, dtype=np.float64)
    tangential_components = np.zeros(n_fa, dtype=np.float64)
    updated: list[FocalAdhesionState] = []

    for i, fa in enumerate(adhesions):
        per_fa_params = replace(
            base_params,
            k_maturity_per_s=eff_k_maturity_list[i],
            k_bind_per_s=eff_k_bind_list[i],
            k_unbind_per_s=eff_k_unbind_list[i],
        )
        axis = (
            traction_axis_xy_per_fa[i]
            if traction_axis_xy_per_fa is not None
            else None
        )
        result_one = step_focal_adhesions_static(
            (fa,),
            centroid_um_xy,
            per_fa_params,
            traction_axis_xy_per_fa=[axis],
        )
        updated.append(result_one.updated_adhesions[0])
        cell_force[i] = result_one.cell_force_nN_xy[0]
        substrate_reaction[i] = result_one.substrate_reaction_nN_xy[0]
        radial_components[i] = result_one.radial_components[0]
        tangential_components[i] = result_one.tangential_components[0]

    updated_tuple = tuple(updated)

    if n_fa > 0:
        aggregate_cell_force = tuple(cell_force.sum(axis=0).tolist())
        aggregate_substrate_reaction = tuple(substrate_reaction.sum(axis=0).tolist())
        max_traction_magnitude = float(
            np.linalg.norm(cell_force, axis=1).max()
        )
        post_traction_sum = cell_force.sum(axis=0)
    else:
        aggregate_cell_force = (0.0, 0.0)
        aggregate_substrate_reaction = (0.0, 0.0)
        max_traction_magnitude = 0.0
        post_traction_sum = np.zeros(2, dtype=np.float64)

    fa_dynamics_result = FocalAdhesionDynamicsResult(
        updated_adhesions=updated_tuple,
        cell_force_nN_xy=cell_force,
        substrate_reaction_nN_xy=substrate_reaction,
        radial_components=radial_components,
        tangential_components=tangential_components,
        diagnostics=FocalAdhesionDynamicsDiagnostics(
            n_adhesions=n_fa,
            aggregate_cell_force_nN_xy=aggregate_cell_force,
            aggregate_substrate_reaction_nN_xy=aggregate_substrate_reaction,
            max_traction_magnitude_nN=max_traction_magnitude,
        ),
    )

    aggregate_change = tuple(
        (post_traction_sum - pre_traction_sum).tolist()
    )

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
