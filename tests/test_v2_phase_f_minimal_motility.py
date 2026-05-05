"""B-tier 테스트 (A-tier-compressed cycle): Phase F minimal cell motility pilot.

12 tests per locked §4 (`docs/v2_phase_f_minimal_cell_motility_pilot_locked.md`,
commit ``7f8c171``) + 2 추가 acceptance tests per Codex `id=1953` / PI
`id=1949` rigid translation 챌린지:

- Composition + validation (4): tests 1-4
- Sign + force coupling (3): tests 5-7
- HB#4 multipliers diagnostic-only guard (1, Y4 + Y10 wrapper seam): test 8
- FA position + ECM evolution (2): tests 9-10
- Fixture-specific evolution (1, Y9 required): test 11
- Exports (1): test 12
- PI acceptance tightening (2, Codex id=1953): tests 13-14
"""

from __future__ import annotations

import numpy as np
import pytest

import acs.v2 as v2_pkg
from acs.v2 import dynamics as v2_dynamics
from acs.v2.active_contour import (
    ActiveContourParameters,
    ActiveContourState,
    regular_polygon_vertices,
)
from acs.v2.dynamics import phase_f_minimal_motility as pf_mod
from acs.v2.dynamics.active_contour import ActiveContourStepError
from acs.v2.dynamics.closed_loop_phase_e import (
    PhaseEStepResult,
    step_closed_loop_phase_e_v2,
)
from acs.v2.dynamics.phase_f_minimal_motility import (
    PhaseFStepResult,
    _step_active_contour_with_external_force,
    step_phase_f_minimal_motility,
)
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState


def _make_test_contour(
    *,
    n_vertices: int = 12,
    radius_um: float = 1.0,
    center=(2.0, 2.0),
    k_a_nN_per_um: float = 0.0,
    target_area_um2: float = np.pi,
    dt_cell_s: float = 1e-3,
    lambda_c_nN: float = 0.0,
    sigma_c_nN_per_um: float = 0.0,
    xi_line_nN_s_per_um2: float = 1e-3,
    effective_height_um: float = 1.0,
) -> ActiveContourState:
    params = ActiveContourParameters(
        k_a_nN_per_um=k_a_nN_per_um,
        target_area_um2=target_area_um2,
        dt_cell_s=dt_cell_s,
        n_vertices=n_vertices,
        lambda_c_nN=lambda_c_nN,
        sigma_c_nN_per_um=sigma_c_nN_per_um,
        xi_line_nN_s_per_um2=xi_line_nN_s_per_um2,
        effective_height_um=effective_height_um,
    ).validate()
    vertices = regular_polygon_vertices(n_vertices, radius_um, center=center)
    return ActiveContourState(
        cell_id="cell-pf",
        params=params,
        vertices_xy_um=vertices,
    )


def _make_test_ecm(
    *, grid_n: int = 8, spacing_um: float = 0.5, isotropic_value: float = 0.5
) -> ECMSubstrateState:
    orientation = np.zeros((grid_n, grid_n, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = isotropic_value
    orientation[..., 1, 1] = isotropic_value
    return ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=spacing_um,
        stiffness_kpa=np.zeros((grid_n, grid_n), dtype=np.float64),
        ligand_density=np.zeros((grid_n, grid_n), dtype=np.float64),
        fiber_density=np.zeros((grid_n, grid_n), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((grid_n, grid_n), dtype=np.float64),
        source="synthetic",
    )


def _make_test_fa(
    adhesion_id: str = "fa-0",
    *,
    position=(2.0, 2.0),
    traction=(1.0, 0.0),
) -> FocalAdhesionState:
    return FocalAdhesionState(
        adhesion_id=adhesion_id,
        cell_id="cell-pf",
        position_um_xy=position,
        age_s=0.0,
        maturity=1.0,
        bound_fraction=1.0,
        state="mature",
        traction_force_nN_xy=traction,
    )


# ---------------------------------------------------------------------------
# Composition + validation (4)
# ---------------------------------------------------------------------------


def test_phase_f_step_returns_phase_f_step_result_dataclass():
    """Test 1: returns PhaseFStepResult with all 4 fields populated +
    identity invariant updated_ecm is phase_e_v2.updated_ecm."""
    contour = _make_test_contour()
    ecm = _make_test_ecm()
    # Place FA at vertex 0 position
    v0 = contour.vertices_xy_um[0]
    fa = _make_test_fa(position=tuple(v0.tolist()), traction=(1.0, 0.0))
    fa_to_vertex = {fa.adhesion_id: 0}

    result = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.0
    )
    assert isinstance(result, PhaseFStepResult)
    assert isinstance(result.contour, ActiveContourState)
    assert isinstance(result.adhesions, tuple)
    assert isinstance(result.phase_e_v2, PhaseEStepResult)
    assert isinstance(result.updated_ecm, ECMSubstrateState)
    assert result.updated_ecm is result.phase_e_v2.updated_ecm


def test_phase_f_step_unknown_fa_id_raises():
    """Test 2 (Y1): fa_to_vertex_index references FA id not in adhesions."""
    contour = _make_test_contour()
    ecm = _make_test_ecm()
    v0 = contour.vertices_xy_um[0]
    fa = _make_test_fa(position=tuple(v0.tolist()))
    fa_to_vertex = {"unknown-fa": 0, fa.adhesion_id: 0}
    with pytest.raises(ValueError, match="unknown FA ids"):
        step_phase_f_minimal_motility(
            contour, (fa,), ecm, fa_to_vertex, k_active=1.0
        )


def test_phase_f_step_vertex_index_out_of_range_raises():
    """Test 3 (Y1): vertex index outside [0, n_vertices) -> ValueError."""
    contour = _make_test_contour(n_vertices=8)
    ecm = _make_test_ecm()
    v0 = contour.vertices_xy_um[0]
    fa = _make_test_fa(position=tuple(v0.tolist()))
    fa_to_vertex = {fa.adhesion_id: 999}
    with pytest.raises(ValueError, match="outside"):
        step_phase_f_minimal_motility(
            contour, (fa,), ecm, fa_to_vertex, k_active=1.0
        )


def test_phase_f_step_dt_violation_raises_active_contour_step_error():
    """Test 4 (Y7): dt safety-margin breach -> P1 sister-consistent
    ActiveContourStepError(failure_kind='dt_violation')."""
    # Force dt_violation: large dt_cell_s + large lambda_c (consistent with
    # sigma_c * H_eff per P1 cortex param invariant) -> dt*rate_max > 0.5
    contour = _make_test_contour(
        dt_cell_s=10.0, lambda_c_nN=10.0, sigma_c_nN_per_um=10.0
    )
    ecm = _make_test_ecm()
    v0 = contour.vertices_xy_um[0]
    fa = _make_test_fa(position=tuple(v0.tolist()))
    fa_to_vertex = {fa.adhesion_id: 0}
    with pytest.raises(ActiveContourStepError) as excinfo:
        step_phase_f_minimal_motility(
            contour, (fa,), ecm, fa_to_vertex, k_active=1.0
        )
    assert excinfo.value.failure_kind == "dt_violation"


# ---------------------------------------------------------------------------
# Sign + force coupling (3)
# ---------------------------------------------------------------------------


def test_phase_f_step_traction_one_zero_moves_vertex_negative_x():
    """Test 5 (Y2): isolated `lambda_c=0, k_a=0`; FA traction (1, 0) ->
    Newton 3rd reaction on vertex = (-1, 0); vertex moves to -x."""
    contour = _make_test_contour(lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    ecm = _make_test_ecm()
    v0_before = contour.vertices_xy_um[0].copy()
    fa = _make_test_fa(position=tuple(v0_before.tolist()), traction=(1.0, 0.0))
    fa_to_vertex = {fa.adhesion_id: 0}

    result = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.0
    )
    v0_after = result.contour.vertices_xy_um[0]
    # Newton 3rd: traction (1, 0) -> reaction (-1, 0) -> vertex moves -x
    assert v0_after[0] < v0_before[0]
    assert abs(v0_after[1] - v0_before[1]) < 1e-12


def test_phase_f_step_zero_traction_and_zero_internal_force_no_movement():
    """Test 6 (Y8 — replaces invalid dt=0 test): lambda_c=0, k_a=0,
    finite positive dt_cell_s, traction (0, 0) -> no motion."""
    contour = _make_test_contour(lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    ecm = _make_test_ecm()
    v0_before = contour.vertices_xy_um[0].copy()
    fa = _make_test_fa(position=tuple(v0_before.tolist()), traction=(0.0, 0.0))
    fa_to_vertex = {fa.adhesion_id: 0}

    result = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.0
    )
    np.testing.assert_allclose(
        result.contour.vertices_xy_um, contour.vertices_xy_um, atol=1e-15
    )


def test_phase_f_step_internal_force_alone_relaxes_per_p1():
    """Test 7: external force zero (no FAs) + nonzero internal -> Phase F
    behaves exactly like P1 step (sister-consistent helper)."""
    from acs.v2.dynamics.active_contour import step as p1_step
    # Setup small nonzero internal force (avoid dt_violation given small
    # vertex spacing → low zeta → large rate_max)
    contour = _make_test_contour(
        lambda_c_nN=0.01, sigma_c_nN_per_um=0.01, k_a_nN_per_um=0.005
    )
    ecm = _make_test_ecm()
    fa_to_vertex: dict[str, int] = {}

    result_pf = step_phase_f_minimal_motility(
        contour, (), ecm, fa_to_vertex, k_active=1.0
    )
    result_p1 = p1_step(contour)
    np.testing.assert_allclose(
        result_pf.contour.vertices_xy_um, result_p1.vertices_xy_um, atol=1e-12
    )


# ---------------------------------------------------------------------------
# HB#4 multipliers diagnostic-only guard (1, Y4 + Y10 wrapper seam)
# ---------------------------------------------------------------------------


def test_phase_f_step_multipliers_diagnostic_only_not_used_for_force_scaling(
    monkeypatch,
):
    """Test 8 (Y10 wrapper seam): monkeypatch step_closed_loop_phase_e_v2
    in Phase F module namespace to return two results with identical
    ECM but neutral vs huge multipliers; vertex displacement IDENTICAL.
    Forward guard against silent magic-number injection (Y4 critical)."""
    contour = _make_test_contour(lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    ecm = _make_test_ecm()
    v0 = contour.vertices_xy_um[0]
    fa = _make_test_fa(position=tuple(v0.tolist()), traction=(1.0, 0.0))
    fa_to_vertex = {fa.adhesion_id: 0}

    # Compute the actual Phase E v2 result to construct fakes with identical ECM
    real_result = step_closed_loop_phase_e_v2(
        (fa,), ecm, dt_s=contour.params.dt_cell_s, k_active=1.0
    )

    from acs.v2.dynamics.ecm_to_fa_bias import ECMToFABiasResult
    from dataclasses import replace as dc_replace

    bias_neutral = real_result.ecm_to_fa_bias
    multipliers_neutral = np.ones((1, 3), dtype=np.float64)
    multipliers_huge = np.full((1, 3), 1e6, dtype=np.float64)
    bias_neutral_fake = dc_replace(
        bias_neutral, multipliers_per_fa=multipliers_neutral
    )
    bias_huge_fake = dc_replace(
        bias_neutral, multipliers_per_fa=multipliers_huge
    )
    fake_neutral = dc_replace(real_result, ecm_to_fa_bias=bias_neutral_fake)
    fake_huge = dc_replace(real_result, ecm_to_fa_bias=bias_huge_fake)

    monkeypatch.setattr(
        pf_mod, "step_closed_loop_phase_e_v2", lambda *a, **k: fake_neutral
    )
    result_neutral = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.0
    )

    monkeypatch.setattr(
        pf_mod, "step_closed_loop_phase_e_v2", lambda *a, **k: fake_huge
    )
    result_huge = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.0
    )

    np.testing.assert_allclose(
        result_neutral.contour.vertices_xy_um,
        result_huge.contour.vertices_xy_um,
        atol=1e-15,
    )


# ---------------------------------------------------------------------------
# FA position + ECM evolution (2)
# ---------------------------------------------------------------------------


def test_phase_f_step_fa_position_follows_attached_vertex():
    """Test 9: after step, FA position matches attached vertex position."""
    contour = _make_test_contour(lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    ecm = _make_test_ecm()
    v0 = contour.vertices_xy_um[0]
    fa = _make_test_fa(position=tuple(v0.tolist()), traction=(1.0, 0.0))
    fa_to_vertex = {fa.adhesion_id: 0}

    result = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.0
    )
    new_vertex = result.contour.vertices_xy_um[0]
    new_fa_position = result.adhesions[0].position_um_xy
    assert new_fa_position[0] == pytest.approx(float(new_vertex[0]), abs=1e-15)
    assert new_fa_position[1] == pytest.approx(float(new_vertex[1]), abs=1e-15)


def test_phase_f_step_phase_e_v2_invoked_and_ecm_carried():
    """Test 10 (Y9 invariant): Phase E v2 result populated + identity
    invariant inherited."""
    contour = _make_test_contour()
    ecm = _make_test_ecm()
    v0 = contour.vertices_xy_um[0]
    fa = _make_test_fa(position=tuple(v0.tolist()))
    fa_to_vertex = {fa.adhesion_id: 0}

    result = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.0
    )
    assert result.phase_e_v2 is not None
    assert isinstance(result.phase_e_v2, PhaseEStepResult)
    assert result.updated_ecm is result.phase_e_v2.updated_ecm


# ---------------------------------------------------------------------------
# Fixture-specific evolution (1, Y9 required per Codex bookkeeping)
# ---------------------------------------------------------------------------


def test_phase_f_step_ecm_orientation_evolves_under_nonzero_traction_fixture():
    """Test 11 (Y9): nonzero traction + isotropic 0.5*I IC + k_active=1.5
    + dt_cell_s=1e-3 -> ECM orientation evolves at float64 level."""
    contour = _make_test_contour(dt_cell_s=1e-3)
    ecm = _make_test_ecm(isotropic_value=0.5)
    v0 = contour.vertices_xy_um[0]
    fa = _make_test_fa(position=tuple(v0.tolist()), traction=(-1.0, 0.0))
    fa_to_vertex = {fa.adhesion_id: 0}

    result = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.5
    )
    delta = (
        np.asarray(result.updated_ecm.orientation_tensor, dtype=np.float64)
        - np.asarray(ecm.orientation_tensor, dtype=np.float64)
    )
    assert float(np.abs(delta).max()) > 0.0


# ---------------------------------------------------------------------------
# Exports (1)
# ---------------------------------------------------------------------------


def test_phase_f_step_exports_through_both_init():
    """Test 12: PhaseFStepResult + step_phase_f_minimal_motility 둘 다 surface
    importable; private helpers NOT exported."""
    expected = {"PhaseFStepResult", "step_phase_f_minimal_motility"}
    forbidden = {
        "_step_active_contour_with_external_force",
        "_validate_fa_to_vertex_index",
    }
    for sym in expected:
        assert hasattr(v2_pkg, sym)
        assert hasattr(v2_dynamics, sym)
        assert getattr(v2_pkg, sym) is getattr(pf_mod, sym)
        assert getattr(v2_dynamics, sym) is getattr(pf_mod, sym)
        assert sym in v2_pkg.__all__
        assert sym in v2_dynamics.__all__
    for sym in forbidden:
        assert sym not in v2_pkg.__all__
        assert sym not in v2_dynamics.__all__


# ---------------------------------------------------------------------------
# PI acceptance tightening (Codex id=1953, PI id=1949 rigid translation 챌린지)
# ---------------------------------------------------------------------------


def test_phase_f_step_attached_vertex_displacement_exceeds_non_attached():
    """Test 13 (Codex id=1953 acceptance #2): attached vertex의 displacement가
    non-attached vertex의 displacement보다 큼 — local pulling/deformation 시연,
    rigid translation 아님."""
    contour = _make_test_contour(
        n_vertices=12,
        radius_um=1.0,
        lambda_c_nN=0.0,  # cortex 약하게 → attached vertex가 유의미하게 움직임
        sigma_c_nN_per_um=0.0,
        k_a_nN_per_um=0.0,
    )
    ecm = _make_test_ecm()
    v0 = contour.vertices_xy_um[0]
    # Single FA at vertex 0; non-attached vertices = vertices 1-11
    fa = _make_test_fa(
        adhesion_id="fa-attached",
        position=tuple(v0.tolist()),
        traction=(1.0, 0.0),
    )
    fa_to_vertex = {fa.adhesion_id: 0}

    result = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.0
    )
    displacements = np.linalg.norm(
        result.contour.vertices_xy_um - contour.vertices_xy_um, axis=1
    )
    attached_disp = displacements[0]
    non_attached_disps = displacements[1:]
    # Attached vertex displacement이 모든 non-attached vertex displacement보다 큼
    assert attached_disp > 0.0
    assert attached_disp > float(non_attached_disps.max()) + 1e-15


def test_phase_f_step_contour_shape_changes_measurably_under_external_force():
    """Test 14 (Codex id=1953 acceptance #3): external force 적용 후 contour
    perimeter / area 가 비영하게 변화 — local deformation 시연."""
    contour = _make_test_contour(
        n_vertices=16,
        radius_um=1.0,
        lambda_c_nN=0.0,
        sigma_c_nN_per_um=0.0,
        k_a_nN_per_um=0.0,
    )
    ecm = _make_test_ecm()
    v0 = contour.vertices_xy_um[0]
    # Small traction so single-step displacement does not self-intersect
    # the polygon (n_vertices=16 packs vertices closer; large translation
    # crosses adjacent edges).
    fa = _make_test_fa(position=tuple(v0.tolist()), traction=(0.1, 0.0))
    fa_to_vertex = {fa.adhesion_id: 0}

    result = step_phase_f_minimal_motility(
        contour, (fa,), ecm, fa_to_vertex, k_active=1.0
    )
    perimeter_before = float(contour.perimeter_um())
    perimeter_after = float(result.contour.perimeter_um())
    area_before = float(contour.area_um2())
    area_after = float(result.contour.area_um2())
    # 변화가 비영해야 — local deformation 시연
    assert abs(perimeter_after - perimeter_before) > 0.0
    assert abs(area_after - area_before) > 0.0
