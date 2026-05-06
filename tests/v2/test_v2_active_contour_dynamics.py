from __future__ import annotations

import numpy as np
import pytest

from acs.v2.active_contour import (
    ActiveContourParameters,
    ActiveContourState,
    compute_rate_max,
    ellipse_polygon_vertices,
    regular_polygon_vertices,
)
from acs.v2.dynamics.active_contour import (
    ActiveContourStepError,
    compute_area_forces,
    compute_cortex_forces,
    energy_area,
    energy_cortex,
    step,
)


def _params(**overrides) -> ActiveContourParameters:
    base = dict(
        k_a_nN_per_um=1.0,
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=8,
        lambda_c_nN=1.0,
        xi_line_nN_s_per_um2=1.0,
    )
    base.update(overrides)
    return ActiveContourParameters(**base).validate()


def _regular_state(n: int, radius_um: float, *, target_area_um2=None, **param_overrides):
    if target_area_um2 is None:
        target_area_um2 = 0.5 * n * radius_um * radius_um * np.sin(2.0 * np.pi / n)
    params = _params(n_vertices=n, target_area_um2=target_area_um2, **param_overrides)
    verts = regular_polygon_vertices(n, radius_um=radius_um)
    return ActiveContourState(cell_id="cell-A", vertices_xy_um=verts, params=params)


def test_cortex_forces_shape_and_zero_when_lambda_zero():
    state = _regular_state(8, 1.5, lambda_c_nN=0.0)
    forces = compute_cortex_forces(state)
    assert forces.shape == (8, 2)
    np.testing.assert_array_equal(forces, np.zeros_like(forces))


def test_area_forces_zero_when_k_a_zero():
    state = _regular_state(8, 1.5, k_a_nN_per_um=0.0)
    forces = compute_area_forces(state)
    assert forces.shape == (8, 2)
    np.testing.assert_array_equal(forces, np.zeros_like(forces))


def test_area_forces_zero_when_at_target_area():
    n, radius = 8, 1.5
    target_area = 0.5 * n * radius * radius * np.sin(2.0 * np.pi / n)
    state = _regular_state(n, radius, target_area_um2=target_area)
    forces = compute_area_forces(state)
    np.testing.assert_allclose(forces, np.zeros_like(forces), atol=1e-12)


def test_cortex_force_pulls_inward_on_regular_polygon():
    """Sanity Gate §5: cortex contractile. Force at each vertex points
    toward the polygon centroid for a regular polygon."""

    state = _regular_state(8, 1.5, lambda_c_nN=1.0, k_a_nN_per_um=0.0)
    forces = compute_cortex_forces(state)
    centroid = state.vertices_xy_um.mean(axis=0)
    inward_direction = centroid - state.vertices_xy_um  # vertex → centroid
    # Each force should have positive component along the inward direction.
    dot = np.einsum("ij,ij->i", forces, inward_direction)
    assert np.all(dot > 0.0)


def test_area_force_pushes_outward_when_below_target():
    state = _regular_state(8, 1.0, target_area_um2=10.0, lambda_c_nN=0.0)
    # current area < target_area
    assert state.area_um2() < state.params.target_area_um2
    forces = compute_area_forces(state)
    centroid = state.vertices_xy_um.mean(axis=0)
    outward_direction = state.vertices_xy_um - centroid
    dot = np.einsum("ij,ij->i", forces, outward_direction)
    assert np.all(dot > 0.0)


def test_area_force_pulls_inward_when_above_target():
    state = _regular_state(8, 2.0, target_area_um2=1.0, lambda_c_nN=0.0)
    assert state.area_um2() > state.params.target_area_um2
    forces = compute_area_forces(state)
    centroid = state.vertices_xy_um.mean(axis=0)
    outward_direction = state.vertices_xy_um - centroid
    dot = np.einsum("ij,ij->i", forces, outward_direction)
    assert np.all(dot < 0.0)


def test_step_zero_force_baseline_no_motion():
    """Test 1: lambda_c=0 and k_a=0 → step returns identical vertices."""

    state = _regular_state(8, 1.5, lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    new_state = step(state)
    np.testing.assert_array_equal(new_state.vertices_xy_um, state.vertices_xy_um)
    assert new_state.step_count == state.step_count + 1


def test_step_at_target_area_only_no_motion():
    """Test 2 (equilibrium): lambda_c=0, k_a>0, A=A_0 → no motion."""

    n, radius = 8, 1.5
    target_area = 0.5 * n * radius * radius * np.sin(2.0 * np.pi / n)
    state = _regular_state(n, radius, target_area_um2=target_area, lambda_c_nN=0.0)
    new_state = step(state)
    np.testing.assert_allclose(
        new_state.vertices_xy_um, state.vertices_xy_um, atol=1e-12
    )


def test_step_returns_new_state_does_not_mutate_input():
    state = _regular_state(8, 1.5, lambda_c_nN=1.0, k_a_nN_per_um=0.0)
    original = state.vertices_xy_um.copy()
    new_state = step(state)
    # Original state is unchanged.
    np.testing.assert_array_equal(state.vertices_xy_um, original)
    # New state is a distinct object with new vertices.
    assert new_state is not state
    assert not np.shares_memory(new_state.vertices_xy_um, state.vertices_xy_um)
    assert new_state.step_count == state.step_count + 1


def test_step_increments_step_count():
    state = _regular_state(8, 1.5, lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    s1 = step(state)
    s2 = step(s1)
    assert state.step_count == 0
    assert s1.step_count == 1
    assert s2.step_count == 2


def test_step_post_step_measurement_boundary_validation():
    state = _regular_state(8, 1.5, lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    # A clean step must export a canonical MeasurementBoundary without
    # auto_orient (regular CCW polygon stays CCW under zero force).
    new_state = step(state)
    boundary = new_state.to_measurement_boundary()
    assert not boundary.orientation_was_reversed


def test_step_aborts_on_dt_violation():
    """When dt is large enough to violate dt·rate_max ≤ 0.5, step raises."""

    state = _regular_state(8, 1.5, dt_cell_s=1.0)  # huge dt for unit drag
    info = compute_rate_max(state)
    assert info["margin_ratio"] > 1.0
    with pytest.raises(ActiveContourStepError) as exc_info:
        step(state)
    assert exc_info.value.failure_kind == "dt_violation"


def test_step_does_not_auto_shrink_dt():
    state = _regular_state(8, 1.5, dt_cell_s=1.0)
    with pytest.raises(ActiveContourStepError):
        step(state)
    # The original state's dt is unchanged (no auto-shrink semantics).
    assert state.params.dt_cell_s == 1.0


def test_cortex_only_energy_decreases_under_step():
    """Test 3 monotonicity sketch: cortex-only step, regular polygon
    shrinks self-similarly and E_c decreases."""

    n, radius = 16, 2.0
    state = _regular_state(
        n, radius, lambda_c_nN=1.0, k_a_nN_per_um=0.0, dt_cell_s=1e-3
    )
    e0 = energy_cortex(state)
    new_state = step(state)
    e1 = energy_cortex(new_state)
    assert e1 < e0
    # Self-similar shrinkage: the new polygon's centroid is the same.
    np.testing.assert_allclose(
        new_state.vertices_xy_um.mean(axis=0),
        state.vertices_xy_um.mean(axis=0),
        atol=1e-9,
    )


def test_area_only_energy_decreases_under_step_when_above_target():
    n, radius = 12, 2.0
    state = _regular_state(
        n,
        radius,
        target_area_um2=2.0,
        lambda_c_nN=0.0,
        k_a_nN_per_um=1.0,
        dt_cell_s=1e-3,
    )
    e0 = energy_area(state)
    new_state = step(state)
    e1 = energy_area(new_state)
    assert e1 < e0


def test_area_only_energy_decreases_under_step_when_below_target():
    n, radius = 12, 0.5
    state = _regular_state(
        n,
        radius,
        target_area_um2=20.0,
        lambda_c_nN=0.0,
        k_a_nN_per_um=1.0,
        dt_cell_s=1e-3,
    )
    e0 = energy_area(state)
    new_state = step(state)
    e1 = energy_area(new_state)
    assert e1 < e0


def test_ellipse_initial_step_runs_and_validates():
    """Test 4 diagnostic smoke: an ellipse can take at least one step
    without contract/measurement-boundary failure."""

    n = 64
    semi_a, semi_b = 1.5, 1.0
    target_area = float(np.pi * semi_a * semi_b)
    params = _params(
        n_vertices=n,
        target_area_um2=target_area,
        lambda_c_nN=1.0,
        k_a_nN_per_um=1.0,
        xi_line_nN_s_per_um2=1.0,
        dt_cell_s=1e-3,
    )
    verts = ellipse_polygon_vertices(n, semi_a_um=semi_a, semi_b_um=semi_b)
    state = ActiveContourState(cell_id="cell-A", vertices_xy_um=verts, params=params)
    new_state = step(state)
    boundary = new_state.to_measurement_boundary()
    assert boundary.projected_area_um2() > 0.0
    assert new_state.step_count == 1
