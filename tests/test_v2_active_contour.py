from __future__ import annotations

import numpy as np
import pytest

from acs.v2.active_contour import (
    ActiveContourParameters,
    ActiveContourParametersError,
    ActiveContourState,
    _AGREEMENT_TOLERANCE_RELATIVE,
    _DT_RATE_SAFETY_MARGIN,
    compute_rate_max,
    ellipse_polygon_vertices,
    regular_polygon_vertices,
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
    return ActiveContourParameters(**base)


def test_minimal_params_validate_with_direct_paths():
    params = _params().validate()
    assert params.lambda_c_resolved_nN == 1.0
    assert params.xi_line_resolved_nN_s_per_um2 == 1.0


def test_cortex_underspecified_rejected():
    params = ActiveContourParameters(
        k_a_nN_per_um=1.0,
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=8,
        xi_line_nN_s_per_um2=1.0,
    )
    with pytest.raises(ActiveContourParametersError) as info:
        params.validate()
    assert info.value.failure_kind == "cortex_underspecified"


def test_drag_underspecified_rejected():
    params = ActiveContourParameters(
        k_a_nN_per_um=1.0,
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=8,
        lambda_c_nN=1.0,
    )
    with pytest.raises(ActiveContourParametersError) as info:
        params.validate()
    assert info.value.failure_kind == "drag_underspecified"


def test_cortex_derived_path_resolves():
    params = _params(
        lambda_c_nN=None,
        sigma_c_nN_per_um=0.5,
        effective_height_um=2.0,
    ).validate()
    assert params.lambda_c_resolved_nN == pytest.approx(1.0)


def test_drag_derived_path_resolves():
    params = _params(
        xi_line_nN_s_per_um2=None,
        xi_areal_nN_s_per_um3=0.25,
        effective_height_um=4.0,
    ).validate()
    assert params.xi_line_resolved_nN_s_per_um2 == pytest.approx(1.0)


def test_dual_path_agreement_accepted():
    params = _params(
        lambda_c_nN=1.0,
        sigma_c_nN_per_um=0.5,
        effective_height_um=2.0,
        xi_line_nN_s_per_um2=1.0,
        xi_areal_nN_s_per_um3=0.5,
    ).validate()
    assert params.lambda_c_resolved_nN == 1.0


def test_dual_path_disagreement_rejected_for_cortex():
    params = _params(
        lambda_c_nN=1.0,
        sigma_c_nN_per_um=0.5,
        effective_height_um=2.5,
    )
    with pytest.raises(ActiveContourParametersError) as info:
        params.validate()
    assert info.value.failure_kind == "cortex_overspecified_disagreement"


def test_dual_path_disagreement_rejected_for_drag():
    params = _params(
        xi_line_nN_s_per_um2=1.0,
        xi_areal_nN_s_per_um3=0.4,
        effective_height_um=3.0,
    )
    with pytest.raises(ActiveContourParametersError) as info:
        params.validate()
    assert info.value.failure_kind == "drag_overspecified_disagreement"


def test_agreement_tolerance_within_named_relative():
    params = _params(
        lambda_c_nN=1.0,
        sigma_c_nN_per_um=0.5,
        effective_height_um=2.0
        * (1.0 + 0.5 * _AGREEMENT_TOLERANCE_RELATIVE),
    ).validate()
    assert params.lambda_c_resolved_nN == pytest.approx(1.0, rel=1e-8)


def test_required_scalars_reject_negative_and_nan():
    for field in ("target_area_um2", "dt_cell_s"):
        params = _params(**{field: -1.0})
        with pytest.raises(ActiveContourParametersError):
            params.validate()
        params_nan = _params(**{field: float("nan")})
        with pytest.raises(ActiveContourParametersError):
            params_nan.validate()


def test_k_a_rejects_negative_and_nan_but_accepts_zero():
    params_neg = _params(k_a_nN_per_um=-1.0)
    with pytest.raises(ActiveContourParametersError):
        params_neg.validate()
    params_nan = _params(k_a_nN_per_um=float("nan"))
    with pytest.raises(ActiveContourParametersError):
        params_nan.validate()
    params_zero = _params(k_a_nN_per_um=0.0)
    params_zero.validate()


def test_lambda_c_rejects_negative_and_nan_but_accepts_zero():
    params_neg = _params(lambda_c_nN=-1.0)
    with pytest.raises(ActiveContourParametersError):
        params_neg.validate()
    params_nan = _params(lambda_c_nN=float("nan"))
    with pytest.raises(ActiveContourParametersError):
        params_nan.validate()
    params_zero = _params(lambda_c_nN=0.0)
    params_zero.validate()


def test_n_vertices_must_be_int_at_least_three():
    with pytest.raises(ActiveContourParametersError):
        _params(n_vertices=2).validate()
    with pytest.raises(ActiveContourParametersError):
        _params(n_vertices=True).validate()  # type: ignore[arg-type]
    with pytest.raises(ActiveContourParametersError):
        _params(n_vertices=4.5).validate()  # type: ignore[arg-type]


def test_state_construction_validates_vertices_shape():
    params = _params(n_vertices=6).validate()
    verts = regular_polygon_vertices(6, radius_um=1.0)
    state = ActiveContourState(cell_id="cell-A", vertices_xy_um=verts, params=params)
    assert state.vertices_xy_um.shape == (6, 2)
    assert state.vertices_xy_um.flags.writeable


def test_state_construction_rejects_count_mismatch():
    params = _params(n_vertices=6).validate()
    with pytest.raises(ActiveContourParametersError) as info:
        ActiveContourState(
            cell_id="cell-A",
            vertices_xy_um=regular_polygon_vertices(8, radius_um=1.0),
            params=params,
        )
    assert info.value.failure_kind == "vertices_count_mismatch"


def test_state_to_measurement_boundary_round_trip():
    params = _params(n_vertices=8).validate()
    verts = regular_polygon_vertices(8, radius_um=1.5)
    state = ActiveContourState(cell_id="cell-A", vertices_xy_um=verts, params=params)
    boundary = state.to_measurement_boundary()
    assert boundary.coordinate_convention == "world_um_y_up"
    assert boundary.source_modality == "active_contour"
    assert boundary.projected_area_um2() == pytest.approx(state.area_um2())


def test_state_area_perimeter_match_shoelace():
    params = _params(n_vertices=6).validate()
    radius = 2.0
    verts = regular_polygon_vertices(6, radius_um=radius)
    state = ActiveContourState(cell_id="cell-A", vertices_xy_um=verts, params=params)
    expected_area = 0.5 * 6 * radius * radius * np.sin(2.0 * np.pi / 6)
    expected_perimeter = 6 * 2.0 * radius * np.sin(np.pi / 6)
    assert state.area_um2() == pytest.approx(expected_area, rel=1e-12)
    assert state.perimeter_um() == pytest.approx(expected_perimeter, rel=1e-12)


def test_compute_rate_max_metadata_keys():
    params = _params(n_vertices=8).validate()
    verts = regular_polygon_vertices(8, radius_um=1.5)
    state = ActiveContourState(cell_id="cell-A", vertices_xy_um=verts, params=params)
    info = compute_rate_max(state)
    expected = {
        "rate_max_s_inv",
        "dt_cell_s",
        "dt_rate_product",
        "margin_ratio",
        "jacobian_row_cortex_nN_per_um",
        "jacobian_row_area_nN_per_um",
        "jacobian_row_total_nN_per_um",
        "zeta_nN_s_per_um",
        "rate_per_vertex_s_inv",
        "min_edge_length_um",
    }
    assert expected.issubset(info.keys())
    assert info["rate_max_s_inv"] > 0.0
    assert info["dt_rate_product"] == pytest.approx(
        info["dt_cell_s"] * info["rate_max_s_inv"]
    )
    assert info["margin_ratio"] == pytest.approx(
        info["dt_rate_product"] / _DT_RATE_SAFETY_MARGIN
    )
    np.testing.assert_allclose(
        info["jacobian_row_total_nN_per_um"],
        info["jacobian_row_cortex_nN_per_um"]
        + info["jacobian_row_area_nN_per_um"],
    )
    np.testing.assert_allclose(
        info["rate_per_vertex_s_inv"],
        info["jacobian_row_total_nN_per_um"] / info["zeta_nN_s_per_um"],
    )


def test_compute_rate_max_zero_force_at_target_area_only_cortex():
    radius = 2.0
    verts = regular_polygon_vertices(8, radius_um=radius)
    expected_area = 0.5 * 8 * radius * radius * np.sin(2.0 * np.pi / 8)
    params = ActiveContourParameters(
        k_a_nN_per_um=1.0,
        target_area_um2=expected_area,
        dt_cell_s=1e-4,
        n_vertices=8,
        lambda_c_nN=0.5,
        xi_line_nN_s_per_um2=1.0,
    ).validate()
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=verts,
        params=params,
    )
    info = compute_rate_max(state)
    # Area force vanishes at the target area, so jacobian_row_area only
    # has the (K_A/A_0) |g_i| sum term. Rate per vertex must be finite.
    assert np.all(np.isfinite(info["rate_per_vertex_s_inv"]))


def test_zero_force_baseline_lambda_c_zero_k_a_zero_validate():
    """Lock §6 Test 1 zero-force baseline: lambda_c=0 AND k_a=0 must validate."""

    params = ActiveContourParameters(
        k_a_nN_per_um=0.0,
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=8,
        lambda_c_nN=0.0,
        xi_line_nN_s_per_um2=1.0,
    ).validate()
    assert params.lambda_c_resolved_nN == 0.0
    assert params.k_a_nN_per_um == 0.0


def test_area_only_lambda_c_zero_validate():
    """Lock §6 Test 2 area-only: lambda_c=0, k_a>0 must validate."""

    params = ActiveContourParameters(
        k_a_nN_per_um=1.0,
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=8,
        lambda_c_nN=0.0,
        xi_line_nN_s_per_um2=1.0,
    ).validate()
    assert params.lambda_c_resolved_nN == 0.0


def test_cortex_only_k_a_zero_validate():
    """Lock §6 Test 3 cortex-only: lambda_c>0, k_a=0 must validate."""

    params = ActiveContourParameters(
        k_a_nN_per_um=0.0,
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=8,
        lambda_c_nN=1.0,
        xi_line_nN_s_per_um2=1.0,
    ).validate()
    assert params.k_a_nN_per_um == 0.0
    assert params.lambda_c_resolved_nN == 1.0


def test_zero_force_compute_rate_max_zero_rate():
    """When both coefficients are exactly zero, rate_max is zero; force-Jacobian rows are zero."""

    params = ActiveContourParameters(
        k_a_nN_per_um=0.0,
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=8,
        lambda_c_nN=0.0,
        xi_line_nN_s_per_um2=1.0,
    ).validate()
    verts = regular_polygon_vertices(8, radius_um=1.5)
    state = ActiveContourState(cell_id="c", vertices_xy_um=verts, params=params)
    info = compute_rate_max(state)
    assert info["rate_max_s_inv"] == 0.0
    assert np.all(info["jacobian_row_cortex_nN_per_um"] == 0.0)
    assert np.all(info["jacobian_row_area_nN_per_um"] == 0.0)
    assert info["margin_ratio"] == 0.0


def test_drag_must_be_strictly_positive_even_when_force_coefficients_zero():
    """Zero drag is rejected because overdamped EOM is undefined; force coefficients can be zero."""

    params_zero_drag = ActiveContourParameters(
        k_a_nN_per_um=0.0,
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=8,
        lambda_c_nN=0.0,
        xi_line_nN_s_per_um2=0.0,
    )
    with pytest.raises(ActiveContourParametersError) as info:
        params_zero_drag.validate()
    assert info.value.failure_kind == "drag_underspecified"


def test_state_construction_calls_params_validate():
    """State construction rejects an unvalidated/invalid params object so it cannot
    silently propagate into compute_rate_max."""

    bad_params = ActiveContourParameters(
        k_a_nN_per_um=-1.0,  # invalid: rejected by validate()
        target_area_um2=10.0,
        dt_cell_s=1e-4,
        n_vertices=8,
        lambda_c_nN=1.0,
        xi_line_nN_s_per_um2=1.0,
    )
    verts = regular_polygon_vertices(8, radius_um=1.0)
    with pytest.raises(ActiveContourParametersError):
        ActiveContourState(cell_id="c", vertices_xy_um=verts, params=bad_params)


def test_jacobian_row_cortex_grows_when_edges_shrink():
    params = _params(n_vertices=6).validate()
    big = regular_polygon_vertices(6, radius_um=2.0)
    small = regular_polygon_vertices(6, radius_um=0.5)
    info_big = compute_rate_max(
        ActiveContourState(cell_id="c", vertices_xy_um=big, params=params)
    )
    info_small = compute_rate_max(
        ActiveContourState(cell_id="c", vertices_xy_um=small, params=params)
    )
    assert (
        info_small["jacobian_row_cortex_nN_per_um"].max()
        > info_big["jacobian_row_cortex_nN_per_um"].max()
    )


def test_compute_rate_max_dt_safety_margin_marker():
    params = _params(n_vertices=6, dt_cell_s=1e-8).validate()
    verts = regular_polygon_vertices(6, radius_um=2.0)
    state = ActiveContourState(cell_id="c", vertices_xy_um=verts, params=params)
    info = compute_rate_max(state)
    assert info["margin_ratio"] < 1.0


def test_ellipse_polygon_vertices_aspect_ratio():
    verts = ellipse_polygon_vertices(64, semi_a_um=1.5, semi_b_um=1.0)
    assert verts.shape == (64, 2)
    extent_x = verts[:, 0].max() - verts[:, 0].min()
    extent_y = verts[:, 1].max() - verts[:, 1].min()
    assert extent_x == pytest.approx(3.0)
    assert extent_y == pytest.approx(2.0)


def test_regular_polygon_vertices_ccw_canonical():
    verts = regular_polygon_vertices(8, radius_um=1.0)
    boundary_state = ActiveContourState(
        cell_id="c", vertices_xy_um=verts, params=_params(n_vertices=8).validate()
    )
    # CCW means positive shoelace.
    assert boundary_state.area_um2() > 0.0
    # Round-trip through MeasurementBoundary must succeed without
    # auto_orient (i.e. already canonical CCW).
    boundary = boundary_state.to_measurement_boundary()
    assert not boundary.orientation_was_reversed


def test_state_rejects_non_finite_vertices():
    params = _params(n_vertices=4).validate()
    bad = regular_polygon_vertices(4, radius_um=1.0)
    bad[0, 0] = float("nan")
    with pytest.raises(ActiveContourParametersError) as info:
        ActiveContourState(cell_id="c", vertices_xy_um=bad, params=params)
    assert info.value.failure_kind == "non_finite_vertices"
