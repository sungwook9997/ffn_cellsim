"""Sanity Gate test harness for P1 alpha active contour.

Implements ``docs/v2_p1_active_contour_sanity_gate.md`` §4 Tests 1-4
on top of :mod:`acs.v2.active_contour_harness`. Tests 1, 2, 3 are
gating: their pass criteria are asserted. Test 4 is diagnostic
(non-gating per gate §4 amendment): its values are recorded but no
tolerance comparison gates the run.
"""

from __future__ import annotations

import json
import os

import h5py
import numpy as np
import pytest

from acs.v2.active_contour import (
    ActiveContourParameters,
    ActiveContourState,
    ellipse_polygon_vertices,
    regular_polygon_vertices,
)
from acs.v2.active_contour_harness import (
    StepDiagnostics,
    compute_finite_n_residual,
    equilibrium_radius_from_cubic,
    run_active_contour_test,
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


def _energy_non_increasing(history: list[StepDiagnostics], *, atol: float = 1e-12) -> bool:
    """Energy is non-increasing each step; strict decrease only while
    ||F||_2 > 0 (gate §4 wording)."""

    for prev, curr in zip(history, history[1:]):
        if curr.force_norm_nN > 0.0:
            if not (curr.energy_cortex_nNum + curr.energy_area_nNum
                    <= prev.energy_cortex_nNum + prev.energy_area_nNum + atol):
                return False
        else:
            if abs(
                (curr.energy_cortex_nNum + curr.energy_area_nNum)
                - (prev.energy_cortex_nNum + prev.energy_area_nNum)
            ) > atol:
                return False
    return True


def test_gate_test_1_zero_force_no_motion(tmp_path):
    n, radius = 8, 1.5
    params = _params(
        n_vertices=n,
        target_area_um2=10.0,
        lambda_c_nN=0.0,
        k_a_nN_per_um=0.0,
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    run = run_active_contour_test(
        "test_1_zero_force",
        state,
        str(tmp_path / "test_1"),
        n_steps=20,
        frame_interval=10,
    )
    assert run.failure is None
    final = run.final_state
    np.testing.assert_array_equal(final.vertices_xy_um, state.vertices_xy_um)
    assert final.step_count == 20
    # Frame, png, html written at step 0, 10, 20.
    assert len(run.frame_paths) == 3
    assert len(run.png_paths) == 3
    assert len(run.html_paths) == 3
    assert run.diagnostic_plot_path is not None
    assert os.path.exists(run.diagnostic_plot_path)
    assert run.summary_path is not None
    assert os.path.exists(run.summary_path)
    assert "PASS" in open(run.summary_path, encoding="utf-8").read()


def test_gate_test_2_area_only_energy_non_increasing_and_drives_to_target(tmp_path):
    n = 16
    target_area = 4.0  # below current area to force inward push
    radius = 2.0
    params = _params(
        n_vertices=n,
        target_area_um2=target_area,
        lambda_c_nN=0.0,
        k_a_nN_per_um=1.0,
        dt_cell_s=1e-3,
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    run = run_active_contour_test(
        "test_2_area_only",
        state,
        str(tmp_path / "test_2"),
        n_steps=200,
        frame_interval=50,
    )
    assert run.failure is None
    assert _energy_non_increasing(run.history)
    # Area drives toward target (closer to A_0 at end than at start).
    assert abs(run.history[-1].area_um2 - target_area) < abs(
        run.history[0].area_um2 - target_area
    )


def test_gate_test_3_cortex_only_energy_non_increasing_and_self_similar(tmp_path):
    n, radius = 16, 2.0
    params = _params(
        n_vertices=n,
        target_area_um2=10.0,
        lambda_c_nN=1.0,
        k_a_nN_per_um=0.0,
        dt_cell_s=1e-3,
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    run = run_active_contour_test(
        "test_3_cortex_only",
        state,
        str(tmp_path / "test_3"),
        n_steps=100,
        frame_interval=25,
    )
    assert run.failure is None
    assert _energy_non_increasing(run.history)
    # Anisotropy stays near zero (regular polygon stays regular under
    # cortex-only relaxation).
    assert all(d.anisotropy < 1e-9 for d in run.history)
    # Centroid is preserved.
    initial_centroid = state.vertices_xy_um.mean(axis=0)
    final_centroid = run.final_state.vertices_xy_um.mean(axis=0)
    np.testing.assert_allclose(initial_centroid, final_centroid, atol=1e-9)


def test_gate_test_4_coupled_ellipse_diagnostic_records(tmp_path):
    n = 64
    semi_a, semi_b = 1.5, 1.0
    target_area = float(np.pi * semi_a * semi_b)
    params = _params(
        n_vertices=n,
        target_area_um2=target_area,
        lambda_c_nN=0.05,
        k_a_nN_per_um=1.0,
        xi_line_nN_s_per_um2=1.0,
        dt_cell_s=1e-3,
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=ellipse_polygon_vertices(
            n, semi_a_um=semi_a, semi_b_um=semi_b
        ),
        params=params,
    )
    r_star = equilibrium_radius_from_cubic(
        lambda_c_nN=params.lambda_c_resolved_nN,
        k_a_nN_per_um=params.k_a_nN_per_um,
        target_area_um2=target_area,
    )
    reference_state = ActiveContourState(
        cell_id="cell-A-ref",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=r_star),
        params=params,
    )
    finite_n_residual = compute_finite_n_residual(reference_state)

    run = run_active_contour_test(
        "test_4_coupled_ellipse",
        state,
        str(tmp_path / "test_4"),
        n_steps=200,
        frame_interval=50,
        extra_metrics={
            "r_star_um": r_star,
            "finite_n_residual_nN": finite_n_residual,
        },
    )
    # Test 4 is diagnostic-only (gate §4 amendment): we record but do
    # not gate on residual / anisotropy thresholds. We still require the
    # run to make it through the configured n_steps without dt or
    # measurement-boundary contract violation.
    assert run.failure is None
    assert run.history[-1].step_index == 200
    assert _energy_non_increasing(run.history)
    assert "finite_n_residual_nN" in run.extra_metrics
    assert run.extra_metrics["r_star_um"] == r_star


def test_failure_path_writes_failure_report_and_truncated_artifacts(tmp_path):
    """When a step fails (e.g. dt too large for the chosen rate), the
    harness records the failure, truncates frame writes, generates a
    diagnostic plot of what ran, and writes failure_report.md."""

    n, radius = 8, 1.5
    # dt is huge; rate_max enforcement raises on the first step.
    params = _params(n_vertices=n, dt_cell_s=1.0)
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    run = run_active_contour_test(
        "test_failure",
        state,
        str(tmp_path / "fail"),
        n_steps=10,
        frame_interval=1,
    )
    assert run.failure is not None
    assert run.failure["kind"] == "dt_violation"
    assert run.failure["step_index"] == 1
    # Step 0 was emitted before the failure on step 1.
    assert len(run.frame_paths) == 1
    assert run.failure_report_path is not None
    assert os.path.exists(run.failure_report_path)
    assert "failure" in open(run.failure_report_path, encoding="utf-8").read().lower()
    assert "FAIL" in open(run.summary_path, encoding="utf-8").read()


def test_frame_dump_files_are_valid_hdf5(tmp_path):
    n, radius = 8, 1.0
    params = _params(
        n_vertices=n,
        lambda_c_nN=0.0,
        k_a_nN_per_um=0.0,
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    run = run_active_contour_test(
        "test_hdf5_validity",
        state,
        str(tmp_path / "hdf5"),
        n_steps=4,
        frame_interval=2,
    )
    assert run.failure is None
    assert len(run.frame_paths) == 3  # step 0, 2, 4
    for frame_path in run.frame_paths:
        with h5py.File(frame_path, "r") as f:
            assert "schema_version" in f.attrs
            assert "meta" in f
            assert "cells" in f


def test_summary_html_includes_status_table_and_extra_metrics(tmp_path):
    n, radius = 8, 1.0
    params = _params(
        n_vertices=n,
        lambda_c_nN=0.0,
        k_a_nN_per_um=0.0,
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    run = run_active_contour_test(
        "test_summary_metrics",
        state,
        str(tmp_path / "summary"),
        n_steps=2,
        frame_interval=1,
        git_commit_hash="abc1234",
        extra_metrics={"finite_n_residual_nN": 0.0, "custom_marker": "alpha"},
    )
    text = open(run.summary_path, encoding="utf-8").read()
    # All gate §8 status table fields are exposed in the rendered HTML.
    for required in (
        "abc1234",
        "Frame thumbnails",
        "Per-step metrics",
        "wall_clock_s",
        "final_residual_nN",
        "finite_N_residual_nN",
        "final_area_um2",
        "target_area_um2",
        "frame_interval",
        "extra_metrics",
        "custom_marker",
        "alpha",
    ):
        assert required in text, f"summary.html missing {required!r}"


def test_metadata_json_persists_status_payload(tmp_path):
    n, radius = 8, 1.0
    params = _params(
        n_vertices=n,
        lambda_c_nN=0.0,
        k_a_nN_per_um=0.0,
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    run = run_active_contour_test(
        "metadata_persistence",
        state,
        str(tmp_path / "meta"),
        n_steps=4,
        frame_interval=2,
        git_commit_hash="commit123",
        extra_metrics={"finite_n_residual_nN": 0.42, "r_star_um": 1.0},
    )
    assert run.metadata_path is not None
    assert os.path.exists(run.metadata_path)
    with open(run.metadata_path, encoding="utf-8") as fh:
        data = json.load(fh)
    expected_keys = {
        "test_name",
        "status",
        "n_steps",
        "final_step_index",
        "frame_interval",
        "wall_clock_s",
        "git_commit_hash",
        "final_residual_nN",
        "finite_N_residual_nN",
        "final_area_um2",
        "target_area_um2",
        "extra_metrics",
        "failure",
    }
    assert expected_keys.issubset(data.keys())
    assert data["status"] == "PASS"
    assert data["frame_interval"] == 2
    assert data["n_steps"] == 4
    assert data["git_commit_hash"] == "commit123"
    assert data["finite_N_residual_nN"] == pytest.approx(0.42)
    assert data["extra_metrics"]["r_star_um"] == 1.0
    assert data["wall_clock_s"] >= 0.0


def test_reference_r_star_png_generated_when_extra_metrics_have_r_star(tmp_path):
    n = 32
    semi_a, semi_b = 1.5, 1.0
    target_area = float(np.pi * semi_a * semi_b)
    params = _params(
        n_vertices=n,
        target_area_um2=target_area,
        lambda_c_nN=0.05,
        k_a_nN_per_um=1.0,
        xi_line_nN_s_per_um2=1.0,
        dt_cell_s=1e-3,
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=ellipse_polygon_vertices(
            n, semi_a_um=semi_a, semi_b_um=semi_b
        ),
        params=params,
    )
    r_star = equilibrium_radius_from_cubic(
        lambda_c_nN=params.lambda_c_resolved_nN,
        k_a_nN_per_um=params.k_a_nN_per_um,
        target_area_um2=target_area,
    )
    run = run_active_contour_test(
        "test_reference_panel",
        state,
        str(tmp_path / "ref"),
        n_steps=20,
        frame_interval=10,
        extra_metrics={"r_star_um": r_star},
    )
    assert run.failure is None
    assert run.reference_plot_path is not None
    assert os.path.exists(run.reference_plot_path)
    summary_text = open(run.summary_path, encoding="utf-8").read()
    assert "Reference equilibrium" in summary_text


def test_reference_r_star_png_skipped_when_no_r_star_in_extras(tmp_path):
    n, radius = 8, 1.0
    params = _params(
        n_vertices=n, lambda_c_nN=0.0, k_a_nN_per_um=0.0
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    run = run_active_contour_test(
        "test_no_ref",
        state,
        str(tmp_path / "no_ref"),
        n_steps=2,
        frame_interval=1,
    )
    assert run.reference_plot_path is None


def test_run_active_contour_test_rejects_bool_n_steps(tmp_path):
    n, radius = 8, 1.0
    params = _params(n_vertices=n, lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    with pytest.raises(TypeError, match="n_steps"):
        run_active_contour_test(
            "test_bool_steps",
            state,
            str(tmp_path / "bool_steps"),
            n_steps=True,  # type: ignore[arg-type]
        )


def test_run_active_contour_test_rejects_float_n_steps(tmp_path):
    n, radius = 8, 1.0
    params = _params(n_vertices=n, lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    with pytest.raises(TypeError, match="n_steps"):
        run_active_contour_test(
            "test_float_steps",
            state,
            str(tmp_path / "float_steps"),
            n_steps=2.5,  # type: ignore[arg-type]
        )


def test_run_active_contour_test_rejects_zero_frame_interval(tmp_path):
    n, radius = 8, 1.0
    params = _params(n_vertices=n, lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    with pytest.raises(ValueError, match="frame_interval"):
        run_active_contour_test(
            "test_zero_K",
            state,
            str(tmp_path / "zero_K"),
            n_steps=2,
            frame_interval=0,
        )


def test_run_active_contour_test_default_frame_interval_is_100(tmp_path):
    """Gate §8: K default 100, stored in metadata."""

    n, radius = 8, 1.0
    params = _params(n_vertices=n, lambda_c_nN=0.0, k_a_nN_per_um=0.0)
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=radius),
        params=params,
    )
    # n_steps < default K → only step 0 + final step are emitted.
    run = run_active_contour_test(
        "test_default_K",
        state,
        str(tmp_path / "default_K"),
        n_steps=4,
    )
    assert run.frame_interval == 100
    with open(run.metadata_path, encoding="utf-8") as fh:
        data = json.load(fh)
    assert data["frame_interval"] == 100


def test_equilibrium_radius_solves_cubic():
    target_area = 10.0
    r_star = equilibrium_radius_from_cubic(
        lambda_c_nN=0.5,
        k_a_nN_per_um=1.0,
        target_area_um2=target_area,
    )
    # Substitute back into K_A·(πR² - A_0)/A_0 + λ_c/R == 0.
    residual = (
        1.0 * (np.pi * r_star * r_star - target_area) / target_area
        + 0.5 / r_star
    )
    assert abs(residual) < 1e-9
    # Sanity: the equilibrium radius should be near sqrt(A_0/π) for
    # small lambda_c.
    assert r_star > 0.0


def test_compute_finite_n_residual_zero_for_zero_force_state():
    n = 8
    params = _params(
        n_vertices=n,
        lambda_c_nN=0.0,
        k_a_nN_per_um=0.0,
    )
    state = ActiveContourState(
        cell_id="cell-A",
        vertices_xy_um=regular_polygon_vertices(n, radius_um=1.0),
        params=params,
    )
    assert compute_finite_n_residual(state) == 0.0
