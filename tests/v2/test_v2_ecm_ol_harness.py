from __future__ import annotations

import json
import os

import h5py
import numpy as np
import pytest

from acs.v2.ecm_open_loop_harness import (
    EcmOlScenario,
    make_default_ecm,
    run_ecm_ol_scenario,
)


def _zero_traction_factory(shape):
    def _factory(step_index: int) -> np.ndarray:  # noqa: ARG001
        return np.zeros(shape, dtype=np.float64)
    return _factory


def _uniform_traction_factory(shape, value: float):
    def _factory(step_index: int) -> np.ndarray:  # noqa: ARG001
        return np.full(shape, value, dtype=np.float64)
    return _factory


def _ramp_traction_factory(shape, low: float, high: float):
    def _factory(step_index: int) -> np.ndarray:  # noqa: ARG001
        return np.linspace(low, high, shape[0] * shape[1], dtype=np.float64).reshape(shape)
    return _factory


def _step_function_traction_factory(shape, value: float, switch_step: int):
    def _factory(step_index: int) -> np.ndarray:
        if step_index < switch_step:
            return np.zeros(shape, dtype=np.float64)
        return np.full(shape, value, dtype=np.float64)
    return _factory


def test_zero_traction_scenario_no_accumulation_change(tmp_path):
    ecm = make_default_ecm(nx=4, ny=4)
    scenario = EcmOlScenario(
        name="zero",
        initial_ecm=ecm,
        traction_factory=_zero_traction_factory((4, 4)),
        n_steps=5,
        dt_s=0.1,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "zero"), frame_interval=2)
    assert run.failure is None
    assert run.final_ecm.accumulated_traction_nNs_per_um2.max() == 0.0
    # 6 step entries (0..5)
    assert len(run.history) == 6
    # Frames at step 0, 2, 4, 5 = 4 frames
    assert len(run.frame_paths) == 4


def test_uniform_traction_scenario_linear_accumulation(tmp_path):
    ecm = make_default_ecm(nx=3, ny=3)
    scenario = EcmOlScenario(
        name="uniform",
        initial_ecm=ecm,
        traction_factory=_uniform_traction_factory((3, 3), 0.5),
        n_steps=10,
        dt_s=0.1,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "uniform"), frame_interval=5)
    assert run.failure is None
    expected_final = 0.5 * 0.1 * 10
    np.testing.assert_allclose(
        run.final_ecm.accumulated_traction_nNs_per_um2,
        np.full((3, 3), expected_final),
    )


def test_spatial_ramp_scenario_preserves_pattern(tmp_path):
    ecm = make_default_ecm(nx=4, ny=2)
    scenario = EcmOlScenario(
        name="ramp",
        initial_ecm=ecm,
        traction_factory=_ramp_traction_factory((4, 2), 0.0, 0.7),
        n_steps=4,
        dt_s=0.25,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "ramp"), frame_interval=2)
    assert run.failure is None
    # 4 steps × 0.25 s × ramp pattern.
    expected_pattern = np.linspace(0.0, 0.7, 8).reshape(4, 2)
    np.testing.assert_allclose(
        run.final_ecm.accumulated_traction_nNs_per_um2,
        expected_pattern * 0.25 * 4,
    )


def test_step_function_scenario_late_activation(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="step_fn",
        initial_ecm=ecm,
        traction_factory=_step_function_traction_factory((2, 2), 0.4, switch_step=6),
        n_steps=10,
        dt_s=0.1,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "step"), frame_interval=2)
    assert run.failure is None
    # Steps 6..10 active, each contributing 0.4 * 0.1 = 0.04. 5 active steps → 0.20.
    expected_final = 0.4 * 0.1 * 5
    np.testing.assert_allclose(
        run.final_ecm.accumulated_traction_nNs_per_um2,
        np.full((2, 2), expected_final),
    )


def test_metadata_json_persists_status_payload(tmp_path):
    ecm = make_default_ecm(nx=3, ny=3)
    scenario = EcmOlScenario(
        name="meta",
        initial_ecm=ecm,
        traction_factory=_uniform_traction_factory((3, 3), 0.2),
        n_steps=4,
        dt_s=0.5,
    )
    run = run_ecm_ol_scenario(
        scenario,
        str(tmp_path / "meta"),
        frame_interval=2,
        git_commit_hash="abc1234",
    )
    assert run.metadata_path is not None
    with open(run.metadata_path) as fh:
        data = json.load(fh)
    expected_keys = {
        "scenario_name",
        "status",
        "n_steps",
        "final_step_index",
        "frame_interval",
        "dt_s",
        "wall_clock_s",
        "git_commit_hash",
        "final_accumulated_max_nNs_per_um2",
        "final_accumulated_mean_nNs_per_um2",
        "failure",
    }
    assert expected_keys.issubset(data.keys())
    assert data["status"] == "PASS"
    assert data["scenario_name"] == "meta"
    assert data["git_commit_hash"] == "abc1234"
    assert data["dt_s"] == 0.5


def test_summary_html_includes_status_table_and_thumbnails(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="summary_test",
        initial_ecm=ecm,
        traction_factory=_uniform_traction_factory((2, 2), 0.3),
        n_steps=3,
        dt_s=0.5,
    )
    run = run_ecm_ol_scenario(
        scenario, str(tmp_path / "summary"), frame_interval=1, git_commit_hash="xyz9876"
    )
    text = open(run.summary_path).read()
    for required in (
        "scenario_name",
        "summary_test",
        "Frame thumbnails",
        "Per-step metrics",
        "wall_clock_s",
        "final_accumulated_max_nNs_per_um2",
        "xyz9876",
    ):
        assert required in text


def test_diagnostic_plot_generated(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="diag",
        initial_ecm=ecm,
        traction_factory=_uniform_traction_factory((2, 2), 0.1),
        n_steps=4,
        dt_s=0.25,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "diag"), frame_interval=2)
    assert run.diagnostic_plot_path is not None
    assert os.path.exists(run.diagnostic_plot_path)


def test_frame_dump_files_are_valid_hdf5(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="hdf5",
        initial_ecm=ecm,
        traction_factory=_uniform_traction_factory((2, 2), 0.5),
        n_steps=2,
        dt_s=0.1,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "hdf5"), frame_interval=1)
    assert run.failure is None
    assert len(run.frame_paths) >= 2
    for frame_path in run.frame_paths:
        with h5py.File(frame_path, "r") as f:
            assert "schema_version" in f.attrs
            assert "meta" in f
            assert "cells" in f
            assert "ecm" in f


def test_invalid_traction_factory_output_recorded_as_failure(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)

    def _bad_factory(step_index: int) -> np.ndarray:
        if step_index >= 3:
            return np.full((2, 2), -1.0, dtype=np.float64)  # negative → reject
        return np.full((2, 2), 0.1, dtype=np.float64)

    scenario = EcmOlScenario(
        name="bad",
        initial_ecm=ecm,
        traction_factory=_bad_factory,
        n_steps=10,
        dt_s=0.1,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "bad"), frame_interval=1)
    assert run.failure is not None
    assert run.failure["kind"] == "ecm_open_loop_violation"
    assert run.failure["failure_kind"] == "traction_negative"
    assert run.failure["step_index"] == 3
    assert run.failure_report_path is not None
    assert os.path.exists(run.failure_report_path)


def test_input_ecm_is_not_mutated_by_harness(tmp_path):
    ecm = make_default_ecm(nx=3, ny=3)
    snapshot = ecm.accumulated_traction_nNs_per_um2.copy()
    scenario = EcmOlScenario(
        name="immut",
        initial_ecm=ecm,
        traction_factory=_uniform_traction_factory((3, 3), 0.5),
        n_steps=3,
        dt_s=0.1,
    )
    run_ecm_ol_scenario(scenario, str(tmp_path / "immut"), frame_interval=1)
    np.testing.assert_array_equal(ecm.accumulated_traction_nNs_per_um2, snapshot)


def test_run_ecm_ol_scenario_rejects_bool_n_steps(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="bad_steps",
        initial_ecm=ecm,
        traction_factory=_zero_traction_factory((2, 2)),
        n_steps=True,  # type: ignore[arg-type]
        dt_s=0.1,
    )
    with pytest.raises(TypeError, match="n_steps"):
        run_ecm_ol_scenario(scenario, str(tmp_path / "bool_steps"), frame_interval=1)


def test_run_ecm_ol_scenario_rejects_zero_frame_interval(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="zero_K",
        initial_ecm=ecm,
        traction_factory=_zero_traction_factory((2, 2)),
        n_steps=2,
        dt_s=0.1,
    )
    with pytest.raises(ValueError, match="frame_interval"):
        run_ecm_ol_scenario(scenario, str(tmp_path / "zero_K"), frame_interval=0)


def test_default_frame_interval_is_100(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="default_K",
        initial_ecm=ecm,
        traction_factory=_uniform_traction_factory((2, 2), 0.1),
        n_steps=4,
        dt_s=0.1,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "default_K"))
    assert run.frame_interval == 100
    with open(run.metadata_path) as fh:
        data = json.load(fh)
    assert data["frame_interval"] == 100


# ---------------------------------------------------------------------------
# Multi-channel scenario factory helpers (caller-supplied test inputs only)
# ---------------------------------------------------------------------------


def _zero_grid_factory(shape):
    def _factory(step_index: int) -> np.ndarray:  # noqa: ARG001
        return np.zeros(shape, dtype=np.float64)
    return _factory


def _uniform_grid_factory(shape, value: float):
    def _factory(step_index: int) -> np.ndarray:  # noqa: ARG001
        return np.full(shape, value, dtype=np.float64)
    return _factory


def _uniform_diagonal_orientation_rate_factory(shape, diag_value: float):
    """Symmetric diagonal-only orientation rate; choose ``diag_value``
    negative when the initial ECM has identity orientation (T_00=T_11=1.0)
    so post-state stays inside ``|T_ij| <= 1.0``."""

    def _factory(step_index: int) -> np.ndarray:  # noqa: ARG001
        rate = np.zeros((*shape, 2, 2), dtype=np.float64)
        rate[..., 0, 0] = diag_value
        rate[..., 1, 1] = diag_value
        return rate
    return _factory


# ---------------------------------------------------------------------------
# Per-channel single-channel scenarios (Codex zero/no-op + per-channel paths)
# ---------------------------------------------------------------------------


def test_stiffness_only_scenario_updates_stiffness_and_preserves_other_fields(tmp_path):
    ecm = make_default_ecm(nx=3, ny=3)
    scenario = EcmOlScenario(
        name="stiffness_only",
        initial_ecm=ecm,
        n_steps=4,
        dt_s=0.5,
        stiffness_rate_factory=_uniform_grid_factory((3, 3), 0.2),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "stiffness_only"), frame_interval=2)
    assert run.failure is None
    expected_final = 0.2 * 0.5 * 4
    np.testing.assert_allclose(
        run.final_ecm.stiffness_kpa, np.full((3, 3), expected_final)
    )
    np.testing.assert_array_equal(
        run.final_ecm.accumulated_traction_nNs_per_um2,
        np.zeros((3, 3)),
    )
    np.testing.assert_array_equal(run.final_ecm.ligand_density, np.zeros((3, 3)))
    np.testing.assert_array_equal(run.final_ecm.fiber_density, np.zeros((3, 3)))


def test_density_only_scenario_updates_both_density_fields(tmp_path):
    ecm = make_default_ecm(nx=3, ny=3)
    scenario = EcmOlScenario(
        name="density_only",
        initial_ecm=ecm,
        n_steps=4,
        dt_s=0.5,
        ligand_density_rate_factory=_uniform_grid_factory((3, 3), 0.1),
        fiber_density_rate_factory=_uniform_grid_factory((3, 3), 0.05),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "density_only"), frame_interval=2)
    assert run.failure is None
    np.testing.assert_allclose(
        run.final_ecm.ligand_density, np.full((3, 3), 0.1 * 0.5 * 4)
    )
    np.testing.assert_allclose(
        run.final_ecm.fiber_density, np.full((3, 3), 0.05 * 0.5 * 4)
    )
    np.testing.assert_array_equal(
        run.final_ecm.stiffness_kpa, np.zeros((3, 3))
    )


def test_density_only_scenario_with_only_ligand_factory_leaves_fiber_unchanged(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="ligand_only_subchannel",
        initial_ecm=ecm,
        n_steps=3,
        dt_s=0.1,
        ligand_density_rate_factory=_uniform_grid_factory((2, 2), 0.2),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "ligand_only"), frame_interval=1)
    assert run.failure is None
    np.testing.assert_allclose(
        run.final_ecm.ligand_density, np.full((2, 2), 0.2 * 0.1 * 3)
    )
    np.testing.assert_array_equal(
        run.final_ecm.fiber_density, np.zeros((2, 2))
    )


def test_orientation_only_scenario_updates_orientation(tmp_path):
    ecm = make_default_ecm(nx=3, ny=3)
    scenario = EcmOlScenario(
        name="orientation_only",
        initial_ecm=ecm,
        n_steps=4,
        dt_s=0.5,
        orientation_rate_factory=_uniform_diagonal_orientation_rate_factory((3, 3), -0.1),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "orient_only"), frame_interval=2)
    assert run.failure is None
    expected_diag = 1.0 + (-0.1) * 0.5 * 4
    np.testing.assert_allclose(
        run.final_ecm.orientation_tensor[..., 0, 0],
        np.full((3, 3), expected_diag),
    )
    np.testing.assert_allclose(
        run.final_ecm.orientation_tensor[..., 1, 1],
        np.full((3, 3), expected_diag),
    )
    np.testing.assert_allclose(
        run.final_ecm.orientation_tensor[..., 0, 1],
        np.zeros((3, 3)),
    )
    np.testing.assert_array_equal(
        run.final_ecm.accumulated_traction_nNs_per_um2, np.zeros((3, 3))
    )


def test_zero_stiffness_rate_no_op_change(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="zero_stiff",
        initial_ecm=ecm,
        n_steps=3,
        dt_s=0.1,
        stiffness_rate_factory=_zero_grid_factory((2, 2)),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "zero_stiff"), frame_interval=1)
    assert run.failure is None
    np.testing.assert_array_equal(
        run.final_ecm.stiffness_kpa, np.zeros((2, 2))
    )


def test_zero_density_rates_no_op_change(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="zero_density",
        initial_ecm=ecm,
        n_steps=3,
        dt_s=0.1,
        ligand_density_rate_factory=_zero_grid_factory((2, 2)),
        fiber_density_rate_factory=_zero_grid_factory((2, 2)),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "zero_density"), frame_interval=1)
    assert run.failure is None
    np.testing.assert_array_equal(run.final_ecm.ligand_density, np.zeros((2, 2)))
    np.testing.assert_array_equal(run.final_ecm.fiber_density, np.zeros((2, 2)))


def test_zero_orientation_rate_no_op_change(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    initial_orientation = ecm.orientation_tensor.copy()
    scenario = EcmOlScenario(
        name="zero_orient",
        initial_ecm=ecm,
        n_steps=3,
        dt_s=0.1,
        orientation_rate_factory=lambda _i: np.zeros((2, 2, 2, 2), dtype=np.float64),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "zero_orient"), frame_interval=1)
    assert run.failure is None
    np.testing.assert_array_equal(
        run.final_ecm.orientation_tensor, initial_orientation
    )


# ---------------------------------------------------------------------------
# All-channel smoke + sequential application order
# ---------------------------------------------------------------------------


def test_all_channel_smoke_scenario_updates_every_field(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="all_channels",
        initial_ecm=ecm,
        n_steps=2,
        dt_s=0.1,
        traction_factory=_uniform_grid_factory((2, 2), 0.3),
        stiffness_rate_factory=_uniform_grid_factory((2, 2), 0.5),
        ligand_density_rate_factory=_uniform_grid_factory((2, 2), 0.2),
        fiber_density_rate_factory=_uniform_grid_factory((2, 2), 0.1),
        orientation_rate_factory=_uniform_diagonal_orientation_rate_factory(
            (2, 2), -0.1
        ),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "all"), frame_interval=1)
    assert run.failure is None
    np.testing.assert_allclose(
        run.final_ecm.accumulated_traction_nNs_per_um2,
        np.full((2, 2), 0.3 * 0.1 * 2),
    )
    np.testing.assert_allclose(
        run.final_ecm.stiffness_kpa, np.full((2, 2), 0.5 * 0.1 * 2)
    )
    np.testing.assert_allclose(
        run.final_ecm.ligand_density, np.full((2, 2), 0.2 * 0.1 * 2)
    )
    np.testing.assert_allclose(
        run.final_ecm.fiber_density, np.full((2, 2), 0.1 * 0.1 * 2)
    )
    expected_diag = 1.0 + (-0.1) * 0.1 * 2
    np.testing.assert_allclose(
        run.final_ecm.orientation_tensor[..., 0, 0],
        np.full((2, 2), expected_diag),
    )


# ---------------------------------------------------------------------------
# Failure capture per channel — failure_kind round-trip
# ---------------------------------------------------------------------------


def test_invalid_stiffness_rate_recorded_as_failure(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)

    def _bad_stiffness(step_index: int) -> np.ndarray:
        if step_index >= 3:
            return np.full((2, 2), -10.0, dtype=np.float64)
        return np.full((2, 2), 0.1, dtype=np.float64)

    scenario = EcmOlScenario(
        name="bad_stiff",
        initial_ecm=ecm,
        n_steps=10,
        dt_s=0.1,
        stiffness_rate_factory=_bad_stiffness,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "bad_stiff"), frame_interval=1)
    assert run.failure is not None
    assert run.failure["kind"] == "ecm_open_loop_violation"
    assert run.failure["failure_kind"] == "stiffness_negative_post_update"
    assert run.failure["step_index"] == 3
    assert run.failure_report_path is not None
    assert os.path.exists(run.failure_report_path)


def test_invalid_ligand_density_rate_recorded_as_failure(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)

    def _bad_ligand(step_index: int) -> np.ndarray:
        if step_index >= 3:
            return np.full((2, 2), 1000.0, dtype=np.float64)
        return np.full((2, 2), 0.01, dtype=np.float64)

    scenario = EcmOlScenario(
        name="bad_ligand",
        initial_ecm=ecm,
        n_steps=10,
        dt_s=0.1,
        ligand_density_rate_factory=_bad_ligand,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "bad_ligand"), frame_interval=1)
    assert run.failure is not None
    assert run.failure["kind"] == "ecm_open_loop_violation"
    assert run.failure["failure_kind"] == "ligand_density_exceeds_one_post_update"
    assert run.failure["step_index"] == 3


def test_invalid_fiber_density_rate_recorded_as_failure(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)

    def _bad_fiber(step_index: int) -> np.ndarray:
        if step_index >= 3:
            return np.full((2, 2), -1000.0, dtype=np.float64)
        return np.full((2, 2), 0.01, dtype=np.float64)

    scenario = EcmOlScenario(
        name="bad_fiber",
        initial_ecm=ecm,
        n_steps=10,
        dt_s=0.1,
        fiber_density_rate_factory=_bad_fiber,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "bad_fiber"), frame_interval=1)
    assert run.failure is not None
    assert run.failure["kind"] == "ecm_open_loop_violation"
    assert run.failure["failure_kind"] == "fiber_density_negative_post_update"
    assert run.failure["step_index"] == 3


def test_invalid_orientation_rate_recorded_as_failure(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)

    def _bad_orientation(step_index: int) -> np.ndarray:
        if step_index >= 2:
            rate = np.zeros((2, 2, 2, 2), dtype=np.float64)
            rate[..., 0, 1] = 0.1
            rate[..., 1, 0] = 0.0
            return rate
        return np.zeros((2, 2, 2, 2), dtype=np.float64)

    scenario = EcmOlScenario(
        name="bad_orient",
        initial_ecm=ecm,
        n_steps=5,
        dt_s=0.1,
        orientation_rate_factory=_bad_orientation,
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "bad_orient"), frame_interval=1)
    assert run.failure is not None
    assert run.failure["kind"] == "ecm_open_loop_violation"
    assert run.failure["failure_kind"] == "orientation_rate_asymmetric"
    assert run.failure["step_index"] == 2


# ---------------------------------------------------------------------------
# Metadata + active channels + per-channel field extrema
# ---------------------------------------------------------------------------


def test_metadata_includes_active_channels_subset(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="channels",
        initial_ecm=ecm,
        n_steps=2,
        dt_s=0.1,
        traction_factory=_uniform_grid_factory((2, 2), 0.1),
        stiffness_rate_factory=_uniform_grid_factory((2, 2), 0.1),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "channels"), frame_interval=1)
    with open(run.metadata_path) as fh:
        data = json.load(fh)
    assert "active_channels" in data
    assert data["active_channels"] == ["traction", "stiffness_rate"]
    assert run.active_channels == ["traction", "stiffness_rate"]


def test_metadata_includes_per_channel_field_extrema(tmp_path):
    ecm = make_default_ecm(nx=3, ny=3)
    scenario = EcmOlScenario(
        name="extrema",
        initial_ecm=ecm,
        n_steps=2,
        dt_s=0.1,
        stiffness_rate_factory=_uniform_grid_factory((3, 3), 0.5),
        ligand_density_rate_factory=_uniform_grid_factory((3, 3), 0.2),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "extrema"), frame_interval=1)
    with open(run.metadata_path) as fh:
        data = json.load(fh)
    expected_keys = {
        "final_stiffness_kpa_max",
        "final_stiffness_kpa_min",
        "final_ligand_density_max",
        "final_ligand_density_min",
        "final_fiber_density_max",
        "final_fiber_density_min",
        "final_orientation_tensor_abs_max",
    }
    assert expected_keys.issubset(data.keys())
    expected_stiffness = 0.5 * 0.1 * 2
    expected_ligand = 0.2 * 0.1 * 2
    assert data["final_stiffness_kpa_max"] == pytest.approx(expected_stiffness)
    assert data["final_stiffness_kpa_min"] == pytest.approx(expected_stiffness)
    assert data["final_ligand_density_max"] == pytest.approx(expected_ligand)
    assert data["final_orientation_tensor_abs_max"] == pytest.approx(1.0)


def test_summary_html_includes_active_channels_and_per_channel_keys(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="html_multi",
        initial_ecm=ecm,
        n_steps=2,
        dt_s=0.1,
        stiffness_rate_factory=_uniform_grid_factory((2, 2), 0.1),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "html_multi"), frame_interval=1)
    text = open(run.summary_path).read()
    for required in (
        "active_channels",
        "final_stiffness_kpa_max",
        "final_ligand_density_max",
        "final_orientation_tensor_abs_max",
        "stiffness_rate",
    ):
        assert required in text


def test_diagnostic_plot_generated_for_multi_channel_run(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="plot_multi",
        initial_ecm=ecm,
        n_steps=3,
        dt_s=0.1,
        stiffness_rate_factory=_uniform_grid_factory((2, 2), 0.1),
        ligand_density_rate_factory=_uniform_grid_factory((2, 2), 0.05),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "plot_multi"), frame_interval=1)
    assert run.diagnostic_plot_path is not None
    assert os.path.exists(run.diagnostic_plot_path)


# ---------------------------------------------------------------------------
# Contract: at-least-one-factory + multi-channel input ECM not mutated
# ---------------------------------------------------------------------------


def test_scenario_requires_at_least_one_factory():
    ecm = make_default_ecm(nx=2, ny=2)
    with pytest.raises(ValueError, match="at least one channel factory"):
        EcmOlScenario(
            name="empty",
            initial_ecm=ecm,
            n_steps=2,
            dt_s=0.1,
        )


def test_input_ecm_not_mutated_by_multi_channel_scenario(tmp_path):
    ecm = make_default_ecm(nx=2, ny=2)
    snapshots = {
        "accumulated_traction_nNs_per_um2": ecm.accumulated_traction_nNs_per_um2.copy(),
        "stiffness_kpa": ecm.stiffness_kpa.copy(),
        "ligand_density": ecm.ligand_density.copy(),
        "fiber_density": ecm.fiber_density.copy(),
        "orientation_tensor": ecm.orientation_tensor.copy(),
    }
    scenario = EcmOlScenario(
        name="immut_multi",
        initial_ecm=ecm,
        n_steps=3,
        dt_s=0.1,
        traction_factory=_uniform_grid_factory((2, 2), 0.5),
        stiffness_rate_factory=_uniform_grid_factory((2, 2), 0.1),
        ligand_density_rate_factory=_uniform_grid_factory((2, 2), 0.05),
        fiber_density_rate_factory=_uniform_grid_factory((2, 2), 0.02),
        orientation_rate_factory=_uniform_diagonal_orientation_rate_factory(
            (2, 2), -0.05
        ),
    )
    run_ecm_ol_scenario(scenario, str(tmp_path / "immut_multi"), frame_interval=1)
    for name, snap in snapshots.items():
        np.testing.assert_array_equal(getattr(ecm, name), snap)


def test_frame_dump_files_valid_hdf5_after_multi_channel_step(tmp_path):
    import h5py

    ecm = make_default_ecm(nx=2, ny=2)
    scenario = EcmOlScenario(
        name="hdf5_multi",
        initial_ecm=ecm,
        n_steps=2,
        dt_s=0.1,
        stiffness_rate_factory=_uniform_grid_factory((2, 2), 0.5),
        ligand_density_rate_factory=_uniform_grid_factory((2, 2), 0.2),
    )
    run = run_ecm_ol_scenario(scenario, str(tmp_path / "hdf5_multi"), frame_interval=1)
    assert run.failure is None
    assert len(run.frame_paths) >= 2
    for frame_path in run.frame_paths:
        with h5py.File(frame_path, "r") as f:
            assert "schema_version" in f.attrs
            assert "ecm" in f
