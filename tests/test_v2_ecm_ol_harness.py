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
