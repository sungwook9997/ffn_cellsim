"""Tests for Phase E v1 Item 5 sweep harness.

Test catalog mirrors
``docs/v2_item_5_sweep_harness_locked.md`` §4 (commit ``c5634d5``)
and the impl-work Sanity Gate
``docs/v2_item_5_sweep_harness_sanity_gate.md`` §8 (commits
``338be1b`` + ``19825e9`` + ``774d963``). 24 tests total: 23
per locked §4 + 1 wording-boundary meta-test (test 24) added
per Codex ``id=1581`` BLOCKER.

Per Codex ``id=1587`` Sanity Gate PASS, this catalog covers:
- Validation tests (11): tests 1-11
- Determinism + metadata (3): tests 12-14
- Evidence-providing (4) Phase E v1 Y2 sister-pattern: tests 15-18
- Boundary + Step 0 (3): tests 19-21
- Composition guard + exports + wording boundary (3): tests 22-24
"""

from __future__ import annotations

import ast
import inspect
import math
import subprocess

import numpy as np
import pytest

from acs.v2 import dynamics as v2_dynamics
import acs.v2 as v2_pkg
from acs.v2.dynamics import closed_loop_phase_e_sweep as ce_sweep
from acs.v2.dynamics.closed_loop_phase_e_sweep import (
    PhaseEV1Item5SweepConfig,
    PhaseEV1Item5SweepMetadata,
    PhaseEV1Item5SweepResult,
    run_phase_e_v1_sensitivity_sweep,
)


def _baseline_config(**overrides) -> PhaseEV1Item5SweepConfig:
    """Helper: build a minimal valid config; overrides for individual tests."""

    defaults = dict(
        spacing_um_values=(1.0, 2.0),
        dt_s_values=(60.0, 120.0),
        n_steps=4,
        domain_size_um_xy=(8.0, 8.0),
        origin_um_xy=(0.0, 0.0),
        traction_scenario="rotating_uniform_single_fa",
        traction_magnitude_nN=2.0,
        fa_position_um_xy=(4.0, 4.0),
        n_revolutions=1,
        git_sha_override="test-sha-fixed",
    )
    defaults.update(overrides)
    return PhaseEV1Item5SweepConfig(**defaults)


def test_item_5_sweep_empty_spacing_values_raises():
    """Test 1 (Y3): empty spacing_um_values must raise."""

    with pytest.raises(ValueError, match="spacing_um_values"):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(spacing_um_values=()))


def test_item_5_sweep_empty_dt_values_raises():
    """Test 2 (Y3): empty dt_s_values must raise."""

    with pytest.raises(ValueError, match="dt_s_values"):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(dt_s_values=()))


def test_item_5_sweep_negative_or_nonfinite_value_raises():
    """Test 3 (Y3): non-positive / nonfinite values raise."""

    with pytest.raises(ValueError):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(spacing_um_values=(-1.0,)))
    with pytest.raises(ValueError):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(dt_s_values=(0.0,)))
    with pytest.raises(ValueError):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(spacing_um_values=(float("nan"),)))


def test_item_5_sweep_indivisible_domain_raises():
    """Test 4 (Y1): indivisible domain raises with explicit error."""

    with pytest.raises(ValueError, match="exactly divisible"):
        run_phase_e_v1_sensitivity_sweep(
            _baseline_config(spacing_um_values=(1.5,))
        )


def test_item_5_sweep_unsupported_traction_scenario_raises():
    """Test 5 (Y20): runtime Literal validation."""

    config = PhaseEV1Item5SweepConfig(
        spacing_um_values=(1.0,),
        dt_s_values=(60.0,),
        n_steps=4,
        domain_size_um_xy=(8.0, 8.0),
        origin_um_xy=(0.0, 0.0),
        traction_scenario="invalid_scenario",
        traction_magnitude_nN=2.0,
        fa_position_um_xy=(4.0, 4.0),
    )
    with pytest.raises(ValueError, match="traction_scenario"):
        run_phase_e_v1_sensitivity_sweep(config)


def test_item_5_sweep_zero_n_revolutions_raises():
    """Test 6 (Y20): n_revolutions == 0 raises."""

    with pytest.raises(ValueError, match="n_revolutions"):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(n_revolutions=0))


def test_item_5_sweep_negative_n_revolutions_raises():
    """Test 7 (Y20): n_revolutions < 0 raises."""

    with pytest.raises(ValueError, match="n_revolutions"):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(n_revolutions=-1))


def test_item_5_sweep_bool_n_revolutions_raises():
    """Test 8 (Y20): n_revolutions of type bool raises (Python bool ⊂ int trap)."""

    with pytest.raises(ValueError, match="n_revolutions"):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(n_revolutions=True))


def test_item_5_sweep_invalid_domain_size_raises():
    """Test 9 (Y21): domain_size invalid (length, finite, positive) raises."""

    with pytest.raises(ValueError, match="domain_size_um_xy"):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(domain_size_um_xy=(8.0,)))
    with pytest.raises(ValueError, match="domain_size_um_xy"):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(domain_size_um_xy=(0.0, 8.0)))


def test_item_5_sweep_invalid_origin_raises():
    """Test 10 (Y21): origin invalid (length, finite) raises."""

    with pytest.raises(ValueError, match="origin_um_xy"):
        run_phase_e_v1_sensitivity_sweep(_baseline_config(origin_um_xy=(0.0,)))
    with pytest.raises(ValueError, match="origin_um_xy"):
        run_phase_e_v1_sensitivity_sweep(
            _baseline_config(origin_um_xy=(float("inf"), 0.0))
        )


def test_item_5_sweep_fa_position_outside_domain_raises():
    """Test 11 (Y21): FA position outside inclusive physical domain raises."""

    with pytest.raises(ValueError, match="fa_position_um_xy"):
        run_phase_e_v1_sensitivity_sweep(
            _baseline_config(fa_position_um_xy=(20.0, 4.0))
        )


def test_item_5_sweep_metadata_uses_git_sha_override_when_provided():
    """Test 12 (Y14): test override path."""

    result = run_phase_e_v1_sensitivity_sweep(
        _baseline_config(git_sha_override="custom-test-sha")
    )
    assert result.metadata.git_sha == "custom-test-sha"


def test_item_5_sweep_metadata_falls_back_to_unknown_on_git_failure(monkeypatch):
    """Test 13 (Y14): subprocess failure → fallback to "unknown"."""

    def _raise(*args, **kwargs):
        raise FileNotFoundError("git not installed (test stub)")

    monkeypatch.setattr(ce_sweep.subprocess, "run", _raise)
    result = run_phase_e_v1_sensitivity_sweep(
        _baseline_config(git_sha_override=None)
    )
    assert result.metadata.git_sha == "unknown"


def test_item_5_sweep_no_timestamp_in_metadata():
    """Test 14 (Y13): no wall-clock timestamp field."""

    result = run_phase_e_v1_sensitivity_sweep(_baseline_config())
    assert not hasattr(result.metadata, "timestamp_iso")
    assert not hasattr(result.metadata, "timestamp")
    assert not hasattr(result.metadata, "wall_time_s")
    metadata_fields = set(PhaseEV1Item5SweepMetadata.__dataclass_fields__.keys())
    forbidden_fields = {"timestamp_iso", "timestamp", "wall_time_s"}
    assert not (metadata_fields & forbidden_fields), (
        f"PhaseEV1Item5SweepMetadata fields {metadata_fields} must not contain "
        f"{forbidden_fields} (Y13 reproducibility)"
    )


def test_item_5_sweep_provides_grid_spacing_sensitivity_evidence():
    """Test 15 (Y19 structural): all locked spacings present + nx/ny match."""

    config = _baseline_config(
        spacing_um_values=(1.0, 2.0, 4.0),
        dt_s_values=(60.0,),
    )
    result = run_phase_e_v1_sensitivity_sweep(config)

    assert set(result.spacing_um_by_run.tolist()) == {1.0, 2.0, 4.0}
    for i in range(len(result.spacing_um_by_run)):
        spacing = float(result.spacing_um_by_run[i])
        assert int(result.nx_by_run[i]) == int(round(8.0 / spacing))
        assert int(result.ny_by_run[i]) == int(round(8.0 / spacing))

    assert result.v_active_um2_traj.shape == (3, config.n_steps + 1)


def test_item_5_sweep_provides_dt_substep_sensitivity_evidence():
    """Test 16 (Y19 analytical): step-1 max_convex_weight monotone non-decreasing in dt."""

    config = _baseline_config(
        spacing_um_values=(1.0,),
        dt_s_values=(30.0, 60.0, 120.0, 300.0),
    )
    result = run_phase_e_v1_sensitivity_sweep(config)

    weights_step1 = result.max_convex_weight_traj[:, 1]
    for i in range(len(weights_step1) - 1):
        assert weights_step1[i] <= weights_step1[i + 1] + 1e-12, (
            f"Step-1 max_convex_weight must be monotone non-decreasing in dt; "
            f"weights_step1={weights_step1.tolist()}"
        )


def test_item_5_sweep_provides_v_active_bound_invariant_evidence():
    """Test 17 (Y18): max_bound_ratio ≤ 1.0 + 1e-12 everywhere.

    1e-12 is IEEE float64 roundoff allowance, NOT a tunable gate
    threshold. The mathematical invariant is ≤ 1.0 per HB#5.
    """

    result = run_phase_e_v1_sensitivity_sweep(_baseline_config())
    for ratio in result.max_bound_ratio:
        assert ratio <= 1.0 + 1e-12, (
            f"max_bound_ratio {ratio} exceeds 1.0 + 1e-12 IEEE roundoff "
            f"allowance; HB#5 mathematical invariant is ≤ 1.0"
        )


def test_item_5_sweep_provides_neutral_multiplier_preservation_evidence():
    """Test 18: HB#4 neutral hard-wired invariant — multipliers all 1.0 everywhere."""

    result = run_phase_e_v1_sensitivity_sweep(_baseline_config())
    np.testing.assert_array_equal(
        result.neutral_multiplier_max_deviation_traj,
        np.zeros_like(result.neutral_multiplier_max_deviation_traj),
    )
    np.testing.assert_array_equal(
        result.max_neutral_multiplier_deviation,
        np.zeros_like(result.max_neutral_multiplier_deviation),
    )


def test_item_5_sweep_max_bound_ratio_zero_when_no_active_cells():
    """Test 19 (Y16): degenerate near-zero traction scenario.

    With traction_magnitude_nN very small but the scenario still
    rotating, scatter produces nonzero traction → V_bound > 0 → ratio
    well-defined. The Y16 zero-active-bound case is internal to
    `_compute_max_bound_ratio` which we test directly via the
    private helper for the synthetic empty case.
    """

    empty_v = np.zeros(5)
    empty_b = np.zeros(5)
    assert ce_sweep._compute_max_bound_ratio(empty_v, empty_b) == 0.0

    v = np.array([0.0, 0.5, 0.0, 1.0, 0.0])
    b = np.array([0.0, 1.0, 0.0, 2.0, 0.0])
    assert ce_sweep._compute_max_bound_ratio(v, b) == 0.5


def test_item_5_sweep_trajectory_shape_n_runs_n_steps_plus_1():
    """Test 20 (Y6): explicit shape verification."""

    config = _baseline_config(
        spacing_um_values=(1.0, 2.0),
        dt_s_values=(60.0, 120.0, 300.0),
        n_steps=10,
    )
    result = run_phase_e_v1_sensitivity_sweep(config)

    n_runs_expected = 2 * 3
    n_steps_plus_1_expected = 10 + 1

    for traj in [
        result.v_active_um2_traj,
        result.v_active_bound_um2_traj,
        result.max_convex_weight_traj,
        result.max_orientation_delta_frobenius_traj,
        result.max_traction_norm_nN_per_um2_traj,
        result.neutral_multiplier_max_deviation_traj,
    ]:
        assert traj.shape == (n_runs_expected, n_steps_plus_1_expected)


def test_item_5_sweep_step_0_initial_against_target_not_zero_baseline():
    """Test 21 (Y17): step 0 = initial-against-target (NOT zero baseline)."""

    config = _baseline_config()
    result = run_phase_e_v1_sensitivity_sweep(config)

    assert (result.v_active_um2_traj[:, 0] > 0.0).all(), (
        "v_active_um2_traj[:, 0] must be > 0 for misaligned 0.5*I IC; "
        "Y17 step 0 is initial-against-target NOT zero baseline"
    )

    np.testing.assert_array_equal(
        result.max_convex_weight_traj[:, 0],
        np.zeros_like(result.max_convex_weight_traj[:, 0]),
    )

    np.testing.assert_array_equal(
        result.max_orientation_delta_frobenius_traj[:, 0],
        np.zeros_like(result.max_orientation_delta_frobenius_traj[:, 0]),
    )


def test_item_5_sweep_no_local_failure_kinds_or_phase_e_error_redefinition():
    """Test 22 (sister-pattern with Phase E v1 test 3): failure-kind discipline.

    Source/meta guard for failure-kind discipline ONLY; wording
    discipline owned by test 24 per Codex id=1581 separation.
    """

    src = inspect.getsource(ce_sweep)
    tree = ast.parse(src)

    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            assert not (
                node.name.startswith("Item5") and node.name.endswith("Error")
            ), (
                f"Item 5 sweep harness must not define wrapper failure "
                f"kinds; found class {node.name}"
            )

    src_no_spaces = src.replace(" ", "")
    assert "failure_kind=" not in src_no_spaces, (
        "Item 5 sweep harness must not introduce failure_kind labels; "
        "sub-call errors propagate as-is"
    )


def test_item_5_sweep_exports_through_both_init():
    """Test 23 (Y12): 4 public exports through both __init__ surfaces."""

    expected_symbols = {
        "PhaseEV1Item5SweepConfig",
        "PhaseEV1Item5SweepMetadata",
        "PhaseEV1Item5SweepResult",
        "run_phase_e_v1_sensitivity_sweep",
    }

    forbidden_helper_names = {
        "_validate_config",
        "_build_rotating_uniform_single_fa",
        "_build_initial_ecm_at_spacing",
        "_detect_git_sha",
        "_compute_max_bound_ratio",
        "_record_step",
        "_default_phase_e_v1_item5_sweep_config",
    }

    for sym in expected_symbols:
        assert hasattr(v2_pkg, sym), f"acs.v2 missing {sym!r}"
        assert hasattr(v2_dynamics, sym), f"acs.v2.dynamics missing {sym!r}"
        assert getattr(v2_pkg, sym) is getattr(ce_sweep, sym), (
            f"acs.v2.{sym} is not "
            f"acs.v2.dynamics.closed_loop_phase_e_sweep.{sym}"
        )
        assert getattr(v2_dynamics, sym) is getattr(ce_sweep, sym), (
            f"acs.v2.dynamics.{sym} is not "
            f"acs.v2.dynamics.closed_loop_phase_e_sweep.{sym}"
        )
        assert sym in v2_pkg.__all__, f"acs.v2.__all__ missing {sym!r}"
        assert sym in v2_dynamics.__all__, (
            f"acs.v2.dynamics.__all__ missing {sym!r}"
        )

    for sym in forbidden_helper_names:
        assert sym not in v2_pkg.__all__, (
            f"acs.v2.__all__ must NOT export private helper {sym!r}"
        )
        assert sym not in v2_dynamics.__all__, (
            f"acs.v2.dynamics.__all__ must NOT export private helper {sym!r}"
        )


def test_item_5_sweep_module_doc_does_not_overclaim_item_5_satisfaction():
    """Test 24 (NEW per Codex id=1581 BLOCKER): wording-boundary discipline.

    Sister-pattern with Phase E v1 test 14 wording-boundary forward
    guard. Asserts module docstring + function docstring + 3
    dataclass docstrings do NOT contain "Item 5 satisfied" /
    "satisfies_item_5" / "satisfies_*_item_5" forbidden substrings.
    Test names in this test file also do NOT match
    test_*_satisfies_item_5_* pattern.
    """

    forbidden_phrases = [
        "Item 5 satisfied",
        "satisfies_item_5",
    ]

    module_doc = ce_sweep.__doc__ or ""
    function_doc = run_phase_e_v1_sensitivity_sweep.__doc__ or ""
    config_doc = PhaseEV1Item5SweepConfig.__doc__ or ""
    metadata_doc = PhaseEV1Item5SweepMetadata.__doc__ or ""
    result_doc = PhaseEV1Item5SweepResult.__doc__ or ""

    for phrase in forbidden_phrases:
        assert phrase not in module_doc, (
            f"closed_loop_phase_e_sweep module docstring must not "
            f"contain forbidden phrase {phrase!r}"
        )
        assert phrase not in function_doc, (
            f"run_phase_e_v1_sensitivity_sweep docstring must not "
            f"contain forbidden phrase {phrase!r}"
        )
        assert phrase not in config_doc, (
            f"PhaseEV1Item5SweepConfig docstring must not contain "
            f"forbidden phrase {phrase!r}"
        )
        assert phrase not in metadata_doc, (
            f"PhaseEV1Item5SweepMetadata docstring must not contain "
            f"forbidden phrase {phrase!r}"
        )
        assert phrase not in result_doc, (
            f"PhaseEV1Item5SweepResult docstring must not contain "
            f"forbidden phrase {phrase!r}"
        )

    test_module = inspect.getmodule(test_item_5_sweep_module_doc_does_not_overclaim_item_5_satisfaction)
    test_function_names = [
        name
        for name, obj in inspect.getmembers(test_module, inspect.isfunction)
        if name.startswith("test_")
    ]
    for name in test_function_names:
        assert "satisfies_item_5" not in name, (
            f"Test name {name!r} must not contain forbidden substring "
            f"'satisfies_item_5'; use 'provides_*_evidence' wording per "
            f"Phase E v1 Y2 sister-pattern"
        )
