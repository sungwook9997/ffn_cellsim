"""Focused B-tier tests for the Phase E v2 Item 5 sweep harness variant."""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import acs.v2 as v2_pkg
from acs.v2 import dynamics as v2_dynamics
from acs.v2.dynamics import closed_loop_phase_e_sweep as ce_sweep
from acs.v2.dynamics.closed_loop_phase_e_sweep import (
    PhaseEV2Item5SweepConfig,
    PhaseEV2Item5SweepMetadata,
    PhaseEV2Item5SweepResult,
    run_phase_e_v2_sensitivity_sweep,
)


def _baseline_v2_config(**overrides) -> PhaseEV2Item5SweepConfig:
    """Build a small valid v2 sweep config for focused tests."""

    defaults = dict(
        spacing_um_values=(1.0,),
        dt_s_values=(60.0,),
        n_steps=2,
        domain_size_um_xy=(4.0, 4.0),
        origin_um_xy=(0.0, 0.0),
        traction_scenario="rotating_uniform_single_fa",
        traction_magnitude_nN=2.0,
        fa_position_um_xy=(2.0, 2.0),
        k_active=1.0,
        n_revolutions=1,
        git_sha_override="test-v2-sweep-sha",
    )
    defaults.update(overrides)
    return PhaseEV2Item5SweepConfig(**defaults)


def test_phase_e_v2_sweep_runs_and_records_active_multiplier_trajectories():
    """V2 sweep returns typed arrays with active multiplier summary fields."""

    config = _baseline_v2_config(
        spacing_um_values=(1.0, 2.0),
        dt_s_values=(30.0, 60.0),
    )
    result = run_phase_e_v2_sensitivity_sweep(config)

    assert isinstance(result, PhaseEV2Item5SweepResult)
    assert isinstance(result.metadata, PhaseEV2Item5SweepMetadata)
    assert result.metadata.config is config
    assert result.metadata.git_sha == "test-v2-sweep-sha"

    n_runs = 4
    n_steps_plus_1 = config.n_steps + 1
    for traj in [
        result.v_active_um2_traj,
        result.v_active_bound_um2_traj,
        result.max_convex_weight_traj,
        result.max_orientation_delta_frobenius_traj,
        result.max_traction_norm_nN_per_um2_traj,
        result.active_multiplier_min_traj,
        result.active_multiplier_mean_traj,
        result.active_multiplier_max_traj,
        result.active_multiplier_max_deviation_traj,
    ]:
        assert traj.shape == (n_runs, n_steps_plus_1)


def test_phase_e_v2_sweep_provides_nonneutral_active_multiplier_evidence():
    """Active HB#4 readout becomes non-neutral after the first ECM update."""

    result = run_phase_e_v2_sensitivity_sweep(_baseline_v2_config())

    # Step 0 reads the isotropic 0.5*I IC after dt=0, so its deviatoric
    # component is zero and the active multiplier is neutral.
    np.testing.assert_allclose(result.active_multiplier_max_deviation_traj[:, 0], 0.0)

    # Step 1 uses the locked v2 active readout on the updated ECM. This is
    # evidence that the harness consumes step_closed_loop_phase_e_v2, not
    # the v1 neutral wrapper.
    assert (result.active_multiplier_max_deviation_traj[:, 1:] > 0.0).any()
    assert (result.active_multiplier_max_traj[:, 1:] > 1.0).any()


def test_phase_e_v2_sweep_invalid_k_active_propagates_upstream_error():
    """Invalid k_active raises the HB#4-active error; no local failure kind."""

    with pytest.raises(Exception) as excinfo:
        run_phase_e_v2_sensitivity_sweep(_baseline_v2_config(k_active=-1.0))

    assert getattr(excinfo.value, "failure_kind", None) == "k_active_invalid"


def test_phase_e_v2_sweep_inherits_v1_decimal_divisibility_validation():
    """The v2 sister harness reuses the exact-decimal spacing validation."""

    with pytest.raises(ValueError, match="exactly divisible"):
        run_phase_e_v2_sensitivity_sweep(
            _baseline_v2_config(spacing_um_values=(1.0 + 1e-11,))
        )


def test_phase_e_v2_sweep_exports_through_both_init_surfaces():
    """Public v2 sweep symbols are exported; private helpers are not."""

    expected_symbols = {
        "PhaseEV2Item5SweepConfig",
        "PhaseEV2Item5SweepMetadata",
        "PhaseEV2Item5SweepResult",
        "run_phase_e_v2_sensitivity_sweep",
    }
    forbidden_helpers = {
        "_validate_v2_config",
        "_as_phase_e_v1_config",
        "_record_step_v2",
    }

    for sym in expected_symbols:
        assert hasattr(v2_pkg, sym)
        assert hasattr(v2_dynamics, sym)
        assert getattr(v2_pkg, sym) is getattr(ce_sweep, sym)
        assert getattr(v2_dynamics, sym) is getattr(ce_sweep, sym)
        assert sym in v2_pkg.__all__
        assert sym in v2_dynamics.__all__

    for sym in forbidden_helpers:
        assert sym not in v2_pkg.__all__
        assert sym not in v2_dynamics.__all__


def test_phase_e_v2_sweep_wording_boundary_does_not_overclaim_item_5():
    """Docstrings and test names stay in evidence-only wording."""

    forbidden_phrases = [
        "Item 5 satisfied",
        "satisfies_item_5",
    ]
    docs = [
        ce_sweep.__doc__ or "",
        run_phase_e_v2_sensitivity_sweep.__doc__ or "",
        PhaseEV2Item5SweepConfig.__doc__ or "",
        PhaseEV2Item5SweepMetadata.__doc__ or "",
        PhaseEV2Item5SweepResult.__doc__ or "",
    ]
    for phrase in forbidden_phrases:
        for doc in docs:
            assert phrase not in doc

    test_module = inspect.getmodule(
        test_phase_e_v2_sweep_wording_boundary_does_not_overclaim_item_5
    )
    test_function_names = [
        name
        for name, obj in inspect.getmembers(test_module, inspect.isfunction)
        if name.startswith("test_")
    ]
    for name in test_function_names:
        assert "satisfies_item_5" not in name
