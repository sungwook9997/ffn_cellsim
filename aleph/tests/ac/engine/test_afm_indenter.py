"""Structural and host-oracle tests for the CUDA-only spherical AFM apparatus."""

import ast
import inspect
import textwrap
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import warp as wp

from aleph.engine.afm_indenter import (
    AFMIndenterForceHook,
    AFMIndenterRuntime,
    AFMIndenterRuntimeSpec,
    sphere_contact_reference,
)
from aleph.engine.connector_joints import _unilateral_contact_force_kernel

_CUDA_DEVICE = next((device for device in wp.get_devices() if device.is_cuda), None)


def _spec() -> AFMIndenterRuntimeSpec:
    return AFMIndenterRuntimeSpec(
        radius_um=2.0,
        surface_contact_gap_um=0.1,
        stiffness_pn_per_um=10.0,
        provenance="test apparatus calibration",
    )


def test_centre_point_representation_is_force_free_at_and_beyond_the_sphere() -> None:
    spec = _spec()
    assert spec.centre_contact_distance_um == pytest.approx(2.1)
    assert sphere_contact_reference(2.1, spec) == 0.0
    assert sphere_contact_reference(2.2, spec) == 0.0


def test_overlap_is_repulsive_and_uses_the_existing_contact_sign() -> None:
    spec = _spec()
    assert sphere_contact_reference(2.0, spec) == pytest.approx(-1.0)


def test_runtime_imports_the_exact_sixth_edge_kernel_and_reaction_channel_zero() -> None:
    source = (Path(__file__).resolve().parents[3] / "engine" / "afm_indenter.py").read_text()
    assert _unilateral_contact_force_kernel.__name__ in source
    assert "reaction[0] = centre_force[0][2]" in source
    assert "_conditional_restore_centre_kernel" in source
    assert "resolved.is_cuda" in source


def test_force_hook_rejects_a_membrane_block_outside_the_global_layout() -> None:
    runtime = SimpleNamespace(n_pairs=5, spec=SimpleNamespace(stiffness_pn_per_um=10.0))
    with pytest.raises(ValueError, match="must lie inside"):
        AFMIndenterForceHook(runtime=runtime, membrane_node_offset=8, total_node_count=10)


def test_every_inner_force_reassembly_receives_the_experiment_hook() -> None:
    from aleph.components.incumbent.driver import make_inner_solve

    tree = ast.parse(textwrap.dedent(inspect.getsource(make_inner_solve)))
    calls = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        and node.func.id == "_accumulate_all"
    ]
    assert len(calls) >= 4
    for call in calls:
        keywords = {keyword.arg for keyword in call.keywords}
        assert "additional_force_accumulators" in keywords


@pytest.mark.parametrize(
    ("field", "value"),
    [("radius_um", 0.0), ("surface_contact_gap_um", -1.0), ("stiffness_pn_per_um", 0.0)],
)
def test_apparatus_magnitudes_are_never_defaulted_or_nonpositive(field: str, value: float) -> None:
    values = dict(radius_um=2.0, surface_contact_gap_um=0.1, stiffness_pn_per_um=10.0)
    values[field] = value
    with pytest.raises(ValueError):
        AFMIndenterRuntimeSpec(**values, provenance="test")


@pytest.mark.skipif(_CUDA_DEVICE is None, reason="Warp CUDA gate runs on the shared workstation")
def test_cuda_contact_is_repulsive_conservative_and_accepted_gated() -> None:
    """One exact contact exercises the reused kernel, reaction reduction and apparatus rollback."""
    device = str(_CUDA_DEVICE)
    position_d = wp.array(
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.9], [0.0, 0.0, 1.0]], dtype=np.float64),
        dtype=wp.vec3d,
        device=device,
    )
    force_d = wp.zeros(3, dtype=wp.vec3d, device=device)
    runtime = AFMIndenterRuntime(
        membrane_position_d=position_d,
        membrane_force_d=force_d,
        initial_centre_um=(0.0, 0.0, 2.0),
        spec=AFMIndenterRuntimeSpec(
            radius_um=1.0,
            surface_contact_gap_um=0.1,
            stiffness_pn_per_um=10.0,
            provenance="CUDA sign/conservation oracle",
        ),
        device=device,
    )
    runtime.accumulate()
    membrane_force = force_d.numpy()
    centre_force = runtime.centre_force_d.numpy()
    reaction = runtime.reaction_d.numpy()
    assert membrane_force[:, 2] == pytest.approx([0.0, 0.0, -1.0])
    assert centre_force[0, 2] == pytest.approx(1.0)
    assert (membrane_force.sum(axis=0) + centre_force[0]) == pytest.approx(np.zeros(3), abs=1e-14)
    assert reaction == pytest.approx([1.0, 1.0])

    global_position_d = wp.array(
        np.array([[9.0, 9.0, 9.0], [8.0, 8.0, 8.0], [0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.9], [0.0, 0.0, 1.0]], dtype=np.float64),
        dtype=wp.vec3d,
        device=device,
    )
    global_force_d = wp.zeros(5, dtype=wp.vec3d, device=device)
    hook = AFMIndenterForceHook(runtime=runtime, membrane_node_offset=2, total_node_count=5)
    hook(global_position_d, global_force_d)
    hooked_force = global_force_d.numpy()
    assert hooked_force[:2] == pytest.approx(np.zeros((2, 3)))
    assert hooked_force[2:, 2] == pytest.approx([0.0, 0.0, -1.0])
    assert runtime.reaction_d.numpy() == pytest.approx([1.0, 1.0])

    runtime.snapshot_candidate()
    runtime.set_candidate_depth(0.2)
    rejected_d = wp.zeros(1, dtype=wp.int32, device=device)
    runtime.rollback(rejected_d)
    assert runtime.centre_d.numpy()[0] == pytest.approx([0.0, 0.0, 2.0])

    runtime.snapshot_candidate()
    runtime.set_candidate_depth(0.2)
    accepted_d = wp.ones(1, dtype=wp.int32, device=device)
    runtime.rollback(accepted_d)
    assert runtime.centre_d.numpy()[0] == pytest.approx([0.0, 0.0, 1.8])
