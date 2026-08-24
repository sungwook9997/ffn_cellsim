r"""Host gate for the segment-anchored crossbridge (P0#2 production motor).

Validates the NumPy oracle :func:`crossbridge_segment_reference` the Warp kernels mirror — the novel physics is
that a head bound to a POINT on a segment loads the two real actin nodes by the barycentric weights ``(1−t):t``
while Newton's 3rd law holds across the segment. CUDA kernels are not launched here (dev-Mac I0-A); on-device
parity is a native gate. The gated attach/step kernels are exercised on the A5000 via the scheduler.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.components.motor.hand import NMIIHandParams
from aleph.components.motor.segment_motor import (
    SegmentMotorRuntime,
    allocate_segment_hand_state,
    crossbridge_segment_reference,
    crossbridge_segment_split_reference,
    step_detach_segment_gated_kernel,
)


def _cuda_device() -> str | None:
    wp.init()
    device = next((item for item in wp.get_devices() if item.is_cuda), None)
    return str(device) if device is not None else None


def test_newton_third_law_across_segment() -> None:
    """f_head + f_seg_a + f_seg_b = 0 exactly (barycentric reaction split loads two real nodes)."""
    r = crossbridge_segment_reference(
        head_pos=np.array([0.0, 0.0, 0.0]), pos_a=np.array([0.0, 0.2, 0.0]), pos_b=np.array([1.0, 0.2, 0.0]),
        bary_t=0.3, abscissa=0.05, walk_dir=np.array([-1.0, 0.0, 0.0]), k_xb=1000.0)
    net = r["f_head"] + r["f_seg_a"] + r["f_seg_b"]
    assert np.abs(net).max() < 1e-9


def test_barycentric_reaction_split() -> None:
    """The reaction on the two segment nodes is −(1−t)·f_head and −t·f_head (t = the attachment fraction)."""
    t = 0.25
    r = crossbridge_segment_reference(
        head_pos=np.array([0.1, 0.0, 0.0]), pos_a=np.array([0.0, 0.15, 0.0]), pos_b=np.array([0.5, 0.15, 0.0]),
        bary_t=t, abscissa=0.02, walk_dir=np.array([1.0, 0.0, 0.0]), k_xb=800.0)
    assert r["f_seg_a"] == pytest.approx(-(1.0 - t) * r["f_head"])
    assert r["f_seg_b"] == pytest.approx(-t * r["f_head"])
    # the endpoint nearer the attachment (larger weight) carries the larger share
    assert np.linalg.norm(r["f_seg_a"]) > np.linalg.norm(r["f_seg_b"])   # (1−t)=0.75 > t=0.25


def test_attachment_advances_along_walk_dir() -> None:
    """x_att = (1−t)·a + t·b + abscissa·walk_dir — the power stroke moves the anchor along the barbed dir."""
    a, b, t, absc = np.array([0.0, 0, 0]), np.array([1.0, 0, 0]), 0.5, 0.1
    wd = np.array([0.0, 1.0, 0.0])
    r = crossbridge_segment_reference(np.zeros(3), a, b, t, absc, wd, k_xb=1000.0)
    assert r["x_att"] == pytest.approx([0.5, 0.1, 0.0])                  # midpoint + 0.1 along +y


def test_directional_crossbridge_load_is_walk_tangent_only() -> None:
    """Directional (defect#3, 2026-07-23) crossbridge: the load is the walk-tangent projection ONLY.

    The directional fix makes the crossbridge force act purely along ``walk_dir`` with magnitude equal to the
    tangential extension (``crossbridge_segment_reference`` — ``f = load_tangential · walk_dir``), so a
    transverse displacement generates NO perpendicular force and does NOT add to the load.  This is the
    invariant that removed the spurious 3D normal force (645→52 pN, defect#3); ``load_full`` is
    ``|load_tangential|`` by construction, collinear or not.  (Supersedes the pre-fix expectation
    ``load_full > load_tangential`` off-axis, which encoded the deleted full-displacement model.)
    """
    # collinear: displacement purely along walk_dir ⇒ tangential == full == k_xb · extension
    r = crossbridge_segment_reference(
        head_pos=np.array([0.0, 0.0, 0.0]), pos_a=np.array([0.005, 0.0, 0.0]), pos_b=np.array([1.005, 0.0, 0.0]),
        bary_t=0.0, abscissa=0.0, walk_dir=np.array([1.0, 0.0, 0.0]), k_xb=1000.0)
    assert r["load_full"] == pytest.approx(abs(r["load_tangential"]), rel=1e-9)
    assert r["load_tangential"] == pytest.approx(5.0, rel=1e-9)
    # off-axis: the transverse component is IGNORED — |F| still equals |tangential|, the load equals the same
    # walk-tangent projection as the collinear case, and the emitted force has zero transverse component.
    r2 = crossbridge_segment_reference(
        head_pos=np.array([0.0, 0.0, 0.0]), pos_a=np.array([0.005, 0.01, 0.0]), pos_b=np.array([1.0, 0.01, 0.0]),
        bary_t=0.0, abscissa=0.0, walk_dir=np.array([1.0, 0.0, 0.0]), k_xb=1000.0)
    assert r2["load_full"] == pytest.approx(abs(r2["load_tangential"]), rel=1e-9)
    assert r2["load_tangential"] == pytest.approx(r["load_tangential"], rel=1e-9)  # transverse offset ignored
    assert r2["f_head"][1] == pytest.approx(0.0, abs=1e-12)
    assert r2["f_head"][2] == pytest.approx(0.0, abs=1e-12)


def test_split_two_array_scatter_preserves_combined_physics() -> None:
    """The two-array (split-ownership) crossbridge is a pure relayout of the single-array formula.

    The head node lives in the ``nmii``-owned actuator array and the two segment nodes in the target-owned port
    array; the arrays are NOT co-located (different lengths + a head index that would collide with a segment
    index if merged). The per-array scatter must equal the single-array reference term-for-term, and Newton's
    3rd law must hold across the two arrays.
    """
    # actuator array: 3 head particles; port array: 4 actin nodes. Head index 1 == segment index 1 on purpose:
    # if the relayout secretly merged the arrays that collision would corrupt the scatter.
    actuator_pos = np.array([[9.0, 9.0, 9.0], [0.0, 0.0, 0.0], [9.0, 9.0, 9.0]], dtype=np.float64)
    port_pos = np.array(
        [[0.0, 0.2, 0.0], [1.0, 0.2, 0.0], [5.0, 5.0, 5.0], [6.0, 6.0, 6.0]], dtype=np.float64)
    head_node, seg_a, seg_b = 1, 0, 1
    t, absc, wd, k_xb = 0.3, 0.05, np.array([-1.0, 0.0, 0.0]), 1000.0

    split = crossbridge_segment_split_reference(
        actuator_pos, port_pos, head_node, seg_a, seg_b, t, absc, wd, k_xb)
    combined = crossbridge_segment_reference(
        actuator_pos[head_node], port_pos[seg_a], port_pos[seg_b], t, absc, wd, k_xb)

    # identical values: the relayout changes only WHERE forces land, not the physics.
    assert split["f_head"] == pytest.approx(combined["f_head"])
    assert split["f_seg_a"] == pytest.approx(combined["f_seg_a"])
    assert split["f_seg_b"] == pytest.approx(combined["f_seg_b"])
    assert split["load_tangential"] == pytest.approx(combined["load_tangential"])
    assert split["load_full"] == pytest.approx(combined["load_full"])

    # head force lands ONLY in the actuator array; segment reactions ONLY in the port array.
    assert split["actuator_force"][head_node] == pytest.approx(combined["f_head"])
    assert np.abs(np.delete(split["actuator_force"], head_node, axis=0)).max() == 0.0
    assert split["port_force"][seg_a] == pytest.approx(combined["f_seg_a"])
    assert split["port_force"][seg_b] == pytest.approx(combined["f_seg_b"])

    # Newton's 3rd law ACROSS the two independent arrays (no shared storage): total force sums to zero.
    net = split["actuator_force"].sum(axis=0) + split["port_force"].sum(axis=0)
    assert np.abs(net).max() < 1e-9


def test_cuda_outer_acceptance_predicates_segment_step() -> None:
    """Rejected outer steps preserve state exactly; accepted steps advance with the derived stall guard."""
    device = _cuda_device()
    if device is None:
        pytest.skip("I0-A: accepted-predicated segment KMC requires a CUDA GPU")
    state = allocate_segment_hand_state(1, device=device)
    state["bound"].fill_(1)
    state["seg_id"].fill_(0)
    state["seg_a"].fill_(0)
    state["seg_b"].fill_(1)
    params = NMIIHandParams()
    params.k_on = 0.0
    params.k_off0 = 0.0
    params.f0 = 1.0
    params.v0 = 0.2
    params.f_stall = 0.5
    params.kappa = 1.0
    params.k_xb = 1.0
    params.r0_head = 0.0
    params.r0_xb = 0.0
    params.capture_radius = 0.1
    loads = wp.zeros(1, dtype=wp.float64, device=device)
    rejected = wp.zeros(1, dtype=wp.int32, device=device)
    accepted = wp.ones(1, dtype=wp.int32, device=device)
    rng_epoch = wp.zeros(1, dtype=wp.int32, device=device)

    args = [
        state["bound"], state["seg_id"], state["seg_a"], state["seg_b"], state["bary_t"],
        state["abscissa"], state["walk_dir"], loads, loads, params, wp.float64(0.1), wp.int32(17),
    ]
    wp.launch(step_detach_segment_gated_kernel, dim=1, inputs=[rejected, rng_epoch, *args], device=device)
    assert float(state["abscissa"].numpy()[0]) == 0.0
    wp.launch(step_detach_segment_gated_kernel, dim=1, inputs=[accepted, rng_epoch, *args], device=device)
    assert float(state["abscissa"].numpy()[0]) == pytest.approx(0.02)


def test_cuda_segment_runtime_query_commit_and_epoch_are_transactional() -> None:
    """The composed runtime queries a segment on GPU and only an accepted commit changes state/RNG epoch."""
    device = _cuda_device()
    if device is None:
        pytest.skip("I0-A: segment runtime transaction gate requires a CUDA GPU")
    pos = wp.array(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.01, 0.0]], dtype=np.float64),
        dtype=wp.vec3d,
        device=device,
    )
    head_node = wp.array(np.array([2], dtype=np.int32), dtype=wp.int32, device=device)
    seg_a = wp.array(np.array([0], dtype=np.int32), dtype=wp.int32, device=device)
    seg_b = wp.array(np.array([1], dtype=np.int32), dtype=wp.int32, device=device)
    polarity = wp.array(np.array([1], dtype=np.int32), dtype=wp.int32, device=device)
    params = NMIIHandParams()
    params.k_on = 1.0e6
    params.k_off0 = 0.0
    params.f0 = 1.0
    params.v0 = 0.1
    params.f_stall = 1.0
    params.kappa = 1.0
    params.k_xb = 1.0
    params.r0_head = 0.0
    params.r0_xb = 0.0
    params.capture_radius = 0.05
    runtime = SegmentMotorRuntime(
        head_node,
        params,
        seg_a,
        seg_b,
        polarity,
        max_segment_length_um=1.0,
        grid_dim=8,
        device=device,
    )
    rejected = wp.zeros(1, dtype=wp.int32, device=device)
    accepted = wp.ones(1, dtype=wp.int32, device=device)

    runtime.commit_kinetics(pos, rejected, 1.0, 123)
    assert int(runtime.state["bound"].numpy()[0]) == 0
    assert int(runtime.rng_epoch.numpy()[0]) == 0

    runtime.commit_kinetics(pos, accepted, 1.0, 123)
    assert int(runtime.state["bound"].numpy()[0]) == 1
    assert int(runtime.state["seg_a"].numpy()[0]) == 0
    assert int(runtime.state["seg_b"].numpy()[0]) == 1
    assert float(runtime.state["bary_t"].numpy()[0]) == pytest.approx(0.5)
    np.testing.assert_allclose(runtime.state["walk_dir"].numpy()[0], [1.0, 0.0, 0.0])
    assert float(runtime.state["abscissa"].numpy()[0]) == pytest.approx(0.1)
    assert int(runtime.rng_epoch.numpy()[0]) == 1

    snapshot = {name: array.numpy().copy() for name, array in runtime.state.items()}
    runtime.commit_kinetics(pos, rejected, 1.0, 123)
    for name, expected in snapshot.items():
        np.testing.assert_array_equal(runtime.state[name].numpy(), expected)
    assert int(runtime.rng_epoch.numpy()[0]) == 1
