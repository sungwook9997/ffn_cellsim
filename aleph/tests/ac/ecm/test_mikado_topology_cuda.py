"""Native Warp-CUDA gate for Mikado chain and material-point candidate topology."""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.components.ecm import MikadoInitConfig, MikadoTopologyBuilder, PopulationSlot

_CUDA_DEVICE = next((device for device in wp.get_devices() if device.is_cuda), None)
pytestmark = pytest.mark.skipif(
    _CUDA_DEVICE is None,
    reason="I0-A: ECM topology native gate requires a CUDA GPU",
)


def _closest_segment_pair(
    p0: np.ndarray,
    p1: np.ndarray,
    q0: np.ndarray,
    q1: np.ndarray,
) -> tuple[float, float, float]:
    """Bit-for-formula host oracle for the device closest-segment function."""
    d1 = p1 - p0
    d2 = q1 - q0
    r = p0 - q0
    a = float(np.dot(d1, d1))
    e = float(np.dot(d2, d2))
    f = float(np.dot(d2, r))
    eps = 1.0e-24
    s = 0.0
    t = 0.0
    if a <= eps and e <= eps:
        return s, t, float(np.linalg.norm(p0 - q0))
    if a <= eps:
        t = float(np.clip(f / e, 0.0, 1.0))
    else:
        c = float(np.dot(d1, r))
        if e <= eps:
            s = float(np.clip(-c / a, 0.0, 1.0))
        else:
            b = float(np.dot(d1, d2))
            denominator = a * e - b * b
            if denominator > eps:
                s = float(np.clip((b * f - c * e) / denominator, 0.0, 1.0))
            t = (b * s + f) / e
            if t < 0.0:
                t = 0.0
                s = float(np.clip(-c / a, 0.0, 1.0))
            elif t > 1.0:
                t = 1.0
                s = float(np.clip((b - c) / a, 0.0, 1.0))
    point_p = p0 + s * d1
    point_q = q0 + t * d2
    return s, t, float(np.linalg.norm(point_p - point_q))


def test_gpu_mikado_initialization_and_material_point_candidates() -> None:
    cfg = MikadoInitConfig(
        box_lo_um=(-2.0, -2.0, -2.0),
        box_hi_um=(2.0, 2.0, 2.0),
        n_fibers=32,
        fiber_length_um=8.0,
        target_segment_um=1.0,
        crosslink_capture_um=1.0,
        pin_faces=("x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi"),
        pin_margin_um=0.25,
        rng_seed=1234,
        max_refinement_level=1,
        persistent_id_base=10_000,
    )
    builder = MikadoTopologyBuilder(cfg, device=str(_CUDA_DEVICE))
    topology = builder.initialize()
    candidates = builder.initial_contact_candidates(topology)
    wp.synchronize_device(str(_CUDA_DEVICE))

    positions = topology.position_d.numpy()
    reference = topology.reference_position_d.numpy()
    segments = topology.segments_d.numpy()
    segment_fiber = topology.segment_fiber_d.numpy()
    segment_active = topology.segment_active_d.numpy().astype(bool)
    bends = topology.bend_triples_d.numpy()
    bend_fiber = topology.bend_fiber_d.numpy()
    bend_active = topology.bend_active_d.numpy().astype(bool)
    segment_prev = topology.segment_prev_d.numpy()
    segment_next = topology.segment_next_d.numpy()
    segment_bend = topology.segment_bend_d.numpy()
    bend_left_segment = topology.bend_left_segment_d.numpy()
    bend_right_segment = topology.bend_right_segment_d.numpy()
    offsets = topology.fiber_offsets_d.numpy()
    rest = topology.segment_rest_length_d.numpy()
    populations = topology.active_population_d.numpy()
    free_populations = topology.free_population_d.numpy()
    node_active = topology.node_active_d.numpy().astype(bool)

    assert positions.shape == (cfg.n_fibers * cfg.node_slots_per_fiber, 3)
    np.testing.assert_array_equal(reference, positions)
    np.testing.assert_array_equal(
        offsets,
        np.arange(cfg.n_fibers + 1) * cfg.node_slots_per_fiber,
    )
    np.testing.assert_array_equal(segments[segment_active, 1], segments[segment_active, 0] + 1)
    np.testing.assert_array_equal(
        topology.node_fiber_d.numpy()[segments[segment_active, 0]],
        segment_fiber[segment_active],
    )
    np.testing.assert_array_equal(
        topology.node_fiber_d.numpy()[segments[segment_active, 1]],
        segment_fiber[segment_active],
    )
    assert np.all(rest[segment_active] > 0.0)
    assert np.all(rest[segment_active] <= cfg.max_segment_length_um * (1.0 + 1.0e-12))
    assert np.all(rest[~segment_active] == 0.0)
    if np.any(bend_active):
        np.testing.assert_array_equal(bends[bend_active, 1], bends[bend_active, 0] + 1)
        np.testing.assert_array_equal(bends[bend_active, 2], bends[bend_active, 0] + 2)
        np.testing.assert_array_equal(
            topology.node_fiber_d.numpy()[bends[bend_active, 1]],
            bend_fiber[bend_active],
        )
    for fiber in range(cfg.n_fibers):
        segment_base = fiber * cfg.segment_slots_per_fiber
        bend_base = fiber * cfg.bend_slots_per_fiber
        initial_segments_per_fiber = cfg.nodes_per_fiber - 1
        initial_bends_per_fiber = max(0, cfg.nodes_per_fiber - 2)
        for local in range(initial_segments_per_fiber):
            segment = segment_base + local
            expected_prev = segment - 1 if local > 0 else -1
            expected_next = segment + 1 if local + 1 < initial_segments_per_fiber else -1
            expected_bend = bend_base + local if local < initial_bends_per_fiber else -1
            assert segment_prev[segment] == expected_prev
            assert segment_next[segment] == expected_next
            assert segment_bend[segment] == expected_bend
        for local in range(initial_bends_per_fiber):
            bend = bend_base + local
            assert bend_left_segment[bend] == segment_base + local
            assert bend_right_segment[bend] == segment_base + local + 1
    assert np.count_nonzero(topology.boundary_face_mask_d.numpy()) > 0
    assert populations[PopulationSlot.FIBER] == topology.n_fiber_capacity
    assert populations[PopulationSlot.NODE] == topology.n_initial_nodes
    assert populations[PopulationSlot.SEGMENT] == topology.n_initial_segments
    assert populations[PopulationSlot.BEND] == topology.n_initial_bends
    assert free_populations[0] == topology.n_node_capacity - topology.n_initial_nodes
    assert free_populations[1] == topology.n_segment_capacity - topology.n_initial_segments
    assert free_populations[2] == topology.n_bend_capacity - topology.n_initial_bends
    node_free = topology.node_free_list_d.numpy()[: free_populations[0]]
    segment_free = topology.segment_free_list_d.numpy()[: free_populations[1]]
    bend_free = topology.bend_free_list_d.numpy()[: free_populations[2]]
    assert np.unique(node_free).size == node_free.size
    assert np.unique(segment_free).size == segment_free.size
    assert np.unique(bend_free).size == bend_free.size
    assert not np.any(node_active[node_free])
    assert not np.any(segment_active[segment_free])
    assert not np.any(bend_active[bend_free])
    assert topology.allocated_array_bytes() > 0
    persistent_fibers = topology.persistent_fiber_id_d.numpy()
    persistent_nodes = topology.persistent_node_id_d.numpy()[node_active]
    persistent_segments = topology.persistent_segment_id_d.numpy()[segment_active]
    all_initial_ids = np.concatenate([persistent_fibers, persistent_nodes, persistent_segments])
    assert np.unique(all_initial_ids).size == cfg.initial_identity_count
    assert all_initial_ids.min() == cfg.persistent_id_base
    assert all_initial_ids.max() == cfg.persistent_id_base + cfg.initial_identity_count - 1
    assert topology.next_persistent_id_d.numpy()[0] == (
        cfg.persistent_id_base + cfg.initial_identity_count
    )
    assert topology.topology_epoch_d.numpy()[0] == 0
    assert topology.topology_dirty_d.numpy()[0] == 0
    assert np.all(topology.persistent_node_id_d.numpy()[~node_active] == -1)
    assert np.all(topology.persistent_segment_id_d.numpy()[~segment_active] == -1)
    assert np.all(topology.refinement_parent_segment_id_d.numpy() == -1)
    np.testing.assert_array_equal(
        topology.refinement_root_segment_id_d.numpy()[segment_active],
        topology.persistent_segment_id_d.numpy()[segment_active],
    )
    assert np.all(topology.refinement_root_segment_id_d.numpy()[~segment_active] == -1)
    assert np.all(topology.refinement_path_d.numpy() == 0)

    assert candidates.initialization_count > 0
    segment_a = candidates.segment_a_d.numpy()
    segment_b = candidates.segment_b_d.numpy()
    fiber_a = candidates.fiber_a_d.numpy()
    fiber_b = candidates.fiber_b_d.numpy()
    u_a = candidates.u_a_d.numpy()
    u_b = candidates.u_b_d.numpy()
    separation = candidates.separation_d.numpy()
    assert np.all(segment_b > segment_a)
    assert np.all(fiber_b != fiber_a)
    assert np.all((u_a >= 0.0) & (u_a <= 1.0))
    assert np.all((u_b >= 0.0) & (u_b <= 1.0))
    assert np.all(separation <= cfg.crosslink_capture_um * (1.0 + 1.0e-12))
    np.testing.assert_array_equal(
        candidates.generation_a_d.numpy(),
        topology.segment_generation_d.numpy()[segment_a],
    )
    np.testing.assert_array_equal(
        candidates.generation_b_d.numpy(),
        topology.segment_generation_d.numpy()[segment_b],
    )
    np.testing.assert_array_equal(
        candidates.persistent_segment_a_d.numpy(),
        topology.persistent_segment_id_d.numpy()[segment_a],
    )
    np.testing.assert_array_equal(
        candidates.persistent_segment_b_d.numpy(),
        topology.persistent_segment_id_d.numpy()[segment_b],
    )
    assert candidates.allocated_array_bytes() > 0

    candidate_by_pair = {
        (int(a), int(b)): (float(ua), float(ub), float(distance))
        for a, b, ua, ub, distance in zip(
            segment_a,
            segment_b,
            u_a,
            u_b,
            separation,
            strict=True,
        )
    }
    assert len(candidate_by_pair) == candidates.initialization_count
    brute_by_pair: dict[tuple[int, int], tuple[float, float, float]] = {}
    live_segments = np.flatnonzero(segment_active)
    for ordinal_a, index_a in enumerate(live_segments):
        for index_b in live_segments[ordinal_a + 1 :]:
            if segment_fiber[index_a] == segment_fiber[index_b]:
                continue
            ua, ub, distance = _closest_segment_pair(
                positions[segments[index_a, 0]],
                positions[segments[index_a, 1]],
                positions[segments[index_b, 0]],
                positions[segments[index_b, 1]],
            )
            if distance <= cfg.crosslink_capture_um:
                brute_by_pair[(int(index_a), int(index_b))] = (ua, ub, distance)
    assert candidate_by_pair.keys() == brute_by_pair.keys()
    for pair, expected in brute_by_pair.items():
        np.testing.assert_allclose(candidate_by_pair[pair], expected, rtol=1.0e-12, atol=1.0e-12)


def test_gpu_mikado_seed_is_topologically_deterministic() -> None:
    cfg = MikadoInitConfig(
        box_lo_um=(-3.0, -2.0, -1.0),
        box_hi_um=(3.0, 2.0, 1.0),
        n_fibers=8,
        fiber_length_um=4.0,
        target_segment_um=0.5,
        crosslink_capture_um=0.2,
        pin_faces=("x_lo",),
        pin_margin_um=0.5,
        rng_seed=99,
        max_refinement_level=1,
        persistent_id_base=700,
    )
    first = MikadoTopologyBuilder(cfg, device=str(_CUDA_DEVICE)).initialize()
    second = MikadoTopologyBuilder(cfg, device=str(_CUDA_DEVICE)).initialize()
    wp.synchronize_device(str(_CUDA_DEVICE))
    for first_array, second_array in (
        (first.position_d, second.position_d),
        (first.segments_d, second.segments_d),
        (first.bend_triples_d, second.bend_triples_d),
        (first.persistent_fiber_id_d, second.persistent_fiber_id_d),
        (first.persistent_node_id_d, second.persistent_node_id_d),
        (first.persistent_segment_id_d, second.persistent_segment_id_d),
        (first.refinement_root_segment_id_d, second.refinement_root_segment_id_d),
        (first.refinement_path_d, second.refinement_path_d),
        (first.segment_prev_d, second.segment_prev_d),
        (first.segment_next_d, second.segment_next_d),
        (first.segment_bend_d, second.segment_bend_d),
        (first.bend_left_segment_d, second.bend_left_segment_d),
        (first.bend_right_segment_d, second.bend_right_segment_d),
        (first.boundary_face_mask_d, second.boundary_face_mask_d),
        (first.node_free_list_d, second.node_free_list_d),
        (first.segment_free_list_d, second.segment_free_list_d),
        (first.bend_free_list_d, second.bend_free_list_d),
    ):
        np.testing.assert_array_equal(first_array.numpy(), second_array.numpy())


def test_two_node_fiber_has_valid_zero_bend_schema() -> None:
    cfg = MikadoInitConfig(
        box_lo_um=(-1.0, -1.0, -1.0),
        box_hi_um=(1.0, 1.0, 1.0),
        n_fibers=2,
        fiber_length_um=0.25,
        target_segment_um=0.5,
        crosslink_capture_um=0.1,
        pin_faces=("z_lo",),
        pin_margin_um=0.1,
        rng_seed=4,
        max_refinement_level=0,
    )
    topology = MikadoTopologyBuilder(cfg, device=str(_CUDA_DEVICE)).initialize()
    wp.synchronize_device(str(_CUDA_DEVICE))
    assert cfg.nodes_per_fiber == 2
    assert topology.n_bend_capacity == 0
    assert topology.bend_triples_d.shape == (0, 3)
    assert topology.active_population_d.numpy()[PopulationSlot.BEND] == 0
    assert topology.free_population_d.numpy()[2] == 0
