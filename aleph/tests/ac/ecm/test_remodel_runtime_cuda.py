"""Native CUDA gates for accepted-step ECM remodeling transactions."""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.components.ecm import (
    ECMRemodelRuntime,
    FiberRemodelEvent,
    MikadoInitConfig,
    MikadoTopologyBuilder,
    PopulationSlot,
    ProposalStatus,
    RemapKind,
    RemodelEvent,
)

_CUDA_DEVICE = next((device for device in wp.get_devices() if device.is_cuda), None)
pytestmark = pytest.mark.skipif(
    _CUDA_DEVICE is None,
    reason="I0-A: ECM remodeling gate requires a CUDA GPU",
)


def _config() -> MikadoInitConfig:
    return MikadoInitConfig(
        box_lo_um=(-10.0, -10.0, -10.0),
        box_hi_um=(10.0, 10.0, 10.0),
        n_fibers=2,
        fiber_length_um=4.0,
        target_segment_um=1.0,
        crosslink_capture_um=0.2,
        pin_faces=("x_lo",),
        pin_margin_um=0.01,
        rng_seed=741,
        max_refinement_level=2,
        persistent_id_base=50_000,
    )


def _put(array: wp.array, values: np.ndarray) -> None:
    source = wp.array(values, dtype=array.dtype, device=str(_CUDA_DEVICE))
    wp.copy(array, source)


def _set_i32(array: wp.array, index: int, value: int) -> None:
    values = array.numpy()
    values[index] = value
    _put(array, values.astype(np.int32, copy=False))


def _set_f64(array: wp.array, index: int, value: float) -> None:
    values = array.numpy()
    values[index] = value
    _put(array, values.astype(np.float64, copy=False))


def _accepted(value: int) -> wp.array:
    return wp.array([value], dtype=wp.int32, device=str(_CUDA_DEVICE))


def _ack_split(runtime: ECMRemodelRuntime, source: int, left: int, right: int) -> None:
    _set_i32(runtime.endpoint_remap.ack_count_d, source, left + right)
    _set_i32(runtime.endpoint_remap.target_a_refcount_d, source, left)
    _set_i32(runtime.endpoint_remap.target_b_refcount_d, source, right)


def _ack_merge(runtime: ECMRemodelRuntime, source: int, count: int) -> None:
    _set_i32(runtime.endpoint_remap.ack_count_d, source, count)
    _set_i32(runtime.endpoint_remap.target_a_refcount_d, source, count)


def _assert_chain_and_bends(topology, fiber: int) -> list[int]:  # type: ignore[no-untyped-def]
    active = topology.segment_active_d.numpy().astype(bool)
    owners = topology.segment_fiber_d.numpy()
    previous = topology.segment_prev_d.numpy()
    following = topology.segment_next_d.numpy()
    segment_bend = topology.segment_bend_d.numpy()
    segments = topology.segments_d.numpy()
    bend_active = topology.bend_active_d.numpy().astype(bool)
    bend_triples = topology.bend_triples_d.numpy()
    bend_left = topology.bend_left_segment_d.numpy()
    bend_right = topology.bend_right_segment_d.numpy()
    heads = np.flatnonzero(active & (owners == fiber) & (previous < 0))
    assert heads.size == 1
    chain: list[int] = []
    segment = int(heads[0])
    while segment >= 0:
        assert segment not in chain
        chain.append(segment)
        segment = int(following[segment])
    assert len(chain) == np.count_nonzero(active & (owners == fiber))
    for ordinal, segment in enumerate(chain):
        assert previous[segment] == (-1 if ordinal == 0 else chain[ordinal - 1])
        assert following[segment] == (-1 if ordinal + 1 == len(chain) else chain[ordinal + 1])
        if ordinal + 1 == len(chain):
            assert segment_bend[segment] == -1
        else:
            right = chain[ordinal + 1]
            bend = int(segment_bend[segment])
            assert bend_active[bend]
            assert bend_left[bend] == segment
            assert bend_right[bend] == right
            np.testing.assert_array_equal(
                bend_triples[bend],
                (segments[segment, 0], segments[segment, 1], segments[right, 1]),
            )
    return chain


def test_refine_reject_ack_guard_and_conservative_coarsen() -> None:
    cfg = _config()
    topology = MikadoTopologyBuilder(cfg, device=str(_CUDA_DEVICE)).initialize()
    runtime = ECMRemodelRuntime(topology, max_refinement_level=cfg.max_refinement_level)
    source = 1
    _set_i32(topology.endpoint_refcount_d, source, 3)
    original_positions = topology.position_d.numpy().copy()
    original_segments = topology.segments_d.numpy().copy()
    original_active = topology.active_population_d.numpy().copy()
    original_free = topology.free_population_d.numpy().copy()
    original_next_id = int(topology.next_persistent_id_d.numpy()[0])
    original_rest = float(topology.segment_rest_length_d.numpy()[source])
    original_s0 = float(topology.segment_material_s0_d.numpy()[source])
    original_s1 = float(topology.segment_material_s1_d.numpy()[source])
    original_nodes = topology.segments_d.numpy()[source].copy()

    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.segment_event_d, source, RemodelEvent.REFINE)
    _set_f64(runtime.proposal.damage_increment_d, source, 0.25)
    runtime.prepare_candidate()
    assert runtime.candidate_valid_d.numpy()[0] == 0
    assert runtime.plan.plan_valid_d.numpy()[0] == 1
    assert runtime.endpoint_remap.kind_d.numpy()[source] == RemapKind.SPLIT
    _ack_split(runtime, source, 1, 1)  # one endpoint deliberately unacknowledged
    runtime.validate_remap_acknowledgements()
    assert runtime.candidate_valid_d.numpy()[0] == 0
    runtime.commit_irreversible(_accepted(1), dt_phys=0.01, rng_seed=1)
    np.testing.assert_array_equal(topology.position_d.numpy(), original_positions)
    np.testing.assert_array_equal(topology.segments_d.numpy(), original_segments)
    np.testing.assert_array_equal(topology.active_population_d.numpy(), original_active)
    assert topology.topology_epoch_d.numpy()[0] == 0

    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.segment_event_d, source, RemodelEvent.REFINE)
    _set_f64(runtime.proposal.damage_increment_d, source, 0.25)
    runtime.prepare_candidate()
    right = int(runtime.endpoint_remap.target_b_segment_d.numpy()[source])
    midpoint = int(runtime.plan.new_node_slot_d.numpy()[source])
    assert right // cfg.segment_slots_per_fiber == 0
    assert midpoint // cfg.node_slots_per_fiber == 0
    _ack_split(runtime, source, 1, 2)
    runtime.validate_remap_acknowledgements()
    assert runtime.candidate_valid_d.numpy()[0] == 1

    runtime.rollback(_accepted(0))
    runtime.commit_irreversible(_accepted(0), dt_phys=0.01, rng_seed=2)
    np.testing.assert_array_equal(topology.position_d.numpy(), original_positions)
    np.testing.assert_array_equal(topology.segments_d.numpy(), original_segments)
    np.testing.assert_array_equal(topology.active_population_d.numpy(), original_active)
    np.testing.assert_array_equal(topology.free_population_d.numpy(), original_free)
    assert topology.next_persistent_id_d.numpy()[0] == original_next_id

    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.segment_event_d, source, RemodelEvent.REFINE)
    _set_f64(runtime.proposal.damage_increment_d, source, 0.25)
    runtime.prepare_candidate()
    right = int(runtime.endpoint_remap.target_b_segment_d.numpy()[source])
    midpoint = int(runtime.plan.new_node_slot_d.numpy()[source])
    _ack_split(runtime, source, 1, 2)
    runtime.validate_remap_acknowledgements()
    runtime.commit_irreversible(_accepted(1), dt_phys=0.01, rng_seed=3)
    wp.synchronize_device(str(_CUDA_DEVICE))

    positions = topology.position_d.numpy()
    np.testing.assert_allclose(
        positions[midpoint],
        0.5 * (original_positions[original_nodes[0]] + original_positions[original_nodes[1]]),
        rtol=0.0,
        atol=1.0e-12,
    )
    rest = topology.segment_rest_length_d.numpy()
    s0 = topology.segment_material_s0_d.numpy()
    s1 = topology.segment_material_s1_d.numpy()
    assert rest[source] + rest[right] == pytest.approx(original_rest)
    assert s0[source] == pytest.approx(original_s0)
    assert s1[right] == pytest.approx(original_s1)
    assert s1[source] == pytest.approx(s0[right])
    assert topology.endpoint_refcount_d.numpy()[source] == 1
    assert topology.endpoint_refcount_d.numpy()[right] == 2
    assert topology.segment_damage_d.numpy()[source] == pytest.approx(0.25)
    assert topology.segment_damage_d.numpy()[right] == pytest.approx(0.25)
    assert topology.active_population_d.numpy()[PopulationSlot.NODE] == original_active[1] + 1
    assert topology.active_population_d.numpy()[PopulationSlot.SEGMENT] == original_active[2] + 1
    assert topology.active_population_d.numpy()[PopulationSlot.BEND] == original_active[3] + 1
    np.testing.assert_array_equal(topology.free_population_d.numpy(), original_free - 1)
    assert topology.next_persistent_id_d.numpy()[0] == original_next_id + 3
    assert topology.topology_epoch_d.numpy()[0] == 1
    assert topology.topology_dirty_d.numpy()[0] == 1
    _assert_chain_and_bends(topology, 0)

    forces = topology.force_d.numpy()
    forces[midpoint] = (4.0, -3.0, 2.0)
    _put(topology.force_d, forces)
    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.segment_event_d, source, RemodelEvent.COARSEN)
    _set_f64(runtime.proposal.damage_increment_d, source, 0.10)
    _set_f64(runtime.proposal.damage_increment_d, right, 0.30)
    runtime.prepare_candidate()
    assert runtime.endpoint_remap.kind_d.numpy()[source] == RemapKind.MERGE_LEFT
    assert runtime.endpoint_remap.kind_d.numpy()[right] == RemapKind.MERGE_RIGHT
    _ack_merge(runtime, source, 1)
    _ack_merge(runtime, right, 2)
    runtime.validate_remap_acknowledgements()
    assert runtime.candidate_valid_d.numpy()[0] == 1
    runtime.commit_irreversible(_accepted(1), dt_phys=0.01, rng_seed=4)
    wp.synchronize_device(str(_CUDA_DEVICE))

    np.testing.assert_array_equal(topology.active_population_d.numpy(), original_active)
    np.testing.assert_array_equal(topology.free_population_d.numpy(), original_free)
    assert topology.segment_active_d.numpy()[right] == 0
    assert topology.node_active_d.numpy()[midpoint] == 0
    np.testing.assert_array_equal(topology.force_d.numpy()[midpoint], np.zeros(3))
    assert topology.endpoint_refcount_d.numpy()[source] == 3
    assert topology.segment_rest_length_d.numpy()[source] == pytest.approx(original_rest)
    assert topology.segment_material_s0_d.numpy()[source] == pytest.approx(original_s0)
    assert topology.segment_material_s1_d.numpy()[source] == pytest.approx(original_s1)
    assert topology.segment_damage_d.numpy()[source] == pytest.approx(0.45)
    assert topology.segment_refinement_level_d.numpy()[source] == 0
    assert topology.refinement_path_d.numpy()[source] == 0
    assert topology.next_persistent_id_d.numpy()[0] == original_next_id + 4
    assert topology.topology_epoch_d.numpy()[0] == 2
    _assert_chain_and_bends(topology, 0)

    for active, free_list, count in (
        (topology.node_active_d, topology.node_free_list_d, topology.free_population_d.numpy()[0]),
        (
            topology.segment_active_d,
            topology.segment_free_list_d,
            topology.free_population_d.numpy()[1],
        ),
        (topology.bend_active_d, topology.bend_free_list_d, topology.free_population_d.numpy()[2]),
    ):
        free = free_list.numpy()[:count]
        assert np.unique(free).size == count
        assert not np.any(active.numpy().astype(bool)[free])
    runtime.mark_candidate_query_refreshed()
    assert topology.topology_dirty_d.numpy()[0] == 0


def test_conflict_arbitration_and_sleep_refcount_guard() -> None:
    cfg = _config()
    topology = MikadoTopologyBuilder(cfg, device=str(_CUDA_DEVICE)).initialize()
    runtime = ECMRemodelRuntime(topology, max_refinement_level=cfg.max_refinement_level)
    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.segment_event_d, 0, RemodelEvent.REFINE)
    _set_i32(runtime.proposal.segment_event_d, 1, RemodelEvent.REFINE)
    runtime.prepare_candidate()
    status = runtime.plan.proposal_status_d.numpy()
    assert status[0] == ProposalStatus.SELECTED
    assert status[1] == ProposalStatus.DEFERRED_CONFLICT
    assert runtime.plan.event_totals_d.numpy()[0] == 1
    runtime.rollback(_accepted(0))

    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.segment_event_d, 0, RemodelEvent.REFINE)
    _set_i32(runtime.proposal.segment_event_d, 2, RemodelEvent.REFINE)
    runtime.prepare_candidate()
    runtime.validate_remap_acknowledgements()
    assert runtime.plan.event_totals_d.numpy()[0] == 2
    assert runtime.candidate_valid_d.numpy()[0] == 1
    runtime.commit_irreversible(_accepted(1), dt_phys=0.01, rng_seed=5)
    assert len(_assert_chain_and_bends(topology, 0)) == cfg.nodes_per_fiber - 1 + 2

    _set_i32(topology.endpoint_refcount_d, 0, 1)
    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.fiber_event_d, 0, FiberRemodelEvent.SLEEP)
    runtime.prepare_candidate()
    runtime.validate_remap_acknowledgements()
    assert runtime.candidate_valid_d.numpy()[0] == 0
    runtime.commit_irreversible(_accepted(1), dt_phys=0.01, rng_seed=6)
    assert topology.fiber_sleep_state_d.numpy()[0] == 0

    _set_i32(topology.endpoint_refcount_d, 0, 0)
    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.fiber_event_d, 0, FiberRemodelEvent.SLEEP)
    runtime.prepare_candidate()
    runtime.validate_remap_acknowledgements()
    assert runtime.candidate_valid_d.numpy()[0] == 1
    runtime.commit_irreversible(_accepted(1), dt_phys=0.01, rng_seed=7)
    assert topology.fiber_sleep_state_d.numpy()[0] == 1
    assert topology.topology_dirty_d.numpy()[0] == 1

    runtime.snapshot_candidate()
    _set_i32(runtime.proposal.fiber_event_d, 0, FiberRemodelEvent.WAKE)
    runtime.prepare_candidate()
    runtime.validate_remap_acknowledgements()
    runtime.commit_irreversible(_accepted(1), dt_phys=0.01, rng_seed=8)
    assert topology.fiber_sleep_state_d.numpy()[0] == 0
    assert runtime.allocated_array_bytes() > 0


def test_two_level_root_path_lineage_round_trip() -> None:
    cfg = _config()
    topology = MikadoTopologyBuilder(cfg, device=str(_CUDA_DEVICE)).initialize()
    runtime = ECMRemodelRuntime(topology, max_refinement_level=cfg.max_refinement_level)
    source = 1
    original_counts = topology.active_population_d.numpy().copy()
    original_rest = float(topology.segment_rest_length_d.numpy()[source])
    original_root = int(topology.refinement_root_segment_id_d.numpy()[source])
    original_next_id = int(topology.next_persistent_id_d.numpy()[0])

    def commit(event: RemodelEvent) -> None:
        runtime.snapshot_candidate()
        _set_i32(runtime.proposal.segment_event_d, source, event)
        runtime.prepare_candidate()
        runtime.validate_remap_acknowledgements()
        assert runtime.candidate_valid_d.numpy()[0] == 1
        runtime.commit_irreversible(_accepted(1), dt_phys=0.01, rng_seed=10)

    commit(RemodelEvent.REFINE)
    first_right = int(topology.segment_next_d.numpy()[source])
    assert topology.refinement_path_d.numpy()[source] == 0
    assert topology.refinement_path_d.numpy()[first_right] == 1
    commit(RemodelEvent.REFINE)
    second_right = int(topology.segment_next_d.numpy()[source])
    assert topology.segment_refinement_level_d.numpy()[source] == 2
    assert topology.segment_refinement_level_d.numpy()[second_right] == 2
    assert topology.refinement_path_d.numpy()[source] == 0
    assert topology.refinement_path_d.numpy()[second_right] == 1
    assert topology.refinement_root_segment_id_d.numpy()[source] == original_root
    assert topology.refinement_root_segment_id_d.numpy()[second_right] == original_root

    commit(RemodelEvent.COARSEN)
    assert topology.segment_refinement_level_d.numpy()[source] == 1
    assert topology.refinement_path_d.numpy()[source] == 0
    assert topology.segment_next_d.numpy()[source] == first_right
    commit(RemodelEvent.COARSEN)
    np.testing.assert_array_equal(topology.active_population_d.numpy(), original_counts)
    assert topology.segment_rest_length_d.numpy()[source] == pytest.approx(original_rest)
    assert topology.segment_refinement_level_d.numpy()[source] == 0
    assert topology.refinement_path_d.numpy()[source] == 0
    assert topology.refinement_root_segment_id_d.numpy()[source] == original_root
    assert topology.next_persistent_id_d.numpy()[0] == original_next_id + 8
    _assert_chain_and_bends(topology, 0)
