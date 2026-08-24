"""Accepted-step CUDA transactions for collagen ECM topology remodeling.

This runtime owns only ECM-internal topology.  It accepts externally generated event
and damage proposals, arbitrates non-overlapping chain edits, publishes a material-point
endpoint remap, and mutates topology only under the scheduler's final accepted predicate.
It contains no collagen modulus, force threshold, clutch kinetics, or acceptance policy.

Sanity Gate:
    * Rejection: proposal planning never mutates authoritative topology, so rejection is
      an O(proposal) device clear rather than a full-network D2D restore.
    * Conservation: refine splits material interval/rest length exactly; coarsen sums
      them and uses contour-length-weighted damage.  Active/free counts are inverse.
    * Identity: every created node/segment receives a fresh monotone int64 identity;
      slot generations increment on reuse/retirement, and root/path lineage survives
      repeated dyadic refine--coarsen cycles.
    * Connectivity: selected edits have a one-segment exclusion halo, preventing write
      races in shared bends and neighbor links.
    * Endpoints: commit is disabled until graph acknowledgements equal the pre-edit
      endpoint refcounts and conserve split/merge target counts.
    * Residency: selection, scans, planning, acknowledgement, commit, and rollback are
      CUDA-resident; no authoritative hot-loop scalar is copied to the host.
"""

from __future__ import annotations

import math

import warp as wp

from aleph.components.ecm.device_schema import ECMTopologyState, PopulationSlot
from aleph.components.ecm.remodel_schema import (
    ECMEndpointRemapView,
    ECMRemodelPlan,
    ECMRemodelProposal,
    FiberRemodelEvent,
    ProposalStatus,
    RemapKind,
    RemodelEvent,
)

__all__ = ["ECMRemodelRuntime"]


@wp.kernel
def _clear_segment_plan_kernel(
    proposal_status: wp.array(dtype=wp.int32),
    selected_event: wp.array(dtype=wp.int32),
    refine_flag: wp.array(dtype=wp.int32),
    coarsen_flag: wp.array(dtype=wp.int32),
    id_weight: wp.array(dtype=wp.int32),
    new_node_slot: wp.array(dtype=wp.int32),
    new_segment_slot: wp.array(dtype=wp.int32),
    new_bend_slot: wp.array(dtype=wp.int32),
    planned_node_id: wp.array(dtype=wp.int64),
    planned_segment_a_id: wp.array(dtype=wp.int64),
    planned_segment_b_id: wp.array(dtype=wp.int64),
    remap_kind: wp.array(dtype=wp.int32),
    remap_source_generation: wp.array(dtype=wp.int32),
    remap_source_id: wp.array(dtype=wp.int64),
    remap_target_a: wp.array(dtype=wp.int32),
    remap_target_b: wp.array(dtype=wp.int32),
    remap_target_a_generation: wp.array(dtype=wp.int32),
    remap_target_b_generation: wp.array(dtype=wp.int32),
    remap_target_a_id: wp.array(dtype=wp.int64),
    remap_target_b_id: wp.array(dtype=wp.int64),
    remap_ack: wp.array(dtype=wp.int32),
    remap_target_a_refcount: wp.array(dtype=wp.int32),
    remap_target_b_refcount: wp.array(dtype=wp.int32),
) -> None:
    segment = wp.tid()
    proposal_status[segment] = wp.int32(ProposalStatus.NONE)
    selected_event[segment] = wp.int32(RemodelEvent.NONE)
    refine_flag[segment] = wp.int32(0)
    coarsen_flag[segment] = wp.int32(0)
    id_weight[segment] = wp.int32(0)
    new_node_slot[segment] = wp.int32(-1)
    new_segment_slot[segment] = wp.int32(-1)
    new_bend_slot[segment] = wp.int32(-1)
    planned_node_id[segment] = wp.int64(-1)
    planned_segment_a_id[segment] = wp.int64(-1)
    planned_segment_b_id[segment] = wp.int64(-1)
    remap_kind[segment] = wp.int32(RemapKind.NONE)
    remap_source_generation[segment] = wp.int32(-1)
    remap_source_id[segment] = wp.int64(-1)
    remap_target_a[segment] = wp.int32(-1)
    remap_target_b[segment] = wp.int32(-1)
    remap_target_a_generation[segment] = wp.int32(-1)
    remap_target_b_generation[segment] = wp.int32(-1)
    remap_target_a_id[segment] = wp.int64(-1)
    remap_target_b_id[segment] = wp.int64(-1)
    remap_ack[segment] = wp.int32(0)
    remap_target_a_refcount[segment] = wp.int32(0)
    remap_target_b_refcount[segment] = wp.int32(0)


@wp.kernel
def _clear_plan_scalars_kernel(
    selected_fiber_event: wp.array(dtype=wp.int32),
    event_totals: wp.array(dtype=wp.int32),
    plan_valid: wp.array(dtype=wp.int32),
    remap_valid: wp.array(dtype=wp.int32),
    candidate_valid: wp.array(dtype=wp.int32),
) -> None:
    index = wp.tid()
    if index < selected_fiber_event.shape[0]:
        selected_fiber_event[index] = wp.int32(FiberRemodelEvent.NONE)
    if index < event_totals.shape[0]:
        event_totals[index] = wp.int32(0)
    if index == 0:
        plan_valid[0] = wp.int32(1)
        remap_valid[0] = wp.int32(0)
        candidate_valid[0] = wp.int32(0)


@wp.kernel
def _validate_raw_event_kernel(
    segment_event: wp.array(dtype=wp.int32),
    segment_active: wp.array(dtype=wp.int32),
    proposal_status: wp.array(dtype=wp.int32),
    plan_valid: wp.array(dtype=wp.int32),
) -> None:
    segment = wp.tid()
    event = segment_event[segment]
    recognized = (
        event == wp.int32(RemodelEvent.NONE)
        or event == wp.int32(RemodelEvent.REFINE)
        or event == wp.int32(RemodelEvent.COARSEN)
    )
    if not recognized or (event != wp.int32(RemodelEvent.NONE) and segment_active[segment] == 0):
        proposal_status[segment] = wp.int32(ProposalStatus.INVALID)
        wp.atomic_min(plan_valid, 0, wp.int32(0))


@wp.kernel
def _validate_damage_kernel(
    segment_active: wp.array(dtype=wp.int32),
    damage_increment: wp.array(dtype=wp.float64),
    plan_valid: wp.array(dtype=wp.int32),
) -> None:
    segment = wp.tid()
    increment = damage_increment[segment]
    if (
        not wp.isfinite(increment)
        or increment < wp.float64(0.0)
        or (segment_active[segment] == 0 and increment != wp.float64(0.0))
    ):
        wp.atomic_min(plan_valid, 0, wp.int32(0))


@wp.kernel
def _select_segment_events_kernel(
    segment_slots_per_fiber: wp.int32,
    max_refinement_level: wp.int32,
    segment_event: wp.array(dtype=wp.int32),
    segment_active: wp.array(dtype=wp.int32),
    segment_fiber: wp.array(dtype=wp.int32),
    segment_prev: wp.array(dtype=wp.int32),
    segment_next: wp.array(dtype=wp.int32),
    segment_level: wp.array(dtype=wp.int32),
    segment_root_id: wp.array(dtype=wp.int64),
    segment_path: wp.array(dtype=wp.int64),
    segments: wp.array2d(dtype=wp.int32),
    boundary_face_mask: wp.array(dtype=wp.int32),
    fiber_sleep_state: wp.array(dtype=wp.int32),
    proposal_status: wp.array(dtype=wp.int32),
    selected_event: wp.array(dtype=wp.int32),
    refine_flag: wp.array(dtype=wp.int32),
    coarsen_flag: wp.array(dtype=wp.int32),
    id_weight: wp.array(dtype=wp.int32),
    plan_valid: wp.array(dtype=wp.int32),
) -> None:
    """Greedily select every valid edit separated by one untouched segment."""
    fiber = wp.tid()
    base = fiber * segment_slots_per_fiber
    head = wp.int32(-1)
    for local in range(segment_slots_per_fiber):
        candidate = base + local
        if segment_active[candidate] != 0 and segment_prev[candidate] < 0:
            head = candidate

    segment = head
    blocked = wp.int32(0)
    for _step in range(segment_slots_per_fiber):
        if segment < 0:
            break
        event = segment_event[segment]
        recognized = (
            event == wp.int32(RemodelEvent.NONE)
            or event == wp.int32(RemodelEvent.REFINE)
            or event == wp.int32(RemodelEvent.COARSEN)
        )
        if not recognized or segment_active[segment] == 0:
            if event != wp.int32(RemodelEvent.NONE):
                proposal_status[segment] = wp.int32(ProposalStatus.INVALID)
                wp.atomic_min(plan_valid, 0, wp.int32(0))
        elif event != wp.int32(RemodelEvent.NONE) and fiber_sleep_state[fiber] != 0:
            proposal_status[segment] = wp.int32(ProposalStatus.INVALID)
            wp.atomic_min(plan_valid, 0, wp.int32(0))
        elif event != wp.int32(RemodelEvent.NONE) and blocked > 0:
            proposal_status[segment] = wp.int32(ProposalStatus.DEFERRED_CONFLICT)
        elif event == wp.int32(RemodelEvent.REFINE):
            if segment_level[segment] >= max_refinement_level:
                proposal_status[segment] = wp.int32(ProposalStatus.INVALID)
                wp.atomic_min(plan_valid, 0, wp.int32(0))
            else:
                proposal_status[segment] = wp.int32(ProposalStatus.SELECTED)
                selected_event[segment] = event
                refine_flag[segment] = wp.int32(1)
                id_weight[segment] = wp.int32(3)
                blocked = wp.int32(1)
        elif event == wp.int32(RemodelEvent.COARSEN):
            right = segment_next[segment]
            valid_pair = right >= 0
            if valid_pair:
                valid_pair = (
                    segment_active[right] != 0
                    and segment_fiber[right] == fiber
                    and segment_level[segment] > 0
                    and segment_level[right] == segment_level[segment]
                    and segment_root_id[segment] >= wp.int64(0)
                    and segment_root_id[right] == segment_root_id[segment]
                    and segment_path[segment] % wp.int64(2) == wp.int64(0)
                    and segment_path[right] == segment_path[segment] + wp.int64(1)
                    and segments[segment, 1] == segments[right, 0]
                    and boundary_face_mask[segments[segment, 1]] == 0
                )
            if not valid_pair:
                proposal_status[segment] = wp.int32(ProposalStatus.INVALID)
                wp.atomic_min(plan_valid, 0, wp.int32(0))
            else:
                proposal_status[segment] = wp.int32(ProposalStatus.SELECTED)
                selected_event[segment] = event
                coarsen_flag[segment] = wp.int32(1)
                id_weight[segment] = wp.int32(1)
                blocked = wp.int32(2)
        if blocked > 0 and selected_event[segment] == wp.int32(RemodelEvent.NONE):
            blocked = blocked - wp.int32(1)
        segment = segment_next[segment]


@wp.kernel
def _select_fiber_events_kernel(
    segment_slots_per_fiber: wp.int32,
    fiber_event: wp.array(dtype=wp.int32),
    fiber_sleep_state: wp.array(dtype=wp.int32),
    segment_active: wp.array(dtype=wp.int32),
    endpoint_refcount: wp.array(dtype=wp.int32),
    selected_event: wp.array(dtype=wp.int32),
    selected_fiber_event: wp.array(dtype=wp.int32),
    event_totals: wp.array(dtype=wp.int32),
    plan_valid: wp.array(dtype=wp.int32),
) -> None:
    fiber = wp.tid()
    event = fiber_event[fiber]
    has_endpoint = wp.int32(0)
    has_topology_event = wp.int32(0)
    base = fiber * segment_slots_per_fiber
    for local in range(segment_slots_per_fiber):
        segment = base + local
        if segment_active[segment] != 0 and endpoint_refcount[segment] > 0:
            has_endpoint = wp.int32(1)
        if selected_event[segment] != wp.int32(RemodelEvent.NONE):
            has_topology_event = wp.int32(1)
    if event == wp.int32(FiberRemodelEvent.NONE):
        pass
    elif event == wp.int32(FiberRemodelEvent.SLEEP):
        if fiber_sleep_state[fiber] != 0 or has_endpoint != 0 or has_topology_event != 0:
            wp.atomic_min(plan_valid, 0, wp.int32(0))
        else:
            selected_fiber_event[fiber] = event
            wp.atomic_add(event_totals, 3, wp.int32(1))
    elif event == wp.int32(FiberRemodelEvent.WAKE):
        if fiber_sleep_state[fiber] == 0 or has_topology_event != 0:
            wp.atomic_min(plan_valid, 0, wp.int32(0))
        else:
            selected_fiber_event[fiber] = event
            wp.atomic_add(event_totals, 3, wp.int32(1))
    else:
        wp.atomic_min(plan_valid, 0, wp.int32(0))


@wp.kernel
def _publish_scan_totals_kernel(
    refine_flag: wp.array(dtype=wp.int32),
    coarsen_flag: wp.array(dtype=wp.int32),
    id_weight: wp.array(dtype=wp.int32),
    refine_rank: wp.array(dtype=wp.int32),
    coarsen_rank: wp.array(dtype=wp.int32),
    id_offset: wp.array(dtype=wp.int32),
    old_free_population: wp.array(dtype=wp.int32),
    event_totals: wp.array(dtype=wp.int32),
    plan_valid: wp.array(dtype=wp.int32),
) -> None:
    if wp.tid() == 0:
        last = refine_flag.shape[0] - wp.int32(1)
        n_refine = refine_rank[last] + refine_flag[last]
        n_coarsen = coarsen_rank[last] + coarsen_flag[last]
        n_ids = id_offset[last] + id_weight[last]
        event_totals[0] = n_refine
        event_totals[1] = n_coarsen
        event_totals[2] = n_ids
        if (
            old_free_population[0] < n_refine
            or old_free_population[1] < n_refine
            or old_free_population[2] < n_refine
        ):
            plan_valid[0] = wp.int32(0)


@wp.kernel
def _assign_fiber_local_resources_kernel(
    node_slots_per_fiber: wp.int32,
    segment_slots_per_fiber: wp.int32,
    bend_slots_per_fiber: wp.int32,
    selected_event: wp.array(dtype=wp.int32),
    node_active: wp.array(dtype=wp.int32),
    segment_active: wp.array(dtype=wp.int32),
    bend_active: wp.array(dtype=wp.int32),
    new_node_slot: wp.array(dtype=wp.int32),
    new_segment_slot: wp.array(dtype=wp.int32),
    new_bend_slot: wp.array(dtype=wp.int32),
    plan_valid: wp.array(dtype=wp.int32),
) -> None:
    """Assign every refine resource from its owning fiber's fixed capacity block."""
    fiber = wp.tid()
    node_base = fiber * node_slots_per_fiber
    segment_base = fiber * segment_slots_per_fiber
    bend_base = fiber * bend_slots_per_fiber
    node_cursor = wp.int32(0)
    segment_cursor = wp.int32(0)
    bend_cursor = wp.int32(0)
    for local in range(segment_slots_per_fiber):
        segment = segment_base + local
        if selected_event[segment] == wp.int32(RemodelEvent.REFINE):
            node_slot = wp.int32(-1)
            segment_slot = wp.int32(-1)
            bend_slot = wp.int32(-1)
            for node_local in range(node_cursor, node_slots_per_fiber):
                candidate = node_base + node_local
                if node_slot < 0 and node_active[candidate] == 0:
                    node_slot = candidate
                    node_cursor = node_local + wp.int32(1)
            for segment_local in range(segment_cursor, segment_slots_per_fiber):
                candidate = segment_base + segment_local
                if segment_slot < 0 and segment_active[candidate] == 0:
                    segment_slot = candidate
                    segment_cursor = segment_local + wp.int32(1)
            for bend_local in range(bend_cursor, bend_slots_per_fiber):
                candidate = bend_base + bend_local
                if bend_slot < 0 and bend_active[candidate] == 0:
                    bend_slot = candidate
                    bend_cursor = bend_local + wp.int32(1)
            if node_slot < 0 or segment_slot < 0 or bend_slot < 0:
                wp.atomic_min(plan_valid, 0, wp.int32(0))
            else:
                new_node_slot[segment] = node_slot
                new_segment_slot[segment] = segment_slot
                new_bend_slot[segment] = bend_slot


@wp.kernel
def _materialize_plan_kernel(
    selected_event: wp.array(dtype=wp.int32),
    id_offset: wp.array(dtype=wp.int32),
    next_persistent_id: wp.array(dtype=wp.int64),
    segment_next: wp.array(dtype=wp.int32),
    segment_generation: wp.array(dtype=wp.int32),
    persistent_segment_id: wp.array(dtype=wp.int64),
    new_segment_slot: wp.array(dtype=wp.int32),
    planned_node_id: wp.array(dtype=wp.int64),
    planned_segment_a_id: wp.array(dtype=wp.int64),
    planned_segment_b_id: wp.array(dtype=wp.int64),
    remap_kind: wp.array(dtype=wp.int32),
    remap_source_generation: wp.array(dtype=wp.int32),
    remap_source_id: wp.array(dtype=wp.int64),
    remap_target_a: wp.array(dtype=wp.int32),
    remap_target_b: wp.array(dtype=wp.int32),
    remap_target_a_generation: wp.array(dtype=wp.int32),
    remap_target_b_generation: wp.array(dtype=wp.int32),
    remap_target_a_id: wp.array(dtype=wp.int64),
    remap_target_b_id: wp.array(dtype=wp.int64),
    plan_valid: wp.array(dtype=wp.int32),
) -> None:
    segment = wp.tid()
    event = selected_event[segment]
    if event == wp.int32(RemodelEvent.REFINE) and plan_valid[0] != 0:
        right_slot = new_segment_slot[segment]
        identity = next_persistent_id[0] + wp.int64(id_offset[segment])
        planned_node_id[segment] = identity
        planned_segment_a_id[segment] = identity + wp.int64(1)
        planned_segment_b_id[segment] = identity + wp.int64(2)
        remap_kind[segment] = wp.int32(RemapKind.SPLIT)
        remap_source_generation[segment] = segment_generation[segment]
        remap_source_id[segment] = persistent_segment_id[segment]
        remap_target_a[segment] = segment
        remap_target_b[segment] = right_slot
        remap_target_a_generation[segment] = segment_generation[segment] + wp.int32(1)
        remap_target_b_generation[segment] = segment_generation[right_slot] + wp.int32(1)
        remap_target_a_id[segment] = identity + wp.int64(1)
        remap_target_b_id[segment] = identity + wp.int64(2)
    elif event == wp.int32(RemodelEvent.COARSEN):
        right = segment_next[segment]
        identity = next_persistent_id[0] + wp.int64(id_offset[segment])
        planned_segment_a_id[segment] = identity
        remap_kind[segment] = wp.int32(RemapKind.MERGE_LEFT)
        remap_source_generation[segment] = segment_generation[segment]
        remap_source_id[segment] = persistent_segment_id[segment]
        remap_target_a[segment] = segment
        remap_target_a_generation[segment] = segment_generation[segment] + wp.int32(1)
        remap_target_a_id[segment] = identity
        remap_kind[right] = wp.int32(RemapKind.MERGE_RIGHT)
        remap_source_generation[right] = segment_generation[right]
        remap_source_id[right] = persistent_segment_id[right]
        remap_target_a[right] = segment
        remap_target_a_generation[right] = segment_generation[segment] + wp.int32(1)
        remap_target_a_id[right] = identity


@wp.kernel
def _validate_remap_kernel(
    remap_kind: wp.array(dtype=wp.int32),
    endpoint_refcount: wp.array(dtype=wp.int32),
    ack_count: wp.array(dtype=wp.int32),
    target_a_refcount: wp.array(dtype=wp.int32),
    target_b_refcount: wp.array(dtype=wp.int32),
    remap_valid: wp.array(dtype=wp.int32),
) -> None:
    segment = wp.tid()
    kind = remap_kind[segment]
    if kind != wp.int32(RemapKind.NONE):
        source_count = endpoint_refcount[segment]
        valid = ack_count[segment] == source_count
        if kind == wp.int32(RemapKind.SPLIT):
            valid = valid and (
                target_a_refcount[segment] + target_b_refcount[segment] == source_count
            )
        else:
            valid = valid and target_a_refcount[segment] == source_count
            valid = valid and target_b_refcount[segment] == 0
        if not valid:
            wp.atomic_min(remap_valid, 0, wp.int32(0))


@wp.kernel
def _publish_candidate_valid_kernel(
    plan_valid: wp.array(dtype=wp.int32),
    remap_valid: wp.array(dtype=wp.int32),
    candidate_valid: wp.array(dtype=wp.int32),
) -> None:
    if wp.tid() == 0:
        candidate_valid[0] = plan_valid[0] * remap_valid[0]


@wp.kernel
def _apply_damage_kernel(
    accepted: wp.array(dtype=wp.int32),
    candidate_valid: wp.array(dtype=wp.int32),
    segment_active: wp.array(dtype=wp.int32),
    damage_increment: wp.array(dtype=wp.float64),
    segment_damage: wp.array(dtype=wp.float64),
) -> None:
    segment = wp.tid()
    if accepted[0] != 0 and candidate_valid[0] != 0 and segment_active[segment] != 0:
        segment_damage[segment] = segment_damage[segment] + damage_increment[segment]


@wp.kernel
def _commit_topology_kernel(  # noqa: PLR0913
    accepted: wp.array(dtype=wp.int32),
    candidate_valid: wp.array(dtype=wp.int32),
    selected_event: wp.array(dtype=wp.int32),
    node_slots_per_fiber: wp.int32,
    new_node_slot: wp.array(dtype=wp.int32),
    new_segment_slot: wp.array(dtype=wp.int32),
    new_bend_slot: wp.array(dtype=wp.int32),
    planned_node_id: wp.array(dtype=wp.int64),
    planned_segment_a_id: wp.array(dtype=wp.int64),
    planned_segment_b_id: wp.array(dtype=wp.int64),
    remap_target_a_refcount: wp.array(dtype=wp.int32),
    remap_target_b_refcount: wp.array(dtype=wp.int32),
    position: wp.array(dtype=wp.vec3d),
    reference_position: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    node_active: wp.array(dtype=wp.int32),
    node_fiber: wp.array(dtype=wp.int32),
    node_local_index: wp.array(dtype=wp.int32),
    node_generation: wp.array(dtype=wp.int32),
    persistent_node_id: wp.array(dtype=wp.int64),
    boundary_face_mask: wp.array(dtype=wp.int32),
    segments: wp.array2d(dtype=wp.int32),
    segment_fiber: wp.array(dtype=wp.int32),
    segment_active: wp.array(dtype=wp.int32),
    segment_generation: wp.array(dtype=wp.int32),
    persistent_segment_id: wp.array(dtype=wp.int64),
    refinement_parent_id: wp.array(dtype=wp.int64),
    refinement_root_id: wp.array(dtype=wp.int64),
    refinement_path: wp.array(dtype=wp.int64),
    segment_prev: wp.array(dtype=wp.int32),
    segment_next: wp.array(dtype=wp.int32),
    segment_bend: wp.array(dtype=wp.int32),
    material_s0: wp.array(dtype=wp.float64),
    material_s1: wp.array(dtype=wp.float64),
    rest_length: wp.array(dtype=wp.float64),
    refinement_level: wp.array(dtype=wp.int32),
    damage: wp.array(dtype=wp.float64),
    endpoint_refcount: wp.array(dtype=wp.int32),
    bend_triples: wp.array2d(dtype=wp.int32),
    bend_fiber: wp.array(dtype=wp.int32),
    bend_active: wp.array(dtype=wp.int32),
    bend_left_segment: wp.array(dtype=wp.int32),
    bend_right_segment: wp.array(dtype=wp.int32),
) -> None:
    segment = wp.tid()
    if accepted[0] == 0 or candidate_valid[0] == 0:
        return
    event = selected_event[segment]
    if event == wp.int32(RemodelEvent.REFINE):
        node_a = segments[segment, 0]
        node_b = segments[segment, 1]
        midpoint = new_node_slot[segment]
        right = new_segment_slot[segment]
        internal_bend = new_bend_slot[segment]
        previous = segment_prev[segment]
        following = segment_next[segment]
        following_bend = segment_bend[segment]
        previous_bend = wp.int32(-1)
        if previous >= 0:
            previous_bend = segment_bend[previous]
        old_id = persistent_segment_id[segment]
        root_id = refinement_root_id[segment]
        path = refinement_path[segment]
        level = refinement_level[segment]
        s0 = material_s0[segment]
        s1 = material_s1[segment]
        half_s = wp.float64(0.5) * (s0 + s1)
        half_rest = wp.float64(0.5) * rest_length[segment]
        old_damage = damage[segment]

        position[midpoint] = wp.float64(0.5) * (position[node_a] + position[node_b])
        reference_position[midpoint] = wp.float64(0.5) * (
            reference_position[node_a] + reference_position[node_b]
        )
        force[midpoint] = wp.vec3d(0.0, 0.0, 0.0)
        node_active[midpoint] = wp.int32(1)
        node_fiber[midpoint] = segment_fiber[segment]
        node_local_index[midpoint] = midpoint % node_slots_per_fiber
        node_generation[midpoint] = node_generation[midpoint] + wp.int32(1)
        persistent_node_id[midpoint] = planned_node_id[segment]
        boundary_face_mask[midpoint] = boundary_face_mask[node_a] & boundary_face_mask[node_b]

        segments[segment, 1] = midpoint
        segment_generation[segment] = segment_generation[segment] + wp.int32(1)
        persistent_segment_id[segment] = planned_segment_a_id[segment]
        refinement_parent_id[segment] = old_id
        refinement_path[segment] = path * wp.int64(2)
        segment_next[segment] = right
        segment_bend[segment] = internal_bend
        material_s1[segment] = half_s
        rest_length[segment] = half_rest
        refinement_level[segment] = level + wp.int32(1)
        endpoint_refcount[segment] = remap_target_a_refcount[segment]

        segments[right, 0] = midpoint
        segments[right, 1] = node_b
        segment_fiber[right] = segment_fiber[segment]
        segment_active[right] = wp.int32(1)
        segment_generation[right] = segment_generation[right] + wp.int32(1)
        persistent_segment_id[right] = planned_segment_b_id[segment]
        refinement_parent_id[right] = old_id
        refinement_root_id[right] = root_id
        refinement_path[right] = path * wp.int64(2) + wp.int64(1)
        segment_prev[right] = segment
        segment_next[right] = following
        segment_bend[right] = following_bend
        material_s0[right] = half_s
        material_s1[right] = s1
        rest_length[right] = half_rest
        refinement_level[right] = level + wp.int32(1)
        damage[right] = old_damage
        endpoint_refcount[right] = remap_target_b_refcount[segment]
        refinement_root_id[segment] = root_id

        if following >= 0:
            segment_prev[following] = right
        if previous_bend >= 0:
            bend_triples[previous_bend, 2] = midpoint
        bend_triples[internal_bend, 0] = node_a
        bend_triples[internal_bend, 1] = midpoint
        bend_triples[internal_bend, 2] = node_b
        bend_fiber[internal_bend] = segment_fiber[segment]
        bend_active[internal_bend] = wp.int32(1)
        bend_left_segment[internal_bend] = segment
        bend_right_segment[internal_bend] = right
        if following_bend >= 0:
            bend_triples[following_bend, 0] = midpoint
            bend_left_segment[following_bend] = right

    elif event == wp.int32(RemodelEvent.COARSEN):
        right = segment_next[segment]
        node_a = segments[segment, 0]
        midpoint = segments[segment, 1]
        node_b = segments[right, 1]
        previous = segment_prev[segment]
        following = segment_next[right]
        internal_bend = segment_bend[segment]
        following_bend = segment_bend[right]
        previous_bend = wp.int32(-1)
        if previous >= 0:
            previous_bend = segment_bend[previous]
        left_rest = rest_length[segment]
        right_rest = rest_length[right]
        total_rest = left_rest + right_rest
        merged_damage = (damage[segment] * left_rest + damage[right] * right_rest) / total_rest
        segments[segment, 1] = node_b
        segment_generation[segment] = segment_generation[segment] + wp.int32(1)
        persistent_segment_id[segment] = planned_segment_a_id[segment]
        refinement_parent_id[segment] = wp.int64(-1)
        refinement_path[segment] = refinement_path[segment] // wp.int64(2)
        segment_next[segment] = following
        segment_bend[segment] = following_bend
        material_s1[segment] = material_s1[right]
        rest_length[segment] = total_rest
        refinement_level[segment] = refinement_level[segment] - wp.int32(1)
        damage[segment] = merged_damage
        endpoint_refcount[segment] = (
            remap_target_a_refcount[segment] + remap_target_a_refcount[right]
        )

        if following >= 0:
            segment_prev[following] = segment
        if previous_bend >= 0:
            bend_triples[previous_bend, 2] = node_b
        if following_bend >= 0:
            bend_triples[following_bend, 0] = node_a
            bend_left_segment[following_bend] = segment

        node_active[midpoint] = wp.int32(0)
        force[midpoint] = wp.vec3d(0.0, 0.0, 0.0)
        node_generation[midpoint] = node_generation[midpoint] + wp.int32(1)
        persistent_node_id[midpoint] = wp.int64(-1)
        boundary_face_mask[midpoint] = wp.int32(0)
        segment_active[right] = wp.int32(0)
        segment_generation[right] = segment_generation[right] + wp.int32(1)
        persistent_segment_id[right] = wp.int64(-1)
        refinement_parent_id[right] = wp.int64(-1)
        refinement_root_id[right] = wp.int64(-1)
        refinement_path[right] = wp.int64(0)
        segment_prev[right] = wp.int32(-1)
        segment_next[right] = wp.int32(-1)
        segment_bend[right] = wp.int32(-1)
        material_s0[right] = wp.float64(0.0)
        material_s1[right] = wp.float64(0.0)
        rest_length[right] = wp.float64(0.0)
        refinement_level[right] = wp.int32(0)
        damage[right] = wp.float64(0.0)
        endpoint_refcount[right] = wp.int32(0)
        bend_active[internal_bend] = wp.int32(0)
        bend_left_segment[internal_bend] = wp.int32(-1)
        bend_right_segment[internal_bend] = wp.int32(-1)


@wp.kernel
def _reset_free_population_kernel(
    accepted: wp.array(dtype=wp.int32),
    candidate_valid: wp.array(dtype=wp.int32),
    event_totals: wp.array(dtype=wp.int32),
    free_population: wp.array(dtype=wp.int32),
) -> None:
    if (
        wp.tid() == 0
        and accepted[0] != 0
        and candidate_valid[0] != 0
        and event_totals[0] + event_totals[1] > 0
    ):
        free_population[0] = wp.int32(0)
        free_population[1] = wp.int32(0)
        free_population[2] = wp.int32(0)


@wp.kernel
def _rebuild_free_list_kernel(
    accepted: wp.array(dtype=wp.int32),
    candidate_valid: wp.array(dtype=wp.int32),
    event_totals: wp.array(dtype=wp.int32),
    population_slot: wp.int32,
    active: wp.array(dtype=wp.int32),
    free_list: wp.array(dtype=wp.int32),
    free_population: wp.array(dtype=wp.int32),
) -> None:
    slot = wp.tid()
    if (
        accepted[0] != 0
        and candidate_valid[0] != 0
        and event_totals[0] + event_totals[1] > 0
        and active[slot] == 0
    ):
        free_index = wp.atomic_add(free_population, population_slot, wp.int32(1))
        free_list[free_index] = slot


@wp.kernel
def _commit_fiber_events_kernel(
    accepted: wp.array(dtype=wp.int32),
    candidate_valid: wp.array(dtype=wp.int32),
    selected_fiber_event: wp.array(dtype=wp.int32),
    fiber_sleep_state: wp.array(dtype=wp.int32),
) -> None:
    fiber = wp.tid()
    if accepted[0] != 0 and candidate_valid[0] != 0:
        event = selected_fiber_event[fiber]
        if event == wp.int32(FiberRemodelEvent.SLEEP):
            fiber_sleep_state[fiber] = wp.int32(1)
        elif event == wp.int32(FiberRemodelEvent.WAKE):
            fiber_sleep_state[fiber] = wp.int32(0)


@wp.kernel
def _commit_counters_kernel(
    accepted: wp.array(dtype=wp.int32),
    candidate_valid: wp.array(dtype=wp.int32),
    event_totals: wp.array(dtype=wp.int32),
    active_population: wp.array(dtype=wp.int32),
    next_persistent_id: wp.array(dtype=wp.int64),
    topology_epoch: wp.array(dtype=wp.int32),
    topology_dirty: wp.array(dtype=wp.int32),
) -> None:
    if wp.tid() == 0 and accepted[0] != 0 and candidate_valid[0] != 0:
        delta = event_totals[0] - event_totals[1]
        active_population[wp.int32(PopulationSlot.NODE)] = (
            active_population[wp.int32(PopulationSlot.NODE)] + delta
        )
        active_population[wp.int32(PopulationSlot.SEGMENT)] = (
            active_population[wp.int32(PopulationSlot.SEGMENT)] + delta
        )
        active_population[wp.int32(PopulationSlot.BEND)] = (
            active_population[wp.int32(PopulationSlot.BEND)] + delta
        )
        next_persistent_id[0] = next_persistent_id[0] + wp.int64(event_totals[2])
        if event_totals[0] + event_totals[1] + event_totals[3] > 0:
            topology_epoch[0] = topology_epoch[0] + wp.int32(1)
            topology_dirty[0] = wp.int32(1)


@wp.kernel
def _clear_proposal_after_kernel(
    clear_when_accepted: wp.int32,
    accepted: wp.array(dtype=wp.int32),
    segment_event: wp.array(dtype=wp.int32),
    damage_increment: wp.array(dtype=wp.float64),
) -> None:
    segment = wp.tid()
    should_clear = accepted[0] != 0
    if clear_when_accepted == 0:
        should_clear = not should_clear
    if should_clear:
        segment_event[segment] = wp.int32(RemodelEvent.NONE)
        damage_increment[segment] = wp.float64(0.0)


@wp.kernel
def _clear_fiber_proposal_after_kernel(
    clear_when_accepted: wp.int32,
    accepted: wp.array(dtype=wp.int32),
    fiber_event: wp.array(dtype=wp.int32),
) -> None:
    fiber = wp.tid()
    should_clear = accepted[0] != 0
    if clear_when_accepted == 0:
        should_clear = not should_clear
    if should_clear:
        fiber_event[fiber] = wp.int32(FiberRemodelEvent.NONE)


@wp.kernel
def _clear_topology_dirty_kernel(topology_dirty: wp.array(dtype=wp.int32)) -> None:
    if wp.tid() == 0:
        topology_dirty[0] = wp.int32(0)


class ECMRemodelRuntime:
    """GPU-resident proposal/plan/commit owner for one ECM topology state."""

    def __init__(self, topology: ECMTopologyState, *, max_refinement_level: int) -> None:
        if max_refinement_level < 0:
            raise ValueError("max_refinement_level must be nonnegative")
        self.topology = topology
        self.device = topology.device
        self.max_refinement_level = int(max_refinement_level)
        self.segment_slots_per_fiber = topology.n_segment_capacity // topology.n_fiber_capacity
        self.node_slots_per_fiber = topology.node_slots_per_fiber
        self.bend_slots_per_fiber = topology.n_bend_capacity // topology.n_fiber_capacity
        initial_segments = topology.initial_nodes_per_fiber - 1
        capacity_ratio = self.segment_slots_per_fiber // initial_segments
        supported_level = capacity_ratio.bit_length() - 1
        if (1 << supported_level) != capacity_ratio or max_refinement_level > supported_level:
            raise ValueError("max_refinement_level exceeds the topology's fixed per-fiber capacity")
        self.proposal = ECMRemodelProposal.allocate(
            n_segments=topology.n_segment_capacity,
            n_fibers=topology.n_fiber_capacity,
            device=self.device,
        )
        self.plan = ECMRemodelPlan.allocate(
            n_segments=topology.n_segment_capacity,
            n_fibers=topology.n_fiber_capacity,
            device=self.device,
        )

    @property
    def candidate_valid_d(self) -> wp.array:
        """Device scalar contributed to the global acceptance ledger."""
        return self.plan.candidate_valid_d

    @property
    def endpoint_remap(self) -> ECMEndpointRemapView:
        """Borrow the graph-facing endpoint remap/ack view."""
        return self.plan.remap

    def _clear_plan(self) -> None:
        p = self.plan
        r = p.remap
        wp.launch(
            _clear_segment_plan_kernel,
            dim=self.topology.n_segment_capacity,
            inputs=[
                p.proposal_status_d,
                p.selected_event_d,
                p.refine_flag_d,
                p.coarsen_flag_d,
                p.id_weight_d,
                p.new_node_slot_d,
                p.new_segment_slot_d,
                p.new_bend_slot_d,
                p.planned_node_id_d,
                p.planned_segment_a_id_d,
                p.planned_segment_b_id_d,
                r.kind_d,
                r.source_generation_d,
                r.source_persistent_id_d,
                r.target_a_segment_d,
                r.target_b_segment_d,
                r.target_a_generation_d,
                r.target_b_generation_d,
                r.target_a_persistent_id_d,
                r.target_b_persistent_id_d,
                r.ack_count_d,
                r.target_a_refcount_d,
                r.target_b_refcount_d,
            ],
            device=self.device,
        )
        wp.launch(
            _clear_plan_scalars_kernel,
            dim=max(self.topology.n_fiber_capacity, 4),
            inputs=[
                p.selected_fiber_event_d,
                p.event_totals_d,
                p.plan_valid_d,
                p.remap_valid_d,
                p.candidate_valid_d,
            ],
            device=self.device,
        )

    def snapshot_candidate(self) -> None:
        """Open a proposal transaction without copying unchanged ECM topology.

        Authoritative arrays are commit-only, so the accepted topology itself is already
        the rollback snapshot.  This resets stale proposal/plan state entirely on CUDA.
        """
        self._clear_plan()
        self.proposal.segment_event_d.zero_()
        self.proposal.damage_increment_d.zero_()
        self.proposal.fiber_event_d.zero_()

    def prepare_candidate(self) -> None:
        """Arbitrate proposals and emit deterministic resource/remap plans on CUDA."""
        t = self.topology
        p = self.plan
        q = self.proposal
        self._clear_plan()
        wp.copy(p.old_free_population_d, t.free_population_d)
        wp.launch(
            _validate_raw_event_kernel,
            dim=t.n_segment_capacity,
            inputs=[
                q.segment_event_d,
                t.segment_active_d,
                p.proposal_status_d,
                p.plan_valid_d,
            ],
            device=self.device,
        )
        wp.launch(
            _validate_damage_kernel,
            dim=t.n_segment_capacity,
            inputs=[t.segment_active_d, q.damage_increment_d, p.plan_valid_d],
            device=self.device,
        )
        wp.launch(
            _select_segment_events_kernel,
            dim=t.n_fiber_capacity,
            inputs=[
                wp.int32(self.segment_slots_per_fiber),
                wp.int32(self.max_refinement_level),
                q.segment_event_d,
                t.segment_active_d,
                t.segment_fiber_d,
                t.segment_prev_d,
                t.segment_next_d,
                t.segment_refinement_level_d,
                t.refinement_root_segment_id_d,
                t.refinement_path_d,
                t.segments_d,
                t.boundary_face_mask_d,
                t.fiber_sleep_state_d,
                p.proposal_status_d,
                p.selected_event_d,
                p.refine_flag_d,
                p.coarsen_flag_d,
                p.id_weight_d,
                p.plan_valid_d,
            ],
            device=self.device,
        )
        wp.launch(
            _select_fiber_events_kernel,
            dim=t.n_fiber_capacity,
            inputs=[
                wp.int32(self.segment_slots_per_fiber),
                q.fiber_event_d,
                t.fiber_sleep_state_d,
                t.segment_active_d,
                t.endpoint_refcount_d,
                p.selected_event_d,
                p.selected_fiber_event_d,
                p.event_totals_d,
                p.plan_valid_d,
            ],
            device=self.device,
        )
        wp.utils.array_scan(p.refine_flag_d, p.refine_rank_d, inclusive=False)
        wp.utils.array_scan(p.coarsen_flag_d, p.coarsen_rank_d, inclusive=False)
        wp.utils.array_scan(p.id_weight_d, p.id_offset_d, inclusive=False)
        wp.launch(
            _publish_scan_totals_kernel,
            dim=1,
            inputs=[
                p.refine_flag_d,
                p.coarsen_flag_d,
                p.id_weight_d,
                p.refine_rank_d,
                p.coarsen_rank_d,
                p.id_offset_d,
                p.old_free_population_d,
                p.event_totals_d,
                p.plan_valid_d,
            ],
            device=self.device,
        )
        wp.launch(
            _assign_fiber_local_resources_kernel,
            dim=t.n_fiber_capacity,
            inputs=[
                wp.int32(self.node_slots_per_fiber),
                wp.int32(self.segment_slots_per_fiber),
                wp.int32(self.bend_slots_per_fiber),
                p.selected_event_d,
                t.node_active_d,
                t.segment_active_d,
                t.bend_active_d,
                p.new_node_slot_d,
                p.new_segment_slot_d,
                p.new_bend_slot_d,
                p.plan_valid_d,
            ],
            device=self.device,
        )
        wp.launch(
            _materialize_plan_kernel,
            dim=t.n_segment_capacity,
            inputs=[
                p.selected_event_d,
                p.id_offset_d,
                t.next_persistent_id_d,
                t.segment_next_d,
                t.segment_generation_d,
                t.persistent_segment_id_d,
                p.new_segment_slot_d,
                p.planned_node_id_d,
                p.planned_segment_a_id_d,
                p.planned_segment_b_id_d,
                p.remap.kind_d,
                p.remap.source_generation_d,
                p.remap.source_persistent_id_d,
                p.remap.target_a_segment_d,
                p.remap.target_b_segment_d,
                p.remap.target_a_generation_d,
                p.remap.target_b_generation_d,
                p.remap.target_a_persistent_id_d,
                p.remap.target_b_persistent_id_d,
                p.plan_valid_d,
            ],
            device=self.device,
        )

    def validate_remap_acknowledgements(self) -> None:
        """Publish endpoint-conservation validity after graph-owned remapping."""
        self.plan.remap_valid_d.fill_(1)
        r = self.plan.remap
        wp.launch(
            _validate_remap_kernel,
            dim=self.topology.n_segment_capacity,
            inputs=[
                r.kind_d,
                self.topology.endpoint_refcount_d,
                r.ack_count_d,
                r.target_a_refcount_d,
                r.target_b_refcount_d,
                self.plan.remap_valid_d,
            ],
            device=self.device,
        )
        wp.launch(
            _publish_candidate_valid_kernel,
            dim=1,
            inputs=[
                self.plan.plan_valid_d,
                self.plan.remap_valid_d,
                self.plan.candidate_valid_d,
            ],
            device=self.device,
        )

    def rollback(self, accepted: wp.array) -> None:
        """Discard rejected proposals; authoritative topology was never mutated."""
        self._clear_proposals_after(accepted, clear_when_accepted=False)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Apply damage/topology/sleep edits only under accepted and valid predicates."""
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        if rng_seed < 0:
            raise ValueError("rng_seed must be nonnegative")
        t = self.topology
        p = self.plan
        r = p.remap
        wp.launch(
            _apply_damage_kernel,
            dim=t.n_segment_capacity,
            inputs=[
                accepted,
                p.candidate_valid_d,
                t.segment_active_d,
                self.proposal.damage_increment_d,
                t.segment_damage_d,
            ],
            device=self.device,
        )
        wp.launch(
            _commit_topology_kernel,
            dim=t.n_segment_capacity,
            inputs=[
                accepted,
                p.candidate_valid_d,
                p.selected_event_d,
                wp.int32(self.node_slots_per_fiber),
                p.new_node_slot_d,
                p.new_segment_slot_d,
                p.new_bend_slot_d,
                p.planned_node_id_d,
                p.planned_segment_a_id_d,
                p.planned_segment_b_id_d,
                r.target_a_refcount_d,
                r.target_b_refcount_d,
                t.position_d,
                t.reference_position_d,
                t.force_d,
                t.node_active_d,
                t.node_fiber_d,
                t.node_local_index_d,
                t.node_generation_d,
                t.persistent_node_id_d,
                t.boundary_face_mask_d,
                t.segments_d,
                t.segment_fiber_d,
                t.segment_active_d,
                t.segment_generation_d,
                t.persistent_segment_id_d,
                t.refinement_parent_segment_id_d,
                t.refinement_root_segment_id_d,
                t.refinement_path_d,
                t.segment_prev_d,
                t.segment_next_d,
                t.segment_bend_d,
                t.segment_material_s0_d,
                t.segment_material_s1_d,
                t.segment_rest_length_d,
                t.segment_refinement_level_d,
                t.segment_damage_d,
                t.endpoint_refcount_d,
                t.bend_triples_d,
                t.bend_fiber_d,
                t.bend_active_d,
                t.bend_left_segment_d,
                t.bend_right_segment_d,
            ],
            device=self.device,
        )
        wp.launch(
            _reset_free_population_kernel,
            dim=1,
            inputs=[accepted, p.candidate_valid_d, p.event_totals_d, t.free_population_d],
            device=self.device,
        )
        for population_slot, active, free_list in (
            (0, t.node_active_d, t.node_free_list_d),
            (1, t.segment_active_d, t.segment_free_list_d),
            (2, t.bend_active_d, t.bend_free_list_d),
        ):
            if int(active.shape[0]) == 0:
                continue
            wp.launch(
                _rebuild_free_list_kernel,
                dim=int(active.shape[0]),
                inputs=[
                    accepted,
                    p.candidate_valid_d,
                    p.event_totals_d,
                    wp.int32(population_slot),
                    active,
                    free_list,
                    t.free_population_d,
                ],
                device=self.device,
            )
        wp.launch(
            _commit_fiber_events_kernel,
            dim=t.n_fiber_capacity,
            inputs=[
                accepted,
                p.candidate_valid_d,
                p.selected_fiber_event_d,
                t.fiber_sleep_state_d,
            ],
            device=self.device,
        )
        wp.launch(
            _commit_counters_kernel,
            dim=1,
            inputs=[
                accepted,
                p.candidate_valid_d,
                p.event_totals_d,
                t.active_population_d,
                t.next_persistent_id_d,
                t.topology_epoch_d,
                t.topology_dirty_d,
            ],
            device=self.device,
        )
        self._clear_proposals_after(accepted, clear_when_accepted=True)

    def _clear_proposals_after(self, accepted: wp.array, *, clear_when_accepted: bool) -> None:
        flag = wp.int32(1 if clear_when_accepted else 0)
        wp.launch(
            _clear_proposal_after_kernel,
            dim=self.topology.n_segment_capacity,
            inputs=[
                flag,
                accepted,
                self.proposal.segment_event_d,
                self.proposal.damage_increment_d,
            ],
            device=self.device,
        )
        wp.launch(
            _clear_fiber_proposal_after_kernel,
            dim=self.topology.n_fiber_capacity,
            inputs=[flag, accepted, self.proposal.fiber_event_d],
            device=self.device,
        )

    def mark_candidate_query_refreshed(self) -> None:
        """Clear the device dirty bit after graph-owned broadphase refresh."""
        wp.launch(
            _clear_topology_dirty_kernel,
            dim=1,
            inputs=[self.topology.topology_dirty_d],
            device=self.device,
        )

    def allocated_array_bytes(self) -> int:
        """Return proposal+plan bytes, excluding topology-owned storage."""
        arrays = (*self.proposal.owned_arrays(), *self.plan.owned_arrays())
        return sum(int(array.capacity) for array in arrays)
