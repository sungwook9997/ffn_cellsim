"""Device-only proposal and endpoint-remap schema for ECM topology transactions.

The arrays in this module contain topology intent, never collagen constitutive laws.
Stress, strain, damage thresholds, clutch chemistry, and acceptance remain external
owners.  An external CUDA producer writes proposals; the ECM planner emits a stable
material-point remap table for graph-owned endpoints before accepted commit.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import warp as wp

from aleph.components.ecm.device_schema import require_cuda_device

__all__ = [
    "ECMEndpointRemapView",
    "ECMRemodelPlan",
    "ECMRemodelProposal",
    "FiberRemodelEvent",
    "ProposalStatus",
    "RemapKind",
    "RemodelEvent",
]


class RemodelEvent(IntEnum):
    """Externally requested per-segment topology operation."""

    NONE = 0
    REFINE = 1
    COARSEN = 2


class FiberRemodelEvent(IntEnum):
    """Externally requested per-fiber query-activity operation."""

    NONE = 0
    SLEEP = 1
    WAKE = 2


class ProposalStatus(IntEnum):
    """Planner disposition for each raw segment proposal."""

    NONE = 0
    SELECTED = 1
    DEFERRED_CONFLICT = 2
    INVALID = 3


class RemapKind(IntEnum):
    """Transformation of a graph-owned endpoint's material coordinate."""

    NONE = 0
    SPLIT = 1
    MERGE_LEFT = 2
    MERGE_RIGHT = 3


@dataclass(frozen=True, slots=True)
class ECMRemodelProposal:
    """CUDA proposal buffers written by damage/contact/KMC owners.

    ``damage_increment_d`` is a nonnegative accepted-step increment.  This schema does
    not decide when damage occurs and intentionally defines no production threshold.
    """

    device: str
    segment_event_d: wp.array
    damage_increment_d: wp.array
    fiber_event_d: wp.array

    @classmethod
    def allocate(
        cls,
        *,
        n_segments: int,
        n_fibers: int,
        device: str | None,
    ) -> ECMRemodelProposal:
        """Allocate empty proposal buffers on the required CUDA device."""
        if n_segments <= 0 or n_fibers <= 0:
            raise ValueError("ECM remodeling capacities must be positive")
        dev = require_cuda_device(device)
        return cls(
            device=dev,
            segment_event_d=wp.zeros(n_segments, dtype=wp.int32, device=dev),
            damage_increment_d=wp.zeros(n_segments, dtype=wp.float64, device=dev),
            fiber_event_d=wp.zeros(n_fibers, dtype=wp.int32, device=dev),
        )

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return proposal buffers for byte accounting."""
        return (self.segment_event_d, self.damage_increment_d, self.fiber_event_d)


@dataclass(frozen=True, slots=True)
class ECMEndpointRemapView:
    """Borrowed CUDA handoff from ECM topology to the connector graph.

    The graph reads source/target identity and generation, remaps each live endpoint,
    then atomically increments ``ack_count_d`` and the applicable target count.  ECM
    commit is disabled unless acknowledgements conserve every source refcount.
    """

    kind_d: wp.array
    source_generation_d: wp.array
    source_persistent_id_d: wp.array
    target_a_segment_d: wp.array
    target_b_segment_d: wp.array
    target_a_generation_d: wp.array
    target_b_generation_d: wp.array
    target_a_persistent_id_d: wp.array
    target_b_persistent_id_d: wp.array
    ack_count_d: wp.array
    target_a_refcount_d: wp.array
    target_b_refcount_d: wp.array


@dataclass(frozen=True, slots=True)
class ECMRemodelPlan:
    """Fixed-capacity CUDA scratch for deterministic remodeling arbitration."""

    proposal_status_d: wp.array
    selected_event_d: wp.array
    refine_flag_d: wp.array
    coarsen_flag_d: wp.array
    id_weight_d: wp.array
    refine_rank_d: wp.array
    coarsen_rank_d: wp.array
    id_offset_d: wp.array
    new_node_slot_d: wp.array
    new_segment_slot_d: wp.array
    new_bend_slot_d: wp.array
    planned_node_id_d: wp.array
    planned_segment_a_id_d: wp.array
    planned_segment_b_id_d: wp.array
    selected_fiber_event_d: wp.array
    old_free_population_d: wp.array
    event_totals_d: wp.array
    plan_valid_d: wp.array
    remap_valid_d: wp.array
    candidate_valid_d: wp.array
    remap: ECMEndpointRemapView

    @classmethod
    def allocate(
        cls,
        *,
        n_segments: int,
        n_fibers: int,
        device: str | None,
    ) -> ECMRemodelPlan:
        """Allocate planner, scan, resource, and graph-handoff arrays."""
        if n_segments <= 0 or n_fibers <= 0:
            raise ValueError("ECM remodeling capacities must be positive")
        dev = require_cuda_device(device)

        def i32(fill: int = 0) -> wp.array:
            return wp.full(n_segments, fill, dtype=wp.int32, device=dev)

        def i64(fill: int = 0) -> wp.array:
            return wp.full(n_segments, fill, dtype=wp.int64, device=dev)

        remap = ECMEndpointRemapView(
            kind_d=i32(),
            source_generation_d=i32(-1),
            source_persistent_id_d=i64(-1),
            target_a_segment_d=i32(-1),
            target_b_segment_d=i32(-1),
            target_a_generation_d=i32(-1),
            target_b_generation_d=i32(-1),
            target_a_persistent_id_d=i64(-1),
            target_b_persistent_id_d=i64(-1),
            ack_count_d=i32(),
            target_a_refcount_d=i32(),
            target_b_refcount_d=i32(),
        )
        return cls(
            proposal_status_d=i32(),
            selected_event_d=i32(),
            refine_flag_d=i32(),
            coarsen_flag_d=i32(),
            id_weight_d=i32(),
            refine_rank_d=i32(),
            coarsen_rank_d=i32(),
            id_offset_d=i32(),
            new_node_slot_d=i32(-1),
            new_segment_slot_d=i32(-1),
            new_bend_slot_d=i32(-1),
            planned_node_id_d=i64(-1),
            planned_segment_a_id_d=i64(-1),
            planned_segment_b_id_d=i64(-1),
            selected_fiber_event_d=wp.zeros(n_fibers, dtype=wp.int32, device=dev),
            old_free_population_d=wp.zeros(3, dtype=wp.int32, device=dev),
            event_totals_d=wp.zeros(4, dtype=wp.int32, device=dev),
            plan_valid_d=wp.ones(1, dtype=wp.int32, device=dev),
            remap_valid_d=wp.ones(1, dtype=wp.int32, device=dev),
            candidate_valid_d=wp.ones(1, dtype=wp.int32, device=dev),
            remap=remap,
        )

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return all scratch/handoff arrays for exact byte accounting."""
        remap = self.remap
        return (
            self.proposal_status_d,
            self.selected_event_d,
            self.refine_flag_d,
            self.coarsen_flag_d,
            self.id_weight_d,
            self.refine_rank_d,
            self.coarsen_rank_d,
            self.id_offset_d,
            self.new_node_slot_d,
            self.new_segment_slot_d,
            self.new_bend_slot_d,
            self.planned_node_id_d,
            self.planned_segment_a_id_d,
            self.planned_segment_b_id_d,
            self.selected_fiber_event_d,
            self.old_free_population_d,
            self.event_totals_d,
            self.plan_valid_d,
            self.remap_valid_d,
            self.candidate_valid_d,
            remap.kind_d,
            remap.source_generation_d,
            remap.source_persistent_id_d,
            remap.target_a_segment_d,
            remap.target_b_segment_d,
            remap.target_a_generation_d,
            remap.target_b_generation_d,
            remap.target_a_persistent_id_d,
            remap.target_b_persistent_id_d,
            remap.ack_count_d,
            remap.target_a_refcount_d,
            remap.target_b_refcount_d,
        )
