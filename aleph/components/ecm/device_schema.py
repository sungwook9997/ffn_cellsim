"""Authoritative Warp-CUDA schema for explicit collagen-network topology.

The schema deliberately contains no spring law, modulus, crosslink chemistry, or solver.
It is the ECM component's device-resident geometry/identity state and the future-safe
storage needed for accepted-step sleeping, damage, and refinement.  Graph-owned
connectors may borrow segment material points but never acquire these arrays.

Sanity Gate:
    * Identity: persistent fiber IDs never encode an array address; node/segment
      generations make stale connector endpoints detectable after accepted remeshing.
    * Population: active counts and allocated capacities are distinct.  Dormant slots
      remain allocated and therefore visible in the byte ledger.
    * Topology: every node, segment, and bend carries an owning fiber index; consecutive
      chain connectivity is initialized by the GPU Mikado builder.
    * Residency: every authoritative array is allocated on a resolved CUDA device.
      Host-visible shapes/capacities are configuration metadata, not per-step state.
    * Accounting: :meth:`allocated_array_bytes` sums Warp allocation capacities exactly;
      HashGrid/workspace and process-lifetime peak bytes remain separate native ledgers.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import warp as wp

__all__ = ["ECMTopologyState", "PopulationSlot", "require_cuda_device"]


class PopulationSlot(IntEnum):
    """Offsets in ``active_population_d`` for topology populations."""

    FIBER = 0
    NODE = 1
    SEGMENT = 2
    BEND = 3


def require_cuda_device(device: str | None) -> str:
    """Resolve ``device`` and reject every non-CUDA runtime."""
    wp.init()
    resolved = wp.get_device(device)
    if not resolved.is_cuda:
        raise RuntimeError(
            "ac.ecm topology requires NVIDIA Warp CUDA (I0-A); no CPU simulation "
            f"path exists. Resolved device {resolved!s}."
        )
    return str(resolved)


@dataclass(frozen=True, slots=True)
class ECMTopologyState:
    """Fixed-capacity, CUDA-owned topology for an explicit collagen Mikado network.

    Crosslink bonds are intentionally absent.  ``endpoint_refcount_d`` only prevents
    remeshing/sleeping under a graph-owned live endpoint.  Boundary tags classify
    topology; the boundary connector remains responsible for reaction and work.
    """

    device: str
    n_fiber_capacity: int
    initial_nodes_per_fiber: int
    node_slots_per_fiber: int
    position_d: wp.array
    reference_position_d: wp.array
    force_d: wp.array
    node_active_d: wp.array
    node_fiber_d: wp.array
    node_local_index_d: wp.array
    node_generation_d: wp.array
    persistent_node_id_d: wp.array
    boundary_face_mask_d: wp.array
    fiber_offsets_d: wp.array
    fiber_active_d: wp.array
    persistent_fiber_id_d: wp.array
    fiber_generation_d: wp.array
    fiber_sleep_state_d: wp.array
    segments_d: wp.array
    segment_fiber_d: wp.array
    segment_active_d: wp.array
    segment_generation_d: wp.array
    persistent_segment_id_d: wp.array
    refinement_parent_segment_id_d: wp.array
    refinement_root_segment_id_d: wp.array
    refinement_path_d: wp.array
    segment_prev_d: wp.array
    segment_next_d: wp.array
    segment_bend_d: wp.array
    segment_material_s0_d: wp.array
    segment_material_s1_d: wp.array
    segment_rest_length_d: wp.array
    segment_refinement_level_d: wp.array
    segment_damage_d: wp.array
    endpoint_refcount_d: wp.array
    bend_triples_d: wp.array
    bend_fiber_d: wp.array
    bend_active_d: wp.array
    bend_left_segment_d: wp.array
    bend_right_segment_d: wp.array
    node_free_list_d: wp.array
    segment_free_list_d: wp.array
    bend_free_list_d: wp.array
    free_population_d: wp.array
    active_population_d: wp.array
    next_persistent_id_d: wp.array
    topology_epoch_d: wp.array
    topology_dirty_d: wp.array

    @classmethod
    def allocate(
        cls,
        *,
        n_fibers: int,
        initial_nodes_per_fiber: int,
        node_slots_per_fiber: int,
        device: str | None = None,
    ) -> ECMTopologyState:
        """Allocate an empty fixed-capacity topology without host population arrays."""
        if isinstance(n_fibers, bool) or not isinstance(n_fibers, int) or n_fibers <= 0:
            raise ValueError("n_fibers must be a positive integer")
        if (
            isinstance(initial_nodes_per_fiber, bool)
            or not isinstance(initial_nodes_per_fiber, int)
            or initial_nodes_per_fiber < 2
        ):
            raise ValueError("initial_nodes_per_fiber must be an integer >= 2")
        if (
            isinstance(node_slots_per_fiber, bool)
            or not isinstance(node_slots_per_fiber, int)
            or node_slots_per_fiber < initial_nodes_per_fiber
        ):
            raise ValueError("node_slots_per_fiber must cover every initial node")
        dev = require_cuda_device(device)
        n_nodes = int(n_fibers) * int(node_slots_per_fiber)
        n_segments = int(n_fibers) * (int(node_slots_per_fiber) - 1)
        n_bends = int(n_fibers) * max(0, int(node_slots_per_fiber) - 2)
        with wp.ScopedDevice(dev):
            return cls(
                device=dev,
                n_fiber_capacity=int(n_fibers),
                initial_nodes_per_fiber=int(initial_nodes_per_fiber),
                node_slots_per_fiber=int(node_slots_per_fiber),
                position_d=wp.zeros(n_nodes, dtype=wp.vec3d),
                reference_position_d=wp.zeros(n_nodes, dtype=wp.vec3d),
                force_d=wp.zeros(n_nodes, dtype=wp.vec3d),
                node_active_d=wp.zeros(n_nodes, dtype=wp.int32),
                node_fiber_d=wp.zeros(n_nodes, dtype=wp.int32),
                node_local_index_d=wp.zeros(n_nodes, dtype=wp.int32),
                node_generation_d=wp.zeros(n_nodes, dtype=wp.int32),
                persistent_node_id_d=wp.full(n_nodes, -1, dtype=wp.int64),
                boundary_face_mask_d=wp.zeros(n_nodes, dtype=wp.int32),
                fiber_offsets_d=wp.zeros(int(n_fibers) + 1, dtype=wp.int32),
                fiber_active_d=wp.zeros(int(n_fibers), dtype=wp.int32),
                persistent_fiber_id_d=wp.zeros(int(n_fibers), dtype=wp.int64),
                fiber_generation_d=wp.zeros(int(n_fibers), dtype=wp.int32),
                fiber_sleep_state_d=wp.zeros(int(n_fibers), dtype=wp.int32),
                segments_d=wp.zeros((n_segments, 2), dtype=wp.int32),
                segment_fiber_d=wp.zeros(n_segments, dtype=wp.int32),
                segment_active_d=wp.zeros(n_segments, dtype=wp.int32),
                segment_generation_d=wp.zeros(n_segments, dtype=wp.int32),
                persistent_segment_id_d=wp.full(n_segments, -1, dtype=wp.int64),
                refinement_parent_segment_id_d=wp.full(n_segments, -1, dtype=wp.int64),
                refinement_root_segment_id_d=wp.full(n_segments, -1, dtype=wp.int64),
                refinement_path_d=wp.zeros(n_segments, dtype=wp.int64),
                segment_prev_d=wp.full(n_segments, -1, dtype=wp.int32),
                segment_next_d=wp.full(n_segments, -1, dtype=wp.int32),
                segment_bend_d=wp.full(n_segments, -1, dtype=wp.int32),
                segment_material_s0_d=wp.zeros(n_segments, dtype=wp.float64),
                segment_material_s1_d=wp.zeros(n_segments, dtype=wp.float64),
                segment_rest_length_d=wp.zeros(n_segments, dtype=wp.float64),
                segment_refinement_level_d=wp.zeros(n_segments, dtype=wp.int32),
                segment_damage_d=wp.zeros(n_segments, dtype=wp.float64),
                endpoint_refcount_d=wp.zeros(n_segments, dtype=wp.int32),
                bend_triples_d=wp.zeros((n_bends, 3), dtype=wp.int32),
                bend_fiber_d=wp.zeros(n_bends, dtype=wp.int32),
                bend_active_d=wp.zeros(n_bends, dtype=wp.int32),
                bend_left_segment_d=wp.full(n_bends, -1, dtype=wp.int32),
                bend_right_segment_d=wp.full(n_bends, -1, dtype=wp.int32),
                node_free_list_d=wp.zeros(n_nodes, dtype=wp.int32),
                segment_free_list_d=wp.zeros(n_segments, dtype=wp.int32),
                bend_free_list_d=wp.zeros(n_bends, dtype=wp.int32),
                free_population_d=wp.zeros(3, dtype=wp.int32),
                active_population_d=wp.zeros(len(PopulationSlot), dtype=wp.int32),
                next_persistent_id_d=wp.zeros(1, dtype=wp.int64),
                topology_epoch_d=wp.zeros(1, dtype=wp.int32),
                topology_dirty_d=wp.zeros(1, dtype=wp.int32),
            )

    @property
    def n_node_capacity(self) -> int:
        """Allocated node slots."""
        return self.n_fiber_capacity * self.node_slots_per_fiber

    @property
    def n_segment_capacity(self) -> int:
        """Allocated consecutive-segment slots."""
        return self.n_fiber_capacity * (self.node_slots_per_fiber - 1)

    @property
    def n_bend_capacity(self) -> int:
        """Allocated consecutive bending-triple slots."""
        return self.n_fiber_capacity * max(0, self.node_slots_per_fiber - 2)

    @property
    def n_initial_nodes(self) -> int:
        """Active nodes immediately after Mikado construction."""
        return self.n_fiber_capacity * self.initial_nodes_per_fiber

    @property
    def n_initial_segments(self) -> int:
        """Active segments immediately after Mikado construction."""
        return self.n_fiber_capacity * (self.initial_nodes_per_fiber - 1)

    @property
    def n_initial_bends(self) -> int:
        """Active bending triples immediately after Mikado construction."""
        return self.n_fiber_capacity * max(0, self.initial_nodes_per_fiber - 2)

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return every authoritative device array owned by the ECM topology."""
        return (
            self.position_d,
            self.reference_position_d,
            self.force_d,
            self.node_active_d,
            self.node_fiber_d,
            self.node_local_index_d,
            self.node_generation_d,
            self.persistent_node_id_d,
            self.boundary_face_mask_d,
            self.fiber_offsets_d,
            self.fiber_active_d,
            self.persistent_fiber_id_d,
            self.fiber_generation_d,
            self.fiber_sleep_state_d,
            self.segments_d,
            self.segment_fiber_d,
            self.segment_active_d,
            self.segment_generation_d,
            self.persistent_segment_id_d,
            self.refinement_parent_segment_id_d,
            self.refinement_root_segment_id_d,
            self.refinement_path_d,
            self.segment_prev_d,
            self.segment_next_d,
            self.segment_bend_d,
            self.segment_material_s0_d,
            self.segment_material_s1_d,
            self.segment_rest_length_d,
            self.segment_refinement_level_d,
            self.segment_damage_d,
            self.endpoint_refcount_d,
            self.bend_triples_d,
            self.bend_fiber_d,
            self.bend_active_d,
            self.bend_left_segment_d,
            self.bend_right_segment_d,
            self.node_free_list_d,
            self.segment_free_list_d,
            self.bend_free_list_d,
            self.free_population_d,
            self.active_population_d,
            self.next_persistent_id_d,
            self.topology_epoch_d,
            self.topology_dirty_d,
        )

    def allocated_array_bytes(self) -> int:
        """Return exact bytes reserved by owned Warp arrays from allocation metadata."""
        return sum(int(array.capacity) for array in self.owned_arrays())

    def capacity_ledger(self) -> dict[str, int]:
        """Return host metadata for the required active/allocated native ledger.

        Active counts remain in ``active_population_d`` and are read only by post-loop
        diagnostics; this method never relabels capacity as active population.
        """
        return {
            "ecm_fiber_capacity": self.n_fiber_capacity,
            "ecm_node_capacity": self.n_node_capacity,
            "ecm_segment_capacity": self.n_segment_capacity,
            "ecm_bend_capacity": self.n_bend_capacity,
            "ecm_initial_active_fibers": self.n_fiber_capacity,
            "ecm_initial_active_nodes": self.n_initial_nodes,
            "ecm_initial_active_segments": self.n_initial_segments,
            "ecm_initial_active_bends": self.n_initial_bends,
            "ecm_initial_dormant_nodes": self.n_node_capacity - self.n_initial_nodes,
            "ecm_initial_dormant_segments": self.n_segment_capacity - self.n_initial_segments,
            "ecm_initial_dormant_bends": self.n_bend_capacity - self.n_initial_bends,
            "ecm_owned_array_bytes": self.allocated_array_bytes(),
        }
