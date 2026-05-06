"""Cell cluster state schema.

Canonical Phase 1 ``CellClusterState`` per
``docs/v2/v2_phase1_plan_consolidated.md`` §5.6. A cluster aggregates
single-cell states with a shared ECM substrate and an explicit junction
list; cells remain separate polygons and never merge geometry on
contact.

Cluster-level dynamics (contact detection, junction creation, contact
inhibition, neighbor reorientation) are explicitly out of P0 scope.
This schema only validates structural consistency: every junction's
cell ids point at known cluster members, the symmetric neighbor graph
matches `SingleCellState.neighbor_cell_ids`, and lineage references
(parent ids) are either resolvable within the cluster or marked as
external.

Sanity Gate scope: non-physics schema module. Full 6-item gate N/A;
this module owns the boundary-case checks (cell id consistency,
junction → member resolution, symmetric neighbor graph, no self-loops,
unique junction ids) and the measurement-protocol consistency that
``role_labels`` use a closed enum so downstream visualization can color
by role without inventing categories.

Magic-Number Block: no tunable numeric. N/A.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional

from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.junction import JunctionState
from acs.v2.single_cell import SingleCellState

ClusterRole = Literal["leader", "follower", "internal", "isolated"]


@dataclass
class CellClusterState:
    """Schema for a small cluster of single cells on a shared ECM substrate."""

    cells: dict[str, SingleCellState]
    ecm: ECMSubstrateState
    junctions: tuple[JunctionState, ...] = ()
    role_labels: dict[str, ClusterRole] = field(default_factory=dict)
    external_parent_ids: tuple[str, ...] = ()

    def validate(self) -> None:
        if not self.cells:
            raise ValueError("CellClusterState.cells must contain at least one cell")
        for cell_id, cell in self.cells.items():
            if not isinstance(cell_id, str) or not cell_id.strip():
                raise ValueError("CellClusterState.cells keys must be non-empty strings")
            if cell.cell_id != cell_id:
                raise ValueError(
                    f"cluster key {cell_id!r} disagrees with cell.cell_id "
                    f"{cell.cell_id!r}"
                )
            cell.validate()

        self.ecm.validate()

        seen_junction_ids: set[str] = set()
        for junction in self.junctions:
            junction.validate()
            if junction.junction_id in seen_junction_ids:
                raise ValueError(
                    f"duplicate junction_id in cluster: {junction.junction_id!r}"
                )
            seen_junction_ids.add(junction.junction_id)
            for jid in (junction.cell_id_a, junction.cell_id_b):
                if jid not in self.cells:
                    raise ValueError(
                        f"junction {junction.junction_id!r} references unknown cell {jid!r}"
                    )

        external_ids = set(self.external_parent_ids)
        for ext_id in self.external_parent_ids:
            if not isinstance(ext_id, str) or not ext_id.strip():
                raise ValueError("external_parent_ids must contain non-empty strings")
        if len(external_ids) != len(self.external_parent_ids):
            raise ValueError("external_parent_ids must be unique")
        cell_ids = set(self.cells.keys())
        if external_ids & cell_ids:
            raise ValueError(
                f"external_parent_ids overlap with cluster cell ids: "
                f"{sorted(external_ids & cell_ids)!r}"
            )

        for cell_id, cell in self.cells.items():
            for neighbor in cell.neighbor_cell_ids:
                if neighbor not in self.cells:
                    raise ValueError(
                        f"cell {cell_id!r} declares neighbor {neighbor!r} "
                        f"not in cluster"
                    )
                neighbor_cell = self.cells[neighbor]
                if cell_id not in neighbor_cell.neighbor_cell_ids:
                    raise ValueError(
                        f"asymmetric neighbor link: {cell_id!r} → {neighbor!r} "
                        f"is not mirrored by {neighbor!r}"
                    )
            if cell.parent_cell_id is not None:
                if (
                    cell.parent_cell_id not in self.cells
                    and cell.parent_cell_id not in external_ids
                ):
                    raise ValueError(
                        f"cell {cell_id!r} parent_cell_id {cell.parent_cell_id!r} "
                        f"is neither in cluster nor declared in external_parent_ids"
                    )

        for cell_id, role in self.role_labels.items():
            if cell_id not in self.cells:
                raise ValueError(
                    f"role_labels references unknown cell {cell_id!r}"
                )
            if role not in {"leader", "follower", "internal", "isolated"}:
                raise ValueError(
                    f"role for {cell_id!r} must be leader, follower, internal, "
                    f"or isolated, got {role!r}"
                )

    def neighbor_graph_edges(self) -> tuple[tuple[str, str], ...]:
        """Return canonicalized neighbor graph edges (lexicographic pairs)."""

        edges: set[tuple[str, str]] = set()
        for cell_id, cell in self.cells.items():
            for neighbor in cell.neighbor_cell_ids:
                a, b = sorted((cell_id, neighbor))
                edges.add((a, b))
        return tuple(sorted(edges))

    @property
    def cell_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self.cells.keys()))

    def get_role(self, cell_id: str) -> Optional[ClusterRole]:
        return self.role_labels.get(cell_id)
