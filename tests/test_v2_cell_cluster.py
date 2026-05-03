from __future__ import annotations

import numpy as np
import pytest

from acs.v2.cell_cluster import CellClusterState
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.junction import JunctionState
from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.single_cell import SingleCellState


def _boundary() -> MeasurementBoundary:
    return MeasurementBoundary.from_array(
        np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )


def _ecm() -> ECMSubstrateState:
    nx, ny = 2, 2
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    return ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=1.0,
        stiffness_kpa=np.ones((nx, ny), dtype=np.float64),
        ligand_density=np.zeros((nx, ny), dtype=np.float64),
        fiber_density=np.zeros((nx, ny), dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
    )


def _cell(cell_id: str, **overrides) -> SingleCellState:
    base = dict(
        cell_id=cell_id,
        time_s=0.0,
        measurement_boundary=_boundary(),
    )
    base.update(overrides)
    return SingleCellState(**base)


def test_minimal_cluster_validates():
    cluster = CellClusterState(
        cells={"cell-1": _cell("cell-1")},
        ecm=_ecm(),
    )
    cluster.validate()
    assert cluster.cell_ids == ("cell-1",)


def test_cluster_key_disagreement_rejected():
    cluster = CellClusterState(
        cells={"cell-1": _cell("cell-2")},
        ecm=_ecm(),
    )
    with pytest.raises(ValueError, match="disagrees"):
        cluster.validate()


def test_empty_cluster_rejected():
    cluster = CellClusterState(cells={}, ecm=_ecm())
    with pytest.raises(ValueError, match="at least one cell"):
        cluster.validate()


def test_junction_unknown_cell_rejected():
    cluster = CellClusterState(
        cells={"cell-1": _cell("cell-1")},
        ecm=_ecm(),
        junctions=(
            JunctionState(
                junction_id="j-1",
                cell_id_a="cell-1",
                cell_id_b="cell-2",
                contact_length_um=1.0,
                age_s=0.0,
                maturity=0.1,
            ),
        ),
    )
    with pytest.raises(ValueError, match="unknown cell"):
        cluster.validate()


def test_duplicate_junction_id_rejected():
    cluster = CellClusterState(
        cells={
            "cell-1": _cell("cell-1", neighbor_cell_ids=("cell-2",)),
            "cell-2": _cell("cell-2", neighbor_cell_ids=("cell-1",)),
        },
        ecm=_ecm(),
        junctions=(
            JunctionState.from_cell_pair(
                junction_id="j-1",
                cell_id_pair=("cell-1", "cell-2"),
                contact_length_um=1.0,
                age_s=0.0,
                maturity=0.1,
            ),
            JunctionState.from_cell_pair(
                junction_id="j-1",
                cell_id_pair=("cell-1", "cell-2"),
                contact_length_um=1.0,
                age_s=0.0,
                maturity=0.1,
            ),
        ),
    )
    with pytest.raises(ValueError, match="duplicate junction_id"):
        cluster.validate()


def test_asymmetric_neighbor_rejected():
    cluster = CellClusterState(
        cells={
            "cell-1": _cell("cell-1", neighbor_cell_ids=("cell-2",)),
            "cell-2": _cell("cell-2"),  # missing back-link
        },
        ecm=_ecm(),
    )
    with pytest.raises(ValueError, match="asymmetric"):
        cluster.validate()


def test_neighbor_to_unknown_cell_rejected():
    cluster = CellClusterState(
        cells={"cell-1": _cell("cell-1", neighbor_cell_ids=("cell-2",))},
        ecm=_ecm(),
    )
    with pytest.raises(ValueError, match="not in cluster"):
        cluster.validate()


def test_symmetric_neighbor_pair_validates_and_edges():
    cluster = CellClusterState(
        cells={
            "cell-1": _cell("cell-1", neighbor_cell_ids=("cell-2",)),
            "cell-2": _cell("cell-2", neighbor_cell_ids=("cell-1",)),
        },
        ecm=_ecm(),
    )
    cluster.validate()
    assert cluster.neighbor_graph_edges() == (("cell-1", "cell-2"),)


def test_parent_id_outside_cluster_requires_external_declaration():
    cluster = CellClusterState(
        cells={
            "cell-1": _cell("cell-1", parent_cell_id="founder"),
        },
        ecm=_ecm(),
    )
    with pytest.raises(ValueError, match="parent_cell_id"):
        cluster.validate()
    cluster_ok = CellClusterState(
        cells={
            "cell-1": _cell("cell-1", parent_cell_id="founder"),
        },
        ecm=_ecm(),
        external_parent_ids=("founder",),
    )
    cluster_ok.validate()


def test_external_parent_overlap_with_cell_rejected():
    cluster = CellClusterState(
        cells={"cell-1": _cell("cell-1")},
        ecm=_ecm(),
        external_parent_ids=("cell-1",),
    )
    with pytest.raises(ValueError, match="overlap"):
        cluster.validate()


def test_role_label_unknown_cell_rejected():
    cluster = CellClusterState(
        cells={"cell-1": _cell("cell-1")},
        ecm=_ecm(),
        role_labels={"cell-7": "leader"},
    )
    with pytest.raises(ValueError, match="role_labels"):
        cluster.validate()


def test_role_invalid_value_rejected():
    cluster = CellClusterState(
        cells={"cell-1": _cell("cell-1")},
        ecm=_ecm(),
        role_labels={"cell-1": "boss"},  # type: ignore[dict-item]
    )
    with pytest.raises(ValueError, match="role"):
        cluster.validate()


def test_role_valid_values_accepted():
    cluster = CellClusterState(
        cells={
            "cell-1": _cell("cell-1", neighbor_cell_ids=("cell-2",)),
            "cell-2": _cell("cell-2", neighbor_cell_ids=("cell-1",)),
        },
        ecm=_ecm(),
        role_labels={"cell-1": "leader", "cell-2": "follower"},
    )
    cluster.validate()
    assert cluster.get_role("cell-1") == "leader"
