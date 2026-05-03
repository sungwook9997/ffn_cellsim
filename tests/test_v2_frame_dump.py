from __future__ import annotations

import os

import h5py
import numpy as np
import pytest

from acs.v2.cell_cluster import CellClusterState
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.junction import JunctionState
from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.output.frame_dump import FrameDumpError, read_frame, write_frame
from acs.v2.protrusion import ProtrusionEvent
from acs.v2.single_cell import SingleCellState


def _ecm(nx: int = 2, ny: int = 2) -> ECMSubstrateState:
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    return ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=1.0,
        stiffness_kpa=np.full((nx, ny), 0.5, dtype=np.float64),
        ligand_density=np.full((nx, ny), 0.3, dtype=np.float64),
        fiber_density=np.full((nx, ny), 0.1, dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
    )


def _square_cell(
    cell_id: str,
    side_um: float = 2.0,
    *,
    polarity: tuple[float, float] | None = (1.0, 0.0),
    height_um: float | None = 5.0,
    neighbor_cell_ids: tuple[str, ...] = (),
) -> SingleCellState:
    boundary = MeasurementBoundary.from_array(
        np.array(
            [
                [0.0, 0.0],
                [side_um, 0.0],
                [side_um, side_um],
                [0.0, side_um],
            ]
        ),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    return SingleCellState(
        cell_id=cell_id,
        time_s=0.0,
        measurement_boundary=boundary,
        polarity_xy=polarity,
        height_um=height_um,
        neighbor_cell_ids=neighbor_cell_ids,
    )


def _hex_cell(cell_id: str, neighbor_cell_ids: tuple[str, ...] = ()) -> SingleCellState:
    angles = np.linspace(0.0, 2.0 * np.pi, 7)[:-1]
    radius = 1.5
    verts = np.column_stack((10.0 + radius * np.cos(angles), 10.0 + radius * np.sin(angles)))
    boundary = MeasurementBoundary.from_array(
        verts,
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    return SingleCellState(
        cell_id=cell_id,
        time_s=0.0,
        measurement_boundary=boundary,
        neighbor_cell_ids=neighbor_cell_ids,
    )


def test_round_trip_single_cell(tmp_path):
    cluster = CellClusterState(
        cells={"cell-1": _square_cell("cell-1")},
        ecm=_ecm(),
    )
    out = tmp_path / "frame_000000.h5"
    write_frame(str(out), cluster, time_s=120.0)
    result = read_frame(str(out))

    assert result.time_s == pytest.approx(120.0)
    assert result.schema_version == 2
    assert tuple(result.cluster.cells.keys()) == ("cell-1",)
    cell = result.cluster.cells["cell-1"]
    assert cell.measurement_boundary.projected_area_um2() == pytest.approx(4.0)
    assert cell.polarity_xy == pytest.approx((1.0, 0.0))
    assert cell.height_um == pytest.approx(5.0)


def test_variable_vertex_counts_round_trip(tmp_path):
    cells = {
        "square": _square_cell("square", side_um=3.0, neighbor_cell_ids=("hex",)),
        "hex": _hex_cell("hex", neighbor_cell_ids=("square",)),
    }
    cluster = CellClusterState(cells=cells, ecm=_ecm())
    out = tmp_path / "frame.h5"
    write_frame(str(out), cluster, time_s=60.0)
    result = read_frame(str(out))

    restored = result.cluster.cells
    assert restored["square"].measurement_boundary.vertices_xy_um.shape[0] == 4
    assert restored["hex"].measurement_boundary.vertices_xy_um.shape[0] == 6
    assert restored["square"].measurement_boundary.projected_area_um2() == pytest.approx(9.0)


def test_round_trip_preserves_protrusions_and_adhesions(tmp_path):
    cell = _square_cell("cell-1")
    cell = SingleCellState(
        cell_id="cell-1",
        time_s=cell.time_s,
        measurement_boundary=cell.measurement_boundary,
        polarity_xy=cell.polarity_xy,
        height_um=cell.height_um,
        protrusions=[
            ProtrusionEvent(
                cell_id="cell-1",
                start_time_s=0.0,
                boundary_angle_rad=0.5,
                length_um=1.2,
                event_id="evt-1",
                event_type="filopodium",
                state="growing",
            )
        ],
        adhesions=[
            FocalAdhesionState(
                adhesion_id="fa-1",
                cell_id="cell-1",
                position_um_xy=(0.5, 0.5),
                age_s=10.0,
                maturity=0.6,
                bound_fraction=0.7,
                state="mature",
                traction_force_nN_xy=(0.1, -0.2),
            )
        ],
    )
    cluster = CellClusterState(cells={"cell-1": cell}, ecm=_ecm())
    out = tmp_path / "frame.h5"
    write_frame(str(out), cluster, time_s=42.0)
    result = read_frame(str(out))

    restored_cell = result.cluster.cells["cell-1"]
    assert len(restored_cell.protrusions) == 1
    assert restored_cell.protrusions[0].event_id == "evt-1"
    assert restored_cell.protrusions[0].event_type == "filopodium"
    assert len(restored_cell.adhesions) == 1
    assert restored_cell.adhesions[0].state == "mature"
    assert restored_cell.adhesions[0].traction_force_nN_xy == pytest.approx((0.1, -0.2))


def test_round_trip_preserves_junctions(tmp_path):
    cells = {
        "alpha": _square_cell("alpha", neighbor_cell_ids=("beta",)),
        "beta": _square_cell("beta", neighbor_cell_ids=("alpha",)),
    }
    junction = JunctionState.from_cell_pair(
        junction_id="j-1",
        cell_id_pair=("beta", "alpha"),
        contact_length_um=1.0,
        age_s=2.0,
        maturity=0.4,
    )
    cluster = CellClusterState(cells=cells, ecm=_ecm(), junctions=(junction,))
    out = tmp_path / "frame.h5"
    write_frame(str(out), cluster, time_s=10.0)
    result = read_frame(str(out))

    assert len(result.cluster.junctions) == 1
    j = result.cluster.junctions[0]
    assert j.cell_id_a == "alpha"
    assert j.cell_id_b == "beta"
    assert j.contact_length_um == pytest.approx(1.0)


def test_round_trip_preserves_ecm_grid(tmp_path):
    cluster = CellClusterState(cells={"c": _square_cell("c")}, ecm=_ecm(nx=3, ny=4))
    out = tmp_path / "frame.h5"
    write_frame(str(out), cluster, time_s=0.0)
    result = read_frame(str(out))

    ecm = result.cluster.ecm
    assert ecm.grid_shape == (3, 4)
    np.testing.assert_allclose(ecm.stiffness_kpa, cluster.ecm.stiffness_kpa)
    np.testing.assert_allclose(ecm.orientation_tensor, cluster.ecm.orientation_tensor)


def test_invalid_time_rejected(tmp_path):
    cluster = CellClusterState(cells={"c": _square_cell("c")}, ecm=_ecm())
    out = tmp_path / "frame.h5"
    with pytest.raises(FrameDumpError, match="time_s"):
        write_frame(str(out), cluster, time_s=float("nan"))


def test_atomic_write_no_partial_target_on_failure(tmp_path, monkeypatch):
    cluster = CellClusterState(cells={"c": _square_cell("c")}, ecm=_ecm())
    out = tmp_path / "frame.h5"

    def boom(src, dst):
        raise RuntimeError("simulated crash")

    monkeypatch.setattr(os, "replace", boom)
    with pytest.raises(RuntimeError, match="simulated crash"):
        write_frame(str(out), cluster, time_s=0.0)
    assert not out.exists(), "target frame file must not exist after failed write"
    leftover_tmps = [p for p in tmp_path.iterdir() if p.name.endswith(".tmp")]
    assert leftover_tmps == [], f"writer must clean up its tmp files, got {leftover_tmps!r}"


def test_missing_schema_version_rejected(tmp_path):
    bogus = tmp_path / "bogus.h5"
    with h5py.File(bogus, "w") as f:
        f.create_group("cells")
    with pytest.raises(FrameDumpError, match="schema_version"):
        read_frame(str(bogus))


def test_missing_cells_group_rejected(tmp_path):
    bogus = tmp_path / "bogus.h5"
    with h5py.File(bogus, "w") as f:
        f.attrs["schema_version"] = 1
        meta = f.create_group("meta")
        meta.attrs["time_s"] = 0.0
    with pytest.raises(FrameDumpError, match="/cells"):
        read_frame(str(bogus))


def test_offsets_inconsistent_with_flat_rejected(tmp_path):
    cluster = CellClusterState(cells={"c": _square_cell("c")}, ecm=_ecm())
    out = tmp_path / "frame.h5"
    write_frame(str(out), cluster, time_s=0.0)
    with h5py.File(out, "r+") as f:
        del f["cells/boundary_offsets"]
        f["cells"].create_dataset(
            "boundary_offsets", data=np.array([0, 99], dtype=np.int64)
        )
    with pytest.raises(FrameDumpError, match="boundary_offsets"):
        read_frame(str(out))


def test_read_missing_file_rejected(tmp_path):
    with pytest.raises(FrameDumpError, match="not found"):
        read_frame(str(tmp_path / "missing.h5"))


def _full_optional_cell() -> SingleCellState:
    boundary = MeasurementBoundary.from_array(
        np.array([[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0]]),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    return SingleCellState(
        cell_id="cell-A",
        time_s=42.0,
        measurement_boundary=boundary,
        height_um=4.0,
        polarity_xy=(1.0, 2.0),
        cell_state="alive",
        cell_age_s=300.0,
        cell_cycle_phase="S",
        division_count=2,
        parent_cell_id="cell-A-parent",
        mechanosignal_yap_taz=0.7,
        neighbor_cell_ids=("cell-B", "cell-C"),
        protrusions=[
            ProtrusionEvent(
                cell_id="cell-A",
                start_time_s=10.0,
                boundary_angle_rad=0.5,
                length_um=2.5,
                event_id="evt-7",
                end_time_s=12.5,
                event_type="filopodium",
                state="ended",
                root_position_um_xy=(0.5, 0.5),
                tip_position_um_xy=(2.0, 1.0),
                direction_um_xy=(0.6, 0.8),
                width_um=0.3,
                width_rad=0.2,
                force_candidate_nN=4.5,
                associated_adhesion_ids=("fa-A", "fa-B"),
                source="detected",
                lifetime_s=2.5,
                confidence=0.9,
            )
        ],
        adhesions=[
            FocalAdhesionState(
                adhesion_id="fa-A",
                cell_id="cell-A",
                position_um_xy=(0.7, 0.4),
                age_s=20.0,
                maturity=0.55,
                bound_fraction=0.8,
                state="mature",
                traction_force_nN_xy=(0.3, -0.2),
                linked_protrusion_id="evt-7",
                source="marker_detected",
            )
        ],
    )


def _other_cell() -> SingleCellState:
    boundary = MeasurementBoundary.from_array(
        np.array([[10.0, 10.0], [12.0, 10.0], [12.0, 11.5], [10.0, 11.5]]),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    return SingleCellState(
        cell_id="cell-B",
        time_s=42.0,
        measurement_boundary=boundary,
        cell_state="alive",
        neighbor_cell_ids=("cell-A",),
    )


def _other_cell_c() -> SingleCellState:
    boundary = MeasurementBoundary.from_array(
        np.array([[20.0, 20.0], [22.0, 20.0], [22.0, 21.0], [20.0, 21.0]]),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    return SingleCellState(
        cell_id="cell-C",
        time_s=42.0,
        measurement_boundary=boundary,
        cell_state="alive",
        neighbor_cell_ids=("cell-A",),
    )


def test_round_trip_full_optional_metadata(tmp_path):
    cluster = CellClusterState(
        cells={
            "cell-A": _full_optional_cell(),
            "cell-B": _other_cell(),
            "cell-C": _other_cell_c(),
        },
        ecm=_ecm(),
        junctions=(
            JunctionState.from_cell_pair(
                junction_id="j-AB",
                cell_id_pair=("cell-A", "cell-B"),
                contact_length_um=2.0,
                age_s=5.0,
                maturity=0.6,
                state="cortex_coupled",
                contact_edge_indices_a=(1, 2, 3),
                contact_edge_indices_b=(0, 1),
                cadherin_proxy=0.7,
                tension_proxy_nN=0.4,
                contact_inhibition_signal=0.2,
            ),
        ),
        role_labels={"cell-A": "leader", "cell-B": "follower"},
        external_parent_ids=("cell-A-parent",),
    )

    out = tmp_path / "frame.h5"
    write_frame(str(out), cluster, time_s=42.0)
    result = read_frame(str(out))
    restored = result.cluster

    cell_a = restored.cells["cell-A"]
    assert cell_a.cell_cycle_phase == "S"
    assert cell_a.parent_cell_id == "cell-A-parent"
    assert cell_a.mechanosignal_yap_taz == pytest.approx(0.7)
    assert cell_a.neighbor_cell_ids == ("cell-B", "cell-C")
    assert cell_a.cell_age_s == pytest.approx(300.0)
    assert cell_a.division_count == 2

    event = cell_a.protrusions[0]
    assert event.event_id == "evt-7"
    assert event.event_type == "filopodium"
    assert event.state == "ended"
    assert event.end_time_s == pytest.approx(12.5)
    assert event.root_position_um_xy == pytest.approx((0.5, 0.5))
    assert event.tip_position_um_xy == pytest.approx((2.0, 1.0))
    assert event.direction_um_xy == pytest.approx((0.6, 0.8))
    assert event.width_um == pytest.approx(0.3)
    assert event.width_rad == pytest.approx(0.2)
    assert event.force_candidate_nN == pytest.approx(4.5)
    assert event.associated_adhesion_ids == ("fa-A", "fa-B")
    assert event.source == "detected"
    assert event.lifetime_s == pytest.approx(2.5)
    assert event.confidence == pytest.approx(0.9)

    fa = cell_a.adhesions[0]
    assert fa.linked_protrusion_id == "evt-7"
    assert fa.source == "marker_detected"
    assert fa.state == "mature"

    assert restored.role_labels == {"cell-A": "leader", "cell-B": "follower"}
    assert restored.external_parent_ids == ("cell-A-parent",)

    assert len(restored.junctions) == 1
    j = restored.junctions[0]
    assert j.contact_edge_indices_a == (1, 2, 3)
    assert j.contact_edge_indices_b == (0, 1)
    assert j.cadherin_proxy == pytest.approx(0.7)
    assert j.tension_proxy_nN == pytest.approx(0.4)
    assert j.contact_inhibition_signal == pytest.approx(0.2)
    assert j.state == "cortex_coupled"


def test_round_trip_optional_none_fields_stay_none(tmp_path):
    cluster = CellClusterState(
        cells={"cell-1": _square_cell("cell-1", polarity=None, height_um=None)},
        ecm=_ecm(),
    )
    out = tmp_path / "frame.h5"
    write_frame(str(out), cluster, time_s=0.0)
    restored = read_frame(str(out)).cluster
    cell = restored.cells["cell-1"]
    assert cell.polarity_xy is None
    assert cell.height_um is None
    assert cell.cell_cycle_phase is None
    assert cell.parent_cell_id is None
    assert cell.mechanosignal_yap_taz is None
    assert cell.neighbor_cell_ids == ()
