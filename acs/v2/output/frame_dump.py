"""HDF5 frame dump writer and reader for v2 single-cell / small-cluster state.

The frame schema is defined in ``docs/v2/v2_phase1_plan_consolidated.md`` §7
and §10 ("frame_dump = schema fidelity, round-trip, complete artifact
metadata"). A frame stores one timepoint of a :class:`CellClusterState`
(or a single :class:`SingleCellState` wrapped in a one-cell cluster)
plus the underlying ECM grid and any associated events.

Boundary vertices are stored CSR-style (one flat ``(N, 2)`` array plus
an ``offsets`` array of length ``n_cells + 1``) so cells can carry
variable vertex counts without padding. The same CSR pattern is used
for variable-length per-row sequences inside event tables
(``associated_adhesion_ids`` per protrusion, ``contact_edge_indices_a``
and ``contact_edge_indices_b`` per junction, and per-cell
``neighbor_cell_ids``).

Optional scalar fields use sentinel values:

- ``NaN`` for optional floats (``mechanosignal_yap_taz``, protrusion
  ``end_time_s``, position/direction xy components, widths, force
  candidate, ``lifetime_s``, ``confidence``, junction
  ``cadherin_proxy``, ``tension_proxy_nN``,
  ``contact_inhibition_signal``).
- Empty bytes string for optional strings (``cell_cycle_phase``,
  ``parent_cell_id``, protrusion ``event_id``, FA
  ``linked_protrusion_id``).

Writer atomicity: ``write_frame`` writes to a unique sibling tmp file
created by :func:`tempfile.mkstemp` (named
``<basename>.<random>.tmp`` next to the final path) and atomically
renames it onto the target with :func:`os.replace`. On any failure the
tmp file is removed in a ``finally``-style block, so the target path
either contains a complete frame or remains untouched. The writer
does not overwrite a pre-existing tmp file (each invocation gets its
own unique tmp).

Sanity Gate (frame_dump):
    1. Dimensional analysis: every grid/array carries an explicit unit
       in its dataset name (``boundary_vertices_flat_um_xy``,
       ``area_um2``, ``perimeter_um``, ``stiffness_kpa``,
       ``accumulated_traction_nNs_per_um2``). dt/CFL is N/A — this
       module is I/O.
    2. Boundary cases: zero-cell cluster, missing ECM, missing
       protrusion/adhesion/junction tables, mismatched cell ids,
       invalid CSR offsets, and missing required datasets are guarded
       by the schema check on read.
    3. Conservation/provenance invariants: full round-trip of every
       canonical ``CellClusterState`` schema field, including slow-
       biology hooks, neighbor graph, role labels, lineage references,
       and all optional event/adhesion/junction metadata. The reader
       rebuilds ``MeasurementBoundary`` via ``from_array`` so the
       canonical orientation/coordinate-convention contract holds for
       every restored cell.
    4. Numerical sanity: float64 throughout for geometry. Integer cell
       indices are stored as int64. Optional scalars use NaN; optional
       strings use the empty bytes sentinel.
    5. Sign/sense check: not applicable — this module performs no
       physics, only persistence. The validation contract on the
       restored state is enforced by the schema modules.
    6. Measurement-protocol consistency (Hard Rule 11): the frame
       persists ``boundary_vertices_flat_um_xy`` and the per-cell
       ``area_um2``/``perimeter_um`` derived from the canonical
       :class:`MeasurementBoundary`, which matches PI top-down
       segmentation modality.

Magic-Number Block: this module declares no tunable numeric. The frame
``schema_version`` is an identifier, not a tunable.
"""

from __future__ import annotations

import math
import os
import tempfile
from dataclasses import dataclass
from typing import Optional

import h5py
import numpy as np

from acs.v2.cell_cluster import CellClusterState
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState
from acs.v2.junction import JunctionState
from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.protrusion import ProtrusionEvent
from acs.v2.single_cell import SingleCellState

_FRAME_SCHEMA_VERSION = 2


class FrameDumpError(IOError):
    """Raised when an HDF5 frame is malformed or fails schema check."""


@dataclass(frozen=True)
class FrameDumpReadResult:
    """Materialized cluster state plus the file's frame metadata."""

    cluster: CellClusterState
    time_s: float
    schema_version: int


def _encode_str(values) -> np.ndarray:
    return np.asarray(list(values), dtype="S")


def _decode_str(values) -> tuple[str, ...]:
    return tuple(v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values)


def _opt_float(value: Optional[float]) -> float:
    return float("nan") if value is None else float(value)


def _from_opt_float(value) -> Optional[float]:
    f = float(value)
    return None if math.isnan(f) else f


def _opt_str(value: Optional[str]) -> str:
    return "" if value is None else str(value)


def _from_opt_str(value: str) -> Optional[str]:
    return None if value == "" else value


def _flatten_str_csr(seqs: list[tuple[str, ...]]) -> tuple[np.ndarray, np.ndarray]:
    flat: list[str] = []
    offsets = [0]
    for seq in seqs:
        flat.extend(seq)
        offsets.append(offsets[-1] + len(seq))
    return _encode_str(flat), np.asarray(offsets, dtype=np.int64)


def _flatten_int_csr(seqs: list[tuple[int, ...]]) -> tuple[np.ndarray, np.ndarray]:
    flat: list[int] = []
    offsets = [0]
    for seq in seqs:
        flat.extend(int(v) for v in seq)
        offsets.append(offsets[-1] + len(seq))
    return np.asarray(flat, dtype=np.int64), np.asarray(offsets, dtype=np.int64)


def _slice_csr_str(flat, offsets, idx: int) -> tuple[str, ...]:
    start = int(offsets[idx])
    end = int(offsets[idx + 1])
    return tuple(_decode_str(flat[start:end]))


def _slice_csr_int(flat, offsets, idx: int) -> tuple[int, ...]:
    start = int(offsets[idx])
    end = int(offsets[idx + 1])
    return tuple(int(v) for v in flat[start:end])


def _build_cell_payload(cluster: CellClusterState) -> dict:
    cluster.validate()
    cell_ids = list(cluster.cell_ids)
    cells = [cluster.cells[cid] for cid in cell_ids]

    flat_vertices: list[np.ndarray] = []
    boundary_offsets = [0]
    centroids = np.zeros((len(cells), 2), dtype=np.float64)
    polarities = np.full((len(cells), 2), np.nan, dtype=np.float64)
    heights = np.full((len(cells),), np.nan, dtype=np.float64)
    areas = np.zeros((len(cells),), dtype=np.float64)
    perimeters = np.zeros((len(cells),), dtype=np.float64)
    state_strings: list[str] = []
    times = np.zeros((len(cells),), dtype=np.float64)
    cell_age_s = np.zeros((len(cells),), dtype=np.float64)
    division_count = np.zeros((len(cells),), dtype=np.int64)
    cell_cycle_phase: list[str] = []
    parent_cell_id: list[str] = []
    mechanosignal = np.full((len(cells),), np.nan, dtype=np.float64)
    neighbor_lists: list[tuple[str, ...]] = []

    for idx, cell in enumerate(cells):
        verts = np.asarray(cell.measurement_boundary.vertices_xy_um, dtype=np.float64)
        flat_vertices.append(verts)
        boundary_offsets.append(boundary_offsets[-1] + verts.shape[0])
        centroids[idx] = verts.mean(axis=0)
        if cell.polarity_xy is not None:
            polarities[idx] = np.asarray(cell.polarity_xy, dtype=np.float64)
        if cell.height_um is not None:
            heights[idx] = float(cell.height_um)
        areas[idx] = cell.measurement_boundary.projected_area_um2()
        perimeters[idx] = cell.measurement_boundary.perimeter_um()
        state_strings.append(cell.cell_state)
        times[idx] = float(cell.time_s)
        cell_age_s[idx] = float(cell.cell_age_s)
        division_count[idx] = int(cell.division_count)
        cell_cycle_phase.append(_opt_str(cell.cell_cycle_phase))
        parent_cell_id.append(_opt_str(cell.parent_cell_id))
        mechanosignal[idx] = _opt_float(cell.mechanosignal_yap_taz)
        neighbor_lists.append(cell.neighbor_cell_ids)

    flat = (
        np.concatenate(flat_vertices, axis=0)
        if flat_vertices
        else np.zeros((0, 2), dtype=np.float64)
    )
    neighbor_flat, neighbor_offsets = _flatten_str_csr(neighbor_lists)

    return dict(
        cell_ids=cell_ids,
        cells=cells,
        boundary_vertices_flat_um_xy=flat,
        boundary_offsets=np.asarray(boundary_offsets, dtype=np.int64),
        centroid_um_xy=centroids,
        polarity_xy=polarities,
        height_um=heights,
        area_um2=areas,
        perimeter_um=perimeters,
        cell_state=state_strings,
        time_s=times,
        cell_age_s=cell_age_s,
        division_count=division_count,
        cell_cycle_phase=cell_cycle_phase,
        parent_cell_id=parent_cell_id,
        mechanosignal_yap_taz=mechanosignal,
        neighbor_flat=neighbor_flat,
        neighbor_offsets=neighbor_offsets,
    )


def _write_cells(group: h5py.Group, payload: dict) -> None:
    group.create_dataset("id", data=_encode_str(payload["cell_ids"]))
    group.create_dataset("state", data=_encode_str(payload["cell_state"]))
    group.create_dataset("centroid_um_xy", data=payload["centroid_um_xy"])
    group.create_dataset("polarity_xy", data=payload["polarity_xy"])
    group.create_dataset("height_um", data=payload["height_um"])
    group.create_dataset(
        "boundary_vertices_flat_um_xy",
        data=payload["boundary_vertices_flat_um_xy"],
    )
    group.create_dataset("boundary_offsets", data=payload["boundary_offsets"])
    group.create_dataset("area_um2", data=payload["area_um2"])
    group.create_dataset("perimeter_um", data=payload["perimeter_um"])
    group.create_dataset("time_s", data=payload["time_s"])
    group.create_dataset("cell_age_s", data=payload["cell_age_s"])
    group.create_dataset("division_count", data=payload["division_count"])
    group.create_dataset("cell_cycle_phase", data=_encode_str(payload["cell_cycle_phase"]))
    group.create_dataset("parent_cell_id", data=_encode_str(payload["parent_cell_id"]))
    group.create_dataset("mechanosignal_yap_taz", data=payload["mechanosignal_yap_taz"])
    group.create_dataset("neighbor_cell_ids_flat", data=payload["neighbor_flat"])
    group.create_dataset("neighbor_cell_ids_offsets", data=payload["neighbor_offsets"])


def _write_protrusions(group: h5py.Group, cluster: CellClusterState) -> None:
    rows: list[ProtrusionEvent] = []
    cell_ids: list[str] = []
    associated_lists: list[tuple[str, ...]] = []
    for cell in cluster.cells.values():
        for event in cell.protrusions:
            rows.append(event)
            cell_ids.append(cell.cell_id)
            associated_lists.append(event.associated_adhesion_ids)

    n = len(rows)

    def col_str(values: list[str]) -> np.ndarray:
        return _encode_str(values) if n else np.array([], dtype="S")

    def col_float(values: list[float]) -> np.ndarray:
        return np.asarray(values, dtype=np.float64) if n else np.array([], dtype=np.float64)

    def opt_xy(values: list[Optional[tuple[float, float]]]) -> tuple[np.ndarray, np.ndarray]:
        xs = [_opt_float(v[0]) if v is not None else float("nan") for v in values]
        ys = [_opt_float(v[1]) if v is not None else float("nan") for v in values]
        return col_float(xs), col_float(ys)

    group.create_dataset("cell_id", data=col_str(cell_ids))
    group.create_dataset("event_id", data=col_str([_opt_str(e.event_id) for e in rows]))
    group.create_dataset("event_type", data=col_str([e.event_type for e in rows]))
    group.create_dataset("state", data=col_str([e.state for e in rows]))
    group.create_dataset("source", data=col_str([e.source for e in rows]))
    group.create_dataset("start_time_s", data=col_float([e.start_time_s for e in rows]))
    group.create_dataset(
        "end_time_s", data=col_float([_opt_float(e.end_time_s) for e in rows])
    )
    group.create_dataset(
        "boundary_angle_rad", data=col_float([e.boundary_angle_rad for e in rows])
    )
    group.create_dataset("length_um", data=col_float([e.length_um for e in rows]))
    root_x, root_y = opt_xy([e.root_position_um_xy for e in rows])
    tip_x, tip_y = opt_xy([e.tip_position_um_xy for e in rows])
    dir_x, dir_y = opt_xy([e.direction_um_xy for e in rows])
    group.create_dataset("root_position_um_x", data=root_x)
    group.create_dataset("root_position_um_y", data=root_y)
    group.create_dataset("tip_position_um_x", data=tip_x)
    group.create_dataset("tip_position_um_y", data=tip_y)
    group.create_dataset("direction_um_x", data=dir_x)
    group.create_dataset("direction_um_y", data=dir_y)
    group.create_dataset(
        "width_um", data=col_float([_opt_float(e.width_um) for e in rows])
    )
    group.create_dataset(
        "width_rad", data=col_float([_opt_float(e.width_rad) for e in rows])
    )
    group.create_dataset(
        "force_candidate_nN", data=col_float([_opt_float(e.force_candidate_nN) for e in rows])
    )
    group.create_dataset(
        "lifetime_s", data=col_float([_opt_float(e.lifetime_s) for e in rows])
    )
    group.create_dataset(
        "confidence", data=col_float([_opt_float(e.confidence) for e in rows])
    )
    associated_flat, associated_offsets = _flatten_str_csr(associated_lists)
    group.create_dataset("associated_adhesion_ids_flat", data=associated_flat)
    group.create_dataset("associated_adhesion_ids_offsets", data=associated_offsets)


def _read_protrusions(group: h5py.Group) -> dict[str, list[ProtrusionEvent]]:
    cell_ids = _decode_str(group["cell_id"][...])
    n = len(cell_ids)
    event_ids = _decode_str(group["event_id"][...])
    event_types = _decode_str(group["event_type"][...])
    states = _decode_str(group["state"][...])
    sources = _decode_str(group["source"][...])
    starts = group["start_time_s"][...]
    ends = group["end_time_s"][...]
    angles = group["boundary_angle_rad"][...]
    lengths = group["length_um"][...]
    root_x = group["root_position_um_x"][...]
    root_y = group["root_position_um_y"][...]
    tip_x = group["tip_position_um_x"][...]
    tip_y = group["tip_position_um_y"][...]
    dir_x = group["direction_um_x"][...]
    dir_y = group["direction_um_y"][...]
    widths_um = group["width_um"][...]
    widths_rad = group["width_rad"][...]
    force = group["force_candidate_nN"][...]
    lifetime = group["lifetime_s"][...]
    confidence = group["confidence"][...]
    associated_flat = group["associated_adhesion_ids_flat"][...]
    associated_offsets = group["associated_adhesion_ids_offsets"][...]
    if len(associated_offsets) != n + 1:
        raise FrameDumpError(
            "associated_adhesion_ids offsets length mismatch with protrusion rows"
        )

    def opt_xy(x, y) -> Optional[tuple[float, float]]:
        if math.isnan(float(x)) or math.isnan(float(y)):
            return None
        return (float(x), float(y))

    out: dict[str, list[ProtrusionEvent]] = {}
    for i, cid in enumerate(cell_ids):
        event = ProtrusionEvent(
            cell_id=cid,
            start_time_s=float(starts[i]),
            boundary_angle_rad=float(angles[i]),
            length_um=float(lengths[i]),
            event_id=_from_opt_str(event_ids[i]),
            end_time_s=_from_opt_float(ends[i]),
            event_type=event_types[i],  # type: ignore[arg-type]
            state=states[i],  # type: ignore[arg-type]
            root_position_um_xy=opt_xy(root_x[i], root_y[i]),
            tip_position_um_xy=opt_xy(tip_x[i], tip_y[i]),
            direction_um_xy=opt_xy(dir_x[i], dir_y[i]),
            width_um=_from_opt_float(widths_um[i]),
            width_rad=_from_opt_float(widths_rad[i]),
            force_candidate_nN=_from_opt_float(force[i]),
            associated_adhesion_ids=_slice_csr_str(associated_flat, associated_offsets, i),
            source=sources[i],  # type: ignore[arg-type]
            lifetime_s=_from_opt_float(lifetime[i]),
            confidence=_from_opt_float(confidence[i]),
        )
        out.setdefault(cid, []).append(event)
    return out


def _write_focal_adhesions(group: h5py.Group, cluster: CellClusterState) -> None:
    rows: list[FocalAdhesionState] = []
    cell_ids: list[str] = []
    for cell in cluster.cells.values():
        for fa in cell.adhesions:
            rows.append(fa)
            cell_ids.append(cell.cell_id)
    n = len(rows)

    def col_str(values: list[str]) -> np.ndarray:
        return _encode_str(values) if n else np.array([], dtype="S")

    def col_float(values: list[float]) -> np.ndarray:
        return np.asarray(values, dtype=np.float64) if n else np.array([], dtype=np.float64)

    group.create_dataset("adhesion_id", data=col_str([fa.adhesion_id for fa in rows]))
    group.create_dataset("cell_id", data=col_str(cell_ids))
    group.create_dataset("state", data=col_str([fa.state for fa in rows]))
    group.create_dataset("source", data=col_str([fa.source for fa in rows]))
    group.create_dataset(
        "linked_protrusion_id",
        data=col_str([_opt_str(fa.linked_protrusion_id) for fa in rows]),
    )
    group.create_dataset(
        "position_um_x", data=col_float([fa.position_um_xy[0] for fa in rows])
    )
    group.create_dataset(
        "position_um_y", data=col_float([fa.position_um_xy[1] for fa in rows])
    )
    group.create_dataset("age_s", data=col_float([fa.age_s for fa in rows]))
    group.create_dataset("maturity", data=col_float([fa.maturity for fa in rows]))
    group.create_dataset(
        "bound_fraction", data=col_float([fa.bound_fraction for fa in rows])
    )
    group.create_dataset(
        "traction_force_nN_x",
        data=col_float([fa.traction_force_nN_xy[0] for fa in rows]),
    )
    group.create_dataset(
        "traction_force_nN_y",
        data=col_float([fa.traction_force_nN_xy[1] for fa in rows]),
    )


def _read_focal_adhesions(group: h5py.Group) -> dict[str, list[FocalAdhesionState]]:
    adhesion_ids = _decode_str(group["adhesion_id"][...])
    cell_ids = _decode_str(group["cell_id"][...])
    states = _decode_str(group["state"][...])
    sources = _decode_str(group["source"][...])
    linked = _decode_str(group["linked_protrusion_id"][...])
    pos_x = group["position_um_x"][...]
    pos_y = group["position_um_y"][...]
    age_s = group["age_s"][...]
    maturity = group["maturity"][...]
    bound_fraction = group["bound_fraction"][...]
    traction_x = group["traction_force_nN_x"][...]
    traction_y = group["traction_force_nN_y"][...]
    out: dict[str, list[FocalAdhesionState]] = {}
    for i, cid in enumerate(cell_ids):
        fa = FocalAdhesionState(
            adhesion_id=adhesion_ids[i],
            cell_id=cid,
            position_um_xy=(float(pos_x[i]), float(pos_y[i])),
            age_s=float(age_s[i]),
            maturity=float(maturity[i]),
            bound_fraction=float(bound_fraction[i]),
            state=states[i],  # type: ignore[arg-type]
            traction_force_nN_xy=(float(traction_x[i]), float(traction_y[i])),
            linked_protrusion_id=_from_opt_str(linked[i]),
            source=sources[i],  # type: ignore[arg-type]
        )
        out.setdefault(cid, []).append(fa)
    return out


def _write_junctions(group: h5py.Group, cluster: CellClusterState) -> None:
    junctions = list(cluster.junctions)
    n = len(junctions)

    def col_str(values: list[str]) -> np.ndarray:
        return _encode_str(values) if n else np.array([], dtype="S")

    def col_float(values: list[float]) -> np.ndarray:
        return np.asarray(values, dtype=np.float64) if n else np.array([], dtype=np.float64)

    group.create_dataset(
        "junction_id", data=col_str([j.junction_id for j in junctions])
    )
    group.create_dataset(
        "cell_id_a", data=col_str([j.cell_id_a for j in junctions])
    )
    group.create_dataset(
        "cell_id_b", data=col_str([j.cell_id_b for j in junctions])
    )
    group.create_dataset("state", data=col_str([j.state for j in junctions]))
    group.create_dataset(
        "contact_length_um",
        data=col_float([j.contact_length_um for j in junctions]),
    )
    group.create_dataset("age_s", data=col_float([j.age_s for j in junctions]))
    group.create_dataset(
        "maturity", data=col_float([j.maturity for j in junctions])
    )
    group.create_dataset(
        "cadherin_proxy",
        data=col_float([_opt_float(j.cadherin_proxy) for j in junctions]),
    )
    group.create_dataset(
        "tension_proxy_nN",
        data=col_float([_opt_float(j.tension_proxy_nN) for j in junctions]),
    )
    group.create_dataset(
        "contact_inhibition_signal",
        data=col_float([_opt_float(j.contact_inhibition_signal) for j in junctions]),
    )
    edges_a_flat, edges_a_off = _flatten_int_csr(
        [j.contact_edge_indices_a for j in junctions]
    )
    edges_b_flat, edges_b_off = _flatten_int_csr(
        [j.contact_edge_indices_b for j in junctions]
    )
    group.create_dataset("contact_edge_indices_a_flat", data=edges_a_flat)
    group.create_dataset("contact_edge_indices_a_offsets", data=edges_a_off)
    group.create_dataset("contact_edge_indices_b_flat", data=edges_b_flat)
    group.create_dataset("contact_edge_indices_b_offsets", data=edges_b_off)


def _read_junctions(group: h5py.Group) -> tuple[JunctionState, ...]:
    junction_ids = _decode_str(group["junction_id"][...])
    cell_ids_a = _decode_str(group["cell_id_a"][...])
    cell_ids_b = _decode_str(group["cell_id_b"][...])
    states = _decode_str(group["state"][...])
    contact_length = group["contact_length_um"][...]
    age = group["age_s"][...]
    maturity = group["maturity"][...]
    cadherin = group["cadherin_proxy"][...]
    tension = group["tension_proxy_nN"][...]
    inhibition = group["contact_inhibition_signal"][...]
    edges_a_flat = group["contact_edge_indices_a_flat"][...]
    edges_a_off = group["contact_edge_indices_a_offsets"][...]
    edges_b_flat = group["contact_edge_indices_b_flat"][...]
    edges_b_off = group["contact_edge_indices_b_offsets"][...]

    n = len(junction_ids)
    if len(edges_a_off) != n + 1 or len(edges_b_off) != n + 1:
        raise FrameDumpError("junction contact edge offsets length mismatch with junction rows")

    out: list[JunctionState] = []
    for i, jid in enumerate(junction_ids):
        out.append(
            JunctionState(
                junction_id=jid,
                cell_id_a=cell_ids_a[i],
                cell_id_b=cell_ids_b[i],
                contact_length_um=float(contact_length[i]),
                age_s=float(age[i]),
                maturity=float(maturity[i]),
                state=states[i],  # type: ignore[arg-type]
                contact_edge_indices_a=_slice_csr_int(edges_a_flat, edges_a_off, i),
                contact_edge_indices_b=_slice_csr_int(edges_b_flat, edges_b_off, i),
                cadherin_proxy=_from_opt_float(cadherin[i]),
                tension_proxy_nN=_from_opt_float(tension[i]),
                contact_inhibition_signal=_from_opt_float(inhibition[i]),
            )
        )
    return tuple(out)


def _write_ecm(group: h5py.Group, ecm: ECMSubstrateState) -> None:
    group.attrs["origin_um_x"] = ecm.origin_um_xy[0]
    group.attrs["origin_um_y"] = ecm.origin_um_xy[1]
    group.attrs["spacing_um"] = ecm.spacing_um
    group.attrs["source"] = ecm.source
    group.create_dataset("stiffness_kpa", data=np.asarray(ecm.stiffness_kpa, dtype=np.float64))
    group.create_dataset(
        "ligand_density", data=np.asarray(ecm.ligand_density, dtype=np.float64)
    )
    group.create_dataset(
        "fiber_density", data=np.asarray(ecm.fiber_density, dtype=np.float64)
    )
    group.create_dataset(
        "orientation_tensor",
        data=np.asarray(ecm.orientation_tensor, dtype=np.float64),
    )
    group.create_dataset(
        "accumulated_traction_nNs_per_um2",
        data=np.asarray(ecm.accumulated_traction_nNs_per_um2, dtype=np.float64),
    )


def _read_ecm(group: h5py.Group) -> ECMSubstrateState:
    return ECMSubstrateState(
        origin_um_xy=(float(group.attrs["origin_um_x"]), float(group.attrs["origin_um_y"])),
        spacing_um=float(group.attrs["spacing_um"]),
        stiffness_kpa=np.asarray(group["stiffness_kpa"][...], dtype=np.float64),
        ligand_density=np.asarray(group["ligand_density"][...], dtype=np.float64),
        fiber_density=np.asarray(group["fiber_density"][...], dtype=np.float64),
        orientation_tensor=np.asarray(group["orientation_tensor"][...], dtype=np.float64),
        accumulated_traction_nNs_per_um2=np.asarray(
            group["accumulated_traction_nNs_per_um2"][...], dtype=np.float64
        ),
        source=str(group.attrs["source"]),  # type: ignore[arg-type]
    )


def _write_cluster_metadata(group: h5py.Group, cluster: CellClusterState) -> None:
    role_keys = sorted(cluster.role_labels.keys())
    role_values = [cluster.role_labels[k] for k in role_keys]
    group.create_dataset("role_label_cell_ids", data=_encode_str(role_keys))
    group.create_dataset("role_label_values", data=_encode_str(role_values))
    group.create_dataset(
        "external_parent_ids", data=_encode_str(cluster.external_parent_ids)
    )


def _read_cluster_metadata(group: h5py.Group) -> dict:
    keys = _decode_str(group["role_label_cell_ids"][...])
    values = _decode_str(group["role_label_values"][...])
    if len(keys) != len(values):
        raise FrameDumpError("role_label_cell_ids and role_label_values length mismatch")
    role_labels = {k: v for k, v in zip(keys, values)}
    external = _decode_str(group["external_parent_ids"][...])
    return {"role_labels": role_labels, "external_parent_ids": tuple(external)}


def write_frame(path: str, cluster: CellClusterState, time_s: float) -> str:
    """Write ``cluster`` at ``time_s`` to ``path`` atomically.

    Returns the absolute path written. Raises :class:`FrameDumpError`
    on schema or I/O failure.
    """

    cluster.validate()
    if not np.isfinite(time_s) or time_s < 0.0:
        raise FrameDumpError("time_s must be finite and non-negative")
    abs_path = os.path.abspath(path)
    directory = os.path.dirname(abs_path) or "."
    os.makedirs(directory, exist_ok=True)

    payload = _build_cell_payload(cluster)

    fd, tmp_path = tempfile.mkstemp(
        prefix=os.path.basename(abs_path) + ".",
        suffix=".tmp",
        dir=directory,
    )
    os.close(fd)
    try:
        with h5py.File(tmp_path, "w") as f:
            f.attrs["schema_version"] = _FRAME_SCHEMA_VERSION
            meta = f.create_group("meta")
            meta.attrs["time_s"] = float(time_s)

            _write_cells(f.create_group("cells"), payload)
            _write_protrusions(f.create_group("protrusions"), cluster)
            _write_focal_adhesions(f.create_group("focal_adhesions"), cluster)
            _write_junctions(f.create_group("junctions"), cluster)
            _write_ecm(f.create_group("ecm"), cluster.ecm)
            _write_cluster_metadata(f.create_group("cluster_metadata"), cluster)
            f.create_group("closure")  # P0: schema-only placeholder
        os.replace(tmp_path, abs_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise
    return abs_path


def read_frame(
    path: str,
    *,
    coordinate_convention: str = "world_um_y_up",
    source_modality: str = "frame_dump",
) -> FrameDumpReadResult:
    """Read a frame and reconstruct its :class:`CellClusterState`."""

    abs_path = os.path.abspath(path)
    if not os.path.exists(abs_path):
        raise FrameDumpError(f"frame file not found: {abs_path!r}")
    with h5py.File(abs_path, "r") as f:
        if "schema_version" not in f.attrs:
            raise FrameDumpError("missing schema_version attribute")
        schema_version = int(f.attrs["schema_version"])
        if "meta" not in f or "time_s" not in f["meta"].attrs:
            raise FrameDumpError("missing /meta/time_s")
        time_s = float(f["meta"].attrs["time_s"])

        if "cells" not in f:
            raise FrameDumpError("missing /cells group")
        cells_group = f["cells"]
        required = {
            "id",
            "state",
            "boundary_vertices_flat_um_xy",
            "boundary_offsets",
            "polarity_xy",
            "height_um",
            "time_s",
            "cell_age_s",
            "division_count",
            "cell_cycle_phase",
            "parent_cell_id",
            "mechanosignal_yap_taz",
            "neighbor_cell_ids_flat",
            "neighbor_cell_ids_offsets",
        }
        missing = required - set(cells_group.keys())
        if missing:
            raise FrameDumpError(f"cells group missing datasets: {sorted(missing)!r}")

        cell_ids = _decode_str(cells_group["id"][...])
        states = _decode_str(cells_group["state"][...])
        flat = cells_group["boundary_vertices_flat_um_xy"][...]
        offsets = cells_group["boundary_offsets"][...]
        polarities = cells_group["polarity_xy"][...]
        heights = cells_group["height_um"][...]
        times = cells_group["time_s"][...]
        ages = cells_group["cell_age_s"][...]
        divisions = cells_group["division_count"][...]
        cycle_phases = _decode_str(cells_group["cell_cycle_phase"][...])
        parent_ids = _decode_str(cells_group["parent_cell_id"][...])
        mechanosignal = cells_group["mechanosignal_yap_taz"][...]
        neighbor_flat = cells_group["neighbor_cell_ids_flat"][...]
        neighbor_offsets = cells_group["neighbor_cell_ids_offsets"][...]

        if len(offsets) != len(cell_ids) + 1:
            raise FrameDumpError(
                f"boundary_offsets must have length n_cells+1, got {len(offsets)} "
                f"for {len(cell_ids)} cells"
            )
        if int(offsets[0]) != 0 or int(offsets[-1]) != flat.shape[0]:
            raise FrameDumpError(
                "boundary_offsets must start at 0 and end at len(boundary_vertices)"
            )
        if len(neighbor_offsets) != len(cell_ids) + 1:
            raise FrameDumpError(
                "neighbor_cell_ids_offsets length must be n_cells+1"
            )

        protrusions_by_cell = _read_protrusions(f["protrusions"])
        adhesions_by_cell = _read_focal_adhesions(f["focal_adhesions"])
        junctions = _read_junctions(f["junctions"])
        ecm = _read_ecm(f["ecm"])
        cluster_metadata = _read_cluster_metadata(f["cluster_metadata"])

    cells: dict[str, SingleCellState] = {}
    for idx, cid in enumerate(cell_ids):
        start = int(offsets[idx])
        end = int(offsets[idx + 1])
        verts = np.asarray(flat[start:end], dtype=np.float64)
        boundary = MeasurementBoundary.from_array(
            verts,
            coordinate_convention=coordinate_convention,  # type: ignore[arg-type]
            source_modality=source_modality,
        )
        polarity_row = polarities[idx]
        polarity_value: Optional[tuple[float, float]] = None
        if np.all(np.isfinite(polarity_row)):
            polarity_value = (float(polarity_row[0]), float(polarity_row[1]))
        height_value = float(heights[idx]) if np.isfinite(heights[idx]) else None
        cell = SingleCellState(
            cell_id=cid,
            time_s=float(times[idx]),
            measurement_boundary=boundary,
            height_um=height_value,
            polarity_xy=polarity_value,
            protrusions=list(protrusions_by_cell.get(cid, ())),
            adhesions=list(adhesions_by_cell.get(cid, ())),
            cell_state=states[idx],  # type: ignore[arg-type]
            cell_age_s=float(ages[idx]),
            cell_cycle_phase=_from_opt_str(cycle_phases[idx]),  # type: ignore[arg-type]
            division_count=int(divisions[idx]),
            parent_cell_id=_from_opt_str(parent_ids[idx]),
            mechanosignal_yap_taz=_from_opt_float(mechanosignal[idx]),
            neighbor_cell_ids=_slice_csr_str(neighbor_flat, neighbor_offsets, idx),
        )
        cells[cid] = cell

    cluster = CellClusterState(
        cells=cells,
        ecm=ecm,
        junctions=junctions,
        role_labels=cluster_metadata["role_labels"],
        external_parent_ids=cluster_metadata["external_parent_ids"],
    )
    cluster.validate()
    return FrameDumpReadResult(
        cluster=cluster, time_s=time_s, schema_version=schema_version
    )
