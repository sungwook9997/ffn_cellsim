from __future__ import annotations

import os

import numpy as np

from acs.v2.cell_cluster import CellClusterState
from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.measurement_boundary import MeasurementBoundary
from acs.v2.output.frame_dump import write_frame
from acs.v2.single_cell import SingleCellState
from acs.v2.viz.stub3d import render_frame_html, render_frame_png


def _ecm() -> ECMSubstrateState:
    nx, ny = 3, 3
    orientation = np.zeros((nx, ny, 2, 2), dtype=np.float64)
    orientation[..., 0, 0] = 1.0
    orientation[..., 1, 1] = 1.0
    return ECMSubstrateState(
        origin_um_xy=(0.0, 0.0),
        spacing_um=2.0,
        stiffness_kpa=np.linspace(0.0, 1.0, nx * ny).reshape(nx, ny),
        ligand_density=np.full((nx, ny), 0.5, dtype=np.float64),
        fiber_density=np.full((nx, ny), 0.2, dtype=np.float64),
        orientation_tensor=orientation,
        accumulated_traction_nNs_per_um2=np.zeros((nx, ny), dtype=np.float64),
    )


def _cell(cell_id: str, *, height_um: float = 5.0) -> SingleCellState:
    boundary = MeasurementBoundary.from_array(
        np.array([[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0]]),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    return SingleCellState(
        cell_id=cell_id,
        time_s=0.0,
        measurement_boundary=boundary,
        height_um=height_um,
    )


def _frame(tmp_path) -> str:
    cluster = CellClusterState(cells={"cell-1": _cell("cell-1")}, ecm=_ecm())
    out = tmp_path / "frame.h5"
    write_frame(str(out), cluster, time_s=12.0)
    return str(out)


def test_render_png_writes_non_empty_file(tmp_path):
    frame_path = _frame(tmp_path)
    png_path = tmp_path / "frame.png"
    out = render_frame_png(frame_path, str(png_path))
    assert os.path.exists(out)
    assert os.path.getsize(out) > 1024  # PNG should not be tiny placeholder
    with open(out, "rb") as fh:
        head = fh.read(8)
    assert head[:4] == b"\x89PNG"


def test_render_html_writes_self_contained_file(tmp_path):
    frame_path = _frame(tmp_path)
    html_path = tmp_path / "frame.html"
    out = render_frame_html(frame_path, str(html_path))
    assert os.path.exists(out)
    text = open(out, encoding="utf-8").read()
    assert text.startswith("<!doctype html>")
    assert "<svg" in text
    assert "cell-1" in text
    assert "v2 stub3d" in text
    assert "schema_version = 2" in text


def test_render_handles_two_cells(tmp_path):
    cell_a = _cell("cell-a")
    boundary_b = MeasurementBoundary.from_array(
        np.array(
            [
                [10.0, 10.0],
                [12.0, 10.0],
                [12.0, 11.0],
                [10.0, 11.0],
            ]
        ),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    cell_b = SingleCellState(
        cell_id="cell-b",
        time_s=0.0,
        measurement_boundary=boundary_b,
        height_um=8.0,
    )
    cluster = CellClusterState(cells={"cell-a": cell_a, "cell-b": cell_b}, ecm=_ecm())
    frame = tmp_path / "two.h5"
    write_frame(str(frame), cluster, time_s=20.0)

    png = render_frame_png(str(frame), str(tmp_path / "two.png"))
    html = render_frame_html(str(frame), str(tmp_path / "two.html"))

    assert os.path.getsize(png) > 1024
    text = open(html, encoding="utf-8").read()
    assert "cell-a" in text
    assert "cell-b" in text


def test_render_html_escapes_cell_id_special_chars(tmp_path):
    boundary = MeasurementBoundary.from_array(
        np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
        coordinate_convention="world_um_y_up",
        source_modality="manual_outline",
    )
    cell = SingleCellState(
        cell_id="<script>alert(1)</script>",
        time_s=0.0,
        measurement_boundary=boundary,
    )
    cluster = CellClusterState(cells={"<script>alert(1)</script>": cell}, ecm=_ecm())
    frame = tmp_path / "esc.h5"
    write_frame(str(frame), cluster, time_s=0.0)

    html = render_frame_html(str(frame), str(tmp_path / "esc.html"))
    text = open(html, encoding="utf-8").read()
    assert "<script>alert(1)</script>" not in text
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in text


def test_render_png_into_nested_dir(tmp_path):
    frame = _frame(tmp_path)
    out_dir = tmp_path / "nested" / "deep"
    png = out_dir / "frame.png"
    render_frame_png(frame, str(png))
    assert png.exists()
