"""Stub 3D visualization for v2 frame dumps.

The renderers here are dev-grade and explicitly **not** the confocal
production path. They read an HDF5 frame written by
:mod:`acs.v2.output.frame_dump` and produce one of two outputs:

- ``render_frame_png(input_path, output_path)``: a top-down 2D view
  with z-axis indication via per-cell height tinting. Intended for
  fast inspection while iterating on schema/data plumbing.
- ``render_frame_html(input_path, output_path)``: a self-contained
  HTML file with an embedded SVG drawing of the same scene plus a
  short caption table. SSH-friendly: no GPU, no JavaScript runtime.

Both renderers consume the canonical :class:`CellClusterState` round-
tripped through ``read_frame`` so they cannot drift from the schema
(per Plan §9 "visualization consumes frame dumps; it must not drive
physics" and Hard Rule 11 "measurement panels separate from 3D
scene"). The scene shows only what the schema declares; quantitative
metric overlays (e.g. PI A/A0) are not embedded in the 3D scene.

Sanity Gate scope: rendering wrapper, no physics. Full 6-item gate
N/A. Boundary-case checks owned here:

- file-not-found / unreadable HDF5 → :class:`FrameDumpError`
  bubbles up from :func:`read_frame` before any plotting work.
- malformed cluster / empty-cell frames are rejected by
  :meth:`CellClusterState.validate` inside ``read_frame`` so the
  renderer never sees a structurally invalid scene.
- ECM grid plotted as substrate background only when the grid has
  non-zero size.

Magic-Number Block: ``_DEFAULT_FIGSIZE_INCHES = (4.0, 4.0)`` and
``_HEIGHT_TINT_SCALE_UM = 10.0`` are presentation defaults documented
as numerical render parameters (no physics tuning, no fitted
constant). They are caller-overridable via function arguments.
"""

from __future__ import annotations

import html as _html
import os
from typing import Optional

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from acs.v2.cell_cluster import CellClusterState  # noqa: E402
from acs.v2.output.frame_dump import read_frame  # noqa: E402

# Presentation defaults; overridable via function arguments.
_DEFAULT_FIGSIZE_INCHES = (4.0, 4.0)
_HEIGHT_TINT_SCALE_UM = 10.0  # used to map height_um -> face alpha [0, 1]
_CELL_FILL_COLOR = "#6699cc"
_CELL_EDGE_COLOR = "#0d3d66"


def _render_axes(
    ax,
    cluster: CellClusterState,
    time_s: float,
    *,
    height_tint_scale_um: float,
) -> None:
    if cluster.ecm.grid_shape[0] > 0 and cluster.ecm.grid_shape[1] > 0:
        nx, ny = cluster.ecm.grid_shape
        spacing = cluster.ecm.spacing_um
        ox, oy = cluster.ecm.origin_um_xy
        extent = (
            ox,
            ox + nx * spacing,
            oy,
            oy + ny * spacing,
        )
        ax.imshow(
            np.asarray(cluster.ecm.stiffness_kpa).T,
            origin="lower",
            extent=extent,
            cmap="Greys",
            alpha=0.4,
            interpolation="nearest",
        )

    for cell_id, cell in cluster.cells.items():
        verts = np.asarray(cell.measurement_boundary.vertices_xy_um, dtype=np.float64)
        if cell.height_um is not None and height_tint_scale_um > 0.0:
            alpha = float(
                np.clip(cell.height_um / height_tint_scale_um, 0.15, 0.95)
            )
        else:
            alpha = 0.55
        ax.fill(
            np.append(verts[:, 0], verts[0, 0]),
            np.append(verts[:, 1], verts[0, 1]),
            color=_CELL_FILL_COLOR,
            alpha=alpha,
            edgecolor=_CELL_EDGE_COLOR,
            linewidth=1.2,
            label=cell_id,
        )

    if cluster.cells:
        all_verts = np.concatenate(
            [
                np.asarray(c.measurement_boundary.vertices_xy_um, dtype=np.float64)
                for c in cluster.cells.values()
            ],
            axis=0,
        )
        margin = max(1.0, 0.1 * float(np.ptp(all_verts)))
        ax.set_xlim(all_verts[:, 0].min() - margin, all_verts[:, 0].max() + margin)
        ax.set_ylim(all_verts[:, 1].min() - margin, all_verts[:, 1].max() + margin)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"v2 stub3d t={time_s:.3f} s · cells={len(cluster.cells)}")


def render_frame_png(
    input_path: str,
    output_path: str,
    *,
    figsize_inches: Optional[tuple[float, float]] = None,
    height_tint_scale_um: float = _HEIGHT_TINT_SCALE_UM,
) -> str:
    """Render a frame to PNG via matplotlib (Agg backend, no display)."""

    result = read_frame(input_path)
    abs_out = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(abs_out) or ".", exist_ok=True)
    fig, ax = plt.subplots(figsize=figsize_inches or _DEFAULT_FIGSIZE_INCHES)
    try:
        _render_axes(
            ax,
            result.cluster,
            result.time_s,
            height_tint_scale_um=height_tint_scale_um,
        )
        fig.tight_layout()
        fig.savefig(abs_out, dpi=120)
    finally:
        plt.close(fig)
    return abs_out


def _polygon_svg(verts: np.ndarray, *, fill: str, alpha: float) -> str:
    points = " ".join(f"{x:.4f},{y:.4f}" for x, y in verts)
    return (
        f'<polygon points="{points}" fill="{fill}" fill-opacity="{alpha:.3f}" '
        f'stroke="{_CELL_EDGE_COLOR}" stroke-width="0.05" />'
    )


def render_frame_html(
    input_path: str,
    output_path: str,
    *,
    height_tint_scale_um: float = _HEIGHT_TINT_SCALE_UM,
) -> str:
    """Render a frame to a self-contained HTML file (SVG + caption)."""

    result = read_frame(input_path)
    cluster = result.cluster
    abs_out = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(abs_out) or ".", exist_ok=True)

    if cluster.cells:
        all_verts = np.concatenate(
            [
                np.asarray(c.measurement_boundary.vertices_xy_um, dtype=np.float64)
                for c in cluster.cells.values()
            ],
            axis=0,
        )
        x_min, x_max = float(all_verts[:, 0].min()), float(all_verts[:, 0].max())
        y_min, y_max = float(all_verts[:, 1].min()), float(all_verts[:, 1].max())
    else:
        x_min, x_max, y_min, y_max = 0.0, 1.0, 0.0, 1.0
    span_x = max(1.0, x_max - x_min)
    span_y = max(1.0, y_max - y_min)
    margin = 0.1 * max(span_x, span_y)
    view_x = x_min - margin
    view_y = y_min - margin
    view_w = (x_max - x_min) + 2.0 * margin
    view_h = (y_max - y_min) + 2.0 * margin

    polygons: list[str] = []
    for cell in cluster.cells.values():
        verts = np.asarray(cell.measurement_boundary.vertices_xy_um, dtype=np.float64)
        flipped = np.column_stack(
            (verts[:, 0], view_y + view_h - (verts[:, 1] - view_y))
        )
        if cell.height_um is not None and height_tint_scale_um > 0.0:
            alpha = float(np.clip(cell.height_um / height_tint_scale_um, 0.15, 0.95))
        else:
            alpha = 0.55
        polygons.append(
            _polygon_svg(flipped, fill=_CELL_FILL_COLOR, alpha=alpha)
        )

    rows = "".join(
        f"<tr><td>{_html.escape(cid)}</td>"
        f"<td>{cell.measurement_boundary.projected_area_um2():.3f}</td>"
        f"<td>{cell.measurement_boundary.perimeter_um():.3f}</td>"
        f"<td>{'-' if cell.height_um is None else f'{cell.height_um:.3f}'}</td></tr>"
        for cid, cell in sorted(cluster.cells.items())
    )

    body = (
        f"<!doctype html><html><head><meta charset=\"utf-8\">"
        f"<title>v2 stub3d t={result.time_s:.3f} s</title>"
        f"<style>body{{font-family:sans-serif;margin:1em;}}"
        f"table{{border-collapse:collapse;margin-top:1em;}}"
        f"td,th{{border:1px solid #999;padding:0.25em 0.5em;}}"
        f"</style></head><body>"
        f"<h1>v2 stub3d frame</h1>"
        f"<p>time = {result.time_s:.3f} s &middot; cells = {len(cluster.cells)} "
        f"&middot; schema_version = {result.schema_version}</p>"
        f"<svg xmlns=\"http://www.w3.org/2000/svg\" "
        f"viewBox=\"{view_x:.4f} {view_y:.4f} {view_w:.4f} {view_h:.4f}\" "
        f"width=\"480\" height=\"480\">"
        f"<rect x=\"{view_x:.4f}\" y=\"{view_y:.4f}\" "
        f"width=\"{view_w:.4f}\" height=\"{view_h:.4f}\" "
        f"fill=\"#f0f0f0\" />"
        f"{''.join(polygons)}"
        f"</svg>"
        f"<table><thead><tr><th>cell_id</th><th>area_um2</th>"
        f"<th>perimeter_um</th><th>height_um</th></tr></thead>"
        f"<tbody>{rows}</tbody></table>"
        f"<p><em>Stub renderer; not confocal-quality.</em></p>"
        f"</body></html>"
    )

    with open(abs_out, "w", encoding="utf-8") as fh:
        fh.write(body)
    return abs_out
