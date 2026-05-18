"""2D / 3D rendering of fiber networks (no KU; supports Unit 1.1 figs)."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from acs_kb.ecm.fiber_network import FiberNetwork
from acs_kb.ecm.cross_links import CrossLink


def _orientation_colors(angles: np.ndarray) -> np.ndarray:
    """Map fiber director angles ∈ [0, π) to HSV-rgb colors."""
    hue = (np.mod(angles, np.pi) / np.pi).astype(float)
    hsv = np.stack([hue, np.full_like(hue, 0.85), np.full_like(hue, 0.85)], axis=1)
    return mcolors.hsv_to_rgb(hsv)


def plot_2d_static(
    net: FiberNetwork,
    cross_links: Sequence[CrossLink] | None = None,
    save_path: str | Path | None = None,
    show: bool = False,
    ax: plt.Axes | None = None,
    title: str | None = None,
) -> plt.Figure:
    """Plot the fiber network with matplotlib (one polyline per fiber).

    Fibers are colored by orientation (HSV). Cross-links shown as
    red dots at the midpoint of each linked bead pair.
    """
    L = net.box_size
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))
    else:
        fig = ax.figure

    colors = _orientation_colors(net.fiber_orientations)
    for f in range(net.bead_positions.shape[0]):
        pts = net.bead_positions[f]
        # Detect periodic wrap: split polyline if a bond jumps > L/2.
        deltas = np.diff(pts, axis=0)
        jumps = np.linalg.norm(deltas, axis=1) > 0.5 * L
        if not jumps.any():
            ax.plot(pts[:, 0], pts[:, 1], color=colors[f], lw=0.8)
        else:
            start = 0
            for k in np.where(jumps)[0]:
                ax.plot(pts[start : k + 1, 0], pts[start : k + 1, 1],
                        color=colors[f], lw=0.8)
                start = k + 1
            ax.plot(pts[start:, 0], pts[start:, 1], color=colors[f], lw=0.8)

    if cross_links:
        flat = net.bead_positions.reshape(-1, 2)
        N = net.bead_positions.shape[1]
        xs, ys = [], []
        for xl in cross_links:
            ia = xl.fiber_a * N + xl.bead_a
            ib = xl.fiber_b * N + xl.bead_b
            # midpoint with minimum image
            d = flat[ib] - flat[ia]
            d -= L * np.round(d / L)
            mid = np.mod(flat[ia] + 0.5 * d, L)
            xs.append(mid[0])
            ys.append(mid[1])
        ax.scatter(xs, ys, s=4, c="red", alpha=0.7, zorder=3, label=f"XL ({len(cross_links)})")
        ax.legend(loc="upper right", fontsize=8)

    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    if title:
        ax.set_title(title)
    fig.tight_layout()

    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150)
    if show:
        plt.show()
    return fig


def plot_3d_pyvista(
    net: FiberNetwork,
    cross_links: Sequence[CrossLink] | None = None,
    save_path: str | Path | None = None,
    off_screen: bool = True,
):
    """Render the 2D network as line segments in a 3D PyVista scene (z=0).

    Sets up the infrastructure that Unit 1.2 (full 3D) will reuse.
    """
    import pyvista as pv  # imported here so matplotlib-only callers don't pay VTK init.

    bp = net.bead_positions  # (F, N, 2)
    F_, N_ = bp.shape[:2]
    pts3 = np.concatenate([bp.reshape(-1, 2), np.zeros((F_ * N_, 1))], axis=1)

    cells = []
    for f in range(F_):
        seg = [N_] + list(range(f * N_, f * N_ + N_))
        cells.extend(seg)
    poly = pv.PolyData(pts3)
    poly.lines = np.array(cells, dtype=np.int64)

    p = pv.Plotter(off_screen=off_screen, window_size=(800, 800))
    p.add_mesh(poly, color="steelblue", line_width=2)

    if cross_links:
        flat = pts3
        L = net.box_size
        line_cells = []
        line_pts = []
        for xl in cross_links:
            ia = xl.fiber_a * N_ + xl.bead_a
            ib = xl.fiber_b * N_ + xl.bead_b
            d = flat[ib, :2] - flat[ia, :2]
            d -= L * np.round(d / L)
            if np.linalg.norm(d) > 0.5 * L:
                continue  # skip wrapping XLs in 3D view
            i0 = len(line_pts)
            line_pts.append(flat[ia])
            line_pts.append(np.concatenate([flat[ia, :2] + d, [0.0]]))
            line_cells.extend([2, i0, i0 + 1])
        if line_pts:
            xl_poly = pv.PolyData(np.array(line_pts))
            xl_poly.lines = np.array(line_cells, dtype=np.int64)
            p.add_mesh(xl_poly, color="red", line_width=1)

    p.view_xy()
    p.add_axes()
    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        p.screenshot(str(save_path))
    p.close()
    return save_path
