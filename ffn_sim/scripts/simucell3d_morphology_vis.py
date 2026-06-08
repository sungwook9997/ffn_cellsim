#!/usr/bin/env python3
"""Render SimuCell3D MCF7 P0 morphology (deformable cell surfaces) to PNG.

Reads SimuCell3D face_data VTK frames (triangulated cell surfaces + per-face
owner cell id) and renders the spheroid with matplotlib, cells distinctly
coloured, across a few timepoints to show the 1->N division morphology.
Part of the ffn_cellsim -> SimuCell3D pipeline (the foundry->engine demo).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402


def parse_face_vtk(path: Path):
    """Minimal VTK-legacy reader: POINTS, triangle CELLS, face_cell_id field."""
    txt = path.read_text().split("\n")
    pts = tris = cid = None
    i = 0
    while i < len(txt):
        line = txt[i]
        if line.startswith("POINTS"):
            n = int(line.split()[1])
            vals: list[float] = []
            i += 1
            while len(vals) < 3 * n:
                vals += [float(x) for x in txt[i].split()]
                i += 1
            pts = np.array(vals[: 3 * n]).reshape(-1, 3)
            continue
        if line.startswith("CELLS"):
            nc = int(line.split()[1])
            rows = []
            i += 1
            while len(rows) < nc:
                parts = txt[i].split()
                i += 1
                if parts:
                    rows.append([int(parts[1]), int(parts[2]), int(parts[3])])
            tris = np.array(rows)
            continue
        if line.startswith("face_cell_id"):
            nc = int(line.split()[2])
            vals = []
            i += 1
            while len(vals) < nc:
                vals += [float(x) for x in txt[i].split()]
                i += 1
            cid = np.array(vals[:nc]).astype(int)
            continue
        i += 1
    return pts, tris, cid


def render_frame(ax, path: Path) -> int:
    pts, tris, cid = parse_face_vtk(path)
    verts = pts[tris] * 1e6  # m -> µm
    uids = np.unique(cid)
    cmap = plt.cm.tab20
    cmap_for = {u: cmap(j % 20) for j, u in enumerate(uids)}
    pc = Poly3DCollection(
        verts,
        facecolors=[cmap_for[c] for c in cid],
        edgecolors=(0, 0, 0, 0.06),
        linewidths=0.1,
    )
    ax.add_collection3d(pc)
    flat = verts.reshape(-1, 3)
    ctr = flat.mean(0)
    rad = np.linalg.norm(flat - ctr, axis=1).max()
    for setlim, c0 in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), ctr):
        setlim(c0 - rad, c0 + rad)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=40)
    ax.set_xlabel("µm", fontsize=8)
    return len(uids)


def main() -> None:
    D = Path("/Users/sw1/simucell3d_ref/simulation_results/mcf7_p0/face_data")
    frames = sorted(D.glob("result_*.vtk"), key=lambda p: int(p.stem.split("_")[1]))
    n = len(frames)
    picks = [frames[0], frames[n // 3], frames[2 * n // 3], frames[-1]]

    fig = plt.figure(figsize=(16, 4.6))
    for k, fp in enumerate(picks):
        ax = fig.add_subplot(1, 4, k + 1, projection="3d")
        ncell = render_frame(ax, fp)
        ax.set_title(f"frame {int(fp.stem.split('_')[1])} — {ncell} cell(s)", fontsize=10)
    fig.suptitle(
        "SimuCell3D MCF7 P0: ffn closures (γ=0.27 mN/m Hosseini cortical tension + "
        "cadherin cohesion) → deformable-cell spheroid morphology",
        fontsize=11,
    )
    out = Path("/Users/sw1/ffn_cellsim/ffn_sim/outputs/simucell3d/mcf7_p0_morphology.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=130)
    print("saved", out)


if __name__ == "__main__":
    main()
