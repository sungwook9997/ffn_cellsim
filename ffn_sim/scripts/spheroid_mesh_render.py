"""Render a spheroid npz as ACTUAL CELL MESHES (lit Poly3DCollection per cell) — full perspective +
an interior cross-section — so the PI can SEE the cells as solid polyhedra (round, all-touching,
flattened junctions), NOT a node scatter. (PI 2026-06-23: "메쉬 형태로 다시 시각화".)

Each cell's triangulated shell is drawn as a shaded surface (Lambert lighting), distinct colour per
cell. The cross-section keeps only triangles on the −y side of the COM so the interior packing /
flattened junctions are visible.

Run: python -m ffn_sim.scripts.spheroid_mesh_render <npz> --png out.png [--frame -1]
"""
from __future__ import annotations
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

UM = 1e6


def _shade(tris, base_rgb, az=-35.0, alt=55.0):
    """Lambert-shaded per-triangle rgba from a base colour + a fixed light direction."""
    a, b, c = tris[:, 0], tris[:, 1], tris[:, 2]
    n = np.cross(b - a, c - a)
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.clip(ln, 1e-30, None)
    azr, altr = np.radians(az), np.radians(alt)
    light = np.array([np.cos(altr) * np.cos(azr), np.cos(altr) * np.sin(azr), np.sin(altr)])
    inten = 0.35 + 0.65 * np.clip(np.abs(n @ light), 0.0, 1.0)
    rgb = np.clip(np.asarray(base_rgb)[None, :] * inten[:, None], 0, 1)
    return np.concatenate([rgb, np.ones((len(rgb), 1))], axis=1)


def render(npz, png, frame=-1):
    d = np.load(npz, allow_pickle=True)
    P = d["frames"][frame].astype(float) * UM           # → µm
    faces = d["faces"].astype(int); cof = d["cof"].astype(int)
    fcell = cof[faces[:, 0]]
    cells = np.unique(cof[cof >= 0])
    com = P[cof >= 0].mean(0)
    cmap = plt.get_cmap("tab20")

    cell_cen = {c: P[cof == c].mean(0) for c in cells}
    fig = plt.figure(figsize=(17, 8.5))
    # left = full spheroid (azim -60); right = interior: only the BACK-half cells (centroid y ≤ com),
    # drawn as full meshes and viewed from the FRONT (+y) so we look straight into the cut → the
    # interior packing + flattened cell-cell junctions are exposed.
    for col, (title, interior) in enumerate([("full spheroid (cell meshes)", False),
                                             ("interior — back-half cells, viewed into the cut", True)]):
        ax = fig.add_subplot(1, 2, col + 1, projection="3d")
        for k, c in enumerate(cells):
            if interior and cell_cen[c][1] > com[1]:     # drop front-half cells → expose interior
                continue
            f = faces[fcell == c]
            tris = P[f]
            base = np.array(cmap(k % 20)[:3])
            coll = Poly3DCollection(tris, facecolors=_shade(tris, base),
                                    edgecolors=(0, 0, 0, 0.12), linewidths=0.18)
            ax.add_collection3d(coll)
        Q = P[cof >= 0]
        rng = np.ptp(Q, axis=0).max() * 0.55
        ctr = Q.mean(0)
        ax.set_xlim(ctr[0] - rng, ctr[0] + rng); ax.set_ylim(ctr[1] - rng, ctr[1] + rng)
        ax.set_zlim(ctr[2] - rng, ctr[2] + rng)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=18, azim=(-60 if not interior else 90))   # interior: look down +y into the cut
        ax.set_title(title, fontsize=12); ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    tag = npz.split("/")[-1].replace(".npz", "")
    fig.suptitle(f"{tag} — {len(cells)} cells (actual meshes)  "
                 f"V/V0={float(d['vv0'][frame]):.3f}", fontsize=13)
    fig.tight_layout(); fig.savefig(png, dpi=130); plt.close(fig)
    print("saved", png)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz"); ap.add_argument("--png", required=True); ap.add_argument("--frame", type=int, default=-1)
    a = ap.parse_args()
    render(a.npz, a.png, a.frame)
