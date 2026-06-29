"""Top-down + side polygon-mesh view of a saved DCM spheroid (.npy), longer xy axes.

Renders each cell as a shaded triangulated surface (SimuCell3D style), coloured by
per-cell V/V0, with an EXTENDED xy range so structure + later spreading is visible
on a common, generous scale. Output → outputs/h_dcm_two_stage/figs/.

Usage: python -m ffn_sim.scripts.dcm_spheroid_mesh_viz --npy /tmp/sph.npy \
           --subdiv 2 --out figs/<name>.png --xy-um 50 --title "..."
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from ffn_sim.archive.hoomd_legacy.cell.dcm import icosphere_mesh

UM = 1e6
OUT = Path("ffn_sim/outputs/h_dcm_two_stage")


def _vg(P, tris):
    a, b, c = P[tris[:, 0]], P[tris[:, 1]], P[tris[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b - a, c - a)).sum() / 6.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npy", required=True)
    ap.add_argument("--subdiv", type=int, default=2)
    ap.add_argument("--r-cell-um", type=float, default=7.5)
    ap.add_argument("--xy-um", type=float, default=50.0,
                    help="HALF-width of the (extended) xy axes [µm].")
    ap.add_argument("--out", default="figs/spheroid_mesh.png")
    ap.add_argument("--title", default="DCM spheroid (subdiv-2)")
    args = ap.parse_args()

    v0i, _, tris0 = icosphere_mesh(args.r_cell_um * 1e-6, args.subdiv)
    nv = v0i.shape[0]
    pos = np.load(args.npy).reshape(-1, nv, 3)
    ncell = pos.shape[0]
    V0 = _vg(v0i, tris0)
    vv = np.array([_vg(pos[c], tris0) / V0 for c in range(ncell)])
    cmap = plt.get_cmap("tab20").colors
    cen = pos.reshape(-1, 3).mean(0) * UM
    h = args.xy_um

    fig = plt.figure(figsize=(14, 6.4))
    for k, (elev, azim, ttl) in enumerate([(90, -90, "TOP-DOWN (xy)"),
                                           (10, -72, "SIDE (xz)")]):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        for c in range(ncell):
            verts = pos[c] * UM
            base = np.array(cmap[c % len(cmap)])
            tri = verts[tris0]
            n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
            nn = np.linalg.norm(n, axis=1, keepdims=True)
            sh = 0.55 + 0.45 * np.clip((n / np.where(nn > 0, nn, 1)) @ np.array([0.3, 0.3, 0.9]), 0, 1)
            rgba = np.zeros((tris0.shape[0], 4)); rgba[:, :3] = base[None, :] * sh[:, None]; rgba[:, 3] = 0.9
            ax.add_collection3d(Poly3DCollection(tri, facecolors=rgba,
                                edgecolors=(0, 0, 0, 0.15), linewidths=0.1))
        ax.set_xlim(cen[0] - h, cen[0] + h); ax.set_ylim(cen[1] - h, cen[1] + h)
        zmax = max(20.0, pos[:, :, 2].max() * UM * 1.1)
        ax.set_zlim(-2, zmax)
        ax.set_box_aspect((2 * h, 2 * h, (zmax + 2) if k == 1 else 0.5 * h))
        ax.view_init(elev=elev, azim=azim)
        ax.set_xlabel("x[µm]", fontsize=8); ax.set_ylabel("y[µm]", fontsize=8)
        ax.set_title(ttl, fontsize=10)
    fig.suptitle(f"{args.title}\nN={ncell} cells, subdiv-{args.subdiv} ({nv} nodes/cell) · "
                 f"per-cell V/V0 = {vv.mean():.3f} ± {vv.std():.3f} (min {vv.min():.3f}, max {vv.max():.3f})",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = OUT / args.out
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=140)
    print(f"wrote {p} | V/V0 {vv.mean():.3f}±{vv.std():.3f}")


if __name__ == "__main__":
    main()
