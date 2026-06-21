"""Mesh-surface visualization of the Warp de-cohesion spread (Phase-C migration).

Renders the cell SURFACE meshes (triangulated faces, colored per cell) — NOT point
clouds — across frames, in two views:
  * top-down (xy): the A/A0 silhouette (footprint) — what spreading grows.
  * side (xz): the height profile — shows maxZ HELD while the footprint widens
    (genuine lateral spread, not the vertical-flatten settle-confound).

Integrity (the viz rules + the landmine lessons):
  * mesh (filled triangles), not points; equal aspect (true proportions, no distortion);
  * axis limits FIXED across frames (from the max extent) so the spread is comparable —
    NOT per-frame autoscale (which would hide it);
  * A/A0 annotation = the driver's TOP-DOWN silhouette value (never basal contact);
  * SI→µm, config (no body-force proxy) stated on the figure.

Input: the npz written by `dcm_warp_decohesion.py --save-frames`.

    python -m ffn_sim.scripts.warp_decoh_mesh_viz --npz frames.npz --out fig.png
"""

from __future__ import annotations

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

UM = 1.0e6


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ncols", type=int, default=5)
    args = ap.parse_args()

    d = np.load(args.npz)
    frames = d["frames"]                 # (F, N, 3) metres
    faces = d["faces"]                   # (M, 3) node indices
    cof = d["cof"]                       # (N,) cell-of-node
    step, aa0, maxZ, vv0 = d["step"], d["aa0"], d["maxZ"], d["vv0"]
    F = frames.shape[0]

    sel = np.unique(np.linspace(0, F - 1, args.ncols).astype(int))
    face_cell = cof[faces[:, 0]]
    fcolors = plt.get_cmap("tab20")(face_cell % 20)
    fcolors[:, 3] = 0.55                 # translucent so the ball's depth reads

    allp = frames * UM
    xylim = float(np.abs(allp[..., :2]).max()) * 1.05      # FIXED across frames
    zmin, zmax = float(allp[..., 2].min()), float(allp[..., 2].max())

    def tri(P, i, j):                    # (M,3,2) triangle polygons for projection (i,j)
        return P[faces][:, :, [i, j]]

    fig, axes = plt.subplots(2, len(sel), figsize=(3.2 * len(sel), 6.8))
    if len(sel) == 1:
        axes = axes.reshape(2, 1)
    for c, fi in enumerate(sel):
        P = frames[fi] * UM
        ax = axes[0, c]
        ax.add_collection(PolyCollection(tri(P, 0, 1), facecolors=fcolors,
                                         edgecolors=(0, 0, 0, 0.12), linewidths=0.1))
        ax.set_xlim(-xylim, xylim); ax.set_ylim(-xylim, xylim); ax.set_aspect("equal")
        ax.set_title(f"step {int(step[fi])}\nA/A0 = {aa0[fi]:.2f}", fontsize=9)
        ax.tick_params(labelsize=6)
        if c == 0:
            ax.set_ylabel("top-down (xy)  µm\n[A/A0 silhouette]", fontsize=8)

        ax = axes[1, c]
        ax.add_collection(PolyCollection(tri(P, 0, 2), facecolors=fcolors,
                                         edgecolors=(0, 0, 0, 0.12), linewidths=0.1))
        ax.set_xlim(-xylim, xylim); ax.set_ylim(zmin - 2.0, zmax + 4.0); ax.set_aspect("equal")
        ax.axhline(0.0, color="saddlebrown", lw=1.2, alpha=0.7)   # substrate plane z0
        ax.set_title(f"maxZ = {maxZ[fi]:.0f} µm\nV/V0 = {vv0[fi]:.2f}", fontsize=9)
        ax.tick_params(labelsize=6)
        if c == 0:
            ax.set_ylabel("side (xz)  µm\n[maxZ held]", fontsize=8)

    fig.suptitle(
        "Warp DCM engine — de-cohesion N=100 wetting cleanball (gbook A5000)\n"
        "cell SURFACE mesh, colored per cell · substrate wetting only, NO body-force proxy, "
        "physiological γ, dt=8e-6\n"
        "top-down (A/A0 silhouette) GROWS while side (maxZ) HELD → genuine lateral spread",
        fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(args.out, dpi=130)
    print(f"wrote {args.out}  ({len(sel)} frames, {faces.shape[0]} triangles/frame, "
          f"A/A0 {aa0[0]:.2f}->{aa0[-1]:.2f}, maxZ {maxZ[0]:.0f}->{maxZ[-1]:.0f}µm)")


if __name__ == "__main__":
    main()
