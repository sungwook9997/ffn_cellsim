"""Render a saved DCM spheroid npz as an MP4 of the ACTUAL CELL MESHES deforming over time — the
cells change shape (turgor + cohesion + contact deform the icospheres into packed polyhedra) while
the spheroid consolidates, and division daughters POP IN at their birth frame. This is the
shape-deformation-during-aggregation video the PI wants (markers/static frames don't show it).

Each frame redraws every active cell's triangulated shell as a lit Poly3DCollection (distinct colour),
with a slow camera rotation; frames are piped to ffmpeg (FFMpegWriter). Parked/unborn daughters
(centroid far from the cluster) are skipped per frame so they appear only once created.

Run: python -m ffn_sim.scripts.spheroid_mesh_mp4 <npz> --out out.mp4 [--max-frames 60] [--fps 12]
"""
from __future__ import annotations
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FFMpegWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

UM = 1e6


def _shade(tris, base_rgb, az=-35.0, alt=55.0):
    a, b, c = tris[:, 0], tris[:, 1], tris[:, 2]
    n = np.cross(b - a, c - a)
    n = n / np.clip(np.linalg.norm(n, axis=1, keepdims=True), 1e-30, None)
    azr, altr = np.radians(az), np.radians(alt)
    light = np.array([np.cos(altr) * np.cos(azr), np.cos(altr) * np.sin(azr), np.sin(altr)])
    inten = 0.4 + 0.6 * np.clip(np.abs(n @ light), 0.0, 1.0)
    rgb = np.clip(np.asarray(base_rgb)[None, :] * inten[:, None], 0, 1)
    return np.concatenate([rgb, np.ones((len(rgb), 1))], axis=1)


def render(npz, out, max_frames=60, fps=12, rotate=True):
    d = np.load(npz, allow_pickle=True)
    frames_all = d["frames"]
    faces = d["faces"].astype(int)
    cof = d["cof"].astype(int)
    fcell = cof[faces[:, 0]]
    step = d["step"]
    aa0 = d["aa0"] if "aa0" in d.files else None
    cells = np.unique(cof[cof >= 0])
    cmap = cm.get_cmap("tab20")
    cellcol = {int(c): np.array(cmap(k % 20)[:3]) for k, c in enumerate(cells)}
    face_of = {int(c): (fcell == c) for c in cells}

    nfr = frames_all.shape[0]
    idx = np.unique(np.linspace(0, nfr - 1, min(max_frames, nfr)).astype(int))

    # fixed axis box from the final active cluster (so the camera is stable)
    Pf = frames_all[-1].astype(float)[cof >= 0] * UM
    ctr = Pf.mean(0); rng = float(np.ptp(Pf, 0).max()) * 0.6

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    writer = FFMpegWriter(fps=fps, bitrate=4000)
    with writer.saving(fig, out, dpi=120):
        for fi, i in enumerate(idx):
            ax.clear()
            P = frames_all[i].astype(float) * UM
            for c in cells:
                ci = int(c)
                nodes = P[cof == c]
                if np.linalg.norm(nodes.mean(0) - ctr) > 500.0:   # parked/unborn → skip
                    continue
                tris = P[faces[face_of[ci]]]
                ax.add_collection3d(Poly3DCollection(
                    tris, facecolors=_shade(tris, cellcol[ci]),
                    edgecolors=(0, 0, 0, 0.10), linewidths=0.12))
            ax.set_xlim(ctr[0] - rng, ctr[0] + rng)
            ax.set_ylim(ctr[1] - rng, ctr[1] + rng)
            ax.set_zlim(ctr[2] - rng, ctr[2] + rng)
            ax.set_box_aspect((1, 1, 1))
            az = -65 + (fi * 1.5 if rotate else 0)
            ax.view_init(elev=14, azim=az)
            ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
            t = f"step {int(step[i])}"
            if aa0 is not None:
                t += f"   A/A0={float(aa0[i]):.3f}"
            nact = int(sum(np.linalg.norm(P[cof == c].mean(0) - ctr) <= 500.0 for c in cells))
            ax.set_title(f"{t}   cells={nact}", fontsize=11)
            writer.grab_frame()
            print(f"  frame {fi+1}/{len(idx)} (step {int(step[i])})", flush=True)
    plt.close(fig)
    print(f"saved {out} ({len(idx)} frames, {len(cells)} cells max)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-frames", type=int, default=60)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--no-rotate", action="store_true")
    a = ap.parse_args()
    render(a.npz, a.out, a.max_frames, a.fps, rotate=not a.no_rotate)
