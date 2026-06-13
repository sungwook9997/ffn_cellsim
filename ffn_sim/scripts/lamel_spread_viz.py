"""Top-down polygon-mesh visualization of a mechanistic-lamellipodium spread.

PI requirement (2026-06-12/14): render the spread as a POLYGON MESH (SimuCell3D
style, every cell a shaded triangulated surface) viewed TOP-DOWN (xy silhouette =
the experimental A/A₀ assay), NOT a point scatter and NOT only z-expansion. Reads
a ``two_stage_n{N}.pkl`` from ``dcm_two_stage_production --lamellipodium`` and emits:

  two_stage_n{N}_lamel_topdown.png  — top-down mesh: spread start vs end + curves
  two_stage_n{N}_lamel_spread.mp4   — top-down mesh spreading movie

Usage:  python -m ffn_sim.scripts.lamel_spread_viz --n 400 [--pkl path]
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.animation import FFMpegWriter

UM = 1e6
OUT = Path("ffn_sim/outputs/h_dcm_two_stage")
_TAB = plt.get_cmap("tab20").colors


def _lambert(verts, tris, base_rgb, ls, alpha):
    tri = verts[tris]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    nn = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.where(nn > 0, n / nn, n)
    shade = 0.55 + 0.45 * np.clip(n @ np.array([0.3, 0.3, 0.9]), 0, 1)
    rgba = np.zeros((tris.shape[0], 4))
    rgba[:, :3] = np.array(base_rgb)[None, :] * shade[:, None]
    rgba[:, 3] = alpha
    return rgba


def draw_topdown(ax, pos, ranges, tris0, ls):
    """Top-down (looking down z) shaded polygon mesh, one colour per cell."""
    for c, (lo, hi) in enumerate(ranges):
        verts = pos[lo:hi]
        rgba = _lambert(verts, tris0, to_rgb(_TAB[c % len(_TAB)]), ls, 0.9)
        ax.add_collection3d(Poly3DCollection(
            verts[tris0] * UM, facecolors=rgba,
            edgecolors=(0, 0, 0, 0.12), linewidths=0.1))
    ax.view_init(elev=90, azim=-90)        # straight down z → xy silhouette


def _xybox(frames, ranges):
    nmem = ranges[-1][1]
    allp = np.concatenate([f[:nmem] for f in frames], axis=0) * UM
    cx, cy = allp[:, 0].mean(), allp[:, 1].mean()
    r = 0.5 * max(np.ptp(allp[:, 0]), np.ptp(allp[:, 1])) * 1.1
    return (cx - r, cx + r), (cy - r, cy + r)


def _setup_top(ax, xlim, ylim):
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_zlim(-20, 200)
    ax.set_box_aspect((xlim[1] - xlim[0], ylim[1] - ylim[0],
                       0.25 * (xlim[1] - xlim[0])))
    ax.view_init(elev=90, azim=-90)
    ax.set_xlabel("x[µm]", fontsize=7); ax.set_ylabel("y[µm]", fontsize=7)
    ax.set_zticks([])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--pkl", default=None)
    args = ap.parse_args()
    path = Path(args.pkl) if args.pkl else OUT / f"two_stage_n{args.n}.pkl"
    d = pickle.load(open(path, "rb"))
    ranges, tris0 = d["ranges"], d["tris0"]
    s = d["spread"]
    A0 = s["A0_topdown_um2"]
    ls = LightSource(azdeg=225, altdeg=65)
    xb, yb = _xybox(s["frames"], ranges)

    aa = [g["topdown_um2"] / A0 for g in s["diags"]]
    vv = [g["VV0_mean"] for g in s["diags"]]
    mz = [g["maxZ_um"] for g in s["diags"]]
    steps = s["steps"]

    fig = plt.figure(figsize=(16, 5.2))
    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    draw_topdown(ax1, s["frames"][0], ranges, tris0, ls); _setup_top(ax1, xb, yb)
    ax1.set_title(f"spread start (top-down)\nA/A0=1.00  V/V0={vv[0]:.2f}", fontsize=9)
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    draw_topdown(ax2, s["frames"][-1], ranges, tris0, ls); _setup_top(ax2, xb, yb)
    ax2.set_title(f"spread end (top-down)\nA/A0={aa[-1]:.2f}  V/V0={vv[-1]:.2f}", fontsize=9)
    ax3 = fig.add_subplot(1, 3, 3)
    ax3.plot(steps, aa, "-o", ms=3, color="C2", label="A/A0 (top-down, the assay)")
    ax3.axhline(1.0, color="0.6", lw=0.8, ls=":")
    ax3.set_ylabel("A/A0", color="C2"); ax3.tick_params(axis="y", colors="C2")
    ax3.set_xlabel("spread step")
    axb = ax3.twinx()
    axb.plot(steps, vv, "-s", ms=3, color="C0", label="V/V0 (volume held?)")
    axb.plot(steps, np.array(mz) / max(mz), "-^", ms=3, color="C3",
             label="maxZ (norm; flattening?)")
    axb.set_ylabel("V/V0 · maxZ(norm)", color="C0")
    axb.axhline(1.0, color="0.6", lw=0.8, ls=":")
    h1, l1 = ax3.get_legend_handles_labels(); h2, l2 = axb.get_legend_handles_labels()
    ax3.legend(h1 + h2, l1 + l2, fontsize=7, loc="best")
    ax3.set_title("A/A0 should RISE while V/V0 stays ~1 (spread, not squash)", fontsize=8)
    fig.suptitle(f"DCM mechanistic lamellipodium spread · N={args.n} · TOP-DOWN mesh · "
                 f"A/A0 {1.0:.2f}→{aa[-1]:.2f}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    png = OUT / f"two_stage_n{args.n}_lamel_topdown.png"
    fig.savefig(png, dpi=140); plt.close(fig)
    print(f"wrote {png}", flush=True)

    # top-down spreading movie
    fig = plt.figure(figsize=(7, 6.6))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    writer = FFMpegWriter(fps=6, bitrate=2800)
    mp4 = OUT / f"two_stage_n{args.n}_lamel_spread.mp4"
    with writer.saving(fig, str(mp4), dpi=115):
        for pos, dg, st in zip(s["frames"], s["diags"], steps):
            ax.clear()
            draw_topdown(ax, pos, ranges, tris0, ls); _setup_top(ax, xb, yb)
            ax.set_title(f"top-down spread · step {st} · A/A0="
                         f"{dg['topdown_um2']/A0:.2f} · V/V0={dg['VV0_mean']:.2f} · "
                         f"maxZ={dg['maxZ_um']:.0f}µm", fontsize=9)
            writer.grab_frame()
    print(f"wrote {mp4}", flush=True)


if __name__ == "__main__":
    main()
