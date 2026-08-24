"""Render a spheroid/spreading run's saved frames as a time montage (top-down xy + side xz),
straight from a --save-frames npz. Actual cell nodes coloured by cell. The picture the PI judges
by — and a cross-section/side row so flattening or peeling is VISIBLE (a 3D trisurf hides both).

Run: python -m aleph.scripts.spheroid_frames_montage <frames.npz> --png out.png [--ncol 6]
"""
from __future__ import annotations
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

UM = 1e6


def montage(npz, png, ncol=6):
    d = np.load(npz, allow_pickle=True)
    F = d["frames"].astype(np.float64)
    cof = d["cof"].astype(int)
    live = cof >= 0
    nf = len(F)
    aa0 = d["aa0"] if "aa0" in d else np.full(nf, np.nan)
    maxZ = d["maxZ"] if "maxZ" in d else np.full(nf, np.nan)
    vv0 = d["vv0"] if "vv0" in d else np.full(nf, np.nan)
    idx = np.linspace(0, nf - 1, min(ncol, nf)).round().astype(int)
    # global limits from last frame for consistent scale
    com = F[-1][live].mean(0)
    span = np.abs((F[-1][live] - com) * UM).max() * 1.1
    fig, axs = plt.subplots(2, len(idx), figsize=(3.2 * len(idx), 6.6))
    if len(idx) == 1:
        axs = axs.reshape(2, 1)
    col = cof[live] % 20
    for j, fi in enumerate(idx):
        P = (F[fi] - com) * UM
        # top-down xy
        axs[0, j].scatter(P[live, 0], P[live, 1], s=2, c=col, cmap="tab20", alpha=0.6)
        axs[0, j].set_aspect("equal"); axs[0, j].set_xlim(-span, span); axs[0, j].set_ylim(-span, span)
        axs[0, j].set_title(f"f{fi}  A/A0={aa0[fi]:.2f}\nV/V0={vv0[fi]:.3f}", fontsize=9)
        axs[0, j].set_xticks([]); axs[0, j].set_yticks([])
        # side xz
        axs[1, j].scatter(P[live, 0], P[live, 2], s=2, c=col, cmap="tab20", alpha=0.6)
        axs[1, j].set_aspect("equal"); axs[1, j].set_xlim(-span, span); axs[1, j].set_ylim(-span, span)
        axs[1, j].set_title(f"side  maxZ={maxZ[fi]:.0f}µm", fontsize=9)
        axs[1, j].set_xticks([]); axs[1, j].set_yticks([])
    axs[0, 0].set_ylabel("top-down (xy)"); axs[1, 0].set_ylabel("side (xz)")
    fig.suptitle(f"{npz.split('/')[-1]}  ({np.unique(cof[live]).size} cells, {nf} frames)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(png, dpi=110); plt.close(fig)
    print("saved", png)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--png", required=True)
    ap.add_argument("--ncol", type=int, default=6)
    a = ap.parse_args()
    montage(a.npz, a.png, a.ncol)
