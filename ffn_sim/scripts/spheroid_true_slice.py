"""TRUE plane-section of a spheroid: cut every cell's triangulated mesh with a horizontal plane
(exact triangle-plane intersection, Ericson) → each cell's real cross-section polygon. Unlike the
crude convex-hull-of-nearby-nodes version, this draws the ACTUAL cell outline at the cut, so flat
shared junctions, true gaps, and real overlaps are faithful (no hull artifact).

Run: python -m ffn_sim.scripts.spheroid_true_slice <npz> --png out.png
"""
from __future__ import annotations
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

UM = 1e6


def _cell_section(P, faces, z0):
    """Intersection segments of a triangle soup with plane z=z0 → ordered polygon (xy)."""
    segs = []
    for tri in faces:
        v = P[tri]                                   # (3,3)
        z = v[:, 2] - z0
        pts = []
        for i in range(3):
            a, b = v[i], v[(i + 1) % 3]
            za, zb = z[i], z[(i + 1) % 3]
            if (za <= 0 < zb) or (zb <= 0 < za):     # edge crosses the plane
                t = za / (za - zb)
                pts.append((a + t * (b - a))[:2])
        if len(pts) == 2:
            segs.append(pts)
    if not segs:
        return None
    P2 = np.array([p for s in segs for p in s])
    c = P2.mean(0)
    ang = np.arctan2(P2[:, 1] - c[1], P2[:, 0] - c[0])
    return P2[np.argsort(ang)]                        # order around centroid → closed polygon


def render(npz, png):
    d = np.load(npz, allow_pickle=True)
    P = d["frames"][-1].astype(float) * UM
    faces = d["faces"].astype(int); cof = d["cof"].astype(int)
    fcell = cof[faces[:, 0]]
    cells = np.unique(cof[cof >= 0])
    com = P[cof >= 0].mean(0)
    cmap = plt.get_cmap("tab20")
    zspan = np.ptp(P[cof >= 0, 2])
    offsets = np.linspace(-0.18 * zspan, 0.18 * zspan, 3)

    fig, axs = plt.subplots(1, 3, figsize=(18, 6.2))
    for ax, off in zip(axs, offsets):
        z0 = com[2] + off
        n = 0
        for k, c in enumerate(cells):
            poly = _cell_section(P[None][0], faces[fcell == c], z0)
            if poly is None or len(poly) < 3:
                continue
            ax.add_patch(Polygon(poly, closed=True, facecolor=cmap(k % 20),
                                 edgecolor="k", linewidth=0.7, alpha=0.9))
            n += 1
        live = P[cof >= 0]
        m = np.abs(live[:, 2] - z0) < 0.10 * zspan
        if m.any():
            ax.set_xlim(live[m, 0].min() - 5, live[m, 0].max() + 5)
            ax.set_ylim(live[m, 1].min() - 5, live[m, 1].max() + 5)
        ax.set_aspect("equal")
        ax.set_title(f"z={off:+.0f}µm  ({n} cells sectioned)", fontsize=12)
        ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    tag = npz.split("/")[-1].replace(".npz", "")
    fig.suptitle(f"{tag} — TRUE mesh plane-sections (exact cell outlines at the cut)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(png, dpi=135); plt.close(fig)
    print("saved", png)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz"); ap.add_argument("--png", required=True)
    a = ap.parse_args()
    render(a.npz, a.png)
