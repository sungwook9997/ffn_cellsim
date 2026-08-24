"""Mid-plane SLAB cross-section of a spheroid — the definitive "solid tissue vs marble pile" check.

For a thin slab at the spheroid mid-height, draw EACH cell's cross-section (2D convex hull of its
nodes in the slab) as a filled coloured polygon. Reading:
  • cells tile TIGHTLY with shared straight edges + no gaps → flattened tissue (real spheroid).
  • circular cells with GAPS between them → loose "marble pile".
Also reports the slab AREA-FILL fraction (covered cell-area / hull-of-all-area) and a per-cell
polygonality proxy. (PI 2026-06-23: confirm the spheroid is solid before proceeding.)

Run: python -m aleph.scripts.spheroid_slab_slice <npz> --png out.png
"""
from __future__ import annotations
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull

UM = 1e6


def render(npz, png, n_slabs=3):
    d = np.load(npz, allow_pickle=True)
    P = d["frames"][-1].astype(float) * UM
    cof = d["cof"].astype(int)
    live = cof >= 0
    cells = np.unique(cof[live])
    com = P[live].mean(0)
    cmap = plt.get_cmap("tab20")
    zspan = np.ptp(P[live, 2])
    delta = 0.10 * zspan / 2.0                       # thin slab half-thickness (~10% of height/2)
    offsets = np.linspace(-0.22 * zspan, 0.22 * zspan, n_slabs)

    fig, axs = plt.subplots(1, n_slabs, figsize=(6.2 * n_slabs, 6.4))
    if n_slabs == 1:
        axs = [axs]
    for ax, off in zip(axs, offsets):
        z0 = com[2] + off
        polys = []; covered = 0.0; total_nodes = 0
        for k, c in enumerate(cells):
            nm = (cof == c) & (np.abs(P[:, 2] - z0) < delta)
            pts = P[nm][:, :2]
            if len(pts) < 3:
                continue
            try:
                h = ConvexHull(pts)
            except Exception:
                continue
            poly = pts[h.vertices]
            polys.append((poly, cmap(k % 20)))
            covered += h.volume                       # 2D hull "volume" = area
            total_nodes += len(pts)
        for poly, color in polys:
            ax.add_patch(Polygon(poly, closed=True, facecolor=color, edgecolor="k",
                                 linewidth=0.6, alpha=0.85))
        allpts = P[live & (np.abs(P[:, 2] - z0) < delta)][:, :2]
        fill = 0.0
        if len(allpts) >= 3:
            envel = ConvexHull(allpts).volume
            fill = covered / envel if envel > 0 else 0.0
            ax.set_xlim(allpts[:, 0].min() - 5, allpts[:, 0].max() + 5)
            ax.set_ylim(allpts[:, 1].min() - 5, allpts[:, 1].max() + 5)
        ax.set_aspect("equal")
        ax.set_title(f"slab z={z0-com[2]:+.0f}µm  ({len(polys)} cells)\n"
                     f"area-fill={fill:.2f} (1.0=gapless tissue)", fontsize=11)
        ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    tag = npz.split("/")[-1].replace(".npz", "")
    fig.suptitle(f"{tag} — mid-plane SLAB cross-sections (cell polygons): tight tiling = solid tissue, "
                 f"gaps = marble pile", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(png, dpi=130); plt.close(fig)
    print("saved", png)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz"); ap.add_argument("--png", required=True)
    a = ap.parse_args()
    render(a.npz, a.png)
