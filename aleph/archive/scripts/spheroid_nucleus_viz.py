"""Standard spheroid cross-section visualization WITH the nucleus rendered inside each cell.

PI directive (2026-06-23): going forward, every spheroid view must show the nucleus (the deformable
core, radius R_nuc = R_nuc_factor·R_cell) inside the cell — not just the membrane. The DCM nucleus is
a virtual core (per-cell centroid + R_nuc resisting membrane intrusion), so we render it as a filled
circle of radius R_nuc at each cell's centroid wherever that centroid sits within the cross-section slab.

Renders an equatorial mid-plane SLAB (render-independent, the trustworthy view per CONTACT_RESOLVED —
NOT the 3D trisurf that hides flat junctions): membrane node-cloud + per-cell convex outline + nucleus
circles. Reads the committed `frames/faces/cof` npz that the de-cohesion driver saves.

Usage:
  python -m aleph.scripts.spheroid_nucleus_viz <npz> [--out fig.png] [--frame -1] [--rcell 7.5e-6] [--rnuc-factor 0.33]
"""
from __future__ import annotations
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection


def render(npz_path, out=None, frame=-1, r_cell=7.5e-6, rnuc_factor=0.33, slab_factor=0.6):
    d = np.load(npz_path, allow_pickle=True)
    pos = np.asarray(d["frames"][frame], dtype=float)            # (N,3) µm-scale SI
    cof = np.asarray(d["cof"], dtype=int)
    um = 1e6
    r_nuc = rnuc_factor * r_cell
    cells = np.unique(cof[cof >= 0])
    com = pos[cof >= 0].mean(0)
    z0 = float(np.median(pos[cof >= 0, 2]))                      # equatorial plane
    dz = slab_factor * r_cell                                    # slab half-thickness
    in_slab = (np.abs(pos[:, 2] - z0) < dz) & (cof >= 0)

    fig, ax = plt.subplots(figsize=(9, 9))
    rng = np.random.default_rng(0)
    cmap = plt.get_cmap("tab20")
    nuc_patches = []
    n_nuc = 0
    try:
        from scipy.spatial import ConvexHull
        have_hull = True
    except Exception:
        have_hull = False

    for ci in cells:
        m = (cof == ci)
        cell_nodes = pos[m]
        cz = cell_nodes[:, 2].mean()
        ccx, ccy = cell_nodes[:, 0].mean(), cell_nodes[:, 1].mean()
        col = cmap(int(ci) % 20)
        # membrane cross-section: in-slab nodes of this cell
        sm = m & in_slab
        P = pos[sm][:, :2] * um
        if P.shape[0] >= 1:
            ax.scatter(P[:, 0], P[:, 1], s=2, color=col, alpha=0.35, linewidths=0)
        if have_hull and P.shape[0] >= 3:
            try:
                h = ConvexHull(P)
                poly = P[h.vertices]
                ax.fill(poly[:, 0], poly[:, 1], color=col, alpha=0.10, lw=0)
                ax.plot(np.append(poly[:, 0], poly[0, 0]), np.append(poly[:, 1], poly[0, 1]),
                        color=col, lw=0.6, alpha=0.55)
            except Exception:
                pass
        # NUCLEUS: filled circle radius R_nuc at the cell centroid, if its centroid is within the slab
        if abs(cz - z0) < r_nuc:
            # the visible nucleus radius shrinks for centroids off the mid-plane (sphere ∩ plane)
            rr = np.sqrt(max(r_nuc**2 - (cz - z0) ** 2, 0.0)) * um
            nuc_patches.append(Circle((ccx * um, ccy * um), rr))
            n_nuc += 1
    if nuc_patches:
        pc = PatchCollection(nuc_patches, facecolor="#2c3e50", edgecolor="#e74c3c",
                             alpha=0.55, lw=0.8)
        ax.add_collection(pc)

    aa0 = float(d["aa0"][frame]) if "aa0" in d else float("nan")
    vv0 = float(d["vv0"][frame]) if "vv0" in d else float("nan")
    pen = float(d["pen_frac"][frame]) if "pen_frac" in d else float("nan")
    step = int(d["step"][frame]) if "step" in d else frame
    ax.set_aspect("equal")
    R = np.linalg.norm(pos[cof >= 0, :2] - com[:2], axis=1).max() * um
    ax.set_xlim((com[0] * um - 1.1 * R, com[0] * um + 1.1 * R))
    ax.set_ylim((com[1] * um - 1.1 * R, com[1] * um + 1.1 * R))
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    ax.set_title(f"{os.path.basename(npz_path)}  step {step}  ({len(cells)} cells)\n"
                 f"equatorial slab ±{dz*um:.1f}µm · {n_nuc} nuclei in plane (red-ring cores, R_nuc={r_nuc*um:.1f}µm)\n"
                 f"A/A0={aa0:.3f}  V/V0={vv0:.3f}  pen={pen:.2f}", fontsize=10)
    fig.tight_layout()
    if out is None:
        out = os.path.join(os.path.dirname(npz_path), "figs",
                           os.path.splitext(os.path.basename(npz_path))[0] + "_nucleus_slab.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"saved {out}  ({len(cells)} cells, {n_nuc} nuclei in slab)")
    return out


def render_montage(npz_path, out=None, frame=-1, n_slices=6, r_cell=7.5e-6, rnuc_factor=0.33):
    """Multi-slice z-montage of the spheroid WITH nuclei — shows the full 3D internal structure
    (a single equatorial slab hides the polar caps). ``n_slices`` z-levels from bottom to top, each
    a thin slab rendered as cells (convex outline) + nucleus circles (sphere-plane intersection)."""
    d = np.load(npz_path, allow_pickle=True)
    pos = np.asarray(d["frames"][frame], dtype=float)
    cof = np.asarray(d["cof"], dtype=int)
    um = 1e6; r_nuc = rnuc_factor * r_cell
    cells = np.unique(cof[cof >= 0])
    live = cof >= 0
    com = pos[live].mean(0)
    zmin, zmax = np.percentile(pos[live, 2], 3), np.percentile(pos[live, 2], 97)
    z_levels = np.linspace(zmin, zmax, n_slices)
    dz = (zmax - zmin) / (2.2 * n_slices)
    R = np.linalg.norm(pos[live, :2] - com[:2], axis=1).max() * um
    try:
        from scipy.spatial import ConvexHull
        have_hull = True
    except Exception:
        have_hull = False
    cmap = plt.get_cmap("tab20")
    ncol = 3; nrow = int(np.ceil(n_slices / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.2 * nrow))
    axes = np.atleast_1d(axes).ravel()
    cell_z = {ci: pos[cof == ci].mean(0) for ci in cells}
    for s, z0 in enumerate(z_levels):
        ax = axes[s]
        in_slab = (np.abs(pos[:, 2] - z0) < dz) & live
        nuc = []
        for ci in cells:
            m = cof == ci
            col = cmap(int(ci) % 20)
            P = pos[m & in_slab][:, :2] * um
            if P.shape[0] >= 1:
                ax.scatter(P[:, 0], P[:, 1], s=2, color=col, alpha=0.3, linewidths=0)
            if have_hull and P.shape[0] >= 3:
                try:
                    h = ConvexHull(P); poly = P[h.vertices]
                    ax.fill(poly[:, 0], poly[:, 1], color=col, alpha=0.10, lw=0)
                    ax.plot(np.append(poly[:, 0], poly[0, 0]), np.append(poly[:, 1], poly[0, 1]), color=col, lw=0.5, alpha=0.5)
                except Exception:
                    pass
            cc = cell_z[ci]
            if abs(cc[2] - z0) < r_nuc:
                rr = np.sqrt(max(r_nuc**2 - (cc[2] - z0) ** 2, 0.0)) * um
                nuc.append(Circle((cc[0] * um, cc[1] * um), rr))
        if nuc:
            ax.add_collection(PatchCollection(nuc, facecolor="#2c3e50", edgecolor="#e74c3c", alpha=0.55, lw=0.7))
        ax.set_aspect("equal")
        ax.set_xlim(com[0] * um - 1.1 * R, com[0] * um + 1.1 * R)
        ax.set_ylim(com[1] * um - 1.1 * R, com[1] * um + 1.1 * R)
        ax.set_title(f"z = {(z0 - com[2]) * um:+.0f} µm  ({len(nuc)} nuclei)", fontsize=9)
        ax.tick_params(labelsize=6)
    for s in range(n_slices, len(axes)):
        axes[s].axis("off")
    aa0 = float(d["aa0"][frame]) if "aa0" in d else float("nan")
    vv0 = float(d["vv0"][frame]) if "vv0" in d else float("nan")
    pen = float(d["pen_frac"][frame]) if "pen_frac" in d else float("nan")
    fig.suptitle(f"{os.path.basename(npz_path)} — {len(cells)}-cell aggregation, z-montage with nuclei\n"
                 f"A/A0={aa0:.3f}  V/V0={vv0:.3f}  pen={pen:.2f}  (clean round all-touching spheroid)",
                 fontsize=11, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    if out is None:
        out = os.path.join(os.path.dirname(npz_path), "figs",
                           os.path.splitext(os.path.basename(npz_path))[0] + "_nucleus_montage.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=125)
    print(f"saved {out}  ({len(cells)} cells, {n_slices} z-slices)")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--out", default=None)
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--rcell", type=float, default=7.5e-6)
    ap.add_argument("--rnuc-factor", type=float, default=0.33, dest="rnuc_factor")
    a = ap.parse_args()
    render(a.npz, a.out, a.frame, a.rcell, a.rnuc_factor)
