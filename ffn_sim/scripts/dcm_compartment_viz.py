"""Render the DCM full-compartment cell(s) so the compartments are VISIBLE:
plasma-membrane/cortex (outer shell) + cytoplasm (interior) + nucleus (inner sphere).

The nucleus in the Warp DCM is a centroid-anchored force (R_nuc = R_nuc_factor·R); it has no stored
mesh, so its position/size are reconstructed from each cell's node centroid + radius. Produces:
  (A) a 3D view: cortex shells (translucent) + nucleus spheres (opaque)
  (B) a mid-plane CROSS-SECTION: cortex outline (mesh-slice) + cytoplasm fill + nucleus disc
so 'why can't I see cytoplasm/plasma-membrane/nucleus' is answered directly.
"""
import argparse, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon


def cell_slice_outline(pos, faces_cell, z_mid):
    """xy points where the cell's triangle edges cross z=z_mid (the cortex cross-section outline)."""
    pts = []
    seen = set()
    for tri in faces_cell:
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            key = (min(a, b), max(a, b))
            if key in seen:
                continue
            seen.add(key)
            za, zb = pos[a, 2], pos[b, 2]
            if (za - z_mid) * (zb - z_mid) < 0:
                t = (z_mid - za) / (zb - za)
                pts.append(pos[a, :2] + t * (pos[b, :2] - pos[a, :2]))
    return np.array(pts) if pts else np.zeros((0, 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--r-nuc-factor", type=float, default=0.25)
    a = ap.parse_args()
    d = np.load(a.npz)
    fr, cof, faces = d["frames"], d["cof"], d["faces"]
    nc = int(cof.max()) + 1
    pos = fr[a.frame] * 1e6  # um
    faces_um = faces
    cens, radii = [], []
    for c in range(nc):
        p = pos[cof == c]
        cen = p.mean(0)
        cens.append(cen); radii.append(np.linalg.norm(p - cen, axis=1).mean())
    cens = np.array(cens); radii = np.array(radii)
    z_mid = float(np.median(cens[:, 2]))
    cmap = plt.cm.tab20(np.linspace(0, 1, max(nc, 2)))

    fig = plt.figure(figsize=(13, 5.6))
    # ---- Panel A: 3D cortex shells (translucent) + nucleus spheres ----
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    for c in range(nc):
        p = pos[cof == c]
        ax.scatter(p[:, 0], p[:, 1], p[:, 2], s=2, color=cmap[c], alpha=0.25)
        # nucleus sphere (centroid + R_nuc)
        rn = a.r_nuc_factor * radii[c]
        u, v = np.mgrid[0:2*np.pi:16j, 0:np.pi:8j]
        xs = cens[c, 0] + rn*np.cos(u)*np.sin(v); ys = cens[c, 1] + rn*np.sin(u)*np.sin(v); zs = cens[c, 2] + rn*np.cos(v)
        ax.plot_surface(xs, ys, zs, color=cmap[c], alpha=0.9, linewidth=0)
    ax.set_title("3D: cortex/membrane shells (translucent)\n+ NUCLEUS spheres (opaque, R_nuc=0.25R)")
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)"); ax.set_zlabel("z (um)")

    # ---- Panel B: mid-plane CROSS-SECTION (membrane / cytoplasm / nucleus) ----
    ax2 = fig.add_subplot(1, 2, 2)
    for c in range(nc):
        fc = faces_um[np.all(np.isin(faces_um, np.where(cof == c)[0]), axis=1)]
        out = cell_slice_outline(pos, fc, z_mid)
        if out.shape[0] >= 3:
            # order the outline points by angle about the centroid for a clean polygon
            ang = np.arctan2(out[:, 1]-cens[c, 1], out[:, 0]-cens[c, 0])
            out = out[np.argsort(ang)]
            ax2.add_patch(Polygon(out, closed=True, facecolor=cmap[c], alpha=0.30, edgecolor=cmap[c], lw=1.8))  # cytoplasm fill + membrane edge
        # nucleus disc at this plane
        rn = a.r_nuc_factor * radii[c]; dz = abs(cens[c, 2] - z_mid)
        if dz < rn:
            ax2.add_patch(Circle((cens[c, 0], cens[c, 1]), np.sqrt(rn**2 - dz**2), facecolor=cmap[c], alpha=0.95, edgecolor="k", lw=0.6))
    ax2.set_aspect("equal"); ax2.autoscale_view()
    lim = np.abs(pos[:, :2]).max() * 1.05
    ax2.set_xlim(pos[:, 0].min()-2, pos[:, 0].max()+2); ax2.set_ylim(pos[:, 1].min()-2, pos[:, 1].max()+2)
    ax2.set_title(f"cross-section @ z={z_mid:.1f}um: MEMBRANE (edge) +\nCYTOPLASM (fill) + NUCLEUS (inner disc)")
    ax2.set_xlabel("x (um)"); ax2.set_ylabel("y (um)")
    fig.suptitle(f"DCM full-compartment cell (N={nc}) — plasma membrane/cortex + cytoplasm + nucleus, GPU-native", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(a.out, dpi=120)
    print(f"saved {a.out}  (N={nc}, z_mid={z_mid:.2f}um, mean R={radii.mean():.2f}um, R_nuc={a.r_nuc_factor*radii.mean():.2f}um)")


if __name__ == "__main__":
    main()
