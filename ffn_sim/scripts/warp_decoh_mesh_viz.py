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


def _rim_mask_cells(frame0: np.ndarray, cof: np.ndarray, z0: float = 0.0,
                    contact_band: float = 1.5) -> np.ndarray:
    """Which CELLS are rim (ECM-contacting) at the rested baseline (frame 0), by the same
    criterion the driver uses: centroid z within ``contact_band·R`` of z0. R is estimated
    from the mesh (mean per-cell node radius). Returns a per-cell bool array."""
    n_cells = int(cof.max()) + 1
    cen = np.zeros((n_cells, 3)); cnt = np.zeros(n_cells)
    np.add.at(cen, cof, frame0); np.add.at(cnt, cof, 1.0)
    cen /= np.maximum(cnt, 1.0)[:, None]
    # per-cell radius = mean node distance from its own centroid
    rad = np.linalg.norm(frame0 - cen[cof], axis=1)
    R = np.zeros(n_cells); rc = np.zeros(n_cells)
    np.add.at(R, cof, rad); np.add.at(rc, cof, 1.0)
    R /= np.maximum(rc, 1.0)
    Rmean = float(R.mean())
    return (cen[:, 2] - z0) <= contact_band * Rmean


def _face_colors(faces, cof, color_by_contact, frame0, alpha=0.55):
    """Per-face RGBA. Default: tab20 per cell. ``color_by_contact``: rim (ECM-contacting)
    cells warm/red, dragged (non-ECM) cells cool/blue — so the collective drag reads."""
    face_cell = cof[faces[:, 0]]
    if color_by_contact:
        rim = _rim_mask_cells(frame0, cof)
        fc = np.where(rim[face_cell][:, None],
                      np.array([0.85, 0.20, 0.15, alpha]),   # rim — crawl front (red)
                      np.array([0.15, 0.40, 0.85, alpha]))    # dragged — non-ECM (blue)
        return fc, rim
    fc = plt.get_cmap("tab20")(face_cell % 20); fc[:, 3] = alpha
    return fc, None


def _write_mp4(frames, faces, fcolors, step, aa0, maxZ, vv0, out, fps, rim_legend):
    """Animate the two-view surface render across ALL frames → mp4 (ffmpeg) or gif fallback."""
    import matplotlib.animation as animation
    allp = frames * UM
    xylim = float(np.abs(allp[..., :2]).max()) * 1.05
    zmin, zmax = float(allp[..., 2].min()), float(allp[..., 2].max())

    def tri(P, i, j):
        return P[faces][:, :, [i, j]]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 5.4))
    pc0 = PolyCollection([], facecolors=fcolors, edgecolors=(0, 0, 0, 0.12), linewidths=0.1)
    pc1 = PolyCollection([], facecolors=fcolors, edgecolors=(0, 0, 0, 0.12), linewidths=0.1)
    ax0.add_collection(pc0); ax1.add_collection(pc1)
    ax0.set_xlim(-xylim, xylim); ax0.set_ylim(-xylim, xylim); ax0.set_aspect("equal")
    ax0.set_xlabel("x (µm)"); ax0.set_ylabel("y (µm)")
    ax1.set_xlim(-xylim, xylim); ax1.set_ylim(zmin - 2.0, zmax + 4.0); ax1.set_aspect("equal")
    ax1.axhline(0.0, color="saddlebrown", lw=1.2, alpha=0.7)
    ax1.set_xlabel("x (µm)"); ax1.set_ylabel("z (µm)")
    if rim_legend:
        from matplotlib.patches import Patch
        ax0.legend(handles=[Patch(color=(0.85, 0.20, 0.15), label="rim (ECM-contacting, crawls)"),
                            Patch(color=(0.15, 0.40, 0.85), label="dragged (no ECM contact)")],
                   loc="upper right", fontsize=7, framealpha=0.9)
    sup = fig.suptitle("")

    def update(fi):
        P = frames[fi] * UM
        pc0.set_verts(tri(P, 0, 1)); pc1.set_verts(tri(P, 0, 2))
        ax0.set_title(f"top-down (xy) — A/A0 = {aa0[fi]:.2f}", fontsize=10)
        ax1.set_title(f"side (xz) — maxZ = {maxZ[fi]:.0f}µm  V/V0 = {vv0[fi]:.2f}", fontsize=10)
        sup.set_text(f"Warp DCM · de-cohesion N=100 lamellipodium+wetting (M2)  ·  step {int(step[fi])}")
        return pc0, pc1

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=False)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    try:
        anim.save(out, writer=animation.FFMpegWriter(fps=fps, bitrate=2400))
        print(f"wrote {out}  ({len(frames)} frames @ {fps}fps)")
    except Exception as e:
        gif = out.rsplit(".", 1)[0] + ".gif"
        anim.save(gif, writer=animation.PillowWriter(fps=fps))
        print(f"ffmpeg unavailable ({e}); wrote {gif} instead")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ncols", type=int, default=5)
    ap.add_argument("--mp4", default=None, help="also render a surface animation to this mp4 path")
    ap.add_argument("--fps", type=int, default=6)
    ap.add_argument("--color-by-contact", action="store_true",
                    help="color rim (ECM-contacting) cells red, dragged non-ECM cells blue")
    args = ap.parse_args()

    d = np.load(args.npz)
    frames = d["frames"]                 # (F, N, 3) metres
    faces = d["faces"]                   # (M, 3) node indices
    cof = d["cof"]                       # (N,) cell-of-node
    step, aa0, maxZ, vv0 = d["step"], d["aa0"], d["maxZ"], d["vv0"]
    F = frames.shape[0]

    sel = np.unique(np.linspace(0, F - 1, args.ncols).astype(int))
    fcolors, rim = _face_colors(faces, cof, args.color_by_contact, frames[0])
    if rim is not None:
        print(f"rim (ECM-contacting) cells: {int(rim.sum())}/{rim.size}  "
              f"dragged (non-ECM): {int((~rim).sum())}")
    if args.mp4:
        _write_mp4(frames, faces, fcolors, step, aa0, maxZ, vv0, args.mp4, args.fps,
                   rim_legend=args.color_by_contact)

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

    _drag = " · rim=red (ECM-contacting), dragged=blue (no ECM)" if rim is not None else ""
    fig.suptitle(
        "Warp DCM engine — de-cohesion N=100 lamellipodium+wetting cleanball (M2, gbook A5000)\n"
        "cell SURFACE mesh" + _drag + " · NO body-force proxy, physiological γ, dt=8e-6\n"
        "top-down (A/A0 silhouette) GROWS while side (maxZ) HELD → genuine lateral spread; "
        "bonded rim cells drag the non-ECM cells outward (collective)",
        fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(args.out, dpi=130)
    print(f"wrote {args.out}  ({len(sel)} frames, {faces.shape[0]} triangles/frame, "
          f"A/A0 {aa0[0]:.2f}->{aa0[-1]:.2f}, maxZ {maxZ[0]:.0f}->{maxZ[-1]:.0f}µm)")


if __name__ == "__main__":
    main()
