"""Mesh-surface visualization of the Warp de-cohesion spread (Phase-C migration).

Renders the cell SURFACE meshes (triangulated faces) — NOT point clouds — across frames,
in two views:
  * top-down (xy): the A/A0 silhouette (footprint) — what aggregation compacts / spread grows.
  * side (xz): the height profile — maxZ HELD while the footprint widens = genuine lateral
    spread (not the vertical-flatten settle-confound).

The npz (from ``dcm_warp_decohesion.py --save-frames``) spans BOTH phases when
``--settle-frames`` was used: phase 0 = AGGREGATE (the spherical cluster compacting into a
rounded spheroid), phase 1 = SPREAD. ``aa0`` is normalised to the rested baseline, so
aggregate frames read >1 falling to 1.0 (compaction) and spread frames climb above 1.0.

Color modes:
  * default       — per cell (tab20), distinguishes cells.
  * --color-by-contact  — rim (ECM-contacting) red, dragged (no ECM) blue: the collective drag.
  * --color-by-junction — per-frame from ``cad``: junction-SWITCHED cells (cad<1, cadherin
    weakened → integrin strengthened) orange, intact cells teal: shows the switch propagating.

Integrity (viz rules + landmine lessons): mesh not points; equal aspect; axis limits FIXED
across frames (never per-frame autoscale); A/A0 = top-down silhouette (never basal contact);
SI→µm; NO body-force proxy.

    python -m ffn_sim.scripts.warp_decoh_mesh_viz --npz frames.npz --out fig.png \
        [--mp4 surf.mp4] [--color-by-contact | --color-by-junction] [--ncols 6] [--fps 8]
"""

from __future__ import annotations

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch

UM = 1.0e6
ALPHA = 0.6
RIM_RGBA = np.array([0.85, 0.20, 0.15, ALPHA])      # rim / ECM-contacting
DRAG_RGBA = np.array([0.15, 0.40, 0.85, ALPHA])     # dragged / non-ECM
SWITCH_RGBA = np.array([0.95, 0.55, 0.10, ALPHA])   # junction-switched (cadherin weak)
INTACT_RGBA = np.array([0.10, 0.65, 0.55, ALPHA])   # intact junction


def _rim_mask_cells(frame0, cof, z0=0.0, contact_band=1.5):
    """Per-cell bool: rim (ECM-contacting) at the baseline by the driver's criterion
    (centroid z within ``contact_band·R`` of z0; R estimated as the mean per-cell node radius)."""
    n_cells = int(cof.max()) + 1
    valid = cof >= 0                              # exclude dormant/parked nodes (cof=-1) — else
    vc = cof[valid]; vp = frame0[valid]           # np.add.at folds them into cell[-1] (review fix #3)
    cen = np.zeros((n_cells, 3)); cnt = np.zeros(n_cells)
    np.add.at(cen, vc, vp); np.add.at(cnt, vc, 1.0)
    cen /= np.maximum(cnt, 1.0)[:, None]
    rad = np.linalg.norm(vp - cen[vc], axis=1)
    R = np.zeros(n_cells); rc = np.zeros(n_cells)
    np.add.at(R, vc, rad); np.add.at(rc, vc, 1.0)
    return (cen[:, 2] - z0) <= contact_band * float((R / np.maximum(rc, 1.0)).mean())


def _make_color_fn(faces_at, cof_at, frames, mode, cad):
    """Return (colors_for(fi), legend_handles). ``colors_for(fi)`` gives the (M_fi,4) per-face
    RGBA for frame fi. Uses per-frame topology (``faces_at``/``cof_at``) so it is correct even
    when remesh changed the mesh between frames."""
    if mode == "contact":
        rim = _rim_mask_cells(frames[-1], cof_at(len(frames) - 1))   # rim at the spread baseline
        def colors_for(fi):
            fcell = cof_at(fi)[faces_at(fi)[:, 0]]
            return np.where(rim[fcell][:, None], RIM_RGBA, DRAG_RGBA)
        leg = [Patch(color=RIM_RGBA[:3], label="rim (ECM-contacting, crawls)"),
               Patch(color=DRAG_RGBA[:3], label="dragged (no ECM contact)")]
        print(f"rim cells: {int(rim.sum())}/{rim.size}  dragged: {int((~rim).sum())}")
        return colors_for, leg
    if mode == "junction":
        if cad is None:
            raise SystemExit("--color-by-junction needs a `cad` array in the npz (run with --junction-switch)")
        def colors_for(fi):
            fcell = cof_at(fi)[faces_at(fi)[:, 0]]
            switched = cad[fi][fcell] < 0.999
            return np.where(switched[:, None], SWITCH_RGBA, INTACT_RGBA)
        leg = [Patch(color=SWITCH_RGBA[:3], label="junction-switched (cadherin↓ integrin↑)"),
               Patch(color=INTACT_RGBA[:3], label="intact junction")]
        return colors_for, leg
    def colors_for(fi):
        fcell = cof_at(fi)[faces_at(fi)[:, 0]]
        fc = plt.get_cmap("tab20")(fcell % 20); fc[:, 3] = ALPHA
        return fc
    return colors_for, None


def _phase_tag(phase, fi):
    return "AGG" if (phase is not None and phase[fi] == 0) else "SPREAD"


def _write_mp4(frames, faces_at, colors_for, step, phase, aa0, maxZ, vv0, pen, out, fps, legend, title):
    import matplotlib.animation as animation
    allp = frames * UM
    xylim = float(np.abs(allp[..., :2]).max()) * 1.05
    zmin, zmax = float(allp[..., 2].min()), float(allp[..., 2].max())

    def tri(P, i, j, fi):
        return P[faces_at(fi)][:, :, [i, j]]

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 5.4))
    pc0 = PolyCollection([], edgecolors=(0, 0, 0, 0.12), linewidths=0.1)
    pc1 = PolyCollection([], edgecolors=(0, 0, 0, 0.12), linewidths=0.1)
    ax0.add_collection(pc0); ax1.add_collection(pc1)
    ax0.set_xlim(-xylim, xylim); ax0.set_ylim(-xylim, xylim); ax0.set_aspect("equal")
    ax0.set_xlabel("x (µm)"); ax0.set_ylabel("y (µm)")
    ax1.set_xlim(-xylim, xylim); ax1.set_ylim(zmin - 2.0, zmax + 4.0); ax1.set_aspect("equal")
    ax1.axhline(0.0, color="saddlebrown", lw=1.2, alpha=0.7)
    ax1.set_xlabel("x (µm)"); ax1.set_ylabel("z (µm)")
    if legend:
        ax0.legend(handles=legend, loc="upper right", fontsize=7, framealpha=0.9)
    sup = fig.suptitle("")

    def update(fi):
        P = frames[fi] * UM
        col = colors_for(fi)
        pc0.set_verts(tri(P, 0, 1, fi)); pc0.set_facecolors(col)
        pc1.set_verts(tri(P, 0, 2, fi)); pc1.set_facecolors(col)
        ax0.set_title(f"[{_phase_tag(phase, fi)}] top-down — A/A0 = {aa0[fi]:.2f}", fontsize=10)
        ax1.set_title(f"side — maxZ = {maxZ[fi]:.0f}µm  V/V0 = {vv0[fi]:.2f}  pen = {pen[fi]:.2f}", fontsize=10)
        sup.set_text(f"{title}  ·  step {int(step[fi])}")
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
    ap.add_argument("--ncols", type=int, default=6)
    ap.add_argument("--mp4", default=None, help="also render a surface animation to this mp4 path")
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--color-by-contact", action="store_true",
                    help="color rim (ECM-contacting) cells red, dragged non-ECM cells blue")
    ap.add_argument("--color-by-junction", action="store_true",
                    help="color junction-switched cells (cad<1) orange, intact teal (per frame)")
    ap.add_argument("--title", default="Warp DCM · de-cohesion (M1-M3)")
    args = ap.parse_args()

    d = np.load(args.npz, allow_pickle=True)
    frames = d["frames"]; faces = d["faces"]; cof = d["cof"]
    step, aa0, maxZ, vv0 = d["step"], d["aa0"], d["maxZ"], d["vv0"]
    phase = d["phase"] if "phase" in d.files else None
    cad = d["cad"] if "cad" in d.files else None
    pen = d["pen_frac"] if "pen_frac" in d.files else np.zeros(frames.shape[0])
    # per-frame topology (remesh): faces/cof vary across frames → object sequences
    faces_seq = d["faces_seq"] if "faces_seq" in d.files else None
    cof_seq = d["cof_seq"] if "cof_seq" in d.files else None
    # object sequences store per-frame int arrays; cast to int (object dtype can't index)
    faces_at = (lambda fi: np.asarray(faces_seq[fi], dtype=np.int64)) if faces_seq is not None else (lambda fi: faces)
    cof_at = (lambda fi: np.asarray(cof_seq[fi], dtype=np.int64)) if cof_seq is not None else (lambda fi: cof)
    F = frames.shape[0]
    print(f"interpenetration max/mean_edge: peak {float(pen.max()):.3f}  final {float(pen[-1]):.3f} "
          f"(0 = no mesh overlap; >~0.3 = cells interpenetrating)"
          + ("  [per-frame topology: remesh ON]" if faces_seq is not None else ""))

    mode = "junction" if args.color_by_junction else ("contact" if args.color_by_contact else "default")
    colors_for, legend = _make_color_fn(faces_at, cof_at, frames, mode, cad)

    if args.mp4:
        _write_mp4(frames, faces_at, colors_for, step, phase, aa0, maxZ, vv0, pen,
                   args.mp4, args.fps, legend, args.title)

    sel = np.unique(np.linspace(0, F - 1, args.ncols).astype(int))
    allp = frames * UM
    xylim = float(np.abs(allp[..., :2]).max()) * 1.05
    zmin, zmax = float(allp[..., 2].min()), float(allp[..., 2].max())

    def tri(P, i, j, fi):
        return P[faces_at(fi)][:, :, [i, j]]

    fig, axes = plt.subplots(2, len(sel), figsize=(3.2 * len(sel), 6.8))
    if len(sel) == 1:
        axes = axes.reshape(2, 1)
    for c, fi in enumerate(sel):
        P = frames[fi] * UM
        col = colors_for(fi)
        ax = axes[0, c]
        ax.add_collection(PolyCollection(tri(P, 0, 1, fi), facecolors=col,
                                         edgecolors=(0, 0, 0, 0.12), linewidths=0.1))
        ax.set_xlim(-xylim, xylim); ax.set_ylim(-xylim, xylim); ax.set_aspect("equal")
        ax.set_title(f"[{_phase_tag(phase, fi)}] step {int(step[fi])}\nA/A0 = {aa0[fi]:.2f}", fontsize=9)
        ax.tick_params(labelsize=6)
        if c == 0:
            ax.set_ylabel("top-down (xy)  µm\n[footprint / A0_rested]", fontsize=8)
            if legend:
                ax.legend(handles=legend, loc="upper left", fontsize=6, framealpha=0.9)
        ax = axes[1, c]
        ax.add_collection(PolyCollection(tri(P, 0, 2, fi), facecolors=col,
                                         edgecolors=(0, 0, 0, 0.12), linewidths=0.1))
        ax.set_xlim(-xylim, xylim); ax.set_ylim(zmin - 2.0, zmax + 4.0); ax.set_aspect("equal")
        ax.axhline(0.0, color="saddlebrown", lw=1.2, alpha=0.7)
        ax.set_title(f"maxZ = {maxZ[fi]:.0f} µm\nV/V0 = {vv0[fi]:.2f}  pen = {pen[fi]:.2f}", fontsize=9)
        ax.tick_params(labelsize=6)
        if c == 0:
            ax.set_ylabel("side (xz)  µm\n[maxZ held]", fontsize=8)

    fig.suptitle(
        args.title + " — cell SURFACE mesh, NO body-force proxy, physiological γ, dt=8e-6\n"
        "AGG phase: spherical cluster COMPACTS (A/A0→1) · SPREAD phase: footprint GROWS at held maxZ",
        fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(args.out, dpi=130)
    print(f"wrote {args.out}  ({len(sel)} frames, {faces.shape[0]} tri/frame, "
          f"A/A0 {aa0[0]:.2f}->{aa0[-1]:.2f}, maxZ {maxZ[0]:.0f}->{maxZ[-1]:.0f}µm)")


if __name__ == "__main__":
    main()
