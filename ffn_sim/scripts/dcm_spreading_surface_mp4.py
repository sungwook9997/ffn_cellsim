"""Render the GPU-friendly DCM spheroid SPREADING as real CELL SURFACES.

WHY THIS SCRIPT (PI 2026-06-11, rightly angry). The prior "spreading"
visualisations were DOT / disc scatter and — worse — the spheroid did NOT
actually attach + spread, it FLOATED. Two fixes landed in the build:
``cell/dcm_gpu_build.build_gpu_dcm_simulation`` now adds a SETTLING body force
(default ``settle_force=4e-10``) so the spheroid SEDIMENTS onto the dish,
CONTACTS, FLATTENS and SPREADS (genuine wetting). This script VISUALISES that as
real triangulated cell SURFACES (à la ``scripts/dcm_surface_mp4.py`` /
``spheroid_confluent_surface.mp4``), NOT dots — so the PI can SEE a 3D ball land
on the substrate, wet it, and spread out over (real-time-equivalent) time.

WHAT IT SHOWS over the frames (the whole point):
    3D ball above/at the dish  →  descends (settling)  →  contacts the floor  →
    flattens  →  spreads its footprint.

HOW. We build the FIXED spheroid (settling ON by default), run a finite gate
(run(0)→run(2000)), then run in CHUNKS capturing tag-ordered node positions per
frame. Every cell's triangulated icosphere shell (the SAME ``tris`` per cell from
``dcm.icosphere_mesh``, offset into each cell's node range) is drawn as a shaded
``Poly3DCollection``. A translucent grey SUBSTRATE PLANE at z=z0 makes the
contact + spreading obvious. Panels: a SIDE (x-z) view (primary — shows descend /
contact / flatten / spread), a 3D perspective, and a top-down (x-y) footprint.
Each frame's title carries step, real-time-equivalent (``common/sim_realtime``),
substrate-contact %, and footprint width.

Run FROM REPO ROOT (CPU on this Mac is fine):
    conda activate ffn_sim
    PYTHONPATH=. python ffn_sim/scripts/dcm_spreading_surface_mp4.py

Outputs:
    ffn_sim/outputs/h_dcm_gpu_lod/figs/spreading_surface.mp4
    ffn_sim/outputs/h_dcm_gpu_lod/figs/spreading_surface_final.png

This script ONLY CALLS the read-only DCM build/integrator modules; it does not
modify any of them.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FFMpegWriter  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from matplotlib.colors import LightSource  # noqa: E402

from ffn_sim.archive.hoomd_legacy.cell.dcm import icosphere_mesh  # noqa: E402
from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import (  # noqa: E402
    ResolvedGpuDCM, build_gpu_dcm_simulation)
from ffn_sim.archive.hoomd_legacy.cell.dcm_confluence import capture_positions  # noqa: E402
from ffn_sim.common.sim_realtime import map_realtime  # noqa: E402

matplotlib.rcParams["animation.ffmpeg_path"] = "/opt/homebrew/bin/ffmpeg"

OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod/figs")
OUT.mkdir(parents=True, exist_ok=True)
MP4 = OUT / "spreading_surface.mp4"
PNG = OUT / "spreading_surface_final.png"

UM = 1.0e6
_TAB = plt.get_cmap("tab20").colors


# ---------------------------------------------------------------------------
# Surface shading / drawing (copied from scripts/dcm_surface_mp4.py — the
# template the PI named; one shaded Poly3DCollection per cell, NOT dots).
# ---------------------------------------------------------------------------
def cell_colour(c: int):
    return _TAB[c % len(_TAB)]


def shade_faces(verts, tris, base_rgb, ls: LightSource, alpha: float = 1.0):
    """Per-face RGBA with simple Lambert shading from a LightSource."""
    tri = verts[tris]                                   # (k,3,3)
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.where(norm > 0, norm, 1.0)
    az = np.deg2rad(ls.azdeg)
    alt = np.deg2rad(ls.altdeg)
    light = np.array([np.cos(alt) * np.cos(az),
                      np.cos(alt) * np.sin(az), np.sin(alt)])
    inten = np.clip(np.abs(n @ light), 0.0, 1.0)
    shade = 0.45 + 0.55 * inten                          # ambient + diffuse
    base = np.array(base_rgb)
    rgba = np.empty((tris.shape[0], 4))
    rgba[:, :3] = np.clip(base[None, :] * shade[:, None], 0, 1)
    rgba[:, 3] = alpha
    return rgba


def draw_surfaces(ax, pos, ranges, tris0, ls, *, alpha=1.0):
    """Draw each cell's triangulated shell as a shaded Poly3DCollection.

    ``tris0`` is the single per-cell icosphere triangle list (LOCAL indices into
    a cell's node block) — every cell shares it, offset by its ``ranges`` block.
    """
    for c, (lo, hi) in enumerate(ranges):
        verts = pos[lo:hi]
        rgba = shade_faces(verts, tris0, cell_colour(c), ls, alpha=alpha)
        polys = verts[tris0]
        coll = Poly3DCollection(polys, facecolors=rgba,
                                edgecolors=(0, 0, 0, 0.10), linewidths=0.12)
        ax.add_collection3d(coll)


def draw_substrate(ax, *, xlim, ylim, z0_um, colour="0.55", alpha=0.30):
    """Draw a translucent grey floor plane at z=z0 so contact/spread is obvious."""
    xs = np.array([[xlim[0], xlim[1]], [xlim[0], xlim[1]]])
    ys = np.array([[ylim[0], ylim[0]], [ylim[1], ylim[1]]])
    zs = np.full_like(xs, z0_um)
    ax.plot_surface(xs, ys, zs, color=colour, alpha=alpha,
                    shade=False, zorder=0)


# ---------------------------------------------------------------------------
# Per-frame physical metrics (maxZ drop / substrate-contact % / footprint width)
# ---------------------------------------------------------------------------
def frame_metrics(pos, *, z0: float, R: float) -> dict:
    """maxZ [µm], substrate-contact fraction [%], footprint width [µm], area [µm²].

    contact = fraction of nodes within the basal contact band (z < z0 + 0.6·R).
    footprint width = lateral (x) extent of basal-band nodes; area = xy
    convex-hull area of basal-contact nodes (same proxy as
    ``dcm_production_spread_mp4._footprint_area``).
    """
    from scipy.spatial import ConvexHull

    maxZ = float(pos[:, 2].max()) * UM
    band_mask = pos[:, 2] < z0 + 0.6 * R
    contact = float(band_mask.mean()) * 100.0
    basal = pos[pos[:, 2] < z0 + 0.5 * R][:, :2]
    if basal.shape[0] >= 3:
        width = float(np.ptp(basal[:, 0])) * UM
        try:
            area = float(ConvexHull(basal).volume) * UM * UM   # 2D hull area
        except Exception:  # noqa: BLE001 — degenerate (collinear) basal set
            area = float(np.ptp(basal[:, 0]) * np.ptp(basal[:, 1])) * UM * UM
    else:
        width = 0.0
        area = 0.0
    return dict(maxZ=maxZ, contact=contact, width=width, area=area)


# ---------------------------------------------------------------------------
# Fixed view box (from the union of all frames so the camera never jumps)
# ---------------------------------------------------------------------------
def view_box(frames, *, z0_um: float, pad: float = 1.08):
    allpos = np.concatenate(frames, axis=0) * UM
    cx = 0.5 * (allpos[:, 0].min() + allpos[:, 0].max())
    cy = 0.5 * (allpos[:, 1].min() + allpos[:, 1].max())
    rx = 0.5 * np.ptp(allpos[:, 0]) * pad
    ry = 0.5 * np.ptp(allpos[:, 1]) * pad
    r = max(rx, ry)
    xlim = (cx - r, cx + r)
    ylim = (cy - r, cy + r)
    zlo = min(z0_um, allpos[:, 2].min()) - 2.0
    zhi = allpos[:, 2].max() * pad
    return xlim, ylim, (zlo, zhi)


def set_lims(ax, xlim, ylim, zlim):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_zlim(*zlim)
    ax.set_box_aspect((xlim[1] - xlim[0], ylim[1] - ylim[0], zlim[1] - zlim[0]))


# ---------------------------------------------------------------------------
# Build + capture
# ---------------------------------------------------------------------------
def build_and_capture(*, n_cells: int, n_frames: int, steps_per_frame: int,
                      gate_steps: int):
    """Build the FIXED (settling-ON) spheroid and capture per-frame node positions.

    Returns ``(frames, metrics, ranges, tris0, p, steps)`` where ``frames`` is a
    list of (N,3) tag-ordered positions and ``metrics`` the per-frame physical
    metrics; ``steps`` is the integrator step at each captured frame.
    """
    p = ResolvedGpuDCM(R_cell=7.5e-6, subdivisions=1)
    # CPU device (this Mac is CPU-only); pick_device(None) falls back to CPU.
    h = build_gpu_dcm_simulation(p, n_cells, active=True)
    sim = h["sim"]
    ranges = h["ranges"]
    n_built = h["n_cells"]
    _verts0, _edges, tris0 = icosphere_mesh(p.R_cell, p.subdivisions)

    z0 = float(p.z_substrate)
    R = float(p.R_cell)

    # finite gate
    sim.run(0)
    t0 = time.time()
    sim.run(gate_steps)
    print(f"[gate] run(0)->run({gate_steps}) in {time.time()-t0:.0f}s; "
          f"n_cells={n_built}, N_nodes={sum(hi-lo for lo, hi in ranges)}",
          flush=True)

    frames, metrics, steps = [], [], []
    pos0 = capture_positions(sim)
    if not np.isfinite(pos0).all():
        raise RuntimeError("non-finite positions after the finite gate")
    frames.append(pos0)
    metrics.append(frame_metrics(pos0, z0=z0, R=R))
    steps.append(gate_steps)
    m0 = metrics[0]
    print(f"  frame  0  step {gate_steps:>6d}  maxZ {m0['maxZ']:5.1f}µm  "
          f"contact {m0['contact']:4.1f}%  width {m0['width']:5.1f}µm",
          flush=True)

    for f in range(1, n_frames):
        t0 = time.time()
        sim.run(steps_per_frame)
        pos = capture_positions(sim)
        if not np.isfinite(pos).all():
            print(f"[warn] non-finite at frame {f}; truncating", flush=True)
            break
        frames.append(pos)
        metrics.append(frame_metrics(pos, z0=z0, R=R))
        steps.append(gate_steps + f * steps_per_frame)
        m = metrics[-1]
        print(f"  frame {f:>2d}  step {steps[-1]:>6d}  maxZ {m['maxZ']:5.1f}µm  "
              f"contact {m['contact']:4.1f}%  width {m['width']:5.1f}µm  "
              f"({time.time()-t0:.0f}s)", flush=True)

    return frames, metrics, ranges, tris0, p, steps


# ---------------------------------------------------------------------------
# Render MP4
# ---------------------------------------------------------------------------
def render(frames, metrics, ranges, tris0, p, steps, *, fps: int):
    z0_um = float(p.z_substrate) * UM
    xlim, ylim, zlim = view_box(frames, z0_um=z0_um)
    ls = LightSource(azdeg=225, altdeg=45)
    n = len(frames)
    accel = map_realtime(1, p.dt).accel_factor

    fig = plt.figure(figsize=(15.5, 6.2))
    axXZ = fig.add_subplot(131, projection="3d")   # primary SIDE view
    ax3D = fig.add_subplot(132, projection="3d")   # 3D perspective
    axXY = fig.add_subplot(133, projection="3d")   # top-down footprint

    writer = FFMpegWriter(fps=fps, metadata=dict(artist="ffn_cellsim"),
                          bitrate=3600)

    def stage(i):
        # narrative stage from the actual descent fraction
        return ("ball above the dish" if i < n * 0.18 else
                "descending / contacting" if i < n * 0.5 else
                "flattening" if i < n * 0.75 else
                "spreading on the dish")

    with writer.saving(fig, str(MP4), dpi=110):
        for i, pos in enumerate(frames):
            for ax in (axXZ, ax3D, axXY):
                ax.clear()
            pm = pos * UM
            m = metrics[i]
            rt = map_realtime(steps[i], p.dt, accel=accel)

            # SIDE (x-z): elev 0, azim -90 → look along +y, x horizontal, z up.
            draw_surfaces(axXZ, pm, ranges, tris0, ls)
            draw_substrate(axXZ, xlim=xlim, ylim=ylim, z0_um=z0_um)
            set_lims(axXZ, xlim, ylim, zlim)
            axXZ.view_init(elev=2, azim=-90)
            axXZ.set_axis_off()
            axXZ.set_title("SIDE (x–z): descend · contact · flatten · spread",
                           fontsize=10)

            # 3D perspective (gentle rotation)
            draw_surfaces(ax3D, pm, ranges, tris0, ls)
            draw_substrate(ax3D, xlim=xlim, ylim=ylim, z0_um=z0_um)
            set_lims(ax3D, xlim, ylim, zlim)
            ax3D.view_init(elev=18, azim=(-60 + 50.0 * i / max(1, n)))
            ax3D.set_axis_off()
            ax3D.set_title("3D perspective", fontsize=10)

            # top-down (x-y): footprint
            draw_surfaces(axXY, pm, ranges, tris0, ls)
            draw_substrate(axXY, xlim=xlim, ylim=ylim, z0_um=z0_um, alpha=0.18)
            set_lims(axXY, xlim, ylim, zlim)
            axXY.view_init(elev=89, azim=-90)
            axXY.set_axis_off()
            axXY.set_title("top-down (x–y): footprint", fontsize=10)

            fig.suptitle(
                f"DCM spheroid WETTING + SPREADING on the substrate — "
                f"{len(ranges)} deformable cells (surfaces)   ·   {stage(i)}\n"
                f"step {steps[i]:>6d}   "
                f"t = {rt.t_sim_human} sim ≈ {rt.t_real_human} real "
                f"(accel {rt.accel_factor:.0e})   ·   "
                f"maxZ = {m['maxZ']:.1f} µm   "
                f"substrate-contact = {m['contact']:.0f}%   "
                f"footprint = {m['width']:.1f} µm",
                fontsize=11)
            fig.tight_layout(rect=(0, 0, 1, 0.92))
            writer.grab_frame()
            if i % 5 == 0:
                print(f"  rendered frame {i}/{n}", flush=True)
    plt.close(fig)
    print(f"wrote {MP4}", flush=True)

    # final-frame PNG (the same three panels, fresh figure)
    figP = plt.figure(figsize=(15.5, 6.2))
    pm = frames[-1] * UM
    m = metrics[-1]
    rt = map_realtime(steps[-1], p.dt, accel=accel)
    for k, (elev, azim, ttl) in enumerate((
            (2, -90, "SIDE (x–z): flattened + spread on the dish"),
            (18, -35, "3D perspective"),
            (89, -90, "top-down (x–y): footprint"))):
        ax = figP.add_subplot(1, 3, k + 1, projection="3d")
        draw_surfaces(ax, pm, ranges, tris0, ls)
        draw_substrate(ax, xlim=xlim, ylim=ylim, z0_um=z0_um,
                       alpha=(0.18 if k == 2 else 0.30))
        set_lims(ax, xlim, ylim, zlim)
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(ttl, fontsize=10)
    figP.suptitle(
        f"DCM spheroid wetting + spreading — FINAL frame   "
        f"({len(ranges)} cells)\n"
        f"step {steps[-1]}   t = {rt.t_sim_human} sim ≈ {rt.t_real_human} real   ·   "
        f"maxZ = {m['maxZ']:.1f} µm   contact = {m['contact']:.0f}%   "
        f"footprint = {m['width']:.1f} µm", fontsize=11)
    figP.tight_layout(rect=(0, 0, 1, 0.90))
    figP.savefig(PNG, dpi=140)
    plt.close(figP)
    print(f"wrote {PNG}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=30,
                    help="cells in the spheroid (CPU surface render; 30 is a clear "
                         "spheroid that renders feasibly on CPU). Drop if too slow.")
    ap.add_argument("--frames", type=int, default=30,
                    help="captured frames (incl. the gate frame).")
    ap.add_argument("--steps-per-frame", type=int, default=600,
                    help="integrator steps between captured frames.")
    ap.add_argument("--gate-steps", type=int, default=2000,
                    help="finite settle-gate steps before frame 0.")
    ap.add_argument("--fps", type=int, default=6)
    args = ap.parse_args()

    t_all = time.time()
    frames, metrics, ranges, tris0, p, steps = build_and_capture(
        n_cells=args.n_cells, n_frames=args.frames,
        steps_per_frame=args.steps_per_frame, gate_steps=args.gate_steps)

    if len(frames) < 2:
        raise RuntimeError("too few finite frames captured to render")

    render(frames, metrics, ranges, tris0, p, steps, fps=args.fps)

    # spreading-reality summary (the validation the PI asked for)
    m0, m1 = metrics[0], metrics[-1]
    rt0 = map_realtime(steps[0], p.dt)
    rt1 = map_realtime(steps[-1], p.dt)
    print("\n=== SPREADING REALITY (first → last captured frame) ===", flush=True)
    print(f"  frames captured : {len(frames)}  "
          f"(steps {steps[0]} → {steps[-1]})", flush=True)
    print(f"  maxZ            : {m0['maxZ']:.1f} → {m1['maxZ']:.1f} µm  "
          f"(Δ {m1['maxZ']-m0['maxZ']:+.1f})", flush=True)
    print(f"  substrate-contact: {m0['contact']:.0f}% → {m1['contact']:.0f}%  "
          f"(Δ {m1['contact']-m0['contact']:+.0f})", flush=True)
    print(f"  footprint width : {m0['width']:.1f} → {m1['width']:.1f} µm  "
          f"(Δ {m1['width']-m0['width']:+.1f})", flush=True)
    print(f"  footprint area  : {m0['area']:.0f} → {m1['area']:.0f} µm²", flush=True)
    print(f"  real-time span  : {rt0.t_real_human} → {rt1.t_real_human} "
          f"(accel {rt1.accel_factor:.0e})", flush=True)
    print(f"  total wall time : {time.time()-t_all:.0f}s", flush=True)


if __name__ == "__main__":
    main()
