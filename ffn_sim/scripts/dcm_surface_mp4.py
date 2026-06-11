"""PART B — render the native-mesh DCM spheroid as ADHERED CELL SURFACES.

The emphasised deliverable: visualise the actual deformable cell SURFACES (each
cell's triangulated shell drawn as a lit 3D ``Poly3DCollection``), NOT a stylised
node/centroid scatter — so the PI can SEE the cells sticking together into a
confluent, space-filling tissue (à la SimuCell3D figures).

What it shows over time (the aggregation the PI wants to watch):
  separated FCC ball  →  cells adhere + deform  →  confluent spheroid.

We re-run a spheroid at the CONFLUENT params (Part A), capture per-frame
tag-ordered mesh-node positions, and for every frame rebuild each cell's surface
from its own triangulation (``per_cell_tris``) as a shaded Poly3DCollection in a
distinct colour, with a slow camera rotation. Output is MP4 (ffmpeg
``FFMpegWriter``) plus a still PNG (full view + an interior cross-section so the
internal packing is visible too).

Run FROM REPO ROOT:
    conda activate ffn_sim
    python scripts/dcm_surface_mp4.py

Outputs:
    ffn_sim/outputs/h_dcm_active/figs/spheroid_confluent_surface.mp4
    ffn_sim/outputs/h_dcm_active/figs/spheroid_confluent_final.png
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FFMpegWriter  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from matplotlib.colors import LightSource  # noqa: E402

from ffn_sim.cell.dcm_native_shell import (  # noqa: E402
    ResolvedNativeDCM, build_native_dcm_simulation, build_native_snapshot)
from ffn_sim.cell.dcm_confluence import (  # noqa: E402
    capture_positions, compute_confluence, per_cell_tris)
from ffn_sim.cell.dcm import _cluster_centers  # noqa: E402
from ffn_sim.common.sim_realtime import map_realtime  # noqa: E402

matplotlib.rcParams["animation.ffmpeg_path"] = "/opt/homebrew/bin/ffmpeg"

OUT = Path("ffn_sim/outputs/h_dcm_active/figs")
OUT.mkdir(parents=True, exist_ok=True)
MP4 = OUT / "spheroid_confluent_surface.mp4"
PNG = OUT / "spheroid_confluent_final.png"

# --- CONFLUENT regime from Part A (scripts/dcm_confluent_tune.py) -------------
# Chosen regime (Part A refinement, config 3): Φ=0.855, contact-node frac=0.839,
# Vcell/Vrest=1.22 (cells keep their volume AND deform to fill space) vs the
# separated-sphere baseline Φ=0.44, contact=0.00.
N_CELLS = 16
CONFLUENT = dict(
    subdivisions=2,         # 162 nodes/cell → smooth cell surfaces for rendering
    spacing_factor=1.8,     # cells start close enough for adhesion to bridge
    c_adh=5.0e-6,           # wide adhesion range → reaches & pulls neighbours
    adh_strength=8.0e8,     # strong cohesion → cells press together
    k_edge=5.0e-5,          # soft cortex → cells DEFORM into polygonal contact
    rep_strength=2.0e8,     # non-penetration → cells keep their volume (no implosion)
    turgor_inflate=1.08,    # moderate turgor → inflated, space-filling
    gamma_node=3.9e-10,
    dt=2.0e-10,
    force_cap=8.0e-8,
)

N_FRAMES = 60
STEPS_PER_FRAME = 120       # 60×120 = 7200 BAOAB steps total relaxation
FPS = 12

_TAB = plt.get_cmap("tab20").colors


def cell_colour(c: int):
    return _TAB[c % len(_TAB)]


def shade_faces(verts, tris, base_rgb, ls: LightSource):
    """Return per-face RGBA with simple Lambert shading from a LightSource."""
    tri = verts[tris]                                   # (k,3,3)
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    norm = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.where(norm > 0, norm, 1.0)
    # Lambert intensity from the light direction.
    az = np.deg2rad(ls.azdeg)
    alt = np.deg2rad(ls.altdeg)
    light = np.array([np.cos(alt) * np.cos(az),
                      np.cos(alt) * np.sin(az), np.sin(alt)])
    inten = np.clip(np.abs(n @ light), 0.0, 1.0)
    shade = 0.45 + 0.55 * inten                          # ambient + diffuse
    base = np.array(base_rgb)
    rgba = np.empty((tris.shape[0], 4))
    rgba[:, :3] = np.clip(base[None, :] * shade[:, None], 0, 1)
    rgba[:, 3] = 1.0
    return rgba


def draw_surfaces(ax, pos, ranges, tris_pc, ls, *, clip=None, alpha=1.0):
    """Draw each cell's triangulated shell as a shaded Poly3DCollection."""
    for c, (lo, hi) in enumerate(ranges):
        verts = pos[lo:hi]
        tris = tris_pc[c]
        if clip is not None:
            # cross-section: keep only triangles whose centroid is on one side.
            cen = verts[tris].mean(axis=1)
            keep = cen[:, clip[0]] <= clip[1]
            tris = tris[keep]
            if tris.shape[0] == 0:
                continue
        rgba = shade_faces(verts, tris, cell_colour(c), ls)
        polys = verts[tris]
        coll = Poly3DCollection(polys, facecolors=rgba, edgecolors=(0, 0, 0, 0.12),
                                linewidths=0.15)
        coll.set_alpha(alpha)
        ax.add_collection3d(coll)


def set_equal(ax, pos, pad=1.1):
    c = pos.mean(0)
    r = np.abs(pos - c).max() * pad
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1))


def main() -> None:
    p = ResolvedNativeDCM(**CONFLUENT)
    h = build_native_dcm_simulation(p, N_CELLS, substrate=False, contact=True)
    sim = h["sim"]
    ranges = h["ranges"]

    centers = _cluster_centers(N_CELLS, p.spacing_factor * p.R_cell,
                               p.z_substrate, p.R_cell, mode=p.cluster)
    _snap, mt, mty, *_rest = build_native_snapshot(p, N_CELLS, centers)
    tris_pc = per_cell_tris(mt, mty, ranges)

    ls = LightSource(azdeg=225, altdeg=45)

    # capture frame 0, then relax frame-by-frame
    frames = [capture_positions(sim)]
    metrics = []
    cr = p.c_adh
    metrics.append(compute_confluence(frames[0], cell_of_tag=h["cell_of_tag"],
                                      ranges=ranges, mesh_tris_per_cell=tris_pc,
                                      contact_range=cr))
    for f in range(1, N_FRAMES):
        sim.run(STEPS_PER_FRAME)
        pos = capture_positions(sim)
        frames.append(pos)
        if not np.isfinite(pos).all():
            print(f"[warn] non-finite at frame {f}; truncating", flush=True)
            frames = frames[:-1]
            break
    for pos in frames[1:]:
        metrics.append(compute_confluence(pos, cell_of_tag=h["cell_of_tag"],
                                          ranges=ranges, mesh_tris_per_cell=tris_pc,
                                          contact_range=cr))

    print(f"captured {len(frames)} frames; "
          f"Phi {metrics[0].packing_fraction:.3f}->{metrics[-1].packing_fraction:.3f} "
          f"contact {metrics[0].contact_fraction:.3f}->{metrics[-1].contact_fraction:.3f}",
          flush=True)

    # fixed view box from the LAST frame (most spread); um units for labels
    UM = 1e6
    last = frames[-1]

    fig = plt.figure(figsize=(7.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    writer = FFMpegWriter(fps=FPS, metadata=dict(artist="ffn_cellsim"),
                          bitrate=3200)

    n = len(frames)
    with writer.saving(fig, str(MP4), dpi=120):
        for i, pos in enumerate(frames):
            ax.clear()
            draw_surfaces(ax, pos * UM, ranges, tris_pc, ls)
            set_equal(ax, last * UM)
            ax.view_init(elev=22, azim=(360.0 * i / n))
            ax.set_axis_off()
            m = metrics[i]
            stage = ("separated FCC ball" if i < n * 0.18 else
                     "cells adhere + deform" if i < n * 0.6 else
                     "confluent spheroid")
            # per-frame real-time-equivalent (t_sim → t_real via accel S; a MAPPING)
            rt_f = map_realtime(i * STEPS_PER_FRAME, p.dt)
            ax.set_title(
                f"DCM spheroid — {N_CELLS} deformable cells  ·  {stage}\n"
                f"step {i*STEPS_PER_FRAME:>5d}   "
                f"t={rt_f.t_sim_human} sim ≈ {rt_f.t_real_human} real "
                f"(accel {rt_f.accel_factor:.0e})\n"
                f"packing Φ={m.packing_fraction:.2f}   "
                f"contact-node frac={m.contact_fraction:.2f}",
                fontsize=10)
            writer.grab_frame()
            if i % 10 == 0:
                print(f"  frame {i}/{n}", flush=True)
    plt.close(fig)
    print(f"wrote {MP4}", flush=True)

    # --- still PNG: final aggregate (full) + interior cross-section ----------
    fig2 = plt.figure(figsize=(13, 6.5))
    axL = fig2.add_subplot(121, projection="3d")
    draw_surfaces(axL, last * UM, ranges, tris_pc, ls)
    set_equal(axL, last * UM)
    axL.view_init(elev=22, azim=35)
    axL.set_axis_off()
    mfin = metrics[-1]
    axL.set_title(f"Confluent aggregate (surfaces)\n"
                  f"Φ={mfin.packing_fraction:.2f}  "
                  f"contact-node frac={mfin.contact_fraction:.2f}", fontsize=11)

    axR = fig2.add_subplot(122, projection="3d")
    mid = (last[:, 1].min() + last[:, 1].max()) * 0.5
    draw_surfaces(axR, last * UM, ranges, tris_pc, ls,
                  clip=(1, mid * UM))   # cut at y = mid → interior visible
    set_equal(axR, last * UM)
    axR.view_init(elev=22, azim=35)
    axR.set_axis_off()
    axR.set_title("Cross-section (interior packing)\n"
                  "deformed cells fill space, no gaps", fontsize=11)
    fig2.suptitle(
        f"Native-mesh DCM — {N_CELLS} cells relaxed into a confluent tissue "
        f"(SimuCell3D-style)", fontsize=13)
    fig2.tight_layout()
    fig2.savefig(PNG, dpi=140)
    plt.close(fig2)
    print(f"wrote {PNG}", flush=True)


if __name__ == "__main__":
    main()
