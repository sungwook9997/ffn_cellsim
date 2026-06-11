"""SEPARATE lamellipodium visualisation for the DCM spheroid (PI 2026-06-11).

PI asked to see the LAMELLIPODIUM on its own so it can be checked. In the DCM
spheroid the lamellipodium is the coarse-grained engine of
``cell/dcm_active.py::ActiveRimTraction`` (= ``DcmActiveRimTractionGPU(Vec)``): a
RIM cell (on the colony periphery — few neighbours — AND with basal nodes touching
the dish) puts down a basal PROTRUSION and pulls OUTWARD (radially from the colony
centroid). The leading-edge basal nodes pull extra; a weak apical belt drags the
body after the edge. This script renders THAT, faithfully, on the now-COHESIVE
spheroid (cells adhered from t=0 — the cohesion fix, scripts/dcm_cohesion_check.py):

  Panel A  top-down (x-y): the colony footprint. Basal lamellipodium nodes
           (z < z_c) drawn as the wetted protrusion layer; RIM (lamellipodiating)
           cells highlighted; OUTWARD traction arrows (the lamellipodial pull) at
           each rim cell — the advancing lamellipodial front around the perimeter.
  Panel B  side (x-z): rounded cell BODIES above z_c vs the thin basal
           lamellipodium layer wetting the dish below z_c (dashed line = z_c).
  Panel C  single rim-cell ANATOMY: one lamellipodiating cell — basal lamellipodium
           nodes vs apical body, LEADING-EDGE basal nodes (extra pull) highlighted,
           and the outward traction vector. The lamellipodium at the single-cell
           scale.

The rim / basal / outward-direction geometry is recomputed here EXACTLY as the
force law computes it (neighbour-crowding rim test, z<z_c basal split, in-plane
centroid→cell outward unit vector) and cross-checked against the live force
object's ``rim_cells`` so the picture matches the physics, not a proxy.

Run (CPU dev on this Mac is fine):
    conda activate ffn_sim
    PYTHONPATH=. python ffn_sim/scripts/dcm_lamellipodium_viz.py --n 60 --steps 9000

Outputs (ffn_sim/outputs/h_dcm_gpu_lod/figs/):
    lamellipodium_anatomy.png   — the 3-panel still (the check figure)
    lamellipodium_front.mp4     — the lamellipodial front developing over spread time
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

from ffn_sim.cell.dcm import icosphere_mesh  # noqa: E402
from ffn_sim.cell.dcm_gpu_build import (  # noqa: E402
    ResolvedGpuDCM, build_gpu_dcm_simulation)
from ffn_sim.cell.dcm_confluence import capture_positions  # noqa: E402
from ffn_sim.common.sim_realtime import map_realtime  # noqa: E402

for cand in ("/opt/homebrew/bin/ffmpeg", "/usr/bin/ffmpeg", "/usr/local/bin/ffmpeg"):
    if Path(cand).exists():
        matplotlib.rcParams["animation.ffmpeg_path"] = cand
        break

OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod/figs")
OUT.mkdir(parents=True, exist_ok=True)
UM = 1.0e6
_TAB = plt.get_cmap("tab20").colors

# lamellipodium colour language (consistent across panels)
C_BODY = "#9ecae1"        # rounded cell body (apical)
C_LAMEL = "#fb6a4a"       # basal lamellipodium protrusion
C_LEAD = "#a50f15"        # leading-edge basal nodes (extra pull)
C_ARROW = "#08519c"       # outward traction vector


def cell_colour(c: int):
    return _TAB[c % len(_TAB)]


# ---------------------------------------------------------------------------
# Lamellipodium geometry — recomputed EXACTLY as ActiveRimTraction.set_forces.
# ---------------------------------------------------------------------------
def lamellipodium_geometry(pos: np.ndarray, ranges, traction) -> dict:
    """Return rim cells, per-cell basal/leading masks, outward dirs, traction mag.

    Mirrors ``cell/dcm_active.py::ActiveRimTraction.set_forces``: rim = active cell
    with neighbour-crowd ≤ max_neigh AND ≥1 basal node (z < z_c); outward dir =
    in-plane unit (cell centroid − colony centroid); leading edge = basal nodes on
    the outward side. ``traction`` supplies R, z0, r_neigh, max_neigh, contact_band,
    f_per_cell (live).
    """
    R, z0 = traction.R, traction.z0
    r_neigh, max_neigh = traction.r_neigh, traction.max_neigh
    zc = z0 + traction.contact_band * R
    n_cells = len(ranges)
    cents = np.array([pos[lo:hi].mean(0) for lo, hi in ranges])
    cluster_cen = cents.mean(0)
    d2 = np.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
    within = d2 < r_neigh ** 2
    np.fill_diagonal(within, False)
    crowd = within.sum(axis=1)

    rim, basal_idx, lead_idx, rhat_of, fmag_of = [], {}, {}, {}, {}
    for c in range(n_cells):
        if crowd[c] > max_neigh:
            continue
        lo, hi = ranges[c]
        cp = pos[lo:hi]
        basal = np.where(cp[:, 2] < zc)[0]
        if basal.size == 0:
            continue
        rxy = cents[c][:2] - cluster_cen[:2]
        rn = float(np.hypot(*rxy))
        if rn < 1e-12:
            continue
        rhat = np.array([rxy[0] / rn, rxy[1] / rn, 0.0])
        off = cp[basal, :2] - cents[c][:2]
        lead = basal[(off @ rhat[:2]) > 0.0]
        rim.append(c)
        basal_idx[c] = basal
        lead_idx[c] = lead
        rhat_of[c] = rhat
        fmag_of[c] = float(traction.f_per_cell.get(c, traction.f_act))
    return dict(rim=rim, basal_idx=basal_idx, lead_idx=lead_idx, rhat=rhat_of,
                fmag=fmag_of, cents=cents, cluster_cen=cluster_cen, zc=zc,
                R=R, z0=z0)


def basal_footprint_area(pos: np.ndarray, z0: float, R: float) -> float:
    """xy convex-hull area [µm²] of basal contact nodes (the wetted footprint)."""
    from scipy.spatial import ConvexHull

    basal = pos[pos[:, 2] < z0 + 0.5 * R][:, :2]
    if basal.shape[0] < 3:
        return 0.0
    try:
        return float(ConvexHull(basal).volume) * UM * UM
    except Exception:  # noqa: BLE001 — collinear/degenerate
        return float(np.ptp(basal[:, 0]) * np.ptp(basal[:, 1])) * UM * UM


# ---------------------------------------------------------------------------
# Build + run, capturing frames AND the live force-object rim state per frame.
# ---------------------------------------------------------------------------
def build_and_capture(*, n_cells, n_frames, steps_per_frame, gate_steps):
    p = ResolvedGpuDCM(R_cell=7.5e-6, subdivisions=1)
    h = build_gpu_dcm_simulation(p, n_cells, active=True)   # cohesive + wetting + active
    sim, ranges, traction = h["sim"], h["ranges"], h["traction"]
    _v, _e, tris0 = icosphere_mesh(p.R_cell, p.subdivisions)

    sim.run(0)
    t0 = time.time()
    sim.run(gate_steps)
    traction.set_forces(sim.timestep)
    print(f"[gate] run({gate_steps}) in {time.time()-t0:.0f}s; "
          f"rim cells={len(traction.rim_cells)} / {h['n_cells']}", flush=True)

    frames, geos, steps = [], [], []

    def grab():
        pos = capture_positions(sim)
        traction.set_forces(sim.timestep)          # refresh rim_cells / f_per_cell
        geo = lamellipodium_geometry(pos, ranges, traction)
        # cross-check our rim against the live force object's rim_cells
        live = set(int(x) for x in traction.rim_cells)
        ours = set(geo["rim"])
        if live and live != ours:
            print(f"  [note] rim set differs from force obj: "
                  f"ours\\live={sorted(ours-live)[:4]} live\\ours={sorted(live-ours)[:4]}",
                  flush=True)
        return pos, geo

    pos0, geo0 = grab()
    frames.append(pos0); geos.append(geo0); steps.append(sim.timestep)
    for f in range(1, n_frames):
        sim.run(steps_per_frame)
        pos = capture_positions(sim)
        if not np.isfinite(pos).all():
            print(f"[warn] non-finite at frame {f}; truncating", flush=True)
            break
        _, geo = grab()
        frames.append(pos); geos.append(geo); steps.append(sim.timestep)
        area = basal_footprint_area(pos, p.z_substrate, p.R_cell)
        print(f"  frame {f:>2d} step {sim.timestep:>6d}  rim={len(geo['rim']):>3d}  "
              f"footprint={area:7.0f}µm²", flush=True)
    # cache the (expensive) sim so render tweaks don't re-simulate
    import pickle
    cache = OUT / f"_lamel_cache_n{len(ranges)}.pkl"
    with open(cache, "wb") as fh:
        pickle.dump(dict(frames=frames, geos=geos, steps=steps,
                         ranges=ranges, tris0=tris0, p=p), fh)
    print(f"[cache] wrote {cache}", flush=True)
    return frames, geos, steps, ranges, tris0, p


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def _arrow_scale(geos, ref_len):
    """A common arrow length so the longest traction ≈ ref_len (µm) in plot units."""
    fmax = max((max(g["fmag"].values()) for g in geos if g["fmag"]), default=1.0)
    return ref_len / max(fmax, 1e-30)


def panelA_full(ax, pos, geo, ranges, *, arrow_k):
    """Top-down done right: per-cell basal + leading nodes from ranges blocks."""
    R = geo["R"]
    # faint: ALL basal nodes (the wetted footprint layer)
    basal_all = pos[pos[:, 2] < geo["zc"]]
    ax.scatter(basal_all[:, 0] * UM, basal_all[:, 1] * UM, s=5, c="0.80",
               edgecolors="none", zorder=1)
    # outline the WETTED FOOTPRINT (xy convex hull of basal nodes) so the spread
    # area the lamellipodial front is advancing is unmistakable.
    try:
        from scipy.spatial import ConvexHull
        xy = basal_all[:, :2] * UM
        if xy.shape[0] >= 3:
            hull = ConvexHull(xy)
            loop = np.append(hull.vertices, hull.vertices[0])
            ax.plot(xy[loop, 0], xy[loop, 1], "-", color="0.45", lw=1.2, zorder=2)
            ax.fill(xy[loop, 0], xy[loop, 1], color="0.85", alpha=0.25, zorder=0)
            ax.plot([], [], color="0.45", lw=1.2,
                    label=f"wetted footprint {hull.volume:.0f}µm²")
    except Exception:  # noqa: BLE001 — degenerate hull
        pass
    for c in geo["rim"]:
        lo, hi = ranges[c]
        cp = pos[lo:hi]
        b = cp[geo["basal_idx"][c]]
        ax.scatter(b[:, 0] * UM, b[:, 1] * UM, s=14, c=C_LAMEL,
                   edgecolors="none", zorder=3)
        if geo["lead_idx"][c].size:
            le = cp[geo["lead_idx"][c]]
            ax.scatter(le[:, 0] * UM, le[:, 1] * UM, s=20, c=C_LEAD,
                       edgecolors="none", zorder=4)
        cen = geo["cents"][c]; rhat = geo["rhat"][c]; fmag = geo["fmag"][c]
        L = fmag * arrow_k
        ax.arrow(cen[0] * UM, cen[1] * UM, rhat[0] * L, rhat[1] * L,
                 head_width=0.30 * R * UM, head_length=0.38 * R * UM,
                 fc=C_ARROW, ec=C_ARROW, length_includes_head=True, zorder=6,
                 width=0.05 * R * UM)
    ax.scatter([geo["cluster_cen"][0] * UM], [geo["cluster_cen"][1] * UM],
               marker="+", s=90, c="k", zorder=7)
    # legend proxies
    ax.scatter([], [], s=14, c=C_LAMEL, label="rim basal lamellipodium")
    ax.scatter([], [], s=20, c=C_LEAD, label="leading edge (extra pull)")
    ax.plot([], [], c=C_ARROW, lw=2, label="outward traction")
    ax.scatter([], [], s=5, c="0.80", label="basal (wetted) nodes")
    ax.set_aspect("equal")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]")
    ax.legend(loc="upper right", fontsize=6, framealpha=0.9)


def panelB_side(ax, pos, geo, ranges):
    """Side (x-z): cell bodies above z_c vs thin basal lamellipodium below z_c."""
    R, z0, zc = geo["R"], geo["z0"], geo["zc"]
    apical = pos[pos[:, 2] >= zc]
    basal = pos[pos[:, 2] < zc]
    ax.scatter(apical[:, 0] * UM, apical[:, 2] * UM, s=5, c=C_BODY,
               edgecolors="none", zorder=2)
    ax.scatter(basal[:, 0] * UM, basal[:, 2] * UM, s=8, c=C_LAMEL,
               edgecolors="none", zorder=3)
    xlo, xhi = pos[:, 0].min() * UM, pos[:, 0].max() * UM
    ax.axhline(z0 * UM, color="0.4", lw=1.4)                       # substrate
    ax.axhline(zc * UM, color="0.4", lw=0.9, ls="--")             # basal/apical split
    ax.text(xhi, zc * UM, " z_c", va="bottom", ha="right", fontsize=7, color="0.4")
    ax.text(xlo, z0 * UM, " substrate", va="bottom", ha="left", fontsize=7, color="0.4")
    ax.scatter([], [], s=5, c=C_BODY, label="cell body (apical)")
    ax.scatter([], [], s=8, c=C_LAMEL, label="basal lamellipodium")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("z [µm]")
    ax.legend(loc="upper right", fontsize=6, framealpha=0.9)


def panelC_single(ax, pos, geo, ranges, tris0):
    """Single rim-cell anatomy (3D): basal lamellipodium / leading edge / body / arrow."""
    if not geo["rim"]:
        ax.text2D(0.5, 0.5, "no rim cells yet", ha="center", transform=ax.transAxes)
        return
    # pick the rim cell with the most leading-edge nodes (clearest lamellipodium)
    c = max(geo["rim"], key=lambda c: geo["lead_idx"][c].size)
    lo, hi = ranges[c]
    cp = pos[lo:hi]
    zc = geo["zc"]
    basal = cp[:, 2] < zc
    lead = np.zeros(len(cp), bool); lead[geo["lead_idx"][c]] = True
    # body surface (faint)
    polys = cp[tris0]
    ax.add_collection3d(Poly3DCollection(polys * UM, facecolors=(0.62, 0.79, 0.88, 0.18),
                                         edgecolors=(0, 0, 0, 0.06), linewidths=0.1))
    body = cp[~basal]
    ax.scatter(body[:, 0] * UM, body[:, 1] * UM, body[:, 2] * UM, s=10, c=C_BODY,
               depthshade=False, label="cell body (apical)")
    bas = cp[basal & ~lead]
    ax.scatter(bas[:, 0] * UM, bas[:, 1] * UM, bas[:, 2] * UM, s=22, c=C_LAMEL,
               depthshade=False, label="basal lamellipodium")
    led = cp[lead]
    if led.size:
        ax.scatter(led[:, 0] * UM, led[:, 1] * UM, led[:, 2] * UM, s=34, c=C_LEAD,
                   depthshade=False, label="leading edge")
    cen = geo["cents"][c]; rhat = geo["rhat"][c]
    L = 1.6 * geo["R"]
    ax.quiver(cen[0] * UM, cen[1] * UM, (geo["z0"] + 0.3 * geo["R"]) * UM,
              rhat[0] * L * UM, rhat[1] * L * UM, 0.0,
              color=C_ARROW, lw=2.4, arrow_length_ratio=0.3)
    ax.set_xlabel("x [µm]", fontsize=7); ax.set_ylabel("y [µm]", fontsize=7)
    ax.set_zlabel("z [µm]", fontsize=7)
    ax.legend(loc="upper left", fontsize=6)
    ax.set_title(f"single rim cell #{c} — lamellipodium anatomy", fontsize=8)


# ---------------------------------------------------------------------------
# SURFACE rendering — real triangulated cell SURFACES with the lamellipodium
# made VISIBLE (PI: "시각화도 같이 그리고 lamellipodium 보이게" — surfaces, not dots).
# ---------------------------------------------------------------------------
def _lambert(verts, tris, base_rgb, ls, alpha):
    """Per-face RGBA, simple Lambert shading (à la dcm_spreading_surface_mp4)."""
    tri = verts[tris]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    nn = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.where(nn > 0, nn, 1.0)
    az, alt = np.deg2rad(ls.azdeg), np.deg2rad(ls.altdeg)
    light = np.array([np.cos(alt) * np.cos(az), np.cos(alt) * np.sin(az), np.sin(alt)])
    shade = 0.5 + 0.5 * np.clip(np.abs(n @ light), 0, 1)
    rgba = np.empty((tris.shape[0], 4))
    rgba[:, :3] = np.clip(np.array(base_rgb)[None, :] * shade[:, None], 0, 1)
    rgba[:, 3] = alpha
    return rgba


def draw_lamellipodium_surfaces(ax, pos, geo, ranges, tris0, ls):
    """Draw every cell as a shaded surface; the LAMELLIPODIUM is made visible as the
    bright basal cap of the rim (lamellipodiating) cells.

    Colour language: BODY cells = cool blue (translucent); RIM cell apical body =
    warm orange; RIM cell BASAL faces (z<z_c) = bright red = THE LAMELLIPODIUM (the
    protrusion wetting the dish). So the red caps spreading on the substrate ARE the
    lamellipodia, unmistakable against the blue body.
    """
    from matplotlib.colors import to_rgb
    rim = set(geo["rim"]); zc = geo["zc"]
    body_rgb, rim_body_rgb, lamel_rgb = to_rgb(C_BODY), to_rgb("#f0a05a"), to_rgb(C_LEAD)
    for c, (lo, hi) in enumerate(ranges):
        verts = pos[lo:hi]
        polys = verts[tris0] * UM
        if c in rim:
            # split faces into basal (lamellipodium) vs apical (rim body) by face-centroid z
            fz = verts[tris0][:, :, 2].mean(axis=1)
            basal_f = fz < zc
            rgba = _lambert(verts, tris0, rim_body_rgb, ls, 0.95)
            rgba[basal_f] = _lambert(verts, tris0, lamel_rgb, ls, 1.0)[basal_f]
            ax.add_collection3d(Poly3DCollection(polys, facecolors=rgba,
                                edgecolors=(0, 0, 0, 0.12), linewidths=0.12))
        else:
            rgba = _lambert(verts, tris0, body_rgb, ls, 0.38)   # translucent so the
            ax.add_collection3d(Poly3DCollection(polys, facecolors=rgba,  # red shows
                                edgecolors=(0, 0, 0, 0.04), linewidths=0.06))
    # outward traction arrows on the lamellipodiating cells (the protrusion pull)
    for c in geo["rim"]:
        cen = geo["cents"][c]; rhat = geo["rhat"][c]
        L = 2.0 * geo["R"]
        ax.quiver(cen[0] * UM, cen[1] * UM, (geo["z0"] + 0.25 * geo["R"]) * UM,
                  rhat[0] * L * UM, rhat[1] * L * UM, 0.0, color=C_ARROW,
                  lw=1.8, arrow_length_ratio=0.35, zorder=10)


def draw_substrate(ax, xlim, ylim, z0_um):
    xs = np.array([[xlim[0], xlim[1]], [xlim[0], xlim[1]]])
    ys = np.array([[ylim[0], ylim[0]], [ylim[1], ylim[1]]])
    ax.plot_surface(xs, ys, np.full_like(xs, z0_um), color="0.6", alpha=0.18, shade=False)


def render_surface_figure(frames, geos, steps, ranges, tris0, p, xlim, ylim):
    """Surface still (3D + side) + MP4 — the lamellipodium VISIBLE on real surfaces."""
    ls = LightSource(azdeg=225, altdeg=45)
    z0_um = p.z_substrate * UM
    zlo = z0_um - 2.0
    zhi = max(np.concatenate(frames)[:, 2]) * UM * 1.08

    def setup3d(ax, elev, azim):
        draw_substrate(ax, xlim, ylim, z0_um)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_zlim(zlo, zhi)
        ax.set_box_aspect((xlim[1]-xlim[0], ylim[1]-ylim[0], (zhi-zlo)))
        ax.view_init(elev=elev, azim=azim)
        ax.set_xlabel("x [µm]", fontsize=7); ax.set_ylabel("y [µm]", fontsize=7)
        ax.set_zlabel("z [µm]", fontsize=7)

    pos, geo = frames[-1], geos[-1]
    rt = map_realtime(int(steps[-1]), float(p.dt))
    fig = plt.figure(figsize=(15, 6.2))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    draw_lamellipodium_surfaces(ax1, pos, geo, ranges, tris0, ls)
    setup3d(ax1, elev=22, azim=-60)
    ax1.set_title("3D perspective — red basal caps = lamellipodia wetting the dish", fontsize=9)
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    draw_lamellipodium_surfaces(ax2, pos, geo, ranges, tris0, ls)
    setup3d(ax2, elev=4, azim=-90)            # near-side view: protrusions at the base
    ax2.set_title("side view — thin lamellipodium protrusions vs rounded bodies", fontsize=9)
    # colour legend
    from matplotlib.patches import Patch
    fig.legend(handles=[
        Patch(facecolor=C_LEAD, label="lamellipodium (rim cell basal cap)"),
        Patch(facecolor="#f0a05a", label="rim cell body"),
        Patch(facecolor=C_BODY, label="interior / body cells"),
    ], loc="lower center", ncol=3, fontsize=8, framealpha=0.9)
    fig.suptitle(
        f"DCM spheroid spreading — LAMELLIPODIA VISIBLE (surfaces) · N={len(ranges)}, "
        f"{len(geo['rim'])} lamellipodiating · step {steps[-1]} · "
        f"{rt.t_sim_human} sim ≈ {rt.t_real_human} real", fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    png = OUT / "lamellipodium_surface.png"
    fig.savefig(png, dpi=145); plt.close(fig)
    print(f"wrote {png}", flush=True)

    # MP4: spreading with lamellipodia visible (3D perspective)
    mp4 = OUT / "lamellipodium_surface.mp4"
    figm = plt.figure(figsize=(8.5, 7))
    axm = figm.add_subplot(1, 1, 1, projection="3d")
    writer = FFMpegWriter(fps=6, bitrate=3000)
    try:
        with writer.saving(figm, str(mp4), dpi=120):
            for pos, geo, st in zip(frames, geos, steps):
                axm.clear()
                draw_lamellipodium_surfaces(axm, pos, geo, ranges, tris0, ls)
                setup3d(axm, elev=20, azim=-60)
                rt = map_realtime(int(st), float(p.dt))
                axm.set_title(f"step {st} · {len(geo['rim'])} lamellipodia · "
                              f"{rt.t_sim_human} sim ≈ {rt.t_real_human} real", fontsize=10)
                writer.grab_frame()
        print(f"wrote {mp4}", flush=True)
    except Exception as e:  # noqa: BLE001
        print(f"[warn] surface MP4 skipped: {e}", flush=True)
    plt.close(figm)


def fixed_box(frames, geo0):
    allp = np.concatenate(frames, axis=0) * UM
    cx = 0.5 * (allp[:, 0].min() + allp[:, 0].max())
    cy = 0.5 * (allp[:, 1].min() + allp[:, 1].max())
    r = 0.5 * max(np.ptp(allp[:, 0]), np.ptp(allp[:, 1])) * 1.12
    return (cx - r, cx + r), (cy - r, cy + r)


# ---------------------------------------------------------------------------
# Compose: still PNG + front MP4
# ---------------------------------------------------------------------------
def render(frames, geos, steps, ranges, tris0, p):
    xlim, ylim = fixed_box(frames, geos[0])
    arrow_k = _arrow_scale(geos, ref_len=2.2 * p.R_cell)   # plot-unit length per N
    z0 = p.z_substrate

    # ---- still: last frame, 3 panels ----
    pos, geo = frames[-1], geos[-1]
    rt = map_realtime(int(steps[-1]), float(p.dt))
    fig = plt.figure(figsize=(16, 5.2))
    axA = fig.add_subplot(1, 3, 1)
    panelA_full(axA, pos, geo, ranges, arrow_k=arrow_k)
    axA.set_xlim(*xlim); axA.set_ylim(*ylim)            # xlim/ylim already in µm
    axA.set_title("A · top-down: lamellipodial front (outward traction)", fontsize=9)
    axB = fig.add_subplot(1, 3, 2)
    panelB_side(axB, pos, geo, ranges)
    axB.set_xlim(*xlim)
    axB.set_title("B · side: basal lamellipodium vs cell body", fontsize=9)
    axC = fig.add_subplot(1, 3, 3, projection="3d")
    panelC_single(axC, pos, geo, ranges, tris0)
    rt_lbl = f"{rt.t_sim_human} sim ≈ {rt.t_real_human} real (accel {rt.accel_factor:.0e})"
    fig.suptitle(
        f"DCM spheroid LAMELLIPODIUM — N={len(ranges)} cells, R_cell={p.R_cell*UM:.1f}µm, "
        f"step {steps[-1]}, {len(geo['rim'])} rim cells lamellipodiating · {rt_lbl}",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    png = OUT / "lamellipodium_anatomy.png"
    fig.savefig(png, dpi=140); plt.close(fig)
    print(f"wrote {png}", flush=True)

    # ---- MP4: the lamellipodial front developing (top-down + side) ----
    mp4 = OUT / "lamellipodium_front.mp4"
    figm = plt.figure(figsize=(12, 5.4))
    axm1 = figm.add_subplot(1, 2, 1)
    axm2 = figm.add_subplot(1, 2, 2)
    writer = FFMpegWriter(fps=8, bitrate=2400)
    try:
        with writer.saving(figm, str(mp4), dpi=120):
            for pos, geo, st in zip(frames, geos, steps):
                axm1.clear(); axm2.clear()
                panelA_full(axm1, pos, geo, ranges, arrow_k=arrow_k)
                axm1.set_xlim(*xlim); axm1.set_ylim(*ylim)     # already in µm
                axm1.set_title("lamellipodial front (top-down)", fontsize=9)
                panelB_side(axm2, pos, geo, ranges)
                axm2.set_xlim(*xlim)
                axm2.set_title("basal lamellipodium (side)", fontsize=9)
                rt = map_realtime(int(st), float(p.dt))
                figm.suptitle(f"step {st} · {len(geo['rim'])} rim cells · "
                              f"{rt.t_sim_human} sim ≈ {rt.t_real_human} real", fontsize=10)
                figm.tight_layout(rect=(0, 0, 1, 0.94))
                writer.grab_frame()
        print(f"wrote {mp4}", flush=True)
    except Exception as e:  # noqa: BLE001 — ffmpeg missing etc.
        print(f"[warn] MP4 skipped: {e}", flush=True)
    plt.close(figm)

    # ---- SURFACE figure + MP4: the lamellipodium VISIBLE on real cell surfaces ----
    render_surface_figure(frames, geos, steps, ranges, tris0, p, xlim, ylim)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--frames", type=int, default=12)
    ap.add_argument("--steps", type=int, default=9000, help="total spread steps")
    ap.add_argument("--gate", type=int, default=4500)
    ap.add_argument("--render-only", action="store_true",
                    help="re-render from the cached sim (no re-simulation)")
    args = ap.parse_args()
    if args.render_only:
        import pickle
        cache = OUT / f"_lamel_cache_n{args.n}.pkl"
        with open(cache, "rb") as fh:
            d = pickle.load(fh)
        print(f"[cache] loaded {cache}", flush=True)
        render(d["frames"], d["geos"], d["steps"], d["ranges"], d["tris0"], d["p"])
        return
    spf = max(1, args.steps // max(1, args.frames - 1))
    frames, geos, steps, ranges, tris0, p = build_and_capture(
        n_cells=args.n, n_frames=args.frames, steps_per_frame=spf, gate_steps=args.gate)
    render(frames, geos, steps, ranges, tris0, p)


if __name__ == "__main__":
    main()
