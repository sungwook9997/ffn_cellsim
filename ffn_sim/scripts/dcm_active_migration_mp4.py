"""Migration MP4 of the ACTIVE-spreading DCM spheroid — make the PULLING VISIBLE.

The morphology MP4 (`dcm_active_morphology_mp4.py`) shows the basal footprint
spreading, but the PI's critique was that you cannot SEE the rim cells MIGRATING
outward (the lamellipodium pulling the cell body + followers dragged by cohesion).
This script renders that explicitly. It re-runs the active spheroid (same
build/run loop) capturing PER-FRAME:

  * per-CELL centroids (so we can draw migration TRAJECTORY TRAILS — fading tails
    over the last K frames — that show each cell crawling outward),
  * the ACTIVE TRACTION VECTORS on the rim cells (direction + magnitude of the
    active pull, read live from ``ActiveRimTraction.cell_net_force``),
  * per-node positions + 3-zone state + junction-switched mask.

It renders an MP4 (ffmpeg writer, NOT gif) with two synchronized views:

  * top-down (x-y): centroid trails + traction arrows on the rim, footprint nodes,
  * side (x-z):     the ball flattening + crawling onto the substrate.

Cells coloured by 3-zone state; junction-switched cells ringed; title tracks
step / N cells / A/A0. The collective migration (leaders out front, followers
dragged behind by tent cohesion) is the thing made visible.

  ~/miniconda3/envs/ffn_sim/bin/python ffn_sim/scripts/dcm_active_migration_mp4.py [--quick]

Output: ffn_sim/outputs/h_dcm_active/figs/active_spheroid_migration.mp4
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

from ffn_sim.archive.hoomd_legacy.cell.dcm_active import (
    ResolvedActiveSpheroid, build_active_spheroid, CellState,
)
from ffn_sim.scripts.dcm_active_spheroid import (
    _positions, _active_node_mask, _footprint_area,
)
from ffn_sim.common.sim_realtime import map_realtime, format_time

_UM = 1e6
_STATE_COLOR = {
    int(CellState.PROLIFERATING): "#27ae60",   # green rim
    int(CellState.QUIESCENT): "#e67e22",        # amber middle
    int(CellState.NECROTIC): "#4d4d4d",         # grey core
}
_STATE_LABEL = {
    int(CellState.PROLIFERATING): "proliferating",
    int(CellState.QUIESCENT): "quiescent",
    int(CellState.NECROTIC): "necrotic",
}


def _centroids(pg, ranges, ids):
    return np.array([pg[ranges[int(c)][0]:ranges[int(c)][1]].mean(0) for c in ids])


def capture(p: ResolvedActiveSpheroid, *, n_active, n_max, dt=3.0e-10,
            aggreg_blocks=6, spread_blocks=24, block=1000, div_every=20000):
    """Build + run the active spheroid; return per-frame snapshots for the MP4.

    Each frame records node positions+state, per-cell centroids, and the live
    active-traction net force per rim cell (for arrow rendering). Division is SLOW
    (div_every large) so the migration, not division stacking, drives the spread.
    """
    h = build_active_spheroid(p, n_active=n_active, n_max=n_max,
                              integrin_substrate=True, belt=True)
    sim = h["sim"]; st = h["st"]; active = h["active"]
    cell_of_node = h["cell_of_node"]; ranges = h["ranges"]; tr = h["traction"]
    z0 = p.z_substrate; R = p.R_cell
    necrosis = h["necrosis"]; pressure = h["pressure"]
    junction = h["junction"]; prolif = h["prolif"]
    necrosis.attach(sim); pressure.attach(sim)
    junction.attach(sim); prolif.attach(sim)

    sim.run(0)
    pg0 = _positions(sim)
    A0 = _footprint_area(pg0, _active_node_mask(cell_of_node), z0, R)
    sim.run(2000)
    if not np.all(np.isfinite(_positions(sim))):
        if dt > 1.1e-10:
            return capture(p, n_active=n_active, n_max=n_max, dt=1.0e-10,
                           aggreg_blocks=aggreg_blocks, spread_blocks=spread_blocks,
                           block=block, div_every=div_every)
        raise RuntimeError("non-finite at finite gate")

    frames = []

    def snap(phase):
        pg = _positions(sim)
        con = cell_of_node
        mask = con >= 0
        A = _footprint_area(pg, mask, z0, R)
        ids = np.where(active)[0]
        cents = _centroids(pg, ranges, ids) if ids.size else np.empty((0, 3))
        # per-cell traction net force (rim only) keyed by cell id
        tractions = {int(c): tr.cell_net_force.get(int(c)) for c in ids
                     if tr.cell_net_force.get(int(c)) is not None}
        frames.append(dict(
            t=int(sim.timestep), dt=float(dt), phase=phase,
            xyz_um=(pg[mask]) * _UM,
            node_state=st.state[con[mask]].copy(),
            node_switched=st.switched[con[mask]].copy(),
            cell_ids=ids.copy(),
            cent_um=cents * _UM,
            cent_state=st.state[ids].copy(),
            cent_switched=st.switched[ids].copy(),
            tractions={c: f.copy() for c, f in tractions.items()},
            rim_cells=tr.rim_cells.copy(),
            n_cells=int(active.sum()), AoverA0=float(A / A0), z0_um=z0 * _UM,
        ))

    def tick():
        necrosis.act(int(sim.timestep)); pressure.act(int(sim.timestep))
        junction.act(int(sim.timestep))

    tick(); snap("aggregation")
    for _ in range(aggreg_blocks):
        sim.run(block); tick(); snap("aggregation")
    for b in range(spread_blocks):
        sim.run(block); tick()
        if int(sim.timestep) % div_every < block:
            prolif.act(int(sim.timestep))
        snap("spreading")
    return frames, h


def _radial_centroid_migration(frames):
    """Mean + max OUTWARD radial centroid displacement of rim cells, agg-end -> end.

    Tracks cells present from the aggregation-end frame to the last frame (the
    stable initial cohort) and reports their radial displacement from the cluster
    centroid. Returns (mean_rim_dr_um, max_dr_um, agg_idx)."""
    agg_idx = max(0, sum(1 for f in frames if f["phase"] == "aggregation") - 1)
    f0, f1 = frames[agg_idx], frames[-1]
    ids0 = list(f0["cell_ids"]); ids1 = list(f1["cell_ids"])
    common = [c for c in ids0 if c in ids1]
    if len(common) < 2:
        return 0.0, 0.0, agg_idx
    c0 = {int(c): f0["cent_um"][ids0.index(c)] for c in common}
    c1 = {int(c): f1["cent_um"][ids1.index(c)] for c in common}
    cc0 = np.mean([c0[c] for c in common], axis=0)
    cc1 = np.mean([c1[c] for c in common], axis=0)
    dr = {}
    for c in common:
        r0 = np.hypot(c0[c][0] - cc0[0], c0[c][1] - cc0[1])
        r1 = np.hypot(c1[c][0] - cc1[0], c1[c][1] - cc1[1])
        dr[c] = r1 - r0
    rim = set(int(x) for x in f1["rim_cells"])
    rim_dr = [dr[c] for c in common if c in rim] or list(dr.values())
    return float(np.mean(rim_dr)), float(np.max(list(dr.values()))), agg_idx


def render(frames, out_path, trail_k=8):
    xy_all = np.vstack([f["xyz_um"] for f in frames])
    xlim = (xy_all[:, 0].min() - 5, xy_all[:, 0].max() + 5)
    ylim = (xy_all[:, 1].min() - 5, xy_all[:, 1].max() + 5)
    zlim = (min(xy_all[:, 2].min(), frames[0]["z0_um"]) - 3, xy_all[:, 2].max() + 5)

    # typical traction magnitude (for arrow scaling) — use the median nonzero |F|
    fmags = [np.linalg.norm(F) for f in frames for F in f["tractions"].values()]
    fmags = [m for m in fmags if m > 0]
    f_ref = np.median(fmags) if fmags else 1e-8
    arrow_um = 6.0  # an f_ref-magnitude pull renders as this many µm long

    fig, (axxy, axxz) = plt.subplots(1, 2, figsize=(13.5, 6.4))
    # REAL-TIME readout (PI 2026-06-11): spreading-front accel factor maps the
    # accelerated-sim time to a real-time equivalent (a MAPPING, not native).
    _dt = float(frames[-1].get("dt", 3.0e-10))
    _rt_end = map_realtime(frames[-1]["t"], _dt)
    fig.suptitle("ACTIVE DCM spheroid — TRACTION-DRIVEN migration "
                 "(centroid trails + active-pull arrows)\n"
                 f" end: t = {_rt_end.t_sim_human} sim  ≈ {_rt_end.t_real_human} real  "
                 f"(accel {_rt_end.accel_factor:.0e}, spreading-front)",
                 fontweight="bold")

    mean_dr, max_dr, agg_idx = _radial_centroid_migration(frames)

    def draw(i):
        f = frames[i]
        for ax in (axxy, axxz):
            ax.clear()
        xyz = f["xyz_um"]; nst = f["node_state"]; nsw = f["node_switched"]
        # faint footprint nodes coloured by state
        for sval, col in _STATE_COLOR.items():
            m = nst == sval
            if m.any():
                axxy.scatter(xyz[m, 0], xyz[m, 1], s=6, c=col, alpha=0.25,
                             edgecolors="none")
                axxz.scatter(xyz[m, 0], xyz[m, 2], s=6, c=col, alpha=0.25,
                             edgecolors="none")
        if nsw.any():
            axxy.scatter(xyz[nsw, 0], xyz[nsw, 1], s=10, facecolors="none",
                         edgecolors="#2980b9", linewidths=0.3, alpha=0.4)

        # --- centroid TRAILS (fading) over the last trail_k frames ---
        lo = max(0, i - trail_k)
        for c in f["cell_ids"]:
            c = int(c)
            xs, ys = [], []
            for j in range(lo, i + 1):
                fj = frames[j]
                idsj = list(fj["cell_ids"])
                if c in idsj:
                    p = fj["cent_um"][idsj.index(c)]
                    xs.append(p[0]); ys.append(p[1])
            if len(xs) >= 2:
                for k in range(len(xs) - 1):
                    a = 0.15 + 0.75 * (k / max(1, len(xs) - 1))
                    axxy.plot(xs[k:k + 2], ys[k:k + 2], "-", color="#34495e",
                              lw=1.4, alpha=a, zorder=2)

        # --- current centroids coloured by state, switched ringed ---
        cent = f["cent_um"]; cst = f["cent_state"]; csw = f["cent_switched"]
        for sval, col in _STATE_COLOR.items():
            m = cst == sval
            if m.any():
                axxy.scatter(cent[m, 0], cent[m, 1], s=120, c=col, alpha=0.95,
                             edgecolors="k", linewidths=0.6, zorder=4,
                             label=_STATE_LABEL[sval])
                axxz.scatter(cent[m, 0], cent[m, 2], s=120, c=col, alpha=0.95,
                             edgecolors="k", linewidths=0.6, zorder=4)
        if csw.any():
            axxy.scatter(cent[csw, 0], cent[csw, 1], s=240, facecolors="none",
                         edgecolors="#2980b9", linewidths=2.0, zorder=5,
                         label="junction-switched")

        # --- TRACTION ARROWS on the rim cells (active pull direction+magnitude) ---
        ids = list(f["cell_ids"])
        for c, F in f["tractions"].items():
            if c not in ids:
                continue
            p = f["cent_um"][ids.index(c)]
            mag = np.linalg.norm(F)
            if mag <= 0:
                continue
            scale = arrow_um * (mag / f_ref)
            ux, uy = F[0] / mag, F[1] / mag
            axxy.annotate("", xy=(p[0] + ux * scale, p[1] + uy * scale),
                          xytext=(p[0], p[1]),
                          arrowprops=dict(arrowstyle="-|>", color="#c0392b",
                                          lw=2.0, alpha=0.9), zorder=6)

        axxz.axhline(f["z0_um"], color="#7f8c8d", lw=2, ls="--", alpha=0.6)
        axxy.set_xlim(*xlim); axxy.set_ylim(*ylim); axxy.set_aspect("equal")
        axxy.set_title(f"top-down: centroid trails + active-pull arrows   "
                       f"A/A0={f['AoverA0']:.2f}")
        axxy.set_xlabel("x [µm]"); axxy.set_ylabel("y [µm]")
        axxz.set_xlim(*xlim); axxz.set_ylim(*zlim)
        axxz.set_title("side (z = crawl onto substrate)")
        axxz.set_xlabel("x [µm]"); axxz.set_ylabel("z [µm]")
        axxy.legend(loc="upper right", fontsize=7, framealpha=0.9)
        # per-frame real-time-equivalent timestamp (t_sim → t_real via accel S)
        rt_f = map_realtime(f["t"], f.get("dt", _dt))
        fig.text(0.5, 0.885,
                 f"{f['phase']}   step {f['t']}   "
                 f"t={rt_f.t_sim_human} sim ≈ {rt_f.t_real_human} real   "
                 f"N={f['n_cells']} cells   "
                 f"| rim centroids migrated outward {mean_dr:+.2f} µm (mean), "
                 f"{max_dr:+.2f} µm (max)",
                 ha="center", fontsize=9.5)
        return []

    anim = FuncAnimation(fig, draw, frames=len(frames), blit=False)
    writer = FFMpegWriter(fps=6, bitrate=2600)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    anim.save(out_path, writer=writer, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}  ({len(frames)} frames)")
    print(f"  rim-cell radial centroid migration: mean {mean_dr:+.3f} µm, "
          f"max {max_dr:+.3f} µm")
    return mean_dr, max_dr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--div-every", type=int, default=20000,
                    help="slow division cadence (steps); large => traction-driven")
    ap.add_argument("--out",
                    default="ffn_sim/outputs/h_dcm_active/figs/active_spheroid_migration.mp4")
    args = ap.parse_args()

    p = ResolvedActiveSpheroid()
    if args.quick:
        n_active, n_max, agg, spr = 10, 16, 3, 12
    else:
        n_active, n_max, agg, spr = 14, 26, 6, 24
    print(f"[migration] capturing active spheroid (n_active={n_active}, "
          f"div_every={args.div_every}) ...")
    frames, _ = capture(p, n_active=n_active, n_max=n_max,
                        aggreg_blocks=agg, spread_blocks=spr,
                        div_every=args.div_every)
    print(f"[migration] {len(frames)} frames captured -> rendering MP4")
    render(frames, args.out)


if __name__ == "__main__":
    main()
