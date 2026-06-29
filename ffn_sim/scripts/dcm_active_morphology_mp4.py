"""Morphology MP4 of the ACTIVE-spreading DCM spheroid (ball -> compact -> spread).

Re-runs the active spheroid (active rim traction + bulk-pressure junction switch + live
division on the native-mesh + bilinear-tent stack) capturing PER-FRAME mesh-node positions
+ per-cell state, and renders an MP4 with two synchronized views:

  * top-down (x-y): the basal FOOTPRINT spreading outward,
  * side (x-z):     the ball flattening onto the substrate.

Mesh nodes are coloured by their cell's 3-zone state (proliferating / quiescent / necrotic);
junction-switched cells (cadherin->integrin) are ringed. Title shows t, N cells, A/A0.

  ~/miniconda3/envs/ffn_sim/bin/python ffn_sim/scripts/dcm_active_morphology_mp4.py [--quick]

Output: ffn_sim/outputs/h_dcm_active/figs/active_spheroid_morphology.mp4
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


def capture(p: ResolvedActiveSpheroid, *, n_active, n_max, dt=3.0e-10,
            aggreg_blocks=6, spread_blocks=24, block=1000):
    """Build + run the active spheroid, returning a list of per-frame snapshots."""
    h = build_active_spheroid(p, n_active=n_active, n_max=n_max,
                              integrin_substrate=True, belt=True)
    sim = h["sim"]; st = h["st"]; active = h["active"]
    cell_of_node = h["cell_of_node"]
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
                           block=block)
        raise RuntimeError("non-finite at finite gate")

    frames = []

    def snap(phase):
        pg = _positions(sim)
        con = cell_of_node                       # node -> cell id (>=0 active)
        mask = con >= 0
        A = _footprint_area(pg, mask, z0, R)
        node_cell = con[mask]
        node_state = st.state[node_cell]
        node_switched = st.switched[node_cell]
        frames.append(dict(
            t=int(sim.timestep), phase=phase,
            xyz_um=(pg[mask] - [0, 0, 0]) * _UM,
            state=node_state.copy(), switched=node_switched.copy(),
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
        if b % 4 == 0 or b == spread_blocks - 1:
            prolif.act(int(sim.timestep))
        snap("spreading")
    return frames


def render(frames, out_path):
    xy_all = np.vstack([f["xyz_um"] for f in frames])
    xlim = (xy_all[:, 0].min() - 5, xy_all[:, 0].max() + 5)
    ylim = (xy_all[:, 1].min() - 5, xy_all[:, 1].max() + 5)
    zlim = (min(xy_all[:, 2].min(), frames[0]["z0_um"]) - 3, xy_all[:, 2].max() + 5)

    fig, (axxy, axxz) = plt.subplots(1, 2, figsize=(13, 6.2))
    fig.suptitle("ACTIVE-spreading DCM spheroid — morphology (ball -> compact -> spread)",
                 fontweight="bold")

    def draw(i):
        f = frames[i]
        for ax in (axxy, axxz):
            ax.clear()
        xyz = f["xyz_um"]; stt = f["state"]; sw = f["switched"]
        for sval, col in _STATE_COLOR.items():
            m = stt == sval
            if m.any():
                axxy.scatter(xyz[m, 0], xyz[m, 1], s=14, c=col, alpha=0.7,
                             edgecolors="none", label=_STATE_LABEL[sval])
                axxz.scatter(xyz[m, 0], xyz[m, 2], s=14, c=col, alpha=0.7,
                             edgecolors="none")
        if sw.any():
            axxy.scatter(xyz[sw, 0], xyz[sw, 1], s=30, facecolors="none",
                         edgecolors="#c0392b", linewidths=0.6,
                         label="junction-switched")
        # substrate line (side view)
        axxz.axhline(f["z0_um"], color="#7f8c8d", lw=2, ls="--", alpha=0.6)
        axxy.set_xlim(*xlim); axxy.set_ylim(*ylim); axxy.set_aspect("equal")
        axxy.set_title(f"top-down (basal footprint)   A/A0={f['AoverA0']:.2f}")
        axxy.set_xlabel("x [µm]"); axxy.set_ylabel("y [µm]")
        axxz.set_xlim(*xlim); axxz.set_ylim(*zlim)
        axxz.set_title(f"side (z = spread on substrate)")
        axxz.set_xlabel("x [µm]"); axxz.set_ylabel("z [µm]")
        axxy.legend(loc="upper right", fontsize=7, framealpha=0.9)
        fig.text(0.5, 0.94, f"{f['phase']}   step {f['t']}   N={f['n_cells']} cells",
                 ha="center", fontsize=10)
        return []

    anim = FuncAnimation(fig, draw, frames=len(frames), blit=False)
    writer = FFMpegWriter(fps=6, bitrate=2400)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    anim.save(out_path, writer=writer, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}  ({len(frames)} frames)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out",
                    default="ffn_sim/outputs/h_dcm_active/figs/active_spheroid_morphology.mp4")
    args = ap.parse_args()

    p = ResolvedActiveSpheroid()
    if args.quick:
        n_active, n_max, agg, spr = 10, 16, 3, 10
    else:
        n_active, n_max, agg, spr = 14, 26, 6, 24
    print(f"[morphology] capturing active spheroid (n_active={n_active}) ...")
    frames = capture(p, n_active=n_active, n_max=n_max,
                     aggreg_blocks=agg, spread_blocks=spr)
    print(f"[morphology] {len(frames)} frames captured -> rendering MP4")
    render(frames, args.out)


if __name__ == "__main__":
    main()
