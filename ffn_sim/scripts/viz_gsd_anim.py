"""Animate a GSD trajectory: cortex shell + growing crosslink/myosin bonds.

Renders every frame as a rotating 3D view — faint actin shell, crosslink bonds
(red), myosin-attach bonds (blue) — annotated with per-frame bond counts and
(optionally) the matched sample's adv / s_grip parsed from the run log. Writes
mp4 (ffmpeg) or falls back to gif (pillow). For the v1500 native run the global
shape is ~static (r/r0≈1); what moves is the binding network assembling.

Usage:
    python ffn_sim/scripts/viz_gsd_anim.py <traj.gsd> <out.mp4|.gif> [run.log]
"""
from __future__ import annotations

import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import gsd.hoomd

PAT = re.compile(r"s_grip=([\d.eE+-]+)nm\s+engaged=(\d+)\s+adv=(\d+)")


def main():
    path, out = sys.argv[1], sys.argv[2]
    log = sys.argv[3] if len(sys.argv) > 3 else None

    samp = []
    if log:
        for line in open(log):
            m = PAT.search(line)
            if m:
                samp.append((float(m.group(1)), int(m.group(2)), int(m.group(3))))

    with gsd.hoomd.open(path, "r") as t:
        frames = [t[i] for i in range(len(t))]
    n = len(frames)
    bnames = list(frames[0].bonds.types)
    xl_ids = [i for i, nm in enumerate(bnames) if "xlink_attach" in nm]
    my_ids = [i for i, nm in enumerate(bnames) if "myosin_attach" in nm]
    tid0 = np.asarray(frames[0].particles.typeid)
    actin_tid = int(np.argmax(np.bincount(tid0)))

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")

    def draw(fi):
        ax.clear()
        fr = frames[fi]
        pos = np.asarray(fr.particles.position) * 1e6  # µm
        tid = np.asarray(fr.particles.typeid)
        bt = np.asarray(fr.bonds.typeid); bg = np.asarray(fr.bonds.group)
        act = tid == actin_tid
        pa = pos[act]
        st = max(1, pa.shape[0] // 25000)
        ax.scatter(pa[::st, 0], pa[::st, 1], pa[::st, 2], s=0.4, c="0.8",
                   alpha=0.25, linewidths=0)
        # crosslink bonds (red) + myosin-attach bonds (blue)
        xlb = bg[np.isin(bt, xl_ids)] if xl_ids else np.empty((0, 2), int)
        myb = bg[np.isin(bt, my_ids)] if my_ids else np.empty((0, 2), int)
        if len(xlb):
            ax.add_collection3d(Line3DCollection(pos[xlb], colors="red",
                                                 linewidths=2.5, alpha=1.0))
        if len(myb):
            ax.add_collection3d(Line3DCollection(pos[myb], colors="blue",
                                                 linewidths=1.2, alpha=0.7))
        lim = np.abs(pos).max()
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=20, azim=4 * fi)  # rotate
        extra = ""
        if fi < len(samp):
            sg, eng, adv = samp[fi]
            extra = f"  adv={adv}  s_grip={sg:.0f}nm  engaged={eng}"
        ax.set_title(f"frame {fi}/{n-1}   xlinks={len(xlb)}  myosin={len(myb)}{extra}",
                     fontsize=9)
        for a in ("x", "y", "z"):
            getattr(ax, f"set_{a}label")(f"{a} (µm)", fontsize=8)

    anim = animation.FuncAnimation(fig, draw, frames=n, interval=600)
    try:
        if out.endswith(".mp4"):
            anim.save(out, writer="ffmpeg", dpi=120, fps=2)
        else:
            anim.save(out, writer="pillow", dpi=110, fps=2)
    except Exception as exc:
        alt = out.rsplit(".", 1)[0] + ".gif"
        print(f"[writer fallback: {exc} -> {alt}]")
        anim.save(alt, writer="pillow", dpi=110, fps=2)
        out = alt
    print(f"wrote {out}  ({n} frames)")


if __name__ == "__main__":
    main()
