"""Render the spreading PROCESS as an animated GIF from captured morphology frames.

Reads outputs/layer2/morphology_frames.npz (written by layer2_morphology_vis.py) and animates the
spheroid: TOP-DOWN (xy footprint) + SIDE (xz cap on the dish), cells coloured by height. This is the
actual spatial process the A/A0 curve summarizes — you watch the footprint grow only modestly while
the aggregate stays a 3D cap (the §C–§F structural-limit, animated). Local render (no sim).

Usage:  python -m ffn_sim.scripts.layer2_morphology_anim [max_frames=50] [fps=8]
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from ffn_sim.archive.hoomd_legacy.spheroid.observables import core_projected_area, effective_radius

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "layer2"
_NPZ = _OUT / "morphology_frames.npz"
_GIF = _OUT / "figs" / "fig_layer2_morphology_anim.gif"


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    max_frames = int(args[0]) if len(args) > 0 else 50
    fps = int(args[1]) if len(args) > 1 else 8
    if not _NPZ.exists():
        print(f"no frames at {_NPZ} — run layer2_morphology_vis first", file=sys.stderr); return 1
    z = np.load(_NPZ, allow_pickle=True)
    times = z["times"]; a0c = float(z["a0_core"]); z_sub = float(z["z_sub"]); r0 = float(z["r0"])
    cond = str(z["cond"])
    n = len(times)
    frames = [z[f"f{i}"] for i in range(n)]
    idx = np.linspace(0, n - 1, min(max_frames, n)).round().astype(int)
    um = 1e6

    pf = frames[idx[-1]]
    ext = np.abs(pf[:, :2] - pf[:, :2].mean(axis=0)).max() * um * 1.12
    zmax = (pf[:, 2].max() - z_sub) * um * 1.12
    rad = r0 * um / 2.0
    s = max((rad ** 2) * 0.45, 0.6)

    fig, (axt, axs) = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle(f"Layer-2 spheroid spreading — the actual process (condition {cond}, "
                 f"R0≈{effective_radius(a0c)*um:.0f} µm)\nTOP-DOWN footprint (left) grows only modestly; "
                 f"SIDE view (right) stays a 3D CAP on the dish — not a flat monolayer", fontsize=11)

    def draw(fi):
        axt.clear(); axs.clear()
        p = frames[fi]
        com = p[:, :2].mean(axis=0)
        x = (p[:, 0]-com[0])*um; y = (p[:, 1]-com[1])*um; zz = (p[:, 2]-z_sub)*um
        aa = core_projected_area(np.asarray(p), 1.6*r0)/a0c
        axt.scatter(x, y, s=s, c=zz, cmap="coolwarm", edgecolors="none", alpha=0.85, vmin=0, vmax=zmax)
        axt.set_xlim(-ext, ext); axt.set_ylim(-ext, ext); axt.set_aspect("equal")
        axt.set_title(f"TOP-DOWN (xy)   t={times[fi]/3600:.0f} h   A/A0={aa:.2f}   N={len(p)}", fontsize=10)
        axt.set_xlabel("x (µm)"); axt.set_ylabel("y (µm)"); axt.grid(alpha=0.2)
        axs.scatter(x, zz, s=s, c=zz, cmap="coolwarm", edgecolors="none", alpha=0.85, vmin=0, vmax=zmax)
        axs.axhline(0, color="saddlebrown", lw=2.5)
        axs.text(0.02, 0.93, "substrate (dish)", color="saddlebrown", transform=axs.transAxes, fontsize=9)
        axs.set_xlim(-ext, ext); axs.set_ylim(-rad, zmax); axs.set_aspect("equal")
        axs.set_title(f"SIDE (xz) — stays a 3D cap   height≈{zz.max():.0f} µm", fontsize=10)
        axs.set_xlabel("x (µm)"); axs.set_ylabel("z above dish (µm)"); axs.grid(alpha=0.2)

    anim = FuncAnimation(fig, draw, frames=idx, interval=1000 // fps)
    _GIF.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(_GIF), writer=PillowWriter(fps=fps))
    print(f"wrote {_GIF} ({len(idx)} frames @ {fps} fps)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
