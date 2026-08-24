"""Render the H.7 single-cell spreading / cell-movement animation.

Reads the per-chunk basal-xy frames npz written by `h7_spreading_compare.py
--frames-out` and renders a side-by-side animated GIF of the two lamellipodium
geometries spreading over time (the basal footprint = the cell's contact with
the substrate, growing/moving chunk by chunk). This is the "cell moving"
visualization (PI 2026-06-07).

Usage:
    python -m aleph.scripts.h7_spreading_anim --frames <frames.npz> --out <anim.gif>
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

_COLORS = {"basal_ring": "#2f6fb0", "polarized_patch": "#e23b3b"}


def render(frames_npz: str, out_gif: str, fps: int = 4) -> str:
    z = np.load(frames_npz, allow_pickle=False)
    geoms = [g if isinstance(g, str) else g.decode() for g in z["geoms"]]
    n_frames = {g: int(n) for g, n in zip(geoms, z["n_frames"])}
    nf = min(n_frames.values())
    frames = {g: [z[f"{g}__{i}"] for i in range(n_frames[g])] for g in geoms}

    # Common axis limits across all frames/geometries.
    allxy = np.vstack([f for g in geoms for f in frames[g] if len(f)])
    lim = float(np.abs(allxy).max()) * 1.08 if len(allxy) else 8.0

    fig, axes = plt.subplots(1, len(geoms), figsize=(5.2 * len(geoms), 5.4),
                             constrained_layout=True)
    if len(geoms) == 1:
        axes = [axes]

    def draw(k: int):
        for ax, g in zip(axes, geoms):
            ax.clear()
            xy = frames[g][k]
            if len(xy):
                ax.scatter(xy[:, 0], xy[:, 1], s=7, c=_COLORS.get(g, "0.4"),
                           alpha=0.7, edgecolors="none")
            ax.set_aspect("equal")
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
            n0 = len(frames[g][0]) or 1
            aa = (len(xy) / n0)
            ax.set_title(f"{g}\nbasal footprint, chunk {k}/{nf - 1} (N={len(xy)})",
                         fontsize=10)
        fig.suptitle("H.7 single-cell lamellipodium spreading (basal footprint over time) "
                     "— isotropic rim vs polarized patch", fontsize=12, fontweight="bold")

    anim = FuncAnimation(fig, draw, frames=nf, interval=1000 // max(1, fps))
    Path(out_gif).parent.mkdir(parents=True, exist_ok=True)
    anim.save(out_gif, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out_gif


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--frames", required=True, help="frames npz from h7_spreading_compare --frames-out")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1]
                    / "outputs" / "h7" / "figs" / "h7_spreading_movement.gif"))
    ap.add_argument("--fps", type=int, default=4)
    args = ap.parse_args()
    out = render(args.frames, args.out, fps=args.fps)
    print(f"[h7-anim] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
