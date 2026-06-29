#!/usr/bin/env python
"""H.3 cortex structure visualisation — dense PNG + multi-angle dynamics MP4s.

Renders the actual cell-cortex structure from a production trajectory
(``outputs/h3/production/lp_{scale}_{device}.npz``, frames (S,F,N,3)):

    fig_h3_cortex_structure.png         all F filaments on the R=10 μm shell,
                                        bending-energy coloured.
    h3_cortex_dynamics_equatorial_x.mp4 time evolution, camera fixed at +x.
    h3_cortex_dynamics_equatorial_y.mp4 time evolution, camera fixed at +y.
    h3_cortex_dynamics_top.mp4          time evolution, top-down view.
    h3_cortex_dynamics_perspective.mp4  time evolution, oblique perspective.

Visualization integrity (CLAUDE.md): SI/μm axes, no truncation, a colour
scale annotated in physical units (E_bend/kT), the R_cell shell radius in
the title. MP4 via matplotlib FFMpegWriter (ffmpeg). Per PI 2026-05-29: no
rotating-camera MP4 (mixes camera motion with time evolution); each
viewpoint gets its own stable-camera MP4 so frame-to-frame comparison is
unambiguous. The earlier ``h3_cortex_rotate.mp4`` (360° orbit) and the
in-frame orbit on ``h3_cortex_dynamics.mp4`` are both retired.

Usage:
    python ffn_sim/scripts/h3_cortex_anim.py --scale full --device gpu
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"
PROD = PKG / "outputs" / "h3" / "production"
FIGS = PKG / "outputs" / "h3" / "figs"


def _bending_energy_per_bond(frame: np.ndarray, k_theta: float, kT: float) -> np.ndarray:
    """Mean per-filament bending energy E_bend/kT mapped to each segment.

    frame: (F, N, 3). Returns (F, N-1) segment colour values (a filament's
    interior-angle energies broadcast to its bonds).
    """
    F, N, _ = frame.shape
    bv = frame[:, 1:, :] - frame[:, :-1, :]
    bn = bv / np.linalg.norm(bv, axis=-1, keepdims=True)
    cos = -np.einsum("...i,...i->...", bn[:, :-1, :], bn[:, 1:, :])  # (F,N-2)
    th = np.arccos(np.clip(cos, -1.0, 1.0))
    E = 0.5 * k_theta * (th - np.pi) ** 2 / kT                        # (F,N-2)
    seg = np.zeros((F, N - 1))
    seg[:, 1:-1] = 0.5 * (E[:, :-1] + E[:, 1:]) if N > 3 else 0.0
    if N >= 3:
        seg[:, 0] = E[:, 0]; seg[:, -1] = E[:, -1]
    return seg


def _segments(frame_um: np.ndarray) -> np.ndarray:
    """(F,N,3) μm → (F*(N-1), 2, 3) line segments for Line3DCollection."""
    F, N, _ = frame_um.shape
    a = frame_um[:, :-1, :].reshape(-1, 3)
    b = frame_um[:, 1:, :].reshape(-1, 3)
    return np.stack([a, b], axis=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scale", default="full")
    ap.add_argument("--device", default="gpu")
    ap.add_argument("--n-filaments", type=int, default=0,
                    help="subset for MP4 (0 = all; dynamics MP4 caps at 400).")
    ap.add_argument("--fps", type=int, default=24)
    args = ap.parse_args()

    with open(CFG) as f:
        cfg = yaml.safe_load(f)
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    R = p.R_cell * 1e6
    k_theta, kT = p.angle_k, p.kT

    frames = np.load(PROD / f"lp_{args.scale}_{args.device}.npz")["frames"]
    S, F, N, _ = frames.shape
    FIGS.mkdir(parents=True, exist_ok=True)
    lim = 1.15 * R

    # colour scale fixed across all renders (98th pct of last frame)
    Elast = _bending_energy_per_bond(frames[-1], k_theta, kT)
    vmax = float(np.percentile(Elast, 98)) or 1.0

    def draw(ax, frame, title):
        ax.clear()
        seg = _segments(frame * 1e6)
        col = _bending_energy_per_bond(frame, k_theta, kT).ravel()
        lc = Line3DCollection(seg, cmap="inferno", linewidths=0.8)
        lc.set_array(col); lc.set_clim(0, vmax)
        ax.add_collection3d(lc)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
        ax.set_xlabel("x [μm]"); ax.set_ylabel("y [μm]"); ax.set_zlabel("z [μm]")
        ax.set_title(title)
        return lc

    # --- PNG: dense full-cortex structure (last/equilibrated frame) ---
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    lc = draw(ax, frames[-1],
              f"H.3 cortex structure — {F} filaments × {N} beads on R={R:.0f} μm shell")
    ax.view_init(elev=18, azim=35)
    cb = fig.colorbar(lc, ax=ax, shrink=0.6, pad=0.1)
    cb.set_label(r"$E_{bend}/k_BT$ (bending / curvature)")
    f_png = FIGS / "fig_h3_cortex_structure.png"
    fig.savefig(f_png, dpi=150, bbox_inches="tight"); plt.close(fig)

    # --- MP4s: time evolution (thermal internal dynamics) from N fixed
    # viewpoints, one MP4 per view. Per PI 2026-05-29: no rotating-camera
    # MP4 (mixes camera motion with time evolution in the same frame).
    # The 360° orbit MP4 from earlier versions of this script is retired.
    nf = F if args.n_filaments <= 0 else min(args.n_filaments, F)
    nf = min(nf, 400)                       # cap for redraw cost
    sel = np.random.default_rng(0).choice(F, size=nf, replace=False)
    sub = frames[:, sel]                    # (S, nf, N, 3)

    cortex_views = [
        ("equatorial_x", "equatorial — +x axis", dict(elev=0, azim=0)),
        ("equatorial_y", "equatorial — +y axis", dict(elev=0, azim=90)),
        ("top",          "top (−z axis)",        dict(elev=88, azim=-90)),
        ("perspective",  "perspective",          dict(elev=22, azim=35)),
    ]
    mp4_paths = []
    for slug, title_suffix, view in cortex_views:
        writer = FFMpegWriter(fps=args.fps, bitrate=3200)
        fig = plt.figure(figsize=(8, 7))
        ax = fig.add_subplot(111, projection="3d")
        out_mp4 = FIGS / f"h3_cortex_dynamics_{slug}.mp4"
        with writer.saving(fig, str(out_mp4), dpi=120):
            for k in range(S):
                draw(ax, sub[k],
                     f"H.3 cortex — {title_suffix} — snapshot {k + 1}/{S} "
                     f"({nf} filaments)")
                ax.view_init(**view)
                writer.grab_frame()
        plt.close(fig)
        mp4_paths.append(out_mp4)
        print(f"WROTE_MP4 {out_mp4}", flush=True)

    print(f"VIZ_WRITTEN {f_png.name} " +
          " ".join(p.name for p in mp4_paths) +
          f" F={F} N={N} S={S} R={R:.1f}um vmax_EkT={vmax:.2f}", flush=True)


if __name__ == "__main__":
    main()
