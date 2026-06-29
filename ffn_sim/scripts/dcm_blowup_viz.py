"""Visualize a dt-instability spreading failure mode (PI 2026-06-12).

The dt-stability sweep showed two DISTINCT large-dt failure modes for N=100
spreading from the converged aggregate (the stable dt is 1e-9):

  * dt=1e-7 (100×) — BLOW-UP: V/V0 runs away to ~1400× (cells balloon to the box
    ceiling). BAOAB integrator unstable.
  * dt=1e-8 (10×)  — COLLAPSE: no explosion, but the turgor volume is not conserved
    (V/V0 drifts to ~0.77) so the aggregate slowly IMPLODES (A/A0 1.0→1.29→0.25,
    Rg 37→14 µm). A large-dt accuracy artifact, not real physics.

This script REPRODUCES one failure mode (the instability is device-independent)
WITH the position frames the sweep pkl does not store, and renders it:

  row 1 — cell node clouds (coloured by cell) at 4 timepoints, per-panel auto-limits
          + extent label so the spatial scale (growing for blow-up, shrinking for
          collapse) is visible.
  row 2 — diagnostics vs step: V/V0 (log; the stable band V/V0≈1 dashed), Rg, maxZ.

Capture (gbook GPU, fast) / render (Mac, has ffmpeg) can be split:
  gbook: python dcm_blowup_viz.py --dt 1e-8 --steps 45000 --capture-pkl F.pkl --no-render
  Mac:   python dcm_blowup_viz.py --render-pkl F.pkl --label "10x (1e-8) COLLAPSE" --tag dt1e-8
Or all-in-one (default): python dcm_blowup_viz.py --dt 1e-7 --steps 7000 --tag dt1e-7
"""

from __future__ import annotations

import argparse
import dataclasses
import pickle
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

OUT = Path("ffn_sim/outputs/h_dcm_two_stage")
FIGS = OUT / "figs"
UM = 1.0e6
S_ACCEL = 6.0e5


def _ext_str(pos):
    return (f"{np.ptp(pos[:, 0]) * UM:.0f}×{np.ptp(pos[:, 1]) * UM:.0f}×"
            f"{np.ptp(pos[:, 2]) * UM:.0f}µm")


def _scatter(ax, pos, ranges, title):
    cmap = plt.get_cmap("tab20")
    col = np.zeros((pos.shape[0], 4))
    for ci, (lo, hi) in enumerate(ranges):
        col[lo:hi] = cmap(ci % 20)
    ax.scatter(pos[:, 0] * UM, pos[:, 1] * UM, pos[:, 2] * UM, c=col, s=2,
               depthshade=True)
    ax.set_title(title, fontsize=8)
    ax.set_xlabel("x µm", fontsize=6); ax.set_ylabel("y µm", fontsize=6)
    ax.tick_params(labelsize=5)


def capture(dt, steps, frames):
    """Reproduce the spreading run at `dt`, return frames + diagnostics."""
    from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import pick_device
    from ffn_sim.scripts.dcm_two_stage_production import spread, centroids
    d = pickle.load(open(OUT / "two_stage_n100.pkl", "rb"))
    R = d["R_cell"]; V0 = d["V0"]; z0 = d["z0"]; n_cells = int(d["n_cells"])
    ranges = d["ranges"]; tris0 = d["tris0"]
    agg_pos = d["agg"]["frames"][-1].copy()
    agg_pos[:, 2] -= (agg_pos[:, 2].min() - z0)
    cxy = centroids(agg_pos, ranges).mean(0)[:2]
    agg_pos[:, 0] -= cxy[0]; agg_pos[:, 1] -= cxy[1]
    dev = pick_device(None)
    p = dataclasses.replace(d["p"], dt=dt)
    print(f"[capture] dt={dt:.0e} steps={steps} frames={frames} "
          f"dev={type(dev).__name__}", flush=True)
    s2, _h = spread(p, n_cells, dev=dev, init_pos=agg_pos.copy(), steps=steps,
                    frames=frames, R=R, z0=z0, V0=V0, tris0=tris0)
    diags = s2["diags"]
    return dict(
        dt=dt, ranges=ranges, frames=s2["frames"], steps=s2["steps"],
        A0_top=s2["A0_topdown_um2"],
        vv0=np.array([g["VV0_mean"] for g in diags]),
        maxz=np.array([g["maxZ_um"] for g in diags]),
        rg=np.array([g["Rg_um"] for g in diags]),
        aa0=np.array([g["topdown_um2"] / s2["A0_topdown_um2"] for g in diags]))


def render(cap, label, tag, mp4=True):
    frames = cap["frames"]; steps = cap["steps"]; ranges = cap["ranges"]
    vv0 = cap["vv0"]; maxz = cap["maxz"]; rg = cap["rg"]; aa0 = cap["aa0"]
    dt = cap["dt"]
    print(f"  {label}: V/V0 {vv0.min():.2f}–{vv0.max():.1f}, Rg {rg.max():.0f}→"
          f"{rg.min():.0f}µm, A/A0 {aa0[0]:.2f}→{aa0[-1]:.2f}", flush=True)

    onset = next((i for i in range(1, len(vv0))
                  if abs(np.log(max(vv0[i], 1e-9))) > 0.15), len(frames) // 3)
    picks = sorted(set([0, onset, (onset + len(frames) - 1) // 2, len(frames) - 1]))
    while len(picks) < 4:
        picks.append(min(len(frames) - 1, picks[-1] + 1))
    picks = picks[:4]

    fig = plt.figure(figsize=(15, 8))
    fig.suptitle(f"dt={dt:.0e} {label} — N=100 spreading reproduced", fontsize=12)
    for k, fi in enumerate(picks):
        ax = fig.add_subplot(2, 4, k + 1, projection="3d")
        _scatter(ax, frames[fi], ranges,
                 f"step {steps[fi]}  A/A0={aa0[fi]:.2f}\n"
                 f"V/V0={vv0[fi]:.2f}  Rg={rg[fi]:.0f}µm\nextent~{_ext_str(frames[fi])}")
    axv = fig.add_subplot(2, 1, 2)
    axv.semilogy(steps, vv0, "-o", color="crimson", ms=3, label="V/V0")
    axv.axhline(1.0, color="green", ls="--", lw=1, label="V/V0=1 (stable / dt=1e-9)")
    axv.set_xlabel(f"spreading step (dt={dt:.0e})")
    axv.set_ylabel("V/V0 (log)", color="crimson")
    axv.tick_params(axis="y", labelcolor="crimson"); axv.legend(loc="upper left", fontsize=8)
    axr = axv.twinx()
    axr.plot(steps, rg, "-s", color="navy", ms=3, label="Rg µm")
    axr.plot(steps, maxz, "-^", color="teal", ms=3, label="maxZ µm")
    axr.set_ylabel("Rg / maxZ µm", color="navy")
    axr.tick_params(axis="y", labelcolor="navy"); axr.legend(loc="upper right", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    png = FIGS / f"blowup_{tag}.png"
    fig.savefig(png, dpi=110); plt.close(fig)
    print(f"  wrote {png}", flush=True)

    if mp4:
        figA = plt.figure(figsize=(6, 6)); axA = figA.add_subplot(111, projection="3d")
        def _draw(i):
            axA.clear()
            _scatter(axA, frames[i], ranges, f"dt={dt:.0e} {label}  step {steps[i]}\n"
                     f"A/A0={aa0[i]:.2f}  V/V0={vv0[i]:.2f}  Rg={rg[i]:.0f}µm")
        anim = animation.FuncAnimation(figA, _draw, frames=len(frames), interval=200)
        mp4p = FIGS / f"blowup_{tag}.mp4"
        try:
            anim.save(mp4p, writer="ffmpeg", fps=5, dpi=90)
            print(f"  wrote {mp4p}", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"  (mp4 skipped: {e})", flush=True)
        plt.close(figA)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dt", type=float, default=1.0e-7)
    ap.add_argument("--steps", type=int, default=7000)
    ap.add_argument("--frames", type=int, default=22)
    ap.add_argument("--label", default="BLOW-UP")
    ap.add_argument("--tag", default="dt1e-7")
    ap.add_argument("--capture-pkl", default="")
    ap.add_argument("--render-pkl", default="")
    ap.add_argument("--no-render", action="store_true")
    args = ap.parse_args()

    if args.render_pkl:
        cap = pickle.load(open(args.render_pkl, "rb"))
        render(cap, args.label, args.tag)
        return
    cap = capture(args.dt, args.steps, args.frames)
    if args.capture_pkl:
        with open(args.capture_pkl, "wb") as fh:
            pickle.dump(cap, fh)
        print(f"  saved frames → {args.capture_pkl}", flush=True)
    if not args.no_render:
        render(cap, args.label, args.tag)


if __name__ == "__main__":
    main()
