"""Visualise the TWO-STAGE DCM spheroid run (PI 2026-06-12) — MacBook side.

gbook computes (scripts/dcm_two_stage_production.py → two_stage_n{N}.pkl); the MacBook
renders. Shows BOTH stages as real cell SURFACES (not dots) + the physical diagnostics
the PI asked about:

  STAGE 1 · AGGREGATION — the loose placement COMPACTS + ROUNDS into a cohesive ball
    (Rg shrinks, asphericity → 0). Rendered as surfaces (initial loose vs final compact).
  STAGE 2 · SPREADING — the aggregate (started from stage 1's FINISHED ball) sediments,
    contacts and spreads; the lamellipodia (rim cell basal caps) are coloured red.
  DIAGNOSTICS — V/V0 (volume held? cells must keep volume while the SHAPE flattens —
    PI's question), Rg, max height, basal footprint over both stages.

Outputs (ffn_sim/outputs/h_dcm_two_stage/figs/):
  two_stage_n{N}.png          — combined still (agg loose→compact, spread, diagnostics)
  two_stage_n{N}_agg.mp4      — aggregation compacting
  two_stage_n{N}_spread.mp4   — spreading with lamellipodia

Run:  PYTHONPATH=. python ffn_sim/scripts/dcm_two_stage_viz.py --n 200
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FFMpegWriter  # noqa: E402
from matplotlib.colors import LightSource  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

from ffn_sim.scripts.dcm_lamellipodium_viz import (  # noqa: E402
    _lambert, draw_lamellipodium_surfaces, draw_substrate, lamellipodium_geometry,
    C_BODY, C_LEAD)

for cand in ("/opt/homebrew/bin/ffmpeg", "/usr/bin/ffmpeg", "/usr/local/bin/ffmpeg"):
    if Path(cand).exists():
        matplotlib.rcParams["animation.ffmpeg_path"] = cand
        break

OUT = Path("ffn_sim/outputs/h_dcm_two_stage/figs")
OUT.mkdir(parents=True, exist_ok=True)
UM = 1.0e6
_TAB = plt.get_cmap("tab20").colors


class _RimStub:
    """Minimal traction-like object so lamellipodium_geometry can recompute the rim
    (few neighbours + basal contact) from the saved rim_params."""

    def __init__(self, rp):
        self.R = rp["R"]; self.z0 = rp["z0"]
        self.r_neigh = rp["r_neigh"]; self.max_neigh = rp["max_neigh"]
        self.contact_band = rp["contact_band"]
        self.f_per_cell = {}; self.f_act = 1.2e-10


def draw_plain_surfaces(ax, pos, ranges, tris0, ls):
    """Aggregation render: every cell a shaded surface coloured by cell id."""
    from matplotlib.colors import to_rgb
    for c, (lo, hi) in enumerate(ranges):
        verts = pos[lo:hi]
        rgba = _lambert(verts, tris0, to_rgb(_TAB[c % len(_TAB)]), ls, 0.92)
        ax.add_collection3d(Poly3DCollection(verts[tris0] * UM, facecolors=rgba,
                            edgecolors=(0, 0, 0, 0.10), linewidths=0.1))


def _box(frames):
    allp = np.concatenate(frames, axis=0) * UM
    cx, cy = allp[:, 0].mean(), allp[:, 1].mean()
    r = 0.5 * max(np.ptp(allp[:, 0]), np.ptp(allp[:, 1])) * 1.15
    zlo, zhi = min(0.0, allp[:, 2].min()) - 2, allp[:, 2].max() * 1.1
    return (cx - r, cx + r), (cy - r, cy + r), (zlo, zhi)


def _setup3d(ax, xlim, ylim, zlim, elev, azim):
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_zlim(*zlim)
    ax.set_box_aspect((xlim[1] - xlim[0], ylim[1] - ylim[0], zlim[1] - zlim[0]))
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("x[µm]", fontsize=7); ax.set_ylabel("y[µm]", fontsize=7)
    ax.set_zlabel("z[µm]", fontsize=7)


def diagnostics_panel(ax, d, R):
    """V/V0, Rg, maxZ, footprint over both stages (agg then spread)."""
    a, s = d["agg"], d["spread"]
    na = len(a["diags"])
    xa = np.arange(na)
    xs = np.arange(len(s["diags"])) + na
    def col(stg, k): return [g[k] for g in stg["diags"]]
    ax.axvspan(-0.5, na - 0.5, color="0.93", zorder=0)
    ax.axvspan(na - 0.5, na + len(s["diags"]) - 0.5, color="0.86", zorder=0)
    ax.plot(xa, col(a, "VV0_mean"), "-o", ms=3, color="C0", label="V/V0 (agg)")
    ax.plot(xs, col(s, "VV0_mean"), "-o", ms=3, color="C0", mfc="w", label="V/V0 (spread)")
    ax.axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax.set_ylabel("V / V0  (volume held)", color="C0")
    ax.set_ylim(0.8, 1.3)
    ax.tick_params(axis="y", colors="C0")
    ax2 = ax.twinx()
    ax2.plot(xa, col(a, "Rg_um"), "-s", ms=3, color="C3", label="Rg (agg→shrinks)")
    ax2.plot(xs, col(s, "Rg_um"), "-s", ms=3, color="C3", mfc="w")
    ax2.plot(np.concatenate([xa, xs]),
             [g["footprint_um2"] ** 0.5 for g in a["diags"] + s["diags"]],
             "-^", ms=3, color="C2", label="√footprint")
    ax2.set_ylabel("Rg [µm] · √footprint [µm]", color="C3")
    ax.set_xlabel("frame  (grey=aggregation · darker=spreading)")
    ax.set_title("volume held (flat ~1) while Rg compacts then footprint spreads", fontsize=8)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=6, loc="upper center")


def render(d, n):
    R, z0 = d["R_cell"], d["z0"]
    ranges, tris0 = d["ranges"], d["tris0"]
    ls = LightSource(azdeg=225, altdeg=45)
    a, s = d["agg"], d["spread"]
    rim = _RimStub(d["rim_params"])

    # ---- combined still ----
    fig = plt.figure(figsize=(17, 5.4))
    # agg initial (loose)
    ax1 = fig.add_subplot(1, 4, 1, projection="3d")
    xb, yb, zb = _box(a["frames"])
    draw_plain_surfaces(ax1, a["frames"][0], ranges, tris0, ls)
    _setup3d(ax1, xb, yb, zb, 18, -60)
    ax1.set_title(f"STAGE 1 start — loose\nRg={a['diags'][0]['Rg_um']:.0f}µm", fontsize=8)
    # agg final (compact)
    ax2 = fig.add_subplot(1, 4, 2, projection="3d")
    draw_plain_surfaces(ax2, a["frames"][-1], ranges, tris0, ls)
    _setup3d(ax2, xb, yb, zb, 18, -60)
    ax2.set_title(f"STAGE 1 end — AGGREGATED\nRg={a['diags'][-1]['Rg_um']:.0f}µm "
                  f"asph={a['diags'][-1]['asphericity']:.2f}", fontsize=8)
    # spread final (lamellipodia)
    ax3 = fig.add_subplot(1, 4, 3, projection="3d")
    xs2, ys2, zs2 = _box(s["frames"])
    geo = lamellipodium_geometry(s["frames"][-1], ranges, rim)
    draw_lamellipodium_surfaces(ax3, s["frames"][-1], geo, ranges, tris0, ls)
    draw_substrate(ax3, xs2, ys2, z0 * UM)
    _setup3d(ax3, xs2, ys2, zs2, 14, -60)
    ax3.set_title(f"STAGE 2 — SPREAD on dish\n{len(geo['rim'])} lamellipodia (red), "
                  f"foot={s['diags'][-1]['footprint_um2']:.0f}µm²", fontsize=8)
    # diagnostics
    ax4 = fig.add_subplot(1, 4, 4)
    diagnostics_panel(ax4, d, R)
    fig.suptitle(f"DCM two-stage spheroid · N={n} cells · aggregation → spreading · "
                 f"{d['device']} · vol held V/V0≈{s['diags'][-1]['VV0_mean']:.2f} "
                 f"(cells keep volume, shape flattens)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    png = OUT / f"two_stage_n{n}.png"
    fig.savefig(png, dpi=140); plt.close(fig)
    print(f"wrote {png}", flush=True)

    # ---- MP4s ----
    _mp4_stage(a, ranges, tris0, ls, None, z0, OUT / f"two_stage_n{n}_agg.mp4",
               "AGGREGATION (cells come together)")
    _mp4_stage(s, ranges, tris0, ls, rim, z0, OUT / f"two_stage_n{n}_spread.mp4",
               "SPREADING (lamellipodia on dish)")


def _mp4_stage(stg, ranges, tris0, ls, rim, z0, path, title):
    xb, yb, zb = _box(stg["frames"])
    fig = plt.figure(figsize=(7.5, 6.8))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    writer = FFMpegWriter(fps=6, bitrate=2800)
    try:
        with writer.saving(fig, str(path), dpi=115):
            for pos, dg, st in zip(stg["frames"], stg["diags"], stg["steps"]):
                ax.clear()
                if rim is None:
                    draw_plain_surfaces(ax, pos, ranges, tris0, ls)
                else:
                    geo = lamellipodium_geometry(pos, ranges, rim)
                    draw_lamellipodium_surfaces(ax, pos, geo, ranges, tris0, ls)
                    draw_substrate(ax, xb, yb, z0 * UM)
                _setup3d(ax, xb, yb, zb, 16, -60)
                ax.set_title(f"{title}\nstep {st} · Rg={dg['Rg_um']:.0f}µm · "
                             f"V/V0={dg['VV0_mean']:.2f} · foot={dg['footprint_um2']:.0f}µm²",
                             fontsize=9)
                writer.grab_frame()
        print(f"wrote {path}", flush=True)
    except Exception as e:  # noqa: BLE001
        print(f"[warn] MP4 {path.name} skipped: {e}", flush=True)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--pkl", default=None, help="explicit pickle path (else by --n)")
    args = ap.parse_args()
    path = Path(args.pkl) if args.pkl else \
        Path(f"ffn_sim/outputs/h_dcm_two_stage/two_stage_n{args.n}.pkl")
    d = pickle.load(open(path, "rb"))
    print(f"loaded {path} (N={d['n_cells']}, {d['device']})", flush=True)
    render(d, d["n_cells"])


if __name__ == "__main__":
    main()
