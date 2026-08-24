"""Time-series MESH morphology (shaded triangulated surfaces, like dcm_spheroid_mesh_viz).

Renders a clean-ball spread pickle at several timepoints as proper per-cell polygon
meshes (Poly3DCollection, shaded by face normal), top-down + side on a common scale,
so the overshoot-then-relax is visible as actual cell shapes spreading out then
retracting/compacting. A/A0 + maxZ curves below mark the timepoints. Output → figs/.
"""
from __future__ import annotations
import argparse
import pickle
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

UM = 1e6
OUT = Path("aleph/outputs/h_dcm_two_stage")


def _mesh_panel(ax, frame, ranges, tris0, elev, azim, h, zmax, topdown):
    cmap = plt.get_cmap("tab20").colors
    nv = int(tris0.max()) + 1
    for c, (a, b) in enumerate(ranges):
        verts = frame[a:b] * UM
        if verts.shape[0] < nv:
            continue
        tri = verts[tris0]
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        nn = np.linalg.norm(n, axis=1, keepdims=True)
        sh = 0.55 + 0.45 * np.clip((n / np.where(nn > 0, nn, 1)) @ np.array([0.3, 0.3, 0.9]), 0, 1)
        base = np.array(cmap[c % 20])
        rgba = np.zeros((tris0.shape[0], 4)); rgba[:, :3] = base[None, :] * sh[:, None]; rgba[:, 3] = 0.92
        ax.add_collection3d(Poly3DCollection(tri, facecolors=rgba,
                            edgecolors=(0, 0, 0, 0.18), linewidths=0.12))
    ax.set_xlim(-h, h); ax.set_ylim(-h, h); ax.set_zlim(-2, zmax)
    ax.set_box_aspect((2 * h, 2 * h, 0.5 * h if topdown else (zmax + 2)))
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("x[µm]", fontsize=7); ax.set_ylabel("y[µm]", fontsize=7)
    ax.tick_params(labelsize=6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", required=True)
    ap.add_argument("--out", default="figs/cleanball_mesh_timeseries.png")
    ap.add_argument("--title", default="N=12 clean ball + wetting")
    args = ap.parse_args()
    with open(args.pkl, "rb") as fh:
        d = pickle.load(fh)
    s2 = d["spread"]; ranges = d["ranges"]; tris0 = np.asarray(d["tris0"])
    A0 = s2.get("A0_topdown_um2") or s2["diags"][0]["topdown_um2"]
    aa = np.array([g["topdown_um2"] / A0 for g in s2["diags"]])
    mz = np.array([g["maxZ_um"] for g in s2["diags"]])
    steps = np.array(s2["steps"]); frames = s2["frames"]
    nmem = int(max(b for (_, b) in ranges))

    ipk = int(np.argmax(aa))
    idxs = sorted(set([0, ipk // 2, ipk, (ipk + len(aa) - 1) // 2, len(aa) - 1]))
    recs, labels = [], []
    for ix in idxs:
        fr = np.asarray(frames[ix]).copy()
        c = fr[:nmem].mean(0); c[2] = 0.0
        recs.append(fr - c)
        ph = ("START" if ix == 0 else "PEAK" if ix == ipk else
              "rising" if ix < ipk else "END" if ix == len(aa) - 1 else "relaxing")
        labels.append(f"{ph} · step {int(steps[ix])}\nA/A0={aa[ix]:.2f}  maxZ={mz[ix]:.0f}µm")
    nc = len(recs)
    h = max(np.abs(r[:nmem, :2]).max() * UM for r in recs) * 1.04
    zmax = max(r[:nmem, 2].max() * UM for r in recs) * 1.08

    fig = plt.figure(figsize=(3.3 * nc, 10))
    gs = fig.add_gridspec(3, nc, height_ratios=[1.05, 1.05, 0.8])
    for k, (fr, lab) in enumerate(zip(recs, labels)):
        axt = fig.add_subplot(gs[0, k], projection="3d")
        _mesh_panel(axt, fr, ranges, tris0, 90, -90, h, zmax, True)
        axt.set_title(lab, fontsize=9)
        axs = fig.add_subplot(gs[1, k], projection="3d")
        _mesh_panel(axs, fr, ranges, tris0, 10, -72, h, zmax, False)
    fig.text(0.012, 0.80, "TOP-DOWN (xy)", rotation=90, va="center", fontsize=10)
    fig.text(0.012, 0.50, "SIDE (xz)", rotation=90, va="center", fontsize=10)

    axc = fig.add_subplot(gs[2, :])
    axc.plot(steps, aa, "o-", color="C3", lw=2, label="A/A0 (footprint)")
    axc.axhline(1.0, ls=":", color="0.6")
    for ix in idxs:
        axc.axvline(steps[ix], color="0.7", lw=0.8, ls="--")
    axc.plot(steps[ipk], aa[ipk], "*", color="k", ms=16)
    axc.set_xlabel("spread step"); axc.set_ylabel("A/A0", color="C3")
    axc.tick_params(axis="y", labelcolor="C3")
    ax2 = axc.twinx(); ax2.plot(steps, mz, "s-", color="C0", ms=3)
    ax2.set_ylabel("maxZ [µm]", color="C0"); ax2.tick_params(axis="y", labelcolor="C0")
    axc.set_title("RISE: wetting spreads the basal footprint (A/A0↑)   →   "
                  "RELAX: cohesion retracts + compacts (A/A0↓ AND maxZ↓)", fontsize=10)

    fig.suptitle(f"{args.title} — MESH time-series\nA/A0 {aa[0]:.2f} → peak {aa[ipk]:.2f} "
                 f"(step {int(steps[ipk])}) → {aa[-1]:.2f}: transient overshoot, then cohesion "
                 "pulls it back toward the compact-ball equilibrium", fontsize=12)
    fig.tight_layout(rect=(0.02, 0, 1, 0.93))
    p = OUT / args.out
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=135)
    print(f"wrote {p}  peak={aa[ipk]:.2f}@{int(steps[ipk])} end={aa[-1]:.2f} "
          f"maxZ {mz[0]:.0f}->{mz[-1]:.0f}")


if __name__ == "__main__":
    main()
