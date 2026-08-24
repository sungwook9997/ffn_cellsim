"""Time-series morphology: why A/A0 overshoots then relaxes (clean-ball wetting spread).

Reads one clean-ball spread pickle and renders the spheroid at several timepoints
(top-down xy + side xz on a COMMON scale) across the trajectory, with the A/A0 and
maxZ curves below and the timepoints marked. Shows the two phases: (1) wetting spreads
the freshly-placed (not-yet-cohered) cells outward at ~constant height → A/A0 rises;
(2) cell-cell cohesion engages and pulls them back together → footprint AND height
shrink → A/A0 relaxes toward the compact-ball equilibrium. Output → figs/.
"""
from __future__ import annotations
import argparse
import pickle
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

UM = 1e6
OUT = Path("aleph/outputs/h_dcm_two_stage")


def _panel(ax, frame, ranges, axes, h, zlim):
    from scipy.spatial import ConvexHull
    i, j = axes
    cmap = plt.get_cmap("tab20").colors
    allp = []
    for c, (a, b) in enumerate(ranges):
        v = frame[a:b] * UM
        if v.shape[0] < 3:
            continue
        allp.append(v[:, [i, j]])
        ax.scatter(v[:, i], v[:, j], s=1.0, color=cmap[c % 20], alpha=0.5, linewidths=0)
    if allp:
        P = np.concatenate(allp)
        try:
            hh = ConvexHull(P); poly = P[hh.vertices]
            ax.fill(poly[:, 0], poly[:, 1], facecolor="0.75", edgecolor="k",
                    lw=1.3, ls="-", alpha=0.35, zorder=0)
        except Exception:
            pass
    ax.set_xlim(-h, h); ax.set_aspect("equal")
    ax.set_ylim(*((-h, h) if j == 1 else zlim))
    ax.tick_params(labelsize=7)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", required=True)
    ap.add_argument("--out", default="figs/cleanball_timeseries.png")
    ap.add_argument("--title", default="N=12 clean ball + wetting")
    args = ap.parse_args()
    with open(args.pkl, "rb") as fh:
        d = pickle.load(fh)
    s2 = d["spread"]; ranges = d["ranges"]
    A0 = s2.get("A0_topdown_um2") or s2["diags"][0]["topdown_um2"]
    aa = np.array([g["topdown_um2"] / A0 for g in s2["diags"]])
    mz = np.array([g["maxZ_um"] for g in s2["diags"]])
    steps = np.array(s2["steps"])
    frames = s2["frames"]
    nmem = int(max(b for (_, b) in ranges))

    ipk = int(np.argmax(aa))
    # 5 timepoints: start, mid-rise, peak, mid-relax, end
    idxs = sorted(set([0, ipk // 2, ipk, (ipk + len(aa) - 1) // 2, len(aa) - 1]))
    labels = []
    recs = []
    for ix in idxs:
        fr = np.asarray(frames[ix]).copy()
        c = fr[:nmem].mean(0); c[2] = 0.0
        recs.append(fr - c)
        phase = "START" if ix == 0 else ("PEAK" if ix == ipk else
                ("rising" if ix < ipk else "relaxing" if ix < len(aa) - 1 else "END"))
        labels.append(f"{phase}\nstep {int(steps[ix])}\nA/A0={aa[ix]:.2f} maxZ={mz[ix]:.0f}µm")
    nc = len(recs)
    H = max(np.abs(r[:, :2]).max() * UM for r in recs) * 1.05
    Zmax = max(r[:, 2].max() * UM for r in recs) * 1.08
    zlim = (-2, Zmax)

    fig = plt.figure(figsize=(3.2 * nc, 9.5))
    gs = fig.add_gridspec(3, nc, height_ratios=[1.0, 1.0, 0.95])
    for k, (fr, lab) in enumerate(zip(recs, labels)):
        axt = fig.add_subplot(gs[0, k]); _panel(axt, fr, ranges, (0, 1), H, zlim)
        axt.set_title(lab, fontsize=9)
        if k == 0:
            axt.set_ylabel("TOP-DOWN xy [µm]", fontsize=9)
        axs = fig.add_subplot(gs[1, k]); _panel(axs, fr, ranges, (0, 2), H, zlim)
        axs.set_xlabel("x [µm]", fontsize=8)
        if k == 0:
            axs.set_ylabel("SIDE xz [µm]", fontsize=9)
    axc = fig.add_subplot(gs[2, :])
    axc.plot(steps, aa, "o-", color="C3", lw=2, label="A/A0 (footprint)")
    axc.axhline(1.0, ls=":", color="0.6")
    for ix in idxs:
        axc.axvline(steps[ix], color="0.7", lw=0.8, ls="--")
    axc.plot(steps[ipk], aa[ipk], "*", color="k", ms=16)
    axc.set_xlabel("spread step"); axc.set_ylabel("A/A0", color="C3")
    axc.tick_params(axis="y", labelcolor="C3")
    ax2 = axc.twinx()
    ax2.plot(steps, mz, "s-", color="C0", ms=3, label="maxZ [µm]")
    ax2.set_ylabel("maxZ [µm]", color="C0"); ax2.tick_params(axis="y", labelcolor="C0")
    axc.set_title("RISE: wetting spreads the un-cohered cells out (A/A0↑, maxZ≈const)   →   "
                  "RELAX: cohesion pulls them back (A/A0↓ AND maxZ↓ = compacting)", fontsize=10)

    fig.suptitle(f"{args.title} — overshoot then relax\nThe 1.87 peak is TRANSIENT: the clean "
                 "ball's cells spread before cohesion engages, then cohesion retracts them "
                 "toward the compact-ball equilibrium (~A/A0 1)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = OUT / args.out
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=135)
    print(f"wrote {p}  peak A/A0={aa[ipk]:.2f}@step{int(steps[ipk])} "
          f"end={aa[-1]:.2f}  maxZ {mz[0]:.0f}->{mz[-1]:.0f}µm")


if __name__ == "__main__":
    main()
