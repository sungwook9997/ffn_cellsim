"""Size-dependent spreading morphology: why A/A0 falls with spheroid radius.

Reads the per-N clean-ball spread pickles (cleanball_n{N}.pkl), takes each spheroid
at its PEAK-spread frame, and renders top-down (xy) + side (xz) on a COMMON scale so
the shape change is directly visible: small spheroids flatten into wide low caps
(A/A0 high), large ones stay tall balls (only the basal cap wets → A/A0 ~ 1). A small
A/A0-vs-R panel ties the morphology to the trend. Output → figs/.
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
R_CELL_UM = 7.5


def _load(N, prefix):
    p = OUT / f"{prefix}{N}.pkl"
    if not p.exists():
        return None
    with open(p, "rb") as fh:
        d = pickle.load(fh)
    s2 = d["spread"]
    A0 = s2.get("A0_topdown_um2") or s2["diags"][0]["topdown_um2"]
    aa = np.array([g["topdown_um2"] / A0 for g in s2["diags"]])
    mz = np.array([g["maxZ_um"] for g in s2["diags"]])
    ipk = int(np.argmax(aa))
    fr = np.asarray(s2["frames"][ipk]).copy()
    ranges = d["ranges"]
    nmem = int(max(b for (_, b) in ranges))
    c = fr[:nmem].mean(0); c[2] = 0.0          # recenter on membrane xy
    fr = fr - c
    return dict(N=N, frame=fr, ranges=ranges, aa_pk=float(aa[ipk]),
                aa_fin=float(aa[-1]), mz_pk=float(mz[ipk]),
                R0=(N ** (1.0 / 3.0)) * R_CELL_UM)


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
        ax.scatter(v[:, i], v[:, j], s=1.2, color=cmap[c % 20], alpha=0.5, linewidths=0)
    if allp:
        P = np.concatenate(allp)
        try:
            hh = ConvexHull(P); poly = P[hh.vertices]
            ax.fill(poly[:, 0], poly[:, 1], facecolor="none", edgecolor="k", lw=1.3, ls="--")
        except Exception:
            pass
    ax.set_xlim(-h, h); ax.set_aspect("equal")
    ax.set_ylim(*( (-h, h) if j == 1 else zlim))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=int, nargs="+", default=[12, 30, 60, 100])
    ap.add_argument("--prefix", default="cleanball_n")
    ap.add_argument("--out", default="figs/cleanball_size_morphology.png")
    args = ap.parse_args()
    data = [_load(N, args.prefix) for N in args.ns]
    data = [d for d in data if d is not None]
    nc = len(data)
    # common scales
    H = max(np.abs(d["frame"][:, :2]).max() * UM for d in data) * 1.05
    Zmax = max(d["frame"][:, 2].max() * UM for d in data) * 1.08
    zlim = (-2, Zmax)

    fig = plt.figure(figsize=(3.4 * nc + 2, 9))
    gs = fig.add_gridspec(3, nc, height_ratios=[1.0, 1.0, 0.9])
    for k, d in enumerate(data):
        axt = fig.add_subplot(gs[0, k])
        _panel(axt, d["frame"], d["ranges"], (0, 1), H, zlim)
        axt.set_title(f"N={d['N']}  R0={d['R0']:.0f}µm\nA/A0 peak={d['aa_pk']:.2f}",
                      fontsize=10)
        if k == 0:
            axt.set_ylabel("TOP-DOWN xy\ny [µm]", fontsize=9)
        axs = fig.add_subplot(gs[1, k])
        _panel(axs, d["frame"], d["ranges"], (0, 2), H, zlim)
        axs.set_title(f"maxZ={d['mz_pk']:.0f}µm", fontsize=9)
        axs.set_xlabel("x [µm]", fontsize=8)
        if k == 0:
            axs.set_ylabel("SIDE xz\nz [µm]", fontsize=9)
    # trend panel spanning the bottom row
    axc = fig.add_subplot(gs[2, :])
    R = np.array([d["R0"] for d in data]); Apk = np.array([d["aa_pk"] for d in data])
    axc.plot(R, Apk, "o-", color="C3", ms=9, lw=2)
    for d in data:
        axc.annotate(f"N={d['N']}", (d["R0"], d["aa_pk"]),
                     textcoords="offset points", xytext=(0, 8), fontsize=8, ha="center")
    axc.axhline(1.0, ls=":", color="0.6")
    axc.set_xlabel("spheroid radius R0 = N^(1/3)·R_cell  [µm]")
    axc.set_ylabel("A/A0 peak")
    axc.set_title("A/A0 falls with R: small spheroids flatten (wide low cap), "
                  "large ones stay tall balls (only basal cap wets)", fontsize=10)
    axc.grid(alpha=0.3)

    fig.suptitle("DCM spheroid spreading vs SIZE — clean FCC ball + wetting (proxy OFF), "
                 "peak-spread morphology on a common scale\n(why A/A0 decreases: the basal "
                 "wetting footprint cannot keep up with the taller bulk as N grows)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = OUT / args.out
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=135)
    print(f"wrote {p}  (common H={H:.0f}µm Zmax={Zmax:.0f}µm)")


if __name__ == "__main__":
    main()
