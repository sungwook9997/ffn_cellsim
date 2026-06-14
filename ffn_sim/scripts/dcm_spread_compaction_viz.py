"""Compaction diagnostic figure for a DCM two-stage spread pickle.

Renders the START-of-spread aggregate vs the FINAL spread state (top-down + side,
per-cell shaded polygon mesh) alongside the A/A0(top-down) trajectory, so the
spheroid-compacts-not-spreads finding is visible at a glance. Reads the
``two_stage_n{N}.pkl`` written by dcm_two_stage_production. Output → figs/.

Usage:
  python -m ffn_sim.scripts.dcm_spread_compaction_viz --pkl /tmp/wetting_proxy_n12.pkl \
      --out figs/n12_wetting_proxy_compaction.png --xy-um 45 \
      --title "N=12 subdiv-2 · wetting+active-proxy spread"
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
OUT = Path("ffn_sim/outputs/h_dcm_two_stage")


def _vg(P, tris):
    a, b, c = P[tris[:, 0]], P[tris[:, 1]], P[tris[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b - a, c - a)).sum() / 6.0)


def _draw2d(ax, frame, ranges, h, axes, ttl):
    """Robust 2D projection (no 3D autoscale quirks): per-cell node cloud + convex
    hull outline. axes=(0,1) top-down xy; axes=(0,2) side xz."""
    from scipy.spatial import ConvexHull
    i, j = axes
    cmap = plt.get_cmap("tab20").colors
    allp = []
    for c, (a, b) in enumerate(ranges):
        verts = frame[a:b] * UM
        if verts.shape[0] < 3:
            continue
        allp.append(verts)
        base = cmap[c % len(cmap)]
        ax.scatter(verts[:, i], verts[:, j], s=2.0, color=base, alpha=0.55,
                   linewidths=0)
    if allp:
        P = np.concatenate(allp)[:, [i, j]]
        try:
            hull = ConvexHull(P)
            poly = P[hull.vertices]
            ax.fill(poly[:, 0], poly[:, 1], facecolor="none", edgecolor="k",
                    lw=1.4, ls="--")
        except Exception:
            pass
    ax.set_xlim(-h, h)
    if j == 1:
        ax.set_ylim(-h, h); ax.set_aspect("equal")
        ax.set_xlabel("x[µm]", fontsize=8); ax.set_ylabel("y[µm]", fontsize=8)
    else:
        ax.set_ylim(-2, max(20.0, frame[:, 2].max() * UM * 1.1)); ax.set_aspect("equal")
        ax.set_xlabel("x[µm]", fontsize=8); ax.set_ylabel("z[µm]", fontsize=8)
    ax.set_title(ttl, fontsize=9)
    ax.grid(alpha=0.25)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", required=True)
    ap.add_argument("--xy-um", type=float, default=45.0)
    ap.add_argument("--out", default="figs/spread_compaction.png")
    ap.add_argument("--title", default="DCM spread")
    args = ap.parse_args()

    with open(args.pkl, "rb") as fh:
        d = pickle.load(fh)
    tris0 = np.asarray(d["tris0"])
    ranges = d["ranges"]
    V0 = d["V0"]
    s2 = d["spread"]
    frames, diags, steps = s2["frames"], s2["diags"], s2["steps"]
    # recenter each frame on its OWN xy-centroid so the top-down panel always frames
    # the cells (the aggregate can sit off-origin at the drop centre).
    nmem = int(max(b for (_, b) in ranges))    # membrane nodes only (exclude pool)

    def _recenter(fr):
        fr = np.asarray(fr).copy()
        c = fr[:nmem].mean(0); c[2] = 0.0
        return fr - c
    f_start, f_end = _recenter(frames[0]), _recenter(frames[-1])
    A0 = s2.get("A0_topdown_um2") or diags[0]["topdown_um2"]
    aa = np.array([g["topdown_um2"] / A0 for g in diags])
    mz = np.array([g["maxZ_um"] for g in diags])
    vv = np.array([g["VV0_mean"] for g in diags])
    cen = np.array([0.0, 0.0, 0.0])
    h = args.xy_um

    fig = plt.figure(figsize=(16, 8.6))
    views = [((0, 1), "TOP-DOWN xy (silhouette = A)"), ((0, 2), "SIDE xz")]
    for col, (frame, tag) in enumerate([(f_start, "START (aggregate)"),
                                        (f_end, "FINAL (after spread)")]):
        for row, (axes, vname) in enumerate(views):
            ax = fig.add_subplot(2, 3, row * 3 + col + 1)
            aval = aa[0] if col == 0 else aa[-1]
            mval = mz[0] if col == 0 else mz[-1]
            _draw2d(ax, frame, ranges, h, axes,
                    f"{tag} · {vname}\nA/A0={aval:.2f} maxZ={mval:.1f}µm")
    # A/A0 + maxZ trajectory panel (spans the right column)
    axc = fig.add_subplot(2, 3, 3)
    axc.plot(steps, aa, "o-", color="C3", label="A/A0 (top-down silhouette)")
    axc.axhline(1.0, ls="--", color="0.5", lw=1)
    axc.set_xlabel("spread step"); axc.set_ylabel("A/A0", color="C3")
    axc.tick_params(axis="y", labelcolor="C3")
    axc.set_title("A/A0 falls below 1 → COMPACTION", fontsize=9)
    ax2 = axc.twinx()
    ax2.plot(steps, mz, "s-", color="C0", label="maxZ [µm]", ms=3)
    ax2.set_ylabel("maxZ [µm]", color="C0"); ax2.tick_params(axis="y", labelcolor="C0")
    axd = fig.add_subplot(2, 3, 6)
    axd.plot(steps, vv, "o-", color="C2")
    axd.axhline(1.0, ls="--", color="0.5", lw=1)
    axd.set_ylim(0.5, 1.5)
    axd.set_xlabel("spread step"); axd.set_ylabel("V/V0 (mean)")
    axd.set_title("V/V0 held ~1.0 (incompressible, STABLE)", fontsize=9)

    fig.suptitle(f"{args.title}\nN={len(ranges)} subdiv-2 · A/A0 {aa[0]:.2f}→{aa[-1]:.2f} "
                 f"(spheroid COMPACTS, does not spread); V/V0 bounded {vv.min():.2f}-{vv.max():.2f}",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = OUT / args.out
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=135)
    print(f"wrote {p} | A/A0 {aa[0]:.3f}->{aa[-1]:.3f}  maxZ {mz[0]:.1f}->{mz[-1]:.1f}µm  "
          f"V/V0 {vv.min():.3f}-{vv.max():.3f}")


if __name__ == "__main__":
    main()
