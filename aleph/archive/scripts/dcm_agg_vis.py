"""Render a DCM AGGREGATE spheroid as polygon SURFACES (PI inspection of stage 1).

Loads an agg-only pickle (agg frames + ranges + tris0) and renders the final
aggregated spheroid as shaded triangle-mesh surfaces — a 3D oblique view + a
top-down (xy) view + the contact/V0 trajectory — so the PI can see whether the
cells properly aggregated into a cohesive ball.

Usage:  python -m aleph.scripts.dcm_agg_vis AGG_ONLY.pkl OUT.png
"""
from __future__ import annotations

import pickle
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, LightSource
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.collections import PolyCollection

UM = 1e6
_TAB = plt.get_cmap("tab20").colors


def _lambert(verts, tris, base_rgb, ls):
    tn = np.cross(verts[tris[:, 1]] - verts[tris[:, 0]],
                  verts[tris[:, 2]] - verts[tris[:, 0]])
    nrm = np.linalg.norm(tn, axis=1, keepdims=True)
    tn = tn / np.where(nrm > 0, nrm, 1)
    shade = 0.45 + 0.55 * np.clip(tn @ np.array([0.3, 0.3, 0.9]), 0, 1)
    return np.array([(base_rgb[0] * s, base_rgb[1] * s, base_rgb[2] * s, 0.95) for s in shade])


def main():
    pkl, out = sys.argv[1], sys.argv[2]
    d = pickle.load(open(pkl, "rb"))
    agg = d["agg"]; ranges = d["ranges"]; tris0 = np.asarray(d["tris0"])
    pos = agg["frames"][-1] * UM
    ls = LightSource(azdeg=315, altdeg=45)
    g = agg["diags"][-1]

    fig = plt.figure(figsize=(16, 6))
    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)

    for c, (a, b) in enumerate(ranges):
        v = pos[a:b]
        rgba = _lambert(v, tris0, to_rgb(_TAB[c % len(_TAB)]), ls)
        ax1.add_collection3d(Poly3DCollection(v[tris0], facecolors=rgba,
                                              edgecolors="k", linewidths=0.05))
        # top-down xy polygon
        ax2.add_collection(PolyCollection(v[tris0][:, :, :2],
                           facecolors=to_rgb(_TAB[c % len(_TAB)]),
                           edgecolors="k", linewidths=0.05, alpha=0.5))
    lim = np.abs(pos).max()
    ax1.set_xlim(-lim, lim); ax1.set_ylim(-lim, lim); ax1.set_zlim(0, 2 * lim)
    ax1.set_title(f"3D surfaces (N={len(ranges)})", fontsize=10)
    ax1.set_xlabel("x µm"); ax1.set_ylabel("y µm")
    ax2.set_aspect("equal"); ax2.set_xlim(-lim, lim); ax2.set_ylim(-lim, lim)
    ax2.set_title("top-down (xy)", fontsize=10); ax2.set_xlabel("x µm"); ax2.set_ylabel("y µm")

    st = agg["steps"]
    ax3.plot(st, [x["contact_frac"] for x in agg["diags"]], "C0-o", label="contact frac")
    ax3.plot(st, [x["VV0_mean"] for x in agg["diags"]], "C1-o", label="V/V0")
    ax3.plot(st, [x["asphericity"] for x in agg["diags"]], "C3-o", label="asphericity")
    ax3.axhline(1.0, color="grey", lw=0.4, ls=":")
    ax3.set_xlabel("agg step"); ax3.set_title("aggregation trajectory", fontsize=10)
    ax3.legend(fontsize=8); ax3.grid(alpha=0.3)

    fig.suptitle(f"DCM AGGREGATE N={len(ranges)} — contact={g['contact_frac']:.2f} "
                 f"asph={g['asphericity']:.3f} V/V0={g['VV0_mean']:.2f} "
                 f"converged={agg['converged']}", fontsize=12)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
