"""Figure + interactive animated HTML for ff_protrusion_into_ecm (native full-cell filopodium into ECM)."""
from __future__ import annotations

import argparse
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from aleph.scripts.ff_viewer_html import build_viewer
from aleph.laws.polymerization_warp import resolve_polymerization, ratchet_velocity_np


def _groups(d):
    """Return the viz-local index slices for cortex / filo / ecm / nucleus, + seg-pair builders."""
    ncs = int(d["n_cortex_sub"]); nff = int(d["n_filo_fib"]); nbf = int(d["nb_filo"])
    nef = int(d["n_ecm_fib"]); nbe = int(d["nb_ecm"])
    filo0 = ncs; ecm0 = ncs + nff * nbf; nuc0 = ecm0 + nef * nbe
    cortex_sl = slice(0, ncs); filo_sl = slice(filo0, ecm0); ecm_sl = slice(ecm0, nuc0); nuc_sl = slice(nuc0, None)

    def seg_pairs(base, nfib, nb):
        pr = []
        for f in range(nfib):
            for k in range(nb - 1):
                pr.append([base + f * nb + k, base + f * nb + k + 1])
        return np.array(pr, np.int64)

    return dict(cortex=cortex_sl, filo=filo_sl, ecm=ecm_sl, nuc=nuc_sl,
                filo_seg=seg_pairs(filo0, nff, nbf), ecm_seg=seg_pairs(ecm0, nef, nbe))


def to_lines(frame, seg):
    """(Nviz,3) frame + (S,2) seg-pairs → (S,2,3) line-segment endpoint pairs."""
    return frame[seg].astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="aleph/outputs/ff/figs/protrusion_ecm.npz")
    ap.add_argument("--out", default="aleph/outputs/ff/figs/protrusion_ecm")
    args = ap.parse_args()
    d = np.load(args.npz, allow_pickle=True)
    g = _groups(d)
    fe, ff = d["frames_ecm"], d["frames_free"]          # (macro, Nviz, 3)
    tip_ecm, tip_free = d["tip_ecm"], d["tip_free"]
    load, vmed_ecm, vmed_free = d["load"], d["vmed_ecm"], d["vmed_free"]
    v0 = float(d["v0"]); dt = float(d["dt_phys"]); macro = fe.shape[0]; t = np.arange(macro) * dt
    poly = resolve_polymerization(G_actin_uM=20.0)

    # ---------- interactive animated HTML ----------
    f0 = fe[0]
    scenes = {"cell + filopodium → ECM (protrusion)": [
        {"name": "cortex (subsample)", "kind": "points", "color": "#3b6ea5", "size": 1.5,
         "verts": f0[g["cortex"]]},
        {"name": "nucleus", "kind": "points", "color": "#c65b7c", "size": 2.0, "verts": f0[g["nuc"]]},
        {"name": "ECM (collagen Mikado)", "kind": "lines", "color": "#7a7a7a",
         "verts": to_lines(f0, g["ecm_seg"]), "frames": [to_lines(fr, g["ecm_seg"]) for fr in fe]},
        {"name": "filopodium (polymerizing)", "kind": "lines", "color": "#ffb000",
         "verts": to_lines(f0, g["filo_seg"]), "frames": [to_lines(fr, g["filo_seg"]) for fr in fe]},
    ], "filopodium FREE (no ECM control)": [
        {"name": "cortex (subsample)", "kind": "points", "color": "#3b6ea5", "size": 1.5,
         "verts": ff[0][g["cortex"]]},
        {"name": "nucleus", "kind": "points", "color": "#c65b7c", "size": 2.0, "verts": ff[0][g["nuc"]]},
        {"name": "filopodium (free)", "kind": "lines", "color": "#20b070",
         "verts": to_lines(ff[0], g["filo_seg"]), "frames": [to_lines(fr, g["filo_seg"]) for fr in ff]},
    ]}
    build_viewer(scenes, out=f"{args.out}.html",
                 title="FF native full-cell — filopodium polymerization protrusion into a collagen ECM")
    print("wrote", f"{args.out}.html")

    # ---------- figure ----------
    fig = plt.figure(figsize=(16, 4.8))
    ax0 = fig.add_subplot(1, 3, 1); ax1 = fig.add_subplot(1, 3, 2); ax2 = fig.add_subplot(1, 3, 3, projection="3d")
    # (a) tip advance
    ax0.plot(t, tip_ecm - tip_ecm[0], "o-", color="#ffb000", ms=3, label="into ECM (emergent load)")
    ax0.plot(t, tip_free - tip_free[0], "s-", color="#20b070", ms=3, label="free (no ECM control)")
    ax0.set_xlabel("time [s]"); ax0.set_ylabel("filopodium tip advance [µm]")
    ax0.set_title("(a) protrusion: free vs ECM-impeded"); ax0.legend(fontsize=8); ax0.grid(alpha=0.3)
    # (b) emergent load + v(f) vs analytic
    axb = ax1.twinx()
    ax1.plot(t, load, "-", color="#d62728", label="emergent ECM load / tip [pN]")
    ax1.set_ylabel("emergent load [pN]", color="#d62728"); ax1.tick_params(axis="y", labelcolor="#d62728")
    axb.plot(t, vmed_ecm, "-", color="#ffb000", label="v (into ECM)")
    axb.plot(t, ratchet_velocity_np(load, poly), "k:", lw=1.6, label="analytic v(load)")
    axb.plot(t, vmed_free, "-", color="#20b070", alpha=0.6, label="v (free) = v0")
    axb.set_ylabel("polymerization v [µm/s]"); axb.set_xlabel("time [s]")
    ax1.set_xlabel("time [s]"); ax1.set_title("(b) emergent load throttles the ratchet v(f)")
    h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = axb.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=7, loc="center right")
    # (c) 3D final snapshot: cortex + nucleus + filo + ECM
    fl = fe[-1]
    ax2.scatter(fl[g["cortex"]][:, 0], fl[g["cortex"]][:, 1], fl[g["cortex"]][:, 2], s=1, c="#3b6ea5", alpha=0.25)
    ax2.scatter(fl[g["nuc"]][:, 0], fl[g["nuc"]][:, 1], fl[g["nuc"]][:, 2], s=2, c="#c65b7c", alpha=0.5)
    for seg in g["ecm_seg"]:
        p = fl[seg]; ax2.plot(p[:, 0], p[:, 1], p[:, 2], color="#7a7a7a", lw=0.4, alpha=0.5)
    for seg in g["filo_seg"]:
        p = fl[seg]; ax2.plot(p[:, 0], p[:, 1], p[:, 2], color="#ffb000", lw=1.3)
    ax2.set_title("(c) native full-cell + filopodium → ECM (final)"); ax2.set_xlabel("x [µm]")
    try:
        ax2.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    meta = json.load(open(args.npz.replace(".npz", ".json"))) if False else None
    fig.suptitle("FF native full-compartment cell (cortex+turgor+membrane+nucleus+cytoplasm) — filopodium "
                 f"polymerization protrusion into a collagen ECM Mikado network  (Δtip free {tip_free[-1]-tip_free[0]:+.2f}µm "
                 f"vs ECM {tip_ecm[-1]-tip_ecm[0]:+.2f}µm; load→{load.max():.1f}pN)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(f"{args.out}.png", dpi=130); print("wrote", f"{args.out}.png")


if __name__ == "__main__":
    main()
