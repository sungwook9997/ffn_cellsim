"""S4 analysis (DCM real-timescale redesign): does cadherin maturation gate the aggregate-sigma
compaction RATE? Overlay porosity(t) / Rg(t) for cluster no-maturation vs maturation tau=30s vs
tau=600s (native N=400, sigma=5mN/m, loose fcc). VERDICT (2026-07-10): NO — the tau=30s and tau=600s
trajectories nearly coincide despite 20x different tau_mature, and the first ~20s is identical across
all three (drag-limited gap-closing). The junction properties modulate the compaction EXTENT (nascent
weak junctions let sigma densify further) but NOT the timescale. => the aggregate-sigma liquid-drop
drive bypasses the junction rate-limiter; the (b) T1-neighbour-exchange escalation is needed for a
maturation-rate-limited (min-hr) compaction. See DCM_CADHERIN_CLUSTER_REARRANGEMENT_DESIGN_2026-07-09 6c.

Usage:  python dcm_cadherin_s4_analyze.py [npz_dir]   (npz named agg_compaction_n400_sig5e-03_adt0.008<TAG>.npz)
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull

NPZ_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/ff_scratch/_prod_out")
FIGDIR = "ffn_sim/outputs/h_dcm_two_stage/figs"
os.makedirs(FIGDIR, exist_ok=True)
ACCEL_DT = 8e-3
CONDS = [("_clu_nomat", "cluster, no maturation", "tab:gray"),
         ("_clu_mat30", "cluster + maturation τ=30s (engages)", "tab:blue"),
         ("_clu_mat600", "cluster + maturation τ=600s (lit)", "tab:red")]


def load(tag):
    f = f"{NPZ_DIR}/agg_compaction_n400_sig5e-03_adt0.008{tag}.npz"
    if not os.path.exists(f):
        return None
    d = np.load(f)
    return d["frames"], d["faces"].astype(int), d["cof"].astype(int), (d["step"] if "step" in d.files else None)


def metrics(frames, faces, cof):
    nc = int(cof.max()) + 1
    poros, rg = [], []
    for p in frames:
        try:
            Vh = float(ConvexHull(p).volume)
            v0, v1, v2 = p[faces[:, 0]], p[faces[:, 1]], p[faces[:, 2]]
            tv = np.einsum("ij,ij->i", v0, np.cross(v1, v2)) / 6.0
            V = np.zeros(nc); np.add.at(V, cof[faces[:, 0]], tv)
            poros.append(1.0 - np.abs(V).sum() / Vh)
        except Exception:
            poros.append(np.nan)
        cen = np.array([p[cof == c].mean(0) for c in range(nc)])
        rg.append(float(np.sqrt(((cen - cen.mean(0)) ** 2).sum(1).mean())) * 1e6)
    return np.array(poros), np.array(rg)


def main():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    print(f"{'condition':40s} {'poros0':>7} {'porosF':>7} {'Δporos':>7} {'t_half[s]':>9}")
    for tag, lab, col in CONDS:
        r = load(tag)
        if r is None:
            print(f"{lab:40s}  (npz not found)"); continue
        frames, faces, cof, step = r
        poros, rg = metrics(frames, faces, cof)
        t = (step * ACCEL_DT) if step is not None else np.arange(len(frames))
        a1.plot(t, poros, "-o", ms=3, color=col, label=lab)
        a2.plot(t, rg, "-o", ms=3, color=col, label=lab)
        p0, pf = poros[0], np.nanmin(poros)
        thalf = np.nan
        if p0 - pf > 1e-3:
            below = np.where(poros <= p0 - 0.5 * (p0 - pf))[0]
            if below.size:
                thalf = float(t[below[0]])
        print(f"{lab:40s} {poros[0]:7.3f} {np.nanmin(poros):7.3f} {p0-pf:7.3f} {thalf:9.2f}")
    a1.set_xlabel("physical time [s]"); a1.set_ylabel("porosity  (1 − ΣVcell/Vhull)")
    a1.set_title("Compaction porosity(t) — does maturation gate the rate?")
    a1.legend(fontsize=8); a1.grid(alpha=0.3)
    a2.set_xlabel("physical time [s]"); a2.set_ylabel("Rg [µm]"); a2.set_title("Aggregate Rg(t)")
    a2.legend(fontsize=8); a2.grid(alpha=0.3)
    fig.tight_layout()
    path = f"{FIGDIR}/dcm_cadherin_s4_compaction_rate.png"
    fig.savefig(path, dpi=130)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
