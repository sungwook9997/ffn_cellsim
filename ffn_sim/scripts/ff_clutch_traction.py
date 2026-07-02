"""FF piece-3 demo — molecular-clutch traction under retrograde flow (Chan-Odde 2008 load-and-fail).

Retrograde actin flow (imposed v_retro) drags the actin against N integrin clutches bound to a fixed substrate;
each clutch stretches → builds force → ruptures at the load-dependent catch-slip rate (fa_clutch_warp KMC) →
re-engages stress-free. The collective TRACTION builds and periodically collapses (load-and-fail) — de-adhesion
emergent from catch-slip, not a latch. Uses the GPU clutch kernels; the actin advance + stress-free re-binding
are host-side. Produces the traction + bound-fraction trace (data + figure).
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import warp as wp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffn_sim.ff.fa_clutch_warp import clutch_catchslip_kmc_kernel, resolve_clutch, clutch_off_rate_np


def run(M=200, v_retro=0.02, dt=0.02, steps=4000, device="cpu", seed0=1):
    """N=M clutches on one retrograde-flowing actin node; return traction(t), bound_frac(t), t."""
    cp = resolve_clutch()
    d = device
    x_act = 0.0                                                   # actin node x (flows +x at v_retro)
    anchor = np.zeros((M, 3))                                     # substrate anchors (x set at bind time)
    bound = np.ones(M, np.int32)
    pos_d = wp.array(np.array([[0.0, 0, 0]]), dtype=wp.vec3d, device=d)
    ac_d = wp.array(np.zeros(M, np.int32), dtype=wp.int32, device=d)   # all clutches → the one actin node
    anch_d = wp.array(anchor, dtype=wp.vec3d, device=d)
    bd_d = wp.array(bound, dtype=wp.int32, device=d)
    traction, bfrac = [], []
    for s in range(steps):
        x_act += v_retro * dt                                    # retrograde flow
        pos_d = wp.array(np.array([[x_act, 0, 0]]), dtype=wp.vec3d, device=d)
        # catch-slip turnover on each clutch's current load (kernel reads |anchor - actin|)
        wp.launch(clutch_catchslip_kmc_kernel, dim=M, inputs=[pos_d, ac_d, anch_d, bd_d, wp.float64(cp.k_int),
                  wp.float64(cp.rest_um), wp.float64(cp.kc0), wp.float64(cp.xc_um), wp.float64(cp.ks0),
                  wp.float64(cp.xs_um), wp.float64(cp.kT), wp.float64(cp.cap_um), wp.float64(cp.k_on),
                  wp.float64(dt), wp.int32(seed0 + s)], device=d)
        bd = bd_d.numpy(); anchor = anch_d.numpy()
        # any clutch that just RE-bound (was 0, now 1) binds stress-free at the current actin x
        just_bound = (bd == 1) & (bound == 0)
        if just_bound.any():
            anchor[just_bound, 0] = x_act - cp.rest_um           # zero stretch (L=rest) at bind
            anch_d = wp.array(anchor, dtype=wp.vec3d, device=d)
        bound = bd.copy()
        stretch = np.maximum(x_act - anchor[:, 0] - cp.rest_um, 0.0)
        traction.append(float(cp.k_int * (stretch * (bd == 1)).sum() / 1000.0))   # nN
        bfrac.append(float(bd.mean()))
    return np.array(traction), np.array(bfrac), np.arange(steps) * dt, cp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--M", type=int, default=200)
    ap.add_argument("--out", default="ffn_sim/outputs/ff/figs/clutch_traction")
    args = ap.parse_args()
    wp.init()
    tr, bf, t, cp = run(M=args.M)
    # steady-window stats (drop the first 20%)
    w = tr[len(tr) // 5:]
    print(f"clutch traction: mean {w.mean():.3f} nN  peak {tr.max():.3f}  bound frac {bf[len(bf)//5:].mean():.2f}  "
          f"F*={cp.F_star_pN:.1f}pN")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    ax[0].plot(t, tr, color="#7d4bcb", lw=0.7)
    ax[0].axhline(w.mean(), color="k", ls="--", lw=0.8, label=f"mean {w.mean():.2f} nN")
    ax[0].set_xlabel("time [s]"); ax[0].set_ylabel("traction [nN]")
    ax[0].set_title(f"(a) clutch traction — stochastic catch-slip, fluctuating ({args.M} clutches)"); ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)
    axb = ax[0].twinx(); axb.plot(t, bf, color="#c98a00", lw=0.5, alpha=0.5); axb.set_ylabel("bound fraction", color="#c98a00")

    # (b) catch-slip off-rate curve (the mechanism)
    f = np.linspace(0, 40, 300)
    ax[1].plot(f, clutch_off_rate_np(f, cp), color="#1a7f37")
    ax[1].axvline(cp.F_star_pN, ls=":", color="k"); ax[1].annotate(f"F*={cp.F_star_pN:.1f} pN\n(catch peak)",
                 (cp.F_star_pN, clutch_off_rate_np(cp.F_star_pN, cp)), (12, 0.85), fontsize=8)
    ax[1].set_xlabel("clutch load F [pN]"); ax[1].set_ylabel("off-rate [1/s]")
    ax[1].set_title("(b) integrin catch-slip off-rate (Kong 2009; ⚠️F*=7 vs KB 30 pN)"); ax[1].grid(alpha=0.3)

    fig.suptitle("FF piece-3/5 — integrin FA-clutch traction (GPU catch-slip; de-adhesion emergent from load)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(f"{args.out}.png", dpi=130); print("wrote", f"{args.out}.png")
    json.dump({"mean_traction_nN": float(w.mean()), "peak_nN": float(tr.max()),
               "bound_frac": float(bf[len(bf)//5:].mean()), "F_star_pN": cp.F_star_pN, "M": args.M},
              open(f"{args.out}.json", "w"), indent=2)


if __name__ == "__main__":
    main()
