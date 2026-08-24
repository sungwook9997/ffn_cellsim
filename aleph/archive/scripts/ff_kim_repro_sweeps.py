"""Aggregate the NATIVE FF ⟷ Slater/Kim 2021 parameter-sweep runs into Kim's Fig-2–5-style TREND figures —
the strong validation (reproducing Kim's parameter DEPENDENCIES, not just one run).

Reads ``outputs/ff/kim_repro/native/kim_repro_*.json`` (contraction / crosslink-density / koff sweeps + the
r^-1 ensemble) and renders four panels:
  (1) peak cortical load vs contraction strength     — Kim: peak stress ∝ cortical contraction ks,c;
  (2) relaxation % and transmission exponent n vs ⟨z⟩ — Kim: denser cross-links → slower relaxation, longer range;
  (3) stress-relaxation curves across unbinding rate koff (fixed time axis) — Kim: faster unbinding → faster relaxation;
  (4) the r^-1 ensemble decay-exponent distribution vs Kim's n=1.

Run:  python -m aleph.scripts.ff_kim_repro_sweeps
Out:  aleph/outputs/ff/kim_repro/figs/kim_sweeps.png
"""

from __future__ import annotations

import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAT = "aleph/outputs/ff/kim_repro/native"
FIGS = "aleph/outputs/ff/kim_repro/figs"


def load(pattern):
    out = []
    for p in sorted(glob.glob(f"{NAT}/{pattern}")):
        try:
            out.append(json.load(open(p)))
        except Exception:
            pass
    return out


def main():
    os.makedirs(FIGS, exist_ok=True)
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.7))

    # (1) contraction sweep — peak stress vs contraction strength
    c = load("kim_repro_contr*.json")
    if c:
        strength = np.array([1.0 - d["params"]["contract_frac"] for d in c])   # more contraction = smaller frac
        peak = np.array([d["peak_load_on_pN"] for d in c])
        o = np.argsort(strength)
        ax[0].plot(strength[o], peak[o], "o-", color="crimson", ms=9)
        ax[0].set_xlabel("contraction strength  (1 − frac)"); ax[0].set_ylabel("peak cortical load [pN]")
        ax[0].set_title("Peak stress rises with contraction\n(Kim: peak ∝ ks,c)"); ax[0].grid(alpha=0.3)

    # (2) crosslink-density sweep — relaxation & transmission vs connectivity ⟨z⟩
    z = load("kim_repro_z*.json")
    if z:
        zz = np.array([d["connectivity_z"] for d in z])
        rel = np.array([100 * d["relax_frac_on"] for d in z])
        nex = np.array([d["n_cyl_final_on"] for d in z])
        o = np.argsort(zz)
        ax[1].plot(zz[o], rel[o], "o-", color="teal", label="relaxation %")
        ax[1].set_xlabel("connectivity ⟨z⟩"); ax[1].set_ylabel("relaxation %", color="teal")
        ax[1].tick_params(axis="y", labelcolor="teal")
        a1b = ax[1].twinx()
        a1b.plot(zz[o], nex[o], "s--", color="purple", label="decay n")
        a1b.axhline(1.0, color="k", ls=":", lw=1)
        a1b.set_ylabel("transmission exponent n", color="purple"); a1b.tick_params(axis="y", labelcolor="purple")
        ax[1].set_title("Denser network: less relaxation,\nlonger-range (Kim connectivity effect)"); ax[1].grid(alpha=0.3)

    # (3) koff sweep — relaxation curves on a fixed time axis. HONEST: under strong contraction the
    #     force-accelerated (Bell) unbinding SATURATES, so relaxation is ~koff-independent here (all overlap) —
    #     Kim's koff-dependence lives at moderate forces, a regime this strong-contraction sweep doesn't probe.
    k = load("kim_repro_koff*.json")
    if k:
        for d in sorted(k, key=lambda x: x["params"]["koff_per_s"]):
            s = np.array(d["on"]["sigma_membrane_Pa"]); t = np.array(d["on"]["t_s"])
            ax[2].plot(t, s / s[0], "o-", ms=4, label=f"k_off={d['params']['koff_per_s']}/s (relax {100*d['relax_frac_on']:.0f}%)")
        ax[2].set_xlabel("time [s]"); ax[2].set_ylabel("normalized cortical load")
        ax[2].set_title("koff sweep — curves OVERLAP:\nBell unbinding SATURATED (strong contraction)"); ax[2].legend(fontsize=7); ax[2].grid(alpha=0.3)

    # (4) r^-1 ensemble — average the σ(r) PROFILES (robust) not the noisy per-seed exponents (which scatter
    #     −0.8..1.2 with unphysical negatives). HONEST: WLC+⟨z⟩3.2 transmits FAR longer-range than linear
    #     (n≈5–12), landing in the r^-1 REGIME, but the exponent is too noisy to claim a clean n=1.
    ens = load("kim_repro_ens*.json") + load("kim_native_wlc.json") + load("kim_repro_geom.json")
    if ens:
        curves = []
        for d in ens:
            pf = d["on"]["profiles"][-1]
            r = np.asarray(pf["r_cyl"]); s = np.asarray(pf["sigma_cyl"])
            ok = (r > 0) & (s > 0)
            if ok.sum() > 4:
                ax[3].loglog(r[ok], s[ok] / s[ok][0], color="crimson", alpha=0.25, lw=1)
                curves.append(np.interp(np.log(r[ok][0] * (r[ok][-1] / r[ok][0]) ** np.linspace(0, 1, 12)),
                                        np.log(r[ok]), np.log(s[ok] / s[ok][0])))
        non = np.array([d["n_cyl_final_on"] for d in ens])
        if curves:
            rg = ens[0]["on"]["profiles"][-1]["r_cyl"]; r0, r1 = min([x for x in rg if x > 0]), max(rg)
            rr = r0 * (r1 / r0) ** np.linspace(0, 1, 12)
            mean_log = np.mean(curves, axis=0)
            ax[3].loglog(rr, np.exp(mean_log), "o-", color="navy", lw=2, label="ensemble-mean σ(r)")
            ax[3].loglog(rr, (rr / rr[0]) ** -1.0, "k--", lw=1.6, label="Kim r$^{-1}$")
        ax[3].set_xlabel("r [µm]"); ax[3].set_ylabel("σ(r)/σ(r₀)")
        ax[3].set_title(f"r$^{{-1}}$ regime (WLC, N={non.size}): long-range,\nbut n noisy ({np.mean(non):.2f}±{np.std(non):.2f}, not clean n=1)")
        ax[3].legend(fontsize=8); ax[3].grid(alpha=0.3, which="both")

    fig.suptitle("FF ⟷ Slater/Kim 2021 — NATIVE parameter sweeps (honest): peak∝contraction ✓, relaxation↓ with ⟨z⟩ ✓, "
                 "ensemble-mean transmission ≈ r$^{-1}$ (single-seed noisy); koff SATURATED in this regime", fontsize=11)
    fig.tight_layout()
    out = f"{FIGS}/kim_sweeps.png"
    fig.savefig(out, dpi=130)
    print(f"[sweeps] contraction N={len(c)}, ⟨z⟩ N={len(z)}, koff N={len(k)}, ensemble N={len(ens)}", flush=True)
    print(f"[out] {out}", flush=True)


if __name__ == "__main__":
    main()
