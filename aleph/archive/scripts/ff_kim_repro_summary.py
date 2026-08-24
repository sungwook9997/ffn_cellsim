"""Consolidate the NATIVE (A5000, 105k-node) Slater/Kim 2021 reproduction runs into one summary figure.

Reads the native run logs + the WLC reference JSON under ``outputs/ff/kim_repro/native/`` and renders:
  (1) stress-decay exponent n by matrix regime vs Kim's quasi-2D r^-1 (n=1) — the headline: linear-spring
      fibers transmit too locally (n≈5-12) at every connectivity, but PHYSIOLOGICAL sub-isostatic ⟨z⟩≈3.2 +
      thermal-WLC strain-stiffening brings n into the Kim r^-1 range (ensemble ON≈0.8), reproducing the
      long-range tensed-fiber force transmission;
  (2) the WLC σ(r) cylindrical-flux profile vs r^-1;
  (3) cortical-load stress relaxation (WLC) — ON relaxes, OFF (elastic) holds.

Run:  python -m aleph.scripts.ff_kim_repro_summary
Out:  aleph/outputs/ff/kim_repro/figs/kim_native_summary.png
"""

from __future__ import annotations

import glob
import json
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAT = "aleph/outputs/ff/kim_repro/native"
FIGS = "aleph/outputs/ff/kim_repro/figs"


def parse_log(path):
    """Pull ⟨z⟩, WLC flag, ON/OFF decay exponent and ON relaxation % from one native run log."""
    txt = open(path).read()
    z = re.search(r"⟨z⟩=([\d.]+)", txt)
    dec = re.search(r"exponent n: ON=([\-\d.]+) OFF=([\-\d.]+)", txt)
    rel = re.search(r"relaxes the cortical load ([\d.]+)%", txt)
    if not (z and dec):
        return None
    return {"z": float(z.group(1)), "wlc": "wlc" in os.path.basename(path),
            "n_on": float(dec.group(1)), "n_off": float(dec.group(2)),
            "relax": float(rel.group(1)) if rel else float("nan"),
            "name": os.path.basename(path)}


def main():
    os.makedirs(FIGS, exist_ok=True)
    runs = [r for r in (parse_log(p) for p in sorted(glob.glob(f"{NAT}/*.log"))) if r]
    # group into the three regimes
    spring_hi = [r for r in runs if not r["wlc"] and r["z"] > 10]        # over-connected linear
    spring_lo = [r for r in runs if not r["wlc"] and r["z"] <= 10]       # physiological linear
    wlc_lo = [r for r in runs if r["wlc"]]                               # physiological + WLC (ensemble)
    groups = [("spring\n⟨z⟩≈31", spring_hi, "steelblue"),
              ("spring\n⟨z⟩≈3.2", spring_lo, "darkorange"),
              ("WLC\n⟨z⟩≈3.2", wlc_lo, "crimson")]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 4.8))

    # (1) decay exponent by regime
    for i, (lbl, g, c) in enumerate(groups):
        non = [r["n_on"] for r in g]; noff = [r["n_off"] for r in g]
        ax1.scatter([i - 0.08] * len(non), non, c=c, marker="o", s=70, label=f"{lbl} ON" if i == 0 else None, zorder=3)
        ax1.scatter([i + 0.08] * len(noff), noff, c=c, marker="s", s=70, alpha=0.5, zorder=3)
        if non:
            ax1.annotate(f"n̄={np.mean(non):.2f}", (i - 0.08, np.mean(non)), fontsize=9,
                         ha="right", va="bottom", color=c, weight="bold")
    ax1.axhline(1.0, color="k", ls="--", lw=1.4, label="Kim quasi-2D r$^{-1}$ (n=1)")
    ax1.set_xticks(range(3)); ax1.set_xticklabels([g[0] for g in groups])
    ax1.set_ylabel("stress-decay exponent n (|σ|~r$^{-n}$)")
    ax1.set_title("Long-range transmission: only physiological\n⟨z⟩ + WLC reaches Kim's r$^{-1}$")
    ax1.legend(fontsize=8, loc="upper center"); ax1.grid(alpha=0.3)
    ax1.text(0.02, 0.02, "○ ON (viscoelastic)   □ OFF (elastic)", transform=ax1.transAxes, fontsize=8)

    # (2) WLC σ(r) profile vs r^-1
    wlcj = f"{NAT}/kim_native_wlc.json"
    if os.path.exists(wlcj):
        d = json.load(open(wlcj))
        pf = d["on"]["profiles"][-1]
        r = np.asarray(pf["r_cyl"]); s = np.asarray(pf["sigma_cyl"])
        ok = (r > 0) & (s > 0)
        ax2.loglog(r[ok], s[ok], "o-", color="crimson", label=f"FF WLC σ(r) (n={pf['n_cyl']:.2f})")
        rr = r[ok]; ax2.loglog(rr, s[ok][0] * (rr / rr[0]) ** -1.0, "k--", label="Kim r$^{-1}$")
        ax2.set_xlabel("r [µm]"); ax2.set_ylabel("σ(r) cylindrical flux [Pa]")
        ax2.set_title("WLC stress transmission vs Kim r$^{-1}$\n(single seed — noisy; ensemble n̄≈0.8)")
        ax2.legend(fontsize=8); ax2.grid(alpha=0.3, which="both")

        # (3) relaxation (WLC)
        on, off = d["on"], d["off"]
        ax3.plot(on["t_s"], on["sigma_membrane_Pa"], "o-", color="crimson", label="turnover ON (viscoelastic)")
        ax3.plot(off["t_s"], off["sigma_membrane_Pa"], "s--", color="steelblue", label="turnover OFF (elastic)")
        ax3.set_xlabel("time [s]"); ax3.set_ylabel("cortical load = mean FA tension [pN]")
        ax3.set_title(f"Stress relaxation (WLC, native)\nON relaxes {100*d['relax_frac_on']:.0f}%, OFF holds")
        ax3.legend(fontsize=8); ax3.grid(alpha=0.3)

    fig.suptitle("FF ⟷ Slater/Kim 2021 — NATIVE (A5000, 105k nodes, ⟨z⟩ physiological): reproduces BOTH the "
                 "viscoelastic relaxation AND the r$^{-1}$ transmission (with WLC strain-stiffening)", fontsize=11)
    fig.tight_layout()
    out = f"{FIGS}/kim_native_summary.png"
    fig.savefig(out, dpi=130)
    print(f"[summary] regimes: spring⟨z⟩31 n_on={[r['n_on'] for r in spring_hi]}, "
          f"spring⟨z⟩3.2 n_on={[r['n_on'] for r in spring_lo]}, "
          f"WLC⟨z⟩3.2 n_on={[round(r['n_on'],2) for r in wlc_lo]} (mean {np.mean([r['n_on'] for r in wlc_lo]):.2f})", flush=True)
    print(f"[out] {out}", flush=True)


if __name__ == "__main__":
    main()
