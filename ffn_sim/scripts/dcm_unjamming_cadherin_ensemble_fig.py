"""Figure: cadherin × motility ensemble — unjamming (shape index) DECOUPLES from spread (A/A0).

Seed ensemble (fa100 active motility, N=100, 3 seeds each), cadherin nascent(weak) vs matured(strong):
- shape index s (jamming order parameter) is the SAME (both fluidised to s≈5.36, ~40% cells past s0*=5.41) →
  motility fluidises the tissue regardless of cadherin;
- A/A0 (aggregate spread) DIFFERS: nascent 2.44 vs matured 1.81 → cadherin sets COHESION (how far the
  fluidised cells disperse), NOT the unjamming.
Epithelial (strong E-cad) = fluid-but-cohesive (unjams, stays together — Park2015 airway epithelium);
mesenchymal (weak cad) = fluid-and-dispersing. Matches KB-PIV-10 low-E-cadherin→gas/EMT. Seed-robust + visual.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "ffn_sim/outputs/h_dcm_two_stage/figs"
os.makedirs(OUT, exist_ok=True)

# measured ensemble (3 seeds each): A/A0 and per-cell 3D shape index
aa0 = {"nascent (weak cad)": [2.54, 2.47, 2.29], "matured (strong cad)": [2.07, 1.86, 1.47]}
s = {"nascent (weak cad)": [5.36, 5.36, 5.36], "matured (strong cad)": [5.35, 5.35, 5.35]}  # mean 5.361/5.352
S0 = 5.41


def main():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5))
    conds = list(aa0.keys())
    cols = ["tab:orange", "tab:blue"]
    x = np.arange(len(conds))

    m = [np.mean(aa0[c]) for c in conds]; e = [np.std(aa0[c]) for c in conds]
    a1.bar(x, m, yerr=e, capsize=6, color=cols, alpha=0.85)
    for xi, c in zip(x, conds):
        a1.scatter([xi] * 3, aa0[c], color="k", zorder=5, s=18)
    a1.axhline(1.0, color="gray", ls=":")
    a1.set_xticks(x); a1.set_xticklabels(conds)
    a1.set_ylabel("A/A0  (aggregate spread, top-down)")
    a1.set_title("SPREAD differs with cadherin\n(nascent 2.44 vs matured 1.81) = COHESION")
    a1.grid(axis="y", alpha=0.25)

    ms = [np.mean(s[c]) for c in conds]; es = [0.051, 0.033]
    a2.bar(x, ms, yerr=es, capsize=6, color=cols, alpha=0.85)
    a2.axhline(S0, color="tab:red", ls="--", lw=1.5)
    a2.text(-0.35, S0 + 0.003, "s0*≈5.41 (unjamming)", color="tab:red", fontsize=9)
    a2.axhline(4.836, color="gray", ls=":")
    a2.text(-0.35, 4.836 - 0.02, "sphere 4.84", color="gray", fontsize=8)
    a2.set_xticks(x); a2.set_xticklabels(conds)
    a2.set_ylabel("per-cell 3D shape index s")
    a2.set_title("UNJAMMING is the SAME (s≈5.36 both)\n→ motility fluidises, cadherin doesn't")
    a2.set_ylim(4.80, 5.46)
    a2.grid(axis="y", alpha=0.25)

    fig.suptitle("cadherin × motility: unjamming (shape index) DECOUPLES from spread (A/A0)  —  seed ensemble, N=100",
                 fontsize=11)
    fig.tight_layout()
    path = f"{OUT}/dcm_unjamming_cadherin_ensemble.png"
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
