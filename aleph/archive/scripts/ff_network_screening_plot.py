"""(a) network-screening validation figure — FF contractile transmission vs reconstituted ground truth.
Left: Ronceray 2016 Table I amplification ratio σ_exp/σ_lin by regime — the DENSE actomyosin (cortex-relevant)
regime is ~unity; the 17-70× amplification is a sparse-active-unit + extensible-fiber effect that saturates. FF
runs at the dense end → no-amplification regime. Right: FF buckling control (σ/σ_dipole vs buckling seed) — robust
~0.42-0.45, buckles when seeded (not the rigid-rod ~0 over-screening artifact). Verdict: FF screening is PHYSICAL.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/Users/sw1/ffn_cellsim/aleph/outputs/ff/figs/ff_network_screening.png"
# Ronceray 2016 Table I
regimes = ["dense 3D\nactomyosin\n(cortex)", "FF network\n(dense, this)", "sparse 2D\non vesicle", "fibrin/platelet\nclot"]
ratios = [1.17, 1.0, 70.0, 17.0]   # FF ~unity vs linear (0.42 vs coherent ≈ 1 vs linear)
cols = ["#2166ac", "#5aae61", "#b2182b", "#b2182b"]
# FF buckling control
seed = np.array([0.00, 0.08, 0.20]); ratio_c = np.array([0.453, 0.417, 0.406]); zbow = np.array([0.0, 0.168, 0.237])

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True})
fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
fig.suptitle("FF contractile-network screening validated PHYSICAL vs reconstituted ground truth · 2026-07-08 · task (a)\n"
             "dense actomyosin transmits at ~unity (no amplification); FF is in that regime; the 25–50× tension gap is DENSITY, not transmission",
             fontsize=10, y=0.99)

# A — Ronceray Table I regimes
a = ax[0]
bars = a.bar(range(len(regimes)), ratios, color=cols, width=0.62)
a.axhline(1.0, color="#333", lw=1.5, ls="--", zorder=1)
a.text(3.4, 1.05, "unity (transmit = linear sum)", fontsize=7.5, color="#333", ha="right", va="bottom")
for i, r in enumerate(ratios):
    a.text(i, r*1.08, f"{r:.2f}×" if r < 2 else f"{r:.0f}×", ha="center", fontsize=9, color=cols[i])
a.set_yscale("log"); a.set_ylim(0.5, 150)
a.set_xticks(range(len(regimes))); a.set_xticklabels(regimes, fontsize=8)
a.set_ylabel("amplification ratio  σ_exp / σ_lin   (log)")
a.set_title("Ronceray 2016 Table I — amplification by regime")
a.text(0.5, 0.93, "amplification = SPARSE units + extensible fiber (saturates → ~1 as density rises)",
       transform=a.transAxes, fontsize=7.5, color="#7a0000", ha="center", va="top")

# B — FF buckling control
a = ax[1]
a.axhspan(0.35, 0.5, color="#fde0b8", alpha=0.5, zorder=0, label="FF ~unity vs linear (0.42 vs coherent)")
a.plot(seed, ratio_c, "-o", color="#5aae61", lw=2.2, ms=9, label="σ/σ_dipole (vs coherent)")
a.plot(seed, zbow, "-s", color="#2166ac", lw=2, ms=8, label="buckling amplitude z_bow [µm]")
a.annotate("planar (no buckle seed)\n→ still 0.45, NOT ~0\n(not rigid-rod artifact)", (0.0, 0.453), (0.05, 0.62),
           fontsize=7.5, color="#123", arrowprops=dict(arrowstyle="->", color="#123"))
a.set_ylim(-0.02, 0.75); a.set_xlabel("out-of-plane buckling seed  seed_z")
a.set_ylabel("σ/σ_dipole   ·   z_bow [µm]")
a.set_title("FF buckling control (f_act=20 pN)"); a.legend(fontsize=7.5, loc="upper right")

fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(OUT, dpi=120, bbox_inches="tight")
print("wrote", OUT)
