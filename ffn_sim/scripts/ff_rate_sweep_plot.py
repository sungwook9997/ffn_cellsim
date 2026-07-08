"""Rate-sweep figure — apparent modulus E_fit vs loading time, full physics (drainage + drained solid +
α-actinin crosslink turnover + membrane reservoir). Data from ff_hertz_validation.py --mode ratesweep
(N_fil=2000 CPU, Lp=1e-7, K_drained=300, koff=0.066/s, f_excess=0.25).

Shows the rate-dependent viscoelastic/poroelastic softening: fast (load ≪ τ_osm≈34s & 1/koff≈15s) → elastic/
undrained (16.5×); slow → drained+remodeled, plateauing at ~5×. The ~5× floor is the cortex inextensibility/
bending + tension-vs-Hertz-inversion effect (≈1100 Pa = TOP of the sharp-tip MCF7 band 200–1000 Pa; the gap to
the colloidal whole-cell 249 Pa is large-vs-small contact geometry).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/Users/sw1/ffn_cellsim/ffn_sim/outputs/ff/figs/ff_rate_sweep.png"
load_t = np.array([0.09, 0.90, 9.00, 89.96, 899.64])      # s (@3% strain)
E_fit = np.array([4109.9, 4035.2, 3393.8, 1116.2, 1273.9])  # Pa
v_lab = ["5", "0.5", "0.05", "0.005", "0.0005"]           # µm/s
TAU_OSM, INV_KOFF = 34.0, 15.0
E_COLLOID, E_SHARP = 249.0, (200.0, 1000.0)
MODEL, EDGE, FILL = "#2166ac", "#e08214", "#fde0b8"

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True})
fig, ax = plt.subplots(figsize=(9.5, 6.0))
fig.suptitle("FF cell — apparent modulus vs loading time (full physics: drainage + drained solid + crosslink turnover)\n"
             "2026-07-08 · rate-dependent softening 16.5×→~5×; the ~5× floor is cortex inextensibility/bending (next lever)",
             fontsize=10.5, y=0.99)

# MCF7 reference bands
ax.axhspan(E_SHARP[0], E_SHARP[1], color=FILL, alpha=0.7, zorder=0, label="MCF7 sharp-tip 200–1000 Pa (Li 2008)")
ax.axhline(E_COLLOID, color=EDGE, lw=2, ls="--", zorder=1, label="MCF7 colloidal whole-cell 249 Pa (Zbiral)")

ax.plot(load_t, E_fit, "-o", color=MODEL, lw=2.2, ms=8, zorder=3, label="model E_fit (full physics)")
for x, y, v in zip(load_t, E_fit, v_lab):
    ax.annotate(f"{v} µm/s\n{y:.0f} Pa", (x, y), (0, 12), textcoords="offset points",
                ha="center", fontsize=7.5, color="#204060")
# transition timescales
ax.axvline(INV_KOFF, color="#7a0000", lw=1.2, ls=":", zorder=2)
ax.text(INV_KOFF, 5200, "1/k_off≈15s\n(turnover)", fontsize=7.5, color="#7a0000", ha="center", va="top")
ax.axvline(TAU_OSM, color="#4d004d", lw=1.2, ls=":", zorder=2)
ax.text(TAU_OSM, 5200, "τ_osm≈34s\n(drainage)", fontsize=7.5, color="#4d004d", ha="center", va="top")
ax.text(0.12, 4400, "fast (Zbiral 5µm/s):\nelastic + undrained\n→ 16.5×", fontsize=8, color="#204060", va="top")
ax.text(300, 1500, "slow: drained + remodeled\n→ ~5× floor", fontsize=8, color="#204060", ha="center")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("loading time @3% strain  [s]   (log; = 2R·strain / v_load)")
ax.set_ylabel("apparent Young's modulus E_fit  [Pa]   (log)")
ax.set_ylim(150, 6000)
ax.legend(fontsize=8, loc="lower left")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT, dpi=120, bbox_inches="tight")
print("wrote", OUT)
