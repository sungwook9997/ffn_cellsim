"""Native-scale confirmation figure — apparent undrained modulus vs relaxation steps, CPU (Nc≈2000 fil) vs
NATIVE (Nc=266,000, A5000). Both drop monotonically through many orders as the cell is allowed to finish
deforming, cross the MCF7 band, and reach the SOFT side (γ-floor). Confirms the small-strain 'stiffness' was an
under-relaxation (fast-press) transient at full resolution too; native equilibrates ~10× slower (denser mesh).
Data: ff_converge_probe.py (CPU) + ff_native_confirm.py long (A5000), strain 3%.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/Users/sw1/ffn_cellsim/ffn_sim/outputs/ff/figs/ff_native_convergence.png"
ns_cpu = np.array([800, 1600, 3200, 6400, 12800]);           E_cpu = np.array([16643, 4127, 262, 20, 19])
ns_nat = np.array([500, 4000, 8000, 16000, 32000, 64000]);   E_nat = np.array([5455562, 216163, 18080, 6588, 1417.5, 181.7])
E_ZB, BAND, FILL, EDGE, CPU, NAT = 249.0, (200.0, 1000.0), "#fde0b8", "#e08214", "#9970ab", "#2166ac"

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True})
fig, ax = plt.subplots(figsize=(9.5, 6.2))
fig.suptitle("FF cell — undrained apparent modulus vs RELAXATION steps: CPU vs NATIVE (A5000, Nc=266k) · 2026-07-08\n"
             "both drop through many orders as the cell finishes deforming → soft equilibrium (γ-floor); the 'stiffness' was an under-relaxation transient",
             fontsize=10, y=0.99)
ax.axhspan(BAND[0], BAND[1], color=FILL, alpha=0.7, zorder=0, label="MCF7 sharp-tip 200–1000 Pa")
ax.axhline(E_ZB, color=EDGE, lw=2, ls="--", zorder=1, label="MCF7 colloidal 249 Pa (Zbiral)")
ax.plot(ns_cpu, E_cpu, "-o", color=CPU, lw=2, ms=8, label="CPU (Nc≈2000 fil)")
ax.plot(ns_nat, E_nat, "-s", color=NAT, lw=2.4, ms=9, label="NATIVE (Nc=266,000, A5000)")
for x, y in zip(ns_nat, E_nat):
    ax.annotate(f"{y/E_ZB:.2g}×", (x, y), (6, 4), textcoords="offset points", fontsize=7.5, color="#123")
ax.annotate("native equilibrium heading\nto the soft γ-floor", (64000, 181.7), (20000, 40), fontsize=8.5,
            color="#123", arrowprops=dict(arrowstyle="->", color="#123"))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("relaxation steps  (log)  — proxy for how completely the cell has deformed")
ax.set_ylabel("undrained apparent modulus E  [Pa]  (log)")
ax.set_ylim(10, 1e7)
ax.legend(fontsize=8.5, loc="upper right")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUT, dpi=120, bbox_inches="tight")
print("wrote", OUT)
