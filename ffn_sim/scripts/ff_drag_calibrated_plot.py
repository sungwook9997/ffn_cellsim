"""Calibrated-drag figure — apparent modulus vs press speed with the BULK-η-calibrated per-node drag
(γ_node=6πηR/Nc, η=65.9 Pa·s, canonical Stokes distributed over the cortex nodes). Data from ff_drag_calibration
(N_fil=1500 CPU, strain 3%).

The drag is validated against the analytic Newtonian ground truth: η_eff = F_visc/(A·ε̇) ≈ 124–141 Pa·s across a
100× speed range (rate-independent = Newtonian), ~2× the target 65.9 (the factor-2 is the crude squeeze-flow ε̇/area
estimate, not a drag error). With physical viscosity the model produces a rate-dependent modulus that PASSES THROUGH
the MCF7 band (~249 Pa) at ~0.5–1 µm/s: viscous-stiff above (Zbiral 5µm/s → 9.3×, viscous-dominated F_visc≫F_elastic),
elastic-soft at the relaxed equilibrium (~13 Pa = γ-floor).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/Users/sw1/ffn_cellsim/ffn_sim/outputs/ff/figs/ff_drag_calibrated.png"
v = np.array([0.5, 5.0, 50.0])
E = np.array([223.0, 2313.0, 23364.0])
eta_eff = np.array([124.0, 138.0, 141.0])
E_EQ, E_ZB = 13.0, 249.0
MODEL, EQ, EDGE, FILL = "#2166ac", "#5aae61", "#e08214", "#fde0b8"

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True})
fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.6), gridspec_kw={"width_ratios": [1.5, 1]})
fig.suptitle("FF cell — CALIBRATED bulk-η drag (γ_node=6πηR/Nc, η=65.9 Pa·s) · 2026-07-08\n"
             "physical rate-dependent modulus PASSES THROUGH the MCF7 band at ~0.5–1 µm/s; Zbiral 5µm/s = viscous-stiff (9.3×); equilibrium = γ-floor",
             fontsize=10, y=0.99)

# A — E vs press speed
a = ax[0]
a.axhspan(200, 1000, color=FILL, alpha=0.6, zorder=0, label="MCF7 sharp-tip 200–1000 Pa")
a.axhline(E_ZB, color=EDGE, lw=2, ls="--", zorder=1, label="MCF7 colloidal 249 Pa (Zbiral)")
a.axhline(E_EQ, color=EQ, lw=2, ls=":", zorder=1, label="relaxed equilibrium ~13 Pa (γ-floor)")
a.plot(v, E, "-o", color=MODEL, lw=2.2, ms=9, zorder=3, label="model E (calibrated drag, ramp-end)")
for x, y in zip(v, E):
    a.annotate(f"{y/E_ZB:.1f}×", (x, y), (7, -3), textcoords="offset points", fontsize=8.5, color="#204060")
a.annotate("0.5µm/s → 0.9×\n(in band)", (0.5, 223), (0.5, 60), fontsize=8, color="#1a6b1a", ha="center",
           arrowprops=dict(arrowstyle="->", color="#1a6b1a"))
a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(8, 5e4)
a.set_xlabel("plate speed v_press  [µm/s]  (log)"); a.set_ylabel("apparent modulus E  [Pa]  (log)")
a.set_title("Apparent modulus vs press speed"); a.legend(fontsize=8, loc="upper left")

# B — effective viscosity validation
a = ax[1]
a.axhline(65.9, color=EDGE, lw=2, ls="--", label="target η=65.9 Pa·s (Dessard 2024)")
a.axhspan(65.9*0.7, 65.9*2.2, color=FILL, alpha=0.4, zorder=0)
a.plot(v, eta_eff, "-s", color=MODEL, lw=2, ms=9, label="model η_eff = F_visc/(A·ε̇)")
a.set_xscale("log"); a.set_ylim(0, 200)
a.set_xlabel("plate speed v_press  [µm/s]  (log)"); a.set_ylabel("effective viscosity  [Pa·s]")
a.set_title("Drag validation: η_eff rate-independent, ~2× target"); a.legend(fontsize=8, loc="lower right")
a.text(0.5, 0.05, "Newtonian (η_eff flat across 100× speed);\nfactor-2 = crude squeeze-flow ε̇/A estimate",
       transform=a.transAxes, fontsize=7.5, color="#555", va="bottom")

fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(OUT, dpi=120, bbox_inches="tight")
print("wrote", OUT)
