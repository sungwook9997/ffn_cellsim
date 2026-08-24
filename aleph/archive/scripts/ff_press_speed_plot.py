"""Press-speed figure — apparent modulus vs plate speed with the PHYSICAL cytoplasm drag (η=65.9 Pa·s, NF2007
mobility). Data from press_speed experiment (N_fil=2000 CPU, strain 3%, undrained/mechanical-only).

Confirms the PI's hypothesis (2026-07-08): pressing faster than the η-limited cell-deformation speed reads a
VISCOUS TRANSIENT (E ∝ v_press), NOT the elastic modulus. The relaxed equilibrium (after dwell) is ~17 Pa =
γ-floor (too soft) — so the previously-reported ~16.5× "stiffness" was a fast-press / under-relaxation transient,
not a real elastic gap. CAVEAT: the per-node drag (single-fiber NF2007 mobility) over-estimates the BULK cytoplasm
viscous response by ~100× (continuum η·ε̇≈44 Pa vs the model's huge transient) → the absolute crossover speed needs
drag calibration to the bulk η; the DIRECTION and the E∝v_press scaling are robust.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/Users/sw1/ffn_cellsim/aleph/outputs/ff/figs/ff_press_speed.png"
v = np.array([0.05, 0.2, 1.0, 5.0])                       # µm/s
E = np.array([2553.0, 10157.0, 50777.0, 251255.0])        # Pa (ramp-end, physical drag)
E_EQ = 17.0                                                # Pa (relaxed equilibrium = γ-floor)
E_ZB = 249.0
MODEL, EQ, EDGE, FILL = "#2166ac", "#5aae61", "#e08214", "#fde0b8"

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True})
fig, ax = plt.subplots(figsize=(9.5, 6.2))
fig.suptitle("FF cell — apparent modulus vs PRESS SPEED (physical η=65.9 Pa·s cytoplasm drag) · 2026-07-08\n"
             "PI confirmed: fast press reads a VISCOUS TRANSIENT (E ∝ v_press), not the elastic modulus; relaxed equilibrium = γ-floor (~17 Pa)",
             fontsize=10, y=0.99)

# reference lines
ax.axhspan(200, 1000, color=FILL, alpha=0.6, zorder=0, label="MCF7 sharp-tip 200–1000 Pa")
ax.axhline(E_ZB, color=EDGE, lw=2, ls="--", zorder=1, label="MCF7 colloidal 249 Pa (Zbiral, 5µm/s)")
ax.axhline(E_EQ, color=EQ, lw=2, ls=":", zorder=1, label="model RELAXED equilibrium ~17 Pa (γ-floor)")

# model viscous-transient points + ~linear (E∝v) guide
ax.plot(v, E, "-o", color=MODEL, lw=2.2, ms=9, zorder=3, label="model E (ramp-end, physical drag)")
vg = np.array([0.03, 6.0]); ax.plot(vg, E[-1]*(vg/5.0), "-", color=MODEL, lw=1, alpha=0.4, zorder=2)
ax.text(0.06, E[-1]*(0.06/5.0)*1.3, "slope 1\n(E ∝ v_press = viscous)", fontsize=7.5, color=MODEL)
for x, y in zip(v, E):
    ax.annotate(f"{y/E_ZB:.0f}×", (x, y), (7, -2), textcoords="offset points", fontsize=8, color="#204060")
ax.annotate("Zbiral 5µm/s →\nmodel 1009× (transient)", (5.0, 251255), (1.2, 251255), fontsize=8, color="#7a0000",
            arrowprops=dict(arrowstyle="->", color="#7a0000"))
ax.text(0.032, 40, "prior 1600-step result (16.5×) was\na point on THIS transient curve, not equilibrium",
        fontsize=7.5, color="#555", va="top")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("plate speed v_press  [µm/s]  (log)")
ax.set_ylabel("apparent Young's modulus E  [Pa]  (log)")
ax.set_ylim(8, 5e5)
ax.legend(fontsize=8, loc="upper left")
fig.text(0.5, 0.008, "CAVEAT: per-node drag (single-fiber NF2007 mobility) over-estimates bulk cytoplasm viscosity ~100× "
         "→ absolute crossover speed needs drag calibration to η; direction + E∝v_press are robust.",
         ha="center", fontsize=7.3, color="#7a0000")
fig.tight_layout(rect=[0, 0.03, 1, 0.93])
fig.savefig(OUT, dpi=120, bbox_inches="tight")
print("wrote", OUT)
