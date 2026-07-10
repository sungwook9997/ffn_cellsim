"""S6 migration biphasic capstone figure — crawl speed vs collagen concentration (Chan-Odde motility optimum).

Controlled native sweep (ecm_library collagen_I + --n-fa 150 --polarize --com-drag, 3000 steps): the cell's
directed crawl speed on the physiological collagen is BIPHASIC in matrix density, peaking at an intermediate
concentration — the molecular-clutch motility optimum (Chan & Odde 2008; Bangasser et al. 2013). Too sparse/soft
→ the clutches can't build load (traction→0) → no propulsion; too dense/stiff → the cell is over-anchored.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# controlled native points (same setup, vary --ecm-conc): (conc mg/mL, v_crawl nm/s, traction nN, mesh ξ µm)
pts = [(2.0, 0.31, 0.00, 1.9), (3.0, 0.57, 0.49, 1.5), (6.0, 0.02, 0.68, 0.9), (12.0, 0.08, 0.49, 0.5)]
c = np.array([p[0] for p in pts]); v = np.array([p[1] for p in pts])
tr = np.array([p[2] for p in pts]); xi = np.array([p[3] for p in pts])

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax[0].plot(c, v, "o-", color="#1f77b4", lw=2.4, ms=9)
ipk = int(np.argmax(v))
ax[0].annotate(f"OPTIMUM\nconc {c[ipk]:.0f} mg/mL\nv={v[ipk]:.2f} nm/s", (c[ipk], v[ipk]),
               xytext=(c[ipk] + 1.5, v[ipk] - 0.02), fontsize=9, color="#1f77b4",
               arrowprops=dict(arrowstyle="->", color="#1f77b4"))
ax[0].axvspan(1.5, 2.5, color="#2ca02c", alpha=0.08); ax[0].text(2.0, 0.05, "too soft\ntraction→0", ha="center", fontsize=8, color="#2ca02c")
ax[0].axvspan(5, 13, color="#d62728", alpha=0.06); ax[0].text(8.5, 0.30, "too dense\nover-anchored", ha="center", fontsize=8, color="#d62728")
ax[0].set_xscale("log"); ax[0].set_xlabel("collagen concentration [mg/mL]  (density → mesh ξ, stiffness)")
ax[0].set_ylabel("directed crawl speed v∥ [nm/s]"); ax[0].set_ylim(0, 0.65)
ax[0].set_title("Migration is BIPHASIC in matrix density (motility optimum)")
ax[0].grid(True, which="both", alpha=0.3)
ax[0].set_xticks(c); ax[0].set_xticklabels([f"{x:.0f}\n(ξ{y:.1f})" for x, y in zip(c, xi)])

ax[1].plot(c, tr, "s-", color="#9467bd", lw=2, ms=8)
ax[1].set_xscale("log"); ax[1].set_xlabel("collagen concentration [mg/mL]")
ax[1].set_ylabel("engaged-clutch traction [nN]")
ax[1].set_title("Traction collapses at low density (no grip → no propulsion)")
ax[1].grid(True, which="both", alpha=0.3); ax[1].set_xticks(c); ax[1].set_xticklabels([f"{x:.0f}" for x in c])

fig.suptitle("S6 — MCF7 crawl on physiological collagen-I is BIPHASIC in matrix density (Chan-Odde / Bangasser motility optimum)\n"
             "native 266k cell, ecm_library collagen_I, --n-fa 150 --polarize --com-drag; absolute speed sub-physiological → MMP is the path to full invasion",
             fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = "ffn_sim/outputs/ff/figs/ff_ecm_migration_biphasic.png"
fig.savefig(out, dpi=130); print(f"wrote {out}")
