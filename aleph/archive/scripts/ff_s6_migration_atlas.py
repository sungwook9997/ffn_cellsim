"""S6 migration-mode ATLAS — the session's synthesis: caveat→remodel, migration biphasic, mechanism progression.

One figure the PI can read at a glance. All values are native results committed this session (2026-07-10/11);
this only re-plots them, it runs nothing. Panels:
  (1) the S6 caveat resolution: coherent inward remodel, follows-not-remodels baseline → two-way collagen-substrate;
  (2) migration is BIPHASIC in matrix density (Chan-Odde motility optimum, peak ~conc 3);
  (3) directionality (remodel coherence) rose steadily across the mechanisms, but native translocation stayed
      sub-physiological — the honest migration limit.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))

# (1) caveat resolution — coherent inward remodel (nm), baseline vs native two-way (100s) vs 400s
labels1 = ["baseline\n(follows-\nnot-remodels)", "two-way\n(100 s)", "two-way\n(400 s)"]
dens = [15.7, 264.9, 391.4]
b = ax[0].bar(labels1, dens, color=["#bbbbbb", "#1f77b4", "#0d3b66"])
ax[0].set_ylabel("coherent inward remodel [nm]")
ax[0].set_title("(1) S6 caveat RESOLVED — collagen IS the substrate")
for r, v in zip(b, dens):
    ax[0].text(r.get_x() + r.get_width() / 2, v + 8, f"{v:.0f}", ha="center", fontsize=10)
ax[0].text(0.5, 0.92, "densification +264.9 nm = 17× baseline\n(traction-driven, native, stable)",
           transform=ax[0].transAxes, ha="center", va="top", fontsize=8, color="#0d3b66")

# (2) migration biphasic vs collagen concentration
c = np.array([2, 3, 6, 12]); v = np.array([0.31, 0.57, 0.02, 0.08])
ax[1].plot(c, v, "o-", color="#d62728", lw=2.4, ms=9)
ax[1].annotate("OPTIMUM\nconc 3, v=0.57", (3, 0.57), xytext=(4.2, 0.5), fontsize=9, color="#d62728",
               arrowprops=dict(arrowstyle="->", color="#d62728"))
ax[1].set_xscale("log"); ax[1].set_xticks(c); ax[1].set_xticklabels([str(x) for x in c])
ax[1].set_xlabel("collagen concentration [mg/mL]"); ax[1].set_ylabel("crawl speed v∥ [nm/s]")
ax[1].set_title("(2) Migration BIPHASIC in matrix density"); ax[1].set_ylim(0, 0.65); ax[1].grid(True, which="both", alpha=0.3)
ax[1].text(0.5, 0.06, "too soft→traction↓ · too dense→anchored\n(Chan-Odde / Bangasser motility optimum)",
           transform=ax[1].transAxes, ha="center", fontsize=8, color="#d62728")

# (3) directionality (coherence) rose across mechanisms, but translocation stayed sub-physiological
mech = ["two-way", "re-grip", "polarize", "n-fa", "mesench."]
coh = [0.05, 0.03, 0.15, 0.29, 0.44]
ax[2].plot(range(len(mech)), coh, "s-", color="#2ca02c", lw=2, ms=8)
ax[2].set_xticks(range(len(mech))); ax[2].set_xticklabels(mech, rotation=20, fontsize=8)
ax[2].set_ylabel("remodel coherence (directionality)"); ax[2].set_title("(3) Directionality ↑, but speed sub-physiological")
ax[2].grid(True, alpha=0.3)
ax[2].axhspan(0, 0.5, color="#2ca02c", alpha=0.04)
ax[2].text(0.5, 0.9, "coherence 0.05→0.44 as the motile\nprogram is added — but v stays 0.2-0.57 nm/s\n(<< physiological 10-30) → MMP invasion open",
           transform=ax[2].transAxes, ha="center", va="top", fontsize=8, color="#2ca02c")

fig.suptitle("S6 migration-mode ATLAS — MCF7 on collagen-I: coherent matrix REMODELLING resolved (native +264.9 nm); "
             "MIGRATION is biphasic in matrix density but sub-physiological (honest limit); mesenchymal + MMP invasion implemented",
             fontsize=10.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = "aleph/outputs/ff/figs/ff_s6_migration_atlas.png"
fig.savefig(out, dpi=130); print(f"wrote {out}")
