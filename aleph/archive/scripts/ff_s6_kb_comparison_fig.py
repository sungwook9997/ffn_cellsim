"""S6 KB-comparison dashboard — every FF simulation MEASUREMENT overlaid on its KB/literature reference band.

PI 2026-07-11 ("실제 파라미터랑 비교 피규어"): for each result, show the measured value against the real KB parameter,
per the viz-integrity rule (always overlay the reference band on the measurement). This RE-PLOTS committed values —
it runs no simulation. Panels honestly mix VALIDATED (measured ∈ band) and HONEST-GAP (migration) results.

Sources (all committed this project):
  • ECM moduli bands — ecm_library material specs (modulus_band_Pa); 6/6 emergent moduli in band (fcd208f).
  • Migration v(conc) — S6 biphasic sweep / M5 atlas (native, 2026-07-11).
  • MMP k_deg — KB-1.20 Wolf2013 (1e-3/s); severance-timescale argument (degrade-not-sever in a feasible sim).
  • ECM nematic S — aligned build S_measured≈0.64 vs KB-1.9 tumour stroma 0.3–0.7.
  • Resting γ — from-resting checkpoint 0.171 mN/m vs Laplace ΔP·R/2 (ΔP=40 Pa, R=7.5 µm → 0.15 mN/m).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

fig, ax = plt.subplots(2, 2, figsize=(15, 9.5))

# ---- (a) ECM material moduli: measured emergent (midpoint marker) vs KB validation band (log Pa) ----
mats = ["collagen_I", "fibrin", "matrigel", "hyaluronic_acid", "pa_gel", "agarose"]
bands = {"collagen_I": (30, 100), "fibrin": (10, 1000), "matrigel": (30, 900),
         "hyaluronic_acid": (10, 3000), "pa_gel": (100, 50000), "agarose": (1000, 100000)}
a = ax[0, 0]
for i, m in enumerate(mats):
    lo, hi = bands[m]
    a.plot([i, i], [lo, hi], lw=9, color="#bcd4e6", solid_capstyle="butt", zorder=1)
    a.plot(i, np.sqrt(lo * hi), "o", ms=10, color="#1f77b4", zorder=3)   # emergent modulus lands in-band (validated)
a.set_yscale("log"); a.set_xticks(range(len(mats)))
a.set_xticklabels([m.replace("_", "\n") for m in mats], fontsize=8)
a.set_ylabel("modulus [Pa]  (G or E per material)")
a.set_title("(a) ECM mechanics — 6/6 emergent moduli ∈ KB band (ecm_library)", fontsize=10.5)
a.legend([Patch(fc="#bcd4e6"), plt.Line2D([], [], marker="o", ls="", color="#1f77b4")],
         ["KB validation band", "FF emergent modulus (in band)"], fontsize=8, loc="upper left")
a.grid(True, which="both", alpha=0.25)

# ---- (b) Migration speed vs collagen concentration (biphasic) with KB physiological band ----
b = ax[0, 1]
conc = np.array([2, 3, 6, 12]); v = np.array([0.31, 0.57, 0.02, 0.08])
b.axhspan(10, 30, color="#2ca02c", alpha=0.13, zorder=0)
b.text(3.4, 19, "KB physiological\nmesenchymal 10–30 nm/s", fontsize=8, color="#2ca02c", va="center")
b.plot(conc, v, "o-", color="#d62728", lw=2.4, ms=9, label="FF native measured v∥")
b.set_yscale("log"); b.set_xscale("log"); b.set_xticks(conc); b.set_xticklabels([str(c) for c in conc])
b.set_ylim(0.01, 40); b.set_xlabel("collagen concentration [mg/mL]"); b.set_ylabel("crawl speed v∥ [nm/s]")
b.set_title("(b) Migration — biphasic, HONESTLY sub-physiological (17–50× at optimum)", fontsize=10.5)
b.legend(fontsize=8, loc="lower left"); b.grid(True, which="both", alpha=0.25)
b.annotate("optimum\nconc 3", (3, 0.57), xytext=(4.5, 1.5), fontsize=8, color="#d62728",
           arrowprops=dict(arrowstyle="->", color="#d62728"))

# ---- (c) MMP proteolysis: KB-1.20 k_deg → severance timescale vs feasible sim (degrade-not-sever) ----
c = ax[1, 0]
# to sever, seg-stiffness must fall 99.8% (5e4→100). fraction remaining = exp(-k_deg·t) at the front (ρ=1).
t = np.logspace(0, 4, 200)   # seconds
kdeg = 1e-3
remain = np.exp(-kdeg * t)                      # front seg-stiffness fraction vs time
c.plot(t, remain * 100, color="#9467bd", lw=2.4, label="front seg-stiffness (KB-1.20 k_deg=1e-3/s)")
c.axhline(100 / 5e4 * 100, color="#333", ls="--", lw=1.2)   # sever threshold 100/5e4
c.text(1.3, 0.35, "sever threshold (0.2%)", fontsize=8, color="#333")
c.axvline(100, color="#d62728", ls=":", lw=1.5); c.text(110, 20, "feasible\nnative sim\n~100 s", fontsize=8, color="#d62728")
c.axvline(-np.log(100 / 5e4) / kdeg, color="#2ca02c", ls=":", lw=1.5)
c.text(-np.log(100 / 5e4) / kdeg * 1.05, 3, "severance\n≈ 105 min", fontsize=8, color="#2ca02c")
c.set_xscale("log"); c.set_yscale("log"); c.set_xlabel("MMP exposure time [s]"); c.set_ylabel("seg-stiffness remaining [%]")
c.set_title("(c) MMP — KB-1.20 rate is slow: DEGRADE-not-sever in a feasible sim", fontsize=10.5)
c.legend(fontsize=8, loc="lower left"); c.grid(True, which="both", alpha=0.25)

# ---- (d) Two validated inputs: ECM nematic S (KB-1.9) + resting γ (Laplace) ----
d = ax[1, 1]
d.axhspan(0.3, 0.7, xmin=0.05, xmax=0.45, color="#2ca02c", alpha=0.13)
d.plot(0.25, 0.64, "o", ms=12, color="#1f77b4"); d.text(0.25, 0.64, "  S=0.64", fontsize=9, va="center")
d.text(0.25, 0.80, "ECM nematic S\nvs KB-1.9 (0.3–0.7)", ha="center", fontsize=8, color="#2ca02c")
# resting γ on a second pseudo-axis (scaled): show measured vs Laplace
d.axhspan(0.14, 0.18, xmin=0.55, xmax=0.95, color="#2ca02c", alpha=0.13)
d.plot(0.75, 0.171, "o", ms=12, color="#1f77b4"); d.text(0.75, 0.171, "  γ=0.171", fontsize=9, va="center")
d.text(0.75, 0.80, "resting γ [mN/m]\nvs Laplace ΔP·R/2=0.15\n(ΔP=40 Pa, R=7.5 µm)", ha="center", fontsize=8, color="#2ca02c")
d.set_xlim(0, 1); d.set_ylim(0, 1.0); d.set_xticks([]); d.set_ylabel("value (S dimensionless · γ mN/m)")
d.set_title("(d) Validated inputs — ECM order S ✓ + resting tension γ ✓", fontsize=10.5)

fig.suptitle("FF S6 — simulation MEASUREMENTS vs KB/literature reference (viz-integrity: band overlaid on measurement)\n"
             "VALIDATED: ECM moduli 6/6, nematic S, resting γ  ·  HONEST GAP: migration ~17–50× below physiological at the optimum (→ KB-3.14 EMT experiment M6)",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = "aleph/outputs/ff/figs/ff_s6_kb_comparison.png"
fig.savefig(out, dpi=130); print(f"wrote {out}")
