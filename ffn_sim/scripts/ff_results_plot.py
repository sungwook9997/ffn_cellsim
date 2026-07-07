"""FF full-compartment quantitative results — record + plot (with literature reality-bands overlaid).
Data from this session's native AFM full-compartment sweep + coupling + field percentiles. Reality bands are
LITERATURE reference (not fits): MCF7 suspended cortical tension ~10 mN/m (Moazzeni 2021); interphase cortical
tension ~0.3-1 mN/m; resting turgor ΔP ~40 Pa (Fischer-Friedrich 2014, HeLa) up to ~few-hundred Pa; whole-cell
parallel-plate compression force ~10-150 nN (Fischer-Friedrich). Log axes noted on the axis label.
"""
import csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_PNG = "/Users/sw1/ffn_cellsim/ffn_sim/outputs/ff/figs/ff_fullcompartment_results.png"
OUT_CSV = "/Users/sw1/ffn_cellsim/ffn_sim/outputs/ff/data/ff_afm_results_2026-07-07.csv"

s   = np.array([0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27])
F   = np.array([0.01,13.93,102.98,331.97,800.52,1580.72,2805.57,4560.76,6974.72,10185.07])  # nN
dP  = np.array([36.4,1283.7,4841.8,10724.7,18966.7,29645.7,42654.7,58294.1,77636.7,100530.1])# Pa
gA  = np.array([0.14,4.82,18.19,40.31,71.31,111.50,160.50,219.47,292.49,379.04])             # mN/m
VV0 = np.array([1.0,0.9983,0.9935,0.9858,0.9752,0.9620,0.9465,0.9289,0.9084,0.8858])
svm95 = np.array([0,16,58,127,227,360,530,748,1024,1379.0])   # Pa
str95 = np.array([0,0.018,0.164,0.581,1.332,2.464,4.08,6.247,9.083,12.691])

MODEL, REF, REFEDGE = "#2166ac", "#fde0b8", "#e08214"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True})
fig, ax = plt.subplots(2, 4, figsize=(18, 8.5))
fig.suptitle("FF full-compartment MCF7 cell — quantitative results vs literature reality bands (2026-07-07)\n"
             "model = blue line/dots · amber band = experimental reference (NOT a fit). Native cortex Nc=494,802 + MT aster + 3000-bead nucleus.",
             fontsize=11, y=0.99)

def band(a, lo, hi, label):
    a.axhspan(lo, hi, color=REF, alpha=0.7, zorder=0); a.axhline(hi, color=REFEDGE, lw=1, ls="--", zorder=1)
    a.axhline(lo, color=REFEDGE, lw=1, ls="--", zorder=1)
    a.text(0.02, hi, label, transform=a.get_yaxis_transform(), fontsize=7.5, color="#9c4d00", va="bottom", ha="left")

# 1 force-indentation
a = ax[0,0]; a.semilogy(s*100, F, "-o", color=MODEL, lw=2, ms=5, label="model F_plate")
band(a, 10, 150, "real parallel-plate ~10–150 nN"); a.set_title("Force–indentation (AFM)"); a.set_xlabel("strain [%]"); a.set_ylabel("plate force [nN]  (log)"); a.legend(fontsize=7, loc="lower right")
a.annotate("~30–100× too STIFF", (20, 3920), (8, 400), color="#b30000", fontsize=8, arrowprops=dict(arrowstyle="->", color="#b30000"))
# 2 apparent cortical tension
a = ax[0,1]; a.semilogy(s*100, np.maximum(gA,1e-2), "-o", color=MODEL, lw=2, ms=5, label="model apparent γ")
band(a, 1, 10, "MCF7/interphase γ ~1–10 mN/m"); a.set_title("Apparent cortical tension"); a.set_xlabel("strain [%]"); a.set_ylabel("γ_apparent [mN/m]  (log)"); a.legend(fontsize=7, loc="lower right")
a.annotate("rest 0.14: ~70× too LOW\n(γ-floor)", (0, 0.14), (4, 0.5), color="#00408b", fontsize=8, arrowprops=dict(arrowstyle="->", color="#00408b"))
# 3 turgor dP
a = ax[0,2]; a.semilogy(s*100, np.maximum(dP,1), "-o", color=MODEL, lw=2, ms=5, label="model ΔP turgor")
band(a, 40, 500, "physiological rest ΔP ~40–500 Pa"); a.set_title("Osmotic turgor pressure"); a.set_xlabel("strain [%]"); a.set_ylabel("ΔP [Pa]  (log)"); a.legend(fontsize=7, loc="lower right")
# 4 V/V0
a = ax[0,3]; a.plot(s*100, VV0, "-o", color=MODEL, lw=2, ms=5); a.set_title("Volume (compaction)"); a.set_xlabel("strain [%]"); a.set_ylabel("V/V0"); a.axhline(1.0, color="#888", lw=1, ls=":")
# 5 σ_vm field p95
a = ax[1,0]; a.plot(s*100, svm95, "-o", color=MODEL, lw=2, ms=5); a.set_title("Cortex von-Mises σ_vm field (p95)"); a.set_xlabel("strain [%]"); a.set_ylabel("σ_vm p95 [Pa]")
a.text(0.5, 0.9, "median ≈ 0 (deviatoric γ-floor);\nconcentrates at contact", transform=a.transAxes, fontsize=7.5, color="#555", va="top")
# 6 areal strain p95
a = ax[1,1]; a.plot(s*100, str95, "-o", color=MODEL, lw=2, ms=5); a.set_title("Cortex areal strain field (p95)"); a.set_xlabel("strain [%]"); a.set_ylabel("areal strain p95")
a.text(0.05, 0.9, "localizes sharply at\nindentation contact", transform=a.transAxes, fontsize=7.5, color="#555", va="top")
# 7 compartment coupling ΔF
a = ax[1,2]; comps = ["nucleus", "membrane", "MT aster"]; dF = [45.26, 0.02, 0.04]
cols = ["#2166ac", "#9aa5b1", "#9aa5b1"]; a.bar(comps, dF, color=cols)
for i, v in enumerate(dF): a.text(i, v + 0.8, f"{v:.2f}%", ha="center", fontsize=8)
a.set_title("Compartment coupling (ΔF_plate @10% removing it)"); a.set_ylabel("|Δ plate force| [%]"); a.set_ylim(0, 52)
a.text(0.5, 0.6, "nucleus REAL (V_cyto);\nmembrane & MT negligible\nunder AFM turgor", transform=a.transAxes, fontsize=7.5, color="#555", ha="center")
# 8 model / reality ratio (how close?)
a = ax[1,3]; labels = ["rest γ\n(model/real)", "AFM force @20%\n(model/real)", "turgor ΔP @10%\n(model/real)"]
ratio = [0.14/10.0, 3920/80.0, 13187/300.0]  # <1 too low, >1 too high
rc = ["#00408b" if r < 1 else "#b30000" for r in ratio]
a.bar(labels, ratio, color=rc); a.axhline(1.0, color="#333", lw=1.5, ls="--"); a.set_yscale("log")
for i, r in enumerate(ratio): a.text(i, r*(1.3 if r>=1 else 0.6), f"{r:.2g}×", ha="center", fontsize=8, color=rc[i])
a.set_title("Model ÷ reality  (1× = match)"); a.set_ylabel("ratio  (log)")
a.text(0.5, 0.04, "blue<1 too low · red>1 too high", transform=a.transAxes, fontsize=7.5, color="#555", ha="center")

fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT_PNG, dpi=110, bbox_inches="tight")
print("wrote", OUT_PNG)

import os; os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["strain", "F_plate_nN", "dP_turgor_Pa", "gamma_apparent_mN_m", "V_over_V0", "svm_p95_Pa", "areal_strain_p95"])
    for i in range(len(s)): w.writerow([s[i], F[i], dP[i], gA[i], VV0[i], svm95[i], str95[i]])
print("wrote", OUT_CSV)
