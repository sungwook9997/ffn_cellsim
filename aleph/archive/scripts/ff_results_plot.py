"""FF full-compartment quantitative results — record + plot vs REAL measured MCF7 AFM curves.

The reality reference is now a measured MCF7 force–indentation / stress–strain LAW (Hertz with the
measured Young's modulus), not a vague band, so at every indentation you read model-vs-real directly.

Real MCF7 anchors (LITERATURE, KB-verified — NOT fits):
  * whole-cell modulus, 10 um colloidal bead (matched large-contact geometry): E ~ 249 Pa (Zbiral 2023; KB-6.1.1)
  * sharp-tip adherent MCF7: E ~ 0.2-1.0 kPa (Li 2008 BBRC 374:609; KB-6.1.1)
  * cytoplasm: G' ~ 33 Pa, eta ~ 56 Pa s (Dessard 2024 Nanoscale Adv; KB-3.B3.1)
  * resting cortical tension: interphase ~0.3-1 mN/m; MCF7 suspended ~10 mN/m (Moazzeni 2021)
Hertz spherical: F(d) = (4/3)(E/(1-nu^2)) sqrt(R) d^1.5 ; nu=0.5 ; R = 7.5 um (cell/probe).
Effective apparent stress: sigma = F/(pi R^2). Linear small-strain reference: sigma = E * eps. Log axes noted.
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_PNG = "/Users/sw1/ffn_cellsim/aleph/outputs/ff/figs/ff_fullcompartment_results.png"
OUT_CSV = "/Users/sw1/ffn_cellsim/aleph/outputs/ff/data/ff_afm_results_2026-07-07.csv"

# ---- model data (this session's native AFM full-compartment sweep) ----
s   = np.array([0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27])
F   = np.array([0.01,13.93,102.98,331.97,800.52,1580.72,2805.57,4560.76,6974.72,10185.07])  # nN
dP  = np.array([36.4,1283.7,4841.8,10724.7,18966.7,29645.7,42654.7,58294.1,77636.7,100530.1])# Pa
gA  = np.array([0.14,4.82,18.19,40.31,71.31,111.50,160.50,219.47,292.49,379.04])             # mN/m
VV0 = np.array([1.0,0.9983,0.9935,0.9858,0.9752,0.9620,0.9465,0.9289,0.9084,0.8858])
svm95 = np.array([0,16,58,127,227,360,530,748,1024,1379.0])   # Pa

# ---- real MCF7 references ----
R   = 7.5e-6         # m
nu  = 0.5
E_bead = 249.0       # Pa  (Zbiral 2023, 10um colloidal, whole-cell)
E_lo, E_hi = 200.0, 1000.0   # Pa (sharp-tip adherent MCF7 band, Li 2008)
d_m = s * R          # indentation depth ~ strain*R  (m)
def hertz(E, d):     # N
    return (4.0/3.0) * (E/(1.0-nu**2)) * np.sqrt(R) * np.power(np.clip(d,0,None), 1.5)
F_real_bead = hertz(E_bead, d_m) * 1e9       # nN
F_real_lo   = hertz(E_lo,   d_m) * 1e9
F_real_hi   = hertz(E_hi,   d_m) * 1e9

sig_model = F * 1e-9 / (np.pi * R**2)         # Pa  (apparent stress F/piR^2)
sig_bead  = E_bead * s                        # Pa  (linear sigma = E*eps)
sig_lo, sig_hi = E_lo * s, E_hi * s

# ---- figure: 2x2, all about how-hard-you-press -> how-much-stress vs REAL MCF7 ----
MODEL, REAL, RLO, FILL = "#2166ac", "#b2182b", "#e08214", "#fde0b8"
plt.rcParams.update({"font.size": 9.5, "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True})
fig, ax = plt.subplots(2, 2, figsize=(13.5, 10.2))
fig.suptitle("FF full-compartment MCF7 cell — press-vs-stress against MEASURED MCF7 AFM law (2026-07-07)\n"
             "blue = model · red = real MCF7 E=249 Pa (Zbiral 2023 colloidal, matched geometry) · amber = MCF7 sharp-tip E=0.2–1.0 kPa (Li 2008). Native cortex Nc=494,802 + MT aster + 3000-bead nucleus.",
             fontsize=11, y=0.995)

def ratio_txt(a, xi, ym, yr, xytext):
    a.annotate(f"{ym/yr:,.0f}× stiffer", (xi, ym), xytext, color="#7a0000", fontsize=8.5,
               arrowprops=dict(arrowstyle="->", color="#7a0000", lw=1.2))

# A — force vs indentation depth
a = ax[0,0]
a.fill_between(d_m*1e6, F_real_lo, F_real_hi, color=FILL, alpha=0.8, zorder=0, label="real MCF7 sharp-tip 0.2–1.0 kPa")
a.semilogy(d_m*1e6, np.maximum(F_real_bead,1e-4), "-", color=REAL, lw=2.2, label="real MCF7 E=249 Pa (colloidal)")
a.semilogy(d_m*1e6, np.maximum(F,1e-4), "-o", color=MODEL, lw=2.2, ms=5, label="model F_plate")
a.set_title("Force vs indentation depth"); a.set_xlabel("indentation depth δ ≈ strain·R  [µm]")
a.set_ylabel("force  [nN]   (log)"); a.legend(fontsize=7.5, loc="lower right")
i = 4  # 12% strain
ratio_txt(a, d_m[i]*1e6, F[i], F_real_bead[i], (0.15, 30))

# B — apparent stress vs strain
a = ax[0,1]
a.fill_between(s*100, sig_lo, sig_hi, color=FILL, alpha=0.8, zorder=0, label="real MCF7 sharp-tip 0.2–1.0 kPa")
a.semilogy(s*100, np.maximum(sig_bead,1e-2), "-", color=REAL, lw=2.2, label="real MCF7 E=249 Pa (σ=E·ε)")
a.semilogy(s*100, np.maximum(sig_model,1e-2), "-o", color=MODEL, lw=2.2, ms=5, label="model σ=F/πR²")
a.set_title("Apparent stress vs compressive strain"); a.set_xlabel("compressive strain  [%]")
a.set_ylabel("apparent stress  [Pa]   (log)"); a.legend(fontsize=7.5, loc="lower right")
ratio_txt(a, s[i]*100, sig_model[i], sig_bead[i], (1, 12))

# C — apparent cortical tension vs strain
a = ax[1,0]
a.axhspan(0.3, 1.0, color=FILL, alpha=0.8, zorder=0)
a.axhline(10.0, color=RLO, lw=1.4, ls="--", zorder=1)
a.text(0.5, 10.0, "MCF7 suspended ~10 mN/m (Moazzeni 2021)", fontsize=7.5, color="#9c4d00", va="bottom")
a.text(0.5, 1.0, "interphase 0.3–1 mN/m", fontsize=7.5, color="#9c4d00", va="bottom")
a.semilogy(s*100, np.maximum(gA,1e-2), "-o", color=MODEL, lw=2.2, ms=5, label="model apparent γ")
a.set_title("Apparent cortical tension vs strain"); a.set_xlabel("compressive strain  [%]")
a.set_ylabel("γ_apparent  [mN/m]   (log)"); a.legend(fontsize=7.5, loc="lower right")
a.annotate("rest 0.14: ~70× too LOW\n(γ-floor, no myosin generation)", (0, 0.14), (2, 0.4),
           color="#00408b", fontsize=8, arrowprops=dict(arrowstyle="->", color="#00408b"))

# D — internal cortex von-Mises stress field p95 vs strain
a = ax[1,1]
a.plot(s*100, svm95, "-o", color=MODEL, lw=2.2, ms=5, label="model σ_vm field p95")
a.set_title("Internal cortex von-Mises stress field (p95)"); a.set_xlabel("compressive strain  [%]")
a.set_ylabel("σ_vm p95  [Pa]"); a.legend(fontsize=7.5, loc="upper left")
a.text(0.97, 0.05, "median ≈ 0 (deviatoric γ-floor)\np95 concentrates at contact rim\n(FEM field — see HTML viewer)",
       transform=a.transAxes, fontsize=8, color="#555", va="bottom", ha="right")

fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUT_PNG, dpi=120, bbox_inches="tight")
print("wrote", OUT_PNG)

# ---- machine-readable data (adds the real MCF7 references at each indentation) ----
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["strain", "indent_depth_um", "F_model_nN", "F_realMCF7_E249_nN", "F_realMCF7_0.2kPa_nN",
                "F_realMCF7_1.0kPa_nN", "sigma_model_Pa", "sigma_realMCF7_E249_Pa",
                "dP_turgor_Pa", "gamma_apparent_mN_m", "V_over_V0", "svm_p95_Pa"])
    for j in range(len(s)):
        w.writerow([s[j], d_m[j]*1e6, F[j], F_real_bead[j], F_real_lo[j], F_real_hi[j],
                    sig_model[j], sig_bead[j], dP[j], gA[j], VV0[j], svm95[j]])
print("wrote", OUT_CSV)

# ---- console: the press->stress table the PI asked for ----
print("\n  strain  depth[um]  F_model[nN]  F_real249[nN]  ratio   sig_model[Pa]  sig_real249[Pa]  ratio")
for j in range(len(s)):
    rf = F[j]/F_real_bead[j] if F_real_bead[j] > 0 else float("nan")
    rs = sig_model[j]/sig_bead[j] if sig_bead[j] > 0 else float("nan")
    print(f"  {s[j]*100:4.0f}%   {d_m[j]*1e6:6.2f}   {F[j]:10.1f}   {F_real_bead[j]:11.3f}  {rf:6.0f}x "
          f"  {sig_model[j]:11.0f}   {sig_bead[j]:13.1f}  {rs:6.0f}x")
