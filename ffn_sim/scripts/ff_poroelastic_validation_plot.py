"""Piece #2 validation figure — Hertz E_eff, undrained → drained, vs the measured MCF7 249 Pa band.

Data from ff_hertz_validation.py (N_fil=2000 CPU, small-strain 2–3%, Sneddon-parabolic inversion, R_cell=7.5µm).
Shows: (A) the F vs δ1^1.5 Hertz fit for each regime with the analytic real-MCF7 (E=249 Pa) line overlaid, and
(B) the fitted apparent modulus E_fit by regime vs the Zbiral band [224,279] Pa. Honest finding: the biphasic
cytoplasm (drainage + drained solid) cuts the undrained 16.5× overshoot to ~5×, but at Zbiral's FAST rate (5µm/s,
τ_load≈0.1s ≪ τ_osm≈30–200s) the cell cannot drain → the model stays ~undrained, and the residual ~5× is the
CORTEX indentation stiffness (the lateral-bulge / crosslinker-stiffness lever), NOT the cytoplasm.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/Users/sw1/ffn_cellsim/ffn_sim/outputs/ff/figs/ff_poroelastic_validation.png"
R, NU = 7.5, 0.5
d1 = np.array([0.150, 0.1875, 0.225])                 # per-contact indentation [µm] at strain 2/2.5/3%
x = d1 ** 1.5
# measured plate force [pN] per regime
F = {
    "undrained (fast rate ≈ Zbiral 5µm/s)": np.array([1187.3, 1598.7, 2142.5]),
    "drained: setpoint 40Pa + K_drained": np.array([596.9, 643.9, 692.3]),
    "drained: Lp + slow ramp + K_drained": np.array([452.0, 484.2, 522.1]),
}
E_fit = {"undrained (fast rate ≈ Zbiral 5µm/s)": 4114.0,
         "drained: setpoint 40Pa + K_drained": 1548.0,
         "drained: Lp + slow ramp + K_drained": 1167.0}
E_TARGET, E_BAND = 249.0, (224.0, 279.0)
COL = {"undrained (fast rate ≈ Zbiral 5µm/s)": "#b2182b",
       "drained: setpoint 40Pa + K_drained": "#5aae61",
       "drained: Lp + slow ramp + K_drained": "#2166ac"}
FILL, EDGE = "#fde0b8", "#e08214"

def F_real(E, d):   # analytic parabolic Hertz for a real cell of modulus E
    return (4.0/3.0) * (E/(1.0-NU**2)) * np.sqrt(R) * d**1.5

plt.rcParams.update({"font.size": 9.5, "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True})
fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4))
fig.suptitle("FF biphasic-cytoplasm validation — Hertz apparent modulus vs measured MCF7 (Zbiral 249 Pa) · 2026-07-08\n"
             "undrained 16.5× → drained ~5×; residual ~5× is the CORTEX (lateral-bulge/crosslinker), not cytoplasm. At Zbiral's fast rate the cell can't drain (τ_load≈0.1s ≪ τ_osm≈30–200s).",
             fontsize=10.5, y=0.99)

# A — F vs δ1^1.5 (Hertz fit); real MCF7 line + band
a = ax[0]
dd = np.linspace(0, d1.max()*1.05, 50)
a.fill_between(dd**1.5, F_real(E_BAND[0], dd), F_real(E_BAND[1], dd), color=FILL, alpha=0.8, zorder=0,
               label=f"real MCF7 band {E_BAND[0]:.0f}–{E_BAND[1]:.0f} Pa")
a.plot(dd**1.5, F_real(E_TARGET, dd), "-", color=EDGE, lw=2, label="real MCF7 E=249 Pa (Zbiral)")
for k, f in F.items():
    a.plot(x, f, "-o", color=COL[k], lw=2, ms=6, label=k)
a.set_title("Hertz fit: plate force vs δ¹·⁵"); a.set_xlabel("δ₁^1.5  [µm^1.5]  (δ₁ = per-contact indent)")
a.set_ylabel("plate force  [pN]"); a.legend(fontsize=7.3, loc="upper left")

# B — fitted apparent modulus by regime vs band
a = ax[1]
labels = list(E_fit.keys()); vals = [E_fit[k] for k in labels]
a.axhspan(E_BAND[0], E_BAND[1], color=FILL, alpha=0.8, zorder=0)
a.axhline(E_TARGET, color=EDGE, lw=1.5, ls="--", zorder=1)
a.text(2.5, E_TARGET, "MCF7 249 Pa (224–279)", fontsize=8, color="#9c4d00", va="bottom", ha="right")
bars = a.bar(range(len(labels)), vals, color=[COL[k] for k in labels], width=0.6)
for i, v in enumerate(vals):
    a.text(i, v*1.03, f"{v:.0f} Pa\n{v/E_TARGET:.1f}×", ha="center", fontsize=8.5)
a.set_yscale("log"); a.set_ylim(150, 6000)
a.set_xticks(range(len(labels))); a.set_xticklabels(["undrained\n(fast=Zbiral)", "drained\nsetpoint", "drained\nLp+slow"], fontsize=8)
a.set_title("Apparent modulus E_fit by regime  (log)"); a.set_ylabel("E_fit  [Pa]  (log)")

fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(OUT, dpi=120, bbox_inches="tight")
print("wrote", OUT)
