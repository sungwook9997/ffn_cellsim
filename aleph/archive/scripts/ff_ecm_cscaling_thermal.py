"""FF ECM — c-scaling gap diagnosis: athermal vs analytic THERMAL semiflexible modulus (KB-anchored).

The collagen G'(c) exponent gap (FF athermal n≈1 vs literature n≈2.0-2.1) is now KB-diagnosed as a THERMAL
semiflexible effect, NOT a connectivity ⟨z⟩(c) effect:
  * KB-1.3 — collagen ⟨z⟩=3-3.5, sub-isostatic, roughly concentration-INDEPENDENT (no ⟨z⟩(c) growth datum) →
    the ⟨z⟩(c) route is ruled out; at fixed ⟨z⟩ the athermal network is density-linear G∝ρ∝c¹.
  * KB-1.30 benchmark #1 — the entropic (MacKintosh) plateau modulus G0 ~ kB·T·ℓp²/ξ⁵ (thermal bending of
    semiflexible segments between crosslinks), with KB-1.7 mesh ξ ~ c^(−1/2) (3D) ⇒ G0 ~ c^(5/2) = c^2.5.
  * KB-1.V.2.1 — the measured literature exponent is n≈2.0-2.1 (37 °C), BETWEEN the athermal c¹ and the pure
    affine-thermal c^2.5 — the real semiflexible crossover.

This module overlays the three on ONE figure (no new dynamics — the athermal points are the committed fixed-⟨z⟩
values, the thermal curve is the analytic MacKintosh scaling anchored at the reference concentration), quantifying
that the gap is thermal and bracketing the literature. Reaching c² in the RUNTIME needs the PI-gated thermal-WLC
force-extension (a production-default physics change) — this figure is the diagnosis, not that change.

Run:  python -m aleph.scripts.ff_ecm_cscaling_thermal
Out:  aleph/outputs/ff/ecm_lib/figs/cscaling_thermal_diagnosis.png + cscaling_thermal.json
"""

from __future__ import annotations

import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "aleph/outputs/ff/ecm_lib"

# Committed athermal G'(c) at FIXED physical ⟨z⟩=3.2 (CSCALING_REGIME_FINDING.md, box 26 µm) — measured, not tuned.
C_ATHERMAL = np.array([1.0, 2.0, 4.0, 7.0])
G_ATHERMAL = np.array([8.3, 18.5, 36.8, 64.1])
C_REF = 1.5                                   # reference concentration [mg/mL]
G_REF = 13.4                                  # measured G'(1.5 mg/mL) [Pa] (matches Yang-Kaufman ≈11)
LIT_N = (2.0, 2.1)                            # KB-1.V.2.1 literature exponent band
THERMAL_N = 2.5                               # KB-1.30 (G0~kBT ℓp²/ξ⁵) + KB-1.7 (ξ~c^-1/2) ⇒ c^(5/2)


def _fit_n(c, g):
    return float(np.polyfit(np.log(c), np.log(g), 1)[0])


def main():
    os.makedirs(f"{OUT}/figs", exist_ok=True)
    n_ath = _fit_n(C_ATHERMAL, G_ATHERMAL)
    cc = np.linspace(0.8, 8.0, 60)
    # anchor both power laws at the reference point (G_REF at C_REF)
    g_thermal = G_REF * (cc / C_REF) ** THERMAL_N
    g_lit_lo = G_REF * (cc / C_REF) ** LIT_N[0]
    g_lit_hi = G_REF * (cc / C_REF) ** LIT_N[1]
    g_ath_line = G_REF * (cc / C_REF) ** n_ath

    fig, ax = plt.subplots(figsize=(8.0, 5.8))
    ax.fill_between(cc, g_lit_lo, g_lit_hi, color="0.7", alpha=0.5, label=f"literature n=2.0-2.1 (KB-1.V.2.1)")
    ax.plot(cc, g_thermal, "-", color="#d62728", lw=2,
            label=f"analytic THERMAL MacKintosh  G0~c^{THERMAL_N:g} (KB-1.30·1.7)")
    ax.plot(cc, g_ath_line, "-", color="#1f77b4", lw=2, label=f"FF athermal fit  n={n_ath:.2f} (fixed ⟨z⟩=3.2)")
    ax.plot(C_ATHERMAL, G_ATHERMAL, "o", color="#1f77b4", ms=8, zorder=5, label="FF athermal measured (native REV)")
    ax.plot([C_REF], [G_REF], "*", color="k", ms=15, zorder=6, label=f"reference G'({C_REF})={G_REF} Pa (Yang-Kaufman)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("collagen concentration  c  [mg/mL]")
    ax.set_ylabel("shear modulus  G'  [Pa]")
    ax.set_title("Collagen G'(c) — the c-scaling gap is THERMAL, not connectivity (KB-anchored)\n"
                 f"athermal c^{n_ath:.2f} (⟨z⟩ fixed, KB-1.3 rules out ⟨z⟩(c)) → literature c^2 → thermal c^2.5; "
                 "runtime c² needs PI-gated thermal-WLC", fontsize=9)
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    path = f"{OUT}/figs/cscaling_thermal_diagnosis.png"
    fig.savefig(path, dpi=140); plt.close(fig)

    out = dict(n_athermal_fit=n_ath, thermal_n=THERMAL_N, literature_n=list(LIT_N),
               reasoning="KB-1.3 ⟨z⟩~3.2 c-independent rules out ⟨z⟩(c); KB-1.30 G0~kBT·ℓp²/ξ⁵ + KB-1.7 ξ~c^-1/2 "
                         "⇒ thermal c^2.5; literature (KB-1.V.2.1) c^2.0-2.1 sits between athermal c^1 and thermal "
                         "c^2.5 = the semiflexible crossover. Runtime resolution = PI-gated thermal-WLC.",
               athermal_points={float(c): float(g) for c, g in zip(C_ATHERMAL, G_ATHERMAL)})
    with open(f"{OUT}/cscaling_thermal.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"athermal fit n={n_ath:.2f} · thermal c^{THERMAL_N} · literature c^{LIT_N}")
    print(f"wrote {path} + cscaling_thermal.json")


if __name__ == "__main__":
    main()
