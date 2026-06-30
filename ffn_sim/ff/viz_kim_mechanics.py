"""Visualize the FF↔Kim absolute shear-modulus result (Stage 6c, task (c)).

The FF cross-linked actin network's elastic shear modulus G [Pa = pN/µm²], measured by affine simple
shear + interior relaxation (ff/kim_network.shear_modulus), with the actin backbone and the crosslink
junction as SEPARATELY-sourced stiffnesses (the prior code conflated them into one knob). Renders:

  (A) crosslink-limited regime (inextensible reshape actin, k_xl ≪ EA/L_seg): G is LINEAR in the
      crosslink stiffness k_xl and tracks the analytic affine form G ≈ f_na·k_xl·ρ_L·ℓc with a STABLE
      non-affine factor f_na (the robust, no-tuning analytic cross-check — analytic is primary).
  (B) the crosslinker-stiffness UNIT-SLIP impact: at the repo's broken k_xl=0.1 pN/µm the network is
      ~10⁵× too soft; at the lit-anchored α-actinin stiffness (Ferrer 2008, 4.6e5 pN/µm, with finite
      Kojima-1994 actin EA) the enthalpic G lands in the cortex-stiff range — overlaid with the
      in-vitro Kim/Gardel entropic band and the in-vivo cortex band.

Writes outputs/ff/figs/kim_shear_modulus.png. Run: python -m ffn_sim.ff.viz_kim_mechanics
"""

from __future__ import annotations

import os

import numpy as np

from ffn_sim.ff.kim_network import EA_ACTIN_PN, shear_modulus

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "ff", "figs")


def render(outdir: str = OUTDIR) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # (A) crosslink-limited linearity (reshape / rigid actin), with analytic affine form
    kxls = np.array([1.0, 3.0, 10.0, 30.0, 100.0])
    G, Gana = [], []
    for k in kxls:
        g, _, m = shear_modulus(C_A_uM=300.0, R_acp=1.5, k_xl=k, n_steps=6000, seed=0)
        G.append(g); Gana.append(m["G_analytic_crosslink"])
    G = np.array(G); Gana = np.array(Gana)
    f_na = float(np.mean(G / Gana))                         # stable non-affine factor

    # (B) regime points — ONE consistent model (reshape / rigid actin) at C_A=300, R=1.5; broken,
    # default, and lit-anchored α-actinin. (Reshape gives the rigid-actin UPPER bound at high k_xl;
    # the finite-EA solver lowers the sourced point ~2× — annotated. broken/default are exact: there
    # k_xl ≪ EA/L_seg so actin is genuinely rigid.)
    G_broken, _, _ = shear_modulus(C_A_uM=300.0, R_acp=1.5, k_xl=0.1, n_steps=6000, seed=0)
    G_default, _, _ = shear_modulus(C_A_uM=300.0, R_acp=1.5, k_xl=10.0, n_steps=6000, seed=0)
    G_sourced, _, _ = shear_modulus(C_A_uM=300.0, R_acp=1.5, k_xl=4.6e5, n_steps=20000, seed=0)

    os.makedirs(outdir, exist_ok=True)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.4))

    axA.loglog(kxls, G, "o", color="navy", ms=8, label="FF measured (reshape actin)")
    axA.loglog(kxls, f_na * Gana, "-", color="grey", lw=2,
               label=f"analytic affine ×f_na  (f_na={f_na:.3f}, stable)")
    axA.set_xlabel("crosslink stiffness k_xl  [pN/µm]"); axA.set_ylabel("shear modulus G  [Pa]")
    axA.set_title("Crosslink-limited regime — G LINEAR in k_xl\n"
                  "matches analytic form G ≈ f_na·k_xl·ρ_L·ℓc (non-affine factor stable)")
    axA.legend(fontsize=9); axA.grid(True, which="both", alpha=0.3)

    labels = ["broken\nk_xl=0.1\n(repo)", "knob\nk_xl=10\n(old default)",
              "sourced α-actinin\nk_xl=4.6e5\n(Ferrer; UB)"]
    vals = [max(G_broken, 1e-3), max(G_default, 1e-3), G_sourced]
    colors = ["crimson", "darkorange", "seagreen"]
    xpos = np.arange(3)
    axB.bar(xpos, vals, color=colors, alpha=0.85, width=0.6)
    axB.set_yscale("log"); axB.set_xticks(xpos); axB.set_xticklabels(labels, fontsize=8.5)
    axB.set_ylabel("shear modulus G  [Pa]")
    axB.axhspan(0.1, 1000.0, color="steelblue", alpha=0.15, label="in-vitro Kim/Gardel (entropic) 0.1–1000 Pa")
    axB.axhspan(100.0, 1000.0, color="seagreen", alpha=0.18, label="in-vivo cortex ~0.1–1 kPa")
    for xi, v in zip(xpos, vals):
        axB.text(xi, v * 1.3, f"{v:.2g} Pa", ha="center", fontsize=9)
    axB.set_title("Crosslinker-stiffness UNIT-SLIP impact\n"
                  "broken 0.1 pN/µm ⇒ cortex ~10⁵× too soft; sourced ⇒ cortex-stiff enthalpic")
    axB.legend(fontsize=8, loc="lower right")

    fig.suptitle("FF↔Kim absolute shear modulus (Stage 6c) — actin (inextensible/EA) + crosslink (k_xl) "
                 "separately sourced; analytic primary", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = os.path.join(outdir, "kim_shear_modulus.png")
    fig.savefig(path, dpi=135); plt.close(fig)
    print(f"wrote {path}")
    print(f"  crosslink-limited: G(k_xl) = {G.tolist()} Pa for k_xl={kxls.tolist()} (non-affine factor f_na={f_na:.4f})")
    print(f"  regime: broken(0.1)={G_broken:.3f} Pa | default(10)={G_default:.2f} Pa | sourced α-actinin+EA={G_sourced:.1f} Pa")
    return path


if __name__ == "__main__":
    render()
