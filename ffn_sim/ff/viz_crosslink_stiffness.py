"""Visualize the crosslink-stiffness γ-floor finding (PI surface, 2026-06-30).

PI approved correcting the crosslinker stiffness link_k 0.1 → Ferrer-2008 4.6e5 pN/µm (α-actinin). This
figure shows what that correction does to the cortical-tension channels, measured by the method-of-
planes at the resting geometry + myosin modulator (gamma_floor.measure_gamma), broken vs corrected:

  * γ_active / γ_myo (the ACTIVE myosin-generated channel) stays FLOORED at ~0.15 pN/µm = 1.5e-4 mN/m
    at BOTH stiffnesses — i.e. the γ-floor is force-magnitude-limited, ROBUST to a 10⁶× crosslink-
    stiffness change. This CLOSES the last root-cause candidate (cross-bridge softness).
  * γ_xl / γ_total (the PASSIVE crosslink elastic channel) is ~0 at the broken value but becomes a
    large SPURIOUS number at the corrected value — stiff crosslinks amplify the tiny residual
    relaxation stretch (~0.6 nm) into ~450 pN tension, and the relaxer occasionally over-stretches
    (numerical instability). This is NOT active tension; it is why the stiff value cannot go into the
    production cortex relaxation until a robust stiff-crosslink relaxer exists.

Writes outputs/ff/figs/crosslink_stiffness_gamma.png. Run: python -m ffn_sim.ff.viz_crosslink_stiffness
"""

from __future__ import annotations

import os

import numpy as np

from ffn_sim.ff.gamma_estimator import SALBREUX_BAND_PN_UM
from ffn_sim.ff.gamma_floor import (
    NMIIA_MINIFIL_STALL_PN,
    CortexParams,
    build_crosslinked_cortex,
    equilibrate,
    measure_gamma,
)

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "ff", "figs")


def _channels(k_corrected, seeds=(0, 1, 2, 3)):
    myo, xl = [], []
    for s in seeds:
        cx = build_crosslinked_cortex(CortexParams(), n_filaments=1000, n_xl=1000, n_myo=100,
                                      rng=np.random.default_rng(s))
        if k_corrected is not None:
            rng2 = np.random.default_rng(1000 + s); is_a = rng2.random(cx.xl_i.size) < 0.30
            cx.xl_k = np.where(is_a, k_corrected[0], k_corrected[1])
        equilibrate(cx, 0.0, n_steps=600, turgor=False, method="implicit")
        g = measure_gamma(cx, NMIIA_MINIFIL_STALL_PN, turgor=False)
        myo.append(g.get("gamma_myo", 0.0)); xl.append(g.get("gamma_xl", 0.0))
    return np.array(myo), np.array(xl)


def render(outdir: str = OUTDIR) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    myo_b, xl_b = _channels(None)                      # broken 0.1
    myo_c, xl_c = _channels((4.6e5, 8.2e5))            # Ferrer corrected

    os.makedirs(outdir, exist_ok=True)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.4))
    band = SALBREUX_BAND_PN_UM

    # (A) ACTIVE myosin channel — floored at both (robust)
    axA.axhspan(band[0], band[1], color="seagreen", alpha=0.2, label=f"Salbreux band {band[0]:.0f}–{band[1]:.0f} pN/µm")
    axA.bar([0, 1], [myo_b.mean(), myo_c.mean()], yerr=[myo_b.std(), myo_c.std()],
            color=["crimson", "seagreen"], alpha=0.85, width=0.5, capsize=5)
    axA.set_yscale("log"); axA.set_xticks([0, 1]); axA.set_xticklabels(["broken\nk_xl=0.1", "Ferrer\nk_xl=4.6e5"])
    axA.set_ylabel("γ_active (myosin channel)  [pN/µm]")
    for xi, v in zip([0, 1], [myo_b.mean(), myo_c.mean()]):
        axA.text(xi, v * 1.5, f"{v:.2f}", ha="center", fontsize=9)
    axA.set_title("ACTIVE γ (myosin) — FLOORED at BOTH stiffnesses\n"
                  "robust to 10⁶× crosslink-stiffness change (force-magnitude limit)")
    axA.legend(fontsize=8)

    # (B) PASSIVE crosslink channel — ~0 broken, spurious + unstable corrected
    axB.bar([0, 1], [abs(xl_b.mean()) + 1e-3, abs(xl_c.mean())], yerr=[xl_b.std(), xl_c.std()],
            color=["crimson", "darkorange"], alpha=0.85, width=0.5, capsize=5)
    axB.set_yscale("log"); axB.set_xticks([0, 1]); axB.set_xticklabels(["broken\nk_xl=0.1", "Ferrer\nk_xl=4.6e5"])
    axB.set_ylabel("γ_xl (passive crosslink channel)  [pN/µm]")
    axB.axhspan(band[0], band[1], color="seagreen", alpha=0.2)
    axB.set_title("PASSIVE γ_xl — spurious + UNSTABLE at the stiff value\n"
                  "(stiff crosslinks amplify residual relax-stretch; occasional blow-up)")

    fig.suptitle("Crosslink-stiffness correction (PI 2026-06-30): the ACTIVE γ-floor is robust; the "
                 "stiff PASSIVE channel is a relaxation artifact", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = os.path.join(outdir, "crosslink_stiffness_gamma.png")
    fig.savefig(path, dpi=135); plt.close(fig)
    print(f"wrote {path}")
    print(f"  γ_active(myo): broken {myo_b.mean():.3f}±{myo_b.std():.3f} | corrected {myo_c.mean():.3f}±{myo_c.std():.3f} pN/µm (FLOORED both)")
    print(f"  γ_xl(passive): broken {xl_b.mean():.3f} | corrected {xl_c.mean():.1f}±{xl_c.std():.1f} pN/µm (spurious+unstable)")
    return path


if __name__ == "__main__":
    render()
