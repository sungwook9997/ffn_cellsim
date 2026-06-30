"""Visualize the on-device Hand KMC turnover (Stage 6k) — the γ-floor turnover-robustness axis.

Sweeps the per-myosin contractile load and, for each, runs ``network_warp.simulate_turnover_on_device``
(Bell detach + Poisson re-attach on the GPU) → the steady ENGAGED fraction + the actomyosin γ measured
on the engaged subset. Renders:

  (A) engaged myosin fraction vs load — the simulated KMC steady state overlaid with the analytic Bell
      law k_on/(k_on+p₀·exp(f/f₀)); the duty-ratio self-limit (engaged↓ as load↑).
  (B) actomyosin γ vs load with the Salbreux band — γ stays ~10³× under band at EVERY load: turnover
      cannot lift the floor (at physiological load engaged≈1, no headroom; higher load only detaches
      motors, lowering the engaged density). The γ-floor is a force-magnitude limit.

Device-agnostic (GPU↔CPU parity); the figure is identical on the A5000. Writes
``outputs/ff/figs/turnover_gamma_floor.png``.

Run: ``python -m ffn_sim.ff.viz_turnover`` (``--n 300 --device cuda:0`` on gbook).
"""

from __future__ import annotations

import os

import numpy as np

from ffn_sim.ff.gamma_estimator import SALBREUX_BAND_PN_UM
from ffn_sim.ff.gamma_floor import (
    CortexParams,
    build_crosslinked_cortex,
    equilibrate,
    measure_gamma,
)
from ffn_sim.ff.hand_kmc import NMIIA_MYOSIN, bell_off_rate
from ffn_sim.ff.network_warp import simulate_turnover_on_device

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "ff", "figs")


def _bell_steady(f):
    return NMIIA_MYOSIN.k_on / (NMIIA_MYOSIN.k_on + bell_off_rate(f, NMIIA_MYOSIN.p0, NMIIA_MYOSIN.f0))


def render(n_filaments: int = 300, loads=(1.0, 5.0, 10.0, 20.0, 30.0, 40.0), n_steps: int = 12000,
           seed: int = 0, device: str = "cpu", outdir: str = OUTDIR) -> str:
    import copy

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    engaged, gam = [], []
    for f in loads:
        cx = build_crosslinked_cortex(CortexParams(), n_filaments=n_filaments, n_xl=n_filaments,
                                      n_myo=max(10, n_filaments // 4), rng=np.random.default_rng(seed))
        equilibrate(cx, 0.0, n_steps=300, method="device", device=device)
        _, bmask, _ = simulate_turnover_on_device(cx, f, n_steps=n_steps, kmc_every=50,
                                                  tau_kmc=0.01, device=device)
        cxe = copy.copy(cx); cxe.myo_i = cx.myo_i[bmask]; cxe.myo_j = cx.myo_j[bmask]
        engaged.append(float(bmask.mean()))
        gam.append(measure_gamma(cxe, f, turgor=False)["gamma_active"])
    loads = np.array(loads); engaged = np.array(engaged); gam = np.array(gam)
    ff = np.linspace(loads.min(), loads.max(), 200)

    os.makedirs(outdir, exist_ok=True)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.2))

    axA.plot(ff, _bell_steady(ff), "-", color="grey", lw=2,
             label="analytic Bell  k_on/(k_on+p₀e^{f/f₀})")
    axA.plot(loads, engaged, "o", color="crimson", ms=8, label="KMC (on-device)")
    axA.set_xlabel("per-myosin load f [pN]"); axA.set_ylabel("engaged fraction")
    axA.set_ylim(0, 1.05)
    axA.set_title("Myosin Hand turnover — engaged fraction self-limits with load\n"
                  "(NF2007 §10.1 Bell; on-device KMC matches the analytic steady state)")
    axA.legend(fontsize=9)

    band = SALBREUX_BAND_PN_UM
    axB.axhspan(band[0] * 1e-3, band[1] * 1e-3, color="seagreen", alpha=0.2,
                label=f"Salbreux band {band[0]*1e-3:.2f}–{band[1]*1e-3:.2f} mN/m")
    axB.semilogy(loads, gam * 1e-3, "s-", color="navy", lw=2, ms=7, label="γ_active (turnover, engaged)")
    axB.set_xlabel("per-myosin load f [pN]"); axB.set_ylabel("γ_active  [mN/m]")
    axB.set_title("Actomyosin γ stays ~10³× under band at every load\n"
                  "turnover cannot lift the force-magnitude floor")
    axB.legend(fontsize=9)

    fig.suptitle("FF γ-floor — TURNOVER robustness axis (Stage 6k, on-device Hand KMC)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = os.path.join(outdir, "turnover_gamma_floor.png")
    fig.savefig(path, dpi=135)
    plt.close(fig)
    print(f"wrote {path}")
    print("  loads[pN]:", list(loads))
    print("  engaged  :", [f"{e:.3f}" for e in engaged])
    print("  γ[mN/m]  :", [f"{g*1e-3:.2e}" for g in gam])
    return path


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--steps", type=int, default=12000)
    ap.add_argument("--device", type=str, default="cpu")
    a = ap.parse_args()
    render(n_filaments=a.n, n_steps=a.steps, device=a.device)
