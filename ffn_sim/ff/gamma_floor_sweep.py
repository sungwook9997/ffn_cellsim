"""γ-floor experiment — prestress sweep + quenched ensemble + figure (Stage 6d, ENGINE.md §4).

Runs the MD-free γ as a function of the myosin prestress f_myo (the controlled variable — NOT tuned
to a target), each point an ensemble of quenched realizations, and renders γ_active(f_myo) with the
Salbreux / MCF7 literature bands and the passive turgor (Young-Laplace) baseline overlaid. The plot
makes the central finding legible: at the lit-anchored NMIIA prestress the actomyosin γ sits
hundreds-of-× under the band (the γ-floor), reproducing the BAOAB-MD result by an independent
MD-free route.

Writes ``outputs/ff/figs/gamma_floor_sweep.png`` + ``outputs/ff/gamma_floor_sweep.npz``.
Run: ``python -m ffn_sim.ff.gamma_floor_sweep``.
"""

from __future__ import annotations

import os

import numpy as np

from ffn_sim.ff.gamma_estimator import MCF7_IQR_PN_UM, SALBREUX_BAND_PN_UM
from ffn_sim.ff.gamma_floor import (
    NMIIA_MINIFIL_STALL_PN,
    TURGOR_DP0,
    gamma_floor_ensemble,
)

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "outputs", "ff")


def run_sweep(f_myo_values=None, *, n_real: int = 8, n_filaments: int = 100, n_xl: int = 400,
              n_myo: int = 200, n_steps: int = 300, base_seed: int = 0, parallel: bool = True):
    """γ_active ensemble at each prestress in ``f_myo_values`` (default: 0 → 10× the lit anchor)."""
    if f_myo_values is None:
        a = NMIIA_MINIFIL_STALL_PN
        f_myo_values = [0.0, 0.2 * a, a, 2 * a, 4 * a, 10 * a, 30 * a]
    rows = []
    for f in f_myo_values:
        ens = gamma_floor_ensemble(f, n_real=n_real, n_filaments=n_filaments, n_xl=n_xl,
                                   n_myo=n_myo, n_steps=n_steps, base_seed=base_seed,
                                   parallel=parallel)
        active = np.array([r["gamma_active"] for r in ens["runs"]])
        rows.append({"f_myo": f, "active_mean": float(active.mean()),
                     "active_std": float(active.std()), "active_samples": active.tolist(),
                     "passive": float(ens["runs"][0]["gamma_passive"])})
        print(f"  f_myo={f:6.2f} pN : γ_active = {active.mean():.4f} ± {active.std():.4f} pN/µm "
              f"({active.mean() * 1e-3:.2e} mN/m)   γ_passive = {rows[-1]['passive']:.1f} pN/µm",
              flush=True)
    return rows


def render(rows, outdir: str = OUTDIR) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.join(outdir, "figs"), exist_ok=True)
    f = np.array([r["f_myo"] for r in rows])
    gm = np.array([r["active_mean"] for r in rows])
    gs = np.array([r["active_std"] for r in rows])
    passive = rows[0]["passive"]

    fig, ax = plt.subplots(figsize=(9, 6.5))
    # literature bands
    ax.axhspan(*SALBREUX_BAND_PN_UM, color="tab:green", alpha=0.15,
               label=f"Salbreux band {SALBREUX_BAND_PN_UM} pN/µm (0.35–0.65 mN/m)")
    ax.axhspan(*MCF7_IQR_PN_UM, color="tab:olive", alpha=0.20,
               label=f"MCF7 IQR {MCF7_IQR_PN_UM} pN/µm (0.18–0.40 mN/m)")
    ax.axhline(passive, color="tab:blue", ls="--",
               label=f"passive turgor γ=ΔP·R/2 = {passive:.0f} pN/µm (dP0={TURGOR_DP0:.0f} Pa)")
    # active sweep (per-realization spread + mean)
    for r in rows:
        ax.scatter([r["f_myo"]] * len(r["active_samples"]), r["active_samples"],
                   s=12, color="tab:red", alpha=0.35)
    ax.plot(f, gm, "-o", color="tab:red", lw=2, label="γ_active (actomyosin, method-of-planes)")
    ax.axvline(NMIIA_MINIFIL_STALL_PN, color="k", ls=":", alpha=0.7,
               label=f"NMIIA per-side stall = {NMIIA_MINIFIL_STALL_PN:.0f} pN (lit anchor)")
    ax.set_yscale("log")
    ax.set_xlabel("myosin prestress f_myo  [pN per link]  (controlled variable — swept, not tuned)")
    ax.set_ylabel("cortical tension γ  [pN/µm]   (1 pN/µm = 1e-3 mN/m)")
    ax.set_title("γ-floor experiment (MD-free, FF engine): actomyosin γ vs myosin prestress\n"
                 "actomyosin floored ~100–1000× under band at the lit prestress; turgor at-band "
                 "→ reproduces the BAOAB-MD γ-floor")
    ax.legend(fontsize=8, loc="center right")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    path = os.path.join(outdir, "figs", "gamma_floor_sweep.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def main():
    print("γ-floor prestress sweep (MD-free FF mechanical solve)...", flush=True)
    rows = run_sweep()
    os.makedirs(OUTDIR, exist_ok=True)
    npz = os.path.join(OUTDIR, "gamma_floor_sweep.npz")
    np.savez(npz, f_myo=[r["f_myo"] for r in rows],
             active_mean=[r["active_mean"] for r in rows],
             active_std=[r["active_std"] for r in rows],
             passive=[r["passive"] for r in rows],
             salbreux=SALBREUX_BAND_PN_UM, mcf7=MCF7_IQR_PN_UM,
             lit_anchor=NMIIA_MINIFIL_STALL_PN)
    fig = render(rows)
    print(f"wrote {npz}\nwrote {fig}", flush=True)


if __name__ == "__main__":
    main()
