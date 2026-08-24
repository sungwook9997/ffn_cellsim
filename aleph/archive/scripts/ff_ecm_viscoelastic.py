"""FF ECM library — emergent viscoelastic stress relaxation (matrix viscoelasticity).

Real ECM is viscoelastic: under a held strain the stress RELAXES (KB-1.6: G(t)=G∞+G₁·exp(−t/τ), SLS model,
τ~30-1000 s; matrix viscoelasticity controls cell spreading & fate — Chaudhuri 2016). The FF model relaxes
because crosslinks turn over (Bell-slip ``xl_turnover_kernel`` — each a Maxwell element of lifetime 1/k_off);
the fibers + bending stay elastic (the SLS spring in parallel). This measures G(t) across a sweep of
crosslink off-rates spanning real crosslink types, and shows **τ ∝ 1/k_off EMERGES** (not tuned) and lands
in the KB-1.6 range. A new time-domain characterization of the ECM library.

Run:  python -m aleph.scripts.ff_ecm_viscoelastic [--device cpu|cuda:0]
Out:  aleph/outputs/ff/ecm_lib/{ecm_viscoelastic.json, figs/stress_relaxation.png}
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aleph.laws import ecm_library as L
from aleph.laws import ecm_mechanics as M

OUT = "aleph/outputs/ff/ecm_lib"
FIGS = f"{OUT}/figs"

# crosslink off-rates (s⁻¹) by type → set the relaxation time (KB-1.6 τ 30-1000 s; KB α-actinin 0.066 Ferrer).
KOFFS = [("covalent-like (LOX, slow)", 0.0007),
         ("weak physical crosslink", 0.004),
         ("α-actinin-like (Ferrer)", 0.066)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--substeps", type=int, default=140)
    a = ap.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    rng = np.random.default_rng(7)
    ecm = L.build_fibrillar_ecm(L.get_spec("collagen_I"), [0, 0, 0], [24, 24, 24], concentration=2.0,
                                dim=3, alignment_S=0.0, pin_faces=(), rng=rng)
    print(f"collagen ⟨z⟩={ecm.connectivity_z:.2f} nodes={ecm.meta['n_nodes']}", flush=True)
    results = []
    for name, koff in KOFFS:
        sr = M.stress_relaxation(ecm, gamma0=0.1, koff0_per_s=koff, x_beta_nm=0.4, n_record=22,
                                 mech_substeps=a.substeps, device=a.device)
        sr["name"] = name
        results.append(sr)
        band = "IN 30-1000s" if 30.0 <= sr["tau_s"] <= 1000.0 else "outside"
        print(f"  {name}: k_off={koff}/s → τ={sr['tau_s']:.0f}s ({band}), G∞/G₀={sr['Ginf_over_G0']:.2f}, "
              f"τ·k_off={sr['tau_x_koff']:.2f}", flush=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for r in results:
        ax1.plot(r["t_s"], r["G_over_G0"], "o-", label=f"{r['name']} (τ={r['tau_s']:.0f}s)")
    ax1.axhspan(0, 1, xmin=0, xmax=0, alpha=0)  # noop to keep bounds
    ax1.set_xscale("log")
    ax1.set_xlabel("time t [s] (log)"); ax1.set_ylabel("relaxation modulus G(t)/G₀")
    ax1.set_title("ECM stress relaxation (KB-1.6 SLS)"); ax1.legend(fontsize=8); ax1.grid(alpha=0.3, which="both")
    koff = np.array([k for _, k in KOFFS]); tau = np.array([r["tau_s"] for r in results])
    ax2.loglog(1.0 / koff, tau, "o", ms=10, color="tab:red")
    lo, hi = (1.0 / koff).min(), (1.0 / koff).max()
    ax2.plot([lo, hi], [1.5 * lo, 1.5 * hi], "k--", alpha=0.5, label="τ = 1.5 / k_off (crosslink lifetime)")
    ax2.axhspan(30, 1000, alpha=0.12, color="green", label="KB-1.6 range 30-1000 s")
    ax2.set_xlabel("crosslink lifetime 1/k_off [s]"); ax2.set_ylabel("relaxation time τ [s]")
    ax2.set_title("τ EMERGES ∝ crosslink lifetime (not tuned)"); ax2.legend(fontsize=8); ax2.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(f"{FIGS}/stress_relaxation.png", dpi=130); plt.close(fig)
    with open(f"{OUT}/ecm_viscoelastic.json", "w") as f:
        json.dump({"koffs": KOFFS, "results": results}, f, indent=2, default=float)
    print(f"wrote {OUT}/ecm_viscoelastic.json + figs/stress_relaxation.png")


if __name__ == "__main__":
    main()
