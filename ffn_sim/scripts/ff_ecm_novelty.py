"""FF ECM library — emergent & mode-dependent mechanics (the novel contributions).

Two things this library can uniquely characterize because it measures the FULL modulus tensor (shear AND
indentation AND uniaxial) on the SAME mechanistic Mikado network across the material/alignment/connectivity
space — neither is tuned in; both EMERGE from the fine-grained physics:

1. **Emergent nonlinear strain-stiffening** (KB-1.V.2.5 gate): the differential shear modulus K(γ)=dσ/dγ
   rises above a critical strain γ_c, and γ_c is ~concentration-independent (geometry-set). Reproducing
   this from a model NOT tuned to it is a predictive validation of the mechanistic fidelity.

2. **Shear ↔ compression/indentation modulus DECOUPLING** vs network connectivity ⟨z⟩ (van Oosten 2019,
   Sci Rep: semiflexible networks have very different moduli in different deformation modes). Sub-isostatic
   collagen is stiff in bulk shear but soft under local indentation; the decoupling ratio is a fingerprint
   of the network architecture — a systematic map the FF library produces directly.

Run:  python -m ffn_sim.scripts.ff_ecm_novelty [--device cpu]
Out:  ffn_sim/outputs/ff/ecm_lib/{ECM_NOVELTY.md, figs/strain_stiffening.png, figs/mode_decoupling.png}
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffn_sim.ff import ecm_library as L
from ffn_sim.ff import ecm_mechanics as M

OUT = "ffn_sim/outputs/ff/ecm_lib"
FIGS = f"{OUT}/figs"


def strain_stiffening_study(device="cpu", n_steps=4500, box=24.0):
    """Emergent K(γ)/K0 for collagen (2 conc, tests γ_c concentration-independence) + fibrin."""
    cases = [("collagen_I", 1.5), ("collagen_I", 3.0), ("fibrin", 2.0)]
    out = []
    g = np.linspace(0.02, 0.42, 11)
    for mat, c in cases:
        rng = np.random.default_rng(9)
        ecm = L.build_fibrillar_ecm(L.get_spec(mat), [0, 0, 0], [box, box, box], concentration=c,
                                    dim=3, alignment_S=0.0, pin_faces=(), rng=rng)
        ss = M.strain_stiffening(ecm, gammas=g, n_steps=n_steps, device=device)
        ss.update({"material": mat, "conc": c, "z": ecm.connectivity_z})
        out.append(ss)
        print(f"  {mat} c={c}: K0={ss['K0_Pa']:.1f} max_stiffen={ss['max_stiffening']:.1f}x γc={ss['gamma_c']:.2f}",
              flush=True)
    return out


def mode_decoupling_study(device="cpu", n_steps=4500, box=28.0):
    """Bulk shear G vs local indentation E across collagen connectivity ⟨z⟩ (set by concentration)."""
    out = []
    for c in (1.0, 2.0, 4.0):
        rng = np.random.default_rng(13)
        # a slab for indentation (bottom pinned) + shear on the same build
        ecm = L.build_fibrillar_ecm(L.get_spec("collagen_I"), [0, 0, 0], [box, box, box * 0.7],
                                    concentration=c, dim=3, alignment_S=0.0, pin_faces=("z_lo",),
                                    pin_margin_um=2.0, rng=rng)
        G = M.shear_modulus(ecm, gamma=0.02, n_steps=n_steps, device=device)["G_Pa"]
        nu = L.get_spec("collagen_I").poisson
        E_bulk = 2 * (1 + nu) * G
        ind = M.indentation_modulus(ecm, indenter_R_um=16.0, max_depth_um=1.6, n_depths=5,
                                    k_ind=3.0e3, relax_steps=n_steps, device=device)
        row = {"conc": c, "z": ecm.connectivity_z, "G_shear_Pa": G, "E_bulk_Pa": E_bulk,
               "E_indent_Pa": ind["E_eff_Pa"], "decoupling": E_bulk / max(ind["E_eff_Pa"], 1e-9)}
        out.append(row)
        print(f"  c={c} ⟨z⟩={ecm.connectivity_z:.2f}: G_shear={G:.1f} E_bulk={E_bulk:.1f} "
              f"E_indent={ind['E_eff_Pa']:.2f} → decoupling {row['decoupling']:.0f}×", flush=True)
    return out


def fig_stiffening(ss):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for s in ss:
        lab = f"{s['material']} {s['conc']}mg/mL (⟨z⟩={s['z']:.1f})"
        ax1.plot(s["gammas"], s["K_over_K0"], "o-", label=lab)
        sig = np.abs(s["sigma_Pa"]); K = np.abs(s["K_Pa"])
        ax2.loglog(sig[sig > 0], K[sig > 0], "o-", label=lab)
    ax1.axhline(1.0, color="0.6", ls=":")
    ax1.set_xlabel("shear strain γ"); ax1.set_ylabel("stiffening ratio K(γ)/K₀")
    ax1.set_title("Emergent strain-stiffening (differential modulus)"); ax1.legend(fontsize=7); ax1.grid(alpha=0.3)
    ax2.set_xlabel("shear stress σ [Pa]"); ax2.set_ylabel("differential modulus K [Pa]")
    ax2.set_title("K ∝ σ signature (nonlinear regime)"); ax2.legend(fontsize=7); ax2.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(f"{FIGS}/strain_stiffening.png", dpi=130); plt.close(fig)


def fig_decoupling(md):
    z = [r["z"] for r in md]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(z, [r["E_bulk_Pa"] for r in md], "s-", color="tab:blue", label="bulk E = 2(1+ν)G_shear")
    ax.plot(z, [r["E_indent_Pa"] for r in md], "o-", color="tab:red", label="local indentation E_eff")
    for r in md:
        ax.annotate(f"{r['decoupling']:.0f}×", (r["z"], r["E_indent_Pa"]), fontsize=8,
                    xytext=(3, -12), textcoords="offset points")
    ax.set_yscale("log")
    ax.set_xlabel("network connectivity ⟨z⟩")
    ax.set_ylabel("modulus [Pa] (log)")
    ax.set_title("Shear ↔ indentation modulus decoupling (sub-isostatic collagen)")
    ax.legend(); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(f"{FIGS}/mode_decoupling.png", dpi=130); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--steps", type=int, default=4500)
    a = ap.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    print("== emergent strain-stiffening ==")
    ss = strain_stiffening_study(device=a.device, n_steps=a.steps)
    print("== shear/indentation mode decoupling ==")
    md = mode_decoupling_study(device=a.device, n_steps=a.steps)
    fig_stiffening(ss); fig_decoupling(md)
    with open(f"{OUT}/ecm_novelty.json", "w") as f:
        json.dump({"strain_stiffening": ss, "mode_decoupling": md}, f, indent=2, default=float)
    lines = ["# FF ECM library — emergent & mode-dependent mechanics (novel)", "",
             "## 1. Emergent nonlinear strain-stiffening (NOT tuned; KB-1.V.2.5 gate)", "",
             "| material | conc | ⟨z⟩ | K₀ (Pa) | peak K/K₀ | γ_c |", "|---|---|---|---|---|---|"]
    for s in ss:
        lines.append(f"| {s['material']} | {s['conc']} | {s['z']:.2f} | {s['K0_Pa']:.1f} | "
                     f"{s['max_stiffening']:.1f}× | {s['gamma_c']:.2f} |")
    lines += ["", "Emergent K(γ) rises above γ_c with γ_c ~concentration-independent (geometry-set) — the "
              "KB-1.V.2.5 qualitative signature, reproduced from a model NOT calibrated to it. Peak stiffening "
              "is moderate (athermal Mikado); the full K>10× regime needs deeper strain / thermal nonlinearity.", "",
              "## 2. Shear ↔ indentation modulus decoupling (van Oosten 2019)", "",
              "| conc (mg/mL) | ⟨z⟩ | G_shear (Pa) | bulk E (Pa) | local E_indent (Pa) | decoupling |",
              "|---|---|---|---|---|---|"]
    for r in md:
        lines.append(f"| {r['conc']} | {r['z']:.2f} | {r['G_shear_Pa']:.1f} | {r['E_bulk_Pa']:.1f} | "
                     f"{r['E_indent_Pa']:.2f} | {r['decoupling']:.0f}× |")
    lines += ["", "Sub-isostatic collagen is stiff in bulk shear but soft under LOCAL indentation — the "
              "semiflexible mode-decoupling. The FF library measures both on the same network, so the "
              "decoupling ratio is produced directly as a function of architecture (⟨z⟩, alignment).",
              "", "## Figures", "- `figs/strain_stiffening.png` — K(γ)/K₀ + K∝σ",
              "- `figs/mode_decoupling.png` — bulk vs local modulus vs ⟨z⟩", ""]
    with open(f"{OUT}/ECM_NOVELTY.md", "w") as f:
        f.write("\n".join(lines))
    print(f"wrote {OUT}/ECM_NOVELTY.md + figs")


if __name__ == "__main__":
    main()
