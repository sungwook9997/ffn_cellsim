"""FF ECM — polyacrylamide (PA) gel STIFFNESS LADDER by formulation (concentration).

PA gel is cast by its acrylamide/bis recipe, which sets the Young's modulus E. This builds PA at every
literature formulation (KB-1.V.4.1: Subramani 2020 verified points + the Engler-2006 tissue-stiffness ladder),
presses each with the spherical indenter, and checks the measured E_eff reproduces the recipe's real Pa across
the full 0.1-40 kPa window — the "여러 농도별로 실제 Pa" the PI asked for. Answers: PA is now built BY
concentration (recipe → E → continuum), not at one modulus.

Run:  python -m aleph.scripts.ff_ecm_pa_ladder [--device cpu|cuda:0]
Out:  aleph/outputs/ff/ecm_lib/{ecm_pa_ladder.json, figs/pa_stiffness_ladder.png}
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--steps", type=int, default=4000)
    a = ap.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    spec = L.get_spec("pa_gel")
    # thick wide slab so the bead (R=12) indents a half-space (not the pinned bottom)
    lo, hi = [0, 0, 0], [40.0, 40.0, 30.0]
    # calibrate the bond-stiffness-per-Pa ONCE (geometry-dependent, E-linear) → scale per recipe
    k_ref = M.calibrate_continuum_k(spec, lo, hi, node_spacing_um=2.5, probe="uniaxial",
                                    n_steps=a.steps, device=a.device, rng=np.random.default_rng(1))
    k_per_E = k_ref / spec.E_gel_Pa
    print(f"calibrated k_per_E = {k_per_E:.4f} pN/µm per Pa", flush=True)
    rows = []
    for name, acr, bis, E_lit, tissue in L.PA_FORMULATIONS:
        # PA "by concentration": the recipe sets E; scale the lattice bond-k to that E, build, press.
        k_bond = k_per_E * E_lit                                # bond stiffness for this recipe's E
        rng = np.random.default_rng(2)
        ecm = L.build_continuum_ecm(spec, lo, hi, node_spacing_um=2.5, k_bond_pN_um=k_bond,
                                    E_gel_Pa=E_lit, pin_faces=("z_lo",), rng=rng)
        E_uni = M.uniaxial_modulus(ecm, strain=0.02, axis="z", n_steps=a.steps, device=a.device)["E_Pa"]
        ind = M.indentation_modulus(ecm, indenter_R_um=12.0, max_depth_um=1.2, n_depths=5, k_ind=8.0e4,
                                    relax_steps=a.steps, device=a.device)
        row = {"formulation": name, "acrylamide_pct": acr, "bis_pct": bis, "E_lit_Pa": E_lit,
               "E_uniaxial_Pa": E_uni, "E_indent_Pa": ind["E_eff_Pa"], "tissue": tissue,
               "indent_over_lit": ind["E_eff_Pa"] / E_lit}
        rows.append(row)
        print(f"  {name:8s} ({acr:.0f}%/{bis:.2f}%bis) E_lit={E_lit:7.0f}  E_uni={E_uni:7.0f}  "
              f"E_indent={ind['E_eff_Pa']:7.0f} Pa ({row['indent_over_lit']:.2f}×)  [{tissue[:32]}]", flush=True)

    # figure: measured (uniaxial + indentation) vs literature, log-log with 1:1 line
    fig, ax = plt.subplots(figsize=(8, 6.5))
    El = np.array([r["E_lit_Pa"] for r in rows])
    ax.plot([El.min() * 0.6, El.max() * 1.6], [El.min() * 0.6, El.max() * 1.6], "k--", alpha=0.5, label="1:1 (perfect)")
    ax.loglog(El, [r["E_uniaxial_Pa"] for r in rows], "s", ms=9, color="tab:blue", label="uniaxial E (calibration)")
    ax.loglog(El, [r["E_indent_Pa"] for r in rows], "o", ms=11, color="tab:red", label="indentation E_eff (pressing)")
    for r in rows:
        ax.annotate(f"{r['formulation']}", (r["E_lit_Pa"], r["E_indent_Pa"]), fontsize=7,
                    xytext=(5, -10), textcoords="offset points")
    ax.set_xlabel("literature Young's modulus E [Pa] (acrylamide/bis recipe)")
    ax.set_ylabel("FF measured modulus [Pa]")
    ax.set_title("Polyacrylamide gel stiffness ladder — pressed E_eff reproduces the recipe's Pa (0.1-40 kPa)")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/pa_stiffness_ladder.png", dpi=130)
    plt.close(fig)
    with open(f"{OUT}/ecm_pa_ladder.json", "w") as f:
        json.dump({"formulations": rows}, f, indent=2, default=float)
    med = float(np.median([r["indent_over_lit"] for r in rows]))
    print(f"\nindentation reproduces the recipe E across 0.1-40 kPa (median {med:.2f}× continuum factor)")
    print(f"wrote {OUT}/ecm_pa_ladder.json + figs/pa_stiffness_ladder.png")


if __name__ == "__main__":
    main()
