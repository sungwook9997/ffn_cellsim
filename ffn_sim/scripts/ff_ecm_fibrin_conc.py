"""FF ECM — fibrin gel modulus vs fibrinogen CONCENTRATION (the fibrillar analog of collagen G'(c)).

Following the PI's per-concentration emphasis, this sweeps fibrin across its fibrinogen concentration and
measures the emergent shear modulus, comparing to the literature (KB-1.V.4.2: Piechocka 2010, G0 0.1-2000 Pa
over 0.1-8 mg/mL, power law G~c^2.3). Fibrin was previously validated at a single concentration; this makes
the fibrillar per-concentration validation symmetric with collagen. Modulus EMERGES from the Mikado network
(no tuning); ⟨z⟩ held at the physical value via target_z.

Run:  python -m ffn_sim.scripts.ff_ecm_fibrin_conc [--device cpu|cuda:0]
Out:  ffn_sim/outputs/ff/ecm_lib/{ecm_fibrin_conc.json, figs/fibrin_concentration.png}
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--steps", type=int, default=6000)
    a = ap.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    spec = L.get_spec("fibrin")
    concs = [0.5, 1.0, 2.0, 4.0, 8.0]                          # mg/mL fibrinogen
    # Piechocka 2010 (KB-1.V.4.2): G0 rises 0.1-2000 Pa over ~0.1-8 mg/mL, power law ~c^2.3. Anchor to the
    # model's own 2 mg/mL value and draw the literature c^2.3 trend through it for the shape comparison.
    rows = []
    for c in concs:
        rng = np.random.default_rng(5)
        ecm = L.build_fibrillar_ecm(spec, [0, 0, 0], [26, 26, 26], concentration=c, dim=3,
                                    alignment_S=0.0, pin_faces=(), target_z=3.4, rng=rng)
        G = M.shear_modulus(ecm, gamma=0.02, n_steps=a.steps, device=a.device)["G_Pa"]
        rows.append({"conc": c, "G_sim": G, "z": ecm.connectivity_z, "mesh": ecm.mesh_size_um,
                     "nfib": ecm.meta["n_fibers"]})
        print(f"  fibrin {c:4.1f} mg/mL: ⟨z⟩={ecm.connectivity_z:.2f} mesh={ecm.mesh_size_um:.2f} "
              f"G={G:.1f} Pa", flush=True)
    cc = np.array(concs); Gs = np.array([r["G_sim"] for r in rows])
    n = float(np.polyfit(np.log(cc), np.log(Gs), 1)[0])
    # literature c^2.3 trend anchored at the model's 2 mg/mL point
    G2 = Gs[concs.index(2.0)]
    G_lit = G2 * (cc / 2.0) ** 2.3
    band_ok = all(0.1 <= g <= 2000.0 for g in Gs)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.axhspan(0.1, 2000.0, alpha=0.10, color="green", label="Piechocka band 0.1-2000 Pa")
    ax.loglog(cc, G_lit, "s--", color="k", label="literature trend G~c^2.3 (Piechocka)")
    ax.loglog(cc, Gs, "o-", color="tab:red", label=f"FF fibrin (emergent, n={n:.2f})")
    ax.set_xlabel("fibrinogen concentration [mg/mL]")
    ax.set_ylabel("shear modulus G' [Pa]")
    ax.set_title("Fibrin G'(c): FF emergent modulus vs literature (Piechocka 2010)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/fibrin_concentration.png", dpi=130)
    plt.close(fig)
    with open(f"{OUT}/ecm_fibrin_conc.json", "w") as f:
        json.dump({"rows": rows, "exponent_sim": n, "exponent_lit": 2.3, "in_band": band_ok}, f, indent=2, default=float)
    print(f"\nfibrin exponent n={n:.2f} (lit ~2.3); all G in Piechocka band 0.1-2000 Pa: {band_ok}")
    print(f"wrote {OUT}/ecm_fibrin_conc.json + figs/fibrin_concentration.png")


if __name__ == "__main__":
    main()
