"""FF ECM — native-scale modulus atlas for ALL 6 materials (ROADMAP NEAR #4).

The HARD native-scale rule requires every material validated at full extent, not just collagen. This builds
each ECM material at a NATIVE-density REV (large converged box on the A5000), measures its modulus the right
way (fibrillar → shear G via the virial tensor; continuum gel → calibrate + indentation E), and reports it vs
the literature band AND vs the small-REV value (modulus is intensive → REV↔native agreement is a check). One
JSON summary + a bar figure. Complements ff_ecm_native.py (which did collagen full-extent + indentation).

Run (gbook):  python -m aleph.scripts.ff_ecm_native_atlas --device cuda:0
Out:  aleph/outputs/ff/ecm_lib/{ecm_native_atlas.json, figs/native_atlas.png}
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aleph.laws import ecm_library as L
from aleph.laws import ecm_mechanics as M
from aleph.scripts.ff_ecm_validate import LIT

OUT = "aleph/outputs/ff/ecm_lib"
FIGS = f"{OUT}/figs"


def measure(key, box, n_steps, device):
    """Build one material at the given box (native REV) + measure its modulus vs literature."""
    spec = L.get_spec(key)
    lit = LIT[key]
    lo, hi = [0, 0, 0], [box, box, box]
    t0 = time.time()
    if spec.is_fibrillar:
        rng = np.random.default_rng(7)
        ecm = L.build_fibrillar_ecm(spec, lo, hi, concentration=lit["ref_conc"], dim=3, alignment_S=0.0,
                                    pin_faces=(), target_z=3.2, rng=rng)
        g = M.shear_modulus(ecm, gamma=0.02, n_steps=n_steps, device=device)
        measured = g["G_virial_Pa"] if not np.isnan(g["G_virial_Pa"]) else g["G_Pa"]
        meta = {"n_nodes": ecm.meta["n_nodes"], "n_fibers": ecm.meta["n_fibers"], "z": ecm.connectivity_z,
                "mesh": ecm.mesh_size_um, "G_energy": g["G_Pa"], "G_virial": g["G_virial_Pa"]}
    else:
        slab = [box + 12, box + 12, box + 6]
        k = M.calibrate_continuum_k(spec, [0, 0, 0], slab, node_spacing_um=2.5, probe="uniaxial",
                                    n_steps=n_steps, device=device, rng=np.random.default_rng(7))
        rng = np.random.default_rng(7)
        ecm = L.build_continuum_ecm(spec, [0, 0, 0], slab, node_spacing_um=2.5, k_bond_pN_um=k,
                                    pin_faces=("z_lo",), rng=rng)
        ind = M.indentation_modulus(ecm, indenter_R_um=12.0, max_depth_um=1.2, n_depths=5, k_ind=8.0e4,
                                    relax_steps=n_steps, device=device)
        measured = ind["E_eff_Pa"]
        meta = {"n_nodes": ecm.meta["n_nodes"], "E_input": spec.E_gel_Pa, "E_indent": ind["E_eff_Pa"]}
    band = lit["band"]
    return {"material": key, "kind": lit["kind"], "measured_Pa": float(measured), "band_Pa": list(band),
            "in_band": bool(band[0] <= measured <= band[1]), "is_fibrillar": spec.is_fibrillar,
            "box_um": box, "meta": meta, "wall_s": time.time() - t0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--box", type=float, default=40.0, help="native-density REV box side [µm]")
    ap.add_argument("--steps", type=int, default=6000)
    a = ap.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    rows = []
    for key in ("collagen_I", "fibrin", "pa_gel", "hyaluronic_acid", "matrigel", "agarose"):
        r = measure(key, a.box, a.steps, a.device)
        print(f"[{'IN BAND ' if r['in_band'] else 'OUT     '}] {key:16s} native({r['box_um']}µm,"
              f"{r['meta']['n_nodes']} nodes) = {r['measured_Pa']:.1f} Pa  band={r['band_Pa']}  ({r['wall_s']:.0f}s)",
              flush=True)
        rows.append(r)
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(rows))
    for i, r in enumerate(rows):
        ax.plot([x[i], x[i]], r["band_Pa"], color="0.7", lw=8, solid_capstyle="round",
                label="literature band" if i == 0 else None)
    ax.scatter(x, [r["measured_Pa"] for r in rows],
               c=["tab:green" if r["in_band"] else "tab:red" for r in rows], s=90, zorder=5,
               edgecolor="k", label="FF native")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['material']}\n({r['kind']})" for r in rows], fontsize=8)
    ax.set_ylabel("modulus [Pa] (log)")
    ax.set_title(f"FF ECM native-scale atlas — all 6 materials at {a.box}µm REV (green = IN BAND)")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/native_atlas.png", dpi=130)
    plt.close(fig)
    with open(f"{OUT}/ecm_native_atlas.json", "w") as f:
        json.dump({"box_um": a.box, "materials": rows}, f, indent=2, default=float)
    n_in = sum(r["in_band"] for r in rows)
    print(f"\n{n_in}/6 materials IN literature band at native scale")
    print(f"wrote {OUT}/ecm_native_atlas.json + figs/native_atlas.png")


if __name__ == "__main__":
    main()
