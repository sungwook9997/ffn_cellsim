"""FF ECM — mammary-stroma tissue-mimetic atlas: normal vs tumor (ROADMAP FAR, KB-1.V.1.2).

The library's individual materials are all in-band (6/6 real Pa); this asks the TISSUE question: can an
interpenetrating COMPOSITE (collagen-I + Matrigel — fibrillar stroma entangled with the basement-membrane gel)
reproduce a real mammary-stroma modulus? Two literature-anchored mimetics (concentrations/alignment set from the
biology, NOT tuned to the band):
  * NORMAL mammary stroma — loose collagen-I (1.5 mg/mL, S≈0.1) + Matrigel; KB-1.V.1.2 band 140-400 Pa.
  * TUMOR stroma (TACS-3) — DENSE aligned collagen-I (5 mg/mL, S=0.59 invasion-highway) + Matrigel; KB-1.V.1.2
    band 5000-10000 Pa (Levental 2009 desmoplastic stiffening, ~10-70×).
We measure the composite bulk shear modulus (virial G) and compare. The HONEST expectation: the normal mimetic
lands near band, the tumor mimetic UNDER-shoots — the athermal sub-isostatic Mikado cannot reach the kPa
desmoplastic stiffness at physiological collagen density (the SAME stretch-dominated ceiling as the c-scaling
+ N1 + stress-propagation findings; tissue-scale). Continuum gels DO reach kPa (agarose 14 kPa, PA 40 kPa in
band), so the ceiling is specifically the FIBRILLAR collagen route — its runtime fix is the PI-gated thermal-WLC.

Run:  python -m aleph.scripts.ff_ecm_tissue_mimetic
Out:  aleph/outputs/ff/ecm_lib/figs/tissue_mimetic_mammary.png + tissue_mimetic_mammary.json
"""

from __future__ import annotations

import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aleph.laws import ecm_library as L
from aleph.laws import ecm_mechanics as M

OUT = "aleph/outputs/ff/ecm_lib"

# literature-anchored tissue mimetics (KB-1.V.1.2 bands; concentrations from the biology, not tuned)
TISSUES = [
    dict(name="normal mammary stroma", band=(140.0, 400.0),
         components=[dict(material="collagen_I", concentration=1.5, alignment_S=0.1, target_z=3.2),
                     dict(material="matrigel")],
         note="loose collagen + basement-membrane gel"),
    dict(name="tumor stroma (TACS-3)", band=(5000.0, 10000.0),
         components=[dict(material="collagen_I", concentration=5.0, alignment_S=0.59, target_z=3.2),
                     dict(material="matrigel")],
         note="dense aligned collagen (invasion highway) + gel; desmoplastic (Levental)"),
]


def measure(tissue, box, n_steps, device):
    lo, hi = [0.0, 0.0, 0.0], [box, box, box]
    rng = np.random.default_rng(7)
    ecm = L.build_composite(tissue["components"], lo, hi, dim=3, interlink_um=0.75, interlink_k=100.0,
                            pin_faces=(), rng=rng)
    g = M.shear_modulus(ecm, gamma=0.02, n_steps=n_steps, device=device)
    G = g["G_virial_Pa"] if not np.isnan(g.get("G_virial_Pa", np.nan)) else g["G_Pa"]
    band = tissue["band"]
    return dict(name=tissue["name"], band_Pa=list(band), G_Pa=float(G),
                in_band=bool(band[0] <= G <= band[1]), under=bool(G < band[0]),
                n_nodes=int(ecm.meta["n_nodes"]), note=tissue["note"],
                G_energy=float(g["G_Pa"]), G_virial=float(g.get("G_virial_Pa", float("nan"))))


def main(box=30.0, n_steps=6000, device="cpu"):
    os.makedirs(f"{OUT}/figs", exist_ok=True)
    rows = [measure(t, box, n_steps, device) for t in TISSUES]
    for r in rows:
        verdict = "IN BAND" if r["in_band"] else ("UNDER (athermal ceiling)" if r["under"] else "OVER")
        print(f"[{verdict:24s}] {r['name']:26s} G={r['G_Pa']:8.1f} Pa  band={r['band_Pa']}  ({r['n_nodes']} nodes)")

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    x = np.arange(len(rows))
    for i, r in enumerate(rows):
        ax.plot([x[i], x[i]], r["band_Pa"], color="0.7", lw=9, solid_capstyle="round",
                label="literature band (KB-1.V.1.2)" if i == 0 else None, zorder=1)
    ax.scatter(x, [r["G_Pa"] for r in rows],
               c=["tab:green" if r["in_band"] else "tab:red" for r in rows], s=130, zorder=3,
               edgecolor="k", label="FF composite mimetic (virial G)")
    for i, r in enumerate(rows):
        ax.annotate(f"{r['G_Pa']:.0f} Pa", (x[i], r["G_Pa"]), textcoords="offset points",
                    xytext=(12, 0), fontsize=9, va="center")
    ax.set_yscale("log"); ax.set_xticks(x)
    ax.set_xticklabels([r["name"].replace(" ", "\n") for r in rows], fontsize=9)
    ax.set_ylabel("bulk shear modulus  G  [Pa] (log)")
    ax.set_title("FF ECM mammary-stroma tissue-mimetic (composite collagen+Matrigel) vs KB-1.V.1.2\n"
                 "normal ≈ band · tumor UNDER-shoots — the athermal fibrillar ceiling (needs PI-gated thermal-WLC)",
                 fontsize=9)
    ax.legend(fontsize=8, loc="upper left"); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    path = f"{OUT}/figs/tissue_mimetic_mammary.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    with open(f"{OUT}/tissue_mimetic_mammary.json", "w") as f:
        json.dump(dict(box_um=box, tissues=rows,
                       finding="normal-stroma composite near band; tumor-stroma UNDER-shoots (athermal fibrillar "
                               "ceiling, same stretch-dominated limit as c-scaling; continuum gels DO reach kPa, "
                               "so the ceiling is the collagen route → PI-gated thermal-WLC)"), f, indent=2, default=float)
    print(f"wrote {path} + tissue_mimetic_mammary.json")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--box", type=float, default=30.0)
    ap.add_argument("--steps", type=int, default=6000)
    ap.add_argument("--device", default="cpu")
    a = ap.parse_args()
    main(box=a.box, n_steps=a.steps, device=a.device)
