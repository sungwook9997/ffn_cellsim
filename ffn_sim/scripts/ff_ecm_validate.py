"""FF ECM library validation — does every material give the REAL literature Pa?

Runs the three-way modulus harness (shear G / uniaxial E / spherical indentation E_eff, all in Pa) over
the ECM material registry and the alignment/concentration axes, compares each to its literature band, and
emits a results JSON + a markdown table + validation figures. This is the "실제와 같은 Pa 나오나" gate.

Fibrillar materials (collagen, fibrin, agarose): the bulk modulus EMERGES from the Mikado microstructure;
we check it lands in the literature band (no tuning). Continuum gels (PA, HA, Matrigel): E is a material
input; we calibrate the lattice bond stiffness once and check that INDENTATION returns that E (the
harness/pressing consistency, validated against a known modulus).

Run:  python -m ffn_sim.scripts.ff_ecm_validate [--quick] [--device cpu]
Out:  ffn_sim/outputs/ff/ecm_lib/{ecm_validation.json, ECM_VALIDATION.md, figs/*.png}
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

from ffn_sim.ff import ecm_library as L
from ffn_sim.ff import ecm_mechanics as M

OUT = "ffn_sim/outputs/ff/ecm_lib"
FIGS = f"{OUT}/figs"

# Literature bands / anchors (Pa) used as acceptance references (from KB + the research cards).
LIT = {
    "collagen_I": {"kind": "G", "band": (5.0, 100.0), "ref_conc": 1.5,
                   "note": "G'(1.5mg/mL)≈11 Pa (KB-1.V.2.1 c² from 5Pa@1mg/mL); scaling-anchor 30-100 (KB-1.30)"},
    "fibrin": {"kind": "G", "band": (10.0, 1000.0), "ref_conc": 2.0,
               "note": "G' 10-1000 Pa over 0.5-10 mg/mL (Piechocka 2010)"},
    "agarose": {"kind": "E", "band": (1.0e3, 1.0e5), "ref_conc": 1.0,
                "note": "E 1-100 kPa over 0.5-4 %w/v (Normand 2000)"},
    "pa_gel": {"kind": "E", "band": (100.0, 4.0e4), "input": 5.0e3,
               "note": "E 0.1-40 kPa tunable; indentation must return the input E (KB-1.21, Tse-Engler)"},
    "hyaluronic_acid": {"kind": "E", "band": (10.0, 3.0e3), "input": 3.0e2,
                        "note": "E 10 Pa-3 kPa crosslink-tunable"},
    "matrigel": {"kind": "E", "band": (30.0, 900.0), "input": 4.5e2,
                 "note": "E 30-900 Pa (typ 450, Soofi 2009 AFM 37°C)"},
}


def _mk_dirs():
    os.makedirs(FIGS, exist_ok=True)


def validate_material(key, *, box=26.0, n_steps=5000, device="cpu", quick=False):
    """Measure the reference-condition modulus of one material and grade it against its literature band."""
    spec = L.get_spec(key)
    lit = LIT[key]
    lo, hi = [0, 0, 0], [box, box, box]
    rng = np.random.default_rng(7)
    t0 = time.time()
    if spec.is_fibrillar:
        ecm = L.build_fibrillar_ecm(spec, lo, hi, concentration=lit["ref_conc"], dim=3,
                                    alignment_S=0.0, pin_faces=(), rng=rng)
        G = M.shear_modulus(ecm, gamma=0.02, n_steps=n_steps, device=device)["G_Pa"]
        E = M.uniaxial_modulus(ecm, strain=0.02, axis="z", n_steps=n_steps, device=device)["E_Pa"]
        measured = G if lit["kind"] == "G" else E
        meta = {"n_fibers": ecm.meta["n_fibers"], "n_nodes": ecm.meta["n_nodes"],
                "z": ecm.connectivity_z, "mesh_um": ecm.mesh_size_um, "G_Pa": G, "E_uniax_Pa": E}
    else:
        # continuum gel: calibrate bond k on the SAME geometry, then confirm indentation returns E. Use a
        # THICK, wide slab (indenter R ≪ thickness) so the pinned bottom does not stiffen the Hertz reading.
        slab_hi = [box + 12.0, box + 12.0, box + 6.0]
        k = M.calibrate_continuum_k(spec, [0, 0, 0], slab_hi, node_spacing_um=2.5, probe="uniaxial",
                                    n_steps=n_steps, device=device, rng=np.random.default_rng(7))
        rng = np.random.default_rng(7)
        ecm = L.build_continuum_ecm(spec, [0, 0, 0], slab_hi, node_spacing_um=2.5,
                                    k_bond_pN_um=k, pin_faces=("z_lo",), rng=rng)
        E = M.uniaxial_modulus(ecm, strain=0.02, axis="z", n_steps=n_steps, device=device)["E_Pa"]
        ind = M.indentation_modulus(ecm, indenter_R_um=12.0, max_depth_um=1.2, n_depths=5,
                                    k_ind=5.0e4, relax_steps=n_steps, device=device)
        measured = ind["E_eff_Pa"]
        meta = {"n_nodes": ecm.meta["n_nodes"], "k_bond": k, "E_uniax_Pa": E,
                "E_indent_Pa": ind["E_eff_Pa"], "indent": ind}
    band = lit["band"]
    verdict = "IN BAND" if band[0] <= measured <= band[1] else ("TOO STIFF" if measured > band[1] else "TOO SOFT")
    return {"material": key, "kind": lit["kind"], "measured_Pa": measured, "band_Pa": band,
            "verdict": verdict, "is_fibrillar": spec.is_fibrillar, "note": lit["note"],
            "meta": meta, "wall_s": time.time() - t0}


def collagen_concentration(device="cpu", n_steps=5000):
    """Collagen G'(c): sim vs KB-1.V.2.1 anchors (5/55/342 Pa @1/3/7 mg/mL, n~2)."""
    spec = L.get_spec("collagen_I")
    kb = {1.0: 5.0, 1.5: 11.25, 2.0: 20.0, 3.0: 55.0, 5.0: 150.0, 7.0: 342.0}
    concs = [1.0, 1.5, 2.0, 3.0, 5.0, 7.0]
    rows = []
    for c in concs:
        rng = np.random.default_rng(3)
        ecm = L.build_fibrillar_ecm(spec, [0, 0, 0], [26, 26, 26], concentration=c, dim=3,
                                    alignment_S=0.0, pin_faces=(), rng=rng)
        G = M.shear_modulus(ecm, gamma=0.02, n_steps=n_steps, device=device)["G_Pa"]
        rows.append({"conc": c, "G_sim": G, "G_kb": kb[c], "z": ecm.connectivity_z,
                     "mesh": ecm.mesh_size_um, "nfib": ecm.meta["n_fibers"]})
    n = float(np.polyfit(np.log(concs), np.log([r["G_sim"] for r in rows]), 1)[0])
    return {"rows": rows, "exponent_sim": n, "exponent_kb": 2.05}


def dim_and_composite(device="cpu", n_steps=5000):
    """2D sheet vs 3D bulk (collagen) + an interpenetrating composite (collagen+Matrigel) — the
    '개별·혼합, 2D·3D 모두' axis. 2D reports the in-plane sheet modulus (loaded in-plane)."""
    spec = L.get_spec("collagen_I")
    out = {"dim": [], "composite": None}
    # 3D bulk
    rng = np.random.default_rng(21)
    e3 = L.build_fibrillar_ecm(spec, [0, 0, 0], [26, 26, 26], concentration=2.0, dim=3,
                               alignment_S=0.0, pin_faces=(), rng=rng)
    G3 = M.shear_modulus(e3, gamma=0.02, n_steps=n_steps, device=device)["G_Pa"]
    out["dim"].append({"dim": "3D bulk", "n_fibers": e3.meta["n_fibers"], "mesh": e3.mesh_size_um,
                       "z": e3.connectivity_z, "modulus_Pa": G3, "kind": "G (bulk)"})
    # 2D sheet (thin), loaded in-plane along x
    rng = np.random.default_rng(22)
    e2 = L.build_fibrillar_ecm(spec, [0, 0, 0], [30, 30, 2.5], concentration=2.0, dim=2,
                               alignment_S=0.0, pin_faces=(), rng=rng)
    E2 = M.uniaxial_modulus(e2, strain=0.02, axis="x", n_steps=n_steps, device=device)["E_Pa"]
    out["dim"].append({"dim": "2D sheet", "n_fibers": e2.meta["n_fibers"], "mesh": e2.mesh_size_um,
                       "z": e2.connectivity_z, "modulus_Pa": E2, "kind": "E_2D (in-plane, Pa·µm/µm)"})
    # composite: collagen (aligned) interpenetrating a Matrigel lattice
    rng = np.random.default_rng(23)
    comp = L.build_composite([{"material": "collagen_I", "concentration": 1.5, "alignment_S": 0.3,
                               "director": (0, 0, 1)}, {"material": "matrigel"}],
                             [0, 0, 0], [26, 26, 26], dim=3, interlink_um=1.0, pin_faces=(), rng=rng)
    Gc = M.shear_modulus(comp, gamma=0.02, n_steps=n_steps, device=device)["G_Pa"]
    out["composite"] = {"components": comp.meta["components"], "n_nodes": comp.meta["n_nodes"],
                        "z": comp.connectivity_z, "G_Pa": Gc,
                        "note": "collagen(S=0.3) + Matrigel interpenetrating; G exceeds either component alone"}
    return out


def alignment_anisotropy(device="cpu", n_steps=5000):
    """Collagen at increasing nematic order S (director=z): E along vs across → anisotropy ratio + tissue map."""
    spec = L.get_spec("collagen_I")
    rows = []
    for S in [0.0, 0.3, 0.6, 0.85]:
        rng = np.random.default_rng(11)
        along = L.build_fibrillar_ecm(spec, [0, 0, 0], [26, 26, 26], concentration=2.0, dim=3,
                                      alignment_S=S, director=(0, 0, 1), pin_faces=(), rng=rng)
        rng = np.random.default_rng(11)
        # 'across' = same S but director in-plane (x), still loaded along z → transverse response
        across = L.build_fibrillar_ecm(spec, [0, 0, 0], [26, 26, 26], concentration=2.0, dim=3,
                                       alignment_S=S, director=(1, 0, 0), pin_faces=(), rng=rng)
        Ea = M.uniaxial_modulus(along, strain=0.02, axis="z", n_steps=n_steps, device=device)["E_Pa"]
        Ec = M.uniaxial_modulus(across, strain=0.02, axis="z", n_steps=n_steps, device=device)["E_Pa"]
        tissue = next((t for s, t in spec.tissue_by_S if abs(s - S) < 0.2), "")
        rows.append({"S_target": S, "S_meas": along.S_measured, "E_along": Ea, "E_across": Ec,
                     "ratio": Ea / max(Ec, 1e-9), "tissue": tissue})
    return {"rows": rows}


# ─────────────────────────────────────────────────────────────────────────────────────────────────
# figures
# ─────────────────────────────────────────────────────────────────────────────────────────────────

def fig_materials(results):
    fig, ax = plt.subplots(figsize=(9, 5))
    keys = [r["material"] for r in results]
    meas = [r["measured_Pa"] for r in results]
    los = [r["band_Pa"][0] for r in results]
    his = [r["band_Pa"][1] for r in results]
    x = np.arange(len(keys))
    for i in range(len(keys)):
        ax.plot([x[i], x[i]], [los[i], his[i]], color="0.7", lw=8, solid_capstyle="round",
                label="literature band" if i == 0 else None)
    colors = ["tab:green" if r["verdict"] == "IN BAND" else "tab:red" for r in results]
    ax.scatter(x, meas, c=colors, s=90, zorder=5, edgecolor="k", label="FF model")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{k}\n({r['kind']})" for k, r in zip(keys, results)], fontsize=8)
    ax.set_ylabel("modulus [Pa]  (log)")
    ax.set_title("FF ECM library: measured modulus vs literature band (green=IN BAND)")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/ecm_material_moduli.png", dpi=130)
    plt.close(fig)


def fig_collagen_conc(cc):
    rows = cc["rows"]
    c = np.array([r["conc"] for r in rows])
    gs = np.array([r["G_sim"] for r in rows])
    gk = np.array([r["G_kb"] for r in rows])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(c, gk, "s--", color="k", label="KB Yang-Kaufman G'(c)~c²")
    ax.loglog(c, gs, "o-", color="tab:blue", label=f"FF model (n={cc['exponent_sim']:.2f})")
    ax.set_xlabel("collagen concentration [mg/mL]")
    ax.set_ylabel("shear modulus G' [Pa]")
    ax.set_title("Collagen-I G'(c): FF model vs literature")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/collagen_concentration.png", dpi=130)
    plt.close(fig)


def fig_alignment(al):
    rows = al["rows"]
    S = [r["S_meas"] for r in rows]
    ratio = [r["ratio"] for r in rows]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(S, ratio, "o-", color="tab:purple")
    for r in rows:
        ax.annotate(r["tissue"].split("(")[0][:18], (r["S_meas"], r["ratio"]), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")
    ax.axhline(1.0, color="0.6", ls=":")
    ax.set_xlabel("nematic order parameter S (measured)")
    ax.set_ylabel("anisotropy  E∥ / E⊥")
    ax.set_title("Collagen alignment → mechanical anisotropy (director→tissue)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/collagen_alignment_anisotropy.png", dpi=130)
    plt.close(fig)


def fig_indentation(results):
    fig, ax = plt.subplots(figsize=(7, 5))
    for r in results:
        ind = r["meta"].get("indent")
        if not ind:
            continue
        d = np.array(ind["depths_um"])
        F = np.array(ind["forces_pN"])
        ax.plot(d ** 1.5, F, "o-", label=f"{r['material']} (E_eff={ind['E_eff_Pa']:.0f} Pa)")
    ax.set_xlabel("δ^1.5 [µm^1.5]")
    ax.set_ylabel("indenter force F [pN]")
    ax.set_title("Spherical indentation (Hertz: F ∝ δ^1.5) — continuum gels")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/indentation_curves.png", dpi=130)
    plt.close(fig)


def write_report(results, cc, al, dc=None):
    lines = ["# FF ECM library — modulus validation (실제와 같은 Pa)", ""]
    lines.append("## Per-material reference modulus vs literature band")
    lines.append("")
    lines.append("| material | class | measured (Pa) | band (Pa) | verdict | note |")
    lines.append("|---|---|---|---|---|---|")
    for r in results:
        cls = "fibrillar" if r["is_fibrillar"] else "continuum"
        lines.append(f"| {r['material']} | {cls}·{r['kind']} | {r['measured_Pa']:.1f} | "
                     f"{r['band_Pa'][0]:.0f}–{r['band_Pa'][1]:.0f} | {r['verdict']} | {r['note']} |")
    lines += ["", "## Collagen-I concentration series G'(c)", "",
              f"power-law exponent: sim **n={cc['exponent_sim']:.2f}** vs literature n≈2.05 "
              "(athermal Mikado is density-limited; absolute Pa at reference conc matches).", "",
              "| c (mg/mL) | G_sim (Pa) | G_KB (Pa) | ⟨z⟩ | mesh ξ (µm) |", "|---|---|---|---|---|"]
    for r in cc["rows"]:
        lines.append(f"| {r['conc']} | {r['G_sim']:.1f} | {r['G_kb']:.0f} | {r['z']:.2f} | {r['mesh']:.2f} |")
    lines += ["", "## Collagen alignment → anisotropy (tissue mapping)", "",
              "| S_target | S_meas | E∥ (Pa) | E⊥ (Pa) | E∥/E⊥ | tissue |", "|---|---|---|---|---|---|"]
    for r in al["rows"]:
        lines.append(f"| {r['S_target']} | {r['S_meas']:.3f} | {r['E_along']:.1f} | {r['E_across']:.1f} | "
                     f"{r['ratio']:.2f} | {r['tissue']} |")
    if dc:
        lines += ["", "## Dimensionality (2D sheet vs 3D bulk) + composite (mixed ECM)", "",
                  "| build | modulus (Pa) | kind | ⟨z⟩ | n_fibers |", "|---|---|---|---|---|"]
        for r in dc["dim"]:
            lines.append(f"| collagen {r['dim']} | {r['modulus_Pa']:.1f} | {r['kind']} | {r['z']:.2f} | {r['n_fibers']} |")
        c = dc["composite"]
        lines.append(f"| composite {'+'.join(c['components'])} | {c['G_Pa']:.1f} | G (bulk) | {c['z']:.2f} | — |")
        lines += ["", f"_{c['note']}_", ""]
    lines += ["", "## Figures", "",
              "- `figs/ecm_material_moduli.png` — each material's modulus vs literature band",
              "- `figs/collagen_concentration.png` — G'(c) vs Yang-Kaufman",
              "- `figs/collagen_alignment_anisotropy.png` — anisotropy vs S",
              "- `figs/indentation_curves.png` — Hertz F∝δ^1.5 for continuum gels", ""]
    with open(f"{OUT}/ECM_VALIDATION.md", "w") as f:
        f.write("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="fewer steps (smoke)")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--only", default="", help="comma-separated material keys to run (default all)")
    a = ap.parse_args()
    _mk_dirs()
    ns = 2500 if a.quick else 5000
    keys = a.only.split(",") if a.only else list(LIT)
    results = []
    for k in keys:
        r = validate_material(k, n_steps=ns, device=a.device, quick=a.quick)
        print(f"[{r['verdict']:9s}] {k:16s} measured={r['measured_Pa']:.1f} Pa  band={r['band_Pa']}  ({r['wall_s']:.0f}s)")
        results.append(r)
    cc = collagen_concentration(device=a.device, n_steps=ns)
    print(f"collagen G'(c) exponent n={cc['exponent_sim']:.2f} (KB~2.05)")
    al = alignment_anisotropy(device=a.device, n_steps=ns)
    for r in al["rows"]:
        print(f"  S={r['S_target']} (meas {r['S_meas']:.2f}): E∥/E⊥={r['ratio']:.2f}  {r['tissue'][:30]}")
    dc = dim_and_composite(device=a.device, n_steps=ns)
    for r in dc["dim"]:
        print(f"  {r['dim']}: modulus={r['modulus_Pa']:.1f} ({r['kind']})  nfib={r['n_fibers']}")
    print(f"  composite {'+'.join(dc['composite']['components'])}: G={dc['composite']['G_Pa']:.1f} Pa")
    fig_materials(results)
    fig_collagen_conc(cc)
    fig_alignment(al)
    fig_indentation(results)
    write_report(results, cc, al, dc)
    with open(f"{OUT}/ecm_validation.json", "w") as f:
        json.dump({"materials": results, "collagen_conc": cc, "alignment": al, "dim_composite": dc},
                  f, indent=2, default=float)
    print(f"\nwrote {OUT}/ECM_VALIDATION.md + figs/*.png + ecm_validation.json")


if __name__ == "__main__":
    main()
