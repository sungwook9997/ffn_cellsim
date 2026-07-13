"""FF ECM library — consolidated validation dashboard (all native results vs their literature/KB anchors).

ONE reviewable figure tying the whole ECM-library program to its literature bands + KB claims (per the PI viz
directive: "각각 실제 파라미터랑 비교 피규어"). Reads the committed native-result JSONs (no re-simulation) and lays
out four panels:
  1. 6-material effective moduli vs literature bands (native atlas)                — real-Pa validation
  2. alignment → matrix directional response: library E∥/E⊥ + contact-guidance R_σ(S) (native)
  3. stress-propagation exponent n(S): collagen vs fibrin vs KB-1.10 (n~1 fibrous / n=3 continuum)
  4. stiffness sensing: Path b Winkler (ascending) + Path a live-ECM (descending) bracket the biphasic
Everything native-confirmed, unfitted; literature bands are acceptance oracles (overlaid, never fit targets).

Run:  python -m ffn_sim.scripts.ff_ecm_summary_dashboard
Out:  ffn_sim/outputs/ff/ecm_lib/figs/ecm_validation_dashboard.png
"""

from __future__ import annotations

import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "ffn_sim/outputs/ff/ecm_lib"
FIGS = f"{OUT}/figs"
E_ANISO_LIB = {0.0: 1.0, 0.02: 1.0, 0.30: 3.1, 0.59: 10.3, 0.83: 63.0}


def _load(name):
    p = f"{OUT}/{name}"
    return json.load(open(p)) if os.path.exists(p) else None


def _agg_by_S(rows, key):
    by = {}
    for r in rows:
        s = round(float(r.get("S_target", r.get("S", 0.0))), 2)
        v = r.get(key)
        if v is not None and np.isfinite(v):
            by.setdefault(s, []).append(float(v))
    S = sorted(by)
    return np.array(S), np.array([np.mean(by[s]) for s in S]), np.array([np.std(by[s]) for s in S])


def panel_moduli(ax):
    d = _load("ecm_native_atlas.json")
    if not d:
        ax.set_visible(False); return
    rows = d["materials"]
    x = np.arange(len(rows))
    for i, r in enumerate(rows):
        ax.plot([x[i], x[i]], r["band_Pa"], color="0.75", lw=7, solid_capstyle="round",
                label="literature band" if i == 0 else None, zorder=1)
    ax.scatter(x, [r["measured_Pa"] for r in rows],
               c=["tab:green" if r["in_band"] else "tab:red" for r in rows], s=70, zorder=3,
               edgecolor="k", label="FF native")
    ax.set_yscale("log"); ax.set_xticks(x)
    ax.set_xticklabels([r["material"].replace("_", "-") for r in rows], fontsize=7, rotation=25, ha="right")
    ax.set_ylabel("effective modulus [Pa]")
    n_in = sum(r["in_band"] for r in rows)
    ax.set_title(f"(1) 6 materials vs literature bands — {n_in}/{len(rows)} IN BAND (real Pa)", fontsize=9)
    ax.legend(fontsize=7, loc="upper left"); ax.grid(True, which="both", alpha=0.25)


def panel_alignment(ax):
    d = _load("contact_guidance_cg_native.json")
    Slad = sorted(E_ANISO_LIB)
    ax.plot(Slad, [E_ANISO_LIB[s] for s in Slad], "--D", color="0.5", ms=6, lw=1.6,
            label="library E∥/E⊥ (validated)")
    if d:
        S, rs, rsd = _agg_by_S(d["rows"], "R_sigma")
        ax.errorbar(S, rs, yerr=rsd, fmt="-o", color="#d62728", lw=2, ms=7, capsize=3,
                    label="contact-guidance R_σ=σ∥/σ⊥ (native)")
    ax.axhline(35, color="tab:purple", ls="-.", lw=1.2, label="Szulczewski ≤35× directional stiffness")
    ax.set_yscale("log"); ax.set_xlabel("nematic alignment order S")
    ax.set_ylabel("directional ratio")
    ax.set_title("(2) alignment → matrix directional response (R_σ EMERGES with S)", fontsize=9)
    ax.legend(fontsize=7, loc="upper left"); ax.grid(True, which="both", alpha=0.25)


def panel_stress_prop(ax):
    refs = [("stress_propagation_sp_native.json", "collagen-I", "#1f77b4"),
            ("stress_propagation_sp_fibrin_native.json", "fibrin", "#2ca02c")]
    for fn, lbl, col in refs:
        d = _load(fn)
        if not d:
            continue
        ne = d["meta"]["n_exp"]
        S = sorted(float(k.split("=")[1]) for k in ne)
        m = [ne[f"S={s:g}"]["mean"] for s in S]
        sd = [ne[f"S={s:g}"]["sd"] for s in S]
        ax.errorbar(S, m, yerr=sd, fmt="-o", color=col, lw=2, ms=7, capsize=3, label=f"{lbl} n(S)")
    ax.axhline(3.0, color="0.5", ls="--", lw=1.4, label="KB-1.10 continuum 3D (n=3)")
    ax.axhline(1.0, color="0.5", ls=":", lw=1.4, label="KB-1.10 fibrous long-range (n~1)")
    ax.set_xlabel("nematic alignment order S"); ax.set_ylabel("stress-decay exponent n  (|σ|~r⁻ⁿ)")
    ax.set_title("(3) stress propagation: n DROPS with alignment (aligned channels farther)", fontsize=9)
    ax.legend(fontsize=7, loc="upper right"); ax.grid(True, alpha=0.25)


def panel_sensing(ax):
    da = _load("ecm_stiffness_ecmnet_native.json")     # Path a (live ECM, descending)
    db = _load("ecm_stiffness_sensing_native.json")    # Path b (Winkler, ascending)
    if db:
        E = np.array([r["E_pa"] for r in db["rows"]]); tr = np.array([r["traction_nN"] for r in db["rows"]])
        ax.plot(E, tr, "-o", color="#1f77b4", lw=2, ms=6, label="Path b Winkler (ascending/catch)")
    if da:
        E = np.array([r["E_pa"] for r in da["rows"]]); tr = np.array([r["traction_nN"] for r in da["rows"]])
        ax.plot(E, tr, "-s", color="#d62728", lw=2, ms=6, label="Path a live ECM (descending/slip)")
    ax.set_xscale("log"); ax.set_xlabel("substrate / ECM modulus E [Pa]")
    ax.set_ylabel("engaged-clutch traction [nN]")
    ax.set_title("(4) stiffness sensing: two paths bracket the motor-clutch biphasic", fontsize=9)
    ax.legend(fontsize=7, loc="upper right"); ax.grid(True, which="both", alpha=0.25)


def panel_wlc(ax):
    d = _load("wlc_cscaling.json")
    if not d:
        ax.set_visible(False); return
    c = np.array(d["concentration_mgml"])
    ax.fill_between([0.8, 8], [5, 5], [100, 100], color="0.85", alpha=0.6, label="lit band 5-100 Pa")
    ax.plot(c, d["G_spring_Pa"], "-o", color="#1f77b4", lw=2, ms=6, label=f"athermal spring n={d['exponent_spring']:.2f}")
    ax.plot(c, d["G_wlc_Pa"], "-s", color="#d62728", lw=2, ms=6, label=f"thermal WLC n={d['exponent_wlc']:.2f}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("collagen c [mg/mL]"); ax.set_ylabel("G' [Pa]")
    ax.set_title("(5) thermal-WLC FIXES the c-exponent (1.11->2.63; confirms KB-1.30 thermal; PROVISIONAL)", fontsize=9)
    ax.legend(fontsize=7, loc="upper left"); ax.grid(True, which="both", alpha=0.25)


def main():
    os.makedirs(FIGS, exist_ok=True)
    fig, axs = plt.subplots(2, 3, figsize=(19.5, 9.5))
    panel_moduli(axs[0, 0]); panel_alignment(axs[0, 1])
    panel_stress_prop(axs[1, 0]); panel_sensing(axs[1, 1])
    panel_wlc(axs[0, 2]); axs[1, 2].set_visible(False)
    fig.suptitle("FF ECM library — native validation dashboard (all results vs literature/KB anchors; "
                 "unfitted, bands = acceptance oracles)\ncollagen-I · fibrin · PA · HA · Matrigel · agarose  ·  "
                 "6/6 real Pa · alignment anisotropy · directional stress channeling · motor-clutch biphasic",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = f"{FIGS}/ecm_validation_dashboard.png"
    fig.savefig(path, dpi=140); plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
