"""FF ECM library — durotaxis stiffness-gradient substrates (spatially-varying Pa).

Cells migrate up stiffness gradients (durotaxis). This builds ECM substrates whose LOCAL Young's modulus
varies with position, then probes E(x) by indenting at a row of positions — validating that the substrate
reproduces the intended gradient in real Pa/µm against the in-vivo durotaxis range (KB-1.V.1.3: physiological
~1 Pa/µm, pathological ~10 Pa/µm, sharp interface ≥100 Pa/µm; Vincent-Engler 2013). A new capability of the
FF ECM library (spatially-graded continuum substrate), self-contained (no cell-crawl needed).

Run:  python -m aleph.scripts.ff_ecm_gradient [--device cpu|cuda:0]
Out:  aleph/outputs/ff/ecm_lib/{ecm_gradient.json, figs/durotaxis_gradient.png}
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


def probe_gradient(ecm, *, axis="x", n_probe=6, R_ind=10.0, relax_steps=3000, device="cpu"):
    """Indent at a row of positions along ``axis`` → local E_eff(x) profile [Pa] and fitted dE/dx [Pa/µm]."""
    lo, hi = ecm.box_lo, ecm.box_hi
    ax = {"x": 0, "y": 1, "z": 2}[axis]
    cy = 0.5 * (lo[1] + hi[1]) if ax == 0 else 0.5 * (lo[0] + hi[0])
    span = hi[ax] - lo[ax]
    xs = np.linspace(lo[ax] + 0.15 * span, hi[ax] - 0.15 * span, n_probe)
    Es, built = [], []
    for x in xs:
        cxy = (x, cy) if ax == 0 else (cy, x)
        ind = M.indentation_modulus(ecm, indenter_R_um=R_ind, max_depth_um=1.4, n_depths=4,
                                    k_ind=5.0e4, relax_steps=relax_steps, center_xy=cxy, device=device)
        Es.append(ind["E_eff_Pa"])
        frac = (x - lo[ax]) / span
        built.append(ecm.meta["E_lo_Pa"] + frac * (ecm.meta["E_hi_Pa"] - ecm.meta["E_lo_Pa"]))
    xs = np.asarray(xs); Es = np.asarray(Es); built = np.asarray(built)
    slope = float(np.polyfit(xs, Es, 1)[0])            # measured dE_eff/dx [Pa/µm]
    return {"x_um": xs.tolist(), "E_eff_Pa": Es.tolist(), "E_built_Pa": built.tolist(),
            "measured_grad_Pa_um": slope, "built_grad_Pa_um": ecm.meta["grad_Pa_per_um"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--steps", type=int, default=3000)
    a = ap.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    spec = L.get_spec("pa_gel")
    # long slab along x (gradient axis); calibrate bond-stiffness-per-modulus once on a uniform gel
    lo, hi = [0, 0, 0], [100.0, 34.0, 26.0]
    k_ref = M.calibrate_continuum_k(spec, lo, hi, node_spacing_um=2.5, probe="uniaxial",
                                    n_steps=a.steps, device=a.device, rng=np.random.default_rng(1))
    k_per_E = k_ref / spec.E_gel_Pa
    print(f"calibrated k_per_E = {k_per_E:.4f} pN/µm per Pa", flush=True)

    # KB-1.V.1.3 durotaxis regimes: physiological 1, pathological 10, sharp interface ~100 Pa/µm
    cases = [("physiological (1 Pa/µm)", 300.0, 300.0 + 1.0 * 100.0),
             ("pathological (10 Pa/µm)", 300.0, 300.0 + 10.0 * 100.0),
             ("sharp interface (~100 Pa/µm)", 300.0, 300.0 + 100.0 * 100.0)]
    results = []
    for name, E_lo, E_hi in cases:
        rng = np.random.default_rng(4)
        ecm = L.build_gradient_ecm(spec, lo, hi, E_lo_Pa=E_lo, E_hi_Pa=E_hi, axis="x",
                                   node_spacing_um=2.5, k_per_E=k_per_E, pin_faces=("z_lo",), rng=rng)
        pr = probe_gradient(ecm, axis="x", n_probe=6, R_ind=10.0, relax_steps=a.steps, device=a.device)
        pr["name"] = name
        # indentation returns ~c·E (continuum factor); the RELATIVE gradient ratio is the validation
        ratio = pr["measured_grad_Pa_um"] / max(pr["built_grad_Pa_um"], 1e-9)
        pr["measured_over_built"] = ratio
        results.append(pr)
        print(f"  {name}: built {pr['built_grad_Pa_um']:.1f} Pa/µm → measured {pr['measured_grad_Pa_um']:.1f} "
              f"Pa/µm ({ratio:.2f}× continuum factor)", flush=True)

    # figure
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["tab:green", "tab:orange", "tab:red"]
    for r, c in zip(results, colors):
        ax.plot(r["x_um"], r["E_eff_Pa"], "o-", color=c, label=f"{r['name']} — meas {r['measured_grad_Pa_um']:.1f} Pa/µm")
        ax.plot(r["x_um"], np.array(r["E_built_Pa"]) * r["measured_over_built"], "--", color=c, alpha=0.5)
    ax.set_xlabel("position along gradient axis x [µm]")
    ax.set_ylabel("local indentation modulus E_eff [Pa]")
    ax.set_yscale("log")
    ax.set_title("Durotaxis stiffness-gradient ECM substrates (KB-1.V.1.3 regimes)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/durotaxis_gradient.png", dpi=130)
    plt.close(fig)
    with open(f"{OUT}/ecm_gradient.json", "w") as f:
        json.dump({"k_per_E": k_per_E, "gradients": results}, f, indent=2, default=float)
    print(f"wrote {OUT}/ecm_gradient.json + figs/durotaxis_gradient.png")


if __name__ == "__main__":
    main()
