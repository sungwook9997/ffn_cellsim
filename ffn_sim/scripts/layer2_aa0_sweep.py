"""L2.3: emergent A/A0(R0) spreading law from the CBM, fit to A/A0 = a + b/R + c/R^2 (G3).

Sweeps initial spheroid size R0 (via cell count), drives active-wetting spreading
(edge-directed traction vs cohesion), measures the steady spread ratio A/A0, and fits the
PI's novel law. The 1/R, 1/R^2 terms should emerge from the edge/bulk cell-fraction (~1/R),
NOT be hard-coded. PI poster A/A0 values are overlay-only at comparison time (we hold only
the functional form). MCF7 is cohesion-dominated, so the spread signal is small — this is a
first-pass harness; clean extraction may need longer runs / ensemble averaging (reported).

Usage: python -m ffn_sim.scripts.layer2_aa0_sweep
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.observables import effective_radius
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2
from ffn_sim.archive.hoomd_legacy.spheroid.spreading import run_edge_spreading
from ffn_sim.validation.oracles.spheroid.aa0_law import aa0_model, fit_aa0

_ROOT = Path(__file__).resolve().parents[1]
_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_OCFG = _ROOT / "validation" / "oracles" / "configs" / "layer2_cbm.yaml"
_FIG_DIR = _ROOT / "outputs" / "layer2" / "figs"


def main(argv: list[str] | None = None) -> int:
    resolved = resolve_layer2(yaml.safe_load(_CFG.read_text()))
    g3 = yaml.safe_load(_OCFG.read_text())["spheroid"]["acceptance"]["g3"]

    n_list = [60, 120, 250, 450, 650]
    n_seeds = 8
    # f_traction kept safely below the 6.5 nN cohesion (further from the detachment threshold
    # = stable spreading, lower realization variance); ensemble-averaged over n_seeds.
    f_tr, Lp = 5.0e-9, 11.0e-6
    print(f"[sweep] f_traction={f_tr*1e9:.1f} nN, Lp={Lp*1e6:.0f} µm, cohesion="
          f"{resolved.D_e*resolved.morse_alpha/2*1e9:.1f} nN, n_seeds={n_seeds}")

    R0s, AA0m, AA0e, dets = [], [], [], []
    for n in n_list:
        r0r, aar, detr = [], [], []
        for s in range(n_seeds):
            res = run_edge_spreading(
                resolved, n_cells=n, f_traction=f_tr, Lp=Lp,
                settle_steps=1_500, spread_steps=8_000, recompute_every=400, seed=1000 + s,
            )
            r0r.append(effective_radius(res["a0"]))
            aar.append(float(res["area_over_a0"][-1]))
            detr.append(res["detached_fraction_final"])
        R0i, m, e = float(np.mean(r0r)), float(np.mean(aar)), float(np.std(aar))
        R0s.append(R0i); AA0m.append(m); AA0e.append(e); dets.append(max(detr))
        print(f"  N={n:4d}  R0={R0i*1e6:6.1f} µm  A/A0={m:.4f}±{e:.4f}  detached={max(detr):.4f}")

    R0 = np.array(R0s); AA0 = np.array(AA0m); AA0err = np.array(AA0e)
    fit = fit_aa0(R0, AA0)
    a, b, c, r2 = fit["a"], fit["b"], fit["c"], fit["r_squared"]
    print(f"\n[fit] A/A0 = {a:.4f} + ({b*1e6:.4f} µm)/R + ({c*1e12:.4f} µm²)/R²   r²={r2:.3f}")

    ok_r2 = r2 >= g3["aa0_fit_r_squared_min"]
    ok_n = len(n_list) >= g3["n_radii_min"]
    print(f"[G3 fit r² ≥ {g3['aa0_fit_r_squared_min']}] {'PASS' if ok_r2 else 'FAIL'}  (r²={r2:.3f})")
    print(f"[G3 ≥{g3['n_radii_min']} radii        ] {'PASS' if ok_n else 'FAIL'}  ({len(n_list)})")
    signal = AA0.max() - AA0.min()
    print(f"[signal A/A0 spread] {signal:.4f}  "
          f"({'measurable' if signal > 0.02 else 'SMALL — needs longer runs / ensemble avg'})")

    try:
        _make_figure(R0, AA0, AA0err, fit, f_tr, signal)
    except Exception as exc:
        print(f"[viz] skipped figure: {exc}")
    return 0 if (ok_r2 and ok_n) else 1


def _make_figure(R0, AA0, AA0err, fit, f_tr, signal) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    Rg = np.linspace(R0.min() * 0.9, R0.max() * 1.05, 200)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.errorbar(R0 * 1e6, AA0, yerr=AA0err, fmt="o", ms=8, color="steelblue", capsize=4,
                zorder=3, label="CBM (mean ± sd, ensemble)")
    ax.plot(Rg * 1e6, aa0_model(Rg, fit["a"], fit["b"], fit["c"]), "-", color="crimson",
            label=f"fit a+b/R+c/R²  (r²={fit['r_squared']:.3f})")
    ax.axhline(fit["a"], color="gray", ls=":", lw=1, label=f"a (baseline) = {fit['a']:.3f}")
    ax.set_xlabel("initial effective radius R₀ (µm)")
    ax.set_ylabel("A / A₀ (spread ratio)")
    ax.set_title(
        f"Layer-2 L2.3 — emergent MCF7 spreading law (f_traction={f_tr*1e9:.0f} nN)\n"
        f"A/A0 = {fit['a']:.3f} + ({fit['b']*1e6:.3f} µm)/R + ({fit['c']*1e12:.3f} µm²)/R²"
        f"   |  signal Δ(A/A0)={signal:.3f}"
        + ("" if signal > 0.02 else "  (small — cohesion-dominated MCF7; longer runs/ensemble for clean law)"),
        fontsize=9,
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = _FIG_DIR / "fig_layer2_aa0_law.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
