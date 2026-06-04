"""L2 Track-2 bridge — predicted aggregate spheroid surface tension from single-cell γ.

Closes the single-cell-γ → aggregate-σ bridge with the Winklbauer 2015 keystone relation,
driven by the IN-BAND passive cortical tension ``g_rigid = 0.57 mN/m`` (KU-3.5). This is the
scientifically productive Track-2 deliverable that is NOT blocked by the active-channel
(g_soft) force-aggregation floor: the bridge runs off the in-band structural channel.

The relation (Winklbauer 2015, J Cell Sci 128:3687, ``references/jcs174623.pdf``, SE311 —
verified directly from the PDF):

    σ ≈ β − β*                                  (PDF lines 403/541/563/572; Brodland DITH)
    β* = β_cont − Γ/2                            (PDF lines 163/180/209; Brodland-Chen 2000)

where β = single-cell cortical tension at the FREE surface (PDF line 178) = our g_rigid, β* =
the residual tension at a cell–cell contact (lowered by the adhesion tension Γ), and σ = the
emergent aggregate (tissue) surface tension. Adhesion lowers β* → so σ rises toward β as the
aggregate becomes more cohesive (β* → 0 ⇒ σ → β, the Roffay "outer = free cortex" identity
the existing oracle's ``aggregate_surface_tension`` encodes — that identity is the
STRONG-ADHESION endpoint of this more general relation).

Because the MCF7-specific β*/β is not measured here, we report σ over the adhesion fraction
φ = β*/β ∈ [0, ~0.625], marking three literature anchors:

  - φ = 0           → σ = β = 0.57 mN/m   (max-adhesion / Roffay σ=γ identity endpoint)
  - φ ≈ 0.25        → σ ≈ ¾β = 0.43 mN/m  (David et al. 2014 typical β*≈¼β, cited in Winklbauer)
  - φ ∈ [0.5,0.625] → σ ∈ [0.21,0.29]     (Roffay outer/interior ratio β/β* ∈ [1.6,2.0]
                                            ⇒ β*/β ∈ [0.5,0.625] ⇒ σ/β = 1−φ ∈ [0.375,0.5])

So the in-band g_rigid predicts an aggregate spheroid σ bracketed in ~[0.21, 0.57] mN/m. This
is compared to (a) the emergent catch-bond CBM measurement (rest σ≈0, strain-engaged ~0.016
mN/m — the stretch-latch, ``emergent_sigma_catch.json``), (b) the single-epithelial-cell γ_CST
≈ 0.013 mN/m (Warmt 2021, MCF-10A, SE307 — a lower magnitude scale), and (c) the tissue
tensiometry PROXY 21 mN/m (MCF10DCIS, Nagle 2022 — different line, flagged proxy).

This is a MEASUREMENT/prediction analysis: it imports the runtime-forbidden bridge oracle
ONLY for the published constants/relations (never a cell-build path).

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.layer2_winklbauer_bridge
Outputs: ffn_sim/outputs/layer2/winklbauer_bridge.{json,png}
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

# Oracle import is PREDICTION-only (published anchors); never a runtime cell-build path.
from ffn_sim.validation.oracles.spheroid import surface_tension_bridge as br

_ROOT = Path(__file__).resolve().parents[1]
_OUT_DIR = _ROOT / "outputs" / "layer2"

#: β = single-cell cortical tension at the free surface = in-band g_rigid (KU-3.5) [mN/m].
BETA_MN_M: float = br.G_RIGID_NATIVE_MN_M  # 0.57

#: David 2014 typical residual-tension fraction β*/β (cited in Winklbauer 2015).
DAVID_PHI: float = 0.25
#: Roffay 2021 outer/interior ratio β/β* → φ = β*/β window.
ROFFAY_PHI_WINDOW: tuple[float, float] = (
    1.0 / br.ROFFAY_SURFACE_INTERIOR_RATIO[1],   # 1/2.0 = 0.5
    1.0 / br.ROFFAY_SURFACE_INTERIOR_RATIO[0],   # 1/1.6 = 0.625
)
#: Warmt 2021 single-epithelial-cell cortical tension γ_CST (MCF-10A, SE307) [mN/m].
WARMT_GAMMA_CST_MN_M: float = 0.0129


def winklbauer_sigma(beta: float, phi: float) -> float:
    """Aggregate surface tension σ = β − β* = β·(1 − φ), φ = β*/β (Winklbauer 2015)."""
    return beta * (1.0 - phi)


def _load_emergent_catch() -> dict | None:
    p = _OUT_DIR / "emergent_sigma_catch.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    cb = d["catch_bond"]
    return {
        "rest_mN_m": d["ensemble"]["sigma_mean_mN_m"],
        "rest_std_mN_m": d["ensemble"]["sigma_std_mN_m"],
        "strain_grid": cb["strain_grid"],
        "sigma_strain_mN_m": cb["sigma_strain_mN_m"],
    }


def make_figure(res: dict, out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    phi = np.array(res["phi_grid"])
    sig = np.array(res["sigma_phi_mN_m"])
    beta = res["beta_mN_m"]

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.0))

    # Panel A: σ(φ) = β(1−φ) with the three literature anchors + KU band.
    ax[0].axhspan(*br.CORTICAL_TENSION_BAND_MN_M, color="tab:green", alpha=0.12,
                  label=f"KU-3.5 γ band {br.CORTICAL_TENSION_BAND_MN_M}")
    ax[0].plot(phi, sig, "-", color="tab:blue", lw=2, label="σ = β(1−φ)  (Winklbauer)")
    ax[0].plot([0.0], [beta], "o", color="tab:green", ms=9,
               label=f"φ=0: σ=β={beta:.2f} (Roffay σ=γ endpoint)")
    ax[0].plot([DAVID_PHI], [winklbauer_sigma(beta, DAVID_PHI)], "s", color="tab:orange", ms=9,
               label=f"φ=¼ (David'14): σ=¾β={winklbauer_sigma(beta, DAVID_PHI):.2f}")
    ax[0].axvspan(*ROFFAY_PHI_WINDOW, color="tab:red", alpha=0.12,
                  label=f"Roffay window σ∈[{res['sigma_roffay_lo']:.2f},{res['sigma_roffay_hi']:.2f}]")
    ax[0].set_xlabel("adhesion fraction  φ = β*/β"); ax[0].set_ylabel("aggregate σ (mN/m)")
    ax[0].set_title("A  Winklbauer bridge: σ = β − β*  from in-band g_rigid\n"
                    f"β = single-cell γ (g_rigid) = {beta:.2f} mN/m")
    ax[0].legend(fontsize=7, loc="upper right"); ax[0].grid(alpha=0.3)
    ax[0].set_ylim(0.0, max(0.7, beta * 1.05))

    # Panel B: predicted σ range vs the emergent CBM measurement + magnitude scales (log).
    labels, vals, colors = [], [], []
    labels.append("Winklbauer\npredicted range"); vals.append(res["sigma_mid_mN_m"]); colors.append("tab:blue")
    cb = res.get("emergent_catch")
    if cb:
        labels.append("emergent catch\n(rest)"); vals.append(max(cb["rest_mN_m"], 1e-4)); colors.append("tab:purple")
        labels.append("emergent catch\n(ε=10% strain)"); vals.append(cb["sigma_strain_mN_m"][-1]); colors.append("tab:cyan")
    labels.append("Warmt single-cell\nγ_CST (MCF-10A)"); vals.append(WARMT_GAMMA_CST_MN_M); colors.append("tab:olive")
    labels.append("tissue proxy\nMCF10DCIS"); vals.append(br.TISSUE_SIGMA_PROXY_MN_M); colors.append("tab:gray")
    x = np.arange(len(labels))
    ax[1].bar(x, vals, 0.6, color=colors)
    # predicted range error bar
    ax[1].errorbar([0], [res["sigma_mid_mN_m"]],
                   yerr=[[res["sigma_mid_mN_m"] - res["sigma_roffay_lo"]],
                         [beta - res["sigma_mid_mN_m"]]], fmt="none", ecolor="k", capsize=6)
    ax[1].axhspan(*br.CORTICAL_TENSION_BAND_MN_M, color="tab:green", alpha=0.12)
    ax[1].set_yscale("log"); ax[1].set_xticks(x); ax[1].set_xticklabels(labels, fontsize=7)
    ax[1].set_ylabel("surface tension (mN/m, log)")
    ax[1].set_title("B  predicted aggregate σ vs emergent measurement + magnitude scales")
    ax[1].grid(alpha=0.3, axis="y")

    fig.suptitle("L2 Track-2 — single-cell γ (g_rigid) → aggregate spheroid σ bridge "
                 "(Winklbauer 2015 σ=β−β*)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def main() -> int:
    beta = BETA_MN_M
    phi_grid = np.linspace(0.0, 0.625, 64)
    sigma_phi = [winklbauer_sigma(beta, float(p)) for p in phi_grid]

    sigma_roffay_lo = winklbauer_sigma(beta, ROFFAY_PHI_WINDOW[1])  # largest φ → smallest σ
    sigma_roffay_hi = winklbauer_sigma(beta, ROFFAY_PHI_WINDOW[0])
    sigma_david = winklbauer_sigma(beta, DAVID_PHI)
    sigma_mid = sigma_david  # the David-typical mid estimate

    res = {
        "beta_mN_m": beta,
        "relation": "sigma = beta - beta* = beta*(1 - phi),  phi = beta*/beta  (Winklbauer 2015, SE311)",
        "phi_grid": phi_grid.tolist(),
        "sigma_phi_mN_m": sigma_phi,
        "anchors": {
            "phi0_sigma_eq_beta_mN_m": beta,
            "david_phi": DAVID_PHI,
            "david_sigma_mN_m": sigma_david,
            "roffay_phi_window": list(ROFFAY_PHI_WINDOW),
            "roffay_sigma_window_mN_m": [sigma_roffay_lo, sigma_roffay_hi],
        },
        "sigma_predicted_range_mN_m": [sigma_roffay_lo, beta],   # [Roffay-window low, σ=γ endpoint]
        "sigma_mid_mN_m": sigma_mid,
        "sigma_roffay_lo": sigma_roffay_lo,
        "sigma_roffay_hi": sigma_roffay_hi,
        "ku35_band_mN_m": list(br.CORTICAL_TENSION_BAND_MN_M),
        "in_band": bool(br.CORTICAL_TENSION_BAND_MN_M[0] <= sigma_mid <= br.CORTICAL_TENSION_BAND_MN_M[1]),
        "comparisons": {
            "warmt_single_cell_gamma_cst_mN_m": WARMT_GAMMA_CST_MN_M,
            "tissue_proxy_MCF10DCIS_mN_m": br.TISSUE_SIGMA_PROXY_MN_M,
        },
        "emergent_catch": _load_emergent_catch(),
    }

    print(f"[bridge] β (single-cell γ = g_rigid) = {beta:.3f} mN/m  (KU-3.5 in-band)")
    print(f"[bridge] σ = β − β* = β(1−φ), φ=β*/β  (Winklbauer 2015, SE311, PDF-verified)")
    print(f"  φ=0    (max adhesion / Roffay σ=γ endpoint): σ = {beta:.3f} mN/m")
    print(f"  φ=¼    (David 2014 typical β*≈¼β):          σ = {sigma_david:.3f} mN/m  (=¾β)")
    print(f"  φ∈[{ROFFAY_PHI_WINDOW[0]:.3f},{ROFFAY_PHI_WINDOW[1]:.3f}] (Roffay window):     "
          f"σ ∈ [{sigma_roffay_lo:.3f}, {sigma_roffay_hi:.3f}] mN/m")
    print(f"  ⇒ predicted aggregate σ ≈ [{sigma_roffay_lo:.2f}, {beta:.2f}] mN/m "
          f"(mid {sigma_mid:.2f}); KU band {tuple(br.CORTICAL_TENSION_BAND_MN_M)} — "
          f"David-mid {'IN' if res['in_band'] else 'OUT of'} band")
    cb = res["emergent_catch"]
    if cb:
        print(f"[compare] emergent catch-bond CBM: rest σ={cb['rest_mN_m']:+.4f} mN/m, "
              f"ε=10% strain σ={cb['sigma_strain_mN_m'][-1]:+.4f} mN/m (stretch-latch ≪ predicted)")
    print(f"[compare] Warmt single-cell γ_CST (MCF-10A) = {WARMT_GAMMA_CST_MN_M:.4f} mN/m; "
          f"tissue proxy (MCF10DCIS) = {br.TISSUE_SIGMA_PROXY_MN_M:.1f} mN/m")

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_json = _OUT_DIR / "winklbauer_bridge.json"
    out_json.write_text(json.dumps(res, indent=2))
    print(f"\n[json] {out_json}")
    try:
        out_png = _OUT_DIR / "winklbauer_bridge.png"
        make_figure(res, out_png)
        print(f"[viz]  {out_png}")
    except Exception as exc:  # best-effort
        print(f"[viz]  skipped figure: {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
