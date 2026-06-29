"""A4′: partial β1 uniformity — sweep the traction localization Lp (edge → uniform) for Lam4.

A4 found full-uniform β1 (Lp ≫ spheroid) is the RIGHT KNOB (it removes/reverses Lam4's
small-size penalty and lifts Lam4 above Pre at large R₀) but OVERSHOOTS: it flips the
size-dependence sign (PI Lam4 still DECREASES, corr −0.91; full-uniform gives +0.95) and adds
large seed variance. A4′ sweeps the traction screening length Lp between the edge value
(11 µm, A1) and the uniform sentinel to locate the PARTIAL uniformity that lifts Lam4's
magnitude / flattens the penalty WITHOUT reversing the slope sign — the physical Lam4.

Only the DISTRIBUTION (Lp) is swept; the magnitude stays the A1 clutch-kinetics value
(Lam4 1.53 nN). Cohesion / substrate / proliferation identical to A1/A4. OVERLAY-ONLY: the Lp
that best matches the PI Lam4 shape is REPORTED, never used to fit the model to PI (the PI
shape is a qualitative target — a decreasing curve with a lifted magnitude — not a fit).

For each Lp we report the emergent corr(R₀,A/A₀) (sign of the size-dependence) and the
median A/A₀ (magnitude), so the Lp→(slope, magnitude) trend is visible and the partial-
uniformity regime (slope still < 0, magnitude lifted) is identified.

Usage: python -m ffn_sim.scripts.layer2_a4prime_partial_uniformity
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.archive.hoomd_legacy.spheroid.ligand_traction import LP_EDGE, resolve_ligand_traction
from ffn_sim.archive.hoomd_legacy.spheroid.observables import effective_radius
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.archive.hoomd_legacy.spheroid.proliferation import run_growth_pooled
from ffn_sim.archive.hoomd_legacy.spheroid.substrate import resolve_substrate

_ROOT = Path(__file__).resolve().parents[1]
_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_FIG_DIR = _ROOT / "outputs" / "layer2" / "figs"
_CKPT = _ROOT / "outputs" / "layer2" / "a4prime_partial_uniformity.checkpoint.json"

# Lp ladder from edge (A1) through partial to ~uniform (A4). µm.
_LP_LADDER_UM = [11.0, 40.0, 120.0, 1000.0]
_N_LIST = [120, 400, 600]          # 3 R0 (small/mid/large) — enough for the slope sign
_N_SEEDS = 3                       # A4 showed high variance at uniform → average more


def main(argv: list[str] | None = None) -> int:
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved)
    sub = resolve_substrate(resolved, adhesion_ratio=1.0)
    f_lam = resolve_ligand_traction("Lam4").f_traction   # 1.53 nN (magnitude held)
    print(f"[A4′] Lam4 f_traction={f_lam*1e9:.2f} nN held; sweeping Lp (edge 11µm → uniform):")

    ck = json.loads(_CKPT.read_text()) if _CKPT.exists() else {}
    for lp_um in _LP_LADDER_UM:
        Lp = lp_um * 1e-6
        for n in _N_LIST:
            key = f"{lp_um:.0f}:{n}"
            if key in ck:
                continue
            r0r, aar, ej = [], [], 0
            for s in range(_N_SEEDS):
                res = run_growth_pooled(
                    resolved, prolif, n_cells_init=n, total_time=2.0 * prolif.cycle_time_mean,
                    epoch_steps=1200, settle_steps=1000, seed=3000 + s, max_cells=3000,
                    cohesion="catch", cad=cad, substrate=sub, f_traction=f_lam, Lp=Lp,
                )
                if res.get("ejected"):
                    ej += 1
                r0r.append(effective_radius(res["a0_core"]))
                aar.append(float(res["area_core_over_a0"][-1]))
            ck[key] = {"Lp_um": lp_um, "R0": float(np.mean(r0r)),
                       "aa": float(np.mean(aar)), "aa_sd": float(np.std(aar)), "ejected": ej}
            _CKPT.write_text(json.dumps(ck, indent=2))
            ejs = f" ⚠EJ{ej}/{_N_SEEDS}" if ej else ""
            print(f"  Lp={lp_um:6.0f}µm N0={n:4d} R0={ck[key]['R0']*1e6:5.1f}µm "
                  f"A/A0={ck[key]['aa']:.2f}±{ck[key]['aa_sd']:.2f}{ejs}", flush=True)

    # per-Lp slope (corr) + magnitude
    print("\n[A4′] Lp → (size-dependence sign, magnitude):")
    rows = []
    for lp_um in _LP_LADDER_UM:
        R0 = np.array([ck[f"{lp_um:.0f}:{n}"]["R0"] for n in _N_LIST])
        AA = np.array([ck[f"{lp_um:.0f}:{n}"]["aa"] for n in _N_LIST])
        corr = float(np.corrcoef(R0, AA)[0, 1])
        rows.append((lp_um, corr, float(np.median(AA)), AA))
        tag = "edge (A1)" if lp_um == 11.0 else ("~uniform (A4)" if lp_um >= 1000 else "partial")
        print(f"  Lp={lp_um:6.0f}µm  corr(R0,A/A0)={corr:+.2f}  med(A/A0)={np.median(AA):.2f}  "
              f"[{AA.min():.2f}–{AA.max():.2f}]  {tag}")
    # identify the partial regime: slope still negative (like PI −0.91) but magnitude lifted vs edge
    edge_med = rows[0][2]
    partial = [r for r in rows if r[1] < 0 and r[2] > edge_med + 0.05]
    print(f"\n[A4′] PI Lam4 target: slope NEGATIVE (corr≈−0.91) + magnitude LIFTED.")
    if partial:
        best = max(partial, key=lambda r: r[2])
        print(f"  → partial-uniformity regime found: Lp≈{best[0]:.0f}µm keeps slope {best[1]:+.2f} "
              f"(<0, like PI) while lifting med(A/A0) {edge_med:.2f}→{best[2]:.2f}.")
    else:
        print(f"  → NO intermediate Lp keeps slope<0 while lifting magnitude in this ladder "
              f"(the lift and the sign-flip are coupled here) — report honestly; finer Lp / more "
              f"seeds (B2/GPU) or a different uniformity model may be needed.")

    try:
        _make_figure(rows)
    except Exception as exc:  # noqa: BLE001
        print(f"[viz] skipped figure: {exc}")
    return 0


def _make_figure(rows) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.8, 5.2))
    lps = [r[0] for r in rows]
    corrs = [r[1] for r in rows]
    meds = [r[2] for r in rows]

    ax.plot(lps, corrs, "o-", color="tab:green")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.axhline(-0.91, color="crimson", ls="--", lw=1.2, label="PI Lam4 slope (corr ≈ −0.91)")
    ax.set_xscale("log")
    ax.set_xlabel("traction localization Lp (µm)  [edge 11 → uniform 1000]")
    ax.set_ylabel("emergent corr(R₀, A/A₀)  (size-dependence sign)")
    ax.set_title("A4′ — partial β1 uniformity: Lp → size-dependence sign\n"
                 "edge (1/R, −) → uniform (penalty reversed, +); PI target stays −", fontsize=9)
    ax.legend(fontsize=8)

    ax2.plot(lps, meds, "s-", color="tab:purple")
    ax2.axhline(meds[0], color="gray", ls=":", lw=1, label=f"edge magnitude ({meds[0]:.2f})")
    ax2.set_xscale("log")
    ax2.set_xlabel("traction localization Lp (µm)")
    ax2.set_ylabel("median A/A₀ (magnitude)")
    ax2.set_title("A4′ — Lp → magnitude (uniform engages interior → lift)\n"
                  "the partial regime = slope still <0 AND magnitude lifted", fontsize=9)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out = _FIG_DIR / "fig_layer2_a4prime_partial_uniformity.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
