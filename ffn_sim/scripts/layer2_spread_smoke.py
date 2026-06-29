"""L2.2 spreading-mechanism smoke + figure: verify f=0 reduces to G1, and the active force
drives the spread/detach response; auto-generates a figure (per the production-driver-auto-viz
rule).

Mechanism check only. With the corrected nN cohesion (Iturri 2020), MCF7 is appropriately
cohesive: traction below the measured detachment force does NOT shred it (correct low-invasion
behaviour). A QUANTITATIVE active-wetting A/A0 sweep needs the per-cell MIGRATION drag anchored
(water-Stokes γ is wrong for motility) + longer edge-driven runs — that lands as the L2.3 driver.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2
from ffn_sim.archive.hoomd_legacy.spheroid.spreading import run_spreading

_ROOT = Path(__file__).resolve().parents[1]
_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_FIG_DIR = _ROOT / "outputs" / "layer2" / "figs"


def main(argv: list[str] | None = None) -> int:
    resolved = resolve_layer2(yaml.safe_load(_CFG.read_text()))
    f_coh = resolved.D_e * resolved.morse_alpha / 2.0  # Morse max |F| = cohesion scale
    print(f"[scale] Morse max |F| = {f_coh*1e9:.2f} nN  (measured MCF7-MCF7 de-adhesion anchor)")

    common = dict(n_cells=120, settle_steps=5_000, spread_steps=15_000, n_samples=4)
    forces_nN = [0.0, 2.0, 4.0, 6.0]
    runs = []
    for fn in forces_nN:
        r = run_spreading(resolved, f_active=fn * 1e-9, **common)
        runs.append((fn, r))
        print(f"  f_active={fn:5.1f} nN  ->  A/A0={r['area_over_a0'][-1]:.3f}  "
              f"detached={r['detached_fraction_final']:.4f}")

    base = runs[0][1]["area_over_a0"][-1]
    ok_boundary = 0.9 <= base <= 1.15
    ok_connected = runs[-1][1]["detached_fraction_final"] < 0.10
    print(f"\n[boundary f=0 -> A/A0~1 ] {'PASS' if ok_boundary else 'FAIL'}  (A/A0={base:.3f})")
    print(f"[cohesion holds (MCF7)   ] {'PASS' if ok_connected else 'FAIL'}  "
          f"(detached={runs[-1][1]['detached_fraction_final']:.3f} @ {forces_nN[-1]:.0f} nN)")

    try:
        _make_figure(resolved, forces_nN, runs, f_coh)
    except Exception as exc:  # best-effort; never fail the gate on a viz error
        print(f"[viz] skipped figure: {exc}")
    return 0 if (ok_boundary and ok_connected) else 1


def _make_figure(resolved, forces_nN, runs, f_coh) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    aa0 = [r["area_over_a0"][-1] for _, r in runs]
    det = [r["detached_fraction_final"] for _, r in runs]
    init = runs[0][1]["pos_settled"] * 1e6
    final = runs[-1][1]["pos_final"] * 1e6

    fig, ax = plt.subplots(2, 2, figsize=(11, 9))
    ax[0, 0].scatter(init[:, 0], init[:, 1], s=16, alpha=0.7, edgecolors="k", linewidths=0.3)
    ax[0, 0].set_aspect("equal"); ax[0, 0].set_title("settled aggregate (f_active = 0)")
    ax[0, 0].set_xlabel("x (µm)"); ax[0, 0].set_ylabel("y (µm)")
    ax[0, 1].scatter(final[:, 0], final[:, 1], s=16, alpha=0.7, edgecolors="k",
                     linewidths=0.3, color="darkorange")
    ax[0, 1].set_aspect("equal")
    ax[0, 1].set_title(f"under active traction (f_active = {forces_nN[-1]:.0f} nN)")
    ax[0, 1].set_xlabel("x (µm)"); ax[0, 1].set_ylabel("y (µm)")

    ax[1, 0].plot(forces_nN, aa0, "o-", color="steelblue")
    ax[1, 0].axhline(1.0, color="gray", ls=":", lw=1)
    ax[1, 0].axvline(f_coh * 1e9, color="crimson", ls="--", lw=2,
                     label=f"cohesion (detachment) ={f_coh*1e9:.1f} nN")
    ax[1, 0].set_xlabel("active traction f_active (nN)"); ax[1, 0].set_ylabel("A / A₀")
    ax[1, 0].set_title("spread area vs traction"); ax[1, 0].legend(fontsize=8)

    ax[1, 1].plot(forces_nN, det, "s-", color="firebrick")
    ax[1, 1].axvline(f_coh * 1e9, color="crimson", ls="--", lw=2)
    ax[1, 1].set_xlabel("active traction f_active (nN)"); ax[1, 1].set_ylabel("detached fraction")
    ax[1, 1].set_title("detachment vs traction"); ax[1, 1].set_ylim(-0.02, 1.02)

    fig.suptitle(
        f"Layer-2 CBM L2.2 — isotropic motility on a cohesive aggregate "
        f"(cohesion {f_coh*1e9:.1f} nN, γ={resolved.gamma_cell:.2f} N·s/m, "
        f"r₀={resolved.morse_r0*1e6:.0f} µm)\n"
        f"FINDING: isotropic self-propulsion does NOT spread a cohesive cluster (detached≈0); "
        f"the spreading driver is EDGE-directed traction (active wetting) — see L2.2-quant",
        fontsize=9,
    )
    fig.tight_layout()
    out = _FIG_DIR / "fig_layer2_l2_2_motility_mechanism.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
