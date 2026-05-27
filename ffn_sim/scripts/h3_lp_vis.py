#!/usr/bin/env python
"""H.3 cortex L_p production visualization (visualize-at-closeout, CLAUDE.md).

Reads the L_p production trajectory written by ``h3_lp_gpu_production.py``
(``outputs/h3/production/lp_{scale}_{device}.npz``) and writes PNG figures
into ``outputs/h3/figs/``:

    fig_h3_lp_tangent_correlation.png   C(s) decay (per-snapshot + ensemble) + WLC fit
    fig_h3_lp_distribution.png          per-snapshot L_p vs KU-1.1 band
    fig_h3_lp_filament_configs.png      sample 3D filament conformations

Visualization integrity (CLAUDE.md): per-realisation thin lines + ensemble mean
overlay, KU-1.1 reference band overlaid, SI/μm units annotated, no axis
truncation. Created 2026-05-28 per PI directive.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.common.filament_math import fit_persistence_length, tangent_correlation

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"
PROD = PKG / "outputs" / "h3" / "production"
FIGS = PKG / "outputs" / "h3" / "figs"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scale", default="medium")
    ap.add_argument("--device", default="gpu")
    args = ap.parse_args()

    with open(CFG) as f:
        cfg = yaml.safe_load(f)
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    l0 = p.rest_length
    lo, hi = p.L_p_band_m

    frames = np.load(PROD / f"lp_{args.scale}_{args.device}.npz")["frames"]  # (S,F,N,3)
    S, F, N, _ = frames.shape
    smax = N // 2
    FIGS.mkdir(parents=True, exist_ok=True)

    Lp, Cs = [], []
    for k in range(S):
        fit = fit_persistence_length(frames[k], rest_length=l0)
        if np.isfinite(fit.L_p_m):
            Lp.append(fit.L_p_m)
        Cs.append(tangent_correlation(frames[k], max_separation=smax))
    Lp = np.asarray(Lp)
    Cs = np.asarray(Cs)
    Lp_mean = float(Lp.mean())
    arclen_um = np.arange(smax + 1) * l0 * 1e6
    Cmean = Cs.mean(axis=0)

    # Fig 1 — tangent correlation C(s) + WLC fit
    fig, ax = plt.subplots(figsize=(7, 5))
    for k in range(S):
        ax.plot(arclen_um, Cs[k], color="0.8", lw=0.5, zorder=1)
    ax.plot(arclen_um, Cmean, "o-", color="C0", lw=2, zorder=3,
            label="ensemble mean C(s)")
    xf = np.linspace(0, arclen_um[-1], 200)
    ax.plot(xf, np.exp(-xf * 1e-6 / Lp_mean), "--", color="C3", lw=2, zorder=4,
            label=f"WLC fit  $L_p$={Lp_mean * 1e6:.1f} μm")
    ax.set_xlabel("arclength  s·ℓ₀  [μm]")
    ax.set_ylabel(r"C(s) = $\langle \hat{t}(0)\cdot\hat{t}(s)\rangle$")
    ax.set_title(f"H.3 cortex tangent correlation ({args.scale} GPU, F={F}, {S} snapshots)")
    ax.legend(); ax.grid(alpha=0.3)
    f1 = FIGS / "fig_h3_lp_tangent_correlation.png"
    fig.savefig(f1, dpi=140, bbox_inches="tight"); plt.close(fig)

    # Fig 2 — per-snapshot L_p vs KU-1.1 band
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.axhspan(lo * 1e6, hi * 1e6, color="C2", alpha=0.2,
               label=f"KU-1.1 band [{lo * 1e6:.1f}, {hi * 1e6:.1f}] μm")
    ax.plot(np.arange(len(Lp)), Lp * 1e6, "o", color="C0", ms=4,
            label="per-snapshot $L_p$")
    ax.axhline(Lp_mean * 1e6, color="C3", lw=2,
               label=f"ensemble mean {Lp_mean * 1e6:.2f} μm")
    ax.set_xlabel("snapshot index")
    ax.set_ylabel("$L_p$  [μm]")
    ax.set_title(f"H.3 per-snapshot persistence length ({args.scale} GPU)")
    ax.legend(); ax.grid(alpha=0.3)
    f2 = FIGS / "fig_h3_lp_distribution.png"
    fig.savefig(f2, dpi=140, bbox_inches="tight"); plt.close(fig)

    # Fig 3 — sample 3D filament conformations (last snapshot)
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    last = frames[-1]
    sel = np.random.default_rng(0).choice(F, size=min(12, F), replace=False)
    for i in sel:
        xyz = last[i] * 1e6
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], lw=1.2)
    ax.set_xlabel("x [μm]"); ax.set_ylabel("y [μm]"); ax.set_zlabel("z [μm]")
    ax.set_title(f"H.3 sample cortex filaments ({len(sel)} of {F}, last snapshot)")
    f3 = FIGS / "fig_h3_lp_filament_configs.png"
    fig.savefig(f3, dpi=140, bbox_inches="tight"); plt.close(fig)

    print(f"FIGS_WRITTEN {f1.name} {f2.name} {f3.name} "
          f"L_p_mean_um={Lp_mean * 1e6:.3f} (band {lo * 1e6:.1f}-{hi * 1e6:.1f})")


if __name__ == "__main__":
    main()
