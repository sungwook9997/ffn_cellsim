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
    ap.add_argument("--exclude-transient", type=int, default=5,
                    help="early eq-transient snapshots to exclude from the "
                         "plateau L_p mean (FULL eq=100k leaves ~5).")
    args = ap.parse_args()

    with open(CFG) as f:
        cfg = yaml.safe_load(f)
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    l0 = p.rest_length
    lo, hi = p.L_p_band_m
    kT = p.kT
    k_theta = p.angle_k
    eq_target = p.equipartition_target_3d_kT
    eq_tol = p.equipartition_rel_tol

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

    # Fig 2 — per-snapshot L_p vs KU-1.1 band, eq-transient marked, plateau mean±sem
    nx = min(args.exclude_transient, max(0, len(Lp) - 2))
    plateau = Lp[nx:]
    p_mean = float(plateau.mean())
    p_sem = float(plateau.std(ddof=1) / np.sqrt(len(plateau))) if len(plateau) > 1 else 0.0
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.axhspan(lo * 1e6, hi * 1e6, color="C2", alpha=0.2,
               label=f"KU-1.1 band [{lo * 1e6:.1f}, {hi * 1e6:.1f}] μm")
    if nx > 0:
        ax.axvspan(-0.5, nx - 0.5, color="0.85", alpha=0.6,
                   label=f"eq-transient (first {nx}, excluded)")
        ax.plot(np.arange(nx), Lp[:nx] * 1e6, "x", color="0.5", ms=6)
    ax.plot(np.arange(nx, len(Lp)), plateau * 1e6, "o", color="C0", ms=4,
            label="plateau per-snapshot $L_p$")
    ax.axhline(p_mean * 1e6, color="C3", lw=2,
               label=f"plateau mean {p_mean * 1e6:.2f} ± {p_sem * 1e6:.2f} μm")
    ax.set_xlabel("snapshot index")
    ax.set_ylabel("$L_p$  [μm]")
    ax.set_title(f"H.3 persistence length ({args.scale} GPU, F={F}, {S} snap)")
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

    # --- bending statistics on the plateau (eq-transient excluded) ---
    pf = frames[nx:]                                   # (Sp, F, N, 3)
    bv = pf[:, :, 1:, :] - pf[:, :, :-1, :]            # bond vectors
    bn = bv / np.linalg.norm(bv, axis=-1, keepdims=True)
    cosang = -np.einsum("...i,...i->...", bn[:, :, :-1, :], bn[:, :, 1:, :])
    theta = np.arccos(np.clip(cosang, -1.0, 1.0)).ravel()   # interior angles
    E_kT = 0.5 * k_theta * (theta - np.pi) ** 2 / kT
    E_mean = float(E_kT.mean())
    E_sem = float(E_kT.std(ddof=1) / np.sqrt(E_kT.size))

    # Fig 4 — per-angle bending energy vs 3D equipartition target (0.9898 kT)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.hist(E_kT, bins=80, density=True, color="C0", alpha=0.6,
            label=r"per-angle $E_{bend}/k_BT$")
    ax.axvspan(eq_target * (1 - eq_tol), eq_target * (1 + eq_tol),
               color="C2", alpha=0.25,
               label=f"equipartition target {eq_target:.4f} ±{eq_tol*100:.0f}%")
    ax.axvline(E_mean, color="C3", lw=2,
               label=fr"$\langle E\rangle$={E_mean:.4f} ± {E_sem:.4f} $k_BT$")
    in_eq = abs(E_mean - eq_target) / eq_target <= eq_tol
    ax.set_xlabel(r"$E_{bend}/k_BT$ per angle"); ax.set_ylabel("density")
    ax.set_title(f"H.3 bending equipartition ({args.scale} GPU) — "
                 f"{'IN-band' if in_eq else 'OUT'}")
    ax.legend(); ax.grid(alpha=0.3)
    f4 = FIGS / "fig_h3_lp_bending_equipartition.png"
    fig.savefig(f4, dpi=140, bbox_inches="tight"); plt.close(fig)

    # Fig 5 — interior-angle PDF vs 3D Boltzmann (sin θ · Boltzmann)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.hist(theta, bins=120, density=True, color="C0", alpha=0.6,
            label="measured angle PDF")
    tg = np.linspace(theta.min(), np.pi - 1e-6, 600)
    w = np.sin(tg) * np.exp(-0.5 * k_theta * (tg - np.pi) ** 2 / kT)
    trapz = getattr(np, "trapezoid", None) or np.trapz
    w = w / trapz(w, tg)
    ax.plot(tg, w, "-", color="C3", lw=2,
            label=r"3D Boltzmann $\propto\sin\theta\,e^{-\frac{k_\theta}{2kT}(\theta-\pi)^2}$")
    ax.set_xlabel(r"interior angle $\theta$ [rad]"); ax.set_ylabel("density")
    ax.set_title(f"H.3 angle distribution vs 3D Boltzmann ({args.scale} GPU)")
    ax.legend(); ax.grid(alpha=0.3)
    f5 = FIGS / "fig_h3_lp_angle_pdf.png"
    fig.savefig(f5, dpi=140, bbox_inches="tight"); plt.close(fig)

    # Fig 6 — bending-energy-coloured 3D filament render (last plateau frame)
    last = pf[-1]
    # per-bead bending energy (interior beads), endpoints get nearest value.
    bvl = last[:, 1:, :] - last[:, :-1, :]
    bnl = bvl / np.linalg.norm(bvl, axis=-1, keepdims=True)
    cosl = -np.einsum("...i,...i->...", bnl[:, :-1, :], bnl[:, 1:, :])
    thl = np.arccos(np.clip(cosl, -1.0, 1.0))             # (F, N-2)
    Ebead = 0.5 * k_theta * (thl - np.pi) ** 2 / kT       # (F, N-2)
    Efull = np.zeros((F, N))
    Efull[:, 1:-1] = Ebead
    Efull[:, 0] = Ebead[:, 0]; Efull[:, -1] = Ebead[:, -1]
    fig = plt.figure(figsize=(8, 6.5))
    ax = fig.add_subplot(111, projection="3d")
    sel = np.random.default_rng(0).choice(F, size=min(20, F), replace=False)
    vmax = float(np.percentile(Efull[sel], 98)) or 1.0
    sc = None
    for i in sel:
        xyz = last[i] * 1e6
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color="0.7", lw=0.6, zorder=1)
        sc = ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=Efull[i],
                        cmap="inferno", vmin=0, vmax=vmax, s=18, zorder=2)
    if sc is not None:
        cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1)
        cb.set_label(r"$E_{bend}/k_BT$ per bead")
    ax.set_xlabel("x [μm]"); ax.set_ylabel("y [μm]"); ax.set_zlabel("z [μm]")
    ax.set_title(f"H.3 cortex filaments, bending-energy coloured "
                 f"({len(sel)} of {F})")
    f6 = FIGS / "fig_h3_lp_bending_3d.png"
    fig.savefig(f6, dpi=140, bbox_inches="tight"); plt.close(fig)

    print(
        f"FIGS_WRITTEN {f1.name} {f2.name} {f3.name} {f4.name} {f5.name} {f6.name} "
        f"plateau_L_p_um={p_mean * 1e6:.3f}±{p_sem * 1e6:.3f} "
        f"(band {lo * 1e6:.1f}-{hi * 1e6:.1f}) "
        f"E_bend_kT={E_mean:.4f}±{E_sem:.4f} (target {eq_target:.4f}±{eq_tol})"
    )


if __name__ == "__main__":
    main()
