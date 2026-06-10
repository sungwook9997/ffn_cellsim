"""Validation + figure for the Slater radial stress(r) estimator.

Reproduces the estimator-audit discipline on a SYNTHETIC field whose answer is
known analytically, and emits the closeout figure (visualize-at-closeout rule).

Synthetic field: ``M`` isotropic radial rays of equal-tension ``T₀`` segments
(force-balanced) → every shell carries the same total radial force ``M·T₀`` →
``σ(r)·A(r) = M·T₀`` constant, i.e. ``σ ∝ r⁻²`` (sphere / 3D, Slater O1) and
``σ ∝ r⁻¹`` (cylinder / 2D). The figure overlays the measured σ(r) on those
reference power laws and shows the σ·A(r) flux-constancy panel.

With ``--with-cell`` it also smoke-runs :func:`radial_stress_from_sim` on a small
physiological-baseline MCF7 cell to exercise the live-sim path (printed only — a
thin cortex SHELL is not a volume-filling r⁻² field, so no scaling is asserted).

Run:
    PYTHONPATH=. python ffn_sim/scripts/h7_radial_stress_validate.py
    PYTHONPATH=. python ffn_sim/scripts/h7_radial_stress_validate.py --with-cell
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ffn_sim.cortex.radial_stress import (
    fit_power_law,
    radial_stress_profile,
)

PKG = Path(__file__).resolve().parents[1]
FIG_DIR = PKG / "outputs" / "h7" / "figs"

T0 = 1.0e-12          # 1 pN per segment
M_RAYS = 200
R_IN, R_OUT = 1.0e-6, 11.0e-6
N_BEADS = 41


def _isotropic_dirs(m: int) -> np.ndarray:
    i = np.arange(m, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / m)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.column_stack(
        [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)]
    )


def _radial_chain_field():
    dirs = _isotropic_dirs(M_RAYS)
    rb = np.linspace(R_IN, R_OUT, N_BEADS)
    rA_list, rB_list = [], []
    for d in dirs:
        for j in range(N_BEADS - 1):
            rA_list.append(rb[j] * d)
            rB_list.append(rb[j + 1] * d)
    rA = np.asarray(rA_list)
    rB = np.asarray(rB_list)
    dd = rB - rA
    u = dd / np.linalg.norm(dd, axis=1)[:, None]
    T = np.full(rA.shape[0], T0)
    r_mid = 0.5 * (rb[:-1] + rb[1:])
    return rA, rB, u, T, r_mid


def _smoke_cell() -> None:
    """Run radial_stress_from_sim on a small built cell (live-sim path)."""
    from copy import deepcopy

    from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
    from ffn_sim.cortex.radial_stress import radial_stress_from_sim

    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["cortex_overrides"] = {"cortex": {"n_filaments": 120, "demo_mode": True}}
    m["compartments"]["nucleus"]["n_beads"] = 300
    cell = build_baseline_cell(manifest=m, device=None)
    res = radial_stress_from_sim(cell.simulation, n_radii=24)
    print("\n[--with-cell] small MCF7 cell radial_stress_from_sim:")
    print(f"  cortical bonds counted : {res['n_bonds']} "
          f"(excluded {res['n_bonds_excluded']} of {res['n_bonds_total']})")
    print(f"  centre [µm]            : "
          f"{np.round(res['center'] * 1e6, 3)}")
    finite = np.isfinite(res["stress"])
    nz = finite & (res["n_crossing"] > 0)
    print(f"  shells with crossings  : {int(np.count_nonzero(nz))}/{len(nz)}")
    if np.any(nz):
        smax = res["stress"][nz]
        print(f"  |σ(r)| range [Pa]      : "
              f"[{np.abs(smax).min():.3g}, {np.abs(smax).max():.3g}]")
        print(f"  power-law fit          : exponent={res['fit']['exponent']:.3f} "
              f"(r²={res['fit']['r2']:.3f}, n={res['fit']['n_points']}) "
              f"— informational only (thin shell, not volume-filling)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--with-cell", action="store_true",
                    help="also smoke radial_stress_from_sim on a built cell")
    args = ap.parse_args()

    rA, rB, u, T, r_mid = _radial_chain_field()
    s_sph = radial_stress_profile(rA, rB, u, T, r_mid, geometry="sphere")
    h_cyl = 5.0e-6
    s_cyl = radial_stress_profile(
        rA, rB, u, T, r_mid, geometry="cylinder", height=h_cyl
    )
    flux_sph = s_sph * 4.0 * np.pi * r_mid ** 2
    flux_cyl = s_cyl * 2.0 * np.pi * r_mid * h_cyl
    fit_sph = fit_power_law(r_mid, s_sph)
    fit_cyl = fit_power_law(r_mid, s_cyl)

    print("Slater radial stress(r) estimator — synthetic validation")
    print(f"  M={M_RAYS} isotropic rays, T₀={T0:.2e} N, "
          f"r∈[{R_IN*1e6:.0f},{R_OUT*1e6:.0f}] µm")
    print(f"  sphere  : exponent={fit_sph['exponent']:.6f} "
          f"(ref −2), r²={fit_sph['r2']:.6f}; "
          f"σ·4πr² const? {np.allclose(flux_sph, M_RAYS*T0, rtol=1e-9)} "
          f"(={flux_sph.mean():.4e} N, target {M_RAYS*T0:.4e})")
    print(f"  cylinder: exponent={fit_cyl['exponent']:.6f} "
          f"(ref −1), r²={fit_cyl['r2']:.6f}; "
          f"σ·2πr·h const? {np.allclose(flux_cyl, M_RAYS*T0, rtol=1e-9)}")

    rm_um = r_mid * 1e6
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))

    # Panel 1: σ(r) vs power-law references (log-log).
    ax[0].loglog(rm_um, s_sph, "o", ms=4, color="#1f77b4",
                 label="measured σ(r), sphere")
    ref_sph = (M_RAYS * T0) / (4.0 * np.pi * r_mid ** 2)
    ax[0].loglog(rm_um, ref_sph, "-", color="#1f77b4", alpha=0.5,
                 label=r"$M T_0/4\pi r^2 \propto r^{-2}$ (O1, 3D)")
    ax[0].loglog(rm_um, s_cyl, "s", ms=4, color="#d62728",
                 label="measured σ(r), cylinder")
    ref_cyl = (M_RAYS * T0) / (2.0 * np.pi * r_mid * h_cyl)
    ax[0].loglog(rm_um, ref_cyl, "-", color="#d62728", alpha=0.5,
                 label=r"$M T_0/2\pi r h \propto r^{-1}$ (O1, 2D)")
    ax[0].set_xlabel("shell radius r [µm]")
    ax[0].set_ylabel("radial stress σ(r) [Pa]")
    ax[0].set_title("σ(r) vs Slater O1 power laws\n"
                    f"fit exponent: sphere {fit_sph['exponent']:.3f}, "
                    f"cyl {fit_cyl['exponent']:.3f}")
    ax[0].legend(fontsize=7, loc="lower left")
    ax[0].grid(True, which="both", alpha=0.25)

    # Panel 2: flux constancy σ·A(r) (linear) — the conservation check.
    ax[1].plot(rm_um, flux_sph * 1e12, "o-", ms=4, color="#1f77b4",
               label="σ·4πr² (sphere)")
    ax[1].plot(rm_um, flux_cyl * 1e12, "s-", ms=4, color="#d62728",
               label="σ·2πr·h (cylinder)")
    ax[1].axhline(M_RAYS * T0 * 1e12, ls="--", color="k", alpha=0.6,
                  label=f"M·T₀ = {M_RAYS*T0*1e12:.1f} pN")
    ax[1].set_xlabel("shell radius r [µm]")
    ax[1].set_ylabel("crossing radial force σ·A(r) [pN]")
    ax[1].set_title("force-balance / conservation check\n"
                    "(flat ⇒ correct area normalisation)")
    ax[1].legend(fontsize=7)
    ax[1].grid(True, alpha=0.25)
    # Headroom so the flat lines are not visually clipped to the band.
    ax[1].set_ylim(0, M_RAYS * T0 * 1e12 * 1.6)

    fig.suptitle("Slater radial stress(r) estimator — synthetic validation "
                 "(ffn_sim.cortex.radial_stress)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = FIG_DIR / "h7_radial_stress_validation.png"
    fig.savefig(out, dpi=140)
    print(f"\nfigure → {out}")

    if args.with_cell:
        _smoke_cell()


if __name__ == "__main__":
    main()
