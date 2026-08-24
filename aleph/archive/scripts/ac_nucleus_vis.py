"""Regenerate the nucleus (I2) analytic-oracle figures — the ac/ visualization pattern (copies fluid-spine).

The single regeneration entry-point for the nucleus track (extend it as I7 oracles land). Every figure
OVERLAYS the closed-form oracle/band on the numeric result, annotates SI units, and never truncates an
axis (the professor's visualization-integrity rules). Pure numpy/matplotlib on the committed host
oracles — runs on the dev Mac (no Warp/CUDA).

    python aleph/scripts/ac_nucleus_vis.py    ->  aleph/outputs/ac/nucleus/figs/*.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from aleph.components.nucleus.geometry import build_oblate_mesh, build_hinges, mesh_volume, oblate_volume
from aleph.components.nucleus.lamina_analytic import (
    LaminaParams,
    calibrate_kappa_tilde,
    dihedral_bending_energy,
    dihedral_bending_forces,
    lamina_tangent_modulus,
    lamina_tension,
    sphere_willmore_energy,
)
from aleph.components.nucleus.linc_analytic import (
    capstan_line_integral,
    oblate_equatorial_radius_at_constant_volume,
)

KAPPA_NE = 0.0828        # pN·µm (KB-DRAFT-3.B-01, ~20 kBT)
E_NUC_PA = 399.0         # Pa (KB-DRAFT-3.B-05, MCF7 in-situ, low-strain anchor)
OUT = Path(__file__).resolve().parents[1] / "outputs" / "ac" / "nucleus" / "figs"

# framework-#6 split — LABELED structural demo params (the numerical split is I0-B2 GAP-PI; the ratio
# is REPORTED not tuned). k_chrom+k_lamin_b = soft small-strain; k_lamin_ac = strain-stiffening.
DEMO = LaminaParams(k_chrom=E_NUC_PA * 0.015 * 0.6, k_lamin_b=E_NUC_PA * 0.015 * 0.4,
                    k_lamin_ac=E_NUC_PA * 0.015 * 3.0, knee_strain=0.10, eps_rupture=0.50,
                    kappa_ne=KAPPA_NE)


def fig_bending() -> Path:
    """Sphere→8πκ (κ̃-calibrated, R-/resolution-invariant Σ_ref) + the FD-gradient force sign arbiter."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    subs = [1, 2, 3, 4]
    sig = [dihedral_bending_energy(v, build_hinges(f), 1.0)
           for v, f in (build_oblate_mesh(3.0, 1.0, s) for s in subs)]
    ax0.axhline(7.50, ls="--", c="crimson", label="icosphere Σ_ref ≈ 7.50 (topological const)")
    ax0.plot(subs, sig, "o-", lw=2, label="raw dihedral sum Σ(1−n̂₁·n̂₂)")
    for r in (1.5, 6.0):
        s3 = dihedral_bending_energy(*(lambda m: (m[0], build_hinges(m[1])))(build_oblate_mesh(r, 1.0, 3)), 1.0)
        ax0.scatter([3], [s3], marker="x", s=80, label=f"R={r} µm (R-invariant)")
    ax0.set_xlabel("icosphere subdivision level  [–]")
    ax0.set_ylabel("dihedral sum  Σ  [–]")
    ax0.set_title(f"Σ_ref grid-/R-invariant → κ̃=8πκ/Σ_ref\n(calibrated E = 8πκ = {sphere_willmore_energy(KAPPA_NE):.4f} pN·µm exactly)")
    ax0.set_ylim(7.3, 7.7); ax0.legend(loc="lower right", fontsize=8)

    rng = np.random.default_rng(1)
    v, f = build_oblate_mesh(3.0, 1.4, 3)
    v = v + 0.04 * rng.standard_normal(v.shape)
    h = build_hinges(f)
    v0, f0 = build_oblate_mesh(3.0, 1.0, 3)
    kt = calibrate_kappa_tilde(KAPPA_NE, v0, build_hinges(f0))
    f_ana = dihedral_bending_forces(v, h, kt)
    eps = 1e-7
    f_fd = np.zeros_like(v)
    for i in range(v.shape[0]):
        for d in range(3):
            vp = v.copy(); vp[i, d] += eps
            vm = v.copy(); vm[i, d] -= eps
            f_fd[i, d] = -(dihedral_bending_energy(vp, h, kt) - dihedral_bending_energy(vm, h, kt)) / (2 * eps)
    lim = np.max(np.abs(f_ana)) * 1.05
    ax1.plot([-lim, lim], [-lim, lim], "crimson", lw=1, label="y = x (f = −∂E/∂x)")
    ax1.scatter(f_ana.ravel(), f_fd.ravel(), s=6, alpha=0.5, label="per-DOF force")
    ax1.set_xlabel("analytic bending force  [pN]")
    ax1.set_ylabel("finite-difference −∂E/∂x  [pN]")
    ax1.set_title(f"FD-gradient sign arbiter (Σf={np.abs(f_ana.sum(0)).max():.1e} pN)")
    ax1.set_xlim(-lim, lim); ax1.set_ylim(-lim, lim); ax1.legend(loc="upper left", fontsize=8)
    fig.suptitle("I2 nucleus — dihedral Helfrich bending (reuses the ff/membrane kernel)")
    fig.tight_layout()
    p = OUT / "i2_bending.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_lamin_knee() -> Path:
    """framework-#6 lamin-split tension σ(ε) + the tangent-modulus crossover at the knee."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    e = np.linspace(0.0, 0.45, 400)
    ax0.plot(e, lamina_tension(e, DEMO), lw=2, label="σ(ε) lamin-split tension")
    ax0.axvline(DEMO.knee_strain, ls=":", c="gray")
    ax0.annotate("lamin-A/C knee", (DEMO.knee_strain, lamina_tension(DEMO.knee_strain, DEMO)),
                 (0.18, lamina_tension(0.05, DEMO)), arrowprops=dict(arrowstyle="->", color="gray"))
    ax0.fill_between(e, 0, lamina_tension(e, DEMO), where=e <= DEMO.knee_strain, alpha=0.15,
                     label="chromatin + lamin-B (soft)")
    ax0.fill_between(e, 0, lamina_tension(e, DEMO), where=e > DEMO.knee_strain, alpha=0.15, color="orange",
                     label="lamin-A/C (strain-stiffening)")
    ax0.set_xlabel("areal strain  ε = (A−A₀)/A₀  [–]")
    ax0.set_ylabel("areal tension  σ  [pN/µm]")
    ax0.set_title("Lamin-split tension (KB-3.B2.2)"); ax0.legend(loc="upper left", fontsize=8)

    ax1.plot(e, lamina_tangent_modulus(e, DEMO), lw=2, label="tangent dσ/dε")
    ax1.axhline(DEMO.k_soft, ls="--", c="C0", label=f"K_soft = {DEMO.k_soft:.3f} (chrom+lamin-B)")
    ax1.axhline(DEMO.k_lamin_ac, ls="--", c="orange", label=f"K_laminAC = {DEMO.k_lamin_ac:.3f}")
    ax1.axvline(DEMO.knee_strain, ls=":", c="gray")
    ax1.set_xlabel("areal strain  ε  [–]")
    ax1.set_ylabel("tangent areal modulus  [pN/µm]")
    ax1.set_title(f"Modulus crossover ×{DEMO.stiffening_ratio:.2f} at the knee (REPORT, not tuned)")
    ax1.set_ylim(0, DEMO.k_lamin_ac * 1.2); ax1.legend(loc="center right", fontsize=8)
    fig.suptitle("I2 nucleus — framework-#6 chromatin→lamin-B→lamin-A/C areal elasticity")
    fig.tight_layout()
    p = OUT / "i2_lamin_knee.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_volume() -> Path:
    """Nucleoplasm volume: mesh→continuum convergence + oblate flatten at CONSTANT volume (I7 kinematics)."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    subs = [2, 3, 4, 5]
    v_cont = 4.0 / 3.0 * np.pi * 3.0 ** 3
    vols = [mesh_volume(*build_oblate_mesh(3.0, 1.0, s)) for s in subs]
    ax0.axhline(v_cont, ls="--", c="crimson", label=f"continuum (4/3)πR³ = {v_cont:.2f} µm³")
    ax0.plot(subs, vols, "o-", lw=2, label="signed-tet mesh volume")
    ax0.set_xlabel("icosphere subdivision level  [–]")
    ax0.set_ylabel("enclosed volume  V  [µm³]")
    ax0.set_title("Divergence-theorem volume → continuum"); ax0.legend(loc="lower right", fontsize=8)

    v0 = v_cont
    aspects = np.linspace(1.0, 3.0, 40)
    a_pred = np.array([oblate_equatorial_radius_at_constant_volume(v0, A) for A in aspects])
    v_recon = 4.0 / 3.0 * np.pi * a_pred ** 2 * (a_pred / aspects)
    ax1.plot(aspects, a_pred, lw=2, label="equatorial a(A) = (3V₀A/4π)^⅓ ∝ A^⅓")
    ax1.plot(aspects, 3.0 * aspects ** (1 / 3), "crimson", ls="--", label="3·A^⅓ (analytic)")
    ax1b = ax1.twinx()
    ax1b.plot(aspects, v_recon / v0, c="green", ls=":", label="V(A)/V₀ (conserved)")
    ax1b.set_ylabel("reconstructed volume  V/V₀  [–]", color="green")
    ax1b.set_ylim(0.9, 1.1); ax1b.axhline(1.0, c="green", lw=0.5)
    ax1.set_xlabel("oblate aspect  A = a/c  [–]")
    ax1.set_ylabel("equatorial radius  a  [µm]")
    ax1.set_title("Flatten at constant nucleoplasm volume (Khatau)"); ax1.legend(loc="upper left", fontsize=8)
    fig.suptitle("I2 nucleus — incompressible nucleoplasm (ν→½) volume conservation")
    fig.tight_layout()
    p = OUT / "i2_volume.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_rupture_capstan() -> Path:
    """EMERGENT rupture on/off past the threshold + the I7-pre-authored capstan ∮T·κ ds = 2T."""
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.2))
    e = np.linspace(0.0, 0.75, 500)
    ax0.plot(e, lamina_tension(e, DEMO), lw=2, label="σ(ε)")
    ax0.axvline(DEMO.eps_rupture, ls="--", c="crimson", label=f"rupture threshold ε={DEMO.eps_rupture} (GAP-PI)")
    ax0.fill_between(e, 0, lamina_tension(e, DEMO).max() * 1.05, where=e > DEMO.eps_rupture, alpha=0.12,
                     color="crimson", label="torn (σ→0)")
    ax0.set_xlabel("areal strain  ε  [–]"); ax0.set_ylabel("areal tension  σ  [pN/µm]")
    ax0.set_ylim(0, lamina_tension(e, DEMO).max() * 1.1)
    ax0.set_title("EMERGENT rupture (on/off; threshold not tuned)"); ax0.legend(loc="upper right", fontsize=8)

    T, R = 5.0, 3.0
    wraps = np.linspace(0.05, np.pi, 40)
    nets = []
    for phi in wraps:
        th = np.linspace(-phi / 2, phi / 2, 3000)
        ds = np.full_like(th, R * (th[1] - th[0]))
        kappa = np.full_like(th, 1.0 / R)
        inward = -np.column_stack([np.sin(th), np.cos(th), np.zeros_like(th)])  # toward centre; bisector=-y
        nets.append(capstan_line_integral(T, kappa, ds, inward, direction=(0.0, -1.0, 0.0)))
    ax1.plot(np.degrees(wraps), nets, lw=2, label="∮ T·κ (n̂·d̂) ds")
    ax1.plot(np.degrees(wraps), 2 * T * np.sin(wraps / 2), "crimson", ls="--",
             label="2T·sin(φ/2) (belt-over-cylinder)")
    ax1.axhline(2 * T, c="gray", ls=":"); ax1.annotate("full half-wrap → 2T", (150, 2 * T), (90, 2 * T * 0.8),
                                                        arrowprops=dict(arrowstyle="->", color="gray"))
    ax1.set_xlabel("cap wrap angle  φ  [deg]"); ax1.set_ylabel("net inward force  [pN]")
    ax1.set_title("Capstan: tangential cap tension → normal pressure (I7 pre-authored)")
    ax1.legend(loc="upper left", fontsize=8)
    fig.suptitle("I2 nucleus — EMERGENT rupture + I7 capstan load-path oracle")
    fig.tight_layout()
    p = OUT / "i2_rupture_capstan.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    figs = [fig_bending(), fig_lamin_knee(), fig_volume(), fig_rupture_capstan()]
    for f in figs:
        print(f"  wrote {f.relative_to(Path(__file__).resolve().parents[2])}")


if __name__ == "__main__":
    main()
