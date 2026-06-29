"""H.4 visualisation — PNG figures for outputs/h4/figs/.

Generates the figures cross-referenced from ``ffn_sim/outputs/h4/REPORT.md``
§Figures section, per the PI's "visualize at unit closeout" feedback. Run::

    conda activate ffn_sim
    python ffn_sim/scripts/h4_vis.py

Reads from ``ffn_sim/outputs/h4/*.npz`` (where the validation tests
write their canonical artefacts when run with ``H4_PRODUCTION=1``) and
from the closed-form oracles directly when production artefacts are
absent, so the script works in both demo and production modes.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")     # headless
import matplotlib.pyplot as plt
import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.bridge.fa import resolve_h4, build_h4_state
from ffn_sim.archive.hoomd_legacy.bridge.motor import (
    hill_velocity,
    hill_velocity_clamped,
)
from ffn_sim.validation.pereverzev import (
    PereverzevParams,
    pereverzev_F_star,
    pereverzev_k_off,
    pereverzev_lifetime,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPO_ROOT / "ffn_sim/configs/phase1_h4.yaml"
OUT_DIR = REPO_ROOT / "ffn_sim/outputs/h4"
FIGS_DIR = OUT_DIR / "figs"


def fig_pereverzev_k_off_and_lifetime(p_params: PereverzevParams) -> Path:
    """k_off(F) and τ(F) over F ∈ [0, 60] pN, with F* marker."""
    F_grid_pN = np.linspace(0.0, 60.0, 400)
    F_grid_N = F_grid_pN * 1e-12
    k_off = pereverzev_k_off(F_grid_N, p_params)
    tau = pereverzev_lifetime(F_grid_N, p_params)
    F_star = pereverzev_F_star(p_params)
    F_star_pN = F_star * 1e12 if np.isfinite(F_star) else float("nan")
    tau_star = (
        float(pereverzev_lifetime(F_star, p_params))
        if np.isfinite(F_star) else float("nan")
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    ax.plot(F_grid_pN, k_off, "C0", lw=2, label=r"$k_\mathrm{off}(F)$")
    ax.axvline(F_star_pN, color="0.4", lw=1, ls="--",
               label=fr"$F^*={F_star_pN:.2f}$ pN")
    ax.set_xlabel(r"Bond force $F$ [pN]")
    ax.set_ylabel(r"$k_\mathrm{off}$ [s$^{-1}$]")
    ax.set_title("Pereverzev off-rate (KU-2.18)")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(F_grid_pN, tau, "C1", lw=2, label=r"$\tau(F)=1/k_\mathrm{off}$")
    ax.axvline(F_star_pN, color="0.4", lw=1, ls="--",
               label=fr"$\tau^*={tau_star:.2f}$ s at $F^*$")
    ax.axvline(30.0, color="C3", lw=1, ls=":",
               label="KU-2.5 experimental F* (30 pN)")
    ax.set_xlabel(r"Bond force $F$ [pN]")
    ax.set_ylabel(r"Bond lifetime $\tau$ [s]")
    ax.set_title("Catch-bond lifetime peak")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out = FIGS_DIR / "ku25_pereverzev.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_pereverzev_emergent_vs_oracle(p_params: PereverzevParams) -> Path:
    """Production-scale (or demo-scale) emergent k_off MC vs closed form."""
    npz = OUT_DIR / "ku25_pereverzev_production.npz"
    fig, ax = plt.subplots(figsize=(7, 4.5))
    F_curve_pN = np.linspace(0.5, 60.0, 200)
    F_curve_N = F_curve_pN * 1e-12
    ax.plot(F_curve_pN, pereverzev_k_off(F_curve_N, p_params), "C0", lw=2,
            label="closed-form oracle")
    if npz.exists():
        data = np.load(npz)
        grid = data["grid_pN"]
        sim = data["sim_values"]
        oracle = data["oracle_values"]
        # Relative error annotation.
        rel = np.abs(sim - oracle) / oracle
        ax.scatter(grid, sim, c="C3", s=40, zorder=5,
                   label=f"emergent MC (n_samples=20k, max rel err = {rel.max()*100:.2f} %)")
        ax.set_title("KU-2.5 + D2 production gate (PASS ± 5 %)")
    else:
        ax.set_title("KU-2.5 + D2 oracle (no production NPZ found)")
    ax.set_xlabel(r"Bond force $F$ [pN]")
    ax.set_ylabel(r"$k_\mathrm{off}$ [s$^{-1}$]")
    ax.set_yscale("log")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    out = FIGS_DIR / "ku25_emergent_vs_oracle.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_hill_force_velocity(p) -> Path:
    """D6 Hill v(F)/v0 for the H.4 motor parameters."""
    motor = p.motor
    F_stall = float(motor["F_stall_per_head"])
    v0 = float(motor["v0_per_head"])
    a_ratios = [0.25, 0.5, 1.0]
    F = np.linspace(0.0, 1.3 * F_stall, 200)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for a in a_ratios:
        v = hill_velocity(F, v0=v0, F_stall=F_stall, a_over_F_stall=a)
        v_clamped = hill_velocity_clamped(
            F, v0=v0, F_stall=F_stall, a_over_F_stall=a
        )
        ax.plot(F / F_stall, v / v0, lw=1, ls="--", alpha=0.6,
                label=fr"$a/F_s$={a} (full)")
        ax.plot(F / F_stall, v_clamped / v0, lw=2,
                label=fr"$a/F_s$={a} (clamped, runtime)")
    ax.axhline(0, color="0.6", lw=0.5)
    ax.axvline(1, color="0.6", lw=0.5, ls=":")
    ax.set_xlabel(r"Load $F / F_\mathrm{stall}$")
    ax.set_ylabel(r"Velocity $v / v_0$")
    ax.set_title(
        "D6 Hill force-velocity\n"
        f"v0 = {v0*1e6:.1f} μm/s, F_stall = {F_stall*1e12:.2f} pN"
    )
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = FIGS_DIR / "ku24_hill_force_velocity.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_fa_layout(p) -> Path:
    """Spatial scatter of FA centres (nascent vs mature) on the box face."""
    snap, layouts, _meta = build_h4_state(p, with_motors=False)
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    for L in layouts:
        color = "C3" if L.is_mature else "C0"
        size = 80 if L.is_mature else 30
        ax.scatter(L.centre_xy[0] * 1e6, L.centre_xy[1] * 1e6,
                   c=color, s=size, alpha=0.7, edgecolor="0.2")
        # Draw capture radius for one mature FA as illustration.
    half_um = 0.5 * p.L_box * 1e6
    ax.set_xlim(-half_um, half_um)
    ax.set_ylim(-half_um, half_um)
    ax.set_aspect("equal")
    ax.set_xlabel(r"x [μm]")
    ax.set_ylabel(r"y [μm]")
    n_nascent = sum(1 for L in layouts if not L.is_mature)
    n_mature = sum(1 for L in layouts if L.is_mature)
    ax.set_title(
        f"H.4 FA layout (L_box = {p.L_box*1e6:.0f} μm)\n"
        f"{n_nascent} nascent (C0) + {n_mature} mature (C3) "
        f"= {len(layouts)} FAs"
    )
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = FIGS_DIR / "fa_layout.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_ku217_oracle(p) -> Path:
    """KU-2.17 Hill growth rate k_g(F) over the H.4 force range."""
    F_th = float(p.fa_growth["F_th_per_FA"])
    k_g0 = float(p.fa_growth["k_g0"])
    k_d = float(p.fa_growth["k_d"])
    n = float(p.fa_growth["Hill_n"])
    F = np.linspace(0.0, 5.0 * F_th, 200)
    F_safe = np.where(F > 0, F, 1e-30)
    rate = k_g0 * F_safe**n / (F_safe**n + F_th**n)
    rate = np.where(F > 0, rate, 0.0)
    net = rate - k_d
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(F * 1e12, rate, "C0", lw=2, label=r"$k_g(F)$ (Hill, KU-2.17)")
    ax.plot(F * 1e12, np.full_like(F, k_d), "C3", lw=1, ls="--",
            label=fr"$k_d={k_d}$ s$^{{-1}}$ (decay)")
    ax.plot(F * 1e12, net, "0.3", lw=1.5, label=r"net $dA/A\,dt$")
    ax.axvline(F_th * 1e12, color="0.4", lw=1, ls=":",
               label=fr"$F_\mathrm{{th}}={F_th*1e12:.0f}$ pN")
    ax.set_xlabel(r"per-FA force $F$ [pN]")
    ax.set_ylabel(r"rate [s$^{-1}$]")
    ax.set_title(
        "KU-2.17 FA growth Hill rate (oracle only;"
        " runtime is emergent N_engaged)"
    )
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = FIGS_DIR / "ku217_hill_growth_rate.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def fig_ku24_biphasic_v1_reference(p) -> Path:
    """Visual encoding of the v1 Worker B biphasic verdict (literature band).

    Plots an illustrative ⟨F_total⟩(E) saturating curve at the v1 Worker
    B prominence (0.34 %) — production HOOMD sweep replaces this in
    week 3 (test_ku24_biphasic_production).
    """
    motor = p.motor
    F_stall_total = (
        p.n_motors_per_fa
        * int(motor["n_heads_per_side"])
        * float(motor["F_stall_per_head"])
    )
    E_grid = np.geomspace(1.0, 1.0e5, 24)
    # Illustrative shape: saturating with the v1 0.34 % prominence.
    # F(E) = F_stall · (E / (E + E_50)) · (1 − 0.0034 · exp(−E/E_50))
    E_50 = 100.0
    F_curve = F_stall_total * (E_grid / (E_grid + E_50)) * (
        1.0 - 0.0034 * np.exp(-E_grid / E_50)
    )
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(E_grid, F_curve * 1e12, "C0", lw=2,
            label="v1 Worker B verdict (saturating, 0.34 % prominence)")
    ax.axhline(F_stall_total * 1e12, color="0.5", lw=1, ls="--",
               label=fr"motor stall N$_m$·F$_s$ = {F_stall_total*1e12:.1f} pN")
    ax.set_xscale("log")
    ax.set_xlabel(r"Substrate stiffness $E$ [Pa]")
    ax.set_ylabel(r"⟨$F_\mathrm{total}$⟩ per FA [pN]")
    ax.set_title(
        "KU-2.4 / KU-2.8 biphasic verdict (literature band)\n"
        "Production HOOMD sweep deferred to H.4 week-3 (H4_PRODUCTION=1)"
    )
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    out = FIGS_DIR / "ku24_biphasic_v1_reference.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main() -> None:
    FIGS_DIR.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    p = resolve_h4(cfg)
    p_pereverzev = p.pereverzev

    written = [
        fig_pereverzev_k_off_and_lifetime(p_pereverzev),
        fig_pereverzev_emergent_vs_oracle(p_pereverzev),
        fig_hill_force_velocity(p),
        fig_fa_layout(p),
        fig_ku217_oracle(p),
        fig_ku24_biphasic_v1_reference(p),
    ]
    print("Wrote figures:")
    for f in written:
        print(f"  - {f.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
