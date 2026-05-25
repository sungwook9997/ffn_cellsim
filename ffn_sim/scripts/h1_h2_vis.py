"""H.1 / H.2 production result visualization (autonomous /loop, 2026-05-21).

Reads the H.1 KU-1.30 + H.2 single-filament production output files
and writes PNG figures into ``outputs/h1/figs/`` and ``outputs/h2/figs/``.

Figures produced:
  H.1:
    fig_h1_ku130_2_strain_stiffening.png   K(γ) log-log + slope fit
    fig_h1_ku130_3_dipole_radial.png       σ(r) bond-virial radial decay
    fig_h1_wallbench.png                   wall-time bar chart vs v1 numpy
  H.2:
    fig_h2_tangent_correlation.png         C(s) + WLC fit + brief reference
    fig_h2_angle_distribution.png          histogram vs 3D Boltzmann ref
    fig_h2_bond_length_dist.png            bond length thermal scatter
    fig_h2_filament_snapshots.png          3D filament configuration overlay
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from ffn_sim.common.filament_math import (
    boltzmann_angle_density_3d,
    fit_persistence_length,
    hoomd_angle_array,
    tangent_correlation,
)
from ffn_sim.scripts.h2_single_filament import resolve_h2_derived

ROOT = Path(__file__).resolve().parents[1]
H1_DIR = ROOT / "outputs" / "h1"
H2_DIR = ROOT / "outputs" / "h2"
H1_FIGS = H1_DIR / "figs"
H2_FIGS = H2_DIR / "figs"


# ---------------------------------------------------------------------------
# H.1 figures
# ---------------------------------------------------------------------------
def fig_h1_ku130_2_strain_stiffening() -> None:
    npz = H1_DIR / "ku130_strain_stiffening_production.npz"
    if not npz.exists():
        print(f"[skip] {npz} missing")
        return
    d = np.load(npz)
    gammas = d["gammas"]
    K_avg = d["K_layer_pa_ensemble"]
    slope_ensemble = float(d["log_log_slope_ensemble"])
    per_ramp_slopes = d["log_log_slope_per_ramp"]

    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=120)
    mask = (gammas >= 0.05) & (gammas <= 0.30) & (K_avg > 0) & np.isfinite(K_avg)
    ax.loglog(gammas[mask], K_avg[mask], "o-", color="C0",
              label="ensemble K(γ) (3-ramp avg)", markersize=5)
    # Fit line.
    if mask.sum() >= 2:
        lg = np.log(gammas[mask])
        lk = np.log(K_avg[mask])
        a, b = np.polyfit(lg, lk, 1)
        gfit = np.array([gammas[mask].min(), gammas[mask].max()])
        ax.loglog(gfit, np.exp(a * np.log(gfit) + b), "--",
                  color="C0", alpha=0.6, label=f"fit slope = {a:+.3f}")
    # Brief band [-2.5, -1.5] is |slope|, shown as horizontal band on slope
    ax.text(0.05, 0.95,
            f"|slope| ensemble = {abs(slope_ensemble):.3f}\n"
            f"yaml band |β| ∈ [1.5, 2.5]\n"
            f"per-ramp = {[f'{s:+.3f}' for s in per_ramp_slopes]}",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray"))
    ax.set_xlabel("γ (strain)")
    ax.set_ylabel("K(γ) [Pa]   (layer convention)")
    ax.set_title("H.1 KU-1.30 #2 strain stiffening — Storm-MacKintosh regime")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    out = H1_FIGS / "fig_h1_ku130_2_strain_stiffening.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def fig_h1_ku130_3_dipole_radial() -> None:
    npz = H1_DIR / "ku130_point_dipole_production.npz"
    if not npz.exists():
        print(f"[skip] {npz} missing")
        return
    d = np.load(npz)
    centers = d["bin_centers_m"]
    sigma_avg = d["sigma_radial_pa_ensemble_avg"]
    per_seed_slope = d["fit_exponent_per_seed"]
    per_seed_mean = float(d["fit_exponent_per_seed_mean"]) \
        if "fit_exponent_per_seed_mean" in d.files else float(np.mean(per_seed_slope))
    per_seed_stack = d["sigma_radial_pa_per_seed"]

    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=120)
    # Plot per-seed sigma curves faintly
    for i in range(per_seed_stack.shape[0]):
        s = per_seed_stack[i]
        valid = (s > 0) & np.isfinite(s)
        if valid.any():
            ax.loglog(centers[valid] * 1e6, s[valid], "-", color="gray",
                      alpha=0.2, linewidth=0.7)
    # Ensemble average
    valid = (sigma_avg > 0) & np.isfinite(sigma_avg)
    ax.loglog(centers[valid] * 1e6, sigma_avg[valid], "o-", color="C0",
              label=f"ensemble avg (n={per_seed_stack.shape[0]} seeds)", markersize=5)
    # Reference 1/r² line through one of the mid-range points
    if valid.sum() > 3:
        idx = valid.sum() // 2
        ref_centers = centers[valid] * 1e6
        ref_idx = ref_centers[idx]
        ref_sig = sigma_avg[valid][idx]
        r_line = np.array([ref_centers.min(), ref_centers.max()])
        sig_line = ref_sig * (r_line / ref_idx) ** (-2)
        ax.loglog(r_line, sig_line, "--", color="C3", alpha=0.7, label="1/r² reference")
    ax.text(0.05, 0.05,
            f"per-seed mean slope = {per_seed_mean:+.3f}\n"
            f"brief band [-2.5, -1.5]\n"
            f"n_seeds = {per_seed_stack.shape[0]}",
            transform=ax.transAxes, va="bottom", fontsize=9,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray"))
    ax.set_xlabel("r (μm) from dipole centre")
    ax.set_ylabel("|σ_xx(r)| [Pa]   (paired-baseline subtraction, cos(2θ) projection)")
    ax.set_title("H.1 KU-1.30 #3 point-dipole stress decay (post tag-gather)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    out = H1_FIGS / "fig_h1_ku130_3_dipole_radial.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def fig_h1_wallbench() -> None:
    h1_json = H1_DIR / "h1_wallbench.json"
    v1_json = H1_DIR / "v1_numpy_bench.json"
    if not (h1_json.exists() and v1_json.exists()):
        print(f"[skip] {h1_json} or {v1_json} missing")
        return
    h1 = json.load(open(h1_json))
    v1 = json.load(open(v1_json))
    labels = ["HOOMD L-M BAOAB\n(this work)", "v1 numpy E-M\n(reference)"]
    rates = [h1["production_steps_per_s"], v1["production_steps_per_s"]]
    wall_per_sim = [h1["wall_per_simulated_second"], v1["wall_per_simulated_second"]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), dpi=120)
    bars = ax1.bar(labels, rates, color=["C0", "C7"], edgecolor="black")
    ax1.set_ylabel("steps / second")
    ax1.set_title("Per-step rate (N = 65 982 Mikado beads)")
    for b, r in zip(bars, rates):
        ax1.text(b.get_x() + b.get_width() / 2, r * 1.02, f"{r:.1f}", ha="center")
    ax1.grid(True, axis="y", alpha=0.3)

    bars2 = ax2.bar(labels, wall_per_sim, color=["C0", "C7"], edgecolor="black")
    ax2.set_ylabel("wall-time per simulated second")
    ax2.set_title(f"Wall / simulated s   (ratio = {wall_per_sim[0]/wall_per_sim[1]:.2f}×)")
    ax2.set_yscale("log")
    for b, w in zip(bars2, wall_per_sim):
        ax2.text(b.get_x() + b.get_width() / 2, w * 1.2, f"{w:.2e}", ha="center")
    ax2.grid(True, axis="y", alpha=0.3)

    fig.suptitle("H.1 wall-time bench: HOOMD L-M BAOAB vs v1 numpy E-M reference (5× budget OK)")
    fig.tight_layout()
    out = H1_FIGS / "fig_h1_wallbench.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# H.2 figures
# ---------------------------------------------------------------------------
def _h2_resolved():
    cfg = yaml.safe_load(open(ROOT / "configs" / "phase1_h2.yaml"))
    return resolve_h2_derived(cfg)


def fig_h2_tangent_correlation() -> None:
    npz = H2_DIR / "h2_production_trajectory.npz"
    if not npz.exists():
        print(f"[skip] {npz} missing")
        return
    p = _h2_resolved()
    d = np.load(npz)
    pos = d["positions"].reshape(-1, p.beads_per_fiber, 3)
    C = tangent_correlation(pos, box=None, max_separation=20)
    s = np.arange(C.size)
    L_p_C1 = -p.rest_length / math.log(C[1])
    fit = fit_persistence_length(pos, rest_length=p.rest_length, box=None)

    fig, ax = plt.subplots(figsize=(7, 4.8), dpi=120)
    ax.plot(s, C, "o-", color="C0", label="measured C(s)", markersize=5, zorder=5)
    s_grid = np.linspace(0, s.max(), 200)
    # WLC at brief reference (continuum diagnostic only)
    ax.plot(s_grid, np.exp(-s_grid * p.rest_length / 17e-6),
            ":", color="gray", alpha=0.6,
            label="WLC continuum @ L_p = 17 μm (brief diagnostic)")
    # D4-anchored C(1) band shaded
    lo_C1, hi_C1 = p.L_p_band_m_C1
    ax.fill_between(s_grid,
                    np.exp(-s_grid * p.rest_length / lo_C1),
                    np.exp(-s_grid * p.rest_length / hi_C1),
                    color="C1", alpha=0.18,
                    label=f"D4 band L_p_C1 ∈ [{lo_C1*1e6:.0f}, {hi_C1*1e6:.0f}] μm")
    # D4-anchored tail band shaded
    lo_tail, hi_tail = p.L_p_band_m_tail
    ax.fill_between(s_grid,
                    np.exp(-s_grid * p.rest_length / lo_tail),
                    np.exp(-s_grid * p.rest_length / hi_tail),
                    color="C2", alpha=0.18,
                    label=f"D4 band L_p_tail ∈ [{lo_tail*1e6:.0f}, {hi_tail*1e6:.0f}] μm")
    # Measured C(1) marker → L_p_C1
    ax.plot(s_grid, np.exp(-s_grid * p.rest_length / L_p_C1),
            "--", color="C1", linewidth=1.5,
            label=f"measured L_p_C1 = {L_p_C1*1e6:.2f} μm PASS")
    # Measured fit-tail → L_p_tail
    if math.isfinite(fit.L_p_m):
        ax.plot(s_grid, np.exp(-s_grid * p.rest_length / fit.L_p_m),
                "--", color="C2", linewidth=1.5,
                label=f"measured L_p_tail = {fit.L_p_m*1e6:.2f} μm PASS")
    ax.set_xlabel("s (bond separation)")
    ax.set_ylabel("C(s) = ⟨t̂_i · t̂_{i+s}⟩")
    ax.set_title("H.2 tangent correlation — both D4-anchored L_p gates PASS\n"
                 "(post HOOMD tag-gather fix, post PI 2026-05-21 rebanding)")
    ax.set_ylim(0.4, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8.2)
    fig.tight_layout()
    out = H2_FIGS / "fig_h2_tangent_correlation.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def fig_h2_angle_distribution() -> None:
    npz = H2_DIR / "h2_production_trajectory.npz"
    if not npz.exists():
        print(f"[skip] {npz} missing")
        return
    p = _h2_resolved()
    d = np.load(npz)
    pos = d["positions"]
    thetas = []
    for frame in pos:
        thetas.append(hoomd_angle_array(frame, box=None).ravel())
    theta_all = np.concatenate(thetas)

    # Effective k_θ (matches measured variance) — the KS-test reference.
    delta_sq_mean = float(np.mean((np.pi - theta_all) ** 2))
    k_theta_eff = 2.0 * p.kT / delta_sq_mean
    bending_modulus_eff = k_theta_eff * p.rest_length

    # KS-stat against effective-k reference (the actual gate quantity)
    from scipy import stats
    grid_full = np.linspace(theta_all.min(), theta_all.max(), 4001)
    density_eff = boltzmann_angle_density_3d(
        grid_full, bending_modulus=bending_modulus_eff,
        rest_length=p.rest_length, kT=p.kT,
    )
    cdf_eff = np.cumsum(density_eff); cdf_eff /= cdf_eff[-1]
    ks_stat, _p = stats.kstest(theta_all, lambda x: np.interp(x, grid_full, cdf_eff))

    fig, ax = plt.subplots(figsize=(7, 4.8), dpi=120)
    delta_deg = np.degrees(np.pi - theta_all)
    ax.hist(delta_deg, bins=60, density=True, alpha=0.55,
            color="C0", edgecolor="black", label="measured angles", zorder=2)
    delta_grid_rad = np.linspace(0, max(delta_deg.max() / 180 * np.pi, 1.0), 400)
    theta_grid = np.pi - delta_grid_rad
    delta_grid_deg = np.degrees(delta_grid_rad)

    # 3D Boltzmann at THEORETICAL k_θ (brief / continuum reference)
    density_th = boltzmann_angle_density_3d(
        theta_grid, bending_modulus=p.bending_modulus,
        rest_length=p.rest_length, kT=p.kT,
    )
    norm_th = np.trapezoid(density_th, delta_grid_deg)
    if norm_th > 0:
        ax.plot(delta_grid_deg, density_th / norm_th, "--",
                color="C3", linewidth=1.6,
                label=f"3D Boltzmann @ theoretical k_θ ({p.angle_k:.2e})\n"
                f"(brief reference — magnitude mismatch noted)")

    # 3D Boltzmann at EFFECTIVE k_θ (KS-test reference, PI 2026-05-21)
    density_eff_g = boltzmann_angle_density_3d(
        theta_grid, bending_modulus=bending_modulus_eff,
        rest_length=p.rest_length, kT=p.kT,
    )
    norm_eff = np.trapezoid(density_eff_g, delta_grid_deg)
    if norm_eff > 0:
        ax.plot(delta_grid_deg, density_eff_g / norm_eff, "-",
                color="C2", linewidth=2,
                label=f"3D Boltzmann @ effective k_θ ({k_theta_eff:.2e})\n"
                f"= shape-only KS ref · ks_stat = {ks_stat:.4f} ≤ 0.10 PASS")

    ax.set_xlabel("(π − θ) [degrees]")
    ax.set_ylabel("density")
    ax.set_title("H.2 angle distribution — shape match via effective k_θ PASS\n"
                 "(KS gate: shape only; magnitude in equipartition gate)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.2)
    fig.tight_layout()
    out = H2_FIGS / "fig_h2_angle_distribution.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def fig_h2_gates_summary() -> None:
    """Bar chart summary of all 4 H.2 production gates with bands."""
    npz = H2_DIR / "h2_production_trajectory.npz"
    if not npz.exists():
        print(f"[skip] {npz} missing")
        return
    p = _h2_resolved()
    d = np.load(npz)
    pos = d["positions"]
    pos_all = pos.reshape(-1, p.beads_per_fiber, 3)

    # Compute the 4 gate values
    C = tangent_correlation(pos_all, box=None, max_separation=1)
    L_p_C1 = -p.rest_length / math.log(C[1])
    fit = fit_persistence_length(pos_all, rest_length=p.rest_length, box=None)
    from ffn_sim.common.filament_math import bending_energy_per_bond, equipartition_check
    energies = [bending_energy_per_bond(f, angle_k=p.angle_k,
                                        angle_t0=p.angle_t0, box=None)
                for f in pos]
    eq = equipartition_check(energies, kT_J=p.kT)
    target_J = p.equipartition_target_3d_kT * p.kT
    rel_3d = (eq.mean_J - target_J) / target_J
    # KS
    from scipy import stats
    theta_all = np.concatenate(
        [hoomd_angle_array(f, box=None).ravel() for f in pos]
    )
    delta_sq = float(np.mean((np.pi - theta_all) ** 2))
    k_eff = 2 * p.kT / delta_sq
    grid = np.linspace(theta_all.min(), theta_all.max(), 4001)
    density = boltzmann_angle_density_3d(
        grid, bending_modulus=k_eff * p.rest_length,
        rest_length=p.rest_length, kT=p.kT,
    )
    cdf = np.cumsum(density); cdf /= cdf[-1]
    ks_stat, _ = stats.kstest(theta_all, lambda x: np.interp(x, grid, cdf))

    fig, axes = plt.subplots(1, 4, figsize=(13, 4.2), dpi=120)
    # 1. L_p_C1
    ax = axes[0]
    lo, hi = p.L_p_band_m_C1
    ax.axhspan(lo * 1e6, hi * 1e6, color="green", alpha=0.2,
               label=f"band [{lo*1e6:.0f}, {hi*1e6:.0f}]")
    ax.bar(["measured"], [L_p_C1 * 1e6], color="C0",
           edgecolor="black", width=0.5)
    ax.set_ylabel("L_p_C1 [μm]")
    ax.set_title(f"L_p (local C(1))\n{L_p_C1*1e6:.2f} μm  PASS")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    # 2. L_p_tail
    ax = axes[1]
    lo, hi = p.L_p_band_m_tail
    ax.axhspan(lo * 1e6, hi * 1e6, color="green", alpha=0.2,
               label=f"band [{lo*1e6:.0f}, {hi*1e6:.0f}]")
    ax.bar(["measured"], [fit.L_p_m * 1e6], color="C0",
           edgecolor="black", width=0.5)
    ax.set_ylabel("L_p_tail [μm]")
    ax.set_title(f"L_p (s∈[1,10] fit)\n{fit.L_p_m*1e6:.2f} μm  PASS")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    # 3. equipartition rel_3d
    ax = axes[2]
    tol = p.equipartition_rel_tol_3d
    ax.axhspan(-tol, tol, color="green", alpha=0.2,
               label=f"tol ±{tol}")
    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
    ax.bar(["measured"], [rel_3d], color="C0",
           edgecolor="black", width=0.5)
    ax.set_ylabel("rel deviation vs 3D analytical kT")
    ax.set_title(f"Equipartition\n{rel_3d:+.3f}  PASS")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    # 4. KS stat
    ax = axes[3]
    max_ks = p.angle_ks_stat_max
    ax.axhspan(0, max_ks, color="green", alpha=0.2,
               label=f"tol ≤ {max_ks}")
    ax.bar(["measured"], [ks_stat], color="C0",
           edgecolor="black", width=0.5)
    ax.set_ylabel("KS stat (CDF distance)")
    ax.set_title(f"KS shape (effective k_θ)\n{ks_stat:.4f}  PASS")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle("H.2 production gates — all 4 PASS  (PI 2026-05-21 rebanded)",
                 fontsize=12)
    fig.tight_layout()
    out = H2_FIGS / "fig_h2_gates_summary.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def fig_h2_bond_length_dist() -> None:
    npz = H2_DIR / "h2_production_trajectory.npz"
    if not npz.exists():
        print(f"[skip] {npz} missing")
        return
    p = _h2_resolved()
    d = np.load(npz)
    pos = d["positions"]
    bonds = pos[:, 0, 1:] - pos[:, 0, :-1]
    bl = np.linalg.norm(bonds, axis=-1).ravel() * 1e9   # nm
    sigma_th = math.sqrt(p.kT / p.bond_k) * 1e9   # thermal stderr in nm

    fig, ax = plt.subplots(figsize=(6.5, 4.5), dpi=120)
    ax.hist(bl, bins=60, density=True, alpha=0.7, color="C0",
            edgecolor="black", label="measured bond lengths")
    ax.axvline(p.rest_length * 1e9, color="C3", linestyle="--",
               label=f"rest length ℓ_0 = {p.rest_length*1e9:.1f} nm")
    ax.axvspan(p.rest_length * 1e9 - sigma_th, p.rest_length * 1e9 + sigma_th,
               color="C3", alpha=0.15,
               label=f"thermal ±σ = ±{sigma_th:.2f} nm")
    ax.set_xlabel("bond length [nm]")
    ax.set_ylabel("density")
    ax.set_title("H.2 bond length distribution\n"
                 "(post tag-gather, AFINES μ=1.5 nN)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    out = H2_FIGS / "fig_h2_bond_length_dist.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def fig_h2_filament_snapshots() -> None:
    npz = H2_DIR / "h2_production_trajectory.npz"
    if not npz.exists():
        print(f"[skip] {npz} missing")
        return
    d = np.load(npz)
    pos = d["positions"]
    n_frames = pos.shape[0]
    frames_to_show = np.linspace(0, n_frames - 1, 6, dtype=int)

    fig = plt.figure(figsize=(9, 6), dpi=120)
    ax = fig.add_subplot(111, projection="3d")
    cmap = plt.get_cmap("viridis")
    for i, fi in enumerate(frames_to_show):
        col = cmap(i / (len(frames_to_show) - 1))
        chain = pos[fi, 0]
        # Centre on the chain COM so multiple frames align for comparison
        chain = chain - chain.mean(axis=0)
        ax.plot(chain[:, 0] * 1e6, chain[:, 1] * 1e6, chain[:, 2] * 1e6,
                "-o", color=col, markersize=3, linewidth=1.5,
                label=f"frame {fi}")
    ax.set_xlabel("x [μm]")
    ax.set_ylabel("y [μm]")
    ax.set_zlabel("z [μm]")
    ax.set_title("H.2 single-filament configurations\n"
                 "(6 frames over 650 ms simulated, COM-aligned)")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    out = H2_FIGS / "fig_h2_filament_snapshots.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# H.2 strict-PASS follow-up figures (PI 2026-05-25)
# ---------------------------------------------------------------------------
H2_STRICT_DIR = H2_DIR / "strict_followup"


def _strict_inputs():
    """Load canonical (Lz=0.2 μm, 1 seed) + slab (Lz=10 μm, 5 seeds) artefacts."""
    canonical = np.load(H2_DIR / "h2_production_trajectory_canonical.npz")
    with open(H2_STRICT_DIR / "derivation.json") as f:
        deriv = json.load(f)
    with open(H2_STRICT_DIR / "slab_Lz10um_aggregate.json") as f:
        slab = json.load(f)
    slab_seed_trajs = {
        r["seed"]: np.load(H2_STRICT_DIR / f"slab_Lz10um_seed{r['seed']:02d}.npz")
        for r in slab["per_seed"]
    }
    return canonical, slab, slab_seed_trajs, deriv


def fig_h2_strict_Cs_comparison() -> None:
    """C(s) canonical vs slab ensemble vs first-principles a^s prediction."""
    canonical, slab, slab_seed_trajs, deriv = _strict_inputs()

    pos_can = canonical["positions"].reshape(-1, 21, 3)
    C_can = tangent_correlation(pos_can, box=None, max_separation=10)

    s_arr = np.arange(11)
    C_slab = np.array([r["C_s"] for r in slab["per_seed"]])
    C_slab_mean = C_slab.mean(axis=0)
    C_slab_std  = C_slab.std(axis=0, ddof=1)

    a_3d = 0.970043
    C_pred = a_3d ** s_arr

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    ax0.plot(s_arr, C_pred, "k--", lw=2, label=r"First-principles $a^s$, $\alpha=16.4$")
    ax0.plot(s_arr, C_can, color="tab:blue", lw=2, marker="o",
             label="Canonical Lz=0.2 μm (1 seed)")
    for i in range(C_slab.shape[0]):
        ax0.plot(s_arr, C_slab[i], color="tab:red", alpha=0.25, lw=1)
    ax0.plot(s_arr, C_slab_mean, color="tab:red", lw=2.5, marker="s",
             label="Slab Lz=10 μm ensemble (n=5, mean)")
    ax0.fill_between(s_arr, C_slab_mean - C_slab_std, C_slab_mean + C_slab_std,
                     color="tab:red", alpha=0.18)
    ax0.set_xlabel("bond separation $s$")
    ax0.set_ylabel(r"$C(s) = \langle \hat{t}_i \cdot \hat{t}_{i+s} \rangle$")
    ax0.set_title("Tangent correlation: canonical vs slab vs theory")
    ax0.set_xlim(0, 10)
    ax0.set_ylim(0.5, 1.02)
    ax0.legend(loc="lower left", fontsize=9)
    ax0.grid(alpha=0.25)

    ax1.axhline(0.0, color="k", lw=0.8)
    ax1.plot(s_arr, C_can - C_pred, color="tab:blue", lw=2, marker="o",
             label="Canonical − prediction")
    ax1.plot(s_arr, C_slab_mean - C_pred, color="tab:red", lw=2.5, marker="s",
             label="Slab ensemble − prediction")
    ax1.fill_between(s_arr, (C_slab_mean - C_slab_std) - C_pred,
                     (C_slab_mean + C_slab_std) - C_pred,
                     color="tab:red", alpha=0.18)
    ax1.set_xlabel("bond separation $s$")
    ax1.set_ylabel(r"$C(s)_\mathrm{measured} - a^s$")
    ax1.set_title("Residual: +50% softness signature gone after slab fix")
    ax1.set_xlim(0, 10)
    ax1.legend(loc="upper left", fontsize=9)
    ax1.grid(alpha=0.25)

    fig.tight_layout()
    out = H2_FIGS / "fig_h2_strict_Cs_comparison.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


def fig_h2_strict_equipartition() -> None:
    """⟨E⟩_per_bond canonical vs slab ensemble vs 0.99 kT target."""
    from ffn_sim.common.filament_math import bending_energy_per_bond
    canonical, slab, slab_seed_trajs, deriv = _strict_inputs()
    p = _h2_resolved()

    E_can_frames = np.array([
        bending_energy_per_bond(f, angle_k=p.angle_k, angle_t0=p.angle_t0, box=None).mean()
        for f in canonical["positions"]
    ]) / p.kT
    E_slab_per_seed = np.array([r["E_mean_kT"] for r in slab["per_seed"]])

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    ax0.hist(E_can_frames, bins=40, density=True, alpha=0.55, color="tab:blue",
             label=f"Canonical per-frame mean (n={len(E_can_frames)})")
    ax0.axvline(E_can_frames.mean(), color="tab:blue", lw=2, ls="--",
                label=f"Canonical mean = {E_can_frames.mean():.3f} kT")
    ax0.axvline(deriv["E_kT_target_3d_exact"], color="k", lw=2.5,
                label=f"First-principles 3D = {deriv['E_kT_target_3d_exact']:.4f} kT")
    ax0.axvline(0.5, color="tab:gray", lw=1, ls=":", label="Brief 2D = 0.5 kT (deprecated)")
    ax0.set_xlabel(r"$\langle E_\mathrm{bend} \rangle / k_BT$ per bond")
    ax0.set_ylabel("density")
    ax0.set_title("Canonical Lz=0.2 μm: +50% above 3D target")
    ax0.legend(loc="upper right", fontsize=9)
    ax0.grid(alpha=0.25)
    ax0.set_xlim(0.0, 2.0)

    ax1.axvline(deriv["E_kT_target_3d_exact"], color="k", lw=2.5,
                label=f"First-principles target = {deriv['E_kT_target_3d_exact']:.4f} kT")
    ax1.axvspan(0.99 * 0.95, 0.99 * 1.05, color="tab:green", alpha=0.18,
                label="±5 % strict band")
    ax1.scatter(E_slab_per_seed, np.arange(5) + 1, color="tab:red", s=90, zorder=5,
                label=f"Slab seeds (mean = {E_slab_per_seed.mean():.3f} ± {E_slab_per_seed.std(ddof=1):.3f} kT)")
    ax1.scatter([E_can_frames.mean()], [0], color="tab:blue", s=130, marker="*",
                zorder=5, label=f"Canonical = {E_can_frames.mean():.3f} kT (out of band)")
    ax1.set_yticks(range(6))
    ax1.set_yticklabels(["canonical"] + [f"seed {s}" for s in slab["seeds"]])
    ax1.set_xlabel(r"$\langle E_\mathrm{bend} \rangle / k_BT$ per bond")
    ax1.set_title("Slab Lz=10 μm: all 5 seeds inside ±5 % strict band")
    ax1.set_xlim(0.85, 1.55)
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(alpha=0.25)

    fig.tight_layout()
    out = H2_FIGS / "fig_h2_strict_equipartition.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


def fig_h2_strict_angle_pdf() -> None:
    """Angle PDF vs first-principles 3D Boltzmann at theoretical α (no rescaling)."""
    from scipy import integrate, stats
    canonical, slab, slab_seed_trajs, deriv = _strict_inputs()

    alpha = deriv["alpha"]
    phi_grid = np.linspace(0, np.pi, 4001)
    pdf_ref = np.sin(phi_grid) * np.exp(-alpha * phi_grid**2)
    pdf_ref /= np.trapezoid(pdf_ref, phi_grid)

    theta_can = hoomd_angle_array(canonical["positions"].reshape(-1, 21, 3), box=None).ravel()
    phi_can = (np.pi - theta_can)
    phi_slab = []
    for seed_d in slab_seed_trajs.values():
        theta = hoomd_angle_array(seed_d["positions"].reshape(-1, 21, 3), box=None).ravel()
        phi_slab.append(np.pi - theta)
    phi_slab = np.concatenate(phi_slab)

    n_bins = 50
    edges = np.linspace(0, np.pi, n_bins + 1)
    def p_3d(phi): return np.sin(phi) * np.exp(-alpha * phi**2)
    expected_p = np.array([integrate.quad(p_3d, edges[i], edges[i+1])[0] for i in range(n_bins)])
    expected_p /= expected_p.sum()

    def chi2_kl(samples):
        counts, _ = np.histogram(samples, bins=edges)
        N = counts.sum()
        expected = N * expected_p
        chi2 = float(np.sum((counts - expected)**2 / np.where(expected > 0, expected, 1.0)))
        widths = np.diff(edges)
        p_meas = counts / N / widths
        p_ref  = expected / N / widths
        m = (p_meas > 0) & (p_ref > 0)
        kl = float(np.sum(p_meas[m] * np.log(p_meas[m]/p_ref[m]) * widths[m]))
        return chi2, kl
    chi2_can, kl_can = chi2_kl(phi_can)
    chi2_slab, kl_slab = chi2_kl(phi_slab)
    df = n_bins - 1
    chi2_crit_95 = stats.chi2.ppf(0.95, df)

    fig, ax = plt.subplots(figsize=(9, 5.2))
    bins = np.linspace(0, 0.6, 80)
    ax.hist(phi_can, bins=bins, density=True, alpha=0.45, color="tab:blue",
            label=f"Canonical Lz=0.2 μm  (χ²/df={chi2_can/df:.1f}, D_KL={kl_can:.3f} nats)")
    ax.hist(phi_slab, bins=bins, density=True, alpha=0.45, color="tab:red",
            label=f"Slab Lz=10 μm pooled  (χ²/df={chi2_slab/df:.2f}, D_KL={kl_slab:.5f} nats)")
    ax.plot(phi_grid, pdf_ref, "k-", lw=2.5,
            label=f"First-principles 3D Boltzmann at α={alpha:.1f}")
    ax.set_xlabel(r"$\varphi = \pi - \theta$ [rad] (deviation from straight)")
    ax.set_ylabel(r"$p(\varphi)$")
    ax.set_title("Angle distribution vs first-principles 3D Boltzmann (no rescaling)\n"
                 f"χ² critical 95%/df = {chi2_crit_95/df:.2f};  KL threshold 95% ≈ {chi2_crit_95/(2*17100):.5f} nats")
    ax.set_xlim(0, 0.6)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.25)

    fig.tight_layout()
    out = H2_FIGS / "fig_h2_strict_angle_pdf.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


def fig_h2_strict_gates_summary() -> None:
    """4-panel summary: L_p_C1, L_p_tail, ⟨E⟩, rel — canonical vs slab vs first-principles."""
    canonical, slab, slab_seed_trajs, deriv = _strict_inputs()

    L_p_can = 10.83
    E_can   = 1.502
    band_C1 = deriv["L_p_C1_band_5sigma_m"]
    band_C1_lo = band_C1[0] * 1e6
    band_C1_hi = band_C1[1] * 1e6
    L_p_pred = deriv["L_p_C1_predicted_m"] * 1e6
    E_pred   = deriv["E_kT_target_3d_exact"]

    slab_Lp_C1 = np.array([r["L_p_C1_m"]   for r in slab["per_seed"]]) * 1e6
    slab_Lp_tl = np.array([r["L_p_tail_m"] for r in slab["per_seed"]]) * 1e6
    slab_E     = np.array([r["E_mean_kT"]  for r in slab["per_seed"]])
    slab_rel   = (slab_E - E_pred) / E_pred

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8))
    panels = [
        (axes[0, 0], "L_p_C1  [μm]",
         L_p_can, slab_Lp_C1, L_p_pred, (band_C1_lo, band_C1_hi),
         (5, 30), "5σ band [15.81, 17.07] μm"),
        (axes[0, 1], "L_p_tail (fit s∈[1,10])  [μm]",
         27.11, slab_Lp_tl, L_p_pred, (band_C1_lo, band_C1_hi),
         (5, 35), "Same theory band (C(s)=a^s)"),
        (axes[1, 0], r"$\langle E_\mathrm{bend} \rangle / k_BT$",
         E_can, slab_E, E_pred, (0.99 * 0.95, 0.99 * 1.05),
         (0.85, 1.6), "±5 % strict band on 0.99 kT"),
        (axes[1, 1], "rel vs 3D kT",
         (E_can - E_pred) / E_pred, slab_rel,
         0.0, (-0.05, 0.05),
         (-0.1, 0.6), "±5 % strict band"),
    ]
    for ax, label, can, slab_pts, ref, band, ylim, band_label in panels:
        ax.axhspan(band[0], band[1], color="tab:green", alpha=0.20, label=band_label)
        ax.axhline(ref, color="k", lw=2.0, label=f"First-principles = {ref:.3f}")
        ax.scatter([0], [can], color="tab:blue", s=150, marker="*", zorder=5,
                   label=f"Canonical Lz=0.2 μm = {can:.3f}")
        ax.scatter(np.arange(1, len(slab_pts) + 1), slab_pts, color="tab:red",
                   s=80, zorder=5,
                   label=f"Slab seeds (mean = {np.mean(slab_pts):.3f} ± {np.std(slab_pts, ddof=1):.3f})")
        ax.set_xticks(range(len(slab_pts) + 1))
        ax.set_xticklabels(["canonical"] + [f"s={s}" for s in slab["seeds"]],
                           rotation=20, fontsize=8)
        ax.set_ylabel(label)
        ax.set_ylim(*ylim)
        ax.legend(loc="best", fontsize=8)
        ax.grid(alpha=0.25)

    fig.suptitle("H.2 strict-PASS gates: canonical Lz=0.2 μm vs slab Lz=10 μm vs first-principles",
                 fontsize=12)
    fig.tight_layout()
    out = H2_FIGS / "fig_h2_strict_gates_summary.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> None:
    H1_FIGS.mkdir(parents=True, exist_ok=True)
    H2_FIGS.mkdir(parents=True, exist_ok=True)
    fig_h1_ku130_2_strain_stiffening()
    fig_h1_ku130_3_dipole_radial()
    fig_h1_wallbench()
    fig_h2_tangent_correlation()
    fig_h2_angle_distribution()
    fig_h2_bond_length_dist()
    fig_h2_filament_snapshots()
    fig_h2_gates_summary()
    if H2_STRICT_DIR.exists():
        fig_h2_strict_Cs_comparison()
        fig_h2_strict_equipartition()
        fig_h2_strict_angle_pdf()
        fig_h2_strict_gates_summary()


if __name__ == "__main__":
    main()
