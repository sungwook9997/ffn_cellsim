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

    fig, ax = plt.subplots(figsize=(6.5, 4.5), dpi=120)
    ax.plot(s, C, "o-", color="C0", label="measured C(s)", markersize=5)
    # WLC reference at brief L_p = 17 μm
    s_grid = np.linspace(0, s.max(), 200)
    ax.plot(s_grid, np.exp(-s_grid * p.rest_length / 17e-6),
            "--", color="C3", alpha=0.7,
            label="WLC continuum @ L_p = 17 μm (brief)")
    # WLC at L_p_C1 (local fit)
    ax.plot(s_grid, np.exp(-s_grid * p.rest_length / L_p_C1),
            "--", color="C1", alpha=0.7,
            label=f"WLC @ L_p_C1 = {L_p_C1*1e6:.2f} μm (local)")
    # WLC at L_p_tail (fit window 1..10)
    if math.isfinite(fit.L_p_m):
        ax.plot(s_grid, np.exp(-s_grid * p.rest_length / fit.L_p_m),
                "--", color="C2", alpha=0.7,
                label=f"WLC @ L_p_tail = {fit.L_p_m*1e6:.2f} μm (s∈[1,10] fit)")
    ax.set_xlabel("s (bond separation)")
    ax.set_ylabel("C(s) = ⟨t̂_i · t̂_{i+s}⟩")
    ax.set_title("H.2 single filament — tangent correlation\n"
                 "(post HOOMD tag-gather fix)")
    ax.set_ylim(0.4, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.5)
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

    fig, ax = plt.subplots(figsize=(6.5, 4.5), dpi=120)
    # Histogram of measured angles (in degrees from π).
    delta_deg = np.degrees(np.pi - theta_all)
    ax.hist(delta_deg, bins=60, density=True, alpha=0.7,
            color="C0", edgecolor="black", label="measured angles")
    # 3D Boltzmann reference (in same coordinate).
    delta_grid_rad = np.linspace(0, max(delta_deg.max() / 180 * np.pi, 1.0), 400)
    theta_grid = np.pi - delta_grid_rad
    density_3d = boltzmann_angle_density_3d(
        theta_grid, bending_modulus=p.bending_modulus,
        rest_length=p.rest_length, kT=p.kT,
    )
    # Normalise to match histogram domain
    # Convert grid back to degrees-from-π
    delta_grid_deg = np.degrees(delta_grid_rad)
    # Integrate density over delta_grid_deg
    norm = np.trapezoid(density_3d, delta_grid_deg)
    if norm > 0:
        ax.plot(delta_grid_deg, density_3d / norm, "--", color="C3", linewidth=2,
                label="3D Boltzmann reference\np(θ) ∝ sin(θ)·exp(-κ(π-θ)²/(2ℓ_0 kT))")
    ax.set_xlabel("(π − θ) [degrees]")
    ax.set_ylabel("density")
    ax.set_title("H.2 angle distribution — measured vs 3D Boltzmann reference\n"
                 "(BAOAB freeze §Open #1 3D-corrected form)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    out = H2_FIGS / "fig_h2_angle_distribution.png"
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


if __name__ == "__main__":
    main()
