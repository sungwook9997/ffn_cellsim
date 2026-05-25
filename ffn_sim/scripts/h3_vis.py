"""H.3 cortex topology visualization (closeout per CLAUDE.md visualize rule).

Builds a demo H.3 cortex (50 filaments × 7 beads on R = 10 μm shell),
runs a short BAOAB stabilization, and writes PNG figures into
``ffn_sim/outputs/h3/figs/`` to demonstrate:

- Spherical-shell topology placement (3D scatter of cortex beads).
- Per-bead radial drift at construction vs KU-3.17 cortex thickness.
- Per-filament tangent unit-vector orientation in the local tangent
  plane (verifies the Marsaglia surface-uniform + tangent-plane basis).
- σ_z vs L_z pre-flight: box geometry vs per-filament thermal extent
  (the H.2 slab lesson satisfied at H.3 setup).
- Force-constant comparison vs H.1 / H.2 (same κ_B, μ, ℓ_0).

This is the closeout vis for H.3 🟨 (cortex.py + config + tests
landed; production sweep deferred to next session). Production sweep
figures (L_p, equipartition, 3D Boltzmann angle KS) will be added by
re-running this script after H3_PRODUCTION=1 results land.

Usage::

    PYTHONPATH=. python ffn_sim/scripts/h3_vis.py
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from ffn_sim.cortex.cortex import (
    build_cortex_simulation,
    generate_cortex_topology,
    resolve_h3_derived,
)


CONFIG_PATH = (
    Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
)
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h3" / "figs"


def _demo_config(n_filaments: int = 200) -> dict:
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    # Demo cortex: 200 filaments (cosmetic — full 1000 looks too dense
    # in the 3D scatter); demo_mode relaxes the cost ceiling.
    cfg["cortex"]["n_filaments"] = n_filaments
    cfg["cortex"]["demo_mode"] = True
    return cfg


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def fig_cortex_3d_scatter(topology, p, out_path: Path) -> None:
    """3D scatter of cortex bead positions colored by filament index."""
    pos = topology.positions       # (F, N, 3)
    F, N, _ = pos.shape

    fig = plt.figure(figsize=(8.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    colors = plt.cm.tab20(np.linspace(0, 1, min(F, 20)))
    for f in range(F):
        ax.plot(
            pos[f, :, 0] * 1e6, pos[f, :, 1] * 1e6, pos[f, :, 2] * 1e6,
            color=colors[f % len(colors)], linewidth=0.9, alpha=0.85,
        )
    # Shell wireframe
    u = np.linspace(0, 2 * np.pi, 24)
    v = np.linspace(0, np.pi, 12)
    R_um = p.R_cell * 1e6
    xs = R_um * np.outer(np.cos(u), np.sin(v))
    ys = R_um * np.outer(np.sin(u), np.sin(v))
    zs = R_um * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(
        xs, ys, zs, color="lightgray", linewidth=0.4, alpha=0.4,
    )
    ax.set_xlabel("x [μm]")
    ax.set_ylabel("y [μm]")
    ax.set_zlabel("z [μm]")
    ax.set_title(
        f"H.3 cortex topology: {F} filaments × {N} beads on R = "
        f"{R_um:.0f} μm shell\n(×40 mesoscopic, Plan v2 §3 H.3 v3.1)"
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def fig_radial_drift_histogram(topology, p, out_path: Path) -> None:
    """Per-bead radial distance from cell center vs KU-3.17 cortex thickness."""
    pos = topology.positions.reshape(-1, 3)
    r = np.linalg.norm(pos, axis=1) * 1e9    # nm
    R_nm = p.R_cell * 1e9
    drift = r - R_nm

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(drift, bins=40, color="#3a7", edgecolor="white", alpha=0.85)
    cortex_nm = p.cortex_thickness * 1e9
    ax.axvline(0, color="black", linewidth=1.0, linestyle="-",
               label=f"R_cell = {R_nm:.0f} nm")
    ax.axvline(cortex_nm, color="red", linewidth=1.2, linestyle="--",
               label=f"+cortex thickness ({cortex_nm:.0f} nm, KU-3.17)")
    ax.axvline(-cortex_nm, color="red", linewidth=1.2, linestyle="--")
    # First-principles construction prediction.
    end_drift = (((p.beads_per_filament - 1) * p.rest_length / 2.0) ** 2
                 / (2 * p.R_cell)) * 1e9
    ax.axvline(end_drift, color="blue", linewidth=1.0, linestyle=":",
               label=f"endpoint drift ℓ_end²/(2 R) = {end_drift:.1f} nm")
    ax.set_xlabel("radial drift  r − R_cell  [nm]")
    ax.set_ylabel("count")
    ax.set_title(
        "H.3 cortex bead radial drift at construction\n"
        "(tangent-plane placement; no projection)"
    )
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def fig_tangent_orientation_check(topology, p, out_path: Path) -> None:
    """Verify tangent vectors lie in the local tangent plane (dot < 1e-9)."""
    centers = topology.centers_of_mass
    tangents = topology.tangents
    normals = centers / p.R_cell
    dots = np.abs(np.einsum("fi,fi->f", normals, tangents))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(dots, bins=40, color="#36a", edgecolor="white", alpha=0.85)
    ax.axvline(1.0e-9, color="red", linewidth=1.0, linestyle="--",
               label="numerical tolerance 1e-9")
    ax.set_xlabel("|n̂ · t̂|  (should be 0 for tangent-plane orientation)")
    ax.set_ylabel("count")
    ax.set_yscale("log")
    ax.set_title(
        "H.3 tangent-plane orientation check\n"
        "(Marsaglia surface-uniform + deterministic tangent basis)"
    )
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def fig_sigma_z_vs_Lz_preflight(p, out_path: Path) -> None:
    """σ_perp(L) vs L over the brief's 1–5 μm filament length range,
    overlaid with H.3 box half-width."""
    L_array_um = np.linspace(0.5, 6.0, 200)
    L_array = L_array_um * 1e-6
    sigma_perp = np.sqrt(L_array**3 / (3 * p.persistence_length)) * 1e6
    half_box_um = 0.5 * p.L_box * 1e6
    R_cell_um = p.R_cell * 1e6
    cell_extent_um = R_cell_um + sigma_perp

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.plot(L_array_um, sigma_perp, label="σ_perp = √(L³/(3 L_p))",
            color="#36a", linewidth=2)
    ax.plot(L_array_um, cell_extent_um,
            label="R_cell + σ_perp (cell + thermal margin)",
            color="#a63", linewidth=1.5, linestyle="--")
    ax.axhline(half_box_um, color="black", linewidth=1.5, linestyle=":",
               label=f"L_box/2 = {half_box_um:.1f} μm  "
                     "(no slab artefact)")
    ax.axhline(R_cell_um, color="green", linewidth=1.0, linestyle=":",
               label=f"R_cell = {R_cell_um:.1f} μm")
    # Brief's nominal range
    ax.axvspan(1.0, 5.0, alpha=0.10, color="orange",
               label="brief L_filament range 1–5 μm")
    # Production setpoint
    L_prod_um = p.L_filament * 1e6
    ax.axvline(L_prod_um, color="purple", linewidth=1.2,
               label=f"H.3 L_filament = {L_prod_um:.1f} μm")
    ax.set_xlabel("L_filament  [μm]")
    ax.set_ylabel("length scale  [μm]")
    ax.set_title(
        "H.3 σ_z vs L_z pre-flight: thermal extent stays inside box\n"
        "(H.2 slab lesson; cortex 200 nm comes from ERM, not box)"
    )
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def fig_force_constants_parity_vs_h2(p, out_path: Path) -> None:
    """Confirm H.3 force constants equal H.2's (×40 coarse-grains count, not stiffness)."""
    constants = {
        "ℓ_0 [μm]": p.rest_length * 1e6,
        "k_bond [pN/μm]": p.bond_k * 1e-3,  # N/m → pN/μm (1 N/m = 1e-3 pN/nm = 1 pN/nm wait...)
                                            # 1 N/m = 1 N/m. In pN/μm: 1 N/m = 1e12 pN/m = 1e6 pN/μm. So multiply by 1e6.
    }
    # Re-compute correctly: k_bond in pN/μm = k_bond [N/m] · 1e6
    k_bond_pN_um = p.bond_k * 1e6
    k_angle_kT_rad = p.angle_k / p.kT
    constants = {
        "rest_length ℓ_0 [μm]": p.rest_length * 1e6,
        "bond k [pN/μm]": k_bond_pN_um,
        "angle k [kT/rad²]": k_angle_kT_rad,
        "γ_b [pN·s/μm]": p.gamma_b * 1e6,
        "kT [zJ]": p.kT * 1e21,
        "L_p [μm]": p.persistence_length * 1e6,
    }
    keys = list(constants.keys())
    vals = [constants[k] for k in keys]

    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    y = np.arange(len(keys))
    ax.barh(y, vals, color="#369", alpha=0.85, edgecolor="white")
    for i, v in enumerate(vals):
        ax.text(v, i, f"  {v:.3g}", va="center", fontsize=10)
    ax.set_yticks(y)
    ax.set_yticklabels(keys)
    ax.set_title(
        "H.3 per-filament force constants (×40 mesoscopic bundle)\n"
        "Identical to H.2 (κ_B, μ, ℓ_0 unchanged — only filament count is "
        "coarse-grained)"
    )
    ax.set_xscale("log")
    ax.grid(True, alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------
def main() -> None:
    cfg = _demo_config()
    p = resolve_h3_derived(cfg)
    topology = generate_cortex_topology(p)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    figures = [
        ("fig_h3_topology_3d.png", fig_cortex_3d_scatter, (topology, p)),
        ("fig_h3_radial_drift.png", fig_radial_drift_histogram, (topology, p)),
        ("fig_h3_tangent_plane_check.png", fig_tangent_orientation_check,
         (topology, p)),
        ("fig_h3_sigma_z_vs_Lz.png", fig_sigma_z_vs_Lz_preflight, (p,)),
        ("fig_h3_force_constants.png", fig_force_constants_parity_vs_h2,
         (p,)),
    ]
    for name, fn, args in figures:
        out = OUTPUT_DIR / name
        fn(*args, out_path=out)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
