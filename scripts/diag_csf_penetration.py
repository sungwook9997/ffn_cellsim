"""Diagnostic — CSF interior penetration & smoothing study (Option e, no
solver code modification).

Builds the same post-calibration spheroid the v11 pilot starts from, extracts
the grid mass field, then applies the colour seed → N-pass box smoothing →
central-difference ∇c pipeline **entirely in numpy** for several smoothing-pass
counts. Produces:

    - interface thickness (cells with 0.05 < c < 0.95) in dx and R₀ units
    - volume fraction of cells with |∇c| above several relative thresholds —
      a direct measurement of the spatial extent of the CSF body force
    - radial profile of c(r) and |∇c|(r) about the spheroid centre

The script does **not** modify solver state, does not change gates, and is
not part of the production codepath. Output: stdout table + JSON dump.

Run:
    python -m scripts.diag_csf_penetration configs/stage1a_pilot.yaml
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from acs.config import load_config
from acs.gpu import init_taichi
from acs.physics.mlsmpm import MLSMPMSolver, SolverConfig


def _solver_cfg_from_yaml(cfg: dict) -> SolverConfig:
    nd = cfg["nondim"]
    nm = cfg["numerics"]
    sim = cfg["simulation"]
    return SolverConfig(
        n_particles=int(sim["n_material_points"]),
        grid_n=int(nm["background_grid_resolution"]),
        domain_star=float(nd["domain_star"]),
        radius_star=float(nd["radius_star"]),
        dt_star=float(nd["dt_star"]),
        K_star=float(nd["K_star"]),
        mu_star=float(nd["mu_star"]),
        tau_star=float(nd["tau_star"]),
        capillary_number=float(nd["capillary_number"]),
        drag_xi_star=float(nd["drag_xi_star"]),
        density_star=float(nd["density_star"]),
        free_surface_threshold=float(nm["free_surface_density_threshold"]),
        seed=int(cfg["run"]["seed"]),
    )


def smooth_box_3x3x3(c: np.ndarray, n_passes: int) -> np.ndarray:
    """3³ box-average (boundary-aware: same edge-clipping the solver does).

    Each pass replaces each cell with the mean over its 3³ neighbourhood,
    skipping out-of-domain neighbours. Same algorithm as
    `MLSMPMSolver._smooth_color_pass`, but in numpy.
    """
    out = c.astype(np.float64).copy()
    if n_passes <= 0:
        return out
    n_g = c.shape[0]
    for _ in range(n_passes):
        s = np.zeros_like(out)
        cnt = np.zeros_like(out)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                for dk in (-1, 0, 1):
                    src = np.roll(np.roll(np.roll(out, -di, 0), -dj, 1), -dk, 2)
                    # Mask out-of-bounds rolls (np.roll is periodic; we want clip).
                    mask = np.ones((n_g, n_g, n_g), dtype=bool)
                    if di == -1:
                        mask[-1, :, :] = False
                    elif di == 1:
                        mask[0, :, :] = False
                    if dj == -1:
                        mask[:, -1, :] = False
                    elif dj == 1:
                        mask[:, 0, :] = False
                    if dk == -1:
                        mask[:, :, -1] = False
                    elif dk == 1:
                        mask[:, :, 0] = False
                    s += np.where(mask, src, 0.0)
                    cnt += mask.astype(np.float64)
        out = s / np.clip(cnt, 1.0, None)
    return out


def grad_central(c: np.ndarray, dx: float) -> np.ndarray:
    """Central-difference ∇c, matching `_color_gradient_from_smoothed`."""
    n_g = c.shape[0]
    g = np.zeros((n_g, n_g, n_g, 3), dtype=np.float64)
    inv_2dx = 1.0 / (2.0 * dx)
    g[1:-1, :, :, 0] = (c[2:, :, :] - c[:-2, :, :]) * inv_2dx
    g[:, 1:-1, :, 1] = (c[:, 2:, :] - c[:, :-2, :]) * inv_2dx
    g[:, :, 1:-1, 2] = (c[:, :, 2:] - c[:, :, :-2]) * inv_2dx
    return g


def radial_profile(field: np.ndarray, dx: float, centre: np.ndarray, n_bins: int = 40,
                   r_max: float | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin a 3D scalar field by radial distance from `centre`.

    Returns (r_centres, mean_per_bin, count_per_bin).
    """
    n_g = field.shape[0]
    ix = np.arange(n_g)
    Ix, Iy, Iz = np.meshgrid(ix, ix, ix, indexing="ij")
    # Cell centres at (i+0.5)*dx
    Xc = (Ix + 0.5) * dx
    Yc = (Iy + 0.5) * dx
    Zc = (Iz + 0.5) * dx
    rr = np.sqrt((Xc - centre[0]) ** 2 + (Yc - centre[1]) ** 2 + (Zc - centre[2]) ** 2)
    if r_max is None:
        r_max = float(rr.max())
    edges = np.linspace(0.0, r_max, n_bins + 1)
    centres_r = 0.5 * (edges[:-1] + edges[1:])
    flat_r = rr.ravel()
    flat_v = field.ravel()
    bin_idx = np.clip(np.digitize(flat_r, edges) - 1, 0, n_bins - 1)
    sums = np.bincount(bin_idx, weights=flat_v, minlength=n_bins)
    counts = np.bincount(bin_idx, minlength=n_bins).astype(np.float64)
    means = np.where(counts > 0, sums / np.clip(counts, 1.0, None), 0.0)
    return centres_r, means, counts


def main(cfg_path: Path) -> dict:
    cfg = load_config(cfg_path)
    init_taichi(
        backend=cfg["gpu"].get("backend", "auto"),
        device_memory_GB=float(cfg["gpu"].get("device_memory_GB", 4)),
    )
    solver_cfg = _solver_cfg_from_yaml(cfg)
    solver = MLSMPMSolver(solver_cfg)
    centre_v = np.full(3, solver_cfg.domain_star * 0.5, dtype=np.float32)
    solver.initialize_sphere(centre_v)
    calib = solver.calibrate_reference_state()

    # Re-scatter mass at the post-calibration configuration to obtain grid_m.
    # _scatter_mass_only first clears grid_m internally.
    solver._scatter_mass_only()
    grid_m = solver.grid_m.to_numpy().astype(np.float64)

    dx = solver_cfg.dx_star
    rho_bulk = solver_cfg.density_star
    R0 = solver_cfg.radius_star

    c_raw = grid_m / (rho_bulk * (dx ** 3))           # matches _seed_color_from_mass
    pass_counts = [0, 1, 2, 3, 4]
    grad_thresholds = [0.001, 0.01, 0.10, 0.25]        # relative to max |∇c|

    domain_voxels = grid_m.size
    centre = np.array([solver_cfg.domain_star * 0.5] * 3, dtype=np.float64)

    summaries: list[dict] = []
    radial_records: dict[str, dict] = {}

    for n_pass in pass_counts:
        c = smooth_box_3x3x3(c_raw, n_pass)
        g = grad_central(c, dx)
        gmag = np.sqrt((g ** 2).sum(axis=-1))
        gmax = float(gmag.max())

        # Interface thickness via cells with 0.05 < c < 0.95.
        interface_mask = (c > 0.05) & (c < 0.95)
        n_interface_cells = int(interface_mask.sum())
        interface_volume_frac = n_interface_cells / domain_voxels
        # Express as a "thickness" = volume / surface-area; surface area ~ 4πR0².
        # In dx units: total spheroid surface ≈ 4π R₀² / dx². Interface thickness
        # in dx ≈ n_interface_cells / (surface in dx² units).
        surface_area_dx2 = 4.0 * np.pi * (R0 / dx) ** 2
        thickness_in_dx = n_interface_cells / max(surface_area_dx2, 1.0)
        thickness_in_R0 = thickness_in_dx * (dx / R0)

        # |∇c| volume fractions above several relative thresholds.
        grad_fractions = {}
        for th in grad_thresholds:
            grad_fractions[f"frac_|grad_c|>{int(th*100)}%_max"] = float(
                (gmag > th * gmax).sum()
            ) / domain_voxels

        # Spheroid-volume-fraction context: how big is the spheroid in the box?
        # Spheroid volume / box volume = (4πR₀³/3)/L³.
        spheroid_volume_frac = (4.0 / 3.0 * np.pi * R0 ** 3) / solver_cfg.domain_star ** 3

        # Radial profiles.
        r_centres, c_radial, _ = radial_profile(c, dx, centre, n_bins=40,
                                                r_max=1.5 * R0)
        _, gmag_radial, _ = radial_profile(gmag, dx, centre, n_bins=40,
                                           r_max=1.5 * R0)

        summary = {
            "n_passes": n_pass,
            "max_|grad_c|": gmax,
            "interface_cells": n_interface_cells,
            "interface_volume_fraction": interface_volume_frac,
            "interface_thickness_in_dx": thickness_in_dx,
            "interface_thickness_in_R0": thickness_in_R0,
            "spheroid_volume_fraction_of_box": spheroid_volume_frac,
            **grad_fractions,
        }
        summaries.append(summary)
        radial_records[f"n_pass={n_pass}"] = {
            "r_over_R0": (r_centres / R0).tolist(),
            "c_radial": c_radial.tolist(),
            "grad_c_mag_radial": gmag_radial.tolist(),
        }

    # ASCII table.
    print(f"\n=== CSF interior-penetration diagnostic ===")
    print(f"config: {cfg_path}")
    print(f"grid_n = {solver_cfg.grid_n}, dx* = {dx:.5f}, R0* = {R0:.4f}, domain* = {solver_cfg.domain_star}")
    print(f"spheroid volume / box volume = {(4.0/3.0*np.pi*R0**3)/solver_cfg.domain_star**3:.4f}")
    print(f"calibration: <J>_well_resolved = {calib['J_mean_well_resolved']:.8f}, "
          f"n_resolved = {calib['n_well_resolved']}/{solver_cfg.n_particles}")
    print()
    header = (
        f"{'passes':>7}  {'max|∇c|':>10}  "
        f"{'interface':>10}  {'thick/dx':>9}  {'thick/R0':>9}  "
        f"{'>1% max':>9}  {'>10% max':>9}  {'>25% max':>9}"
    )
    print(header)
    print("-" * len(header))
    for s in summaries:
        print(
            f"{s['n_passes']:>7}  "
            f"{s['max_|grad_c|']:>10.3f}  "
            f"{s['interface_cells']:>10d}  "
            f"{s['interface_thickness_in_dx']:>9.2f}  "
            f"{s['interface_thickness_in_R0']:>9.4f}  "
            f"{s['frac_|grad_c|>1%_max']*100:>8.2f}%  "
            f"{s['frac_|grad_c|>10%_max']*100:>8.2f}%  "
            f"{s['frac_|grad_c|>25%_max']*100:>8.2f}%"
        )
    print()
    print("Interpretation:")
    print("  - thick/R0 close to (n_passes+1)·dx/R0 means smoothing stays surface-localised.")
    print("  - frac_|grad_c|>1%_max ≫ surface band (~6πR0²/L³ for thin shell)")
    print("    means the CSF body force has interior penetration.")
    print()

    # JSON dump.
    out = {
        "config": str(cfg_path),
        "grid_n": solver_cfg.grid_n,
        "dx_star": dx,
        "R0_star": R0,
        "domain_star": solver_cfg.domain_star,
        "spheroid_volume_fraction_of_box": (4.0 / 3.0 * np.pi * R0 ** 3) / solver_cfg.domain_star ** 3,
        "thin_shell_volume_fraction_estimate": 4.0 * np.pi * R0 ** 2 * dx / solver_cfg.domain_star ** 3,
        "calibration_J_well_resolved": calib["J_mean_well_resolved"],
        "summaries": summaries,
        "radial_profiles": radial_records,
    }
    out_path = Path("results/stage1a_pilot/csf_penetration_diag.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"JSON → {out_path}")
    return out


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python -m scripts.diag_csf_penetration <config.yaml>", file=sys.stderr)
        raise SystemExit(2)
    main(Path(sys.argv[1]))
