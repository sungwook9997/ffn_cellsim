"""H.7 manifold lamellipodium REGION MASKS — driver + viz (PI direction (b)).

The next increment of the (b) spatial-substrate wiring after the FA-traction
field (``h7_manifold_traction.py``): it carves the lamellipodium regions
(``basal_ring`` isotropic-spreading collar; ``polarized_patch`` migrating
leading-edge cap) onto the surface manifold as **patch masks + local frames**
(``cortex/manifold_regions.py``), and demonstrates that the masks co-register
with the H.5/H.7 lamellipodium WAVE *particle* placement (every WAVE bead's home
patch lies inside its mask — the wiring proof).

Geometry only: no HOOMD, no cell build, no forces — the manifold is the shared
spatial substrate, the regions are pure geometric selections, the WAVE positions
come from the lamellipodium layout generators. The co-registration check is a
GATE (exits nonzero on failure).

Outputs:
  * ``ffn_sim/outputs/h7/figs/h7_manifold_regions.png`` — south-cap Lambert
    projection of the two region masks + per-patch in-plane growth directions,
    with the lamellipodium WAVE home patches overlaid.
  * ``ffn_sim/outputs/h7/figs/h7_manifold_regions.json`` — scalar summary
    (area fractions, patch counts, co-registration gate, resolution sweep).

Usage:
    python -m ffn_sim.scripts.h7_manifold_regions
    python -m ffn_sim.scripts.h7_manifold_regions --subdivisions 4 --R-cell 7.5e-6
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import yaml

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ffn_sim.cortex.manifold_regions import (  # noqa: E402
    basal_ring_region,
    polarized_patch_region,
)
from ffn_sim.common.surface_manifold import SurfaceManifold  # noqa: E402

_OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"
_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "phase1_h5.yaml"


def _resolve_h5(R_cell: float, n_WAVE: int):
    """Resolve H.5 lamellipodium params for the MCF7 single-cell box."""
    from ffn_sim.cell.lamellipodium import resolve_h5_lamellipodium

    with open(_CONFIG) as f:
        cfg = yaml.safe_load(f)
    L_box = 3.0 * R_cell
    cfg["lamellipodium"]["n_WAVE"] = n_WAVE
    cfg["lamellipodium"]["Y_max"] = 0.45 * L_box
    return resolve_h5_lamellipodium(cfg, L_box=L_box, dt=13.0e-9)


def _lambert_azimuthal_south(centroids: np.ndarray) -> np.ndarray:
    """South-cap Lambert azimuthal equal-area projection (mirrors traction viz).

    The basal contact sits on the south cap; equal-area so patch areas are not
    distorted. North-hemisphere patches map to NaN (the basal regions never use
    them).
    """
    rho = np.linalg.norm(centroids, axis=1)
    rho = np.where(rho > 0.0, rho, 1.0)
    zc = centroids[:, 2] / rho  # cos(colatitude from +z); south pole zc = −1
    denom = 1.0 + zc  # project from the NORTH pole: k = sqrt(2/(1+zc))
    denom = np.where(denom > 1.0e-12, denom, np.nan)
    k = np.sqrt(2.0 / denom)
    X = k * (centroids[:, 0] / rho)
    Y = k * (centroids[:, 1] / rho)
    return np.stack([X, Y], axis=1)


def _make_figure(
    manifold: SurfaceManifold,
    basal,
    polar,
    *,
    wave_basal_home: np.ndarray,
    wave_polar_home: np.ndarray,
    R_cell: float,
    out_png: Path,
) -> None:
    proj = _lambert_azimuthal_south(manifold.tri_centroids)
    south = manifold.tri_centroids[:, 2] < 0.0  # only the basal hemisphere

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 6.4), constrained_layout=True)

    for ax, reg, home, title, cmap_c in (
        (axes[0], basal, wave_basal_home, "basal_ring (isotropic spreading)", "tab:blue"),
        (axes[1], polar, wave_polar_home, "polarized_patch (migrating front)", "tab:red"),
    ):
        # All basal-hemisphere patches (context).
        ok = south & np.isfinite(proj[:, 0])
        ax.scatter(proj[ok, 0], proj[ok, 1], s=6, c="0.85", label="patch (other)")
        # Region patches.
        rmask = reg.mask & ok
        ax.scatter(
            proj[rmask, 0], proj[rmask, 1], s=18, c=cmap_c,
            label=f"{reg.name} ({reg.patch_ids.size} patches)",
        )
        # In-plane growth directions (project the tangent vector tip into 2D).
        tip = manifold.tri_centroids[rmask] + 0.18 * R_cell * reg.in_plane_dir[
            np.isin(reg.patch_ids, np.flatnonzero(rmask))
        ]
        tip_proj = _lambert_azimuthal_south(tip)
        base_proj = proj[rmask]
        good = np.isfinite(tip_proj[:, 0]) & np.isfinite(base_proj[:, 0])
        ax.quiver(
            base_proj[good, 0], base_proj[good, 1],
            (tip_proj[good, 0] - base_proj[good, 0]),
            (tip_proj[good, 1] - base_proj[good, 1]),
            angles="xy", scale_units="xy", scale=1.0, width=0.003,
            color="0.25", alpha=0.7,
        )
        # WAVE home patches (co-registration overlay).
        if home.size:
            hp = proj[home]
            hgood = np.isfinite(hp[:, 0])
            ax.scatter(
                hp[hgood, 0], hp[hgood, 1], s=26, facecolors="none",
                edgecolors="k", linewidths=0.9, label="WAVE home patch",
            )
        ax.set_title(title)
        ax.set_xlabel("Lambert X (south cap)")
        ax.set_ylabel("Lambert Y")
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        f"H.7 manifold lamellipodium region masks  (R_cell={R_cell*1e6:.1f} µm, "
        f"{manifold.n_tri} patches)  — arrows = per-patch in-plane growth dir",
        fontsize=11,
    )
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def _resolution_sweep(R_cell: float) -> dict:
    """Area-fraction of each region across the canonical {320,1280,5120} grid."""
    sweep = {"subdivisions": [2, 3, 4], "basal_ring": [], "polarized_patch": []}
    for s in sweep["subdivisions"]:
        m = SurfaceManifold.icosphere(subdivisions=s, radius=R_cell)
        sweep["basal_ring"].append(
            basal_ring_region(
                m, R_cell=R_cell, rest_length=0.5e-6, collar_half_angle=0.25,
            ).area_fraction
        )
        sweep["polarized_patch"].append(
            polarized_patch_region(
                m, R_cell=R_cell, half_angle_azimuth=math.pi / 6.0,
            ).area_fraction
        )
    return sweep


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--R-cell", type=float, default=7.5e-6, help="cortex radius [m]")
    ap.add_argument("--subdivisions", type=int, default=3, help="icosphere level")
    ap.add_argument("--n-wave", type=int, default=120, help="WAVE count for overlay")
    args = ap.parse_args()

    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    R_cell = float(args.R_cell)
    manifold = SurfaceManifold.icosphere(
        subdivisions=int(args.subdivisions), radius=R_cell
    )
    p = _resolve_h5(R_cell, int(args.n_wave))

    # --- basal_ring region + WAVE co-registration ---
    from ffn_sim.cell.lamellipodium_basal_ring import (
        generate_basal_ring_lamellipodium_layout,
    )

    br_layout = generate_basal_ring_lamellipodium_layout(
        p, wave_tag_start=0, R_cell=R_cell,
    )
    br_collar = 3.0 * manifold.mean_edge_length / R_cell
    basal = basal_ring_region(
        manifold, R_cell=R_cell, rest_length=p.rest_length,
        cap_depth=br_layout.geometry.cap_depth, collar_half_angle=br_collar,
    )
    br_home = manifold.nearest_patch(br_layout.layout.wave_positions)
    br_coreg = bool(np.all(basal.mask[br_home])) if br_home.size else True

    # --- polarized_patch region + WAVE co-registration ---
    from ffn_sim.cell.lamellipodium_polarized_patch import (
        generate_polarized_patch_layout,
    )

    pp_layout = generate_polarized_patch_layout(
        p, wave_tag_start=0, R_cell=R_cell, polarization=(1.0, 0.0, 0.0),
        half_angle_azimuth=math.pi / 9.0, half_angle_linear=math.pi / 9.0,
    )
    polar = polarized_patch_region(
        manifold, R_cell=R_cell, polarization=(1.0, 0.0, 0.0),
        half_angle_azimuth=math.pi / 6.0, half_angle_linear=math.pi / 6.0,
    )
    pp_home = manifold.nearest_patch(pp_layout.layout.wave_positions)
    pp_coreg = bool(np.all(polar.mask[pp_home])) if pp_home.size else True

    out_png = _OUT_DIR / "h7_manifold_regions.png"
    _make_figure(
        manifold, basal, polar,
        wave_basal_home=br_home, wave_polar_home=pp_home,
        R_cell=R_cell, out_png=out_png,
    )

    sweep = _resolution_sweep(R_cell)
    summary = {
        "R_cell_m": R_cell,
        "subdivisions": int(args.subdivisions),
        "n_tri": manifold.n_tri,
        "basal_ring": {
            "n_patches": int(basal.patch_ids.size),
            "area_fraction": basal.area_fraction,
            "theta_c_rad": basal.meta["theta_c_rad"],
            "collar_half_angle_rad": basal.meta["collar_half_angle_rad"],
            "wave_coregistration_pass": br_coreg,
        },
        "polarized_patch": {
            "n_patches": int(polar.patch_ids.size),
            "area_fraction": polar.area_fraction,
            "cap_radius_rad": polar.meta["cap_radius_rad"],
            "wave_coregistration_pass": pp_coreg,
        },
        "resolution_sweep": sweep,
        "figure": str(out_png),
    }
    out_json = _OUT_DIR / "h7_manifold_regions.json"
    out_json.write_text(json.dumps(summary, indent=2))

    print(json.dumps(summary, indent=2))
    gate_ok = br_coreg and pp_coreg
    print(
        f"\nco-registration gate: basal_ring={'PASS' if br_coreg else 'FAIL'}, "
        f"polarized_patch={'PASS' if pp_coreg else 'FAIL'}"
    )
    print(f"figure → {out_png}")
    return 0 if gate_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
