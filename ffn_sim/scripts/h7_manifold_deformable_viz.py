"""H.7 deformable surface-manifold viz — shell slaved to a spread-cell cloud + regions.

Visualises the (b) spatial-substrate increments together: the surface manifold
fitted (``fit_to_cloud``) to a NON-spherical adherent/spread bead cloud — the
slaved-to-beads deformable shell — with the lamellipodium region masks
(``basal_ring`` collar + ``polarized_patch`` cap, ``cortex/manifold_regions.py``)
carved onto the DEFORMED shell. Shows that the contact/region infrastructure
tracks the real cell shape, not just the construction icosphere.

Geometry only (no HOOMD, no sim): a synthetic spread-cell cloud (oblate, basally
flattened) stands in for an adherent cortex; the same machinery applies to a real
cortex snapshot via ``manifold.fit_to_cloud(actin_positions)``.

Outputs:
  * ``ffn_sim/outputs/h7/figs/h7_manifold_deformable.png`` — 3D shells (icosphere
    vs fitted-to-spread-cloud) coloured by region + region-area-fraction bars.
  * ``ffn_sim/outputs/h7/figs/h7_manifold_deformable.json`` — region metrics for
    the spherical vs deformed shell.

Usage:
    python -m ffn_sim.scripts.h7_manifold_deformable_viz --subdivisions 3
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

from ffn_sim.cortex.manifold_regions import (  # noqa: E402
    basal_ring_region,
    polarized_patch_region,
)
from ffn_sim.cortex.surface_manifold import SurfaceManifold  # noqa: E402

_FIG_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"
_R_CELL = 7.5e-6


def _spread_cell_cloud(R, *, n=6000, oblate=0.55, basal_flatten=0.85, seed=3):
    """Synthetic adherent/spread cortex cloud: oblate + basal (south) flattening.

    Oblate squashes z by ``oblate``; ``basal_flatten`` pushes the south cap up
    toward a flat contact plane (a spread cell rests/flattens on the substrate).
    """
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    p = v * np.array([R, R, oblate * R])
    # Basal flattening: lift the lowest cap toward z_floor (clip the south pole).
    z_floor = -oblate * R * basal_flatten
    p[:, 2] = np.maximum(p[:, 2], z_floor)
    return p


def _region_colors(manifold, R_cell):
    br = basal_ring_region(manifold, R_cell=R_cell, rest_length=0.5e-6,
                           collar_half_angle=0.22)
    pp = polarized_patch_region(manifold, R_cell=R_cell, polarization=(1.0, 0.0, 0.0),
                                half_angle_azimuth=math.pi / 6.0)
    colors = np.full((manifold.n_tri, 4), 0.0)
    colors[:] = (0.82, 0.82, 0.82, 1.0)        # other patches: gray
    colors[br.mask] = (0.12, 0.47, 0.71, 1.0)  # basal_ring: blue
    colors[pp.mask] = (0.84, 0.19, 0.15, 1.0)  # polarized_patch: red (overrides)
    return colors, br, pp


def _draw_shell(ax, manifold, colors, title):
    verts = manifold.verts
    tris = manifold.tris
    polys = verts[tris] * 1e6  # µm
    pc = Poly3DCollection(polys, facecolors=colors, edgecolors=(0, 0, 0, 0.12),
                          linewidths=0.2)
    ax.add_collection3d(pc)
    lim = np.max(np.abs(verts)) * 1e6 * 1.05
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)"); ax.set_zlabel("z (µm)")
    ax.set_title(title, fontsize=10)
    ax.view_init(elev=-22, azim=-60)  # look up at the basal (south) cap


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--subdivisions", type=int, default=3)
    ap.add_argument("--R-cell", type=float, default=_R_CELL)
    args = ap.parse_args()
    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    R = float(args.R_cell)

    # Sphere (construction) manifold.
    m_sphere = SurfaceManifold.icosphere(args.subdivisions, R)
    col_s, br_s, pp_s = _region_colors(m_sphere, R)

    # Deformable shell fitted to a spread-cell cloud.
    cloud = _spread_cell_cloud(R)
    m_fit = SurfaceManifold.icosphere(args.subdivisions, R)
    m_fit.fit_to_cloud(cloud)
    col_f, br_f, pp_f = _region_colors(m_fit, R)

    fig = plt.figure(figsize=(16.0, 6.6), constrained_layout=True)
    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    _draw_shell(ax1, m_sphere, col_s,
                f"(a) icosphere (construction)\nbasal_ring {br_s.patch_ids.size} · "
                f"polar {pp_s.patch_ids.size}")
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    _draw_shell(ax2, m_fit, col_f,
                f"(b) deformable shell ⟵ spread-cell cloud\nbasal_ring "
                f"{br_f.patch_ids.size} · polar {pp_f.patch_ids.size}")

    ax3 = fig.add_subplot(1, 3, 3)
    x = np.arange(2)
    w = 0.35
    ax3.bar(x - w / 2, [br_s.area_fraction, pp_s.area_fraction], w,
            label="icosphere", color="0.6")
    ax3.bar(x + w / 2, [br_f.area_fraction, pp_f.area_fraction], w,
            label="deformed shell", color="tab:orange")
    ax3.set_xticks(x); ax3.set_xticklabels(["basal_ring", "polarized_patch"])
    ax3.set_ylabel("region area fraction")
    ax3.set_title("(c) region footprint: sphere vs deformed")
    ax3.legend(fontsize=8)

    fig.suptitle(
        "H.7 deformable surface-manifold: shell slaved to a spread-cell cloud, "
        "lamellipodium regions tracked on the deformed shape",
        fontsize=12,
    )
    out_png = _FIG_DIR / "h7_manifold_deformable.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

    summary = {
        "R_cell_m": R, "subdivisions": args.subdivisions, "n_tri": m_sphere.n_tri,
        "cloud_extent_um": {
            "x": float(np.ptp(cloud[:, 0]) * 1e6),
            "y": float(np.ptp(cloud[:, 1]) * 1e6),
            "z": float(np.ptp(cloud[:, 2]) * 1e6),
        },
        "sphere": {"total_area_m2": m_sphere.total_area(),
                   "basal_ring_area_frac": br_s.area_fraction,
                   "polarized_patch_area_frac": pp_s.area_fraction},
        "deformed": {"total_area_m2": m_fit.total_area(),
                     "basal_ring_area_frac": br_f.area_fraction,
                     "polarized_patch_area_frac": pp_f.area_fraction},
        "figure": str(out_png),
    }
    (_FIG_DIR / "h7_manifold_deformable.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"figure → {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
