"""S1 FLAT basal-SURFACE gate (KU-3.5, PI directive (a) 2026-06-09).

The S-layer: a FLAT 2D triangulated surface mesh (the adherent ventral surface) that
POSITIONS FA particles + carries patch connectivity ("표면 연락 일부") — NO force, not
a bond (γ-invisible). Adherent cells flatten their ventral surface against the
substrate (PI 2026-06-09), so the basal surface is flat (not a sphere cap) — also
makes the lamellipodium natural to add. Build-time controls: FA ON the flat surface
(z=z_basal), inside the footprint disk, triangulation connected, and the VG-1
resolution-invariance master gate (disk area independent of ring resolution).
Auto-viz: top triangulation + FA + side view (planarity).

Run:  python ffn_sim/scripts/h7_basal_surface_gate.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffn_sim.cell.basal_surface import (
    basal_surface_report,
    build_flat_basal_surface,
    disk_area_resolution_invariance,
)

_OUT = _HERE.parents[1] / "outputs" / "h7"

FOOT_R = 5.0e-6
Z_BASAL = -7.0e-6
N_FA = 60
N_RINGS = 10


def run() -> dict:
    surf = build_flat_basal_surface(
        footprint_radius=FOOT_R, z_basal=Z_BASAL, n_rings=N_RINGS, n_fa=N_FA,
        rng=np.random.default_rng(7),
    )
    rep = basal_surface_report(surf)
    inv = disk_area_resolution_invariance(footprint_radius=FOOT_R, z_basal=Z_BASAL)
    rep["resolution_invariance"] = inv
    rep["verdict"] = "PASS" if (rep["verdict"] == "PASS" and inv["invariant"]) else "REVIEW"
    rep["gate"] = "S1 flat basal_surface (positioning + connectivity, NO force)"
    rep["params"] = {"footprint_radius_m": FOOT_R, "z_basal_m": Z_BASAL,
                     "n_fa": N_FA, "n_rings": N_RINGS}
    rep["note"] = (
        "S-layer = FLAT ventral surface mesh (adherent cells flatten the ventral "
        "surface; PI 2026-06-09). Positions FA ON the 2D surface + patch connectivity; "
        "NO in-plane/edge force (γ-invisible), per the 7b19276 cortex-as-mesh "
        "rejection. Real forces = the F-layer filament network (B1/B2) anchored to "
        "these FA — wired in B3. Flat surface also eases the lamellipodium."
    )
    _figure(surf, rep)
    return rep


def _figure(surf, rep) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 5.2))
    v = np.asarray(surf.verts) * 1e6
    tris = np.asarray(surf.tris)
    # top view: triangulation + FA.
    ax0.triplot(v[:, 0], v[:, 1], tris, color="#1f77b4", lw=0.4, alpha=0.6)
    fa = np.asarray(surf.fa_positions) * 1e6
    ax0.scatter(fa[:, 0], fa[:, 1], s=18, color="#d62728", zorder=3, label="FA on surface")
    th = np.linspace(0, 2 * np.pi, 200)
    ax0.plot(FOOT_R * 1e6 * np.cos(th), FOOT_R * 1e6 * np.sin(th), "k--", lw=1.0,
             label="footprint")
    ax0.set_aspect("equal", adjustable="datalim")
    ax0.set_xlabel("x [µm]"); ax0.set_ylabel("y [µm]")
    c = rep["controls"]
    ax0.set_title(f"FLAT ventral surface: {c['n_tri']} tris, {c['n_fa_placed']} FA\n"
                  f"area {c['area_fraction_of_disk']:.1%} of disk")
    ax0.legend(fontsize=8)

    # side view x-z: planarity at z_basal.
    ax1.scatter(v[:, 0], v[:, 2], s=3, color="#1f77b4", alpha=0.5, label="surface verts")
    ax1.scatter(fa[:, 0], fa[:, 2], s=18, color="#d62728", label="FA")
    ax1.set_xlabel("x [µm]"); ax1.set_ylabel("z [µm]")
    inv = rep["resolution_invariance"]
    ax1.set_title(f"side view (FLAT at z_basal); res-invariant area "
                  f"(rel spread {inv['rel_spread']:.1%})")
    ax1.legend(fontsize=8)

    fig.suptitle(f"S1 flat basal_surface gate — {rep['verdict']} | planar="
                 f"{c['planar']} within_disk={c['within_disk']} connected={c['connected']}",
                 fontsize=10)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_basal_surface_gate.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    rep = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_basal_surface_gate.json"
    with open(path, "w") as fh:
        json.dump(rep, fh, indent=2, default=str)
    c = rep["controls"]; inv = rep["resolution_invariance"]
    print(f"[S1 flat basal_surface gate] {rep['verdict']}")
    print(f"  {c['n_tri']} tris / {c['n_vert']} verts, FA placed={c['n_fa_placed']}, "
          f"area={c['area_fraction_of_disk']:.1%} of disk")
    print(f"  planar={c['planar']} (max z dev {c['max_z_dev']:.1e}) "
          f"within_disk={c['within_disk']} connected={c['connected']}")
    print(f"  VG-1 resolution-invariance: invariant={inv['invariant']} "
          f"(rel spread {inv['rel_spread']:.1%})")
    print(f"  json: {path}")
    return 0 if rep["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
