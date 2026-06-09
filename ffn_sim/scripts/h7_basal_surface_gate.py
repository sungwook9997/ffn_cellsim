"""S1 basal-SURFACE gate (KU-3.5, PI directive (a) 2026-06-09).

The S-layer: a 2D surface mesh (surface_manifold) that POSITIONS FA particles on the
basal cap + carries patch connectivity ("표면 연락 일부") — NO force, not a bond
(γ-invisible). Build-time controls: FA ON the surface (radial=R), inside the
south-pole basal cap, basal patches connected, and the VG-1 resolution-invariance
master gate (cap area fraction independent of subdivision level). Auto-viz: 3D cap +
FA + side view.

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
    build_basal_surface,
    cap_fraction_resolution_invariance,
)

_OUT = _HERE.parents[1] / "outputs" / "h7"

R = 7.5e-6
FOOT_R = 5.0e-6
N_FA = 60
SUBDIV = 3


def run() -> dict:
    surf = build_basal_surface(
        R=R, footprint_radius=FOOT_R, n_fa=N_FA, subdivisions=SUBDIV,
        rng=np.random.default_rng(7),
    )
    rep = basal_surface_report(surf, footprint_radius=FOOT_R)
    inv = cap_fraction_resolution_invariance(R=R, footprint_radius=FOOT_R)
    rep["resolution_invariance"] = inv
    rep["verdict"] = "PASS" if (rep["verdict"] == "PASS" and inv["invariant"]) else "REVIEW"
    rep["gate"] = "S1 basal_surface (positioning + connectivity, NO force)"
    rep["params"] = {"R_m": R, "footprint_radius_m": FOOT_R, "n_fa": N_FA,
                     "subdivisions": SUBDIV}
    rep["note"] = (
        "S-layer = surface_manifold positions FA ON the 2D surface + patch "
        "connectivity; NO in-plane/edge force (γ-invisible), per the 7b19276 "
        "cortex-as-mesh rejection bright line. Real forces are the F-layer "
        "explicit filament network (B1/B2) anchored to these FA — wired in B3."
    )
    _figure(surf, rep)
    return rep


def _figure(surf, rep) -> None:
    m = surf.manifold
    fig = plt.figure(figsize=(12, 5.2))
    ax0 = fig.add_subplot(1, 2, 1, projection="3d")
    cent = np.asarray(m.tri_centroids) * 1e6
    # all patches (faint) + basal patches (blue) + FA (red).
    ax0.scatter(cent[:, 0], cent[:, 1], cent[:, 2], s=2, color="0.85", alpha=0.3)
    b = np.asarray(surf.basal_tris)
    ax0.scatter(cent[b, 0], cent[b, 1], cent[b, 2], s=6, color="#1f77b4",
                alpha=0.6, label="basal patches")
    fa = np.asarray(surf.fa_positions) * 1e6
    ax0.scatter(fa[:, 0], fa[:, 1], fa[:, 2], s=18, color="#d62728",
                label="FA on surface")
    ax0.set_title(f"S-layer: {rep['controls']['n_fa_placed']} FA on the basal cap "
                  f"(θ={rep['controls']['theta_cap_deg']:.0f}°)")
    ax0.legend(fontsize=8)
    ax0.set_xlabel("x [µm]"); ax0.set_ylabel("y [µm]"); ax0.set_zlabel("z [µm]")

    # side view x-z: FA on the surface, south cap.
    ax1 = fig.add_subplot(1, 2, 2)
    ax1.scatter(cent[:, 0], cent[:, 2], s=2, color="0.85", alpha=0.3)
    ax1.scatter(cent[b, 0], cent[b, 2], s=6, color="#1f77b4", alpha=0.6)
    ax1.scatter(fa[:, 0], fa[:, 2], s=18, color="#d62728")
    ax1.set_aspect("equal", adjustable="datalim")
    ax1.set_xlabel("x [µm]"); ax1.set_ylabel("z [µm]")
    inv = rep["resolution_invariance"]
    ax1.set_title(f"side view (south cap); res-invariant cap frac "
                  f"(rel spread {inv['rel_spread']:.1%})")

    fig.suptitle(f"S1 basal_surface gate — {rep['verdict']} | on_surface="
                 f"{rep['controls']['on_surface']} in_cap={rep['controls']['in_basal_cap']} "
                 f"connected={rep['controls']['basal_connected']}", fontsize=10)
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
    print(f"[S1 basal_surface gate] {rep['verdict']}")
    print(f"  basal patches={c['n_basal_tris']}, FA placed={c['n_fa_placed']}, "
          f"theta_cap={c['theta_cap_deg']:.1f}°, cap area frac={c['cap_area_fraction']:.4f}")
    print(f"  on_surface={c['on_surface']} in_basal_cap={c['in_basal_cap']} "
          f"basal_connected={c['basal_connected']}")
    print(f"  VG-1 resolution-invariance: invariant={inv['invariant']} "
          f"(rel spread {inv['rel_spread']:.1%}; analytic frac {inv['analytic_cap_fraction']:.4f})")
    print(f"  json: {path}")
    return 0 if rep["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
