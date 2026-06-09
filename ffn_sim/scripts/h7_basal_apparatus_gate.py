"""B3 basal-APPARATUS gate (KU-3.5, PI directive (a) 2026-06-09).

Integrates the two layers: the F-layer explicit filament network is built ON the
S-layer flat ventral surface, with the LONG cables (ventral stress fibers) spanning
FA→FA (ends AT the FA particles = force-free sf_anchor) and the SHORT infill filling
the basal actin slab. Build-time controls: cables anchored FA→FA, anchors force-free,
filaments within the basal slab [z_basal, z_basal+band], per-filament force-free
chains, and the combined mesh percolates (B2: giant>=0.9, z in [3,3.5]). The surface
carries NO force (γ-invisible). Auto-viz: top (surface + cables + infill + FA) + side.

Run:  python ffn_sim/scripts/h7_basal_apparatus_gate.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffn_sim.cell.basal_surface import build_flat_basal_surface
from ffn_sim.cell.basal_mesh import (
    basal_apparatus_report,
    basal_connectivity_report,
    build_basal_filament_network,
    connect_basal_mesh,
)
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers

_OUT = _HERE.parents[1] / "outputs" / "h7"

ELL0 = 0.5e-6
FOOT_R = 5.0e-6
Z_BASAL = -7.0e-6
N_FA = 60
N_CABLES = 12
N_INFILL = 500


def run() -> dict:
    surf = build_flat_basal_surface(
        footprint_radius=FOOT_R, z_basal=Z_BASAL, n_rings=10, n_fa=N_FA,
        rng=np.random.default_rng(7),
    )
    app = build_basal_filament_network(
        surf, ell0=ELL0, n_cables=N_CABLES, n_infill=N_INFILL,
        rng=np.random.default_rng(8),
    )
    rep = basal_apparatus_report(app, ell0=ELL0)
    # combined-mesh connectivity (B2 on the integrated layout).
    _h3 = _HERE.parents[1] / "configs" / "phase1_h3.yaml"
    p_xl = resolve_crosslinkers(yaml.safe_load(open(_h3)), dt=1e-9)
    seed = connect_basal_mesh(app.layout, p_xl, footprint_radius=FOOT_R,
                              z_struct=3.3, rng=np.random.default_rng(9))
    conn = basal_connectivity_report(seed)
    rep["connectivity"] = conn
    rep["verdict"] = "PASS" if (rep["verdict"] == "PASS" and conn["verdict"] == "PASS") else "REVIEW"
    rep["gate"] = "B3 basal_apparatus (F-layer ON S-layer, cables FA→FA)"
    rep["params"] = {"ell0_m": ELL0, "footprint_radius_m": FOOT_R, "z_basal_m": Z_BASAL,
                     "n_fa": N_FA, "n_cables": N_CABLES, "n_infill": N_INFILL}
    rep["note"] = (
        "2-LAYER INTEGRATION. F-layer explicit filaments on the S-layer flat surface; "
        "cables span FA→FA (force-free sf_anchor); infill fills the basal actin slab; "
        "combined mesh percolates. Surface = positioning + connectivity ONLY (no "
        "force, γ-invisible). Next: B4 sf_myosin_ NMII on the apparatus + B5 "
        "equilibrated Kumar single-SF tension + Balaban traction."
    )
    _figure(surf, app, rep)
    return rep


def _figure(surf, app, rep) -> None:
    lay = app.layout
    pos = np.asarray(lay.positions_flat) * 1e6
    starts = np.asarray(lay.filament_starts, dtype=np.int64)
    nb = np.asarray(lay.n_beads_per_filament, dtype=np.int64)
    is_c = np.asarray(lay.is_formin, dtype=bool)
    fa = np.asarray(app.fa_positions) * 1e6

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 5.4))
    sv = np.asarray(surf.verts) * 1e6
    ax0.triplot(sv[:, 0], sv[:, 1], np.asarray(surf.tris), color="0.85", lw=0.3, alpha=0.5)
    for f in range(starts.shape[0]):
        s, n = int(starts[f]), int(nb[f])
        seg = pos[s:s + n]
        if is_c[f]:
            ax0.plot(seg[:, 0], seg[:, 1], color="#d62728", lw=1.6, alpha=0.9, zorder=3)
        else:
            ax0.plot(seg[:, 0], seg[:, 1], color="#2ca02c", lw=0.5, alpha=0.4, zorder=1)
    ax0.scatter(fa[:, 0], fa[:, 1], s=22, color="k", marker="x", zorder=4, label="FA")
    ax0.set_aspect("equal", adjustable="datalim")
    ax0.set_xlabel("x [µm]"); ax0.set_ylabel("y [µm]")
    c = rep["controls"]
    ax0.set_title(f"basal apparatus: {c['n_cables']} cables (red, FA→FA) + "
                  f"{c['n_infill']} infill (green)\nsurface triangulation (grey), FA (×)")
    ax0.legend(fontsize=8)

    ax1.scatter(pos[:, 0], pos[:, 2], s=2, color="#1f77b4", alpha=0.4)
    ax1.axhline(Z_BASAL * 1e6, color="k", ls="--", lw=1.0, label="z_basal (surface)")
    ax1.axhline((Z_BASAL + app.band_thickness) * 1e6, color="0.5", ls=":", lw=1.0,
                label="z_basal + band")
    ax1.set_xlabel("x [µm]"); ax1.set_ylabel("z [µm]")
    cn = rep["connectivity"]["controls"]
    ax1.set_title(f"side: actin slab on the flat surface\n"
                  f"connectivity giant={cn['giant_fraction']['value']:.3f} "
                  f"z={cn['coordination_z']['value']:.2f}")
    ax1.legend(fontsize=8)

    fig.suptitle(f"B3 basal_apparatus gate — {rep['verdict']} | cables_anchored="
                 f"{c['cables_anchored']} in_slab={c['in_slab']['ok']} "
                 f"chains_ff={c['chains_force_free']['ok']}", fontsize=10)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_basal_apparatus_gate.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    rep = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_basal_apparatus_gate.json"
    with open(path, "w") as fh:
        json.dump(rep, fh, indent=2, default=str)
    c = rep["controls"]; cn = rep["connectivity"]["controls"]
    print(f"[B3 basal_apparatus gate] {rep['verdict']}")
    print(f"  {c['n_cables']} cables (FA→FA) + {c['n_infill']} infill, "
          f"{c['n_beads_total']} beads, {c['n_anchor_bonds']} anchors")
    print(f"  cables_anchored={c['cables_anchored']} anchors_ff={c['anchors_force_free']['ok']} "
          f"in_slab={c['in_slab']['ok']} chains_ff={c['chains_force_free']['ok']} "
          f"(strain {c['chains_force_free']['max_per_filament_strain']:.1e})")
    print(f"  combined connectivity: giant={cn['giant_fraction']['value']:.3f} "
          f"z={cn['coordination_z']['value']:.3f} L/lc={cn['L_over_lc']['value']:.2f}")
    print(f"  json: {path}")
    return 0 if rep["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
