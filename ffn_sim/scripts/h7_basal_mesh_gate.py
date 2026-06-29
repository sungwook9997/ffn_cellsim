"""B1 basal-mesh GEOMETRY gate (KU-3.5, PI directive (a) 2026-06-09).

Builds the connected-basal-mesh GEOMETRY layout (bimodal cables + infill on a flat
basal disk) and checks the build-time controls: PLANAR (in the basal band), within
the footprint disk, force-free chains, long-axis-aligned cables, isotropic infill.
NO force, NO HOOMD state — geometry only (the connectivity / FA-anchor / NMII / Kumar
gates are increments B2-B5). Auto-viz: top view (cables vs infill) + side view
(planarity).

Run:  python ffn_sim/scripts/h7_basal_mesh_gate.py
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

import yaml

from ffn_sim.archive.hoomd_legacy.cell.basal_mesh import (
    basal_connectivity_report,
    basal_mesh_build_report,
    connect_basal_mesh,
    generate_basal_mesh_layout,
)
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers

_OUT = _HERE.parents[1] / "outputs" / "h7"

# Representative MCF7 basal footprint (R_cell 7.5 µm; ℓ0 0.5 µm cortex segment).
ELL0 = 0.5e-6
FOOT_R = 5.0e-6
Z_BASAL = -7.0e-6
BAND = 200.0e-9
N_FIL = 600
LONG_AXIS = np.array([1.0, 0.0, 0.0])


def run() -> dict:
    lay = generate_basal_mesh_layout(
        n_filaments=N_FIL, ell0=ELL0, footprint_radius=FOOT_R, z_basal=Z_BASAL,
        band_thickness=BAND, long_axis=LONG_AXIS, seed=73,
    )
    rep = basal_mesh_build_report(
        lay, ell0=ELL0, footprint_radius=FOOT_R, z_basal=Z_BASAL,
        band_thickness=BAND, long_axis=LONG_AXIS,
    )
    rep["gate"] = "B1+B2 basal_mesh geometry + connectivity"
    rep["params"] = {
        "n_filaments": N_FIL, "ell0_m": ELL0, "footprint_radius_m": FOOT_R,
        "z_basal_m": Z_BASAL, "band_thickness_m": BAND,
    }

    # ---- B2: connect the filament network (F-layer connectivity) ----
    _h3 = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"
    p_xl = resolve_crosslinkers(yaml.safe_load(open(_h3)), dt=1e-9)
    seed = connect_basal_mesh(
        lay, p_xl, footprint_radius=FOOT_R, z_struct=3.3, bundle_mult=2,
        rng=np.random.default_rng(1),
    )
    b2 = basal_connectivity_report(seed)
    rep["b2_connectivity"] = b2
    # combined verdict: B1 geometry AND B2 connectivity must pass.
    rep["verdict"] = "PASS" if (rep["verdict"] == "PASS" and b2["verdict"] == "PASS") else "REVIEW"

    rep["note"] = (
        "2-LAYER (F-layer = explicit filament network). B1 geometry + B2 "
        "bridge-different-filament connectivity (giant>=0.9 / z in [3,3.5] / "
        "L/lc>=5.9). The S-layer (surface_manifold positions FA particles + patch "
        "connectivity + soft normal confinement; NO force-bearing edges per the "
        "7b19276 cortex-as-mesh rejection) and B3 anchor-on-surface / B4 sf_myosin "
        "NMII / B5 equilibrated Kumar+Balaban gate are the following increments."
    )
    _figure(lay, rep)
    return rep


def _figure(lay, rep) -> None:
    pos = np.asarray(lay.positions_flat) * 1e6
    is_f = np.asarray(lay.is_formin, dtype=bool)
    starts = np.asarray(lay.filament_starts, dtype=np.int64)
    nb = np.asarray(lay.n_beads_per_filament, dtype=np.int64)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 5.2))
    # top view x-y: cables (long, red) over infill (short, grey)
    for f in range(starts.shape[0]):
        s, n = int(starts[f]), int(nb[f])
        seg = pos[s:s + n]
        if is_f[f]:
            ax0.plot(seg[:, 0], seg[:, 1], color="#d62728", lw=1.3, alpha=0.85, zorder=3)
        else:
            ax0.plot(seg[:, 0], seg[:, 1], color="0.6", lw=0.6, alpha=0.5, zorder=1)
    th = np.linspace(0, 2 * np.pi, 200)
    ax0.plot(FOOT_R * 1e6 * np.cos(th), FOOT_R * 1e6 * np.sin(th),
             "k--", lw=1.0, label="footprint")
    ax0.set_xlabel("x [µm]"); ax0.set_ylabel("y [µm]")
    ax0.set_aspect("equal", adjustable="datalim")
    c = rep["controls"]
    ax0.set_title(f"basal mesh top view: {c['n_cables']} cables (red) + "
                  f"{c['n_infill']} infill (grey)\n"
                  f"cables |cos|={c['cables_aligned']['mean_abs_cos']:.3f}, "
                  f"infill |cos|={c['infill_isotropic']['mean_abs_cos']:.3f}")
    ax0.legend(fontsize=8)

    # side view x-z: planarity in the basal band
    ax1.scatter(pos[:, 0], pos[:, 2], s=2, color="#1f77b4", alpha=0.5)
    ax1.axhline(Z_BASAL * 1e6, color="k", ls="--", lw=1.0, label="z_basal")
    ax1.axhline((Z_BASAL - BAND) * 1e6, color="0.5", ls=":", lw=1.0,
                label="z_basal − band")
    ax1.set_xlabel("x [µm]"); ax1.set_ylabel("z [µm]")
    ax1.set_title(f"side view (planar in band): planar={c['planar_ok']}, "
                  f"force-free strain={c['force_free']['max_strain']:.1e}")
    ax1.legend(fontsize=8)

    fig.suptitle(f"B1 basal_mesh GEOMETRY gate — {rep['verdict']}", fontsize=11)
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_basal_mesh_gate.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    rep = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_basal_mesh_gate.json"
    with open(path, "w") as fh:
        json.dump(rep, fh, indent=2, default=str)
    c = rep["controls"]
    print(f"[B1 basal_mesh geometry gate] {rep['verdict']}")
    print(f"  {c['n_filaments']} filaments = {c['n_cables']} cables + {c['n_infill']} infill, "
          f"{c['n_beads_total']} beads")
    print(f"  planar={c['planar_ok']} within_footprint={c['within_footprint']} "
          f"force_free={c['force_free']['ok']} (strain {c['force_free']['max_strain']:.1e})")
    print(f"  cables_aligned={c['cables_aligned']['ok']} "
          f"(|cos|={c['cables_aligned']['mean_abs_cos']:.3f}); "
          f"infill_isotropic={c['infill_isotropic']['ok']} "
          f"(|cos|={c['infill_isotropic']['mean_abs_cos']:.3f})")
    b2 = rep.get("b2_connectivity", {}).get("controls", {})
    if b2:
        print(f"  [B2 connectivity] n_xl={b2['n_xl']} "
              f"giant={b2['giant_fraction']['value']:.3f} (ok={b2['giant_fraction']['ok']}) "
              f"z={b2['coordination_z']['value']:.3f} (ok={b2['coordination_z']['ok']}) "
              f"L/lc={b2['L_over_lc']['value']:.2f} (ok={b2['L_over_lc']['ok']})")
    print(f"  json: {path}")
    return 0 if rep["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
