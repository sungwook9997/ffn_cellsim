"""B4 basal-NMII gate (KU-3.5, PI directive (a) 2026-06-09).

Places sf_myosin_ NMII minifilaments ON the integrated basal apparatus (B3), reusing
the cortex Stam-Hocky/Hill machinery via the ①a prefix split (sf_myosin_ → sf_
γ-denylist, never contaminating cortical active-γ) and the ①b actin-pool
generalization (heads bind the apparatus filament beads). Build-time controls:
minifilaments assembled, force-free intra bonds, ≥80% of heads binding-eligible
(perp-to-segment ≤ capture_perp), sf_-denylisted, distinct from cortex_myosin_. The
dynamic contraction (Kumar single-SF tension) is the equilibrated B5 gate.
Auto-viz: apparatus + minifilament backbones + heads.

Run:  python ffn_sim/scripts/h7_basal_nmii_gate.py
"""

from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ffn_sim.archive.hoomd_legacy.cell.basal_surface import build_flat_basal_surface
from ffn_sim.archive.hoomd_legacy.cell.basal_mesh import (
    build_basal_filament_network,
    place_sf_myosin_on_apparatus,
    sf_myosin_placement_report,
)
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin

_OUT = _HERE.parents[1] / "outputs" / "h7"

ELL0 = 0.5e-6
FOOT_R = 5.0e-6
Z_BASAL = -7.0e-6
N_FA = 60
N_CABLES = 12
N_INFILL = 500
N_MOTORS = 40


def run() -> dict:
    surf = build_flat_basal_surface(
        footprint_radius=FOOT_R, z_basal=Z_BASAL, n_rings=10, n_fa=N_FA,
        rng=np.random.default_rng(7),
    )
    app = build_basal_filament_network(
        surf, ell0=ELL0, n_cables=N_CABLES, n_infill=N_INFILL,
        rng=np.random.default_rng(8),
    )
    cfg = deepcopy(yaml.safe_load(open(_HERE.parents[1] / "configs" / "phase1_h3.yaml")))
    cfg["cortex"]["myosin"]["prefix"] = "sf_myosin_"
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = N_MOTORS
    p_sf = resolve_cortex_myosin(cfg, dt=1e-9)
    myo = place_sf_myosin_on_apparatus(app, p_sf, rng=np.random.default_rng(11))
    rep = sf_myosin_placement_report(app, p_sf, myo)
    rep["gate"] = "B4 basal sf_myosin_ NMII placement"
    rep["params"] = {"n_motors": N_MOTORS, "prefix": p_sf.prefix,
                     "n_backbone": p_sf.n_backbone, "n_heads_per_side": p_sf.n_heads_per_side}
    rep["note"] = (
        "sf_myosin_ NMII placed ON the basal apparatus (prefix split ①a + actin-pool "
        "①b). Force-free + ≥80% heads binding-eligible + sf_-denylisted (cortical "
        "active-γ untouched). Dynamic contraction → Kumar single-SF tension is the "
        "equilibrated B5 gate (needs prelude + PI N_filaments/k_actin + gbook GPU)."
    )
    _figure(app, myo, p_sf, rep)
    return rep


def _figure(app, myo, p_sf, rep) -> None:
    lay = app.layout
    pos = np.asarray(lay.positions_flat) * 1e6
    starts = np.asarray(lay.filament_starts, dtype=np.int64)
    nb = np.asarray(lay.n_beads_per_filament, dtype=np.int64)
    is_c = np.asarray(lay.is_formin, dtype=bool)
    N = p_sf.n_backbone
    mp = np.asarray(myo.positions) * 1e6  # (M, N+2H, 3)

    fig, ax = plt.subplots(figsize=(7.5, 7.0))
    for f in range(starts.shape[0]):
        s, n = int(starts[f]), int(nb[f])
        seg = pos[s:s + n]
        col = "#d62728" if is_c[f] else "#2ca02c"
        lw = 1.4 if is_c[f] else 0.4
        ax.plot(seg[:, 0], seg[:, 1], color=col, lw=lw, alpha=0.5, zorder=1)
    # minifilament backbones (blue) + heads (orange dots).
    for m in range(mp.shape[0]):
        bb = mp[m, :N]
        ax.plot(bb[:, 0], bb[:, 1], color="#1f77b4", lw=2.0, zorder=3)
        hd = mp[m, N:]
        ax.scatter(hd[:, 0], hd[:, 1], s=4, color="#ff7f0e", zorder=4)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x [µm]"); ax.set_ylabel("y [µm]")
    hn = rep["controls"]["heads_near_actin"]
    ax.set_title(f"B4 basal NMII — {rep['verdict']}\n"
                 f"{rep['params']['n_motors']} sf_myosin_ minifilaments (blue) + heads "
                 f"(orange) on cables (red)/infill (green)\n"
                 f"{hn['eligible_fraction']*100:.0f}% heads binding-eligible "
                 f"(median perp {hn['median_perp_m']*1e9:.0f} nm ≤ cap {hn['capture_perp_m']*1e9:.0f} nm)")
    (_OUT / "figs").mkdir(parents=True, exist_ok=True)
    fig.savefig(_OUT / "figs" / "h7_basal_nmii_gate.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    rep = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_basal_nmii_gate.json"
    with open(path, "w") as fh:
        json.dump(rep, fh, indent=2, default=str)
    c = rep["controls"]; hn = c["heads_near_actin"]; ff = c["force_free"]
    print(f"[B4 basal NMII gate] {rep['verdict']}")
    print(f"  {c['n_motors']} sf_myosin_ minifilaments, {c['n_beads_per_motor']} beads/motor, "
          f"assembled={c['assembled']}")
    print(f"  force_free={ff['ok']} (backbone strain {ff['backbone_strain']:.1e})")
    print(f"  heads binding-eligible={hn['ok']} (frac {hn['eligible_fraction']:.3f}, "
          f"median perp {hn['median_perp_m']*1e9:.0f} nm, cap {hn['capture_perp_m']*1e9:.0f} nm)")
    print(f"  sf_denylisted={c['sf_denylisted']} distinct_from_cortex_myosin={c['distinct_from_cortex_myosin']}")
    print(f"  json: {path}")
    return 0 if rep["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
