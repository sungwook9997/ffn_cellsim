"""H.7 Option-A — ManifoldIndex validation: VG-6 candidate identity + staleness/drift gate.

Option A (PI 2026-06-09) wires the geometry manifold into the cortex as a SEARCH/COORDINATE
accelerator. Before any binder is rewired, the 2026-06-09 adversarial audit
(``H7_2D_MESH_AB_DECISION_2026-06-09.md``) demands two correctness proofs at production scale:

  VG-6  (broad-phase neutrality): at refresh time, the index's within-reach pair set is
        IDENTICAL to the global cKDTree query — wiring it in changes no connectivity.

  STALENESS GATE (the NEW hole the audit surfaced): a lazily-refreshed (stale) bead→patch
        map silently DROPS within-reach pairs once beads drift. This sweeps the missed-pair
        fraction vs bead drift, confirms it is ~0 below the geometry-derived drift threshold
        (frac·circumradius) and grows past it, and so derives the safe refresh policy a
        binder MUST honour. (NOT a γ test — this index is orthogonal to the active-γ floor.)

Geometry only: builds the real connected-cortex bead cloud (seed, no sim) for honest
positions; no mechanics, no integrator, no production-config change.

Usage:
    python -m ffn_sim.scripts.h7_manifold_index_validate --n-filaments 1000 --subdiv 3
"""
from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import load_manifest, _deep_merge
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cortex.connected_mesh import build_connected_cortex, mesoscale_reach
from ffn_sim.cortex.surface_manifold import SurfaceManifold
from ffn_sim.cortex.manifold_index import ManifoldIndex, global_pairs_within_reach

_UM = 1.0e6
_NM = 1.0e9


def _cortex_cloud(n_filaments, seed):
    """Real connected-cortex actin bead positions (seed only, no sim)."""
    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    base = load_manifest(manifest["base_cortex_config"])
    cfg = _deep_merge(base, manifest.get("cortex_overrides"))
    cfg.setdefault("cortex", {})["R_cell"] = float(manifest["R_cell"])
    cfg["cortex"]["n_filaments"] = int(n_filaments)
    cfg["cortex"]["demo_mode"] = True
    p_cortex = resolve_h3_derived(cfg)
    p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
    rng = np.random.default_rng(seed)
    handles = build_connected_cortex(
        p_cortex, p_xl, n_filaments=int(n_filaments),
        with_baoab=False, with_simulation=False, rng=rng,
    )
    pos = np.asarray(handles.layout.positions_flat, dtype=np.float64)
    reach = mesoscale_reach(p_cortex.R_cell, int(n_filaments))
    return pos, float(p_cortex.R_cell), reach


def run(*, n_filaments, subdiv, seed, drift_steps, out_json, out_png):
    pos0, R, reach = _cortex_cloud(n_filaments, seed)
    n = pos0.shape[0]
    manifold = SurfaceManifold.icosphere(subdivisions=subdiv, radius=R)
    idx = ManifoldIndex(manifold, reach=reach)
    idx.refresh(pos0)

    print(f"  n_beads={n}  R={R*_UM:.2f}µm  reach={reach*_NM:.0f}nm  "
          f"n_tri={manifold.n_tri}  k-ring={idx.k}  "
          f"circumradius={manifold.max_circumradius*_NM:.0f}nm  "
          f"drift_thresh={idx.drift_threshold*_NM:.0f}nm", flush=True)

    # ---- VG-6: candidate-set identity at refresh time ----
    truth0 = global_pairs_within_reach(pos0, reach)
    idx_pairs0 = idx.pairs_within_reach(pos0)
    missed0 = truth0 - idx_pairs0
    extra0 = idx_pairs0 - truth0  # narrow-phase removes extras → should be empty
    vg6_pass = (len(missed0) == 0)
    print(f"  [VG-6] truth pairs={len(truth0)}  index pairs={len(idx_pairs0)}  "
          f"missed={len(missed0)}  extra={len(extra0)}  → {'PASS' if vg6_pass else 'FAIL'}",
          flush=True)

    # ---- STALENESS GATE: miss-fraction vs drift (stale bead_tri) ----
    # refresh ONCE at pos0, then displace beads by increasing isotropic drift, keep the
    # STALE index, and measure the fraction of TRUE within-reach pairs (on the drifted
    # cloud) the stale broad-phase pool MISSES.
    rng = np.random.default_rng(seed + 7)
    unit = rng.normal(size=pos0.shape)
    unit /= np.linalg.norm(unit, axis=1, keepdims=True).clip(min=1e-30)
    drift_max_nm = np.linspace(0.0, 3.0 * idx.drift_threshold * _NM, drift_steps)
    rows = []
    for dmax_nm in drift_max_nm:
        dmax = dmax_nm / _NM
        # per-bead drift magnitude uniform in [0, dmax] along a random direction
        mag = rng.uniform(0.0, dmax, size=n) if dmax > 0 else np.zeros(n)
        pos1 = pos0 + unit * mag[:, None]
        idx.refresh(pos0)  # reset to the STALE (pos0) home-patch map
        truth1 = global_pairs_within_reach(pos1, reach)   # truth on the drifted cloud
        stale_pairs = idx.pairs_within_reach(pos1)        # stale broad phase + narrow on pos1
        missed = truth1 - stale_pairs
        miss_frac = len(missed) / max(1, len(truth1))
        actual_drift = idx.max_drift_since_refresh(pos1)
        rows.append({
            "drift_max_nm": float(dmax_nm),
            "actual_max_drift_nm": float(actual_drift * _NM),
            "n_truth_pairs": int(len(truth1)),
            "n_missed": int(len(missed)),
            "miss_fraction": float(miss_frac),
            "past_gate": bool(actual_drift >= idx.drift_threshold),
        })
        print(f"  [drift≤{dmax_nm:6.0f}nm] truth={len(truth1):5d} missed={len(missed):5d} "
              f"miss={miss_frac*100:5.1f}%  {'(past gate)' if rows[-1]['past_gate'] else ''}",
              flush=True)

    # miss fraction strictly below the gate (safe region)
    below = [r["miss_fraction"] for r in rows if not r["past_gate"]]
    max_miss_below_gate = max(below) if below else 0.0
    verdict = {
        "vg6_pass": vg6_pass,
        "vg6_missed": len(missed0),
        "vg6_extra": len(extra0),
        "drift_threshold_nm": idx.drift_threshold * _NM,
        "max_miss_fraction_below_gate": max_miss_below_gate,
        "max_miss_fraction_overall": max(r["miss_fraction"] for r in rows),
        "gate_holds": bool(max_miss_below_gate < 0.01),  # <1% missed inside the safe region
    }
    print("=" * 78, flush=True)
    print(f"  VG-6 candidate identity: {'PASS' if vg6_pass else 'FAIL'}", flush=True)
    print(f"  staleness gate: max miss below drift-threshold "
          f"({verdict['drift_threshold_nm']:.0f}nm) = {max_miss_below_gate*100:.2f}%  "
          f"→ {'HOLDS (<1%)' if verdict['gate_holds'] else 'LEAKS'}", flush=True)
    print(f"  worst-case miss (3× threshold drift) = "
          f"{verdict['max_miss_fraction_overall']*100:.1f}%  (confirms the audit hole: "
          f"stale index drops pairs once past the gate)", flush=True)
    print("=" * 78, flush=True)

    out = {
        "scale": {"n_beads": n, "n_filaments": n_filaments, "R_um": R * _UM,
                  "reach_nm": reach * _NM, "n_tri": manifold.n_tri, "k_ring": idx.k,
                  "circumradius_nm": manifold.max_circumradius * _NM},
        "note": "Option-A ManifoldIndex validation; geometry only; NOT a γ test; no config change.",
        "verdict": verdict,
        "drift_sweep": rows,
    }
    Path(out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(out_json).write_text(json.dumps(out, indent=2))
    print(f"  json → {out_json}", flush=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    ax.plot([r["actual_max_drift_nm"] for r in rows],
            [r["miss_fraction"] * 100 for r in rows], "o-", color="#c0392b",
            label="missed within-reach pairs (stale index)")
    ax.axvline(verdict["drift_threshold_nm"], color="green", ls="--", lw=1.5,
               label=f"refresh gate ({verdict['drift_threshold_nm']:.0f}nm = "
                     f"{idx.refresh_drift_frac:.1f}·circumradius)")
    ax.axhline(1.0, color="grey", ls=":", lw=1, label="1% tolerance")
    ax.set_xlabel("max bead drift since refresh (nm)")
    ax.set_ylabel("missed within-reach binding pairs (%)")
    ax.set_title(f"Option-A staleness gate (n_beads={n}, n_tri={manifold.n_tri}) — "
                 f"VG-6 {'PASS' if vg6_pass else 'FAIL'}")
    ax.legend(fontsize=8)
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    print(f"  figure → {out_png}", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=1000)
    ap.add_argument("--subdiv", type=int, default=3, help="icosphere subdivision level (3=1280 tri)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--drift-steps", type=int, default=10)
    ap.add_argument("--out-json", type=str,
                    default="ffn_sim/outputs/h7/production/h7_manifold_index_validate.json")
    ap.add_argument("--out-png", type=str,
                    default="ffn_sim/outputs/h7/figs/h7_manifold_index_validate.png")
    args = ap.parse_args()
    return run(n_filaments=args.n_filaments, subdiv=args.subdiv, seed=args.seed,
               drift_steps=args.drift_steps, out_json=args.out_json, out_png=args.out_png)


if __name__ == "__main__":
    raise SystemExit(main())
