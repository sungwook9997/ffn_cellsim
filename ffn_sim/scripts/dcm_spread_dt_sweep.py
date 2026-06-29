"""Spreading dt-stability sweep (PI 2026-06-12).

QUESTION: can a LARGER BAOAB timestep cover the same real-time CHEAPLY without
losing stability or changing the physics? The spreading stage currently uses
dt=1e-9 s ("strong adhesion → small dt"), so ~10 min real = 1e6 steps (~5 h wall).
If dt can go 10× / 100× larger and stay BAOAB-stable, the SAME real-time costs
1e5 / 1e4 steps (~33 / ~3 min wall) — long, physiological spreading becomes feasible.

METHOD (clean, no re-aggregation): reuse the CONVERGED N=100 aggregate from
``outputs/h_dcm_two_stage/two_stage_n100.pkl`` (stage-1 final, dish-placed exactly
as the production driver does) and run ONLY stage-2 spreading at each dt, each to
the SAME real-time target ``t_real = S·dt·steps`` (S=6e5). Because the frames are
sampled at the same real-time fractions, frame i is directly comparable ACROSS dt:
if a bigger dt is stable + converged, its A/A0(t_real) must track the small-dt one.

OBSERVABLE: top-down silhouette A/A0 (the experimental assay observable — xy convex
hull of ALL nodes / its t=0 value), NOT the basal contact footprint. Arrest is OFF.

STABILITY = no non-finite truncation (spread() truncates on NaN/Inf) AND V/V0 held.

  dt=1e-9 (1e6 steps) ~ dt=1e-8 (1e5, 10×) ~ dt=1e-7 (1e4, 100×) — all ~10 min real.

Run (gbook, the 10× and 100× tests):
    PYTHONPATH=. python ffn_sim/scripts/dcm_spread_dt_sweep.py --dts 1e-8,1e-7
Add the dt=1e-9 reference (slow, ~5 h): --dts 1e-7,1e-8,1e-9
"""

from __future__ import annotations

import argparse
import dataclasses
import pickle
import time
from pathlib import Path

import numpy as np

from ffn_sim.archive.hoomd_legacy.cell.dcm import icosphere_mesh  # noqa: F401  (parity w/ driver imports)
from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import pick_device
from ffn_sim.scripts.dcm_two_stage_production import spread, centroids

OUT = Path("ffn_sim/outputs/h_dcm_two_stage")
S_ACCEL = 6.0e5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aggregate", default=str(OUT / "two_stage_n100.pkl"),
                    help="two-stage pkl whose stage-1 final aggregate is reused")
    ap.add_argument("--dts", default="1e-8,1e-7",
                    help="comma-sep spreading dt [s] to test (10×/100× of 1e-9)")
    ap.add_argument("--t-real-min", type=float, default=10.0,
                    help="real-time target per dt [min]; steps = t_real/(S·dt)")
    ap.add_argument("--frames", type=int, default=16)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    d = pickle.load(open(args.aggregate, "rb"))
    R = d["R_cell"]; V0 = d["V0"]; z0 = d["z0"]; n_cells = int(d["n_cells"])
    ranges = d["ranges"]; tris0 = d["tris0"]
    # converged aggregate (stage-1 final), dish-placed EXACTLY as the driver main():
    agg_pos = d["agg"]["frames"][-1].copy()
    agg_pos[:, 2] -= (agg_pos[:, 2].min() - z0)
    cxy = centroids(agg_pos, ranges).mean(0)[:2]
    agg_pos[:, 0] -= cxy[0]; agg_pos[:, 1] -= cxy[1]

    dev = pick_device(None)
    is_gpu = type(dev).__name__.endswith("GPU")
    t_real = args.t_real_min * 60.0
    dts = [float(x) for x in args.dts.split(",")]
    print(f"[spread dt-sweep] N={n_cells} dev={'GPU' if is_gpu else 'CPU'} "
          f"t_real_target={args.t_real_min}min S={S_ACCEL:.0e} dts={dts}", flush=True)

    results: dict[float, dict] = {}
    t_start = time.time()
    for dt in dts:
        steps = int(round(t_real / (S_ACCEL * dt)))
        p = dataclasses.replace(d["p"], dt=dt)
        print(f"\n=== dt={dt:.0e}s  steps={steps:,}  "
              f"(t_real={S_ACCEL * dt * steps / 60:.2f}min) ===", flush=True)
        s2, _h = spread(p, n_cells, dev=dev, init_pos=agg_pos.copy(),
                        steps=steps, frames=args.frames, R=R, z0=z0, V0=V0,
                        tris0=tris0)
        diags = s2["diags"]; A0t = s2["A0_topdown_um2"]
        finite_ok = len(diags) == args.frames     # spread() truncates on non-finite
        vv0 = [g["VV0_mean"] for g in diags]
        aa0 = [g["topdown_um2"] / A0t for g in diags]
        results[dt] = dict(steps=steps, diags=diags, A0_top=A0t,
                           finite_ok=bool(finite_ok), wall_s=s2["wall_s"])
        print(f"  -> dt={dt:.0e}: stable={finite_ok} "
              f"(frames {len(diags)}/{args.frames})  A/A0 {aa0[0]:.3f}→{aa0[-1]:.3f} "
              f"(top-down)  V/V0 {min(vv0):.3f}–{max(vv0):.3f}  wall={s2['wall_s']}s",
              flush=True)
        with open(OUT / f"spread_dt_sweep{args.tag}.pkl", "wb") as fh:   # incremental
            pickle.dump(dict(n_cells=n_cells, S=S_ACCEL, t_real_min=args.t_real_min,
                             R_cell=R, z0=z0, V0=V0, results=results), fh)

    print(f"\n=== SUMMARY (top-down A/A0, ~{args.t_real_min}min real each) ===",
          flush=True)
    for dt in dts:
        r = results[dt]
        aa0 = [g["topdown_um2"] / r["A0_top"] for g in r["diags"]]
        print(f"  dt={dt:.0e}  steps={r['steps']:>9,}  wall={r['wall_s']:>7.1f}s  "
              f"stable={str(r['finite_ok']):>5}  A/A0_final={aa0[-1]:.3f}", flush=True)
    print(f"\ntotal wall {time.time() - t_start:.0f}s; "
          f"wrote {OUT / f'spread_dt_sweep{args.tag}.pkl'}", flush=True)


if __name__ == "__main__":
    main()
