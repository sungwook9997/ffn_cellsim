"""TWO-STAGE DCM spheroid production (PI 2026-06-12) — aggregation THEN spreading.

PI directive: spreading must NOT start from a just-placed cell lattice. It must take a
spheroid that has gone through real AGGREGATION (cells pulled together into a compact
ball) and place THAT on the dish — else the rim traction is ill-defined. So:

  STAGE 1 · AGGREGATION (free-float, no dish): N cells start loosely spaced
    (spacing ~2.4·R, within the cadherin reach) and AGGREGATE under cohesion + turgor
    ALONE into a compact, rounded spheroid (Rg shrinks, contact fraction rises, the
    cluster rounds up). No substrate, no settling, no active traction.

  STAGE 2 · SPREADING (on the dish): the aggregated spheroid is translated onto the
    substrate (bottom at z0) and run WITH substrate adhesion + settling + active rim
    traction — it sediments, contacts, wets and spreads as a cohesive aggregate.

Both stages save per-frame node positions + physical diagnostics so the MacBook can
visualise them (gbook computes, MacBook renders). Diagnostics include per-cell VOLUME
conservation V/V0 (PI asked: is volume held, or do cells just flatten? — turgor must
hold V while the SHAPE flattens), radius of gyration, max height, basal footprint.

Designed for the gbook RTX A5000 (build auto-selects GPU; CPU fallback for dev). Save
output is a pickle the viz script loads.

Run (gbook, one size):
    PYTHONPATH=. python ffn_sim/scripts/dcm_two_stage_production.py --n 200 \
        --agg-steps 30000 --spread-steps 30000
Run (CPU dev smoke, tiny):
    PYTHONPATH=. python ffn_sim/scripts/dcm_two_stage_production.py --n 30 \
        --agg-steps 3000 --spread-steps 3000 --frames 6
"""

from __future__ import annotations

import argparse
import dataclasses
import pickle
import time
from pathlib import Path

import numpy as np

from ffn_sim.cell.dcm import icosphere_mesh
from ffn_sim.cell.dcm_gpu_build import (
    AggregationDrive, ResolvedGpuDCM, build_gpu_dcm_simulation, pick_device)
from ffn_sim.cell.dcm_confluence import capture_positions
from ffn_sim.common.sim_realtime import map_realtime

OUT = Path("ffn_sim/outputs/h_dcm_two_stage")
OUT.mkdir(parents=True, exist_ok=True)
UM = 1.0e6


# ---------------------------------------------------------------------------
# physical diagnostics (per frame)
# ---------------------------------------------------------------------------
def cell_volume(verts, tris):
    """Enclosed volume of one triangulated cell shell (divergence theorem)."""
    t = verts[tris]
    return abs(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum()) / 6.0


def centroids(pos, ranges):
    return np.array([pos[lo:hi].mean(0) for lo, hi in ranges])


def diagnostics(pos, ranges, tris0, *, R, z0, V0):
    """V/V0 (mean/min/max), Rg, asphericity, maxZ, basal footprint area [µm²]."""
    vols = np.array([cell_volume(pos[lo:hi], tris0) for lo, hi in ranges])
    vr = vols / V0
    cen = centroids(pos, ranges)
    cc = cen.mean(0)
    r = np.linalg.norm(cen - cc, axis=1)
    Rg = float(np.sqrt((r ** 2).mean())) if len(cen) else 0.0
    # asphericity from the gyration tensor eigenvalues (0 = sphere)
    if len(cen) >= 4:
        d = cen - cc
        G = (d[:, :, None] * d[:, None, :]).mean(0)
        ev = np.sort(np.linalg.eigvalsh(G))[::-1]
        asph = float((ev[0] - 0.5 * (ev[1] + ev[2])) / max(ev.sum(), 1e-30))
    else:
        asph = 0.0
    maxZ = float(pos[:, 2].max())
    # basal footprint (xy hull of substrate-contact nodes)
    from scipy.spatial import ConvexHull
    basal = pos[pos[:, 2] < z0 + 0.5 * R][:, :2]
    if basal.shape[0] >= 3:
        try:
            foot = float(ConvexHull(basal).volume)
        except Exception:  # noqa: BLE001
            foot = float(np.ptp(basal[:, 0]) * np.ptp(basal[:, 1]))
    else:
        foot = 0.0
    return dict(VV0_mean=float(vr.mean()), VV0_min=float(vr.min()),
                VV0_max=float(vr.max()), Rg_um=Rg * UM, asphericity=asph,
                maxZ_um=maxZ * UM, footprint_um2=foot * UM * UM)


def run_stage(sim, *, n_frames, steps_per_frame, ranges, tris0, R, z0, V0, label, dt):
    """Run a stage in chunks, capturing per-frame positions + diagnostics."""
    sim.run(0)
    frames, diags, steps = [], [], []
    pos = capture_positions(sim)
    frames.append(pos.copy())
    diags.append(diagnostics(pos, ranges, tris0, R=R, z0=z0, V0=V0))
    steps.append(int(sim.timestep))
    print(f"  [{label}] frame 0 step {sim.timestep}: "
          f"V/V0={diags[0]['VV0_mean']:.3f} Rg={diags[0]['Rg_um']:.1f}µm "
          f"maxZ={diags[0]['maxZ_um']:.1f}µm foot={diags[0]['footprint_um2']:.0f}µm²",
          flush=True)
    t0 = time.time()
    for f in range(1, n_frames):
        sim.run(steps_per_frame)
        pos = capture_positions(sim)
        if not np.isfinite(pos).all():
            print(f"  [{label}] NON-FINITE at frame {f} — truncating", flush=True)
            break
        frames.append(pos.copy())
        diags.append(diagnostics(pos, ranges, tris0, R=R, z0=z0, V0=V0))
        steps.append(int(sim.timestep))
        dd = diags[-1]
        print(f"  [{label}] frame {f} step {sim.timestep}: "
              f"V/V0={dd['VV0_mean']:.3f} Rg={dd['Rg_um']:.1f}µm "
              f"maxZ={dd['maxZ_um']:.1f}µm foot={dd['footprint_um2']:.0f}µm²", flush=True)
    wall = time.time() - t0
    return dict(frames=frames, diags=diags, steps=steps, label=label,
                wall_s=round(wall, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200, help="cell count")
    ap.add_argument("--agg-steps", type=int, default=30000)
    ap.add_argument("--spread-steps", type=int, default=30000)
    ap.add_argument("--frames", type=int, default=12)
    ap.add_argument("--agg-spacing", type=float, default=2.6,
                    help="initial cell spacing for aggregation (·R; loose start)")
    ap.add_argument("--k-agg", type=float, default=1.5e-5,
                    help="aggregation-confinement stiffness N/m (0 disables). The "
                         "centripetal pull that compacts the loose cluster into a ball.")
    ap.add_argument("--r-cell-um", type=float, default=7.5,
                    help="cell radius µm (coarse-grain to a tissue patch for big R)")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    R = args.r_cell_um * 1e-6
    V0 = (4.0 / 3.0) * np.pi * R ** 3
    z0 = 0.0
    _v, _e, tris0 = icosphere_mesh(R, 1)
    dev = pick_device(None)
    is_gpu = type(dev).__name__.endswith("GPU")
    spf_a = max(1, args.agg_steps // max(1, args.frames - 1))
    spf_s = max(1, args.spread_steps // max(1, args.frames - 1))
    print(f"[two-stage] N={args.n} R_cell={args.r_cell_um}µm device={'GPU' if is_gpu else 'CPU'} "
          f"agg={args.agg_steps} spread={args.spread_steps}\n", flush=True)

    # ---- STAGE 1: AGGREGATION (loose start, no dish, no active) ----
    p_agg = dataclasses.replace(ResolvedGpuDCM(seed=args.seed), R_cell=R,
                                spacing_factor=args.agg_spacing)
    h1 = build_gpu_dcm_simulation(p_agg, args.n, device=dev, active=False,
                                  with_substrate=False, settle_force=0.0)
    ranges = h1["ranges"]
    # AGGREGATION DRIVE: pull cells toward the initial cluster centre so the loose
    # placement COMPACTS + ROUNDS into a real cohesive spheroid (the short-range tent
    # cannot pull cells across a gap on its own). Appended BEFORE the first run.
    snap0 = h1["sim"].state.get_snapshot()
    agg_center = np.asarray(snap0.particles.position).mean(0)
    if args.k_agg > 0.0:
        h1["sim"].operations.integrator.forces.append(
            AggregationDrive(center=agg_center, k_agg=args.k_agg))
    print(f"STAGE 1 · AGGREGATION (n_cells={h1['n_cells']}, free-float, "
          f"k_agg={args.k_agg:.1e}):", flush=True)
    s1 = run_stage(h1["sim"], n_frames=args.frames, steps_per_frame=spf_a,
                   ranges=ranges, tris0=tris0, R=R, z0=z0, V0=V0, label="agg", dt=p_agg.dt)
    agg_pos = s1["frames"][-1].copy()
    # translate the aggregated ball onto the dish: bottom node at z0
    agg_pos[:, 2] -= (agg_pos[:, 2].min() - z0)
    # also recentre xy over the origin so it sits centred on the dish
    cen_xy = centroids(agg_pos, ranges).mean(0)[:2]
    agg_pos[:, 0] -= cen_xy[0]; agg_pos[:, 1] -= cen_xy[1]
    print(f"  -> aggregated: Rg {s1['diags'][0]['Rg_um']:.1f}→{s1['diags'][-1]['Rg_um']:.1f}µm, "
          f"asph {s1['diags'][0]['asphericity']:.3f}→{s1['diags'][-1]['asphericity']:.3f}, "
          f"V/V0 {s1['diags'][-1]['VV0_mean']:.3f}\n", flush=True)

    # ---- STAGE 2: SPREADING from the aggregate (on the dish, active) ----
    p = dataclasses.replace(ResolvedGpuDCM(seed=args.seed), R_cell=R)   # cohesive defaults
    h2 = build_gpu_dcm_simulation(p, args.n, device=dev, active=True,
                                  with_substrate=True, init_pos=agg_pos)
    print(f"STAGE 2 · SPREADING (from aggregate, active rim traction):", flush=True)
    s2 = run_stage(h2["sim"], n_frames=args.frames, steps_per_frame=spf_s,
                   ranges=h2["ranges"], tris0=tris0, R=R, z0=z0, V0=V0, label="spread", dt=p.dt)

    rt_a = map_realtime(s1["steps"][-1], float(p_agg.dt))
    rt_s = map_realtime(s2["steps"][-1], float(p.dt))
    tr = h2["traction"]
    out = dict(
        n_cells=int(h2["n_cells"]), R_cell=R, V0=V0, z0=z0,
        tris0=tris0, ranges=h2["ranges"], device="GPU" if is_gpu else "CPU",
        agg=s1, spread=s2, p_agg=p_agg, p=p,
        # rim-detection params so the viz can colour the lamellipodiating cells
        # exactly as the traction does (rim = few neighbours + basal contact).
        rim_params=dict(r_neigh=float(tr.r_neigh), max_neigh=int(tr.max_neigh),
                        contact_band=float(tr.contact_band), R=R, z0=z0),
        realtime=dict(agg=rt_a.t_real_human, spread=rt_s.t_real_human,
                      accel=rt_a.accel_factor))
    path = OUT / f"two_stage_n{args.n}.pkl"
    with open(path, "wb") as fh:
        pickle.dump(out, fh)
    print(f"\nwrote {path} (agg {s1['wall_s']}s + spread {s2['wall_s']}s)", flush=True)
    print(f"VOLUME held: agg V/V0={s1['diags'][-1]['VV0_mean']:.3f}, "
          f"spread V/V0={s2['diags'][-1]['VV0_mean']:.3f} "
          f"(≈1 + turgor inflation; flat = cells keep volume, shape flattens)", flush=True)


if __name__ == "__main__":
    main()
