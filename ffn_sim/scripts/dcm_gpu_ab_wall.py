"""A/B wall-time benchmark: original-loop vs vectorized-active vs fused (gbook A5000).

The full-stack profiler (``dcm_gpu_profile.py``) showed the dominant per-step cost
is the active-traction Python per-cell loop (~22 ms/step, ~57% of ~39 ms total at
N=200), NOT the gpu_local context (fusing 4 forces saves only ~0.7 ms). This script
measures the HONEST end-to-end wall speedup of the two fixes on a real run:

  A. SEPARATE + original loop active  (DcmActiveRimTractionGPU)         — baseline
  B. SEPARATE + vectorized active     (DcmActiveRimTractionGPUVec)      — the fix
  C. FUSED single callback            (DcmFusedForceGPU, vec active)    — fix+fuse

Each builds the SAME N, runs the SAME steps after a warmup, times ``sim.run`` with
a device sync, and prints ms/step + speedup vs baseline.

Run on gbook:
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python scripts/dcm_gpu_ab_wall.py \
        --n-cells 200 --steps 2000 --warmup 200
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.cell.dcm_gpu_build import (
    ResolvedGpuDCM,
    build_gpu_dcm_simulation,
)
from ffn_sim.cell.dcm_gpu_forces import DcmFusedForceGPU

_OUT = Path("ffn_sim/outputs/h_dcm_gpu_lod/profile")
_ON_GPU = False


def _sync():
    if _ON_GPU:
        import cupy as cp
        cp.cuda.runtime.deviceSynchronize()


def _time_run(h, steps, warmup):
    sim = h["sim"]
    sim.run(warmup)
    _sync()
    t0 = time.perf_counter()
    sim.run(steps)
    _sync()
    wall = time.perf_counter() - t0
    finite = bool(np.all(np.isfinite(_pos(sim))))
    return wall, finite


def _pos(sim):
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    return pos[np.argsort(tag)]


def _build_fused(p, n_cells, device):
    """Separate bond+turgor REPLACED by one fused custom force (+ native bond)."""
    from ffn_sim.integrator.baoab import make_baoab_updater
    from ffn_sim.cell.dcm_gpu_build import (
        build_gpu_dcm_snapshot, pick_device,
    )
    b = build_gpu_dcm_snapshot(p, n_cells)
    snap = b["snap"]
    nv, mean_edge = b["nv"], b["mean_edge"]
    cell_of_node, ranges = b["cell_of_node"], b["ranges"]
    faces, face_cell = b["faces"], b["face_cell"]

    dev = pick_device(device)
    sim = hoomd.Simulation(device=dev, seed=p.seed)
    sim.create_state_from_snapshot(snap)
    ig = md.Integrator(dt=p.dt)
    bond = md.bond.Harmonic()
    bond.params["dcm_edge"] = dict(k=p.k_edge, r0=mean_edge)
    ig.forces.append(bond)

    area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
    W_cs = p.W_cs_Jm2 * p.ligand_density * area_per_node
    V0 = (4.0 / 3.0) * np.pi * p.R_cell ** 3
    cad_mult = np.ones(n_cells, dtype=np.float64)
    integrin_gain = np.ones(n_cells, dtype=np.float64)
    active_mask = np.ones(n_cells, dtype=bool)
    r_contact = p.r_contact_factor * mean_edge

    fused = DcmFusedForceGPU(
        n_cells=n_cells,
        faces=faces, face_cell=face_cell, V0=V0, turgor_dP0=p.turgor_dP0,
        K_vol=p.K_vol,
        cell_of_node=cell_of_node, r_contact=r_contact, c_adh=p.c_adh,
        rep_strength=p.rep_strength, adh_strength=p.adh_strength,
        patch_area=area_per_node, contact_force_cap=p.contact_force_cap,
        cad_mult=cad_mult,
        z0=p.z_substrate, W_cs=W_cs, adh_range=p.R_cell, k_sub=p.k_sub_Nm,
        with_active=True, ranges=ranges, active=active_mask,
        int_mult=integrin_gain, R_cell=p.R_cell, f_act=1.2e-10, f_cap=6.0e-10,
        ramp_steps=4000, contact_band=0.5, neighbour_factor=2.6,
        max_neighbours=9, belt_factor=0.25)
    ig.forces.append(fused)
    sim.operations.integrator = ig
    sim.run(0)
    gamma = {"dcm_mem": p.gamma_node, "dcm_inert": p.gamma_node}
    _, updater = make_baoab_updater(kT=p.kT, gamma=gamma, dt=p.dt, seed=p.seed + 1)
    sim.operations.updaters.append(updater)
    return dict(sim=sim, fused=fused)


def main() -> None:
    global _ON_GPU
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=200)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--warmup", type=int, default=200)
    ap.add_argument("--cpu", action="store_true")
    args = ap.parse_args()
    _OUT.mkdir(parents=True, exist_ok=True)

    p = ResolvedGpuDCM()
    device = hoomd.device.CPU(notice_level=0) if args.cpu else None

    # A. baseline (original loop active)
    hA = build_gpu_dcm_simulation(p, args.n_cells, device=device,
                                  active=True, fast_active=False)
    _ON_GPU = isinstance(hA["sim"].device, hoomd.device.GPU)
    wallA, okA = _time_run(hA, args.steps, args.warmup)

    # B. vectorized active (separate forces)
    hB = build_gpu_dcm_simulation(p, args.n_cells, device=device,
                                  active=True, fast_active=True)
    wallB, okB = _time_run(hB, args.steps, args.warmup)

    # C. fused single callback (vec active inside)
    hC = _build_fused(p, args.n_cells, device)
    wallC, okC = _time_run(hC, args.steps, args.warmup)

    steps = args.steps
    msA, msB, msC = 1e3 * wallA / steps, 1e3 * wallB / steps, 1e3 * wallC / steps
    dev = type(hA["sim"].device).__name__
    print(f"\n=== A/B WALL  N={args.n_cells} nodes={hA['sim'].state.N_particles} "
          f"{dev}  {steps} steps (after {args.warmup} warmup) ===")
    print(f"  A. separate + ORIGINAL loop active : {msA:8.3f} ms/step  "
          f"(wall {wallA:7.2f}s, finite={okA})   [baseline 1.00x]")
    print(f"  B. separate + VECTORIZED active    : {msB:8.3f} ms/step  "
          f"(wall {wallB:7.2f}s, finite={okB})   [{wallA / wallB:.2f}x]")
    print(f"  C. FUSED single callback (vec)     : {msC:8.3f} ms/step  "
          f"(wall {wallC:7.2f}s, finite={okC})   [{wallA / wallC:.2f}x]")
    print(f"\n  vectorize-active speedup  = {wallA / wallB:.2f}x")
    print(f"  +fuse extra over vec      = {wallB / wallC:.2f}x")
    print(f"  total (fused+vec) speedup = {wallA / wallC:.2f}x")

    out = dict(device=dev, n_cells=args.n_cells,
               nodes=int(hA["sim"].state.N_particles), steps=steps,
               warmup=args.warmup,
               A_original_loop=dict(ms_per_step=msA, wall_s=wallA, finite=okA),
               B_vectorized=dict(ms_per_step=msB, wall_s=wallB, finite=okB,
                                 speedup_vs_A=wallA / wallB),
               C_fused=dict(ms_per_step=msC, wall_s=wallC, finite=okC,
                            speedup_vs_A=wallA / wallC,
                            speedup_over_B=wallB / wallC))
    fp = _OUT / f"ab_wall_n{args.n_cells}_{dev}.json"
    fp.write_text(json.dumps(out, indent=2))
    print(f"\n[json] {fp}")


if __name__ == "__main__":
    main()
