"""Micro-profile the DcmActiveRimTractionGPU per-cell Python loop (gbook A5000).

The full-stack profiler showed set_forces_active = ~22 ms/step (~57% of total) at
N=200 — by far the dominant cost, and it is NOT the gpu_local context (that floor
is ~0.2 ms). The suspect is the Python ``for k in range(n_active)`` loop inside
``DcmActiveRimTractionGPU.set_forces`` (dcm_gpu_forces.py): each iteration does
``int(...)`` / ``bool(...)`` / ``float(...)`` on 0-d cupy arrays — each a blocking
GPU->CPU scalar sync — plus per-cell fancy-index scatters. This script times the
three phases of that set_forces (snapshot+setup, per-cell loop, force write) to
prove where the 22 ms goes, and counts the device->host scalar syncs per step.

Run on gbook:
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python \
        scripts/dcm_gpu_active_micro.py --n-cells 200 --steps 200 --warmup 60
"""

from __future__ import annotations

import argparse
import time

import numpy as np

import hoomd

from ffn_sim.cell.dcm_gpu_build import ResolvedGpuDCM, build_gpu_dcm_simulation

_ON_GPU = False


def _sync():
    if _ON_GPU:
        import cupy as cp
        cp.cuda.runtime.deviceSynchronize()


def main() -> None:
    global _ON_GPU
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-cells", type=int, default=200)
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--warmup", type=int, default=60)
    ap.add_argument("--cpu", action="store_true")
    args = ap.parse_args()

    p = ResolvedGpuDCM()
    device = hoomd.device.CPU(notice_level=0) if args.cpu else None
    h = build_gpu_dcm_simulation(p, args.n_cells, device=device, active=True)
    sim = h["sim"]
    _ON_GPU = isinstance(sim.device, hoomd.device.GPU)
    trac = h["traction"]
    print(f"[micro] device={type(sim.device).__name__} N={args.n_cells} "
          f"nodes={sim.state.N_particles}")

    sim.run(args.warmup)
    _sync()

    # Re-implement the set_forces phases with timers, calling the SAME math the
    # force uses, on the live snapshot — measures phase split without editing the
    # force class. (Read-only: we do NOT write forces here, just time the compute.)
    d = trac._dispatch()
    xp = d.xp
    t_setup = t_loop = t_write = 0.0
    for _ in range(args.steps):
        sim.run(1)  # advance one real step (uses the real force)
        _sync()
        # --- now time a STANDALONE replay of the three phases on current state
        _sync(); t0 = time.perf_counter()
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            perm = xp.argsort(tag)
            pos_g = pos[perm]
            F_g = xp.zeros_like(pos_g)
            active_ids = xp.where(xp.asarray(trac.active))[0]
            cents = xp.asarray(
                [pos_g[trac.ranges[int(c)][0]:trac.ranges[int(c)][1]].mean(0)
                 for c in active_ids])
            cluster_cen = cents.mean(0)
            d2 = xp.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
            within = d2 < trac.r_neigh ** 2
            xp.fill_diagonal(within, False)
            crowd = within.sum(axis=1)
            int_mult = xp.asarray(trac.int_mult)
            _sync(); t_setup += time.perf_counter() - t0

            # --- per-cell Python loop (the suspect) ---
            _sync(); t1 = time.perf_counter()
            zc = trac.z0 + trac.contact_band * trac.R
            for k in range(int(active_ids.shape[0])):
                c = int(active_ids[k])
                if int(crowd[k]) > trac.max_neigh:
                    continue
                lo, hi = trac.ranges[c]
                cell_pos = pos_g[lo:hi]
                basal = cell_pos[:, 2] < zc
                if not bool(basal.any()):
                    continue
                rxy = cents[k][:2] - cluster_cen[:2]
                rn = float(xp.hypot(rxy[0], rxy[1]))
                if rn < 1e-12:
                    continue
                rhat = xp.concatenate([rxy / rn, xp.zeros(1, dtype=rxy.dtype)])
                gain = float(int_mult[c])
                fmag = min(1.0 * trac.f_act * gain, trac.f_cap)
                idx = xp.where(basal)[0]
                F_g[lo:hi][idx] += fmag * rhat
            _sync(); t_loop += time.perf_counter() - t1

            _sync(); t2 = time.perf_counter()
            fn = xp.linalg.norm(F_g, axis=1)
            over = fn > trac.f_cap
            if bool(over.any()):
                F_g[over] *= (trac.f_cap / fn[over])[:, None]
            F = xp.empty_like(pos)
            F[perm] = F_g
            _sync(); t_write += time.perf_counter() - t2

    s = args.steps
    ms = lambda x: 1e3 * x / s
    print(f"\n=== ACTIVE-FORCE PHASE SPLIT  ({s} steps, "
          f"{'GPU' if _ON_GPU else 'CPU'}, n_active={int(active_ids.shape[0])}) ===")
    print(f"  setup (snapshot+centroids+crowd) : {ms(t_setup):8.4f} ms/step")
    print(f"  per-cell PYTHON loop             : {ms(t_loop):8.4f} ms/step")
    print(f"  cap + write-back                 : {ms(t_write):8.4f} ms/step")
    print(f"  TOTAL (replay)                   : "
          f"{ms(t_setup + t_loop + t_write):8.4f} ms/step")
    # device->host syncs per loop iteration: int(c) + int(crowd) + bool(basal)
    # + float(rn) + float(gain) + bool(over) ≈ 5-6 per active cell per step.
    print(f"\n  est. blocking GPU->CPU scalar syncs in loop "
          f"≈ {5 * int(active_ids.shape[0])}/step "
          f"(5 per active cell x {int(active_ids.shape[0])} cells)")


if __name__ == "__main__":
    main()
