"""Speed: native device-resident constrained BAOAB vs cupy ConstrainedLM-BAOAB.

Constraint chains only (no external forces) → isolates the constrained-integrator
cost (Fixman + predictor + M-SHAKE), which on the cupy path is the launch-bound
38k-tiny-op hotspot the native single-thread-per-chain kernels target.

    PYTHONPATH="$PWD/build:$PWD/python:$HOME/ffn_cellsim" python bench_constrained.py --n-chains 2500
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import hoomd

from ffn_sim.integrator.constrained_baoab import ConstrainedLeimkuhlerMatthewsBAOAB
from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater


def _build(F, m, L, r0, rng):
    N = F * (m + 1)
    pos = np.zeros((N, 3))
    chains, cpairs = [], []
    for f in range(F):
        base = f * (m + 1)
        c = (rng.random(3) - 0.5) * (L * 1e-3)
        ch = []
        for b in range(m + 1):
            pos[base + b] = c
            ch.append(base + b)
            d = rng.normal(size=3); d /= np.linalg.norm(d)
            c = c + r0 * d
        chains.append(np.array(ch, dtype=np.int64))
        for a in range(m):
            cpairs.append((base + a, base + a + 1))
    return N, pos, chains, np.array(cpairs, dtype=np.int64)


def _snap(N, L, pos):
    s = hoomd.Snapshot()
    if s.communicator.rank == 0:
        s.configuration.box = [L, L, L, 0, 0, 0]
        s.particles.N = N
        s.particles.types = ["A"]
        s.particles.position[:] = pos
        s.particles.typeid[:] = np.zeros(N, dtype=np.int32)
    return s


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-chains", type=int, default=2500)
    ap.add_argument("--m", type=int, default=6)
    ap.add_argument("--steps", type=int, default=3000)
    args = ap.parse_args()
    F, m = args.n_chains, args.m
    L, r0, gamma, dt, kT, seed = 1.0e6, 100.0, 1.0, 1.0e-3, 1.0, 7
    dev = hoomd.device.GPU()
    rng = np.random.default_rng(11)
    N, pos, chains, cpairs = _build(F, m, L, r0, rng)
    cl = np.full(cpairs.shape[0], r0)
    chains_flat = np.concatenate(chains).astype(np.int32)

    def mk(native):
        sim = hoomd.Simulation(device=dev, seed=1)
        sim.create_state_from_snapshot(_snap(N, L, pos))
        ig = hoomd.md.Integrator(dt=dt); ig.forces = []; ig.methods = []
        sim.operations.integrator = ig
        if native:
            sim.operations.updaters.append(NativeConstrainedBaoabUpdater(
                dt=dt, kT=kT, inv_gamma_by_tag=[1.0 / gamma] * N, chains_tag=chains_flat,
                n_chains=F, bonds_per_chain=m, rest_length=r0, seed=seed,
                trigger=hoomd.trigger.Periodic(1)))
        else:
            act = ConstrainedLeimkuhlerMatthewsBAOAB(
                kT=kT, gamma={"A": gamma}, dt=dt, constraint_pairs=cpairs,
                constraint_lengths=cl, chains=chains, seed=seed)
            sim.operations.updaters.append(
                hoomd.update.CustomUpdater(action=act, trigger=hoomd.trigger.Periodic(1)))
        sim.run(0); sim.run(200)  # warm
        return sim

    res = {}
    for label, native in (("native", True), ("cupy", False)):
        sim = mk(native)
        t0 = time.perf_counter(); sim.run(args.steps); wall = time.perf_counter() - t0
        res[label] = args.steps / wall
        print(f"[bench] {label:6s} N={N} ({F} chains x {m+1}) steps={args.steps} "
              f"{res[label]:.1f} steps/s ({1e6*wall/args.steps:.1f} us/step)")
    print(f"[bench] constrained-integrator speedup native/cupy = {res['native']/res['cupy']:.2f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
