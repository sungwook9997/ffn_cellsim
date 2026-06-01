"""Controlled N-scaling benchmark for the cupy-ported constrained BAOAB Action.

Validation-ladder step 4 (GPU_MAIN_PORT_2026-05-31.md §3b): "re-profile GPU
steps/s at mesoscale AND a larger N to confirm the sync is gone and GPU scales
with N." The cortex demo build packs filaments and needs long LJ-overlap
equilibration that breaks at large n_fil — so this benchmark isolates the thing
that was actually ported (the constrained Action: BAOAB predictor + M-SHAKE +
Fixman) on a clean grid of NON-overlapping rigid filaments with NO pair force
(hence no neighbour list, no LJ-overlap blow-up). Each filament is n_beads with
rigid backbone bonds + a soft harmonic angle, exactly the per-step solver work
the port touches.

Reports CPU vs GPU steps/s as the filament count F grows. The crossover (if
any) is where the per-step cupy kernel-launch latency is amortised by enough
data that the GPU's parallelism wins — the regime the native KU-3.5 run needs.

Run on gbook (A5000):
    ~/miniconda3/envs/ffn_sim/bin/python -m ffn_sim.scripts.gpu_scaling_constrained
"""
from __future__ import annotations

import argparse
import time

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.integrator.constrained_baoab import make_constrained_baoab_updater


def _make_device(kind: str):
    if kind == "gpu":
        return hoomd.device.GPU(notice_level=0)
    return hoomd.device.CPU(notice_level=0)


def build(device_kind: str, *, n_fil: int, n_beads: int = 10, L0: float = 1.0,
          kT: float = 1.0, gamma: float = 1.0, kth: float = 2.0, dt: float = 1e-3,
          seed: int = 7):
    N = n_fil * n_beads
    # Lay filaments out on a cubic grid, well separated (no pair force anyway,
    # but keeps them visually disjoint and the box large vs each filament).
    side = int(np.ceil(n_fil ** (1.0 / 3.0)))
    spacing = (n_beads + 4) * L0
    box = spacing * (side + 2)
    pos = np.zeros((N, 3))
    for f in range(n_fil):
        ix, iy, iz = f % side, (f // side) % side, f // (side * side)
        origin = np.array([ix, iy, iz], dtype=float) * spacing
        for k in range(n_beads):
            pos[f * n_beads + k] = origin + np.array([k * L0, 0.0, 0.0])
    pos -= pos.mean(0)

    snap = hoomd.Snapshot()
    snap.particles.N = N
    snap.particles.types = ["A"]
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = 0
    snap.configuration.box = [box, box, box, 0, 0, 0]
    snap.bonds.N = n_fil * (n_beads - 1)
    snap.bonds.types = ["b"]
    bg, ag = [], []
    for f in range(n_fil):
        base = f * n_beads
        for k in range(n_beads - 1):
            bg.append([base + k, base + k + 1])
        for k in range(n_beads - 2):
            ag.append([base + k, base + k + 1, base + k + 2])
    snap.bonds.group[:] = bg
    snap.bonds.typeid[:] = 0
    snap.angles.N = n_fil * (n_beads - 2)
    snap.angles.types = ["a"]
    snap.angles.group[:] = ag
    snap.angles.typeid[:] = 0

    sim = hoomd.Simulation(device=_make_device(device_kind), seed=1)
    sim.create_state_from_snapshot(snap)
    angle = md.angle.Harmonic()
    angle.params["a"] = dict(k=kth, t0=np.pi)
    sim.operations.integrator = md.Integrator(dt=dt, forces=[angle], methods=[])

    pairs = np.array(bg, dtype=np.int64)
    chains = [np.arange(f * n_beads, (f + 1) * n_beads) for f in range(n_fil)]
    act, upd = make_constrained_baoab_updater(
        kT=kT, gamma={"A": gamma}, dt=dt,
        constraint_pairs=pairs,
        constraint_lengths=np.full(pairs.shape[0], L0),
        chains=chains, seed=seed,
    )
    sim.operations.updaters.append(upd)
    sim.run(0)
    return sim, act, N


def timed(sim, steps):
    sim.run(steps)  # warm the step path / JIT cupy kernels
    t0 = time.time()
    sim.run(steps)
    return steps / (time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fils", type=int, nargs="+",
                    default=[50, 200, 1000, 5000])
    ap.add_argument("--n-beads", type=int, default=10)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--devices", nargs="+", default=["cpu", "gpu"])
    args = ap.parse_args()

    print(f"=== constrained-Action N-scaling (n_beads={args.n_beads}, "
          f"steps={args.steps}, hoomd_gpu_build={hoomd.version.gpu_enabled}) ===",
          flush=True)
    print(f"{'n_fil':>7} {'N_part':>8} {'CPU s/s':>10} {'GPU s/s':>10} "
          f"{'GPU/CPU':>8} {'drift_cpu':>11} {'drift_gpu':>11}", flush=True)

    for nf in args.fils:
        row = {"cpu": (None, None), "gpu": (None, None)}
        npart = None
        for dev in args.devices:
            if dev == "gpu":
                try:
                    import cupy  # noqa: F401
                except Exception as exc:
                    print(f"  [GPU skip] {exc}", flush=True)
                    continue
            sim, act, npart = build(dev, n_fil=nf, n_beads=args.n_beads)
            sps = timed(sim, args.steps)
            row[dev] = (sps, act.max_constraint_drift)
            del sim, act
        cps, cdr = row["cpu"]
        gps, gdr = row["gpu"]
        ratio = (gps / cps) if (cps and gps) else float("nan")
        print(f"{nf:>7} {npart:>8} "
              f"{(cps or float('nan')):>10.1f} {(gps or float('nan')):>10.1f} "
              f"{ratio:>8.2f} {(cdr if cdr is not None else float('nan')):>11.2e} "
              f"{(gdr if gdr is not None else float('nan')):>11.2e}", flush=True)


if __name__ == "__main__":
    main()
