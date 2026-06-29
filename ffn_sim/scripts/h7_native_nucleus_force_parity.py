"""Parity + host-sync-removal for NucleusConfinementGPU (compartment-force port).

Builds a nucleus-bead cloud spanning both bilinear regimes, attaches the original
cpu_local_snapshot NucleusConfinement AND the device-resident NucleusConfinementGPU,
and (1) compares their per-particle force/energy (must match to ~float64 eps modulo
the centroid reduction order), (2) times set_forces for each to show the per-step
host-sync the GPU twin removes. Read-only; wires nothing.

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_nucleus_force_parity --n 4000 --steps 600
"""

from __future__ import annotations

import argparse
import math
import time

import hoomd
import hoomd.md as md
import numpy as np

from ffn_sim.archive.hoomd_legacy.cell.nucleus import resolve_nucleus, NucleusConfinement
from ffn_sim.archive.hoomd_legacy.cell.nucleus_confinement_gpu import NucleusConfinementGPU


def _cloud(n, R_nuc, d_knee, seed):
    rng = np.random.default_rng(seed)
    # radii spanning interior (|d|<d_knee) and exterior (|d|>d_knee) regimes
    r = R_nuc + (rng.random(n) - 0.5) * 6.0 * d_knee
    u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    return (r[:, None] * u).astype(np.float64)


def _sim(pos, dt=1e-6):
    n = pos.shape[0]
    snap = hoomd.Snapshot()
    if snap.communicator.rank == 0:
        L = 40.0e-6
        snap.configuration.box = [L, L, L, 0, 0, 0]
        snap.particles.N = n
        snap.particles.types = ["nucleus_bead"]
        snap.particles.position[:] = pos
        snap.particles.typeid[:] = np.zeros(n, dtype=np.int32)
    sim = hoomd.Simulation(device=hoomd.device.GPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.integrator = md.Integrator(dt=dt)
    return sim


def _force_of(sim, force):
    sim.operations.integrator.forces = [force]
    sim.run(0)
    return np.asarray(force.forces, np.float64).copy(), np.asarray(force.energies, np.float64).copy()


def _time(sim, force, steps):
    sim.operations.integrator.forces = [force]
    sim.run(40)
    t0 = time.perf_counter(); sim.run(steps)
    return 1e6 * (time.perf_counter() - t0) / steps


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=4000)
    ap.add_argument("--steps", type=int, default=600)
    ap.add_argument("--seed", type=int, default=3)
    args = ap.parse_args()

    cfg = {"nucleus": {"E_nuc": 5.0e3, "ratio_lamin": 3.0, "knee_strain": 0.10,
                       "critical_pore_area_um2": 7.0}}
    R_nuc = 3.0e-6
    p = resolve_nucleus(cfg, R_nuc=R_nuc, n_beads=args.n)
    pos = _cloud(args.n, R_nuc, p.d_knee, args.seed)
    rng_tags = (0, args.n)

    sim = _sim(pos)
    f_cpu = NucleusConfinement(p, rng_tags)
    f_gpu = NucleusConfinementGPU(p, rng_tags)
    Fc, Uc = _force_of(sim, f_cpu)
    Fg, Ug = _force_of(sim, f_gpu)

    dF = float(np.max(np.abs(Fc - Fg)))
    sF = float(np.max(np.abs(Fc))) + 1e-300
    dU = float(np.max(np.abs(Uc - Ug)))
    sU = float(np.max(np.abs(Uc))) + 1e-300
    okF = dF < 1e-10 * sF
    okU = dU < 1e-10 * sU
    print(f"nucleus-force parity: n={args.n}", flush=True)
    print(f"[force]  max|cpu-gpu|={dF:.3e} rel={dF/sF:.2e}  {'PASS' if okF else 'FAIL'}", flush=True)
    print(f"[energy] max|cpu-gpu|={dU:.3e} rel={dU/sU:.2e}  {'PASS' if okU else 'FAIL'}", flush=True)

    t_cpu = _time(sim, f_cpu, args.steps)
    t_gpu = _time(sim, f_gpu, args.steps)
    print(f"[set_forces] cpu_local_snapshot {t_cpu:7.1f} us/step  ->  "
          f"gpu_local (cupy) {t_gpu:7.1f} us/step  ({t_cpu/max(t_gpu,1e-9):.2f}x, "
          f"saves {t_cpu-t_gpu:.1f} us host-sync)", flush=True)
    ok = okF and okU
    print("NucleusConfinementGPU parity " + ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
