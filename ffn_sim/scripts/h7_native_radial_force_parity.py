"""Parity + speed: native C++ FFNRadialShellForce vs the cpu compartment forces.

The cupy *_gpu twins removed the cpu_local_snapshot host-sync but stayed launch-
bound at the ~233 us/force gpu_local_snapshot floor (h7_native_force_diag). The
native ForceCompute (ArrayHandle, no snapshot) should both MATCH the physics
(modulo atomic-reduction order) and drop to ~tens of us. This compares all 3 laws
vs their originals (force + energy) and times set_forces. Read-only.

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_radial_force_parity --n 12840
"""

from __future__ import annotations

import argparse
import time

import hoomd
import hoomd.md as md
import numpy as np

from ffn_sim.cell.nucleus import resolve_nucleus, NucleusConfinement
from ffn_sim.cell.membrane_surface import resolve_membrane_surface, MembraneSurfaceTension
from ffn_sim.cortex.enclosed_volume import resolve_enclosed_volume, EnclosedVolumePressure
from ffn_hoomd_plugin import NativeRadialShellForce


def _cloud(n, R, spread, seed):
    rng = np.random.default_rng(seed)
    r = R + (rng.random(n) - 0.5) * 2.0 * spread
    u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    return (r[:, None] * u).astype(np.float64)


def _sim(pos):
    n = pos.shape[0]
    snap = hoomd.Snapshot()
    if snap.communicator.rank == 0:
        L = 60.0e-6
        snap.configuration.box = [L, L, L, 0, 0, 0]
        snap.particles.N = n
        snap.particles.types = ["shell"]
        snap.particles.position[:] = pos
        snap.particles.typeid[:] = np.zeros(n, dtype=np.int32)
    sim = hoomd.Simulation(device=hoomd.device.GPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.integrator = md.Integrator(dt=1e-6)
    return sim


def _force_of(sim, force):
    sim.operations.integrator.forces = [force]
    sim.run(0)
    return (np.asarray(force.forces, np.float64).copy(),
            np.asarray(force.energies, np.float64).copy())


def _time(sim, force, steps):
    sim.operations.integrator.forces = [force]
    sim.run(40)
    t0 = time.perf_counter(); sim.run(steps)
    return 1e6 * (time.perf_counter() - t0) / steps


def _one(name, cpu_force, nat_force, pos, steps):
    sim = _sim(pos)
    Fc, Uc = _force_of(sim, cpu_force)
    Fn, Un = _force_of(sim, nat_force)
    dF = float(np.max(np.abs(Fc - Fn))); sF = float(np.max(np.abs(Fc))) + 1e-300
    dU = float(np.max(np.abs(Uc - Un))); sU = float(np.max(np.abs(Uc))) + 1e-300
    ok = (dF < 1e-8 * sF) and (dU < 1e-8 * sU)
    t_cpu = _time(sim, cpu_force, steps)
    t_nat = _time(sim, nat_force, steps)
    print(f"[{name:9s}] force rel={dF/sF:.1e} energy rel={dU/sU:.1e}  {'PASS' if ok else 'FAIL'}"
          f"  |  cpu {t_cpu:7.1f}us -> native {t_nat:6.1f}us ({t_cpu/max(t_nat,1e-9):.1f}x)", flush=True)
    return ok, t_nat


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12840)
    ap.add_argument("--steps", type=int, default=500)
    ap.add_argument("--seed", type=int, default=3)
    args = ap.parse_args()
    n, R_cell, R_nuc, rng = args.n, 7.5e-6, 3.0e-6, (0, args.n)
    p_nuc = resolve_nucleus(
        {"nucleus": {"E_nuc": 5e3, "ratio_lamin": 3.0, "knee_strain": 0.10,
                     "critical_pore_area_um2": 7.0}}, R_nuc=R_nuc, n_beads=n)
    p_mem = resolve_membrane_surface({"membrane_surface": {"gamma_mem": 1.0e-4}}, R_cell=R_cell)
    p_ev = resolve_enclosed_volume({}, R_cell=R_cell)
    cl_nuc = _cloud(n, R_nuc, 3.0 * p_nuc.d_knee, args.seed)
    cl_shell = _cloud(n, R_cell, 0.05 * R_cell, args.seed)

    results, tnat = [], 0.0
    for name, cf, nf, pos in (
        ("nucleus", NucleusConfinement(p_nuc, rng), NativeRadialShellForce.from_nucleus(p_nuc, rng), cl_nuc),
        ("membrane", MembraneSurfaceTension(p_mem, rng), NativeRadialShellForce.from_membrane(p_mem, rng), cl_shell),
        ("turgor", EnclosedVolumePressure(p_ev, rng), NativeRadialShellForce.from_turgor(p_ev, rng), cl_shell),
    ):
        ok, t = _one(name, cf, nf, pos, args.steps)
        results.append(ok); tnat += t
    print(f"[total] 3 native compartment forces ~= {tnat:.0f} us/step (was ~1860 us cupy / ~3700 us cpu)",
          flush=True)
    print("native radial-shell force parity " + ("PASS" if all(results) else "FAIL"), flush=True)
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
