"""Parity + host-sync-removal for ALL 3 compartment-force GPU ports (Stage 2).

For each of NucleusConfinement / MembraneSurfaceTension / EnclosedVolumePressure:
builds a shell-bead cloud, attaches the original (cpu_local_snapshot) force AND its
device-resident GPU twin, compares per-particle force/energy (must match to ~eps),
and times set_forces (the per-step host-sync the twin removes). Read-only.

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_compartment_force_parity --n 12840
"""

from __future__ import annotations

import argparse
import time

import hoomd
import hoomd.md as md
import numpy as np

from ffn_sim.cell.nucleus import resolve_nucleus, NucleusConfinement
from ffn_sim.cell.nucleus_confinement_gpu import NucleusConfinementGPU
from ffn_sim.cell.membrane_surface import resolve_membrane_surface, MembraneSurfaceTension
from ffn_sim.cell.membrane_surface_gpu import MembraneSurfaceTensionGPU
from ffn_sim.cortex.enclosed_volume import resolve_enclosed_volume, EnclosedVolumePressure
from ffn_sim.cortex.enclosed_volume_gpu import EnclosedVolumePressureGPU


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


def _one(name, cpu_force, gpu_force, pos, steps):
    sim = _sim(pos)
    Fc, Uc = _force_of(sim, cpu_force)
    Fg, Ug = _force_of(sim, gpu_force)
    dF = float(np.max(np.abs(Fc - Fg))); sF = float(np.max(np.abs(Fc))) + 1e-300
    dU = float(np.max(np.abs(Uc - Ug))); sU = float(np.max(np.abs(Uc))) + 1e-300
    ok = (dF < 1e-10 * sF) and (dU < 1e-10 * sU)
    t_cpu = _time(sim, cpu_force, steps)
    t_gpu = _time(sim, gpu_force, steps)
    print(f"[{name:9s}] force rel={dF/sF:.1e} energy rel={dU/sU:.1e}  "
          f"{'PASS' if ok else 'FAIL'}  |  cpu {t_cpu:7.1f}us -> gpu {t_gpu:7.1f}us "
          f"({t_cpu/max(t_gpu,1e-9):.2f}x, -{t_cpu-t_gpu:.0f}us)", flush=True)
    return ok, t_cpu - t_gpu


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12840)
    ap.add_argument("--steps", type=int, default=500)
    ap.add_argument("--seed", type=int, default=3)
    args = ap.parse_args()
    n, R_cell = args.n, 7.5e-6
    R_nuc = 3.0e-6
    rng = (0, n)

    p_nuc = resolve_nucleus(
        {"nucleus": {"E_nuc": 5e3, "ratio_lamin": 3.0, "knee_strain": 0.10,
                     "critical_pore_area_um2": 7.0}}, R_nuc=R_nuc, n_beads=n)
    p_mem = resolve_membrane_surface(
        {"membrane_surface": {"gamma_mem": 1.0e-4}}, R_cell=R_cell)  # T_m in KU-3.B1 band
    p_ev = resolve_enclosed_volume({}, R_cell=R_cell)

    cloud_nuc = _cloud(n, R_nuc, 3.0 * p_nuc.d_knee, args.seed)
    cloud_shell = _cloud(n, R_cell, 0.05 * R_cell, args.seed)  # V != V0 -> nonzero elastic

    results, saved = [], 0.0
    for name, cf, gf, pos in (
        ("nucleus", NucleusConfinement(p_nuc, rng), NucleusConfinementGPU(p_nuc, rng), cloud_nuc),
        ("membrane", MembraneSurfaceTension(p_mem, rng), MembraneSurfaceTensionGPU(p_mem, rng), cloud_shell),
        ("turgor", EnclosedVolumePressure(p_ev, rng), EnclosedVolumePressureGPU(p_ev, rng), cloud_shell),
    ):
        ok, dt = _one(name, cf, gf, pos, args.steps)
        results.append(ok); saved += max(dt, 0.0)

    print(f"[total] 3 compartment forces host-sync saved ~= {saved:.0f} us/step "
          f"(full-cell per-step)", flush=True)
    print("compartment-force GPU parity " + ("PASS" if all(results) else "FAIL"), flush=True)
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
