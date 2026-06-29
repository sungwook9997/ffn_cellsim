"""Baseline runtime: the EXISTING single-cell DCM per-step path (for the Warp compare).

Same physics + same icosphere as ``dcm_warp_hybrid`` (radial-shell turgor via
``EnclosedVolumePressure`` + cortex edge springs via ``md.bond.Harmonic`` + the
L-M BAOAB Action), run on the existing HOOMD path. The turgor + BAOAB Action read
``cpu_local_snapshot`` each step (the per-step host touch the Warp hybrid removes),
so this is the apples-to-apples "existing production runtime" number to compare the
Warp GPU-resident per-step throughput against.

    python -m ffn_sim.dcm.bench_dcm_baseline --device cuda:0 --steps 2000 --subdiv 3
"""

from __future__ import annotations

import argparse
import time

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.dcm.geometry import icosphere_mesh, ResolvedDCM
from ffn_sim.archive.hoomd_legacy.cortex.enclosed_volume import resolve_enclosed_volume, EnclosedVolumePressure
from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater


def run_baseline(*, steps: int, subdiv: int = 3, device: str = "cpu",
                 kT: float = 0.0, dt: float = 1.0e-7, warmup: int = 50) -> dict:
    p = ResolvedDCM(subdivisions=subdiv)
    verts, edges, tris = icosphere_mesh(p.R_cell, subdiv)
    N = verts.shape[0]
    r0 = float(np.linalg.norm(verts[edges[:, 0]] - verts[edges[:, 1]], axis=1).mean())

    dev = hoomd.device.GPU(notice_level=0) if device.startswith("cuda") else hoomd.device.CPU(notice_level=0)
    snap = hoomd.Snapshot()
    L = 60.0e-6
    snap.configuration.box = [L, L, L, 0, 0, 0]
    snap.particles.N = N
    snap.particles.types = ["node"]
    snap.particles.position[:] = verts
    snap.particles.typeid[:] = np.zeros(N, dtype=np.int32)
    snap.bonds.N = edges.shape[0]
    snap.bonds.types = ["edge"]
    snap.bonds.typeid[:] = np.zeros(edges.shape[0], dtype=np.int32)
    snap.bonds.group[:] = edges.astype(np.int32)

    sim = hoomd.Simulation(device=dev, seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.tuners.clear()

    integrator = md.Integrator(dt=dt)
    p_ev = resolve_enclosed_volume({"turgor_dP0": p.turgor_dP0}, R_cell=p.R_cell)
    turgor = EnclosedVolumePressure(p_ev, (0, N))
    bond = md.bond.Harmonic()
    bond.params["edge"] = dict(k=p.k_edge, r0=r0)
    integrator.forces = [turgor, bond]
    integrator.methods = []
    sim.operations.integrator = integrator
    _action, updater = make_baoab_updater(kT=kT, gamma={"node": p.gamma_node}, dt=dt)
    sim.operations.updaters.append(updater)

    sim.run(warmup)
    t0 = time.perf_counter()
    sim.run(steps)
    elapsed = time.perf_counter() - t0
    return {"device": device, "N": N, "subdiv": subdiv, "steps": steps,
            "elapsed_s": elapsed, "steps_per_s": steps / elapsed, "dt": dt}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--subdiv", type=int, default=3)
    ap.add_argument("--kT", type=float, default=0.0)
    args = ap.parse_args()
    import json
    print(json.dumps(run_baseline(steps=args.steps, subdiv=args.subdiv,
                                  device=args.device, kT=args.kT), indent=2))


if __name__ == "__main__":
    main()
