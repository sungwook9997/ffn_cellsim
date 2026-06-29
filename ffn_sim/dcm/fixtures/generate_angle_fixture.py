"""Generate + commit the HOOMD reference fixture for the harmonic-angle Warp parity test.

Builds chains of beads with consecutive triplets as angles (vertex = middle bead),
bent away from the rest angle, and runs HOOMD's native ``md.angle.Harmonic`` to get
the per-bead net force + energy — the parity TARGET for the angle parity test.

    python ffn_sim/warp_port/fixtures/generate_angle_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

import hoomd
import hoomd.md as md

HERE = os.path.dirname(os.path.abspath(__file__))

NCHAIN, LCHAIN = 400, 8
N = NCHAIN * LCHAIN
K = 20.0           # bending constant
T0 = 2.6           # rest angle [rad] (~149°), away from π so θ != t0 and not collinear
BOND = 1.0
L = 200.0
SEED = 11


def main() -> None:
    rng = np.random.default_rng(SEED)
    pos = np.zeros((N, 3))
    angles = []
    for c in range(NCHAIN):
        base = c * LCHAIN
        start = (rng.random(3) - 0.5) * (L * 0.2)
        p = start.copy()
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        for b in range(LCHAIN):
            pos[base + b] = p
            # turn the direction by a random moderate angle each step so the
            # chain is bent (θ around ~120-160°, never collinear)
            turn = rng.normal(size=3) * 0.4
            d = d + turn
            d /= np.linalg.norm(d)
            p = p + BOND * d
            if b < LCHAIN - 2:
                angles.append((base + b, base + b + 1, base + b + 2))
    angles = np.array(angles, dtype=np.int32)
    A = angles.shape[0]

    snap = hoomd.Snapshot()
    snap.configuration.box = [L, L, L, 0, 0, 0]
    snap.particles.N = N
    snap.particles.types = ["bead"]
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = np.zeros(N, dtype=np.int32)
    snap.angles.N = A
    snap.angles.types = ["A"]
    snap.angles.typeid[:] = np.zeros(A, dtype=np.int32)
    snap.angles.group[:] = angles

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.tuners.clear()
    integrator = md.Integrator(dt=1e-4)
    angle = md.angle.Harmonic()
    angle.params["A"] = dict(k=K, t0=T0)
    integrator.forces.append(angle)
    sim.operations.integrator = integrator
    sim.run(0)

    F = np.asarray(angle.forces, dtype=np.float64).copy()
    U = np.asarray(angle.energies, dtype=np.float64).copy()

    out = os.path.join(HERE, "network_angle_ref.npz")
    np.savez(
        out,
        N=N, A=A, k=K, t0=T0, L=L,
        pos=pos, angles=angles, ref_force=F, ref_energy=U,
    )
    print(f"wrote {out}  N={N} A={A}  max|F|={np.abs(F).max():.4e}  "
          f"sum|U|={U.sum():.4e}")


if __name__ == "__main__":
    main()
