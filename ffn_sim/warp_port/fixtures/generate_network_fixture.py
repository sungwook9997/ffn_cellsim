"""Generate + commit the HOOMD reference fixture for the harmonic-bond Warp parity test.

Builds chains of beads (consecutive bonds), perturbs them off rest length, and runs
HOOMD's native ``md.bond.Harmonic`` to get the per-bead net force + energy — the
bit/tolerance-parity TARGET for ``tests/warp_port/test_network_warp_parity.py``.

    python ffn_sim/warp_port/fixtures/generate_network_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

import hoomd
import hoomd.md as md

HERE = os.path.dirname(os.path.abspath(__file__))

NCHAIN, LCHAIN = 400, 8     # 400 chains of 8 beads -> 3200 beads, 2800 bonds
N = NCHAIN * LCHAIN
K = 50.0                    # spring constant
R0 = 1.0                    # rest length
L = 200.0                   # box (>> chain extent → bonds never wrap)
SEED = 7


def main() -> None:
    rng = np.random.default_rng(SEED)
    pos = np.zeros((N, 3))
    bonds = []
    for c in range(NCHAIN):
        base = c * LCHAIN
        start = (rng.random(3) - 0.5) * (L * 0.2)
        p = start.copy()
        for b in range(LCHAIN):
            pos[base + b] = p
            # next bead ~R0 away in a random direction, perturbed ±20% so r != r0
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            p = p + R0 * (1.0 + 0.2 * rng.normal()) * d
            if b < LCHAIN - 1:
                bonds.append((base + b, base + b + 1))
    bonds = np.array(bonds, dtype=np.int32)
    B = bonds.shape[0]

    snap = hoomd.Snapshot()
    snap.configuration.box = [L, L, L, 0, 0, 0]
    snap.particles.N = N
    snap.particles.types = ["bead"]
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = np.zeros(N, dtype=np.int32)
    snap.bonds.N = B
    snap.bonds.types = ["A"]
    snap.bonds.typeid[:] = np.zeros(B, dtype=np.int32)
    snap.bonds.group[:] = bonds

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.tuners.clear()
    integrator = md.Integrator(dt=1e-4)
    bond = md.bond.Harmonic()
    bond.params["A"] = dict(k=K, r0=R0)
    integrator.forces.append(bond)
    sim.operations.integrator = integrator
    sim.run(0)

    F = np.asarray(bond.forces, dtype=np.float64).copy()
    U = np.asarray(bond.energies, dtype=np.float64).copy()

    out = os.path.join(HERE, "network_bond_ref.npz")
    np.savez(
        out,
        N=N, B=B, k=K, r0=R0, L=L,
        pos=pos, bonds=bonds, ref_force=F, ref_energy=U,
    )
    print(f"wrote {out}  N={N} B={B}  max|F|={np.abs(F).max():.4e}  "
          f"sum|U|={U.sum():.4e}")


if __name__ == "__main__":
    main()
