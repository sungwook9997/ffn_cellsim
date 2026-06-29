"""Generate + commit the HOOMD reference fixture for the LJ-excluded-volume parity test.

Builds a dense, BOND-FREE bead cluster on a cubic lattice (spacing ≈ σ so the WCA
repulsion is active) with small jitter, and runs HOOMD's native ``md.pair.LJ`` with
``mode='shift'`` and ``r_cut = 2^(1/6) σ`` (the cortex excluded-volume convention,
see ``cortex/cortex.py``) to get the per-particle net force + energy — the
tolerance-parity TARGET for ``tests/dcm/test_network_warp_parity.py``.

Bond-free on purpose: with no bonds the nlist has no bonded-pair exclusions, so
"every pair within r_cut" is unambiguously what HOOMD computes and the all-pairs
Warp kernel is exact. (Bonded-pair exclusion is an nlist/topology concern, separate
from the WCA force law ported here.) Lattice spacing < r_cut puts the 6 axial
neighbours inside the cutoff; the jitter pushes some pairs just outside, exercising
the cutoff boundary; beads fill the periodic box, exercising min-image wrap.

    python ffn_sim/warp_port/fixtures/generate_lj_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

import hoomd
import hoomd.md as md

HERE = os.path.dirname(os.path.abspath(__file__))

NSIDE = 8                       # 8³ = 512 beads
A = 1.0                         # lattice spacing (≈ σ → WCA active)
N = NSIDE ** 3
EPS = 1.0                       # ε_LJ
SIGMA = 1.0                     # σ_LJ
RCUT = 2.0 ** (1.0 / 6.0) * SIGMA   # WCA cutoff at the LJ minimum
L = NSIDE * A                   # periodic box; L > 2·r_cut so min-image is unique
JITTER = 0.06                   # ± lattice jitter (keeps min separation ~0.88 σ)
SEED = 11


def main() -> None:
    rng = np.random.default_rng(SEED)
    grid = np.arange(NSIDE) * A
    gx, gy, gz = np.meshgrid(grid, grid, grid, indexing="ij")
    pos = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1).astype(np.float64)
    pos += (rng.random((N, 3)) - 0.5) * (2.0 * JITTER)
    pos = np.mod(pos, L)                 # wrap into [0, L)
    pos -= L * (pos > 0.5 * L)           # fold to [-L/2, L/2) (HOOMD box convention)

    snap = hoomd.Snapshot()
    snap.configuration.box = [L, L, L, 0, 0, 0]
    snap.particles.N = N
    snap.particles.types = ["A"]
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = np.zeros(N, dtype=np.int32)
    snap.bonds.N = 0                     # bond-free → no nlist exclusions

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.tuners.clear()
    integrator = md.Integrator(dt=1e-4)
    nlist = md.nlist.Tree(buffer=0.5 * SIGMA)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    lj.params[("A", "A")] = dict(epsilon=EPS, sigma=SIGMA)
    lj.r_cut[("A", "A")] = RCUT
    lj.mode = "shift"
    integrator.forces.append(lj)
    sim.operations.integrator = integrator
    sim.run(0)

    # Read back the wrapped positions HOOMD actually used (so the Warp kernel sees
    # the identical configuration), then its per-particle force + energy.
    snap_out = sim.state.get_snapshot()
    pos_used = np.asarray(snap_out.particles.position, dtype=np.float64).copy()
    F = np.asarray(lj.forces, dtype=np.float64).copy()
    U = np.asarray(lj.energies, dtype=np.float64).copy()

    out = os.path.join(HERE, "network_lj_ref.npz")
    np.savez(
        out,
        N=N, epsilon=EPS, sigma=SIGMA, r_cut=RCUT, L=L,
        pos=pos_used, ref_force=F, ref_energy=U,
    )
    nnz = int((np.abs(F).sum(axis=1) > 0).sum())
    print(f"wrote {out}  N={N}  max|F|={np.abs(F).max():.4e}  "
          f"sum_U={U.sum():.4e}  beads_with_force={nnz}")


if __name__ == "__main__":
    main()
