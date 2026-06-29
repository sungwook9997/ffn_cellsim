"""Generate + commit the HOOMD reference fixture for the lamellipodium parity test.

Runs the REAL ``cell.dcm_lamellipodium.LamellipodialTractionTether`` in a HOOMD CPU
sim on a minimal 2-cell basal config: each cell is a ring of membrane nodes with a
ring of ``actin_lamel`` beads just OUTWARD of it, so each cell's leading (outward)
basal nodes tether to a nearest outward actin bead. The per-node net force is the
parity TARGET for ``tests/warp_port/test_lamellipodium_warp_parity.py``.

    python ffn_sim/warp_port/fixtures/generate_lamellipodium_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.cell.dcm_lamellipodium import LamellipodialTractionTether

HERE = os.path.dirname(os.path.abspath(__file__))

N_MEM = 12                    # membrane nodes per cell (ring)
N_ACT = 16                    # actin beads per cell (outward ring)
R_MEM = 1.0e-6
R_ACT = 1.4e-6
Z_BASAL = 0.0
BASAL_BAND = 0.5e-6
K_TETHER = 1.0e-4
FORCE_CAP = 5.0e-8
TETHER_RADIUS = 1.0e-6
LEAD_FRAC = 0.0
CENTERS = np.array([[-2.0e-6, 0.0, 0.0], [2.0e-6, 0.0, 0.0]])
SEED = 29


def _ring(n, r, center, rng):
    th = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    p = np.stack([r * np.cos(th), r * np.sin(th), np.zeros(n)], axis=1) + center
    p += rng.uniform(-0.02e-6, 0.02e-6, size=p.shape)  # break exact ties
    p[:, 2] = Z_BASAL
    return p


def main() -> None:
    rng = np.random.default_rng(SEED)
    mem = np.vstack([_ring(N_MEM, R_MEM, c, rng) for c in CENTERS])
    act = np.vstack([_ring(N_ACT, R_ACT, c, rng) for c in CENTERS])
    pos = np.vstack([mem, act])
    n_mem = mem.shape[0]
    N = pos.shape[0]
    typeid = np.concatenate([
        np.zeros(n_mem, np.int32), np.ones(act.shape[0], np.int32)])
    cell_of_node = np.concatenate([
        np.zeros(N_MEM, np.int64), np.ones(N_MEM, np.int64)])   # membrane block only
    rim_cells = np.array([0, 1], dtype=np.int64)

    snap = hoomd.Snapshot()
    L = 40.0e-6
    snap.configuration.box = [L, L, L, 0, 0, 0]
    snap.particles.N = N
    snap.particles.types = ["mem", "actin"]
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = typeid

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.tuners.clear()
    integrator = md.Integrator(dt=1e-6)
    tether = LamellipodialTractionTether(
        cell_of_node=cell_of_node, rim_cells=rim_cells,
        actin_typeid=1, mem_typeid=0, z_basal=Z_BASAL, basal_band=BASAL_BAND,
        k_tether=K_TETHER, force_cap=FORCE_CAP, tether_radius=TETHER_RADIUS,
        lead_frac=LEAD_FRAC)
    integrator.forces.append(tether)
    sim.operations.integrator = integrator
    sim.run(0)

    F = np.asarray(tether.forces, dtype=np.float64).copy()
    n_teth = int(tether.n_tethered)
    assert n_teth > 0, "no nodes tethered — fixture exercises nothing"

    out = os.path.join(HERE, "lamellipodium_ref.npz")
    np.savez(
        out,
        N=N, n_mem=n_mem, actin_typeid=1, mem_typeid=0,
        z_basal=Z_BASAL, basal_band=BASAL_BAND, k_tether=K_TETHER,
        force_cap=FORCE_CAP, tether_radius=TETHER_RADIUS, lead_frac=LEAD_FRAC,
        pos=pos, typeid=typeid, cell_of_node=cell_of_node, rim_cells=rim_cells,
        ref_force=F,
    )
    print(f"wrote {out}  N={N}  n_tethered={n_teth}  "
          f"max|F|={np.abs(F).max():.4e}  nodes-with-force={int((np.abs(F).sum(1)>0).sum())}")


if __name__ == "__main__":
    main()
