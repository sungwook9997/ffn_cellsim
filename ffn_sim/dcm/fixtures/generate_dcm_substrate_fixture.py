"""Generate + commit the HOOMD reference fixture for the DCM substrate Warp parity test.

Runs the REAL ``cell.dcm_gpu_forces.DcmSubstrateForceGPU`` (z-well + rigid-dish floor)
and ``DcmSubstrateWettingGPU`` (in-plane area-maximizing wetting) on CPU, on ONE
icosahedral cell straddling the substrate plane z0 so that ALL branches are exercised:

  * z-well: top nodes dz>rng (F=0), mid nodes |dz|<=rng (harmonic), bottom nodes
    dz<-rng (capped push-up) and z<z0 (rigid floor).
  * wetting: low (basal) faces engaged w>0, high faces w=0; non-degenerate xy areas.

The per-node net force of each is the parity TARGET for
``tests/warp_port/test_dcm_substrate_warp_parity.py`` (guard-rail 2: graded vs THIS
committed HOOMD output, not a re-authored oracle).

    python ffn_sim/warp_port/fixtures/generate_dcm_substrate_fixture.py
"""

from __future__ import annotations

import math
import os

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_forces import DcmSubstrateForceGPU, DcmSubstrateWettingGPU

HERE = os.path.dirname(os.path.abspath(__file__))

R = 1.0e-6
Z0 = 0.0                 # substrate plane; cell centred ON it → straddles z0
SEED = 23
# z-well params (DcmSubstrateForceGPU)
W_CS = 2.0e-16           # J / node (well depth)
ADH_RANGE = 0.5e-6       # m (well/contact range)
K_FLOOR = 1.0            # N/m (rigid dish, production default)
# wetting params (DcmSubstrateWettingGPU)
W_CS_JM2 = 2.85e-3       # J/m² (lit MCF7, Gil-Redondo 2023)
FORCE_CAP = 5.0e-8       # N (production default; no node caps at these scales — clean parity)


def _icosahedron(radius, center):
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    verts = np.array([
        (-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
        (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
        (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1),
    ], dtype=np.float64)
    verts /= np.linalg.norm(verts[0])
    verts = verts * radius + np.asarray(center, dtype=np.float64)
    faces = np.array([
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ], dtype=np.int64)
    return verts, faces


def _sim(pos: np.ndarray):
    n = pos.shape[0]
    snap = hoomd.Snapshot()
    L = 20.0e-6
    snap.configuration.box = [L, L, L, 0, 0, 0]
    snap.particles.N = n
    snap.particles.types = ["node"]
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = np.zeros(n, dtype=np.int32)
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.tuners.clear()
    return sim


def main() -> None:
    rng = np.random.default_rng(SEED)
    pos, faces = _icosahedron(R, (0.0, 0.0, Z0))
    pos = pos + rng.uniform(-0.05e-6, 0.05e-6, size=pos.shape)
    N = pos.shape[0]

    # --- z-well reference ---
    sim_w = _sim(pos)
    integ_w = md.Integrator(dt=1e-6)
    well = DcmSubstrateForceGPU(z0=Z0, W_cs=W_CS, adh_range=ADH_RANGE,
                                k_sub=None, k_floor=K_FLOOR)
    integ_w.forces.append(well)
    sim_w.operations.integrator = integ_w
    sim_w.run(0)
    F_well = np.asarray(well.forces, dtype=np.float64).copy()

    # --- wetting reference ---
    sim_t = _sim(pos)
    integ_t = md.Integrator(dt=1e-6)
    wet = DcmSubstrateWettingGPU(faces=faces, z0=Z0, W_cs_Jm2=W_CS_JM2,
                                 adh_range=ADH_RANGE, force_cap=FORCE_CAP)
    integ_t.forces.append(wet)
    sim_t.operations.integrator = integ_t
    sim_t.run(0)
    F_wet = np.asarray(wet.forces, dtype=np.float64).copy()

    out = os.path.join(HERE, "dcm_substrate_ref.npz")
    np.savez(
        out,
        N=N, z0=Z0, W_cs=W_CS, adh_range=ADH_RANGE, k_floor=K_FLOOR,
        W_cs_Jm2=W_CS_JM2, force_cap=FORCE_CAP,
        pos=pos, faces=faces.astype(np.int64),
        ref_force_well=F_well, ref_force_wetting=F_wet,
    )
    # branch-coverage audit (so the fixture is known to exercise all paths)
    dz = pos[:, 2] - Z0
    n_well = int((np.abs(dz) <= ADH_RANGE).sum())
    n_deep = int((dz < -ADH_RANGE).sum())
    n_floor = int((pos[:, 2] < Z0).sum())
    zc = pos[faces].mean(axis=1)[:, 2]
    w = np.clip(1.0 - (zc - Z0) / ADH_RANGE, 0.0, 1.0)
    print(f"wrote {out}  N={N} faces={faces.shape[0]}")
    print(f"  z-well: |dz|<=rng {n_well}, deep(dz<-rng) {n_deep}, floor(z<z0) {n_floor}; "
          f"max|F_well|={np.abs(F_well).max():.4e}")
    print(f"  wetting: faces engaged(w>0) {int((w > 0).sum())}/{faces.shape[0]}, "
          f"w=1 {int((w >= 1).sum())}; max|F_wet|={np.abs(F_wet).max():.4e}")


if __name__ == "__main__":
    main()
