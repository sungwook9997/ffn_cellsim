"""Generate + commit the HOOMD reference fixture for the DCM turgor Warp parity test.

Runs the REAL ``cell.dcm.DcmTurgorForce`` (exact divergence-theorem volume → osmotic
face-normal pressure) in a HOOMD CPU sim on a 2-cell triangulated mesh (two
icosahedra), perturbed so each cell's volume deviates from V0 → non-zero ΔP → non-
zero force. The per-node net force is the parity TARGET for
``tests/dcm/test_dcm_turgor_warp_parity.py``.

    python ffn_sim/warp_port/fixtures/generate_dcm_turgor_fixture.py
"""

from __future__ import annotations

import math
import os

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.cell.dcm import DcmTurgorForce

HERE = os.path.dirname(os.path.abspath(__file__))

R = 1.0e-6
SEP = 3.0e-6
TURGOR_DP0 = 40.0
K_VOL = 1.0e3
SEED = 19


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


def _exact_volume(pos, faces):
    v0, v1, v2 = pos[faces[:, 0]], pos[faces[:, 1]], pos[faces[:, 2]]
    cross = np.cross(v1 - v0, v2 - v0)
    return float(np.einsum("ij,ij->i", v0, cross).sum() / 6.0)


def main() -> None:
    rng = np.random.default_rng(SEED)
    v0, f0 = _icosahedron(R, (-SEP / 2.0, 0.0, 0.0))
    v1, f1 = _icosahedron(R, (SEP / 2.0, 0.0, 0.0))
    nv = v0.shape[0]
    V0 = abs(_exact_volume(v0, f0))          # rest volume = unperturbed cell volume

    pos = np.vstack([v0, v1])
    pos += rng.uniform(-0.08e-6, 0.08e-6, size=pos.shape)   # perturb → Vc ≠ V0
    faces = np.vstack([f0, f1 + nv]).astype(np.int64)
    face_cell = np.concatenate([
        np.zeros(f0.shape[0], np.int64), np.ones(f1.shape[0], np.int64)])
    N = pos.shape[0]

    snap = hoomd.Snapshot()
    L = 20.0e-6
    snap.configuration.box = [L, L, L, 0, 0, 0]
    snap.particles.N = N
    snap.particles.types = ["node"]
    snap.particles.position[:] = pos
    snap.particles.typeid[:] = np.zeros(N, dtype=np.int32)

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)
    sim.operations.tuners.clear()
    integrator = md.Integrator(dt=1e-6)
    turgor = DcmTurgorForce(
        faces=faces, face_cell=face_cell, n_cells=2,
        V0=V0, turgor_dP0=TURGOR_DP0, K_vol=K_VOL)
    integrator.forces.append(turgor)
    sim.operations.integrator = integrator
    sim.run(0)

    F = np.asarray(turgor.forces, dtype=np.float64).copy()

    out = os.path.join(HERE, "dcm_turgor_ref.npz")
    np.savez(
        out,
        N=N, n_cells=2, V0=V0, turgor_dP0=TURGOR_DP0, K_vol=K_VOL,
        pos=pos, faces=faces, face_cell=face_cell, ref_force=F,
    )
    print(f"wrote {out}  N={N} faces={faces.shape[0]}  V0={V0:.4e}  "
          f"max|F|={np.abs(F).max():.4e}  nodes-with-force={int((np.abs(F).sum(1)>0).sum())}")


if __name__ == "__main__":
    main()
