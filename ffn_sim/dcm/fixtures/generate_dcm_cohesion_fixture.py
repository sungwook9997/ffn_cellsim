"""Generate + commit the HOOMD reference fixture for the DCM cohesion parity test.

Runs the REAL ``cell.dcm_contact.DcmTentContact`` (node-node inter-cell repulsion +
adhesive tent) in a HOOMD CPU sim on two overlapping clouds of nodes (two cells)
positioned so inter-cell node pairs span BOTH the repulsion (d < r_contact) and the
adhesion (r_contact ≤ d < c_adh) branches. The per-node net force is the parity
TARGET for ``tests/warp_port/test_dcm_cohesion_warp_parity.py``.

    python ffn_sim/warp_port/fixtures/generate_dcm_cohesion_fixture.py
"""

from __future__ import annotations

import os

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.cell.dcm_contact import DcmTentContact

HERE = os.path.dirname(os.path.abspath(__file__))

NPC = 40                      # nodes per cell
R = 0.45e-6                   # cloud radius
SEP = 0.7e-6                  # center separation (< 2R → overlap → rep + adh span)
R_CONTACT = 0.5e-6
C_ADH = 1.0e-6
REP = 2.0e8
ADH = 8.0e7
PATCH_AREA = 1.0e-13
FORCE_CAP = 5.0e-8
SEED = 23


def main() -> None:
    rng = np.random.default_rng(SEED)
    c0 = rng.normal(size=(NPC, 3)); c0 /= np.linalg.norm(c0, axis=1)[:, None]
    c1 = rng.normal(size=(NPC, 3)); c1 /= np.linalg.norm(c1, axis=1)[:, None]
    # radial jitter so distances span the rep + adh ranges
    r0 = R * (0.8 + 0.4 * rng.random(NPC))
    r1 = R * (0.8 + 0.4 * rng.random(NPC))
    v0 = c0 * r0[:, None] + np.array([-SEP / 2.0, 0.0, 0.0])
    v1 = c1 * r1[:, None] + np.array([SEP / 2.0, 0.0, 0.0])
    pos = np.vstack([v0, v1])
    N = pos.shape[0]
    cell_of_node = np.concatenate([
        np.zeros(NPC, np.int64), np.ones(NPC, np.int64)])
    cad_mult = np.array([1.2, 0.7], dtype=np.float64)

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
    coh = DcmTentContact(
        cell_of_node=cell_of_node, r_contact=R_CONTACT, c_adh=C_ADH,
        rep_strength=REP, adh_strength=ADH, patch_area=PATCH_AREA,
        force_cap=FORCE_CAP, cad_mult=cad_mult)
    integrator.forces.append(coh)
    sim.operations.integrator = integrator
    sim.run(0)

    F = np.asarray(coh.forces, dtype=np.float64).copy()

    # branch-coverage audit
    n_rep = n_adh = 0
    for i in range(N):
        for j in range(N):
            if i < j and cell_of_node[i] != cell_of_node[j]:
                d = float(np.linalg.norm(pos[i] - pos[j]))
                if d < R_CONTACT:
                    n_rep += 1
                elif d < C_ADH:
                    n_adh += 1
    assert n_rep > 0 and n_adh > 0, f"branches not covered: rep={n_rep} adh={n_adh}"

    out = os.path.join(HERE, "dcm_cohesion_ref.npz")
    np.savez(
        out,
        N=N, r_contact=R_CONTACT, c_adh=C_ADH, rep=REP, adh=ADH,
        patch_area=PATCH_AREA, force_cap=FORCE_CAP,
        pos=pos, cell_of_node=cell_of_node, cad_mult=cad_mult, ref_force=F,
    )
    print(f"wrote {out}  N={N}  max|F|={np.abs(F).max():.4e}  "
          f"nodes-with-force={int((np.abs(F).sum(1)>0).sum())}  "
          f"branches rep={n_rep} adh={n_adh}")


if __name__ == "__main__":
    main()
