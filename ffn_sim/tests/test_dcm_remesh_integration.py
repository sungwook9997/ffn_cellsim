"""Phase 2 step 2b integration: the live DcmRemeshUpdater holds the mesh contract.

Builds a real single-cell DCM sim on CPU with node-vs-face contact + a dormant node
pool, drives it with the remesh CustomUpdater at low cadence under a spreading
forcing, and asserts that after the run the turgor/contact face arrays are still a
closed manifold, stay length-consistent, the enclosed volume is positive (outward
winding intact), and positions are finite (no LJ-explosion). CPU + short — the GPU
A/A0 spreading validation runs on gbook.
"""

import numpy as np
import pytest

from ffn_sim.cell.dcm_gpu_build import build_gpu_dcm_simulation, ResolvedGpuDCM
from ffn_sim.cell.dcm_remesh import mesh_edges, enclosed_volume
from ffn_sim.cell.dcm_remesh_updater import attach_remesh_updater


def _closed_manifold(faces):
    _e, ef, _ea, bnd = mesh_edges(faces)
    return (not bnd.any()) and (ef >= 0).all()


def test_remesh_updater_preserves_mesh_contract_under_spreading():
    import hoomd
    p = ResolvedGpuDCM()
    h = build_gpu_dcm_simulation(
        p, n_cells=1, device=hoomd.device.CPU(notice_level=0),
        active=True, with_substrate=True, node_face_contact=True, n_pool=40)
    sim = h["sim"]
    turgor, contact = h["turgor"], h["contact"]
    faces0 = np.array(turgor.faces, copy=True)
    V0 = h["V0"]

    # remesh at low cadence; l_min tight enough that spreading stretch triggers SPLITs
    action, _cu = attach_remesh_updater(h, period=250, max_ops=12)

    sim.run(3000)

    faces = np.asarray(turgor.faces)
    # contract: turgor and contact share the SAME topology
    assert np.array_equal(np.asarray(contact.faces), faces)
    assert np.array_equal(np.asarray(contact.face_cell), np.asarray(turgor.face_cell))
    assert faces.shape[0] == np.asarray(turgor.face_cell).shape[0]
    # mesh stayed a closed sphere
    assert _closed_manifold(faces)
    assert np.unique(faces).size - mesh_edges(faces)[0].shape[0] + faces.shape[0] == 2
    # positions finite (no explosion) and turgor volume positive (outward winding)
    with sim.state.cpu_local_snapshot as snap:
        pos = np.asarray(snap.particles.position).copy()
        tag = np.asarray(snap.particles.tag).copy()
    assert np.isfinite(pos).all()
    pos_g = pos[np.argsort(tag)]
    assert enclosed_volume(pos_g, faces) > 0
    # the volume should not have collapsed or blown up (turgor holds it near V0)
    assert 0.2 * V0 < enclosed_volume(pos_g, faces) < 5.0 * V0
    # the updater actually ran (acts fired); ops may or may not have triggered
    assert action.n_acts >= 1
