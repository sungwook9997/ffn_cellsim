"""DcmActiveMotilitySPP (Stage-1 aggregation self-propulsion) — correctness + parity.

On a CUDA-less machine only the CPU dispatch path is exercisable; this pins the CPU
path (the GPU parity reference) and the active-matter-aggregation invariants the PI
requires:

  1. NET force on each live cell == f_active · p_c  (the self-propulsion is per-cell,
     distributed over its nodes, NOT toward any centre).
  2. EVERY live cell gets a non-zero active force — no cell is silently frozen
     (the all-cells-full-physics rule; contrast the rim traction which excludes
     interior cells).
  3. Dormant / non-membrane nodes get zero force.
  4. The polarity updater keeps |p_c| == 1 exactly (rotational diffusion, no drift).

Run:  ~/miniconda3/envs/ffn_sim/bin/python -m pytest \
          ffn_sim/tests/test_dcm_motility_gpu_parity.py -q
"""

from __future__ import annotations

import numpy as np

from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_build import ResolvedGpuDCM, build_gpu_dcm_simulation
from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_forces import (
    DcmActiveMotilitySPP, DcmPolarityUpdater, on_gpu)


def _force_local(force):
    with force.cpu_local_force_arrays as arr:
        return np.asarray(arr.force).copy()


def _build(n=12):
    p = ResolvedGpuDCM(subdivisions=1, seed=7)
    h = build_gpu_dcm_simulation(p, n, active=False, with_substrate=False,
                                 settle_force=0.0)
    return p, h


def test_cpu_path_selected():
    _, h = _build()
    assert not on_gpu(h["sim"])


def test_net_force_per_cell_is_factive_times_polarity():
    """Sum of the motility force over each cell's nodes == f_active · p_c."""
    p, h = _build(12)
    sim, ranges = h["sim"], h["ranges"]
    f_active = 1.6e-10
    mot = DcmActiveMotilitySPP(cell_of_node=h["cell_of_node"], ranges=ranges,
                               n_cells=h["n_cells"], f_active=f_active, f_cap=5.0e-9,
                               ramp_steps=1, seed=7)
    sim.operations.integrator.forces.append(mot)
    sim.run(0)
    mot.set_forces(100)                          # past the 1-step ramp
    F = _force_local(mot)                         # (N,3) LOCAL order
    # map local node -> cell via the LOCAL tags (cell_of_node is global/tag order)
    with sim.state.cpu_local_snapshot as lsnap:
        tag = np.asarray(lsnap.particles.tag)
    cell_of_local = np.asarray(h["cell_of_node"])[tag]
    nv = ranges[0][1] - ranges[0][0]
    nonzero_cells = 0
    for c in range(h["n_cells"]):
        mask = cell_of_local == c
        if not mask.any():
            continue
        net = F[mask].sum(axis=0)
        expect = f_active * mot.p[c]              # NET self-propulsion
        assert np.allclose(net, expect, atol=1e-18), \
            f"cell {c}: net {net} != f_active*p {expect}"
        # per-node share is uniform f_active/nv * p_c (within cap, never hit here)
        assert np.allclose(F[mask], (f_active / nv) * mot.p[c], atol=1e-18)
        if np.linalg.norm(net) > 0:
            nonzero_cells += 1
    # ALL cells full physics: every live cell propelled, none frozen
    assert nonzero_cells == h["n_cells"], \
        f"only {nonzero_cells}/{h['n_cells']} cells propelled — a cell is frozen"
    assert np.all(np.isfinite(F))


def test_dormant_nodes_zero():
    """Non-membrane / dormant nodes get zero active force (typeid mask)."""
    p, h = _build(12)
    sim = h["sim"]
    mot = DcmActiveMotilitySPP(cell_of_node=h["cell_of_node"], ranges=h["ranges"],
                               n_cells=h["n_cells"], f_active=1.6e-10, ramp_steps=1,
                               seed=7)
    sim.operations.integrator.forces.append(mot)
    sim.run(0)
    mot.set_forces(100)
    F = _force_local(mot)
    with sim.state.cpu_local_snapshot as lsnap:
        tid = np.asarray(lsnap.particles.typeid)
    if (tid != 0).any():
        assert np.allclose(F[tid != 0], 0.0)


def test_polarity_updater_unit_norm():
    """The OU reorientation keeps every polarity a unit vector (no drift)."""
    p, h = _build(12)
    mot = DcmActiveMotilitySPP(cell_of_node=h["cell_of_node"], ranges=h["ranges"],
                               n_cells=h["n_cells"], f_active=1.6e-10, ramp_steps=1,
                               seed=7)
    assert np.allclose(np.linalg.norm(mot.p, axis=1), 1.0)
    pol = DcmPolarityUpdater(motility=mot, D_r_eff=1.0e4, dt_reorient=2.0e-7, seed=7)
    p_before = mot.p.copy()
    for k in range(50):
        pol.act(k)
    assert np.allclose(np.linalg.norm(mot.p, axis=1), 1.0), "polarity not unit-norm"
    # the headings actually moved (decorrelated), not frozen
    assert not np.allclose(mot.p, p_before)
