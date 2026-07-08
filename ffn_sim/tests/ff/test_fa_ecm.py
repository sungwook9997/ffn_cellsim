r"""FF FA↔ECM-fiber coupling (ff/fa_ecm) — the S4 primitive: a clutch binds a LIVE collagen node and traction
REMODELS the matrix. Analytic ground truth (Hookean two-sided spring + Newton pair). CPU device.
"""

import numpy as np
import pytest
import warp as wp

from ffn_sim.ff.fa_ecm import attach_clutches_to_ecm, clutch_ecm_spring_kernel

D = "cpu"


def _run(cell_xyz, actin, ecm_xyz, ecm_node, k_int, rest):
    cell_pos = wp.array(np.asarray(cell_xyz, float), dtype=wp.vec3d, device=D)
    ecm_pos = wp.array(np.asarray(ecm_xyz, float), dtype=wp.vec3d, device=D)
    ac = wp.array(np.asarray(actin, np.int32), dtype=wp.int32, device=D)
    en = wp.array(np.asarray(ecm_node, np.int32), dtype=wp.int32, device=D)
    cf = wp.zeros(len(cell_xyz), dtype=wp.vec3d, device=D)
    ef = wp.zeros(len(ecm_xyz), dtype=wp.vec3d, device=D)
    wp.launch(clutch_ecm_spring_kernel, dim=len(actin),
              inputs=[cell_pos, ac, ecm_pos, en, wp.float64(k_int), wp.float64(rest), cf, ef], device=D)
    return cf.numpy(), ef.numpy()


def test_hookean_force_and_newton_pair():
    """f = k_int·(L−rest)/L·(ECM−actin): +f on the actin (traction), −f on the fiber (reaction), and the pair
    sums to ZERO (no spurious net force injected into the cell+ECM system)."""
    k_int, rest, d = 1000.0, 0.05, 0.20
    cf, ef = _run([[0.0, 0, 0]], [0], [[d, 0, 0]], [0], k_int, rest)
    fmag = k_int * (d - rest)                                  # along +x
    assert cf[0] == pytest.approx([fmag, 0, 0])                # traction pulls actin toward the fiber (+x)
    assert ef[0] == pytest.approx([-fmag, 0, 0])               # reaction pulls the fiber toward the cell (−x)
    assert np.allclose(cf.sum(0) + ef.sum(0), 0.0, atol=1e-9)  # Newton: total force on the coupled system = 0


def test_traction_remodels_toward_the_cell():
    """The reaction on the collagen node points from the fiber TOWARD the actin — i.e. a loaded clutch recruits
    the fiber to the cell (the sign that gives the emergent inward matrix remodel)."""
    cf, ef = _run([[0.0, 0, 0]], [0], [[0.3, 0.0, 0.0]], [0], 800.0, 0.05)
    # fiber at +x, actin at origin → reaction on fiber must have −x component (moves fiber toward the cell)
    assert ef[0][0] < 0.0
    assert np.dot(ef[0], np.array([0.0, 0, 0]) - np.array([0.3, 0, 0])) > 0   # points fiber→actin


def test_unbound_clutch_exerts_nothing():
    cf, ef = _run([[0.0, 0, 0]], [0], [[0.2, 0, 0]], [-1], 1000.0, 0.05)      # ecm_node<0 = no fiber bound
    assert np.allclose(cf, 0.0) and np.allclose(ef, 0.0)


def test_attach_binds_nearest_within_capture():
    actin = np.array([[0.0, 0, 0]])
    ecm = np.array([[0.1, 0, 0], [5.0, 0, 0]])                 # near fiber node 0, far node 1
    assert attach_clutches_to_ecm(actin, ecm, capture_um=0.5)[0] == 0     # binds the near node
    assert attach_clutches_to_ecm(actin, ecm, capture_um=0.05)[0] == -1   # none within a tight capture
    assert attach_clutches_to_ecm(actin, np.zeros((0, 3)), capture_um=0.5)[0] == -1   # no ECM → unbound


def test_attach_then_force_are_consistent():
    """End-to-end: attach picks a fiber, the kernel then loads that exact node (bound index feeds the spring)."""
    actin = np.array([[0.0, 0, 0]]); ecm = np.array([[0.15, 0, 0], [9.0, 0, 0]])
    en = attach_clutches_to_ecm(actin, ecm, capture_um=0.5)
    assert en[0] == 0
    cf, ef = _run(actin, [0], ecm, en, 1000.0, 0.05)
    assert ef[0][0] < 0 and abs(ef[1][0]) == 0.0              # only the bound fiber node (0) feels the reaction
