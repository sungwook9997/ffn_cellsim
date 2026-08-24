"""Algebra gates for the overlapping vector-block contact Schwarz solve."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from aleph.components.incumbent.contact_schwarz import additive_pair_schwarz_oracle


def test_isolated_contact_pair_is_the_exact_six_dof_block_solve() -> None:
    unit = np.array([1.0, 2.0, -1.0])
    unit /= np.linalg.norm(unit)
    coupling = 7.0 * np.outer(unit, unit) + 2.0 * np.eye(3)
    diagonal = np.stack([coupling + 3.0 * np.eye(3), coupling + 5.0 * np.eye(3)])
    residual = np.array([[2.0, -1.0, 4.0], [-3.0, 5.0, 1.0]])
    edges = np.array([[0, 1]], dtype=np.int32)

    measured = additive_pair_schwarz_oracle(residual, diagonal, edges, coupling[None, :, :])
    system = np.block([[diagonal[0], -coupling], [-coupling, diagonal[1]]])
    expected = np.linalg.solve(system, residual.reshape(-1)).reshape(2, 3)
    np.testing.assert_allclose(measured, expected, rtol=3.0e-15, atol=3.0e-15)


def test_overlapping_pair_operator_is_symmetric_positive_definite() -> None:
    edges = np.array([[0, 1], [1, 2]], dtype=np.int32)
    directions = np.array([[1.0, 0.3, -0.2], [-0.4, 1.0, 0.5]])
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    blocks = np.stack([4.0 * np.outer(directions[0], directions[0]),
                       6.0 * np.outer(directions[1], directions[1])])
    diagonal = np.repeat((2.0 * np.eye(3))[None, :, :], 3, axis=0)
    diagonal[0] += blocks[0]
    diagonal[1] += blocks[0] + blocks[1]
    diagonal[2] += blocks[1]

    columns = []
    for column in range(9):
        basis = np.zeros((3, 3))
        basis.reshape(-1)[column] = 1.0
        columns.append(additive_pair_schwarz_oracle(basis, diagonal, edges, blocks).reshape(-1))
    operator = np.column_stack(columns)
    np.testing.assert_allclose(operator, operator.T, rtol=0.0, atol=1.0e-16)
    assert np.linalg.eigvalsh(operator).min() > 0.0


def test_cuda_isolated_wca_crosslink_pair_matches_dense_block() -> None:
    """CUDA assembles coincident WCA+xlink curvature and solves the exact degree-one pair block."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.contact_schwarz import ContactPairSchwarz
    from aleph.components.solid.steric_warp import StericForce
    from aleph.components.solid.wca_analytic import epsilon_from_contact_stiffness

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device_obj is None:
        pytest.skip("I0-A: contact Schwarz runtime gate requires a CUDA GPU")
    device = str(device_obj)
    sigma = 1.0
    distance = 0.95
    k_ev = 10.0
    k_xlink = 3.0
    regularization = 2.0
    positions = np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]], dtype=np.float64)
    residual = np.array([[2.0, -1.0, 0.5], [-3.0, 4.0, 1.5]], dtype=np.float64)
    steric = StericForce(
        np.array([0, 1], dtype=np.int32), sigma, k_ev, device,
        force_cap=None, grid_dim=(8, 8, 8),
    )
    cell = SimpleNamespace(
        device=device,
        n_total=2,
        n_fibers=0,
        srest_d=wp.zeros(0, dtype=wp.float64, device=device),
        steric=steric,
        xl_d=wp.array(np.array([[0, 1]], np.int32), dtype=wp.int32, device=device),
        kxl_d=wp.array(np.array([k_xlink]), dtype=wp.float64, device=device),
        r0xl_d=wp.array(np.array([distance]), dtype=wp.float64, device=device),
    )
    workspace = ContactPairSchwarz(cell)
    measured_d = workspace.direction(
        wp.array(positions, dtype=wp.vec3d, device=device),
        wp.array(residual, dtype=wp.vec3d, device=device),
        wp.array(np.array([regularization]), dtype=wp.float64, device=device),
        wp.ones(1, dtype=wp.int32, device=device),
        wp.ones(1, dtype=wp.int32, device=device),
    )
    wp.synchronize_device(device)

    epsilon = epsilon_from_contact_stiffness(k_ev, sigma)
    sr6 = (sigma / distance) ** 6
    radial = 24.0 * epsilon / distance**2 * (26.0 * sr6**2 - 7.0 * sr6)
    coupling = np.diag([radial + k_xlink, 0.0, 0.0])
    diagonal = regularization * np.eye(3) + coupling
    system = np.block([[diagonal, -coupling], [-coupling, diagonal]])
    expected = np.linalg.solve(system, residual.reshape(-1)).reshape(2, 3)
    np.testing.assert_allclose(measured_d.numpy(), expected, rtol=3.0e-12, atol=3.0e-12)


def test_cuda_two_single_node_fibers_match_dense_translation_block() -> None:
    """The fiber-cluster CUDA path reduces to the same exact block for one-node fibers."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.contact_schwarz import FiberContactSchwarz
    from aleph.components.solid.steric_warp import StericForce
    from aleph.components.solid.wca_analytic import epsilon_from_contact_stiffness

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device_obj is None:
        pytest.skip("I0-A: fiber contact Schwarz runtime gate requires a CUDA GPU")
    device = str(device_obj)
    sigma, distance, k_ev, k_xlink, regularization = 1.0, 0.95, 10.0, 3.0, 2.0
    positions = np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]], dtype=np.float64)
    residual = np.array([[2.0, -1.0, 0.5], [-3.0, 4.0, 1.5]], dtype=np.float64)
    steric = StericForce(
        np.array([0, 1], dtype=np.int32), sigma, k_ev, device,
        force_cap=None, grid_dim=(8, 8, 8),
    )
    cell = SimpleNamespace(
        device=device, n_total=2, n_fibers=2, steric=steric,
        foff_d=wp.array(np.array([0, 1, 2], np.int32), dtype=wp.int32, device=device),
        soff_d=wp.array(np.array([0, 0, 0], np.int32), dtype=wp.int32, device=device),
        srest_d=wp.zeros(0, dtype=wp.float64, device=device),
        xl_d=wp.array(np.array([[0, 1]], np.int32), dtype=wp.int32, device=device),
        kxl_d=wp.array(np.array([k_xlink]), dtype=wp.float64, device=device),
        r0xl_d=wp.array(np.array([distance]), dtype=wp.float64, device=device),
    )
    workspace = FiberContactSchwarz(cell)
    measured_d = workspace.direction(
        wp.array(positions, dtype=wp.vec3d, device=device),
        wp.array(residual, dtype=wp.vec3d, device=device),
        wp.array(np.array([regularization]), dtype=wp.float64, device=device),
        wp.ones(1, dtype=wp.int32, device=device),
        wp.ones(1, dtype=wp.int32, device=device),
    )
    wp.synchronize_device(device)

    epsilon = epsilon_from_contact_stiffness(k_ev, sigma)
    sr6 = (sigma / distance) ** 6
    radial = 24.0 * epsilon / distance**2 * (26.0 * sr6**2 - 7.0 * sr6)
    coupling = np.diag([radial + k_xlink, 0.0, 0.0])
    diagonal = regularization * np.eye(3) + coupling
    expected = np.linalg.solve(
        np.block([[diagonal, -coupling], [-coupling, diagonal]]), residual.reshape(-1),
    ).reshape(2, 3)
    np.testing.assert_allclose(measured_d.numpy(), expected, rtol=3.0e-12, atol=3.0e-12)


def test_cuda_rigid_fiber_contact_direction_is_finite_and_constraint_tangent() -> None:
    """The 12-DoF cluster path produces a finite nonzero NF2007-tangent direction."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.contact_schwarz import RigidFiberContactSchwarz
    from aleph.components.solid.steric_warp import StericForce

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device_obj is None:
        pytest.skip("I0-A: rigid fiber contact Schwarz runtime gate requires a CUDA GPU")
    device = str(device_obj)
    positions = np.array([
        [0.0, 0.0, 0.0], [0.0, 0.0, 1.0],
        [0.95, 0.0, 0.0], [2.0, 0.0, 1.0],
    ], dtype=np.float64)
    residual = np.array([
        [2.0, -1.0, 0.5], [1.0, 0.5, -0.25],
        [-3.0, 4.0, 1.5], [-0.5, -1.0, 0.25],
    ], dtype=np.float64)
    segments = np.array([np.linalg.norm(positions[1] - positions[0]),
                         np.linalg.norm(positions[3] - positions[2])])
    steric = StericForce(
        np.array([0, 0, 1, 1], dtype=np.int32), 1.0, 10.0, device,
        force_cap=None, grid_dim=(8, 8, 8),
    )
    cell = SimpleNamespace(
        device=device, n_total=4, n_fibers=2, n_xl=0, steric=steric,
        foff_d=wp.array(np.array([0, 2, 4], np.int32), dtype=wp.int32, device=device),
        soff_d=wp.array(np.array([0, 1, 2], np.int32), dtype=wp.int32, device=device),
        srest_d=wp.array(segments, dtype=wp.float64, device=device),
        xl_d=wp.array(np.zeros((0, 2), np.int32), dtype=wp.int32, device=device),
        kxl_d=wp.zeros(0, dtype=wp.float64, device=device),
        r0xl_d=wp.zeros(0, dtype=wp.float64, device=device),
    )
    workspace = RigidFiberContactSchwarz(cell)
    finite_d = wp.ones(1, dtype=wp.int32, device=device)
    measured = workspace.direction(
        wp.array(positions, dtype=wp.vec3d, device=device),
        wp.array(residual, dtype=wp.vec3d, device=device),
        wp.array(np.array([2.0]), dtype=wp.float64, device=device),
        wp.ones(1, dtype=wp.int32, device=device), finite_d,
    ).numpy()
    assert int(finite_d.numpy()[0]) == 1
    assert np.isfinite(measured).all()
    assert np.linalg.norm(measured) > 0.0
    for i, j in ((0, 1), (2, 3)):
        assert abs(float(np.dot(positions[j] - positions[i], measured[j] - measured[i]))) < 2.0e-12
