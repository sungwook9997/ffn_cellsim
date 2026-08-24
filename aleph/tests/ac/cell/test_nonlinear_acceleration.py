"""CUDA and algebra gates for source-independent nonlinear mechanics accelerators."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest


def _cuda_device(wp: object) -> str:
    device = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device is None:
        pytest.skip("I0-A: nonlinear accelerator gates require a CUDA GPU")
    return str(device)


def test_cuda_block_descent_matches_live_diagonal_oracle() -> None:
    """The direct block direction is ``P M^-1 P F`` and is not mislabeled as a converged solve."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import ProjectedAnalyticCG

    wp.init()
    device = _cuda_device(wp)
    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.1, 0.0], [1.0, 0.0, 0.1]], dtype=np.float64)
    rhs = np.array([[1.0, -0.5, 0.25], [-0.2, 0.7, -0.3], [0.4, -0.1, 0.8]], dtype=np.float64)
    triples = np.array([[0, 1, 2]], dtype=np.int32)
    alpha = 2.5
    regularization = 0.75
    cell = SimpleNamespace(
        device=device,
        n_total=3,
        n_fibers=0,
        srest_d=wp.zeros(0, dtype=wp.float64, device=device),
        n_tri=1,
        tri_d=wp.array(triples, dtype=wp.int32, device=device),
        alpha_d=wp.array(np.array([alpha]), dtype=wp.float64, device=device),
        n_xl=0,
        steric=None,
        myosin=None,
        nucleus=None,
        membrane=None,
    )
    workspace = ProjectedAnalyticCG(cell, max_iterations=1)
    finite_d = wp.ones(1, dtype=wp.int32, device=device)
    measured_d = workspace.preconditioned_direction(
        wp.array(positions, dtype=wp.vec3d, device=device),
        wp.array(rhs, dtype=wp.vec3d, device=device),
        wp.array(np.array([regularization]), dtype=wp.float64, device=device),
        finite_d,
    )
    expected = rhs / np.array([
        regularization + alpha,
        regularization + 4.0 * alpha,
        regularization + alpha,
    ])[:, None]
    np.testing.assert_allclose(measured_d.numpy(), expected, rtol=2.0e-15, atol=2.0e-15)
    assert int(finite_d.numpy()[0]) == 1


def test_cuda_anderson_depth_one_matches_dense_formula_and_singular_fallback() -> None:
    """Depth-one Anderson evaluates its FP64 secant formula and has an exact base-map fallback."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.nonlinear_acceleration import AndersonDepthOne

    wp.init()
    device = _cuda_device(wp)
    accelerator = AndersonDepthOne(2, device)
    active_d = wp.ones(1, dtype=wp.int32, device=device)
    finite_d = wp.ones(1, dtype=wp.int32, device=device)

    pos0 = np.array([[0.0, 0.0, 0.0], [1.0, -1.0, 0.5]], dtype=np.float64)
    residual0 = np.array([[0.3, -0.2, 0.1], [-0.4, 0.2, 0.5]], dtype=np.float64)
    first_d = accelerator.direction(
        wp.array(pos0, dtype=wp.vec3d, device=device),
        wp.array(residual0, dtype=wp.vec3d, device=device),
        active_d,
        finite_d,
    )
    np.testing.assert_array_equal(first_d.numpy(), residual0)

    pos1 = np.array([[0.1, 0.2, -0.1], [0.8, -0.7, 0.4]], dtype=np.float64)
    residual1 = np.array([[0.2, -0.1, 0.4], [-0.1, 0.3, 0.2]], dtype=np.float64)
    second_d = accelerator.direction(
        wp.array(pos1, dtype=wp.vec3d, device=device),
        wp.array(residual1, dtype=wp.vec3d, device=device),
        active_d,
        finite_d,
    )
    delta_f = residual1 - residual0
    mapping0 = pos0 + residual0
    mapping1 = pos1 + residual1
    gamma = float(np.vdot(delta_f, residual1) / np.vdot(delta_f, delta_f))
    expected = mapping1 - gamma * (mapping1 - mapping0) - pos1
    np.testing.assert_allclose(second_d.numpy(), expected, rtol=3.0e-15, atol=3.0e-15)
    assert float(accelerator.coefficient.numpy()[0]) == pytest.approx(gamma, rel=3.0e-15)

    # Same fixed-point residual gives Delta f=0. The coefficient is exactly zero and the base direction wins.
    pos2 = pos1 + 0.05
    third_d = accelerator.direction(
        wp.array(pos2, dtype=wp.vec3d, device=device),
        wp.array(residual1, dtype=wp.vec3d, device=device),
        active_d,
        finite_d,
    )
    np.testing.assert_allclose(third_d.numpy(), residual1, rtol=0.0, atol=4.0e-17)
    assert float(accelerator.coefficient.numpy()[0]) == 0.0


def test_cuda_rkc1_recurrence_matches_chebyshev_stage() -> None:
    """The device recurrence evaluates ``Yj=2Yj-1-Yj-2+2dtF(Yj-1)`` componentwise."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.nonlinear_acceleration import rkc1_recurrence_kernel

    wp.init()
    device = _cuda_device(wp)
    previous = np.array([[0.1, -0.2, 0.3], [0.4, 0.5, -0.6]], dtype=np.float64)
    current = np.array([[0.2, -0.1, 0.4], [0.3, 0.7, -0.2]], dtype=np.float64)
    force = np.array([[-1.0, 0.5, -0.25], [0.2, -0.4, 0.8]], dtype=np.float64)
    dt_mu = 0.125
    measured_d = wp.zeros(2, dtype=wp.vec3d, device=device)
    wp.launch(
        rkc1_recurrence_kernel,
        dim=2,
        inputs=[
            wp.array(previous, dtype=wp.vec3d, device=device),
            wp.array(current, dtype=wp.vec3d, device=device),
            wp.array(force, dtype=wp.vec3d, device=device),
            wp.array(np.array([dt_mu]), dtype=wp.float64, device=device),
            wp.ones(1, dtype=wp.int32, device=device),
            wp.ones(1, dtype=wp.int32, device=device),
            measured_d,
        ],
        device=device,
    )
    expected = 2.0 * current - previous + 2.0 * dt_mu * force
    np.testing.assert_allclose(measured_d.numpy(), expected, rtol=0.0, atol=2.0e-16)


def test_cuda_backtracking_cannot_fake_mechanical_convergence() -> None:
    """A tiny accepted displacement cannot converge while its explicit-equivalent force remains large."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.inner_mechanics import update_convergence_kernel

    wp.init()
    device = _cuda_device(wp)
    max_displacement_d = wp.array(np.array([0.5]), dtype=wp.float64, device=device)
    max_constraint_d = wp.array(np.array([0.1]), dtype=wp.float64, device=device)
    max_force_d = wp.array(np.array([2.0]), dtype=wp.float64, device=device)
    dt_mu_d = wp.array(np.array([1.0]), dtype=wp.float64, device=device)
    finite_d = wp.ones(1, dtype=wp.int32, device=device)
    active_d = wp.ones(1, dtype=wp.int32, device=device)
    converged_d = wp.zeros(1, dtype=wp.int32, device=device)
    iterations_d = wp.zeros(1, dtype=wp.int32, device=device)

    wp.launch(
        update_convergence_kernel,
        dim=1,
        inputs=[
            max_displacement_d, max_constraint_d, max_force_d, dt_mu_d,
            wp.float64(1.0), wp.int32(7), finite_d, active_d, converged_d, iterations_d,
        ],
        device=device,
    )
    assert int(active_d.numpy()[0]) == 1
    assert int(converged_d.numpy()[0]) == 0

    max_force_d.fill_(0.5)
    wp.launch(
        update_convergence_kernel,
        dim=1,
        inputs=[
            max_displacement_d, max_constraint_d, max_force_d, dt_mu_d,
            wp.float64(1.0), wp.int32(8), finite_d, active_d, converged_d, iterations_d,
        ],
        device=device,
    )
    assert int(active_d.numpy()[0]) == 0
    assert int(converged_d.numpy()[0]) == 1
    assert int(iterations_d.numpy()[0]) == 8
