"""Acceptance algebra and CUDA parity gates for analytic projected implicit mechanics."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from aleph.components.incumbent.implicit_mechanics_analytic import (
    bending_stiffness,
    central_spring_tangent,
    projected_implicit_solve,
)


def test_bending_stiffness_is_exact_symmetric_psd_stencil() -> None:
    triples = np.array([[0, 1, 2]], dtype=np.int32)
    stiffness = bending_stiffness(3, triples, np.array([2.5]))
    expected_scalar = 2.5 * np.outer([1.0, -2.0, 1.0], [1.0, -2.0, 1.0])
    np.testing.assert_array_equal(stiffness, np.kron(expected_scalar, np.eye(3)))
    np.testing.assert_array_equal(stiffness, stiffness.T)
    roundoff_floor = -64.0 * np.finfo(np.float64).eps * np.linalg.norm(stiffness, ord=2)
    assert np.linalg.eigvalsh(stiffness).min() >= roundoff_floor
    np.testing.assert_allclose(stiffness.reshape(3, 3, 3, 3).sum(axis=2), 0.0, atol=0.0)


def test_central_spring_tangent_has_axial_and_only_nonnegative_transverse_curvature() -> None:
    tensile = central_spring_tangent([2.0, 0.0, 0.0], stiffness=6.0, rest_length=1.0)
    np.testing.assert_allclose(tensile, np.diag([6.0, 3.0, 3.0]), rtol=0.0, atol=0.0)

    compressed = central_spring_tangent([0.5, 0.0, 0.0], stiffness=6.0, rest_length=1.0)
    np.testing.assert_allclose(compressed, np.diag([6.0, 0.0, 0.0]), rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(
        central_spring_tangent([2.0, 0.0, 0.0], stiffness=6.0, rest_length=2.0, tension_only=True),
        np.zeros((3, 3)),
    )


def test_dense_projected_solve_satisfies_tangent_spd_system() -> None:
    stiffness = np.diag([9.0, 4.0, 1.0])
    normal = np.array([1.0, 2.0, -1.0])
    normal /= np.linalg.norm(normal)
    projector = np.eye(3) - np.outer(normal, normal)
    force = np.array([2.0, -3.0, 5.0])
    displacement = projected_implicit_solve(force, stiffness, projector, regularization=0.75)
    np.testing.assert_allclose(
        (0.75 * np.eye(3) + projector @ stiffness @ projector) @ displacement,
        projector @ force,
        rtol=2.0e-14,
        atol=2.0e-14,
    )
    assert abs(float(normal @ displacement)) <= 2.0e-15


def test_omitted_regularization_uses_declared_stiffness_and_machine_floor() -> None:
    from aleph.components.incumbent.implicit_mechanics import omitted_regularization_base

    nucleus = SimpleNamespace(
        mean_edge_um=2.0, kappa_tilde=8.0, k_soft=12.0, k_ac=25.0, k_vol=1000.0)
    membrane = SimpleNamespace(
        mean_edge_um=0.5, kappa_tilde=0.125, gamma_mem=10.0, with_area_tension=True)
    cell = SimpleNamespace(mechanical_kmax=1.0e6, nucleus=nucleus, membrane=membrane)
    assert omitted_regularization_base(cell) == 1000.0

    floor_only = SimpleNamespace(mechanical_kmax=1.0e6, nucleus=None, membrane=None)
    assert omitted_regularization_base(floor_only) == pytest.approx(
        np.sqrt(np.finfo(np.float64).eps) * 1.0e6)


def test_cuda_projected_cg_matches_dense_bending_oracle() -> None:
    """The fixed-launch all-device CG returns the same 3-node bending solve as dense algebra."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import ProjectedAnalyticCG

    wp.init()
    cuda_device = next((device for device in wp.get_devices() if device.is_cuda), None)
    if cuda_device is None:
        pytest.skip("I0-A: analytic implicit runtime gate requires a CUDA GPU")
    device = str(cuda_device)

    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.1, 0.0], [1.0, 0.0, 0.1]], dtype=np.float64)
    rhs = np.array([[1.0, -0.5, 0.25], [-0.2, 0.7, -0.3], [0.4, -0.1, 0.8]], dtype=np.float64)
    triples = np.array([[0, 1, 2]], dtype=np.int32)
    alpha = np.array([2.5], dtype=np.float64)
    regularization = 0.75
    cell = SimpleNamespace(
        device=device,
        n_total=3,
        n_fibers=0,
        srest_d=wp.zeros(0, dtype=wp.float64, device=device),
        n_tri=1,
        tri_d=wp.array(triples, dtype=wp.int32, device=device),
        alpha_d=wp.array(alpha, dtype=wp.float64, device=device),
        n_xl=0,
        steric=None,
        myosin=None,
        nucleus=None,
        membrane=None,
    )
    solver = ProjectedAnalyticCG(cell, max_iterations=8)
    finite_d = wp.ones(1, dtype=wp.int32, device=device)
    measured_d = solver.solve(
        wp.array(positions, dtype=wp.vec3d, device=device),
        wp.array(rhs, dtype=wp.vec3d, device=device),
        wp.array(np.array([regularization]), dtype=wp.float64, device=device),
        finite_d,
    )
    wp.synchronize_device(device)

    expected = np.linalg.solve(
        regularization * np.eye(9) + bending_stiffness(3, triples, alpha),
        rhs.reshape(-1),
    ).reshape(3, 3)
    np.testing.assert_allclose(measured_d.numpy(), expected, rtol=2.0e-12, atol=2.0e-12)
    assert int(finite_d.numpy()[0]) == 1
    assert int(solver.converged.numpy()[0]) == 1
    assert 1 <= int(solver.iterations.numpy()[0]) <= 8


def test_cuda_wca_radial_tangent_matches_closed_form() -> None:
    """The hash-grid tangent keeps exact positive radial WCA curvature and Newton symmetry."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import add_wca_stiffness_kernel
    from aleph.components.solid.steric_warp import StericForce
    from aleph.components.solid.wca_analytic import epsilon_from_contact_stiffness

    wp.init()
    cuda_device = next((device for device in wp.get_devices() if device.is_cuda), None)
    if cuda_device is None:
        pytest.skip("I0-A: WCA implicit-tangent gate requires a CUDA GPU")
    device = str(cuda_device)
    sigma = 1.0
    distance = 0.95
    epsilon = epsilon_from_contact_stiffness(10.0, sigma)
    positions = wp.array(
        np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]], dtype=np.float64),
        dtype=wp.vec3d,
        device=device,
    )
    vector = wp.array(np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]), dtype=wp.vec3d, device=device)
    steric = StericForce(
        np.array([0, 1], dtype=np.int32), sigma, 10.0, device,
        force_cap=None, grid_dim=(8, 8, 8),
    )
    steric._accumulate_pos(positions, wp.zeros(2, dtype=wp.vec3d, device=device))
    measured_d = wp.zeros(2, dtype=wp.vec3d, device=device)
    wp.launch(
        add_wca_stiffness_kernel,
        dim=2,
        inputs=[
            steric.grid.id, steric._qpts, positions, vector, steric.fiber_id, steric.active,
            wp.float32(steric.r_c), wp.float64(sigma), wp.float64(epsilon), wp.float64(-1.0), measured_d,
        ],
        device=device,
    )
    wp.synchronize_device(device)

    sr6 = (sigma / distance) ** 6
    radial_k = 24.0 * epsilon / distance**2 * (26.0 * sr6**2 - 7.0 * sr6)
    np.testing.assert_allclose(
        measured_d.numpy(), [[radial_k, 0.0, 0.0], [-radial_k, 0.0, 0.0]],
        rtol=2.0e-12, atol=2.0e-12,
    )


def test_cuda_cg_matches_dense_projected_system() -> None:
    """The matrix-free ``P K P`` path matches an independently constructed dense NF2007 projector."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import ProjectedAnalyticCG

    wp.init()
    cuda_device = next((device for device in wp.get_devices() if device.is_cuda), None)
    if cuda_device is None:
        pytest.skip("I0-A: projected implicit runtime gate requires a CUDA GPU")
    device = str(cuda_device)

    positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.1, 0.0], [0.9, 0.35, 0.1]], dtype=np.float64)
    triples = np.array([[0, 1, 2]], dtype=np.int32)
    alpha = np.array([3.25], dtype=np.float64)
    raw_force = np.array([[1.0, -0.5, 0.25], [-0.2, 0.7, -0.3], [0.4, -0.1, 0.8]], dtype=np.float64)
    jacobian = np.zeros((2, 9), dtype=np.float64)
    for segment in range(2):
        g = positions[segment + 1] - positions[segment]
        jacobian[segment, 3 * segment:3 * segment + 3] = -2.0 * g
        jacobian[segment, 3 * (segment + 1):3 * (segment + 1) + 3] = 2.0 * g
    projector = np.eye(9) - jacobian.T @ np.linalg.solve(jacobian @ jacobian.T, jacobian)
    rhs = (projector @ raw_force.reshape(-1)).reshape(3, 3)
    regularization = 1.5

    cell = SimpleNamespace(
        device=device,
        n_total=3,
        n_fibers=1,
        foff_d=wp.array(np.array([0, 3], dtype=np.int32), dtype=wp.int32, device=device),
        soff_d=wp.array(np.array([0, 2], dtype=np.int32), dtype=wp.int32, device=device),
        srest_d=wp.array(np.linalg.norm(np.diff(positions, axis=0), axis=1), dtype=wp.float64, device=device),
        n_tri=1,
        tri_d=wp.array(triples, dtype=wp.int32, device=device),
        alpha_d=wp.array(alpha, dtype=wp.float64, device=device),
        n_xl=0,
        steric=None,
        myosin=None,
        nucleus=None,
        membrane=None,
    )
    solver = ProjectedAnalyticCG(cell, max_iterations=32)
    finite_d = wp.ones(1, dtype=wp.int32, device=device)
    measured_d = solver.solve(
        wp.array(positions, dtype=wp.vec3d, device=device),
        wp.array(rhs, dtype=wp.vec3d, device=device),
        wp.array(np.array([regularization]), dtype=wp.float64, device=device),
        finite_d,
    )
    wp.synchronize_device(device)

    expected = projected_implicit_solve(
        raw_force, bending_stiffness(3, triples, alpha), projector, regularization=regularization,
    ).reshape(3, 3)
    np.testing.assert_allclose(measured_d.numpy(), expected, rtol=2.0e-11, atol=2.0e-11)
    np.testing.assert_allclose(jacobian @ measured_d.numpy().reshape(-1), 0.0, rtol=0.0, atol=2.0e-12)
    assert int(finite_d.numpy()[0]) == 1
    assert int(solver.converged.numpy()[0]) == 1


def test_cuda_line_search_rejects_unconverged_linear_candidate() -> None:
    """A smaller line-search residual cannot authorize displacement from unconverged PCG."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import (
        commit_better_line_search_trial_kernel,
        conditional_copy_vec3_kernel,
        decide_better_line_search_trial_kernel,
        finalize_line_search_kernel,
    )

    wp.init()
    cuda_device = next((device for device in wp.get_devices() if device.is_cuda), None)
    if cuda_device is None:
        pytest.skip("I0-A: implicit residual-guard gate requires a CUDA GPU")
    device = str(cuda_device)
    trial_d = wp.array(np.array([[1.0, 0.0, 0.0]]), dtype=wp.vec3d, device=device)
    best_d = wp.array(np.array([[2.0, 0.0, 0.0]]), dtype=wp.vec3d, device=device)
    trial_residual_d = wp.array(np.array([1.0]), dtype=wp.float64, device=device)
    best_residual_d = wp.array(np.array([2.0]), dtype=wp.float64, device=device)
    trial_finite_d = wp.ones(1, dtype=wp.int32, device=device)
    best_finite_d = wp.ones(1, dtype=wp.int32, device=device)
    converged_d = wp.zeros(1, dtype=wp.int32, device=device)
    active_d = wp.ones(1, dtype=wp.int32, device=device)
    take_d = wp.zeros(1, dtype=wp.int32, device=device)
    best_index_d = wp.array(np.array([-2], np.int32), dtype=wp.int32, device=device)
    output_finite_d = wp.zeros(1, dtype=wp.int32, device=device)
    accepted_d = wp.zeros(1, dtype=wp.int32, device=device)
    explicit_d = wp.zeros(1, dtype=wp.int32, device=device)
    stationary_d = wp.zeros(1, dtype=wp.int32, device=device)
    scale_counts_d = wp.zeros(3, dtype=wp.int32, device=device)

    decision_inputs = [
        trial_residual_d, best_residual_d, trial_finite_d, best_finite_d,
        converged_d, active_d, take_d,
    ]
    wp.launch(decide_better_line_search_trial_kernel, dim=1, inputs=decision_inputs, device=device)
    wp.launch(conditional_copy_vec3_kernel, dim=1, inputs=[trial_d, best_d, take_d], device=device)
    wp.launch(
        commit_better_line_search_trial_kernel,
        dim=1,
        inputs=[
            trial_residual_d, trial_finite_d, wp.int32(2), take_d,
            best_residual_d, best_finite_d, best_index_d,
        ],
        device=device,
    )
    wp.launch(
        finalize_line_search_kernel,
        dim=1,
        inputs=[best_finite_d, best_index_d, output_finite_d, accepted_d,
                explicit_d, stationary_d, scale_counts_d],
        device=device,
    )
    np.testing.assert_array_equal(best_d.numpy(), [[2.0, 0.0, 0.0]])
    assert int(accepted_d.numpy()[0]) == 0
    assert int(stationary_d.numpy()[0]) == 1
    np.testing.assert_array_equal(scale_counts_d.numpy(), [0, 0, 0])

    converged_d.fill_(1)
    wp.launch(decide_better_line_search_trial_kernel, dim=1, inputs=decision_inputs, device=device)
    wp.launch(conditional_copy_vec3_kernel, dim=1, inputs=[trial_d, best_d, take_d], device=device)
    wp.launch(
        commit_better_line_search_trial_kernel,
        dim=1,
        inputs=[
            trial_residual_d, trial_finite_d, wp.int32(2), take_d,
            best_residual_d, best_finite_d, best_index_d,
        ],
        device=device,
    )
    wp.launch(
        finalize_line_search_kernel,
        dim=1,
        inputs=[best_finite_d, best_index_d, output_finite_d, accepted_d,
                explicit_d, stationary_d, scale_counts_d],
        device=device,
    )
    np.testing.assert_array_equal(best_d.numpy(), [[1.0, 0.0, 0.0]])
    assert int(accepted_d.numpy()[0]) == 1
    assert int(stationary_d.numpy()[0]) == 1
    np.testing.assert_array_equal(scale_counts_d.numpy(), [0, 0, 1])
    assert int(output_finite_d.numpy()[0]) == 1


def test_segment_crossbridge_tangent_matches_dense_barycentric_oracle_on_cuda() -> None:
    """The implicit three-node segment tangent matches the explicit barycentric force map."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import add_segment_crossbridge_stiffness_kernel

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device_obj is None:
        pytest.skip("I0-A: segment implicit-tangent gate requires a CUDA GPU")
    device = str(device_obj)
    pos = np.array([[0.0, 0.0, 0.0], [1.0, 0.2, 0.0], [2.0, 0.2, 0.0]], dtype=np.float64)
    vector = np.array([[0.3, -0.2, 0.1], [-0.1, 0.4, 0.2], [0.2, -0.3, 0.5]], dtype=np.float64)
    t = 0.25
    abscissa = 0.1
    walk = np.array([[1.0, 0.0, 0.0]], dtype=np.float64)
    stiffness = 3.0
    rest = 0.1
    out = wp.zeros(3, dtype=wp.vec3d, device=device)
    wp.launch(
        add_segment_crossbridge_stiffness_kernel,
        dim=1,
        inputs=[
            wp.array(pos, dtype=wp.vec3d, device=device),
            wp.array(vector, dtype=wp.vec3d, device=device),
            wp.array(np.array([0], np.int32), dtype=wp.int32, device=device),
            wp.ones(1, dtype=wp.int32, device=device),
            wp.array(np.array([1], np.int32), dtype=wp.int32, device=device),
            wp.array(np.array([2], np.int32), dtype=wp.int32, device=device),
            wp.array(np.array([t]), dtype=wp.float64, device=device),
            wp.array(np.array([abscissa]), dtype=wp.float64, device=device),
            wp.array(walk, dtype=wp.vec3d, device=device),
            wp.float64(stiffness),
            wp.float64(rest),
            out,
        ],
        device=device,
    )
    attach = (1.0 - t) * pos[1] + t * pos[2] + abscissa * walk[0]
    delta = attach - pos[0]
    length = np.linalg.norm(delta)
    unit = delta / length
    transverse = max(stiffness * (length - rest) / length, 0.0)
    tangent = stiffness * np.outer(unit, unit) + transverse * (np.eye(3) - np.outer(unit, unit))
    relative = vector[0] - (1.0 - t) * vector[1] - t * vector[2]
    action = tangent @ relative
    expected = np.stack([action, -(1.0 - t) * action, -t * action])
    np.testing.assert_allclose(out.numpy(), expected, rtol=2.0e-13, atol=2.0e-13)


@pytest.mark.parametrize("straight", [False, True])
def test_cuda_angle_tangent_matches_dense_psd_oracle(straight: bool) -> None:
    """NMII angle tangents match Gauss--Newton and the two-mode straight-backbone limit."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import add_angle_gauss_newton_stiffness_kernel

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device_obj is None:
        pytest.skip("I0-A: NMII angle implicit-tangent gate requires a CUDA GPU")
    device = str(device_obj)
    if straight:
        pos = np.array([[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float64)
        theta0 = np.pi
    else:
        pos = np.array([[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 2.0, 0.0]], dtype=np.float64)
        theta0 = 0.5 * np.pi
    vector = np.array([[0.3, -0.2, 0.7], [-0.1, 0.4, -0.2], [0.2, -0.3, 0.5]], dtype=np.float64)
    stiffness = 3.25
    out_d = wp.zeros(3, dtype=wp.vec3d, device=device)
    wp.launch(
        add_angle_gauss_newton_stiffness_kernel,
        dim=1,
        inputs=[
            wp.array(pos, dtype=wp.vec3d, device=device),
            wp.array(vector, dtype=wp.vec3d, device=device),
            wp.array(np.array([[0, 1, 2]], np.int32), dtype=wp.int32, device=device),
            wp.float64(stiffness),
            wp.float64(theta0),
            out_d,
        ],
        device=device,
    )

    r1 = pos[0] - pos[1]
    r2 = pos[2] - pos[1]
    n1 = np.linalg.norm(r1)
    n2 = np.linalg.norm(r2)
    u1 = r1 / n1
    u2 = r2 / n2
    if straight:
        projection = np.eye(3) - np.outer(u2, u2)
        bend = projection @ ((vector[0] - vector[1]) / n1 + (vector[2] - vector[1]) / n2)
        action_i = stiffness * projection @ bend / n1
        action_k = stiffness * projection @ bend / n2
        expected = np.stack([action_i, -(action_i + action_k), action_k])
    else:
        cosine = float(u1 @ u2)
        sine = np.sqrt(1.0 - cosine * cosine)
        grad_i = -(u2 - cosine * u1) / (n1 * sine)
        grad_k = -(u1 - cosine * u2) / (n2 * sine)
        grad_j = -(grad_i + grad_k)
        gradient = np.concatenate([grad_i, grad_j, grad_k])
        expected = (stiffness * np.outer(gradient, gradient) @ vector.reshape(-1)).reshape(3, 3)
    np.testing.assert_allclose(out_d.numpy(), expected, rtol=3.0e-13, atol=3.0e-13)
    np.testing.assert_allclose(expected.sum(axis=0), 0.0, rtol=0.0, atol=3.0e-15)
    assert float(np.vdot(vector, expected)) >= -64.0 * np.finfo(np.float64).eps


def test_cuda_fiber_translation_coarse_block_matches_dense_oracle() -> None:
    """The topology-aware additive block solves each disjoint fiber translation mode componentwise."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import add_fiber_translation_coarse_kernel

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device_obj is None:
        pytest.skip("I0-A: fiber coarse-preconditioner gate requires a CUDA GPU")
    device = str(device_obj)
    residual = np.array([
        [1.0, 2.0, -1.0], [3.0, -1.0, 5.0],
        [-2.0, 4.0, 1.0], [6.0, 2.0, 3.0],
    ], dtype=np.float64)
    external = np.array([
        [2.0, 4.0, 5.0], [6.0, 2.0, 5.0],
        [1.0, 3.0, 2.0], [3.0, 1.0, 6.0],
    ], dtype=np.float64)
    initial = np.full((4, 3), 0.25, dtype=np.float64)
    measured_d = wp.array(initial, dtype=wp.vec3d, device=device)
    wp.launch(
        add_fiber_translation_coarse_kernel,
        dim=2,
        inputs=[
            wp.array(residual, dtype=wp.vec3d, device=device),
            wp.array(np.array([0, 2, 4], np.int32), dtype=wp.int32, device=device),
            wp.array(external, dtype=wp.vec3d, device=device),
            measured_d,
        ],
        device=device,
    )
    expected = initial.copy()
    for begin, end in ((0, 2), (2, 4)):
        expected[begin:end] += residual[begin:end].sum(axis=0) / external[begin:end].sum(axis=0)
    np.testing.assert_allclose(measured_d.numpy(), expected, rtol=0.0, atol=2.0e-15)


def test_cuda_fiber_trace_majorizer_is_isotropic_component_trace() -> None:
    """The coarse Jacobi denominator retains the full Cartesian trace on each fiber."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import restrict_fiber_trace_majorizer_kernel

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device_obj is None:
        pytest.skip("I0-A: fiber trace-majorizer gate requires a CUDA GPU")
    device = str(device_obj)
    fine = np.array([
        [1.0, 2.0, 3.0], [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0], [10.0, 11.0, 12.0],
    ], dtype=np.float64)
    measured_d = wp.zeros(2, dtype=wp.vec3d, device=device)
    wp.launch(
        restrict_fiber_trace_majorizer_kernel,
        dim=2,
        inputs=[
            wp.array(fine, dtype=wp.vec3d, device=device),
            wp.array(np.array([0, 2, 4], np.int32), dtype=wp.int32, device=device),
            measured_d,
        ],
        device=device,
    )
    traces = np.array([fine[:2].sum(), fine[2:].sum()])
    expected = np.repeat(traces[:, None], 3, axis=1)
    np.testing.assert_allclose(measured_d.numpy(), expected, rtol=0.0, atol=0.0)


def test_cuda_fiber_cholesky_block_matches_dense_bending_oracle() -> None:
    """The production topology block exactly solves bending plus the live external Cartesian diagonal."""
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import (
        apply_fiber_block_preconditioner_kernel,
        build_fiber_block_cholesky_kernel,
    )

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device_obj is None:
        pytest.skip("I0-A: fiber Cholesky-preconditioner gate requires a CUDA GPU")
    device = str(device_obj)
    external = np.array([[2.0, 4.0, 3.0], [5.0, 2.0, 6.0], [3.0, 7.0, 2.0]], dtype=np.float64)
    residual = np.array([[1.0, 2.0, -1.0], [3.0, -1.0, 5.0], [-2.0, 4.0, 1.0]], dtype=np.float64)
    alpha = 2.5
    factor_d = wp.zeros(3 * 3 * 3, dtype=wp.float64, device=device)
    finite_d = wp.ones(1, dtype=wp.int32, device=device)
    offsets_d = wp.array(np.array([0, 3], np.int32), dtype=wp.int32, device=device)
    wp.launch(
        build_fiber_block_cholesky_kernel,
        dim=1,
        inputs=[
            wp.array(external, dtype=wp.vec3d, device=device),
            offsets_d,
            wp.array(np.array([0, 1], np.int32), dtype=wp.int32, device=device),
            wp.array(np.array([alpha]), dtype=wp.float64, device=device),
            wp.int32(3),
            factor_d,
            finite_d,
        ],
        device=device,
    )
    measured_d = wp.zeros(3, dtype=wp.vec3d, device=device)
    wp.launch(
        apply_fiber_block_preconditioner_kernel,
        dim=1,
        inputs=[
            wp.array(residual, dtype=wp.vec3d, device=device),
            offsets_d,
            wp.int32(3),
            factor_d,
            wp.zeros(9, dtype=wp.float64, device=device),
            finite_d,
            measured_d,
        ],
        device=device,
    )
    stencil = np.array([1.0, -2.0, 1.0])
    expected = np.column_stack([
        np.linalg.solve(np.diag(external[:, component]) + alpha * np.outer(stencil, stencil),
                        residual[:, component])
        for component in range(3)
    ])
    np.testing.assert_allclose(measured_d.numpy(), expected, rtol=3.0e-13, atol=3.0e-13)
    assert int(finite_d.numpy()[0]) == 1


def test_rigid_strain_coarse_basis_is_orthonormal_and_rigid_null() -> None:
    """The l<=2 coarse space is L2-orthonormal; its rigid modes are exact zero-energy of a spring network."""
    from aleph.components.incumbent.implicit_mechanics_analytic import rigid_strain_coarse_basis

    rng = np.random.default_rng(11)
    n = 40
    pos = rng.normal(size=(n, 3))
    pos /= np.linalg.norm(pos, axis=1, keepdims=True)
    pos *= 7.5
    basis = rigid_strain_coarse_basis(pos, n_modes=12)
    assert basis.shape == (12, n, 3)
    flat = basis.reshape(12, -1).T
    np.testing.assert_allclose(flat.T @ flat, np.eye(12), rtol=0.0, atol=1.0e-12)

    stiffness = np.zeros((3 * n, 3 * n))
    for _ in range(300):
        i, j = rng.integers(0, n, 2)
        if i == j:
            continue
        unit = pos[j] - pos[i]
        unit /= np.linalg.norm(unit)
        block = 5.0 * np.outer(unit, unit)
        for a, sa in ((i, 1.0), (j, -1.0)):
            for b, sb in ((i, 1.0), (j, -1.0)):
                stiffness[3 * a:3 * a + 3, 3 * b:3 * b + 3] += sa * sb * block
    # 3 translations + 3 rotations are exact zero-energy modes of any central-force network.
    for mode in range(6):
        vector = flat[:, mode]
        assert abs(float(vector @ stiffness @ vector)) <= 1.0e-9
    # The 6 constant-strain modes stretch bonds and carry positive energy.
    for mode in range(6, 12):
        vector = flat[:, mode]
        assert float(vector @ stiffness @ vector) > 1.0e-3
    # Determinism: a second call reproduces the sign-fixed basis bit-for-bit.
    np.testing.assert_array_equal(basis, rigid_strain_coarse_basis(pos, n_modes=12))


def test_two_level_coarse_correction_matches_exact_galerkin() -> None:
    """The additive coarse correction is exactly ``B (B^T A B)^-1 B^T r`` and zeroes the coarse residual."""
    from aleph.components.incumbent.implicit_mechanics_analytic import (
        rigid_strain_coarse_basis,
        two_level_coarse_correction,
    )

    rng = np.random.default_rng(5)
    n = 24
    pos = rng.normal(size=(n, 3))
    pos /= np.linalg.norm(pos, axis=1, keepdims=True)
    pos *= 7.5
    basis = rigid_strain_coarse_basis(pos, n_modes=12)
    root = rng.normal(size=(3 * n, 3 * n))
    operator = root @ root.T + 0.5 * np.eye(3 * n)  # SPD
    residual = rng.normal(size=(n, 3))
    correction = two_level_coarse_correction(operator, basis, residual)

    flat = basis.reshape(12, -1).T
    coarse_operator = flat.T @ operator @ flat
    reference = (flat @ np.linalg.solve(coarse_operator, flat.T @ residual.reshape(-1))).reshape(n, 3)
    np.testing.assert_allclose(correction, reference, rtol=0.0, atol=1.0e-12)
    # Galerkin defining property: the coarse-space residual after correction is zero.
    coarse_residual = flat.T @ (residual.reshape(-1) - operator @ correction.reshape(-1))
    np.testing.assert_allclose(coarse_residual, 0.0, rtol=0.0, atol=1.0e-9)


def test_warp_global_coarse_matches_two_level_reference_and_preserves_solution() -> None:
    """The device coarse-deflation preconditioner matches the numpy Galerkin oracle and leaves the fixed point.

    Device-agnostic kernel parity (CUDA on the A5000, CPU-Warp on the dev host): a preconditioner may only
    change the CG path, never the accepted solution, so ``coarse_modes=12`` must return the same displacement
    as the dense bending oracle, and its isolated coarse term must equal ``two_level_coarse_correction``.
    """
    wp = pytest.importorskip("warp")
    from aleph.components.incumbent.implicit_mechanics import ProjectedAnalyticCG
    from aleph.components.incumbent.implicit_mechanics_analytic import (
        bending_stiffness,
        rigid_strain_coarse_basis,
        two_level_coarse_correction,
    )

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    device = str(device_obj) if device_obj is not None else "cpu"

    rng = np.random.default_rng(2)
    n = 30
    pos = rng.normal(size=(n, 3))
    pos /= np.linalg.norm(pos, axis=1, keepdims=True)
    pos *= 7.5
    triples = np.array([[i, i + 1, i + 2] for i in range(0, n - 2, 3)], dtype=np.int32)
    alpha = np.full(triples.shape[0], 2.5, dtype=np.float64)
    regularization = 0.75
    rhs = rng.normal(size=(n, 3))

    cell = SimpleNamespace(
        device=device, n_total=n, n_fibers=0,
        srest_d=wp.zeros(0, dtype=wp.float64, device=device),
        n_tri=triples.shape[0], tri_d=wp.array(triples, dtype=wp.int32, device=device),
        alpha_d=wp.array(alpha, dtype=wp.float64, device=device),
        n_xl=0, steric=None, myosin=None, nucleus=None, membrane=None,
        pos_d=wp.array(pos, dtype=wp.vec3d, device=device),
    )
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    rhs_d = wp.array(rhs, dtype=wp.vec3d, device=device)
    reg_d = wp.array(np.array([regularization]), dtype=wp.float64, device=device)

    operator = regularization * np.eye(3 * n) + bending_stiffness(n, triples, alpha)  # P = I (no fibers)
    expected = np.linalg.solve(operator, rhs.reshape(-1)).reshape(n, 3)

    solver = ProjectedAnalyticCG(cell, max_iterations=128, coarse_modes=12)
    assert solver.coarse_modes == 12
    finite_d = wp.ones(1, dtype=wp.int32, device=device)
    measured = solver.solve(pos_d, rhs_d, reg_d, finite_d)
    wp.synchronize_device(device)
    np.testing.assert_allclose(measured.numpy(), expected, rtol=1.0e-9, atol=1.0e-9)
    assert int(finite_d.numpy()[0]) == 1
    assert int(solver.converged.numpy()[0]) == 1

    # Isolate the coarse term: a huge fine diagonal drives the node-Jacobi contribution to ~0, so the
    # preconditioned vector is the pure global-coarse correction and must equal the numpy Galerkin oracle.
    probe = rng.normal(size=(n, 3))
    isolate = ProjectedAnalyticCG(cell, max_iterations=1, coarse_modes=12)
    isolate._build_preconditioner(pos_d, reg_d, finite_d)
    isolate.preconditioner.fill_(1.0e300)
    isolate.r = wp.array(probe, dtype=wp.vec3d, device=device)
    isolate._precondition(pos_d, finite_d)
    wp.synchronize_device(device)
    basis = rigid_strain_coarse_basis(pos, n_modes=12)
    reference = two_level_coarse_correction(operator, basis, probe)
    np.testing.assert_allclose(isolate.z.numpy(), reference, rtol=1.0e-9, atol=1.0e-11)
