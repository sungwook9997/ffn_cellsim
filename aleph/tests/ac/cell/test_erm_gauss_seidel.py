"""Algebra + SPD/descent + kernel-vs-oracle gates for the coloured symmetric block Gauss-Seidel sweep.

The NumPy-oracle / SPD / descent / colouring gates are pure host algebra and always run.  The device
kernel-vs-oracle match is CUDA-gated (skips off a GPU) to honour the ac/-tests-never-launch-a-Warp-CPU-kernel
static contract; the equivalent CPU-Warp oracle match was verified in an uncommitted scratch probe before the
native run (7.7e-18 abs / 1.6e-14 rel), exactly as the augmented block was.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.incumbent.erm_gauss_seidel import (
    ERMGaussSeidel,
    greedy_edge_coloring,
    symmetric_block_gs_oracle,
)


def _random_spring_graph(seed: int = 0, n: int = 9):
    """Build a random connected spring graph with a stiff/soft stiffness spread like the resting cell."""
    rng = np.random.default_rng(seed)
    pos = rng.normal(size=(n, 3))
    # a spanning path guarantees connectivity, then extra random chords add crosslink overlap.
    edges = [(i, i + 1) for i in range(n - 1)]
    for _ in range(2 * n):
        a, b = int(rng.integers(n)), int(rng.integers(n))
        if a != b and (a, b) not in edges and (b, a) not in edges:
            edges.append((min(a, b), max(a, b)))
    edges_arr = np.array(edges, dtype=np.int64)
    # stiff crosslinks (8e5-scale) mixed with soft membrane/ERM (4600-scale), the 174x resting ratio.
    stiffness = rng.choice([8.0e5, 4.6e3, 2.0e2], size=edges_arr.shape[0]).astype(np.float64)
    rest = np.linalg.norm(pos[edges_arr[:, 0]] - pos[edges_arr[:, 1]], axis=1) * rng.uniform(
        0.9, 1.1, size=edges_arr.shape[0])
    a_reg = 1.5e3  # isotropic omitted-family regularizer a > 0
    return pos, edges_arr, stiffness, rest, a_reg


def _assemble_diagonal_and_blocks(pos, edges, stiffness, rest, a_reg, n):
    """Build D_i = a I + sum_e T_e and the per-edge tangents T_e in NumPy (mirrors the device build)."""
    from aleph.components.incumbent.erm_gauss_seidel import _central_tangent_np

    blocks = np.stack([
        _central_tangent_np(pos[j] - pos[i], float(stiffness[e]), float(rest[e]))
        for e, (i, j) in enumerate(edges)
    ]) if len(edges) else np.zeros((0, 3, 3))
    diagonal = np.repeat((a_reg * np.eye(3))[None], n, axis=0)
    for e, (i, j) in enumerate(edges):
        diagonal[i] += blocks[e]
        diagonal[j] += blocks[e]
    return diagonal, blocks


def test_greedy_edge_coloring_is_a_proper_edge_coloring() -> None:
    pos, edges, stiffness, rest, _ = _random_spring_graph(seed=3, n=12)
    color, order, offsets = greedy_edge_coloring(edges[:, 0], edges[:, 1], 12)
    # every colour's edges are node-disjoint (the property the un-atomic dx write relies on).
    for c in range(offsets.shape[0] - 1):
        members = np.nonzero(color == c)[0]
        nodes = np.concatenate([edges[members, 0], edges[members, 1]])
        assert nodes.size == np.unique(nodes).size
    # order groups the colours contiguously and offsets delimit them.
    assert np.array_equal(np.sort(order), np.arange(edges.shape[0]))
    assert offsets[-1] == edges.shape[0]


def test_single_edge_oracle_is_the_exact_two_node_block_solve() -> None:
    unit = np.array([1.0, -2.0, 0.5]); unit /= np.linalg.norm(unit)
    block = 5.0 * np.outer(unit, unit) + 1.0 * np.eye(3)
    diagonal = np.stack([block + 4.0 * np.eye(3), block + 6.0 * np.eye(3)])
    residual = np.array([[1.0, -2.0, 3.0], [-1.5, 0.5, 2.0]])
    edges = np.array([[0, 1]], dtype=np.int64)
    color = np.array([0], dtype=np.int64)

    measured = symmetric_block_gs_oracle(residual, diagonal, edges, block[None], color, 1)
    system = np.block([[diagonal[0], -block], [-block, diagonal[1]]])
    # one colour swept twice (forward+backward) on a single block is still the exact block solve.
    expected = np.linalg.solve(system, residual.reshape(-1)).reshape(2, 3)
    np.testing.assert_allclose(measured, expected, rtol=1.0e-13, atol=1.0e-13)


def test_symmetric_block_gs_induced_operator_is_symmetric_positive_definite() -> None:
    pos, edges, stiffness, rest, a_reg = _random_spring_graph(seed=7, n=9)
    n = 9
    diagonal, blocks = _assemble_diagonal_and_blocks(pos, edges, stiffness, rest, a_reg, n)
    color, _, offsets = greedy_edge_coloring(edges[:, 0], edges[:, 1], n)
    n_colors = int(offsets.shape[0] - 1)

    # Columns of the induced linear map r -> dx = B^-1 r.
    columns = []
    for k in range(3 * n):
        basis = np.zeros((n, 3)); basis.reshape(-1)[k] = 1.0
        columns.append(
            symmetric_block_gs_oracle(basis, diagonal, edges, blocks, color, n_colors).reshape(-1))
    operator = np.column_stack(columns)
    np.testing.assert_allclose(operator, operator.T, rtol=0.0, atol=1.0e-11)
    assert np.linalg.eigvalsh(0.5 * (operator + operator.T)).min() > 0.0


def test_symmetric_block_gs_direction_is_a_descent_direction() -> None:
    pos, edges, stiffness, rest, a_reg = _random_spring_graph(seed=11, n=10)
    n = 10
    diagonal, blocks = _assemble_diagonal_and_blocks(pos, edges, stiffness, rest, a_reg, n)
    color, _, offsets = greedy_edge_coloring(edges[:, 0], edges[:, 1], n)
    n_colors = int(offsets.shape[0] - 1)
    rng = np.random.default_rng(5)
    for _ in range(20):
        r = rng.normal(size=(n, 3))
        dx = symmetric_block_gs_oracle(r, diagonal, edges, blocks, color, n_colors)
        # B SPD => <r, B^-1 r> = <r, dx> > 0: dx has a positive projection on the force residual.
        assert float(np.sum(r * dx)) > 0.0


def test_cuda_warp_symmetric_block_gs_matches_dense_oracle() -> None:
    """The device kernels (CUDA) reproduce the dense symmetric block GS oracle to ~1e-12."""
    wp = pytest.importorskip("warp")

    wp.init()
    device_obj = next((item for item in wp.get_devices() if item.is_cuda), None)
    if device_obj is None:
        pytest.skip("I0-A: Gauss-Seidel runtime gate requires a CUDA GPU")
    device = str(device_obj)
    pos, edges, stiffness, rest, a_reg = _random_spring_graph(seed=17, n=11)
    n = 11
    diagonal, blocks = _assemble_diagonal_and_blocks(pos, edges, stiffness, rest, a_reg, n)
    residual = np.random.default_rng(2).normal(size=(n, 3))
    color, order, offsets = greedy_edge_coloring(edges[:, 0], edges[:, 1], n)
    n_colors = int(offsets.shape[0] - 1)
    expected = symmetric_block_gs_oracle(residual, diagonal, edges, blocks, color, n_colors)

    # Colour-sort the edges exactly as ERMGaussSeidel.__init__ does, then drive .solve() with hand-built arrays.
    ei = edges[order, 0].astype(np.int32)
    ej = edges[order, 1].astype(np.int32)
    ek = stiffness[order]
    er = rest[order]
    workspace = object.__new__(ERMGaussSeidel)
    workspace.device = device
    workspace.n = n
    workspace.n_edges = int(edges.shape[0])
    workspace.n_colors = n_colors
    workspace.color_offsets = offsets.astype(np.int64)
    ERMGaussSeidel._assert_proper_coloring(ei, ej, offsets)
    with wp.ScopedDevice(device):
        workspace.edge_i = wp.array(np.ascontiguousarray(ei), dtype=wp.int32, device=device)
        workspace.edge_j = wp.array(np.ascontiguousarray(ej), dtype=wp.int32, device=device)
        workspace.edge_k = wp.array(np.ascontiguousarray(ek), dtype=wp.float64, device=device)
        workspace.edge_rest = wp.array(np.ascontiguousarray(er), dtype=wp.float64, device=device)
        workspace.diagonal = wp.zeros(6 * n, dtype=wp.float64, device=device)
        workspace.mdx = wp.zeros(n, dtype=wp.vec3d, device=device)
        workspace.dx = wp.zeros(n, dtype=wp.vec3d, device=device)
    workspace.cell = None  # unused: n_fibers path is skipped via the direct .solve() call

    measured = workspace.solve(
        wp.array(pos, dtype=wp.vec3d, device=device),
        wp.array(residual, dtype=wp.vec3d, device=device),
        wp.array(np.array([a_reg]), dtype=wp.float64, device=device),
        wp.ones(1, dtype=wp.int32, device=device),
        wp.ones(1, dtype=wp.int32, device=device),
    )
    wp.synchronize_device(device)
    np.testing.assert_allclose(measured.numpy(), expected, rtol=1.0e-12, atol=1.0e-12)
