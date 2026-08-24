"""I1a immersed-boundary transfer gates — Peskin spread/interp adjointness + net-force projection.

Pure NumPy (no Warp/CUDA). Certifies the fluid<->solid coupling operators are a matched transpose
pair so the FSI momentum transfer conserves total force and does no spurious work:
  * partition of unity (interp reproduces a constant field);
  * spread conserves the total (no net force created) -> net-force projection;
  * adjointness <Spread(u),v>_grid h^d == <u,Interp(v)>_node -> closed-system COM drift ~ 0.
Gated against the RAW Peskin operators (matching ff/biot_fluid_warp.py::_peskin4), not the FSI
num/den mass-lumping wrapper.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.fluid.ibm_reference import interp, peskin4, spread, spread_matrix


def test_peskin_partition_of_unity_1d() -> None:
    """The 4-point delta sums to 1 over the integer grid for ANY sub-cell node position."""
    grid = np.arange(-3, 5).astype(float)  # integer cell centres, h=1
    for xf in np.linspace(0.0, 1.0, 11):
        assert peskin4(grid - xf).sum() == pytest.approx(1.0, abs=1e-12)


def test_peskin_second_moment_1d() -> None:
    """Peskin's delta reproduces the first moment (sum r*phi == 0): centred, no drift."""
    grid = np.arange(-3, 5).astype(float)
    for xf in np.linspace(0.0, 1.0, 6):
        w = peskin4(grid - xf)
        assert np.sum((grid - xf) * w) == pytest.approx(0.0, abs=1e-12)


def _rng_nodes(n: int, box: float, dim: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.uniform(2.0, box - 2.0, size=(n, dim))


def test_interp_reproduces_constant_field() -> None:
    """Partition of unity => interpolating a constant grid field returns that constant at every node."""
    shape = (12, 12, 12)
    dx = 0.5
    origin = np.zeros(3)
    nodes = _rng_nodes(40, 12 * dx, 3, 0)
    field = np.full(shape, 3.7)
    vals = interp(field, nodes, origin, dx)
    assert np.allclose(vals, 3.7, atol=1e-12)


def test_spread_conserves_total_force() -> None:
    """sum_cells Spread(f) h^d == sum_nodes f: the transfer creates no net force (projection gate)."""
    shape = (16, 16, 16)
    dx = 0.4
    origin = np.zeros(3)
    rng = np.random.default_rng(1)
    nodes = _rng_nodes(60, 16 * dx, 3, 1)
    f = rng.standard_normal(60)
    field = spread(f, nodes, origin, dx, shape)
    total_grid = field.sum() * dx**3
    assert total_grid == pytest.approx(f.sum(), rel=1e-12)


def test_spread_interp_are_adjoint() -> None:
    """<Spread(u), v>_grid h^d == <u, Interp(v)>_node to round-off (Interp = h^d Spread^T)."""
    shape = (10, 10, 10)
    dx = 0.5
    origin = np.zeros(3)
    rng = np.random.default_rng(2)
    nodes = _rng_nodes(30, 10 * dx, 3, 2)
    u = rng.standard_normal(30)              # node field
    v = rng.standard_normal(shape)           # grid field

    lhs = float(np.sum(spread(u, nodes, origin, dx, shape) * v) * dx**3)
    rhs = float(np.sum(u * interp(v, nodes, origin, dx)))
    assert lhs == pytest.approx(rhs, rel=1e-11, abs=1e-12)


def test_spread_matrix_transpose_is_interp() -> None:
    """The dense spread matrix W satisfies Interp(v) == h^d W^T v exactly (per-node)."""
    shape = (9, 9, 9)
    dx = 0.6
    origin = np.zeros(3)
    nodes = _rng_nodes(20, 9 * dx, 3, 3)
    W = spread_matrix(nodes, origin, dx, shape)
    rng = np.random.default_rng(4)
    v = rng.standard_normal(shape)
    from_matrix = dx**3 * (W.T @ v.ravel())
    from_interp = interp(v, nodes, origin, dx)
    assert np.allclose(from_matrix, from_interp, atol=1e-12)


def test_adjoint_holds_with_edge_clipped_stencil() -> None:
    """The transpose survives edge clipping: a node near the box face still adjoints exactly."""
    shape = (8, 8, 8)
    dx = 0.5
    origin = np.zeros(3)
    nodes = np.array([[0.2, 3.0, 3.0], [7 * dx - 0.1, 2.0, 2.0]])  # both hug a face
    rng = np.random.default_rng(5)
    u = rng.standard_normal(len(nodes))
    v = rng.standard_normal(shape)
    lhs = float(np.sum(spread(u, nodes, origin, dx, shape) * v) * dx**3)
    rhs = float(np.sum(u * interp(v, nodes, origin, dx)))
    assert lhs == pytest.approx(rhs, rel=1e-11, abs=1e-12)
