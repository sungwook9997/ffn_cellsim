"""Peskin immersed-boundary spread/interp REFERENCE — the adjoint-transfer oracle for I1a coupling.

Host-side acceptance oracle (pure NumPy, NO Warp/simulation runtime). The fluid<->solid coupling
moves node quantities to the field (SPREAD, Lagrangian->grid) and field quantities to nodes (INTERP,
grid->Lagrangian). For the FSI transfer to conserve momentum and do no spurious work, these two must
be a matched transpose pair (Peskin's ``Interp = h^d * Spread^T``). This module reimplements the exact
4-point Peskin regularized delta used by the Warp kernel ``ff/biot_fluid_warp.py::_peskin4`` /
``ibm_spread_kernel`` / ``ibm_interp_kernel`` so the adjoint identity + net-force projection can be
gated on the dev Mac.

Weights (per node ``X``, per grid cell centre ``x``; ``r_axis = (x_axis - X_axis)/h``):

    phi(r) = _peskin4(r)  in [0,2] support, sum_cells prod_axes phi = 1  (partition of unity)
    W_spread[node, cell] = (1/h^d) prod_axes phi          (deposits a density -> integrates to node value)
    W_interp[node, cell] =            prod_axes phi = h^d W_spread   (samples the field at the node)

so ``Interp = h^d Spread^T`` exactly, cell-by-cell, EVEN when the 4^d stencil is clipped at the box
edge (the clip is identical in both directions). Gates (test_i1a_ibm_gates.py):
  * partition of unity: sum of interp weights == 1 per node (interp reproduces a constant field);
  * spread conserves the total: sum_cells Spread(f) h^d == sum_nodes f (no net force created);
  * adjointness: <Spread(u), v>_grid h^d == <u, Interp(v)>_node to round-off (=> closed-system COM drift
    ~ 0 when the same operator carries the -alpha p I traction and its fluid reaction).

Build the gate against the RAW spread/interp (NOT the FSI ``num/den`` mass-lumping wrapper in
``network_warp.py``, which is nonlinear and breaks the transpose).

Reference: Peskin 2002, "The immersed boundary method", Acta Numerica 11:479.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["peskin4", "spread", "interp", "spread_matrix"]


def peskin4(r: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """The 4-point Peskin regularized delta ``phi(r)`` (support ``|r| <= 2``), matching ``_peskin4``.

    Args:
        r: Normalized offset ``(x_grid - X_node)/h`` (any shape).

    Returns:
        ``phi(r)`` with the shape of ``r``; zero for ``|r| > 2``.
    """
    a = np.abs(np.asarray(r, dtype=np.float64))
    out = np.zeros_like(a)
    m1 = a <= 1.0
    m2 = (a > 1.0) & (a <= 2.0)
    out[m1] = (3.0 - 2.0 * a[m1] + np.sqrt(1.0 + 4.0 * a[m1] - 4.0 * a[m1] ** 2)) / 8.0
    inner = np.clip(-7.0 + 12.0 * a[m2] - 4.0 * a[m2] ** 2, 0.0, None)
    out[m2] = (5.0 - 2.0 * a[m2] - np.sqrt(inner)) / 8.0
    return out


def _stencil_weights(
    nodes: npt.NDArray[np.float64], origin: npt.NDArray[np.float64], dx: float, shape: tuple[int, ...]
) -> list[tuple[int, tuple[np.ndarray, ...], np.ndarray]]:
    """For each node return ``(node_idx, cell_multi_index, prod_phi)`` over its clipped 4^d stencil."""
    dim = len(shape)
    entries: list[tuple[int, tuple[np.ndarray, ...], np.ndarray]] = []
    offsets = np.array(np.meshgrid(*([np.arange(-1, 3)] * dim), indexing="ij")).reshape(dim, -1).T
    for j, X in enumerate(nodes):
        base = np.floor((X - origin) / dx).astype(int)
        cells = base[None, :] + offsets  # (4^d, dim)
        in_box = np.all((cells >= 0) & (cells < np.array(shape)), axis=1)
        cells = cells[in_box]
        centres = origin[None, :] + cells * dx
        phi = np.prod(peskin4((centres - X[None, :]) / dx), axis=1)
        entries.append((j, tuple(cells.T), phi))
    return entries


def spread(
    values: npt.NDArray[np.float64],
    nodes: npt.NDArray[np.float64],
    origin: npt.ArrayLike,
    dx: float,
    shape: tuple[int, ...],
) -> npt.NDArray[np.float64]:
    """Spread node values to a grid density: ``field[cell] = sum_node (1/h^d) prod phi * value``.

    Args:
        values: Per-node scalar values, shape ``(N,)``.
        nodes: Node coordinates, shape ``(N, dim)``.
        origin: Grid origin (cell-centre of index 0), shape ``(dim,)``.
        dx: Uniform spacing.
        shape: Grid shape.

    Returns:
        Grid density field, shape ``shape`` (integrates to ``sum(values)``).
    """
    origin = np.asarray(origin, dtype=np.float64)
    field = np.zeros(shape, dtype=np.float64)
    inv_hd = 1.0 / dx ** len(shape)
    for j, cell_idx, phi in _stencil_weights(np.atleast_2d(nodes), origin, dx, shape):
        np.add.at(field, cell_idx, inv_hd * phi * values[j])
    return field


def interp(
    field: npt.NDArray[np.float64],
    nodes: npt.NDArray[np.float64],
    origin: npt.ArrayLike,
    dx: float,
) -> npt.NDArray[np.float64]:
    """Sample a grid field at node positions: ``value[node] = sum_cell prod phi * field[cell]``.

    Args:
        field: Grid field, shape ``grid_shape``.
        nodes: Node coordinates, shape ``(N, dim)``.
        origin: Grid origin.
        dx: Uniform spacing.

    Returns:
        Per-node sampled values, shape ``(N,)``.
    """
    origin = np.asarray(origin, dtype=np.float64)
    nodes = np.atleast_2d(nodes)
    out = np.zeros(len(nodes), dtype=np.float64)
    for j, cell_idx, phi in _stencil_weights(nodes, origin, dx, field.shape):
        out[j] = float(np.sum(phi * field[cell_idx]))
    return out


def spread_matrix(
    nodes: npt.NDArray[np.float64], origin: npt.ArrayLike, dx: float, shape: tuple[int, ...]
) -> npt.NDArray[np.float64]:
    """Dense ``(n_cells, N)`` spread matrix ``W`` (for the adjoint gate); ``field = W @ values``.

    ``interp`` is then ``value = h^d W^T @ field``, so ``W`` alone certifies the transpose relation.
    """
    origin = np.asarray(origin, dtype=np.float64)
    nodes = np.atleast_2d(nodes)
    n_cells = int(np.prod(shape))
    W = np.zeros((n_cells, len(nodes)), dtype=np.float64)
    inv_hd = 1.0 / dx ** len(shape)
    for j, cell_idx, phi in _stencil_weights(nodes, origin, dx, shape):
        flat = np.ravel_multi_index(cell_idx, shape)
        W[flat, j] = inv_hd * phi
    return W
