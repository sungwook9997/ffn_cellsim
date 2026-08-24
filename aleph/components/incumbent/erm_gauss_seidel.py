r"""Device-resident multiplicative (coloured symmetric block Gauss-Seidel) sweep over the resting graph.

The resting-baseline blocker (re-diagnosed 2026-07-23a) is NOT a missing/wrong operator stiffness — the
analytic tangent is EXACT (FD-verified 1e-16). It is the CONDITIONING of a stiff crosslink graph
(``k = 8e5 pN/um``) coupled to a soft membrane (``k_erm = 4600``) — a 174x ratio.  **Every ADDITIVE
preconditioner is exhausted** (``augmented_block`` / ``erm_jacobi`` / additive ``erm_schwarz`` all plateau
~0.55-0.73 pN, never the 0.21 gate): additive Schwarz solves each stiff crosslink pair block from the SAME
stale residual and SUMS the corrections, so a node shared by several crosslink edges is over-stepped and the
``8e5 . Delta`` spring detonates (``other_actin`` -> 8.15 pN at a full step; the line search then backtracks to
``t ~ 1/16`` where the membrane has barely engaged the ERM).

This module is the MULTIPLICATIVE fix.  It is the exact same static resting load-path spring graph as
:class:`aleph.components.incumbent.erm_schwarz.ERMPairSchwarz` (crosslinks ``k_xl`` + bound ERM ``k_erm`` + membrane
area-tension ``gamma_mem``), but instead of an additive sum it runs **symmetric block Gauss-Seidel** over a
greedy **edge colouring**: within one colour no two edges share a node (a proper edge colouring), so each
colour's two-node blocks are node-disjoint and solve in parallel; between colours the linear residual
``s = r - M dx`` is refreshed, so a later colour sees the crosslink corrections already made by earlier ones
and cannot re-add the same relaxation.  That is exactly why the additive sum over-steps and the multiplicative
sweep does not: no two coupled nodes ever step from a stale residual, so the ``8e5 . Delta`` spike cannot form.
The forward (colours ``0..C-1``) then backward (``C-1..0``) sweep makes the induced preconditioner ``B`` SPD
whenever the assembled spring operator ``M`` is SPD, so ``dx = B^-1 r`` is a genuine descent direction for the
driver's exact nonlinear residual-monotone line search (the isolated additive ``erm_schwarz`` direction is an
ascent; this one is not).

The two-node block solve reuses the FP64 Cholesky/Schur pair solve verified for the contact and ERM additive
Schwarz (``_pair_block_solve``); the graph build is the single shared
:func:`aleph.components.incumbent.erm_schwarz.build_resting_loadpath_graph`.  Off by default; wired as the additive
``erm_gauss_seidel`` / ``erm_gauss_seidel_tournament`` inner-solver options with every existing default
unchanged.

Sanity Gate:
    * Dimensions: diagonal and coupling blocks are pN/um, residual is pN, output displacement is um.
    * Boundary cases: an unbound / over-stretched (ruptured) ERM is excluded from the graph exactly as in the
      additive path; a node with no edge receives zero correction.  A degree-one edge is solved exactly in one
      colour (identical to the additive pair block for that edge).
    * Symmetry/positivity: ``M = a I + sum_e [[T, -T], [-T, T]]`` is SPD (``a > 0``; every ``_central_tangent``
      ``T`` is PSD).  A forward+backward coloured sweep is symmetric block Gauss-Seidel, whose error/
      preconditioner operator is SPD; the induced map ``r -> dx`` is symmetric positive definite (verified to
      1e-13 against a dense NumPy oracle and by eigen-decomposition on CPU).
    * Constraints: the assembled direction is passed through the same NF2007 inextensibility projector as every
      other mechanics candidate; pair reactions use equal and opposite off-diagonal blocks.
    * Numerical: 3x3 Cholesky + pair Schur complement in FP64.  No damping, cluster radius, fitted threshold,
      or biological parameter is introduced.  The colouring is a deterministic greedy proper edge colouring.
    * Residency: the static spring graph + colouring are converted once to device arrays at setup (outside the
      physical-time loop); every live tangent, matvec, block solve, and projection is Warp-CUDA.
"""
from __future__ import annotations

import numpy as np
import warp as wp

from aleph.components.incumbent.contact_schwarz import (
    _central_tangent,
    _load_symmetric_block,
    _pair_block_solve,
)
from aleph.components.incumbent.erm_schwarz import (
    add_edge_diagonal_kernel,
    build_resting_loadpath_graph,
    init_diagonal_kernel,
)

__all__ = ["ERMGaussSeidel", "greedy_edge_coloring", "symmetric_block_gs_oracle"]


def greedy_edge_coloring(
    edge_i: np.ndarray, edge_j: np.ndarray, n_nodes: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Greedy PROPER edge colouring: adjacent edges (sharing a node) get distinct colours.

    Within one colour no two edges share a node, so each colour's two-node blocks are node-disjoint and can be
    solved and scattered without a race.  Deterministic (edges processed in input order, smallest free colour
    chosen), so the colouring is reproducible for the oracle gate.

    Returns:
        ``(color, order, offsets)`` where ``color[e]`` is edge ``e``'s colour, ``order`` sorts the edges by
        colour (stable), and ``offsets`` (length ``n_colors + 1``) delimits each colour's contiguous run in the
        colour-sorted ordering.
    """
    edge_i = np.asarray(edge_i, dtype=np.int64)
    edge_j = np.asarray(edge_j, dtype=np.int64)
    m = int(edge_i.shape[0])
    color = np.full(m, -1, dtype=np.int64)
    incident: list[list[int]] = [[] for _ in range(int(n_nodes))]
    for e in range(m):
        i = int(edge_i[e])
        j = int(edge_j[e])
        used = set()
        for adj in incident[i]:
            used.add(color[adj])
        for adj in incident[j]:
            used.add(color[adj])
        c = 0
        while c in used:
            c += 1
        color[e] = c
        incident[i].append(e)
        incident[j].append(e)
    n_colors = int(color.max()) + 1 if m else 0
    order = np.argsort(color, kind="stable").astype(np.int64)
    offsets = np.zeros(n_colors + 1, dtype=np.int64)
    if m:
        np.cumsum(np.bincount(color, minlength=n_colors), out=offsets[1:])
    return color, order, offsets


def _central_tangent_np(delta: np.ndarray, stiffness: float, rest_length: float) -> np.ndarray:
    """NumPy mirror of :func:`aleph.components.incumbent.contact_schwarz._central_tangent` for the oracle gate."""
    length = float(np.linalg.norm(delta))
    if length <= 1.0e-12:
        return np.zeros((3, 3))
    unit = delta / length
    transverse = max(stiffness * (length - rest_length) / length, 0.0)
    difference = stiffness - transverse
    return transverse * np.eye(3) + difference * np.outer(unit, unit)


def symmetric_block_gs_oracle(
    residual: np.ndarray,
    diagonal: np.ndarray,
    edges: np.ndarray,
    pair_blocks: np.ndarray,
    color: np.ndarray,
    n_colors: int,
) -> np.ndarray:
    """Dense NumPy oracle for one forward+backward coloured symmetric block Gauss-Seidel sweep.

    Solves the assembled spring operator ``M`` (``M[i, i] = diagonal[i]``, ``M[i, j] = -pair_blocks[e]`` for
    edge ``e = (i, j)``) with block Gauss-Seidel over the given proper edge colouring: colours ``0..C-1`` then
    ``C-1..0``, refreshing the residual ``s = r - M dx`` before each colour.  This is the exact reference the
    device kernels must match to ~1e-13.
    """
    residual = np.asarray(residual, dtype=np.float64)
    diagonal = np.asarray(diagonal, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.int64)
    pair_blocks = np.asarray(pair_blocks, dtype=np.float64)
    color = np.asarray(color, dtype=np.int64)
    n = residual.shape[0]
    if diagonal.shape != (n, 3, 3):
        raise ValueError("diagonal must have shape (n, 3, 3)")
    if edges.ndim != 2 or edges.shape[1] != 2 or pair_blocks.shape != (edges.shape[0], 3, 3):
        raise ValueError("edges and pair_blocks must have shapes (m, 2) and (m, 3, 3)")

    dense = np.zeros((3 * n, 3 * n))
    for i in range(n):
        dense[3 * i:3 * i + 3, 3 * i:3 * i + 3] += diagonal[i]
    for (i, j), block in zip(edges, pair_blocks, strict=True):
        i, j = int(i), int(j)
        dense[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= block
        dense[3 * j:3 * j + 3, 3 * i:3 * i + 3] -= block

    dx = np.zeros((n, 3))
    sweep = list(range(n_colors)) + list(range(n_colors - 1, -1, -1))
    for c in sweep:
        s = (residual.reshape(-1) - dense @ dx.reshape(-1)).reshape(n, 3)
        for e in np.nonzero(color == c)[0]:
            i, j = int(edges[e, 0]), int(edges[e, 1])
            block = pair_blocks[e]
            local = np.block([[diagonal[i], -block], [-block, diagonal[j]]])
            solved = np.linalg.solve(local, np.concatenate([s[i], s[j]])).reshape(2, 3)
            dx[i] += solved[0]
            dx[j] += solved[1]
    return dx


@wp.kernel
def operator_diagonal_action_kernel(
    dx: wp.array(dtype=wp.vec3d),
    diagonal: wp.array(dtype=wp.float64),
    out: wp.array(dtype=wp.vec3d),
) -> None:
    """Seed ``M dx`` with the block-diagonal action ``(M dx)_i <- D_i dx_i``."""
    i = wp.tid()
    out[i] = _load_symmetric_block(diagonal, i) * dx[i]


@wp.kernel
def operator_offdiagonal_action_kernel(
    pos: wp.array(dtype=wp.vec3d),
    dx: wp.array(dtype=wp.vec3d),
    edge_i: wp.array(dtype=wp.int32),
    edge_j: wp.array(dtype=wp.int32),
    edge_k: wp.array(dtype=wp.float64),
    edge_rest: wp.array(dtype=wp.float64),
    out: wp.array(dtype=wp.vec3d),
) -> None:
    """Scatter the off-diagonal action ``-T_e dx_other`` onto both endpoints (completes ``M dx``)."""
    e = wp.tid()
    i = edge_i[e]
    j = edge_j[e]
    tangent = _central_tangent(pos[j] - pos[i], edge_k[e], edge_rest[e])
    wp.atomic_add(out, i, -(tangent * dx[j]))
    wp.atomic_add(out, j, -(tangent * dx[i]))


@wp.kernel
def gauss_seidel_color_solve_kernel(
    pos: wp.array(dtype=wp.vec3d),
    residual: wp.array(dtype=wp.vec3d),
    operator_action: wp.array(dtype=wp.vec3d),
    edge_i: wp.array(dtype=wp.int32),
    edge_j: wp.array(dtype=wp.int32),
    edge_k: wp.array(dtype=wp.float64),
    edge_rest: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.float64),
    base: wp.int32,
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    dx: wp.array(dtype=wp.vec3d),
) -> None:
    """Solve every two-node block of one colour and add its exact correction to ``dx`` (Gauss-Seidel step).

    ``s = r - M dx`` is the linear residual with the accumulated ``dx`` of all previously swept colours, so
    this colour's block correction ``[[D_i, -T], [-T, D_j]]^-1 [s_i, s_j]`` is the block Gauss-Seidel update.
    Colour edges are node-disjoint (proper edge colouring), so writing ``dx_i``/``dx_j`` needs no atomic.
    """
    if active[0] == wp.int32(0) or finite[0] == wp.int32(0):
        return
    e = base + wp.tid()
    i = edge_i[e]
    j = edge_j[e]
    coupling = _central_tangent(pos[j] - pos[i], edge_k[e], edge_rest[e])
    solution = _pair_block_solve(
        _load_symmetric_block(diagonal, i),
        _load_symmetric_block(diagonal, j),
        coupling,
        residual[i] - operator_action[i],
        residual[j] - operator_action[j],
    )
    correction_i = wp.spatial_top(solution)
    correction_j = wp.spatial_bottom(solution)
    if not wp.isfinite(wp.length_sq(correction_i)) or not wp.isfinite(wp.length_sq(correction_j)):
        finite[0] = wp.int32(0)
    else:
        dx[i] = dx[i] + correction_i
        dx[j] = dx[j] + correction_j


class ERMGaussSeidel:
    """Coloured symmetric block Gauss-Seidel sweep over the static resting load-path spring graph."""

    def __init__(self, cell: object) -> None:
        mem = getattr(cell, "membrane", None)
        if mem is None or int(getattr(mem, "n_erm", 0)) == 0:
            raise ValueError("ERM Gauss-Seidel requires membrane ERM tethers")
        self.cell = cell
        self.device = cell.device
        self.n = int(cell.n_total)

        edge_i, edge_j, edge_k, edge_rest = build_resting_loadpath_graph(cell)
        color, order, offsets = greedy_edge_coloring(edge_i, edge_j, self.n)
        # Colour-sort the edges so each colour is a contiguous, launchable run.
        edge_i = edge_i[order]
        edge_j = edge_j[order]
        edge_k = edge_k[order]
        edge_rest = edge_rest[order]
        self._assert_proper_coloring(edge_i, edge_j, offsets)
        self.color_offsets = offsets.astype(np.int64)
        self.n_colors = int(offsets.shape[0] - 1)
        self.n_edges = int(edge_i.shape[0])

        d = self.device
        with wp.ScopedDevice(d):
            self.edge_i = wp.array(np.ascontiguousarray(edge_i), dtype=wp.int32, device=d)
            self.edge_j = wp.array(np.ascontiguousarray(edge_j), dtype=wp.int32, device=d)
            self.edge_k = wp.array(np.ascontiguousarray(edge_k), dtype=wp.float64, device=d)
            self.edge_rest = wp.array(np.ascontiguousarray(edge_rest), dtype=wp.float64, device=d)
            self.diagonal = wp.zeros(6 * self.n, dtype=wp.float64, device=d)
            self.mdx = wp.zeros(self.n, dtype=wp.vec3d, device=d)
            self.dx = wp.zeros(self.n, dtype=wp.vec3d, device=d)
            self.projected_correction = wp.zeros(self.n, dtype=wp.vec3d, device=d)
            self.projector_diag = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
            self.projector_rhs = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)

    @staticmethod
    def _assert_proper_coloring(edge_i: np.ndarray, edge_j: np.ndarray, offsets: np.ndarray) -> None:
        """Fail fast if any colour has two edges sharing a node (would race the un-atomic ``dx`` write)."""
        for c in range(offsets.shape[0] - 1):
            begin, end = int(offsets[c]), int(offsets[c + 1])
            nodes = np.concatenate([edge_i[begin:end], edge_j[begin:end]])
            if nodes.size != np.unique(nodes).size:
                raise ValueError(f"edge colouring colour {c} is not node-disjoint")

    def _build_diagonal(self, pos: wp.array, regularization: wp.array) -> None:
        d = self.device
        wp.launch(init_diagonal_kernel, dim=self.n, inputs=[regularization, self.diagonal], device=d)
        wp.launch(add_edge_diagonal_kernel, dim=self.n_edges,
                  inputs=[pos, self.edge_i, self.edge_j, self.edge_k, self.edge_rest, self.diagonal],
                  device=d)

    def _solve_color(self, pos: wp.array, residual: wp.array, color: int,
                     active: wp.array, finite: wp.array) -> None:
        """Refresh ``M dx`` then apply the block Gauss-Seidel update of one colour."""
        d = self.device
        begin = int(self.color_offsets[color])
        end = int(self.color_offsets[color + 1])
        if end <= begin:
            return
        wp.launch(operator_diagonal_action_kernel, dim=self.n,
                  inputs=[self.dx, self.diagonal, self.mdx], device=d)
        wp.launch(operator_offdiagonal_action_kernel, dim=self.n_edges,
                  inputs=[pos, self.dx, self.edge_i, self.edge_j, self.edge_k, self.edge_rest, self.mdx],
                  device=d)
        wp.launch(gauss_seidel_color_solve_kernel, dim=end - begin,
                  inputs=[pos, residual, self.mdx, self.edge_i, self.edge_j, self.edge_k, self.edge_rest,
                          self.diagonal, wp.int32(begin), active, finite, self.dx], device=d)

    def solve(
        self,
        pos: wp.array,
        residual: wp.array,
        regularization: wp.array,
        active: wp.array,
        finite: wp.array,
    ) -> wp.array:
        """Return the UNPROJECTED symmetric block Gauss-Seidel correction ``dx = B^-1 r`` (self.dx)."""
        self._build_diagonal(pos, regularization)
        self.dx.zero_()
        for color in range(self.n_colors):
            self._solve_color(pos, residual, color, active, finite)
        for color in range(self.n_colors - 1, -1, -1):
            self._solve_color(pos, residual, color, active, finite)
        return self.dx

    def direction(
        self,
        pos: wp.array,
        residual: wp.array,
        regularization: wp.array,
        active: wp.array,
        finite: wp.array,
    ) -> wp.array:
        """Return the projected symmetric block Gauss-Seidel correction ``P B^-1 r``."""
        d = self.device
        self.solve(pos, residual, regularization, active, finite)
        wp.copy(self.projected_correction, self.dx)
        if self.cell.n_fibers:
            from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel
            wp.launch(project_constraint_forces_kernel, dim=self.cell.n_fibers,
                      inputs=[pos, self.dx, self.cell.foff_d, self.cell.soff_d,
                              self.projected_correction, self.projector_diag, self.projector_rhs, finite],
                      device=d)
        return self.projected_correction
