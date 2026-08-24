r"""Device-resident overlapping vector-block Schwarz over the resting load-path spring graph.

The resting-baseline blocker (diagnosed 2026-07-22d) is NOT a missing/wrong operator stiffness — the analytic
tangent is EXACT (FD-verified, ``s_operator = s_fd``, ratio 1.0000). It is the CONDITIONING of the coupled
**soft-membrane / stiff-cortex** system: the exact projected-Newton step ``(aI+PKP)^-1 PF`` overshoots (native
8×) and a scalar-damped Newton stalls at 0.78 pN (>> the 0.21 gate).

This mirrors :class:`ac.cell.contact_schwarz.ContactPairSchwarz` (which resolved the analogous stiff-WCA +
crosslink plateau), but over the **static** resting load-path spring graph rather than the live WCA graph.  The
graph is the union of three central-spring families that couple the soft membrane to the stiff cortex:

* actin **crosslinks** (``k ~ 8e5 pN/um`` — the stiffest springs; two crosslinked cortex nodes must move
  consistently),
* membrane↔cortex **ERM** tethers (``k_erm = 4600``),
* membrane **area-tension** edges (``gamma_mem``).

Each edge ``(i, j)`` owns one six-degree-of-freedom subdomain and solves its exact two-node vector block

``[[D_i, -H_ij], [-H_ij, D_j]] [dx_i, dx_j]^T = [r_i, r_j]^T``,

with ``D_i`` the full 3×3 diagonal (regularizer + every spring tangent at ``i``) and ``H_ij`` that edge's
tangent.  **Solving crosslinks as their OWN blocks — not folding them into a node diagonal — is what a first
attempt (ERM-only pairs, crosslinks in the diagonal) missed:** a crosslinked cortex node reached by one ERM
block and its crosslink partner reached by another block then stepped inconsistently and detonated the ``8e5``
spring.  Symmetric ``degree^-1/2`` partition weights make the additive operator ``sum R_e^T W_e A_e^-1 W_e R_e``
SPD without a relaxation parameter, so the one-sweep direction ``P M^-1 PF`` is a projected descent for the
driver's exact nonlinear residual-monotone line search.

STATUS (2026-07-22e, see ``docs/v2_audit/RESTING_BASELINE_DIAGNOSIS_2026-07-22e.md``): this is a PARTIAL
improvement, NOT a gate closure. In the ``erm_tournament`` (fiber-block + this) it lifts the native resting
plateau 0.776 -> 0.620 pN, but does not reach the 0.21 gate; the *isolated* ``erm_schwarz`` direction is an
ascent for the nonlinear residual. Root reason: the resting ERM-cortex nodes are mostly NON-crosslinked
(backbone-inextensibility-dominated), so a spring-pair Schwarz under-captures them — the gate needs a
backbone-aware fiber-block extended with the ERM-tethered membrane DOFs, plus a membrane tangential relaxation.
Kept as an additive, defaults-unchanged option and as the harness for that next attempt.

Sanity Gate:
    * Dimensions: diagonal and coupling blocks are pN/um, residual is pN, output displacement is um.
    * Boundary cases: an unbound / compressed / super-rupture ERM is excluded from the graph (same physics as
      ``erm_tether_force_kernel`` / ``add_erm_stiffness_kernel``); at the preloaded resting state every ERM is
      extended and sub-rupture.
    * Symmetry/positivity: each local matrix is a principal block of the SPD regularized operator ``aI+PKP``;
      symmetric partition weights preserve SPD under additive overlap; every ``_central_tangent`` is PSD.
    * Constraints: pair reactions use equal and opposite off-diagonal blocks; the assembled direction is passed
      through the same NF2007 inextensibility projector as every other mechanics candidate.
    * Numerical: 3×3 Cholesky + pair Schur complement in FP64.  No damping, cluster radius, fitted threshold,
      or biological parameter is introduced.
    * Residency: the static spring graph is converted once to device edge arrays at setup (outside the
      physical-time loop); every live tangent, block solve, overlap accumulation, and projection is Warp-CUDA.
"""
from __future__ import annotations

import numpy as np
import warp as wp

from aleph.components.incumbent.contact_schwarz import (
    _central_tangent,
    _load_symmetric_block,
    _pair_block_solve,
    _store_symmetric_block,
)

__all__ = ["ERMPairSchwarz", "build_resting_loadpath_graph"]


def build_resting_loadpath_graph(
    cell: object,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the static resting load-path spring graph as edge arrays.

    The graph is the union of the three central-spring families that couple the soft membrane to the stiff
    cortex at the resting operating point: actin **crosslinks** (``k ~ 8e5 pN/um``), bound membrane<->cortex
    **ERM** tethers (``k_erm``), and membrane **area-tension** edges (``gamma_mem``).  Every consumer (the
    additive :class:`ERMPairSchwarz` and the multiplicative
    :class:`aleph.components.incumbent.erm_gauss_seidel.ERMGaussSeidel`) reads exactly this graph, so the two share one
    topology definition and cannot drift.

    Returns:
        ``(edge_i, edge_j, edge_k, edge_rest)`` as contiguous ``int32``/``int32``/``float64``/``float64``
        arrays; ``edge_k`` is per-edge stiffness [pN/um] and ``edge_rest`` the per-edge rest length [um].
    """
    mem = getattr(cell, "membrane", None)
    if mem is None or int(getattr(mem, "n_erm", 0)) == 0:
        raise ValueError("resting load-path graph requires membrane ERM tethers")
    ei_list, ej_list, ek_list, er_list = [], [], [], []
    xl = np.asarray(cell.xl_d.numpy(), dtype=np.int64).reshape(-1, 2)
    if xl.size:
        ei_list.append(xl[:, 0]); ej_list.append(xl[:, 1])
        ek_list.append(np.asarray(cell.kxl_d.numpy(), dtype=np.float64))
        er_list.append(np.asarray(cell.r0xl_d.numpy(), dtype=np.float64))
    erm_m = np.asarray(mem.erm_m_d.numpy(), dtype=np.int64)
    erm_c = np.asarray(mem.erm_c_d.numpy(), dtype=np.int64)
    bound = np.asarray(mem.erm_bound_d.numpy(), dtype=np.int64).astype(bool)
    if bound.any():
        ei_list.append(erm_m[bound]); ej_list.append(erm_c[bound])
        ek_list.append(np.full(int(bound.sum()), float(mem.k_erm), dtype=np.float64))
        er_list.append(np.asarray(mem.erm_rest_d.numpy(), dtype=np.float64)[bound])
    if getattr(mem, "with_area_tension", False) and int(getattr(mem, "n_faces", 0)):
        faces = mem.faces_d.numpy()
        medges = np.unique(np.sort(np.concatenate(
            [faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=0), axis=1), axis=0).astype(np.int64)
        pos_build = cell.pos_d.numpy()
        ei_list.append(medges[:, 0]); ej_list.append(medges[:, 1])
        ek_list.append(np.full(medges.shape[0], float(mem.gamma_mem), dtype=np.float64))
        er_list.append(np.linalg.norm(pos_build[medges[:, 0]] - pos_build[medges[:, 1]], axis=1))
    edge_i = np.ascontiguousarray(np.concatenate(ei_list).astype(np.int32))
    edge_j = np.ascontiguousarray(np.concatenate(ej_list).astype(np.int32))
    edge_k = np.ascontiguousarray(np.concatenate(ek_list))
    edge_rest = np.ascontiguousarray(np.concatenate(er_list))
    return edge_i, edge_j, edge_k, edge_rest


@wp.func
def _atomic_add_symmetric_block(
    storage: wp.array(dtype=wp.float64), node: wp.int32, matrix: wp.mat33d,
) -> None:
    """Atomically accumulate one symmetric 3×3 tangent into the 6-float packed diagonal store."""
    base = wp.int32(6) * node
    wp.atomic_add(storage, base, matrix[0, 0])
    wp.atomic_add(storage, base + wp.int32(1), matrix[1, 1])
    wp.atomic_add(storage, base + wp.int32(2), matrix[2, 2])
    wp.atomic_add(storage, base + wp.int32(3), matrix[0, 1])
    wp.atomic_add(storage, base + wp.int32(4), matrix[0, 2])
    wp.atomic_add(storage, base + wp.int32(5), matrix[1, 2])


@wp.kernel
def init_diagonal_kernel(
    regularization: wp.array(dtype=wp.float64), diagonal: wp.array(dtype=wp.float64),
) -> None:
    """Seed every node diagonal with the isotropic omitted-family regularizer ``a I``."""
    i = wp.tid()
    a = regularization[0]
    _store_symmetric_block(diagonal, i, a * wp.identity(n=3, dtype=wp.float64))


@wp.kernel
def add_edge_diagonal_kernel(
    pos: wp.array(dtype=wp.vec3d),
    edge_i: wp.array(dtype=wp.int32),
    edge_j: wp.array(dtype=wp.int32),
    edge_k: wp.array(dtype=wp.float64),
    edge_rest: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.float64),
) -> None:
    """Add each spring tangent to both endpoint diagonals (crosslink / ERM / membrane-edge, one graph)."""
    e = wp.tid()
    i = edge_i[e]
    j = edge_j[e]
    tangent = _central_tangent(pos[j] - pos[i], edge_k[e], edge_rest[e])
    _atomic_add_symmetric_block(diagonal, i, tangent)
    _atomic_add_symmetric_block(diagonal, j, tangent)


@wp.kernel
def apply_edge_schwarz_kernel(
    pos: wp.array(dtype=wp.vec3d),
    edge_i: wp.array(dtype=wp.int32),
    edge_j: wp.array(dtype=wp.int32),
    edge_k: wp.array(dtype=wp.float64),
    edge_rest: wp.array(dtype=wp.float64),
    residual: wp.array(dtype=wp.vec3d),
    diagonal: wp.array(dtype=wp.float64),
    degree: wp.array(dtype=wp.float64),
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    correction: wp.array(dtype=wp.vec3d),
) -> None:
    """Solve and scatter every weighted two-node spring subdomain (one thread per edge)."""
    e = wp.tid()
    if active[0] == wp.int32(0) or finite[0] == wp.int32(0):
        return
    i = edge_i[e]
    j = edge_j[e]
    weight_i = wp.float64(1.0) / wp.sqrt(degree[i])
    weight_j = wp.float64(1.0) / wp.sqrt(degree[j])
    solution = _pair_block_solve(
        _load_symmetric_block(diagonal, i),
        _load_symmetric_block(diagonal, j),
        _central_tangent(pos[j] - pos[i], edge_k[e], edge_rest[e]),
        weight_i * residual[i],
        weight_j * residual[j],
    )
    correction_i = weight_i * wp.spatial_top(solution)
    correction_j = weight_j * wp.spatial_bottom(solution)
    if not wp.isfinite(wp.length_sq(correction_i)) or not wp.isfinite(wp.length_sq(correction_j)):
        finite[0] = wp.int32(0)
    else:
        wp.atomic_add(correction, i, correction_i)
        wp.atomic_add(correction, j, correction_j)


class ERMPairSchwarz:
    """Static load-path spring-graph additive Schwarz mirroring :class:`ContactPairSchwarz`."""

    def __init__(self, cell: object) -> None:
        mem = getattr(cell, "membrane", None)
        if mem is None or int(getattr(mem, "n_erm", 0)) == 0:
            raise ValueError("ERM Schwarz requires membrane ERM tethers")
        self.cell = cell
        self.device = cell.device
        self.n = int(cell.n_total)

        # Build the static resting load-path spring graph: crosslinks + bound ERM tethers + membrane edges.
        edge_i, edge_j, edge_k, edge_rest = build_resting_loadpath_graph(cell)
        degree = (np.bincount(edge_i, minlength=self.n)
                  + np.bincount(edge_j, minlength=self.n)).astype(np.float64)
        degree = np.maximum(degree, 1.0)   # an edge-free node is never solved; keep its weight finite

        d = self.device
        with wp.ScopedDevice(d):
            self.edge_i = wp.array(np.ascontiguousarray(edge_i), dtype=wp.int32, device=d)
            self.edge_j = wp.array(np.ascontiguousarray(edge_j), dtype=wp.int32, device=d)
            self.edge_k = wp.array(np.ascontiguousarray(edge_k), dtype=wp.float64, device=d)
            self.edge_rest = wp.array(np.ascontiguousarray(edge_rest), dtype=wp.float64, device=d)
            self.degree = wp.array(np.ascontiguousarray(degree), dtype=wp.float64, device=d)
            self.diagonal = wp.zeros(6 * self.n, dtype=wp.float64, device=d)
            self.correction = wp.zeros(self.n, dtype=wp.vec3d, device=d)
            self.projected_correction = wp.zeros(self.n, dtype=wp.vec3d, device=d)
            self.projector_diag = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
            self.projector_rhs = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
        self.n_edges = int(edge_i.shape[0])

    def direction(
        self,
        pos: wp.array,
        residual: wp.array,
        regularization: wp.array,
        active: wp.array,
        finite: wp.array,
    ) -> wp.array:
        """Return the projected degree-weighted additive spring-block correction ``P M^-1 r``."""
        d = self.device
        wp.launch(init_diagonal_kernel, dim=self.n, inputs=[regularization, self.diagonal], device=d)
        wp.launch(add_edge_diagonal_kernel, dim=self.n_edges,
                  inputs=[pos, self.edge_i, self.edge_j, self.edge_k, self.edge_rest, self.diagonal], device=d)
        self.correction.zero_()
        wp.launch(apply_edge_schwarz_kernel, dim=self.n_edges,
                  inputs=[pos, self.edge_i, self.edge_j, self.edge_k, self.edge_rest, residual,
                          self.diagonal, self.degree, active, finite, self.correction], device=d)
        wp.copy(self.projected_correction, self.correction)
        if self.cell.n_fibers:
            from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel
            wp.launch(project_constraint_forces_kernel, dim=self.cell.n_fibers,
                      inputs=[pos, self.correction, self.cell.foff_d, self.cell.soff_d,
                              self.projected_correction, self.projector_diag, self.projector_rhs, finite],
                      device=d)
        return self.projected_correction
