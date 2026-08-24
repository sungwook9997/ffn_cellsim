r"""Device-resident overlapping vector-block Schwarz for the live WCA contact graph.

The full-native preload plateau is dominated by nodes where a stiff radial WCA contact and a Hookean
crosslink act on the same small neighborhood.  Scalar/componentwise Jacobi loses the Cartesian coupling and
the off-node reaction.  This module instead assigns one six-degree-of-freedom subdomain to every live WCA
pair and solves its exact two-node vector block,

``[[D_i, -H_ij], [-H_ij, D_j]] [dx_i, dx_j]^T = [r_i, r_j]^T``.

``D_i`` is the full 3x3 diagonal block of the current WCA + crosslink tangent plus the already-derived omitted
family regularizer. ``H_ij`` combines the WCA tangent and every crosslink joining that same node pair.  Pair
subdomains overlap at graph junctions.  Symmetric ``degree^-1/2`` partition weights make the additive operator
``sum R_e^T W_e A_e^-1 W_e R_e`` positive definite on the represented contact nodes without a relaxation
parameter.  A degree-one contact pair is therefore solved exactly.  The resulting direction is passed through
the NF2007 tangent projector and the driver's exact nonlinear residual-monotone line search.

Sanity Gate:
    * Dimensions: diagonal and pair blocks are pN/um, residual is pN, and output displacement is um.
    * Boundary cases: nodes outside the WCA graph receive zero correction; capped WCA pairs have zero radial
      derivative exactly as in the production analytic tangent, while any coincident crosslink remains live.
    * Symmetry/positivity: every local matrix is a principal block of the SPD regularized contact operator;
      symmetric partition weights preserve SPD under additive overlap.
    * Conservation/constraints: pair reactions use equal and opposite off-diagonal blocks; the assembled
      direction is projected by the same NF2007 operator as every other mechanics candidate.
    * Numerical: 3x3 Cholesky and its pair Schur complement use FP64.  No damping, cluster radius, fitted
      threshold, or biological parameter is introduced.
    * Residency: live graph traversal, block assembly, solves, overlap accumulation, and projection are Warp
      CUDA resident.  Static crosslink topology is converted once to a device CSR during solver setup, outside
      the physical-time loop.
"""

from __future__ import annotations

import numpy as np
import warp as wp

from aleph.components.solid.steric_warp import pos_to_f32

__all__ = [
    "ContactPairSchwarz",
    "FiberContactSchwarz",
    "RigidFiberContactSchwarz",
    "additive_pair_schwarz_oracle",
]

_TWO_POW_1_6 = wp.constant(wp.float64(2.0 ** (1.0 / 6.0)))
_SQRT_EPS64 = wp.constant(wp.float64(np.sqrt(np.finfo(np.float64).eps)))
_VEC6D = wp.types.vector(length=6, dtype=wp.float64)
_VEC12D = wp.types.vector(length=12, dtype=wp.float64)
_MAT66D = wp.types.matrix(shape=(6, 6), dtype=wp.float64)


def additive_pair_schwarz_oracle(
    residual: np.ndarray,
    diagonal: np.ndarray,
    edges: np.ndarray,
    pair_blocks: np.ndarray,
) -> np.ndarray:
    """Dense NumPy oracle for the symmetric degree-weighted additive pair solve."""
    residual = np.asarray(residual, dtype=np.float64)
    diagonal = np.asarray(diagonal, dtype=np.float64)
    edges = np.asarray(edges, dtype=np.int64)
    pair_blocks = np.asarray(pair_blocks, dtype=np.float64)
    if residual.ndim != 2 or residual.shape[1] != 3:
        raise ValueError("residual must have shape (n, 3)")
    if diagonal.shape != (residual.shape[0], 3, 3):
        raise ValueError("diagonal must have shape (n, 3, 3)")
    if edges.ndim != 2 or edges.shape[1] != 2 or pair_blocks.shape != (edges.shape[0], 3, 3):
        raise ValueError("edges and pair_blocks must have shapes (m, 2) and (m, 3, 3)")
    degree = np.bincount(edges.reshape(-1), minlength=residual.shape[0]).astype(np.float64)
    correction = np.zeros_like(residual)
    for edge, block in zip(edges, pair_blocks, strict=True):
        i, j = map(int, edge)
        weights = np.array([degree[i] ** -0.5, degree[j] ** -0.5])
        local = np.block([[diagonal[i], -block], [-block, diagonal[j]]])
        local_rhs = np.concatenate([weights[0] * residual[i], weights[1] * residual[j]])
        solved = np.linalg.solve(local, local_rhs).reshape(2, 3)
        correction[i] += weights[0] * solved[0]
        correction[j] += weights[1] * solved[1]
    return correction


@wp.func
def _central_tangent(
    delta: wp.vec3d,
    stiffness: wp.float64,
    rest_length: wp.float64,
) -> wp.mat33d:
    length = wp.length(delta)
    if length <= wp.float64(1.0e-12):
        return wp.mat33d()
    unit = delta / length
    transverse = wp.max(stiffness * (length - rest_length) / length, wp.float64(0.0))
    difference = stiffness - transverse
    return transverse * wp.identity(n=3, dtype=wp.float64) + difference * wp.outer(unit, unit)


@wp.func
def _wca_tangent(
    delta: wp.vec3d,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
) -> wp.mat33d:
    zero = wp.float64(0.0)
    length = wp.length(delta)
    cutoff = _TWO_POW_1_6 * sigma
    if length <= wp.float64(1.0e-12) or length >= cutoff:
        return wp.mat33d()
    ratio = sigma / length
    sr6 = ratio * ratio * ratio
    sr6 = sr6 * sr6
    force_magnitude = (wp.float64(24.0) * epsilon / length) * (
        wp.float64(2.0) * sr6 * sr6 - sr6)
    if force_cap > zero and force_magnitude >= force_cap:
        return wp.mat33d()
    radial = (wp.float64(24.0) * epsilon / (length * length)) * (
        wp.float64(26.0) * sr6 * sr6 - wp.float64(7.0) * sr6)
    unit = delta / length
    return radial * wp.outer(unit, unit)


@wp.func
def _load_symmetric_block(storage: wp.array(dtype=wp.float64), node: wp.int32) -> wp.mat33d:
    base = wp.int32(6) * node
    xx = storage[base]
    yy = storage[base + wp.int32(1)]
    zz = storage[base + wp.int32(2)]
    xy = storage[base + wp.int32(3)]
    xz = storage[base + wp.int32(4)]
    yz = storage[base + wp.int32(5)]
    return wp.mat33d(xx, xy, xz, xy, yy, yz, xz, yz, zz)


@wp.func
def _store_symmetric_block(
    storage: wp.array(dtype=wp.float64),
    node: wp.int32,
    matrix: wp.mat33d,
) -> None:
    base = wp.int32(6) * node
    storage[base] = matrix[0, 0]
    storage[base + wp.int32(1)] = matrix[1, 1]
    storage[base + wp.int32(2)] = matrix[2, 2]
    storage[base + wp.int32(3)] = matrix[0, 1]
    storage[base + wp.int32(4)] = matrix[0, 2]
    storage[base + wp.int32(5)] = matrix[1, 2]


@wp.func
def _cholesky_solve_3x3(matrix: wp.mat33d, rhs: wp.vec3d) -> wp.vec3d:
    """Solve one guaranteed-SPD 3x3 system by explicit FP64 Cholesky."""
    l00 = wp.sqrt(matrix[0, 0])
    l10 = matrix[1, 0] / l00
    l20 = matrix[2, 0] / l00
    l11 = wp.sqrt(matrix[1, 1] - l10 * l10)
    l21 = (matrix[2, 1] - l20 * l10) / l11
    l22 = wp.sqrt(matrix[2, 2] - l20 * l20 - l21 * l21)

    y0 = rhs[0] / l00
    y1 = (rhs[1] - l10 * y0) / l11
    y2 = (rhs[2] - l20 * y0 - l21 * y1) / l22
    x2 = y2 / l22
    x1 = (y1 - l21 * x2) / l11
    x0 = (y0 - l10 * x1 - l20 * x2) / l00
    return wp.vec3d(x0, x1, x2)


@wp.func
def _pair_block_solve(
    diagonal_i: wp.mat33d,
    diagonal_j: wp.mat33d,
    coupling: wp.mat33d,
    rhs_i: wp.vec3d,
    rhs_j: wp.vec3d,
) -> wp.spatial_vectord:
    inverse_rhs_i = _cholesky_solve_3x3(diagonal_i, rhs_i)
    column0 = _cholesky_solve_3x3(
        diagonal_i, wp.vec3d(coupling[0, 0], coupling[1, 0], coupling[2, 0]))
    column1 = _cholesky_solve_3x3(
        diagonal_i, wp.vec3d(coupling[0, 1], coupling[1, 1], coupling[2, 1]))
    column2 = _cholesky_solve_3x3(
        diagonal_i, wp.vec3d(coupling[0, 2], coupling[1, 2], coupling[2, 2]))
    inverse_coupling = wp.mat33d(
        column0[0], column1[0], column2[0],
        column0[1], column1[1], column2[1],
        column0[2], column1[2], column2[2],
    )
    schur = diagonal_j - coupling * inverse_coupling
    solution_j = _cholesky_solve_3x3(schur, rhs_j + coupling * inverse_rhs_i)
    solution_i = inverse_rhs_i + inverse_coupling * solution_j
    return wp.spatial_vectord(
        solution_i[0], solution_i[1], solution_i[2],
        solution_j[0], solution_j[1], solution_j[2],
    )


@wp.kernel
def build_contact_diagonal_kernel(
    grid: wp.uint64,
    query_points: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    fiber_id: wp.array(dtype=wp.int32),
    active_nodes: wp.array(dtype=wp.int32),
    radius: wp.float32,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
    xlink_offset: wp.array(dtype=wp.int32),
    xlink_neighbor: wp.array(dtype=wp.int32),
    xlink_stiffness: wp.array(dtype=wp.float64),
    xlink_rest: wp.array(dtype=wp.float64),
    regularization: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.float64),
    contact_degree: wp.array(dtype=wp.int32),
) -> None:
    """Assemble the full Cartesian diagonal and WCA overlap degree with one owner thread per node."""
    i = wp.tid()
    identity = wp.identity(n=3, dtype=wp.float64)
    block = regularization[0] * identity
    degree = wp.int32(0)
    for edge in range(xlink_offset[i], xlink_offset[i + wp.int32(1)]):
        j = xlink_neighbor[edge]
        block = block + _central_tangent(pos[j] - pos[i], xlink_stiffness[edge], xlink_rest[edge])
    if active_nodes[i] >= wp.int32(0):
        query = wp.hash_grid_query(grid, query_points[i], radius)
        j = wp.int32(0)
        while wp.hash_grid_query_next(query, j):
            if j != i and active_nodes[j] >= wp.int32(0) and fiber_id[j] != fiber_id[i]:
                delta = pos[i] - pos[j]
                length = wp.length(delta)
                if length > wp.float64(1.0e-12) and length < _TWO_POW_1_6 * sigma:
                    block = block + _wca_tangent(delta, sigma, epsilon, force_cap)
                    degree = degree + wp.int32(1)
    _store_symmetric_block(diagonal, i, block)
    contact_degree[i] = degree


@wp.kernel
def apply_contact_pair_schwarz_kernel(
    grid: wp.uint64,
    query_points: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    residual: wp.array(dtype=wp.vec3d),
    fiber_id: wp.array(dtype=wp.int32),
    active_nodes: wp.array(dtype=wp.int32),
    radius: wp.float32,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
    xlink_offset: wp.array(dtype=wp.int32),
    xlink_neighbor: wp.array(dtype=wp.int32),
    xlink_stiffness: wp.array(dtype=wp.float64),
    xlink_rest: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.float64),
    contact_degree: wp.array(dtype=wp.int32),
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    correction: wp.array(dtype=wp.vec3d),
) -> None:
    """Solve and scatter every unique weighted two-node WCA subdomain."""
    i = wp.tid()
    if active[0] == wp.int32(0) or finite[0] == wp.int32(0) or active_nodes[i] < wp.int32(0):
        return
    query = wp.hash_grid_query(grid, query_points[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(query, j):
        if j > i and active_nodes[j] >= wp.int32(0) and fiber_id[j] != fiber_id[i]:
            delta = pos[i] - pos[j]
            length = wp.length(delta)
            if length > wp.float64(1.0e-12) and length < _TWO_POW_1_6 * sigma:
                coupling = _wca_tangent(delta, sigma, epsilon, force_cap)
                for edge in range(xlink_offset[i], xlink_offset[i + wp.int32(1)]):
                    if xlink_neighbor[edge] == j:
                        coupling = coupling + _central_tangent(
                            pos[j] - pos[i], xlink_stiffness[edge], xlink_rest[edge])
                weight_i = wp.float64(1.0) / wp.sqrt(wp.float64(contact_degree[i]))
                weight_j = wp.float64(1.0) / wp.sqrt(wp.float64(contact_degree[j]))
                solution = _pair_block_solve(
                    _load_symmetric_block(diagonal, i),
                    _load_symmetric_block(diagonal, j),
                    coupling,
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


class ContactPairSchwarz:
    """Live WCA-pair additive Schwarz workspace for one assembled cell."""

    def __init__(self, cell: object) -> None:
        if cell.steric is None:
            raise ValueError("contact Schwarz requires the production WCA force")
        self.cell = cell
        self.device = cell.device
        self.n = int(cell.n_total)

        # Static topology conversion only.  The arrays are immutable during the physical run and every live
        # position/tangent calculation remains on CUDA.
        pairs = np.asarray(cell.xl_d.numpy(), dtype=np.int32)
        stiffness = np.asarray(cell.kxl_d.numpy(), dtype=np.float64)
        rest = np.asarray(cell.r0xl_d.numpy(), dtype=np.float64)
        if pairs.size:
            source = np.concatenate([pairs[:, 0], pairs[:, 1]])
            neighbor = np.concatenate([pairs[:, 1], pairs[:, 0]])
            directed_stiffness = np.concatenate([stiffness, stiffness])
            directed_rest = np.concatenate([rest, rest])
            order = np.argsort(source, kind="stable")
            source = source[order]
            neighbor = neighbor[order]
            directed_stiffness = directed_stiffness[order]
            directed_rest = directed_rest[order]
            counts = np.bincount(source, minlength=self.n)
        else:
            neighbor = np.zeros(0, dtype=np.int32)
            directed_stiffness = np.zeros(0, dtype=np.float64)
            directed_rest = np.zeros(0, dtype=np.float64)
            counts = np.zeros(self.n, dtype=np.int64)
        offset = np.empty(self.n + 1, dtype=np.int64)
        offset[0] = 0
        np.cumsum(counts, out=offset[1:])
        if offset[-1] > np.iinfo(np.int32).max:
            raise OverflowError("directed crosslink CSR exceeds int32 indexing")

        d = self.device
        with wp.ScopedDevice(d):
            self.xlink_offset = wp.array(offset.astype(np.int32), dtype=wp.int32, device=d)
            self.xlink_neighbor = wp.array(neighbor.astype(np.int32), dtype=wp.int32, device=d)
            self.xlink_stiffness = wp.array(directed_stiffness, dtype=wp.float64, device=d)
            self.xlink_rest = wp.array(directed_rest, dtype=wp.float64, device=d)
            self.diagonal = wp.zeros(6 * self.n, dtype=wp.float64, device=d)
            self.contact_degree = wp.zeros(self.n, dtype=wp.int32, device=d)
            self.correction = wp.zeros(self.n, dtype=wp.vec3d, device=d)
            self.projected_correction = wp.zeros(self.n, dtype=wp.vec3d, device=d)
            self.projector_diag = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)
            self.projector_rhs = wp.zeros(cell.srest_d.shape[0], dtype=wp.float64, device=d)

    def direction(
        self,
        pos: wp.array,
        residual: wp.array,
        regularization: wp.array,
        active: wp.array,
        finite: wp.array,
    ) -> wp.array:
        """Return the projected degree-weighted additive pair-block correction."""
        steric = self.cell.steric
        d = self.device
        wp.launch(pos_to_f32, dim=steric.n, inputs=[pos, steric._qpts], device=d)
        steric.grid.build(steric._qpts, float(steric.r_c))
        wp.launch(
            build_contact_diagonal_kernel,
            dim=self.n,
            inputs=[
                steric.grid.id, steric._qpts, pos, steric.fiber_id, steric.active,
                wp.float32(steric.r_c), wp.float64(steric.sigma), wp.float64(steric.epsilon),
                wp.float64(steric.f_cap), self.xlink_offset, self.xlink_neighbor,
                self.xlink_stiffness, self.xlink_rest, regularization, self.diagonal,
                self.contact_degree,
            ],
            device=d,
        )
        self.correction.zero_()
        wp.launch(
            apply_contact_pair_schwarz_kernel,
            dim=self.n,
            inputs=[
                steric.grid.id, steric._qpts, pos, residual, steric.fiber_id, steric.active,
                wp.float32(steric.r_c), wp.float64(steric.sigma), wp.float64(steric.epsilon),
                wp.float64(steric.f_cap), self.xlink_offset, self.xlink_neighbor,
                self.xlink_stiffness, self.xlink_rest, self.diagonal, self.contact_degree,
                active, finite, self.correction,
            ],
            device=d,
        )
        wp.copy(self.projected_correction, self.correction)
        if self.cell.n_fibers:
            from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel

            wp.launch(
                project_constraint_forces_kernel,
                dim=self.cell.n_fibers,
                inputs=[
                    pos, self.correction, self.cell.foff_d, self.cell.soff_d,
                    self.projected_correction, self.projector_diag, self.projector_rhs, finite,
                ],
                device=d,
            )
        return self.projected_correction


@wp.kernel
def build_fiber_contact_diagonal_kernel(
    grid: wp.uint64,
    query_points: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    fiber_offset: wp.array(dtype=wp.int32),
    fiber_id: wp.array(dtype=wp.int32),
    active_nodes: wp.array(dtype=wp.int32),
    radius: wp.float32,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
    xlink_offset: wp.array(dtype=wp.int32),
    xlink_neighbor: wp.array(dtype=wp.int32),
    xlink_stiffness: wp.array(dtype=wp.float64),
    xlink_rest: wp.array(dtype=wp.float64),
    regularization: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.float64),
    contact_degree: wp.array(dtype=wp.int32),
) -> None:
    """Galerkin diagonal on the exact per-fiber translation basis."""
    fiber = wp.tid()
    begin = fiber_offset[fiber]
    end = fiber_offset[fiber + wp.int32(1)]
    identity = wp.identity(n=3, dtype=wp.float64)
    block = regularization[0] * wp.float64(end - begin) * identity
    degree = wp.int32(0)
    for i in range(begin, end):
        for edge in range(xlink_offset[i], xlink_offset[i + wp.int32(1)]):
            j = xlink_neighbor[edge]
            if fiber_id[j] != fiber:
                block = block + _central_tangent(
                    pos[j] - pos[i], xlink_stiffness[edge], xlink_rest[edge])
        query = wp.hash_grid_query(grid, query_points[i], radius)
        j = wp.int32(0)
        while wp.hash_grid_query_next(query, j):
            if j != i and active_nodes[j] >= wp.int32(0) and fiber_id[j] != fiber:
                delta = pos[i] - pos[j]
                length = wp.length(delta)
                if length > wp.float64(1.0e-12) and length < _TWO_POW_1_6 * sigma:
                    block = block + _wca_tangent(delta, sigma, epsilon, force_cap)
                    degree = degree + wp.int32(1)
    _store_symmetric_block(diagonal, fiber, block)
    contact_degree[fiber] = degree


@wp.kernel
def restrict_fiber_residual_kernel(
    residual: wp.array(dtype=wp.vec3d),
    fiber_offset: wp.array(dtype=wp.int32),
    coarse_residual: wp.array(dtype=wp.vec3d),
) -> None:
    """Apply the transpose of the piecewise-constant fiber translation basis."""
    fiber = wp.tid()
    value = wp.vec3d(0.0, 0.0, 0.0)
    for node in range(fiber_offset[fiber], fiber_offset[fiber + wp.int32(1)]):
        value = value + residual[node]
    coarse_residual[fiber] = value


@wp.kernel
def apply_fiber_contact_pair_schwarz_kernel(
    grid: wp.uint64,
    query_points: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    coarse_residual: wp.array(dtype=wp.vec3d),
    fiber_id: wp.array(dtype=wp.int32),
    active_nodes: wp.array(dtype=wp.int32),
    radius: wp.float32,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
    xlink_offset: wp.array(dtype=wp.int32),
    xlink_neighbor: wp.array(dtype=wp.int32),
    xlink_stiffness: wp.array(dtype=wp.float64),
    xlink_rest: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=wp.float64),
    contact_degree: wp.array(dtype=wp.int32),
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    coarse_correction: wp.array(dtype=wp.vec3d),
) -> None:
    """Solve each node contact as an overlapping two-fiber translation subdomain."""
    i = wp.tid()
    if active[0] == wp.int32(0) or finite[0] == wp.int32(0) or active_nodes[i] < wp.int32(0):
        return
    fiber_i = fiber_id[i]
    query = wp.hash_grid_query(grid, query_points[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(query, j):
        if j > i and active_nodes[j] >= wp.int32(0) and fiber_id[j] != fiber_i:
            delta = pos[i] - pos[j]
            length = wp.length(delta)
            if length > wp.float64(1.0e-12) and length < _TWO_POW_1_6 * sigma:
                fiber_j = fiber_id[j]
                coupling = _wca_tangent(delta, sigma, epsilon, force_cap)
                for edge in range(xlink_offset[i], xlink_offset[i + wp.int32(1)]):
                    if xlink_neighbor[edge] == j:
                        coupling = coupling + _central_tangent(
                            pos[j] - pos[i], xlink_stiffness[edge], xlink_rest[edge])
                weight_i = wp.float64(1.0) / wp.sqrt(wp.float64(contact_degree[fiber_i]))
                weight_j = wp.float64(1.0) / wp.sqrt(wp.float64(contact_degree[fiber_j]))
                solution = _pair_block_solve(
                    _load_symmetric_block(diagonal, fiber_i),
                    _load_symmetric_block(diagonal, fiber_j),
                    coupling,
                    weight_i * coarse_residual[fiber_i],
                    weight_j * coarse_residual[fiber_j],
                )
                correction_i = weight_i * wp.spatial_top(solution)
                correction_j = weight_j * wp.spatial_bottom(solution)
                if not wp.isfinite(wp.length_sq(correction_i)) or not wp.isfinite(wp.length_sq(correction_j)):
                    finite[0] = wp.int32(0)
                else:
                    wp.atomic_add(coarse_correction, fiber_i, correction_i)
                    wp.atomic_add(coarse_correction, fiber_j, correction_j)


@wp.kernel
def prolong_fiber_correction_kernel(
    coarse_correction: wp.array(dtype=wp.vec3d),
    fiber_offset: wp.array(dtype=wp.int32),
    correction: wp.array(dtype=wp.vec3d),
) -> None:
    """Prolong one uniform translation correction to each node of its fiber."""
    fiber = wp.tid()
    value = coarse_correction[fiber]
    for node in range(fiber_offset[fiber], fiber_offset[fiber + wp.int32(1)]):
        correction[node] = value


class FiberContactSchwarz(ContactPairSchwarz):
    """Overlapping vector-block Schwarz on per-fiber translation clusters."""

    def __init__(self, cell: object) -> None:
        super().__init__(cell)
        d = self.device
        with wp.ScopedDevice(d):
            self.fiber_diagonal = wp.zeros(6 * cell.n_fibers, dtype=wp.float64, device=d)
            self.fiber_contact_degree = wp.zeros(cell.n_fibers, dtype=wp.int32, device=d)
            self.coarse_residual = wp.zeros(cell.n_fibers, dtype=wp.vec3d, device=d)
            self.coarse_correction = wp.zeros(cell.n_fibers, dtype=wp.vec3d, device=d)

    def direction(
        self,
        pos: wp.array,
        residual: wp.array,
        regularization: wp.array,
        active: wp.array,
        finite: wp.array,
    ) -> wp.array:
        """Return the projected two-fiber translation-cluster correction."""
        steric = self.cell.steric
        d = self.device
        wp.launch(pos_to_f32, dim=steric.n, inputs=[pos, steric._qpts], device=d)
        steric.grid.build(steric._qpts, float(steric.r_c))
        wp.launch(
            build_fiber_contact_diagonal_kernel,
            dim=self.cell.n_fibers,
            inputs=[
                steric.grid.id, steric._qpts, pos, self.cell.foff_d, steric.fiber_id, steric.active,
                wp.float32(steric.r_c), wp.float64(steric.sigma), wp.float64(steric.epsilon),
                wp.float64(steric.f_cap), self.xlink_offset, self.xlink_neighbor,
                self.xlink_stiffness, self.xlink_rest, regularization, self.fiber_diagonal,
                self.fiber_contact_degree,
            ],
            device=d,
        )
        wp.launch(
            restrict_fiber_residual_kernel,
            dim=self.cell.n_fibers,
            inputs=[residual, self.cell.foff_d, self.coarse_residual],
            device=d,
        )
        self.coarse_correction.zero_()
        wp.launch(
            apply_fiber_contact_pair_schwarz_kernel,
            dim=self.n,
            inputs=[
                steric.grid.id, steric._qpts, pos, self.coarse_residual, steric.fiber_id,
                steric.active, wp.float32(steric.r_c), wp.float64(steric.sigma),
                wp.float64(steric.epsilon), wp.float64(steric.f_cap), self.xlink_offset,
                self.xlink_neighbor, self.xlink_stiffness, self.xlink_rest, self.fiber_diagonal,
                self.fiber_contact_degree, active, finite, self.coarse_correction,
            ],
            device=d,
        )
        self.correction.zero_()
        wp.launch(
            prolong_fiber_correction_kernel,
            dim=self.cell.n_fibers,
            inputs=[self.coarse_correction, self.cell.foff_d, self.correction],
            device=d,
        )
        wp.copy(self.projected_correction, self.correction)
        if self.cell.n_fibers:
            from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel

            wp.launch(
                project_constraint_forces_kernel,
                dim=self.cell.n_fibers,
                inputs=[
                    pos, self.correction, self.cell.foff_d, self.cell.soff_d,
                    self.projected_correction, self.projector_diag, self.projector_rhs, finite,
                ],
                device=d,
            )
        return self.projected_correction


@wp.func
def _rigid_basis_column(offset: wp.vec3d, column: wp.int32) -> wp.vec3d:
    """Return one column of ``[I, omega x offset]``."""
    if column == wp.int32(0):
        return wp.vec3d(1.0, 0.0, 0.0)
    if column == wp.int32(1):
        return wp.vec3d(0.0, 1.0, 0.0)
    if column == wp.int32(2):
        return wp.vec3d(0.0, 0.0, 1.0)
    if column == wp.int32(3):
        return wp.vec3d(0.0, -offset[2], offset[1])
    if column == wp.int32(4):
        return wp.vec3d(offset[2], 0.0, -offset[0])
    return wp.vec3d(-offset[1], offset[0], 0.0)


@wp.func
def _rigid_basis_transpose(offset: wp.vec3d, force: wp.vec3d) -> _VEC6D:
    torque = wp.cross(offset, force)
    return _VEC6D(force[0], force[1], force[2], torque[0], torque[1], torque[2])


@wp.func
def _rigid_galerkin_block(
    left_offset: wp.vec3d,
    tangent: wp.mat33d,
    right_offset: wp.vec3d,
) -> _MAT66D:
    block = _MAT66D()
    for column in range(6):
        response = tangent * _rigid_basis_column(right_offset, wp.int32(column))
        projected = _rigid_basis_transpose(left_offset, response)
        for row in range(6):
            block[row, column] = projected[row]
    return block


@wp.func
def _cholesky_solve_6x6(matrix: _MAT66D, rhs: _VEC6D) -> _VEC6D:
    factor = _MAT66D()
    for row in range(6):
        for column in range(row + 1):
            value = matrix[row, column]
            for inner in range(column):
                value = value - factor[row, inner] * factor[column, inner]
            if row == column:
                factor[row, column] = wp.sqrt(value)
            else:
                factor[row, column] = value / factor[column, column]
    work = _VEC6D()
    for row in range(6):
        value = rhs[row]
        for column in range(row):
            value = value - factor[row, column] * work[column]
        work[row] = value / factor[row, row]
    solution = _VEC6D()
    for reverse in range(6):
        row = wp.int32(5) - wp.int32(reverse)
        value = work[row]
        for column in range(row + 1, 6):
            value = value - factor[column, row] * solution[column]
        solution[row] = value / factor[row, row]
    return solution


@wp.func
def _rigid_pair_block_solve(
    diagonal_i: _MAT66D,
    diagonal_j: _MAT66D,
    coupling_ij: _MAT66D,
    rhs_i: _VEC6D,
    rhs_j: _VEC6D,
) -> _VEC12D:
    inverse_rhs_i = _cholesky_solve_6x6(diagonal_i, rhs_i)
    inverse_coupling = _MAT66D()
    for column in range(6):
        coupling_column = _VEC6D()
        for row in range(6):
            coupling_column[row] = coupling_ij[row, column]
        solved_column = _cholesky_solve_6x6(diagonal_i, coupling_column)
        for row in range(6):
            inverse_coupling[row, column] = solved_column[row]
    coupling_transpose = wp.transpose(coupling_ij)
    schur = diagonal_j - coupling_transpose * inverse_coupling
    solution_j = _cholesky_solve_6x6(
        schur, rhs_j + coupling_transpose * inverse_rhs_i)
    solution_i = inverse_rhs_i + inverse_coupling * solution_j
    result = _VEC12D()
    for component in range(6):
        result[component] = solution_i[component]
        result[component + 6] = solution_j[component]
    return result


@wp.kernel
def compute_fiber_centers_kernel(
    pos: wp.array(dtype=wp.vec3d),
    fiber_offset: wp.array(dtype=wp.int32),
    centers: wp.array(dtype=wp.vec3d),
) -> None:
    fiber = wp.tid()
    begin = fiber_offset[fiber]
    end = fiber_offset[fiber + wp.int32(1)]
    center = wp.vec3d(0.0, 0.0, 0.0)
    for node in range(begin, end):
        center = center + pos[node]
    centers[fiber] = center / wp.float64(end - begin)


@wp.kernel
def build_rigid_fiber_contact_diagonal_kernel(
    grid: wp.uint64,
    query_points: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    centers: wp.array(dtype=wp.vec3d),
    fiber_offset: wp.array(dtype=wp.int32),
    fiber_id: wp.array(dtype=wp.int32),
    active_nodes: wp.array(dtype=wp.int32),
    radius: wp.float32,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
    xlink_offset: wp.array(dtype=wp.int32),
    xlink_neighbor: wp.array(dtype=wp.int32),
    xlink_stiffness: wp.array(dtype=wp.float64),
    xlink_rest: wp.array(dtype=wp.float64),
    regularization: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=_MAT66D),
    contact_degree: wp.array(dtype=wp.int32),
) -> None:
    fiber = wp.tid()
    begin = fiber_offset[fiber]
    end = fiber_offset[fiber + wp.int32(1)]
    center = centers[fiber]
    identity = wp.identity(n=3, dtype=wp.float64)
    block = _MAT66D()
    degree = wp.int32(0)
    rotational_scale = wp.float64(0.0)
    for i in range(begin, end):
        offset_i = pos[i] - center
        rotational_scale = rotational_scale + wp.dot(offset_i, offset_i)
        block = block + regularization[0] * _rigid_galerkin_block(offset_i, identity, offset_i)
        for edge in range(xlink_offset[i], xlink_offset[i + wp.int32(1)]):
            j = xlink_neighbor[edge]
            if fiber_id[j] != fiber:
                tangent = _central_tangent(
                    pos[j] - pos[i], xlink_stiffness[edge], xlink_rest[edge])
                block = block + _rigid_galerkin_block(offset_i, tangent, offset_i)
        query = wp.hash_grid_query(grid, query_points[i], radius)
        j = wp.int32(0)
        while wp.hash_grid_query_next(query, j):
            if j != i and active_nodes[j] >= wp.int32(0) and fiber_id[j] != fiber:
                delta = pos[i] - pos[j]
                length = wp.length(delta)
                if length > wp.float64(1.0e-12) and length < _TWO_POW_1_6 * sigma:
                    tangent = _wca_tangent(delta, sigma, epsilon, force_cap)
                    block = block + _rigid_galerkin_block(offset_i, tangent, offset_i)
                    degree = degree + wp.int32(1)
    # A perfectly straight fiber has one unobservable spin about its own axis.  A machine-resolution pivot
    # keeps Cholesky defined without changing nodal motion (that spin prolongs to exactly zero on a line).
    rotational_floor = _SQRT_EPS64 * regularization[0] * rotational_scale
    for component in range(3, 6):
        block[component, component] = block[component, component] + rotational_floor
    diagonal[fiber] = block
    contact_degree[fiber] = degree


@wp.kernel
def restrict_rigid_fiber_residual_kernel(
    pos: wp.array(dtype=wp.vec3d),
    residual: wp.array(dtype=wp.vec3d),
    centers: wp.array(dtype=wp.vec3d),
    fiber_offset: wp.array(dtype=wp.int32),
    coarse_residual: wp.array(dtype=_VEC6D),
) -> None:
    fiber = wp.tid()
    value = _VEC6D()
    center = centers[fiber]
    for node in range(fiber_offset[fiber], fiber_offset[fiber + wp.int32(1)]):
        value = value + _rigid_basis_transpose(pos[node] - center, residual[node])
    coarse_residual[fiber] = value


@wp.kernel
def apply_rigid_fiber_contact_schwarz_kernel(
    grid: wp.uint64,
    query_points: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    centers: wp.array(dtype=wp.vec3d),
    coarse_residual: wp.array(dtype=_VEC6D),
    fiber_id: wp.array(dtype=wp.int32),
    active_nodes: wp.array(dtype=wp.int32),
    radius: wp.float32,
    sigma: wp.float64,
    epsilon: wp.float64,
    force_cap: wp.float64,
    xlink_offset: wp.array(dtype=wp.int32),
    xlink_neighbor: wp.array(dtype=wp.int32),
    xlink_stiffness: wp.array(dtype=wp.float64),
    xlink_rest: wp.array(dtype=wp.float64),
    diagonal: wp.array(dtype=_MAT66D),
    contact_degree: wp.array(dtype=wp.int32),
    active: wp.array(dtype=wp.int32),
    finite: wp.array(dtype=wp.int32),
    coarse_correction: wp.array(dtype=_VEC6D),
) -> None:
    i = wp.tid()
    if active[0] == wp.int32(0) or finite[0] == wp.int32(0) or active_nodes[i] < wp.int32(0):
        return
    fiber_i = fiber_id[i]
    query = wp.hash_grid_query(grid, query_points[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(query, j):
        if j > i and active_nodes[j] >= wp.int32(0) and fiber_id[j] != fiber_i:
            delta = pos[i] - pos[j]
            length = wp.length(delta)
            if length > wp.float64(1.0e-12) and length < _TWO_POW_1_6 * sigma:
                fiber_j = fiber_id[j]
                tangent = _wca_tangent(delta, sigma, epsilon, force_cap)
                for edge in range(xlink_offset[i], xlink_offset[i + wp.int32(1)]):
                    if xlink_neighbor[edge] == j:
                        tangent = tangent + _central_tangent(
                            pos[j] - pos[i], xlink_stiffness[edge], xlink_rest[edge])
                coupling = _rigid_galerkin_block(
                    pos[i] - centers[fiber_i], tangent, pos[j] - centers[fiber_j])
                weight_i = wp.float64(1.0) / wp.sqrt(wp.float64(contact_degree[fiber_i]))
                weight_j = wp.float64(1.0) / wp.sqrt(wp.float64(contact_degree[fiber_j]))
                solution = _rigid_pair_block_solve(
                    diagonal[fiber_i], diagonal[fiber_j], coupling,
                    weight_i * coarse_residual[fiber_i], weight_j * coarse_residual[fiber_j])
                correction_i = _VEC6D()
                correction_j = _VEC6D()
                for component in range(6):
                    correction_i[component] = weight_i * solution[component]
                    correction_j[component] = weight_j * solution[component + 6]
                if not wp.isfinite(wp.dot(correction_i, correction_i)) or not wp.isfinite(
                        wp.dot(correction_j, correction_j)):
                    finite[0] = wp.int32(0)
                else:
                    wp.atomic_add(coarse_correction, fiber_i, correction_i)
                    wp.atomic_add(coarse_correction, fiber_j, correction_j)


@wp.kernel
def prolong_rigid_fiber_correction_kernel(
    pos: wp.array(dtype=wp.vec3d),
    centers: wp.array(dtype=wp.vec3d),
    coarse_correction: wp.array(dtype=_VEC6D),
    fiber_offset: wp.array(dtype=wp.int32),
    correction: wp.array(dtype=wp.vec3d),
) -> None:
    fiber = wp.tid()
    value = coarse_correction[fiber]
    translation = wp.vec3d(value[0], value[1], value[2])
    rotation = wp.vec3d(value[3], value[4], value[5])
    center = centers[fiber]
    for node in range(fiber_offset[fiber], fiber_offset[fiber + wp.int32(1)]):
        correction[node] = translation + wp.cross(rotation, pos[node] - center)


class RigidFiberContactSchwarz(ContactPairSchwarz):
    """Overlapping 12-DoF pair Schwarz on rigid translation+rotation fiber clusters."""

    def __init__(self, cell: object) -> None:
        super().__init__(cell)
        d = self.device
        with wp.ScopedDevice(d):
            self.centers = wp.zeros(cell.n_fibers, dtype=wp.vec3d, device=d)
            self.rigid_diagonal = wp.zeros(cell.n_fibers, dtype=_MAT66D, device=d)
            self.rigid_contact_degree = wp.zeros(cell.n_fibers, dtype=wp.int32, device=d)
            self.rigid_residual = wp.zeros(cell.n_fibers, dtype=_VEC6D, device=d)
            self.rigid_correction = wp.zeros(cell.n_fibers, dtype=_VEC6D, device=d)

    def direction(
        self,
        pos: wp.array,
        residual: wp.array,
        regularization: wp.array,
        active: wp.array,
        finite: wp.array,
    ) -> wp.array:
        """Return the projected translation+rotation contact-cluster correction."""
        steric = self.cell.steric
        d = self.device
        wp.launch(pos_to_f32, dim=steric.n, inputs=[pos, steric._qpts], device=d)
        steric.grid.build(steric._qpts, float(steric.r_c))
        wp.launch(compute_fiber_centers_kernel, dim=self.cell.n_fibers,
                  inputs=[pos, self.cell.foff_d, self.centers], device=d)
        wp.launch(
            build_rigid_fiber_contact_diagonal_kernel,
            dim=self.cell.n_fibers,
            inputs=[
                steric.grid.id, steric._qpts, pos, self.centers, self.cell.foff_d,
                steric.fiber_id, steric.active, wp.float32(steric.r_c), wp.float64(steric.sigma),
                wp.float64(steric.epsilon), wp.float64(steric.f_cap), self.xlink_offset,
                self.xlink_neighbor, self.xlink_stiffness, self.xlink_rest, regularization,
                self.rigid_diagonal, self.rigid_contact_degree,
            ], device=d)
        wp.launch(restrict_rigid_fiber_residual_kernel, dim=self.cell.n_fibers,
                  inputs=[pos, residual, self.centers, self.cell.foff_d, self.rigid_residual], device=d)
        self.rigid_correction.zero_()
        wp.launch(
            apply_rigid_fiber_contact_schwarz_kernel,
            dim=self.n,
            inputs=[
                steric.grid.id, steric._qpts, pos, self.centers, self.rigid_residual,
                steric.fiber_id, steric.active, wp.float32(steric.r_c), wp.float64(steric.sigma),
                wp.float64(steric.epsilon), wp.float64(steric.f_cap), self.xlink_offset,
                self.xlink_neighbor, self.xlink_stiffness, self.xlink_rest, self.rigid_diagonal,
                self.rigid_contact_degree, active, finite, self.rigid_correction,
            ], device=d)
        self.correction.zero_()
        wp.launch(prolong_rigid_fiber_correction_kernel, dim=self.cell.n_fibers,
                  inputs=[pos, self.centers, self.rigid_correction, self.cell.foff_d, self.correction], device=d)
        wp.copy(self.projected_correction, self.correction)
        from aleph.components.incumbent.inner_mechanics import project_constraint_forces_kernel

        wp.launch(
            project_constraint_forces_kernel,
            dim=self.cell.n_fibers,
            inputs=[pos, self.correction, self.cell.foff_d, self.cell.soff_d,
                    self.projected_correction, self.projector_diag, self.projector_rhs, finite],
            device=d)
        return self.projected_correction
