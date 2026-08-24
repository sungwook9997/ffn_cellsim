"""Fiber-quotient inter-fiber coarse correction (Path A, matrix-free) for ``ProjectedAnalyticCG``.

This is the one coarse space the resting-cortex CG machinery lacks (design:
``aleph/docs/v2_audit/FIBER_QUOTIENT_COARSE_PLAN_2026-07-24.md``, validated CPU prototype:
``aleph/scripts/ac_fiber_quotient_coarse_proto.py``). Node-Jacobi, per-fiber block Cholesky, the
per-fiber *translation* Rayleigh coarse, and the global ``l<=2`` 12-mode coarse all resolve intra-fiber
and global-smooth modes but NONE couples distinct fibers at the coarse level, so the collective
inter-fiber near-rigid modes excited by scattered motor loads are left un-damped -> the ``max|PF| ~ 1.16``
GATE-A plateau.

The fiber-quotient coarse collapses each filament to its per-fiber rigid-body modes (3 translations + the
<=3 non-null rotations; rank ``r_f == 5`` for a straight cortical filament) and solves the crosslink-coupled
Galerkin coarse system ``A_c x = P^T r`` additively into the preconditioner. ``P`` is block-sparse (each
fine node touches only its own fiber's coarse columns) and is NEVER materialized dense (the dense ``(K,N,3)``
basis would be ~120 TB at the native 70,686 fibers; the packed per-fiber ``Q`` is ~8.5 MB).

Path A (this module): the coarse operator is applied MATRIX-FREE as ``A_c y = P^T A (P y)`` reusing the exact
fine operator ``A = aI + Pi K Pi`` (``ProjectedAnalyticCG._operator``) — provably the true ``P^T A P`` with the
correct inextensibility projector and every live connector tangent. An inner fixed-budget PCG on the
``N_c``-length coarse system, preconditioned by a per-fiber ``r_f x r_f`` block-Jacobi factor
``M_f = Q_f^T D_ext,f Q_f`` (congruence of the already-assembled ``coarse_external_diagonal``, SPD, zero
operator-applies to build), produces the coarse solution; it is prolonged back and ADDED to the smoother
output (``M^-1 = S^-1 + P A_c^-1 P^T`` stays SPD -> PCG remains valid).

This is preconditioning only: it changes the CG convergence path, never the fixed point or any residual gate.
Everything is GPU-resident Warp; the host only builds the geometry-derived prolongator once at setup (same
status as ``rigid_strain_coarse_basis``).

Sanity Gate:
    * G1 rigid fidelity / rank (CPU, ``_g1_selfcheck``): ``Q_f`` reproduces an exact per-fiber rigid motion
      to ``< 1e-9`` and a straight fiber gives ``r_f == 5`` (3 translations + 2 transverse rotations; the
      axial rotation of a collinear chain is correctly dropped as the numerical null column).
    * SPD: ``M_f = Q_f^T D_ext,f Q_f`` with ``D_ext >= a > 0`` diagonal and ``Q_f`` full column rank is SPD;
      a non-positive Cholesky pivot can only be a non-finite input and falls back to the identity block
      (that fiber then contributes an un-preconditioned inner-CG step, never a spurious correction).
    * Additivity: the coarse correction is ADDED to the smoother output, so the composite preconditioner is
      SPD and outer PCG stays valid.
    * No tuned constant: ``rank_tol`` is a numerical rank gauge relative to the matrix scale; the inner-CG
      budget is a solver compute budget, not a physics knob; there is no relaxation parameter.
"""

from __future__ import annotations

import logging
from typing import Callable

import numpy as np
import warp as wp

_LOG = logging.getLogger(__name__)

__all__ = [
    "build_fiber_rigid_prolongator",
    "fiber_rigid_restrict_kernel",
    "fiber_rigid_prolong_add_kernel",
    "build_fq_block_jacobi_kernel",
    "apply_fq_block_jacobi_kernel",
    "FiberQuotientCoarse",
    "assemble_fq_crosslink_kernel",
    "set_fq_diagonal_kernel",
    "factor_block66_kernel",
    "apply_block66_kernel",
    "fq_restrict_vec6_kernel",
    "fq_prolong_add_vec6_kernel",
    "FiberQuotientCoarsePathB",
]


# ----------------------------------------------------------------------------------------------------
# Host: per-fiber rigid prolongator (numpy; built once at setup, NOT in the hot loop).
# ----------------------------------------------------------------------------------------------------


def build_fiber_rigid_prolongator(
    pos_np: np.ndarray, foff_np: np.ndarray, *, rank_tol: float = 1.0e-9, pad_to: str | int | None = None,
) -> dict:
    r"""Build the packed per-fiber rigid-body prolongator ``P`` (the fiber-quotient coarse space).

    Each fiber ``f`` (nodes ``foff[f] .. foff[f+1]``) contributes its local rigid modes
    ``u_i = t + w x (x_i - c_f)`` — the ``(3 L_f, 6)`` candidate ``[T | R]`` about the fiber centroid — which
    is orthonormalized by column-pivoted QR; columns below ``rank_tol`` relative to the largest pivot (the
    axial rotation of a straight chain) are dropped, giving ``r_f = 5`` for a straight fiber. Different
    fibers touch disjoint node DOFs, so ``P`` is exactly block-sparse and full column rank; it is packed, not
    materialized dense.

    Args:
        pos_np: ``(N, 3)`` node positions (only rows referenced by ``foff`` are read; the actin cortex nodes).
        foff_np: ``(n_fibers + 1,)`` CSR fiber -> node offsets (``cell.foff_d``).
        rank_tol: numerical rank tolerance relative to the largest QR pivot (a gauge, not a physics constant).
        pad_to: ``None`` -> variable per-fiber column count ``r_f`` (Path A). ``"max"`` or an int -> pad every
            fiber's ``Q_f`` to a UNIFORM ``block_rank`` columns (zero-filling the dropped axial-rotation
            column), giving fixed-size coarse blocks for the Path-B BSR assembly. Zero-padded columns
            contribute nothing to restrict/prolong and get only the ``a`` regularizer on the diagonal (SPD).

    Returns:
        dict with keys:
          ``q_packed`` (float64): ``Q_f`` packed so entry
            ``q_packed[q_off[f] + (l * s_f + c) * 3 + comp] == Q_f[3 l + comp, c]`` where ``s_f`` is the
            per-fiber column stride (``r_f`` unpadded, ``block_rank`` padded).
          ``q_off`` (int32, ``n_fibers + 1``): prefix sum of ``3 L_f s_f`` (packed offsets).
          ``coarse_off`` (int32, ``n_fibers + 1``): prefix sum of ``s_f`` (coarse-DOF offsets);
            ``n_coarse = coarse_off[-1]``.
          ``rank`` (int32, ``n_fibers``): TRUE ``r_f`` per fiber (before padding).
          ``n_coarse`` (int): total coarse DOF ``N_c``.
          ``max_rank`` (int): ``max_f r_f``.
          ``block_rank`` (int): uniform column stride (== ``max_rank`` when padded, else 0 for variable).
          ``ranks_unique`` (list[int]): distinct observed ranks (diagnostic).
    """
    from scipy.linalg import qr as _qr

    pos = np.ascontiguousarray(pos_np, dtype=np.float64)
    foff = np.ascontiguousarray(foff_np, dtype=np.int64)
    n_fibers = int(foff.shape[0] - 1)

    q_keeps: list[np.ndarray] = []
    lengths = np.zeros(n_fibers, dtype=np.int64)
    ranks = np.zeros(n_fibers, dtype=np.int64)
    eye3 = np.eye(3)

    for f in range(n_fibers):
        begin = int(foff[f])
        end = int(foff[f + 1])
        x = pos[begin:end]
        length = x.shape[0]
        lengths[f] = length
        rel = x - x.mean(axis=0)
        cand = np.zeros((3 * length, 6), dtype=np.float64)
        for ax in range(3):
            trans = np.zeros((length, 3))
            trans[:, ax] = 1.0
            cand[:, ax] = trans.reshape(-1)
            cand[:, 3 + ax] = np.cross(eye3[ax], rel).reshape(-1)
        # Column-pivoted QR robustly detects and drops the near-null (axial-rotation) column.
        q, r, _ = _qr(cand, mode="economic", pivoting=True)
        diag = np.abs(np.diag(r))
        keep = diag > rank_tol * max(1.0, float(diag.max()))
        q_keeps.append(np.ascontiguousarray(q[:, keep], dtype=np.float64))  # (3L, r_f)
        ranks[f] = int(q_keeps[-1].shape[1])

    max_rank = int(ranks.max()) if n_fibers else 0
    if pad_to is None:
        block_rank = 0
    elif pad_to == "max":
        block_rank = max_rank
    else:
        block_rank = int(pad_to)
        if block_rank < max_rank:
            raise ValueError(f"pad_to={block_rank} < observed max_rank={max_rank}")

    q_blocks: list[np.ndarray] = []
    q_off = np.zeros(n_fibers + 1, dtype=np.int64)
    coarse_off = np.zeros(n_fibers + 1, dtype=np.int64)
    for f in range(n_fibers):
        q_keep = q_keeps[f]
        length = int(lengths[f])
        r_f = int(ranks[f])
        stride = block_rank if block_rank else r_f
        if stride > r_f:  # zero-pad the dropped column(s) so blocks are uniform
            q_keep = np.hstack([q_keep, np.zeros((3 * length, stride - r_f))])
        packed = np.ascontiguousarray(q_keep).reshape(length, 3, stride).transpose(0, 2, 1).reshape(-1)
        q_blocks.append(packed)
        q_off[f + 1] = q_off[f] + packed.size
        coarse_off[f + 1] = coarse_off[f] + stride

    q_packed = (np.concatenate(q_blocks) if q_blocks else np.zeros(0, dtype=np.float64))
    return {
        "q_packed": np.ascontiguousarray(q_packed, dtype=np.float64),
        "q_off": np.ascontiguousarray(q_off, dtype=np.int32),
        "coarse_off": np.ascontiguousarray(coarse_off, dtype=np.int32),
        "rank": np.ascontiguousarray(ranks, dtype=np.int32),
        "n_coarse": int(coarse_off[-1]),
        "max_rank": max_rank,
        "block_rank": int(block_rank),
        "ranks_unique": sorted(set(ranks.tolist())),
    }


# ----------------------------------------------------------------------------------------------------
# Device kernels: rigid restrict / prolong (one thread per fiber).
# ----------------------------------------------------------------------------------------------------


@wp.kernel
def fiber_rigid_restrict_kernel(
    residual: wp.array(dtype=wp.vec3d),
    foff: wp.array(dtype=wp.int32),
    q_packed: wp.array(dtype=wp.float64),
    q_off: wp.array(dtype=wp.int32),
    coarse_off: wp.array(dtype=wp.int32),
    rank: wp.array(dtype=wp.int32),
    coarse_rhs: wp.array(dtype=wp.float64),
) -> None:
    r"""``coarse_rhs[coarse_off[f]+c] = sum_l dot(Q_f[l,c], residual[foff[f]+l])`` (``P^T r``)."""
    f = wp.tid()
    r_f = rank[f]
    node_begin = foff[f]
    node_count = foff[f + 1] - node_begin
    qbase = q_off[f]
    cbase = coarse_off[f]
    for c in range(r_f):
        acc = wp.float64(0.0)
        for l in range(node_count):
            idx = qbase + (l * r_f + c) * 3
            rr = residual[node_begin + l]
            acc = acc + q_packed[idx] * rr[0] + q_packed[idx + 1] * rr[1] + q_packed[idx + 2] * rr[2]
        coarse_rhs[cbase + c] = acc


@wp.kernel
def fiber_rigid_prolong_add_kernel(
    coarse_x: wp.array(dtype=wp.float64),
    foff: wp.array(dtype=wp.int32),
    q_packed: wp.array(dtype=wp.float64),
    q_off: wp.array(dtype=wp.int32),
    coarse_off: wp.array(dtype=wp.int32),
    rank: wp.array(dtype=wp.int32),
    fine_out: wp.array(dtype=wp.vec3d),
) -> None:
    r"""``fine_out[foff[f]+l] += sum_c coarse_x[coarse_off[f]+c] * Q_f[l,c]`` (``+= P x``)."""
    f = wp.tid()
    r_f = rank[f]
    node_begin = foff[f]
    node_count = foff[f + 1] - node_begin
    qbase = q_off[f]
    cbase = coarse_off[f]
    for l in range(node_count):
        vx = wp.float64(0.0)
        vy = wp.float64(0.0)
        vz = wp.float64(0.0)
        for c in range(r_f):
            coeff = coarse_x[cbase + c]
            idx = qbase + (l * r_f + c) * 3
            vx = vx + q_packed[idx] * coeff
            vy = vy + q_packed[idx + 1] * coeff
            vz = vz + q_packed[idx + 2] * coeff
        node = node_begin + l
        old = fine_out[node]
        fine_out[node] = wp.vec3d(old[0] + vx, old[1] + vy, old[2] + vz)


# ----------------------------------------------------------------------------------------------------
# Device kernels: per-fiber block-Jacobi preconditioner for the inner coarse CG.
# ----------------------------------------------------------------------------------------------------


@wp.kernel
def build_fq_block_jacobi_kernel(
    external_diagonal: wp.array(dtype=wp.vec3d),
    foff: wp.array(dtype=wp.int32),
    q_packed: wp.array(dtype=wp.float64),
    q_off: wp.array(dtype=wp.int32),
    coarse_off: wp.array(dtype=wp.int32),
    rank: wp.array(dtype=wp.int32),
    max_rank: wp.int32,
    block_factor: wp.array(dtype=wp.float64),
) -> None:
    r"""Factor the per-fiber block-Jacobi preconditioner ``M_f = Q_f^T D_ext,f Q_f`` (tiny Cholesky).

    ``external_diagonal`` is the already-assembled ``coarse_external_diagonal`` (regularizer ``a`` plus the
    crosslink / WCA / crossbridge / LINC / ERM self couplings that do NOT cancel under a rigid fiber motion,
    with NO intra-fiber bending — exactly the diagonal the rigid modes see). The congruence
    ``Q_f^T diag(D_ext) Q_f`` is SPD because ``D_ext >= a > 0`` and ``Q_f`` is full column rank; a
    non-positive pivot can only be non-finite input and latches the identity block. One thread per fiber.
    """
    f = wp.tid()
    r_f = rank[f]
    node_begin = foff[f]
    node_count = foff[f + 1] - node_begin
    qbase = q_off[f]
    stride = max_rank * max_rank
    fbase = f * stride
    # Form the dense r_f x r_f block M into the leading (stride-max_rank) slot.
    for a in range(r_f):
        for b in range(r_f):
            s = wp.float64(0.0)
            for l in range(node_count):
                d = external_diagonal[node_begin + l]
                ia = qbase + (l * r_f + a) * 3
                ib = qbase + (l * r_f + b) * 3
                s = s + q_packed[ia] * d[0] * q_packed[ib]
                s = s + q_packed[ia + 1] * d[1] * q_packed[ib + 1]
                s = s + q_packed[ia + 2] * d[2] * q_packed[ib + 2]
            block_factor[fbase + a * max_rank + b] = s
    # In-place lower-Cholesky of the leading r_f x r_f block.
    for row in range(r_f):
        for column in range(row + 1):
            value = block_factor[fbase + row * max_rank + column]
            for inner in range(column):
                value = value - (
                    block_factor[fbase + row * max_rank + inner]
                    * block_factor[fbase + column * max_rank + inner]
                )
            if row == column:
                if value <= wp.float64(1.0e-300) or not wp.isfinite(value):
                    # Fall back to identity: this fiber's inner-CG step is then un-preconditioned.
                    for i in range(r_f):
                        for j in range(r_f):
                            if i == j:
                                block_factor[fbase + i * max_rank + j] = wp.float64(1.0)
                            else:
                                block_factor[fbase + i * max_rank + j] = wp.float64(0.0)
                    return
                block_factor[fbase + row * max_rank + column] = wp.sqrt(value)
            else:
                block_factor[fbase + row * max_rank + column] = (
                    value / block_factor[fbase + column * max_rank + column])


@wp.kernel
def apply_fq_block_jacobi_kernel(
    coarse_r: wp.array(dtype=wp.float64),
    coarse_off: wp.array(dtype=wp.int32),
    rank: wp.array(dtype=wp.int32),
    max_rank: wp.int32,
    block_factor: wp.array(dtype=wp.float64),
    coarse_z: wp.array(dtype=wp.float64),
) -> None:
    r"""Apply ``coarse_z = M_f^-1 coarse_r`` by per-fiber forward/back substitution. One thread per fiber."""
    f = wp.tid()
    r_f = rank[f]
    cbase = coarse_off[f]
    stride = max_rank * max_rank
    fbase = f * stride
    for i in range(r_f):
        coarse_z[cbase + i] = coarse_r[cbase + i]
    # Forward solve L y = rhs.
    for row in range(r_f):
        value = coarse_z[cbase + row]
        for column in range(row):
            value = value - block_factor[fbase + row * max_rank + column] * coarse_z[cbase + column]
        coarse_z[cbase + row] = value / block_factor[fbase + row * max_rank + row]
    # Back solve L^T z = y.
    for reverse in range(r_f):
        row = r_f - 1 - reverse
        value = coarse_z[cbase + row]
        for column in range(row + 1, r_f):
            value = value - block_factor[fbase + column * max_rank + row] * coarse_z[cbase + column]
        coarse_z[cbase + row] = value / block_factor[fbase + row * max_rank + row]


# ----------------------------------------------------------------------------------------------------
# Device kernels: fixed-budget inner PCG on the N_c-length coarse system (float64 vectors).
# ----------------------------------------------------------------------------------------------------


@wp.kernel
def _fq_dot_kernel(
    a: wp.array(dtype=wp.float64), b: wp.array(dtype=wp.float64), out: wp.array(dtype=wp.float64),
) -> None:
    i = wp.tid()
    wp.atomic_add(out, 0, a[i] * b[i])


@wp.kernel
def _fq_alpha_kernel(
    rz: wp.array(dtype=wp.float64), p_ap: wp.array(dtype=wp.float64), alpha: wp.array(dtype=wp.float64),
) -> None:
    denominator = p_ap[0]
    if denominator <= wp.float64(1.0e-300) or not wp.isfinite(denominator):
        alpha[0] = wp.float64(0.0)
    else:
        alpha[0] = rz[0] / denominator


@wp.kernel
def _fq_update_kernel(
    x: wp.array(dtype=wp.float64), r: wp.array(dtype=wp.float64), p: wp.array(dtype=wp.float64),
    ap: wp.array(dtype=wp.float64), alpha: wp.array(dtype=wp.float64),
) -> None:
    i = wp.tid()
    a = alpha[0]
    x[i] = x[i] + a * p[i]
    r[i] = r[i] - a * ap[i]


@wp.kernel
def _fq_beta_kernel(
    rz_new: wp.array(dtype=wp.float64), rz_old: wp.array(dtype=wp.float64), beta: wp.array(dtype=wp.float64),
) -> None:
    old = rz_old[0]
    new = rz_new[0]
    if old <= wp.float64(1.0e-300) or new <= wp.float64(0.0) or not wp.isfinite(new):
        beta[0] = wp.float64(0.0)
    else:
        beta[0] = new / old
        rz_old[0] = new


@wp.kernel
def _fq_direction_kernel(
    z: wp.array(dtype=wp.float64), beta: wp.array(dtype=wp.float64), p: wp.array(dtype=wp.float64),
) -> None:
    i = wp.tid()
    p[i] = z[i] + beta[0] * p[i]


# ----------------------------------------------------------------------------------------------------
# Driver: owns the device workspace, builds the block-Jacobi factor, applies the additive coarse.
# ----------------------------------------------------------------------------------------------------


class FiberQuotientCoarse:
    """Path-A matrix-free fiber-quotient coarse correction for ``ProjectedAnalyticCG``.

    Owns the packed prolongator ``P``, the coarse workspace, and the inner fixed-budget PCG. It is wired
    additively into the preconditioner: :meth:`build` factors the per-fiber block-Jacobi preconditioner once
    per preconditioner build (zero operator-applies) and :meth:`apply` restricts the residual, runs the inner
    PCG (coarse SpMV ``= P^T A P`` via the passed fine operator), and adds the prolonged coarse solution to the
    smoother output.
    """

    def __init__(self, cell: object, *, iterations: int = 40) -> None:
        if iterations <= 0:
            raise ValueError("fq_coarse_iterations must be positive")
        self.device = cell.device
        self.n = int(cell.n_total)
        self.n_fibers = int(cell.n_fibers)
        self.iterations = int(iterations)

        data = build_fiber_rigid_prolongator(cell.pos_d.numpy(), cell.foff_d.numpy())
        self.n_coarse = int(data["n_coarse"])
        self.max_rank = int(data["max_rank"])
        self.ranks_unique = data["ranks_unique"]

        d = self.device
        self.foff_d = cell.foff_d
        self.q_packed_d = wp.array(data["q_packed"], dtype=wp.float64, device=d)
        self.q_off_d = wp.array(data["q_off"], dtype=wp.int32, device=d)
        self.coarse_off_d = wp.array(data["coarse_off"], dtype=wp.int32, device=d)
        self.rank_d = wp.array(data["rank"], dtype=wp.int32, device=d)

        nc = max(self.n_coarse, 1)
        self.coarse_rhs = wp.zeros(nc, dtype=wp.float64, device=d)
        self.coarse_x = wp.zeros(nc, dtype=wp.float64, device=d)
        self.coarse_r = wp.zeros(nc, dtype=wp.float64, device=d)
        self.coarse_z = wp.zeros(nc, dtype=wp.float64, device=d)
        self.coarse_p = wp.zeros(nc, dtype=wp.float64, device=d)
        self.coarse_ap = wp.zeros(nc, dtype=wp.float64, device=d)
        self.block_factor = wp.zeros(
            max(self.n_fibers * self.max_rank * self.max_rank, 1), dtype=wp.float64, device=d)
        # Dedicated fine workspaces for the matrix-free SpMV (prolong -> operator -> restrict).
        self.fine = wp.zeros(self.n, dtype=wp.vec3d, device=d)
        self.fine_action = wp.zeros(self.n, dtype=wp.vec3d, device=d)
        # Inner-CG scalars (one-element device arrays; no host readback in the loop).
        self.rz = wp.zeros(1, dtype=wp.float64, device=d)
        self.rz_new = wp.zeros(1, dtype=wp.float64, device=d)
        self.p_ap = wp.zeros(1, dtype=wp.float64, device=d)
        self.alpha = wp.zeros(1, dtype=wp.float64, device=d)
        self.beta = wp.zeros(1, dtype=wp.float64, device=d)

    def build(self, external_diagonal: wp.array, finite: wp.array) -> None:
        """Factor the per-fiber block-Jacobi preconditioner (once per preconditioner build)."""
        if self.n_coarse == 0:
            return
        wp.launch(
            build_fq_block_jacobi_kernel, dim=self.n_fibers,
            inputs=[external_diagonal, self.foff_d, self.q_packed_d, self.q_off_d,
                    self.coarse_off_d, self.rank_d, wp.int32(self.max_rank), self.block_factor],
            device=self.device)

    def _spmv(self, pos: wp.array, coarse_in: wp.array, coarse_out: wp.array,
              operator: Callable, regularization: wp.array, finite: wp.array) -> None:
        """Matrix-free coarse operator ``coarse_out = P^T A (P coarse_in)`` (exact fine operator ``A``)."""
        d = self.device
        self.fine.zero_()
        wp.launch(fiber_rigid_prolong_add_kernel, dim=self.n_fibers,
                  inputs=[coarse_in, self.foff_d, self.q_packed_d, self.q_off_d,
                          self.coarse_off_d, self.rank_d, self.fine], device=d)
        operator(pos, self.fine, self.fine_action, regularization, finite)
        wp.launch(fiber_rigid_restrict_kernel, dim=self.n_fibers,
                  inputs=[self.fine_action, self.foff_d, self.q_packed_d, self.q_off_d,
                          self.coarse_off_d, self.rank_d, coarse_out], device=d)

    def apply(self, pos: wp.array, residual: wp.array, preconditioned: wp.array,
              operator: Callable, regularization: wp.array, finite: wp.array) -> None:
        r"""Add the coarse correction ``P A_c^-1 P^T residual`` into ``preconditioned`` (additive, SPD).

        Restrict ``residual`` onto the coarse space, run the fixed-budget block-Jacobi PCG on
        ``A_c x = P^T residual`` (SpMV ``= P^T A P``), and prolong-add ``x`` into ``preconditioned``.
        """
        if self.n_coarse == 0:
            return
        d = self.device
        nc = self.n_coarse
        nf = self.n_fibers
        # coarse_rhs = P^T residual; x0 = 0; r0 = coarse_rhs.
        wp.launch(fiber_rigid_restrict_kernel, dim=nf,
                  inputs=[residual, self.foff_d, self.q_packed_d, self.q_off_d,
                          self.coarse_off_d, self.rank_d, self.coarse_rhs], device=d)
        self.coarse_x.zero_()
        wp.copy(self.coarse_r, self.coarse_rhs)
        # z0 = M^-1 r0; p0 = z0; rz0 = <r0, z0>.
        wp.launch(apply_fq_block_jacobi_kernel, dim=nf,
                  inputs=[self.coarse_r, self.coarse_off_d, self.rank_d, wp.int32(self.max_rank),
                          self.block_factor, self.coarse_z], device=d)
        wp.copy(self.coarse_p, self.coarse_z)
        self.rz.zero_()
        wp.launch(_fq_dot_kernel, dim=nc, inputs=[self.coarse_r, self.coarse_z, self.rz], device=d)
        for _ in range(self.iterations):
            self._spmv(pos, self.coarse_p, self.coarse_ap, operator, regularization, finite)
            self.p_ap.zero_()
            wp.launch(_fq_dot_kernel, dim=nc, inputs=[self.coarse_p, self.coarse_ap, self.p_ap], device=d)
            wp.launch(_fq_alpha_kernel, dim=1, inputs=[self.rz, self.p_ap, self.alpha], device=d)
            wp.launch(_fq_update_kernel, dim=nc,
                      inputs=[self.coarse_x, self.coarse_r, self.coarse_p, self.coarse_ap, self.alpha],
                      device=d)
            wp.launch(apply_fq_block_jacobi_kernel, dim=nf,
                      inputs=[self.coarse_r, self.coarse_off_d, self.rank_d, wp.int32(self.max_rank),
                              self.block_factor, self.coarse_z], device=d)
            self.rz_new.zero_()
            wp.launch(_fq_dot_kernel, dim=nc, inputs=[self.coarse_r, self.coarse_z, self.rz_new], device=d)
            wp.launch(_fq_beta_kernel, dim=1, inputs=[self.rz_new, self.rz, self.beta], device=d)
            wp.launch(_fq_direction_kernel, dim=nc,
                      inputs=[self.coarse_z, self.beta, self.coarse_p], device=d)
        # preconditioned += P coarse_x.
        wp.launch(fiber_rigid_prolong_add_kernel, dim=nf,
                  inputs=[self.coarse_x, self.foff_d, self.q_packed_d, self.q_off_d,
                          self.coarse_off_d, self.rank_d, preconditioned], device=d)


# ====================================================================================================
# PATH B — explicit crosslink-weighted A_c (warp.sparse BSR) + deep block-Jacobi CG (the STRONG solve).
#
# Path A applies A_c matrix-free (P^T A P via the fine operator), so each inner-CG iteration costs a full
# fine operator-apply; only ~40-200 inner iters are affordable and the per-fiber block-Jacobi inner-CG cannot
# converge the ill-conditioned 424k inter-fiber A_c in that budget (native plateau max|PF|~1.36). Path B
# assembles A_c ONCE per preconditioner build as a crosslink-weighted 6x6-block BSR matrix (Galerkin
# A_c = a I + P^T K_xl P; rigid modes are inextensible so backbone/bending project to ~0). The BSR SpMV
# (warp.sparse.bsr_mv) is then CHEAP, so the SAME block-Jacobi inner CG can afford ~600-1500 iters — which a
# downscaled CPU study shows recovers the EXACT-coarse fine convergence (fine PCG 57 iters exact vs 60 at
# inner budget 600) where the smoother-only fine solve plateaus. (A 2-level rigid-aggregate coarse was tried
# on the same downscaled operator and did NOT converge — the sub-isostatic web has too many floppy modes for a
# single aggregate-rigid level — so the deep, cheap, contrast-agnostic block-Jacobi CG is used instead.)
#
# Uniform 6-column (padded) coarse blocks let A_c be a fixed-block-size BSR: native cortical fibers are curved
# so all six rigid modes are non-null (block_rank == max_rank == 6); a rank-5 straight fiber's dropped axial
# rotation is a zero column that only receives the ``a`` regularizer on its diagonal (SPD, harmless).
# ====================================================================================================

from warp.types import matrix as _wp_matrix  # noqa: E402
from warp.types import vector as _wp_vector  # noqa: E402

_VEC6 = _wp_vector(length=6, dtype=wp.float64)
_MAT66 = _wp_matrix(shape=(6, 6), dtype=wp.float64)
_MAT36 = _wp_matrix(shape=(3, 6), dtype=wp.float64)


@wp.func
def _fq_central_tangent33(delta: wp.vec3d, stiffness: wp.float64, rest_length: wp.float64) -> wp.mat33d:
    """Full 3x3 SPD central-spring tangent ``tau I + (k - tau) u u^T`` (matches ``_central_diagonal``)."""
    length = wp.length(delta)
    if length <= wp.float64(1.0e-12):
        return wp.mat33d(wp.float64(0.0))
    u = delta / length
    tau = wp.max(stiffness * (length - rest_length) / length, wp.float64(0.0))
    d = stiffness - tau
    return tau * wp.identity(n=3, dtype=wp.float64) + d * wp.outer(u, u)


@wp.func
def _fq_load_node_block(
    q_packed: wp.array(dtype=wp.float64), qbase: wp.int32, local: wp.int32, block_rank: wp.int32,
) -> _MAT36:
    """Load ``Q_f[local]`` (the 3 x block_rank Cartesian rows for one node) into a 3x6 matrix."""
    m = _MAT36()
    for comp in range(3):
        for c in range(block_rank):
            m[comp, c] = q_packed[qbase + (local * block_rank + c) * 3 + comp]
    return m


@wp.func
def _fq_block_congruence(qp: _MAT36, tangent: wp.mat33d, qq: _MAT36, sign: wp.float64) -> _MAT66:
    """Return ``sign * QP^T T QQ`` (6x6), the projected connector tangent block."""
    out = _MAT66()
    for a in range(6):
        for b in range(6):
            s = wp.float64(0.0)
            for i in range(3):
                for j in range(3):
                    s = s + qp[i, a] * tangent[i, j] * qq[j, b]
            out[a, b] = sign * s
    return out


@wp.kernel
def assemble_fq_crosslink_kernel(
    pos: wp.array(dtype=wp.vec3d),
    links: wp.array(dtype=wp.int32, ndim=2),
    stiffness: wp.array(dtype=wp.float64),
    rest: wp.array(dtype=wp.float64),
    node_fiber: wp.array(dtype=wp.int32),
    foff: wp.array(dtype=wp.int32),
    q_packed: wp.array(dtype=wp.float64),
    q_off: wp.array(dtype=wp.int32),
    block_rank: wp.int32,
    n_actin: wp.int32,
    rows: wp.array(dtype=wp.int32),
    cols: wp.array(dtype=wp.int32),
    vals: wp.array(dtype=_MAT66),
) -> None:
    r"""Scatter one crosslink's four projected 6x6 tangent blocks as BSR triplets (one thread per link).

    For link ``e = (i on fiber f, j on fiber g)`` with 3x3 tangent ``T``: ``+Q_f[i]^T T Q_f[i]`` -> ``(f,f)``,
    ``+Q_g[j]^T T Q_g[j]`` -> ``(g,g)``, ``-Q_f[i]^T T Q_g[j]`` -> ``(f,g)``, its transpose -> ``(g,f)``.
    """
    t = wp.tid()
    base = t * 4
    i = links[t, 0]
    j = links[t, 1]
    if i < 0 or j < 0 or i >= n_actin or j >= n_actin:
        for k in range(4):  # emit harmless zero blocks at (0,0) to keep the triplet arrays valid
            rows[base + k] = 0
            cols[base + k] = 0
            vals[base + k] = _MAT66()
        return
    f = node_fiber[i]
    g = node_fiber[j]
    lp = i - foff[f]
    lq = j - foff[g]
    tangent = _fq_central_tangent33(pos[j] - pos[i], stiffness[t], rest[t])
    qp = _fq_load_node_block(q_packed, q_off[f], lp, block_rank)
    qq = _fq_load_node_block(q_packed, q_off[g], lq, block_rank)
    sff = _fq_block_congruence(qp, tangent, qp, wp.float64(1.0))
    sgg = _fq_block_congruence(qq, tangent, qq, wp.float64(1.0))
    sfg = _fq_block_congruence(qp, tangent, qq, wp.float64(-1.0))
    rows[base + 0] = f; cols[base + 0] = f; vals[base + 0] = sff
    rows[base + 1] = g; cols[base + 1] = g; vals[base + 1] = sgg
    rows[base + 2] = f; cols[base + 2] = g; vals[base + 2] = sfg
    rows[base + 3] = g; cols[base + 3] = f; vals[base + 3] = wp.transpose(sfg)


@wp.kernel
def set_fq_diagonal_kernel(
    reg: wp.array(dtype=wp.float64), base_off: wp.int32,
    rows: wp.array(dtype=wp.int32), cols: wp.array(dtype=wp.int32), vals: wp.array(dtype=_MAT66),
) -> None:
    """Emit the ``a I_6`` regularizer block on each fiber's diagonal (keeps A_c SPD, one thread per fiber)."""
    f = wp.tid()
    a = reg[0]
    m = _MAT66()
    for c in range(6):
        m[c, c] = a
    idx = base_off + f
    rows[idx] = f
    cols[idx] = f
    vals[idx] = m


@wp.kernel
def factor_block66_kernel(
    diag: wp.array(dtype=_MAT66), factor: wp.array(dtype=_MAT66), finite: wp.array(dtype=wp.int32),
) -> None:
    """Lower-Cholesky factor of each fiber's 6x6 A_c diagonal block (block-Jacobi preconditioner)."""
    f = wp.tid()
    lo = diag[f]  # copy; the Cholesky overwrites the lower triangle in place
    for row in range(6):
        for column in range(row + 1):
            value = lo[row, column]
            for inner in range(column):
                value = value - lo[row, inner] * lo[column, inner]
            if row == column:
                if value <= wp.float64(1.0e-300) or not wp.isfinite(value):
                    factor[f] = wp.identity(n=6, dtype=wp.float64)  # fall back to identity for this fiber
                    return
                lo[row, column] = wp.sqrt(value)
            else:
                lo[row, column] = value / lo[column, column]
    factor[f] = lo


@wp.kernel
def apply_block66_kernel(
    coarse_r: wp.array(dtype=_VEC6), factor: wp.array(dtype=_MAT66), coarse_z: wp.array(dtype=_VEC6),
) -> None:
    """Apply ``z = M_f^-1 r`` by 6x6 forward/back substitution (one thread per fiber)."""
    f = wp.tid()
    lo = factor[f]
    r = coarse_r[f]
    y = _VEC6()
    for row in range(6):  # forward: L y = r
        value = r[row]
        for column in range(row):
            value = value - lo[row, column] * y[column]
        y[row] = value / lo[row, row]
    z = _VEC6()
    for reverse in range(6):  # back: L^T z = y
        row = 5 - reverse
        value = y[row]
        for column in range(row + 1, 6):
            value = value - lo[column, row] * z[column]
        z[row] = value / lo[row, row]
    coarse_z[f] = z


@wp.kernel
def fq_restrict_vec6_kernel(
    residual: wp.array(dtype=wp.vec3d), foff: wp.array(dtype=wp.int32),
    q_packed: wp.array(dtype=wp.float64), q_off: wp.array(dtype=wp.int32), block_rank: wp.int32,
    coarse_rhs: wp.array(dtype=_VEC6),
) -> None:
    r"""``coarse_rhs[f][c] = sum_l dot(Q_f[l,c], residual[foff[f]+l])`` as a vec6 (one thread per fiber)."""
    f = wp.tid()
    node_begin = foff[f]
    node_count = foff[f + 1] - node_begin
    qbase = q_off[f]
    v = _VEC6()
    for c in range(block_rank):
        acc = wp.float64(0.0)
        for l in range(node_count):
            idx = qbase + (l * block_rank + c) * 3
            rr = residual[node_begin + l]
            acc = acc + q_packed[idx] * rr[0] + q_packed[idx + 1] * rr[1] + q_packed[idx + 2] * rr[2]
        v[c] = acc
    coarse_rhs[f] = v


@wp.kernel
def fq_prolong_add_vec6_kernel(
    coarse_x: wp.array(dtype=_VEC6), foff: wp.array(dtype=wp.int32),
    q_packed: wp.array(dtype=wp.float64), q_off: wp.array(dtype=wp.int32), block_rank: wp.int32,
    fine_out: wp.array(dtype=wp.vec3d),
) -> None:
    r"""``fine_out[foff[f]+l] += sum_c coarse_x[f][c] * Q_f[l,c]`` (one thread per fiber)."""
    f = wp.tid()
    node_begin = foff[f]
    node_count = foff[f + 1] - node_begin
    qbase = q_off[f]
    cx = coarse_x[f]
    for l in range(node_count):
        vx = wp.float64(0.0)
        vy = wp.float64(0.0)
        vz = wp.float64(0.0)
        for c in range(block_rank):
            coeff = cx[c]
            idx = qbase + (l * block_rank + c) * 3
            vx = vx + q_packed[idx] * coeff
            vy = vy + q_packed[idx + 1] * coeff
            vz = vz + q_packed[idx + 2] * coeff
        node = node_begin + l
        old = fine_out[node]
        fine_out[node] = wp.vec3d(old[0] + vx, old[1] + vy, old[2] + vz)


@wp.kernel
def _fq6_dot_kernel(a: wp.array(dtype=_VEC6), b: wp.array(dtype=_VEC6),
                    out: wp.array(dtype=wp.float64)) -> None:
    i = wp.tid()
    wp.atomic_add(out, 0, wp.dot(a[i], b[i]))


@wp.kernel
def _fq6_update_kernel(x: wp.array(dtype=_VEC6), r: wp.array(dtype=_VEC6), p: wp.array(dtype=_VEC6),
                       ap: wp.array(dtype=_VEC6), alpha: wp.array(dtype=wp.float64)) -> None:
    i = wp.tid()
    a = alpha[0]
    x[i] = x[i] + a * p[i]
    r[i] = r[i] - a * ap[i]


@wp.kernel
def _fq6_direction_kernel(z: wp.array(dtype=_VEC6), beta: wp.array(dtype=wp.float64),
                          p: wp.array(dtype=_VEC6)) -> None:
    i = wp.tid()
    p[i] = z[i] + beta[0] * p[i]


@wp.kernel
def _fq6_rsub_kernel(rhs: wp.array(dtype=_VEC6), ax: wp.array(dtype=_VEC6),
                     r: wp.array(dtype=_VEC6)) -> None:
    """``r = rhs - ax`` (the coarse Richardson residual)."""
    i = wp.tid()
    r[i] = rhs[i] - ax[i]


@wp.kernel
def _fq6_axpy_dev_kernel(x: wp.array(dtype=_VEC6), z: wp.array(dtype=_VEC6),
                         omega: wp.array(dtype=wp.float64)) -> None:
    """``x += omega * z`` with a DEVICE-resident damping scalar (no host readback in the loop)."""
    i = wp.tid()
    x[i] = x[i] + omega[0] * z[i]


@wp.kernel
def _fq6_matvec_diag_kernel(diag: wp.array(dtype=_MAT66), v: wp.array(dtype=_VEC6),
                            out: wp.array(dtype=_VEC6)) -> None:
    """``out[f] = M_f v[f]`` with ``M_f`` the fiber's 6x6 A_c diagonal block (for the M-inner-product Rayleigh)."""
    f = wp.tid()
    out[f] = diag[f] * v[f]


@wp.kernel
def _fq6_normalize_kernel(z: wp.array(dtype=_VEC6), nrm2: wp.array(dtype=wp.float64),
                          out: wp.array(dtype=_VEC6)) -> None:
    """``out = z / ||z||_2`` (guarded; a zero vector is left zero) — the power-iteration renormalization."""
    f = wp.tid()
    n = nrm2[0]
    if n > wp.float64(1.0e-300):
        out[f] = z[f] / wp.sqrt(n)
    else:
        out[f] = z[f]


@wp.kernel
def _fq_set_omega_kernel(num: wp.array(dtype=wp.float64), den: wp.array(dtype=wp.float64),
                         safety: wp.float64, omega: wp.array(dtype=wp.float64)) -> None:
    r"""``omega = safety / lambda_max`` with the M-Rayleigh estimate ``lambda_max = num/den = <v,A_c v>/<v,M v>``.

    Rayleigh quotient in the M-inner-product converges to ``lambda_max(M^-1 A_c)`` from below, so
    ``safety < 1`` keeps ``omega < 1/lambda_max < 2/lambda_max`` — the unconditional-SPD Richardson range.
    """
    n = num[0]
    d = den[0]
    if n > wp.float64(1.0e-300) and d > wp.float64(1.0e-300) and wp.isfinite(n) and wp.isfinite(d):
        omega[0] = safety * d / n
    else:
        omega[0] = wp.float64(0.0)


class FiberQuotientCoarsePathB:
    """Path-B strong fiber-quotient coarse: BSR-assembled A_c + deep block-Jacobi CG.

    Assembles the crosslink-weighted Galerkin coarse operator ``A_c = a I + P^T K_xl P`` as a 6x6-block BSR
    matrix once per preconditioner build, then approximates ``A_c^-1 P^T r`` with a fixed, generous budget of
    block-Jacobi conjugate gradient (cheap ``bsr_mv`` SpMV) and prolongs the result additively into the
    smoother output. Preconditioner-only; SPD; default OFF.
    """

    def __init__(self, cell: object, *, inner_iterations: int = 1000) -> None:
        if inner_iterations <= 0:
            raise ValueError("fq_coarse_iterations must be positive")
        import warp.sparse as wsp
        self._wsp = wsp
        self.device = cell.device
        self.n = int(cell.n_total)
        self.n_actin = int(cell.n_actin)
        self.n_fibers = int(cell.n_fibers)
        self.iterations = int(inner_iterations)

        # The full rigid space is 6-D, so the BSR block size is FIXED at 6 (curved native fibers are all
        # rank 6; a straight fiber's dropped axial rotation is a zero-padded 6th column, SPD via the a I_6
        # diagonal). Padding to 6 (not max_rank) keeps the block type _MAT66/_VEC6 valid for every config.
        data = build_fiber_rigid_prolongator(cell.pos_d.numpy(), cell.foff_d.numpy(), pad_to=6)
        self.block_rank = int(data["block_rank"])  # == 6
        self.max_rank = int(data["max_rank"])
        self.ranks_unique = data["ranks_unique"]
        self.n_coarse = int(data["n_coarse"])
        if self.max_rank > 6:
            raise ValueError(f"fiber rigid rank {self.max_rank} > 6 is impossible")

        d = self.device
        self.foff_d = cell.foff_d
        self.node_fiber_d = cell.node_fiber_d
        self.q_packed_d = wp.array(data["q_packed"], dtype=wp.float64, device=d)
        self.q_off_d = wp.array(data["q_off"], dtype=wp.int32, device=d)
        # crosslink connectivity (the inter-fiber coupling that builds A_c's off-diagonal graph)
        self.xl_d = getattr(cell, "xl_d", None)
        self.kxl_d = getattr(cell, "kxl_d", None)
        self.r0xl_d = getattr(cell, "r0xl_d", None)
        self.n_xl = int(getattr(cell, "n_xl", 0))

        nf = self.n_fibers
        self.n_trip = 4 * self.n_xl + nf
        self.rows_d = wp.zeros(self.n_trip, dtype=wp.int32, device=d)
        self.cols_d = wp.zeros(self.n_trip, dtype=wp.int32, device=d)
        self.vals_d = wp.zeros(self.n_trip, dtype=_MAT66, device=d)
        self.A_c = wsp.bsr_zeros(nf, nf, _MAT66, device=d)
        self.diag_blocks = wp.zeros(nf, dtype=_MAT66, device=d)
        self.bj_factor = wp.zeros(nf, dtype=_MAT66, device=d)
        # vec6 coarse workspaces (one vec6 per fiber; padded block_rank == 6 components)
        self.coarse_rhs = wp.zeros(nf, dtype=_VEC6, device=d)
        self.coarse_x = wp.zeros(nf, dtype=_VEC6, device=d)
        self.coarse_r = wp.zeros(nf, dtype=_VEC6, device=d)
        self.coarse_z = wp.zeros(nf, dtype=_VEC6, device=d)
        self.coarse_p = wp.zeros(nf, dtype=_VEC6, device=d)
        self.coarse_ap = wp.zeros(nf, dtype=_VEC6, device=d)
        self.rz = wp.zeros(1, dtype=wp.float64, device=d)
        self.rz_new = wp.zeros(1, dtype=wp.float64, device=d)
        self.p_ap = wp.zeros(1, dtype=wp.float64, device=d)
        self.alpha = wp.zeros(1, dtype=wp.float64, device=d)
        self.beta = wp.zeros(1, dtype=wp.float64, device=d)
        # Power-iteration workspace for the block-Jacobi Richardson damping estimate (apply_richardson).
        self.pw_v = wp.zeros(nf, dtype=_VEC6, device=d)
        self.pw_av = wp.zeros(nf, dtype=_VEC6, device=d)
        self.pw_mv = wp.zeros(nf, dtype=_VEC6, device=d)
        self.pw_z = wp.zeros(nf, dtype=_VEC6, device=d)
        self.pw_num = wp.zeros(1, dtype=wp.float64, device=d)
        self.pw_den = wp.zeros(1, dtype=wp.float64, device=d)
        self.pw_nrm = wp.zeros(1, dtype=wp.float64, device=d)
        self.omega_dev = wp.zeros(1, dtype=wp.float64, device=d)
        self.richardson_power_iters = 12

    def build(self, pos: wp.array, regularization: wp.array, finite: wp.array) -> None:
        """Assemble A_c at the live position and factor its block-Jacobi preconditioner (per build)."""
        d = self.device
        self.vals_d.zero_()
        if self.n_xl and self.xl_d is not None:
            wp.launch(assemble_fq_crosslink_kernel, dim=self.n_xl,
                      inputs=[pos, self.xl_d, self.kxl_d, self.r0xl_d, self.node_fiber_d, self.foff_d,
                              self.q_packed_d, self.q_off_d, wp.int32(self.block_rank),
                              wp.int32(self.n_actin), self.rows_d, self.cols_d, self.vals_d], device=d)
        wp.launch(set_fq_diagonal_kernel, dim=self.n_fibers,
                  inputs=[regularization, wp.int32(4 * self.n_xl), self.rows_d, self.cols_d, self.vals_d],
                  device=d)
        self._wsp.bsr_set_from_triplets(self.A_c, self.rows_d, self.cols_d, self.vals_d)
        self._wsp.bsr_get_diag(self.A_c, out=self.diag_blocks)
        wp.launch(factor_block66_kernel, dim=self.n_fibers,
                  inputs=[self.diag_blocks, self.bj_factor, finite], device=d)
        self._estimate_richardson_omega(safety=0.9)

    def _estimate_richardson_omega(self, *, safety: float) -> None:
        r"""Device-resident power iteration -> ``omega = safety / lambda_max(M_bj^-1 A_c)`` (SPD Richardson).

        ``M_bj^-1 A_c`` is self-adjoint in the ``M``-inner-product; the M-Rayleigh quotient
        ``<v, A_c v> / <v, M v>`` (``M`` = the A_c diagonal blocks) converges to ``lambda_max`` from below, so
        ``omega = safety/lambda_max`` (``safety < 1``) lands strictly inside the unconditional-SPD Richardson
        window ``(0, 2/lambda_max)``. All scalars stay on device (no host readback)."""
        d = self.device
        nf = self.n_fibers
        if nf == 0:
            return
        self.pw_v.fill_(_VEC6(1.0))
        for _ in range(self.richardson_power_iters):
            self._wsp.bsr_mv(self.A_c, self.pw_v, self.pw_av)              # av = A_c v
            wp.launch(_fq6_matvec_diag_kernel, dim=nf,
                      inputs=[self.diag_blocks, self.pw_v, self.pw_mv], device=d)  # mv = M v
            self.pw_num.zero_(); self.pw_den.zero_()
            wp.launch(_fq6_dot_kernel, dim=nf, inputs=[self.pw_v, self.pw_av, self.pw_num], device=d)
            wp.launch(_fq6_dot_kernel, dim=nf, inputs=[self.pw_v, self.pw_mv, self.pw_den], device=d)
            wp.launch(apply_block66_kernel, dim=nf,
                      inputs=[self.pw_av, self.bj_factor, self.pw_z], device=d)  # z = M^-1 av
            self.pw_nrm.zero_()
            wp.launch(_fq6_dot_kernel, dim=nf, inputs=[self.pw_z, self.pw_z, self.pw_nrm], device=d)
            wp.launch(_fq6_normalize_kernel, dim=nf, inputs=[self.pw_z, self.pw_nrm, self.pw_v], device=d)
        wp.launch(_fq_set_omega_kernel, dim=1,
                  inputs=[self.pw_num, self.pw_den, wp.float64(float(safety)), self.omega_dev], device=d)

    def apply(self, pos: wp.array, residual: wp.array, preconditioned: wp.array,
              regularization: wp.array, finite: wp.array) -> None:
        r"""Add ``P A_c^-1 P^T residual`` (deep block-Jacobi CG) into ``preconditioned`` (additive, SPD)."""
        d = self.device
        nf = self.n_fibers
        br = wp.int32(self.block_rank)
        wp.launch(fq_restrict_vec6_kernel, dim=nf,
                  inputs=[residual, self.foff_d, self.q_packed_d, self.q_off_d, br, self.coarse_rhs],
                  device=d)
        self.coarse_x.zero_()
        wp.copy(self.coarse_r, self.coarse_rhs)
        wp.launch(apply_block66_kernel, dim=nf,
                  inputs=[self.coarse_r, self.bj_factor, self.coarse_z], device=d)
        wp.copy(self.coarse_p, self.coarse_z)
        self.rz.zero_()
        wp.launch(_fq6_dot_kernel, dim=nf, inputs=[self.coarse_r, self.coarse_z, self.rz], device=d)
        for _ in range(self.iterations):
            self._wsp.bsr_mv(self.A_c, self.coarse_p, self.coarse_ap)
            self.p_ap.zero_()
            wp.launch(_fq6_dot_kernel, dim=nf, inputs=[self.coarse_p, self.coarse_ap, self.p_ap], device=d)
            wp.launch(_fq_alpha_kernel, dim=1, inputs=[self.rz, self.p_ap, self.alpha], device=d)
            wp.launch(_fq6_update_kernel, dim=nf,
                      inputs=[self.coarse_x, self.coarse_r, self.coarse_p, self.coarse_ap, self.alpha],
                      device=d)
            wp.launch(apply_block66_kernel, dim=nf,
                      inputs=[self.coarse_r, self.bj_factor, self.coarse_z], device=d)
            self.rz_new.zero_()
            wp.launch(_fq6_dot_kernel, dim=nf, inputs=[self.coarse_r, self.coarse_z, self.rz_new], device=d)
            wp.launch(_fq_beta_kernel, dim=1, inputs=[self.rz_new, self.rz, self.beta], device=d)
            wp.launch(_fq6_direction_kernel, dim=nf,
                      inputs=[self.coarse_z, self.beta, self.coarse_p], device=d)
        wp.launch(fq_prolong_add_vec6_kernel, dim=nf,
                  inputs=[self.coarse_x, self.foff_d, self.q_packed_d, self.q_off_d, br, preconditioned],
                  device=d)

    def apply_richardson(self, pos: wp.array, residual: wp.array, preconditioned: wp.array,
                         regularization: wp.array, finite: wp.array, *, sweeps: int) -> None:
        r"""Add ``P B_k P^T residual`` (fixed damped block-Jacobi RICHARDSON, LINEAR + SPD) into ``preconditioned``.

        Unlike :meth:`apply` (whose inner **CG** makes the ``residual -> correction`` map a Krylov polynomial
        that depends on ``residual`` — nonlinear, so unusable as a fixed SPD preconditioner block), this runs a
        FIXED ``sweeps`` count of the stationary iteration ``x += omega M_bj^-1 (b_c - A_c x)`` on the assembled
        ``A_c`` with the ``omega = 0.9/lambda_max`` estimated in :meth:`build`. The map ``b_c -> x`` is then the
        fixed polynomial ``B_k = [I - (I - omega M^-1 A_c)^k] A_c^-1`` — LINEAR in the residual and, with the
        symmetric block-Jacobi ``M`` and ``0 < omega < 2/lambda_max``, SPD. This is the coarsest inter-fiber
        solve for the arclength-multigrid V-cycle (which requires a fixed SPD ``M^-1``). ``build`` must have
        assembled ``A_c`` + the block-Jacobi factor + the ``omega`` estimate first.
        """
        d = self.device
        nf = self.n_fibers
        br = wp.int32(self.block_rank)
        wp.launch(fq_restrict_vec6_kernel, dim=nf,
                  inputs=[residual, self.foff_d, self.q_packed_d, self.q_off_d, br, self.coarse_rhs],
                  device=d)
        self.coarse_x.zero_()
        for _ in range(int(sweeps)):
            self._wsp.bsr_mv(self.A_c, self.coarse_x, self.coarse_ap)   # coarse_ap = A_c x
            wp.launch(_fq6_rsub_kernel, dim=nf,
                      inputs=[self.coarse_rhs, self.coarse_ap, self.coarse_r], device=d)  # r = b - A_c x
            wp.launch(apply_block66_kernel, dim=nf,
                      inputs=[self.coarse_r, self.bj_factor, self.coarse_z], device=d)    # z = M_bj^-1 r
            wp.launch(_fq6_axpy_dev_kernel, dim=nf,
                      inputs=[self.coarse_x, self.coarse_z, self.omega_dev], device=d)    # x += omega z
        wp.launch(fq_prolong_add_vec6_kernel, dim=nf,
                  inputs=[self.coarse_x, self.foff_d, self.q_packed_d, self.q_off_d, br, preconditioned],
                  device=d)


# ----------------------------------------------------------------------------------------------------
# G1 Sanity Gate (CPU-only; no CUDA required). Run:
#   ~/miniconda3/envs/ffn_sim/bin/python -m aleph.components.incumbent.fiber_quotient_coarse
# ----------------------------------------------------------------------------------------------------


def _dense_prolongator(pos_np: np.ndarray, foff_np: np.ndarray, data: dict) -> np.ndarray:
    """Reconstruct the dense ``(3N, N_c)`` ``P`` from the packed arrays (CPU gate only)."""
    foff = np.asarray(foff_np, dtype=np.int64)
    n = int(pos_np.shape[0])
    q_packed = data["q_packed"]
    q_off = np.asarray(data["q_off"], dtype=np.int64)
    coarse_off = np.asarray(data["coarse_off"], dtype=np.int64)
    rank = np.asarray(data["rank"], dtype=np.int64)
    block_rank = int(data.get("block_rank", 0))
    P = np.zeros((3 * n, int(data["n_coarse"])), dtype=np.float64)
    for f in range(foff.shape[0] - 1):
        begin = int(foff[f])
        length = int(foff[f + 1]) - begin
        stride = block_rank if block_rank else int(rank[f])  # packed column stride (padded or variable)
        qbase = int(q_off[f])
        cbase = int(coarse_off[f])
        for l in range(length):
            for c in range(stride):
                idx = qbase + (l * stride + c) * 3
                for comp in range(3):
                    P[3 * (begin + l) + comp, cbase + c] = q_packed[idx + comp]
    return P


def _g1_selfcheck(n_fibers: int = 24, nodes_per_fiber: int = 7, seed: int = 20260724) -> bool:
    """G1: straight-fiber rank == 5 and an exact per-fiber rigid motion reproduced by ``P`` to < 1e-9."""
    rng = np.random.default_rng(seed)
    length, m = nodes_per_fiber, n_fibers
    seg = 1.0 / (length - 1)
    centers = rng.uniform(0.0, 12.0, size=(m, 3))
    dirs = rng.normal(size=(m, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    s = (np.arange(length) - (length - 1) / 2.0)[None, :, None] * seg
    pos = (centers[:, None, :] + s * dirs[:, None, :]).reshape(m * length, 3)
    foff = np.arange(0, m * length + 1, length, dtype=np.int32)

    data = build_fiber_rigid_prolongator(pos, foff)
    ranks = np.asarray(data["rank"])
    rank_ok = bool(np.all(ranks == 5))

    P = _dense_prolongator(pos, foff, data)
    # Exact rigid motion of fiber 0: translation + rotation about its centroid.
    f0 = np.arange(foff[0], foff[1])
    rel = pos[f0] - pos[f0].mean(axis=0)
    v = np.zeros((pos.shape[0], 3))
    v[f0] = np.array([0.3, -0.2, 0.1]) + np.cross(np.array([0.05, 0.11, -0.07]), rel)
    v = v.reshape(-1)
    proj = P @ np.linalg.lstsq(P, v, rcond=None)[0]
    rigid_res = float(np.linalg.norm(proj - v) / np.linalg.norm(v))
    rigid_ok = rigid_res < 1e-9

    print("=" * 84)
    print("G1 fiber-quotient prolongator Sanity Gate (CPU, numpy)")
    print("=" * 84)
    print(f"  fibers={m}  nodes/fiber={length}  N_c={data['n_coarse']}  max_rank={data['max_rank']}")
    print(f"  observed unique ranks = {data['ranks_unique']}")
    print(f"  [{'PASS' if rank_ok else 'FAIL'}] every straight fiber has coarse rank EXACTLY 5")
    print(f"  [{'PASS' if rigid_ok else 'FAIL'}] exact fiber rigid motion reproduced by P "
          f"(residual {rigid_res:.3e} < 1e-9)")
    ok = rank_ok and rigid_ok
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    return ok


def _pathb_selfcheck(n_fibers: int = 300, nodes_per_fiber: int = 7, seed: int = 424116) -> bool:
    """Path-B CPU gate (Warp CPU device, no CUDA): assembly correctness (G3) + strong-solve confidence.

    Builds a downscaled CURVED-fiber synthetic (block_rank == 6, like native) with stiff sub-isostatic
    crosslinks (contrast 8e5), assembles A_c with the Path-B warp.sparse kernels, and checks:
      * G3: the assembled A_c matches the numpy Galerkin ``a I + P^T K_xl P`` (the connector-block scatter
        equals the true projected operator).
      * SPD + symmetry.
      * STRONG solve: the deep block-Jacobi CG drives the relative coarse residual ``||A_c x - b|| / ||b||``
        below 1e-6 at a modest budget, where a small budget (the Path-A regime) leaves it large.
    """
    import warp.sparse as wsp

    rng = np.random.default_rng(seed)
    m, length, a, kx = n_fibers, nodes_per_fiber, 1.0, 8.0e5
    seg = 1.0 / (length - 1)
    centers = rng.uniform(0.0, 20.0, size=(m, 3))
    axis = rng.normal(size=(m, 3)); axis /= np.linalg.norm(axis, axis=1, keepdims=True)
    perp = np.cross(axis, rng.normal(size=(m, 3))); perp /= np.linalg.norm(perp, axis=1, keepdims=True)
    s = (np.arange(length) - (length - 1) / 2.0)[None, :, None] * seg
    bow = (0.15 * (seg * (length - 1)) * (1.0 - (2.0 * np.arange(length) / (length - 1) - 1.0) ** 2))
    # curved fibers (axial line + transverse bow) => all 6 rigid modes non-null => block_rank == 6 (native).
    pos = (centers[:, None, :] + s * axis[:, None, :] + bow[None, :, None] * perp[:, None, :]).reshape(m * length, 3)
    n = m * length
    foff = np.arange(0, n + 1, length, dtype=np.int32)
    node_fiber = np.repeat(np.arange(m), length).astype(np.int32)
    ii = rng.integers(0, n, 5 * m); jj = rng.integers(0, n, 5 * m)
    keep = node_fiber[ii] != node_fiber[jj]; ii, jj = ii[keep], jj[keep]; nxl = ii.size
    xl = np.stack([ii, jj], 1).astype(np.int32)
    r0xl = np.linalg.norm(pos[jj] - pos[ii], axis=1) * 0.9
    kxl = np.full(nxl, kx)

    data = build_fiber_rigid_prolongator(pos, foff, pad_to=6)
    P = _dense_prolongator(pos, foff, data); n_c = int(data["n_coarse"])
    K = np.zeros((3 * n, 3 * n))
    for t in range(nxl):
        dl = pos[jj[t]] - pos[ii[t]]; ln = np.linalg.norm(dl); u = dl / ln
        tau = max(kxl[t] * (ln - r0xl[t]) / ln, 0.0)
        tng = tau * np.eye(3) + (kxl[t] - tau) * np.outer(u, u)
        p, q = int(ii[t]), int(jj[t])
        for (rr, cc, sg) in [(p, p, 1), (q, q, 1), (p, q, -1), (q, p, -1)]:
            K[3 * rr:3 * rr + 3, 3 * cc:3 * cc + 3] += sg * tng
    Ac_ref = a * np.eye(n_c) + P.T @ K @ P

    class _Cell:
        device = "cpu"; n_total = n; n_actin = n; n_fibers = m
        pos_d = wp.array(pos, dtype=wp.vec3d, device="cpu")
        foff_d = wp.array(foff, dtype=wp.int32, device="cpu")
        node_fiber_d = wp.array(node_fiber, dtype=wp.int32, device="cpu")
        xl_d = wp.array(xl, dtype=wp.int32, device="cpu")
        kxl_d = wp.array(kxl, dtype=wp.float64, device="cpu")
        r0xl_d = wp.array(r0xl, dtype=wp.float64, device="cpu")
        n_xl = nxl

    fq = FiberQuotientCoarsePathB(_Cell(), inner_iterations=40)
    reg = wp.array([a], dtype=wp.float64, device="cpu"); finite = wp.ones(1, dtype=wp.int32, device="cpu")
    fq.build(_Cell.pos_d, reg, finite)

    xr = rng.standard_normal(n_c)
    xd = wp.array(xr.reshape(m, 6), dtype=_VEC6, device="cpu"); yd = wp.zeros(m, dtype=_VEC6, device="cpu")
    wsp.bsr_mv(fq.A_c, xd, yd)
    g3 = float(np.linalg.norm(yd.numpy().reshape(-1) - Ac_ref @ xr) / np.linalg.norm(Ac_ref @ xr))
    eig = float(np.linalg.eigvalsh(Ac_ref).min())
    sym = float(np.abs(Ac_ref - Ac_ref.T).max() / np.abs(Ac_ref).max())

    r3 = rng.standard_normal((n, 3)); r3d = wp.array(r3, dtype=wp.vec3d, device="cpu")
    z = wp.zeros(n, dtype=wp.vec3d, device="cpu")
    rhs = P.T @ r3.reshape(-1)
    residuals = {}
    for budget in (40, 600):
        fq.iterations = budget
        z.zero_(); fq.apply(_Cell.pos_d, r3d, z, reg, finite)
        cx = fq.coarse_x.numpy().reshape(-1)
        residuals[budget] = float(np.linalg.norm(Ac_ref @ cx - rhs) / np.linalg.norm(rhs))

    g3_ok = g3 < 1e-10
    spd_ok = eig > 0.0 and sym < 1e-10
    strong_ok = residuals[600] < 1e-6 and residuals[40] > residuals[600]
    print("=" * 84)
    print("Path-B fiber-quotient A_c Sanity Gate (Warp CPU device, no CUDA)")
    print("=" * 84)
    print(f"  fibers={m}  N_c={n_c}  block_rank={data['block_rank']}  ranks={data['ranks_unique']}  "
          f"crosslinks={nxl}  contrast k/a={kx / a:.0e}")
    print(f"  [{'PASS' if g3_ok else 'FAIL'}] G3 assembled A_c == a I + P^T K_xl P  (rel err {g3:.2e})")
    print(f"  [{'PASS' if spd_ok else 'FAIL'}] A_c SPD (min eig {eig:.3e} > 0), symmetric ({sym:.1e})")
    print(f"  [{'PASS' if strong_ok else 'FAIL'}] deep block-Jacobi CG STRONG: "
          f"coarse residual budget40 {residuals[40]:.2e} -> budget600 {residuals[600]:.2e} (< 1e-6)")
    ok = g3_ok and spd_ok and strong_ok
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    import sys

    passed = _g1_selfcheck()
    print()
    passed = _pathb_selfcheck() and passed
    sys.exit(0 if passed else 1)
