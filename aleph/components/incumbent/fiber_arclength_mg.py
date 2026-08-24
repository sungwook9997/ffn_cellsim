r"""Matrix-free per-fiber arclength geometric-multigrid building blocks (2:1 R/P + O(L) line smoother).

These are the first two matrix-free pieces of the fiber-arclength geometric multigrid recommended in
``aleph/docs/v2_audit/FINE_MESH_MULTIGRID_DESIGN_2026-07-24.md`` (§2.2 the hierarchy, §2.4 the line
smoother, §4.1 pieces 2+3). They are **additive and default-off**: nothing here is wired into
``ProjectedAnalyticCG`` yet — the Lead wires the V-cycle (§2.3) and native-validates on the §4.4 ladder. With
the multigrid switched off the existing coarse-mesh solver path is bit-identical (this module is not imported
by the incumbent solve).

What is delivered here (the two memory + speed wins that retire the ``O(L^2)`` per-fiber Cholesky block):

1. **2:1 per-fiber arclength restriction / prolongation** (§2.2). Each fiber's nodes ``foff[f]..foff[f+1]`` are
   a contiguous 1-D chain (``cell.foff_d``). :func:`build_arclength_restriction` builds, once at setup, a 2:1
   arclength semicoarsening that keeps even-indexed fine nodes as coarse nodes (injection) and linearly
   interpolates the dropped odd nodes from their two kept neighbours; a fiber with ``L_f`` fine nodes yields
   ``ceil(L_f/2)`` coarse nodes. The device kernels :func:`arclength_prolong_add_kernel` (``P``, coarse->fine
   linear interpolation) and :func:`arclength_restrict_kernel` (``R = P^T``, fine->coarse full-weighting) apply
   it matrix-free, one thread per fine node, using the SAME per-fine-node ``(coarse_index, weight)`` table so
   ``R`` is EXACTLY ``P^T`` (the Galerkin requirement ``A_c = R A P`` with a symmetric ``A_c``). The stencil is
   purely index-based on the arclength ordering (geometry-free), so it applies unchanged to the heterogeneous
   mixed formin+Arp2/3 cortex — every fiber coarsens to its own floor and terminates at the shared
   fiber-quotient level (§2.2). Analogue of ``fiber_quotient_coarse.fiber_rigid_restrict_kernel`` /
   ``fiber_rigid_prolong_add_kernel`` but with the 2-weight linear stencil instead of the rigid ``Q_f`` block.

2. **O(L) line (arclength) smoother** (§2.4). :func:`fiber_line_smooth_kernel` runs one damped block-Jacobi
   line relaxation ``x += omega * M_f^-1 (b - A x)`` per fiber, where ``M_f`` is a per-fiber SPD **tridiagonal**
   surrogate of the fiber's own along-arclength sub-operator ``D_external + D2^T diag(alpha) D2`` — the exact
   ``O(L^2)``-dense operator that ``build_fiber_block_cholesky_kernel`` factors, read from the SAME live
   ``coarse_external_diagonal`` (per-node vec3d) and ``alpha_d``/``toff_d`` (per-triple bending). The NF2007
   discrete bending ``D2^T diag(alpha) D2`` is a 4th-difference (pentadiagonal) operator; ``M_f`` keeps its
   tridiagonal band and **row-sum-lumps** the two bandwidth-2 entries onto the diagonal. That lumping is exact
   on both spectral ends — the surrogate annihilates the constant/translation near-null (zero row sums, so the
   smooth mode is preserved, ``no-op`` on it) and reproduces the true oscillatory eigenvalue ``16 alpha`` at
   ``theta = pi`` (so the high-arclength-frequency error is damped) — while a Thomas solve makes one sweep
   ``O(L)`` in time and storage, retiring the ``O(L^2)``-storage / ``O(L^3)``-factor dense Cholesky (the
   §9/§3 2.85 GB block). A symmetric direct tridiagonal solve keeps the smoother SPD, so a symmetric V-cycle
   stays a valid SPD preconditioner (§2.3).

How this plugs into a V-cycle next (the Lead wires ``ProjectedAnalyticCG``; §2.3, §4.3): each level ``l`` holds
its own ``build_arclength_restriction`` table (call it recursively on the coarsened ``foff``). ``Vcycle(l, b)``
pre-smooths with :func:`fiber_line_smooth_kernel` (the residual ``b - A x`` reuses the exact matrix-free
``ProjectedAnalyticCG._operator``), restricts the residual with :func:`arclength_restrict_kernel`, recurses,
prolong-adds the coarse correction with :func:`arclength_prolong_add_kernel`, and post-smooths. The Galerkin
coarse operator is applied matrix-free as ``R (A (P .))`` reusing ``_operator`` (never assembled at the fine
level). The recursion terminates at the fiber-quotient graph (``FiberQuotientCoarsePathB``, this package's
``fiber_quotient_coarse``) as the coarsest solver — where each fiber is ~1 point and the residual is
inter-fiber near-rigid. Used as the SPD preconditioner ``M^-1`` inside the existing outer PCG, this is
preconditioning only: it cannot move the fixed point or any residual gate (CLAUDE.md no-gate-loosening).

Sanity Gate:
    * Dimensions: R/P weights are dimensionless arclength interpolation coefficients; the smoother diagonal is
      ``[pN/um]`` (external diagonal + bending), the correction ``[um] = omega * M^-1 [pN] * ... `` — ``omega``
      is dimensionless, ``M^-1`` has units ``[um/pN]``, residual ``[pN]``.
    * Transpose consistency (Galerkin): ``<R x, y>_c == <x, P y>_f`` to round-off, by sharing one weight table.
    * Linear-interpolation exactness: ``P`` reproduces any arclength-linear coarse field exactly on interior
      odd nodes (the near-null / smooth-bending capture that makes multigrid h-independent, §2.5).
    * Near-null: the row-sum-lumped tridiagonal bending has zero row sums => it annihilates the constant mode
      exactly, so the line smoother is a no-op on the translation near-null and frequency-selective.
    * SPD / boundary cases: a straight/degenerate fiber's ``M_f`` diagonal is ``external_diagonal >= a > 0``; a
      non-positive Thomas pivot can only be non-finite input and latches the device finite predicate (never a
      spurious correction). ``L_f in {1, 2}`` fibers (no triples) reduce to a diagonal solve and do not crash.
    * No tuned constant: ``omega`` is the standard damped-Jacobi / line-relaxation smoothing weight (a solver
      compute parameter, grid-invariant, not chosen to pass a biological gate); the 2:1 coarsening ratio and
      the linear stencil are geometric, not physics knobs.
"""

from __future__ import annotations

import logging
from typing import Callable

import numpy as np
import warp as wp

_LOG = logging.getLogger(__name__)

__all__ = [
    "build_arclength_restriction",
    "arclength_restrict_kernel",
    "arclength_prolong_add_kernel",
    "fiber_line_smooth_kernel",
    "fiber_perfiber_operator_kernel",
    "restrict_diagonal_sq_kernel",
    "FiberArclengthMultigrid",
]


# ----------------------------------------------------------------------------------------------------
# Host: per-fiber 2:1 arclength restriction / prolongation table (built once at setup, NOT in hot loop).
# ----------------------------------------------------------------------------------------------------


def build_arclength_restriction(foff_np: np.ndarray) -> dict:
    r"""Build one 2:1 arclength coarsening level (fine -> coarse) for every fiber.

    Each fiber ``f`` with fine nodes ``foff[f]..foff[f+1]`` (local index ``k = 0..L_f-1``) coarsens 2:1 keeping
    the even-indexed fine nodes as coarse nodes (coarse local ``j`` sits at fine local ``2 j``), giving
    ``L_c = ceil(L_f/2)`` coarse nodes. The prolongation ``P`` (coarse -> fine) is 1-D linear interpolation
    along arclength:

      * even fine node ``k = 2 j``      -> weight 1 from coarse ``j`` (injection / endpoint-preserving);
      * odd  fine node ``k = 2 j + 1``  -> weight 1/2 from coarse ``j`` and 1/2 from coarse ``j + 1``
        (linear interpolation); at the trailing boundary of an even-``L_f`` fiber the last odd node has no
        right coarse neighbour and takes weight 1 from its left coarse node (constant/endpoint extrapolation).

    The restriction is ``R = P^T`` (full-weighting), built from the SAME per-fine-node ``(coarse_index,
    weight)`` slots so ``R`` is EXACTLY the transpose of ``P`` (Galerkin ``A_c = R A P`` is then symmetric). The
    table is index-based on the arclength ordering only (no geometry), so it applies unchanged to heterogeneous
    fibers of any length (formin 41-node, Arp2/3 4-node, ...).

    Args:
        foff_np: ``(n_fibers + 1,)`` CSR fiber -> fine-node offsets (``cell.foff_d``). Fine nodes are assumed
            contiguous per fiber in arclength order (the invariant everywhere in ``ac/cell``).

    Returns:
        dict with keys (all numpy, upload with :func:`warp.array`):
          ``foff_coarse`` (int32, ``n_fibers + 1``): CSR fiber -> coarse-node offsets (the coarse ``foff`` for
            recursing to the next level).
          ``n_coarse`` (int): total coarse nodes ``sum_f ceil(L_f/2)``.
          ``c0`` (int32, ``n_fine``): first coarse-node index each fine node draws from (``-1`` never occurs;
            every fine node has >= 1 weight).
          ``w0`` (float64, ``n_fine``): weight for ``c0``.
          ``c1`` (int32, ``n_fine``): second coarse-node index (``-1`` if the fine node is an injection or a
            boundary node with a single weight).
          ``w1`` (float64, ``n_fine``): weight for ``c1`` (0 where ``c1 == -1``).
    """
    foff = np.ascontiguousarray(foff_np, dtype=np.int64)
    n_fibers = int(foff.shape[0] - 1)
    n_fine = int(foff[-1]) if n_fibers else 0

    coarse_off = np.zeros(n_fibers + 1, dtype=np.int64)
    c0 = np.full(n_fine, -1, dtype=np.int32)
    c1 = np.full(n_fine, -1, dtype=np.int32)
    w0 = np.zeros(n_fine, dtype=np.float64)
    w1 = np.zeros(n_fine, dtype=np.float64)

    for f in range(n_fibers):
        begin = int(foff[f])
        length = int(foff[f + 1]) - begin
        length_c = (length + 1) // 2  # ceil(L_f / 2)
        coarse_off[f + 1] = coarse_off[f] + length_c
        cbase = int(coarse_off[f])
        for k in range(length):
            gi = begin + k
            if k % 2 == 0:  # even -> injection from coarse j = k/2
                c0[gi] = cbase + (k // 2)
                w0[gi] = 1.0
            else:  # odd -> linear interpolation from j and j+1
                jl = (k - 1) // 2
                jr = jl + 1
                if jr < length_c:
                    c0[gi] = cbase + jl
                    w0[gi] = 0.5
                    c1[gi] = cbase + jr
                    w1[gi] = 0.5
                else:  # trailing boundary node of an even-L fiber: no right coarse neighbour
                    c0[gi] = cbase + jl
                    w0[gi] = 1.0

    return {
        "foff_coarse": np.ascontiguousarray(coarse_off, dtype=np.int32),
        "n_coarse": int(coarse_off[-1]),
        "c0": c0,
        "w0": w0,
        "c1": c1,
        "w1": w1,
    }


# ----------------------------------------------------------------------------------------------------
# Device kernels: 2:1 arclength restrict / prolong (one thread per FINE node).
# ----------------------------------------------------------------------------------------------------


@wp.kernel
def arclength_restrict_kernel(
    fine: wp.array(dtype=wp.vec3d),
    c0: wp.array(dtype=wp.int32),
    w0: wp.array(dtype=wp.float64),
    c1: wp.array(dtype=wp.int32),
    w1: wp.array(dtype=wp.float64),
    coarse: wp.array(dtype=wp.vec3d),
) -> None:
    r"""``coarse = R fine`` with ``R = P^T`` (full-weighting). Caller zeroes ``coarse`` first; one thread/fine node."""
    i = wp.tid()
    value = fine[i]
    j0 = c0[i]
    if j0 >= 0:
        wp.atomic_add(coarse, j0, w0[i] * value)
    j1 = c1[i]
    if j1 >= 0:
        wp.atomic_add(coarse, j1, w1[i] * value)


@wp.kernel
def arclength_prolong_add_kernel(
    coarse: wp.array(dtype=wp.vec3d),
    c0: wp.array(dtype=wp.int32),
    w0: wp.array(dtype=wp.float64),
    c1: wp.array(dtype=wp.int32),
    w1: wp.array(dtype=wp.float64),
    fine: wp.array(dtype=wp.vec3d),
) -> None:
    r"""``fine += P coarse`` (1-D linear arclength interpolation). One thread per fine node."""
    i = wp.tid()
    acc = wp.vec3d(0.0, 0.0, 0.0)
    j0 = c0[i]
    if j0 >= 0:
        acc = acc + w0[i] * coarse[j0]
    j1 = c1[i]
    if j1 >= 0:
        acc = acc + w1[i] * coarse[j1]
    fine[i] = fine[i] + acc


# ----------------------------------------------------------------------------------------------------
# Device kernel: O(L) per-fiber line (arclength) smoother — damped block-Jacobi tridiagonal Thomas relaxation.
# ----------------------------------------------------------------------------------------------------


@wp.kernel
def fiber_line_smooth_kernel(
    residual: wp.array(dtype=wp.vec3d),
    x: wp.array(dtype=wp.vec3d),
    foff: wp.array(dtype=wp.int32),
    toff: wp.array(dtype=wp.int32),
    alpha: wp.array(dtype=wp.float64),
    external_diagonal: wp.array(dtype=wp.vec3d),
    omega: wp.float64,
    diag_bend: wp.array(dtype=wp.float64),
    off_sub: wp.array(dtype=wp.float64),
    cprime: wp.array(dtype=wp.vec3d),
    dprime: wp.array(dtype=wp.vec3d),
    finite: wp.array(dtype=wp.int32),
) -> None:
    r"""One damped line relaxation ``x += omega * M_f^-1 residual`` per fiber (``M_f`` SPD tridiagonal).

    ``M_f`` is the fiber's own along-arclength sub-operator ``D_external + D2^T diag(alpha) D2`` reduced to its
    tridiagonal band by **row-sum lumping** the two bandwidth-2 NF2007 bending entries onto the diagonal (per
    local triple ``t`` on nodes ``t, t+1, t+2`` with stencil ``[1, -2, 1]``: the full pentadiagonal contributes
    ``+alpha`` at the ``(t, t+2)`` corners; lumping moves each onto its diagonal, giving the per-triple
    tridiagonal stencil diag ``+2a, +4a, +2a`` and sub/super ``-2a, -2a``). Row sums stay zero, so the constant
    near-null is annihilated (no-op on the smooth mode) and the ``theta = pi`` eigenvalue is the exact
    ``16 alpha`` (high-frequency damping). The three Cartesian components share the scalar bending band but keep
    their own ``external_diagonal`` component, solved in lockstep as one vec3d Thomas sweep (``O(L)`` per
    fiber). Disjoint per-fiber node ranges => race-free. One thread per fiber.

    ``diag_bend``, ``off_sub`` (float64, size >= n_total) and ``cprime``, ``dprime`` (vec3d, size >= n_total)
    are per-fiber scratch (indexed on the fiber's own disjoint node slice). ``finite`` latches to 0 on a
    non-positive Thomas pivot (non-finite input) so the sweep contributes no spurious correction.
    """
    f = wp.tid()
    begin = foff[f]
    node_count = foff[f + 1] - begin
    if node_count <= 0:
        return
    triple_begin = toff[f]

    # ---- Assemble the row-sum-lumped tridiagonal bending band (component-independent scalar). ----
    for k in range(node_count):
        diag_bend[begin + k] = wp.float64(0.0)
        off_sub[begin + k] = wp.float64(0.0)
    for t in range(node_count - 2):
        a = alpha[triple_begin + t]
        two_a = wp.float64(2.0) * a
        four_a = wp.float64(4.0) * a
        diag_bend[begin + t] = diag_bend[begin + t] + two_a
        diag_bend[begin + t + 1] = diag_bend[begin + t + 1] + four_a
        diag_bend[begin + t + 2] = diag_bend[begin + t + 2] + two_a
        # sub-diagonal entry [row, row-1]; symmetric so super[row-1, row] == sub[row, row-1].
        off_sub[begin + t + 1] = off_sub[begin + t + 1] - two_a
        off_sub[begin + t + 2] = off_sub[begin + t + 2] - two_a

    # ---- vec3d Thomas solve of the SPD tridiagonal M_f (3 componentwise systems in lockstep). ----
    # Diagonal b[k] = external_diagonal[k] + diag_bend[k] (vec3d); sub[k] = off_sub[k] (scalar, [k, k-1]);
    # super[k] = off_sub[k+1] (scalar, [k, k+1], symmetric). RHS = residual.
    ext0 = external_diagonal[begin]
    db0 = diag_bend[begin]
    b0x = ext0[0] + db0
    b0y = ext0[1] + db0
    b0z = ext0[2] + db0
    if b0x <= wp.float64(1.0e-300) or b0y <= wp.float64(1.0e-300) or b0z <= wp.float64(1.0e-300):
        finite[0] = wp.int32(0)
        return
    if not (wp.isfinite(b0x) and wp.isfinite(b0y) and wp.isfinite(b0z)):
        finite[0] = wp.int32(0)
        return
    sup0 = wp.float64(0.0)
    if node_count > 1:
        sup0 = off_sub[begin + 1]
    cprime[begin] = wp.vec3d(sup0 / b0x, sup0 / b0y, sup0 / b0z)
    r0 = residual[begin]
    dprime[begin] = wp.vec3d(r0[0] / b0x, r0[1] / b0y, r0[2] / b0z)

    for k in range(1, node_count):
        node = begin + k
        ext = external_diagonal[node]
        db = diag_bend[node]
        sub = off_sub[node]  # entry [k, k-1]
        cp_prev = cprime[node - 1]
        dp_prev = dprime[node - 1]
        mx = ext[0] + db - sub * cp_prev[0]
        my = ext[1] + db - sub * cp_prev[1]
        mz = ext[2] + db - sub * cp_prev[2]
        if mx <= wp.float64(1.0e-300) or my <= wp.float64(1.0e-300) or mz <= wp.float64(1.0e-300):
            finite[0] = wp.int32(0)
            return
        if not (wp.isfinite(mx) and wp.isfinite(my) and wp.isfinite(mz)):
            finite[0] = wp.int32(0)
            return
        sup = wp.float64(0.0)
        if k < node_count - 1:
            sup = off_sub[node + 1]  # entry [k, k+1] == off_sub[k+1]
        cprime[node] = wp.vec3d(sup / mx, sup / my, sup / mz)
        rr = residual[node]
        dprime[node] = wp.vec3d(
            (rr[0] - sub * dp_prev[0]) / mx,
            (rr[1] - sub * dp_prev[1]) / my,
            (rr[2] - sub * dp_prev[2]) / mz,
        )

    # ---- Back substitution (dprime becomes the solved correction e) + damped update x += omega * e. ----
    last = begin + node_count - 1
    e_next = dprime[last]
    xv = x[last]
    x[last] = wp.vec3d(xv[0] + omega * e_next[0], xv[1] + omega * e_next[1], xv[2] + omega * e_next[2])
    for reverse in range(node_count - 1):
        k = node_count - 2 - reverse
        node = begin + k
        cp = cprime[node]
        dp = dprime[node]
        e = wp.vec3d(dp[0] - cp[0] * e_next[0], dp[1] - cp[1] * e_next[1], dp[2] - cp[2] * e_next[2])
        xv2 = x[node]
        x[node] = wp.vec3d(xv2[0] + omega * e[0], xv2[1] + omega * e[1], xv2[2] + omega * e[2])
        e_next = e


# ----------------------------------------------------------------------------------------------------
# Device kernel: assembled per-fiber coarse operator (ext diagonal + full NF2007 pentadiagonal bending).
# ----------------------------------------------------------------------------------------------------


@wp.kernel
def fiber_perfiber_operator_kernel(
    x: wp.array(dtype=wp.vec3d),
    foff: wp.array(dtype=wp.int32),
    toff: wp.array(dtype=wp.int32),
    alpha: wp.array(dtype=wp.float64),
    ext: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.vec3d),
) -> None:
    r"""``out = A_perfiber x`` per fiber: ``ext ⊙ x + (D2^T diag(alpha) D2) x`` (matrix-free, ``O(L)``).

    This is the assembled **per-fiber block** of the fine operator ``aI + P K P`` at a coarse arclength level:
    the node-diagonal ``ext`` (the live ``coarse_external_diagonal`` restricted with the squared interpolation
    weights — ``regularization a`` + crosslink self-coupling + external turgor/ERM/LINC) plus the FULL NF2007
    pentadiagonal bending ``sum_t alpha_t D2_t^T D2_t`` (per local triple ``t`` on nodes ``t, t+1, t+2`` with the
    ``[1, -2, 1]`` second-difference stencil). It drops ONLY the inter-fiber crosslink OFF-diagonal (captured
    ADDITIVELY by the fiber-quotient coarsest level, never nested here) and the inextensibility projector ``P``
    (which the sibling line smoother also omits — so this operator is consistent with the smoother it drives).

    Using it as the coarse-level operator (levels ``l >= 1``) replaces the full-fine ``R ∘ _operator ∘ P`` apply
    (a 2.9M-node ``_operator`` at EVERY level) with a cheap ``O(n_nodes[l])`` per-fiber apply — the whole point
    of a multigrid coarse level. SPD: ``ext > 0`` (>= regularization ``a`` > 0) and ``D2^T diag(alpha) D2`` is
    PSD (``alpha >= 0``), so ``A_perfiber`` is SPD and the coarse Richardson/line-smoother solve stays a valid
    SPD preconditioner component. One thread per fiber (disjoint node ranges => race-free).
    """
    f = wp.tid()
    begin = foff[f]
    node_count = foff[f + 1] - begin
    if node_count <= 0:
        return
    for k in range(node_count):
        node = begin + k
        e = ext[node]
        xv = x[node]
        out[node] = wp.vec3d(e[0] * xv[0], e[1] * xv[1], e[2] * xv[2])
    triple_begin = toff[f]
    for t in range(node_count - 2):
        a = alpha[triple_begin + t]
        n0 = begin + t
        n1 = begin + t + 1
        n2 = begin + t + 2
        x0 = x[n0]
        x1 = x[n1]
        x2 = x[n2]
        d2x = x0[0] - wp.float64(2.0) * x1[0] + x2[0]
        d2y = x0[1] - wp.float64(2.0) * x1[1] + x2[1]
        d2z = x0[2] - wp.float64(2.0) * x1[2] + x2[2]
        ad2 = wp.vec3d(a * d2x, a * d2y, a * d2z)
        out[n0] = out[n0] + ad2
        o1 = out[n1]
        out[n1] = wp.vec3d(o1[0] - wp.float64(2.0) * ad2[0], o1[1] - wp.float64(2.0) * ad2[1],
                           o1[2] - wp.float64(2.0) * ad2[2])
        out[n2] = out[n2] + ad2


# ----------------------------------------------------------------------------------------------------
# Device kernels: V-cycle glue (Galerkin-diagonal restriction of the external diagonal + vec3d saxpy/copy).
# ----------------------------------------------------------------------------------------------------


@wp.kernel
def restrict_diagonal_sq_kernel(
    fine_diag: wp.array(dtype=wp.vec3d),
    c0: wp.array(dtype=wp.int32),
    w0: wp.array(dtype=wp.float64),
    c1: wp.array(dtype=wp.int32),
    w1: wp.array(dtype=wp.float64),
    coarse_diag: wp.array(dtype=wp.vec3d),
) -> None:
    r"""Galerkin diagonal restriction ``coarse_diag[c] = sum_i P_ic^2 fine_diag[i]`` (caller zeroes coarse).

    ``fine_diag`` is a per-node componentwise diagonal (the smoother's external diagonal ``D_ext`` at the
    current level). The exact diagonal of the Galerkin coarse of a diagonal operator is ``P^T diag(d) P``'s
    diagonal, i.e. the SQUARED interpolation weights times the fine diagonal — positive whenever ``fine_diag``
    is, so the coarse smoother stays SPD. One thread per fine node (same stencil as ``arclength_restrict``).
    """
    i = wp.tid()
    d = fine_diag[i]
    j0 = c0[i]
    if j0 >= 0:
        s0 = w0[i] * w0[i]
        wp.atomic_add(coarse_diag, j0, wp.vec3d(s0 * d[0], s0 * d[1], s0 * d[2]))
    j1 = c1[i]
    if j1 >= 0:
        s1 = w1[i] * w1[i]
        wp.atomic_add(coarse_diag, j1, wp.vec3d(s1 * d[0], s1 * d[1], s1 * d[2]))


@wp.kernel
def _mg_sub_kernel(
    a: wp.array(dtype=wp.vec3d), b: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3d),
) -> None:
    """``out = a - b`` (the level residual ``b_l - A_l x_l``)."""
    i = wp.tid()
    out[i] = a[i] - b[i]


@wp.kernel
def _mg_add_kernel(
    a: wp.array(dtype=wp.vec3d), b: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3d),
) -> None:
    """``out = a + b`` (the additive coarsest correction: line-smoother solve + fiber-quotient rigid)."""
    i = wp.tid()
    out[i] = a[i] + b[i]


@wp.kernel
def _mg_copy_kernel(src: wp.array(dtype=wp.vec3d), dst: wp.array(dtype=wp.vec3d)) -> None:
    """``dst[i] = src[i]`` over the launch range (used to move the leading actin block between buffers)."""
    i = wp.tid()
    dst[i] = src[i]


@wp.kernel
def _mg_scale_kernel(vec: wp.array(dtype=wp.vec3d), scale: wp.array(dtype=wp.float64)) -> None:
    """``vec *= scale[0]`` — pre-scale the smoother residual by the device damping ``omega`` (M^-1 is linear, so
    ``x += M^-1 (omega r) == x += omega M^-1 r``), keeping the fiber line smoother's constant-omega signature."""
    i = wp.tid()
    vec[i] = scale[0] * vec[i]


@wp.kernel
def _mg_dot_kernel(a: wp.array(dtype=wp.vec3d), b: wp.array(dtype=wp.vec3d),
                   out: wp.array(dtype=wp.float64)) -> None:
    """``out[0] += <a, b>`` (device reduction for the smoother-omega power iteration)."""
    i = wp.tid()
    wp.atomic_add(out, 0, wp.dot(a[i], b[i]))


@wp.kernel
def _mg_normalize_kernel(z: wp.array(dtype=wp.vec3d), nrm2: wp.array(dtype=wp.float64),
                         out: wp.array(dtype=wp.vec3d)) -> None:
    """``out = z / ||z||_2`` (guarded) — power-iteration renormalization."""
    i = wp.tid()
    n = nrm2[0]
    if n > wp.float64(1.0e-300):
        out[i] = z[i] / wp.sqrt(n)
    else:
        out[i] = z[i]


@wp.kernel
def _mg_set_omega_kernel(num: wp.array(dtype=wp.float64), den: wp.array(dtype=wp.float64),
                         safety: wp.float64, omega_cap: wp.float64,
                         omega: wp.array(dtype=wp.float64)) -> None:
    r"""``omega = min(omega_cap, safety / lambda_max)`` with the A-Rayleigh estimate ``lambda_max = num/den``.

    ``M^-1 A`` is self-adjoint in the A-inner-product, so ``num/den = <Av, M^-1 Av>/<v, Av>`` converges to
    ``lambda_max(M^-1 A)`` from below; ``safety < 1`` keeps ``omega < 1/lambda_max < 2/lambda_max`` (the
    unconditional-SPD damped-block-Jacobi range). The ``omega_cap`` prevents an over-large step when
    ``lambda_max`` is small (a well-conditioned level)."""
    n = num[0]
    d = den[0]
    if n > wp.float64(1.0e-300) and d > wp.float64(1.0e-300) and wp.isfinite(n) and wp.isfinite(d):
        omega[0] = wp.min(omega_cap, safety * d / n)
    else:
        omega[0] = omega_cap


# ----------------------------------------------------------------------------------------------------
# Host helpers: coarsen the per-fiber bending band (toff / alpha) for a 2:1 arclength level.
# ----------------------------------------------------------------------------------------------------


def _coarsen_bending(
    foff_prev: np.ndarray, toff_prev: np.ndarray, alpha_prev: np.ndarray, foff_coarse: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Build the coarse per-fiber triple offsets + bending coefficients for one 2:1 arclength level.

    A coarse fiber of ``L_c`` nodes has ``max(L_c - 2, 0)`` bending triples. Linear-interpolation Galerkin of
    the raw NF2007 bending sub-block ``alpha * (D2^T D2)`` on a 2:1 grid scales the coefficient by ``1/4`` (the
    fine second difference of a piecewise-linear coarse field is half the coarse second difference at the kept
    nodes and zero at the dropped nodes, so ``P^T (D2^T alpha D2) P = (alpha/4) D2c^T D2c``). Heterogeneous
    ``alpha`` is sampled at the fine triple centred on the coarse triple's middle node (``2 j_c + 1``), clamped
    to the fiber's own triple range — a surrogate for the smoother, exact for the uniform-``alpha`` case.
    """
    foff_prev = np.asarray(foff_prev, dtype=np.int64)
    foff_coarse = np.asarray(foff_coarse, dtype=np.int64)
    toff_prev = np.asarray(toff_prev, dtype=np.int64)
    alpha_prev = np.asarray(alpha_prev, dtype=np.float64)
    n_fibers = int(foff_prev.shape[0] - 1)
    toff_c = np.zeros(n_fibers + 1, dtype=np.int64)
    alpha_c_list: list[float] = []
    for f in range(n_fibers):
        length_p = int(foff_prev[f + 1] - foff_prev[f])
        length_c = int(foff_coarse[f + 1] - foff_coarse[f])
        n_tri_c = max(length_c - 2, 0)
        n_tri_p = max(length_p - 2, 0)
        tp_base = int(toff_prev[f])
        for jc in range(n_tri_c):
            if n_tri_p > 0:
                tp = min(max(2 * jc + 1, 0), n_tri_p - 1)
                alpha_c_list.append(0.25 * float(alpha_prev[tp_base + tp]))
            else:
                alpha_c_list.append(0.0)
        toff_c[f + 1] = toff_c[f] + n_tri_c
    alpha_c = np.asarray(alpha_c_list, dtype=np.float64) if alpha_c_list else np.zeros(0, dtype=np.float64)
    return np.ascontiguousarray(toff_c, dtype=np.int32), np.ascontiguousarray(alpha_c, dtype=np.float64)


# ====================================================================================================
# Driver: matrix-free fiber-arclength geometric-multigrid V-cycle (SPD preconditioner for ProjectedAnalyticCG).
# ====================================================================================================


class FiberArclengthMultigrid:
    r"""Matrix-free fiber-arclength geometric-multigrid V-cycle (design memo §2.2-2.4, §4.3).

    Builds a per-fiber 2:1 arclength semicoarsening hierarchy (each level's ``R = P^T`` from
    :func:`build_arclength_restriction`) down to short fibers, then attaches
    :class:`FiberQuotientCoarsePathB` as the coarsest inter-fiber level. :meth:`apply` runs a **symmetric**
    additively-corrected V(``nu1``, ``nu2``) cycle used as the SPD preconditioner ``M^-1`` inside the existing
    outer PCG:

      * pre-smooth with the ``O(L)`` :func:`fiber_line_smooth_kernel` (damped block-Jacobi line relaxation),
      * matrix-free residual via the exact fine ``_operator`` (the ``aI + P K P`` apply),
      * restrict (``R``) -> recurse -> prolong-add (``P``) -> post-smooth (mirrored ``nu2 == nu1`` sweeps),

    with every coarse operator applied MATRIX-FREE as ``R ∘ _operator ∘ P`` (never assembled at the fine
    level). Bramble's symmetric construction (``z += omega M^-1 (r - A z)`` pre, symmetric coarse correction,
    equal post-smooth) keeps ``M^-1`` SPD, so the outer PCG stays valid — preconditioning only, it cannot move
    the fixed point or any residual gate (CLAUDE.md no-gate-loosening). The ``finite`` latch is carried through
    every operator apply (a non-finite pivot in the line smoother / coarsest Cholesky latches it).

    It preconditions only the actin fiber block (the leading ``n_actin`` nodes owned by ``foff``); the caller
    node-Jacobi-preconditions the disjoint remaining DOFs, so the composite ``M^-1`` is block-diagonal SPD.

    Sanity Gate:
        * SPD/symmetry: ``M^-1`` built column-by-column is symmetric to round-off and has positive eigenvalues
          (symmetric V-cycle + SPD smoother + ``R == P^T`` + SPD coarsest solve).
        * h-independence: one V-cycle's residual-energy reduction factor is roughly constant across fiber
          lengths (the multigrid signature).
        * Off-by-construction: nothing here is imported by the incumbent solve; the class is instantiated only
          under the ``multigrid`` flag / ``AC_MG`` env, and the OFF preconditioner path is byte-untouched.
    """

    def __init__(
        self, cell: object, *, nu1: int = 2, nu2: int = 2, omega: float = 0.8,
        min_fiber_nodes: int = 4, max_levels: int = 6, attach_fiber_quotient: bool = True,
        coarsest_sweeps: int = 20, fq_inner_iterations: int = 1000,
        assembled_coarse: bool = False,
    ) -> None:
        if nu1 <= 0 or nu2 <= 0:
            raise ValueError("nu1/nu2 must be positive")
        if nu1 != nu2:
            raise ValueError("symmetric V-cycle requires nu1 == nu2 (SPD preconditioner)")
        if not omega > 0.0:
            raise ValueError("omega must be positive")
        if min_fiber_nodes < 3:
            raise ValueError("min_fiber_nodes must be >= 3 (a triple needs 3 nodes)")
        self.device = cell.device
        self.n_total = int(cell.n_total)
        self.n_actin = int(getattr(cell, "n_actin", 0)) or int(cell.foff_d.numpy()[-1])
        self.n_fibers = int(cell.n_fibers)
        self.nu1 = int(nu1)
        self.nu2 = int(nu2)
        self.omega = float(omega)
        self.coarsest_sweeps = int(coarsest_sweeps)
        # SPEED: when True, coarse levels (l >= 1) apply the cheap assembled per-fiber operator
        # (fiber_perfiber_operator_kernel, O(n_nodes[l])) instead of the full-fine R ∘ _operator ∘ P (a 2.9M-node
        # apply at EVERY level). Level 0 always uses the exact fine operator. Inter-fiber coupling stays with the
        # additive fiber-quotient coarsest level. Default False keeps the currently-native-validated V-cycle
        # byte-for-byte; the CPU h-independence/SPD gates cover the True path.
        self.assembled_coarse = bool(assembled_coarse)

        foff0 = np.ascontiguousarray(cell.foff_d.numpy(), dtype=np.int64)
        toff0 = np.ascontiguousarray(cell.toff_d.numpy(), dtype=np.int64)
        alpha0 = np.ascontiguousarray(cell.alpha_d.numpy(), dtype=np.float64)

        # ---- Host: build the level hierarchy (foff / toff / alpha per level; R/P tables per transition). ----
        level_foff = [foff0]
        level_toff = [np.ascontiguousarray(toff0, dtype=np.int32)]
        level_alpha = [np.ascontiguousarray(alpha0, dtype=np.float64)]
        self._tables_host: list[dict | None] = [None]  # tables[k] = transition (k-1) -> k
        cur_foff, cur_toff, cur_alpha = foff0, level_toff[0], level_alpha[0]
        while (cur_foff.shape[0] > 1 and int(np.diff(cur_foff).max()) > min_fiber_nodes
               and len(self._tables_host) <= max_levels):
            table = build_arclength_restriction(cur_foff)
            foff_c = np.ascontiguousarray(table["foff_coarse"], dtype=np.int64)
            toff_c, alpha_c = _coarsen_bending(cur_foff, cur_toff, cur_alpha, foff_c)
            self._tables_host.append(table)
            level_foff.append(foff_c)
            level_toff.append(toff_c)
            level_alpha.append(alpha_c)
            cur_foff, cur_toff, cur_alpha = foff_c, toff_c, alpha_c
        self.n_levels = len(self._tables_host) - 1  # coarsest arclength level index
        self.level_lengths = [int(np.diff(f).max()) if f.shape[0] > 1 else 0 for f in level_foff]

        d = self.device
        self.n_nodes = [int(f[-1]) for f in level_foff]  # nodes per level (n_nodes[0] == n_actin)

        # ---- Device: per-level operator/smoother arrays. ----
        self.foff_d = [wp.array(np.ascontiguousarray(f, np.int32), dtype=wp.int32, device=d) for f in level_foff]
        self.toff_d = [wp.array(t, dtype=wp.int32, device=d) for t in level_toff]
        self.alpha_d = [wp.array(a, dtype=wp.float64, device=d) for a in level_alpha]
        self.tables_d = [None]
        for k in range(1, self.n_levels + 1):
            tb = self._tables_host[k]
            self.tables_d.append({
                "c0": wp.array(tb["c0"], dtype=wp.int32, device=d),
                "w0": wp.array(tb["w0"], dtype=wp.float64, device=d),
                "c1": wp.array(tb["c1"], dtype=wp.int32, device=d),
                "w1": wp.array(tb["w1"], dtype=wp.float64, device=d),
            })

        # V-cycle state (level 0 working vectors are n_actin; the operator I/O buffers t1[0]/t2[0] are n_total).
        def _z(n: int) -> wp.array:
            return wp.zeros(max(n, 1), dtype=wp.vec3d, device=d)

        def _zf(n: int) -> wp.array:
            return wp.zeros(max(n, 1), dtype=wp.float64, device=d)

        self.x = [_z(n) for n in self.n_nodes]
        self.b = [_z(n) for n in self.n_nodes]
        self.res = [_z(n) for n in self.n_nodes]
        self.t1 = [_z(self.n_total if l == 0 else self.n_nodes[l]) for l in range(self.n_levels + 1)]
        self.t2 = [_z(self.n_total if l == 0 else self.n_nodes[l]) for l in range(self.n_levels + 1)]
        # Smoother external diagonal per level: ext[0] is aliased to the live coarse_external_diagonal in
        # build(); ext[l>=1] are the Galerkin-restricted diagonals rebuilt every build().
        self.ext = [None] + [_z(self.n_nodes[l]) for l in range(1, self.n_levels + 1)]
        self.diag_bend = [_zf(n) for n in self.n_nodes]
        self.off_sub = [_zf(n) for n in self.n_nodes]
        self.cprime = [_z(n) for n in self.n_nodes]
        self.dprime = [_z(n) for n in self.n_nodes]

        # ---- Coarsest inter-fiber level: assembled fiber-quotient A_c (Path B), only if the cell couples. ----
        self.fq = None
        can_fq = (
            attach_fiber_quotient and int(getattr(cell, "n_xl", 0)) > 0
            and getattr(cell, "xl_d", None) is not None
            and getattr(cell, "node_fiber_d", None) is not None)
        if can_fq:
            from aleph.components.incumbent.fiber_quotient_coarse import FiberQuotientCoarsePathB
            self.fq = FiberQuotientCoarsePathB(cell, inner_iterations=fq_inner_iterations)

        # Smoother-damping power-iteration workspace + device omega (estimated per build). The line smoother is
        # damped block-Jacobi, stable/SPD only for omega < 2/lambda_max(M_line^-1 A); a fixed omega is unsafe
        # when inter-fiber coupling is strong, so omega is estimated at setup (design §2.4) and kept on device.
        self.pw_v = _z(self.n_actin)
        self.pw_av = _z(self.n_actin)
        self.pw_w = _z(self.n_actin)
        self.pw_num = _zf(1)
        self.pw_den = _zf(1)
        self.pw_nrm = _zf(1)
        self.omega_dev = wp.zeros(1, dtype=wp.float64, device=d)
        self.omega_dev.fill_(wp.float64(self.omega))  # safe default until the first build() estimate
        self.smoother_power_iters = 10

        self._pos = None
        self._reg = None
        self._finite = None
        self.operator: Callable | None = None
        _LOG.info(
            "fiber-arclength multigrid: %d arclength levels (lengths %s), n_actin=%d, coarsest=%s",
            self.n_levels, self.level_lengths, self.n_actin,
            "fiber-quotient A_c" if self.fq is not None else f"{self.coarsest_sweeps} line-smooth sweeps")

    # -- build: restrict the live external diagonal down the levels + estimate omega + assemble the coarsest A_c. --
    def build(self, pos: wp.array, coarse_external_diagonal: wp.array, regularization: wp.array,
              finite: wp.array, operator: Callable) -> None:
        """Refresh the position-dependent smoother diagonals + omega estimate + coarsest ``A_c`` (per build)."""
        d = self.device
        self._pos = pos
        self._reg = regularization
        self._finite = finite
        self.operator = operator
        self.ext[0] = coarse_external_diagonal
        for l in range(1, self.n_levels + 1):
            tb = self.tables_d[l]
            self.ext[l].zero_()
            wp.launch(restrict_diagonal_sq_kernel, dim=self.n_nodes[l - 1],
                      inputs=[self.ext[l - 1], tb["c0"], tb["w0"], tb["c1"], tb["w1"], self.ext[l]], device=d)
        self._estimate_smoother_omega(safety=0.9)
        if self.fq is not None:
            self.fq.build(pos, regularization, finite)

    def _estimate_smoother_omega(self, *, safety: float) -> None:
        r"""Device power iteration -> ``omega = min(self.omega, safety/lambda_max(M_line^-1 A_0))`` (SPD smoother).

        Power iteration with the line-smoother solve as the operator drives ``pw_v`` to the dominant eigenvector
        of ``M_line^-1 A``; the A-inner-product Rayleigh quotient ``<A v, M^-1 A v> / <v, A v>`` converges to
        ``lambda_max`` from below, so ``safety < 1`` lands ``omega`` strictly inside the unconditional-SPD
        damped-block-Jacobi window. All scalars stay device-resident (no host readback in the physical-time
        loop). ``M_line^-1 y`` is one line-smoother sweep from ``x = 0`` at ``omega = 1``."""
        d = self.device
        if self.n_fibers == 0 or self.n_actin == 0:
            return
        self.pw_v.fill_(wp.vec3d(1.0, 1.0, 1.0))
        for _ in range(self.smoother_power_iters):
            self._apply_level_operator(0, self.pw_v, self.pw_av)          # av = A_0 v
            self.pw_w.zero_()
            wp.launch(fiber_line_smooth_kernel, dim=self.n_fibers,        # w = M_line^-1 av (x=0, omega=1)
                      inputs=[self.pw_av, self.pw_w, self.foff_d[0], self.toff_d[0], self.alpha_d[0],
                              self.ext[0], wp.float64(1.0), self.diag_bend[0], self.off_sub[0],
                              self.cprime[0], self.dprime[0], self._finite], device=d)
            self.pw_num.zero_(); self.pw_den.zero_()
            wp.launch(_mg_dot_kernel, dim=self.n_actin, inputs=[self.pw_av, self.pw_w, self.pw_num], device=d)
            wp.launch(_mg_dot_kernel, dim=self.n_actin, inputs=[self.pw_v, self.pw_av, self.pw_den], device=d)
            self.pw_nrm.zero_()
            wp.launch(_mg_dot_kernel, dim=self.n_actin, inputs=[self.pw_w, self.pw_w, self.pw_nrm], device=d)
            wp.launch(_mg_normalize_kernel, dim=self.n_actin,
                      inputs=[self.pw_w, self.pw_nrm, self.pw_v], device=d)
        wp.launch(_mg_set_omega_kernel, dim=1,
                  inputs=[self.pw_num, self.pw_den, wp.float64(float(safety)),
                          wp.float64(self.omega), self.omega_dev], device=d)

    # -- matrix-free level operator: xout(level l) = A_l xin = R_{0..l} ( _operator ( P_{l..0} xin ) ). --
    def _apply_level_operator(self, l: int, xin: wp.array, xout: wp.array) -> None:
        d = self.device
        if l >= 1 and self.assembled_coarse:
            # Cheap assembled coarse operator (SPEED): ext[l] ⊙ x + pentadiagonal bending, O(n_nodes[l]).
            # Replaces the full-fine prolong→_operator→restrict; inter-fiber coupling is the additive fq level.
            wp.launch(fiber_perfiber_operator_kernel, dim=self.n_fibers,
                      inputs=[xin, self.foff_d[l], self.toff_d[l], self.alpha_d[l], self.ext[l], xout],
                      device=d)
            return
        if l == 0:
            self.t1[0].zero_()
            wp.launch(_mg_copy_kernel, dim=self.n_actin, inputs=[xin, self.t1[0]], device=d)
        else:
            wp.copy(self.t1[l], xin)
            for m in range(l, 0, -1):
                self.t1[m - 1].zero_()
                tb = self.tables_d[m]
                wp.launch(arclength_prolong_add_kernel, dim=self.n_nodes[m - 1],
                          inputs=[self.t1[m], tb["c0"], tb["w0"], tb["c1"], tb["w1"], self.t1[m - 1]], device=d)
        self.operator(self._pos, self.t1[0], self.t2[0], self._reg, self._finite)
        if l == 0:
            wp.launch(_mg_copy_kernel, dim=self.n_actin, inputs=[self.t2[0], xout], device=d)
        else:
            src = self.t2[0]
            for m in range(1, l + 1):
                self.t2[m].zero_()
                tb = self.tables_d[m]
                wp.launch(arclength_restrict_kernel, dim=self.n_nodes[m - 1],
                          inputs=[src, tb["c0"], tb["w0"], tb["c1"], tb["w1"], self.t2[m]], device=d)
                src = self.t2[m]
            wp.copy(xout, self.t2[l])

    def _smooth(self, l: int, sweeps: int, x_is_zero: bool) -> None:
        """``sweeps`` damped line-relaxation sweeps ``x_l += omega M_l^-1 (b_l - A_l x_l)`` at level ``l``."""
        d = self.device
        for s in range(sweeps):
            if x_is_zero and s == 0:
                wp.copy(self.res[l], self.b[l])
            else:
                self._apply_level_operator(l, self.x[l], self.res[l])
                wp.launch(_mg_sub_kernel, dim=self.n_nodes[l],
                          inputs=[self.b[l], self.res[l], self.res[l]], device=d)
            # Pre-scale the residual by the device omega (M^-1 linear ⇒ omega M^-1 r = M^-1 (omega r)); the
            # smoother then runs at the constant omega = 1 (its committed signature) with the SAFE damping.
            wp.launch(_mg_scale_kernel, dim=self.n_nodes[l], inputs=[self.res[l], self.omega_dev], device=d)
            wp.launch(fiber_line_smooth_kernel, dim=self.n_fibers,
                      inputs=[self.res[l], self.x[l], self.foff_d[l], self.toff_d[l], self.alpha_d[l],
                              self.ext[l], wp.float64(1.0), self.diag_bend[l], self.off_sub[l],
                              self.cprime[l], self.dprime[l], self._finite], device=d)

    def _prolong_to_fine(self, l: int, source: wp.array) -> None:
        """Prolong a level-``l`` vector all the way to the fine actin block of ``t1[0]`` (non-actin zeroed)."""
        d = self.device
        if l == 0:
            self.t1[0].zero_()
            wp.launch(_mg_copy_kernel, dim=self.n_actin, inputs=[source, self.t1[0]], device=d)
            return
        wp.copy(self.t1[l], source)
        for m in range(l, 0, -1):
            self.t1[m - 1].zero_()
            tb = self.tables_d[m]
            wp.launch(arclength_prolong_add_kernel, dim=self.n_nodes[m - 1],
                      inputs=[self.t1[m], tb["c0"], tb["w0"], tb["c1"], tb["w1"], self.t1[m - 1]], device=d)

    def _restrict_from_fine(self, l: int, dest: wp.array) -> None:
        """Restrict the fine actin field in ``t2[0]`` down to level ``l`` and copy it into ``dest``."""
        d = self.device
        if l == 0:
            wp.launch(_mg_copy_kernel, dim=self.n_actin, inputs=[self.t2[0], dest], device=d)
            return
        src = self.t2[0]
        for m in range(1, l + 1):
            self.t2[m].zero_()
            tb = self.tables_d[m]
            wp.launch(arclength_restrict_kernel, dim=self.n_nodes[m - 1],
                      inputs=[src, tb["c0"], tb["w0"], tb["c1"], tb["w1"], self.t2[m]], device=d)
            src = self.t2[m]
        wp.copy(dest, self.t2[l])

    def _coarsest_solve(self) -> None:
        r"""Solve the coarsest arclength level ``A_L = R_{0..L} A P_{L..0}`` with damped line-smoother Richardson.

        This is the Galerkin-CONSISTENT terminal solve of the arclength hierarchy (``coarsest_sweeps`` SPD line
        relaxations). The fiber-quotient rigid inter-fiber space is a DIFFERENT (non-arclength) coarse space, so
        it is NOT nested here — nesting a non-Galerkin-consistent coarse solve inside the recursion over-corrects
        and compounds through the post-smooths. It is instead injected ADDITIVELY at the fine level in
        :meth:`apply` (``M^-1 = V_arclength + P A_c^-1 P^T``, the proven-SPD additive form used by the incumbent
        fiber-quotient coarse and the ``FIBER_QUOTIENT_COARSE_PLAN`` scope).
        """
        L = self.n_levels
        self.x[L].zero_()
        self._smooth(L, self.coarsest_sweeps, x_is_zero=True)

    def _vcycle(self, l: int) -> None:
        d = self.device
        if l == self.n_levels:
            self._coarsest_solve()
            return
        self.x[l].zero_()
        self._smooth(l, self.nu1, x_is_zero=True)
        # r_l = b_l - A_l x_l
        self._apply_level_operator(l, self.x[l], self.res[l])
        wp.launch(_mg_sub_kernel, dim=self.n_nodes[l],
                  inputs=[self.b[l], self.res[l], self.res[l]], device=d)
        # b_{l+1} = R r_l
        self.b[l + 1].zero_()
        tb = self.tables_d[l + 1]
        wp.launch(arclength_restrict_kernel, dim=self.n_nodes[l],
                  inputs=[self.res[l], tb["c0"], tb["w0"], tb["c1"], tb["w1"], self.b[l + 1]], device=d)
        self._vcycle(l + 1)
        # x_l += P e_{l+1}
        wp.launch(arclength_prolong_add_kernel, dim=self.n_nodes[l],
                  inputs=[self.x[l + 1], tb["c0"], tb["w0"], tb["c1"], tb["w1"], self.x[l]], device=d)
        self._smooth(l, self.nu2, x_is_zero=False)

    def apply(self, pos: wp.array, residual: wp.array, preconditioned: wp.array,
              operator: Callable, regularization: wp.array, finite: wp.array) -> None:
        r"""OVERWRITE the actin block of ``preconditioned`` with the V-cycle ``M^-1 residual`` (SPD).

        ``residual`` and ``preconditioned`` are full ``n_total`` vectors; only the leading ``n_actin`` fiber
        block is read/written (the caller node-Jacobi-preconditions the disjoint remainder). The stored ``pos``
        / ``regularization`` / ``finite`` and the fine ``operator`` are used by the recursive Galerkin coarse
        operators.
        """
        d = self.device
        self._pos = pos
        self._reg = regularization
        self._finite = finite
        self.operator = operator
        wp.launch(_mg_copy_kernel, dim=self.n_actin, inputs=[residual, self.b[0]], device=d)
        self._vcycle(0)
        wp.launch(_mg_copy_kernel, dim=self.n_actin, inputs=[self.x[0], preconditioned], device=d)
        # Additive fine-level fiber-quotient rigid inter-fiber correction: M^-1 = V_arclength + P A_c^-1 P^T.
        # Both terms are SPD ⇒ the composite is SPD; applied ONCE at the fine level (not nested in the arclength
        # recursion) so a non-Galerkin-consistent A_c cannot over-correct/compound through the post-smooths.
        # apply_richardson ADDS its correction into the actin block of `preconditioned` (the incumbent form).
        if self.fq is not None:
            self.fq.apply_richardson(pos, residual, preconditioned, regularization, finite,
                                     sweeps=self.coarsest_sweeps)
