"""Standalone CPU prototype: fiber-quotient inter-fiber coarse correction for the GATE-A stall.

Context
-------
The resting cortex static solve (GATE A) plateaus at ``max|PF| ~ 1.16 pN``. A solver-expert
consultation traced this to low-frequency *inter-fiber* error modes: on a very stiff crosslinked
network (crosslink tangent ~8e5 pN/um >> everything), ~4,420 spatially-distributed contractile
motor loads excite collective rearrangements where whole filaments move nearly rigidly, coupled
only through crosslinks. The production preconditioners in
``aleph/components/incumbent/implicit_mechanics.py`` are:

  * per-node Jacobi                         (``apply_preconditioner_kernel``)
  * per-fiber dense Cholesky block          (``build_fiber_block_cholesky_kernel``)
  * per-fiber DIAGONAL translation Rayleigh (``add_fiber_translation_coarse_kernel``) -- independent
    per fiber, NO inter-fiber coupling in the coarse solve
  * a <=12-mode GLOBAL rigid+strain (l<=2) two-level Galerkin coarse
    (``rigid_strain_coarse_basis`` + ``assemble/restrict/solve/prolong`` kernels)

All four resolve intra-fiber modes and the 12 global l<=2 modes, but NONE couples fibers to each
other at the coarse level, so the collective inter-fiber near-null space is left un-damped -> plateau.

This prototype reproduces the essential difficulty in a small synthetic network and shows that a
FIBER-QUOTIENT coarse correction (collapse each ~L-node filament to its per-fiber rigid-body modes,
build the crosslink-weighted Galerkin coarse operator ``A_c = P^T K P`` on the SPARSE per-fiber
prolongator, solve, prolong) closes the plateau that per-node / per-fiber / global-12 cannot.

Solvers compared (all share the per-fiber block-Cholesky smoother except A0):
  A0  node Jacobi PCG                              -> plateaus (baseline)
  A1  per-fiber block-Cholesky PCG                 -> plateaus (mirrors production block)
  A2  A1 + diagonal per-fiber translation Rayleigh -> plateaus (mirrors add_fiber_translation_coarse)
  A3  A1 + global 12-mode l<=2 Galerkin coarse     -> plateaus (mirrors rigid_strain_coarse_basis)
  B   A1 + FIBER-QUOTIENT inter-fiber Galerkin coarse  -> CONVERGES (the missing piece)

Every coarse correction is added ADDITIVELY (M^-1 = S^-1 + C A_c^-1 C^T), which stays SPD so PCG is
valid -- exactly the additive design used on device. This is preconditioning only; it moves the
convergence path, never the fixed point.

Sanity Gate
-----------
  * Dimensions / SPD: K = sum_e k_e (u_e u_e^T) block (SPSD, k_e>0) + a*I (a>0) is symmetric PD;
    verified numerically on a small twin (dense eigvalsh, min eig > 0).
  * Rigid-mode content: each fiber's coarse block reproduces an exact rigid motion of that fiber to
    ~1e-10; a straight fiber yields rank EXACTLY 5 (3 translations + 2 rotations; the axial rotation
    is correctly dropped as a null column) -- matches the consultation's "3 translations + 2 nonzero
    rotations for an approximately straight chain".
  * Near-null capture: the softest eigenvectors of K (the collective floppy modes) are captured by
    the fiber-quotient range to ~1e-6 but NOT by the global-12 space (projection residual O(1)) --
    this is the quantitative reason B converges where A3 plateaus.
  * Galerkin symmetry: A_c = P^T K P is symmetric to round-off.
  * No tuned constant: the coarse builders take only geometry + fiber topology (no relaxation
    parameter, no fitting knob); the rank tolerance is a numerical gauge relative to the matrix
    scale, not a physics constant. All physical constants are declared in ``ModelParams`` and printed.

Run:  ~/miniconda3/envs/ffn_sim/bin/python aleph/scripts/ac_fiber_quotient_coarse_proto.py
Pure CPU numpy/scipy. Imports NOTHING from the Warp runtime; modifies no simulation code.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

# ----------------------------------------------------------------------------------------------------
# Model definition (all constants physical / geometric; no solver knob is tuned to make a gate pass).
# ----------------------------------------------------------------------------------------------------


@dataclass
class ModelParams:
    """Synthetic stiff-crosslinked fiber network capturing the GATE-A difficulty.

    Attributes:
        n_fibers: number of short filaments (each a nearly-rigid braced chain).
        nodes_per_fiber: nodes per filament (7 mirrors the native cortical discretization).
        fiber_length: filament arclength [arb. length unit].
        box: cubic domain edge; sets crosslink density.
        k_backbone: axial (consecutive) central-force spring stiffness -> fiber stretch stiffness.
        k_bend: NF2007-style bending-stencil stiffness (``alpha (D2^T D2)``) -> transverse rigidity of
            the straight chain (mirrors ``bending_stiffness`` in implicit_mechanics_analytic). Axial +
            bending make each straight fiber a nearly rigid worm-like rod whose free near-null space is
            EXACTLY the 5 rigid modes (3 translations + 2 transverse rotations).
        k_crosslink: inter-fiber crosslink spring stiffness (the STIFF scale, contrast vs a).
        regularizer: the soft ``a`` in ``aI + PKP``; the ONLY restoring term on the collective
            near-rigid (floppy) inter-fiber modes -> contrast = k_crosslink / regularizer.
        crosslink_cutoff: KD-tree radius for candidate inter-fiber crosslinks.
        crosslink_fraction: fraction of the fiber-quotient isostatic count (~5*n_fibers) to keep,
            <1 -> a genuinely floppy (sub-isostatic) assembly with many collective soft modes.
        n_load_pairs: scattered opposing point-load pairs (the motor analog; the hard part).
        smooth_load: amplitude of the easy smooth body load.
        pair_load: magnitude of each concentrated opposing pair.
        seed: RNG seed (reproducible).
    """

    n_fibers: int = 600
    nodes_per_fiber: int = 7
    fiber_length: float = 1.0
    box: float = 12.0
    k_backbone: float = 1.0e6
    k_bend: float = 1.0e6
    k_crosslink: float = 1.0e6
    regularizer: float = 1.0
    crosslink_cutoff: float = 0.55
    crosslink_fraction: float = 0.65
    n_load_pairs: int = 2000
    smooth_load: float = 5.0
    pair_load: float = 50.0
    seed: int = 20260724

    @property
    def contrast(self) -> float:
        return self.k_crosslink / self.regularizer


@dataclass
class Network:
    """Assembled synthetic network."""

    pos: np.ndarray                      # (N, 3) node positions
    fiber_off: np.ndarray                # (n_fibers+1,) CSR fiber->node offsets
    node_fiber: np.ndarray               # (N,) fiber id per node
    springs_i: np.ndarray                # (S,) endpoint a
    springs_j: np.ndarray                # (S,) endpoint b
    springs_k: np.ndarray                # (S,) stiffness
    springs_u: np.ndarray                # (S, 3) unit direction
    is_crosslink: np.ndarray             # (S,) bool -- inter-fiber crosslink vs backbone
    bend_triples: np.ndarray             # (T, 3) consecutive node triples per fiber
    bend_alpha: np.ndarray               # (T,) bending stencil stiffness
    params: ModelParams
    n_crosslinks: int = 0
    isostatic: int = 0

    @property
    def n_nodes(self) -> int:
        return self.pos.shape[0]

    @property
    def n_dof(self) -> int:
        return 3 * self.pos.shape[0]


def build_network(params: ModelParams) -> Network:
    """Build many short nearly-rigid fibers, stiff-crosslinked sub-isostatically into one floppy web."""
    rng = np.random.default_rng(params.seed)
    m, L = params.n_fibers, params.nodes_per_fiber
    seg = params.fiber_length / (L - 1)

    # Fiber centers + random orientations; nodes are an EXACTLY straight chain so the fiber's rigid
    # near-null space is exactly 5-D (the axial rotation of collinear nodes is the trivial zero mode);
    # transverse rigidity comes from the NF2007 bending stencil (below), not from geometry.
    centers = rng.uniform(0.0, params.box, size=(m, 3))
    dirs = rng.normal(size=(m, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    s = (np.arange(L) - (L - 1) / 2.0)[None, :, None] * seg           # (1, L, 1)
    pos = (centers[:, None, :] + s * dirs[:, None, :]).reshape(m * L, 3)

    node_fiber = np.repeat(np.arange(m), L)
    fiber_off = np.arange(0, m * L + 1, L)
    base = np.arange(m) * L

    # Backbone: consecutive AXIAL central-force springs (stretch stiffness).
    bi = np.concatenate([base + p for p in range(L - 1)])
    bj = np.concatenate([base + p + 1 for p in range(L - 1)])

    # Backbone: NF2007 bending stencil over consecutive triples (transverse rigidity, keeps fiber straight).
    tri = np.stack([np.concatenate([base + p for p in range(L - 2)]),
                    np.concatenate([base + p + 1 for p in range(L - 2)]),
                    np.concatenate([base + p + 2 for p in range(L - 2)])], axis=1)
    bend_alpha = np.full(tri.shape[0], params.k_bend, dtype=np.float64)

    # Crosslinks: KD-tree candidate pairs on DIFFERENT fibers, subsampled sub-isostatic, then unioned
    # with an MST over fibers so the whole web is one connected component (else the demo decouples).
    tree = cKDTree(pos)
    pairs = np.array(sorted(tree.query_pairs(params.crosslink_cutoff)), dtype=np.int64)
    inter = node_fiber[pairs[:, 0]] != node_fiber[pairs[:, 1]]
    pairs = pairs[inter]
    isostatic = 5 * m - 6                                            # fiber-quotient Maxwell count
    target = int(params.crosslink_fraction * isostatic)
    if pairs.shape[0] > target:
        keep = rng.choice(pairs.shape[0], size=target, replace=False)
        pairs = pairs[keep]
    mi, mj = _spanning_crosslinks(pos, node_fiber, fiber_off, m)
    xi = np.concatenate([pairs[:, 0], mi])
    xj = np.concatenate([pairs[:, 1], mj])
    key = np.minimum(xi, xj).astype(np.int64) * (m * L) + np.maximum(xi, xj)
    _, uniq = np.unique(key, return_index=True)
    xi, xj = xi[uniq], xj[uniq]

    springs_i = np.concatenate([bi, xi])
    springs_j = np.concatenate([bj, xj])
    is_crosslink = np.concatenate([np.zeros(bi.size, bool), np.ones(xi.size, bool)])
    stiff = np.where(is_crosslink, params.k_crosslink, params.k_backbone)

    delta = pos[springs_j] - pos[springs_i]
    length = np.linalg.norm(delta, axis=1, keepdims=True)
    unit = delta / length

    return Network(
        pos=pos, fiber_off=fiber_off, node_fiber=node_fiber,
        springs_i=springs_i, springs_j=springs_j, springs_k=stiff, springs_u=unit,
        is_crosslink=is_crosslink, bend_triples=tri, bend_alpha=bend_alpha, params=params,
        n_crosslinks=int(xi.size), isostatic=isostatic,
    )


def _spanning_crosslinks(
    pos: np.ndarray, node_fiber: np.ndarray, fiber_off: np.ndarray, m: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return one crosslink per MST edge over fiber centroids (guarantees a single connected web)."""
    centroids = np.stack([pos[fiber_off[f]:fiber_off[f + 1]].mean(axis=0) for f in range(m)])
    mst = minimum_spanning_tree(cdist(centroids, centroids)).tocoo()  # m-1 spanning edges
    xi, xj = [], []
    for fa, fb in zip(mst.row, mst.col):
        na = np.arange(fiber_off[fa], fiber_off[fa + 1])
        nb = np.arange(fiber_off[fb], fiber_off[fb + 1])
        d = np.linalg.norm(pos[na][:, None, :] - pos[nb][None, :, :], axis=2)
        a, b = np.unravel_index(np.argmin(d), d.shape)
        xi.append(int(na[a]))
        xj.append(int(nb[b]))
    return np.asarray(xi, np.int64), np.asarray(xj, np.int64)


# ----------------------------------------------------------------------------------------------------
# Stiffness assembly (rank-1 central-force spring tangents at the reference geometry + regularizer).
# ----------------------------------------------------------------------------------------------------


def assemble_stiffness(net: Network) -> tuple[sp.csr_matrix, np.ndarray]:
    """Return SPD ``K = sum_e k_e (u_e u_e^T) + sum_t alpha_t (D2^T D2) + a I`` and the external diag.

    Central-force springs contribute the SPSD rank-1 tangent ``k_e [u;-u][u;-u]^T``; bending triples
    contribute the NF2007 stencil ``alpha (D2^T D2)`` (componentwise ``w w^T``, ``w=[1,-2,1]``). The
    external diagonal is ``a`` plus only the CROSSLINK self blocks -- exactly the ``external_diagonal``
    the production per-fiber block / translation coarse consume (backbone + bending stay in the block).
    """
    a = net.params.regularizer
    n = net.n_dof
    ai, aj = np.arange(3), np.arange(3)
    rr, cc = np.meshgrid(ai, aj, indexing="ij")

    # Central-force spring blocks.
    uu = net.springs_u[:, :, None] * net.springs_u[:, None, :]       # (S, 3, 3)
    blk = net.springs_k[:, None, None] * uu                          # (S, 3, 3)

    def emit(base_r, base_c, sign):
        rows = (3 * base_r[:, None, None] + rr[None]).reshape(-1)
        cols = (3 * base_c[:, None, None] + cc[None]).reshape(-1)
        vals = (sign * blk).reshape(-1)
        return rows, cols, vals

    i, j = net.springs_i, net.springs_j
    parts = [emit(i, i, +1.0), emit(j, j, +1.0), emit(i, j, -1.0), emit(j, i, -1.0)]

    # Bending stencil: for each triple and Cartesian component, alpha * w w^T on the 3 nodes.
    w = np.array([1.0, -2.0, 1.0])
    ww = np.outer(w, w)                                              # (3, 3)
    tri = net.bend_triples                                          # (T, 3)
    bblk = net.bend_alpha[:, None, None] * ww[None]                 # (T, 3, 3)
    li, lj = np.meshgrid(np.arange(3), np.arange(3), indexing="ij")
    for comp in range(3):
        gdof = 3 * tri + comp                                       # (T, 3)
        parts.append((gdof[:, li.reshape(-1)].reshape(-1),
                      gdof[:, lj.reshape(-1)].reshape(-1),
                      bblk.reshape(tri.shape[0], 9).reshape(-1)))

    rows = np.concatenate([p[0] for p in parts])
    cols = np.concatenate([p[1] for p in parts])
    vals = np.concatenate([p[2] for p in parts])
    K = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    K = K + a * sp.identity(n, format="csr")

    # External diagonal (vec3 per node): a + crosslink self blocks only.
    ext = np.full((net.n_nodes, 3), a, dtype=np.float64)
    xmask = net.is_crosslink
    diag = net.springs_k[xmask, None] * (net.springs_u[xmask] ** 2)  # (Sx, 3) diag of k uu^T
    np.add.at(ext, net.springs_i[xmask], diag)
    np.add.at(ext, net.springs_j[xmask], diag)
    return K, ext


def build_rhs(net: Network) -> np.ndarray:
    """Smooth body load (easy) + many scattered opposing concentrated pairs (the motor-analog hard part)."""
    p = net.params
    rng = np.random.default_rng(p.seed + 1)
    pos = net.pos
    f = np.zeros((net.n_nodes, 3))
    # (i) smooth body load: a divergence-free-ish smooth field, projects onto easy modes.
    kx = 2.0 * np.pi / p.box
    f[:, 0] += p.smooth_load * np.sin(kx * pos[:, 1])
    f[:, 1] += p.smooth_load * np.sin(kx * pos[:, 2])
    f[:, 2] += p.smooth_load * np.sin(kx * pos[:, 0])
    # (ii) concentrated opposing pairs between nearby DIFFERENT fibers (the contractile motor analog).
    tree = cKDTree(pos)
    pairs = np.array(sorted(tree.query_pairs(p.crosslink_cutoff * 1.5)), dtype=np.int64)
    inter = net.node_fiber[pairs[:, 0]] != net.node_fiber[pairs[:, 1]]
    pairs = pairs[inter]
    sel = rng.choice(pairs.shape[0], size=min(p.n_load_pairs, pairs.shape[0]), replace=False)
    pairs = pairs[sel]
    d = rng.normal(size=(pairs.shape[0], 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    np.add.at(f, pairs[:, 0], +p.pair_load * d)
    np.add.at(f, pairs[:, 1], -p.pair_load * d)
    return f.reshape(-1)


# ----------------------------------------------------------------------------------------------------
# Preconditioners.
# ----------------------------------------------------------------------------------------------------


def node_jacobi(K: sp.csr_matrix) -> Callable[[np.ndarray], np.ndarray]:
    """Per-node scalar Jacobi: ``z = r / diag(K)`` (mirrors apply_preconditioner_kernel)."""
    d = K.diagonal()
    d = np.where(d > 0.0, d, 1.0)
    return lambda r: r / d


def fiber_block_cholesky(K: sp.csr_matrix, net: Network) -> Callable[[np.ndarray], np.ndarray]:
    """Exact dense solve of each fiber's principal block (mirrors build_fiber_block_cholesky_kernel).

    The block = fiber backbone (axial + bending, full) + crosslink self diagonal + a*I -- the fiber's
    OWN DOFs with inter-fiber crosslink coupling dropped to the diagonal, exactly the production
    external-diagonal block. Solving it resolves intra-fiber + that fiber's local rigid adjustment,
    but cannot coordinate the collective inter-fiber rearrangement. Vectorized: fibers are contiguous
    equal-size (3L) DOF blocks, so the 600 solves batch into one einsum.
    """
    m, L = net.params.n_fibers, net.params.nodes_per_fiber
    d = 3 * L
    Kc = K.tocsc()
    inv = np.empty((m, d, d))
    for f in range(m):
        s = 3 * int(net.fiber_off[f])
        inv[f] = np.linalg.inv(Kc[s:s + d, s:s + d].toarray())
    return lambda r: np.einsum("fij,fj->fi", inv, r.reshape(m, d)).reshape(-1)


def fiber_translation_rayleigh(net: Network, ext: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """Diagonal per-fiber translation Rayleigh coarse (mirrors add_fiber_translation_coarse_kernel).

    Independent per fiber, 3 translation modes, diagonal (Rayleigh) approximation, NO inter-fiber
    coupling: ``corr_f[c] = (sum_i r_i[c]) / (sum_i ext_i[c])`` added to every node of fiber f.
    """
    m, L = net.params.n_fibers, net.params.nodes_per_fiber
    den = ext.reshape(m, L, 3).sum(axis=1)                          # (m, 3)
    den = np.where(den > 0.0, den, 1.0)

    def apply(r: np.ndarray) -> np.ndarray:
        num = r.reshape(m, L, 3).sum(axis=1)                        # (m, 3)
        corr = num / den
        return np.broadcast_to(corr[:, None, :], (m, L, 3)).reshape(-1).copy()

    return apply


def global_rigid_strain_basis(pos: np.ndarray, n_modes: int = 12) -> np.ndarray:
    """Orthonormal (3N, K) global l<=2 rigid+strain basis (reproduces rigid_strain_coarse_basis).

    3 translations + 3 rotations + 6 constant symmetric strains, QR-orthonormalized. This is the
    EXISTING global coarse space; the prototype uses it to show it does NOT close the plateau.
    """
    n = pos.shape[0]
    rel = pos - pos.mean(axis=0)
    cols = []
    eye = np.eye(3)
    for ax in range(3):                       # translations
        v = np.zeros((n, 3)); v[:, ax] = 1.0; cols.append(v)
    for ax in range(3):                       # rotations
        cols.append(np.cross(eye[ax], rel))
    for a, b in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]:  # strains
        E = np.zeros((3, 3)); E[a, b] = 1.0; E[b, a] = 1.0; cols.append(rel @ E)
    M = np.stack(cols[:n_modes], axis=0).reshape(n_modes, 3 * n).T
    q, _ = np.linalg.qr(M)
    return q  # (3N, K), orthonormal columns


def fiber_quotient_prolongator(net: Network, rank_tol: float = 1.0e-9) -> tuple[sp.csr_matrix, dict]:
    """Build the SPARSE per-fiber rigid-body prolongator P (the NEW inter-fiber coarse space).

    Each fiber contributes its local rigid modes ``u_i = t + w x (x_i - c_f)`` (3 translations + 3
    rotations). Per fiber the 3L x 6 candidate matrix is orthonormalized (thin QR) and rank-deficient
    columns (the axial rotation of a straight chain) are dropped, giving rank r_f = 5 for a straight
    fiber. Column blocks of different fibers touch disjoint node DOFs so P is exactly block-sparse and
    full column rank; it is NEVER materialized as a dense (K, N, 3) array.

    Returns:
        P: scipy CSR (3N, Nc); Nc = sum_f r_f.
        info: per-fiber ranks + diagnostics.
    """
    pos, off = net.pos, net.fiber_off
    rows, cols, vals = [], [], []
    ranks = []
    col_base = 0
    for f in range(net.params.n_fibers):
        nodes = np.arange(off[f], off[f + 1])
        x = pos[nodes]
        c = x.mean(axis=0)
        rel = x - c
        L = nodes.size
        cand = np.zeros((3 * L, 6))
        for ax in range(3):
            t = np.zeros((L, 3)); t[:, ax] = 1.0
            cand[:, ax] = t.reshape(-1)
            cand[:, 3 + ax] = np.cross(np.eye(3)[ax], rel).reshape(-1)
        # Column-pivoted QR to detect & drop the near-null (axial-rotation) column robustly.
        q, r, piv = _pivoted_qr(cand)
        diag = np.abs(np.diag(r))
        keep = diag > rank_tol * max(1.0, diag.max())
        Q = q[:, keep]                       # (3L, r_f) orthonormal fiber rigid modes
        rf = Q.shape[1]
        ranks.append(rf)
        dofs = (3 * nodes[:, None] + np.arange(3)).reshape(-1)   # (3L,)
        rr, cc = np.nonzero(Q)
        rows.append(dofs[rr])
        cols.append(col_base + cc)
        vals.append(Q[rr, cc])
        col_base += rf
    P = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(net.n_dof, col_base)).tocsr()
    info = {
        "n_coarse_dof": col_base,
        "ranks": np.asarray(ranks),
        "nnz_per_col": P.nnz / max(1, col_base),
        "dense_KN3_bytes": net.params.n_fibers * net.n_nodes * 3 * 8,  # what the dense basis WOULD cost
        "sparse_bytes": P.data.nbytes + P.indices.nbytes + P.indptr.nbytes,
    }
    return P, info


def _pivoted_qr(a: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Column-pivoted QR via scipy (Q, R, piv) for robust rank detection."""
    from scipy.linalg import qr as sqr
    q, r, piv = sqr(a, mode="economic", pivoting=True)
    return q, r, piv


class GalerkinCoarse:
    """Additive two-level Galerkin coarse correction ``z = C A_c^{-1} C^T r`` with ``A_c = C^T K C``.

    ``C`` may be sparse (fiber-quotient) or dense (global l<=2). SPD because K is SPD and C is full
    column rank; adding it to an SPD smoother keeps the composite SPD -> PCG stays valid.
    """

    def __init__(self, K: sp.csr_matrix, C) -> None:
        self.sparse = sp.issparse(C)
        if self.sparse:
            self.C = C.tocsr()
            Ac = (self.C.T @ (K @ self.C)).tocsc()
            self.Ac = Ac
            self.solve = spla.factorized(Ac)          # sparse LU factor (SuperLU)
        else:
            self.C = np.asarray(C)
            Ac = self.C.T @ (K @ self.C)
            self.Ac = Ac
            self._chol = np.linalg.cholesky(Ac)

    def apply(self, r: np.ndarray) -> np.ndarray:
        rc = self.C.T @ r
        if self.sparse:
            yc = self.solve(rc)
        else:
            y = np.linalg.solve(self._chol, rc)
            yc = np.linalg.solve(self._chol.T, y)
        return self.C @ yc

    def symmetry(self) -> float:
        Ac = self.Ac.toarray() if self.sparse else self.Ac
        return float(np.abs(Ac - Ac.T).max() / max(1.0, np.abs(Ac).max()))


def additive(*precs: Callable[[np.ndarray], np.ndarray]) -> Callable[[np.ndarray], np.ndarray]:
    """Additive combination of SPD preconditioners: ``M^-1 r = sum_k P_k r``."""
    return lambda r: sum(p(r) for p in precs)


# ----------------------------------------------------------------------------------------------------
# Preconditioned CG (records true relative residual + max nodal residual, mirrors max|PF|).
# ----------------------------------------------------------------------------------------------------


def pcg(
    K: sp.csr_matrix, b: np.ndarray, Minv: Callable[[np.ndarray], np.ndarray],
    *, tol: float = 1.0e-8, maxiter: int = 1500,
) -> dict:
    """Left-preconditioned CG. Returns history, iters-to-tol, and residual at the device budget (32)."""
    x = np.zeros_like(b)
    r = b.copy()
    bnorm = np.linalg.norm(b)
    z = Minv(r)
    p = z.copy()
    rz = float(r @ z)
    hist = [1.0]
    iters_to_tol = None
    k = 0
    for k in range(1, maxiter + 1):
        Ap = K @ p
        pAp = float(p @ Ap)
        if pAp <= 0.0:
            break
        alpha = rz / pAp
        x += alpha * p
        r -= alpha * Ap
        rel = float(np.linalg.norm(r) / bnorm)
        hist.append(rel)
        if rel <= tol:
            iters_to_tol = k
            break
        z = Minv(r)
        rz_new = float(r @ z)
        beta = rz_new / rz
        p = z + beta * p
        rz = rz_new
    hist = np.asarray(hist)
    max_nodal = float(np.linalg.norm(r.reshape(-1, 3), axis=1).max())
    return {
        "history": hist,
        "iters": k,
        "iters_to_tol": iters_to_tol,
        "relres": float(hist[-1]),
        "res_at_32": float(hist[min(32, hist.size - 1)]),
        "max_nodal": max_nodal,
    }


# ----------------------------------------------------------------------------------------------------
# Sanity Gate.
# ----------------------------------------------------------------------------------------------------


def sanity_gate(params: ModelParams) -> list[str]:
    """Run the Sanity Gate checks on a small twin; return human-readable PASS/FAIL lines."""
    out = []
    small = ModelParams(**{**params.__dict__})
    small.n_fibers = 24
    small.crosslink_fraction = 0.6
    net = build_network(small)
    K, ext = assemble_stiffness(net)
    Kd = K.toarray()

    # (1) symmetry + SPD (dense eigvalsh on the small twin).
    asym = np.abs(Kd - Kd.T).max() / max(1.0, np.abs(Kd).max())
    evals = np.linalg.eigvalsh(Kd)
    spd = evals.min()
    out.append(f"[{'PASS' if asym < 1e-10 else 'FAIL'}] K symmetric (rel asym {asym:.2e})")
    out.append(f"[{'PASS' if spd > 0.0 else 'FAIL'}] K SPD (min eig {spd:.3e} > 0; "
               f"= regularizer a={small.regularizer:g} as expected)")

    # (2) fiber-quotient rigid-mode content + straight-fiber rank == 5.
    P, info = fiber_quotient_prolongator(net)
    ranks = info["ranks"]
    rank5 = np.all(ranks == 5)
    out.append(f"[{'PASS' if rank5 else 'FAIL'}] every straight fiber has coarse rank EXACTLY 5 "
               f"(3 trans + 2 rot); observed unique ranks {sorted(set(ranks.tolist()))}")
    # exact rigid motion of fiber 0 must lie in range(P).
    f0 = np.arange(net.fiber_off[0], net.fiber_off[1])
    rel = net.pos[f0] - net.pos[f0].mean(axis=0)
    v = np.zeros((net.n_nodes, 3))
    v[f0] = np.array([0.3, -0.2, 0.1]) + np.cross(np.array([0.05, 0.11, -0.07]), rel)
    v = v.reshape(-1)
    Pd = P.toarray()
    proj = Pd @ np.linalg.lstsq(Pd, v, rcond=None)[0]
    rig_res = np.linalg.norm(proj - v) / np.linalg.norm(v)
    out.append(f"[{'PASS' if rig_res < 1e-9 else 'FAIL'}] exact fiber rigid motion reproduced by P "
               f"(residual {rig_res:.2e})")

    # (3) near-null capture: softest eigvecs of K in range(P) but NOT in global-12.
    n_soft = 40
    idx = np.argsort(evals)[:n_soft]
    Vs = np.linalg.eigh(Kd)[1][:, idx]           # (3N, n_soft) softest modes
    G = global_rigid_strain_basis(net.pos, 12)
    def capture(B):
        Q, _ = np.linalg.qr(B if B.ndim == 2 else B)
        res = Vs - Q @ (Q.T @ Vs)
        return np.linalg.norm(res, axis=0).mean()
    cap_fq = capture(Pd)
    cap_g = capture(G)
    out.append(f"[{'PASS' if cap_fq < 1e-4 else 'FAIL'}] fiber-quotient captures softest {n_soft} "
               f"modes (mean residual {cap_fq:.2e})")
    out.append(f"[{'PASS' if cap_g > 0.5 else 'FAIL'}] global-12 does NOT capture them "
               f"(mean residual {cap_g:.2e}) -> this is why B beats A3")

    # (4) Galerkin symmetry.
    gc = GalerkinCoarse(K, P)
    out.append(f"[{'PASS' if gc.symmetry() < 1e-10 else 'FAIL'}] A_c = P^T K P symmetric "
               f"(rel asym {gc.symmetry():.2e})")

    # (5) no tuned constant: coarse builder signature carries only geometry/topology.
    import inspect
    sig = set(inspect.signature(fiber_quotient_prolongator).parameters) - {"net", "rank_tol"}
    out.append(f"[{'PASS' if not sig else 'FAIL'}] fiber-quotient builder has no physics knob "
               f"(only net + numerical rank_tol; extra params {sig or 'none'})")
    return out


# ----------------------------------------------------------------------------------------------------
# Driver.
# ----------------------------------------------------------------------------------------------------


def main() -> None:
    params = ModelParams()
    t0 = time.time()
    net = build_network(params)
    K, ext = assemble_stiffness(net)
    b = build_rhs(net)
    build_s = time.time() - t0

    print("=" * 96)
    print("FIBER-QUOTIENT INTER-FIBER COARSE -- synthetic GATE-A stall prototype")
    print("=" * 96)
    print(f"fibers={params.n_fibers}  nodes/fiber={params.nodes_per_fiber}  "
          f"nodes={net.n_nodes}  DOF={net.n_dof}")
    print(f"backbone springs + crosslinks={net.n_crosslinks}  "
          f"(fiber-quotient isostatic ~5F-6={net.isostatic}; "
          f"kept fraction {net.n_crosslinks / net.isostatic:.2f} -> "
          f"{'SUB-isostatic (floppy)' if net.n_crosslinks < net.isostatic else 'over-constrained'})")
    print(f"stiffness contrast k_crosslink/a = {params.contrast:.0e}   "
          f"(k_backbone={params.k_backbone:g}, k_crosslink={params.k_crosslink:g}, a={params.regularizer:g})")
    print(f"loads: smooth body + {params.n_load_pairs} scattered opposing pairs "
          f"(|pair|={params.pair_load:g})")
    print(f"assembly time {build_s:.2f}s")
    print()

    # Build preconditioners.
    tp = time.time()
    jac = node_jacobi(K)
    block = fiber_block_cholesky(K, net)
    transl = fiber_translation_rayleigh(net, ext)
    Gbasis = global_rigid_strain_basis(net.pos, 12)
    global_coarse = GalerkinCoarse(K, Gbasis)
    P, pinfo = fiber_quotient_prolongator(net)
    fq_coarse = GalerkinCoarse(K, P)
    prec_s = time.time() - tp
    print(f"fiber-quotient coarse: {pinfo['n_coarse_dof']} coarse DOF, "
          f"{pinfo['nnz_per_col']:.1f} nnz/col, sparse P = "
          f"{pinfo['sparse_bytes'] / 1e6:.2f} MB   vs dense (K,N,3) basis = "
          f"{pinfo['dense_KN3_bytes'] / 1e9:.1f} GB (infeasible)")
    print(f"preconditioner build time {prec_s:.2f}s")
    print()

    solvers = {
        "A0  node Jacobi": jac,
        "A1  per-fiber block Cholesky": block,
        "A2  A1 + fiber-translation Rayleigh": additive(block, transl),
        "A3  A1 + GLOBAL 12-mode l<=2 coarse": additive(block, global_coarse.apply),
        "B   A1 + FIBER-QUOTIENT coarse": additive(block, fq_coarse.apply),
    }

    results = {}
    for name, Minv in solvers.items():
        t = time.time()
        results[name] = pcg(K, b, Minv, tol=1e-8, maxiter=1500)
        results[name]["wall"] = time.time() - t

    # Table.
    print("-" * 96)
    print(f"{'solver':40s} {'res@32':>10s} {'iters->1e-8':>12s} "
          f"{'final rel':>11s} {'max nodal':>11s} {'verdict':>9s}")
    print("-" * 96)
    for name, r in results.items():
        it = str(r["iters_to_tol"]) if r["iters_to_tol"] else f">{r['iters']}"
        verdict = "converged" if r["iters_to_tol"] else "PLATEAU"
        print(f"{name:40s} {r['res_at_32']:10.2e} {it:>12s} "
              f"{r['relres']:11.2e} {r['max_nodal']:11.2e} {verdict:>9s}")
    print("-" * 96)
    print()

    # Verdict.
    a3 = results["A3  A1 + GLOBAL 12-mode l<=2 coarse"]
    b_ = results["B   A1 + FIBER-QUOTIENT coarse"]
    a3_plateau = a3["iters_to_tol"] is None
    b_conv = b_["iters_to_tol"] is not None
    speedup = (a3["res_at_32"] / max(b_["res_at_32"], 1e-300))
    print("VERDICT")
    if a3_plateau and b_conv:
        print(f"  PASS -- the existing per-node/per-fiber/global-12 preconditioners PLATEAU "
              f"(A3 res@32 = {a3['res_at_32']:.2e}, max nodal {a3['max_nodal']:.2e}),")
        print(f"         the FIBER-QUOTIENT inter-fiber coarse CONVERGES to 1e-8 in "
              f"{b_['iters_to_tol']} iters (max nodal {b_['max_nodal']:.2e}).")
        print(f"         residual@32 improved {speedup:.1e}x. The missing ingredient is exactly the "
              f"crosslink-weighted inter-fiber Galerkin coarse.")
    else:
        print(f"  INCONCLUSIVE -- A3 plateau={a3_plateau}, B converged={b_conv}; "
              f"adjust contrast/crosslink_fraction.")
    print()

    # Contrast robustness: the real crosslink is ~8e5 pN/um >> everything. Show the block-smoother
    # plateau WORSENS with contrast (its floor is set by the un-damped stiff-coupled inter-fiber
    # modes) while the fiber-quotient two-level stays converged in ~constant iterations -- the
    # signature of a contrast-independent coarse space, and the reason it will close the real stall.
    print("CONTRAST ROBUSTNESS  (A1 block-smoother floor vs B fiber-quotient; sweep k_crosslink/a)")
    print("-" * 96)
    print(f"{'contrast k/a':>14s} {'A1 res@64':>12s} {'A1 res@256':>12s} "
          f"{'B iters->1e-8':>14s} {'B res@64':>12s}")
    print("-" * 96)
    for contrast in (1e3, 1e4, 1e5, 1e6, 1e7):
        sp_ = ModelParams(**{**params.__dict__})
        sp_.k_backbone = sp_.k_bend = sp_.k_crosslink = float(contrast)
        sp_.regularizer = 1.0
        sn = build_network(sp_)
        Ks, es = assemble_stiffness(sn)
        bs = build_rhs(sn)
        blk = fiber_block_cholesky(Ks, sn)
        Ps, _ = fiber_quotient_prolongator(sn)
        fq = GalerkinCoarse(Ks, Ps)
        ra = pcg(Ks, bs, blk, tol=1e-8, maxiter=256)
        rb = pcg(Ks, bs, additive(blk, fq.apply), tol=1e-8, maxiter=256)
        bit = str(rb["iters_to_tol"]) if rb["iters_to_tol"] else ">256"
        print(f"{contrast:14.0e} {ra['history'][min(64, ra['history'].size - 1)]:12.2e} "
              f"{ra['history'][min(256, ra['history'].size - 1)]:12.2e} "
              f"{bit:>14s} {rb['history'][min(64, rb['history'].size - 1)]:12.2e}")
    print("-" * 96)
    print()

    # Sanity Gate.
    print("SANITY GATE")
    for line in sanity_gate(params):
        print("  " + line)
    print()
    print(f"total wall {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
