"""Cytosim fiber inextensibility — constraint projector + reshape (Nédélec & Foethke 2007 §5.3).

NF2007 enforces fixed segment length as a HARD constraint, **not** an axial penalty spring (the
paper explicitly rejects the spring; ENGINE.md §2). For each fiber the ``p`` segment-length
constraints are

    C_k = (m_{k+1} − m_k)² − (L/p)² = 0,      k ∈ [0, p−1]                  [NF2007 §5.3, p11]

Two operations implement this, both grounded verbatim in the paper:

1. **Constraint projector** ``P = I − Jᵀ(JJᵀ)⁻¹J`` (symmetric, idempotent — an orthogonal
   projection). It removes the part of the force that would change any segment length, so the
   resolved motion stays tangent to the constraint manifold. With segment vectors ``g_k =
   m_{k+1}−m_k``, the Jacobian ``J_{kj} = ∂C_k/∂x_j`` has block ``−2 g_k`` at node ``k`` and
   ``+2 g_k`` at node ``k+1``; ``JJᵀ`` is then symmetric **tridiagonal** (consecutive constraints
   share one node) — cheap to invert per fiber (the paper notes ``JJᵀ`` is banded symmetric).

       (Jf)_k       = 2 g_k · (f_{k+1} − f_k)
       (JJᵀ)_{kk}   = 8 |g_k|²
       (JJᵀ)_{k,k+1}= −4 g_k · g_{k+1}
       Pf           = f − Jᵀ (JJᵀ)⁻¹ (Jf)

   The projector is the *force/velocity* side: it keeps the dynamics on the manifold to first
   order, but finite steps still drift, so —

2. **Reshape** restores the constraints EXACTLY after the points have moved (NF2007 §5.3, p11):
   "moving the points m₀,…,m_k in the direction of m_{k+1}−m_k and m_{k+1},…,m_p in the opposite
   direction, to restore |m_{k+1}−m_k| = L/p while conserving the center of gravity of the fiber."
   Done sequentially over segments; iterating converges (each pass corrects the residual coupling).

This is the inextensibility layer for the FF engine; the bending force lives in ``forces_warp``.
Pure-numpy, per-fiber (fibers are short ⇒ dense/tridiagonal solves are cheap); Warp kernelisation
is deferred to the cortex-scale assembly (Stage 6c) if profiling needs it.
"""

from __future__ import annotations

import numpy as np

from ffn_sim.ff.fiber_network import FiberNetwork


def _fiber_node_slice(net: FiberNetwork, f: int) -> slice:
    """Node index range owned by fiber ``f``."""
    return slice(int(net.fiber_offsets[f]), int(net.fiber_offsets[f + 1]))


def _project_one_fiber(pos: np.ndarray, force: np.ndarray, seg_rest: np.ndarray) -> np.ndarray:
    """Apply ``P = I − Jᵀ(JJᵀ)⁻¹J`` to ``force`` for one fiber.

    Args:
        pos: (n, 3) model points of the fiber [length].
        force: (n, 3) force on those points [force].
        seg_rest: (n−1,) segment rest lengths (unused by P itself — kept for signature symmetry
            with reshape; the projector depends only on the current geometry).

    Returns:
        (n, 3) projected force ``Pf``.
    """
    n = pos.shape[0]
    if n < 2:
        return force.copy()
    g = pos[1:] - pos[:-1]                       # (p, 3) segment vectors g_k
    p = g.shape[0]
    # Jf : (p,)  (Jf)_k = 2 g_k · (f_{k+1} − f_k)
    Jf = 2.0 * np.einsum("kd,kd->k", g, force[1:] - force[:-1])
    # JJᵀ : symmetric tridiagonal (p, p)
    diag = 8.0 * np.einsum("kd,kd->k", g, g)                 # 8|g_k|²
    if p > 1:
        off = -4.0 * np.einsum("kd,kd->k", g[:-1], g[1:])    # −4 g_k·g_{k+1}, length p−1
    else:
        off = np.zeros(0)
    JJt = np.diag(diag)
    if p > 1:
        JJt += np.diag(off, 1) + np.diag(off, -1)
    lam = np.linalg.solve(JJt, Jf)                           # (JJᵀ)⁻¹ Jf, size p
    # Jᵀλ : node i gets −2 g_i λ_i (constraint i, i<p) + 2 g_{i−1} λ_{i−1} (constraint i−1, i>0)
    JtLam = np.zeros_like(force)
    JtLam[:-1] += -2.0 * g * lam[:, None]                    # constraint k contributes to node k
    JtLam[1:] += 2.0 * g * lam[:, None]                      # constraint k contributes to node k+1
    return force - JtLam


def project_constraint_forces(net: FiberNetwork, force: np.ndarray) -> np.ndarray:
    """Project a force field onto the inextensibility constraint manifold, per fiber (NF2007 §5.3).

    Args:
        net: the fiber network (topology + current ``pos``).
        force: (N, 3) force on all model points [force].

    Returns:
        (N, 3) projected force ``Pf`` (the component that changes segment lengths removed).
    """
    force = np.asarray(force, dtype=np.float64).reshape(net.n_nodes, 3)
    out = force.copy()
    off = net.fiber_offsets
    seg_per_fiber = np.diff(off) - 1
    seg_cum = np.concatenate([[0], np.cumsum(seg_per_fiber)])
    for f in range(net.n_fibers):
        sl = _fiber_node_slice(net, f)
        s0, s1 = int(seg_cum[f]), int(seg_cum[f + 1])
        out[sl] = _project_one_fiber(net.pos[sl], force[sl], net.seg_rest[s0:s1])
    return out


def _reshape_one_fiber(pos: np.ndarray, seg_rest: np.ndarray, n_iter: int) -> np.ndarray:
    """Restore |m_{k+1}−m_k| = seg_rest_k exactly, conserving COG (NF2007 §5.3 reshape).

    Sequential over segments per the paper; ``n_iter`` passes converge the residual coupling
    (each correction perturbs neighbours). Center of gravity of the fiber is conserved per segment
    by splitting the correction between the two groups inversely to their sizes.
    """
    x = pos.copy()
    n = x.shape[0]
    if n < 2:
        return x
    p = n - 1
    for _ in range(n_iter):
        for k in range(p):
            g = x[k + 1] - x[k]
            d = float(np.linalg.norm(g))
            if d < 1e-300:
                continue
            u = g / d
            e = d - float(seg_rest[k])           # +e = too long → groups move together
            n_a = k + 1                           # group A = points 0..k
            n_b = p - k                           # group B = points k+1..p
            tot = n_a + n_b
            d_a = e * n_b / tot                   # A moves +d_a·u (toward B); COG: n_a d_a = n_b d_b
            d_b = e * n_a / tot                   # B moves −d_b·u
            x[: k + 1] += d_a * u
            x[k + 1:] -= d_b * u
    return x


def reshape(net: FiberNetwork, *, n_iter: int = 4, in_place: bool = False) -> np.ndarray:
    """Reshape every fiber to restore its segment rest lengths exactly (NF2007 §5.3).

    Args:
        net: the fiber network (mutated in place iff ``in_place``).
        n_iter: sequential-restore passes per fiber (more = tighter constraint satisfaction).
        in_place: if True, write the result back into ``net.pos`` and return it.

    Returns:
        (N, 3) reshaped positions.
    """
    out = net.pos.copy()
    off = net.fiber_offsets
    seg_per_fiber = np.diff(off) - 1
    seg_cum = np.concatenate([[0], np.cumsum(seg_per_fiber)])
    for f in range(net.n_fibers):
        sl = _fiber_node_slice(net, f)
        s0, s1 = int(seg_cum[f]), int(seg_cum[f + 1])
        out[sl] = _reshape_one_fiber(net.pos[sl], net.seg_rest[s0:s1], n_iter)
    if in_place:
        net.pos = out
    return out


def segment_lengths(net: FiberNetwork) -> np.ndarray:
    """Current segment lengths (S,) [length] — the constraint quantities, for diagnostics/tests."""
    s = net.segments
    if s.shape[0] == 0:
        return np.zeros(0)
    return np.linalg.norm(net.pos[s[:, 1]] - net.pos[s[:, 0]], axis=1)


def make_projected_force_fn(net: FiberNetwork, force_fn):
    """Wrap a force callable so it returns the constraint-projected force ``P·force_fn(pos)``.

    ``force_fn(pos_flat) → force_flat`` (e.g. from :func:`ff.forces_warp.make_bending_force_fn`).
    The returned callable updates ``net.pos`` from its argument (so P uses the current geometry),
    then projects. Use inside a relaxation/implicit loop together with periodic :func:`reshape`.
    """
    N = net.n_nodes

    def projected(pos_flat: np.ndarray) -> np.ndarray:
        pos = np.asarray(pos_flat, dtype=np.float64).reshape(N, 3)
        net.pos = pos
        F = np.asarray(force_fn(pos_flat), dtype=np.float64).reshape(N, 3)
        return project_constraint_forces(net, F).reshape(-1)

    return projected
