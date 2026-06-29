"""Cytosim fiber mechanics — Warp force kernels (grounded: Nédélec & Foethke 2007, NJP 9:427).

**Bending** (Stage 6b, the load-bearing kernel). For each interior model-point ``m_i`` the
consecutive triplet {m_{i-1}, m_i, m_{i+1}} carries the force triplet {−F, +2F, −F} with

    F = α (m_{i-1} − 2 m_i + m_{i+1}),    α = κ (p/L)³ = κ / seg³            [NF2007 p9, fig.6]

where κ = bending modulus = k_B·T·L_p = E·I, p = #segments, L = fiber length, seg = L/p. The
discrete bending energy whose negative gradient this force is:

    E_bend = (α/2) Σ_triples |m_{i-1} − 2 m_i + m_{i+1}|²                     [J]

(the assembled stiffness is the symmetric banded matrix α·E with interior 4th-difference stencil
[−1, 4, −6, 4, −1] and reduced end rows; row/col sums = 0 → internal torque only, fig.6).
Validation anchor: the buckling threshold equals Euler's π²κ/L²  [NF2007 p9].

Non-extensibility is **NOT** a spring here — NF2007 enforces fixed segment length as a hard
constraint + "reshape" projection (p3, §5.3), explicitly rejecting an axial penalty spring; that
constraint/implicit layer is a separate step (wired through ``dcm.dcm_warp_implicit``). This module
is the bending force only, kernel-pattern-matched to ``dcm/network_warp.py`` (one thread per
element, ``wp.atomic_add`` accumulation, float64).
"""

from __future__ import annotations

import numpy as np
import warp as wp

from ffn_sim.ff.fiber_network import FiberNetwork

wp.init()


@wp.kernel
def cytosim_bending_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) model-point positions [m]
    triples: wp.array(dtype=wp.int32, ndim=2),  # (T, 3) consecutive node indices (i-1, i, i+1)
    alpha: wp.array(dtype=wp.float64),          # (T,) α = κ/seg³ of the triple's fiber [N/m]
    force: wp.array(dtype=wp.vec3d),            # (N,) out (zeroed; atomic accumulate)
):
    t = wp.tid()
    a = triples[t, 0]
    b = triples[t, 1]
    c = triples[t, 2]
    # discrete curvature (second difference) d = m_{i-1} − 2 m_i + m_{i+1}
    d = pos[a] - wp.float64(2.0) * pos[b] + pos[c]
    f = alpha[t] * d                            # F = α d
    wp.atomic_add(force, a, -f)                 # triplet {−F, +2F, −F}
    wp.atomic_add(force, b, wp.float64(2.0) * f)
    wp.atomic_add(force, c, -f)


def _per_triple_alpha(net: FiberNetwork) -> np.ndarray:
    """α = κ/seg³ for each bending triple (seg = mean segment length of the triple's fiber).

    NF2007: α = κ (p/L)³ with seg = L/p, so α = κ/seg³.
    """
    off = net.fiber_offsets
    seg_per_fiber = np.diff(off) - 1                       # p_f = #segments per fiber
    seg_off = np.concatenate([[0], np.cumsum(seg_per_fiber)])
    seg_mean = np.ones(net.n_fibers, dtype=np.float64)
    for f in range(net.n_fibers):
        s0, s1 = int(seg_off[f]), int(seg_off[f + 1])
        if s1 > s0:
            seg_mean[f] = float(net.seg_rest[s0:s1].mean())
    alpha_fiber = net.kappa / seg_mean**3                  # κ/seg³  [N/m]
    centers = net.bend_triples[:, 1]                       # the m_i of each triple
    fib_of_node = np.searchsorted(off, centers, side="right") - 1
    return alpha_fiber[fib_of_node].astype(np.float64)


def bending_force(net: FiberNetwork, *, device: str = "cpu") -> np.ndarray:
    """Cytosim bending force on the model points, (N, 3) [N] (NF2007 p9).

    Reads ``net.pos``. For an integration loop use :func:`make_bending_force_fn` instead.
    """
    return _launch_bending(net.pos, net, _per_triple_alpha(net), device=device)


def make_bending_force_fn(net: FiberNetwork, *, device: str = "cpu"):
    """Return ``force_fn(pos_flat) → force_flat`` for ``dcm.dcm_warp_implicit`` / a relaxation loop.

    Topology + α are uploaded once; each call evaluates the bending force at the given positions.
    """
    triples = np.ascontiguousarray(net.bend_triples, dtype=np.int32)
    alpha = np.ascontiguousarray(_per_triple_alpha(net), dtype=np.float64)
    triples_d = wp.array(triples, dtype=wp.int32, device=device)
    alpha_d = wp.array(alpha, dtype=wp.float64, device=device)
    N = net.n_nodes
    T = triples.shape[0]

    def force_fn(pos_flat: np.ndarray) -> np.ndarray:
        pos = np.ascontiguousarray(np.asarray(pos_flat, dtype=np.float64).reshape(N, 3))
        pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
        force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
        if T > 0:
            wp.launch(cytosim_bending_kernel, dim=T,
                      inputs=[pos_d, triples_d, alpha_d, force_d], device=device)
            wp.synchronize_device(device)
        return force_d.numpy().astype(np.float64).reshape(-1)

    return force_fn


def _launch_bending(pos: np.ndarray, net: FiberNetwork, alpha: np.ndarray,
                    *, device: str) -> np.ndarray:
    triples = np.ascontiguousarray(net.bend_triples, dtype=np.int32)
    N = net.n_nodes
    T = triples.shape[0]
    pos_d = wp.array(np.ascontiguousarray(pos, dtype=np.float64), dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    if T > 0:
        triples_d = wp.array(triples, dtype=wp.int32, device=device)
        alpha_d = wp.array(np.ascontiguousarray(alpha, dtype=np.float64), dtype=wp.float64,
                           device=device)
        wp.launch(cytosim_bending_kernel, dim=T,
                  inputs=[pos_d, triples_d, alpha_d, force_d], device=device)
        wp.synchronize_device(device)
    return force_d.numpy().astype(np.float64)


def bending_energy(net: FiberNetwork) -> float:
    """Discrete bending energy E = (α/2) Σ |m_{i-1} − 2 m_i + m_{i+1}|²  [J] (NF2007)."""
    tri = net.bend_triples
    if tri.shape[0] == 0:
        return 0.0
    alpha = _per_triple_alpha(net)
    d = net.pos[tri[:, 0]] - 2.0 * net.pos[tri[:, 1]] + net.pos[tri[:, 2]]
    return float(0.5 * np.sum(alpha * np.sum(d * d, axis=1)))
