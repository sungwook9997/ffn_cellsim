"""IPC incremental-potential ENERGY kernels — the scalar E(x) the Newton line-search needs.

#1 Step 1 (DCM large-dt plan 2026-07-07). The projected-Newton IPC step minimizes the incremental
potential; its Armijo line-search needs a SCALAR energy whose gradient is EXACTLY the applied stiff
force (barrier + cortex-edge + turgor). If the energy's nearest-face pick differs from the force
kernel's, the line-search rejects valid steps — the subtle trap. These kernels mirror the force
kernels term-by-term:

  barrier  E_c = κ·area·b(d)            (OUTSIDE, d<d̂)   |  ½·κ·area·d²  (INSIDE, feasibilization)
           matches nearest_face_ipc_kernel: force = κ·area·(−b'(d)) outward  →  F = −∇E.  (adhesion is a
           separate soft/lagged force, NOT in this energy.)
  edge     E_e = ½·k_edge·(L−r0)²        matches _bond_accumulate: F = −k_edge·(L−r0)·û  = −∇E.
  turgor   E_c = (K_vol/2V0)(V_c−V0)² − dP0·(V_c−V0)   matches dcm_turgor_force: −dE/dV = dP0+K_vol(V0−V)/V0.

Gate 1 (in `fd_energy_gradient_check`): central-difference ∇E along random directions must equal
−F_stiff to ~1e-5 before the Newton driver (Step 2) may use this energy.
"""
from __future__ import annotations

import numpy as np
import warp as wp

from ffn_sim.dcm.dcm_contact_warp import closest_bary
from ffn_sim.dcm.dcm_contact_implicit_warp import _ipc_b


@wp.kernel
def ipc_barrier_energy_kernel(
    grid: wp.uint64,
    qpts: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    radius: wp.float32,
    rep: wp.float64,                          # κ (= derived rep), same as the force kernel
    d_hat: wp.float64,                         # activation gap d̂ (= c_rep)
    e_out: wp.array(dtype=wp.float64),         # (1,) scalar accumulator
):
    """Barrier energy on the SINGLE nearest other-cell face — the identical localization as
    ``nearest_face_ipc_kernel`` so ∇(this) = −(that force). Adhesion is excluded (soft/lagged)."""
    ni = wp.tid()
    z = wp.float64(0.0)
    c1 = cof[ni]
    if c1 < wp.int32(0):
        return
    p = pos[ni]
    best_d = wp.float64(1.0e300)
    best_i = wp.int32(-1)
    q = wp.hash_grid_query(grid, qpts[ni], radius)
    fj = wp.int32(0)
    while wp.hash_grid_query_next(q, fj):
        if fcell[fj] != c1:
            a = pos[faces[fj, 0]]
            b = pos[faces[fj, 1]]
            c = pos[faces[fj, 2]]
            bary = closest_bary(p, a, b, c)
            cpa = a * bary[0] + b * bary[1] + c * bary[2]
            d = wp.length(p - cpa)
            if d < best_d:
                best_d = d
                best_i = fj
    if best_i < wp.int32(0):
        return
    a = pos[faces[best_i, 0]]
    b = pos[faces[best_i, 1]]
    c = pos[faces[best_i, 2]]
    bary = closest_bary(p, a, b, c)
    cpa = a * bary[0] + b * bary[1] + c * bary[2]
    r_vec = p - cpa
    min_d = wp.length(r_vec)
    fnv = wp.cross(b - a, c - a)
    nrm = wp.length(fnv)
    area = wp.float64(0.5) * nrm
    if nrm <= z or min_d <= z:
        return
    sign = wp.dot(r_vec, fnv) / nrm
    e = z
    if sign > z:                                # OUTSIDE — barrier cushion
        if min_d < d_hat:
            e = rep * area * _ipc_b(min_d, d_hat)
    else:                                       # INSIDE — feasibilization ½·κ·area·d² (force = κ·area·d)
        e = wp.float64(0.5) * rep * area * min_d * min_d
    if e != z:
        wp.atomic_add(e_out, 0, e)


@wp.kernel
def edge_energy_kernel(
    pos: wp.array(dtype=wp.vec3d),
    edges: wp.array(dtype=wp.int32, ndim=2),
    k: wp.float64, r0: wp.array(dtype=wp.float64),
    e_out: wp.array(dtype=wp.float64),
):
    """½·k_edge·(L−r0)² per edge — matches ``_bond_accumulate`` (F = −k·(L−r0)·û)."""
    e = wp.tid()
    i = edges[e, 0]
    j = edges[e, 1]
    L = wp.length(pos[i] - pos[j])
    dL = L - r0[e]
    wp.atomic_add(e_out, 0, wp.float64(0.5) * k * dL * dL)


@wp.kernel
def turgor_energy_kernel(
    Vc: wp.array(dtype=wp.float64),            # per-cell current volume
    V0: wp.array(dtype=wp.float64),            # per-cell rest volume
    K_vol: wp.float64, dP0: wp.float64,
    e_out: wp.array(dtype=wp.float64),
):
    """(K_vol/2V0)(V_c−V0)² − dP0·(V_c−V0) per cell — matches dcm_turgor_force (−dE/dV = dP)."""
    c = wp.tid()
    v0 = V0[c]
    if v0 <= wp.float64(0.0):
        return
    u = Vc[c] - v0
    wp.atomic_add(e_out, 0, (K_vol / (wp.float64(2.0) * v0)) * u * u - dP0 * u)


@wp.kernel
def inertial_energy_kernel(
    x: wp.array(dtype=wp.vec3d),
    xn: wp.array(dtype=wp.vec3d),
    a: wp.float64,                             # a = γ_node/dt (implicit-Euler regulariser)
    cof: wp.array(dtype=wp.int32),
    e_out: wp.array(dtype=wp.float64),
):
    """½·a·|x−xn|² per LIVE node — the implicit-Euler inertial term of the incremental potential
    Φ(x)=½a|x−xn|²+U(x). ∇=a(x−xn), so ∇Φ = a(x−xn)−F = G (the Newton residual)."""
    i = wp.tid()
    if cof[i] < wp.int32(0):
        return
    d = x[i] - xn[i]
    wp.atomic_add(e_out, 0, wp.float64(0.5) * a * wp.dot(d, d))


@wp.kernel
def soft_linear_energy_kernel(
    x: wp.array(dtype=wp.vec3d),
    xn: wp.array(dtype=wp.vec3d),
    fsoft: wp.array(dtype=wp.vec3d),           # lagged soft force F_soft(xn), frozen over the step
    cof: wp.array(dtype=wp.int32),
    e_out: wp.array(dtype=wp.float64),
):
    """−F_soft(xn)·(x−xn) per LIVE node — the linear potential of the LAGGED soft/explicit drivers
    (cadherin, wetting, ECM, gravity). ∇=−F_soft(xn) (constant force), so its gradient reproduces
    the explicit RHS force exactly, keeping Φ's gradient == the full residual for the line-search."""
    i = wp.tid()
    if cof[i] < wp.int32(0):
        return
    d = x[i] - xn[i]
    wp.atomic_add(e_out, 0, -wp.dot(fsoft[i], d))
