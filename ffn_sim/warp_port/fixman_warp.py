"""Warp port of the Fixman metric pseudo-force (Phase B, Fixman piece).

Ports the uniform fast path of
``integrator.constrained_baoab.fixman_logdet_and_force`` (and its native CUDA twin
``fixman_kernel.cu``): ONE THREAD PER CHAIN builds the m×m tridiagonal mobility-
metric Gram matrix G, computes ln det G and the diagonal + first-super-diagonal
entries of G⁻¹ via the **tridiagonal continuant (θ/φ) recurrences**, and assembles
the analytic ∇ln det G force.

Per chain (m bonds, m+1 beads; bond vectors b_a = r_{a+1} − r_a, mobility M):
    G_aa     =  4 |b_a|² (M_a + M_{a+1})
    G_{a,a+1}= −4 M_{a+1} (b_a · b_{a+1})
    U_F      =  FIXMAN_SIGN · ½ kT · ln det G   (FIXMAN_SIGN = +1)
    F_{p_a}  = +∇  ,  F_{p_{a+1}} = −∇  with ∇ = FIXMAN_SIGN·½kT·(d ln det G / d b_a)

Continuant (0-based; d_a=G_aa, e_a=G_{a,a+1}):
    θ_0=1, θ_1=d_0, θ_k = d_{k-1}θ_{k-1} − e_{k-2}²θ_{k-2}            det G = θ_m
    φ_m=1, φ_{m-1}=d_{m-1}, φ_k = d_kφ_{k+1} − e_k²φ_{k+2}
    (G⁻¹)_aa     = θ_a φ_{a+1} / det
    (G⁻¹)_{a,a+1}= −e_a θ_a φ_{a+2} / det

Parity is TOLERANCE-based (not bitwise), exactly as the native plugin documents:
the Python reference uses LAPACK ``slogdet``/``inv`` while this uses the continuant
recurrence — they differ in the last ULPs. Gate = the native plugin's own contract
(force rel < 1e-7, U_F rel < 1e-9, no non-positive metric). The force is what feeds
the dynamics.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()

FIXMAN_SIGN = 1.0  # matches constrained_baoab.FIXMAN_SIGN
MMAX = 64
LMAX = MMAX + 1
vecL = wp.types.vector(length=LMAX, dtype=wp.float64)


@wp.kernel
def fixman_chain_kernel(
    pos: wp.array(dtype=wp.vec3d),                 # (N,) ro
    inv_gamma: wp.array(dtype=wp.float64),         # (N,) ro per-bead mobility
    chains: wp.array(dtype=wp.int32, ndim=2),      # (F, m+1)
    force: wp.array(dtype=wp.vec3d),               # (N,) out (zeroed; chain beads set)
    logdet: wp.array(dtype=wp.float64),            # (F,) out ln det G
    bad_sign: wp.array(dtype=wp.int32),            # (F,) out 1 if det G <= 0
    m: wp.int32,
    half_kT: wp.float64,
    sign: wp.float64,
    Lx: wp.float64,
    Ly: wp.float64,
    Lz: wp.float64,
):
    f = wp.tid()

    bx = vecL()
    by = vecL()
    bz = vecL()
    Mv = vecL()
    for k in range(m + 1):
        Mv[k] = inv_gamma[chains[f, k]]
    for a in range(m):
        ra = pos[chains[f, a]]
        rb = pos[chains[f, a + 1]]
        dx = rb[0] - ra[0]
        dy = rb[1] - ra[1]
        dz = rb[2] - ra[2]
        dx = dx - Lx * wp.rint(dx / Lx)
        dy = dy - Ly * wp.rint(dy / Ly)
        dz = dz - Lz * wp.rint(dz / Lz)
        bx[a] = dx
        by[a] = dy
        bz[a] = dz

    diag = vecL()
    off = vecL()
    for a in range(m):
        b2 = bx[a] * bx[a] + by[a] * by[a] + bz[a] * bz[a]
        diag[a] = wp.float64(4.0) * b2 * (Mv[a] + Mv[a + 1])
    for a in range(m - 1):
        bdot = bx[a] * bx[a + 1] + by[a] * by[a + 1] + bz[a] * bz[a + 1]
        off[a] = wp.float64(-4.0) * Mv[a + 1] * bdot

    # continuant: leading (θ) and trailing (φ) minor determinants
    theta = vecL()
    phi = vecL()
    theta[0] = wp.float64(1.0)
    theta[1] = diag[0]
    for k in range(2, m + 1):
        theta[k] = diag[k - 1] * theta[k - 1] - off[k - 2] * off[k - 2] * theta[k - 2]
    phi[m] = wp.float64(1.0)
    phi[m - 1] = diag[m - 1]
    for k in range(m - 2, -1, -1):
        phi[k] = diag[k] * phi[k + 1] - off[k] * off[k] * phi[k + 2]

    det = theta[m]
    if det <= wp.float64(0.0):
        bad_sign[f] = wp.int32(1)
        logdet[f] = wp.float64(0.0)
        return

    logdet[f] = wp.log(det)
    inv_det = wp.float64(1.0) / det

    gdiag = vecL()
    gsuper = vecL()
    for a in range(m):
        gdiag[a] = theta[a] * phi[a + 1] * inv_det
    for a in range(m - 1):
        gsuper[a] = -off[a] * theta[a] * phi[a + 2] * inv_det

    gradx = vecL()
    grady = vecL()
    gradz = vecL()
    for a in range(m):
        # diagonal term: Ginv[a,a]·8·(M_a+M_{a+1})·b_a
        cf = gdiag[a] * wp.float64(8.0) * (Mv[a] + Mv[a + 1])
        cx = cf * bx[a]
        cy = cf * by[a]
        cz = cf * bz[a]
        # coupling to previous bond: 2·Ginv[a-1,a]·(-4 M_a)·b_{a-1}
        if a >= 1:
            cp = wp.float64(2.0) * gsuper[a - 1] * (wp.float64(-4.0) * Mv[a])
            cx = cx + cp * bx[a - 1]
            cy = cy + cp * by[a - 1]
            cz = cz + cp * bz[a - 1]
        # coupling to next bond: 2·Ginv[a,a+1]·(-4 M_{a+1})·b_{a+1}
        if a <= m - 2:
            cn = wp.float64(2.0) * gsuper[a] * (wp.float64(-4.0) * Mv[a + 1])
            cx = cx + cn * bx[a + 1]
            cy = cy + cn * by[a + 1]
            cz = cz + cn * bz[a + 1]
        gradx[a] = sign * half_kT * cx
        grady[a] = sign * half_kT * cy
        gradz[a] = sign * half_kT * cz

    # F_{p_k} = +grad_k (k≤m-1)  −grad_{k-1} (k≥1); beads disjoint across chains
    for k in range(m + 1):
        fx = wp.float64(0.0)
        fy = wp.float64(0.0)
        fz = wp.float64(0.0)
        if k <= m - 1:
            fx = fx + gradx[k]
            fy = fy + grady[k]
            fz = fz + gradz[k]
        if k >= 1:
            fx = fx - gradx[k - 1]
            fy = fy - grady[k - 1]
            fz = fz - gradz[k - 1]
        force[chains[f, k]] = wp.vec3d(fx, fy, fz)


def run_fixman_warp(
    *,
    pos: np.ndarray,        # (N, 3)
    chains: np.ndarray,     # (F, m+1) int32
    inv_gamma: np.ndarray,  # (N,)
    kT: float,
    box_L: tuple[float, float, float],
    device: str = "cpu",
) -> dict:
    """Fixman force + per-chain ln det G in Warp (one thread/chain).

    Returns ``force`` (N,3), ``U_F`` (scalar = FIXMAN_SIGN·½kT·Σ ln det G),
    per-chain ``logdet`` (F,), and ``bad_sign`` count — matching
    ``fixman_logdet_and_force`` + the native nonpositive-metric flag.
    """
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    chains = np.ascontiguousarray(chains, dtype=np.int32)
    inv_gamma = np.ascontiguousarray(inv_gamma, dtype=np.float64)
    N = pos.shape[0]
    F, Np1 = chains.shape
    m = Np1 - 1
    if m < 2:
        raise ValueError(f"Fixman uniform fast path requires m>=2 bonds, got {m}")
    if m > MMAX:
        raise ValueError(f"chain has {m} bonds > MMAX={MMAX}")
    Lx, Ly, Lz = box_L

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    invg_d = wp.array(inv_gamma, dtype=wp.float64, device=device)
    chains_d = wp.array(chains, dtype=wp.int32, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    logdet_d = wp.zeros(F, dtype=wp.float64, device=device)
    bad_d = wp.zeros(F, dtype=wp.int32, device=device)

    wp.launch(
        fixman_chain_kernel, dim=F,
        inputs=[pos_d, invg_d, chains_d, force_d, logdet_d, bad_d,
                wp.int32(m), wp.float64(0.5 * kT), wp.float64(FIXMAN_SIGN),
                wp.float64(Lx), wp.float64(Ly), wp.float64(Lz)],
        device=device,
    )
    wp.synchronize_device(device)
    logdet = logdet_d.numpy().astype(np.float64)
    bad = int(bad_d.numpy().sum())
    U_F = FIXMAN_SIGN * 0.5 * kT * float(logdet.sum())
    return {
        "force": force_d.numpy().astype(np.float64),
        "U_F": U_F,
        "logdet": logdet,
        "bad_sign": bad,
    }
