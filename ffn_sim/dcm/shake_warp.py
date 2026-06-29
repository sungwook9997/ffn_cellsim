"""Warp port of Matrix-SHAKE chain constraint projection (Phase B, M-SHAKE).

Ports the rigid uniform fast path of
``integrator.constrained_baoab.shake_project_chains`` (and its native CUDA twin
``shake_kernel.cu``): ONE THREAD PER CHAIN, per-chain serial Newton + tridiagonal
Thomas solve in fixed-length local vectors, so the float-op order matches the
Python einsum / ``_thomas_batched`` exactly → bit-identical (the native header
makes the same claim for its kernel).

Algorithm per chain (m bonds, m+1 beads, reference bond vectors d0 fixed):
  repeat (Newton):
    s_a   = min_image(pos[p_a] - pos[p_{a+1}])         current bond vectors
    g_a   = s_a·s_a − ℓ²                                constraint violation
    if max_a|g_a|/ℓ² ≤ tol: converged
    sd_a  = s_a·d0_a
    tridiagonal J Δλ = −g  (Thomas), J built from s·d0 cross-terms
    Δr[p_k] = M_{p_k} (λ_{k−1} d0_{k−1} − λ_k d0_k)     mobility-weighted projection
    λ_total += λ

Only the RIGID equality constraint (``compression_release=False``) is ported here
— the bit-identical core. The relaxed unilateral / Euler-buckling Gate-B mode is
a separate later piece. Orthorhombic box only (matches the cortex constraint path
and ``_min_image_orthorhombic``). ``wp.rint`` == numpy ``np.round`` (verified).
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()

# Max bonds per chain (matches native FFN_SHAKE_MMAX); +1 bead. The per-chain
# local vectors are sized to LMAX so a thread holds its whole chain in local mem.
MMAX = 64
LMAX = MMAX + 1
vecL = wp.types.vector(length=LMAX, dtype=wp.float64)


@wp.kernel
def shake_chain_kernel(
    pos: wp.array(dtype=wp.vec3d),                 # (N,) rw projected in place
    ref: wp.array(dtype=wp.vec3d),                 # (N,) ro start-of-step (d0)
    inv_mass: wp.array(dtype=wp.float64),          # (N,) ro mobility M_i
    chains: wp.array(dtype=wp.int32, ndim=2),      # (F, m+1) bead row indices
    lam_out: wp.array(dtype=wp.float64, ndim=2),   # (F, m) rw accumulated λ
    nonconv: wp.array(dtype=wp.int32),             # (1,) atomic nonconverged count
    m: wp.int32,
    rest_length: wp.float64,
    Lx: wp.float64,
    Ly: wp.float64,
    Lz: wp.float64,
    tol: wp.float64,
    max_iter: wp.int32,
):
    f = wp.tid()
    L2 = rest_length * rest_length

    d0x = vecL()
    d0y = vecL()
    d0z = vecL()
    Mv = vecL()

    for k in range(m + 1):
        Mv[k] = inv_mass[chains[f, k]]
    for a in range(m):
        ra = ref[chains[f, a]]
        rb = ref[chains[f, a + 1]]
        dx = ra[0] - rb[0]
        dy = ra[1] - rb[1]
        dz = ra[2] - rb[2]
        dx = dx - Lx * wp.rint(dx / Lx)
        dy = dy - Ly * wp.rint(dy / Ly)
        dz = dz - Lz * wp.rint(dz / Lz)
        d0x[a] = dx
        d0y[a] = dy
        d0z[a] = dz

    converged = wp.int32(0)
    for _it in range(max_iter):
        sx = vecL()
        sy = vecL()
        sz = vecL()
        g = vecL()
        sd = vecL()
        gmax_abs = wp.float64(0.0)
        for a in range(m):
            pa = pos[chains[f, a]]
            pb = pos[chains[f, a + 1]]
            dx = pa[0] - pb[0]
            dy = pa[1] - pb[1]
            dz = pa[2] - pb[2]
            dx = dx - Lx * wp.rint(dx / Lx)
            dy = dy - Ly * wp.rint(dy / Ly)
            dz = dz - Lz * wp.rint(dz / Lz)
            sx[a] = dx
            sy[a] = dy
            sz[a] = dz
            ga = dx * dx + dy * dy + dz * dz - L2
            g[a] = ga
            sd[a] = dx * d0x[a] + dy * d0y[a] + dz * d0z[a]
            ag = wp.abs(ga)
            if ag > gmax_abs:
                gmax_abs = ag
        if gmax_abs / L2 <= tol:
            converged = wp.int32(1)
            break

        diag = vecL()
        sub = vecL()
        sup = vecL()
        rhs = vecL()
        for a in range(m):
            diag[a] = wp.float64(-2.0) * (Mv[a] + Mv[a + 1]) * sd[a]
            rhs[a] = -g[a]
            sub[a] = wp.float64(0.0)
            sup[a] = wp.float64(0.0)
        if m > 1:
            for a in range(1, m):
                sub[a] = wp.float64(2.0) * Mv[a] * (
                    sx[a] * d0x[a - 1] + sy[a] * d0y[a - 1] + sz[a] * d0z[a - 1])
            for a in range(0, m - 1):
                sup[a] = wp.float64(2.0) * Mv[a + 1] * (
                    sx[a] * d0x[a + 1] + sy[a] * d0y[a + 1] + sz[a] * d0z[a + 1])

        # Thomas tridiagonal solve  J λ = rhs  (sub[0], sup[m-1] unused)
        cc = vecL()
        dp = vecL()
        lam = vecL()
        cc[0] = sup[0] / diag[0]
        dp[0] = rhs[0] / diag[0]
        for k in range(1, m):
            den = diag[k] - sub[k] * cc[k - 1]
            cc[k] = sup[k] / den
            dp[k] = (rhs[k] - sub[k] * dp[k - 1]) / den
        lam[m - 1] = dp[m - 1]
        for k in range(m - 2, -1, -1):
            lam[k] = dp[k] - cc[k] * lam[k + 1]

        for a in range(m):
            lam_out[f, a] = lam_out[f, a] + lam[a]

        # Δr[p_k] = M_k (λ_{k-1} d0_{k-1} − λ_k d0_k); subtract-then-add order
        # matches numpy's disp[:-1] -= ... ; disp[1:] += ... .
        for k in range(m + 1):
            ddx = wp.float64(0.0)
            ddy = wp.float64(0.0)
            ddz = wp.float64(0.0)
            if k <= m - 1:
                ddx = ddx - Mv[k] * lam[k] * d0x[k]
                ddy = ddy - Mv[k] * lam[k] * d0y[k]
                ddz = ddz - Mv[k] * lam[k] * d0z[k]
            if k >= 1:
                ddx = ddx + Mv[k] * lam[k - 1] * d0x[k - 1]
                ddy = ddy + Mv[k] * lam[k - 1] * d0y[k - 1]
                ddz = ddz + Mv[k] * lam[k - 1] * d0z[k - 1]
            idx = chains[f, k]
            p = pos[idx]
            pos[idx] = wp.vec3d(p[0] + ddx, p[1] + ddy, p[2] + ddz)

    if converged == 0:
        wp.atomic_add(nonconv, 0, 1)


def run_shake_warp(
    *,
    pred_pos: np.ndarray,   # (N, 3) post-predictor positions
    ref_pos: np.ndarray,    # (N, 3) start-of-step positions (d0 reference)
    chains: np.ndarray,     # (F, m+1) int32 bead row indices (uniform)
    inv_mass: np.ndarray,   # (N,) mobility
    rest_length: float,
    box_L: tuple[float, float, float],
    tol: float = 1.0e-10,
    max_iter: int = 100,
    device: str = "cpu",
) -> dict:
    """Project chains onto the rigid-bond manifold in Warp (one thread/chain).

    Returns the projected positions (N,3), accumulated per-bond Lagrange
    multipliers (F,m), and the nonconverged-chain count — matching
    ``shake_project_chains(..., return_lambdas=True)`` + the native nonconv out.
    """
    pred_pos = np.ascontiguousarray(pred_pos, dtype=np.float64)
    ref_pos = np.ascontiguousarray(ref_pos, dtype=np.float64)
    chains = np.ascontiguousarray(chains, dtype=np.int32)
    inv_mass = np.ascontiguousarray(inv_mass, dtype=np.float64)
    F, Np1 = chains.shape
    m = Np1 - 1
    if m > MMAX:
        raise ValueError(f"chain has {m} bonds > MMAX={MMAX}")
    Lx, Ly, Lz = box_L

    pos_d = wp.array(pred_pos, dtype=wp.vec3d, device=device)
    ref_d = wp.array(ref_pos, dtype=wp.vec3d, device=device)
    invm_d = wp.array(inv_mass, dtype=wp.float64, device=device)
    chains_d = wp.array(chains, dtype=wp.int32, device=device)
    lam_d = wp.zeros((F, m), dtype=wp.float64, device=device)
    nonconv_d = wp.zeros(1, dtype=wp.int32, device=device)

    wp.launch(
        shake_chain_kernel, dim=F,
        inputs=[pos_d, ref_d, invm_d, chains_d, lam_d, nonconv_d,
                wp.int32(m), wp.float64(rest_length),
                wp.float64(Lx), wp.float64(Ly), wp.float64(Lz),
                wp.float64(tol), wp.int32(max_iter)],
        device=device,
    )
    wp.synchronize_device(device)
    return {
        "pos": pos_d.numpy().astype(np.float64),
        "lam": lam_d.numpy().astype(np.float64),
        "nonconverged": int(nonconv_d.numpy()[0]),
    }
