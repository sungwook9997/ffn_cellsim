"""FF-specific implicit overdamped solver (NF2007 Eq 2) — FF-scale, SEPARATE from the DCM solver.

Why separate (PI 2026-07-06): the DCM ``dcm_warp_implicit.implicit_overdamped_step`` is a matrix-free CG on a
finite-difference Jacobian-vector product. That is VERIFIED to blow up at FF force/position scales (positions
~µm, elastic forces ~1e5-1e6 pN): the FD JVP is unreliable, the CG bails at 1 iteration, and the step degenerates
to explicit → diverges at large dt (tested: ε from 1e-9 to 1e-2 all identical blow-up). It cannot be salvaged by
a tolerance/ε calibration.

This solver uses the NF2007 approach instead — the ONE the FF engine is supposed to have (ENGINE.md §2, the
"key contribution" / Stage-6b piece): assemble the **ANALYTIC** elastic stiffness matrix and solve a real sparse
linear system. No finite-difference Jacobian.

Semi-implicit (IMEX) overdamped Euler:  γ(x_{n+1}−x_n)/dt = F(x_{n+1}) ≈ F(x_n) + A·Δx,  A = ∂F_elastic/∂x
    ⇒ (γ/dt·I − A)·Δx = F_total(x_n),   K ≔ −A (elastic stiffness, PSD)  ⇒  (γ/dt·I + K)·Δx = F_total(x_n)

- **K assembled analytically**: bending = α·DᵀD (banded 4th-difference `[1,−4,6,−4,1]` per fiber, per
  coordinate — the discrete-bending Hessian); crosslinks = k·ûûᵀ rank-1 spring stiffness per bond.
- **M = γ/dt·I + K is SPD**, assembled at the reference config → **prefactored ONCE** (sparse LU); each step is
  one cheap back-substitution with a fresh RHS. Unconditionally stable ⇒ dt bounded by ACCURACY not the stiff
  α-actinin CFL — the NF2007 10⁴×-class speed-up that makes minutes-scale spreading/migration reachable.
- **Active/soft forces** (turgor, membrane, nucleus, substrate, clutch, protrusion, gravity, volume) enter
  ``F_total`` EXPLICITLY (they are soft; keeping them out of K keeps M SPD and the factorization reusable).

Inextensibility is handled by the caller's periodic reshape (NF2007 §5.3), as in the explicit path.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import factorized


def assemble_elastic_stiffness(net, alpha_per_triple: np.ndarray, xl_ij: np.ndarray,
                               k_xl: np.ndarray, pos: np.ndarray) -> sp.csr_matrix:
    """Analytic elastic stiffness K = K_bend + K_xl (3N×3N, PSD). K_bend = α·DᵀD (bending Hessian, per fiber,
    per coordinate); K_xl = k·ûûᵀ per crosslink (leading-order spring stiffness at the reference geometry)."""
    N = int(net.n_nodes)
    n3 = 3 * N
    tri = np.ascontiguousarray(net.bend_triples, np.int64)     # (T,3): a,b,c per interior triple
    al = np.ascontiguousarray(alpha_per_triple, np.float64)
    rows, cols, vals = [], [], []

    # ---- bending: per triple (a,b,c), per coordinate d, add α·(2nd-diff)ᵀ(2nd-diff) ----
    # force pattern {−F,+2F,−F} with F=α(m_a−2m_b+m_c) ⇒ K rows: a:[α,−2α,α] b:[−2α,4α,−2α] c:[α,−2α,α]
    blk = np.array([[1.0, -2.0, 1.0], [-2.0, 4.0, -2.0], [1.0, -2.0, 1.0]])
    for t in range(tri.shape[0]):
        a, b, c = tri[t]
        aa = al[t]
        nodes = (a, b, c)
        for r in range(3):
            for cc in range(3):
                v = aa * blk[r, cc]
                if v == 0.0:
                    continue
                for d in range(3):
                    rows.append(3 * nodes[r] + d); cols.append(3 * nodes[cc] + d); vals.append(v)

    # ---- crosslinks: per bond (i,j), K block [[P,−P],[−P,P]] with P = k·ûûᵀ ----
    xl = np.ascontiguousarray(xl_ij, np.int64)
    kx = np.ascontiguousarray(k_xl, np.float64)
    for e in range(xl.shape[0]):
        i, j = xl[e]
        rvec = pos[j] - pos[i]
        L = float(np.linalg.norm(rvec))
        if L < 1e-9:
            continue
        u = rvec / L
        P = kx[e] * np.outer(u, u)                             # 3×3 rank-1
        for d1 in range(3):
            for d2 in range(3):
                p = P[d1, d2]
                if p == 0.0:
                    continue
                ii, jj = 3 * i + d1, 3 * j + d1
                ci, cj = 3 * i + d2, 3 * j + d2
                rows += [ii, ii, jj, jj]
                cols += [ci, cj, cj, ci]
                vals += [p, -p, p, -p]

    K = sp.coo_matrix((vals, (rows, cols)), shape=(n3, n3)).tocsr()
    return K


class FFImplicitStepper:
    """Prefactored FF implicit overdamped stepper. Build once (assembles K at the reference config, forms
    M = γ/dt·I + K, LU-factorises it); then call :meth:`step` each timestep with the current TOTAL force."""

    def __init__(self, net, alpha_per_triple, xl_ij, k_xl, pos, *, gamma: float, dt: float):
        self.N = int(net.n_nodes)
        self.gamma = float(gamma)
        self.dt = float(dt)
        K = assemble_elastic_stiffness(net, alpha_per_triple, xl_ij, k_xl, pos)
        M = (self.gamma / self.dt) * sp.identity(3 * self.N, format="csr") + K
        self.K = K
        self.M = M.tocsc()
        self._solve = factorized(self.M)                       # LU prefactor (reused every step)

    def step(self, x_flat: np.ndarray, F_total_flat: np.ndarray) -> np.ndarray:
        """One implicit overdamped step: solve (γ/dt·I + K)·Δx = F_total, return x + Δx. Unconditionally stable
        in the elastic modes → dt can be large. ``F_total_flat`` is the full force (elastic + active) at x."""
        dx = self._solve(np.ascontiguousarray(F_total_flat, np.float64))
        return x_flat + dx



# ------------------------------------------------------------------------------------------------------
# CURRENT-config solver (the full NF2007 stability: A assembled at x_n each step ⇒ NSD ⇒ unconditionally
# stable). Vectorised assembly (fast enough to re-do each step) + CG on the ASSEMBLED SPD matrix (reliable,
# unlike the FD Jacobian-vector product). Optional Newton iterations + rigid-mode-safe.
# ------------------------------------------------------------------------------------------------------

_BEND_BLOCK = np.array([[1.0, -2.0, 1.0], [-2.0, 4.0, -2.0], [1.0, -2.0, 1.0]])   # α·DᵀD pattern


def assemble_K_current(pos: np.ndarray, bend_triples: np.ndarray, alpha: np.ndarray,
                       xl_ij: np.ndarray, k_xl: np.ndarray, N: int) -> sp.csr_matrix:
    """Analytic elastic stiffness K = −∂F_elastic/∂x at the CURRENT config ``pos`` (vectorised). Bending =
    α·DᵀD (config-independent); crosslinks = k·ûûᵀ from the current bond directions. PSD."""
    n3 = 3 * N
    # ---- bending (vectorised): T triples × 9 (r,c) × 3 coords ----
    tri = np.ascontiguousarray(bend_triples, np.int64)
    al = np.ascontiguousarray(alpha, np.float64)
    T = tri.shape[0]
    ridx, cidx = np.meshgrid(np.arange(3), np.arange(3), indexing="ij")
    ridx = ridx.ravel(); cidx = cidx.ravel(); bval = _BEND_BLOCK.ravel()               # (9,)
    nr = tri[:, ridx]; nc = tri[:, cidx]                                                 # (T,9)
    v9 = al[:, None] * bval[None, :]                                                     # (T,9)
    coord = np.arange(3)
    rows_b = (3 * nr[:, :, None] + coord[None, None, :]).ravel()
    cols_b = (3 * nc[:, :, None] + coord[None, None, :]).ravel()
    vals_b = np.broadcast_to(v9[:, :, None], (T, 9, 3)).ravel()
    # ---- crosslinks (vectorised): E bonds × 9 (d1,d2) × 4 blocks ----
    xl = np.ascontiguousarray(xl_ij, np.int64)
    i = xl[:, 0]; j = xl[:, 1]
    r = pos[j] - pos[i]
    L = np.linalg.norm(r, axis=1)
    good = L > 1e-9
    i, j, r, L = i[good], j[good], r[good], L[good]
    kx = np.ascontiguousarray(k_xl, np.float64)[good]
    u = r / L[:, None]                                                                   # (E,3)
    P = kx[:, None, None] * u[:, :, None] * u[:, None, :]                                # (E,3,3) = k·ûûᵀ
    d1, d2 = np.meshgrid(np.arange(3), np.arange(3), indexing="ij")
    d1 = d1.ravel(); d2 = d2.ravel()                                                     # (9,)
    Pf = P.reshape(-1, 9)                                                                # (E,9)
    ri, rj = 3 * i[:, None] + d1[None, :], 3 * j[:, None] + d1[None, :]                  # (E,9)
    ci, cj = 3 * i[:, None] + d2[None, :], 3 * j[:, None] + d2[None, :]
    rows_x = np.concatenate([ri.ravel(), ri.ravel(), rj.ravel(), rj.ravel()])
    cols_x = np.concatenate([ci.ravel(), cj.ravel(), cj.ravel(), ci.ravel()])
    vals_x = np.concatenate([Pf.ravel(), -Pf.ravel(), Pf.ravel(), -Pf.ravel()])
    rows = np.concatenate([rows_b, rows_x]); cols = np.concatenate([cols_b, cols_x])
    vals = np.concatenate([vals_b, vals_x])
    return sp.coo_matrix((vals, (rows, cols)), shape=(n3, n3)).tocsr()


def implicit_step_current(x_flat, force_fn, bend_triples, alpha, xl_ij, k_xl, *, gamma, dt,
                          n_newton=1, cg_tol=1e-6, cg_maxiter=300, vol_g=None, k_vol=0.0):
    """One NF2007 implicit overdamped step with the CURRENT-config Jacobian (unconditionally stable). Solves
    (γ/dt·I + K(x))·Δx = F_total(x) via CG on the assembled SPD matrix — no finite-difference JVP. ``force_fn``
    returns the FULL force (elastic + active) at a given x. Newton (n_newton>1) re-linearises for big steps.

    ``vol_g`` (3N,) + ``k_vol``>0: also treat the OSMOTIC VOLUME constraint implicitly as the rank-1 stiffness
    ``k_vol·g·gᵀ`` (g = ∂V/∂x). Without this the stiff osmotic (Π_in0≈5e5) caps the stable dt; with it, the
    implicit step is stable at large dt → the real wall-fast, minutes-scale speed-up. The operator stays SPD
    (k_vol>0), so CG converges."""
    from scipy.sparse.linalg import cg, LinearOperator
    N = x_flat.size // 3
    n3 = 3 * N
    a = gamma / dt
    x = np.ascontiguousarray(x_flat, np.float64).copy()
    x0 = x.copy()
    info = {"cg_iters": 0, "newton": 0}
    for _ in range(n_newton):
        pos = x.reshape(N, 3)
        K = assemble_K_current(pos, bend_triples, alpha, xl_ij, k_xl, N)
        M = (a * sp.identity(n3, format="csr") + K).tocsr()
        F = np.asarray(force_fn(x), dtype=np.float64).reshape(-1)
        rhs = F - a * (x - x0)                              # residual RHS (0 net at x0 for n_newton=1)
        if vol_g is not None and k_vol > 0.0:              # M + k_vol·g·gᵀ  (osmotic volume, implicit rank-1)
            g = np.ascontiguousarray(vol_g, np.float64).reshape(-1)
            op = LinearOperator((n3, n3), matvec=lambda v: M @ v + k_vol * (g @ v) * g, dtype=np.float64)
        else:
            op = M
        it = [0]
        dx, _ = cg(op, rhs, rtol=cg_tol, maxiter=cg_maxiter, callback=lambda *_a: it.__setitem__(0, it[0] + 1))
        x = x + dx
        info["cg_iters"] += it[0]; info["newton"] += 1
    return x, info



# ------------------------------------------------------------------------------------------------------
# GPU-RESIDENT version (cupy sparse + CG) — the PI's GPU-only mandate. Same math as the host path, but K is
# assembled and CG-solved ON THE DEVICE (cupy), and the force / positions stay device-resident (Warp ↔ cupy via
# the CUDA array interface, no host round-trip). This is what makes native (≈495k nodes) + large dt feasible.
# ------------------------------------------------------------------------------------------------------

def assemble_K_current_cupy(pos, bend_triples, alpha, xl_ij, k_xl, N, diag_extra=None):
    """Analytic elastic stiffness K (bending 4th-diff + crosslink k·ûûᵀ) assembled ON THE GPU from cupy arrays.
    ``pos`` (N,3), ``bend_triples`` (T,3), ``xl_ij`` (E,2), ``alpha`` (T,), ``k_xl`` (E,) are all cupy. Returns a
    cupy CSR (3N×3N). ``diag_extra`` (3N,) optionally adds a per-DOF diagonal stiffness — used to make the CLUTCH
    (k_int·I on bound basal nodes) and SUBSTRATE (k_plane on basal z) IMPLICIT too, so dt is no longer capped by
    their explicit CFL (it rises to the next-softest force ≈ the nucleus)."""
    import cupy as cp
    import cupyx.scipy.sparse as csp
    n3 = 3 * N
    tri = bend_triples; al = alpha; T = tri.shape[0]
    ridx, cidx = cp.meshgrid(cp.arange(3), cp.arange(3), indexing="ij")
    ridx = ridx.ravel(); cidx = cidx.ravel()
    bval = cp.asarray(_BEND_BLOCK).ravel()
    nr = tri[:, ridx]; nc = tri[:, cidx]
    v9 = al[:, None] * bval[None, :]
    coord = cp.arange(3)
    rows_b = (3 * nr[:, :, None] + coord[None, None, :]).ravel()
    cols_b = (3 * nc[:, :, None] + coord[None, None, :]).ravel()
    vals_b = cp.broadcast_to(v9[:, :, None], (T, 9, 3)).ravel()
    i = xl_ij[:, 0]; j = xl_ij[:, 1]
    r = pos[j] - pos[i]; L = cp.linalg.norm(r, axis=1)
    good = L > 1e-9
    i = i[good]; j = j[good]; u = (r[good] / L[good, None]); kx = k_xl[good]
    P = kx[:, None, None] * u[:, :, None] * u[:, None, :]
    d1, d2 = cp.meshgrid(cp.arange(3), cp.arange(3), indexing="ij"); d1 = d1.ravel(); d2 = d2.ravel()
    Pf = P.reshape(-1, 9)
    ri = 3 * i[:, None] + d1[None, :]; rj = 3 * j[:, None] + d1[None, :]
    ci = 3 * i[:, None] + d2[None, :]; cj = 3 * j[:, None] + d2[None, :]
    rows_x = cp.concatenate([ri.ravel(), ri.ravel(), rj.ravel(), rj.ravel()])
    cols_x = cp.concatenate([ci.ravel(), cj.ravel(), cj.ravel(), ci.ravel()])
    vals_x = cp.concatenate([Pf.ravel(), -Pf.ravel(), Pf.ravel(), -Pf.ravel()])
    rows = cp.concatenate([rows_b, rows_x]); cols = cp.concatenate([cols_b, cols_x])
    vals = cp.concatenate([vals_b, vals_x])
    K = csp.coo_matrix((vals, (rows, cols)), shape=(n3, n3)).tocsr()
    if diag_extra is not None:
        K = K + csp.diags(cp.asarray(diag_extra, cp.float64), format="csr")
    return K


def ff_implicit_step_gpu(x, force_fn, bend_triples, alpha, xl_ij, k_xl, *, gamma, dt,
                         vol_g=None, k_vol=0.0, diag_extra=None, cg_tol=1e-6, cg_maxiter=400):
    """One GPU-resident NF2007 implicit overdamped step. All arrays cupy, device-resident. ``force_fn(x_cp)``
    returns the full force as a cupy (3N,) array (Warp kernels → cupy view, no host). Solves
    (γ/dt·I + K(x) + k_vol·g·gᵀ)·Δx = F(x) with cupy CG. Returns x+Δx (cupy).

    Un-preconditioned CG is used deliberately: a diagonal (Jacobi) preconditioner was measured to INFLATE
    iterations here (K is dominated by the rank-1 crosslink blocks k·ûûᵀ with k_xl≈4.6e5, which are not
    diagonally dominant, so 1/diag mis-scales them — Thread-C R-CG, 2026-07-07). Merging the MT aster costs a
    bounded ~1.5× iters (native 486k: 113→171, well under maxiter); a block-Jacobi / IC(0) preconditioner is the
    flagged follow-up if native CG throughput becomes a production bottleneck."""
    import cupy as cp
    import cupyx.scipy.sparse as csp
    from cupyx.scipy.sparse.linalg import cg, LinearOperator
    N = x.size // 3; n3 = 3 * N; a = gamma / dt
    K = assemble_K_current_cupy(x.reshape(N, 3), bend_triples, alpha, xl_ij, k_xl, N, diag_extra=diag_extra)
    M = (a * csp.identity(n3, format="csr", dtype=cp.float64) + K).tocsr()
    F = force_fn(x)
    if vol_g is not None and k_vol > 0.0:
        g = vol_g
        op = LinearOperator((n3, n3), matvec=lambda v: M @ v + k_vol * float(g @ v) * g, dtype=cp.float64)
    else:
        op = M
    dx, _ = cg(op, F, rtol=cg_tol, maxiter=cg_maxiter)
    return x + dx


__all__ = ["assemble_elastic_stiffness", "FFImplicitStepper", "assemble_K_current", "implicit_step_current",
           "assemble_K_current_cupy", "ff_implicit_step_gpu"]
