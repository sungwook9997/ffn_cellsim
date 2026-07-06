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


__all__ = ["assemble_elastic_stiffness", "FFImplicitStepper"]
