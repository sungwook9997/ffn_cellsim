"""SC2 deflated-PCG gate (CPU, numpy) — parity + iteration reduction + failure guards.

Mimics the FF implicit operator's two SC0-measured conditioning targets on a small SPD system:
  A = a·I + K (bending/crosslink SPD) with (1) an a_com rank-3 rigid-translation softening (near-singular mode)
  and (2) a stiff heterogeneous "basal" diagonal. Verifies the two-level deflated PCG (rigid-mode deflation +
  stiff-diag Jacobi smoother) converges to the SAME solution as a direct solve, in FEWER iterations than plain CG.
"""
from __future__ import annotations

import numpy as np
import pytest

from aleph.laws.solver.deflated_pcg import deflated_pcg, build_rigid_translation_basis


def _chain_laplacian(n):
    """SPD-ish graph Laplacian on a chain + small pin (so it is SPD, mimics bending/crosslink stiffness)."""
    L = np.zeros((n, n))
    for i in range(n - 1):
        L[i, i] += 1.0; L[i + 1, i + 1] += 1.0; L[i, i + 1] -= 1.0; L[i + 1, i] -= 1.0
    return L


def _build_operator(N=30, a=100.0, a_com=0.5, k_stiff=5.0e3, n_basal=4, seed=0):
    """Return (matvec, A_dense, a, diag_extra) for a 3N SPD system with the two conditioning targets."""
    rng = np.random.default_rng(seed)
    n3 = 3 * N
    Lx = _chain_laplacian(N)
    K = np.zeros((n3, n3))
    for d in range(3):                                   # bending/crosslink stiffness per coordinate
        K[d::3][:, d::3] += 20.0 * Lx
    K += 0.5 * (lambda M: M + M.T)(rng.standard_normal((n3, n3))) * 0.0   # (kept SPD; no random asymmetry)
    diag_extra = np.zeros(n3)
    basal = np.arange(n_basal)                           # stiff "basal" nodes → large diagonal (k_int/k_plane)
    for i in basal:
        diag_extra[3 * i:3 * i + 3] += k_stiff
    A = a * np.eye(n3) + K + np.diag(diag_extra)
    # a_com rank-3 softening: lower the rigid-translation eigenvalue a→a_com (symmetric, keeps SPD since a_com>0)
    for d in range(3):
        e = np.zeros(n3); e[d::3] = 1.0 / np.sqrt(N)     # unit rigid-translation mode in axis d
        A = A - (a - a_com) * np.outer(e, e)

    def matvec(v):
        return A @ v

    return matvec, A, a, diag_extra


def _plain_cg_iters(matvec, b, tol=1e-8, maxiter=2000):
    """Reference unpreconditioned CG iteration count to the same tol."""
    x = np.zeros_like(b); r = b - matvec(x); p = r.copy(); rr = r @ r; bn = np.linalg.norm(b); it = 0
    while np.sqrt(rr) > tol * bn and it < maxiter:
        Ap = matvec(p); alpha = rr / (p @ Ap); x += alpha * p; r -= alpha * Ap
        rr_new = r @ r; p = r + (rr_new / rr) * p; rr = rr_new; it += 1
    return x, it


def test_deflated_pcg_parity_and_speedup():
    N, a, a_com, k_stiff = 30, 100.0, 0.5, 5.0e3
    matvec, A, aval, diag_extra = _build_operator(N=N, a=a, a_com=a_com, k_stiff=k_stiff)
    n3 = 3 * N
    rng = np.random.default_rng(1)
    b = rng.standard_normal(n3)
    x_direct = np.linalg.solve(A, b)                     # ground truth

    # deflation basis = rigid translation (the a_com target); smoother = 1/(a + diag_extra) (the stiff-diag target)
    W = build_rigid_translation_basis(N, N, np, extra_cols=None)
    D_inv = 1.0 / (aval + diag_extra)

    x_def, info = deflated_pcg(matvec, b, np, D_inv=D_inv, W=W, rtol=1e-10, maxiter=2000, check_every=1)
    _, it_plain = _plain_cg_iters(matvec, b, tol=1e-10)

    # (1) PARITY: deflated PCG solves the SAME system as the direct solve
    assert np.linalg.norm(x_def - x_direct) / np.linalg.norm(x_direct) < 1e-6, info
    assert info["converged"] and info["breakdown"] is None
    # (2) SPEEDUP: fewer iterations than unpreconditioned CG (the whole point)
    assert info["iters"] < it_plain, f"deflated {info['iters']} not < plain {it_plain}"
    # the near-singular a_com mode makes plain CG slow; deflation should cut it hard
    assert info["iters"] <= 0.7 * it_plain


def test_deflation_removes_the_rigid_mode_penalty():
    """a_com softening inflates plain CG; deflating the rigid mode should restore ~the no-softening iters."""
    N, a, k_stiff = 24, 100.0, 0.0                       # no stiff diag → isolate the a_com effect
    b = np.random.default_rng(2).standard_normal(3 * N)

    mv_soft, A_soft, aval, de = _build_operator(N=N, a=a, a_com=0.2, k_stiff=k_stiff)
    mv_stiff, A_stiff, _, _ = _build_operator(N=N, a=a, a_com=a, k_stiff=k_stiff)   # a_com=a ⇒ no softening
    _, it_soft = _plain_cg_iters(mv_soft, b, tol=1e-10)
    _, it_stiff = _plain_cg_iters(mv_stiff, b, tol=1e-10)
    assert it_soft > it_stiff                            # softening genuinely inflates plain CG

    W = build_rigid_translation_basis(N, N, np)
    _, info = deflated_pcg(mv_soft, b, np, D_inv=1.0 / (aval + de), W=W, rtol=1e-10, maxiter=2000, check_every=1)
    # deflating the rigid mode brings the SOFT operator's iters down to ~the STIFF (unsoftened) baseline
    assert info["iters"] <= it_stiff + 3, (info["iters"], it_stiff)


def test_guards_fire_on_indefinite_operator():
    """A non-SPD (indefinite) operator must trip a guard, not silently return garbage."""
    n = 30
    A = np.diag(np.concatenate([np.ones(n - 1), [-5.0]]))   # one negative eigenvalue
    b = np.ones(n)
    _, info = deflated_pcg(lambda v: A @ v, b, np, D_inv=None, W=None, rtol=1e-10, maxiter=200)
    assert info["breakdown"] in ("pAp_floor", "rz_nonpos") or not info["converged"]


def test_structured_rigid_equals_dense_W():
    """The structured rigid-mode coarse apply (rigid_ncortex[, g_col]) must equal the dense-W apply — same
    solution and (nearly) same iteration count — but without materialising the dense 3N×k basis per iteration."""
    N, a, a_com, k_stiff = 30, 100.0, 0.5, 5.0e3
    matvec, A, aval, de = _build_operator(N=N, a=a, a_com=a_com, k_stiff=k_stiff)
    b = np.random.default_rng(7).standard_normal(3 * N)
    D_inv = 1.0 / (aval + de)
    g = np.random.default_rng(8).standard_normal(3 * N)          # a dense extra column (stand-in for osmotic g)

    W = build_rigid_translation_basis(N, N, np, extra_cols=[g])
    x_dense, i_dense = deflated_pcg(matvec, b, np, D_inv=D_inv, W=W, rtol=1e-10, maxiter=2000, check_every=1)
    x_str, i_str = deflated_pcg(matvec, b, np, D_inv=D_inv, rigid_ncortex=N, g_col=g,
                                rtol=1e-10, maxiter=2000, check_every=1)
    assert np.linalg.norm(x_str - x_dense) / np.linalg.norm(x_dense) < 1e-8, (i_str, i_dense)
    assert abs(i_str["iters"] - i_dense["iters"]) <= 1
    assert np.linalg.norm(x_str - np.linalg.solve(A, b)) / np.linalg.norm(b) < 1e-6


def test_no_deflation_equals_plain_pcg_parity():
    """W=None, D_inv=None ⇒ deflated_pcg reduces to plain CG and still solves an SPD system."""
    N = 20
    matvec, A, aval, de = _build_operator(N=N, a=50.0, a_com=50.0, k_stiff=0.0)  # no softening/stiff
    b = np.random.default_rng(3).standard_normal(3 * N)
    x_def, info = deflated_pcg(matvec, b, np, D_inv=None, W=None, rtol=1e-10, maxiter=2000, check_every=1)
    assert np.linalg.norm(x_def - np.linalg.solve(A, b)) / np.linalg.norm(b) < 1e-6
    assert info["converged"]
