"""Two-level deflated PCG for the FF implicit overdamped operator (SC2, FF_SINGLE_CELL_OPTIMIZATION_PLAN §4.4).

SC0 (native profile) found the implicit step is CG-DOMINATED (86-89% of wall) and identified TWO conditioning
targets that inflate the unpreconditioned CG iteration count:

  1. the ``a_com`` modal rigid-COM drag pulls the 3 rigid-translation eigenvalues from ``a=γ/dt`` down to
     ``a_com=6πηR/Nc/dt ≪ a`` (near-singular) — CG stalls on it (measured 71→135 iters);
  2. the stiff, heterogeneous basal ``k_int``/``k_plane`` diagonal widens the spectrum (measured +47-52% iters);
     it is diagonally dominant on the basal dofs, so a *selective* Jacobi smoother helps there (unlike the
     crosslink rank-1 blocks, where full Jacobi INFLATES).

This module attacks both with an SPD two-level preconditioner used inside standard PCG:

    M⁻¹ r = D⁻¹ ⊙ r  +  W (Wᵀ A W)⁻¹ Wᵀ r
            └ smoother ┘   └── coarse deflation correction ──┘

- **D⁻¹ = 1/(a + diag_extra)** — Jacobi on the *diagonally-dominant stiff part only* (interior dofs get the
  harmless uniform ``1/a``; basal dofs get ``1/(a+k_int/k_plane)``). It deliberately does NOT try to
  precondition the crosslink k·ûûᵀ blocks (those are not diagonally dominant → Jacobi inflates).
- **W** = the closed-form troublesome modes (rigid translation over the cortex nodes, + the osmotic volume
  gradient ``g``). Because W is known analytically, the coarse solve ``E⁻¹`` (E = Wᵀ A W, tiny k×k) is exact,
  cheap, and robust — no heuristic eigenvector estimation. This is the principled replacement for the
  unstable ``--bulk-drag`` uniform-γ hack (which created a near-null-space and diverged).

Both terms are SPD (``D⁻¹>0``; ``W E⁻¹ Wᵀ`` SPD since A SPD, W full-rank) ⇒ the sum is an SPD preconditioner ⇒
standard PCG converges. Guards (pAp floor, non-finite, residual growth, maxiter) match the production solver.

The operator ``A`` is passed as a matvec so this works identically for the assembled-CSR path and any future
matrix-free operator, and for the crawl / AFM / ECM variants (they differ only in which forces enter A).
"""
from __future__ import annotations


def build_rigid_translation_basis(N: int, n_cortex: int, xp, *, extra_cols=None):
    """Closed-form deflation basis W (3N × k): unit rigid translation of the cortex nodes in x, y, z
    (+ optional extra columns, e.g. the osmotic volume gradient). Columns are individually normalised.

    Args:
        N: total node count; ``3N`` dofs.
        n_cortex: cortex nodes ``[0, n_cortex)`` participate in the rigid-translation modes (the a_com target).
        xp: the array module (``cupy`` on device, ``numpy`` for CPU tests).
        extra_cols: optional list of (3N,) arrays to append (e.g. ``vol_g``); each is normalised, zero-cols dropped.

    Returns:
        W (3N × k) in ``xp``, columns unit-norm and (numerically) A-independent analytic modes.
    """
    n3 = 3 * N
    cols = []
    for d in range(3):
        w = xp.zeros(n3, dtype=xp.float64)
        idx = xp.arange(n_cortex) * 3 + d
        w[idx] = 1.0
        nrm = float(xp.linalg.norm(w))
        if nrm > 0.0:
            cols.append(w / nrm)
    for c in (extra_cols or []):
        c = xp.asarray(c, dtype=xp.float64).reshape(-1)
        nrm = float(xp.linalg.norm(c))
        if nrm > 1e-300:
            cols.append(c / nrm)
    return xp.stack(cols, axis=1) if cols else None


def _structured_coarse(matvec, xp, N, n_cortex, g_col):
    """Structured coarse operator for the rigid-translation (+ optional osmotic ``g``) deflation space, applied
    WITHOUT a dense 3N×k matrix. The rigid modes ``u_d`` are uniform ``1/√Nc`` on the cortex axis-``d`` dofs, so
    ``Wᵀr`` = 3 strided cortex axis-sums (+ one ``g·r`` dot) and ``Wc`` = 3 strided broadcasts (+ one ``g`` axpy) —
    O(Nc) per iteration, not a 47 MB GEMV (the naive dense apply measured ~2× wall-SLOWER despite 1.7-1.9× fewer
    iters). Returns ``coarse(r) -> W (WᵀAW)⁻¹ Wᵀ r``.

    Setup materialises the k basis columns ONCE (k matvecs for AW + the k×k coarse matrix E); only the small
    ``g`` column is kept dense thereafter.
    """
    n3 = 3 * N
    inv_sq = 1.0 / (n_cortex ** 0.5)
    g_hat = None
    if g_col is not None:
        g_col = xp.asarray(g_col, dtype=xp.float64).reshape(-1)
        gn = float(xp.linalg.norm(g_col))
        g_hat = (g_col / gn) if gn > 1e-300 else None
    k = 3 + (1 if g_hat is not None else 0)

    def wt_apply(r):                                             # Wᵀ r  → (k,) on device
        s = [r[d::3][:n_cortex].sum() * inv_sq for d in range(3)]
        if g_hat is not None:
            s.append(g_hat @ r)
        return xp.stack(s)

    def wc_apply(c):                                            # W c  → (3N,) on device
        out = xp.zeros(n3, dtype=xp.float64)
        for d in range(3):
            out[d::3][:n_cortex] = c[d] * inv_sq
        if g_hat is not None:
            out = out + c[3] * g_hat
        return out

    # E = WᵀAW: materialise the k columns once for the k setup matvecs (one-time transient), then E via wt_apply.
    U_cols = []
    for d in range(3):
        u = xp.zeros(n3, dtype=xp.float64); u[d::3][:n_cortex] = inv_sq; U_cols.append(u)
    if g_hat is not None:
        U_cols.append(g_hat)
    AW = xp.stack([matvec(u) for u in U_cols], axis=1)          # (3N × k)  — k setup matvecs
    E = xp.stack([wt_apply(AW[:, j]) for j in range(k)], axis=1)   # (k × k) = WᵀAW (SPD)
    E_inv = xp.linalg.inv(E)

    def coarse(r):
        return wc_apply(E_inv @ wt_apply(r))
    return coarse


def deflated_pcg(matvec, b, xp, *, D_inv=None, W=None, rigid_ncortex=None, g_col=None,
                 rtol=1e-6, maxiter=400, callback=None, atol=0.0, check_every=10):
    """Preconditioned CG for SPD ``A`` (given as ``matvec``) with the two-level preconditioner above.

    DEVICE-SCALAR: the CG scalars (rz, pAp, α, β) stay on-device (0-d arrays) — no per-iteration host sync — and
    the residual norm / guards are read only every ``check_every`` iterations (the ENGINE_ACCELERATION_PLAN
    "read residual on a cadence, never drop the divergence exit" rule). This is what makes the iteration
    reduction from deflation actually turn into a WALL-TIME reduction (a per-iter ``float()`` sync erases it).

    Args:
        matvec: callable ``v -> A v`` (both (3N,) in ``xp``). Must be the SAME operator being solved.
        b: right-hand side (3N,) in ``xp``.
        xp: array module (``cupy`` on device / ``numpy`` for CPU tests).
        D_inv: (3N,) Jacobi smoother = ``1/(a + diag_extra)`` (diagonally-dominant part). ``None`` → identity.
        W: (3N × k) deflation basis (rigid translation [+ g]). ``None`` → no coarse correction (plain PCG).
        rtol, atol: stop when ``||r|| <= max(rtol*||b||, atol)``.
        maxiter: cap; check_every: residual/guard readback cadence; callback(it) fired once per iteration.

    Returns:
        (x, info) — x (3N,), info dict {iters, rel_resid, converged, breakdown}.
    """
    # ---- coarse operator: structured rigid-mode apply (preferred) or dense W (tests/generality) ----
    coarse = None
    if rigid_ncortex is not None:                                    # STRUCTURED: O(Nc) reductions/broadcasts (no 47MB GEMV)
        N = b.shape[0] // 3
        coarse = _structured_coarse(matvec, xp, N, int(rigid_ncortex), g_col)
    elif W is not None and W.shape[1] > 0:                           # DENSE fallback (small systems / arbitrary basis)
        k = W.shape[1]
        AW = xp.stack([matvec(W[:, j]) for j in range(k)], axis=1)   # (3N × k)
        E = W.T @ AW                                                  # (k × k) SPD
        E_inv = xp.linalg.inv(E)
        Wt = xp.ascontiguousarray(W.T)                               # (k × 3N) contiguous → coalesced Wᵀr GEMV

        def coarse(r):
            return W @ (E_inv @ (Wt @ r))                            # W (WᵀAW)⁻¹ Wᵀ r  (all on-device)

    def apply_Minv(r):
        z = (D_inv * r) if D_inv is not None else r
        if coarse is not None:
            z = z + coarse(r)
        return z

    x = xp.zeros_like(b)
    r = b - matvec(x)                                                # x0 = 0 ⇒ r = b
    bnorm = float(xp.linalg.norm(b))                                 # one host read at setup
    stop = max(rtol * bnorm, atol)
    stop2 = stop * stop
    z = apply_Minv(r)
    p = z.copy()
    rz = r @ z                                                       # 0-d device scalar (NO float())
    it = 0
    breakdown = None
    converged = bnorm <= stop
    while not converged and it < maxiter:
        Ap = matvec(p)
        pAp = p @ Ap                                                 # 0-d device scalar
        alpha = rz / pAp
        x = x + alpha * p
        r = r - alpha * Ap
        it += 1
        if callback is not None:
            callback(it)
        if it % check_every == 0 or it >= maxiter:                  # guards + convergence on a cadence
            rr = float(r @ r); pApf = float(pAp)
            if pApf != pApf or rr != rr:                            # non-finite guard
                breakdown = "nonfinite"; break
            if not (pApf > 1e-300):                                 # pAp floor
                breakdown = "pAp_floor"; break
            if rr <= stop2:
                converged = True; break
        z = apply_Minv(r)
        rz_new = r @ z
        beta = rz_new / rz
        p = z + beta * p
        rz = rz_new
    rel = float(xp.linalg.norm(b - matvec(x))) / (bnorm if bnorm > 0 else 1.0)
    if not converged and breakdown is None and rel <= rtol:
        converged = True
    return x, {"iters": it, "rel_resid": rel, "converged": bool(converged), "breakdown": breakdown}
