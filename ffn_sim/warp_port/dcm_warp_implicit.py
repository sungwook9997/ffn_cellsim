"""Matrix-free linearly-implicit overdamped integrator for the Warp DCM engine (Phase C accel I1).

The explicit BAOAB step is CFL-capped by the stiffest force (`dt < ~γ/k_stiffest`), so it cannot
reach biological timescales (see PHASE_C_IMPLICIT_IMEX_ACCEL_PLAN_2026-06-22.md). This module is
the implicit path: it integrates the SAME overdamped law `γ ẋ = F(x)` with a linearly-implicit
(semi-implicit) Euler step whose stability is unconditional — dt is then bounded by ACCURACY, not
stability, unlocking 100–1000× larger steps.

Per-node drag γ is diagonal (a scalar here). The implicit-Euler step for `γ ẋ = F`:

    γ (xₙ₊₁ − xₙ)/dt = F(xₙ₊₁)
    ⇒ linearise F(xₙ₊₁) ≈ F(xₙ) + J·Δx ,  J = ∂F/∂x ,  K = −J (stiffness)
    ⇒ (γ/dt · I + K) Δx = F(xₙ)            (the SPD-dominant linear solve)
    ⇒ xₙ₊₁ = xₙ + Δx

We never form K. **Matrix-free**: `K·v = −J·v ≈ −[F(x+εv) − F(x)]/ε` — one extra force evaluation
per CG iteration (the existing Warp force kernels ARE the operator). Solved with conjugate
gradient. This file is **integrator-agnostic about the force**: pass any `force_fn(pos)→force`
(numpy or device-backed); I1 validates the solver itself on a controlled spring system, I2 wires
it to the DCM stiff force.

kT=0 (deterministic overdamped) — no thermostat, so no implicit-Langevin subtlety (that is I5).
"""

from __future__ import annotations

import numpy as np


def implicit_overdamped_step(pos: np.ndarray, force_fn, gamma: float, dt: float,
                             *, eps: float = 1e-9, cg_tol: float = 1e-8, cg_maxiter: int = 200):
    """One linearly-implicit overdamped Euler step. ``force_fn(pos)->(N,3)`` is the (stiff) force.

    Solves (γ/dt I + K)Δx = F(xₙ) by matrix-free CG (K·v via a finite-difference JVP of force_fn),
    returns (pos+Δx, info). γ scalar (diagonal drag). Pure-numpy CG here (the validation/reference
    prototype); the production all-device Warp CG is the I-later optimisation."""
    from scipy.sparse.linalg import cg, LinearOperator

    x0 = np.ascontiguousarray(pos, dtype=np.float64)
    n3 = x0.size
    F0 = np.asarray(force_fn(x0), dtype=np.float64).reshape(-1)
    a = gamma / dt

    def matvec(v):
        V = v.reshape(x0.shape)
        # K·v = −J·v ≈ −[F(x+εv) − F(x)]/ε   (scale ε to v for conditioning)
        s = eps / (np.linalg.norm(v) / max(np.sqrt(n3), 1.0) + 1e-30)
        Fp = np.asarray(force_fn(x0 + s * V), dtype=np.float64).reshape(-1)
        Kv = -(Fp - F0) / s
        return a * v + Kv

    A = LinearOperator((n3, n3), matvec=matvec, dtype=np.float64)
    dx, status = cg(A, F0, rtol=cg_tol, maxiter=cg_maxiter)
    x1 = x0 + dx.reshape(x0.shape)
    return x1, {"cg_status": int(status), "dx_norm": float(np.linalg.norm(dx))}


def explicit_overdamped_step(pos: np.ndarray, force_fn, gamma: float, dt: float):
    """One explicit overdamped Euler step xₙ₊₁ = xₙ + (dt/γ) F(xₙ) — the CFL-capped reference."""
    x0 = np.ascontiguousarray(pos, dtype=np.float64)
    return x0 + (dt / gamma) * np.asarray(force_fn(x0), dtype=np.float64), {}


# ---------------------------------------------------------------------------
# I1 validation — coupled stiff spring chain (K tridiagonal ⇒ CG genuinely iterates)
# ---------------------------------------------------------------------------
def _spring_demo():
    """Relax a pinned stiff bond-chain. Show: explicit needs dt < CFL; implicit is STABLE at
    100× that dt and reaches the SAME equilibrium (the acceleration claim, concretely)."""
    rng = np.random.default_rng(0)
    N, k, L, gamma = 24, 1.0e3, 1.0, 1.0       # stiff springs (k=1e3), unit drag
    eq = (np.arange(N) * L).astype(float)       # even-spaced equilibrium (pinned ends)
    eq3 = np.zeros((N, 3)); eq3[:, 0] = eq

    def force_fn(p):
        p = np.asarray(p, dtype=np.float64).reshape(N, 3)
        x = p[:, 0].copy()
        F = np.zeros((N, 3))
        f = np.zeros(N)
        d = (x[1:] - x[:-1]) - L                # bond stretch
        f[:-1] += k * d
        f[1:] -= k * d
        f[0] = 0.0; f[-1] = 0.0                  # pinned ends
        F[:, 0] = f
        return F

    x0 = eq3.copy(); x0[1:-1, 0] += rng.normal(scale=0.3, size=N - 2)   # perturb interior
    dt_cfl = gamma / (4.0 * k)                   # ~explicit stability limit (λ_max≈4k)
    err0 = np.linalg.norm(x0 - eq3)

    def run(stepper, dt, nsteps):
        x = x0.copy()
        for s in range(nsteps):
            x, _ = stepper(x, force_fn, gamma, dt)
            if not np.isfinite(x).all() or np.linalg.norm(x - eq3) > 1e6:
                return x, s + 1, False           # diverged (non-finite OR runaway)
        return x, nsteps, True

    T = 0.5                                      # > slow-mode relax time (γN²/π²k ≈ 0.06)
    def err(x): return float(np.linalg.norm(x - eq3))
    dt_e = 0.5 * dt_cfl
    xe, ne, oke = run(explicit_overdamped_step, dt_e, int(T / dt_e))
    dt_big = 100 * dt_cfl
    xeb, neb, okeb = run(explicit_overdamped_step, dt_big, int(T / dt_big))
    xi, ni, oki = run(implicit_overdamped_step, dt_big, int(T / dt_big))

    print(f"N={N} k={k:.0e} γ={gamma}  dt_CFL≈{dt_cfl:.2e}  initial |x−eq|={err0:.3f}")
    print(f"  explicit @0.5·CFL  (dt={dt_e:.1e}, {ne} steps):  stable={oke}  |x−eq|={err(xe):.2e}  → relaxed")
    print(f"  explicit @100·CFL  (dt={dt_big:.1e}, {neb} steps): stable={okeb}  |x−eq|={err(xeb):.2e}  ← BLOW-UP")
    print(f"  IMPLICIT @100·CFL  (dt={dt_big:.1e}, {ni} steps):  stable={oki}  |x−eq|={err(xi):.2e}  ← STABLE+relaxed")
    print(f"  implicit≈explicit equilibrium: |xi−xe|={np.linalg.norm(xi-xe):.2e}  (×{int(dt_big/dt_e)} fewer steps)")
    return dict(dt_cfl=dt_cfl, explicit_ok=oke, explicit_big_ok=okeb, implicit_ok=oki,
                implicit_eq_err=err(xi), implicit_vs_explicit=float(np.linalg.norm(xi - xe)),
                speedup_steps=int(dt_big / dt_e))


if __name__ == "__main__":
    _spring_demo()
