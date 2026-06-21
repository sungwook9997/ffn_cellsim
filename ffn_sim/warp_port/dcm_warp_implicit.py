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
    iters = [0]
    def _cb(xk): iters[0] += 1
    dx, status = cg(A, F0, rtol=cg_tol, maxiter=cg_maxiter, callback=_cb)
    x1 = x0 + dx.reshape(x0.shape)
    return x1, {"cg_status": int(status), "cg_iters": iters[0], "dx_norm": float(np.linalg.norm(dx))}


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


def make_dcm_stiff_force(device="cpu", subdiv=2):
    """Build a single-cell DCM **stiff** force callable `stiff_force(pos_np)->force_np` from the
    real Warp kernels (turgor + cortex edges — the stiff terms that set the CFL). Used by I2 to
    validate the implicit step against the actual DCM physics (not just a toy spring)."""
    import warp as wp
    from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM
    from ffn_sim.warp_port.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
    from ffn_sim.warp_port.dcm_warp_hybrid import _bond_accumulate
    from ffn_sim.warp_port.dcm_warp_hybrid_multicell import _dp_from_vol, _zero_vec

    p = ResolvedDCM(subdivisions=subdiv)
    verts, edges, faces = icosphere_mesh(p.R_cell, subdiv)
    N, nf, ne = verts.shape[0], faces.shape[0], edges.shape[0]
    R0 = float(np.linalg.norm(verts, axis=1).mean()); V0 = (4/3)*np.pi*R0**3
    fcell = np.zeros(nf, np.int32)
    pos_d = wp.array(verts, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    faces_d = wp.array(faces.astype(np.int32), dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell, dtype=wp.int32, device=device)
    edges_d = wp.array(edges.astype(np.int32), dtype=wp.int32, device=device)
    r0_d = wp.array(np.linalg.norm(verts[edges[:,0]]-verts[edges[:,1]], axis=1), dtype=wp.float64, device=device)
    Vc_d = wp.zeros(1, dtype=wp.float64, device=device); dP_d = wp.zeros(1, dtype=wp.float64, device=device)
    gamma = 6.0*np.pi*65.9*p.R_cell/N

    def stiff_force(pos_np):
        pos_d.assign(np.ascontiguousarray(pos_np.reshape(N,3)))
        wp.launch(_zero_vec, dim=N, inputs=[force_d], device=device)
        Vc_d.zero_()
        wp.launch(dcm_volume_kernel, dim=nf, inputs=[pos_d, faces_d, fcell_d, Vc_d], device=device)
        wp.launch(_dp_from_vol, dim=1, inputs=[Vc_d, wp.float64(V0), wp.float64(p.turgor_dP0), wp.float64(7.73e5), dP_d], device=device)
        wp.launch(dcm_turgor_force_kernel, dim=nf, inputs=[pos_d, faces_d, fcell_d, dP_d, force_d], device=device)
        wp.launch(_bond_accumulate, dim=ne, inputs=[pos_d, edges_d, wp.float64(p.k_edge), r0_d, force_d], device=device)
        wp.synchronize_device(device)
        return force_d.numpy().astype(np.float64)

    return stiff_force, verts.copy(), gamma, V0, faces, fcell


def _dcm_stiff_demo(device="cpu"):
    """I2: relax a SQUASHED single DCM cell under the real turgor+edge stiff force. Explicit @8e-6
    (production dt) vs explicit @100× (blow-up) vs IMPLICIT @100× (stable, same rest shape)."""
    stiff_force, verts, gamma, V0, faces, fcell = make_dcm_stiff_force(device)
    N = verts.shape[0]
    def vol(x):
        v0,v1,v2 = x[faces[:,0]],x[faces[:,1]],x[faces[:,2]]
        return abs(float(np.einsum('ij,ij->i', v0, np.cross(v1-v0,v2-v0)).sum()/6.0))
    x0 = verts.copy(); x0[:,2] *= 0.7            # squash to 70% in z (off-equilibrium)
    dt_e = 8.0e-6                                 # production explicit dt (known stable)
    dt_big = 100*dt_e

    def run(stepper, dt, nsteps):
        x = x0.copy()
        for _ in range(nsteps):
            x,_ = stepper(x, stiff_force, gamma, dt)
            if not np.isfinite(x).all() or vol(x) > 50*V0: return x, False
        return x, True
    T = 0.02
    xe, oke = run(explicit_overdamped_step, dt_e, int(T/dt_e))
    xeb, okeb = run(explicit_overdamped_step, dt_big, int(T/dt_big))
    xi, oki = run(implicit_overdamped_step, dt_big, int(T/dt_big))
    print(f"single DCM cell (N={N}), squashed z×0.7  V0={V0:.2e}  γ={gamma:.2e}  dt_exp={dt_e:.0e}")
    print(f"  explicit @dt   ({int(T/dt_e)} steps): stable={oke}  V/V0={vol(xe)/V0:.3f}")
    print(f"  explicit @100× ({int(T/dt_big)} steps): stable={okeb}  ← expect BLOW-UP")
    print(f"  IMPLICIT @100× ({int(T/dt_big)} steps): stable={oki}  V/V0={vol(xi)/V0:.3f}  |xi−xe|={np.linalg.norm(xi-xe):.2e}")
    return dict(explicit_ok=oke, explicit_big_ok=okeb, implicit_ok=oki,
                vv0_implicit=vol(xi)/V0, implicit_vs_explicit=float(np.linalg.norm(xi-xe)))


def _dt_ramp(device="cpu"):
    """I3: ramp dt (1/10/100/1000/10000× the explicit dt) on the real DCM stiff force. Measure
    stability, equilibrium V/V0, error vs the converged-explicit reference, CG iters/step, and the
    step-count speedup. Finds the accuracy ceiling (where the equilibrium stops matching), NOT a
    stability limit (implicit is unconditionally stable)."""
    stiff_force, verts, gamma, V0, faces, fcell = make_dcm_stiff_force(device)
    N = verts.shape[0]
    def vol(x):
        v0,v1,v2 = x[faces[:,0]],x[faces[:,1]],x[faces[:,2]]
        return abs(float(np.einsum('ij,ij->i', v0, np.cross(v1-v0,v2-v0)).sum()/6.0))
    x0 = verts.copy(); x0[:,2] *= 0.7
    dt_e = 8.0e-6; T = 0.02

    # converged-explicit reference equilibrium
    xref = x0.copy()
    for _ in range(int(T/dt_e)):
        xref,_ = explicit_overdamped_step(xref, stiff_force, gamma, dt_e)
    print(f"reference (explicit dt={dt_e:.0e}, {int(T/dt_e)} steps): V/V0={vol(xref)/V0:.4f}")
    print(f"{'dt/dt_e':>8} {'steps':>7} {'stable':>7} {'V/V0':>7} {'|x-ref|':>9} {'cgit/step':>9} {'speedup':>8}")
    for mult in (1, 10, 100, 1000, 10000):
        dt = dt_e*mult; ns = max(1, int(T/dt)); x = x0.copy(); tot_it=0; ok=True
        for _ in range(ns):
            x, info = implicit_overdamped_step(x, stiff_force, gamma, dt)
            tot_it += info["cg_iters"]
            if not np.isfinite(x).all() or vol(x) > 50*V0: ok=False; break
        e = float(np.linalg.norm(x-xref)) if ok else float('inf')
        print(f"{mult:>8} {ns:>7} {str(ok):>7} {vol(x)/V0 if ok else 0:>7.4f} {e:>9.2e} "
              f"{tot_it/max(ns,1):>9.1f} {int(T/dt_e)/ns:>7.0f}x")


if __name__ == "__main__":
    print("== I1 spring solver =="); _spring_demo()
    print("== I2 DCM stiff force =="); _dcm_stiff_demo()
    print("== I3 dt ramp =="); _dt_ramp()
