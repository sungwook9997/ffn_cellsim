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

import math
import os

import numpy as np


def implicit_overdamped_step(pos: np.ndarray, force_fn, gamma: float, dt: float,
                             *, eps: float = 1e-9, cg_tol: float = 1e-8, cg_maxiter: int = 200,
                             n_newton: int = 1, newton_tol: float = 1e-10, precond: bool = False):
    """One implicit overdamped Euler step via Newton (I4). The implicit-Euler nonlinear equation is
    G(x) = γ(x−xₙ)/dt − F(x) = 0; Newton solves J_G Δx = −G with J_G = γ/dt·I − ∂F/∂x = γ/dt·I + K,
    matrix-free (K·v = −[F(x+εv)−F(x)]/ε), CG per Newton iter. ``n_newton=1`` = the linearly-implicit
    step (I1/I2); ``n_newton>1`` corrects the large-Δx nonlinearity that drifts V/V0 at huge dt (I3
    finding). Returns (xₙ₊₁, info). γ scalar diagonal drag. scipy-CG prototype (all-device = I-opt)."""
    from scipy.sparse.linalg import cg, LinearOperator

    x0 = np.ascontiguousarray(pos, dtype=np.float64)
    n3 = x0.size
    a = gamma / dt
    x = x0.copy()
    rng = np.random.default_rng(0)
    tot_cg = 0; n_used = 0
    for _ in range(n_newton):
        Fx = np.asarray(force_fn(x), dtype=np.float64).reshape(-1)
        G = a * (x.reshape(-1) - x0.reshape(-1)) - Fx        # residual (=0 at the implicit solution)
        gnorm = np.linalg.norm(G)
        n_used += 1
        if gnorm < newton_tol * max(a, 1.0):
            break

        def matvec(v):
            V = v.reshape(x.shape)
            s = eps / (np.linalg.norm(v) / max(np.sqrt(n3), 1.0) + 1e-30)
            Fp = np.asarray(force_fn(x + s * V), dtype=np.float64).reshape(-1)
            return a * v - (Fp - Fx) / s                     # (γ/dt I + K) v
        A = LinearOperator((n3, n3), matvec=matvec, dtype=np.float64)
        # Jacobi preconditioner (OPT-IN, default OFF). MEASURED to HURT here: a Hutchinson diagonal
        # from the noisy forward-diff JVP is unreliable → worse conditioning (2-cell contact: plain
        # CG 3 iters vs this 157). Near equilibrium the γ/dt regulariser already conditions the
        # operator well (plain CG ~3-5 iters), so no preconditioner is needed. A correct one needs
        # the ANALYTIC diagonal stiffness (per-term contact/turgor/edge) — future work if far-from-
        # equilibrium steps ever need it. Contact does NOT dominate the conditioning here.
        M = None
        if precond:
            diagA = np.zeros(n3)
            for _pi in range(4):
                v = rng.choice(np.array([-1.0, 1.0]), size=n3)
                diagA += v * matvec(v)
            diagA = np.maximum(diagA / 4.0, a)
            M = LinearOperator((n3, n3), matvec=lambda r: r / diagA, dtype=np.float64)
        it = [0]
        dx, _ = cg(A, -G, rtol=cg_tol, maxiter=cg_maxiter, M=M,
                   callback=lambda xk: it.__setitem__(0, it[0] + 1))
        tot_cg += it[0]
        x = x + dx.reshape(x.shape)
    return x, {"cg_iters": tot_cg, "newton_iters": n_used, "g_norm": float(gnorm)}


import warp as wp
wp.init()


@wp.kernel
def _vdot(a: wp.array(dtype=wp.vec3d), b: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.float64)):
    i = wp.tid(); wp.atomic_add(out, 0, wp.dot(a[i], b[i]))

@wp.kernel
def _vaxpy(y: wp.array(dtype=wp.vec3d), alpha: wp.float64, x: wp.array(dtype=wp.vec3d)):  # y += αx
    i = wp.tid(); y[i] = y[i] + alpha * x[i]

@wp.kernel
def _vxpby(out: wp.array(dtype=wp.vec3d), x: wp.array(dtype=wp.vec3d), beta: wp.float64,
           p: wp.array(dtype=wp.vec3d)):  # out = x + β·p
    i = wp.tid(); out[i] = x[i] + beta * p[i]

@wp.kernel
def _vcopy(dst: wp.array(dtype=wp.vec3d), src: wp.array(dtype=wp.vec3d)):
    i = wp.tid(); dst[i] = src[i]

@wp.kernel
def _vaxpy_active(y: wp.array(dtype=wp.vec3d), alpha: wp.float64, x: wp.array(dtype=wp.vec3d),
                  cof: wp.array(dtype=wp.int32)):  # y += αx for LIVE nodes only (cof>=0)
    i = wp.tid()
    if cof[i] >= wp.int32(0):
        y[i] = y[i] + alpha * x[i]


@wp.kernel
def _vaxpy_active_capped(y: wp.array(dtype=wp.vec3d), x: wp.array(dtype=wp.vec3d),
                         cof: wp.array(dtype=wp.int32), cap: wp.float64):
    """D8: y += Δx for LIVE nodes, but with each node's step clamped to ‖Δx‖ ≤ cap. The
    frozen-grid implicit step evaluates contact on xₙ, so a node solved to move > the thin
    contact shell can TUNNEL through a face in one big step (the pen-runaway). Clamping the
    per-node displacement to the contact-shell scale (cap = c_rep, geometric — not tuned)
    means no node can cross the shell un-checked; the residual motion is taken next step."""
    i = wp.tid()
    if cof[i] >= wp.int32(0):
        dx = x[i]
        n = wp.length(dx)
        if n > cap and n > wp.float64(0.0):
            dx = dx * (cap / n)
        y[i] = y[i] + dx

@wp.kernel
def _operator(out: wp.array(dtype=wp.vec3d), a: wp.float64, v: wp.array(dtype=wp.vec3d),
              Fp: wp.array(dtype=wp.vec3d), Fx: wp.array(dtype=wp.vec3d), inv_s: wp.float64):
    # out = a·v + K·v,  K·v = −(F(x+s v) − F(x))/s = −(Fp − Fx)·inv_s
    i = wp.tid(); out[i] = a * v[i] - (Fp[i] - Fx[i]) * inv_s


def device_cg(stiff_into, x_d, a, b_d, scratch, *, tol=1e-8, maxiter=200, eps=1e-9, device="cpu"):
    """All-device matrix-free CG for (a·I + K)Δx = b. Big vectors stay on the GPU; only the CG
    scalars (dot products) cross to host. ``stiff_into(pos_d, out_d)`` writes F(pos) on device;
    K·v via a perturbed stiff eval. ``scratch`` = dict of pre-allocated device vec3d buffers."""
    N = x_d.shape[0]
    r, p, Ap, dx, Fx, Fp, xp, sca = (scratch[k] for k in ("r", "p", "Ap", "dx", "Fx", "Fp", "xp", "sca"))

    def dot(u, v):
        sca.zero_(); wp.launch(_vdot, dim=N, inputs=[u, v, sca], device=device)
        wp.synchronize_device(device); return float(sca.numpy()[0])

    stiff_into(x_d, Fx)                                  # F(x) once
    _dbg = os.environ.get("CG_DEBUG")

    dx.zero_()
    wp.launch(_vcopy, dim=N, inputs=[r, b_d], device=device)     # r = b - A·0 = b
    wp.launch(_vcopy, dim=N, inputs=[p, r], device=device)
    rs = dot(r, r); bnorm = max(dot(b_d, b_d) ** 0.5, 1e-30); rs0 = rs; it = 0
    if _dbg:
        print(f"    [cg-dbg] N={N} a={a:.3e} bnorm={bnorm:.3e} rs0={rs:.3e}", flush=True)
    # Already at (or below) tolerance — F(x)≈0 ⇒ Δx≈0. Skip CG; a zero step is correct and
    # avoids the forward-difference JVP catastrophically cancelling on Fp−Fx≈0 (a near-equilibrium
    # configuration is exactly where that noise turns the matrix-free operator non-SPD).
    if rs ** 0.5 <= tol * bnorm or not math.isfinite(rs):
        return dx, 0
    for it in range(1, maxiter + 1):
        # JVP perturbation scale from the current search direction. Guard a degenerate/diverged
        # p (vn=0 or non-finite) BEFORE forming 1/s — an unfinite vn underflows s→0 and the old
        # code crashed on the 1.0/s division. A zero/garbage direction means we stop here.
        pp = dot(p, p)
        vn = (pp / N) ** 0.5
        if not math.isfinite(vn) or vn <= 1e-300:
            if _dbg:
                print(f"    [cg-dbg] it={it} degenerate vn={vn:.3e} — break", flush=True)
            break
        s = eps / vn
        wp.launch(_vxpby, dim=N, inputs=[xp, x_d, wp.float64(s), p], device=device)   # xp = x + s p
        stiff_into(xp, Fp)
        wp.launch(_operator, dim=N, inputs=[Ap, wp.float64(a), p, Fp, Fx, wp.float64(1.0 / s)], device=device)
        pAp = dot(p, Ap)
        # A = (γ/dt)·I + K. The diagonal a·‖p‖² is exact and strictly positive; the physical
        # stiffness K (penalty/turgor/edge/bending Hessians) is PSD at a stable equilibrium, so
        # pᵀAp ≥ a·‖p‖² ALWAYS. Near equilibrium the forward-difference Jacobian K·p collapses to
        # rounding noise (catastrophic cancellation on Fp−Fx≈0) and pᵀKp can go spuriously negative
        # → a runaway α=rs/pAp that injected the giant displacement we saw at 100× dt. Floor pAp at
        # the exact diagonal (ignore the FD noise) — physically justified, no scale-tuned constant.
        pAp_diag = a * pp
        if not math.isfinite(pAp) or pAp < pAp_diag:
            if _dbg and pAp < pAp_diag:
                print(f"    [cg-dbg] it={it} K-noise floor pAp={pAp:.3e}<diag={pAp_diag:.3e}", flush=True)
            pAp = pAp_diag
        alpha = rs / pAp
        wp.launch(_vaxpy, dim=N, inputs=[dx, wp.float64(alpha), p], device=device)    # x += α p
        wp.launch(_vaxpy, dim=N, inputs=[r, wp.float64(-alpha), Ap], device=device)   # r -= α Ap
        rs_new = dot(r, r)
        if _dbg:
            print(f"    [cg-dbg] it={it} rs={rs:.3e} pAp={pAp:.3e} -> rs_new={rs_new:.3e}", flush=True)
        # Divergence guard. For this non-symmetric, K-dominated operator (large dt ⇒ a=γ/dt ≪‖K‖),
        # diagonal-floored CG can diverge — the residual grows monotonically and ran away to ~1e128
        # within a single solve. Stop at the first runaway (4× the initial residual) and return the
        # current iterate. A diverging solve here is the SIGNAL that this dt exceeds the stable
        # ceiling for this stack — the caller's V/V0 will then flag it, instead of an opaque NaN.
        if not math.isfinite(rs_new) or rs_new > 4.0 * rs0:
            if _dbg:
                print(f"    [cg-dbg] it={it} DIVERGING rs_new={rs_new:.3e} > 4·rs0={4 * rs0:.3e} — break", flush=True)
            break
        if (rs_new ** 0.5) < tol * bnorm:
            rs = rs_new
            break
        beta = rs_new / rs
        wp.launch(_vxpby, dim=N, inputs=[p, r, wp.float64(beta), p], device=device)   # p = r + β p
        rs = rs_new
    return dx, it


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


def make_dcm_stiff_into(device="cpu", subdiv=2):
    """Device-native single-cell DCM stiff force: ``stiff_into(pos_d, out_d)`` runs turgor+edges
    on device buffers (no numpy round-trip) for the all-device CG. Returns the callable + handles."""
    from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM
    from ffn_sim.warp_port.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
    from ffn_sim.warp_port.dcm_warp_hybrid import _bond_accumulate
    from ffn_sim.warp_port.dcm_warp_hybrid_multicell import _dp_from_vol, _zero_vec
    p = ResolvedDCM(subdivisions=subdiv)
    verts, edges, faces = icosphere_mesh(p.R_cell, subdiv)
    N, nf, ne = verts.shape[0], faces.shape[0], edges.shape[0]
    R0 = float(np.linalg.norm(verts, axis=1).mean()); V0 = (4/3)*np.pi*R0**3
    faces_d = wp.array(faces.astype(np.int32), dtype=wp.int32, device=device)
    fcell_d = wp.zeros(nf, dtype=wp.int32, device=device)
    edges_d = wp.array(edges.astype(np.int32), dtype=wp.int32, device=device)
    r0_d = wp.array(np.linalg.norm(verts[edges[:,0]]-verts[edges[:,1]], axis=1), dtype=wp.float64, device=device)
    Vc_d = wp.zeros(1, dtype=wp.float64, device=device); dP_d = wp.zeros(1, dtype=wp.float64, device=device)
    gamma = 6.0*np.pi*65.9*p.R_cell/N

    def stiff_into(pos_d, out_d):
        wp.launch(_zero_vec, dim=N, inputs=[out_d], device=device)
        Vc_d.zero_()
        wp.launch(dcm_volume_kernel, dim=nf, inputs=[pos_d, faces_d, fcell_d, Vc_d], device=device)
        wp.launch(_dp_from_vol, dim=1, inputs=[Vc_d, wp.float64(V0), wp.float64(p.turgor_dP0), wp.float64(7.73e5), dP_d], device=device)
        wp.launch(dcm_turgor_force_kernel, dim=nf, inputs=[pos_d, faces_d, fcell_d, dP_d, out_d], device=device)
        wp.launch(_bond_accumulate, dim=ne, inputs=[pos_d, edges_d, wp.float64(p.k_edge), r0_d, out_d], device=device)
    return stiff_into, verts, gamma, V0, faces, N


def _device_cg_demo(device="cpu"):
    """(b) validate the all-device Warp CG implicit step == the scipy implicit step (same answer),
    on the real DCM stiff force. Proves the host-transfer-free CG path is correct before wiring it."""
    stiff_into, verts, gamma, V0, faces, N = make_dcm_stiff_into(device)
    stiff_np, _, _, _, _, _ = make_dcm_stiff_force(device)   # numpy wrapper for the scipy ref
    def vol(x):
        v0,v1,v2 = x[faces[:,0]],x[faces[:,1]],x[faces[:,2]]
        return abs(float(np.einsum('ij,ij->i', v0, np.cross(v1-v0,v2-v0)).sum()/6.0))
    x0 = verts.copy(); x0[:,2] *= 0.95         # MILD perturbation (well-conditioned single step)
    dt = 8.0e-6*10; a = gamma/dt               # dt×10 (CG converges cleanly → fair device==scipy check)
    scratch = {k: wp.zeros(N, dtype=wp.vec3d, device=device) for k in ("r","p","Ap","dx","Fx","Fp","xp","Fb")}
    scratch["sca"] = wp.zeros(1, dtype=wp.float64, device=device)
    # one DEVICE implicit-Euler step (n_newton=1): b=F(x0), solve (aI+K)dx=b
    x_d = wp.array(np.ascontiguousarray(x0), dtype=wp.vec3d, device=device)
    stiff_into(x_d, scratch["Fb"])
    dx_d, iters = device_cg(stiff_into, x_d, a, scratch["Fb"], scratch, device=device)
    x_dev = x0 + dx_d.numpy()
    # one SCIPY implicit step (reference)
    x_sci, info = implicit_overdamped_step(x0, stiff_np, gamma, dt, n_newton=1)
    print(f"(b) device-CG vs scipy-CG implicit step (DCM stiff, dt=10×, N={N}):")
    print(f"  device  : V/V0={vol(x_dev)/V0:.5f}  CG iters={iters}")
    print(f"  scipy   : V/V0={vol(x_sci)/V0:.5f}  CG iters={info['cg_iters']}")
    print(f"  |x_dev − x_sci| = {np.linalg.norm(x_dev-x_sci):.2e}  (expect ~0 → device CG correct)")


def make_two_cell_contact_stiff(device="cpu", subdiv=1):
    """TWO touching cells with node-face CONTACT — the testbed that actually exercises the contact
    stiffness (rep≈2e8) that dominates the implicit conditioning. Returns numpy force_fn +
    analytic per-node diagonal-stiffness estimator (contact rep·area + edges) for a Jacobi PCG."""
    import warp as wp
    from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM
    from ffn_sim.warp_port.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
    from ffn_sim.warp_port.dcm_warp_hybrid import _bond_accumulate
    from ffn_sim.warp_port.dcm_warp_hybrid_multicell import _dp_from_vol, _zero_vec
    from ffn_sim.warp_port.dcm_neighbor_warp import (pos_to_f32, face_centroids_f32,
        cohesion_grid_kernel, contact_grid_kernel)
    p = ResolvedDCM(subdivisions=subdiv); R = p.R_cell
    v1, e1, f1 = icosphere_mesh(R, subdiv); npc = v1.shape[0]
    me = float(np.linalg.norm(v1[e1[:,0]]-v1[e1[:,1]], axis=1).mean())
    c_rep, c_adh, r_contact = 0.3*me, 0.8*me, 0.3*me
    rep, adh, A = 2.0e8, 1.0e7, 4*np.pi*R**2/npc
    verts = np.concatenate([v1 + [1.95*R,0,0], v1 + [-1.95*R,0,0]])  # two cells just touching
    faces = np.concatenate([f1, f1+npc]); edges = np.concatenate([e1, e1+npc])
    cof = np.array([0]*npc+[1]*npc, np.int64); fcell = np.array([0]*f1.shape[0]+[1]*f1.shape[0], np.int64)
    N, nf, ne = verts.shape[0], faces.shape[0], edges.shape[0]
    V0 = (4/3)*np.pi*(np.linalg.norm(v1,axis=1).mean())**3
    coh_q = c_adh; con_q = c_adh + 0.7*(3*me/2.9)
    gamma = 6*np.pi*65.9*R/npc
    pos_d=wp.array(verts,dtype=wp.vec3d,device=device); force_d=wp.zeros(N,dtype=wp.vec3d,device=device)
    cof_d=wp.array(cof.astype(np.int32),dtype=wp.int32,device=device)
    faces_d=wp.array(faces.astype(np.int32),dtype=wp.int32,device=device); fcell_d=wp.array(fcell.astype(np.int32),dtype=wp.int32,device=device)
    edges_d=wp.array(edges.astype(np.int32),dtype=wp.int32,device=device)
    r0_d=wp.array(np.linalg.norm(verts[edges[:,0]]-verts[edges[:,1]],axis=1),dtype=wp.float64,device=device)
    nf32=wp.zeros(N,dtype=wp.vec3,device=device); cf32=wp.zeros(nf,dtype=wp.vec3,device=device)
    Vc=wp.zeros(2,dtype=wp.float64,device=device); dP=wp.zeros(2,dtype=wp.float64,device=device)
    ng=wp.HashGrid(32,32,32,device=device); fg=wp.HashGrid(32,32,32,device=device)
    fc_cap = 5e-8
    def force_fn(pos_np):
        pos_d.assign(np.ascontiguousarray(pos_np.reshape(N,3)))
        wp.launch(_zero_vec,dim=N,inputs=[force_d],device=device)
        wp.launch(pos_to_f32,dim=N,inputs=[pos_d,nf32],device=device); ng.build(points=nf32,radius=coh_q)
        wp.launch(face_centroids_f32,dim=nf,inputs=[pos_d,faces_d,cf32],device=device); fg.build(points=cf32,radius=con_q)
        wp.launch(cohesion_grid_kernel,dim=N,inputs=[ng.id,nf32,pos_d,cof_d,wp.float32(coh_q),wp.float64(r_contact),wp.float64(c_adh),wp.float64(rep),wp.float64(adh),wp.float64(A),wp.float64(fc_cap),force_d],device=device)
        Vc.zero_(); wp.launch(dcm_volume_kernel,dim=nf,inputs=[pos_d,faces_d,fcell_d,Vc],device=device)
        wp.launch(_dp_from_vol,dim=2,inputs=[Vc,wp.float64(V0),wp.float64(p.turgor_dP0),wp.float64(7.73e5),dP],device=device)
        wp.launch(dcm_turgor_force_kernel,dim=nf,inputs=[pos_d,faces_d,fcell_d,dP,force_d],device=device)
        wp.launch(contact_grid_kernel,dim=N,inputs=[fg.id,nf32,pos_d,cof_d,faces_d,fcell_d,wp.float32(con_q),wp.float64(rep),wp.float64(adh),wp.float64(c_rep),wp.float64(c_adh),force_d],device=device)
        wp.launch(_bond_accumulate,dim=ne,inputs=[pos_d,edges_d,wp.float64(p.k_edge),r0_d,force_d],device=device)
        wp.synchronize_device(device); return force_d.numpy().astype(np.float64)
    return force_fn, verts, gamma, V0, faces, N


def _precond_contact_demo(device="cpu"):
    """Does the all-important CONTACT stiffness wreck CG conditioning, and does a preconditioner
    help? Compress two touching cells, take ONE implicit step at dt×100, count CG iters with
    NO preconditioner vs Jacobi(Hutchinson)."""
    force_fn, verts, gamma, V0, faces, N = make_two_cell_contact_stiff(device)
    x0 = verts.copy(); x0[:N//2,0] -= 0.3e-6; x0[N//2:,0] += 0.3e-6   # push the two cells together
    dt = 8e-6*100
    xa, ia = implicit_overdamped_step(x0, force_fn, gamma, dt, n_newton=1, cg_maxiter=600, precond=False)
    xb, ib = implicit_overdamped_step(x0, force_fn, gamma, dt, n_newton=1, cg_maxiter=600, precond=True)
    print(f"two-cell CONTACT (rep=2e8), 1 implicit step @dt×100 (N={N}):")
    print(f"  NO preconditioner   : CG iters = {ia['cg_iters']}   |Δx|={np.linalg.norm(xa-x0):.2e}")
    print(f"  Jacobi (Hutchinson) : CG iters = {ib['cg_iters']}   |Δx|={np.linalg.norm(xb-x0):.2e}")
    print(f"  → confirms whether CONTACT stiffness dominates conditioning + if Jacobi helps")


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
    for nN in (1, 6):
        tag = "linearly-implicit (n_newton=1)" if nN == 1 else "NEWTON (n_newton=6)"
        print(f"-- {tag} --")
        print(f"{'dt/dt_e':>8} {'steps':>7} {'V/V0':>7} {'|x-ref|':>9} {'cgit/step':>9} {'newton/step':>11}")
        for mult in (100, 1000, 10000):
            dt = dt_e*mult; ns = max(1, int(T/dt)); x = x0.copy(); tot_it=0; tot_nw=0; ok=True
            for _ in range(ns):
                x, info = implicit_overdamped_step(x, stiff_force, gamma, dt, n_newton=nN)
                tot_it += info["cg_iters"]; tot_nw += info["newton_iters"]
                if not np.isfinite(x).all() or vol(x) > 50*V0: ok=False; break
            e = float(np.linalg.norm(x-xref)) if ok else float('inf')
            print(f"{mult:>8} {ns:>7} {vol(x)/V0 if ok else 0:>7.4f} {e:>9.2e} "
                  f"{tot_it/max(ns,1):>9.1f} {tot_nw/max(ns,1):>11.1f}")


if __name__ == "__main__":
    print("== I1 spring solver =="); _spring_demo()
    print("== I2 DCM stiff force =="); _dcm_stiff_demo()
    print("== I3 dt ramp =="); _dt_ramp()
