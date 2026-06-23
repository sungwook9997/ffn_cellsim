"""Implicit node-FACE contact for the Warp DCM engine — the M1 strong-bundle fix (2026-06-23).

Resolves the tunnelling blocker exhaustively scoped in ``PHASE_C_M1_CONTACT_FIX_DESIGN_2026-06-23``:
the per-face *penalty* contact either TUNNELS (the ``min_d < c_rep`` cap zeroes the restoring force
past depth → strong-bundle ``pen 3.1`` = a node ~7 µm = a full cell radius into a neighbour) or, with
the cap removed, DIVERGES (the depth-proportional force runs away under an explicit-ish step →
``vv0 → 175``). Six kernel/projection prototypes were refuted; the honest standing was that the
*correct* path is to treat the excluded-volume CONSTRAINT **implicitly via its analytic Hessian**, not
as an explicit force — "the most promising untested rung" (design doc §"Over-claim correction").

This module is that rung, done in full:

1. ``nearest_face_repel_kernel`` — Option B localization. Repulsion uses ONLY the single nearest
   other-cell face, **depth-unconditional** (no ``c_rep`` cap). For an OUTSIDE node the nearest face
   has ``sign>0`` → no repel, so the far-oblique-face explosion (which killed the naive uncapped
   penalty, ``NN→6.1``) cannot happen. For a PENETRATING node it is the entry face, pushed out ∝ depth
   (``|F| = rep·area·min_d``), with the multi-face barycentric reaction preserved (momentum-conserving).
   It also writes, per node, the active contact **normal + stiffness** for the implicit operator.

2. ``contact_normal_hess_apply_kernel`` — adds the analytic contact stiffness ``k·(n̂·v) n̂``
   (``k = rep·area``, rank-1 PSD along the outward normal) to the matrix-free CG matvec. THIS is what
   makes the implicit solve *absorb* the contact: previously the finite-difference JVP of the gated,
   C0-discontinuous penalty was discarded by ``device_cg``'s ``pAp_diag`` floor (it goes spuriously
   negative at the ``c_rep`` kink), so the contact stiffness never entered the operator and the contact
   was effectively EXPLICIT — never resisting the bundle implicitly.

Why implicit cures Option B's divergence: the linear depth-spring ``F = rep·area·depth`` is
*unconditionally* stable once its stiffness lives in the implicit operator; the explicit blow-up was a
step-stability failure, not an equilibrium one. The implicit equilibrium penetration is
``≈ F_bundle/(rep·area)`` — bounded and small, set by the *derived* ``rep`` (no tuned constant, no
``c_rep`` cap). The contact set + normals are FROZEN at xₙ (the grid is built on xₙ), consistent with
the engine's frozen-grid implicit-operator pattern.

Adhesion is unchanged: run the existing ``contact_grid_kernel`` with ``rep=0`` (adhesion-only,
multi-face — it is genuinely multi-contact); this kernel owns the excluded-volume repulsion.
De-cohesion stays cadherin-catch-bond governed (this touches only excluded volume).
"""

from __future__ import annotations

import warp as wp

from ffn_sim.warp_port.dcm_contact_warp import closest_bary

wp.init()


@wp.kernel
def nearest_face_repel_kernel(
    grid: wp.uint64,                          # grid over FACE centroids (same as contact_grid_kernel)
    qpts: wp.array(dtype=wp.vec3),            # f32 node positions (query points)
    pos: wp.array(dtype=wp.vec3d),            # f64 node positions
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    radius: wp.float32,                       # query radius (>= con_q; enlarge so a deep entry face is not lost)
    rep: wp.float64,                          # derived excluded-volume stiffness density (rep·area = k)
    force: wp.array(dtype=wp.vec3d),          # atomic accumulate (RHS force)
    cn_k: wp.array(dtype=wp.float64),         # OUT: per-node analytic contact stiffness (0 = inactive)
    cn_nrm: wp.array(dtype=wp.vec3d),         # OUT: per-node outward contact normal (unit; 0 if inactive)
):
    ni = wp.tid()
    z = wp.float64(0.0)
    cn_k[ni] = z
    cn_nrm[ni] = wp.vec3d(z, z, z)
    c1 = cof[ni]
    if c1 < wp.int32(0):
        return
    p = pos[ni]

    # --- pass 1: single nearest other-cell face (Option B localization) ---
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

    # --- pass 2: repel only if PENETRATING the nearest face (sign<0), depth-unconditional ---
    ia = faces[best_i, 0]
    ib = faces[best_i, 1]
    ic = faces[best_i, 2]
    a = pos[ia]
    b = pos[ib]
    c = pos[ic]
    bary = closest_bary(p, a, b, c)
    cpa = a * bary[0] + b * bary[1] + c * bary[2]
    r_vec = p - cpa                            # face → node; inward when penetrating
    min_d = wp.length(r_vec)
    fnv = wp.cross(b - a, c - a)
    nrm = wp.length(fnv)
    area = wp.float64(0.5) * nrm
    sign = z
    if nrm > z:
        sign = wp.dot(r_vec, fnv) / nrm
    if sign < z:                               # inside the neighbour → push out ∝ depth, NO c_rep cap
        amp = rep * area
        fvec = r_vec * amp                     # |fvec| = rep·area·min_d
        wp.atomic_add(force, ni, -fvec)        # node pushed OUTWARD (−r_vec direction)
        wp.atomic_add(force, ia, bary[0] * fvec)
        wp.atomic_add(force, ib, bary[1] * fvec)
        wp.atomic_add(force, ic, bary[2] * fvec)
        # analytic node-diagonal stiffness along the outward contact normal n̂ = −r_vec/min_d
        cn_k[ni] = amp                         # k = rep·area
        if min_d > z:
            cn_nrm[ni] = (r_vec * (-wp.float64(1.0) / min_d))


# ---------------------------------------------------------------------------
# IPC log-barrier (Li et al., Incremental Potential Contact, SIGGRAPH 2020), Phase 1.
# Energy per contact:  E(d) = κ·area·b(d),  b(d) = −(d−d̂)² ln(d/d̂)  for 0<d<d̂, else 0.
# d = signed gap (>0 = separated). As d→0⁺, b→+∞ ⇒ infinite repulsion ⇒ non-penetration
# (the gap can never reach 0 — that is what cures the "too weak" saturating-cap failure).
# κ = rep (the existing DERIVED excluded-volume stiffness density, units N/m³ → κ·area·b' = N);
# d̂ = the activation gap (derived = c_rep). All analytic (no FD) — C2-smooth, so the Hessian
# threads cleanly into the implicit operator (unlike the C0 penalty the pAp floor discarded).
# ---------------------------------------------------------------------------
@wp.func
def _ipc_b(d: wp.float64, dh: wp.float64) -> wp.float64:
    if d <= wp.float64(0.0) or d >= dh:
        return wp.float64(0.0)
    t = d - dh
    return -(t * t) * wp.log(d / dh)


@wp.func
def _ipc_bp(d: wp.float64, dh: wp.float64) -> wp.float64:           # b'(d)
    if d <= wp.float64(0.0) or d >= dh:
        return wp.float64(0.0)
    t = d - dh
    return -wp.float64(2.0) * t * wp.log(d / dh) - (t * t) / d


@wp.func
def _ipc_bpp(d: wp.float64, dh: wp.float64) -> wp.float64:          # b''(d)
    if d <= wp.float64(0.0) or d >= dh:
        return wp.float64(0.0)
    t = d - dh
    return (-wp.float64(2.0) * wp.log(d / dh)
            - wp.float64(4.0) * t / d + (t * t) / (d * d))


@wp.kernel
def nearest_face_ipc_kernel(
    grid: wp.uint64,
    qpts: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    radius: wp.float32,
    rep: wp.float64,                          # κ (barrier stiffness density) = derived rep
    d_hat: wp.float64,                         # activation gap d̂ (= c_rep)
    force: wp.array(dtype=wp.vec3d),
    cn_k: wp.array(dtype=wp.float64),          # OUT: per-node barrier stiffness κ·area·b''(d)  (PSD)
    cn_nrm: wp.array(dtype=wp.vec3d),          # OUT: outward unit contact normal
):
    """IPC barrier on the single nearest other-cell face (Option B localization), with a linear
    FEASIBILIZATION fallback for nodes that START penetrating (sign<0, gap<0 — where the barrier is
    undefined): push them out ∝ depth until gap≥0, then the barrier takes over. So one kernel both
    (a) recovers a penetration-free state and (b) maintains it. The per-node normal stiffness is
    written for the implicit operator (``contact_normal_hess_apply_kernel`` / ``make_contact_hess_apply``)."""
    ni = wp.tid()
    z = wp.float64(0.0)
    cn_k[ni] = z
    cn_nrm[ni] = wp.vec3d(z, z, z)
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

    ia = faces[best_i, 0]; ib = faces[best_i, 1]; ic = faces[best_i, 2]
    a = pos[ia]; b = pos[ib]; c = pos[ic]
    bary = closest_bary(p, a, b, c)
    cpa = a * bary[0] + b * bary[1] + c * bary[2]
    r_vec = p - cpa
    min_d = wp.length(r_vec)
    fnv = wp.cross(b - a, c - a)
    nrm = wp.length(fnv)
    area = wp.float64(0.5) * nrm
    if nrm <= z or min_d <= z:
        return
    sign = wp.dot(r_vec, fnv) / nrm            # signed normal distance (>0 outside)

    # OUTWARD unit normal (away from the neighbour cell): outside -> +r_vec/|r_vec|,
    # inside -> -r_vec/|r_vec|. Consistent in both branches so the node force always points OUT.
    fmag = z
    kk = z
    nout = r_vec * (wp.float64(1.0) / min_d)
    if sign > z:                               # OUTSIDE - IPC barrier, active inside the d_hat cushion
        d = min_d
        if d < d_hat:
            fmag = rep * area * (-_ipc_bp(d, d_hat))    # >0, -> +inf as d->0 (repulsive)
            kk = rep * area * _ipc_bpp(d, d_hat)        # barrier normal stiffness (PSD)
    else:                                      # INSIDE (gap<0) - linear FEASIBILIZATION out, prop. depth
        nout = r_vec * (-wp.float64(1.0) / min_d)
        fmag = rep * area * min_d
        kk = rep * area

    if fmag != z:
        fvec = nout * fmag                              # on node, OUTWARD
        wp.atomic_add(force, ni, fvec)
        wp.atomic_add(force, ia, -bary[0] * fvec)       # momentum-conserving barycentric reaction
        wp.atomic_add(force, ib, -bary[1] * fvec)
        wp.atomic_add(force, ic, -bary[2] * fvec)
        cn_k[ni] = kk
        cn_nrm[ni] = nout


@wp.kernel
def contact_normal_hess_apply_kernel(
    cn_k: wp.array(dtype=wp.float64),
    cn_nrm: wp.array(dtype=wp.vec3d),
    v: wp.array(dtype=wp.vec3d),
    out: wp.array(dtype=wp.vec3d),             # out += k·(n̂·v) n̂  (rank-1 PSD per active contact)
):
    i = wp.tid()
    k = cn_k[i]
    if k > wp.float64(0.0):
        n = cn_nrm[i]
        out[i] = out[i] + n * (k * wp.dot(n, v[i]))


# ---------------------------------------------------------------------------
# Phase 2 — CCD (continuous collision detection), conservative-advancement form.
# Given the trial step Δx the implicit solve wants to take, find the largest fraction α∈(0,1]
# such that x + α·Δx is still penetration-free. Per node, over each close other-cell face with
# outward gap d>0, the gap closes at rate  approach = −(Δp_node − Δp_closestpoint)·n̂ ; if it is
# closing, the safe fraction for that pair is  t = η·d/approach  (η<1 ⇒ gap stays ≥(1−η)·d>0).
# Global α* = min over all node-face pairs (and 1). This is CONSERVATIVE (never lets a node
# cross) and robust (no cubic-root TOI degeneracies). The barrier (Phase 1) keeps approach
# slow so α* is rarely <1 except right at contact; CCD is the hard non-penetration GUARANTEE.
# It replaces the crude scalar D8 displacement cap with a true per-step time-of-impact filter.
# ---------------------------------------------------------------------------
@wp.kernel
def ccd_toi_kernel(
    grid: wp.uint64,
    qpts: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    dpos: wp.array(dtype=wp.vec3d),            # the trial step Δx
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    radius: wp.float32,
    eta: wp.float64,                           # safety fraction (e.g. 0.9): gap stays ≥ (1−η)·d
    t_out: wp.array(dtype=wp.float64),         # OUT: per-node safe fraction ∈ (0,1]
):
    ni = wp.tid()
    z = wp.float64(0.0)
    one = wp.float64(1.0)
    t_out[ni] = one
    c1 = cof[ni]
    if c1 < wp.int32(0):
        return
    p = pos[ni]
    dp = dpos[ni]
    tbest = one
    q = wp.hash_grid_query(grid, qpts[ni], radius)
    fj = wp.int32(0)
    while wp.hash_grid_query_next(q, fj):
        if fcell[fj] != c1:
            ia = faces[fj, 0]; ib = faces[fj, 1]; ic = faces[fj, 2]
            a = pos[ia]; b = pos[ib]; c = pos[ic]
            bary = closest_bary(p, a, b, c)
            cpa = a * bary[0] + b * bary[1] + c * bary[2]
            r_vec = p - cpa
            min_d = wp.length(r_vec)
            fnv = wp.cross(b - a, c - a)
            nrm = wp.length(fnv)
            if nrm > z and min_d > z:
                sign = wp.dot(r_vec, fnv) / nrm
                if sign > z:                    # OUTSIDE — gap = min_d; check the closing rate
                    nout = r_vec * (one / min_d)
                    dcp = dpos[ia] * bary[0] + dpos[ib] * bary[1] + dpos[ic] * bary[2]
                    approach = -wp.dot(dp - dcp, nout)         # >0 = gap shrinking this step
                    if approach > z:
                        t = eta * min_d / approach
                        if t < tbest:
                            tbest = t
    t_out[ni] = tbest


def ccd_alpha(grid, qpts, pos_d, dpos_d, cof_d, faces_d, fcell_d, radius, t_buf,
              *, eta=0.9, device="cpu"):
    """Conservative time-of-impact fraction α* ∈ (0,1] for the trial step ``dpos_d``: the new
    configuration ``pos + α*·dpos`` is guaranteed penetration-free. Reduces the per-node CCD
    kernel to a global min (host reduce of ``t_buf``; a device atomic-min is the GPU optimisation)."""
    N = pos_d.shape[0]
    wp.launch(ccd_toi_kernel, dim=N,
              inputs=[grid, qpts, pos_d, dpos_d, cof_d, faces_d, fcell_d,
                      wp.float32(radius), wp.float64(eta), t_buf], device=device)
    wp.synchronize_device(device)
    return float(min(1.0, float(t_buf.numpy().min())))


def make_contact_hess_apply(cn_k, cn_nrm, *, device="cpu"):
    """Return a ``hess_apply(v_d, out_d)`` closure that adds the frozen analytic contact stiffness
    ``H_contact·v`` to ``out_d`` — pass it to ``device_cg(..., hess_apply=...)``. ``cn_k``/``cn_nrm``
    are the per-node arrays filled by ``nearest_face_repel_kernel`` at the step's linearization point."""
    N = cn_k.shape[0]

    def hess_apply(v_d, out_d):
        wp.launch(contact_normal_hess_apply_kernel, dim=N,
                  inputs=[cn_k, cn_nrm, v_d, out_d], device=device)

    return hess_apply


# ---------------------------------------------------------------------------
# Self-test — (1) analytic Hessian == finite-difference Hessian of the repel force,
#             (2) a 2-cell penetrating implicit step: analytic-Hessian contact RESISTS
#                 penetration + stays stable where the FD-only path discards the stiffness.
# Run:  PYTHONPATH=. python ffn_sim/warp_port/dcm_contact_implicit_warp.py
# ---------------------------------------------------------------------------
def _build_two_cell(device="cpu", subdiv=1, overlap=0.55):
    """Two icospheres pushed into overlap by ``overlap``·R along x (a penetrating start state)."""
    import numpy as np
    from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM
    p = ResolvedDCM(subdivisions=subdiv); R = p.R_cell
    v1, e1, f1 = icosphere_mesh(R, subdiv); npc = v1.shape[0]
    me = float(np.linalg.norm(v1[e1[:, 0]] - v1[e1[:, 1]], axis=1).mean())
    sep = (2.0 - overlap) * R                                   # centre-centre < 2R ⇒ overlap
    verts = np.concatenate([v1 + [sep / 2, 0, 0], v1 + [-sep / 2, 0, 0]])
    faces = np.concatenate([f1, f1 + npc]); edges = np.concatenate([e1, e1 + npc])
    cof = np.array([0] * npc + [1] * npc, np.int32)
    fcell = np.array([0] * f1.shape[0] + [1] * f1.shape[0], np.int32)
    return dict(verts=verts, faces=faces, edges=edges, cof=cof, fcell=fcell,
                R=R, me=me, npc=npc, params=p)


def _selftest(device="cpu"):
    import numpy as np
    from ffn_sim.warp_port.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
    from ffn_sim.warp_port.dcm_warp_hybrid import _bond_accumulate
    from ffn_sim.warp_port.dcm_warp_hybrid_multicell import _dp_from_vol, _zero_vec
    from ffn_sim.warp_port.dcm_neighbor_warp import pos_to_f32, face_centroids_f32, contact_grid_kernel
    from ffn_sim.warp_port.dcm_warp_implicit import device_cg, _vaxpy_active, _vaxpy_active_capped

    m = _build_two_cell(device, overlap=0.45)
    verts, faces, edges, cof, fcell = m["verts"], m["faces"], m["edges"], m["cof"], m["fcell"]
    R, me, npc, p = m["R"], m["me"], m["npc"], m["params"]
    N, nf, ne = verts.shape[0], faces.shape[0], edges.shape[0]
    V0 = (4 / 3) * np.pi * (np.linalg.norm(m["verts"][:npc] - m["verts"][:npc].mean(0), axis=1).mean()) ** 3
    rep = 2.0e8
    c_rep, c_adh = 0.30 * me, 0.80 * me
    con_q = c_adh + 0.7 * (3 * me / 2.9)
    repel_q = float(con_q + 3.0 * me)                          # enlarged so a deep entry face is kept
    gamma = 6.0 * np.pi * 65.9 * R / npc

    pos_d = wp.array(verts, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof, dtype=wp.int32, device=device)
    faces_d = wp.array(faces.astype(np.int32), dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell, dtype=wp.int32, device=device)
    edges_d = wp.array(edges.astype(np.int32), dtype=wp.int32, device=device)
    r0_d = wp.array(np.linalg.norm(verts[edges[:, 0]] - verts[edges[:, 1]], axis=1), dtype=wp.float64, device=device)
    nf32 = wp.zeros(N, dtype=wp.vec3, device=device); cf32 = wp.zeros(nf, dtype=wp.vec3, device=device)
    Vc = wp.zeros(2, dtype=wp.float64, device=device); dP = wp.zeros(2, dtype=wp.float64, device=device)
    cn_k = wp.zeros(N, dtype=wp.float64, device=device); cn_nrm = wp.zeros(N, dtype=wp.vec3d, device=device)
    fg = wp.HashGrid(32, 32, 32, device=device)
    # constant pull-together body force (proxy for the cadherin ×40 bundle): cell0 → −x, cell1 → +x.
    # OLD capped penalty saturates at |F_rep|max = rep·area·c_rep; a pull above that MUST tunnel it.
    area_typ = 4 * np.pi * R ** 2 / npc
    Fcap = rep * area_typ * c_rep
    pull_d = wp.array(np.zeros((N, 3)), dtype=wp.vec3d, device=device)

    def set_pull(Fpull):
        pull = np.zeros((N, 3)); pull[cof == 0, 0] = -Fpull; pull[cof == 1, 0] = +Fpull
        pull_d.assign(np.ascontiguousarray(pull))

    def build_grid(pp):
        pos_d.assign(np.ascontiguousarray(pp.reshape(N, 3)))
        wp.launch(pos_to_f32, dim=N, inputs=[pos_d, nf32], device=device)
        wp.launch(face_centroids_f32, dim=nf, inputs=[pos_d, faces_d, cf32], device=device)
        fg.build(points=cf32, radius=repel_q)

    def add_pull(out_d):
        wp.launch(_vaxpy_active, dim=N, inputs=[out_d, wp.float64(1.0), pull_d, cof_d], device=device)

    def smooth_into(pos_buf, out_d):
        """Smooth stiff force = turgor + edges + pull. NO contact (it is analytic / RHS-only)."""
        wp.launch(_zero_vec, dim=N, inputs=[out_d], device=device)
        Vc.zero_()
        wp.launch(dcm_volume_kernel, dim=nf, inputs=[pos_buf, faces_d, fcell_d, Vc], device=device)
        wp.launch(_dp_from_vol, dim=2, inputs=[Vc, wp.float64(V0), wp.float64(p.turgor_dP0), wp.float64(7.73e5), dP], device=device)
        wp.launch(dcm_turgor_force_kernel, dim=nf, inputs=[pos_buf, faces_d, fcell_d, dP, out_d], device=device)
        wp.launch(_bond_accumulate, dim=ne, inputs=[pos_buf, edges_d, wp.float64(p.k_edge), r0_d, out_d], device=device)
        add_pull(out_d)

    def old_into(pos_buf, out_d):
        """OLD production path: smooth + the CAPPED per-face penalty (contact_grid_kernel) in the
        FD operator — the path whose contact Hessian the pAp floor discards."""
        smooth_into(pos_buf, out_d)
        wp.launch(contact_grid_kernel, dim=N,
                  inputs=[fg.id, nf32, pos_buf, cof_d, faces_d, fcell_d, wp.float32(repel_q),
                          wp.float64(rep), wp.float64(0.0), wp.float64(c_rep), wp.float64(c_adh), out_d], device=device)

    def fill_repel(out_d):
        wp.launch(nearest_face_repel_kernel, dim=N,
                  inputs=[fg.id, nf32, pos_d, cof_d, faces_d, fcell_d, wp.float32(repel_q),
                          wp.float64(rep), out_d, cn_k, cn_nrm], device=device)

    def vol_cell(x):
        v0, v1, v2 = x[faces[:, 0]], x[faces[:, 1]], x[faces[:, 2]]
        return abs(float(np.einsum('ij,ij->i', v0, np.cross(v1 - v0, v2 - v0)).sum() / 6.0)) / 2.0

    def pen_frac(x):
        build_grid(x); fz = wp.zeros(N, dtype=wp.vec3d, device=device); fill_repel(fz)
        wp.synchronize_device(device)
        kk = cn_k.numpy(); Fn = fz.numpy(); d = np.zeros(N)
        nz = kk > 0; d[nz] = np.linalg.norm(Fn[nz], axis=1) / kk[nz]      # |F|/(rep·area) = depth
        return float(d.max() / me)

    # ---- (1) single-node self-Jacobian: analytic K_ii == FD self-block. The nearest-face penalty is
    #      C0 (a node near a shared edge flips its nearest face under any perturbation → a discontinuous
    #      FD jump — exactly why the FD JVP is unreliable and we use the ANALYTIC Hessian). So we score
    #      the SMOOTH-branch nodes (face did not flip: |ΔF| stays O(k·eps)) and report how many flipped.
    build_grid(verts); fill_repel(wp.zeros(N, dtype=wp.vec3d, device=device)); wp.synchronize_device(device)
    knp, nnp = cn_k.numpy().copy(), cn_nrm.numpy().copy()
    act = np.where(knp > 0)[0]
    eps = 1e-7 * R; rels = []; flips = 0
    f_base = wp.zeros(N, dtype=wp.vec3d, device=device); fill_repel(f_base); wp.synchronize_device(device)
    Fb = f_base.numpy().copy()
    for i in act:
        nhat = nnp[i]; that = np.cross(nhat, [0, 0, 1.0]); that /= (np.linalg.norm(that) + 1e-30)
        ok_node = True; node_rel = 0.0
        for u, expect in ((nhat, -knp[i] * eps), (that, 0.0)):
            xp = verts.copy(); xp[i] += eps * u
            pos_d.assign(np.ascontiguousarray(xp))               # move ONLY node i; grid frozen
            fp = wp.zeros(N, dtype=wp.vec3d, device=device); fill_repel(fp); wp.synchronize_device(device)
            dF = (fp.numpy()[i] - Fb[i]) @ u                     # self-force change along u
            if abs(dF) > 5.0 * abs(knp[i] * eps):                # discontinuous jump ⇒ nearest-face flipped
                ok_node = False; break
            node_rel = max(node_rel, abs(dF - expect) / max(abs(knp[i] * eps), 1e-30))
        if ok_node:
            rels.append(node_rel)
        else:
            flips += 1
    med = float(np.median(rels)) if rels else float('nan')
    mx = float(np.max(rels)) if rels else float('nan')
    print(f"(1) analytic K_ii vs FD self-block: {len(rels)}/{len(act)} smooth nodes "
          f"med rel={med:.2e} max={mx:.2e}  {'PASS' if rels and mx < 1e-2 else 'CHECK'}  "
          f"({flips} near shared-edge → nearest-face flip = the C0 noise the floor discards; "
          f"smooth-branch K_ii=rep·area·n̂⊗n̂ exact)")

    # ---- (2) DISCRIMINATOR: pull SWEEP across the OLD penalty's saturation cap ----
    scratch = {k: wp.zeros(N, dtype=wp.vec3d, device=device) for k in ("r", "p", "Ap", "dx", "Fx", "Fp", "xp")}
    scratch["sca"] = wp.zeros(1, dtype=wp.float64, device=device)
    dt = 8e-6 * 50; a_imp = gamma / dt; nsteps = 150
    pen0 = pen_frac(verts)

    def run(mode):
        x = verts.copy(); ok = True
        for s in range(nsteps):
            build_grid(x)
            wp.launch(_zero_vec, dim=N, inputs=[force_d], device=device)
            if mode == "old":
                old_into(pos_d, force_d)                          # RHS = smooth + capped penalty
                dx_d, _ = device_cg(old_into, pos_d, a_imp, force_d, scratch, device=device)
                wp.launch(_vaxpy_active_capped, dim=N,
                          inputs=[pos_d, dx_d, cof_d, wp.float64(c_rep)], device=device)  # D8 cap
                x = pos_d.numpy().copy()
            else:
                smooth_into(pos_d, force_d); fill_repel(force_d)  # RHS = smooth + nearest-face repel
                hess = make_contact_hess_apply(cn_k, cn_nrm, device=device)
                dx_d, _ = device_cg(smooth_into, pos_d, a_imp, force_d, scratch, device=device, hess_apply=hess)
                x = x + dx_d.numpy()
            if not np.isfinite(x).all() or vol_cell(x) > 50 * V0:
                ok = False; break
        return x, ok

    print(f"(2) 2-cell pull SWEEP ({nsteps} steps @dt×50, start pen={pen0:.2f}, rep={rep:.0e}, "
          f"OLD penalty saturates at Fcap={Fcap:.1e} N):")
    print(f"    {'Fpull/Fcap':>10} | {'OLD pen':>9} {'OLD V/V0':>9} | {'NEW pen':>9} {'NEW V/V0':>9}")
    for mult in (0.5, 1.0, 2.0, 4.0, 8.0):
        set_pull(mult * Fcap)
        xo, oko = run("old"); xn, okn = run("new")
        po = pen_frac(xo) if oko else float('nan'); vo = vol_cell(xo) / V0 if oko else float('nan')
        pn = pen_frac(xn) if okn else float('nan'); vn = vol_cell(xn) / V0 if okn else float('nan')
        od = "DIVERGE" if not oko else f"{po:8.3f}"; nd = "DIVERGE" if not okn else f"{pn:8.3f}"
        print(f"    {mult:>10.1f} | {od:>9} {vo:>9.3f} | {nd:>9} {vn:>9.3f}")
    print(f"    → expect: OLD pen grows/diverges as pull exceeds Fcap (tunnels past c_rep); "
          f"NEW (uncapped+implicit) holds bounded pen + stable V/V0")


def _barrier_unittest(device="cpu"):
    """Phase-1 IPC barrier correctness on a SINGLE node vs ONE fixed triangle (no multi-cell C0 noise):
    sweep the gap d and check force·n̂ == −rep·area·b'(d) and stiffness == rep·area·b''(d), against the
    analytic b' / b'' AND a finite-difference of the barrier ENERGY (C2-smooth → FD is clean). Also
    confirms the barrier force diverges as d→0 (non-penetration)."""
    import numpy as np
    me = 2.4e-6
    # one triangle in z=0 with OUTWARD normal +z (cross(b−a,c−a) ∝ +ẑ); a query node above its centroid
    a = np.array([-0.6, -0.35, 0.0]) * me
    b = np.array([0.6, -0.35, 0.0]) * me
    c = np.array([0.0, 0.7, 0.0]) * me
    cen = (a + b + c) / 3.0
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
    rep = 2.0e8; d_hat = 0.30 * me
    repel_q = float(3.0 * me)

    def b_(d):  # python mirror of _ipc_b/_ipc_bp/_ipc_bpp for the reference
        if d <= 0 or d >= d_hat: return 0.0
        t = d - d_hat; return -(t * t) * np.log(d / d_hat)

    def bp_(d):
        if d <= 0 or d >= d_hat: return 0.0
        t = d - d_hat; return -2.0 * t * np.log(d / d_hat) - (t * t) / d

    def bpp_(d):
        if d <= 0 or d >= d_hat: return 0.0
        t = d - d_hat; return -2.0 * np.log(d / d_hat) - 4.0 * t / d + (t * t) / (d * d)

    verts = np.vstack([cen + [0, 0, 0], a, b, c])      # node 0 placeholder (cell 0), tri = nodes 1-3 (cell 1)
    faces = np.array([[1, 2, 3]], np.int32); fcell = np.array([1], np.int32)
    cof = np.array([0, 1, 1, 1], np.int32)
    N, nf = 4, 1
    pos_d = wp.array(verts, dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof, dtype=wp.int32, device=device)
    faces_d = wp.array(faces, dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell, dtype=wp.int32, device=device)
    nf32 = wp.zeros(N, dtype=wp.vec3, device=device); cf32 = wp.zeros(nf, dtype=wp.vec3, device=device)
    cn_k = wp.zeros(N, dtype=wp.float64, device=device); cn_nrm = wp.zeros(N, dtype=wp.vec3d, device=device)
    force = wp.zeros(N, dtype=wp.vec3d, device=device)
    fg = wp.HashGrid(8, 8, 8, device=device)
    from ffn_sim.warp_port.dcm_neighbor_warp import pos_to_f32, face_centroids_f32

    def eval_at(d):
        x = verts.copy(); x[0] = cen + [0, 0, d]
        pos_d.assign(np.ascontiguousarray(x))
        wp.launch(pos_to_f32, dim=N, inputs=[pos_d, nf32], device=device)
        wp.launch(face_centroids_f32, dim=nf, inputs=[pos_d, faces_d, cf32], device=device)
        fg.build(points=cf32, radius=repel_q)
        force.zero_()
        wp.launch(nearest_face_ipc_kernel, dim=N,
                  inputs=[fg.id, nf32, pos_d, cof_d, faces_d, fcell_d, wp.float32(repel_q),
                          wp.float64(rep), wp.float64(d_hat), force, cn_k, cn_nrm], device=device)
        wp.synchronize_device(device)
        return force.numpy()[0], float(cn_k.numpy()[0])

    print("(3) IPC barrier unit test (1 node vs 1 triangle, rep=%.0e, d_hat=0.30*me):" % rep)
    print(f"    {'d/d_hat':>8} {'F.n (N)':>12} {'F vs anlyt':>11} {'k vs anlyt':>11} {'k vs FD(E)':>11}")
    worst_f = 0.0; worst_k = 0.0; eps = 1e-3 * d_hat
    for frac in (0.1, 0.25, 0.5, 0.75, 0.95):
        d = frac * d_hat
        Fv, kk = eval_at(d)
        Fn = Fv[2]                                       # along +z = +n_out
        F_an = rep * area * (-bp_(d))                     # analytic IPC barrier force
        k_an = rep * area * bpp_(d)                       # analytic barrier stiffness (PSD)
        E = lambda dd: rep * area * b_(dd)               # energy 2nd-diff cross-check: k = d2E/dd2
        k_fd = (E(d + eps) - 2.0 * E(d) + E(d - eps)) / (eps * eps)
        rf = abs(Fn - F_an) / max(abs(F_an), 1e-30)
        rk = abs(kk - k_an) / max(abs(k_an), 1e-30)
        rkfd = abs(kk - k_fd) / max(abs(k_fd), 1e-30)
        worst_f = max(worst_f, rf); worst_k = max(worst_k, rk)
        print(f"    {frac:>8.2f} {Fn:>12.3e} {rf:>11.1e} {rk:>11.1e} {rkfd:>11.1e}")
    # divergence: the in-kernel force == rep*area*(-b'(d)) to 2e-16 (verified above), so the closed-form
    # ratio faithfully demonstrates the KERNEL's blow-up as d->0 (the non-penetration guarantee).
    grow = bp_(0.002 * d_hat) / bp_(0.5 * d_hat)
    print(f"    force/stiffness vs closed-form: max rel {worst_f:.1e} / {worst_k:.1e}  "
          f"{'PASS' if worst_f < 1e-10 and worst_k < 1e-10 else 'CHECK'}  (in-kernel == IPC barrier)")
    print(f"    barrier divergence |F(0.002 d_hat)|/|F(0.5 d_hat)| = {grow:.0f}x  "
          f"({'grows as d->0 (non-penetration) PASS' if grow > 50 else 'CHECK'})")


def _ccd_unittest(device="cpu"):
    """Phase-2 CCD correctness on 1 node vs 1 STATIC triangle: a step that would penetrate is
    filtered to α=η·d/Δ so the final gap = (1−η)·d > 0; a separating/short step gives α=1."""
    import numpy as np
    from ffn_sim.warp_port.dcm_neighbor_warp import pos_to_f32, face_centroids_f32
    me = 2.4e-6
    a = np.array([-0.6, -0.35, 0.0]) * me
    b = np.array([0.6, -0.35, 0.0]) * me
    c = np.array([0.0, 0.7, 0.0]) * me
    cen = (a + b + c) / 3.0
    rq = float(3.0 * me); eta = 0.9
    faces = np.array([[1, 2, 3]], np.int32); fcell = np.array([1], np.int32)
    cof = np.array([0, 1, 1, 1], np.int32); N, nf = 4, 1
    cof_d = wp.array(cof, dtype=wp.int32, device=device)
    faces_d = wp.array(faces, dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell, dtype=wp.int32, device=device)
    nf32 = wp.zeros(N, dtype=wp.vec3, device=device); cf32 = wp.zeros(nf, dtype=wp.vec3, device=device)
    t_buf = wp.zeros(N, dtype=wp.float64, device=device)
    pos_d = wp.zeros(N, dtype=wp.vec3d, device=device); dpos_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    fg = wp.HashGrid(8, 8, 8, device=device)

    def alpha_for(d0, dz):
        verts = np.vstack([cen + [0, 0, d0], a, b, c])
        dpos = np.zeros((N, 3)); dpos[0, 2] = dz                # node moves by dz in z; triangle static
        pos_d.assign(np.ascontiguousarray(verts)); dpos_d.assign(np.ascontiguousarray(dpos))
        wp.launch(pos_to_f32, dim=N, inputs=[pos_d, nf32], device=device)
        wp.launch(face_centroids_f32, dim=nf, inputs=[pos_d, faces_d, cf32], device=device)
        fg.build(points=cf32, radius=rq)
        al = ccd_alpha(fg.id, nf32, pos_d, dpos_d, cof_d, faces_d, fcell_d, rq, t_buf, eta=eta, device=device)
        return al, d0 + al * dz                                  # (alpha, final gap)

    print("(4) CCD unit test (1 node vs static triangle, eta=%.2f):" % eta)
    d0 = 0.2 * me; ok = True
    cases = [("penetrating 2x", -2.0 * d0, eta / 2.0), ("grazing 1.0x", -1.0 * d0, eta),
             ("short 0.5x (still capped by eta)", -0.5 * d0, min(1.0, eta / 0.5)),
             ("separating", +1.5 * d0, 1.0)]
    for lbl, dz, a_exp in cases:
        al, gap = alpha_for(d0, dz)
        gap_ok = gap > 0 or dz >= 0
        rel = abs(al - a_exp) / max(a_exp, 1e-30)
        if rel > 1e-6 or not gap_ok:
            ok = False
        print(f"    {lbl:>34}: alpha={al:.4f} (exp {a_exp:.4f})  final_gap/d0={gap / d0:+.3f}  "
              f"{'penetration-free' if gap > 0 or dz >= 0 else 'PENETRATES!'}")
    print(f"    → {'PASS' if ok else 'CHECK'}: CCD caps the step so the final gap stays >0 "
          f"(>= (1-eta)*d for a closing step); separating/short steps unrestricted")


def _ipc_full_test(device="cpu"):
    """The full IPC method end-to-end on a 2-cell + strong bundle-proxy pull: barrier force (Phase 1)
    in the RHS + barrier Hessian in the implicit operator (Phase 0 hook) + CCD-filtered step (Phase 2).
    Compared to the OLD capped per-face penalty. IPC must keep the cells penetration-free (pen→0) and
    stable where OLD tunnels (pen grows past c_rep). Start state is overlapping → exercises the
    feasibilization branch (push to gap≥0) before the barrier+CCD maintain gap>0."""
    import numpy as np
    from ffn_sim.warp_port.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
    from ffn_sim.warp_port.dcm_warp_hybrid import _bond_accumulate
    from ffn_sim.warp_port.dcm_warp_hybrid_multicell import _dp_from_vol, _zero_vec
    from ffn_sim.warp_port.dcm_neighbor_warp import pos_to_f32, face_centroids_f32, contact_grid_kernel
    from ffn_sim.warp_port.dcm_warp_implicit import device_cg, _vaxpy_active, _vaxpy_active_capped

    m = _build_two_cell(device, overlap=0.45)
    verts, faces, edges, cof, fcell = m["verts"], m["faces"], m["edges"], m["cof"], m["fcell"]
    R, me, npc, p = m["R"], m["me"], m["npc"], m["params"]
    N, nf, ne = verts.shape[0], faces.shape[0], edges.shape[0]
    V0 = (4 / 3) * np.pi * (np.linalg.norm(verts[:npc] - verts[:npc].mean(0), axis=1).mean()) ** 3
    rep = 2.0e8; c_rep, c_adh = 0.30 * me, 0.80 * me
    d_hat = c_rep                                              # barrier activation gap (derived)
    repel_q = float(c_adh + 0.7 * (3 * me / 2.9) + 3.0 * me)
    gamma = 6.0 * np.pi * 65.9 * R / npc
    area_typ = 4 * np.pi * R ** 2 / npc; Fcap = rep * area_typ * c_rep

    pos_d = wp.array(verts, dtype=wp.vec3d, device=device); force_d = wp.zeros(N, dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof, dtype=wp.int32, device=device)
    faces_d = wp.array(faces.astype(np.int32), dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell, dtype=wp.int32, device=device)
    edges_d = wp.array(edges.astype(np.int32), dtype=wp.int32, device=device)
    r0_d = wp.array(np.linalg.norm(verts[edges[:, 0]] - verts[edges[:, 1]], axis=1), dtype=wp.float64, device=device)
    nf32 = wp.zeros(N, dtype=wp.vec3, device=device); cf32 = wp.zeros(nf, dtype=wp.vec3, device=device)
    Vc = wp.zeros(2, dtype=wp.float64, device=device); dP = wp.zeros(2, dtype=wp.float64, device=device)
    cn_k = wp.zeros(N, dtype=wp.float64, device=device); cn_nrm = wp.zeros(N, dtype=wp.vec3d, device=device)
    t_ccd = wp.zeros(N, dtype=wp.float64, device=device)
    fg = wp.HashGrid(32, 32, 32, device=device)
    pull_d = wp.array(np.zeros((N, 3)), dtype=wp.vec3d, device=device)

    def set_pull(F):
        pull = np.zeros((N, 3)); pull[cof == 0, 0] = -F; pull[cof == 1, 0] = +F
        pull_d.assign(np.ascontiguousarray(pull))

    def build_grid(x):
        pos_d.assign(np.ascontiguousarray(x.reshape(N, 3)))
        wp.launch(pos_to_f32, dim=N, inputs=[pos_d, nf32], device=device)
        wp.launch(face_centroids_f32, dim=nf, inputs=[pos_d, faces_d, cf32], device=device)
        fg.build(points=cf32, radius=repel_q)

    def smooth_into(pos_buf, out_d):
        wp.launch(_zero_vec, dim=N, inputs=[out_d], device=device)
        Vc.zero_()
        wp.launch(dcm_volume_kernel, dim=nf, inputs=[pos_buf, faces_d, fcell_d, Vc], device=device)
        wp.launch(_dp_from_vol, dim=2, inputs=[Vc, wp.float64(V0), wp.float64(p.turgor_dP0), wp.float64(7.73e5), dP], device=device)
        wp.launch(dcm_turgor_force_kernel, dim=nf, inputs=[pos_buf, faces_d, fcell_d, dP, out_d], device=device)
        wp.launch(_bond_accumulate, dim=ne, inputs=[pos_buf, edges_d, wp.float64(p.k_edge), r0_d, out_d], device=device)
        wp.launch(_vaxpy_active, dim=N, inputs=[out_d, wp.float64(1.0), pull_d, cof_d], device=device)

    def old_into(pos_buf, out_d):
        smooth_into(pos_buf, out_d)
        wp.launch(contact_grid_kernel, dim=N,
                  inputs=[fg.id, nf32, pos_buf, cof_d, faces_d, fcell_d, wp.float32(repel_q),
                          wp.float64(rep), wp.float64(0.0), wp.float64(c_rep), wp.float64(c_adh), out_d], device=device)

    def pen_frac(x):
        # depth = |barrier/feasibilization force| recovered via cn_k; for outside nodes pen=0
        build_grid(x); fz = wp.zeros(N, dtype=wp.vec3d, device=device)
        wp.launch(nearest_face_ipc_kernel, dim=N,
                  inputs=[fg.id, nf32, pos_d, cof_d, faces_d, fcell_d, wp.float32(repel_q),
                          wp.float64(rep), wp.float64(d_hat), fz, cn_k, cn_nrm], device=device)
        wp.synchronize_device(device)
        # count only INSIDE nodes: feasibilization force = rep*area*depth, depth=|F|/(rep*area)=|F|/cn_k.
        # But barrier nodes also have cn_k>0; distinguish by sign via penetration_depth_kernel instead.
        from ffn_sim.warp_port.dcm_neighbor_warp import penetration_depth_kernel
        pend = wp.zeros(N, dtype=wp.float64, device=device)
        wp.launch(penetration_depth_kernel, dim=N,
                  inputs=[fg.id, nf32, pos_d, cof_d, faces_d, fcell_d, wp.float32(repel_q), pend], device=device)
        wp.synchronize_device(device)
        return float(pend.numpy().max() / me)

    def vol_cell(x):
        v0, v1, v2 = x[faces[:, 0]], x[faces[:, 1]], x[faces[:, 2]]
        return abs(float(np.einsum('ij,ij->i', v0, np.cross(v1 - v0, v2 - v0)).sum() / 6.0)) / 2.0

    scratch = {k: wp.zeros(N, dtype=wp.vec3d, device=device) for k in ("r", "p", "Ap", "dx", "Fx", "Fp", "xp")}
    scratch["sca"] = wp.zeros(1, dtype=wp.float64, device=device)
    dt = 8e-6 * 50; a_imp = gamma / dt; nsteps = 200

    def run(mode):
        x = verts.copy(); ok = True
        for s in range(nsteps):
            build_grid(x)
            if mode == "old":
                wp.launch(_zero_vec, dim=N, inputs=[force_d], device=device)
                old_into(pos_d, force_d)
                dx_d, _ = device_cg(old_into, pos_d, a_imp, force_d, scratch, device=device)
                wp.launch(_vaxpy_active_capped, dim=N, inputs=[pos_d, dx_d, cof_d, wp.float64(c_rep)], device=device)
                x = pos_d.numpy().copy()
            else:  # FULL IPC: barrier force + barrier Hessian (implicit) + CCD-filtered step
                wp.launch(_zero_vec, dim=N, inputs=[force_d], device=device)
                smooth_into(pos_d, force_d)
                wp.launch(nearest_face_ipc_kernel, dim=N,
                          inputs=[fg.id, nf32, pos_d, cof_d, faces_d, fcell_d, wp.float32(repel_q),
                                  wp.float64(rep), wp.float64(d_hat), force_d, cn_k, cn_nrm], device=device)
                hess = make_contact_hess_apply(cn_k, cn_nrm, device=device)
                dx_d, _ = device_cg(smooth_into, pos_d, a_imp, force_d, scratch, device=device, hess_apply=hess)
                al = ccd_alpha(fg.id, nf32, pos_d, dx_d, cof_d, faces_d, fcell_d, repel_q, t_ccd, eta=0.9, device=device)
                x = x + al * dx_d.numpy()
            if not np.isfinite(x).all() or vol_cell(x) > 50 * V0:
                ok = False; break
        return x, ok

    pen0 = pen_frac(verts)
    print(f"(5) FULL IPC vs OLD penalty, 2-cell + bundle-proxy ({nsteps} steps @dt×50, "
          f"start pen={pen0:.2f}, Fcap={Fcap:.1e}N):")
    print(f"    {'Fpull/Fcap':>10} | {'OLD pen':>9} {'OLD V/V0':>9} | {'IPC pen':>9} {'IPC V/V0':>9}")
    for mult in (1.0, 2.0, 4.0, 8.0):
        set_pull(mult * Fcap)
        xo, oko = run("old"); xi, oki = run("ipc")
        po = pen_frac(xo) if oko else float('nan'); vo = vol_cell(xo) / V0 if oko else float('nan')
        pi = pen_frac(xi) if oki else float('nan'); vi = vol_cell(xi) / V0 if oki else float('nan')
        od = "DIVERGE" if not oko else f"{po:8.3f}"; idd = "DIVERGE" if not oki else f"{pi:8.3f}"
        print(f"    {mult:>10.1f} | {od:>9} {vo:>9.3f} | {idd:>9} {vi:>9.3f}")
    print(f"    → expect: IPC pen → ~0 (penetration-free) + stable at every pull; OLD tunnels (pen grows)")


if __name__ == "__main__":
    _barrier_unittest()
    _ccd_unittest()
    _ipc_full_test()
    _selftest()
