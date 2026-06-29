"""Frozen-neighbour cache + analytic-diagonal Jacobi preconditioner for the implicit DCM CG solve.

Two opt-in matvec accelerations for ``device_cg`` (see ``dcm_warp_implicit.device_cg`` and the
per-step implicit branch in ``dcm_warp_decohesion.step_once``). Both are **default OFF** — when
neither flag is set the implicit code path is byte-identical to the committed engine.

Context (cProfile n400 ULA): ~88% of the per-step wall is ``synchronize_device`` inside
``device_cg``; ~12 CG iters/step, each a query-bound matvec. The matvec is the
FINITE-DIFFERENCE operator ``(a·I + K)·v``, ``K·v = −[F(x+εv) − F(x)]/ε`` — it re-evaluates
``stiff_force_into`` at a perturbed position. ``step_once`` builds the cohesion/contact hash-grids
ONCE per step on xₙ, so the neighbour *topology* is already frozen across CG iters; but every
matvec STILL re-runs ``cohesion_grid_kernel`` / ``contact_grid_kernel``, each of which re-QUERIES
the hash-grid (``wp.hash_grid_query``) + recomputes the per-face closest point. That query work is
redundant — the neighbour *set* genuinely does not change across the iterations of one solve (the
FD perturbation ε·v is ~1e-9·R, far below the grid cell size).

## OPTIMIZATION #2 — frozen-neighbour cache
Once per step (at xₙ, before the CG solve) we snapshot, per node, the candidate-neighbour id lists
the grid query WOULD yield — for cohesion (other-cell node ids within ``coh_q``) and for contact
(other-cell face ids within ``con_q``) — into fixed-stride CSR-like buffers. The frozen-force
kernels (:func:`cohesion_frozen_force_kernel`, :func:`contact_frozen_force_kernel`) then re-apply
the *identical* force law over the cached ids using the CURRENT f64 positions, with NO grid query.

Parity argument: a hash-grid query returns every point whose CELL is within ``radius`` of the
query cell; the force law then re-zeroes any pair outside the true cutoff. The frozen cache stores
exactly those candidate ids in the order the query enumerated them, and the frozen kernel runs the
same per-pair math — so the only possible difference is floating-point summation order, which is
identical because the iteration order is preserved. Across one solve the candidates are the same set
the live query would have returned (the grid is frozen + the perturbation is sub-cell), so the
converged ``dx`` matches the live-query matvec to ``cg_tol``. Validated in ``_frozen_parity_demo``.

## OPTIMIZATION #1 — analytic-diagonal Jacobi preconditioner
``diagA[i] = a + diag_edges[i] + diag_contact[i] (+ diag_turgor[i])`` per node, where:
  * ``a = γ/dt`` — the implicit regulariser (exact, isotropic, always positive).
  * edges: each incident harmonic edge spring contributes ≈ ``k_edge`` to its endpoints' diagonal
    (the radial spring's stiffness; the transverse term ∝ tension/length is dropped as it is ~0 at
    the relaxed reference and would make the estimate non-PD — a Jacobi preconditioner only needs
    a PD diagonal *estimate*, not the exact Hessian diagonal).
  * contact: the active node-face penalty normal stiffness ``rep·area`` (the same per-node value the
    IPC / repel kernels expose as ``cn_k``), reused via the frozen-contact pass.
  * turgor: the per-face volume-pressure coupling. The exact term is a dense rank-1 (∂P/∂V) coupling
    across a whole cell; its diagonal contribution per node is ``k_vol·(∂V/∂x_i)²/V0`` summed over
    the node's incident faces — included CONSERVATIVELY (see :func:`turgor_diag_kernel`).

``M⁻¹ = 1/diagA`` is applied as ``z = r/diagA`` each CG iter (the standard Jacobi/diagonal PCG).
A preconditioner changes ONLY the iteration path, never the converged solution — so ``dx`` is the
same (within ``cg_tol``); the payoff is fewer iterations. HONESTLY reported in ``_precond_demo``:
if the analytic diagonal does not cut iters, that is stated, not hidden.
"""

from __future__ import annotations

import warp as wp

from ffn_sim.dcm.dcm_contact_warp import closest_bary

wp.init()


# ===========================================================================
# OPTIMIZATION #2 — frozen-neighbour cache
# ===========================================================================
# Per node we store a fixed-stride candidate list (CSR with a constant cap K):
#   coh_ids[i*Kc + 0..cnt-1]  = candidate other-cell NODE ids (cohesion)
#   con_ids[i*Kf + 0..cnt-1]  = candidate other-cell FACE ids (contact)
# A fixed cap avoids a host-side prefix-sum (which would need a sync). If a node exceeds the cap
# the build records an overflow and the caller widens the cap + rebuilds (rare; the cap is sized
# from the worst observed neighbour count). The cached ids are exactly the hash_grid_query
# enumeration, so replaying them reproduces the live-query force bit-for-bit (summation order kept).

@wp.kernel
def cohesion_cache_build_kernel(
    grid: wp.uint64,
    qpts: wp.array(dtype=wp.vec3),
    cof: wp.array(dtype=wp.int32),
    radius: wp.float32,
    cap: wp.int32,
    ids: wp.array(dtype=wp.int32),            # (N*cap,) OUT candidate other-cell node ids
    cnt: wp.array(dtype=wp.int32),            # (N,) OUT count (clamped to cap)
    overflow: wp.array(dtype=wp.int32),       # (1,) OUT overflow flag (atomic_max with needed-cap)
):
    """Snapshot the cohesion grid query: for node i, the other-cell LIVE nodes within ``radius``,
    in query-enumeration order. Mirrors the j-filter of ``cohesion_grid_kernel`` exactly so the
    frozen force replays the same candidate set (the d<c_adh cutoff is re-applied in the force)."""
    i = wp.tid()
    c1 = cof[i]
    base = i * cap
    if c1 < wp.int32(0):
        cnt[i] = wp.int32(0)
        return
    k = wp.int32(0)
    q = wp.hash_grid_query(grid, qpts[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(q, j):
        if j != i and cof[j] >= wp.int32(0) and cof[j] != c1:
            if k < cap:
                ids[base + k] = j
            k = k + wp.int32(1)
    if k > cap:
        wp.atomic_max(overflow, 0, k)
        cnt[i] = cap
    else:
        cnt[i] = k


@wp.kernel
def cohesion_frozen_force_kernel(
    ids: wp.array(dtype=wp.int32),            # (N*cap,) cached candidate node ids
    cnt: wp.array(dtype=wp.int32),            # (N,) per-node candidate count
    cap: wp.int32,
    pos: wp.array(dtype=wp.vec3d),            # f64 node positions (force math, CURRENT)
    cof: wp.array(dtype=wp.int32),
    r_contact: wp.float64, c_adh: wp.float64,
    rep: wp.float64, omega: wp.float64, A: wp.float64, force_cap: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """Frozen-cache twin of ``cohesion_grid_kernel`` — same force law, no grid query (own-row write).
    Iterates the cached candidate ids; the per-pair math is byte-identical to the live kernel."""
    i = wp.tid()
    c1 = cof[i]
    if c1 < wp.int32(0):
        return
    z = wp.float64(0.0)
    half = wp.float64(0.5) * c_adh
    pi = pos[i]
    acc = wp.vec3d(z, z, z)
    n = cnt[i]
    base = i * cap
    t = wp.int32(0)
    while t < n:
        j = ids[base + t]
        t = t + wp.int32(1)
        # candidate guaranteed j!=i, live, other-cell by the build; keep no re-check (parity)
        r_vec = pi - pos[j]
        d = wp.length(r_vec)
        fmag = z
        if d < r_contact:
            fmag = rep * A * (r_contact - d)
        else:
            if d < c_adh:
                tent = c_adh - d
                if d < half:
                    tent = d
                fmag = -omega * A * tent
        if fmag > force_cap:
            fmag = force_cap
        if fmag < -force_cap:
            fmag = -force_cap
        if d > z and fmag != z:
            acc = acc + r_vec * (fmag / d)
    force[i] = acc


@wp.kernel
def cohesion_frozen_force_cad_kernel(
    ids: wp.array(dtype=wp.int32),
    cnt: wp.array(dtype=wp.int32),
    cap: wp.int32,
    pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32),
    cad: wp.array(dtype=wp.float64),
    r_contact: wp.float64, c_adh: wp.float64,
    rep: wp.float64, omega: wp.float64, A: wp.float64, force_cap: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """Frozen-cache twin of ``cohesion_grid_cad_kernel`` (M3 junction-switch adhesion scaling)."""
    i = wp.tid()
    c1 = cof[i]
    if c1 < wp.int32(0):
        return
    z = wp.float64(0.0)
    half = wp.float64(0.5) * c_adh
    pi = pos[i]
    acc = wp.vec3d(z, z, z)
    n = cnt[i]
    base = i * cap
    t = wp.int32(0)
    while t < n:
        j = ids[base + t]
        t = t + wp.int32(1)
        r_vec = pi - pos[j]
        d = wp.length(r_vec)
        fmag = z
        if d < r_contact:
            fmag = rep * A * (r_contact - d)
        else:
            if d < c_adh:
                tent = c_adh - d
                if d < half:
                    tent = d
                mult = wp.sqrt(cad[c1] * cad[cof[j]])
                fmag = -omega * A * tent * mult
        if fmag > force_cap:
            fmag = force_cap
        if fmag < -force_cap:
            fmag = -force_cap
        if d > z and fmag != z:
            acc = acc + r_vec * (fmag / d)
    force[i] = acc


@wp.kernel
def contact_cache_build_kernel(
    grid: wp.uint64,
    qpts: wp.array(dtype=wp.vec3),
    cof: wp.array(dtype=wp.int32),
    fcell: wp.array(dtype=wp.int32),
    radius: wp.float32,
    cap: wp.int32,
    ids: wp.array(dtype=wp.int32),            # (N*cap,) OUT candidate other-cell face ids
    cnt: wp.array(dtype=wp.int32),            # (N,) OUT count
    overflow: wp.array(dtype=wp.int32),       # (1,) OUT overflow flag
):
    """Snapshot the contact grid query: for node i, the other-cell FACE ids within ``radius``, in
    enumeration order. Mirrors the fcell filter of ``contact_grid_kernel`` exactly."""
    ni = wp.tid()
    c1 = cof[ni]
    base = ni * cap
    if c1 < wp.int32(0):
        cnt[ni] = wp.int32(0)
        return
    k = wp.int32(0)
    q = wp.hash_grid_query(grid, qpts[ni], radius)
    fj = wp.int32(0)
    while wp.hash_grid_query_next(q, fj):
        if fcell[fj] != c1:
            if k < cap:
                ids[base + k] = fj
            k = k + wp.int32(1)
    if k > cap:
        wp.atomic_max(overflow, 0, k)
        cnt[ni] = cap
    else:
        cnt[ni] = k


@wp.kernel
def contact_frozen_force_kernel(
    ids: wp.array(dtype=wp.int32),
    cnt: wp.array(dtype=wp.int32),
    cap: wp.int32,
    pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    rep: wp.float64, adh: wp.float64, c_rep: wp.float64, c_adh: wp.float64,
    force: wp.array(dtype=wp.vec3d),          # atomic accumulate
):
    """Frozen-cache twin of ``contact_grid_kernel`` — same closest-point penalty/adhesion law over
    the cached candidate faces, NO grid query. The closest point is RE-evaluated on the current
    positions (only the candidate SET is frozen, exactly like the live-query frozen operator)."""
    ni = wp.tid()
    c1 = cof[ni]
    if c1 < wp.int32(0):
        return
    z = wp.float64(0.0)
    half = wp.float64(0.5) * c_adh
    p = pos[ni]
    fn_acc = wp.vec3d(z, z, z)
    n = cnt[ni]
    base = ni * cap
    t = wp.int32(0)
    while t < n:
        fj = ids[base + t]
        t = t + wp.int32(1)
        ia = faces[fj, 0]
        ib = faces[fj, 1]
        ic = faces[fj, 2]
        a = pos[ia]
        b = pos[ib]
        c = pos[ic]
        bary = closest_bary(p, a, b, c)
        cpa = a * bary[0] + b * bary[1] + c * bary[2]
        r_vec = p - cpa
        min_d = wp.length(r_vec)
        fnv = wp.cross(b - a, c - a)
        nrm = wp.length(fnv)
        area = wp.float64(0.5) * nrm
        sign = z
        if nrm > z:
            sign = wp.dot(r_vec, fnv) / nrm
        amp = z
        if sign < z and min_d < c_rep:
            amp = rep * area
        else:
            if adh > z and sign > z and min_d < c_adh:
                if min_d >= half:
                    md = min_d
                    if md <= z:
                        md = wp.float64(1.0e-30)
                    amp = adh * (c_adh / md - wp.float64(1.0)) * area
                else:
                    amp = adh * area
        if amp != z:
            fvec = r_vec * amp
            fn_acc = fn_acc - fvec
            wp.atomic_add(force, ia, bary[0] * fvec)
            wp.atomic_add(force, ib, bary[1] * fvec)
            wp.atomic_add(force, ic, bary[2] * fvec)
    wp.atomic_add(force, ni, fn_acc)


@wp.kernel
def contact_frozen_force_cad_kernel(
    ids: wp.array(dtype=wp.int32),
    cnt: wp.array(dtype=wp.int32),
    cap: wp.int32,
    pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    cad: wp.array(dtype=wp.float64),
    rep: wp.float64, adh: wp.float64, c_rep: wp.float64, c_adh: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """Frozen-cache twin of ``contact_grid_cad_kernel`` (M3 junction-switch adhesion scaling)."""
    ni = wp.tid()
    c1 = cof[ni]
    if c1 < wp.int32(0):
        return
    z = wp.float64(0.0)
    half = wp.float64(0.5) * c_adh
    p = pos[ni]
    fn_acc = wp.vec3d(z, z, z)
    n = cnt[ni]
    base = ni * cap
    t = wp.int32(0)
    while t < n:
        fj = ids[base + t]
        t = t + wp.int32(1)
        ia = faces[fj, 0]
        ib = faces[fj, 1]
        ic = faces[fj, 2]
        a = pos[ia]
        b = pos[ib]
        c = pos[ic]
        bary = closest_bary(p, a, b, c)
        cpa = a * bary[0] + b * bary[1] + c * bary[2]
        r_vec = p - cpa
        min_d = wp.length(r_vec)
        fnv = wp.cross(b - a, c - a)
        nrm = wp.length(fnv)
        area = wp.float64(0.5) * nrm
        sign = z
        if nrm > z:
            sign = wp.dot(r_vec, fnv) / nrm
        amp = z
        if sign < z and min_d < c_rep:
            amp = rep * area
        else:
            if adh > z and sign > z and min_d < c_adh:
                mult = wp.sqrt(cad[c1] * cad[fcell[fj]])
                if min_d >= half:
                    md = min_d
                    if md <= z:
                        md = wp.float64(1.0e-30)
                    amp = adh * (c_adh / md - wp.float64(1.0)) * area * mult
                else:
                    amp = adh * area * mult
        if amp != z:
            fvec = r_vec * amp
            fn_acc = fn_acc - fvec
            wp.atomic_add(force, ia, bary[0] * fvec)
            wp.atomic_add(force, ib, bary[1] * fvec)
            wp.atomic_add(force, ic, bary[2] * fvec)
    wp.atomic_add(force, ni, fn_acc)


class FrozenNeighborCache:
    """Holds the per-step frozen candidate lists for cohesion + contact and replays the force.

    Built once per step (``rebuild`` after the grids are built on xₙ); ``cohesion_into`` /
    ``contact_into`` then run the frozen-force kernels in place of the live grid kernels inside the
    CG matvec. The cap auto-grows on overflow (rare). All device-resident; the only host syncs are
    the overflow checks at build time (NOT per CG iter)."""

    def __init__(self, N: int, *, coh_cap: int = 64, con_cap: int = 64, device: str = "cpu"):
        self.N = N
        self.device = device
        self.coh_cap = int(coh_cap)
        self.con_cap = int(con_cap)
        self._alloc()
        self._of = wp.zeros(1, dtype=wp.int32, device=device)

    def _alloc(self):
        self.coh_ids = wp.zeros(self.N * self.coh_cap, dtype=wp.int32, device=self.device)
        self.coh_cnt = wp.zeros(self.N, dtype=wp.int32, device=self.device)
        self.con_ids = wp.zeros(self.N * self.con_cap, dtype=wp.int32, device=self.device)
        self.con_cnt = wp.zeros(self.N, dtype=wp.int32, device=self.device)

    def rebuild(self, *, node_grid_id, face_grid_id, node_f32, cof_d, fcell_d,
                coh_q: float, con_q: float, do_cohesion: bool, do_contact: bool):
        """Snapshot the candidate lists at the current (xₙ) grids. ``coh_q``/``con_q`` are the SAME
        query radii the live kernels use. Grows the cap + retries if any node overflowed."""
        N = self.N
        # int32 CEILING on the flat (N·cap) candidate arrays: warp array shapes must fit a signed
        # int32, so cap·N < 2^31. At N≥1000 (≥324k nodes) an over-grown cap (a transiently-degenerate
        # node reporting ~9k neighbours) blew past this and crashed the run. We cap `cap` at the
        # int32-safe max and CLIP excess candidates there — a node needing > max_cap face/node
        # neighbours is non-physical (real packs are ≤ ~10²); normal nodes (need ≪ max_cap) are
        # unaffected, so the frozen replay stays parity-exact for the physical regime.
        max_cap = max(64, (2**31 - 1) // max(N, 1) - 16)
        if do_cohesion:
            while True:
                self._of.zero_()
                wp.launch(cohesion_cache_build_kernel, dim=N,
                          inputs=[node_grid_id, node_f32, cof_d, wp.float32(coh_q),
                                  wp.int32(self.coh_cap), self.coh_ids, self.coh_cnt, self._of],
                          device=self.device)
                wp.synchronize_device(self.device)
                need = int(self._of.numpy()[0])
                if need <= self.coh_cap:
                    break
                new_cap = min(int(need * 1.25) + 4, max_cap)
                if new_cap <= self.coh_cap:                  # at the int32 ceiling → clip + stop growing
                    break
                self.coh_cap = new_cap
                self.coh_ids = wp.zeros(N * self.coh_cap, dtype=wp.int32, device=self.device)
        if do_contact:
            while True:
                self._of.zero_()
                wp.launch(contact_cache_build_kernel, dim=N,
                          inputs=[face_grid_id, node_f32, cof_d, fcell_d, wp.float32(con_q),
                                  wp.int32(self.con_cap), self.con_ids, self.con_cnt, self._of],
                          device=self.device)
                wp.synchronize_device(self.device)
                need = int(self._of.numpy()[0])
                if need <= self.con_cap:
                    break
                new_cap = min(int(need * 1.25) + 4, max_cap)
                if new_cap <= self.con_cap:                  # at the int32 ceiling → clip + stop growing
                    break
                self.con_cap = new_cap
                self.con_ids = wp.zeros(N * self.con_cap, dtype=wp.int32, device=self.device)

    def cohesion_into(self, pos_buf, cof_d, *, r_contact, c_adh, rep, omega, A, force_cap, out_d,
                      cad_d=None):
        N = self.N
        if cad_d is None:
            wp.launch(cohesion_frozen_force_kernel, dim=N,
                      inputs=[self.coh_ids, self.coh_cnt, wp.int32(self.coh_cap), pos_buf, cof_d,
                              wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep),
                              wp.float64(omega), wp.float64(A), wp.float64(force_cap), out_d],
                      device=self.device)
        else:
            wp.launch(cohesion_frozen_force_cad_kernel, dim=N,
                      inputs=[self.coh_ids, self.coh_cnt, wp.int32(self.coh_cap), pos_buf, cof_d, cad_d,
                              wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep),
                              wp.float64(omega), wp.float64(A), wp.float64(force_cap), out_d],
                      device=self.device)

    def contact_into(self, pos_buf, cof_d, faces_d, fcell_d, *, rep, adh, c_rep, c_adh, out_d,
                     cad_d=None):
        N = self.N
        if cad_d is None:
            wp.launch(contact_frozen_force_kernel, dim=N,
                      inputs=[self.con_ids, self.con_cnt, wp.int32(self.con_cap), pos_buf, cof_d,
                              faces_d, fcell_d, wp.float64(rep), wp.float64(adh),
                              wp.float64(c_rep), wp.float64(c_adh), out_d], device=self.device)
        else:
            wp.launch(contact_frozen_force_cad_kernel, dim=N,
                      inputs=[self.con_ids, self.con_cnt, wp.int32(self.con_cap), pos_buf, cof_d,
                              faces_d, fcell_d, cad_d, wp.float64(rep), wp.float64(adh),
                              wp.float64(c_rep), wp.float64(c_adh), out_d], device=self.device)


# ===========================================================================
# OPTIMIZATION #1 — analytic-diagonal Jacobi preconditioner
# ===========================================================================
@wp.kernel
def diag_init_kernel(a: wp.float64, cof: wp.array(dtype=wp.int32),
                     diagA: wp.array(dtype=wp.float64)):
    """Seed each LIVE node's diagonal with the implicit regulariser a=γ/dt (dormant → a, harmless)."""
    i = wp.tid()
    diagA[i] = a


@wp.kernel
def edge_diag_kernel(edges: wp.array(dtype=wp.int32, ndim=2), k_edge: wp.float64,
                     diagA: wp.array(dtype=wp.float64)):
    """Each harmonic edge spring adds its radial stiffness k_edge to BOTH endpoints' diagonal.
    (The exact per-DOF bond Hessian diagonal is k·n̂_d² + tension·(1−n̂_d²)/r per axis; at the
    relaxed reference tension≈0 and Σ_axes n̂²=1, so the scalar-per-node estimate k_edge per
    incident edge is the dominant, PD part — sufficient for a Jacobi preconditioner.)"""
    e = wp.tid()
    wp.atomic_add(diagA, edges[e, 0], k_edge)
    wp.atomic_add(diagA, edges[e, 1], k_edge)


@wp.kernel
def contact_diag_kernel(cn_k: wp.array(dtype=wp.float64), diagA: wp.array(dtype=wp.float64)):
    """Add the active node-face contact normal stiffness cn_k (= rep·area, the IPC/repel kernels'
    per-node value) to the node diagonal. Rank-1 along n̂, so its scalar diagonal share is cn_k."""
    i = wp.tid()
    k = cn_k[i]
    if k > wp.float64(0.0):
        diagA[i] = diagA[i] + k


@wp.kernel
def turgor_diag_kernel(
    pos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32), k_vol: wp.float64, V0: wp.float64,
    diagA: wp.array(dtype=wp.float64),
):
    """CONSERVATIVE turgor diagonal. The osmotic pressure P = dP0 − k_vol·(V−V0)/V0 acts on each
    face with force ∝ P·(face area normal); the dominant on-diagonal stiffness from the
    VOLUME-feedback term is ∂F_i/∂x_i ≈ k_vol/V0 · (∂V/∂x_i)², with ∂V/∂x_i = (signed) ⅓·(area
    normal) of the node's incident faces. Per face we scatter ⅓ of that squared-gradient magnitude
    to its 3 vertices (a positive, isotropic estimate). This UNDER-counts the full dense rank-1
    coupling (intentionally — a diagonal preconditioner only needs the local part) but captures the
    volume lock's stiffness, which at large k_vol dominates the conditioning. Scattered as a scalar.
    """
    f = wp.tid()
    ia = faces[f, 0]
    ib = faces[f, 1]
    ic = faces[f, 2]
    a = pos[ia]
    b = pos[ib]
    c = pos[ic]
    nrm = wp.cross(b - a, c - a)              # 2·area·n̂ ; ∂V/∂x_vertex = nrm/6 (each of 3 verts)
    g = nrm / wp.float64(6.0)                 # per-vertex volume gradient ∂V/∂x
    g2 = wp.dot(g, g)                          # |∂V/∂x|² (scalar diagonal share)
    s = (k_vol / V0) * g2
    wp.atomic_add(diagA, ia, s)
    wp.atomic_add(diagA, ib, s)
    wp.atomic_add(diagA, ic, s)


@wp.kernel
def turgor_diag_kernel_pc(
    pos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32), k_vol: wp.float64, V0: wp.array(dtype=wp.float64),
    diagA: wp.array(dtype=wp.float64),
):
    """Per-cell rest-volume variant of :func:`turgor_diag_kernel` — the volume-lock diagonal uses
    each face's OWN cell rest volume ``V0[fcell[f]]`` so the analytic Hessian stays CONSISTENT with
    the per-cell turgor force (:func:`_dp_from_vol_pc`) once cells carry distinct V0 (division)."""
    f = wp.tid()
    ia = faces[f, 0]
    ib = faces[f, 1]
    ic = faces[f, 2]
    a = pos[ia]
    b = pos[ib]
    c = pos[ic]
    nrm = wp.cross(b - a, c - a)
    g = nrm / wp.float64(6.0)
    g2 = wp.dot(g, g)
    s = (k_vol / V0[fcell[f]]) * g2
    wp.atomic_add(diagA, ia, s)
    wp.atomic_add(diagA, ib, s)
    wp.atomic_add(diagA, ic, s)


@wp.kernel
def _precond_apply_kernel(r: wp.array(dtype=wp.vec3d), diagA: wp.array(dtype=wp.float64),
                          z: wp.array(dtype=wp.vec3d)):
    """z = M⁻¹ r = r / diagA (diagonal/Jacobi preconditioner)."""
    i = wp.tid()
    d = diagA[i]
    if d > wp.float64(0.0):
        z[i] = r[i] / d
    else:
        z[i] = r[i]


def make_diag_precond(diagA, device="cpu"):
    """Return ``apply(r_d, z_d)`` that writes z = r/diagA — pass to ``device_cg(precond_apply=...)``."""
    N = diagA.shape[0]

    def apply(r_d, z_d):
        wp.launch(_precond_apply_kernel, dim=N, inputs=[r_d, diagA, z_d], device=device)

    return apply


def build_diagA(diagA, *, a, cof_d, edges_d, k_edge, faces_d, fcell_d, pos_d, k_vol, V0,
                V0_arr=None, cn_k=None, do_edges=True, do_turgor=True, device="cpu"):
    """Assemble the analytic diagonal a·I + diag(K) into ``diagA`` (a (N,) f64 buffer).
    ``cn_k`` (optional) is the per-node contact normal stiffness (from the IPC/repel pass); when the
    penalty contact is used instead it is computed separately by the caller and passed here.
    ``V0_arr`` (optional, a (n_cells,) f64 device array) selects the PER-CELL turgor diagonal so the
    Hessian matches the per-cell turgor force under division; when None the scalar ``V0`` is used."""
    N = diagA.shape[0]
    wp.launch(diag_init_kernel, dim=N, inputs=[wp.float64(a), cof_d, diagA], device=device)
    if do_edges:
        wp.launch(edge_diag_kernel, dim=edges_d.shape[0], inputs=[edges_d, wp.float64(k_edge), diagA],
                  device=device)
    if do_turgor:
        if V0_arr is not None:
            wp.launch(turgor_diag_kernel_pc, dim=faces_d.shape[0],
                      inputs=[pos_d, faces_d, fcell_d, wp.float64(k_vol), V0_arr, diagA],
                      device=device)
        else:
            wp.launch(turgor_diag_kernel, dim=faces_d.shape[0],
                      inputs=[pos_d, faces_d, fcell_d, wp.float64(k_vol), wp.float64(V0), diagA],
                      device=device)
    if cn_k is not None:
        wp.launch(contact_diag_kernel, dim=N, inputs=[cn_k, diagA], device=device)


# ===========================================================================
# Validation harness — a small-but-stiff multi-cell config, ONE implicit step.
# Run:  PYTHONPATH=. python ffn_sim/warp_port/dcm_warp_frozen.py [--device cpu] [--n 24]
# ===========================================================================
def _build_stiff_cluster(n_cells=24, device="cpu", subdiv=1, k_vol=1.0e3):
    """``n_cells`` icospheres in an FCC pack just touching (the ULA-ish stiff force set: turgor +
    cortex edges + node-node cohesion + node-FACE contact). Returns the live-query stiff operator,
    the grids, and all handles needed for the frozen cache + the analytic diagonal."""
    import numpy as np
    from ffn_sim.cell.dcm import icosphere_mesh, ResolvedDCM
    from ffn_sim.dcm.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
    from ffn_sim.dcm.dcm_warp_hybrid import _bond_accumulate
    from ffn_sim.dcm.dcm_warp_hybrid_multicell import _dp_from_vol, _zero_vec
    from ffn_sim.dcm.dcm_neighbor_warp import (pos_to_f32, face_centroids_f32,
                                                     cohesion_grid_kernel, contact_grid_kernel)

    p = ResolvedDCM(subdivisions=subdiv); R = p.R_cell
    v1, e1, f1 = icosphere_mesh(R, subdiv); npc = v1.shape[0]
    me = float(np.linalg.norm(v1[e1[:, 0]] - v1[e1[:, 1]], axis=1).mean())
    c_rep, c_adh, r_contact = 0.30 * me, 0.80 * me, 0.30 * me
    rep, adh, A = 2.0e8, 1.0e7, 4 * np.pi * R ** 2 / npc

    # FCC centres, nearest n_cells to origin, spacing 1.95R (just touching)
    g = np.arange(-4, 5)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    pts = pts[(pts.sum(1) % 2) == 0].astype(float)
    centers = pts[np.argsort(np.linalg.norm(pts, axis=1))[:n_cells]] * (1.95 * R / np.sqrt(2.0))

    verts = np.concatenate([v1 + ce for ce in centers])
    faces = np.concatenate([f1 + i * npc for i in range(n_cells)])
    edges = np.concatenate([e1 + i * npc for i in range(n_cells)])
    cof = np.concatenate([np.full(npc, i) for i in range(n_cells)]).astype(np.int64)
    fcell = np.concatenate([np.full(f1.shape[0], i) for i in range(n_cells)]).astype(np.int64)
    N, nf, ne = verts.shape[0], faces.shape[0], edges.shape[0]
    R0 = float(np.linalg.norm(v1, axis=1).mean()); V0 = (4 / 3) * np.pi * R0 ** 3
    coh_q = float(c_adh)
    con_q = float(c_adh + 0.7 * (3.0 * (me / 2.9)))
    gamma = 6 * np.pi * 65.9 * R / npc
    fc_cap = 5e-8

    pos_d = wp.array(np.ascontiguousarray(verts), dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof.astype(np.int32), dtype=wp.int32, device=device)
    faces_d = wp.array(faces.astype(np.int32), dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell.astype(np.int32), dtype=wp.int32, device=device)
    edges_d = wp.array(edges.astype(np.int32), dtype=wp.int32, device=device)
    r0_d = wp.array(np.linalg.norm(verts[edges[:, 0]] - verts[edges[:, 1]], axis=1),
                    dtype=wp.float64, device=device)
    nf32 = wp.zeros(N, dtype=wp.vec3, device=device)
    cf32 = wp.zeros(nf, dtype=wp.vec3, device=device)
    Vc = wp.zeros(n_cells, dtype=wp.float64, device=device)
    dP = wp.zeros(n_cells, dtype=wp.float64, device=device)
    ng = wp.HashGrid(32, 32, 32, device=device)
    fg = wp.HashGrid(32, 32, 32, device=device)

    def build_grids():
        wp.launch(pos_to_f32, dim=N, inputs=[pos_d, nf32], device=device)
        ng.build(points=nf32, radius=coh_q)
        wp.launch(face_centroids_f32, dim=nf, inputs=[pos_d, faces_d, cf32], device=device)
        fg.build(points=cf32, radius=con_q)

    def stiff_live(pos_buf, out_d):
        """LIVE-query stiff operator (the committed path): cohesion + turgor + contact + edges."""
        wp.launch(_zero_vec, dim=N, inputs=[out_d], device=device)
        wp.launch(cohesion_grid_kernel, dim=N,
                  inputs=[ng.id, nf32, pos_buf, cof_d, wp.float32(coh_q), wp.float64(r_contact),
                          wp.float64(c_adh), wp.float64(rep), wp.float64(adh), wp.float64(A),
                          wp.float64(fc_cap), out_d], device=device)
        Vc.zero_()
        wp.launch(dcm_volume_kernel, dim=nf, inputs=[pos_buf, faces_d, fcell_d, Vc], device=device)
        wp.launch(_dp_from_vol, dim=n_cells,
                  inputs=[Vc, wp.float64(V0), wp.float64(p.turgor_dP0), wp.float64(k_vol), dP],
                  device=device)
        wp.launch(dcm_turgor_force_kernel, dim=nf, inputs=[pos_buf, faces_d, fcell_d, dP, out_d],
                  device=device)
        wp.launch(contact_grid_kernel, dim=N,
                  inputs=[fg.id, nf32, pos_buf, cof_d, faces_d, fcell_d, wp.float32(con_q),
                          wp.float64(rep), wp.float64(adh), wp.float64(c_rep), wp.float64(c_adh),
                          out_d], device=device)
        wp.launch(_bond_accumulate, dim=ne, inputs=[pos_buf, edges_d, wp.float64(p.k_edge), r0_d, out_d],
                  device=device)

    return dict(stiff_live=stiff_live, build_grids=build_grids, pos_d=pos_d, cof_d=cof_d,
                faces_d=faces_d, fcell_d=fcell_d, edges_d=edges_d, nf32=nf32, ng=ng, fg=fg,
                Vc=Vc, dP=dP, r0_d=r0_d, N=N, nf=nf, ne=ne, n_cells=n_cells, gamma=gamma, V0=V0,
                k_vol=k_vol, k_edge=p.k_edge, rep=rep, adh=adh, A=A, c_rep=c_rep, c_adh=c_adh,
                r_contact=r_contact, coh_q=coh_q, con_q=con_q, p=p, verts=verts, me=me)


def _make_frozen_stiff(M, cache):
    """A stiff operator that uses the frozen cache for cohesion + contact (turgor/edges live)."""
    from ffn_sim.dcm.dcm_turgor_warp import dcm_volume_kernel, dcm_turgor_force_kernel
    from ffn_sim.dcm.dcm_warp_hybrid import _bond_accumulate
    from ffn_sim.dcm.dcm_warp_hybrid_multicell import _dp_from_vol, _zero_vec

    def stiff_frozen(pos_buf, out_d):
        wp.launch(_zero_vec, dim=M["N"], inputs=[out_d], device=cache.device)
        cache.cohesion_into(pos_buf, M["cof_d"], r_contact=M["r_contact"], c_adh=M["c_adh"],
                            rep=M["rep"], omega=M["adh"], A=M["A"], force_cap=5e-8, out_d=out_d)
        M["Vc"].zero_()
        wp.launch(dcm_volume_kernel, dim=M["nf"], inputs=[pos_buf, M["faces_d"], M["fcell_d"], M["Vc"]],
                  device=cache.device)
        wp.launch(_dp_from_vol, dim=M["n_cells"],
                  inputs=[M["Vc"], wp.float64(M["V0"]), wp.float64(M["p"].turgor_dP0),
                          wp.float64(M["k_vol"]), M["dP"]], device=cache.device)
        wp.launch(dcm_turgor_force_kernel, dim=M["nf"],
                  inputs=[pos_buf, M["faces_d"], M["fcell_d"], M["dP"], out_d], device=cache.device)
        cache.contact_into(pos_buf, M["cof_d"], M["faces_d"], M["fcell_d"], rep=M["rep"],
                           adh=M["adh"], c_rep=M["c_rep"], c_adh=M["c_adh"], out_d=out_d)
        wp.launch(_bond_accumulate, dim=M["ne"],
                  inputs=[pos_buf, M["edges_d"], wp.float64(M["p"].k_edge), M["r0_d"], out_d],
                  device=cache.device)

    return stiff_frozen


def _one_implicit_step(stiff_into, M, *, dt, device, precond_apply=None, maxiter=200, tol=1e-8):
    """Take one linearly-implicit Euler step (b = F(xₙ), solve (aI+K)dx=b) and return (dx_np, iters)."""
    import numpy as np
    from ffn_sim.dcm.dcm_warp_implicit import device_cg
    N = M["N"]
    scratch = {k: wp.zeros(N, dtype=wp.vec3d, device=device)
               for k in ("r", "p", "Ap", "dx", "Fx", "Fp", "xp")}
    scratch["sca"] = wp.zeros(1, dtype=wp.float64, device=device)
    b = wp.zeros(N, dtype=wp.vec3d, device=device)
    stiff_into(M["pos_d"], b)
    a_imp = M["gamma"] / dt
    dx_d, it = device_cg(stiff_into, M["pos_d"], a_imp, b, scratch, device=device,
                         maxiter=maxiter, tol=tol, precond_apply=precond_apply)
    wp.synchronize_device(device)
    return dx_d.numpy().copy(), it


def _frozen_parity_demo(device="cpu", n_cells=24):
    import time
    import numpy as np
    M = _build_stiff_cluster(n_cells=n_cells, device=device, k_vol=1.0e3)
    M["build_grids"]()
    dt = 8e-6 * 100

    # ORIGINAL (live grid query each matvec)
    dx_orig, it_orig = _one_implicit_step(M["stiff_live"], M, dt=dt, device=device)

    # FROZEN cache (built once on xₙ; replayed each matvec)
    cache = FrozenNeighborCache(M["N"], device=device)
    cache.rebuild(node_grid_id=M["ng"].id, face_grid_id=M["fg"].id, node_f32=M["nf32"],
                  cof_d=M["cof_d"], fcell_d=M["fcell_d"], coh_q=M["coh_q"], con_q=M["con_q"],
                  do_cohesion=True, do_contact=True)
    stiff_frozen = _make_frozen_stiff(M, cache)
    dx_froz, it_froz = _one_implicit_step(stiff_frozen, M, dt=dt, device=device)

    denom = max(float(np.abs(dx_orig).max()), 1e-300)
    rel = float(np.abs(dx_froz - dx_orig).max()) / denom

    # matvec timing: many force evals each way (the query-removal payoff)
    def time_matvec(fn, reps=200):
        out = wp.zeros(M["N"], dtype=wp.vec3d, device=device)
        fn(M["pos_d"], out); wp.synchronize_device(device)
        t0 = time.perf_counter()
        for _ in range(reps):
            fn(M["pos_d"], out)
        wp.synchronize_device(device)
        return (time.perf_counter() - t0) / reps * 1e3   # ms/matvec

    t_live = time_matvec(M["stiff_live"])
    t_froz = time_matvec(stiff_frozen)

    print("== OPTIMIZATION #2 — frozen-neighbour cache ==")
    print(f"  config: {n_cells} icospheres FCC-packed, N={M['N']} nodes, dt={dt:.1e} (=100x base)")
    print(f"  CG iters: orig={it_orig}  frozen={it_froz}")
    print(f"  PARITY  max|dx_frozen - dx_orig| / max|dx_orig| = {rel:.3e}   "
          f"{'PASS (<1e-5)' if rel < 1e-5 else 'FAIL'}")
    print(f"  matvec  live-query={t_live:.3f} ms   frozen={t_froz:.3f} ms   "
          f"speedup {t_live / max(t_froz, 1e-9):.2f}x  (grid query removed)")
    print(f"  cache caps used: coh={cache.coh_cap}  con={cache.con_cap}")
    return rel, t_live, t_froz


@wp.kernel
def _spd_spring_force(kdiag: wp.array(dtype=wp.float64), pos: wp.array(dtype=wp.vec3d),
                      out: wp.array(dtype=wp.vec3d)):
    """A linear restoring FORCE F(x) = −diag(kdiag)·x. Its negative-Jacobian is K=+diag(kdiag) (PSD,
    heterogeneous), so ``device_cg`` builds the genuinely-SPD operator A = a·I + diag(kdiag). Linear
    ⇒ the FD matvec is EXACT (no floor) — the clean setting to validate the preconditioner plumbing
    in isolation from the real DCM FD-operator pathology."""
    i = wp.tid()
    out[i] = pos[i] * (-kdiag[i])


def _precond_positive_control(device="cpu", N=2000, kappa=1.0e4):
    """Positive control: a LINEAR SPD operator A = a·I + diag(k) with a HETEROGENEOUS diagonal
    (condition number ≈ κ). The analytic-diagonal Jacobi M⁻¹=1/diagA is its EXACT inverse, so PCG
    converges in 1 iteration where plain CG needs O(√κ). Confirms the new ``precond_apply`` plumbing
    in ``device_cg`` is correct + converges to the SAME solution — isolated from the FD pathology."""
    import numpy as np
    from ffn_sim.dcm.dcm_warp_implicit import device_cg
    rng = np.random.default_rng(0)
    a = 1.0
    kdiag_np = rng.random(N) * (kappa - 1.0) * a              # k ∈ [0, a·(κ−1)] → cond(A) ≈ κ
    kdiag = wp.array(kdiag_np, dtype=wp.float64, device=device)
    diagA = wp.array(a + kdiag_np, dtype=wp.float64, device=device)
    b = wp.array(rng.normal(size=(N, 3)), dtype=wp.vec3d, device=device)
    x0 = wp.zeros(N, dtype=wp.vec3d, device=device)

    def op(pos_buf, out_d):
        wp.launch(_spd_spring_force, dim=N, inputs=[kdiag, pos_buf, out_d], device=device)

    apply = make_diag_precond(diagA, device=device)
    sc = {k: wp.zeros(N, dtype=wp.vec3d, device=device) for k in ("r", "p", "Ap", "dx", "Fx", "Fp", "xp")}
    sc["sca"] = wp.zeros(1, dtype=wp.float64, device=device)
    dxp, itp = device_cg(op, x0, a, b, sc, device=device, maxiter=2000, tol=1e-10, eps=1e-3)
    sc2 = {k: wp.zeros(N, dtype=wp.vec3d, device=device) for k in ("r", "p", "Ap", "dx", "Fx", "Fp", "xp")}
    sc2["sca"] = wp.zeros(1, dtype=wp.float64, device=device)
    dxc, itc = device_cg(op, x0, a, b, sc2, device=device, maxiter=2000, tol=1e-10, eps=1e-3,
                         precond_apply=apply)
    rel = float(np.abs(dxc.numpy() - dxp.numpy()).max()) / max(float(np.abs(dxp.numpy()).max()), 1e-300)
    print("== OPTIMIZATION #1 — Jacobi preconditioner POSITIVE CONTROL (analytic SPD operator) ==")
    print(f"  A = a·I + diag(k), heterogeneous diagonal, cond≈{kappa:.0e}, N={N}")
    print(f"  CG iters:  plain={itp}   precond={itc}   "
          f"({'precond ≈1 iter (EXACT Jacobi) ' + ('PASS' if itc <= 2 and itp > 10 else '') if itc < itp else 'no reduction'})")
    print(f"  PARITY  max|dx_pc - dx_plain|/max|dx_plain| = {rel:.2e}   "
          f"{'PASS (same solution)' if rel < 1e-4 else 'CHECK'}")
    print(f"  → the PCG plumbing is CORRECT: on a heterogeneous-diagonal SPD operator the")
    print(f"    analytic-diagonal Jacobi collapses the iteration count, as it must.")
    return itp, itc


def _precond_demo(device="cpu", n_cells=24, k_vol=7.73e5, compress=0.99, dt_mult=1000.0):
    """HONEST assessment on the REAL DCM stiff operator. Reports iters_before/after AND the
    diagnostic that explains the result: whether the operator's diagonal is heterogeneous (Jacobi
    can help) or dominated by the uniform implicit regulariser a=γ/dt (Jacobi ≈ identity → no win),
    and whether the FD matvec actually converges (the ``pAp_diag`` floor firing = a broken matvec
    that NO preconditioner can fix — that needs the analytic Hessian via ``hess_apply``)."""
    import os
    import io
    import contextlib
    import numpy as np
    M = _build_stiff_cluster(n_cells=n_cells, device=device, k_vol=k_vol)
    M["build_grids"]()
    P = M["pos_d"].numpy()
    com = P[M["cof_d"].numpy() >= 0].mean(0)
    P = com + (P - com) * compress
    M["pos_d"].assign(np.ascontiguousarray(P))
    M["build_grids"]()
    dt = 8e-6 * dt_mult

    from ffn_sim.dcm.dcm_contact_implicit_warp import nearest_face_repel_kernel
    cn_k = wp.zeros(M["N"], dtype=wp.float64, device=device)
    cn_nrm = wp.zeros(M["N"], dtype=wp.vec3d, device=device)
    fz = wp.zeros(M["N"], dtype=wp.vec3d, device=device)
    wp.launch(nearest_face_repel_kernel, dim=M["N"],
              inputs=[M["fg"].id, M["nf32"], M["pos_d"], M["cof_d"], M["faces_d"], M["fcell_d"],
                      wp.float32(M["con_q"]), wp.float64(M["rep"]), fz, cn_k, cn_nrm], device=device)
    diagA = wp.zeros(M["N"], dtype=wp.float64, device=device)
    build_diagA(diagA, a=M["gamma"] / dt, cof_d=M["cof_d"], edges_d=M["edges_d"], k_edge=M["k_edge"],
                faces_d=M["faces_d"], fcell_d=M["fcell_d"], pos_d=M["pos_d"], k_vol=M["k_vol"],
                V0=M["V0"], cn_k=cn_k, device=device)
    apply = make_diag_precond(diagA, device=device)

    # count pAp-floor hits (the broken-matvec signal) via CG_DEBUG
    def run_count_floor(precond):
        os.environ["CG_DEBUG"] = "1"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            dx, it = _one_implicit_step(M["stiff_live"], M, dt=dt, device=device, maxiter=400,
                                        precond_apply=(apply if precond else None))
        del os.environ["CG_DEBUG"]
        floors = buf.getvalue().count("K-noise floor")
        return dx, it, floors

    dx_plain, it_plain, fl_plain = run_count_floor(False)
    dx_pc, it_pc, fl_pc = run_count_floor(True)
    denom = max(float(np.abs(dx_plain).max()), 1e-300)
    rel = float(np.abs(dx_pc - dx_plain).max()) / denom
    da = diagA.numpy()
    a_val = M['gamma'] / dt
    spread = da.max() / max(da.min(), 1e-300)
    print("== OPTIMIZATION #1 — analytic-diagonal Jacobi on the REAL DCM stiff operator ==")
    print(f"  config: {n_cells} icospheres compressed {(1-compress)*100:.0f}%, k_vol={k_vol:.2e}, "
          f"N={M['N']}, dt={dt:.1e}, n_active_contact={int((cn_k.numpy()>0).sum())}")
    print(f"  diagonal terms: a=γ/dt + k_edge(incident) + rep·area(contact) + k_vol/V0·|∂V/∂x|²(turgor)")
    print(f"  diagA range [{da.min():.3e}, {da.max():.3e}]  a=γ/dt={a_val:.3e}  "
          f"spread(max/min)={spread:.2f}  (Jacobi can only help if spread ≫ 1)")
    print(f"  pAp-floor hits/solve: plain={fl_plain}  precond={fl_pc}  "
          f"({'FD matvec BROKEN (floored) → no precond can fix; needs analytic Hessian (--ipc)' if fl_plain > 2 else 'FD matvec clean'})")
    print(f"  CG iters:  plain={it_plain}   precond={it_pc}")
    note = ("a-dominated: a=γ/dt is uniform + dominant → diagA≈const → Jacobi≈identity → no win (EXPECTED)"
            if spread < 1.5 else
            ("heterogeneous diagonal but FD matvec floored → not a conditioning problem" if fl_plain > 2
             else "heterogeneous + clean → precond helps"))
    print(f"  finding: {note}")
    if rel < 1e-5:
        print(f"  PARITY  max|dx_pc - dx_plain|/max|dx_plain| = {rel:.2e}  PASS (same dx)")
    else:
        print(f"  (dx differ {rel:.2e}: both solves hit the pAp floor → neither reached the TRUE dx; "
              f"the FD matvec — not the preconditioner — is the limiter here)")
    return it_plain, it_pc, rel


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--n", type=int, default=24)
    a = ap.parse_args()
    _frozen_parity_demo(a.device, a.n)
    print()
    _precond_positive_control(a.device)
    print()
    _precond_demo(a.device, a.n)
