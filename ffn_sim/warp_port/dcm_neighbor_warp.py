"""Hash-grid neighbour-list acceleration for the multi-cell DCM inter-cell forces.

The committed ``dcm_cohesion_warp`` (O(N^2)) and ``dcm_contact_warp`` (O(N*M),
closest-point per pair) are brute-force — written for PARITY, and the multi-cell
A5000 bottleneck. These grid kernels reproduce the SAME force laws but iterate
only the local neighbours from a Warp ``HashGrid`` -> O(N*k):

  * cohesion: grid over node positions; node i queries node neighbours within c_adh.
  * contact : grid over FACE CENTROIDS; node i queries faces whose centroid is
    within c_adh + face_reach (so no face with a closest-point in range is missed),
    then the exact ``closest_bary`` law on those candidates only.

The grid is built on a float32 copy of the positions each step (positions ~1e-6 m
=> cell index floor(p/radius) is O(1), float32-safe); the force math uses the f64
positions, so the grid only prunes far pairs that contribute exactly 0 — the
result equals the brute-force force to summation-order tolerance.
"""

from __future__ import annotations

import warp as wp

from ffn_sim.warp_port.dcm_contact_warp import closest_bary

wp.init()


@wp.kernel
def pos_to_f32(pos: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3)):
    i = wp.tid()
    p = pos[i]
    out[i] = wp.vec3(wp.float32(p[0]), wp.float32(p[1]), wp.float32(p[2]))


@wp.kernel
def face_centroids_f32(pos: wp.array(dtype=wp.vec3d),
                       faces: wp.array(dtype=wp.int32, ndim=2),
                       out: wp.array(dtype=wp.vec3)):
    f = wp.tid()
    a = pos[faces[f, 0]]
    b = pos[faces[f, 1]]
    c = pos[faces[f, 2]]
    cx = (a[0] + b[0] + c[0]) / wp.float64(3.0)
    cy = (a[1] + b[1] + c[1]) / wp.float64(3.0)
    cz = (a[2] + b[2] + c[2]) / wp.float64(3.0)
    out[f] = wp.vec3(wp.float32(cx), wp.float32(cy), wp.float32(cz))


@wp.kernel
def cohesion_grid_kernel(
    grid: wp.uint64,
    qpts: wp.array(dtype=wp.vec3),          # f32 node positions (grid query points)
    pos: wp.array(dtype=wp.vec3d),          # f64 node positions (force math)
    cof: wp.array(dtype=wp.int32),
    radius: wp.float32,                      # query radius (= c_adh) in f32
    r_contact: wp.float64, c_adh: wp.float64,
    rep: wp.float64, omega: wp.float64, A: wp.float64, force_cap: wp.float64,
    force: wp.array(dtype=wp.vec3d),         # (N,) own-row write
):
    i = wp.tid()
    c1 = cof[i]
    if c1 < wp.int32(0):
        return
    z = wp.float64(0.0)
    half = wp.float64(0.5) * c_adh
    pi = pos[i]
    acc = wp.vec3d(z, z, z)
    q = wp.hash_grid_query(grid, qpts[i], radius)
    j = wp.int32(0)
    while wp.hash_grid_query_next(q, j):
        if j != i and cof[j] >= wp.int32(0) and cof[j] != c1:
            r_vec = pi - pos[j]
            d = wp.length(r_vec)
            fmag = z
            if d < r_contact:                                  # REPULSION
                fmag = rep * A * (r_contact - d)
            else:
                if d < c_adh:                                  # ADHESION tent
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
def contact_grid_kernel(
    grid: wp.uint64,                         # grid over FACE centroids
    qpts: wp.array(dtype=wp.vec3),           # f32 node positions (query points)
    pos: wp.array(dtype=wp.vec3d),           # f64 node positions
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    radius: wp.float32,                      # = c_adh + face_reach (f32)
    rep: wp.float64, adh: wp.float64, c_rep: wp.float64, c_adh: wp.float64,
    force: wp.array(dtype=wp.vec3d),         # atomic accumulate
):
    ni = wp.tid()
    c1 = cof[ni]
    if c1 < wp.int32(0):
        return
    z = wp.float64(0.0)
    half = wp.float64(0.5) * c_adh
    p = pos[ni]
    fn_acc = wp.vec3d(z, z, z)
    q = wp.hash_grid_query(grid, qpts[ni], radius)
    fj = wp.int32(0)
    while wp.hash_grid_query_next(q, fj):
        if fcell[fj] != c1:
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
def gather_lead_pos(pos: wp.array(dtype=wp.vec3d), lead_idx: wp.array(dtype=wp.int32),
                    out: wp.array(dtype=wp.vec3d)):
    """Device gather lead_rp[t] = pos[lead_idx[t]] (per-step; leading-node SET is
    clutch-gripped/low-cadence, but their POSITIONS move every step)."""
    t = wp.tid()
    out[t] = pos[lead_idx[t]]


@wp.kernel
def lamellipodium_tether_accum(
    lead_rp: wp.array(dtype=wp.vec3d), lead_ccx: wp.array(dtype=wp.float64),
    lead_ccy: wp.array(dtype=wp.float64), lead_ox: wp.array(dtype=wp.float64),
    lead_oy: wp.array(dtype=wp.float64), lead_proj: wp.array(dtype=wp.float64),
    lead_idx: wp.array(dtype=wp.int32), actin: wp.array(dtype=wp.vec3d),
    n_actin: wp.int32, k_tether: wp.float64, force_cap: wp.float64,
    tether_radius: wp.float64, force: wp.array(dtype=wp.vec3d),
):
    """ACCUMULATING variant of lamellipodium_tether_kernel (atomic_add, not overwrite)
    so it composes with turgor/contact/cohesion/edges in the shared force array."""
    t = wp.tid()
    rp = lead_rp[t]
    ccx = lead_ccx[t]
    ccy = lead_ccy[t]
    ox = lead_ox[t]
    oy = lead_oy[t]
    npj = lead_proj[t]
    min_dist = wp.float64(1.0e300)
    jmin = wp.int32(-1)
    for a in range(n_actin):
        ap = actin[a]
        dx = ap[0] - rp[0]
        dy = ap[1] - rp[1]
        dz = ap[2] - rp[2]
        dist = wp.sqrt(dx * dx + dy * dy + dz * dz)
        act_proj = (ap[0] - ccx) * ox + (ap[1] - ccy) * oy
        if dist <= tether_radius and act_proj > npj:
            if dist < min_dist:
                min_dist = dist
                jmin = wp.int32(a)
    if jmin >= wp.int32(0):
        ap = actin[jmin]
        fvx = k_tether * (ap[0] - rp[0])
        fvy = k_tether * (ap[1] - rp[1])
        fvz = k_tether * (ap[2] - rp[2])
        fmag = wp.sqrt(fvx * fvx + fvy * fvy + fvz * fvz)
        if fmag > force_cap:
            s = force_cap / fmag
            fvx = fvx * s
            fvy = fvy * s
            fvz = fvz * s
        wp.atomic_add(force, lead_idx[t], wp.vec3d(fvx, fvy, fvz))
