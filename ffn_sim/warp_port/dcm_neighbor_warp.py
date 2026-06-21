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
def edge_midpoints_f32(pos: wp.array(dtype=wp.vec3d),
                       edges: wp.array(dtype=wp.int32, ndim=2),
                       out: wp.array(dtype=wp.vec3)):
    e = wp.tid()
    a = pos[edges[e, 0]]
    b = pos[edges[e, 1]]
    h = wp.float64(0.5)
    out[e] = wp.vec3(wp.float32(h * (a[0] + b[0])),
                     wp.float32(h * (a[1] + b[1])),
                     wp.float32(h * (a[2] + b[2])))


@wp.func
def _clamp01(x: wp.float64) -> wp.float64:
    if x < wp.float64(0.0):
        return wp.float64(0.0)
    if x > wp.float64(1.0):
        return wp.float64(1.0)
    return x


@wp.kernel
def edge_edge_contact_kernel(
    grid: wp.uint64,                         # grid over EDGE midpoints
    mids: wp.array(dtype=wp.vec3),           # f32 edge midpoints (query points)
    pos: wp.array(dtype=wp.vec3d),
    edges: wp.array(dtype=wp.int32, ndim=2),
    edge_cell: wp.array(dtype=wp.int32),     # per-edge owner cell
    radius: wp.float32,                      # query radius (= c_rep + max edge half-length)
    c_rep: wp.float64, rep: wp.float64, A: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """A2 edge-edge excluded volume: the case node-vs-face MISSES — a triangle edge of one
    cell sliding between the (sparse) nodes of another. One thread per edge e1: find the
    nearest OTHER-cell edge e2 (segment-segment closest point, Ericson), and if the segments
    pass within ``c_rep`` push e1's two endpoints away from e2 by ``rep·A·(c_rep−d)`` split
    by the closest-point barycentric. Each edge pushes only ITS OWN endpoints (e2's thread
    pushes e2) → symmetric, no double count, momentum-clean. Same-cell edges are skipped."""
    e1 = wp.tid()
    c1 = edge_cell[e1]
    if c1 < wp.int32(0):
        return
    p1 = pos[edges[e1, 0]]
    q1 = pos[edges[e1, 1]]
    d1 = q1 - p1
    a = wp.dot(d1, d1)
    z = wp.float64(0.0)
    eps = wp.float64(1.0e-30)
    best_d = wp.float64(1.0e300)
    best_s = z
    best_dir = wp.vec3d(z, z, z)
    q = wp.hash_grid_query(grid, mids[e1], radius)
    e2 = wp.int32(0)
    while wp.hash_grid_query_next(q, e2):
        if e2 != e1 and edge_cell[e2] != c1 and edge_cell[e2] >= wp.int32(0):
            p2 = pos[edges[e2, 0]]
            q2 = pos[edges[e2, 1]]
            d2 = q2 - p2
            r = p1 - p2
            e = wp.dot(d2, d2)
            f = wp.dot(d2, r)
            s = z
            t = z
            if a <= eps and e <= eps:
                s = z
                t = z
            elif a <= eps:
                t = _clamp01(f / e)
            else:
                cc = wp.dot(d1, r)
                if e <= eps:
                    s = _clamp01(-cc / a)
                else:
                    b = wp.dot(d1, d2)
                    denom = a * e - b * b
                    if denom > eps or denom < -eps:
                        s = _clamp01((b * f - cc * e) / denom)
                    t = (b * s + f) / e
                    if t < z:
                        t = z
                        s = _clamp01(-cc / a)
                    elif t > wp.float64(1.0):
                        t = wp.float64(1.0)
                        s = _clamp01((b - cc) / a)
            cp1 = p1 + d1 * s
            cp2 = p2 + d2 * t
            sep = cp1 - cp2
            dist = wp.length(sep)
            if dist < best_d:
                best_d = dist
                best_s = s
                best_dir = sep
    if best_d < c_rep:
        dlen = best_d
        if dlen < eps:
            dlen = eps
        amp = rep * A * (c_rep - best_d)
        fvec = best_dir * (amp / dlen)               # push e1 AWAY from e2
        wp.atomic_add(force, edges[e1, 0], fvec * (wp.float64(1.0) - best_s))
        wp.atomic_add(force, edges[e1, 1], fvec * best_s)


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
def cohesion_grid_cad_kernel(
    grid: wp.uint64,
    qpts: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32),
    cad: wp.array(dtype=wp.float64),         # (n_cells,) per-cell cadherin multiplier
    radius: wp.float32,
    r_contact: wp.float64, c_adh: wp.float64,
    rep: wp.float64, omega: wp.float64, A: wp.float64, force_cap: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """M3 junction-switch variant of :func:`cohesion_grid_kernel`: the cell-cell
    ADHESION (only) is scaled by ``mult = sqrt(cad[c1]·cad[c2])`` so a crowd-switched
    cell (``cad`` lowered to ``cadherin_weak_factor``) de-coheres from its neighbours.
    Repulsion is NOT modulated (excluded volume is unconditional)."""
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
            if d < r_contact:                                  # REPULSION (unmodulated)
                fmag = rep * A * (r_contact - d)
            else:
                if d < c_adh:                                  # ADHESION tent (× cad)
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
def contact_grid_cad_kernel(
    grid: wp.uint64,
    qpts: wp.array(dtype=wp.vec3),
    pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    cad: wp.array(dtype=wp.float64),         # (n_cells,) per-cell cadherin multiplier
    radius: wp.float32,
    rep: wp.float64, adh: wp.float64, c_rep: wp.float64, c_adh: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """M3 junction-switch variant of :func:`contact_grid_kernel`: the node-face
    ADHESION (only) is scaled by ``mult = sqrt(cad[c1]·cad[c_face])``. Repulsion
    (excluded volume) is NOT modulated."""
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
                amp = rep * area                               # REPULSION (unmodulated)
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


@wp.kernel
def penetration_depth_kernel(
    grid: wp.uint64,                         # grid over FACE centroids
    qpts: wp.array(dtype=wp.vec3),           # f32 node positions (query points)
    pos: wp.array(dtype=wp.vec3d),           # f64 node positions
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    radius: wp.float32,
    pen: wp.array(dtype=wp.float64),         # (N,) out: penetration depth (0 = outside)
):
    """DIAGNOSTIC (no force): per node, the depth it has penetrated INTO another cell.

    Finds the node's nearest face on any OTHER cell (Ericson closest-point, same geometry
    as ``contact_grid_kernel``); if the node sits on the INNER side of that face
    (``sign = (node−closestpoint)·faceNormal < 0``), it is inside that cell and the
    penetration depth is the closest-point distance. ``pen[ni]`` = that depth (0 if outside
    or no neighbour). The driver reduces ``max(pen)/mean_edge`` per frame — interpenetration
    is then tracked as a number (it is INVISIBLE in the surface render), never eyeballed."""
    ni = wp.tid()
    c1 = cof[ni]
    if c1 < wp.int32(0):
        pen[ni] = wp.float64(0.0)
        return
    z = wp.float64(0.0)
    p = pos[ni]
    best_d = wp.float64(1.0e300)
    best_sign = wp.float64(1.0)
    q = wp.hash_grid_query(grid, qpts[ni], radius)
    fj = wp.int32(0)
    while wp.hash_grid_query_next(q, fj):
        if fcell[fj] != c1:
            a = pos[faces[fj, 0]]
            b = pos[faces[fj, 1]]
            c = pos[faces[fj, 2]]
            bary = closest_bary(p, a, b, c)
            cpa = a * bary[0] + b * bary[1] + c * bary[2]
            r_vec = p - cpa
            d = wp.length(r_vec)
            if d < best_d:
                best_d = d
                fnv = wp.cross(b - a, c - a)
                nrm = wp.length(fnv)
                if nrm > z:
                    best_sign = wp.dot(r_vec, fnv) / nrm
                else:
                    best_sign = wp.float64(1.0)
    if best_d < wp.float64(1.0e299) and best_sign < z:
        pen[ni] = best_d
    else:
        pen[ni] = z


@wp.kernel
def cadherin_bond_force_kernel(
    bonds: wp.array(dtype=wp.vec2i),         # (M,) node-index pairs (i in cell A, j in cell B)
    n_bonds: wp.int32,
    pos: wp.array(dtype=wp.vec3d),
    k_trans: wp.float64, r0_trans: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """E1 explicit cadherin trans-dimer force: ATTRACTIVE-only harmonic tether. For a
    stretched bond (L > r0_trans) the ectodomain bridge pulls its two membrane nodes
    together with ``F = k_trans·(L − r0_trans)``; a slack bond (L ≤ r0) is force-free (a
    floppy tether doesn't push). Excluded-volume repulsion is the cohesion/contact kernels'
    job (run adhesion-OFF in cadherin mode), so this is the SOLE cell-cell adhesion — and
    de-cohesion is emergent: the host breaks each bond at the Rakshit catch-slip rate."""
    t = wp.tid()
    if t >= n_bonds:
        return
    e = bonds[t]
    i = e[0]
    j = e[1]
    rij = pos[j] - pos[i]
    L = wp.length(rij)
    if L > r0_trans and L > wp.float64(1.0e-30):
        amp = k_trans * (L - r0_trans) / L
        fvec = rij * amp                          # on i toward j (pull together)
        wp.atomic_add(force, i, fvec)
        wp.atomic_add(force, j, -fvec)


@wp.kernel
def edge_neighbor_sum_kernel(
    vals: wp.array(dtype=wp.vec3d), edges: wp.array(dtype=wp.int32, ndim=2),
    vsum: wp.array(dtype=wp.vec3d), vcnt: wp.array(dtype=wp.float64),
):
    """Accumulate each node's 1-ring neighbour sum of ``vals`` over the mesh edges (zero
    vsum/vcnt first). Symmetric: edge (i,j) adds vals[j]→i and vals[i]→j."""
    e = wp.tid()
    i = edges[e, 0]; j = edges[e, 1]
    wp.atomic_add(vsum, i, vals[j]); wp.atomic_add(vcnt, i, wp.float64(1.0))
    wp.atomic_add(vsum, j, vals[i]); wp.atomic_add(vcnt, j, wp.float64(1.0))


@wp.kernel
def umbrella_kernel(
    vals: wp.array(dtype=wp.vec3d), vsum: wp.array(dtype=wp.vec3d),
    vcnt: wp.array(dtype=wp.float64), out: wp.array(dtype=wp.vec3d),
):
    """Discrete umbrella Laplacian: out[i] = mean(neighbours of vals) − vals[i] (0 if isolated)."""
    i = wp.tid()
    n = vcnt[i]
    if n > wp.float64(0.0):
        out[i] = vsum[i] / n - vals[i]
    else:
        out[i] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))


@wp.kernel
def bending_apply_kernel(
    bilap: wp.array(dtype=wp.vec3d), cof: wp.array(dtype=wp.int32),
    k_bend: wp.float64, force: wp.array(dtype=wp.vec3d),
):
    """B5 thin-plate (Helfrich-like) bending: F = −k_bend·Δ²r (biharmonic = umbrella-of-umbrella).
    Resists curvature VARIATION (flat and uniformly-curved are both low-energy) — true bending,
    distinct from the area-gradient surface tension. Skips dormant nodes (cof<0)."""
    i = wp.tid()
    if cof[i] >= wp.int32(0):
        wp.atomic_add(force, i, bilap[i] * (-k_bend))


@wp.func
def _tri_area_grads(a: wp.vec3d, b: wp.vec3d, c: wp.vec3d):
    """Per-vertex area gradient of triangle (a,b,c): ∂A/∂a = ½ n̂×(c−b), etc. Returns the three
    gradient vectors packed; magnitude ½|opposite edge|, in-plane, pointing to INCREASE area."""
    nrm = wp.cross(b - a, c - a)
    L = wp.length(nrm)
    z = wp.float64(0.0)
    if L <= wp.float64(1.0e-30):
        return wp.vec3d(z, z, z), wp.vec3d(z, z, z), wp.vec3d(z, z, z)
    nh = nrm / L
    h = wp.float64(0.5)
    ga = wp.cross(nh, c - b) * h
    gb = wp.cross(nh, a - c) * h
    gc = wp.cross(nh, b - a) * h
    return ga, gb, gc


@wp.kernel
def surface_tension_kernel(
    pos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
    gamma: wp.float64, force: wp.array(dtype=wp.vec3d),
):
    """B4 membrane surface tension: F = −γ·∂A/∂r per face, scattered to its 3 vertices — the
    area-minimising (rounding) tension of the cortex/membrane (SimuCell3D shell energy term)."""
    f = wp.tid()
    ia = faces[f, 0]; ib = faces[f, 1]; ic = faces[f, 2]
    ga, gb, gc = _tri_area_grads(pos[ia], pos[ib], pos[ic])
    wp.atomic_add(force, ia, ga * (-gamma))
    wp.atomic_add(force, ib, gb * (-gamma))
    wp.atomic_add(force, ic, gc * (-gamma))


@wp.kernel
def face_area_accum_kernel(
    pos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32), acell: wp.array(dtype=wp.float64),
):
    """Scatter each face area into its owner cell (zero acell first) → per-cell surface area."""
    f = wp.tid()
    nrm = wp.cross(pos[faces[f, 1]] - pos[faces[f, 0]], pos[faces[f, 2]] - pos[faces[f, 0]])
    wp.atomic_add(acell, fcell[f], wp.float64(0.5) * wp.length(nrm))


@wp.kernel
def global_area_force_kernel(
    pos: wp.array(dtype=wp.vec3d), faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32), acell: wp.array(dtype=wp.float64),
    a0cell: wp.array(dtype=wp.float64), k_a: wp.float64, force: wp.array(dtype=wp.vec3d),
):
    """B4 global area constraint: F = −k_a·(A_cell − A0_cell)·∂A/∂r per face → resists the
    cell's TOTAL surface area drifting from A0 (the membrane-reservoir / areal-incompressibility
    term that complements turgor's volume constraint)."""
    f = wp.tid()
    c = fcell[f]
    dev = k_a * (acell[c] - a0cell[c])
    ia = faces[f, 0]; ib = faces[f, 1]; ic = faces[f, 2]
    ga, gb, gc = _tri_area_grads(pos[ia], pos[ib], pos[ic])
    wp.atomic_add(force, ia, ga * (-dev))
    wp.atomic_add(force, ib, gb * (-dev))
    wp.atomic_add(force, ic, gc * (-dev))


@wp.kernel
def cell_centroid_accum_kernel(
    pos: wp.array(dtype=wp.vec3d), cof: wp.array(dtype=wp.int32),
    csum: wp.array(dtype=wp.vec3d), ccnt: wp.array(dtype=wp.float64),
):
    """Scatter membrane-node positions into per-cell sum + count (zero csum/ccnt first)."""
    i = wp.tid()
    c = cof[i]
    if c >= wp.int32(0):
        wp.atomic_add(csum, c, pos[i])
        wp.atomic_add(ccnt, c, wp.float64(1.0))


@wp.kernel
def nucleus_force_kernel(
    pos: wp.array(dtype=wp.vec3d), cof: wp.array(dtype=wp.int32),
    csum: wp.array(dtype=wp.vec3d), ccnt: wp.array(dtype=wp.float64),
    R_nuc: wp.float64, d_knee: wp.float64, k_chrom: wp.float64, k_lamin: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """E2 nucleus as a deformable stiff core (H.9 KU-3.B2). The DCM shell is hollow, so the
    nucleus is the radial confinement a membrane node feels when it intrudes within R_nuc of
    its cell centroid (the cell thinning below the nuclear size): an OUTWARD bilinear push —
    soft chromatin (k_chrom) up to the lamin knee (d_knee = knee_strain·R_nuc), stiff lamin-A
    shell (k_chrom+k_lamin) beyond. k_chrom = 4π·E_nuc·R_nuc/npc (continuum bridge, N-invariant);
    k_lamin = (ratio_lamin−1)·k_chrom. Force-free at r = R_nuc; never pulls inward (excluded
    volume of the nuclear core) → sets the cell-thinning / maxZ floor as it spreads."""
    i = wp.tid()
    c = cof[i]
    if c < wp.int32(0):
        return
    cen = csum[c] / ccnt[c]
    rv = pos[i] - cen
    r = wp.length(rv)
    if r < R_nuc and r > wp.float64(1.0e-30):
        d = R_nuc - r                         # intrusion depth into the nucleus (>0)
        if d <= d_knee:
            fmag = k_chrom * d
        else:
            fmag = k_chrom * d_knee + (k_chrom + k_lamin) * (d - d_knee)
        wp.atomic_add(force, i, rv * (fmag / r))   # push OUTWARD


@wp.kernel
def ecm_clutch_force_kernel(
    node_idx: wp.array(dtype=wp.int32),      # (M,) basal node of each engaged clutch
    anchor: wp.array(dtype=wp.vec3d),        # (M,) fixed substrate site it gripped
    n_clutch: wp.int32,
    k_fa: wp.float64,
    pos: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
):
    """C6 integrin-ECM clutch force: harmonic spring from the basal node to the substrate site
    it gripped, ``F = k_fa·(anchor − r_node)``. The anchor is the rigid dish ligand (no
    reaction), so the clutch RESISTS the node sliding off its grip — transmitting the
    lamellipodium/membrane pull to the substrate as traction. The host breaks each clutch at
    the Pereverzev catch-slip rate (catch near F*, slip at F≫F_s) → traction-limited spread."""
    t = wp.tid()
    if t >= n_clutch:
        return
    i = node_idx[t]
    d = anchor[t] - pos[i]
    wp.atomic_add(force, i, d * k_fa)


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


@wp.kernel
def lamellipodium_tether_multicell(
    lead_rp: wp.array(dtype=wp.vec3d), lead_ccx: wp.array(dtype=wp.float64),
    lead_ccy: wp.array(dtype=wp.float64), lead_ox: wp.array(dtype=wp.float64),
    lead_oy: wp.array(dtype=wp.float64), lead_proj: wp.array(dtype=wp.float64),
    lead_idx: wp.array(dtype=wp.int32), lead_cell: wp.array(dtype=wp.int32),
    actin: wp.array(dtype=wp.vec3d), actin_cell: wp.array(dtype=wp.int32),
    n_actin: wp.int32, k_tether: wp.float64, force_cap: wp.float64,
    tether_radius: wp.float64, force: wp.array(dtype=wp.vec3d),
):
    """Multi-cell traction tether (M2): pull each rim cell's leading basal node
    toward the nearest OUTWARD clutch-anchored actin bead **of its OWN cell**.

    Two fixes over :func:`lamellipodium_tether_accum` (M2b §E): (a) a **same-cell
    mask** (``actin_cell[a] == lead_cell[t]``) so a node only tethers to its own
    cell's front — required once the actin pool is shared across all rim cells;
    (b) **no Newton-3 reaction** on the actin bead — the host pins the clutch-gripped
    actin as a rigid anchor (slip ≤ tether_cap/k_clutch = 0.1µm ≤ 0.2·ℓ₀), so the
    bead is the fixed point and is never integrated. ``F = +k_tether·(r_actin − r_node)``
    capped at ``force_cap`` (= tether_cap), accumulated into the shared force array."""
    t = wp.tid()
    rp = lead_rp[t]
    ccx = lead_ccx[t]
    ccy = lead_ccy[t]
    ox = lead_ox[t]
    oy = lead_oy[t]
    npj = lead_proj[t]
    lc = lead_cell[t]
    min_dist = wp.float64(1.0e300)
    jmin = wp.int32(-1)
    for a in range(n_actin):
        if actin_cell[a] == lc:
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
