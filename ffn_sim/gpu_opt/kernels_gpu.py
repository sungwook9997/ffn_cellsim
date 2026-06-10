"""cupy port of kernels_cpu.py — 4 particle/mesh force kernels for GPU (RTX A5000).

Same function names/signatures as the NumPy reference. Each kernel accepts numpy or
cupy arrays in and returns a cupy (N,3) float64 array. The per-element math is kept
identical to the reference (parity rtol 1e-9). K4 replaces the O(N^2) reference loop
with a fully vectorized uniform-grid neighbour list (no RawKernels, no per-particle
Python loops — only a fixed 27-iteration loop over neighbour-cell offsets).
"""

from __future__ import annotations

import cupy as cp
import cupyx


# ---------------------------------------------------------------------------
# K1. mesh_pressure_forces — closed-mesh volume-conserving normal pressure
# ---------------------------------------------------------------------------
def mesh_pressure_forces(pos, tris, tri_group, n_groups, V0, p0, K):
    """Outward facet-normal pressure on closed triangulated mesh-groups. Returns (N,3)."""
    p = cp.asarray(pos)
    t = cp.asarray(tris)
    tg = cp.asarray(tri_group)
    v0, v1, v2 = p[t[:, 0]], p[t[:, 1]], p[t[:, 2]]
    a = v1 - v0
    b = v2 - v0
    # manual cross product (same FP ops as np.cross: multiply, multiply, subtract)
    cx = a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1]
    cy = a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2]
    cz = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    cross = cp.stack((cx, cy, cz), axis=1)                   # 2·area·n̂
    vol_contrib = cp.einsum("ij,ij->i", v0, cross) / 6.0
    Vg = cp.zeros(n_groups, dtype=cp.float64)
    cupyx.scatter_add(Vg, tg, vol_contrib)
    P = p0 + K * (V0 - Vg) / V0
    Pf = P[tg]
    f_facet = (Pf[:, None] * cross / 2.0) / 3.0
    F = cp.zeros_like(p)
    for col in range(3):
        cupyx.scatter_add(F, t[:, col], f_facet)
    return F


# ---------------------------------------------------------------------------
# K2. plane_well_forces — capped-harmonic well about the plane z = z0
# ---------------------------------------------------------------------------
def plane_well_forces(pos, z0, k, w, mask=None):
    """z-only capped-harmonic well. |dz|<=w: F_z=-k·dz; dz<-w: F_z=+k·w; else 0."""
    p = cp.asarray(pos)
    dz = p[:, 2] - z0
    F = cp.zeros_like(p)
    within = cp.abs(dz) <= w
    deep = dz < -w                                           # disjoint from `within`
    F[:, 2] = cp.where(within, -k * dz, cp.where(deep, k * w, 0.0))
    if mask is not None:
        m = cp.asarray(mask)
        F = cp.where(m[:, None], F, 0.0)
    return F


# ---------------------------------------------------------------------------
# K3. bond_spring_forces — per-bond Hookean springs between index pairs
# ---------------------------------------------------------------------------
def bond_spring_forces(pos, bonds, r0, k):
    """Per-bond Hookean spring; ± scatter to the two endpoints. Returns (N,3)."""
    p = cp.asarray(pos)
    F = cp.zeros_like(p)
    bnd = cp.asarray(bonds)
    if bnd.shape[0] == 0:
        return F
    r0c = cp.asarray(r0)
    a, b = bnd[:, 0], bnd[:, 1]
    dr = p[a] - p[b]
    r = cp.sqrt(cp.sum(dr * dr, axis=1))
    r_safe = cp.where(r > 1e-18, r, 1.0)
    fmag = k * (r - r0c)
    fvec = (-fmag / r_safe)[:, None] * dr
    cupyx.scatter_add(F, a, fvec)
    cupyx.scatter_add(F, b, -fvec)
    return F


# ---------------------------------------------------------------------------
# K4. group_pair_forces — inter-group pairwise soft-core + attractive well
#     (uniform-grid neighbour list; replaces the O(N^2) reference loop)
# ---------------------------------------------------------------------------
def group_pair_forces(pos, group_id, sigma, r_cut, k_core, f_well0):
    """Pairwise potential between particles of DIFFERENT groups within r_cut.

        r < sigma          : F_rep = k_core·(sigma − r)
        sigma <= r < r_cut : F_att = −f_well0·(r_cut − r)/(r_cut − sigma)
        r >= r_cut         : 0
    group_id (N,) int, −1 = inactive. Exact Newton-3 ± scatter. Returns (N,3).
    """
    p = cp.asarray(pos)
    gid = cp.asarray(group_id)
    F = cp.zeros_like(p)
    active = cp.flatnonzero(gid >= 0)
    if int(active.size) < 2:
        return F
    ap = p[active]
    ag = gid[active]
    n = int(active.size)

    # --- bin active particles on a uniform grid of edge r_cut -------------
    mn = ap.min(axis=0)
    cell = cp.floor((ap - mn) / r_cut).astype(cp.int64)
    dims = cell.max(axis=0) + 1
    nx, ny, nz = (int(dims[0]), int(dims[1]), int(dims[2]))
    cid = (cell[:, 0] * ny + cell[:, 1]) * nz + cell[:, 2]

    order = cp.argsort(cid)
    cid_sorted = cid[order]
    # NOTE: no dense per-cell start/end tables — those would be O(nx*ny*nz),
    # which scales with domain volume, not particle count (OOM on sparse
    # inputs). Instead each offset looks up its neighbour cell directly in
    # the length-n sorted cell-id array via searchsorted (O(n log n)).

    dims_arr = cp.asarray([nx, ny, nz], dtype=cp.int64)

    # --- candidate pairs from the 27 neighbour-cell offsets ---------------
    i_parts, j_parts = [], []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                off = cp.asarray([dx, dy, dz], dtype=cp.int64)
                ncoord = cell + off
                valid = cp.logical_and(
                    (ncoord >= 0).all(axis=1), (ncoord < dims_arr).all(axis=1))
                ncid = (ncoord[:, 0] * ny + ncoord[:, 1]) * nz + ncoord[:, 2]
                ncid = cp.where(valid, ncid, 0)
                st = cp.searchsorted(cid_sorted, ncid, side="left")
                en = cp.searchsorted(cid_sorted, ncid, side="right")
                st = cp.where(valid, st, 0)
                cnt = cp.where(valid, en - st, 0)
                total = int(cnt.sum())
                if total == 0:
                    continue
                # Variable-length segment expansion WITHOUT cp.repeat(x, ndarray)
                # (released cupy <= 14.0.x rejects ndarray `repeats`): for each
                # output slot k, the source particle is the segment whose
                # exclusive-cumsum bin contains k.
                cum = cp.cumsum(cnt)
                k = cp.arange(total, dtype=cp.int64)
                src = cp.searchsorted(cum, k, side="right")
                pair_i = src                                 # == arange(n)[src]
                within = k - (cum[src] - cnt[src])
                j_sortedpos = st[src] + within
                pair_j = order[j_sortedpos]
                i_parts.append(pair_i)
                j_parts.append(pair_j)

    if not i_parts:
        return F
    pi = cp.concatenate(i_parts)
    pj = cp.concatenate(j_parts)

    # --- filter: i<j dedup, different groups, within cutoff ---------------
    keep = pi < pj
    pi, pj = pi[keep], pj[keep]
    keep = ag[pi] != ag[pj]
    pi, pj = pi[keep], pj[keep]
    d = ap[pj] - ap[pi]
    r = cp.sqrt(cp.sum(d * d, axis=1))
    keep = cp.logical_and(r < r_cut, r > 1e-18)
    if int(keep.sum()) == 0:
        return F
    pi, pj, d, r = pi[keep], pj[keep], d[keep], r[keep]

    # --- per-pair math, identical to the reference -------------------------
    rep = r < sigma
    fmag = cp.where(
        rep,
        k_core * cp.clip(sigma - r, 0.0, sigma),
        -f_well0 * (r_cut - r) / (r_cut - sigma),
    )
    fvec = fmag[:, None] * (d / r[:, None])
    cupyx.scatter_add(F, active[pj], fvec)                   # +f on j
    cupyx.scatter_add(F, active[pi], -fvec)                  # −f on i (Newton-3 exact)
    return F


# ---------------------------------------------------------------------------
# K5. tent_contact_forces — SimuCell3D bilinear-tent cell-cell contact
#     (uniform-grid neighbour list; same per-pair law as kernels_cpu.K5)
# ---------------------------------------------------------------------------
def tent_contact_forces(pos, cell_id, r_contact, c_adh, rep_strength,
                        adh_strength, patch_area, force_cap, cad_mult=None):
    """SimuCell3D bilinear-tent inter-cell contact, cupy uniform-grid version.

    Per-pair force law is identical to ``kernels_cpu.tent_contact_forces`` (and to
    ``cell/dcm_contact.py::DcmTentContact``); the O(N²) inter-node search is
    replaced by the same uniform-grid neighbour list as K4 (bin on a grid of edge
    ``r_search = max(r_contact, c_adh)``, 27-cell stencil, ``i<j`` dedup). Returns
    per-node force (N,3). cad_mult (n_cells,) or None.
    """
    p = cp.asarray(pos)
    cid = cp.asarray(cell_id)
    F = cp.zeros_like(p)
    active = cp.flatnonzero(cid >= 0)
    if int(active.size) < 2:
        return F
    ap = p[active]
    ac = cid[active]
    r_search = max(float(r_contact), float(c_adh))
    half = 0.5 * c_adh
    n = int(active.size)

    # --- bin active nodes on a uniform grid of edge r_search --------------
    mn = ap.min(axis=0)
    cell = cp.floor((ap - mn) / r_search).astype(cp.int64)
    dims = cell.max(axis=0) + 1
    nx, ny, nz = (int(dims[0]), int(dims[1]), int(dims[2]))
    cid_grid = (cell[:, 0] * ny + cell[:, 1]) * nz + cell[:, 2]

    order = cp.argsort(cid_grid)
    cid_sorted = cid_grid[order]
    dims_arr = cp.asarray([nx, ny, nz], dtype=cp.int64)

    # --- candidate pairs from the 27 neighbour-cell offsets ---------------
    i_parts, j_parts = [], []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                off = cp.asarray([dx, dy, dz], dtype=cp.int64)
                ncoord = cell + off
                valid = cp.logical_and(
                    (ncoord >= 0).all(axis=1), (ncoord < dims_arr).all(axis=1))
                ncid = (ncoord[:, 0] * ny + ncoord[:, 1]) * nz + ncoord[:, 2]
                ncid = cp.where(valid, ncid, 0)
                st = cp.searchsorted(cid_sorted, ncid, side="left")
                en = cp.searchsorted(cid_sorted, ncid, side="right")
                st = cp.where(valid, st, 0)
                cnt = cp.where(valid, en - st, 0)
                total = int(cnt.sum())
                if total == 0:
                    continue
                cum = cp.cumsum(cnt)
                k = cp.arange(total, dtype=cp.int64)
                src = cp.searchsorted(cum, k, side="right")
                pair_i = src
                within = k - (cum[src] - cnt[src])
                j_sortedpos = st[src] + within
                pair_j = order[j_sortedpos]
                i_parts.append(pair_i)
                j_parts.append(pair_j)

    if not i_parts:
        return F
    pi = cp.concatenate(i_parts)
    pj = cp.concatenate(j_parts)

    # --- filter: i<j dedup, different cells, within search radius ---------
    keep = pi < pj
    pi, pj = pi[keep], pj[keep]
    keep = ac[pi] != ac[pj]
    pi, pj = pi[keep], pj[keep]
    # r_vec points FROM node pj TO node pi (matches kernels_cpu: ap[ii]-ap[jj]).
    d_vec = ap[pi] - ap[pj]
    r = cp.sqrt(cp.sum(d_vec * d_vec, axis=1))
    keep = cp.logical_and(r < r_search, r > 1e-18)
    if int(keep.sum()) == 0:
        return F
    pi, pj, d_vec, r = pi[keep], pj[keep], d_vec[keep], r[keep]

    # --- per-pair tent law, identical to the reference --------------------
    r_safe = cp.where(r > 1e-18, r, 1e-18)
    rhat = d_vec / r_safe[:, None]
    ov = r_contact - r
    rep_m = ov > 0.0
    tent = cp.where(r >= half, c_adh - r, r)
    if cad_mult is not None:
        cm = cp.asarray(cad_mult)
        mult = cp.sqrt(cm[ac[pi]] * cm[ac[pj]])
    else:
        mult = 1.0
    adh_m = cp.logical_and(~rep_m, r < c_adh) if adh_strength > 0.0 \
        else cp.zeros_like(rep_m)
    fmag = cp.where(
        rep_m,
        rep_strength * patch_area * ov,
        cp.where(adh_m, -(adh_strength * patch_area * tent * mult), 0.0),
    )
    fmag = cp.clip(fmag, -force_cap, force_cap)
    fvec = fmag[:, None] * rhat
    cupyx.scatter_add(F, active[pi], fvec)                   # + on pi
    cupyx.scatter_add(F, active[pj], -fvec)                  # − on pj (Newton-3)
    return F
