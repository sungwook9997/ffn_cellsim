"""SimuCell3D node-vs-face penalty contact — CPU reference (Phase 1).

Replaces the node-NODE tent contact (``dcm_contact.py`` / ``kernels_*.tent_contact``)
with SimuCell3D's node-vs-FACE signed-distance penalty (CONTACT_MODEL_INDEX=0,
``contact_node_face_via_spring.cpp``). A node of cell c1 is tested against the
closest point of a triangular face of cell c2≠c1; the force is measured in the
FACE plane, so two shells cannot slip between each other's nodes (the failure mode
of node-node center-line contact — the H.7 "LJ explosion").

Per-pair law (re-derived in docs/v2_audit/SIMUCELL3D_INTEGRATION_2026-06-11.md §2.2):
    cpa, bary = closest_point_on_triangle(n.pos; v0,v1,v2)
    r_vec = n.pos − cpa            # face→node, |r_vec| = min_d  (UNNORMALIZED)
    sign  = r_vec · n̂_face
    REPULSION  (sign<0, min_d<c_rep):  F = r_vec · (ξ · A_face)     → |F| = ξ·A·d
    ADHESION   (sign>0, min_d<c_adh):  F = r_vec · force_amp
        hardening (min_d ≥ c_adh/2): force_amp = ω·A·(c_adh/min_d − 1) → |F| = ω·A·(c_adh−d)
        softening (min_d <  c_adh/2): force_amp = ω·A                  → |F| = ω·A·d
    Newton-3: face's 3 nodes += F·bary ;  node −= F
The 1/min_d in the hardening ``force_amp`` is cancelled by |r_vec|=min_d, giving a
symmetric bilinear tent peaking at d=c_adh/2 (= paper Eq.2 traction-separation ×A).

This module is the numpy GROUND TRUTH for the cupy/gpu_local port + its sanity gate
(``tests/test_dcm_face_contact.py``: push a node into a static triangle, verify the
implemented force equals −dE/dx of the stored contact energy). SI units.
"""

from __future__ import annotations

import numpy as np


def closest_point_on_triangle(p: np.ndarray, a: np.ndarray, b: np.ndarray,
                              c: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Closest point of ``p`` on triangle (a,b,c) — Ericson 7-region (RTCD §5.1.5).

    Returns (closest_point (3,), barycentric (3,) for (a,b,c)). Exact barycentric
    in every region (vertex / edge / interior).
    """
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = float(ab @ ap)
    d2 = float(ac @ ap)
    if d1 <= 0.0 and d2 <= 0.0:                       # vertex a
        return a.copy(), np.array([1.0, 0.0, 0.0])
    bp = p - b
    d3 = float(ab @ bp)
    d4 = float(ac @ bp)
    if d3 >= 0.0 and d4 <= d3:                         # vertex b
        return b.copy(), np.array([0.0, 1.0, 0.0])
    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:          # edge ab
        v = d1 / (d1 - d3)
        return a + v * ab, np.array([1.0 - v, v, 0.0])
    cp = p - c
    d5 = float(ab @ cp)
    d6 = float(ac @ cp)
    if d6 >= 0.0 and d5 <= d6:                         # vertex c
        return c.copy(), np.array([0.0, 0.0, 1.0])
    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:          # edge ac
        w = d2 / (d2 - d6)
        return a + w * ac, np.array([1.0 - w, 0.0, w])
    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:   # edge bc
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return b + w * (c - b), np.array([0.0, 1.0 - w, w])
    denom = 1.0 / (va + vb + vc)                       # interior
    v = vb * denom
    w = vc * denom
    return a + ab * v + ac * w, np.array([1.0 - v - w, v, w])


def face_contact_pair(n_pos, v0, v1, v2, *, rep_strength, adh_strength,
                      c_rep, c_adh, cad_mult=1.0):
    """One node-vs-face contact. Returns (F_node (3,), bary (3,), face_idx_force F).

    F is the SimuCell3D ``F = r_vec·amp`` (node gets −F; the 3 face nodes get
    +F·bary). Returns (F_node, F_on_face_nodes (3,3), energy) where F_node = −F.
    """
    cpa, bary = closest_point_on_triangle(n_pos, v0, v1, v2)
    r_vec = n_pos - cpa                               # face→node, |·|=min_d
    min_d = float(np.linalg.norm(r_vec))
    fn = np.cross(v1 - v0, v2 - v0)                    # unnormalised face normal
    area = 0.5 * float(np.linalg.norm(fn))
    nrm = float(np.linalg.norm(fn))
    n_hat = fn / nrm if nrm > 0 else np.zeros(3)
    sign = float(r_vec @ n_hat)
    F = np.zeros(3)
    energy = 0.0
    if sign < 0.0 and min_d < c_rep:                  # REPULSION (penetration)
        amp = rep_strength * area
        F = r_vec * amp
        energy = 0.5 * rep_strength * area * min_d ** 2
    elif adh_strength > 0.0 and sign > 0.0 and min_d < c_adh:   # ADHESION
        half = 0.5 * c_adh
        mult = adh_strength * cad_mult
        if min_d >= half:                             # hardening
            force_amp = mult * (c_adh / max(min_d, 1e-30) - 1.0) * area
            energy = 0.5 * area * mult * (0.5 * c_adh ** 2 - (c_adh - min_d) ** 2)
        else:                                         # softening
            force_amp = mult * area
            energy = 0.5 * mult * area * min_d ** 2
        F = r_vec * force_amp
    F_node = -F
    F_face = np.outer(bary, F)                          # (3,3): +F·bary on 3 nodes
    return F_node, F_face, energy


def node_face_contact_forces(pos, cell_of_node, faces, face_cell, *,
                             rep_strength, adh_strength, c_rep, c_adh,
                             cad_mult=None):
    """Brute-force numpy reference: per-node force (N,3) from node-vs-face contact.

    pos (N,3); cell_of_node (N,) int (−1 = dormant); faces (M,3) node-id triplets;
    face_cell (M,) owner cell of each face. Every active node is tested against
    every face whose owner cell differs, within the relevant cutoff. O(N·M) — the
    reference for the cupy USPG port; not for production N.
    """
    pos = np.asarray(pos, float)
    cof = np.asarray(cell_of_node)
    faces = np.asarray(faces)
    fcell = np.asarray(face_cell)
    F = np.zeros_like(pos)
    cutoff = max(c_rep, c_adh)
    for ni in range(pos.shape[0]):
        c1 = int(cof[ni])
        if c1 < 0:
            continue
        npos = pos[ni]
        for fi in range(faces.shape[0]):
            if int(fcell[fi]) == c1:
                continue                               # same cell — skip
            a, b, c = pos[faces[fi]]
            # cheap AABB pre-filter
            if (npos < np.minimum.reduce([a, b, c]) - cutoff).any() or \
               (npos > np.maximum.reduce([a, b, c]) + cutoff).any():
                continue
            cm = 1.0 if cad_mult is None else float(
                np.sqrt(cad_mult[c1] * cad_mult[int(fcell[fi])]))
            F_node, F_face, _ = face_contact_pair(
                npos, a, b, c, rep_strength=rep_strength,
                adh_strength=adh_strength, c_rep=c_rep, c_adh=c_adh, cad_mult=cm)
            F[ni] += F_node
            for k in range(3):
                F[faces[fi, k]] += F_face[k]
    return F


# ---------------------------------------------------------------------------
# Vectorized, array-module-parameterized port (numpy OR cupy) — the production
# kernel. xp=numpy validates bit-parity vs the brute-force reference locally;
# xp=cupy is the gpu_local GPU path. Same closest-point + law; the O(N·M) search
# is replaced by a SYNC-FREE uniform grid over FACE CENTROIDS (nodes query their
# 27-voxel neighbourhood; r_search = max(c_rep,c_adh)+max_edge so every in-range
# face's centroid is co-resident). Candidates with min_d ≥ cutoff get amp=0, so
# the force matches the brute-force reference regardless of the candidate set.
# ---------------------------------------------------------------------------
def _closest_point_vec(xp, p, a, b, c):
    """Vectorized Ericson 7-region closest point. p,a,b,c are (K,3). Returns
    (cpa (K,3), bary (K,3)). Region barycentrics are computed with guarded
    denominators and selected in REVERSE priority so the vertex regions win."""
    ab = b - a; ac = c - a; ap = p - a
    d1 = xp.sum(ab * ap, axis=1); d2 = xp.sum(ac * ap, axis=1)
    bp = p - b; d3 = xp.sum(ab * bp, axis=1); d4 = xp.sum(ac * bp, axis=1)
    cp = p - c; d5 = xp.sum(ab * cp, axis=1); d6 = xp.sum(ac * cp, axis=1)
    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4
    K = p.shape[0]
    z = xp.zeros(K); o = xp.ones(K)

    def sdiv(num, den):
        return num / xp.where(den != 0, den, 1.0)

    v_ab = sdiv(d1, d1 - d3)
    w_ac = sdiv(d2, d2 - d6)
    w_bc = sdiv(d4 - d3, (d4 - d3) + (d5 - d6))
    denom = sdiv(o, va + vb + vc); v_in = vb * denom; w_in = vc * denom
    bary = xp.stack([o - v_in - w_in, v_in, w_in], axis=1)            # interior
    bary = xp.where(((va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0))[:, None],
                    xp.stack([z, o - w_bc, w_bc], axis=1), bary)       # edge bc
    bary = xp.where(((vb <= 0) & (d2 >= 0) & (d6 <= 0))[:, None],
                    xp.stack([o - w_ac, z, w_ac], axis=1), bary)       # edge ac
    bary = xp.where(((d6 >= 0) & (d5 <= d6))[:, None],
                    xp.stack([z, z, o], axis=1), bary)                 # vertex c
    bary = xp.where(((vc <= 0) & (d1 >= 0) & (d3 <= 0))[:, None],
                    xp.stack([o - v_ab, v_ab, z], axis=1), bary)       # edge ab
    bary = xp.where(((d3 >= 0) & (d4 <= d3))[:, None],
                    xp.stack([z, o, z], axis=1), bary)                 # vertex b
    bary = xp.where(((d1 <= 0) & (d2 <= 0))[:, None],
                    xp.stack([o, z, z], axis=1), bary)                 # vertex a
    cpa = a * bary[:, 0:1] + b * bary[:, 1:2] + c * bary[:, 2:3]
    return cpa, bary


def node_face_contact_forces_vec(xp, pos, cell_of_node, faces, face_cell, *,
                                 rep_strength, adh_strength, c_rep, c_adh,
                                 max_edge, cad_mult=None):
    """Per-node contact force (N,3), array module ``xp`` (numpy or cupy)."""
    pos = xp.asarray(pos, dtype=xp.float64)
    cof = xp.asarray(cell_of_node)
    faces = xp.asarray(faces)
    fcell = xp.asarray(face_cell)
    F = xp.zeros_like(pos)
    M = int(faces.shape[0])
    active = xp.flatnonzero(cof >= 0)
    n = int(active.shape[0])
    if n < 1 or M < 1:
        return F
    # face geometry (recomputed: force ∝ A_face on a deforming mesh)
    v0 = pos[faces[:, 0]]; v1 = pos[faces[:, 1]]; v2 = pos[faces[:, 2]]
    fn = xp.cross(v1 - v0, v2 - v0)
    fnorm = xp.sqrt(xp.sum(fn * fn, axis=1))
    area = 0.5 * fnorm
    n_hat = fn / xp.where(fnorm[:, None] > 0, fnorm[:, None], 1.0)
    fc = (v0 + v1 + v2) / 3.0
    ap_ = pos[active]; ac_ = cof[active]
    r_search = max(float(c_rep), float(c_adh)) + float(max_edge)

    # --- SYNC-FREE uniform grid over FACE CENTROIDS, node queries 27 voxels --
    mn = xp.minimum(ap_.min(axis=0), fc.min(axis=0))
    fcl = xp.floor((fc - mn) / r_search).astype(xp.int64)
    ndl = xp.floor((ap_ - mn) / r_search).astype(xp.int64)
    dims = xp.maximum(fcl.max(axis=0), ndl.max(axis=0)) + 1
    ny, nz = dims[1], dims[2]
    fcid = (fcl[:, 0] * ny + fcl[:, 1]) * nz + fcl[:, 2]
    order = xp.argsort(fcid); fcid_sorted = fcid[order]
    offs = xp.asarray([(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1)
                       for dz in (-1, 0, 1)], dtype=xp.int64)
    st_list, cnt_list = [], []
    for oo in range(27):
        nc = ndl + offs[oo]
        valid = xp.logical_and((nc >= 0).all(axis=1), (nc < dims).all(axis=1))
        ncid = (nc[:, 0] * ny + nc[:, 1]) * nz + nc[:, 2]
        ncid = xp.where(valid, ncid, 0)
        st = xp.searchsorted(fcid_sorted, ncid, side="left")
        en = xp.searchsorted(fcid_sorted, ncid, side="right")
        st_list.append(st)
        cnt_list.append(xp.where(valid, en - st, 0))
    st_all = xp.concatenate(st_list); cnt_all = xp.concatenate(cnt_list)
    src_node = xp.tile(xp.arange(n, dtype=xp.int64), 27)
    cum = xp.cumsum(cnt_all)
    total = int(cum[-1])
    if total == 0:
        return F
    kk = xp.arange(total, dtype=xp.int64)
    seg = xp.searchsorted(cum, kk, side="right")
    within = kk - (cum[seg] - cnt_all[seg])
    ni = src_node[seg]                          # node (index into active)
    fj = order[st_all[seg] + within]            # face index
    keep = ac_[ni] != fcell[fj]                 # different cell only
    ni, fj = ni[keep], fj[keep]
    if int(ni.shape[0]) == 0:
        return F

    # --- closest point + branch + tent/repulsion law ----------------------
    cpa, bary = _closest_point_vec(xp, ap_[ni], v0[fj], v1[fj], v2[fj])
    r_vec = ap_[ni] - cpa
    min_d = xp.sqrt(xp.sum(r_vec * r_vec, axis=1))
    sign = xp.sum(r_vec * n_hat[fj], axis=1)
    Af = area[fj]
    half = 0.5 * c_adh
    if cad_mult is not None:
        cm = xp.asarray(cad_mult)
        mult = xp.sqrt(cm[ac_[ni]] * cm[fcell[fj]])
    else:
        mult = 1.0
    rep_m = (sign < 0.0) & (min_d < c_rep)
    adh_m = ((sign > 0.0) & (min_d < c_adh)) if adh_strength > 0.0 \
        else xp.zeros_like(rep_m)
    md_safe = xp.where(min_d > 0, min_d, 1.0)
    amp_adh = xp.where(min_d >= half,
                       adh_strength * mult * (c_adh / md_safe - 1.0) * Af,
                       adh_strength * mult * Af)
    amp = xp.where(rep_m, rep_strength * Af, xp.where(adh_m, amp_adh, 0.0))
    Fp = r_vec * amp[:, None]                    # F = r_vec·amp
    cupyx_scatter(xp, F, active[ni], -Fp)        # node -= F
    for k in range(3):
        cupyx_scatter(xp, F, faces[fj, k], bary[:, k:k + 1] * Fp)   # face += F·bary
    return F


def cupyx_scatter(xp, out, idx, val):
    """xp-agnostic scatter-add (numpy ``add.at`` / cupyx ``scatter_add``)."""
    if xp.__name__ == "numpy":
        xp.add.at(out, idx, val)
    else:
        import cupyx
        cupyx.scatter_add(out, idx, val)
