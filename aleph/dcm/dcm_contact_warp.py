"""Warp port of the DCM SimuCell3D node-vs-face contact (Phase B, DCM forces).

Mirrors ``cell/dcm_face_contact.py`` (the committed numpy GROUND TRUTH) and the
``gpu_opt/nodeface_rawkernel.py`` CUDA kernel: a node of cell c1 is tested against
the closest point of each triangular face of a DIFFERENT cell c2≠c1; the penalty is
measured in the face plane (signed-distance), so two shells cannot slip between each
other's nodes.

Per-pair law (SimuCell3D, docs/v2_audit/_historical/SIMUCELL3D_INTEGRATION_2026-06-11.md §2.2):
    cpa, bary = closest_point_on_triangle(p; a,b,c)   # Ericson 7-region
    r_vec = p − cpa,  min_d = |r_vec|,  sign = r_vec·n̂_face,  A = ½|n_face|
    REPULSION (sign<0, min_d<c_rep):   amp = ξ·A
    ADHESION  (sign>0, min_d<c_adh):   amp = ω·m·(c_adh/min_d − 1)·A   if min_d ≥ c_adh/2
                                       amp = ω·m·A                      if min_d <  c_adh/2
        m = √(cad[c1]·cad[c2])  (cadherin multiplier; 1 if absent)
    F = r_vec·amp ;  node −= F ;  face's 3 nodes += bary·F   (Newton-3)

One thread per node, looping ALL faces (all-pairs O(N·M) — exact vs the brute-force
``node_face_contact_forces`` reference; the production kernel uses the voxel grid in
``nodeface_rawkernel.py``, an acceleration structure that changes speed, not forces:
any face the grid drops has min_d ≥ cutoff → amp=0). Node's own force and the
face-node scatter are both ``wp.atomic_add`` (a node is also a face vertex of other
threads' faces, so plain writes would race). SI units.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


@wp.func
def closest_bary(p: wp.vec3d, a: wp.vec3d, b: wp.vec3d, c: wp.vec3d) -> wp.vec3d:
    """Barycentric (a,b,c) of the closest point of p on triangle (a,b,c) — Ericson
    7-region (RTCD §5.1.5), branch-for-branch identical to the scalar reference
    ``dcm_face_contact.closest_point_on_triangle`` (same priority: vertex a → b →
    edge ab → vertex c → edge ac → edge bc → interior)."""
    z = wp.float64(0.0)
    o = wp.float64(1.0)
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = wp.dot(ab, ap)
    d2 = wp.dot(ac, ap)
    if d1 <= z and d2 <= z:                       # vertex a
        return wp.vec3d(o, z, z)
    bp = p - b
    d3 = wp.dot(ab, bp)
    d4 = wp.dot(ac, bp)
    if d3 >= z and d4 <= d3:                       # vertex b
        return wp.vec3d(z, o, z)
    vc = d1 * d4 - d3 * d2
    if vc <= z and d1 >= z and d3 <= z:            # edge ab
        v = d1 / (d1 - d3)
        return wp.vec3d(o - v, v, z)
    cp = p - c
    d5 = wp.dot(ab, cp)
    d6 = wp.dot(ac, cp)
    if d6 >= z and d5 <= d6:                       # vertex c
        return wp.vec3d(z, z, o)
    vb = d5 * d2 - d1 * d6
    if vb <= z and d2 >= z and d6 <= z:            # edge ac
        w = d2 / (d2 - d6)
        return wp.vec3d(o - w, z, w)
    va = d3 * d6 - d5 * d4
    if va <= z and (d4 - d3) >= z and (d5 - d6) >= z:   # edge bc
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return wp.vec3d(z, o - w, w)
    denom = o / (va + vb + vc)                      # interior
    v = vb * denom
    w = vc * denom
    return wp.vec3d(o - v - w, v, w)


@wp.kernel
def node_face_contact_kernel(
    pos: wp.array(dtype=wp.vec3d),              # (N,) ro node positions
    cof: wp.array(dtype=wp.int32),              # (N,) cell_of_node (−1 = dormant)
    faces: wp.array(dtype=wp.int32, ndim=2),    # (M, 3) node-id triplets
    fcell: wp.array(dtype=wp.int32),            # (M,) owner cell of each face
    cad: wp.array(dtype=wp.float64),            # (n_cell,) cadherin mult (or dummy)
    use_cad: wp.int32,                          # 1 → apply √(cad·cad), 0 → 1.0
    m: wp.int32,                                # face count (inner-loop bound)
    rep: wp.float64,                            # ξ repulsion strength
    adh: wp.float64,                            # ω adhesion strength
    c_rep: wp.float64,
    c_adh: wp.float64,
    force: wp.array(dtype=wp.vec3d),            # (N,) out (zeroed; atomic accumulate)
):
    ni = wp.tid()
    c1 = cof[ni]
    if c1 < wp.int32(0):                        # dormant node — no contact
        return
    p = pos[ni]
    z = wp.float64(0.0)
    half = wp.float64(0.5) * c_adh
    fn_acc = wp.vec3d(z, z, z)                  # node's own force, register-accumulated
    for fj in range(m):
        if fcell[fj] != c1:                     # different cell only
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
            fn = wp.cross(b - a, c - a)
            nrm = wp.length(fn)
            area = wp.float64(0.5) * nrm
            sign = z
            if nrm > z:
                sign = wp.dot(r_vec, fn) / nrm
            amp = z
            if sign < z and min_d < c_rep:                       # REPULSION
                amp = rep * area
            else:
                if adh > z and sign > z and min_d < c_adh:       # ADHESION
                    mult = adh
                    if use_cad == wp.int32(1):
                        mult = adh * wp.sqrt(cad[c1] * cad[fcell[fj]])
                    if min_d >= half:                            # hardening
                        md = min_d
                        if md <= z:
                            md = wp.float64(1.0e-30)
                        amp = mult * (c_adh / md - wp.float64(1.0)) * area
                    else:                                        # softening
                        amp = mult * area
            if amp != z:
                fvec = r_vec * amp
                fn_acc = fn_acc - fvec                           # node −= F
                wp.atomic_add(force, ia, bary[0] * fvec)         # face += bary·F
                wp.atomic_add(force, ib, bary[1] * fvec)
                wp.atomic_add(force, ic, bary[2] * fvec)
    wp.atomic_add(force, ni, fn_acc)


def run_node_face_contact_warp(
    *,
    pos: np.ndarray,            # (N, 3)
    cell_of_node: np.ndarray,   # (N,) int, −1 = dormant
    faces: np.ndarray,          # (M, 3) int
    face_cell: np.ndarray,      # (M,) int
    rep_strength: float,
    adh_strength: float,
    c_rep: float,
    c_adh: float,
    cad_mult: np.ndarray | None = None,   # (n_cell,) or None
    device: str = "cpu",
) -> dict:
    """Node-vs-face contact force (N,3) in Warp (one thread/node, all-pairs over
    faces), matching ``dcm_face_contact.node_face_contact_forces``."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    cof = np.ascontiguousarray(cell_of_node, dtype=np.int32)
    faces = np.ascontiguousarray(faces, dtype=np.int32)
    fcell = np.ascontiguousarray(face_cell, dtype=np.int32)
    N = pos.shape[0]
    M = faces.shape[0]
    use_cad = 1 if cad_mult is not None else 0
    cad = (np.ascontiguousarray(cad_mult, dtype=np.float64)
           if cad_mult is not None else np.zeros(1, dtype=np.float64))

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof, dtype=wp.int32, device=device)
    faces_d = wp.array(faces, dtype=wp.int32, device=device)
    fcell_d = wp.array(fcell, dtype=wp.int32, device=device)
    cad_d = wp.array(cad, dtype=wp.float64, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)

    wp.launch(
        node_face_contact_kernel, dim=N,
        inputs=[pos_d, cof_d, faces_d, fcell_d, cad_d, wp.int32(use_cad),
                wp.int32(M), wp.float64(rep_strength), wp.float64(adh_strength),
                wp.float64(c_rep), wp.float64(c_adh), force_d],
        device=device,
    )
    wp.synchronize_device(device)
    return {"force": force_d.numpy().astype(np.float64)}
