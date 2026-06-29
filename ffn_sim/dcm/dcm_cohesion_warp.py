"""Warp port of the DCM node-node cell-cell cohesion force (Phase B, cohesion).

Ports ``cell.dcm_contact.DcmTentContact`` — the inter-cell node-node contact:
non-penetration REPULSION below the contact radius plus a symmetric bilinear
ADHESION tent (cadherin-scale cohesion) that pulls neighbouring cells together.

Per inter-cell pair (i in cell c1, j in cell c2≠c1), r_vec = pos_i − pos_j, d = |r_vec|:
    REPULSION (d < r_contact):           fmag = +ξ·A·(r_contact − d)     (push apart)
    ADHESION  (r_contact ≤ d < c_adh):   tent = (c_adh − d) if d ≥ c_adh/2 else d
                                         fmag = −ω·A·m·tent              (pull together)
        m = √(cad[c1]·cad[c2])  (cadherin multiplier; 1 if absent)
    fmag clamped to ±force_cap ;  F_i += fmag·r̂ ,  F_j −= fmag·r̂  (Newton-3)

Atomic-free Newton-3: ONE THREAD PER NODE i sums fmag·r̂_ij over all j≠i (different
cell). The symmetric pair is computed independently by thread j with r_vec_ji =
−r_vec_ij and the same fmag, giving exactly −(fmag·r̂_ij) on node j — so writing only
the thread's own row reproduces F_i += fvec / F_j −= fvec with no atomics. Brute-force
O(N²) (the reference's KD-tree cutoff only drops pairs that contribute 0 anyway).
Per-node summation order differs from the reference's ``np.add.at`` → tolerance-class.
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()


@wp.kernel
def dcm_cohesion_kernel(
    pos: wp.array(dtype=wp.vec3d),
    cof: wp.array(dtype=wp.int32),             # (N,) cell_of_node (−1 dormant)
    cad: wp.array(dtype=wp.float64),           # (n_cells,) cadherin mult (or dummy)
    use_cad: wp.int32,
    N: wp.int32,
    r_contact: wp.float64,
    c_adh: wp.float64,
    rep: wp.float64,
    omega: wp.float64,
    A: wp.float64,
    force_cap: wp.float64,
    force: wp.array(dtype=wp.vec3d),           # (N,) out (own-row write, no atomics)
):
    i = wp.tid()
    ci = cof[i]
    if ci < wp.int32(0):
        return
    z = wp.float64(0.0)
    pi = pos[i]
    fx = z
    fy = z
    fz = z
    half = wp.float64(0.5) * c_adh
    for j in range(N):
        cj = cof[j]
        if j != i and cj >= wp.int32(0) and cj != ci:
            pj = pos[j]
            rx = pi[0] - pj[0]
            ry = pi[1] - pj[1]
            rz = pi[2] - pj[2]
            d = wp.sqrt(rx * rx + ry * ry + rz * rz)
            d_safe = d
            if d_safe <= wp.float64(1.0e-18):
                d_safe = wp.float64(1.0e-18)
            fmag = z
            ov = r_contact - d
            if ov > z:
                fmag = rep * A * ov
            else:
                if omega > z and d < c_adh:
                    if d >= half:
                        tent = c_adh - d
                    else:
                        tent = d
                    amag = omega * A * tent
                    if use_cad == wp.int32(1):
                        amag = amag * wp.sqrt(cad[ci] * cad[cj])
                    fmag = -amag
            if fmag > force_cap:
                fmag = force_cap
            if fmag < -force_cap:
                fmag = -force_cap
            rhx = rx / d_safe
            rhy = ry / d_safe
            rhz = rz / d_safe
            fx = fx + fmag * rhx
            fy = fy + fmag * rhy
            fz = fz + fmag * rhz
    force[i] = wp.vec3d(fx, fy, fz)


def run_dcm_cohesion_warp(
    *,
    pos: np.ndarray,
    cell_of_node: np.ndarray,
    r_contact: float,
    c_adh: float,
    rep_strength: float,
    adh_strength: float,
    patch_area: float,
    force_cap: float = 5.0e-8,
    cad_mult: np.ndarray | None = None,
    device: str = "cpu",
) -> dict:
    """DCM node-node cohesion force (N,3) in Warp, matching dcm_contact.DcmTentContact."""
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    cof = np.ascontiguousarray(cell_of_node, dtype=np.int32)
    N = pos.shape[0]
    use_cad = 1 if cad_mult is not None else 0
    cad = (np.ascontiguousarray(cad_mult, dtype=np.float64)
           if cad_mult is not None else np.zeros(1, dtype=np.float64))

    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    cof_d = wp.array(cof, dtype=wp.int32, device=device)
    cad_d = wp.array(cad, dtype=wp.float64, device=device)
    force_d = wp.zeros(N, dtype=wp.vec3d, device=device)

    wp.launch(
        dcm_cohesion_kernel, dim=N,
        inputs=[pos_d, cof_d, cad_d, wp.int32(use_cad), wp.int32(N),
                wp.float64(r_contact), wp.float64(c_adh), wp.float64(rep_strength),
                wp.float64(adh_strength), wp.float64(patch_area),
                wp.float64(force_cap), force_d],
        device=device,
    )
    wp.synchronize_device(device)
    return {"force": force_d.numpy().astype(np.float64)}
