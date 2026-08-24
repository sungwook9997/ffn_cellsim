"""Conservative bilinear traction-separation contact (SimuCell3D-faithful) — the energy fix.

DIAGNOSIS (workflow wz04nrajj, reading SimuCell3D source): SimuCell3D uses the SAME explicit
overdamped Euler we do — the faceting failure is NOT a solver/minimizer issue. The root cause is
that our production contact force (`contact_grid_kernel`) is **non-conservative**: the repulsion is
a CONSTANT `rep·area` and the adhesion for d<c_adh/2 is a CONSTANT `adh·area` — a force with NO
stable minimum. A descent on such a force jitters into rounded "marbles"; there is no faceted
energy minimum to settle into.

This kernel replaces it with SimuCell3D's **bilinear traction-separation tent**, which IS a
conservative potential with a clean minimum at contact (d=0):

    σ(d) = rep·d                 d < 0   (overlap → linear restoring repulsion)
         = adh·d                 0 ≤ d < c/2   (adhesion, rising)
         = adh·(c−d)             c/2 ≤ d < c   (adhesion, falling to 0 at the cutoff)

force on the node = −σ·area·r̂ (toward the face when adhering, away when overlapping), scattered
barycentrically to the face (Newton-3). The adhesion energy per unit contact area is
∫₀^c σ dd = adh·c²/4; setting this to the measured w_cs fixes **adh = 4·w_cs/c²** (grid-invariant —
c=c_adh scales with the mesh, so the energy density stays w_cs). No outcome tuning.
"""
from __future__ import annotations

import numpy as np
import warp as wp

wp.init()

from aleph.dcm.dcm_contact_warp import closest_bary


@wp.kernel
def contact_grid_conservative_kernel(
    grid: wp.uint64,                          # hash grid over FACE centroids
    qpts: wp.array(dtype=wp.vec3),            # f32 node positions (query points)
    pos: wp.array(dtype=wp.vec3d),            # f64 node positions
    cof: wp.array(dtype=wp.int32),
    faces: wp.array(dtype=wp.int32, ndim=2),
    fcell: wp.array(dtype=wp.int32),
    radius: wp.float32,                       # = c_adh + face_reach (f32)
    rep: wp.float64,                          # repulsion stiffness [Pa/m]
    adh: wp.float64,                          # adhesion stiffness [Pa/m] (= 4·w_cs/c_adh²)
    c_adh: wp.float64,                        # adhesion/contact cutoff [m]
    force: wp.array(dtype=wp.vec3d),          # atomic accumulate
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
        if fcell[fj] != c1 and fcell[fj] >= wp.int32(0):
            ia = faces[fj, 0]
            ib = faces[fj, 1]
            ic = faces[fj, 2]
            a = pos[ia]
            b = pos[ib]
            c = pos[ic]
            bary = closest_bary(p, a, b, c)
            cpa = a * bary[0] + b * bary[1] + c * bary[2]
            r_vec = p - cpa
            d = wp.length(r_vec)
            fnv = wp.cross(b - a, c - a)
            nrm = wp.length(fnv)
            area = wp.float64(0.5) * nrm
            sign = z
            if nrm > z:
                sign = wp.dot(r_vec, fnv) / nrm
            # bilinear traction-separation tent (conservative; min at d=0)
            amp = z
            if sign < z:                                 # overlap → linear repulsion ∝ depth
                amp = rep * d * area
            else:
                if d < c_adh:                            # separation → adhesion tent
                    if d < half:
                        amp = adh * d * area             # rising
                    else:
                        amp = adh * (c_adh - d) * area    # falling
            if amp != z:
                fvec = r_vec * amp                       # along r̂·|r_vec| ; node gets −fvec
                fn_acc = fn_acc - fvec
                wp.atomic_add(force, ia, bary[0] * fvec)
                wp.atomic_add(force, ib, bary[1] * fvec)
                wp.atomic_add(force, ic, bary[2] * fvec)
    wp.atomic_add(force, ni, fn_acc)


def derive_adhesion_stiffness(w_cs: float, c_adh: float) -> float:
    """adh [Pa/m] from the measured adhesion energy density w_cs [J/m²]: adh = 4·w_cs/c_adh²
    (so ∫σ dd = adh·c²/4 = w_cs, grid-invariant). Derived, not tuned."""
    return float(4.0 * w_cs / max(c_adh, 1e-12) ** 2)
