"""Differential interfacial tension (foam / DAH) for the DCM cell-cell contact — the missing
faceting drive (PI 2026-06-30).

DIAGNOSIS (workflow wf_2088b933 + the N=48 regime/ω̃ tests): the cells stayed rounded "bag of
marbles" at ANY adhesion strength (contact saturated ~0.19 even at ω̃≈2.4 = 10× lit adhesion),
because the existing adhesion is a *perpendicular* node-face traction that only closes the gap — it
never makes the contact SPREAD — and the cortical surface tension was UNIFORM (same γ on free and
contact faces), so there was no energetic drive to convert free surface into contact area.

THE FIX (this module): the foam / Differential-Adhesion-Hypothesis interfacial tension. The shared
interface between two adhered cells has tension γ₁₂ = γ₁ + γ₂ − w_adh (two cortices minus the
adhesion energy density). Represented per-face (each of the two apposed faces carries half the
interface), a cell-cell **contact** face gets a REDUCED tension

    γ_contact = max(0, γ_free − w_adh/2)        (free faces keep γ_free)

so the system lowers its energy by GROWING contact area at the expense of free area → the contact
spreads (Young-Dupré) → cells deform into a space-filling faceted foam. This is the standard
deformable-cell / wet-foam formulation (Steinberg DAH; Manning PNAS 2010; SimuCell3D Runser 2024;
the Douezan spreading coefficient S = w_cs − 2γ the codebase already references).

w_adh is the MEASURED MCF7 cell-cell adhesion energy density (w_cs = 2.85e-3 J/m²); γ_free is the
physiological cortical tension. **No outcome tuning** — the wetting regime is whatever S = w_cs − 2γ
yields at the lit values. Contact detection reuses the hash-grid-over-face-centroids pattern
(``face_contact_count_kernel``); area gradients reuse ``_tri_area_grads``.
"""

from __future__ import annotations

import numpy as np
import warp as wp

wp.init()

from ffn_sim.dcm.dcm_neighbor_warp import _tri_area_grads


@wp.kernel
def differential_surface_tension_kernel(
    grid: wp.uint64,                         # hash grid over FACE centroids (built each step)
    fcent: wp.array(dtype=wp.vec3),          # (M,) f32 face centroids (query points)
    pos: wp.array(dtype=wp.vec3d),           # (N,) node positions
    faces: wp.array(dtype=wp.int32, ndim=2), # (M,3)
    fcell: wp.array(dtype=wp.int32),         # (M,) owner cell of each face
    radius: wp.float32,                      # contact range (≈ c_adh + face reach)
    gamma_free: wp.float64,                  # cortical tension on media-exposed faces
    w_adh_half: wp.float64,                  # w_adh/2 — tension reduction on a contact face
    force: wp.array(dtype=wp.vec3d),         # (N,) out (atomic scatter)
):
    """Per face: F = −γ_eff·∂A/∂r, with γ_eff REDUCED on cell-cell contact faces (foam/DAH).

    A face is a cell-cell contact iff another cell's face centroid lies within ``radius``. Contact
    faces use γ_eff = max(0, γ_free − w_adh/2); free faces use γ_free. The differential
    (γ_free − γ_contact = w_adh/2) is the spreading drive that flattens cells into a faceted tissue.
    """
    f = wp.tid()
    c1 = fcell[f]
    if c1 < wp.int32(0):
        return
    # cell-cell contact detection (apposed to a DIFFERENT cell's face within radius)
    found = wp.int32(0)
    q = wp.hash_grid_query(grid, fcent[f], radius)
    fj = wp.int32(0)
    while wp.hash_grid_query_next(q, fj):
        if fcell[fj] != c1 and fcell[fj] >= wp.int32(0):
            found = wp.int32(1)
    g = gamma_free
    if found == wp.int32(1):
        g = gamma_free - w_adh_half
        if g < wp.float64(0.0):
            g = wp.float64(0.0)                # complete-wetting floor (no negative/runaway tension)
    ia = faces[f, 0]
    ib = faces[f, 1]
    ic = faces[f, 2]
    ga, gb, gc = _tri_area_grads(pos[ia], pos[ib], pos[ic])
    wp.atomic_add(force, ia, ga * (-g))
    wp.atomic_add(force, ib, gb * (-g))
    wp.atomic_add(force, ic, gc * (-g))


def douezan_spreading(w_cs: float, gamma: float) -> float:
    """Douezan spreading coefficient S = w_cs − 2γ. S>0 ⇒ cells wet/spread on each other (faceting);
    S<0 ⇒ non-wetting (rounded). The regime is SET by the lit values, not tuned."""
    return float(w_cs - 2.0 * gamma)
