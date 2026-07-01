"""Aggregate-level (Foty-Steinberg) liquid-drop surface tension for the DCM spheroid — the MISSING
global densification driver (PI 8h-goal 2026-07-02).

DIAGNOSIS (SPHEROID_STAGED_FORCES_2026-07-01 §2e): dynamic loose→compact compaction does NOT emerge from any
junction lever (reach + contraction + differential γ) — those are all PAIRWISE or PER-CELL. Reach+contraction
give Rg −0.27% but porosity RISES (no global densification). The honest conclusion: compaction needs an
AGGREGATE-level surface tension (Foty-Steinberg tissue surface tension σ, several mN/m) that treats the whole
spheroid as a liquid drop and minimises its OUTER (aggregate-medium) surface — pulling the boundary inward,
eliminating voids — which no per-cell force does.

THE MECHANISM (this module): the tissue behaves as a liquid drop of surface tension σ. The drop's Laplace
pressure ΔP_agg = 2σ/R_agg acts INWARD on the aggregate envelope (the media-exposed / FREE faces). We apply it
as a radial-inward force toward the aggregate centroid on every free face, weighted by its area:

    F_node = (ΔP_agg · A_face / 3) · (−d̂)        d̂ = (face_centroid − c_agg)/|·|   (free faces only)

Free faces = the aggregate-medium interface (a face with NO other-cell face within the contact radius — reuses
the same hash-grid-over-face-centroids detection as the differential-tension kernel). Contact (internal) faces
are skipped (they are handled by adhesion / differential tension / turgor). This is GENUINELY aggregate-level:
c_agg and R_agg are global (the whole drop), so a LOOSE aggregate is pulled together (which per-cell γ, being
local, cannot do). σ is the lit Foty-Steinberg tissue surface tension (swept as a controlled variable over the
~1–20 mN/m band; not tuned to a compaction target). Combined with the incompressible cells (turgor) + non-
penetration (IPC), the envelope contraction becomes void-elimination at constant cell volume = compaction.
"""

from __future__ import annotations

import numpy as np
import warp as wp

wp.init()


@wp.kernel
def aggregate_laplace_kernel(
    grid: wp.uint64,                          # hash grid over FACE centroids (built each step)
    fcent: wp.array(dtype=wp.vec3),           # (M,) f32 face centroids (query points)
    pos: wp.array(dtype=wp.vec3d),            # (N,) node positions
    faces: wp.array(dtype=wp.int32, ndim=2),  # (M,3)
    fcell: wp.array(dtype=wp.int32),          # (M,) owner cell of each face
    radius: wp.float32,                       # contact range (≈ c_adh + face reach) for free/contact test
    cagg: wp.vec3d,                           # aggregate centroid (global)
    dP_agg: wp.float64,                       # aggregate Laplace pressure 2σ/R_agg [Pa]
    force: wp.array(dtype=wp.vec3d),          # (N,) out (atomic scatter)
):
    """Inward aggregate Laplace pressure on the FREE (media-exposed) faces → the liquid-drop that
    contracts + densifies the spheroid. Free faces only; toward the global aggregate centroid."""
    f = wp.tid()
    c1 = fcell[f]
    if c1 < wp.int32(0):
        return
    # free-face detection: a face is FREE iff NO other cell's face centroid is within radius
    found = wp.int32(0)
    q = wp.hash_grid_query(grid, fcent[f], radius)
    fj = wp.int32(0)
    while wp.hash_grid_query_next(q, fj):
        if fcell[fj] != c1 and fcell[fj] >= wp.int32(0):
            found = wp.int32(1)
    if found == wp.int32(1):
        return                                # contact (internal) face — not the aggregate envelope
    ia = faces[f, 0]
    ib = faces[f, 1]
    ic = faces[f, 2]
    a = pos[ia]
    b = pos[ib]
    c = pos[ic]
    fc = (a + b + c) / wp.float64(3.0)
    d = fc - cagg
    r = wp.length(d)
    if r < wp.float64(1e-9):
        return
    nin = -d / r                              # inward, toward the aggregate centroid
    area = wp.float64(0.5) * wp.length(wp.cross(b - a, c - a))
    fmag = dP_agg * area / wp.float64(3.0)
    fv = nin * fmag
    wp.atomic_add(force, ia, fv)
    wp.atomic_add(force, ib, fv)
    wp.atomic_add(force, ic, fv)


def aggregate_centroid_radius(pos_np: np.ndarray, fcell_np: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    """Global aggregate centroid + effective radius R_agg (mean node distance from the centroid) [µm].

    Uses ALL live nodes (the whole drop). R_agg is the Laplace-pressure length scale ΔP=2σ/R_agg.
    """
    c = pos_np.mean(axis=0)
    R = float(np.linalg.norm(pos_np - c, axis=1).mean())
    return c.astype(np.float64), max(R, 1e-9)
