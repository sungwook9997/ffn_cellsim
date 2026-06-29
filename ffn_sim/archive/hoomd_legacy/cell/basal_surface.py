"""Basal SURFACE layer (S1; KU-3.5, PI directive (a) 2026-06-09).

Two-layer basal contractile apparatus (``docs/v2_audit/BASAL_MESH_DESIGN_2026-06-09.md``):

* **S-layer (THIS module)** — a **FLAT 2D triangulated SURFACE MESH** (SimuCell3D-style)
  used as a *positioning + connectivity* substrate: it gives the focal-adhesion
  particles a real 2D surface to be PLACED on (FA is otherwise hard to attach — the
  spatial-disjoint integrin problem), and it carries SOME surface connectivity (the
  triangulation's patch/vertex adjacency, "표면 연락 일부"). It carries **NO in-plane /
  edge force** — a force-bearing surface edge is the 2nd unsanctioned coarse-graining
  and over-constrains single-filament buckling, the Gate-B lever (REJECTED,
  ``H7_CORTEX_AS_MESH_2026-06-07.md`` / commit ``7b19276``). Geometry / positioning /
  contact ONLY (not a HOOMD bond → invisible to the γ estimator).
* **F-layer (``cell/basal_mesh.py``)** — the explicit fine-grained filament network
  built ON this surface, anchored to the FA particles, where the REAL FA mechanics /
  traction are measured ("실제 FA 메카닉스는 그 표면 위 필라멘트 네트워크에서 본다").

Why FLAT (PI 2026-06-09)
------------------------
An ADHERENT / spread cell flattens its VENTRAL surface against the substrate — the
basal contact zone where FA and ventral stress fibers live is approximately a FLAT
disk (the cell body / apical surface stays rounded above). A curved sphere-cap is the
SUSPENDED-cell geometry, wrong for the traction / SF scenario. A flat ventral surface
also makes the lamellipodium (a flat basal protrusion) natural to add later. The
filaments (``cell/basal_mesh.py``, also planar) thus live consistently ON this flat
surface — there is no flat-vs-curved choice once the ventral surface is flattened.

The flat surface is a Delaunay triangulation of disk-sampled points (concentric rings
for even coverage) embedded at ``z = z_basal`` with outward (cell-interior, +z) normal.
Disk area ``π·footprint_radius²`` is resolution-invariant by construction (VG-1).

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. Checks in ``ffn_sim/tests/test_basal_surface.py``.*

1. **Dimensional** — ``footprint_radius`` [m]; positions [m]; areas [m²].
2. **Boundary** — ``footprint_radius ≤ 0`` → ValueError; ``n_fa ≤ 0`` → ValueError;
   ``n_fa`` above the patch count → clamp.
3. **Conservation / topology** — the triangulation is a SINGLE connected component
   (a Delaunay disk is connected); area → ``π·footprint_radius²`` with resolution.
4. **Numerical** — positions finite float64.
5. **Sign / sense (FLAT, in the disk)** — every vertex/FA at ``z = z_basal`` exactly
   (planar); every FA inside the footprint disk (radial ≤ footprint_radius); outward
   normal is ``+z`` (cell-interior; substrate is below at lower z).
6. **Measurement-protocol (resolution-invariance, VG-1)** — triangulated disk area is
   invariant to ring resolution (→ ``π·footprint_radius²``); NO force, not a bond →
   γ-invisible.

References
----------
- ``docs/v2_audit/H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md`` (manifold as
  geometry) + ``H7_CORTEX_AS_MESH_2026-06-07.md`` (force-bearing-surface rejection).
- ``docs/v2_audit/BASAL_MESH_DESIGN_2026-06-09.md`` (the 2-layer plan).
- ``cell/basal_mesh.py`` (the F-layer filament network placed ON this surface).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(slots=True)
class FlatBasalSurface:
    """A flat triangulated ventral surface + FA placement (S-layer geometry).

    Attributes:
        verts: (n_vert, 3) float64 — surface vertices at z = z_basal.
        tris: (n_tri, 3) int64 — triangle vertex indices (Delaunay).
        tri_adj: (n_tri, 3) int64 — edge-neighbour triangle indices (−1 = boundary).
        tri_centroids: (n_tri, 3) float64 — triangle centroids (at z_basal).
        normal: (3,) float64 — outward (cell-interior, +z) unit normal.
        footprint_radius: disk radius [m].
        z_basal: basal plane height [m].
        area: triangulated disk area [m²] (→ π·footprint_radius²).
        fa_positions: (n_fa, 3) float64 — FA-anchor positions ON the surface.
        fa_home_tri: (n_fa,) int64 — the triangle each FA sits on.
    """

    verts: np.ndarray
    tris: np.ndarray
    tri_adj: np.ndarray
    tri_centroids: np.ndarray
    normal: np.ndarray
    footprint_radius: float
    z_basal: float
    area: float
    fa_positions: np.ndarray
    fa_home_tri: np.ndarray

    @property
    def n_tri(self) -> int:
        return int(self.tris.shape[0])

    @property
    def n_vert(self) -> int:
        return int(self.verts.shape[0])


def _disk_ring_points(footprint_radius: float, n_rings: int) -> np.ndarray:
    """Concentric-ring point cloud on a disk (even coverage), shape (P, 2)."""
    if n_rings < 1:
        raise ValueError(f"n_rings must be ≥ 1; got {n_rings}")
    pts = [np.array([[0.0, 0.0]])]
    dr = footprint_radius / n_rings
    for k in range(1, n_rings + 1):
        rk = k * dr
        # ~even arc spacing ≈ dr between ring points.
        n_k = max(6, int(round(2.0 * math.pi * rk / dr)))
        ang = np.linspace(0.0, 2.0 * math.pi, n_k, endpoint=False)
        pts.append(np.stack([rk * np.cos(ang), rk * np.sin(ang)], axis=1))
    return np.concatenate(pts, axis=0)


def _triangulate_disk(points2d: np.ndarray):
    """Delaunay-triangulate 2D disk points → (tris, tri_neighbors)."""
    from scipy.spatial import Delaunay
    d = Delaunay(points2d)
    return (
        np.asarray(d.simplices, dtype=np.int64),
        np.asarray(d.neighbors, dtype=np.int64),  # −1 for boundary edges
    )


def build_flat_basal_surface(
    *,
    footprint_radius: float,
    z_basal: float,
    n_rings: int = 10,
    n_fa: int = 60,
    rng: np.random.Generator | None = None,
) -> FlatBasalSurface:
    """Build the FLAT ventral surface mesh + FA particles placed ON it.

    Args:
        footprint_radius: ventral contact-disk radius [m].
        z_basal: basal plane height [m] (the substrate-contact z).
        n_rings: concentric-ring resolution of the disk triangulation.
        n_fa: number of FA-anchor particles to place on the surface.
        rng: RNG for FA placement.

    Returns:
        A :class:`FlatBasalSurface`.
    """
    if not (math.isfinite(footprint_radius) and footprint_radius > 0.0):
        raise ValueError(f"footprint_radius must be finite > 0; got {footprint_radius!r}")
    if not math.isfinite(z_basal):
        raise ValueError(f"z_basal must be finite; got {z_basal!r}")
    if n_fa <= 0:
        raise ValueError(f"n_fa must be > 0; got {n_fa}")
    if rng is None:
        rng = np.random.default_rng(0)

    pts2d = _disk_ring_points(footprint_radius, n_rings)
    tris, tri_adj = _triangulate_disk(pts2d)
    verts = np.column_stack([pts2d, np.full(pts2d.shape[0], z_basal)])

    v0 = verts[tris[:, 0]]
    v1 = verts[tris[:, 1]]
    v2 = verts[tris[:, 2]]
    tri_centroids = (v0 + v1 + v2) / 3.0
    # flat-triangle areas (xy only since planar).
    cross_z = (v1[:, 0] - v0[:, 0]) * (v2[:, 1] - v0[:, 1]) - \
              (v1[:, 1] - v0[:, 1]) * (v2[:, 0] - v0[:, 0])
    area = float(0.5 * np.abs(cross_z).sum())
    normal = np.array([0.0, 0.0, 1.0])   # outward = into the cell (+z)

    # FA placed at triangle centroids (subset, no replacement), ON the surface.
    m = min(int(n_fa), tris.shape[0])
    chosen = rng.permutation(tris.shape[0])[:m].astype(np.int64)
    fa_positions = tri_centroids[chosen].copy()
    fa_home_tri = chosen

    return FlatBasalSurface(
        verts=verts,
        tris=tris,
        tri_adj=tri_adj,
        tri_centroids=tri_centroids,
        normal=normal,
        footprint_radius=float(footprint_radius),
        z_basal=float(z_basal),
        area=area,
        fa_positions=fa_positions,
        fa_home_tri=fa_home_tri,
    )


def _is_connected(tris: np.ndarray, tri_adj: np.ndarray) -> bool:
    """BFS over edge-adjacency → single connected component?"""
    n = tris.shape[0]
    if n == 0:
        return False
    seen = {0}
    stack = [0]
    while stack:
        t = stack.pop()
        for nb in tri_adj[t]:
            nb = int(nb)
            if nb >= 0 and nb not in seen:
                seen.add(nb)
                stack.append(nb)
    return len(seen) == n


def basal_surface_report(
    surf: FlatBasalSurface, *, planar_tol: float = 1e-12,
) -> dict[str, Any]:
    """S1 build-time gate metrics (FLAT, in the disk, connected).

    Controls:
      * planar          — every vertex/FA at z = z_basal (within ``planar_tol``).
      * within_disk     — every FA radial (x,y) ≤ footprint_radius (+ε).
      * connected       — the triangulation is one edge-connected component.
      * area_matches     — triangulated area ≈ π·footprint_radius² (disk coverage).
    """
    z_dev = float(np.max(np.abs(surf.verts[:, 2] - surf.z_basal))) if surf.n_vert else 0.0
    fa = np.asarray(surf.fa_positions, dtype=np.float64)
    planar = bool(z_dev <= planar_tol
                  and (fa.size == 0 or np.all(np.abs(fa[:, 2] - surf.z_basal) <= planar_tol)))
    radial = np.linalg.norm(fa[:, :2], axis=1) if fa.size else np.zeros(0)
    within = bool(np.all(radial <= surf.footprint_radius + 1e-12)) if fa.size else False
    connected = _is_connected(surf.tris, surf.tri_adj)
    disk_area = math.pi * surf.footprint_radius ** 2
    area_frac = surf.area / disk_area if disk_area > 0 else 0.0
    area_ok = 0.9 <= area_frac <= 1.0 + 1e-9  # Delaunay hull underfills the circle a bit
    controls = {
        "n_vert": surf.n_vert,
        "n_tri": surf.n_tri,
        "n_fa_placed": int(fa.shape[0]),
        "planar": planar,
        "max_z_dev": z_dev,
        "within_disk": within,
        "connected": connected,
        "area_m2": surf.area,
        "area_fraction_of_disk": area_frac,
        "area_matches": area_ok,
    }
    ok = planar and within and connected and area_ok and fa.shape[0] > 0
    return {"verdict": "PASS" if ok else "REVIEW", "controls": controls}


def disk_area_resolution_invariance(
    *, footprint_radius: float, z_basal: float = 0.0,
    ring_levels: tuple[int, ...] = (6, 10, 16),
) -> dict[str, Any]:
    """VG-1 master gate: triangulated disk area is invariant to ring resolution.

    Builds the flat surface at several ring resolutions and checks the area
    converges to the analytic disk area ``π·footprint_radius²`` (the Delaunay hull
    underfills the circle by a deficit that SHRINKS with resolution → the area is
    resolution-invariant, not drifting with physics).
    """
    analytic = math.pi * footprint_radius ** 2
    areas = {}
    for n in ring_levels:
        surf = build_flat_basal_surface(
            footprint_radius=footprint_radius, z_basal=z_basal, n_rings=n, n_fa=1,
        )
        areas[n] = surf.area
    fracs = np.array([a / analytic for a in areas.values()])
    spread = float(fracs.max() - fracs.min())
    return {
        "analytic_disk_area": analytic,
        "area_by_rings": {int(k): float(v) for k, v in areas.items()},
        "fraction_by_rings": {int(k): float(v / analytic) for k, v in areas.items()},
        "rel_spread": spread,
        "invariant": bool(spread < 0.10),  # converges to the disk; hull deficit only
    }
