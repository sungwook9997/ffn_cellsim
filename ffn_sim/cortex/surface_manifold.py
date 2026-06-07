"""Deformable-capable triangulated SPHERE manifold — GEOMETRY / broad-phase ONLY.

A standalone, simulator-free triangulated cell-surface manifold scaffold for the
H.7 surface-manifold design
(``docs/v2_audit/H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX_2026-06-07.md`` §1-§7;
``docs/v2_audit/H1_H5_H7_MANIFOLD_CONTACT_ARCHITECTURE_2026-06-07.md``).

Per that design the manifold is a **coordinate / locality / soft-confinement
substrate ONLY** — it carries NO mechanics. This module is the *geometry half*:
it supplies the cell-surface shape, per-triangle local tangent/normal frames,
a bead→triangle map (``nearest_patch``), and a geodesic k-ring neighbourhood
(``patch_kring``) for broad-phase candidate search. It contains **no HOOMD, no
forces, no tension, no edge/area springs, and no γ** — only ``numpy`` and
``scipy.spatial`` (for the bead→triangle nearest-centroid query).

The triangulation is an **icosphere** (recursive 4-into-1 subdivision of a unit
icosahedron, projected to the sphere) at radius ``R``. Subdivision level ``s``
gives ``20·4ˢ`` triangles and ``10·4ˢ + 2`` vertices — the canonical
resolution-invariance reference grid {320, 1280, 5120} is levels {2, 3, 4}.

The class stores vertices in a mutable ``verts`` array and exposes ``set_verts``
so the manifold can be re-fit / advected to a deforming bead cloud (the
"slaved to beads" invariant); the connectivity (``tris``, ``tri_adj``) is fixed
under deformation. That makes it *deformable-capable* without ever owning an
elastic DOF.

Sanity Gate (per CLAUDE.md Sanity Gate Protocol; checked in
``manifold_search_benchmark.py`` and ``test_surface_manifold.py``):

* **Dimensional.** ``R``, ``verts``, ``tri_centroids`` are lengths [m];
  ``tri_normals``/``tri_e1``/``tri_e2`` are dimensionless unit vectors;
  ``mean_edge_length`` is a length [m]. ``nearest_patch`` maps a length-3 point
  to a dimensionless triangle index. ``patch_kring`` is purely combinatorial
  (no units). No quantity here enters a force budget.
* **Boundary.** A degenerate request (``subdivisions < 0``) raises. ``k < 0``
  raises. ``patch_kring(t, 0)`` returns ``{t}`` (a triangle is its own 0-ring).
  An empty ``points`` array maps to an empty index array. The base icosahedron
  (level 0) has exactly 20 triangles / 12 vertices / 30 edges
  (Euler: V−E+F = 12−30+20 = 2 — a closed genus-0 surface).
* **Conservation / topology.** Every triangle has **exactly 3** edge-neighbours
  (a closed triangulated sphere is 3-regular in its face-adjacency graph);
  ``tri_adj`` is symmetric (t in adj[u] ⟺ u in adj[t]); the face-adjacency
  graph is a single connected component. ``Σ`` triangle areas → ``4πR²`` as
  ``s → ∞`` (flat-triangle area underestimates the curved sphere; the relative
  deficit shrinks monotonically with resolution — asserted in the benchmark as a
  geometry check, never used as a physical observable).
* **Normal orientation (sign-sense).** Every face normal points **outward**
  (``n̂ · ĉ > 0`` where ``ĉ`` is the unit centroid direction) — enforced at
  construction by orienting the vertex winding. ``(n̂, e1, e2)`` is a
  right-handed orthonormal frame (``e1 × e2 ≈ n̂``).
* **Numerical.** Pure float64 geometry; no time integration, no CFL. The
  bead→triangle query is an exact nearest-centroid ``cKDTree`` lookup. Frames are
  re-orthonormalised so round-off cannot accumulate a non-unit basis.
* **Measurement-protocol consistency.** ``mean_edge_length`` is the *unique*
  geometric scale the broad-phase k-ring count is derived from
  (``k = ceil(reach / mean_edge_length)`` + a safety ring) — it is a geometry
  property of the mesh, NOT a tunable. This is the hook that keeps candidate sets
  resolution-invariant (the manifold supplies geometry; the *physics reach* is
  supplied by the caller).
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

__all__ = [
    "SurfaceManifold",
    "icosphere",
    "n_tri_for_subdivisions",
    "n_vert_for_subdivisions",
]


# ---------------------------------------------------------------------------
# Icosphere construction (pure geometry)
# ---------------------------------------------------------------------------
def n_tri_for_subdivisions(subdivisions: int) -> int:
    """Triangle count of an icosphere at the given subdivision level.

    Args:
        subdivisions: Recursive subdivision level ``s ≥ 0``.

    Returns:
        ``20 · 4ˢ`` — the number of triangular faces.
    """
    if subdivisions < 0:
        raise ValueError(f"subdivisions must be ≥ 0; got {subdivisions}")
    return 20 * (4 ** subdivisions)


def n_vert_for_subdivisions(subdivisions: int) -> int:
    """Vertex count of an icosphere at the given subdivision level.

    Args:
        subdivisions: Recursive subdivision level ``s ≥ 0``.

    Returns:
        ``10 · 4ˢ + 2`` — the number of vertices (Euler-consistent with
        ``20·4ˢ`` faces on a closed genus-0 surface).
    """
    if subdivisions < 0:
        raise ValueError(f"subdivisions must be ≥ 0; got {subdivisions}")
    return 10 * (4 ** subdivisions) + 2


def _base_icosahedron() -> tuple[np.ndarray, np.ndarray]:
    """Return (verts (12,3) unit, faces (20,3) int) of a regular icosahedron.

    The 12 vertices are the cyclic permutations of ``(0, ±1, ±φ)`` with
    ``φ`` the golden ratio, normalised to the unit sphere.
    """
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    raw = np.array(
        [
            (-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
            (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
            (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1),
        ],
        dtype=np.float64,
    )
    verts = raw / np.linalg.norm(raw, axis=1, keepdims=True)
    faces = np.array(
        [
            (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
            (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
            (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
            (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
        ],
        dtype=np.int64,
    )
    return verts, faces


def icosphere(
    subdivisions: int, radius: float
) -> tuple[np.ndarray, np.ndarray]:
    """Build a triangulated icosphere by recursive 4-into-1 subdivision.

    Each subdivision step splits every triangle into four by inserting a vertex
    at each edge midpoint and re-projecting all vertices to the sphere. The
    result is a closed, genus-0, (nearly) geodesic triangulation.

    Args:
        subdivisions: Recursive subdivision level ``s ≥ 0`` (``20·4ˢ`` faces).
        radius: Sphere radius ``R`` [m] (``> 0``).

    Returns:
        ``(verts, faces)`` where ``verts`` has shape ``(n_vert, 3)`` [m] and
        ``faces`` has shape ``(n_tri, 3)`` int64 (vertex indices, outward
        winding).

    Raises:
        ValueError: If ``subdivisions < 0`` or ``radius`` is not finite > 0.
    """
    if subdivisions < 0:
        raise ValueError(f"subdivisions must be ≥ 0; got {subdivisions}")
    if not (np.isfinite(radius) and radius > 0.0):
        raise ValueError(f"radius must be finite > 0; got {radius!r}")

    verts, faces = _base_icosahedron()
    verts = list(map(tuple, verts))  # list of unit (x,y,z) tuples
    midpoint_cache: dict[tuple[int, int], int] = {}

    def _midpoint(i: int, j: int) -> int:
        key = (i, j) if i < j else (j, i)
        cached = midpoint_cache.get(key)
        if cached is not None:
            return cached
        a = np.asarray(verts[i])
        b = np.asarray(verts[j])
        m = a + b
        m /= np.linalg.norm(m)  # re-project to the unit sphere
        idx = len(verts)
        verts.append(tuple(m))
        midpoint_cache[key] = idx
        return idx

    for _ in range(subdivisions):
        new_faces: list[tuple[int, int, int]] = []
        for a, b, c in faces:
            ab = _midpoint(a, b)
            bc = _midpoint(b, c)
            ca = _midpoint(c, a)
            new_faces.extend(
                [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
            )
        faces = np.asarray(new_faces, dtype=np.int64)

    v = np.asarray(verts, dtype=np.float64) * float(radius)
    return v, faces


# ---------------------------------------------------------------------------
# The manifold
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class SurfaceManifold:
    """Triangulated cell-surface manifold — geometry / broad-phase substrate.

    A finite triangulated sphere at radius ``R`` that supplies, per the H.7
    surface-manifold design, ONLY: (a) the cell-surface shape (``verts``,
    ``tris``); (b) per-triangle local frames ``(n̂, e1, e2)`` and centroids;
    (c) a bead→triangle map (:meth:`nearest_patch`); and (d) geodesic k-ring
    neighbourhoods (:meth:`patch_kring`) for broad-phase candidate search.

    It carries **no mechanics**: no edge/area/bending springs, no tension, no
    force, no γ — and must never be given any. The vertices are mutable
    (:meth:`set_verts`) so the manifold can be slaved to a deforming bead cloud;
    the *connectivity* (``tris``, ``tri_adj``) is fixed under deformation.

    Attributes:
        R: Construction sphere radius [m].
        verts: Vertex positions, shape ``(n_vert, 3)`` [m] (mutable).
        tris: Triangle vertex indices, shape ``(n_tri, 3)`` int64 (outward
            winding).
        tri_adj: Edge-neighbour adjacency, shape ``(n_tri, 3)`` int64 — the
            three triangles sharing an edge with each triangle (3-regular on a
            closed sphere). Used for geodesic BFS / k-ring.
        tri_centroids: Triangle centroids, shape ``(n_tri, 3)`` [m].
        tri_normals: Outward unit normals, shape ``(n_tri, 3)``.
        tri_e1: First tangent unit vector per triangle, shape ``(n_tri, 3)``.
        tri_e2: Second tangent unit vector per triangle, shape ``(n_tri, 3)``
            (``e1 × e2 = n̂``, right-handed).
        tri_areas: Flat-triangle areas, shape ``(n_tri,)`` [m²].
        mean_edge_length: Mean triangle-edge length [m] (a reporting scale).
        min_centroid_step: Minimum centroid-to-edge-neighbour-centroid distance
            [m] — the WORST-CASE (smallest) geodesic advance per k-ring hop, used
            to DERIVE a conservative k-ring count.
        max_circumradius: Maximum triangle circumradius [m] — the farthest a bead
            can sit from its home-patch centroid; doubled into the k-ring
            coverage bound (centre bead + candidate bead both at a far corner).
    """

    R: float
    verts: np.ndarray
    tris: np.ndarray
    tri_adj: np.ndarray
    tri_centroids: np.ndarray = field(default=None)  # type: ignore[assignment]
    tri_normals: np.ndarray = field(default=None)  # type: ignore[assignment]
    tri_e1: np.ndarray = field(default=None)  # type: ignore[assignment]
    tri_e2: np.ndarray = field(default=None)  # type: ignore[assignment]
    tri_areas: np.ndarray = field(default=None)  # type: ignore[assignment]
    mean_edge_length: float = 0.0
    min_centroid_step: float = 0.0
    max_circumradius: float = 0.0
    _centroid_tree: cKDTree = field(default=None, repr=False)  # type: ignore[assignment]

    # -- construction -------------------------------------------------------
    @classmethod
    def icosphere(cls, subdivisions: int, radius: float) -> "SurfaceManifold":
        """Build an icosphere manifold at ``radius`` and the given subdivision.

        Args:
            subdivisions: Recursive subdivision level ``s ≥ 0`` (``20·4ˢ`` faces;
                {2, 3, 4} are the canonical {320, 1280, 5120} resolution grid).
            radius: Sphere radius ``R`` [m] (``> 0``).

        Returns:
            A fully-populated :class:`SurfaceManifold`.
        """
        verts, tris = icosphere(subdivisions, radius)
        tri_adj = _build_tri_adjacency(tris)
        m = cls(R=float(radius), verts=verts, tris=tris, tri_adj=tri_adj)
        m._recompute_geometry()
        return m

    # -- deformation (slaved-to-beads invariant) ----------------------------
    def set_verts(self, verts: np.ndarray) -> None:
        """Replace vertex positions (e.g. re-fit / advect to a deformed cloud).

        The connectivity (``tris``, ``tri_adj``) is unchanged; only per-triangle
        geometry (centroids, frames, areas, edge length, the centroid tree) is
        recomputed. This is the "manifold slaved to the bead cloud" update — it
        owns no elastic DOF and never relaxes toward its own minimum.

        Args:
            verts: New vertex positions, shape ``(n_vert, 3)`` [m].
        """
        verts = np.asarray(verts, dtype=np.float64)
        if verts.shape != self.verts.shape:
            raise ValueError(
                f"verts shape {verts.shape} != current {self.verts.shape}; "
                "set_verts changes positions only, not connectivity"
            )
        self.verts = verts
        self._recompute_geometry()

    def _recompute_geometry(self) -> None:
        """(Re)compute centroids, outward frames, areas, edge length, KD-tree."""
        v0 = self.verts[self.tris[:, 0]]
        v1 = self.verts[self.tris[:, 1]]
        v2 = self.verts[self.tris[:, 2]]
        centroids = (v0 + v1 + v2) / 3.0

        # Face normal from the winding; orient OUTWARD via the centroid ray.
        cross = np.cross(v1 - v0, v2 - v0)
        cross_norm = np.linalg.norm(cross, axis=1, keepdims=True)
        normals = cross / np.maximum(cross_norm, 1.0e-300)
        cdir = centroids / np.maximum(
            np.linalg.norm(centroids, axis=1, keepdims=True), 1.0e-300
        )
        flip = np.sum(normals * cdir, axis=1) < 0.0
        normals[flip] *= -1.0

        # Right-handed tangent frame: e1 along edge v0→v1 projected ⟂ n̂, e2 = n̂×e1.
        e1 = (v1 - v0) - np.sum((v1 - v0) * normals, axis=1, keepdims=True) * normals
        e1_norm = np.linalg.norm(e1, axis=1, keepdims=True)
        # Degenerate fallback: pick any vector ⟂ n̂ (only if v0,v1 ~ coincident).
        bad = (e1_norm[:, 0] < 1.0e-30)
        if bad.any():
            e1[bad] = _arbitrary_perp(normals[bad])
            e1_norm = np.linalg.norm(e1, axis=1, keepdims=True)
        e1 = e1 / np.maximum(e1_norm, 1.0e-300)
        e2 = np.cross(normals, e1)
        e2 = e2 / np.maximum(np.linalg.norm(e2, axis=1, keepdims=True), 1.0e-300)

        # 0.5·|cross| is the flat-triangle area.
        areas = 0.5 * cross_norm[:, 0]

        # Mean edge length over the three edges of every triangle.
        e_ab = np.linalg.norm(v1 - v0, axis=1)
        e_bc = np.linalg.norm(v2 - v1, axis=1)
        e_ca = np.linalg.norm(v0 - v2, axis=1)
        mean_edge = float(np.mean(np.concatenate([e_ab, e_bc, e_ca])))

        # --- k-ring coverage geometry (worst-case, for a conservative bound) ---
        # Per-hop geodesic advance is the centroid→edge-neighbour-centroid step;
        # use its MINIMUM (smallest advance = most hops needed → conservative).
        nb = self.tri_adj
        valid = nb >= 0
        steps = np.linalg.norm(
            centroids[:, None, :] - centroids[np.where(valid, nb, 0)], axis=2
        )
        steps = steps[valid]
        min_step = float(steps.min()) if steps.size else mean_edge
        # Farthest a bead can sit from its home-patch centroid = max circumradius
        # (max vertex-to-centroid distance over all triangles).
        circ = np.maximum.reduce([
            np.linalg.norm(v0 - centroids, axis=1),
            np.linalg.norm(v1 - centroids, axis=1),
            np.linalg.norm(v2 - centroids, axis=1),
        ])
        max_circ = float(circ.max())

        self.tri_centroids = centroids
        self.tri_normals = normals
        self.tri_e1 = e1
        self.tri_e2 = e2
        self.tri_areas = areas
        self.mean_edge_length = mean_edge
        self.min_centroid_step = min_step
        self.max_circumradius = max_circ
        self._centroid_tree = cKDTree(centroids)

    # -- introspection ------------------------------------------------------
    @property
    def n_tri(self) -> int:
        """Number of triangular faces."""
        return int(self.tris.shape[0])

    @property
    def n_vert(self) -> int:
        """Number of vertices."""
        return int(self.verts.shape[0])

    def total_area(self) -> float:
        """Sum of flat-triangle areas [m²] (→ ``4πR²`` as resolution → ∞)."""
        return float(np.sum(self.tri_areas))

    def frame(self, tri: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the local ``(n̂, e1, e2)`` orthonormal frame of a triangle.

        Args:
            tri: Triangle index.

        Returns:
            ``(n_hat, e1, e2)`` each a unit 3-vector; ``e1 × e2 = n_hat``.
        """
        return (
            self.tri_normals[tri].copy(),
            self.tri_e1[tri].copy(),
            self.tri_e2[tri].copy(),
        )

    # -- bead → triangle map ------------------------------------------------
    def nearest_patch(self, points: np.ndarray) -> np.ndarray:
        """Map each point (bead) to its nearest-centroid triangle (patch id).

        This is the broad-phase "home patch" of each bead. It is a nearest-
        centroid lookup (exact, via ``cKDTree``) — adequate at the mesoscale
        where the cortex shell band (``h_cortex = 200 nm``) is far thinner than a
        patch, so a bead's home patch is unambiguous.

        Args:
            points: Query positions, shape ``(n, 3)`` [m].

        Returns:
            ``(n,)`` int64 triangle indices. An empty input gives an empty array.
        """
        points = np.asarray(points, dtype=np.float64)
        if points.size == 0:
            return np.empty((0,), dtype=np.int64)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError(f"points must be (n, 3); got {points.shape}")
        _, idx = self._centroid_tree.query(points, k=1)
        return np.asarray(idx, dtype=np.int64)

    # -- geodesic k-ring ----------------------------------------------------
    def patch_kring(self, tri: int, k: int) -> np.ndarray:
        """Geodesic k-ring neighbourhood of a triangle (BFS over ``tri_adj``).

        Returns every triangle within ``k`` edge-adjacency hops of ``tri``
        (inclusive of ``tri`` itself at hop 0). This is the patch-local
        candidate-search neighbourhood: beads homed in ``patch_kring(t, k)`` are
        the broad-phase candidate pool for a bead on triangle ``t``.

        ``k`` must be chosen (by the caller) so the ring covers the physical
        reach ball — see :meth:`kring_for_reach`. Resolution-invariance of the
        candidate set follows from that coverage condition, NOT from any tuning
        here.

        Args:
            tri: Centre triangle index.
            k: Number of adjacency hops (``≥ 0``). ``k = 0`` → ``{tri}``.

        Returns:
            Sorted int64 array of triangle indices within ``k`` hops.
        """
        if k < 0:
            raise ValueError(f"k must be ≥ 0; got {k}")
        if not (0 <= tri < self.n_tri):
            raise ValueError(f"tri {tri} out of range [0, {self.n_tri})")
        visited = {int(tri)}
        frontier = deque([(int(tri), 0)])
        while frontier:
            t, hop = frontier.popleft()
            if hop >= k:
                continue
            for nb in self.tri_adj[t]:
                nb = int(nb)
                if nb < 0 or nb in visited:
                    continue
                visited.add(nb)
                frontier.append((nb, hop + 1))
        return np.array(sorted(visited), dtype=np.int64)

    def kring_for_reach(self, reach: float, safety_rings: int = 1) -> int:
        """DERIVE the k-ring hop count that covers a physical reach ball.

        The coverage condition (companion design doc: *patch geodesic radius
        ≥ physical reach*) is solved from the mesh's OWN geometry, bead-to-bead:

            ``k · min_centroid_step − 2 · max_circumradius ≥ reach``

        i.e. ``k = ceil((reach + 2·max_circumradius) / min_centroid_step)`` plus
        ``safety_rings``. The terms are:

        * ``min_centroid_step`` — the smallest centroid→neighbour-centroid step,
          a **lower bound** on the geodesic advance per face-hop, so ``k`` hops
          guarantee at least ``k·min_centroid_step`` of geodesic radius (using
          the minimum is conservative — never undercovers).
        * ``2·max_circumradius`` — a home-patch bead and a candidate bead can each
          sit a full circumradius from their patch centroids, so the bead-to-bead
          reach must be paid on both ends.

        This is a **DERIVED, grid-aware** quantity: ``reach`` is supplied by the
        caller (the physics: ``max_bind_dist`` or the ``√(A/n)`` bridge reach),
        the geometry (``min_centroid_step``, ``max_circumradius``) is a property
        of the mesh. It is **never tuned to pass a gate**. As resolution rises,
        ``min_centroid_step`` falls and ``k`` grows automatically (the "k scales
        with resolution" requirement) — this is exactly what keeps the patch
        candidate set resolution-invariant (a coarser mesh covers the same
        physical ball with fewer, larger hops; a finer mesh with more, smaller
        hops; the covered *physical radius* is the same).

        Args:
            reach: Physical candidate-search radius [m] (e.g. ``max_bind_dist``
                or the ``√(A/n)`` bridge reach).
            safety_rings: Extra rings added for ball-coverage robustness against
                chord-vs-arc and discrete-step slack (``≥ 0``; default 1).

        Returns:
            The k-ring hop count (an int ``≥ 1``).
        """
        if not (np.isfinite(reach) and reach > 0.0):
            raise ValueError(f"reach must be finite > 0; got {reach!r}")
        if safety_rings < 0:
            raise ValueError(f"safety_rings must be ≥ 0; got {safety_rings}")
        if self.min_centroid_step <= 0.0:
            raise ValueError("coverage geometry not computed")
        k = math.ceil(
            (reach + 2.0 * self.max_circumradius) / self.min_centroid_step
        )
        return int(k) + int(safety_rings)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _build_tri_adjacency(tris: np.ndarray) -> np.ndarray:
    """Edge-neighbour adjacency of a closed triangle mesh, shape ``(n_tri, 3)``.

    Two faces are neighbours iff they share an (undirected) edge. On a closed
    genus-0 surface every edge is shared by exactly two faces, so each face has
    exactly three neighbours; the result is dense and 3-regular.

    Args:
        tris: Triangle vertex indices, shape ``(n_tri, 3)`` int.

    Returns:
        ``(n_tri, 3)`` int64 neighbour indices (entries are ``-1`` only for an
        open boundary edge, which does not occur on the closed icosphere).

    Raises:
        ValueError: If any edge is shared by more than two faces (non-manifold).
    """
    n_tri = int(tris.shape[0])
    edge_to_faces: dict[tuple[int, int], list[int]] = {}
    for f in range(n_tri):
        a, b, c = int(tris[f, 0]), int(tris[f, 1]), int(tris[f, 2])
        for i, j in ((a, b), (b, c), (c, a)):
            key = (i, j) if i < j else (j, i)
            edge_to_faces.setdefault(key, []).append(f)

    adj = np.full((n_tri, 3), -1, dtype=np.int64)
    slot = np.zeros(n_tri, dtype=np.int64)
    for faces in edge_to_faces.values():
        if len(faces) > 2:
            raise ValueError(
                f"edge shared by {len(faces)} faces (non-manifold mesh)"
            )
        if len(faces) == 2:
            f0, f1 = faces
            adj[f0, slot[f0]] = f1
            slot[f0] += 1
            adj[f1, slot[f1]] = f0
            slot[f1] += 1
    return adj


def _arbitrary_perp(normals: np.ndarray) -> np.ndarray:
    """Return a unit vector perpendicular to each normal (degenerate fallback).

    Uses the minimum-|component| world axis crossed with the normal, mirroring
    ``cortex.py:_tangent_plane_basis`` so a degenerate triangle still yields a
    valid frame.
    """
    n = normals.shape[0]
    abs_n = np.abs(normals)
    min_axis = np.argmin(abs_n, axis=1)
    ref = np.zeros_like(normals)
    ref[np.arange(n), min_axis] = 1.0
    perp = np.cross(ref, normals)
    return perp / np.maximum(
        np.linalg.norm(perp, axis=1, keepdims=True), 1.0e-300
    )
