"""The plasma membrane: a closed triangulated surface built on the device, by kernels.

WHAT THIS IS, AND WHY IT IS NOT ``surface.build_surface``.  ``aleph.world.surface`` builds a surface
from vertices and faces a caller already holds, on the host, in NumPy.  That is the right shape for a
few hundred vertices and the wrong one for a hundred and sixty thousand: the host icosphere dedups
midpoints through a Python dict, one edge at a time, and at subdivision 7 that is ~1.5 million
dictionary operations before a single byte reaches the card.  PI 2026-08-20 settled the general form —
construction is a Warp CUDA kernel, never a host build uploaded — and this module is that form for the
closed-surface populations.

THE HASH IS WHAT HAD TO GO, AND AN EDGE LIST IS WHAT REPLACES IT.  Loop subdivision needs each edge's
midpoint to be created exactly once and found again from all the faces that share it; a dict keyed on
the vertex pair is the obvious way and it is inherently serial.  Carry an explicit EDGE LIST beside the
faces and the lookup disappears into arithmetic::

    level n:    V vertices,  E edges,  F faces
    level n+1:  V + E,       2E + 3F,  4F

    the midpoint of edge e IS vertex  V + e          — no search, no dedup, no collision
    the two children of edge e ARE    2e, 2e+1
    the three interior edges of face f ARE  2E + 3f + {0,1,2}

Every index a level needs is a closed-form function of the level below, so both halves of a
subdivision are one parallel launch over edges and one over faces, with no atomics and no ordering
dependence.  The cost is one extra ``(E, 2)`` array, which is smaller than the faces it accompanies.

WHY THE BASE ICOSAHEDRON IS STILL BUILT ON THE HOST.  It has twelve vertices.  Uploading a twelve-vertex
seed is not "building the population on the host"; it is a constant, and it is taken from
``surface.icosphere(0)`` rather than retyped so the two modules cannot drift into two different
icosahedra.  Everything from level 1 up — 100% of the vertices this module exists to make — is written
by kernels into device memory that no host array ever mirrors.

EQUIVALENCE TO THE HOST BUILDER IS TESTED, NOT ASSERTED.  ``surface.icosphere`` normalises once at the
END rather than at every level, which is a different vertex distribution from per-level projection, and
this module reproduces the host one: the base vertices all share a norm, so scaling them to the unit
sphere commutes with taking midpoints, and the two agree exactly.  ``_demo`` checks that at subdivision
2 against ``surface.icosphere`` — as a multiset, because the edge-list scheme numbers vertices in a
different order than the dict does and vertex ORDER is not part of what a mesh is.

THE SUBDIVISION LEVEL IS DERIVED, NEVER TYPED.  A level is not a physiological quantity; a mesh SPACING
is.  :func:`subdivisions_for_mesh` returns the coarsest level whose mean triangle edge falls at or below
a declared spacing, so the caller declares the length it can source and the level falls out.  At
``R = 7.5 µm`` and the sourced 50-100 nm cortical/membrane-skeleton mesh this returns 7 (70.6 nm edge;
level 6 is 141 nm, outside the band), which is the level ``WORLD_PORT_PLAN`` and PI 2026-07-22 both
name — but it is reached from the spacing, so a different radius or a different spacing gives a
different answer rather than the same typed 7.

WHAT IS DELIBERATELY ABSENT.  No bending rigidity, no surface tension, no area reservoir, no pressure,
no ERM tether.  PHASE 1 is geometry and nothing else.  ``mobility`` is left at its allocated zero for
the same reason: a mobility is a drag law's parameter, and a builder that filled it in would be
deciding physics under the name of construction.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — vertex positions and mesh spacing [µm]; ``area0`` [µm²]; enclosed volume [µm³];
    subdivision level and every index are dimensionless.
  * boundary — a negative subdivision level, a non-positive radius or spacing, and a spacing coarser
    than the base icosahedron are refused; ``radius_um`` and ``mesh_um`` have NO defaults and the build
    refuses without them, the ratified pattern from ``surface.build_surface``; an edge used by other
    than exactly two faces is counted on the device and raises rather than producing a partial hinge
    list.
  * conservation/invariant — Euler's ``V - E + F = 2`` is checked at every level from the closed-form
    counts before a byte is allocated, and the built claims tile the arena's live prefix
    (``assert_partitioned``).  Hinge count is E = 3F/2, not F.
  * CFL/precision — no integration.  float64 for every position, area and volume.  The random-free
    construction means the build is bit-reproducible for a given level and radius.
  * sign sense — winding decides the sign of every face's pressure traction, so the discrete
    divergence-theorem volume is computed from the mesh's own triangles and a non-positive value
    raises.  The child-face orders below preserve the parent's winding; ``_demo`` checks the volume is
    positive and that inverting the base makes it negative.
  * measurement protocol — the volume and total area are summed on the HOST from a per-face device
    array rather than by device atomics: an atomic float64 reduction is order-dependent and would make
    the rest volume irreproducible in its last bits.  Both readbacks are build-time, between accepted
    steps by construction, since nothing has stepped yet.

engine units: length µm, area µm², volume µm³.  Runtime: NVIDIA Warp on CUDA.  This module is
CPU-importable — it allocates and launches nothing until a device is given — but it will not build on
one: ``build_closed_surface`` requires a CUDA arena and raises otherwise.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import warp as wp

from aleph.world.arena import Claim, Kind, WorldArena
from aleph.world.surface import icosphere

__all__ = [
    "ClosedSurface",
    "build_closed_surface",
    "build_membrane",
    "icosphere_counts",
    "mean_edge_um",
    "subdivisions_for_mesh",
]


# ── counts, spacings and the level that follows from them ───────────────────────────────────────
def icosphere_counts(subdivisions: int) -> tuple[int, int, int]:
    """``(V, E, F)`` of an icosphere at ``subdivisions``, in closed form.

    The recursion ``V' = V + E``, ``E' = 2E + 3F``, ``F' = 4F`` from the base ``(12, 30, 20)``, which
    is the same recursion the kernels below implement.  Having it available as arithmetic is what lets
    every array be sized, and every Euler check run, before a device allocation exists.

    Raises:
        ValueError: on a negative level.
    """
    if subdivisions < 0:
        raise ValueError(f"subdivisions must be nonnegative; got {subdivisions}")
    v, e, f = 12, 30, 20
    for _ in range(subdivisions):
        v, e, f = v + e, 2 * e + 3 * f, 4 * f
        if v - e + f != 2:
            raise AssertionError(f"Euler characteristic broke at V={v} E={e} F={f}")
    return v, e, f


def mean_edge_um(radius_um: float, subdivisions: int) -> float:
    """Mean triangle edge length [µm] of an icosphere of ``radius_um`` at ``subdivisions``.

    From area rather than from geometry: the sphere's ``4πR²`` is shared by ``F`` near-equilateral
    triangles, so ``a = sqrt(4A_face / sqrt(3))``.  This is the length a mesh SPACING is compared
    against, and it is a mean — an icosphere's edges vary by a few percent about it.
    """
    _, _, f = icosphere_counts(subdivisions)
    area_face = 4.0 * math.pi * radius_um * radius_um / f
    return math.sqrt(4.0 * area_face / math.sqrt(3.0))


def subdivisions_for_mesh(radius_um: float, mesh_um: float, *, max_level: int = 10) -> int:
    """The COARSEST subdivision whose mean edge is at or below ``mesh_um``, for ``radius_um``.

    This is the whole reason no subdivision level is typed anywhere in this package.  A level is a
    discretisation index and carries no physiology; a mesh spacing is a sourced length.  Declare the
    length, get the level.

    Args:
        radius_um: sphere radius [µm].
        mesh_um: the largest acceptable mean triangle edge [µm] — a declared, sourced length.
        max_level: refuse rather than run away if the spacing is unreachably fine.  Level 10 is
            ~10.5 M vertices, already past a whole cell's node budget, so reaching it means the
            spacing is wrong rather than that the ceiling is low.

    Raises:
        ValueError: on a non-positive radius or spacing, or a spacing needing more than ``max_level``.
    """
    if not (math.isfinite(radius_um) and radius_um > 0.0):
        raise ValueError(f"radius_um must be finite and positive; got {radius_um!r}")
    if not (math.isfinite(mesh_um) and mesh_um > 0.0):
        raise ValueError(f"mesh_um must be finite and positive; got {mesh_um!r}")
    for n in range(max_level + 1):
        if mean_edge_um(radius_um, n) <= mesh_um:
            return n
    raise ValueError(
        f"a mesh spacing of {mesh_um} um at R={radius_um} um needs more than subdivision {max_level} "
        f"({icosphere_counts(max_level)[0]:,} vertices, already past a whole cell's budget). Either the "
        "spacing is not the one that was sourced, or this population needs a representation that is "
        "not a uniform icosphere."
    )


# ── the built population ────────────────────────────────────────────────────────────────────────
@dataclass(frozen=True, slots=True)
class ClosedSurface:
    """One closed-surface population: ONE claim per primitive kind, and device-resident topology.

    Three claims for the whole population, not three per element — the arena's own words are that "a
    population is a NAMED, CONTIGUOUS, HALF-OPEN ID RANGE", and a surface was already built that way.
    The strand builder was not, which is the change :mod:`aleph.world.build.cortex` carries.

    Attributes:
        population: the name the ranges are attributed to.
        nodes / faces / hinges: the three claims.
        subdivisions: the level built — DERIVED from a mesh spacing, never typed.
        radius_um / centre_um: the sphere this is.
        mesh_um_requested: the spacing the level was derived from [µm].
        mesh_um_realised: the mean triangle edge actually built [µm].  **The one a law must read.**
        mesh_provenance: a free-text label for where ``mesh_um_requested`` came from.  Carried into the
            record so a PI-GAP travels with the number instead of beside it.
        face_idx: ``(F, 3)`` device int32, outward-wound, GLOBAL arena node ids.
        area0_um2: ``(F,)`` device float64 rest areas [µm²].
        hinge_idx: ``(H, 4)`` device int32 dihedrals ``(a, b, c, d)`` with ``b-c`` the shared edge,
            GLOBAL ids.  ``a < d`` always, so the list does not depend on kernel scheduling order.
        volume0_um3: enclosed rest volume [µm³], discrete, by the divergence theorem.
        area0_total_um2: total rest area [µm²].
        topology_bytes: device bytes this population's own arrays hold, EXCLUDING the arena's shared
            node arrays.  Reported so bytes-per-node can be split into shared and owned.
    """

    population: str
    nodes: Claim
    faces: Claim
    hinges: Claim
    subdivisions: int
    radius_um: float
    centre_um: tuple[float, float, float]
    mesh_um_requested: float
    mesh_um_realised: float
    mesh_provenance: str
    face_idx: wp.array = field(repr=False)
    area0_um2: wp.array = field(repr=False)
    hinge_idx: wp.array = field(repr=False)
    volume0_um3: float
    area0_total_um2: float
    topology_bytes: int

    @property
    def n_vertices(self) -> int:
        """Number of vertices."""
        return self.nodes.count

    def record(self) -> dict[str, object]:
        """Host-side summary for a run record.  Touches no device array."""
        return {
            "population": self.population,
            "kind": "closed_surface",
            "subdivisions": self.subdivisions,
            "n_vertices": self.n_vertices,
            "n_faces": self.faces.count,
            "n_hinges": self.hinges.count,
            "radius_um": self.radius_um,
            "centre_um": list(self.centre_um),
            "mesh_um_requested": self.mesh_um_requested,
            "mesh_um_realised": self.mesh_um_realised,
            "mesh_provenance": self.mesh_provenance,
            "volume0_um3": self.volume0_um3,
            "area0_total_um2": self.area0_total_um2,
            "topology_bytes": self.topology_bytes,
            "claims": {
                "node": [self.nodes.lo, self.nodes.hi],
                "face": [self.faces.lo, self.faces.hi],
                "angle4": [self.hinges.lo, self.hinges.hi],
            },
        }


# ── kernels ─────────────────────────────────────────────────────────────────────────────────────
@wp.func
def _edge_child(edge: wp.array(dtype=wp.int32, ndim=2), e: wp.int32, endpoint: wp.int32) -> wp.int32:
    """Index of the child of edge ``e`` that still touches ``endpoint``.

    Edge ``e = (u, v)`` splits into ``2e = (u, m)`` and ``2e+1 = (m, v)``, so the child touching an
    endpoint is decided by which end that endpoint is.  This is the entire replacement for the
    midpoint hash: a face finds its children's edges by asking which side of a parent edge it is on.
    """
    if edge[e, 0] == endpoint:
        return 2 * e
    return 2 * e + 1


@wp.kernel
def _copy_vertices(src: wp.array(dtype=wp.vec3d), dst: wp.array(dtype=wp.vec3d)):
    """Carry a level's vertices into the next level's buffer; midpoints append after them."""
    i = wp.tid()
    dst[i] = src[i]


@wp.kernel
def _split_edges(
    vert: wp.array(dtype=wp.vec3d),
    edge: wp.array(dtype=wp.int32, ndim=2),
    n_vert: wp.int32,
    vert_out: wp.array(dtype=wp.vec3d),
    edge_out: wp.array(dtype=wp.int32, ndim=2),
):
    """One midpoint and two child edges per parent edge.

    The midpoint is the plain arithmetic mean and is NOT projected onto the sphere here: the host
    builder projects once at the end, and matching it is what makes the two comparable.  Projection is
    scale-commuting across a level because every base vertex shares a norm, so the final projection in
    :func:`_place_vertices` recovers exactly the host's vertex set.
    """
    e = wp.tid()
    u = edge[e, 0]
    v = edge[e, 1]
    m = n_vert + e
    vert_out[m] = (vert[u] + vert[v]) * wp.float64(0.5)
    edge_out[2 * e + 0, 0] = u
    edge_out[2 * e + 0, 1] = m
    edge_out[2 * e + 1, 0] = m
    edge_out[2 * e + 1, 1] = v


@wp.kernel
def _split_faces(
    face_v: wp.array(dtype=wp.int32, ndim=2),
    face_e: wp.array(dtype=wp.int32, ndim=2),
    edge: wp.array(dtype=wp.int32, ndim=2),
    n_vert: wp.int32,
    n_edge: wp.int32,
    edge_out: wp.array(dtype=wp.int32, ndim=2),
    face_v_out: wp.array(dtype=wp.int32, ndim=2),
    face_e_out: wp.array(dtype=wp.int32, ndim=2),
):
    """Four children and three interior edges per parent face, winding preserved.

    ``face_e[f, k]`` is the edge joining ``face_v[f, k]`` to ``face_v[f, (k+1) % 3]``; the children
    keep that convention, which is what lets the next level call :func:`_edge_child` at all.  The four
    child orders ``(a, m_ab, m_ca) (b, m_bc, m_ab) (c, m_ca, m_bc) (m_ab, m_bc, m_ca)`` are the host
    builder's, and each traverses its parent's sense — an inverted child would flip the sign of its own
    pressure traction, which is the failure this convention exists to avoid.
    """
    f = wp.tid()
    a = face_v[f, 0]
    b = face_v[f, 1]
    c = face_v[f, 2]
    e_ab = face_e[f, 0]
    e_bc = face_e[f, 1]
    e_ca = face_e[f, 2]
    m_ab = n_vert + e_ab
    m_bc = n_vert + e_bc
    m_ca = n_vert + e_ca

    i0 = 2 * n_edge + 3 * f
    i1 = i0 + 1
    i2 = i0 + 2
    edge_out[i0, 0] = m_ab
    edge_out[i0, 1] = m_bc
    edge_out[i1, 0] = m_bc
    edge_out[i1, 1] = m_ca
    edge_out[i2, 0] = m_ca
    edge_out[i2, 1] = m_ab

    f0 = 4 * f
    face_v_out[f0 + 0, 0] = a
    face_v_out[f0 + 0, 1] = m_ab
    face_v_out[f0 + 0, 2] = m_ca
    face_e_out[f0 + 0, 0] = _edge_child(edge, e_ab, a)
    face_e_out[f0 + 0, 1] = i2
    face_e_out[f0 + 0, 2] = _edge_child(edge, e_ca, a)

    face_v_out[f0 + 1, 0] = b
    face_v_out[f0 + 1, 1] = m_bc
    face_v_out[f0 + 1, 2] = m_ab
    face_e_out[f0 + 1, 0] = _edge_child(edge, e_bc, b)
    face_e_out[f0 + 1, 1] = i0
    face_e_out[f0 + 1, 2] = _edge_child(edge, e_ab, b)

    face_v_out[f0 + 2, 0] = c
    face_v_out[f0 + 2, 1] = m_ca
    face_v_out[f0 + 2, 2] = m_bc
    face_e_out[f0 + 2, 0] = _edge_child(edge, e_ca, c)
    face_e_out[f0 + 2, 1] = i1
    face_e_out[f0 + 2, 2] = _edge_child(edge, e_bc, c)

    face_v_out[f0 + 3, 0] = m_ab
    face_v_out[f0 + 3, 1] = m_bc
    face_v_out[f0 + 3, 2] = m_ca
    face_e_out[f0 + 3, 0] = i0
    face_e_out[f0 + 3, 1] = i1
    face_e_out[f0 + 3, 2] = i2


@wp.kernel
def _place_vertices(
    vert: wp.array(dtype=wp.vec3d),
    radius: wp.float64,
    centre: wp.vec3d,
    node_lo: wp.int32,
    range_id: wp.int32,
    pos: wp.array(dtype=wp.vec3d),
    strand_id: wp.array(dtype=wp.int32),
    rng_id: wp.array(dtype=wp.int32),
):
    """Project onto the sphere, scale, translate, and write into the arena's claimed node range.

    ``strand_id`` is left at ``-1``: a surface vertex belongs to no filament, and writing a fabricated
    strand index would make a later "which strand is this node on" query answer confidently and wrongly.
    """
    v = wp.tid()
    pos[node_lo + v] = wp.normalize(vert[v]) * radius + centre
    strand_id[node_lo + v] = -1
    rng_id[node_lo + v] = range_id


@wp.kernel
def _face_metrics(
    pos: wp.array(dtype=wp.vec3d),
    face_v: wp.array(dtype=wp.int32, ndim=2),
    node_lo: wp.int32,
    face_idx: wp.array(dtype=wp.int32, ndim=2),
    area: wp.array(dtype=wp.float64),
    vol6: wp.array(dtype=wp.float64),
):
    """Globalise the face indices and record each face's rest area and volume contribution.

    ``vol6`` holds ``a · (b × c)`` per face; the ``/6`` and the sum happen on the host, deterministically.
    The sum is translation-invariant for a closed mesh, so an off-origin centre needs no correction.
    """
    f = wp.tid()
    ia = node_lo + face_v[f, 0]
    ib = node_lo + face_v[f, 1]
    ic = node_lo + face_v[f, 2]
    face_idx[f, 0] = ia
    face_idx[f, 1] = ib
    face_idx[f, 2] = ic
    a = pos[ia]
    b = pos[ib]
    c = pos[ic]
    area[f] = wp.float64(0.5) * wp.length(wp.cross(b - a, c - a))
    vol6[f] = wp.dot(a, wp.cross(b, c))


@wp.kernel
def _gather_opposites(
    face_v: wp.array(dtype=wp.int32, ndim=2),
    face_e: wp.array(dtype=wp.int32, ndim=2),
    count: wp.array(dtype=wp.int32),
    opp: wp.array(dtype=wp.int32, ndim=2),
):
    """For every (face, edge) incidence, record the vertex opposite that edge.

    The slot is taken with an atomic, so WHICH of the two opposites lands first is scheduling-dependent
    — :func:`_write_hinges` removes that dependence by ordering the pair, so the built hinge list is
    reproducible even though this kernel is not.
    """
    f = wp.tid()
    for k in range(3):
        e = face_e[f, k]
        slot = wp.atomic_add(count, e, 1)
        if slot < 2:
            opp[e, slot] = face_v[f, (k + 2) % 3]


@wp.kernel
def _write_hinges(
    edge: wp.array(dtype=wp.int32, ndim=2),
    count: wp.array(dtype=wp.int32),
    opp: wp.array(dtype=wp.int32, ndim=2),
    node_lo: wp.int32,
    hinge: wp.array(dtype=wp.int32, ndim=2),
    bad: wp.array(dtype=wp.int32),
):
    """One dihedral per edge, ordered so the build does not depend on kernel scheduling.

    An edge with other than two incident faces is counted into ``bad`` rather than written: on a closed
    mesh it cannot happen, and if it does the surface is open or non-manifold and its hinge list would
    be quietly incomplete — which is a bending law with missing terms, not a warning.
    """
    e = wp.tid()
    if count[e] != 2:
        wp.atomic_add(bad, 0, 1)
        return
    o0 = opp[e, 0]
    o1 = opp[e, 1]
    hinge[e, 0] = node_lo + wp.min(o0, o1)
    hinge[e, 1] = node_lo + edge[e, 0]
    hinge[e, 2] = node_lo + edge[e, 1]
    hinge[e, 3] = node_lo + wp.max(o0, o1)


# ── the builder ─────────────────────────────────────────────────────────────────────────────────
def _base_topology() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The 12-vertex icosahedron plus its edge list and face-edge incidence, on the host.

    Twelve vertices and thirty edges: a constant, not a population.  Taken from
    ``surface.icosphere(0)`` rather than retyped, so this module and ``surface`` cannot drift into two
    different icosahedra — a divergence that would show up only as a small mismatch in a parity number
    much later.
    """
    verts, faces = icosphere(0)
    edge_id: dict[tuple[int, int], int] = {}
    face_e = np.empty((faces.shape[0], 3), np.int32)
    for f, tri in enumerate(faces):
        for k in range(3):
            u, v = int(tri[k]), int(tri[(k + 1) % 3])
            key = (min(u, v), max(u, v))
            if key not in edge_id:
                edge_id[key] = len(edge_id)
            face_e[f, k] = edge_id[key]
    edge = np.empty((len(edge_id), 2), np.int32)
    for (u, v), e in edge_id.items():
        edge[e] = (u, v)
    return verts, edge, faces.astype(np.int32), face_e


def build_closed_surface(
    arena: WorldArena,
    population: str,
    *,
    radius_um: float | None = None,
    mesh_um: float | None = None,
    mesh_provenance: str = "",
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> ClosedSurface:
    """Claim ONE range per kind and build a closed icosphere population into ``arena`` on the device.

    Args:
        arena: the world to claim from.  Must hold a CUDA allocation — there is no host build path.
        population: the name the ranges are attributed to, e.g. ``"membrane"``.
        radius_um: sphere radius [µm].  **REQUIRED, no default** — a surface has a size and a size is
            a declared physiological axis.
        mesh_um: the largest acceptable mean triangle edge [µm].  **REQUIRED, no default.**  The
            subdivision level is derived from it, so this is the only discretisation number a caller
            supplies and it is a sourced length rather than an index.
        mesh_provenance: where ``mesh_um`` came from.  Free text, carried into :meth:`ClosedSurface.record`
            so a PI-GAP label travels inside the artifact instead of beside it.
        centre_um: translation applied after scaling [µm].  A coordinate origin, not a physiological
            value, which is why this one may default.

    Returns:
        The built :class:`ClosedSurface`.

    Raises:
        ValueError: on a missing or non-positive radius or spacing, or a spacing no icosphere reaches.
        RuntimeError: if the arena holds no device allocation.
        AssertionError: if any edge is not shared by exactly two faces — an open or non-manifold mesh.
    """
    if radius_um is None:
        raise ValueError(
            "radius_um has no default and must be declared. A surface has a size, and a size is a "
            "physiological axis with a scope — not a number chosen in a builder."
        )
    if mesh_um is None:
        raise ValueError(
            "mesh_um has no default and must be declared. The subdivision level is DERIVED from it; a "
            "default spacing would recreate a typed subdivision level under a new name, which is the "
            "thing this signature exists to prevent."
        )
    radius = float(radius_um)
    if not (math.isfinite(radius) and radius > 0.0):
        raise ValueError(f"radius_um must be finite and positive; got {radius_um!r}")
    if not arena.node_arrays:
        raise RuntimeError(
            f"build_closed_surface({population}) needs a CUDA arena — this one was built with "
            "device=None, the bookkeeping-only mode. Construction is a Warp kernel (PI 2026-08-20); "
            "there is no host build path to fall back to."
        )

    level = subdivisions_for_mesh(radius, float(mesh_um))
    n_vert, n_edge, n_face = icosphere_counts(level)

    device = arena.node_arrays["position"].device
    with wp.ScopedDevice(device):
        base_v, base_e, base_fv, base_fe = _base_topology()
        vert = wp.array(base_v, dtype=wp.vec3d)
        edge = wp.array(base_e, dtype=wp.int32)
        face_v = wp.array(base_fv, dtype=wp.int32)
        face_e = wp.array(base_fe, dtype=wp.int32)
        v, e, f = 12, 30, 20

        for _ in range(level):
            v2, e2, f2 = v + e, 2 * e + 3 * f, 4 * f
            vert2 = wp.zeros(v2, dtype=wp.vec3d)
            edge2 = wp.zeros((e2, 2), dtype=wp.int32)
            face_v2 = wp.zeros((f2, 3), dtype=wp.int32)
            face_e2 = wp.zeros((f2, 3), dtype=wp.int32)
            wp.launch(_copy_vertices, dim=v, inputs=[vert, vert2])
            wp.launch(_split_edges, dim=e, inputs=[vert, edge, v, vert2, edge2])
            wp.launch(_split_faces, dim=f, inputs=[face_v, face_e, edge, v, e, edge2, face_v2, face_e2])
            vert, edge, face_v, face_e = vert2, edge2, face_v2, face_e2
            v, e, f = v2, e2, f2

        if (v, e, f) != (n_vert, n_edge, n_face):  # pragma: no cover — arithmetic, not a runtime path
            raise AssertionError(f"built ({v}, {e}, {f}) but the closed form says {(n_vert, n_edge, n_face)}")

        nodes = arena.claim(population, Kind.NODE, n_vert)
        faces_claim = arena.claim(population, Kind.FACE, n_face)
        hinges_claim = arena.claim(population, Kind.ANGLE4, n_edge)
        range_id = len(arena.claims(kind=Kind.NODE)) - 1

        wp.launch(
            _place_vertices,
            dim=n_vert,
            inputs=[
                vert, wp.float64(radius), wp.vec3d(*(float(x) for x in centre_um)),
                nodes.lo, range_id,
                arena.node_arrays["position"], arena.node_arrays["strand_id"],
                arena.node_arrays["range_id"],
            ],
        )

        face_idx = wp.zeros((n_face, 3), dtype=wp.int32)
        area0 = wp.zeros(n_face, dtype=wp.float64)
        vol6 = wp.zeros(n_face, dtype=wp.float64)
        wp.launch(
            _face_metrics,
            dim=n_face,
            inputs=[arena.node_arrays["position"], face_v, nodes.lo, face_idx, area0, vol6],
        )

        count = wp.zeros(n_edge, dtype=wp.int32)
        opp = wp.full((n_edge, 2), -1, dtype=wp.int32)
        hinge_idx = wp.zeros((n_edge, 4), dtype=wp.int32)
        bad = wp.zeros(1, dtype=wp.int32)
        wp.launch(_gather_opposites, dim=n_face, inputs=[face_v, face_e, count, opp])
        wp.launch(_write_hinges, dim=n_edge, inputs=[edge, count, opp, nodes.lo, hinge_idx, bad])

        n_bad = int(bad.numpy()[0])

    if n_bad:
        raise AssertionError(
            f"{population}: {n_bad} of {n_edge} edges are not shared by exactly two faces. The surface "
            "is open or non-manifold, so its hinge list is incomplete — a bending law with missing "
            "terms, which is worse than one that is wrong everywhere."
        )

    # Host-side, deterministic: an atomic float64 reduction would make the rest volume depend on the
    # order the blocks happened to retire, and the rest volume is what the pressure law reads.
    volume = float(np.sum(vol6.numpy(), dtype=np.float64) / 6.0)
    area_total = float(np.sum(area0.numpy(), dtype=np.float64))
    if volume <= 0.0:
        raise ValueError(
            f"{population} encloses {volume:.6g} um^3 — it is wound INWARD. An inverted face flips the "
            "sign of its own pressure traction, so this is asserted at build rather than discovered as "
            "a force pointing the wrong way."
        )

    topology_bytes = int(face_idx.size * 4 + area0.size * 8 + hinge_idx.size * 4)
    return ClosedSurface(
        population=population, nodes=nodes, faces=faces_claim, hinges=hinges_claim,
        subdivisions=level, radius_um=radius, centre_um=tuple(float(x) for x in centre_um),
        mesh_um_requested=float(mesh_um), mesh_um_realised=mean_edge_um(radius, level),
        mesh_provenance=mesh_provenance,
        face_idx=face_idx, area0_um2=area0, hinge_idx=hinge_idx,
        volume0_um3=volume, area0_total_um2=area_total, topology_bytes=topology_bytes,
    )


def build_membrane(
    arena: WorldArena,
    *,
    radius_um: float | None = None,
    mesh_um: float | None = None,
    mesh_provenance: str = "",
    centre_um: tuple[float, float, float] = (0.0, 0.0, 0.0),
    population: str = "membrane",
) -> ClosedSurface:
    """The plasma membrane, as a closed surface population named ``"membrane"``.

    The membrane is the cell's OUTER radius — PI 2026-07-17, *"measurement basis was the membrane"* —
    so ``radius_um`` here is the cell radius itself and the cortex is a shell just inside it, never the
    other way round.  This wrapper exists to hold that sentence next to the population name; it adds no
    behaviour to :func:`build_closed_surface`.

    Args:
        arena: the world to claim from.
        radius_um: the CELL radius [µm].  **REQUIRED** — a sourced axis, e.g. the MCF7 volume anchor in
            ``aleph.laws.cell_geometry``.
        mesh_um: largest acceptable mean triangle edge [µm].  **REQUIRED.**
        mesh_provenance: where ``mesh_um`` came from.
        centre_um: cell centre [µm].
        population: the claim name; only change it to build a second membrane.

    Returns:
        The built :class:`ClosedSurface`.
    """
    return build_closed_surface(
        arena, population, radius_um=radius_um, mesh_um=mesh_um,
        mesh_provenance=mesh_provenance, centre_um=centre_um,
    )


# ── self-check ──────────────────────────────────────────────────────────────────────────────────
def _host_reference_subdivide(level: int) -> tuple[np.ndarray, np.ndarray]:
    """The edge-list recursion, on the host, as an ORACLE for the index algebra.

    This is NOT a build path and must never become one — it exists so the arithmetic in
    :func:`_split_edges` and :func:`_split_faces` can be checked against
    ``surface.icosphere`` without a card, which is the only part of this module a device is not
    required to test.  It is a few hundred vertices in ``_demo`` and nothing else calls it.
    """
    verts, edge, face_v, face_e = _base_topology()
    verts = verts.astype(np.float64)
    for _ in range(level):
        v, e, f = verts.shape[0], edge.shape[0], face_v.shape[0]
        mid = 0.5 * (verts[edge[:, 0]] + verts[edge[:, 1]])
        verts = np.concatenate([verts, mid], axis=0)

        edge2 = np.empty((2 * e + 3 * f, 2), np.int32)
        m = v + np.arange(e, dtype=np.int32)
        edge2[0:2 * e:2] = np.stack([edge[:, 0], m], 1)
        edge2[1:2 * e:2] = np.stack([m, edge[:, 1]], 1)

        a, b, c = face_v[:, 0], face_v[:, 1], face_v[:, 2]
        e_ab, e_bc, e_ca = face_e[:, 0], face_e[:, 1], face_e[:, 2]
        m_ab, m_bc, m_ca = v + e_ab, v + e_bc, v + e_ca
        base = 2 * e + 3 * np.arange(f, dtype=np.int32)
        edge2[base + 0] = np.stack([m_ab, m_bc], 1)
        edge2[base + 1] = np.stack([m_bc, m_ca], 1)
        edge2[base + 2] = np.stack([m_ca, m_ab], 1)

        def child(ee, p):
            return np.where(edge[ee, 0] == p, 2 * ee, 2 * ee + 1).astype(np.int32)

        face_v2 = np.empty((4 * f, 3), np.int32)
        face_e2 = np.empty((4 * f, 3), np.int32)
        for slot, (tri, edges) in enumerate((
            ((a, m_ab, m_ca), (child(e_ab, a), base + 2, child(e_ca, a))),
            ((b, m_bc, m_ab), (child(e_bc, b), base + 0, child(e_ab, b))),
            ((c, m_ca, m_bc), (child(e_ca, c), base + 1, child(e_bc, c))),
            ((m_ab, m_bc, m_ca), (base + 0, base + 1, base + 2)),
        )):
            face_v2[slot::4] = np.stack(tri, 1)
            face_e2[slot::4] = np.stack(edges, 1)
        edge, face_v, face_e = edge2, face_v2, face_e2
    return verts / np.linalg.norm(verts, axis=1, keepdims=True), face_v


def _demo() -> None:
    """Self-check: the counts, the derived level, and the index algebra against the host icosphere."""
    from aleph.world.surface import enclosed_volume

    # The closed-form counts are Euler-consistent and reproduce the documented vertex series.
    assert icosphere_counts(0) == (12, 30, 20)
    assert [icosphere_counts(n)[0] for n in range(5)] == [12, 42, 162, 642, 2562]
    assert icosphere_counts(7) == (163_842, 491_520, 327_680)
    assert icosphere_counts(7)[1] == 3 * icosphere_counts(7)[2] // 2, "a closed mesh has 3F/2 edges"

    # The level is DERIVED from a spacing, and the derivation is what picks 7 at the cell radius.
    assert subdivisions_for_mesh(7.5, 0.100) == 7, "coarsest level inside the sourced 50-100 nm band"
    assert mean_edge_um(7.5, 6) > 0.100 >= mean_edge_um(7.5, 7), "6 is outside the band, 7 is inside"
    assert 0.050 <= mean_edge_um(7.5, 7) <= 0.100, "subdivision 7 lands INSIDE the band, not past it"
    # A different radius gives a different level from the same spacing — which is the whole point.
    assert subdivisions_for_mesh(5.1, 0.100) == 6, "the nuclear radius needs one level less"
    try:
        subdivisions_for_mesh(7.5, 1e-4)
    except ValueError as exc:
        assert "past a whole cell's budget" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an unreachable spacing must refuse")

    # The index algebra the kernels implement, checked against the host builder it must reproduce.
    for level in (1, 2, 3):
        v_ref, f_ref = icosphere(level)
        v_new, f_new = _host_reference_subdivide(level)
        assert v_new.shape == v_ref.shape and f_new.shape == f_ref.shape
        # Vertex ORDER differs — the edge list numbers midpoints by edge, the dict by face traversal —
        # so the meshes are compared as the sets of points and triangles they are.
        key = lambda p: np.lexsort(np.round(p.T, 9))  # noqa: E731
        assert np.allclose(v_new[key(v_new)], v_ref[key(v_ref)], atol=1e-12), \
            f"level {level}: the edge-list vertices are not the host's"
        assert np.isclose(enclosed_volume(v_new, f_new), enclosed_volume(v_ref, f_ref), rtol=1e-12)
        assert enclosed_volume(v_new, f_new) > 0.0, "child faces preserve the parent's outward winding"
        # Every edge used by exactly two faces: the closure the hinge kernel counts on the device.
        seen: dict[tuple[int, int], int] = {}
        for tri in f_new:
            for k in range(3):
                u, w = int(tri[k]), int(tri[(k + 1) % 3])
                seen[(min(u, w), max(u, w))] = seen.get((min(u, w), max(u, w)), 0) + 1
        assert set(seen.values()) == {2}, f"level {level} is not closed"
        assert len(seen) == icosphere_counts(level)[1]

    # Inverting the base inverts the volume — the control for the winding assertion above.
    v2, f2 = icosphere(2)
    assert enclosed_volume(v2, f2[:, ::-1]) < 0.0

    # Both lengths refuse to default, and a bookkeeping-only arena refuses to build.
    arena = WorldArena(capacity={Kind.NODE: 10, Kind.FACE: 10, Kind.ANGLE4: 10})
    for kwargs in ({"radius_um": 7.5}, {"mesh_um": 0.1}, {}):
        try:
            build_closed_surface(arena, "membrane", **kwargs)
        except ValueError as exc:
            assert "no default" in str(exc)
        else:  # pragma: no cover
            raise AssertionError(f"build_closed_surface{kwargs} must refuse")
    try:
        build_closed_surface(arena, "membrane", radius_um=7.5, mesh_um=0.1)
    except RuntimeError as exc:
        assert "needs a CUDA arena" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("a device=None arena must refuse to build")

    print(
        "membrane self-check OK — subdiv 7 at R=7.5 um: "
        f"{icosphere_counts(7)[0]:,} verts / {icosphere_counts(7)[2]:,} faces / "
        f"{icosphere_counts(7)[1]:,} hinges, mean edge {mean_edge_um(7.5, 7) * 1e3:.1f} nm"
    )


if __name__ == "__main__":
    _demo()
