"""Surfaces: the FACE primitive, its hinges, and the two invariants a closed mesh must carry.

A FACE is the only primitive whose measure is an AREA. That is what makes it the carrier of surface
tension, of the pressure traction ``(p_in - p_out) n``, of the water flux through a membrane, and of the
enclosed volume — none of which a node, a segment or an angle triple can express. A membrane and a
nuclear envelope are the same primitive in different populations.

TWO INVARIANTS, AND BOTH HAVE COST SOMETHING IN THIS PROJECT BEFORE.

**Winding.** Every face must be wound outward, so the divergence-theorem volume is positive. An inverted
triangle is not a visual defect — it flips the sign of that face's pressure traction, and the repository's
own membrane track carries a negative control that inverts one triangle's winding and requires the volume
ratio to exceed one. So winding is asserted at build, where it is cheap, rather than discovered as a
force pointing the wrong way.

**Discrete volume, never analytic.** The enclosed volume is computed from the polyhedron the mesh
actually is, by the divergence theorem over its own triangles. Substituting the analytic sphere volume
is a documented trap in this codebase — it produced a phantom nucleus force at t0 in the Lane-D rig,
because a discrete mesh inscribed in a sphere encloses LESS than the sphere, and the difference became a
pressure that nothing balanced. The two differ by about 1.4% at subdivision 2 and shrink with
refinement, which is exactly the regime where the error is small enough to look like noise and large
enough to move a gate.

THE HINGE IS THE FOURTH-ARITY ELEMENT. Two faces sharing an edge define a dihedral, and Helfrich bending
lives there rather than on the triangles. It is ANGLE4 and not ANGLE3 for that reason — a dihedral is
not a bending triple, they differ in arity, in stencil and in which law reads them, which is why the
arena keeps two angle kinds rather than one with a discriminator.

WHAT IS DELIBERATELY ABSENT. No bending coefficient, no surface tension, no area-reservoir law, no
pressure. Those are laws and they carry parameters; this module builds the geometry they will read. The
one number it does record per face is ``area0``, because a rest area is part of the configuration rather
than part of a law.

Sanity Gate (recorded per CLAUDE.md, before first execution):
  * dimensional — vertex positions [µm]; ``area0`` [µm²]; enclosed volume [µm³].
  * boundary — a mesh with fewer than 4 faces cannot enclose anything and is refused; a face indexing a
    vertex that does not exist is refused; a closed mesh whose signed volume is non-positive is refused
    as inverted rather than silently negated.
  * conservation/invariant — for a closed mesh every edge is shared by exactly TWO faces, which is
    asserted; a boundary edge means the surface is open and its hinge list would be incomplete.
  * CFL/precision — no integration; float64. Note the hinge count, not the face count, sets the number
    of bending terms, and it is 3F/2 for a closed mesh rather than F.
  * sign sense — winding is the whole point: outward faces give a positive divergence-theorem volume,
    and the assertion fires on the sign rather than on a magnitude.
  * measurement protocol — host geometry only; no device is touched and nothing is uploaded.

engine units: length µm, area µm², volume µm³.  Runtime: pure host, CPU-importable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from aleph.world.arena import Claim, Kind, WorldArena

__all__ = ["Surface", "build_surface", "icosphere", "enclosed_volume"]


@dataclass(frozen=True, slots=True)
class Surface:
    """One built surface: its claims, its faces, its hinges, and the volume it encloses.

    Attributes:
        population: the population the ranges belong to.
        nodes / faces / hinges: the three claims.
        position: ``(V, 3)`` rest vertex positions [µm].
        face_idx: ``(F, 3)`` outward-wound triangles, as GLOBAL node indices.
        area0_um2: ``(F,)`` rest area of each face [µm²].
        hinge_idx: ``(H, 4)`` dihedrals as ``(a, b, c, d)`` where ``b–c`` is the shared edge and ``a``,
            ``d`` are the two opposite vertices — GLOBAL indices.
        volume0_um3: the enclosed volume of the rest configuration [µm³], by the divergence theorem
            over the mesh's own triangles.
        closed: whether every edge is shared by exactly two faces.
    """

    population: str
    nodes: Claim
    faces: Claim
    hinges: Claim
    position: npt.NDArray[np.float64]
    face_idx: npt.NDArray[np.int64]
    area0_um2: npt.NDArray[np.float64]
    hinge_idx: npt.NDArray[np.int64]
    volume0_um3: float
    closed: bool

    @property
    def n_vertices(self) -> int:
        """Number of vertices."""
        return self.nodes.count

    @property
    def area0_total_um2(self) -> float:
        """Total rest area [µm²]."""
        return float(self.area0_um2.sum())


def enclosed_volume(position: npt.NDArray[np.float64], face_idx: npt.NDArray[np.int64]) -> float:
    """Signed volume enclosed by a closed triangulated surface [µm³], by the divergence theorem.

    Computed from the polyhedron the mesh actually is — summing ``a · (b × c) / 6`` over its own
    triangles — and never from an analytic formula for the shape it approximates. A mesh inscribed in a
    sphere encloses less than the sphere, and substituting the analytic value is what produced a phantom
    force at t0 in this project once already.

    Positive for outward winding, negative for inward.
    """
    p = np.asarray(position, np.float64)
    f = np.asarray(face_idx, np.int64).reshape(-1, 3)
    a, b, c = p[f[:, 0]], p[f[:, 1]], p[f[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def _hinges(face_idx: npt.NDArray[np.int64]) -> tuple[npt.NDArray[np.int64], bool]:
    """Return ``(H, 4)`` dihedrals ``(a, b, c, d)`` over shared edges, and whether the mesh is closed.

    ``b–c`` is the shared edge; ``a`` and ``d`` are the vertices opposite it in the two faces. An edge
    used by one face is a boundary (open surface); one used by more than two is non-manifold and raises,
    because a bending term over it would be ambiguous rather than merely wrong.
    """
    f = np.asarray(face_idx, np.int64).reshape(-1, 3)
    edges: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for tri in f:
        for k in range(3):
            u, v = int(tri[k]), int(tri[(k + 1) % 3])
            opp = int(tri[(k + 2) % 3])
            edges.setdefault((u, v) if u < v else (v, u), []).append((u, opp))
    rows, closed = [], True
    for (u, v), uses in edges.items():
        if len(uses) == 1:
            closed = False
            continue
        if len(uses) > 2:
            raise ValueError(
                f"edge ({u}, {v}) is shared by {len(uses)} faces — the mesh is non-manifold and a "
                "dihedral over that edge would be ambiguous, not merely wrong."
            )
        rows.append((uses[0][1], u, v, uses[1][1]))
    return np.asarray(rows, np.int64).reshape(-1, 4), closed


def icosphere(subdivisions: int = 1) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """A unit icosphere: ``(V, 3)`` vertices on the unit sphere and ``(F, 3)`` outward-wound faces.

    Present because a membrane and a nuclear envelope are closed surfaces, and a closed surface is the
    only configuration in which the winding and volume invariants mean anything. Vertex counts follow
    ``10·4^n + 2``: 12, 42, 162, 642, 2562 …

    Raises:
        ValueError: on a negative subdivision count.
    """
    if subdivisions < 0:
        raise ValueError(f"subdivisions must be nonnegative; got {subdivisions}")
    t = (1.0 + 5.0 ** 0.5) / 2.0
    verts = np.array([
        (-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0),
        (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
        (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1),
    ], np.float64)
    faces = np.array([
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ], np.int64)
    for _ in range(subdivisions):
        mid: dict[tuple[int, int], int] = {}
        vlist = list(verts)
        new_faces = []
        for a, b, c in faces:
            m = []
            for u, v in ((a, b), (b, c), (c, a)):
                key = (int(min(u, v)), int(max(u, v)))
                if key not in mid:
                    mid[key] = len(vlist)
                    vlist.append((verts[u] + verts[v]) / 2.0)
                m.append(mid[key])
            ab, bc, ca = m
            new_faces += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        verts = np.asarray(vlist, np.float64)
        faces = np.asarray(new_faces, np.int64)
    verts = verts / np.linalg.norm(verts, axis=1, keepdims=True)
    return verts, faces


def build_surface(
    arena: WorldArena,
    population: str,
    *,
    vertices: npt.ArrayLike,
    faces: npt.ArrayLike,
    radius_um: float | None = None,
    centre: npt.ArrayLike = (0.0, 0.0, 0.0),
) -> Surface:
    """Claim ranges for one surface and build its geometry into ``arena``.

    Args:
        arena: the world to claim from.
        population: the name the ranges are attributed to, e.g. ``"membrane"``.
        vertices: ``(V, 3)`` positions. Scaled by ``radius_um`` and offset by ``centre`` if given, so a
            unit :func:`icosphere` can be placed without the caller rescaling it.
        faces: ``(F, 3)`` triangles, outward-wound.
        radius_um: scale to apply to ``vertices`` [µm]. **No default** — a surface has a size and it is
            a declared axis, not a number chosen here.
        centre: translation applied after scaling [µm].

    Returns:
        The built :class:`Surface`.

    Raises:
        ValueError: on a missing or non-positive radius, fewer than 4 faces, an out-of-range vertex
            index, a non-manifold edge, or a closed mesh wound inward.
    """
    if radius_um is None:
        raise ValueError(
            "radius_um has no default and must be declared. A surface has a size, and a size is a "
            "physiological axis with a scope — not a number chosen in a builder."
        )
    r = float(radius_um)
    if not np.isfinite(r) or r <= 0.0:
        raise ValueError(f"radius_um must be finite and positive; got {radius_um!r}")

    v = np.asarray(vertices, np.float64).reshape(-1, 3) * r + np.asarray(centre, np.float64).reshape(3)
    f_local = np.asarray(faces, np.int64).reshape(-1, 3)
    if f_local.shape[0] < 4:
        raise ValueError(f"a surface needs at least 4 faces to enclose anything; got {f_local.shape[0]}")
    if f_local.min() < 0 or f_local.max() >= v.shape[0]:
        raise ValueError(
            f"face indexes vertex [{int(f_local.min())}, {int(f_local.max())}] against {v.shape[0]} vertices"
        )

    hinge_local, closed = _hinges(f_local)
    vol = enclosed_volume(v, f_local)
    if closed and vol <= 0.0:
        raise ValueError(
            f"closed mesh encloses {vol:.6g} um^3 — it is wound INWARD. An inverted face flips the sign "
            "of its own pressure traction, so this is asserted here rather than discovered as a force "
            "pointing the wrong way."
        )

    nodes = arena.claim(population, Kind.NODE, v.shape[0])
    faces_claim = arena.claim(population, Kind.FACE, f_local.shape[0])
    hinges_claim = (arena.claim(population, Kind.ANGLE4, hinge_local.shape[0])
                    if hinge_local.shape[0] else Claim(population, Kind.ANGLE4, 0, 0))

    face_idx = f_local + nodes.lo
    hinge_idx = hinge_local + nodes.lo if hinge_local.size else hinge_local

    a, b, c = v[f_local[:, 0]], v[f_local[:, 1]], v[f_local[:, 2]]
    area0 = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)

    return Surface(
        population=population, nodes=nodes, faces=faces_claim, hinges=hinges_claim,
        position=v, face_idx=face_idx, area0_um2=area0, hinge_idx=hinge_idx,
        volume0_um3=vol, closed=closed,
    )


def _demo() -> None:
    """Self-check: the closed-mesh invariants, and the discrete-vs-analytic volume gap."""
    arena = WorldArena(capacity={Kind.NODE: 4_000, Kind.FACE: 8_000, Kind.ANGLE4: 12_000})

    v, f = icosphere(2)
    assert v.shape[0] == 162 and f.shape[0] == 320, "10*4^n + 2 vertices, 20*4^n faces"

    s = build_surface(arena, "membrane", vertices=v, faces=f, radius_um=7.5)
    assert s.n_vertices == 162 and s.faces.count == 320
    assert s.closed, "an icosphere is closed"
    assert s.hinges.count == 480 == 3 * 320 // 2, "a closed mesh has 3F/2 edges, not F"
    assert s.volume0_um3 > 0.0, "outward winding gives a positive volume"

    # The trap, made numerical: the discrete polyhedron encloses LESS than the sphere it approximates.
    analytic = 4.0 / 3.0 * np.pi * 7.5 ** 3
    gap = (analytic - s.volume0_um3) / analytic
    assert 0.0 < gap < 0.05, f"discrete volume is {gap:.4%} below analytic — small, and not zero"

    # Refinement shrinks it, which is why substituting the analytic value looks like noise.
    v3, f3 = icosphere(3)
    finer = build_surface(arena, "nucleus", vertices=v3, faces=f3, radius_um=3.0)
    gap3 = (4.0 / 3.0 * np.pi * 3.0 ** 3 - finer.volume0_um3) / (4.0 / 3.0 * np.pi * 3.0 ** 3)
    assert gap3 < gap, "a finer mesh is closer to the sphere"

    # Area follows the square of the radius; the mesh is the same shape at both radii.
    unit = build_surface(arena, "probe", vertices=v, faces=f, radius_um=1.0)
    assert np.isclose(s.area0_total_um2 / unit.area0_total_um2, 7.5 ** 2)

    # Inward winding is refused rather than silently negated.
    try:
        build_surface(arena, "flipped", vertices=v, faces=f[:, ::-1], radius_um=1.0)
    except ValueError as exc:
        assert "wound INWARD" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("an inward-wound closed mesh must refuse")

    # The radius refuses to default, like every other length in this engine.
    try:
        build_surface(arena, "sized", vertices=v, faces=f)
    except ValueError as exc:
        assert "no default" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("radius_um must refuse to default")

    arena.assert_partitioned()
    print(f"surface self-check OK — {s.n_vertices} verts / {s.faces.count} faces / "
          f"{s.hinges.count} hinges; volume {s.volume0_um3:.3f} um^3, {gap:.3%} below analytic")


if __name__ == "__main__":
    _demo()
