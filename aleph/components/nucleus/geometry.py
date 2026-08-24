"""Closed-mesh geometry helpers for the deformable nucleus — pure NumPy (ZERO Warp import).

Host-side geometry the analytic oracles and the fluid mask-provider share: an oblate-capable
triangulated envelope, dihedral-hinge adjacency, the divergence-theorem signed volume + its exact
vertex gradient, and areas/normals. These are independent re-derivations (not re-imports of the
kernels under test) so they can serve as the FD-gradient ground truth for the Warp force kernels.

Conventions: FF units — µm, pN, s. An oblate spheroid has equatorial radius ``a = R_eq`` and polar
semi-axis ``c = R_eq / aspect`` (aspect = a/c ≥ 1 flattens along z); volume = (4/3)π a² c.
Icosphere winding is outward, so the signed volume is positive.
"""
from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.laws.surface_manifold import icosphere

__all__ = [
    "build_oblate_mesh",
    "build_hinges",
    "mesh_volume",
    "mesh_volume_gradient",
    "mesh_area",
    "oblate_volume",
    "face_normals_areas",
]


def build_oblate_mesh(
    r_eq: float, aspect: float = 1.0, subdivisions: int = 3, centre=(0.0, 0.0, 0.0)
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """Oblate triangulated envelope by scaling an icosphere's z-axis by ``1/aspect``.

    Args:
        r_eq: Equatorial radius ``a`` [µm] (> 0).
        aspect: Oblate aspect ratio ``a/c`` ≥ 1 (1 = sphere; 1.5–3 = adherent-flattened nucleus).
        subdivisions: Icosphere subdivision level (3 → 1280 faces / 642 nodes; nucleus needs less
            resolution than the cortex — it is a continuum shell).
        centre: Envelope centre ``centre_nuc`` [µm].

    Returns:
        ``(verts, faces)`` — verts (Nv,3) float64 [µm], faces (Nf,3) int64 outward-wound.

    Raises:
        ValueError: If ``r_eq`` ≤ 0 or ``aspect`` < 1.
    """
    if not (np.isfinite(r_eq) and r_eq > 0.0):
        raise ValueError(f"r_eq must be finite > 0; got {r_eq!r}")
    if aspect < 1.0:
        raise ValueError(f"aspect must be ≥ 1 (oblate flattens along z); got {aspect!r}")
    v, f = icosphere(subdivisions, r_eq)
    v = np.array(v, dtype=np.float64)
    v[:, 2] /= float(aspect)                       # flatten the pole → oblate
    v += np.asarray(centre, np.float64)[None, :]
    return np.ascontiguousarray(v), np.ascontiguousarray(f, np.int64)


def build_hinges(faces: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
    """Dihedral-hinge adjacency ``(Nh,4)=[i,j,k,l]`` from DIRECTED half-edges (winding is the
    correctness-critical detail — mirrors ``ff.membrane_surface.build_membrane_hinges``, re-derived
    here to keep the analytic layer Warp-free).

    For a manifold edge {i,j} the two half-edges i→j and j→i live in triangles T1=(i,j,k) and
    T2=(j,i,l), so the two flap vertices carry each triangle's outward winding ⇒ consistent normals.
    Closed mesh: Nh = 3·Nf/2.
    """
    opp: dict[tuple[int, int], int] = {}
    for tri in faces:
        v0, v1, v2 = int(tri[0]), int(tri[1]), int(tri[2])
        opp[(v0, v1)] = v2
        opp[(v1, v2)] = v0
        opp[(v2, v0)] = v1
    hinges: list[tuple[int, int, int, int]] = []
    seen: set[tuple[int, int]] = set()
    for (a, b), _ in opp.items():
        key = (a, b) if a < b else (b, a)
        if key in seen or (b, a) not in opp:
            continue
        seen.add(key)
        i, j = key
        hinges.append((i, j, opp[(i, j)], opp[(j, i)]))
    return np.ascontiguousarray(hinges, np.int64)


def face_normals_areas(
    verts: npt.NDArray[np.float64], faces: npt.NDArray[np.int64]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Per-face outward unit normals (Nf,3) and areas (Nf,) [µm²]."""
    p0 = verts[faces[:, 0]]
    p1 = verts[faces[:, 1]]
    p2 = verts[faces[:, 2]]
    cr = np.cross(p1 - p0, p2 - p0)
    a2 = np.linalg.norm(cr, axis=1)
    n = cr / (a2[:, None] + 1e-300)
    return n, 0.5 * a2


def mesh_area(verts: npt.NDArray[np.float64], faces: npt.NDArray[np.int64]) -> float:
    """Total surface area Σ triangle areas [µm²]."""
    _, areas = face_normals_areas(verts, faces)
    return float(areas.sum())


def mesh_volume(verts: npt.NDArray[np.float64], faces: npt.NDArray[np.int64]) -> float:
    """Enclosed volume of the closed mesh via the divergence theorem [µm³]:
    ``V = (1/6) Σ_faces x_i · (x_j × x_k)`` (signed; positive for outward winding). Translation-
    invariant for a CLOSED mesh (the shift terms telescope to zero over the closed surface)."""
    p0 = verts[faces[:, 0]]
    p1 = verts[faces[:, 1]]
    p2 = verts[faces[:, 2]]
    return float(np.einsum("ij,ij->i", p0, np.cross(p1, p2)).sum() / 6.0)


def mesh_volume_gradient(
    verts: npt.NDArray[np.float64], faces: npt.NDArray[np.int64]
) -> npt.NDArray[np.float64]:
    """Exact vertex gradient ``∂V/∂x_i`` (Nv,3) [µm²]. For vertex i on face (i,j,k),
    ``∂V/∂x_i = (1/6)(x_j × x_k)``; summed over all incident faces. This is the area-weighted
    outward normal and is the force direction of the nucleoplasm volume constraint."""
    grad = np.zeros_like(verts)
    p0 = verts[faces[:, 0]]
    p1 = verts[faces[:, 1]]
    p2 = verts[faces[:, 2]]
    np.add.at(grad, faces[:, 0], np.cross(p1, p2) / 6.0)
    np.add.at(grad, faces[:, 1], np.cross(p2, p0) / 6.0)
    np.add.at(grad, faces[:, 2], np.cross(p0, p1) / 6.0)
    return grad


def oblate_volume(r_eq: float, aspect: float) -> float:
    """Continuum oblate-spheroid volume ``(4/3)π a² c`` with ``a=r_eq``, ``c=r_eq/aspect`` [µm³]."""
    c = r_eq / aspect
    return float(4.0 / 3.0 * np.pi * r_eq * r_eq * c)
