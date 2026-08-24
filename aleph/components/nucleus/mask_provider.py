"""Nucleus relative-no-flux boundary provider — the §1.4 fluid-track interface (pure NumPy, NO Warp).

CRITICAL INTERFACE CONTRACT (AC_PARALLEL_SESSIONS §1.4). The fluid track OWNS
``ac/fluid/domain.py :: Domain``, which exposes ``set_nucleus_boundary(provider)`` where ``provider``
yields the inner relative-no-flux mask on the conservative field grid. I1a ships a STATIC-sphere
placeholder; THIS track ships ``DeformableNucleusMaskProvider`` conforming to the SAME signature — the
moving oblate mesh (about ``centre_nuc``) is I1a's inner boundary, NOT a static sphere R_nuc
(NEW_ENGINE_BUILD_PLAN I2 row: "mask the LIVE deformable oblate nucleus mesh"). Neither track edits the
other's file: we code against the ``NucleusBoundaryProvider`` Protocol + the local ``StaticSphere...``
test-double (which mirrors what I1a ships), so ac/fluid/ can adopt this drop-in.

The relative-no-flux BC: Darcy discharge relative to the moving solid boundary has zero normal
component across the envelope (q·n̂ = 0 in the boundary frame); the nucleoplasm exchanges no bulk pore
fluid with the cytoplasm on the poroelastic timescale. The provider supplies (i) the inside mask (grid
cells the nucleus occludes) and (ii) the live boundary faces + outward normals for the flux stencil.
Fluid–solid velocity coupling (the mesh moves with the mechanical solve) is wired at integration; here
the geometry is exact and re-queryable each outer step via ``update``.

Point-in-closed-mesh uses the signed-solid-angle winding number (Van Oosterom–Strackee) — robust for a
DEFORMED (non-ellipsoidal) envelope, exact for the oblate baseline, and orientation-checked against the
mesh winding. A fast analytic ellipsoid test is provided as the smooth-baseline cross-check.
"""
from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np
import numpy.typing as npt

from aleph.components.nucleus.geometry import mesh_volume

__all__ = [
    "NucleusBoundaryProvider",
    "StaticSphereMaskProvider",
    "DeformableNucleusMaskProvider",
    "inside_ellipsoid",
    "solid_angle_winding_number",
]


@runtime_checkable
class NucleusBoundaryProvider(Protocol):
    """The contract ``Domain.set_nucleus_boundary(provider)`` consumes each outer step."""

    def inside_mask(self, grid_points: npt.NDArray[np.float64]) -> npt.NDArray[np.bool_]:
        """Boolean mask (N,) — True where a grid point lies INSIDE the nucleus (relative no-flux)."""
        ...

    def boundary(self) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Current boundary: face centroids (Nf,3) [µm] and outward unit normals (Nf,3)."""
        ...


def inside_ellipsoid(
    grid_points: npt.NDArray[np.float64], centre, a: float, c: float
) -> npt.NDArray[np.bool_]:
    """Analytic oblate-ellipsoid occupancy ``(x²+y²)/a² + z²/c² ≤ 1`` (equatorial radius a, polar
    semi-axis c). Exact for the smooth oblate baseline; the winding-number cross-check target."""
    d = np.asarray(grid_points, dtype=np.float64) - np.asarray(centre, np.float64)[None, :]
    return (d[:, 0] ** 2 + d[:, 1] ** 2) / (a * a) + d[:, 2] ** 2 / (c * c) <= 1.0


def solid_angle_winding_number(
    grid_points: npt.NDArray[np.float64],
    verts: npt.NDArray[np.float64],
    faces: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    """Winding number w=Ω_total/4π of each point w.r.t. the closed oriented mesh (Van Oosterom–
    Strackee signed solid angle per triangle). w≈1 inside, w≈0 outside — robust for any closed
    orientable surface, no convexity assumption. Returns (N,) float."""
    gp = np.asarray(grid_points, dtype=np.float64)
    v0 = verts[faces[:, 0]]
    v1 = verts[faces[:, 1]]
    v2 = verts[faces[:, 2]]
    omega = np.zeros(gp.shape[0], dtype=np.float64)
    for n in range(gp.shape[0]):
        a = v0 - gp[n]
        b = v1 - gp[n]
        cc = v2 - gp[n]
        la = np.linalg.norm(a, axis=1)
        lb = np.linalg.norm(b, axis=1)
        lc = np.linalg.norm(cc, axis=1)
        num = np.einsum("ij,ij->i", a, np.cross(b, cc))          # a·(b×c) = signed 6·tet volume
        den = la * lb * lc + np.einsum("ij,ij->i", a, b) * lc \
            + np.einsum("ij,ij->i", a, cc) * lb + np.einsum("ij,ij->i", b, cc) * la
        omega[n] = 2.0 * np.sum(np.arctan2(num, den))
    return omega / (4.0 * np.pi)


class StaticSphereMaskProvider:
    """The I1a placeholder (test-double): a static sphere ``‖x−centre‖ ≤ radius`` no-flux inclusion.
    Mirrors what fluid-spine ships so this track can be tested against the frozen contract before
    ac/fluid/domain.py exists."""

    def __init__(self, centre, radius: float):
        self.centre = np.asarray(centre, np.float64)
        self.radius = float(radius)

    def inside_mask(self, grid_points: npt.NDArray[np.float64]) -> npt.NDArray[np.bool_]:
        d = np.asarray(grid_points, dtype=np.float64) - self.centre[None, :]
        return np.einsum("ij,ij->i", d, d) <= self.radius ** 2

    def boundary(self) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        # No mesh — the sphere boundary is analytic; return an empty face set (I1a uses the analytic BC).
        empty = np.zeros((0, 3), np.float64)
        return empty, empty


class DeformableNucleusMaskProvider:
    """The LIVE deformable oblate-mesh nucleus as I1a's relative-no-flux inner boundary (the I2 real
    provider). Conforms to ``NucleusBoundaryProvider``: ``inside_mask`` (winding-number occupancy of
    the current mesh) + ``boundary`` (face centroids + outward normals for the flux stencil).
    ``update(verts)`` refreshes the geometry when the mechanical solve moves the envelope — the mask
    tracks the moving surface, never a frozen R_nuc.

    Args:
        verts: (Nv,3) envelope node positions [µm] (the current deformed mesh).
        faces: (Nf,3) outward-wound triangles (topology is fixed; only ``verts`` moves).
    """

    def __init__(self, verts: npt.NDArray[np.float64], faces: npt.NDArray[np.int64]):
        self.faces = np.ascontiguousarray(faces, np.int64)
        self.update(verts)

    def update(self, verts: npt.NDArray[np.float64]) -> None:
        """Refresh with the current envelope nodes (call once per outer step; topology is invariant)."""
        self.verts = np.ascontiguousarray(verts, np.float64)
        if mesh_volume(self.verts, self.faces) < 0.0:
            raise ValueError("mesh winding is inward (negative volume): normals would point inward")

    def inside_mask(self, grid_points: npt.NDArray[np.float64]) -> npt.NDArray[np.bool_]:
        """True where a grid cell lies inside the CURRENT envelope (winding number ≥ ½)."""
        return solid_angle_winding_number(grid_points, self.verts, self.faces) >= 0.5

    def boundary(self) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Face centroids (Nf,3) + outward unit normals (Nf,3) of the current envelope."""
        p0 = self.verts[self.faces[:, 0]]
        p1 = self.verts[self.faces[:, 1]]
        p2 = self.verts[self.faces[:, 2]]
        centroids = (p0 + p1 + p2) / 3.0
        nrm = np.cross(p1 - p0, p2 - p0)
        nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-300
        return centroids, nrm
