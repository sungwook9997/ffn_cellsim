r"""Pore-pressure traction on the live plasma-membrane mesh (Warp-CUDA).

The Biot body coupling ``-alpha * grad(p)`` is zero for a spatially uniform resting pressure.  That is correct
in the bulk, but it does **not** remove the boundary traction exerted by the enclosed fluid on its membrane.
For every outward-wound live membrane triangle this module applies

    f_face = (p_inside - p_ext) * A * n_out,

distributed by the linear-triangle shape functions (one third to each vertex).  In the engine units
``1 Pa == 1 pN/um^2`` exactly, so pressure times triangle area is pN.  This is the mechanistic surface term
whose spherical continuum balance is Young-Laplace ``gamma = DeltaP * R / 2``; it does not resurrect the
retired scalar ``turgor_kernel`` and it reads the spatial live pressure field rather than a lumped pressure.

The inside pressure is reconstructed at the live triangle centroid from two mask-aware affine fits centred at
the first two inward finite-volume locations (``dx/2`` and ``dx``).  The one-sided extrapolation
``p_surface = 2 p(dx/2) - p(dx)`` evaluates the physical interior limit, reproduces linear fields exactly, and
is second-order for a smooth trace. OUTSIDE/NUCLEUS values cannot leak into either fit. A face with a
rank-deficient stencil increments a device diagnostic and receives no invented fallback load.

Sanity Gate:
    * Dimensions: ``DeltaP [pN/um^2] * A [um^2] = force [pN]``.
    * Boundary cases: ``DeltaP=0`` gives bit-zero force; positive inside pressure pushes outward; reversing
      pressure reverses every force.  A closed surface has zero net force under uniform pressure.
    * Conservation/work: the assembled nodal load is the weak linear-triangle form of ``DeltaP dV``;
      finite-difference volume work and ``sum_i f_i dot x_i = 3 DeltaP V`` arbitrate the sign.
    * Numerical: interpolation is mask-aware and the unresolved-face count must be zero.  The pressure
      stiffness scale used by the explicit mechanics CFL is derived as ``DeltaP * mean_edge`` [pN/um].
    * Residency: pressure, mask, live positions, face forces and diagnostics remain on CUDA in the inner loop.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.components.fluid.surface_trace import membrane_surface_pressure_trace

__all__ = [
    "MembranePressureTraction",
    "membrane_pressure_traction_kernel",
    "uniform_pressure_force_reference",
    "signed_volume",
]

@wp.kernel
def membrane_pressure_traction_kernel(
    pos: wp.array(dtype=wp.vec3d),
    faces: wp.array(dtype=wp.int32, ndim=2),
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    origin: wp.vec3d,
    dx: wp.float64,
    p_ext: wp.float64,
    unresolved_faces: wp.array(dtype=wp.int32),
    force: wp.array(dtype=wp.vec3d),
) -> None:
    """Add live spatial pressure traction to the three GLOBAL vertices of each membrane triangle."""
    t = wp.tid()
    i0 = faces[t, 0]
    i1 = faces[t, 1]
    i2 = faces[t, 2]
    x0 = pos[i0]
    x1 = pos[i1]
    x2 = pos[i2]
    area2 = wp.cross(x1 - x0, x2 - x0)  # 2 A n_out for outward-wound faces
    a2 = wp.length(area2)
    if a2 <= wp.float64(0.0):
        wp.atomic_add(unresolved_faces, 0, 1)
        return

    n_out = area2 / a2
    centroid = (x0 + x1 + x2) / wp.float64(3.0)

    trace = membrane_surface_pressure_trace(p, mask, origin, dx, centroid, n_out)
    if trace[1] == wp.float64(0.0):
        wp.atomic_add(unresolved_faces, 0, 1)
        return
    p_inside = trace[0]
    # A*n = area2/2; linear triangle shape function integral is 1/3 per vertex -> area2/6.
    f_vertex = (p_inside - p_ext) * area2 / wp.float64(6.0)
    wp.atomic_add(force, i0, f_vertex)
    wp.atomic_add(force, i1, f_vertex)
    wp.atomic_add(force, i2, f_vertex)


@dataclass
class MembranePressureTraction:
    """Device force primitive coupling one ``FieldGrid`` to the live membrane triangles."""

    grid: object
    faces_d: wp.array
    n_faces: int
    p_ext: float = 0.0

    def __post_init__(self) -> None:
        device = wp.get_device(self.grid.device)
        if not device.is_cuda:
            raise RuntimeError("MembranePressureTraction is Warp-CUDA-only")
        self.device = str(device)
        with wp.ScopedDevice(device):
            self.unresolved_faces_d = wp.zeros(1, dtype=wp.int32, device=device)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Add ``(p_inside-p_ext) n`` traction to the live membrane without any host readback."""
        g = self.grid
        origin = wp.vec3d(float(g.origin[0]), float(g.origin[1]), float(g.origin[2]))
        wp.launch(
            membrane_pressure_traction_kernel,
            dim=self.n_faces,
            inputs=[pos, self.faces_d, g.p, g.mask, origin, wp.float64(g.dx), wp.float64(self.p_ext),
                    self.unresolved_faces_d],
            outputs=[force],
            device=self.device,
        )

    def reset_diagnostics(self) -> None:
        """Explicitly start a new gate/run window; force assembly never erases an earlier failure."""
        self.unresolved_faces_d.zero_()


def signed_volume(verts: npt.NDArray[np.float64], faces: npt.NDArray[np.integer]) -> float:
    """Signed volume of an outward-wound closed triangular mesh [um^3]."""
    x = np.asarray(verts, dtype=np.float64)
    f = np.asarray(faces, dtype=np.int64)
    return float(np.einsum("ij,ij->i", x[f[:, 0]], np.cross(x[f[:, 1]], x[f[:, 2]])).sum() / 6.0)


def uniform_pressure_force_reference(
    verts: npt.NDArray[np.float64], faces: npt.NDArray[np.integer], delta_p: float,
) -> npt.NDArray[np.float64]:
    """Pure-NumPy weak triangle load for a uniform pressure jump (analytic acceptance oracle)."""
    x = np.asarray(verts, dtype=np.float64)
    f = np.asarray(faces, dtype=np.int64)
    area2 = np.cross(x[f[:, 1]] - x[f[:, 0]], x[f[:, 2]] - x[f[:, 0]])
    per_vertex = float(delta_p) * area2 / 6.0
    out = np.zeros_like(x)
    np.add.at(out, f[:, 0], per_vertex)
    np.add.at(out, f[:, 1], per_vertex)
    np.add.at(out, f[:, 2], per_vertex)
    return out
