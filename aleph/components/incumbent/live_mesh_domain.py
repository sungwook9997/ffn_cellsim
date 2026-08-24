"""GPU-resident live-mesh classifiers for the Active Cell fluid domain.

The conservative Biot/RAD grid is the volume enclosed by the *current* plasma-membrane mesh minus the
*current* nuclear-envelope mesh.  A production classifier therefore cannot use the historical static spheres,
nor may it download the mechanical mesh to NumPy for a host winding-number pass.  This module keeps a
``wp.Mesh`` acceleration structure for each closed surface, refreshes it from the authoritative float64
mechanical positions on the same CUDA device, and classifies every finite-volume cell with Warp's winding-
number point query.

``wp.Mesh`` stores float32 points, while the mechanical state is float64.  The conversion is a device-to-device
copy and is used only for the topological inside/outside query.  At the production grid spacing (0.5 um), the
float32 coordinate round-off at the 10.5 um domain extent is < 2e-6 of one cell.  The native parity gate still
requires exact cell-for-cell agreement with the independent float64 host winding oracle at rest and after a
prescribed deformation; a disagreement is a numerical finding, never a tolerance to loosen.

Sanity Gate:
    * Dimensions: positions, ``origin``, ``dx`` and ``max_dist`` are all [um]; the winding sign is
      dimensionless.  No material magnitude enters classification.
    * Boundary cases: a watertight outward-wound surface gives sign < 0 inside and sign > 0 outside.  The
      topological 0.5 winding threshold is Warp's closed-surface definition, not a fitted cutoff. Grid centres
      on the surface use the closed-set policy (inside); they are counted separately with an epsilon derived
      from float32 round-off, not silently absorbed into a biological tolerance.
    * Conservation: these kernels only classify.  ``Domain.remap`` remains the sole owner of moving-face
      content transport when a cell changes class.
    * Numerical: ``max_dist`` is the finite-volume box diagonal, derived from grid geometry.  There is no
      hard-coded query radius or device id.
    * Residency: refresh = float64->float32 D2D copy + device BVH refit + device query.  There is no ``.numpy()``
      or authoritative host geometry in the physical-time loop. Query failures are device-counted and must be
      zero at the native gate.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import warp as wp

__all__ = [
    "LiveMeshMembraneMaskProvider",
    "LiveMeshNucleusMaskProvider",
    "copy_live_mesh_points_kernel",
    "classify_live_membrane_kernel",
    "classify_live_nucleus_kernel",
]

_FLUID = wp.constant(1)
_OUTSIDE = wp.constant(0)
_NUCLEUS = wp.constant(2)


def _validate_closed_surface(faces: np.ndarray, n_verts: int, reference_verts: np.ndarray | None) -> None:
    """Reject non-manifold, open, degenerate, or inward-wound topology before creating the query mesh."""
    if np.any(faces[:, 0] == faces[:, 1]) or np.any(faces[:, 1] == faces[:, 2]) \
            or np.any(faces[:, 2] == faces[:, 0]):
        raise ValueError("live boundary contains a triangle with repeated vertices")
    directed = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=0)
    undirected = np.sort(directed, axis=1)
    unique, counts = np.unique(undirected, axis=0, return_counts=True)
    if unique.shape[0] == 0 or np.any(counts != 2):
        raise ValueError("live boundary must be a watertight 2-manifold (every edge shared by two faces)")
    # Every shared edge must occur once in each direction; equal directions expose inconsistent winding.
    edge_balance: dict[tuple[int, int], int] = {}
    for a, b in directed:
        lo, hi = (int(a), int(b)) if a < b else (int(b), int(a))
        edge_balance[(lo, hi)] = edge_balance.get((lo, hi), 0) + (1 if a < b else -1)
    if any(balance != 0 for balance in edge_balance.values()):
        raise ValueError("live boundary triangle winding is inconsistent across a shared edge")
    if reference_verts is not None:
        x = np.ascontiguousarray(reference_verts, dtype=np.float64)
        if x.shape != (n_verts, 3):
            raise ValueError("reference_verts must have shape (n_verts, 3)")
        area2 = np.linalg.norm(np.cross(x[faces[:, 1]] - x[faces[:, 0]],
                                       x[faces[:, 2]] - x[faces[:, 0]]), axis=1)
        if np.any(area2 <= 0.0):
            raise ValueError("live boundary contains a zero-area reference triangle")
        volume6 = float(np.einsum("ij,ij->i", x[faces[:, 0]],
                                  np.cross(x[faces[:, 1]], x[faces[:, 2]])).sum())
        if volume6 <= 0.0:
            raise ValueError("live boundary must be consistently outward-wound (positive signed volume)")


@wp.kernel
def copy_live_mesh_points_kernel(
    pos_global: wp.array(dtype=wp.vec3d), node_off: wp.int32, points_local: wp.array(dtype=wp.vec3),
) -> None:
    """Copy one closed-surface node block from the float64 mechanics state to the float32 query mesh (D2D)."""
    i = wp.tid()
    p = pos_global[node_off + i]
    points_local[i] = wp.vec3(wp.float32(p[0]), wp.float32(p[1]), wp.float32(p[2]))


@wp.kernel
def classify_live_membrane_kernel(
    mesh_id: wp.uint64,
    origin: wp.vec3,
    dx: wp.float32,
    max_dist: wp.float32,
    surface_epsilon: wp.float32,
    query_failures: wp.array(dtype=wp.int32),
    surface_ties: wp.array(dtype=wp.int32),
    query_failure_mask: wp.array3d(dtype=wp.int32),
    surface_tie_mask: wp.array3d(dtype=wp.int32),
    mask: wp.array3d(dtype=wp.int32),
) -> None:
    """Write FLUID inside the current membrane mesh and OUTSIDE beyond it."""
    i, j, k = wp.tid()
    query_failure_mask[i, j, k] = 0
    surface_tie_mask[i, j, k] = 0
    point = origin + dx * wp.vec3(wp.float32(i), wp.float32(j), wp.float32(k))
    query = wp.mesh_query_point_sign_winding_number(
        mesh_id, point, max_dist, wp.float32(2.0), wp.float32(0.5))
    if not query.result:
        wp.atomic_add(query_failures, 0, 1)
        query_failure_mask[i, j, k] = 1
        mask[i, j, k] = _OUTSIDE
        return
    closest = wp.mesh_eval_position(mesh_id, query.face, query.u, query.v)
    on_surface = wp.length(closest - point) <= surface_epsilon
    if on_surface:
        wp.atomic_add(surface_ties, 0, 1)
        surface_tie_mask[i, j, k] = 1
    if query.sign < wp.float32(0.0) or on_surface:
        mask[i, j, k] = _FLUID
    else:
        mask[i, j, k] = _OUTSIDE


@wp.kernel
def classify_live_nucleus_kernel(
    mesh_id: wp.uint64,
    origin: wp.vec3,
    dx: wp.float32,
    max_dist: wp.float32,
    surface_epsilon: wp.float32,
    query_failures: wp.array(dtype=wp.int32),
    surface_ties: wp.array(dtype=wp.int32),
    query_failure_mask: wp.array3d(dtype=wp.int32),
    surface_tie_mask: wp.array3d(dtype=wp.int32),
    mask: wp.array3d(dtype=wp.int32),
) -> None:
    """Overwrite membrane-interior cells enclosed by the current nuclear envelope as NUCLEUS."""
    i, j, k = wp.tid()
    query_failure_mask[i, j, k] = 0
    surface_tie_mask[i, j, k] = 0
    if mask[i, j, k] != _FLUID:
        return
    point = origin + dx * wp.vec3(wp.float32(i), wp.float32(j), wp.float32(k))
    query = wp.mesh_query_point_sign_winding_number(
        mesh_id, point, max_dist, wp.float32(2.0), wp.float32(0.5))
    if not query.result:
        wp.atomic_add(query_failures, 0, 1)
        query_failure_mask[i, j, k] = 1
        return
    closest = wp.mesh_eval_position(mesh_id, query.face, query.u, query.v)
    on_surface = wp.length(closest - point) <= surface_epsilon
    if on_surface:
        wp.atomic_add(surface_ties, 0, 1)
        surface_tie_mask[i, j, k] = 1
    if query.sign < wp.float32(0.0) or on_surface:
        mask[i, j, k] = _NUCLEUS


class _LiveMeshGeometry:
    """Device-resident classifier for one live, closed triangular surface.

    Args:
        pos_global: Authoritative combined mechanical position array (float64, CUDA).
        node_off: First vertex of this surface in ``pos_global``.
        n_verts: Number of contiguous surface vertices.
        faces_local: Outward-wound local triangle indices in ``[0, n_verts)``.
    """

    def __init__(
        self,
        pos_global: wp.array,
        node_off: int,
        n_verts: int,
        faces_local: npt.NDArray[np.integer],
        reference_verts: npt.NDArray[np.floating] | None = None,
    ) -> None:
        device = pos_global.device
        if not device.is_cuda:
            raise RuntimeError("live mesh classification is Warp-CUDA-only; CPU classification is forbidden")
        if pos_global.dtype != wp.vec3d or pos_global.ndim != 1 or not pos_global.is_contiguous:
            raise TypeError("pos_global must be a contiguous 1-D wp.vec3d array")
        if int(n_verts) <= 0:
            raise ValueError("a live closed surface must contain at least one vertex")
        if int(node_off) < 0 or int(node_off) + int(n_verts) > int(pos_global.shape[0]):
            raise ValueError("live surface node block lies outside pos_global")
        faces = np.ascontiguousarray(faces_local, dtype=np.int32)
        if faces.ndim != 2 or faces.shape[1] != 3:
            raise ValueError("faces_local must have shape (n_faces, 3)")
        if faces.shape[0] == 0:
            raise ValueError("a live closed surface must contain at least one triangle")
        if faces.size and (int(faces.min()) < 0 or int(faces.max()) >= int(n_verts)):
            raise ValueError("faces_local contains an index outside the live surface node block")
        _validate_closed_surface(faces, int(n_verts), None if reference_verts is None else np.asarray(reference_verts))
        self.pos_global = pos_global
        self.node_off = int(node_off)
        self.n_verts = int(n_verts)
        self.device = device
        with wp.ScopedDevice(device):
            self.points_d = wp.zeros(self.n_verts, dtype=wp.vec3, device=device)
            self.indices_d = wp.array(faces.reshape(-1), dtype=wp.int32, device=device)
            self.query_failures_d = wp.zeros(1, dtype=wp.int32, device=device)
            self.surface_ties_d = wp.zeros(1, dtype=wp.int32, device=device)
            wp.launch(
                copy_live_mesh_points_kernel,
                dim=self.n_verts,
                inputs=[self.pos_global, wp.int32(self.node_off)],
                outputs=[self.points_d],
                device=device,
            )
            self.mesh = wp.Mesh(
                points=self.points_d,
                indices=self.indices_d,
                support_winding_number=True,
            )
        self._grid_shape: tuple[int, int, int] | None = None
        self.query_failure_mask_d: wp.array | None = None
        self.surface_tie_mask_d: wp.array | None = None

    def _ensure_grid_scratch(self, grid: object) -> None:
        """Allocate per-cell diagnostics once, before the first physical-time remap."""
        shape = tuple(int(n) for n in grid.shape)
        if self._grid_shape is None:
            with wp.ScopedDevice(self.device):
                self.query_failure_mask_d = wp.zeros(shape, dtype=wp.int32, device=self.device)
                self.surface_tie_mask_d = wp.zeros(shape, dtype=wp.int32, device=self.device)
            self._grid_shape = shape
        elif shape != self._grid_shape:
            raise ValueError("a live mesh provider cannot be rebound to a differently shaped field grid")

    def _refresh(self) -> None:
        """Refresh the query mesh from live mechanics state without leaving the CUDA device."""
        with wp.ScopedDevice(self.device):
            wp.launch(
                copy_live_mesh_points_kernel,
                dim=self.n_verts,
                inputs=[self.pos_global, wp.int32(self.node_off)],
                outputs=[self.points_d],
                device=self.device,
            )
            self.mesh.refit()

    @staticmethod
    def _query_geometry(grid: object) -> tuple[wp.vec3, wp.float32, wp.float32, wp.float32]:
        origin_np = np.asarray(grid.origin, dtype=np.float64)
        extent = (np.asarray(grid.shape, dtype=np.float64) - 1.0) * float(grid.dx)
        # The box diagonal is an exact upper bound on point-to-surface distance because the live surface is
        # required to stay inside the field box.
        max_dist = float(np.linalg.norm(extent))
        # Eight float32 rounding operations bound coordinate construction + closest-point subtraction/norm.
        surface_epsilon = 8.0 * float(np.finfo(np.float32).eps) * max_dist
        return (
            wp.vec3(float(origin_np[0]), float(origin_np[1]), float(origin_np[2])),
            wp.float32(grid.dx),
            wp.float32(max_dist),
            wp.float32(surface_epsilon),
        )

    def _classify_membrane(self, grid: object, mask: wp.array) -> None:
        if grid.device != self.device:
            raise ValueError("field grid and live membrane mesh must reside on the same CUDA device")
        self._ensure_grid_scratch(grid)
        self._refresh()
        origin, dx, max_dist, surface_epsilon = self._query_geometry(grid)
        with wp.ScopedDevice(self.device):
            wp.launch(
                classify_live_membrane_kernel,
                dim=grid.shape,
                inputs=[self.mesh.id, origin, dx, max_dist, surface_epsilon,
                        self.query_failures_d, self.surface_ties_d,
                        self.query_failure_mask_d, self.surface_tie_mask_d],
                outputs=[mask],
                device=self.device,
            )

    def _classify_nucleus(self, grid: object, mask: wp.array) -> None:
        if grid.device != self.device:
            raise ValueError("field grid and live nucleus mesh must reside on the same CUDA device")
        self._ensure_grid_scratch(grid)
        self._refresh()
        origin, dx, max_dist, surface_epsilon = self._query_geometry(grid)
        with wp.ScopedDevice(self.device):
            wp.launch(
                classify_live_nucleus_kernel,
                dim=grid.shape,
                inputs=[self.mesh.id, origin, dx, max_dist, surface_epsilon,
                        self.query_failures_d, self.surface_ties_d,
                        self.query_failure_mask_d, self.surface_tie_mask_d],
                outputs=[mask],
                device=self.device,
            )

    def reset_diagnostics(self) -> None:
        """Explicitly start a new gate/run diagnostic window; classifications never clear failures."""
        with wp.ScopedDevice(self.device):
            self.query_failures_d.zero_()
            self.surface_ties_d.zero_()


class LiveMeshMembraneMaskProvider(_LiveMeshGeometry):
    """Outer live mesh implementing only the ``MembraneMaskProvider`` protocol."""

    def __init__(self, pos_global: wp.array, node_off: int, n_verts: int, faces_local: np.ndarray,
                 reference_verts: np.ndarray | None = None) -> None:
        super().__init__(pos_global, node_off, n_verts, faces_local, reference_verts)

    def classify_membrane(self, grid: object, mask: wp.array) -> None:
        self._classify_membrane(grid, mask)


class LiveMeshNucleusMaskProvider(_LiveMeshGeometry):
    """Inner live mesh implementing only the ``NucleusMaskProvider`` protocol."""

    def __init__(self, pos_global: wp.array, node_off: int, n_verts: int, faces_local: np.ndarray,
                 reference_verts: np.ndarray | None = None) -> None:
        super().__init__(pos_global, node_off, n_verts, faces_local, reference_verts)

    def classify_nucleus(self, grid: object, mask: wp.array) -> None:
        self._classify_nucleus(grid, mask)
