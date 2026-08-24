"""The live membrane-minus-nucleus fluid domain (Warp-CUDA) — OWNED by the fluid-spine track (§1.4).

``Domain`` classifies the conservative field grid into OUTSIDE / FLUID / NUCLEUS cells from the LIVE
outer membrane and the LIVE nuclear-envelope inclusion, and remaps fluid content conservatively as both
boundaries move. It exposes the FROZEN INTERFACE CONTRACT the nucleus track (Session B) codes against:

    Domain.set_nucleus_boundary(provider)

where ``provider`` yields the inner relative-no-flux mask. I1a ships a STATIC-sphere provider placeholder
(``StaticSphereNucleusMaskProvider``); the nucleus track ships a ``DeformableNucleusMaskProvider``
conforming to the SAME ``NucleusMaskProvider`` protocol (its moving oblate mesh is I1a's inner no-flux
boundary). Neither co-edits the other's file: the fluid track owns ``domain.py``, the nucleus track owns
its provider. A fixed cubic no-flux box is an analytic oracle only (``fv_reference``), never the cell.

Conservation under domain motion: a cell changing class transports its content across the moving boundary
(the ``S p*(L,t) L'(t)`` moving-face channel of ``manufactured.moving_interval_content_rate``); the remap
records exactly that transport so total content only changes by genuine boundary flux + moving-face
transport, never a numerical leak. Runtime CUDA-only (I0-A). Host acceptance = the analytic moving-interval
identity (already gated in ``test_manufactured_oracle``) + the native moving-boundary gate the lead runs.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

import warp as wp  # noqa: E402  (Warp-CUDA runtime; HOOMD never imported)

from aleph.components.fluid.field_grid import FieldGrid

__all__ = [
    "NucleusMaskProvider",
    "MembraneMaskProvider",
    "StaticSphereNucleusMaskProvider",
    "StaticSphereMembraneProvider",
    "Domain",
]

_FLUID = wp.constant(1)
_OUTSIDE = wp.constant(0)
_NUCLEUS = wp.constant(2)


@runtime_checkable
class MembraneMaskProvider(Protocol):
    """Classifies each grid cell as inside (FLUID) or outside (OUTSIDE) the live outer membrane."""

    def classify_membrane(self, grid: FieldGrid, mask: wp.array) -> None:
        """Write OUTSIDE where the cell centre is beyond the membrane, FLUID where inside (into ``mask``)."""


@runtime_checkable
class NucleusMaskProvider(Protocol):
    """§1.4 contract: mark cells enclosed by the live nuclear envelope as NUCLEUS (relative no-flux).

    Session B's ``DeformableNucleusMaskProvider`` conforms to this exact signature; I1a ships the static
    sphere below. Called AFTER the membrane classification — overwrites FLUID cells inside the nucleus.
    """

    def classify_nucleus(self, grid: FieldGrid, mask: wp.array) -> None:
        """Set ``mask`` to NUCLEUS where the (moving) nuclear envelope encloses the cell centre."""


@wp.kernel
def _sphere_membrane_kernel(
    origin: wp.vec3d, dx: wp.float64, centre: wp.vec3d, radius: wp.float64, mask: wp.array3d(dtype=wp.int32)
) -> None:
    """FLUID inside the sphere of ``radius`` about ``centre``, OUTSIDE beyond it."""
    i, j, k = wp.tid()
    x = origin[0] + wp.float64(i) * dx
    y = origin[1] + wp.float64(j) * dx
    z = origin[2] + wp.float64(k) * dx
    d = wp.vec3d(x, y, z) - centre
    if wp.length(d) <= radius:
        mask[i, j, k] = _FLUID
    else:
        mask[i, j, k] = _OUTSIDE


@wp.kernel
def _sphere_nucleus_kernel(
    origin: wp.vec3d, dx: wp.float64, centre: wp.vec3d, radius: wp.float64, mask: wp.array3d(dtype=wp.int32)
) -> None:
    """Overwrite FLUID cells inside the nucleus sphere with NUCLEUS (relative no-flux inclusion)."""
    i, j, k = wp.tid()
    if mask[i, j, k] != _FLUID:
        return
    x = origin[0] + wp.float64(i) * dx
    y = origin[1] + wp.float64(j) * dx
    z = origin[2] + wp.float64(k) * dx
    d = wp.vec3d(x, y, z) - centre
    if wp.length(d) <= radius:
        mask[i, j, k] = _NUCLEUS


@wp.kernel
def _remap_transport_kernel(
    p: wp.array3d(dtype=wp.float64),
    p_old: wp.array3d(dtype=wp.float64),
    mask_old: wp.array3d(dtype=wp.int32),
    mask_new: wp.array3d(dtype=wp.int32),
    p_bar: wp.array3d(dtype=wp.float64),
    storage_S: wp.float64,
    cell_vol: wp.float64,
    moving_face_content: wp.array(dtype=wp.float64),
) -> None:
    """Conservatively remap content as cells change class; accumulate the moving-face transport.

    * newly FLUID (membrane swept over it): seed ``p`` from the mean pressure of its FLUID neighbours (or
      ``p_bar`` if none) and record ``+S p V`` as content that entered across the moving membrane.  Every
      seed reads the immutable ``p_old`` and only neighbours that remain FLUID, so multi-cell boundary motion
      is deterministic and free of in-place read/write races;
    * newly non-FLUID (boundary retreated / nucleus advanced): record ``-S p V`` as content that left, and
      hold ``p`` (the cell is no longer a DOF).
    Cells that keep their class are untouched. The net recorded transport is exactly the moving-face
    channel of the analytic conservation identity, so total content changes only by real transport + flux.
    """
    i, j, k = wp.tid()
    nx = mask_new.shape[0]
    ny = mask_new.shape[1]
    nz = mask_new.shape[2]
    was_fluid = mask_old[i, j, k] == _FLUID
    now_fluid = mask_new[i, j, k] == _FLUID
    if was_fluid == now_fluid:
        return
    if now_fluid and not was_fluid:
        # Seed only from old-fluid neighbours that survive in the new domain.  Reading p_old makes this a
        # deterministic gather even when several adjacent cells become fluid in the same remap.
        s = wp.float64(0.0)
        cnt = 0
        if i > 0 and mask_old[i - 1, j, k] == _FLUID and mask_new[i - 1, j, k] == _FLUID:
            s += p_old[i - 1, j, k]
            cnt += 1
        if i < nx - 1 and mask_old[i + 1, j, k] == _FLUID and mask_new[i + 1, j, k] == _FLUID:
            s += p_old[i + 1, j, k]
            cnt += 1
        if j > 0 and mask_old[i, j - 1, k] == _FLUID and mask_new[i, j - 1, k] == _FLUID:
            s += p_old[i, j - 1, k]
            cnt += 1
        if j < ny - 1 and mask_old[i, j + 1, k] == _FLUID and mask_new[i, j + 1, k] == _FLUID:
            s += p_old[i, j + 1, k]
            cnt += 1
        if k > 0 and mask_old[i, j, k - 1] == _FLUID and mask_new[i, j, k - 1] == _FLUID:
            s += p_old[i, j, k - 1]
            cnt += 1
        if k < nz - 1 and mask_old[i, j, k + 1] == _FLUID and mask_new[i, j, k + 1] == _FLUID:
            s += p_old[i, j, k + 1]
            cnt += 1
        seeded = p_bar[i, j, k]
        if cnt > 0:
            seeded = s / wp.float64(cnt)
        p[i, j, k] = seeded
        wp.atomic_add(moving_face_content, 0, storage_S * seeded * cell_vol)
    else:
        # leaving the fluid domain
        wp.atomic_add(moving_face_content, 0, -storage_S * p_old[i, j, k] * cell_vol)


class StaticSphereMembraneProvider:
    """I1a placeholder outer membrane = a static sphere (the live Helfrich mesh replaces it at wiring)."""

    def __init__(self, centre: tuple[float, float, float], radius: float) -> None:
        self.centre = tuple(float(c) for c in centre)
        self.radius = float(radius)

    def classify_membrane(self, grid: FieldGrid, mask: wp.array) -> None:
        with wp.ScopedDevice(grid.device):
            wp.launch(
                _sphere_membrane_kernel,
                dim=grid.shape,
                inputs=[wp.vec3d(*grid.origin.tolist()), wp.float64(grid.dx),
                        wp.vec3d(*self.centre), wp.float64(self.radius)],
                outputs=[mask],
            )


class StaticSphereNucleusMaskProvider:
    """I1a placeholder nucleus = a static sphere (§1.4). Session B swaps in the deformable oblate mesh."""

    def __init__(self, centre: tuple[float, float, float], radius: float) -> None:
        self.centre = tuple(float(c) for c in centre)
        self.radius = float(radius)

    def classify_nucleus(self, grid: FieldGrid, mask: wp.array) -> None:
        with wp.ScopedDevice(grid.device):
            wp.launch(
                _sphere_nucleus_kernel,
                dim=grid.shape,
                inputs=[wp.vec3d(*grid.origin.tolist()), wp.float64(grid.dx),
                        wp.vec3d(*self.centre), wp.float64(self.radius)],
                outputs=[mask],
            )


class Domain:
    """The live membrane-minus-nucleus fluid domain on a :class:`FieldGrid` (CUDA).

    Args:
        grid: The device field grid it classifies.
        membrane: Provider for the outer membrane (I1a: static sphere; wiring: live Helfrich mesh).
        nucleus: Optional nucleus provider; set later via :meth:`set_nucleus_boundary` (§1.4).
    """

    def __init__(
        self, grid: FieldGrid, membrane: MembraneMaskProvider, nucleus: NucleusMaskProvider | None = None
    ) -> None:
        self.grid = grid
        self.membrane = membrane
        self._nucleus = nucleus
        with wp.ScopedDevice(grid.device):
            # Persistent scratch: physical-time remap performs only D2D copy + kernels, never allocation/readback.
            self._mask_old_d = wp.empty(grid.shape, dtype=wp.int32, device=grid.device)
            self._p_old_d = wp.empty(grid.shape, dtype=wp.float64, device=grid.device)
            self._moving_face_content_d = wp.zeros(1, dtype=wp.float64, device=grid.device)

    def set_nucleus_boundary(self, provider: NucleusMaskProvider) -> None:
        """§1.4 contract: install the inner relative-no-flux provider (static sphere or deformable mesh)."""
        if not isinstance(provider, NucleusMaskProvider):
            raise TypeError("provider must implement NucleusMaskProvider.classify_nucleus(grid, mask)")
        self._nucleus = provider

    def classify(self) -> None:
        """(Re)build the cell classification on device from the live membrane + nucleus boundaries."""
        self.membrane.classify_membrane(self.grid, self.grid.mask)
        if self._nucleus is not None:
            self._nucleus.classify_nucleus(self.grid, self.grid.mask)

    def remap(self) -> wp.array:
        """Reclassify after boundary motion and conservatively transport content on CUDA.

        Returns the persistent device scalar holding the net content that crossed the moving boundaries this
        outer step (the analytic ``moving_face`` channel).  Read it only after the physical-time loop, or clone
        it device-to-device when retaining per-step history.  No authoritative D2H read occurs here (NG-6).
        """
        g = self.grid
        with wp.ScopedDevice(g.device):
            wp.copy(self._mask_old_d, g.mask)
            wp.copy(self._p_old_d, g.p)
            self.classify()  # writes the NEW classification into g.mask
            self._moving_face_content_d.zero_()
            # NOTE: storage_S/p_bar are the substrate's; passed via attributes set by the scheduler.
            wp.launch(
                _remap_transport_kernel,
                dim=g.shape,
                inputs=[g.p, self._p_old_d, self._mask_old_d, g.mask, g.p_bar,
                        wp.float64(self._storage_S), wp.float64(g.dx**g.dim)],
                outputs=[self._moving_face_content_d],
            )
            return self._moving_face_content_d

    # the scheduler binds the storage closure so remap() can weight content correctly
    _storage_S: float = 1.0

    def bind_storage(self, storage_S: float) -> None:
        """Bind the storativity used to weight remapped content (called once by the scheduler)."""
        self._storage_S = float(storage_S)
