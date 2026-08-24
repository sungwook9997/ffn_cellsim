"""Mixed q-p pore-fluid velocity (Warp-CUDA kernel SOURCE) — I1b (the fluid FLOWS).

GPU-resident reconstruction of the Darcy discharge and the absolute pore-fluid velocity on the same
masked moving domain as I1a:

    q   = -(k/mu) (grad p - rho_f b)          relative discharge (cell-centred, masked central diff)
    v_f = v_s + q / phi                        absolute pore-fluid velocity (I1c advects tracers with v_f)

Myosin enters only through the SOLID balance (it deforms the skeleton -> div(v_s) -> the fluid source in
I1a); it is NOT re-added here as a second direct Darcy body force (double-count guard, FSI_WIRING_DESIGN).
The only genuine fluid body force is ``rho_f b`` (gravity), physiologically tiny at cell scale but kept.

Runtime CUDA-only (I0-A). Host acceptance = ``darcy_analytic`` (flux linearity, ``v_f - v_s = q/phi``,
hydrostatic zero-discharge, radial conservation). ``phi`` is an I0-B1b GAP (params_i0b1b.yaml) — the
kinematic gate is phi-agnostic; the NATIVE run is blocked on a sourced ``phi``.
"""

from __future__ import annotations

import warp as wp  # noqa: E402  (Warp-CUDA runtime; HOOMD never imported)

from aleph.components.fluid.field_grid import FieldGrid

__all__ = ["darcy_discharge_kernel", "pore_velocity_kernel", "DarcyVelocity"]

_FLUID = wp.constant(1)


@wp.kernel
def darcy_discharge_kernel(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    mobility: wp.float64,
    dx: wp.float64,
    rho_f_b: wp.vec3d,
    q_out: wp.array3d(dtype=wp.vec3d),
) -> None:
    """Cell-centred Darcy discharge ``q = -mobility (grad p - rho_f b)`` on fluid cells (0 elsewhere).

    Masked central difference: a face to a non-fluid/edge neighbour is no-flux, so that axis falls back to
    a one-sided stencil against the centre (zero contribution) — consistent with the I1a update stencil.
    """
    i, j, k = wp.tid()
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]
    if mask[i, j, k] != _FLUID:
        q_out[i, j, k] = wp.vec3d(0.0, 0.0, 0.0)
        return
    pc = p[i, j, k]
    xm = pc
    if i > 0 and mask[i - 1, j, k] == _FLUID:
        xm = p[i - 1, j, k]
    xp = pc
    if i < nx - 1 and mask[i + 1, j, k] == _FLUID:
        xp = p[i + 1, j, k]
    ym = pc
    if j > 0 and mask[i, j - 1, k] == _FLUID:
        ym = p[i, j - 1, k]
    yp = pc
    if j < ny - 1 and mask[i, j + 1, k] == _FLUID:
        yp = p[i, j + 1, k]
    zm = pc
    if k > 0 and mask[i, j, k - 1] == _FLUID:
        zm = p[i, j, k - 1]
    zp = pc
    if k < nz - 1 and mask[i, j, k + 1] == _FLUID:
        zp = p[i, j, k + 1]
    inv2dx = wp.float64(1.0) / (wp.float64(2.0) * dx)
    grad = wp.vec3d((xp - xm) * inv2dx, (yp - ym) * inv2dx, (zp - zm) * inv2dx)
    q_out[i, j, k] = -mobility * (grad - rho_f_b)


@wp.kernel
def pore_velocity_kernel(
    v_s: wp.array3d(dtype=wp.vec3d),
    q: wp.array3d(dtype=wp.vec3d),
    mask: wp.array3d(dtype=wp.int32),
    phi: wp.float64,
    v_f_out: wp.array3d(dtype=wp.vec3d),
) -> None:
    """Absolute pore-fluid velocity ``v_f = v_s + q/phi`` (fluid cells only)."""
    i, j, k = wp.tid()
    if mask[i, j, k] != _FLUID:
        v_f_out[i, j, k] = wp.vec3d(0.0, 0.0, 0.0)
        return
    v_f_out[i, j, k] = v_s[i, j, k] + q[i, j, k] / phi


class DarcyVelocity:
    """Reconstructs q and v_f on a :class:`FieldGrid` (CUDA).

    Args:
        grid: The device field grid (reads ``p``, ``mask``).
        mobility: Darcy mobility ``k/mu`` (I0-B1).
        phi: Porosity (I0-B1b GAP — must be sourced before the native run).
        rho_f_b: Fluid body force ``rho_f b`` (Pa/um), default gravity-free.
    """

    def __init__(
        self,
        grid: FieldGrid,
        *,
        mobility: float,
        phi: float,
        rho_f_b: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> None:
        if not (0.0 < phi <= 1.0):
            raise ValueError("phi must be in (0, 1] (I0-B1b) — a sourced value is required for native runs")
        self.grid = grid
        self.mobility = float(mobility)
        self.phi = float(phi)
        self.rho_f_b = tuple(float(c) for c in rho_f_b)
        with wp.ScopedDevice(grid.device):
            self.q = wp.zeros(grid.shape, dtype=wp.vec3d)
            self.v_f = wp.zeros(grid.shape, dtype=wp.vec3d)

    def discharge(self) -> None:
        """Compute the Darcy discharge field ``q`` from the current pressure (device)."""
        g = self.grid
        with wp.ScopedDevice(g.device):
            wp.launch(
                darcy_discharge_kernel,
                dim=g.shape,
                inputs=[g.p, g.mask, wp.float64(self.mobility), wp.float64(g.dx), wp.vec3d(*self.rho_f_b)],
                outputs=[self.q],
            )

    def pore_velocity(self, v_s: wp.array) -> None:
        """Reconstruct ``v_f = v_s + q/phi`` from a solid-velocity grid field ``v_s`` (device)."""
        g = self.grid
        with wp.ScopedDevice(g.device):
            wp.launch(
                pore_velocity_kernel,
                dim=g.shape,
                inputs=[v_s, self.q, g.mask, wp.float64(self.phi)],
                outputs=[self.v_f],
            )
