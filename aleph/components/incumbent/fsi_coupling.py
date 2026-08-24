#!/usr/bin/env python3
r"""FSI div(v_s) -> fluid coupling INFRA (lead-owned, net-new) — activates the stubbed alpha*div(v_s) term.

The conservative Biot p/mass balance the fluid track solves is

    S dp/dt = mobility * lap(p) - alpha * div(v_s) + s_water ,

where ``div(v_s)`` is the SOLID skeleton's volumetric dilatation rate. ``ac.fluid.field_grid.FieldGrid``
allocates ``div_vs`` as ZEROS and the resting driver never fills it, so the fluid<-solid back-coupling is
STUBBED at 0: a deforming skeleton does NOT drive the pore fluid. This module is the missing wire.

It computes ``div(v_s)`` on the Eulerian field grid from the LAGRANGIAN solid node velocities using the same
regularized-delta (Peskin) machinery as ``PressureCoupling`` (the transpose-adjoint spread), then writes it
into ``grid.div_vs`` so the very next ``BiotSubstrate.step`` sees it. This is pure INFRA under ``ac/cell/``:
it imports the fluid ``FieldGrid`` read-only and edits NO track module and NO ``ff/`` file.

Method (standard immersed-boundary velocity interpolation, then Eulerian divergence):
  1. FLUID-mask-aware, per-node-normalized spread of the nodal solid velocity to a grid field:
     ``M = sum_n delta(x-X_n) v_s,n V_n`` and ``W = sum_n delta(x-X_n) V_n`` (atomics); normalization over
     each node's retained stencil preserves its full volume at the immersed boundary and is the exact
     transpose of ``PressureCoupling``; the grid solid velocity is ``V_s = M / W`` where W>0;
  2. masked central-difference divergence ``div(V_s)`` on FLUID cells (a face to a non-fluid/edge neighbour
     is one-sided => zero contribution, matching the I1a Darcy stencil), written to ``grid.div_vs``.

Sign/physics sanity: a skeleton CONTRACTING (nodes moving toward the centroid) has ``div(v_s) < 0``, so
``-alpha*div(v_s) > 0`` is a POSITIVE pressure source — squeezing the skeleton pressurises the pore fluid.
A rigid translation has ``div(v_s) = 0`` (no spurious source). Dormant particles (``active < 0``, e.g. myosin
backbone beads) are excluded so only the load-bearing actin skeleton drives the fluid.

Runtime: Warp-CUDA only (I0-A). Runs on the gbook A5000.
"""

from __future__ import annotations

import warp as wp

from aleph.components.fluid.field_grid import FieldGrid

__all__ = [
    "SolidDilatationCoupling",
    "displacement_rate_kernel",
    "spread_vs_kernel",
    "normalize_vs_kernel",
    "divergence_kernel",
]

_FLUID = wp.constant(1)


@wp.kernel
def displacement_rate_kernel(
    before: wp.array(dtype=wp.vec3d),
    after: wp.array(dtype=wp.vec3d),
    inv_dt_phys: wp.float64,
    accepted: wp.array(dtype=wp.int32),
    velocity: wp.array(dtype=wp.vec3d),
) -> None:
    """Convert one accepted outer-step displacement to solid velocity [um/s] on device."""
    i = wp.tid()
    if accepted[0] == 1:
        velocity[i] = (after[i] - before[i]) * inv_dt_phys
    else:
        velocity[i] = wp.vec3d(0.0, 0.0, 0.0)


@wp.func
def _peskin4(r: wp.float64) -> wp.float64:
    """4-point Peskin regularized delta (identical to ``biot_substrate._peskin4`` / ``ibm_reference``)."""
    a = wp.abs(r)
    if a <= wp.float64(1.0):
        return (wp.float64(3.0) - wp.float64(2.0) * a
                + wp.sqrt(wp.float64(1.0) + wp.float64(4.0) * a - wp.float64(4.0) * a * a)) / wp.float64(8.0)
    if a <= wp.float64(2.0):
        inner = wp.max(wp.float64(0.0), wp.float64(-7.0) + wp.float64(12.0) * a - wp.float64(4.0) * a * a)
        return (wp.float64(5.0) - wp.float64(2.0) * a - wp.sqrt(inner)) / wp.float64(8.0)
    return wp.float64(0.0)


@wp.kernel
def spread_vs_kernel(
    node_pos: wp.array(dtype=wp.vec3d),
    node_vs: wp.array(dtype=wp.vec3d),
    node_volume: wp.array(dtype=wp.float64),
    active: wp.array(dtype=wp.int32),
    mask: wp.array3d(dtype=wp.int32),
    origin: wp.vec3d,
    dx: wp.float64,
    momentum: wp.array3d(dtype=wp.vec3d),
    weight: wp.array3d(dtype=wp.float64),
) -> None:
    """Mask-normalized Peskin spread into ``momentum``/``weight`` for the mass-weighted grid velocity.

    Each node's retained FLUID weights sum to one, including when its stencil crosses the membrane.
    ``active[n] < 0`` (dormant myosin particles) are skipped so only the actin skeleton spreads.
    """
    n = wp.tid()
    if active[n] < 0:
        return
    nx = momentum.shape[0]
    ny = momentum.shape[1]
    nz = momentum.shape[2]
    Xp = node_pos[n]
    Vn = node_volume[n]
    vs = node_vs[n]
    bi = int(wp.floor((Xp[0] - origin[0]) / dx))
    bj = int(wp.floor((Xp[1] - origin[1]) / dx))
    bk = int(wp.floor((Xp[2] - origin[2]) / dx))
    weight_sum = wp.float64(0.0)
    for di in range(-1, 3):
        ii = bi + di
        if ii < 0 or ii >= nx:
            continue
        wx = _peskin4((origin[0] + wp.float64(ii) * dx - Xp[0]) / dx)
        for dj in range(-1, 3):
            jj = bj + dj
            if jj < 0 or jj >= ny:
                continue
            wy = _peskin4((origin[1] + wp.float64(jj) * dx - Xp[1]) / dx)
            for dk in range(-1, 3):
                kk = bk + dk
                if kk < 0 or kk >= nz:
                    continue
                if mask[ii, jj, kk] == _FLUID:
                    wz = _peskin4((origin[2] + wp.float64(kk) * dx - Xp[2]) / dx)
                    weight_sum += wx * wy * wz
    if weight_sum <= wp.float64(0.0):
        return
    for di in range(-1, 3):
        ii = bi + di
        if ii < 0 or ii >= nx:
            continue
        wx = _peskin4((origin[0] + wp.float64(ii) * dx - Xp[0]) / dx)
        for dj in range(-1, 3):
            jj = bj + dj
            if jj < 0 or jj >= ny:
                continue
            wy = _peskin4((origin[1] + wp.float64(jj) * dx - Xp[1]) / dx)
            for dk in range(-1, 3):
                kk = bk + dk
                if kk < 0 or kk >= nz:
                    continue
                if mask[ii, jj, kk] != _FLUID:
                    continue
                wz = _peskin4((origin[2] + wp.float64(kk) * dx - Xp[2]) / dx)
                wgt = wx * wy * wz / weight_sum
                wp.atomic_add(momentum, ii, jj, kk, vs * (wgt * Vn))
                wp.atomic_add(weight, ii, jj, kk, wgt * Vn)


@wp.kernel
def normalize_vs_kernel(
    momentum: wp.array3d(dtype=wp.vec3d),
    weight: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    w_floor: wp.float64,
    vs_grid: wp.array3d(dtype=wp.vec3d),
) -> None:
    """Grid solid velocity ``V_s = momentum / weight`` on cells with enough spread mass; else 0."""
    i, j, k = wp.tid()
    w = weight[i, j, k]
    if w > w_floor:
        vs_grid[i, j, k] = momentum[i, j, k] / w
    else:
        vs_grid[i, j, k] = wp.vec3d(0.0, 0.0, 0.0)


@wp.kernel
def divergence_kernel(
    vs_grid: wp.array3d(dtype=wp.vec3d),
    mask: wp.array3d(dtype=wp.int32),
    dx: wp.float64,
    div_out: wp.array3d(dtype=wp.float64),
) -> None:
    """Masked central-difference ``div(V_s)`` on FLUID cells (one-sided => 0 contribution at a non-fluid face)."""
    i, j, k = wp.tid()
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]
    if mask[i, j, k] != _FLUID:
        div_out[i, j, k] = wp.float64(0.0)
        return
    vc = vs_grid[i, j, k]
    xm = vc[0]
    if i > 0 and mask[i - 1, j, k] == _FLUID:
        xm = vs_grid[i - 1, j, k][0]
    xp = vc[0]
    if i < nx - 1 and mask[i + 1, j, k] == _FLUID:
        xp = vs_grid[i + 1, j, k][0]
    ym = vc[1]
    if j > 0 and mask[i, j - 1, k] == _FLUID:
        ym = vs_grid[i, j - 1, k][1]
    yp = vc[1]
    if j < ny - 1 and mask[i, j + 1, k] == _FLUID:
        yp = vs_grid[i, j + 1, k][1]
    zm = vc[2]
    if k > 0 and mask[i, j, k - 1] == _FLUID:
        zm = vs_grid[i, j, k - 1][2]
    zp = vc[2]
    if k < nz - 1 and mask[i, j, k + 1] == _FLUID:
        zp = vs_grid[i, j, k + 1][2]
    inv2dx = wp.float64(1.0) / (wp.float64(2.0) * dx)
    div_out[i, j, k] = (xp - xm) * inv2dx + (yp - ym) * inv2dx + (zp - zm) * inv2dx


class SolidDilatationCoupling:
    """Computes ``div(v_s)`` on the field grid from solid node motion and writes it into ``grid.div_vs``.

    Args:
        grid: the fluid ``FieldGrid`` (its ``div_vs`` array is the coupling target; read-only otherwise).
        w_floor: minimum spread weight before a cell's velocity is trusted (else V_s := 0), [um^3].
    """

    def __init__(self, grid: FieldGrid, *, w_floor: float = 0.0) -> None:
        self.grid = grid
        self.w_floor = float(w_floor)
        with wp.ScopedDevice(grid.device):
            self._momentum = wp.zeros(grid.shape, dtype=wp.vec3d)
            self._weight = wp.zeros(grid.shape, dtype=wp.float64)
            self._vs_grid = wp.zeros(grid.shape, dtype=wp.vec3d)

    def allocate_velocity(self, n_nodes: int) -> wp.array:
        """Allocate persistent solid-velocity scratch for one coupled mechanical state."""
        with wp.ScopedDevice(self.grid.device):
            return wp.zeros(int(n_nodes), dtype=wp.vec3d, device=self.grid.device)

    def update_from_displacement(
        self,
        before: wp.array,
        after: wp.array,
        dt_phys: float,
        node_volume: wp.array,
        active: wp.array,
        velocity: wp.array,
        accepted: wp.array,
    ) -> None:
        """Compute ``(after-before)/dt_phys`` and update ``grid.div_vs`` without a host roundtrip."""
        if dt_phys <= 0.0:
            raise ValueError("dt_phys must be positive")
        with wp.ScopedDevice(self.grid.device):
            wp.launch(
                displacement_rate_kernel,
                dim=after.shape[0],
                inputs=[before, after, wp.float64(1.0 / dt_phys), accepted],
                outputs=[velocity],
                device=self.grid.device,
            )
            self.update(after, velocity, node_volume, active)

    def update(self, node_pos: wp.array, node_vs: wp.array, node_volume: wp.array,
               active: wp.array) -> None:
        """Spread ``node_vs`` -> grid, take divergence, and store it in ``grid.div_vs`` (device-resident)."""
        g = self.grid
        origin = wp.vec3d(float(g.origin[0]), float(g.origin[1]), float(g.origin[2]))
        with wp.ScopedDevice(g.device):
            self._momentum.zero_()
            self._weight.zero_()
            wp.launch(spread_vs_kernel, dim=node_pos.shape[0],
                      inputs=[node_pos, node_vs, node_volume, active, g.mask, origin, wp.float64(g.dx)],
                      outputs=[self._momentum, self._weight])
            wp.launch(normalize_vs_kernel, dim=g.shape,
                      inputs=[self._momentum, self._weight, g.mask, wp.float64(self.w_floor)],
                      outputs=[self._vs_grid])
            wp.launch(divergence_kernel, dim=g.shape,
                      inputs=[self._vs_grid, g.mask, wp.float64(g.dx)], outputs=[g.div_vs])

    def div_vs_host(self):
        """Download the current div(v_s) field (out-of-hot-loop diagnostic)."""
        return self.grid.div_vs.numpy()
