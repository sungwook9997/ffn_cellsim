"""Conservative Biot p/mass substrate (Warp-CUDA kernel SOURCE) + PressureCoupling force primitive — I1a.

GPU-resident implementation of the fluid-content balance on the masked moving cell domain:

    S dp/dt + alpha div(v_s) + div(q) = s_water ,   q = -(k/mu) grad p = -mobility grad p
=>  S dp/dt = mobility * laplacian_masked(p) - alpha div(v_s) + s_water            (p = p_ext+p_bar+p_excess)

The explicit update kernel is the GATHER form of the conservative finite-volume stencil the pure-NumPy
oracle ``fv_reference.darcy_divergence`` verifies on the dev Mac: each fluid-fluid face contributes
``(p_nbr - p_c)`` symmetrically to both cells (discretely conservative without atomics); OUTSIDE /
NUCLEUS / grid-edge faces are no-flux; the membrane hydraulic-flux BC is folded into ``s_water`` by
``boundary.MembraneFluxBC``. Runtime is CUDA-only (I0-A); construction requires the FieldGrid's CUDA
device. The native gate runs constant-state, impermeable, source/dilatation, live-surface-flux, Green-spread,
and Terzaghi drained-slab checks on CUDA.  The drained-slab path uses the same conservative cell-centred
operator with the standard half-cell Dirichlet ghost contribution exposed below.

``PressureCoupling`` is the §1.4 inner-force-assembly primitive: ``accumulate(state, out_force)`` adds the
fluid->solid body force ``-alpha grad(p)`` (the ``-alpha p I`` traction) at each solid node, interpolated
from the grid with the Peskin stencil (the transpose-adjoint of the spread — gated by ``ibm_reference``).
The lead-owned integrator sums it with ``MyosinForce`` (I3) and ``StericForce`` (I2b).
"""

from __future__ import annotations

import warp as wp  # noqa: E402  (Warp-CUDA runtime; HOOMD never imported — contract test enforced)

from aleph.components.fluid.field_grid import FieldGrid

__all__ = [
    "biot_pmass_update_kernel",
    "biot_pmass_dirichlet_x_update_kernel",
    "BiotSubstrate",
    "PressureCoupling",
]

_FLUID = wp.constant(1)  # must match fv_reference.FLUID / field_grid mask codes


@wp.kernel
def biot_pmass_update_kernel(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    s_water: wp.array3d(dtype=wp.float64),
    div_vs: wp.array3d(dtype=wp.float64),
    mobility: wp.float64,
    dx: wp.float64,
    storage_S: wp.float64,
    alpha: wp.float64,
    dt: wp.float64,
    p_new: wp.array3d(dtype=wp.float64),
) -> None:
    """One explicit conservative p/mass step (gather form of ``fv_reference.darcy_divergence``).

    Non-fluid cells are held (they are not DOFs). Only FLUID|FLUID faces carry Darcy flux; OUTSIDE
    (membrane) and NUCLEUS faces are no-flux at the kernel level — the membrane hydraulic flux enters
    through ``s_water``, the nucleus is relative no-flux.
    """
    i, j, k = wp.tid()
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]
    if mask[i, j, k] != _FLUID:
        p_new[i, j, k] = p[i, j, k]
        return
    pc = p[i, j, k]
    acc = wp.float64(0.0)  # sum over fluid-fluid faces of (p_nbr - pc)
    if i > 0 and mask[i - 1, j, k] == _FLUID:
        acc += p[i - 1, j, k] - pc
    if i < nx - 1 and mask[i + 1, j, k] == _FLUID:
        acc += p[i + 1, j, k] - pc
    if j > 0 and mask[i, j - 1, k] == _FLUID:
        acc += p[i, j - 1, k] - pc
    if j < ny - 1 and mask[i, j + 1, k] == _FLUID:
        acc += p[i, j + 1, k] - pc
    if k > 0 and mask[i, j, k - 1] == _FLUID:
        acc += p[i, j, k - 1] - pc
    if k < nz - 1 and mask[i, j, k + 1] == _FLUID:
        acc += p[i, j, k + 1] - pc
    inv_dx2 = wp.float64(1.0) / (dx * dx)
    neg_div_q = mobility * inv_dx2 * acc  # -div(q) over fluid faces
    rate = (neg_div_q - alpha * div_vs[i, j, k] + s_water[i, j, k]) / storage_S
    p_new[i, j, k] = pc + dt * rate


@wp.kernel
def biot_pmass_dirichlet_x_update_kernel(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    mobility: wp.float64,
    dx: wp.float64,
    storage_S: wp.float64,
    dt: wp.float64,
    p_low: wp.float64,
    p_high: wp.float64,
    p_new: wp.array3d(dtype=wp.float64),
) -> None:
    """Conservative diffusion with drained Dirichlet values on the low/high x faces.

    The edge-cell centre is ``dx/2`` from its boundary.  A ghost value ``2*p_bc-p_edge`` therefore contributes
    ``2*(p_bc-p_edge)`` to the dimensionless face sum, exactly matching
    :class:`fv_reference.BiotFVReference` and the classic doubly-drained Terzaghi slab.
    """
    i, j, k = wp.tid()
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]
    if mask[i, j, k] != _FLUID:
        p_new[i, j, k] = p[i, j, k]
        return
    pc = p[i, j, k]
    acc = wp.float64(0.0)
    if i > 0 and mask[i - 1, j, k] == _FLUID:
        acc += p[i - 1, j, k] - pc
    if i < nx - 1 and mask[i + 1, j, k] == _FLUID:
        acc += p[i + 1, j, k] - pc
    if j > 0 and mask[i, j - 1, k] == _FLUID:
        acc += p[i, j - 1, k] - pc
    if j < ny - 1 and mask[i, j + 1, k] == _FLUID:
        acc += p[i, j + 1, k] - pc
    if k > 0 and mask[i, j, k - 1] == _FLUID:
        acc += p[i, j, k - 1] - pc
    if k < nz - 1 and mask[i, j, k + 1] == _FLUID:
        acc += p[i, j, k + 1] - pc
    if i == 0:
        acc += wp.float64(2.0) * (p_low - pc)
    if i == nx - 1:
        acc += wp.float64(2.0) * (p_high - pc)
    p_new[i, j, k] = pc + dt * mobility * acc / (storage_S * dx * dx)


@wp.func
def _peskin4(r: wp.float64) -> wp.float64:
    """4-point Peskin regularized delta (mirrors ``ff/biot_fluid_warp._peskin4`` and ``ibm_reference.peskin4``)."""
    a = wp.abs(r)
    if a <= wp.float64(1.0):
        return (wp.float64(3.0) - wp.float64(2.0) * a + wp.sqrt(wp.float64(1.0) + wp.float64(4.0) * a - wp.float64(4.0) * a * a)) / wp.float64(8.0)
    if a <= wp.float64(2.0):
        inner = wp.max(wp.float64(0.0), wp.float64(-7.0) + wp.float64(12.0) * a - wp.float64(4.0) * a * a)
        return (wp.float64(5.0) - wp.float64(2.0) * a - wp.sqrt(inner)) / wp.float64(8.0)
    return wp.float64(0.0)


@wp.kernel
def _pressure_gradient_masked_kernel(
    p: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    dx: wp.float64,
    gradp: wp.array3d(dtype=wp.vec3d),
) -> None:
    """Central/one-sided ``grad(p)`` on fluid cells, reproducing affine fields at a masked boundary."""
    i, j, k = wp.tid()
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]
    if mask[i, j, k] != _FLUID:
        gradp[i, j, k] = wp.vec3d(0.0, 0.0, 0.0)
        return
    pc = p[i, j, k]
    invdx = wp.float64(1.0) / dx
    gx = wp.float64(0.0)
    gy = wp.float64(0.0)
    gz = wp.float64(0.0)
    has_xm = i > 0 and mask[i - 1, j, k] == _FLUID
    has_xp = i < nx - 1 and mask[i + 1, j, k] == _FLUID
    has_ym = j > 0 and mask[i, j - 1, k] == _FLUID
    has_yp = j < ny - 1 and mask[i, j + 1, k] == _FLUID
    has_zm = k > 0 and mask[i, j, k - 1] == _FLUID
    has_zp = k < nz - 1 and mask[i, j, k + 1] == _FLUID
    if has_xm and has_xp and has_ym and has_yp and has_zm and has_zp:
        gx = (p[i + 1, j, k] - p[i - 1, j, k]) * (wp.float64(0.5) * invdx)
        gy = (p[i, j + 1, k] - p[i, j - 1, k]) * (wp.float64(0.5) * invdx)
        gz = (p[i, j, k + 1] - p[i, j, k - 1]) * (wp.float64(0.5) * invdx)
    else:
        # A component-wise one-sided difference is not affine-complete on voxel slivers where one axis has no
        # direct FLUID neighbour.  Fit the local 3-D affine polynomial instead.  This is the same mechanistic
        # reconstruction used for the membrane trace, evaluated only on boundary cells; full-interior cells
        # retain the six-read central fast path above.
        normal = wp.mat44d()
        rhs = wp.vec4d(0.0, 0.0, 0.0, 0.0)
        sample_count = wp.float64(0.0)
        for di in range(-1, 3):
            ii = i + di
            if ii < 0 or ii >= nx:
                continue
            for dj in range(-1, 3):
                jj = j + dj
                if jj < 0 or jj >= ny:
                    continue
                for dk in range(-1, 3):
                    kk = k + dk
                    if kk < 0 or kk >= nz or mask[ii, jj, kk] != _FLUID:
                        continue
                    basis = wp.vec4d(
                        wp.float64(1.0), wp.float64(di), wp.float64(dj), wp.float64(dk))
                    normal += wp.outer(basis, basis)
                    rhs += basis * p[ii, jj, kk]
                    sample_count += wp.float64(1.0)
        det = wp.determinant(normal)
        scale = wp.max(wp.float64(1.0), sample_count)
        scale2 = scale * scale
        determinant_floor = wp.float64(64.0) * wp.float64(2.220446049250313e-16) * scale2 * scale2
        if wp.abs(det) > determinant_floor:
            coefficients = wp.mul(wp.inverse(normal), rhs)
            gx = coefficients[1] * invdx
            gy = coefficients[2] * invdx
            gz = coefficients[3] * invdx
        else:
            # Rank loss can only occur in an under-resolved one-cell-thick fragment.  Keep the available
            # directional derivative and let the native affine gate reject any non-negligible occurrence.
            if has_xp:
                gx = (p[i + 1, j, k] - pc) * invdx
            elif has_xm:
                gx = (pc - p[i - 1, j, k]) * invdx
            if has_yp:
                gy = (p[i, j + 1, k] - pc) * invdx
            elif has_ym:
                gy = (pc - p[i, j - 1, k]) * invdx
            if has_zp:
                gz = (p[i, j, k + 1] - pc) * invdx
            elif has_zm:
                gz = (pc - p[i, j, k - 1]) * invdx
    gradp[i, j, k] = wp.vec3d(gx, gy, gz)


@wp.kernel
def _interp_grad_force_kernel(
    node_pos: wp.array(dtype=wp.vec3d),
    node_volume: wp.array(dtype=wp.float64),
    gradp: wp.array3d(dtype=wp.vec3d),
    mask: wp.array3d(dtype=wp.int32),
    origin: wp.vec3d,
    dx: wp.float64,
    alpha: wp.float64,
    out_force: wp.array(dtype=wp.vec3d),
) -> None:
    """Peskin-interp grad(p) to a node and add the fluid->solid force -alpha*V_node*grad(p).

    Uses the SAME FLUID-masked, per-node-normalized delta/stencil as the spread (transpose-adjoint), so the
    transfer conserves momentum even for cortex nodes whose support crosses the immersed membrane
    (``ibm_reference`` plus the native CUDA work gate certify the identity). Force sign realises the
    ``sigma_total = sigma_eff - alpha p I`` coupling on the solid balance; myosin stays on the solid and is
    NOT re-added as a Darcy force.
    """
    n = wp.tid()
    nx = gradp.shape[0]
    ny = gradp.shape[1]
    nz = gradp.shape[2]
    Xp = node_pos[n]
    bi = int(wp.floor((Xp[0] - origin[0]) / dx))
    bj = int(wp.floor((Xp[1] - origin[1]) / dx))
    bk = int(wp.floor((Xp[2] - origin[2]) / dx))
    acc = wp.vec3d(0.0, 0.0, 0.0)
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
                wz = _peskin4((origin[2] + wp.float64(kk) * dx - Xp[2]) / dx)
                if mask[ii, jj, kk] == _FLUID:
                    weight = wx * wy * wz
                    acc += gradp[ii, jj, kk] * weight
                    weight_sum += weight
    if weight_sum > wp.float64(0.0):
        wp.atomic_add(out_force, n, acc * (-alpha * node_volume[n] / weight_sum))


class BiotSubstrate:
    """Explicit conservative Biot p/mass solver on a :class:`FieldGrid` (CUDA).

    Args:
        grid: The device field grid (owns p / mask / s_water / div_vs).
        mobility: Darcy mobility ``k/mu`` (I0-B1 derived, consistent with ``c_v``).
        storage_S: Storativity ``S = 1/M`` (I0-B1 derived from the ``c_v`` anchor).
        alpha: Biot-Willis coupling (I0-B1 ratified 1.0).
    """

    def __init__(self, grid: FieldGrid, *, mobility: float, storage_S: float, alpha: float = 1.0) -> None:
        self.grid = grid
        self.mobility = float(mobility)
        self.storage_S = float(storage_S)
        self.alpha = float(alpha)

    @property
    def c_v(self) -> float:
        """Consolidation coefficient ``mobility / S`` implied by the closure."""
        return self.mobility / self.storage_S

    def cfl_dt(self, safety: float = 0.9) -> float:
        return self.grid.cfl_dt(self.mobility, self.storage_S, safety)

    def step(self, dt: float) -> None:
        """Advance the pressure field one explicit step (double-buffer swap; no host roundtrip)."""
        g = self.grid
        with wp.ScopedDevice(g.device):
            wp.launch(
                biot_pmass_update_kernel,
                dim=g.shape,
                inputs=[g.p, g.mask, g.compose_water_sources(), g.div_vs,
                        wp.float64(self.mobility), wp.float64(g.dx), wp.float64(self.storage_S),
                        wp.float64(self.alpha), wp.float64(dt)],
                outputs=[g.p_new],
            )
        g.p, g.p_new = g.p_new, g.p  # ping-pong

    def step_dirichlet_x(self, dt: float, p_low: float, p_high: float) -> None:
        """Advance with drained pressure values on both x faces (Terzaghi/native-oracle path)."""
        g = self.grid
        with wp.ScopedDevice(g.device):
            wp.launch(
                biot_pmass_dirichlet_x_update_kernel,
                dim=g.shape,
                inputs=[
                    g.p, g.mask, wp.float64(self.mobility), wp.float64(g.dx),
                    wp.float64(self.storage_S), wp.float64(dt), wp.float64(p_low), wp.float64(p_high),
                ],
                outputs=[g.p_new],
                device=g.device,
            )
        g.p, g.p_new = g.p_new, g.p


class PressureCoupling:
    """§1.4 inner-force primitive: the fluid->solid ``-alpha p I`` body force at each solid node.

    ``accumulate(state, out_force)`` computes grad(p) on the masked grid, Peskin-interpolates it to the
    node positions in ``state``, and ADDS ``-alpha V_node grad(p)`` into ``out_force`` (the global
    per-node force array the lead-owned integrator sums with the other primitives). It never overwrites.
    """

    def __init__(self, substrate: BiotSubstrate) -> None:
        self.substrate = substrate
        g = substrate.grid
        with wp.ScopedDevice(g.device):
            self._gradp = wp.zeros(g.shape, dtype=wp.vec3d)

    def accumulate(self, state: object, out_force: wp.array) -> None:
        """Add the fluid->solid pressure force to ``out_force``.

        Args:
            state: Carries ``node_pos`` (wp.array vec3d, um) and ``node_volume`` (wp.array float64, um^3),
                the solid nodes and their representative control volumes.
            out_force: Global per-node force accumulator (wp.array vec3d, pN); ADDED to, not overwritten.
        """
        g = self.substrate.grid
        origin = wp.vec3d(float(g.origin[0]), float(g.origin[1]), float(g.origin[2]))
        with wp.ScopedDevice(g.device):
            wp.launch(_pressure_gradient_masked_kernel, dim=g.shape,
                      inputs=[g.p, g.mask, wp.float64(g.dx)], outputs=[self._gradp])
            wp.launch(_interp_grad_force_kernel, dim=state.node_pos.shape[0],
                      inputs=[state.node_pos, state.node_volume, self._gradp, g.mask, origin,
                              wp.float64(g.dx), wp.float64(self.substrate.alpha)],
                      outputs=[out_force])
