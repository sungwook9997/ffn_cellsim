r"""Device-resident Biot poroelastic pore-pressure field (Warp) — the FEM→FEM+CFD going-forward fluid.

The 0-D volumetric poroelastic (validated, `network_warp.py` biphasic block) becomes a spatially-resolved
pore-pressure field ``p(x,t)`` on a device grid, solved by Biot/Terzaghi consolidation ∂p/∂t = c_v ∇²p
(Darcy pore-flow through the saturated fiber skeleton; ``c_v`` = the consolidation coefficient = the
Moeendarbary 2013 diffusivity D≈40–60 µm²/s). This is the LARGE-build device-resident field solver — the
physics is validated (``scripts/engine_reval.py`` P7.F.1 1-D Terzaghi, P7.F.2 3-D Green's function, P7.F.3
IBM coupling, P7.F.4 coupled two-way FSI); this module is the Warp-kernel version, parity-gated against
those numpy prototypes.

Units: FF is µm, pN, s. p in pN/µm² (= Pa); D in µm²/s; dx in µm.

Going-forward (LARGE build, PI-gated): couple to the fiber nodes via immersed-boundary spread/interpolate
(node divergence → grid source; ∇p → node Darcy force), and RETIRE the lumped 6πηR node drag in the same
commit (the η double-count trap — the single biggest correctness hazard, per the roadmap).

Sanity Gate (this module): the explicit stencil is stable iff ``dt ≤ dx²/(6D)`` (3-D CFL); no-flux
(Neumann) box conserves total pore mass ∫p dV; parity with the numpy Green's-function prototype (P7.F.2)
to <1e-10. Verified in ``tests/ff/test_biot_fluid_warp.py`` / engine_reval P7.F.5.
"""
from __future__ import annotations

import numpy as np
import warp as wp

wp.init()


@wp.kernel
def biot_diffusion_kernel(
    p: wp.array3d(dtype=wp.float64),
    p_new: wp.array3d(dtype=wp.float64),
    D: wp.float64,
    dt: wp.float64,
    inv_dx2: wp.float64,
    n: wp.int32,
):
    """One explicit backward-difference-free (forward-Euler) Biot diffusion step on a cubic grid.

    Interior: p_new = p + dt·D·∇²p (7-point Laplacian). Faces: no-flux (Neumann) → the Laplacian's
    outward term is dropped so the face value is held (∂p/∂n = 0). Stable for dt ≤ dx²/(6D)."""
    i, j, k = wp.tid()
    # no-flux box: hold the six faces (Neumann) — matches the numpy prototype's lap[face]=0
    if i == 0 or i == n - 1 or j == 0 or j == n - 1 or k == 0 or k == n - 1:
        p_new[i, j, k] = p[i, j, k]
        return
    lap = (p[i + 1, j, k] + p[i - 1, j, k]
           + p[i, j + 1, k] + p[i, j - 1, k]
           + p[i, j, k + 1] + p[i, j, k - 1]
           - wp.float64(6.0) * p[i, j, k]) * inv_dx2
    p_new[i, j, k] = p[i, j, k] + dt * D * lap


@wp.kernel
def biot_gradient_kernel(
    pgrid: wp.array3d(dtype=wp.vec3d),      # pressure packed in component 0
    gradp: wp.array3d(dtype=wp.vec3d),      # output ∇p [Pa/µm]
    inv_2dx: wp.float64,                    # 1/(2·dx)
    n: wp.int32,
):
    """Central-difference gradient ∇p of the scalar field packed in ``pgrid[...,0]`` → a vec3d ``gradp``
    field [Pa/µm], for the fiber (MT) −V_seg·∇p buoyancy force. Faces: one-sided → 0 (no gradient across the
    no-flux box), matching ``biot_diffusion_kernel``'s Neumann faces so ∇p is not spuriously large there."""
    i, j, k = wp.tid()
    if i == 0 or i == n - 1 or j == 0 or j == n - 1 or k == 0 or k == n - 1:
        gradp[i, j, k] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
        return
    gx = (pgrid[i + 1, j, k][0] - pgrid[i - 1, j, k][0]) * inv_2dx
    gy = (pgrid[i, j + 1, k][0] - pgrid[i, j - 1, k][0]) * inv_2dx
    gz = (pgrid[i, j, k + 1][0] - pgrid[i, j, k - 1][0]) * inv_2dx
    gradp[i, j, k] = wp.vec3d(gx, gy, gz)


class BiotField:
    """A device-resident 3-D Biot pore-pressure field with an explicit consolidation solver."""

    def __init__(self, n: int, dx: float, D: float, device: str = "cpu"):
        self.n = int(n)
        self.dx = float(dx)
        self.D = float(D)
        self.device = device
        self.dt_max = dx * dx / (6.0 * D)                    # 3-D explicit CFL
        self._p = wp.zeros((n, n, n), dtype=wp.float64, device=device)
        self._p2 = wp.zeros((n, n, n), dtype=wp.float64, device=device)

    def set_field(self, p0: np.ndarray) -> None:
        assert p0.shape == (self.n, self.n, self.n)
        self._p = wp.array(np.ascontiguousarray(p0, np.float64), dtype=wp.float64, device=self.device)
        self._p2 = wp.zeros_like(self._p)

    def get_field(self) -> np.ndarray:
        return self._p.numpy()

    def step(self, dt: float) -> None:
        """Advance one step (double-buffered). ``dt`` must satisfy the CFL ``dt ≤ dt_max``."""
        if dt > self.dt_max * (1.0 + 1e-9):
            raise ValueError(f"dt={dt} exceeds 3-D CFL dt_max={self.dt_max}")
        wp.launch(biot_diffusion_kernel, dim=(self.n, self.n, self.n),
                  inputs=[self._p, self._p2, wp.float64(self.D), wp.float64(dt),
                          wp.float64(1.0 / (self.dx * self.dx)), wp.int32(self.n)],
                  device=self.device)
        self._p, self._p2 = self._p2, self._p                # swap buffers

    def run(self, t_end: float, dt: float | None = None) -> np.ndarray:
        dt = dt if dt is not None else 0.2 * self.dt_max
        for _ in range(int(t_end / dt)):
            self.step(dt)
        return self.get_field()

    def pore_mass(self) -> float:
        """∫ p dV — conserved by the no-flux box (a mass-conservation invariant)."""
        return float(self._p.numpy().sum() * self.dx ** 3)


# ── Immersed-boundary two-way coupling (Peskin 4-point) — device kernels ────────────────────────────
@wp.func
def _peskin4(r: wp.float64) -> wp.float64:
    """Peskin 4-point regularized delta (1-D factor), argument r = (x_grid − X_node)/h."""
    a = wp.abs(r)
    if a <= wp.float64(1.0):
        return (wp.float64(3.0) - wp.float64(2.0) * a
                + wp.sqrt(wp.float64(1.0) + wp.float64(4.0) * a - wp.float64(4.0) * a * a)) / wp.float64(8.0)
    if a <= wp.float64(2.0):
        rad = wp.float64(-7.0) + wp.float64(12.0) * a - wp.float64(4.0) * a * a
        if rad < wp.float64(0.0):
            rad = wp.float64(0.0)
        return (wp.float64(5.0) - wp.float64(2.0) * a - wp.sqrt(rad)) / wp.float64(8.0)
    return wp.float64(0.0)


@wp.kernel
def ibm_spread_kernel(
    node_pos: wp.array(dtype=wp.vec3d),
    node_src: wp.array(dtype=wp.vec3d),
    grid_src: wp.array(dtype=wp.vec3d, ndim=3),
    origin: wp.vec3d,
    dx: wp.float64,
    n: wp.int32,
):
    """SPREAD each node's source to the grid via the 3-D Peskin δ_h over its 4³ support:
    grid_src(x) += node_src · δ_h(x−X) with δ_h = ∏ peskin4/h (units 1/h³). ∫ grid_src dV = node_src
    (partition of unity → momentum conservation)."""
    t = wp.tid()
    Xp = node_pos[t]
    Fp = node_src[t]
    bi = wp.int32(wp.floor((Xp[0] - origin[0]) / dx))
    bj = wp.int32(wp.floor((Xp[1] - origin[1]) / dx))
    bk = wp.int32(wp.floor((Xp[2] - origin[2]) / dx))
    inv = wp.float64(1.0) / dx
    for di in range(-1, 3):
        i = bi + di
        if i < 0 or i >= n:
            continue
        wx = _peskin4(((origin[0] + wp.float64(i) * dx) - Xp[0]) * inv) * inv
        for dj in range(-1, 3):
            j = bj + dj
            if j < 0 or j >= n:
                continue
            wy = _peskin4(((origin[1] + wp.float64(j) * dx) - Xp[1]) * inv) * inv
            for dk in range(-1, 3):
                k = bk + dk
                if k < 0 or k >= n:
                    continue
                wz = _peskin4(((origin[2] + wp.float64(k) * dx) - Xp[2]) * inv) * inv
                w = wx * wy * wz
                wp.atomic_add(grid_src, i, j, k, Fp * w)


@wp.kernel
def ibm_interp_kernel(
    node_pos: wp.array(dtype=wp.vec3d),
    grid_field: wp.array(dtype=wp.vec3d, ndim=3),
    node_val: wp.array(dtype=wp.vec3d),
    origin: wp.vec3d,
    dx: wp.float64,
    n: wp.int32,
):
    """INTERPOLATE the grid field to each node via the same δ_h: node_val = Σ grid·δ_h·h³. The adjoint
    of spread (Peskin) → the coupling conserves momentum (no spurious force)."""
    t = wp.tid()
    Xp = node_pos[t]
    bi = wp.int32(wp.floor((Xp[0] - origin[0]) / dx))
    bj = wp.int32(wp.floor((Xp[1] - origin[1]) / dx))
    bk = wp.int32(wp.floor((Xp[2] - origin[2]) / dx))
    inv = wp.float64(1.0) / dx
    h3 = dx * dx * dx
    acc = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
    for di in range(-1, 3):
        i = bi + di
        if i < 0 or i >= n:
            continue
        wx = _peskin4(((origin[0] + wp.float64(i) * dx) - Xp[0]) * inv) * inv
        for dj in range(-1, 3):
            j = bj + dj
            if j < 0 or j >= n:
                continue
            wy = _peskin4(((origin[1] + wp.float64(j) * dx) - Xp[1]) * inv) * inv
            for dk in range(-1, 3):
                k = bk + dk
                if k < 0 or k >= n:
                    continue
                wz = _peskin4(((origin[2] + wp.float64(k) * dx) - Xp[2]) * inv) * inv
                acc = acc + grid_field[i, j, k] * (wx * wy * wz * h3)
    node_val[t] = acc


def ibm_spread(node_pos: np.ndarray, node_src: np.ndarray, n: int, dx: float, origin,
               device: str = "cpu") -> np.ndarray:
    """Spread node vector sources → a (n,n,n,3) grid source array (numpy convenience)."""
    npd = wp.array(np.ascontiguousarray(node_pos, np.float64), dtype=wp.vec3d, device=device)
    nsd = wp.array(np.ascontiguousarray(node_src, np.float64), dtype=wp.vec3d, device=device)
    grid = wp.zeros((n, n, n), dtype=wp.vec3d, device=device)
    wp.launch(ibm_spread_kernel, dim=node_pos.shape[0],
              inputs=[npd, nsd, grid, wp.vec3d(*origin), wp.float64(dx), wp.int32(n)], device=device)
    return grid.numpy()


def ibm_interp(node_pos: np.ndarray, grid_field: np.ndarray, dx: float, origin,
               device: str = "cpu") -> np.ndarray:
    """Interpolate a (n,n,n,3) grid field → node values (numpy convenience)."""
    n = grid_field.shape[0]
    npd = wp.array(np.ascontiguousarray(node_pos, np.float64), dtype=wp.vec3d, device=device)
    gfd = wp.array(np.ascontiguousarray(grid_field, np.float64), dtype=wp.vec3d, device=device)
    outd = wp.zeros(node_pos.shape[0], dtype=wp.vec3d, device=device)
    wp.launch(ibm_interp_kernel, dim=node_pos.shape[0],
              inputs=[npd, gfd, outd, wp.vec3d(*origin), wp.float64(dx), wp.int32(n)], device=device)
    return outd.numpy()
