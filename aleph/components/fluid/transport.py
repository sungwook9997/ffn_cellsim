"""Conservative G-actin monomer transport (Warp-CUDA kernel SOURCE) — I1c.

GPU-resident conservative reaction-advection-diffusion of the monomer density ``u = phi c`` on the SAME
masked moving cell domain as I1a/I1b:

    d(phi c)/dt + div(phi c v_f - phi D_c grad c) = R_polymer            (advect with v_f, NOT the discharge q)

The explicit update kernel is the GATHER form of the conservative UPWIND-advection + central-diffusion
stencil the pure-NumPy oracle ``transport_reference.RADTransportReference`` verifies on the dev Mac (each
face flux is a symmetric function of the two cells, so it is discretely conservative). Non-fluid/edge faces
are no-flux — the membrane is IMPERMEABLE to monomer (ENGINE_ARCHITECTURE_PLAN §②); water crosses in I1a,
monomer does not. This retires the fixed-scalar ``G_actin`` in ``ff/polymerization_warp.py`` (a spatially
uniform, infinite pool) in favour of a transported field (INTEGRATION.md patch-note).

Reaction ``R`` is IBM-coupled to the filaments: barbed-end consumption (``-k_on c delta`` at barbed-end
cells) SUBTRACTS from ``u`` and ADDS to the bound-polymer count; pointed-end/depoly release does the
reverse. Every polymer length change maps to an exact monomer count so ``integral phi c dV + N_polymer`` is
conserved to machine precision (the non-vacuous I1c invariant). Barbed/pointed events are spread/interp'd
with the SAME Peskin operator as the pressure coupling (adjoint, ``ibm_reference``).

Runtime CUDA-only (I0-A). Host acceptance = ``transport_analytic`` (FRAP D_c|k|^2, advection-follows-v_f,
treadmill) + ``transport_reference`` (conservation to round-off, positivity). ``D_c`` is an I0-B1c draft
(~2-6 um^2/s, provisional); ``c_0`` (initial MCF7 free-monomer pool) is a GAP -> the native run is blocked.
"""

from __future__ import annotations

import warp as wp  # noqa: E402  (Warp-CUDA runtime; HOOMD never imported)

from aleph.components.fluid.field_grid import FieldGrid

__all__ = ["rad_transport_kernel", "MonomerField"]

_FLUID = wp.constant(1)


@wp.func
def _upwind(w: wp.float64, u_lo: wp.float64, u_hi: wp.float64) -> wp.float64:
    """Upwind density for a face with normal velocity ``w`` (lo cell = -side, hi cell = +side)."""
    if w > wp.float64(0.0):
        return u_lo
    return u_hi


@wp.kernel
def rad_transport_kernel(
    u: wp.array3d(dtype=wp.float64),
    mask: wp.array3d(dtype=wp.int32),
    v_f: wp.array3d(dtype=wp.vec3d),
    reaction: wp.array3d(dtype=wp.float64),
    d_c: wp.float64,
    dx: wp.float64,
    dt: wp.float64,
    u_new: wp.array3d(dtype=wp.float64),
) -> None:
    """One explicit conservative RAD step for ``u = phi c`` (gather form of ``RADTransportReference``).

    Only FLUID|FLUID faces carry flux (upwind advection + central diffusion); OUTSIDE (membrane) and
    NUCLEUS faces are no-flux (monomer-impermeable). ``reaction`` is the per-cell net source R (release +,
    consumption -); the matching bound-pool update is applied by the IBM barbed/pointed events (host-free,
    exact-exchange), audited by the conservation gate.
    """
    i, j, k = wp.tid()
    nx = mask.shape[0]
    ny = mask.shape[1]
    nz = mask.shape[2]
    if mask[i, j, k] != _FLUID:
        u_new[i, j, k] = u[i, j, k]
        return
    uc = u[i, j, k]
    vc = v_f[i, j, k]
    inv_dx = wp.float64(1.0) / dx
    net_out = wp.float64(0.0)  # sum of outward face fluxes / dx  (= div F)

    # -x face
    if i > 0 and mask[i - 1, j, k] == _FLUID:
        w = wp.float64(0.5) * (v_f[i - 1, j, k][0] + vc[0])
        f = w * _upwind(w, u[i - 1, j, k], uc) - d_c * (uc - u[i - 1, j, k]) * inv_dx  # +x-oriented flux at the -x face
        net_out += -f * inv_dx  # outward (-x) = -(+x flux)
    # +x face
    if i < nx - 1 and mask[i + 1, j, k] == _FLUID:
        w = wp.float64(0.5) * (vc[0] + v_f[i + 1, j, k][0])
        f = w * _upwind(w, uc, u[i + 1, j, k]) - d_c * (u[i + 1, j, k] - uc) * inv_dx
        net_out += f * inv_dx
    # -y face
    if j > 0 and mask[i, j - 1, k] == _FLUID:
        w = wp.float64(0.5) * (v_f[i, j - 1, k][1] + vc[1])
        f = w * _upwind(w, u[i, j - 1, k], uc) - d_c * (uc - u[i, j - 1, k]) * inv_dx
        net_out += -f * inv_dx
    # +y face
    if j < ny - 1 and mask[i, j + 1, k] == _FLUID:
        w = wp.float64(0.5) * (vc[1] + v_f[i, j + 1, k][1])
        f = w * _upwind(w, uc, u[i, j + 1, k]) - d_c * (u[i, j + 1, k] - uc) * inv_dx
        net_out += f * inv_dx
    # -z face
    if k > 0 and mask[i, j, k - 1] == _FLUID:
        w = wp.float64(0.5) * (v_f[i, j, k - 1][2] + vc[2])
        f = w * _upwind(w, u[i, j, k - 1], uc) - d_c * (uc - u[i, j, k - 1]) * inv_dx
        net_out += -f * inv_dx
    # +z face
    if k < nz - 1 and mask[i, j, k + 1] == _FLUID:
        w = wp.float64(0.5) * (vc[2] + v_f[i, j, k + 1][2])
        f = w * _upwind(w, uc, u[i, j, k + 1]) - d_c * (u[i, j, k + 1] - uc) * inv_dx
        net_out += f * inv_dx

    u_new[i, j, k] = uc + dt * (-net_out + reaction[i, j, k])


class MonomerField:
    """Conservative G-actin monomer transport on a :class:`FieldGrid` (CUDA).

    Args:
        grid: The device field grid (shares mask + geometry with the Biot substrate).
        phi: Porosity (I0-B1b GAP; ``u = phi c``).
        d_c: Monomer diffusivity (I0-B1c draft ~2-6 um^2/s).
        c0: Initial uniform free-monomer concentration (I0-B1c GAP — MCF7-specific; required for native).
    """

    def __init__(self, grid: FieldGrid, *, phi: float, d_c: float, c0: float | None = None) -> None:
        if not (0.0 < phi <= 1.0):
            raise ValueError("phi must be in (0, 1] (I0-B1b)")
        self.grid = grid
        self.phi = float(phi)
        self.d_c = float(d_c)
        self.c0 = c0
        with wp.ScopedDevice(grid.device):
            init = 0.0 if c0 is None else float(phi * c0)
            self.u = wp.full(grid.shape, init, dtype=wp.float64)       # density u = phi c
            self.u_new = wp.zeros(grid.shape, dtype=wp.float64)
            self.reaction = wp.zeros(grid.shape, dtype=wp.float64)     # IBM-assembled barbed/pointed source
        self.bound_monomer = 0.0  # host mirror of the polymer-bound count (for the conservation gate)

    def cfl_dt(self, v_max: float, safety: float = 0.5) -> float:
        """Combined advection+diffusion CFL (matches ``transport_reference.rad_cfl_dt``)."""
        d = self.grid.dim
        return safety / (2.0 * d * self.d_c / self.grid.dx**2 + v_max / self.grid.dx)

    def step(self, dt: float, v_f: wp.array) -> None:
        """Advance one conservative RAD step advecting with ``v_f`` (device; double-buffer swap)."""
        g = self.grid
        with wp.ScopedDevice(g.device):
            wp.launch(
                rad_transport_kernel,
                dim=g.shape,
                inputs=[self.u, g.mask, v_f, self.reaction,
                        wp.float64(self.d_c), wp.float64(g.dx), wp.float64(dt)],
                outputs=[self.u_new],
            )
        self.u, self.u_new = self.u_new, self.u

    def total_actin(self) -> float:
        """Total actin ``sum u dV + bound_monomer`` (host readback — conservation gate, not the hot loop)."""
        u = self.u.numpy()
        m = self.grid.mask.numpy()
        return float(u[m == 1].sum() * self.grid.dx**self.grid.dim + self.bound_monomer)
