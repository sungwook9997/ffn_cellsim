"""Discrete conservative RAD REFERENCE for the I1c G-actin monomer transport.

Host-side acceptance oracle (pure NumPy, NO Warp/simulation runtime) — the DISCRETE companion to
``transport_analytic``. It implements the exact conservative reaction-advection-diffusion stencil the
Warp-CUDA kernel (``transport.py``) ports to the device, so the HARD invariant a closed form cannot
express — total actin ``integral phi c dV + N_polymer`` conserved to machine precision — is gated on the
dev Mac. It is NOT a runtime and never a CPU simulation fallback (I0-A).

Conserved density ``u = phi c`` (monomer per bulk volume). With uniform porosity ``phi`` the flux is

    F = u v_f  -  D_c grad u          (advective at v_f + Fickian; phi D_c grad c = D_c grad u)

Conservative FV: an UPWIND advective face flux + a central diffusive face flux, each computed once per
face and scatter-subtracted into both cells (interior fluxes telescope). Non-fluid/edge faces are no-flux
(the membrane is impermeable to monomer — ENGINE_ARCHITECTURE_PLAN §②). The reaction ``R`` EXCHANGES with a
scalar bound-polymer pool: every monomer removed from the field (barbed consumption) is added to
``bound_monomer`` and vice versa (pointed release), so ``sum_i u_i V_i + bound_monomer`` is invariant to
round-off regardless of advection/diffusion/reaction.

Sanity Gate (self-tested in tests/ac/fluid/test_i1c_transport_oracle.py):
  * total actin conserved to round-off (closed domain, arbitrary v_f + reaction);
  * positivity: a non-negative field stays non-negative under the CFL;
  * pure diffusion (v_f=0, R=0) recovers a cosine bleach at rate D_c|k|^2 (FRAP reads D_c);
  * pure advection (D_c=0, R=0) transports a blob at v_f (centroid speed == v_f, != q);
  * uniform-state preservation under (divergence-free) motion.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

__all__ = ["rad_cfl_dt", "RADTransportReference"]


def rad_cfl_dt(dx: float, d_c: float, v_max: float, dim: int, safety: float = 0.5) -> float:
    """Combined advection+diffusion explicit CFL ``safety / (2 d D_c/dx^2 + v_max/dx)``.

    Args:
        dx: Grid spacing.
        d_c: Monomer diffusivity.
        v_max: Max |v_f| component (advective speed).
        dim: Spatial dimension.
        safety: Fraction of the stability limit (< 1).

    Returns:
        Stable explicit step for upwind advection + central diffusion.
    """
    diff_rate = 2.0 * dim * d_c / (dx * dx)
    adv_rate = v_max / dx
    denom = diff_rate + adv_rate
    if denom <= 0.0:
        raise ValueError("need D_c > 0 or v_max > 0 for a finite CFL")
    return safety / denom


@dataclass
class RADTransportReference:
    """Explicit conservative reaction-advection-diffusion stepper for ``u = phi c`` (NumPy oracle).

    Attributes:
        shape: Grid shape.
        dx: Uniform spacing.
        phi: Porosity (uniform; I0-B1b GAP — swept as an oracle parameter).
        d_c: Monomer diffusivity (I0-B1c draft ~2-6 um^2/s).
        mask: FLUID/OUTSIDE/NUCLEUS classification (defaults all-fluid); non-fluid faces are no-flux.
    """

    shape: tuple[int, ...]
    dx: float
    phi: float
    d_c: float
    mask: npt.NDArray[np.int_] | None = None
    bound_monomer: float = 0.0
    _FLUID: int = field(default=1, init=False, repr=False)

    def __post_init__(self) -> None:
        if not (0.0 < self.phi <= 1.0):
            raise ValueError("phi must be in (0, 1]")
        if self.mask is None:
            self.mask = np.full(self.shape, self._FLUID, dtype=np.int_)
        self.dim = len(self.shape)

    def concentration(self, u: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Free-monomer concentration ``c = u / phi``."""
        return u / self.phi

    def field_content(self, u: npt.NDArray[np.float64]) -> float:
        """Field monomer content ``sum_i u_i V_i`` over the fluid cells."""
        return float(np.sum(u[self.mask == self._FLUID]) * self.dx**self.dim)

    def total_actin(self, u: npt.NDArray[np.float64]) -> float:
        """Total actin ``sum_i u_i V_i + bound_monomer`` (the conserved quantity)."""
        return self.field_content(u) + self.bound_monomer

    def _flux_divergence(
        self, u: npt.NDArray[np.float64], v_f: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Conservative div(F), F = u v_f - D_c grad u, with upwind advection + no-flux masked faces."""
        fluid = self.mask == self._FLUID
        div = np.zeros_like(u)
        for axis in range(self.dim):
            n = self.shape[axis]
            lo = np.take(fluid, np.arange(n - 1), axis=axis)
            hi = np.take(fluid, np.arange(1, n), axis=axis)
            both = lo & hi  # only fluid-fluid faces carry flux

            u_lo = np.take(u, np.arange(n - 1), axis=axis)
            u_hi = np.take(u, np.arange(1, n), axis=axis)
            vlo = np.take(v_f[..., axis], np.arange(n - 1), axis=axis)
            vhi = np.take(v_f[..., axis], np.arange(1, n), axis=axis)
            w = 0.5 * (vlo + vhi)  # face-normal velocity (this axis)

            u_up = np.where(w > 0.0, u_lo, u_hi)  # upwind density
            f_adv = w * u_up
            f_diff = -self.d_c * (u_hi - u_lo) / self.dx
            f = np.where(both, f_adv + f_diff, 0.0)

            idx_lo = [slice(None)] * self.dim
            idx_hi = [slice(None)] * self.dim
            idx_lo[axis] = slice(0, n - 1)
            idx_hi[axis] = slice(1, n)
            div[tuple(idx_lo)] += f / self.dx
            div[tuple(idx_hi)] -= f / self.dx
        div[~fluid] = 0.0
        return div

    def step(
        self,
        u: npt.NDArray[np.float64],
        dt: float,
        v_f: npt.NDArray[np.float64] | None = None,
        reaction: npt.NDArray[np.float64] | None = None,
    ) -> npt.NDArray[np.float64]:
        """One explicit conservative RAD step; updates ``bound_monomer`` in place, returns the new ``u``.

        Args:
            u: Current density field ``phi c``.
            dt: Time step (<= :func:`rad_cfl_dt`).
            v_f: Cell-centred pore-fluid velocity, shape ``(*grid, dim)``; ``None`` = quiescent.
            reaction: Net field source ``R`` per cell (release positive, consumption negative). Its
                integral is transferred to/from ``bound_monomer`` so total actin is conserved.

        Returns:
            The updated density field.
        """
        if v_f is None:
            v_f = np.zeros((*self.shape, self.dim))
        div = self._flux_divergence(u, v_f)
        r = np.zeros_like(u) if reaction is None else np.asarray(reaction, dtype=np.float64)
        u_new = u + dt * (-div + r)
        u_new = np.where(self.mask == self._FLUID, u_new, u)
        # monomers added to the field come OUT of the bound pool (and vice versa): exact exchange.
        self.bound_monomer -= float(np.sum(r[self.mask == self._FLUID]) * self.dx**self.dim * dt)
        return u_new
