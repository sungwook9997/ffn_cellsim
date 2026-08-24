"""Discrete conservative finite-volume REFERENCE for the I1a Biot p/mass solver.

Host-side acceptance oracle (pure NumPy, NO Warp, NO simulation runtime). This is the DISCRETE
companion to the closed-form oracles (``consolidation_analytic`` / ``greens_analytic`` /
``manufactured``): where those give the CONTINUUM ground truth, this implements the exact
conservative cell-centred finite-volume STENCIL that the Warp-CUDA kernel (``biot_substrate.py``)
ports to the device. It lets the discrete conservation identities a closed form cannot express —
constant-state preservation, impermeable-limit mass conservation, global content == integrated
boundary flux to machine precision, discrete negative-semidefiniteness (pressure-work sign), and
refinement convergence to Terzaghi/Green — be gated on the dev Mac where CUDA is unavailable.

It is NOT a runtime and NEVER a CPU simulation fallback (I0-A: Warp-CUDA is the only runtime). It
exists solely so the ALGORITHM is CPU-green before the lead opens its A5000 gate slot; the Warp
kernel is the authoritative implementation of the identical stencil.

Governing balance (I1a, reduced pressure form; solid velocity ``v_s`` frozen at the outer step):

    S dp/dt + alpha div(v_s) + div(q) = s_water ,   q = -(k/mu) grad p = -mobility grad p

so, moving Darcy into the divergence,

    S dp/dt = -div(q) - alpha div(v_s) + s_water
            = mobility laplacian(p) - alpha div(v_s) + s_water      (constant mobility, interior).

Conservative discretisation: cell-centred ``p`` on a uniform structured grid (spacing ``dx``, cell
volume ``dx**d``, face area ``dx**(d-1)``). Each internal Darcy face flux
``q_face = -mobility (p_{k+1}-p_k)/dx`` is computed ONCE and scatter-subtracted into its two adjacent
cells, so interior fluxes telescope exactly and the global content balance closes to round-off:

    d/dt sum_i S p_i V_i = -sum_{boundary faces} q.n A + sum_i (s_water - alpha div v_s)_i V_i .

Domain masking (the moving membrane-minus-nucleus cell): a face carries Darcy flux only between two
FLUID cells; a FLUID|OUTSIDE face carries the prescribed membrane hydraulic flux; a FLUID|NUCLEUS
face carries zero (relative no-flux). This is the CPU shadow of ``domain.py`` cut-cell classification.

Sanity Gate (self-tested in tests/ac/fluid/test_i1a_solver_gates.py):
  * constant-state: a uniform field with zero source and no solid motion does not drift (also under
    a masked/cut domain — the remap adds no spurious source);
  * impermeable limit (all no-flux faces, zero source): sum_i S p_i V_i is conserved to round-off;
  * content balance: d/dt(total content) equals the integrated boundary flux + interior source;
  * pressure-work sign: the no-flux, source-free diffusion operator is negative-semidefinite
    (d/dt of 0.5 S integral p^2 <= 0), i.e. the discrete Laplacian is SPD;
  * CFL: the explicit step is stable at dt <= S dx^2 / (2 d mobility) = dx^2 / (2 d c_v);
  * convergence: the explicit scheme reproduces Terzaghi U(T_v) and the Green kernel spread at 2nd
    order in dx.

References: Biot 1941 J Appl Phys 12:155; LeVeque 2002, Finite Volume Methods for Hyperbolic
Problems (conservative FV); Terzaghi 1943; Carslaw & Jaeger 1959. c_v anchor: Moeendarbary 2013
Nat Mater 12:253 (KB-3.B3.2).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

__all__ = [
    "FLUID",
    "OUTSIDE",
    "NUCLEUS",
    "cfl_dt",
    "darcy_divergence",
    "boundary_flux_integral",
    "BiotFVReference",
]

# Cell-classification codes for the cut-cell domain mask.
OUTSIDE = 0  # beyond the membrane (not part of the cell interior)
FLUID = 1    # membrane-minus-nucleus fluid-saturated interior
NUCLEUS = 2  # inside the nuclear-envelope inclusion (relative no-flux)


def cfl_dt(dx: float, mobility: float, storage_S: float, dim: int, safety: float = 0.9) -> float:
    """Stable explicit time step ``dt <= safety * S dx^2 / (2 d mobility)``.

    Since ``c_v = mobility / S`` this is ``safety * dx^2 / (2 d c_v)`` — the classic explicit
    diffusion CFL. The Warp kernel must register the SAME limit (the k_EV/kmax bookkeeping analogue).

    Args:
        dx: Uniform grid spacing (length).
        mobility: Darcy mobility ``k/mu`` (length^2 / (pressure*time)).
        storage_S: Storativity ``S = 1/M`` (per pressure).
        dim: Spatial dimension.
        safety: Fraction of the stability limit (< 1).

    Returns:
        The maximum stable explicit step (time).
    """
    if dx <= 0.0 or mobility <= 0.0 or storage_S <= 0.0:
        raise ValueError("dx, mobility, storage_S must be positive")
    return safety * storage_S * dx * dx / (2.0 * dim * mobility)


def _face_darcy_flux(p: npt.NDArray[np.float64], axis: int, dx: float, mobility: float) -> npt.NDArray[np.float64]:
    """Internal-face Darcy flux ``q = -mobility (p_{k+1}-p_k)/dx`` along one axis (per unit area).

    Returns an array whose length along ``axis`` is one less than ``p`` (one value per internal face).
    """
    dp = np.diff(p, axis=axis)
    return -mobility * dp / dx


def darcy_divergence(
    p: npt.NDArray[np.float64],
    mobility: float,
    dx: float,
    mask: npt.NDArray[np.int_] | None = None,
    membrane_flux: npt.NDArray[np.float64] | None = None,
) -> npt.NDArray[np.float64]:
    """Conservative discrete ``div(q)`` with ``q = -mobility grad p`` on the masked cut domain.

    Face rule (per internal face between cells ``k`` and ``k+1`` along each axis):
      * FLUID | FLUID     -> Darcy flux ``-mobility (p_{k+1}-p_k)/dx``;
      * FLUID | OUTSIDE   -> prescribed membrane flux (from ``membrane_flux``; 0 if not supplied);
      * FLUID | NUCLEUS   -> 0 (relative no-flux);
      * neither FLUID     -> irrelevant (0).
    Grid-edge faces are treated as no-flux here (a fixed box); Dirichlet drainage is applied by the
    caller via ghost cells (see :meth:`BiotFVReference.step`). Each face flux is scatter-subtracted
    into both neighbours, so interior fluxes cancel exactly (discrete conservation).

    Args:
        p: Cell-centred pressure field, shape ``grid_shape``.
        mobility: Darcy mobility ``k/mu``.
        dx: Grid spacing.
        mask: Cell classification (:data:`FLUID` / :data:`OUTSIDE` / :data:`NUCLEUS`); ``None`` = all fluid.
        membrane_flux: Optional signed flux (+ = out of the fluid cell, in +axis sense) prescribed on
            FLUID|OUTSIDE faces; must broadcast to each axis' internal-face shape (a scalar is fine).

    Returns:
        ``div(q)`` per cell (same shape as ``p``); zero in non-fluid cells.
    """
    if mask is None:
        mask = np.full(p.shape, FLUID, dtype=np.int_)
    fluid = mask == FLUID
    div = np.zeros_like(p, dtype=np.float64)

    for axis in range(p.ndim):
        q = _face_darcy_flux(p, axis, dx, mobility)  # internal faces along this axis

        lo = np.take(fluid, np.arange(p.shape[axis] - 1), axis=axis)      # cell k
        hi = np.take(fluid, np.arange(1, p.shape[axis]), axis=axis)       # cell k+1
        lo_code = np.take(mask, np.arange(p.shape[axis] - 1), axis=axis)
        hi_code = np.take(mask, np.arange(1, p.shape[axis]), axis=axis)

        both_fluid = lo & hi
        # FLUID|OUTSIDE membrane faces: prescribe the hydraulic flux (signed by outward normal).
        memb = np.zeros_like(q)
        if membrane_flux is not None:
            mf = np.broadcast_to(np.asarray(membrane_flux, dtype=np.float64), q.shape)
            out_hi = lo & (hi_code == OUTSIDE)   # fluid at k, outside at k+1 -> outward normal +axis
            out_lo = (lo_code == OUTSIDE) & hi   # outside at k, fluid at k+1 -> outward normal -axis
            memb = np.where(out_hi, mf, memb)
            memb = np.where(out_lo, -mf, memb)
        # nucleus faces and all other non-both-fluid faces carry zero Darcy flux.
        q_eff = np.where(both_fluid, q, memb)

        # scatter into the two neighbouring cells: cell k loses q, cell k+1 gains q.
        idx_lo = [slice(None)] * p.ndim
        idx_hi = [slice(None)] * p.ndim
        idx_lo[axis] = slice(0, p.shape[axis] - 1)
        idx_hi[axis] = slice(1, p.shape[axis])
        div[tuple(idx_lo)] += q_eff / dx
        div[tuple(idx_hi)] -= q_eff / dx

    div[~fluid] = 0.0
    return div


def boundary_flux_integral(
    p: npt.NDArray[np.float64],
    mobility: float,
    dx: float,
    mask: npt.NDArray[np.int_],
    membrane_flux: npt.NDArray[np.float64] | None = None,
) -> float:
    """Net Darcy content flux OUT of the fluid domain across all non-fluid-fluid faces.

    For a masked interior domain (no Dirichlet grid edges) this must equal ``sum_i div(q)_i V_i``
    (the discrete divergence theorem) — the ``content == integrated boundary flux`` gate.

    Args:
        p: Pressure field.
        mobility: Darcy mobility.
        dx: Grid spacing.
        mask: Cell classification.
        membrane_flux: Optional prescribed membrane face flux (see :func:`darcy_divergence`).

    Returns:
        Total outward content flux (content per unit time).
    """
    face_area = dx ** (p.ndim - 1)
    fluid = mask == FLUID
    total = 0.0
    for axis in range(p.ndim):
        n_face = p.shape[axis] - 1
        lo = np.take(fluid, np.arange(n_face), axis=axis)
        hi = np.take(fluid, np.arange(1, p.shape[axis]), axis=axis)
        lo_code = np.take(mask, np.arange(n_face), axis=axis)
        hi_code = np.take(mask, np.arange(1, p.shape[axis]), axis=axis)
        # membrane faces (FLUID|OUTSIDE) carry the prescribed outward flux mf regardless of which side
        # the fluid is on; nucleus faces carry 0; fluid-fluid faces are interior and telescope to 0.
        is_membrane = (lo & (hi_code == OUTSIDE)) | ((lo_code == OUTSIDE) & hi)
        if membrane_flux is not None:
            mf = np.broadcast_to(np.asarray(membrane_flux, dtype=np.float64), is_membrane.shape)
            total += float(np.sum(np.where(is_membrane, mf, 0.0)) * face_area)
    return total


@dataclass
class BiotFVReference:
    """Stateful explicit reference stepper for the I1a conservative p/mass balance.

    Cell-centred pressure on a uniform ``dim``-D grid. The step is the explicit forward-Euler update
    of ``S dp/dt = -div(q) - alpha div(v_s) + s_water`` with the conservative :func:`darcy_divergence`.
    Domain edges are homogeneous no-flux by default; per-side Dirichlet drainage (Terzaghi) is applied
    with ghost cells. This is an ACCEPTANCE ORACLE, not the runtime.

    Attributes:
        shape: Grid shape (cells per axis).
        dx: Uniform spacing.
        mobility: Darcy mobility ``k/mu``.
        storage_S: Storativity ``S = 1/M``.
        mask: Cell classification; defaults to all-fluid.
        dirichlet: Optional dict ``{(axis, side): value}`` (side 0 = low face, 1 = high face) giving a
            fixed pressure on that grid-edge face (drained boundary).
    """

    shape: tuple[int, ...]
    dx: float
    mobility: float
    storage_S: float
    mask: npt.NDArray[np.int_] | None = None
    dirichlet: dict[tuple[int, int], float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.mask is None:
            self.mask = np.full(self.shape, FLUID, dtype=np.int_)
        self.dim = len(self.shape)

    @property
    def c_v(self) -> float:
        """Consolidation coefficient ``c_v = mobility / S`` implied by the closure."""
        return self.mobility / self.storage_S

    def cfl_dt(self, safety: float = 0.9) -> float:
        """Stable explicit step for this grid (see :func:`cfl_dt`)."""
        return cfl_dt(self.dx, self.mobility, self.storage_S, self.dim, safety)

    def total_content(self, p: npt.NDArray[np.float64]) -> float:
        """Total fluid content ``sum_i S p_i V_i`` over the FLUID cells."""
        cell_vol = self.dx**self.dim
        return float(self.storage_S * np.sum(p[self.mask == FLUID]) * cell_vol)

    def _apply_dirichlet_divergence(self, p: npt.NDArray[np.float64], div: npt.NDArray[np.float64]) -> float:
        """Add drained-edge (Dirichlet) face fluxes to ``div`` in place; return their outward integral.

        Ghost-cell drainage with the face dx/2 outside the edge-cell centre and the ghost centre dx
        outside: the face value equals ``p_bc`` when ``p_ghost = 2 p_bc - p_edge``. The OUTWARD Darcy
        content flux is then symmetric on both sides, ``q_out = 2 mobility (p_edge - p_bc)/dx``, and the
        edge cell picks up ``div += q_out/dx`` — equivalent to the standard Dirichlet Laplacian stencil
        ``(p_nbr - 3 p_edge + 2 p_bc)/dx^2``.
        """
        face_area = self.dx ** (self.dim - 1)
        outward = 0.0
        for (axis, side), p_bc in self.dirichlet.items():
            idx = [slice(None)] * self.dim
            idx[axis] = 0 if side == 0 else self.shape[axis] - 1
            edge = p[tuple(idx)]
            q_out = self.mobility * 2.0 * (edge - p_bc) / self.dx  # outward-normal Darcy flux (per area)
            div[tuple(idx)] += q_out / self.dx
            outward += float(np.sum(q_out) * face_area)
        return outward

    def rhs(
        self,
        p: npt.NDArray[np.float64],
        solid_dilatation_rate: npt.NDArray[np.float64] | float = 0.0,
        s_water: npt.NDArray[np.float64] | float = 0.0,
        alpha: float = 1.0,
        membrane_flux: npt.NDArray[np.float64] | None = None,
    ) -> npt.NDArray[np.float64]:
        """Right-hand side ``dp/dt`` of the explicit balance (per cell).

        Args:
            p: Current pressure field.
            solid_dilatation_rate: ``div(v_s)`` per cell (frozen outer-step solid source).
            s_water: Interior water source per cell (membrane BC is separate via ``membrane_flux``).
            alpha: Biot-Willis coupling (I0-B1: ratified 1.0; swept as an oracle here).
            membrane_flux: FLUID|OUTSIDE face flux (see :func:`darcy_divergence`).

        Returns:
            ``dp/dt`` per fluid cell; zero elsewhere.
        """
        div_q = darcy_divergence(p, self.mobility, self.dx, self.mask, membrane_flux)
        self._apply_dirichlet_divergence(p, div_q)
        rate = (-div_q - alpha * np.asarray(solid_dilatation_rate) + np.asarray(s_water)) / self.storage_S
        rate = np.where(self.mask == FLUID, rate, 0.0)
        return rate

    def step(
        self,
        p: npt.NDArray[np.float64],
        dt: float,
        **rhs_kwargs: object,
    ) -> npt.NDArray[np.float64]:
        """One explicit forward-Euler step; returns the updated field (does not mutate ``p``)."""
        return p + dt * self.rhs(p, **rhs_kwargs)  # type: ignore[arg-type]

    def diffusion_energy(self, p: npt.NDArray[np.float64]) -> float:
        """Discrete field energy ``0.5 S sum_i p_i^2 V_i`` (the Lyapunov functional).

        Under no-flux boundaries and zero source this must be non-increasing every step
        (pressure-work sign / SPD-Laplacian gate).
        """
        cell_vol = self.dx**self.dim
        return float(0.5 * self.storage_S * np.sum(p[self.mask == FLUID] ** 2) * cell_vol)
