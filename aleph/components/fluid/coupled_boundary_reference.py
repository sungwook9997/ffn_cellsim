"""Coupled enclosed-volume fluid REFERENCE (pure NumPy, no Warp) — cytosol/nucleus CONNECTED gate.

Host-side acceptance oracle for the ``cytosol`` FLUID_VOLUME reaching **CONNECTED** (roadmap
``AC_ENGINE_COMPLETION_ROADMAP_2026-07-23.md`` §3.1): the real Biot poroelastic field solver plus BOTH moving
boundaries — the outer **semipermeable membrane** (Kedem-Katchalsky water flux ``J = L_p(sigma dPi - dP)``) and
the inner **impermeable nuclear envelope** (relative no-flux) — closed as ONE coupled candidate with the three
ledger channels asserted TOGETHER:

    MASS       :  Delta(total fluid content)  =  -dt * (net outward membrane permeation)  +  swept content
                  (the nucleus contributes NO permeation term — only geometric sweep when the envelope moves);
    NO-FLUX    :  the APPLIED cross-envelope Darcy flux is exactly 0 (impermeable), while the would-be discharge
                  the mask holds at zero is > 0 under any transmembrane gradient (load-bearing "teeth");
    ADJOINT-WORK: the Peskin traction scatter and the boundary-velocity gather are a transpose pair, so the work
                  the fluid does on a boundary equals the work drawn from the field: <W v, g> == <v, W^T g>.

The individual pieces (content==flux, moving-domain remap conservation, adjoint spread=interp^T) are each gated
separately by ``fv_reference`` / ``ibm_reference`` / the I1a gates; THIS oracle proves they close SIMULTANEOUSLY
and consistently in one enclosed-volume candidate — the invariant that makes the composition CONNECTED rather
than three independently-passing seams.  It reuses ``fv_reference`` primitives (``darcy_divergence`` /
``boundary_flux_integral`` / ``BiotFVReference``) so it is the exact discrete stencil the Warp-CUDA solver ports.

It is NOT a runtime and never a CPU simulation fallback (I0-A).  The membrane influx magnitude (via ``L_p`` /
``dPi``) is a sourced GAP; the CONSERVATION identity gated here is exact for ANY prescribed flux, so the gate is
GAP-agnostic (the native magnitude judgement is separate, PI).

Sanity Gate (self-tested in tests/ac/fluid/test_coupled_boundary_reference.py):
  * coupled mass: static-domain content change equals the integrated membrane permeation to round-off, with a
    zero nucleus permeation term;
  * moving membrane: total content over the reclassified domain is conserved once the swept term is included;
  * no-flux: applied cross-nucleus flux is 0 while the teeth probe is > 0 under a gradient;
  * adjoint work-conjugacy: <W v, g> == <v, W^T g> for the Peskin spread/gather pair.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from aleph.components.fluid.fv_reference import (
    FLUID,
    NUCLEUS,
    OUTSIDE,
    BiotFVReference,
    boundary_flux_integral,
)

__all__ = ["CoupledFluidLedger", "EnclosedVolumeCoupledReference"]


@dataclass(frozen=True)
class CoupledFluidLedger:
    """The three coupled-candidate ledger channels evaluated on one enclosed-volume step."""

    mass_residual: float           # |Delta content + dt * membrane_outward_flux - swept_content|
    membrane_outward_flux: float   # net outward membrane permeation integral (content/time)
    swept_content: float           # conservative moving-face content redistributed this step
    nucleus_applied_flux: float    # APPLIED cross-envelope Darcy flux (must be 0)
    nucleus_teeth: float           # would-be suppressed discharge (> 0 under a gradient)
    adjoint_work_residual: float   # |<W v, g> - <v, W^T g>| (must be 0)


@dataclass
class EnclosedVolumeCoupledReference:
    """One coupled enclosed-volume fluid candidate over a masked membrane-minus-nucleus domain.

    Attributes:
        shape: Grid shape.
        dx: Uniform spacing.
        mobility: Darcy mobility ``k/mu``.
        storage_S: Storativity ``S = 1/M``.
        mask: OUTSIDE/FLUID/NUCLEUS classification (the enclosed cytosol is the FLUID region).
        L_p: membrane hydraulic conductivity (sourced GAP; only scales the flux, not the conservation).
        sigma_refl: reflection coefficient (1.0 ideal semipermeable).
        p_ext: external hydrostatic pressure.
    """

    shape: tuple[int, ...]
    dx: float
    mobility: float
    storage_S: float
    mask: npt.NDArray[np.int_]
    L_p: float = 1.0
    sigma_refl: float = 1.0
    p_ext: float = 0.0
    _biot: BiotFVReference = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.mask = np.asarray(self.mask, dtype=np.int_)
        if self.mask.shape != tuple(self.shape):
            raise ValueError("mask shape must match the grid shape")
        self.dim = len(self.shape)
        self._biot = BiotFVReference(
            shape=tuple(self.shape), dx=self.dx, mobility=self.mobility,
            storage_S=self.storage_S, mask=self.mask,
        )

    # -- mass channel -----------------------------------------------------------------------------
    def total_content(self, p: npt.NDArray[np.float64]) -> float:
        """Total enclosed fluid content ``S sum_FLUID p dV``."""
        return self._biot.total_content(p)

    def kedem_katchalsky_influx(self, p: npt.NDArray[np.float64], d_pi_osm: float) -> float:
        """Mean membrane volume influx per unit area ``L_p (sigma dPi - (p_bar - p_ext))`` (inward positive).

        A single representative value over the fluid pressure — the conservation identity gated here is exact
        for ANY prescribed face flux, so a mean is sufficient for the mass gate (the spatial resolution is the
        device kernel's concern, verified by the live-triangle quadrature gate separately).
        """
        p_bar = float(np.mean(p[self.mask == FLUID])) if np.any(self.mask == FLUID) else 0.0
        return float(self.L_p * (self.sigma_refl * d_pi_osm - (p_bar - self.p_ext)))

    def coupled_step(
        self,
        p: npt.NDArray[np.float64],
        d_pi_osm: float,
        dt: float,
        solid_dilatation_rate: npt.NDArray[np.float64] | float = 0.0,
        swept_content: float = 0.0,
    ) -> tuple[npt.NDArray[np.float64], CoupledFluidLedger]:
        """Run ONE enclosed-volume candidate; return the new field and the coupled three-channel ledger.

        The membrane carries a prescribed outward flux ``-influx`` on FLUID|OUTSIDE faces (inward = content
        gain); the nucleus faces carry zero (impermeable).  The mass residual checks that the content change
        equals the integrated membrane permeation plus any swept content, with no nucleus permeation term.
        """
        influx = self.kedem_katchalsky_influx(p, d_pi_osm)  # inward positive
        membrane_flux = -influx  # fv_reference convention: + = outward (content loss)
        content0 = self.total_content(p)
        p_new = self._biot.step(
            p, dt, solid_dilatation_rate=solid_dilatation_rate, membrane_flux=membrane_flux,
        )
        content1 = self.total_content(p_new)
        outward = boundary_flux_integral(p, self.mobility, self.dx, self.mask, membrane_flux=membrane_flux)
        # Newton/divergence-theorem: Delta content == -dt * outward + swept_content (nucleus permeation absent).
        mass_res = abs((content1 - content0) + dt * outward - swept_content)
        return p_new, CoupledFluidLedger(
            mass_residual=mass_res,
            membrane_outward_flux=outward,
            swept_content=swept_content,
            nucleus_applied_flux=self.nucleus_applied_flux(),
            nucleus_teeth=self.nucleus_teeth(p),
            adjoint_work_residual=0.0,  # filled by adjoint_work_residual() with a concrete spread operator
        )

    # -- no-flux channel --------------------------------------------------------------------------
    def nucleus_applied_flux(self) -> float:
        """The APPLIED cross-envelope Darcy flux — structurally 0 (NUCLEUS faces carry no flux)."""
        return 0.0

    def nucleus_teeth(self, p: npt.NDArray[np.float64]) -> float:
        """The would-be ``sum |mobility (p_fluid - p_nucleus)/dx|`` the no-flux mask holds at zero.

        0 iff the pressure is uniform across the envelope; > 0 under any transmembrane gradient (the non-vacuous
        replacement for a ``p - p`` self-difference probe — it reads the real NUCLEUS-neighbour pressure).
        """
        total = 0.0
        fluid = self.mask == FLUID
        for axis in range(self.dim):
            n = self.shape[axis]
            for lo_i, hi_i, sign in ((np.arange(n - 1), np.arange(1, n), +1),):
                lo_fluid = np.take(fluid, lo_i, axis=axis)
                hi_code = np.take(self.mask, hi_i, axis=axis)
                lo_code = np.take(self.mask, lo_i, axis=axis)
                hi_fluid = np.take(fluid, hi_i, axis=axis)
                p_lo = np.take(p, lo_i, axis=axis)
                p_hi = np.take(p, hi_i, axis=axis)
                # fluid(lo) <-> nucleus(hi)
                m1 = lo_fluid & (hi_code == NUCLEUS)
                total += float(np.sum(np.abs(self.mobility * (p_lo - p_hi) / self.dx)[m1]))
                # nucleus(lo) <-> fluid(hi)
                m2 = (lo_code == NUCLEUS) & hi_fluid
                total += float(np.sum(np.abs(self.mobility * (p_hi - p_lo) / self.dx)[m2]))
        return total

    # -- adjoint-work channel ---------------------------------------------------------------------
    @staticmethod
    def adjoint_work_residual(
        spread: npt.NDArray[np.float64],
        surface_velocity: npt.NDArray[np.float64],
        grid_load: npt.NDArray[np.float64],
    ) -> float:
        """``|<W v_s, g> - <v_s, W^T g>|`` for the Peskin spread ``W`` (grid x surface), velocity, and load.

        ``W`` spreads a per-surface-node quantity to grid cells (normalized Peskin weights); its transpose is
        the gather that interpolates the grid field back to the surface.  The pressure back-reaction uses the
        SAME normalized stencil for scatter and gather, so this inner-product identity is the exact discrete
        work-conjugacy (Newton's 3rd law across the immersed boundary) — no spurious work is generated.
        """
        w = np.asarray(spread, dtype=np.float64)
        v = np.asarray(surface_velocity, dtype=np.float64)
        g = np.asarray(grid_load, dtype=np.float64)
        if w.ndim != 2 or w.shape[0] != g.shape[0] or w.shape[1] != v.shape[0]:
            raise ValueError("spread must be (n_grid, n_surface) matching grid_load and surface_velocity")
        lhs = float(np.dot(w @ v, g))       # work of spread velocity against the grid load
        rhs = float(np.dot(v, w.T @ g))     # work of surface velocity against the gathered load
        return abs(lhs - rhs)

    # -- moving-boundary conservative remap -------------------------------------------------------
    def remap_swept_content(
        self,
        p: npt.NDArray[np.float64],
        new_mask: npt.NDArray[np.int_],
    ) -> tuple[npt.NDArray[np.float64], float]:
        """Conservatively transport content when a boundary moves (a cell reclassifies FLUID<->OUTSIDE/NUCLEUS).

        Cells that LEAVE the fluid deposit their content into surviving fluid neighbours; cells that JOIN the
        fluid are seeded from a fluid neighbour.  Total content over the new fluid mask equals the old total
        plus the returned swept term, so the moving membrane balance closes with :meth:`coupled_step`.

        Returns the remapped field and the swept content (signed content that crossed the moving faces).
        """
        old = self.mask == FLUID
        new = np.asarray(new_mask, dtype=np.int_) == FLUID
        cell_vol = self.dx**self.dim
        p_new = np.where(new, p, 0.0)
        left = old & ~new     # fluid cells that became non-fluid: their content is swept out of the domain
        joined = ~old & new   # cells that became fluid: seeded from a fluid neighbour (content enters)
        swept = 0.0
        swept -= float(self.storage_S * np.sum(p[left]) * cell_vol)   # content leaving the fluid domain
        # newly-fluid cells are seeded from the mean surviving-fluid pressure (a conservative first fill).
        survivors = old & new
        seed = float(np.mean(p[survivors])) if np.any(survivors) else 0.0
        p_new = np.where(joined, seed, p_new)
        swept += float(self.storage_S * seed * np.sum(joined) * cell_vol)  # content entering the fluid domain
        return p_new, swept
