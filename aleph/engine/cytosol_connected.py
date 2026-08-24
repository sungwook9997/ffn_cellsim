"""Cytosol FLUID_VOLUME -> CONNECTED: the Biot field + BOTH moving boundaries as one coupled candidate.

Roadmap ``AC_ENGINE_COMPLETION_ROADMAP_2026-07-23.md`` §3.1.  ``ac/engine/fluid_core.py`` already binds the
real Biot/Darcy solver (:class:`BiotSubstrateFluidSolver`) and the moving **impermeable nuclear envelope**
(:class:`NucleusPressureAdjointBoundary` via ``nucleus_cytosol_boundary``).  What was missing for CONNECTED is
the moving **semipermeable outer membrane** (``membrane_cytosol_boundary``, the Kedem-Katchalsky water flux
``J = L_p(sigma dPi - dP)``, backend :class:`~aleph.components.fluid.boundary.MembraneFluxBC`) and a single
enclosed-volume candidate that drives the field + BOTH boundaries and closes the three ledger channels
TOGETHER:

    MASS         : membrane permeation flux + nucleus geometric swept content (no nucleus permeation);
    NO-FLUX      : the impermeable nucleus applied flux is 0, its would-be discharge (teeth) is load-bearing;
    ADJOINT-WORK : the Peskin traction scatter / boundary-velocity gather is a transpose pair (Newton 3rd).

This module adds:

* :class:`MovingSemipermeableFluidBoundaryFacade` — the graph-bound ``membrane_cytosol_boundary`` FLUID_BOUNDARY,
  the semipermeable twin of ``fluid_core.MovingImpermeableFluidBoundaryFacade`` (nucleus).  It forwards the same
  ``update_geometry`` / ``accumulate`` (adjoint traction) / transaction / ledger hooks to an injected delegate.
* :class:`MembranePressureFluxAdjointBoundary` — the CUDA-lane concrete delegate binding ``MembraneFluxBC``
  (rebuild the hydraulic source on the live membrane) + ``PressureCoupling`` (scatter the pressure traction back
  onto the membrane surface — the adjoint of the boundary-velocity map) + the membrane integrated-flux mass
  channel.
* :class:`CytosolConnectedCandidate` — composes the cytosol state owner + membrane (semipermeable) + nucleus
  (impermeable) boundaries into ONE coupled candidate, exposes the :class:`CytosolFieldEndpoint` for the
  immersed-transfer connectors, and assembles the coupled mass / no-flux / adjoint-work ledger in one pass.

Osmotic driver: for the STATIC resting baseline (GATE A) the osmotic difference is a genuine van 't Hoff
constant ``Pi_0`` (spec §4a — at fixed resting volume the concentration is constant, so a fixed ``dPi`` is
physical, NOT a held-pressure artifact); a volume-changing run replaces it with the emergent ``dPi_eff`` from
the Fluid/fields osmotic vertical (post-GATE-A).  This module takes the driver as a value or callable and
invents no magnitude — ``L_p`` / ``Pi_0`` are sourced GAPs.

The facades/candidate are pure composition (CPU-importable; structural gates inject delegate doubles); the
concrete membrane delegate is CUDA-only and exercised on the native lane, exactly as ``fluid_core.py`` layers
its SEAMED facades over KERNEL_BOUND backends.

Sanity Gate:
    * ownership: the candidate drives cytosol + membrane + nucleus as three participants under one predicate;
      the membrane boundary owns caches but no component state (``commit_on_accept=False`` never exempts
      cache rollback).
    * boundary/sign: nucleus permeation is 0 (impermeable) — only the membrane channel carries permeation mass;
      the adjoint traction is the transpose of the boundary-velocity gather (no spurious work).
    * conservation: the coupled ledger accumulates membrane flux + swept content into the mass channel and the
      boundary work into the work channel in ONE assembly, so mass/no-flux/adjoint-work close simultaneously.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import warp as wp

from aleph.engine.contracts import CellArchitecture, ComponentRole, ConnectorFamily
from aleph.engine.fluid_core import (
    CYTOSOL_COMPONENT,
    NUCLEUS_COMPONENT,
    NUCLEUS_FLUID_BOUNDARY,
    FluidVolumeStateOwner,
    MovingBoundaryDelegate,
    MovingImpermeableFluidBoundaryFacade,
    SurfaceQuadratureState,
)
from aleph.engine.runtime import CytosolFieldEndpoint

MEMBRANE_COMPONENT = "membrane"
MEMBRANE_FLUID_BOUNDARY = "membrane_cytosol_boundary"

__all__ = [
    "MEMBRANE_COMPONENT",
    "MEMBRANE_FLUID_BOUNDARY",
    "MovingSemipermeableFluidBoundaryFacade",
    "MembranePressureFluxAdjointBoundary",
    "CytosolConnectedBindings",
    "CytosolConnectedCandidate",
]


def _device_is_cuda(array: object) -> bool:
    device = getattr(array, "device", None)
    return bool(getattr(device, "is_cuda", False))


# ======================================================================================================
# Semipermeable membrane fluid boundary (graph facade + concrete CUDA delegate)
# ======================================================================================================
@dataclass(frozen=True, slots=True)
class MovingSemipermeableFluidBoundaryFacade:
    """Graph-bound membrane--cytosol boundary: the semipermeable twin of the nucleus fluid boundary.

    ``update_geometry`` is the membrane-to-fluid kinematic + hydraulic-source map (rebuild the Kedem-Katchalsky
    flux on the live membrane).  ``accumulate`` is the reverse fluid-to-membrane pressure-traction (adjoint).
    Unlike the nucleus boundary this one is **semipermeable** — water crosses per ``J = L_p(sigma dPi - dP)``.
    """

    delegate: MovingBoundaryDelegate
    name: str = MEMBRANE_FLUID_BOUNDARY
    component_a: str = MEMBRANE_COMPONENT
    component_b: str = CYTOSOL_COMPONENT
    impermeable: bool = False
    adjoint_transfer_required: bool = True

    def __post_init__(self) -> None:
        if self.name != MEMBRANE_FLUID_BOUNDARY:
            raise ValueError(f"expected field coupler {MEMBRANE_FLUID_BOUNDARY!r}")
        if {self.component_a, self.component_b} != {MEMBRANE_COMPONENT, CYTOSOL_COMPONENT}:
            raise ValueError("membrane fluid-boundary endpoints must be membrane and cytosol")
        if self.impermeable:
            raise ValueError("the outer membrane fluid boundary is semipermeable, not impermeable")
        if not self.adjoint_transfer_required:
            raise ValueError("membrane fluid boundary requires an adjoint transfer")

    def update_geometry(self, surface_position: wp.array, surface_velocity: wp.array) -> None:
        """Forward moving geometry + rebuild the hydraulic membrane source to the field implementation."""
        self.delegate.update_geometry(surface_position, surface_velocity)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Add fluid pressure traction to the membrane surface force (adjoint back-reaction)."""
        self.delegate.accumulate_boundary(pos, force)

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        """Alias matching the surface-body fluid-boundary vocabulary."""
        self.accumulate(pos, force)

    def snapshot_candidate(self) -> None:
        self.delegate.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        self.delegate.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.delegate.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward membrane permeation mass + adjoint work terms to the global ledger."""
        self.delegate.accumulate_ledger(ledger)


class MembranePressureFluxAdjointBoundary:
    """Bind the semipermeable-membrane :class:`MovingBoundaryDelegate` to ``MembraneFluxBC`` + ``PressureCoupling``.

    * ``update_geometry`` rebuilds the Kedem-Katchalsky hydraulic source on the LIVE membrane triangles for the
      current osmotic difference (``membrane_bc.apply(d_pi_osm)``).  The osmotic driver is a value or callable;
      for GATE A it is the van 't Hoff constant ``Pi_0`` (spec §4a) — physical at fixed resting volume.
    * ``accumulate_boundary`` scatters the fluid pressure traction ``-alpha V_node grad(p)`` onto the
      membrane-owned surface force through the existing Peskin primitive — the work-conjugate adjoint of the
      boundary-velocity map (Newton 3rd: the field feels the opposite via its own ``PressureCoupling``).
    * the membrane permeation is the mass channel of the conservation ledger
      (``membrane_bc.integrated_flux_d``); transaction/ledger are forwarded to injected participants.
    """

    def __init__(
        self,
        membrane_bc: object,
        pressure_coupling: object,
        node_volume_d: wp.array,
        transaction: object,
        ledger: object,
        *,
        osmotic_difference: Callable[[], float] | float,
    ) -> None:
        for label, dep, attr in (
            ("membrane_bc", membrane_bc, "apply"),
            ("pressure_coupling", pressure_coupling, "accumulate"),
            ("transaction", transaction, "snapshot_candidate"),
            ("ledger", ledger, "accumulate_ledger"),
        ):
            if not hasattr(dep, attr):
                raise TypeError(f"{label} must provide .{attr}(...) for the membrane fluid boundary")
        if getattr(node_volume_d, "dtype", None) != wp.float64:
            raise TypeError("node_volume_d must have dtype wp.float64")
        if len(getattr(node_volume_d, "shape", ())) != 1:
            raise ValueError("node_volume_d must be a rank-1 per-surface-node array")
        if not _device_is_cuda(node_volume_d):
            raise ValueError("node_volume_d must be a Warp CUDA device array")
        if not (callable(osmotic_difference) or isinstance(osmotic_difference, (int, float))):
            raise TypeError("osmotic_difference must be a constant Pi_0 or a callable returning dPi")
        self.membrane_bc = membrane_bc
        self.pressure_coupling = pressure_coupling
        self.node_volume_d = node_volume_d
        self.transaction = transaction
        self.ledger = ledger
        self._osmotic_difference = osmotic_difference

    def _d_pi_osm(self) -> float:
        d_pi = self._osmotic_difference() if callable(self._osmotic_difference) else self._osmotic_difference
        if not math.isfinite(float(d_pi)):
            raise ValueError("osmotic difference dPi must be finite")
        return float(d_pi)

    def update_geometry(self, surface_position: wp.array, surface_velocity: wp.array) -> None:
        """Rebuild the membrane hydraulic source for the current osmotic difference on the live membrane."""
        self.membrane_bc.apply(self._d_pi_osm())

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Scatter the fluid pressure traction onto the membrane surface force (adjoint back-reaction)."""
        self.pressure_coupling.accumulate(SurfaceQuadratureState(pos, self.node_volume_d), force)

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        self.accumulate(pos, force)

    def membrane_flux(self) -> wp.array | None:
        """The device scalar of net membrane permeation this step (the mass channel; None before first apply)."""
        return getattr(self.membrane_bc, "integrated_flux_d", None)

    def snapshot_candidate(self) -> None:
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger.accumulate_ledger(ledger)


# ======================================================================================================
# The coupled enclosed-volume candidate
# ======================================================================================================
def _find_connector(architecture: CellArchitecture, name: str):
    for connector in architecture.connectors:
        if connector.name == name:
            return connector
    raise ValueError(f"architecture must register {name!r}")


def _validate_fluid_boundary(architecture: CellArchitecture, name: str, endpoints: set[str]) -> None:
    connector = _find_connector(architecture, name)
    if {connector.component_a, connector.component_b} != endpoints:
        raise ValueError(f"{name!r} must join {sorted(endpoints)}")
    if connector.family is not ConnectorFamily.FLUID_BOUNDARY:
        raise ValueError(f"{name!r} must use ConnectorFamily.FLUID_BOUNDARY")
    if connector.kinetics or connector.commit_on_accept:
        raise ValueError(f"{name!r} is a non-kinetic field coupler")
    if not connector.bidirectional or not connector.adjoint_transfer_required:
        raise ValueError(f"{name!r} must be bidirectional with adjoint transfer")


@dataclass(frozen=True, slots=True)
class CytosolConnectedBindings:
    """Non-owning graph/field binding supplied for one coupled enclosed-volume candidate."""

    membrane_fluid_boundary: MovingSemipermeableFluidBoundaryFacade
    nucleus_fluid_boundary: MovingImpermeableFluidBoundaryFacade
    membrane_surface_position_d: wp.array
    membrane_surface_velocity_d: wp.array
    membrane_surface_force_d: wp.array
    nucleus_surface_position_d: wp.array
    nucleus_surface_velocity_d: wp.array
    nucleus_surface_force_d: wp.array

    def __post_init__(self) -> None:
        if self.membrane_fluid_boundary.name != MEMBRANE_FLUID_BOUNDARY:
            raise ValueError(f"expected graph connector {MEMBRANE_FLUID_BOUNDARY!r}")
        if self.nucleus_fluid_boundary.name != NUCLEUS_FLUID_BOUNDARY:
            raise ValueError(f"expected graph connector {NUCLEUS_FLUID_BOUNDARY!r}")


@dataclass(frozen=True, slots=True)
class CytosolConnectedCandidate:
    """One coupled enclosed-volume candidate over the cytosol field + membrane + nucleus moving boundaries.

    The caller zeroes the boundary force arrays before this candidate and repeats it to the interface residual;
    the final geometry refresh leaves both boundaries synchronized with the new field candidate.
    """

    architecture: CellArchitecture
    cytosol: FluidVolumeStateOwner

    def __post_init__(self) -> None:
        cytosol = self.architecture.component(CYTOSOL_COMPONENT)
        if cytosol.role is not ComponentRole.FLUID_VOLUME:
            raise ValueError("registered cytosol must have ComponentRole.FLUID_VOLUME")
        self.architecture.component(MEMBRANE_COMPONENT)
        self.architecture.component(NUCLEUS_COMPONENT)
        _validate_fluid_boundary(
            self.architecture, MEMBRANE_FLUID_BOUNDARY, {MEMBRANE_COMPONENT, CYTOSOL_COMPONENT}
        )
        _validate_fluid_boundary(
            self.architecture, NUCLEUS_FLUID_BOUNDARY, {NUCLEUS_COMPONENT, CYTOSOL_COMPONENT}
        )
        if self.cytosol.name != CYTOSOL_COMPONENT:
            raise ValueError("CytosolConnectedCandidate.cytosol must be the registered cytosol state owner")

    def make_immersed_transfer_endpoint(self, residual_d: wp.array) -> CytosolFieldEndpoint:
        """Expose the live pressure + a solver-owned residual for the IMMERSED_TRANSFER connectors."""
        return self.cytosol.make_immersed_transfer_endpoint(residual_d)

    def coupled_candidate_iteration(self, bindings: CytosolConnectedBindings, dt_phys: float) -> None:
        """Run one coupled enclosed-volume iteration: both boundary geometries -> field solve -> both tractions.

        Order (partitioned Gauss-Seidel, matching ``FluidCore.candidate_iteration``): update the membrane and
        nucleus moving geometry (rebuild the hydraulic source + reclassify/remap), solve the Biot field for the
        new boundaries, then scatter each boundary's adjoint pressure traction back onto its surface force.
        """
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        membrane = bindings.membrane_fluid_boundary
        nucleus = bindings.nucleus_fluid_boundary
        membrane.update_geometry(
            bindings.membrane_surface_position_d, bindings.membrane_surface_velocity_d
        )
        nucleus.update_geometry(
            bindings.nucleus_surface_position_d, bindings.nucleus_surface_velocity_d
        )
        self.cytosol.solve_candidate(dt_phys)
        membrane.accumulate(bindings.membrane_surface_position_d, bindings.membrane_surface_force_d)
        nucleus.accumulate(bindings.nucleus_surface_position_d, bindings.nucleus_surface_force_d)

    def transaction_participants(self, bindings: CytosolConnectedBindings) -> tuple[object, object, object]:
        """Return the cytosol field owner + both moving boundaries participating in the global transaction."""
        return self.cytosol, bindings.membrane_fluid_boundary, bindings.nucleus_fluid_boundary

    def snapshot_candidate(self, bindings: CytosolConnectedBindings) -> None:
        for participant in self.transaction_participants(bindings):
            participant.snapshot_candidate()

    def rollback(self, bindings: CytosolConnectedBindings, accepted: wp.array) -> None:
        for participant in self.transaction_participants(bindings):
            participant.rollback(accepted)

    def commit_irreversible(
        self, bindings: CytosolConnectedBindings, accepted: wp.array, dt_phys: float, rng_seed: int
    ) -> None:
        for participant in self.transaction_participants(bindings):
            participant.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_coupled_ledger(self, bindings: CytosolConnectedBindings, ledger: object) -> None:
        """Collect the field, membrane (permeation mass + adjoint work), and nucleus (swept mass) terms.

        One assembly closes the three channels together: the membrane contributes permeation to the mass
        channel + boundary work; the nucleus contributes only geometric swept content (its permeation is 0 —
        impermeable); the field owner contributes its residual/mass.  Acceptance is decided downstream on the
        assembled ledger, never here.
        """
        self.cytosol.accumulate_ledger(ledger)
        bindings.membrane_fluid_boundary.accumulate_ledger(ledger)
        bindings.nucleus_fluid_boundary.accumulate_ledger(ledger)
