"""The ``immersed_transfer`` connector family: any solid component immersed in the cytosol field.

WHAT THIS CLOSES.  Seven edges of :func:`~aleph.engine.contracts.reference_cell_architecture` carry
``ConnectorFamily.IMMERSED_TRANSFER``, every one of them ``<solid> -> cytosol``::

    surface_porous_transfer          cortex                 -> cytosol
    sf_cytosol_transfer              sf_arc                 -> cytosol
    mt_cytosol_transfer              microtubule            -> cytosol
    if_cytosol_transfer              intermediate_filament  -> cytosol
    lamellipodium_cytosol_transfer   lamellipodium          -> cytosol
    filopodium_cytosol_transfer      filopodium             -> cytosol
    nmii_cytosol_transfer            nmii                   -> cytosol

One of the seven had a runtime and six did not, which is a sixth of the whole 36-connector census.  The
one that existed — ``CortexCytosolPorousTransfer``, formerly in :mod:`aleph.engine.interior_column_slice`
— was already generic in substance: it scatters ``-alpha * V_node * grad(p)`` with the Peskin-4 stencil
onto whatever solid view it is handed, and nothing in that arithmetic is cortical.  What pinned it to the
cortex was a pair of hard-coded name checks in its constructor.  So this module is that class with the
guard rewritten rather than six new implementations of one physics.

**The guard is now derived, not repeated.**  The constructor takes the edge NAME, looks the
:class:`~aleph.engine.contracts.ConnectorContract` up in the architecture, and reads its endpoints from
there.  A caller cannot mis-state which components it joins, because it is not asked; and a name whose
declared family is not ``IMMERSED_TRANSFER`` is refused.  The previous guard could only ever protect the
one pair it named, and would have had to be copied — and eventually mis-copied — six times.

WHICH PHYSICS THIS IS, stated plainly because the decision was frozen and is now open.
``darcy-or-stokes-for-the-cytosol`` is an inherited open question, and this code already answers it in
one direction: the coupling implemented here is **Biot poroelastic (Darcy-type)** — a pressure gradient
scattered onto the solid as ``-alpha * V * grad(p)``, with the work-conjugate reaction reaching the fluid
through the ``div(v_s)`` source the moving solid imprints on the Biot update.  That is not a Stokes drag
on an immersed body.  Extending the family does not choose the answer; it carries the answer that was
already running for ``cortex`` to the other six solids, which is the only self-consistent thing to do
until the question is decided.  **Ratifying or changing it is a PI decision, and it is one decision for
all seven edges, not seven.**

WHY IT OWNS NO REVERSIBLE STATE.  The authoritative pressure field belongs to the ``cytosol`` component,
so the transaction hooks here are structural no-ops: the field's own transaction snapshots and
reject-restores it.  The hooks exist because
:class:`~aleph.engine.world.CellWorldTransaction` treats the transaction API as all-or-nothing — a
connector implementing two of the three is a ``TypeError``, not a partial participant.

Sanity Gate (before first execution; :mod:`aleph.tests.ac.engine.test_immersed_transfer`):
  * dimensions: ``node_volume_d`` is a per-node control volume in um^3; the injected coupling produces
    pN on the solid's force array.  This module introduces no constant and no unit conversion of its own.
  * boundary cases: a solid with zero nodes scatters nothing (the injected coupling's own launch is
    dimension-zero); a name outside the family, or absent from the architecture, raises at construction
    rather than binding wrong.
  * conservation invariant: the scatter is ``wp.atomic_add`` and never overwrites — the caller owns force
    zeroing, exactly as for every other accumulate in the composed candidate solve.  Gather and scatter
    are the transpose pair, which is why the contract carries ``bidirectional`` and
    ``adjoint_transfer_required``; a one-sided scatter is what the balance gate is there to catch.
  * sign sense: the traction is ``-alpha * V * grad(p)`` — a solid in a pressure gradient is pushed from
    high to low pressure, so the sign is carried by the injected coupling and is not re-derived here.
  * precision / measurement protocol: no host readback and no magnitude is produced by this module; the
    only thing it can be asked for is which calls it made, which the recording lists carry for the gate.

CPU-importable: it allocates nothing and launches nothing itself, delegating every launch to the injected
pressure coupling, so the structural gate drives it with recording doubles.
"""

from __future__ import annotations

import warp as wp

from aleph.engine.contracts import CellArchitecture, ConnectorFamily, reference_cell_architecture
from aleph.engine.fluid_core import SurfaceQuadratureState

__all__ = [
    "IMMERSED_TRANSFER_EDGES",
    "ImmersedPorousTransfer",
    "build_immersed_porous_transfer",
    "immersed_transfer_edges",
]


def immersed_transfer_edges(architecture: CellArchitecture | None = None) -> tuple[str, ...]:
    """Return every declared ``immersed_transfer`` edge name, in architecture order.

    Read from the contract rather than listed here, so a new solid component that declares a cytosol
    transfer is covered without editing this module — and so this docstring's list cannot go stale.
    """
    arch = architecture or reference_cell_architecture()
    return tuple(
        c.name for c in arch.connectors if c.family is ConnectorFamily.IMMERSED_TRANSFER
    )


#: Convenience snapshot of the family at import time (the reference architecture is a constant).
IMMERSED_TRANSFER_EDGES = immersed_transfer_edges()


class ImmersedPorousTransfer:
    """One ``immersed_transfer`` edge: the cytosol pressure traction on an immersed solid component.

    Binds a real ``PressureCoupling`` primitive (it launches the masked pressure-gradient kernel and the
    Peskin-4 interpolated scatter, putting ``-alpha * V_node * grad(p)`` on the solid's OWN force array).
    The fluid feels the work-conjugate reaction through the Biot ``div(v_s)`` source the moving solid
    imprints, so gather and scatter are the transpose pair (Newton's 3rd law).

    A NON-kinetic field coupler: the authoritative pressure field belongs to the ``cytosol`` component, so
    this connector owns no reversible per-candidate state and its transaction hooks are structural no-ops.
    """

    #: fidelity marker (a real per-node Peskin pressure scatter, never a lumped whole-body pressure).
    per_node_pressure_traction = True
    aggregate_or_lumped = False

    def __init__(
        self,
        *,
        name: str,
        pressure_coupling: object,
        node_volume_d: wp.array,
        architecture: CellArchitecture | None = None,
    ) -> None:
        arch = architecture or reference_cell_architecture()
        contract = next((c for c in arch.connectors if c.name == name), None)
        if contract is None:
            raise ValueError(f"connector {name!r} is not declared in the architecture")
        if contract.family is not ConnectorFamily.IMMERSED_TRANSFER:
            raise ValueError(
                f"connector {name!r} is family {contract.family.value!r}, not "
                f"{ConnectorFamily.IMMERSED_TRANSFER.value!r}; an immersed transfer cannot stand in for it"
            )
        if not callable(getattr(pressure_coupling, "accumulate", None)):
            raise TypeError("immersed transfer requires a PressureCoupling exposing accumulate(state, force)")

        self.name = contract.name
        # Endpoints come from the contract, never from the caller: a runtime cannot disagree with the
        # declaration it is bound under.
        self.component_a = contract.component_a
        self.component_b = contract.component_b
        self.bidirectional = contract.bidirectional
        self.adjoint_transfer_required = contract.adjoint_transfer_required
        self._pressure_coupling = pressure_coupling
        self._node_volume_d = node_volume_d
        self.transfer_calls: list[object] = []
        self.ledger_calls: list[object] = []

    def accumulate_transfer(self, solid: object, cytosol: object) -> None:
        """Scatter the cytosol pressure traction ``-alpha*V*grad(p)`` onto the solid-owned force array.

        ``solid`` is a non-owning view of the immersed component (``position_d`` / ``force_d``); ``cytosol``
        is the public field endpoint, present for signature parity — the live field the traction reads is
        the one the injected pressure coupling already holds.  The scatter never overwrites; the caller
        owns force zeroing.
        """
        self.transfer_calls.append((solid, cytosol))
        self._pressure_coupling.accumulate(
            SurfaceQuadratureState(solid.position_d, self._node_volume_d), solid.force_d
        )

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Mechanics-hook alias: scatter the pressure traction onto ``force`` at solid node positions ``pos``."""
        self._pressure_coupling.accumulate(SurfaceQuadratureState(pos, self._node_volume_d), force)

    def snapshot_candidate(self) -> None:
        """No reversible connector state: the authoritative pressure field belongs to the cytosol component."""

    def rollback(self, accepted: wp.array) -> None:
        """No-op: the cytosol field's own transaction reject-restores the pressure (no connector cache here)."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """No-op: a non-kinetic field coupler advances no irreversible state."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Record the ledger handle (the adjoint-work term reduces device-side on the native lane)."""
        self.ledger_calls.append(ledger)


def build_immersed_porous_transfer(
    *,
    name: str,
    pressure_coupling: object,
    node_volume_d: wp.array,
    architecture: CellArchitecture | None = None,
) -> ImmersedPorousTransfer:
    """Build one immersed solid<->cytosol porous transfer for the declared edge ``name``.

    Args:
        name: the declared connector edge, e.g. ``"sf_cytosol_transfer"``.  Must carry family
            ``IMMERSED_TRANSFER`` in the architecture; its endpoints are read from there.
        pressure_coupling: the real ``PressureCoupling`` primitive over the live Biot field.
        node_volume_d: per-node immersed control volume [um^3] for the solid's nodes, in node order.
        architecture: optional composition (defaults to :func:`reference_cell_architecture`).

    Returns:
        An :class:`ImmersedPorousTransfer` ready to bind at the ``name`` connector slot.
    """
    return ImmersedPorousTransfer(
        name=name, pressure_coupling=pressure_coupling, node_volume_d=node_volume_d,
        architecture=architecture,
    )
