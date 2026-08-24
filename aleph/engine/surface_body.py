"""Representation-neutral runtime seam for the membrane/cortex Surface Body.

This module fixes ownership and coupling before any new surface constitutive law is introduced:

* ``membrane`` and ``cortex`` remain distinct state-owning components with independent CUDA arrays;
* :class:`SurfaceBody` is a coupled-mechanics facade, never a third state owner;
* the ``membrane_erm_cortex`` connector is supplied by the graph for each step, so the facade cannot own or
  privately commit ERM state;
* a registered membrane--cytosol boundary is required before construction, rather than hiding membrane
  pressure/flux inside the cortex porous-transfer edge; and
* mechanics and transaction hooks receive existing device arrays/contributors.  No host-side physics path is
created here.

The Cortex owner may wrap the native 70,686-filament implementation or a separately validated
condensation candidate.  This seam does not authorize the latter to replace the authoritative
full-population baseline.

The adapter at the bottom can reuse the landed :class:`ac.cell.compartments.MembraneCompartment` Helfrich and
area kernels when that legacy object is rebuilt against a *local* membrane array with ``node_off == 0`` and
``n_erm == 0``.  ERM is deliberately excluded because its authoritative state belongs to the connector graph.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Protocol, runtime_checkable

import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
)
from aleph.engine.runtime import (
    CytosolFieldEndpoint,
    ImmersedTransferConnector,
    LedgerContributor,
    MechanicsContributor,
    TransactionParticipant,
)

MEMBRANE_COMPONENT = "membrane"
CORTEX_COMPONENT = "cortex"
CYTOSOL_COMPONENT = "cytosol"
ERM_CONNECTOR = "membrane_erm_cortex"
MEMBRANE_FLUID_BOUNDARY = "membrane_cytosol_boundary"
SURFACE_POROUS_TRANSFER = "surface_porous_transfer"


def _device_is_cuda(array: object) -> bool:
    """Return whether an array-like object declares a CUDA device."""
    device = getattr(array, "device", None)
    return bool(getattr(device, "is_cuda", False))


def storage_key(array: object) -> tuple[object, ...]:
    """Return a conservative storage-identity key for Warp arrays and structural test doubles.

    Public because array-ownership exclusivity is a whole-composition property, not a membrane/cortex one:
    :func:`aleph.engine.cortex_state.assert_component_state_disjoint` applies the same identity across
    every registered component (``COMPARTMENT_VALIDATION_TRACKS_2026-07-25.md`` T2 step iii).
    """
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    device = getattr(array, "device", None)
    return ("device-ptr", str(device), int(ptr))


#: Backwards-compatible private alias — this module and its tests used the underscore name first.
_storage_key = storage_key


def _validate_vec3_device_array(array: object, *, label: str) -> None:
    """Validate the minimum CUDA-array contract without reading device state."""
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != wp.vec3d:
        raise TypeError(f"{label} must have dtype wp.vec3d")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != 1:
        raise ValueError(f"{label} must be a one-dimensional vec3 array")
    if int(shape[0]) <= 0:
        raise ValueError(f"{label} must contain at least one surface vertex")


def _validate_topology_pair_array(array: object, *, label: str, cols: int) -> int:
    """Validate a CUDA ``(rows, cols)`` int32 connectivity array; return the row count.

    Used for the cortex link/segment table ``(L, 2)`` and the Cytosim bending triples ``(T, 3)``.
    No device data is read; only the array metadata contract is enforced on the host.
    """
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != wp.int32:
        raise TypeError(f"{label} must have dtype wp.int32")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != 2 or int(shape[1]) != int(cols):
        raise ValueError(f"{label} must be a two-dimensional (rows, {cols}) int32 array")
    return int(shape[0])


def _validate_scalar_f64_array(array: object, *, label: str, rows: int) -> None:
    """Validate a CUDA one-dimensional float64 parameter array of a fixed length."""
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != wp.float64:
        raise TypeError(f"{label} must have dtype wp.float64")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != 1:
        raise ValueError(f"{label} must be a one-dimensional float64 array")
    if int(shape[0]) != int(rows):
        raise ValueError(f"{label} length {int(shape[0])} must match its topology row count {int(rows)}")


@runtime_checkable
class SurfacePairConnector(Protocol):
    """Graph-owned membrane/cortex connector view used by the coupled surface solve."""

    name: str
    component_a: str
    component_b: str

    def accumulate_pair(
        self,
        pos_a: wp.array,
        force_a: wp.array,
        pos_b: wp.array,
        force_b: wp.array,
    ) -> None:
        """Add equal-and-opposite connector mechanics to two independently owned states."""

    def snapshot_candidate(self) -> None:
        """Snapshot graph-owned connector state before candidate mutation."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore connector state under the global rejected predicate."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit connector kinetics only under the final accepted predicate."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute connector force/work/topology terms without deciding acceptance."""


@runtime_checkable
class MembraneFluidBoundary(Protocol):
    """Field-owned membrane boundary view for pressure, flux, and moving-boundary work."""

    name: str
    component_a: str
    component_b: str

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        """Add current fluid traction to the membrane-owned force array."""

    def snapshot_candidate(self) -> None:
        """Snapshot boundary state before a coupled candidate step."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore boundary state under the global rejected predicate."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit boundary transport only after global acceptance."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute pressure/flux/work terms without deciding acceptance."""


@dataclass(frozen=True, slots=True)
class SurfacePorousEndpointView:
    """Non-owning cortical surface view consumed by the porous-fluid transfer."""

    component: str
    position_d: wp.array
    force_d: wp.array

    def __post_init__(self) -> None:
        if self.component != CORTEX_COMPONENT:
            raise ValueError("surface porous endpoint must be the registered cortex")
        _validate_vec3_device_array(self.position_d, label="cortex.position_d")
        _validate_vec3_device_array(self.force_d, label="cortex.force_d")
        if self.position_d.shape != self.force_d.shape:
            raise ValueError("cortex porous endpoint position and force shapes must match")
        if str(self.position_d.device) != str(self.force_d.device):
            raise ValueError("cortex porous endpoint arrays must share one CUDA device")
        if _storage_key(self.position_d) == _storage_key(self.force_d):
            raise ValueError("cortex porous endpoint position and force must not alias")


@dataclass(frozen=True, slots=True)
class SurfaceStepBindings:
    """Non-owning views supplied by the connector/field registries for one candidate step."""

    erm_connector: SurfacePairConnector
    membrane_fluid_boundary: MembraneFluidBoundary
    cortex_cytosol_transfer: ImmersedTransferConnector
    cytosol: CytosolFieldEndpoint

    def __post_init__(self) -> None:
        if self.erm_connector.name != ERM_CONNECTOR:
            raise ValueError(f"expected graph connector {ERM_CONNECTOR!r}")
        if {self.erm_connector.component_a, self.erm_connector.component_b} != {
            MEMBRANE_COMPONENT,
            CORTEX_COMPONENT,
        }:
            raise ValueError("ERM runtime endpoints must be membrane and cortex")
        if self.membrane_fluid_boundary.name != MEMBRANE_FLUID_BOUNDARY:
            raise ValueError(f"expected field coupler {MEMBRANE_FLUID_BOUNDARY!r}")
        if {
            self.membrane_fluid_boundary.component_a,
            self.membrane_fluid_boundary.component_b,
        } != {MEMBRANE_COMPONENT, CYTOSOL_COMPONENT}:
            raise ValueError("membrane fluid-boundary endpoints must be membrane and cytosol")
        if self.cortex_cytosol_transfer.name != SURFACE_POROUS_TRANSFER:
            raise ValueError(f"expected field coupler {SURFACE_POROUS_TRANSFER!r}")
        if {
            self.cortex_cytosol_transfer.component_a,
            self.cortex_cytosol_transfer.component_b,
        } != {CORTEX_COMPONENT, CYTOSOL_COMPONENT}:
            raise ValueError("surface porous-transfer endpoints must be cortex and cytosol")
        if self.cytosol.component != CYTOSOL_COMPONENT:
            raise ValueError("surface porous transfer requires the registered cytosol field")


@dataclass(frozen=True, slots=True)
class SurfaceComponentStateOwner:
    """One registered surface component and the CUDA arrays it exclusively owns.

    The mechanics, transaction, and ledger delegates may reuse existing Warp implementations.  This wrapper
    does not compute forces on the host and does not concatenate another component's arrays into its state.
    """

    name: str
    position_d: wp.array
    force_d: wp.array
    mechanics: MechanicsContributor
    transaction: TransactionParticipant
    ledger: LedgerContributor

    def __post_init__(self) -> None:
        if self.name not in {MEMBRANE_COMPONENT, CORTEX_COMPONENT}:
            raise ValueError("surface state owner name must be 'membrane' or 'cortex'")
        _validate_vec3_device_array(self.position_d, label=f"{self.name}.position_d")
        _validate_vec3_device_array(self.force_d, label=f"{self.name}.force_d")
        if self.position_d.shape != self.force_d.shape:
            raise ValueError(f"{self.name} position and force arrays must have identical shape")
        if str(self.position_d.device) != str(self.force_d.device):
            raise ValueError(f"{self.name} position and force arrays must share a CUDA device")
        if _storage_key(self.position_d) == _storage_key(self.force_d):
            raise ValueError(f"{self.name} position and force arrays must not alias")
        n_vertices = getattr(self.mechanics, "n_vertices", None)
        if n_vertices is not None and int(n_vertices) != int(self.position_d.shape[0]):
            raise ValueError(
                f"{self.name} mechanics expects {int(n_vertices)} vertices but state owns "
                f"{int(self.position_d.shape[0])}"
            )

    def accumulate(self) -> None:
        """Launch this component's existing mechanics against its own CUDA arrays."""
        self.mechanics.accumulate(self.position_d, self.force_d)

    def snapshot_candidate(self) -> None:
        """Forward the global transaction snapshot hook to the state owner."""
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Forward the scheduler-owned predicate without inspecting it on the host."""
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Forward accepted-step internal-state kinetics without advancing a private clock."""
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward device ledger contributions without owning acceptance."""
        self.ledger.accumulate_ledger(ledger)


def _find_connector(architecture: CellArchitecture, name: str) -> ConnectorContract:
    """Return a named build-time connector or raise a targeted construction error."""
    for connector in architecture.connectors:
        if connector.name == name:
            return connector
    raise ValueError(f"architecture must register {name!r}")


def _validate_architecture(architecture: CellArchitecture) -> None:
    """Validate the graph edges required by the first Surface Body runtime seam."""
    membrane = architecture.component(MEMBRANE_COMPONENT)
    cortex = architecture.component(CORTEX_COMPONENT)
    architecture.component(CYTOSOL_COMPONENT)
    if membrane.role is not ComponentRole.SURFACE_BODY:
        raise ValueError("registered membrane must have ComponentRole.SURFACE_BODY")
    if cortex.role is not ComponentRole.SURFACE_BODY:
        raise ValueError("registered cortex must have ComponentRole.SURFACE_BODY")
    if not membrane.owns_geometry or not cortex.owns_geometry:
        raise ValueError("membrane and cortex must independently own geometry")

    erm = _find_connector(architecture, ERM_CONNECTOR)
    if erm.family is not ConnectorFamily.ERM:
        raise ValueError("membrane_erm_cortex must use ConnectorFamily.ERM")
    if {erm.component_a, erm.component_b} != {MEMBRANE_COMPONENT, CORTEX_COMPONENT}:
        raise ValueError("membrane_erm_cortex must join membrane and cortex")
    if not erm.kinetics or not erm.commit_on_accept:
        raise ValueError("membrane_erm_cortex must be kinetic and accepted-step committed")
    if not erm.adjoint_transfer_required or not erm.bidirectional:
        raise ValueError("membrane_erm_cortex must be bidirectional with adjoint transfer")

    fluid_boundary = _find_connector(architecture, MEMBRANE_FLUID_BOUNDARY)
    if {fluid_boundary.component_a, fluid_boundary.component_b} != {
        MEMBRANE_COMPONENT,
        CYTOSOL_COMPONENT,
    }:
        raise ValueError("membrane_cytosol_boundary must join membrane and cytosol")
    # ``FLUID_BOUNDARY`` is the planned family; IMMERSED_TRANSFER is the non-invasive Wave-1 bridge until the
    # common enum grows that value.  Compare the string so this module accepts the future StrEnum member.
    if fluid_boundary.family.value not in {"fluid_boundary", ConnectorFamily.IMMERSED_TRANSFER.value}:
        raise ValueError("membrane_cytosol_boundary must be a fluid-boundary/immersed-transfer coupler")
    if fluid_boundary.kinetics or fluid_boundary.commit_on_accept:
        raise ValueError("membrane_cytosol_boundary is a non-kinetic field coupler")
    if not fluid_boundary.adjoint_transfer_required or not fluid_boundary.bidirectional:
        raise ValueError("membrane_cytosol_boundary must be bidirectional with adjoint transfer")

    porous = _find_connector(architecture, SURFACE_POROUS_TRANSFER)
    if porous.family is not ConnectorFamily.IMMERSED_TRANSFER:
        raise ValueError("surface_porous_transfer must use ConnectorFamily.IMMERSED_TRANSFER")
    if {porous.component_a, porous.component_b} != {CORTEX_COMPONENT, CYTOSOL_COMPONENT}:
        raise ValueError("surface_porous_transfer must join cortex and cytosol")
    if porous.kinetics or porous.commit_on_accept:
        raise ValueError("surface_porous_transfer must be a non-kinetic field coupler")
    if not porous.adjoint_transfer_required or not porous.bidirectional:
        raise ValueError("surface_porous_transfer must be bidirectional with adjoint transfer")


@dataclass(frozen=True, slots=True)
class SurfaceBody:
    """Coupled facade over two state owners; ERM remains an external graph binding.

    The facade owns no position, force, ERM, kinetic epoch, or accepted-time array.  It validates the build and
    supplies deterministic mechanics/transaction ordering over objects owned by the component, connector, and
    field registries.
    """

    architecture: CellArchitecture
    membrane: SurfaceComponentStateOwner
    cortex: SurfaceComponentStateOwner

    def __post_init__(self) -> None:
        _validate_architecture(self.architecture)
        if self.membrane.name != MEMBRANE_COMPONENT:
            raise ValueError("SurfaceBody.membrane must be the registered membrane state owner")
        if self.cortex.name != CORTEX_COMPONENT:
            raise ValueError("SurfaceBody.cortex must be the registered cortex state owner")
        if self.membrane is self.cortex:
            raise ValueError("membrane and cortex must be distinct state owners")
        if _storage_key(self.membrane.position_d) == _storage_key(self.cortex.position_d):
            raise ValueError("membrane and cortex must have distinct position kinematics")
        if _storage_key(self.membrane.force_d) == _storage_key(self.cortex.force_d):
            raise ValueError("membrane and cortex must have distinct force ownership")

    def accumulate_mechanics(self, bindings: SurfaceStepBindings) -> None:
        """Assemble component, graph-connector, and field-boundary loads into separate force arrays.

        The caller owns force-array zeroing.  Every delegate is expected to launch Warp kernels; this method
        only wires arrays and never reads or computes authoritative mechanics on the host.
        """
        self.membrane.accumulate()
        self.cortex.accumulate()
        bindings.erm_connector.accumulate_pair(
            self.membrane.position_d,
            self.membrane.force_d,
            self.cortex.position_d,
            self.cortex.force_d,
        )
        bindings.membrane_fluid_boundary.accumulate_boundary(
            self.membrane.position_d,
            self.membrane.force_d,
        )
        bindings.cortex_cytosol_transfer.accumulate_transfer(
            SurfacePorousEndpointView(
                CORTEX_COMPONENT,
                self.cortex.position_d,
                self.cortex.force_d,
            ),
            bindings.cytosol,
        )

    def transaction_participants(
        self, bindings: SurfaceStepBindings
    ) -> tuple[
        TransactionParticipant,
        TransactionParticipant,
        SurfacePairConnector,
        MembraneFluidBoundary,
        ImmersedTransferConnector,
    ]:
        """Expose both owners, ERM, and the moving fluid boundary to the scheduler."""
        return (
            self.membrane,
            self.cortex,
            bindings.erm_connector,
            bindings.membrane_fluid_boundary,
            bindings.cortex_cytosol_transfer,
        )

    def snapshot_candidate(self, bindings: SurfaceStepBindings) -> None:
        """Snapshot every surface transaction participant without creating facade-owned state."""
        for participant in self.transaction_participants(bindings):
            participant.snapshot_candidate()

    def rollback(self, bindings: SurfaceStepBindings, accepted: wp.array) -> None:
        """Forward one device predicate to every surface transaction participant."""
        for participant in self.transaction_participants(bindings):
            participant.rollback(accepted)

    def commit_irreversible(
        self,
        bindings: SurfaceStepBindings,
        accepted: wp.array,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        """Commit participant kinetics only through the scheduler-owned final predicate."""
        for participant in self.transaction_participants(bindings):
            participant.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, bindings: SurfaceStepBindings, ledger: object) -> None:
        """Collect component, connector, and fluid-boundary ledgers without deciding acceptance."""
        self.membrane.accumulate_ledger(ledger)
        self.cortex.accumulate_ledger(ledger)
        bindings.erm_connector.accumulate_ledger(ledger)
        bindings.membrane_fluid_boundary.accumulate_ledger(ledger)
        bindings.cortex_cytosol_transfer.accumulate_ledger(ledger)


class _LegacyMembraneCompartment(Protocol):
    """Narrow structural subset of the landed MembraneCompartment used by the adapter."""

    node_off: int
    n_verts: int
    n_erm: int

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch existing membrane Helfrich/area kernels."""


@dataclass(frozen=True, slots=True)
class LegacyMembraneCompartmentAdapter:
    """Reuse a local, ERM-free ``MembraneCompartment`` as membrane mechanics only.

    Build the legacy compartment with ``node_off=0`` and no cortex candidates, then give this adapter to a
    ``SurfaceComponentStateOwner``.  The adapter deliberately refuses a monolithic/global offset or embedded
    ERM population: accepting either would collapse the two state owners or steal graph-owned connector state.
    """

    compartment: _LegacyMembraneCompartment

    def __post_init__(self) -> None:
        if int(self.compartment.node_off) != 0:
            raise ValueError("legacy membrane adapter requires a local array with node_off == 0")
        if int(self.compartment.n_erm) != 0:
            raise ValueError("legacy membrane adapter cannot own ERM state; register it in ConnectorGraph")
        if int(self.compartment.n_verts) <= 0:
            raise ValueError("legacy membrane compartment must contain surface vertices")

    @property
    def n_vertices(self) -> int:
        """Number of vertices required by the wrapped local membrane law."""
        return int(self.compartment.n_verts)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Delegate to the existing Warp Helfrich/area force implementation."""
        self.compartment.accumulate(pos, force)


@dataclass(frozen=True, slots=True)
class CortexFilamentMechanics:
    """Cortex mechanics delegate binding the native ``ff`` link + Cytosim-bending Warp kernels.

    This is the cortex-side counterpart of :class:`LegacyMembraneCompartmentAdapter`.  It launches the
    already-audited ``ff.network_warp.link_spring_kernel`` (Hookean axial backbone + crosslinker springs) and
    ``ff.forces_warp.cytosim_bending_kernel`` (NF2007 discrete bending) against the *cortex-owned* position/force
    arrays supplied by the :class:`SurfaceComponentStateOwner`.  It never allocates or touches membrane arrays,
    never concatenates cortex nodes into a whole-cell array, and never reads authoritative state on the host.

    The connectivity/parameter tables (``links_d``, ``link_k_d``, ``link_r0_d``, ``bend_triples_d``,
    ``bend_alpha_d``) are cortex-local device state.  Because SF/protrusion filaments are *separate* components,
    this table must contain only the cortex's own filament segments and crosslinks; joints to other components go
    through the connector graph, never through a shared row here.

    ``link_kernel``/``bending_kernel`` and ``launch`` are injected so the CUDA-unit structural test can substitute
    a recorder on a host without CUDA.  :meth:`bind_native` wires the real kernels for production.
    """

    device: str
    n_vertices: int
    links_d: wp.array
    link_k_d: wp.array
    link_r0_d: wp.array
    bend_triples_d: wp.array
    bend_alpha_d: wp.array
    link_kernel: object
    bending_kernel: object
    launch: Callable[..., object] = wp.launch
    n_links: int = field(init=False)
    n_triples: int = field(init=False)

    def __post_init__(self) -> None:
        if int(self.n_vertices) <= 0:
            raise ValueError("cortex mechanics must own at least one vertex")
        n_links = _validate_topology_pair_array(self.links_d, label="cortex.links_d", cols=2)
        _validate_scalar_f64_array(self.link_k_d, label="cortex.link_k_d", rows=n_links)
        _validate_scalar_f64_array(self.link_r0_d, label="cortex.link_r0_d", rows=n_links)
        n_triples = _validate_topology_pair_array(self.bend_triples_d, label="cortex.bend_triples_d", cols=3)
        _validate_scalar_f64_array(self.bend_alpha_d, label="cortex.bend_alpha_d", rows=n_triples)
        if n_links == 0 and n_triples == 0:
            raise ValueError("cortex mechanics must bind at least one link or bending triple")
        object.__setattr__(self, "n_links", n_links)
        object.__setattr__(self, "n_triples", n_triples)

    @classmethod
    def bind_native(
        cls,
        *,
        device: str,
        n_vertices: int,
        links_d: wp.array,
        link_k_d: wp.array,
        link_r0_d: wp.array,
        bend_triples_d: wp.array,
        bend_alpha_d: wp.array,
        launch: Callable[..., object] = wp.launch,
    ) -> "CortexFilamentMechanics":
        """Wire the production ``ff`` link/bending kernels (imported lazily to keep the seam light)."""
        from aleph.laws.forces_warp import cytosim_bending_kernel
        from aleph.laws.network_warp import link_spring_kernel

        return cls(
            device=device,
            n_vertices=int(n_vertices),
            links_d=links_d,
            link_k_d=link_k_d,
            link_r0_d=link_r0_d,
            bend_triples_d=bend_triples_d,
            bend_alpha_d=bend_alpha_d,
            link_kernel=link_spring_kernel,
            bending_kernel=cytosim_bending_kernel,
            launch=launch,
        )

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the cortex axial/crosslink and bending kernels onto the cortex force array.

        The two kernels ``wp.atomic_add`` into ``force``; the caller owns force zeroing.  Argument order matches
        the existing ``ff.ecm_mechanics`` launch sites so the reused kernels are called with their published
        signatures and never a re-derived variant.
        """
        if self.n_links:
            self.launch(
                self.link_kernel,
                dim=self.n_links,
                inputs=[pos, self.links_d, self.link_k_d, self.link_r0_d, force],
                device=self.device,
            )
        if self.n_triples:
            self.launch(
                self.bending_kernel,
                dim=self.n_triples,
                inputs=[pos, self.bend_triples_d, self.bend_alpha_d, force],
                device=self.device,
            )


@dataclass(frozen=True, slots=True)
class MembranePressureFluidBoundary:
    """``membrane_cytosol_boundary`` delegate binding the live membrane pressure-traction kernel.

    Wraps :class:`ac.cell.membrane_pressure.MembranePressureTraction`, which launches
    ``membrane_pressure_traction_kernel`` — the weak linear-triangle form of ``(p_inside - p_ext) n`` sampled from
    the graph-owned cytosol ``FieldGrid`` on the live membrane triangles.  This is the mechanistic surface term
    whose spherical continuum balance is Young--Laplace ``gamma = DeltaP*R/2``; it is the fluid boundary load, not
    a cortex body force, so binding it here keeps pressure off the cortex skeleton and prevents double counting.

    It is a *non-kinetic* field coupler: the authoritative pressure field belongs to the cytosol component, so this
    adapter owns no reversible per-candidate state and its snapshot/rollback/commit hooks are structural no-ops.
    The one device diagnostic it does hold (unresolved-face count) is monotone by contract — force assembly never
    erases an earlier failure — so it is deliberately not rewound on rollback.
    """

    traction: object
    name: str = MEMBRANE_FLUID_BOUNDARY
    component_a: str = MEMBRANE_COMPONENT
    component_b: str = CYTOSOL_COMPONENT

    def __post_init__(self) -> None:
        if self.name != MEMBRANE_FLUID_BOUNDARY:
            raise ValueError(f"pressure boundary must register as {MEMBRANE_FLUID_BOUNDARY!r}")
        if {self.component_a, self.component_b} != {MEMBRANE_COMPONENT, CYTOSOL_COMPONENT}:
            raise ValueError("pressure boundary must join membrane and cytosol")
        if not callable(getattr(self.traction, "accumulate", None)):
            raise TypeError("pressure boundary requires a traction primitive exposing accumulate(pos, force)")

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        """Add live spatial pressure traction to the membrane-owned force array (no host readback)."""
        self.traction.accumulate(pos, force)

    def snapshot_candidate(self) -> None:
        """No reversible boundary state is owned here; the pressure field belongs to cytosol."""

    def rollback(self, accepted: wp.array) -> None:
        """No-op: the monotone unresolved-face diagnostic is intentionally not rewound."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """No-op: pressure traction is a non-kinetic field coupler."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Optional device diagnostic hand-off; the coupler decides no acceptance."""
        sink = getattr(ledger, "accept_membrane_pressure_boundary", None)
        if callable(sink):
            sink(self.traction)


__all__ = [
    "CORTEX_COMPONENT",
    "CYTOSOL_COMPONENT",
    "CortexFilamentMechanics",
    "ERM_CONNECTOR",
    "LegacyMembraneCompartmentAdapter",
    "MEMBRANE_COMPONENT",
    "MEMBRANE_FLUID_BOUNDARY",
    "MembraneFluidBoundary",
    "MembranePressureFluidBoundary",
    "SurfaceBody",
    "SurfaceComponentStateOwner",
    "SurfacePairConnector",
    "SurfacePorousEndpointView",
    "SurfaceStepBindings",
    "SURFACE_POROUS_TRANSFER",
]
