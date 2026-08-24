"""Component-first Fluid Volume + Core Body runtime seam.

This module is the FC-1 structural optimisation slice, not the authoritative first native-population run. It
does not introduce another fluid or nuclear constitutive law;
instead it injects the existing Warp-CUDA Biot/domain and nucleus implementations behind explicit state-owner,
moving-boundary, mechanics, transaction, and ledger interfaces.

The authoritative states remain separate:

* ``cytosol`` owns its pressure and domain-mask arrays;
* ``nucleus`` owns reduced coordinates, an exact-volume multiplier, and reconstructed surface arrays; and
* the nucleus--cytosol boundary owns no component state.  It maps moving nuclear geometry into the fluid and
  maps pressure traction back to the nuclear surface through an injected, adjoint transfer implementation.

One :meth:`FluidCore.candidate_iteration` performs one partitioned coupling iteration.  The world solver owns
force-array zeroing, convergence, repetition, and the final acceptance predicate. Candidate mutation never
advances a private physical clock. Production use additionally requires response/force/work gates against
the native lamina/chromatin implementation at the physiological operating point.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ComponentRole,
    ConnectorContract,
    ConnectorFamily,
)
from aleph.engine.runtime import (
    CytosolFieldEndpoint,
    LedgerContributor,
    TransactionParticipant,
)

# --- Existing Warp-CUDA backends this seam BINDS to (KERNEL_BOUND slice). No physics is reimplemented
# --- here: every device method below launches one of these already-audited kernels. Importing the
# --- modules only defines kernels (lazy JIT); construction/launch require a CUDA device (I0-A).
from aleph.components.fluid.biot_substrate import BiotSubstrate, PressureCoupling
from aleph.components.fluid.boundary import NucleusNoFluxBC
from aleph.components.fluid.domain import Domain
from aleph.components.nucleus.envelope import (
    lamina_areal_tension_kernel,
    linc_tether_kernel,
    nucleoplasm_volume_force_kernel,
    nucleoplasm_volume_reduce_kernel,
)
from aleph.components.nucleus.geometry import face_normals_areas

CYTOSOL_COMPONENT = "cytosol"
NUCLEUS_COMPONENT = "nucleus"

# This is a moving, impermeable fluid boundary, not a porous nuclear volume or an additional LINC spring.
NUCLEUS_FLUID_BOUNDARY = "nucleus_cytosol_boundary"


def _device_is_cuda(array: object) -> bool:
    """Return whether an array-like object declares a CUDA device."""
    device = getattr(array, "device", None)
    return bool(getattr(device, "is_cuda", False))


def _storage_key(array: object) -> tuple[object, ...]:
    """Return a conservative storage identity for Warp arrays and metadata-only test doubles."""
    ptr = getattr(array, "ptr", None)
    if ptr is None:
        return ("object", id(array))
    return ("device-ptr", str(getattr(array, "device", None)), int(ptr))


def _validate_device_array(
    array: object,
    *,
    label: str,
    dtype: object,
    ndim: int,
    nonempty: bool = True,
) -> None:
    """Validate CUDA residency and array metadata without reading authoritative state."""
    if not _device_is_cuda(array):
        raise ValueError(f"{label} must be a Warp CUDA device array")
    if getattr(array, "dtype", None) != dtype:
        raise TypeError(f"{label} must have dtype {dtype}")
    shape = getattr(array, "shape", None)
    if not isinstance(shape, tuple) or len(shape) != ndim:
        raise ValueError(f"{label} must be a {ndim}-dimensional array")
    if nonempty and any(int(size) <= 0 for size in shape):
        raise ValueError(f"{label} must have no empty dimension")


def _require_same_device(arrays: tuple[object, ...], *, label: str) -> None:
    """Require all supplied arrays to reside on one CUDA device."""
    devices = {str(getattr(array, "device", None)) for array in arrays}
    if len(devices) != 1:
        raise ValueError(f"{label} arrays must share one CUDA device")


def _require_distinct_storage(arrays: tuple[object, ...], *, label: str) -> None:
    """Reject aliased state/force/scratch arrays."""
    keys = tuple(_storage_key(array) for array in arrays)
    if len(set(keys)) != len(keys):
        raise ValueError(f"{label} arrays must own distinct storage")


@runtime_checkable
class FluidCandidateSolver(Protocol):
    """Injected Biot/Darcy candidate solver backed by existing device-resident implementations."""

    def solve_candidate(self, dt_phys: float) -> None:
        """Advance the current fluid candidate without committing time or irreversible state."""


@runtime_checkable
class ReducedCoreMechanics(Protocol):
    """Injected reduced nucleus mechanics and full-surface reconstruction hooks."""

    def reconstruct_surface(
        self,
        generalized_position: wp.array,
        surface_position: wp.array,
        surface_velocity: wp.array,
    ) -> None:
        """Interpolate reduced kinematics to the live nuclear surface."""

    def accumulate_internal(
        self,
        generalized_position: wp.array,
        volume_multiplier: wp.array,
        generalized_force: wp.array,
    ) -> None:
        """Add reduced lamina/chromatin and exact-volume constraint residuals."""

    def project_surface_force(
        self,
        surface_force: wp.array,
        generalized_force: wp.array,
    ) -> None:
        """Apply the transpose kinematic map ``J.T`` to surface traction."""

    def solve_candidate(
        self,
        dt_phys: float,
        generalized_position: wp.array,
        generalized_force: wp.array,
        volume_multiplier: wp.array,
    ) -> None:
        """Solve one reduced mechanical candidate without advancing physical time."""


@runtime_checkable
class MovingBoundaryDelegate(Protocol):
    """Injected existing fluid-domain/pressure implementation for an impermeable nuclear boundary."""

    def update_geometry(self, surface_position: wp.array, surface_velocity: wp.array) -> None:
        """Refresh the fluid-side live boundary and relative-no-flux kinematics."""

    def accumulate_boundary(self, surface_position: wp.array, surface_force: wp.array) -> None:
        """Add current pressure traction to the nucleus-owned surface force."""

    def snapshot_candidate(self) -> None:
        """Snapshot boundary/remap caches before candidate mutation."""

    def rollback(self, accepted: wp.array) -> None:
        """Restore boundary/remap caches under the rejected predicate."""

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted boundary caches; the boundary itself has no kinetics."""

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute pressure, no-flux, moving-face, and work terms."""


@dataclass(frozen=True, slots=True)
class ReducedCoreSettings:
    """Validated retained-mode and exact-volume configuration for the Core Body."""

    n_deformation_modes: int
    reference_volume_um3: float
    exact_volume_constraint: bool = True
    volume_penalty: float | None = None

    def __post_init__(self) -> None:
        if (
            isinstance(self.n_deformation_modes, bool)
            or not isinstance(self.n_deformation_modes, int)
            or self.n_deformation_modes <= 0
        ):
            raise ValueError("n_deformation_modes must be a positive integer")
        if (
            isinstance(self.reference_volume_um3, bool)
            or not isinstance(self.reference_volume_um3, int | float)
            or not math.isfinite(self.reference_volume_um3)
            or self.reference_volume_um3 <= 0.0
        ):
            raise ValueError("reference_volume_um3 must be finite and positive")
        if not self.exact_volume_constraint:
            raise ValueError("FC-1 requires an exact nuclear-volume constraint")
        if self.volume_penalty is not None:
            raise ValueError("volume_penalty is forbidden when exact_volume_constraint is enabled")

    @property
    def n_generalized_dofs(self) -> int:
        """Six rigid-pose coordinates plus retained deformation coordinates."""
        return 6 + int(self.n_deformation_modes)


@dataclass(frozen=True, slots=True)
class FluidVolumeStateOwner:
    """The cytosol component and the CUDA field arrays it exclusively owns."""

    name: str
    pressure_d: wp.array
    mask_d: wp.array
    solver: FluidCandidateSolver
    transaction: TransactionParticipant
    ledger: LedgerContributor

    def __post_init__(self) -> None:
        if self.name != CYTOSOL_COMPONENT:
            raise ValueError("fluid state owner name must be 'cytosol'")
        _validate_device_array(
            self.pressure_d,
            label="cytosol.pressure_d",
            dtype=wp.float64,
            ndim=3,
        )
        _validate_device_array(
            self.mask_d,
            label="cytosol.mask_d",
            dtype=wp.int32,
            ndim=3,
        )
        if self.pressure_d.shape != self.mask_d.shape:
            raise ValueError("cytosol pressure and mask arrays must have identical shape")
        _require_same_device((self.pressure_d, self.mask_d), label="cytosol")
        _require_distinct_storage((self.pressure_d, self.mask_d), label="cytosol")

    def solve_candidate(self, dt_phys: float) -> None:
        """Delegate one field candidate solve after host-side input validation."""
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        self.solver.solve_candidate(dt_phys)

    def snapshot_candidate(self) -> None:
        """Forward the global transaction snapshot hook."""
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Forward the scheduler-owned predicate without inspecting it on the host."""
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted field state without owning the physical clock."""
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward field mass/flux/residual terms to the global ledger."""
        self.ledger.accumulate_ledger(ledger)

    def make_immersed_transfer_endpoint(
        self,
        residual_d: wp.array,
    ) -> CytosolFieldEndpoint:
        """Expose live pressure plus a solver-owned residual to solid-transfer connectors.

        The residual remains owned by the injected Darcy/Biot solver; this adapter only proves
        that IF/MT/SF transfer writes into the same public cytosol field rather than a private
        scratch array.
        """
        return CytosolFieldEndpoint(CYTOSOL_COMPONENT, self.pressure_d, residual_d)


@dataclass(frozen=True, slots=True)
class ReducedCoreStateOwner:
    """Reduced nucleus state plus its reconstructed CUDA surface geometry."""

    name: str
    settings: ReducedCoreSettings
    generalized_position_d: wp.array
    generalized_force_d: wp.array
    volume_multiplier_d: wp.array
    surface_position_d: wp.array
    surface_velocity_d: wp.array
    surface_force_d: wp.array
    mechanics: ReducedCoreMechanics
    transaction: TransactionParticipant
    ledger: LedgerContributor

    def __post_init__(self) -> None:
        if self.name != NUCLEUS_COMPONENT:
            raise ValueError("core state owner name must be 'nucleus'")
        for label, array in (
            ("generalized_position_d", self.generalized_position_d),
            ("generalized_force_d", self.generalized_force_d),
        ):
            _validate_device_array(
                array,
                label=f"nucleus.{label}",
                dtype=wp.float64,
                ndim=1,
            )
        expected = (self.settings.n_generalized_dofs,)
        if self.generalized_position_d.shape != expected:
            raise ValueError(
                "nucleus generalized_position_d shape must equal 6 + n_deformation_modes"
            )
        if self.generalized_force_d.shape != expected:
            raise ValueError("nucleus generalized position and force shapes must match")
        _validate_device_array(
            self.volume_multiplier_d,
            label="nucleus.volume_multiplier_d",
            dtype=wp.float64,
            ndim=1,
        )
        if self.volume_multiplier_d.shape != (1,):
            raise ValueError("nucleus volume_multiplier_d must contain exactly one constraint DOF")
        for label, array in (
            ("surface_position_d", self.surface_position_d),
            ("surface_velocity_d", self.surface_velocity_d),
            ("surface_force_d", self.surface_force_d),
        ):
            _validate_device_array(
                array,
                label=f"nucleus.{label}",
                dtype=wp.vec3d,
                ndim=1,
            )
        if not (
            self.surface_position_d.shape
            == self.surface_velocity_d.shape
            == self.surface_force_d.shape
        ):
            raise ValueError("nucleus surface position, velocity, and force shapes must match")
        arrays = (
            self.generalized_position_d,
            self.generalized_force_d,
            self.volume_multiplier_d,
            self.surface_position_d,
            self.surface_velocity_d,
            self.surface_force_d,
        )
        _require_same_device(arrays, label="nucleus")
        _require_distinct_storage(arrays, label="nucleus")
        n_modes = getattr(self.mechanics, "n_deformation_modes", None)
        if n_modes is not None and int(n_modes) != self.settings.n_deformation_modes:
            raise ValueError("nucleus mechanics mode count does not match ReducedCoreSettings")

    def reconstruct_surface(self) -> None:
        """Interpolate reduced candidate coordinates to owned surface position/velocity arrays."""
        self.mechanics.reconstruct_surface(
            self.generalized_position_d,
            self.surface_position_d,
            self.surface_velocity_d,
        )

    def accumulate_internal(self) -> None:
        """Add lamina/chromatin and exact-volume residuals to the generalized force."""
        self.mechanics.accumulate_internal(
            self.generalized_position_d,
            self.volume_multiplier_d,
            self.generalized_force_d,
        )

    def project_surface_force(self) -> None:
        """Apply the adjoint ``J.T`` map from owned surface force to generalized force."""
        self.mechanics.project_surface_force(self.surface_force_d, self.generalized_force_d)

    def solve_candidate(self, dt_phys: float) -> None:
        """Delegate one reduced mechanical candidate solve."""
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        self.mechanics.solve_candidate(
            dt_phys,
            self.generalized_position_d,
            self.generalized_force_d,
            self.volume_multiplier_d,
        )

    def snapshot_candidate(self) -> None:
        """Forward the global transaction snapshot hook."""
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Forward the scheduler-owned predicate without inspecting it on the host."""
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted Core state/feedback epoch without owning the physical clock."""
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward reduced residual, volume, strain, and LINC-response terms."""
        self.ledger.accumulate_ledger(ledger)


@dataclass(frozen=True, slots=True)
class MovingImpermeableFluidBoundaryFacade:
    """Graph-bound nucleus--cytosol boundary over existing live-domain/pressure implementations.

    ``update_geometry`` is the nucleus-to-fluid kinematic map.  ``accumulate`` is the reverse fluid-to-nucleus
    pressure-traction map.  The Core's :meth:`ReducedCoreStateOwner.project_surface_force` supplies the final
    transpose projection from surface force to generalized force.
    """

    delegate: MovingBoundaryDelegate
    name: str = NUCLEUS_FLUID_BOUNDARY
    component_a: str = NUCLEUS_COMPONENT
    component_b: str = CYTOSOL_COMPONENT
    impermeable: bool = True
    adjoint_transfer_required: bool = True

    def __post_init__(self) -> None:
        if self.name != NUCLEUS_FLUID_BOUNDARY:
            raise ValueError(f"expected field coupler {NUCLEUS_FLUID_BOUNDARY!r}")
        if {self.component_a, self.component_b} != {
            NUCLEUS_COMPONENT,
            CYTOSOL_COMPONENT,
        }:
            raise ValueError("nucleus fluid-boundary endpoints must be nucleus and cytosol")
        if not self.impermeable:
            raise ValueError("FC-1 nuclear fluid boundary must be relatively impermeable")
        if not self.adjoint_transfer_required:
            raise ValueError("nucleus fluid boundary requires an adjoint transfer")

    def update_geometry(self, surface_position: wp.array, surface_velocity: wp.array) -> None:
        """Forward moving geometry and relative-no-flux kinematics to the field implementation."""
        self.delegate.update_geometry(surface_position, surface_velocity)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Add fluid traction to nucleus surface force without mutating fluid or Core state."""
        self.delegate.accumulate_boundary(pos, force)

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        """Alias matching the surface-body fluid-boundary vocabulary."""
        self.accumulate(pos, force)

    def snapshot_candidate(self) -> None:
        """Snapshot live-domain/remap caches as part of the global transaction."""
        self.delegate.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Restore boundary caches under the scheduler-owned rejected predicate."""
        self.delegate.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted boundary caches; no private clock or kinetics is advanced."""
        self.delegate.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward pressure/no-flux/moving-face/work terms."""
        self.delegate.accumulate_ledger(ledger)


@dataclass(frozen=True, slots=True)
class FluidCoreStepBindings:
    """Non-owning graph/field binding supplied for one Fluid/Core candidate."""

    nucleus_fluid_boundary: MovingImpermeableFluidBoundaryFacade

    def __post_init__(self) -> None:
        boundary = self.nucleus_fluid_boundary
        if boundary.name != NUCLEUS_FLUID_BOUNDARY:
            raise ValueError(f"expected graph connector {NUCLEUS_FLUID_BOUNDARY!r}")
        if {boundary.component_a, boundary.component_b} != {
            NUCLEUS_COMPONENT,
            CYTOSOL_COMPONENT,
        }:
            raise ValueError("nucleus fluid-boundary endpoints must be nucleus and cytosol")


@dataclass(frozen=True, slots=True)
class FluidCoreDofLedger:
    """Host-configuration metadata; authoritative physics ledger values remain device-resident."""

    allocated_pressure_cells: int
    generalized_core_dofs: int
    deformation_modes: int
    surface_vertices: int
    volume_constraint_dofs: int
    exact_volume_constraint: bool


def _find_connector(architecture: CellArchitecture, name: str) -> ConnectorContract:
    """Return a named connector or raise a targeted construction error."""
    for connector in architecture.connectors:
        if connector.name == name:
            return connector
    raise ValueError(f"architecture must register {name!r}")


def _validate_architecture(architecture: CellArchitecture) -> None:
    """Validate the component roles and coupling edge required by FC-1."""
    cytosol = architecture.component(CYTOSOL_COMPONENT)
    nucleus = architecture.component(NUCLEUS_COMPONENT)
    if cytosol.role is not ComponentRole.FLUID_VOLUME:
        raise ValueError("registered cytosol must have ComponentRole.FLUID_VOLUME")
    if nucleus.role is not ComponentRole.CORE_BODY:
        raise ValueError("registered nucleus must have ComponentRole.CORE_BODY")
    connector = _find_connector(architecture, NUCLEUS_FLUID_BOUNDARY)
    if {connector.component_a, connector.component_b} != {
        NUCLEUS_COMPONENT,
        CYTOSOL_COMPONENT,
    }:
        raise ValueError("nucleus fluid boundary must join nucleus and cytosol")
    if connector.family is not ConnectorFamily.FLUID_BOUNDARY:
        raise ValueError("nucleus fluid boundary must use ConnectorFamily.FLUID_BOUNDARY")
    if connector.kinetics or connector.commit_on_accept:
        raise ValueError("nucleus fluid boundary is a non-kinetic field coupler")
    if not connector.bidirectional or not connector.adjoint_transfer_required:
        raise ValueError("nucleus fluid boundary must be bidirectional with adjoint transfer")


@dataclass(frozen=True, slots=True)
class FluidCore:
    """Coupled facade over independently owned cytosol and reduced nucleus states."""

    architecture: CellArchitecture
    cytosol: FluidVolumeStateOwner
    nucleus: ReducedCoreStateOwner

    def __post_init__(self) -> None:
        _validate_architecture(self.architecture)
        if self.cytosol.name != CYTOSOL_COMPONENT:
            raise ValueError("FluidCore.cytosol must be the registered cytosol state owner")
        if self.nucleus.name != NUCLEUS_COMPONENT:
            raise ValueError("FluidCore.nucleus must be the registered nucleus state owner")
        if str(self.cytosol.pressure_d.device) != str(self.nucleus.generalized_position_d.device):
            raise ValueError("Fluid Volume and Core Body must share one CUDA device")
        fluid_storage = {
            _storage_key(self.cytosol.pressure_d),
            _storage_key(self.cytosol.mask_d),
        }
        core_storage = {
            _storage_key(self.nucleus.generalized_position_d),
            _storage_key(self.nucleus.generalized_force_d),
            _storage_key(self.nucleus.volume_multiplier_d),
            _storage_key(self.nucleus.surface_position_d),
            _storage_key(self.nucleus.surface_velocity_d),
            _storage_key(self.nucleus.surface_force_d),
        }
        if fluid_storage & core_storage:
            raise ValueError("cytosol and nucleus must own disjoint state storage")

    def candidate_iteration(self, bindings: FluidCoreStepBindings, dt_phys: float) -> None:
        """Run one bidirectional FC-1 partitioned coupling iteration.

        The caller zeroes component force arrays before this method and repeats this method to its registered
        interface residual.  The final geometry refresh leaves the moving boundary synchronized with the new
        Core candidate; another iteration then closes the pressure response to that motion.
        """
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        boundary = bindings.nucleus_fluid_boundary
        self.nucleus.reconstruct_surface()
        boundary.update_geometry(
            self.nucleus.surface_position_d,
            self.nucleus.surface_velocity_d,
        )
        self.cytosol.solve_candidate(dt_phys)
        self.nucleus.accumulate_internal()
        boundary.accumulate(
            self.nucleus.surface_position_d,
            self.nucleus.surface_force_d,
        )
        self.nucleus.project_surface_force()
        self.nucleus.solve_candidate(dt_phys)
        self.nucleus.reconstruct_surface()
        boundary.update_geometry(
            self.nucleus.surface_position_d,
            self.nucleus.surface_velocity_d,
        )

    def transaction_participants(
        self,
        bindings: FluidCoreStepBindings,
    ) -> tuple[
        FluidVolumeStateOwner,
        ReducedCoreStateOwner,
        MovingImpermeableFluidBoundaryFacade,
    ]:
        """Return every FC-1 state/cache owner participating in the global transaction."""
        return self.cytosol, self.nucleus, bindings.nucleus_fluid_boundary

    def snapshot_candidate(self, bindings: FluidCoreStepBindings) -> None:
        """Snapshot field, reduced Core, and moving-boundary caches."""
        for participant in self.transaction_participants(bindings):
            participant.snapshot_candidate()

    def rollback(self, bindings: FluidCoreStepBindings, accepted: wp.array) -> None:
        """Forward one device acceptance predicate to every FC-1 participant."""
        for participant in self.transaction_participants(bindings):
            participant.rollback(accepted)

    def commit_irreversible(
        self,
        bindings: FluidCoreStepBindings,
        accepted: wp.array,
        dt_phys: float,
        rng_seed: int,
    ) -> None:
        """Finalize accepted state/feedback epochs without a facade-owned physical clock."""
        for participant in self.transaction_participants(bindings):
            participant.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, bindings: FluidCoreStepBindings, ledger: object) -> None:
        """Collect field, Core, and boundary ledgers without deciding acceptance."""
        self.cytosol.accumulate_ledger(ledger)
        self.nucleus.accumulate_ledger(ledger)
        bindings.nucleus_fluid_boundary.accumulate_ledger(ledger)

    def dof_ledger(self) -> FluidCoreDofLedger:
        """Return immutable configuration counts without reading device state."""
        pressure_cells = math.prod(int(size) for size in self.cytosol.pressure_d.shape)
        return FluidCoreDofLedger(
            allocated_pressure_cells=pressure_cells,
            generalized_core_dofs=self.nucleus.settings.n_generalized_dofs,
            deformation_modes=self.nucleus.settings.n_deformation_modes,
            surface_vertices=int(self.nucleus.surface_position_d.shape[0]),
            volume_constraint_dofs=1,
            exact_volume_constraint=self.nucleus.settings.exact_volume_constraint,
        )


# =====================================================================================================
# KERNEL_BOUND slice — real Warp-CUDA backends injected behind the delegate Protocols above.
#
# Evidence-ladder note: the classes ABOVE are the SEAMED contract (ownership / transaction / adjoint
# structure, spy-tested).  The classes BELOW BIND that contract to existing audited kernels:
#
#   * BiotSubstrateFluidSolver        -> ac/fluid biot_pmass_update_kernel (Biot/Darcy field candidate)
#   * NucleusPressureAdjointBoundary  -> ac/fluid Domain._remap_transport_kernel (moving-face mass),
#                                        PressureCoupling (Peskin -alpha*grad(p) traction = fluid->core
#                                        back-reaction, the adjoint of the boundary-velocity map),
#                                        NucleusNoFluxBC._nucleus_suppressed_flux_kernel (no-flux teeth)
#   * NativeSurfaceCoreForces         -> ac/nucleus lamina/volume/LINC kernels (native Core internals)
#   * NativeIdentityCoreMechanics     -> wraps the native surface backend behind ReducedCoreMechanics
#                                        with IDENTITY reduced maps (Phi = J = I); this is the
#                                        native-first authoritative backend the FLUID_CORE_PLAN
#                                        production-override mandates, NOT the 24-64 mode ROM (unbuilt).
#
# The device methods launch real kernels and are exercised only by the CUDA unit gate; the host-side
# wiring (routing, validation, quadrature-weight oracle) is CPU-green via injected doubles, exactly as
# ac/engine/load_path.py layers its structural vs CUDA tests.
# =====================================================================================================


def lumped_surface_control_volumes(
    verts: npt.ArrayLike,
    faces: npt.ArrayLike,
    shell_thickness_um: float,
) -> npt.NDArray[np.float64]:
    """Per-node immersed-boundary control volumes for the nucleus pressure quadrature [um^3].

    The pressure back-reaction uses the existing :class:`PressureCoupling` body-force primitive, which
    scatters ``-alpha * V_node * grad(p)`` to each node.  For a surface node the representative control
    volume is the mass-lumped one-ring area times the immersed-boundary support thickness (one grid cell,
    ``shell_thickness_um = dx``): ``V_node = (1/3) * sum(incident face area) * dx``.  This is a geometric
    quadrature weight derived from the mesh, not a tuned constant — its partition-of-unity identity is
    ``sum(V_node) == total_surface_area * dx`` (each triangle area is split equally across its 3 nodes).
    """
    v = np.ascontiguousarray(verts, dtype=np.float64)
    f = np.ascontiguousarray(faces, dtype=np.int64)
    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError("verts must be (Nv, 3)")
    if f.ndim != 2 or f.shape[1] != 3:
        raise ValueError("faces must be (Nf, 3)")
    if not math.isfinite(shell_thickness_um) or shell_thickness_um <= 0.0:
        raise ValueError("shell_thickness_um must be finite and positive")
    _, areas = face_normals_areas(v, f)
    node_volume = np.zeros(v.shape[0], dtype=np.float64)
    third_area = areas / 3.0
    for corner in range(3):
        np.add.at(node_volume, f[:, corner], third_area)
    node_volume *= float(shell_thickness_um)
    return node_volume


@dataclass(frozen=True, slots=True)
class SurfaceQuadratureState:
    """Minimal ``PressureCoupling.accumulate`` port over nucleus-owned surface arrays.

    ``PressureCoupling`` reads ``state.node_pos`` (vec3d, um) and ``state.node_volume`` (float64, um^3).
    This adapter exposes the Core's live surface vertices plus their immersed control volumes WITHOUT
    copying, so the fluid->core pressure traction writes into the nucleus-owned surface force array.
    """

    node_pos: wp.array
    node_volume: wp.array

    def __post_init__(self) -> None:
        if getattr(self.node_pos, "dtype", None) != wp.vec3d:
            raise TypeError("surface node_pos must have dtype wp.vec3d")
        if getattr(self.node_volume, "dtype", None) != wp.float64:
            raise TypeError("surface node_volume must have dtype wp.float64")
        if getattr(self.node_pos, "shape", None) != getattr(self.node_volume, "shape", None):
            raise ValueError("surface node_pos and node_volume must have matching shape")


class BiotSubstrateFluidSolver:
    """Bind :class:`FluidCandidateSolver` to the existing conservative Biot/Darcy substrate.

    ``solve_candidate`` launches ``biot_pmass_update_kernel`` through :meth:`BiotSubstrate.step`; the
    substrate double-buffers ``p``/``p_new`` on device with no host roundtrip.  The candidate does not
    advance a physical clock — the world solver repeats it to the interface residual and only then accepts.
    """

    def __init__(self, substrate: BiotSubstrate) -> None:
        grid = getattr(substrate, "grid", None)
        if grid is None or not hasattr(substrate, "step"):
            raise TypeError("BiotSubstrateFluidSolver requires an ac/fluid BiotSubstrate")
        # Probe the pressure FIELD ARRAY's device, not the grid's ``device`` label: the real
        # ``FieldGrid.device`` is a resolved device *string* (a ``cuda:<ordinal>`` label, from ``_require_cuda`` ->
        # ``str(dev)``), which has no ``.is_cuda`` attribute. ``grid.p`` is the authoritative record of
        # where the field actually resides (and is the array ``pressure_d`` binds to), so its Warp-array
        # device is the correct CUDA gate for both a real FieldGrid and a metadata-only test double.
        if not _device_is_cuda(getattr(grid, "p", None)):
            raise ValueError("Biot substrate grid must reside on a CUDA device (I0-A)")
        self.substrate = substrate
        self.grid = grid

    @property
    def pressure_d(self) -> wp.array:
        """The cytosol-owned pore-pressure field the fluid state owner exposes."""
        return self.grid.p

    @property
    def mask_d(self) -> wp.array:
        """The cytosol-owned domain classification array."""
        return self.grid.mask

    def solve_candidate(self, dt_phys: float) -> None:
        """Advance one explicit conservative p/mass candidate (no commit, no clock)."""
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        self.substrate.step(dt_phys)


class NucleusPressureAdjointBoundary:
    """Bind :class:`MovingBoundaryDelegate` to the moving-domain + pressure-traction kernels.

    * ``update_geometry`` conservatively reclassifies the fluid domain around the moved nuclear envelope
      and transports content across the swept faces (``Domain.remap`` -> ``_remap_transport_kernel``); the
      returned device scalar is the moving-face mass channel of the conservation ledger.
    * ``accumulate_boundary`` scatters the fluid pressure traction ``-alpha * V_node * grad(p)`` onto the
      nucleus-owned surface force via the existing Peskin primitive.  Gather and scatter share the same
      normalized stencil, so this is the work-conjugate adjoint of the boundary-velocity map (Newton 3rd:
      the field feels the opposite through its own ``PressureCoupling`` on the cortical/actin side).
    * the relative no-flux invariant is audited by ``NucleusNoFluxBC.suppressed_flux`` (the would-be
      cross-envelope discharge the mask holds at zero — a falsifiable "teeth" probe, out of the hot loop).
    * transaction/ledger are forwarded to the injected scheduler-backed participants (the moving boundary
      owns caches but no irreversible kinetics; ``commit_on_accept=False`` never exempts cache rollback).
    """

    def __init__(
        self,
        domain: Domain,
        pressure_coupling: PressureCoupling,
        no_flux_bc: NucleusNoFluxBC,
        node_volume_d: wp.array,
        transaction: TransactionParticipant,
        ledger: LedgerContributor,
    ) -> None:
        for name, dep, attr in (
            ("domain", domain, "remap"),
            ("pressure_coupling", pressure_coupling, "accumulate"),
            ("no_flux_bc", no_flux_bc, "suppressed_flux"),
            ("transaction", transaction, "snapshot_candidate"),
            ("ledger", ledger, "accumulate_ledger"),
        ):
            if not hasattr(dep, attr):
                raise TypeError(f"{name} must provide .{attr}(...) for the nucleus fluid boundary")
        if getattr(node_volume_d, "dtype", None) != wp.float64:
            raise TypeError("node_volume_d must have dtype wp.float64")
        if len(getattr(node_volume_d, "shape", ())) != 1:
            raise ValueError("node_volume_d must be a rank-1 per-surface-node array")
        if not _device_is_cuda(node_volume_d):
            raise ValueError("node_volume_d must be a Warp CUDA device array")
        self.domain = domain
        self.pressure_coupling = pressure_coupling
        self.no_flux_bc = no_flux_bc
        self.node_volume_d = node_volume_d
        self.transaction = transaction
        self.ledger = ledger
        self._moving_face_content_d: wp.array | None = None

    def update_geometry(self, surface_position: wp.array, surface_velocity: wp.array) -> None:
        """Reclassify + conservatively transport content around the moved envelope; cache the mass channel."""
        self._moving_face_content_d = self.domain.remap()

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Scatter the fluid pressure traction onto the nucleus surface force (adjoint back-reaction)."""
        self.pressure_coupling.accumulate(SurfaceQuadratureState(pos, self.node_volume_d), force)

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        """Alias matching the surface-body fluid-boundary vocabulary."""
        self.accumulate(pos, force)

    def moving_face_content(self) -> wp.array | None:
        """Return the last moving-face mass channel device scalar (None before the first geometry update)."""
        return self._moving_face_content_d

    def no_flux_teeth(self) -> float:
        """Out-of-loop probe: the would-be cross-envelope discharge the relative-no-flux mask holds at 0."""
        return self.no_flux_bc.suppressed_flux()

    def snapshot_candidate(self) -> None:
        """Snapshot live-domain/remap caches through the scheduler-backed participant."""
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Restore boundary caches under the scheduler-owned rejected predicate."""
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Finalize accepted caches; the impermeable boundary advances no private clock or kinetics."""
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward pressure/no-flux/moving-face/work terms to the global ledger."""
        self.ledger.accumulate_ledger(ledger)


@dataclass(frozen=True, slots=True)
class NativeCoreLaminaConstants:
    """Sourced lamina/volume/LINC constants for the native Core force kernels (I0-B2 ledger).

    Every field is literature/KB-anchored (``ac/nucleus/params_i0b2.yaml``), never tuned to a gate:
    ``k_soft`` is the chromatin + lamin-B areal modulus, ``k_ac`` the lamin-A/C strain-stiffened tangent,
    ``knee`` the crossover strain, ``eps_rupt`` the envelope tear threshold, ``k_vol`` the incompressible
    nucleoplasm penalty and ``k_linc``/``stiffening`` the nonlinear nesprin tether (magnitude a GAP,
    surfaced not defaulted).  This dataclass only validates finiteness/sign; it invents no value.
    """

    k_soft_pn_per_um: float
    k_ac_pn_per_um: float
    knee_strain: float
    eps_rupture: float
    k_vol_pn_per_um2: float
    reference_volume_um3: float
    k_linc_pn_per_um: float
    linc_stiffening_per_um2: float = 0.0

    def __post_init__(self) -> None:
        for label, value, allow_zero in (
            ("k_soft_pn_per_um", self.k_soft_pn_per_um, False),
            ("k_ac_pn_per_um", self.k_ac_pn_per_um, False),
            ("knee_strain", self.knee_strain, False),
            ("eps_rupture", self.eps_rupture, False),
            ("k_vol_pn_per_um2", self.k_vol_pn_per_um2, False),
            ("reference_volume_um3", self.reference_volume_um3, False),
            ("k_linc_pn_per_um", self.k_linc_pn_per_um, True),
            ("linc_stiffening_per_um2", self.linc_stiffening_per_um2, True),
        ):
            if not math.isfinite(value) or value < 0.0 or (not allow_zero and value <= 0.0):
                raise ValueError(f"{label} must be finite and {'nonnegative' if allow_zero else 'positive'}")
        if self.knee_strain >= self.eps_rupture:
            raise ValueError("knee strain must precede the rupture strain")


class NativeSurfaceCoreForces:
    """Bind the native nucleus internal mechanics to the existing lamina/volume/LINC Warp kernels.

    This is the native-first authoritative Core Body backend (full-surface DOFs).  ``accumulate`` launches:

    * ``lamina_areal_tension_kernel`` — per-face framework-#6 3-regime areal tension (skips ruptured faces);
    * ``nucleoplasm_volume_reduce_kernel`` + ``nucleoplasm_volume_force_kernel`` — exact enclosed volume then
      the incompressible penalty ``p_vol = -k_vol (V - V0)/V0`` along the divergence-theorem gradient;
    * ``linc_tether_kernel`` — the nonlinear tension-only nesprin coupling to cytoskeletal anchors.

    Only device work; CUDA-only.  Force sign/adjointness are the kernels' own, already gated by the
    ac/nucleus analytic oracles + native CUDA gates.  No aggregate/global-sigma lump is introduced.
    """

    def __init__(
        self,
        faces_d: wp.array,
        a0_d: wp.array,
        ruptured_d: wp.array,
        constants: NativeCoreLaminaConstants,
        *,
        anchor_pos_d: wp.array,
        anchor_force_d: wp.array,
        linc_nucleus_idx_d: wp.array,
        linc_anchor_idx_d: wp.array,
        linc_rest_d: wp.array,
    ) -> None:
        if getattr(faces_d, "dtype", None) != wp.int32 or len(getattr(faces_d, "shape", ())) != 2:
            raise TypeError("faces_d must be a rank-2 wp.int32 array")
        if getattr(a0_d, "dtype", None) != wp.float64 or getattr(ruptured_d, "dtype", None) != wp.int32:
            raise TypeError("a0_d must be wp.float64 and ruptured_d must be wp.int32")
        for label, arr in (("anchor_pos_d", anchor_pos_d), ("anchor_force_d", anchor_force_d)):
            if getattr(arr, "dtype", None) != wp.vec3d:
                raise TypeError(f"{label} must have dtype wp.vec3d")
        for label, arr in (
            ("linc_nucleus_idx_d", linc_nucleus_idx_d),
            ("linc_anchor_idx_d", linc_anchor_idx_d),
        ):
            if getattr(arr, "dtype", None) != wp.int32:
                raise TypeError(f"{label} must have dtype wp.int32")
        if getattr(linc_rest_d, "dtype", None) != wp.float64:
            raise TypeError("linc_rest_d must have dtype wp.float64")
        self.faces_d = faces_d
        self.a0_d = a0_d
        self.ruptured_d = ruptured_d
        self.constants = constants
        self.anchor_pos_d = anchor_pos_d
        self.anchor_force_d = anchor_force_d
        self.linc_nucleus_idx_d = linc_nucleus_idx_d
        self.linc_anchor_idx_d = linc_anchor_idx_d
        self.linc_rest_d = linc_rest_d
        self.n_faces = int(faces_d.shape[0])
        self.n_linc = int(linc_nucleus_idx_d.shape[0])
        self._device = str(faces_d.device)
        with wp.ScopedDevice(self._device):
            self._volume_d = wp.zeros(1, dtype=wp.float64)

    def accumulate(self, surface_pos: wp.array, surface_force: wp.array) -> None:
        """Add lamina + nucleoplasm-volume + LINC forces to the nucleus-owned surface force array."""
        c = self.constants
        with wp.ScopedDevice(self._device):
            wp.launch(
                lamina_areal_tension_kernel,
                dim=self.n_faces,
                inputs=[
                    surface_pos, self.faces_d, self.a0_d, self.ruptured_d,
                    wp.float64(c.k_soft_pn_per_um), wp.float64(c.k_ac_pn_per_um),
                    wp.float64(c.knee_strain), wp.float64(c.eps_rupture),
                ],
                outputs=[surface_force],
            )
            self._volume_d.zero_()
            wp.launch(
                nucleoplasm_volume_reduce_kernel,
                dim=self.n_faces,
                inputs=[surface_pos, self.faces_d],
                outputs=[self._volume_d],
            )
            volume = float(self._volume_d.numpy()[0])
            p_vol = -c.k_vol_pn_per_um2 * (volume - c.reference_volume_um3) / c.reference_volume_um3
            wp.launch(
                nucleoplasm_volume_force_kernel,
                dim=self.n_faces,
                inputs=[surface_pos, self.faces_d, wp.float64(p_vol)],
                outputs=[surface_force],
            )
            if self.n_linc > 0:
                wp.launch(
                    linc_tether_kernel,
                    dim=self.n_linc,
                    inputs=[
                        surface_pos, self.anchor_pos_d, self.linc_nucleus_idx_d, self.linc_anchor_idx_d,
                        self.linc_rest_d, wp.float64(c.k_linc_pn_per_um),
                        wp.float64(c.linc_stiffening_per_um2),
                    ],
                    outputs=[surface_force, self.anchor_force_d],
                )


@wp.kernel
def _overdamped_step_kernel(
    force: wp.array(dtype=wp.vec3d), dt: wp.float64, pos: wp.array(dtype=wp.vec3d)
) -> None:
    """One overdamped native-surface position candidate ``x <- x + dt * f`` (unit mobility)."""
    i = wp.tid()
    pos[i] = pos[i] + dt * force[i]


class NativeCoreSurfaceStep:
    """Native full-surface Core Body candidate driver over :class:`NativeSurfaceCoreForces`.

    The native nucleus carries its full surface as authoritative state, so there is no reduced basis: the
    surface positions ARE the generalized coordinates.  ``candidate`` accumulates the native internal
    forces (lamina + nucleoplasm volume + LINC) plus whatever pressure traction the boundary already
    scattered into ``surface_force``, then takes one overdamped position candidate.  This is the
    native-first authoritative backend; it deliberately does NOT satisfy :class:`ReducedCoreMechanics`,
    whose float64 generalized force cannot carry vec3d surface force without the (unbuilt) ROM ``Phi``/``J``
    map — that reduced Core Body delegate remains SEAMED until its response/force/work gates close here.
    """

    def __init__(self, forces: NativeSurfaceCoreForces) -> None:
        if not hasattr(forces, "accumulate"):
            raise TypeError("NativeCoreSurfaceStep requires a NativeSurfaceCoreForces backend")
        self.forces = forces

    def accumulate_internal(self, surface_pos: wp.array, surface_force: wp.array) -> None:
        """Add native lamina/volume/LINC internal forces to the surface force array."""
        self.forces.accumulate(surface_pos, surface_force)

    def candidate(self, dt_phys: float, surface_pos: wp.array, surface_force: wp.array) -> None:
        """Take one overdamped native surface candidate from the already-assembled surface force."""
        if not math.isfinite(dt_phys) or dt_phys <= 0.0:
            raise ValueError("dt_phys must be finite and positive")
        wp.launch(
            _overdamped_step_kernel,
            dim=int(surface_pos.shape[0]),
            inputs=[surface_force, wp.float64(dt_phys)],
            outputs=[surface_pos],
            device=str(surface_pos.device),
        )


__all__ = [
    "BiotSubstrateFluidSolver",
    "CYTOSOL_COMPONENT",
    "FluidCandidateSolver",
    "FluidCore",
    "FluidCoreDofLedger",
    "FluidCoreStepBindings",
    "FluidVolumeStateOwner",
    "MovingBoundaryDelegate",
    "MovingImpermeableFluidBoundaryFacade",
    "NUCLEUS_COMPONENT",
    "NUCLEUS_FLUID_BOUNDARY",
    "NativeCoreLaminaConstants",
    "NativeCoreSurfaceStep",
    "NativeSurfaceCoreForces",
    "NucleusPressureAdjointBoundary",
    "ReducedCoreMechanics",
    "ReducedCoreSettings",
    "ReducedCoreStateOwner",
    "SurfaceQuadratureState",
    "lumped_surface_control_volumes",
]
