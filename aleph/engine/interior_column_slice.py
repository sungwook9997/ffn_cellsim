r"""GATE-B third dynamic-runtime vertical slice: the interior column membrane -> cortex -> cytosol -> nucleus.

The first two dynamic slices grew ACTIVE surface tension from binding EVENTS along the cortex:
:mod:`aleph.engine.cortex_motor_slice` (NMII head binding) and
:mod:`aleph.engine.erm_cortex_slice` (ERM membrane--cortex tether).  Both stay on the SURFACE.  This
slice drives :class:`~aleph.engine.transaction.CellTransaction` INWARD, over one radial column of the
composed cell graph::

    membrane ─(membrane_cytosol_boundary)─ cytosol ─(nucleus_cytosol_boundary)─ nucleus
        │                                    │
        └────────── cortex ─(surface_porous_transfer)─┘

so the four load-bearing interior compartments run TOGETHER as accepted-step participants.  It introduces NO
new physics: every mechanical kernel it launches is one already bound by ``ac/engine/fluid_core.py`` +
``ac/engine/cytosol_connected.py`` (the Biot p/mass field, the Peskin ``-alpha*V*grad(p)`` pressure traction,
the Kedem--Katchalsky membrane water flux, the moving-domain conservative remap) or by ``ac/nucleus`` (the
per-face lamina areal tension, the nucleoplasm incompressible-volume penalty, the nesprin LINC tether, the
device-gated envelope-rupture commit).  The slice is the COMPOSITION: it re-expresses the incumbent
single-owner coupled fluid step as one :class:`CellTransaction` whose participants snapshot / roll back /
commit-on-accept under ONE device-resident acceptance predicate, exactly as
:class:`~aleph.engine.erm_cortex_slice.ERMCortexSlice` does for the ERM pair.

Per accepted physical step the transaction:

1. **snapshots** every participant's authoritative reversible state as one candidate — the four surface/field
   owners (membrane + cortex surface positions, the Biot ``p``/``mask``/``div_vs``/``s_*`` field, the nucleus
   surface positions + per-face ``ruptured`` flag) and the three moving-boundary caches;
2. **proposes events** — a no-op: every interior-column coupling is a NON-kinetic field coupler
   (``kinetics=False``); the only irreversible kinetic is nucleus envelope rupture, committed in step 5;
3. runs the caller-owned **inner mechanical solve** — during which :meth:`InteriorColumnSlice.accumulate`
   scatters the coupling loads as bidirectional adjoint pairs (Newton's 3rd law): the cytosol pressure traction
   onto the cortex (``surface_porous_transfer``) and onto the membrane / nuclear envelope (the two
   ``*_cytosol_boundary`` connectors), plus the nucleus's own lamina/volume/LINC internal forces — while the
   Biot field relaxes to the moved boundaries and the compartments relax to convergence;
4. accumulates the coupled MASS / NO-FLUX / ADJOINT-WORK ledger channels (membrane permeation + nucleus
   swept content; the impermeable-nucleus no-flux teeth; the transpose pressure-traction work);
5. **commits** the accepted-gated nucleus envelope rupture on the converged geometry (the single irreversible
   interior-column kinetic) — or rolls the whole column back bit-exactly on a rejected step;
6. **advances the device event clock** iff the step was accepted.

No fixed prestress is imposed anywhere: the pressure the field carries and the tension the envelope develops
are the reactions of the coupled candidate the world solver relaxes.  Every CUDA op flows through an injected
launcher/copy or an injected fluid primitive, so the module is CPU-importable and the structural gates drive it
with recording doubles; the CUDA-lane builders assemble the real Warp-resident runtimes on the gbook A5000
(:mod:`aleph.scripts.ac_gate_b_interior_column_native`).

Sanity Gate:
    * ownership: membrane / cortex / cytosol / nucleus are four DISJOINT state-owners; the three fluid
      connectors are the ONLY interior mechanical couplings (co-location in the global array is never a
      connection — the split is logical, Newton's 3rd law still holds within it).
    * boundary/sign: every coupling is bidirectional + adjoint (the Peskin gather/scatter is a transpose
      pair, so no spurious work); the impermeable nucleus carries zero permeation (mask-structural no-flux)
      while the semipermeable membrane carries the only water flux.
    * conservation: the field + surface + ruptured state are snapshot and reject-gated bit-restored; the
      rupture commit and the clock advance are device-gated on the accepted predicate — a rejected step
      advances no interior state and no clock (no host ``bool()``/``numpy()`` on the predicate).
    * REQUIRED-PARAM: the nucleus lamina/volume/LINC moduli + the rupture strain are the open I0-B2 GAP-PI
      (``ac/nucleus/params_i0b2.yaml``); the osmotic driver ``Pi_0`` / ``L_p`` are the membrane-boundary GAP.
      The builders REQUIRE them (no convenient defaults); the native script prints a loud PROVISIONAL banner.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from math import isfinite
from typing import Any

import warp as wp

from aleph.engine.actor import CellActor
from aleph.engine.contracts import CellArchitecture, reference_cell_architecture
from aleph.engine.cytosol_connected import (
    MEMBRANE_FLUID_BOUNDARY,
    MovingSemipermeableFluidBoundaryFacade,
)
from aleph.engine.erm_cortex_connector import (
    _restore_float64_if_rejected_kernel,
    _restore_int32_if_rejected_kernel,
)
from aleph.engine.erm_cortex_slice import _restore_vec3_if_rejected_kernel
from aleph.engine.immersed_transfer import ImmersedPorousTransfer, build_immersed_porous_transfer
from aleph.engine.fluid_core import (
    CYTOSOL_COMPONENT,
    NUCLEUS_COMPONENT,
    NUCLEUS_FLUID_BOUNDARY,
    FluidVolumeStateOwner,
    MovingImpermeableFluidBoundaryFacade,
    SurfaceQuadratureState,
)
from aleph.engine.surface_body import (
    CORTEX_COMPONENT,
    MEMBRANE_COMPONENT,
    SURFACE_POROUS_TRANSFER,
    SurfacePorousEndpointView,
)
from aleph.engine.transaction import CellTransaction
from aleph.engine.world import CellWorldTransaction
from aleph.components.fluid.scheduler import (
    _conditional_restore_f64_3d_kernel,
    _conditional_restore_i32_3d_kernel,
)
from aleph.components.nucleus.envelope import conditional_rupture_update_kernel

__all__ = [
    "CORE_COMPONENT",
    "CoreSurfaceStateOwner",
    "CoreSurfaceTransaction",
    "FluidFieldTransaction",
    "InteriorColumnSlice",
    "NoOpSurfaceMechanics",
    "RecordingLedgerStub",
    "assert_disjoint_blocks",
    "build_core_surface_owner",
    "build_cortex_cytosol_transfer",
    "build_fluid_field_owner",
    "component_block_views",
    "build_interior_column_slice",
]

CORE_COMPONENT = NUCLEUS_COMPONENT

LaunchFn = Callable[..., None]
CopyFn = Callable[..., None]


def _default_launch(kernel: object, *, dim: object, inputs: list, outputs: list | None = None,
                    device: object | None = None) -> None:
    """Production launcher: forward to ``wp.launch`` on the owning CUDA device."""
    if outputs is None:
        wp.launch(kernel, dim=dim, inputs=inputs, device=device)
    else:
        wp.launch(kernel, dim=dim, inputs=inputs, outputs=outputs, device=device)


def _default_copy(dst: object, src: object) -> None:
    """Production D2D snapshot/restore: forward to ``wp.copy``."""
    wp.copy(dst, src)


def _require_param(name: str, value: float | None, *, gap: bool = False) -> float:
    """Return a supplied positive-finite param, or raise a REQUIRED-PARAM error naming the GAP."""
    suffix = " (KB/PI GAP — source it or surface to PI, never default)" if gap else ""
    if value is None or not isfinite(float(value)) or float(value) <= 0.0:
        raise ValueError(f"REQUIRED-PARAM {name!r} must be a supplied positive-finite value{suffix}; got {value!r}")
    return float(value)


# ── Card-5: each interior-column component ADDRESSES its own block of the incumbent's global arrays ──
#
# The incumbent's inner solve relaxes ONE global ``pos_d``, so a component cannot hold a private
# allocation without splitting that solve — which is new physics and out of the feature freeze.  What it
# CAN hold is its own contiguous block, and that is the difference between a component that owns a
# compartment and one that claims every node in the cell.  It is a strictly smaller claim than a private
# allocation and must never be reported as one (``STATE.md`` (c) 4).
def assert_disjoint_blocks(blocks: dict[str, tuple[int, int]], *, n_total: int) -> None:
    """Assert every named ``(node_off, n_nodes)`` block lies inside the global array and touches no other.

    This is where a wrong offset becomes an exception.  It is otherwise silent: an offset into the wrong
    compartment still indexes something, the kernel still launches, and the number it produces is wrong in
    a way no downstream gate reads as an error.

    Args:
        blocks: component name -> ``(node_off, n_nodes)`` in the global node array.
        n_total: the global array's node count.

    Raises:
        ValueError: If a block is empty, leaves the global array, or shares a node with another.
    """
    ordered = sorted(blocks.items(), key=lambda item: int(item[1][0]))
    for name, (node_off, n_nodes) in ordered:
        if int(n_nodes) <= 0:
            raise ValueError(f"component {name!r} must address at least one node; got n_nodes={n_nodes!r}")
        if int(node_off) < 0 or int(node_off) + int(n_nodes) > int(n_total):
            raise ValueError(
                f"component {name!r} block [{node_off}, {int(node_off) + int(n_nodes)}) falls outside the "
                f"global array of {int(n_total)} nodes — an offset into the wrong array still indexes "
                "something, which is why this is refused rather than clamped"
            )
    for (a_name, (a_off, a_n)), (b_name, (b_off, _b_n)) in zip(ordered, ordered[1:]):
        if int(a_off) + int(a_n) > int(b_off):
            raise ValueError(
                f"components {a_name!r} and {b_name!r} address overlapping nodes "
                f"([{a_off}, {int(a_off) + int(a_n)}) vs starting at {b_off}) — co-location in an array is "
                "NEVER a connection"
            )


def component_block_views(
    position_d: wp.array,
    force_d: wp.array,
    *,
    node_off: int,
    n_nodes: int,
) -> tuple[wp.array, wp.array]:
    """Return one component's ``(position, force)`` SLICE VIEWS of the incumbent's global arrays.

    **These are views, not private allocations.**  The returned arrays share their bytes with the global
    arrays the incumbent's inner solve relaxes; what they do NOT share is addressing — the view is exactly
    ``n_nodes`` long, so a kernel launched over it cannot reach another component's node, and the
    component's snapshot / reject-restore covers its own block instead of the whole cell.  A run artifact
    must record this as ``arrays_private: false`` with the block it addresses; reporting a view as a
    private allocation is the precise claim ``STATE.md`` (c) 4 refuses.

    Args:
        position_d / force_d: the built cell's global node arrays.
        node_off / n_nodes: the component's contiguous block, validated by :func:`assert_disjoint_blocks`.

    Returns:
        ``(position_view, force_view)``, each ``n_nodes`` long, addressing the same nodes in both.
    """
    stop = int(node_off) + int(n_nodes)
    return position_d[int(node_off):stop], force_d[int(node_off):stop]


class NoOpSurfaceMechanics:
    """No-op mechanics for a native owner whose forces the driver's own inner solve already assembles.

    In the native interior-column driver the membrane Helfrich/area, the cortex network, the nucleus
    lamina/volume/LINC, and the fluid pressure coupling are ALL assembled by the driver's ``_accumulate_all``
    (the incumbent single-owner coupled fluid step).  Re-adding any of them through the slice would DOUBLE
    COUNT, so an owner wired for that driver is given this no-op mechanics — exactly as the ERM native driver
    routes its membrane through ``_SliceMembraneAdapter`` and gives the surface owner a no-op mechanics.
    """

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """No-op: the driver's inner solve owns the compartment mechanics (no double count)."""


class RecordingLedgerStub:
    """Minimal :class:`LedgerContributor` for an interior-column participant (records the ledger handle)."""

    def __init__(self) -> None:
        self.ledger_calls: list[object] = []

    def accumulate_ledger(self, ledger: object) -> None:
        """Record the ledger handle without deciding acceptance (device reductions land on the native lane)."""
        self.ledger_calls.append(ledger)


# ── Cytosol Biot p/mass field: the concrete snapshot / reject-restore participant (KERNEL_BOUND). ──────
class FluidFieldTransaction:
    """Snapshot + reject-gated restore of the cytosol Biot field — the missing FLUID_VOLUME participant.

    ``ac/engine/fluid_core.FluidVolumeStateOwner`` forwards its transaction hooks to an INJECTED
    :class:`TransactionParticipant` that the seam did not previously construct — that concrete
    snapshot/restore owner is exactly what stood between the cytosol CUDA_UNIT and CONNECTED (the cytosol
    binding study, §(c) item 1).  This is it.  It D2D-snapshots the reversible field arrays the Biot candidate
    mutates and, on a REJECTED outer step, restores them IN-DEVICE via the already-audited
    ``ac/fluid/scheduler._conditional_restore_{f64,i32}_3d_kernel`` (the same kernels
    ``PhysicalScheduler.outer_step`` uses).  The pressure field is REVERSIBLE state (snapshot/rollback); the
    field owns no irreversible kinetic, so ``commit_irreversible`` is a no-op (the moving-boundary swept-content
    caches + the RNG epoch commit under the same predicate through the boundary participants and the clock).

    ``grid`` is read LIVE each call so the Biot double-buffer ping-pong (``grid.p, grid.p_new = ...``) is
    transparent: the reject-gated restore writes the snapshot back into whichever buffer ``grid.p`` currently
    names, by value, so the restored field content is bit-exact regardless of parity.
    """

    def __init__(
        self,
        *,
        grid: object,
        snap_p_d: wp.array,
        snap_mask_d: wp.array,
        snap_div_vs_d: wp.array,
        snap_s_membrane_d: wp.array,
        snap_s_total_d: wp.array,
        shape: tuple[int, ...],
        launch: LaunchFn = _default_launch,
        copy: CopyFn = _default_copy,
        device: object | None = None,
    ) -> None:
        self._grid = grid
        self._snap_p_d = snap_p_d
        self._snap_mask_d = snap_mask_d
        self._snap_div_vs_d = snap_div_vs_d
        self._snap_s_membrane_d = snap_s_membrane_d
        self._snap_s_total_d = snap_s_total_d
        self._shape = tuple(int(n) for n in shape)
        self._launch = launch
        self._copy = copy
        self._device = device
        self.launched_kernels: list[object] = []
        self.commit_calls: list[tuple[Any, ...]] = []

    def _f64_pairs(self) -> tuple[tuple[wp.array, wp.array], ...]:
        """(live, snap) float64 field pairs, read live so the ``p``/``p_new`` swap is transparent."""
        return (
            (self._grid.p, self._snap_p_d),
            (self._grid.div_vs, self._snap_div_vs_d),
            (self._grid.s_membrane, self._snap_s_membrane_d),
            (self._grid.s_total, self._snap_s_total_d),
        )

    def snapshot_candidate(self) -> None:
        """D2D-snapshot the reversible Biot field (``snap <- live``) before the candidate field solve."""
        for live, snap in self._f64_pairs():
            self._copy(snap, live)
        self._copy(self._snap_mask_d, self._grid.mask)

    def rollback(self, accepted: wp.array) -> None:
        """Reject-gated restore of every field array; an accepted step keeps the converged candidate field."""
        self.launched_kernels = []
        for live, snap in self._f64_pairs():
            self._launch(_conditional_restore_f64_3d_kernel, dim=self._shape,
                         inputs=[live, snap, accepted], device=self._device)
            self.launched_kernels.append(_conditional_restore_f64_3d_kernel)
        self._launch(_conditional_restore_i32_3d_kernel, dim=self._shape,
                     inputs=[self._grid.mask, self._snap_mask_d, accepted], device=self._device)
        self.launched_kernels.append(_conditional_restore_i32_3d_kernel)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """No-op: the pressure field is reversible state (the swept-content + epoch commit elsewhere)."""
        self.commit_calls.append((accepted, dt_phys, rng_seed))


def build_fluid_field_owner(
    *,
    grid: object,
    solver: object,
    device: str,
    launch: LaunchFn = _default_launch,
    copy: CopyFn = _default_copy,
    ledger: object | None = None,
) -> FluidVolumeStateOwner:
    """Wire a native ``cytosol`` :class:`FluidVolumeStateOwner` over a live Biot ``FieldGrid`` (CUDA lane).

    ``solver`` is the injected :class:`~aleph.engine.fluid_core.FluidCandidateSolver` (production:
    :class:`~aleph.engine.fluid_core.BiotSubstrateFluidSolver`, which launches the real
    ``biot_pmass_update_kernel``).  This allocates the field snapshot twins and wires the concrete
    :class:`FluidFieldTransaction`; the owner exposes ``grid.p`` / ``grid.mask`` as its authoritative arrays.
    """
    shape = tuple(int(n) for n in grid.shape)

    def zf(dtype: object) -> wp.array:
        return wp.zeros(shape, dtype=dtype, device=device)

    transaction = FluidFieldTransaction(
        grid=grid, snap_p_d=zf(wp.float64), snap_mask_d=zf(wp.int32), snap_div_vs_d=zf(wp.float64),
        snap_s_membrane_d=zf(wp.float64), snap_s_total_d=zf(wp.float64), shape=shape,
        launch=launch, copy=copy, device=device,
    )
    return FluidVolumeStateOwner(
        name=CYTOSOL_COMPONENT, pressure_d=grid.p, mask_d=grid.mask,
        solver=solver, transaction=transaction,
        ledger=ledger if ledger is not None else RecordingLedgerStub(),
    )


def build_cortex_cytosol_transfer(
    *,
    pressure_coupling: object,
    node_volume_d: wp.array,
) -> ImmersedPorousTransfer:
    """Assemble the ``surface_porous_transfer`` connector over a real :class:`PressureCoupling` (CUDA lane).

    The cortex edge is one of SEVEN in the ``immersed_transfer`` family and shares its runtime with the
    other six (:mod:`aleph.engine.immersed_transfer`); nothing in ``-alpha*V*grad(p)`` is cortical.  This
    stays as the cortex-named entry point the interior-column driver already calls.
    """
    return build_immersed_porous_transfer(
        name=SURFACE_POROUS_TRANSFER, pressure_coupling=pressure_coupling, node_volume_d=node_volume_d,
    )


# ── nucleus full-surface Core Body owner + its snapshot/restore + accepted-gated rupture commit. ───────
class CoreSurfaceTransaction:
    """Nucleus surface snapshot / reject-restore + the accepted-gated envelope-rupture commit.

    The native nucleus carries its full deformable surface as authoritative state (the ROM is unbuilt — the
    surface positions ARE the generalized coordinates).  This snapshots the reversible surface position (vec3d)
    + the per-face ``ruptured`` flag (int32); a rejected step reject-restores both in-device.  On an ACCEPTED
    step :meth:`commit_irreversible` launches the real ``ac/nucleus.conditional_rupture_update_kernel`` — the
    ONE irreversible interior-column kinetic (a face tears when its areal strain exceeds ``eps_rupt``), gated on
    the device accepted predicate, exactly as ``NucleusCompartment.rupture_step`` drives it in the incumbent.

    ``eps_rupt`` is the open I0-B2 GAP-PI rupture strain (``params_i0b2.yaml``: "do NOT keep the ff proxy
    0.50"); it is REQUIRED here (no default) and the native script supplies it with a loud PROVISIONAL banner.
    """

    def __init__(
        self,
        *,
        position_d: wp.array,
        position_snap_d: wp.array,
        faces_d: wp.array,
        a0_d: wp.array,
        ruptured_d: wp.array,
        ruptured_snap_d: wp.array,
        eps_rupture: float,
        launch: LaunchFn = _default_launch,
        copy: CopyFn = _default_copy,
        device: object | None = None,
    ) -> None:
        self._position_d = position_d
        self._position_snap_d = position_snap_d
        self._faces_d = faces_d
        self._a0_d = a0_d
        self._ruptured_d = ruptured_d
        self._ruptured_snap_d = ruptured_snap_d
        self._eps_rupture = _require_param("eps_rupture", eps_rupture, gap=True)
        self._n_faces = int(faces_d.shape[0])
        self._launch = launch
        self._copy = copy
        self._device = device
        self.launched_kernels: list[object] = []
        self.commit_calls: list[tuple[Any, ...]] = []

    def snapshot_candidate(self) -> None:
        """D2D-snapshot the reversible surface position + the per-face ``ruptured`` flag."""
        self._copy(self._position_snap_d, self._position_d)
        self._copy(self._ruptured_snap_d, self._ruptured_d)

    def rollback(self, accepted: wp.array) -> None:
        """Reject-gated restore of the surface position (vec3d) + ``ruptured`` flag (int32)."""
        self.launched_kernels = []
        self._launch(_restore_vec3_if_rejected_kernel, dim=int(self._position_d.shape[0]),
                     inputs=[accepted, self._position_d, self._position_snap_d], device=self._device)
        self.launched_kernels.append(_restore_vec3_if_rejected_kernel)
        self._launch(_restore_int32_if_rejected_kernel, dim=self._n_faces,
                     inputs=[accepted, self._ruptured_d, self._ruptured_snap_d], device=self._device)
        self.launched_kernels.append(_restore_int32_if_rejected_kernel)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Accepted-gated per-face envelope rupture on the converged surface (device-gated; the one kinetic)."""
        self.commit_calls.append((accepted, dt_phys, rng_seed))
        self.launched_kernels = []
        self._launch(
            conditional_rupture_update_kernel, dim=self._n_faces,
            inputs=[self._position_d, self._faces_d, self._a0_d, wp.float64(self._eps_rupture),
                    accepted, self._ruptured_d], device=self._device,
        )
        self.launched_kernels.append(conditional_rupture_update_kernel)


@dataclass(slots=True)
class CoreSurfaceStateOwner:
    """The ``nucleus`` Core Body component and the CUDA surface arrays it exclusively owns (full-surface).

    ``mechanics`` may be the real :class:`~aleph.engine.fluid_core.NativeSurfaceCoreForces` (lamina +
    nucleoplasm volume + LINC) for a standalone column, or :class:`NoOpSurfaceMechanics` when the native driver
    owns the nucleus force (no double count).  ``transaction`` is the :class:`CoreSurfaceTransaction`
    snapshot/restore + rupture commit; ``ledger`` records the ledger handle.
    """

    name: str
    position_d: wp.array
    force_d: wp.array
    mechanics: Any
    transaction: CoreSurfaceTransaction
    ledger: Any

    def __post_init__(self) -> None:
        if self.name != CORE_COMPONENT:
            raise ValueError(f"core surface owner name must be {CORE_COMPONENT!r}")

    def accumulate(self) -> None:
        """Launch this component's mechanics against its own surface arrays (no-op under the native driver)."""
        self.mechanics.accumulate(self.position_d, self.force_d)

    def snapshot_candidate(self) -> None:
        """Forward the global transaction snapshot hook to the surface transaction."""
        self.transaction.snapshot_candidate()

    def rollback(self, accepted: wp.array) -> None:
        """Forward the scheduler-owned predicate without inspecting it on the host."""
        self.transaction.rollback(accepted)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Forward the accepted-gated envelope-rupture commit without owning the physical clock."""
        self.transaction.commit_irreversible(accepted, dt_phys, rng_seed)

    def accumulate_ledger(self, ledger: object) -> None:
        """Forward the reduced residual / volume / strain / no-flux terms without deciding acceptance."""
        self.ledger.accumulate_ledger(ledger)


def build_core_surface_owner(
    *,
    position_d: wp.array,
    force_d: wp.array,
    faces_d: wp.array,
    a0_d: wp.array,
    ruptured_d: wp.array,
    eps_rupture: float,
    device: str,
    mechanics: Any | None = None,
    launch: LaunchFn = _default_launch,
    copy: CopyFn = _default_copy,
    ledger: Any | None = None,
) -> CoreSurfaceStateOwner:
    """Wire a native ``nucleus`` :class:`CoreSurfaceStateOwner` from a built deformable-mesh envelope (CUDA lane).

    ``position_d``/``force_d`` are the live (global-index) surface arrays; ``faces_d``/``a0_d``/``ruptured_d``
    are the per-face topology + reference areas + tear flags from the built ``NucleusCompartment``.  Allocates
    the position + ruptured snapshot twins and wires the :class:`CoreSurfaceTransaction`.  ``mechanics``
    defaults to :class:`NoOpSurfaceMechanics` (the native driver owns the nucleus force — no double count).
    """
    n_nodes = int(position_d.shape[0])
    n_faces = int(faces_d.shape[0])
    transaction = CoreSurfaceTransaction(
        position_d=position_d, position_snap_d=wp.zeros(n_nodes, dtype=wp.vec3d, device=device),
        faces_d=faces_d, a0_d=a0_d, ruptured_d=ruptured_d,
        ruptured_snap_d=wp.zeros(n_faces, dtype=wp.int32, device=device),
        eps_rupture=eps_rupture, launch=launch, copy=copy, device=device,
    )
    return CoreSurfaceStateOwner(
        name=CORE_COMPONENT, position_d=position_d, force_d=force_d,
        mechanics=mechanics if mechanics is not None else NoOpSurfaceMechanics(),
        transaction=transaction, ledger=ledger if ledger is not None else RecordingLedgerStub(),
    )


# ── The wired slice: one accepted physical step over the membrane → cytosol → nucleus interior column. ──
@dataclass(slots=True)
class InteriorColumnSlice:
    """One radial column of the cell graph wired for the accepted-step transaction.

    Holds the four state-owners (membrane + cortex surface owners, the cytosol Biot field owner, the nucleus
    full-surface owner) and the three fluid coupling connectors (``surface_porous_transfer`` cortex↔cytosol,
    ``membrane_cytosol_boundary`` semipermeable, ``nucleus_cytosol_boundary`` impermeable), the device event
    clock, and the top-level :class:`CellTransaction`.  Constructed from already-assembled parts, so it is
    CPU-importable and the structural gates build it with recording doubles; :func:`build_interior_column_slice`
    assembles the real Warp-resident parts on the CUDA lane.

    ``cytosol_endpoint`` is the public :class:`CytosolFieldEndpoint` the porous transfer reads; the membrane /
    nucleus surface velocity arrays feed the moving boundaries' ``update_geometry``.
    """

    transaction: CellTransaction
    membrane: Any
    cortex: Any
    cytosol: Any
    nucleus: Any
    cortex_transfer: ImmersedPorousTransfer
    membrane_boundary: MovingSemipermeableFluidBoundaryFacade
    nucleus_boundary: MovingImpermeableFluidBoundaryFacade
    cytosol_endpoint: Any
    clock: Any
    membrane_surface_velocity_d: Any = None
    nucleus_surface_velocity_d: Any = None
    ledger: Any | None = None
    tol_sq_d: wp.array | None = None

    def accumulate(self) -> None:
        """One coupled candidate-force pass the caller's inner solve invokes each iteration (after zeroing).

        Refreshes both moving-boundary geometries, then scatters every interior coupling load as a bidirectional
        adjoint pair (Newton's 3rd law): the two ``*_cytosol_boundary`` pressure tractions onto the membrane and
        nuclear-envelope surfaces, the ``surface_porous_transfer`` pressure traction onto the cortex, plus the
        compartments' own internal mechanics.  No fixed prestress is injected; each coupling force is the
        reaction of the coupled candidate the world solver relaxes.  (In the native no-double-count driver the
        compartment mechanics are no-ops — the driver's ``_accumulate_all`` owns them.)
        """
        self.membrane_boundary.update_geometry(self.membrane.position_d, self.membrane_surface_velocity_d)
        self.nucleus_boundary.update_geometry(self.nucleus.position_d, self.nucleus_surface_velocity_d)
        self.membrane.accumulate()
        self.cortex.accumulate()
        self.nucleus.accumulate()
        self.membrane_boundary.accumulate(self.membrane.position_d, self.membrane.force_d)
        self.nucleus_boundary.accumulate(self.nucleus.position_d, self.nucleus.force_d)
        self.cortex_transfer.accumulate_transfer(
            SurfacePorousEndpointView(CORTEX_COMPONENT, self.cortex.position_d, self.cortex.force_d),
            self.cytosol_endpoint,
        )

    def accumulate_ledger(self, ledger: object) -> None:
        """Collect the coupled MASS / NO-FLUX / ADJOINT-WORK channels from every interior-column participant."""
        for participant in (self.membrane, self.cortex, self.cytosol, self.nucleus,
                            self.cortex_transfer, self.membrane_boundary, self.nucleus_boundary):
            participant.accumulate_ledger(ledger)

    def step(
        self,
        solve: Callable[[], None],
        *,
        dt_phys: float,
        accepted_d: wp.array | None = None,
        ledger: Any | None = None,
        tol_sq_d: wp.array | None = None,
    ) -> None:
        """Run exactly one accepted physical step of the dynamic interior-column runtime.

        Args:
            solve: the caller-owned inner mechanical solve (zero forces → :meth:`accumulate` + Biot field
                relaxation → integrate to convergence).  Runs after events are proposed and owns all physics.
            dt_phys: the outer physical timestep [s]; finite and positive.
            accepted_d: the device acceptance predicate.  If ``None`` the balance gate decides (a ``ledger``
                with a ``tol_sq_d`` must then be available).
            ledger / tol_sq_d: optional balance-gate inputs; default to the slice's stored ones.

        The nucleus rupture commit reads the converged geometry directly, so no per-face load pre-pass is
        needed.  The predicate stays device-resident throughout.
        """
        ledger = ledger if ledger is not None else self.ledger
        tol_sq_d = tol_sq_d if tol_sq_d is not None else self.tol_sq_d
        self.transaction.step(
            dt_phys=dt_phys, solve=solve, accepted_d=accepted_d, ledger=ledger, tol_sq_d=tol_sq_d,
        )

    def ruptured_face_count(self) -> int:
        """Host read of the ruptured-face count for OUT-OF-LOOP diagnostics only (device→host readback).

        A telemetry helper the caller may call BETWEEN steps to watch envelope tearing under load; it must
        never be called inside the inner solve (I0-A: no authoritative GPU→CPU roundtrip in the hot loop).
        """
        return int(self.nucleus.transaction._ruptured_d.numpy().sum())


def build_interior_column_slice(
    *,
    membrane_owner: Any,
    cortex_owner: Any,
    cytosol_owner: Any,
    nucleus_owner: Any,
    cortex_transfer: ImmersedPorousTransfer,
    membrane_boundary: MovingSemipermeableFluidBoundaryFacade,
    nucleus_boundary: MovingImpermeableFluidBoundaryFacade,
    cytosol_endpoint: Any,
    base_seed: int,
    device: str,
    membrane_surface_velocity_d: Any = None,
    nucleus_surface_velocity_d: Any = None,
    clock: Any | None = None,
    ledger: Any | None = None,
    tol_sq_d: wp.array | None = None,
    population_ledgers: tuple[Any, ...] = (),
    architecture: CellArchitecture | None = None,
    additional_component_bindings: tuple[tuple[str, Any], ...] = (),
    additional_connector_bindings: tuple[tuple[str, Any], ...] = (),
) -> InteriorColumnSlice:
    """Assemble the dynamic interior-column slice: four owners + three fluid connectors under one clock.

    Binds the membrane/cortex surface owners, the cytosol Biot field owner, and the nucleus full-surface owner
    as the ``membrane``/``cortex``/``cytosol``/``nucleus`` components, and the three fluid connectors as
    ``surface_porous_transfer`` / ``membrane_cytosol_boundary`` / ``nucleus_cytosol_boundary`` on the reference
    architecture (which already carries all four components + all three edges).  Then wires a
    :class:`CellWorldTransaction` (``require_complete=False`` — a vertical slice, not the whole cell) and a
    :class:`CellTransaction`.  The three fluid connectors are the ONLY interior mechanical couplings.

    Args:
        membrane_owner / cortex_owner: native ``membrane``/``cortex`` :class:`SurfaceComponentStateOwner`.
        cytosol_owner: the ``cytosol`` :class:`FluidVolumeStateOwner` (see :func:`build_fluid_field_owner`).
        nucleus_owner: the ``nucleus`` :class:`CoreSurfaceStateOwner` (see :func:`build_core_surface_owner`).
        cortex_transfer / membrane_boundary / nucleus_boundary: the three graph fluid connectors.
        cytosol_endpoint: the public :class:`CytosolFieldEndpoint` the porous transfer reads.
        base_seed: the run's fixed host RNG seed (nonnegative int).
        device: the resolved CUDA device string (never a hard-coded id).
        clock / ledger / tol_sq_d: optional; a :class:`WarpEventClock` is built from ``base_seed`` if omitted.
        population_ledgers: optional disjoint-population ledgers asserted at build + each accepted step.
        architecture: optional composition; defaults to :func:`reference_cell_architecture`.
        additional_component_bindings / additional_connector_bindings: experiment-scoped bindings added
            before the world transaction captures its participant set. This is how an apparatus joins the
            same accepted-step predicate without changing the reference-cell census.

    Returns:
        An :class:`InteriorColumnSlice` whose :meth:`~InteriorColumnSlice.step` runs one accepted physical step.
    """
    if architecture is None:
        architecture = reference_cell_architecture()

    actor = CellActor(architecture)
    actor.bind_component(MEMBRANE_COMPONENT, membrane_owner)
    actor.bind_component(CORTEX_COMPONENT, cortex_owner)
    actor.bind_component(CYTOSOL_COMPONENT, cytosol_owner)
    actor.bind_component(NUCLEUS_COMPONENT, nucleus_owner)
    actor.bind_connector(SURFACE_POROUS_TRANSFER, cortex_transfer)
    actor.bind_connector(MEMBRANE_FLUID_BOUNDARY, membrane_boundary)
    actor.bind_connector(NUCLEUS_FLUID_BOUNDARY, nucleus_boundary)
    for name, runtime in additional_component_bindings:
        actor.bind_component(name, runtime)
    for name, runtime in additional_connector_bindings:
        actor.bind_connector(name, runtime)
    world = CellWorldTransaction(actor, require_complete=False)

    if clock is None:
        from aleph.engine.events import make_event_clock  # CUDA-lane import (allocates device scalars).
        clock = make_event_clock(base_seed=base_seed, device=device)
    transaction = CellTransaction(world, clock, population_ledgers=tuple(population_ledgers))
    return InteriorColumnSlice(
        transaction=transaction, membrane=membrane_owner, cortex=cortex_owner, cytosol=cytosol_owner,
        nucleus=nucleus_owner, cortex_transfer=cortex_transfer, membrane_boundary=membrane_boundary,
        nucleus_boundary=nucleus_boundary, cytosol_endpoint=cytosol_endpoint, clock=clock,
        membrane_surface_velocity_d=membrane_surface_velocity_d,
        nucleus_surface_velocity_d=nucleus_surface_velocity_d, ledger=ledger, tol_sq_d=tol_sq_d,
    )


# ── the native interior-column recipe, in ONE place ───────────────────────────────────────────────────

@dataclass(frozen=True, slots=True)
class InteriorColumnParts:
    """The four owners + three connectors of the interior column, built over a native cell's own arrays.

    Attributes:
        membrane / cortex / nucleus: :class:`SurfaceComponentStateOwner` / :class:`CoreSurfaceStateOwner`
            addressing DISJOINT blocks of the incumbent's global node arrays.
        cytosol: the Biot field owner (holds no node array, so it is not party to the disjointness check).
        cortex_transfer / membrane_boundary / nucleus_boundary: the three fluid connectors.
        cytosol_endpoint: the public endpoint the porous transfer reads.
        blocks: ``{component: (node_offset, n_nodes)}`` — the addressing this build certified.
        arrays_private: **always False.** These are SLICE VIEWS of one global array, not private
            allocations. A private allocation without splitting the incumbent's inner solve would make
            the component a passenger — distinct arrays the solve never updates, which is exactly the T2
            negative control. What the views buy is exclusive ADDRESSING. `STATE.md` (c) 4 stands.
    """

    membrane: Any
    cortex: Any
    cytosol: Any
    nucleus: Any
    cortex_transfer: Any
    membrane_boundary: Any
    nucleus_boundary: Any
    cytosol_endpoint: Any
    blocks: dict[str, tuple[int, int]]
    arrays_private: bool = False


def build_native_interior_column_parts(
    cell: Any,
    *,
    device: str,
    osmotic_pi0_pa: float,
    participant: Any | None = None,
    ledger_stub: Any | None = None,
) -> InteriorColumnParts:
    """Build membrane / cortex / cytosol / nucleus owners + their three fluid connectors from a built cell.

    This was 60 lines inside ``ac_gate_b_interior_column_native.py``. It moved here so the composed native
    world and that driver call the SAME recipe: a second copy of a build this delicate — local face tables,
    per-component node-volume weights, disjoint block views — would drift, and the drift would look like a
    physics difference between two lanes.

    Args:
        cell: a built incumbent cell with ``grid``, ``substrate``, ``pressure``, ``membrane_bc``, ``domain``,
            ``nucleus`` and ``membrane`` present (``with_pressure``/``with_membrane``/``with_nucleus``).
        device: resolved CUDA device string.
        osmotic_pi0_pa: van 't Hoff Π₀ [Pa]. **Required, no default** — it is a PI-GAP datum and a default
            here would let a convenient value reach production silently.
        participant: transaction participant the moving-boundary delegates forward to; a no-op when the
            caller's own inner solve owns the moving-domain remap and pressure caches.
        ledger_stub: energy-ledger sink for the two boundary delegates.

    Raises:
        RuntimeError: if a required compartment is absent — an interior column without its fluid grid is
            not a weaker column, it is a different object.
    """
    import numpy as _np

    from aleph.components.fluid.boundary import NucleusNoFluxBC
    from aleph.engine.cortex_state import assert_component_state_disjoint
    from aleph.engine.cytosol_connected import (
        MembranePressureFluxAdjointBoundary,
        MovingSemipermeableFluidBoundaryFacade,
    )
    from aleph.engine.erm_cortex_slice import build_native_surface_owner
    from aleph.engine.fluid_core import (
        BiotSubstrateFluidSolver,
        MovingImpermeableFluidBoundaryFacade,
        NucleusPressureAdjointBoundary,
        lumped_surface_control_volumes,
    )

    for name in ("grid", "substrate", "pressure", "membrane_bc", "domain", "nucleus", "membrane"):
        if getattr(cell, name, None) is None:
            raise RuntimeError(
                f"the interior column needs the {name!r} compartment; build the cell with "
                f"with_pressure/with_membrane/with_nucleus=True"
            )

    blocks = {
        "cortex": (0, int(cell.n_actin)),
        "nucleus": (int(cell.nucleus.node_off), int(cell.nucleus.n_verts)),
        "membrane": (int(cell.membrane.node_off), int(cell.membrane.n_verts)),
    }
    assert_disjoint_blocks(blocks, n_total=int(cell.n_total))

    views = {
        name: component_block_views(cell.pos_d, cell.f_d, node_off=off, n_nodes=n)
        for name, (off, n) in blocks.items()
    }

    # The nucleus rupture kernel indexes `position_d[faces[t, k]]`, so a LOCAL position view needs a LOCAL
    # face table. `cell.nucleus.faces_d` stays untouched — the incumbent still reads it with global indices.
    pos_np = cell.pos_d.numpy()
    dx = float(cell.grid.dx)

    def _local(surface):
        faces = surface.faces_d.numpy().astype(_np.int64) - int(surface.node_off)
        verts = pos_np[int(surface.node_off):int(surface.node_off) + int(surface.n_verts)]
        return verts, faces

    m_verts, m_faces = _local(cell.membrane)
    n_verts, n_faces = _local(cell.nucleus)

    membrane_owner = build_native_surface_owner(
        name="membrane", position_d=views["membrane"][0], force_d=views["membrane"][1], device=device)
    cortex_owner = build_native_surface_owner(
        name="cortex", position_d=views["cortex"][0], force_d=views["cortex"][1], device=device)
    cytosol_owner = build_fluid_field_owner(
        grid=cell.grid, solver=BiotSubstrateFluidSolver(cell.substrate), device=device)
    nucleus_owner = build_core_surface_owner(
        position_d=views["nucleus"][0], force_d=views["nucleus"][1],
        faces_d=wp.array(_np.ascontiguousarray(n_faces, _np.int32), dtype=wp.int32, device=device),
        a0_d=cell.nucleus.a0_d, ruptured_d=cell.nucleus.ruptured_d,
        eps_rupture=float(cell.nucleus.eps_rupt), device=device)

    # Definition of done, checked on the built objects rather than declared. Overlapping views would pass a
    # pointer-identity test and are refused here by byte span.
    assert_component_state_disjoint([membrane_owner, cortex_owner, nucleus_owner])

    participant = participant if participant is not None else _NoOpTransactionParticipant()
    ledger_stub = ledger_stub if ledger_stub is not None else RecordingLedgerStub()

    membrane_node_vol = wp.array(
        lumped_surface_control_volumes(m_verts, m_faces, dx), dtype=wp.float64, device=device)
    nucleus_node_vol = wp.array(
        lumped_surface_control_volumes(n_verts, n_faces, dx), dtype=wp.float64, device=device)
    cortex_node_vol = wp.full(int(cell.n_actin), dx ** 3, dtype=wp.float64, device=device)

    cortex_transfer = build_cortex_cytosol_transfer(
        pressure_coupling=cell.pressure, node_volume_d=cortex_node_vol)
    membrane_boundary = MovingSemipermeableFluidBoundaryFacade(
        delegate=MembranePressureFluxAdjointBoundary(
            cell.membrane_bc, cell.pressure, membrane_node_vol, participant, ledger_stub,
            osmotic_difference=float(osmotic_pi0_pa),
        ))
    nucleus_boundary = MovingImpermeableFluidBoundaryFacade(
        delegate=NucleusPressureAdjointBoundary(
            cell.domain, cell.pressure,
            NucleusNoFluxBC(cell.grid, mobility=float(cell.substrate.mobility)),
            nucleus_node_vol, participant, ledger_stub))

    residual_d = wp.zeros(cell.grid.shape, dtype=wp.float64, device=device)
    return InteriorColumnParts(
        membrane=membrane_owner, cortex=cortex_owner, cytosol=cytosol_owner, nucleus=nucleus_owner,
        cortex_transfer=cortex_transfer, membrane_boundary=membrane_boundary,
        nucleus_boundary=nucleus_boundary,
        cytosol_endpoint=cytosol_owner.make_immersed_transfer_endpoint(residual_d),
        blocks=blocks,
    )


class _NoOpTransactionParticipant:
    """Forwarding target for the moving-boundary delegates when the caller's solve owns those caches."""

    def snapshot_candidate(self) -> None: ...
    def rollback(self, accepted: wp.array) -> None: ...
    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None: ...
