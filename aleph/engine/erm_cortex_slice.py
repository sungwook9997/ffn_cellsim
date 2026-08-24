r"""GATE-B second dynamic-runtime vertical slice: emergent membrane--cortex ERM tether tension.

The first dynamic slice (:mod:`aleph.engine.cortex_motor_slice`) grew active cortical tension from NMII
head-binding EVENTS.  This slice extends the actomyosin cortex OUTWARD toward the plasma membrane: it drives
:class:`~aleph.engine.transaction.CellTransaction` over one vertical column of the composed cell graph::

    membrane  ──(membrane_erm_cortex : ERM connector)──  cortex

Per accepted physical step the transaction:

1. **snapshots** the two surface state-owners and the ERM connector's ``bound``/``rest``/``epoch`` SoA;
2. **proposes events** — a near-no-op (:meth:`ErmCortexConnector.propose_events`): ERM pairs are FIXED by
   radial pairing, and a shed tether rebinds geometrically-within-capture INSIDE the KMC (there is no
   point-to-segment attach query — this is why ERM is the smaller, simpler next participant after NMII);
3. runs the caller-owned **inner mechanical solve** — during which :meth:`ERMCortexSlice.accumulate` scatters
   the two-array Bell tether load ``-f`` into the membrane force array and ``+f`` into the cortex force array
   (Newton's 3rd law across two never-merged owner arrays), and the membrane relaxes onto the cortex;
4. **commits** the accepted-gated ERM Bell KMC on the CONVERGED geometry: a bound tether SLIP-detaches under
   its live tensile load (Bell 1978; ERM is a pure slip bond) and a shed pair rebinds within the capture
   radius — or rolls the SoA back bit-exactly on a rejected step;
5. **advances the device event clock** iff the step was accepted.

No fixed rupture threshold is imposed: the membrane--cortex tether tension EMERGES from the bound-tether
population the KMC maintains, and under load the Bell shedding IS the mechanistic bleb-onset (it replaces the
hard ``f_rupt``).  Every CUDA op flows through an injected launcher/copy so the module is CPU-importable; the
CUDA-lane builders assemble the real Warp-resident runtimes on the gbook A5000.

Sanity Gate:
    * ownership: membrane and cortex are the two state-owners; the ``membrane_erm_cortex`` connector is the
      ONLY mechanical coupling (co-location is never a connection).  The ERM SoA is owned by the connector.
    * boundary/sign: the two-array tether is momentum-conserving; the KMC and clock advance are device-gated
      on the accepted predicate — a rejected step advances no ERM state and no clock (no host read of it).
    * REQUIRED-PARAM: ``k_on``/``k_off0``/``bell_force_pn``/``capture_radius_um`` are PI/KB GAPs (no MCF7 ERM
      Bell datum in the Contract-Graph, 2026-07-21 audit) — :class:`ERMSliceParams` refuses to construct
      without them; ``k_erm`` is the sourced Braunger 4.6e3 pN/µm.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Any, Callable

import warp as wp

from aleph.components.incumbent.erm_tether import ERMBellKinetics
from aleph.engine.actor import CellActor
from aleph.engine.contracts import CellArchitecture, reference_cell_architecture
from aleph.engine.erm_cortex_connector import ErmCortexConnector
from aleph.engine.surface_body import (
    CORTEX_COMPONENT,
    ERM_CONNECTOR,
    MEMBRANE_COMPONENT,
    SurfaceComponentStateOwner,
)
from aleph.engine.transaction import CellTransaction
from aleph.engine.world import CellWorldTransaction

__all__ = [
    "ERMCortexSlice",
    "ERMSliceParams",
    "build_erm_cortex_connector",
    "build_erm_cortex_slice",
    "build_native_surface_owner",
]

LaunchFn = Callable[..., None]
CopyFn = Callable[..., None]


def _default_launch(kernel: object, *, dim: int, inputs: list, outputs: list | None = None,
                    device: object | None = None) -> None:
    """Production launcher: forward to ``wp.launch`` on the owning CUDA device."""
    if outputs is None:
        wp.launch(kernel, dim=dim, inputs=inputs, device=device)
    else:
        wp.launch(kernel, dim=dim, inputs=inputs, outputs=outputs, device=device)


def _default_copy(dst: object, src: object) -> None:
    """Production D2D snapshot/restore: forward to ``wp.copy``."""
    wp.copy(dst, src)


@wp.kernel
def _restore_vec3_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.vec3d),
    snap: wp.array(dtype=wp.vec3d),
) -> None:
    """Restore a vec3 array ``live <- snap`` iff the outer step was rejected (device-only, no host branch)."""
    if accepted[0] == 0:
        t = wp.tid()
        live[t] = snap[t]


# ── REQUIRED-PARAM discipline: the 4 ERM Bell-kinetic GAPs are supplied or the config refuses to build. ──
#: The Contract-Graph has NO MCF7 ERM Bell on/off datum (2026-07-21 audit); these four are supplied by the
#: caller (with an explicit provenance banner), never defaulted to a convenient number here.  ``k_erm`` is
#: NOT in this set — it is the sourced Braunger 4.6e3 pN/µm single-bond stiffness (PI-ratified 2026-07-21).
_ERM_GAP_PARAMS = frozenset({"k_on", "k_off0", "bell_force_pn", "capture_radius_um"})


def _require_param(name: str, value: float | None) -> float:
    """Return a supplied positive-finite param, or raise a REQUIRED-PARAM error naming the GAP."""
    gap = " (KB/PI GAP — source it or surface to PI, never default)" if name in _ERM_GAP_PARAMS else ""
    if value is None or not isfinite(value) or value <= 0.0:
        raise ValueError(f"REQUIRED-PARAM {name!r} must be a supplied positive-finite value{gap}; got {value!r}")
    return float(value)


@dataclass(frozen=True, slots=True)
class ERMSliceParams:
    """Sourced-or-GAP stiffness + Bell on/off constants for the ERM membrane--cortex clutch (engine units).

    ``k_erm`` is the sourced single-ezrin--F-actin linker stiffness (Braunger 2014 JBC single-bond estimate,
    4.6 pN/nm = 4.6e3 pN/µm; PI-ratified 2026-07-21, KB-3.B1.6 registration pending).  The four Bell kinetics
    (``k_on``/``k_off0``/``bell_force_pn``/``capture_radius_um``) are a PI/KB GAP — no MCF7 ERM on/off datum
    exists in the Contract-Graph — so they are REQUIRED (no defaults); the caller supplies provisional proxies
    with a loud banner (see :mod:`aleph.scripts.ac_gate_b_erm_cortex_native`).  ``bell_force_pn`` is the Bell
    characteristic force ``F0 = kBT/x_beta`` — use :func:`aleph.components.motor.bell_kinetics_analytic.bell_f0_from_x_beta`
    when an ERM bond length ``x_beta`` is available.  ERM is a pure SLIP bond (the off-rate rises with load).
    """

    k_erm: float
    k_on: float
    k_off0: float
    bell_force_pn: float
    capture_radius_um: float
    source: str
    rebind_rest_policy: str = "formation_length"

    def __post_init__(self) -> None:
        for name in ("k_erm", "k_on", "k_off0", "bell_force_pn", "capture_radius_um"):
            _require_param(name, getattr(self, name))
        if not self.source.strip():
            raise ValueError("ERM slice params require an explicit literature/provenance source")
        if self.rebind_rest_policy != "formation_length":
            raise ValueError("only the explicit ERM rebind_rest_policy='formation_length' is implemented")

    def to_bell_kinetics(self) -> ERMBellKinetics:
        """Build the source-gated :class:`ERMBellKinetics` contract this slice's KMC consumes."""
        return ERMBellKinetics(
            k_on_s=self.k_on, k_off0_s=self.k_off0, bell_force_pn=self.bell_force_pn,
            capture_radius_um=self.capture_radius_um, source=self.source,
            rebind_rest_policy=self.rebind_rest_policy,
        )


# ── Native-lane surface state-owner delegates (the ERM demo's membrane/cortex forces are DRIVER-owned). ──
class _NoOpSurfaceMechanics:
    """No-op mechanics for a native surface owner whose forces the driver's own inner solve assembles.

    In the native ERM driver the membrane Helfrich/area + cortex network + turgor forces are assembled by the
    driver's ``_accumulate_all`` (which the membrane adapter augments with the slice's ERM tether force);
    adding them again through the owner would DOUBLE-COUNT, so the owner mechanics is deliberately a no-op —
    the same reason the NMII driver routes its motor through the driver's ``_accumulate_all``, not the slice.
    """

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """No-op: the driver's inner solve owns the membrane/cortex mechanics (no double count)."""


class _SurfaceStateLedger:
    """Minimal :class:`LedgerContributor` for a native surface owner (records the ledger handle)."""

    def __init__(self) -> None:
        self.ledger_calls: list[object] = []

    def accumulate_ledger(self, ledger: object) -> None:
        """Record the ledger handle without deciding acceptance."""
        self.ledger_calls.append(ledger)


class _SurfaceStateTransaction:
    """Reject-gated position snapshot/restore for one native surface state-owner (no kinetic epoch).

    The membrane/cortex positions are relaxed by the driver's inner solve; this delegate snapshots them so a
    rejected outer step restores the pre-step geometry in-device.  It advances no epoch on commit — the ERM
    kinetic epoch is the connector's.
    """

    def __init__(
        self,
        *,
        component_name: str,
        position_d: wp.array,
        position_snap_d: wp.array,
        launch: LaunchFn = _default_launch,
        copy: CopyFn = _default_copy,
        device: object | None = None,
    ) -> None:
        self.component_name = component_name
        self._position_d = position_d
        self._position_snap_d = position_snap_d
        self._launch = launch
        self._copy = copy
        self._device = device

    def snapshot_candidate(self) -> None:
        """D2D-snapshot the surface position (``snap <- live``) before the candidate solve."""
        self._copy(self._position_snap_d, self._position_d)

    def rollback(self, accepted: wp.array) -> None:
        """Reject-gated restore of the surface position; an accepted step keeps the converged geometry."""
        self._launch(_restore_vec3_if_rejected_kernel, dim=int(self._position_d.shape[0]),
                     inputs=[accepted, self._position_d, self._position_snap_d], device=self._device)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """No-op: a surface owner carries no kinetic epoch in this slice (the ERM epoch is the connector's)."""


def build_native_surface_owner(
    *,
    name: str,
    position_d: wp.array,
    force_d: wp.array,
    device: str,
    launch: LaunchFn = _default_launch,
    copy: CopyFn = _default_copy,
) -> SurfaceComponentStateOwner:
    """Wire a native ``membrane``/``cortex`` :class:`SurfaceComponentStateOwner` (CUDA lane; allocates a snap).

    The owner's mechanics is a no-op (the driver's inner solve owns the surface forces — see
    :class:`_NoOpSurfaceMechanics`); its transaction snapshots/reject-restores the position, and its ledger is
    a recording stub.  ``position_d``/``force_d`` are the live arrays the owner exposes to the slice's tether
    scatter (in the native driver both surfaces alias one global array; the split is logical).
    """
    n = int(position_d.shape[0])
    transaction = _SurfaceStateTransaction(
        component_name=name, position_d=position_d,
        position_snap_d=wp.zeros(n, dtype=wp.vec3d, device=device),
        launch=launch, copy=copy, device=device,
    )
    return SurfaceComponentStateOwner(
        name=name, position_d=position_d, force_d=force_d,
        mechanics=_NoOpSurfaceMechanics(), transaction=transaction, ledger=_SurfaceStateLedger(),
    )


def build_erm_cortex_connector(
    *,
    membrane_idx_d: wp.array,
    cortex_idx_d: wp.array,
    bound_d: wp.array,
    rest_d: wp.array,
    membrane_pos_d: wp.array,
    cortex_pos_d: wp.array,
    params: ERMSliceParams,
    device: str,
    rng_epoch_d: wp.array | None = None,
    detach_events_d: wp.array | None = None,
    attach_events_d: wp.array | None = None,
) -> ErmCortexConnector:
    """Assemble an :class:`ErmCortexConnector` from a live ERM SoA (CUDA lane; allocates the snapshot twins).

    The live pair indices + ``bound``/``rest`` SoA come from the built cell's radial-paired tethers (e.g. a
    native :class:`ac.cell.compartments.MembraneCompartment`); this allocates the ``*_snap`` twins, the RNG
    epoch, and the Bell detach/attach event counters, and binds the source-gated ERM Bell constants.
    """
    n_erm = int(bound_d.shape[0])

    def z(n: int, dtype: object) -> wp.array:
        return wp.zeros(n, dtype=dtype, device=device)

    return ErmCortexConnector(
        membrane_idx_d=membrane_idx_d, cortex_idx_d=cortex_idx_d, bound_d=bound_d, rest_d=rest_d,
        rng_epoch_d=rng_epoch_d if rng_epoch_d is not None else z(1, wp.int32),
        bound_snap_d=z(n_erm, wp.int32), rest_snap_d=z(n_erm, wp.float64), rng_epoch_snap_d=z(1, wp.int32),
        detach_events_d=detach_events_d if detach_events_d is not None else z(1, wp.int32),
        attach_events_d=attach_events_d if attach_events_d is not None else z(1, wp.int32),
        membrane_pos_d=membrane_pos_d, cortex_pos_d=cortex_pos_d,
        k_erm=params.k_erm, k_on=params.k_on, k_off0=params.k_off0,
        bell_force_pn=params.bell_force_pn, capture_radius_um=params.capture_radius_um, device=device,
    )


# ── The wired slice: one accepted physical step over membrane ⟷ cortex through the ERM tether. ─────────
@dataclass(slots=True)
class ERMCortexSlice:
    """One vertical column of the cell graph wired for the accepted-step transaction.

    Holds the three participants (membrane + cortex state-owners and the ``membrane_erm_cortex`` connector),
    the device event clock, and the top-level :class:`CellTransaction`.  Constructed from already-assembled
    parts, so it is CPU-importable and the structural gates build it with recording doubles;
    :func:`build_erm_cortex_slice` assembles the real Warp-resident parts on the CUDA lane.

    Unlike the NMII slice there is NO separate ``compute_loads`` pass: the ERM Bell KMC recomputes each
    tether's tensile load from the converged geometry INSIDE its own accepted-gated kernel (ERM is simpler).
    """

    transaction: CellTransaction
    membrane: Any
    cortex: Any
    connector: ErmCortexConnector
    clock: Any
    ledger: Any | None = None
    tol_sq_d: wp.array | None = None

    def accumulate(self) -> None:
        """One candidate-force pass the caller's inner solve invokes each iteration (after zeroing forces).

        Scatters the two-array Bell tether load — ``-f`` into the membrane force array and ``+f`` into the
        cortex force array (Newton's 3rd law across two never-merged arrays).  No fixed prestress is injected;
        the tension is the reaction of whatever tethers the KMC has kept bound so far.
        """
        self.connector.accumulate_pair(
            self.membrane.position_d, self.membrane.force_d,
            self.cortex.position_d, self.cortex.force_d,
        )

    def step(
        self,
        solve: Callable[[], None],
        *,
        dt_phys: float,
        accepted_d: wp.array | None = None,
        ledger: Any | None = None,
        tol_sq_d: wp.array | None = None,
    ) -> None:
        """Run exactly one accepted physical step of the dynamic membrane--cortex ERM runtime.

        Args:
            solve: the caller-owned inner mechanical solve (zero forces → :meth:`accumulate` + surface passive
                relaxation → integrate to convergence).  Runs after events are proposed and owns all physics.
            dt_phys: the outer physical timestep [s]; finite and positive.
            accepted_d: the device acceptance predicate.  If ``None`` the balance gate decides (a ``ledger``
                with a ``tol_sq_d`` must then be available).
            ledger / tol_sq_d: optional balance-gate inputs; default to the slice's stored ones.

        The ERM Bell KMC commit reads the converged geometry directly, so no per-tether load pre-pass is
        needed.  The predicate stays device-resident throughout.
        """
        ledger = ledger if ledger is not None else self.ledger
        tol_sq_d = tol_sq_d if tol_sq_d is not None else self.tol_sq_d
        self.transaction.step(
            dt_phys=dt_phys, solve=solve, accepted_d=accepted_d, ledger=ledger, tol_sq_d=tol_sq_d,
        )

    def bound_tether_count(self) -> int:
        """Host read of the bound-tether population for OUT-OF-LOOP diagnostics only (device→host readback).

        A viz/telemetry helper the caller may call BETWEEN steps to watch the bound fraction settle under load;
        it must never be called inside the inner solve (I0-A: no authoritative GPU→CPU roundtrip in the loop).
        """
        return int(self.connector.bound_d.numpy().sum())


def build_erm_cortex_slice(
    *,
    membrane_owner: Any,
    cortex_owner: Any,
    connector: ErmCortexConnector,
    base_seed: int,
    device: str,
    clock: Any | None = None,
    ledger: Any | None = None,
    tol_sq_d: wp.array | None = None,
    population_ledgers: tuple[Any, ...] = (),
    architecture: CellArchitecture | None = None,
) -> ERMCortexSlice:
    """Assemble the dynamic ERM membrane--cortex slice: two surface owners + the ERM connector under one clock.

    Binds the membrane and cortex state-owners as the ``membrane``/``cortex`` components and the ERM connector
    as the ``membrane_erm_cortex`` connector on the reference architecture (which already carries the ERM
    edge), then wires a :class:`CellWorldTransaction` (``require_complete=False`` — this is a vertical slice,
    not the whole cell) and a :class:`CellTransaction`.  The ERM connector is the ONLY mechanical coupling.

    Args:
        membrane_owner / cortex_owner: the native ``membrane``/``cortex`` state-owners (see
            :func:`build_native_surface_owner`), or recording doubles for the CPU structural gates.
        connector: the :class:`ErmCortexConnector` owning the ERM SoA + Bell kinetics.
        base_seed: the run's fixed host RNG seed (nonnegative int).
        device: the resolved CUDA device string (never a hard-coded id).
        clock / ledger / tol_sq_d: optional; a :class:`WarpEventClock` is built from ``base_seed`` if omitted.
        population_ledgers: optional disjoint-population ledgers asserted at build + each accepted step.
        architecture: optional composition; defaults to :func:`reference_cell_architecture`.

    Returns:
        An :class:`ERMCortexSlice` whose :meth:`~ERMCortexSlice.step` runs one accepted physical step.
    """
    if architecture is None:
        architecture = reference_cell_architecture()

    actor = CellActor(architecture)
    actor.bind_component(MEMBRANE_COMPONENT, membrane_owner)
    actor.bind_component(CORTEX_COMPONENT, cortex_owner)
    actor.bind_connector(ERM_CONNECTOR, connector)
    world = CellWorldTransaction(actor, require_complete=False)

    if clock is None:
        from aleph.engine.events import make_event_clock  # CUDA-lane import (allocates device scalars).
        clock = make_event_clock(base_seed=base_seed, device=device)
    transaction = CellTransaction(world, clock, population_ledgers=tuple(population_ledgers))
    return ERMCortexSlice(
        transaction=transaction, membrane=membrane_owner, cortex=cortex_owner, connector=connector,
        clock=clock, ledger=ledger, tol_sq_d=tol_sq_d,
    )
