r"""``build_native_composed_cell_world`` — the REAL-owner analog of the census-double bridge.

WHAT THIS CLOSES.  :func:`~aleph.engine.composition.build_composed_cell_world` binds concrete runtime
objects into one ``(actor, pipeline, transaction)`` composed world, but its only concrete-object caller,
``scripts/ac_composed_world_dump.py``, feeds CPU census DOUBLES (``_CensusRuntime`` / ``_FacadeSpy``) — host
count carriers with the transaction API but **no geometry and no kernel** (see
``docs/v2_audit/cell_engine/NATIVE_COMPOSITION_SCOPE_2026-07-25.md``).  Each of the four advanced components
(cortex, ``sf_arc``, ECM, interior column) has its OWN per-slice native driver over a SEPARATE ``CellActor``;
nothing binds the REAL private-array owners into ONE composed ``CellTransaction`` that steps them together.

This module is the missing native assembly lane.  It composes the owners that already own **private device
arrays** and need **no frozen-driver edit** — ``sf_arc`` (:class:`SFArcStateOwner` over
:class:`aleph.engine.sf_mechanics.SFFilamentMechanics`, KERNEL_BOUND), ECM
(:class:`aleph.engine.ecm_world.ECMStateOwner`, owner-driven), and the head-resolved NMII actuator +
its ``nmii_cortex_motor`` connector (:mod:`aleph.engine.cortex_motor_slice`) — with the **cortex as a
bind-target PORT** (a :class:`~aleph.engine.nmii_actuator.FilamentMotorPortView`, never a driven
participant for the motor edge, exactly as ``cortex_motor_slice`` treats it).  They are registered into the
REAL :func:`build_composed_cell_world` runtimes as **real objects** — each launches its OWN Warp kernels over
its OWN disjoint population — and wrapped in ONE top-level
:class:`~aleph.engine.transaction.CellTransaction` (a
:class:`~aleph.engine.world.CellWorldTransaction` + a
:class:`~aleph.engine.events.WarpEventClock` + the per-component
:class:`~aleph.engine.population.PopulationLedger`\ s, with the disjoint no-double-count invariant
asserted at build and after every accepted step).  One accepted physical step then drives every owner under
one device-resident acceptance predicate.

HONESTLY EXCLUDED — the interior fluid column (membrane / cytosol / nucleus) is **NOT** composed here.  Its
couplings (``pressure`` / ``membrane_pressure`` / ``nucleus``) are computed inside the incumbent driver's
``_accumulate_all`` and moving them into engine connectors edits the **feature-frozen** ``ac/cell/driver.py``
— that is the PI-gated Card-5 strangler transition
(``INTERIOR_COLUMN_CONNECTED_PLAN_2026-07-25.md``), not this autonomous slice.  The five non-composed facade
types (surface_body / fluid / MT / IF / protrusion) are bound as documented dispatch-coverage placeholders so
the exact-once structural gate still fires; the composed STEP drives only the three real owners' own kernels.

TWO ENTRY POINTS (host-wiring vs CUDA lane):
  * :func:`compose_native_cell_world` — pure wiring.  Given the already-built real owners + clock + ledgers,
    it registers them into :func:`build_composed_cell_world` and wraps the returned world in one
    :class:`CellTransaction`.  It allocates no device memory and launches no kernel, so it is CPU-importable
    and the structural gate drives it with owner doubles (the CUDA lane injects the reals).
  * :func:`build_native_composed_cell_world` — the CUDA-lane builder.  It CONSTRUCTS the real ``sf_arc`` and
    ECM owners (via :func:`build_sf_arc_state_owner` / :func:`build_ecm_state_owner`) and the
    ``nmii_cortex_motor`` connector, then calls :func:`compose_native_cell_world`.  The NMII actuator
    state-owner + cortex port come from the caller (they depend on the cell-layer minifilament builder), so
    this stays engine-layer and never imports ``ac/cell``.

engine units: length µm, force pN, stiffness pN/µm.  Runtime: NVIDIA Warp on CUDA only (I0-A).

Sanity Gate (self-tested in tests/ac/engine/test_composed_native_smoke.py):
  * ownership: the composed world binds the REAL owners at the ``sf_arc`` / ``ecm`` / ``nmii`` component slots
    and the real connector at ``nmii_cortex_motor`` — verbatim, never wrapped in a census double; cortex is a
    PORT, never a transaction participant.
  * conservation / no-double-count: the cortex / ``sf_arc`` / ECM / NMII population ledgers are disjoint at
    build and re-asserted every accepted step (:func:`assert_disjoint_populations`).
  * one clock / one predicate: exactly one :class:`CellTransaction` wraps the deduplicated owners; a single
    device predicate reaches every ``rollback`` + ``commit_irreversible`` and the clock advances once.
  * honest exclusion: the interior-column fluid connectors are NOT bound (the PI-gated Card-5 seam).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import warp as wp

from aleph.engine.composition import (
    ComposedCellRuntimes,
    ComposedCellWorld,
    FacadeDispatch,
    build_composed_cell_world,
)
from aleph.engine.contracts import CellArchitecture, reference_cell_architecture
from aleph.engine.load_path import BorrowedSegmentActorView
from aleph.engine.dispatch import canonical_facade_claims
from aleph.engine.ecm_world import ECM_COMPONENT
from aleph.engine.nmii_actuator import CORTEX_COMPONENT, NMII_COMPONENT
from aleph.engine.population import PopulationLedger, assert_disjoint_populations
from aleph.engine.sf_mechanics import (
    SF_COMPONENT,
    SFFilamentMechanics,
    SFInternalArcJointConnector,
    build_sf_filament_mechanics,
    build_sf_internal_arc_connector,
    build_sf_mechanics_topology,
)
from aleph.engine.transaction import CellTransaction

__all__ = [
    "SF_COMPONENT",
    "ECM_COMPONENT",
    "NMII_COMPONENT",
    "NMII_CORTEX_MOTOR",
    "NMII_SF_MOTOR",
    "SF_CYTOSOL_TRANSFER",
    "NMII_CYTOSOL_TRANSFER",
    "FOCAL_ADHESION_COMPONENT",
    "INTERIOR_COLUMN_FLUID_CONNECTORS",
    "SFArcStateOwner",
    "ComposedNativeCell",
    "build_sf_arc_state_owner",
    "compose_native_cell_world",
    "build_native_composed_cell_world",
]

NMII_CORTEX_MOTOR = "nmii_cortex_motor"
SF_CORTEX_TRANSIENT = "sf_cortex_transient"
DORSAL_ARC_CROSSLINK = "dorsal_arc_crosslink"
NMII_SF_MOTOR = "nmii_sf_motor"
SF_CYTOSOL_TRANSFER = "sf_cytosol_transfer"
NMII_CYTOSOL_TRANSFER = "nmii_cytosol_transfer"
FOCAL_ADHESION_COMPONENT = "focal_adhesion"

#: The interior fluid column's three coupling connectors — HONESTLY EXCLUDED here (PI-gated Card-5 driver
#: seam; ``INTERIOR_COLUMN_CONNECTED_PLAN_2026-07-25.md``).  This slice never binds them.
INTERIOR_COLUMN_FLUID_CONNECTORS = (
    "surface_porous_transfer",
    "membrane_cytosol_boundary",
    "nucleus_cytosol_boundary",
)


# ── device kernels (import-safe on a CUDA-free host; compiled lazily only on a CUDA launch) ───────────────
@wp.kernel
def _zero_vec3d_kernel(force: wp.array(dtype=wp.vec3d)) -> None:
    """Zero one per-node force array (caller owns force zeroing before each candidate accumulate)."""
    force[wp.tid()] = wp.vec3d(0.0, 0.0, 0.0)


@wp.kernel
def _restore_vec3_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.vec3d),
    snap: wp.array(dtype=wp.vec3d),
) -> None:
    """Restore ``live <- snap`` iff the outer step was REJECTED (``accepted[0] == 0``; device-only branch)."""
    if accepted[0] == wp.int32(0):
        t = wp.tid()
        live[t] = snap[t]


@wp.kernel
def _bump_epoch_if_accepted_kernel(
    accepted: wp.array(dtype=wp.int32),
    epoch: wp.array(dtype=wp.int32),
) -> None:
    """Advance the accepted-step epoch iff ``accepted[0] != 0`` (device-only; no host predicate read)."""
    if accepted[0] != wp.int32(0):
        epoch[0] = epoch[0] + wp.int32(1)


# ── the REAL sf_arc state-owner (a transaction participant over the SF-owned passive rod-cable) ───────────
@dataclass(frozen=True, slots=True)
class SFArcStateOwner:
    """The ``sf_arc`` component as a composed-world transaction participant over its OWN device arrays.

    :class:`aleph.engine.sf_mechanics.SFFilamentMechanics` is the KERNEL_BOUND passive mechanics
    delegate (it LAUNCHES ``link_spring_kernel`` + ``cytosim_bending_kernel``); this owner adds the
    accepted-step transaction API the composed :class:`CellTransaction` requires — a device D2D snapshot of
    the SF node positions, a reject-gated restore, and an accepted-step epoch bump — over the SF-owned
    ``position_d`` / ``force_d``, exactly the way :class:`ECMStateOwner` / :class:`NMIIActuatorStateOwner`
    delegate their own accepted-step hooks.  No cortex node is ever referenced; ``sf_arc`` owns a DISJOINT
    filament population (PI 2026-07-22), so it snapshots/rolls-back/commits its OWN arrays only.

    The active contractile prestress is NOT lumped here — it enters through the separate ``nmii_sf_motor`` MOTOR
    connector, which landed on 2026-07-25 in :mod:`aleph.engine.sf_motor_slice` (this owner is reused there
    verbatim as the SF participant).  This owner itself carries only the passive backbone/bending laws + the
    transaction discipline; folding the SF motor edge into :func:`build_native_composed_cell_world` alongside the
    cortex motor edge additionally requires head-exclusivity between the two MOTOR connectors (one head may bind
    only one target), which is a separate composition step.
    """

    name: str
    device: str
    n_filaments: int
    n_nodes: int
    position_d: wp.array
    force_d: wp.array
    position_snap_d: wp.array
    epoch_d: wp.array
    epoch_snap_d: wp.array
    mechanics: SFFilamentMechanics
    arc_connector: SFInternalArcJointConnector | None = None
    launch: Callable[..., object] = wp.launch
    copy: Callable[..., object] = wp.copy

    def __post_init__(self) -> None:
        if self.name != SF_COMPONENT:
            raise ValueError(f"SF state owner must be {SF_COMPONENT!r}")
        if int(self.mechanics.n_nodes) != int(self.n_nodes):
            raise ValueError("SF mechanics node count must match the owned node count")

    def accumulate(self) -> None:
        """Launch the SF passive backbone/bending (+ internal arc crosslink) over the SF-owned arrays.

        Both kernels ``wp.atomic_add`` into ``force_d``; the caller owns force zeroing (the composed
        candidate solve zeroes every owner's force array before this pass).
        """
        self.mechanics.accumulate(self.position_d, self.force_d)
        if self.arc_connector is not None and self.arc_connector.n_joints:
            self.arc_connector.accumulate(self.position_d, self.force_d)

    def snapshot_candidate(self) -> None:
        """Snapshot the SF node positions + epoch before candidate mechanics (device D2D copy)."""
        self.copy(self.position_snap_d, self.position_d)
        self.copy(self.epoch_snap_d, self.epoch_d)

    def rollback(self, accepted: wp.array) -> None:
        """Restore the rejected candidate's SF positions in-device (accepted step keeps the converged ones)."""
        self.launch(
            _restore_vec3_if_rejected_kernel, dim=self.n_nodes,
            inputs=[accepted, self.position_d, self.position_snap_d], device=self.device,
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Advance the SF accepted-step epoch under the final predicate (device-gated)."""
        self.launch(
            _bump_epoch_if_accepted_kernel, dim=1, inputs=[accepted, self.epoch_d], device=self.device,
        )

    def accumulate_ledger(self, ledger: object) -> None:
        """Reduce this body's OWN force-array resultant into the balance gate's ``reaction`` channel.

        The ``sf_arc`` array is one of the two never-merged bodies of the ``nmii_sf_motor`` cut, and it is
        the one carrying the Dirichlet FA anchors, so it takes the reaction side (the assignment is fixed
        per lane; the gate is symmetric in the two).  Every SF-internal force cancels inside this
        resultant, leaving exactly what the connector scattered in from the ``nmii`` body — which the NMII
        owner reduces into the other channel, so the two must cancel by Newton's third law.  That is a
        check of the ADJOINT WIRING and not of convergence: it holds at any configuration and fails when a
        scatter is one-sided, sign-flipped, or double counted.

        A ledger exposing no such accessor is left alone (the delegate seam is exercised by recording
        doubles that own no accumulator), and the resulting one-sided gate does not fail silently: an
        unfilled channel makes ``|reaction + traction|`` equal to the OTHER channel's magnitude, which the
        relative tolerance cannot admit, so every step is rejected rather than waved through.
        """
        add_body_force = getattr(ledger, "add_body_force", None)
        if callable(add_body_force):
            add_body_force(self.force_d, side="reaction")

    def geometry_arrays(self) -> tuple[wp.array, wp.array]:
        """Return ``(position_d, force_d)`` for the composed candidate solve to drive/zero."""
        return (self.position_d, self.force_d)


def build_sf_arc_state_owner(
    population: Any,
    *,
    k_axial_pn_per_um: float,
    device: str,
    with_internal_arc: bool = True,
    launch: Callable[..., object] = wp.launch,
    topology: Any | None = None,
) -> SFArcStateOwner:
    """Build the REAL ``sf_arc`` state-owner from a disjoint :class:`SFArcPopulation` (CUDA lane).

    Uploads the SF node positions, allocates the SF-owned force + snapshot + epoch arrays, binds the real
    :class:`SFFilamentMechanics` (link + Cytosim bending kernels), and — when the population has internal
    dorsal↔arc joints — the genuinely KERNEL_BOUND :class:`SFInternalArcJointConnector`.  Sizes come from the
    passed population; a biological count is never lowered here (the sourced SF inventory arrives from
    ``CellState`` in a biology phase).

    Args:
        population: a built :class:`aleph.engine.sf_population.SFArcPopulation` (host NumPy).
        k_axial_pn_per_um: the actin axial backbone stiffness [pN/µm] — a REQUIRED modelling GAP (NF2007
            inextensible; no sourced ``EA_actin``); forwarded to :func:`build_sf_mechanics_topology`.
        device: the resolved CUDA device string (never a hard-coded id).
        with_internal_arc: bind the internal dorsal↔arc crosslink connector when the population has joints.
        launch: injected launcher (defaults to ``wp.launch``; the CUDA-free gate substitutes a recorder).
        topology: an ALREADY-built :class:`SFMechanicsTopology` to reuse.  Pass it when the caller also needs
            the topology (e.g. to build the ``sf_arc`` motor port from the same per-segment metadata) so the
            segment ordering has exactly ONE source of truth; omit it to build one here.

    Returns:
        The wired :class:`SFArcStateOwner` (a composed-world transaction participant).
    """
    import numpy as np

    if topology is None:
        topology = build_sf_mechanics_topology(population, k_axial_pn_per_um=k_axial_pn_per_um)
    elif float(topology.k_axial_pn_per_um) != float(k_axial_pn_per_um):
        raise ValueError(
            f"the supplied SF topology was built with k_axial={topology.k_axial_pn_per_um!r} pN/µm but this "
            f"owner was asked for {k_axial_pn_per_um!r}; pass one value, not two"
        )
    mechanics = build_sf_filament_mechanics(topology, device=device, launch=launch)
    arc_connector = None
    if with_internal_arc and topology.n_arc_joints:
        arc_connector = build_sf_internal_arc_connector(topology, device=device, launch=launch)

    pos = np.ascontiguousarray(population.pos, dtype=np.float64)
    n_nodes = int(pos.shape[0])
    position_d = wp.array(pos, dtype=wp.vec3d, device=device)
    return SFArcStateOwner(
        name=SF_COMPONENT, device=str(device),
        n_filaments=int(population.ledger.active_count), n_nodes=n_nodes,
        position_d=position_d,
        force_d=wp.zeros(n_nodes, dtype=wp.vec3d, device=device),
        position_snap_d=wp.zeros(n_nodes, dtype=wp.vec3d, device=device),
        epoch_d=wp.zeros(1, dtype=wp.int32, device=device),
        epoch_snap_d=wp.zeros(1, dtype=wp.int32, device=device),
        mechanics=mechanics, arc_connector=arc_connector, launch=launch,
    )


# ── dispatch-coverage facades (the exact-once structural gate; NOT transaction participants) ──────────────
def _excluded_facade_noop(*_args: object, **_kwargs: object) -> None:
    """Coverage placeholder for a facade type NOT composed in this native slice.

    The composed STEP drives only the three real owners' own kernels (through the caller solve); this no-op
    exists solely so :func:`build_composed_cell_world`'s exact-once dispatch-coverage gate still fires for the
    honestly-excluded facade types (surface_body / fluid / MT / IF / protrusion).  It launches nothing.
    """
    return None


@dataclass(slots=True)
class _NativeMechanicsFacade:
    """One dispatch facade whose canonical method drives a REAL owner's mechanics (or a coverage no-op).

    For the three composed owner types (``StressFiberActor`` → ``sf_arc``, ``ECMWorld`` → ``ecm``,
    ``NMIIActuator`` → ``nmii``) ``run`` launches the real owner's own kernels; for the five excluded types it
    is :func:`_excluded_facade_noop`.  All three canonical method names are exposed so one class serves every
    facade type; :func:`build_composed_cell_world` invokes exactly the one method the manifest claims.
    """

    run: Callable[[], None]
    label: str

    def accumulate_mechanics(self, *args: object, **kwargs: object) -> None:
        self.run()

    def accumulate_candidate(self, *args: object, **kwargs: object) -> None:
        self.run()

    def candidate_iteration(self, *args: object, **kwargs: object) -> None:
        self.run()


# ── the composed native cell: one CellActor + one CellTransaction over the real owners ────────────────────
@dataclass(slots=True)
class ComposedNativeCell:
    """A composed native cell: the real owners bound into one world + one accepted-step transaction.

    Args:
        world: the :class:`ComposedCellWorld` (actor + exact-once pipeline + ``CellWorldTransaction``).
        transaction: the ONE top-level :class:`CellTransaction` (world + clock + population ledgers).
        component_owners / connector_runtimes: the REAL owners bound (by architecture name).
        cortex_port: the cortex bind-target PORT (never a transaction participant), or ``None``.
        clock: the device event clock (advanced once per accepted step).
        population_ledgers: the disjoint per-component ledgers (asserted every accepted step).
    """

    world: ComposedCellWorld
    transaction: CellTransaction
    component_owners: dict[str, object]
    connector_runtimes: dict[str, object]
    cortex_port: object | None
    clock: Any
    population_ledgers: tuple[PopulationLedger, ...]
    _force_arrays: tuple[wp.array, ...] = field(default_factory=tuple)
    _load_compute: Callable[[], None] | None = None
    _launch: Callable[..., object] = wp.launch
    #: The incumbent inner solve that owns the interior column's mechanics (Card-5), or ``None``. It runs
    #: FIRST in the candidate pass because it zeroes the global force array; see ``solve_candidate``.
    _interior_solve: Callable[[], None] | None = None
    #: Diagnostics for the `nmii_sf_motor` head partition, or ``{'bound': False}``.
    sf_motor_report: dict = field(default_factory=lambda: {"bound": False})
    #: (component name, its OverdampedStep, its position array) for each component that advances
    #: its OWN positions. Empty means the world is force-assembly only, which is what it was.
    _relaxers: tuple[tuple[str, object, wp.array], ...] = field(default_factory=tuple)

    def borrowed_pairs_for(self, component: str) -> tuple[tuple[Any, Any], ...]:
        """Stiffness every BOUND connector scatters into ``component``'s force array, in its indices.

        A component's own topology cannot see this: `sf_cortex_transient` borrows
        ``sf_owner.position_d``/``force_d`` through :class:`BorrowedSegmentActorView` and owns no
        array here, so a topology scan is structurally blind to its 4.6e5 pN/µm. Feed the result to
        ``build_overdamped_step(..., borrowed_pairs=...)`` — that argument has no default precisely so
        this cannot be skipped.

        CANDIDATE pairs, not engaged ones: engagement is transient and a bound that changes with
        state is not a bound.  Each joint is charged to BOTH endpoints of its segment at full ``k``,
        which is conservative — the scatter weights are ``(1−u, u) ≤ 1``.

        Raises:
            ValueError: if a connector is bound whose stiffness cannot be read. Reporting `()` there
                would hand back an optimistic bound, which is the defect this exists to prevent.
        """
        owner = self.component_owners[component]
        pairs: list[tuple[Any, Any]] = []
        for name, connector in self.connector_runtimes.items():
            joints = getattr(connector, "_joints", connector)
            for side in ("actor_a", "actor_b"):
                actor = getattr(joints, side, None)
                if actor is None or getattr(actor, "force_d", None) is not owner.force_d:
                    continue
                letter = side[-1]
                elements = getattr(joints, f"element_{letter}_id_d", None)
                stiffness = getattr(joints, "stiffness_d", None)
                segments = getattr(actor, "segments_d", None)
                if elements is None or stiffness is None or segments is None:
                    raise ValueError(
                        f"connector {name!r} scatters into {component!r}'s force array but its "
                        f"stiffness cannot be read, so no bound covering it can be derived"
                    )
                pairs.append(
                    (segments.numpy().reshape(-1, 2)[elements.numpy()], stiffness.numpy())
                )
        return tuple(pairs)

    def attach_relaxer(self, component: str, step: object, position_d: wp.array) -> None:
        """Let ``component`` advance its OWN positions inside the accepted step.

        OPT-IN, and deliberately not done at build: a world with no relaxer behaves exactly as it did
        before this existed, so the three instruments that call ``solve_candidate()`` bare are
        untouched and the relax can be measured against its own absence in the same run.

        Raises:
            KeyError: if ``component`` owns nothing here — a relaxer on an unbound name would move
                nothing and report success.
            ValueError: if ``step`` was derived for a different node count than the array it is
                being pointed at, or if the component is already registered.
        """
        owner = self.component_owners[component]
        if any(name == component for name, _, _ in self._relaxers):
            raise ValueError(f"{component!r} already advances its own positions")
        if int(getattr(step, "n_nodes", -1)) != int(position_d.shape[0]):
            raise ValueError(
                f"step derived for {getattr(step, 'n_nodes', None)} nodes, array has "
                f"{int(position_d.shape[0])} — a mismatched relaxer would move the wrong nodes"
            )
        if int(owner.force_d.shape[0]) != int(position_d.shape[0]):
            raise ValueError(f"{component!r} force and position arrays disagree in length")
        self._relaxers = (*self._relaxers, (component, step, position_d))

    def subcycles_for(self, dt_phys: float) -> int:
        """How many relax substeps this world's stiffest relaxing component demands for ``dt_phys``.

        DERIVED, never typed: ``ceil(dt / dt_max)`` from each component's own Gershgorin bound.  For
        ``sf_arc`` at native this is large — α-actinin's 4.6e5 pN/µm is the fastest mode in the graph
        — and that is the measurement, not a problem to tune away.  Reporting 1 when the bound says
        otherwise would be integrating an unstable configuration quietly.
        """
        if not self._relaxers or dt_phys <= 0.0:
            return 1
        worst = min(float(step.dt_max_s) for _, step, _ in self._relaxers)
        return max(int(math.ceil(dt_phys / worst)), 1)

    def _solve_and_relax(self, dt_phys: float) -> None:
        """Reassemble forces and advance the relaxing components' own positions, ``n`` times.

        Forces MUST be re-evaluated at each substep: holding them fixed across the interval is one
        large explicit step wearing a subcycle's name, and the derived bound would no longer apply.

        The trailing ``solve_candidate`` is not optional bookkeeping.  ``accumulate_ledger`` reduces
        ``force_d``, so without a final assembly AT THE MOVED GEOMETRY the acceptance predicate is
        bit-identical with and without the relax — a divergent move could not be rejected, and the
        gate would pass by being blind rather than by being satisfied.
        """
        n = self.subcycles_for(dt_phys)
        sub_dt = dt_phys / n
        for _ in range(n):
            self.solve_candidate()
            for name, step, position_d in self._relaxers:
                step.advance(position_d, self.component_owners[name].force_d, sub_dt)
        self.solve_candidate()

    @property
    def participants(self) -> tuple[object, ...]:
        """The deduplicated accepted-step participants (the real state-owners + the motor connector)."""
        return self.transaction.world.participants

    def solve_candidate(self) -> None:
        """One composed candidate-force pass: every real owner launches its OWN kernels, disjointly.

        Zeroes each owner's force array, then dispatches the exact-once pipeline (the real facades launch the
        ``sf_arc`` / ECM passive mechanics and the NMII internal mechanics + split crossbridge into the cortex
        port), then computes the per-head motor loads on the converged geometry before commit consumes them.
        Cross-component coupling happens ONLY through the connector scattering into the borrowed cortex
        ``force_d`` (Newton's 3rd law, two never-merged arrays) — never a shared node array.
        """
        # ORDER IS THE CORRECTNESS. The interior solve's first act is to zero the ONE global force array, so
        # it must run BEFORE the engine's own zeroing — otherwise whichever runs second erases the other's
        # forces, silently and with no exception. The engine's zero list excludes the cortex port for exactly
        # this reason (see `compose_native_cell_world`), and the motor's scatter then ADDS onto the cortex
        # block the interior solve just filled.
        if self._interior_solve is not None:
            self._interior_solve()
        for force in self._force_arrays:
            self._launch(_zero_vec3d_kernel, dim=int(force.shape[0]), inputs=[force], device=str(force.device))
        self.world.pipeline.run_candidate()
        if self._load_compute is not None:
            self._load_compute()

    def step_once(
        self,
        *,
        dt_phys: float,
        accepted_d: wp.array | None = None,
        ledger: Any | None = None,
        tol_sq_d: wp.array | None = None,
        rates: object | None = None,
        neighbors: object | None = None,
    ) -> None:
        """Run exactly one accepted physical step of the composed cell under ONE device predicate.

        ``solve`` is a zero-argument closure, not a changed ``solve_candidate`` signature: three
        instruments call ``solve_candidate()`` bare and outside any transaction, and a parameter
        would ``TypeError`` all three on the GPU host while a default ``dt`` would be a chosen
        constant.  With no relaxer registered the closure IS ``solve_candidate``, unchanged.
        """
        solve = self.solve_candidate if not self._relaxers else (
            lambda: self._solve_and_relax(dt_phys)
        )
        self.transaction.step(
            dt_phys=dt_phys, solve=solve, accepted_d=accepted_d,
            ledger=ledger, tol_sq_d=tol_sq_d, rates=rates, neighbors=neighbors,
        )


def compose_native_cell_world(
    *,
    sf_owner: Any,
    ecm_owner: Any,
    nmii_actuator: Any,
    nmii_cortex_connector: Any,
    cortex_port: Any,
    clock: Any,
    population_ledgers: tuple[PopulationLedger, ...],
    sf_internal_arc_connector: Any | None = None,
    sf_cortex_connector: Any | None = None,
    sf_motor_connector: Any | None = None,
    solid_cytosol_transfers: dict[str, Any] | None = None,
    fa_owner: Any | None = None,
    interior_column: Any | None = None,
    interior_solve: Callable[[], None] | None = None,
    architecture: CellArchitecture | None = None,
    launch: Callable[..., object] = wp.launch,
) -> ComposedNativeCell:
    """Register the REAL owners into :func:`build_composed_cell_world` and wrap ONE :class:`CellTransaction`.

    Pure wiring: allocates no device memory and launches no kernel, so it is CPU-importable and the
    structural gate drives it with owner doubles (the CUDA lane injects the reals).  The three real owners are
    bound VERBATIM at the ``sf_arc`` / ``ecm`` / ``nmii`` component slots and the real connector at
    ``nmii_cortex_motor`` — never wrapped in a census double.  The **cortex is a PORT**: it is NOT bound as a
    component, so it is never a transaction participant; it enters only through the connector's adjoint force
    scatter into ``cortex_port.force_d``.  The interior fluid column is NOT bound (PI-gated Card-5 seam).

    Args:
        sf_owner: the real :class:`SFArcStateOwner` (bound at ``sf_arc``).
        ecm_owner: the real :class:`ECMStateOwner` (bound at ``ecm``).
        nmii_actuator: the real :class:`NMIIActuatorStateOwner` (bound at ``nmii``).
        nmii_cortex_connector: the real ``nmii_cortex_motor`` MOTOR connector (a participant + event proposer).
        cortex_port: the cortex :class:`FilamentMotorPortView` bind target (PORT only, never a participant).
        clock: the device event clock (a :class:`WarpEventClock` on the CUDA lane, a spy on the gate).
        population_ledgers: the disjoint per-component ledgers (asserted at build + every accepted step).
        sf_internal_arc_connector: optional internal dorsal↔arc connector, driven in ``solve_candidate``
            (the SF owner already drives it; passing it separately is not required).
        architecture: optional composition (defaults to :func:`reference_cell_architecture`).
        launch: injected launcher for the composed candidate zeroing (defaults to ``wp.launch``).

    Returns:
        A :class:`ComposedNativeCell` whose :meth:`ComposedNativeCell.step_once` runs one accepted step.
    """
    architecture = architecture or reference_cell_architecture()

    component_owners: dict[str, object] = {
        SF_COMPONENT: sf_owner,
        ECM_COMPONENT: ecm_owner,
        NMII_COMPONENT: nmii_actuator,
    }
    if fa_owner is not None:
        component_owners[FOCAL_ADHESION_COMPONENT] = fa_owner
    connector_runtimes: dict[str, object] = {NMII_CORTEX_MOTOR: nmii_cortex_connector}
    # `dorsal_arc_crosslink` was BUILT all along and never REGISTERED.  `build_sf_arc_state_owner` binds it
    # by default and `SFArcStateOwner.accumulate` already launches its `link_spring_kernel`, so the joint
    # force has been in every composed run — attributed to the `sf_arc` COMPONENT, because the actor was
    # never told the connector exists.  The 2026-08-11 device-run census made the gap visible
    # (`outputs/ac/connector_devicerun/`): it reported the edge `NOT_BOUND` while its kernel was running.
    # Registering it changes NO physics — same object, same launch, now under its declared name — and the
    # census re-run is the check.  Two slots were already reserved for this and left unfilled: the
    # `DORSAL_ARC_CROSSLINK` constant above, and this function's own `sf_internal_arc_connector` parameter,
    # which was accepted and documented but read by nothing.
    if sf_internal_arc_connector is not None and getattr(sf_internal_arc_connector, "n_joints", 0):
        connector_runtimes[DORSAL_ARC_CROSSLINK] = sf_internal_arc_connector
    # `sf_cortex_transient` is bound only when the sourced capture radius actually selected joints;
    # the builder returns None otherwise, and None must not become an empty connector.
    if sf_cortex_connector is not None:
        connector_runtimes[SF_CORTEX_TRANSIENT] = sf_cortex_connector
    # `nmii_sf_motor`: an INTER-component MOTOR edge, so like `sf_cortex_transient` nothing else would
    # call it. Job 94 measured exactly this gap — a bound connector that is never accumulated carries
    # no force and the cut still reads as wired.
    if sf_motor_connector is not None:
        connector_runtimes[NMII_SF_MOTOR] = sf_motor_connector
    for _edge, _transfer in (solid_cytosol_transfers or {}).items():
        connector_runtimes[_edge] = _transfer

    # CARD-5 (PI-approved 2026-08-11 15:24 KST). The interior fluid column — membrane / cortex / cytosol /
    # nucleus and the three fluid connectors — was "honestly excluded" from this world because its coupling
    # is computed inside the feature-frozen `driver.py::_accumulate_all`. It enters here as built parts plus
    # the incumbent solve that owns their mechanics, with D2's `omit=` keeping the two from double-counting.
    #
    # CORTEX IS PROMOTED FROM PORT TO COMPONENT. It stays the motor's bind target — same memory, two access
    # paths — but it is now also a transaction participant with its own snapshot/rollback, which is what the
    # remaining cortex-endpoint connectors need in order to be bindable at all.
    if interior_column is not None:
        from aleph.engine.interior_column_slice import (
            CYTOSOL_COMPONENT, MEMBRANE_COMPONENT, MEMBRANE_FLUID_BOUNDARY,
            NUCLEUS_COMPONENT, NUCLEUS_FLUID_BOUNDARY, SURFACE_POROUS_TRANSFER,
        )
        component_owners[MEMBRANE_COMPONENT] = interior_column.membrane
        component_owners[CORTEX_COMPONENT] = interior_column.cortex
        component_owners[CYTOSOL_COMPONENT] = interior_column.cytosol
        component_owners[NUCLEUS_COMPONENT] = interior_column.nucleus
        connector_runtimes[SURFACE_POROUS_TRANSFER] = interior_column.cortex_transfer
        connector_runtimes[MEMBRANE_FLUID_BOUNDARY] = interior_column.membrane_boundary
        connector_runtimes[NUCLEUS_FLUID_BOUNDARY] = interior_column.nucleus_boundary
        if interior_solve is None:
            raise ValueError(
                "an interior column without its solve would be a PASSENGER: four owners whose arrays "
                "nothing relaxes, which is precisely the T2 negative control. Pass the incumbent "
                "inner solve (with `omit=` for whatever the engine now supplies)."
            )

    # Real dispatch facades drive the real owners' own kernels; excluded facade types are coverage no-ops.
    def _sf_run() -> None:
        sf_owner.accumulate()
        # `sf_cortex_transient` is an INTER-component edge, so unlike the internal arc joint (which the
        # SF owner's own accumulate launches) nothing else would call it.  Registering it in
        # `connector_runtimes` makes the actor know it exists; it does not put it in the candidate
        # pass.  Job 94 measured exactly that gap — the cut stopped reporting `unwired` and still
        # cancelled nothing, because a bound connector that is never accumulated carries no force.
        if sf_cortex_connector is not None:
            sf_cortex_connector.accumulate()
        if sf_motor_connector is not None:
            sf_motor_connector.accumulate_candidate(nmii_actuator.geometry(), sf_motor_connector.port)
        _sf_transfer = (solid_cytosol_transfers or {}).get(SF_CYTOSOL_TRANSFER)
        if _sf_transfer is not None:
            _sf_transfer.accumulate(sf_owner.position_d, sf_owner.force_d)

    def _ecm_run() -> None:
        ecm_owner.accumulate()

    def _nmii_run() -> None:
        nmii_actuator.accumulate_internal()
        _nmii_transfer = (solid_cytosol_transfers or {}).get(NMII_CYTOSOL_TRANSFER)
        if _nmii_transfer is not None:
            _nmii_transfer.accumulate(nmii_actuator.position_d, nmii_actuator.force_d)
        nmii_cortex_connector.accumulate_candidate(nmii_actuator.geometry(), cortex_port)

    from aleph.engine.ecm_world import ECMWorld
    from aleph.engine.nmii_actuator import NMIIActuator
    from aleph.engine.stress_fiber import StressFiberActor

    real_runs: dict[type[object], Callable[[], None]] = {
        StressFiberActor: _sf_run,
        ECMWorld: _ecm_run,
        NMIIActuator: _nmii_run,
    }
    facade_dispatch: dict[type[object], FacadeDispatch] = {}
    for claim in canonical_facade_claims():
        run = real_runs.get(claim.owner_type, _excluded_facade_noop)
        facade_dispatch[claim.owner_type] = FacadeDispatch(
            _NativeMechanicsFacade(run, claim.owner_type.__name__), args=(),
        )

    runtimes = ComposedCellRuntimes(
        component_owners=component_owners,
        connector_runtimes=connector_runtimes,
        facade_dispatch=facade_dispatch,
    )
    # require_complete=False: this is a bring-up composed world (three owners + the motor connector), not the
    # whole 13/32 cell.  The exact-once dispatch-coverage gate still fires (all eight facades present).
    world = build_composed_cell_world(runtimes, architecture, require_complete=False)

    # ONE top-level transaction: the CellWorldTransaction + the device clock + the disjoint population ledgers.
    population_ledgers = tuple(population_ledgers)
    assert_disjoint_populations(population_ledgers)  # build-time no-double-count (also re-run every step)
    transaction = CellTransaction(world.transaction, clock, population_ledgers=population_ledgers)

    # Force arrays the composed candidate zeroes before each accumulate (owner-private + the cortex port).
    #
    # With an interior column bound the cortex port is REMOVED from this list, and the omission is the whole
    # correctness of the ordering. The incumbent's `_accumulate_all` zeroes the ONE global force array — which
    # contains the cortex block — as its first act. Zeroing the port again afterwards would erase the interior
    # solve's cortex physics; zeroing it before is redundant. So the sequence is: interior solve (zeroes the
    # global array, fills membrane/cortex/nucleus/cytosol), then zero the engine's OWN private arrays, then
    # the engine pipeline, whose motor scatter ADDS onto the cortex block rather than replacing it.
    engine_private = (sf_owner, ecm_owner, nmii_actuator)
    zeroed_owners = engine_private if interior_column is not None else (*engine_private, cortex_port)
    force_arrays: list[wp.array] = []
    for owner in zeroed_owners:
        force = getattr(owner, "force_d", None)
        if force is not None:
            force_arrays.append(force)

    load_compute: Callable[[], None] | None = None
    compute_loads = getattr(nmii_cortex_connector, "compute_loads", None)
    if callable(compute_loads):
        def load_compute() -> None:
            compute_loads(nmii_actuator.geometry(), cortex_port)

    return ComposedNativeCell(
        world=world,
        transaction=transaction,
        component_owners=component_owners,
        connector_runtimes=connector_runtimes,
        cortex_port=cortex_port,
        clock=clock,
        population_ledgers=population_ledgers,
        _force_arrays=tuple(force_arrays),
        _load_compute=load_compute,
        _launch=launch,
        _interior_solve=interior_solve,
    )


def build_native_composed_cell_world(
    *,
    device: str,
    base_seed: int,
    sf_population: Any,
    sf_k_axial_pn_per_um: float,
    ecm_topology: Any,
    ecm_settings: Any,
    ecm_spec: Any | None,
    nmii_actuator: Any,
    cortex_port: Any,
    nmii_cortex_connector: Any,
    cortex_capacity: int = 70_686,
    sf_motor_params: Any | None = None,
    cortex_positions_um: Any | None = None,
    cytosol_grid_dx_um: float | None = None,
    cortex_active_ids: range | None = None,
    ecm_id_base: int = 2_000_000,
    nmii_id_base: int = 3_000_000,
    interior_column: Any | None = None,
    interior_solve: Callable[[], None] | None = None,
    architecture: CellArchitecture | None = None,
) -> ComposedNativeCell:
    """CUDA-lane builder: construct the real ``sf_arc`` + ECM owners + ledgers + clock, then compose.

    Builds the ``sf_arc`` state-owner (:func:`build_sf_arc_state_owner`) and the ECM state-owner
    (:func:`aleph.engine.ecm_world.build_ecm_state_owner`) — each owns its own private device arrays and
    launches its own kernels — plus the per-component population ledgers and a
    :class:`WarpEventClock`, then calls :func:`compose_native_cell_world`.  The NMII actuator state-owner, the
    cortex bind-target port, and the ``nmii_cortex_motor`` connector are supplied by the caller (they depend
    on the cell-layer minifilament/straddle builder, so they are built in the gbook driver and injected here —
    keeping this module engine-layer, never importing ``ac/cell``).

    Populations are DISJOINT by construction: cortex ``[0, cortex_capacity)``, ``sf_arc`` its own
    ``id_base=1_000_000`` block (from the population ledger), ECM ``[ecm_id_base, …)``, NMII
    ``[nmii_id_base, …)`` — asserted at build and after every accepted step.

    Args:
        device: the resolved CUDA device string.
        base_seed: the run's fixed host RNG seed (nonnegative int) for the :class:`WarpEventClock`.
        sf_population: a built :class:`SFArcPopulation`.
        sf_k_axial_pn_per_um: the SF axial backbone stiffness GAP [pN/µm].
        ecm_topology: a built :class:`ECMTopologyState` (device SoA).
        ecm_settings: an :class:`ECMWorldSettings` (boundary anchor mode).
        ecm_spec: the collagen material spec (or ``None`` for the builder default).
        nmii_actuator: the real :class:`NMIIActuatorStateOwner` (caller-built).
        cortex_port: the cortex :class:`FilamentMotorPortView` (caller-built).
        nmii_cortex_connector: the real ``nmii_cortex_motor`` connector (caller-built).
        cortex_capacity / cortex_active_ids: the cortex population block + seeded active slice.
        ecm_id_base / nmii_id_base: the ECM / NMII host-mirror ledger blocks (disjoint from cortex + sf_arc).
        architecture: optional composition (defaults to :func:`reference_cell_architecture`).

    Returns:
        A :class:`ComposedNativeCell` ready for :meth:`ComposedNativeCell.step_once`.
    """
    from aleph.engine.ecm_world import build_ecm_state_owner
    from aleph.engine.events import make_event_clock
    from aleph.engine.load_path import ActorRecord

    # 1. REAL sf_arc owner (private device arrays; launches link + Cytosim-bending kernels).
    sf_owner = build_sf_arc_state_owner(
        sf_population, k_axial_pn_per_um=sf_k_axial_pn_per_um, device=device,
    )

    # 2. REAL ECM owner (its own device SoA; launches the collagen constitutive force).
    n_segments = int(ecm_topology.n_segment_capacity)
    ecm_owner = build_ecm_state_owner(
        ecm_topology,
        settings=ecm_settings,
        actor=ActorRecord(ECM_COMPONENT, 0, 0, 0, 0, n_segments),
        transaction=_ECMOwnerTransaction(),
        ledger=_ECMOwnerLedger(),
        spec=ecm_spec,
        device=device,
    )

    # 2b. sf_cortex_transient — the only unwired power port whose pairing rule is already SOURCED
    # (alpha-actinin capture radius 0.06 um + link_k 4.6e5 pN/um, Ferrer 2008 / PI 2026-06-30), so its
    # topology is DERIVED from geometry rather than chosen.  The 2026-08-11 interface residual measured
    # this cut carrying force on one side with nothing crossing it; this is what crosses.
    #
    # It returns None when the geometry supplies no pair inside the sourced radius, and None must stay
    # None: a connector with no joints cannot fail to balance, and registering one would manufacture
    # exactly the vacuous PASS the residual reading exists to expose.  The radius is never widened to
    # produce pairs — that would be tuning a sourced constant to an outcome.
    sf_cortex_connector = None
    if cortex_port is not None:
        from aleph.engine.load_path import ActorRegistry
        from aleph.engine.alpha_actinin_kinetics import build_alpha_actinin_kinetics
        from aleph.engine.sf_cortex_transient import (
            build_sf_cortex_transient_connector,
            pair_sf_to_cortex_segments,
        )

        sf_links = sf_owner.mechanics.links_d.numpy()
        cortex_segments = np.stack(
            [cortex_port.segment_node_a_d.numpy(), cortex_port.segment_node_b_d.numpy()], axis=1,
        )
        pairing = pair_sf_to_cortex_segments(
            sf_positions=sf_owner.position_d.numpy(), sf_segments=sf_links,
            cortex_positions=cortex_port.position_d.numpy(), cortex_segments=cortex_segments,
        )
        sf_seg_record = ActorRecord(SF_COMPONENT, 10, 0, 0, 0, max(int(sf_links.shape[0]), 1))
        cortex_seg_record = ActorRecord(
            CORTEX_COMPONENT, 11, 0, 0, 0, max(int(cortex_segments.shape[0]), 1)
        )
        # BOUND (2026-08-11). The edge declares kinetics=True, so it needs a rate law or it is a
        # PERMANENT WELD — the charter forbids that by name and `build_sf_cortex_transient_connector`
        # refuses a None kinetics. What was missing was a DELEGATE, never a datum: every alpha-actinin
        # constant is sourced (k_on 10/s, k_off0 0.066/s, Bell f0, capture 0.06 um, link_k 4.6e5
        # pN/um; Ferrer 2008, PI-approved 2026-06-30), so writing it cost no decision.
        #
        # None stays None when the geometry supplies no pair: a connector with no joints cannot fail
        # to balance, and registering one would manufacture the vacuous PASS the interface-residual
        # reading exists to expose.
        sf_cortex_connector = build_sf_cortex_transient_connector(
            pairing=pairing,
            sf_actor=BorrowedSegmentActorView(
                sf_seg_record, sf_owner.position_d, sf_owner.force_d,
                wp.array(np.ascontiguousarray(sf_links, dtype=np.int32), dtype=wp.int32, device=device),
            ),
            cortex_actor=BorrowedSegmentActorView(
                cortex_seg_record, cortex_port.position_d, cortex_port.force_d,
                wp.array(np.ascontiguousarray(cortex_segments, dtype=np.int32), dtype=wp.int32,
                         device=device),
            ),
            registry=ActorRegistry((sf_seg_record, cortex_seg_record)),
            sf_record=sf_seg_record, cortex_record=cortex_seg_record,
            kinetics=build_alpha_actinin_kinetics(),
            architecture=architecture,
        )

    # 3. disjoint per-component population ledgers (host mirror of the device free-lists).
    cortex_ledger = PopulationLedger("cortex", id_base=0, capacity=int(cortex_capacity))
    cortex_ledger.seed_active(cortex_active_ids if cortex_active_ids is not None else range(0, 128))
    ecm_ledger = PopulationLedger(ECM_COMPONENT, id_base=int(ecm_id_base), capacity=max(n_segments, 1))
    ecm_ledger.seed_active(range(int(ecm_id_base), int(ecm_id_base) + max(n_segments, 1)))
    nmii_n = int(nmii_actuator.n_minifilaments)
    nmii_ledger = PopulationLedger(NMII_COMPONENT, id_base=int(nmii_id_base), capacity=max(nmii_n, 1))
    nmii_ledger.seed_active(range(int(nmii_id_base), int(nmii_id_base) + max(nmii_n, 1)))
    population_ledgers = (cortex_ledger, sf_population.ledger, ecm_ledger, nmii_ledger)

    # 4. ONE device event clock, advanced only on an accepted step.
    clock = make_event_clock(base_seed=int(base_seed), device=device)

    # `nmii_sf_motor` — the one power port the 2026-08-11 interface residual reported UNWIRED, and the
    # reason `sf_arc` carries no load. The connector class is already TARGET-AGNOSTIC (its own docstring
    # names this edge); what was missing was the `sf_arc` PORT and, behind it, head exclusivity. One
    # myosin head engages ONE actin filament, so the two motor edges get DISJOINT head sets, partitioned
    # by which target each minifilament's centre is nearer — a distance ORDERING, so no radius is chosen.
    sf_motor_connector = None
    sf_motor_report: dict[str, object] = {"bound": False}
    if sf_motor_params is not None and cortex_positions_um is not None:
        from aleph.engine.cortex_motor_slice import (
            CortexMotorConnector, SegmentAttachQuery, allocate_segment_connector_state,
        )
        from aleph.engine.nmii_sf_motor import (
            HeadMaskedAttachQuery, build_sf_motor_port,
            partition_minifilaments_by_proximity, serves_mask,
        )
        from aleph.engine.sf_mechanics import build_sf_mechanics_topology

        topology = build_sf_mechanics_topology(sf_population, k_axial_pn_per_um=sf_k_axial_pn_per_um)
        sf_port = build_sf_motor_port(sf_owner, topology, device=device)
        n_heads = int(nmii_actuator.n_heads)
        partition = partition_minifilaments_by_proximity(
            head_node_d=nmii_actuator.head_node_d,
            head_positions_um=nmii_actuator.position_d.numpy()[nmii_actuator.head_node_d.numpy()],
            cortex_positions_um=cortex_positions_um,
            sf_positions_um=sf_owner.position_d.numpy(),
            heads_per_minifilament=n_heads // max(int(nmii_actuator.n_minifilaments), 1),
        )
        # A connector with no heads cannot fail to balance; registering one would manufacture the
        # vacuous PASS the interface-residual instrument exists to expose. None stays None.
        if partition.sf_heads.size:
            state = allocate_segment_connector_state(
                n_heads, int(sf_port.segment_node_a_d.shape[0]), device=device)
            sf_motor_connector = CortexMotorConnector(
                state=state, head_node_d=nmii_actuator.head_node_d,
                params=sf_motor_params.to_hand_params(), catch_slip=None,
                query=HeadMaskedAttachQuery(
                    inner=SegmentAttachQuery(
                        state=state, head_node_d=nmii_actuator.head_node_d, port=sf_port,
                        capture_radius=sf_motor_params.capture_radius,
                        max_segment_length_um=float(np.max(np.asarray(topology.link_r0))),
                        device=device),
                    serves_d=serves_mask(partition.sf_heads, n_heads=n_heads, device=device)),
                actuator_view=nmii_actuator.geometry, port=sf_port,
                r0_bind_d=wp.zeros(n_heads, dtype=wp.float64, device=device),
                r0_bind_snap_d=wp.zeros(n_heads, dtype=wp.float64, device=device),
                detach_kinetics=sf_motor_params.detach_kinetics, device=device,
                name=NMII_SF_MOTOR, component_b=SF_COMPONENT,
            )
        sf_motor_report = {
            "bound": sf_motor_connector is not None,
            "sf_heads": int(partition.sf_heads.size),
            "cortex_heads": int(partition.cortex_heads.size),
            "sf_head_fraction": partition.sf_fraction,
            "tie_broken_to_cortex": int(partition.tie_broken_to_cortex),
        }

    # `focal_adhesion` — the component, from the SF population's OWN declared FA sites. The count is
    # DERIVED (one clutch slot per `fa_sites` entry the SF builder placed), never typed here, and the
    # population is fixed at build because a slot is the EXISTENCE of a potential adhesion — its
    # bound/unbound history is the state, not the allocation. Deallocating on unbind would make the
    # population non-conserved and break the disjoint-ID re-assert.
    fa_owner = None
    fa_sites = np.asarray(getattr(sf_population, "fa_sites", np.zeros(0, np.int64))).reshape(-1)
    if fa_sites.size:
        from aleph.engine.focal_adhesion import build_focal_adhesion_state_owner

        fa_owner = build_focal_adhesion_state_owner(
            cluster_ids=np.arange(fa_sites.size, dtype=np.int64), device=device,
        )

    # `sf_cytosol_transfer` and `nmii_cytosol_transfer` — two more of the SEVEN edges in the
    # `immersed_transfer` family, sharing the runtime the cortex edge already uses. Nothing in
    # `-alpha*V*grad(p)` is cortical (that module says so), so what was missing per edge was only its
    # component's immersed control volumes.
    #
    # The volume convention is READ FROM THE CORTEX EDGE, not invented here: a filament-network node is
    # given the uniform immersed cell `dx**3` of the fluid grid it sits in. sf_arc and the NMII backbone
    # are filament networks in that same grid, so they take the same convention. A different rule for
    # these two would make the three edges physically inconsistent while all three still "bound".
    solid_cytosol_transfers: dict[str, Any] = {}
    if interior_column is not None and cytosol_grid_dx_um is not None:
        from aleph.engine.immersed_transfer import build_immersed_porous_transfer

        cell_volume_um3 = float(cytosol_grid_dx_um) ** 3
        for edge, owner in ((SF_CYTOSOL_TRANSFER, sf_owner), (NMII_CYTOSOL_TRANSFER, nmii_actuator)):
            positions = getattr(owner, "position_d", None)
            if positions is None:
                continue
            solid_cytosol_transfers[edge] = build_immersed_porous_transfer(
                name=edge,
                pressure_coupling=interior_column.cortex_transfer._pressure_coupling,
                node_volume_d=wp.full(int(positions.shape[0]), cell_volume_um3,
                                      dtype=wp.float64, device=device),
                architecture=architecture,
            )

    composed = compose_native_cell_world(
        sf_owner=sf_owner,
        fa_owner=fa_owner,
        solid_cytosol_transfers=solid_cytosol_transfers,
        ecm_owner=ecm_owner,
        nmii_actuator=nmii_actuator,
        nmii_cortex_connector=nmii_cortex_connector,
        cortex_port=cortex_port,
        clock=clock,
        population_ledgers=population_ledgers,
        # Built above by `build_sf_arc_state_owner` (with_internal_arc defaults True) and, until
        # 2026-08-11, never handed on — so the composed world's actor did not know the edge existed.
        sf_internal_arc_connector=sf_owner.arc_connector,
        sf_cortex_connector=sf_cortex_connector,
        sf_motor_connector=sf_motor_connector,
        interior_column=interior_column,
        interior_solve=interior_solve,
        architecture=architecture,
    )
    composed.sf_motor_report = sf_motor_report
    return composed


class _ECMOwnerTransaction:
    """Accepted-step transaction delegate for the composed ECM owner (mechanism slice — no remesh commit).

    The composed native gate exercises the ECM owner's force pass + the accepted-step commit/refresh path;
    damage/sleep/refinement/remesh commits are a later ECM biology phase, so each hook is a device-free no-op
    that never advances a private clock or mutates authoritative state.  The owner's own
    ``commit_irreversible`` still re-derives ``k_seg``/``α`` after this returns (a value-idempotent no-op when
    the committed topology did not change).
    """

    def snapshot_candidate(self) -> None:
        return None

    def rollback(self, accepted: wp.array) -> None:
        return None

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        return None


class _ECMOwnerLedger:
    """Ledger delegate for the composed ECM owner (no reduction in this mechanism slice)."""

    def accumulate_ledger(self, ledger: object) -> None:
        return None
