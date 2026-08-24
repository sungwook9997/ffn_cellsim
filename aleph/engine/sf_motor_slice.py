r"""``nmii_sf_motor``: the stress fiber's ACTIVE contractile tension as an emergent myosin-event outcome.

WHAT THIS CLOSES.  ``sf_arc`` reached ``KERNEL_BOUND`` for its PASSIVE rod-cable laws and
:mod:`aleph.engine.sf_mechanics` recorded the exact remaining hole in
:data:`~aleph.engine.sf_mechanics.SF_CONNECTOR_BINDING_STATUS`: the active prestress *"is the reaction of
the discrete NMII heads, which enters through a separate MOTOR connector (the SF counterpart of
``nmii_cortex_motor``) … left honestly SEAMED because it needs the NMII actuator + a live SF bind-target port."*
``nmii_sf_motor`` has been a declared edge of the reference graph since the architecture landed
(:func:`~aleph.engine.nmii_actuator.nmii_contract_proposal`, MOTOR, ``nmii`` ⟷ ``sf_arc``, bidirectional +
adjoint) with no runtime behind it.  This module supplies that runtime.

The slice drives :class:`~aleph.engine.transaction.CellTransaction` over::

    sf_arc F-actin  ──(nmii_sf_motor : MOTOR connector)──  NMII actuator (explicit heads on the sarcomeres)

and it is the SF analogue of :mod:`aleph.engine.cortex_motor_slice`, reusing that lane's
:class:`~aleph.engine.cortex_motor_slice.FilamentMotorConnector` verbatim (target-agnostic since
2026-07-25) rather than re-deriving any motor physics.  No lumped ``k_SF`` bundle spring exists anywhere: the
SF axial tension EMERGES from whichever heads the accepted-step KMC has bound.

TWO OWNERSHIP DIFFERENCES FROM THE CORTEX SLICE (both make this the stricter case):

1. **``sf_arc`` is a real transaction PARTICIPANT, not a bind-target port only.**  The cortex slice treats the
   cortex purely as a port (the incumbent driver relaxes it), so it has exactly two participants.  Here the SF
   component OWNS its arrays and LAUNCHES its own passive mechanics, so the slice has THREE participants —
   ``sf_arc`` (:class:`~aleph.engine.composed_native.SFArcStateOwner`), ``nmii``, and the connector — each
   snapshotting / rolling back / committing its own state under one predicate.
2. **The split ownership is physical, not logical.**  The cortex lane has the actuator and the port ALIAS one
   global ``cell.pos_d`` (honest, but a single array with global indices).  Here the minifilament particles live
   in a separate engine-owned array from the SF nodes, so the crossbridge's two-array adjoint scatter
   (``+f`` on the head, ``−(1−t)f`` / ``−t·f`` on the two SF nodes) is a genuine inter-component transfer —
   Newton's 3rd law across two never-merged arrays, which is what the connector contract asks for.

WHY THE HEADS ACTUALLY REACH THE ACTIN.  A bipolar minifilament can only straddle an ANTI-PARALLEL PAIR.  The SF
sarcomere is built as that pair with the lateral separation and axial overlap DERIVED from the minifilament
(``lateral = 2·head_offset``, ``overlap = backbone_length``;
:func:`~aleph.engine.sf_population.build_sf_arc_population`), and
:func:`~aleph.engine.sf_nmii_population.build_sf_straddle_nmii_population` places the motors on the known
pairs.  Measured on the host for the default population: every head's perpendicular distance to the filament
line it sits on is ≈0 µm (median 0.000, max 0.055) versus 0.200 µm for the legacy co-located sarcomere — i.e.
binding no longer sits at the very edge of the capture reach, and the ``+``/``−`` head groups land on the two
DIFFERENT filaments (bipolar dot ≈ −0.92 ⇒ the pair contracts inward).

PARAMS ARE PI-GAPs (report-not-tune).  ``k_on``, ``f_stall``, ``v0``, ``kappa``, ``k_xb``, the catch-slip
constants, the SF axial ``k_axial`` and every NMII mechanical stiffness are unresolved KB/PI GAPs (cards N1–N9,
S1).  Both config paths REQUIRE them (:class:`~aleph.engine.cortex_motor_slice.NMIIMotorParams` and
:func:`build_sf_nmii_actuator_state`) and never default: this slice demonstrates the MECHANISM (tension rises
with the bound-head population, from zero), not a quantitative magnitude.

engine units: length µm, force pN, stiffness pN/µm.  Runtime: NVIDIA Warp on CUDA only (I0-A) — the device
builders below run on the gbook A5000; the module itself is CPU-importable so the structural gates drive it with
recording doubles.

Sanity Gate (self-tested in tests/ac/engine/test_sf_motor_slice.py):
  * ownership: the three participants are the ``sf_arc`` state owner, the ``nmii`` state owner, and the
    ``nmii_sf_motor`` connector; the connector is bound at the ``nmii_sf_motor`` graph slot with
    ``component_b == "sf_arc"``; a cortex port handed to it raises (no silent mis-binding).
  * no-double-count: the ``sf_arc`` / ``nmii`` population ledgers are disjoint at build and re-asserted after
    every accepted step; the minifilament particle array is NOT the SF node array (distinct storage).
  * boundary/sign: the split crossbridge is momentum-conserving (``f_head + f_seg_a + f_seg_b = 0``); the KMC
    and the clock advance only under the device-resident accepted predicate (no host ``bool()``).
  * acceptance: that momentum conservation is not merely asserted — each body reduces its OWN force array
    into one half of :class:`~aleph.engine.ledger.GlobalCellLedger`'s balance gate, whose device flag
    IS the acceptance predicate when ``accepted_d`` is omitted, with the tolerance derived on-device from
    the ledger's own resultants (PI D8).  The gate tests the adjoint wiring, not the convergence, and the
    driver's positive control shows it rejecting a deliberately unbalanced step.
  * conservation: a rejected step restores the binding SoA and the SF/NMII positions bit-exactly.
  * REQUIRED-PARAM: every kinetic/mechanical GAP must be supplied; the builders raise otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Any, Callable

import numpy as np
import warp as wp

from aleph.engine.actor import CellActor
from aleph.engine.composed_native import SFArcStateOwner, build_sf_arc_state_owner
from aleph.engine.contracts import CellArchitecture, reference_cell_architecture
from aleph.engine.cortex_motor_slice import (
    FilamentMotorConnector,
    NMIIMotorParams,
    SegmentAttachQuery,
    allocate_segment_connector_state,
    build_nmii_actuator_state,
)
from aleph.engine.nmii_actuator import (
    NMII_COMPONENT,
    NMII_SF_MOTOR,
    SF_COMPONENT,
    FilamentMotorPortView,
    NMIIActuatorStateOwner,
    NMIIBendingStiffness,
    proposed_nmii_architecture,
)
from aleph.engine.population import PopulationLedger
from aleph.engine.protrusion import proposed_protrusion_architecture
from aleph.engine.sf_mechanics import SFMechanicsTopology, build_sf_mechanics_topology
from aleph.engine.sf_nmii_population import SFNMIIPopulation, build_sf_straddle_nmii_population
from aleph.engine.sf_population import SFArcPopulation
from aleph.engine.transaction import CellTransaction
from aleph.engine.world import CellWorldTransaction
from aleph.components.motor.minifilament_topology import MinifilamentTopology
from aleph.components.motor.segment_motor import SegmentDetachKinetics

__all__ = [
    "NMII_SF_MOTOR",
    "SFMotorSlice",
    "build_sf_motor_port",
    "build_sf_motor_slice",
    "build_sf_nmii_actuator_state",
    "sf_sarcomere_geometry_for",
]


def _require(name: str, value: float | None, *, allow_zero: bool = False) -> float:
    """Return a supplied positive-finite (or nonnegative-finite) value, or raise naming the GAP."""
    bad = value is None or not isfinite(value) or (value < 0.0 if allow_zero else value <= 0.0)
    if bad:
        floor = "nonnegative-finite" if allow_zero else "positive-finite"
        raise ValueError(
            f"REQUIRED-PARAM {name!r} must be a supplied {floor} value (KB/PI GAP — source it or surface to "
            f"PI, never default); got {value!r}"
        )
    return float(value)


def sf_sarcomere_geometry_for(topology: MinifilamentTopology) -> dict[str, float]:
    """Return the sarcomere geometry DERIVED from the minifilament that must sit inside it.

    Both values follow from the motor's own rigid geometry, so neither is a free constant (Magic-Number Block):

    * ``sarcomere_lateral_um = 2·head_offset_um`` — the ``+`` heads sit at ``+head_offset`` and the ``−`` heads
      at ``−head_offset`` across the backbone, so the two anti-parallel filaments must be exactly that far
      apart for both head groups to lie ON the actin (perpendicular crossbridge residual → 0).
    * ``sarcomere_overlap_um = backbone_length_um`` — the heads are distributed along the whole backbone
      contour, so the shared span must be at least that long for every head to face actin.

    Feed the result straight into :func:`~aleph.engine.sf_population.build_sf_arc_population`.
    """
    return {
        "sarcomere_lateral_um": 2.0 * float(topology.head_offset_um),
        "sarcomere_overlap_um": float(topology.backbone_length_um),
    }


def build_sf_motor_port(
    *,
    position_d: wp.array,
    force_d: wp.array,
    topology: SFMechanicsTopology,
    device: str,
) -> FilamentMotorPortView:
    """Build the ``sf_arc`` bind-target port over the SF-owned live segment topology (CUDA lane).

    The port exposes the SF component's OWN ``position_d``/``force_d`` (the arrays its passive mechanics and the
    outer relax already integrate) plus the per-segment topology.  Every array comes from the ONE loop that
    emitted the axial links (:func:`build_sf_mechanics_topology`), so segment order, barbed polarity, persistent
    filament identity and material arclength cannot drift apart.  The α-actinin ``arc_joints`` are deliberately
    excluded: they are crosslinks, not filament segments, and a myosin head must not bind them.

    What the motor kernels ACTUALLY read is ``position_d``, ``force_d``, ``segment_node_a_d``/``_b_d`` and
    ``segment_polarity_d``.  ``persistent_filament_id_d``, ``material_s0_d``/``s1_d`` and ``topology_epoch_d``
    are identity/remap metadata that NO kernel consumes today (the same is true of the landed cortex port); they
    are carried because they are the handles a future SF remesh must remap stale ``seg_id`` bindings against, and
    because the filament id makes the disjoint-population invariant checkable from a dump.  Filling them with
    real values rather than placeholders costs nothing here since the emitting loop already has them.

    Polarity convention (``refresh_segment_barbed_kernel``): links are emitted node-index-increasing, so a
    segment inherits its filament's barbed flag unchanged — ``+1`` ⇒ barbed end on the ``node_b`` side, ``−1``
    ⇒ on the ``node_a`` side.  A sarcomere's two filaments carry OPPOSITE flags (barbed ends outward), which is
    exactly what makes the straddling minifilament contractile.

    Args:
        position_d: the ``sf_arc``-owned node positions (``wp.vec3d``, CUDA).
        force_d: the ``sf_arc``-owned node force array the crossbridge reaction scatters into.
        topology: the built SF mechanics topology (host) carrying the per-segment metadata.
        device: the resolved CUDA device string (never a hard-coded id).

    Returns:
        The validated :class:`FilamentMotorPortView` for ``sf_arc``.

    Raises:
        ValueError: if the topology carries no axial segment or its per-segment metadata is missing/mismatched
            (an older topology built before the metadata was emitted).
    """
    n_seg = int(topology.n_links)
    if n_seg <= 0:
        raise ValueError("sf_arc motor port needs at least one axial segment")
    for label, array in (("seg_filament_id", topology.seg_filament_id),
                         ("seg_polarity", topology.seg_polarity),
                         ("seg_s0", topology.seg_s0), ("seg_s1", topology.seg_s1)):
        if int(np.asarray(array).size) != n_seg:
            raise ValueError(
                f"SF topology {label} has {int(np.asarray(array).size)} rows but {n_seg} axial segments — "
                f"rebuild it with build_sf_mechanics_topology so the per-segment metadata is emitted in link "
                f"order (a re-walk of the population would risk drifting out of that order)"
            )
    links = np.ascontiguousarray(topology.links.reshape(-1, 2))
    return FilamentMotorPortView(
        component=SF_COMPONENT,
        position_d=position_d,
        force_d=force_d,
        segment_node_a_d=wp.array(np.ascontiguousarray(links[:, 0], dtype=np.int32),
                                  dtype=wp.int32, device=device),
        segment_node_b_d=wp.array(np.ascontiguousarray(links[:, 1], dtype=np.int32),
                                  dtype=wp.int32, device=device),
        segment_polarity_d=wp.array(np.ascontiguousarray(topology.seg_polarity, dtype=np.int32),
                                    dtype=wp.int32, device=device),
        persistent_filament_id_d=wp.array(np.ascontiguousarray(topology.seg_filament_id, dtype=np.int64),
                                          dtype=wp.int64, device=device),
        material_s0_d=wp.array(np.ascontiguousarray(topology.seg_s0, dtype=np.float64),
                               dtype=wp.float64, device=device),
        material_s1_d=wp.array(np.ascontiguousarray(topology.seg_s1, dtype=np.float64),
                               dtype=wp.float64, device=device),
        topology_epoch_d=wp.zeros(1, dtype=wp.int32, device=device),
    )


def build_sf_nmii_actuator_state(
    population: SFNMIIPopulation,
    *,
    k_backbone_pn_per_um: float,
    k_head_spring_pn_per_um: float,
    k_xb_pn_per_um: float,
    backbone_persistence_length_um: float,
    device: str,
) -> NMIIActuatorStateOwner:
    """Upload a straddle-placed SF minifilament population into an ``nmii`` state owner (CUDA lane).

    The particle array allocated here is the ``nmii`` component's OWN storage — it is NOT the ``sf_arc`` node
    array, so the crossbridge connector performs a genuine two-array adjoint transfer.  The F6 angle-harmonic
    stiffnesses are DERIVED (never fit) by :class:`NMIIBendingStiffness` from the supplied persistence length
    and crossbridge stiffness; the axial/arm stiffnesses are passed through as REQUIRED params.

    Args:
        population: the host straddle-placed minifilament population.
        k_backbone_pn_per_um: backbone rod stiffness [pN/µm] (REQUIRED — GAP).
        k_head_spring_pn_per_um: head↔backbone arm stiffness [pN/µm] (REQUIRED — GAP).
        k_xb_pn_per_um: crossbridge stiffness [pN/µm] (REQUIRED — GAP); feeds the derived arm-orientation
            stiffness so the F6 lever is consistent with the crossbridge it transmits.
        backbone_persistence_length_um: minifilament backbone ``L_p`` [µm] (REQUIRED — GAP; a mature-bipolar
            value is unresolved, so the caller must say what it supplied and why).
        device: the resolved CUDA device string.

    Returns:
        The validated :class:`NMIIActuatorStateOwner` for the ``nmii`` component.
    """
    k_backbone = _require("k_backbone_pn_per_um", k_backbone_pn_per_um)
    k_head_spring = _require("k_head_spring_pn_per_um", k_head_spring_pn_per_um)
    k_xb = _require("k_xb_pn_per_um", k_xb_pn_per_um)
    lp = _require("backbone_persistence_length_um", backbone_persistence_length_um)

    topo = population.topology
    r0_backbone = float(topo.segment_length_um)
    r0_head = float(topo.head_offset_um)
    n_mf = int(population.n_minifilaments)
    n_heads = int(population.n_heads)
    n_particles = int(population.n_particles)
    if int(population.backbone_angles.shape[0]) == 0:
        raise ValueError(
            "the minifilament topology produced no backbone-angle triple (n_bb < 3), which the F6 "
            "BackboneArmMechanics requires; supply a topology with at least three backbone beads"
        )

    def up(array: np.ndarray, dtype: object) -> wp.array:
        return wp.array(np.ascontiguousarray(array), dtype=dtype, device=device)

    bending = NMIIBendingStiffness(
        persistence_length_um=lp, segment_len_um=r0_backbone,
        k_xb_pn_per_um=k_xb, r0_head_um=r0_head,
    )
    return build_nmii_actuator_state(
        n_minifilaments=n_mf, n_heads=n_heads,
        position_d=up(population.pos.astype(np.float64), wp.vec3d),
        minifilament_offset_d=up(population.minifilament_offset.astype(np.int32), wp.int32),
        active_minifilament_d=up(np.ones(n_mf, np.int32), wp.int32),
        minifilament_id_d=up(np.asarray(population.ledger.block, np.int64)[0]
                             + np.arange(n_mf, dtype=np.int64), wp.int64),
        particle_role_d=up(np.zeros(n_particles, np.int32), wp.int32),
        head_node_d=up(population.head_node.astype(np.int32), wp.int32),
        head_id_d=up(np.arange(n_heads, dtype=np.int64), wp.int64),
        head_side_d=up(population.head_side.astype(np.int32), wp.int32),
        backbone_bonds_d=up(population.backbone_bonds.reshape(-1, 2).astype(np.int32), wp.int32),
        head_bonds_d=up(population.head_bonds.reshape(-1, 2).astype(np.int32), wp.int32),
        backbone_angles_d=up(population.backbone_angles.reshape(-1, 3).astype(np.int32), wp.int32),
        head_arm_angles_d=up(population.head_arm_angles.reshape(-1, 3).astype(np.int32), wp.int32),
        k_backbone_pn_per_um=k_backbone, r0_backbone_um=r0_backbone,
        k_head_spring_pn_per_um=k_head_spring, r0_head_um=r0_head,
        bending=bending, device=device,
    )


@dataclass(slots=True)
class SFMotorSlice:
    """The wired ``sf_arc`` ⟷ ``nmii`` dynamic slice: THREE participants under one accepted-step transaction.

    Holds the SF state owner (which launches its own passive rod-cable mechanics), the NMII state owner, the
    ``nmii_sf_motor`` connector, the SF bind-target port, the device clock, and the one
    :class:`CellTransaction`.  Constructed from already-assembled parts so it is CPU-importable and the
    structural gates build it with recording doubles; :func:`build_sf_motor_slice` assembles the real
    Warp-resident parts on the CUDA lane.
    """

    transaction: CellTransaction
    sf_owner: SFArcStateOwner
    actuator_state: NMIIActuatorStateOwner
    connector: FilamentMotorConnector
    port: FilamentMotorPortView
    clock: Any
    ledger: Any | None = None
    tol_sq_d: wp.array | None = None
    balance_gamma_n_d: wp.array | None = None
    _launch: Callable[..., object] = wp.launch

    @property
    def participants(self) -> tuple[object, ...]:
        """The deduplicated accepted-step participants (SF owner + NMII owner + the motor connector)."""
        return self.transaction.world.participants

    def zero_forces(self) -> None:
        """Zero both owners' force arrays before a candidate accumulate (the caller owns force zeroing)."""
        from aleph.engine.composed_native import _zero_vec3d_kernel

        for force in (self.sf_owner.force_d, self.actuator_state.force_d):
            self._launch(_zero_vec3d_kernel, dim=int(force.shape[0]), inputs=[force],
                         device=str(force.device))

    def accumulate(self) -> None:
        """One candidate-force pass: SF passive mechanics + NMII internal mechanics + the split crossbridge.

        Adds (1) the ``sf_arc`` axial-backbone/bending (and internal dorsal↔arc joint) forces into the SF-owned
        force array, (2) the minifilament's own backbone/head-arm mechanics into the ``nmii``-owned array, and
        (3) the split crossbridge — ``+f`` on the head, ``−(1−t)f``/``−t·f`` on the two SF nodes (Newton's 3rd
        law across two never-merged arrays).  No fixed prestress is injected anywhere: the SF axial tension is
        the reaction of whatever heads the KMC has bound so far.
        """
        self.sf_owner.accumulate()
        self.actuator_state.accumulate_internal()
        self.connector.accumulate_candidate(self.actuator_state.geometry(), self.port)

    def step(
        self,
        solve: Callable[[], None],
        *,
        dt_phys: float,
        accepted_d: wp.array | None = None,
        ledger: Any | None = None,
        tol_sq_d: wp.array | None = None,
        balance_gamma_n_d: wp.array | None = None,
    ) -> None:
        """Run exactly one accepted physical step of the dynamic SF-motor runtime.

        Args:
            solve: the caller-owned inner mechanical solve (zero forces → :meth:`accumulate` → integrate the SF
                nodes and the minifilament particles toward convergence).  Runs after events are proposed and
                owns all physics.
            dt_phys: the outer physical timestep [s]; finite and positive.
            accepted_d: the device acceptance predicate.  If ``None`` the balance gate decides — then a
                ``ledger`` with a ``tol_sq_d`` must be available.  Passing ``None`` is what makes this lane
                accept on PHYSICS (the two bodies' force resultants cancelling across the connector) rather
                than on a constant array of ones.
            ledger / tol_sq_d / balance_gamma_n_d: optional balance-gate inputs; default to the slice's
                stored ones.  ``balance_gamma_n_d`` carries the D8 float64 accumulation ratio, from which
                the transaction derives ``tol_sq_d`` on the device after this step's own contributions land.

        The per-head Hill/Bell loads are computed on the CONVERGED geometry (after ``solve``, before the KMC
        commit consumes them).  The predicate stays device-resident throughout.

        The ledger's accumulators are ZEROED here before the pass: the balance channels are per-step
        resultants, so carrying the previous step's sums into this step's gate would compare a residual
        against the wrong configuration's scale.
        """
        ledger = ledger if ledger is not None else self.ledger
        tol_sq_d = tol_sq_d if tol_sq_d is not None else self.tol_sq_d
        balance_gamma_n_d = (
            balance_gamma_n_d if balance_gamma_n_d is not None else self.balance_gamma_n_d)

        def _inner() -> None:
            solve()
            self.connector.compute_loads(self.actuator_state.geometry(), self.port)

        if ledger is not None:
            ledger.reset()
        self.transaction.step(
            dt_phys=dt_phys, solve=_inner, accepted_d=accepted_d, ledger=ledger, tol_sq_d=tol_sq_d,
            balance_gamma_n_d=balance_gamma_n_d,
        )

    def bound_head_count(self) -> int:
        """Host read of the bound-head population for OUT-OF-LOOP diagnostics only (device→host readback).

        A viz/telemetry helper for BETWEEN steps; never call it inside the inner solve (I0-A forbids an
        authoritative GPU→CPU roundtrip in the physical-time loop).
        """
        return int(self.connector.state.bound_d.numpy().sum())


def build_sf_motor_slice(
    *,
    sf_population: SFArcPopulation,
    minifilament_topology: MinifilamentTopology,
    params: NMIIMotorParams,
    sf_k_axial_pn_per_um: float,
    nmii_k_backbone_pn_per_um: float,
    nmii_k_head_spring_pn_per_um: float,
    nmii_backbone_persistence_length_um: float,
    base_seed: int,
    device: str,
    nmii_id_base: int = 3_000_000,
    extra_population_ledgers: tuple[PopulationLedger, ...] = (),
    clock: Any | None = None,
    ledger: Any | None = None,
    tol_sq_d: wp.array | None = None,
    balance_gamma_n_d: wp.array | None = None,
    architecture: CellArchitecture | None = None,
) -> SFMotorSlice:
    """Assemble the CUDA-lane dynamic SF-motor slice: three participants under one accepted-step clock.

    Builds, in order: the SF mechanics topology + the real ``sf_arc`` state owner (its own device arrays and
    kernels); the SF bind-target port over those arrays; the straddle-placed NMII minifilament population and
    its ``nmii`` state owner (a SEPARATE particle array); the ``nmii_sf_motor``
    :class:`FilamentMotorConnector` with its attach-query accelerator; then one :class:`CellActor` binding
    ``sf_arc`` + ``nmii`` + ``nmii_sf_motor``, one :class:`CellWorldTransaction` (``require_complete=False`` —
    this is a vertical slice, not the whole cell) and one :class:`CellTransaction` carrying the disjoint
    population ledgers.

    Args:
        sf_population: a built ``sf_arc`` population WITH the straddle sarcomere geometry (see
            :func:`sf_sarcomere_geometry_for`); a legacy co-located population raises in the placement builder.
        minifilament_topology: the NMII reference topology the sarcomere geometry was derived from.
        params: the sourced-or-GAP head kinetic/crossbridge constants (all REQUIRED).
        sf_k_axial_pn_per_um: SF actin axial backbone stiffness [pN/µm] (REQUIRED — modelling GAP).
        nmii_k_backbone_pn_per_um / nmii_k_head_spring_pn_per_um: minifilament internal stiffnesses (REQUIRED).
        nmii_backbone_persistence_length_um: minifilament backbone ``L_p`` [µm] (REQUIRED — GAP).
        base_seed: the run's fixed host RNG seed (nonnegative int).
        device: the resolved CUDA device string (never a hard-coded id).
        nmii_id_base: inclusive lower bound of the ``nmii`` population ID block (disjoint from ``sf_arc``).
        extra_population_ledgers: further ledgers (e.g. cortex, ECM) to assert disjoint every accepted step.
        clock / ledger / tol_sq_d / balance_gamma_n_d: optional; a :class:`WarpEventClock` is built from
            ``base_seed`` if omitted.  Supplying a ledger + ``tol_sq_d`` + ``balance_gamma_n_d`` is what lets
            :meth:`SFMotorSlice.step` be called with ``accepted_d=None``, i.e. with the two bodies' force
            resultants deciding acceptance on the device instead of a constant predicate.
        architecture: optional composition; defaults to the protrusion+NMII-promoted reference architecture.

    Returns:
        A :class:`SFMotorSlice` whose :meth:`SFMotorSlice.step` runs one accepted physical step.
    """
    from aleph.engine.events import make_event_clock

    if architecture is None:
        architecture = proposed_nmii_architecture(proposed_protrusion_architecture(reference_cell_architecture()))

    # 1. the REAL sf_arc owner (its own device arrays; launches link_spring + cytosim_bending itself).
    topology = build_sf_mechanics_topology(
        sf_population, k_axial_pn_per_um=_require("sf_k_axial_pn_per_um", sf_k_axial_pn_per_um))
    sf_owner = build_sf_arc_state_owner(
        sf_population, k_axial_pn_per_um=float(topology.k_axial_pn_per_um), device=device,
        topology=topology)   # ONE topology feeds both the mechanics and the port (single segment-order SoT)

    # 2. the SF bind-target port over those same SF-owned arrays (the target side of the MOTOR edge).
    port = build_sf_motor_port(
        position_d=sf_owner.position_d, force_d=sf_owner.force_d, topology=topology, device=device)

    # 3. the straddle-placed NMII population -> the nmii state owner (a SEPARATE particle array).
    nmii_population = build_sf_straddle_nmii_population(
        sf_population, minifilament_topology, id_base=int(nmii_id_base))
    actuator_state = build_sf_nmii_actuator_state(
        nmii_population,
        k_backbone_pn_per_um=nmii_k_backbone_pn_per_um,
        k_head_spring_pn_per_um=nmii_k_head_spring_pn_per_um,
        k_xb_pn_per_um=params.k_xb,
        backbone_persistence_length_um=nmii_backbone_persistence_length_um,
        device=device,
    )

    # 4. the nmii_sf_motor connector (target-agnostic FilamentMotorConnector bound to the SF port).
    n_heads = int(actuator_state.n_heads)
    connector_state = allocate_segment_connector_state(n_heads, int(topology.n_links), device=device)
    max_segment_um = float(np.max(topology.link_r0)) if topology.n_links else 0.0
    query = SegmentAttachQuery(
        state=connector_state, head_node_d=actuator_state.head_node_d, port=port,
        capture_radius=params.capture_radius, max_segment_length_um=max_segment_um, device=device,
    )
    catch_slip = (
        params.to_catch_slip_params()
        if params.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP else None
    )
    connector = FilamentMotorConnector(
        state=connector_state, head_node_d=actuator_state.head_node_d, params=params.to_hand_params(),
        catch_slip=catch_slip, query=query, actuator_view=actuator_state.geometry, port=port,
        r0_bind_d=wp.zeros(n_heads, dtype=wp.float64, device=device),
        r0_bind_snap_d=wp.zeros(n_heads, dtype=wp.float64, device=device),
        detach_kinetics=params.detach_kinetics,
        name=NMII_SF_MOTOR, component_b=SF_COMPONENT, device=device,
    )

    # 5. ONE actor + ONE world + ONE transaction over the three participants.
    actor = CellActor(architecture)
    actor.bind_component(SF_COMPONENT, sf_owner)
    actor.bind_component(NMII_COMPONENT, actuator_state)
    actor.bind_connector(NMII_SF_MOTOR, connector)
    world = CellWorldTransaction(actor, require_complete=False)
    if clock is None:
        clock = make_event_clock(base_seed=int(base_seed), device=device)
    transaction = CellTransaction(
        world, clock,
        population_ledgers=(sf_population.ledger, nmii_population.ledger, *extra_population_ledgers),
    )
    return SFMotorSlice(
        transaction=transaction, sf_owner=sf_owner, actuator_state=actuator_state, connector=connector,
        port=port, clock=clock, ledger=ledger, tol_sq_d=tol_sq_d,
        balance_gamma_n_d=balance_gamma_n_d,
    )
