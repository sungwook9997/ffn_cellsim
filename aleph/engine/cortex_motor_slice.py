r"""GATE-B first dynamic-runtime vertical slice: emergent cortical NMII tension.

GATE A closed the STATIC converged resting baseline on the clean turgor-pressurised cortex.  Its verdict:
the resting myosin cortical tension is **not** a static prestress — it is a DYNAMIC steady state that must
EMERGE from explicit NMII head binding/unbinding EVENTS driven by the accepted-step transaction.  This module
is the smallest real slice that shows exactly that.

It drives :class:`~aleph.engine.transaction.CellTransaction` over one vertical column of the composed
cell graph::

    cortex actin  ──(nmii_cortex_motor : MOTOR connector)──  NMII actuator (explicit heads)

Per accepted physical step the transaction:

1. **snapshots** the NMII component state and the connector's per-head binding SoA (one candidate);
2. **proposes attach events** — :class:`CortexMotorConnector.propose_events` runs the injected point-to-segment
   attach-query accelerator over the live cortex port and re-binds the live segment topology (INTEGRATION
   items 2 & 3), writing only candidate buffers (no force, no host state);
3. runs the caller-owned **inner mechanical solve** — during which :meth:`CortexMotorSlice.accumulate` scatters
   the split crossbridge force ``+f`` into the ``nmii``-owned head array and ``−(1−t)f``/``−t·f`` into the two
   cortex port nodes (Newton's 3rd law across two never-merged arrays) plus the minifilament's own internal
   backbone/head-arm mechanics — then the caller relaxes the cortex to convergence;
4. computes the per-head Hill (tangential) / Bell (full-|F|) loads on the converged geometry;
5. **commits** the accepted-predicated KMC: a free head ATTACHES by the ``k_on`` Poisson rate, a bound head
   Hill-steps and DETACHES by the PHYSIOLOGICAL NMII catch-slip hazard (Kovacs 2007, default) — or rolls the
   binding SoA back bit-exactly on a rejected step;
6. **advances the device event clock** iff the step was accepted.

No fixed prestress is added anywhere; the cortical tension EMERGES from the bound-head population the KMC grows
over steps.  Every CUDA op flows through an injected launcher/copy/query so the module is CPU-importable and the
structural gates drive it with recording doubles; the CUDA-lane :func:`build_cortex_motor_slice` assembles the
real Warp-resident runtimes on the gbook A5000.

Sanity Gate:
    * ownership: cortex is a bind-target PORT (never a participant); the two transaction participants are the
      NMII state-owner and the ``nmii_cortex_motor`` connector; co-location is never a connection.
    * boundary/sign: the split crossbridge is momentum-conserving (``f_head + f_seg_a + f_seg_b = 0``); the
      KMC and clock advance are device-gated on the accepted predicate — a rejected step advances no binding
      state and no clock (no host ``bool()``/``numpy()`` on the predicate).
    * conservation: the connector's binding SoA is snapshot/rolled-back bit-exactly; the emergent tension is
      the reaction of the bound-head crossbridges, never an imposed setpoint.
    * REQUIRED-PARAM: ``k_on`` (attach) and the catch-slip constants are PI/KB GAPs — supplied, never invented;
      :class:`CortexMotorParams` refuses to construct without them.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Any, Callable, Protocol, runtime_checkable

import warp as wp

from aleph.engine.actor import CellActor
from aleph.engine.contracts import CellArchitecture, reference_cell_architecture
from aleph.engine.nmii_actuator import (
    BackboneArmMechanics,
    CORTEX_COMPONENT,
    FilamentMotorPortView,
    NMII_COMPONENT,
    NMII_CONNECTOR_EVENT_CHANNELS,
    NMII_CORTEX_MOTOR,
    NMII_LEDGER_CHANNELS,
    NMII_STATE_EVENT_CHANNELS,
    NMIIActuatorStateOwner,
    NMIIActuatorView,
    NMIIBendingStiffness,
    SegmentConnectorState,
    SegmentMotorConnectorRuntime,
    proposed_nmii_architecture,
)
from aleph.engine.transaction import CellTransaction
from aleph.engine.world import CellWorldTransaction
from aleph.components.motor.bell_kinetics_analytic import KBT_PN_UM
from aleph.components.motor.hand import NMIICatchSlipParams, NMIIHandParams
from aleph.components.motor.segment_motor import (
    SegmentDetachKinetics,
    _increment_epoch_if_accepted_kernel,
    _refresh_bound_walk_dir_kernel,
    attach_segment_gated_r0bind_kernel,
    compute_head_loads_segment_split_r0bind_kernel,
    crossbridge_segment_split_r0bind_kernel,
    refresh_segment_barbed_kernel,
    step_detach_segment_gated_catch_slip_kernel,
    step_detach_segment_gated_kernel,
)
from aleph.components.motor.segment_query import SegmentQuery, _gather_head_pos_kernel

__all__ = [
    "AttachQueryAccelerator",
    "CortexMotorConnector",
    "CortexMotorParams",
    "CortexMotorSlice",
    "FilamentMotorConnector",
    "NMIIMotorParams",
    "SegmentAttachQuery",
    "allocate_segment_connector_state",
    "build_cortex_motor_slice",
    "build_nmii_actuator_state",
]

# ── Injected CUDA seam: one launch/copy indirection keeps the wiring exercisable on a CUDA-free host. ──
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


# ── Reject-gated restore kernels for the NMII state-owner's reversible mechanical/candidate arrays. ────
# These select the REJECTED branch in-device (``accepted[0]==0``): a rejected outer step restores the
# pre-candidate value, an accepted step keeps the converged candidate.  The connector's binding SoA uses the
# landed unconditional-copy rollback (solve never mutates it — only the gated commit does); the state-owner's
# position IS relaxed by solve, so it needs the gated restore instead.
@wp.kernel
def _restore_vec3_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.vec3d),
    snap: wp.array(dtype=wp.vec3d),
) -> None:
    """Restore ``live <- snap`` iff the outer step was rejected (device-only, no host branch)."""
    if accepted[0] == 0:
        t = wp.tid()
        live[t] = snap[t]


@wp.kernel
def _restore_int32_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.int32),
    snap: wp.array(dtype=wp.int32),
) -> None:
    """Restore an int32 array ``live <- snap`` iff the outer step was rejected."""
    if accepted[0] == 0:
        t = wp.tid()
        live[t] = snap[t]


@wp.kernel
def _restore_int64_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.int64),
    snap: wp.array(dtype=wp.int64),
) -> None:
    """Restore an int64 array ``live <- snap`` iff the outer step was rejected."""
    if accepted[0] == 0:
        t = wp.tid()
        live[t] = snap[t]


# ── REQUIRED-PARAM discipline: GAP inputs are supplied or the config refuses to build. ────────────────
#: Parameters that are a KB/PI GAP in ``ac/motor/params_i0b3.yaml`` (value: null / "GAP — PI"), NEVER
#: defaulted to a convenient number here.  The magnitude ones additionally gate the emergent-tension
#: comparison; the Lead sources them (or surfaces to PI) before the native GATE-B run.
_GAP_PARAMS = frozenset(
    {"k_on", "f_stall", "v0", "kappa", "k_xb", "r0_head", "k_catch0", "x_catch", "k_slip0", "x_slip"}
)


def _require_param(name: str, value: float | None, *, allow_zero: bool = False) -> float:
    """Return a supplied finite param, or raise a REQUIRED-PARAM error naming the GAP.

    Never substitutes a default: an unset (``None``) or non-finite value halts with a message the Lead
    surfaces to PI, per the no-magic-number / physiological-baseline rules.  Most params must be positive;
    ``allow_zero`` permits a physical zero (e.g. ``r0_xb`` — a bound head sits on the actin site, rest ≈ 0).
    """
    gap = " (KB/PI GAP — source it or surface to PI, never default)" if name in _GAP_PARAMS else ""
    bad = value is None or not isfinite(value) or (value < 0.0 if allow_zero else value <= 0.0)
    if bad:
        floor = "nonnegative-finite" if allow_zero else "positive-finite"
        raise ValueError(f"REQUIRED-PARAM {name!r} must be a supplied {floor} value{gap}; got {value!r}")
    return float(value)


@dataclass(frozen=True, slots=True)
class CortexMotorParams:
    """Sourced-or-GAP kinetic + mechanical constants for the cortex NMII motor slice (engine units).

    Every physical value is REQUIRED (no convenient defaults).  The subset in :data:`_GAP_PARAMS` is a
    ``params_i0b3.yaml`` GAP — PI: ``k_on`` (per-head attach rate, ``null``, provisional 50/s unaudited),
    ``f_stall``/``v0``/``kappa``/``k_xb``/``r0_head`` (magnitude GAPs), and the catch-slip constants
    (``k_catch0``/``x_catch``/``k_slip0``/``x_slip``, constrained by Kovacs 5×/12× + NM2B duty 0.2–0.3 but no
    direct single-molecule fit).  ``k_off0``/``f0`` are the provisional-sourced Bell slip inputs (only used in
    ``SLIP`` mode); ``kT`` is the physical 37 °C constant.  The slip path is retained but ``CATCH_SLIP`` is the
    physiological default, so the catch-slip constants are then mandatory.
    """

    # Hill mechanics (magnitudes are PI GAPs; passed, never invented).
    v0: float
    f_stall: float
    kappa: float
    k_xb: float
    r0_head: float
    r0_xb: float
    capture_radius: float
    # attach (GAP — PI).
    k_on: float
    # Bell slip inputs (provisional-sourced); only consumed when detach_kinetics == SLIP.
    k_off0: float
    f0: float
    # Pereverzev catch-slip constants (GAP — PI); mandatory when detach_kinetics == CATCH_SLIP.
    k_catch0: float | None = None
    x_catch: float | None = None
    k_slip0: float | None = None
    x_slip: float | None = None
    kT: float = KBT_PN_UM
    detach_kinetics: SegmentDetachKinetics = SegmentDetachKinetics.CATCH_SLIP

    def __post_init__(self) -> None:
        object.__setattr__(self, "detach_kinetics", SegmentDetachKinetics(self.detach_kinetics))
        for name in ("v0", "f_stall", "kappa", "k_xb", "r0_head", "capture_radius",
                     "k_on", "k_off0", "f0", "kT"):
            _require_param(name, getattr(self, name))
        _require_param("r0_xb", self.r0_xb, allow_zero=True)  # rest ≈ 0: a bound head sits on the actin site.
        if self.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP:
            for name in ("k_catch0", "x_catch", "k_slip0", "x_slip"):
                _require_param(name, getattr(self, name))

    def to_hand_params(self) -> NMIIHandParams:
        """Build the device-resident :class:`NMIIHandParams` (CUDA lane; constructs a Warp struct)."""
        params = NMIIHandParams()
        params.k_on = wp.float64(self.k_on)
        params.k_off0 = wp.float64(self.k_off0)
        params.f0 = wp.float64(self.f0)
        params.v0 = wp.float64(self.v0)
        params.f_stall = wp.float64(self.f_stall)
        params.kappa = wp.float64(self.kappa)
        params.k_xb = wp.float64(self.k_xb)
        params.r0_head = wp.float64(self.r0_head)
        params.r0_xb = wp.float64(self.r0_xb)
        params.capture_radius = wp.float64(self.capture_radius)
        return params

    def to_catch_slip_params(self) -> NMIICatchSlipParams:
        """Build the device-resident :class:`NMIICatchSlipParams` (CUDA lane; ``CATCH_SLIP`` only)."""
        if self.detach_kinetics is not SegmentDetachKinetics.CATCH_SLIP:
            raise ValueError("catch-slip params requested for a non-CATCH_SLIP config")
        cs = NMIICatchSlipParams()
        cs.k_catch0 = wp.float64(_require_param("k_catch0", self.k_catch0))
        cs.x_catch = wp.float64(_require_param("x_catch", self.x_catch))
        cs.k_slip0 = wp.float64(_require_param("k_slip0", self.k_slip0))
        cs.x_slip = wp.float64(_require_param("x_slip", self.x_slip))
        cs.kT = wp.float64(self.kT)
        return cs


# ── NMII component state-owner delegates (concrete, minimal, faithful). ───────────────────────────────
class _ActuatorStateTransaction:
    """Concrete :class:`~aleph.engine.nmii_actuator.NMIITransaction` for the NMII state-owner.

    Snapshots the reversible/candidate arrays the state-owner declares mutable (position, ATP cycle/consumed,
    active population, RNG epoch); on a REJECTED step restores them in-device via the reject-gated kernels; on
    an ACCEPTED step advances the actuator RNG epoch.  Position is the only array the inner solve actually
    relaxes in this slice — the ATP/population arrays are static until the biology phase adds turnover events,
    but they are covered so the transaction contract is honoured.
    """

    component_name = NMII_COMPONENT
    event_channels = NMII_STATE_EVENT_CHANNELS

    def __init__(
        self,
        *,
        position_d: wp.array,
        active_minifilament_d: wp.array,
        atp_cycle_state_d: wp.array,
        atp_consumed_d: wp.array,
        actuator_epoch_d: wp.array,
        position_snap_d: wp.array,
        active_snap_d: wp.array,
        atp_state_snap_d: wp.array,
        atp_consumed_snap_d: wp.array,
        epoch_snap_d: wp.array,
        launch: LaunchFn = _default_launch,
        copy: CopyFn = _default_copy,
        device: object | None = None,
    ) -> None:
        self._vec3_pairs = ((position_d, position_snap_d),)
        self._int32_pairs = (
            (active_minifilament_d, active_snap_d),
            (atp_cycle_state_d, atp_state_snap_d),
            (actuator_epoch_d, epoch_snap_d),
        )
        self._int64_pairs = ((atp_consumed_d, atp_consumed_snap_d),)
        self._epoch_d = actuator_epoch_d
        self._mutable = (
            position_d, active_minifilament_d, atp_cycle_state_d, atp_consumed_d, actuator_epoch_d,
        )
        self._launch = launch
        self._copy = copy
        self._device = device
        self.commit_calls: list[tuple[Any, ...]] = []

    def owned_arrays(self) -> tuple[wp.array, ...]:
        """Return the authoritative mutable arrays covered by snapshot/rollback (storage-key checked)."""
        return self._mutable

    def snapshot_candidate(self) -> None:
        """D2D-snapshot every covered array (``snap <- live``) before candidate mutation."""
        for live, snap in (*self._vec3_pairs, *self._int32_pairs, *self._int64_pairs):
            self._copy(snap, live)

    def rollback(self, accepted: wp.array) -> None:
        """Reject-gated restore of every covered array; the accepted branch keeps the converged candidate."""
        for live, snap in self._vec3_pairs:
            self._launch(_restore_vec3_if_rejected_kernel, dim=int(live.shape[0]),
                         inputs=[accepted, live, snap], device=self._device)
        for live, snap in self._int32_pairs:
            self._launch(_restore_int32_if_rejected_kernel, dim=int(live.shape[0]),
                         inputs=[accepted, live, snap], device=self._device)
        for live, snap in self._int64_pairs:
            self._launch(_restore_int64_if_rejected_kernel, dim=int(live.shape[0]),
                         inputs=[accepted, live, snap], device=self._device)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Advance the actuator RNG epoch iff the step was accepted (device-gated)."""
        self.commit_calls.append((accepted, dt_phys, rng_seed))
        self._launch(_increment_epoch_if_accepted_kernel, dim=1,
                     inputs=[accepted, self._epoch_d], device=self._device)


class _ActuatorStateLedger:
    """Concrete :class:`~aleph.engine.nmii_actuator.NMIILedgerContributor` for the NMII state-owner.

    Declares the four motor accounting channels and records the ledger handle; the force/work/ATP/population
    reductions land with the two-array kernel on the native lane (same seam the connector uses).
    """

    ledger_channels = NMII_LEDGER_CHANNELS

    def __init__(self) -> None:
        self.ledger_calls: list[object] = []

    def accumulate_ledger(self, ledger: object) -> None:
        """Record the ledger handle without deciding acceptance."""
        self.ledger_calls.append(ledger)


# ── Attach-query accelerator seam. ────────────────────────────────────────────────────────────────────
@runtime_checkable
class AttachQueryAccelerator(Protocol):
    """Fills the connector's accepted-step attach scratch from the live target port (INTEGRATION item 2)."""

    def fill_attach_scratch(self, actuator: NMIIActuatorView, port: FilamentMotorPortView) -> None:
        """Write ``query_seg_id_d`` / ``query_t_d`` / ``query_barbed_d`` for every head (candidate only)."""


class SegmentAttachQuery:
    """Production attach-query accelerator over a device HashGrid of the cortex port's segment midpoints.

    Runs :class:`~aleph.components.motor.segment_query.SegmentQuery` over the TARGET port's live segments (grid from
    segment midpoints) and gathers head query points from the ``nmii``-owned actuator array (split ownership —
    the heads are NOT co-located in the port array), writing ``seg_id``/``t``/``barbed`` directly into the
    connector's attach scratch (no extra copy).  CUDA-lane only (it allocates a ``wp.HashGrid``); the module is
    importable on a CUDA-free host and the structural gates inject a recording spy instead.
    """

    def __init__(
        self,
        *,
        state: SegmentConnectorState,
        head_node_d: wp.array,
        port: FilamentMotorPortView,
        capture_radius: float,
        max_segment_length_um: float,
        grid_dim: int = 128,
        device: str | None = None,
    ) -> None:
        if capture_radius <= 0.0 or max_segment_length_um <= 0.0:
            raise ValueError("capture_radius and max_segment_length_um must be positive")
        self._state = state
        self._head_node_d = head_node_d
        self._capture_radius = float(capture_radius)
        # A segment is discoverable through its midpoint whenever the head is within capture of its body.
        self._query_radius = float(capture_radius) + 0.5 * float(max_segment_length_um)
        self._device = device
        self._query = SegmentQuery(
            port.segment_node_a_d, port.segment_node_b_d, state.seg_barbed_d, state.n_heads,
            grid_dim=grid_dim, device=device,
        )
        # Route the query outputs into the connector-owned attach scratch (no copy, INTEGRATION item 2).
        self._query.seg_id = state.query_seg_id_d
        self._query.t = state.query_t_d
        self._query.barbed = state.query_barbed_d
        self._head_pos = wp.zeros(state.n_heads, dtype=wp.vec3d, device=device)

    def fill_attach_scratch(self, actuator: NMIIActuatorView, port: FilamentMotorPortView) -> None:
        """Refresh the live barbed field, gather split-ownership head points, and query the cortex segments."""
        wp.launch(
            refresh_segment_barbed_kernel, dim=int(port.segment_node_a_d.shape[0]),
            inputs=[port.position_d, port.segment_node_a_d, port.segment_node_b_d,
                    port.segment_polarity_d, self._state.seg_barbed_d], device=self._device,
        )
        wp.launch(
            _gather_head_pos_kernel, dim=self._state.n_heads,
            inputs=[self._head_node_d, actuator.position_d], outputs=[self._head_pos], device=self._device,
        )
        self._query.build(port.position_d, self._query_radius)
        self._query.query(self._head_pos, port.position_d, self._query_radius, self._capture_radius)


# ── The graph-owned nmii_cortex_motor connector with the GATE-B dynamic-runtime additions. ────────────
class FilamentMotorConnector(SegmentMotorConnectorRuntime):
    """A graph-owned NMII MOTOR connector against ONE live filament port + the two GATE-B additions.

    Extends :class:`SegmentMotorConnectorRuntime` with the two pieces the dynamic slice needs, without
    changing the landed crossbridge scatter / snapshot / rollback / attach KMC:

    * :meth:`propose_events` — the 5th scheduler method.  Re-binds the live target segment topology used by the
      accepted attach kernel (INTEGRATION item 3) and fills the accepted-step attach scratch via the injected
      attach-query accelerator (INTEGRATION item 2).  Reads geometry only; writes candidate buffers only.
    * :meth:`commit_irreversible` — defaults the DETACH hazard to the PHYSIOLOGICAL NMII catch-slip law
      (Kovacs 2007) instead of the pure Bell slip, so the bound-head duty RISES under load for sustained
      tension.  ``SegmentDetachKinetics.SLIP`` selects the landed base-class behaviour.

    TARGET-AGNOSTIC (2026-07-25).  Every kernel this class launches addresses the target only through the
    :class:`FilamentMotorPortView` interface (``position_d`` / ``force_d`` / ``segment_node_a_d`` /
    ``segment_node_b_d`` / ``segment_polarity_d``), so the same connector serves any declared NMII MOTOR edge —
    ``nmii_cortex_motor`` (the default, kept byte-identical), ``nmii_sf_motor``, and later the protrusion
    edges.  ``name``/``component_b`` are the only cortex-specific values and are now parameters; the base class
    still validates the pair against the architecture's ``_REQUIRED_CONNECTORS``, so a mismatched combination
    raises instead of silently binding the wrong target.  :data:`CortexMotorConnector` remains as the original
    name for the landed cortex call sites.
    """

    def __init__(
        self,
        *,
        state: SegmentConnectorState,
        head_node_d: wp.array,
        params: object,
        catch_slip: NMIICatchSlipParams | None,
        query: AttachQueryAccelerator,
        actuator_view: Callable[[], NMIIActuatorView],
        port: FilamentMotorPortView,
        r0_bind_d: wp.array,
        r0_bind_snap_d: wp.array,
        detach_kinetics: SegmentDetachKinetics = SegmentDetachKinetics.CATCH_SLIP,
        event_channels: frozenset[str] = NMII_CONNECTOR_EVENT_CHANNELS,
        ledger_channels: frozenset[str] = NMII_LEDGER_CHANNELS,
        name: str = NMII_CORTEX_MOTOR,
        component_b: str = CORTEX_COMPONENT,
        launch: LaunchFn | None = None,
        copy: CopyFn | None = None,
        device: object | None = None,
    ) -> None:
        super_kwargs: dict[str, Any] = {}
        if launch is not None:
            super_kwargs["launch"] = launch
        if copy is not None:
            super_kwargs["copy"] = copy
        if port.component != component_b:
            raise ValueError(
                f"MOTOR connector {name!r} targets {component_b!r} but was handed a "
                f"{port.component!r} port (co-location is never a connection — bind the target's own port)"
            )
        super().__init__(
            name=name, component_b=component_b, state=state, head_node_d=head_node_d,
            params=params, event_channels=event_channels, ledger_channels=ledger_channels,
            device=device, **super_kwargs,
        )
        self.detach_kinetics = SegmentDetachKinetics(detach_kinetics)
        if self.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP and catch_slip is None:
            raise ValueError(
                f"REQUIRED-PARAM: {name} CATCH_SLIP detach needs NMIICatchSlipParams "
                "(k_catch0/x_catch/k_slip0/x_slip are GAP — PI; pass detach_kinetics=SLIP for the "
                "provisional-sourced Bell slip smoke path)"
            )
        self._catch_slip = catch_slip
        self._query = query
        self._actuator_view = actuator_view
        self._port = port
        # Per-head zero-strain crossbridge reference (attach-unstrained): captured at the accepted attach and
        # snapshot/rolled-back with the rest of the binding SoA.  Owned here (not in the read-only
        # SegmentConnectorState) so the fix stays in the GATE-B lane; the r0bind kernels consume it.
        self._r0_bind_d = r0_bind_d
        self._r0_bind_snap_d = r0_bind_snap_d
        # INTEGRATION item 3: bind the live target segment endpoints up front (refreshed each propose_events).
        self.bind_target_topology(port.segment_node_a_d, port.segment_node_b_d)

    def binding_field_arrays(self) -> tuple[wp.array, ...]:
        """Return every live array the crossbridge force reads — what an OBSERVER must freeze and restore.

        A path integral is a property of the force FIELD, so it has to be walked with the discrete state
        held fixed, and an observer that alters that state (detaching every head to isolate the passive
        channel, or reinstating a previous step's binding to reproduce the field the step actually used)
        must put back exactly what it changed.  Enumerating the set here rather than in each driver stops
        the two from drifting: a driver restoring one array fewer than the force reads would silently
        measure a different field on the way back, and the restore check would still pass.

        Returns:
            The binding SoA arrays plus the per-head ``r0_bind`` crossbridge reference, in a fixed order.
        """
        return (
            self.state.bound_d, self.state.seg_id_d, self.state.seg_a_d, self.state.seg_b_d,
            self.state.bary_t_d, self.state.abscissa_d, self.state.walk_dir_d, self._r0_bind_d,
        )

    def accumulate_candidate(self, actuator: NMIIActuatorView, port: FilamentMotorPortView) -> None:
        """Refresh live walk polarity, then scatter the attach-UNSTRAINED two-array crossbridge.

        Identical to the base connector except the crossbridge uses the per-head ``r0_bind`` reference instead
        of the scalar ``r0_xb``, so a freshly-bound head injects ZERO passive placement force (the coarse-mesh
        ``k_xb·⟨offset, ŵ⟩`` artifact is removed) and the load develops only from the power stroke.
        """
        self.launched_kernels = []
        self._launch(refresh_segment_barbed_kernel, dim=int(port.segment_node_a_d.shape[0]),
                     inputs=[port.position_d, port.segment_node_a_d, port.segment_node_b_d,
                             port.segment_polarity_d, self.state.seg_barbed_d], device=self._device)
        self.launched_kernels.append(refresh_segment_barbed_kernel)
        self._launch(_refresh_bound_walk_dir_kernel, dim=self.state.n_heads,
                     inputs=[port.position_d, self.state.bound_d, self.state.seg_id_d,
                             port.segment_node_a_d, port.segment_node_b_d, port.segment_polarity_d,
                             self.state.walk_dir_d], device=self._device)
        self.launched_kernels.append(_refresh_bound_walk_dir_kernel)
        self._launch(crossbridge_segment_split_r0bind_kernel, dim=self.state.n_heads,
                     inputs=[actuator.position_d, actuator.force_d, port.position_d, port.force_d,
                             self.head_node_d, self.state.bound_d, self.state.seg_a_d, self.state.seg_b_d,
                             self.state.bary_t_d, self.state.abscissa_d, self.state.walk_dir_d,
                             self.params.k_xb, self._r0_bind_d], device=self._device)
        self.launched_kernels.append(crossbridge_segment_split_r0bind_kernel)

    def compute_loads(self, actuator: NMIIActuatorView, port: FilamentMotorPortView) -> None:
        """Fill the per-head Hill/Bell loads from the converged geometry using the ``r0_bind`` reference."""
        self.launched_kernels = []
        self._launch(compute_head_loads_segment_split_r0bind_kernel, dim=self.state.n_heads,
                     inputs=[actuator.position_d, port.position_d, self.head_node_d, self.state.bound_d,
                             self.state.seg_a_d, self.state.seg_b_d, self.state.bary_t_d,
                             self.state.abscissa_d, self.state.walk_dir_d, self.params.k_xb,
                             self._r0_bind_d, self.state.loads_hill_d, self.state.loads_bell_d],
                     device=self._device)
        self.launched_kernels.append(compute_head_loads_segment_split_r0bind_kernel)

    def snapshot_candidate(self) -> None:
        """Snapshot the base binding SoA plus the per-head ``r0_bind`` reference (one candidate)."""
        super().snapshot_candidate()
        self._copy(self._r0_bind_snap_d, self._r0_bind_d)

    def rollback(self, accepted: wp.array) -> None:
        """Bit-restore the base binding SoA plus ``r0_bind`` (the gated commit never mutates on a rejected step)."""
        super().rollback(accepted)
        self._copy(self._r0_bind_d, self._r0_bind_snap_d)

    def propose_events(self, rates: object, dt_phys: float, rng_seed: int, neighbors: object) -> None:
        """Fill the accepted-step attach scratch from the live cortex port (no force, no host state).

        Re-binds the live segment endpoints (INTEGRATION item 3) then runs the injected point-to-segment query
        (INTEGRATION item 2): each free head's candidate segment/barycentric-fraction/barbed-polarity is written
        into the connector-owned ``query_*`` scratch that :meth:`commit_irreversible`'s attach kernel consumes
        under the accepted predicate.  A rejected outer step re-proposes identically (device RNG folds the
        accepted-step epoch), so the candidate is reproducible.
        """
        self.bind_target_topology(self._port.segment_node_a_d, self._port.segment_node_b_d)
        self._query.fill_attach_scratch(self._actuator_view(), self._port)

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Accepted-predicated attach-UNSTRAINED → Hill-step/DETACH → epoch KMC (detach defaults to catch-slip).

        The attach kernel is the ``r0bind`` variant: on a successful bind it records the per-head zero-strain
        reference ``r0_bind[h] = dot(x_att − x_head, walk_dir)`` from the CONVERGED geometry, so the crossbridge
        is unstrained at the binding instant (the coarse-mesh placement artifact is removed).  ``CATCH_SLIP``
        (default) uses the physiological Kovacs 2007 catch-slip detach; ``SLIP`` uses the pure Bell slip.  Both
        detach laws read the ``r0bind``-corrected loads, and the epoch advances only on acceptance.
        """
        self.launched_kernels = []
        seed_attach = wp.int32(int(rng_seed) & 0x7FFFFFFF)
        seed_detach = wp.int32((int(rng_seed) ^ 0x5BD1E995) & 0x7FFFFFFF)
        actuator_pos = self._actuator_view().position_d
        port_pos = self._port.position_d
        self._launch(
            attach_segment_gated_r0bind_kernel, dim=self.state.n_heads,
            inputs=[accepted, self.state.rng_epoch_d, self.state.bound_d, self.state.seg_id_d,
                    self.state.seg_a_d, self.state.seg_b_d, self.state.bary_t_d, self.state.abscissa_d,
                    self.state.walk_dir_d, self._r0_bind_d, self.state.query_seg_id_d, self.state.query_t_d,
                    self.state.query_barbed_d, self._seg_node_a(), self._seg_node_b(),
                    actuator_pos, port_pos, self.head_node_d, self.params,
                    wp.float64(dt_phys), seed_attach], device=self._device,
        )
        self.launched_kernels.append(attach_segment_gated_r0bind_kernel)
        if self.detach_kinetics is SegmentDetachKinetics.SLIP:
            self._launch(
                step_detach_segment_gated_kernel, dim=self.state.n_heads,
                inputs=[accepted, self.state.rng_epoch_d, self.state.bound_d, self.state.seg_id_d,
                        self.state.seg_a_d, self.state.seg_b_d, self.state.bary_t_d, self.state.abscissa_d,
                        self.state.walk_dir_d, self.state.loads_hill_d, self.state.loads_bell_d, self.params,
                        wp.float64(dt_phys), seed_detach], device=self._device,
            )
            self.launched_kernels.append(step_detach_segment_gated_kernel)
        else:
            self._launch(
                step_detach_segment_gated_catch_slip_kernel, dim=self.state.n_heads,
                inputs=[accepted, self.state.rng_epoch_d, self.state.bound_d, self.state.seg_id_d,
                        self.state.seg_a_d, self.state.seg_b_d, self.state.bary_t_d, self.state.abscissa_d,
                        self.state.walk_dir_d, self.state.loads_hill_d, self.state.loads_bell_d, self.params,
                        self._catch_slip, wp.float64(dt_phys), seed_detach], device=self._device,
            )
            self.launched_kernels.append(step_detach_segment_gated_catch_slip_kernel)
        self._launch(
            _increment_epoch_if_accepted_kernel, dim=1,
            inputs=[accepted, self.state.rng_epoch_d], device=self._device,
        )
        self.launched_kernels.append(_increment_epoch_if_accepted_kernel)


#: The original name of the (now target-agnostic) connector, kept for the landed cortex call sites.  It is the
#: SAME class — a cortex binding is just its default ``name``/``component_b``, so cortex behaviour is unchanged.
CortexMotorConnector = FilamentMotorConnector

#: The head/crossbridge constants are properties of the NMII HEAD, not of the filament it pulls on, so the same
#: config serves every MOTOR edge.  This alias exists so a non-cortex call site (e.g. ``nmii_sf_motor``) reads
#: correctly; it is the SAME class as :class:`CortexMotorParams`, with the same REQUIRED-PARAM discipline.
NMIIMotorParams = CortexMotorParams


# ── The wired slice: one accepted physical step over cortex ⟷ NMII. ───────────────────────────────────
@dataclass(slots=True)
class CortexMotorSlice:
    """One vertical column of the cell graph wired for the accepted-step transaction.

    Holds the two participants (the NMII state-owner and the ``nmii_cortex_motor`` connector), the cortex
    bind-target port, the device event clock, and the top-level :class:`CellTransaction`.  Constructed from
    already-assembled parts, so it is CPU-importable and the structural gates build it with recording doubles;
    :func:`build_cortex_motor_slice` assembles the real Warp-resident parts on the CUDA lane.
    """

    transaction: CellTransaction
    actuator_state: NMIIActuatorStateOwner
    connector: CortexMotorConnector
    port: FilamentMotorPortView
    clock: Any
    ledger: Any | None = None
    tol_sq_d: wp.array | None = None

    def accumulate(self) -> None:
        """One candidate-force pass the caller's inner solve invokes each iteration (after zeroing forces).

        Adds (1) the minifilament's own internal backbone/head-arm mechanics into the ``nmii``-owned force
        array and (2) the split crossbridge — ``+f`` into the head array, ``−(1−t)f``/``−t·f`` into the two
        cortex port nodes (Newton's 3rd law across two never-merged arrays).  No fixed prestress is injected;
        the crossbridge force is the reaction of whatever heads the KMC has bound so far.
        """
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
    ) -> None:
        """Run exactly one accepted physical step of the dynamic cortical-motor runtime.

        Args:
            solve: the caller-owned inner mechanical solve (zero forces → :meth:`accumulate` + cortex passive
                relaxation → integrate to convergence).  Runs after events are proposed and owns all physics.
            dt_phys: the outer physical timestep [s]; finite and positive.
            accepted_d: the device acceptance predicate.  If ``None`` the balance gate decides — then a
                ``ledger`` with a ``tol_sq_d`` must be available.
            ledger / tol_sq_d: optional balance-gate inputs; default to the slice's stored ones.

        The per-head Hill/Bell loads are computed on the converged geometry (after ``solve``, before the KMC
        commit consumes them).  The predicate stays device-resident throughout.
        """
        ledger = ledger if ledger is not None else self.ledger
        tol_sq_d = tol_sq_d if tol_sq_d is not None else self.tol_sq_d

        def _inner() -> None:
            solve()
            # Per-head loads on the converged geometry, before commit_irreversible consumes them.
            self.connector.compute_loads(self.actuator_state.geometry(), self.port)

        self.transaction.step(
            dt_phys=dt_phys, solve=_inner, accepted_d=accepted_d, ledger=ledger, tol_sq_d=tol_sq_d,
        )

    def bound_head_count(self) -> int:
        """Host read of the bound-head population for OUT-OF-LOOP diagnostics only (device→host readback).

        This is a viz/telemetry helper the caller may call BETWEEN steps to watch the emergent bound fraction
        rise; it must never be called inside the inner solve (I0-A: no authoritative GPU→CPU roundtrip in the
        hot loop).
        """
        return int(self.connector.state.bound_d.numpy().sum())


# ── CUDA-lane builders (allocate device memory; called on the gbook A5000, not on the dev Mac). ────────
def allocate_segment_connector_state(
    n_heads: int, n_segments: int, *, device: str | None = None,
) -> SegmentConnectorState:
    """Allocate a free (all-unbound) :class:`SegmentConnectorState` + its snapshot twins + query scratch.

    Args:
        n_heads: number of NMII heads owned by this MOTOR connector.
        n_segments: number of live actin segments in the target cortex port.
        device: CUDA device string (never a hard-coded id in production; the caller passes the resolved one).
    """
    if n_heads <= 0 or n_segments <= 0:
        raise ValueError("n_heads and n_segments must be positive")

    def z(n: int, dtype: object) -> wp.array:
        return wp.zeros(n, dtype=dtype, device=device)

    def full(n: int, value: int, dtype: object) -> wp.array:
        return wp.full(n, value, dtype=dtype, device=device)

    return SegmentConnectorState(
        n_heads=int(n_heads), n_segments=int(n_segments),
        bound_d=z(n_heads, wp.int32), seg_id_d=full(n_heads, -1, wp.int32),
        seg_a_d=full(n_heads, -1, wp.int32), seg_b_d=full(n_heads, -1, wp.int32),
        bary_t_d=z(n_heads, wp.float64), abscissa_d=z(n_heads, wp.float64),
        walk_dir_d=z(n_heads, wp.vec3d), rng_epoch_d=z(1, wp.int32),
        loads_hill_d=z(n_heads, wp.float64), loads_bell_d=z(n_heads, wp.float64),
        seg_barbed_d=z(n_segments, wp.vec3d),
        query_seg_id_d=full(n_heads, -1, wp.int32), query_t_d=z(n_heads, wp.float64),
        query_barbed_d=z(n_heads, wp.vec3d),
        bound_snap_d=z(n_heads, wp.int32), seg_id_snap_d=full(n_heads, -1, wp.int32),
        seg_a_snap_d=full(n_heads, -1, wp.int32), seg_b_snap_d=full(n_heads, -1, wp.int32),
        bary_t_snap_d=z(n_heads, wp.float64), abscissa_snap_d=z(n_heads, wp.float64),
        walk_dir_snap_d=z(n_heads, wp.vec3d), rng_epoch_snap_d=z(1, wp.int32),
    )


def build_nmii_actuator_state(
    *,
    n_minifilaments: int,
    n_heads: int,
    position_d: wp.array,
    minifilament_offset_d: wp.array,
    active_minifilament_d: wp.array,
    minifilament_id_d: wp.array,
    particle_role_d: wp.array,
    head_node_d: wp.array,
    head_id_d: wp.array,
    head_side_d: wp.array,
    backbone_bonds_d: wp.array,
    head_bonds_d: wp.array,
    backbone_angles_d: wp.array,
    head_arm_angles_d: wp.array,
    k_backbone_pn_per_um: float,
    r0_backbone_um: float,
    k_head_spring_pn_per_um: float,
    r0_head_um: float,
    bending: NMIIBendingStiffness,
    device: str | None = None,
) -> NMIIActuatorStateOwner:
    """Assemble a validated NMII state-owner from a straddle-placed minifilament population (CUDA lane).

    The caller (the Lead) supplies the minifilament geometry + backbone/head bond/angle topology from the
    landed minifilament builder + straddle placement; this allocates the per-head force/load/ATP/epoch arrays
    and the snapshot twins, binds the real :class:`BackboneArmMechanics`, and wires the concrete state
    transaction/ledger delegates.  It never lowers a biological count — sizes come from the passed population.
    """
    n_particles = int(position_d.shape[0])
    force_d = wp.zeros(n_particles, dtype=wp.vec3d, device=device)
    head_load_hill_d = wp.zeros(n_heads, dtype=wp.float64, device=device)
    head_load_bell_d = wp.zeros(n_heads, dtype=wp.float64, device=device)
    atp_cycle_state_d = wp.zeros(n_heads, dtype=wp.int32, device=device)
    atp_consumed_d = wp.zeros(1, dtype=wp.int64, device=device)
    actuator_epoch_d = wp.zeros(1, dtype=wp.int32, device=device)

    mechanics = BackboneArmMechanics(
        n_particles=n_particles, n_heads=int(n_heads),
        backbone_bonds_d=backbone_bonds_d, head_bonds_d=head_bonds_d,
        backbone_angles_d=backbone_angles_d, head_arm_angles_d=head_arm_angles_d,
        k_backbone_pn_per_um=k_backbone_pn_per_um, r0_backbone_um=r0_backbone_um,
        k_head_spring_pn_per_um=k_head_spring_pn_per_um, r0_head_um=r0_head_um,
        bending=bending, device=device,
    )
    transaction = _ActuatorStateTransaction(
        position_d=position_d, active_minifilament_d=active_minifilament_d,
        atp_cycle_state_d=atp_cycle_state_d, atp_consumed_d=atp_consumed_d,
        actuator_epoch_d=actuator_epoch_d,
        position_snap_d=wp.zeros(n_particles, dtype=wp.vec3d, device=device),
        active_snap_d=wp.zeros(int(n_minifilaments), dtype=wp.int32, device=device),
        atp_state_snap_d=wp.zeros(n_heads, dtype=wp.int32, device=device),
        atp_consumed_snap_d=wp.zeros(1, dtype=wp.int64, device=device),
        epoch_snap_d=wp.zeros(1, dtype=wp.int32, device=device),
        device=device,
    )
    return NMIIActuatorStateOwner(
        n_minifilaments=int(n_minifilaments), n_heads=int(n_heads),
        position_d=position_d, force_d=force_d, minifilament_offset_d=minifilament_offset_d,
        active_minifilament_d=active_minifilament_d, minifilament_id_d=minifilament_id_d,
        particle_role_d=particle_role_d, head_node_d=head_node_d, head_id_d=head_id_d,
        head_side_d=head_side_d, head_load_hill_d=head_load_hill_d, head_load_bell_d=head_load_bell_d,
        atp_cycle_state_d=atp_cycle_state_d, atp_consumed_d=atp_consumed_d,
        actuator_epoch_d=actuator_epoch_d, mechanics=mechanics, transaction=transaction,
        ledger=_ActuatorStateLedger(),
    )


def build_cortex_motor_slice(
    *,
    actuator_state: NMIIActuatorStateOwner,
    port: FilamentMotorPortView,
    params: CortexMotorParams,
    base_seed: int,
    device: str,
    max_segment_length_um: float,
    connector_state: SegmentConnectorState | None = None,
    clock: Any | None = None,
    ledger: Any | None = None,
    tol_sq_d: wp.array | None = None,
    population_ledgers: tuple[Any, ...] = (),
    architecture: CellArchitecture | None = None,
) -> CortexMotorSlice:
    """Assemble the CUDA-lane dynamic cortex-motor slice: two participants under one accepted-step clock.

    Binds the NMII state-owner as the ``nmii`` component and a :class:`CortexMotorConnector` as the
    ``nmii_cortex_motor`` connector on a graph that already carries the protrusion + NMII contracts, then wires
    a :class:`CellWorldTransaction` (``require_complete=False`` — this is the first vertical slice, not the
    whole cell) and a :class:`CellTransaction`.  The cortex is a bind-target PORT, never a participant.

    Args:
        actuator_state: the NMII state-owner (see :func:`build_nmii_actuator_state`).
        port: the live cortex actin bind-target port (its own ``seg_node_a/b`` + polarity + material coords).
        params: the sourced-or-GAP kinetic/mechanical constants (:class:`CortexMotorParams`).
        base_seed: the run's fixed host RNG seed (nonnegative int).
        device: the resolved CUDA device string (never a hard-coded id).
        max_segment_length_um: longest cortex segment [µm]; sets the attach-query grid radius.
        connector_state: optional pre-allocated binding SoA; allocated from the head/segment counts if omitted.
        clock / ledger / tol_sq_d: optional; a :class:`WarpEventClock` is built from ``base_seed`` if omitted.
        population_ledgers: optional disjoint-population ledgers asserted at build + each accepted step.
        architecture: optional composition; defaults to the protrusion+NMII-promoted reference architecture.

    Returns:
        A :class:`CortexMotorSlice` whose :meth:`~CortexMotorSlice.step` runs one accepted physical step.
    """
    from aleph.engine.events import make_event_clock  # CUDA-lane import (allocates device scalars).

    if architecture is None:
        # Local: `protrusion` is 1,107 lines this lane needs only to name the DEFAULT architecture, and
        # importing it at module scope put a lamellipodium the cortex-motor run never builds into every
        # native run's import set (2026-07-28).
        from aleph.engine.protrusion import proposed_protrusion_architecture

        architecture = proposed_nmii_architecture(proposed_protrusion_architecture(reference_cell_architecture()))

    n_heads = int(actuator_state.n_heads)
    n_segments = int(port.segment_node_a_d.shape[0])
    if connector_state is None:
        connector_state = allocate_segment_connector_state(n_heads, n_segments, device=device)
    elif int(connector_state.n_heads) != n_heads:
        raise ValueError("connector_state head count must match the NMII state-owner")

    hand_params = params.to_hand_params()
    catch_slip = (
        params.to_catch_slip_params()
        if params.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP else None
    )
    query = SegmentAttachQuery(
        state=connector_state, head_node_d=actuator_state.head_node_d, port=port,
        capture_radius=params.capture_radius, max_segment_length_um=max_segment_length_um, device=device,
    )
    # Per-head attach-unstrained crossbridge reference + its snapshot twin (owned by the connector).
    r0_bind_d = wp.zeros(n_heads, dtype=wp.float64, device=device)
    r0_bind_snap_d = wp.zeros(n_heads, dtype=wp.float64, device=device)
    connector = CortexMotorConnector(
        state=connector_state, head_node_d=actuator_state.head_node_d, params=hand_params,
        catch_slip=catch_slip, query=query, actuator_view=actuator_state.geometry, port=port,
        r0_bind_d=r0_bind_d, r0_bind_snap_d=r0_bind_snap_d,
        detach_kinetics=params.detach_kinetics, device=device,
    )

    actor = CellActor(architecture)
    actor.bind_component(NMII_COMPONENT, actuator_state)
    actor.bind_connector(NMII_CORTEX_MOTOR, connector)
    world = CellWorldTransaction(actor, require_complete=False)

    if clock is None:
        clock = make_event_clock(base_seed=base_seed, device=device)
    transaction = CellTransaction(world, clock, population_ledgers=tuple(population_ledgers))
    return CortexMotorSlice(
        transaction=transaction, actuator_state=actuator_state, connector=connector, port=port,
        clock=clock, ledger=ledger, tol_sq_d=tol_sq_d,
    )
