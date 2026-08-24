"""CPU-green structural + numeric gates for the ``nmii_sf_motor`` dynamic SF-motor vertical slice.

Proves — with recording CUDA doubles and the landed NumPy crossbridge oracle, no device — that
:class:`~aleph.engine.sf_motor_slice.SFMotorSlice` binds the SF MOTOR edge the graph has declared since
the architecture landed, and that it is a real THREE-participant transaction rather than a cortex copy:

* the connector is bound at the ``nmii_sf_motor`` slot with ``component_b == "sf_arc"``, and handing it a
  cortex port RAISES (no silent mis-binding of the wrong target);
* the ``sf_arc`` state owner is a real participant that LAUNCHES its own passive rod-cable mechanics, so one
  candidate pass runs SF passive → NMII internal → the split crossbridge, in that order;
* the split crossbridge conserves momentum across the two never-merged arrays
  (``f_head + f_seg_a + f_seg_b = 0``) and is UNSTRAINED at the binding instant (``r0_bind``);
* the sarcomere geometry the motor needs is DERIVED from the minifilament, not chosen;
* the cortex connector's default binding is unchanged by the target-agnostic generalisation.

The CUDA-lane numeric behaviour (bound-fraction growth, emergent SF tension) is the Lead's native gate.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest
import warp as wp

from aleph.engine.actor import CellActor
from aleph.engine.composed_native import SFArcStateOwner
from aleph.engine.contracts import reference_cell_architecture
from aleph.engine.cortex_motor_slice import (
    CortexMotorConnector,
    FilamentMotorConnector,
    NMIIMotorParams,
    _ActuatorStateLedger,
    _ActuatorStateTransaction,
)
from aleph.engine.nmii_actuator import (
    CORTEX_COMPONENT,
    NMII_COMPONENT,
    NMII_CORTEX_MOTOR,
    NMII_SF_MOTOR,
    SF_COMPONENT,
    BackboneArmMechanics,
    FilamentMotorPortView,
    NMIIActuatorStateOwner,
    NMIIBendingStiffness,
    SegmentConnectorState,
    proposed_nmii_architecture,
)
from aleph.engine.protrusion import proposed_protrusion_architecture
from aleph.engine.sf_mechanics import SFFilamentMechanics, build_sf_mechanics_topology
from aleph.engine.sf_motor_slice import (
    SFMotorSlice,
    build_sf_motor_port,
    sf_sarcomere_geometry_for,
)
from aleph.engine.sf_population import build_sf_arc_population
from aleph.engine.world import CellWorldTransaction
from aleph.components.motor.minifilament_topology import MinifilamentTopology
from aleph.components.motor.segment_motor import (
    SegmentDetachKinetics,
    crossbridge_segment_split_r0bind_reference,
    r0_bind_at_attach_reference,
)

_TOPOLOGY = MinifilamentTopology(
    n_bb=14, n_heads_per_side=10, backbone_length_um=0.301, head_offset_um=0.200)


# ── CUDA metadata doubles (no host physics is evaluated by the structural gates) ──────────────────────
@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    ptr: int
    shape: tuple[int, ...]
    dtype: object
    device: _FakeDevice = _FakeDevice()


def _fa(ptr: int, n: int, dtype: object) -> _FakeArray:
    return _FakeArray(ptr, (n,), dtype)


def _fa2(ptr: int, rows: int, cols: int, dtype: object) -> _FakeArray:
    return _FakeArray(ptr, (rows, cols), dtype)


@dataclass(slots=True)
class _RecordingLauncher:
    calls: list = field(default_factory=list)

    def __call__(self, kernel, *, dim, inputs, outputs=None, device=None) -> None:
        self.calls.append((kernel, int(dim), tuple(inputs), tuple(outputs or ()), device))

    @property
    def kernels(self) -> list:
        return [call[0] for call in self.calls]


@dataclass(slots=True)
class _RecordingCopy:
    calls: list = field(default_factory=list)

    def __call__(self, dst, src) -> None:
        self.calls.append((dst, src))


@dataclass(slots=True)
class _SpyQuery:
    calls: list = field(default_factory=list)

    def fill_attach_scratch(self, actuator, port) -> None:
        self.calls.append((actuator, port))


class _HandParamsDouble:
    k_xb = wp.float64(1000.0)
    r0_xb = wp.float64(0.0)
    k_on = wp.float64(50.0)


_N_NODES_SF = 12
_N_SEG_SF = 6
_N_HEADS = 4


def _sf_owner(launch: _RecordingLauncher, copy: _RecordingCopy) -> SFArcStateOwner:
    """A real :class:`SFArcStateOwner` over CUDA-metadata doubles (its kernels flow through the recorder)."""
    mechanics = SFFilamentMechanics(
        device="cuda:0", n_nodes=_N_NODES_SF,
        links_d=_fa2(400, _N_SEG_SF, 2, wp.int32),
        link_k_d=_fa(401, _N_SEG_SF, wp.float64),
        link_r0_d=_fa(402, _N_SEG_SF, wp.float64),
        bend_triples_d=_fa2(403, 4, 3, wp.int32),
        bend_alpha_d=_fa(404, 4, wp.float64),
        link_kernel=object(), bending_kernel=object(), launch=launch,
    )
    return SFArcStateOwner(
        name=SF_COMPONENT, device="cuda:0", n_filaments=2, n_nodes=_N_NODES_SF,
        position_d=_fa(410, _N_NODES_SF, wp.vec3d), force_d=_fa(411, _N_NODES_SF, wp.vec3d),
        position_snap_d=_fa(412, _N_NODES_SF, wp.vec3d),
        epoch_d=_fa(413, 1, wp.int32), epoch_snap_d=_fa(414, 1, wp.int32),
        mechanics=mechanics, arc_connector=None, launch=launch, copy=copy,
    )


def _actuator(launch: _RecordingLauncher, copy: _RecordingCopy) -> NMIIActuatorStateOwner:
    position_d = _fa(100, 10, wp.vec3d)
    active_d = _fa(103, 2, wp.int32)
    atp_state_d = _fa(111, _N_HEADS, wp.int32)
    atp_consumed_d = _fa(112, 1, wp.int64)
    epoch_d = _fa(113, 1, wp.int32)
    mechanics = BackboneArmMechanics(
        n_particles=10, n_heads=_N_HEADS,
        backbone_bonds_d=_fa(120, 2, wp.int32), head_bonds_d=_fa(121, _N_HEADS, wp.int32),
        backbone_angles_d=_fa(122, 1, wp.int32), head_arm_angles_d=_fa(123, _N_HEADS, wp.int32),
        k_backbone_pn_per_um=1.0e3, r0_backbone_um=0.023, k_head_spring_pn_per_um=1.0e2, r0_head_um=0.200,
        bending=NMIIBendingStiffness(persistence_length_um=1.0, segment_len_um=0.023,
                                     k_xb_pn_per_um=1.0e3, r0_head_um=0.200),
        launch=launch,
    )
    transaction = _ActuatorStateTransaction(
        position_d=position_d, active_minifilament_d=active_d, atp_cycle_state_d=atp_state_d,
        atp_consumed_d=atp_consumed_d, actuator_epoch_d=epoch_d,
        position_snap_d=_fa(200, 10, wp.vec3d), active_snap_d=_fa(203, 2, wp.int32),
        atp_state_snap_d=_fa(211, _N_HEADS, wp.int32), atp_consumed_snap_d=_fa(212, 1, wp.int64),
        epoch_snap_d=_fa(213, 1, wp.int32), launch=launch, copy=copy,
    )
    return NMIIActuatorStateOwner(
        n_minifilaments=2, n_heads=_N_HEADS, position_d=position_d, force_d=_fa(101, 10, wp.vec3d),
        minifilament_offset_d=_fa(102, 3, wp.int32), active_minifilament_d=active_d,
        minifilament_id_d=_fa(104, 2, wp.int64), particle_role_d=_fa(105, 10, wp.int32),
        head_node_d=_fa(106, _N_HEADS, wp.int32), head_id_d=_fa(107, _N_HEADS, wp.int64),
        head_side_d=_fa(108, _N_HEADS, wp.int32), head_load_hill_d=_fa(109, _N_HEADS, wp.float64),
        head_load_bell_d=_fa(110, _N_HEADS, wp.float64), atp_cycle_state_d=atp_state_d,
        atp_consumed_d=atp_consumed_d, actuator_epoch_d=epoch_d,
        mechanics=mechanics, transaction=transaction, ledger=_ActuatorStateLedger(),
    )


def _port(component: str, base: int) -> FilamentMotorPortView:
    return FilamentMotorPortView(
        component=component, position_d=_fa(base, _N_NODES_SF, wp.vec3d),
        force_d=_fa(base + 1, _N_NODES_SF, wp.vec3d),
        segment_node_a_d=_fa(base + 2, _N_SEG_SF, wp.int32),
        segment_node_b_d=_fa(base + 3, _N_SEG_SF, wp.int32),
        segment_polarity_d=_fa(base + 4, _N_SEG_SF, wp.int32),
        persistent_filament_id_d=_fa(base + 5, _N_SEG_SF, wp.int64),
        material_s0_d=_fa(base + 6, _N_SEG_SF, wp.float64),
        material_s1_d=_fa(base + 7, _N_SEG_SF, wp.float64),
        topology_epoch_d=_fa(base + 8, 1, wp.int32),
    )


def _connector_state() -> SegmentConnectorState:
    return SegmentConnectorState(
        n_heads=_N_HEADS, n_segments=_N_SEG_SF,
        bound_d=_fa(800, _N_HEADS, wp.int32), seg_id_d=_fa(801, _N_HEADS, wp.int32),
        seg_a_d=_fa(802, _N_HEADS, wp.int32), seg_b_d=_fa(803, _N_HEADS, wp.int32),
        bary_t_d=_fa(804, _N_HEADS, wp.float64), abscissa_d=_fa(805, _N_HEADS, wp.float64),
        walk_dir_d=_fa(806, _N_HEADS, wp.vec3d), rng_epoch_d=_fa(807, 1, wp.int32),
        loads_hill_d=_fa(808, _N_HEADS, wp.float64), loads_bell_d=_fa(809, _N_HEADS, wp.float64),
        seg_barbed_d=_fa(810, _N_SEG_SF, wp.vec3d), query_seg_id_d=_fa(811, _N_HEADS, wp.int32),
        query_t_d=_fa(812, _N_HEADS, wp.float64), query_barbed_d=_fa(813, _N_HEADS, wp.vec3d),
        bound_snap_d=_fa(820, _N_HEADS, wp.int32), seg_id_snap_d=_fa(821, _N_HEADS, wp.int32),
        seg_a_snap_d=_fa(822, _N_HEADS, wp.int32), seg_b_snap_d=_fa(823, _N_HEADS, wp.int32),
        bary_t_snap_d=_fa(824, _N_HEADS, wp.float64), abscissa_snap_d=_fa(825, _N_HEADS, wp.float64),
        walk_dir_snap_d=_fa(826, _N_HEADS, wp.vec3d), rng_epoch_snap_d=_fa(827, 1, wp.int32),
    )


def _sf_connector(launch, copy, query, actuator, port) -> FilamentMotorConnector:
    return FilamentMotorConnector(
        state=_connector_state(), head_node_d=actuator.head_node_d, params=_HandParamsDouble(),
        catch_slip=object(), query=query, actuator_view=actuator.geometry, port=port,
        r0_bind_d=_fa(850, _N_HEADS, wp.float64), r0_bind_snap_d=_fa(851, _N_HEADS, wp.float64),
        detach_kinetics=SegmentDetachKinetics.CATCH_SLIP,
        name=NMII_SF_MOTOR, component_b=SF_COMPONENT, launch=launch, copy=copy,
    )


def _architecture():
    return proposed_nmii_architecture(proposed_protrusion_architecture(reference_cell_architecture()))


def _slice() -> tuple[SFMotorSlice, _RecordingLauncher, _RecordingCopy, _SpyQuery]:
    launch, copy, query = _RecordingLauncher(), _RecordingCopy(), _SpyQuery()
    sf_owner = _sf_owner(launch, copy)
    actuator = _actuator(launch, copy)
    port = _port(SF_COMPONENT, 300)
    connector = _sf_connector(launch, copy, query, actuator, port)
    actor = CellActor(_architecture())
    actor.bind_component(SF_COMPONENT, sf_owner)
    actor.bind_component(NMII_COMPONENT, actuator)
    actor.bind_connector(NMII_SF_MOTOR, connector)
    world = CellWorldTransaction(actor, require_complete=False)

    @dataclass(slots=True)
    class _SpyClock:
        base_seed: int = 7
        epoch: int = 0
        advances: list = field(default_factory=list)

        def advance(self, accepted, dt_phys) -> None:
            self.advances.append((accepted, dt_phys))

    from aleph.engine.transaction import CellTransaction

    transaction = CellTransaction(world, _SpyClock(), population_ledgers=())
    return (SFMotorSlice(transaction=transaction, sf_owner=sf_owner, actuator_state=actuator,
                         connector=connector, port=port, clock=transaction.clock, _launch=launch),
            launch, copy, query)


# ── the gates ─────────────────────────────────────────────────────────────────────────────────────────
def test_sarcomere_geometry_is_derived_from_the_minifilament() -> None:
    """Magic-Number Block: both sarcomere lengths follow from the motor's own rigid geometry."""
    geometry = sf_sarcomere_geometry_for(_TOPOLOGY)
    assert geometry["sarcomere_lateral_um"] == pytest.approx(2.0 * _TOPOLOGY.head_offset_um)
    assert geometry["sarcomere_overlap_um"] == pytest.approx(_TOPOLOGY.backbone_length_um)


def test_connector_binds_the_sf_edge_and_refuses_a_cortex_port() -> None:
    """Ownership: the SF MOTOR edge targets ``sf_arc``; a wrong-target port is refused, never mis-bound."""
    launch, copy, query = _RecordingLauncher(), _RecordingCopy(), _SpyQuery()
    actuator = _actuator(launch, copy)
    connector = _sf_connector(launch, copy, query, actuator, _port(SF_COMPONENT, 300))
    assert connector.name == NMII_SF_MOTOR
    assert connector.component_a == NMII_COMPONENT
    assert connector.component_b == SF_COMPONENT
    assert connector.individual_head_state and connector.per_head_bell and connector.hill_force_velocity
    assert connector.aggregate_or_two_anchor is False

    with pytest.raises(ValueError, match="was handed a"):
        _sf_connector(launch, copy, query, actuator, _port(CORTEX_COMPONENT, 500))


def test_cortex_connector_default_binding_is_unchanged() -> None:
    """Additive: the target-agnostic generalisation leaves the landed cortex binding byte-identical."""
    launch, copy, query = _RecordingLauncher(), _RecordingCopy(), _SpyQuery()
    actuator = _actuator(launch, copy)
    cortex = CortexMotorConnector(
        state=_connector_state(), head_node_d=actuator.head_node_d, params=_HandParamsDouble(),
        catch_slip=object(), query=query, actuator_view=actuator.geometry,
        port=_port(CORTEX_COMPONENT, 500),
        r0_bind_d=_fa(850, _N_HEADS, wp.float64), r0_bind_snap_d=_fa(851, _N_HEADS, wp.float64),
        detach_kinetics=SegmentDetachKinetics.CATCH_SLIP, launch=launch, copy=copy,
    )
    assert CortexMotorConnector is FilamentMotorConnector
    assert cortex.name == NMII_CORTEX_MOTOR
    assert cortex.component_b == CORTEX_COMPONENT


def test_slice_has_three_participants_including_the_sf_owner() -> None:
    """The SF component is a real participant here, not a port-only bind target as in the cortex slice."""
    slice_, _, _, _ = _slice()
    participants = slice_.participants
    assert len(participants) == 3
    assert slice_.sf_owner in participants
    assert slice_.actuator_state in participants
    assert slice_.connector in participants
    actor = slice_.transaction.world.actor
    assert actor.component_runtime(SF_COMPONENT) is slice_.sf_owner
    assert actor.component_runtime(NMII_COMPONENT) is slice_.actuator_state
    assert actor.connector_runtime(NMII_SF_MOTOR) is slice_.connector
    # the graph contract for this edge: kinetic, accepted-step committed, bidirectional + adjoint.
    contract = actor.connector_contract(NMII_SF_MOTOR)
    assert contract.kinetics and contract.commit_on_accept
    assert contract.bidirectional and contract.adjoint_transfer_required
    assert {contract.component_a, contract.component_b} == {NMII_COMPONENT, SF_COMPONENT}


def test_candidate_pass_launches_sf_passive_then_nmii_internal_then_crossbridge() -> None:
    """One accumulate = SF rod-cable + minifilament internal + the split crossbridge, in that order."""
    slice_, launch, _, _ = _slice()
    launch.calls.clear()
    slice_.accumulate()

    kernels = launch.kernels
    assert len(kernels) >= 6
    # SF passive first: its own link + bending kernels over the SF-owned arrays.
    assert kernels[0] is slice_.sf_owner.mechanics.link_kernel
    assert kernels[1] is slice_.sf_owner.mechanics.bending_kernel
    sf_link_call = launch.calls[0]
    assert sf_link_call[2][0] is slice_.sf_owner.position_d          # SF pos in
    assert sf_link_call[1] == _N_SEG_SF
    # then the crossbridge, which must address BOTH arrays (split ownership, never merged).
    crossbridge = launch.calls[-1]
    inputs = crossbridge[2]
    assert slice_.actuator_state.position_d in inputs and slice_.actuator_state.force_d in inputs
    assert slice_.port.position_d in inputs and slice_.port.force_d in inputs
    # split ownership: the actuator's and the target's arrays are DISTINCT storage, never one merged array.
    assert slice_.actuator_state.position_d.ptr != slice_.port.position_d.ptr
    assert slice_.actuator_state.force_d.ptr != slice_.port.force_d.ptr


def test_split_crossbridge_conserves_momentum_and_is_unstrained_at_bind() -> None:
    """Newton's 3rd law across the two never-merged arrays + zero force at the binding instant (r0_bind)."""
    rng = np.random.default_rng(4)
    actuator_pos = rng.normal(size=(3, 3))
    port_pos = rng.normal(size=(4, 3))
    walk_dir = np.array([0.0, 0.0, 1.0])
    head_node, seg_a, seg_b, bary_t, k_xb = 1, 2, 3, 0.35, 1000.0

    r0_bind = r0_bind_at_attach_reference(
        actuator_pos=actuator_pos, port_pos=port_pos, head_node=head_node,
        seg_a=seg_a, seg_b=seg_b, bary_t=bary_t, walk_dir=walk_dir)
    at_bind = crossbridge_segment_split_r0bind_reference(
        actuator_pos, port_pos, head_node, seg_a, seg_b, bary_t, 0.0, walk_dir, k_xb, r0_bind)
    assert np.allclose(at_bind["f_head"], 0.0, atol=1e-12)           # unstrained at bind

    stroked = crossbridge_segment_split_r0bind_reference(
        actuator_pos, port_pos, head_node, seg_a, seg_b, bary_t, 0.004, walk_dir, k_xb, r0_bind)
    total = stroked["f_head"] + stroked["f_seg_a"] + stroked["f_seg_b"]
    assert np.allclose(total, 0.0, atol=1e-12)                       # momentum-conserving split
    assert np.linalg.norm(stroked["f_head"]) == pytest.approx(k_xb * 0.004, rel=1e-9)


def test_port_builder_refuses_a_topology_without_per_segment_metadata() -> None:
    """The port needs identity/polarity emitted IN LINK ORDER; a stale topology raises before any upload."""
    population = build_sf_arc_population(
        n_ventral=2, n_dorsal=0, n_arc=0, n_cap=0, n_per_fiber=5,
        **sf_sarcomere_geometry_for(_TOPOLOGY))
    topology = build_sf_mechanics_topology(population, k_axial_pn_per_um=1000.0)
    assert topology.seg_polarity.size == topology.n_links > 0        # the fresh builder emits it

    topology.seg_polarity = np.zeros(0, np.int32)                    # simulate a pre-metadata topology
    with pytest.raises(ValueError, match="rebuild it with build_sf_mechanics_topology"):
        build_sf_motor_port(position_d=_fa(1, 4, wp.vec3d), force_d=_fa(2, 4, wp.vec3d),
                            topology=topology, device="cuda:0")


def test_segment_metadata_matches_the_link_order_and_the_barbed_convention() -> None:
    """MEASURED: per-segment polarity/identity/arclength line up with the links they describe."""
    population = build_sf_arc_population(
        n_ventral=2, n_dorsal=1, n_arc=1, n_cap=1, n_per_fiber=6,
        **sf_sarcomere_geometry_for(_TOPOLOGY))
    topology = build_sf_mechanics_topology(population, k_axial_pn_per_um=1000.0)

    assert topology.seg_filament_id.size == topology.n_links
    assert set(np.unique(topology.seg_polarity).tolist()) <= {-1, 1}  # never 0 (the kernel treats 0 as +1)
    # arclength runs monotonically within a filament and restarts at each new one.
    for fid in np.unique(topology.seg_filament_id):
        mask = topology.seg_filament_id == fid
        s0, s1 = topology.seg_s0[mask], topology.seg_s1[mask]
        assert s0[0] == pytest.approx(0.0)
        assert np.all(np.diff(s0) > 0.0)
        assert np.allclose(s1 - s0, topology.link_r0[mask])
    # every segment id belongs to the sf_arc population's own ID block (no cortex id leaks in).
    lo, hi = population.ledger.block
    assert int(topology.seg_filament_id.min()) >= lo
    assert int(topology.seg_filament_id.max()) < hi


def test_required_params_have_no_defaults() -> None:
    """REQUIRED-PARAM discipline: the head config refuses to construct with an unset GAP."""
    with pytest.raises(ValueError, match="REQUIRED-PARAM"):
        NMIIMotorParams(
            v0=0.12, f_stall=0.5, kappa=0.5, k_xb=1000.0, r0_head=0.2, r0_xb=0.0,
            capture_radius=0.21, k_on=None, k_off0=0.35, f0=7.13,
            k_catch0=0.35, x_catch=1e-3, k_slip0=0.35, x_slip=6e-4)
