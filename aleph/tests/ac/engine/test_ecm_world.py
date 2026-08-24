"""Structural/runtime gates for the collagen ECM World seam."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import pytest
import warp as wp

from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorContract,
    ConnectorFamily,
    ConnectorScope,
    reference_cell_architecture,
)
from aleph.engine.ecm_world import (
    ECM_CROSSLINK_CONNECTOR,
    FA_SERIES_GROUP,
    MEMBRANE_ECM_CONTACT,
    BoundaryAnchorMode,
    CompositeLoadPathECMClutchAdapter,
    ECMBoundaryAnchorFacade,
    ECMStateOwner,
    ECMWorld,
    ECMWorldSettings,
    ECMWorldStepBindings,
    FarFieldDirichletBoundaryRuntime,
    LegacyLocalMikadoMechanicsAdapter,
    MembraneContactSurfaceView,
)
from aleph.engine.load_path import (
    ActorRecord,
    BorrowedSegmentActorView,
    ElementKind,
    EndpointRole,
)
from aleph.engine.ecm_world import BoundaryAnchorDelegate
from aleph.engine.runtime import LedgerContributor, TransactionParticipant


def _cuda_available() -> bool:
    """Return whether a real Warp CUDA device is present (GPU lane only)."""
    try:
        wp.init()
        return bool(wp.get_device().is_cuda)
    except Exception:  # pragma: no cover - depends on hardware
        return False


def _dirichlet_runtime(
    *,
    ptr_base: int = 700,
    device: _FakeDevice | None = None,
    n_pinned: int = 4,
) -> FarFieldDirichletBoundaryRuntime:
    if device is None:
        device = _FakeDevice()
    return FarFieldDirichletBoundaryRuntime(
        pinned_index_d=_FakeArray(ptr_base, (n_pinned,), wp.int32, device),
        target_position_d=_FakeArray(ptr_base + 1, (n_pinned,), wp.vec3d, device),
        reaction_d=_FakeArray(ptr_base + 2, (n_pinned,), wp.vec3d, device),
        committed_reaction_d=_FakeArray(ptr_base + 3, (n_pinned,), wp.vec3d, device),
        committed_target_d=_FakeArray(ptr_base + 4, (n_pinned,), wp.vec3d, device),
        boundary_work_d=_FakeArray(ptr_base + 5, (1,), wp.float64, device),
    )


@dataclass(slots=True)
class _BoundaryLedgerSink:
    calls: list[tuple[object, object]] = field(default_factory=list)

    def add_far_field_reaction(self, reaction: object, work: object) -> None:
        self.calls.append((reaction, work))


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:0"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """Metadata-only CUDA array double; authoritative physics never executes here."""

    ptr: int
    shape: tuple[int, ...]
    dtype: object
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _MechanicsSpy:
    order: list[str]
    n_nodes: int = 8
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.order.append("ecm.internal")
        self.calls.append((pos, force))


@dataclass(slots=True)
class _TransactionSpy:
    calls: list[tuple[object, ...]] = field(default_factory=list)

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))


@dataclass(slots=True)
class _LedgerSpy:
    calls: list[object] = field(default_factory=list)

    def accumulate_ledger(self, ledger: object) -> None:
        self.calls.append(ledger)


@dataclass(slots=True)
class _CrosslinkSpy:
    order: list[str]
    name: str = ECM_CROSSLINK_CONNECTOR
    component_a: str = "ecm"
    component_b: str = "ecm"
    adjoint_transfer_required: bool = True
    mechanics_calls: list[tuple[object, object]] = field(default_factory=list)
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_internal(self, pos: wp.array, force: wp.array) -> None:
        self.order.append("ecm.crosslink")
        self.mechanics_calls.append((pos, force))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


@dataclass(slots=True)
class _BoundaryDelegateSpy:
    order: list[str]
    mechanics_calls: list[tuple[object, object]] = field(default_factory=list)
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_boundary(self, pos: wp.array, force: wp.array) -> None:
        self.order.append("ecm.boundary")
        self.mechanics_calls.append((pos, force))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


@dataclass(slots=True)
class _ContactSpy:
    order: list[str]
    name: str = MEMBRANE_ECM_CONTACT
    component_a: str = "membrane"
    component_b: str = "ecm"
    adjoint_transfer_required: bool = True
    mechanics_calls: list[tuple[object, object]] = field(default_factory=list)
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate_contact(self, membrane: object, ecm: object) -> None:
        self.order.append("membrane.ecm_contact")
        self.mechanics_calls.append((membrane, ecm))

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


@dataclass(slots=True)
class _ClutchSpy:
    order: list[str]
    ecm_endpoint: object
    mechanical_group: str = FA_SERIES_GROUP
    equal_opposite_required: bool = True
    adjoint_transfer_required: bool = True
    work_ledger_required: bool = True
    mechanics_calls: int = 0
    transaction_calls: list[tuple[object, ...]] = field(default_factory=list)
    ledger_calls: list[object] = field(default_factory=list)

    def accumulate(self) -> None:
        self.order.append("fa.composite")
        self.mechanics_calls += 1

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: wp.array) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


@dataclass(slots=True)
class _JointSpy:
    actor_a: object
    actor_b: object
    calls: list[tuple[object, ...]] = field(default_factory=list)

    def accumulate(self) -> None:
        self.calls.append(("accumulate",))

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted: object) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: object, dt_phys: float, rng_seed: int) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))


def _owner(
    order: list[str],
    *,
    ptr_base: int = 100,
    device: _FakeDevice | None = None,
    mechanics_nodes: int = 8,
) -> tuple[ECMStateOwner, _MechanicsSpy, _TransactionSpy, _LedgerSpy]:
    if device is None:
        device = _FakeDevice()
    mechanics = _MechanicsSpy(order, n_nodes=mechanics_nodes)
    transaction = _TransactionSpy()
    ledger = _LedgerSpy()
    owner = ECMStateOwner(
        name="ecm",
        settings=ECMWorldSettings(BoundaryAnchorMode.FAR_FIELD_DIRICHLET),
        actor=ActorRecord("ecm", 7, 2, 13, 3, 6),
        position_d=_FakeArray(ptr_base, (8,), wp.vec3d, device),
        force_d=_FakeArray(ptr_base + 1, (8,), wp.vec3d, device),
        segments_d=_FakeArray(ptr_base + 2, (6, 2), wp.int32, device),
        fiber_offsets_d=_FakeArray(ptr_base + 3, (3,), wp.int32, device),
        fiber_sleep_state_d=_FakeArray(ptr_base + 4, (2,), wp.int32, device),
        segment_refinement_level_d=_FakeArray(ptr_base + 5, (6,), wp.int32, device),
        segment_damage_d=_FakeArray(ptr_base + 6, (6,), wp.float64, device),
        endpoint_refcount_d=_FakeArray(ptr_base + 7, (6,), wp.int32, device),
        topology_epoch_d=_FakeArray(ptr_base + 8, (1,), wp.int32, device),
        mechanics=mechanics,
        transaction=transaction,
        ledger=ledger,
    )
    return owner, mechanics, transaction, ledger


def _world_and_bindings() -> tuple[
    ECMWorld,
    ECMWorldStepBindings,
    list[str],
    tuple[_MechanicsSpy, _TransactionSpy, _LedgerSpy],
    _CrosslinkSpy,
    _BoundaryDelegateSpy,
    _ContactSpy,
    _ClutchSpy,
]:
    order: list[str] = []
    owner, mechanics, transaction, ledger = _owner(order)
    crosslink = _CrosslinkSpy(order)
    boundary_delegate = _BoundaryDelegateSpy(order)
    boundary = ECMBoundaryAnchorFacade(
        boundary_delegate,
        BoundaryAnchorMode.FAR_FIELD_DIRICHLET,
    )
    clutch = _ClutchSpy(order, owner.collagen_endpoint_view())
    contact = _ContactSpy(order)
    membrane = MembraneContactSurfaceView(
        "membrane",
        _FakeArray(400, (5,), wp.vec3d),
        _FakeArray(401, (5,), wp.vec3d),
    )
    bindings = ECMWorldStepBindings(clutch, crosslink, boundary, membrane, contact)
    return (
        ECMWorld(reference_cell_architecture(), owner),
        bindings,
        order,
        (mechanics, transaction, ledger),
        crosslink,
        boundary_delegate,
        contact,
        clutch,
    )


def _architecture_replacing_connector(connector: ConnectorContract) -> CellArchitecture:
    base = reference_cell_architecture()
    connectors = tuple(
        connector if item.name == connector.name else item for item in base.connectors
    )
    return CellArchitecture(base.components, connectors)


def test_reference_graph_and_owner_register_a_live_collagen_environment() -> None:
    world, bindings, *_ = _world_and_bindings()
    ecm_contract = world.architecture.component("ecm")
    crosslink_contract = next(
        item for item in world.architecture.connectors if item.name == ECM_CROSSLINK_CONNECTOR
    )

    assert ecm_contract.owns_geometry
    assert crosslink_contract.family is ConnectorFamily.FIBER_CROSSLINK
    assert crosslink_contract.scope is ConnectorScope.INTERNAL
    assert crosslink_contract.chemistry_card == "collagen_crosslink"
    assert bindings.alpha2beta1_clutch.mechanical_group == FA_SERIES_GROUP
    assert isinstance(world.ecm, TransactionParticipant)
    assert isinstance(world.ecm, LedgerContributor)


def test_load_path_joint_adapter_borrows_live_ecm_and_closes_composite_api() -> None:
    owner, *_ = _owner([])
    sf_record = ActorRecord("sf_arc", 3, 1, 9, 1, 1)
    sf = BorrowedSegmentActorView(
        sf_record,
        _FakeArray(900, (2,), wp.vec3d),
        _FakeArray(901, (2,), wp.vec3d),
        _FakeArray(902, (1, 2), wp.int32),
    )
    ecm = BorrowedSegmentActorView(
        owner.actor,
        owner.position_d,
        owner.force_d,
        owner.segments_d,
    )
    joint = _JointSpy(sf, ecm)
    ledger_spy = _LedgerSpy()
    adapter = CompositeLoadPathECMClutchAdapter(
        joint, owner.collagen_endpoint_view(), ledger_spy
    )
    accepted = object()
    ledger = object()
    adapter.accumulate()
    adapter.snapshot_candidate()
    adapter.rollback(accepted)
    adapter.commit_irreversible(accepted, 0.05, 4)
    adapter.accumulate_ledger(ledger)

    assert adapter.mechanical_group == FA_SERIES_GROUP
    assert joint.calls == [
        ("accumulate",),
        ("snapshot",),
        ("rollback", accepted),
        ("commit", accepted, 0.05, 4),
    ]
    assert ledger_spy.calls == [ledger]


def test_state_owner_requires_cuda_disjoint_shape_consistent_arrays() -> None:
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        _owner([], device=cpu)

    owner, mechanics, transaction, ledger = _owner([])
    with pytest.raises(ValueError, match="distinct storage"):
        ECMStateOwner(
            name="ecm",
            settings=owner.settings,
            actor=owner.actor,
            position_d=owner.position_d,
            force_d=replace(owner.force_d, ptr=owner.position_d.ptr),
            segments_d=owner.segments_d,
            fiber_offsets_d=owner.fiber_offsets_d,
            fiber_sleep_state_d=owner.fiber_sleep_state_d,
            segment_refinement_level_d=owner.segment_refinement_level_d,
            segment_damage_d=owner.segment_damage_d,
            endpoint_refcount_d=owner.endpoint_refcount_d,
            topology_epoch_d=owner.topology_epoch_d,
            mechanics=mechanics,
            transaction=transaction,
            ledger=ledger,
        )

    with pytest.raises(ValueError, match="one entry per collagen segment"):
        ECMStateOwner(
            name="ecm",
            settings=owner.settings,
            actor=owner.actor,
            position_d=owner.position_d,
            force_d=owner.force_d,
            segments_d=owner.segments_d,
            fiber_offsets_d=owner.fiber_offsets_d,
            fiber_sleep_state_d=owner.fiber_sleep_state_d,
            segment_refinement_level_d=replace(
                owner.segment_refinement_level_d,
                shape=(5,),
            ),
            segment_damage_d=owner.segment_damage_d,
            endpoint_refcount_d=owner.endpoint_refcount_d,
            topology_epoch_d=owner.topology_epoch_d,
            mechanics=mechanics,
            transaction=transaction,
            ledger=ledger,
        )


def test_state_owner_checks_actor_and_mechanics_population_metadata() -> None:
    with pytest.raises(ValueError, match="actor element count"):
        owner, mechanics, transaction, ledger = _owner([])
        ECMStateOwner(
            name="ecm",
            settings=owner.settings,
            actor=replace(owner.actor, n_elements=5),
            position_d=owner.position_d,
            force_d=owner.force_d,
            segments_d=owner.segments_d,
            fiber_offsets_d=owner.fiber_offsets_d,
            fiber_sleep_state_d=owner.fiber_sleep_state_d,
            segment_refinement_level_d=owner.segment_refinement_level_d,
            segment_damage_d=owner.segment_damage_d,
            endpoint_refcount_d=owner.endpoint_refcount_d,
            topology_epoch_d=owner.topology_epoch_d,
            mechanics=mechanics,
            transaction=transaction,
            ledger=ledger,
        )

    with pytest.raises(ValueError, match="expects 7 nodes"):
        _owner([], mechanics_nodes=7)


def test_collagen_endpoint_is_a_nonowning_generation_aware_segment_material_point() -> None:
    owner, *_ = _owner([])
    endpoint = owner.collagen_endpoint_view()
    port = endpoint.port_ref(4, 0.25)

    assert endpoint.position_d is owner.position_d
    assert endpoint.force_d is owner.force_d
    assert endpoint.segments_d is owner.segments_d
    assert endpoint.topology_epoch_d is owner.topology_epoch_d
    assert endpoint.endpoint_refcount_d is owner.endpoint_refcount_d
    assert port.component == "ecm"
    assert port.element_kind is ElementKind.SEGMENT
    assert port.element_id == 4
    assert port.segment_u == pytest.approx(0.25)
    assert port.role is EndpointRole.COLLAGEN_LIGAND
    assert port.actor_generation == owner.actor.actor_generation
    assert port.entity_generation == owner.actor.entity_generation

    with pytest.raises(ValueError, match="outside"):
        endpoint.port_ref(6, 0.5)
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        endpoint.port_ref(0, 1.1)


def test_mechanics_wiring_visits_each_contribution_once_in_order() -> None:
    """The far-field boundary must close LAST, after every force contributor.

    This assertion previously encoded the opposite order (boundary third), which is the defect
    itself: ``_far_field_reaction_ledger_kernel`` records ``reaction = -f`` and then zeroes the
    pinned node, so running it before ``membrane_contact`` and ``alpha2beta1_clutch`` both excluded
    those two channels from the recorded traction AND left the pinned nodes carrying a net force
    from them, silently violating the Dirichlet pin.
    """
    world, bindings, order, owner_spies, crosslink, boundary, contact, clutch = (
        _world_and_bindings()
    )
    world.accumulate_mechanics(bindings)

    assert order == [
        "ecm.internal",
        "ecm.crosslink",
        "membrane.ecm_contact",
        "fa.composite",
        "ecm.boundary",
    ]
    assert order[-1] == "ecm.boundary", "the reaction ledger must see every contribution"
    assert owner_spies[0].calls == [(world.ecm.position_d, world.ecm.force_d)]
    assert crosslink.mechanics_calls == [(world.ecm.position_d, world.ecm.force_d)]
    assert boundary.mechanics_calls == [(world.ecm.position_d, world.ecm.force_d)]
    membrane, endpoint = contact.mechanics_calls[0]
    assert membrane is bindings.membrane
    assert endpoint.owner is world.ecm
    assert clutch.mechanics_calls == 1


def test_bindings_require_same_endpoint_owner_and_physiological_anchor_mode() -> None:
    world, bindings, *_ = _world_and_bindings()
    foreign_owner, *_ = _owner([], ptr_base=500)
    foreign_clutch = _ClutchSpy([], foreign_owner.collagen_endpoint_view())
    foreign_bindings = replace(bindings, alpha2beta1_clutch=foreign_clutch)
    with pytest.raises(ValueError, match="this ECM state owner"):
        world.accumulate_mechanics(foreign_bindings)

    boundary_delegate = _BoundaryDelegateSpy([])
    impedance = ECMBoundaryAnchorFacade(
        boundary_delegate,
        BoundaryAnchorMode.CONTINUUM_IMPEDANCE,
    )
    with pytest.raises(ValueError, match="match the physiological"):
        world.accumulate_mechanics(replace(bindings, boundary_anchor=impedance))


def test_binding_metadata_enforces_crosslink_and_clutch_conservation_contracts() -> None:
    world, bindings, *_ = _world_and_bindings()
    with pytest.raises(ValueError, match="equal/opposite"):
        ECMWorldStepBindings(
            replace(bindings.alpha2beta1_clutch, equal_opposite_required=False),
            bindings.crosslink_connector,
            bindings.boundary_anchor,
            bindings.membrane,
            bindings.membrane_contact,
        )
    with pytest.raises(ValueError, match=ECM_CROSSLINK_CONNECTOR):
        ECMWorldStepBindings(
            bindings.alpha2beta1_clutch,
            replace(bindings.crosslink_connector, name="private_link"),
            bindings.boundary_anchor,
            bindings.membrane,
            bindings.membrane_contact,
        )
    with pytest.raises(ValueError, match=MEMBRANE_ECM_CONTACT):
        replace(bindings, membrane_contact=replace(bindings.membrane_contact, name="private_contact"))
    assert world.ecm.settings.anchor_mode is BoundaryAnchorMode.FAR_FIELD_DIRICHLET


def test_one_transaction_covers_remodelling_crosslinks_boundary_contact_and_clutch() -> None:
    world, bindings, _, owner_spies, crosslink, boundary, contact, clutch = (
        _world_and_bindings()
    )
    accepted_d = object()

    world.snapshot_candidate(bindings)
    world.rollback(bindings, accepted_d)
    world.commit_irreversible(bindings, accepted_d, dt_phys=0.05, rng_seed=19)

    expected = [
        ("snapshot",),
        ("rollback", accepted_d),
        ("commit", accepted_d, 0.05, 19),
    ]
    assert owner_spies[1].calls == expected
    assert crosslink.transaction_calls == expected
    assert boundary.transaction_calls == expected
    assert contact.transaction_calls == expected
    assert clutch.transaction_calls == expected


@pytest.mark.parametrize(
    ("dt_phys", "rng_seed"),
    [(0.0, 1), (-0.1, 1), (float("nan"), 1), (float("inf"), 1), (0.1, -1)],
)
def test_invalid_commit_inputs_fail_before_any_participant_mutates(
    dt_phys: float,
    rng_seed: int,
) -> None:
    world, bindings, _, owner_spies, crosslink, boundary, contact, clutch = (
        _world_and_bindings()
    )
    with pytest.raises(ValueError):
        world.commit_irreversible(bindings, object(), dt_phys, rng_seed)
    assert owner_spies[1].calls == []
    assert crosslink.transaction_calls == []
    assert boundary.transaction_calls == []
    assert contact.transaction_calls == []
    assert clutch.transaction_calls == []


def test_ledgers_include_component_connectors_boundary_reaction_and_allocation_metadata() -> None:
    world, bindings, _, owner_spies, crosslink, boundary, contact, clutch = (
        _world_and_bindings()
    )
    ledger = object()
    world.accumulate_ledger(bindings, ledger)

    assert owner_spies[2].calls == [ledger]
    assert crosslink.ledger_calls == [ledger]
    assert boundary.ledger_calls == [ledger]
    assert contact.ledger_calls == [ledger]
    assert clutch.ledger_calls == [ledger]
    assert world.dof_ledger() == replace(
        world.dof_ledger(),
        allocated_nodes=8,
        allocated_segments=6,
        allocated_fibers=2,
        sleep_state_slots=2,
        refinement_state_slots=6,
        damage_state_slots=6,
        endpoint_refcount_slots=6,
        topology_epoch_slots=1,
    )


def test_settings_and_boundary_reject_nonphysiological_or_unledgered_builds() -> None:
    with pytest.raises(ValueError, match="physiological far-field anchor"):
        ECMWorldSettings(
            BoundaryAnchorMode.FAR_FIELD_DIRICHLET,
            physiological_anchor_required=False,
        )
    with pytest.raises(ValueError, match="remodelling and damage"):
        ECMWorldSettings(
            BoundaryAnchorMode.FAR_FIELD_DIRICHLET,
            damage_enabled=False,
        )
    with pytest.raises(ValueError, match="reaction and boundary work"):
        ECMBoundaryAnchorFacade(
            _BoundaryDelegateSpy([]),
            BoundaryAnchorMode.FAR_FIELD_DIRICHLET,
            reaction_ledger_required=False,
        )


def test_architecture_rejects_wrong_collagen_chemistry_and_nonaccepted_kinetics() -> None:
    owner, *_ = _owner([])
    wrong_chemistry = ConnectorContract(
        ECM_CROSSLINK_CONNECTOR,
        ConnectorFamily.FIBER_CROSSLINK,
        "ecm",
        "ecm",
        kinetics=True,
        commit_on_accept=True,
        scope=ConnectorScope.INTERNAL,
        chemistry_card="actin_crosslink",
        generation_required=True,
        remap_on_accept=True,
        blocks_sleep_refine=True,
    )
    with pytest.raises(ValueError, match="collagen_crosslink chemistry"):
        ECMWorld(_architecture_replacing_connector(wrong_chemistry), owner)

    wrong_clutch = ConnectorContract(
        "integrin_collagen_clutch",
        ConnectorFamily.FA_CLUTCH,
        "focal_adhesion",
        "ecm",
        kinetics=True,
        commit_on_accept=True,
        chemistry_card="alpha2beta1_collagen",
        mechanical_group=FA_SERIES_GROUP,
        adjoint_transfer_required=False,
    )
    with pytest.raises(ValueError, match="bidirectional with adjoint"):
        ECMWorld(_architecture_replacing_connector(wrong_clutch), owner)


@dataclass(slots=True)
class _LegacyMikadoSpy:
    node_off: int = 0
    n_nodes: int = 8
    owns_crosslinks: bool = False
    owns_boundary_anchor: bool = False
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        self.calls.append((pos, force))


def test_legacy_adapter_reuses_only_local_connector_free_boundary_free_mechanics() -> None:
    legacy = _LegacyMikadoSpy()
    adapter = LegacyLocalMikadoMechanicsAdapter(legacy)
    pos = _FakeArray(1, (8,), wp.vec3d)
    force = _FakeArray(2, (8,), wp.vec3d)
    adapter.accumulate(pos, force)
    assert adapter.n_nodes == 8
    assert legacy.calls == [(pos, force)]

    with pytest.raises(ValueError, match="node_off == 0"):
        LegacyLocalMikadoMechanicsAdapter(_LegacyMikadoSpy(node_off=3))
    with pytest.raises(ValueError, match="cannot own graph crosslink"):
        LegacyLocalMikadoMechanicsAdapter(_LegacyMikadoSpy(owns_crosslinks=True))
    with pytest.raises(ValueError, match="cannot hide"):
        LegacyLocalMikadoMechanicsAdapter(_LegacyMikadoSpy(owns_boundary_anchor=True))


def test_far_field_runtime_owns_disjoint_cuda_ledger_and_satisfies_the_boundary_seam() -> None:
    runtime = _dirichlet_runtime()

    assert runtime.n_pinned == 4
    assert runtime.device == "cuda:0"
    # KERNEL_BOUND: the concrete runtime is the real boundary delegate and transacts/ledgers.
    assert isinstance(runtime, BoundaryAnchorDelegate)
    assert isinstance(runtime, TransactionParticipant)
    assert isinstance(runtime, LedgerContributor)

    # It plugs straight into the existing far-field anchor facade seam.
    facade = ECMBoundaryAnchorFacade(runtime, BoundaryAnchorMode.FAR_FIELD_DIRICHLET)
    assert facade.delegate is runtime


def test_far_field_runtime_rejects_host_free_and_aliased_boundary_state() -> None:
    cpu = _FakeDevice(alias="cpu", is_cuda=False)
    with pytest.raises(ValueError, match="CUDA device array"):
        _dirichlet_runtime(device=cpu)

    # An unanchored (empty) far-field set is non-physiological and must be rejected.
    with pytest.raises(ValueError, match="no empty dimension"):
        _dirichlet_runtime(n_pinned=0)

    base = _dirichlet_runtime()
    with pytest.raises(ValueError, match="distinct storage"):
        replace(base, committed_reaction_d=replace(base.reaction_d))
    with pytest.raises(ValueError, match="one entry per pinned"):
        replace(base, reaction_d=replace(base.reaction_d, shape=(3,)))
    with pytest.raises(ValueError, match="one accepted work scalar"):
        replace(base, boundary_work_d=replace(base.boundary_work_d, shape=(2,)))


def test_far_field_runtime_forwards_reaction_ledger_without_a_device_read() -> None:
    runtime = _dirichlet_runtime()
    sink = _BoundaryLedgerSink()
    facade = ECMBoundaryAnchorFacade(runtime, BoundaryAnchorMode.FAR_FIELD_DIRICHLET)

    facade.accumulate_ledger(sink)

    assert sink.calls == [(runtime.committed_reaction_d, runtime.boundary_work_d)]
    # A ledger sink without the boundary hook is tolerated (schema gap), never a host readback.
    runtime.accumulate_ledger(object())


def test_far_field_runtime_dispatches_to_the_concrete_global_ledger_type() -> None:
    """The typed GlobalCellLedger is the recognised sink; a recording double falls to the duck path.

    The concrete reduction launch is a CUDA_UNIT gate (test_ledger.py); here we prove — CPU-green — that the
    boundary's forwarding branches on the concrete typed ledger, so the ``add_far_field_reaction`` slot is no
    longer a duck-typed gap.
    """
    from aleph.engine.ledger import GlobalCellLedger

    device = _FakeDevice()
    ledger = GlobalCellLedger(
        reaction_resultant_d=_FakeArray(500, (1,), wp.vec3d, device),
        traction_resultant_d=_FakeArray(501, (1,), wp.vec3d, device),
        force_resultant_d=_FakeArray(502, (1,), wp.vec3d, device),
        work_d=_FakeArray(503, (1,), wp.float64, device),
        mass_d=_FakeArray(504, (1,), wp.float64, device),
        topology_counts_d=_FakeArray(505, (1,), wp.int32, device),
        balance_residual_sq_d=_FakeArray(506, (1,), wp.float64, device),
        balance_ok_d=_FakeArray(507, (1,), wp.int32, device),
    )
    assert isinstance(ledger, GlobalCellLedger)          # boundary takes the typed branch for this sink
    assert not isinstance(_BoundaryLedgerSink(), GlobalCellLedger)  # a recorder uses the tolerant duck path
    assert callable(ledger.add_far_field_reaction)


@pytest.mark.parametrize(
    ("dt_phys", "rng_seed"),
    [(0.0, 1), (-1.0, 1), (float("nan"), 1), (float("inf"), 1), (0.1, -1)],
)
def test_far_field_commit_rejects_invalid_inputs_before_any_launch(
    dt_phys: float,
    rng_seed: int,
) -> None:
    runtime = _dirichlet_runtime()
    # These raise in Python before any wp.launch, so the guard is CPU-green off-CUDA.
    with pytest.raises(ValueError):
        runtime.commit_irreversible(object(), dt_phys, rng_seed)


@pytest.mark.skipif(not _cuda_available(), reason="CUDA_UNIT reaction gate runs on the GPU lane")
def test_far_field_reaction_kernel_preserves_then_removes_pinned_force() -> None:  # pragma: no cover
    import numpy as np

    device = "cuda:0"
    force = wp.array(
        np.array([[1.0, -2.0, 3.0], [0.0, 0.0, 0.0], [-4.0, 5.0, -6.0]], dtype=np.float64),
        dtype=wp.vec3d,
        device=device,
    )
    runtime = FarFieldDirichletBoundaryRuntime(
        pinned_index_d=wp.array(np.array([0, 2], dtype=np.int32), dtype=wp.int32, device=device),
        target_position_d=wp.zeros(2, dtype=wp.vec3d, device=device),
        reaction_d=wp.zeros(2, dtype=wp.vec3d, device=device),
        committed_reaction_d=wp.zeros(2, dtype=wp.vec3d, device=device),
        committed_target_d=wp.zeros(2, dtype=wp.vec3d, device=device),
        boundary_work_d=wp.zeros(1, dtype=wp.float64, device=device),
    )

    runtime.accumulate_boundary(force, force)

    reaction = runtime.reaction_d.numpy()
    residual = force.numpy()
    # Reaction is the removed constraint reaction R = -f (sign gate).
    assert np.allclose(reaction[0], [-1.0, 2.0, -3.0])
    assert np.allclose(reaction[1], [4.0, -5.0, 6.0])
    # The pinned nodes are zeroed; the free node is untouched (boundary gate).
    assert np.allclose(residual[0], [0.0, 0.0, 0.0])
    assert np.allclose(residual[2], [0.0, 0.0, 0.0])
    assert np.allclose(residual[1], [0.0, 0.0, 0.0])
