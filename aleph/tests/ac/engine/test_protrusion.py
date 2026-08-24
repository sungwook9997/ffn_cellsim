"""Structural gates for explicit lamellipodium and filopodium actors."""

from __future__ import annotations

from dataclasses import dataclass, replace

import pytest
import warp as wp

from aleph.engine.contracts import ComponentRole, ConnectorFamily, reference_cell_architecture
from aleph.engine.protrusion import (
    FASCIN_INTERFILAMENT_SPACING_UM,
    FASCIN_K_CROSSLINK_PN_PER_UM,
    FILO_CORTEX,
    FILO_CYTOSOL,
    FILO_MEMBRANE,
    FILO_NASCENT_FA,
    FILOPODIUM_COMPONENT,
    FILOPODIUM_EVENT_CHANNELS,
    FILOPODIUM_REPRESENTATION,
    FilopodiumFascinBundleMechanics,
    fascin_bundle_pair_reference,
    LAM_CORTEX,
    LAM_CYTOSOL,
    LAM_MEMBRANE,
    LAM_NASCENT_FA,
    LAMELLIPODIUM_COMPONENT,
    LAMELLIPODIUM_EVENT_CHANNELS,
    LAMELLIPODIUM_REPRESENTATION,
    FilopodiumStateOwner,
    LAM_BRANCH_K_THETA_PN_UM,
    LAM_BRANCH_THETA0_RAD,
    LamellipodiumBranchAngleMechanics,
    LamellipodiumStateOwner,
    ProtrusionActors,
    ProtrusionStepBindings,
    RefinementClosureReport,
    proposed_protrusion_architecture,
    protrusion_contract_proposal,
)
from aleph.components.weave.branch_angle_warp import branch_angle_kernel


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda:test"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    ptr: int
    shape: tuple[int, ...]
    dtype: object
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _MechanicsSpy:
    calls: list[tuple[object, object]]

    def accumulate(self, pos: object, force: object) -> None:
        self.calls.append((pos, force))


@dataclass(slots=True)
class _LaunchSpy:
    """CPU capturing double for ``wp.launch`` (no CUDA on the dev Mac)."""

    launches: list[dict[str, object]]

    def __call__(
        self,
        kernel: object,
        *,
        dim: int,
        inputs: list[object],
        outputs: list[object],
        device: object,
    ) -> None:
        self.launches.append(
            {
                "kernel": kernel,
                "dim": dim,
                "inputs": tuple(inputs),
                "outputs": tuple(outputs),
                "device": device,
            }
        )


@dataclass(slots=True)
class _TransactionSpy:
    component_name: str
    event_channels: frozenset[str]
    arrays: tuple[object, ...]
    calls: list[tuple[object, ...]]

    def owned_arrays(self) -> tuple[object, ...]:
        return self.arrays

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted: object) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: object, dt_phys: float, rng_seed: int) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))


@dataclass(slots=True)
class _LedgerSpy:
    calls: list[object]

    def accumulate_ledger(self, ledger: object) -> None:
        self.calls.append(ledger)


@dataclass(slots=True)
class _ConnectorSpy:
    name: str
    component_a: str
    component_b: str
    mechanics_calls: list[str]
    transaction_calls: list[tuple[object, ...]]
    ledger_calls: list[object]

    def accumulate_actor(self, actor: object) -> None:
        self.mechanics_calls.append(actor.component)

    def snapshot_candidate(self) -> None:
        self.transaction_calls.append(("snapshot",))

    def rollback(self, accepted: object) -> None:
        self.transaction_calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: object, dt_phys: float, rng_seed: int) -> None:
        self.transaction_calls.append(("commit", accepted, dt_phys, rng_seed))

    def accumulate_ledger(self, ledger: object) -> None:
        self.ledger_calls.append(ledger)


def _arrays(start: int, *, n_node: int = 4, n_fil: int = 2) -> dict[str, _FakeArray]:
    """Create disjoint fake CUDA storage for common protrusion ownership metadata."""
    return {
        "position": _FakeArray(start, (n_node,), wp.vec3d),
        "force": _FakeArray(start + 1, (n_node,), wp.vec3d),
        "offset": _FakeArray(start + 2, (n_fil + 1,), wp.int32),
        "active_filament": _FakeArray(start + 3, (n_fil,), wp.int32),
        "active_node": _FakeArray(start + 4, (n_node,), wp.int32),
        "persistent_id": _FakeArray(start + 5, (n_fil,), wp.int64),
        "parent_id": _FakeArray(start + 6, (n_fil,), wp.int64),
        "level": _FakeArray(start + 7, (n_fil,), wp.int32),
        "generation": _FakeArray(start + 8, (1,), wp.int32),
        "epoch": _FakeArray(start + 9, (1,), wp.int32),
    }


def _lamellipodium() -> tuple[LamellipodiumStateOwner, _MechanicsSpy, _TransactionSpy, _LedgerSpy]:
    a = _arrays(100)
    branch_triples = _FakeArray(110, (1, 3), wp.int32)
    branch_active = _FakeArray(111, (1,), wp.int32)
    barbed = _FakeArray(112, (2,), wp.int32)
    mutable = (
        a["position"],
        a["offset"],
        a["active_filament"],
        a["active_node"],
        a["persistent_id"],
        a["parent_id"],
        a["level"],
        branch_triples,
        branch_active,
        barbed,
        a["generation"],
        a["epoch"],
    )
    mechanics = _MechanicsSpy([])
    transaction = _TransactionSpy(
        LAMELLIPODIUM_COMPONENT,
        LAMELLIPODIUM_EVENT_CHANNELS,
        mutable,
        [],
    )
    ledger = _LedgerSpy([])
    owner = LamellipodiumStateOwner(
        n_filament_capacity=2,
        n_branch_capacity=1,
        position_d=a["position"],
        force_d=a["force"],
        fiber_offset_d=a["offset"],
        active_filament_d=a["active_filament"],
        active_node_d=a["active_node"],
        persistent_filament_id_d=a["persistent_id"],
        refinement_parent_id_d=a["parent_id"],
        refinement_level_d=a["level"],
        branch_triples_d=branch_triples,
        branch_active_d=branch_active,
        barbed_state_d=barbed,
        actor_generation_d=a["generation"],
        topology_epoch_d=a["epoch"],
        mechanics=mechanics,
        transaction=transaction,
        ledger=ledger,
    )
    return owner, mechanics, transaction, ledger


def _lamellipodium_kernel_bound() -> tuple[
    LamellipodiumStateOwner,
    LamellipodiumBranchAngleMechanics,
    _LaunchSpy,
    _FakeArray,
    _FakeArray,
]:
    """Build a lamellipodium whose mechanics launches the real branch_angle_kernel via a spy launcher."""
    a = _arrays(300)
    branch_triples = _FakeArray(310, (1, 3), wp.int32)
    branch_active = _FakeArray(311, (1,), wp.int32)
    barbed = _FakeArray(312, (2,), wp.int32)
    mutable = (
        a["position"],
        a["offset"],
        a["active_filament"],
        a["active_node"],
        a["persistent_id"],
        a["parent_id"],
        a["level"],
        branch_triples,
        branch_active,
        barbed,
        a["generation"],
        a["epoch"],
    )
    launcher = _LaunchSpy([])
    mechanics = LamellipodiumBranchAngleMechanics(
        n_branch=1,
        branch_triples_d=branch_triples,
        branch_active_d=branch_active,
        launch=launcher,
    )
    transaction = _TransactionSpy(
        LAMELLIPODIUM_COMPONENT,
        LAMELLIPODIUM_EVENT_CHANNELS,
        mutable,
        [],
    )
    owner = LamellipodiumStateOwner(
        n_filament_capacity=2,
        n_branch_capacity=1,
        position_d=a["position"],
        force_d=a["force"],
        fiber_offset_d=a["offset"],
        active_filament_d=a["active_filament"],
        active_node_d=a["active_node"],
        persistent_filament_id_d=a["persistent_id"],
        refinement_parent_id_d=a["parent_id"],
        refinement_level_d=a["level"],
        branch_triples_d=branch_triples,
        branch_active_d=branch_active,
        barbed_state_d=barbed,
        actor_generation_d=a["generation"],
        topology_epoch_d=a["epoch"],
        mechanics=mechanics,
        transaction=transaction,
        ledger=_LedgerSpy([]),
    )
    return owner, mechanics, launcher, branch_triples, branch_active


def _filopodium() -> tuple[FilopodiumStateOwner, _MechanicsSpy, _TransactionSpy, _LedgerSpy]:
    a = _arrays(200)
    bundle = _FakeArray(210, (2,), wp.int32)
    polarity = _FakeArray(211, (2,), wp.int32)
    pairs = _FakeArray(212, (1, 2), wp.int32)
    crosslink_active = _FakeArray(213, (1,), wp.int32)
    barbed = _FakeArray(214, (2,), wp.int32)
    mutable = (
        a["position"],
        a["offset"],
        a["active_filament"],
        a["active_node"],
        a["persistent_id"],
        a["parent_id"],
        a["level"],
        bundle,
        polarity,
        pairs,
        crosslink_active,
        barbed,
        a["generation"],
        a["epoch"],
    )
    mechanics = _MechanicsSpy([])
    transaction = _TransactionSpy(
        FILOPODIUM_COMPONENT,
        FILOPODIUM_EVENT_CHANNELS,
        mutable,
        [],
    )
    ledger = _LedgerSpy([])
    owner = FilopodiumStateOwner(
        n_filament_capacity=2,
        n_crosslink_capacity=1,
        position_d=a["position"],
        force_d=a["force"],
        fiber_offset_d=a["offset"],
        active_filament_d=a["active_filament"],
        active_node_d=a["active_node"],
        persistent_filament_id_d=a["persistent_id"],
        refinement_parent_id_d=a["parent_id"],
        refinement_level_d=a["level"],
        bundle_id_d=bundle,
        polarity_d=polarity,
        crosslink_pairs_d=pairs,
        crosslink_active_d=crosslink_active,
        barbed_state_d=barbed,
        actor_generation_d=a["generation"],
        topology_epoch_d=a["epoch"],
        mechanics=mechanics,
        transaction=transaction,
        ledger=ledger,
    )
    return owner, mechanics, transaction, ledger


def _connector(name: str, a: str, b: str) -> _ConnectorSpy:
    return _ConnectorSpy(name, a, b, [], [], [])


def _bindings() -> ProtrusionStepBindings:
    return ProtrusionStepBindings(
        (
            _connector(LAM_MEMBRANE, LAMELLIPODIUM_COMPONENT, "membrane"),
            _connector(LAM_CORTEX, LAMELLIPODIUM_COMPONENT, "cortex"),
            _connector(LAM_CYTOSOL, LAMELLIPODIUM_COMPONENT, "cytosol"),
            _connector(LAM_NASCENT_FA, LAMELLIPODIUM_COMPONENT, "focal_adhesion"),
            _connector(FILO_MEMBRANE, FILOPODIUM_COMPONENT, "membrane"),
            _connector(FILO_CORTEX, FILOPODIUM_COMPONENT, "cortex"),
            _connector(FILO_CYTOSOL, FILOPODIUM_COMPONENT, "cytosol"),
            _connector(FILO_NASCENT_FA, FILOPODIUM_COMPONENT, "focal_adhesion"),
        )
    )


def _facade() -> tuple[
    ProtrusionActors,
    LamellipodiumStateOwner,
    FilopodiumStateOwner,
    tuple[object, ...],
]:
    lam, lam_mech, lam_tx, lam_ledger = _lamellipodium()
    filo, filo_mech, filo_tx, filo_ledger = _filopodium()
    architecture = proposed_protrusion_architecture(reference_cell_architecture())
    return (
        ProtrusionActors(architecture, lam, filo),
        lam,
        filo,
        (lam_mech, lam_tx, lam_ledger, filo_mech, filo_tx, filo_ledger),
    )


def test_proposal_adds_two_explicit_state_owning_actors_idempotently() -> None:
    base = reference_cell_architecture()
    proposal = protrusion_contract_proposal()
    assert {component.name for component in proposal.components} == {
        LAMELLIPODIUM_COMPONENT,
        FILOPODIUM_COMPONENT,
    }
    assert all(component.role is ComponentRole.ACTIVE_LOAD_PATH for component in proposal.components)
    assert "Arp2/3 branched F-actin" in proposal.components[0].representation
    assert "bundled F-actin" in proposal.components[1].representation
    extended = proposed_protrusion_architecture(base)
    twice = proposed_protrusion_architecture(extended)
    assert twice == extended
    for component in proposal.components:
        assert sum(item.name == component.name for item in extended.components) == 1
    for connector in proposal.connectors:
        assert sum(item.name == connector.name for item in extended.connectors) == 1
    proposed = {connector.name: connector for connector in proposal.connectors}
    assert proposed[LAM_MEMBRANE].family is ConnectorFamily.CONTACT
    assert proposed[FILO_MEMBRANE].family is ConnectorFamily.CONTACT


def test_two_actor_facade_preserves_distinct_geometry_and_topology() -> None:
    facade, lam, filo, _ = _facade()
    assert facade.lamellipodium is lam
    assert facade.filopodium is filo
    assert lam.position_d.ptr != filo.position_d.ptr
    assert lam.geometry().representation == LAMELLIPODIUM_REPRESENTATION
    assert filo.geometry().representation == FILOPODIUM_REPRESENTATION
    assert lam.branch_triples_d.shape == (1, 3)
    assert filo.crosslink_pairs_d.shape == (1, 2)


def test_lumped_active_stress_relabelling_is_rejected() -> None:
    lam, *_ = _lamellipodium()
    with pytest.raises(ValueError, match="explicit Arp2/3"):
        replace(lam, representation="lumped active-stress patch")
    filo, *_ = _filopodium()
    with pytest.raises(ValueError, match="explicit bundled"):
        replace(filo, representation="lumped active-stress cable")


def test_transactions_cover_required_events_and_authoritative_arrays() -> None:
    lam, _, transaction, _ = _lamellipodium()
    missing_event = replace(
        transaction,
        event_channels=transaction.event_channels - {"severing"},
    )
    with pytest.raises(ValueError, match="missing event channels"):
        replace(lam, transaction=missing_event)
    missing_array = replace(transaction, arrays=transaction.arrays[:-1])
    with pytest.raises(ValueError, match="all mutable topology arrays"):
        replace(lam, transaction=missing_array)


def test_external_connector_bindings_require_all_named_edges_and_endpoints() -> None:
    bindings = _bindings()
    assert len(bindings.for_actor(LAMELLIPODIUM_COMPONENT)) == 4
    assert len(bindings.for_actor(FILOPODIUM_COMPONENT)) == 4
    with pytest.raises(ValueError, match="missing"):
        ProtrusionStepBindings(bindings.connectors[:-1])
    bad = replace(bindings.connectors[0], component_b="cytosol")
    with pytest.raises(ValueError, match="incorrect component endpoints"):
        ProtrusionStepBindings((bad, *bindings.connectors[1:]))


def test_mechanics_connector_transaction_and_ledger_forwarding() -> None:
    facade, _, _, spies = _facade()
    lam_mech, lam_tx, lam_ledger, filo_mech, filo_tx, filo_ledger = spies
    bindings = _bindings()
    accepted = object()
    ledger = object()

    facade.accumulate_mechanics(bindings)
    assert len(lam_mech.calls) == len(filo_mech.calls) == 1
    for connector in bindings.connectors:
        expected_component = (
            LAMELLIPODIUM_COMPONENT if connector.name.startswith("lamellipodium")
            else FILOPODIUM_COMPONENT
        )
        assert connector.mechanics_calls == [expected_component]

    facade.snapshot_candidate(bindings)
    facade.rollback(bindings, accepted)
    facade.commit_irreversible(bindings, accepted, dt_phys=0.025, rng_seed=71)
    expected = [("snapshot",), ("rollback", accepted), ("commit", accepted, 0.025, 71)]
    assert lam_tx.calls == filo_tx.calls == expected
    assert all(connector.transaction_calls == expected for connector in bindings.connectors)

    facade.accumulate_ledger(bindings, ledger)
    assert lam_ledger.calls == filo_ledger.calls == [ledger]
    assert all(connector.ledger_calls == [ledger] for connector in bindings.connectors)


def test_refinement_mapping_closes_identity_force_work_and_topology() -> None:
    RefinementClosureReport(0, 0, 2.0e-10, 3.0e-11).assert_closed(
        force_tolerance_pn=1.0e-9,
        work_tolerance_pn_um=1.0e-10,
    )
    with pytest.raises(ValueError, match="persistent filament identity"):
        RefinementClosureReport(1, 0, 0.0, 0.0).assert_closed(
            force_tolerance_pn=0.0,
            work_tolerance_pn_um=0.0,
        )
    with pytest.raises(ValueError, match="active graph topology"):
        RefinementClosureReport(0, 1, 0.0, 0.0).assert_closed(
            force_tolerance_pn=0.0,
            work_tolerance_pn_um=0.0,
        )
    with pytest.raises(ValueError, match="resultant force"):
        RefinementClosureReport(0, 0, 1.0e-3, 0.0).assert_closed(
            force_tolerance_pn=1.0e-4,
            work_tolerance_pn_um=0.0,
        )
    with pytest.raises(ValueError, match="virtual work"):
        RefinementClosureReport(0, 0, 0.0, 1.0e-3).assert_closed(
            force_tolerance_pn=0.0,
            work_tolerance_pn_um=1.0e-4,
        )


def test_lamellipodium_binds_the_real_branch_angle_kernel_through_the_seam() -> None:
    owner, _, launcher, branch_triples, branch_active = _lamellipodium_kernel_bound()

    owner.accumulate()

    assert len(launcher.launches) == 1
    launch = launcher.launches[0]
    # The seam carries the real Warp kernel, not a lumped active-stress proxy.
    assert launch["kernel"] is branch_angle_kernel
    assert launch["dim"] == owner.n_branch_capacity
    pos, triples, active, theta0, k_theta = launch["inputs"]
    # No private state path: the kernel reads the actor-owned geometry/topology arrays...
    assert pos.ptr == owner.position_d.ptr
    assert triples.ptr == branch_triples.ptr
    assert active.ptr == branch_active.ptr
    # ...and scatters force straight into the actor-owned accumulator (no private buffer).
    assert launch["outputs"][0].ptr == owner.force_d.ptr
    # Constants are the sourced Faessler 2020 anchors (no magic numbers).
    assert theta0 == LAM_BRANCH_THETA0_RAD
    assert k_theta == LAM_BRANCH_K_THETA_PN_UM
    assert launch["device"] is owner.position_d.device


def test_branch_mechanics_rejects_off_source_constants_and_a_fake_kernel() -> None:
    triples = _FakeArray(410, (1, 3), wp.int32)
    active = _FakeArray(411, (1,), wp.int32)
    with pytest.raises(ValueError, match="theta0 must be the sourced"):
        LamellipodiumBranchAngleMechanics(
            n_branch=1,
            branch_triples_d=triples,
            branch_active_d=active,
            theta0_rad=1.5,
        )
    with pytest.raises(ValueError, match="k_theta must be the equipartition"):
        LamellipodiumBranchAngleMechanics(
            n_branch=1,
            branch_triples_d=triples,
            branch_active_d=active,
            k_theta_pn_um=0.5,
        )
    with pytest.raises(ValueError, match="real branch_angle_kernel"):
        LamellipodiumBranchAngleMechanics(
            n_branch=1,
            branch_triples_d=triples,
            branch_active_d=active,
            kernel=object(),
        )


def test_lamellipodium_rejects_branch_mechanics_bound_to_foreign_topology() -> None:
    a = _arrays(500)
    branch_triples = _FakeArray(510, (1, 3), wp.int32)
    branch_active = _FakeArray(511, (1,), wp.int32)
    barbed = _FakeArray(512, (2,), wp.int32)
    foreign_triples = _FakeArray(520, (1, 3), wp.int32)
    foreign_active = _FakeArray(521, (1,), wp.int32)
    mutable = (
        a["position"],
        a["offset"],
        a["active_filament"],
        a["active_node"],
        a["persistent_id"],
        a["parent_id"],
        a["level"],
        branch_triples,
        branch_active,
        barbed,
        a["generation"],
        a["epoch"],
    )
    # A delegate that carries a private copy of the topology (not the actor's arrays) is rejected.
    mechanics = LamellipodiumBranchAngleMechanics(
        n_branch=1,
        branch_triples_d=foreign_triples,
        branch_active_d=foreign_active,
        launch=_LaunchSpy([]),
    )
    with pytest.raises(ValueError, match="actor-owned branch_triples_d"):
        LamellipodiumStateOwner(
            n_filament_capacity=2,
            n_branch_capacity=1,
            position_d=a["position"],
            force_d=a["force"],
            fiber_offset_d=a["offset"],
            active_filament_d=a["active_filament"],
            active_node_d=a["active_node"],
            persistent_filament_id_d=a["persistent_id"],
            refinement_parent_id_d=a["parent_id"],
            refinement_level_d=a["level"],
            branch_triples_d=branch_triples,
            branch_active_d=branch_active,
            barbed_state_d=barbed,
            actor_generation_d=a["generation"],
            topology_epoch_d=a["epoch"],
            mechanics=mechanics,
            transaction=_TransactionSpy(
                LAMELLIPODIUM_COMPONENT, LAMELLIPODIUM_EVENT_CHANNELS, mutable, []
            ),
            ledger=_LedgerSpy([]),
        )


def test_actor_facade_accepts_the_idempotently_extended_reference_graph() -> None:
    lam, *_ = _lamellipodium()
    filo, *_ = _filopodium()
    architecture = proposed_protrusion_architecture(reference_cell_architecture())
    facade = ProtrusionActors(architecture, lam, filo)
    assert facade.architecture.component(LAMELLIPODIUM_COMPONENT).owns_geometry
    assert facade.architecture.component(FILOPODIUM_COMPONENT).owns_geometry


# ── Filopodium fascin bundle: PI-GAP stiffness slot bound to a real crosslink kernel ─────────────


def test_fascin_geometry_is_sourced_and_stiffness_is_an_unset_slot() -> None:
    # bundle GEOMETRY is sourced (Courson & Rock 2010 ~8 nm); the crosslink STIFFNESS is NOT FOUND.
    assert FASCIN_INTERFILAMENT_SPACING_UM == pytest.approx(0.008)
    assert FASCIN_K_CROSSLINK_PN_PER_UM is None


def test_fascin_mechanics_requires_a_supplied_positive_stiffness() -> None:
    pairs = _FakeArray(600, (1, 2), wp.int32)
    active = _FakeArray(601, (1,), wp.int32)
    for bad in (0.0, -5.0, float("nan")):
        with pytest.raises(ValueError, match="UNSET PI-GAP slot"):
            FilopodiumFascinBundleMechanics(
                n_crosslink=1,
                crosslink_pairs_d=pairs,
                crosslink_active_d=active,
                k_fascin_pn_per_um=bad,
            )
    with pytest.raises(ValueError, match="rest spacing must be positive"):
        FilopodiumFascinBundleMechanics(
            n_crosslink=1,
            crosslink_pairs_d=pairs,
            crosslink_active_d=active,
            k_fascin_pn_per_um=100.0,
            rest_spacing_um=0.0,
        )


def test_fascin_mechanics_binds_the_real_crosslink_kernel_through_the_seam() -> None:
    pairs = _FakeArray(610, (1, 2), wp.int32)
    active = _FakeArray(611, (1,), wp.int32)
    launcher = _LaunchSpy([])
    mechanics = FilopodiumFascinBundleMechanics(
        n_crosslink=1,
        crosslink_pairs_d=pairs,
        crosslink_active_d=active,
        k_fascin_pn_per_um=250.0,   # PI-supplied stiffness activates the mechanism
        launch=launcher,
    )
    pos = _FakeArray(612, (2,), wp.vec3d)
    force = _FakeArray(613, (2,), wp.vec3d)
    mechanics.accumulate(pos, force)

    assert len(launcher.launches) == 1
    launch = launcher.launches[0]
    assert launch["kernel"] is mechanics.kernel
    assert launch["dim"] == 1
    p, cl_pairs, cl_active, k_fascin, rest = launch["inputs"]
    # no private state path: the kernel reads the actor-owned crosslink arrays and scatters into actor force.
    assert p.ptr == pos.ptr
    assert cl_pairs.ptr == pairs.ptr
    assert cl_active.ptr == active.ptr
    assert k_fascin == 250.0
    assert rest == FASCIN_INTERFILAMENT_SPACING_UM
    assert launch["outputs"][0].ptr == force.ptr


def test_fascin_mechanics_rejects_a_fake_kernel_and_foreign_topology() -> None:
    pairs = _FakeArray(620, (1, 2), wp.int32)
    active = _FakeArray(621, (1,), wp.int32)
    with pytest.raises(ValueError, match="real _fascin_bundle_kernel"):
        FilopodiumFascinBundleMechanics(
            n_crosslink=1,
            crosslink_pairs_d=pairs,
            crosslink_active_d=active,
            k_fascin_pn_per_um=100.0,
            kernel=object(),
        )
    mechanics = FilopodiumFascinBundleMechanics(
        n_crosslink=1,
        crosslink_pairs_d=pairs,
        crosslink_active_d=active,
        k_fascin_pn_per_um=100.0,
        launch=_LaunchSpy([]),
    )
    foreign_pairs = _FakeArray(630, (1, 2), wp.int32)
    with pytest.raises(ValueError, match="actor-owned crosslink_pairs_d"):
        mechanics.assert_bound_to(
            crosslink_pairs_d=foreign_pairs,
            crosslink_active_d=active,
            n_crosslink=1,
        )


def test_filopodium_owner_accepts_the_fascin_backend_and_rejects_a_private_copy() -> None:
    owner, *_ = _filopodium()
    # the actor-owned crosslink arrays bind cleanly as a fascin mechanics backend.
    good = FilopodiumFascinBundleMechanics(
        n_crosslink=owner.n_crosslink_capacity,
        crosslink_pairs_d=owner.crosslink_pairs_d,
        crosslink_active_d=owner.crosslink_active_d,
        k_fascin_pn_per_um=180.0,
        launch=_LaunchSpy([]),
    )
    installed = replace(owner, mechanics=good)
    assert isinstance(installed.mechanics, FilopodiumFascinBundleMechanics)
    # a fascin backend carrying a PRIVATE copy of the crosslink topology is rejected at construction.
    private_pairs = _FakeArray(640, (1, 2), wp.int32)
    private_active = _FakeArray(641, (1,), wp.int32)
    private = FilopodiumFascinBundleMechanics(
        n_crosslink=owner.n_crosslink_capacity,
        crosslink_pairs_d=private_pairs,
        crosslink_active_d=private_active,
        k_fascin_pn_per_um=180.0,
        launch=_LaunchSpy([]),
    )
    with pytest.raises(ValueError, match="actor-owned crosslink_pairs_d"):
        replace(owner, mechanics=private)


def test_fascin_bundle_pair_reference_is_equal_and_opposite_and_signed() -> None:
    fi, fj, load = fascin_bundle_pair_reference(
        [0.0, 0.0, 0.0], [0.02, 0.0, 0.0],
        k_fascin_pn_per_um=100.0, rest_spacing_um=0.008,
    )
    # stretched beyond rest: force on i points toward j (+x), reaction on j is its negative.
    assert fi[0] > 0.0
    assert fj == [-fi[0], -fi[1], -fi[2]]
    assert load == pytest.approx(abs(100.0 * (0.02 - 0.008)))
    # rejects an unset stiffness (PI-GAP guard).
    with pytest.raises(ValueError, match="k_fascin_pn_per_um"):
        fascin_bundle_pair_reference(
            [0.0, 0.0, 0.0], [0.02, 0.0, 0.0],
            k_fascin_pn_per_um=0.0, rest_spacing_um=0.008,
        )
