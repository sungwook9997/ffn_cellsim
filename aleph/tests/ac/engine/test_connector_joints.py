"""Structural gates for the last eleven connector runtimes — two force laws, not eleven.

CUDA-free: the joint SoA is a recording double, so these gates check the BINDING CONTRACT (which edges
the runtime serves, where its identity comes from, which pair law each edge gets, what it refuses) and
never physics.  The scatter itself is
:class:`~aleph.engine.load_path.LoadPathJointRuntime`'s and runs first on the GPU.

Two gates carry most of the weight.

:func:`test_a_contact_edge_gets_the_unilateral_law_and_nothing_else_does` is the one that would catch
the expensive mistake.  A contact edge silently handed the bilateral kernel is a connector that PULLS —
a membrane towed inward by the cortex, an adhesion nobody declared — and every downstream shape
assertion would still pass, because the force array has the right dtype and the pair still sums to zero.

:func:`test_the_served_set_is_read_from_the_contract_not_listed_here` is the lesson four families
already taught on 2026-08-09: counting declared edges as implementation units was wrong arithmetic in
BOTH directions.  ``contact`` is one family whose five edges are one law; ``motor`` + ``actin_anchor`` +
``fa_clutch`` are three families whose remaining edges are another single law that was already written.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from aleph.engine.actor import CellActor
from aleph.engine.connector_joints import (
    AdhesionSeriesJointSpec,
    ConnectorJointSpec,
    ContractJointConnector,
    adhesion_series_groups,
    connector_joint_edges,
    contact_pair_reference,
    force_kernel_for,
    resolve_adhesion_series_joint,
    resolve_connector_joint,
)
from aleph.engine.contracts import (
    CellArchitecture,
    ConnectorContract,
    ConnectorFamily,
    ConnectorScope,
    reference_cell_architecture,
)
from aleph.engine.load_path import (
    ActorRecord,
    ActorRegistry,
    ElementKind,
    EndpointRole,
    JointKind,
    JointState,
    resolve_fa_series_group,
)

ARCH = reference_cell_architecture()
BY_NAME = {c.name: c for c in ARCH.connectors}

#: Every component that appears as a joint endpoint in this file, one actor each.
_COMPONENTS = (
    "membrane", "cortex", "ecm", "nucleus", "microtubule",
    "lamellipodium", "filopodium", "sf_arc",
)
REGISTRY = ActorRegistry(
    tuple(
        ActorRecord(
            component=name, actor_id=index, actor_generation=0,
            entity_id=index, entity_generation=0, n_elements=8,
        )
        for index, name in enumerate(_COMPONENTS)
    )
)
_ACTOR_ID = {name: index for index, name in enumerate(_COMPONENTS)}


def _port(component: str, element: int = 0, u: float = 0.5) -> object:
    """A valid segment port on ``component``.

    The ROLE is a placeholder and is deliberately not asserted anywhere: ``EndpointRole`` declares
    members for stress-fibre and collagen endpoints only, and adding one for a membrane quadrature or a
    cortical dynein site is a contract change and therefore a PI decision.  The resolver carries the
    caller's role through rather than asserting a role it has no declaration for.
    """
    from aleph.engine.load_path import PortRef

    index = _ACTOR_ID[component]
    return PortRef(
        component=component, actor_id=index, actor_generation=0,
        entity_id=index, entity_generation=0,
        element_kind=ElementKind.SEGMENT, element_id=element,
        local_coordinates=(u, 0.0, 0.0, 0.0),
        role=EndpointRole.TRANSVERSE_ARC_MATERIAL,
    )


@dataclass(slots=True)
class _RecordingJoints:
    """Stands in for ``LoadPathJointRuntime``; records which transaction hooks it was driven through."""

    capacity: int = 3
    calls: list[str] = field(default_factory=list)

    def accumulate(self) -> None:
        self.calls.append("accumulate")

    def snapshot_candidate(self) -> None:
        self.calls.append("snapshot_candidate")

    def rollback(self, accepted: object) -> None:
        self.calls.append("rollback")

    def commit_irreversible(self, accepted: object, dt_phys: float, rng_seed: int) -> None:
        self.calls.append("commit_irreversible")


@dataclass(slots=True)
class _RecordingKinetics:
    calls: list[object] = field(default_factory=list)

    def propose_events(self, joints: object, *args: object, **kwargs: object) -> str:
        self.calls.append(joints)
        return "proposed"


def _connector(name: str, kinetics: object | None = None) -> ContractJointConnector:
    """Build the connector for ``name`` around a recording SoA — no device allocation."""
    contract = BY_NAME[name]
    return ContractJointConnector(
        name=contract.name,
        component_a=contract.component_a,
        component_b=contract.component_b,
        joint_runtime=_RecordingJoints(),
        kinetics_declared=contract.kinetics,
        chemistry_card=contract.chemistry_card,
        bidirectional=contract.bidirectional,
        adjoint_transfer_required=contract.adjoint_transfer_required,
        kinetics=kinetics,
    )


# --- which edges this runtime serves, and why that is derived ------------------------------------


def test_the_served_set_is_read_from_the_contract_not_listed_here() -> None:
    """Every inter-component CONTACT edge, plus every MOTOR edge that is not an NMII crossbridge."""
    served = set(connector_joint_edges())
    expected = {
        c.name for c in ARCH.connectors
        if c.scope is ConnectorScope.INTER_COMPONENT
        and (
            c.family is ConnectorFamily.CONTACT
            or (c.family is ConnectorFamily.MOTOR and "nmii" not in (c.component_a, c.component_b))
        )
    }
    assert served == expected
    assert served == {
        "membrane_cortex_contact", "membrane_ecm_contact", "nucleus_cortex_contact",
        "lamellipodium_membrane_contact", "filopodium_membrane_tip", "mt_cortex_capture",
    }


def test_the_nmii_motor_edges_belong_to_the_actuator_and_are_not_claimed_here() -> None:
    """An ABSENCE gate: four MOTOR edges are head-resolved Stam-Hocky crossbridges, not this law.

    Claiming them here would put a second runtime on an edge that already has one, which
    ``CellActor.bind_connector`` refuses — but only after the wrong module had already been written.
    """
    served = set(connector_joint_edges())
    nmii_motors = {
        c.name for c in ARCH.connectors
        if c.family is ConnectorFamily.MOTOR and "nmii" in (c.component_a, c.component_b)
    }
    assert len(nmii_motors) == 4
    assert served.isdisjoint(nmii_motors)


def test_together_with_the_series_groups_every_open_edge_is_covered() -> None:
    """Six single edges + two new series (two edges each) + the medium traction = the open eleven."""
    single = set(connector_joint_edges())
    nascent = {
        c.name for c in ARCH.connectors
        if c.mechanical_group in {"lamellipodium_nascent_series", "filopodium_nascent_series"}
    }
    assert len(nascent) == 4
    assert len(single | nascent) == 10  # membrane_medium_traction is an ENVIRONMENT_BOUNDARY facade


# --- the pair law, which is the mistake that would not otherwise be caught ------------------------


def test_a_contact_edge_gets_the_unilateral_law_and_nothing_else_does() -> None:
    """A contact that pulls is an undeclared adhesion; the choice is made from the family, not the caller."""
    for name in connector_joint_edges():
        contract = BY_NAME[name]
        kernel = force_kernel_for(contract)
        if contract.family is ConnectorFamily.CONTACT:
            assert kernel is not None, f"{name} must not fall back to the bilateral spring"
        else:
            assert kernel is None, f"{name} is a bound crossbridge and resists in both directions"


def test_the_unilateral_law_never_pulls_and_vanishes_at_the_physiological_gap() -> None:
    """Sign sense, as an executable statement: overlap repels, separation is force-free."""
    d0, k = 0.030, 500.0
    assert contact_pair_reference(0.020, d0, k) < 0.0, "overlapping pairs must be pushed apart"
    assert contact_pair_reference(d0, d0, k) == 0.0, "force-free exactly at the contact distance"
    assert contact_pair_reference(0.050, d0, k) == 0.0, "a separated contact may not pull"
    # Linear in the overlap, so a stiffness is the only magnitude the caller supplies.
    assert contact_pair_reference(0.010, d0, k) == pytest.approx(2.0 * contact_pair_reference(0.020, d0, k))


def test_the_two_contact_physics_are_told_apart_by_kinetics_not_by_the_chemistry_card() -> None:
    """Three steric edges, two Brownian-ratchet ones — read from ``kinetics``, which cannot be mislabelled."""
    kinds = {}
    for name in connector_joint_edges():
        contract = BY_NAME[name]
        if contract.family is not ConnectorFamily.CONTACT:
            continue
        spec = ConnectorJointSpec(
            joint_id=0, edge=name,
            port_a=_port(contract.component_a), port_b=_port(contract.component_b),
            stiffness_pn_per_um=100.0, rest_um=0.03,
        )
        kinds[name] = resolve_connector_joint(spec, REGISTRY, ARCH).kind
    steric = {n for n, k in kinds.items() if k is JointKind.STERIC_CONTACT}
    ratchet = {n for n, k in kinds.items() if k is JointKind.BROWNIAN_RATCHET_CONTACT}
    assert steric == {"membrane_cortex_contact", "membrane_ecm_contact", "nucleus_cortex_contact"}
    assert ratchet == {"lamellipodium_membrane_contact", "filopodium_membrane_tip"}
    assert all(not BY_NAME[n].kinetics for n in steric)
    assert all(BY_NAME[n].kinetics for n in ratchet)


def test_the_capture_edge_resolves_as_a_bound_crossbridge() -> None:
    contract = BY_NAME["mt_cortex_capture"]
    spec = ConnectorJointSpec(
        joint_id=0, edge="mt_cortex_capture",
        port_a=_port(contract.component_a), port_b=_port(contract.component_b),
        stiffness_pn_per_um=50.0, rest_um=0.02,
    )
    assert resolve_connector_joint(spec, REGISTRY, ARCH).kind is JointKind.CORTICAL_DYNEIN_CAPTURE


# --- identity and refusals -------------------------------------------------------------------------


def test_endpoints_come_from_the_contract_not_the_caller() -> None:
    for name in connector_joint_edges():
        contract = BY_NAME[name]
        connector = _connector(name)
        assert (connector.component_a, connector.component_b) == (
            contract.component_a, contract.component_b
        )
        assert connector.chemistry_card == contract.chemistry_card
        assert connector.adjoint_transfer_required is contract.adjoint_transfer_required
        assert connector.mechanical_group is None


def test_a_port_on_the_wrong_component_is_refused() -> None:
    spec = ConnectorJointSpec(
        joint_id=0, edge="nucleus_cortex_contact",
        port_a=_port("membrane"), port_b=_port("cortex"),
        stiffness_pn_per_um=100.0, rest_um=0.03,
    )
    with pytest.raises(ValueError, match="endpoint A must belong to 'nucleus'"):
        resolve_connector_joint(spec, REGISTRY, ARCH)


@pytest.mark.parametrize(
    "edge",
    ["nmii_cortex_motor", "sf_cortex_transient", "dorsal_arc_crosslink", "surface_porous_transfer"],
)
def test_an_edge_this_runtime_does_not_serve_is_refused(edge: str) -> None:
    """A crosslink, an internal joint and a field transfer all have their own runtimes already."""
    contract = BY_NAME[edge]
    spec = ConnectorJointSpec(
        joint_id=0, edge=edge,
        port_a=_port("cortex"), port_b=_port("cortex"),
        stiffness_pn_per_um=1.0, rest_um=0.01,
    )
    with pytest.raises(ValueError, match="not served by this runtime"):
        resolve_connector_joint(spec, REGISTRY, ARCH)
    assert contract is not None


def test_an_undeclared_edge_name_is_refused() -> None:
    spec = ConnectorJointSpec(
        joint_id=0, edge="contact_that_does_not_exist",
        port_a=_port("cortex"), port_b=_port("membrane"),
        stiffness_pn_per_um=1.0, rest_um=0.01,
    )
    with pytest.raises(ValueError, match="not declared in the architecture"):
        resolve_connector_joint(spec, REGISTRY, ARCH)


# --- kinetics: refuse rather than silently propose nothing ------------------------------------------


def test_a_non_kinetic_contact_has_no_events_and_says_so() -> None:
    connector = _connector("membrane_cortex_contact")
    with pytest.raises(RuntimeError, match="declares kinetics=False"):
        connector.propose_events()


def test_a_kinetic_edge_without_a_rate_law_refuses_rather_than_proposing_nothing() -> None:
    """Silently proposing nothing is indistinguishable from a permanent weld."""
    connector = _connector("lamellipodium_membrane_contact")
    with pytest.raises(RuntimeError, match="permanent weld"):
        connector.propose_events()


def test_a_bound_rate_law_receives_the_joint_soa() -> None:
    kinetics = _RecordingKinetics()
    connector = _connector("mt_cortex_capture", kinetics=kinetics)
    assert connector.propose_events() == "proposed"
    assert len(kinetics.calls) == 1


def test_every_transaction_hook_reaches_the_joint_soa() -> None:
    connector = _connector("membrane_ecm_contact")
    joints = connector._joints  # noqa: SLF001 - the gate is about what the facade forwards
    connector.snapshot_candidate()
    connector.rollback(accepted=object())
    connector.commit_irreversible(accepted=object(), dt_phys=1.0e-4, rng_seed=7)
    connector.accumulate()
    assert joints.calls == [
        "snapshot_candidate", "rollback", "commit_irreversible", "accumulate",
    ]


def test_the_facade_mechanics_aliases_all_dispatch_one_launch() -> None:
    """One force launch per call however the owning facade's protocol spells it."""
    for alias in ("accumulate_contact", "accumulate_actor", "accumulate_motor"):
        connector = _connector("membrane_ecm_contact")
        getattr(connector, alias)(object(), object())
        assert connector._joints.calls == ["accumulate"]  # noqa: SLF001


# --- the adhesion series ----------------------------------------------------------------------------


def test_every_declared_series_chains_through_one_geometry_less_middle() -> None:
    """What makes a series a series is structural, and is now what is checked."""
    groups = adhesion_series_groups()
    assert set(groups) == {
        "alpha2beta1_collagen_series", "lamellipodium_nascent_series", "filopodium_nascent_series",
    }
    for group in groups:
        resolution = resolve_fa_series_group(ARCH, group)
        assert resolution.middle_component == "focal_adhesion"
        assert resolution.ligand_component == "ecm"
        assert resolution.mechanical_joint_count == 1, "each series is ONE composite spring"
    assert {resolve_fa_series_group(ARCH, g).actin_component for g in groups} == {
        "sf_arc", "lamellipodium", "filopodium",
    }


def test_a_group_whose_edges_do_not_chain_is_refused() -> None:
    """Two edges sharing no component are two unrelated springs filed under one name."""
    broken = CellArchitecture(
        components=ARCH.components,
        connectors=(
            ConnectorContract(
                "anchor_x", ConnectorFamily.ACTIN_ANCHOR, "sf_arc", "focal_adhesion", True, True,
                chemistry_card="actin_talin_integrin", mechanical_group="broken_series",
            ),
            ConnectorContract(
                "clutch_x", ConnectorFamily.FA_CLUTCH, "membrane", "ecm", True, True,
                chemistry_card="alpha2beta1_collagen", mechanical_group="broken_series",
            ),
        ),
    )
    with pytest.raises(ValueError, match="share exactly one"):
        resolve_fa_series_group(broken, "broken_series")


def test_the_series_spring_runs_actin_to_collagen_and_has_no_fa_port() -> None:
    """``focal_adhesion`` owns no geometry, so there is no third endpoint to supply."""
    spec = AdhesionSeriesJointSpec(
        joint_id=0, group="lamellipodium_nascent_series", fa_cluster_id=4,
        actin_port=_port("lamellipodium"), ligand_port=_port("ecm"),
        stiffness_pn_per_um=250.0, rest_um=0.05,
    )
    joint = resolve_adhesion_series_joint(spec, REGISTRY, ARCH)
    assert joint.kind is JointKind.NASCENT_FA_SERIES
    assert joint.port_a.component == "lamellipodium"
    assert joint.port_b.component == "ecm"
    assert joint.fa_cluster_id == 4
    assert set(joint.semantic_edges) == {
        "lamellipodium_nascent_fa", "lamellipodium_nascent_clutch"
    }, "one spring registered under BOTH semantic edges"


def test_a_series_port_on_the_wrong_actin_population_is_refused() -> None:
    """The filopodial series may not be built over lamellipodial actin."""
    spec = AdhesionSeriesJointSpec(
        joint_id=0, group="filopodium_nascent_series", fa_cluster_id=0,
        actin_port=_port("lamellipodium"), ligand_port=_port("ecm"),
        stiffness_pn_per_um=250.0, rest_um=0.05,
    )
    with pytest.raises(ValueError, match="actin-side port must belong to 'filopodium'"):
        resolve_adhesion_series_joint(spec, REGISTRY, ARCH)


def test_one_composite_runtime_binds_both_series_edges_and_a_second_is_refused() -> None:
    """Each nascent adhesion gets its OWN group; two runtimes in one group would be two springs."""
    actor = CellActor(ARCH)
    composite = ContractJointConnector(
        name="lamellipodium_nascent_series",
        component_a="lamellipodium", component_b="ecm",
        joint_runtime=_RecordingJoints(), kinetics_declared=True,
        chemistry_card="actin_talin_integrin_nascent + alpha2beta1_collagen_nascent",
        bidirectional=True, adjoint_transfer_required=True,
        mechanical_group="lamellipodium_nascent_series",
    )
    actor.bind_connector("lamellipodium_nascent_fa", composite)
    actor.bind_connector("lamellipodium_nascent_clutch", composite)
    assert actor.connector_runtime("lamellipodium_nascent_clutch") is composite

    rival = CellActor(ARCH)
    rival.bind_connector("filopodium_nascent_fa", composite)
    with pytest.raises(ValueError, match="must bind one composite runtime"):
        rival.bind_connector("filopodium_nascent_clutch", ContractJointConnector(
            name="filopodium_nascent_series",
            component_a="filopodium", component_b="ecm",
            joint_runtime=_RecordingJoints(), kinetics_declared=True,
            chemistry_card="x", bidirectional=True, adjoint_transfer_required=True,
            mechanical_group="filopodium_nascent_series",
        ))


def test_the_nascent_series_does_not_join_the_mature_group() -> None:
    """Folding them together would make one spring serve three different adhesions."""
    mature = {c.name for c in ARCH.mechanical_group("alpha2beta1_collagen_series")}
    assert mature == {"fa_actin_anchor", "integrin_collagen_clutch"}


# --- what this runtime deliberately does NOT own -----------------------------------------------------


def test_the_connector_owns_no_geometry_and_no_force_array() -> None:
    """An ABSENCE gate.

    The pair force lands in the two ENDPOINT components' own arrays.  A connector that also carried a
    force array would contribute a resultant into the balance ledger that its endpoints already
    contributed — the two channels agreeing for a reason unrelated to the adjoint wiring under test,
    which is the same trap ``test_focal_adhesion`` guards with its own absence assertions.
    """
    connector = _connector("nucleus_cortex_contact")
    for forbidden in ("position_d", "force_d", "owned_arrays", "geometry"):
        assert not hasattr(connector, forbidden), f"a connector must not own {forbidden}"


def test_the_initial_state_gate_is_the_default_and_below_it_nothing_carries_force() -> None:
    spec = ConnectorJointSpec(
        joint_id=0, edge="membrane_cortex_contact",
        port_a=_port("membrane"), port_b=_port("cortex"),
        stiffness_pn_per_um=100.0, rest_um=0.03,
    )
    assert spec.initial_state is JointState.ACTIN_ENGAGED
    assert int(JointState.FREE) < int(JointState.ACTIN_ENGAGED)
