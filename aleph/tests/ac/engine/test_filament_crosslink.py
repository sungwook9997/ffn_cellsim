"""Structural gates for the two-array adjoint filament-crosslink runtime — 3 families, 5 edges.

CUDA-free: the joint runtime is a recording double, so these gates check the BINDING CONTRACT and the
kinetics seam, never physics.  The two-array adjoint scatter is
``load_path._segment_joint_force_kernel``'s and it runs first on the GPU.

The gate that carries the most weight is :func:`test_the_internal_edge_is_refused_because_it_is_other_physics`.
The family has four edges and TWO physics: ``dorsal_arc_crosslink`` has both endpoints in one array and
binds ``link_spring_kernel`` directly, while the other three join components that own disjoint arrays and
need the two-array adjoint runtime.  Letting the wrong one bind here would silently produce a joint whose
two endpoints index the same array — a weld that looks like a connector.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pytest

from aleph.engine.contracts import ConnectorFamily, ConnectorScope, reference_cell_architecture
from aleph.engine.load_path import JointKind, JointState
from aleph.engine.filament_crosslink import (
    CROSSLINK_FAMILIES,
    FILAMENT_CROSSLINK_EDGES,
    FilamentCrosslinkConnector,
    FilamentCrosslinkJointSpec,
    filament_crosslink_edges,
)

CORTEX = "cortex"


@dataclass(slots=True)
class _RecordingJointRuntime:
    """Stands in for ``LoadPathJointRuntime``; records which hooks the connector delegated."""

    capacity: int = 4
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


def _connector(name: str, kinetics: object | None = None) -> tuple[FilamentCrosslinkConnector, _RecordingJointRuntime]:
    joints = _RecordingJointRuntime()
    return FilamentCrosslinkConnector(name=name, joint_runtime=joints, kinetics=kinetics), joints


def test_the_runtime_serves_five_edges_across_three_families_read_from_the_contract() -> None:
    arch = reference_cell_architecture()
    declared = tuple(
        c.name for c in arch.connectors
        if c.family in CROSSLINK_FAMILIES and c.scope is ConnectorScope.INTER_COMPONENT
    )
    assert filament_crosslink_edges() == declared
    assert FILAMENT_CROSSLINK_EDGES == declared
    assert len(declared) == 5
    assert {f.value for f in CROSSLINK_FAMILIES} == {"transient_actin", "plectin", "spectraplakin"}


def test_actin_anchor_is_absent_and_that_is_a_recorded_gap_not_an_omission() -> None:
    """Its nascent edges end on ``focal_adhesion``, which owns no geometry — nothing to pull against."""
    assert ConnectorFamily.ACTIN_ANCHOR not in CROSSLINK_FAMILIES
    arch = reference_cell_architecture()
    fa = arch.component("focal_adhesion")
    assert fa.owns_geometry is False


def test_every_edge_binds_the_one_runtime_with_endpoints_from_its_contract() -> None:
    arch = reference_cell_architecture()
    by_name = {c.name: c for c in arch.connectors}
    for name in FILAMENT_CROSSLINK_EDGES:
        connector, _ = _connector(name)
        contract = by_name[name]
        assert connector.name == name
        assert (connector.component_a, connector.component_b) == (contract.component_a, contract.component_b)
        assert connector.chemistry_card == contract.chemistry_card


def test_five_edges_are_one_mechanics_but_four_chemistries() -> None:
    """Why the rate law is injected rather than baked in."""
    cards = {name: _connector(name)[0].chemistry_card for name in FILAMENT_CROSSLINK_EDGES}
    assert cards["sf_cortex_transient"] == cards["lamellipodium_cortex_seam"] == "transient_actin_crosslink"
    assert cards["filopodium_cortex_root"] == "formin_fascin_root_coupling"
    assert cards["if_sf_plectin"] == "plectin_if_actin"
    assert cards["mt_sf_spectraplakin"] == "spectraplakin_mt_actin"
    assert len(set(cards.values())) == 4


def test_the_internal_edge_is_refused_because_it_is_other_physics() -> None:
    """``dorsal_arc_crosslink`` has both endpoints in ONE array; binding it here would be a weld."""
    with pytest.raises(ValueError, match="INTERNAL"):
        _connector("dorsal_arc_crosslink")


def test_an_edge_from_another_family_is_refused() -> None:
    with pytest.raises(ValueError, match="this runtime serves"):
        _connector("nmii_cortex_motor")


def test_an_undeclared_edge_is_refused() -> None:
    with pytest.raises(ValueError, match="not declared in the architecture"):
        _connector("cortex_seam_that_does_not_exist")


def test_propose_events_refuses_when_no_rate_law_is_bound() -> None:
    """All three declare kinetics=True; a transient crosslink without kinetics is a permanent weld."""
    connector, _ = _connector("if_sf_plectin")
    with pytest.raises(RuntimeError, match="permanent weld"):
        connector.propose_events()


def test_propose_events_delegates_to_the_bound_rate_law() -> None:
    kinetics = _RecordingKinetics()
    connector, joints = _connector("lamellipodium_cortex_seam", kinetics=kinetics)
    assert connector.propose_events() == "proposed"
    assert kinetics.calls == [joints], "the rate law proposes onto the joint runtime's candidate state"


def test_mechanics_and_transaction_hooks_all_delegate_to_the_joint_runtime() -> None:
    connector, joints = _connector("filopodium_cortex_root")
    connector.snapshot_candidate()
    connector.accumulate()
    connector.rollback(object())
    connector.commit_irreversible(object(), dt_phys=0.05, rng_seed=3)
    assert joints.calls == ["snapshot_candidate", "accumulate", "rollback", "commit_irreversible"]
    assert connector.n_joints == joints.capacity


def test_joint_kinds_are_appended_never_renumbered() -> None:
    """The value is written into device SoA, so existing members must keep their numbers.

    A PREFIX check, not a closed list.  It was written as ``== [0, 1, 2, 3, 4]`` and so refused an
    APPEND — the one change the invariant it is named for explicitly permits — while a renumbering that
    kept the length would have passed it just as happily.  What must hold is that the members already
    written into a device SoA keep their values and that the sequence stays gap-free and ordered;
    2026-08-09 appended four (contact, ratchet, dynein capture, nascent series).
    """
    assert [int(k) for k in JointKind] == list(range(len(JointKind)))
    assert len(JointKind) >= 5
    assert (int(JointKind.FA_COMPOSITE_SERIES), int(JointKind.DORSAL_ARC_INTERNAL)) == (0, 1)
    assert int(JointKind.TRANSIENT_ACTIN_CROSSLINK) == 2
    assert int(JointKind.PLECTIN_IF_ACTIN) == 3
    assert int(JointKind.SPECTRAPLAKIN_MT_ACTIN) == 4


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"stiffness_pn_per_um": 0.0}, "stiffness"),
        ({"stiffness_pn_per_um": float("nan")}, "stiffness"),
        ({"rest_um": -1.0}, "rest length"),
        ({"joint_id": -1}, "joint ID"),
    ],
)
def test_spec_refuses_unphysical_values(kwargs: dict, match: str) -> None:
    base = dict(
        joint_id=0, edge="sf_cortex_transient", port_a=object(), port_b=object(),
        stiffness_pn_per_um=12.0, rest_um=0.05, initial_state=JointState.ACTIN_ENGAGED,
    )
    with pytest.raises(ValueError, match=match):
        FilamentCrosslinkJointSpec(**(base | kwargs))
