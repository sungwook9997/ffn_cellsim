"""Structural gates for ``ecm_crosslink`` — the collagen bond population that was declared and never built.

CUDA-free: launches and copies go to recorders, so these gates check the BINDING CONTRACT and the
accepted-step discipline, never physics.  The gated spring kernel runs first on the GPU.

Two gates carry the weight.  :func:`test_it_satisfies_the_protocol_ecm_world_already_validates` is the
point of the module — ``ECMWorld`` has held a ``crosslink_connector`` field and validated it for as long
as the tier-(a) ECM row has read ``zero crosslink bonds``; nothing satisfied the Protocol.  And
:func:`test_the_population_is_fixed_so_a_dissociated_bond_is_gated_not_removed` records why the kernel
skips unbound slots instead of the builder compacting them: compaction would change the bond count
between steps and make the ECM population ledger non-conserved.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest

from aleph.engine.collagen_crosslink import (
    CollagenCrosslinkConnector,
    _bound_link_spring_kernel,
    _commit_bond_if_accepted_kernel,
    _restore_float64_if_rejected_kernel,
    _restore_int32_if_rejected_kernel,
    build_collagen_crosslink_connector,
)
from aleph.engine.ecm_world import ECM_COMPONENT, ECM_CROSSLINK_CONNECTOR, ECMInternalCrosslink

BONDS = np.array([[0, 1], [1, 2], [2, 3]], dtype=np.int32)
K = np.array([12.0, 12.0, 12.0])
R0 = np.array([0.1, 0.1, 0.1])


@dataclass(frozen=True, slots=True)
class _FakeArray:
    tag: str


@dataclass(slots=True)
class _LaunchRecorder:
    calls: list[tuple[object, int, list, object]] = field(default_factory=list)

    def __call__(self, kernel: object, *, dim: int, inputs: list, device: object | None = None) -> None:
        self.calls.append((kernel, dim, inputs, device))

    def kernels(self) -> list[object]:
        return [c[0] for c in self.calls]


@dataclass(slots=True)
class _CopyRecorder:
    pairs: list[tuple[object, object]] = field(default_factory=list)

    def __call__(self, dst: object, src: object) -> None:
        self.pairs.append((dst, src))


def _connector(capacity: int = 3) -> tuple[CollagenCrosslinkConnector, _LaunchRecorder, _CopyRecorder]:
    launch, copy = _LaunchRecorder(), _CopyRecorder()
    tags = [
        "bonds", "bound", "cand", "k", "r0", "age", "bound_snap", "cand_snap", "age_snap",
    ]
    arrays = {t: _FakeArray(t) for t in tags}
    connector = CollagenCrosslinkConnector(
        name=ECM_CROSSLINK_CONNECTOR, component_a=ECM_COMPONENT, component_b=ECM_COMPONENT,
        device="cuda:0", capacity=capacity,
        bonds_d=arrays["bonds"], bound_d=arrays["bound"], candidate_bound_d=arrays["cand"],
        k_d=arrays["k"], r0_d=arrays["r0"], age_d=arrays["age"],
        bound_snap_d=arrays["bound_snap"], candidate_bound_snap_d=arrays["cand_snap"],
        age_snap_d=arrays["age_snap"], launch=launch, copy=copy,
    )
    return connector, launch, copy


def test_it_satisfies_the_protocol_ecm_world_already_validates() -> None:
    """The field and its validation existed; what was missing was anything that satisfied them."""
    connector, _, _ = _connector()
    assert isinstance(connector, ECMInternalCrosslink)
    assert connector.name == ECM_CROSSLINK_CONNECTOR
    assert {connector.component_a, connector.component_b} == {ECM_COMPONENT}
    assert connector.adjoint_transfer_required is True


def test_it_refuses_the_three_things_ecm_world_would_refuse() -> None:
    for kwargs, match in (
        ({"name": "membrane_ecm_contact"}, "ecm_crosslink"),
        ({"component_b": "membrane"}, "INTERNAL"),
        ({"adjoint_transfer_required": False}, "adjoint"),
    ):
        base = dict(
            name=ECM_CROSSLINK_CONNECTOR, component_a=ECM_COMPONENT, component_b=ECM_COMPONENT,
            device="cuda:0", capacity=0,
            bonds_d=None, bound_d=None, candidate_bound_d=None, k_d=None, r0_d=None, age_d=None,
            bound_snap_d=None, candidate_bound_snap_d=None, age_snap_d=None,
        )
        with pytest.raises(ValueError, match=match):
            CollagenCrosslinkConnector(**(base | kwargs))


def test_the_population_is_fixed_so_a_dissociated_bond_is_gated_not_removed() -> None:
    """Compaction would change the bond count between steps and break population conservation."""
    connector, launch, _ = _connector()
    connector.accumulate_internal(_FakeArray("pos"), _FakeArray("force"))
    assert launch.kernels() == [_bound_link_spring_kernel]
    kernel, dim, inputs, _ = launch.calls[0]
    assert dim == 3, "every slot is launched; the kernel skips the unbound ones"
    assert [a.tag for a in inputs] == ["pos", "bonds", "bound", "k", "r0", "force"]


def test_snapshot_covers_bond_state_and_age_and_nothing_immutable() -> None:
    connector, _, copy = _connector()
    connector.snapshot_candidate()
    assert [(d.tag, s.tag) for d, s in copy.pairs] == [
        ("bound_snap", "bound"), ("cand_snap", "cand"), ("age_snap", "age"),
    ]


def test_rejected_step_restores_state_and_age_under_the_device_predicate() -> None:
    connector, launch, _ = _connector()
    accepted = _FakeArray("accepted")
    connector.rollback(accepted)
    assert launch.kernels() == [
        _restore_int32_if_rejected_kernel,
        _restore_int32_if_rejected_kernel,
        _restore_float64_if_rejected_kernel,
    ]
    assert [c[2][1].tag for c in launch.calls] == ["bound", "cand", "age"]
    assert all(c[2][0] is accepted for c in launch.calls)


def test_commit_promotes_the_proposal_only_under_acceptance() -> None:
    connector, launch, _ = _connector()
    accepted = _FakeArray("accepted")
    connector.commit_irreversible(accepted, dt_phys=0.05, rng_seed=11)
    assert launch.kernels() == [_commit_bond_if_accepted_kernel]
    inputs = launch.calls[0][2]
    assert inputs[0] is accepted
    assert float(inputs[1]) == pytest.approx(0.05)
    assert [a.tag for a in inputs[2:]] == ["bound", "cand", "bound_snap", "age"]


def test_an_uncrosslinked_network_is_legal_and_every_hook_is_a_no_op() -> None:
    """Zero bonds is the configuration the tier-(a) native ECM row was measured on."""
    connector, launch, copy = _connector(capacity=0)
    connector.accumulate_internal(_FakeArray("pos"), _FakeArray("force"))
    connector.snapshot_candidate()
    connector.rollback(_FakeArray("accepted"))
    connector.commit_irreversible(_FakeArray("accepted"), dt_phys=0.05, rng_seed=0)
    assert launch.calls == []
    assert copy.pairs == []
    assert connector.bound_count() == 0


def test_ledger_contributes_the_topology_count_not_a_force() -> None:
    """The force resultant belongs to the ECM owner's own array; this reports bond topology."""
    recorded: list[tuple[str, int]] = []

    class _Ledger:
        def add_topology_count(self, name: str, count: int) -> None:
            recorded.append((name, count))

    connector, _, _ = _connector()
    connector.accumulate_ledger(_Ledger())
    assert recorded == [(ECM_CROSSLINK_CONNECTOR, 3)]
    connector.accumulate_ledger(object())  # a ledger with no such accessor is left alone


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"bonds": np.zeros((3, 3), dtype=np.int32)}, r"\(J, 2\)"),
        ({"bonds": np.array([[0, 0]], dtype=np.int32), "stiffness_pn_per_um": np.array([1.0]),
          "rest_um": np.array([0.1])}, "node to itself"),
        ({"stiffness_pn_per_um": np.array([0.0, 12.0, 12.0])}, "stiffness"),
        ({"rest_um": np.array([-1.0, 0.1, 0.1])}, "rest length"),
        ({"stiffness_pn_per_um": np.array([12.0])}, "one entry per bond"),
    ],
)
def test_builder_refuses_a_malformed_bond_table_before_it_allocates(kwargs: dict, match: str) -> None:
    base = dict(bonds=BONDS, stiffness_pn_per_um=K, rest_um=R0, device="cuda:0")
    with pytest.raises(ValueError, match=match):
        build_collagen_crosslink_connector(**(base | kwargs))
