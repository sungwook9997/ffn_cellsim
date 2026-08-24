"""Structural gates for the ``focal_adhesion`` clutch-graph owner — the component that had no owner.

CUDA-free by construction: every array is a metadata-only double and every launch/copy goes to a recorder,
so these gates check the OWNERSHIP CONTRACT (what is owned, what is snapshotted, what the device predicate
gates, and what is deliberately ABSENT) and never physics.  The kernel bodies run for the first time on the
GPU, in the composed native driver.

What each gate is for:
  * the transaction API is all-or-nothing — :class:`~aleph.engine.world.CellWorldTransaction` raises
    ``TypeError`` on a partial one, so a missing hook must fail here rather than at world construction;
  * the owner must carry NO force/position array and NO ``accumulate_ledger`` — the FA spring belongs to the
    ONE composite joint runtime, and an owner that also carried it would turn one series clutch into two in
    parallel.  That absence is a contract, so it is asserted, not left to a reader's memory;
  * ``capacity == 0`` (a cell with no engaged adhesion) is legal and every hook is a no-op — the Sanity
    Gate's boundary case;
  * a rejected step restores every owned array, and an accepted step promotes the proposal — which is the
    charter's rule that a kinetic transition commits only on an accepted physical step.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest
import warp as wp

from aleph.engine.focal_adhesion import (
    FA_COMPONENT,
    FA_COMPOSITE_EDGES,
    FocalAdhesionStateOwner,
    _bump_epoch_if_accepted_kernel,
    _commit_clutch_if_accepted_kernel,
    _restore_float64_if_rejected_kernel,
    _restore_int32_if_rejected_kernel,
    build_focal_adhesion_state_owner,
)
from aleph.engine.load_path import JointState

CAPACITY = 5


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """Metadata-only CUDA array double; it never executes authoritative physics."""

    tag: str
    dtype: object = wp.int32


@dataclass(slots=True)
class _LaunchRecorder:
    calls: list[tuple[object, int, list, object]] = field(default_factory=list)

    def __call__(self, kernel: object, *, dim: int, inputs: list, device: object | None = None) -> None:
        self.calls.append((kernel, dim, inputs, device))

    def kernels(self) -> list[object]:
        return [call[0] for call in self.calls]


@dataclass(slots=True)
class _CopyRecorder:
    pairs: list[tuple[object, object]] = field(default_factory=list)

    def __call__(self, dst: object, src: object) -> None:
        self.pairs.append((dst, src))


def _owner(capacity: int = CAPACITY) -> tuple[FocalAdhesionStateOwner, _LaunchRecorder, _CopyRecorder]:
    launch, copy = _LaunchRecorder(), _CopyRecorder()
    owner = FocalAdhesionStateOwner(
        name=FA_COMPONENT,
        device="cuda:0",
        capacity=capacity,
        active_d=_FakeArray("active"),
        cluster_id_d=_FakeArray("cluster_id"),
        state_d=_FakeArray("state"),
        candidate_state_d=_FakeArray("candidate_state"),
        age_d=_FakeArray("age", wp.float64),
        rng_epoch_d=_FakeArray("rng_epoch"),
        state_snap_d=_FakeArray("state_snap"),
        candidate_state_snap_d=_FakeArray("candidate_state_snap"),
        age_snap_d=_FakeArray("age_snap", wp.float64),
        rng_epoch_snap_d=_FakeArray("rng_epoch_snap"),
        epoch_d=_FakeArray("epoch"),
        launch=launch,
        copy=copy,
    )
    return owner, launch, copy


def test_transaction_api_is_complete_so_the_world_admits_it_as_a_participant() -> None:
    """All three hooks present: a partial API is a ``TypeError`` at world construction, not a skip."""
    owner, _, _ = _owner()
    for hook in ("snapshot_candidate", "rollback", "commit_irreversible"):
        assert callable(getattr(owner, hook, None)), f"missing transaction hook {hook!r}"


def test_owner_carries_no_mechanics_and_no_ledger_because_the_spring_is_the_composite_joint() -> None:
    """The absence is the contract: ``owns_geometry=False`` and one series clutch, not two in parallel."""
    owner, _, _ = _owner()
    assert not hasattr(owner, "geometry_arrays")
    assert not hasattr(owner, "accumulate_ledger")
    assert not hasattr(owner, "force_d")
    assert not hasattr(owner, "position_d")


def test_the_two_semantic_edges_are_named_so_a_driver_cannot_bind_two_runtimes() -> None:
    assert FA_COMPOSITE_EDGES == ("fa_actin_anchor", "integrin_collagen_clutch")


def test_snapshot_copies_every_owned_mutable_array_and_nothing_else() -> None:
    owner, _, copy = _owner()
    owner.snapshot_candidate()
    assert [(dst.tag, src.tag) for dst, src in copy.pairs] == [
        ("state_snap", "state"),
        ("candidate_state_snap", "candidate_state"),
        ("age_snap", "age"),
        ("rng_epoch_snap", "rng_epoch"),
    ]


def test_rollback_restores_every_owned_array_under_the_device_predicate() -> None:
    owner, launch, _ = _owner()
    accepted = _FakeArray("accepted")
    owner.rollback(accepted)

    assert launch.kernels() == [
        _restore_int32_if_rejected_kernel,
        _restore_int32_if_rejected_kernel,
        _restore_int32_if_rejected_kernel,
        _restore_float64_if_rejected_kernel,
    ]
    restored = [call[2][1].tag for call in launch.calls]
    assert restored == ["state", "candidate_state", "rng_epoch", "age"]
    for kernel, dim, inputs, device in launch.calls:
        assert dim == CAPACITY
        assert device == "cuda:0"
        assert inputs[0] is accepted, "every restore must be gated on the device predicate"


def test_commit_promotes_the_proposal_and_bumps_the_epoch_under_the_predicate() -> None:
    owner, launch, _ = _owner()
    accepted = _FakeArray("accepted")
    owner.commit_irreversible(accepted, dt_phys=0.05, rng_seed=7)

    assert launch.kernels() == [_commit_clutch_if_accepted_kernel, _bump_epoch_if_accepted_kernel]
    commit_inputs = launch.calls[0][2]
    assert commit_inputs[0] is accepted
    assert float(commit_inputs[1]) == pytest.approx(0.05)
    assert [a.tag for a in commit_inputs[2:]] == [
        "active", "state", "candidate_state", "state_snap", "age",
    ]
    assert launch.calls[1][1] == 1, "the accepted-step epoch is one scalar, not one per clutch"


def test_empty_clutch_population_is_legal_and_every_hook_is_a_no_op() -> None:
    """A cell with no engaged adhesion is a boundary case, not an error."""
    owner, launch, copy = _owner(capacity=0)
    owner.snapshot_candidate()
    owner.rollback(_FakeArray("accepted"))
    assert copy.pairs == []
    assert launch.calls == []

    owner.commit_irreversible(_FakeArray("accepted"), dt_phys=0.05, rng_seed=0)
    assert launch.kernels() == [_bump_epoch_if_accepted_kernel], "the epoch advances even with no clutches"


def test_owner_refuses_a_name_that_is_not_the_declared_component() -> None:
    with pytest.raises(ValueError, match="focal_adhesion"):
        FocalAdhesionStateOwner(
            name="cortex", device="cuda:0", capacity=0,
            active_d=_FakeArray("a"), cluster_id_d=_FakeArray("b"), state_d=_FakeArray("c"),
            candidate_state_d=_FakeArray("d"), age_d=_FakeArray("e"), rng_epoch_d=_FakeArray("f"),
            state_snap_d=_FakeArray("g"), candidate_state_snap_d=_FakeArray("h"),
            age_snap_d=_FakeArray("i"), rng_epoch_snap_d=_FakeArray("j"), epoch_d=_FakeArray("k"),
        )


@pytest.mark.parametrize(
    "ids, match",
    [
        (np.zeros((2, 2), dtype=np.int32), "1-D"),
        (np.array([0.5, 1.5]), "nonnegative integers"),
        (np.array([0, -1], dtype=np.int32), "nonnegative integers"),
    ],
)
def test_builder_rejects_malformed_cluster_ids_before_it_allocates(ids: np.ndarray, match: str) -> None:
    """Validation precedes allocation, so these raise with no device touched."""
    with pytest.raises(ValueError, match=match):
        build_focal_adhesion_state_owner(cluster_ids=ids, device="cuda:0")


def test_initial_state_default_matches_the_composite_joint_spec_default() -> None:
    """Owner and joint must not disagree at step zero; both default to ``ACTIN_ENGAGED``."""
    from aleph.engine.load_path import CompositeFAJointSpec

    assert CompositeFAJointSpec.__dataclass_fields__["initial_state"].default is JointState.ACTIN_ENGAGED
