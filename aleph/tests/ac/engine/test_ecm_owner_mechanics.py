r"""Owner-driven gate for the ``ecm`` collagen constitutive force pass (CPU-safe; native forbidden).

Proves the advance from KERNEL_BOUND → owner-driven CONNECTED for the ``ecm`` component: the
:class:`~aleph.engine.ecm_world.ECMStateOwner` now CONSTRUCTS + OWNS the concrete
:class:`~aleph.engine.ecm_mechanics.ECMConstitutiveForce` and drives it through its own per-step
force hook (:meth:`ECMStateOwner.accumulate`) and its accepted-step commit
(:meth:`ECMStateOwner.commit_irreversible` re-derives the constitutive parameters after a remesh).

Every proof is CUDA-free (metadata-only ``_FakeArray`` doubles + a recording launcher — no Warp-CPU
production-kernel launch), and additivity is checked: a mechanics delegate without ``refresh_parameters``
(the build-once cortex/SF style or a plain spy) is left completely untouched, so existing components and
``test_ecm_world.py`` stay bit-identical.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import warp as wp

from aleph.engine.ecm_mechanics import ECMConstitutiveForce
from aleph.engine.ecm_world import (
    ECM_COMPONENT,
    BoundaryAnchorMode,
    ECMStateOwner,
    ECMWorldSettings,
    build_ecm_state_owner,
)
from aleph.engine.load_path import ActorRecord


# ── metadata-only CUDA doubles + recording launcher (no host kernel execution) ────────────────────────
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


@dataclass(slots=True)
class _LaunchRecorder:
    calls: list[dict[str, object]] = field(default_factory=list)

    def __call__(self, kernel, *, dim, inputs, device) -> None:
        self.calls.append({"kernel": kernel, "dim": dim, "inputs": inputs, "device": device})


@dataclass(slots=True)
class _TransactionSpy:
    calls: list[tuple[object, ...]] = field(default_factory=list)

    def snapshot_candidate(self) -> None:
        self.calls.append(("snapshot",))

    def rollback(self, accepted: object) -> None:
        self.calls.append(("rollback", accepted))

    def commit_irreversible(self, accepted: object, dt_phys: float, rng_seed: int) -> None:
        self.calls.append(("commit", accepted, dt_phys, rng_seed))


@dataclass(slots=True)
class _LedgerSpy:
    calls: list[object] = field(default_factory=list)

    def accumulate_ledger(self, ledger: object) -> None:
        self.calls.append(ledger)


@dataclass(slots=True)
class _PlainMechanicsSpy:
    """A build-once mechanics delegate WITHOUT a ``refresh_parameters`` hook (additivity control)."""

    n_nodes: int = 8
    calls: list[tuple[object, object]] = field(default_factory=list)

    def accumulate(self, pos: object, force: object) -> None:
        self.calls.append((pos, force))


_LINK_KERNEL = object()
_BENDING_KERNEL = object()
_STIFFNESS_KERNEL = object()
_ALPHA_KERNEL = object()

# EA = E_fibril·π r_f² (GAP value from the collagen card, for the double); κ = k_B·T·L_p (SOURCED).
_EA_PN = 1.1e6 * 3.141592653589793 * 0.05 ** 2
_KAPPA = 4.28e-3 * 17.0


def _constitutive_force(
    recorder: _LaunchRecorder,
    *,
    n_nodes: int = 8,
    n_segments: int = 4,
    n_bends: int = 3,
) -> ECMConstitutiveForce:
    """A real :class:`ECMConstitutiveForce` bound to fakes + a recording launcher (no CUDA)."""
    return ECMConstitutiveForce(
        device="cuda:0",
        n_nodes=n_nodes,
        segments_d=_FakeArray(500, (n_segments, 2), wp.int32),
        seg_rest_d=_FakeArray(501, (n_segments,), wp.float64),
        segment_active_d=_FakeArray(502, (n_segments,), wp.int32),
        bend_triples_d=_FakeArray(503, (n_bends, 3), wp.int32),
        bend_left_segment_d=_FakeArray(504, (n_bends,), wp.int32),
        bend_right_segment_d=_FakeArray(505, (n_bends,), wp.int32),
        bend_active_d=_FakeArray(506, (n_bends,), wp.int32),
        k_seg_d=_FakeArray(507, (n_segments,), wp.float64),
        alpha_d=_FakeArray(508, (n_bends,), wp.float64),
        ea_pn=_EA_PN,
        kappa_pn_um2=_KAPPA,
        link_kernel=_LINK_KERNEL,
        bending_kernel=_BENDING_KERNEL,
        stiffness_kernel=_STIFFNESS_KERNEL,
        alpha_kernel=_ALPHA_KERNEL,
        launch=recorder,
    )


def _owner_arrays(device: _FakeDevice, *, ptr_base: int = 100) -> dict[str, _FakeArray]:
    """The nine CUDA arrays an :class:`ECMStateOwner` owns (8 nodes / 6 segments / 2 fibers)."""
    return {
        "position_d": _FakeArray(ptr_base, (8,), wp.vec3d, device),
        "force_d": _FakeArray(ptr_base + 1, (8,), wp.vec3d, device),
        "segments_d": _FakeArray(ptr_base + 2, (6, 2), wp.int32, device),
        "fiber_offsets_d": _FakeArray(ptr_base + 3, (3,), wp.int32, device),
        "fiber_sleep_state_d": _FakeArray(ptr_base + 4, (2,), wp.int32, device),
        "segment_refinement_level_d": _FakeArray(ptr_base + 5, (6,), wp.int32, device),
        "segment_damage_d": _FakeArray(ptr_base + 6, (6,), wp.float64, device),
        "endpoint_refcount_d": _FakeArray(ptr_base + 7, (6,), wp.int32, device),
        "topology_epoch_d": _FakeArray(ptr_base + 8, (1,), wp.int32, device),
    }


def _owner_with_constitutive_force(
    recorder: _LaunchRecorder,
) -> tuple[ECMStateOwner, ECMConstitutiveForce, _TransactionSpy, _LedgerSpy]:
    device = _FakeDevice()
    arrays = _owner_arrays(device)
    force_pass = _constitutive_force(recorder)  # n_nodes == 8 == owner position capacity
    transaction = _TransactionSpy()
    ledger = _LedgerSpy()
    owner = ECMStateOwner(
        name=ECM_COMPONENT,
        settings=ECMWorldSettings(BoundaryAnchorMode.FAR_FIELD_DIRICHLET),
        actor=ActorRecord("ecm", 7, 2, 13, 3, 6),  # n_elements == 6 == owner segment capacity
        mechanics=force_pass,
        transaction=transaction,
        ledger=ledger,
        **arrays,
    )
    return owner, force_pass, transaction, ledger


# ── owner OWNS the real constitutive force pass and drives it through its own force hook ───────────────
def test_owner_owns_and_launches_the_constitutive_force_over_owned_arrays() -> None:
    recorder = _LaunchRecorder()
    owner, force_pass, _, _ = _owner_with_constitutive_force(recorder)

    # The mechanics slot is the genuine ECMConstitutiveForce, not a spy/facade.
    assert owner.mechanics is force_pass
    assert isinstance(owner.mechanics, ECMConstitutiveForce)

    owner.accumulate()

    # The owner's force hook launched BOTH ff kernels over the owner-owned pos/force arrays.
    assert len(recorder.calls) == 2
    link_call, bend_call = recorder.calls
    assert link_call["kernel"] is _LINK_KERNEL and link_call["dim"] == force_pass.n_segments
    assert bend_call["kernel"] is _BENDING_KERNEL and bend_call["dim"] == force_pass.n_bends
    for call in recorder.calls:
        assert call["inputs"][0] is owner.position_d   # reads owner pos first
        assert call["inputs"][-1] is owner.force_d      # scatters onto owner force last


def test_accepted_commit_refreshes_constitutive_parameters_after_the_transaction() -> None:
    recorder = _LaunchRecorder()
    owner, force_pass, transaction, _ = _owner_with_constitutive_force(recorder)
    accepted = object()

    owner.commit_irreversible(accepted, dt_phys=0.05, rng_seed=19)

    # The accepted-step transaction is forwarded first, THEN the two derivation kernels re-run
    # (k_seg = EA/seg_rest, α = κ/seg³) so a committed remesh never leaves stale parameters.
    assert transaction.calls == [("commit", accepted, 0.05, 19)]
    assert len(recorder.calls) == 2
    k_call, a_call = recorder.calls
    assert k_call["kernel"] is _STIFFNESS_KERNEL and k_call["inputs"][-1] is force_pass.k_seg_d
    assert a_call["kernel"] is _ALPHA_KERNEL and a_call["inputs"][-1] is force_pass.alpha_d


def test_commit_does_not_refresh_inside_accumulate_only_after_commit() -> None:
    recorder = _LaunchRecorder()
    owner, _, _, _ = _owner_with_constitutive_force(recorder)

    owner.accumulate()
    assert all(c["kernel"] in (_LINK_KERNEL, _BENDING_KERNEL) for c in recorder.calls)
    assert not any(c["kernel"] in (_STIFFNESS_KERNEL, _ALPHA_KERNEL) for c in recorder.calls)


# ── additivity: a delegate without refresh_parameters is untouched (existing components unchanged) ─────
def test_refresh_is_a_no_op_for_a_delegate_without_the_hook() -> None:
    device = _FakeDevice()
    arrays = _owner_arrays(device)
    transaction = _TransactionSpy()
    plain = _PlainMechanicsSpy(n_nodes=8)
    owner = ECMStateOwner(
        name=ECM_COMPONENT,
        settings=ECMWorldSettings(BoundaryAnchorMode.FAR_FIELD_DIRICHLET),
        actor=ActorRecord("ecm", 7, 2, 13, 3, 6),
        mechanics=plain,
        transaction=transaction,
        ledger=_LedgerSpy(),
        **arrays,
    )
    accepted = object()

    # commit still forwards to the transaction; refresh_mechanics_parameters is a silent no-op.
    owner.commit_irreversible(accepted, dt_phys=0.1, rng_seed=3)
    owner.refresh_mechanics_parameters()
    assert transaction.calls == [("commit", accepted, 0.1, 3)]
    assert plain.calls == []  # never launched anything for a build-once delegate


# ── the builder constructs the force pass from the topology + card and wires the mechanics seam ────────
@dataclass(frozen=True, slots=True)
class _FakeTopology:
    """Minimal ECMTopologyState double exposing only the arrays the builder maps into the owner."""

    position_d: _FakeArray
    force_d: _FakeArray
    segments_d: _FakeArray
    fiber_offsets_d: _FakeArray
    fiber_sleep_state_d: _FakeArray
    segment_refinement_level_d: _FakeArray
    segment_damage_d: _FakeArray
    endpoint_refcount_d: _FakeArray
    topology_epoch_d: _FakeArray


def test_builder_constructs_force_pass_from_topology_and_wires_owner_mechanics() -> None:
    device = _FakeDevice()
    arrays = _owner_arrays(device)
    topology = _FakeTopology(**arrays)
    recorder = _LaunchRecorder()
    force_pass = _constitutive_force(recorder)

    seen: dict[str, object] = {}

    def _fake_force_builder(topo, card, *, device=None, launch=None) -> ECMConstitutiveForce:
        seen["topology"] = topo
        seen["card"] = card
        return force_pass

    owner = build_ecm_state_owner(
        topology,
        settings=ECMWorldSettings(BoundaryAnchorMode.FAR_FIELD_DIRICHLET),
        actor=ActorRecord("ecm", 7, 2, 13, 3, 6),
        transaction=_TransactionSpy(),
        ledger=_LedgerSpy(),
        force_builder=_fake_force_builder,
    )

    # The builder fed the topology + defaulted the SOURCED collagen card to the force factory,
    assert seen["topology"] is topology
    assert seen["card"] == "collagen_I"
    # wired the real force pass into the mechanics seam,
    assert owner.mechanics is force_pass
    # and mapped every owned array straight from the topology SoA (no copy).
    assert owner.position_d is topology.position_d
    assert owner.force_d is topology.force_d
    assert owner.segments_d is topology.segments_d
    assert owner.topology_epoch_d is topology.topology_epoch_d
    assert owner.name == ECM_COMPONENT
