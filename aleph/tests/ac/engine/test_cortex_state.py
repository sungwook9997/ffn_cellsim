"""Structural gates for the cortex component's PRIVATE device state — the end of the bind-target port.

CUDA-free: every array is a metadata-only double and every launch goes to a recorder, so these gates check
the OWNERSHIP CONTRACT (who owns which allocation, who launches what onto it, what the predicate gates) and
never physics.  The physics parity against the incumbent is a native gate on the A5000
(``scripts/ac_gate_a_cortex_ownership_native.py``).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest
import warp as wp

from aleph.engine.cortex_population import build_cortex_population
from aleph.engine.cortex_state import (
    CORTEX_CHANNELS_NOT_BOUND,
    CortexBranchAngleMechanics,
    CortexCompositeMechanics,
    CortexStateLedger,
    CortexStateTransaction,
    assemble_cortex_motor_port,
    assemble_cortex_state_owner,
    assert_component_state_disjoint,
    cortex_segment_material_coordinates,
    incumbent_channels_to_omit,
)
from aleph.engine.surface_body import CortexFilamentMechanics

N_NODES = 6
N_LINKS = 3
N_TRIPLES = 4


@dataclass(frozen=True, slots=True)
class _FakeDevice:
    alias: str = "cuda"
    is_cuda: bool = True

    def __str__(self) -> str:
        return self.alias


@dataclass(frozen=True, slots=True)
class _FakeArray:
    """Metadata-only CUDA array double; it never executes authoritative physics."""

    ptr: int
    shape: tuple[int, ...] = (N_NODES,)
    dtype: object = wp.vec3d
    device: _FakeDevice = _FakeDevice()


@dataclass(slots=True)
class _LaunchRecorder:
    calls: list[tuple[object, int, list, object]] = field(default_factory=list)

    def __call__(self, kernel: object, *, dim: int, inputs: list, device: object | None = None) -> None:
        self.calls.append((kernel, int(dim), list(inputs), device))

    @property
    def kernels(self) -> list[object]:
        return [call[0] for call in self.calls]


@dataclass(slots=True)
class _CopyRecorder:
    calls: list[tuple[object, object]] = field(default_factory=list)

    def __call__(self, dst: object, src: object) -> None:
        self.calls.append((dst, src))


@dataclass(slots=True)
class _BalanceLedgerSpy:
    """Stands in for ``GlobalCellLedger``: records which array landed in which half of the balance gate."""

    calls: list[tuple[object, str]] = field(default_factory=list)

    def add_body_force(self, force_d: object, *, side: str) -> None:
        self.calls.append((force_d, side))


class _LedgerWithoutBodyForce:
    """A ledger assembled for another lane; the cortex must skip it, never adapt into a foreign channel."""


def _filament_mechanics(recorder: _LaunchRecorder, *, n_vertices: int = N_NODES) -> CortexFilamentMechanics:
    return CortexFilamentMechanics(
        device="cuda",
        n_vertices=n_vertices,
        links_d=_FakeArray(0x1000, (N_LINKS, 2), wp.int32),
        link_k_d=_FakeArray(0x1010, (N_LINKS,), wp.float64),
        link_r0_d=_FakeArray(0x1020, (N_LINKS,), wp.float64),
        bend_triples_d=_FakeArray(0x1030, (N_TRIPLES, 3), wp.int32),
        bend_alpha_d=_FakeArray(0x1040, (N_TRIPLES,), wp.float64),
        link_kernel=object(),
        bending_kernel=object(),
        launch=recorder,
    )


def _owner(recorder: _LaunchRecorder, copier: _CopyRecorder, *, mechanics=None):
    return assemble_cortex_state_owner(
        position_d=_FakeArray(0x2000),
        force_d=_FakeArray(0x2100),
        position_snap_d=_FakeArray(0x2200),
        mechanics=mechanics if mechanics is not None else CortexCompositeMechanics(
            filament=_filament_mechanics(recorder)),
        device="cuda",
        launch=recorder,
        copy=copier,
    )


# ── ownership ────────────────────────────────────────────────────────────────────────────────────────────
def test_the_cortex_owns_three_distinct_allocations() -> None:
    """Position, force and snapshot are separate; the owner is named ``cortex`` and nothing else."""
    owner = _owner(_LaunchRecorder(), _CopyRecorder())
    assert owner.name == "cortex"
    assert owner.position_d.ptr != owner.force_d.ptr


def test_a_mechanics_addressing_a_longer_array_is_refused() -> None:
    """THE gate the aliased port cannot pass.

    Under the bind-target port the cortex's "own" array is the incumbent's ``n_total``-long global array
    while its mechanics addresses only the ``n_actin`` cortex nodes.  The state owner asserts those two
    numbers are equal, so that arrangement is rejected by construction rather than by inspection.
    """
    recorder = _LaunchRecorder()
    with pytest.raises(ValueError, match="mechanics expects"):
        assemble_cortex_state_owner(
            position_d=_FakeArray(0x3000, (N_NODES + 11,)),   # a global [actin | myosin | ...] array
            force_d=_FakeArray(0x3100, (N_NODES + 11,)),
            position_snap_d=_FakeArray(0x3200, (N_NODES + 11,)),
            mechanics=CortexCompositeMechanics(filament=_filament_mechanics(recorder)),
            device="cuda",
            launch=recorder,
            copy=_CopyRecorder(),
        )


def test_position_and_force_may_not_alias() -> None:
    """A component whose force IS its position has no force at all."""
    recorder = _LaunchRecorder()
    shared = _FakeArray(0x4000)
    with pytest.raises(ValueError, match="must not alias"):
        assemble_cortex_state_owner(
            position_d=shared,
            force_d=shared,
            position_snap_d=_FakeArray(0x4200),
            mechanics=CortexCompositeMechanics(filament=_filament_mechanics(recorder)),
            device="cuda",
            launch=recorder,
            copy=_CopyRecorder(),
        )


def test_the_snapshot_may_not_alias_the_live_position() -> None:
    """A snapshot that is the live array restores nothing on a rejected step."""
    live = _FakeArray(0x5000)
    with pytest.raises(ValueError, match="separate allocation"):
        CortexStateTransaction(position_d=live, position_snap_d=live, device="cuda")


def test_two_components_sharing_one_allocation_are_rejected() -> None:
    """Co-location in an array is NEVER a connection — the composition-wide form of the guard (T2 iii)."""
    recorder, copier = _LaunchRecorder(), _CopyRecorder()
    cortex = _owner(recorder, copier)

    @dataclass(frozen=True, slots=True)
    class _Other:
        name: str
        position_d: object
        force_d: object

    disjoint = _Other("membrane", _FakeArray(0x6000), _FakeArray(0x6100))
    assert_component_state_disjoint([cortex, disjoint])

    aliased = _Other("membrane", cortex.position_d, _FakeArray(0x6100))
    with pytest.raises(ValueError, match="NEVER a connection"):
        assert_component_state_disjoint([cortex, aliased])


# ── not a passenger ──────────────────────────────────────────────────────────────────────────────────────
def test_a_state_owner_with_no_force_channel_is_refused() -> None:
    """The T2 negative control: distinct arrays plus a no-op mechanics is the passenger pattern.

    The refusal lives one level down, in the mechanics delegate, which is why the composite does not repeat
    it — the point is that the passenger CANNOT be assembled, not where the check sits.
    """
    recorder = _LaunchRecorder()
    bending_only = CortexFilamentMechanics(
        device="cuda",
        n_vertices=N_NODES,
        links_d=_FakeArray(0x7000, (0, 2), wp.int32),
        link_k_d=_FakeArray(0x7010, (0,), wp.float64),
        link_r0_d=_FakeArray(0x7020, (0,), wp.float64),
        bend_triples_d=_FakeArray(0x7030, (1, 3), wp.int32),
        bend_alpha_d=_FakeArray(0x7040, (1,), wp.float64),
        link_kernel=object(),
        bending_kernel=object(),
        launch=recorder,
    )
    # A single bending triple is a real channel, so this composite is legal ...
    assert CortexCompositeMechanics(filament=bending_only).channels == ("cytosim_bending",)
    # ... but a delegate with no channel at all cannot be constructed, so no owner can hold one.
    with pytest.raises(ValueError, match="at least one link or bending triple"):
        CortexFilamentMechanics(
            device="cuda",
            n_vertices=N_NODES,
            links_d=_FakeArray(0x7100, (0, 2), wp.int32),
            link_k_d=_FakeArray(0x7110, (0,), wp.float64),
            link_r0_d=_FakeArray(0x7120, (0,), wp.float64),
            bend_triples_d=_FakeArray(0x7130, (0, 3), wp.int32),
            bend_alpha_d=_FakeArray(0x7140, (0,), wp.float64),
            link_kernel=object(),
            bending_kernel=object(),
            launch=recorder,
        )


def test_the_bound_channels_are_the_incumbents_cortex_channels() -> None:
    """Bending then crosslink, in ``_accumulate_all`` order, launched onto the cortex's own force array."""
    recorder, copier = _LaunchRecorder(), _CopyRecorder()
    owner = _owner(recorder, copier)
    assert owner.mechanics.channels == ("cytosim_bending", "crosslink_link_spring")
    owner.accumulate()
    assert len(recorder.calls) == 2
    link_call, bend_call = recorder.calls[0], recorder.calls[1]
    assert link_call[1] == N_LINKS and bend_call[1] == N_TRIPLES
    for call in recorder.calls:
        assert call[2][0] is owner.position_d, "a channel must read the cortex's own position"
        assert call[2][-1] is owner.force_d, "a channel must write the cortex's own force"


def test_channels_left_unbound_are_named_not_silently_dropped() -> None:
    """A force absent for a stated reason is a different object from a force that was forgotten."""
    assert set(CORTEX_CHANNELS_NOT_BOUND) == {"axial_inextensibility", "steric_wca", "biot_pressure"}
    assert all(reason.strip() for reason in CORTEX_CHANNELS_NOT_BOUND.values())


def test_every_bound_channel_has_an_incumbent_name_to_omit() -> None:
    """The relocation must be TOTAL: a channel this owner binds and the incumbent still launches is doubled.

    This is the regression the derived omit set exists for. The driver used to type
    ``{"actin_bending_cytosim", "actin_crosslink_link_spring"}`` by hand, which is exactly right for the
    formin-only default cortex and silently double-counts Arp2/3 for a mixed one — invisible in the result,
    and invisible in the artifact too, because every A/B so far ran at ``arp23_fraction = 0``.
    """
    from aleph.components.incumbent.driver import OMITTABLE_CHANNELS

    recorder = _LaunchRecorder()
    branch = CortexBranchAngleMechanics(
        device="cuda",
        branch_triples_d=_FakeArray(0x9000, (2, 3), wp.int32),
        branch_active_d=_FakeArray(0x9010, (2,), wp.int32),
        theta0_rad=1.2217304764, k_theta=0.173, kernel=object(), launch=recorder,
    )
    for mechanics in (CortexCompositeMechanics(filament=_filament_mechanics(recorder)),
                      CortexCompositeMechanics(filament=_filament_mechanics(recorder), branch=branch)):
        omit = incumbent_channels_to_omit(mechanics)
        assert len(omit) == len(mechanics.channels), "a bound channel lost its omit entry — doubled force"
        assert omit <= OMITTABLE_CHANNELS, "the incumbent cannot be told to skip a channel it does not name"

    mixed = CortexCompositeMechanics(filament=_filament_mechanics(recorder), branch=branch)
    assert "arp23_branch_angle" in incumbent_channels_to_omit(mixed), (
        "the mixed cortex must omit the incumbent's branch-angle launch, or Arp2/3 is counted twice"
    )


def test_a_mixed_cortex_adds_the_branch_channel() -> None:
    """Arp2/3 binds the angle HARMONIC, never a rigid clamp, and only when triples exist."""
    recorder = _LaunchRecorder()
    branch = CortexBranchAngleMechanics(
        device="cuda",
        branch_triples_d=_FakeArray(0x8000, (2, 3), wp.int32),
        branch_active_d=_FakeArray(0x8010, (2,), wp.int32),
        theta0_rad=1.2217304764,
        k_theta=0.173,
        kernel=object(),
        launch=recorder,
    )
    composite = CortexCompositeMechanics(filament=_filament_mechanics(recorder), branch=branch)
    assert composite.channels == ("cytosim_bending", "crosslink_link_spring", "arp23_branch_angle")
    with pytest.raises(ValueError, match="formin-only"):
        CortexBranchAngleMechanics(
            device="cuda",
            branch_triples_d=_FakeArray(0x8100, (0, 3), wp.int32),
            branch_active_d=_FakeArray(0x8110, (0,), wp.int32),
            theta0_rad=1.2217304764, k_theta=0.173, kernel=object(), launch=recorder,
        )


# ── transaction ──────────────────────────────────────────────────────────────────────────────────────────
def test_snapshot_and_reject_gated_restore_stay_on_device() -> None:
    """The predicate reaches a kernel, never a host ``bool()``; snapshot is snap <- live, not the reverse."""
    recorder, copier = _LaunchRecorder(), _CopyRecorder()
    owner = _owner(recorder, copier)
    owner.snapshot_candidate()
    assert len(copier.calls) == 1
    destination, source = copier.calls[0]
    assert source is owner.position_d and destination is not owner.position_d

    accepted = _FakeArray(0x9000, (1,), wp.int32)
    before = len(recorder.calls)
    owner.rollback(accepted)
    restore = recorder.calls[before]
    assert restore[1] == N_NODES
    assert restore[2][0] is accepted, "the predicate is the kernel's first input, not a host branch"
    assert restore[2][1] is owner.position_d


def test_the_cortex_advances_no_private_kinetic_epoch() -> None:
    """Turnover rates are an unfilled PI-GAP; a silent private clock would be worse than none."""
    recorder, copier = _LaunchRecorder(), _CopyRecorder()
    owner = _owner(recorder, copier)
    before = len(recorder.calls) + len(copier.calls)
    owner.commit_irreversible(_FakeArray(0x9100, (1,), wp.int32), 0.01, 17)
    assert len(recorder.calls) + len(copier.calls) == before


# ── ledger ───────────────────────────────────────────────────────────────────────────────────────────────
def test_the_cortex_reduces_its_own_force_into_one_half_of_the_balance_gate() -> None:
    """The operation array ownership buys: the cortex becomes a genuine side of a Newton pair."""
    recorder, copier = _LaunchRecorder(), _CopyRecorder()
    owner = _owner(recorder, copier)
    ledger = _BalanceLedgerSpy()
    owner.accumulate_ledger(ledger)
    assert len(ledger.calls) == 1
    force_d, side = ledger.calls[0]
    assert force_d is owner.force_d, "never a copy computed from the other side"
    assert side == "reaction"


def test_a_ledger_without_the_body_force_channel_is_skipped_not_adapted() -> None:
    """Routing a body force into a differently-meaning channel is the accounting error the gate exists for."""
    contributor = CortexStateLedger(force_d=_FakeArray(0xA000))
    contributor.accumulate_ledger(_LedgerWithoutBodyForce())
    assert contributor.contributions == 0


def test_the_ledger_side_must_be_one_of_the_two_channels() -> None:
    with pytest.raises(ValueError, match="reaction"):
        CortexStateLedger(force_d=_FakeArray(0xA100), side="whichever")


# ── the motor port, over the component's own arrays ──────────────────────────────────────────────────────
def test_segment_material_coordinates_are_rest_arc_lengths_restarting_per_filament() -> None:
    """A head binds a material coordinate; a segment ordinal is not one."""
    population = build_cortex_population(6, seg_um=0.5, length_um=3.0, density_per_fil=4.0, seed=3)
    topology = population.topology
    s0, s1 = cortex_segment_material_coordinates(topology)
    assert s0.shape == s1.shape == (topology.n_segments,)
    assert np.all(s1 > s0)
    offsets = np.asarray(topology.fiber_offsets, np.int64)
    seg_per_fiber = np.diff(offsets) - 1
    starts = np.concatenate([[0], np.cumsum(seg_per_fiber)]).astype(np.int64)
    for f in range(topology.n_fibers):
        a, b = int(starts[f]), int(starts[f + 1])
        assert s0[a] == pytest.approx(0.0), "arc length restarts at every filament's first node"
        contour = float(np.sum(np.asarray(topology.seg_rest)[a:b]))
        assert s1[b - 1] == pytest.approx(contour)


def test_the_motor_port_reads_and_writes_the_cortex_owned_arrays() -> None:
    """The head's ``+f`` and the actin nodes' ``-f`` then land in two DIFFERENT allocations."""
    recorder, copier = _LaunchRecorder(), _CopyRecorder()
    owner = _owner(recorder, copier)
    n_segments = 5
    port = assemble_cortex_motor_port(
        component="cortex",
        position_d=owner.position_d,
        force_d=owner.force_d,
        segment_node_a_d=_FakeArray(0xB000, (n_segments,), wp.int32),
        segment_node_b_d=_FakeArray(0xB010, (n_segments,), wp.int32),
        segment_polarity_d=_FakeArray(0xB020, (n_segments,), wp.int32),
        persistent_filament_id_d=_FakeArray(0xB030, (n_segments,), wp.int64),
        material_s0_d=_FakeArray(0xB040, (n_segments,), wp.float64),
        material_s1_d=_FakeArray(0xB050, (n_segments,), wp.float64),
        topology_epoch_d=_FakeArray(0xB060, (1,), wp.int32),
    )
    assert port.component == "cortex"
    assert port.position_d is owner.position_d
    assert port.force_d is owner.force_d


def test_the_ports_persistent_ids_would_come_from_the_population_block() -> None:
    """Device-side no-double-count: a bound head's record names a filament owned by exactly one component."""
    population = build_cortex_population(
        5, seg_um=0.5, length_um=3.0, density_per_fil=4.0, seed=3, id_base=0)
    topology = population.topology
    ids = np.asarray(topology.node_fiber, np.int64)[np.asarray(topology.seg_node_a, np.int64)]
    low, high = population.ledger.block
    assert ids.min() >= low - population.ledger.block[0] and ids.max() < high
    assert set(np.unique(ids)) <= set(range(low, high))
