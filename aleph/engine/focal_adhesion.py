"""The ``focal_adhesion`` component: a geometry-less stochastic clutch graph over its OWN device state.

WHAT THIS CLOSES.  ``focal_adhesion`` was the one component in
:func:`~aleph.engine.contracts.reference_cell_architecture` with no state owner at all — every other
declared component had one, so a composed world could never reach ``require_complete=True`` however many
other slices landed.  Its contract is::

    ComponentContract("focal_adhesion", ACTIVE_LOAD_PATH, "stochastic clutch graph", "joint KMC",
                      owns_geometry=False)

``owns_geometry=False`` is the whole design.  This owner holds **no positions and no force array**: the
mechanical FA spring is the composite :class:`~aleph.engine.load_path.LoadPathJointRuntime` that
:class:`~aleph.engine.actor.CellActor` binds ONCE under the ``alpha2beta1_collagen_series`` mechanical
group and registers under BOTH semantic edges (``fa_actin_anchor`` and ``integrin_collagen_clutch``).
What this owner owns is the part that is not a spring — the per-clutch **bound state**, its age, and its
RNG epoch — i.e. the "joint KMC" half of the contract.  Splitting it the other way (an owner that also
carries a force array) would put the same spring in two objects and turn one series clutch into two in
parallel, which is exactly the failure the composite mechanical group exists to prevent.

WHY THE STATE IS DOUBLE-BUFFERED.  A clutch bind/unbind is a kinetic event, and the charter admits a
kinetic transition ONLY on an accepted physical step, under one device-resident transaction.  So the KMC
proposal writes ``candidate_state_d`` and never ``state_d``; :meth:`commit_irreversible` promotes the
candidate and ages the bond under the device predicate, and :meth:`rollback` discards it.  Neither reads
the predicate on the host.  A single-buffered design cannot express "proposed but not yet accepted" and
would let a rejected step leave a bond bound.

Sanity Gate (before first execution; :mod:`aleph.tests.ac.engine.test_focal_adhesion`):
  * dimensions: ``age_d`` seconds, ``dt_phys`` seconds; every other owned array is a dimensionless index,
    enum or counter.  This owner carries NO force, length or stiffness — those live on the joint runtime,
    in pN, um and pN/um respectively.
  * boundary cases: ``capacity == 0`` is legal (a cell with no engaged adhesion) and every hook is then a
    no-op launch of dimension zero, not a raise; ``active_d == 0`` masks a slot out of both KMC and ageing.
  * conservation invariant: the clutch population is FIXED at build.  ``commit_irreversible`` may change a
    slot's STATE, never the slot count, so the FA population ledger block cannot grow into another
    component's ID range.  Asserted per step by the composed world's disjoint-population re-assert.
  * accepted-step invariant: after a REJECTED step every owned array is bit-identical to the snapshot;
    after an ACCEPTED step ``state_d == candidate_state_d`` and ``age_d`` advanced by exactly ``dt_phys``
    on bound slots only.  Both are the test's assertions, driven on recording doubles.
  * precision: state/epoch are ``int32`` (a step count, never a physical magnitude); ``age_d`` is
    ``float64`` because it accumulates ``dt_phys`` over a run and a ``float32`` bond age loses the last
    digit within one physical second at the production ``dt``.
  * sign sense: ``age_d`` is monotone non-decreasing under accepted steps and is RESET to zero on any
    slot whose committed state differs from the state it held at snapshot — a bond that just formed is
    zero seconds old, not the age of the bond it replaced.
  * measurement protocol: nothing here is a physics magnitude.  The only quantity this owner can be asked
    for is a COUNT (bound slots per state), read between accepted steps, never inside the loop.

This module allocates device memory but launches nothing at import, and every launch goes through the
injected ``launch``/``copy`` seam, so it is CPU-importable and the structural gate drives it with
recording doubles exactly as :class:`~aleph.engine.composed_native.SFArcStateOwner` is driven.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import warp as wp

from aleph.engine.load_path import JointState

FA_COMPONENT = "focal_adhesion"

#: The two semantic edges the ONE composite series joint is registered under (``ConnectorContract``
#: ``mechanical_group="alpha2beta1_collagen_series"``).  Named here so a driver cannot bind a second
#: runtime to one of them by accident; :meth:`CellActor.bind_connector` is what enforces it.
FA_COMPOSITE_EDGES = ("fa_actin_anchor", "integrin_collagen_clutch")


@wp.kernel
def _restore_int32_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.int32),
    snap: wp.array(dtype=wp.int32),
) -> None:
    """Restore ``live <- snap`` iff the outer step was REJECTED (``accepted[0] == 0``; device-only branch)."""
    if accepted[0] == wp.int32(0):
        t = wp.tid()
        live[t] = snap[t]


@wp.kernel
def _restore_float64_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.float64),
    snap: wp.array(dtype=wp.float64),
) -> None:
    """Restore ``live <- snap`` iff the outer step was REJECTED (``accepted[0] == 0``; device-only branch)."""
    if accepted[0] == wp.int32(0):
        t = wp.tid()
        live[t] = snap[t]


@wp.kernel
def _commit_clutch_if_accepted_kernel(
    accepted: wp.array(dtype=wp.int32),
    dt_phys: wp.float64,
    active: wp.array(dtype=wp.int32),
    state: wp.array(dtype=wp.int32),
    candidate_state: wp.array(dtype=wp.int32),
    state_snap: wp.array(dtype=wp.int32),
    age: wp.array(dtype=wp.float64),
) -> None:
    """Promote the proposed clutch state and age the bond, iff the outer step was ACCEPTED.

    A slot whose committed state differs from the state it held at snapshot is a bond that changed in
    THIS step, so its age restarts at zero rather than inheriting the previous bond's age.  Inactive
    slots are left untouched: they carry no bond to age and no transition to commit.
    """
    if accepted[0] != wp.int32(0):
        t = wp.tid()
        if active[t] != wp.int32(0):
            proposed = candidate_state[t]
            if proposed != state_snap[t]:
                state[t] = proposed
                age[t] = wp.float64(0.0)
            else:
                state[t] = proposed
                age[t] = age[t] + dt_phys


@wp.kernel
def _bump_epoch_if_accepted_kernel(
    accepted: wp.array(dtype=wp.int32),
    epoch: wp.array(dtype=wp.int32),
) -> None:
    """Advance the accepted-step epoch iff ``accepted[0] != 0`` (device-only; no host predicate read)."""
    if accepted[0] != wp.int32(0):
        epoch[0] = epoch[0] + wp.int32(1)


@dataclass(slots=True)
class FocalAdhesionStateOwner:
    """The ``focal_adhesion`` component as a composed-world transaction participant.

    Owns the per-clutch KMC state and NOTHING mechanical — see the module docstring for why the spring
    stays on the composite joint runtime.  Implements the three transaction hooks
    :class:`~aleph.engine.world.CellWorldTransaction` requires (all-or-nothing: a partial transaction API
    is a ``TypeError`` at world construction, not a silent skip).

    It deliberately does NOT implement ``accumulate_ledger``.  That hook reduces a body's own force-array
    resultant into the balance gate, and this body has no force array; contributing a zero would make the
    gate's two channels look balanced for a reason that has nothing to do with the adjoint wiring under
    test.  The FA load path enters the ledger through the composite joint's two endpoint bodies.
    """

    name: str
    device: str
    capacity: int
    active_d: wp.array
    cluster_id_d: wp.array
    state_d: wp.array
    candidate_state_d: wp.array
    age_d: wp.array
    rng_epoch_d: wp.array
    state_snap_d: wp.array
    candidate_state_snap_d: wp.array
    age_snap_d: wp.array
    rng_epoch_snap_d: wp.array
    epoch_d: wp.array
    launch: Callable[..., object] = wp.launch
    copy: Callable[..., object] = wp.copy

    def __post_init__(self) -> None:
        if self.name != FA_COMPONENT:
            raise ValueError(f"focal-adhesion state owner must be {FA_COMPONENT!r}")
        if self.capacity < 0:
            raise ValueError("focal-adhesion clutch capacity must be nonnegative")

    def snapshot_candidate(self) -> None:
        """Snapshot every owned array before candidate kinetics (device D2D copies, no host read)."""
        if not self.capacity:
            return
        self.copy(self.state_snap_d, self.state_d)
        self.copy(self.candidate_state_snap_d, self.candidate_state_d)
        self.copy(self.age_snap_d, self.age_d)
        self.copy(self.rng_epoch_snap_d, self.rng_epoch_d)

    def rollback(self, accepted: wp.array) -> None:
        """Restore the rejected candidate's clutch state in-device (accepted steps keep the proposal)."""
        if not self.capacity:
            return
        for live, snap in (
            (self.state_d, self.state_snap_d),
            (self.candidate_state_d, self.candidate_state_snap_d),
            (self.rng_epoch_d, self.rng_epoch_snap_d),
        ):
            self.launch(
                _restore_int32_if_rejected_kernel, dim=self.capacity,
                inputs=[accepted, live, snap], device=self.device,
            )
        self.launch(
            _restore_float64_if_rejected_kernel, dim=self.capacity,
            inputs=[accepted, self.age_d, self.age_snap_d], device=self.device,
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Promote the proposed clutch states and age the bonds under the final device predicate."""
        if self.capacity:
            self.launch(
                _commit_clutch_if_accepted_kernel, dim=self.capacity,
                inputs=[
                    accepted, wp.float64(float(dt_phys)), self.active_d, self.state_d,
                    self.candidate_state_d, self.state_snap_d, self.age_d,
                ],
                device=self.device,
            )
        self.launch(
            _bump_epoch_if_accepted_kernel, dim=1,
            inputs=[accepted, self.epoch_d], device=self.device,
        )

    def bound_counts(self) -> dict[str, int]:
        """OUT-OF-LOOP diagnostic: clutch slots per state.  Never call inside the physical-time loop."""
        if not self.capacity:
            return {state.name: 0 for state in JointState}
        states = self.state_d.numpy()
        active = self.active_d.numpy()
        return {
            state.name: int(np.count_nonzero((states == int(state)) & (active != 0)))
            for state in JointState
        }


def build_focal_adhesion_state_owner(
    *,
    cluster_ids: np.ndarray,
    device: str,
    initial_state: JointState = JointState.ACTIN_ENGAGED,
    launch: Callable[..., object] = wp.launch,
    copy: Callable[..., object] = wp.copy,
) -> FocalAdhesionStateOwner:
    """Allocate the FA clutch graph's PRIVATE device state, one slot per composite series joint.

    The clutch population is fixed at build and every slot starts ACTIVE: a slot is the *existence* of a
    potential adhesion, and its bound/unbound history is the state, not the allocation.  Deallocating on
    unbind would make the population non-conserved and break the disjoint-ID re-assert.

    Args:
        cluster_ids: FA cluster ID per clutch slot, shape ``(capacity,)``.  One composite
            :class:`~aleph.engine.load_path.CompositeFAJointSpec` per slot, in the same order.
        device: the resolved CUDA device string.
        initial_state: the state every slot starts in.  Defaults to the same
            :class:`~aleph.engine.load_path.JointState` the composite joint spec defaults to, so the owner
            and its joint cannot disagree at step zero.
        launch: injected kernel launcher (the structural gate passes a recording double).
        copy: injected device-copy (likewise).

    Returns:
        A :class:`FocalAdhesionStateOwner` ready to bind at the ``focal_adhesion`` component slot.

    Raises:
        ValueError: if ``cluster_ids`` is not a 1-D nonnegative integer array.
    """
    ids = np.asarray(cluster_ids)
    if ids.ndim != 1:
        raise ValueError("cluster_ids must be a 1-D array, one entry per clutch slot")
    if ids.size and (not np.issubdtype(ids.dtype, np.integer) or int(ids.min()) < 0):
        raise ValueError("cluster_ids must be nonnegative integers")
    capacity = int(ids.size)

    def i32(fill: int) -> wp.array:
        return wp.array(np.full(capacity, fill, dtype=np.int32), dtype=wp.int32, device=device)

    def f64(fill: float) -> wp.array:
        return wp.array(np.full(capacity, fill, dtype=np.float64), dtype=wp.float64, device=device)

    state = int(initial_state)
    return FocalAdhesionStateOwner(
        name=FA_COMPONENT,
        device=device,
        capacity=capacity,
        active_d=i32(1),
        cluster_id_d=wp.array(ids.astype(np.int32), dtype=wp.int32, device=device),
        state_d=i32(state),
        candidate_state_d=i32(state),
        age_d=f64(0.0),
        rng_epoch_d=i32(0),
        state_snap_d=i32(state),
        candidate_state_snap_d=i32(state),
        age_snap_d=f64(0.0),
        rng_epoch_snap_d=i32(0),
        epoch_d=wp.zeros(1, dtype=wp.int32, device=device),
        launch=launch,
        copy=copy,
    )
