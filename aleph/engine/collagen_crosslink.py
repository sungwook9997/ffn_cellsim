"""``ecm_crosslink``: the collagen--collagen bond population, kinetic and internal to the ECM.

WHAT THIS CLOSES.  :class:`~aleph.engine.ecm_world.ECMInternalCrosslink` is a Protocol that
:class:`~aleph.engine.ecm_world.ECMWorld` already holds a field for and already validates — the name, the
internal endpoints and the adjoint requirement are all checked at world construction.  What did not exist
was anything satisfying it, which is why the tier-(a) native ECM row reads ``zero crosslink bonds``: the
constitutive fibre force was measured on a network whose bonds were declared and never built.

INTERNAL, so ONE array.  ``ecm_crosslink`` is ``ecm -> ecm``: both endpoints are ECM-owned nodes in the
single ECM position array, so the bond binds a ``(J, 2)`` pair table directly rather than the two-array
adjoint runtime the inter-component crosslinks need
(:mod:`aleph.engine.filament_crosslink`).  Equal-and-opposite is then automatic — the same kernel adds
``+f`` at one end and ``-f`` at the other of the same array.

KINETIC, unlike the other internal joint in the tree.  ``dorsal_arc_crosslink`` is a fixed internal
connector; ``ecm_crosslink`` declares ``kinetics=True`` and ``commit_on_accept=True``, so bonds associate
and dissociate.  That is why this cannot reuse :class:`~aleph.engine.sf_mechanics.SFInternalArcJointConnector`
verbatim: it needs the double-buffered bond state that only commits on an accepted physical step, the
same discipline :class:`~aleph.engine.focal_adhesion.FocalAdhesionStateOwner` carries for clutches.

WHY THE FORCE KERNEL IS GATED RATHER THAN THE BOND LIST COMPACTED.  A dissociated crosslink must carry no
force, and the two ways to arrange that are to drop it from the launch or to zero its contribution.  The
population is FIXED at build here and the kernel skips unbound slots, because compaction would change the
bond count between steps and make the ECM population ledger non-conserved — the same reason the FA clutch
graph keeps its slots.  ``proposals`` therefore write ``candidate_bound``; nothing rewrites the table.

Sanity Gate (before first execution; :mod:`aleph.tests.ac.engine.test_collagen_crosslink`):
  * dimensions: stiffness pN/um, rest length um, ``age`` seconds.  No constant is introduced here; the
    collagen crosslink stiffness arrives on the builder and is a chemistry-card quantity, not a default.
  * boundary cases: zero bonds is legal (an uncrosslinked mikado network is a real configuration and is
    what the tier-(a) row measured) and every hook is then a no-op; a bond whose two endpoints are the
    same node is refused at build, since a self-bond has no direction and would divide by zero.
  * conservation invariant: the population is fixed — commit changes a bond's STATE, never the slot
    count.  Equal-and-opposite holds identically because both scatters are into one array.
  * accepted-step invariant: after a REJECTED step bond state and age are bit-identical to the snapshot;
    after an ACCEPTED step ``bound == candidate_bound`` and ages advance on bound slots only, with a
    slot that just changed state restarting at zero.
  * sign sense: a stretched bond pulls its two nodes together; an unbound bond (``bound == 0``) is
    skipped entirely, so a dissociated crosslink cannot pull.  A bond exactly at rest contributes zero.
  * precision: state/epoch ``int32``; age and the mechanical scalars ``float64``, matching every other
    joint SoA in the tree.
  * measurement protocol: the only quantity this connector can be asked for is a bond COUNT, read
    between accepted steps.  It produces no magnitude.

CPU-importable: allocation happens in the builder and every launch goes through the injected seam.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import warp as wp

from aleph.engine.ecm_world import ECM_COMPONENT, ECM_CROSSLINK_CONNECTOR

__all__ = ["CollagenCrosslinkConnector", "build_collagen_crosslink_connector"]


@wp.kernel
def _bound_link_spring_kernel(
    pos: wp.array(dtype=wp.vec3d),
    bonds: wp.array(dtype=wp.int32, ndim=2),
    bound: wp.array(dtype=wp.int32),
    k: wp.array(dtype=wp.float64),
    r0: wp.array(dtype=wp.float64),
    force: wp.array(dtype=wp.vec3d),
) -> None:
    """State-gated Hookean bond over ONE array: a dissociated crosslink contributes nothing.

    Equal-and-opposite is identical here rather than adjoint-by-construction: both scatters land in the
    same force array, so the pair sums to zero exactly whatever the configuration.
    """
    t = wp.tid()
    if bound[t] == wp.int32(0):
        return
    i = bonds[t, 0]
    j = bonds[t, 1]
    d = pos[j] - pos[i]
    length = wp.length(d)
    if length <= wp.float64(0.0):
        return
    f = k[t] * (length - r0[t]) / length * d
    wp.atomic_add(force, i, f)
    wp.atomic_add(force, j, -f)


@wp.kernel
def _restore_int32_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.int32),
    snap: wp.array(dtype=wp.int32),
) -> None:
    """Restore ``live <- snap`` iff the outer step was REJECTED (device-only branch)."""
    if accepted[0] == wp.int32(0):
        t = wp.tid()
        live[t] = snap[t]


@wp.kernel
def _restore_float64_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.float64),
    snap: wp.array(dtype=wp.float64),
) -> None:
    """Restore ``live <- snap`` iff the outer step was REJECTED (device-only branch)."""
    if accepted[0] == wp.int32(0):
        t = wp.tid()
        live[t] = snap[t]


@wp.kernel
def _commit_bond_if_accepted_kernel(
    accepted: wp.array(dtype=wp.int32),
    dt_phys: wp.float64,
    bound: wp.array(dtype=wp.int32),
    candidate_bound: wp.array(dtype=wp.int32),
    bound_snap: wp.array(dtype=wp.int32),
    age: wp.array(dtype=wp.float64),
) -> None:
    """Promote the proposed bond state and age the bond, iff the outer step was ACCEPTED.

    A slot whose committed state differs from the state it held at snapshot changed in THIS step, so its
    age restarts: a bond that just associated is zero seconds old, and a dissociated slot carries no age.
    """
    if accepted[0] != wp.int32(0):
        t = wp.tid()
        proposed = candidate_bound[t]
        bound[t] = proposed
        if proposed != bound_snap[t]:
            age[t] = wp.float64(0.0)
        elif proposed != wp.int32(0):
            age[t] = age[t] + dt_phys


@dataclass(slots=True)
class CollagenCrosslinkConnector:
    """The ``ecm_crosslink`` bond population — internal to ``ecm``, kinetic, fixed capacity."""

    name: str
    component_a: str
    component_b: str
    device: str
    capacity: int
    bonds_d: wp.array
    bound_d: wp.array
    candidate_bound_d: wp.array
    k_d: wp.array
    r0_d: wp.array
    age_d: wp.array
    bound_snap_d: wp.array
    candidate_bound_snap_d: wp.array
    age_snap_d: wp.array
    adjoint_transfer_required: bool = True
    launch: Callable[..., object] = wp.launch
    copy: Callable[..., object] = wp.copy

    def __post_init__(self) -> None:
        if self.name != ECM_CROSSLINK_CONNECTOR:
            raise ValueError(f"collagen crosslink connector must be {ECM_CROSSLINK_CONNECTOR!r}")
        if {self.component_a, self.component_b} != {ECM_COMPONENT}:
            raise ValueError("ecm_crosslink is INTERNAL — both endpoints must be ecm")
        if not self.adjoint_transfer_required:
            raise ValueError("ECM crosslink connector requires an adjoint force scatter")

    def accumulate_internal(self, pos: wp.array, force: wp.array) -> None:
        """Add the bound crosslink forces to the ECM-owned force array (caller owns zeroing)."""
        if not self.capacity:
            return
        self.launch(
            _bound_link_spring_kernel, dim=self.capacity,
            inputs=[pos, self.bonds_d, self.bound_d, self.k_d, self.r0_d, force],
            device=self.device,
        )

    def snapshot_candidate(self) -> None:
        """Snapshot the committed bond topology before candidate proposals (device D2D)."""
        if not self.capacity:
            return
        self.copy(self.bound_snap_d, self.bound_d)
        self.copy(self.candidate_bound_snap_d, self.candidate_bound_d)
        self.copy(self.age_snap_d, self.age_d)

    def rollback(self, accepted: wp.array) -> None:
        """Discard rejected association/dissociation proposals under the device predicate."""
        if not self.capacity:
            return
        for live, snap in ((self.bound_d, self.bound_snap_d),
                           (self.candidate_bound_d, self.candidate_bound_snap_d)):
            self.launch(
                _restore_int32_if_rejected_kernel, dim=self.capacity,
                inputs=[accepted, live, snap], device=self.device,
            )
        self.launch(
            _restore_float64_if_rejected_kernel, dim=self.capacity,
            inputs=[accepted, self.age_d, self.age_snap_d], device=self.device,
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """Commit association/dissociation only after global acceptance."""
        if not self.capacity:
            return
        self.launch(
            _commit_bond_if_accepted_kernel, dim=self.capacity,
            inputs=[accepted, wp.float64(float(dt_phys)), self.bound_d, self.candidate_bound_d,
                    self.bound_snap_d, self.age_d],
            device=self.device,
        )

    def accumulate_ledger(self, ledger: object) -> None:
        """Contribute the bond topology count; the force resultant is the ECM owner's own array."""
        add_topology = getattr(ledger, "add_topology_count", None)
        if callable(add_topology):
            add_topology(self.name, self.capacity)

    def bound_count(self) -> int:
        """OUT-OF-LOOP diagnostic: currently associated bonds.  Never call inside the physical-time loop."""
        return 0 if not self.capacity else int(np.count_nonzero(self.bound_d.numpy()))


def build_collagen_crosslink_connector(
    *,
    bonds: np.ndarray,
    stiffness_pn_per_um: np.ndarray,
    rest_um: np.ndarray,
    device: str,
    initially_bound: np.ndarray | None = None,
    launch: Callable[..., object] = wp.launch,
    copy: Callable[..., object] = wp.copy,
) -> CollagenCrosslinkConnector:
    """Upload the collagen bond table and bind the state-gated internal spring (CUDA lane).

    Args:
        bonds: ``(J, 2)`` ECM node-index pairs, one row per crosslink slot.
        stiffness_pn_per_um: per-bond stiffness ``(J,)`` [pN/um].  A chemistry-card quantity — required,
            never defaulted here.
        rest_um: per-bond zero-force separation ``(J,)`` [um].
        device: the resolved CUDA device string.
        initially_bound: optional ``(J,)`` 0/1 initial association.  Defaults to all bound, which is the
            configuration the mikado builder produces before any kinetics runs.
        launch: injected kernel launcher (the structural gate passes a recording double).
        copy: injected device-copy (likewise).

    Returns:
        A :class:`CollagenCrosslinkConnector` satisfying
        :class:`~aleph.engine.ecm_world.ECMInternalCrosslink`.

    Raises:
        ValueError: on a malformed table, a non-positive stiffness, a negative rest length, or a
            self-bond (two identical endpoints have no direction).
    """
    pairs = np.asarray(bonds)
    if pairs.ndim != 2 or (pairs.size and pairs.shape[1] != 2):
        raise ValueError("bonds must be a (J, 2) index table")
    capacity = int(pairs.shape[0])
    k = np.asarray(stiffness_pn_per_um, dtype=np.float64).reshape(-1)
    r0 = np.asarray(rest_um, dtype=np.float64).reshape(-1)
    if k.size != capacity or r0.size != capacity:
        raise ValueError("stiffness and rest length need one entry per bond")
    if capacity:
        if int(pairs.min()) < 0:
            raise ValueError("bond node indices must be nonnegative")
        if np.any(pairs[:, 0] == pairs[:, 1]):
            raise ValueError("a crosslink cannot bond a node to itself: it has no direction")
        if not np.all(np.isfinite(k)) or np.any(k <= 0.0):
            raise ValueError("bond stiffness must be finite and positive")
        if not np.all(np.isfinite(r0)) or np.any(r0 < 0.0):
            raise ValueError("bond rest length must be finite and nonnegative")

    bound = (
        np.ones(capacity, dtype=np.int32) if initially_bound is None
        else np.asarray(initially_bound, dtype=np.int32).reshape(-1)
    )
    if bound.size != capacity:
        raise ValueError("initially_bound needs one entry per bond")

    def i32(values: np.ndarray) -> wp.array:
        return wp.array(np.ascontiguousarray(values, dtype=np.int32), dtype=wp.int32, device=device)

    def f64(values: np.ndarray) -> wp.array:
        return wp.array(np.ascontiguousarray(values, dtype=np.float64), dtype=wp.float64, device=device)

    return CollagenCrosslinkConnector(
        name=ECM_CROSSLINK_CONNECTOR,
        component_a=ECM_COMPONENT,
        component_b=ECM_COMPONENT,
        device=device,
        capacity=capacity,
        bonds_d=wp.array(np.ascontiguousarray(pairs, dtype=np.int32), dtype=wp.int32, ndim=2, device=device),
        bound_d=i32(bound),
        candidate_bound_d=i32(bound),
        k_d=f64(k),
        r0_d=f64(r0),
        age_d=f64(np.zeros(capacity)),
        bound_snap_d=i32(bound),
        candidate_bound_snap_d=i32(bound),
        age_snap_d=f64(np.zeros(capacity)),
        launch=launch,
        copy=copy,
    )
