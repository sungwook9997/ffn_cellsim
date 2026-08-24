r"""The ``cortex`` component's PRIVATE device state — the end of the bind-target port.

WHAT CHANGES.  Until now the engine's cortex endpoint was a *view* onto the incumbent's global array:
``scripts/ac_gate_b_cortex_motor_native.py:198,227`` hands ``cell.pos_d``/``cell.f_d`` out under the cortex's
name, and :mod:`aleph.engine.cortex_motor_slice` records in its own Sanity Gate that "cortex is a
bind-target PORT (never a participant)".  Under that arrangement the cortex owns nothing: its "force array"
is the same allocation the myosin, nucleus and membrane accumulate into, so no adjoint pair terminating on
the cortex can be closed — the two halves are the same numbers.  Nine of the 36 connectors terminate here
(:data:`~aleph.engine.cortex_population.CORTEX_CONNECTORS`), which is why the execution plan puts cortex
array ownership on the critical path (``AC_EXECUTION_PLAN_2026-07-25.md`` §10, ``…TRACKS…`` T2).

This module gives the cortex its own ``position_d``/``force_d`` and makes it a real
:class:`~aleph.engine.surface_body.SurfaceComponentStateOwner`: it launches its OWN mechanics on its OWN
arrays, snapshots and reject-restores its OWN position, and reduces its OWN force resultant into one half of
the balance gate — which is the operation that was impossible before, because one array cannot be both halves
of a Newton pair.

NOT A PASSENGER.  The T2 negative control is explicit that distinct arrays with a no-op mechanics is the
failure mode to catch, so this owner refuses to be built without a real force channel: the composite
mechanics binds the SAME two ``ff`` kernels ``ac.cell.driver._accumulate_all`` launches on the cortex
(``cytosim_bending_kernel``, ``link_spring_kernel``) plus the Arp2/3 branch-angle harmonic when the cortex is
mixed, and :class:`CortexCompositeMechanics` raises if every channel is empty.  Deliberately NOT bound here:
axial inextensibility (an NF2007 solver constraint, not an accumulated force — a spring would double-count
it), steric WCA and Biot pressure, which address populations wider than the cortex and belong to the
composed world rather than to this component.  They are named in :data:`CORTEX_CHANNELS_NOT_BOUND` so the
omission is visible in the artifact instead of being discovered later as a missing force.

WHAT THIS DOES NOT DO.  It does not modify the incumbent, and it must not be composed into the same solve as
the incumbent's cortex without an omit-mask — that path double-counts and is PI-gated (``…TRACKS…`` §6 D2).
The supported use today is the same one ``sf_arc`` and ``ecm`` already have: a component that owns private
arrays and can be assembled into a bring-up world with no frozen-driver edit at all.

engine units: length µm, force pN, stiffness pN/µm.  The dataclasses take injected arrays and an injected
launcher so the structural gates drive them on a CUDA-free host; the ``build_*`` functions allocate and are
the gbook A5000 lane.

Sanity Gate (self-tested in ``tests/ac/engine/test_cortex_state.py``):
    * ownership: :func:`assert_component_state_disjoint` refuses any two registered components sharing a
      position or force allocation, and the owner refuses arrays longer than the population it was built
      from — the aliased port's array is ``n_total`` long, so that single length check is what the old
      arrangement cannot pass.
    * boundary/sign: the reject-restore is device-gated on ``accepted[0] == 0`` with no host ``bool()``; an
      accepted step keeps the converged candidate, a rejected one restores bit-exactly.
    * conservation: the ledger contribution is the component's OWN force array reduced into ONE side of the
      balance gate, never a copy computed from the other side.
    * measurement-protocol: the bound channels are the same kernels, in the same argument order and units,
      that the incumbent launches on the cortex, so a native node-by-node parity check is a comparison of
      one model with itself rather than of two models.
    * numerical: no tolerance is introduced here; acceptance tolerances stay with the ledger (D8).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable, Sequence

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.engine.cortex_population import (
    CORTEX_COMPONENT,
    CortexPopulation,
    CortexTopology,
)
from aleph.engine.surface_body import (
    CortexFilamentMechanics,
    SurfaceComponentStateOwner,
    storage_key,
)

__all__ = [
    "CORTEX_CHANNELS_NOT_BOUND",
    "CortexBranchAngleMechanics",
    "CortexCompositeMechanics",
    "CortexStateLedger",
    "CortexStateTransaction",
    "assemble_cortex_motor_port",
    "assemble_cortex_state_owner",
    "assert_component_state_disjoint",
    "build_cortex_motor_port",
    "build_cortex_state_owner",
    "cortex_segment_material_coordinates",
    "incumbent_channels_to_omit",
    "INCUMBENT_CHANNEL_NAME",
]

#: Cortex force channels the incumbent launches that this owner deliberately does NOT bind, with the reason.
#: Recorded rather than silently omitted: a force that is absent for a stated reason is a different thing
#: from a force that was forgotten, and only the artifact can tell them apart.
CORTEX_CHANNELS_NOT_BOUND: dict[str, str] = {
    "axial_inextensibility": (
        "NF2007 constraint carried by the solver, not an accumulated force — binding a spring here would "
        "double-count the backbone"
    ),
    "steric_wca": (
        "the WCA neighbour search addresses every solid node in the cell, not just the cortex; it belongs "
        "to the composed world"
    ),
    "biot_pressure": (
        "-alpha*grad(p)*V is a field->solid coupling owned by the cytosol component and delivered through "
        "surface_porous_transfer, not a cortex-internal channel"
    ),
}

#: This owner's channel name -> the name `ac/cell/driver.py::_accumulate_all` knows the SAME kernel by.
#: The two vocabularies exist because each names the channel from its own side; the map is what keeps a
#: relocation total.  Extend this dict in the SAME change that binds a new channel — a channel bound here
#: and missing here is exactly the silent double-count `omit=` was added to prevent.
INCUMBENT_CHANNEL_NAME: dict[str, str] = {
    "cytosim_bending": "actin_bending_cytosim",
    "crosslink_link_spring": "actin_crosslink_link_spring",
    "arp23_branch_angle": "arp23_branch_angle",
}


def incumbent_channels_to_omit(mechanics: "CortexCompositeMechanics") -> frozenset[str]:
    """The incumbent channels a caller MUST omit once this owner computes the cortex.

    DERIVED, never typed.  The caller cannot know which channels are bound without asking: the branch-angle
    delegate appears only for a mixed formin+Arp2/3 cortex (``n_branch > 0``), so a hand-written omit set is
    correct for the default cortex and silently doubles the Arp2/3 force for the mixed one — the force is
    launched once by this owner and once by the incumbent that was never told to stop.  Deriving it makes the
    omit set a function of what is actually bound, so the two cannot drift apart.

    Args:
        mechanics: The composite mechanics of a built cortex owner (``owner.mechanics``).

    Returns:
        Incumbent channel names, suitable for ``_accumulate_all(omit=)``.

    Raises:
        KeyError: If a bound channel has no incumbent name.  Refusing is the point: an unmapped channel would
            otherwise be dropped from the omit set and double-counted, which is invisible in the result.
    """
    return frozenset(INCUMBENT_CHANNEL_NAME[name] for name in mechanics.channels)


# ── Reject-gated restore: selects the REJECTED branch in-device so no host `bool()` ever sees the predicate.
# The same five-line primitive exists in `cortex_motor_slice` and `erm_cortex_slice`; it is repeated here
# rather than imported so a foundational component module does not depend on a slice module (the dependency
# would run the wrong way and would drag the whole motor stack into every cortex build).
@wp.kernel
def _restore_vec3_if_rejected_kernel(
    accepted: wp.array(dtype=wp.int32),
    live: wp.array(dtype=wp.vec3d),
    snap: wp.array(dtype=wp.vec3d),
) -> None:
    """Restore ``live <- snap`` iff the outer step was rejected (device-only, no host branch)."""
    if accepted[0] == 0:
        t = wp.tid()
        live[t] = snap[t]


def _default_launch(kernel: object, *, dim: int, inputs: list, device: object | None = None) -> None:
    """Production launcher: forward to ``wp.launch`` on the owning CUDA device."""
    wp.launch(kernel, dim=dim, inputs=inputs, device=device)


def _default_copy(dst: object, src: object) -> None:
    """Production D2D snapshot/restore: forward to ``wp.copy``."""
    wp.copy(dst, src)


@dataclass(frozen=True, slots=True)
class CortexBranchAngleMechanics:
    """The Arp2/3 70° branch-angle harmonic over cortex-owned branch triples.

    The mixed formin+Arp2/3 cortex weaves ``(m_after, branch_vertex, daughter_arm)`` triples; this launches
    the SAME ``branch_angle_kernel`` with the SAME Faessler-2020 anchors the incumbent uses
    (``ac.cell.driver._accumulate_all``), so the channel is the incumbent's channel and not a re-derivation.
    A formin-only cortex weaves zero triples, in which case this delegate is not constructed at all.

    Per the charter's worked example the branch is an angle HARMONIC with thermal fluctuation, never a rigid
    72°/70° clamp.
    """

    device: str
    branch_triples_d: wp.array
    branch_active_d: wp.array
    theta0_rad: float
    k_theta: float
    kernel: object
    launch: Callable[..., object] = _default_launch
    n_branch: int = field(init=False)

    def __post_init__(self) -> None:
        shape = getattr(self.branch_triples_d, "shape", None)
        if not isinstance(shape, tuple) or len(shape) != 2 or int(shape[1]) != 3:
            raise ValueError("cortex branch triples must be a (B, 3) int32 device array")
        n_branch = int(shape[0])
        if n_branch <= 0:
            raise ValueError("do not construct the branch delegate for a formin-only cortex (zero triples)")
        active_shape = getattr(self.branch_active_d, "shape", None)
        if not isinstance(active_shape, tuple) or int(active_shape[0]) != n_branch:
            raise ValueError("cortex branch_active length must match the branch triple count")
        if not np.isfinite(self.theta0_rad) or not np.isfinite(self.k_theta) or self.k_theta <= 0.0:
            raise ValueError("branch-angle anchors must be finite with a positive stiffness")
        object.__setattr__(self, "n_branch", n_branch)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Launch the branch-angle harmonic into the cortex-owned force array."""
        self.launch(
            self.kernel,
            dim=self.n_branch,
            inputs=[pos, self.branch_triples_d, self.branch_active_d,
                    self.theta0_rad, self.k_theta, force],
            device=self.device,
        )


@dataclass(frozen=True, slots=True)
class CortexCompositeMechanics:
    """Every cortex-internal force channel, launched on the cortex's own arrays in incumbent order.

    Order matters only for float64 atomic arrival, not for the sum, but it is kept identical to
    ``_accumulate_all`` (bending, crosslink, branch) so a native parity comparison differs by summation
    order alone.  ``n_vertices`` is exposed because
    :class:`~aleph.engine.surface_body.SurfaceComponentStateOwner` asserts the mechanics addresses
    exactly the owned node count — the check that catches a concatenated array.

    The passenger refusal (distinct arrays + a no-op mechanics, the T2 negative control) is NOT repeated
    here: :class:`~aleph.engine.surface_body.CortexFilamentMechanics` already refuses to construct with
    zero links AND zero bending triples, so a composite with no channel at all cannot be reached.  A second
    copy of that check would be unreachable code no test could exercise.
    """

    filament: CortexFilamentMechanics
    branch: CortexBranchAngleMechanics | None = None

    @property
    def n_vertices(self) -> int:
        """Node count this mechanics addresses; must equal the owner's array length."""
        return int(self.filament.n_vertices)

    @property
    def channels(self) -> tuple[str, ...]:
        """Names of the bound force channels, for the run artifact."""
        bound = []
        if self.filament.n_triples:
            bound.append("cytosim_bending")
        if self.filament.n_links:
            bound.append("crosslink_link_spring")
        if self.branch is not None:
            bound.append("arp23_branch_angle")
        return tuple(bound)

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Sum every bound cortex channel into the cortex force array (caller owns zeroing)."""
        self.filament.accumulate(pos, force)
        if self.branch is not None:
            self.branch.accumulate(pos, force)


@dataclass(slots=True)
class CortexStateTransaction:
    """Snapshot / reject-restore of the cortex's OWN position array.

    The cortex carries no kinetic epoch of its own: nucleation, severing and turnover arrive with sourced
    rates in a biology phase, and every kinetic edge terminating on the cortex today
    (``membrane_erm_cortex``, ``nmii_cortex_motor``, …) owns its own epoch on the connector side.  So
    :meth:`commit_irreversible` is a structural no-op and says so, rather than silently advancing a private
    clock the transaction does not know about.
    """

    position_d: wp.array
    position_snap_d: wp.array
    device: object | None = None
    launch: Callable[..., object] = _default_launch
    copy: Callable[..., None] = _default_copy

    def __post_init__(self) -> None:
        if storage_key(self.position_d) == storage_key(self.position_snap_d):
            raise ValueError("the cortex snapshot must be a separate allocation from the live position")
        live_shape = getattr(self.position_d, "shape", None)
        snap_shape = getattr(self.position_snap_d, "shape", None)
        if live_shape != snap_shape:
            raise ValueError("cortex snapshot and live position must have identical shape")

    def snapshot_candidate(self) -> None:
        """D2D-snapshot the cortex position (``snap <- live``) before the candidate solve."""
        self.copy(self.position_snap_d, self.position_d)

    def rollback(self, accepted: wp.array) -> None:
        """Reject-gated restore; an accepted step keeps the converged candidate geometry."""
        self.launch(
            _restore_vec3_if_rejected_kernel,
            dim=int(self.position_d.shape[0]),
            inputs=[accepted, self.position_d, self.position_snap_d],
            device=self.device,
        )

    def commit_irreversible(self, accepted: wp.array, dt_phys: float, rng_seed: int) -> None:
        """No-op: the cortex owns no kinetic epoch yet (turnover rates are an unfilled PI-GAP)."""


@dataclass(slots=True)
class CortexStateLedger:
    """Reduce the cortex's OWN force array into one half of the device balance gate.

    This is the operation array ownership buys.  ``GlobalCellLedger.add_body_force`` requires the two
    channels to be *independently accumulated* halves of one Newton pair; while the cortex shared the
    incumbent's array, its resultant and the other body's resultant were reductions of the same numbers, so
    the gate could not test the adjoint wiring at all.  With a private array the cortex is a genuine side.

    ``side`` is arbitrary but fixed per lane, and it is the caller's job to ensure the far side of each
    connector occupies the other channel.  A ledger without ``add_body_force`` (a recording double, or a
    ledger assembled for a different lane) is skipped rather than adapted — silently routing a body force
    into a differently-meaning channel is precisely the accounting error this gate exists to catch.
    """

    force_d: wp.array
    side: str = "reaction"
    contributions: int = field(default=0, init=False)

    def __post_init__(self) -> None:
        if self.side not in ("reaction", "traction"):
            raise ValueError("cortex ledger side must be 'reaction' or 'traction'")

    def accumulate_ledger(self, ledger: object) -> None:
        """Hand the cortex's own force array to the balance gate; decide no acceptance here."""
        sink = getattr(ledger, "add_body_force", None)
        if not callable(sink):
            return
        sink(self.force_d, side=self.side)
        self.contributions += 1


def _byte_span(array: object) -> tuple[str, int, int] | None:
    """The device byte interval ``[start, end)`` an array addresses, or ``None`` if it is not spannable.

    Pointer identity is not enough once any owner holds a SLICE VIEW.  ``wp.array.__getitem__`` returns a
    real array whose ``ptr`` is the base pointer plus the slice offset, so two views that OVERLAP carry two
    different pointers while addressing the same nodes — the defect a careless split introduces, and one
    that ``storage_key`` equality cannot see.  Spanning the bytes makes overlap the thing being tested
    rather than a proxy for it.

    Only rank-1 arrays are spanned; every component state array here is rank-1, and anything else falls
    back to pointer identity, which is what the check did for everything before.
    """
    ptr = getattr(array, "ptr", None)
    shape = getattr(array, "shape", None)
    strides = getattr(array, "strides", None)
    if ptr is None or not shape or not strides or len(shape) != 1 or len(strides) != 1:
        return None
    n, stride = int(shape[0]), int(strides[0])
    if n <= 0 or stride <= 0:
        return None
    return (str(getattr(array, "device", None)), int(ptr), int(ptr) + n * stride)


def assert_component_state_disjoint(owners: Sequence[object]) -> None:
    """Assert no two registered components address one position or force allocation.

    T2 step (iii): "enforce ``storage_key`` exclusivity across *all* registered components".  The pairwise
    membrane/cortex guard already in :class:`~aleph.engine.surface_body.SurfaceBody` only covers those
    two and — measured 2026-07-25 — had never executed outside tests, because the facade was never
    constructed for real.  This is the same check over an arbitrary component set, so a third component
    joining the composition cannot quietly re-alias.

    Exclusivity is tested on the BYTES a component addresses, not on the pointer it was handed: a slice
    view of another component's buffer has its own pointer (see :func:`_byte_span`).  Passing this says
    the components address disjoint nodes; it does NOT say each holds a private allocation, and no caller
    may report it as one — that distinction is what ``STATE.md`` (c) 4 is about.

    Args:
        owners: objects exposing ``name``, ``position_d`` and ``force_d``.

    Raises:
        ValueError: If two components share an allocation or address overlapping bytes of one, or if a
            component's own position and force alias each other (which would make its force a
            displacement).
    """
    seen: dict[tuple[object, ...], tuple[str, str]] = {}
    spans: list[tuple[str, str, tuple[str, int, int]]] = []
    for owner in owners:
        name = str(getattr(owner, "name", owner.__class__.__name__))
        arrays = (("position_d", getattr(owner, "position_d", None)),
                  ("force_d", getattr(owner, "force_d", None)))
        keys = []
        for label, array in arrays:
            if array is None:
                raise ValueError(f"component {name!r} exposes no {label}; it is not a state owner")
            key = storage_key(array)
            keys.append(key)
            previous = seen.get(key)
            if previous is not None:
                raise ValueError(
                    f"component {name!r}.{label} shares one allocation with {previous[0]!r}."
                    f"{previous[1]} — co-location in an array is NEVER a connection, so this composition "
                    "would report a weld as a coupling"
                )
            seen[key] = (name, label)
            span = _byte_span(array)
            if span is None:
                continue
            for other_name, other_label, other in spans:
                if other[0] == span[0] and other[1] < span[2] and span[1] < other[2]:
                    raise ValueError(
                        f"component {name!r}.{label} overlaps {other_name!r}.{other_label} in device "
                        f"memory — two views of one buffer are one allocation wearing two pointers, and "
                        "co-location in an array is NEVER a connection"
                    )
            spans.append((name, label, span))
        if keys[0] == keys[1]:
            raise ValueError(f"component {name!r} position and force must not alias")


def cortex_segment_material_coordinates(
    topology: CortexTopology,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Per-segment arc-length material coordinates ``(s0, s1)`` along each owning filament [µm].

    A motor head binds a *material coordinate on a polar filament*, so the port's ``material_s0``/``s1``
    must be an arc length, not a segment ordinal.  The coordinate runs from 0 at each filament's first node
    and accumulates that filament's NF2007 segment rest lengths, so it is a REST-configuration material
    coordinate: it is fixed by construction and does not drift as the cortex deforms, which is what makes it
    usable as a persistent identity across a topology epoch.

    Returns:
        ``(s0, s1)``, each ``(S,)`` float64, with ``s1 > s0`` for every segment.
    """
    off = np.asarray(topology.fiber_offsets, np.int64)
    rest = np.asarray(topology.seg_rest, np.float64)
    seg_per_fiber = np.diff(off) - 1
    seg_start = np.concatenate([[0], np.cumsum(seg_per_fiber)]).astype(np.int64)
    s0 = np.zeros(int(rest.shape[0]), np.float64)
    s1 = np.zeros(int(rest.shape[0]), np.float64)
    for f in range(int(seg_per_fiber.shape[0])):
        a, b = int(seg_start[f]), int(seg_start[f + 1])
        if b <= a:
            continue
        cumulative = np.concatenate([[0.0], np.cumsum(rest[a:b])])
        s0[a:b] = cumulative[:-1]
        s1[a:b] = cumulative[1:]
    return s0, s1


def assemble_cortex_state_owner(
    *,
    position_d: wp.array,
    force_d: wp.array,
    position_snap_d: wp.array,
    mechanics: CortexCompositeMechanics,
    device: object | None = None,
    launch: Callable[..., object] = _default_launch,
    copy: Callable[..., None] = _default_copy,
    ledger_side: str = "reaction",
) -> SurfaceComponentStateOwner:
    """Wire already-allocated cortex arrays into a real state owner (no allocation, host-testable).

    Separated from :func:`build_cortex_state_owner` so the structural gates can drive the whole ownership
    contract with recording doubles on a CUDA-free host, exactly as the SF lane does.

    Raises:
        ValueError: via :class:`~aleph.engine.surface_body.SurfaceComponentStateOwner` if the arrays
            alias, differ in shape or device, or if the mechanics addresses a different node count — the
            last being the check the aliased global array cannot pass.
    """
    return SurfaceComponentStateOwner(
        name=CORTEX_COMPONENT,
        position_d=position_d,
        force_d=force_d,
        mechanics=mechanics,
        transaction=CortexStateTransaction(
            position_d=position_d, position_snap_d=position_snap_d,
            device=device, launch=launch, copy=copy,
        ),
        ledger=CortexStateLedger(force_d=force_d, side=ledger_side),
    )


def assemble_cortex_motor_port(
    *,
    component: str,
    position_d: wp.array,
    force_d: wp.array,
    segment_node_a_d: wp.array,
    segment_node_b_d: wp.array,
    segment_polarity_d: wp.array,
    persistent_filament_id_d: wp.array,
    material_s0_d: wp.array,
    material_s1_d: wp.array,
    topology_epoch_d: wp.array,
) -> object:
    """Wrap cortex-owned arrays as the ``nmii_cortex_motor`` filament port (no allocation).

    Imported lazily so a cortex build does not pull the NMII actuator stack.  The point of routing the port
    through the component's own arrays is that the head's ``+f`` and the two actin nodes' ``-(1-t)f``/``-t f``
    then land in two DIFFERENT allocations, which is the only arrangement in which the split is physical
    rather than logical.
    """
    from aleph.engine.nmii_actuator import FilamentMotorPortView

    return FilamentMotorPortView(
        component=component,
        position_d=position_d,
        force_d=force_d,
        segment_node_a_d=segment_node_a_d,
        segment_node_b_d=segment_node_b_d,
        segment_polarity_d=segment_polarity_d,
        persistent_filament_id_d=persistent_filament_id_d,
        material_s0_d=material_s0_d,
        material_s1_d=material_s1_d,
        topology_epoch_d=topology_epoch_d,
    )


# ── CUDA-lane builders (allocate device memory; the gbook A5000, never the dev Mac) ──────────────────────
def _upload(array: npt.NDArray, dtype: object, device: object | None) -> wp.array:
    return wp.array(np.ascontiguousarray(array), dtype=dtype, device=device)


def build_cortex_state_owner(
    topology: CortexTopology,
    *,
    device: str,
    launch: Callable[..., object] = _default_launch,
    copy: Callable[..., None] = _default_copy,
    ledger_side: str = "reaction",
) -> SurfaceComponentStateOwner:
    """Allocate the cortex's PRIVATE arrays and bind its real mechanics (CUDA lane).

    Allocates exactly ``topology.n_nodes`` positions, forces and one snapshot — never ``n_total`` — so the
    owner's arrays cannot address another component's node even by accident.  The kernels bound are the
    incumbent's own (``ff.forces_warp.cytosim_bending_kernel``, ``ff.network_warp.link_spring_kernel``, and
    ``ac.weave.branch_angle_warp.branch_angle_kernel`` for a mixed cortex).

    Args:
        topology: the cortex-local tables from
            :func:`~aleph.engine.cortex_population.build_cortex_population`.
        device: the resolved CUDA device string (never a hard-coded ordinal).
        launch / copy: injected CUDA seam; defaults forward to Warp.
        ledger_side: which half of the balance gate this component occupies.

    Returns:
        A :class:`~aleph.engine.surface_body.SurfaceComponentStateOwner` named ``"cortex"``.
    """
    topology.assert_local()
    n = topology.n_nodes
    position_d = _upload(topology.pos, wp.vec3d, device)
    force_d = wp.zeros(n, dtype=wp.vec3d, device=device)
    position_snap_d = wp.zeros(n, dtype=wp.vec3d, device=device)

    filament = CortexFilamentMechanics.bind_native(
        device=str(device),
        n_vertices=n,
        links_d=_upload(topology.links.reshape(-1, 2), wp.int32, device),
        link_k_d=_upload(topology.link_k, wp.float64, device),
        link_r0_d=_upload(topology.link_r0, wp.float64, device),
        bend_triples_d=_upload(topology.bend_triples.reshape(-1, 3), wp.int32, device),
        bend_alpha_d=_upload(topology.bend_alpha, wp.float64, device),
        launch=launch,
    )
    branch = None
    if topology.n_branch:
        from aleph.components.weave.branch_angle import ARP23_K_THETA, ARP23_THETA0_RAD
        from aleph.components.weave.branch_angle_warp import branch_angle_kernel

        branch = CortexBranchAngleMechanics(
            device=str(device),
            branch_triples_d=_upload(topology.branch_triples.reshape(-1, 3), wp.int32, device),
            branch_active_d=_upload(topology.branch_active, wp.int32, device),
            theta0_rad=float(ARP23_THETA0_RAD),
            k_theta=float(ARP23_K_THETA),
            kernel=branch_angle_kernel,
            launch=launch,
        )
    return assemble_cortex_state_owner(
        position_d=position_d,
        force_d=force_d,
        position_snap_d=position_snap_d,
        mechanics=CortexCompositeMechanics(filament=filament, branch=branch),
        device=device,
        launch=launch,
        copy=copy,
        ledger_side=ledger_side,
    )


def build_cortex_motor_port(
    population: CortexPopulation,
    owner: SurfaceComponentStateOwner,
    *,
    device: str,
) -> object:
    """Upload the cortex segment table as the ``nmii_cortex_motor`` port over the OWNER's arrays (CUDA lane).

    The persistent filament id is the population ledger's GLOBAL unique id
    (``id_base + node_fiber[segment's first node]``), so a head's binding record names a filament that
    belongs to exactly one component — the device-side form of the no-double-count invariant.  The material
    coordinates are rest arc lengths from :func:`cortex_segment_material_coordinates`, not segment ordinals.
    """
    topology = population.topology
    id_base = int(population.ledger.block[0])
    seg_fid = (id_base + np.asarray(topology.node_fiber, np.int64)[
        np.asarray(topology.seg_node_a, np.int64)]).astype(np.int64)
    s0, s1 = cortex_segment_material_coordinates(topology)
    return assemble_cortex_motor_port(
        component=CORTEX_COMPONENT,
        position_d=owner.position_d,
        force_d=owner.force_d,
        segment_node_a_d=_upload(topology.seg_node_a, wp.int32, device),
        segment_node_b_d=_upload(topology.seg_node_b, wp.int32, device),
        segment_polarity_d=_upload(topology.seg_polarity, wp.int32, device),
        persistent_filament_id_d=_upload(seg_fid, wp.int64, device),
        material_s0_d=_upload(s0, wp.float64, device),
        material_s1_d=_upload(s1, wp.float64, device),
        topology_epoch_d=wp.zeros(1, dtype=wp.int32, device=device),
    )


def cortex_ownership_census(
    population: CortexPopulation,
    owner: SurfaceComponentStateOwner,
    *,
    other_owners: Iterable[object] = (),
) -> dict[str, object]:
    """Assemble the ownership facts a run artifact must carry, and check them while assembling.

    Every field is MEASURED from the built objects, never declared: the array lengths and the distinct
    allocation identities are read off the arrays themselves, so an artifact claiming ownership that the
    build does not have cannot be produced.
    """
    owners = [owner, *other_owners]
    assert_component_state_disjoint(owners)
    mechanics = owner.mechanics
    return {
        "component": owner.name,
        "population": population.census,
        "owned_position_len": int(owner.position_d.shape[0]),
        "owned_force_len": int(owner.force_d.shape[0]),
        "arrays_private": True,
        "position_storage": str(storage_key(owner.position_d)),
        "force_storage": str(storage_key(owner.force_d)),
        "bound_force_channels": list(getattr(mechanics, "channels", ())),
        "channels_not_bound": dict(CORTEX_CHANNELS_NOT_BOUND),
        "connector_domains": {
            name: {"domain": domain.domain, "size": int(domain.size)}
            for name, domain in sorted(population.endpoints.items())
        },
        "provenance": dict(population.provenance),
    }
