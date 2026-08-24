r"""Segment-anchored, accepted-predicated NMII KMC + crossbridge (P0#2 production motor kernels).

The NG-1 motor anchors a bound head to a single actin NODE; production actin is a segmented filament, so a head
binds a POINT ON A SEGMENT (from :mod:`aleph.components.motor.segment_query`). This module is the device-resident
motor Codex's canonical outer loop wires per Codex's fixed order — KMC runs ONCE per ACCEPTED physical step:

    mechanics converge → outer accept decision → head load → point-to-segment query
      → attach (accepted-predicated) → walk-dir hand-off → Hill step / Bell detach (accepted-predicated)

Every state mutation is predicated on the device ``accepted`` scalar (the scheduler's final outer-acceptance
predicate, after fluid/query/remap validation):
a REJECTED outer step advances NO motor state (bound / segment anchor / abscissa / walk_dir), so the scheduler's
rollback boundary is honoured with no host branch (I0-A: no D2H in the hot loop).

Segment anchor state (per head): ``bound``, ``seg_a`` / ``seg_b`` (the bound segment's two node indices),
``bary_t`` (the barycentric fraction of the attachment on that segment), ``abscissa`` (the walked displacement),
``walk_dir`` (the barbed-end direction from the query — the I4 polarity hand-off). The crossbridge anchors at

    x_att = (1 − t)·pos[seg_a] + t·pos[seg_b] + abscissa·walk_dir

and its reaction is split between the two segment nodes by the barycentric weights ``(1 − t)`` and ``t`` so
Newton's 3rd law holds across the segment (``f_head = −(f_seg_a + f_seg_b)`` exactly). The Hill/Bell laws and
the F4 tangential/full-|F| load split are identical to :mod:`aleph.components.motor.hand` — only the anchor changes.

The NumPy :func:`crossbridge_segment_reference` is the host oracle (dev-Mac I0-A: kernels not launched here).
"""

from __future__ import annotations

from enum import StrEnum

import numpy as np
import numpy.typing as npt
import warp as wp

from aleph.components.motor.hand import (
    NMIICatchSlipParams,
    NMIIHandParams,
    attach_prob,
    bell_off_rate,
    catch_slip_off_rate,
    detach_prob,
    hill_velocity,
)
from aleph.components.motor.segment_query import SegmentQuery

__all__ = [
    "allocate_segment_hand_state",
    "attach_segment_gated_kernel",
    "attach_segment_gated_r0bind_kernel",
    "step_detach_segment_gated_kernel",
    "step_detach_segment_gated_catch_slip_kernel",
    "SegmentDetachKinetics",
    "crossbridge_segment_kernel",
    "compute_head_loads_segment_kernel",
    "crossbridge_segment_split_kernel",
    "compute_head_loads_segment_split_kernel",
    "crossbridge_segment_split_r0bind_kernel",
    "compute_head_loads_segment_split_r0bind_kernel",
    "refresh_segment_barbed_kernel",
    "SegmentMotorRuntime",
    "crossbridge_segment_reference",
    "crossbridge_segment_split_reference",
    "crossbridge_segment_split_r0bind_reference",
    "r0_bind_at_attach_reference",
]


class SegmentDetachKinetics(StrEnum):
    """Selectable per-head detachment law for the segment-anchored NMII KMC (both mechanistic Bell-Evans).

    ``SLIP`` is the pure Bell slip (:func:`aleph.components.motor.hand.bell_off_rate`) — off-rate rises
    monotonically with load; correct for the passive slip linkers.  ``CATCH_SLIP`` is the Pereverzev
    two-pathway catch-slip (:func:`aleph.components.motor.hand.catch_slip_off_rate`) — the physiological
    NMII head-actin bond whose off-rate FALLS with load up to the peak-lifetime force ``F*`` (Kovacs
    2007), so the engaged fraction RISES under resistive load for sustained cortical tension.  This is
    a mode SELECTION over two already-landed device laws — it introduces no new physics and never
    removes the slip path.
    """

    SLIP = "slip"
    CATCH_SLIP = "catch_slip"


def allocate_segment_hand_state(n_hands: int, device: str | None = None) -> dict:
    """All-free segment-anchored hand population (bound=0, seg_a/seg_b=-1, bary_t=0, abscissa=0, walk_dir=0)."""
    return {
        "bound": wp.zeros(n_hands, dtype=wp.int32, device=device),
        "seg_id": wp.full(n_hands, -1, dtype=wp.int32, device=device),
        "seg_a": wp.full(n_hands, -1, dtype=wp.int32, device=device),
        "seg_b": wp.full(n_hands, -1, dtype=wp.int32, device=device),
        "bary_t": wp.zeros(n_hands, dtype=wp.float64, device=device),
        "abscissa": wp.zeros(n_hands, dtype=wp.float64, device=device),
        "walk_dir": wp.zeros(n_hands, dtype=wp.vec3d, device=device),
    }


@wp.kernel
def attach_segment_gated_kernel(
    accepted: wp.array(dtype=wp.int32),
    rng_epoch: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    bound_seg_id: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    q_seg_id: wp.array(dtype=wp.int32),
    q_t: wp.array(dtype=wp.float64),
    q_barbed: wp.array(dtype=wp.vec3d),
    seg_node_a: wp.array(dtype=wp.int32),
    seg_node_b: wp.array(dtype=wp.int32),
    params: NMIIHandParams,
    tau: wp.float64,
    rng_seed: wp.int32,
):
    """Accepted-predicated attach: a free head binds the queried segment (walk-dir = its barbed polarity).

    No-op if the outer step was rejected (``accepted[0]==0``), if the head is already bound, or if the query
    found no segment within capture (``q_seg_id < 0``). Otherwise binds with the Poisson attach probability.
    """
    if accepted[0] == 0:
        return
    h = wp.tid()
    if bound[h] == 1:
        return
    sid = q_seg_id[h]
    if sid < 0:
        return
    state = wp.rand_init(rng_seed ^ rng_epoch[0], h)
    if wp.float64(wp.randf(state)) < attach_prob(tau, params.k_on):
        bound[h] = 1
        bound_seg_id[h] = sid
        seg_a[h] = seg_node_a[sid]
        seg_b[h] = seg_node_b[sid]
        bary_t[h] = q_t[h]
        walk_dir[h] = q_barbed[h]
        abscissa[h] = wp.float64(0.0)


@wp.kernel
def step_detach_segment_gated_kernel(
    accepted: wp.array(dtype=wp.int32),
    rng_epoch: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    bound_seg_id: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    loads_hill: wp.array(dtype=wp.float64),
    loads_bell: wp.array(dtype=wp.float64),
    params: NMIIHandParams,
    tau: wp.float64,
    rng_seed: wp.int32,
):
    """Accepted-predicated step+detach: advance abscissa by τ·v_Hill(tangential), detach by Bell(full |F|).

    Identical kinetics to :func:`aleph.components.motor.hand.step_detach_kernel` (F4 split) but no-op on a rejected
    outer step, and resets the SEGMENT anchor (seg_a/seg_b) on detach.
    """
    if accepted[0] == 0:
        return
    h = wp.tid()
    if bound[h] == 0:
        return
    # DERIVED stall-crossing guard: cap Δa by the Hill-stall bound (f_stall − F)/k_xb so the along-actin
    # tension cannot overshoot f_stall within a tick (the same overshoot fix as hand.step_detach_kernel).
    da = tau * hill_velocity(loads_hill[h], params.v0, params.f_stall, params.kappa)
    da_max = wp.max(wp.float64(0.0), (params.f_stall - loads_hill[h]) / params.k_xb)
    abscissa[h] = abscissa[h] + wp.min(da, da_max)
    p_off = bell_off_rate(loads_bell[h], params.k_off0, params.f0)
    state = wp.rand_init(rng_seed ^ rng_epoch[0], h)
    if wp.float64(wp.randf(state)) < detach_prob(tau, p_off):
        bound[h] = 0
        bound_seg_id[h] = wp.int32(-1)
        seg_a[h] = wp.int32(-1)
        seg_b[h] = wp.int32(-1)
        bary_t[h] = wp.float64(0.0)
        abscissa[h] = wp.float64(0.0)
        walk_dir[h] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))


@wp.kernel
def step_detach_segment_gated_catch_slip_kernel(
    accepted: wp.array(dtype=wp.int32),
    rng_epoch: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    bound_seg_id: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    loads_hill: wp.array(dtype=wp.float64),
    loads_bell: wp.array(dtype=wp.float64),
    params: NMIIHandParams,
    cs: NMIICatchSlipParams,
    tau: wp.float64,
    rng_seed: wp.int32,
):
    r"""Accepted-predicated step+detach with the PHYSIOLOGICAL NMII catch-slip off-rate (Kovacs 2007).

    Identical Hill stepping, DERIVED stall-crossing guard, and segment-anchor reset as
    :func:`step_detach_segment_gated_kernel`, but the ``loads_bell`` full-|F| detachment hazard is the
    Pereverzev two-pathway catch-slip law (:func:`aleph.components.motor.hand.catch_slip_off_rate`) instead of the
    pure Bell slip.  This is the Kovacs 2007 fix: resistive load SLOWS detachment (catch) up to the
    peak-lifetime force ``F*``, so the engaged fraction RISES under load for sustained cortical tension, then
    falls past ``F*`` (slip).  ``params`` supplies only the Hill mechanics (``v0``/``f_stall``/``kappa``/
    ``k_xb``); its ``k_off0``/``f0`` are unused here — the off-rate comes entirely from ``cs``.  The engaged
    fraction still EMERGES from ``k_on``/``k_off`` (P2), never an imposed duty.  No-op on a rejected outer step
    (``accepted[0]==0``), so the scheduler's rollback boundary is honoured with no host branch.
    """
    if accepted[0] == 0:
        return
    h = wp.tid()
    if bound[h] == 0:
        return
    # DERIVED stall-crossing guard: cap Δa by the Hill-stall bound (f_stall − F)/k_xb (same as the slip kernel).
    da = tau * hill_velocity(loads_hill[h], params.v0, params.f_stall, params.kappa)
    da_max = wp.max(wp.float64(0.0), (params.f_stall - loads_hill[h]) / params.k_xb)
    abscissa[h] = abscissa[h] + wp.min(da, da_max)
    p_off = catch_slip_off_rate(loads_bell[h], cs.k_catch0, cs.x_catch, cs.k_slip0, cs.x_slip, cs.kT)
    state = wp.rand_init(rng_seed ^ rng_epoch[0], h)
    if wp.float64(wp.randf(state)) < detach_prob(tau, p_off):
        bound[h] = 0
        bound_seg_id[h] = wp.int32(-1)
        seg_a[h] = wp.int32(-1)
        seg_b[h] = wp.int32(-1)
        bary_t[h] = wp.float64(0.0)
        abscissa[h] = wp.float64(0.0)
        walk_dir[h] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))


@wp.func
def _segment_attachment(pos_a: wp.vec3d, pos_b: wp.vec3d, t: wp.float64, abscissa: wp.float64,
                        walk_dir: wp.vec3d) -> wp.vec3d:
    """Crossbridge anchor on a segment, advanced by the power stroke: ``(1−t)·a + t·b + abscissa·walk_dir``."""
    return (wp.float64(1.0) - t) * pos_a + t * pos_b + abscissa * walk_dir


@wp.kernel
def crossbridge_segment_kernel(
    pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    k_xb: wp.float64,
    r0_xb: wp.float64,
    force: wp.array(dtype=wp.vec3d),
):
    """Segment-anchored power-stroke crossbridge; the reaction splits ``(1−t):t`` between the segment nodes.

    ``F = k_xb (|r| − r0_xb) r̂`` from the head to the advanced segment attachment. Newton's 3rd law holds
    across the segment: the head gets ``+f`` and the two segment nodes get ``−(1−t)·f`` and ``−t·f`` (their sum
    is ``−f``), so a barycentric attachment loads the two real actin nodes exactly.
    """
    h = wp.tid()
    if bound[h] == 0:
        return
    a = seg_a[h]
    b = seg_b[h]
    if a < 0 or b < 0:
        return
    t = bary_t[h]
    hn = head_node[h]
    x_att = _segment_attachment(pos[a], pos[b], t, abscissa[h], walk_dir[h])
    d = x_att - pos[hn]
    # DIRECTIONAL (tangential) crossbridge (2026-07-23 defect#3 fix): the myosin power stroke acts ALONG the
    # actin (``walk_dir``, unit); the perpendicular head<->actin offset is carried by the compliant lever ARM,
    # not the crossbridge. Matches the Hill/Bell load model (``dot(d, walk_dir)``). The prior isotropic 3-D
    # spring made an off-actin head carry a spurious NORMAL force ``k_xb*offset`` (native 645 pN at 0.6 µm).
    f = (k_xb * (wp.dot(d, walk_dir[h]) - r0_xb)) * walk_dir[h]
    wp.atomic_add(force, hn, f)
    wp.atomic_add(force, a, -(wp.float64(1.0) - t) * f)
    wp.atomic_add(force, b, -t * f)


@wp.kernel
def compute_head_loads_segment_kernel(
    pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    k_xb: wp.float64,
    r0_xb: wp.float64,
    loads: wp.array(dtype=wp.float64),
    loads_full: wp.array(dtype=wp.float64),
):
    """F4/F5 split for segment anchors: tangential (walk-dir) load → Hill; full |F| load → Bell."""
    h = wp.tid()
    if bound[h] == 0 or seg_a[h] < 0:
        loads[h] = wp.float64(0.0)
        loads_full[h] = wp.float64(0.0)
        return
    x_att = _segment_attachment(pos[seg_a[h]], pos[seg_b[h]], bary_t[h], abscissa[h], walk_dir[h])
    d = x_att - pos[head_node[h]]
    loads[h] = k_xb * (wp.dot(d, walk_dir[h]) - r0_xb)
    # Directional crossbridge: the applied force is tangential, so |F| == |tangential load| (Bell sees the
    # tangential magnitude, not the 3-D deformation — no spurious detachment from a geometric offset).
    loads_full[h] = wp.abs(loads[h])


@wp.kernel
def crossbridge_segment_split_kernel(
    actuator_pos: wp.array(dtype=wp.vec3d),
    actuator_force: wp.array(dtype=wp.vec3d),
    port_pos: wp.array(dtype=wp.vec3d),
    port_force: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    k_xb: wp.float64,
    r0_xb: wp.float64,
):
    """Two-array adjoint crossbridge for the split-ownership MOTOR connector (physics-preserving relayout).

    Identical power-stroke formula to :func:`crossbridge_segment_kernel`, but the head node is addressed in the
    ``nmii``-owned ``actuator_pos``/``actuator_force`` arrays while the two actin segment nodes are addressed in
    the *target*-owned ``port_pos``/``port_force`` arrays.  The head node index (``head_node``) indexes the
    actuator array; ``seg_a``/``seg_b`` index the port array.  The arrays are NEVER merged: there is no shared
    ``pos``/``force`` and no co-location assumption.  Newton's 3rd law still holds *across* the two arrays —
    the head gets ``+f`` (into ``actuator_force``) and the segment endpoints get ``−(1−t)·f`` and ``−t·f``
    (into ``port_force``), whose resultant is exactly ``−f``.
    """
    h = wp.tid()
    if bound[h] == 0:
        return
    a = seg_a[h]
    b = seg_b[h]
    if a < 0 or b < 0:
        return
    t = bary_t[h]
    hn = head_node[h]
    x_att = _segment_attachment(port_pos[a], port_pos[b], t, abscissa[h], walk_dir[h])
    d = x_att - actuator_pos[hn]
    # DIRECTIONAL (tangential) crossbridge — same 2026-07-23 fix as crossbridge_segment_kernel, across the two
    # split-ownership arrays (Newton's 3rd law still holds: head +f into actuator, segment reactions into port).
    f = (k_xb * (wp.dot(d, walk_dir[h]) - r0_xb)) * walk_dir[h]
    wp.atomic_add(actuator_force, hn, f)
    wp.atomic_add(port_force, a, -(wp.float64(1.0) - t) * f)
    wp.atomic_add(port_force, b, -t * f)


@wp.kernel
def compute_head_loads_segment_split_kernel(
    actuator_pos: wp.array(dtype=wp.vec3d),
    port_pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    k_xb: wp.float64,
    r0_xb: wp.float64,
    loads: wp.array(dtype=wp.float64),
    loads_full: wp.array(dtype=wp.float64),
):
    """Split-ownership F4/F5 per-head load: head from ``actuator_pos``, segment from ``port_pos``.

    Same tangential(walk-dir)→Hill / full-|F|→Bell split as :func:`compute_head_loads_segment_kernel`, but the
    head node reads the ``nmii``-owned array and the segment nodes read the target-owned array; the two arrays
    are not merged.
    """
    h = wp.tid()
    if bound[h] == 0 or seg_a[h] < 0:
        loads[h] = wp.float64(0.0)
        loads_full[h] = wp.float64(0.0)
        return
    x_att = _segment_attachment(port_pos[seg_a[h]], port_pos[seg_b[h]], bary_t[h], abscissa[h], walk_dir[h])
    d = x_att - actuator_pos[head_node[h]]
    loads[h] = k_xb * (wp.dot(d, walk_dir[h]) - r0_xb)
    # Directional crossbridge: |F| == |tangential load| (same fix as compute_head_loads_segment_kernel).
    loads_full[h] = wp.abs(loads[h])


@wp.kernel
def attach_segment_gated_r0bind_kernel(
    accepted: wp.array(dtype=wp.int32),
    rng_epoch: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    bound_seg_id: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    r0_bind: wp.array(dtype=wp.float64),
    q_seg_id: wp.array(dtype=wp.int32),
    q_t: wp.array(dtype=wp.float64),
    q_barbed: wp.array(dtype=wp.vec3d),
    seg_node_a: wp.array(dtype=wp.int32),
    seg_node_b: wp.array(dtype=wp.int32),
    actuator_pos: wp.array(dtype=wp.vec3d),
    port_pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    params: NMIIHandParams,
    tau: wp.float64,
    rng_seed: wp.int32,
):
    r"""Accepted-predicated attach that ALSO captures the per-head zero-strain reference ``r0_bind``.

    Identical binding rule to :func:`attach_segment_gated_kernel`, but on a successful bind it records

        r0_bind[h] = dot(x_att − x_head, walk_dir)          (with abscissa=0, so x_att = (1−t)·a + t·b)

    the initial along-walk projection of the head→attachment offset.  The crossbridge/load kernels then use
    ``r0_bind[h]`` as the crossbridge rest, so the PASSIVE force is exactly 0 at the instant of binding and the
    ACTIVE force develops ONLY from the subsequent power stroke (``abscissa``) and relative sliding.  This is the
    dynamic-path realisation of the ``⟨offset, ŵ⟩`` correction ``resting_setpoint`` already applies to the seed
    path (RESTING_BASELINE_DIAGNOSIS): on the coarse cortex mesh a head binds the NEAREST segment up to the
    capture radius away, so a scalar ``r0_xb=0`` (head-on-actin) would inject a spurious ``k_xb·⟨offset, ŵ⟩``
    placement force; ``r0_bind`` removes that discretization artifact resolution-independently.

    The head node reads the ``nmii``-owned ``actuator_pos``; the segment nodes read the target-owned
    ``port_pos`` — the two arrays are addressed independently (split ownership), never merged.
    """
    if accepted[0] == 0:
        return
    h = wp.tid()
    if bound[h] == 1:
        return
    sid = q_seg_id[h]
    if sid < 0:
        return
    state = wp.rand_init(rng_seed ^ rng_epoch[0], h)
    if wp.float64(wp.randf(state)) < attach_prob(tau, params.k_on):
        a = seg_node_a[sid]
        b = seg_node_b[sid]
        t = q_t[h]
        bound[h] = 1
        bound_seg_id[h] = sid
        seg_a[h] = a
        seg_b[h] = b
        bary_t[h] = t
        walk_dir[h] = q_barbed[h]
        abscissa[h] = wp.float64(0.0)
        # zero-strain reference at bind (abscissa=0): the along-walk projection of the head→attachment offset.
        x_att = (wp.float64(1.0) - t) * port_pos[a] + t * port_pos[b]
        r0_bind[h] = wp.dot(x_att - actuator_pos[head_node[h]], q_barbed[h])


@wp.kernel
def crossbridge_segment_split_r0bind_kernel(
    actuator_pos: wp.array(dtype=wp.vec3d),
    actuator_force: wp.array(dtype=wp.vec3d),
    port_pos: wp.array(dtype=wp.vec3d),
    port_force: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    k_xb: wp.float64,
    r0_bind: wp.array(dtype=wp.float64),
):
    r"""Two-array adjoint crossbridge with a PER-HEAD zero-strain reference (attach-unstrained).

    Identical power-stroke formula and Newton-3rd split as :func:`crossbridge_segment_split_kernel`, except the
    scalar rest ``r0_xb`` is replaced by the per-head ``r0_bind[h]`` captured at attach:

        f = k_xb · (dot(x_att − x_head, walk_dir) − r0_bind[h]) · walk_dir

    so the passive force is 0 at the binding instant and builds only from the walked ``abscissa`` and relative
    sliding — the coarse-mesh placement force ``k_xb·⟨offset, ŵ⟩`` is removed.  ``r0_bind`` is a physical zero-
    strain length (grounded in the bind geometry), never a tuned magnitude.
    """
    h = wp.tid()
    if bound[h] == 0:
        return
    a = seg_a[h]
    b = seg_b[h]
    if a < 0 or b < 0:
        return
    t = bary_t[h]
    hn = head_node[h]
    x_att = _segment_attachment(port_pos[a], port_pos[b], t, abscissa[h], walk_dir[h])
    d = x_att - actuator_pos[hn]
    f = (k_xb * (wp.dot(d, walk_dir[h]) - r0_bind[h])) * walk_dir[h]
    wp.atomic_add(actuator_force, hn, f)
    wp.atomic_add(port_force, a, -(wp.float64(1.0) - t) * f)
    wp.atomic_add(port_force, b, -t * f)


@wp.kernel
def compute_head_loads_segment_split_r0bind_kernel(
    actuator_pos: wp.array(dtype=wp.vec3d),
    port_pos: wp.array(dtype=wp.vec3d),
    head_node: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    seg_a: wp.array(dtype=wp.int32),
    seg_b: wp.array(dtype=wp.int32),
    bary_t: wp.array(dtype=wp.float64),
    abscissa: wp.array(dtype=wp.float64),
    walk_dir: wp.array(dtype=wp.vec3d),
    k_xb: wp.float64,
    r0_bind: wp.array(dtype=wp.float64),
    loads: wp.array(dtype=wp.float64),
    loads_full: wp.array(dtype=wp.float64),
):
    """Split F4/F5 per-head load with the per-head zero-strain reference ``r0_bind`` (attach-unstrained).

    Same tangential(walk-dir)→Hill / full-|F|→Bell split as :func:`compute_head_loads_segment_split_kernel`,
    but the crossbridge rest is ``r0_bind[h]`` (0 tangential load at the binding instant), so the Hill stepping
    and Bell detachment see the true POWER-STROKE load, not the coarse-mesh placement artifact.
    """
    h = wp.tid()
    if bound[h] == 0 or seg_a[h] < 0:
        loads[h] = wp.float64(0.0)
        loads_full[h] = wp.float64(0.0)
        return
    x_att = _segment_attachment(port_pos[seg_a[h]], port_pos[seg_b[h]], bary_t[h], abscissa[h], walk_dir[h])
    d = x_att - actuator_pos[head_node[h]]
    loads[h] = k_xb * (wp.dot(d, walk_dir[h]) - r0_bind[h])
    loads_full[h] = wp.abs(loads[h])


@wp.kernel
def refresh_segment_barbed_kernel(
    pos: wp.array(dtype=wp.vec3d),
    seg_node_a: wp.array(dtype=wp.int32),
    seg_node_b: wp.array(dtype=wp.int32),
    seg_polarity: wp.array(dtype=wp.int32),
    seg_barbed: wp.array(dtype=wp.vec3d),
) -> None:
    """Refresh each segment's live unit direction toward its filament's barbed end."""
    s = wp.tid()
    direction = pos[seg_node_b[s]] - pos[seg_node_a[s]]
    if seg_polarity[s] < 0:
        direction = -direction
    length = wp.length(direction)
    if length <= wp.float64(1.0e-12):
        seg_barbed[s] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
    else:
        seg_barbed[s] = direction / length


@wp.kernel
def _refresh_bound_walk_dir_kernel(
    pos: wp.array(dtype=wp.vec3d),
    bound: wp.array(dtype=wp.int32),
    bound_seg_id: wp.array(dtype=wp.int32),
    seg_node_a: wp.array(dtype=wp.int32),
    seg_node_b: wp.array(dtype=wp.int32),
    seg_polarity: wp.array(dtype=wp.int32),
    walk_dir: wp.array(dtype=wp.vec3d),
) -> None:
    """Make a bound head follow its live segment without refreshing every unbound-query segment."""
    h = wp.tid()
    if bound[h] == 0 or bound_seg_id[h] < 0:
        return
    sid = bound_seg_id[h]
    direction = pos[seg_node_b[sid]] - pos[seg_node_a[sid]]
    if seg_polarity[sid] < 0:
        direction = -direction
    length = wp.length(direction)
    if length <= wp.float64(1.0e-12):
        walk_dir[h] = wp.vec3d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
    else:
        walk_dir[h] = direction / length


@wp.kernel
def _increment_epoch_if_accepted_kernel(
    accepted: wp.array(dtype=wp.int32),
    rng_epoch: wp.array(dtype=wp.int32),
) -> None:
    """Advance the motor RNG epoch only for a committed outer physical transaction."""
    if accepted[0] != 0:
        rng_epoch[0] = rng_epoch[0] + 1


class SegmentMotorRuntime:
    """Production segment-anchor state, force primitive, query accelerator, and accepted KMC commit.

    The runtime owns only GPU arrays after construction.  ``commit_kinetics`` is called from the scheduler's
    final irreversible hook with its authoritative outer-acceptance scalar; rejected attempts may refresh
    derived query scratch but cannot mutate binding, abscissa, or attachment ownership.
    """

    def __init__(
        self,
        head_node: wp.array,
        params: NMIIHandParams,
        seg_node_a: wp.array,
        seg_node_b: wp.array,
        seg_polarity: wp.array,
        max_segment_length_um: float,
        *,
        detach_kinetics: SegmentDetachKinetics = SegmentDetachKinetics.SLIP,
        catch_slip: NMIICatchSlipParams | None = None,
        grid_dim: int = 128,
        device: str | None = None,
    ) -> None:
        if int(seg_node_a.shape[0]) != int(seg_node_b.shape[0]):
            raise ValueError("segment endpoint arrays must have equal length")
        if int(seg_node_a.shape[0]) != int(seg_polarity.shape[0]):
            raise ValueError("segment polarity must have one entry per segment")
        if max_segment_length_um <= 0.0:
            raise ValueError("max_segment_length_um must be positive")
        # Detach-law selection (both mechanistic Bell-Evans).  Default SLIP keeps every existing caller's
        # behaviour bit-identical; CATCH_SLIP is the physiological NMII head-actin bond (Kovacs 2007) and
        # then REQUIRES the two-pathway constants — it is never silently defaulted.
        self.detach_kinetics = SegmentDetachKinetics(detach_kinetics)
        if self.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP and catch_slip is None:
            raise ValueError(
                "SegmentDetachKinetics.CATCH_SLIP requires catch_slip=NMIICatchSlipParams(...) "
                "(k_catch0/x_catch/k_slip0/x_slip/kT are GAP — PI-sourced, never invented)"
            )
        self.catch_slip = catch_slip
        self.device = device
        self.head_node = head_node
        self.params = params
        self.seg_node_a = seg_node_a
        self.seg_node_b = seg_node_b
        self.seg_polarity = seg_polarity
        self.n_heads = int(head_node.shape[0])
        self.n_seg = int(seg_node_a.shape[0])
        self.state = allocate_segment_hand_state(self.n_heads, device=device)
        self.loads = wp.zeros(self.n_heads, dtype=wp.float64, device=device)
        self.loads_full = wp.zeros(self.n_heads, dtype=wp.float64, device=device)
        self.seg_barbed = wp.zeros(self.n_seg, dtype=wp.vec3d, device=device)
        self.capture_radius = float(params.capture_radius)
        # A segment is discoverable through its midpoint whenever the head is within capture of its body.
        self.query_radius = self.capture_radius + 0.5 * float(max_segment_length_um)
        self.query = SegmentQuery(
            seg_node_a,
            seg_node_b,
            self.seg_barbed,
            self.n_heads,
            grid_dim=grid_dim,
            device=device,
        )
        self.accepted_diagnostic = wp.ones(1, dtype=wp.int32, device=device)
        self.rng_epoch = wp.zeros(1, dtype=wp.int32, device=device)

    def refresh_query_geometry(self, pos: wp.array) -> None:
        """Refresh all segment barbed directions once before the accepted-boundary attachment query."""
        wp.launch(
            refresh_segment_barbed_kernel,
            dim=self.n_seg,
            inputs=[pos, self.seg_node_a, self.seg_node_b, self.seg_polarity, self.seg_barbed],
            device=self.device,
        )

    def refresh_bound_geometry(self, pos: wp.array) -> None:
        """Refresh only currently bound heads' live segment directions during inner mechanics."""
        wp.launch(
            _refresh_bound_walk_dir_kernel,
            dim=self.n_heads,
            inputs=[
                pos, self.state["bound"], self.state["seg_id"], self.seg_node_a, self.seg_node_b,
                self.seg_polarity, self.state["walk_dir"],
            ],
            device=self.device,
        )

    def accumulate(self, pos: wp.array, force: wp.array) -> None:
        """Add live segment-anchored crossbridge forces with barycentric endpoint reactions."""
        self.refresh_bound_geometry(pos)
        wp.launch(
            crossbridge_segment_kernel,
            dim=self.n_heads,
            inputs=[
                pos,
                self.head_node,
                self.state["bound"],
                self.state["seg_a"],
                self.state["seg_b"],
                self.state["bary_t"],
                self.state["abscissa"],
                self.state["walk_dir"],
                self.params.k_xb,
                self.params.r0_xb,
                force,
            ],
            device=self.device,
        )

    def compute_loads(self, pos: wp.array) -> None:
        """Refresh live geometry, then compute Hill-tangential and Bell-full crossbridge loads."""
        self.refresh_bound_geometry(pos)
        wp.launch(
            compute_head_loads_segment_kernel,
            dim=self.n_heads,
            inputs=[
                pos,
                self.head_node,
                self.state["bound"],
                self.state["seg_a"],
                self.state["seg_b"],
                self.state["bary_t"],
                self.state["abscissa"],
                self.state["walk_dir"],
                self.params.k_xb,
                self.params.r0_xb,
                self.loads,
                self.loads_full,
            ],
            device=self.device,
        )

    def commit_kinetics(self, pos: wp.array, accepted: wp.array, tau: float, rng_seed: int) -> None:
        """Run load→query→attach→Hill-step/Bell-detach once at the final outer commit boundary."""
        self.compute_loads(pos)
        self.refresh_query_geometry(pos)
        self.query.build(pos, self.query_radius)
        self.query.query_nodes(
            self.head_node,
            pos,
            query_radius=self.query_radius,
            capture_radius=self.capture_radius,
        )
        seed_attach = wp.int32(rng_seed & 0x7FFFFFFF)
        seed_detach = wp.int32((rng_seed ^ 0x5BD1E995) & 0x7FFFFFFF)
        wp.launch(
            attach_segment_gated_kernel,
            dim=self.n_heads,
            inputs=[
                accepted,
                self.rng_epoch,
                self.state["bound"],
                self.state["seg_id"],
                self.state["seg_a"],
                self.state["seg_b"],
                self.state["bary_t"],
                self.state["abscissa"],
                self.state["walk_dir"],
                self.query.seg_id,
                self.query.t,
                self.query.barbed,
                self.seg_node_a,
                self.seg_node_b,
                self.params,
                wp.float64(tau),
                seed_attach,
            ],
            device=self.device,
        )
        if self.detach_kinetics is SegmentDetachKinetics.CATCH_SLIP:
            wp.launch(
                step_detach_segment_gated_catch_slip_kernel,
                dim=self.n_heads,
                inputs=[
                    accepted,
                    self.rng_epoch,
                    self.state["bound"],
                    self.state["seg_id"],
                    self.state["seg_a"],
                    self.state["seg_b"],
                    self.state["bary_t"],
                    self.state["abscissa"],
                    self.state["walk_dir"],
                    self.loads,
                    self.loads_full,
                    self.params,
                    self.catch_slip,
                    wp.float64(tau),
                    seed_detach,
                ],
                device=self.device,
            )
        else:
            wp.launch(
                step_detach_segment_gated_kernel,
                dim=self.n_heads,
                inputs=[
                    accepted,
                    self.rng_epoch,
                    self.state["bound"],
                    self.state["seg_id"],
                    self.state["seg_a"],
                    self.state["seg_b"],
                    self.state["bary_t"],
                    self.state["abscissa"],
                    self.state["walk_dir"],
                    self.loads,
                    self.loads_full,
                    self.params,
                    wp.float64(tau),
                    seed_detach,
                ],
                device=self.device,
            )
        wp.launch(
            _increment_epoch_if_accepted_kernel,
            dim=1,
            inputs=[accepted, self.rng_epoch],
            device=self.device,
        )


# ── NumPy host oracle (bit-for-formula twin of crossbridge_segment_kernel) ────────────────────────
def crossbridge_segment_reference(
    head_pos: npt.NDArray, pos_a: npt.NDArray, pos_b: npt.NDArray, bary_t: float,
    abscissa: float, walk_dir: npt.NDArray, k_xb: float, r0_xb: float = 0.0,
) -> dict[str, npt.NDArray]:
    """Return ``{f_head, f_seg_a, f_seg_b, x_att, load_tangential, load_full}`` for one segment-anchored head.

    The barycentric reaction split ``f_seg_a = −(1−t)·f``, ``f_seg_b = −t·f`` satisfies Newton's 3rd law
    (``f_head + f_seg_a + f_seg_b = 0``) — the acceptance the host gate checks.
    """
    t = float(bary_t)
    x_att = (1.0 - t) * pos_a + t * pos_b + abscissa * walk_dir
    d = x_att - head_pos
    # DIRECTIONAL (tangential) crossbridge (2026-07-23 defect#3 fix) — the host twin of crossbridge_segment_kernel:
    # force acts along the unit walk_dir with magnitude = the tangential load; |F| == |tangential load|.
    load_tangential = k_xb * (float(np.dot(d, walk_dir)) - r0_xb)
    f = load_tangential * np.asarray(walk_dir, dtype=np.float64)
    return {
        "f_head": f, "f_seg_a": -(1.0 - t) * f, "f_seg_b": -t * f, "x_att": x_att,
        "load_tangential": load_tangential,
        "load_full": abs(load_tangential),
    }


def crossbridge_segment_split_reference(
    actuator_pos: npt.NDArray, port_pos: npt.NDArray, head_node: int, seg_a: int, seg_b: int,
    bary_t: float, abscissa: float, walk_dir: npt.NDArray, k_xb: float, r0_xb: float = 0.0,
) -> dict[str, npt.NDArray]:
    """Two-array split oracle: head force scatters into ``actuator_force``, segment reactions into ``port_force``.

    Proves the relayout is physics-preserving: the returned per-array scatters are numerically identical to the
    single-array :func:`crossbridge_segment_reference` (the head reads ``actuator_pos[head_node]`` and the
    segment reads ``port_pos[seg_a/seg_b]`` — arrays that need NOT be co-located), and Newton's 3rd law holds
    across the two arrays (``f_head + f_seg_a + f_seg_b = 0``).  The result buffers show where each force lands.
    """
    combined = crossbridge_segment_reference(
        actuator_pos[int(head_node)], port_pos[int(seg_a)], port_pos[int(seg_b)],
        bary_t, abscissa, walk_dir, k_xb, r0_xb,
    )
    actuator_force = np.zeros_like(actuator_pos, dtype=np.float64)
    port_force = np.zeros_like(port_pos, dtype=np.float64)
    actuator_force[int(head_node)] += combined["f_head"]
    port_force[int(seg_a)] += combined["f_seg_a"]
    port_force[int(seg_b)] += combined["f_seg_b"]
    return {
        "actuator_force": actuator_force, "port_force": port_force,
        "f_head": combined["f_head"], "f_seg_a": combined["f_seg_a"], "f_seg_b": combined["f_seg_b"],
        "x_att": combined["x_att"], "load_tangential": combined["load_tangential"],
        "load_full": combined["load_full"],
    }


def r0_bind_at_attach_reference(
    actuator_pos: npt.NDArray, port_pos: npt.NDArray, head_node: int, seg_a: int, seg_b: int,
    bary_t: float, walk_dir: npt.NDArray,
) -> float:
    """Host twin of the attach-time capture: ``r0_bind = dot(x_att − x_head, walk_dir)`` with abscissa=0.

    ``x_att = (1−t)·port_pos[seg_a] + t·port_pos[seg_b]`` (the barycentric attachment at the binding instant).
    Using this as the crossbridge rest makes the passive force exactly 0 at the instant of binding, regardless
    of the head-to-segment placement offset on a coarse mesh.
    """
    t = float(bary_t)
    x_att = (1.0 - t) * port_pos[int(seg_a)] + t * port_pos[int(seg_b)]
    return float(np.dot(x_att - actuator_pos[int(head_node)], np.asarray(walk_dir, dtype=np.float64)))


def crossbridge_segment_split_r0bind_reference(
    actuator_pos: npt.NDArray, port_pos: npt.NDArray, head_node: int, seg_a: int, seg_b: int,
    bary_t: float, abscissa: float, walk_dir: npt.NDArray, k_xb: float, r0_bind: float,
) -> dict[str, npt.NDArray]:
    """Two-array split oracle with a PER-HEAD zero-strain reference ``r0_bind`` (attach-unstrained).

    Bit-for-formula twin of :func:`crossbridge_segment_split_r0bind_kernel`: the scalar ``r0_xb`` of
    :func:`crossbridge_segment_split_reference` is replaced by ``r0_bind`` (the attach-time projection).  With
    ``r0_bind = r0_bind_at_attach_reference(...)`` and ``abscissa=0`` the returned force is the zero vector
    (unstrained at bind); a positive ``abscissa`` (the power stroke) then yields ``|f| = k_xb·abscissa`` along
    ``walk_dir``.  Newton's 3rd law still holds across the two arrays (``f_head + f_seg_a + f_seg_b = 0``).
    """
    t = float(bary_t)
    x_att = (1.0 - t) * port_pos[int(seg_a)] + t * port_pos[int(seg_b)] + abscissa * walk_dir
    d = x_att - actuator_pos[int(head_node)]
    load_tangential = k_xb * (float(np.dot(d, walk_dir)) - float(r0_bind))
    f = load_tangential * np.asarray(walk_dir, dtype=np.float64)
    actuator_force = np.zeros_like(actuator_pos, dtype=np.float64)
    port_force = np.zeros_like(port_pos, dtype=np.float64)
    actuator_force[int(head_node)] += f
    port_force[int(seg_a)] += -(1.0 - t) * f
    port_force[int(seg_b)] += -t * f
    return {
        "actuator_force": actuator_force, "port_force": port_force,
        "f_head": f, "f_seg_a": -(1.0 - t) * f, "f_seg_b": -t * f, "x_att": x_att,
        "load_tangential": load_tangential, "load_full": abs(load_tangential),
    }
