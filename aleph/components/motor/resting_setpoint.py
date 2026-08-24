r"""Resting bound-myosin SETPOINT — the fine-grained-faithful resting cortical-tension source (candidate 1).

The resting-baseline diagnosis (``docs/v2_audit/RESTING_BASELINE_DIAGNOSIS_2026-07-23c.md``) proved that the
resting membrane turgor (``ΔP·A_node`` per membrane node) can ONLY be balanced by CORTICAL TENSION transmitted
through the ERM, and that neither passive crosslink pre-straining (§3, non-conservative) nor a force-free ERM
(§2, transmits nothing) supplies it.  The physiological source of that resting tension is BOUND MYOSIN: at rest
a fraction of the NMII heads are engaged and carry an isometric load, loading the actin network with active
hoop tension.  This module seeds exactly that state onto the EXISTING head-resolved Stam-Hocky motor
(``ac.motor.segment_motor.SegmentMotorRuntime``) — it is NOT a lumped tension term.

Faithful, not lumped:
  * A seeded head is a real bound crossbridge in the production segment runtime (``bound``/``seg_a``/``seg_b``/
    ``bary_t``/``abscissa``/``walk_dir``).  Its force is the SAME ``crossbridge_segment_kernel`` power-stroke
    spring every running motor uses; ``compute_head_loads_segment_kernel`` reads back the same tangential/Bell
    loads.  Nothing here adds a new force law.
  * The resting per-head tension is realized THROUGH the power stroke: a head bound at the nearest point on its
    actin segment has its head→anchor offset PERPENDICULAR to the walk direction (nearest-point property), so
    the tangential crossbridge load is exactly ``k_xb·(abscissa − r0_xb) + k_xb·⟨offset, ŵ⟩`` with the offset
    term ≈ 0.  Setting ``abscissa = f_head/k_xb + r0_xb − ⟨offset, ŵ⟩`` makes each seeded head carry the PI-GAP
    per-head force ``f_head`` tangentially, by construction and verifiable against the reference kernel.

The two magnitudes are PI-GAPs, surfaced never invented:
  * ``fraction`` — the resting bound-head fraction (duty ratio at rest).  Absent from the Contract-Graph.
  * ``per_head_force_pn`` — the resting isometric per-head tension.  Absent from the Contract-Graph.

:class:`RestingBoundMyosinSetpoint` validates a complete, source-labelled contract exactly like
:class:`~aleph.components.incumbent.erm_tether.ERMBellKinetics`; the cell build refuses a partial/defaulted setpoint and,
when the setpoint is UNSET, leaves the motor at its unbound-at-rest t0 and surfaces the PI-GAP.  The physiological
values do NOT live here; a TEST fraction proves the MECHANISM, never the native gate.

Host planning (:func:`plan_resting_bound_heads`) is pure NumPy and out-of-hot-loop — a build-time construction
like the ERM/LINC pairing — so it is dev-Mac testable.  :func:`apply_resting_bound_heads` writes the plan into
the device segment-hand state once (host→device at build; no authoritative per-step host state thereafter, I0-A).

Units: FF µm·pN·s.  ``k_xb`` [pN/µm], ``f_head`` [pN], ``abscissa`` [µm].
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from aleph.components.motor.segment_query import closest_point_on_segment

__all__ = [
    "RestingBoundMyosinSetpoint",
    "RestingBoundPlan",
    "plan_resting_bound_heads",
    "resting_tangential_load_reference",
    "apply_resting_bound_heads",
]


@dataclass(frozen=True, slots=True)
class RestingBoundMyosinSetpoint:
    """Complete, source-gated resting bound-myosin contract.  No field has a physiological default.

    Attributes:
        fraction: resting bound-head fraction (duty ratio at rest), ``0 < fraction ≤ 1``.  PI GAP.
        per_head_force_pn: resting isometric per-head tangential tension [pN], finite and positive.  PI GAP.
        source: literature/provenance label; production use requires it to cover BOTH magnitudes.
        capture_radius_um: bind only heads whose nearest actin segment is within this reach [µm].  When
            ``None`` the motor's own ``params.capture_radius`` is used (the same attach gate the KMC uses).
    """

    fraction: float
    per_head_force_pn: float
    source: str
    capture_radius_um: float | None = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.fraction) or not (0.0 < self.fraction <= 1.0):
            raise ValueError("resting bound-myosin fraction must be finite in (0, 1]")
        if not math.isfinite(self.per_head_force_pn) or self.per_head_force_pn <= 0.0:
            raise ValueError("resting bound-myosin per-head force must be finite and positive [pN]")
        if not self.source.strip():
            raise ValueError("resting bound-myosin setpoint requires an explicit literature/provenance source")
        if self.capture_radius_um is not None and (
            not math.isfinite(self.capture_radius_um) or self.capture_radius_um <= 0.0
        ):
            raise ValueError("resting bound-myosin capture_radius_um, when given, must be finite and positive")


@dataclass(frozen=True, slots=True)
class RestingBoundPlan:
    """Host plan for the heads bound at rest (all arrays indexed by the K selected heads).

    Attributes:
        head_index: (K,) local head indices (0..n_heads) to mark bound.
        seg_id: (K,) the bound segment index (into the seg_node_a/seg_node_b tables).
        seg_a / seg_b: (K,) the bound segment's two GLOBAL actin node indices.
        bary_t: (K,) barycentric fraction of the attachment along the segment (0=a, 1=b).
        walk_dir: (K, 3) unit barbed-end walk direction of the bound segment.
        abscissa: (K,) walked power-stroke displacement [µm] realizing ``f_head`` tangentially.
        diagnostics: build-time bookkeeping (counts, realized loads) — never enters the runtime state.
    """

    head_index: npt.NDArray[np.int64]
    seg_id: npt.NDArray[np.int64]
    seg_a: npt.NDArray[np.int64]
    seg_b: npt.NDArray[np.int64]
    bary_t: npt.NDArray[np.float64]
    walk_dir: npt.NDArray[np.float64]
    abscissa: npt.NDArray[np.float64]
    diagnostics: dict

    @property
    def n_bound(self) -> int:
        return int(self.head_index.shape[0])


def _segment_barbed_directions(
    seg_a_pos: npt.NDArray[np.float64],
    seg_b_pos: npt.NDArray[np.float64],
    seg_polarity: npt.NDArray[np.int64],
) -> npt.NDArray[np.float64]:
    """Per-segment unit barbed direction (matches ``refresh_segment_barbed_kernel``: flip when polarity<0)."""
    direction = seg_b_pos - seg_a_pos
    flip = np.asarray(seg_polarity, dtype=np.int64) < 0
    direction[flip] = -direction[flip]
    norm = np.linalg.norm(direction, axis=1, keepdims=True)
    out = np.zeros_like(direction)
    ok = norm[:, 0] > 1e-12
    out[ok] = direction[ok] / norm[ok]
    return out


def _nearest_segment_per_head(
    head_pos: npt.NDArray[np.float64],
    seg_a_pos: npt.NDArray[np.float64],
    seg_b_pos: npt.NDArray[np.float64],
    *,
    k_candidates: int = 12,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Nearest actin segment for every head → (seg_id, bary_t, point-to-segment dist).

    A KD-tree on segment midpoints selects the ``k_candidates`` nearest bodies per head; the true
    point-to-segment nearest among them is then found with the canonical ``closest_point_on_segment`` (the same
    geometry the device ``SegmentQuery`` uses).  Falls back to the exact brute-force scan on small inputs / when
    SciPy is unavailable, so the result is identical to :func:`segment_query.point_to_segment_query_reference`.
    """
    n_head = head_pos.shape[0]
    n_seg = seg_a_pos.shape[0]
    seg_id = np.full(n_head, -1, np.int64)
    bary_t = np.zeros(n_head, np.float64)
    dist = np.full(n_head, np.inf, np.float64)
    if n_head == 0 or n_seg == 0:
        return seg_id, bary_t, dist

    midpoints = 0.5 * (seg_a_pos + seg_b_pos)
    use_bruteforce = n_seg <= 64
    tree = None
    if not use_bruteforce:
        try:
            from scipy.spatial import cKDTree

            tree = cKDTree(midpoints)
        except Exception:  # noqa: BLE001 — deterministic brute-force fallback (no SciPy)
            use_bruteforce = n_head * n_seg <= 20_000_000
            if not use_bruteforce:
                raise RuntimeError(
                    "SciPy cKDTree is required to plan a native-scale resting bound-myosin setpoint"
                ) from None

    for h in range(n_head):
        if use_bruteforce:
            cand = range(n_seg)
        else:
            kk = min(k_candidates, n_seg)
            _, idx = tree.query(head_pos[h], k=kk)
            cand = np.atleast_1d(idx).astype(np.int64)
        best_d, best_t, best_s = np.inf, 0.0, -1
        for s in cand:
            t, _, d = closest_point_on_segment(head_pos[h], seg_a_pos[s], seg_b_pos[s])
            if d < best_d:
                best_d, best_t, best_s = d, t, int(s)
        seg_id[h] = best_s
        bary_t[h] = best_t
        dist[h] = best_d
    return seg_id, bary_t, dist


def plan_resting_bound_heads(
    pos: npt.NDArray[np.float64],
    head_node: npt.NDArray[np.int64],
    seg_node_a: npt.NDArray[np.int64],
    seg_node_b: npt.NDArray[np.int64],
    seg_polarity: npt.NDArray[np.int64],
    setpoint: RestingBoundMyosinSetpoint,
    *,
    k_xb: float,
    r0_xb: float,
    capture_radius_um: float,
) -> RestingBoundPlan:
    """Plan which heads bind at rest and with what power-stroke abscissa (pure NumPy, deterministic).

    Selection is RNG-free: the ``round(fraction·n_heads)`` targets are strided evenly over the ELIGIBLE heads
    (those with a nearest actin segment within ``capture_radius`` AND a physical, non-negative abscissa), so the
    plan is reproducible without a seed.  When fewer heads are eligible than requested, every eligible head is
    bound and the shortfall is recorded — the fraction is never faked up.

    Each bound head's abscissa is set so its TANGENTIAL crossbridge load equals ``per_head_force_pn`` exactly
    (verified by :func:`resting_tangential_load_reference`).  Because binding uses the nearest point on the
    segment, the head→anchor offset is ~perpendicular to the walk direction, so the abscissa is ≈ ``f_head/k_xb``
    — a physical (non-negative) power-stroke displacement, not a construction pre-strain of the actin bonds.

    Args:
        pos: (N, 3) all node positions [µm] (actin + myosin + …); global-indexed.
        head_node: (n_heads,) GLOBAL indices of the myosin head particles.
        seg_node_a / seg_node_b: (S,) GLOBAL actin node indices of each segment's endpoints.
        seg_polarity: (S,) ±1 barbed-direction sign per segment.
        setpoint: the source-gated resting contract (fraction + per-head force).
        k_xb: crossbridge stiffness [pN/µm] (the motor's ``params.k_xb``).
        r0_xb: crossbridge rest length [µm] (the motor's ``params.r0_xb``, ~0).
        capture_radius_um: bind only heads within this point-to-segment reach [µm].

    Returns:
        A :class:`RestingBoundPlan` (K bound heads) + a diagnostics dict.
    """
    pos = np.ascontiguousarray(pos, np.float64)
    head_node = np.ascontiguousarray(head_node, np.int64)
    seg_node_a = np.ascontiguousarray(seg_node_a, np.int64)
    seg_node_b = np.ascontiguousarray(seg_node_b, np.int64)
    seg_polarity = np.ascontiguousarray(seg_polarity, np.int64)
    n_heads = int(head_node.shape[0])
    k_xb = float(k_xb)
    r0_xb = float(r0_xb)
    if k_xb <= 0.0:
        raise ValueError("k_xb must be positive to realize a resting per-head tension")

    empty_plan = RestingBoundPlan(
        head_index=np.zeros(0, np.int64), seg_id=np.zeros(0, np.int64), seg_a=np.zeros(0, np.int64),
        seg_b=np.zeros(0, np.int64), bary_t=np.zeros(0, np.float64), walk_dir=np.zeros((0, 3), np.float64),
        abscissa=np.zeros(0, np.float64),
        diagnostics={"n_heads": n_heads, "n_target": 0, "n_eligible": 0, "n_bound": 0},
    )
    if n_heads == 0 or seg_node_a.shape[0] == 0:
        return empty_plan

    head_pos = pos[head_node]
    seg_a_pos = pos[seg_node_a]
    seg_b_pos = pos[seg_node_b]
    barbed = _segment_barbed_directions(seg_a_pos, seg_b_pos, seg_polarity)

    seg_id, bary_t, dist = _nearest_segment_per_head(head_pos, seg_a_pos, seg_b_pos)

    f_head = float(setpoint.per_head_force_pn)
    # attachment point on the nearest segment (abscissa=0) and its walk direction
    within = (seg_id >= 0) & (dist <= float(capture_radius_um))
    seg_pt = np.zeros((n_heads, 3), np.float64)
    wdir = np.zeros((n_heads, 3), np.float64)
    valid_seg = seg_id.copy()
    valid_seg[~within] = 0
    seg_pt = (1.0 - bary_t)[:, None] * seg_a_pos[valid_seg] + bary_t[:, None] * seg_b_pos[valid_seg]
    wdir = barbed[valid_seg]
    # tangential geometric offset (≈0 by the nearest-point perpendicularity, exact for interior t)
    geo = np.einsum("ij,ij->i", seg_pt - head_pos, wdir)
    abscissa = f_head / k_xb + r0_xb - geo
    walk_ok = np.linalg.norm(wdir, axis=1) > 1e-12
    eligible = within & walk_ok & (abscissa >= 0.0)

    elig_idx = np.nonzero(eligible)[0]
    n_eligible = int(elig_idx.shape[0])
    n_target = int(np.floor(setpoint.fraction * n_heads + 0.5))
    if n_target > 0 and n_eligible > 0:
        if n_eligible <= n_target:
            chosen = elig_idx
        else:  # evenly strided midpoint sampling over the eligible heads (deterministic, seed-free)
            pick = np.floor((np.arange(n_target, dtype=np.float64) + 0.5) * n_eligible / n_target).astype(np.int64)
            chosen = elig_idx[np.unique(pick)]
    else:
        chosen = np.zeros(0, np.int64)

    plan = RestingBoundPlan(
        head_index=chosen.astype(np.int64),
        seg_id=seg_id[chosen].astype(np.int64),
        seg_a=seg_node_a[seg_id[chosen]].astype(np.int64),
        seg_b=seg_node_b[seg_id[chosen]].astype(np.int64),
        bary_t=bary_t[chosen].astype(np.float64),
        walk_dir=np.ascontiguousarray(wdir[chosen], np.float64),
        abscissa=np.ascontiguousarray(abscissa[chosen], np.float64),
        diagnostics={
            "n_heads": n_heads,
            "n_target": n_target,
            "n_eligible": n_eligible,
            "n_bound": int(chosen.shape[0]),
            "requested_fraction": float(setpoint.fraction),
            "realized_fraction": float(chosen.shape[0]) / float(n_heads) if n_heads else 0.0,
            "per_head_force_pn": f_head,
            "k_xb_pn_per_um": k_xb,
            "abscissa_mean_um": float(abscissa[chosen].mean()) if chosen.size else 0.0,
            "abscissa_max_um": float(abscissa[chosen].max()) if chosen.size else 0.0,
            "capture_radius_um": float(capture_radius_um),
            "eligible_shortfall": bool(n_eligible < n_target),
        },
    )
    return plan


def resting_tangential_load_reference(
    pos: npt.NDArray[np.float64],
    head_node: npt.NDArray[np.int64],
    plan: RestingBoundPlan,
    *,
    k_xb: float,
    r0_xb: float,
) -> npt.NDArray[np.float64]:
    """Per-bound-head tangential crossbridge load [pN] for the seeded plan (matches the device load kernel).

    Reuses the ``compute_head_loads_segment_kernel`` formula ``k_xb·(⟨x_att − head, ŵ⟩ − r0_xb)`` with
    ``x_att = (1−t)·a + t·b + abscissa·ŵ``.  Every entry must equal ``per_head_force_pn`` to numerical
    tolerance — the acceptance that the abscissa seed realized the requested resting tension.
    """
    if plan.n_bound == 0:
        return np.zeros(0, np.float64)
    head_pos = np.ascontiguousarray(pos, np.float64)[np.asarray(head_node, np.int64)[plan.head_index]]
    a_pos = np.ascontiguousarray(pos, np.float64)[plan.seg_a]
    b_pos = np.ascontiguousarray(pos, np.float64)[plan.seg_b]
    x_att = (1.0 - plan.bary_t)[:, None] * a_pos + plan.bary_t[:, None] * b_pos + plan.abscissa[:, None] * plan.walk_dir
    d = x_att - head_pos
    return float(k_xb) * (np.einsum("ij,ij->i", d, plan.walk_dir) - float(r0_xb))


def apply_resting_bound_heads(state: dict, plan: RestingBoundPlan, device: str | None = None) -> int:
    """Write the resting bound-head plan into a device segment-hand ``state`` (host→device once, at build).

    Overwrites the ``bound``/``seg_id``/``seg_a``/``seg_b``/``bary_t``/``abscissa``/``walk_dir`` device arrays
    of :func:`ac.motor.segment_motor.allocate_segment_hand_state` with the seeded state.  All non-selected heads
    stay unbound (the allocate default).  Returns the number of heads bound.

    Warp is imported lazily so the pure-NumPy planner/reference above stay importable on a non-CUDA host.
    """
    import warp as wp

    n_heads = int(state["bound"].shape[0])
    bound = np.zeros(n_heads, np.int32)
    seg_id = np.full(n_heads, -1, np.int32)
    seg_a = np.full(n_heads, -1, np.int32)
    seg_b = np.full(n_heads, -1, np.int32)
    bary_t = np.zeros(n_heads, np.float64)
    abscissa = np.zeros(n_heads, np.float64)
    walk_dir = np.zeros((n_heads, 3), np.float64)
    if plan.n_bound:
        hi = plan.head_index
        bound[hi] = 1
        seg_id[hi] = plan.seg_id.astype(np.int32)
        seg_a[hi] = plan.seg_a.astype(np.int32)
        seg_b[hi] = plan.seg_b.astype(np.int32)
        bary_t[hi] = plan.bary_t
        abscissa[hi] = plan.abscissa
        walk_dir[hi] = plan.walk_dir
    with wp.ScopedDevice(device):
        state["bound"].assign(np.ascontiguousarray(bound, np.int32))
        state["seg_id"].assign(np.ascontiguousarray(seg_id, np.int32))
        state["seg_a"].assign(np.ascontiguousarray(seg_a, np.int32))
        state["seg_b"].assign(np.ascontiguousarray(seg_b, np.int32))
        state["bary_t"].assign(np.ascontiguousarray(bary_t, np.float64))
        state["abscissa"].assign(np.ascontiguousarray(abscissa, np.float64))
        state["walk_dir"].assign(np.ascontiguousarray(walk_dir, np.float64))
    return int(plan.n_bound)
