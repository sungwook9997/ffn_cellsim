r"""Host + device gates for the I4-weave actin↔motor JOINT (WovenCell → production point-to-segment NMII).

The joint (:mod:`aleph.components.weave.woven_cell`) turns the ONE woven actin network into the three device arrays
the production motor binds to — ``seg_node_a`` / ``seg_node_b`` / ``seg_polarity`` — so a head grabs a point on
a REAL filament and its ``walk_dir`` is overwritten from that segment's barbed polarity (not the bipolar-
geometry default, which is the OFF/regression path). These gates certify:

  1. topology — flattened adjacent-node segments, one per actin bond, polarity broadcast per fiber;
  2. barbed-end DIRECTION correctness — the joint orients every segment toward its filament barbed end, and
     agrees with the node-anchored :mod:`aleph.components.weave.walk_dir` polarity (sign always; exact on a straight
     fiber). A polarity flip flips every one of a fiber's segment directions (the reattach hand-off);
  3. reaction symmetry + transmission — a head bound to a woven segment closes Newton's 3rd law across the two
     nodes, and a collinear isometric hold transmits the crossbridge tension 1:1 (Hill == Bell == applied);
  4. (CUDA-gated, native) end-to-end: the joint feeds a ``SegmentMotorRuntime``; an accepted commit binds a
     head with ``walk_dir`` == the woven segment's barbed direction — the polarity overwrite on device.

Dev-Mac I0-A: the CUDA kernels are not launched here (host oracles carry the same formulas); the device gate
skips without a GPU and runs on the lead's gbook A5000.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from aleph.laws.architecture_spec import CORTEX
from aleph.components.motor.segment_motor import crossbridge_segment_reference
from aleph.components.weave.regions import RegionSpec
from aleph.components.weave.walk_dir import barbed_end_node, walk_dir_from_polarity
from aleph.components.weave.woven_cell import (
    actin_segment_topology,
    segment_barbed_directions,
    weave_cell,
)


def _small_cortex(n: int = 120) -> RegionSpec:
    fil = dataclasses.replace(CORTEX.filament, n_filaments=n)
    return RegionSpec(arch=dataclasses.replace(CORTEX, filament=fil))


# ── 1. topology ──────────────────────────────────────────────────────────────────────────────────
def test_topology_flattens_adjacent_nodes_and_broadcasts_polarity() -> None:
    """Two fibers (3 + 4 nodes) → 2 + 3 adjacent-node segments; barbed sign broadcast per fiber."""
    seg_a, seg_b, seg_pol = actin_segment_topology(
        np.array([0, 3, 7], dtype=np.int64), np.array([1, -1], dtype=np.int64))
    assert np.array_equal(seg_a, np.array([0, 1, 3, 4, 5], dtype=np.int32))
    assert np.array_equal(seg_b, seg_a + 1)                     # always adjacent nodes
    assert np.array_equal(seg_pol, np.array([1, 1, -1, -1, -1], dtype=np.int32))


def test_topology_rejects_singleton_fiber_and_bad_polarity() -> None:
    with pytest.raises(ValueError, match="two contiguous nodes"):
        actin_segment_topology(np.array([0, 1, 4], np.int64), np.array([1, 1], np.int64))
    with pytest.raises(ValueError, match="one entry per fiber"):
        actin_segment_topology(np.array([0, 3, 7], np.int64), np.array([1], np.int64))


def test_wovencell_method_matches_module_and_counts_segments() -> None:
    """The WovenCell method is the module joint applied to the network; segment count = Σ(len_fiber − 1)."""
    wc = weave_cell([_small_cortex()], rng=np.random.default_rng(0))
    seg_a, seg_b, seg_pol = wc.actin_segment_topology()
    ref = actin_segment_topology(wc.fiber_offsets, wc.polarity)
    for got, exp in zip((seg_a, seg_b, seg_pol), ref):
        assert np.array_equal(got, exp)
    assert seg_a.shape[0] == int(np.sum(np.diff(wc.fiber_offsets) - 1))
    # each segment's two nodes belong to the SAME fiber (a segment never spans a fiber boundary)
    assert np.array_equal(wc.node_fiber[seg_a], wc.node_fiber[seg_b])


# ── 2. barbed-end direction correctness (the polarity gate) ───────────────────────────────────────
def test_segment_barbed_directions_are_unit_and_reproduce_the_device_formula() -> None:
    """Host oracle == refresh_segment_barbed_kernel formula: sign·(b−a) normalised; |dir| = 1."""
    pos = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 0, 0], [0, 3, 0]], float)
    seg_a = np.array([0, 1, 3], np.int32)
    seg_b = np.array([1, 2, 4], np.int32)
    seg_pol = np.array([1, -1, 1], np.int32)
    got = segment_barbed_directions(pos, seg_a, seg_b, seg_pol)
    expect = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0]], float)   # seg 1 flipped by polarity −1
    np.testing.assert_allclose(got, expect, atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(got, axis=1), 1.0, atol=1e-12)


def test_degenerate_zero_length_segment_is_passive() -> None:
    """A zero-length segment yields a zero direction (a passive head — no directed force)."""
    pos = np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]], float)
    got = segment_barbed_directions(pos, np.array([0], np.int32), np.array([1], np.int32), np.array([1], np.int32))
    assert np.allclose(got[0], 0.0)


def test_barbed_direction_points_toward_the_filament_barbed_end_on_the_woven_cell() -> None:
    """Every woven segment's barbed dir has POSITIVE projection on the node→barbed-end direction (correct sign).

    Ties the segment-anchored polarity (this joint) to the node-anchored walk_dir oracle: a head that grabs the
    segment's first node and walks the segment barbed dir must move toward that fiber's barbed end.
    """
    wc = weave_cell([_small_cortex()], rng=np.random.default_rng(1))
    seg_a, seg_b, seg_pol = wc.actin_segment_topology()
    seg_dir = wc.segment_barbed_directions()
    barbed = barbed_end_node(wc.fiber_offsets, wc.polarity)        # (F,) global barbed node per fiber
    fib = wc.node_fiber[seg_a]
    for s in range(seg_a.shape[0]):
        node_dir = walk_dir_from_polarity(int(seg_a[s]), int(barbed[fib[s]]), wc.pos)
        if np.linalg.norm(node_dir) < 1e-9:      # anchor coincides with the barbed end → no reference direction
            continue
        assert float(np.dot(seg_dir[s], node_dir)) > 0.0, f"segment {s} points away from its barbed end"


def test_straight_fiber_segment_dir_equals_node_anchored_walk_dir_exactly() -> None:
    """On a STRAIGHT fiber the local segment barbed dir == the global node→barbed walk_dir (they coincide)."""
    pos = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], float)     # one straight fiber, 4 nodes
    offsets = np.array([0, 4], np.int64)
    for pol, expect in ((np.array([1]), np.array([1.0, 0, 0])), (np.array([-1]), np.array([-1.0, 0, 0]))):
        seg_a, seg_b, seg_pol = actin_segment_topology(offsets, pol)
        seg_dir = segment_barbed_directions(pos, seg_a, seg_b, seg_pol)
        barbed = barbed_end_node(offsets, pol)
        for s in range(seg_a.shape[0]):
            np.testing.assert_allclose(seg_dir[s], expect, atol=1e-12)   # segment-local dir always defined
            node_dir = walk_dir_from_polarity(int(seg_a[s]), int(barbed[0]), pos)
            if np.linalg.norm(node_dir) < 1e-9:      # anchor == barbed end (polarity −1, first segment): skip
                continue
            np.testing.assert_allclose(seg_dir[s], node_dir, atol=1e-12)


def test_polarity_flip_flips_every_segment_direction_reattach_handoff() -> None:
    """Rebinding to a filament of OPPOSITE polarity flips walk_dir — the joint tracks the actual actin."""
    pos = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], float)
    a, b, pol_p = actin_segment_topology(np.array([0, 3], np.int64), np.array([1], np.int64))
    _, _, pol_m = actin_segment_topology(np.array([0, 3], np.int64), np.array([-1], np.int64))
    dir_p = segment_barbed_directions(pos, a, b, pol_p)
    dir_m = segment_barbed_directions(pos, a, b, pol_m)
    np.testing.assert_allclose(dir_m, -dir_p, atol=1e-12)


# ── 3. reaction symmetry + transmission at the joint ──────────────────────────────────────────────
def test_bound_head_closes_newtons_third_law_across_the_woven_segment() -> None:
    """A head bound to a woven segment with walk_dir from the joint: f_head + f_seg_a + f_seg_b = 0 exactly."""
    pos = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], float)
    seg_a, seg_b, seg_pol = actin_segment_topology(np.array([0, 3], np.int64), np.array([1], np.int64))
    seg_dir = segment_barbed_directions(pos, seg_a, seg_b, seg_pol)
    head = np.array([0.5, 0.05, 0.0])                       # near the first segment, slightly off-axis
    r = crossbridge_segment_reference(head, pos[seg_a[0]], pos[seg_b[0]], bary_t=0.5,
                                      abscissa=0.03, walk_dir=seg_dir[0], k_xb=1000.0)
    net = r["f_head"] + r["f_seg_a"] + r["f_seg_b"]
    assert np.abs(net).max() < 1e-9


def test_collinear_isometric_hold_transmits_one_to_one() -> None:
    """Collinear head displacement along the woven barbed dir: Hill load == Bell load == k_xb·|Δ| (transmission 1)."""
    pos = np.array([[0, 0, 0], [1, 0, 0]], float)
    seg_a, seg_b, seg_pol = actin_segment_topology(np.array([0, 2], np.int64), np.array([1], np.int64))
    seg_dir = segment_barbed_directions(pos, seg_a, seg_b, seg_pol)   # +x
    k_xb, delta = 1000.0, 0.02
    # head sits `delta` behind the attachment along −walk_dir so the crossbridge is a pure +walk_dir stretch
    x_att = pos[seg_a[0]].astype(float)                               # bary_t=0, abscissa=0 → attach at node a
    head = x_att - delta * seg_dir[0]
    r = crossbridge_segment_reference(head, pos[seg_a[0]], pos[seg_b[0]], bary_t=0.0,
                                      abscissa=0.0, walk_dir=seg_dir[0], k_xb=k_xb, r0_xb=0.0)
    assert r["load_tangential"] == pytest.approx(k_xb * delta, rel=1e-9)
    assert r["load_full"] == pytest.approx(r["load_tangential"], rel=1e-9)     # transmission ratio = 1
    # the head is pulled toward the barbed (+walk_dir) direction — a directed contractile pull on actin
    assert float(np.dot(r["f_head"], seg_dir[0])) > 0.0


# ── 4. CUDA-gated end-to-end: the joint feeds the production runtime, polarity overwritten on device ──
def _cuda_device() -> str | None:
    import warp as wp

    wp.init()
    device = next((item for item in wp.get_devices() if item.is_cuda), None)
    return str(device) if device is not None else None


def test_cuda_woven_joint_binds_head_with_segment_polarity_end_to_end() -> None:
    """The WovenCell joint → SegmentMotorRuntime; an accepted commit binds a head whose walk_dir == the woven
    segment's barbed direction (the I4-weave polarity overwrite, proven on device)."""
    device = _cuda_device()
    if device is None:
        pytest.skip("I0-A: end-to-end woven-joint bind gate requires a CUDA GPU")
    import warp as wp

    from aleph.components.motor.hand import NMIIHandParams
    from aleph.components.motor.segment_motor import SegmentMotorRuntime

    # one straight actin fiber (nodes 0..2, barbed at the last node) + one head node (3) just off its midpoint
    node_pos = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [0.5, 0.01, 0.0]], np.float64)
    offsets = np.array([0, 3], np.int64)
    seg_a, seg_b, seg_pol = actin_segment_topology(offsets, np.array([1], np.int64))
    barbed_expected = segment_barbed_directions(node_pos, seg_a, seg_b, seg_pol)[0]   # +x, host oracle

    pos = wp.array(node_pos, dtype=wp.vec3d, device=device)
    head_node = wp.array(np.array([3], np.int32), dtype=wp.int32, device=device)
    params = NMIIHandParams()
    params.k_on = 1.0e6         # attach with probability ≈ 1 this tick
    params.k_off0 = 0.0
    params.f0 = 1.0
    params.v0 = 0.0             # no walk → isolate the binding/polarity overwrite
    params.f_stall = 1.0
    params.kappa = 1.0
    params.k_xb = 1.0
    params.r0_head = 0.0
    params.r0_xb = 0.0
    params.capture_radius = 0.05
    runtime = SegmentMotorRuntime(
        head_node, params,
        wp.array(seg_a, dtype=wp.int32, device=device),
        wp.array(seg_b, dtype=wp.int32, device=device),
        wp.array(seg_pol, dtype=wp.int32, device=device),
        max_segment_length_um=1.0, grid_dim=8, device=device)
    accepted = wp.ones(1, dtype=wp.int32, device=device)
    runtime.commit_kinetics(pos, accepted, 1.0, 7)
    assert int(runtime.state["bound"].numpy()[0]) == 1
    np.testing.assert_allclose(runtime.state["walk_dir"].numpy()[0], barbed_expected, atol=1e-9)
