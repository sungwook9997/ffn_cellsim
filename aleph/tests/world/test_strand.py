"""Strand composition — the N / N-1 / N-2 rule, the arc coordinate, and the refusals.

Host geometry only, so these run anywhere. The properties under test are the ones a later force law
will depend on without re-checking: that the topology stays inside its own claim, that the arc
coordinate tiles the contour without a gap, and that the length the physics sees is the REALISED one
rather than the requested one.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.world.arena import Kind, WorldArena
from aleph.world.strand import Strand, assert_inside_arena, build_strand
from aleph.world import strand as strand_mod


@pytest.fixture()
def arena() -> WorldArena:
    return WorldArena(capacity={Kind.NODE: 5_000, Kind.SEGMENT: 5_000, Kind.ANGLE3: 5_000})


def test_self_check_passes() -> None:
    """The module's own ``_demo`` is the primary check."""
    strand_mod._demo()


def test_composition_rule_is_n_nodes_n_minus_1_segments_n_minus_2_angles(arena: WorldArena) -> None:
    """A strand is nodes + segments + angles in that exact relation — the first composition."""
    s = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=3.0, seg_um=0.05)
    assert s.n_nodes == 61
    assert s.segments.count == s.n_nodes - 1
    assert s.angles.count == s.n_nodes - 2


def test_geometry_matches_the_requested_contour(arena: WorldArena) -> None:
    """The built filament is as long as it was asked to be, end to end."""
    s = build_strand(arena, "cortex", start=(1.0, 2.0, 3.0), direction=(0, 0, 2.0),
                     contour_um=2.0, seg_um=0.1)
    assert np.allclose(s.position[0], (1.0, 2.0, 3.0))
    assert np.isclose(np.linalg.norm(s.position[-1] - s.position[0]), 2.0)
    assert np.allclose(s.seg_rest_um, 0.1), "direction magnitude must not scale the geometry"


def test_realised_segment_length_is_reported_when_rounding_moves_it(arena: WorldArena) -> None:
    """A contour that does not divide by the requested length changes it. The realised value is what a
    downstream law must read: alpha = kappa/seg^3 turns a 1% length error into 3%."""
    s = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=3.0, seg_um=0.055)
    assert s.segments.count == 55
    assert s.seg_um_realised != s.seg_um_requested
    assert s.seg_um_error == pytest.approx(abs(3.0 / 55 - 0.055) / 0.055)
    assert np.allclose(s.seg_rest_um, s.seg_um_realised), "rest lengths follow the realised value"


def test_exact_division_leaves_the_length_untouched(arena: WorldArena) -> None:
    """When it does divide, nothing moves and the error is exactly zero — not merely small."""
    s = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=3.0, seg_um=0.05)
    assert s.seg_um_error == 0.0
    assert s.contour_um_realised == pytest.approx(s.contour_um_requested)


def test_arc_coordinate_tiles_the_contour_without_a_gap(arena: WorldArena) -> None:
    """The material arc coordinate is what a motor binds at, so a gap in it would be a stretch of
    filament nothing can bind to."""
    s = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=1.0, seg_um=0.03)
    assert np.isclose(s.seg_arc_um[0, 0], 0.0)
    assert np.isclose(s.seg_arc_um[-1, 1], s.contour_um_realised)
    assert np.allclose(s.seg_arc_um[:-1, 1], s.seg_arc_um[1:, 0]), "no gap between consecutive segments"
    assert np.all(np.diff(s.seg_arc_um, axis=1) > 0.0), "every span is positive"


def test_topology_indexes_are_global_and_inside_the_claim(arena: WorldArena) -> None:
    """A strand must not reach outside its own node range — sharing a node across populations is a
    weld, not a connection."""
    build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0), contour_um=1.0, seg_um=0.1)
    second = build_strand(arena, "sf_arc", start=(0, 5, 0), direction=(1, 0, 0),
                          contour_um=1.0, seg_um=0.1)
    assert second.nodes.lo > 0, "the second strand starts after the first"
    assert int(second.seg_node.min()) >= second.nodes.lo
    assert int(second.angle_idx.max()) < second.nodes.hi
    assert_inside_arena(second)


def test_out_of_range_topology_is_caught(arena: WorldArena) -> None:
    """The assertion must fire if a builder ever writes outside its claim."""
    s = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=1.0, seg_um=0.1)
    bad = Strand(**{**{f: getattr(s, f) for f in s.__slots__},
                    "seg_node": s.seg_node + s.nodes.hi})
    with pytest.raises(AssertionError, match="weld"):
        assert_inside_arena(bad)


@pytest.mark.parametrize("polarity,expect_last", [(+1, True), (-1, False)])
def test_tip_is_derived_from_polarity_not_passed_separately(
    arena: WorldArena, polarity: int, expect_last: bool
) -> None:
    """Deriving the tip means the two cannot disagree."""
    s = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=1.0, seg_um=0.1, polarity=polarity)
    assert s.tip_node == (s.nodes.hi - 1 if expect_last else s.nodes.lo)


def test_polarity_must_be_plus_or_minus_one(arena: WorldArena) -> None:
    with pytest.raises(ValueError, match="polarity"):
        build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=1.0, seg_um=0.1, polarity=0)


@pytest.mark.parametrize("kwargs", [{"contour_um": 3.0}, {"seg_um": 0.05}, {}])
def test_lengths_refuse_to_default(arena: WorldArena, kwargs: dict) -> None:
    """The ratified pattern applied at its first opportunity: a default segment length would recreate
    the seg_um = 0.5 the provenance audit ranks #2 load-bearing with a source of "none"."""
    with pytest.raises(ValueError, match="no default"):
        build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0), **kwargs)


@pytest.mark.parametrize("bad", [0.0, -1.0, float("inf"), float("nan")])
def test_non_positive_or_non_finite_length_is_refused(arena: WorldArena, bad: float) -> None:
    with pytest.raises(ValueError, match="finite and positive"):
        build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=1.0, seg_um=bad)


def test_too_coarse_to_bend_is_refused(arena: WorldArena) -> None:
    """Two nodes carry no bending triple, so the object has no bending stiffness at all — a different
    physical object, not a coarser one."""
    with pytest.raises(ValueError, match="at least 3"):
        build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0),
                     contour_um=1.0, seg_um=1.0)


def test_degenerate_direction_is_refused(arena: WorldArena) -> None:
    with pytest.raises(ValueError, match="degenerate"):
        build_strand(arena, "cortex", start=(0, 0, 0), direction=(0, 0, 0),
                     contour_um=1.0, seg_um=0.1)


def test_claims_stay_contiguous_across_many_strands(arena: WorldArena) -> None:
    """The live prefix is what makes launching at dim=n_live correct, so it must survive real use."""
    for i in range(20):
        build_strand(arena, "cortex", start=(0, float(i), 0), direction=(1, 0, 0),
                     contour_um=1.0, seg_um=0.1)
    arena.assert_partitioned()
    c = arena.census()
    assert c["populations"]["cortex"]["node"] == 20 * 11
    assert c["populations"]["cortex"]["segment"] == 20 * 10
    assert c["populations"]["cortex"]["angle3"] == 20 * 9
    assert arena.n_live(Kind.NODE) == 20 * 11
