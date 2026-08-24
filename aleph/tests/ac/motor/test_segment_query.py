r"""Host gate for the point-to-segment binding query geometry (P0#2 primitive).

Validates the NumPy reference (:func:`point_to_segment_query_reference`) that the Warp device kernels mirror
bit-for-formula. The CUDA kernels are NOT launched here (dev-Mac I0-A; the foundation contract forbids a
Warp-CPU launch in ac tests); the on-device parity is a native gate on the A5000.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.motor.segment_query import closest_point_on_segment, point_to_segment_query_reference


def test_closest_point_on_segment_endpoints_and_interior() -> None:
    """t clamps to [0,1]; interior projection is exact; distance is the perpendicular for an interior foot."""
    a, b = np.array([0.0, 0, 0]), np.array([1.0, 0, 0])
    t, c, d = closest_point_on_segment(np.array([0.5, 0.3, 0]), a, b)
    assert t == pytest.approx(0.5) and c == pytest.approx([0.5, 0, 0]) and d == pytest.approx(0.3)
    # beyond endpoint b → clamps to t=1 (the endpoint), distance to b
    t2, c2, d2 = closest_point_on_segment(np.array([2.0, 0, 0]), a, b)
    assert t2 == pytest.approx(1.0) and c2 == pytest.approx([1.0, 0, 0]) and d2 == pytest.approx(1.0)
    # degenerate zero-length segment → t=0, distance to a
    t3, c3, d3 = closest_point_on_segment(np.array([0.0, 1.0, 0]), a, a)
    assert t3 == 0.0 and d3 == pytest.approx(1.0)


def test_query_5field_contract() -> None:
    """The query returns seg_id / t / attach / barbed / dist; seg_id = −1 outside capture."""
    seg_a = np.array([[0, 0, 0], [1, 0, 0]], float)
    seg_b = np.array([[1, 0, 0], [1, 1, 0]], float)
    barbed = np.array([[1, 0, 0], [0, 1, 0]], float)          # per-segment barbed unit dir (walk direction)
    heads = np.array([[0.5, 0.05, 0], [1.0, 0.5, 0], [2.0, 2.0, 0]], float)
    r = point_to_segment_query_reference(heads, seg_a, seg_b, barbed, capture_radius=0.1)
    assert list(r["seg_id"]) == [0, 1, -1]                    # near seg0, on seg1, far → unbound
    assert r["t"] == pytest.approx([0.5, 0.5, 1.0])
    assert r["attach"][0] == pytest.approx([0.5, 0, 0])       # the actual crossbridge anchor point
    assert r["barbed"][1] == pytest.approx([0, 1, 0])         # walk direction from the bound segment's polarity
    assert r["dist"][:2] == pytest.approx([0.05, 0.0])


def test_query_binds_nearest_of_several() -> None:
    """A head between two segments binds the closer one (nearest point-to-segment, not nearest midpoint)."""
    # a long segment whose MIDPOINT is far but whose body is close, vs a short one whose midpoint is nearer
    seg_a = np.array([[-5, 0.02, 0], [0, 0.5, 0]], float)
    seg_b = np.array([[5, 0.02, 0], [0.2, 0.5, 0]], float)
    barbed = np.array([[1, 0, 0], [1, 0, 0]], float)
    r = point_to_segment_query_reference(np.array([[0.0, 0.0, 0]]), seg_a, seg_b, barbed, capture_radius=1.0)
    assert r["seg_id"][0] == 0                                # binds the long segment's body (dist 0.02), not seg1 (0.5)
    assert r["dist"][0] == pytest.approx(0.02)
