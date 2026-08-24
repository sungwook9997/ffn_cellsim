"""Self-test of the walk_dir barbed-end-polarity hand-off oracle (pure NumPy — no Warp/CUDA)."""

from __future__ import annotations

import numpy as np
import pytest

from aleph.components.weave.walk_dir import barbed_end_node, fill_walk_dir, walk_dir_from_polarity


def test_barbed_end_node_from_polarity() -> None:
    """+1 polarity -> barbed end at the LAST node; -1 -> the FIRST node."""
    off = np.array([0, 5, 9], np.int64)                # two fibers: [0..4], [5..8]
    pol = np.array([1, -1], np.int64)
    barbed = barbed_end_node(off, pol)
    assert barbed[0] == 4                               # last node of fiber 0
    assert barbed[1] == 5                               # first node of fiber 1


def test_walk_dir_points_to_barbed_end() -> None:
    """A bound head's walk_dir is the unit vector from its anchor toward the filament barbed end."""
    pos = np.array([[0.0, 0, 0], [1, 0, 0], [2, 0, 0]])
    wd = walk_dir_from_polarity(anchor_node=0, barbed_node=2, pos=pos)
    assert np.allclose(wd, [1.0, 0.0, 0.0])
    assert np.linalg.norm(wd) == pytest.approx(1.0)


def test_unbound_head_is_passive() -> None:
    """An unbound head (anchor -1) or a head at the barbed end gets zero walk_dir (passive)."""
    pos = np.array([[0.0, 0, 0], [1, 0, 0]])
    assert np.allclose(walk_dir_from_polarity(-1, 1, pos), 0.0)
    assert np.allclose(walk_dir_from_polarity(1, 1, pos), 0.0)   # anchor == barbed -> undefined -> passive


def test_reattach_flips_walk_dir_on_opposite_polarity() -> None:
    """Rebinding a head to a filament of OPPOSITE polarity flips its walk_dir (tracks the actual actin)."""
    # one straight filament along +x; barbed at last node (pol +1) then at first node (pol -1)
    pos = np.array([[0.0, 0, 0], [1, 0, 0], [2, 0, 0]])
    off = np.array([0, 3], np.int64)
    node_fiber = np.array([0, 0, 0], np.int64)
    barbed_plus = barbed_end_node(off, np.array([1], np.int64))    # node 2
    barbed_minus = barbed_end_node(off, np.array([-1], np.int64))  # node 0
    bound = np.array([True]); anchor = np.array([1], np.int64)
    wd_plus = fill_walk_dir(bound, anchor, node_fiber, barbed_plus, pos)[0]
    wd_minus = fill_walk_dir(bound, anchor, node_fiber, barbed_minus, pos)[0]
    assert np.allclose(wd_plus, [1, 0, 0])
    assert np.allclose(wd_minus, [-1, 0, 0])            # flipped


def test_bipolar_contraction_sign() -> None:
    """Two anti-parallel actin filaments with barbed ends OUTWARD -> the two heads walk outward; the reaction
    slides the filaments so their pointed ends converge (net contraction), matching the powerstroke sign."""
    # filament A: nodes at x in [-2,-1,0], barbed at x=-2 (outward, left); filament B: x in [0,1,2], barbed at x=2
    pos = np.array([[-2.0, 0, 0], [-1, 0, 0], [0, 0, 0],   # A (0,1,2), barbed = node 0 (x=-2)
                    [0.0, 0, 0], [1, 0, 0], [2, 0, 0]])    # B (3,4,5), barbed = node 5 (x=+2)
    node_fiber = np.array([0, 0, 0, 1, 1, 1], np.int64)
    barbed = np.array([0, 5], np.int64)                 # A barbed outward (left), B barbed outward (right)
    # one head on each filament, near the centre
    bound = np.array([True, True]); anchor = np.array([1, 4], np.int64)
    wd = fill_walk_dir(bound, anchor, node_fiber, barbed, pos)
    assert wd[0][0] < 0                                 # head on A walks LEFT (toward A's outward barbed)
    assert wd[1][0] > 0                                 # head on B walks RIGHT (toward B's outward barbed)
    # heads pull actin the OPPOSITE way (reaction) -> A moves right, B moves left -> pointed ends converge
    # => net contraction of the minifilament straddling them (sign check).
    assert wd[0][0] * wd[1][0] < 0                      # anti-parallel walk directions (bipolar)
