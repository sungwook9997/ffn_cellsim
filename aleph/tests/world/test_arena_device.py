"""Arena device transfer — CUDA only, and the refusals that are checkable anywhere.

The transfer itself needs a card, so it is gated. What is NOT gated is the set of refusals, because
those are the part that must hold on the machine where the edits are made: a bookkeeping-only arena must
say so rather than fail obscurely, and a claim that is not this arena's must be rejected by range.
"""

from __future__ import annotations

import numpy as np
import pytest
import warp as wp

from aleph.world.arena import Claim, Kind, WorldArena
from aleph.world.strand import build_strand

CUDA = wp.get_cuda_device_count() > 0
cuda_only = pytest.mark.skipif(not CUDA, reason="requires a CUDA device")


def test_bookkeeping_only_arena_says_so() -> None:
    a = WorldArena(capacity={Kind.NODE: 10})
    c = a.claim("cortex", Kind.NODE, 4)
    with pytest.raises(RuntimeError, match="bookkeeping-only mode"):
        a.upload_nodes(c, "position", np.zeros((4, 3)))


def test_non_node_claim_is_refused() -> None:
    a = WorldArena(capacity={Kind.NODE: 10, Kind.BOND: 10})
    a.node_arrays = {"position": object()}          # bypass the device check to reach the kind check
    bond = a.claim("ecm", Kind.BOND, 2)
    with pytest.raises(ValueError, match="needs a NODE claim"):
        a.upload_nodes(bond, "position", np.zeros((2, 3)))


def test_claim_past_the_live_prefix_is_refused() -> None:
    a = WorldArena(capacity={Kind.NODE: 100})
    a.node_arrays = {"position": object()}
    a.claim("cortex", Kind.NODE, 10)
    with pytest.raises(ValueError, match="past the live prefix"):
        a.upload_nodes(Claim("ghost", Kind.NODE, 50, 60), "position", np.zeros((10, 3)))


@cuda_only
def test_upload_writes_only_its_own_range() -> None:
    """The property that makes one array safe to share: a write to one range must leave every other
    range untouched, byte for byte."""
    dev = wp.get_device("cuda:0")
    a = WorldArena(capacity={Kind.NODE: 64}, device=str(dev))
    first = a.claim("cortex", Kind.NODE, 20)
    second = a.claim("membrane", Kind.NODE, 12)

    a.upload_nodes(first, "position", np.full((20, 3), 7.0))
    a.upload_nodes(second, "position", np.full((12, 3), -3.0))

    assert np.allclose(a.download_nodes(first, "position"), 7.0)
    assert np.allclose(a.download_nodes(second, "position"), -3.0)
    rest = a.node_arrays["position"].numpy()[second.hi:]
    assert np.allclose(rest, 0.0), "capacity beyond the live prefix must stay untouched"


@cuda_only
def test_a_strand_round_trips_through_the_device() -> None:
    dev = wp.get_device("cuda:0")
    a = WorldArena(capacity={Kind.NODE: 500, Kind.SEGMENT: 500, Kind.ANGLE3: 500}, device=str(dev))
    s = build_strand(a, "cortex", start=(1.0, 2.0, 3.0), direction=(0, 1, 0),
                     contour_um=1.0, seg_um=0.1)
    a.upload_nodes(s.nodes, "position", s.position)
    assert np.allclose(a.download_nodes(s.nodes, "position"), s.position)


@cuda_only
def test_wrong_length_is_refused_rather_than_partially_written() -> None:
    dev = wp.get_device("cuda:0")
    a = WorldArena(capacity={Kind.NODE: 32}, device=str(dev))
    c = a.claim("cortex", Kind.NODE, 10)
    with pytest.raises(ValueError, match="values for a range"):
        a.upload_nodes(c, "position", np.zeros((4, 3)))
    assert np.allclose(a.download_nodes(c, "position"), 0.0), "nothing was written"


@cuda_only
def test_byte_spans_of_two_claims_do_not_overlap() -> None:
    """Spans, not pointers — a range at index 0 shares the whole array's pointer."""
    dev = wp.get_device("cuda:0")
    a = WorldArena(capacity={Kind.NODE: 64}, device=str(dev))
    first = a.claim("cortex", Kind.NODE, 20)
    second = a.claim("membrane", Kind.NODE, 12)
    lo0, hi0 = a.byte_span("position", first)
    lo1, hi1 = a.byte_span("position", second)
    assert hi0 == lo1 and lo0 < hi0 <= lo1 < hi1
    a.assert_disjoint_spans("position")
