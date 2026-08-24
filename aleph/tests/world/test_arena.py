"""Arena invariants — the three properties the layer exists to hold.

These run without a device on purpose: the arena's bookkeeping is the part that must be correct before
any allocation happens, and a test that needs a card cannot run in the place most of these edits are
made. The device path gets one test of its own, which asserts the REFUSAL rather than the allocation —
that a non-CUDA device raises is a contract, and it is checkable anywhere.
"""

from __future__ import annotations

import pytest

from aleph.world.arena import Claim, Kind, WorldArena
from aleph.world import arena as arena_mod


def test_self_check_passes() -> None:
    """The module's own ``_demo`` is the primary check; a failure there is a failure here."""
    arena_mod._demo()


def test_claims_tile_the_live_prefix_exactly() -> None:
    """Live is the claim total, and the claims tile ``[0, n_live)`` with no gap — that is what makes
    launching at ``dim=n_live`` against the full array correct without slicing."""
    a = WorldArena(capacity={Kind.NODE: 1_000})
    first = a.claim("cortex", Kind.NODE, 600)
    second = a.claim("membrane", Kind.NODE, 300)

    assert first.hi == second.lo, "the second claim must begin where the first ends"
    assert a.n_live(Kind.NODE) == first.count + second.count
    assert a.n_live(Kind.NODE) < a.capacity[Kind.NODE], "live must not be capacity"
    a.assert_partitioned()


def test_unclaimed_kind_is_live_zero_not_live_capacity() -> None:
    """A reserved-but-unclaimed kind contributes nothing to any launch dimension."""
    a = WorldArena(capacity={Kind.NODE: 100, Kind.BOND: 50})
    a.claim("cortex", Kind.NODE, 10)
    assert a.n_live(Kind.BOND) == 0
    assert a.capacity[Kind.BOND] == 50


def test_population_is_derived_from_the_range_not_declared() -> None:
    """A bond stores two node ids; its component pair is read off the ranges. Beyond the live prefix
    there is no owner — an unclaimed index must not inherit the last claim."""
    a = WorldArena(capacity={Kind.NODE: 100})
    a.claim("cortex", Kind.NODE, 60)
    a.claim("membrane", Kind.NODE, 30)

    assert a.population_of(Kind.NODE, 0) == "cortex"
    assert a.population_of(Kind.NODE, 59) == "cortex"
    assert a.population_of(Kind.NODE, 60) == "membrane"
    assert a.population_of(Kind.NODE, 95) is None


@pytest.mark.parametrize("count", [0, -1])
def test_empty_claim_is_refused(count: int) -> None:
    """"Present but silent" and "failed to build" are different facts, and an empty range erases the
    difference downstream."""
    a = WorldArena(capacity={Kind.BOND: 10})
    with pytest.raises(ValueError, match="must be positive"):
        a.claim("ecm", Kind.BOND, count)


def test_over_capacity_claim_names_both_numbers() -> None:
    """The remedy is to raise the reservation, never to lower a population: the count comes from
    physiological density x geometry, so the error must not read as an invitation to shrink it."""
    a = WorldArena(capacity={Kind.NODE: 100})
    a.claim("cortex", Kind.NODE, 90)
    with pytest.raises(ValueError) as excinfo:
        a.claim("membrane", Kind.NODE, 20)
    message = str(excinfo.value)
    assert "exceeds capacity" in message
    assert "90" in message and "20" in message and "100" in message


def test_unnamed_claim_is_refused() -> None:
    """A range with no population cannot be attributed in a census or a ledger cut."""
    a = WorldArena(capacity={Kind.NODE: 10})
    with pytest.raises(ValueError, match="population name"):
        a.claim("", Kind.NODE, 1)


def test_partition_check_catches_a_gap() -> None:
    """The assertion is true by construction while claiming is contiguous — which is why it must fail
    loudly the moment a code path claims out of order."""
    a = WorldArena(capacity={Kind.NODE: 100})
    a._claims.append(Claim("ghost", Kind.NODE, 10, 20))
    a._live[Kind.NODE] = 20
    with pytest.raises(AssertionError, match="not contiguous"):
        a.assert_partitioned()


def test_partition_check_catches_a_live_count_disagreement() -> None:
    """A live count that outruns its claims would make every launch read past the last real element."""
    a = WorldArena(capacity={Kind.NODE: 100})
    a.claim("cortex", Kind.NODE, 10)
    a._live[Kind.NODE] = 40
    with pytest.raises(AssertionError, match="disagrees with the claims"):
        a.assert_partitioned()


def test_claims_of_different_kinds_never_overlap() -> None:
    """Ranges are per-kind, so identical index bounds in NODE and SEGMENT are not a collision."""
    a = WorldArena(capacity={Kind.NODE: 10, Kind.SEGMENT: 10})
    node = a.claim("cortex", Kind.NODE, 10)
    segment = a.claim("cortex", Kind.SEGMENT, 10)
    assert (node.lo, node.hi) == (segment.lo, segment.hi)
    assert not node.overlaps(segment), "same bounds, different kind — not an overlap"


def test_census_reports_live_against_capacity() -> None:
    """The floorplan view's payload: host-side, touching no device array, safe every accepted step."""
    a = WorldArena(capacity={Kind.NODE: 1_000, Kind.BOND: 100})
    a.claim("cortex", Kind.NODE, 600)
    a.claim("cortex", Kind.BOND, 40)
    a.claim("membrane", Kind.NODE, 300)

    c = a.census()
    assert c["live"]["node"] == 900 and c["capacity"]["node"] == 1_000
    assert c["utilisation"]["node"] == pytest.approx(0.9)
    assert c["populations"]["cortex"] == {"node": 600, "bond": 40}
    assert c["populations"]["membrane"] == {"node": 300}
    assert c["n_claims"] == 3
    assert c["device"] is None and c["node_arrays"] == []


def test_non_cuda_device_raises_rather_than_falling_back() -> None:
    """There is no CPU simulation path and there must never be one, so resolving to a non-CUDA device
    is an error and not a quieter mode."""
    with pytest.raises(RuntimeError, match="not CUDA"):
        WorldArena(capacity={Kind.NODE: 8}, device="cpu")


def test_byte_span_refuses_a_non_node_claim() -> None:
    """Only node arrays are spanned here; asking for a bond's span is a caller error, not an empty
    answer."""
    a = WorldArena(capacity={Kind.BOND: 10})
    bond = a.claim("ecm", Kind.BOND, 4)
    with pytest.raises(ValueError, match="for NODE claims"):
        a.byte_span("position", bond)
