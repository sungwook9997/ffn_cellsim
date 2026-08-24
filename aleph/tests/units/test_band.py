"""Controls for `ALEPH-PORT-4001` — cited acceptance bands.

Every test here is named in the ledger entry's §8 or §9. `tests/ports/test_port_discipline.py`
checks that correspondence in both directions: an entry may not name a control nobody wrote, and
this file may not quietly rename one the entry points at.

Nothing in here imports from `Project_Aleph`. The oracle is the definition of set membership, per
the ledger's §7 — if this port agreed with the source and disagreed with the definition, the
definition wins.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aleph.units import LiteratureBand, OutOfBandError, band_guard

CITED = LiteratureBand(
    lo=1e-2, hi=1.0, unit="pN/um", citation="ALEPH-PORT-4001 test fixture -- not a physical claim"
)


# ── §8 positive controls ──────────────────────────────────────────────────────────────────────
def test_membership_matches_the_definition() -> None:
    """`contains` agrees with `lo <= x <= hi` at every one of 2,001 points across both bounds."""
    grid = np.linspace(-1.0, 2.0, 2001)
    for x in grid:
        expected = CITED.lo <= float(x) <= CITED.hi
        assert CITED.contains(float(x)) is expected, x
        if expected:
            assert CITED.guard("x", float(x)) == float(x)
        else:
            with pytest.raises(OutOfBandError):
                CITED.guard("x", float(x))


def test_guard_returns_the_value_untouched() -> None:
    """I1 — the scalar guard returns an equal float, the array guard the *same object*."""
    assert CITED.guard("sigma", 0.5) == 0.5

    values = np.array([0.02, 0.5, 0.99])
    returned = CITED.guard_all("sigma", values)
    assert returned is values, "guard_all must not hand back a different array"
    assert np.array_equal(returned, np.array([0.02, 0.5, 0.99]))


# ── §9 negative controls, each of which must be able to fail ──────────────────────────────────
def test_out_of_band_raises_and_does_not_clamp() -> None:
    """I2 — one ULP above `hi` is refused, and the message carries the citation."""
    just_over = math.nextafter(CITED.hi, math.inf)
    assert just_over > CITED.hi

    with pytest.raises(OutOfBandError) as excinfo:
        CITED.guard("sigma_cortex", just_over)

    message = str(excinfo.value)
    assert "ALEPH-PORT-4001 test fixture" in message, "the citation must reach the reader"
    assert "pN/um" in message
    # If the guard is ever changed to clamp, it returns hi instead of raising and this test is the
    # thing that notices.
    assert "sigma_cortex" in message


def test_a_band_without_a_citation_is_refused() -> None:
    """I3 — refused at construction, not at first use."""
    for empty in ("", "   ", "\n\t"):
        with pytest.raises(ValueError, match="citation"):
            LiteratureBand(lo=0.0, hi=1.0, unit="pN", citation=empty)


def test_nan_is_outside_even_an_infinite_band() -> None:
    """I4 — the case the two obvious spellings of the interval test disagree about.

    `lo <= x <= hi` says False for NaN; `not (x < lo or x > hi)` says True. A one-sided band is the
    hardest version, because everything else about it admits almost everything.
    """
    one_sided = LiteratureBand(
        lo=0.0, hi=math.inf, unit="pN", citation="strictly positive, by construction"
    )
    assert one_sided.contains(1e300) is True
    assert one_sided.contains(float("nan")) is False
    assert one_sided.contains(float("-inf")) is False
    with pytest.raises(OutOfBandError):
        one_sided.guard("f", float("nan"))

    with pytest.raises(OutOfBandError):
        CITED.guard_all("sigma", np.array([0.5, np.nan, 0.5]))


def test_a_band_infinite_on_both_sides_is_refused() -> None:
    """A band that admits every finite value is not a check, and saying so is cheaper than a note."""
    with pytest.raises(ValueError, match="infinite on both sides"):
        LiteratureBand(lo=-math.inf, hi=math.inf, unit="pN", citation="anything")
    with pytest.raises(ValueError, match="NaN"):
        LiteratureBand(lo=float("nan"), hi=1.0, unit="pN", citation="anything")


def test_malformed_and_unitless_bands_are_refused() -> None:
    with pytest.raises(ValueError, match="malformed"):
        LiteratureBand(lo=1.0, hi=0.0, unit="pN", citation="anything")
    with pytest.raises(ValueError, match="unit"):
        LiteratureBand(lo=0.0, hi=1.0, unit=" ", citation="anything")


# ── vacuity: a control that cannot fail is decoration ─────────────────────────────────────────
def test_the_positive_control_would_notice_a_broken_band() -> None:
    """Widen the band and the definition comparison must disagree — so it is load-bearing."""
    widened = LiteratureBand(
        lo=CITED.lo, hi=CITED.hi * 2.0, unit=CITED.unit, citation=CITED.citation
    )
    disagreements = [
        x for x in np.linspace(-1.0, 2.0, 2001) if widened.contains(float(x)) is not (CITED.lo <= float(x) <= CITED.hi)
    ]
    assert disagreements, (
        "the membership comparison in test_membership_matches_the_definition would pass against a "
        "wrong band, which would make it decoration"
    )


def test_exclusive_bounds_exclude_their_endpoints() -> None:
    """`inclusive=False` is unexercised in practice (ledger §14.5) but is not untested."""
    strict = LiteratureBand(lo=0.0, hi=1.0, unit="1", citation="strict bound", inclusive=False)
    assert strict.contains(0.0) is False
    assert strict.contains(1.0) is False
    assert strict.contains(0.5) is True


def test_band_guard_is_the_inline_spelling() -> None:
    assert band_guard("x", 0.5, 0.0, 1.0, "1", "inline") == 0.5
    with pytest.raises(OutOfBandError):
        band_guard("x", 1.5, 0.0, 1.0, "1", "inline")
    with pytest.raises(ValueError, match="citation"):
        band_guard("x", 0.5, 0.0, 1.0, "1", "")


def test_array_guard_reports_how_many_and_which() -> None:
    with pytest.raises(OutOfBandError) as excinfo:
        CITED.guard_all("sigma", np.array([0.5, 5.0, 0.5, -3.0]))
    message = str(excinfo.value)
    assert "2 of 4" in message
    assert "worst" in message


def test_a_wrong_dtype_raises_rather_than_becoming_nan() -> None:
    """Ledger §5: coercion would turn a dtype bug into an out-of-band report and hide which it was."""
    with pytest.raises((TypeError, ValueError)):
        CITED.guard_all("sigma", np.array(["not", "numbers"]))
