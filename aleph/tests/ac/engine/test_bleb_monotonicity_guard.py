"""A turning point in the swept range destroys R's premise, not just its value.

``R = gamma(k_xb -> 0) / gamma(reference)`` is a LIMIT claim, and every reading of it — "the fraction of
cortical tension that survives a compliant crossbridge" — presumes gamma DESCENDS toward that limit.  If
gamma turns around inside the sweep, the softest point stops being an extreme and dividing by it answers
no question anyone asked.

WHAT THIS CAUGHT.  On 2026-07-29 the sweep read 4.3256 -> 3.4807 -> 2.8839 -> **12.2425** pN/µm at
``k_xb`` 1000 -> 100 -> 10 -> 3.  Every point was STATIONARY, the parity control agreed bit-identically,
and the span requirement was met — so the analyser would have reported **R = 2.83** with a plateau note
and no objection.  Every guard it had asked whether the points had SETTLED; none asked whether the shape
they traced was the shape R assumes.

The threshold is in combined sems rather than a bare inequality because the softest points carry the
largest sems, and a noise-sized rise must not read as a turning point.
"""

from __future__ import annotations

import math

import pytest

from aleph.scripts.ac_bleb_dose_response import (
    MONOTONIC_SIGMA,
    _monotonicity_breaks,
    refuse_unless_admissible,
)


def _pt(k_xb: float, mean: float, sem: float) -> dict:
    return {
        "name": f"kxb_{k_xb:g}", "k_xb": k_xb, "build": "aaaaaaaa", "steps_run": 800,
        "closure_sha256": "", "stationarity": {
            "gamma_total_pn_per_um": {"verdict": "STATIONARY", "mean": mean, "sem": sem,
                                      "drift_over_std": 0.2},
            "bound_fraction": {"verdict": "STATIONARY", "mean": 0.98, "sem": 1e-4},
        },
    }


#: The measured sweep that motivated the guard.
MEASURED = [_pt(1000, 4.3256, 0.00060), _pt(100, 3.4807, 0.00172),
            _pt(10, 2.8839, 0.00312), _pt(3, 12.2425, 0.06743)]


def test_a_descending_sweep_has_no_break() -> None:
    """The three-point sweep that preceded k_xb=3 must stay admissible — no retroactive invalidation."""
    assert _monotonicity_breaks(MEASURED[:3]) == []


def test_the_measured_turning_point_is_found_with_its_ratio() -> None:
    """4.25x at k_xb=3 over k_xb=10 — the number that must reach the reader."""
    [(softer, stiffer, ratio)] = _monotonicity_breaks(MEASURED)
    assert softer["k_xb"] == 3 and stiffer["k_xb"] == 10
    assert ratio == pytest.approx(12.2425 / 2.8839, rel=1e-9)


def test_a_turning_point_BLOCKS_the_ratio_rather_than_annotating_it() -> None:
    """Non-monotonicity is a refusal, not a footnote: the premise is gone, so the number is meaningless."""
    reasons = refuse_unless_admissible(MEASURED, reference=1000.0)
    assert any("NOT MONOTONIC" in r for r in reasons)
    assert any("do not divide across it" in r for r in reasons)


def test_a_noise_sized_rise_is_not_a_turning_point() -> None:
    """The softest points carry the largest sems; a 1-sigma wobble must not read as a reversal."""
    sem = 0.01
    rise = (MONOTONIC_SIGMA - 0.5) * math.hypot(sem, sem)
    assert _monotonicity_breaks([_pt(100, 3.0, sem), _pt(10, 3.0 + rise, sem)]) == []


def test_a_rise_just_past_the_threshold_is_one() -> None:
    """Pinning the boundary so it cannot drift once someone dislikes a refusal."""
    sem = 0.01
    rise = (MONOTONIC_SIGMA + 0.5) * math.hypot(sem, sem)
    assert len(_monotonicity_breaks([_pt(100, 3.0, sem), _pt(10, 3.0 + rise, sem)])) == 1


def test_a_point_with_no_mean_cannot_manufacture_or_mask_a_break() -> None:
    """An unsettled neighbour is already blocked elsewhere; it must not also fabricate a turning point."""
    pts = [_pt(100, 3.0, 0.01), _pt(10, 99.0, 0.01)]
    pts[1]["stationarity"]["gamma_total_pn_per_um"]["mean"] = None
    assert _monotonicity_breaks(pts) == []


def test_every_adjacent_pair_is_examined_not_just_the_softest_end() -> None:
    """A reversal in the middle of a long sweep is the same defect and must not be skipped."""
    breaks = _monotonicity_breaks(
        [_pt(1000, 4.0, 0.001), _pt(100, 9.0, 0.001), _pt(10, 2.0, 0.001)])
    assert [b[0]["k_xb"] for b in breaks] == [100]
