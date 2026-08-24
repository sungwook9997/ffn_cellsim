"""The equilibration search can keep a window that still opens on a transient — and then it self-reinforces.

WHAT THIS CAUGHT.  `equilibration_point` maximises ``n_eff = T / (2·tau_int)``, on the assumption its own
docstring states: "discarding too little leaves the transient in and INFLATES the correlation time".  On
2026-07-29 that assumption failed on a real full-native run.  A γ series that rose 0.8 → 3.4 in 0.5 s and
was then flat for 7 s got ``start_index = 0``: ``T`` is maximal there, and the ACF estimator did not
penalise the trend enough to compensate.  The "steady-state mean" was therefore taken ACROSS THE RISE —
3.373 instead of the plateau's 3.415 — which is the single failure this module exists to prevent.

AND IT HID ITSELF.  The drift verdict is ``|drift| < drift_sigma · std``.  Keeping the transient inflated
``std`` twentyfold (0.0106 → 0.2538), so a BIGGER transient makes the drift test EASIER to pass.  The run
was reported STATIONARY at 0.94 σ.  A guard normalised by the same inflated ``std`` would inherit that
loop, which is why this one is normalised by the trailing tenth instead — transient-free by construction.

WHY THE STATISTIC IS A SEGMENT AND NOT A POINT.  The first draft compared the first SAMPLE against the
window's scatter.  It rejected all three real runs, because a cut placed at the knee necessarily leaves its
first sample out in the tail.  A transient is a property of a segment.  Leading tenth vs trailing tenth
reads 0.0 / 0.8 / 33.5 on the three measured runs — the two good cuts keep, the contaminated one rejects.

The threshold is ONE standard deviation: a definition, not a number chosen after seeing the data.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.engine.observe.stationarity import (
    StationarityVerdict,
    assess_stationarity,
    window_opens_on_a_transient,
)

DRIVING_TAU_S = 1.4558568761100406


def _rise_then_plateau(n: int = 152, rise_frac: float = 0.06, seed: int = 0) -> np.ndarray:
    """A series shaped like the run that broke the detector: sharp rise, then flat with small noise."""
    rng = np.random.default_rng(seed)
    k = max(2, int(n * rise_frac))
    rise = np.linspace(0.8, 3.41, k)
    flat = 3.41 + rng.normal(0.0, 0.011, n - k)
    return np.concatenate([rise, flat])


def _plateau_only(n: int = 152, seed: int = 1) -> np.ndarray:
    """An already-equilibrated series with the same scatter and no transient."""
    return 3.41 + np.random.default_rng(seed).normal(0.0, 0.011, n)


class TestGuard:
    def test_a_window_opening_on_a_rise_is_flagged(self) -> None:
        """The failure case: keeping index 0 of a rise-then-plateau series."""
        contaminated, ratio = window_opens_on_a_transient(_rise_then_plateau(), 0)
        assert contaminated and ratio > 1.0

    def test_an_equilibrated_window_is_not_flagged(self) -> None:
        """The guard must not fire on the thing it is meant to permit."""
        contaminated, ratio = window_opens_on_a_transient(_plateau_only(), 0)
        assert not contaminated and ratio <= 1.0

    def test_discarding_the_rise_clears_the_flag(self) -> None:
        """Cutting past the transient makes the same series admissible — the guard tracks the CUT."""
        series = _rise_then_plateau()
        k = max(2, int(series.size * 0.06))
        assert window_opens_on_a_transient(series, k)[0] is False

    def test_it_is_normalised_by_the_trailing_tenth_not_the_whole_window(self) -> None:
        """The whole window's std is inflated BY the transient; normalising by it would inherit the loop.

        Scaling the plateau noise down makes the transient far more significant.  A guard normalised by the
        full-window std would barely move, because that std is dominated by the rise either way.
        """
        mild = _rise_then_plateau(seed=2)
        sharp = mild.copy()
        k = max(2, int(mild.size * 0.06))
        sharp[k:] = 3.41 + (sharp[k:] - 3.41) * 0.1        # same rise, one-tenth the plateau scatter
        assert window_opens_on_a_transient(sharp, 0)[1] > window_opens_on_a_transient(mild, 0)[1] * 5

    def test_a_series_too_short_to_split_is_not_judged(self) -> None:
        """Under 20 samples there are no tenths to compare; the guard must abstain, not guess."""
        assert window_opens_on_a_transient(np.arange(12.0), 0) == (False, 0.0)

    def test_a_constant_series_does_not_divide_by_zero(self) -> None:
        """Zero scatter in the trailing tenth is a boundary, not an exception."""
        assert window_opens_on_a_transient(np.full(60, 2.5), 0) == (False, 0.0)


class TestVerdict:
    def test_a_contaminated_window_returns_DRIFTING_WITH_NO_MEAN(self) -> None:
        """The point of the guard: nothing downstream may quote the mean taken across a rise."""
        report = assess_stationarity(
            _rise_then_plateau(), 0.05, observable="gamma_total_pn_per_um",
            driving_correlation_time_s=DRIVING_TAU_S, min_windows=5.0, drift_sigma=1.0)
        assert report.verdict is StationarityVerdict.DRIFTING
        assert report.mean is None and report.sem is None

    def test_the_reason_names_the_self_reinforcing_loop(self) -> None:
        """A reader who sees only 'DRIFTING' will re-run it longer; they need to know std inflates too."""
        report = assess_stationarity(
            _rise_then_plateau(), 0.05, observable="g",
            driving_correlation_time_s=DRIVING_TAU_S, min_windows=5.0, drift_sigma=1.0)
        reason = report.notes.get("reason", "")
        assert "OPENS on a transient" in reason and "EASIER to pass" in reason

    def test_an_equilibrated_series_still_gets_its_mean(self) -> None:
        """Guarding must not cost the verdict it was added to protect."""
        report = assess_stationarity(
            _plateau_only(n=400), 0.05, observable="g",
            driving_correlation_time_s=DRIVING_TAU_S, min_windows=5.0, drift_sigma=1.0)
        assert report.verdict is StationarityVerdict.STATIONARY
        assert report.mean == pytest.approx(3.41, abs=0.01)
