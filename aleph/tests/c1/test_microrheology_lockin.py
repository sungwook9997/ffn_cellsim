"""Sanity gates for the lock-in detection and the phase crossing — against signals with known answers.

The oscillatory method's whole result is a phase, and a phase is the easiest thing in this codebase to
get wrong by a sign or a quarter turn while still producing a smooth, plausible, monotone curve. So the
detector is checked on synthetic signals whose amplitude and lag are known exactly, including the
multi-tone case the driver actually uses (orthogonality of distinct tones over whole periods is the
only reason one run can carry a whole spectrum).

Host arithmetic; no Warp kernels, no field, no physics.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aleph.scripts.c1_biot_microrheology import _lock_in, _phase_crossing


def _grid(t_end: float, n: int = 20000) -> np.ndarray:
    return np.linspace(0.0, t_end, n, endpoint=False)


def test_a_pure_tone_returns_its_amplitude_and_zero_lag() -> None:
    w = 3.7
    t = _grid(8 * 2 * math.pi / w)
    amp, ph = _lock_in(t, 2.5 * np.sin(w * t), w)
    assert amp == pytest.approx(2.5, rel=1e-3)
    assert ph == pytest.approx(0.0, abs=1e-3)


@pytest.mark.parametrize("lag_deg", [15.0, 45.0, 90.0, 135.0])
def test_a_lagged_tone_returns_its_lag(lag_deg: float) -> None:
    """sin(wt - phi) must read back as a POSITIVE lag phi — a driven medium lags, never leads."""
    w, phi = 2.1, math.radians(lag_deg)
    t = _grid(12 * 2 * math.pi / w)
    amp, ph = _lock_in(t, 1.3 * np.sin(w * t - phi), w)
    assert amp == pytest.approx(1.3, rel=1e-3)
    assert math.degrees(ph) == pytest.approx(lag_deg, abs=0.2)


def test_distinct_tones_do_not_leak_into_each_other() -> None:
    """The multi-tone drive only works if projections at distinct tones are orthogonal."""
    w1, w2 = 1.0, 4.0                       # whole periods of both fit in the window
    t = _grid(8 * 2 * math.pi / w1)
    y = 1.0 * np.sin(w1 * t - 0.3) + 0.4 * np.sin(w2 * t - 1.1)
    a1, p1 = _lock_in(t, y, w1)
    a2, p2 = _lock_in(t, y, w2)
    assert a1 == pytest.approx(1.0, rel=2e-3) and math.degrees(p1) == pytest.approx(17.19, abs=0.3)
    assert a2 == pytest.approx(0.4, rel=2e-3) and math.degrees(p2) == pytest.approx(63.03, abs=0.3)


def test_a_first_order_response_crosses_45_degrees_at_its_corner() -> None:
    """The declared inversion: for a first-order low-pass, phase = 45 deg exactly at omega = 1/tau.

    This is the arithmetic the driver performs on the measured phases, checked against the model whose
    corner it is defined by — so a wrong interpolation cannot hide behind a plausible-looking sweep.
    """
    tau = 0.037
    w = np.array([0.25, 0.5, 1.0, 2.0, 4.0]) / tau
    phases = np.arctan(w * tau)
    w45 = _phase_crossing(w, phases, math.pi / 4.0)
    assert w45 == pytest.approx(1.0 / tau, rel=2e-2)


def test_no_crossing_inside_the_band_returns_none() -> None:
    """A band that never reaches 45 degrees must not report its endpoint as a corner."""
    w = np.array([1.0, 2.0, 3.0])
    assert _phase_crossing(w, np.radians(np.array([1.0, 2.0, 3.0])), math.pi / 4.0) is None


def test_the_crossing_is_bracketed_by_the_declared_tone_spread() -> None:
    """The tone multiples must straddle the corner, or every run reports None by construction."""
    from aleph.scripts.c1_biot_microrheology import TONE_MULTIPLES
    phases = np.arctan(np.array(TONE_MULTIPLES))          # corner at multiple 1.0
    assert phases[0] < math.pi / 4.0 < phases[-1]


# ── the amended observable, and the reason it was amended ────────────────────────────────────────────
#: Run 34's measured phase profiles, D=50, dx=0.2, tones at (0.25 .. 4.0) x D/a^2.
RUN34_PHASES_DEG = {
    0.8: [23.2, 30.0, 37.3, 43.1, 43.6],
    1.1: [24.5, 32.8, 42.8, 54.1, 65.5],
    1.5: [23.3, 30.1, 37.4, 43.1, 43.5],
    2.0: [23.0, 29.8, 37.0, 42.7, 42.7],
}


def test_a_phase_profile_that_stops_below_45_yields_no_crossing() -> None:
    """Three of run 34's four points never reach 45 degrees; the detector must say so.

    Reporting the band edge as a corner would be the failure mode: a number where there is none.
    """
    from aleph.scripts.c1_biot_microrheology import _phase_crossing
    w = np.array([0.25, 0.5, 1.0, 2.0, 4.0])
    for a in (0.8, 1.5, 2.0):
        assert _phase_crossing(w, np.radians(RUN34_PHASES_DEG[a]), math.pi / 4.0) is None


def test_the_self_similarity_guard_rejects_the_one_point_that_reported_a_number() -> None:
    """Run 34's central finding, pinned: the only point whose inversion fired is the corrupted one.

    Tones are multiples of D/a^2, so the dimensionless response cannot depend on `a` — three profiles
    agree to three figures and one does not. The guard uses the MEDIAN profile as reference so a
    single corrupted point cannot drag the reference onto itself.
    """
    from aleph.scripts.c1_biot_microrheology import SELF_SIMILARITY_MAX_DEG
    prof = np.array([RUN34_PHASES_DEG[a] for a in (0.8, 1.1, 1.5, 2.0)], dtype=float)
    dev = np.max(np.abs(prof - np.median(prof, axis=0)), axis=1)
    assert dev[0] <= SELF_SIMILARITY_MAX_DEG
    assert dev[2] <= SELF_SIMILARITY_MAX_DEG
    assert dev[3] <= SELF_SIMILARITY_MAX_DEG
    assert dev[1] > SELF_SIMILARITY_MAX_DEG, "the a=1.1 outlier must be rejected"
    assert dev[1] == pytest.approx(21.95, abs=0.01)


def test_the_amplitude_crossing_finds_the_minus_3dB_point() -> None:
    """The replacement observable, against a roll-off whose -3 dB point is known exactly."""
    from aleph.scripts.c1_biot_microrheology import AMPLITUDE_LEVEL, _amplitude_crossing
    tau = 0.02
    w = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0]) / tau
    amp = 1.0 / np.sqrt(1.0 + (w * tau) ** 2)             # -3 dB at omega = 1/tau by construction
    w3 = _amplitude_crossing(w, amp / amp[0], AMPLITUDE_LEVEL)
    assert w3 == pytest.approx(1.0 / tau, rel=5e-2)


def test_a_flat_response_has_no_minus_3dB_point() -> None:
    """A band that never rolls off must return None rather than its last frequency."""
    from aleph.scripts.c1_biot_microrheology import AMPLITUDE_LEVEL, _amplitude_crossing
    w = np.array([1.0, 2.0, 4.0])
    assert _amplitude_crossing(w, np.array([1.0, 1.0, 1.0]), AMPLITUDE_LEVEL) is None
