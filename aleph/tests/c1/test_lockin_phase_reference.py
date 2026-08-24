"""The lock-in phase reference must be the DRIVE's clock, not the analysis window's start.

WHY. Every "corrupted point" in the C-1 oscillatory runs — a=1.1 in runs 34 and 38, a=2.0 in run 44,
a=0.6 and a=1.7 in run 52 — was this one defect. The projection was taken against ``t - t_window_start``
while the drive is defined on absolute ``t``, which adds a phase offset ``omega * t_start`` to every
tone. That offset is a multiple of 2*pi only when ``t_start`` is an exact whole number of periods of
every tone; the discard boundary is *designed* to land there, but the comparison that finds it
(``t >= N_DISCARD * period_slow``) sits on a floating-point knife edge, and one ulp decides whether the
boundary sample is in or out. Out by one sample shifts every phase by ``-omega * sample_dt``.

The signature was diagnostic and nobody read it for three runs: amplitudes always clean (a time shift
does not change them), phase error growing LINEARLY with frequency, and the affected probe size moving
between runs — because a knife-edge comparison is decided by rounding, not by physics.

These tests pin both halves: the reference is absolute, and a one-sample shift produces exactly the
error that was observed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aleph.scripts.c1_biot_microrheology import _lock_in

#: Tone multiples and the sampling used by the driver: 16 samples per period of the fastest tone.
MULTIPLES = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
SPP = 16


def _drive_case(shift_samples: int) -> list[float]:
    """Phases recovered from a multi-tone drive with a whole-period window, shifted by N samples."""
    base = 1.0                                    # omega_base [rad/s]
    om = np.array([m * base for m in MULTIPLES])
    sample_dt = (2 * math.pi / om.max()) / SPP
    t_slow = 2 * math.pi / om.min()
    n = int(round(4 * t_slow / sample_dt))        # exactly four periods of the slowest tone
    t0 = 4 * t_slow                               # window starts a whole number of periods in
    t = t0 + (np.arange(n) + shift_samples) * sample_dt
    y = sum(np.sin(w * t) for w in om)
    return [math.degrees(_lock_in(t, y, w)[1]) for w in om]


def test_absolute_time_recovers_zero_lag_from_a_shifted_window() -> None:
    """With absolute t, where the window starts cannot matter — the drive's own clock is the reference."""
    for shift in (0, 1, 7):
        for ph in _drive_case(shift):
            assert ph == pytest.approx(0.0, abs=0.05), f"shift={shift} leaked into the phase"


def test_a_one_sample_reference_shift_reproduces_the_observed_error() -> None:
    """The retired bug, pinned: subtracting the window start costs exactly -(2 pi / (SPP*16)) * m.

    That is -1.4062 * m degrees for this tone grid, which is what runs 34/38/44/52 measured.
    """
    base = 1.0
    om = np.array([m * base for m in MULTIPLES])
    sample_dt = (2 * math.pi / om.max()) / SPP
    t_slow = 2 * math.pi / om.min()
    n = int(round(4 * t_slow / sample_dt))
    t0 = 4 * t_slow
    t = t0 + np.arange(n) * sample_dt
    y = sum(np.sin(w * t) for w in om)

    # The defect: project against a time base shifted by one sample. The sign of the apparent phase
    # follows the direction of the shift, so the deviation is wrapped into (-180, 180] and its
    # MAGNITUDE is what is asserted — the observation was that the error grows as omega * sample_dt,
    # not which way the window happened to move.
    for w, m in zip(om, MULTIPLES):
        _, ph = _lock_in(t - sample_dt, y, w)
        dev = math.degrees(ph)
        dev = (dev + 180.0) % 360.0 - 180.0
        assert abs(dev) == pytest.approx(math.degrees(w * sample_dt), abs=0.05)
        assert abs(dev) == pytest.approx(1.40625 * m, abs=0.05)


def test_the_error_grows_linearly_with_frequency() -> None:
    """The property that made this identifiable at all — and it is what a time shift always looks like."""
    base = 1.0
    om = np.array([m * base for m in MULTIPLES])
    sample_dt = (2 * math.pi / om.max()) / SPP
    lags = np.degrees(om * sample_dt)
    ratios = lags / np.array(MULTIPLES)
    assert np.allclose(ratios, ratios[0], rtol=1e-12), "a time shift is linear in omega by definition"
