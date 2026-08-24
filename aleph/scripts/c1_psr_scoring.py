r"""C-1 / PSR scoring — relaxation time from a decay curve, and the exponent from a sweep.

The C-1 gate contract (``aleph/docs/v2_audit/C1_PROBE_SCALE_RELAXATION_2026-08-10.md``) scores **an
exponent, not a magnitude**: ``tau(x) ~ x^p`` over a declared sweep.  This module is the only place that
arithmetic lives, so the driver that runs the physics and the driver that runs a synthetic control score
through the same code — a scorer written twice is two scorers.

Host arithmetic only: no Warp, no device state, no physics.  It is therefore testable on a CUDA-free host,
which is the point — the physics is measured on the GPU, the *scoring* is checked against curves whose
answer is known in closed form.

Sanity Gate (self-tested in ``aleph/tests/c1/test_psr_scoring.py``):
    * dimensions: ``t`` [s], ``tau`` [s], the amplitude is dimensionless (normalised to its own t=0).
    * boundary: a curve that never reaches ``1/e`` returns ``None`` rather than the last sample — an
      un-decayed run is not a fast one, and returning its end time would read as exactly that.
    * sign sense: amplitudes must be non-increasing to within ``tol``; a rising curve is a defect
      (energy entering a relaxation), so it raises instead of being fitted.
    * numerical: the ``1/e`` crossing is linearly interpolated in ``log`` amplitude, which is exact for a
      single exponential and second-order otherwise; no tolerance is introduced by the fit itself.
    * measurement protocol: :func:`power_law_exponent` returns ``r2`` beside ``p`` so a caller can apply
      the fit-quality floor the contract requires it to declare BEFORE the run.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

__all__ = [
    "ONE_OVER_E",
    "Q_L2",
    "Q_PEAK",
    "gaussian_spreading_tau",
    "one_over_e_time",
    "power_law_exponent",
]

#: Decay exponent ``q`` of the OBSERVER, in ``A(t) = (a^2/(a^2 + 2Dt))^q``.
#:
#: These two are the whole multi-method argument in miniature.  One diffusing Gaussian blob, read two
#: ways, gives two different relaxation times from the SAME physics:
#:
#:   * an observer that integrates the whole field (an L2 norm) sees ``q = 3/4``, because
#:     ``||p||_2 ~ sigma^(-3/2)`` in three dimensions;
#:   * an observer that reads one point (the peak) sees ``q = 3/2``, because the amplitude of a
#:     mass-conserving Gaussian falls as ``sigma^(-3)``.
#:
#: The ``1/e`` times therefore differ by ``(e^(4/3)-1)/(e^(2/3)-1) = 2.9477x`` while the EXPONENT in ``a``
#: is 2 for both.  That is the C-1 thesis, provable in closed form before any cell is simulated: the
#: scaling transfers across measurement methods, the magnitude does not.
Q_L2: float = 0.75
Q_PEAK: float = 1.5

#: The level the relaxation time is read at.  Declared as a constant rather than written inline at the
#: two call sites, because "the time to decay by 1/e" and "the time to decay by half" are different
#: numbers and a run that mixes them is unfittable in a way nothing downstream would notice.
ONE_OVER_E: float = 1.0 / math.e


def gaussian_spreading_tau(a_um: float, d_um2_s: float, q: float = Q_L2) -> float:
    r"""Closed-form ``1/e`` time of a Gaussian pressure blob, AS SEEN BY AN OBSERVER OF EXPONENT ``q``.

    ACCEPTANCE ORACLE, not a runtime mechanism.  A Gaussian ``p(r,0) = P0 exp(-r^2 / 2a^2)`` diffusing
    under ``dp/dt = D grad^2 p`` stays Gaussian with ``sigma^2(t) = a^2 + 2Dt``.  Any observer whose
    reading falls as a power of the width — ``A(t) = (a^2/(a^2+2Dt))^q`` — therefore has

        ``tau = a^2 (e^(1/q) - 1) / (2 D)``

    so the run can declare an EXPONENT **and a PREFACTOR** in advance, per observer.  With
    :data:`Q_L2` = 3/4 this is ``1.396834 a^2/D``; with :data:`Q_PEAK` = 3/2 it is ``0.473867 a^2/D``.
    **The ratio 2.9477 is the method disagreement, in closed form, with the physics held identical.**

    Valid in free space.  A finite no-flux box only ever makes the measured time LONGER, so the caller
    must keep ``tau`` well under the box's own diffusive time and record that it did.

    Args:
        a_um: Gaussian width (the probe length scale) [µm].
        d_um2_s: Poroelastic diffusivity [µm²/s].
        q: The observer's decay exponent in the blob width. Defaults to the L2-norm observer.

    Returns:
        The ``1/e`` time of that observer's normalised reading [s].
    """
    if a_um <= 0.0 or d_um2_s <= 0.0:
        raise ValueError("a and D must be positive to have a diffusive timescale")
    if q <= 0.0:
        raise ValueError("the observer exponent q must be positive, or its reading does not decay")
    return a_um * a_um * (math.exp(1.0 / q) - 1.0) / (2.0 * d_um2_s)


def one_over_e_time(
    t_s: npt.ArrayLike,
    amplitude: npt.ArrayLike,
    *,
    rise_tol: float = 1e-9,
) -> float | None:
    r"""Return the time at which ``amplitude`` first falls to ``1/e`` of its initial value [s].

    Interpolated linearly in ``log(amplitude)``, which is exact when the decay is a single exponential
    and second-order accurate when it is not.  Deliberately NOT a fitted rate: the poroelastic decay of a
    spreading blob is a power law in ``(a^2 + 2Dt)``, not an exponential, so fitting an exponential to it
    would report a number the physics does not have while looking more precise than the crossing.

    Args:
        t_s: Sample times [s], strictly increasing, starting at 0.
        amplitude: Amplitude at each sample; normalised internally by its first element.
        rise_tol: Fractional rise tolerated before the curve is rejected as non-monotone.

    Returns:
        The ``1/e`` crossing time, or ``None`` if the curve never reaches it.

    Raises:
        ValueError: If the inputs disagree in length, ``t`` is not increasing, the first amplitude is
            not positive, or the curve rises by more than ``rise_tol`` (a relaxation that gains
            amplitude is a defect, not a slow relaxation).
    """
    t = np.asarray(t_s, dtype=float)
    y = np.asarray(amplitude, dtype=float)
    if t.shape != y.shape or t.ndim != 1 or t.size < 2:
        raise ValueError("t and amplitude must be one-dimensional arrays of equal length >= 2")
    if not np.all(np.diff(t) > 0.0):
        raise ValueError("sample times must be strictly increasing")
    if not (y[0] > 0.0):
        raise ValueError("the initial amplitude must be positive to normalise against")
    y = y / y[0]
    if np.any(np.diff(y) > rise_tol):
        raise ValueError(
            "amplitude rises during the decay — a relaxation that gains amplitude is a defect "
            "(check the CFL bound and the source term), not a slow relaxation"
        )
    below = np.nonzero(y <= ONE_OVER_E)[0]
    if below.size == 0:
        return None
    i = int(below[0])
    if i == 0:
        return float(t[0])
    y0, y1 = y[i - 1], y[i]
    if y0 <= 0.0 or y1 <= 0.0:
        return float(t[i])
    # linear in log y: exact for a single exponential, second order otherwise
    f = (math.log(y0) - math.log(ONE_OVER_E)) / (math.log(y0) - math.log(y1))
    return float(t[i - 1] + f * (t[i] - t[i - 1]))


def power_law_exponent(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
) -> tuple[float, float, float]:
    r"""Fit ``y = C x^p`` in log-log and return ``(p, C, r2)``.

    The exponent is the C-1 scored scalar.  ``r2`` comes back with it so the caller can apply the
    fit-quality floor it declared before the run — a slope through points that do not lie on a line is
    a number, and without ``r2`` nothing downstream can tell that it is not a power law.

    Args:
        x: Sweep variable, all strictly positive (a probe radius, a diffusivity, a grid spacing).
        y: Measured times, all strictly positive.

    Returns:
        ``(p, C, r2)`` — exponent, prefactor in the same units as ``y``, and the coefficient of
        determination of the log-log fit.

    Raises:
        ValueError: If fewer than three points are supplied (two points always fit a line exactly, so
            ``r2`` would be a vacuous 1.0), or if any value is non-positive.
    """
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    if xs.shape != ys.shape or xs.ndim != 1:
        raise ValueError("x and y must be one-dimensional arrays of equal length")
    if xs.size < 3:
        raise ValueError(
            "a power-law fit needs at least three points; two always fit a line exactly and would "
            "report r2 = 1.0 with no evidence behind it"
        )
    if np.any(xs <= 0.0) or np.any(ys <= 0.0):
        raise ValueError("log-log fitting requires strictly positive x and y")
    lx, ly = np.log(xs), np.log(ys)
    p, log_c = np.polyfit(lx, ly, 1)
    residual = ly - (p * lx + log_c)
    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    return float(p), float(math.exp(log_c)), float(r2)
