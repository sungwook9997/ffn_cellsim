"""The creep observable must not reference the initial flux, because the initial flux is the mesh.

The first creep observable was the 1/e crossing of dQ/dt, and it failed its own grid-invariance check
by a factor of 21 (0.03405 / 0.01830 / 0.01444 s at dx = 0.30 / 0.20 / 0.15). The reason is physics,
not numerics: the flux out of a patch held at fixed potential does NOT decay to zero. For a sphere in
an infinite diffusive medium it tends to 4*pi*D*a, approached as

    flux(t) = 4 pi D a (1 + a / sqrt(pi D t))

so the transient excess falls as t^(-1/2) and there is no intrinsic 1/e time at all. The old
observable measured "one part in e of the INITIAL flux", and the initial flux is the gradient at the
patch edge, which scales as 1/dx.

These gates check the replacement against that closed form: the time at which the excess equals the
steady value, tau = a^2/(pi D).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aleph.scripts.c1_biot_creep import CREEP_ORACLE_C, FIT_SKIP_FRACTION


def _analytic_flux(t: np.ndarray, a: float, d: float) -> np.ndarray:
    """flux(t) for a sphere of radius a held at unit potential in an infinite medium."""
    return 4.0 * math.pi * d * a * (1.0 + a / np.sqrt(math.pi * d * t))


def _tau_from_excess(tm: np.ndarray, flux: np.ndarray) -> float | None:
    """The driver's inversion, reimplemented here so the gate tests the arithmetic, not the driver.

    Fit ``flux = A + B t^(-1/2)``; the crossover ``B/sqrt(t) = A`` gives ``tau = (B/A)^2``.
    """
    skip = max(2, int(FIT_SKIP_FRACTION * flux.size))
    if flux.size - skip < 8 or not np.all(flux[skip:] > 0):
        return None
    b_slope, a_int = np.polyfit(1.0 / np.sqrt(tm[skip:]), flux[skip:], 1)
    if not (a_int > 0.0 and b_slope > 0.0):
        return None
    return float((b_slope / a_int) ** 2)


def test_the_oracle_prefactor_is_one_over_pi() -> None:
    assert CREEP_ORACLE_C == pytest.approx(1.0 / math.pi, rel=1e-15)
    assert CREEP_ORACLE_C == pytest.approx(0.318310, rel=1e-5)


@pytest.mark.parametrize(("a", "d"), [(0.5, 50.0), (1.0, 50.0), (2.0, 12.5)])
def test_the_analytic_flux_gives_back_a2_over_pi_D(a: float, d: float) -> None:
    """The inversion recovers tau = a^2/(pi D) from the closed-form flux it was derived from."""
    tau_true = a * a / (math.pi * d)
    t = np.geomspace(0.02 * tau_true, 40.0 * tau_true, 4000)
    tau = _tau_from_excess(t, _analytic_flux(t, a, d))
    assert tau is not None
    assert tau == pytest.approx(tau_true, rel=0.05)


def test_the_observable_does_not_move_when_the_initial_flux_does() -> None:
    """The whole point: truncating the early curve — which is what refining dx changes — must not move it.

    A finer mesh resolves a larger initial flux, so the recorded curve starts higher and earlier. The
    old 1/e-of-initial observable moved by 2.36x under exactly this change; this one must not.
    """
    a, d = 1.0, 50.0
    tau_true = a * a / (math.pi * d)
    full = np.geomspace(0.005 * tau_true, 40.0 * tau_true, 6000)
    taus = [_tau_from_excess(full[k:], _analytic_flux(full[k:], a, d)) for k in (0, 500, 1000)]
    assert all(t is not None for t in taus)
    assert max(taus) / min(taus) == pytest.approx(1.0, abs=0.02)


def test_a_flux_that_never_exceeds_its_steady_value_returns_none() -> None:
    """No transient means no creep time; returning the first sample would be a number from nothing."""
    t = np.linspace(1.0, 10.0, 100)
    assert _tau_from_excess(t, np.full_like(t, 3.0)) is None
