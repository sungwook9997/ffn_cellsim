"""Conservative RAD transport — closed-form ground truth for the I1c G-actin monomer field.

Host-side acceptance oracles (pure NumPy/SciPy, NO Warp/simulation runtime). I1c advances the monomer
economy on the SAME moving cell domain as I1a/I1b:

    d(phi c)/dt + div(phi c v_f - phi D_c grad c) = R_polymer          (advect with v_f, NOT the discharge q)

with barbed-end consumption + pointed-end/depoly release wired so total actin is conserved. These closed
forms gate three physics facts (the discrete conservation to machine precision is gated by the companion
``transport_reference.RADTransportReference``):

  * FRAP: a bleached cosine pattern recovers at rate ``D_c |k|^2`` (reads back D_c cleanly); the literature
    Soumpasis circular-bleach recovery curve is provided for the figure overlay.
  * ADVECTION follows ``v_f``: a passive front travels at the ABSOLUTE pore-fluid velocity ``v_f``, NOT the
    Darcy discharge ``q = phi(v_f - v_s)`` — using ``q`` would give the wrong front speed (the classic
    solute-velocity bug the build plan warns against).
  * TREADMILL balance: at steady treadmilling barbed consumption == pointed release, so ``integral R dV = 0``
    and the free-monomer pool is stationary.

Sanity Gate (self-tested in tests/ac/fluid/test_i1c_transport_oracle.py): the cosine recovery rate equals
the diffusion eigenvalue; Soumpasis is monotone 0->1 with the right half-time; the front position is linear
in t at slope v_f; treadmill source integrates to zero.

References: Axelrod 1976 Biophys J 16:1055; Soumpasis 1983 Biophys J 41:95 (FRAP); Bear 1972 (advection).
D_c: I0-B1c draft ~2-6 um^2/s (KB-DRAFT-7-03, crowded cytoplasm — NOT free-solution 50-80). c_0: MCF7 GAP.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.special import iv

__all__ = [
    "frap_cosine_recovery_rate",
    "frap_soumpasis_recovery",
    "advection_front_position",
    "treadmill_net_source",
]


def frap_cosine_recovery_rate(wavevector: npt.ArrayLike, d_c: float) -> float:
    """Recovery rate ``lam = D_c |k|^2`` of a bleached cosine concentration mode (reads back D_c).

    Args:
        wavevector: Component wavenumbers ``k_i`` (per length).
        d_c: Monomer diffusivity (length^2/time).

    Returns:
        Recovery/relaxation rate ``lam`` (per time).
    """
    k = np.asarray(wavevector, dtype=np.float64)
    return float(d_c * np.sum(k * k))


def frap_soumpasis_recovery(t: npt.ArrayLike, spot_radius: float, d_c: float) -> npt.NDArray[np.float64]:
    """Soumpasis circular-bleach fluorescence recovery fraction ``f(t)`` (literature FRAP curve).

        f(t) = exp(-2 tau/t) [ I0(2 tau/t) + I1(2 tau/t) ] ,   tau = w^2 / (4 D_c)

    f(0)=0, f(inf)=1; the half-recovery time is ``t_1/2 ~= 0.88 w^2 / (4 D_c)``.

    Args:
        t: Time(s) since bleach (> 0).
        spot_radius: Bleach spot radius ``w`` (length).
        d_c: Monomer diffusivity.

    Returns:
        Recovery fraction ``f(t)`` in [0, 1], same shape as ``t``.
    """
    tt = np.asarray(t, dtype=np.float64)
    tau = spot_radius * spot_radius / (4.0 * d_c)
    with np.errstate(over="ignore", divide="ignore"):
        x = 2.0 * tau / np.where(tt > 0.0, tt, np.nan)
        f = np.exp(-x) * (iv(0, x) + iv(1, x))
    return np.where(tt > 0.0, f, 0.0)


def advection_front_position(x0: float, v_f: float, t: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Passive-front position ``x(t) = x0 + v_f t`` — travels at the ABSOLUTE pore-fluid velocity.

    NOT ``x0 + q t``: the tracer rides ``v_f = v_s + q/phi``, not the relative discharge ``q``. For a moving
    skeleton (``v_s != 0``) these differ by ``v_s + q/phi - q = v_s + q(1-phi)/phi``.

    Args:
        x0: Initial front position (length).
        v_f: Absolute pore-fluid velocity (length/time).
        t: Time(s).

    Returns:
        Front position(s) ``x(t)``.
    """
    return x0 + v_f * np.asarray(t, dtype=np.float64)


def treadmill_net_source(barbed_consumption: float, pointed_release: float) -> float:
    """Net monomer source ``integral R dV = pointed_release - barbed_consumption`` at steady treadmill.

    At balanced treadmilling the two are equal so the net is 0 (the free pool is stationary while monomers
    cycle barbed->bound->pointed->free). Both are counts per unit time (>= 0).

    Args:
        barbed_consumption: Total barbed-end consumption rate (monomers/time), >= 0.
        pointed_release: Total pointed-end/depoly release rate (monomers/time), >= 0.

    Returns:
        The net field source ``integral R dV`` (monomers/time); 0 at treadmill balance.
    """
    if barbed_consumption < 0.0 or pointed_release < 0.0:
        raise ValueError("consumption and release rates must be non-negative")
    return pointed_release - barbed_consumption
