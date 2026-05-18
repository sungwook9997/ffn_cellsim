"""Pereverzev two-pathway catch-slip off-rate (KU-2.5).

The integrin–ligand bond dissociates through two competing pathways: a
slip pathway (off-rate grows with force) and a catch pathway (off-rate
falls with force). The two-pathway Pereverzev 2005 form gives a closed
expression for the total off-rate::

    k_off(F) = k_s · exp(F / F_s) + k_c · exp(−F / F_c)

with bond lifetime ``τ(F) = 1 / k_off(F)``. ``τ(F)`` is non-monotonic
with a single maximum at::

    F* = (F_s · F_c) / (F_s + F_c) · ln(k_c · F_s / (k_s · F_c))

(set ``dτ/dF = 0``, equivalent to ``dk_off/dF = 0``; require the
argument of the log positive — i.e. ``k_c · F_s > k_s · F_c`` —
otherwise the bond is slip-only and there is no catch peak).

Sanity Gate
-----------
1. Dimensional analysis: ``F``, ``F_s``, ``F_c`` in N; ``k_s``, ``k_c``
   in s⁻¹; ``k_off`` in s⁻¹; ``τ`` in s. The ``F/F_s`` and ``F/F_c``
   exponents are dimensionless.
2. Boundary cases: ``F = 0`` → ``k_off = k_s + k_c`` (finite);
   ``F → ∞`` → slip term dominates → ``τ → 0``;
   negative ``F`` clamped to 0 (compression on an integrin bond is not
   physical — clutches operate in tension only).
3. Conservation: rate equation only; no energy term.
4. Numerical sanity: float64; ``exp`` clamped against the slip blow-up
   at unphysically large forces (``F > 200 pN`` is non-physical for
   integrin tension and gives ``exp(200/F_s)`` overflow at F_s ≈ 30 pN).
5. Sign / sense: slip term grows with ``F`` (force destabilises bond);
   catch term falls with ``F`` (force initially stabilises bond).
6. Measurement protocol: the regression test samples lifetimes via
   Gillespie at each fixed ``F`` and matches ``⟨τ_simulated⟩`` against
   ``1/k_off(F)``; the catch peak is located by argmax over a fine
   force grid and compared to the closed-form ``F*`` (not to the
   experimentally quoted F* ≈ 30 pN, which arises from a different fit;
   see REPORT.md notes).

KU-2.18 defaults (Phase 1)
--------------------------
- ``k_off^slip  = 0.5 s⁻¹``,  ``F_s = 30 pN``
- ``k_off^catch = 0.4 s⁻¹``,  ``F_c = 7 pN``

With those values the closed-form ``F*`` is ≈ 7 pN (not the 30 pN
quoted experimentally in KU-2.5 from Kong 2009 α5β1 fits). The
discrepancy is logged in the Unit 2.1 REPORT.md as a known
parameter-fit gap, not a code bug.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# Cap exponent argument to avoid float overflow on absurdly large forces
# (~3·10^2). exp(700) is the float64 ceiling; exp(50) ≈ 5·10^21 is more
# than enough headroom for any realistic integrin force.
_EXP_CAP = 50.0


@dataclass(frozen=True, slots=True)
class CatchSlipParams:
    """Pereverzev two-pathway parameters (KU-2.5 / KU-2.18)."""

    k_off_slip: float = 0.5
    F_s: float = 30.0e-12
    k_off_catch: float = 0.4
    F_c: float = 7.0e-12


DEFAULT_PARAMS = CatchSlipParams()


def catch_slip_off_rate(
    F: float | np.ndarray, params: CatchSlipParams = DEFAULT_PARAMS
) -> np.ndarray:
    """Pereverzev two-pathway off-rate ``k_off(F)``  [s⁻¹]  (KU-2.5)."""
    F_arr = np.maximum(np.asarray(F, dtype=np.float64), 0.0)
    a_s = np.clip(F_arr / params.F_s, -_EXP_CAP, _EXP_CAP)
    a_c = np.clip(F_arr / params.F_c, -_EXP_CAP, _EXP_CAP)
    return params.k_off_slip * np.exp(a_s) + params.k_off_catch * np.exp(-a_c)


def catch_slip_lifetime(
    F: float | np.ndarray, params: CatchSlipParams = DEFAULT_PARAMS
) -> np.ndarray:
    """Bond lifetime ``τ(F) = 1 / k_off(F)``  [s]."""
    return 1.0 / catch_slip_off_rate(F, params)


def lifetime_peak_force(params: CatchSlipParams = DEFAULT_PARAMS) -> float:
    """Closed-form catch peak ``F*``  [N]; ``nan`` if the bond is slip-only.

    Derived by setting ``dk_off/dF = 0``::

        F* = F_s F_c / (F_s + F_c) · ln(k_c F_s / (k_s F_c))

    The bond is catch-then-slip only if ``k_c F_s > k_s F_c`` (so the
    log argument exceeds 1 and ``F*`` is positive); otherwise the slip
    pathway dominates at all forces and there is no interior maximum.
    """
    p = params
    arg = (p.k_off_catch * p.F_s) / (p.k_off_slip * p.F_c)
    if arg <= 1.0:
        return float("nan")
    return (p.F_s * p.F_c) / (p.F_s + p.F_c) * math.log(arg)


def lifetime_peak_value(params: CatchSlipParams = DEFAULT_PARAMS) -> float:
    """Closed-form maximum bond lifetime ``τ(F*)``  [s]; ``nan`` if slip-only."""
    F_star = lifetime_peak_force(params)
    if math.isnan(F_star):
        return float("nan")
    return float(catch_slip_lifetime(F_star, params))
