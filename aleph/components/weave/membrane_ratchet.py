r"""Elastic Brownian ratchet — barbed-end polymerization vs membrane load (host NumPy oracle) — I4 / §3A.b.

The lamellipodium protrudes because dendritic-network barbed ends rectify their thermal (and elastic)
fluctuations against the plasma membrane: a monomer intercalates only when a thermal gap of one subunit
(``delta``) opens between the barbed end and the load, so the polymerization velocity is ATTENUATED by the
Boltzmann work ``f*delta`` needed to open that gap under an opposing force ``f`` (Peskin-Odell-Oster 1993
thermal ratchet; Mogilner-Oster 1996 elastic ratchet). The single-barbed-end force-velocity law is

    V(f) = delta * ( k_on_c * exp(-f*delta / kT) - k_off )                 [um/s]

  * FREE (f = 0):  V0 = delta*(k_on_c - k_off)  — the load-free elongation speed.
  * STALL (V = 0): f_s = (kT/delta) * ln(k_on_c / k_off)  — the per-filament stall force (finite, sourced
    entirely by kT, delta, and the on/off ratio — NOT a tuned magnitude).
  * SUPER-STALL (f -> inf): V -> -delta*k_off — the barbed end depolymerizes; the law is reported, never
    clamped (physiological-baseline rule: the model may recede under an over-load, that is real).

The load-bearing property is MONOTONE force-velocity: raising the opposing load slows protrusion, and — the
Mogilner-Oster ADAPTIVE insight the lamellipodium loop exploits — a branched network shares one membrane load
across many barbed ends, so autocatalytic Arp2/3 nucleation (``ac.weave.nucleation``) recruits more filaments,
drops the PER-FILAMENT load ``W/N``, and holds velocity up. This module is the closed-form force-velocity
oracle; ``ac.weave.lamellipodium.AdaptiveLamellipodium`` drives it over time, and the native realization is the
lead's Warp ratchet body-force (the historical ``ff.motility_warp`` ratchet, reused under the network).

Magnitudes: ``delta`` (actin monomer axial half-rise) is a GROUNDED geometric constant; ``kT`` is the sourced
thermal energy at 310 K (reused from ``ac.weave.branch_angle``). The elongation ``k_on_c = k_on_barbed * c``
and dissociation ``k_off`` are the I0-B4 barbed-end kinetics GAP (params_i0b4.yaml consumes the I1c monomer
field ``c`` read-only) — passed in, NEVER chosen to hit a protrusion/retrograde band (§6.2 report-not-tune).
The SHAPE gates below are magnitude-independent (they sweep k_on_c/k_off/f as oracle variables), so the
analytic build is unblocked exactly as the branch-angle gate is.

engine units: length um, force pN, velocity um/s, energy pN*um; kT_310K = 4.28e-3 pN*um.

Sanity Gate (self-tested in tests/ac/weave/test_membrane_ratchet_oracle.py):
  * dimensional: V in um/s; the ratchet argument f*delta/kT is dimensionless (pN*um / (pN*um)).
  * boundary: V(0) = V0 = delta*(k_on_c-k_off); V -> -delta*k_off as f -> inf (super-stall depolymerization).
  * monotonicity: V(f) is strictly decreasing in the opposing load f (force-velocity, never increasing).
  * stall: V(f_s) = 0 at f_s = (kT/delta)*ln(k_on_c/k_off); requires k_on_c > k_off (a growing end) else the
    filament recedes even free and stall is undefined (raised, not floored to a convenient value).
  * ensemble/adaptive: for fixed total load W, per-filament load W/N decreases in N, so the ensemble
    protrusion velocity is monotone INCREASING in the number of load-sharing barbed ends (Mogilner-Oster).

Reference: Peskin, Odell & Oster 1993 Biophys J 65:316 (thermal ratchet); Mogilner & Oster 1996 Biophys J
71:3030 (elastic ratchet, force-velocity of a branched array); actin axial rise ~2.7 nm/subunit (Pollard).
"""

from __future__ import annotations

import numpy as np

from aleph.components.weave.branch_angle import KT_310K_PN_UM

__all__ = [
    "ACTIN_MONOMER_DELTA_UM",
    "KT_310K_PN_UM",
    "ratchet_velocity",
    "free_polymerization_velocity",
    "stall_force",
    "per_filament_load",
    "ensemble_protrusion_velocity",
]

# Actin monomer axial half-rise: one subunit advances a two-start filament by ~2.7 nm along its axis
# (13 subunits / 37 nm helical repeat ~ 2.75 nm/subunit; the Mogilner-Oster ratchet step size). A GROUNDED
# geometric constant (like ARP23_THETA0_RAD) — the ratchet gap width, not a tuning knob.
ACTIN_MONOMER_DELTA_UM: float = 0.0027


def ratchet_velocity(
    load_pN: float | np.ndarray,
    k_on_c: float,
    k_off: float,
    *,
    delta_um: float = ACTIN_MONOMER_DELTA_UM,
    kT: float = KT_310K_PN_UM,
) -> np.ndarray:
    """Single-barbed-end force-velocity ``V = delta*(k_on_c*exp(-f*delta/kT) - k_off)`` [um/s].

    Args:
        load_pN: opposing membrane load ``f`` on the barbed end [pN]; must be >= 0 (a compressive load).
        k_on_c: elongation rate ``k_on_barbed * c`` [1/s] (I0-B4 GAP; consumes the I1c monomer field c).
        k_off: barbed-end dissociation rate [1/s] (I0-B4 GAP).
        delta_um: monomer axial step / ratchet gap [um] (grounded geometric constant).
        kT: thermal energy [pN*um] at the physiological setpoint.

    Returns:
        Polymerization velocity [um/s]; positive = protruding, negative = super-stall depolymerization.
    """
    load = np.asarray(load_pN, float)
    if np.any(load < 0.0):
        raise ValueError("opposing load must be >= 0 (a compressive membrane load)")
    if k_on_c < 0.0 or k_off < 0.0:
        raise ValueError("k_on_c and k_off must be >= 0")
    if delta_um <= 0.0 or kT <= 0.0:
        raise ValueError("delta_um and kT must be > 0")
    return delta_um * (k_on_c * np.exp(-load * delta_um / kT) - k_off)


def free_polymerization_velocity(
    k_on_c: float,
    k_off: float,
    *,
    delta_um: float = ACTIN_MONOMER_DELTA_UM,
) -> float:
    """Load-free elongation velocity ``V0 = delta*(k_on_c - k_off)`` [um/s] (the ratchet at f = 0)."""
    return float(ratchet_velocity(0.0, k_on_c, k_off, delta_um=delta_um))


def stall_force(
    k_on_c: float,
    k_off: float,
    *,
    delta_um: float = ACTIN_MONOMER_DELTA_UM,
    kT: float = KT_310K_PN_UM,
) -> float:
    """Per-filament stall force ``f_s = (kT/delta)*ln(k_on_c/k_off)`` [pN] where ``V(f_s) = 0``.

    Requires a growing end (``k_on_c > k_off``); a receding end (``k_on_c <= k_off``) has no protruding
    regime and stall is undefined — raised rather than floored to a convenient value (report-not-tune).
    """
    if delta_um <= 0.0 or kT <= 0.0:
        raise ValueError("delta_um and kT must be > 0")
    if k_off <= 0.0:
        raise ValueError("k_off must be > 0 for a finite stall force")
    if k_on_c <= k_off:
        raise ValueError(
            "no stall force: k_on_c <= k_off means the barbed end recedes even at zero load "
            "(not a protruding filament)"
        )
    return float((kT / delta_um) * np.log(k_on_c / k_off))


def per_filament_load(total_load_pN: float, n_barbed: int) -> float:
    """Share a total membrane load ``W`` equally across ``n_barbed`` pushing barbed ends: ``W / N`` [pN]."""
    if total_load_pN < 0.0:
        raise ValueError("total load must be >= 0")
    if n_barbed <= 0:
        raise ValueError("n_barbed must be >= 1 (at least one pushing barbed end)")
    return float(total_load_pN / n_barbed)


def ensemble_protrusion_velocity(
    total_load_pN: float,
    n_barbed: int,
    k_on_c: float,
    k_off: float,
    *,
    delta_um: float = ACTIN_MONOMER_DELTA_UM,
    kT: float = KT_310K_PN_UM,
) -> float:
    """Protrusion velocity of ``n_barbed`` load-sharing ends under a total load ``W`` [um/s].

    Each barbed end feels ``W / N`` (equal load sharing) and moves at the ratchet velocity for that load; the
    rigidly-coupled membrane advances at that common speed. Because ``W / N`` decreases in ``N``, the velocity
    is monotone increasing in the number of barbed ends — the Mogilner-Oster adaptive force-velocity that
    autocatalytic branching exploits (more branches -> lower per-filament load -> faster protrusion).
    """
    w = per_filament_load(total_load_pN, n_barbed)
    return float(ratchet_velocity(w, k_on_c, k_off, delta_um=delta_um, kT=kT))
