r"""Arp2/3 dendritic nucleation + capping — flux-limited KMC rate law (host NumPy oracle) — I4.

Arp2/3 autocatalytic dendritic nucleation at existing filament sides, **flux-limited by the I1c G-actin
monomer field** and gated by an **NPF (WASP/WAVE) activity field** at the membrane (NEW_ENGINE_BUILD_PLAN
§3A.b). The per-site branch-nucleation rate is

    k_branch(x) = k_arp0 * npf(x) * c(x) / (c(x) + K_m)                 [1/s]   (Michaelis-Menten in monomer)

so nucleation SATURATES at high monomer and — the load-bearing property — goes to ZERO as the transported
monomer pool ``c`` is depleted (the emergent network cannot grow faster than G-actin is advected in by
``v_f``; pressure-only cannot carry monomer — [[project-field-actuation-and-actin-gf-gap]]). NPF activity is
a declared membrane field (a SEEDED input, ablatable), not a tuned per-branch bias. Capping protein
terminates growth at rate ``k_cap`` (a competing first-order KMC channel, via ``hand_kmc`` semantics), so the
mean filament length at steady state is set by the growth/cap balance.

This oracle is the closed-form rate law the device branch/cap KMC (``ac.weave.crosslink_kmc_warp`` shares the
same attach/detach primitives) reproduces; it CONSUMES the I1c monomer concentration READ-ONLY (a scalar or
sampled array ``c`` — this module never runs the Warp field). ``k_arp0``, ``K_m``, ``k_cap``, and the NPF
areal density are I0-B4 GAPs (params_i0b4.yaml) — surfaced to PI, never chosen to hit a protrusion band.

engine units: rate 1/s, length um, concentration (monomer) uM (or any consistent unit; only c/(c+K_m) enters).

Sanity Gate (self-tested in tests/ac/weave/test_nucleation_oracle.py):
  * dimensional: k_branch in 1/s; c/(c+K_m) in [0,1]; k_branch=0 when npf=0 or c=0 (flux limit + NPF gate).
  * monotonicity: k_branch strictly increasing in c and in npf; saturates -> k_arp0*npf as c -> inf.
  * flux limit: depleting c drives k_branch -> 0 (the emergence cannot outrun monomer supply).
  * capping: steady-state mean length ~ v_grow / k_cap (longer filaments when capping is rarer); barbed-end
    elongation consumes monomer at ``k_on_barbed * c`` (mass-coupled to I1c, conserved there).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "branch_nucleation_rate",
    "capping_prob",
    "steady_state_length",
    "barbed_consumption_rate",
    "activate_dormant_daughters",
]


def branch_nucleation_rate(
    c: float | npt.NDArray[np.float64],
    npf: float | npt.NDArray[np.float64],
    k_arp0: float,
    k_m: float,
) -> npt.NDArray[np.float64]:
    """Arp2/3 branch-nucleation rate ``k_arp0 * npf * c/(c+K_m)`` [1/s], flux-limited by monomer ``c``.

    Args:
        c: local free-monomer concentration (from the I1c field) [uM]; >= 0.
        npf: local NPF (WASP/WAVE) activity in [0, 1] (a declared membrane field, SEEDED).
        k_arp0: max autocatalytic branch rate at saturating monomer + full NPF [1/s] (I0-B4 GAP).
        k_m: Michaelis monomer half-saturation [uM] (I0-B4 GAP).

    Returns:
        Per-site nucleation rate [1/s]; 0 where c = 0 (flux limit) or npf = 0 (no NPF).
    """
    c = np.asarray(c, float)
    npf = np.asarray(npf, float)
    if k_arp0 < 0.0 or k_m <= 0.0:
        raise ValueError("k_arp0 must be >= 0 and k_m > 0")
    if np.any(c < 0.0):
        raise ValueError("monomer concentration c must be >= 0")
    if np.any((npf < 0.0) | (npf > 1.0)):
        raise ValueError("npf activity must be in [0, 1]")
    return k_arp0 * npf * c / (c + k_m)


def capping_prob(tau: float, k_cap: float) -> float:
    """Per-tick capping probability ``1 - exp(-tau k_cap)`` (Poisson; same form as the hand-KMC detach)."""
    if tau < 0.0 or k_cap < 0.0:
        raise ValueError("tau and k_cap must be >= 0")
    return float(1.0 - np.exp(-tau * k_cap))


def steady_state_length(v_grow: float, k_cap: float, seg_um: float) -> float:
    """Mean barbed-end length at the growth/cap balance ``v_grow / k_cap`` [um] (floored at one segment).

    A barbed end elongates at ``v_grow`` [um/s] until a capping event (rate ``k_cap``); the expected pre-cap
    growth is ``v_grow / k_cap``. Longer filaments emerge when capping is rarer — the capping rate SETS the
    dendritic mesh size, it is not a scripted length.
    """
    if v_grow < 0.0 or k_cap <= 0.0 or seg_um <= 0.0:
        raise ValueError("v_grow >= 0, k_cap > 0, seg_um > 0 required")
    return float(max(v_grow / k_cap, seg_um))


def barbed_consumption_rate(c: float, k_on_barbed: float) -> float:
    """Monomer consumption at a growing barbed end ``k_on_barbed * c`` [monomer/s] (mass-coupled to I1c).

    This is the reaction ``R_polymer`` sink the I1c conservative transport applies to ``u = phi c``; every
    consumed monomer becomes one bound polymer subunit so total actin is conserved there (non-vacuous).
    """
    if k_on_barbed < 0.0 or c < 0.0:
        raise ValueError("k_on_barbed and c must be >= 0")
    return float(k_on_barbed * c)


def activate_dormant_daughters(
    active_mask: npt.NDArray[np.bool_],
    rate: npt.NDArray[np.float64],
    tau: float,
    rng: np.random.Generator,
) -> npt.NDArray[np.bool_]:
    """Activate dormant daughter filaments by a Poisson KMC pass (N-FIXED: no node allocation).

    A dormant daughter (``active_mask[f] == False``) with per-site nucleation ``rate[f]`` [1/s] turns active
    with probability ``1 - exp(-tau*rate)`` in the tick. This is the N-fixed encoding of §3A.b / the growing-N
    reconciliation: branch nucleation ACTIVATES pre-allocated dormant daughters (the dendritic topology
    emerges) — it never allocates new nodes (that is the I9+ true growing-N deferral).

    Args:
        active_mask: (F,) current active flags (mutated copy returned).
        rate: (F,) per-filament nucleation rate [1/s] (0 for already-active or unreachable).
        tau: KMC tick [s].
        rng: NumPy generator (device RNG mirror).

    Returns:
        (F,) updated active mask.
    """
    active_mask = np.asarray(active_mask, bool).copy()
    rate = np.asarray(rate, float)
    if active_mask.shape != rate.shape:
        raise ValueError("active_mask and rate must have the same shape")
    dormant = ~active_mask
    p = 1.0 - np.exp(-tau * np.clip(rate, 0.0, None))
    fired = dormant & (rng.random(active_mask.shape[0]) < p)
    active_mask[fired] = True
    return active_mask
