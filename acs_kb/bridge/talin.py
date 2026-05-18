"""Talin Bell-Evans unfolding (KU-2.6, Phase 1 1-state simplification).

Talin's 13 rod domains (R1–R13) each unfold under tension at a Bell-Evans
rate (del Rio 2009 Science, Yao 2014 Sci Rep, Tapia-Rojo 2019 Sci Adv).
For Phase 1 Unit 2.2 we collapse the 13-state ladder to a single binary
domain per FA: folded (0) → unfolded (1). The unfolding rate is::

    k_unfold(F) = k_{u,0} · exp(F · Δx* / k_B T)

with the *force on a representative talin* coming from the per-clutch
force at the FA (see ``motor_clutch.step``). Once unfolded, the domain
stays unfolded — Phase 1 ignores refolding (KU-2.6 pitfall: refolding is
non-negligible at low force; deferred to Phase 2 along with the full
13-state ladder).

A single unfolded talin exposes one vinculin binding site (VBS), which
:mod:`acs_kb.bridge.vinculin` consumes (KU-2.7 input).

Sanity Gate
-----------
1. Dimensional analysis: ``F`` [N] · ``Δx*`` [m] / ``k_B T`` [J] is
   dimensionless. ``k_u0`` and ``k_unfold`` are in s⁻¹. No CFL.
2. Boundary cases: ``F = 0`` ⇒ ``k_unfold = k_u0`` (resting rate);
   ``F → ∞`` is clipped by ``_EXP_CAP = 50`` to avoid float overflow
   (corresponds to F ≈ 143 pN at default Δx*, well past biological).
   ``fa.talin_unfolded_domains ≥ n_domains`` short-circuits to no-op.
3. Conservation: stochastic state-transition only; no energy term in
   Phase 1 (the work done by the force during unfolding is bookkept in
   the motor-clutch energy balance, not here).
4. Numerical sanity: float64; pure NumPy/math; one ``rng.random()``
   draw per (FA, dt) step.
5. Sign / sense: positive ``F`` (tensile on talin) increases
   ``k_unfold`` ⇒ stabilises the unfolded state under load. Negative
   ``F`` is non-physical for tension and is clamped to 0.
6. Measurement protocol: the unit test sweeps fixed ``F`` and confirms
   the mean time-to-unfold matches ``1/k_unfold(F)`` from Gillespie
   sampling, exactly at the analytical proof point.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from acs_kb.bridge.types import FocalAdhesion

# k_B · T at 310 K (37 °C), SI Joules — matches phase1_unit1.yaml `kT`.
_KT_310K = 4.28e-21
# Cap exponent argument; exp(50) ≈ 5·10²¹ — far above any plausible rate.
_EXP_CAP = 50.0


@dataclass(frozen=True, slots=True)
class TalinParams:
    """KU-2.6 single-state Bell-Evans parameters."""

    k_u0: float = 0.01            # s⁻¹  resting unfold rate
    dx_star: float = 1.5e-9       # m    transition state distance
    kT: float = _KT_310K          # J    thermal energy at 310 K
    n_domains: int = 1            # Phase 1 1-state binary

    @classmethod
    def from_config(cls, cfg: dict) -> "TalinParams":
        b = cfg["bridge"] if "bridge" in cfg else cfg
        t = b["talin"]
        return cls(
            k_u0=float(t["k_u0"]),
            dx_star=float(t["dx_star"]),
            kT=float(t.get("kT", _KT_310K)),
            n_domains=int(t.get("n_domains", 1)),
        )


DEFAULT_PARAMS = TalinParams()


def talin_unfold_rate(
    force: float | np.ndarray, params: TalinParams = DEFAULT_PARAMS
) -> np.ndarray:
    """Bell-Evans unfolding rate ``k_u(F) = k_{u,0} exp(F Δx*/kT)`` [s⁻¹].

    Negative forces are clamped to 0 (tension only). The exponent is
    capped at ``_EXP_CAP`` to avoid overflow at unphysically large F.
    """
    F = np.maximum(np.asarray(force, dtype=np.float64), 0.0)
    arg = np.clip(F * params.dx_star / params.kT, -_EXP_CAP, _EXP_CAP)
    return params.k_u0 * np.exp(arg)


def talin_unfold_step(
    fa: FocalAdhesion,
    force_per_talin: float,
    dt: float,
    rng: np.random.Generator,
    params: TalinParams = DEFAULT_PARAMS,
) -> int:
    """Stochastic 1-state Bell-Evans unfold for the given FA.

    Phase 1 has a *single binary* talin per FA; this function either
    flips it from folded (0) to unfolded (1) with probability
    ``1 − exp(−k_unfold(F) · dt)``, or leaves it alone. No refolding.

    Returns
    -------
    n_new_unfolded : int
        Number of domains that transitioned this step (0 or 1).
    """
    if fa.talin_unfolded_domains >= params.n_domains:
        return 0
    k_u = float(talin_unfold_rate(force_per_talin, params))
    if k_u * dt <= 0.0:
        return 0
    p_unfold = -math.expm1(-k_u * dt)              # numerically stable 1 − exp(−x)
    if rng.random() < p_unfold:
        fa.talin_unfolded_domains += 1
        return 1
    return 0
