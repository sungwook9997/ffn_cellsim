r"""Ensemble-stall EMERGENCE — the analytic reference for the per-minifilament contractile force.

Host-side acceptance oracle (pure NumPy, NO Warp). This is the gate the build plan calls out explicitly:
"ensemble stall emerges from bound heads (NOT imposed ``N_side * F_head``)" (I3, §3 line 243; §3 Session D).

The whole point of the head-resolved primitive (P2) is that the minifilament's contractile force is NOT a
lumped constant and NOT the naive product ``N_side * F_head`` (all heads engaged). It EMERGES from summing the
individual bound heads, whose engaged count is itself set by the Bell/attach kinetics:

    at ensemble isometric (v_slide -> 0) each BOUND head steps to its own stall F_head (v_head=0), so
        F_ensemble = n_bound * F_head
    with n_bound ~ Binomial(N_side, phi_b),  phi_b = k_on / (k_on + k_off0 * exp(F_head / f0))

Hence the mean-field emergent stall is

    E[F_ensemble] = N_side * phi_b(F_head) * F_head        <=  N_side * F_head    (equality only if phi_b=1)

The strict inequality is the "Bell-emergent engaged fraction" that SURVIVES at the settled state (P2): fewer
than all heads are bound because the slip off-rate rises under load, so the emergent stall is BELOW the naive
product. This module returns BOTH the mean-field value and a stochastic (per-realisation) ensemble so the gate
can check the stochastic mean converges to the mean field — i.e. the force is assembled from heads, not set.

⚠ report-not-tune (§6.2): where the emergent stall lands (e.g. vs the KB-3.18 50-100 pN per-minifilament band)
is a FINDING for the GAP params (F_head, N_side, k_on GAP to PI). NEVER add heads/density to move it into a band.
⚠ engaged-fraction caveat: with the archived k_on=50/s >> k_off0=0.35/s and f0~7 pN >> F_head~0.5-2 pN, phi_b
is ~0.99 at per-head stall load, so the emergent stall is only slightly below the naive product HERE — but it
is COMPUTED from kinetics, not assumed. If a future isoform has F_head ~ f0, the self-limiting bites hard. The
oracle exposes this dependence; it does not bake in 0.99. (The measured Kovacs duty ~0.1 differs from this
k_on/k_off equilibrium 0.99 — flagged in params_i0b3.yaml::k_on / duty as a PI reconciliation item.)

Sanity Gate (of the oracle itself — tests/ac/motor/test_ensemble_stall_oracle.py):
  * emergence: mean-field E[F] = N_side * phi_b * F_head, and phi_b in (0,1] is computed from kinetics;
  * strict bound: E[F] <= N_side * F_head, with equality iff phi_b == 1 (k_off0 -> 0 or k_on -> inf);
  * stochastic mean -> mean-field as n_realisations -> inf (law of large numbers; within a few SEM);
  * self-limiting monotonicity: E[F] decreases if k_off0 rises or f0 falls (more shedding under load);
  * NOT-imposed arbiter: changing k_on/k_off0 changes E[F] (a hard-coded N_side*F_head would not move).

References:
  Stam-Hocky 2015 Biophys J 108:1997 (minifilament); Billington 2013 JBC (N per filament); Kovacs 2003 (duty).
  KB-3.18 (per-minifilament 50-100 pN band). Bell 1978 (slip). Freedman 2017 (AFINES ensemble).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.components.motor.bell_kinetics_analytic import engaged_fraction_steady

__all__ = [
    "ensemble_stall_meanfield",
    "ensemble_stall_stochastic",
]


def ensemble_stall_meanfield(
    n_side: int,
    f_stall_head: float,
    k_on: float,
    k_off0: float,
    f0: float,
) -> dict[str, float]:
    """Mean-field emergent per-side ensemble stall force from the Bell-governed engaged head count.

    Args:
        n_side: explicit heads per half-filament ``N_side`` (I0-B3 GAP: 10 AFINES vs 28-30 Billington).
        f_stall_head: per-head isometric stall ``F_head`` [pN] (I0-B3 GAP: 0.5 vs 2.0 pN).
        k_on: per-head attach rate [1/s].
        k_off0: per-head zero-force off-rate [1/s].
        f0: Bell characteristic force [pN].

    Returns:
        dict with:
          ``phi_b`` — steady-state engaged fraction at per-head load ``F_head`` (in (0, 1]);
          ``n_engaged`` — expected engaged head count ``N_side * phi_b``;
          ``f_ensemble`` — emergent ensemble stall ``N_side * phi_b * F_head`` [pN];
          ``f_naive`` — the naive (imposed) product ``N_side * F_head`` [pN], for the strict-bound check.
    """
    if n_side < 1:
        raise ValueError("n_side must be >= 1")
    if f_stall_head <= 0.0:
        raise ValueError("f_stall_head must be positive")
    phi_b = float(engaged_fraction_steady(k_on, f_stall_head, k_off0, f0))
    n_engaged = n_side * phi_b
    return {
        "phi_b": phi_b,
        "n_engaged": n_engaged,
        "f_ensemble": n_engaged * f_stall_head,
        "f_naive": n_side * f_stall_head,
    }


def ensemble_stall_stochastic(
    n_side: int,
    f_stall_head: float,
    k_on: float,
    k_off0: float,
    f0: float,
    n_realizations: int,
    rng: np.random.Generator,
) -> npt.NDArray[np.float64]:
    """Stochastic per-realisation ensemble stall ``F = n_bound * F_head``, ``n_bound ~ Binomial(N_side, phi_b)``.

    Each of ``N_side`` heads is independently bound with probability ``phi_b(F_head)``; the contractile force
    is the sum over bound heads (all at their per-head stall at isometric). The ensemble MEAN of the returned
    array converges to :func:`ensemble_stall_meanfield`'s ``f_ensemble`` — demonstrating that the force is
    ASSEMBLED from heads, not imposed.

    Args:
        n_side: heads per half-filament.
        f_stall_head: per-head stall [pN].
        k_on, k_off0, f0: attach rate / zero-force off-rate / Bell force (see :func:`ensemble_stall_meanfield`).
        n_realizations: number of independent minifilament draws.
        rng: NumPy random generator (caller-seeded for reproducibility).

    Returns:
        ``(n_realizations,)`` array of per-realisation ensemble stall forces [pN].
    """
    if n_realizations < 1:
        raise ValueError("n_realizations must be >= 1")
    phi_b = float(engaged_fraction_steady(k_on, f_stall_head, k_off0, f0))
    n_bound = rng.binomial(n_side, phi_b, size=n_realizations)
    return n_bound.astype(np.float64) * f_stall_head
