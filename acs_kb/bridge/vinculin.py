"""Vinculin recruitment and allosteric clutch reinforcement (KU-2.7).

Vinculin binds to the cryptic VBS exposed when force-unfolded talin
domains open (KU-2.6). Bound vinculin allosterically increases the
effective integrin–actin spring constant::

    k_int^eff = k_int^bare · (1 + α · N_vin)            (Han 2021 eLife)

and the vinculin count itself follows a recruitment / dissociation
balance::

    dN_vin/dt = k_rec · n_unfolded · N_vin^free − k_diss · N_vin

with ``N_vin^free`` the cytoplasmic free-vinculin pool (held constant
in Phase 1 — KU-2.7 pitfall flags this as inaccurate when many FAs
compete; deferred to Phase 2).

Sanity Gate
-----------
1. Dimensional analysis: ``k_rec`` [s⁻¹ · vin⁻¹], ``k_diss`` [s⁻¹],
   ``n_unfolded`` [domains], ``N_vin^free`` [vin], so the ODE has
   units [vin/s] on both sides. ``α`` is dimensionless. No CFL — the
   ``dt`` for the Euler-step here must be ≪ ``1/max(k_rec n N_free,
   k_diss)`` which is checked at config-resolve time.
2. Boundary cases: ``n_unfolded = 0`` ⇒ recruitment rate 0, so N_vin
   decays exponentially to 0. ``N_vin = 0`` ⇒ effective stiffness
   reduces cleanly to the bare value. Negative N_vin clamped to 0.
3. Conservation: per-FA only — the global vinculin pool is NOT updated
   in Phase 1 (constant ``N_vin^free``). Phase 2 must close this loop.
4. Numerical sanity: float64; **Poisson tau-leap** integration of the
   recruitment / dissociation events — an explicit Euler on an integer
   N_vin field would quantise per-step increments below 0.5 to zero
   and trap the count at 0 (caught by
   ``test_vinculin_recruit_steady_state``). The Poisson tau-leap is
   exact in mean and correct in variance for short ``dt``, and degrades
   to Gillespie at larger rates.
5. Sign / sense: positive ``n_unfolded`` *adds* vinculin; positive
   ``k_diss`` *removes* it; ``α > 0`` strengthens the clutch under load.
6. Measurement protocol: the regression test runs a force-clamp on a
   single FA, lets talin unfold, and verifies that the steady-state
   vinculin count matches the closed-form ``k_rec n N_free / k_diss``
   within a 10 % tolerance (sample mean over the second half).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from acs_kb.bridge.types import FocalAdhesion


@dataclass(frozen=True, slots=True)
class VinculinParams:
    """KU-2.7 recruitment + allostery parameters (Han 2021, Franz 2023)."""

    k_rec: float = 0.1        # s⁻¹ per unfolded domain × vinculin⁻¹·free
    k_diss: float = 0.05      # s⁻¹  dissociation rate
    N_free: float = 100.0     # cytoplasmic free vinculin pool (per-FA, Phase 1)
    alpha: float = 0.02       # allostery factor (k_int_eff = k_int (1 + α N_vin))

    @classmethod
    def from_config(cls, cfg: dict) -> "VinculinParams":
        b = cfg["bridge"] if "bridge" in cfg else cfg
        v = b["vinculin"]
        return cls(
            k_rec=float(v["k_rec"]),
            k_diss=float(v["k_diss"]),
            N_free=float(v["N_free"]),
            alpha=float(v["alpha"]),
        )

    def steady_state_count(self, n_unfolded: int) -> float:
        """Closed-form steady-state ``N_vin = k_rec n N_free / k_diss``."""
        return self.k_rec * n_unfolded * self.N_free / max(self.k_diss, 1e-30)


DEFAULT_PARAMS = VinculinParams()


def vinculin_recruit_step(
    fa: FocalAdhesion,
    dt: float,
    params: VinculinParams = DEFAULT_PARAMS,
    rng: np.random.Generator | None = None,
) -> None:
    """Advance ``fa.vinculin_count`` by ``dt`` under Poisson tau-leap.

    Per timestep:

      - recruitment events ~ Poisson(``k_rec · n_unfolded · N_free · dt``)
      - dissociation events ~ Poisson(``k_diss · N_vin · dt``)

    The two are independent (tau-leap on uncoupled birth / death).
    Expected dynamics reproduce the continuum ODE
    ``dN_vin/dt = k_rec n N_free − k_diss N_vin`` with steady state
    ``N_vin* = k_rec n N_free / k_diss``. ``rng`` is required (the
    earlier "deterministic explicit Euler" variant rounded per-step Δ
    below 0.5 to zero and left N_vin pinned at 0 — regression test
    ``test_vinculin_recruit_steady_state``).
    """
    if rng is None:
        rng = np.random.default_rng()
    n_unfolded = fa.talin_unfolded_domains
    N_vin = int(fa.vinculin_count)
    rate_in = params.k_rec * n_unfolded * params.N_free
    rate_out = params.k_diss * N_vin
    n_in = int(rng.poisson(rate_in * dt)) if rate_in > 0.0 else 0
    n_out = int(rng.poisson(rate_out * dt)) if rate_out > 0.0 else 0
    fa.vinculin_count = max(0, N_vin + n_in - n_out)


def effective_clutch_stiffness(
    fa: FocalAdhesion,
    bare_k_int: float,
    params: VinculinParams = DEFAULT_PARAMS,
) -> float:
    """``k_int^eff = k_int^bare (1 + α N_vin)``  [N/m]  (KU-2.7 / Han 2021)."""
    return bare_k_int * (1.0 + params.alpha * float(fa.vinculin_count))
