"""Chan-Odde motor-clutch dynamics on a single FA (KU-2.4, KU-2.18).

This module implements the canonical Chan & Odde 2008 model with the
Pereverzev catch-slip extension (Bangasser 2013, Bangasser 2017):

* Each focal adhesion holds ``N_total`` integrin clutches; at any time a
  subset is engaged (bound to the substrate).
* Engaged clutches anchor to the substrate at the position the actin
  bundle had when the clutch bound (lab frame). While engaged the
  clutch tension is::

      F_i = k_int · (x_actin − x_anchor_i)         [N]

  where ``k_int`` is the integrin spring stiffness (KU-2.18). All
  engaged clutches share a *common* substrate-side displacement
  ``x_sub`` that solves the quasi-static force balance::

      k_sub · x_sub = Σ k_int · (x_actin − a_i − x_sub)
      ⇒ x_sub = k_int (N_eng x_actin − Σ a_i) / (k_sub + N_eng k_int)

  i.e. the clutches and substrate act as a parallel-spring system. This
  is the closed-form solution one gets when the substrate inertia is
  negligible compared to bond on/off rates — true for any soft-tissue
  PAA gel.
* Myosin retrograde flow loads the FA::

      v_actin = v_unloaded · max(0, 1 − F_total / (N_motors · F_stall))
      x_actin ← x_actin + v_actin · dt

* Each engaged clutch breaks at rate ``k_off(F_i)`` from the Pereverzev
  form (KU-2.5); each disengaged clutch rebinds at rate ``k_on``
  (KU-2.18).  Stochastic decision per timestep: ``p_event = 1 −
  exp(−k · dt)``.

Sanity Gate
-----------
1. Dimensional analysis: ``k_int``, ``k_sub`` in N/m; ``x``, ``a`` in
   m; ``F`` in N; ``v`` in m/s; ``k_on``, ``k_off`` in s⁻¹. ``dt`` must
   be ≪ the bond-event timescale (``1 / max(k_on, k_off(F))``) for the
   ``1 − exp(−k·dt)`` probability to be valid; the Sanity Gate checks
   this against the slip blow-up at ``F = F_s``.
2. Boundary cases: ``N_eng = 0`` → ``x_sub = 0`` (no clutches, no
   substrate load); all clutches engaged → force-balance still well
   defined; ``v_actin`` clamped at zero (stall reversal is not
   physical in this model).
3. Conservation: stochastic dynamics, no conservation invariant beyond
   probability normalisation (each clutch is either engaged or not,
   never both).
4. Numerical sanity: float64; vectorised over clutches; ``rng`` is
   ``numpy.random.Generator`` and is the *single* source of stochastic
   decisions (seed-controlled in tests).
5. Sign / sense: actin retrograde flow is the +x direction; engaged
   clutches resist the flow (positive ``F_i``); the substrate displaces
   in the same direction by the smaller amount ``x_sub < x_actin``.
6. Measurement protocol: ``step`` returns a diagnostics ``dict`` of
   ``(force_total, v_actin, n_engaged)`` measured *after* the
   stochastic update of this timestep; the biphasic regression test
   averages ``force_total`` over the second half of the time series to
   discard the engagement transient.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from acs_kb.bridge.catch_bond import (
    DEFAULT_PARAMS as DEFAULT_CATCH_PARAMS,
    CatchSlipParams,
    catch_slip_off_rate,
)
from acs_kb.bridge.substrate_stub import LinearElasticSubstrate
from acs_kb.bridge.types import FocalAdhesion


@dataclass(frozen=True, slots=True)
class MotorClutchParams:
    """KU-2.4 / KU-2.18 defaults for a single FA."""

    n_clutches: int = 50
    n_motors: int = 50
    k_on: float = 1.0                  # s⁻¹
    k_int: float = 1.0e-3              # N/m   (1 pN/nm, KU-2.18)
    v_unloaded: float = 100.0e-9       # m/s   (100 nm/s)
    F_stall_per_motor: float = 2.0e-12 # N     (2 pN, KU-2.18)
    bond: CatchSlipParams = DEFAULT_CATCH_PARAMS


class MotorClutchFA:
    """Single-FA Chan-Odde motor-clutch stepper (KU-2.4).

    Holds a :class:`FocalAdhesion` plus the simulator-side state
    (per-clutch substrate-frame anchors and the common ``x_sub``) that
    is not part of the public week-5 contract.

    Parameters
    ----------
    fa : FocalAdhesion
        The FA being driven. Mutated in place by :meth:`step`.
    substrate : LinearElasticSubstrate
        Phase 1 stub; will become an ECM adapter in Unit 2.2+.
    params : MotorClutchParams, optional
        Override KU-2.18 defaults if needed.
    """

    __slots__ = ("fa", "substrate", "params", "_anchors", "_x_sub", "_total_force")

    def __init__(
        self,
        fa: FocalAdhesion,
        substrate: LinearElasticSubstrate,
        params: MotorClutchParams | None = None,
    ) -> None:
        self.fa = fa
        self.substrate = substrate
        self.params = params if params is not None else MotorClutchParams()
        # Per-clutch anchor on the substrate (substrate frame). ``nan``
        # means disengaged; resets on bind.
        self._anchors = np.full(fa.n_clutches_total, np.nan, dtype=np.float64)
        self._x_sub = 0.0
        self._total_force = 0.0

    # ------------------------------------------------------------------ #
    # Quasi-static force balance                                         #
    # ------------------------------------------------------------------ #

    def _solve_force_balance(self) -> tuple[np.ndarray, float, float]:
        """Compute (per-clutch forces, x_sub, F_total) at the current state.

        Solves the parallel-spring force balance::

            x_sub = k_int (N_eng x_actin − Σ a_i) / (k_sub + N_eng k_int)
            F_i   = k_int (x_actin − a_i − x_sub)

        Returns
        -------
        F : np.ndarray, shape (n_clutches_total,)
            Force per clutch; 0 for disengaged.
        x_sub : float
            Substrate displacement under the current engagement state.
        F_total : float
            Σ F_i, equal to ``k_sub · x_sub`` (force balance).
        """
        eng = self.fa.clutches_engaged
        n_eng = int(eng.sum())
        F = np.zeros(self.fa.n_clutches_total, dtype=np.float64)
        if n_eng == 0:
            return F, 0.0, 0.0

        k_int = self.params.k_int
        k_sub = self.substrate.stiffness
        x_actin = self.fa.actin_position
        anchors_eng = self._anchors[eng]
        sum_a = float(anchors_eng.sum())
        denom = k_sub + n_eng * k_int
        x_sub = k_int * (n_eng * x_actin - sum_a) / denom
        F_eng = k_int * (x_actin - anchors_eng - x_sub)
        F[eng] = F_eng
        return F, float(x_sub), float(F_eng.sum())

    # ------------------------------------------------------------------ #
    # Main step                                                          #
    # ------------------------------------------------------------------ #

    def step(self, dt: float, rng: np.random.Generator) -> dict[str, Any]:
        """Advance the FA by ``dt`` seconds.

        The ordering is:
          1. Force balance at the *current* engagement state.
          2. Advance actin position with the motor-stall velocity
             computed from the current ``F_total``.
          3. Stochastic bind / unbind events (per-clutch).
          4. Re-solve force balance and write back to the FA dataclass.

        Returns
        -------
        diagnostics : dict
            ``force_total`` (N), ``v_actin`` (m/s), ``n_engaged`` (int),
            ``x_sub`` (m). Measured *after* the stochastic update of
            this timestep.
        """
        p = self.params
        fa = self.fa

        # 1. Force balance at current engagement state.
        F_now, _, F_total_now = self._solve_force_balance()

        # 2. Motor-stall actin advance. Clamp at stall (no reverse flow).
        F_max = p.n_motors * p.F_stall_per_motor
        v_actin = p.v_unloaded * max(0.0, 1.0 - F_total_now / F_max)
        fa.actin_position += v_actin * dt
        fa.age += dt

        # 3. Stochastic bind / unbind.
        eng = fa.clutches_engaged
        # 3a. Disengage: each engaged clutch fails with p = 1 − exp(−k_off·dt).
        if eng.any():
            k_off = catch_slip_off_rate(F_now[eng], p.bond)
            p_break = 1.0 - np.exp(-k_off * dt)
            break_draws = rng.random(p_break.shape[0])
            broke = break_draws < p_break
            if broke.any():
                # Map back to global clutch indices.
                eng_idx = np.flatnonzero(eng)
                broke_idx = eng_idx[broke]
                fa.clutches_engaged[broke_idx] = False
                self._anchors[broke_idx] = np.nan
        # 3b. Engage: each disengaged clutch binds with p = 1 − exp(−k_on·dt).
        free = ~fa.clutches_engaged
        if free.any():
            p_bind = 1.0 - np.exp(-p.k_on * dt)
            bind_draws = rng.random(int(free.sum()))
            binds = bind_draws < p_bind
            if binds.any():
                free_idx = np.flatnonzero(free)
                bind_idx = free_idx[binds]
                fa.clutches_engaged[bind_idx] = True
                # Anchor in the substrate frame: clutch attaches at the
                # current actin position minus the current substrate
                # offset, so its instantaneous force is zero.
                self._anchors[bind_idx] = fa.actin_position - self._x_sub

        # 4. Re-solve force balance under the new engagement state.
        F_new, x_sub_new, F_total_new = self._solve_force_balance()
        fa.clutch_forces[:] = F_new
        self._x_sub = x_sub_new
        self._total_force = F_total_new

        return {
            "force_total": F_total_new,
            "v_actin": v_actin,
            "n_engaged": int(fa.clutches_engaged.sum()),
            "x_sub": x_sub_new,
        }

    # ------------------------------------------------------------------ #
    # Diagnostic helpers                                                 #
    # ------------------------------------------------------------------ #

    @property
    def x_sub(self) -> float:
        """Last-solved substrate displacement (m)."""
        return self._x_sub

    @property
    def total_force(self) -> float:
        """Last-solved Σ F_i over engaged clutches (N)."""
        return self._total_force


def run_steady_state(
    mc: MotorClutchFA,
    *,
    dt: float,
    n_steps: int,
    rng: np.random.Generator,
    burn_in_fraction: float = 0.5,
) -> dict[str, np.ndarray | float]:
    """Run a single FA for ``n_steps`` and report steady-state averages.

    ``burn_in_fraction`` of the leading samples is discarded so the
    transient engagement (starting with zero engaged clutches) does not
    contaminate the mean. Returns a dict of arrays of length
    ``n_steps - burn_in`` plus the means.
    """
    f_total = np.empty(n_steps, dtype=np.float64)
    v = np.empty(n_steps, dtype=np.float64)
    n_eng = np.empty(n_steps, dtype=np.int32)
    for i in range(n_steps):
        d = mc.step(dt, rng)
        f_total[i] = d["force_total"]
        v[i] = d["v_actin"]
        n_eng[i] = d["n_engaged"]
    burn_in = int(burn_in_fraction * n_steps)
    return {
        "force_total": f_total,
        "v_actin": v,
        "n_engaged": n_eng,
        "mean_force_total": float(f_total[burn_in:].mean()),
        "mean_v_actin": float(v[burn_in:].mean()),
        "mean_n_engaged": float(n_eng[burn_in:].mean()),
    }
