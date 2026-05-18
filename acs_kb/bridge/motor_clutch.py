"""Chan-Odde motor-clutch dynamics on a single FA (KU-2.4, KU-2.18).

This module implements the canonical Chan & Odde 2008 model with the
Pereverzev catch-slip extension (Bangasser 2013, Bangasser 2017) and,
from Unit 2.2 onward, the talin / vinculin / FA-growth maturation
machinery (KU-2.6 / KU-2.7 / KU-2.17). The maturation substeps are
opt-in via constructor parameters; passing ``None`` for any of them
falls back to the Unit 2.1 nascent-FA behaviour (which is what most of
the unit tests exercise).

* Each focal adhesion holds ``N_total`` integrin clutches; at any time a
  subset is engaged (bound to the substrate).
* Engaged clutches anchor to the substrate at the position the actin
  bundle had when the clutch bound (lab frame). While engaged the
  clutch tension is::

      F_i = k_int^eff · (x_actin − x_anchor_i)         [N]

  with ``k_int^eff = k_int^bare · (1 + α · N_vin)`` (KU-2.7 vinculin
  allostery); in Unit 2.1 there is no vinculin and the effective
  stiffness collapses to ``k_int^bare``.
* All engaged clutches share a *common* substrate-side displacement
  ``x_sub`` that solves the quasi-static force balance::

      k_sub · x_sub = Σ k_int^eff · (x_actin − a_i − x_sub)

* Myosin retrograde flow loads the FA::

      v_actin = v_unloaded · max(0, 1 − F_total / (N_motors · F_stall))

* Bond off-rate from Pereverzev catch-slip (KU-2.5).
* Substrate object exposes ``.stiffness`` and ``.compute_displacement``
  (a :class:`acs_kb.bridge.ecm_adapter.SubstrateProtocol`). The stub
  and the ECM adapter are interchangeable.

Sanity Gate
-----------
Inherits Unit 2.1's six checks; Unit 2.2 adds:

7. Talin substep is gated by ``talin_params is not None`` so Unit 2.1
   tests continue to pass; when active, the unfold rate uses the *mean
   per-engaged-clutch* force, matching KU-2.6's per-domain Bell-Evans.
8. Vinculin substep deterministically integrates the recruitment ODE
   with explicit Euler (timescale 1/k_rec ≈ 10 s ≫ dt).
9. FA growth substep resizes ``fa.clutches_engaged`` and
   ``fa.clutch_forces`` in place; ``self._anchors`` is padded/trimmed
   to keep its size in sync with ``fa.n_clutches_total``.
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
from acs_kb.bridge.fa_growth import (
    FAGrowthParams,
    fa_growth_step,
)
from acs_kb.bridge.substrate_stub import LinearElasticSubstrate
from acs_kb.bridge.talin import (
    TalinParams,
    talin_unfold_step,
)
from acs_kb.bridge.types import FocalAdhesion
from acs_kb.bridge.vinculin import (
    VinculinParams,
    effective_clutch_stiffness,
    vinculin_recruit_step,
)


@dataclass(frozen=True, slots=True)
class MotorClutchParams:
    """KU-2.4 / KU-2.18 motor-clutch parameters for a single FA.

    Defaults mirror KU-2.18 so the dataclass is usable standalone, but
    the *authoritative* source for a simulation run is
    ``configs/phase1_unit2_{1,2}.yaml`` via :meth:`from_config`.
    """

    n_clutches: int = 50
    n_motors: int = 50
    k_on: float = 1.0                  # s⁻¹
    k_int: float = 1.0e-3              # N/m   (1 pN/nm, KU-2.18)
    v_unloaded: float = 100.0e-9       # m/s   (100 nm/s)
    F_stall_per_motor: float = 2.0e-12 # N     (2 pN, KU-2.18)
    bond: CatchSlipParams = DEFAULT_CATCH_PARAMS

    @classmethod
    def from_config(cls, cfg: dict) -> "MotorClutchParams":
        b = cfg["bridge"] if "bridge" in cfg else cfg
        m = b["motor_clutch"]
        return cls(
            n_clutches=int(m["n_clutches"]),
            n_motors=int(m["n_motors"]),
            k_on=float(m["k_on"]),
            k_int=float(m["k_int"]),
            v_unloaded=float(m["v_unloaded"]),
            F_stall_per_motor=float(m["F_stall_per_motor"]),
            bond=CatchSlipParams.from_config({"bridge": b}),
        )


class MotorClutchFA:
    """Single-FA Chan-Odde motor-clutch stepper (KU-2.4 + KU-2.6/2.7/2.17).

    Parameters
    ----------
    fa : FocalAdhesion
        The FA being driven. Mutated in place by :meth:`step`.
    substrate : LinearElasticSubstrate or ECMAdapter
        Anything that satisfies the ``SubstrateProtocol`` (``.stiffness``
        + ``.compute_displacement``). Unit 2.1 uses the stub; Unit 2.2's
        regression suite verifies the ECM adapter is interchangeable.
    params : MotorClutchParams, optional
        Override KU-2.18 defaults if needed.
    talin_params, vinculin_params, fa_growth_params : optional
        Phase-1 Unit-2.2 maturation params. Pass ``None`` (default) to
        skip the corresponding substep — keeps Unit 2.1 behaviour intact.
    """

    __slots__ = (
        "fa", "substrate", "params",
        "talin_params", "vinculin_params", "fa_growth_params",
        "_anchors", "_x_sub", "_total_force",
    )

    def __init__(
        self,
        fa: FocalAdhesion,
        substrate: LinearElasticSubstrate,
        params: MotorClutchParams | None = None,
        *,
        talin_params: TalinParams | None = None,
        vinculin_params: VinculinParams | None = None,
        fa_growth_params: FAGrowthParams | None = None,
    ) -> None:
        self.fa = fa
        self.substrate = substrate
        self.params = params if params is not None else MotorClutchParams()
        self.talin_params = talin_params
        self.vinculin_params = vinculin_params
        self.fa_growth_params = fa_growth_params
        self._anchors = np.full(fa.n_clutches_total, np.nan, dtype=np.float64)
        self._x_sub = 0.0
        self._total_force = 0.0

    # ------------------------------------------------------------------ #
    # Quasi-static force balance                                         #
    # ------------------------------------------------------------------ #

    def _effective_k_int(self) -> float:
        """``k_int^eff = k_int^bare (1 + α N_vin)`` (KU-2.7)."""
        if self.vinculin_params is None:
            return self.params.k_int
        return effective_clutch_stiffness(self.fa, self.params.k_int,
                                          self.vinculin_params)

    def _solve_force_balance(
        self, k_int: float | None = None
    ) -> tuple[np.ndarray, float, float]:
        """Compute (per-clutch forces, x_sub, F_total) under ``k_int^eff``.

        ``k_int`` defaults to the current vinculin-modulated stiffness.
        Returns zeros if no clutches are engaged.
        """
        if k_int is None:
            k_int = self._effective_k_int()
        eng = self.fa.clutches_engaged
        n_eng = int(eng.sum())
        F = np.zeros(self.fa.n_clutches_total, dtype=np.float64)
        if n_eng == 0:
            return F, 0.0, 0.0
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
    # Anchor housekeeping for FA growth                                  #
    # ------------------------------------------------------------------ #

    def _resync_anchors_to_fa_size(self) -> None:
        """Pad / trim ``self._anchors`` to match ``fa.n_clutches_total``.

        FA growth (KU-2.17) resizes the clutch arrays in place. Newly
        added clutches are disengaged (anchor = NaN); shrinking drops
        disengaged clutches first, so any trimmed anchors are already
        NaN — losing them is safe.
        """
        n_target = self.fa.n_clutches_total
        n_have = self._anchors.size
        if n_target == n_have:
            return
        if n_target > n_have:
            extend = n_target - n_have
            self._anchors = np.concatenate([
                self._anchors, np.full(extend, np.nan, dtype=np.float64)
            ])
        else:
            self._anchors = self._anchors[:n_target].copy()

    # ------------------------------------------------------------------ #
    # Main step                                                          #
    # ------------------------------------------------------------------ #

    def step(self, dt: float, rng: np.random.Generator) -> dict[str, Any]:
        """Advance the FA by ``dt`` seconds.

        Ordering (Unit 2.2):
          1. Force balance at the current engagement + vinculin state.
          2. Motor-stall actin advance using the current F_total.
          3. Stochastic clutch bind / unbind using F from step 1.
          4. Talin unfold step (KU-2.6, opt-in).
          5. Vinculin recruit step (KU-2.7, opt-in).
          6. FA growth step (KU-2.17, opt-in) — may resize clutch arrays.
          7. Re-solve force balance with the updated state and write back.

        Returns a diagnostics dict with motor-clutch + maturation observables.
        """
        p = self.params
        fa = self.fa

        # 1. Force balance at current engagement state.
        F_now, _, F_total_now = self._solve_force_balance()

        # 2. Motor-stall actin advance.
        F_max = p.n_motors * p.F_stall_per_motor
        v_actin = p.v_unloaded * max(0.0, 1.0 - F_total_now / F_max)
        fa.actin_position += v_actin * dt
        fa.age += dt

        # 3. Stochastic bind / unbind.
        eng = fa.clutches_engaged
        if eng.any():
            k_off = catch_slip_off_rate(F_now[eng], p.bond)
            p_break = 1.0 - np.exp(-k_off * dt)
            break_draws = rng.random(p_break.shape[0])
            broke = break_draws < p_break
            if broke.any():
                eng_idx = np.flatnonzero(eng)
                broke_idx = eng_idx[broke]
                fa.clutches_engaged[broke_idx] = False
                self._anchors[broke_idx] = np.nan
        free = ~fa.clutches_engaged
        if free.any():
            p_bind = 1.0 - np.exp(-p.k_on * dt)
            bind_draws = rng.random(int(free.sum()))
            binds = bind_draws < p_bind
            if binds.any():
                free_idx = np.flatnonzero(free)
                bind_idx = free_idx[binds]
                fa.clutches_engaged[bind_idx] = True
                self._anchors[bind_idx] = fa.actin_position - self._x_sub

        # 4. Talin unfold (KU-2.6) — Phase 1 1-state, per-FA force = mean
        #    over engaged clutches' last forces.
        if self.talin_params is not None:
            if eng.any():
                eng_last_F = F_now[eng]
                f_per_talin = float(eng_last_F.mean()) if eng_last_F.size else 0.0
            else:
                f_per_talin = 0.0
            talin_unfold_step(fa, f_per_talin, dt, rng, self.talin_params)

        # 5. Vinculin recruit (KU-2.7).
        if self.vinculin_params is not None:
            vinculin_recruit_step(fa, dt, self.vinculin_params, rng)

        # 6. FA growth (KU-2.17) — uses pre-bind/unbind F_total_now.
        if self.fa_growth_params is not None:
            fa_growth_step(fa, F_total_now, dt, self.fa_growth_params)
            self._resync_anchors_to_fa_size()

        # 7. Re-solve force balance under the new engagement + vinculin state.
        F_new, x_sub_new, F_total_new = self._solve_force_balance()
        fa.clutch_forces[:] = F_new
        self._x_sub = x_sub_new
        self._total_force = F_total_new

        return {
            "force_total": F_total_new,
            "v_actin": v_actin,
            "n_engaged": int(fa.clutches_engaged.sum()),
            "x_sub": x_sub_new,
            "talin_unfolded": int(fa.talin_unfolded_domains),
            "vinculin_count": int(fa.vinculin_count),
            "fa_area": float(fa.area),
            "n_clutches_total": int(fa.n_clutches_total),
            "k_int_eff": self._effective_k_int(),
        }

    # ------------------------------------------------------------------ #
    # Diagnostic helpers                                                 #
    # ------------------------------------------------------------------ #

    @property
    def x_sub(self) -> float:
        return self._x_sub

    @property
    def total_force(self) -> float:
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

    Records both Unit 2.1 observables (force_total, v_actin, n_engaged)
    and Unit 2.2 maturation observables (vinculin_count, talin_unfolded,
    fa_area, n_clutches_total). The means are taken over the second
    half by default to discard the engagement / recruitment transient.
    """
    f_total = np.empty(n_steps, dtype=np.float64)
    v = np.empty(n_steps, dtype=np.float64)
    n_eng = np.empty(n_steps, dtype=np.int32)
    n_vin = np.empty(n_steps, dtype=np.int32)
    n_talin = np.empty(n_steps, dtype=np.int32)
    area = np.empty(n_steps, dtype=np.float64)
    n_clutch = np.empty(n_steps, dtype=np.int32)
    for i in range(n_steps):
        d = mc.step(dt, rng)
        f_total[i] = d["force_total"]
        v[i] = d["v_actin"]
        n_eng[i] = d["n_engaged"]
        n_vin[i] = d["vinculin_count"]
        n_talin[i] = d["talin_unfolded"]
        area[i] = d["fa_area"]
        n_clutch[i] = d["n_clutches_total"]
    burn_in = int(burn_in_fraction * n_steps)
    return {
        "force_total": f_total,
        "v_actin": v,
        "n_engaged": n_eng,
        "vinculin_count": n_vin,
        "talin_unfolded": n_talin,
        "fa_area": area,
        "n_clutches_total": n_clutch,
        "mean_force_total": float(f_total[burn_in:].mean()),
        "mean_v_actin": float(v[burn_in:].mean()),
        "mean_n_engaged": float(n_eng[burn_in:].mean()),
        "mean_vinculin_count": float(n_vin[burn_in:].mean()),
        "mean_fa_area": float(area[burn_in:].mean()),
        "mean_n_clutches_total": float(n_clutch[burn_in:].mean()),
    }
