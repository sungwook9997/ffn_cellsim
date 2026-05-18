"""Hill-function focal-adhesion growth (KU-2.17, Walcott & Sun 2010 PNAS).

FA area grows under sustained tension above a threshold and disassembles
otherwise. The continuum approximation of the four-state maturation
(KU-2.2) collapses into a single ODE::

    dA/dt = k_g(F) · A − k_d · A
    k_g(F) = k_g^0 · F^n / (F^n + F_th^n)            (Hill, exponent n)

Below ``F_th`` the growth rate vanishes ⇒ net decay at rate ``k_d``.
Above ``F_th`` the growth rate saturates at ``k_g^0`` ⇒ net rate
``k_g^0 − k_d`` (positive ⇒ growth, KU-2.17 defaults make this 0.05 s⁻¹).

Phase 1 maps the growing area back to a proportional clutch count::

    n_clutches_total = round( (A / A_nascent) · N_clutches_nascent )

with ``A_nascent`` the FA's initial area and ``N_clutches_nascent`` its
initial clutch count (KU-2.18 defaults: A = 1 μm², N = 50). When new
clutches are added their ``clutches_engaged`` slots are False and
``clutch_forces`` are 0 (Unit 2.1 contract is preserved).

Sanity Gate
-----------
1. Dimensional analysis: ``A`` [m²], ``F`` [N], ``k_g, k_d`` [s⁻¹],
   ``F_th`` [N]. ``F^n / (F^n + F_th^n)`` is dimensionless. No CFL —
   1/k_g ≈ 10 s ≫ dt = 1e-4 s, explicit Euler safe.
2. Boundary cases: ``F = 0`` ⇒ k_g = 0 ⇒ pure decay. ``A → 0`` ⇒ the
   FA dissolves and ``n_clutches_total → 0``; the helper guards
   against ``n_clutches_total = 0`` by clamping to a minimum of 1 so
   the dataclass arrays stay well-defined (the FA can still re-grow).
3. Conservation: not a conservation law — phenomenological maturation
   dynamics. The clutch count and area are kept consistent.
4. Numerical sanity: float64; pure NumPy; explicit Euler is fine for
   the slow timescale.
5. Sign / sense: positive ``F`` above ``F_th`` ⇒ positive ``dA/dt`` ⇒
   FA grows; below ``F_th`` ⇒ shrinks. ``k_d`` is strictly positive
   (FA cannot grow with no force).
6. Measurement protocol: the regression test runs 100 FAs in parallel
   for 5 min sim time and checks the resulting size distribution has
   three modes (NA / FC / FA) at ≈ 0.1 / 1 / 2+ μm (KU-2.2 / KU-2.17).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from acs_kb.bridge.types import FocalAdhesion


@dataclass(frozen=True, slots=True)
class FAGrowthParams:
    """KU-2.17 Hill-function parameters (Walcott & Sun 2010)."""

    k_g0: float = 0.1            # s⁻¹  saturating growth rate
    F_th: float = 1.0e-9         # N    1 nN per-FA threshold (KU-2.17)
    n_hill: float = 2.0          # Hill coefficient (KU-2.17, n ≈ 2–4)
    k_d: float = 0.05            # s⁻¹  decay rate
    A_nascent: float = 1.0e-12   # m²   nascent FA area (KU-2.18 / KU-2.2)
    n_clutches_nascent: int = 50 # clutches at A_nascent

    @classmethod
    def from_config(cls, cfg: dict) -> "FAGrowthParams":
        b = cfg["bridge"] if "bridge" in cfg else cfg
        g = b["fa_growth"]
        return cls(
            k_g0=float(g["k_g0"]),
            F_th=float(g["F_th"]),
            n_hill=float(g["n_hill"]),
            k_d=float(g["k_d"]),
            A_nascent=float(g["A_nascent"]),
            n_clutches_nascent=int(g["n_clutches_nascent"]),
        )


DEFAULT_PARAMS = FAGrowthParams()


def hill_growth_rate(force: float, params: FAGrowthParams = DEFAULT_PARAMS) -> float:
    """``k_g(F) = k_g^0 F^n / (F^n + F_th^n)``  [s⁻¹]."""
    if force <= 0.0:
        return 0.0
    Fn = force ** params.n_hill
    Fth_n = params.F_th ** params.n_hill
    return params.k_g0 * Fn / (Fn + Fth_n)


def _resize_focal_adhesion(fa: FocalAdhesion, n_new: int) -> None:
    """Resize the clutch arrays in-place to match a new total clutch count.

    Newly added clutches start disengaged with zero force; removed
    clutches are taken from the disengaged pool first, then from the
    engaged pool if necessary (their loss is recorded as a forced
    disengage; talin / vinculin are not yet per-clutch state in
    Phase 1 so no further bookkeeping is required).
    """
    n_old = fa.n_clutches_total
    if n_new == n_old:
        return
    n_new = max(1, n_new)
    if n_new > n_old:
        # Append disengaged slots.
        extend = n_new - n_old
        fa.clutches_engaged = np.concatenate([
            fa.clutches_engaged, np.zeros(extend, dtype=bool)
        ])
        fa.clutch_forces = np.concatenate([
            fa.clutch_forces, np.zeros(extend, dtype=np.float64)
        ])
    else:
        # Shrink: drop disengaged clutches first.
        keep = n_new
        engaged_idx = np.flatnonzero(fa.clutches_engaged)
        free_idx = np.flatnonzero(~fa.clutches_engaged)
        n_eng = engaged_idx.size
        if n_eng >= keep:
            # Keep the first `keep` engaged clutches; everything else dropped.
            keep_idx = engaged_idx[:keep]
        else:
            keep_idx = np.concatenate([engaged_idx, free_idx[:keep - n_eng]])
        keep_idx = np.sort(keep_idx)
        fa.clutches_engaged = fa.clutches_engaged[keep_idx].copy()
        fa.clutch_forces = fa.clutch_forces[keep_idx].copy()
    fa.n_clutches_total = n_new


def fa_growth_step(
    fa: FocalAdhesion,
    total_force: float,
    dt: float,
    params: FAGrowthParams = DEFAULT_PARAMS,
) -> None:
    """Advance ``fa.area`` (and the clutch count) by one Hill-growth step.

    Parameters
    ----------
    fa : FocalAdhesion
        Mutated in place: ``area`` and ``n_clutches_total`` (with the
        backing arrays) are updated.
    total_force : float
        Σ of engaged-clutch forces on the FA [N]. The Hill response is
        a per-FA function of *total* traction, KU-2.17.
    dt : float
        Timestep [s].
    """
    kg = hill_growth_rate(total_force, params)
    dA = (kg - params.k_d) * fa.area * dt
    new_area = max(0.1 * params.A_nascent, fa.area + dA)
    fa.area = new_area
    # Clutch count scales with area (KU-2.2 / KU-2.18).
    n_target = int(round(params.n_clutches_nascent * new_area / params.A_nascent))
    if n_target != fa.n_clutches_total:
        _resize_focal_adhesion(fa, n_target)
