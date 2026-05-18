"""E-cadherin trans-bond Bell-Evans slip-only model (KU-4.2, KU-4.17).

Phase 1 Unit 4.1 ships a **slip-only** Bell-Evans off-rate for the
cadherin trans-dimer; the full catch-bond detail (KU-4.2, Buckley 2014
Science) is reserved for Phase 2 once the slip-only baseline is
validated against the Maître two-cell rounding-angle gate.

The off-rate is

.. math::
    k_{\\rm off}(F) = k_{\\rm off}^0 \\, \\exp\\!\\bigl(F\\,\\Delta x^* / k_B T\\bigr)

with KU-4.17 defaults

.. list-table::
   :header-rows: 1

   * - Parameter
     - Symbol
     - Default
   * - Zero-force off-rate
     - :math:`k_{\\rm off}^0`
     - ``0.5 s⁻¹``
   * - Slip distance
     - :math:`\\Delta x^*`
     - ``4 nm`` (Buckley 2014 Fig 4)
   * - Encounter on-rate
     - :math:`k_{\\rm on}`
     - ``1 s⁻¹``
   * - Bonds per contact
     - :math:`N_{\\rm cad}`
     - ``100`` (mature junction)
   * - Mean per-bond load
     - :math:`\\langle F_{\\rm bond}\\rangle`
     - ``30 pN`` (KU-4.17)
   * - Temperature
     - :math:`T`
     - ``310 K``, :math:`k_B T = 4.28\\times 10^{-21} J`

Sanity Gate
-----------

1. **Dimensional**: :math:`F\\,\\Delta x^*/k_BT` is dimensionless;
   :math:`k_{\\rm off}` in :math:`s^{-1}`; survival probability over
   :math:`\\Delta t` is :math:`p_{\\rm survive} = \\exp(-k_{\\rm off}\\Delta t)`,
   also dimensionless. No CFL bound — discrete stochastic update.
2. **Boundary cases**:
   - ``F = 0`` ⇒ ``k_off = k_off0`` (zero-force slip rate, KU-4.17 anchor).
   - ``F → ∞`` ⇒ ``k_off → ∞`` (slip bond ruptures); guarded against
     overflow by clamping the exponent to ``700`` so ``exp(·)`` stays in
     ``float64`` range.
   - ``Δt → 0`` ⇒ ``p_break, p_engage → 0``; well-defined.
   - ``Δt → large`` ⇒ ``p → 1``; ``1 - exp(-x)`` is monotone and bounded
     in :math:`[0, 1]`.
   - ``n_engaged = 0`` ⇒ per-bond force is set to ``0`` (no carrier);
     only the ``k_on`` arm fires that step.
3. **Conservation**: ``n_bonds_engaged ∈ [0, n_bonds_total]`` at every
   step by construction (``np.clip`` after the two-binomial update).
   ``bond_forces.shape == (n_bonds_engaged,)`` invariant is restored on
   exit.
4. **Numerical sanity**: ``float64`` throughout; stochastic transitions
   use ``rng.binomial`` (exact for independent two-state events). No
   accumulating roundoff because the state is integer-valued.
5. **Sign / sense**: ``F ≥ 0`` (tensile axial load on a bond pulled to
   separate the two cells); the exponent monotonically increases
   :math:`k_{\\rm off}`, i.e. higher tension → faster break — the
   defining feature of a slip bond. ``bond_forces`` entries are
   non-negative.

References: Bell 1978 *Science*; Evans & Ritchie 1997 *Biophys J*;
Buckley *et al.* 2014 *Science* (cadherin slip vs. catch bond — slip
dominant arm at low force is used here).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from acs_kb.junction.types import EcadherinJunction

# KU-4.17 Phase 1 Unit 4.1 defaults. All quantities SI.
KU417_DEFAULTS: dict[str, float] = {
    "k_off0": 0.5,             # s^-1   zero-force off-rate
    "dx_star": 4.0e-9,         # m      Bell slip distance (Buckley 2014)
    "k_on": 1.0,               # s^-1   encounter on-rate
    "N_cad_per_contact": 100,  # bonds  KU-4.17 mature junction
    "F_per_bond_average": 30.0e-12,  # N  mean per-bond axial load
    "kT": 4.28e-21,            # J      k_B T at 310 K
}

# Numerical guard: exp(700) ≈ 1e304, the float64 ceiling. Force values
# above this exponent would overflow, so we clamp the dimensionless
# Bell exponent F·Δx*/kT to the same bound. At KU-4.17 defaults this
# corresponds to F ≈ 700 · kT / dx_star ≈ 7.5 × 10⁻¹⁰ N = 0.75 nN per
# *single* bond, well above any physically meaningful per-bond load
# (the bond ruptures long before that — k_off would be 10³⁰⁰ s⁻¹).
_EXP_ARG_MAX: float = 700.0


def bond_off_rate(
    F: np.ndarray | float,
    *,
    k_off0: float = KU417_DEFAULTS["k_off0"],
    dx_star: float = KU417_DEFAULTS["dx_star"],
    kT: float = KU417_DEFAULTS["kT"],
) -> np.ndarray | float:
    """Bell-Evans slip-only off-rate :math:`k_{\\rm off}(F)`.

    Parameters
    ----------
    F : array-like or float
        Axial tensile force per bond (N), broadcastable.
    k_off0, dx_star, kT : float
        KU-4.17 / KU-4.2 constants; see module defaults.

    Returns
    -------
    k_off : same type as ``F``
        Off-rate in s⁻¹. Always positive.
    """
    arg = np.asarray(F, dtype=np.float64) * dx_star / kT
    arg = np.clip(arg, -_EXP_ARG_MAX, _EXP_ARG_MAX)
    return float(k_off0) * np.exp(arg)


def update_bonds(
    junction: EcadherinJunction,
    F_total: float,
    dt: float,
    rng: np.random.Generator,
    *,
    k_off0: float = KU417_DEFAULTS["k_off0"],
    dx_star: float = KU417_DEFAULTS["dx_star"],
    k_on: float = KU417_DEFAULTS["k_on"],
    kT: float = KU417_DEFAULTS["kT"],
) -> EcadherinJunction:
    """Advance a junction's bond population by one stochastic timestep.

    The bond population is treated as two compartments: ``n_engaged``
    bound bonds and ``n_total - n_engaged`` free bonds in the reservoir.
    Transitions in one step :math:`\\Delta t`:

    * Each engaged bond breaks with probability
      :math:`p_{\\rm break} = 1 - \\exp(-k_{\\rm off}(F_{\\rm bond}) \\Delta t)`,
      where :math:`F_{\\rm bond} = F_{\\rm total} / n_{\\rm engaged}` is
      the mean-field per-bond load (Phase 1 simplification; per-bond
      variance is a Phase 2 extension). When ``n_engaged = 0``, the
      per-bond load is ``0`` and only the on-rate arm fires.
    * Each free bond engages with probability
      :math:`p_{\\rm engage} = 1 - \\exp(-k_{\\rm on} \\Delta t)`.

    Both transition counts are drawn from a binomial distribution
    (``rng.binomial``), which is exact for independent events. The order
    of the two draws (break first, then engage) does not affect the
    result because we compute :math:`p_{\\rm break}` against the
    pre-step ``n_engaged`` and :math:`p_{\\rm engage}` against the
    pre-step free count, then combine — no double-counting is possible.

    On exit the junction's ``bond_forces`` array is rebuilt to length
    ``n_engaged_new`` with each entry equal to
    :math:`F_{\\rm total} / n_{\\rm engaged}^{\\rm new}` (or all zeros
    if ``n_engaged_new = 0``).

    Parameters
    ----------
    junction : EcadherinJunction
        Modified in place; returned for chained calls.
    F_total : float
        Externally applied axial force across the junction (N), ≥ 0.
    dt : float
        Timestep (s), > 0.
    rng : np.random.Generator
        Source of stochastic decisions; seed at the call site for
        reproducibility.
    k_off0, dx_star, k_on, kT : float
        KU-4.17 defaults; override for sensitivity sweeps.

    Returns
    -------
    EcadherinJunction
        Same instance, with ``n_bonds_engaged``, ``bond_forces``, and
        ``age`` updated.
    """
    if dt <= 0:
        raise ValueError(f"dt must be positive (got {dt}).")
    if F_total < 0:
        raise ValueError(f"F_total must be non-negative (got {F_total}).")
    n_total = int(junction.n_bonds_total)
    n_engaged = int(junction.n_bonds_engaged)
    if not (0 <= n_engaged <= n_total):
        raise ValueError(
            f"n_engaged={n_engaged} out of [0, n_total={n_total}]; "
            "junction state corrupted upstream."
        )

    if n_engaged > 0:
        F_per_bond = float(F_total) / float(n_engaged)
        k_off = bond_off_rate(F_per_bond, k_off0=k_off0, dx_star=dx_star, kT=kT)
        p_break = 1.0 - np.exp(-float(k_off) * float(dt))
    else:
        F_per_bond = 0.0
        p_break = 0.0

    p_engage = 1.0 - np.exp(-float(k_on) * float(dt))
    p_engage = float(np.clip(p_engage, 0.0, 1.0))
    p_break = float(np.clip(p_break, 0.0, 1.0))

    n_broken = int(rng.binomial(n_engaged, p_break)) if n_engaged > 0 else 0
    n_free = n_total - n_engaged
    n_new = int(rng.binomial(n_free, p_engage)) if n_free > 0 else 0

    n_engaged_new = int(np.clip(n_engaged - n_broken + n_new, 0, n_total))

    if n_engaged_new > 0:
        F_per_bond_new = float(F_total) / float(n_engaged_new)
        bond_forces = np.full(n_engaged_new, F_per_bond_new, dtype=np.float64)
    else:
        bond_forces = np.zeros(0, dtype=np.float64)

    junction.n_bonds_engaged = n_engaged_new
    junction.bond_forces = bond_forces
    junction.age = float(junction.age) + float(dt)
    return junction
