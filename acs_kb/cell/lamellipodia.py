"""Lamellipodium structural prep for Phase 1 Unit 3.2 (KU-3.6, KU-3.7, KU-3.12).

This module ships the **interface** that Unit 3.2 will fill in:

* :class:`Lamellipodium` dataclass — the per-protrusion state Worker
  C's Unit 3.2 work will own.
* :func:`brownian_ratchet_velocity` — the literal KU-3.6 force-velocity
  formula. Implemented now (it is a closed-form scalar equation) so
  the acceptance test can run; it has no time-stepping dependency on
  Worker B's clutch.
* :func:`apply_tip_load_from_clutch` — typed coupling to Worker B's
  motor-clutch output. Phase 1 Unit 3.1 ships a **mock** implementation
  that takes the per-filament load directly from a caller-supplied
  callable; Unit 3.2 will swap the callable for the real Worker B FA
  population sampler once the freeze in `acs_kb/bridge/types.py` is in
  place.
* :func:`advance_lamellipodium` — a single quasi-static advance step.
  Stubbed (returns the Lamellipodium unchanged with a NotImplementedError-
  raising hook) so Unit 3.2 can fill in the time integration without
  touching the dataclass shape.

The contract here is what Unit 3.2 will rely on; this file is
deliberately scope-frozen at "schema + KU-3.6 closed-form + mock
motor-clutch interface".

Sanity Gate
-----------
1. Dimensional analysis: [δ] = m, [k_on c_G] = 1/s, [k_off] = 1/s,
   [F δ / k_B T] = N·m / J = dimensionless. v_p has units m/s. ✓.
   No time integration in the formula → no CFL.
2. Boundary cases:
   - F → 0 ⇒ v_p → δ (k_on c_G − k_off) (free-protrusion limit).
   - F → ∞ ⇒ v_p → −δ k_off (depolymerisation-dominated).
   - F · δ ≫ k_B T ⇒ exp factor underflows safely (numpy).
   - δ ≤ 0, k_on ≤ 0, c_G ≤ 0 ⇒ ValueError.
3. Conservation: closed-form scalar — no state to conserve.
4. Numerical sanity: float64; exp underflow returns 0 naturally.
5. Sign / sense: positive v_p = protrusion (membrane forward).
   Negative v_p = retraction (depolymerisation outpacing).
6. Measurement protocol: the KU-3.6 VALIDATION value is the *free*
   protrusion v_p ≈ 30 nm/s (per the KU's literal target — see the
   inline note in :func:`brownian_ratchet_velocity` about the
   physiological c_G that yields this). The test exercises the
   formula at F = 0 against that value with ±20 % tolerance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Protocol, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    # Worker B's motor-clutch output. Stub Protocol for typing only —
    # Unit 3.2 will replace this with the concrete import from
    # acs_kb.bridge.types once Worker B freezes the FA dataclass.
    class _ClutchPopulationSampler(Protocol):
        def sample_tip_load(self, fa_id: int) -> float:
            """Return the instantaneous per-filament tip load (N) at this FA."""
            ...


# ---------------------------------------------------------------- #
# KU-3.6 closed-form polymerisation velocity                       #
# ---------------------------------------------------------------- #

# Boltzmann constant × T = 310 K = 4.28e-21 J  (KU-1.1, used here for
# the exponential load factor F δ / k_BT). Imported as a module-level
# constant rather than from acs_kb.common to keep this stub self-
# contained for Unit 3.2 swap-in.
_KBT_310K = 4.28e-21  # J  (k_B · 310 K — KU-1.1 anchor)

# Per-filament constants from KU-3.6 (Mogilner & Oster 1996 BioPhys J).
# c_G is left as a caller-supplied parameter because the free-protrusion
# VALIDATION number (≈ 30 nm/s) corresponds to a physiological local
# G-actin concentration that is lower than the KU's bulk 10 μM (the
# KU notes itself that the formula with bulk c_G gives ~300 nm/s; the
# "30 nm/s" target uses the leading-edge depleted concentration).
DELTA_ACTIN_MONOMER_M = 2.7e-9            # KU-3.6 actin monomer half-length δ [m]
K_ON_ACTIN_BARBED_PER_MICROMOLAR_S = 11.0  # KU-3.6 barbed-end on-rate [μM⁻¹·s⁻¹]
K_OFF_ACTIN_BARBED_PER_S = 1.0             # KU-3.6 off-rate [s⁻¹]


def brownian_ratchet_velocity(
    *,
    tip_force_N: float,
    g_actin_concentration_uM: float = 1.0,
    delta_m: float = DELTA_ACTIN_MONOMER_M,
    k_on_per_uM_s: float = K_ON_ACTIN_BARBED_PER_MICROMOLAR_S,
    k_off_per_s: float = K_OFF_ACTIN_BARBED_PER_S,
    kBT_J: float = _KBT_310K,
) -> float:
    """KU-3.6 Brownian-ratchet per-filament polymerisation velocity.

    .. math::

        v_p = \\delta \\Big[k_{on}^0 c_G \\exp\\!\\big(-\\tfrac{F \\delta}{k_B T}\\big) - k_{off}\\Big]

    All parameters carry SI / KU-3.6 units (δ in m, c_G in μM, k_on
    in μM⁻¹·s⁻¹). The default ``g_actin_concentration_uM = 1.0`` is
    the physiological depleted-leading-edge concentration that
    reproduces the KU-3.6 VALIDATION free-protrusion v_p ≈ 30 nm/s;
    the bulk 10 μM value the KU mentions in passing yields ~300 nm/s
    (see the docstring at the module level).

    Returns the per-filament velocity in m/s. Positive = protrusion.
    """
    if delta_m <= 0 or g_actin_concentration_uM <= 0:
        raise ValueError("delta_m and g_actin_concentration_uM must be > 0")
    if k_on_per_uM_s < 0 or k_off_per_s < 0 or kBT_J <= 0:
        raise ValueError("rates and kBT must be non-negative / positive")

    load_factor = np.exp(-tip_force_N * delta_m / kBT_J)
    on_rate = k_on_per_uM_s * g_actin_concentration_uM * load_factor
    return float(delta_m * (on_rate - k_off_per_s))


# ---------------------------------------------------------------- #
# Lamellipodium dataclass                                          #
# ---------------------------------------------------------------- #

@dataclass(slots=True)
class Lamellipodium:
    """Phase 1 Unit 3.2 lamellipodium state (scope-frozen by Unit 3.1).

    Attributes
    ----------
    id : int
        Stable identifier for the protrusion within a cell.
    cell_id : int
        Owning cell's :class:`acs_kb.cell.cell.Cell` id.
    base_position : np.ndarray, shape (2,)
        Cortex anchor in m (the bead on the cortex this lamellipodium
        protrudes from).
    tip_position : np.ndarray, shape (2,)
        Current tip position in m (the membrane-leading-edge point).
    polarity : np.ndarray, shape (2,)
        Outward unit normal that the lamellipodium protrudes along
        (KU-3.7 leading-edge direction). Phase 1 keeps this fixed
        per lamellipodium; Unit 3.3+ may update from cortex shape.
    n_filaments_at_tip : int
        Effective number of actin filaments at the leading edge
        (KU-3.7: ~100 per μm of leading-edge width × the
        ``width_um`` of this protrusion).
    width_um : float
        Tangential width of the lamellipodium in μm (KU-3.7 range
        5–20 μm; Phase 1 Unit 3.1 default 5).
    g_actin_uM : float
        Local G-actin concentration at the tip (μM). KU-3.6 VALIDATION
        anchor — default 1.0 reproduces 30 nm/s free protrusion.
    retrograde_flow_nm_per_s : float
        F-actin retrograde flow speed at the tip (KU-3.12). Set by
        Unit 3.2 via the motor-clutch coupling; default 0 = free
        protrusion.
    last_tip_force_N : float
        Snapshot of the most-recent per-filament tip load (N).
        Populated by :func:`apply_tip_load_from_clutch`. Defaults to 0.
    """

    id: int
    cell_id: int
    base_position: np.ndarray
    tip_position: np.ndarray
    polarity: np.ndarray
    n_filaments_at_tip: int
    width_um: float = 5.0
    g_actin_uM: float = 1.0
    retrograde_flow_nm_per_s: float = 0.0
    last_tip_force_N: float = 0.0
    extra: dict = field(default_factory=dict)

    @classmethod
    def at_cortex_bead(
        cls,
        *,
        lamellipodium_id: int,
        cell_id: int,
        bead_position: np.ndarray,
        outward_normal: np.ndarray,
        n_filaments_at_tip: int = 5,
        width_um: float = 5.0,
    ) -> "Lamellipodium":
        """Construct a lamellipodium anchored at a cortex bead.

        The tip starts coincident with the base (no protrusion yet);
        :func:`advance_lamellipodium` extends it.
        """
        base = np.asarray(bead_position, dtype=np.float64)
        normal = np.asarray(outward_normal, dtype=np.float64)
        n = np.linalg.norm(normal)
        if n <= 0.0:
            raise ValueError("outward_normal must be non-zero")
        return cls(
            id=int(lamellipodium_id),
            cell_id=int(cell_id),
            base_position=base,
            tip_position=base.copy(),
            polarity=normal / n,
            n_filaments_at_tip=int(n_filaments_at_tip),
            width_um=float(width_um),
        )

    def protrusion_length_m(self) -> float:
        """Current radial protrusion from base to tip (m)."""
        return float(np.linalg.norm(self.tip_position - self.base_position))


# ---------------------------------------------------------------- #
# Motor-clutch coupling (mock for Unit 3.1; Unit 3.2 replaces)     #
# ---------------------------------------------------------------- #

# Callable signature Unit 3.2 will swap to Worker B's real sampler.
# Inputs: (Lamellipodium) -> per-filament tip load in N.
ClutchTipLoadSampler = Callable[["Lamellipodium"], float]


def zero_load_mock(_lam: "Lamellipodium") -> float:
    """Phase 1 mock: returns zero tip load (free protrusion).

    Use this in tests to exercise the KU-3.6 VALIDATION free-protrusion
    velocity without depending on Worker B's FA population. Unit 3.2
    swaps this for a real
    ``functools.partial(sample_from_population, fa_pop)``.
    """
    return 0.0


def apply_tip_load_from_clutch(
    lam: Lamellipodium,
    sampler: ClutchTipLoadSampler,
) -> Lamellipodium:
    """Set ``lam.last_tip_force_N`` from a clutch sampler (mock-safe)."""
    lam.last_tip_force_N = float(sampler(lam))
    return lam


# ---------------------------------------------------------------- #
# Time advance — Unit 3.2 will implement; Unit 3.1 leaves a stub   #
# ---------------------------------------------------------------- #

class _Unit32NotYetImplemented(NotImplementedError):
    """Marker so Unit 3.2 can grep and replace the stub safely."""


def advance_lamellipodium(
    _lam: Lamellipodium,
    *,
    sampler: ClutchTipLoadSampler,
    dt: float,
) -> Lamellipodium:
    """Stub. Unit 3.2 will implement the full advance step.

    Outline (NOT YET):
      1. Sample per-filament tip load via ``sampler`` (Worker B
         motor-clutch population).
      2. Compute v_p = :func:`brownian_ratchet_velocity` for that load
         and the lam's ``g_actin_uM``.
      3. Subtract ``retrograde_flow_nm_per_s`` (KU-3.12) to get the
         net protrusion velocity v_net.
      4. Update tip position: ``tip += v_net · dt · polarity``,
         clamped to width_um × KU-3.7 max length (3 μm Phase 1).
      5. Return the updated lamellipodium.

    For Unit 3.1 the function raises ``_Unit32NotYetImplemented`` so
    accidental callers cannot silently no-op a time step.
    """
    raise _Unit32NotYetImplemented(
        "advance_lamellipodium is stubbed in Unit 3.1; implement in Unit 3.2 "
        "after Worker B's motor-clutch population freeze."
    )
