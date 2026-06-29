"""B1 — global CFL dt reconciliation across all active cell subsystems.

PI-ratified item B1 (PI_DECISION_QUEUE_2026-05-30 §B1, audit §1#15).

Motivation
----------
``build_cortex_full_simulation`` builds the HOOMD integrator at the
*cortex* CFL timestep ``ResolvedH3.dt_cfl = α · min(τ_xl, τ_stretch,
τ_bend)``. That dt is correct for a cortex-only system. But once stiff
optional subsystems are wired in — in particular H.4 focal adhesions
(integrin bond ``k_int_bare``) — the stiffest bond in the *whole* system
may have a relaxation timescale ``τ = γ / k`` shorter than any cortex
timescale. Integrating that bond at the cortex dt then violates the CFL
condition ``dt ≤ α · τ``; the overdamped Langevin step
``|dr| = |F| dt / γ`` overshoots and the bond energy grows step-on-step
until the int32 image-flag guard in ``baoab.py`` halts the run.

This module is the **build-layer** guard that catches that case BEFORE
the production loop starts. It does NOT touch the frozen integrator: it
only *computes* a binding dt and (optionally, opt-in) lowers the
integrator dt to it.

CFL timescale per element
-------------------------
For an overdamped (Brownian / L-M BAOAB-limit) bead of drag ``γ``
tethered by a harmonic spring of stiffness ``k`` [N/m], the relaxation
time is

    τ = γ / k                                        [s]

and the stable timestep is ``dt ≤ safety · τ`` where ``safety`` is the
cortex ``cfl_safety_factor`` (already ratified, KU-1.26). The binding
global timestep is

    dt_min = safety · min over all active stiff translational springs (γ_e / k_e)

Drag values mirror exactly the per-type ``gamma_map`` that
``cell.build_cortex_full_simulation`` hands the BAOAB updater:

============  ===============================  ==============================
subsystem     stiff spring(s) [N/m]            drag γ [N·s/m]
============  ===============================  ==============================
cortex        bond_k (→ tau_min folds          gamma_b
              stretch/bend/xl)
xlinks        k_intra, k_attach                gamma_b (cortex)
myosin        k_head_spring, k_backbone        gamma_b (cortex)
erm           k_ERM                            gamma_b (cortex)
fa            k_int_bare                       gamma_integrin
lamellipodium bond_k, k_wave_pin               6πη·bead_radius
enclosed_vol  k_eff = 3 K_vol S²/(N² V0)       gamma_b (cortex)
membrane      (tension reservoir; report-only) —
turnover      (rate-based; no new spring)      —
============  ===============================  ==============================

For the cortex we trust ``ResolvedH3.tau_min`` (the module's own
sanity-gated CFL timescale, already = min(τ_xl, τ_stretch, τ_bend)). For
every other subsystem we compute τ = γ/k explicitly from the raw
stiffness fields and the matching drag, because those dataclasses do not
pre-derive a combined τ. The membrane (a bounded tension/ratchet load,
not a stiff Hookean translational bond the inner BAOAB step must resolve)
is reported for transparency only and is NOT folded into dt_min.

Event-rate CFL (batched stochastic updaters)
--------------------------------------------
Several subsystems carry ``bond_event_rate_max`` / ``k_off_max`` (1/s), a
*kinetic* ceiling rather than a mechanical spring. The batched updaters
already self-limit their batch window against it (cell.py FA batch-window
reconciliation). We report it in the breakdown for transparency but do
NOT fold it into ``dt_min`` — it constrains batch stride, not the inner
BAOAB step.

SHAKE / constraint timescale separation (constrained=True)
----------------------------------------------------------
When the cortex actin backbone is run as a rigid M-SHAKE constraint
(``constrained=True``), the backbone *stretch* bond is no longer a soft
spring the timestep must resolve — the constraint removes that fast
degree of freedom. The binding soft timescale is then the cortex bend
plus the next-stiffest spring (xl / fa / myosin). When a constrained dt
is supplied, :func:`compute_global_cfl_dt` asserts it is within the
remaining soft CFL bound (``dt ≤ safety · τ_soft_min``); if a soft
subsystem is stiffer than the constrained dt assumes, it raises so the
caller can lower the constrained dt rather than silently overshoot.

Sanity Gate
-----------
*Per CLAUDE.md hard rule. STATIC + RUNTIME checks in
``ffn_sim/tests/test_dt_equilibration.py``.*

1. **Dimensional**: every τ = γ/k → [N·s/m]/[N/m] = [s]; dt_min = α·τ
   → [s]. ``safety`` dimensionless. Asserted finite-positive.
2. **Boundary**: with only cortex active, ``dt_min`` equals the cortex
   ``cfl_safety_factor · tau_min`` (= ``ResolvedH3.dt_cfl``). Adding any
   subsystem can only *lower* dt_min (monotone min).
3. **Monotonicity**: ``dt_min ≤ safety · τ_e`` for every element e.
4. **Numerical**: a zero/negative k or γ raises (degenerate subsystem).
5. **Sign**: dt_min > 0 always.

References
----------
- ``ffn_sim/cortex/cortex.py`` ``resolve_h3_derived`` (dt_cfl derivation).
- ``ffn_sim/ecm/equilibrate.py`` (the overlap-drain prelude this pairs with).
- PI_DECISION_QUEUE_2026-05-30 §B1; MECHANISM_AUDIT_2026-05-30 §1#15.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)


def _stokes_drag(eta: float, R: float) -> float:
    """γ = 6π η R  [N·s/m]  (Stokes drag for a sphere of radius R)."""
    return 6.0 * math.pi * eta * R


@dataclass(slots=True)
class CFLElement:
    """One CFL-constraining element in the global dt breakdown.

    Attributes
    ----------
    subsystem
        Owning subsystem (``"cortex"``, ``"fa"``, ...).
    label
        Element name (``"k_int_bare"``, ``"tau_min"``, ...).
    k
        Translational spring stiffness [N/m], or ``math.nan`` when the
        element is reported via a pre-derived τ rather than a raw k.
    gamma
        Drag [N·s/m], or ``math.nan`` when reported via a pre-derived τ.
    tau
        Relaxation timescale [s] (= γ/k, or the module's own derived τ).
    dt_bound
        ``safety · tau`` [s].
    is_constraint_removed
        True if this soft spring is removed by the SHAKE constraint
        (cortex stretch when ``constrained=True``) and so does NOT
        constrain the soft timestep.
    is_informational
        True if reported for transparency only (membrane tension /
        event rates) and NOT folded into ``dt_min``.
    """

    subsystem: str
    label: str
    k: float
    gamma: float
    tau: float
    dt_bound: float
    is_constraint_removed: bool = False
    is_informational: bool = False


@dataclass(slots=True)
class GlobalCFLResult:
    """Result of :func:`compute_global_cfl_dt`.

    Attributes
    ----------
    dt_min
        Binding global CFL timestep [s] = ``safety · min τ`` over all
        active, non-removed, non-informational elements.
    safety
        CFL safety factor used (the cortex ``cfl_safety_factor``).
    binding
        The :class:`CFLElement` that set ``dt_min`` (the stiffest one).
    elements
        Per-element breakdown (includes removed / informational ones,
        flagged).
    event_rates
        ``{subsystem: bond_event_rate_max}`` reported for transparency.
    """

    dt_min: float
    safety: float
    binding: CFLElement
    elements: list[CFLElement] = field(default_factory=list)
    event_rates: dict[str, float] = field(default_factory=dict)

    def format_breakdown(self) -> str:
        """Return a multi-line human-readable breakdown for logging."""
        lines = [
            f"global CFL dt reconciliation: dt_min = {self.dt_min:.4e} s "
            f"(safety = {self.safety:g}); binding element = "
            f"{self.binding.subsystem}.{self.binding.label} "
            f"(tau = {self.binding.tau:.4e} s)"
        ]
        for e in self.elements:
            tag = ""
            if e.is_constraint_removed:
                tag = " [SHAKE-removed]"
            elif e.is_informational:
                tag = " [info-only]"
            if math.isfinite(e.k) and math.isfinite(e.gamma):
                lines.append(
                    f"  {e.subsystem:>14s}.{e.label:<18s} "
                    f"k={e.k:.3e} N/m  gamma={e.gamma:.3e}  "
                    f"tau={e.tau:.3e} s  dt<={e.dt_bound:.3e} s{tag}"
                )
            else:
                lines.append(
                    f"  {e.subsystem:>14s}.{e.label:<18s} "
                    f"tau={e.tau:.3e} s  dt<={e.dt_bound:.3e} s{tag}"
                )
        for sub, rate in self.event_rates.items():
            lines.append(
                f"  {sub:>14s}.{'bond_event_rate_max':<18s} "
                f"rate={rate:.3e} 1/s  (batch-stride only; not in dt_min)"
            )
        return "\n".join(lines)


def _tau_from_k_gamma(k: float, gamma: float, *, subsystem: str, label: str) -> float:
    """τ = γ/k with finite-positive guards (Sanity Gate §4)."""
    if not (math.isfinite(k) and k > 0.0):
        raise ValueError(
            f"{subsystem}.{label}: stiffness k must be finite-positive; got {k!r}."
        )
    if not (math.isfinite(gamma) and gamma > 0.0):
        raise ValueError(
            f"{subsystem}.{label}: drag gamma must be finite-positive; got {gamma!r}."
        )
    return gamma / k


def compute_global_cfl_dt(
    p_cortex: Any,
    *,
    p_xlinks: Any | None = None,
    p_myosin: Any | None = None,
    p_lamellipodium: Any | None = None,
    p_fa: Any | None = None,
    p_erm: Any | None = None,
    p_enclosed_volume: Any | None = None,
    p_turnover: Any | None = None,
    p_membrane: Any | None = None,
    constrained: bool = False,
    constrained_dt: float | None = None,
    safety: float | None = None,
) -> GlobalCFLResult:
    """Compute the binding global CFL timestep across active subsystems.

    Parameters
    ----------
    p_cortex
        Resolved cortex params (:class:`ffn_sim.archive.hoomd_legacy.cortex.cortex.ResolvedH3`).
        Always required — supplies the safety factor, the cortex
        ``gamma_b`` drag, and the cortex ``tau_min``.
    p_xlinks, p_myosin, p_lamellipodium, p_fa, p_erm, p_enclosed_volume,
    p_turnover, p_membrane
        Optional resolved subsystem params. Each non-``None`` one
        contributes its stiff translational spring(s).
    constrained
        If True, the cortex actin backbone is a rigid M-SHAKE constraint;
        its soft *stretch* spring is removed from the binding min and the
        constrained dt is asserted separated from the remaining soft τ's.
    constrained_dt
        The requested constrained timestep (when ``constrained=True``).
    safety
        Override the CFL safety factor. Default ``p_cortex.cfl_safety_factor``.

    Returns
    -------
    GlobalCFLResult

    Raises
    ------
    ValueError
        On a degenerate k/γ, or if the SHAKE timescale-separation
        assertion fails when ``constrained``.
    """
    if safety is None:
        safety = float(p_cortex.cfl_safety_factor)
    if not (math.isfinite(safety) and safety > 0.0):
        raise ValueError(f"cfl_safety_factor must be finite-positive; got {safety!r}.")

    gamma_cortex = float(p_cortex.gamma_b)
    eta = float(p_cortex.water_viscosity)

    elements: list[CFLElement] = []
    event_rates: dict[str, float] = {}

    def add(subsystem: str, label: str, *, k: float, gamma: float,
            tau: float | None = None, constraint_removed: bool = False,
            informational: bool = False) -> None:
        if tau is None:
            tau = _tau_from_k_gamma(k, gamma, subsystem=subsystem, label=label)
        elif not (math.isfinite(tau) and tau > 0.0):
            raise ValueError(
                f"{subsystem}.{label}: derived tau must be finite-positive; got {tau!r}."
            )
        elements.append(
            CFLElement(
                subsystem=subsystem, label=label, k=k, gamma=gamma, tau=tau,
                dt_bound=safety * tau, is_constraint_removed=constraint_removed,
                is_informational=informational,
            )
        )

    # --- Cortex (always) -------------------------------------------------
    # ResolvedH3.tau_min = min(tau_xl, tau_stretch, tau_bend), the cortex
    # module's own sanity-gated CFL timescale. When constrained, the
    # backbone stretch DOF is SHAKE-removed; the cortex soft contribution
    # becomes tau_bend (bend survives a bond-length constraint), and the
    # stretch contribution is flagged removed.
    if constrained:
        add("cortex", "tau_bend", k=p_cortex.bond_k, gamma=gamma_cortex,
            tau=p_cortex.tau_bend)
        add("cortex", "tau_stretch", k=p_cortex.bond_k, gamma=gamma_cortex,
            tau=p_cortex.tau_stretch, constraint_removed=True)
    else:
        add("cortex", "tau_min", k=p_cortex.bond_k, gamma=gamma_cortex,
            tau=p_cortex.tau_min)

    # --- Crosslinkers (ride on cortex gamma_b) ---------------------------
    if p_xlinks is not None and getattr(p_xlinks, "n_xl", 0) > 0:
        add("xlinks", "k_intra", k=p_xlinks.k_intra, gamma=gamma_cortex)
        add("xlinks", "k_attach", k=p_xlinks.k_attach, gamma=gamma_cortex)
        if getattr(p_xlinks, "k_off_max", 0.0):
            event_rates["xlinks"] = float(p_xlinks.k_off_max)

    # --- Myosin (ride on cortex gamma_b) ---------------------------------
    if p_myosin is not None and getattr(p_myosin, "n_motors_per_cell", 0) > 0:
        add("myosin", "k_head_spring", k=p_myosin.k_head_spring, gamma=gamma_cortex)
        if getattr(p_myosin, "k_backbone", 0.0):
            add("myosin", "k_backbone", k=p_myosin.k_backbone, gamma=gamma_cortex)
        if getattr(p_myosin, "k_off_max_at_zero_load", 0.0):
            event_rates["myosin"] = float(p_myosin.k_off_max_at_zero_load)

    # --- ERM tether (rides on cortex gamma_b) ----------------------------
    if p_erm is not None:
        add("erm", "k_ERM", k=p_erm.k_ERM, gamma=gamma_cortex)

    # --- Focal adhesions (the B1 motivating case) ------------------------
    n_fa = 0
    if p_fa is not None:
        n_fa = (getattr(p_fa, "n_nascent_per_cell", 0)
                + getattr(p_fa, "n_mature_per_cell", 0))
    if p_fa is not None and n_fa > 0:
        add("fa", "k_int_bare", k=p_fa.k_int_bare, gamma=p_fa.gamma_integrin)
        if getattr(p_fa, "bond_event_rate_max", 0.0):
            event_rates["fa"] = float(p_fa.bond_event_rate_max)

    # --- Lamellipodium (γ = 6πη·bead_radius, matches cell.py) ------------
    if p_lamellipodium is not None and getattr(p_lamellipodium, "n_WAVE", 0) > 0:
        gamma_lamel = _stokes_drag(eta, p_lamellipodium.bead_radius)
        add("lamellipodium", "bond_k", k=p_lamellipodium.bond_k, gamma=gamma_lamel)
        if getattr(p_lamellipodium, "k_wave_pin", 0.0):
            add("lamellipodium", "k_wave_pin", k=p_lamellipodium.k_wave_pin,
                gamma=gamma_lamel)

    # --- Enclosed-volume pressure (rides on cortex gamma_b) --------------
    # k_eff = 3 K_vol S² / (N² V0); S = 4π R². N (shell bead count) is the
    # cortex actin count = n_filaments · beads_per_filament. Mirrors
    # EnclosedVolumePressure.effective_bead_stiffness without needing the
    # constructed force object.
    if p_enclosed_volume is not None:
        n_shell = int(p_cortex.n_filaments) * int(p_cortex.beads_per_filament)
        S = 4.0 * math.pi * p_enclosed_volume.R_cell ** 2
        k_eff = (3.0 * p_enclosed_volume.K_vol * S * S
                 / (float(n_shell) ** 2 * p_enclosed_volume.V0))
        if k_eff > 0.0:
            add("enclosed_volume", "k_eff", k=k_eff, gamma=gamma_cortex)

    # --- Membrane load (tension reservoir; informational only) ----------
    # The membrane is a bounded tension/Brownian-ratchet load, not a stiff
    # Hookean translational bond stepped inside BAOAB; its diagnostic tip
    # stiffness needs a runtime contact count. Report γ_mem for
    # transparency; do NOT fold into dt_min.
    if p_membrane is not None:
        gamma_mem = float(getattr(p_membrane, "gamma_mem", float("nan")))
        elements.append(
            CFLElement(
                subsystem="membrane", label="gamma_mem(tension)",
                k=gamma_mem, gamma=float("nan"), tau=float("nan"),
                dt_bound=float("nan"), is_informational=True,
            )
        )

    # --- Turnover (rate-based; no new stiff translational spring) --------
    if p_turnover is not None and getattr(p_turnover, "k_sev_max", 0.0):
        event_rates["turnover"] = float(p_turnover.k_sev_max)

    # --- Bind the min over active (non-removed, non-informational) -------
    active = [e for e in elements
              if not e.is_constraint_removed and not e.is_informational]
    if not active:
        raise ValueError("No active CFL elements (cortex must always be present).")
    binding = min(active, key=lambda e: e.tau)
    dt_min = safety * binding.tau

    result = GlobalCFLResult(
        dt_min=dt_min, safety=safety, binding=binding,
        elements=elements, event_rates=event_rates,
    )

    # --- SHAKE timescale-separation assertion (Sanity Gate, constrained) -
    if constrained and constrained_dt is not None:
        if constrained_dt > dt_min * (1.0 + 1e-9):
            raise ValueError(
                "SHAKE timescale separation violated: requested constrained_dt "
                f"= {constrained_dt:.4e} s exceeds the binding soft CFL bound "
                f"dt_min = {dt_min:.4e} s (binding element "
                f"{binding.subsystem}.{binding.label}, tau={binding.tau:.4e} s). "
                "Lower constrained_dt or surface to PI. Breakdown:\n"
                + result.format_breakdown()
            )

    return result
