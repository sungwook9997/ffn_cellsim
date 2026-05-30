"""B2 — cortex/cell equilibration prelude (soft-start + BAOAB drain).

PI-ratified item B2 (PI_DECISION_QUEUE_2026-05-30 §B2, audit §1#11;
H4_FA_INTEGRATION_DESIGN §S9).

Why a prelude is needed for the integrated cell
-----------------------------------------------
``build_cortex_full_simulation`` composes the cortex shell, optional
crosslinkers / myosin / lamellipodium, and — critically for H.4 — the
focal-adhesion integrins + immobile substrate ligands. At *construction*:

* The cortex shell sits at ``r ≈ R_cell`` (a sphere).
* The substrate ligands sit on the ``z = 0`` plane.

These two sets are spatially **disjoint**: no integrin is within clutch
capture radius of a substrate ligand, and the cortex bottom is one cell
radius above the substrate. For FA to engage, the cell must first
*settle onto* the substrate so clutch bonds can form. Starting the
production loop directly would either (a) form no clutches (cell
floating) or (b), if a clutch is pre-seeded across the gap, snap a
hugely-stretched spring that explodes the first BAOAB step.

Separately — exactly as in the H.1 ECM case (``ecm/equilibrate.py``) —
construction can place excluded-volume (WCA) pairs inside the hard-core
where ``|F|`` is many decades above thermal; a direct BAOAB step then
overshoots and trips the int32 image guard.

This module is the cell-layer analogue of
``ecm.equilibrate.equilibrate_no_shear``: a two-phase soft-start that
(1) drains construction overlaps with capped-displacement steepest-
descent, then (2) runs true L-M BAOAB at the integrator dt to thermalise
— BEFORE the production / measurement loop.

It does NOT modify the frozen integrator. It reuses the existing
``ecm.equilibrate._softstart_step`` (clipped-Brownian phase) and the
BAOAB updater already attached by the builder, so it introduces no new
integrator and no new force kernel.

Relationship to ``ecm.equilibrate.equilibrate_no_shear``
--------------------------------------------------------
``equilibrate_no_shear`` is H.1-owned and assumes a single bead type
with one ``rest_length``/``gamma_b``. The cell has multiple particle
types with different drags. The clipped-Brownian step caps displacement
per particle using a single representative ``gamma``; we pass the cortex
``gamma_b`` (the overwhelming majority of beads) and a displacement cap
sized to the cortex rest length, which is conservative for the lighter
FA beads. The substrate-ligand pin updater (when FA is on) re-pins
ligands every step, so their soft-start displacement is harmless. We
therefore *reuse* ``_softstart_step`` directly rather than duplicate it,
and never edit the H.1-owned module.

Sanity Gate
-----------
*Per CLAUDE.md hard rule. Checks in
``ffn_sim/tests/test_dt_equilibration.py``.*

1. **Boundary**: ``n_softstart ≥ 0`` and ``n_baoab ≥ 0``. With both 0
   the helper is a no-op returning the construction-state max force.
2. **Numerical**: every read of ``net_force`` is finite-checked; a
   non-finite force raises (the int32 guard should fire first).
3. **Non-increase**: the post-prelude ``max|F|`` must not exceed the
   construction ``max|F|`` (the prelude drains, it does not inject
   energy). Asserted in the test; reported in the diagnostics here.
4. **Strict opt-in**: ``strict_force_target`` reproduces the H.1 raise
   semantics; default is warn-only because the integrated cell can carry
   a legitimate residual from cross-subsystem attractive pairs / the FA
   gap.

References
----------
- ``ffn_sim/ecm/equilibrate.py`` (the H.1 soft-start this mirrors).
- ``ffn_sim/cell/dt_reconcile.py`` (B1 — the companion dt guard).
- PI_DECISION_QUEUE_2026-05-30 §B2; H4_FA_INTEGRATION_DESIGN §S9.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

import hoomd

from ffn_sim.ecm.equilibrate import _softstart_step

logger = logging.getLogger(__name__)


# Default prelude budget. Mirrors the H.1 100 + 900 split rationale
# (ecm/equilibrate.py) but is OPT-IN here: the builder default is 0
# (no prelude, no behaviour change).
DEFAULT_N_SOFTSTART: int = 100
DEFAULT_N_BAOAB: int = 900

# Soft-phase per-step displacement cap as a fraction of the cortex rest
# length (same derivation as ecm/equilibrate: large enough to drain a
# worst-case overlap in O(few) steps, small enough not to tunnel a bead
# through a neighbour).
DEFAULT_SOFTSTART_MAX_STEP_FRAC: float = 0.05


def _max_abs_net_force(sim: hoomd.Simulation) -> float:
    """Return max |net_force| over all particles (finite-checked)."""
    with sim.state.cpu_local_snapshot as snap:
        F = np.asarray(snap.particles.net_force)
        if not np.all(np.isfinite(F)):
            raise FloatingPointError(
                "Non-finite net_force during cell equilibration; the int32 "
                "image guard in baoab.py should have raised first — surface "
                "to PI."
            )
        return float(np.abs(F).max())


def equilibrate_cell(
    handles: dict[str, Any],
    *,
    n_softstart: int = DEFAULT_N_SOFTSTART,
    n_baoab: int = DEFAULT_N_BAOAB,
    softstart_max_step_frac: float = DEFAULT_SOFTSTART_MAX_STEP_FRAC,
    rest_length: float,
    gamma_b: float,
    strict_force_target: float | None = None,
) -> dict[str, float]:
    """Soft-start + BAOAB equilibration prelude on a built cell.

    Operates on the ``handles`` dict returned by
    :func:`ffn_sim.cell.cell.build_cortex_full_simulation`. Drains
    construction overlaps (and lets the cell settle onto the substrate
    when FA is wired) BEFORE the production loop, without touching the
    frozen integrator.

    Parameters
    ----------
    handles
        Builder output. Must contain ``"sim"``; for the BAOAB phase it
        must contain ``"baoab_updater"`` (the attached CustomUpdater).
        When that is ``None`` only the soft phase runs.
    n_softstart
        Clipped-Brownian steepest-descent steps (default 100). 0 skips.
    n_baoab
        True L-M BAOAB steps at the integrator dt (default 900). 0 skips.
    softstart_max_step_frac
        Soft-phase per-step displacement cap as a fraction of
        ``rest_length``.
    rest_length, gamma_b
        Cortex ``ResolvedH3.rest_length`` (ℓ₀) and ``gamma_b`` — size the
        soft-phase displacement cap and the clipped-Brownian step.
    strict_force_target
        If not ``None``, raise ``RuntimeError`` when ``max|F|`` after the
        soft phase exceeds this value (the strict H.1 contract). Default
        ``None`` → warn-only.

    Returns
    -------
    dict[str, float]
        ``{"max_force_before", "max_force_after_soft",
           "max_force_after_baoab", "n_softstart", "n_baoab"}`` — same
        schema as ``equilibrate_no_shear``.

    Raises
    ------
    ValueError
        On a bad budget / displacement cap / gamma.
    KeyError
        If ``handles`` lacks ``"sim"``.
    RuntimeError
        If ``strict_force_target`` is set and the soft phase does not
        drain below it.
    """
    if n_softstart < 0:
        raise ValueError(f"n_softstart must be ≥ 0; got {n_softstart}.")
    if n_baoab < 0:
        raise ValueError(f"n_baoab must be ≥ 0; got {n_baoab}.")
    if "sim" not in handles:
        raise KeyError(
            "handles must contain 'sim' (build_cortex_full_simulation output)."
        )

    sim: hoomd.Simulation = handles["sim"]
    baoab_updater = handles.get("baoab_updater")

    L_box = float(sim.state.box.Lx)
    max_step = softstart_max_step_frac * rest_length
    if not (0.0 < max_step < 0.5 * L_box):
        raise ValueError(
            f"softstart_max_step_frac · rest_length = {max_step:.3e} m must lie "
            f"in (0, 0.5·L_box={0.5 * L_box:.3e}); got "
            f"frac={softstart_max_step_frac}, ℓ₀={rest_length:.3e}."
        )
    if not (gamma_b > 0.0):
        raise ValueError(f"gamma_b must be > 0; got {gamma_b!r}.")

    dt = float(sim.operations.integrator.dt)

    # No-op fast path: both phases off → report construction force only.
    if n_softstart == 0 and n_baoab == 0:
        sim.run(0)
        mf = _max_abs_net_force(sim)
        return {
            "max_force_before": mf,
            "max_force_after_soft": mf,
            "max_force_after_baoab": mf,
            "n_softstart": 0.0,
            "n_baoab": 0.0,
        }

    # Detach BAOAB updater during the soft phase (it advances positions
    # itself; running BAOAB on top would double-step).
    had_baoab = baoab_updater is not None
    if had_baoab and baoab_updater in list(sim.operations.updaters):
        sim.operations.updaters.remove(baoab_updater)

    # Seed forces at the construction state.
    sim.run(0)
    max_force_before = _max_abs_net_force(sim)

    # ----- Phase A: clipped-Brownian soft start -----
    max_force_after_soft = max_force_before
    if n_softstart > 0:
        for _ in range(n_softstart):
            _softstart_step(sim, max_step=max_step, gamma=gamma_b, dt=dt)
            sim.run(1)
        max_force_after_soft = _max_abs_net_force(sim)
        if strict_force_target is not None and max_force_after_soft > strict_force_target:
            # Re-attach before raising so the caller's handles stay valid.
            if had_baoab:
                sim.operations.updaters.append(baoab_updater)
            raise RuntimeError(
                "Cell soft-start did not drain construction force below the "
                f"strict target. max|F| after {n_softstart} steps = "
                f"{max_force_after_soft:.3e} N; target = {strict_force_target:.3e} N. "
                "Bump n_softstart, lower softstart_max_step_frac, or surface to PI."
            )
        if max_force_after_soft > max_force_before:
            logger.warning(
                "Cell soft-start max|F| rose from %.3e to %.3e N (expected to "
                "drain). Check construction overlaps / FA gap.",
                max_force_before, max_force_after_soft,
            )

    # ----- Phase B: BAOAB drain -----
    max_force_after_baoab = max_force_after_soft
    if had_baoab and n_baoab > 0:
        sim.operations.updaters.append(baoab_updater)
        sim.run(n_baoab)
        max_force_after_baoab = _max_abs_net_force(sim)
    elif (not had_baoab) and n_baoab > 0:
        logger.warning(
            "equilibrate_cell: n_baoab=%d requested but no baoab_updater in "
            "handles; BAOAB phase skipped.", n_baoab,
        )

    return {
        "max_force_before": max_force_before,
        "max_force_after_soft": max_force_after_soft,
        "max_force_after_baoab": max_force_after_baoab,
        "n_softstart": float(n_softstart),
        "n_baoab": float(n_baoab),
    }
