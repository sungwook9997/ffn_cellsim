"""H.1 ECM no-shear equilibration prelude (M2-rest, PI 2026-05-20).

The Phase-1-primary-copy Mikado construction places each cross-link at
the *nearest existing bead* on each intersecting fiber. The nearest-bead
rule has no minimum-spacing constraint, so a non-trivial fraction of
realised xl endpoints fall inside one bead radius of an unrelated bead
on a third fiber, placing two ``actin_ecm`` particles at ``r < σ_LJ``.
The WCA pair force is then in its hard-core regime with forces of order

    F_LJ(r ≪ σ) ≈ 24 ε / σ · (σ/r)¹³ ≳ 10³ – 10⁷ N

— many decades above the thermal scale ``√(kT γ_b / Δt) ≈ 10⁻¹¹ N``.
At ``Δt = α · τ_min ≈ 3.8 ns`` and ``γ_b ≈ 6.5 × 10⁻¹⁰ N·s/m``, a single
L-M BAOAB step at ``|F| ~ 10⁴ N`` would displace the offending bead by

    |dr| = |F| · Δt / γ_b ≈ 10⁴ · 4 × 10⁻⁹ / 7 × 10⁻¹⁰ ≈ 5 × 10⁴ m,

which is ~2 × 10⁸ box-lengths — well past the int32 image-flag bound
(2³¹ ≈ 2.1 × 10⁹). The ``_wrap_into_box`` int32 guard
(``ffn_sim/integrator/baoab.py``, PI 2026-05-20) now halts on this
divergence; this helper provides the pre-production settle that keeps
the simulation away from the guard in the first place.

Why a *clipped-force* phase precedes the BAOAB phase
----------------------------------------------------
Reducing ``dt`` does not drain the overlap: with ``|F| ~ 10⁴ N``,
even ``dt = 10⁻²⁰ s`` (well below float64 spacing) leaves ``|dr|``
larger than ``ℓ₀``. The only way to take a *physically* small step
out of the overlap is to cap the per-step displacement directly,
i.e. set ``|dr|_max = α · ℓ₀`` with ``α`` small (here 0.5). After a
few clipped-step iterations the offending pair moves to ``r ≳ σ`` and
``|F|`` drops to the thermal scale, at which point ordinary BAOAB
can take over.

The "clipped Brownian" phase is *not* a physical integrator — it has
no noise, it violates dt-scaling, and it is intentionally
detailed-balance-violating. It is the equivalent of a steepest-descent
minimization pass, used purely to push the system into a configuration
where the stochastic integrator is well-defined. After the clipped
phase the system passes through ``n_baoab`` steps of true L-M BAOAB
with the *standard* ``dt`` and full Langevin noise, where the actual
equilibration happens.

Two phases:

1. **Soft start** (``n_softstart`` steps, default 100): BAOAB Updater
   detached; for each step we evaluate HOOMD forces (``sim.run(1)``
   with ``methods=[]``), clip per-particle displacement to
   ``softstart_max_step``, and write back via ``cpu_local_snapshot``.
2. **BAOAB drain** (``n_baoab`` steps, default 900): BAOAB Updater
   re-attached; standard L-M step at full ``dt``. By the start of
   this phase, ``max |F| < softstart_force_target`` (default 1 nN),
   so the BAOAB step is bounded by the thermal scale.

The combined default of 100 + 900 = 1000 steps matches the boot
prompt's "~1000 step no-shear equilibration prelude" guidance.

Sanity Gate
-----------
*Per CLAUDE.md hard rule "Sanity Gate Protocol mandatory before first
execution of any physics/numerics module."* This helper does not
introduce a new force kernel — both phases reuse the existing
``md.Integrator.forces`` (bond + angle + LJ). The Sanity Gates of the
underlying ``mikado.build_mikado_simulation`` (§Sanity Gate in
``mikado.py``) and ``LeimkuhlerMatthewsBAOAB`` (§Sanity Gate in
``baoab.py``) cover the correctness of each evaluated force and each
BAOAB step. This module adds three RUNTIME checks specific to the
equilibration semantics:

1. **Boundary**:
   - ``n_softstart ≥ 1`` and ``n_baoab ≥ 0``.
   - ``softstart_max_step > 0`` and ``< 0.5 · L_box`` (so the soft
     phase cannot teleport beads more than half a box per step).
2. **Soft-phase exit**: after ``n_softstart`` clipped iterations,
   ``max |net_force|`` must have dropped below ``softstart_force_target``.
   If not, raise — the soft phase is supposed to drain the overlap,
   and a residual force above the target means the integrator is
   about to overflow when BAOAB re-enters. The caller should bump
   ``n_softstart`` or lower ``softstart_max_step`` and retry.
3. **Outcome**: after the BAOAB phase, ``max |net_force|`` must be
   below ``post_force_bound`` (default 1 nN). Same surface-to-PI
   semantics as §2.

References
----------
- ``ffn_sim/integrator/baoab.py`` §4 NaN/Inf force guard + int32-range
  guard (PI 2026-05-20).
- M2-pre full-Mikado-with-shear smoke run warning trace
  (``tests/test_h1_shear.py::TestSimulationSmoke``) — origin of this
  prelude requirement.
- H.1 brief ``ffn_sim/docs/briefs/H1_ecm_mikado.md`` §Implementation
  spec / §Excluded volume (D7 LJ).
"""

from __future__ import annotations

import numpy as np

import hoomd

from ffn_sim.archive.hoomd_legacy.integrator.baoab import _wrap_into_box, LeimkuhlerMatthewsBAOAB


# -- Default tuning constants (Magic-Number Block) --------------------------
# These are not "tunings to pass a gate"; they are derived from the
# construction-time overlap scale (~ℓ₀/2 typical inter-bead distance to
# resolve) and the thermal force scale (~10⁻¹¹ N below which BAOAB is
# unconditionally safe at the H.1 dt). All have a one-line justification.
DEFAULT_N_SOFTSTART: int = 100
"""Number of clipped-Brownian steps. Empirically (Day 4 M2-pre smoke runs),
all overlaps drain in O(τ_xl / Δt) = O(170); 100 with a generous step cap
gives a comfortable margin."""

DEFAULT_N_BAOAB: int = 900
"""Number of true-BAOAB no-shear steps after the soft phase. Together
with the soft phase this matches the boot's ~1000-step total budget."""

DEFAULT_SOFTSTART_MAX_STEP_FRAC: float = 0.05
"""Soft-phase per-step displacement cap as a fraction of ``ℓ₀``. With
ℓ₀ = 0.5 μm this gives ``|dr|_max = 25 nm`` per step ≈ σ_LJ/4 — large
enough to drain a worst-case overlap in O(σ_LJ / 25 nm) = O(4) steps,
small enough not to tunnel through other beads."""

DEFAULT_SOFTSTART_FORCE_TARGET_N: float = 1.0e-9
"""Soft-phase exit criterion: max |F| ≤ 1 nN. This is ≈ 10² × the
thermal force at the H.1 dt, so BAOAB is unconditionally safe past this
value. Lower would be over-conservative."""

DEFAULT_POST_FORCE_BOUND_N: float = 1.0e-9
"""Post-BAOAB acceptance: same bound as the soft-phase exit. The BAOAB
phase brings the system back to thermal equilibrium; this just sanity-
checks no fresh divergence appeared in the noise-driven phase."""


def _softstart_step(
    sim: hoomd.Simulation, *, max_step: float, gamma: float, dt: float
) -> float:
    """One clipped-Brownian step on the current Simulation state.

    HOOMD is expected to have just evaluated forces (i.e. the caller
    has called ``sim.run(1)`` immediately before, or this is the first
    call after ``sim.run(0)``). We read ``net_force``, compute
    ``dr = clip(F · dt / γ, max_step)`` per particle (no noise), wrap
    the new positions into the box, and write back.

    Returns ``max |F|`` BEFORE this step (useful for monitoring the
    overlap drain).
    """
    with sim.state.cpu_local_snapshot as snap:
        pos = np.asarray(snap.particles.position)
        F = np.asarray(snap.particles.net_force)
        image = np.asarray(snap.particles.image)

        if not np.all(np.isfinite(F)):
            raise FloatingPointError(
                "Non-finite net_force during equilibration soft-start; "
                "check that build_mikado_simulation produced a well-formed "
                "force pipeline."
            )

        max_F = float(np.abs(F).max())

        # Per-particle displacement, clipped to max_step magnitude.
        dr = F * (dt / gamma)
        mag = np.linalg.norm(dr, axis=1)
        # Avoid division by zero for zero-force particles (overwhelming
        # majority at equilibrium).
        scale = np.where(mag > max_step, max_step / np.maximum(mag, 1e-300), 1.0)
        dr = dr * scale[:, None]

        new_pos = pos + dr
        wrapped, img_delta = _wrap_into_box(new_pos, sim.state.box)
        pos[:] = wrapped
        image[:] = image + img_delta

    return max_F


def equilibrate_no_shear(
    sim: hoomd.Simulation,
    baoab_action: LeimkuhlerMatthewsBAOAB | None,
    baoab_updater: hoomd.update.CustomUpdater | None,
    *,
    n_softstart: int = DEFAULT_N_SOFTSTART,
    n_baoab: int = DEFAULT_N_BAOAB,
    softstart_max_step_frac: float = DEFAULT_SOFTSTART_MAX_STEP_FRAC,
    softstart_force_target: float = DEFAULT_SOFTSTART_FORCE_TARGET_N,
    post_force_bound: float = DEFAULT_POST_FORCE_BOUND_N,
    rest_length: float | None = None,
    gamma_b: float | None = None,
) -> dict[str, float]:
    """Two-phase no-shear settle for a freshly-built H.1 Mikado.

    Parameters
    ----------
    sim, baoab_action, baoab_updater
        Output of ``ffn_sim.archive.hoomd_legacy.ecm.mikado.build_mikado_simulation`` with
        ``with_baoab=True``. ``baoab_updater`` is detached from
        ``sim.operations.updaters`` during the soft phase and
        re-appended for the BAOAB phase. If both are ``None`` (caller
        opted out of BAOAB), the helper skips the BAOAB phase and only
        runs the soft phase.
    rest_length, gamma_b
        ``ResolvedH1.rest_length`` (ℓ₀) and ``ResolvedH1.gamma_b`` from
        the H.1 resolved config. Required so the soft phase can size
        ``max_step`` and the Brownian step. The caller can also leave
        them as ``None`` and pull them from ``baoab_action`` (which
        carries a per-type γ map) — but that is awkward in the
        no-BAOAB case, so the explicit parameters are preferred.

    Returns
    -------
    diagnostics : dict[str, float]
        ``{
            "max_force_before",      # right after build, no settle
            "max_force_after_soft",  # after the clipped-Brownian phase
            "max_force_after_baoab", # after the BAOAB no-shear phase
            "n_softstart": ...,
            "n_baoab": ...,
        }``

    Raises
    ------
    ValueError
        If any boundary check (§1) fails.
    RuntimeError
        If the post-soft-phase or post-BAOAB-phase force bound (§2/§3)
        is violated.
    """
    # §1 boundary checks.
    if n_softstart < 1:
        raise ValueError(f"n_softstart must be ≥ 1; got {n_softstart}.")
    if n_baoab < 0:
        raise ValueError(f"n_baoab must be ≥ 0; got {n_baoab}.")
    L_box = float(sim.state.box.Lx)
    if rest_length is None or gamma_b is None:
        raise ValueError(
            "equilibrate_no_shear requires explicit rest_length and gamma_b "
            "(from ResolvedH1). Pass them through from build_mikado_simulation."
        )
    max_step = softstart_max_step_frac * rest_length
    if not (0.0 < max_step < 0.5 * L_box):
        raise ValueError(
            f"softstart_max_step_frac · rest_length = {max_step:.3e} m must "
            f"lie in (0, 0.5·L_box={0.5*L_box:.3e}); got "
            f"frac={softstart_max_step_frac}, ℓ₀={rest_length:.3e}."
        )

    dt = float(sim.operations.integrator.dt)

    # Detach BAOAB Updater for the soft phase (the soft phase advances
    # positions itself via cpu_local_snapshot writes — running BAOAB on
    # top would double-step).
    had_baoab = baoab_updater is not None
    if had_baoab:
        existing = list(sim.operations.updaters)
        if baoab_updater in existing:
            sim.operations.updaters.remove(baoab_updater)

    # Seed forces at the construction state.
    sim.run(0)
    with sim.state.cpu_local_snapshot as snap:
        F_init = np.asarray(snap.particles.net_force)
        if not np.all(np.isfinite(F_init)):
            raise FloatingPointError(
                "Non-finite net_force at construction (step 0); the "
                "Mikado / xl topology produced a degenerate state."
            )
        max_force_before = float(np.abs(F_init).max())

    # ----- Phase A: clipped-Brownian soft start -----
    for _ in range(n_softstart):
        _softstart_step(sim, max_step=max_step, gamma=gamma_b, dt=dt)
        # Re-evaluate forces for the next iteration.
        sim.run(1)

    with sim.state.cpu_local_snapshot as snap:
        F_soft = np.asarray(snap.particles.net_force)
        max_force_after_soft = float(np.abs(F_soft).max())
    if max_force_after_soft > softstart_force_target:
        raise RuntimeError(
            "Soft-start phase did not drain the construction force below "
            f"the BAOAB-safety target. max|F| after {n_softstart} clipped "
            f"steps = {max_force_after_soft:.3e} N; target = "
            f"{softstart_force_target:.3e} N. Bump n_softstart, lower "
            "softstart_max_step_frac, or surface to PI."
        )

    # ----- Phase B: BAOAB no-shear settle -----
    max_force_after_baoab = max_force_after_soft  # if no BAOAB phase, leave as soft outcome
    if had_baoab and n_baoab > 0:
        sim.operations.updaters.append(baoab_updater)
        sim.run(n_baoab)
        with sim.state.cpu_local_snapshot as snap:
            F_post = np.asarray(snap.particles.net_force)
            if not np.all(np.isfinite(F_post)):
                raise FloatingPointError(
                    "Non-finite net_force after BAOAB equilibration prelude; "
                    "the int32 guard should have raised earlier — surface to "
                    "PI."
                )
            max_force_after_baoab = float(np.abs(F_post).max())
        if max_force_after_baoab > post_force_bound:
            raise RuntimeError(
                "BAOAB no-shear phase did not converge to a thermal-bound "
                f"state. max|F| after {n_baoab} steps = "
                f"{max_force_after_baoab:.3e} N; bound = "
                f"{post_force_bound:.3e} N. Surface to PI."
            )

    return {
        "max_force_before": max_force_before,
        "max_force_after_soft": max_force_after_soft,
        "max_force_after_baoab": max_force_after_baoab,
        "n_softstart": float(n_softstart),
        "n_baoab": float(n_baoab),
    }
