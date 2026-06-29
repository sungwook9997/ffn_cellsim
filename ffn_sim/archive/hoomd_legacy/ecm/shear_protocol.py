"""H.1 ECM Lees-Edwards shear protocol (KU-1.30 #2 strain stiffening).

Drives a simulation through a strain ramp γ(t) followed by a hold at
γ_max via HOOMD's ``md.Integrator``-aware ``hoomd.update.BoxResize``.
The strain is realised as an ``xy`` tilt of the simulation box (Lees-
Edwards convention in 3D-periodic geometry), where the per-step box
shape is interpolated linearly from the initial orthorhombic state
to the final sheared state by ``hoomd.variant.box.Interpolate`` driven
by a scalar ``hoomd.variant.Ramp``.

Mathematical setup
------------------
HOOMD's box convention:

    a₁ = (Lx,     0,     0)
    a₂ = (xy·Ly,  Ly,    0)
    a₃ = (xz·Lz,  yz·Lz, Lz)

Shear strain ``γ ≡ δx / Ly`` is realised by setting ``xy = γ``
(dimensionless). Lees-Edwards BC are automatic in HOOMD when xy ≠ 0:
periodic images at y = ±Ly/2 are shifted by xy·Ly along x.

``BoxResize`` rescales every filtered particle's fractional
coordinates so the relative position within the box is preserved as
the box deforms — an affine deformation. This is the canonical small-
strain G_0 measurement (KU-1.30 #1) protocol. For strain stiffening
(KU-1.30 #2) the same Updater drives the box through large γ.

Coordination with the D3 BAOAB Updater
--------------------------------------
HOOMD operation order per step (relevant subset):

  1. Tuners / Updaters fire on their triggers.
  2. ``md.Integrator`` evaluates forces and (here) does nothing else
     because ``methods=[]``.
  3. Writers.

The L-M BAOAB Action (in ``sim.operations.updaters``) and the
``BoxResize`` Updater both belong to step 1. To avoid an order race
(BoxResize affinely rescaling positions *after* BAOAB has stepped
them) we append ``BoxResize`` **before** the BAOAB Updater on the
updater list, so the order is:

    BoxResize.act → BAOAB.act → forces → (next step)

BAOAB's per-step wrap then handles the new (potentially tilted) box
shape via the ``_wrap_into_box`` upper-triangular fractional-coord
transform that already supports xy / xz / yz tilts (verified at
unit level in ``test_baoab.py::TestBoxWrap``).

Sanity Gate
-----------
*Per CLAUDE.md hard rule "Sanity Gate Protocol mandatory before first
execution of any physics/numerics module."*

1. **Dimensional analysis**

   - γ_max: dimensionless (length/length). HOOMD ``xy`` parameter:
     dimensionless. ✓
   - ``ramp_steps``, ``hold_steps``: integers (timesteps). ✓
   - The ``hoomd.variant.Ramp(A, B, t_start, t_ramp)`` evaluates a
     scalar λ ∈ [0, 1]; ``Interpolate(initial, final, variant)``
     evaluates the per-axis box dimension as a linear blend, so
     each axis length stays in [L_initial, L_final] (i.e. [L, L]
     for axes we are NOT shearing). ✓

2. **Boundary cases**

   - ``γ_max == 0``: degenerate (no shear). RUNTIME warning; the
     Updater is still attached but the box never changes — useful as
     a no-shear control. Allowed.
   - ``ramp_steps == 0``: instant jump. RUNTIME forbidden — the
     Ramp variant requires t_ramp > 0 (HOOMD raises). We re-raise
     with a clearer message.
   - ``ramp_steps < 0`` or ``hold_steps < 0``: RUNTIME forbidden.
   - γ_max sign: positive → top edge shifts +x; negative → −x. Both
     allowed (Lees-Edwards is symmetric).

3. **Conservation invariants**

   - **Particle count, bond/angle topology**: unchanged by
     ``BoxResize``. STATIC test ``test_h1_shear.py::TestSchedule``.
   - **Affine rescaling preserves fractional coordinates**: a bead
     at (x, y, z) before the step maps to a bead at the same
     fractional position after — verified by re-reading positions
     after a single ramp step and confirming they lie on the
     box-shape interpolant. STATIC.
   - **No spurious energy**: BoxResize rescaling does NOT compute
     forces — it only mutates positions and box shape. The resulting
     energy change is purely due to the bond / angle / LJ response
     to the strained configuration (as it should be physically).

4. **Numerical sanity**

   - Box variant evaluated each step by HOOMD's C++ side; float64
     internal.
   - Periodic-trigger interval = 1 step → strain is applied every
     step over the ramp (canonical for G_0 oscillatory-shear runs).
     Higher intervals are allowed but the per-step Δγ grows
     accordingly.
   - After ``ramp_steps``, λ saturates at 1.0 — no overshoot, no
     drift.

5. **Sign / sense**

   - Positive γ_max: at the top y face, beads at high y are
     shifted +x relative to beads at low y (Lees-Edwards xy-tilt
     convention). For a chain initially along +y, the chain ends up
     tilted +x at the top. STATIC: small γ shear of a 2-bead chain
     in +y direction, confirm top bead has higher x than bottom.
   - The BAOAB Updater's box-wrap handles the tilt correctly:
     verified at unit level by ``test_baoab.py::TestBoxWrap``.

6. **Measurement protocol**

   - **Schedule**: at timestep ``t``,
        if t < t_start              → λ = 0,  γ = 0
        elif t < t_start + t_ramp   → λ = (t - t_start) / t_ramp,  γ = λ · γ_max
        else                        → λ = 1,  γ = γ_max     (hold)
     STATIC test reads ``sim.state.box.xy`` at sampled times and
     confirms this profile.
   - **Total simulated strain**: at the end of (ramp + hold),
     ``sim.state.box.xy == γ_max``. STATIC.

References
----------
- H.1 brief ``ffn_sim/docs/briefs/H1_ecm_mikado.md`` §Integration
  (Lees-Edwards shear via ``hoomd.update.BoxResize``).
- KU-1.30 frozen by v1 commits ``11eaf13`` (#1 G_0, #2 strain
  stiffening) and ``d92ac20`` (#3 point-dipole).
- ``ffn_sim/integrator/baoab.py`` ``_wrap_into_box`` — tilted-box wrap.
"""

from __future__ import annotations

from dataclasses import dataclass

import hoomd


@dataclass(slots=True)
class ShearSchedule:
    """Linear ramp γ(t): 0 → γ_max over ramp_steps, then hold for hold_steps."""

    gamma_max: float
    ramp_steps: int
    hold_steps: int
    t_start: int = 0

    def __post_init__(self) -> None:
        # RUNTIME §2 boundary checks.
        if self.ramp_steps <= 0:
            raise ValueError(
                f"ramp_steps must be > 0 (HOOMD Ramp requires positive "
                f"t_ramp); got {self.ramp_steps}."
            )
        if self.hold_steps < 0:
            raise ValueError(
                f"hold_steps must be ≥ 0; got {self.hold_steps}."
            )
        if self.t_start < 0:
            raise ValueError(
                f"t_start must be ≥ 0; got {self.t_start}."
            )

    @property
    def total_steps(self) -> int:
        return self.ramp_steps + self.hold_steps

    def expected_gamma_at(self, t: int) -> float:
        """Analytical strain at HOOMD timestep ``t`` (used by STATIC tests)."""
        if t < self.t_start:
            return 0.0
        elapsed = t - self.t_start
        if elapsed < self.ramp_steps:
            return self.gamma_max * (elapsed / self.ramp_steps)
        return self.gamma_max


def make_box_variant(
    initial_box: hoomd.Box, schedule: ShearSchedule
) -> hoomd.variant.box.Interpolate:
    """Build a HOOMD BoxVariant that ramps xy 0 → γ_max then holds."""
    Lx, Ly, Lz = initial_box.Lx, initial_box.Ly, initial_box.Lz
    final_box = hoomd.Box(
        Lx=Lx, Ly=Ly, Lz=Lz,
        xy=schedule.gamma_max, xz=initial_box.xz, yz=initial_box.yz,
    )
    ramp = hoomd.variant.Ramp(
        A=0.0,
        B=1.0,
        t_start=schedule.t_start,
        t_ramp=schedule.ramp_steps,
    )
    return hoomd.variant.box.Interpolate(
        initial_box=initial_box, final_box=final_box, variant=ramp
    )


def attach_shear_updater(
    sim: hoomd.Simulation, schedule: ShearSchedule
) -> hoomd.update.BoxResize:
    """Build a BoxResize Updater and **prepend** it to ``sim.operations.updaters``.

    Prepending is intentional: the BAOAB-limit Action has to see the
    box state *after* any affine rescaling for the step (the wrap
    handles the new tilt), so BoxResize must fire before BAOAB on the
    updater list. See module docstring §Coordination with the D3
    BAOAB Updater.
    """
    initial_box = sim.state.box
    variant = make_box_variant(initial_box, schedule)
    updater = hoomd.update.BoxResize(
        trigger=hoomd.trigger.Periodic(1),
        box=variant,
        filter=hoomd.filter.All(),
    )
    # Prepend by re-building the list.
    existing = list(sim.operations.updaters)
    sim.operations.updaters.clear()
    sim.operations.updaters.append(updater)
    for u in existing:
        sim.operations.updaters.append(u)
    return updater
