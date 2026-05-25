"""H.3 Cell composition class (v2 — replaces deleted v1 `acs_kb/cell/cell.py`).

Phase 1 H.3 brief Deliverables: `ffn_sim/cell/cell.py` — new v2 `Cell`
class, HOOMD particle-group-backed, owns:

* cortex (`ffn_sim.cortex.cortex` topology + ResolvedH3)
* ERM tether (`ffn_sim.cortex.erm.ERMHarmonic` custom force)
* dynamic crosslinkers (`ffn_sim.cortex.crosslinkers.XlinkBondUpdater`)
* lamellipodium slot (H.5 deliverable — None for H.3)
* FA / adhesion slot (H.4 Sub-owned — None for H.3 from Main side)
* myosin slot (H.3 next-iteration deliverable — None for now)

PER CLAUDE.md the v1 `acs_kb/cell/cell.py` (single-chain Cortex-coupled)
was DELETED in the 2026-05-20 rename and MUST NOT be resurrected.
This v2 inverts the v1 architecture: HOOMD particle dynamics is the
runtime, cell composition is bookkeeping over particle groups + force
computes + Updaters.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC checks in
``ffn_sim/tests/test_cell.py``.*

1. **Dimensional analysis**
   - All particle-count / bond-count diagnostics are integer counts.
     No unit-bearing scalars introduced here (Cell wraps already-
     resolved Cortex / ERM / Crosslinker dataclasses; the underlying
     SI dimensions are validated in those modules' own sanity gates).

2. **Boundary cases**
   - `Cell.build(...)` with all add-ons None: builds a bare cortex sim
     (equivalent to ``cortex.build_cortex_simulation``).
   - `Cell.with_erm` enabled but the ERM CFL violated at the cortex
     ``dt_cfl``: ``attach_erm_to_simulation`` raises (gate inherited
     from erm.py); ``Cell.build`` propagates the error.
   - `Cell.with_crosslinkers` enabled but `n_xl = 0`: the Updater
     idles each tick; no-op.
   - Empty cortex (`n_filaments = 0`): caught by `resolve_h3_derived`
     ValueError; never reaches the Cell layer.

3. **Conservation invariants**
   - Particle count: ``n_cortex_actin + 2·n_xl + n_myosin_beads +
     n_lamellipodium_beads + n_fa_integrins``. Each subsystem
     contributes to ``Cell.bead_count_summary()``.
   - Bond count: cortex backbone + xlink_intra + (dynamic attach,
     variable) + future myosin rigid-body + FA bonds. Diagnostic
     ``Cell.bond_count_summary()`` reports current totals.
   - Tag conventions: tags 0..n_cortex_actin-1 = cortex; subsequent
     ranges follow the strict insertion order (xlink_head next, then
     myosin_backbone + myosin_head, then lamellipodium, then FA).
     Diagnostic ``Cell.tag_ranges()`` exposes the map.

4. **Numerical sanity**
   - All resolved parameter dataclasses (``ResolvedH3``,
     ``ResolvedCrosslinkers``, ``ResolvedERM``) carry their own
     finite-positive guards; Cell.build does NOT re-validate (would
     double-cover).
   - The HOOMD Simulation is built with the cortex.dt_cfl integrator;
     ERM CFL gate at attach time prevents stiff-spring runaway.

5. **Sign / sense**
   - Cell.build inherits all force-sign contracts from cortex.py
     (bond / angle / LJ), crosslinkers.py (Bell-Evans slip), and
     erm.py (radial inward/outward harmonic).

6. **Measurement protocol**
   - `Cell.bead_count_summary()` returns a dict with per-subsystem
     bead counts. STATIC test asserts the sum equals
     ``sim.state.N_particles``.
   - `Cell.tag_ranges()` returns a dict mapping subsystem name → (start,
     end) tag range. STATIC test asserts ranges are non-overlapping
     and contiguous.
   - `Cell.diagnostics()` returns a compact summary used by REPORT.md
     generation.

References
----------
- Brief: ``ffn_sim/docs/briefs/H3_cortex.md`` Deliverables table row
  `ffn_sim/cell/cell.py`.
- v1 archive: ``~/ActiveCellSim/acs_kb/cell/cell.py`` (NOT used in v2;
  preserved for historical reference only).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd

from ffn_sim.cortex.cortex import (
    CortexTopology,
    ResolvedH3,
    build_cortex_simulation,
    build_cortex_state,
    generate_cortex_topology,
)
from ffn_sim.cortex.crosslinkers import (
    ResolvedCrosslinkers,
    XlinkBondUpdater,
    XlinkLayout,
    build_cortex_xlink_simulation,
    extend_cortex_state_with_xlinks,
    generate_xlink_layout,
    make_xlink_updater,
)
from ffn_sim.cortex.erm import (
    ERMHarmonic,
    ResolvedERM,
    attach_erm_to_simulation,
)


@dataclass(slots=True)
class CellBuildOptions:
    """Per-build feature flags for the Cell composition.

    Each flag is wired only when the corresponding resolved dataclass
    is also passed to :meth:`Cell.build`. The options bundle is a
    small DTO so the test harness can construct a Cell with arbitrary
    combinations of subsystems enabled.
    """

    with_erm: bool = False
    with_crosslinkers: bool = False
    with_baoab: bool = True
    # Future H.3 단계 4 / H.5 / H.7 hooks (None for 단계 3).
    with_myosin: bool = False
    with_lamellipodium: bool = False
    with_fa: bool = False


@dataclass(slots=True)
class Cell:
    """v2 Cell composition: cortex + add-ons over a single HOOMD Simulation.

    Construct via :meth:`Cell.build`, not the dataclass directly — the
    factory wires the HOOMD Simulation, force computes, Updaters, and
    tag-range bookkeeping.
    """

    # Subsystem resolved configs
    p_cortex: ResolvedH3
    p_xlinks: ResolvedCrosslinkers | None
    p_erm: ResolvedERM | None

    # Build flags applied
    options: CellBuildOptions

    # HOOMD wiring
    simulation: hoomd.Simulation
    topology: CortexTopology
    xlink_layout: XlinkLayout | None
    erm_force: ERMHarmonic | None
    baoab_action: Any | None
    baoab_updater: Any | None
    xlink_action: XlinkBondUpdater | None
    xlink_updater: Any | None

    # Diagnostics
    n_cortex_actin: int
    n_xlink_heads: int

    # H.3 단계 3 hooks (None — populated by 단계 4 / H.5 / H.4 integrations).
    myosin: Any | None = None
    lamellipodium: Any | None = None
    fa: Any | None = None

    extras: dict[str, Any] = field(default_factory=dict)

    # ---------------------------------------------------------------
    # Factory
    # ---------------------------------------------------------------
    @classmethod
    def build(
        cls,
        p_cortex: ResolvedH3,
        *,
        p_xlinks: ResolvedCrosslinkers | None = None,
        p_erm: ResolvedERM | None = None,
        options: CellBuildOptions | None = None,
        device: hoomd.device.Device | None = None,
        rng: np.random.Generator | None = None,
    ) -> "Cell":
        """Compose a Cell sim from resolved subsystem configs.

        Parameters
        ----------
        p_cortex : ResolvedH3
            Required. Cortex topology + force constants.
        p_xlinks : ResolvedCrosslinkers, optional
            If provided AND options.with_crosslinkers, builds the
            crosslinker layer + D2 Bell-Evans Updater.
        p_erm : ResolvedERM, optional
            If provided AND options.with_erm, attaches the ERMHarmonic
            custom force (with CFL gate against cortex.dt_cfl).
        options : CellBuildOptions, optional
            Feature flags. Defaults to all-off (bare cortex + BAOAB).
        device, rng : HOOMD device, numpy RNG (both optional).

        Returns
        -------
        Cell
            Fully-wired composition. Access ``cell.simulation`` for the
            HOOMD Simulation handle.
        """
        opts = options or CellBuildOptions()

        # Path A: no xlinks → use the simpler cortex builder
        # Path B: xlinks → use build_cortex_xlink_simulation
        n_cortex_actin = p_cortex.n_filaments * p_cortex.beads_per_filament

        if opts.with_crosslinkers and p_xlinks is not None and p_xlinks.n_xl > 0:
            (sim, baoab_updater, baoab_action,
             xlink_updater, xlink_action, topology, xl_layout) = (
                build_cortex_xlink_simulation(
                    p_cortex, p_xlinks,
                    device=device, with_baoab=opts.with_baoab, rng=rng,
                )
            )
            n_xlink_heads = 2 * p_xlinks.n_xl
        else:
            sim, baoab_updater, baoab_action, topology, _ = (
                build_cortex_simulation(
                    p_cortex,
                    device=device,
                    with_baoab=opts.with_baoab,
                    with_crosslinkers=False,
                    rng=rng,
                )
            )
            xlink_updater = None
            xlink_action = None
            xl_layout = None
            n_xlink_heads = 0

        erm_force = None
        if opts.with_erm:
            if p_erm is None:
                raise ValueError(
                    "options.with_erm=True requires p_erm to be provided."
                )
            erm_force = attach_erm_to_simulation(
                sim, p_erm,
                actin_cortex_tag_range=(0, n_cortex_actin),
                gamma_b=p_cortex.gamma_b,
                cfl_safety_factor=p_cortex.cfl_safety_factor,
                cfl_strict=True,
            )

        return cls(
            p_cortex=p_cortex,
            p_xlinks=p_xlinks,
            p_erm=p_erm,
            options=opts,
            simulation=sim,
            topology=topology,
            xlink_layout=xl_layout,
            erm_force=erm_force,
            baoab_action=baoab_action,
            baoab_updater=baoab_updater,
            xlink_action=xlink_action,
            xlink_updater=xlink_updater,
            n_cortex_actin=n_cortex_actin,
            n_xlink_heads=n_xlink_heads,
        )

    # ---------------------------------------------------------------
    # Diagnostics
    # ---------------------------------------------------------------
    def bead_count_summary(self) -> dict[str, int]:
        """Per-subsystem bead counts (sums to ``simulation.state.N_particles``)."""
        return {
            "cortex_actin": int(self.n_cortex_actin),
            "xlink_head": int(self.n_xlink_heads),
            "myosin_backbone": 0 if self.myosin is None else self.myosin.n_backbone,
            "myosin_head": 0 if self.myosin is None else self.myosin.n_heads,
            "lamellipodium": 0 if self.lamellipodium is None else self.lamellipodium.n_beads,
            "fa_integrin": 0 if self.fa is None else self.fa.n_integrins,
        }

    def tag_ranges(self) -> dict[str, tuple[int, int]]:
        """Map subsystem name → [start, end) tag range. Non-overlapping
        and contiguous in insertion order."""
        offsets = {"cortex_actin": (0, self.n_cortex_actin)}
        cur = self.n_cortex_actin
        offsets["xlink_head"] = (cur, cur + self.n_xlink_heads)
        cur += self.n_xlink_heads
        # 단계 4 / H.5 / H.4 hooks (zero-width for now).
        offsets["myosin"] = (cur, cur)
        offsets["lamellipodium"] = (cur, cur)
        offsets["fa_integrin"] = (cur, cur)
        return offsets

    def diagnostics(self) -> dict[str, Any]:
        """Compact summary for REPORT.md generation."""
        return {
            "particle_total": int(self.simulation.state.N_particles),
            "bead_counts": self.bead_count_summary(),
            "tag_ranges": self.tag_ranges(),
            "cortex_R_cell_m": self.p_cortex.R_cell,
            "cortex_L_box_m": self.p_cortex.L_box,
            "dt_cfl_s": self.p_cortex.dt_cfl,
            "with_erm": self.options.with_erm,
            "with_crosslinkers": self.options.with_crosslinkers,
            "with_baoab": self.options.with_baoab,
            "with_myosin": self.options.with_myosin,
            "with_lamellipodium": self.options.with_lamellipodium,
            "with_fa": self.options.with_fa,
        }

    # ---------------------------------------------------------------
    # Convenience accessors
    # ---------------------------------------------------------------
    @property
    def n_particles(self) -> int:
        return int(self.simulation.state.N_particles)

    def run(self, n_steps: int) -> None:
        """Convenience pass-through to ``self.simulation.run(n_steps)``."""
        self.simulation.run(n_steps)
