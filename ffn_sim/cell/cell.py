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

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.md as md

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
from ffn_sim.cortex.myosin import (
    CortexMyosinLayout,
    MyosinStepUpdater,
    ResolvedCortexMyosin,
    cortex_myosin_attach_bin_names,
    cortex_myosin_attach_bin_rest_lengths,
    extend_state_with_cortex_myosin,
    generate_cortex_myosin_layout,
    make_cortex_myosin_updater,
    register_cortex_myosin_bond_params,
)
from ffn_sim.cell.lamellipodium import (
    ArpBranchingUpdater,
    BarbedEndElongationUpdater,
    CappingUpdater,
    LamellipodiumLayout,
    LamellipodiumState,
    ResolvedH5,
    WaveMembranePin,
    attach_lamellipodium_to_simulation,
    extend_cortex_snapshot_with_lamellipodium,
)
from ffn_sim.integrator.baoab import make_baoab_updater


def build_cortex_full_simulation(
    p_cortex: ResolvedH3,
    *,
    p_xlinks: ResolvedCrosslinkers | None = None,
    p_myosin: ResolvedCortexMyosin | None = None,
    p_lamellipodium: ResolvedH5 | None = None,
    device: hoomd.device.Device | None = None,
    with_baoab: bool = True,
    constrained: bool = False,
    constrained_dt: float | None = None,
    rng: np.random.Generator | None = None,
):
    """End-to-end builder for cortex + (optional) xlinks + myosin + lamellipodium.

    Performs the full state composition (no ERM — attach separately via
    ``attach_erm_to_simulation``):

    1. ``build_cortex_state(p_cortex, with_crosslinkers=False)`` for the
       cortex actin backbone.
    2. If ``p_xlinks`` and ``p_xlinks.n_xl > 0``:
       ``generate_xlink_layout`` + ``extend_cortex_state_with_xlinks``.
    3. If ``p_myosin`` and ``p_myosin.n_motors_per_cell > 0``:
       ``generate_cortex_myosin_layout`` +
       ``extend_state_with_cortex_myosin`` (appended AFTER xlinks so the
       motor tag block sits at the end of the tag space).
    4. If ``p_lamellipodium`` and ``p_lamellipodium.n_WAVE > 0``:
       :func:`ffn_sim.cell.lamellipodium.extend_cortex_snapshot_with_lamellipodium`
       (appended AFTER myosin so the WAVE / mother-actin tag block sits at
       the end of the tag space — H.5 단계 2 composition).
    5. HOOMD ``Simulation`` + ``bond.Harmonic`` (with all registered
       bond-type params for the four subsystems) + ``angle.Harmonic``
       (cortex-angle + optional lamel_branch_angle) + ``pair.LJ`` (WCA
       repulsive on intra-subsystem pairs; DISABLED on cross-subsystem
       pairs so binding / branching can occur).
    6. ``md.Integrator(dt=p_cortex.dt_cfl)`` with all forces.
    7. ``methods=[]`` (BAOAB Updater contract).
    8. If ``with_baoab``: BAOAB Updater attached with per-type γ_b
       inferred from ``p_cortex.gamma_b`` (and per-type Stokes drag for
       the lamellipodium beads when present).
    9. If xlinks enabled: ``XlinkBondUpdater`` attached
       (``trigger=Periodic(batch_steps)``).
    10. If myosin enabled: ``MyosinStepUpdater`` attached.
    11. If lamellipodium enabled: ``BarbedEndElongationUpdater`` +
        ``ArpBranchingUpdater`` + ``CappingUpdater`` attached, plus the
        ``WaveMembranePin`` custom force on the integrator.

    Returns
    -------
    A dict (handles) with keys:
        ``sim`` (hoomd.Simulation),
        ``topology`` (CortexTopology),
        ``xlink_layout`` (XlinkLayout | None),
        ``myosin_layout`` (CortexMyosinLayout | None),
        ``lamellipodium_layout`` (LamellipodiumLayout | None),
        ``lamellipodium_state`` (LamellipodiumState | None),
        ``wave_pin_force`` (WaveMembranePin | None),
        ``baoab_updater``, ``baoab_action``,
        ``xlink_updater``, ``xlink_action``,
        ``myosin_updater``, ``myosin_action``,
        ``elong_action``, ``elong_updater``,
        ``branch_action``, ``branch_updater``,
        ``cap_action``, ``cap_updater``,
        ``n_cortex_actin`` (int), ``n_xlink_heads`` (int),
        ``n_myosin_particles`` (int), ``n_wave_particles`` (int),
        ``n_lamellipodium_actin`` (int).
    """
    # 1. Cortex base — use the topology returned by build_cortex_state so the
    # myosin/xlink layouts see the SAME actin positions the snapshot has.
    # (Previous dual call to generate_cortex_topology + build_cortex_state
    # with a shared rng mutated state between calls → divergent topologies,
    # which made the actin-aware myosin placement land on the wrong beads.)
    cortex_snap, topology, _ = build_cortex_state(
        p_cortex, with_crosslinkers=False, rng=rng,
    )
    n_cortex_actin = p_cortex.n_filaments * p_cortex.beads_per_filament
    snap = cortex_snap

    # 2. Optional xlinks
    xlink_layout = None
    n_xlink_heads = 0
    enable_xl = (
        p_xlinks is not None and p_xlinks.n_xl > 0
    )
    if enable_xl:
        cortex_positions = topology.positions.reshape(n_cortex_actin, 3)
        cortex_filament_idx = np.repeat(
            np.arange(p_cortex.n_filaments, dtype=np.int64),
            p_cortex.beads_per_filament,
        )
        # Use xlink's own seed-derived rng (not the shared rng), so topology
        # generation choices don't reshuffle xlink placement and trip stability.
        xlink_layout = generate_xlink_layout(
            cortex_positions, cortex_filament_idx, p_xlinks,
            n_cortex_beads=n_cortex_actin,
            rng=np.random.default_rng(p_xlinks.seed),
        )
        snap = extend_cortex_state_with_xlinks(snap, xlink_layout, p_xlinks)
        n_xlink_heads = 2 * p_xlinks.n_xl

    # 3. Optional myosin
    myosin_layout = None
    n_myosin_particles = 0
    enable_myo = (
        p_myosin is not None and p_myosin.n_motors_per_cell > 0
    )
    if enable_myo:
        motor_tag_start = int(snap.particles.N)
        myosin_layout = generate_cortex_myosin_layout(
            p_myosin, p_cortex.R_cell,
            motor_tag_start=motor_tag_start,
            rng=np.random.default_rng(p_myosin.seed),
            # Actin-aware placement (KU-3.5 fix 2026-05-29, PI Option C):
            # minifilaments sit at random cortex actin beads, backbone along
            # local actin tangent, heads in tangent-plane lateral. Removes the
            # radial-offset bug that left heads ~824 nm from any actin.
            cortex_positions=topology.positions.reshape(-1, 3),
            cortex_tangents=topology.tangents,
            beads_per_filament=p_cortex.beads_per_filament,
        )
        snap = extend_state_with_cortex_myosin(snap, myosin_layout, p_myosin)
        n_myosin_particles = (
            p_myosin.n_motors_per_cell * p_myosin.n_particles_per_motor
        )

    # 4. Optional lamellipodium (H.5 단계 2 composition).  Appended AFTER
    # myosin so the WAVE / mother-actin block sits at the END of the
    # tag space, leaving the cortex / xlink / myosin tag ranges
    # unchanged.  Must run BEFORE create_state_from_snapshot because the
    # new particle / bond / angle TYPES (wave_particle, actin_lamel,
    # lamel_wave_anchor, lamel_actin_bond, lamel_branch_bond,
    # lamel_branch_angle) cannot be added once HOOMD has initialised
    # the state from a snapshot.
    lamellipodium_layout = None
    lamellipodium_state = None
    n_wave_particles = 0
    n_lamellipodium_actin = 0
    enable_lamel = (
        p_lamellipodium is not None and p_lamellipodium.n_WAVE > 0
    )
    if enable_lamel:
        snap, lamellipodium_layout, lamellipodium_state = (
            extend_cortex_snapshot_with_lamellipodium(
                snap, p_lamellipodium,
                wave_tag_start=int(snap.particles.N),
                rng=np.random.default_rng(p_lamellipodium.seed),
            )
        )
        n_wave_particles = int(p_lamellipodium.n_WAVE)
        # One mother actin seed per WAVE at construction; runtime
        # elongation / branching events grow this count via
        # sim.state.set_snapshot in the D2 Updaters.
        n_lamellipodium_actin = int(p_lamellipodium.n_WAVE)

    # 5. HOOMD Simulation + state
    sim = hoomd.Simulation(
        device=device or hoomd.device.CPU(), seed=p_cortex.seed
    )
    sim.create_state_from_snapshot(snap)

    # 5. Forces
    bond = md.bond.Harmonic()
    # constrained mode: the actin backbone stretch becomes a rigid distance
    # constraint (M-SHAKE) → zero its harmonic force; keep r0 for the
    # constraint length. Myosin/xlink bonds stay harmonic.
    bond.params["cortex-bond"] = dict(
        k=(0.0 if constrained else p_cortex.bond_k), r0=p_cortex.rest_length)
    if enable_xl:
        from ffn_sim.cortex.crosslinkers import (
            xlink_attach_bin_names, xlink_attach_bin_rest_lengths,
        )
        avg_intra_r0 = float(
            p_xlinks.alpha_fraction * p_xlinks.alpha_length
            + (1.0 - p_xlinks.alpha_fraction) * p_xlinks.filamin_length
        )
        bond.params["xlink_intra"] = dict(k=p_xlinks.k_intra, r0=avg_intra_r0)
        xl_bin_r0 = xlink_attach_bin_rest_lengths(
            p_xlinks.n_bins, p_xlinks.max_bind_dist
        )
        for i, name in enumerate(xlink_attach_bin_names(p_xlinks.n_bins)):
            bond.params[name] = dict(
                k=p_xlinks.k_attach, r0=float(xl_bin_r0[i])
            )
    if enable_myo:
        register_cortex_myosin_bond_params(bond, p_myosin)
    if enable_lamel:
        # Soft WAVE-mother anchor (10× softer than WAVE plane pin so it
        # never tightens CFL — see lamellipodium.py module docstring).
        bond.params["lamel_wave_anchor"] = dict(
            k=p_lamellipodium.k_wave_pin * 0.1, r0=p_lamellipodium.rest_length,
        )
        bond.params["lamel_actin_bond"] = dict(
            k=p_lamellipodium.bond_k, r0=p_lamellipodium.rest_length,
        )
        bond.params["lamel_branch_bond"] = dict(
            k=p_lamellipodium.bond_k, r0=p_lamellipodium.rest_length,
        )

    angle = md.angle.Harmonic()
    angle.params["cortex-angle"] = dict(k=p_cortex.angle_k, t0=p_cortex.angle_t0)
    if enable_lamel:
        angle.params["lamel_branch_angle"] = dict(
            k=p_lamellipodium.angle_branch_k, t0=p_lamellipodium.angle_branch_t0,
        )

    # Exclude bonded pairs (and angle-1-3 neighbors) from LJ.  Necessary
    # because some intra-subsystem bond lengths are SHORTER than the WCA
    # cutoff: myosin backbone segment = 700/13 ≈ 54 nm vs r_cut = 67 nm;
    # α-actinin xlink_intra = 35 nm vs r_cut = 67 nm.  Without exclusion,
    # bonded WCA neighbors would compete with the harmonic bond, giving
    # huge LJ energy at construction (verified empirically: ~1e-12 J for
    # the demo 3-way cortex → BAOAB runaway within 100 steps).
    nlist = md.nlist.Tree(buffer=0.5 * p_cortex.lj_sigma,
                          exclusions=("bond", "1-3"))
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)

    def _enable_pair(a: str, b: str, *, repulsive: bool = True) -> None:
        if (a, b) in lj.params:
            # already wired
            pass
        lj.params[(a, b)] = dict(
            epsilon=p_cortex.lj_epsilon, sigma=p_cortex.lj_sigma
        )
        lj.r_cut[(a, b)] = (
            p_cortex.lj_r_cut if (repulsive and p_cortex.lj_enabled) else 0.0
        )

    # LJ wiring rationale: WCA enabled ONLY within each subsystem
    # (actin × actin, xlink × xlink, myosin × myosin).  Inter-subsystem
    # pairs are DISABLED so that:
    # (a) myosin heads can approach cortex actin for D2 binding (no
    #     WCA blocking the close-range contact)
    # (b) xlink heads can approach cortex actin similarly
    # (c) at construction time, randomly-placed myosin backbones can sit
    #     near (or even briefly overlap) cortex actin without numerical
    #     LJ blow-up (myosin sits ABOVE actin in cortex anatomy; the
    #     intra-cortex packing is biology, not steric repulsion at this
    #     coarse-graining scale).
    # Same convention applied to myosin-xlink pairs.
    # All inter-subsystem r_cut = 0 ⇒ NO LJ contribution (HOOMD convention).
    _enable_pair("actin_cortex", "actin_cortex", repulsive=True)
    if enable_xl:
        _enable_pair("xlink_head", "xlink_head", repulsive=True)
        _enable_pair("xlink_head", "actin_cortex", repulsive=False)
    if enable_myo:
        _enable_pair("cortex_myosin_backbone", "cortex_myosin_backbone", repulsive=True)
        _enable_pair("cortex_myosin_head", "cortex_myosin_head", repulsive=True)
        _enable_pair("cortex_myosin_backbone", "cortex_myosin_head", repulsive=True)
        # All myosin × non-myosin pairs DISABLED (see rationale above).
        _enable_pair("cortex_myosin_backbone", "actin_cortex", repulsive=False)
        _enable_pair("cortex_myosin_head", "actin_cortex", repulsive=False)
        if enable_xl:
            _enable_pair("cortex_myosin_backbone", "xlink_head", repulsive=False)
            _enable_pair("cortex_myosin_head", "xlink_head", repulsive=False)
    if enable_lamel:
        # Intra-lamellipodium WCA: lamellipodial actin × itself; WAVE
        # particles are pinned to the membrane plane so we leave the
        # WAVE × WAVE pair at r_cut = 0 (no steric blockade between
        # adjacent WAVE / NPF — the membrane is the steric barrier).
        _enable_pair("actin_lamel", "actin_lamel", repulsive=True)
        _enable_pair("wave_particle", "wave_particle", repulsive=False)
        # All lamellipodium × non-lamellipodium pairs DISABLED — the
        # cortex shell and the leading-edge lamellipodium sit in
        # different spatial domains (cortex at radius R_cell, WAVE plane
        # at y = Y_max).  Suppressing WCA between subsystems matches
        # the (myosin × non-myosin = DISABLED) convention from the
        # myosin block: no steric runaway from incidental construction-
        # time proximity, and the D2 Updaters' search-radius bonds do
        # not have to fight WCA repulsion.
        _enable_pair("actin_lamel", "actin_cortex", repulsive=False)
        _enable_pair("wave_particle", "actin_cortex", repulsive=False)
        _enable_pair("actin_lamel", "wave_particle", repulsive=False)
        if enable_xl:
            _enable_pair("actin_lamel", "xlink_head", repulsive=False)
            _enable_pair("wave_particle", "xlink_head", repulsive=False)
        if enable_myo:
            _enable_pair("actin_lamel", "cortex_myosin_backbone", repulsive=False)
            _enable_pair("actin_lamel", "cortex_myosin_head", repulsive=False)
            _enable_pair("wave_particle", "cortex_myosin_backbone", repulsive=False)
            _enable_pair("wave_particle", "cortex_myosin_head", repulsive=False)

    lj.mode = "shift"

    dt_used = constrained_dt if (constrained and constrained_dt) else p_cortex.dt_cfl
    ig = md.Integrator(dt=dt_used)
    ig.forces.append(bond)
    ig.forces.append(angle)
    ig.forces.append(lj)
    # WAVE membrane plane pin — registered now so the BAOAB gamma_map
    # check at attach time sees the WAVE force in the integrator.
    wave_pin_force = None
    if enable_lamel:
        wave_pin_force = WaveMembranePin(
            p_lamellipodium,
            wave_tag_start=lamellipodium_layout.wave_tag_start,
            n_WAVE=p_lamellipodium.n_WAVE,
        )
        ig.forces.append(wave_pin_force)
    sim.operations.integrator = ig

    baoab_updater = None
    baoab_action = None
    if with_baoab:
        gamma_map: dict[str, float] = {"actin_cortex": p_cortex.gamma_b}
        if enable_xl:
            gamma_map["xlink_head"] = p_cortex.gamma_b
        if enable_myo:
            gamma_map["cortex_myosin_backbone"] = p_cortex.gamma_b
            gamma_map["cortex_myosin_head"] = p_cortex.gamma_b
        if enable_lamel:
            # WAVE + actin_lamel inherit the cortex bead Stokes drag at
            # the lamellipodium's own bead_radius (yaml-anchored, same
            # 30 nm ×40 bundle as the cortex bead by default but kept
            # subsystem-local so KU-5.x can sweep independently).
            gamma_lamel = float(
                6.0 * math.pi * 6.913e-4 * p_lamellipodium.bead_radius
            )
            gamma_map["wave_particle"] = gamma_lamel
            gamma_map["actin_lamel"] = gamma_lamel
        if constrained:
            # Rigid actin backbone (M-SHAKE + Fixman); cortex filaments are
            # contiguous N-bead blocks [f·N, (f+1)·N). Myosin/xlink/ERM beads
            # get the predictor step only (soft forces, no constraint).
            from ffn_sim.integrator.constrained_baoab import (
                make_constrained_baoab_updater,
            )
            N = p_cortex.beads_per_filament
            Ffil = p_cortex.n_filaments
            chains = [np.arange(f * N, (f + 1) * N, dtype=np.int64) for f in range(Ffil)]
            cpairs = np.asarray(topology.bond_groups, dtype=np.int64)
            baoab_action, baoab_updater = make_constrained_baoab_updater(
                kT=p_cortex.kT, gamma=gamma_map, dt=dt_used,
                constraint_pairs=cpairs,
                constraint_lengths=np.full(cpairs.shape[0], p_cortex.rest_length),
                chains=chains, seed=p_cortex.seed,
                shake_tol=1.0e-9, shake_max_iter=200,
            )
        else:
            baoab_action, baoab_updater = make_baoab_updater(
                kT=p_cortex.kT, gamma=gamma_map,
                dt=p_cortex.dt_cfl, seed=p_cortex.seed,
            )
        sim.operations.updaters.append(baoab_updater)

    xlink_updater = None
    xlink_action = None
    if enable_xl:
        xlink_action, xlink_updater = make_xlink_updater(
            p=p_xlinks, layout=xlink_layout, kT=p_cortex.kT,
            n_cortex_actin=n_cortex_actin,
        )
        sim.operations.updaters.append(xlink_updater)

    myosin_updater = None
    myosin_action = None
    if enable_myo:
        myosin_action, myosin_updater = make_cortex_myosin_updater(
            p_myo=p_myosin, layout=myosin_layout, kT=p_cortex.kT,
            n_cortex_actin=n_cortex_actin,
            # KU-3.5 segment-projection binding (Option C): pass actin
            # bond topology so the updater can match heads to segments,
            # not just bead centers (bead-only is a 500 nm grid artefact).
            cortex_bond_groups=topology.bond_groups,
        )
        sim.operations.updaters.append(myosin_updater)

    elong_action = None
    elong_updater = None
    branch_action = None
    branch_updater = None
    cap_action = None
    cap_updater = None
    if enable_lamel:
        elong_action = BarbedEndElongationUpdater(
            p=p_lamellipodium, lamel_state=lamellipodium_state,
        )
        elong_updater = hoomd.update.CustomUpdater(
            action=elong_action,
            trigger=hoomd.trigger.Periodic(p_lamellipodium.batch_steps),
        )
        branch_action = ArpBranchingUpdater(
            p=p_lamellipodium, lamel_state=lamellipodium_state,
            wave_tag_start=lamellipodium_layout.wave_tag_start,
            n_WAVE=p_lamellipodium.n_WAVE,
        )
        branch_updater = hoomd.update.CustomUpdater(
            action=branch_action,
            trigger=hoomd.trigger.Periodic(p_lamellipodium.batch_steps),
        )
        cap_action = CappingUpdater(
            p=p_lamellipodium, lamel_state=lamellipodium_state,
        )
        cap_updater = hoomd.update.CustomUpdater(
            action=cap_action,
            trigger=hoomd.trigger.Periodic(p_lamellipodium.batch_steps),
        )
        sim.operations.updaters.append(elong_updater)
        sim.operations.updaters.append(branch_updater)
        sim.operations.updaters.append(cap_updater)

    return {
        "sim": sim,
        "topology": topology,
        "xlink_layout": xlink_layout,
        "myosin_layout": myosin_layout,
        "lamellipodium_layout": lamellipodium_layout,
        "lamellipodium_state": lamellipodium_state,
        "wave_pin_force": wave_pin_force,
        "baoab_updater": baoab_updater,
        "baoab_action": baoab_action,
        "xlink_updater": xlink_updater,
        "xlink_action": xlink_action,
        "myosin_updater": myosin_updater,
        "myosin_action": myosin_action,
        "elong_action": elong_action,
        "elong_updater": elong_updater,
        "branch_action": branch_action,
        "branch_updater": branch_updater,
        "cap_action": cap_action,
        "cap_updater": cap_updater,
        "n_cortex_actin": n_cortex_actin,
        "n_xlink_heads": n_xlink_heads,
        "n_myosin_particles": n_myosin_particles,
        "n_wave_particles": n_wave_particles,
        "n_lamellipodium_actin": n_lamellipodium_actin,
    }


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
    p_myosin: ResolvedCortexMyosin | None
    p_lamellipodium: ResolvedH5 | None

    # Build flags applied
    options: CellBuildOptions

    # HOOMD wiring
    simulation: hoomd.Simulation
    topology: CortexTopology
    xlink_layout: XlinkLayout | None
    myosin_layout: CortexMyosinLayout | None
    erm_force: ERMHarmonic | None
    baoab_action: Any | None
    baoab_updater: Any | None
    xlink_action: XlinkBondUpdater | None
    xlink_updater: Any | None
    myosin_action: MyosinStepUpdater | None
    myosin_updater: Any | None
    # H.5 lamellipodium handles (None when options.with_lamellipodium=False)
    lamellipodium_layout: LamellipodiumLayout | None
    lamellipodium_state: LamellipodiumState | None
    wave_pin_force: WaveMembranePin | None
    elong_action: BarbedEndElongationUpdater | None
    elong_updater: Any | None
    branch_action: ArpBranchingUpdater | None
    branch_updater: Any | None
    cap_action: CappingUpdater | None
    cap_updater: Any | None

    # Diagnostics
    n_cortex_actin: int
    n_xlink_heads: int
    n_myosin_particles: int
    n_wave_particles: int
    n_lamellipodium_actin: int

    # Future hooks (H.5 lamellipodium populated when wired; H.4 FA when integrated).
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
        p_myosin: ResolvedCortexMyosin | None = None,
        p_lamellipodium: ResolvedH5 | None = None,
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
        p_myosin : ResolvedCortexMyosin, optional
            If provided AND options.with_myosin, builds the bipolar
            minifilament + D5 Stam-Hocky / D6 Hill / D2 Bell-Evans wiring.
        p_lamellipodium : ResolvedH5, optional
            If provided AND options.with_lamellipodium, extends the
            cortex snapshot with WAVE + mother-actin beads and attaches
            the D1 Bieling/Funk Updaters (BarbedEndElongation, ArpBranching,
            Capping) plus the WAVE membrane plane pin.  ``p_lamellipodium.dt``
            must equal ``p_cortex.dt_cfl`` (lamellipodium shares the host
            BAOAB integrator).
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
        n_cortex_actin = p_cortex.n_filaments * p_cortex.beads_per_filament

        # Pre-validate add-on options that require resolved configs so we
        # surface a clear error BEFORE any HOOMD state is created.
        if opts.with_myosin and p_myosin is None:
            raise ValueError(
                "options.with_myosin=True requires p_myosin to be provided."
            )
        if opts.with_lamellipodium and p_lamellipodium is None:
            raise ValueError(
                "options.with_lamellipodium=True requires p_lamellipodium "
                "to be provided."
            )

        # Initialise lamellipodium handles (overwritten by the unified
        # builder when with_lamellipodium=True).
        lamellipodium_layout: LamellipodiumLayout | None = None
        lamellipodium_state: LamellipodiumState | None = None
        wave_pin_force: WaveMembranePin | None = None
        elong_action: BarbedEndElongationUpdater | None = None
        elong_updater = None
        branch_action: ArpBranchingUpdater | None = None
        branch_updater = None
        cap_action: CappingUpdater | None = None
        cap_updater = None
        n_wave_particles = 0
        n_lamellipodium_actin = 0

        # Route through build_cortex_full_simulation whenever myosin OR
        # lamellipodium is requested (both need the unified snapshot-
        # construction path so the extra particle / bond / angle TYPES
        # are registered BEFORE HOOMD initialises the state — types
        # cannot be added post-create_state).  Otherwise dispatch to
        # the simpler cortex / cortex+xlink builders to preserve the
        # existing 단계 1-3 behavior verbatim.
        if opts.with_myosin or opts.with_lamellipodium:
            handles = build_cortex_full_simulation(
                p_cortex,
                p_xlinks=p_xlinks if opts.with_crosslinkers else None,
                p_myosin=p_myosin if opts.with_myosin else None,
                p_lamellipodium=p_lamellipodium if opts.with_lamellipodium else None,
                device=device, with_baoab=opts.with_baoab, rng=rng,
            )
            sim = handles["sim"]
            topology = handles["topology"]
            xl_layout = handles["xlink_layout"]
            myosin_layout = handles["myosin_layout"]
            baoab_updater = handles["baoab_updater"]
            baoab_action = handles["baoab_action"]
            xlink_updater = handles["xlink_updater"]
            xlink_action = handles["xlink_action"]
            myosin_updater = handles["myosin_updater"]
            myosin_action = handles["myosin_action"]
            n_xlink_heads = handles["n_xlink_heads"]
            n_myosin_particles = handles["n_myosin_particles"]
            lamellipodium_layout = handles["lamellipodium_layout"]
            lamellipodium_state = handles["lamellipodium_state"]
            wave_pin_force = handles["wave_pin_force"]
            elong_action = handles["elong_action"]
            elong_updater = handles["elong_updater"]
            branch_action = handles["branch_action"]
            branch_updater = handles["branch_updater"]
            cap_action = handles["cap_action"]
            cap_updater = handles["cap_updater"]
            n_wave_particles = handles["n_wave_particles"]
            n_lamellipodium_actin = handles["n_lamellipodium_actin"]
        elif opts.with_crosslinkers and p_xlinks is not None and p_xlinks.n_xl > 0:
            (sim, baoab_updater, baoab_action,
             xlink_updater, xlink_action, topology, xl_layout) = (
                build_cortex_xlink_simulation(
                    p_cortex, p_xlinks,
                    device=device, with_baoab=opts.with_baoab, rng=rng,
                )
            )
            n_xlink_heads = 2 * p_xlinks.n_xl
            myosin_layout = None
            myosin_updater = None
            myosin_action = None
            n_myosin_particles = 0
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
            myosin_layout = None
            myosin_updater = None
            myosin_action = None
            n_myosin_particles = 0

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
            p_myosin=p_myosin,
            p_lamellipodium=p_lamellipodium,
            options=opts,
            simulation=sim,
            topology=topology,
            xlink_layout=xl_layout,
            myosin_layout=myosin_layout,
            erm_force=erm_force,
            baoab_action=baoab_action,
            baoab_updater=baoab_updater,
            xlink_action=xlink_action,
            xlink_updater=xlink_updater,
            myosin_action=myosin_action,
            myosin_updater=myosin_updater,
            lamellipodium_layout=lamellipodium_layout,
            lamellipodium_state=lamellipodium_state,
            wave_pin_force=wave_pin_force,
            elong_action=elong_action,
            elong_updater=elong_updater,
            branch_action=branch_action,
            branch_updater=branch_updater,
            cap_action=cap_action,
            cap_updater=cap_updater,
            n_cortex_actin=n_cortex_actin,
            n_xlink_heads=n_xlink_heads,
            n_myosin_particles=n_myosin_particles,
            n_wave_particles=n_wave_particles,
            n_lamellipodium_actin=n_lamellipodium_actin,
        )

    # ---------------------------------------------------------------
    # Diagnostics
    # ---------------------------------------------------------------
    def bead_count_summary(self) -> dict[str, int]:
        """Per-subsystem bead counts (sums to ``simulation.state.N_particles``).

        ``lamellipodium`` reports the CONSTRUCTION-TIME bead count
        (``2 · n_WAVE`` = WAVE + one mother seed per WAVE).  The
        BarbedEndElongation / ArpBranching Updaters grow this block at
        runtime via ``sim.state.set_snapshot``; query
        ``simulation.state.N_particles`` for the current total.
        """
        if self.p_myosin is not None and self.myosin_layout is not None:
            n_myosin_backbone = (
                self.p_myosin.n_motors_per_cell * self.p_myosin.n_backbone
            )
            n_myosin_head = (
                self.p_myosin.n_motors_per_cell * 2 * self.p_myosin.n_heads_per_side
            )
        else:
            n_myosin_backbone = 0
            n_myosin_head = 0
        return {
            "cortex_actin": int(self.n_cortex_actin),
            "xlink_head": int(self.n_xlink_heads),
            "myosin_backbone": int(n_myosin_backbone),
            "myosin_head": int(n_myosin_head),
            "wave_particle": int(self.n_wave_particles),
            "lamellipodium_actin": int(self.n_lamellipodium_actin),
            "fa_integrin": 0 if self.fa is None else self.fa.n_integrins,
        }

    def tag_ranges(self) -> dict[str, tuple[int, int]]:
        """Map subsystem name → [start, end) tag range. Non-overlapping
        and contiguous in insertion order: cortex → xlink → myosin
        backbone → myosin heads → wave_particle → lamellipodium_actin → FA.

        Lamellipodium ranges reflect the CONSTRUCTION-TIME tag block;
        runtime elongation / branching events append actin_lamel beads
        at tags ``≥ lamellipodium_actin[1]`` (see
        ``LamellipodiumState.actin_next_tag``).
        """
        offsets = {"cortex_actin": (0, self.n_cortex_actin)}
        cur = self.n_cortex_actin
        offsets["xlink_head"] = (cur, cur + self.n_xlink_heads)
        cur += self.n_xlink_heads
        if self.p_myosin is not None and self.n_myosin_particles > 0:
            N = self.p_myosin.n_backbone
            H = self.p_myosin.n_heads_per_side
            M = self.p_myosin.n_motors_per_cell
            # Myosin particles are interleaved per-motor:
            # [backbone (N), +heads (H), −heads (H)] × M motors.
            # We report the WHOLE myosin block as one range, then split
            # backbone vs head sub-ranges as zero-width markers at the
            # block boundary (since they interleave per motor, not as
            # contiguous global sub-blocks).
            myo_end = cur + self.n_myosin_particles
            offsets["myosin"] = (cur, myo_end)
            cur = myo_end
        else:
            offsets["myosin"] = (cur, cur)
        offsets["wave_particle"] = (cur, cur + self.n_wave_particles)
        cur += self.n_wave_particles
        offsets["lamellipodium_actin"] = (cur, cur + self.n_lamellipodium_actin)
        cur += self.n_lamellipodium_actin
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
            "n_xlinks_realised": (
                int(self.xlink_layout.head_tag_pairs.shape[0])
                if self.xlink_layout is not None else 0
            ),
            "n_myosin_motors_realised": (
                int(self.myosin_layout.positions.shape[0])
                if self.myosin_layout is not None else 0
            ),
            "n_wave_particles": int(self.n_wave_particles),
            "n_lamellipodium_actin_realised": int(self.n_lamellipodium_actin),
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
