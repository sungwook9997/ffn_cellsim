"""H.3 Cell composition class (v2 — replaces deleted v1 `acs_kb/cell/cell.py`).

Phase 1 H.3 brief Deliverables: `ffn_sim/cell/cell.py` — new v2 `Cell`
class, HOOMD particle-group-backed, owns:

* cortex (`ffn_sim.cortex.cortex` topology + ResolvedH3)
* ERM tether (`ffn_sim.cortex.erm.ERMHarmonic` custom force)
* dynamic crosslinkers (`ffn_sim.cortex.crosslinkers.XlinkBondUpdater`)
* lamellipodium slot (H.5 deliverable — None for H.3)
* FA / adhesion slot (H.4 α restart S0/S1/S2 — wired behind `p_fa`)
* myosin slot (H.3 next-iteration deliverable — None for now)

PER CLAUDE.md the v1 `acs_kb/cell/cell.py` (single-chain Cortex-coupled)
was DELETED in the 2026-05-20 rename and MUST NOT be resurrected.
This v2 inverts the v1 architecture: HOOMD particle dynamics is the
runtime, cell composition is bookkeeping over particle groups + force
computes + Updaters.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC checks in
``ffn_sim/tests/test_cell.py``; FA S0/S1/S2 checks in
``ffn_sim/tests/test_fa_integration.py``.*

1. **Dimensional analysis**
   - All particle-count / bond-count diagnostics are integer counts.
     No unit-bearing scalars introduced here (Cell wraps already-
     resolved Cortex / ERM / Crosslinker / FA dataclasses; the underlying
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
     n_lamellipodium_beads + n_fa_integrins + n_substrate_ligands``. Each
     subsystem contributes to ``Cell.bead_count_summary()``.
   - Bond count: cortex backbone + xlink_intra + (dynamic attach,
     variable) + myosin rigid-body + FA clutch bonds. Diagnostic
     ``Cell.bond_count_summary()`` reports current totals.
   - Tag conventions: tags 0..n_cortex_actin-1 = cortex; subsequent
     ranges follow the strict insertion order (xlink_head next, then
     myosin, then lamellipodium, then FA integrin + substrate ligand).
     Diagnostic ``Cell.tag_ranges()`` exposes the map.

4. **Numerical sanity**
   - All resolved parameter dataclasses (``ResolvedH3``,
     ``ResolvedCrosslinkers``, ``ResolvedERM``, ``ResolvedH4``) carry
     their own finite-positive guards; Cell.build does NOT re-validate.
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
- FA brief: ``ffn_sim/docs/briefs/H4_FA_INTEGRATION_DESIGN.md`` (S0/S1/S2).
- v1 archive: ``~/ActiveCellSim/acs_kb/cell/cell.py`` (NOT used in v2;
  preserved for historical reference only).
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

import hoomd
import hoomd.md as md

_LOG = logging.getLogger(__name__)

from ffn_sim.cortex.cortex import (
    CortexTopology,
    ResolvedH3,
    build_cortex_simulation,
    build_cortex_state,
    build_variable_length_cortex_state,
    generate_bimodal_cortex_layout,
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
    seed_connected_mesh_xlinks,
)
from ffn_sim.cortex.erm import (
    ERMHarmonic,
    ResolvedERM,
    attach_erm_to_simulation,
)
from ffn_sim.cortex.enclosed_volume import (
    EnclosedVolumePressure,
    ResolvedEnclosedVolume,
    attach_enclosed_volume_to_simulation,
)
from ffn_sim.cell.nucleus import (  # H.9 (additive, default-off; Template-2)
    NucleusConfinement,
    ResolvedNucleus,
    attach_nucleus_confinement,
    build_nucleus_beads,
)
from ffn_sim.cell.membrane_surface import (  # H.8 (additive, default-off)
    ResolvedMembraneSurface,
    attach_membrane_surface,
)
from ffn_sim.cell.cytoplasm import (  # H.10 Tier-1 (additive, default-off; FDT-safe)
    CytoplasmTier1,
    apply_cytoplasm_drag,
)
from ffn_sim.ecm.substrate import (  # Track A (additive, default-off)
    ResolvedSubstrate,
    attach_substrate_spring,
)
from ffn_sim.cortex.turnover import (
    ActinTurnoverUpdater,
    ResolvedTurnover,
    make_turnover_updater,
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
    ResolvedH5,
    WaveMembranePin,
    ArpBranchingUpdater,
    BarbedEndElongationUpdater,
    CappingUpdater,
    LamellipodiumLayout,
    LamellipodiumState,
    attach_lamellipodium_to_simulation,
    extend_cortex_snapshot_with_lamellipodium,
)
# Plasma-membrane load (KU-5.2 leading-edge Brownian-ratchet brake; additive +
# default-off). The p_membrane=None path never touches it, so the pre-membrane
# builder stays bit-for-bit identical.
from ffn_sim.cell.membrane import (
    ResolvedMembrane,
    attach_membrane_to_simulation,
)
from ffn_sim.integrator.baoab import make_baoab_updater
from ffn_sim.integrator.baoab_device import make_baoab_updater_device

# H.4 FA integration (α restart — S0/S1/S2 vertical slice, additive +
# default-off). The p_fa=None path never touches any FA code, so the
# pre-FA builder stays bit-for-bit identical.
from ffn_sim.bridge.fa import (
    ResolvedH4,
    _build_fa_layout,
)
from ffn_sim.bridge.integrin_bonds import make_integrin_updater


# H.4 FA topology type / bond names (S0/S1/S2). Externalised so the
# substrate-anchor + clutch wiring agree on the strings; the off-path
# (p_fa=None) registers none of them.
FA_TYPE_INTEGRIN: str = "integrin"
FA_TYPE_SUBSTRATE_LIGAND: str = "ligand"
FA_BOND_INTEGRIN_LIGAND: str = "integrin_ligand"   # S1 catch bond (reuses bridge name)
FA_BOND_ACTIN_CLUTCH: str = "fa_actin_clutch"      # S2 NEW load-path bond


class SubstrateLigandPin(hoomd.custom.Action):
    """S0 substrate immobility — post-BAOAB position reset for ligand tags.

    The substrate ligand layer at z=0 is the mechanical *ground*: it must
    not move under the BAOAB heat bath. The frozen BAOAB integrator
    (``ffn_sim/integrator/baoab.py``) integrates **every** particle tag and
    exposes no integration-group-exclusion hook, so — per the H.4 design
    brief S0 ("integration-group exclusion or stiff pin … whichever the
    existing BAOAB Action supports cleanly") — this Action is appended to
    ``sim.operations.updaters`` AFTER the BAOAB Updater and resets each
    ligand tag to its construction position every step. HOOMD runs
    updaters in append order, so the reset overwrites the BAOAB
    displacement: the ligand is bit-exactly immobile (no magic-number pin
    stiffness, no CFL interaction, no edit to the frozen integrator). This
    is the literal "omit their tags from the integrated set" mechanism —
    BAOAB still *touches* the rows (it must, to satisfy its dense-tag
    contract + gamma_map coverage), but a pinned ligand's net per-step
    motion is zero.

    Parameters
    ----------
    ligand_tags : np.ndarray
        Global HOOMD tags of the substrate ligand particles.
    anchor_positions : np.ndarray, shape (len(ligand_tags), 3)
        Construction-time (z=0) positions to hold each ligand at.
    """

    def __init__(
        self,
        *,
        ligand_tags: np.ndarray,
        anchor_positions: np.ndarray,
    ) -> None:
        super().__init__()
        self._ligand_tags = np.asarray(ligand_tags, dtype=np.int64)
        self._anchor = np.asarray(anchor_positions, dtype=np.float64).copy()
        if self._anchor.shape != (self._ligand_tags.shape[0], 3):
            raise ValueError(
                "anchor_positions must have shape (n_ligand, 3); got "
                f"{self._anchor.shape} for {self._ligand_tags.shape[0]} tags."
            )
        self._anchor_by_tag = dict(
            zip(self._ligand_tags.tolist(), self._anchor)
        )
        self._sim_ref: hoomd.Simulation | None = None

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim_ref
        if sim is None:
            return
        with sim.state.cpu_local_snapshot as snap:
            pos = np.asarray(snap.particles.position)
            tag = np.asarray(snap.particles.tag)
            rows = np.flatnonzero(np.isin(tag, self._ligand_tags))
            for r in rows:
                pos[r] = self._anchor_by_tag[int(tag[r])]


@dataclass(slots=True)
class FAIntegration:
    """Bookkeeping for the S0/S1/S2 FA slice wired into a cell build.

    Returned in the builder handles dict (key ``fa_integration``) and
    stored on ``Cell.fa`` so diagnostics / tag-ranges report the real
    integrin block instead of the dead zero-width slot.
    """

    layouts: list                      # list[FALayout] with GLOBAL tags
    integrin_tag_start: int            # first integrin global tag
    n_integrins: int                   # total integrin particles
    ligand_tag_start: int              # first ligand global tag
    n_substrate_ligands: int           # total substrate ligand particles
    n_clutch_bonds: int                # static fa_actin_clutch bonds at construction
    clutch_pairs: np.ndarray           # shape (n_clutch_bonds, 2) [integrin, actin] global tags
    clutch_bin_names: list             # per-r0-bin fa_actin_clutch bond type names
    clutch_bin_r0: np.ndarray          # per-bin rest length [m]
    clutch_max_separation: float       # max integrin↔actin separation realised [m]


def _extend_snapshot_with_fa(
    snap,
    p_fa: ResolvedH4,
    *,
    cortex_positions: np.ndarray,
    n_cortex_actin: int,
    clutch_k: float,
    capture_radius: float,
):
    """Append integrin + immobile substrate-ligand particles + FA bonds (S0/S1/S2).

    Mirrors the lamellipodium snapshot-extension pattern: MUST run BEFORE
    ``create_state_from_snapshot`` because the new particle types
    (``integrin``, ``ligand``) and bond types (``integrin_ligand``,
    ``fa_actin_clutch``) cannot be added after HOOMD initialises the state.

    Tag layout — the FA block (integrins then substrate ligands) is APPENDED
    LAST, in natural append order after cortex / myosin / xlink /
    lamellipodium (mirroring the lamellipodium block). Integrins occupy the
    contiguous range ``[N0, N0 + n_int)`` and ligands ``[N0 + n_int, ...)``,
    where ``N0`` is the pre-FA particle count. NOTHING already present is
    re-tagged, so every other subsystem's absolute-tag Updater bookkeeping
    (myosin head→actin, xlink head→actin, lamellipodium) is unchanged — this
    is the S5 tag-space unification (2026-05-30) that lets FA coexist with
    those subsystems. The reused
    :class:`ffn_sim.bridge.integrin_bonds.IntegrinBondUpdater` was
    generalized (S5) to resolve per-integrin state by an explicit
    tag->local-index map (built from ``FALayout.integrin_tag_start``) and to
    read positions by global tag (a HOOMD snapshot is tag-ordered), so it no
    longer requires integrins at ``[0, n_int)``. The FA-off path never calls
    this function, so the pre-FA builder stays bit-for-bit identical.

    * S0: ligands placed at z=0 (the mechanical ground; held immobile at
      runtime by :class:`SubstrateLigandPin`).
    * S1: ``integrin_ligand`` bond TYPE registered (Pereverzev catch bonds
      are added dynamically by the IntegrinBondUpdater; none at
      construction, matching ``build_h4_state``).
    * S2: ``fa_actin_clutch`` bonds created statically between each integrin
      and its nearest cortex-actin bead within ``capture_radius`` (the
      literal load path; dynamic load-and-fail kinetics is the S5 TODO).

    Returns ``(snap, FAIntegration)``.
    """
    import gsd.hoomd

    # 1. Build the FA layout. ``_build_fa_layout`` numbers integrins
    #    [0, n_int) and ligands [n_int, n_int + n_lig) — a self-contained
    #    0-based tag space. We KEEP the integrins at [0, n_int) (the updater
    #    contract) and move everything else up by n_int.
    rng = np.random.default_rng(p_fa.seed)
    # FLOOR-BLOCKER FIX (a) [KU-3.5 v4]: the cortex shell is origin-centred
    # (south pole at z ~ -R_cell) while the substrate plane is at z ~ 0, so with
    # the default full-face FA scatter the integrin-cortex gap is ~R_cell >>
    # capture_radius and ~0 clutch bonds form (measured 9/52, n_integrin_bound=0,
    # r/r0=1.0). Build the contact footprint = cortex beads in the SOUTH CAP
    # (z <= z_min + capture_radius); _build_fa_layout then seeds each FA under one
    # cap bead so its integrins land within capture of a real cortex bead
    # (measured 52/52, all seeds; max fan-in <= 9, verified stable). No cortex
    # translation (myosin/xlink anchors untouched), all-derived geometry. The
    # standalone build_h4_state path passes no footprint -> flat-substrate
    # behaviour is bit-for-bit unchanged.
    _cortex_xyz_seed = np.asarray(cortex_positions, dtype=np.float64).reshape(
        n_cortex_actin, 3
    )
    _z_min = float(_cortex_xyz_seed[:, 2].min())
    _z_max = float(_cortex_xyz_seed[:, 2].max())
    # South-cap depth: expose ~2*n_int anchor beads so the clutch spreads thin
    # (<=2 bonds/bead via the fan-in cap below). This keeps each cortex bead's
    # bonded-exclusion count -- including HOOMD's "1-3" exclusion expansion
    # through neighbouring clutch beads + dynamic myosin/xlink binding -- within
    # the neighbour-list exclusion capacity (the "Too many bonds to process
    # exclusions" crash). Derived from the geometry (no magic number); floored at
    # capture_radius so the nearest-bead clutch still forms.
    _n_int_expected = (
        (p_fa.n_nascent_per_cell + p_fa.n_mature_per_cell) * p_fa.n_total_per_fa
    )
    _cap_depth = max(
        capture_radius,
        2.0 * (_z_max - _z_min) * (_n_int_expected / max(1, n_cortex_actin)),
    )
    _cap_mask = _cortex_xyz_seed[:, 2] <= (_z_min + _cap_depth)
    _anchor_xyz = _cortex_xyz_seed[_cap_mask]
    layouts, integrin_pos_local, ligand_pos_local = _build_fa_layout(
        p_fa, rng, contact_footprint=(_anchor_xyz, capture_radius),
    )

    N0 = int(snap.particles.N)
    n_int = int(integrin_pos_local.shape[0])
    n_lig = int(ligand_pos_local.shape[0])
    # S5 tag-space unification (2026-05-30): the FA block is APPENDED LAST,
    # in natural append order (after cortex / myosin / xlink / lamellipodium),
    # mirroring how the lamellipodium block is appended. Integrins occupy the
    # contiguous range [N0, N0 + n_int); substrate ligands [N0 + n_int, ...).
    # NOTHING that already exists is re-tagged — the cortex / myosin / xlink /
    # lamellipodium tag ranges (and their Updaters' absolute-tag bookkeeping)
    # are untouched. The reused IntegrinBondUpdater no longer requires
    # integrins at [0, n_int): it resolves per-integrin state by an explicit
    # tag->local map built from FALayout.integrin_tag_start (the S5 fix in
    # bridge/integrin_bonds.py), and reads positions by global tag (snapshot
    # is tag-ordered).
    integrin_tag_start = N0
    ligand_tag_start = N0 + n_int

    # Re-tag the layout to GLOBAL tags: _build_fa_layout numbered integrins
    # [0, n_int) and ligands [n_int, n_int + n_lig) in a self-contained 0-based
    # space; shift both blocks up by N0 so they land in the appended region.
    for lay in layouts:
        lay.integrin_tag_start = lay.integrin_tag_start + N0
        lay.ligand_tag = (lay.ligand_tag - n_int) + ligand_tag_start

    # 2. Particle-type registration (append new types after existing ones).
    old_types = list(snap.particles.types)
    new_types = list(old_types)
    for t in (FA_TYPE_INTEGRIN, FA_TYPE_SUBSTRATE_LIGAND):
        if t not in new_types:
            new_types.append(t)
    int_typeid = new_types.index(FA_TYPE_INTEGRIN)
    lig_typeid = new_types.index(FA_TYPE_SUBSTRATE_LIGAND)

    old_pos = np.asarray(snap.particles.position, dtype=np.float64)
    old_typeid = np.asarray(snap.particles.typeid, dtype=np.uint32)
    old_mass = np.asarray(snap.particles.mass, dtype=np.float64)

    # Order: [old particles (0..N0)] [integrins (N0..N0+n_int)]
    #        [ligands (N0+n_int..)]  — FA block appended LAST (S5).
    pos_all = np.concatenate(
        [old_pos, integrin_pos_local, ligand_pos_local], axis=0
    )
    typeid_all = np.concatenate(
        [
            old_typeid,
            np.full(n_int, int_typeid, dtype=np.uint32),
            np.full(n_lig, lig_typeid, dtype=np.uint32),
        ]
    )
    mass_all = np.concatenate(
        [old_mass, np.ones(n_int, dtype=np.float64), np.ones(n_lig, dtype=np.float64)]
    )

    # 3. S2 static clutch bonds: each integrin → nearest cortex-actin bead
    #    within capture_radius. cortex_positions is the (n_cortex_actin, 3)
    #    flat actin array; cortex beads keep their ORIGINAL global tags
    #    [0, n_cortex_actin) because the FA block is appended last (no shift).
    #
    #    Geometric reality (honest, documented): the substrate-ligand plane
    #    (z≈0) and the cortex shell (r≈R_cell) are SPATIALLY DISJOINT, so
    #    the nearest-actin separation can be up to capture_radius (~µm) —
    #    far larger than ℓ₀. To keep the clutch FORCE-FREE at construction
    #    (so the short BAOAB warm-up is stable WITHOUT the PI-gated
    #    equilibration prelude), each clutch bond gets r0 = its EXACT
    #    initial separation via a per-bond r0 bin (one HOOMD bond type per
    #    realised separation; type family ``fa_actin_clutch[_b{i}]``). Same
    #    per-r0-bin pattern the cortex crosslinkers use, taken to the exact
    #    one-bond-per-bin limit: ½k·0² = 0 J per clutch at construction.
    cortex_xyz = np.asarray(cortex_positions, dtype=np.float64).reshape(
        n_cortex_actin, 3
    )
    clutch_pairs_list: list[tuple[int, int]] = []
    clutch_r0_list: list[float] = []
    # Per-cortex-bead clutch fan-in cap. With contact-footprint seeding (FLOOR
    # FIX a) many integrins argmin to the SAME south-cap bead; an unbounded
    # fan-in pushes that bead's bonded-exclusion count (clutch + cortex backbone
    # + myosin/xlink) past HOOMD's neighbour-list exclusion capacity (the
    # runtime "Too many bonds to process exclusions for particle with tag ..."
    # crash, seed-dependent). This cap is a HOOMD library-capacity bound, NOT a
    # physics parameter: each integrin instead clutches its nearest WITHIN-CAPTURE
    # cortex bead that still has spare capacity. The clutch stays force-free
    # (r0 = the EXACT realised separation to whichever bead it bonds).
    _MAX_CLUTCH_FANIN = 1
    _fanin = np.zeros(n_cortex_actin, dtype=np.int64)
    for i in range(n_int):
        r_int = integrin_pos_local[i]
        d = np.linalg.norm(cortex_xyz - r_int, axis=1)
        # nearest beads first; bond the closest within-capture bead with capacity.
        for j in np.argsort(d):
            if d[j] > capture_radius:
                break  # no within-capture bead has spare capacity → unclutched
            if _fanin[j] < _MAX_CLUTCH_FANIN:
                # integrin global tag = i + N0 (appended); cortex bead j keeps
                # its original global tag j (no shift, FA block is last).
                clutch_pairs_list.append((i + N0, int(j)))
                clutch_r0_list.append(float(d[j]))
                _fanin[j] += 1
                break
    clutch_pairs = (
        np.asarray(clutch_pairs_list, dtype=np.int64).reshape(-1, 2)
        if clutch_pairs_list
        else np.empty((0, 2), dtype=np.int64)
    )
    clutch_r0_arr = np.asarray(clutch_r0_list, dtype=np.float64)
    n_clutch = int(clutch_pairs.shape[0])
    clutch_max_sep = float(clutch_r0_arr.max()) if n_clutch > 0 else 0.0

    # 4. Bond-type registration. Carry existing bonds across (SHIFTED by
    #    +n_int because all old particles moved up), then append the FA
    #    bond types:
    #    - integrin_ligand (S1): registered TYPE only; Pereverzev catch
    #      bonds populated dynamically by IntegrinBondUpdater.
    #    - fa_actin_clutch (S2): one type per clutch bond, each carrying its
    #      exact construction separation as r0. Canonical family name is
    #      ``fa_actin_clutch``; extra bonds are ``fa_actin_clutch_b1 ..``.
    old_bond_types = list(snap.bonds.types)
    old_bond_N = int(snap.bonds.N)
    # FA block is appended LAST (S5), so existing bonds keep their tags — NO
    # shift (was +n_int when integrins were prepended).
    old_bond_group = (
        np.asarray(snap.bonds.group, dtype=np.int64).reshape(old_bond_N, 2)
        if old_bond_N > 0
        else np.empty((0, 2), dtype=np.int64)
    ).astype(np.uint32)
    old_bond_typeid = (
        np.asarray(snap.bonds.typeid, dtype=np.uint32)
        if old_bond_N > 0
        else np.empty((0,), dtype=np.uint32)
    )
    new_bond_types = list(old_bond_types)
    if FA_BOND_INTEGRIN_LIGAND not in new_bond_types:
        new_bond_types.append(FA_BOND_INTEGRIN_LIGAND)

    clutch_bin_names: list[str] = []
    clutch_bin_r0_list: list[float] = []
    if FA_BOND_ACTIN_CLUTCH not in new_bond_types:
        new_bond_types.append(FA_BOND_ACTIN_CLUTCH)
    if n_clutch == 0:
        clutch_bin_names = [FA_BOND_ACTIN_CLUTCH]
        clutch_bin_r0_list = [0.0]
    else:
        for i in range(n_clutch):
            name = (
                FA_BOND_ACTIN_CLUTCH if i == 0
                else f"{FA_BOND_ACTIN_CLUTCH}_b{i}"
            )
            clutch_bin_names.append(name)
            clutch_bin_r0_list.append(float(clutch_r0_arr[i]))
            if name not in new_bond_types:
                new_bond_types.append(name)

    if n_clutch > 0:
        clutch_typeids = np.array(
            [new_bond_types.index(nm) for nm in clutch_bin_names],
            dtype=np.uint32,
        )
        bond_group_all = np.concatenate(
            [old_bond_group, clutch_pairs.astype(np.uint32)], axis=0
        )
        bond_typeid_all = np.concatenate([old_bond_typeid, clutch_typeids])
    else:
        bond_group_all = old_bond_group
        bond_typeid_all = old_bond_typeid
    clutch_bin_r0 = np.asarray(clutch_bin_r0_list, dtype=np.float64)

    # 5. Write a fresh GSD frame: particles + bonds + (passthrough) angles.
    out = gsd.hoomd.Frame()
    out.particles.N = int(pos_all.shape[0])
    out.particles.types = new_types
    out.particles.typeid = typeid_all
    out.particles.position = pos_all
    out.particles.mass = mass_all

    out.bonds.N = int(bond_group_all.shape[0])
    out.bonds.types = new_bond_types
    out.bonds.typeid = bond_typeid_all.astype(np.uint32)
    out.bonds.group = bond_group_all.astype(np.uint32)

    n_ang = int(snap.angles.N)
    if n_ang > 0:
        out.angles.N = n_ang
        out.angles.types = list(snap.angles.types)
        out.angles.typeid = np.asarray(snap.angles.typeid, dtype=np.uint32)
        # FA block appended last (S5): angle groups keep their original tags
        # (was +n_int when integrins were prepended).
        out.angles.group = (
            np.asarray(snap.angles.group, dtype=np.int64)
        ).astype(np.uint32)

    out.configuration.box = list(snap.configuration.box)

    fa_info = FAIntegration(
        layouts=layouts,
        integrin_tag_start=integrin_tag_start,
        n_integrins=n_int,
        ligand_tag_start=ligand_tag_start,
        n_substrate_ligands=n_lig,
        n_clutch_bonds=n_clutch,
        clutch_pairs=clutch_pairs,
        clutch_bin_names=clutch_bin_names,
        clutch_bin_r0=clutch_bin_r0,
        clutch_max_separation=clutch_max_sep,
    )
    return out, fa_info


@dataclass(slots=True)
class NucleusIntegration:
    """Bookkeeping for the H.9 nucleus block appended into a cell build.

    Returned in the builder handles dict (key ``nucleus_integration``). The
    nucleus is a ``nucleus_bead`` particle cloud APPENDED LAST (after cortex /
    myosin / xlink / lamellipodium / FA), so its tags occupy the contiguous
    range ``[nucleus_tag_start, nucleus_tag_start + n_beads)`` and NOTHING
    already present is re-tagged.
    """

    nucleus_tag_start: int             # first nucleus_bead global tag
    n_beads: int                       # nucleus bead count
    centroid: np.ndarray               # (3,) seed centroid [m]
    R_nuc: float                       # m   nuclear radius
    gamma_nuc: float                   # N·s/m per-bead Stokes drag


def _extend_snapshot_with_nucleus(
    snap,
    p_nuc: ResolvedNucleus,
    *,
    centroid: np.ndarray,
    gamma_nuc: float,
):
    """Append a ``nucleus_bead`` particle cloud to a snapshot (H.9, Template-2).

    Mirrors :func:`ffn_sim.cell.lamellipodium.extend_cortex_snapshot_with_lamellipodium`
    and :func:`_extend_snapshot_with_fa`: writes a fresh GSD frame extending
    ``snap`` with ``p_nuc.n_beads`` ``nucleus_bead`` particles seeded by the
    PURE :func:`ffn_sim.cell.nucleus.build_nucleus_beads` helper. MUST run
    BEFORE ``create_state_from_snapshot`` because the new ``nucleus_bead``
    particle type cannot be added once HOOMD initialises the state.

    Tag layout — the nucleus block is APPENDED LAST (after cortex / myosin /
    xlink / lamellipodium / FA), occupying ``[N0, N0 + n_beads)`` where ``N0``
    is the pre-nucleus particle count. The nucleus adds NO bonds/angles and
    references NO existing tags, so every earlier subsystem's absolute-tag
    bookkeeping (myosin head→actin, xlink head→actin, lamellipodium, FA
    clutch/integrin Updaters) is UNCHANGED — this is the CRITICAL tag-APPEND
    invariant (inserting mid-sequence would shift downstream tags and break the
    absolute-tag Updaters). Existing bonds / angles are carried over verbatim
    (no shift, because the nucleus is appended after them).

    The nucleus-off path never calls this function, so the pre-H.9 builder
    stays bit-for-bit identical.

    Returns ``(snap_new, NucleusIntegration)``.
    """
    import gsd.hoomd

    N0 = int(snap.particles.N)
    built = build_nucleus_beads(
        centroid, p_nuc.R_nuc, p_nuc.n_beads,
        gamma_nuc=gamma_nuc, type_name="nucleus_bead",
        seed=0, fill=True,
    )
    nuc_pos = built["positions"]
    n_beads = int(nuc_pos.shape[0])

    # --- Particle types: existing + nucleus_bead (appended after existing) ---
    old_types = list(snap.particles.types)
    new_types = list(old_types)
    if "nucleus_bead" not in new_types:
        new_types.append("nucleus_bead")
    nuc_typeid = new_types.index("nucleus_bead")

    old_pos = np.asarray(snap.particles.position, dtype=np.float64).reshape(N0, 3)
    old_typeid = np.asarray(snap.particles.typeid, dtype=np.uint32)
    old_mass = np.asarray(snap.particles.mass, dtype=np.float64)

    out = gsd.hoomd.Frame()
    out.particles.N = N0 + n_beads
    out.particles.types = new_types
    out.particles.position = np.concatenate([old_pos, nuc_pos], axis=0)
    out.particles.typeid = np.concatenate(
        [old_typeid, np.full(n_beads, nuc_typeid, dtype=np.uint32)]
    )
    out.particles.mass = np.concatenate(
        [old_mass, np.ones(n_beads, dtype=np.float64)]
    )

    # --- Bonds / angles: carry over verbatim (nucleus appended last → no
    #     existing tag shifts; nucleus contributes no bonds/angles). ---
    n_bond = int(snap.bonds.N)
    if n_bond > 0:
        out.bonds.N = n_bond
        out.bonds.types = list(snap.bonds.types)
        out.bonds.typeid = np.asarray(snap.bonds.typeid, dtype=np.uint32)
        out.bonds.group = np.asarray(snap.bonds.group, dtype=np.uint32)
    elif len(list(snap.bonds.types)) > 0:
        # Preserve a registered-but-empty bond-type registry.
        out.bonds.types = list(snap.bonds.types)

    n_ang = int(snap.angles.N)
    if n_ang > 0:
        out.angles.N = n_ang
        out.angles.types = list(snap.angles.types)
        out.angles.typeid = np.asarray(snap.angles.typeid, dtype=np.uint32)
        out.angles.group = np.asarray(snap.angles.group, dtype=np.uint32)
    elif len(list(snap.angles.types)) > 0:
        out.angles.types = list(snap.angles.types)

    out.configuration.box = list(snap.configuration.box)

    nuc_info = NucleusIntegration(
        nucleus_tag_start=N0,
        n_beads=n_beads,
        centroid=np.asarray(centroid, dtype=np.float64).reshape(3).copy(),
        R_nuc=float(p_nuc.R_nuc),
        gamma_nuc=float(gamma_nuc),
    )
    return out, nuc_info


def build_cortex_full_simulation(
    p_cortex: ResolvedH3,
    *,
    p_xlinks: ResolvedCrosslinkers | None = None,
    p_myosin: ResolvedCortexMyosin | None = None,
    p_lamellipodium: ResolvedH5 | None = None,
    lamellipodium_geometry: str = "flat_plane",
    lamellipodium_polarization: "np.ndarray | tuple[float, float, float] | None" = None,
    p_fa: "ResolvedH4 | None" = None,
    fa_clutch_capture_radius: float | None = None,
    fa_clutch_k: float | None = None,
    p_enclosed_volume: "ResolvedEnclosedVolume | None" = None,
    p_membrane_surface: "ResolvedMembraneSurface | None" = None,
    p_nucleus: "ResolvedNucleus | None" = None,
    nucleus_centroid: "np.ndarray | tuple[float, float, float] | None" = None,
    nucleus_gamma: float | None = None,
    p_cytoplasm: "CytoplasmTier1 | None" = None,
    p_substrate: "ResolvedSubstrate | None" = None,
    p_turnover: "ResolvedTurnover | None" = None,
    p_membrane: "ResolvedMembrane | None" = None,
    membrane_with_reaction_force: bool = False,
    device: hoomd.device.Device | None = None,
    with_baoab: bool = True,
    constrained: bool = False,
    constrained_dt: float | None = None,
    p_erm: "ResolvedERM | None" = None,
    reconcile_dt: bool = False,
    equilibrate: bool = False,
    equilibrate_steps: int = 0,
    equilibrate_softstart_steps: int = 100,
    connected_mesh: bool = False,
    faithful_connected_mesh: bool = False,
    cm_z_struct: float = 3.7,
    cm_bundle_mult: int = 2,
    cm_reach: float | None = None,
    cm_formin_fraction: float = 0.12,
    cm_L_long_mean: float = 5.0e-6,
    cm_L_short_mean: float | None = None,
    cm_arp_branch_fraction: float = 0.7,
    cm_branch_angle_deg: float = 70.0,
    rng: np.random.Generator | None = None,
):
    """End-to-end builder for cortex + (optional) xlinks + myosin + lamellipodium + FA.

    CONNECTED-MESH (fast hybrid, 2026-06-04): when ``connected_mesh=True`` the
    crosslinkers are built by :func:`seed_connected_mesh_xlinks` on the EXISTING
    fixed-N cortex topology — BRIDGE-DIFFERENT-FILAMENT (0 same-filament staples)
    + per-filament degree-cap z≈3-4 + bundling + SEEDED at construction (adhered
    baseline) — instead of the random-anchor ``generate_xlink_layout`` that left
    the prior mesh fragmented (z=1.3, giant 7 %, 56 % staples). Myosin is
    UNCHANGED (fixed-N stepping kernel intact); the bimodal variable-N backbone
    is deferred (``cortex/connected_mesh.py`` has the full faithful build). The
    crosslinker count ``n_xl`` becomes an OUTPUT of the seeding. ``cm_reach``
    defaults to the mesoscale spacing √(A_shell/n_fil). DEFAULT-OFF: when False
    the build is bit-for-bit identical to the random-anchor path.

    FAITHFUL-CONNECTED-MESH (2026-06-07): when ``faithful_connected_mesh=True``
    the cortex BASE is the bimodal variable-length layout
    (:func:`generate_bimodal_cortex_layout` — short Arp2/3 infill + long formin
    backbone + Arp2/3 70° branches, projected onto the 200 nm shell band) built
    via :func:`build_variable_length_cortex_state`, INSTEAD of the uniform
    fixed-N :func:`build_cortex_state`. This mirrors
    ``cortex/connected_mesh.py::build_connected_cortex`` (the canonical faithful
    reference) but wired into the full-cell composition. It IMPLIES the connected
    seed (treated as ``connected_mesh``-ON): the bridge mesh is seeded on the
    bimodal beads with the Arp2/3 ``branch_bonds`` counted as prior connectivity
    (``prior_bead_bonds``). The bimodal-layout params ``cm_formin_fraction`` /
    ``cm_L_long_mean`` / ``cm_L_short_mean`` / ``cm_arp_branch_fraction`` /
    ``cm_branch_angle_deg`` carry ``build_connected_cortex``'s defaults. The
    ``cm_z_struct`` / ``cm_bundle_mult`` / ``cm_reach`` seed knobs are reused.
    DEFAULT-OFF: when False the cortex base is the uniform path, bit-for-bit
    identical to today.

    Performs the full state composition (no ERM — attach separately via
    ``attach_erm_to_simulation``):

    1. ``build_cortex_state(p_cortex, with_crosslinkers=False)`` for the
       cortex actin backbone.
    2. If ``p_xlinks`` and ``p_xlinks.n_xl > 0``:
       ``generate_xlink_layout`` + ``extend_cortex_state_with_xlinks``.
    3. If ``p_myosin`` and ``p_myosin.n_motors_per_cell > 0``:
       ``generate_cortex_myosin_layout`` +
       ``extend_state_with_cortex_myosin`` (appended AFTER xlinks).
    4. If ``p_lamellipodium`` and ``p_lamellipodium.n_WAVE > 0``:
       :func:`extend_cortex_snapshot_with_lamellipodium` (appended AFTER
       myosin).
    4b. If ``p_fa`` (H.4 α restart S0/S1/S2):
       :func:`_extend_snapshot_with_fa` appends integrin + immobile
       substrate-ligand particles + the fa_actin_clutch load-path bonds
       (appended LAST so all earlier tag ranges are unchanged).
    5. HOOMD ``Simulation`` + ``bond.Harmonic`` + ``angle.Harmonic`` +
       ``pair.LJ`` (WCA repulsive intra-subsystem; cross-subsystem off).
    6. ``md.Integrator(dt=p_cortex.dt_cfl)`` with all forces.
    7. ``methods=[]`` (BAOAB Updater contract).
    8. If ``with_baoab``: BAOAB Updater with per-type γ_b.
    9-11. Subsystem Updaters (xlink / myosin / lamellipodium D2 batched).
    12. If FA on: ``IntegrinBondUpdater`` (S1) + ``SubstrateLigandPin``
        (S0), attached LAST so the pin reset runs after the BAOAB step.

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

    When ``p_fa`` is provided (H.4 α restart S0/S1/S2) the dict also
    carries: ``fa_integration`` (FAIntegration), ``fa_layout``
    (list[FALayout]), ``integrin_action`` / ``integrin_updater`` (S1
    Pereverzev catch), ``ligand_pin_action`` / ``ligand_pin_updater`` (S0
    substrate immobility), ``n_fa_integrins``, ``n_substrate_ligands``,
    ``n_fa_clutch_bonds`` (int). When ``p_fa`` is None these are None / 0
    and the build is bit-for-bit identical to the pre-FA builder.
    ``fa_clutch_capture_radius`` / ``fa_clutch_k`` override the S2 clutch
    geometry / stiffness (default: ``p_fa.capture_radius_R_FA`` /
    ``p_fa.k_int_bare``) — used to bridge the construction substrate↔cortex
    gap without touching the governed resolve_h4 defaults.
    """
    # 1. Cortex base.
    # FAITHFUL path: bimodal variable-length layout (Arp2/3 short infill + formin
    # long backbone + 70° branches), mirroring connected_mesh.build_connected_cortex.
    # The faithful base IMPLIES the connected seed (treat connected_mesh as ON).
    if faithful_connected_mesh:
        if rng is None:
            rng = np.random.default_rng(p_cortex.seed)
        topology = generate_bimodal_cortex_layout(
            p_cortex,
            formin_fraction=cm_formin_fraction,
            L_long_mean=cm_L_long_mean,
            L_short_mean=cm_L_short_mean,
            arp_branch_fraction=cm_arp_branch_fraction,
            branch_angle_deg=cm_branch_angle_deg,
            n_filaments=p_cortex.n_filaments,
            project_to_shell=True,
            rng=rng,
        )
        n_cortex_actin = int(topology.positions_flat.shape[0])
        cortex_snap = build_variable_length_cortex_state(p_cortex, topology)
        snap = cortex_snap
        connected_mesh = True  # faithful base requires the connected seed
    else:
        # UNIFORM path — use the topology returned by build_cortex_state so the
        # myosin/xlink layouts see the SAME actin positions the snapshot has.
        cortex_snap, topology, _ = build_cortex_state(
            p_cortex, with_crosslinkers=False, rng=rng,
        )
        n_cortex_actin = p_cortex.n_filaments * p_cortex.beads_per_filament
        snap = cortex_snap

    # Unified cortex bead-position view (variable-length layout exposes
    # ``positions_flat``; the fixed-N CortexTopology exposes ``positions`` which
    # reshapes to the same (n_cortex_actin, 3) for the uniform path → bit-identical).
    cortex_xyz = (
        topology.positions_flat
        if hasattr(topology, "positions_flat")
        else topology.positions.reshape(-1, 3)
    )

    # 2. Optional xlinks
    xlink_layout = None
    n_xlink_heads = 0
    _cm_seed = None                  # ConnectedMeshSeed (connected_mesh path)
    _cm_filament_idx = None          # per-actin-bead filament index (bridge rule)
    enable_xl = (
        p_xlinks is not None and p_xlinks.n_xl > 0
    )
    if enable_xl:
        cortex_positions = cortex_xyz
        # Per-bead filament index for the bridge-different-filament rule. The
        # bimodal (faithful) layout has unequal per-filament bead counts, so use
        # its own ``filament_idx``; the uniform path uses the contiguous repeat.
        if faithful_connected_mesh:
            cortex_filament_idx = topology.filament_idx
        else:
            cortex_filament_idx = np.repeat(
                np.arange(p_cortex.n_filaments, dtype=np.int64),
                p_cortex.beads_per_filament,
            )
        if connected_mesh:
            # Seed the CONNECTED percolated mesh on the cortex beads
            # (bridge-different-filament + bundling + adhered seed). On the
            # faithful base the Arp2/3 70° branch bonds count as PRIOR
            # connectivity (prior_bead_bonds) so the total coordination
            # (branches + crosslinks) lands at z_struct — mirroring
            # connected_mesh.build_connected_cortex. On the uniform fast-hybrid
            # base there are no branch bonds, so prior_bead_bonds is None.
            from dataclasses import replace as _replace
            import math as _math
            _reach = cm_reach if cm_reach is not None else _math.sqrt(
                4.0 * _math.pi * p_cortex.R_cell ** 2 / p_cortex.n_filaments
            )
            p_xlinks = _replace(p_xlinks, max_bind_dist=_reach)
            _prior_bonds = (
                topology.branch_bonds if faithful_connected_mesh else None
            )
            _cm_seed = seed_connected_mesh_xlinks(
                cortex_positions, cortex_filament_idx, p_cortex.n_filaments,
                p_xlinks, z_struct=cm_z_struct, bundle_mult=cm_bundle_mult,
                reach=_reach, R_cell=p_cortex.R_cell,
                n_cortex_beads=n_cortex_actin,
                prior_bead_bonds=_prior_bonds,
                rng=np.random.default_rng(p_xlinks.seed),
            )
            xlink_layout = _cm_seed.layout
            _cm_filament_idx = cortex_filament_idx
            snap = extend_cortex_state_with_xlinks(
                snap, xlink_layout, p_xlinks,
                seeded_attach=_cm_seed.seeded_attach,
            )
            n_xlink_heads = 2 * _cm_seed.n_xl
        else:
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
        # FAITHFUL + grip-walk is not yet supported: the grip-walk WALKING
        # geometry (myosin._tag_to_fil_pos / _bead_tag) hard-assumes the fixed-N
        # ``fil·beads_per_filament + pos`` bead-tag map, which mis-maps the
        # variable-length bimodal cortex (unequal per-filament bead counts). The
        # default ``binned_r0`` stepping mode does NOT use that map (its
        # actin-aware seeding fallback already consumes the per-bead
        # cortex_filament_idx), so it works on the bimodal base. grip_walk on the
        # variable-length bimodal cortex is now supported: the bead-tag↔(fil,pos)
        # map below is built from topology.filament_starts (2026-06-07 port).
        motor_tag_start = int(snap.particles.N)
        # Per-bead → filament index for the myosin actin-aware seeding. VARIABLE-
        # LENGTH-aware: a bimodal cortex (VariableLengthCortexLayout) has unequal
        # per-filament bead counts, so the uniform // beads_per_filament map is
        # wrong; build the index from n_beads_per_filament when present.
        if hasattr(topology, "n_beads_per_filament"):
            _myo_fil_idx = np.repeat(
                np.arange(len(topology.n_beads_per_filament), dtype=np.int64),
                topology.n_beads_per_filament,
            )
        else:
            _myo_fil_idx = np.repeat(
                np.arange(p_cortex.n_filaments, dtype=np.int64),
                p_cortex.beads_per_filament,
            )
        myosin_layout = generate_cortex_myosin_layout(
            p_myosin, p_cortex.R_cell,
            motor_tag_start=motor_tag_start,
            rng=np.random.default_rng(p_myosin.seed),
            cortex_positions=cortex_xyz,
            cortex_tangents=topology.tangents,
            beads_per_filament=p_cortex.beads_per_filament,
            cortex_filament_idx=_myo_fil_idx,
        )
        snap = extend_state_with_cortex_myosin(snap, myosin_layout, p_myosin)
        n_myosin_particles = (
            p_myosin.n_motors_per_cell * p_myosin.n_particles_per_motor
        )

    # 4. Optional lamellipodium (H.5 단계 2 composition).  Appended AFTER
    # myosin so the WAVE / mother-actin block sits at the END of the
    # tag space, leaving the cortex / xlink / myosin tag ranges
    # unchanged.  Must run BEFORE create_state_from_snapshot because the
    # new particle / bond / angle TYPES cannot be added once HOOMD has
    # initialised the state from a snapshot.
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
                geometry=lamellipodium_geometry,
                R_cell=(p_cortex.R_cell if lamellipodium_geometry != "flat_plane" else None),
                polarization=lamellipodium_polarization,
            )
        )
        n_wave_particles = int(p_lamellipodium.n_WAVE)
        n_lamellipodium_actin = int(p_lamellipodium.n_WAVE)

    # 4b. Optional FA (H.4 α restart — S0/S1/S2 + S5 tag-space unification).
    # _extend_snapshot_with_fa APPENDS the integrin + substrate-ligand block
    # LAST (natural append order, after cortex / myosin / xlink /
    # lamellipodium), leaving every earlier tag range unchanged. The reused
    # IntegrinBondUpdater was generalized (S5, bridge/integrin_bonds.py) to
    # resolve per-integrin state by an explicit tag->local map instead of
    # assuming integrins at [0, n_int), so FA now coexists with the other
    # subsystems. MUST run BEFORE create_state_from_snapshot (new particle /
    # bond types cannot be added post-create). ADDITIVE + DEFAULT-OFF: when
    # p_fa is None nothing here executes and the snapshot / forces / updaters
    # are bit-for-bit identical to the pre-FA builder.
    fa_integration = None
    enable_fa = p_fa is not None
    # S5 tag-space unification (2026-05-30): FA now coexists with myosin /
    # xlink / lamellipodium. The FA block is appended LAST (see
    # _extend_snapshot_with_fa), so NO existing particle is re-tagged and the
    # myosin / xlink / lamellipodium Updaters' absolute cortex-actin tag
    # bookkeeping stays valid. The reused IntegrinBondUpdater was generalized
    # to resolve per-integrin state by an explicit tag->local map (so
    # integrins need not be at [0, n_int)). The previous NotImplementedError
    # guard for FA + (myosin | xlink | lamellipodium) is therefore removed.
    if enable_fa:
        # Clutch geometry: each integrin → nearest cortex actin bead within
        # the capture radius (default = KU-anchored Plan H.4
        # capture_radius_R_FA = 1.5 μm). Clutch stiffness reuses k_int_bare
        # (KU-2.7, 1 pN/nm) — the actin-side half of the same molecular
        # clutch.
        #
        # HONEST geometric finding: at the real cell geometry the substrate
        # plane (z≈0) and the cortex shell (r≈R_cell ~10 μm) are spatially
        # DISJOINT, so the nearest integrin↔cortex-actin separation is
        # several μm — larger than capture_radius_R_FA. With the physical
        # radius ZERO clutch bonds form at construction (the literal load
        # path needs the cell brought onto the substrate, i.e. the PI-gated
        # equilibration / cell-positioning of S5/B2). The
        # ``fa_clutch_capture_radius`` / ``fa_clutch_k`` overrides let a
        # caller / test bridge the construction gap and exercise the real
        # clutch topology without changing the governed resolve_h4 default.
        # Dynamic load-and-fail clutch kinetics is S5 (TODO).
        snap, fa_integration = _extend_snapshot_with_fa(
            snap, p_fa,
            cortex_positions=cortex_xyz,
            n_cortex_actin=n_cortex_actin,
            clutch_k=(
                fa_clutch_k if fa_clutch_k is not None else p_fa.k_int_bare
            ),
            capture_radius=(
                fa_clutch_capture_radius
                if fa_clutch_capture_radius is not None
                else p_fa.capture_radius_R_FA
            ),
        )

    # 4c. Optional nucleus (H.9 — KU-3.B2; EXTEND Template-2). The
    # nucleus_bead cloud is APPENDED LAST (after cortex / myosin / xlink /
    # lamellipodium / FA), leaving every earlier tag range UNCHANGED — the
    # CRITICAL tag-APPEND invariant (inserting mid-sequence would shift
    # downstream tags and break the absolute-tag myosin/xlink/FA Updaters).
    # MUST run BEFORE create_state_from_snapshot (the nucleus_bead particle
    # type cannot be added post-create). ADDITIVE + DEFAULT-OFF: when
    # p_nucleus is None nothing here executes and the snapshot is bit-for-bit
    # identical to the pre-H.9 builder. The seed centroid defaults to the
    # cortex shell centroid (origin-centred per build_cortex_state); the host
    # may override via nucleus_centroid.
    nucleus_integration = None
    enable_nucleus = p_nucleus is not None
    if enable_nucleus:
        if nucleus_centroid is None:
            nuc_centroid = cortex_xyz.mean(axis=0)
        else:
            nuc_centroid = np.asarray(
                nucleus_centroid, dtype=np.float64
            ).reshape(3)
        gamma_nuc = (
            nucleus_gamma if nucleus_gamma is not None else p_cortex.gamma_b
        )
        snap, nucleus_integration = _extend_snapshot_with_nucleus(
            snap, p_nucleus, centroid=nuc_centroid, gamma_nuc=gamma_nuc,
        )

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
    # FAITHFUL path: Arp2/3 70° branch link (mother↔daughter), emitted as the
    # ``arp_branch`` bond type by build_variable_length_cortex_state when the
    # bimodal layout has branches. Same stiffness/rest-length as the backbone
    # bond (mirrors connected_mesh.build_connected_cortex). It is NOT an M-SHAKE
    # constraint pair (only cortex-bond backbone bonds are), so it stays harmonic
    # at full bond_k even in constrained mode. Register only when present in the
    # snapshot, else md.bond.Harmonic would demand params for an absent type.
    if faithful_connected_mesh and int(topology.branch_bonds.shape[0]) > 0:
        bond.params["arp_branch"] = dict(
            k=p_cortex.bond_k, r0=p_cortex.rest_length)
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
    if enable_fa:
        # S1 integrin-ligand catch bond: k = k_int_bare, r0 = integrin_r0.
        # No bonds at construction (IntegrinBondUpdater adds them); the type
        # params must still be registered so the bond force is defined.
        bond.params[FA_BOND_INTEGRIN_LIGAND] = dict(
            k=p_fa.k_int_bare, r0=p_fa.integrin_r0,
        )
        # S2 fa_actin_clutch: per-clutch-bond exact-r0 types (force-free at
        # construction). k reuses the integrin bare clutch stiffness.
        clutch_k = (
            fa_clutch_k if fa_clutch_k is not None else p_fa.k_int_bare
        )
        for name, r0 in zip(
            fa_integration.clutch_bin_names, fa_integration.clutch_bin_r0,
        ):
            bond.params[name] = dict(k=clutch_k, r0=float(r0))

    angle = md.angle.Harmonic()
    angle.params["cortex-angle"] = dict(k=p_cortex.angle_k, t0=p_cortex.angle_t0)
    # FAITHFUL path: Arp2/3 dendritic 70° branch angle, emitted as the
    # ``arp_branch_angle`` angle type by build_variable_length_cortex_state when
    # the bimodal layout has branch angles. Same bending stiffness as the
    # backbone angle, held at branch_angle_deg (mirrors connected_mesh.
    # build_connected_cortex). Register only when present in the snapshot.
    if faithful_connected_mesh and int(topology.branch_angles.shape[0]) > 0:
        angle.params["arp_branch_angle"] = dict(
            k=p_cortex.angle_k, t0=math.radians(cm_branch_angle_deg),
        )
    if enable_lamel:
        angle.params["lamel_branch_angle"] = dict(
            k=p_lamellipodium.angle_branch_k, t0=p_lamellipodium.angle_branch_t0,
        )

    # Exclude bonded pairs (and angle-1-3 neighbors) from LJ.  Necessary
    # because some intra-subsystem bond lengths are SHORTER than the WCA
    # cutoff: myosin backbone segment = 700/13 ≈ 54 nm vs r_cut = 67 nm;
    # α-actinin xlink_intra = 35 nm vs r_cut = 67 nm.  Without exclusion,
    # bonded WCA neighbors would compete with the harmonic bond, giving
    # huge LJ energy at construction.
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

    # LJ wiring rationale: WCA enabled ONLY within each subsystem.
    # Inter-subsystem pairs are DISABLED so binding / branching / clutch
    # contact can occur and randomly-placed construction-time positions do
    # not blow up. All inter-subsystem r_cut = 0 ⇒ NO LJ (HOOMD convention).
    _enable_pair("actin_cortex", "actin_cortex", repulsive=True)
    if enable_xl:
        # CONNECTED-MESH: head↔head WCA OFF — bundling places many crosslinker
        # heads near the same bridge anchors; a hard core would diverge at
        # construction for coincident heads (their excluded volume is negligible
        # for mechanics, the intra+attach bonds constrain them).
        _enable_pair("xlink_head", "xlink_head",
                     repulsive=not connected_mesh)
        _enable_pair("xlink_head", "actin_cortex", repulsive=False)
    if enable_myo:
        _enable_pair("cortex_myosin_backbone", "cortex_myosin_backbone", repulsive=True)
        _enable_pair("cortex_myosin_head", "cortex_myosin_head", repulsive=True)
        _enable_pair("cortex_myosin_backbone", "cortex_myosin_head", repulsive=True)
        _enable_pair("cortex_myosin_backbone", "actin_cortex", repulsive=False)
        _enable_pair("cortex_myosin_head", "actin_cortex", repulsive=False)
        if enable_xl:
            _enable_pair("cortex_myosin_backbone", "xlink_head", repulsive=False)
            _enable_pair("cortex_myosin_head", "xlink_head", repulsive=False)
    if enable_lamel:
        _enable_pair("actin_lamel", "actin_lamel", repulsive=True)
        _enable_pair("wave_particle", "wave_particle", repulsive=False)
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
    if enable_fa:
        # FA WCA: intra-FA steric repulsion only (integrin × integrin,
        # ligand × ligand). All FA × non-FA pairs AND integrin × ligand are
        # DISABLED (r_cut = 0), mirroring the myosin / lamellipodium
        # convention: the fa_actin_clutch bonded integrin↔actin pair is
        # already excluded from WCA by exclusions=("bond","1-3"), and
        # disabling cross-subsystem WCA prevents construction-time LJ
        # runaway from the spatially-disjoint substrate (z≈0) vs cortex
        # (r≈R_cell) layers. integrin↔ligand is left WCA-off so the
        # IntegrinBondUpdater's dynamic catch bond can close the gap to z=0
        # without fighting steric repulsion.
        _enable_pair(FA_TYPE_INTEGRIN, FA_TYPE_INTEGRIN, repulsive=True)
        _enable_pair(
            FA_TYPE_SUBSTRATE_LIGAND, FA_TYPE_SUBSTRATE_LIGAND, repulsive=True
        )
        _enable_pair(FA_TYPE_INTEGRIN, FA_TYPE_SUBSTRATE_LIGAND, repulsive=False)
        _enable_pair(FA_TYPE_INTEGRIN, "actin_cortex", repulsive=False)
        _enable_pair(FA_TYPE_SUBSTRATE_LIGAND, "actin_cortex", repulsive=False)
        if enable_xl:
            _enable_pair(FA_TYPE_INTEGRIN, "xlink_head", repulsive=False)
            _enable_pair(FA_TYPE_SUBSTRATE_LIGAND, "xlink_head", repulsive=False)
        if enable_myo:
            _enable_pair(FA_TYPE_INTEGRIN, "cortex_myosin_backbone", repulsive=False)
            _enable_pair(FA_TYPE_INTEGRIN, "cortex_myosin_head", repulsive=False)
            _enable_pair(FA_TYPE_SUBSTRATE_LIGAND, "cortex_myosin_backbone", repulsive=False)
            _enable_pair(FA_TYPE_SUBSTRATE_LIGAND, "cortex_myosin_head", repulsive=False)
        if enable_lamel:
            _enable_pair(FA_TYPE_INTEGRIN, "actin_lamel", repulsive=False)
            _enable_pair(FA_TYPE_INTEGRIN, "wave_particle", repulsive=False)
            _enable_pair(FA_TYPE_SUBSTRATE_LIGAND, "actin_lamel", repulsive=False)
            _enable_pair(FA_TYPE_SUBSTRATE_LIGAND, "wave_particle", repulsive=False)

    if enable_nucleus:
        # H.9 nucleus WCA: the nucleus is an INTERNAL bead cloud held by its
        # own NucleusConfinement custom force, NOT a steric body. ALL
        # nucleus_bead pairs (including nucleus×nucleus) are DISABLED
        # (r_cut = 0), mirroring the cross-subsystem disabled convention: the
        # confinement law alone shapes the cloud, and registering every
        # type-pair is required because md.pair.LJ demands params for EVERY
        # present type pair. Cross pairs are r_cut=0 so the nucleus cloud
        # (seeded near the cell centroid, r ≈ 0) does not sterically fight the
        # cortex shell at r ≈ R_cell.
        _enable_pair("nucleus_bead", "nucleus_bead", repulsive=False)
        _enable_pair("nucleus_bead", "actin_cortex", repulsive=False)
        if enable_xl:
            _enable_pair("nucleus_bead", "xlink_head", repulsive=False)
        if enable_myo:
            _enable_pair("nucleus_bead", "cortex_myosin_backbone", repulsive=False)
            _enable_pair("nucleus_bead", "cortex_myosin_head", repulsive=False)
        if enable_lamel:
            _enable_pair("nucleus_bead", "actin_lamel", repulsive=False)
            _enable_pair("nucleus_bead", "wave_particle", repulsive=False)
        if enable_fa:
            _enable_pair("nucleus_bead", FA_TYPE_INTEGRIN, repulsive=False)
            _enable_pair("nucleus_bead", FA_TYPE_SUBSTRATE_LIGAND, repulsive=False)

    lj.mode = "shift"

    dt_used = constrained_dt if (constrained and constrained_dt) else p_cortex.dt_cfl

    # B1 — global CFL dt reconciliation (additive, default-off).
    # When reconcile_dt is True, compute the binding CFL timestep across ALL
    # active stiff subsystems (cortex + xlinks/myosin/fa/erm/...) and lower
    # the integrator dt to min(requested dt_used, dt_min). This catches the
    # case where a stiff optional bond (e.g. FA k_int_bare) has a relaxation
    # timescale shorter than any cortex timescale, which would otherwise
    # violate CFL at the cortex dt and trip the int32 image guard. When False
    # (default) dt_used is unchanged and the builder is bit-for-bit identical
    # to the pre-B1 version. Does NOT touch the frozen integrator — it only
    # computes/lowers the dt scalar passed to md.Integrator + the BAOAB
    # Action.
    cfl_result = None
    if reconcile_dt:
        from ffn_sim.cell.dt_reconcile import compute_global_cfl_dt
        cfl_result = compute_global_cfl_dt(
            p_cortex,
            p_xlinks=p_xlinks, p_myosin=p_myosin,
            p_lamellipodium=p_lamellipodium, p_fa=p_fa, p_erm=p_erm,
            p_enclosed_volume=p_enclosed_volume, p_turnover=p_turnover,
            p_membrane=p_membrane,
            constrained=constrained,
            constrained_dt=(dt_used if constrained else None),
        )
        if cfl_result.dt_min < dt_used:
            _LOG.info(
                "B1 dt reconciliation: lowering integrator dt from %.4e s to "
                "%.4e s.\n%s", dt_used, cfl_result.dt_min,
                cfl_result.format_breakdown(),
            )
            dt_used = cfl_result.dt_min
        else:
            _LOG.info(
                "B1 dt reconciliation: requested dt %.4e s already within the "
                "global CFL bound %.4e s; unchanged.\n%s",
                dt_used, cfl_result.dt_min, cfl_result.format_breakdown(),
            )

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

    # Optional enclosed-volume / intracellular-pressure term (KU-3.1).
    # ADDITIVE + DEFAULT-OFF: when p_enclosed_volume is None this whole
    # block is skipped and the builder is bit-for-bit identical to the
    # pre-enclosed-volume version. When provided, attach the
    # EnclosedVolumePressure custom force over the cortex-actin shell tags
    # [0, n_cortex_actin). The force participates in HOOMD net_force (read
    # by the BAOAB Updater). Its effective per-bead stiffness is far softer
    # than bond/angle/xl, so the CFL gate (parity with erm.py) never fires
    # at the cortex dt_cfl — but it is enforced for safety.
    enclosed_volume_force = None
    if p_enclosed_volume is not None:
        enclosed_volume_force = attach_enclosed_volume_to_simulation(
            sim, p_enclosed_volume,
            shell_tag_range=(0, n_cortex_actin),
            n_shell=n_cortex_actin,
            gamma_b=p_cortex.gamma_b,
            cfl_safety_factor=p_cortex.cfl_safety_factor,
        )

    # H.8 membrane surface (ADDITIVE + DEFAULT-OFF): when p_membrane_surface is
    # None this block is skipped and the builder is bit-for-bit identical. When
    # provided, attach the MembraneSurfaceTension custom force over the same
    # cortex-actin shell tags [0, n_cortex_actin) (Template-1, mirrors
    # enclosed_volume). Supplies whole-cell surface tension toward the composite
    # KU-3.5 (gamma_total = gamma_membrane + gamma_cortex). Soft k_eff -> the CFL
    # gate (parity with erm/enclosed_volume) essentially never fires.
    membrane_surface_force = None
    if p_membrane_surface is not None:
        membrane_surface_force = attach_membrane_surface(
            sim, p_membrane_surface,
            shell_tag_range=(0, n_cortex_actin),
            n_shell=n_cortex_actin,
            gamma_b=p_cortex.gamma_b,
            cfl_safety_factor=p_cortex.cfl_safety_factor,
        )

    # H.9 nucleus confinement (ADDITIVE + DEFAULT-OFF; Template-2 step 3).
    # When p_nucleus is None this block is skipped and the builder is
    # bit-for-bit identical. When provided, attach the NucleusConfinement
    # custom force over the appended nucleus_bead tags [nucleus_tag_start,
    # nucleus_tag_start + n_beads) — attached AFTER the integrator is set
    # (the helper requires it) and BEFORE the BAOAB Updater (parity with
    # enclosed_volume / membrane_surface). The two-regime radial force
    # participates in HOOMD net_force (read by the BAOAB Updater). CFL gate
    # (k_hi = k_chrom + k_lamin; τ = γ_nuc / k_hi) is enforced strict at the
    # host integrator dt for parity with erm / enclosed_volume — a violation
    # surfaces to the host rather than silently using an unstable dt.
    nucleus_force = None
    if enable_nucleus:
        nuc_start = nucleus_integration.nucleus_tag_start
        nucleus_force = attach_nucleus_confinement(
            sim, p_nucleus,
            nucleus_tags=(nuc_start, nuc_start + nucleus_integration.n_beads),
            gamma_nuc=nucleus_integration.gamma_nuc,
            cfl_safety_factor=p_cortex.cfl_safety_factor,
            cfl_strict=True,
        )

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
            gamma_lamel = float(
                6.0 * math.pi * 6.913e-4 * p_lamellipodium.bead_radius
            )
            gamma_map["wave_particle"] = gamma_lamel
            gamma_map["actin_lamel"] = gamma_lamel
        if enable_fa:
            # Integrin / substrate-ligand Stokes drag from resolve_h4
            # (γ = 6π η R). The ligand drag is moot at runtime (its position
            # is pinned by SubstrateLigandPin every step), but BAOAB
            # requires every present particle type in gamma_map and asserts
            # γ > 0 finite, so both are supplied.
            gamma_map[FA_TYPE_INTEGRIN] = p_fa.gamma_integrin
            gamma_map[FA_TYPE_SUBSTRATE_LIGAND] = p_fa.gamma_ligand
        if enable_nucleus:
            # H.9 nucleus_bead Stokes drag (= 6π η R_bead, supplied by the
            # host like the cortex gamma_b). BAOAB asserts every PRESENT
            # particle type has a gamma_map entry (γ > 0 finite) or HALTS,
            # and the appended nucleus_bead cloud is now a present type.
            gamma_map["nucleus_bead"] = nucleus_integration.gamma_nuc
        # H.10 cytoplasm Tier-1 (per-type effective-viscosity drag): scale the
        # immersed types' Stokes drag by eta_eff/eta_water once gamma_map is fully
        # assembled. p_cytoplasm is None => bit-for-bit identical (water). FDT-safe:
        # BAOAB reads gamma per-type for the noise amplitude, so D=kT/gamma_b stays
        # exact under any per-type eta_eff (H.10 design; cell/cytoplasm.py).
        gamma_map = apply_cytoplasm_drag(gamma_map, p_cytoplasm)
        if constrained:
            # Rigid actin backbone (M-SHAKE + Fixman); cortex filaments are
            # contiguous bead blocks. UNIFORM: each filament is N beads at
            # [f·N, (f+1)·N). VARIABLE-LENGTH (faithful): per-filament bead
            # counts differ, so derive each chain from filament_starts +
            # n_beads_per_filament. The constraint pairs are the backbone bonds
            # only (bond_groups, NOT the Arp2/3 branch bonds); ℓ0 spacing is
            # uniform so constraint_lengths stays rest_length in both paths.
            from ffn_sim.integrator.constrained_baoab import (
                make_constrained_baoab_updater,
            )
            if hasattr(topology, "filament_starts"):
                chains = [
                    np.arange(
                        int(topology.filament_starts[f]),
                        int(topology.filament_starts[f])
                        + int(topology.n_beads_per_filament[f]),
                        dtype=np.int64,
                    )
                    for f in range(len(topology.n_beads_per_filament))
                ]
            else:
                N = p_cortex.beads_per_filament
                Ffil = p_cortex.n_filaments
                chains = [
                    np.arange(f * N, (f + 1) * N, dtype=np.int64)
                    for f in range(Ffil)
                ]
            cpairs = np.asarray(topology.bond_groups, dtype=np.int64)
            baoab_action, baoab_updater = make_constrained_baoab_updater(
                kT=p_cortex.kT, gamma=gamma_map, dt=dt_used,
                constraint_pairs=cpairs,
                constraint_lengths=np.full(cpairs.shape[0], p_cortex.rest_length),
                chains=chains, seed=p_cortex.seed,
                shake_tol=1.0e-9, shake_max_iter=200,
            )
        else:
            # GPU-main port (2026-06-07): device-resident BAOAB on GPU
            # (gpu_local_snapshot + cupy — removes the per-step device→host sync that
            # capped throughput on the unconstrained GPU path). OPT-IN (default OFF →
            # frozen cpu BAOAB, GPU path UNCHANGED = zero regression while the device
            # path is validated on the full cortex+myosin system). Enable with
            # FFN_GPU_DEVICE_BAOAB=1. The device factory is signature-identical +
            # Phase-1-validated for the solver (diffusion/equipartition/drift ~1e-13);
            # full-system GPU validation pending. The constrained path already auto-adapts.
            import os as _os
            _use_device = (
                isinstance(sim.device, hoomd.device.GPU)
                and _os.environ.get("FFN_GPU_DEVICE_BAOAB") == "1"
            )
            _baoab_factory = make_baoab_updater_device if _use_device else make_baoab_updater
            baoab_action, baoab_updater = _baoab_factory(
                kT=p_cortex.kT, gamma=gamma_map,
                dt=dt_used, seed=p_cortex.seed,
            )
        sim.operations.updaters.append(baoab_updater)

    xlink_updater = None
    xlink_action = None
    if enable_xl:
        xlink_action, xlink_updater = make_xlink_updater(
            p=p_xlinks, layout=xlink_layout, kT=p_cortex.kT,
            n_cortex_actin=n_cortex_actin,
            # CONNECTED-MESH: bridge-different-filament rule on rebind +
            # seeded-bound initial state (adhered baseline).
            actin_filament_idx=_cm_filament_idx,
            head_to_actin0=(_cm_seed.head_to_actin0
                            if _cm_seed is not None else None),
        )
        sim.operations.updaters.append(xlink_updater)

    myosin_updater = None
    myosin_action = None
    if enable_myo:
        myosin_action, myosin_updater = make_cortex_myosin_updater(
            p_myo=p_myosin, layout=myosin_layout, kT=p_cortex.kT,
            n_cortex_actin=n_cortex_actin,
            cortex_bond_groups=topology.bond_groups,
            # Grip-walk geometry (used iff p_myosin.stepping_mode=="grip_walk";
            # harmless otherwise). The bead-tag ↔ (filament, pos) map + the bead
            # spacing ℓ₀ one walked sub-bead step is measured in. UNIFORM cortex
            # uses the fixed-N map (beads_per_filament); the VARIABLE-LENGTH bimodal
            # faithful cortex passes filament_starts + n_beads_per_filament so the
            # walk maps tags correctly (2026-06-07 port).
            ell0_cortex=p_cortex.rest_length,
            cortex_beads_per_filament=p_cortex.beads_per_filament,
            cortex_filament_starts=(
                topology.filament_starts
                if hasattr(topology, "filament_starts") else None
            ),
            cortex_n_beads_per_filament=(
                topology.n_beads_per_filament
                if hasattr(topology, "n_beads_per_filament") else None
            ),
        )
        sim.operations.updaters.append(myosin_updater)

    # Optional actin turnover (cofilin severing + pointed-end re-annealing).
    # ADDITIVE + DEFAULT-OFF: when p_turnover is None this block is skipped and
    # the builder is bit-for-bit identical to the pre-turnover version. When
    # provided, attach the ActinTurnoverUpdater (D2 batched) which mutates the
    # cortex-bond + cortex-angle topology only (particle count invariant; the
    # BAOAB Action stays the sole position integrator). It is appended AFTER the
    # BAOAB updater so its bond-topology rewrite happens after the position step
    # each tick (matching the other bond-mutating subsystem updaters); it
    # touches a DISJOINT bond type (cortex-bond) from xlink / myosin / FA, so
    # update order among them does not matter.
    turnover_updater = None
    turnover_action = None
    if p_turnover is not None:
        turnover_action, turnover_updater = make_turnover_updater(
            p=p_turnover,
            rest_length=p_cortex.rest_length,
            kT=p_cortex.kT,
            bond_k=p_cortex.bond_k,
        )
        sim.operations.updaters.append(turnover_updater)

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

    # Optional plasma-membrane load (KU-5.2 leading-edge Brownian-ratchet
    # brake; mechanism audit 2026-05-30 §5 #6 / §6 #7; ffn_sim/cell/membrane.py).
    # ADDITIVE + DEFAULT-OFF: when p_membrane is None this whole block is
    # skipped and the builder is bit-for-bit identical to the pre-membrane
    # version. The membrane supplies the per-barbed-end reaction load F the
    # lamellipodium's Bell-Evans elongation/capping laws consume (today F=0 →
    # laws inert → KU-5.1 ~10,000× under). It is a MembraneLoad Action wired
    # to feed elong_action / cap_action, attached AFTER the lamellipodium
    # updaters but BEFORE FA so its load is current; HOOMD runs updaters in
    # append order. Requires the lamellipodium to be present (it drives the
    # barbed-end rates) — a membrane with no lamellipodium has nothing to
    # load, so it is a no-op guarded by enable_lamel.
    membrane_load_action = None
    membrane_updater = None
    membrane_reaction_force = None
    if p_membrane is not None and enable_lamel:
        membrane_handles = attach_membrane_to_simulation(
            sim, p_membrane,
            lamel_state=lamellipodium_state,
            barbed_load_consumers=[elong_action, cap_action],
            batch_steps=p_lamellipodium.batch_steps,
            gamma_b=p_cortex.gamma_b,
            with_reaction_force=membrane_with_reaction_force,
            cfl_safety_factor=p_cortex.cfl_safety_factor,
        )
        membrane_load_action = membrane_handles["membrane_load"]
        membrane_updater = membrane_handles["membrane_updater"]
        membrane_reaction_force = membrane_handles["membrane_reaction_force"]

    # FA Updaters (S0 + S1). Attached LAST so the SubstrateLigandPin reset
    # runs AFTER the BAOAB position step each tick (HOOMD runs updaters in
    # append order) — that is what makes the substrate ligands immobile
    # without editing the frozen integrator.
    integrin_action = None
    integrin_updater = None
    ligand_pin_action = None
    ligand_pin_updater = None
    substrate_spring_force = None
    if enable_fa:
        # S1: Pereverzev catch-slip dynamic integrin↔ligand bonds (reused
        # as-is from bridge/integrin_bonds.py). The PI-gated Pereverzev→Kong
        # migration is NOT done here — the Pereverzev catch parameters are
        # used verbatim.
        #
        # B1 dt reconciliation (in-builder half; the GLOBAL dt = min(τ)
        # reconciliation is PI-gated per the design brief §6). The
        # IntegrinBondUpdater's D2 batch-CFL contract is
        # ``integrin_batch_steps · dt · k_off_max ≤ 1e-3``, and the updater
        # checks it against ``p_fa.dt`` — the dt resolve_h4 computed in
        # ISOLATION for a standalone FA sim. But the integrated cell runs on
        # the cortex integrator's ``dt_used`` (= cortex dt_cfl ≈ 1.5e-9 s),
        # which is ~10³× SMALLER than the standalone FA dt. Feeding the
        # standalone FA dt into the updater both (a) mis-states the batch
        # window (the real per-batch wall-time is batch_steps · dt_used) and
        # (b) trips the 1e-3 ceiling at the small-demo FA scale. So pass the
        # updater a dt-reconciled copy of p_fa: dt := the host dt_used, and
        # integrin_batch_steps recomputed so the SAME physical batch window
        # (≈ resolve_h4's batch_steps · fa.dt) is preserved while honoring
        # the 1e-3 ceiling at the host dt. This is a derived, grid-invariant
        # reconciliation — no magic number, no governed-param edit, no
        # change to integrin_bonds.py or resolve_h4. (Dynamic clutch
        # kinetics + the global multi-timescale dt policy is S5/B1, gated.)
        BATCH_CFL_CEILING = 1.0e-3   # D2 contract (gate_bell_evans_batch_cfl)
        target_batch_window = p_fa.integrin_batch_steps * p_fa.dt
        reconciled_batch_steps = max(1, int(round(target_batch_window / dt_used)))
        # Clamp so events_per_batch = batch · dt_used · k_off_max ≤ ceiling.
        max_batch_at_host_dt = int(
            BATCH_CFL_CEILING / (dt_used * p_fa.bond_event_rate_max)
        )
        if max_batch_at_host_dt >= 1:
            reconciled_batch_steps = min(
                reconciled_batch_steps, max_batch_at_host_dt
            )
        p_fa_host = replace(
            p_fa, dt=dt_used, integrin_batch_steps=reconciled_batch_steps,
        )
        integrin_action, integrin_updater = make_integrin_updater(
            sim=sim, p=p_fa_host, layouts=fa_integration.layouts, bond=bond,
        )
        sim.operations.updaters.append(integrin_updater)

        # S0: pin substrate ligands to their construction z=0 positions
        # every step (post-BAOAB reset → bit-exact immobile ground).
        lig_tags = np.arange(
            fa_integration.ligand_tag_start,
            fa_integration.ligand_tag_start + fa_integration.n_substrate_ligands,
            dtype=np.int64,
        )
        read_snap = sim.state.get_snapshot()
        anchor_pos = np.asarray(
            read_snap.particles.position, dtype=np.float64
        )[lig_tags].copy()
        if p_substrate is not None:
            # Track A compliant substrate (ADDITIVE, default-off): a tunable
            # k_sub anchor spring over the ligand tags REPLACES the rigid pin
            # (mutually exclusive — substrate.py docstring). Default
            # (p_substrate is None) keeps the bit-for-bit rigid SubstrateLigandPin.
            substrate_spring_force = attach_substrate_spring(
                sim, p_substrate.k_sub, lig_tags, anchor_pos,
                p_fa.gamma_ligand,
                cfl_safety_factor=p_cortex.cfl_safety_factor,
            )
        else:
            ligand_pin_action = SubstrateLigandPin(
                ligand_tags=lig_tags, anchor_positions=anchor_pos,
            )
            ligand_pin_updater = hoomd.update.CustomUpdater(
                action=ligand_pin_action, trigger=hoomd.trigger.Periodic(1),
            )
            sim.operations.updaters.append(ligand_pin_updater)

    # Assemble the handles dict (also surfaces B1 results: cfl_result is None
    # when reconcile_dt is False; dt_used is the integrator dt actually used).
    handles = {
        "sim": sim,
        "topology": topology,
        "xlink_layout": xlink_layout,
        "myosin_layout": myosin_layout,
        "lamellipodium_layout": lamellipodium_layout,
        "lamellipodium_state": lamellipodium_state,
        "wave_pin_force": wave_pin_force,
        "enclosed_volume_force": enclosed_volume_force,
        "membrane_surface_force": membrane_surface_force,
        # H.9 nucleus (KU-3.B2) — None / 0 when p_nucleus is None.
        "nucleus_force": nucleus_force,
        "nucleus_integration": nucleus_integration,
        "n_nucleus_beads": (
            nucleus_integration.n_beads if nucleus_integration else 0
        ),
        "baoab_updater": baoab_updater,
        "baoab_action": baoab_action,
        "xlink_updater": xlink_updater,
        "xlink_action": xlink_action,
        "myosin_updater": myosin_updater,
        "myosin_action": myosin_action,
        "turnover_updater": turnover_updater,
        "turnover_action": turnover_action,
        "elong_action": elong_action,
        "elong_updater": elong_updater,
        "branch_action": branch_action,
        "branch_updater": branch_updater,
        "cap_action": cap_action,
        "cap_updater": cap_updater,
        # Plasma-membrane load (KU-5.2) — None when p_membrane is None or no
        # lamellipodium is present (additive + default-off).
        "membrane_load_action": membrane_load_action,
        "membrane_updater": membrane_updater,
        "membrane_reaction_force": membrane_reaction_force,
        "n_cortex_actin": n_cortex_actin,
        "n_xlink_heads": n_xlink_heads,
        "n_myosin_particles": n_myosin_particles,
        "n_wave_particles": n_wave_particles,
        "n_lamellipodium_actin": n_lamellipodium_actin,
        # H.4 FA (S0/S1/S2) handles — None / 0 when p_fa is None.
        "fa_integration": fa_integration,
        "fa_layout": (fa_integration.layouts if fa_integration else None),
        "integrin_action": integrin_action,
        "integrin_updater": integrin_updater,
        "ligand_pin_action": ligand_pin_action,
        "ligand_pin_updater": ligand_pin_updater,
        "substrate_spring_force": substrate_spring_force,
        "n_fa_integrins": (fa_integration.n_integrins if fa_integration else 0),
        "n_substrate_ligands": (
            fa_integration.n_substrate_ligands if fa_integration else 0
        ),
        "n_fa_clutch_bonds": (
            fa_integration.n_clutch_bonds if fa_integration else 0
        ),
        # B1 — global CFL dt reconciliation results (None / cortex dt_cfl on
        # the default off-path). dt_used is the dt actually passed to the
        # integrator + BAOAB Action.
        "cfl_result": cfl_result,
        "dt_used": dt_used,
        # B2 — equilibration prelude diagnostics. None unless the prelude ran
        # (equilibrate=True or equilibrate_steps>0). Populated just below.
        "equilibrate_diagnostics": None,
    }

    # B2 — optional equilibration prelude (additive, default-off). When
    # neither equilibrate nor equilibrate_steps is requested this block is
    # skipped entirely and the builder is bit-for-bit identical to the
    # pre-B2 version (no equilibrate_cell call, no sim.run). When requested,
    # drain construction overlaps / let the cell settle BEFORE the production
    # loop, then stash the diagnostics dict under "equilibrate_diagnostics".
    if equilibrate or equilibrate_steps > 0:
        from ffn_sim.cell.equilibration import equilibrate_cell

        handles["equilibrate_diagnostics"] = equilibrate_cell(
            handles,
            n_softstart=equilibrate_softstart_steps,
            n_baoab=equilibrate_steps,
            rest_length=p_cortex.rest_length,
            gamma_b=p_cortex.gamma_b,
        )

    return handles


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

    # A1 full-cell compartment configs (2026-06-05). None when not wired; set
    # when Cell.build forwards them to build_cortex_full_simulation. The raw
    # builder handles dict (incl. compartment force objects) is in
    # ``extras["handles"]`` for downstream measurement (e.g. gamma estimators).
    p_enclosed_volume: Any | None = None
    p_membrane_surface: Any | None = None
    p_nucleus: Any | None = None
    p_cytoplasm: Any | None = None
    p_substrate: Any | None = None
    p_turnover: Any | None = None
    p_membrane: Any | None = None

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
        lamellipodium_geometry: str = "flat_plane",
        lamellipodium_polarization: "np.ndarray | tuple[float, float, float] | None" = None,
        p_fa: "ResolvedH4 | None" = None,
        fa_clutch_capture_radius: float | None = None,
        fa_clutch_k: float | None = None,
        # --- A1 full-cell unification (2026-06-05) ---
        # Compartment + advanced run-mode params forwarded verbatim to
        # build_cortex_full_simulation so the physiological baseline (A2
        # manifest) is reachable through this single public entry point.
        # All default-None / default-off: when none are passed, routing and
        # numerics are bit-for-bit identical to the pre-A1 Cell.build.
        p_enclosed_volume: "ResolvedEnclosedVolume | None" = None,
        p_membrane_surface: "ResolvedMembraneSurface | None" = None,
        p_nucleus: "ResolvedNucleus | None" = None,
        nucleus_centroid: "np.ndarray | tuple[float, float, float] | None" = None,
        nucleus_gamma: float | None = None,
        p_cytoplasm: "CytoplasmTier1 | None" = None,
        p_substrate: "ResolvedSubstrate | None" = None,
        p_turnover: "ResolvedTurnover | None" = None,
        p_membrane: "ResolvedMembrane | None" = None,
        membrane_with_reaction_force: bool = False,
        constrained: bool = False,
        constrained_dt: float | None = None,
        reconcile_dt: bool = False,
        equilibrate: bool = False,
        equilibrate_steps: int = 0,
        equilibrate_softstart_steps: int = 100,
        connected_mesh: bool = False,
        cm_z_struct: float = 3.7,
        cm_bundle_mult: int = 2,
        cm_reach: float | None = None,
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
            the D1 Bieling/Funk Updaters plus the WAVE membrane plane pin.
        p_fa : ResolvedH4, optional
            If provided AND options.with_fa, wires the H.4 α restart
            S0/S1/S2 slice: a fixed substrate-ligand layer at z=0
            (immobile via SubstrateLigandPin), the integrin–ligand catch
            bond (IntegrinBondUpdater), and the fa_actin_clutch load-path
            bonds. ``fa_clutch_capture_radius`` / ``fa_clutch_k`` override
            the S2 clutch geometry / stiffness (see
            ``build_cortex_full_simulation``).
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
        if opts.with_fa and p_fa is None:
            raise ValueError(
                "options.with_fa=True requires p_fa (ResolvedH4) to be "
                "provided."
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
        # lamellipodium OR FA is requested (all need the unified snapshot-
        # construction path so the extra particle / bond / angle TYPES
        # are registered BEFORE HOOMD initialises the state — types cannot
        # be added post-create_state).  Otherwise dispatch to the simpler
        # cortex / cortex+xlink builders to preserve the existing 단계 1-3
        # behavior verbatim.
        #
        # FA gates on the with_fa flag AND a resolved p_fa, mirroring the
        # per-subsystem flag+config convention. When with_fa is False OR
        # p_fa is None the FA path is fully inert (additive + default-off).
        enable_fa = opts.with_fa and p_fa is not None
        fa_integration = None
        full_handles: dict[str, Any] | None = None
        # A1 (2026-06-05): the simple cortex / cortex+xlink builders cannot wire
        # compartments or advanced run-modes, so route through the unified
        # build_cortex_full_simulation whenever ANY of those is requested (in
        # addition to the original myosin/lamellipodium/FA triggers). When none
        # are set, the legacy dispatch below is preserved bit-for-bit — the
        # simpler builders register fewer particle TYPES, so collapsing them
        # unconditionally would change snapshots and break determinism.
        needs_full = (
            opts.with_myosin or opts.with_lamellipodium or enable_fa
            or any(p is not None for p in (
                p_enclosed_volume, p_membrane_surface, p_nucleus,
                p_cytoplasm, p_substrate, p_turnover, p_membrane,
            ))
            or constrained or reconcile_dt or equilibrate or connected_mesh
        )
        if needs_full:
            handles = build_cortex_full_simulation(
                p_cortex,
                p_xlinks=p_xlinks if opts.with_crosslinkers else None,
                p_myosin=p_myosin if opts.with_myosin else None,
                p_lamellipodium=p_lamellipodium if opts.with_lamellipodium else None,
                lamellipodium_geometry=lamellipodium_geometry,
                lamellipodium_polarization=lamellipodium_polarization,
                p_fa=p_fa if enable_fa else None,
                fa_clutch_capture_radius=fa_clutch_capture_radius,
                fa_clutch_k=fa_clutch_k,
                p_enclosed_volume=p_enclosed_volume,
                p_membrane_surface=p_membrane_surface,
                p_nucleus=p_nucleus,
                nucleus_centroid=nucleus_centroid,
                nucleus_gamma=nucleus_gamma,
                p_cytoplasm=p_cytoplasm,
                p_substrate=p_substrate,
                p_turnover=p_turnover,
                p_membrane=p_membrane,
                membrane_with_reaction_force=membrane_with_reaction_force,
                p_erm=p_erm,
                device=device, with_baoab=opts.with_baoab,
                constrained=constrained, constrained_dt=constrained_dt,
                reconcile_dt=reconcile_dt,
                equilibrate=equilibrate, equilibrate_steps=equilibrate_steps,
                equilibrate_softstart_steps=equilibrate_softstart_steps,
                connected_mesh=connected_mesh, cm_z_struct=cm_z_struct,
                cm_bundle_mult=cm_bundle_mult, cm_reach=cm_reach,
                rng=rng,
            )
            full_handles = handles
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
            fa_integration = handles["fa_integration"]
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
            fa=fa_integration,
            p_enclosed_volume=p_enclosed_volume,
            p_membrane_surface=p_membrane_surface,
            p_nucleus=p_nucleus,
            p_cytoplasm=p_cytoplasm,
            p_substrate=p_substrate,
            p_turnover=p_turnover,
            p_membrane=p_membrane,
            extras={"handles": full_handles} if full_handles is not None else {},
        )

    # ---------------------------------------------------------------
    # Diagnostics
    # ---------------------------------------------------------------
    def bead_count_summary(self) -> dict[str, int]:
        """Per-subsystem bead counts (sums to ``simulation.state.N_particles``).

        ``lamellipodium`` reports the CONSTRUCTION-TIME bead count.  The
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
            "fa_integrin": 0 if self.fa is None else int(self.fa.n_integrins),
            "substrate_ligand": (
                0 if self.fa is None else int(self.fa.n_substrate_ligands)
            ),
        }

    def tag_ranges(self) -> dict[str, tuple[int, int]]:
        """Map subsystem name → [start, end) tag range. Non-overlapping
        and contiguous.

        FA OFF (default): cortex → xlink → myosin → wave_particle →
        lamellipodium_actin → fa_integrin (zero-width slot). Identical to
        the pre-FA Cell (regression-clean).

        FA ON (S5 tag-space unification, 2026-05-30): the FA block is
        APPENDED LAST, so the map is cortex → xlink → myosin → wave_particle
        → lamellipodium_actin → fa_integrin (real range) → substrate_ligand.
        Because nothing is re-tagged, FA coexists with myosin / xlink /
        lamellipodium; the integrin / substrate-ligand offsets match the
        actual appended snapshot layout.

        Lamellipodium ranges reflect the CONSTRUCTION-TIME tag block;
        runtime elongation / branching events append actin_lamel beads at
        tags ``≥ lamellipodium_actin[1]``.
        """
        # Single append-order layout (S5 tag-space unification, 2026-05-30):
        # cortex → xlink → myosin → wave → lamellipodium → FA (integrin then
        # substrate_ligand). The FA block is appended LAST, so this is the
        # SAME order whether FA is on or off — when FA is off the fa_integrin
        # slot is zero-width and substrate_ligand is omitted (bit-for-bit the
        # pre-FA Cell map). When FA is on, fa_integrin is a real range followed
        # by substrate_ligand.
        offsets = {"cortex_actin": (0, self.n_cortex_actin)}
        cur = self.n_cortex_actin
        offsets["xlink_head"] = (cur, cur + self.n_xlink_heads)
        cur += self.n_xlink_heads
        if self.p_myosin is not None and self.n_myosin_particles > 0:
            myo_end = cur + self.n_myosin_particles
            offsets["myosin"] = (cur, myo_end)
            cur = myo_end
        else:
            offsets["myosin"] = (cur, cur)
        offsets["wave_particle"] = (cur, cur + self.n_wave_particles)
        cur += self.n_wave_particles
        offsets["lamellipodium_actin"] = (cur, cur + self.n_lamellipodium_actin)
        cur += self.n_lamellipodium_actin
        if self.fa is not None:
            n_int = int(self.fa.n_integrins)
            n_lig = int(self.fa.n_substrate_ligands)
            offsets["fa_integrin"] = (cur, cur + n_int)
            cur += n_int
            offsets["substrate_ligand"] = (cur, cur + n_lig)
            cur += n_lig
        else:
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
            "n_fa_integrins": (0 if self.fa is None else int(self.fa.n_integrins)),
            "n_substrate_ligands": (
                0 if self.fa is None else int(self.fa.n_substrate_ligands)
            ),
            "n_fa_clutch_bonds": (
                0 if self.fa is None else int(self.fa.n_clutch_bonds)
            ),
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
