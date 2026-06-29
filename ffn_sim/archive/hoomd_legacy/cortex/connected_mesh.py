"""Connected percolated cortex construction (CORTEX REBUILD 2026-06-04).

End-to-end builder that assembles a CONNECTED SPANNING MESH cortex — the fix for
the fragmented mesh (z = 1.3, giant 7 %, 56 % same-filament staples;
``outputs/h3/production/cortex_network.json``) that underlies the γ-floor
(Kadzik-Munro 2026: connectivity is the prerequisite for force transmission).

It composes the two rebuilt construction layers (PI directive 2026-06-04):

1. ``cortex.generate_bimodal_cortex_layout`` — BIMODAL-EXPONENTIAL filament
   length (short Arp2/3 infill + long formin = connecting backbone, Fritzsche
   2016/2017) + DISORDERED ISOTROPIC orientation (Li-Gao-Xu 2022) + 200 nm
   shell-band projection (KU-3.17).
2. ``crosslinkers.seed_connected_mesh_xlinks`` — BRIDGE-DIFFERENT-FILAMENT
   crosslinkers (no same-filament staples; Kim 2007), per-filament degree-capped
   at the literature coordination z ≈ 3-4 (Kadzik-Munro), with BUNDLING (Flormann
   2024: stiffness from bundling) → L/lc ≥ 5.9, SEEDED at construction so the
   mesh starts connected (adhered baseline).

Acceptance (the rebuild's gates, measured from the HOOMD-built seeded state):
  z (distinct-neighbour coordination) ∈ [3.0, 3.5],
  giant-component fraction ≥ 0.9,
  L/lc (crosslinks per filament) ≥ 5.9.

The dynamic ``XlinkBondUpdater`` is wired with the bridge rule + the seeded bound
state, so runtime turnover MAINTAINS connectivity (never re-staples).

Run the connectivity check / viz via ``ffn_sim.scripts.viz_cortex_network``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.cortex.cortex import (
    ResolvedH3,
    generate_bimodal_cortex_layout,
    build_variable_length_cortex_state,
)
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import (
    ResolvedCrosslinkers,
    ConnectedMeshSeed,
    seed_connected_mesh_xlinks,
    extend_cortex_state_with_xlinks,
    make_xlink_updater,
    xlink_attach_bin_names,
    xlink_attach_bin_rest_lengths,
)
from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater


def mesoscale_reach(R_cell: float, n_filaments: int) -> float:
    """Mesoscale crosslinker partner-search radius = the inter-filament spacing
    ``√(A_shell / n_fil)`` — the geometric dual of the ×40 areal coarse-graining
    (DERIVED, grid-aware; bracketed by [60 nm·√(N_native/N_eff), √(A/n_fil)]).
    """
    return math.sqrt(4.0 * math.pi * R_cell ** 2 / max(1, n_filaments))


@dataclass(slots=True)
class ConnectedCortexHandles:
    sim: hoomd.Simulation
    layout: Any                 # VariableLengthCortexLayout (bimodal)
    seed: ConnectedMeshSeed
    baoab_updater: Any
    baoab_action: Any
    xlink_updater: Any
    xlink_action: Any
    n_cortex_actin: int
    n_xlink_heads: int
    reach: float


def build_connected_cortex(
    p_cortex: ResolvedH3,
    p_xl: ResolvedCrosslinkers,
    *,
    formin_fraction: float = 0.12,
    L_long_mean: float = 5.0e-6,
    L_short_mean: float | None = None,
    z_struct: float = 2.8,
    bundle_mult: int = 3,
    arp_branch_fraction: float = 0.7,
    branch_angle_deg: float = 70.0,
    reach: float | None = None,
    n_filaments: int | None = None,
    device: hoomd.device.Device | None = None,
    with_baoab: bool = True,
    with_simulation: bool = True,
    equilibrate: bool = False,
    n_softstart: int = 200,
    n_baoab: int = 0,
    rng: np.random.Generator | None = None,
) -> ConnectedCortexHandles:
    """Build the connected percolated cortex (bimodal + seeded bridge mesh).

    Returns a :class:`ConnectedCortexHandles`.  The crosslinker count ``n_xl`` is
    an OUTPUT of the seeding (set by z_struct·bundle_mult and the geometry), not
    a config input.  ``p_xl.max_bind_dist`` is overridden to ``reach`` so the
    attach bins cover the mesoscale span and the dynamic query searches the
    mesoscale neighbourhood.
    """
    from dataclasses import replace as _replace

    if rng is None:
        rng = np.random.default_rng(p_cortex.seed)
    F = int(n_filaments) if n_filaments is not None else int(p_cortex.n_filaments)
    if reach is None:
        reach = mesoscale_reach(p_cortex.R_cell, F)
    # Mesoscale-consistent attach reach (bins + dynamic query span the mesh).
    p_xl = _replace(p_xl, max_bind_dist=reach)

    # 1. Bimodal cortex topology + Arp2/3 70° branches (faithful nucleator).
    layout = generate_bimodal_cortex_layout(
        p_cortex, formin_fraction=formin_fraction,
        L_long_mean=L_long_mean, L_short_mean=L_short_mean,
        arp_branch_fraction=arp_branch_fraction,
        branch_angle_deg=branch_angle_deg,
        n_filaments=F, project_to_shell=True, rng=rng,
    )
    n_cortex_actin = int(layout.positions_flat.shape[0])
    filament_idx = layout.filament_idx

    # 2. Seed the connected bridge mesh (bridge-different-filament + bundling),
    # counting the Arp2/3 branch bonds as prior connectivity so the TOTAL
    # coordination (branches + crosslinks) lands at z_struct.
    seed = seed_connected_mesh_xlinks(
        layout.positions_flat, filament_idx, F, p_xl,
        z_struct=z_struct, bundle_mult=bundle_mult, reach=reach,
        R_cell=p_cortex.R_cell, n_cortex_beads=n_cortex_actin,
        prior_bead_bonds=layout.branch_bonds, rng=rng,
    )

    # 3. Build HOOMD frame: bimodal cortex + seeded xlink heads + attach bonds.
    cortex_snap = build_variable_length_cortex_state(p_cortex, layout)
    snap = extend_cortex_state_with_xlinks(
        cortex_snap, seed.layout, p_xl, seeded_attach=seed.seeded_attach,
    )
    n_xlink_heads = 2 * seed.n_xl

    handles = ConnectedCortexHandles(
        sim=None, layout=layout, seed=seed,
        baoab_updater=None, baoab_action=None,
        xlink_updater=None, xlink_action=None,
        n_cortex_actin=n_cortex_actin, n_xlink_heads=n_xlink_heads, reach=reach,
    )
    if not with_simulation:
        handles.sim = snap  # caller wants the frame only
        return handles

    sim = hoomd.Simulation(
        device=device or hoomd.device.CPU(notice_level=0), seed=p_cortex.seed
    )
    sim.create_state_from_snapshot(snap)

    # ---- Forces ----
    bond = md.bond.Harmonic()
    bond.params["cortex-bond"] = dict(k=p_cortex.bond_k, r0=p_cortex.rest_length)
    if int(layout.branch_bonds.shape[0]) > 0:
        # Arp2/3 branch link: same stiffness as the backbone bond, rest length
        # = one segment (the daughter base sits one ℓ0 from the mother bead).
        bond.params["arp_branch"] = dict(k=p_cortex.bond_k, r0=p_cortex.rest_length)
    avg_intra_r0 = float(
        p_xl.alpha_fraction * p_xl.alpha_length
        + (1.0 - p_xl.alpha_fraction) * p_xl.filamin_length
    )
    bond.params["xlink_intra"] = dict(k=p_xl.k_intra, r0=avg_intra_r0)
    bin_r0 = xlink_attach_bin_rest_lengths(p_xl.n_bins, p_xl.max_bind_dist)
    for i, name in enumerate(xlink_attach_bin_names(p_xl.n_bins)):
        bond.params[name] = dict(k=p_xl.k_attach, r0=float(bin_r0[i]))

    angle = md.angle.Harmonic()
    angle.params["cortex-angle"] = dict(k=p_cortex.angle_k, t0=p_cortex.angle_t0)
    if int(layout.branch_angles.shape[0]) > 0:
        # Arp2/3 70° dendritic branch angle (Garlick ground-truth), same
        # bending stiffness as the backbone angle.
        angle.params["arp_branch_angle"] = dict(
            k=p_cortex.angle_k, t0=math.radians(branch_angle_deg))

    nlist = md.nlist.Tree(buffer=0.5 * p_cortex.lj_sigma)
    lj = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
    lj.params[("actin_cortex", "actin_cortex")] = dict(
        epsilon=p_cortex.lj_epsilon, sigma=p_cortex.lj_sigma)
    lj.r_cut[("actin_cortex", "actin_cortex")] = (
        p_cortex.lj_r_cut if p_cortex.lj_enabled else 0.0)
    # head↔head: NO LJ.  In the dense connected mesh many crosslinker heads sit
    # near the same bridge anchors (bundling); a head↔head hard core would
    # diverge at construction for coincident heads.  Crosslinker heads are
    # tiny relative to the mesh and their mutual excluded volume is negligible
    # for connectivity / mechanics — the intra + attach bonds constrain them.
    lj.params[("xlink_head", "xlink_head")] = dict(
        epsilon=p_cortex.lj_epsilon, sigma=p_cortex.lj_sigma)
    lj.r_cut[("xlink_head", "xlink_head")] = 0.0
    # head↔actin: NO LJ (heads must approach actin to bind via attach bonds).
    lj.params[("xlink_head", "actin_cortex")] = dict(
        epsilon=p_cortex.lj_epsilon, sigma=p_cortex.lj_sigma)
    lj.r_cut[("xlink_head", "actin_cortex")] = 0.0
    lj.mode = "shift"

    ig = md.Integrator(dt=p_cortex.dt_cfl)
    ig.forces.append(bond)
    ig.forces.append(angle)
    ig.forces.append(lj)
    sim.operations.integrator = ig

    baoab_updater = baoab_action = None
    if with_baoab:
        baoab_action, baoab_updater = make_baoab_updater(
            kT=p_cortex.kT,
            gamma={"actin_cortex": p_cortex.gamma_b,
                   "xlink_head": p_cortex.gamma_b},
            dt=p_cortex.dt_cfl, seed=p_cortex.seed,
        )
        sim.operations.updaters.append(baoab_updater)

    # Dynamic xlink updater: bridge-different-filament rule + seeded bound state.
    xlink_action, xlink_updater = make_xlink_updater(
        p=p_xl, layout=seed.layout, kT=p_cortex.kT,
        n_cortex_actin=n_cortex_actin,
        actin_filament_idx=filament_idx,
        head_to_actin0=seed.head_to_actin0,
    )
    sim.operations.updaters.append(xlink_updater)

    handles.sim = sim
    handles.baoab_updater = baoab_updater
    handles.baoab_action = baoab_action
    handles.xlink_updater = xlink_updater
    handles.xlink_action = xlink_action

    # Optional B2 soft-start: drain construction WCA overlaps (crossing
    # filaments within the band place some bead pairs inside the hard core) with
    # capped-displacement steepest descent BEFORE any production BAOAB, so the
    # first dt = 13 ns step doesn't overshoot.  The seeded connectivity is
    # unaffected (positions only relax locally).
    if equilibrate and with_baoab:
        from ffn_sim.archive.hoomd_legacy.cell.equilibration import equilibrate_cell
        equilibrate_cell(
            {"sim": sim, "baoab_updater": baoab_updater},
            n_softstart=n_softstart, n_baoab=n_baoab,
            rest_length=p_cortex.rest_length, gamma_b=p_cortex.gamma_b,
        )
    return handles
