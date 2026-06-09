"""Two-cell DOUBLET assembler for the explicit cadherin adherens junction (GATE-J).

H.7 platform (2026-06-09). ``cell.build`` / ``build_baseline_cell`` are SINGLE-shell:
they assemble ONE cell in one box. The explicit E-cadherin trans-dimer junction
(:mod:`ffn_sim.junction.cadherin`) is intrinsically a TWO-cell load path, so its
activation (GATE-J) needs a two-cell build: two cortex shells in one box, facing
each other across an interface, with cadherins seeded on each cell's FACING cap and
the dynamic catch-bond binder forming ``cadherin_trans`` dimers ACROSS the interface
(A↔B only — never intra-cell).

This module supplies that assembler, :func:`build_cell_doublet`, as a SEPARATE
build path (it does NOT entangle the single-cell ``Cell.build``). It builds each
cell's explicit cortex (the always-assembled actomyosin spine — backbone + angle +
α-actinin crosslinkers) via :func:`ffn_sim.cortex.cortex.build_cortex_state`,
offsets cell B along +x by ``2·R_cell + interface_gap`` so the two shells touch at
an interface, merges the two frames into one (re-indexing bonds/angles), seeds
cadherin tips on each cell's facing cap, creates ONE simulation with the shared
cortex forces + the L-M BAOAB integrator, and attaches the cadherin binder.

Scope (GATE-J build-time): a CORTEX doublet — the core actomyosin spine per cell +
the cadherin interface. The full per-cell internal compartment stack
(nucleus/MT/IF/membrane/…) on each cell of a doublet is a follow-on (the
single-cell extenders would each need to run per cell); GATE-J only needs the
two cortices + the junction. The cortical-γ no-contamination control holds on the
cortex doublet (cadherin_ denylisted → γ over cortex bonds is unchanged ON/OFF).

References
----------
- :mod:`ffn_sim.junction.cadherin` — the explicit trans-dimer catch-bond junction
  (resolver + extender + ``CadherinTransJunctionUpdater`` binder + attach helper).
- :func:`ffn_sim.cortex.cortex.build_cortex_state` — single-cell cortex frame.
- ``ffn_sim/cell/cell.py`` — the single-cell build whose snapshot-merge + BAOAB
  gamma_map idioms this mirrors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.cortex.cortex import build_cortex_state
from ffn_sim.integrator.baoab import make_baoab_updater
from ffn_sim.junction.cadherin import (
    CADHERIN_TRANS_BOND,
    attach_cadherin_junction,
    extend_snapshot_with_cadherins,
    resolve_cadherin_junction,
)


@dataclass(slots=True)
class CellDoublet:
    """A built two-cell cadherin-junction doublet.

    Attributes:
        simulation: the HOOMD Simulation holding both cells + the junction.
        p_cortex: resolved cortex params (shared by both cells).
        p_cadherin: resolved cadherin-junction params.
        cell_a_cadherin_tags / cell_b_cadherin_tags: global tags of each cell's
            cadherin particles (the interface).
        n_cortex_per_cell: cortex bead count per cell.
        cadherin_binder: the attached CadherinTransJunctionUpdater (or None).
        extras: handles (forces, layout, etc.).
    """

    simulation: Any
    p_cortex: Any
    p_cadherin: Any
    cell_a_cadherin_tags: np.ndarray
    cell_b_cadherin_tags: np.ndarray
    n_cortex_per_cell: int
    cadherin_binder: Any = None
    p_junctional_actin: Any = None
    extras: dict[str, Any] = field(default_factory=dict)


def _merge_two_cortex_frames(frame, offset_x: float):
    """Return a fresh gsd frame = cortex frame + an x-offset duplicate of it.

    Cell A is ``frame`` as-is; cell B is the SAME frame translated by
    ``+offset_x`` along x. Bonds/angles of B are re-indexed by ``n_a``. Returns
    ``(merged_frame, n_a)`` where ``n_a`` is cell A's particle count.
    """
    import gsd.hoomd

    pos = np.asarray(frame.particles.position, dtype=np.float64).reshape(-1, 3)
    n_a = int(frame.particles.N)
    types = list(frame.particles.types)
    tid = np.asarray(frame.particles.typeid, dtype=np.uint32).reshape(-1)

    pos_b = pos.copy()
    pos_b[:, 0] += float(offset_x)
    merged_pos = np.vstack([pos, pos_b])
    merged_tid = np.concatenate([tid, tid])

    out = gsd.hoomd.Frame()
    out.particles.N = 2 * n_a
    out.particles.types = types
    out.particles.typeid = merged_tid
    out.particles.position = merged_pos
    out.particles.mass = np.ones(2 * n_a, dtype=np.float64)

    # bonds: A as-is, B re-indexed by n_a.
    bg = (
        np.asarray(frame.bonds.group, dtype=np.int64).reshape(-1, 2)
        if int(frame.bonds.N) > 0 else np.empty((0, 2), dtype=np.int64)
    )
    bt = (
        np.asarray(frame.bonds.typeid, dtype=np.uint32).reshape(-1)
        if int(frame.bonds.N) > 0 else np.empty((0,), dtype=np.uint32)
    )
    btypes = list(frame.bonds.types) if frame.bonds.types else []
    if bg.shape[0] > 0:
        merged_bg = np.vstack([bg, bg + n_a]).astype(np.uint32)
        merged_bt = np.concatenate([bt, bt]).astype(np.uint32)
        out.bonds.N = int(merged_bg.shape[0])
        out.bonds.types = btypes
        out.bonds.group = merged_bg
        out.bonds.typeid = merged_bt

    # angles: A as-is, B re-indexed by n_a.
    ag = (
        np.asarray(frame.angles.group, dtype=np.int64).reshape(-1, 3)
        if int(frame.angles.N) > 0 else np.empty((0, 3), dtype=np.int64)
    )
    at = (
        np.asarray(frame.angles.typeid, dtype=np.uint32).reshape(-1)
        if int(frame.angles.N) > 0 else np.empty((0,), dtype=np.uint32)
    )
    atypes = list(frame.angles.types) if frame.angles.types else []
    if ag.shape[0] > 0:
        merged_ag = np.vstack([ag, ag + n_a]).astype(np.uint32)
        merged_at = np.concatenate([at, at]).astype(np.uint32)
        out.angles.N = int(merged_ag.shape[0])
        out.angles.types = atypes
        out.angles.group = merged_ag
        out.angles.typeid = merged_at
    elif atypes:
        out.angles.N = 0
        out.angles.types = atypes

    # box: large enough for both shells side-by-side + margin.
    L = 2.0 * (float(np.abs(merged_pos).max()) + 1.0e-6)
    out.configuration.box = [L, L, L, 0.0, 0.0, 0.0]
    return out, n_a


def _interface_cadherin_seeds(
    pos_a: np.ndarray, pos_b: np.ndarray, *,
    n_cad: int, x_iface: float, r0_trans: float,
) -> dict:
    """Seed MATCHED facing cadherin tips + their cortex anchors at the interface.

    For each of ``n_cad`` interface points (taken from cell A's +x facing-cap
    cortex beads' transverse (y, z)), place an A cadherin tip just LEFT of the
    interface midplane ``x_iface`` and a B cadherin tip just RIGHT, at the SAME
    (y, z), separated by ``r0_trans`` — so the trans-dimer is FORCE-FREE at its
    rest length AND within the binder's capture radius. Each tip is anchored to
    its cell's nearest facing-cap cortex bead (so the cadherin rides the surface,
    not free-drifting).

    Returns a dict: ``a_tips``/``b_tips`` (n_cad, 3); ``anchor_a``/``anchor_b``
    (n_cad,) cortex-bead LOCAL indices (into pos_a / pos_b).
    """
    n_cad = min(int(n_cad), pos_a.shape[0], pos_b.shape[0])
    # A facing cap = beads with the largest +x (nearest the interface).
    a_cap = np.argpartition(pos_a[:, 0], -n_cad)[-n_cad:]
    yz = pos_a[a_cap, 1:]                       # transverse interface points
    # B facing cap = beads with the smallest +x (B's −x face toward A).
    b_cap_pool = np.argpartition(pos_b[:, 0], n_cad * 4)[: n_cad * 4]
    # match each A interface point to the nearest B facing-cap bead in (y, z).
    from scipy.spatial import cKDTree
    tree = cKDTree(pos_b[b_cap_pool, 1:])
    _d, j = tree.query(yz, k=1)
    anchor_b = b_cap_pool[np.atleast_1d(j)]
    anchor_a = a_cap

    # A-tip and B-tip share the SAME transverse (y, z) and sit ±r0_trans/2 across
    # the interface midplane → their separation is EXACTLY r0_trans (the trans-
    # dimer rest length), so a seeded cadherin_trans bond is force-free AND within
    # the binder's capture radius. (The B anchor cortex bead carries the match
    # residual in (y,z) on the soft anchor bond, not the trans-dimer.)
    half = 0.5 * float(r0_trans)
    a_tips = np.column_stack([np.full(n_cad, x_iface - half), yz])
    b_tips = np.column_stack([np.full(n_cad, x_iface + half), yz])
    return {
        "a_tips": a_tips, "b_tips": b_tips,
        "anchor_a": anchor_a, "anchor_b": anchor_b,
    }


def build_cell_doublet(
    path_or_name: str = "mcf7_baseline.yaml",
    *,
    manifest: dict | None = None,
    interface_gap: float | None = None,
    n_cad_per_cell: int | None = None,
    device: Any | None = None,
    seed: int = 1,
    with_baoab: bool = True,
    with_junctional_actin: bool = False,
    run_binder_batches: int = 0,
) -> CellDoublet:
    """Assemble a two-cortex-cell cadherin-junction doublet (GATE-J build).

    Args:
        path_or_name: manifest path/name (for the cortex params).
        manifest: optional pre-loaded manifest (overrides path_or_name).
        interface_gap: surface-to-surface gap [m] between the two shells
            (default = 0.5·r_bind so facing cadherin tips fall within capture).
        n_cad_per_cell: cadherins seeded per cell (default from the resolved
            junction's ``n_cad_per_cell``).
        device: HOOMD device (default CPU).
        seed: RNG / sim seed.
        with_baoab: attach the L-M BAOAB integrator (else a no-op integrator
            is still set so the binder can read dt).
        run_binder_batches: if > 0, run the binder this many batches so
            trans-dimers form (the gate uses this to count A↔B dimers).

    Returns:
        A :class:`CellDoublet`.
    """
    from ffn_sim.cell.manifest import load_manifest, resolve_baseline

    mani = manifest if manifest is not None else load_manifest(path_or_name)
    rb = resolve_baseline(mani, allow_no_nucleus=True)
    p_cortex = rb.p_cortex

    # ---- cadherin junction params (needs contact-zone + a cadherin count) ----
    cad_block = mani.get("optional_subsystems", {}).get("cadherin_junction", {})
    cad_cfg = dict(cad_block) if isinstance(cad_block, dict) else {}
    cad_cfg["enabled"] = True
    if n_cad_per_cell is not None:
        cad_cfg["n_cad_per_cell"] = int(n_cad_per_cell)
    # Engagement length: a fraction of the cortex bead spacing (the interface
    # contact zone over which a trans-dimer's force reaches f0). Geometry-derived.
    contact_zone = float(cad_cfg.get("contact_zone_width", p_cortex.rest_length))
    # Capture radius >= rest length so the binder re-captures a seeded dimer (and
    # the matched A↔B tips at separation r0_trans fall within it). Default 2× zone.
    cad_cfg.setdefault("r_bind", 2.0 * contact_zone)
    dtc = float(rb.dtc)
    p_cad = resolve_cadherin_junction(
        {"cadherin": cad_cfg}, kT=p_cortex.kT, dt=dtc,
        contact_zone_width=contact_zone,
    )

    # ---- build one cortex shell, duplicate + offset → doublet frame ----
    frame_a, _topo, _xl = build_cortex_state(
        p_cortex, with_crosslinkers=True, rng=np.random.default_rng(seed),
    )
    R_cell = float(p_cortex.R_cell)
    gap = interface_gap if interface_gap is not None else 0.5 * p_cad.r_bind
    offset_x = 2.0 * R_cell + float(gap)
    snap, n_a = _merge_two_cortex_frames(frame_a, offset_x)

    # ---- seed MATCHED facing cadherins at the interface (+ cortex anchors) ----
    pos = np.asarray(snap.particles.position, dtype=np.float64).reshape(-1, 3)
    pos_a = pos[:n_a]
    pos_b = pos[n_a:]
    n_cad = int(p_cad.n_cad_per_cell)
    x_iface = R_cell + 0.5 * float(gap)
    seeds = _interface_cadherin_seeds(
        pos_a, pos_b, n_cad=n_cad, x_iface=x_iface, r0_trans=p_cad.r0_trans,
    )
    anchor_a_local = np.asarray(seeds["anchor_a"], dtype=np.int64)         # into A
    anchor_b_global = np.asarray(seeds["anchor_b"], dtype=np.int64) + n_a  # into merged
    cad_info = extend_snapshot_with_cadherins(
        snap, p_cad,
        cell_a_surface_points=seeds["a_tips"], cell_b_surface_points=seeds["b_tips"],
    )
    a_cad_tags = np.asarray(cad_info["cell_a_tags"], dtype=np.int64)
    b_cad_tags = np.asarray(cad_info["cell_b_tags"], dtype=np.int64)

    # ---- anchor each cadherin to its cell's nearest facing-cap cortex bead ----
    # cadherin_anchor (cadherin_ prefix → γ-denylisted). r0 = construction sep →
    # force-free; keeps the cadherin riding the surface instead of free-drifting.
    cad_pos = np.asarray(snap.particles.position, dtype=np.float64).reshape(-1, 3)
    anchor_pairs = np.concatenate([
        np.column_stack([a_cad_tags, anchor_a_local]),
        np.column_stack([b_cad_tags, anchor_b_global]),
    ], axis=0).astype(np.int64)
    anchor_r0 = np.linalg.norm(
        cad_pos[anchor_pairs[:, 0]] - cad_pos[anchor_pairs[:, 1]], axis=1
    )
    CAD_ANCHOR_BOND = "cadherin_anchor"
    bond_types = list(snap.bonds.types)
    for nm in (CADHERIN_TRANS_BOND, CAD_ANCHOR_BOND):
        if nm not in bond_types:
            bond_types.append(nm)
    anchor_tid = bond_types.index(CAD_ANCHOR_BOND)
    trans_tid = bond_types.index(CADHERIN_TRANS_BOND)
    # ---- SEED the trans-dimers pre-bound (engaged-junction baseline) ----
    # Physiological baseline = an ENGAGED adherens junction. Matched A↔B tips sit
    # r0_trans apart → each seeded cadherin_trans bond is force-free at build. The
    # binder then MAINTAINS them (catch-slip break/rebind), exactly the cortex
    # seeded-adhered pattern (per-batch binding from scratch is k_on·Δt≪1, so a
    # bare junction would take ~1e3 batches to engage — invalid as a baseline).
    trans_pairs = np.column_stack([a_cad_tags, b_cad_tags]).astype(np.int64)
    old_bg = (
        np.asarray(snap.bonds.group, dtype=np.int64).reshape(-1, 2)
        if int(snap.bonds.N) > 0 else np.empty((0, 2), dtype=np.int64)
    )
    old_bt = (
        np.asarray(snap.bonds.typeid, dtype=np.uint32).reshape(-1)
        if int(snap.bonds.N) > 0 else np.empty((0,), dtype=np.uint32)
    )
    merged_bg = np.vstack([old_bg, anchor_pairs, trans_pairs]).astype(np.uint32)
    merged_bt = np.concatenate([
        old_bt,
        np.full(anchor_pairs.shape[0], anchor_tid, dtype=np.uint32),
        np.full(trans_pairs.shape[0], trans_tid, dtype=np.uint32),
    ])
    snap.bonds.types = bond_types
    snap.bonds.N = int(merged_bg.shape[0])
    snap.bonds.group = merged_bg
    snap.bonds.typeid = merged_bt
    # per-cadherin anchor rest lengths (force-free) for the bond registration.
    _anchor_r0_mean = float(np.mean(anchor_r0)) if anchor_r0.size else 0.0
    _n_trans_seeded = int(trans_pairs.shape[0])

    # ---- (optional) junctional-actin belt: couple each cell's cadherins to its
    #      OWN cortex via the α-catenin/vinculin clutch (BEFORE create_state). ----
    p_ja = None
    if with_junctional_actin:
        from ffn_sim.junction.junctional_actin import (
            extend_snapshot_with_junctional_actin,
            register_junctional_actin_bond_params,
            resolve_junctional_actin,
        )
        ja_block = mani.get("optional_subsystems", {}).get("junctional_actin", {})
        ja_cfg = dict(ja_block) if isinstance(ja_block, dict) else {}
        ja_cfg["enabled"] = True
        p_ja = resolve_junctional_actin(
            {"junctional_actin": ja_cfg}, kT=p_cortex.kT, dt=dtc,
            n_cadherin=int(p_cad.n_cad_per_cell),
        )
        if not p_ja.is_anchored:
            raise RuntimeError(
                "junctional_actin enabled but its catch-set constants are not "
                "anchored — supply the PI-candidate set in the manifest "
                "optional_subsystems.junctional_actin block (k_couple, k_anchor, "
                "k_catch0, x_catch, k_slip0, x_slip, k_on, max_couple_dist, "
                "anchor_r0). See junctional_actin.PI_DECISIONS."
            )
        # cell A cadherins → cell A cortex [0, n_a); cell B → [n_a, 2n_a).
        snap = extend_snapshot_with_junctional_actin(
            snap, p_ja, cadherin_tags=a_cad_tags,
            cortex_positions=pos_a, cortex_tag_start=0,
        )
        snap = extend_snapshot_with_junctional_actin(
            snap, p_ja, cadherin_tags=b_cad_tags,
            cortex_positions=pos_b, cortex_tag_start=n_a,
        )

    # ---- create the simulation + cortex forces ----
    sim = hoomd.Simulation(device=device or hoomd.device.CPU(), seed=seed)
    sim.create_state_from_snapshot(snap)

    bond = md.bond.Harmonic()
    for bt in sim.state.bond_types:
        if bt == CADHERIN_TRANS_BOND:
            # cadherin_trans is carried by attach_cadherin_junction's OWN Harmonic
            # (below); zero it here so the two Harmonics do NOT double-apply it.
            bond.params[bt] = dict(k=0.0, r0=0.0)
        elif bt == "cadherin_anchor":
            # cadherin → cortex anchor (holds the cadherin on the surface); soft
            # spring at the mean construction separation (force-negligible).
            bond.params[bt] = dict(k=p_cortex.bond_k, r0=_anchor_r0_mean)
        elif bt == "cortex-bond":
            bond.params[bt] = dict(k=p_cortex.bond_k, r0=p_cortex.rest_length)
        else:
            # crosslinker / other cortex bond types: cortex backbone params as a
            # safe default (α-actinin xlink bonds; GATE-J observable is the
            # junction, not xl tension).
            bond.params[bt] = dict(k=p_cortex.bond_k, r0=p_cortex.rest_length)
    # H.junctional_actin: override the junc_actin_anchor + junc_actin_couple_b{i}
    # params with their real (anchored-candidate) stiffness/rest-length.
    if p_ja is not None:
        register_junctional_actin_bond_params(bond, p_ja)
    angle = md.angle.Harmonic()
    for at in sim.state.angle_types:
        angle.params[at] = dict(k=p_cortex.angle_k, t0=p_cortex.angle_t0)

    integrator = md.Integrator(dt=dtc)
    integrator.forces.append(bond)
    if len(sim.state.angle_types) > 0:
        integrator.forces.append(angle)
    sim.operations.integrator = integrator

    if with_baoab:
        gamma_map = {"actin_cortex": p_cortex.gamma_b}
        if "cadherin" in sim.state.particle_types:
            # cadherin tips diffuse with the surface; cortex-bead drag order.
            gamma_map["cadherin"] = p_cortex.gamma_b
        if "junc_actin" in sim.state.particle_types:
            gamma_map["junc_actin"] = p_cortex.gamma_b
        _baoab_action, _baoab_updater = make_baoab_updater(
            kT=p_cortex.kT, gamma=gamma_map, dt=dtc, seed=seed + 7,
        )
        sim.operations.updaters.append(_baoab_updater)

    # ---- attach the cadherin binder (forms A↔B trans-dimers) ----
    a_cad = np.asarray(cad_info["cell_a_tags"], dtype=np.int64)
    b_cad = np.asarray(cad_info["cell_b_tags"], dtype=np.int64)
    harmonic, updater = attach_cadherin_junction(
        sim, p_cad, cell_a_tags=a_cad, cell_b_tags=b_cad, seed=seed + 11,
        gamma_cad=p_cortex.gamma_b, cfl_strict=False,
    )
    binder = getattr(sim, "_cadherin_binder", None)

    if run_binder_batches > 0 and binder is not None:
        sim.run(int(run_binder_batches) * int(p_cad.batch_steps))

    return CellDoublet(
        simulation=sim,
        p_cortex=p_cortex,
        p_cadherin=p_cad,
        cell_a_cadherin_tags=a_cad,
        cell_b_cadherin_tags=b_cad,
        n_cortex_per_cell=n_a,
        cadherin_binder=binder,
        p_junctional_actin=p_ja,
        extras={"handles": {"bond_force": harmonic, "cadherin_updater": updater,
                            "cadherin_info": cad_info,
                            "n_trans_seeded": _n_trans_seeded,
                            "junc_actin": getattr(snap, "_junc_actin", None),
                            "n_junc_actin_heads": int(
                                (np.asarray(sim.state.get_snapshot().particles.typeid)
                                 == list(sim.state.particle_types).index("junc_actin")).sum()
                            ) if "junc_actin" in sim.state.particle_types else 0}},
    )
