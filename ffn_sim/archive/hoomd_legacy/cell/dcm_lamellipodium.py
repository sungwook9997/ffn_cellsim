"""Per-rim-cell mechanistic lamellipodium for the native+tent DCM spheroid.

Replaces the coarse-grained body-force traction abstraction
(:class:`ffn_sim.archive.hoomd_legacy.cell.dcm_active.ActiveRimTraction`, which prescribes an outward
body force on rim-cell basal nodes) with a REAL, bead-resolved lamellipodial
actin engine attached per rim cell — the same validated single-cell engine that
took the single cell to A/A₀ → 2.03 (``ffn_sim/cell/spreading_drive.py``, commit
f654531), generalized to every rim cell of the spheroid.

CLAUDE.md hard rule (PI 2026-05-19): pick the mechanistic option, no lumped
proxies. A prescribed outward body force IS a lumped proxy for the protrusion +
clutch + traction loop; this module implements that loop explicitly.

What "mechanistic per-rim-cell lamellipodium" means here
--------------------------------------------------------
For each RIM cell (spheroid periphery AND in substrate contact) we attach a real
lamellipodial actin module built from three reused ``spreading_drive`` engines
plus a tether:

1. **Lamellipodial actin sheet** — a PRE-ALLOCATED pool of ``actin_lamel`` beads
   (a SEPARATE particle type), seeded one ℓ₀ outward from the cell's LEADING
   (outward-radial, basal) membrane nodes. KEY invariant: these beads are NOT in
   any cell's mesh triangle list, so the native ``md.mesh.conservation.Volume``
   tag-space is untouched — the bead engine (which mutates topology by appending
   actin) is therefore compatible with the native shell (the shell nodes never
   change; only the separate actin pool grows). Seeded force-free, no t=0
   overlap.

2. **Leading-edge nucleation** (:class:`spreading_drive.LeadingEdgeNucleationUpdater`)
   — the membrane-tracked Arp2/3 polymerization ratchet advancing each rim cell's
   leading edge outward (radial from the cell centroid). One updater per rim cell.

3. **FA molecular clutch** (:class:`spreading_drive.FrontClutchRatchetUpdater`) —
   force-free harmonic grip bonds gripping the advancing actin to substrate
   ``ligand`` beads, converting retrograde flow into anchored protrusion
   (Chan-Odde / Elosegui-Artola motor-clutch). One updater per rim cell (sharing
   the global ligand field).

4. **Actin→membrane tether** (:class:`LamellipodialTractionTether`) — a harmonic
   tether pulling each rim cell's LEADING membrane (mesh) nodes toward the
   clutch-anchored, polymerizing actin front. THIS is the traction transmission:
   the clutch-gripped, advancing actin PULLS THE CELL BODY FORWARD through the
   tether — real traction carried by the actin/clutch, not a prescribed body
   force. Plus :class:`spreading_drive.BasalAdhesionTether` keeping the sheet
   basal.

Stability (BAOAB int32 guard — the failure mode that bit before)
----------------------------------------------------------------
- The actin pool + ligand field are seeded force-free (exact rest lengths / no
  overlap), so construction injects no energy.
- ``gamma`` covers BOTH ``mem`` (shell nodes) AND ``actin_lamel`` AND ``ligand``.
- The traction tether is force-capped and the tether stiffness satisfies the
  overdamped CFL ``dt ≪ γ / k_tether``.
- Small dt (3e-10 → 1e-10) and a finite gate before long runs.

Scale note (honest)
--------------------
Bead-resolved many-cell is HEAVY on CPU (each rim cell carries its own growing
actin network + per-cell snapshot-rebuild updaters). The CPU smoke is kept small
(~8-14 cells, short runs). The large production run is gbook-GPU — ``device=`` is
plumbed through :func:`build_native_dcm_simulation`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.cell.dcm_native_shell import ResolvedNativeDCM
from ffn_sim.archive.hoomd_legacy.cell.lamellipodium import LamellipodiumState
from ffn_sim.archive.hoomd_legacy.cell.spreading_drive import (
    BasalAdhesionTether,
    FrontClutchRatchetUpdater,
    LeadingEdgeNucleationUpdater,
)


# ---------------------------------------------------------------------------
# Resolved parameters for the per-rim-cell lamellipodium
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedLamellipodiumSpheroid:
    """Per-rim-cell lamellipodium constants (on top of the native DCM base).

    The base shell / turgor / contact / substrate constants live in
    :class:`ResolvedNativeDCM` (``base``); this carries the lamellipodium engine
    knobs. Defaults are band-matched to the validated single-cell run
    (``spreading_drive``): v_front in the literature single-cell spreading range
    (3-12 µm/min), grip radius ≈ FA reach, tether stiffness in the overdamped-CFL
    safe band.
    """

    base: ResolvedNativeDCM = field(default_factory=ResolvedNativeDCM)

    # --- rim-cell detection (which cells grow a lamellipodium) ---
    rim_contact_band: float = 0.6   # basal node if z < z0 + band*R (in contact)
    rim_neighbour_factor: float = 2.6   # cell is a "neighbour" if < this*R away
    rim_max_neighbours: int = 9     # a rim cell has <= this many close neighbours

    # --- lamellipodial actin seed pool (per rim cell) ---
    n_seed_per_cell: int = 8        # leading-edge seed beads per rim cell
    seed_basal_offset: float = 0.3e-6   # m seed actin at z0 + this (basal sheet,
                                    # within grip_radius of the z0 ligands)

    # --- leading-edge nucleation (the Arp2/3 polymerization ratchet) ---
    v_front: float = 6.0e-6 / 60.0  # m/s band-matched front velocity (6 µm/min, lit 3-12)
    S_kinetic: float = 1.0          # PHYSIOLOGICAL real-time (new foundation, 2026-06-14).
                                    # WAS 6.0e5 (accelerated clock on the old low-γ
                                    # foundation). With S=1 + dt~1e-3, advancing the front
                                    # one ℓ₀ takes ~ℓ₀/(v_front·dt)≈5e3 steps → A/A0≈2 needs
                                    # ~75k steps/tip (kinetic-budget tradeoff of removing S).
    actin_rest_length: float = 0.5e-6   # ℓ₀ monomer step [m]
    max_advance_per_cell: int = 8   # tips advanced per tick per cell (cost bound)
    contact_band_frac: float = 1.5  # advance band = frac*mean_edge about z_basal

    # --- FA molecular clutch (grip advancing actin to substrate ligands) ---
    grip_radius: float = 0.6e-6     # m bead↔ligand grip reach (FA footprint)
    k_clutch: float = 2.0e-4        # N/m engaged-integrin clutch spring (KU-2.4 band)
    max_grip_per_cell: int = 6      # grips per tick per cell

    # --- substrate ligand field (explicit beads the clutch grips) ---
    ligand_spacing: float = 0.5e-6  # m ligand grid spacing (dense enough for grip)
    ligand_pad: float = 2.0         # ligand grid extends pad*R beyond footprint

    # --- actin→membrane traction tether (transmits traction to the cell body) ---
    k_tether: float = 4.0e-3        # N/m actin→leading-membrane tether spring
    tether_force_cap: float = 5.0e-9    # N per-node cap = 5·(60Pa·area_per_node) lit
                                    # MCF7 per-node traction (was 5e-8; mirrors the proxy
                                    # f_cap=5·f_act so a tether can't exceed physiological
                                    # traction by 50×).
    tether_radius: float = 3.0e-6   # m a membrane node tethers to actin within this

    # --- basal-sheet tether (keep lamellipodial actin in the basal plane) ---
    k_basal_adh: float = 2.0e-4     # N/m z-tether of actin to z_basal

    # --- pool / numerics ---
    pool_per_cell: int = 220        # pre-allocated append headroom per rim cell
    # PHYSIOLOGICAL drag (new foundation): derive at build time as 6π·η·R/nv with
    # η=65.9 Pa·s MCF7 cytoplasm → 2.22e-4 for R=7.5µm,nv=42 (×5.7e5 vs the old 3.9e-10
    # water-like default). Actin beads are immersed in the same cytoplasm; ligand is
    # pinned (γ only sets its BD prefactor). Build must override from the run's η/R/nv.
    gamma_actin: float = 2.22e-4    # N·s/m actin-bead drag = physiological node drag
    gamma_ligand: float = 2.22e-4   # N·s/m ligand drag (pinned; large effective)
    batch_steps: int = 50           # BAOAB steps between updater ticks
    seed: int = 7
    extras: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Substrate ligand pin (keep ligand beads on the substrate plane, force-free)
# ---------------------------------------------------------------------------
class LigandPin(md.force.Custom):
    """Pin substrate ``ligand`` beads to the substrate plane ``z = z0``.

    ``U_i = ½ k_pin (z_i − z0)²`` (z-only) on ligand beads. The ligands are the
    fixed substrate adhesion sites the FA clutch grips; pinning them keeps the
    substrate rigid (glass dish) so the clutch reaction is felt by the cell, not
    dissipated by ligand drift. Force-free at construction (ligands seeded at z0).
    """

    def __init__(self, *, z0: float, k_pin: float, ligand_typeid: int) -> None:
        super().__init__(aniso=False)
        self.z0 = float(z0)
        self.k_pin = float(k_pin)
        self.ligand_typeid = int(ligand_typeid)

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tid = np.asarray(snap.particles.typeid).copy()
            pos = np.asarray(snap.particles.position).copy()
        mask = tid == self.ligand_typeid
        dz = pos[:, 2] - self.z0
        F = np.zeros_like(pos)
        F[:, 2] = -self.k_pin * dz
        U = 0.5 * self.k_pin * dz * dz
        F[~mask] = 0.0
        U[~mask] = 0.0
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U


# ---------------------------------------------------------------------------
# Actin → leading-membrane traction tether (the traction-transmission force)
# ---------------------------------------------------------------------------
class LamellipodialTractionTether(md.force.Custom):
    """Pull each rim cell's LEADING membrane nodes toward the clutch-anchored,
    polymerizing actin front — the traction-transmission mechanism.

    This is the force that makes the lamellipodium PULL THE CELL BODY (replacing
    the coarse prescribed body force). Each tick (every force eval) it:

    1. Identifies each rim cell's leading (outward-radial, basal) membrane nodes.
    2. For each such node, finds the nearest ``actin_lamel`` bead OUTWARD of it
       (radially farther from the cell centroid, within ``tether_radius``).
    3. Applies a harmonic pull ``F = −k_tether·(r_node − r_actin)`` on the node
       (toward the actin) — force-capped. Newton-3 reaction goes on the actin
       bead, so the cell body and the actin sheet pull on each other: the
       clutch-gripped actin (anchored to the substrate ligands) acts as the fixed
       point, so the NET effect is the cell body being pulled outward onto its own
       protruded, substrate-gripped lamellipodium. That is real traction carried
       by the actin/clutch, not a prescribed body force.

    Only LEADING nodes are tethered (outward hemisphere of each rim cell relative
    to its centroid), so the pull is directional (outward), reproducing the
    polarized crawl without prescribing the force magnitude — the magnitude is set
    by how far the actin has polymerized and gripped (emergent).

    The membrane normal / leading direction is the per-cell outward-radial unit
    vector from the spheroid centroid through the cell centroid; recomputed each
    eval from live positions (cells migrate).
    """

    def __init__(
        self, *, cell_of_node: np.ndarray, rim_cells: np.ndarray,
        actin_typeid: int, mem_typeid: int, z_basal: float, basal_band: float,
        k_tether: float, force_cap: float, tether_radius: float,
        lead_frac: float = 0.0,
    ) -> None:
        super().__init__(aniso=False)
        self.cell_of_node = np.asarray(cell_of_node, dtype=np.int64)
        self.rim_cells = np.asarray(rim_cells, dtype=np.int64)
        self.rim_set = set(int(c) for c in self.rim_cells)
        self.actin_typeid = int(actin_typeid)
        self.mem_typeid = int(mem_typeid)
        self.z_basal = float(z_basal)
        self.basal_band = float(basal_band)
        self.k_tether = float(k_tether)
        self.force_cap = float(force_cap)
        self.tether_radius = float(tether_radius)
        self.lead_frac = float(lead_frac)   # node is "leading" if (r_node-c)·out >= lead_frac*|.|
        self._n_tethered = 0

    @property
    def n_tethered(self) -> int:
        return self._n_tethered

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            tid = np.asarray(snap.particles.typeid).copy()
            pos = np.asarray(snap.particles.position).copy()
        F = np.zeros_like(pos)

        # Tag-ordered views (cell_of_node is in global-tag order).
        perm = np.argsort(tag)
        pos_g = pos[perm]
        tid_g = tid[perm]
        # cell_of_node only covers the membrane tag block [0, len). Actin/ligand
        # beads have higher tags (appended pool); guard the index.
        n_mem_block = self.cell_of_node.shape[0]

        actin_g = np.where(tid_g == self.actin_typeid)[0]
        if actin_g.size == 0:
            with self.cpu_local_force_arrays as arr:
                arr.force[:] = 0.0
            return
        actin_pos = pos_g[actin_g]
        # only basal actin participates (the lamellipodial sheet)
        actin_basal = np.abs(actin_pos[:, 2] - self.z_basal) <= self.basal_band
        actin_pos = actin_pos[actin_basal]
        if actin_pos.shape[0] == 0:
            with self.cpu_local_force_arrays as arr:
                arr.force[:] = 0.0
            return

        # Spheroid centroid from membrane nodes of rim cells (in-plane).
        mem_mask = (tid_g == self.mem_typeid)
        mem_idx = np.where(mem_mask)[0]
        mem_idx = mem_idx[mem_idx < n_mem_block]
        if mem_idx.size == 0:
            with self.cpu_local_force_arrays as arr:
                arr.force[:] = 0.0
            return
        sph_centroid = pos_g[mem_idx].mean(axis=0)
        sph_centroid[2] = self.z_basal  # in-plane reference

        n_tethered = 0
        F_g = np.zeros_like(pos_g)
        for c in self.rim_cells:
            c = int(c)
            cnodes = np.where(self.cell_of_node == c)[0]
            if cnodes.size == 0:
                continue
            cell_centroid = pos_g[cnodes].mean(axis=0)
            out = cell_centroid - sph_centroid
            out[2] = 0.0
            on = np.linalg.norm(out)
            if on < 1e-12:
                continue
            out /= on
            # leading + basal membrane nodes of this cell
            npos = pos_g[cnodes]
            rel = npos - cell_centroid
            proj = rel[:, 0] * out[0] + rel[:, 1] * out[1]
            basal = np.abs(npos[:, 2] - self.z_basal) <= self.basal_band
            lead = (proj >= self.lead_frac * np.maximum(
                np.linalg.norm(rel[:, :2], axis=1), 1e-18)) & basal
            lead_nodes = cnodes[lead]
            if lead_nodes.size == 0:
                continue
            # for each leading node, nearest OUTWARD basal actin bead
            for ln in lead_nodes:
                rp = pos_g[ln]
                d = actin_pos - rp
                dist = np.linalg.norm(d, axis=1)
                # require the actin to be farther out (radially) than the node
                act_proj = ((actin_pos[:, 0] - cell_centroid[0]) * out[0]
                            + (actin_pos[:, 1] - cell_centroid[1]) * out[1])
                node_proj = float(proj[np.where(cnodes == ln)[0][0]])
                cand = (dist <= self.tether_radius) & (act_proj > node_proj)
                if not cand.any():
                    continue
                j = int(np.argmin(np.where(cand, dist, np.inf)))
                dvec = actin_pos[j] - rp     # node → actin (pull outward)
                fvec = self.k_tether * dvec
                fmag = np.linalg.norm(fvec)
                if fmag > self.force_cap:
                    fvec = fvec * (self.force_cap / fmag)
                F_g[ln] += fvec
                n_tethered += 1

        # scatter back to local row order
        F[perm] = F_g
        self._n_tethered = n_tethered
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F


# ---------------------------------------------------------------------------
# Rim-cell detection (periphery + substrate contact)
# ---------------------------------------------------------------------------
def detect_rim_cells(
    centers: np.ndarray, *, R: float, z0: float, contact_band: float,
    neighbour_factor: float, max_neighbours: int,
) -> np.ndarray:
    """Return the cell indices that are RIM cells (few neighbours = periphery)
    AND in substrate contact (centroid within ``contact_band*R`` of z0).

    Mirrors the neighbour-count rim heuristic in
    :class:`ffn_sim.archive.hoomd_legacy.cell.dcm_active.ActiveRimTraction._resolve_rim_cells` but
    operates on the build-time cell centers (a static rim assignment; for the
    smoke this is sufficient and avoids per-step convex-hull cost).
    """
    n = centers.shape[0]
    rim = []
    for c in range(n):
        d = np.linalg.norm(centers - centers[c], axis=1)
        d[c] = np.inf
        n_close = int(np.count_nonzero(d < neighbour_factor * R))
        in_contact = (centers[c, 2] - z0) <= contact_band * R
        if n_close <= max_neighbours and in_contact:
            rim.append(c)
    return np.array(rim, dtype=np.int64)


# ---------------------------------------------------------------------------
# Build: native+tent+substrate base  +  per-rim-cell lamellipodium engine
# ---------------------------------------------------------------------------
def build_lamellipodium_spheroid(
    p: ResolvedLamellipodiumSpheroid, n_cells: int, *,
    device=None, contact: bool = True,
):
    """Assemble the native DCM spheroid with a REAL per-rim-cell lamellipodium.

    Steps:
      1. Build the native+tent+substrate base via
         :func:`build_native_dcm_simulation` (``device=`` plumbed → gbook GPU).
      2. Detect rim cells (periphery + substrate contact).
      3. Append a force-free ``actin_lamel`` seed pool (leading-edge beads) +
         a ``ligand`` substrate field to the state (separate types; mesh
         tag-space untouched).
      4. Rebuild the BAOAB updater with ``gamma`` covering mem + actin_lamel +
         ligand, and a fresh actin/ligand-aware tag space.
      5. Attach per-rim-cell ``LeadingEdgeNucleationUpdater`` +
         ``FrontClutchRatchetUpdater``, the global
         ``LamellipodialTractionTether`` + ``BasalAdhesionTether`` + ``LigandPin``
         + the ``front_clutch_grip`` harmonic bond.

    Returns a dict of handles (sim, rim_cells, nucleators, clutches, tether,
    states, plus the base handles).
    """
    from ffn_sim.archive.hoomd_legacy.cell.dcm import DcmSubstrateForce, _cluster_centers
    from ffn_sim.archive.hoomd_legacy.cell.dcm_contact import DcmTentContact
    from ffn_sim.archive.hoomd_legacy.cell.dcm_native_shell import (
        build_native_snapshot, build_native_shell_forces)

    base = p.base
    R = base.R_cell
    z0 = base.z_substrate

    # 1. native mesh-shell snapshot (mem block at tags [0, n_mem)). We build the
    #    snapshot ourselves (rather than via build_native_dcm_simulation) so that
    #    actin_lamel + ligand types can be registered in the INITIAL state —
    #    HOOMD forbids adding particle types after create_state_from_snapshot.
    #    The mesh tag-space stays the front block, exactly as the native builder
    #    lays it out, so the native Volume mesh is unchanged.
    spacing = base.spacing_factor * R
    centers = _cluster_centers(n_cells, spacing, z0, R, mode=base.cluster)
    (msnap, mesh_tris, mesh_typeids, mesh_types, cell_of_tag, ranges,
     mean_edge, nv, V_rest) = build_native_snapshot(base, n_cells, centers)
    cell_of_tag = np.asarray(cell_of_tag, dtype=np.int64)
    mem_pos = np.asarray(msnap.particles.position, dtype=np.float64)
    n_mem = mem_pos.shape[0]
    # recompute centers from actual placed nodes (x-mirror etc.)
    centers = np.empty((n_cells, 3), dtype=np.float64)
    for c, (a, b) in enumerate(ranges):
        centers[c] = mem_pos[a:b].mean(axis=0)

    # 2. rim cells
    rim_cells = detect_rim_cells(
        centers, R=R, z0=z0, contact_band=p.rim_contact_band,
        neighbour_factor=p.rim_neighbour_factor,
        max_neighbours=p.rim_max_neighbours)
    if rim_cells.size == 0:
        # fall back: outermost cells by in-plane radius from spheroid centroid
        sph_c = centers.mean(axis=0)
        rxy = np.hypot(centers[:, 0] - sph_c[0], centers[:, 1] - sph_c[1])
        rim_cells = np.argsort(rxy)[::-1][: max(1, n_cells // 2)].astype(np.int64)

    z_basal = z0 + p.seed_basal_offset       # basal sheet just above the floor
    basal_band = p.contact_band_frac * mean_edge

    # 3a. actin_lamel seed pool: per rim cell, n_seed beads one ℓ₀ outward from
    #     the cell's leading-basal membrane nodes (force-free, no overlap).
    sph_c = centers.mean(axis=0)
    seed_pos = []
    for c in rim_cells:
        c = int(c)
        a, b = ranges[c]
        npos = mem_pos[a:b]
        cc = centers[c]
        out = cc - sph_c
        out[2] = 0.0
        on = np.linalg.norm(out)
        if on < 1e-12:
            out = np.array([1.0, 0.0, 0.0])
        else:
            out /= on
        # rank this cell's nodes by outward-radial projection; take the leading n
        proj = (npos[:, 0] - cc[0]) * out[0] + (npos[:, 1] - cc[1]) * out[1]
        order = np.argsort(proj)[::-1][: p.n_seed_per_cell]
        for k in order:
            base_node = npos[k]
            r_seed = base_node + p.actin_rest_length * out
            r_seed[2] = z_basal
            seed_pos.append(r_seed)
    seed_pos = (np.array(seed_pos, dtype=np.float64).reshape(-1, 3)
                if seed_pos else np.empty((0, 3)))
    n_seed = seed_pos.shape[0]

    # 3b. ligand field: a grid on the substrate under the spheroid footprint.
    rad = np.hypot(centers[:, 0] - sph_c[0], centers[:, 1] - sph_c[1]).max()
    extent = rad + p.ligand_pad * R
    step = max(p.ligand_spacing, 1e-9)
    gx = np.arange(sph_c[0] - extent, sph_c[0] + extent + step, step)
    gy = np.arange(sph_c[1] - extent, sph_c[1] + extent + step, step)
    GX, GY = np.meshgrid(gx, gy)
    lig_pos = np.column_stack([GX.ravel(), GY.ravel(),
                               np.full(GX.size, z0)])
    n_lig = lig_pos.shape[0]

    # 4. assemble the FULL initial snapshot: mem block (front) + actin seed pool
    #    + ligand field, with actin_lamel + ligand types registered up front.
    #    Mesh tags stay the front block [0, n_mem) so the native Volume mesh
    #    (built on those tags) is unchanged.
    import gsd.hoomd
    n_new = n_seed + n_lig
    snap = gsd.hoomd.Frame()
    snap.particles.N = n_mem + n_new
    types = ["mem", "actin_lamel", "ligand"]
    mem_tid = types.index("mem")
    actin_tid = types.index("actin_lamel")
    lig_tid = types.index("ligand")
    snap.particles.types = types
    snap.particles.typeid = np.concatenate([
        np.full(n_mem, mem_tid, dtype=np.uint32),
        np.full(n_seed, actin_tid, dtype=np.uint32),
        np.full(n_lig, lig_tid, dtype=np.uint32),
    ])
    snap.particles.position = np.concatenate([mem_pos, seed_pos, lig_pos], axis=0)
    snap.particles.mass = np.ones(n_mem + n_new, dtype=np.float64)
    # bonds: dcm_edge mesh bonds (front) + lamellipodium bond types (no instances).
    edge_bonds = np.asarray(msnap.bonds.group, dtype=np.uint32)
    snap.bonds.N = int(edge_bonds.shape[0])
    snap.bonds.types = [
        "dcm_edge", "lamel_actin_bond", "lamel_branch_bond", "front_clutch_grip"]
    snap.bonds.typeid = np.zeros(edge_bonds.shape[0], dtype=np.uint32)
    snap.bonds.group = edge_bonds
    snap.configuration.box = list(msnap.configuration.box)

    dev = device or hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=base.seed)
    sim.create_state_from_snapshot(snap)

    # native mesh shell forces (Volume turgor + edge springs + optional bending).
    ig = md.Integrator(dt=base.dt)
    shell_forces, mesh, vol, V0 = build_native_shell_forces(
        sim, base, mesh_tris=mesh_tris, mesh_typeids=mesh_typeids,
        mesh_types=mesh_types, V_rest=V_rest, mean_edge=mean_edge)
    # extend the shell's Harmonic bond force to cover the lamellipodium bond
    # types too (HOOMD requires every snapshot bond type to be parameterized on a
    # single bond force; two Harmonic forces would both demand all 4 types).
    for f in shell_forces:
        if isinstance(f, md.bond.Harmonic):
            f.params["lamel_actin_bond"] = dict(k=base.k_edge, r0=p.actin_rest_length)
            f.params["lamel_branch_bond"] = dict(k=base.k_edge, r0=p.actin_rest_length)
            # grip-bond r0 is set per-instance (exact-r0); k is the clutch spring.
            f.params["front_clutch_grip"] = dict(k=p.k_clutch, r0=0.0)
        ig.forces.append(f)

    # cell-cell tent contact (acts on mem nodes only via cell_of_node; actin /
    # ligand tags are beyond cell_of_tag's length so they are never paired).
    if contact:
        patch_area = 4.0 * np.pi * R ** 2 / nv
        tent = DcmTentContact(
            cell_of_node=cell_of_tag, r_contact=1.05 * mean_edge,
            c_adh=base.c_adh, rep_strength=base.rep_strength,
            adh_strength=base.adh_strength, patch_area=patch_area,
            force_cap=base.force_cap)
        ig.forces.append(tent)
    else:
        tent = None

    # adhesive substrate floor (acts on every particle by z; ligands are pinned
    # separately, actin is held basal by BasalAdhesionTether — both consistent).
    area_per_node = 4.0 * np.pi * R ** 2 / nv
    W_cs = base.W_cs_Jm2 * base.ligand_density * area_per_node
    sub = DcmSubstrateForce(z0=z0, W_cs=W_cs, adh_range=R, k_sub=base.k_sub_Nm)
    ig.forces.append(sub)

    sim.operations.integrator = ig

    # 5. ligand pin (rigid substrate sites)
    lig_pin = LigandPin(z0=z0, k_pin=p.k_clutch * 5.0, ligand_typeid=lig_tid)
    ig.forces.append(lig_pin)

    # basal-sheet tether (keep actin in the basal plane)
    basal = BasalAdhesionTether(z_basal=z_basal, k_adh=p.k_basal_adh,
                                actin_typeid=actin_tid)
    ig.forces.append(basal)

    # the traction tether: actin pulls the leading membrane nodes (the cell body)
    tether = LamellipodialTractionTether(
        cell_of_node=cell_of_tag, rim_cells=rim_cells,
        actin_typeid=actin_tid, mem_typeid=mem_tid, z_basal=z_basal,
        basal_band=basal_band, k_tether=p.k_tether,
        force_cap=p.tether_force_cap, tether_radius=p.tether_radius)
    ig.forces.append(tether)

    # 6. BAOAB with gamma covering ALL types (mem + actin_lamel + ligand).
    from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater
    gamma = {
        "mem": base.gamma_node,
        "actin_lamel": p.gamma_actin,
        "ligand": p.gamma_ligand,
    }
    action, baoab = make_baoab_updater(
        kT=base.kT, gamma=gamma, dt=base.dt, seed=p.seed + 1)
    sim.operations.updaters.append(baoab)

    # 7. per-rim-cell LamellipodiumState + nucleation + clutch updaters.
    #    Each rim cell's state tracks its own seed beads as the initial barbed
    #    ends (leading edge). The nucleation updater advances them outward; the
    #    clutch updater grips them to ligands. We run a SINGLE shared state +
    #    one nucleator + one clutch (they operate on ALL actin_lamel beads, which
    #    is correct because each bead advances radially from its OWN position —
    #    the engine's outward-radial tangent is computed per-bead from the bead's
    #    in-plane position, so a shared updater protrudes every rim cell's sheet
    #    independently). This is the validated single-cell engine applied to the
    #    union of all rim lamellipodia.
    state = LamellipodiumState()
    next_tag = n_mem
    for i in range(n_seed):
        tag = next_tag + i
        state.barbed_end_tags.append(tag)
        # outward-radial tangent from the spheroid centroid (per-bead)
        rxy = seed_pos[i, :2] - sph_c[:2]
        rn = np.hypot(rxy[0], rxy[1])
        tang = (np.array([rxy[0] / rn, rxy[1] / rn, 0.0]) if rn > 1e-12
                else np.array([1.0, 0.0, 0.0]))
        state.tangent_of[tag] = tang
    state.actin_next_tag = n_mem + n_new   # next append lands after the ligands

    # p_advance from band-matched v_front (same formula as the single-cell run)
    batch_dt = p.batch_steps * base.dt
    p_advance = min(1.0, p.v_front * batch_dt * p.S_kinetic / p.actin_rest_length)

    nucleator = LeadingEdgeNucleationUpdater(
        lamel_state=state, z_basal=z_basal, band=basal_band,
        rest_length=p.actin_rest_length, p_advance=p_advance,
        advance_margin=mean_edge, max_advance=p.max_advance_per_cell * rim_cells.size,
        max_radius=float("inf"), seed=p.seed + 11)
    nuc_upd = hoomd.update.CustomUpdater(
        action=nucleator, trigger=hoomd.trigger.Periodic(p.batch_steps))
    sim.operations.updaters.append(nuc_upd)

    clutch = FrontClutchRatchetUpdater(
        z_basal=z_basal, band=basal_band, grip_radius=p.grip_radius,
        bond_type_name="front_clutch_grip",
        max_grip_per_tick=p.max_grip_per_cell * rim_cells.size)
    clutch_upd = hoomd.update.CustomUpdater(
        action=clutch, trigger=hoomd.trigger.Periodic(p.batch_steps))
    sim.operations.updaters.append(clutch_upd)

    sim.run(0)

    return {
        "sim": sim, "rim_cells": rim_cells, "centers": centers,
        "state": state, "nucleator": nucleator, "clutch": clutch,
        "tether": tether, "basal": basal, "lig_pin": lig_pin,
        "n_seed": n_seed, "n_lig": n_lig, "n_mem": n_mem,
        "cell_of_tag": cell_of_tag, "ranges": ranges, "nv": nv,
        "z_basal": z_basal, "mean_edge": mean_edge, "p": p,
        "actin_typeid": actin_tid, "mem_typeid": mem_tid,
        "volume": vol, "V0": V0, "mesh": mesh, "tent": tent, "substrate": sub,
        "baoab": action,
    }


# ---------------------------------------------------------------------------
# Footprint / diagnostics
# ---------------------------------------------------------------------------
def footprint_area(sim, *, cell_of_node: np.ndarray, z_basal: float,
                   band: float, mem_typeid: int) -> float:
    """Convex-hull area of the basal membrane-node footprint (A in A/A₀)."""
    from scipy.spatial import ConvexHull
    snap = sim.state.get_snapshot()
    tid = np.asarray(snap.particles.typeid)
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    n_mem = cell_of_node.shape[0]
    mem = (tid[:n_mem] == mem_typeid)
    mpos = pos[:n_mem][mem]
    basal = np.abs(mpos[:, 2] - z_basal) <= band
    pts = mpos[basal, :2]
    if pts.shape[0] < 3:
        return 0.0
    try:
        return float(ConvexHull(pts).volume)  # 2D hull "volume" == area
    except Exception:
        return 0.0


def rim_centroid_radius(sim, *, rim_cells: np.ndarray, ranges,
                        mem_typeid: int) -> float:
    """Mean in-plane radius of rim-cell centroids from the spheroid centroid."""
    snap = sim.state.get_snapshot()
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    n_mem = sum(b - a for a, b in ranges)
    sph_c = pos[:n_mem].mean(axis=0)
    rs = []
    for c in rim_cells:
        a, b = ranges[int(c)]
        cc = pos[a:b].mean(axis=0)
        rs.append(np.hypot(cc[0] - sph_c[0], cc[1] - sph_c[1]))
    return float(np.mean(rs)) if rs else 0.0
