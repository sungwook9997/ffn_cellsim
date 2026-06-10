"""DCM tier — PROLIFERATING deformable-cell monolayer (contact-inhibited spread).

PI pivot (2026-06-11, all-DCM graduation timeline). This module extends the
Deformable-Cell-Model (``cell/dcm.py`` — icosphere elastic shells + exact
triangulated turgor + adhesive substrate) with **mid-run cell division** so a
monolayer's footprint grows not only by mechanical wetting (which saturates near
A/A0 ~ 1.2) but by **multiplying the cell count** — the dominant route real
epithelial sheets spread. Division is **contact-inhibited**: only low-density RIM
cells (few same-sheet neighbours) divide; high-density CORE cells are quiescent.

Architecture for HOOMD's fixed-type-set constraint
--------------------------------------------------
HOOMD cannot add particle TYPES (or, cleanly, tags) after ``create_state``.  So:

* **One shared membrane type ``'dcm_mem'``** for ALL cells (the type set is fixed
  at build — never grows). Cell identity is carried in a mutable per-node integer
  ``cell_of_node`` (-1 = dormant), NOT in the particle type.
* **Pre-allocate a pool of ``N_max`` icosphere shells.** ``N0`` start ACTIVE in a
  2D hex cluster on the adhesive substrate; the remaining ``N_max-N0`` are DORMANT
  — parked far outside the active region with near-zero interaction. No particles
  (or tags) are ever added → BAOAB sees a fixed, dense tag space (no tag-space
  mutation path, regression-clean).
* **Cell-cell interaction is a CUSTOM ``md.force.Custom``** (``CellCellAdhesion``)
  keyed on ``cell_of_node`` — excluded-volume (soft WCA-like core) + cadherin-scale
  adhesive well applied ONLY between nodes of DIFFERENT ACTIVE cells. Intra-cell
  mechanics are the edge bonds; dormant nodes are skipped. (Per-cell-type LJ is
  impossible with one type — this is the mechanistic replacement.)
* **Turgor** reuses ``DcmTurgorForce`` over the ACTIVE cells only; its
  ``faces`` / ``face_cell`` / ``n_cells`` are rebuilt (in place) when a dormant
  cell is activated. Dormant shells carry no faces → no turgor → they stay parked.

Division (``ProliferationUpdater``, slow batched ``Action``)
------------------------------------------------------------
Each tick: compute every active cell's xy-centroid, count OTHER active centroids
within ``neighbour_radius`` = local density. RIM cells (density below a threshold)
are contact-inhibited-released and divide with ``p_div`` per tick: ACTIVATE one
dormant shell, place its icosphere nodes adjacent to the parent (force-free, with a
surface GAP — no LJ-core overlap; adhesion pulls it in gently), outward from the
cluster centroid. Update ``cell_of_node``, the active mask, and the turgor face
groups. CORE (high-density) cells do NOT divide.

Stability (BAOAB int32 guard cost many crashes before)
------------------------------------------------------
dt small (5e-10, drop to 2e-10 if needed); ``gamma`` dict MUST include ``'dcm_mem'``
(the only membrane type). New daughter cells placed with a surface gap (no t=0
overlap). All new bonds/positions force-free at creation. The custom adhesion uses a
moderate, CAPPED soft core so a residual small overlap can never blow up.

Units SI. Defaults are MCF7 literature bands (R 7.5 µm, turgor 133 Pa).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md

# Reuse the read-only DCM primitives.
from ffn_sim.cell.dcm import (
    _KT_310,
    icosphere_mesh,
    DcmTurgorForce,
    DcmSubstrateForce,
    _cluster_centers,
)

MEM_TYPE = "dcm_mem"  # single shared membrane particle type (fixed type set)


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedProlifDCM:
    """Literature-band constants for the proliferating monolayer (trusted)."""

    R_cell: float = 7.5e-6          # m  MCF7 radius
    subdivisions: int = 1          # icosphere level (1 = 42 nodes/cell, for speed)
    turgor_dP0: float = 133.0      # Pa  baseline osmotic turgor (physiological)
    K_vol: float = 5.0e3           # Pa  osmotic bulk modulus
    k_edge: float = 5.0e-4         # N/m membrane edge spring (cortical elasticity)
    gamma_node: float = 3.9e-10    # N·s/m per-node Stokes drag (cortex gamma_b)
    # Adhesion energy densities (J/m²), trusted literature bands.
    W_cs_Jm2: float = 0.012        # J/m² cell-substrate adhesion (spreading driver)
    W_cc_Jm2: float = 0.3e-3       # J/m² cell-cell adhesion (sheet cohesion)
    eps_rep_kT: float = 8.0        # excluded-volume (soft-core) strength [kT]
    # --- proliferation knobs ---
    n0: int = 7                    # initial ACTIVE cells (one 2D-hex cluster)
    n_max: int = 40               # pre-allocated pool (active + dormant)
    spacing_factor: float = 2.3    # initial cell-center spacing = factor·R (gap, no overlap)
    neighbour_factor: float = 2.6  # neighbour search radius = factor·R (local density)
    rim_max_neighbours: int = 4    # ≤ this many active neighbours ⇒ RIM (may divide)
    p_div_per_tick: float = 0.5    # per-tick division prob for an eligible rim cell
    daughter_gap_factor: float = 0.35  # surface gap of daughter = factor·R (force-free)
    # --- substrate / numerics ---
    k_sub_Nm: float | None = None  # N/m substrate spring (None = rigid dish)
    ligand_density: float = 1.0    # relative ECM ligand density (adhesion ∝ this)
    z_substrate: float = 0.0       # m substrate plane
    dormant_park: float = 60.0     # dormant shells parked this many R away (±)
    kT: float = _KT_310
    dt: float = 5.0e-10            # s  (small for BAOAB stability under adhesion)
    seed: int = 7
    extras: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Custom cell-cell interaction (one type → cannot use per-cell-type LJ)
# ---------------------------------------------------------------------------
class CellCellAdhesion(md.force.Custom):
    """Excluded-volume + cadherin adhesion between nodes of DIFFERENT active cells.

    With a single membrane type the per-cell-type LJ trick of ``cell/dcm.py`` is
    impossible, so this force reads the *mutable* ``cell_of_node`` (-1 = dormant)
    and applies, for every close node pair (i, j) belonging to two DIFFERENT
    ACTIVE cells:

      * a SOFT excluded-volume core (capped harmonic push-apart) for r < sigma,
        keeping the cap finite so a residual t=0 overlap can never blow BAOAB up;
      * a cadherin-scale adhesive WELL for sigma ≤ r < r_cut (linear-to-zero at
        r_cut), pulling neighbouring cells into contact and giving the sheet its
        cohesion.

    Intra-cell pairs (same cell) are handled by the edge bonds + turgor and are
    SKIPPED here. Dormant nodes (cell == -1) are skipped entirely. A HOOMD
    neighbour list is rebuilt each call from the live snapshot (active-node subset
    only) — F is small (≤ N_max·42 nodes) so an O(N · k) cell-list pass is cheap
    relative to the batched division cadence.

    Force law (per pair, magnitude along r̂, repulsion +, attraction −):
        r < sigma            : F_rep = k_core·(sigma − r),  capped at k_core·sigma
        sigma ≤ r < r_cut    : F_adh = −F_adh0·(r_cut − r)/(r_cut − sigma)
        r ≥ r_cut            : 0
    where F_adh0 = 2·W_cc / (r_cut − sigma) gives well depth ≈ W_cc per pair.
    """

    def __init__(self, *, cell_of_node: np.ndarray, sigma: float, r_cut: float,
                 eps_rep: float, W_cc: float) -> None:
        super().__init__(aniso=False)
        # Shared, MUTABLE array (the updater edits it in place on division).
        self.cell_of_node = cell_of_node                 # (N,) int, -1 = dormant
        self.sigma = float(sigma)
        self.r_cut = float(r_cut)
        # Soft-core stiffness from a target overlap energy ~ eps_rep over sigma.
        self.k_core = float(2.0 * eps_rep / (sigma ** 2))
        self.F_adh0 = float(2.0 * W_cc / (r_cut - sigma)) if r_cut > sigma else 0.0
        self._build_count = 0

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        n = pos.shape[0]
        perm = np.argsort(tag)                 # local row → global (tag) order
        pos_g = pos[perm]                      # pos_g[global_idx] = position
        cell = self.cell_of_node               # already in global-idx order

        F_g = np.zeros_like(pos_g)
        active = np.where(cell >= 0)[0]
        if active.size > 1:
            apos = pos_g[active]
            acell = cell[active]
            # Pairs within r_cut on the ACTIVE subset (≤ ~1680 nodes), via a
            # scipy KD-tree (vectorised C, exact) — cheap relative to a per-step
            # python loop and called every force eval.
            pairs = _within_cutoff_pairs(apos, self.r_cut)
            if pairs.size:
                ii = pairs[:, 0]
                jj = pairs[:, 1]
                # Only DIFFERENT active cells.
                diff = acell[ii] != acell[jj]
                ii, jj = ii[diff], jj[diff]
                if ii.size:
                    d = apos[ii] - apos[jj]
                    r = np.linalg.norm(d, axis=1)
                    r_safe = np.where(r > 1e-18, r, 1e-18)
                    rhat = d / r_safe[:, None]
                    fmag = np.zeros_like(r)
                    rep = r < self.sigma
                    overlap = np.clip(self.sigma - r[rep], 0.0, self.sigma)
                    fmag[rep] = self.k_core * overlap          # capped at k_core·sigma
                    adh = (~rep) & (r < self.r_cut)
                    fmag[adh] = -self.F_adh0 * (self.r_cut - r[adh]) / (
                        self.r_cut - self.sigma)
                    fvec = fmag[:, None] * rhat                # on node ii (+ = apart)
                    gi = active[ii]
                    gj = active[jj]
                    np.add.at(F_g, gi, fvec)
                    np.add.at(F_g, gj, -fvec)
        self._build_count += 1

        F = np.empty_like(pos)
        F[perm] = F_g
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F
            arr.potential_energy[:] = np.zeros(n, dtype=np.float64)


def _within_cutoff_pairs(pos: np.ndarray, r_cut: float) -> np.ndarray:
    """Return (P,2) array of i<j index pairs with |pos_i − pos_j| < r_cut.

    Uses scipy's ``cKDTree.query_pairs`` (vectorised C; exact) — called on every
    force evaluation, so the per-step cost must stay small. Falls back to an empty
    array for fewer than two points.
    """
    n = pos.shape[0]
    if n < 2:
        return np.empty((0, 2), dtype=np.int64)
    from scipy.spatial import cKDTree
    tree = cKDTree(pos)
    pairs = tree.query_pairs(r=r_cut, output_type="ndarray")  # (P,2), i<j
    if pairs.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    return pairs.astype(np.int64)


# ---------------------------------------------------------------------------
# Proliferation updater (contact-inhibited rim division)
# ---------------------------------------------------------------------------
class ProliferationUpdater(hoomd.custom.Action):
    """Contact-inhibited mid-run division: rim cells divide, core cells quiescent.

    State (shared, mutated in place so the forces see it next step):
      * ``cell_of_node`` (N,) int, -1 = dormant.
      * ``active`` (N_max,) bool per-cell active flag.
      * ``ranges`` list[(a,b)] global-node range of each pre-allocated shell.
    On each act():
      1. read positions; compute active cells' xy-centroids;
      2. RIM = cell on the convex-hull BOUNDARY of the active-cell centroids (xy)
         — a free-edge cell with open space on one side (the physiological
         contact-inhibition trigger). Interior CORE cells are contact-inhibited.
         Hull membership is scale-invariant as the cluster compacts on settling,
         unlike a fixed neighbour-radius count; degenerate (collinear / N<3)
         configs fall back to a neighbour-count ≤ ``rim_max_neighbours`` test.
      3. each rim cell divides with ``p_div`` if a dormant shell is available;
      4. ACTIVATE a dormant shell: place its (force-free) icosphere nodes at the
         parent centroid + an outward offset (away from the cluster centroid)
         at distance ``2R + gap`` (surface gap, no LJ-core overlap), resting on
         the substrate (lowest node at z0); set ``cell_of_node`` + ``active``;
         rebuild the turgor face groups in place.
    """

    def __init__(self, *, p, cell_of_node, active, ranges, verts0, tris0, nv,
                 turgor, adhesion, rng_seed=12345):
        super().__init__()
        self.p = p
        self.cell_of_node = cell_of_node
        self.active = active
        self.ranges = ranges
        self.verts0 = np.asarray(verts0, dtype=np.float64)
        self.tris0 = np.asarray(tris0, dtype=np.int64)
        self.nv = int(nv)
        self.turgor = turgor
        self.adhesion = adhesion
        self.R = float(p.R_cell)
        self.z0 = float(p.z_substrate)
        self.neighbour_radius = float(p.neighbour_factor * p.R_cell)
        self.rim_max = int(p.rim_max_neighbours)
        self.p_div = float(p.p_div_per_tick)
        self.gap = float(p.daughter_gap_factor * p.R_cell)
        self._rng = np.random.default_rng(rng_seed)
        self._sim = None
        self.n_divisions = 0

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    # -- helpers -----------------------------------------------------------
    def _centroids(self, pos_g):
        cents = {}
        for c in np.where(self.active)[0]:
            a, b = self.ranges[c]
            cents[int(c)] = pos_g[a:b].mean(axis=0)
        return cents

    def _rebuild_turgor_faces(self):
        """Rebuild turgor faces/face_cell over ACTIVE cells only, in place."""
        faces = []
        face_cell = []
        active_ids = np.where(self.active)[0]
        # Compact active-cell volume index 0..n_active-1 (so np.add.at is dense).
        remap = {int(c): k for k, c in enumerate(active_ids)}
        for c in active_ids:
            a, _ = self.ranges[int(c)]
            faces.append(self.tris0 + a)
            face_cell.append(np.full(self.tris0.shape[0], remap[int(c)],
                                     dtype=np.int64))
        self.turgor.faces = np.concatenate(faces, axis=0)
        self.turgor.face_cell = np.concatenate(face_cell)
        self.turgor.n_cells = int(active_ids.size)

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim
        with sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        perm = np.argsort(tag)
        pos_g = pos[perm]

        active_ids = np.where(self.active)[0]
        if active_ids.size == 0 or not np.any(~self.active):
            return  # nothing active, or pool exhausted

        cents = self._centroids(pos_g)
        cids = np.array(sorted(cents.keys()))
        carr = np.array([cents[int(c)] for c in cids])  # (na, 3)
        cluster_cen = carr[:, :2].mean(axis=0)
        xy = carr[:, :2]

        # RIM detection = contact inhibition geometry. A cell is a free-edge RIM
        # cell if it lies on the convex-hull BOUNDARY of the active-cell centroids
        # (xy) — i.e. it has open space on one side, the physiological trigger for
        # division. Interior (CORE) cells are fully surrounded → contact-inhibited.
        # Convex-hull membership is scale-invariant (robust to the cluster
        # compacting as it settles, which a fixed neighbour-radius count is not).
        # Fallback for degenerate (collinear / N<3) configs: the local-density
        # neighbour count within ``neighbour_radius``.
        rim = np.zeros(xy.shape[0], dtype=bool)
        hull_ok = False
        if xy.shape[0] >= 3:
            try:
                from scipy.spatial import ConvexHull
                hull = ConvexHull(xy)
                rim[np.unique(hull.vertices)] = True
                hull_ok = True
            except Exception:  # noqa: BLE001 (collinear → QhullError)
                hull_ok = False
        if not hull_ok:
            d2 = np.sum((xy[:, None, :] - xy[None, :, :]) ** 2, axis=2)
            within = (d2 < self.neighbour_radius ** 2)
            np.fill_diagonal(within, False)
            rim = within.sum(axis=1) <= self.rim_max
        divided_any = False
        for k in np.where(rim)[0]:
            if not np.any(~self.active):
                break  # pool exhausted
            if self._rng.random() >= self.p_div:
                continue
            parent = int(cids[k])
            dormant = np.where(~self.active)[0]
            if dormant.size == 0:
                break
            daughter = int(dormant[0])

            # outward direction = parent centroid → away from cluster centroid
            pc = carr[k, :2]
            outward = pc - cluster_cen
            nrm = np.linalg.norm(outward)
            if nrm < 1e-12:
                # core/degenerate: pick a random outward direction
                ang = self._rng.random() * 2 * np.pi
                outward = np.array([np.cos(ang), np.sin(ang)])
                nrm = 1.0
            outward = outward / nrm
            offset = (2.0 * self.R + self.gap) * outward   # surface gap, no overlap
            parent_cen3 = cents[parent]
            new_center = np.array([
                parent_cen3[0] + offset[0],
                parent_cen3[1] + offset[1],
                self.z0 + self.R,                          # rest on substrate
            ])
            # force-free placement of the (undeformed) icosphere nodes
            new_nodes = self.verts0 + new_center
            a, b = self.ranges[daughter]
            pos_g[a:b] = new_nodes
            self.cell_of_node[a:b] = daughter
            self.active[daughter] = True
            self.n_divisions += 1
            divided_any = True

        if divided_any:
            # write the relocated daughter nodes back (force-free), then rebuild
            # the turgor groups so the new cell is pressurised next step.
            pos_back = np.empty_like(pos)
            pos_back[perm] = pos_g
            with sim.state.cpu_local_snapshot as snap:
                snp = np.asarray(snap.particles.position)
                snp[:] = pos_back
            self._rebuild_turgor_faces()


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------
def build_prolif_snapshot(p: ResolvedProlifDCM):
    """Pre-allocate N_max icosphere shells: N0 active in a 2D-hex cluster on the
    substrate, the rest dormant + parked far away. Single membrane type.

    Returns (gsd Frame, cell_of_node (N,), ranges, active (N_max,), mean_edge, nv,
             verts0, tris0).
    """
    import gsd.hoomd

    verts0, edges, tris0 = icosphere_mesh(p.R_cell, p.subdivisions)
    nv = verts0.shape[0]
    mean_edge = float(np.linalg.norm(
        verts0[edges[:, 0]] - verts0[edges[:, 1]], axis=1).mean())

    # Active cluster centers (2D hex monolayer on the substrate).
    spacing = p.spacing_factor * p.R_cell
    centers0 = _cluster_centers(p.n0, spacing, p.z_substrate, p.R_cell, mode="2d")

    all_pos = []
    bond_groups = []
    cell_of_node = []
    ranges = []
    active = np.zeros(p.n_max, dtype=bool)
    tag = 0
    park = p.dormant_park * p.R_cell
    for c in range(p.n_max):
        ranges.append((tag, tag + nv))
        if c < p.n0:
            all_pos.append(verts0 + centers0[c])
            cell_of_node.append(np.full(nv, c, dtype=np.int64))
            active[c] = True
        else:
            # Park dormant shell far outside the active region (so even if a stray
            # interaction were ever evaluated it is out of every cutoff), in a
            # loose grid so dormant shells never overlap each other.
            d = c - p.n0
            px = park + (d % 8) * (3.0 * p.R_cell)
            py = park + (d // 8) * (3.0 * p.R_cell)
            all_pos.append(verts0 + np.array([px, py, park]))
            cell_of_node.append(np.full(nv, -1, dtype=np.int64))  # dormant
        bond_groups.append(edges + tag)
        tag += nv

    pos = np.concatenate(all_pos, axis=0)
    cell_of_node = np.concatenate(cell_of_node)
    bonds = np.concatenate(bond_groups, axis=0)

    # Box must contain the parked dormant shells too.
    box_edge = float(np.abs(pos).max() * 2.0 + 8.0 * p.R_cell)
    snap = gsd.hoomd.Frame()
    snap.particles.N = pos.shape[0]
    snap.particles.types = [MEM_TYPE]
    snap.particles.position = pos
    snap.particles.typeid = np.zeros(pos.shape[0], dtype=np.uint32)
    snap.particles.mass = np.ones(pos.shape[0])
    snap.bonds.N = bonds.shape[0]
    snap.bonds.types = ["dcm_edge"]
    snap.bonds.typeid = np.zeros(bonds.shape[0], dtype=np.uint32)
    snap.bonds.group = bonds.astype(np.uint32)
    snap.configuration.box = [box_edge, box_edge, box_edge, 0, 0, 0]
    return snap, cell_of_node, ranges, active, mean_edge, nv, verts0, tris0


def build_prolif_simulation(p: ResolvedProlifDCM, *, device=None):
    """Assemble the proliferating-DCM monolayer on the BAOAB integrator.

    Forces: edge bonds (cortical elasticity), CellCellAdhesion (custom),
    DcmTurgorForce (active cells only), DcmSubstrateForce (adhesive floor).
    A ProliferationUpdater drives contact-inhibited rim division.
    Returns a dict of handles.
    """
    from ffn_sim.integrator.baoab import make_baoab_updater

    (snap, cell_of_node, ranges, active, mean_edge, nv, verts0, tris0
     ) = build_prolif_snapshot(p)
    dev = device or hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=p.seed)
    sim.create_state_from_snapshot(snap)

    ig = md.Integrator(dt=p.dt)

    # Edge springs (cortical elasticity); r0 = mean edge length.
    bond = md.bond.Harmonic()
    bond.params["dcm_edge"] = dict(k=p.k_edge, r0=mean_edge)
    ig.forces.append(bond)

    # Per-node membrane patch area → adhesion energies per node / pair.
    area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
    W_cc = p.W_cc_Jm2 * area_per_node                       # J per cell-cell pair
    W_cs = p.W_cs_Jm2 * p.ligand_density * area_per_node    # J per substrate node

    # Custom cell-cell interaction (one type → no per-type LJ possible).
    sigma = 0.9 * mean_edge        # node excluded-volume diameter
    r_cut = 2.0 * sigma            # adhesive tail reach
    eps_rep = p.eps_rep_kT * p.kT
    adhesion = CellCellAdhesion(cell_of_node=cell_of_node, sigma=sigma,
                                r_cut=r_cut, eps_rep=eps_rep, W_cc=W_cc)
    ig.forces.append(adhesion)

    # Turgor over the ACTIVE cells only (exact triangulated face-normal pressure).
    V0 = (4.0 / 3.0) * np.pi * p.R_cell ** 3
    faces = []
    face_cell = []
    active_ids = np.where(active)[0]
    remap = {int(c): k for k, c in enumerate(active_ids)}
    for c in active_ids:
        a, _ = ranges[int(c)]
        faces.append(tris0 + a)
        face_cell.append(np.full(tris0.shape[0], remap[int(c)], dtype=np.int64))
    faces = np.concatenate(faces, axis=0)
    face_cell = np.concatenate(face_cell)
    turgor = DcmTurgorForce(faces=faces, face_cell=face_cell,
                            n_cells=int(active_ids.size), V0=V0,
                            turgor_dP0=p.turgor_dP0, K_vol=p.K_vol)
    ig.forces.append(turgor)

    # Adhesive substrate floor (wetting). Reach = R so a whole bottom hemisphere
    # wets; capped-harmonic well = BAOAB-stable. Dormant shells parked far above
    # (z ~ park·R) feel nothing (out of reach), so they stay parked.
    substrate = DcmSubstrateForce(z0=p.z_substrate, W_cs=W_cs,
                                  adh_range=p.R_cell, k_sub=p.k_sub_Nm)
    ig.forces.append(substrate)

    sim.operations.integrator = ig
    sim.run(0)

    gamma = {MEM_TYPE: p.gamma_node}
    action, updater = make_baoab_updater(kT=p.kT, gamma=gamma, dt=p.dt,
                                         seed=p.seed + 1)
    sim.operations.updaters.append(updater)

    prolif = ProliferationUpdater(
        p=p, cell_of_node=cell_of_node, active=active, ranges=ranges,
        verts0=verts0, tris0=tris0, nv=nv, turgor=turgor, adhesion=adhesion,
        rng_seed=p.seed + 2)

    return {
        "sim": sim, "cell_of_node": cell_of_node, "ranges": ranges,
        "active": active, "nv": nv, "mean_edge": mean_edge, "V0": V0,
        "turgor": turgor, "adhesion": adhesion, "substrate": substrate,
        "prolif": prolif, "baoab": action, "p": p, "verts0": verts0,
        "tris0": tris0,
    }
