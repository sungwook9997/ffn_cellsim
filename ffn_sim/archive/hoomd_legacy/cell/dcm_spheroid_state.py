"""DCM tier — 3D spheroid STATE physics: necrosis (3-zone), bulk pressure, junction switch.

PI pivot (2026-06-11, all-DCM graduation timeline). The platform has so far only
*described* three spheroid phenomena; this module makes them physically real on a
proper-morphology 3D Deformable-Cell-Model ball and exposes them to the visualiser:

  1. NECROSIS / 3-ZONE (``NecrosisUpdater``) — each cell's DEPTH below the cluster
     surface ``d = R_cluster - r_cell`` sets its state by the O₂-penetration
     architecture (references/analysis/_necrosis): rim (d<40 µm) PROLIFERATING,
     shell (40–150 µm) QUIESCENT, core (d>150 µm) NECROTIC (division off + dead).
     The core forms DURING a days-long aggregation (modelled as an aggregation
     phase BEFORE spreading), so the spheroid already carries a necrotic core when
     it lands on the substrate and spreads.
  2. BULK PRESSURE (``PressureProbe``) — per-cell compressive-stress proxy =
     neighbour-cell crowding within a contact radius, mapped to the literature
     [0.5, 5] kPa band (references/analysis/_junction). Builds toward the core.
  3. JUNCTION SWITCH (``JunctionSwitchUpdater``) — when a cell's bulk pressure
     exceeds the onset threshold it WEAKENS cell-cell (cadherin) adhesion and
     STRENGTHENS cell-substrate (integrin) adhesion for THAT cell — the
     cadherin→integrin clutch switch that lets unjamming rim cells spread outward.

Architecture reuses ``cell/dcm_prolif.py`` (single shared membrane type +
pre-allocated dormant pool + exact triangulated turgor + capped adhesive substrate)
but (a) seeds a 3D FCC BALL instead of a 2D monolayer, and (b) replaces the global
adhesion constants with PER-CELL modulated versions so the junction switch can act
on individual cells. No existing file is modified.

Stability (BAOAB int32 guard): single shared type, gapped 3D placement (no LJ-core
overlap), capped soft-core + capped substrate well, gamma dict covers the one type,
small dt. Coarse-graining: each shell = a tissue patch, R_patch scaled up so a
~40–60-cell ball spans the necrotic-onset scale (see _necrosis doc + figure caption).

Units SI.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.cell.dcm import _KT_310, icosphere_mesh, DcmTurgorForce, _cluster_centers
from ffn_sim.archive.hoomd_legacy.cell.dcm_prolif import MEM_TYPE, _within_cutoff_pairs

_UM = 1.0e6


class CellState(IntEnum):
    """Per-cell radial zone state (necrosis 3-zone)."""

    PROLIFERATING = 0   # rim, O2/nutrient replete, may divide
    QUIESCENT = 1       # viable hypoxic shell, no division
    NECROTIC = 2        # anoxic core, dead
    DORMANT = -1        # unused pool shell (not part of the spheroid)


# ---------------------------------------------------------------------------
# Resolved parameters (literature bands, trusted; see references/analysis/_*)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedSpheroidState:
    """3D-spheroid state-physics constants. COARSE-GRAINED patch scale.

    R_patch is scaled up from the 7.5-µm single-cell radius so a ~40–60-cell ball
    spans the necrotic-onset diameter (~300–500 µm) on CPU — each shell is a tissue
    patch. Depth thresholds stay at their literature values (40 / 150 µm).
    """

    # --- coarse-grained morphology ---
    R_patch: float = 30.0e-6        # m  coarse-grained patch radius (= "cell" here)
    subdivisions: int = 1           # icosphere level (1 = 42 nodes/cell, for speed)
    n_cells: int = 55               # cells in the 3D ball (≤ pool)
    spacing_factor: float = 2.05    # initial centre spacing = factor·R (gapped)

    # --- membrane / turgor mechanics (scaled patch) ---
    turgor_dP0: float = 133.0       # Pa baseline osmotic turgor
    K_vol: float = 5.0e3            # Pa osmotic bulk modulus
    k_edge: float = 5.0e-4          # N/m membrane edge spring
    gamma_node: float = 3.9e-10     # N·s/m per-node Stokes drag
    eps_rep_kT: float = 8.0         # excluded-volume soft-core strength [kT]

    # --- adhesion energy densities (J/m², literature bands) ---
    W_cs_Jm2: float = 0.05          # J/m² cell-substrate adhesion (spreading driver)
    W_cc_Jm2: float = 0.3e-3        # J/m² cell-cell adhesion (spheroid cohesion)

    # --- 3-zone necrosis depth thresholds (µm; _necrosis doc) ---
    d_prolif_um: float = 40.0       # d < this -> PROLIFERATING (rim)
    d_necrotic_um: float = 150.0    # d >= this -> NECROTIC (fresh aggregate threshold)
    d_necrotic_mature_um: float = 90.0  # mature-spheroid necrotic threshold (deeper core)

    # --- bulk pressure proxy (crowding -> kPa; _junction doc) ---
    contact_factor: float = 2.2     # neighbour contact radius = factor·R_patch
    crowd_lo: float = 4.0           # crowd at/below -> P_min
    crowd_hi: float = 12.0          # crowd at/above -> P_max (full FCC shell)
    P_min_kPa: float = 0.0
    P_max_kPa: float = 6.0

    # --- junction switch (cadherin -> integrin) ---
    P_switch_kPa: float = 0.5       # bulk-pressure onset for the switch
    cadherin_weak_factor: float = 0.3   # cell-cell W_cc multiplier when switched
    integrin_strong_factor: float = 3.0  # cell-substrate W_cs multiplier when switched

    # --- substrate / numerics ---
    z_substrate: float = 0.0
    dormant_park: float = 80.0      # dormant shells parked this many R away
    kT: float = _KT_310
    dt: float = 3.0e-10             # s (small for BAOAB stability under adhesion)
    seed: int = 7
    extras: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Per-cell modulated cell-cell adhesion (junction switch weakens cadherin per cell)
# ---------------------------------------------------------------------------
class ModulatedCellCellAdhesion(md.force.Custom):
    """Excluded-volume + cadherin adhesion between nodes of DIFFERENT active cells,
    with a PER-CELL adhesion multiplier so the junction switch can weaken the
    cadherin (cell-cell) adhesion of individual cells.

    Same soft-core + linear-well law as ``dcm_prolif.CellCellAdhesion`` but the
    adhesive well depth of a pair (i, j) is scaled by the geometric mean of the two
    cells' ``cad_mult`` entries (1.0 = full cadherin; <1 = switched/weakened). The
    repulsive core is NOT scaled (excluded volume always holds → BAOAB-safe).
    """

    def __init__(self, *, cell_of_node: np.ndarray, cad_mult: np.ndarray,
                 sigma: float, r_cut: float, eps_rep: float, W_cc: float) -> None:
        super().__init__(aniso=False)
        self.cell_of_node = cell_of_node          # (N,) int, -1 = dormant
        self.cad_mult = cad_mult                  # (n_max,) per-cell cadherin factor
        self.sigma = float(sigma)
        self.r_cut = float(r_cut)
        self.k_core = float(2.0 * eps_rep / (sigma ** 2))
        self.F_adh0 = float(2.0 * W_cc / (r_cut - sigma)) if r_cut > sigma else 0.0

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        n = pos.shape[0]
        perm = np.argsort(tag)
        pos_g = pos[perm]
        cell = self.cell_of_node

        F_g = np.zeros_like(pos_g)
        active = np.where(cell >= 0)[0]
        if active.size > 1:
            apos = pos_g[active]
            acell = cell[active]
            pairs = _within_cutoff_pairs(apos, self.r_cut)
            if pairs.size:
                ii, jj = pairs[:, 0], pairs[:, 1]
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
                    fmag[rep] = self.k_core * overlap
                    adh = (~rep) & (r < self.r_cut)
                    # per-pair cadherin multiplier = geomean of the two cells'.
                    ci = acell[ii][adh]
                    cj = acell[jj][adh]
                    mult = np.sqrt(self.cad_mult[ci] * self.cad_mult[cj])
                    fmag[adh] = -self.F_adh0 * mult * (
                        self.r_cut - r[adh]) / (self.r_cut - self.sigma)
                    fvec = fmag[:, None] * rhat
                    np.add.at(F_g, active[ii], fvec)
                    np.add.at(F_g, active[jj], -fvec)

        F = np.empty_like(pos)
        F[perm] = F_g
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F
            arr.potential_energy[:] = np.zeros(n, dtype=np.float64)


# ---------------------------------------------------------------------------
# Per-node modulated adhesive substrate (junction switch strengthens integrin per cell)
# ---------------------------------------------------------------------------
class ModulatedSubstrateForce(md.force.Custom):
    """Capped-harmonic adhesive substrate well with a PER-CELL integrin multiplier.

    Identical capped-harmonic law to ``dcm.DcmSubstrateForce`` but the well depth
    (and hence stiffness) of each node is scaled by its cell's ``int_mult`` entry
    (1.0 = baseline integrin; >1 = switched/strengthened). The cap is recomputed
    per node so the steric push-up stays finite (BAOAB-safe) at any multiplier.
    """

    def __init__(self, *, z0: float, W_cs: float, adh_range: float,
                 cell_of_node: np.ndarray, int_mult: np.ndarray) -> None:
        super().__init__(aniso=False)
        self.z0 = float(z0)
        self.W_cs = float(W_cs)
        self.rng = float(adh_range)
        self.k_well0 = 2.0 * self.W_cs / (self.rng ** 2)   # baseline stiffness
        self.cell_of_node = cell_of_node
        self.int_mult = int_mult                            # (n_max,) per-cell

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        n = pos.shape[0]
        perm = np.argsort(tag)
        pos_g = pos[perm]
        cell = self.cell_of_node
        mult = np.ones(n)
        act = cell >= 0
        mult[act] = self.int_mult[cell[act]]
        k = self.k_well0 * mult                              # per-node stiffness

        dz = pos_g[:, 2] - self.z0
        F_g = np.zeros_like(pos_g)
        U = np.zeros(n, dtype=np.float64)
        within = np.abs(dz) <= self.rng
        F_g[within, 2] = -k[within] * dz[within]
        U[within] = 0.5 * k[within] * dz[within] ** 2 - self.W_cs * mult[within]
        deep = dz < -self.rng
        F_g[deep, 2] = k[deep] * self.rng
        U[deep] = k[deep] * self.rng * (-(dz[deep]) - 0.5 * self.rng) \
            - self.W_cs * mult[deep]

        F = np.empty_like(pos)
        F[perm] = F_g
        Uo = np.empty(n)
        Uo[perm] = U
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F
            arr.potential_energy[:] = Uo


# ---------------------------------------------------------------------------
# State helpers (shared, read by the updaters + the visualiser)
# ---------------------------------------------------------------------------
class SpheroidStateArrays:
    """Mutable per-cell state shared between forces, updaters, and the driver.

    All arrays are length n_max and indexed by cell id. The updaters mutate them in
    place so the modulated forces see the new values next step; the driver/visualiser
    reads snapshots of them over time.
    """

    def __init__(self, n_max: int):
        self.state = np.full(n_max, int(CellState.DORMANT), dtype=np.int64)
        self.depth_um = np.zeros(n_max)          # d = R_cluster - r_cell  [µm]
        self.pressure_kPa = np.zeros(n_max)      # bulk compressive proxy [kPa]
        self.switched = np.zeros(n_max, dtype=bool)
        self.cad_mult = np.ones(n_max)           # cell-cell adhesion multiplier
        self.int_mult = np.ones(n_max)           # cell-substrate adhesion multiplier


def _active_centroids(pos_g, ranges, active_ids):
    """Return (ids array, centroids (k,3)) for the given active cell ids."""
    cents = np.array([pos_g[ranges[int(c)][0]:ranges[int(c)][1]].mean(axis=0)
                      for c in active_ids])
    return np.asarray(active_ids), cents


# ---------------------------------------------------------------------------
# Necrosis 3-zone updater (depth below the cluster surface)
# ---------------------------------------------------------------------------
class NecrosisUpdater(hoomd.custom.Action):
    """Assign each cell a 3-zone state from its DEPTH below the spheroid surface.

    Each act(): compute active centroids, the cluster centroid, each cell's radius
    r_cell from the centroid, and the cluster surface radius R_cluster (max r_cell +
    one patch radius). Depth d = R_cluster - r_cell. State by depth (literature
    bands): d < d_prolif -> PROLIFERATING; d_prolif <= d < d_necrotic -> QUIESCENT;
    d >= d_necrotic -> NECROTIC. Necrotic cells have division switched off (their
    state gates the ProliferationUpdater) and are marked dead. State is monotone
    toward death within a run: once NECROTIC a cell stays NECROTIC (irreversible).
    """

    def __init__(self, *, p: ResolvedSpheroidState, st: SpheroidStateArrays,
                 ranges, active: np.ndarray):
        super().__init__()
        self.p = p
        self.st = st
        self.ranges = ranges
        self.active = active
        self.R = float(p.R_patch)
        # Maturation ramp (0->1): as the spheroid ages during aggregation, nutrient
        # depletion deepens hypoxia → the effective necrotic depth threshold SHRINKS
        # from d_necrotic toward d_necrotic_mature, so the necrotic core grows
        # gradually over the aggregation phase (PI: necrosis develops over the
        # days-long aggregation, not instantly). Driven by the driver via
        # set_maturation(); 0 = fresh aggregate, 1 = mature spheroid.
        self.maturation = 0.0
        # Once the spheroid lands on the substrate it flattens, so the radial-depth
        # metric (built for a sphere) no longer describes the tissue interior. The
        # 3-zone assignment is therefore FROZEN at the aggregation→spreading
        # transition: the necrotic core / quiescent shell / proliferating rim are
        # carried forward as-is (state latched), and depth is no longer recomputed.
        self.frozen = False
        self._sim = None

    def set_maturation(self, m: float) -> None:
        self.maturation = float(np.clip(m, 0.0, 1.0))

    def freeze(self) -> None:
        self.frozen = True

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        if self.frozen:
            return  # 3-zone state latched at the spreading transition (see __init__)
        sim = self._sim
        with sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        pos_g = pos[np.argsort(tag)]
        active_ids = np.where(self.active)[0]
        if active_ids.size == 0:
            return
        ids, cents = _active_centroids(pos_g, self.ranges, active_ids)
        cluster_cen = cents.mean(axis=0)
        r_cell = np.linalg.norm(cents - cluster_cen, axis=1)
        # Surface = outermost cell CENTROID radius (so the outermost cells sit at
        # depth 0 = the proliferating rim; depth grows inward toward the core).
        R_cluster = r_cell.max()
        depth = (R_cluster - r_cell)               # m, >= 0
        depth_um = depth * _UM
        d_pro = self.p.d_prolif_um
        # Mature spheroid: necrotic threshold shrinks toward d_necrotic_mature so the
        # anoxic core grows during aggregation. m=0 -> d_necrotic (fresh); m=1 ->
        # d_necrotic_mature (deep core). Linear interpolation in the maturation knob.
        d_nec = (self.p.d_necrotic_um
                 + self.maturation
                 * (self.p.d_necrotic_mature_um - self.p.d_necrotic_um))
        for k, c in enumerate(ids):
            c = int(c)
            self.st.depth_um[c] = depth_um[k]
            if self.st.state[c] == int(CellState.NECROTIC):
                continue                           # irreversible
            if depth_um[k] >= d_nec:
                self.st.state[c] = int(CellState.NECROTIC)
            elif depth_um[k] < d_pro:
                self.st.state[c] = int(CellState.PROLIFERATING)
            else:
                self.st.state[c] = int(CellState.QUIESCENT)


# ---------------------------------------------------------------------------
# Bulk pressure probe (neighbour crowding -> kPa band)
# ---------------------------------------------------------------------------
class PressureProbe(hoomd.custom.Action):
    """Per-cell bulk compressive-stress proxy from neighbour-cell crowding.

    Each act(): count, for every active cell, the number of OTHER active cells whose
    centroid is within ``contact_factor·R_patch``; map that crowd count to the
    literature [P_min, P_max] kPa band by a clipped linear ramp between crowd_lo and
    crowd_hi. Core cells (full coordination shell) saturate near P_max; free-surface
    rim cells sit near P_min — so pressure builds toward the core, as observed.
    """

    def __init__(self, *, p: ResolvedSpheroidState, st: SpheroidStateArrays,
                 ranges, active: np.ndarray):
        super().__init__()
        self.p = p
        self.st = st
        self.ranges = ranges
        self.active = active
        self.r_contact = float(p.contact_factor * p.R_patch)
        self._sim = None

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim
        with sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        pos_g = pos[np.argsort(tag)]
        active_ids = np.where(self.active)[0]
        if active_ids.size == 0:
            return
        ids, cents = _active_centroids(pos_g, self.ranges, active_ids)
        d2 = np.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
        within = d2 < self.r_contact ** 2
        np.fill_diagonal(within, False)
        crowd = within.sum(axis=1).astype(float)
        frac = np.clip((crowd - self.p.crowd_lo)
                       / (self.p.crowd_hi - self.p.crowd_lo), 0.0, 1.0)
        P = self.p.P_min_kPa + (self.p.P_max_kPa - self.p.P_min_kPa) * frac
        for k, c in enumerate(ids):
            self.st.pressure_kPa[int(c)] = P[k]


# ---------------------------------------------------------------------------
# Junction switch (cadherin -> integrin clutch) per cell
# ---------------------------------------------------------------------------
class JunctionSwitchUpdater(hoomd.custom.Action):
    """Flip the cadherin->integrin clutch for cells above the pressure threshold.

    Each act(): for every active, NON-necrotic cell whose bulk pressure exceeds
    ``P_switch_kPa``, set ``switched=True``, WEAKEN its cadherin (cad_mult ->
    cadherin_weak_factor) and STRENGTHEN its integrin (int_mult ->
    integrin_strong_factor). The change is read live by the modulated forces, so a
    switched cell loses cell-cell cohesion and gains substrate traction → high-
    pressure rim cells unjam and spread outward. The switch is latched (does not
    revert) within a run; necrotic core cells are excluded (their switch is moot).
    """

    def __init__(self, *, p: ResolvedSpheroidState, st: SpheroidStateArrays,
                 active: np.ndarray):
        super().__init__()
        self.p = p
        self.st = st
        self.active = active
        self._sim = None

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        for c in np.where(self.active)[0]:
            c = int(c)
            if self.st.state[c] == int(CellState.NECROTIC):
                continue
            if self.st.switched[c]:
                continue
            if self.st.pressure_kPa[c] > self.p.P_switch_kPa:
                self.st.switched[c] = True
                self.st.cad_mult[c] = self.p.cadherin_weak_factor
                self.st.int_mult[c] = self.p.integrin_strong_factor


# ---------------------------------------------------------------------------
# Necrosis-gated rim proliferation (division only from a PROLIFERATING rim cell)
# ---------------------------------------------------------------------------
class GatedProliferationUpdater(hoomd.custom.Action):
    """Contact-inhibited rim division GATED by the necrosis 3-zone state.

    A cell divides only if it is (a) ACTIVE, (b) state == PROLIFERATING (rim, by
    depth) and (c) a free-edge cell (on the 3D convex hull of active centroids, the
    contact-inhibition geometry). The daughter is an undeformed icosphere placed
    force-free outward from the cluster centroid at 2R+gap (surface gap, no LJ-core
    overlap), then the turgor faces are rebuilt. Necrotic/quiescent cells never
    divide. Used in the SPREADING phase so the proliferating rim drives footprint
    growth while the core stays dead.
    """

    def __init__(self, *, p: ResolvedSpheroidState, st: SpheroidStateArrays,
                 cell_of_node, active, ranges, verts0, tris0, nv, turgor,
                 p_div=0.4, gap_factor=0.35, rng_seed=99):
        super().__init__()
        self.p = p
        self.st = st
        self.cell_of_node = cell_of_node
        self.active = active
        self.ranges = ranges
        self.verts0 = np.asarray(verts0, dtype=np.float64)
        self.tris0 = np.asarray(tris0, dtype=np.int64)
        self.nv = int(nv)
        self.turgor = turgor
        self.R = float(p.R_patch)
        self.z0 = float(p.z_substrate)
        self.p_div = float(p_div)
        self.gap = float(gap_factor * p.R_patch)
        self._rng = np.random.default_rng(rng_seed)
        self._sim = None
        self.n_divisions = 0

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def _rebuild_turgor_faces(self):
        faces, face_cell = [], []
        active_ids = np.where(self.active)[0]
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
        if not np.any(~self.active):
            return  # pool exhausted
        with sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        perm = np.argsort(tag)
        pos_g = pos[perm]
        active_ids = np.where(self.active)[0]
        if active_ids.size < 4:
            return
        ids, cents = _active_centroids(pos_g, self.ranges, active_ids)
        cluster_cen = cents.mean(axis=0)
        # free-edge cells = 3D convex-hull vertices.
        rim = np.zeros(ids.shape[0], dtype=bool)
        try:
            from scipy.spatial import ConvexHull
            hull = ConvexHull(cents)
            rim[np.unique(hull.vertices)] = True
        except Exception:  # noqa: BLE001
            return
        divided = False
        for k in np.where(rim)[0]:
            c = int(ids[k])
            if self.st.state[c] != int(CellState.PROLIFERATING):
                continue
            if not np.any(~self.active):
                break
            if self._rng.random() >= self.p_div:
                continue
            daughter = int(np.where(~self.active)[0][0])
            outward = cents[k] - cluster_cen
            nrm = np.linalg.norm(outward)
            if nrm < 1e-12:
                outward = self._rng.standard_normal(3)
                nrm = np.linalg.norm(outward)
            outward = outward / nrm
            new_center = cents[k] + (2.0 * self.R + self.gap) * outward
            new_center[2] = max(new_center[2], self.z0 + self.R)
            a, b = self.ranges[daughter]
            pos_g[a:b] = self.verts0 + new_center
            self.cell_of_node[a:b] = daughter
            self.active[daughter] = True
            self.st.state[daughter] = int(CellState.PROLIFERATING)
            self.st.cad_mult[daughter] = 1.0
            self.st.int_mult[daughter] = 1.0
            self.n_divisions += 1
            divided = True
        if divided:
            pos_back = np.empty_like(pos)
            pos_back[perm] = pos_g
            with sim.state.cpu_local_snapshot as snap:
                np.asarray(snap.particles.position)[:] = pos_back
            self._rebuild_turgor_faces()


# ---------------------------------------------------------------------------
# Builder — 3D ball, single shared membrane type, pre-allocated pool
# ---------------------------------------------------------------------------
def build_spheroid_snapshot(p: ResolvedSpheroidState, n_max: int):
    """Pre-allocate n_max shells: n_cells active in a 3D FCC ball, rest dormant.

    The ball is built ABOVE the substrate (so the aggregation phase compacts it
    free of the floor); the driver lowers it onto the substrate for the spreading
    phase. Returns (Frame, cell_of_node, ranges, active, mean_edge, nv, verts0, tris0).
    """
    import gsd.hoomd

    verts0, edges, tris0 = icosphere_mesh(p.R_patch, p.subdivisions)
    nv = verts0.shape[0]
    mean_edge = float(np.linalg.norm(
        verts0[edges[:, 0]] - verts0[edges[:, 1]], axis=1).mean())

    spacing = p.spacing_factor * p.R_patch
    # 3D ball, lifted so its base sits a few R above the substrate during aggregation.
    centers0 = _cluster_centers(p.n_cells, spacing, p.z_substrate, p.R_patch,
                                mode="3d")
    centers0[:, 2] += 2.5 * p.R_patch     # lift a little clear of the floor for aggregation

    all_pos, bond_groups, cell_of_node, ranges = [], [], [], []
    active = np.zeros(n_max, dtype=bool)
    tag = 0
    park = p.dormant_park * p.R_patch
    for c in range(n_max):
        ranges.append((tag, tag + nv))
        if c < p.n_cells:
            all_pos.append(verts0 + centers0[c])
            cell_of_node.append(np.full(nv, c, dtype=np.int64))
            active[c] = True
        else:
            d = c - p.n_cells
            px = park + (d % 8) * (3.0 * p.R_patch)
            py = park + (d // 8) * (3.0 * p.R_patch)
            all_pos.append(verts0 + np.array([px, py, park]))
            cell_of_node.append(np.full(nv, -1, dtype=np.int64))
        bond_groups.append(edges + tag)
        tag += nv

    pos = np.concatenate(all_pos, axis=0)
    cell_of_node = np.concatenate(cell_of_node)
    bonds = np.concatenate(bond_groups, axis=0)
    box_edge = float(np.abs(pos).max() * 2.0 + 8.0 * p.R_patch)
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


def build_spheroid_simulation(p: ResolvedSpheroidState, n_max: int, *, device=None):
    """Assemble the 3D-spheroid state simulation on the BAOAB integrator.

    Forces: edge bonds, ModulatedCellCellAdhesion (per-cell cadherin), turgor,
    ModulatedSubstrateForce (per-cell integrin). State updaters (necrosis, pressure,
    junction switch, gated proliferation) are returned UNATTACHED so the driver can
    schedule them per phase. Returns a dict of handles.
    """
    from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater

    (snap, cell_of_node, ranges, active, mean_edge, nv, verts0, tris0
     ) = build_spheroid_snapshot(p, n_max)
    dev = device or hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=p.seed)
    sim.create_state_from_snapshot(snap)

    st = SpheroidStateArrays(n_max)
    # initial state: all active cells QUIESCENT (zone assigned on first necrosis tick).
    st.state[np.where(active)[0]] = int(CellState.QUIESCENT)

    ig = md.Integrator(dt=p.dt)
    bond = md.bond.Harmonic()
    bond.params["dcm_edge"] = dict(k=p.k_edge, r0=mean_edge)
    ig.forces.append(bond)

    area_per_node = 4.0 * np.pi * p.R_patch ** 2 / nv
    W_cc = p.W_cc_Jm2 * area_per_node
    W_cs = p.W_cs_Jm2 * area_per_node

    sigma = 0.9 * mean_edge
    r_cut = 2.0 * sigma
    eps_rep = p.eps_rep_kT * p.kT
    adhesion = ModulatedCellCellAdhesion(
        cell_of_node=cell_of_node, cad_mult=st.cad_mult, sigma=sigma,
        r_cut=r_cut, eps_rep=eps_rep, W_cc=W_cc)
    ig.forces.append(adhesion)

    V0 = (4.0 / 3.0) * np.pi * p.R_patch ** 3
    faces, face_cell = [], []
    active_ids = np.where(active)[0]
    remap = {int(c): k for k, c in enumerate(active_ids)}
    for c in active_ids:
        a, _ = ranges[int(c)]
        faces.append(tris0 + a)
        face_cell.append(np.full(tris0.shape[0], remap[int(c)], dtype=np.int64))
    turgor = DcmTurgorForce(faces=np.concatenate(faces, axis=0),
                            face_cell=np.concatenate(face_cell),
                            n_cells=int(active_ids.size), V0=V0,
                            turgor_dP0=p.turgor_dP0, K_vol=p.K_vol)
    ig.forces.append(turgor)

    substrate = ModulatedSubstrateForce(
        z0=p.z_substrate, W_cs=W_cs, adh_range=p.R_patch,
        cell_of_node=cell_of_node, int_mult=st.int_mult)
    ig.forces.append(substrate)

    sim.operations.integrator = ig
    sim.run(0)

    gamma = {MEM_TYPE: p.gamma_node}
    baoab_action, baoab_updater = make_baoab_updater(
        kT=p.kT, gamma=gamma, dt=p.dt, seed=p.seed + 1)
    sim.operations.updaters.append(baoab_updater)

    necrosis = NecrosisUpdater(p=p, st=st, ranges=ranges, active=active)
    pressure = PressureProbe(p=p, st=st, ranges=ranges, active=active)
    junction = JunctionSwitchUpdater(p=p, st=st, active=active)
    prolif = GatedProliferationUpdater(
        p=p, st=st, cell_of_node=cell_of_node, active=active, ranges=ranges,
        verts0=verts0, tris0=tris0, nv=nv, turgor=turgor, rng_seed=p.seed + 3)

    return {
        "sim": sim, "st": st, "cell_of_node": cell_of_node, "ranges": ranges,
        "active": active, "nv": nv, "mean_edge": mean_edge, "V0": V0,
        "turgor": turgor, "adhesion": adhesion, "substrate": substrate,
        "necrosis": necrosis, "pressure": pressure, "junction": junction,
        "prolif": prolif, "baoab": baoab_action, "p": p,
        "verts0": verts0, "tris0": tris0,
    }


def lower_to_substrate(sim, ranges, active, R_patch, z0):
    """Translate the whole active cluster down so its lowest node rests at z0+ (gap).

    Called at the aggregation→spreading transition: the compacted ball is dropped
    onto the adhesive substrate (a small gap is kept so the capped well, not a
    t=0 overlap, pulls it down → BAOAB-safe).
    """
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    perm = np.argsort(tag)
    pos_g = pos[perm]
    act_nodes = np.concatenate([np.arange(ranges[int(c)][0], ranges[int(c)][1])
                                for c in np.where(active)[0]])
    zmin = pos_g[act_nodes, 2].min()
    dz = (z0 + 0.5 * R_patch) - zmin          # leave a half-R gap above the floor
    pos_g[act_nodes, 2] += dz
    pos_back = np.empty_like(pos)
    pos_back[perm] = pos_g
    with sim.state.cpu_local_snapshot as snap:
        np.asarray(snap.particles.position)[:] = pos_back
