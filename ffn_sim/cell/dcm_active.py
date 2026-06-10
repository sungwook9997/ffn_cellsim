"""DCM many-cell spheroid — ACTIVE spreading drivers on the native+tent stack.

PI (2026-06-11): the validated native-mesh + bilinear-tent spheroid capstone
(`cell/dcm_native_shell.py`, A/A0=2.36 at 14 cells) is *passive* wetting — turgor
flattens the ball, the adhesive floor pins the base, contact holds cells together.
The ACTIVE spreading machinery the PI specified is NOT in that stack yet. This
module adds it, ADDITIVELY, on the validated stack (no existing file touched):

  (#5 in the master stack) **ACTIVE RIM TRACTION** — a coarse-grained
  lamellipodium + actomyosin contraction-belt + FA molecular-clutch engine at the
  DCM CELL scale. Basal RIM cells (spheroid periphery AND in substrate contact)
  exert an ACTIVE OUTWARD (radial-from-centroid, in-plane) traction on their basal
  nodes — the integrated lamellipodial protrusion + clutch grip pulling the rim
  out. This is an energy INPUT (active), not passive wetting. Optionally a weak
  inward apical contraction-belt line tension drags the cell body after the edge.

  (#6 in the master stack) **BULK-PRESSURE JUNCTION SWITCH** — wires the validated
  `dcm_spheroid_state` updaters (PressureProbe -> JunctionSwitchUpdater) onto the
  tent contact's `cad_mult` seam + the active traction's per-cell integrin gain:
  a cell whose crowding-pressure proxy exceeds P_switch (0.5 kPa onset) WEAKENS its
  cadherin (cad_mult down on DcmTentContact) and STRENGTHENS its integrin (active
  traction + substrate gain up) — the cadherin->integrin clutch switch that lets
  pressure-loaded rim cells unjam and spread (the pressure-release mechanism).

  **LIVE PROLIFERATION** — rim-biased contact-inhibited division at low cadence on
  the pre-allocated native-shell cell pool (no mesh split: a dormant pool shell is
  *activated* as the daughter, its triangles already live in the fixed mesh-tag
  space; the native Volume force conserves it once its V0 is set).

WHY COARSE-GRAIN THE TRACTION (vs the bead-resolved `spreading_drive.py` engine):
  the single-cell engine mutates particle/bond TOPOLOGY (appends actin beads,
  adds clutch grip bonds) which is incompatible with the many-cell native mesh
  (one fixed `hoomd.mesh.Mesh` tag-space shared by every cell's Volume force —
  appending beads renumbers tags and breaks every cell's triangle list). So the
  active traction is coarse-grained to a per-NODE body force on the existing mesh
  nodes (no topology mutation), with its MAGNITUDE calibrated against the
  single-cell engine's FA-clutch traction (see `ActiveRimTraction` docstring).

Stability (BAOAB int32 guard, the failure mode that bit before): the active
traction is CAPPED per node (``f_cap``) and RAMPED in over ``ramp_steps`` so it
never shocks the mesh; division activates a force-free gapped daughter; gamma dict
covers the single membrane type; small dt. Build -> run(0) -> run(2000) finite
BEFORE any long run (the smoke does this).

Units SI.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md

from ffn_sim.cell.dcm_native_shell import (
    ResolvedNativeDCM,
    build_native_dcm_simulation,
    native_icosphere,
    mesh_volume_of,
)
from ffn_sim.cell.dcm_spheroid_state import (
    CellState,
    SpheroidStateArrays,
    ResolvedSpheroidState,
    PressureProbe,
    JunctionSwitchUpdater,
    NecrosisUpdater,
)

_UM = 1.0e6


# ---------------------------------------------------------------------------
# Resolved parameters for the ACTIVE spheroid (extends the native-shell bands)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedActiveSpheroid:
    """Active-spreading DCM spheroid constants (native+tent bands + active layer).

    The mechanical base (R_cell, K_bulk, edge/turgor, tent contact, substrate) is
    the VALIDATED native-shell band (``ResolvedNativeDCM`` defaults). The added
    fields are the active-traction / junction-switch / proliferation layer.
    """

    # --- validated native-shell base (passed through to ResolvedNativeDCM) ---
    R_cell: float = 7.5e-6
    subdivisions: int = 2
    K_bulk: float = 2500.0
    turgor_inflate: float = 1.05
    k_edge: float = 1.0e-3
    gamma_node: float = 3.9e-10
    W_cs_Jm2: float = 0.5e-3
    rep_strength: float = 1.0e8
    adh_strength: float = 1.0e8
    c_adh: float = 5.0e-7
    contact_force_cap: float = 5.0e-8
    spacing_factor: float = 2.3
    dt: float = 3.0e-10
    seed: int = 7

    # --- ACTIVE RIM TRACTION (#5) ---------------------------------------------
    # f_act: per-basal-node outward active traction [N]. CALIBRATION (vs the
    # single-cell spreading_drive engine): that engine's FA molecular clutch
    # carries k_int=1e-3 N/m at clutch extensions ~50-200 nm => per-clutch
    # traction 5e-11..2e-10 N. The DCM substrate adhesive well's max pull per node
    # is ~6e-10 N (W_cs=0.5 mJ/m^2 over a 162-node shell). We set f_act=1.2e-10 N
    # so a rim node's active outward pull is (a) in the single-cell clutch-traction
    # band and (b) ~0.2x the passive substrate well force => it is a genuine
    # OVER-DRIVE of passive wetting (advances the rim beyond the wetting plateau)
    # while staying ~400x below the contact force_cap => BAOAB-safe.
    f_act: float = 1.2e-10
    f_cap: float = 6.0e-10          # N per-node active-force cap (BAOAB guard)
    ramp_steps: int = 4000          # steps to ramp f_act 0 -> full (no shock)
    rim_contact_band: float = 0.5   # basal node if z < z0 + band*R (in contact)
    rim_neighbour_factor: float = 2.6  # rim cell if < this*R from few neighbours
    rim_max_neighbours: int = 9     # a rim cell has <= this many close neighbours
    belt_factor: float = 0.25       # apical contraction-belt tension = belt*f_act (inward)
    integrin_switch_gain: float = 3.0  # switched cells' active traction x this
    # --- WHOLE-CELL MIGRATION (centroid translocation, not just basal splay) ---
    # migrate_factor: the per-NODE outward body force applied to ALL of a rim
    # cell's nodes = migrate_factor * f_act. UNLIKE the basal lamellipodial pull
    # (which deforms/splays the basal lip but, once the inward belt is summed over
    # the ~131 apical nodes, gives a near-ZERO net cell force => the centroid does
    # NOT translocate — the diagnosis the PI flagged), this term is a genuine NET
    # outward force on the whole cell body: the lamellipodium PULLS THE CELL BODY
    # forward (the molecular clutch transmits traction through the cytoskeleton to
    # the nucleus/body, not just to the membrane edge). With it, the cell migrates
    # as a unit and the tent cohesion drags the followers. Calibrated so the
    # whole-cell net (n_node * migrate*f_act ~ 162 * 0.6 * 1.2e-10 ~ 1.2e-8 N) is
    # below the contact cap (5e-8 N) and per node 0.6*f_act=7.2e-11 N is in the
    # single-cell FA-clutch band => BAOAB-safe. Set 0 to recover basal-only splay.
    migrate_factor: float = 0.6
    lead_bias: float = 1.5          # leading-edge basal nodes pull extra (x this) — directional crawl

    # --- BULK-PRESSURE JUNCTION SWITCH (#6) -----------------------------------
    contact_factor: float = 2.6     # neighbour contact radius = factor*R (centroids)
    crowd_lo: float = 4.0
    crowd_hi: float = 12.0
    P_min_kPa: float = 0.0
    P_max_kPa: float = 6.0
    P_switch_kPa: float = 0.5
    cadherin_weak_factor: float = 0.3
    integrin_strong_factor: float = 3.0

    # --- LIVE PROLIFERATION ----------------------------------------------------
    # SLOW vs spreading (PI 2026-06-11): the cell-cycle timescale (~hours) is FAR
    # longer than the spreading/migration timescale (~minutes), so over a single
    # spread phase a rim cell divides at most ~0-2 times. p_div is lowered and the
    # driver's div cadence raised so the active TRACTION (centroid migration), not
    # division-stacking, is the visible A/A0 driver. The earlier p_div=0.18 +
    # div_every=4000 stacked ~6 division checks => 4-9 divisions => A/A0 was
    # division-dominated (~9), masking the ~2-3 traction band.
    p_div: float = 0.04             # per-eligible-rim-cell division prob / cadence (slow)
    div_gap_factor: float = 0.4     # daughter placed 2R + gap*R outward (no overlap)

    # --- necrosis 3-zone (depth-from-surface) ---------------------------------
    d_prolif_um: float = 40.0
    d_necrotic_um: float = 150.0

    z_substrate: float = 0.0
    kT: float = 4.28e-21
    extras: dict[str, Any] = field(default_factory=dict)

    def to_native(self, *, dt: float | None = None) -> ResolvedNativeDCM:
        """Project the mechanical base onto a ``ResolvedNativeDCM``."""
        return ResolvedNativeDCM(
            R_cell=self.R_cell, subdivisions=self.subdivisions, K_bulk=self.K_bulk,
            turgor_inflate=self.turgor_inflate, k_edge=self.k_edge,
            gamma_node=self.gamma_node, W_cs_Jm2=self.W_cs_Jm2,
            rep_strength=self.rep_strength, adh_strength=self.adh_strength,
            c_adh=self.c_adh, force_cap=self.contact_force_cap,
            spacing_factor=self.spacing_factor, cluster="3d",
            dt=dt if dt is not None else self.dt, seed=self.seed,
            z_substrate=self.z_substrate, kT=self.kT)


# ---------------------------------------------------------------------------
# (#5) ACTIVE RIM TRACTION — coarse-grained lamellipodium + belt + FA clutch
# ---------------------------------------------------------------------------
class ActiveRimTraction(md.force.Custom):
    """Active OUTWARD traction on basal rim-cell nodes (coarse-grained engine).

    Represents the integrated single-cell spreading machinery
    (``spreading_drive.py``) at the DCM CELL scale, as a per-node body force on the
    existing mesh nodes (NO topology mutation — see module docstring). Each step:

      1. Identify RIM cells: active cells whose centroid is on the cluster
         periphery (few close neighbours) AND that have at least one basal node in
         substrate contact (z < z0 + band*R). These are the cells with a free edge
         on the substrate — the only ones that can put down a lamellipodium.
      2. For each rim cell, its BASAL nodes (z < z0 + band*R) get an ACTIVE OUTWARD
         traction f_act * r_hat, where r_hat is the in-plane (xy) unit vector from
         the CLUSTER centroid through the cell centroid — the lamellipodial
         protrusion + actomyosin contraction-belt pull + FA-clutch grip pulling the
         rim outward. Per-cell the traction is scaled by ``int_mult`` (the junction
         switch's integrin gain) so a pressure-switched cell pulls harder.
      3. (optional) a weak INWARD apical line tension (belt_factor*f_act) on each
         rim cell's APICAL nodes (z high) drags the cell BODY after the protruding
         edge — the actomyosin contraction belt so the cell follows, not just the
         basal lip.

    The force is RAMPED (0 -> full over ramp_steps) and CAPPED per node (f_cap) so
    it never shocks the BAOAB integrator (the int32-blowup guard). It is an active
    energy input: potential_energy is reported 0 (a non-conservative body force).

    Args:
        cell_of_node: (N,) int32, mutable; cell id of each node (-1 dormant).
        ranges: list[(lo,hi)] per cell id, node-tag range (fixed pool).
        active: (n_max,) bool, mutable; which pool cells are active.
        int_mult: (n_max,) float, per-cell integrin gain (junction switch raises it).
        R_cell, z0: geometry [m].
        f_act, f_cap: base traction + per-node cap [N].
        ramp_steps: steps to ramp the traction in.
        contact_band, neighbour_factor, max_neighbours, integrin_switch_gain,
        belt_factor: rim-detection + belt knobs (see ResolvedActiveSpheroid).
    """

    def __init__(self, *, cell_of_node: np.ndarray, ranges, active: np.ndarray,
                 int_mult: np.ndarray, R_cell: float, z0: float, f_act: float,
                 f_cap: float, ramp_steps: int, contact_band: float,
                 neighbour_factor: float, max_neighbours: int,
                 integrin_switch_gain: float, belt_factor: float,
                 migrate_factor: float = 0.6, lead_bias: float = 1.5) -> None:
        super().__init__(aniso=False)
        self.cell_of_node = cell_of_node
        self.ranges = ranges
        self.active = active
        self.int_mult = int_mult
        self.R = float(R_cell)
        self.z0 = float(z0)
        self.f_act = float(f_act)
        self.f_cap = float(f_cap)
        self.ramp_steps = max(1, int(ramp_steps))
        self.contact_band = float(contact_band)
        self.r_neigh = float(neighbour_factor * R_cell)
        self.max_neigh = int(max_neighbours)
        self.switch_gain = float(integrin_switch_gain)
        self.belt = float(belt_factor)
        self.migrate = float(migrate_factor)
        self.lead_bias = float(lead_bias)
        # diagnostics (read by the driver)
        self.rim_cells: np.ndarray = np.empty(0, dtype=np.int64)
        self.f_per_cell: dict[int, float] = {}
        # per-rim-cell (centroid, outward-unit-vector, net-cell-force) for arrow viz
        self.cell_centroids: dict[int, np.ndarray] = {}
        self.cell_rhat: dict[int, np.ndarray] = {}
        self.cell_net_force: dict[int, np.ndarray] = {}
        self._ramp = 0.0

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        n = pos.shape[0]
        perm = np.argsort(tag)
        pos_g = pos[perm]
        F_g = np.zeros_like(pos_g)

        ramp = float(min(1.0, timestep / self.ramp_steps))
        self._ramp = ramp
        active_ids = np.where(self.active)[0]
        if active_ids.size < 2 or ramp <= 0.0:
            F = np.empty_like(pos)
            F[perm] = F_g
            with self.cpu_local_force_arrays as arr:
                arr.force[:] = F
                arr.potential_energy[:] = np.zeros(n, dtype=np.float64)
            return

        # per-cell centroid
        cents = np.array([pos_g[self.ranges[int(c)][0]:self.ranges[int(c)][1]].mean(0)
                          for c in active_ids])
        cluster_cen = cents.mean(0)
        # neighbour count by centroid distance (rim = few neighbours)
        d2 = np.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
        within = d2 < self.r_neigh ** 2
        np.fill_diagonal(within, False)
        crowd = within.sum(axis=1)

        rim_list = []
        f_per_cell = {}
        cell_centroids = {}
        cell_rhat = {}
        cell_net = {}
        zc = self.z0 + self.contact_band * self.R
        for k, c in enumerate(active_ids):
            c = int(c)
            if crowd[k] > self.max_neigh:
                continue  # interior cell (well-coordinated) — not a rim cell
            lo, hi = self.ranges[c]
            cell_pos = pos_g[lo:hi]
            basal = cell_pos[:, 2] < zc
            if not basal.any():
                continue  # not in substrate contact — cannot lamellipodiate
            rim_list.append(c)
            # in-plane outward direction (cluster centroid -> cell centroid)
            rxy = cents[k][:2] - cluster_cen[:2]
            rn = float(np.hypot(rxy[0], rxy[1]))
            if rn < 1e-12:
                continue  # cell at the very centre — no defined outward dir
            rhat = np.array([rxy[0] / rn, rxy[1] / rn, 0.0])
            gain = self.int_mult[c]
            fmag = ramp * self.f_act * gain
            fmag = min(fmag, self.f_cap)
            f_per_cell[c] = fmag
            cell_F = np.zeros_like(cell_pos)  # this cell's added force (for diag)
            idx = np.where(basal)[0]

            # LEGACY basal-only splay law (migrate_factor=0 AND lead_bias=1.0): the
            # original scheme — uniform outward fmag on every basal node + inward
            # belt on apical. Kept bit-exact so the GPU twin's parity gate (a
            # separate frozen-file port) still holds against this configuration.
            legacy = (self.migrate == 0.0 and self.lead_bias == 1.0)
            if legacy:
                cell_F[idx] += fmag * rhat
            else:
                # (1) WHOLE-CELL MIGRATION: a NET outward body force on EVERY node so
                # the centroid translocates as a unit (clutch transmits traction to
                # the cell body, not just the membrane lip). This is what makes the
                # cell MIGRATE rather than only splay its basal footprint.
                if self.migrate > 0.0:
                    fm = ramp * self.migrate * self.f_act * gain
                    cell_F += fm * rhat
                # (2) LAMELLIPODIAL leading edge: basal nodes on the OUTWARD
                # (leading) side pull extra — the protruding lamellipodium / clutch
                # grip; a directional crawl bias.
                off_xy = cell_pos[idx, :2] - cents[k][:2]
                proj = off_xy @ rhat[:2]
                lead = idx[proj > 0.0]
                trail = idx[proj <= 0.0]
                cell_F[lead] += self.lead_bias * fmag * rhat
                cell_F[trail] += 0.3 * fmag * rhat   # weak trailing-edge grip

            # (3) apical INWARD contraction belt — the actomyosin belt that drags
            # the cell body after the edge; kept SMALL so it does NOT cancel the
            # net migration (the old belt over-cancelled => zero centroid drift).
            if self.belt > 0.0:
                apic = np.where(~basal)[0]
                if apic.size:
                    fb = min(self.belt * fmag, self.f_cap)
                    cell_F[apic] += -fb * rhat

            F_g[lo:hi] += cell_F
            cell_centroids[c] = cents[k].copy()
            cell_rhat[c] = rhat.copy()
            cell_net[c] = cell_F.sum(axis=0)

        # cap per node (BAOAB int32 guard) — never let a node exceed f_cap
        fn = np.linalg.norm(F_g, axis=1)
        over = fn > self.f_cap
        if over.any():
            F_g[over] *= (self.f_cap / fn[over])[:, None]

        self.rim_cells = np.array(rim_list, dtype=np.int64)
        self.f_per_cell = f_per_cell
        self.cell_centroids = cell_centroids
        self.cell_rhat = cell_rhat
        self.cell_net_force = cell_net

        F = np.empty_like(pos)
        F[perm] = F_g
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F
            arr.potential_energy[:] = np.zeros(n, dtype=np.float64)


# ---------------------------------------------------------------------------
# (#6) Bulk-pressure + junction-switch updaters, wired to native+tent geometry
# ---------------------------------------------------------------------------
class SpheroidPressureProbe(hoomd.custom.Action):
    """Per-cell bulk-pressure proxy (neighbour crowding -> kPa), native geometry.

    Same crowding->kPa law as ``dcm_spheroid_state.PressureProbe`` but reads the
    native-shell ``ranges``/``active`` pool (mutable) so it sees daughters that the
    proliferation updater activates. Writes ``st.pressure_kPa`` (read by the
    junction switch).
    """

    def __init__(self, *, p: ResolvedActiveSpheroid, st: SpheroidStateArrays,
                 ranges, active: np.ndarray):
        super().__init__()
        self.p = p
        self.st = st
        self.ranges = ranges
        self.active = active
        self.r_contact = float(p.contact_factor * p.R_cell)
        self._sim = None

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        with self._sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        pos_g = pos[np.argsort(tag)]
        active_ids = np.where(self.active)[0]
        if active_ids.size == 0:
            return
        cents = np.array([pos_g[self.ranges[int(c)][0]:self.ranges[int(c)][1]].mean(0)
                          for c in active_ids])
        d2 = np.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
        within = d2 < self.r_contact ** 2
        np.fill_diagonal(within, False)
        crowd = within.sum(axis=1).astype(float)
        frac = np.clip((crowd - self.p.crowd_lo)
                       / (self.p.crowd_hi - self.p.crowd_lo), 0.0, 1.0)
        P = self.p.P_min_kPa + (self.p.P_max_kPa - self.p.P_min_kPa) * frac
        for k, c in enumerate(active_ids):
            self.st.pressure_kPa[int(c)] = P[k]


class SpheroidJunctionSwitch(hoomd.custom.Action):
    """Pressure-gated cadherin->integrin switch on the native+tent spheroid.

    For each active, non-necrotic cell whose bulk pressure exceeds ``P_switch``:
    set ``switched=True``, WEAKEN cadherin (``st.cad_mult`` down — read live by the
    ``DcmTentContact`` via its ``cad_mult`` seam) and STRENGTHEN integrin
    (``st.int_mult`` up — read live by ``ActiveRimTraction`` + the substrate gain).
    Latched within a run; necrotic core cells excluded. This is the pressure-release
    mechanism: loaded rim cells drop cell-cell cohesion and grip the ECM -> unjam.
    """

    def __init__(self, *, p: ResolvedActiveSpheroid, st: SpheroidStateArrays,
                 active: np.ndarray):
        super().__init__()
        self.p = p
        self.st = st
        self.active = active

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


class SpheroidNecrosis(hoomd.custom.Action):
    """3-zone depth-from-surface necrosis on the native+tent pool (live cells).

    depth d = R_cluster - r_cell (cluster-centroid radius); d < d_prolif ->
    PROLIFERATING (rim), d_prolif <= d < d_necrotic -> QUIESCENT, d >= d_necrotic
    -> NECROTIC (irreversible; gates division off). Writes ``st.state`` +
    ``st.depth_um``. R_cluster = max centroid radius + one R_cell (outer surface).
    """

    def __init__(self, *, p: ResolvedActiveSpheroid, st: SpheroidStateArrays,
                 ranges, active: np.ndarray):
        super().__init__()
        self.p = p
        self.st = st
        self.ranges = ranges
        self.active = active
        self._sim = None

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        with self._sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        pos_g = pos[np.argsort(tag)]
        active_ids = np.where(self.active)[0]
        if active_ids.size == 0:
            return
        cents = np.array([pos_g[self.ranges[int(c)][0]:self.ranges[int(c)][1]].mean(0)
                          for c in active_ids])
        cc = cents.mean(0)
        r = np.linalg.norm(cents - cc, axis=1)
        R_cluster = r.max() + self.p.R_cell
        depth_um = (R_cluster - r) * _UM
        for k, c in enumerate(active_ids):
            c = int(c)
            self.st.depth_um[c] = depth_um[k]
            if self.st.state[c] == int(CellState.NECROTIC):
                continue
            if depth_um[k] >= self.p.d_necrotic_um:
                self.st.state[c] = int(CellState.NECROTIC)
            elif depth_um[k] < self.p.d_prolif_um:
                self.st.state[c] = int(CellState.PROLIFERATING)
            else:
                self.st.state[c] = int(CellState.QUIESCENT)


# ---------------------------------------------------------------------------
# LIVE PROLIFERATION — pre-allocated native-shell pool (no mesh split)
# ---------------------------------------------------------------------------
class NativeRimProliferation(hoomd.custom.Action):
    """Rim-biased contact-inhibited division on the native-shell cell pool.

    No mesh split: each pool cell already owns nv mesh nodes + a ``cell{c}`` mesh
    triangle type with a live ``Volume`` entry. A division ACTIVATES a dormant pool
    cell as the daughter — copies an undeformed x-mirrored icosphere into the
    daughter's node range, placed force-free OUTWARD from a dividing rim cell
    (2R + gap, surface gap => no contact-core overlap => BAOAB-safe), flips
    ``active[daughter]=True`` and resets its per-cell state. The daughter's mesh
    ``Volume`` entry already exists (set at build for the whole pool) so turgor
    conserves it immediately — nothing in the fixed mesh-tag space changes.

    Gating: a cell divides only if ACTIVE, state==PROLIFERATING (rim by necrosis
    depth) AND a free-edge cell (3D convex-hull vertex of active centroids =
    contact-inhibition geometry). Necrotic/quiescent cells never divide.
    """

    def __init__(self, *, p: ResolvedActiveSpheroid, st: SpheroidStateArrays,
                 cell_of_node: np.ndarray, active: np.ndarray, ranges,
                 verts0: np.ndarray, nv: int, p_div: float, gap_factor: float,
                 rng_seed: int = 99):
        super().__init__()
        self.p = p
        self.st = st
        self.cell_of_node = cell_of_node
        self.active = active
        self.ranges = ranges
        self.verts0 = np.asarray(verts0, dtype=np.float64)
        self.nv = int(nv)
        self.R = float(p.R_cell)
        self.z0 = float(p.z_substrate)
        self.p_div = float(p_div)
        self.gap = float(gap_factor * p.R_cell)
        self._rng = np.random.default_rng(rng_seed)
        self._sim = None
        self.n_divisions = 0

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        if not np.any(~self.active):
            return  # pool exhausted
        with self._sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        perm = np.argsort(tag)
        pos_g = pos[perm]
        active_ids = np.where(self.active)[0]
        if active_ids.size < 4:
            return
        cents = np.array([pos_g[self.ranges[int(c)][0]:self.ranges[int(c)][1]].mean(0)
                          for c in active_ids])
        cluster_cen = cents.mean(0)
        rim = np.zeros(active_ids.shape[0], dtype=bool)
        try:
            from scipy.spatial import ConvexHull
            hull = ConvexHull(cents)
            rim[np.unique(hull.vertices)] = True
        except Exception:  # noqa: BLE001
            return
        divided = False
        for k in np.where(rim)[0]:
            c = int(active_ids[k])
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
            lo, hi = self.ranges[daughter]
            pos_g[lo:hi] = self.verts0 + new_center
            self.cell_of_node[lo:hi] = daughter
            self.active[daughter] = True
            self.st.state[daughter] = int(CellState.PROLIFERATING)
            self.st.cad_mult[daughter] = 1.0
            self.st.int_mult[daughter] = 1.0
            self.st.switched[daughter] = False
            self.n_divisions += 1
            divided = True
        if divided:
            pos_back = np.empty_like(pos)
            pos_back[perm] = pos_g
            with self._sim.state.cpu_local_snapshot as snap:
                np.asarray(snap.particles.position)[:] = pos_back


# ---------------------------------------------------------------------------
# Per-cell integrin-modulated substrate (junction switch strengthens grip)
# ---------------------------------------------------------------------------
class IntegrinModulatedSubstrate(md.force.Custom):
    """Adhesive substrate floor with a PER-CELL integrin gain (capped-harmonic).

    Same capped-harmonic well as ``dcm.DcmSubstrateForce`` but each node's well
    depth/stiffness is scaled by its cell's ``int_mult`` (1.0 baseline; >1 after
    the junction switch). So a pressure-switched cell GRIPS the substrate harder —
    the integrin half of the cadherin->integrin switch. Cap recomputed per node so
    the steric push-up stays finite at any gain (BAOAB-safe). Dormant nodes (cell
    -1) keep baseline gain (parked far from the floor — never in reach anyway).
    """

    def __init__(self, *, z0: float, W_cs: float, adh_range: float,
                 cell_of_node: np.ndarray, int_mult: np.ndarray) -> None:
        super().__init__(aniso=False)
        self.z0 = float(z0)
        self.W_cs = float(W_cs)
        self.rng = float(adh_range)
        self.k_well0 = 2.0 * self.W_cs / (self.rng ** 2)
        self.cell_of_node = cell_of_node
        self.int_mult = int_mult

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
        k = self.k_well0 * mult

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
# Full assembly — native+tent base + active traction + state updaters + pool
# ---------------------------------------------------------------------------
def build_active_spheroid(p: ResolvedActiveSpheroid, n_active: int, n_max: int,
                          *, device=None, belt: bool = True,
                          integrin_substrate: bool = True):
    """Assemble the ACTIVE native+tent spheroid with a pre-allocated pool.

    Builds an ``n_max``-cell native+tent spheroid (the validated stack) where the
    first ``n_active`` cells start active in a 3D ball and the rest are parked far
    off (dormant pool for division). Adds: ``ActiveRimTraction`` (#5), a per-cell
    ``cad_mult`` on the tent contact (#6 cadherin seam), and returns the
    pressure/junction/necrosis/proliferation updaters UNATTACHED so the driver
    schedules them. Optionally swaps the substrate for the integrin-modulated one.

    Returns a dict of handles (sim, traction, tent, st, active, ranges, updaters…).
    """
    # The native builder packs ALL n_max cells in a ball; we want n_active in the
    # ball + the rest parked. So we build n_max cells then PARK the surplus far off
    # and mark them dormant. (Simplest reuse of the validated builder + mesh.)
    base = p.to_native()
    h = build_native_dcm_simulation(base, n_max, substrate=True, contact=True,
                                    device=device)
    sim = h["sim"]
    ranges = h["ranges"]
    cell_of_tag = h["cell_of_tag"]      # (N,) cell id per node (immutable here)
    nv = h["nv"]
    tent = h["tent"]
    vol = h["volume"]
    V0 = h["V0"]

    # x-mirrored rest shell (the daughter template, identical to the native build)
    verts0, _, tris0 = native_icosphere(p.R_cell, p.subdivisions)

    # Pool bookkeeping: active[c], mutable cell_of_node (the tent/traction read it).
    active = np.zeros(n_max, dtype=bool)
    active[:n_active] = True
    cell_of_node = cell_of_tag.astype(np.int64).copy()  # node -> cell (mutable)

    st = SpheroidStateArrays(n_max)
    st.state[:n_active] = int(CellState.QUIESCENT)  # zone assigned on first necrosis tick

    # RE-PACK the active cells into a COMPACT n_active ball resting on the
    # substrate + PARK the surplus (dormant) cells far off.
    #
    # Why re-pack (the bug fix): the native builder packs ALL n_max cells into one
    # FCC ball sized for n_max, so the first n_active cells are scattered THROUGH a
    # tall n_max-ball and their lowest node floats ~12 µm above z0 (out of the
    # R-range substrate well) — the active subset never touches the floor and the
    # footprint measure collapses. We instead place the n_active active cells into
    # their OWN compact n_active FCC ball (via the validated `_cluster_centers`)
    # and LOWER it so the lowest node rests one node-gap above z0 (a small gap, not
    # a t=0 overlap → BAOAB-safe; the capped well then pulls it into contact).
    from ffn_sim.cell.dcm import _cluster_centers
    spacing = p.spacing_factor * p.R_cell
    centers = _cluster_centers(n_active, spacing, p.z_substrate, p.R_cell, mode="3d")

    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag).copy()
        pos = np.asarray(snap.particles.position).copy()
    perm = np.argsort(tag)
    pos_g = pos[perm]
    for c in range(n_active):
        lo, hi = ranges[c]
        pos_g[lo:hi] = verts0 + centers[c]
    park0 = 60.0 * p.R_cell
    for d, c in enumerate(range(n_active, n_max)):
        lo, hi = ranges[c]
        px = park0 + (d % 8) * (3.0 * p.R_cell)
        py = park0 + (d // 8) * (3.0 * p.R_cell)
        pos_g[lo:hi] = verts0 + np.array([px, py, park0])
        cell_of_node[lo:hi] = -1
    # lower the active cluster so its lowest node sits 0.3 R above the floor.
    act_nodes = np.concatenate([np.arange(ranges[c][0], ranges[c][1])
                                for c in range(n_active)])
    zmin = pos_g[act_nodes, 2].min()
    pos_g[act_nodes, 2] += (p.z_substrate + 0.3 * p.R_cell) - zmin
    pos_back = np.empty_like(pos)
    pos_back[perm] = pos_g
    with sim.state.cpu_local_snapshot as snap:
        np.asarray(snap.particles.position)[:] = pos_back

    # Wire the per-cell cadherin seam onto the tent contact + the mutable
    # cell_of_node (so dormant nodes don't interact, daughters do once activated).
    tent.cell_of_node = cell_of_node
    tent.cad_mult = st.cad_mult

    ig = sim.operations.integrator

    # Active rim traction (#5).
    traction = ActiveRimTraction(
        cell_of_node=cell_of_node, ranges=ranges, active=active,
        int_mult=st.int_mult, R_cell=p.R_cell, z0=p.z_substrate,
        f_act=p.f_act, f_cap=p.f_cap, ramp_steps=p.ramp_steps,
        contact_band=p.rim_contact_band, neighbour_factor=p.rim_neighbour_factor,
        max_neighbours=p.rim_max_neighbours,
        integrin_switch_gain=p.integrin_switch_gain, belt_factor=p.belt_factor,
        migrate_factor=p.migrate_factor, lead_bias=p.lead_bias)
    ig.forces.append(traction)

    # Optionally replace the substrate with the integrin-modulated one (#6 grip).
    if integrin_substrate:
        old_sub = h.get("substrate")
        if old_sub is not None and old_sub in ig.forces:
            ig.forces.remove(old_sub)
        area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
        W_cs = p.W_cs_Jm2 * area_per_node
        sub_mod = IntegrinModulatedSubstrate(
            z0=p.z_substrate, W_cs=W_cs, adh_range=p.R_cell,
            cell_of_node=cell_of_node, int_mult=st.int_mult)
        ig.forces.append(sub_mod)
        h["substrate"] = sub_mod

    # State updaters (UNATTACHED; the driver schedules them per phase).
    pressure = SpheroidPressureProbe(p=p, st=st, ranges=ranges, active=active)
    junction = SpheroidJunctionSwitch(p=p, st=st, active=active)
    necrosis = SpheroidNecrosis(p=p, st=st, ranges=ranges, active=active)
    prolif = NativeRimProliferation(
        p=p, st=st, cell_of_node=cell_of_node, active=active, ranges=ranges,
        verts0=verts0, nv=nv, p_div=p.p_div, gap_factor=p.div_gap_factor,
        rng_seed=p.seed + 3)

    return {
        "sim": sim, "traction": traction, "tent": tent, "substrate": h["substrate"],
        "volume": vol, "V0": V0, "st": st, "active": active,
        "cell_of_node": cell_of_node, "ranges": ranges, "nv": nv,
        "verts0": verts0, "tris0": tris0, "n_active": n_active, "n_max": n_max,
        "pressure": pressure, "junction": junction, "necrosis": necrosis,
        "prolif": prolif, "p": p,
    }
