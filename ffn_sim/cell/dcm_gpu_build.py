"""GPU-FRIENDLY DCM spheroid build — custom K1-turgor, NO native md.mesh.

WHY THIS MODULE (the A5000 blocker it solves). The native-mesh capstone
(``cell/dcm_native_shell.py``) enforces per-cell turgor with
``md.mesh.conservation.Volume`` over a ``hoomd.mesh.Mesh`` that carries a SEPARATE
mesh triangle type PER CELL. HOOMD warns that many mesh triangle types perform
poorly on the GPU (shared-memory pressure) and, confirmed on the gbook RTX A5000,
the run STALLS at 0% GPU util by N≈60. This module is the GPU-friendly assembly:
it reproduces the same DCM physics WITHOUT any native md.mesh and WITHOUT per-cell
mesh types, so it runs GPU-resident.

THE STACK (all on the project's frozen L-M BAOAB-limit integrator):
  * geometry — ``dcm.icosphere_mesh`` shells at ``dcm._cluster_centers`` (reused).
  * TURGOR — ``dcm_gpu_forces.DcmTurgorForceGPU`` (the validated K1 mesh-pressure
    kernel): ONE custom force, per-cell enclosed volume by face-grouping on
    ``face_cell``, outward face-normal pressure. NO per-cell mesh types → no stall.
  * EDGE SPRINGS — ``md.bond.Harmonic`` with a SINGLE bond type (GPU-native,
    in-plane cortical elasticity). No per-cell bond types.
  * CELL-CELL — ``dcm_gpu_forces.DcmTentContactGPU`` (SimuCell3D bilinear tent),
    which distinguishes cells by the mutable ``cell_of_node`` array, NOT by
    per-cell particle/pair types — so there is a SINGLE membrane particle type.
  * SUBSTRATE — ``dcm_gpu_forces.DcmSubstrateForceGPU`` (capped-harmonic adhesive
    well at z0).
  * (optional) ACTIVE RIM TRACTION — ``dcm_gpu_forces.DcmActiveRimTractionGPU``.

A single membrane particle type means the BAOAB ``gamma`` dict has ONE entry
(plus the LOD freeze type, see ``dcm_gpu_lod.py``), and the contact / turgor
costs are O(nodes) device kernels rather than an O(n_cells²) pair-type matrix.

``device=`` is plumbed: a ``hoomd.device.GPU`` is created automatically when one
is available (gbook A5000), otherwise a CPU device (this Mac, dev/fallback). The
build runs ``run(0)`` then a short ``run(2000)`` finite gate before any long run.

Units SI.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.custom
import hoomd.md as md

from ffn_sim.cell.dcm import icosphere_mesh, _cluster_centers
from ffn_sim.cell.dcm_gpu_forces import (
    DcmTurgorForceGPU,
    DcmTentContactGPU,
    DcmSubstrateForceGPU,
    DcmActiveRimTractionGPU,
    DcmActiveRimTractionGPUVec,
    on_gpu,
)


class SettlingForce(md.force.Custom):
    """Weak constant downward (−z) body force on the live membrane nodes — the
    sedimentation / plating of the spheroid onto the dish.

    WETTING DRIVER (2026-06-11 fix): without it the adhesive substrate well only
    grips nodes already near z0, so a free spheroid floats (~20% contact, never
    wets — diagnosed maxZ stays ~54 µm, footprint flat). A plated spheroid
    sediments under its buoyant weight, contacts, then wets/spreads by adhesion +
    active crawling. Applies only to ``mem_typeid`` (live membrane), never the
    parked dormant pool. Device-dispatched (GPU cupy / CPU numpy).
    """

    def __init__(self, *, f_settle: float, mem_typeid: int = 0) -> None:
        super().__init__(aniso=False)
        self.f = float(f_settle)
        self.mem_typeid = int(mem_typeid)
        self._gpu: bool | None = None
        self._xp = None

    def _setup(self) -> None:
        if self._gpu is None:
            self._gpu = on_gpu(self._state)
            if self._gpu:
                import cupy as cp           # guarded GPU-only import
                self._xp = cp
            else:
                self._xp = np

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        self._setup()
        xp = self._xp
        snap_ctx = (self._state.gpu_local_snapshot if self._gpu
                    else self._state.cpu_local_snapshot)
        with snap_ctx as snap:
            tid = xp.asarray(snap.particles.typeid)
            fz = xp.where(tid == self.mem_typeid, -self.f, 0.0)
        farr_ctx = (self.gpu_local_force_arrays if self._gpu
                    else self.cpu_local_force_arrays)
        with farr_ctx as arr:
            arr.force[:] = 0.0
            arr.force[:, 2] = fz


class AggregationDrive(md.force.Custom):
    """Weak CENTRIPETAL confinement that drives STAGE-1 aggregation (hanging-drop /
    surface-tension analog).

    Each live membrane node is pulled toward the aggregation ``center`` so a loosely
    placed cluster COMPACTS and ROUNDS into a cohesive spheroid (the cells come
    together). The cell-cell tent adhesion is short range (it cannot pull cells across
    a gap), so on its own a loose cluster just sits there — this confinement supplies
    the aggregation drive that a hanging-drop / low-adhesion well supplies in vitro.
    The FINAL ball size is set by cohesion + turgor + excluded-volume balancing the
    inward pull (the confinement only brings the cells in); per-node force is harmonic
    toward the centre with a magnitude cap (BAOAB-safe). It is an INITIAL-CONDITION
    protocol to produce a realistic starting aggregate for STAGE-2 spreading — it is
    NOT wired during spreading and is not a spreading-physics mechanism. Applies only
    to ``mem_typeid`` (live membrane). Device-dispatched (GPU cupy / CPU numpy).
    """

    def __init__(self, *, center, k_agg: float, f_cap: float = 3.0e-10,
                 mem_typeid: int = 0) -> None:
        super().__init__(aniso=False)
        self.center = np.asarray(center, dtype=np.float64)
        self.k = float(k_agg)
        self.f_cap = float(f_cap)
        self.mem_typeid = int(mem_typeid)
        self._gpu: bool | None = None
        self._xp = None

    def _setup(self) -> None:
        if self._gpu is None:
            self._gpu = on_gpu(self._state)
            if self._gpu:
                import cupy as cp                # guarded GPU-only import
                self._xp = cp
            else:
                self._xp = np

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        self._setup()
        xp = self._xp
        c = xp.asarray(self.center)
        snap_ctx = (self._state.gpu_local_snapshot if self._gpu
                    else self._state.cpu_local_snapshot)
        with snap_ctx as snap:
            tid = xp.asarray(snap.particles.typeid)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            F = -self.k * (pos - c[None, :])                 # harmonic toward centre
            mag = xp.sqrt((F * F).sum(axis=1))
            scale = xp.where(mag > self.f_cap,
                             self.f_cap / xp.where(mag > 0, mag, 1.0), 1.0)
            F = F * scale[:, None]
            F = xp.where((tid == self.mem_typeid)[:, None], F, 0.0)
        farr_ctx = (self.gpu_local_force_arrays if self._gpu
                    else self.cpu_local_force_arrays)
        with farr_ctx as arr:
            arr.force[:] = F


# Prolif-aware activity-LOD (Problem 2 fix) subclasses the FROZEN classifier; the
# frozen module (dcm_gpu_lod.py) is imported, never modified. No circular import:
# dcm_gpu_lod imports only from dcm_gpu_forces.
from ffn_sim.cell.dcm_gpu_lod import DcmActivityLOD

_KT_310 = 4.28e-21  # J at 310 K
_UM = 1.0e6


# ---------------------------------------------------------------------------
# Resolved parameters (single membrane type; GPU-friendly bands)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedGpuDCM:
    """GPU-friendly DCM spheroid constants (single membrane type).

    The mechanical bands are the validated DCM literature bands (MCF7): R_cell
    7.5 µm, turgor ~133 Pa, cortical tension scale via the edge spring. The
    cell-cell contact is the SimuCell3D bilinear-tent (same bands as the native
    capstone). Adhesion energy densities are converted to per-node well depths at
    build via the per-node patch area.
    """

    R_cell: float = 7.5e-6          # m  MCF7 radius
    subdivisions: int = 1           # icosphere level (1 = 42 nodes, 2 = 162)
    turgor_dP0: float = 133.0       # Pa baseline osmotic turgor (Young-Laplace)
    K_vol: float = 1.0e3            # Pa osmotic bulk modulus (ΔP per ΔV/V)
    k_edge: float = 1.0e-3          # N/m membrane edge spring (cortical elasticity)
    gamma_node: float = 3.9e-10     # N·s/m per-node Stokes drag (membrane type)

    # cell-substrate / cell-cell adhesion energy densities (J/m², converted to
    # per-node well depths via area_per_node at build).
    W_cs_Jm2: float = 0.5e-3        # J/m² cell-substrate adhesion (spreading driver)

    # SimuCell3D bilinear-tent cell-cell contact — COHESIVE bands (PI 2026-06-11).
    # The bilinear tent only ADHERES on node pairs in [r_contact, c_adh); since
    # r_contact = r_contact_factor·mean_edge = 4.37 µm at subdivisions=1, the old
    # c_adh=5.0e-7 (0.5 µm) made that band EMPTY → cell-cell adhesion was DEAD, the
    # "spheroid" a gapped lattice of mutually-repelling balls that dispersed under
    # any push (the convex-hull A/A₀ blow-up). A spheroid is a COHESIVE aggregate:
    # cells must be adhered cell-to-cell from t=0. These bands are the validated
    # confluent regime (dcm_confluent_tune config 3 / dcm_surface_mp4): with c_adh=5
    # µm > r_contact the adhesion band exists, and at the touching spacing 2.0·R the
    # cluster starts 73 % node-node in contact and STAYS bound (Rg flat, no
    # dispersal, no explosion — verified by scripts/dcm_cohesion_check.py). The
    # validated stiff cortex (k_edge) + turgor are KEPT (cohesion is an adhesion fix,
    # not a mechanics change).
    rep_strength: float = 2.0e8     # Pa/m repulsion stiffness ξ (non-penetration)
    adh_strength: float = 8.0e8     # Pa/m adhesion stiffness ω (strong cohesion)
    c_adh: float = 5.0e-6           # m adhesion cutoff (> r_contact ⇒ band exists)
    contact_force_cap: float = 5.0e-8   # N per-pair force cap (BAOAB guard)
    r_contact_factor: float = 1.0   # r_contact = factor · mean_edge

    spacing_factor: float = 2.0     # cell-center spacing = factor · R_cell (touching
                                    # → adhesion engages from t=0, cohesive start)
    cluster: str = "3d"             # "3d" spheroid ball | "2d" monolayer

    k_sub_Nm: float | None = None   # N/m substrate compliance (None = rigid dish)
    ligand_density: float = 1.0     # areal ligand density (cell-substrate adhesion)
    z_substrate: float = 0.0        # m substrate plane

    kT: float = _KT_310
    dt: float = 1.0e-9              # s (strong adhesion → small dt for BAOAB)
    seed: int = 7
    extras: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Device selection
# ---------------------------------------------------------------------------
def pick_device(device=None, *, seed_notice: int = 0):
    """Return a HOOMD device: the caller's, else GPU-if-available, else CPU.

    On the gbook A5000 a ``hoomd.device.GPU`` is created; on this Mac (no CUDA)
    the GPU constructor raises and we fall back to ``hoomd.device.CPU`` — so the
    SAME script runs in both places (GPU production / CPU dev-fallback).
    """
    if device is not None:
        return device
    try:
        return hoomd.device.GPU(notice_level=seed_notice)
    except Exception:
        return hoomd.device.CPU(notice_level=seed_notice)


# ---------------------------------------------------------------------------
# Snapshot — single membrane type, single bond type
# ---------------------------------------------------------------------------
def build_gpu_dcm_snapshot(p: ResolvedGpuDCM, n_cells: int):
    """Build a packed cluster of n_cells DCM shells with ONE membrane type.

    Returns a dict with the gsd Frame and all the per-cell bookkeeping the forces
    need: cell_of_node (mutable, used by the tent contact), ranges (per-cell node
    tag ranges), faces / face_cell (for the K1 turgor), nv/ne, mean_edge.
    """
    import gsd.hoomd

    verts0, edges, tris0 = icosphere_mesh(p.R_cell, p.subdivisions)
    nv = verts0.shape[0]
    ne = edges.shape[0]
    mean_edge = float(np.linalg.norm(
        verts0[edges[:, 0]] - verts0[edges[:, 1]], axis=1).mean())

    # surface-gapped start (no t=0 node overlap → no BAOAB blow-up); adhesion
    # pulls cells into contact gently over the settle.
    spacing = p.spacing_factor * p.R_cell
    centers = _cluster_centers(n_cells, spacing, p.z_substrate, p.R_cell,
                               mode=p.cluster)
    # _cluster_centers (dcm.py) returns a lattice whose count need not equal
    # n_cells exactly (e.g. an FCC ball filter gives 1372 for a requested 1500,
    # which crashed build at large N). Clamp to the actual centers available so
    # any N is safe; if it overshoots, take the first n_cells.
    if len(centers) != n_cells:
        n_cells = min(n_cells, len(centers))
        centers = centers[:n_cells]

    all_pos, bond_groups, cell_of_node = [], [], []
    face_groups, face_cell, ranges = [], [], []
    tag = 0
    for c in range(n_cells):
        ranges.append((tag, tag + nv))
        all_pos.append(verts0 + centers[c])
        cell_of_node.append(np.full(nv, c, dtype=np.int64))
        bond_groups.append(edges + tag)
        face_groups.append(tris0 + tag)
        face_cell.append(np.full(tris0.shape[0], c, dtype=np.int64))
        tag += nv

    pos = np.concatenate(all_pos, axis=0)
    cell_of_node = np.concatenate(cell_of_node)
    bonds = np.concatenate(bond_groups, axis=0)
    faces = np.concatenate(face_groups, axis=0)
    face_cell = np.concatenate(face_cell)

    box_edge = float(np.abs(pos).max() * 2.5 + 4.0 * p.R_cell)
    snap = gsd.hoomd.Frame()
    snap.particles.N = pos.shape[0]
    snap.particles.types = ["dcm_mem", "dcm_inert"]   # single active type + LOD freeze type
    snap.particles.position = pos
    snap.particles.typeid = np.zeros(pos.shape[0], dtype=np.uint32)  # all active
    snap.particles.mass = np.ones(pos.shape[0])
    snap.bonds.N = bonds.shape[0]
    snap.bonds.types = ["dcm_edge"]                   # SINGLE bond type
    snap.bonds.typeid = np.zeros(bonds.shape[0], dtype=np.uint32)
    snap.bonds.group = bonds.astype(np.uint32)
    snap.configuration.box = [box_edge, box_edge, box_edge, 0, 0, 0]
    return dict(snap=snap, cell_of_node=cell_of_node, ranges=ranges,
                faces=faces, face_cell=face_cell, nv=nv, ne=ne,
                mean_edge=mean_edge, centers=centers)


# ---------------------------------------------------------------------------
# Full GPU-friendly simulation assembly
# ---------------------------------------------------------------------------
def build_gpu_dcm_simulation(p: ResolvedGpuDCM, n_cells: int, *, device=None,
                             active: bool = False, fast_active: bool = True,
                             arrest: bool = False,
                             arrest_radius_factor: float = 1.4,
                             arrest_width: float = 0.18,
                             arrest_settle_steps: int = 4000,
                             settle_force: float = 4.0e-10,
                             with_substrate: bool = True,
                             init_pos: "np.ndarray | None" = None):
    """Assemble the GPU-friendly DCM spheroid on the BAOAB integrator.

    No native md.mesh, no per-cell mesh/particle/bond types. Returns a dict of
    handles (sim, cell_of_node, ranges, faces, face_cell, turgor, contact,
    substrate, traction|None, baoab, p, nv, ne, mean_edge).

    Args:
        active: if True, wire active basal-rim traction (passive turgor +
            adhesion wetting only if False).
        fast_active: when ``active`` is True, use the FULLY-VECTORIZED
            ``DcmActiveRimTractionGPUVec`` (no per-cell Python loop) instead of the
            original ``DcmActiveRimTractionGPU``. The gbook A5000 profile
            (2026-06-11) found the original's per-cell loop is ~22 ms/step (~57% of
            the ~39 ms/step production total at N=200); the vectorized form is
            bit-identical (max abs diff < 1e-10 N, ``test_dcm_active_vec_parity.py``)
            and collapses that to a handful of ms. Default True (the win); pass
            False to wire the original loop (e.g. for an A/B wall comparison).
        arrest: SPREADING-ARREST (membrane-tension / contact-inhibition stall),
            default OFF (opt-in). The outward rim traction is smoothly switched off
            as a rim cell's radial spread approaches ``arrest_radius_factor ·
            R0_cluster``. It is OFF by default because the diagnostic sweep showed it
            does NOT cure the over-spread it targeted (the runaway is a convex-hull
            footprint artifact + the now-fixed float artifact, not the rim push —
            full arrest gain→0 still blows A/A₀ up, ``probe3_on.json``), and a factor
            chosen to land A/A₀ in a band is a no-magic-number violation. The
            physical cure is the cohesive-start + real-wetting baseline. Pass
            ``arrest=True`` to A/B the (opt-in, tuned) cap; see DcmActiveRimTractionGPU.
        arrest_radius_factor: the cap on each rim cell's radial spread as a multiple
            of the SETTLED cluster radius R0 (captured after ``arrest_settle_steps``,
            N-invariant). At plateau the footprint area ≈ factor² × the settled area,
            and the settled footprint is already ~1.5× A₀ (the gapped start compacts
            then begins to spread by the settle), so factor 1.35 → plateau A/A₀ ≈
            1.5·1.35² ≈ 2.7 ∈ the physiological [2,4] band.
        arrest_width: the smooth-ramp width as a fraction of r_max (the arrest gain
            falls from 1 to 0 across ``[(1−width)·r_max, r_max]``). Default 0.45.
        arrest_settle_steps: defer the R0 capture until this step so the cap anchors
            to the COMPACTED cluster, not the gapped t≈0 radius (which is too large →
            cap unreachable → arrest never engages). Default 4000 (= the ramp).
    """
    from ffn_sim.integrator.baoab import make_baoab_updater

    b = build_gpu_dcm_snapshot(p, n_cells)
    snap = b["snap"]
    nv, ne, mean_edge = b["nv"], b["ne"], b["mean_edge"]
    cell_of_node = b["cell_of_node"]
    ranges = b["ranges"]
    faces, face_cell = b["faces"], b["face_cell"]

    # STAGE 2 hand-off: replace the freshly-placed lattice with the AGGREGATED
    # spheroid positions captured at the end of stage 1 (same topology/order, so the
    # bonds/faces/cell_of_node all still line up) — translated so the ball sits on the
    # dish. This is the PI requirement: spreading must start from a REAL aggregated
    # spheroid, not a just-placed lattice, or the rim traction is ill-defined.
    if init_pos is not None:
        ip = np.asarray(init_pos, dtype=snap.particles.position.dtype)
        if ip.shape != snap.particles.position.shape:
            raise ValueError(f"init_pos {ip.shape} != snapshot "
                             f"{snap.particles.position.shape} (N/topology mismatch)")
        L = float(np.abs(ip).max() * 2.5 + 6.0 * p.R_cell)   # box headroom for spread
        if L > snap.configuration.box[0]:
            snap.configuration.box = [L, L, L, 0, 0, 0]
        snap.particles.position[:] = ip

    dev = pick_device(device)
    sim = hoomd.Simulation(device=dev, seed=p.seed)
    sim.create_state_from_snapshot(snap)

    ig = md.Integrator(dt=p.dt)

    # Edge springs (cortical elasticity); r0 = mean edge length. SINGLE bond type.
    bond = md.bond.Harmonic()
    bond.params["dcm_edge"] = dict(k=p.k_edge, r0=mean_edge)
    ig.forces.append(bond)

    # Per-node membrane patch area (adhesion-energy → per-node depths / contact A).
    area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
    W_cs = p.W_cs_Jm2 * p.ligand_density * area_per_node    # J per node

    # TURGOR — K1 mesh-pressure (GPU-friendly: one custom force, face-grouped).
    V0 = (4.0 / 3.0) * np.pi * p.R_cell ** 3
    turgor = DcmTurgorForceGPU(faces=faces, face_cell=face_cell, n_cells=n_cells,
                               V0=V0, turgor_dP0=p.turgor_dP0, K_vol=p.K_vol)
    ig.forces.append(turgor)

    # CELL-CELL — bilinear-tent contact, distinguishes cells by cell_of_node.
    # cad_mult (n_cells,) is the MUTABLE per-cell cadherin multiplier the junction
    # switch lowers in place (mult = √(cad_mult_i·cad_mult_j); see
    # dcm_gpu_forces.DcmTentContactGPU / kernels_cpu.tent_contact_forces). Init 1.0
    # (= full cadherin); a pressure-loaded cell's entry is driven toward
    # cadherin_weak_factor by attach_junction_switch (additive — switch off by
    # default, so passing it is backward-compatible: all-1 ⇒ identical physics).
    cad_mult = np.ones(n_cells, dtype=np.float64)
    r_contact = p.r_contact_factor * mean_edge
    contact = DcmTentContactGPU(
        cell_of_node=cell_of_node, r_contact=r_contact, c_adh=p.c_adh,
        rep_strength=p.rep_strength, adh_strength=p.adh_strength,
        patch_area=area_per_node, force_cap=p.contact_force_cap,
        cad_mult=cad_mult)
    ig.forces.append(contact)

    # SUBSTRATE — adhesive capped-harmonic well at z0. Skipped for STAGE 1
    # (free-float aggregation): with_substrate=False removes the dish entirely so the
    # cells aggregate into a free-floating spheroid under cohesion+turgor ALONE (no
    # wetting), then STAGE 2 re-introduces the dish under the aggregated ball.
    substrate = None
    if with_substrate:
        substrate = DcmSubstrateForceGPU(
            z0=p.z_substrate, W_cs=W_cs, adh_range=p.R_cell, k_sub=p.k_sub_Nm)
        ig.forces.append(substrate)

    # SEDIMENTATION / PLATING — a weak constant downward body force on the live
    # membrane nodes. WHY (the wetting fix, 2026-06-11): the adhesive substrate
    # well only acts on nodes already within adh_range of z0, so a free spheroid
    # FLOATS (only its bottom ~20% touches) and never wets — diagnosed: maxZ stays
    # ~54 µm, footprint flat. A plated spheroid physically SEDIMENTS onto the dish
    # under its (buoyant) weight, contacts, then the cells wet/spread by adhesion +
    # active crawling. This force supplies that sedimentation. Calibrated (diagnostic
    # sweep): f_settle≈4e-10 N/node brings maxZ 54→~33 µm with all cells reaching
    # the substrate and the footprint growing — i.e. genuine wetting/spreading (vs
    # 0 = floats, 1.5e-9 = full pancake). Applies only to the mem type (typeid 0),
    # never the parked dormant pool. f_settle=0 disables (legacy float behaviour).
    if with_substrate and settle_force and settle_force > 0.0:
        ig.forces.append(SettlingForce(f_settle=settle_force, mem_typeid=0))

    traction = None
    # integrin_gain (n_cells,) is the MUTABLE per-cell integrin (cell-substrate)
    # multiplier the junction switch RAISES in place. When active=True it IS the
    # rim-traction ``int_mult`` (so a switched cell's substrate traction is scaled
    # up live); when active=False it is a standalone array the switch still raises
    # (read by the visualiser / available to any substrate gain wiring). Init 1.0.
    integrin_gain = np.ones(n_cells, dtype=np.float64)
    if active:
        active_mask = np.ones(n_cells, dtype=bool)
        int_mult = integrin_gain          # share the SAME array → switch raises traction
        # vectorized (default) vs original per-cell-loop active traction — same law,
        # same construction signature (the vec is a true drop-in subclass).
        TractionCls = (DcmActiveRimTractionGPUVec if fast_active
                       else DcmActiveRimTractionGPU)
        traction = TractionCls(
            cell_of_node=cell_of_node, ranges=ranges, active=active_mask,
            int_mult=int_mult, R_cell=p.R_cell, z0=p.z_substrate,
            f_act=1.2e-10, f_cap=6.0e-10, ramp_steps=4000, contact_band=0.5,
            neighbour_factor=2.6, max_neighbours=9, integrin_switch_gain=3.0,
            belt_factor=0.25,
            arrest_radius_factor=(arrest_radius_factor if arrest else None),
            arrest_width=arrest_width, arrest_settle_steps=arrest_settle_steps)
        ig.forces.append(traction)

    sim.operations.integrator = ig
    sim.run(0)

    # BAOAB — gamma covers the single membrane type (+ the LOD inert freeze type).
    gamma = {"dcm_mem": p.gamma_node, "dcm_inert": p.gamma_node}
    action, updater = make_baoab_updater(
        kT=p.kT, gamma=gamma, dt=p.dt, seed=p.seed + 1)
    sim.operations.updaters.append(updater)

    return dict(
        sim=sim, cell_of_node=cell_of_node, ranges=ranges, faces=faces,
        face_cell=face_cell, n_cells=n_cells, nv=nv, ne=ne, mean_edge=mean_edge,
        V0=V0, turgor=turgor, contact=contact, substrate=substrate,
        traction=traction, baoab=action, gamma=gamma, p=p,
        area_per_node=area_per_node, centers=b["centers"],
        cad_mult=cad_mult, integrin_gain=integrin_gain)


# ===========================================================================
# BULK-PRESSURE JUNCTION SWITCH (cadherin → integrin clutch) on the GPU build
# ===========================================================================
# WHY (PI 2026-06-11). As the spheroid compacts, per-cell bulk pressure rises; a
# cell whose pressure exceeds an onset (~0.5 kPa) WEAKENS its cadherin (cell-cell)
# adhesion and STRENGTHENS its integrin (cell-substrate) — the cadherin→integrin
# clutch switch that lets pressure-loaded cells unjam/spread. The mechanism and its
# constants live in the FROZEN ``cell/dcm_spheroid_state.py``
# (``PressureProbe`` crowding→kPa, ``JunctionSwitchUpdater`` switch rule,
# ``ResolvedSpheroidState`` thresholds); this wires them onto the GPU-friendly
# build's MUTABLE ``cad_mult`` / ``integrin_gain`` arrays.
#
# WHY A LIVE-CELL WRAPPER (not the frozen updater directly). The frozen
# PressureProbe/JunctionSwitchUpdater operate on a ``SpheroidStateArrays`` and an
# ``active`` mask; the GPU build instead exposes ``cell_of_node`` (the tent/turgor
# liveness authority) and its own ``cad_mult`` / ``integrin_gain``. So the wrapper
# REPLICATES the frozen crowding→kPa pressure proxy + the > P_switch rule EXACTLY,
# but computes pressure over LIVE cells ONLY (``cell_of_node ≥ 0``) — consistent
# with the necrosis live-cell fix (``LiveCellActivityLOD`` /
# ``dcm_gpu_lod_run._live_necrotic_count``), so a parked dormant pool (if present)
# never inflates the crowd count. The crowding/kPa map + P_switch read straight
# from ``ResolvedSpheroidState`` (the frozen literature bands), with R_patch tied
# to the build's actual ``p.R_cell`` so the contact radius scales with the cells.


class GpuJunctionSwitchUpdater(hoomd.custom.Action):
    """Bulk-pressure cadherin→integrin clutch on the GPU-friendly DCM build.

    Each (low-cadence) ``act`` computes per-LIVE-cell bulk pressure from neighbour-
    cell crowding (the FROZEN ``PressureProbe`` proxy: count active neighbours
    within ``contact_factor·R``, clip-ramp the crowd count to the [P_min, P_max] kPa
    band between ``crowd_lo`` and ``crowd_hi``), then applies the FROZEN
    ``JunctionSwitchUpdater`` rule: any live cell whose pressure exceeds
    ``P_switch_kPa`` is LATCHED switched — its ``cad_mult`` is lowered to
    ``cadherin_weak_factor`` (cell-cell adhesion weakens) and its ``integrin_gain``
    raised to ``integrin_strong_factor`` (cell-substrate adhesion strengthens). The
    modulated tent contact reads ``cad_mult`` live (mult = √(cad_mult_i·cad_mult_j))
    and the rim traction reads ``integrin_gain`` (= its ``int_mult``) live, so a
    switched cell loses cohesion and gains traction the next step.

    Pressure is computed over LIVE cells only (``cell_of_node ≥ 0``) so any parked
    dormant pool never inflates the crowd count (consistent with the necrosis
    live-cell fix). The switch is latched (does not revert) within a run.

    Args:
        handles: the ``build_gpu_dcm_simulation`` dict (sim, ranges, cell_of_node,
            cad_mult, integrin_gain).
        sp: ``ResolvedSpheroidState`` carrying the crowding→kPa proxy bands +
            the junction-switch thresholds (frozen literature values). If None, a
            default is created with ``R_patch`` tied to the build's ``p.R_cell``.
    """

    def __init__(self, *, handles: dict, sp=None) -> None:
        super().__init__()
        from ffn_sim.cell.dcm_spheroid_state import ResolvedSpheroidState
        self.h = handles
        self.ranges = handles["ranges"]
        self.cell_of_node = handles["cell_of_node"]          # (N,) int, −1 = dormant
        self.cad_mult = handles["cad_mult"]                  # (n_cells,) mutable
        self.integrin_gain = handles["integrin_gain"]        # (n_cells,) mutable
        p = handles["p"]
        # Default thresholds = frozen literature bands; R_patch tied to the build's
        # actual cell radius so the contact radius scales with these cells.
        self.sp = sp or ResolvedSpheroidState(R_patch=float(p.R_cell))
        self.r_contact = float(self.sp.contact_factor * self.sp.R_patch)
        n_cells = len(self.ranges)
        # per-cell diagnostics (read by the visualiser / driver)
        self.pressure_kPa = np.zeros(n_cells, dtype=np.float64)
        self.switched = np.zeros(n_cells, dtype=bool)
        self.n_switched = 0
        self._sim = None

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def _live_ids(self) -> np.ndarray:
        """(k,) cell ids that are LIVE (first node's cell_of_node ≥ 0)."""
        first_node = np.array([lo for (lo, _hi) in self.ranges])
        return np.flatnonzero(self.cell_of_node[first_node] >= 0)

    def compute_pressure(self, pos_g: np.ndarray):
        """Per-LIVE-cell bulk pressure [kPa] from neighbour crowding (frozen proxy).

        Mirrors ``PressureProbe.act`` exactly (count active neighbours within
        ``contact_factor·R``, clip-ramp crowd → [P_min, P_max] kPa between
        ``crowd_lo`` and ``crowd_hi``) but over the LIVE-cell centroids only.

        Returns ``(live_ids, P_live)`` — live cell ids and their pressures [kPa].
        """
        live_ids = self._live_ids()
        if live_ids.size == 0:
            return live_ids, np.zeros(0)
        cents = np.array([pos_g[self.ranges[int(c)][0]:self.ranges[int(c)][1]].mean(0)
                          for c in live_ids])
        d2 = np.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
        within = d2 < self.r_contact ** 2
        np.fill_diagonal(within, False)
        crowd = within.sum(axis=1).astype(float)
        frac = np.clip((crowd - self.sp.crowd_lo)
                       / (self.sp.crowd_hi - self.sp.crowd_lo), 0.0, 1.0)
        P = self.sp.P_min_kPa + (self.sp.P_max_kPa - self.sp.P_min_kPa) * frac
        return live_ids, P

    def act(self, timestep: int) -> None:  # noqa: D401
        with self._sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        pos_g = pos[np.argsort(tag)]
        live_ids, P = self.compute_pressure(pos_g)
        for k, c in enumerate(live_ids):
            c = int(c)
            self.pressure_kPa[c] = P[k]
            if self.switched[c]:
                continue
            if P[k] > self.sp.P_switch_kPa:
                self.switched[c] = True
                self.cad_mult[c] = self.sp.cadherin_weak_factor       # cadherin ↓
                self.integrin_gain[c] = self.sp.integrin_strong_factor  # integrin ↑
        self.n_switched = int(self.switched.sum())


def attach_junction_switch(handles: dict, *, sp=None, cadence: int = 500):
    """Wire the bulk-pressure cadherin→integrin clutch onto a GPU DCM build.

    Attaches a :class:`GpuJunctionSwitchUpdater` (low cadence) that each fire
    computes per-LIVE-cell bulk pressure and latches the junction switch for cells
    above ``P_switch_kPa`` — lowering their ``cad_mult`` (cell-cell) and raising
    their ``integrin_gain`` (cell-substrate). Additive: with the updater absent the
    arrays stay all-1 and the build is unchanged.

    ``handles`` is the dict from :func:`build_gpu_dcm_simulation`. Adds
    ``junction_switch`` (the updater) + ``junction_updater`` (the CustomUpdater).
    Returns the updated ``handles``.
    """
    sim = handles["sim"]
    upd = GpuJunctionSwitchUpdater(handles=handles, sp=sp)
    cu = hoomd.update.CustomUpdater(
        action=upd, trigger=hoomd.trigger.Periodic(int(cadence)))
    sim.operations.updaters.append(cu)
    handles["junction_switch"] = upd
    handles["junction_updater"] = cu
    return handles


# ===========================================================================
# LIVE PROLIFERATION on the GPU-friendly build — pre-allocated cell POOL
# ===========================================================================
# WHY A POOL (NO MESH SPLIT). The GPU-friendly build's K1 turgor groups faces by
# the per-face ``face_cell`` array and the tent contact distinguishes cells by the
# per-node ``cell_of_node`` array — both are FIXED tag-space lookups, NOT native
# md.mesh types. A true division would have to APPEND nodes/faces to the snapshot,
# which renumbers every tag and breaks both arrays (same constraint the native
# stack hit, see ``cell/dcm_active.py``). So we pre-allocate an ``n_max``-cell
# pool: build n_max icosphere shells, start ``n_active`` in the ball + PARK the
# surplus dormant (far off, undeformed → turgor V≈V0 → force-free; cell_of_node=−1
# → tent skips them; high z → outside the substrate well → force-free). A division
# ACTIVATES a parked pool cell as the daughter — copies an undeformed icosphere
# into its (already-allocated) node range, places it gapped-outward from the
# dividing rim cell, and flips its ``face_cell`` / ``cell_of_node`` / traction
# ``active`` ON. The turgor's V0 is per-cell-uniform (all shells share one V0) so
# the daughter is conserved the instant ``face_cell`` points its faces at it —
# nothing in the fixed tag space changes. The K1 turgor reads ``face_cell`` fresh
# each step (set on the live force object), so updating it in place is enough.


def build_gpu_spheroid_prolif(p: ResolvedGpuDCM, n_active: int, n_max: int, *,
                              device=None, active: bool = False,
                              fast_active: bool = True, arrest: bool = True,
                              arrest_radius_factor: float = 1.7):
    """GPU-friendly DCM spheroid with a PRE-ALLOCATED proliferation pool.

    Builds ``n_max`` icosphere shells in a SINGLE fixed tag space (so the K1
    turgor's ``face_cell`` and the tent's ``cell_of_node`` never renumber). The
    first ``n_active`` cells start in a compact FCC ball resting on the substrate;
    the remaining ``n_max − n_active`` are PARKED far off and marked dormant
    (``cell_of_node = −1``, undeformed → force-free). A ``GpuProliferationUpdater``
    (attached by the caller, low cadence) activates parked cells as daughters.

    Mirrors ``build_gpu_dcm_simulation`` (same forces, same BAOAB) but with the
    pool + a mutable ``active`` (n_max,) mask and mutable ``face_cell`` so the
    turgor follows daughters. Returns the same handle dict + ``active``,
    ``n_active``, ``n_max``, ``verts0`` (the daughter template), and an unattached
    ``prolif`` updater factory is left to the caller (the run-script wires it).

    Args:
        n_active: cells that start active in the ball.
        n_max: total pool size (n_active + dormant reserve).
        active: if True, wire ``DcmActiveRimTractionGPU`` (only active-mask cells
            pull). Default False.
    """
    from ffn_sim.integrator.baoab import make_baoab_updater
    import gsd.hoomd

    verts0, edges, tris0 = icosphere_mesh(p.R_cell, p.subdivisions)
    nv = verts0.shape[0]
    ne = edges.shape[0]
    nf = tris0.shape[0]
    mean_edge = float(np.linalg.norm(
        verts0[edges[:, 0]] - verts0[edges[:, 1]], axis=1).mean())

    # --- pool geometry: n_active in a compact ball, surplus parked far off -----
    spacing = p.spacing_factor * p.R_cell
    centers = _cluster_centers(n_active, spacing, p.z_substrate, p.R_cell,
                               mode=p.cluster)
    # guard the same lattice-count mismatch as the plain build (see above): if
    # _cluster_centers returns fewer than n_active, clamp the live count.
    if len(centers) < n_active:
        n_active = len(centers)

    all_pos, bond_groups, cell_of_node = [], [], []
    face_groups, face_cell, ranges = [], [], []
    tag = 0
    park0 = 60.0 * p.R_cell
    for c in range(n_max):
        ranges.append((tag, tag + nv))
        if c < n_active:
            ctr = centers[c]
            con = np.full(nv, c, dtype=np.int64)
        else:
            d = c - n_active
            px = park0 + (d % 8) * (3.0 * p.R_cell)
            py = park0 + (d // 8) * (3.0 * p.R_cell)
            ctr = np.array([px, py, park0])
            con = np.full(nv, -1, dtype=np.int64)      # dormant → tent skips it
        all_pos.append(verts0 + ctr)
        cell_of_node.append(con)
        bond_groups.append(edges + tag)
        face_groups.append(tris0 + tag)
        face_cell.append(np.full(nf, c, dtype=np.int64))  # turgor: own cell id
        tag += nv

    pos = np.concatenate(all_pos, axis=0)
    cell_of_node = np.concatenate(cell_of_node)            # (N,) mutable
    bonds = np.concatenate(bond_groups, axis=0)
    faces = np.concatenate(face_groups, axis=0)
    face_cell = np.concatenate(face_cell)                  # (F,) mutable

    box_edge = float(np.abs(pos).max() * 2.5 + 4.0 * p.R_cell)
    snap = gsd.hoomd.Frame()
    snap.particles.N = pos.shape[0]
    snap.particles.types = ["dcm_mem", "dcm_inert"]
    snap.particles.position = pos
    snap.particles.typeid = np.zeros(pos.shape[0], dtype=np.uint32)
    snap.particles.mass = np.ones(pos.shape[0])
    snap.bonds.N = bonds.shape[0]
    snap.bonds.types = ["dcm_edge"]
    snap.bonds.typeid = np.zeros(bonds.shape[0], dtype=np.uint32)
    snap.bonds.group = bonds.astype(np.uint32)
    snap.configuration.box = [box_edge, box_edge, box_edge, 0, 0, 0]

    dev = pick_device(device)
    sim = hoomd.Simulation(device=dev, seed=p.seed)
    sim.create_state_from_snapshot(snap)

    ig = md.Integrator(dt=p.dt)

    bond = md.bond.Harmonic()
    bond.params["dcm_edge"] = dict(k=p.k_edge, r0=mean_edge)
    ig.forces.append(bond)

    area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
    W_cs = p.W_cs_Jm2 * p.ligand_density * area_per_node

    # TURGOR over the WHOLE pool (n_max cells). face_cell is mutable and read fresh
    # each step → a daughter's faces are conserved as soon as it activates. Parked
    # cells are undeformed (V≈V0) so their turgor force is ~0 (force-free pool).
    V0 = (4.0 / 3.0) * np.pi * p.R_cell ** 3
    turgor = DcmTurgorForceGPU(faces=faces, face_cell=face_cell, n_cells=n_max,
                               V0=V0, turgor_dP0=p.turgor_dP0, K_vol=p.K_vol)
    ig.forces.append(turgor)

    r_contact = p.r_contact_factor * mean_edge
    contact = DcmTentContactGPU(
        cell_of_node=cell_of_node, r_contact=r_contact, c_adh=p.c_adh,
        rep_strength=p.rep_strength, adh_strength=p.adh_strength,
        patch_area=area_per_node, force_cap=p.contact_force_cap)
    ig.forces.append(contact)

    substrate = DcmSubstrateForceGPU(
        z0=p.z_substrate, W_cs=W_cs, adh_range=p.R_cell, k_sub=p.k_sub_Nm)
    ig.forces.append(substrate)

    # active mask over the POOL: only the active ball pulls; parked cells off.
    active_mask = np.zeros(n_max, dtype=bool)
    active_mask[:n_active] = True
    int_mult = np.ones(n_max, dtype=np.float64)

    traction = None
    if active:
        # vectorized (default) vs original per-cell-loop active traction (same law).
        TractionCls = (DcmActiveRimTractionGPUVec if fast_active
                       else DcmActiveRimTractionGPU)
        traction = TractionCls(
            cell_of_node=cell_of_node, ranges=ranges, active=active_mask,
            int_mult=int_mult, R_cell=p.R_cell, z0=p.z_substrate,
            f_act=1.2e-10, f_cap=6.0e-10, ramp_steps=4000, contact_band=0.5,
            neighbour_factor=2.6, max_neighbours=9, integrin_switch_gain=3.0,
            belt_factor=0.25,
            arrest_radius_factor=(arrest_radius_factor if arrest else None),
            arrest_width=arrest_width, arrest_settle_steps=arrest_settle_steps)
        ig.forces.append(traction)

    sim.operations.integrator = ig
    sim.run(0)

    gamma = {"dcm_mem": p.gamma_node, "dcm_inert": p.gamma_node}
    baoab_action, updater = make_baoab_updater(
        kT=p.kT, gamma=gamma, dt=p.dt, seed=p.seed + 1)
    sim.operations.updaters.append(updater)

    return dict(
        sim=sim, cell_of_node=cell_of_node, ranges=ranges, faces=faces,
        face_cell=face_cell, n_cells=n_max, n_max=n_max, n_active=n_active,
        nv=nv, ne=ne, mean_edge=mean_edge, V0=V0, turgor=turgor, contact=contact,
        substrate=substrate, traction=traction, baoab=baoab_action, gamma=gamma,
        p=p, area_per_node=area_per_node, centers=centers, verts0=verts0,
        active=active_mask, int_mult=int_mult)


class GpuProliferationUpdater(hoomd.custom.Action):
    """Rim-biased, contact-inhibited, necrosis-gated division on the GPU pool.

    NO mesh split (see the module-level pool note): each ``act`` activates a parked
    pool cell as the daughter of a dividing rim cell. SLOW relative to spreading
    (PI 2026-06-11: the cell-cycle timescale ≫ the spreading timescale, so over one
    spread phase a rim cell divides ~0–2 times) — small ``p_div`` + a Periodic
    trigger at ``div_every`` so division BOOSTS A/A₀ physically without stacking /
    exploding it.

    Each fired ``act`` (cadence gated by ``div_every``):
      1. read centroids of ACTIVE pool cells;
      2. RIM = 3D convex-hull vertices of the active centroids (free-edge /
         contact-inhibition geometry) AND in substrate contact AND non-necrotic
         (necrosis via the LOD classifier's ``cell_inert`` core mask, if wired —
         a frozen/necrotic core cell never divides);
      3. each eligible rim cell divides with prob ``p_div``: pop the first parked
         pool cell, copy an undeformed icosphere into its node range placed
         ``2R + gap`` OUTWARD from the parent (gapped → no t=0 contact-core overlap
         → BAOAB-safe), lift it ≥ z0+R, then flip ON its ``cell_of_node`` (tent),
         ``face_cell`` (turgor), and traction ``active`` mask. Daughter is
         force-free at activation (undeformed, no overlap).

    The turgor's ``face_cell`` and the tent's ``cell_of_node`` are the SAME mutable
    arrays the live forces read each step, so editing them in place is the whole
    activation — no snapshot rebuild, no tag renumber.

    Args:
        handles: the ``build_gpu_spheroid_prolif`` dict (sim, active, ranges,
            cell_of_node, face_cell, turgor, traction, verts0, nv, p).
        p_div: per-eligible-rim-cell division prob per cadence (SMALL — slow).
        div_every: trigger period in steps (LARGE — slow vs spreading).
        gap_factor: daughter placed 2R + gap_factor·R outward (no overlap).
        cell_inert: optional (n_max,) bool LOD core mask — inert ⇒ no division.
        rng_seed: RNG seed.
    """

    def __init__(self, *, handles: dict, p_div: float, gap_factor: float = 0.4,
                 cell_inert: np.ndarray | None = None,
                 necrotic_depth_um: float = 150.0, rng_seed: int = 99) -> None:
        super().__init__()
        self.h = handles
        self.sim = handles["sim"]
        self.active = handles["active"]              # (n_max,) bool, mutable
        self.ranges = handles["ranges"]
        self.cell_of_node = handles["cell_of_node"]  # (N,) int, mutable (tent)
        # turgor reads face_cell off the live force object → edit THAT array.
        self.face_cell = handles["turgor"].face_cell  # (F,) int, mutable
        self.nf_per_cell = handles["faces"].shape[0] // handles["n_max"]
        self.traction = handles.get("traction")
        self.verts0 = np.asarray(handles["verts0"], dtype=np.float64)
        self.nv = int(handles["nv"])
        p = handles["p"]
        self.R = float(p.R_cell)
        self.z0 = float(p.z_substrate)
        self.contact_band = 0.6
        self.p_div = float(p_div)
        self.gap = float(gap_factor * p.R_cell)
        # NOTE: the pool-wide LOD ``cell_inert`` mask is NOT used to gate division.
        # It is classified over ALL n_max pool cells including the PARKED dormant
        # cells (at z≈60·R, far off), which corrupt the cluster centroid +
        # per-cell depth-from-surface so the live cells spuriously read "deep" →
        # necrotic → every live cell inert → division never fires. We compute the
        # necrosis gate over the LIVE subset only (see ``_live_necrotic``).
        self.cell_inert = cell_inert  # kept for back-compat / introspection only
        self.necrotic_depth_m = float(necrotic_depth_um) * 1.0e-6
        self.n_necrotic_live = 0     # last live-cell necrotic count (diagnostic)
        self._rng = np.random.default_rng(rng_seed)
        self.n_divisions = 0
        self._sim = None

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def _face_range(self, c: int):
        lo = c * self.nf_per_cell
        return lo, lo + self.nf_per_cell

    def _live_mask(self) -> np.ndarray:
        """(n_max,) bool: cell c is LIVE iff its first node's cell_of_node ≥ 0.

        ``cell_of_node`` is the liveness AUTHORITY (the tent + turgor governing
        signal) — NOT the traction ``active`` mask, which the LOD classifier shares
        and overwrites for isolated parked cells. Reading liveness here keeps the
        proliferation accounting correct under LOD.
        """
        first_node = np.array([self.ranges[c][0] for c in range(len(self.ranges))])
        return self.cell_of_node[first_node] >= 0

    @staticmethod
    def live_necrotic_mask(cents_live: np.ndarray, cluster_cen: np.ndarray,
                           necrotic_depth_m: float) -> np.ndarray:
        """Necrosis mask over the LIVE centroids ONLY (radial-depth proxy).

        Mirrors the LOD classifier's necrosis rule (depth from the cluster
        surface > ``necrotic_depth``), but computed on the LIVE-cell centroids
        and their OWN centroid — so the parked dormant pool cells never enter the
        surface-radius / depth estimate. At small spheroids (R_surface ≪ 150 µm)
        every depth is below threshold → necrotic mask all-False → necrosis = 0.

        Args:
            cents_live: (n_live, 3) live-cell centroids.
            cluster_cen: (3,) live-cell cluster centroid.
            necrotic_depth_m: depth-from-surface necrosis onset [m].

        Returns:
            (n_live,) bool — True where the live cell is in the necrotic core.
        """
        if cents_live.shape[0] == 0:
            return np.zeros(0, dtype=bool)
        r_cell = np.linalg.norm(cents_live - cluster_cen, axis=1)
        r_surface = float(r_cell.max())
        depth = r_surface - r_cell
        return depth > necrotic_depth_m

    def act(self, timestep: int) -> None:  # noqa: D401
        live = self._live_mask()
        free = np.where(~live)[0]
        if free.size == 0:
            return  # pool exhausted
        live_ids = np.where(live)[0]
        if live_ids.size < 4:
            return
        with self._sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        perm = np.argsort(tag)
        pos_g = pos[perm]
        cents = np.array([pos_g[self.ranges[int(c)][0]:self.ranges[int(c)][1]].mean(0)
                          for c in live_ids])
        cluster_cen = cents.mean(0)

        # NECROSIS gate over the LIVE subset ONLY (parked pool excluded). Indexed
        # by position in ``live_ids``, parallel to ``cents`` / ``rim_k``.
        necrotic_live = self.live_necrotic_mask(
            cents, cluster_cen, self.necrotic_depth_m)
        self.n_necrotic_live = int(necrotic_live.sum())

        # RIM = convex-hull vertices of the LIVE centroids (free-edge geometry).
        try:
            from scipy.spatial import ConvexHull
            hull = ConvexHull(cents)
            rim_k = np.unique(hull.vertices)
        except Exception:  # noqa: BLE001 — degenerate (coplanar/too few) → no division
            return

        zc = self.z0 + self.contact_band * self.R
        free_ptr = 0                                   # next free pool slot to claim
        divided = False
        for k in rim_k:
            c = int(live_ids[k])
            if free_ptr >= free.size:
                break  # pool exhausted this act
            if bool(necrotic_live[k]):
                continue  # necrotic core (live-cell depth gate) — no division
            lo, hi = self.ranges[c]
            if not bool((pos_g[lo:hi][:, 2] < zc).any()):
                continue  # not in substrate contact — no free basal edge
            if self._rng.random() >= self.p_div:
                continue
            daughter = int(free[free_ptr])             # claim a distinct free slot
            free_ptr += 1
            outward = cents[k] - cluster_cen
            nrm = float(np.linalg.norm(outward))
            if nrm < 1e-12:
                outward = self._rng.standard_normal(3)
                nrm = float(np.linalg.norm(outward))
            outward = outward / nrm
            new_center = cents[k] + (2.0 * self.R + self.gap) * outward
            new_center[2] = max(new_center[2], self.z0 + self.R)
            dlo, dhi = self.ranges[daughter]
            pos_g[dlo:dhi] = self.verts0 + new_center      # undeformed → force-free
            self.cell_of_node[dlo:dhi] = daughter          # tent: ON (+ liveness)
            flo, fhi = self._face_range(daughter)
            self.face_cell[flo:fhi] = daughter             # turgor: ON
            self.active[daughter] = True                   # traction: ON
            if self.traction is not None:
                self.traction.int_mult[daughter] = 1.0
            self.n_divisions += 1
            divided = True
        if divided:
            pos_back = np.empty_like(pos)
            pos_back[perm] = pos_g
            with self._sim.state.cpu_local_snapshot as snap:
                np.asarray(snap.particles.position)[:] = pos_back


# ===========================================================================
# PROLIF-AWARE activity-LOD — classify over LIVE cells ONLY
# ===========================================================================
# WHY (Problem 2, PI 2026-06-11 night). The plain ``attach_activity_lod`` wires a
# ``DcmActivityLOD`` (cell/dcm_gpu_lod.py, FROZEN) that classifies over ALL
# ``n_cells`` ranges. On the PROLIFERATION pool that is ``n_max`` — including the
# PARKED dormant reserve (``cell_of_node = −1``, sitting at z≈60·R, far off the
# cluster). Those parked cells corrupt every per-cell statistic the classifier
# uses: the cluster centroid is dragged toward the parked corner, and — fatally —
# the surface radius ``r_surface = max(‖cent − cluster_cen‖)`` is set by a parked
# cell ~hundreds of µm away, so EVERY live cell reads as "deep" (depth ≫ necrotic
# onset) → spuriously NECROTIC → INERT → its traction ``active`` flag is cleared.
# Result: the prolif build classifies far too few cells active (measured 18/30 vs
# the plain active build's 28/30) and its rim traction is much weaker → A/A₀
# collapses to ~1.19 vs the plain active ~1.55, even with proliferation OFF.
#
# THE FIX (additive, frozen LOD untouched). Subclass ``DcmActivityLOD`` and run the
# EXACT SAME classification rule, but over the LIVE subset only (``cell_of_node ≥
# 0`` — the necrosis fix already applied to the run-script's counting). Live
# centroid / live surface-radius / live crowding / live substrate-contact → no
# parked-pool contamination. The full-pool masks are then written so parked cells
# stay INERT + non-active (they must, they are dormant) and live cells get the
# correct active/inert split. With proliferation OFF (or 0 divisions) this makes
# the prolif build's active set + A/A₀ MATCH the plain active build at the same N.


class LiveCellActivityLOD(DcmActivityLOD):
    """Activity-LOD classifier that classifies over LIVE pool cells ONLY.

    Identical classification PHYSICS to :class:`~ffn_sim.cell.dcm_gpu_lod.DcmActivityLOD`
    (periphery-OR-substrate-contact AND non-necrotic ⇒ active; everything else
    frozen), but the centroid, cluster surface radius, crowding, substrate-contact
    and necrosis-depth are all computed over the LIVE cells of the proliferation
    pool (``cell_of_node ≥ 0``) — never the parked dormant reserve. Parked cells
    are forced INERT + non-active (they are dormant: undeformed, force-free, off
    the substrate). Live cells get the correct rim-active / core-inert split.

    This is the prolif-pool counterpart of ``attach_activity_lod``; it exists so
    that, with proliferation OFF or 0 divisions, the prolif build reproduces the
    plain ``build_gpu_spheroid --active`` active-cell set + A/A₀ (Problem 2 fix).
    The frozen LOD module is not modified — this only OVERRIDES ``act``.
    """

    def act(self, timestep: int) -> None:  # noqa: D401
        if self._last >= 0 and (timestep - self._last) < self.cfg.cadence:
            return
        self._last = timestep
        cfg = self.cfg
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        perm = np.argsort(tag)
        pos_g = pos[perm]                            # pos_g[global tag] = position

        # LIVE = cells whose first node carries a non-negative cell_of_node (the
        # tent + turgor liveness authority; the parked reserve is −1). All
        # statistics below are computed over LIVE cells ONLY so the parked pool
        # (z≈60·R, far off) never corrupts the centroid / surface radius / depth.
        con = self.cell_of_node
        first_node = np.array([lo for (lo, _hi) in self.ranges])
        live = con[first_node] >= 0
        live_ids = np.flatnonzero(live)

        # default: everything inert / non-active (the right state for parked cells)
        active_full = np.zeros(self.n_cells, dtype=bool)
        necrotic_full = np.zeros(self.n_cells, dtype=bool)

        if live_ids.size >= 1:
            cents = np.array([pos_g[self.ranges[int(c)][0]:self.ranges[int(c)][1]].mean(0)
                              for c in live_ids])
            cluster_cen = cents.mean(0)

            # crowding among LIVE centroids only (rim = few live neighbours)
            d2 = np.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
            within = d2 < (cfg.neighbour_factor * self.R) ** 2
            np.fill_diagonal(within, False)
            crowd = within.sum(axis=1)

            # substrate contact (any basal node within the contact band)
            zc = self.z0 + cfg.contact_band * self.R
            in_contact = np.array(
                [bool((pos_g[self.ranges[int(c)][0]:self.ranges[int(c)][1]][:, 2]
                       < zc).any()) for c in live_ids])

            # necrotic = deep from the LIVE-cluster surface (radial-depth proxy).
            # r_surface is the max LIVE centroid distance → no parked-pool inflation.
            r_cell = np.linalg.norm(cents - cluster_cen, axis=1)
            r_surface = float(r_cell.max()) if live_ids.size else 0.0
            depth = r_surface - r_cell
            necrotic = depth > (cfg.necrotic_depth_um * 1.0e-6)

            periphery = crowd <= cfg.max_neighbours
            active = (periphery | in_contact) & (~necrotic)

            active_full[live_ids] = active
            necrotic_full[live_ids] = necrotic

        inert_full = ~active_full

        # capture anchors for cells that are NEWLY inert so they pin at their
        # CURRENT position (live cells that just went inert; parked cells too —
        # harmless, they are already at rest).
        newly_inert = inert_full & (~self.cell_inert)
        for c in np.flatnonzero(newly_inert):
            lo, hi = self.ranges[int(c)]
            self.anchor[lo:hi] = pos_g[lo:hi]
        self.cell_inert[:] = inert_full
        if self.traction_active is not None:
            self.traction_active[:] = active_full

        # diagnostics: report over the LIVE cells (consistent with the run-script's
        # live-cell necrosis counting), not the full n_max pool.
        n_live = int(live_ids.size)
        self.n_active = int(active_full.sum())
        self.n_inert = max(0, n_live - self.n_active)
        self.n_necrotic = int(necrotic_full.sum())
        self.active_frac = self.n_active / max(1, n_live)
        self.necrotic_frac = self.n_necrotic / max(1, n_live)


def attach_activity_lod_prolif(handles: dict, cfg=None):
    """Wire the LIVE-cell activity-LOD onto a ``build_gpu_spheroid_prolif`` build.

    Prolif-pool counterpart of ``dcm_gpu_lod.attach_activity_lod``: same overdamped
    freeze force + the SAME classification physics, but the classifier is the
    :class:`LiveCellActivityLOD` (classifies over LIVE cells only, so the parked
    dormant reserve cannot mis-mark live cells necrotic/inert — Problem 2 fix).

    ``handles`` is the dict from ``build_gpu_spheroid_prolif``. Adds ``lod_freeze``,
    ``lod_updater``, ``cell_inert``, ``anchor``, ``lod_kfreeze`` (mirrors
    ``attach_activity_lod``), sharing the mutable masks with the freeze force and
    flipping the traction ``active`` mask. Returns the updated ``handles``.
    """
    import hoomd

    from ffn_sim.cell.dcm_gpu_lod import (
        DcmLodFreezeForce,
        ResolvedLOD as _ResolvedLOD,
    )

    cfg = cfg or _ResolvedLOD()
    sim = handles["sim"]
    p = handles["p"]
    n_cells = handles["n_cells"]          # = n_max (the full pool)
    cell_of_node = handles["cell_of_node"]
    ranges = handles["ranges"]
    N = sum(hi - lo for (lo, hi) in ranges)

    cell_inert = np.zeros(n_cells, dtype=bool)
    anchor = np.zeros((N, 3), dtype=np.float64)

    # overdamped CFL: k_freeze·dt/γ < 1. Clamp to 0.5·γ/dt for margin (same as the
    # frozen attach_activity_lod).
    gamma = p.gamma_node
    k_max = 0.5 * gamma / p.dt
    k_freeze = min(cfg.k_freeze, k_max)

    freeze = DcmLodFreezeForce(
        cell_of_node=cell_of_node, cell_inert=cell_inert, anchor=anchor,
        k_freeze=k_freeze)
    sim.operations.integrator.forces.append(freeze)
    sim.run(0)

    traction = handles.get("traction")
    traction_active = traction.active if traction is not None else None

    lod = LiveCellActivityLOD(
        cfg=cfg, cell_of_node=cell_of_node, ranges=ranges, n_cells=n_cells,
        R_cell=p.R_cell, z0=p.z_substrate, cell_inert=cell_inert, anchor=anchor,
        traction_active=traction_active)
    updater = hoomd.update.CustomUpdater(
        action=lod, trigger=hoomd.trigger.Periodic(cfg.cadence))
    sim.operations.updaters.append(updater)

    handles["lod_freeze"] = freeze
    handles["lod_updater"] = lod
    handles["cell_inert"] = cell_inert
    handles["anchor"] = anchor
    handles["lod_kfreeze"] = k_freeze
    return handles
