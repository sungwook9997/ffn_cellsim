"""GPU-resident, leak-free mechanistic lamellipodium for the GPU-friendly DCM.

This is the production graft of the per-rim-cell lamellipodium (the engine that
spreads a cell by protrusion + FA clutch + traction, NOT by a prescribed body
force) onto ``cell/dcm_gpu_build.build_gpu_dcm_simulation`` — the build that runs
GPU-resident (no native md.mesh). It replaces the coarse ``DcmActiveRimTraction*``
body-force proxy, which provably contracts rather than spreads (a body force
cannot extend a cell's basal area; diagnosis 2026-06-14, PI-confirmed).

Why a NEW module (vs ``cell/dcm_lamellipodium.py``). The legacy engine grows the
actin sheet by ``sim.state.set_snapshot`` topology appends (LeadingEdgeNucleation)
and ``md.bond`` appends (FrontClutch). Both force a HOOMD State rebuild every
batch — the documented rebuild memory leak (tens of MB/rebuild → GBs) and a GPU→
CPU sync that kills GPU residency. This module is the dormant-POOL twin (same
pattern as the Phase-2 remesh ``DcmRemeshUpdater`` and the proliferation pool):

  * a PRE-ALLOCATED dormant ``actin_lamel`` bead pool is added to the snapshot up
    front (``cell_of_node = −1`` so the node-face / tent contact skips it; typeid
    ``actin_lamel`` so the settling body force skips it; parked far off → force-
    free). The front ADVANCES by ACTIVATING a parked bead at a new force-free
    position (a local-snapshot position write), never by a topology append. No
    rebuild ⇒ no leak ⇒ GPU-resident.

  * the FA molecular clutch is a per-bead ANCHOR SPRING, not a ``md.bond`` to an
    explicit ligand: a bead is "born gripped" at its activation site (the substrate
    point it polymerised onto), and :class:`ActinClutchAnchorGPU` holds it there
    (``F = −k_clutch·(r − anchor)``, capped) — the engaged integrin clutch (KU-2.4
    k_int) pinning the front to the rigid dish. This folds the explicit ``ligand``
    grid + ``LigandPin`` + grip-``md.bond`` of the legacy engine into one anchor
    array (the dense-ligand limit: a ligand exactly under each gripped bead), so
    there is NO bond topology to mutate.

  * traction is transmitted by :class:`LamellipodialTractionTetherGPU` exactly as
    in the legacy engine: each rim cell's LEADING (outward-radial, basal) membrane
    node is pulled (capped harmonic) toward the nearest OUTWARD clutch-anchored
    actin bead of its OWN cell. The anchored actin is the fixed point, so the cell
    body is pulled outward onto its own substrate-gripped lamellipodium — real
    traction carried by the actin/clutch, not a prescribed force.

The ratchet (mechanistic, not a body force). Each low-cadence tick
(:class:`GpuLamellipodiumAdvance`): for each leading membrane node, if the node
has CAUGHT UP to its current front bead (within ``catch_factor·ℓ₀``) it advances
the front — activating a new bead one ℓ₀ further out along the node's outward
radial, with probability ``p = v_front·Δt_tick·S/ℓ₀`` (the band-matched
polymerisation rate, lit single-cell front 3–12 µm/min). So the front velocity is
the one constitutive input (as in continuum spreading models / the validated
engine); everything downstream — footprint shape, traction magnitude, how far the
cell flattens — is emergent BAOAB physics. The front only advances where the
membrane has caught up = the clutch ratchet.

Faithful-vs-simplified (surfaced to PI). Kept explicit (NOT lumped): actin beads as
real BAOAB DOFs, the clutch anchor spring (KU-2.4), the traction tether, the
polymerisation-rate-gated ratchet. Simplified vs ``dcm_lamellipodium.py``:
(1) "born gripped" — a bead anchors at activation rather than nucleate-then-grip a
tick later (it is activated basal, in the contact band, where it would grip next
tick anyway); (2) per-leading-node radial lanes rather than free Arp2/3 branching
(the front is a membrane-tracked ratchet, the same mesoscale proxy the legacy
``LeadingEdgeNucleationUpdater`` uses). These are mesoscale construction choices,
not a return to the lumped body force — protrusion, clutch, and traction are all
explicit. Units SI.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

import hoomd
import hoomd.md as md

from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_forces import on_gpu, DeviceDispatch


def _scatter_add(xp, a, idx, v) -> None:
    """``a[idx] += v`` accumulating duplicate indices, on numpy OR cupy.

    numpy has ``np.add.at``; cupy's equivalent is ``cupyx.scatter_add`` (plain
    fancy-index ``+=`` does NOT accumulate duplicates on either backend).
    """
    if xp is np:
        np.add.at(a, idx, v)
    else:
        import cupyx
        cupyx.scatter_add(a, idx, v)


# ---------------------------------------------------------------------------
# Resolved lamellipodium parameters (on top of ResolvedGpuDCM)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedGpuLamellipodium:
    """Per-rim-cell lamellipodium constants for the GPU-friendly DCM build.

    Defaults band-match the validated single-cell engine (``spreading_drive``):
    v_front 3–12 µm/min, clutch stiffness in the KU-2.4 integrin band, tether in
    the overdamped-CFL safe band, per-node traction cap = 5·(60 Pa·area_per_node)
    (lit MCF7 traction stress, Gil-Redondo 2023). Drag for the actin type is
    DERIVED at build from the cytoplasm viscosity (physiological-baseline rule), so
    it is not carried here.
    """

    # --- protrusion (polymerisation ratchet) ---
    v_front: float = 6.0e-6 / 60.0   # m/s band-matched front velocity (lit 3–12 µm/min)
    actin_rest_length: float = 0.5e-6   # ℓ₀ monomer step [m]
    S_kinetic: float = 1.0           # physiological real-time clock (new foundation)
    catch_factor: float = 1.2        # advance a lane when the node is within
                                     # catch_factor·ℓ₀ of its front bead (the ratchet)
    seed_offset: float = 1.0         # initial seed bead is seed_offset·ℓ₀ outward

    # --- rim-cell detection (the spreading interface = basal contact cells) ---
    rim_contact_band: float = 1.5    # rim cell if centroid z < z0 + band·R (basal-
                                     # contact; ~bottom 1-2 layers of a touching ball)
    rim_neighbour_factor: float = 0.0   # unused (kept for signature compat)
    rim_max_neighbours: int = 0      # unused

    # --- leading-node classification (which membrane nodes are pulled) ---
    lead_frac: float = 0.0           # node leading if outward-proj >= lead_frac·|rel_xy|
    basal_band_factor: float = 0.3   # basal band half-width = factor·R about z_basal
                                     # (the substrate CONTACT zone — must be ≲ tether
                                     # reach so a z_basal-seeded bead is re-findable by
                                     # its leading node, else the lane re-seeds forever)
    seed_basal_offset: float = 0.3e-6   # actin sheet at z0 + this (basal)

    # --- FA molecular clutch (anchor spring) ---
    # The per-bead anchor is the engaged FA ENSEMBLE (~N_int integrin clutches in a
    # mature adhesion, k_single ≈ 2e-4 N/m KU-2.4) attaching the gripped actin to
    # the RIGID glass dish — so the effective anchor is ~N_int·k_single ≈ 5e-2 N/m,
    # MUCH stiffer than the traction tether (k_tether=4e-3). It MUST be: a soft
    # anchor would let the tether's Newton-3 reaction drag the actin inward (both at
    # equal drag → they meet in the middle, no net protrusion). A stiff anchor makes
    # the substrate-gripped actin the FIXED point the membrane is pulled OUT to.
    # CFL ok: dt<2γ/k ⇒ at γ=2.22e-4, k=5e-2 → dt<8.9e-3 (spread dt 1e-4 safe).
    k_clutch: float = 5.0e-2         # N/m FA-ensemble clutch anchor (≈N_int·KU-2.4)
    clutch_cap: float = 1.0e-7       # N per-bead anchor cap (>> tether_cap so it holds)

    # --- actin → leading-membrane traction tether ---
    k_tether: float = 4.0e-3         # N/m actin → leading-membrane tether spring
    tether_cap: float = 5.0e-9       # N per-node cap = 5·(60 Pa·area_per_node) lit MCF7
    tether_radius: float = 3.0e-6    # m a membrane node tethers to actin within this

    # --- pool / cadence ---
    pool_per_cell: int = 60          # dormant actin beads pre-allocated per rim cell
    batch_steps: int = 50            # BAOAB steps between ratchet ticks
    seed: int = 7
    extras: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Shared mutable lamellipodium state (read by the forces, written by the updater)
# ---------------------------------------------------------------------------
class LamelState:
    """Mutable per-actin-bead bookkeeping shared between the updater and forces.

    The dormant actin pool occupies the CONTIGUOUS global-tag block
    ``[a0, a0 + n_pool)`` (appended last in the snapshot). All arrays are indexed
    by pool position ``i`` (global tag ``a0 + i``):

      * ``active[i]``  — bead i is an activated (live) lamellipodial DOF.
      * ``cell_of[i]`` — which rim cell bead i serves (−1 dormant).
      * ``anchor[i]``  — (3,) clutch anchor point (the substrate site it gripped).

    The activation pointer ``n_used`` is a simple contiguous high-water mark (we
    never deactivate a lamellipodial bead within a run — trailing beads stay as the
    adhered sheet), so ``active = arange < n_used`` and the pool is exhausted when
    ``n_used == n_pool``.
    """

    def __init__(self, *, a0: int, n_pool: int) -> None:
        self.a0 = int(a0)
        self.n_pool = int(n_pool)
        self.n_used = 0
        self.cell_of = np.full(n_pool, -1, dtype=np.int64)
        self.anchor = np.zeros((n_pool, 3), dtype=np.float64)
        self.pool_exhausted = 0

    @property
    def active_mask(self) -> np.ndarray:
        m = np.zeros(self.n_pool, dtype=bool)
        m[: self.n_used] = True
        return m

    def activate(self, cell: int, anchor_xyz: np.ndarray) -> int | None:
        """Activate the next dormant bead for ``cell`` anchored at ``anchor_xyz``.

        Returns the pool index activated, or None if the pool is exhausted.
        """
        if self.n_used >= self.n_pool:
            self.pool_exhausted += 1
            return None
        i = self.n_used
        self.cell_of[i] = int(cell)
        self.anchor[i] = anchor_xyz
        self.n_used += 1
        return i


# ---------------------------------------------------------------------------
# Snapshot helper: append the dormant actin pool + register the actin type
# ---------------------------------------------------------------------------
def add_actin_pool_to_snapshot(snap, cell_of_node: np.ndarray, *,
                               n_actin: int, R_cell: float, z_substrate: float):
    """Append a dormant ``actin_lamel`` bead pool to a gsd Frame in place.

    Adds ``actin_lamel`` to ``snap.particles.types`` (typeid 2), appends ``n_actin``
    parked beads (parked in a compact grid well above/beside the cluster so they
    are outside every contact cutoff and the substrate well → force-free), and
    extends ``cell_of_node`` with −1 for each (so node-face / tent contact and the
    settling force all skip them). Returns ``(cell_of_node_ext, a0)`` where ``a0``
    is the first actin global tag.
    """
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    tid = np.asarray(snap.particles.typeid, dtype=np.uint32)
    a0 = pos.shape[0]

    types = list(snap.particles.types)
    if "actin_lamel" not in types:
        types.append("actin_lamel")
    actin_tid = types.index("actin_lamel")

    # park in a compact grid high above the cluster (force-free; activated later).
    side = int(np.ceil(n_actin ** (1.0 / 3.0)))
    gx, gy, gz = np.meshgrid(np.arange(side), np.arange(side), np.arange(side))
    grid = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])[:n_actin]
    park0 = np.array([np.abs(pos[:, 0]).max() + 10.0 * R_cell,
                      np.abs(pos[:, 1]).max() + 10.0 * R_cell,
                      z_substrate + 12.0 * R_cell])
    actin_pos = park0 + grid * (2.0 * R_cell)

    new_pos = np.concatenate([pos, actin_pos], axis=0)
    new_tid = np.concatenate([tid, np.full(n_actin, actin_tid, dtype=np.uint32)])
    cof_ext = np.concatenate([np.asarray(cell_of_node, np.int64),
                              np.full(n_actin, -1, dtype=np.int64)])

    # grow the box if the parked pool sticks out
    L = float(np.abs(new_pos).max() * 2.2 + 4.0 * R_cell)
    if L > snap.configuration.box[0]:
        snap.configuration.box = [L, L, L, 0, 0, 0]

    snap.particles.N = int(new_pos.shape[0])
    snap.particles.types = types
    snap.particles.position = new_pos
    snap.particles.typeid = new_tid
    snap.particles.mass = np.ones(new_pos.shape[0], dtype=np.float64)
    return cof_ext, a0


# ---------------------------------------------------------------------------
# Rim-cell detection (periphery + substrate contact)
# ---------------------------------------------------------------------------
def detect_rim_cells(centers: np.ndarray, *, R: float, z0: float,
                     contact_band: float, neighbour_factor: float = 0.0,
                     max_neighbours: int = 0) -> np.ndarray:
    """Cell indices that form the spreading interface = the BASAL (substrate-
    contacting) cells (centroid within ``contact_band·R`` of z0).

    The spreading wetting front is the spheroid-substrate contact disk; every cell
    that touches the dish is a candidate lamellipodium cell. We DO NOT additionally
    require "few neighbours" (the legacy periphery heuristic) — in a packed
    substrate-touching ball the bottom-contact cells have MANY neighbours, so that
    gate wrongly excluded exactly the cells that crawl. Periphery-vs-centre is then
    self-selected by the updater/tether: a central contact cell's outward-radial is
    ~ill-defined (its centroid ≈ the spheroid axis) so it grows few leading nodes,
    while a peripheral contact cell crawls outward. ``neighbour_factor`` /
    ``max_neighbours`` are accepted for signature compatibility but unused.

    Falls back to the lowest third of cells by z if none are within the band (so a
    not-quite-touching aggregate still seeds a front).
    """
    n = centers.shape[0]
    basal = np.flatnonzero((centers[:, 2] - z0) <= contact_band * R)
    if basal.size == 0:
        return np.argsort(centers[:, 2])[: max(1, n // 3)].astype(np.int64)
    return basal.astype(np.int64)


# ---------------------------------------------------------------------------
# FORCE 1 — FA molecular-clutch anchor spring (holds the gripped front)
# ---------------------------------------------------------------------------
class ActinClutchAnchorGPU(md.force.Custom):
    """Hold each active ``actin_lamel`` bead at its clutch anchor (the substrate
    site it polymerised onto) — the engaged FA molecular clutch.

    ``F_i = −k_clutch·(r_i − anchor_i)`` (full 3D, capped) for every active bead;
    zero for dormant beads. The anchor (xy at the activation site, z = substrate)
    is the rigid-dish ligand the integrin clutch grips, so the bead is the FIXED
    point the traction tether pulls the membrane toward. Device-dispatched
    (cupy/GPU, numpy/CPU). Reads :class:`LamelState` (``n_used``, ``anchor``) fresh
    each step; the actin pool is the contiguous tag block ``[a0, a0+n_pool)``.
    """

    def __init__(self, *, state: LamelState, k_clutch: float, force_cap: float):
        super().__init__(aniso=False)
        self.st = state
        self.k = float(k_clutch)
        self.cap = float(force_cap)
        self._d: DeviceDispatch | None = None

    def _dispatch(self) -> DeviceDispatch:
        if self._d is None:
            self._d = DeviceDispatch(self)
        return self._d

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        d = self._dispatch()
        xp = d.xp
        n_used = self.st.n_used
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            n = int(pos.shape[0])
            perm = xp.argsort(tag)
            F_g = xp.zeros((n, 3), dtype=xp.float64)
            if n_used > 0:
                pos_g = pos[perm]
                a0 = self.st.a0
                rows = slice(a0, a0 + n_used)
                anchor = xp.asarray(self.st.anchor[:n_used])
                dvec = pos_g[rows] - anchor               # bead → anchor offset
                f = -self.k * dvec
                fmag = xp.sqrt((f * f).sum(axis=1, keepdims=True))
                fmag_safe = xp.where(fmag > 0.0, fmag, 1.0)
                scale = xp.where(fmag > self.cap, self.cap / fmag_safe, 1.0)
                F_g[rows] = f * scale
            F = xp.empty_like(pos)
            F[perm] = F_g
        with d.force_arrays() as arr:
            arr.force[:] = F


# ---------------------------------------------------------------------------
# FORCE 2 — actin → leading-membrane traction tether (the traction transmission)
# ---------------------------------------------------------------------------
class LamellipodialTractionTetherGPU(md.force.Custom):
    """Pull each rim cell's LEADING basal membrane node toward the nearest OUTWARD
    clutch-anchored actin bead of its OWN cell — the traction transmission.

    Device-dispatched, fully vectorised. Each step:
      1. cell centroids (segment-mean over ``cell_of_node`` for membrane nodes) and
         the spheroid in-plane centroid (mean of rim-cell membrane nodes);
      2. a membrane node is LEADING if it is in a rim cell, basal (within the basal
         band of ``z_basal``), and its outward-radial projection ≥ ``lead_frac``;
      3. for each leading node, the nearest active actin bead of the SAME cell that
         is OUTWARD of it (radially farther) and within ``tether_radius``;
      4. ``F = +k_tether·(r_actin − r_node)`` on the node (capped), Newton-3 −F on
         the actin bead — the anchored actin pulls the cell body outward.

    Reads :class:`LamelState` (``n_used``, ``cell_of``) for the active actin and its
    per-bead cell; the actin pool is the contiguous tag block ``[a0, a0+n_pool)``.
    """

    def __init__(self, *, state: LamelState, cell_of_node: np.ndarray,
                 rim_cells: np.ndarray, mem_typeid: int, n_mem: int,
                 z_basal: float, basal_band: float, lead_frac: float,
                 k_tether: float, force_cap: float, tether_radius: float):
        super().__init__(aniso=False)
        self.st = state
        self.cell_of_node = np.asarray(cell_of_node, dtype=np.int64)
        self.rim_cells = np.asarray(rim_cells, dtype=np.int64)
        self.mem_typeid = int(mem_typeid)
        self.n_mem = int(n_mem)
        self.z_basal = float(z_basal)
        self.basal_band = float(basal_band)
        self.lead_frac = float(lead_frac)
        self.k = float(k_tether)
        self.cap = float(force_cap)
        self.radius = float(tether_radius)
        self._n_tethered = 0
        self._d: DeviceDispatch | None = None

    @property
    def n_tethered(self) -> int:
        return self._n_tethered

    def _dispatch(self) -> DeviceDispatch:
        if self._d is None:
            self._d = DeviceDispatch(self)
        return self._d

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        d = self._dispatch()
        xp = d.xp
        n_used = self.st.n_used
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            n = int(pos.shape[0])
            perm = xp.argsort(tag)
            F_g = xp.zeros((n, 3), dtype=xp.float64)
            if n_used > 0:
                pos_g = pos[perm]
                self._compute(xp, pos_g, n_used, F_g)
            F = xp.empty_like(pos)
            F[perm] = F_g
        with d.force_arrays() as arr:
            arr.force[:] = F

    def _compute(self, xp, pos_g, n_used, F_g):
        a0 = self.st.a0
        nm = self.n_mem
        cof = xp.asarray(self.cell_of_node)        # (n_mem,) cell id per membrane node
        rim = xp.asarray(self.rim_cells)
        n_cells = int(cof.max()) + 1 if cof.size else 0

        # cell centroids (segment mean over membrane nodes)
        mem_pos = pos_g[:nm]
        cnt = xp.zeros(n_cells, dtype=xp.float64)
        cx = xp.zeros(n_cells, dtype=xp.float64)
        cy = xp.zeros(n_cells, dtype=xp.float64)
        valid = cof >= 0
        idx = cof[valid]
        _scatter_add(xp, cnt, idx, 1.0)
        _scatter_add(xp, cx, idx, mem_pos[valid, 0])
        _scatter_add(xp, cy, idx, mem_pos[valid, 1])
        cnt_safe = xp.where(cnt > 0, cnt, 1.0)
        cen_x = cx / cnt_safe
        cen_y = cy / cnt_safe

        # spheroid in-plane centroid from rim-cell centroids
        rim_mask_cell = xp.zeros(n_cells, dtype=bool)
        rim_mask_cell[rim] = True
        sx = float(cen_x[rim_mask_cell].mean()) if rim.size else 0.0
        sy = float(cen_y[rim_mask_cell].mean()) if rim.size else 0.0

        # per membrane node: outward dir, projection, basal, leading, in-rim
        node_cell = cof
        in_rim = rim_mask_cell[xp.where(node_cell >= 0, node_cell, 0)] & (node_cell >= 0)
        ccx = cen_x[xp.where(node_cell >= 0, node_cell, 0)]
        ccy = cen_y[xp.where(node_cell >= 0, node_cell, 0)]
        outx = ccx - sx
        outy = ccy - sy
        on = xp.sqrt(outx * outx + outy * outy)
        on_safe = xp.where(on > 1e-18, on, 1.0)
        outx = outx / on_safe
        outy = outy / on_safe
        relx = mem_pos[:, 0] - ccx
        rely = mem_pos[:, 1] - ccy
        proj = relx * outx + rely * outy
        rel_xy = xp.sqrt(relx * relx + rely * rely)
        basal = xp.abs(mem_pos[:, 2] - self.z_basal) <= self.basal_band
        leading = in_rim & basal & (proj >= self.lead_frac * xp.maximum(rel_xy, 1e-18))

        lead_idx = xp.where(leading)[0]
        if lead_idx.size == 0:
            self._n_tethered = 0
            return

        # active actin: positions, cells, anchors
        act_pos = pos_g[a0:a0 + n_used]            # (n_used,3)
        act_cell = xp.asarray(self.st.cell_of[:n_used])
        lp = mem_pos[lead_idx]                      # (L,3)
        lcell = node_cell[lead_idx]                 # (L,)
        lout_x = outx[lead_idx]; lout_y = outy[lead_idx]
        lproj = proj[lead_idx]

        # (L, n_used) distances and masks — brute force (cheap: L~1e3, n_used~few e3)
        dx = act_pos[None, :, 0] - lp[:, None, 0]
        dy = act_pos[None, :, 1] - lp[:, None, 1]
        dz = act_pos[None, :, 2] - lp[:, None, 2]
        dist = xp.sqrt(dx * dx + dy * dy + dz * dz)
        same_cell = act_cell[None, :] == lcell[:, None]
        # actin radial projection from the leading node's cell centroid
        act_proj = ((act_pos[None, :, 0] - lp[:, None, 0]) * lout_x[:, None]
                    + (act_pos[None, :, 1] - lp[:, None, 1]) * lout_y[:, None])
        outward = act_proj > 0.0                    # actin farther out than the node
        ok = same_cell & outward & (dist <= self.radius)
        big = xp.asarray(np.float64(1e30))
        masked = xp.where(ok, dist, big)
        jmin = xp.argmin(masked, axis=1)            # (L,) nearest outward same-cell actin
        has = masked[xp.arange(lead_idx.shape[0]), jmin] < big

        sel = xp.where(has)[0]
        if sel.size == 0:
            self._n_tethered = 0
            return
        ln = lead_idx[sel]
        aj = jmin[sel]
        dvec = act_pos[aj] - lp[sel]                # node → actin (pull outward)
        f = self.k * dvec
        fmag = xp.sqrt((f * f).sum(axis=1, keepdims=True))
        fmag_safe = xp.where(fmag > 0.0, fmag, 1.0)
        scale = xp.where(fmag > self.cap, self.cap / fmag_safe, 1.0)
        f = f * scale
        # scatter onto nodes (+f) and onto actin beads (−f, Newton-3)
        _scatter_add(xp, F_g, ln, f)
        _scatter_add(xp, F_g, a0 + aj, -f)
        self._n_tethered = int(sel.shape[0])


# ---------------------------------------------------------------------------
# UPDATER — the membrane-tracked polymerisation ratchet (dormant-pool activation)
# ---------------------------------------------------------------------------
class GpuLamellipodiumAdvance(hoomd.custom.Action):
    """Advance each rim cell's lamellipodial front by activating dormant actin beads.

    Host-side, LOW cadence (every ``batch_steps``). Each tick, for each rim cell, it
    classifies the cell's LEADING basal membrane nodes (outward-radial), and per
    leading node:
      * if the node has NO front bead yet → SEED one ℓ₀ outward (activate);
      * else if the node has CAUGHT UP to its front bead (within ``catch_factor·ℓ₀``)
        and ``rng < p_advance`` → ADVANCE: activate a new bead one ℓ₀ further out.
    Activation = pop a dormant pool bead, write its position to the activation site
    (basal, z = z_basal) via the local snapshot, and record the clutch anchor (that
    site) + serving cell in :class:`LamelState`. No topology change ⇒ no rebuild.

    ``p_advance = v_front·(batch_steps·dt)·S/ℓ₀`` (clamped) is the band-matched
    polymerisation rate; the front only advances where the membrane caught up = the
    clutch ratchet.
    """

    def __init__(self, *, state: LamelState, cell_of_node: np.ndarray,
                 rim_cells: np.ndarray, mem_typeid: int, n_mem: int,
                 z_basal: float, basal_band: float, lead_frac: float,
                 rest_length: float, catch_factor: float, p_advance: float,
                 tether_radius: float, seed_offset: float, seed: int = 7):
        super().__init__()
        self.st = state
        self.cell_of_node = np.asarray(cell_of_node, dtype=np.int64)
        self.rim_cells = np.asarray(rim_cells, dtype=np.int64)
        self.mem_typeid = int(mem_typeid)
        self.n_mem = int(n_mem)
        self.z_basal = float(z_basal)
        self.basal_band = float(basal_band)
        self.lead_frac = float(lead_frac)
        self.l0 = float(rest_length)
        self.catch = float(catch_factor) * self.l0
        self.p_advance = float(min(max(p_advance, 0.0), 1.0))
        self.radius = float(tether_radius)
        self.seed_offset = float(seed_offset) * self.l0
        self._rng = np.random.default_rng(seed)
        self._sim = None
        self.n_seeded = 0
        self.n_advanced = 0

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        st = self.st
        if st.n_used >= st.n_pool:
            return  # pool exhausted (logged via st.pool_exhausted by activate())
        sim = self._sim
        with sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        perm = np.argsort(tag)
        pos_g = pos[perm]
        nm = self.n_mem
        a0 = st.a0
        cof = self.cell_of_node

        # spheroid in-plane centroid from rim-cell membrane nodes
        rim_node_mask = np.isin(cof, self.rim_cells)
        if not rim_node_mask.any():
            return
        sph = pos_g[:nm][rim_node_mask].mean(axis=0)

        # active actin positions/cells (numpy view)
        n_used = st.n_used
        act_pos = pos_g[a0:a0 + n_used] if n_used else np.empty((0, 3))
        act_cell = st.cell_of[:n_used] if n_used else np.empty(0, np.int64)

        activations = []   # (pool_idx-to-be, position) recorded; positions written below
        new_positions = {}

        for c in self.rim_cells:
            c = int(c)
            cnodes = np.where(cof == c)[0]
            if cnodes.size == 0:
                continue
            npos = pos_g[cnodes]
            cc = npos.mean(axis=0)
            out = np.array([cc[0] - sph[0], cc[1] - sph[1], 0.0])
            on = np.linalg.norm(out)
            if on < 1e-18:
                continue
            out /= on
            rel = npos - cc
            proj = rel[:, 0] * out[0] + rel[:, 1] * out[1]
            rel_xy = np.hypot(rel[:, 0], rel[:, 1])
            basal = np.abs(npos[:, 2] - self.z_basal) <= self.basal_band
            lead = basal & (proj >= self.lead_frac * np.maximum(rel_xy, 1e-18))
            lead_nodes = cnodes[lead]
            if lead_nodes.size == 0:
                continue

            # this cell's active actin
            cmask = act_cell == c
            cact = act_pos[cmask] if n_used else np.empty((0, 3))

            node_proj_c = (npos[:, 0] - cc[0]) * out[0] + (npos[:, 1] - cc[1]) * out[1]
            # this cell's active-actin outward projection from the cell centroid
            if cact.shape[0] > 0:
                proj_a = (cact[:, 0] - cc[0]) * out[0] + (cact[:, 1] - cc[1]) * out[1]
            for ki, ln in enumerate(lead_nodes):
                rp = pos_g[ln]
                node_proj = float(node_proj_c[np.where(cnodes == ln)[0][0]])
                # front bead = cell's actin within tether reach whose outward
                # projection is the LARGEST (outermost) — NOT strictly-outward of
                # the node (that re-seeds once the node catches up). None ⇒ seed.
                fb_proj = None
                if cact.shape[0] > 0:
                    dist = np.linalg.norm(cact - rp, axis=1)
                    near = dist <= self.radius
                    if near.any():
                        cand_proj = np.where(near, proj_a, -np.inf)
                        jm = int(np.argmax(cand_proj))
                        fb_proj = float(proj_a[jm])
                if fb_proj is None:
                    # SEED a front bead one seed_offset outward of the node
                    site = np.array([rp[0] + self.seed_offset * out[0],
                                     rp[1] + self.seed_offset * out[1], self.z_basal])
                    i = st.activate(c, site)
                    if i is None:
                        return
                    new_positions[a0 + i] = site
                    self.n_seeded += 1
                elif (fb_proj - node_proj) < self.catch and self._rng.uniform() < self.p_advance:
                    # ADVANCE: node caught up to the front → polymerise one ℓ₀ further
                    site = np.array([cc[0] + (fb_proj + self.l0) * out[0],
                                     cc[1] + (fb_proj + self.l0) * out[1], self.z_basal])
                    i = st.activate(c, site)
                    if i is None:
                        return
                    new_positions[a0 + i] = site
                    self.n_advanced += 1
                if st.n_used >= st.n_pool:
                    break
            if st.n_used >= st.n_pool:
                break

        if not new_positions:
            return
        # write the activated bead positions back (tag order → local row order)
        pos_back = pos_g.copy()
        for grow, site in new_positions.items():
            pos_back[grow] = site
        pos_local = np.empty_like(pos)
        pos_local[perm] = pos_back
        with sim.state.cpu_local_snapshot as snap:
            np.asarray(snap.particles.position)[:] = pos_local
            # zero the activated beads' velocity (no stale parked velocity)
            vel = np.asarray(snap.particles.velocity)
            for grow in new_positions:
                vel[perm[grow]] = 0.0
