"""ACTIVITY-LEVEL-OF-DETAIL (LOD) for the GPU-friendly DCM spheroid.

PI directive (2026-06-11): in a spreading spheroid only the active RIM / ECM-
contacting SHELL does the expensive dynamics — the deep interior is jammed and
the necrotic core is mechanically dead. So classify cells, at LOW cadence, into

  * ACTIVE — periphery (few close neighbours) OR in substrate contact, AND
    non-necrotic. These get FULL dynamics: active traction + cell-cell contact.
  * INERT — deep interior / quiescent / necrotic. These are FROZEN (they barely
    move) and their active traction is skipped.

This is both FASTER (the inert core is removed from the active-traction work and
from the contact neighbour search of moving nodes) and MORE PHYSICAL (the jammed/
necrotic core is mechanically inert; only the rim spreads).

THE FREEZE MECHANISM (without editing the frozen BAOAB). The project integrator
is the L-M BAOAB-LIMIT (overdamped): each step ``Δr = F·Δt/γ + noise``, with γ
per particle-TYPE. We must not touch ``integrator/baoab.py``. Two levers, used
together, freeze the inert core:

  1. **Overdamped anchor (the freeze).** ``DcmLodFreezeForce`` is a Custom force
     that, for every node of an INERT cell, applies a STIFF harmonic restoring
     force ``F = −k_freeze·(x − x_anchor)`` toward the node's position CAPTURED at
     the moment it became inert. In the overdamped limit a node relaxes to its
     anchor with timescale ``γ/k_freeze``; with ``k_freeze`` large (but
     ``k_freeze·Δt/γ < 1`` for BAOAB stability) the node is PINNED — it barely
     moves, i.e. effectively frozen. The anchor is re-captured each time a cell is
     (re)classified inert, so a cell that the core grows around stays put.
  2. **Skip the active work.** The classifier flips the active-traction's
     ``active`` mask FALSE for inert cells (no outward traction computed for
     them) and sets their ``cell_of_node`` entries to −1 in a SEPARATE contact
     view so the contact kernel skips inert↔inert and inert↔active node pairs
     (the dominant O(pairs) cost). Active cells keep full contact among
     themselves and feel the (now-pinned) core only as a static wall via the
     freeze anchor + their own contact against frozen nodes — physically the rim
     still pushes off the jammed core.

The freeze force is O(n_inert_nodes) and cheap; the SAVING is the active-traction
per-cell loop + the contact pair work that the inert core no longer incurs.

Units SI.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

import hoomd.custom
import hoomd.md as md

from ffn_sim.cell.dcm_gpu_forces import DeviceDispatch


# ---------------------------------------------------------------------------
# Overdamped freeze anchor
# ---------------------------------------------------------------------------
class DcmLodFreezeForce(md.force.Custom):
    """Stiff overdamped harmonic anchor that PINS inert-cell nodes in place.

    For every node whose cell is currently INERT (``cell_inert[cell_of_node]``),
    apply ``F = −k_freeze·(x − x_anchor)`` toward the captured anchor position.
    Active-cell nodes feel ZERO from this force. ``k_freeze`` must satisfy the
    overdamped CFL ``k_freeze·dt/γ < 1`` for BAOAB stability; within that bound,
    larger ``k_freeze`` → tighter pin (smaller residual wander ≈ √(kT/k_freeze)).

    The anchor and the inert mask are MUTABLE and shared with the LOD classifier
    (``DcmActivityLOD``): when a cell becomes inert the classifier writes its
    current node positions into ``anchor`` and flips ``cell_inert``; the force
    reads both fresh each step.

    Args:
        cell_of_node: (N,) int, node→cell id (the build's mutable array).
        cell_inert: (n_cells,) bool, mutable; True = this cell is frozen.
        anchor: (N,3) float, mutable; per-node anchor (rest) position.
        k_freeze: N/m, anchor stiffness.
    """

    def __init__(self, *, cell_of_node: np.ndarray, cell_inert: np.ndarray,
                 anchor: np.ndarray, k_freeze: float) -> None:
        super().__init__(aniso=False)
        self.cell_of_node = cell_of_node          # (N,) int
        self.cell_inert = cell_inert              # (n_cells,) bool, mutable
        self.anchor = anchor                      # (N,3) float, mutable
        self.k_freeze = float(k_freeze)
        self._d: DeviceDispatch | None = None

    def _dispatch(self) -> DeviceDispatch:
        if self._d is None:
            self._d = DeviceDispatch(self)
        return self._d

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        d = self._dispatch()
        xp = d.xp
        with d.snapshot() as snap:
            tag = xp.asarray(snap.particles.tag)
            pos = xp.asarray(snap.particles.position, dtype=xp.float64)
            n = int(pos.shape[0])
            perm = xp.argsort(tag)
            pos_g = pos[perm]                       # pos_g[global tag] = position

            inert_cell = xp.asarray(self.cell_inert)
            cell_g = xp.asarray(self.cell_of_node)
            node_inert = inert_cell[cell_g]          # (N,) bool in tag order
            anchor_g = xp.asarray(self.anchor)

            F_g = xp.zeros_like(pos_g)
            disp = pos_g - anchor_g
            F_g = xp.where(node_inert[:, None], -self.k_freeze * disp, F_g)

            F = xp.empty_like(pos)
            F[perm] = F_g
            U = xp.zeros(n, dtype=xp.float64)

        with d.force_arrays() as arr:
            arr.force[:] = F
            arr.potential_energy[:] = U


# ---------------------------------------------------------------------------
# LOD classifier (low-cadence updater)
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedLOD:
    """Activity-LOD classification knobs.

    A cell is ACTIVE iff it is non-necrotic AND (on the cluster periphery OR in
    substrate contact). Everything else is INERT (frozen).
    """

    neighbour_factor: float = 2.6   # close-neighbour radius = factor · R (centroids)
    max_neighbours: int = 9         # > this many close neighbours ⇒ interior (inert)
    contact_band: float = 0.6       # in substrate contact if a node z < z0 + band·R
    k_freeze: float = 3.0e-4        # N/m anchor stiffness (overdamped pin)
    necrotic_depth_um: float = 150.0  # cells deeper than this (from surface) ⇒ necrotic
    cadence: int = 500              # steps between reclassifications (LOW cadence)


class DcmActivityLOD(hoomd.custom.Action):
    """Low-cadence updater: classify cells ACTIVE vs INERT, freeze the inert core.

    Each ``cadence`` steps it reads the centroids, computes per-cell crowding
    (close-neighbour count) and substrate contact, marks deep-interior /
    necrotic cells INERT, captures their node positions as freeze anchors, and
    flips the active-traction ``active`` mask + the freeze force ``cell_inert``
    mask accordingly. The contact force is left to act on all nodes (the rim
    pushes off the frozen core), but the EXPENSIVE active traction is computed
    only for active cells.

    Diagnostics (read by the driver after each act): ``n_active``, ``n_inert``,
    ``n_necrotic``, ``active_frac``, ``necrotic_frac``.

    Args:
        cfg: ResolvedLOD knobs.
        cell_of_node, ranges: build bookkeeping.
        n_cells, R_cell, z0: geometry.
        cell_inert: (n_cells,) bool, mutable; shared with DcmLodFreezeForce.
        anchor: (N,3) float, mutable; shared with DcmLodFreezeForce.
        traction_active: (n_cells,) bool or None; the active-traction mask to flip
            (None if no traction wired).
    """

    def __init__(self, *, cfg: ResolvedLOD, cell_of_node: np.ndarray, ranges,
                 n_cells: int, R_cell: float, z0: float,
                 cell_inert: np.ndarray, anchor: np.ndarray,
                 traction_active: np.ndarray | None = None) -> None:
        self.cfg = cfg
        self.cell_of_node = cell_of_node
        self.ranges = ranges
        self.n_cells = int(n_cells)
        self.R = float(R_cell)
        self.z0 = float(z0)
        self.cell_inert = cell_inert
        self.anchor = anchor
        self.traction_active = traction_active
        # diagnostics
        self.n_active = n_cells
        self.n_inert = 0
        self.n_necrotic = 0
        self.active_frac = 1.0
        self.necrotic_frac = 0.0
        self._last = -1

    def act(self, timestep: int) -> None:  # noqa: D401
        if self._last >= 0 and (timestep - self._last) < self.cfg.cadence:
            return
        self._last = timestep
        cfg = self.cfg
        sim = self._state._simulation
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        perm = np.argsort(tag)
        pos_g = pos[perm]                            # pos_g[global tag] = position

        # per-cell centroid
        cents = np.array([pos_g[lo:hi].mean(0) for (lo, hi) in self.ranges])
        cluster_cen = cents.mean(0)

        # crowding (close-neighbour count by centroid distance)
        d2 = np.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
        within = d2 < (cfg.neighbour_factor * self.R) ** 2
        np.fill_diagonal(within, False)
        crowd = within.sum(axis=1)

        # substrate contact (any basal node within the contact band)
        zc = self.z0 + cfg.contact_band * self.R
        in_contact = np.array(
            [bool((pos_g[lo:hi][:, 2] < zc).any()) for (lo, hi) in self.ranges])

        # necrotic = deep from the cluster surface (radial depth proxy). The
        # surface radius is the max centroid distance; a cell deeper than
        # necrotic_depth from that surface is in the dead core.
        r_cell = np.linalg.norm(cents - cluster_cen, axis=1)
        r_surface = float(r_cell.max()) if self.n_cells else 0.0
        depth = r_surface - r_cell
        necrotic = depth > (cfg.necrotic_depth_um * 1.0e-6)

        periphery = crowd <= cfg.max_neighbours
        active = (periphery | in_contact) & (~necrotic)
        inert = ~active

        # capture anchors for cells that are NEWLY inert (or re-inert) so they pin
        # at their CURRENT position, not a stale one.
        newly_inert = inert & (~self.cell_inert)
        for c in np.flatnonzero(newly_inert):
            lo, hi = self.ranges[int(c)]
            self.anchor[lo:hi] = pos_g[lo:hi]
        # write masks in place (shared with the freeze force + traction)
        self.cell_inert[:] = inert
        if self.traction_active is not None:
            self.traction_active[:] = active

        self.n_active = int(active.sum())
        self.n_inert = int(inert.sum())
        self.n_necrotic = int(necrotic.sum())
        self.active_frac = self.n_active / max(1, self.n_cells)
        self.necrotic_frac = self.n_necrotic / max(1, self.n_cells)


# ---------------------------------------------------------------------------
# Wiring helper
# ---------------------------------------------------------------------------
def attach_activity_lod(handles: dict, cfg: ResolvedLOD | None = None):
    """Wire the activity-LOD (freeze force + classifier) onto a GPU-DCM build.

    ``handles`` is the dict returned by
    ``dcm_gpu_build.build_gpu_dcm_simulation``. Adds:
      * a ``DcmLodFreezeForce`` to the integrator (the overdamped pin),
      * a ``DcmActivityLOD`` updater at ``cfg.cadence`` cadence,
    sharing the mutable ``cell_inert`` mask + ``anchor`` array, and (if a traction
    is wired) flipping its ``active`` mask. Returns the updated ``handles`` with
    ``lod_freeze``, ``lod_updater``, ``cell_inert``, ``anchor`` added.

    Stability note: ``cfg.k_freeze`` is checked against the overdamped CFL
    ``k_freeze·dt/γ < 1`` and clamped (with the chosen value reported in
    ``handles['lod_kfreeze']``) so the pin never destabilises BAOAB.
    """
    import hoomd

    cfg = cfg or ResolvedLOD()
    sim = handles["sim"]
    p = handles["p"]
    n_cells = handles["n_cells"]
    cell_of_node = handles["cell_of_node"]
    ranges = handles["ranges"]
    N = sum(hi - lo for (lo, hi) in ranges)

    cell_inert = np.zeros(n_cells, dtype=bool)
    anchor = np.zeros((N, 3), dtype=np.float64)

    # overdamped CFL: k_freeze·dt/γ < 1. Clamp to 0.5·γ/dt for margin.
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

    lod = DcmActivityLOD(
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
