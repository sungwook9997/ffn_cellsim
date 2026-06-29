"""DCM remeshing — Phase 2, step 2b: the low-cadence HOOMD CustomUpdater.

Applies the SWAP/SPLIT/COLLAPSE mutations (``cell/dcm_remesh.py``) to the LIVE
mesh during a run, keeping every edge in [l_min, 3·l_min] so large spreading never
stretches a triangle into a sliver (the H.7 "LJ explosion" — the node-face contact
force ∝ A_face). It is the mesh-side twin of ``ProliferationUpdater``: a host-side
``hoomd.custom.Action`` driven at LOW cadence by a ``CustomUpdater`` trigger,
mutating in place the arrays the live custom forces read —

  * ``DcmTurgorForceGPU.faces`` / ``.face_cell``  (the K1 turgor groups),
  * ``FaceContactForceGPU.faces`` / ``.face_cell`` (the node-vs-face penalty),
  * the SHARED ``cell_of_node`` array (−1 = dormant pool node).

The fixed-tag-space constraint (HOOMD + BAOAB cannot add particles after
``create_state``) is handled by the pre-allocated DORMANT NODE POOL: a SPLIT
activates one parked dormant node at the new edge midpoint, a COLLAPSE returns one
to the pool. The face arrays are plain numpy attributes (NOT tag space) so they
resize freely. No State-rebuild ⇒ no HOOMD rebuild leak.

⚠ CONTRACT / SCOPE. This updater assumes ``turgor.face_cell`` and
``contact.face_cell`` share the SAME (global owner-cell) indexing — true while all
cells are active / no mid-run division has remapped the turgor groups to a
compact-active index. Remesh + division co-running is future work (it needs the
``_rebuild_turgor_faces`` compact remap reconciled with the per-face owner id).
Edge-spring bonds (``md.bond.Harmonic``) are NOT remeshed here — they are the
Phase-3-replaced surface-tension proxy and lie in fixed bond tag space; a
split-activated node carries no edge spring (it is held by turgor + the surrounding
faces). Surfaced to PI as the known step-2 limitation.
"""

from __future__ import annotations

import numpy as np

import hoomd

from ffn_sim.dcm.dcm_remesh import remesh_pass


class DcmRemeshUpdater(hoomd.custom.Action):
    """Low-cadence local remeshing of the live DCM shells.

    Args:
        turgor: ``DcmTurgorForceGPU`` — its ``faces`` / ``face_cell`` are mutated.
        contact: ``FaceContactForceGPU`` or None (node-node tent ⇒ pass None; the
            face arrays do not apply, only the turgor is remeshed).
        cell_of_node: the SHARED (N,) per-node cell id array (−1 = dormant); mutated
            in place so the contact force sees activations/deactivations next step.
        l_min: lower edge-length bound (l_max = 3·l_min); typically ~0.5·mean_edge so
            the start mesh sits mid-band.
        sliver_q: SWAP fires below this face quality (SimuCell3D 0.2).
        max_ops: cap on mutations per ``act`` (bounds host cost at low cadence).
        park: (3,) position to send deactivated pool nodes (default: derived from the
            simulation box corner at attach).
    """

    def __init__(self, *, turgor, contact, cell_of_node, l_min, sliver_q=0.2,
                 max_ops=24, park=None):
        super().__init__()
        self.turgor = turgor
        self.contact = contact
        self.cell_of_node = cell_of_node
        self.l_min = float(l_min)
        self.sliver_q = float(sliver_q)
        self.max_ops = int(max_ops)
        self.park = None if park is None else np.asarray(park, float)
        self._sim = None
        self.totals = dict(swap=0, split=0, collapse=0, pool_exhausted=0)
        self.n_acts = 0

    def attach(self, simulation):  # noqa: D401
        super().attach(simulation)
        self._sim = simulation
        if self.park is None:
            L = float(simulation.state.box.Lx)
            self.park = np.array([0.45 * L, 0.45 * L, 0.45 * L])

    def act(self, timestep: int) -> None:  # noqa: D401
        sim = self._sim
        self.n_acts += 1
        with sim.state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        perm = np.argsort(tag)               # global-tag order (faces index this)
        pos_g = pos[perm]

        cof_before = np.asarray(self.cell_of_node).copy()
        pos2, faces2, cof2, fc2, counts = remesh_pass(
            pos_g, self.turgor.faces, cof_before, self.l_min,
            face_cell=self.turgor.face_cell, max_ops=self.max_ops,
            sliver_q=self.sliver_q, park=self.park)

        if counts["swap"] + counts["split"] + counts["collapse"] == 0:
            return                            # nothing changed — skip the writeback

        # push the new topology to the live forces (numpy attrs, read next step).
        self.turgor.faces = faces2
        self.turgor.face_cell = fc2
        if self.contact is not None:
            self.contact.faces = faces2
            self.contact.face_cell = fc2
        # the contact force holds the SAME cell_of_node ref → mutate in place.
        np.asarray(self.cell_of_node)[:] = cof2

        # write back moved / activated / parked node positions; zero the velocity of
        # nodes whose pool membership changed (activated split nodes, parked collapse
        # nodes) so no stale dormant velocity is injected (overdamped ⇒ re-thermalised
        # within a step anyway).
        changed_cof = np.flatnonzero(cof2 != cof_before)
        pos_back = np.empty_like(pos)
        pos_back[perm] = pos2
        with sim.state.cpu_local_snapshot as snap:
            np.asarray(snap.particles.position)[:] = pos_back
            if changed_cof.size:
                vel = np.asarray(snap.particles.velocity)
                vel[perm[changed_cof]] = 0.0

        for k in self.totals:
            self.totals[k] += counts[k]


def attach_remesh_updater(handles: dict, *, period: int, l_min=None,
                          sliver_q=0.2, max_ops=24, park=None):
    """Wire a :class:`DcmRemeshUpdater` onto a ``build_gpu_dcm_simulation`` handle dict.

    ``period`` is the LOW remesh cadence (steps between passes). ``l_min`` defaults to
    0.5·mean_edge (the start mesh mid-band). Returns (action, CustomUpdater).
    """
    sim = handles["sim"]
    turgor = handles["turgor"]
    contact = handles.get("contact")
    # only a FaceContactForceGPU carries remeshable face arrays; a node-node tent does
    # not (pass None → turgor-only remesh).
    if contact is not None and not hasattr(contact, "faces"):
        contact = None
    if l_min is None:
        l_min = 0.5 * float(handles["mean_edge"])
    action = DcmRemeshUpdater(
        turgor=turgor, contact=contact, cell_of_node=handles["cell_of_node"],
        l_min=l_min, sliver_q=sliver_q, max_ops=max_ops, park=park)
    cu = hoomd.update.CustomUpdater(
        trigger=hoomd.trigger.Periodic(int(period)), action=action)
    sim.operations.updaters.append(cu)
    return action, cu
