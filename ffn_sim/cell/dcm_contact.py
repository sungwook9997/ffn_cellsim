"""DCM Upgrade B — SimuCell3D bilinear-tent cell-cell contact.

Replaces the crude node-LJ / capped-harmonic ``CellCellAdhesion`` (which
over-dispersed cells → A/A0 blew up to ~17) with the SimuCell3D node-vs-surface
penalty contact (CONTACT_MODEL_INDEX = 0), ported as a ``md.force.Custom`` Action
on the BAOAB path.

Force law (Lead-verified re-derivation of the C++ ``contact_node_face_via_spring``;
see ``docs/v2_audit/SIMUCELL3D_INTEGRATION_2026-06-11.md`` §2.2 step 4):

  Let ``r_vec`` point from the other cell's surface point to this node, and
  ``d = |r_vec|`` the separation. The code's scalar ``force_amp ∝ (c/d − 1)``
  LOOKS like 1/d, but the actual force ``F = r_vec · force_amp`` carries the
  UNnormalised ``r_vec`` whose magnitude ``d`` cancels the 1/d → a clean TENT:

    * ADHESION (cohesion), d ∈ [0, c_adh]:  symmetric bilinear tent, PEAK at
      d = c_adh/2, zero at d = 0 and d = c_adh:
          hardening (c/2 ≤ d < c):  |F| = ω · A · (c_adh − d)
          softening (0 ≤ d < c/2):  |F| = ω · A · d
      Direction: attractive (pulls the two cells together).
    * REPULSION (non-penetration), overlap d_ov > 0 below the contact radius:
          |F| = ξ · A · d_ov     (linear penalty, pushes apart)

  ω = ``adh_strength`` (Pa/m), ξ = ``rep_strength`` (Pa/m), A = contact patch
  area (here a per-node membrane patch; the true model uses the face area).

GEOMETRY USED HERE (stated honestly, per the brief's "state which you used"):
  **node ↔ other-cell NODE** neighbour search (scipy cKDTree), NOT the full
  closest-point-on-triangle (Ericson) node↔face. This is the explicitly-permitted
  tractable first cut: the SAME bilinear tent law is applied with ``r_vec`` the
  inter-node vector and a per-node patch area. The repulsion branch uses a contact
  radius ``r_contact`` (≈ node spacing) as the effective membrane thickness so two
  shells cannot interpenetrate through the gaps between nodes for moderate overlap.
  (The node↔face closest-point kernel is the documented Phase-1 upgrade for exact
  non-penetration; it is not needed to fix the A/A0=17 over-dispersion, which is
  caused by the contact LAW, not the geometry.)

Inter-cell EXCLUSION: only node pairs from DIFFERENT cells interact (uses the
mutable ``cell_of_node`` like ``CellCellAdhesion`` so the proliferation updater
can keep editing it). The ``cad_mult`` junction-switch seam is preserved via the
optional ``cad_mult`` per-cell array (``mult = √(cad_mult_i · cad_mult_j)`` on the
adhesion branch).

Forces are CAPPED (``force_cap``) so a residual t=0 overlap can never blow BAOAB's
int32 step counter. Newton-3 pairwise (equal-and-opposite).

Units SI.
"""

from __future__ import annotations

import numpy as np

import hoomd
import hoomd.md as md


class DcmTentContact(md.force.Custom):
    """SimuCell3D bilinear-tent cell-cell contact (node↔node first cut).

    Args:
        cell_of_node: (N,) int32, mutable; cell id of each node, -1 = dormant.
        r_contact: m, effective contact radius (≈ node spacing). Overlap is
            measured as ``r_contact − d`` for the repulsion branch.
        c_adh: m, adhesion cutoff (SimuCell3D ``c_adh``); tent peaks at c_adh/2.
        rep_strength: Pa/m, repulsion stiffness ξ.
        adh_strength: Pa/m, adhesion stiffness ω (0 ⇒ repulsion-only).
        patch_area: m², contact patch area A per node (= 4πR²/n_nodes).
        force_cap: N, per-pair force-magnitude cap (BAOAB int32 guard).
        cad_mult: optional (n_cells,) float; per-cell adhesion multiplier for the
            junction switch (``mult = √(cad_mult_i·cad_mult_j)``). None ⇒ all 1.
    """

    def __init__(self, *, cell_of_node: np.ndarray, r_contact: float, c_adh: float,
                 rep_strength: float, adh_strength: float, patch_area: float,
                 force_cap: float = 5.0e-8, cad_mult: np.ndarray | None = None) -> None:
        super().__init__(aniso=False)
        self.cell_of_node = cell_of_node            # (N,) int, -1 = dormant
        self.r_contact = float(r_contact)
        self.c_adh = float(c_adh)
        self.rep = float(rep_strength)
        self.omega = float(adh_strength)
        self.A = float(patch_area)
        self.force_cap = float(force_cap)
        self.cad_mult = cad_mult                    # (n_cells,) or None
        # Search radius: the larger of the contact (repulsion) and adhesion ranges.
        self.r_search = max(self.r_contact, self.c_adh)

    def set_forces(self, timestep: int) -> None:  # noqa: D401
        with self._state.cpu_local_snapshot as snap:
            tag = np.asarray(snap.particles.tag).copy()
            pos = np.asarray(snap.particles.position).copy()
        n = pos.shape[0]
        perm = np.argsort(tag)              # local row → global (tag) order
        pos_g = pos[perm]
        cell = self.cell_of_node           # global-idx order

        F_g = np.zeros_like(pos_g)
        active = np.where(cell >= 0)[0]
        if active.size > 1:
            apos = pos_g[active]
            acell = cell[active]
            pairs = _within_cutoff_pairs(apos, self.r_search)
            if pairs.size:
                ii, jj = pairs[:, 0], pairs[:, 1]
                diff = acell[ii] != acell[jj]      # inter-cell only
                ii, jj = ii[diff], jj[diff]
                if ii.size:
                    # r_vec points from node jj (other cell surface) to node ii.
                    rvec = apos[ii] - apos[jj]
                    d = np.linalg.norm(rvec, axis=1)
                    d_safe = np.where(d > 1e-18, d, 1e-18)
                    rhat = rvec / d_safe[:, None]

                    fmag = np.zeros_like(d)         # +out (repulsion), −in (adhesion)

                    # REPULSION: overlap below contact radius (non-penetration).
                    ov = self.r_contact - d
                    rep_m = ov > 0.0
                    fmag[rep_m] = self.rep * self.A * ov[rep_m]   # + = push apart

                    # ADHESION: bilinear tent over [0, c_adh] on the NON-overlapping
                    # pairs (d ≥ r_contact and d < c_adh).
                    if self.omega > 0.0:
                        adh_m = (~rep_m) & (d < self.c_adh)
                        da = d[adh_m]
                        half = 0.5 * self.c_adh
                        tent = np.where(da >= half,
                                        self.c_adh - da,    # hardening branch
                                        da)                 # softening branch
                        amag = self.omega * self.A * tent   # ≥ 0 magnitude
                        if self.cad_mult is not None:
                            mi = self.cad_mult[acell[ii][adh_m]]
                            mj = self.cad_mult[acell[jj][adh_m]]
                            amag = amag * np.sqrt(mi * mj)
                        fmag[adh_m] = -amag             # − = pull together

                    # Cap the per-pair force magnitude (BAOAB int32 guard).
                    np.clip(fmag, -self.force_cap, self.force_cap, out=fmag)

                    fvec = fmag[:, None] * rhat        # on node ii (+ = away from jj)
                    gi = active[ii]
                    gj = active[jj]
                    np.add.at(F_g, gi, fvec)
                    np.add.at(F_g, gj, -fvec)          # Newton-3

        F = np.empty_like(pos)
        F[perm] = F_g
        with self.cpu_local_force_arrays as arr:
            arr.force[:] = F
            arr.potential_energy[:] = np.zeros(n, dtype=np.float64)


def _within_cutoff_pairs(pos: np.ndarray, r_cut: float) -> np.ndarray:
    """(P,2) int i<j pairs with |pos_i − pos_j| < r_cut via scipy cKDTree."""
    nn = pos.shape[0]
    if nn < 2:
        return np.empty((0, 2), dtype=np.int64)
    from scipy.spatial import cKDTree
    tree = cKDTree(pos)
    pairs = tree.query_pairs(r=r_cut, output_type="ndarray")
    if pairs.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    return pairs.astype(np.int64)
