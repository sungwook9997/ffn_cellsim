"""Host-side cell division / proliferation for the Warp DCM loop (C7).

MITOTIC-ROUNDING, VOLUME-CONSERVING cytokinesis (PI 2026-06-24). The earlier scheme activated a
FRESH FULL-SIZE, full-rest-volume daughter icosphere 2.4·R OUTWARD of the mother and dropped it cold
into the dense pack — so (1) the daughter overlapped neighbour cells, (2) turgor saw it at full rest
volume and immediately inflated it, and (3) every division ADDED a whole cell-volume (non-conserving).
The net was an outward "shell" of daughters that punched through the ULA wall and inflated A/A0 / maxZ
(a placement artifact, NOT real spreading — the original core stayed compact). Diagnosis + design:
``docs/v2_audit`` division-mesh-collision workflow.

This replaces that with real cytokinesis (the SimuCell3D / biology pattern, minus the in-place mesh
cleave which needs the remesh co-run that is not yet landed):

  * a dividing cell ROUNDS UP and SPLITS into TWO daughters that TOGETHER occupy the mother volume —
    each is a HALF-VOLUME icosphere (radius R·½^⅓) placed at the mother centroid, offset ±½·sep along
    the Hertwig long axis (PCA largest eigenvector of the mother's node cloud). No outward dump.
  * each daughter's per-cell rest volume ``V0_cell`` is set to V0/2, so turgor sees ~zero pressure
    mismatch at birth (no inflation kick) — VOLUME-CONSERVING at the division instant.
  * a time-consistent REGROWTH ramp (:meth:`grow`) re-inflates each daughter's ``V0_cell`` back to the
    full V0 over one cell cycle; turgor chases that ramping setpoint, so the daughters grow GRADUALLY
    (quasi-statically) and push neighbours apart instead of a cold full-size insert.

Faces/edges are pre-allocated at build (a parked-cell pool, nodes cof=−1), so division still touches
only ``pos`` + ``cof`` + ``V0_cell`` — no device topology resync. Mirrors the rim+ConvexHull gating of
``cell/dcm_active.ProliferationUpdater``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class DivisionParams:
    # per-rim-cell division probability per tick. Anchored to the HOOMD reference
    # (cell/dcm_active.py:157 = 0.04): cell-cycle ~hours ≫ spreading ~minutes ⇒ 0-2 divisions per
    # spread. The reference WARNS a high p_div makes A/A0 division-dominated (~9), masking the
    # ~2-3 traction band — so do NOT raise this for a spreading/de-cohesion verdict (review fix #2).
    p_div: float = 0.04
    # cytokinesis furrow separation between the two daughters, as a multiple of the HALF-cell radius
    # (R·½^⅓). A GEOMETRIC quantity (not tuned to an outcome): 0 = co-centred (max sister overlap,
    # most compact vs neighbours), 2 = just-touching (no sister overlap, longest dumbbell). 1.5 keeps
    # neighbour disturbance small while leaving only a shallow sister overlap the contact resolves.
    div_sep_factor: float = 1.5
    batch_steps: int = 2000        # division updater cadence
    seed: int = 23


class DivisionHost:
    """Owns the active/parked cell mask and proliferates rim cells by mitotic-rounding split.

    Usage (division mode):
        div = DivisionHost(verts0=icosphere_verts, npc=npc, n_total=n_total, n_active0=n_active,
                           R=R, z0=z0, V0_full=V0)
        if s % div.batch_steps == 0:
            if div.update(P, cof, V0_cell, can_divide=can_div):   # mutates P, cof, V0_cell in place
                pos_d.assign(P); cof_d.assign(cof.int32); V0_cell_d.assign(V0_cell)
                lam.cof=cof; cad.cof=cof; ...
            div.grow(V0_cell, dV0)                                # regrowth ramp (every batch)
            V0_cell_d.assign(V0_cell)
    """

    def __init__(self, *, verts0: np.ndarray, npc: int, n_total: int, n_active0: int,
                 R: float, z0: float, V0_full: float | None = None,
                 params: DivisionParams | None = None):
        self.p = params or DivisionParams()
        self.verts0 = np.asarray(verts0, dtype=np.float64)   # template icosphere (centred at 0, radius R)
        self.npc = int(npc)
        self.n_total = int(n_total)
        self.R = float(R)
        self.z0 = float(z0)
        # full rest volume (per-cell osmotic setpoint target). Default from the template radius.
        self.V0_full = float(V0_full) if V0_full is not None \
            else (4.0 / 3.0) * np.pi * (float(np.linalg.norm(self.verts0, axis=1).mean())) ** 3
        self.half_scale = 0.5 ** (1.0 / 3.0)                 # half-volume → radius scale (½^⅓ ≈ 0.7937)
        self.batch_steps = self.p.batch_steps
        self._rng = np.random.default_rng(self.p.seed)
        self.n_active0 = int(n_active0)
        self.n_divisions = 0

    def _cell_active(self, cof: np.ndarray) -> np.ndarray:
        """(n_total,) bool: a cell is active iff its first node is live (cof ≥ 0)."""
        first = cof[np.arange(self.n_total) * self.npc]
        return first >= 0

    def _long_axis(self, Pm: np.ndarray) -> np.ndarray:
        """Hertwig long axis = unit eigenvector of the largest eigenvalue of the node covariance."""
        d = Pm - Pm.mean(0)
        cov = d.T @ d
        w, V = np.linalg.eigh(cov)
        n = V[:, int(np.argmax(w))]
        nn = float(np.linalg.norm(n))
        if nn < 1e-12 or not np.isfinite(nn):
            n = self._rng.standard_normal(3); nn = float(np.linalg.norm(n))
        return n / nn

    def update(self, P: np.ndarray, cof: np.ndarray, V0_cell: np.ndarray,
               can_divide: np.ndarray | None = None,
               force_cell: int | None = None) -> bool:
        """One division tick — MITOTIC-ROUNDING volume-conserving split.

        Mutates ``P`` (daughter + reset mother node positions), ``cof`` (activates the daughter), and
        ``V0_cell`` (halves both daughters' rest volume) IN PLACE; returns True if any cell divided.
        ``can_divide`` (per-cell bool, C8 necrosis) gates to proliferating cells only. ``force_cell``
        (test hook) forces exactly that cell id to divide once, bypassing the rim+probability gate."""
        active = self._cell_active(cof)
        active_ids = np.flatnonzero(active)
        if not np.any(~active):
            return False                                  # no free parked cell

        if force_cell is not None:                        # deterministic test path: one forced split
            candidates = [int(force_cell)] if active[int(force_cell)] else []
        else:
            if active_ids.size < 4:
                return False                              # need a hull to find the rim
            cents = np.array([P[c * self.npc:(c + 1) * self.npc].mean(0) for c in active_ids])
            rim = np.zeros(active_ids.size, dtype=bool)
            try:
                from scipy.spatial import ConvexHull
                # QJ joggle so a near-coplanar cluster doesn't QhullError and silently disable division.
                rim[np.unique(ConvexHull(cents, qhull_options="QJ").vertices)] = True
            except Exception:
                rim[:] = True
            candidates = []
            for k in np.where(rim)[0]:
                cid = int(active_ids[k])
                if can_divide is not None and not can_divide[cid]:
                    continue                              # C8: only proliferating-zone cells divide
                if self._rng.random() >= self.p.p_div:
                    continue
                candidates.append(cid)

        divided = False
        for mother in candidates:
            free = np.flatnonzero(~active)
            if free.size == 0:
                break
            daughter = int(free[0])
            mlo, mhi = mother * self.npc, (mother + 1) * self.npc
            Pm = P[mlo:mhi]
            C = Pm.mean(0)
            n = self._long_axis(Pm)
            r_half = self.R * self.half_scale
            half = 0.5 * self.p.div_sep_factor * r_half               # ±offset of each daughter centre
            tmpl = self.verts0 * self.half_scale                      # half-VOLUME icosphere (radius r_half)
            c_mom = C - half * n
            c_dau = C + half * n
            c_mom[2] = max(c_mom[2], self.z0 + r_half)                # rest on / above the dish
            c_dau[2] = max(c_dau[2], self.z0 + r_half)
            # round-up + split: mother is RESET to a half-volume sphere at c_mom; daughter activated at c_dau
            P[mlo:mhi] = tmpl + c_mom
            dlo, dhi = daughter * self.npc, (daughter + 1) * self.npc
            P[dlo:dhi] = tmpl + c_dau
            cof[dlo:dhi] = daughter
            V0_cell[mother] = 0.5 * self.V0_full                      # volume-conserving: each daughter V0/2
            V0_cell[daughter] = 0.5 * self.V0_full
            active[daughter] = True
            self.n_divisions += 1
            divided = True
        return divided

    def grow(self, V0_cell: np.ndarray, dV0: float, cof: np.ndarray | None = None) -> None:
        """Time-consistent post-mitotic REGROWTH: ramp every ACTIVE cell's rest volume toward the full
        V0 by ``dV0`` (capped at V0). Cells already at V0 are unchanged; freshly born half-cells climb
        back to V0 over the cell cycle (``dV0`` is sized by the caller from the cell-cycle time + the
        time-acceleration S, so it is NOT a free knob). Turgor chases this setpoint → gradual inflation."""
        if dV0 <= 0.0:
            return
        if cof is not None:
            active = self._cell_active(cof)
        else:
            active = np.ones(self.n_total, dtype=bool)
        grow_mask = active & (V0_cell < self.V0_full)
        V0_cell[grow_mask] = np.minimum(V0_cell[grow_mask] + dV0, self.V0_full)
