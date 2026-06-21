"""Host-side cell division / proliferation for the Warp DCM loop (C7).

Ports the H.9/dcm_active proliferation: a pre-allocated CELL POOL (parked undeformed icospheres,
nodes cof=−1 so cohesion/contact/measure skip them and turgor sees V=V0 → force-free) from which
a DIVISION activates a daughter. At low cadence each RIM cell (on the active-cluster convex hull)
divides with probability ``p_div``: a parked cell is activated — its node block's cof is set to the
daughter id and its undeformed icosphere is placed one cell-diameter OUTWARD of the mother (z floored
to rest on the dish). Faces/edges are pre-allocated (build), so division touches only pos + cof — no
device topology resync. Mirrors cell/dcm_active.ProliferationUpdater (rim+ConvexHull, gapped daughter).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class DivisionParams:
    p_div: float = 0.5             # per-rim-cell division probability per tick (dimensionless)
    div_gap: float = 0.4           # daughter placed at (2 + div_gap)·R outward of the mother
    batch_steps: int = 2000        # division updater cadence
    seed: int = 23


class DivisionHost:
    """Owns the active/parked cell mask and proliferates rim cells into parked daughters.

    Usage (division mode):
        div = DivisionHost(verts0=icosphere_verts, npc=npc, n_total=n_total, n_active0=n_active,
                           R=R, z0=z0)
        if s % div.batch_steps == 0:
            if div.update(P, cof):          # mutates P (daughter pos) + cof in place
                pos_d.assign(P); cof_d.assign(cof.int32); lam.cof=cof; cad.cof=cof; ...
    """

    def __init__(self, *, verts0: np.ndarray, npc: int, n_total: int, n_active0: int,
                 R: float, z0: float, params: DivisionParams | None = None):
        self.p = params or DivisionParams()
        self.verts0 = np.asarray(verts0, dtype=np.float64)   # template icosphere (centred at 0)
        self.npc = int(npc)
        self.n_total = int(n_total)
        self.R = float(R)
        self.z0 = float(z0)
        self.batch_steps = self.p.batch_steps
        self._rng = np.random.default_rng(self.p.seed)
        self.n_active0 = int(n_active0)
        self.n_divisions = 0

    def _cell_active(self, cof: np.ndarray) -> np.ndarray:
        """(n_total,) bool: a cell is active iff its first node is live (cof ≥ 0)."""
        first = cof[np.arange(self.n_total) * self.npc]
        return first >= 0

    def update(self, P: np.ndarray, cof: np.ndarray, can_divide: np.ndarray | None = None) -> bool:
        """One division tick. Mutates P and cof IN PLACE; returns True if any cell divided.
        ``can_divide`` (per-cell bool, from C8 necrosis) further gates which cells may divide —
        only PROLIFERATING (rim, nutrient-supplied) cells; necrotic/quiescent never divide."""
        active = self._cell_active(cof)
        active_ids = np.flatnonzero(active)
        if active_ids.size < 4 or not np.any(~active):
            return False                                  # need a hull + a free parked cell
        cents = np.array([P[c * self.npc:(c + 1) * self.npc].mean(0) for c in active_ids])
        cluster_cen = cents.mean(0)
        rim = np.zeros(active_ids.size, dtype=bool)
        try:
            from scipy.spatial import ConvexHull
            rim[np.unique(ConvexHull(cents).vertices)] = True
        except Exception:
            return False
        divided = False
        for k in np.where(rim)[0]:
            free = np.flatnonzero(~active)
            if free.size == 0:
                break
            if can_divide is not None and not can_divide[int(active_ids[k])]:
                continue                                  # C8: only proliferating-zone cells divide
            if self._rng.random() >= self.p.p_div:
                continue
            daughter = int(free[0])
            outward = cents[k] - cluster_cen
            nrm = float(np.linalg.norm(outward))
            if nrm < 1e-12:
                outward = self._rng.standard_normal(3); nrm = float(np.linalg.norm(outward))
            outward = outward / nrm
            new_center = cents[k] + (2.0 + self.p.div_gap) * self.R * outward
            new_center[2] = max(new_center[2], self.z0 + self.R)   # rest on / above the dish
            lo, hi = daughter * self.npc, (daughter + 1) * self.npc
            P[lo:hi] = self.verts0 + new_center
            cof[lo:hi] = daughter
            active[daughter] = True
            self.n_divisions += 1
            divided = True
        return divided
