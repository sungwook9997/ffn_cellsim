"""Host-side 3-zone necrosis (depth-from-surface) for the Warp DCM loop (C8).

Ports cell/dcm_active.SpheroidNecrosis: a spheroid's diffusion-limited O2/glucose supply means
cells deep in the core starve. Modelled by depth from the cluster surface (the standard
depth-proxy; depth>150µm matches the ~100-200µm O2 penetration length — a continuous
Schaller-Meyer-Hermann field is the future oracle). Each tick assigns every live cell a zone:

    d = R_cluster − r_cell          (R_cluster = max active centroid radius + R_cell)
    d < d_prolif            → PROLIFERATING (rim; the only cells allowed to divide)
    d_prolif ≤ d < d_necrotic → QUIESCENT (alive, no division)
    d ≥ d_necrotic           → NECROTIC (irreversible; loses turgor regulation, no division)

Outputs per cell: a zone code, a turgor multiplier (necrotic core softens → ``turgor_necrotic``),
and a ``can_divide`` mask (proliferating only) the DivisionHost respects. Inert at small N (no
cell is >150µm deep until the spheroid is large), which is physically correct.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PROLIFERATING, QUIESCENT, NECROTIC = 0, 1, 2


@dataclass
class NecrosisParams:
    d_prolif_um: float = 40.0
    d_necrotic_um: float = 150.0
    turgor_necrotic: float = 0.3       # necrotic-core turgor multiplier (lost pressure regulation)
    batch_steps: int = 2000


class NecrosisHost:
    """Assigns the 3 depth zones and the per-cell turgor multiplier + division mask each tick."""

    def __init__(self, *, cof: np.ndarray, n_cells: int, npc: int, R: float,
                 params: NecrosisParams | None = None):
        self.p = params or NecrosisParams()
        self.cof = np.asarray(cof, dtype=np.int64)
        self.n_cells = int(n_cells)
        self.npc = int(npc)
        self.R = float(R)
        self.batch_steps = self.p.batch_steps
        self.zone = np.full(n_cells, QUIESCENT, dtype=np.int32)
        self.turgor_mult = np.ones(n_cells, dtype=np.float64)
        self.can_divide = np.zeros(n_cells, dtype=bool)
        self.depth_um = np.zeros(n_cells, dtype=np.float64)

    def update(self, P: np.ndarray, cof: np.ndarray) -> None:
        """Recompute zones from the current geometry. Necrosis is IRREVERSIBLE (a cell already
        NECROTIC stays necrotic even if the surface later recedes)."""
        P = np.asarray(P, dtype=np.float64)
        active = cof[np.arange(self.n_cells) * self.npc] >= 0
        ids = np.flatnonzero(active)
        if ids.size == 0:
            return
        cents = np.array([P[c * self.npc:(c + 1) * self.npc].mean(0) for c in ids])
        sph = cents.mean(0)
        r = np.linalg.norm(cents - sph, axis=1)
        R_cluster = float(r.max()) + self.R
        depth_um = (R_cluster - r) * 1e6
        for k, c in enumerate(ids):
            self.depth_um[c] = depth_um[k]
            if self.zone[c] == NECROTIC:
                continue                                   # irreversible
            if depth_um[k] >= self.p.d_necrotic_um:
                self.zone[c] = NECROTIC
                self.turgor_mult[c] = self.p.turgor_necrotic
            elif depth_um[k] < self.p.d_prolif_um:
                self.zone[c] = PROLIFERATING
            else:
                self.zone[c] = QUIESCENT
        self.can_divide = active & (self.zone == PROLIFERATING)

    def counts(self):
        return {"prolif": int((self.zone == PROLIFERATING).sum()),
                "quiescent": int((self.zone == QUIESCENT).sum()),
                "necrotic": int((self.zone == NECROTIC).sum())}
