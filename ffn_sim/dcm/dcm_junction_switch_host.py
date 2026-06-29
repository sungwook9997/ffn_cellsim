"""Host-side bulk-pressure junction switch (cadherin → integrin clutch) for the Warp loop.

Phase C migration · M3 — ports the de-cohesion junction switch off the HOOMD
``GpuJunctionSwitchUpdater`` / ``attach_junction_switch`` (``cell/dcm_gpu_build.py`` on
branch ``dcm/decohesion``) onto the GPU-resident Warp loop.

Mechanism (faithful to the frozen proxy): at low cadence the host computes each live
cell's bulk pressure from neighbour CROWDING — count cells whose centroid is within
``contact_factor·R``, then clip-ramp that crowd count to the literature ``[P_min, P_max]``
kPa band between ``crowd_lo`` and ``crowd_hi``. Any cell whose pressure exceeds
``P_switch_kPa`` is LATCHED switched: its cadherin multiplier ``cad_mult`` drops to
``cadherin_weak_factor`` (cell-cell adhesion weakens → de-cohesion) and its substrate
``integrin_gain`` rises to ``integrin_strong_factor`` (cell-substrate traction strengthens).
The device kernels read these per-cell arrays live: ``cohesion_grid_cad_kernel`` /
``contact_grid_cad_kernel`` scale cell-cell adhesion by ``sqrt(cad_i·cad_j)``, and
``dcm_wetting_scatter_integrin_kernel`` scales the per-cell substrate wetting by ``integrin``.

Additive / default-off: with the switch absent both arrays stay all-1 and the build is
unchanged. Per the landmine register this is a SECONDARY (settle-dominated) effect.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class JunctionParams:
    """Frozen literature thresholds (``ResolvedSpheroidState`` defaults, dcm/decohesion)."""
    contact_factor: float = 2.2          # neighbour contact radius = factor·R
    crowd_lo: float = 4.0                # crowd ≤ this → P_min
    crowd_hi: float = 12.0               # crowd ≥ this → P_max (full FCC shell)
    P_min_kPa: float = 0.0
    P_max_kPa: float = 6.0
    P_switch_kPa: float = 0.5            # bulk-pressure onset for the switch
    cadherin_weak_factor: float = 0.3    # cell-cell multiplier when switched
    integrin_strong_factor: float = 3.0  # cell-substrate multiplier when switched
    cadence: int = 500


class JunctionSwitchHost:
    """Owns the mutable per-cell ``cad_mult`` / ``integrin_gain`` arrays and latches the
    crowd-pressure junction switch at low cadence.

    Usage in the driver loop:
        js = JunctionSwitchHost(cof=cof_a, n_cells=..., R=...)
        ...
        if s % js.cadence == 0:
            if js.update(pos_d.numpy()):       # returns True if any new cell switched
                cad_d.assign(js.cad_mult); integrin_d.assign(js.integrin_gain)
    """

    def __init__(self, *, cof: np.ndarray, n_cells: int, R: float,
                 params: JunctionParams | None = None):
        self.p = params or JunctionParams()
        self.cof = np.asarray(cof, dtype=np.int64)
        self.n_cells = int(n_cells)
        self.R = float(R)
        self.r_contact = self.p.contact_factor * self.R
        self.cadence = self.p.cadence
        self.cad_mult = np.ones(n_cells, dtype=np.float64)
        self.integrin_gain = np.ones(n_cells, dtype=np.float64)
        self.switched = np.zeros(n_cells, dtype=bool)
        self.pressure_kPa = np.zeros(n_cells, dtype=np.float64)
        self.n_switched = 0

    def _cell_centroids(self, P: np.ndarray) -> np.ndarray:
        c = np.zeros((self.n_cells, 3), dtype=np.float64)
        cnt = np.zeros(self.n_cells, dtype=np.float64)
        valid = self.cof >= 0
        np.add.at(c, self.cof[valid], P[valid])
        np.add.at(cnt, self.cof[valid], 1.0)
        return c / np.maximum(cnt, 1.0)[:, None]

    def update(self, P: np.ndarray) -> bool:
        """One switch tick: compute per-cell crowd pressure and latch cells above
        ``P_switch_kPa``. Returns True if any NEW cell switched this tick (so the
        caller re-uploads the device arrays)."""
        P = np.asarray(P, dtype=np.float64)
        cents = self._cell_centroids(P)
        # crowd = number of OTHER cell centroids within r_contact (frozen proxy)
        d2 = np.sum((cents[:, None, :] - cents[None, :, :]) ** 2, axis=2)
        within = d2 < self.r_contact ** 2
        np.fill_diagonal(within, False)
        crowd = within.sum(axis=1).astype(np.float64)
        frac = np.clip((crowd - self.p.crowd_lo) / (self.p.crowd_hi - self.p.crowd_lo),
                       0.0, 1.0)
        self.pressure_kPa = self.p.P_min_kPa + (self.p.P_max_kPa - self.p.P_min_kPa) * frac
        newly = (~self.switched) & (self.pressure_kPa > self.p.P_switch_kPa)
        if not newly.any():
            self.n_switched = int(self.switched.sum())
            return False
        self.switched |= newly
        self.cad_mult[newly] = self.p.cadherin_weak_factor
        self.integrin_gain[newly] = self.p.integrin_strong_factor
        self.n_switched = int(self.switched.sum())
        return True
