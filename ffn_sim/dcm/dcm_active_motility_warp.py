"""Per-cell active self-propulsion (cell motility) — the active-matter unjamming lever.

DIAGNOSIS (DCM_SPREADING_INVESTIGATION_PLAN_2026-07-10 §6d): the DCM aggregate is JAMMED — cells do not
slide past each other, so neither rate-gated compaction nor spreading/dispersal works, regardless of
cadherin. The physical lever to UNJAM a dense cell packing is ACTIVE MOTILITY: each cell self-propels in a
slowly-reorienting direction (persistent random walk — the standard active-matter model of migrating cells,
Szabo 2006 / Henkes 2011). Active self-propulsion above a threshold FLUIDISES a jammed tissue → cells
rearrange (T1) → the aggregate can flow (compact-rate-gated) or disperse (spread) if de-cohesion allows.

This module adds a per-cell self-propulsion force F_active·p̂_c (p̂_c = the cell's polarity, host-updated by
rotational diffusion). It is applied per node so the whole-cell propulsion = F_active·p̂_c (node force =
F_active·p̂_c / N_nodes_c). Default OFF; a controlled-variable magnitude (physiological single-cell traction
~10-100 nN, KB-2.12), never tuned to a spreading target. This is a COARSE motility proxy (the fine-grained
driver is the lamellipodium+clutch, which is too weak) used to TEST whether motility is the unjamming lever;
if it fluidises, invest in the fine-grained version.
"""
from __future__ import annotations

import numpy as np
import warp as wp


@wp.kernel
def active_self_propulsion_kernel(
    cof: wp.array(dtype=wp.int32),
    pol: wp.array(dtype=wp.vec3d),        # (n_cells,) per-cell unit polarity
    fnode: wp.array(dtype=wp.float64),    # (n_cells,) per-node propulsion magnitude = F_active / N_nodes_c
    force: wp.array(dtype=wp.vec3d),      # (N,) out
):
    """Each live node adds its cell's self-propulsion F_active·p̂_c / N_nodes_c → whole-cell force F_active·p̂_c."""
    i = wp.tid()
    c = cof[i]
    if c < wp.int32(0):
        return
    force[i] = force[i] + pol[c] * fnode[c]


class ActiveMotilityHost:
    """Per-cell polarity with rotational diffusion (persistent random walk). Host-updates p̂_c each batch;
    the device kernel applies F_active·p̂_c per cell every step."""

    def __init__(self, *, n_cells: int, cof: np.ndarray, f_active_N: float,
                 persistence_s: float = 600.0, planar: bool = True, seed: int = 7):
        """f_active_N = whole-cell propulsion force [N] (physiological single-cell traction ~1e-8..1e-7 N).
        persistence_s = polarity reorientation time τ_p (cell directional persistence, ~10 min lit).
        planar = keep polarity in-plane (xy) so cells crawl ON the substrate (default), else 3D."""
        self.n_cells = int(n_cells)
        self.f_active_N = float(f_active_N)
        self.tau_p = float(persistence_s)
        self.planar = bool(planar)
        self._rng = np.random.default_rng(seed)
        cof = np.asarray(cof)
        # per-cell node count → per-node force so the WHOLE-CELL propulsion is F_active (mesh-invariant)
        self._Nc = np.array([max(1, int((cof == c).sum())) for c in range(self.n_cells)], dtype=np.float64)
        # initial random unit polarity
        self.pol = self._rand_dirs(self.n_cells)
        self._pol_d = None
        self._fnode_d = None

    def _rand_dirs(self, n):
        v = self._rng.normal(size=(n, 3))
        if self.planar:
            v[:, 2] = 0.0
        nrm = np.linalg.norm(v, axis=1, keepdims=True)
        nrm[nrm < 1e-12] = 1.0
        return v / nrm

    def update(self, dt_batch: float):
        """Rotational diffusion: p̂ ← normalize(p̂ + sqrt(2·D_r·dt)·ξ), D_r = 1/τ_p (persistent random walk)."""
        if self.tau_p <= 0.0:
            return
        sigma = np.sqrt(2.0 * (1.0 / self.tau_p) * max(dt_batch, 0.0))
        xi = self._rng.normal(size=(self.n_cells, 3))
        if self.planar:
            xi[:, 2] = 0.0
        p = self.pol + sigma * xi
        if self.planar:
            p[:, 2] = 0.0
        nrm = np.linalg.norm(p, axis=1, keepdims=True)
        nrm[nrm < 1e-12] = 1.0
        self.pol = p / nrm

    def upload(self, device):
        import warp as wp
        fnode = (self.f_active_N / self._Nc).astype(np.float64)
        if self._pol_d is None:
            self._pol_d = wp.zeros(self.n_cells, dtype=wp.vec3d, device=device)
            self._fnode_d = wp.array(fnode, dtype=wp.float64, device=device)
        self._pol_d.assign(self.pol.astype(np.float64))
        return self._pol_d, self._fnode_d
