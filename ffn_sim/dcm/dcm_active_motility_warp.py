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


@wp.kernel
def active_drift_translate_kernel(
    cof: wp.array(dtype=wp.int32),
    pol: wp.array(dtype=wp.vec3d),        # (n_cells,) per-cell unit polarity
    dx: wp.float64,                       # v0 · dt_step: the per-step drift distance (v0-mode)
    pos: wp.array(dtype=wp.vec3d),        # (N,) in/out — translated in place
):
    """OPERATOR-SPLIT active drift (timescale attack, 2026-07-11). Translate each live node by dx·p̂_c BEFORE
    the implicit elastic/contact relax, instead of adding a self-propulsion FORCE that the backward-Euler solve
    over-damps. Lie-Trotter split of the overdamped dynamics γẋ = F_active + F_passive: substep-1 (drift) is
    exact — x* = xₙ + (F_active/γ)·dt = xₙ + v0·dt·p̂ — so the SPV drift is preserved at ANY dt (the force-based
    application froze at large dt because the implicit solve equilibrates F_active against contact each step,
    skipping the non-equilibrium T1 creep). substep-2 (relax) is the existing implicit IPC-Newton solve, whose
    backward-Euler anchor xₙ is captured AFTER this translation → relaxes from x*. dx tiny (v0·dt) so the
    per-step overlap IPC must resolve is small at every dt."""
    i = wp.tid()
    c = cof[i]
    if c < wp.int32(0):
        return
    pos[i] = pos[i] + pol[c] * dx


class ActiveMotilityHost:
    """Per-cell polarity with rotational diffusion (persistent random walk). Host-updates p̂_c each batch;
    the device kernel applies F_active·p̂_c per cell every step."""

    def __init__(self, *, n_cells: int, cof: np.ndarray, f_active_N: float = 0.0,
                 persistence_s: float = 600.0, planar: bool = True, seed: int = 7,
                 v0_um_s: float = 0.0, gamma_node: float = 0.0):
        """Two ways to set the self-propulsion (the SPV v0 axis):
          - v0_um_s > 0 (PREFERRED, faithful to SPV): target a physiological migration SPEED v0 [µm/s]; the
            per-node force = gamma_node·v0 so each cell self-propels at v0 (overdamped v=F/γ). Bounds the
            velocity → bounds CFL → maps the phase diagram cleanly at any v0. gamma_node = the physiological
            per-node Stokes drag (6π·η·R/npc) from the caller. (physiological cell speed ~0.5-2 µm/min.)
          - f_active_N > 0 (legacy): a fixed whole-cell propulsion FORCE [N] (mesh-invariant per node); CFL
            grows with F when the drag gives super-physiological speed.
        persistence_s = polarity reorientation time τ_p (~10 min lit). planar = crawl in-plane (on substrate)."""
        self.n_cells = int(n_cells)
        self.f_active_N = float(f_active_N)
        self.v0 = float(v0_um_s) * 1e-6            # µm/s → m/s
        self.gamma_node = float(gamma_node)
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
        if self.v0 > 0.0 and self.gamma_node > 0.0:
            # target speed v0: per-node force = gamma_node·v0 → each node (hence the cell) moves at v0
            fnode = np.full(self.n_cells, self.gamma_node * self.v0, dtype=np.float64)
        else:
            fnode = (self.f_active_N / self._Nc).astype(np.float64)
        if self._pol_d is None:
            self._pol_d = wp.zeros(self.n_cells, dtype=wp.vec3d, device=device)
            self._fnode_d = wp.array(fnode, dtype=wp.float64, device=device)
        self._pol_d.assign(self.pol.astype(np.float64))
        return self._pol_d, self._fnode_d
