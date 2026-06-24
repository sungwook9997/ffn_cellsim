"""Host-side lamellipodium ratchet + leading-node geometry for the Warp DCM loop.

Phase C migration · M2 — the de-cohesion per-cell lamellipodium, ported off the
HOOMD ``cpu_local_snapshot`` updater (``cell/dcm_lamellipodium_gpu.py`` on branch
``dcm/decohesion``) onto the GPU-resident Warp loop.

Design (M2b §F — the leanest faithful port): **actin = host-managed advancing anchor
points, NOT BAOAB DOFs.** The clutch is ~rigid relative to the tether
(``k_clutch 5e-2 ≫ k_tether 4e-3``, tether capped 5e-9 ⇒ actin slip ≤
tether_cap/k_clutch = 0.1µm = 0.2·ℓ₀), so the gripped actin bead is the fixed point
the membrane is pulled toward. This module owns that anchor pool on the host: it
advances the front at low cadence (``batch_steps``) and rewrites the device anchor
array + per-leading-node geometry arrays that the every-step device tether
(:func:`dcm_neighbor_warp.lamellipodium_tether_multicell`) reads.

What is faithful vs the HOOMD reference:
  * rim-cell detection (``detect_rim_cells``): centroid z within ``contact_band·R`` of z0.
  * two-centroid leading-node predicate (per-cell centroid for proj, spheroid centroid
    for the outward dir; z dropped → 2D radial).
  * ratchet SEED/ADVANCE with ``rng = np.random.default_rng(7)`` and
    ``p_advance = min(1, v_front·batch_steps·dt·S/ℓ₀)``.
The only simplifications (documented, immaterial for a spreading verdict): actin is
pinned (no DOF integration, so no clutch-spring force / no Newton-3 reaction), and
there is no tag permutation (Warp node row i == node i; actin is a separate array).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class LamelParams:
    """Resolved lamellipodium params (defaults = ``ResolvedGpuLamellipodium``, M2b §G)."""
    v_front: float = 6.0e-6 / 60.0      # m/s (6 µm/min; lit 3–12)
    l0: float = 0.5e-6                  # m actin rest length ℓ₀
    S_kinetic: float = 1.0
    catch_factor: float = 1.2           # ×ℓ₀
    seed_offset: float = 1.0            # ×ℓ₀
    rim_contact_band: float = 1.5       # ×R
    lead_frac: float = 0.0
    basal_band_factor: float = 0.3      # ×R
    seed_basal_offset: float = 0.3e-6   # m
    # B2: node-to-plane substrate clutch. Off (legacy) → the advancing actin anchor floats at
    # z0+seed_basal_offset (an abstract actin-bead the membrane is tethered to). On → the anchor
    # sits ON the dish plane (z=z0): the lamellipodial protrusion's nascent adhesion GRIPS a
    # substrate ligand at the advancing front and pulls the leading node toward it — the SAME
    # node-to-plane geometry as the ECM integrin clutch (unified-actin-architecture). The
    # advancing-front + tether spring are unchanged; only the anchor's z is the dish, so the
    # traction is transmitted to the rigid substrate (node-to-plane) instead of a floating bead.
    # (Full Pereverzev catch-slip on the lamellipodial clutch is a documented follow-up.)
    substrate_clutch: bool = False
    k_tether: float = 4.0e-3            # N/m
    tether_cap: float = 5.0e-9          # N (5·60Pa·area_node, Gil-Redondo 2023)
    tether_radius: float = 3.0e-6       # m
    pool_per_cell: int = 60
    batch_steps: int = 50
    seed: int = 7


def detect_rim_cells(centers: np.ndarray, z0: float, R: float,
                     contact_band: float) -> np.ndarray:
    """Rim cells = those whose centroid z is within ``contact_band·R`` of the dish z0
    (``dcm_lamellipodium_gpu.py:272``). A substrate-touching ball's bottom cells crawl;
    periphery-vs-centre then self-selects via the outward-radial dir (a central contact
    cell's outward dir is ill-defined so it grows few leading nodes). Falls back to the
    lowest third by z if none are within the band."""
    n = centers.shape[0]
    basal = np.flatnonzero((centers[:, 2] - z0) <= contact_band * R)
    if basal.size == 0:
        return np.argsort(centers[:, 2])[: max(1, n // 3)].astype(np.int64)
    return basal.astype(np.int64)


class LamellipodiumHost:
    """Owns the shared actin anchor pool and advances the front each tick.

    Usage in the driver loop:
        lam = LamellipodiumHost(pos0=pos_a, cof=cof_a, n_cells=..., z0=..., R=..., dt=...)
        ...
        if s % lam.batch_steps == 0:           # low cadence
            lam.update(pos_d.numpy())           # ratchet + rebuild geometry (host)
            lam.upload(device)                  # push anchors + lead arrays to device
        # every step: gather_lead_pos(pos_d, lam.lead_idx_d) then the tether kernel
    """

    def __init__(self, *, pos0: np.ndarray, cof: np.ndarray, n_cells: int,
                 z0: float, R: float, dt: float, params: LamelParams | None = None,
                 use_gpu_ratchet: bool = False, device: str = "cpu"):
        self.p = params or LamelParams()
        self.cof = np.asarray(cof, dtype=np.int64)
        self.n_cells = int(n_cells)
        self.z0 = float(z0)
        self.R = float(R)
        self.batch_steps = self.p.batch_steps
        # Opt-in GPU geometry path (additive; default OFF → the host path below is
        # byte-identical to before this change). When ON, the deterministic geometry
        # (per-cell centroid, rim set, spheroid centroid, per-leading-node out/proj/rp)
        # is computed by the Warp kernels in :mod:`dcm_lamellipodium_ratchet_warp` reading
        # the LIVE device positions — no per-batch ``pos_d.numpy()`` of all N nodes — and
        # the SAME host ratchet (:meth:`_ratchet_from_geo`) consumes it, so the seed/advance
        # RNG stream and the device anchor arrays are identical to the host path.
        self.use_gpu_ratchet = bool(use_gpu_ratchet)
        self.device = str(device)

        centers = self._cell_centroids(np.asarray(pos0, dtype=np.float64))
        self.rim_cells = detect_rim_cells(centers, z0, R, self.p.rim_contact_band)
        self.n_rim = int(self.rim_cells.size)
        self.n_pool = self.n_rim * self.p.pool_per_cell

        # B2: node-to-plane clutch anchors the protrusion ON the dish (z0); legacy floats it above.
        self.z_basal = z0 if self.p.substrate_clutch else z0 + self.p.seed_basal_offset
        self.basal_band = self.p.basal_band_factor * R
        self.catch = self.p.catch_factor * self.p.l0
        self.seed_off = self.p.seed_offset * self.p.l0
        # band-matched polymerisation probability per tick (computed once)
        self.p_advance = min(1.0, max(0.0,
            self.p.v_front * (self.batch_steps * dt) * self.p.S_kinetic / self.p.l0))
        self._rng = np.random.default_rng(self.p.seed)

        # host anchor pool: position + serving cell (−1 dormant), contiguous high-water mark
        self.actin_xyz = np.zeros((self.n_pool, 3), dtype=np.float64)
        self.actin_cell = np.full(self.n_pool, -1, dtype=np.int64)
        self.n_used = 0
        self.n_seeded = 0
        self.n_advanced = 0

        # per-leading-node geometry (rebuilt each tick) — numpy staging
        self._lead = _empty_lead()
        # device arrays (allocated lazily on first upload)
        self._dev = None

    def _cell_centroids(self, P: np.ndarray) -> np.ndarray:
        c = np.zeros((self.n_cells, 3), dtype=np.float64)
        cnt = np.zeros(self.n_cells, dtype=np.float64)
        valid = self.cof >= 0
        np.add.at(c, self.cof[valid], P[valid])
        np.add.at(cnt, self.cof[valid], 1.0)
        cnt = np.where(cnt > 0, cnt, 1.0)
        return c / cnt[:, None]

    # -- the ratchet + geometry rebuild (HOOMD GpuLamellipodiumAdvance.act + tether geo) --
    def update(self, P: np.ndarray | None = None, *, pos_d=None) -> None:
        """One ratchet tick: SEED/ADVANCE the front, then rebuild leading-node geometry.

        Host path (default): pass ``P`` = current node positions (device → host, (N,3));
        the deterministic geometry is built in numpy here and the ratchet runs on it.

        GPU path (``use_gpu_ratchet=True``): the deterministic geometry is built by the
        Warp kernels in :mod:`dcm_lamellipodium_ratchet_warp` reading the LIVE device
        positions — pass ``pos_d`` (a ``wp.array(vec3d)`` of positions on device) so no
        per-batch ``pos_d.numpy()`` of all N nodes happens; only the O(rim) geometry
        crosses the bus. The SAME host ratchet (:meth:`_ratchet_from_geo`) then consumes
        it, so the RNG seed/advance stream and the device anchor arrays are identical.
        ``P`` may still be passed (and is ignored) for a drop-in call signature."""
        if self.use_gpu_ratchet:
            geo = self._geometry_gpu(P=P, pos_d=pos_d)
        else:
            geo = self._geometry_host(np.asarray(P, dtype=np.float64))
        # ``self.rim_cells`` is refreshed inside the geometry builder (byte-matching the
        # original update, which assigned it before the empty-mask early-return).
        if geo is None:
            self._lead = _empty_lead()
            return
        self._ratchet_from_geo(geo)

    # -- deterministic geometry: host numpy reference --------------------------------
    def _geometry_host(self, P: np.ndarray):
        """Build the per-leading-node geometry in numpy (the reference path). Returns the
        SAME dict shape as :meth:`_geometry_gpu` so :meth:`_ratchet_from_geo` is shared."""
        cof = self.cof
        # Recompute rim from the LIVE geometry each tick (review fix #5): rim was frozen at init,
        # so division daughters (new cof ids the driver sets in self.cof) never crawled. Live
        # detection is ~identical when there is no division (centroids stable) and lets daughters
        # join the rim. (Actin pool stays sized from the initial n_rim; it caps gracefully if rim grows.)
        rim_cells = detect_rim_cells(self._cell_centroids(P), self.z0, self.R,
                                     self.p.rim_contact_band)
        self.rim_cells = rim_cells           # refresh before any early-return (orig order)
        rim_node_mask = np.isin(cof, rim_cells)
        if not rim_node_mask.any():
            return None
        sph = P[rim_node_mask].mean(axis=0)
        sx, sy = float(sph[0]), float(sph[1])

        lead_node_id, lead_cc, lead_out, lead_proj, lead_rp, lead_cell = (
            [], [], [], [], [], [])
        for c in rim_cells:
            c = int(c)
            cnodes = np.where(cof == c)[0]
            if cnodes.size == 0:
                continue
            npos = P[cnodes]
            cc = npos.mean(axis=0)
            out = np.array([cc[0] - sx, cc[1] - sy, 0.0])
            on = np.linalg.norm(out)
            if on < 1e-18:
                continue
            out /= on
            rel = npos - cc
            proj = rel[:, 0] * out[0] + rel[:, 1] * out[1]
            rel_xy = np.hypot(rel[:, 0], rel[:, 1])
            basal = np.abs(npos[:, 2] - self.z_basal) <= self.basal_band
            lead = basal & (proj >= self.p.lead_frac * np.maximum(rel_xy, 1e-18))
            lead_local = np.where(lead)[0]
            for ki in lead_local:
                ln = int(cnodes[ki])
                lead_node_id.append(ln)
                lead_cc.append((cc[0], cc[1]))
                lead_out.append((out[0], out[1]))
                lead_proj.append(float(proj[ki]))
                lead_rp.append(P[ln].copy())
                lead_cell.append(c)
        return {
            "rim_cells": rim_cells,
            "lead_node_id": np.asarray(lead_node_id, dtype=np.int64),
            "lead_ccx": np.asarray([p[0] for p in lead_cc], dtype=np.float64),
            "lead_ccy": np.asarray([p[1] for p in lead_cc], dtype=np.float64),
            "lead_ox": np.asarray([p[0] for p in lead_out], dtype=np.float64),
            "lead_oy": np.asarray([p[1] for p in lead_out], dtype=np.float64),
            "lead_proj": np.asarray(lead_proj, dtype=np.float64),
            "lead_rp": (np.asarray(lead_rp, dtype=np.float64).reshape(-1, 3)
                        if lead_rp else np.empty((0, 3))),
            "lead_cell": np.asarray(lead_cell, dtype=np.int64),
        }

    # -- deterministic geometry: GPU Warp kernels ------------------------------------
    def _geometry_gpu(self, *, P=None, pos_d=None):
        """Build the same geometry on the GPU via
        :func:`dcm_lamellipodium_ratchet_warp.compute_lamellipodium_geometry_gpu`,
        reading device positions (no full ``pos_d.numpy()``)."""
        from ffn_sim.warp_port.dcm_lamellipodium_ratchet_warp import (
            compute_lamellipodium_geometry_gpu)
        g = compute_lamellipodium_geometry_gpu(
            pos=(np.asarray(P, dtype=np.float64) if pos_d is None else None),
            cof=self.cof, n_cells=self.n_cells,
            z0=self.z0, R=self.R, contact_band=self.p.rim_contact_band,
            z_basal=self.z_basal, basal_band=self.basal_band, lead_frac=self.p.lead_frac,
            device=self.device, pos_d=pos_d)
        self.rim_cells = g["rim_cells"]      # refresh before any early-return (orig order)
        if g["sph"] is None:
            return None
        return g  # already the shared geo dict shape (lead_node_id / lead_* / rim_cells)

    # -- the RNG seed/advance ratchet (shared by host + GPU geometry) ----------------
    def _ratchet_from_geo(self, geo: dict) -> None:
        """SEED/ADVANCE the actin front from a per-leading-node geometry dict, then write
        ``self._lead``. Identical logic to the original host inner loop — the only change
        is that ``(cc, out, proj, rp)`` come PRE-COMPUTED (host numpy OR GPU kernels)
        instead of being recomputed inline, so both paths share one RNG stream and produce
        byte-identical anchor pools. Leading nodes MUST arrive grouped by cell (rim order),
        node id ascending within a cell — both geometry builders guarantee this."""
        lead_node_id = geo["lead_node_id"]
        ccx = geo["lead_ccx"]; ccy = geo["lead_ccy"]
        ox = geo["lead_ox"]; oy = geo["lead_oy"]
        node_proj_all = geo["lead_proj"]
        rp_all = geo["lead_rp"]
        cell_all = geo["lead_cell"]

        lead_idx_all, l_ccx, l_ccy, l_ox, l_oy, l_proj, l_cell = (
            [], [], [], [], [], [], [])
        pool_full = self.n_used >= self.n_pool

        L = lead_node_id.size
        ki = 0
        while ki < L:
            c = int(cell_all[ki])
            ccc0 = float(ccx[ki]); ccc1 = float(ccy[ki])
            out0 = float(ox[ki]); out1 = float(oy[ki])
            # this cell's active actin + its outward projection from the cell centroid
            cmask = self.actin_cell[: self.n_used] == c
            cact = self.actin_xyz[: self.n_used][cmask] if self.n_used else np.empty((0, 3))
            proj_a = ((cact[:, 0] - ccc0) * out0 + (cact[:, 1] - ccc1) * out1
                      if cact.shape[0] else np.empty(0))

            # consume all leading nodes of this cell (contiguous block)
            while ki < L and int(cell_all[ki]) == c:
                ln = int(lead_node_id[ki])
                rp = rp_all[ki]
                node_proj = float(node_proj_all[ki])
                fb_proj = None
                if cact.shape[0] > 0:
                    dist = np.linalg.norm(cact - rp, axis=1)
                    near = dist <= self.p.tether_radius
                    if near.any():
                        fb_proj = float(np.max(np.where(near, proj_a, -np.inf)))

                if not pool_full:
                    if fb_proj is None:
                        site = np.array([rp[0] + self.seed_off * out0,
                                         rp[1] + self.seed_off * out1, self.z_basal])
                        self._activate(c, site)
                        self.n_seeded += 1
                        cmask = self.actin_cell[: self.n_used] == c
                        cact = self.actin_xyz[: self.n_used][cmask]
                        proj_a = ((cact[:, 0] - ccc0) * out0 + (cact[:, 1] - ccc1) * out1)
                    elif (fb_proj - node_proj) < self.catch and self._rng.uniform() < self.p_advance:
                        site = np.array([ccc0 + (fb_proj + self.p.l0) * out0,
                                         ccc1 + (fb_proj + self.p.l0) * out1, self.z_basal])
                        self._activate(c, site)
                        self.n_advanced += 1
                        cmask = self.actin_cell[: self.n_used] == c
                        cact = self.actin_xyz[: self.n_used][cmask]
                        proj_a = ((cact[:, 0] - ccc0) * out0 + (cact[:, 1] - ccc1) * out1)
                    if self.n_used >= self.n_pool:
                        pool_full = True

                lead_idx_all.append(ln)
                l_ccx.append(ccc0); l_ccy.append(ccc1)
                l_ox.append(out0); l_oy.append(out1)
                l_proj.append(node_proj); l_cell.append(c)
                ki += 1

        self._lead = {
            "idx": np.asarray(lead_idx_all, dtype=np.int32),
            "ccx": np.asarray(l_ccx, dtype=np.float64),
            "ccy": np.asarray(l_ccy, dtype=np.float64),
            "ox": np.asarray(l_ox, dtype=np.float64),
            "oy": np.asarray(l_oy, dtype=np.float64),
            "proj": np.asarray(l_proj, dtype=np.float64),
            "cell": np.asarray(l_cell, dtype=np.int32),
        }

    def _activate(self, cell: int, site: np.ndarray) -> None:
        i = self.n_used
        self.actin_cell[i] = cell
        self.actin_xyz[i] = site
        self.n_used += 1

    @property
    def n_lead(self) -> int:
        return int(self._lead["idx"].size)

    # -- device sync --------------------------------------------------------------
    def upload(self, device: str):
        """Push the anchor pool + leading-node geometry to device arrays. Returns a
        dict of warp arrays (allocated once, reused/resized as the pool/lead grow)."""
        import warp as wp
        L = self.n_lead
        nu = self.n_used
        if self._dev is None or self._dev["cap_lead"] < L or self._dev["cap_act"] < nu:
            # geometric (×2) growth like CadherinBondHost/EcmClutchHost (review fix #10) — the
            # append-only actin pool grew every tick, so max(L,old) reallocated the whole geometry
            # block each cadence; doubling amortises it to O(log) reallocations.
            cap_lead = max(L, 1, 2 * (self._dev["cap_lead"] if self._dev else 0))
            cap_act = max(nu, 1, 2 * (self._dev["cap_act"] if self._dev else 0))
            self._dev = {
                "cap_lead": cap_lead, "cap_act": cap_act,
                "lead_rp": wp.zeros(cap_lead, dtype=wp.vec3d, device=device),
                "lead_idx": wp.zeros(cap_lead, dtype=wp.int32, device=device),
                "lead_ccx": wp.zeros(cap_lead, dtype=wp.float64, device=device),
                "lead_ccy": wp.zeros(cap_lead, dtype=wp.float64, device=device),
                "lead_ox": wp.zeros(cap_lead, dtype=wp.float64, device=device),
                "lead_oy": wp.zeros(cap_lead, dtype=wp.float64, device=device),
                "lead_proj": wp.zeros(cap_lead, dtype=wp.float64, device=device),
                "lead_cell": wp.zeros(cap_lead, dtype=wp.int32, device=device),
                "actin": wp.zeros(cap_act, dtype=wp.vec3d, device=device),
                "actin_cell": wp.zeros(cap_act, dtype=wp.int32, device=device),
            }
        d = self._dev
        if L > 0:
            d["lead_idx"].assign(self._lead["idx"])
            d["lead_ccx"].assign(self._lead["ccx"])
            d["lead_ccy"].assign(self._lead["ccy"])
            d["lead_ox"].assign(self._lead["ox"])
            d["lead_oy"].assign(self._lead["oy"])
            d["lead_proj"].assign(self._lead["proj"])
            d["lead_cell"].assign(self._lead["cell"])
        if nu > 0:
            d["actin"].assign(np.ascontiguousarray(self.actin_xyz[:nu]))
            d["actin_cell"].assign(self.actin_cell[:nu].astype(np.int32))
        d["n_lead"] = L
        d["n_used"] = nu
        return d


def _empty_lead() -> dict:
    z32 = np.empty(0, dtype=np.int32)
    z64 = np.empty(0, dtype=np.float64)
    return {"idx": z32, "ccx": z64, "ccy": z64, "ox": z64, "oy": z64,
            "proj": z64, "cell": z32}
