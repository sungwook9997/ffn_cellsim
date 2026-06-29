"""Host-side explicit FILOPODIA protrusion management for the Warp DCM loop (B3).

A filopodium is a thin actin-bundle finger: it nucleates at a cell-surface node,
elongates outward at a polymerization velocity, and its TIP probes the environment
and adheres on contact. This host owns that lifecycle (the same host-managed
protrusion pattern as :class:`dcm_lamellipodium_host.LamellipodiumHost` — tips are
advancing anchor points, NOT BAOAB DOFs); the per-step tip-adhesion spring force is
applied by the device kernels in :mod:`dcm_filopodia_warp`.

Mechanism (fine-grained, faithful to the unified actin-architecture view):
  1. **SEED** — at a low cadence (``batch_steps``), nucleate new filopodia from
     surface nodes that face OUTWARD (away from the cell centroid), at a Poisson-like
     per-node rate. The tip starts at the surface node, growing along the outward
     surface normal direction.
  2. **EXTEND** — each tick, advance every un-adhered tip outward by
     ``v_poly · batch_steps · dt`` (Brownian-ratchet polymerization, KB-3.6), up to a
     maximum length ``L_max`` (KB-3.8 filopodial reach), then it caps/retracts.
  3. **PROBE + ADHERE** — for each un-adhered tip, test (a) the substrate (z ≤ z0 +
     grip band) → form a node-to-plane clutch at the dish, and (b) the nearest FACE of
     ANOTHER cell within ``tip_capture`` (Ericson closest-point, the SAME node-FACE
     geometry as the contact kernel) → form a tip-face adhesion at that face's
     barycentric grip point.
  4. **DETACH** — break a tip adhesion when its spring stretches past a multiple of
     the rest length (distance-based; documented choice — see below) or when its
     base/face cell dies (cof<0). A detached tip is freed back to the dormant pool.

The host exposes two device-array bundles via :meth:`upload`:
  * face adhesions: ``(base_idx, face_id, bary)`` for ``filopodia_tip_face_force_kernel``.
  * substrate clutches: ``(base_idx, anchor)`` for ``filopodia_tip_plane_force_kernel``.

Detach policy (DOCUMENTED CHOICE). We use a **distance-based** break, not a
catch-slip rate: a filopodial tip complex is a transient probe, and a length-based
rupture (``L > break_factor·L_grip``) is the leanest faithful model that keeps tips
from dragging the base unboundedly — the force is capped AND the adhesion releases
once over-stretched. The sibling ECM clutch uses Pereverzev catch-slip because it
models a load-bearing integrin focal adhesion with a measured off-rate law; the
filopodial tip-complex off-rate is NOT in the KB (see PROVISIONAL list), so we do
not invent a rate. Distance break is the conservative, magic-number-free choice.

Sanity Gate
-----------
Dimensional analysis:
  * ``v_poly`` [m/s] · (``batch_steps``·``dt``) [s] = advance step [m] — consistent.
    ``p_extend`` and ``p_seed`` are unitless probabilities clamped to [0,1].
  * ``L_max`` [m], ``tip_capture`` [m], ``grip_band`` [m] are lengths; ``break_factor``
    unitless. Tip stretch ``L`` [m] compared to ``break_factor·L_grip`` [m] — consistent.
Boundary cases:
  * no outward surface node (degenerate centroid) → no seed, empty arrays (kernels
    no-op on n=0). A 1-cell geometry can only form SUBSTRATE adhesions (no other cell
    face) — face arrays stay empty, exercised by the self-test.
  * pool exhausted → seeding caps gracefully (high-water mark), no overflow.
  * a tip that reaches neither substrate nor a face stays free and keeps extending
    until ``L_max``, then retracts (no force until it adheres).
Sign-sense:
  * outward dir = (node − cell_centroid) normalised → tips grow AWAY from the cell.
  * face grip is the closest point on the OTHER cell's face (cof != base cell), so
    the resulting spring (kernel) is attractive toward the neighbour — cells approach.
  * substrate anchor z is clamped to z0 (the dish), so the clutch pulls DOWN/along to
    the dish, never into free space.

SI units throughout.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class FilopodiaParams:
    """Resolved filopodia params. KB anchors: KB-3.8 (filopodium structure / reach),
    KB-3.6 (Brownian-ratchet polymerization velocity), KB-3.21 (shared G-actin pool).

    Values with no KB record are marked PROVISIONAL and chosen mechanistically (a
    mesh/scale bridge or a sibling-module analogue), NEVER tuned to a spreading
    outcome — see the module docstring + the B3 report.
    """
    # --- protrusion kinematics (KB anchored) ---
    v_poly: float = 0.1e-6              # m/s  filopodial tip polymerization velocity.
    #   KB-3.6: single-filament v0 ~30 nm/s; bundled filopodia extend ~50-200 nm/s.
    #   0.1 µm/s = 100 nm/s sits in that band (Mallavarapu-Mitchison growth phase).
    l0: float = 0.2e-6                 # m    initial protrusion length at seed (~1 bundle
    #   segment; KB-3.8 filopodium width 100-300 nm sets the natural seed scale).
    L_max: float = 5.0e-6              # m    KB-3.8: filopodial reach r_filo ~5 µm
    #   (length range 1-10 µm; 5 µm = the cited reach radius).
    # --- tip adhesion geometry (mesh/scale bridge; PROVISIONAL magnitudes) ---
    tip_capture: float = 0.5e-6        # m    tip↔face / tip↔dish grip range. PROVISIONAL:
    #   set to the mesh contact scale (≈ c_adh of the engine) so a tip grips a
    #   neighbour face once it reaches the contact band — a ×40 scale bridge, not a fit.
    grip_band: float = 0.5e-6          # m    z-band above z0 within which a tip grips the
    #   dish. PROVISIONAL: same mesh contact scale as tip_capture.
    k_tip: float = 1.0e-3              # N/m  filopodial tip-complex spring stiffness.
    #   PROVISIONAL: matched to the ECM integrin clutch k_fa (KU-2.7 molecular clutch
    #   ~1e-3 N/m) — a filopodial tip adhesion is a small integrin/cadherin cluster of
    #   the same molecular order. No filopodial tip-complex stiffness datum in the KB.
    force_cap: float = 5.0e-9          # N    per-adhesion force cap. PROVISIONAL: matches
    #   the lamellipodium tether_cap (5·60Pa·area_node, Gil-Redondo 2023) so a single
    #   protrusion cannot inject a non-physical impulse; same physiological force scale.
    break_factor: float = 2.0          # ×L_grip  distance-based detach threshold.
    #   PROVISIONAL: a tip over-stretched past 2× its grip length releases (transient
    #   probe). No filopodial tip-adhesion off-rate in the KB (see report).
    # --- population (PROVISIONAL — no number-per-cell / lifetime KB record) ---
    p_seed: float = 0.05               # per outward-node per tick nucleation probability.
    #   PROVISIONAL: low cadence so only a sparse rim of filopodia forms (filopodia are
    #   sparse fingers, not a continuous front like the lamellipodium). Not outcome-tuned.
    pool_per_cell: int = 30            # max concurrent filopodia per cell. PROVISIONAL:
    #   order of the KB-3.8 "10-30 fil/bundle" count reused as a per-cell finger budget;
    #   no measured filopodia-per-cell datum in the KB. Caps memory, graceful.
    batch_steps: int = 50              # extension/probe cadence (matches lamellipodium/ECM).
    seed: int = 13


def _surface_outward(P: np.ndarray, cof: np.ndarray, n_cells: int):
    """Per-node outward unit direction = (node − its cell centroid) normalised, plus
    the per-cell centroid. Nodes with a degenerate (~0) outward dir are flagged invalid.

    Returns ``(out_dir (N,3), centroid (n_cells,3), valid (N,) bool)``."""
    c = np.zeros((n_cells, 3), dtype=np.float64)
    cnt = np.zeros(n_cells, dtype=np.float64)
    live = cof >= 0
    np.add.at(c, cof[live], P[live])
    np.add.at(cnt, cof[live], 1.0)
    cnt = np.where(cnt > 0, cnt, 1.0)
    cen = c / cnt[:, None]
    out = np.zeros_like(P)
    valid = np.zeros(P.shape[0], dtype=bool)
    rel = P - cen[np.clip(cof, 0, n_cells - 1)]
    nrm = np.linalg.norm(rel, axis=1)
    ok = live & (nrm > 1e-18)
    out[ok] = rel[ok] / nrm[ok, None]
    valid[ok] = True
    return out, cen, valid


class FilopodiaHost:
    """Owns the dynamic filopodia set: seed, extend, probe/adhere, detach.

    Usage in the driver loop (additive, constructed only when enabled):
        filo = FilopodiaHost(cof=cof_a, n_cells=n_cells, faces=faces_a,
                             fcell=fcell_a, z0=z0, R=R, dt=dt)
        if s % filo.batch_steps == 0:
            filo.update(pos_d.numpy().astype(np.float64))
            filo.upload(device)
        # every step: launch filopodia_tip_face_force_kernel + _plane_force_kernel
        #             over filo._dev["n_face"] / filo._dev["n_plane"].
    """

    def __init__(self, *, cof: np.ndarray, n_cells: int, faces: np.ndarray,
                 fcell: np.ndarray, z0: float, R: float, dt: float,
                 params: FilopodiaParams | None = None,
                 use_gpu_probe: bool = False, device: str = "cpu"):
        self.p = params or FilopodiaParams()
        # opt-in GPU PROBE: when True, the O(free×n_faces) nearest-other-cell-face scan
        # in update() runs as one Warp hash-grid kernel launch (dcm_filopodia_probe_warp)
        # instead of the per-tip numpy scan. Default OFF = byte-identical host behaviour.
        self.use_gpu_probe = bool(use_gpu_probe)
        self.device = device
        self.cof = np.asarray(cof, dtype=np.int64)
        self.n_cells = int(n_cells)
        self.faces = np.asarray(faces, dtype=np.int64)
        self.fcell = np.asarray(fcell, dtype=np.int64)
        self.z0 = float(z0)
        self.R = float(R)
        self.dt = float(dt)
        self.batch_steps = self.p.batch_steps
        self.dt_batch = self.batch_steps * self.dt
        self._rng = np.random.default_rng(self.p.seed)

        n_nodes = self.cof.shape[0]
        self.n_pool = self.n_cells * self.p.pool_per_cell
        # advance per tick (band-matched, clamp ≤ L_max growth)
        self.adv = self.p.v_poly * self.dt_batch

        # --- dynamic filopodia pool (struct-of-arrays) ---
        # base_idx : surface node the finger grew from (−1 = dormant slot)
        # tip      : current tip position (3,)
        # out_dir  : outward growth direction (3,)
        # length   : current protrusion length from base
        # state    : 0 free/extending, 1 face-adhered, 2 substrate-adhered
        # tgt_face : adhered face id (state==1), else −1
        # bary     : frozen grip barycentric on tgt_face (state==1)
        # anchor   : fixed dish ligand site (state==2)
        self.base_idx = np.full(self.n_pool, -1, dtype=np.int64)
        self.tip = np.zeros((self.n_pool, 3), dtype=np.float64)
        self.out_dir = np.zeros((self.n_pool, 3), dtype=np.float64)
        self.length = np.zeros(self.n_pool, dtype=np.float64)
        self.state = np.zeros(self.n_pool, dtype=np.int64)
        self.tgt_face = np.full(self.n_pool, -1, dtype=np.int64)
        self.bary = np.zeros((self.n_pool, 3), dtype=np.float64)
        self.anchor = np.zeros((self.n_pool, 3), dtype=np.float64)
        self.alive = np.zeros(self.n_pool, dtype=bool)
        self.n_per_cell = np.zeros(self.n_cells, dtype=np.int64)

        # diagnostics
        self.n_seeded = 0
        self.n_extended = 0
        self.n_face_adhered = 0
        self.n_sub_adhered = 0
        self.n_detached = 0
        self._dev = None

    # -- closest point of p on triangle (a,b,c): Ericson, returns (cpa, bary) --
    @staticmethod
    def _closest_point_bary(p, a, b, c):
        ab = b - a; ac = c - a; ap = p - a
        d1 = ab @ ap; d2 = ac @ ap
        if d1 <= 0 and d2 <= 0:
            return a.copy(), np.array([1.0, 0.0, 0.0])
        bp = p - b
        d3 = ab @ bp; d4 = ac @ bp
        if d3 >= 0 and d4 <= d3:
            return b.copy(), np.array([0.0, 1.0, 0.0])
        vc = d1 * d4 - d3 * d2
        if vc <= 0 and d1 >= 0 and d3 <= 0:
            v = d1 / (d1 - d3)
            return a + v * ab, np.array([1.0 - v, v, 0.0])
        cp = p - c
        d5 = ab @ cp; d6 = ac @ cp
        if d6 >= 0 and d5 <= d6:
            return c.copy(), np.array([0.0, 0.0, 1.0])
        vb = d5 * d2 - d1 * d6
        if vb <= 0 and d2 >= 0 and d6 <= 0:
            w = d2 / (d2 - d6)
            return a + w * ac, np.array([1.0 - w, 0.0, w])
        va = d3 * d6 - d5 * d4
        if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
            w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
            return b + w * (c - b), np.array([0.0, 1.0 - w, w])
        denom = 1.0 / (va + vb + vc)
        v = vb * denom; w = vc * denom
        return a + ab * v + ac * w, np.array([1.0 - v - w, v, w])

    def update(self, P: np.ndarray) -> None:
        """One protrusion tick: DETACH stale, EXTEND free tips, PROBE/ADHERE, SEED new."""
        P = np.asarray(P, dtype=np.float64)
        cof = self.cof
        out_dir, cen, node_valid = _surface_outward(P, cof, self.n_cells)

        # --- 1. DETACH: dead base/face cell, or over-stretched spring ---
        for i in np.flatnonzero(self.alive):
            b = self.base_idx[i]
            if b < 0 or cof[b] < 0:
                self._free(i); continue
            if self.state[i] == 1:                       # face adhesion
                fj = self.tgt_face[i]
                if fj < 0 or self.fcell[fj] < 0 or self.fcell[fj] == cof[b]:
                    self._free(i); continue
                tri = self.faces[fj]
                grip = (P[tri[0]] * self.bary[i, 0] + P[tri[1]] * self.bary[i, 1]
                        + P[tri[2]] * self.bary[i, 2])
                stretch = np.linalg.norm(grip - P[b])
                if stretch > self.p.break_factor * max(self.length[i], self.p.l0):
                    self._free(i); continue
            elif self.state[i] == 2:                     # substrate clutch
                stretch = np.linalg.norm(self.anchor[i] - P[b])
                if stretch > self.p.break_factor * max(self.length[i], self.p.l0):
                    self._free(i); continue

        # --- 2. EXTEND free tips along their outward dir, capped at L_max ---
        free = np.flatnonzero(self.alive & (self.state == 0))
        for i in free:
            b = self.base_idx[i]
            new_len = self.length[i] + self.adv
            if new_len > self.p.L_max:                   # cap/retract: free the finger
                self._free(i); continue
            self.length[i] = new_len
            self.tip[i] = P[b] + self.out_dir[i] * new_len
            self.n_extended += 1

        # --- 3. PROBE + ADHERE free tips (substrate first, then nearest other-cell face) ---
        # (a) SUBSTRATE pass (always host): tips that reach the dish band grip it and are removed
        #     from the face-probe set. This is the substrate-first branch the GPU never touches.
        free = np.flatnonzero(self.alive & (self.state == 0))
        probe_idx = []                                # free, non-substrate tips → face probe
        for i in free:
            tip = self.tip[i]
            if (tip[2] - self.z0) <= self.p.grip_band:
                self.state[i] = 2
                site = tip.copy(); site[2] = self.z0     # grip the dish ligand at z0
                self.anchor[i] = site
                self.n_sub_adhered += 1
            else:
                probe_idx.append(i)

        # (b) FACE pass: nearest OTHER-cell face within tip_capture (node-FACE closest point).
        if self.use_gpu_probe and probe_idx:
            # GPU PROBE: one hash-grid kernel over ALL free non-substrate tips, replacing the
            # per-tip O(n_faces) numpy scan. Same acceptance + transitions as the host branch.
            from ffn_sim.dcm.dcm_filopodia_probe_warp import probe_faces_gpu
            pidx = np.asarray(probe_idx, dtype=np.int64)
            tips_xyz = self.tip[pidx]
            owncell = cof[self.base_idx[pidx]].astype(np.int32)
            bf, bb, _bd = probe_faces_gpu(
                free_tips_xyz=tips_xyz, free_tip_owncell=owncell,
                pos=P, faces=self.faces, fcell=self.fcell,
                tip_capture=self.p.tip_capture, R=self.R, device=self.device)
            for k, i in enumerate(pidx):
                best_f = int(bf[k])
                if best_f >= 0:                          # kernel already gated by tip_capture
                    self.state[i] = 1
                    self.tgt_face[i] = best_f
                    self.bary[i] = bb[k]
                    self.n_face_adhered += 1
        elif probe_idx:
            # HOST PROBE (default): precompute ALL face centroids ONCE (was recomputed per free
            # tip → O(free·n_faces) gathers, the dominant per-batch host cost). Cell-independent,
            # so hoisting is exact (each tip just indexes fc_all[cand]). No physics change.
            fc_all = (P[self.faces[:, 0]] + P[self.faces[:, 1]] + P[self.faces[:, 2]]) / 3.0
            _cand_cache = {}                          # other-cell face set is SAME for a cell's tips
            for i in probe_idx:
                tip = self.tip[i]
                b = self.base_idx[i]
                cb = cof[b]
                cand = _cand_cache.get(cb)
                if cand is None:
                    cand = np.flatnonzero((self.fcell != cb) & (self.fcell >= 0))
                    _cand_cache[cb] = cand
                if cand.size == 0:
                    continue
                # cheap pre-prune by face-centroid distance, then exact closest-point
                dcent = np.linalg.norm(fc_all[cand] - tip, axis=1)
                near = cand[dcent <= (self.p.tip_capture + self.R)]
                best_d = np.inf; best_f = -1; best_bary = None
                for fj in near:
                    t = self.faces[fj]
                    cpa, bw = self._closest_point_bary(tip, P[t[0]], P[t[1]], P[t[2]])
                    dd = np.linalg.norm(cpa - tip)
                    if dd < best_d:
                        best_d = dd; best_f = int(fj); best_bary = bw
                if best_f >= 0 and best_d <= self.p.tip_capture:
                    self.state[i] = 1
                    self.tgt_face[i] = best_f
                    self.bary[i] = best_bary
                    self.n_face_adhered += 1

        # --- 4. SEED new filopodia from outward surface nodes (low cadence) ---
        if self.alive.sum() < self.n_pool:
            outward_nodes = np.flatnonzero(node_valid)
            # only nodes whose cell still has pool headroom + outward (away from centroid)
            fire = outward_nodes[self._rng.random(outward_nodes.size) < self.p.p_seed]
            for b in fire:
                c = int(cof[b])
                if self.n_per_cell[c] >= self.p.pool_per_cell:
                    continue
                slot = self._alloc()
                if slot < 0:
                    break
                self.base_idx[slot] = int(b)
                self.out_dir[slot] = out_dir[b]
                self.length[slot] = self.p.l0
                self.tip[slot] = P[b] + out_dir[b] * self.p.l0
                self.state[slot] = 0
                self.tgt_face[slot] = -1
                self.alive[slot] = True
                self.n_per_cell[c] += 1
                self.n_seeded += 1

    def _alloc(self) -> int:
        free = np.flatnonzero(~self.alive)
        return int(free[0]) if free.size else -1

    def _free(self, i: int) -> None:
        if self.alive[i]:
            c = int(self.cof[self.base_idx[i]]) if self.base_idx[i] >= 0 else -1
            if 0 <= c < self.n_cells and self.n_per_cell[c] > 0:
                self.n_per_cell[c] -= 1
            self.n_detached += 1
        self.alive[i] = False
        self.base_idx[i] = -1
        self.state[i] = 0
        self.tgt_face[i] = -1
        self.length[i] = 0.0

    @property
    def n_alive(self) -> int:
        return int(self.alive.sum())

    @property
    def n_face(self) -> int:
        return int((self.alive & (self.state == 1)).sum())

    @property
    def n_plane(self) -> int:
        return int((self.alive & (self.state == 2)).sum())

    # -- device sync (geometric growth like the sibling hosts) --------------------
    def upload(self, device):
        """Push the two adhesion bundles to device arrays (allocated once, ×2 grow).

        Returns a dict with face arrays ``(face_base, face_id, face_bary)`` sized
        ``n_face`` and substrate arrays ``(plane_base, plane_anchor)`` sized ``n_plane``."""
        import warp as wp
        fmask = self.alive & (self.state == 1)
        pmask = self.alive & (self.state == 2)
        nf = int(fmask.sum()); npl = int(pmask.sum())
        cap_need = max(nf, npl, 1)
        if self._dev is None or self._dev["cap"] < cap_need:
            cap = max(cap_need, (self._dev["cap"] * 2 if self._dev else 1))
            self._dev = {
                "cap": cap,
                "face_base": wp.zeros(cap, dtype=wp.int32, device=device),
                "face_id": wp.zeros(cap, dtype=wp.int32, device=device),
                "face_bary": wp.zeros(cap, dtype=wp.vec3d, device=device),
                "plane_base": wp.zeros(cap, dtype=wp.int32, device=device),
                "plane_anchor": wp.zeros(cap, dtype=wp.vec3d, device=device),
            }
        d = self._dev
        if nf:
            d["face_base"].assign(self.base_idx[fmask].astype(np.int32))
            d["face_id"].assign(self.tgt_face[fmask].astype(np.int32))
            d["face_bary"].assign(np.ascontiguousarray(self.bary[fmask]))
        if npl:
            d["plane_base"].assign(self.base_idx[pmask].astype(np.int32))
            d["plane_anchor"].assign(np.ascontiguousarray(self.anchor[pmask]))
        d["n_face"] = nf
        d["n_plane"] = npl
        return d


# ---------------------------------------------------------------------------
# Self-test: tiny 2-cell geometry — seed / extend / adhere + finite force check.
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import warp as wp

    from ffn_sim.dcm.dcm_filopodia_warp import (
        filopodia_tip_face_force_kernel, filopodia_tip_plane_force_kernel)

    wp.init()
    device = "cpu"
    rng = np.random.default_rng(0)

    # Build two small icosahedron-ish cells (octahedra) sitting on the dish, near each
    # other so cell-0's outward filopodia can reach cell-1's faces.
    def octa(center, r):
        v = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0],
                      [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=np.float64) * r + center
        f = np.array([[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
                      [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]], dtype=np.int64)
        return v, f

    R = 1.0e-6
    z0 = 0.0
    c0 = np.array([0.0, 0.0, R])
    c1 = np.array([1.6e-6, 0.0, R])       # 1.6 µm apart → faces within ~0.6 µm gap
    v0, f0 = octa(c0, R)
    v1, f1 = octa(c1, R)
    pos = np.vstack([v0, v1])
    faces = np.vstack([f0, f1 + 6])
    cof = np.array([0] * 6 + [1] * 6, dtype=np.int64)
    fcell = np.array([0] * 8 + [1] * 8, dtype=np.int64)
    n_cells = 2

    # Lower cell-0 bottom node into the dish band so a substrate adhesion can form.
    pos[5] = c0 + np.array([0, 0, -R]); pos[5, 2] = 0.3e-6   # bottom vertex near z0

    p = FilopodiaParams(p_seed=1.0, batch_steps=1, tip_capture=0.8e-6,
                        grip_band=0.6e-6, v_poly=0.5e-6, L_max=3.0e-6)
    filo = FilopodiaHost(cof=cof, n_cells=n_cells, faces=faces, fcell=fcell,
                         z0=z0, R=R, dt=1.0e-3, params=p)

    print("== filopodia self-test (2-cell octahedra) ==")
    for tick in range(8):
        filo.update(pos)
        print(f" tick {tick}: alive={filo.n_alive:3d} face_adh={filo.n_face:2d} "
              f"sub_adh={filo.n_plane:2d} seeded={filo.n_seeded} "
              f"extended={filo.n_extended} detached={filo.n_detached}")

    assert filo.n_seeded > 0, "no filopodia seeded"
    assert filo.n_alive > 0, "no live filopodia after ticks"
    assert (filo.n_face + filo.n_plane) > 0, "no tip adhesions formed (face or substrate)"

    # Upload + launch both force kernels; check finite, attractive forces.
    d = filo.upload(device)
    pos_d = wp.array(pos, dtype=wp.vec3d, device=device)
    force_d = wp.zeros(pos.shape[0], dtype=wp.vec3d, device=device)
    faces_d = wp.array(faces.astype(np.int32), dtype=wp.int32, ndim=2, device=device)

    if d["n_face"] > 0:
        wp.launch(filopodia_tip_face_force_kernel, dim=d["n_face"],
                  inputs=[d["face_base"], d["face_id"], d["face_bary"], wp.int32(d["n_face"]),
                          faces_d, pos_d, wp.float64(p.k_tip), wp.float64(p.force_cap),
                          force_d], device=device)
    if d["n_plane"] > 0:
        wp.launch(filopodia_tip_plane_force_kernel, dim=d["n_plane"],
                  inputs=[d["plane_base"], d["plane_anchor"], wp.int32(d["n_plane"]),
                          wp.float64(p.k_tip), wp.float64(p.force_cap), pos_d, force_d],
                  device=device)
    wp.synchronize_device(device)
    F = force_d.numpy()
    assert np.all(np.isfinite(F)), "non-finite force produced"
    fnorm = np.linalg.norm(F, axis=1)
    print(f" force: max|F|={fnorm.max():.3e} N  total|F|sum={fnorm.sum():.3e} N  "
          f"(cap={p.force_cap:.1e})")
    assert fnorm.max() <= p.force_cap * 1.000001, "force exceeds cap"
    assert fnorm.max() > 0.0, "force kernels produced zero force despite adhesions"

    # Sign check: a face-adhered base node is pulled TOWARD the neighbour's face
    # (i.e. its force has a component toward the other cell's centroid).
    if d["n_face"] > 0:
        fb = d["face_base"].numpy()[: d["n_face"]]
        for k, b in enumerate(fb):
            cb = cof[b]
            other_cen = c1 if cb == 0 else c0
            toward = (other_cen - pos[b])
            toward /= (np.linalg.norm(toward) + 1e-18)
            proj = F[b] @ toward
            print(f"  face-base node {b} (cell {cb}): F·toward_neighbour = {proj:.3e} N")
        # at least one face adhesion pulls toward the neighbour
        assert any((F[b] @ ((c1 if cof[b] == 0 else c0) - pos[b])) > 0 for b in fb), \
            "face adhesion force is not attractive toward the neighbour"

    print("PASS: seed / extend / adhere + finite attractive tip forces verified.")
