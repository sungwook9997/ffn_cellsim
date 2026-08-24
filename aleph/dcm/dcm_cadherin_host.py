"""Host-side explicit cadherin catch-bond ensemble for the Warp DCM loop (E1).

Replaces the LUMPED cell-cell adhesion (the cohesion tent attraction + the M3 binary
crowd-pressure junction switch) with the architecturally-correct FINE-GRAINED mechanism
(CLAUDE.md hard rule; PI 2026-06-21): cell-cell adhesion is an ensemble of EXPLICIT
E-cadherin trans-dimer bonds between apposed membrane nodes of different cells, each with
the faithful Rakshit-2012 sliding-rebinding catch-slip off-rate. De-cohesion is then
EMERGENT — under the spreading traction the bonds load past the catch peak (f0≈29 pN) into
the slip regime, rupture faster than they reform, and the cells let go. No state-switch.

Design (the leanest faithful realization on the existing mesh, no extra particles): one
trans-dimer bond per apposed membrane-node PAIR (the membrane carries the cadherins); each
node hosts at most one trans-dimer ("one trans-dimer per cadherin"). Bonds are managed on
the host at the binder cadence ``batch_steps`` (KDTree partner search, per-batch stochastic
break/form); a device kernel applies the harmonic trans-dimer force every step. The
force-dependent off-rate uses ``validation.cadherin_sliding_rebinding.effective_k_off``
(exact 3×3 generator) precomputed into a force→k_off lookup table (the closed form is too
heavy per-bond-per-step; the table is faithful to interpolation error).

Repulsion / excluded volume is NOT here — it stays in the cohesion+contact kernels (run
with adhesion OFF in cadherin mode). Cadherin bonds are attractive-only (a stretched
ectodomain tether pulls; a slack one, L≤r0, is force-free).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aleph.dcm.dcm_cadherin_cluster import cluster_bd_step
from aleph.validation.cadherin_sliding_rebinding import (
    effective_k_off, RAKSHIT_W2A, CadherinCatchParams)


@dataclass
class CadherinParams:
    """Resolved trans-dimer params (defaults = resolve_cadherin_junction for MCF7:
    Iturri F_detach 6.5 nN / Rakshit f0, contact_zone_width 0.5 µm)."""
    k_trans: float = 5.84e-5       # N/m  trans-dimer ectodomain stiffness
    r0_trans: float = 0.5e-6       # m    rest length (force-free bridge)
    r_bind: float = 0.5e-6         # m    capture radius (unbound A↔B → trans-dimer)
    k_on: float = 27.96            # s⁻¹  reformation rate (= k_off(0), rest-symmetric)
    # binder cadence. The resolver's CFL-derived 8 made the host bond-management a near-per-step
    # GPU→CPU sync, defeating the Warp port's whole purpose (GPU-main rule). Raised to 50 (matching
    # the ecm-clutch / lamellipodium host-hybrid cadence): the EXACT survival probability
    # 1−exp(−k_off·Δt_batch) is correct for ANY Δt_batch when k_off is ~constant over the batch,
    # and in the overdamped slow spread the bond force (hence k_off) barely changes over 50 steps
    # (4e-4 s) — so the kinetics stay faithful while the host sync is amortised (the accepted
    # low-cadence host-hybrid, like remesh). Override via --cad-batch.
    batch_steps: int = 50
    seed: int = 7
    # ×40 mesoscale FORCE bridge (sanctioned coarse-graining). A node-pair bond is a BUNDLE of
    # ~bundle_n cadherins over the contact patch: it transmits bundle_n × the single-molecule
    # force, but its catch-slip k_off is still evaluated at the PER-MOLECULE load F/bundle_n (the
    # Rakshit SHAPE + f0 stay molecular). bundle_n=1 → legacy single-molecule force (back-compat).
    # Set so a typical junction lands in KB-4.11's audited 1–10 nN/junction (LeDuc2010, Buckley2014).
    bundle_n: float = 1.0
    # Active actomyosin junctional CONTRACTION [N per single trans-dimer, bundle-scaled at launch] —
    # the Stage-2 compaction motor the passive catch-bond lacks (RhoA/ROCK-gated junctional NMII that
    # actively pulls bonded cells together). 0 = off (passive-only, legacy). SWEEP as a controlled
    # variable anchored to junctional actomyosin tension (per-motor ~5-15 pN × engaged motors); never
    # tune to a compaction target. Applied as f_contract·bundle_n along each bond (always contracting).
    f_contract: float = 0.0
    # #1 Step 5: multiscale KMC sub-cycling. The exact survival 1−exp(−k·Δt) needs k·Δt SMALL (k_off
    # ~constant over the increment); at a large mechanics accel_dt, dt_batch=batch_steps·accel_dt can
    # be ≫ the ~36 ms bond lifetime (1/k_on), so a single KMC pass overshoots the turnover (all loaded
    # bonds break + all free nodes bind in one saturated step → wrong junction remodelling). Sub-cycle
    # instead: advance the KMC in n_sub micro-steps of δt_cad = 1/(micro_M·k_on) — a DERIVED fraction of
    # the physical cadherin timescale (grid-invariant, not tuned). n_sub=1 whenever dt_batch ≤ δt_cad
    # (base dt) → byte-identical to the single-pass version. k_on here is the CURRENT (S_accel-scaled)
    # rate, so the biology-time-acceleration wrinkle stays consistent automatically.
    subcycle: bool = True
    micro_M: int = 8               # micro-steps per cadherin timescale 1/k_on (δt_cad = 1/(M·k_on))
    micro_cap: int = 512           # hard ceiling on n_sub/call (logged if hit; no silent truncation)
    # Junction MATURATION (2026-07-09; default OFF, back-compat byte-identical). A nascent trans-dimer
    # turns over at the fast single-molecule catch-slip rate (rest k_off≈k_on≈28/s, ~36 ms), but a
    # PERSISTING junction MATURES (E-cadherin clustering + α-catenin/vinculin reinforcement, KB-4.3) and
    # its lifetime rises to 5–30 min (KB-4.11, LeDuc2010/Buckley2014). A bond accumulates AGE; the
    # matured fraction m=1−exp(−age/τ_mature) blends its off-rate from the nascent catch-slip toward a
    # mature floor: k_off_eff = k_off_catchslip(F)·[1 − m·(1 − r)], r = k_off_mature/k_off_nascent_rest =
    # (1/mature_lifetime)/k_off(0). Force-dependence is retained (a mature bond still slips under strong
    # load). This raises tissue REARRANGEMENT viscosity over the maturation timescale → slow (min–hr)
    # compaction matching real spheroids; the fast ~100 s mechanical rounding is the immature limit.
    # τ_mature + mature_lifetime are LIT-ANCHORED (KB-4.11); the resulting compaction rate is MEASURED,
    # not tuned to a target (magic-number rule).
    mature: bool = False
    tau_mature: float = 600.0      # s  maturation timescale (KB-4.11 5–30 min; mid ≈10 min)
    mature_lifetime: float = 600.0 # s  mature junction lifetime (KB-4.11) → k_off_mature = 1/this
    # LOAD-SHARING CLUSTER (2026-07-09; default OFF, back-compat byte-identical). A junction is a
    # cluster of ``bundle_n`` PARALLEL trans-dimers sharing the extension: every engaged molecule bears
    # the same per-molecule load F1=k_trans·(L−r0), unbinds at eps(F1), and an empty slot rebinds at
    # k_on; the junction is lost only when ALL molecules are simultaneously unbound (m→0). This replaces
    # the lumped "whole bundle breaks at the single-molecule rate" (which fixed the junction lifetime at
    # ~1/k_off and starved maturation) with the fine-grained cluster whose collective lifetime T(n_b)≫
    # 1/eps is EMERGENT (CLAUDE.md hard rule). ``bundle_n`` IS the cluster size n_b here. See
    # dcm_cadherin_cluster + DCM_CADHERIN_CLUSTER_REARRANGEMENT_DESIGN_2026-07-09.
    cluster: bool = False
    # S2 maturation-capacity: a nascent junction NUCLEATES small (n_nascent parallel dimers, turns over
    # fast → rearrangement-permissive) and, if it PERSISTS, RECRUITS toward the full cluster over
    # tau_mature (E-cadherin clustering + α-catenin/vinculin reinforcement, KB-4.3): N_b(age) = n_nascent
    # + (bundle_n − n_nascent)·(1 − exp(−age/tau_mature)). As capacity grows the empty slots refill at
    # k_on and the cluster lifetime rises steeply → the junction LOCKS. The RATE at which junctions lock
    # is 1/tau_mature (lit-anchored, robust) — NOT the super-exponentially-sensitive raw lifetime (see
    # DCM_CADHERIN_CLUSTER_REARRANGEMENT_DESIGN §2). Active only when cluster AND mature are both on.
    # n_nascent = the post-nucleation clustered dimer count (KB-4.3 nascent puncta); CONTROLLED VARIABLE,
    # characterised in G2/G3, flagged for a firm KB anchor — never tuned to a compaction target.
    n_nascent: int = 4
    catch: CadherinCatchParams = None   # set in __post_init__ to RAKSHIT_W2A


def _build_koff_table(catch: CadherinCatchParams, f_max_pN: float = 300.0, n: int = 512):
    """Force→k_off(F) lookup [grid F in N, k_off in s⁻¹] from the exact Rakshit generator."""
    fs = np.linspace(0.0, f_max_pN * 1e-12, n)
    koff = np.array([effective_k_off(float(f), catch) for f in fs])
    return fs, koff


class CadherinBondHost:
    """Owns the dynamic trans-dimer bond set between apposed membrane nodes of different cells.

    Usage in the driver loop (cadherin mode; cohesion/contact run adhesion-OFF):
        cad = CadherinBondHost(cof=cof_a, n_cells=..., dt=dt)
        ...
        if s % cad.batch_steps == 0:
            cad.update(pos_d.numpy()); cad.upload(device)   # break/form on host, push to device
        # every step: cadherin_bond_force_kernel(bonds_d, r0, k_trans, pos_d, force_d)
    """

    def __init__(self, *, cof: np.ndarray, n_cells: int, dt: float,
                 params: CadherinParams | None = None):
        self.p = params or CadherinParams()
        if self.p.catch is None:
            self.p.catch = RAKSHIT_W2A
        self.cof = np.asarray(cof, dtype=np.int64)
        self.n_cells = int(n_cells)
        self.dt = float(dt)
        self.batch_steps = self.p.batch_steps
        self.dt_batch = self.batch_steps * self.dt
        self.subcycle = bool(self.p.subcycle)
        self.micro_M = int(self.p.micro_M)
        self.micro_cap = int(self.p.micro_cap)
        self._subcycle_logged = False
        self._rng = np.random.default_rng(self.p.seed)
        self._fs, self._koff = _build_koff_table(self.p.catch)
        # dynamic bond set: (M,2) node-index pairs (i in cell A, j in cell B); each node ≤1 bond
        self.bonds = np.zeros((0, 2), dtype=np.int64)
        # load-sharing cluster: per-bond engaged-molecule count m (parallel to self.bonds); a new
        # junction nucleates full (m=n_b). n_b = bundle_n as an int cluster size (cluster mode only).
        self.cluster = bool(self.p.cluster)
        self.n_b = max(1, int(round(float(self.p.bundle_n))))   # full (mature) cluster capacity
        self.m = np.zeros((0,), dtype=np.int64)
        # S2 maturation-capacity: maturation is a property of the sustained CONTACT (apposition), NOT a
        # single trans-dimer's uninterrupted lifetime — a nascent junction's lifetime (~0.24 s even with
        # load sharing) is « τ_mature, so a persisting CONTACT must accumulate maturation THROUGH bond
        # turnover (break→reform). Track per-NODE contact_age [s] (incremented while the node is apposed
        # to / bonded with another cell, reset when it leaves contact); N_b for a bond derives from its
        # endpoints' contact_age. So the capacity matures over τ_mature regardless of which specific
        # trans-dimer is currently engaged → the lock RATE is 1/τ_mature (robust).
        self.cluster_mature = bool(self.p.cluster and self.p.mature)
        self.n_nascent = max(1, min(int(self.p.n_nascent), self.n_b))
        self.contact_age = np.zeros(self.cof.size, dtype=np.float64)
        self.n_formed = 0
        self.n_broken = 0
        self._dev = None

    def _koff_of(self, F: np.ndarray) -> np.ndarray:
        return np.interp(F, self._fs, self._koff)

    def _nb_of_age(self, age: np.ndarray) -> np.ndarray:
        """Maturing cluster capacity N_b(age) = n_nascent + (n_b − n_nascent)·(1 − exp(−age/τ_mature)),
        rounded to an integer molecule count. Constant n_b when maturation is off (S1)."""
        if not self.cluster_mature:
            return np.full(np.shape(age), self.n_b, dtype=np.int64)
        grow = 1.0 - np.exp(-np.asarray(age) / max(self.p.tau_mature, 1e-30))
        nb = self.n_nascent + (self.n_b - self.n_nascent) * grow
        return np.clip(np.rint(nb), self.n_nascent, self.n_b).astype(np.int64)

    def _apposed_mask(self, P: np.ndarray) -> np.ndarray:
        """(N,) bool: a node is IN CONTACT iff within r_bind of a live DIFFERENT-cell node, OR currently
        bonded (a loaded trans-dimer can stretch past r_bind but is still an apposed contact). This is
        the sustained-contact signal that accumulates maturation across bond turnover."""
        from scipy.spatial import cKDTree
        cof = self.cof
        app = np.zeros(cof.size, dtype=bool)
        idx = np.flatnonzero(cof >= 0)
        if idx.size >= 2:
            pairs = cKDTree(P[idx]).query_pairs(r=self.p.r_bind, output_type="ndarray")
            if pairs.size:
                a = idx[pairs[:, 0]]; b = idx[pairs[:, 1]]
                diff = cof[a] != cof[b]                    # only across a cell-cell interface
                a, b = a[diff], b[diff]
                app[a] = True; app[b] = True
        if self.bonds.shape[0]:                            # a bonded node is apposed regardless of stretch
            app[self.bonds[:, 0]] = True; app[self.bonds[:, 1]] = True
        return app

    def _n_subcycle(self):
        """(n_sub, dt_sub) for advancing the KMC over dt_batch. δt_cad = 1/(micro_M·k_ref) with
        k_ref = max(k_on, koff at the MAX FORMABLE stretch k_trans·(r_bind−r0)) — the fastest rate a
        *persisting* bond can have (bonds beyond r_bind can't form, and koff→millions there = instant
        rupture, irrelevant to the steady population). Resolving to k_ref keeps koff·δt ≤ 1/micro_M for
        every formable bond, so 1−exp(−koff·δt)≈koff·δt stays linear and the steady fraction is
        dt-invariant. n_sub=1 whenever dt_batch is already ≤ δt_cad (base dt) → byte-identical."""
        kon = float(self.p.k_on)
        if not (self.subcycle and self.dt_batch > 0.0 and kon > 0.0):
            return 1, float(self.dt_batch)
        f_maxform = float(self.p.k_trans) * max(0.0, float(self.p.r_bind) - float(self.p.r0_trans))
        k_ref = max(kon, float(np.interp(f_maxform, self._fs, self._koff)))
        n_sub = int(np.ceil(self.dt_batch * float(self.micro_M) * k_ref))
        n_sub = max(1, min(n_sub, self.micro_cap))
        dt_sub = float(self.dt_batch) / float(n_sub)
        if n_sub > 1 and not self._subcycle_logged:
            print(f"  [cad-subcycle] dt_batch={self.dt_batch:.2e}s → {n_sub} micro-steps @ "
                  f"δt_cad={dt_sub*1e3:.3f}ms (k_ref={k_ref:.0f}/s = max[k_on, slip@r_bind])"
                  f"{' (CAPPED — fastest bonds under-resolved)' if n_sub >= self.micro_cap else ''}",
                  flush=True)
            self._subcycle_logged = True
        return n_sub, dt_sub

    def update(self, P: np.ndarray) -> None:
        """One binder update over ``dt_batch``. #1 Step 5: SUB-CYCLED at the cadherin timescale — the
        break+form tick is applied n_sub times at δt_cad=1/(micro_M·k_on) so a large mechanics dt does
        not overshoot the ~36ms turnover. n_sub=1 at base dt → byte-identical (same RNG draw order)."""
        P = np.asarray(P, dtype=np.float64)
        n_sub, dt_sub = self._n_subcycle()
        for _ in range(n_sub):
            self._tick(P, dt_sub)

    def _tick(self, P: np.ndarray, dt: float) -> None:
        """One break+form pass advancing the KMC by ``dt`` at frozen positions ``P``."""
        cof = self.cof
        # --- MATURATION: age the sustained contacts (per node), reset those out of contact (S2) ---
        if self.cluster_mature:
            app = self._apposed_mask(P)
            self.contact_age[app] += dt
            self.contact_age[~app] = 0.0
        # --- BREAK (catch-slip, force-dependent) ---
        if self.bonds.shape[0]:
            i = self.bonds[:, 0]; j = self.bonds[:, 1]
            # a bond whose node was deactivated (remesh collapse, cof<0) is dropped
            live = (cof[i] >= 0) & (cof[j] >= 0)
            i, j = i[live], j[live]
            L = np.linalg.norm(P[i] - P[j], axis=1)
            if self.cluster:
                # load-sharing cluster: per-molecule load F1=k_trans·(L−r0) (parallel springs at one
                # extension) drives each molecule's catch-slip unbind; empty slots rebind at k_on; the
                # junction dies only when all molecules are simultaneously unbound (m→0). The capacity
                # N_b matures over τ_mature via the endpoints' CONTACT age (S2), so a persisting contact
                # recruits + locks even across bond turnover; a fresh contact turns over at the nascent size.
                f1 = self.p.k_trans * np.maximum(0.0, L - self.p.r0_trans)
                p_off = 1.0 - np.exp(-self._koff_of(f1) * dt)
                p_on = 1.0 - np.exp(-self.p.k_on * dt)
                nb = (self._nb_of_age(np.minimum(self.contact_age[i], self.contact_age[j]))
                      if self.cluster_mature else self.n_b)
                m = cluster_bd_step(self.m[live], nb, p_off, p_on, self._rng)
                keep = m > 0
                self.n_broken += int((~keep).sum())
                self.bonds = np.stack([i[keep], j[keep]], axis=1)
                self.m = m[keep]
            else:
                F = self.p.k_trans * np.maximum(0.0, L - self.p.r0_trans)
                p_break = 1.0 - np.exp(-self._koff_of(F) * dt)
                keep = self._rng.random(i.size) >= p_break
                self.n_broken += int((~keep).sum())
                self.bonds = np.stack([i[keep], j[keep]], axis=1)
        # --- FORM (rest-symmetric on-rate over apposed unbonded different-cell pairs) ---
        bonded = np.zeros(P.shape[0], dtype=bool)
        if self.bonds.shape[0]:
            bonded[self.bonds[:, 0]] = True
            bonded[self.bonds[:, 1]] = True
        free = np.flatnonzero((cof >= 0) & (~bonded))
        if free.size >= 2:
            self._form(P, free, bonded, dt)

    def _form(self, P, free, bonded, dt: float):
        from scipy.spatial import cKDTree
        cof = self.cof
        pts = P[free]
        tree = cKDTree(pts)
        pairs = tree.query_pairs(r=self.p.r_bind, output_type="ndarray")
        if pairs.size == 0:
            return
        a = free[pairs[:, 0]]; b = free[pairs[:, 1]]
        diff = cof[a] != cof[b]                       # only ACROSS a cell-cell interface
        a, b = a[diff], b[diff]
        if a.size == 0:
            return
        p_on = 1.0 - np.exp(-self.p.k_on * dt)
        fire = self._rng.random(a.size) < p_on
        a, b = a[fire], b[fire]
        # greedily accept pairs keeping the "≤1 trans-dimer per node" invariant
        new = []
        for ia, ib in zip(a.tolist(), b.tolist()):
            if not bonded[ia] and not bonded[ib]:
                bonded[ia] = True; bonded[ib] = True
                new.append((ia, ib))
        if new:
            self.n_formed += len(new)
            new_arr = np.array(new, dtype=np.int64)
            self.bonds = np.concatenate([self.bonds, new_arr], axis=0)
            if self.cluster:
                # nucleate at the CONTACT's current capacity: a fresh contact → n_nascent; a re-forming
                # bond on an already-matured contact → its grown N_b (the clustered cadherins re-engage).
                if self.cluster_mature:
                    m0 = self._nb_of_age(
                        np.minimum(self.contact_age[new_arr[:, 0]], self.contact_age[new_arr[:, 1]]))
                else:
                    m0 = np.full(len(new), self.n_b, dtype=np.int64)
                self.m = np.concatenate([self.m, m0])

    @property
    def n_bonds(self) -> int:
        if getattr(self, "_n_gpu", None) is not None:   # GPU-native path: device bond count
            return int(self._n_gpu)
        return int(self.bonds.shape[0])

    def upload(self, device):
        """Push the current bond pairs to a device int32 (M,2) array (realloc as M grows). In cluster
        mode also push the per-bond engaged count m (the CPU-device / host-update path, so the cluster
        force kernel can read cad.m_device())."""
        import warp as wp
        M = self.n_bonds
        if self._dev is None or self._dev["cap"] < M:
            cap = max(M, 1, (self._dev["cap"] * 2 if self._dev else 0))
            self._dev = {"cap": cap,
                         "bonds": wp.zeros(cap, dtype=wp.vec2i, device=device)}
            if self.cluster:
                self._dev["m"] = wp.zeros(cap, dtype=wp.int32, device=device)
        if self.cluster and "m" not in self._dev:
            self._dev["m"] = wp.zeros(self._dev["cap"], dtype=wp.int32, device=device)
        if M:
            self._dev["bonds"].assign(self.bonds.astype(np.int32))
            if self.cluster:
                self._dev["m"].assign(self.m.astype(np.int32))
        self._dev["n"] = M
        return self._dev

    # ---------------------------------------------------------------- GPU-native path -----
    def _ensure_gpu(self, N: int, device):
        """Allocate the device buffers for the GPU break/form path (bonds ping-pong + scratch)."""
        import warp as wp
        from aleph.dcm.dcm_cadherin_gpu import build_koff_device
        cap = int(N)                                   # ≤1 bond/node → ≤N/2 bonds; N is a safe cap
        if self._dev is None or self._dev.get("cap", 0) < cap or "bondsA" not in self._dev:
            koff_d, fs0, df, nk = build_koff_device(self._fs, self._koff, device)
            n0 = int(self._dev["n"]) if self._dev else 0
            bA = wp.zeros(cap, dtype=wp.vec2i, device=device)
            if self._dev is not None and self._dev.get("n"):    # carry any existing bonds over
                bA.assign(np.resize(self.bonds.astype(np.int32), (cap, 2)))
            self._dev = {
                "cap": cap, "n": n0,
                "bondsA": bA, "bondsB": wp.zeros(cap, dtype=wp.vec2i, device=device),
                "ageA": wp.zeros(cap, dtype=wp.float64, device=device),   # per-bond age [s], ping-pongs
                "ageB": wp.zeros(cap, dtype=wp.float64, device=device),   # with bonds (maturation)
                "which": "A", "bonds": bA,
                "count": wp.zeros(1, dtype=wp.int32, device=device),
                "bonded": wp.zeros(N, dtype=wp.int32, device=device),
                "partner": wp.zeros(N, dtype=wp.int32, device=device),
                "koff_d": koff_d, "fs0": fs0, "df": df, "nk": nk,
                "grid": wp.HashGrid(48, 48, 48, device=device),   # OWN grid, built at r_bind
            }
        return self._dev

    def update_gpu(self, pos_d, cof_d, node_f32, N: int, batch_index: int, device) -> None:
        """GPU-native binder tick — break + mutual-nearest form on the device, advancing ``dt_batch``
        of physical time. #1 Step 5: SUB-CYCLED at the cadherin timescale — when dt_batch ≫ 1/k_on
        (large mechanics dt) the KMC is resolved in ``n_sub`` micro-steps of δt_cad=1/(micro_M·k_on)
        rather than one saturated pass (which would break every loaded bond + bind every free node in a
        single 4 s jump). Positions are frozen over the sub-cycle (operator split: the fast bonds relax
        at xₙ), so the r_bind grid is built ONCE. n_sub=1 at base dt → byte-identical (same RNG salt)."""
        import warp as wp
        from aleph.dcm.dcm_cadherin_gpu import (cad_break_kernel, cad_partner_kernel, cad_form_kernel)
        d = self._ensure_gpu(N, device)
        d["grid"].build(points=node_f32, radius=float(self.p.r_bind))   # positions frozen → build once
        n_sub, dt_sub = self._n_subcycle()
        p_on = 1.0 - np.exp(-float(self.p.k_on) * dt_sub)
        # maturation: r = k_off_mature / k_off_nascent_rest (=koff at F=0 ≈ k_on). LIT-anchored lifetimes.
        mature_on = 1 if bool(self.p.mature) else 0
        koff_rest = float(self._koff[0]) if self._koff[0] > 0 else float(self.p.k_on)
        mature_factor = ((1.0 / float(self.p.mature_lifetime)) / koff_rest) if mature_on else 1.0
        for isub in range(n_sub):
            cur = d["bondsA"] if d["which"] == "A" else d["bondsB"]
            other = d["bondsB"] if d["which"] == "A" else d["bondsA"]
            cur_age = d["ageA"] if d["which"] == "A" else d["ageB"]
            other_age = d["ageB"] if d["which"] == "A" else d["ageA"]
            n = int(d["n"])
            d["count"].zero_(); d["bonded"].zero_()
            tick = batch_index if n_sub == 1 else (batch_index * self.micro_cap + isub)   # base dt: original salt
            sb = wp.int32(((2 * tick) * N) % 2000000000)          # disjoint RNG salt ranges: break vs form
            sf = wp.int32(((2 * tick + 1) * N) % 2000000000)
            wp.launch(cad_break_kernel, dim=max(n, 1), inputs=[
                cur, wp.int32(n), pos_d, cof_d,
                wp.float64(self.p.k_trans), wp.float64(self.p.r0_trans), wp.float64(dt_sub),
                d["koff_d"], wp.float64(d["fs0"]), wp.float64(d["df"]), wp.int32(d["nk"]),
                wp.int32(self.p.seed), sb,
                cur_age, other_age, wp.int32(mature_on),
                wp.float64(self.p.tau_mature), wp.float64(mature_factor),
                other, d["count"], d["bonded"]], device=device)
            wp.launch(cad_partner_kernel, dim=N, inputs=[
                d["grid"].id, node_f32, pos_d, cof_d, d["bonded"], wp.float64(self.p.r_bind),
                d["partner"]], device=device)
            wp.launch(cad_form_kernel, dim=N, inputs=[
                d["partner"], wp.int32(N), wp.float64(p_on), wp.int32(self.p.seed), sf,
                wp.int32(d["cap"]), other, d["count"], other_age], device=device)
            new_n = min(int(d["count"].numpy()[0]), d["cap"])   # scalar sync per micro-step
            self.n_broken += max(0, n - new_n)
            self.n_formed += max(0, new_n - n)
            d["which"] = "B" if d["which"] == "A" else "A"
            d["bonds"] = other
            d["n"] = new_n
        self._n_gpu = int(d["n"])

    # -------------------------------------------------- GPU-native CLUSTER path (S3) -----
    def _ensure_gpu_cluster(self, N: int, device):
        """Device buffers for the load-sharing cluster path: bonds + engaged-count m ping-pong, per-node
        contact_age (persistent maturation clock) + apposed mask, partner/bonded scratch, koff table, grid."""
        import warp as wp
        from aleph.dcm.dcm_cadherin_gpu import build_koff_device
        cap = int(N)
        if self._dev is None or self._dev.get("cap", 0) < cap or "mA" not in self._dev:
            koff_d, fs0, df, nk = build_koff_device(self._fs, self._koff, device)
            bA = wp.zeros(cap, dtype=wp.vec2i, device=device)
            mA = wp.zeros(cap, dtype=wp.int32, device=device)
            n0 = 0
            if self._dev is not None and self._dev.get("n"):     # carry existing bonds/m over
                n0 = int(self._dev["n"])
                bA.assign(np.resize(self.bonds.astype(np.int32), (cap, 2)))
                mA.assign(np.resize(self.m.astype(np.int32), (cap,)) if self.m.size else np.zeros(cap, np.int32))
            self._dev = {
                "cap": cap, "n": n0, "cluster": True,
                "bondsA": bA, "bondsB": wp.zeros(cap, dtype=wp.vec2i, device=device),
                "mA": mA, "mB": wp.zeros(cap, dtype=wp.int32, device=device),
                "which": "A", "bonds": bA, "m": mA,
                "contact_age": wp.zeros(N, dtype=wp.float64, device=device),   # persistent per-node
                "apposed": wp.zeros(N, dtype=wp.int32, device=device),
                "bonded": wp.zeros(N, dtype=wp.int32, device=device),
                "partner": wp.zeros(N, dtype=wp.int32, device=device),
                "count": wp.zeros(1, dtype=wp.int32, device=device),
                "koff_d": koff_d, "fs0": fs0, "df": df, "nk": nk,
                "grid": wp.HashGrid(48, 48, 48, device=device),
            }
        return self._dev

    def update_gpu_cluster(self, pos_d, cof_d, node_f32, N: int, batch_index: int, device) -> None:
        """GPU-native load-sharing cluster binder tick (S3): per-node contact-age maturation + per-bond
        engaged-count birth-death break + mutual-nearest form. Statistical parity with the host
        _tick cluster path (Binomial via per-molecule Bernoulli sums; different RNG stream)."""
        import warp as wp
        from aleph.dcm.dcm_cadherin_gpu import (
            cad_break_cluster_kernel, cad_form_cluster_kernel, cad_partner_kernel,
            cad_apposition_kernel, cad_bonded_appose_kernel, contact_age_update_kernel)
        d = self._ensure_gpu_cluster(N, device)
        d["grid"].build(points=node_f32, radius=float(self.p.r_bind))
        n_sub, dt_sub = self._n_subcycle()
        p_on = 1.0 - np.exp(-float(self.p.k_on) * dt_sub)
        mature_on = 1 if (self.p.cluster and self.p.mature) else 0
        n_nascent = max(1, min(int(self.p.n_nascent), self.n_b))
        for isub in range(n_sub):
            cur_b = d["bondsA"] if d["which"] == "A" else d["bondsB"]
            oth_b = d["bondsB"] if d["which"] == "A" else d["bondsA"]
            cur_m = d["mA"] if d["which"] == "A" else d["mB"]
            oth_m = d["mB"] if d["which"] == "A" else d["mA"]
            n = int(d["n"])
            # --- maturation: age sustained contacts (apposed OR bonded), reset the rest ---
            if mature_on != 0:                            # contact_age is unused when maturation off
                d["apposed"].zero_()
                wp.launch(cad_apposition_kernel, dim=N, inputs=[
                    d["grid"].id, node_f32, cof_d, wp.float64(self.p.r_bind), d["apposed"]], device=device)
                if n > 0:
                    wp.launch(cad_bonded_appose_kernel, dim=n, inputs=[cur_b, wp.int32(n), d["apposed"]],
                              device=device)
                wp.launch(contact_age_update_kernel, dim=N, inputs=[
                    d["apposed"], wp.float64(dt_sub), d["contact_age"]], device=device)
            # --- break (cluster birth-death) ---
            d["count"].zero_(); d["bonded"].zero_()
            tick = batch_index if n_sub == 1 else (batch_index * self.micro_cap + isub)
            sb = wp.int32(((2 * tick) * N) % 2000000000)
            sf = wp.int32(((2 * tick + 1) * N) % 2000000000)
            wp.launch(cad_break_cluster_kernel, dim=max(n, 1), inputs=[
                cur_b, wp.int32(n), cur_m, pos_d, cof_d, d["contact_age"],
                wp.float64(self.p.k_trans), wp.float64(self.p.r0_trans), wp.float64(dt_sub),
                d["koff_d"], wp.float64(d["fs0"]), wp.float64(d["df"]), wp.int32(d["nk"]),
                wp.float64(self.p.k_on), wp.int32(mature_on), wp.float64(self.p.tau_mature),
                wp.int32(n_nascent), wp.int32(self.n_b), wp.int32(self.p.seed), sb,
                oth_b, oth_m, d["count"], d["bonded"]], device=device)
            # --- form (mutual-nearest; nucleate at the contact's current capacity) ---
            wp.launch(cad_partner_kernel, dim=N, inputs=[
                d["grid"].id, node_f32, pos_d, cof_d, d["bonded"], wp.float64(self.p.r_bind),
                d["partner"]], device=device)
            wp.launch(cad_form_cluster_kernel, dim=N, inputs=[
                d["partner"], wp.int32(N), d["contact_age"], wp.int32(mature_on),
                wp.float64(self.p.tau_mature), wp.int32(n_nascent), wp.int32(self.n_b),
                wp.float64(p_on), wp.int32(self.p.seed), sf, wp.int32(d["cap"]),
                oth_b, d["count"], oth_m], device=device)
            new_n = min(int(d["count"].numpy()[0]), d["cap"])
            self.n_broken += max(0, n - new_n)
            self.n_formed += max(0, new_n - n)
            d["which"] = "B" if d["which"] == "A" else "A"
            d["bonds"] = oth_b
            d["m"] = oth_m
            d["n"] = new_n
        self._n_gpu = int(d["n"])

    def m_device(self):
        """Current per-bond engaged-count device array (for the cluster force kernel)."""
        return self._dev["m"] if (self._dev is not None and "m" in self._dev) else None

    def bonds_now(self) -> np.ndarray:
        """Current (M,2) bond node-pairs — device download in GPU mode, else the host array."""
        if getattr(self, "_n_gpu", None) is not None and self._dev is not None:
            n = int(self._n_gpu)
            if n <= 0:
                return np.zeros((0, 2), np.int64)
            return self._dev["bonds"].numpy()[:n].astype(np.int64)
        return self.bonds
