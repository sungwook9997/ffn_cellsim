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

from ffn_sim.validation.cadherin_sliding_rebinding import (
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
        self.n_formed = 0
        self.n_broken = 0
        self._dev = None

    def _koff_of(self, F: np.ndarray) -> np.ndarray:
        return np.interp(F, self._fs, self._koff)

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
        # --- BREAK (catch-slip, force-dependent) ---
        if self.bonds.shape[0]:
            i = self.bonds[:, 0]; j = self.bonds[:, 1]
            # a bond whose node was deactivated (remesh collapse, cof<0) is dropped
            live = (cof[i] >= 0) & (cof[j] >= 0)
            i, j = i[live], j[live]
            L = np.linalg.norm(P[i] - P[j], axis=1)
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
            self.bonds = np.concatenate([self.bonds, np.array(new, dtype=np.int64)], axis=0)

    @property
    def n_bonds(self) -> int:
        if getattr(self, "_n_gpu", None) is not None:   # GPU-native path: device bond count
            return int(self._n_gpu)
        return int(self.bonds.shape[0])

    def upload(self, device):
        """Push the current bond pairs to a device int32 (M,2) array (realloc as M grows)."""
        import warp as wp
        M = self.n_bonds
        if self._dev is None or self._dev["cap"] < M:
            cap = max(M, 1, (self._dev["cap"] * 2 if self._dev else 0))
            self._dev = {"cap": cap,
                         "bonds": wp.zeros(cap, dtype=wp.vec2i, device=device)}
        if M:
            self._dev["bonds"].assign(self.bonds.astype(np.int32))
        self._dev["n"] = M
        return self._dev

    # ---------------------------------------------------------------- GPU-native path -----
    def _ensure_gpu(self, N: int, device):
        """Allocate the device buffers for the GPU break/form path (bonds ping-pong + scratch)."""
        import warp as wp
        from ffn_sim.dcm.dcm_cadherin_gpu import build_koff_device
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
        from ffn_sim.dcm.dcm_cadherin_gpu import (cad_break_kernel, cad_partner_kernel, cad_form_kernel)
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

    def bonds_now(self) -> np.ndarray:
        """Current (M,2) bond node-pairs — device download in GPU mode, else the host array."""
        if getattr(self, "_n_gpu", None) is not None and self._dev is not None:
            n = int(self._n_gpu)
            if n <= 0:
                return np.zeros((0, 2), np.int64)
            return self._dev["bonds"].numpy()[:n].astype(np.int64)
        return self.bonds
