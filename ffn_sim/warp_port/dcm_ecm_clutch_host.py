"""Host-side explicit integrin-ECM catch-slip clutch for the Warp DCM loop (C6).

Replaces the conservative substrate-WETTING proxy (the lumped in-plane spreading force) with
the real molecular-clutch mechanism (CLAUDE.md hard rule): basal membrane nodes engage fixed
substrate ligand sites through explicit integrin focal-adhesion clutches whose lifetime
follows the faithful Pereverzev two-pathway catch-slip off-rate
``k_off(F)=k_s·exp(F/F_s)+k_c·exp(−F/F_c)`` (``validation.pereverzev``, the SAME law H.4's
IntegrinBondUpdater and cell/dcm_ecm.py's FaCatchSlipUpdater use). A clutch is a harmonic
spring from the basal node to the substrate site it gripped; as the lamellipodium (M2) and
membrane tension pull the node, the clutch transmits that load to the dish as TRACTION, and
catch-slip governs whether it holds (catch, near F*) or releases (slip, F≫F_s). Spreading is
then traction-limited and fully mechanistic — protrusion (M2) + clutch (C6), no wetting body force.

Host-managed at the FA cadence (engage/break, KDTree-free — basal nodes only); a device kernel
applies the clutch force every step. ×40 scale bridge: a node-clutch is the integrin BUNDLE over
a basal patch, so its engage range is the mesh contact scale and k_fa is set so the slip force F_s
is reached at the engage limit; the Pereverzev SHAPE (F*, F_s, F_c) stays molecular.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ffn_sim.validation.pereverzev import pereverzev_k_off, PereverzevParams


@dataclass
class EcmClutchParams:
    """Resolved FA clutch params (dcm_ecm.ResolvedDcmEcm defaults; Pereverzev KU-2.18, F*≈7pN)."""
    k_fa: float = 1.0e-3           # N/m  molecular integrin clutch stiffness (KU-2.7) — bridged below
    fa_capture: float = 0.5e-6     # m    basal-node↔site engage / rebind range (mesh-bridged)
    k_on: float = 1.0              # s⁻¹  FA rebind on-rate (KU-2.18)
    fa_batch_steps: int = 50       # catch-slip updater cadence
    seed: int = 11
    # Pereverzev two-pathway (KU-2.18): k_off = k_s·exp(F/F_s) + k_c·exp(−F/F_c)
    k_off_slip: float = 0.5        # s⁻¹  k_s
    F_s: float = 30.0e-12          # N    slip force scale
    k_off_catch: float = 0.4       # s⁻¹  k_c
    F_c: float = 7.0e-12           # N    catch force scale
    # FA-patch FORCE bridge. A basal node-clutch is one focal-adhesion PATCH (a BUNDLE of ~bundle_n
    # integrins), so it transmits bundle_n × the single-integrin force; its Pereverzev k_off is still
    # evaluated at the PER-INTEGRIN load F/bundle_n (F_s, F_c stay molecular). bundle_n=1 → legacy
    # single-integrin force (back-compat). Set so per-clutch lands in KB-2.12's audited per-FA
    # 1–10 nN (~5 nN; Plotnikov2012, Trichet2012) and per-cell in 10–100 nN.
    bundle_n: float = 1.0
    # C4: substrate ligand-coating density (Bare/Pre/Lam4). More ECM ligand → a basal node finds a
    # binding site faster → higher effective engagement on-rate. Scales k_on (engagement), so the
    # steady-state engaged-clutch count / traction rises with ligand density. ligand_density=1 =
    # baseline (Bare). Per-condition values are EXPERIMENT-anchored to the Lam4>Pre>Bare ordering
    # (the platform's PI MCF7-on-pV4D4/col-I conditions) — NOT tuned to a spreading outcome.
    ligand_density: float = 1.0


class EcmClutchHost:
    """Owns the dynamic integrin-ECM clutch set: basal node ↔ fixed substrate site.

    Usage (ecm-clutch mode; wetting OFF):
        ecm = EcmClutchHost(cof=cof_a, n_cells=..., z0=z0, R=R, dt=dt, c_adh=c_adh)
        if s % ecm.fa_batch_steps == 0:
            ecm.update(pos_d.numpy()); ecm.upload(device)
        # every step: ecm_clutch_force_kernel(node_idx_d, anchor_d, n, k_fa, pos_d, force_d)
    """

    def __init__(self, *, cof: np.ndarray, n_cells: int, z0: float, R: float, dt: float,
                 c_adh: float, params: EcmClutchParams | None = None):
        self.p = params or EcmClutchParams()
        self.cof = np.asarray(cof, dtype=np.int64)
        self.n_cells = int(n_cells)
        self.z0 = float(z0)
        self.R = float(R)
        self.dt = float(dt)
        self.fa_batch_steps = self.p.fa_batch_steps
        self.dt_batch = self.fa_batch_steps * self.dt
        self._rng = np.random.default_rng(self.p.seed)
        self._pp = PereverzevParams(k_s=self.p.k_off_slip, F_s=self.p.F_s,
                                    k_c=self.p.k_off_catch, F_c=self.p.F_c)
        # ×40 scale bridge: a basal node = integrin BUNDLE over a patch. Engage range = the mesh
        # contact scale (c_adh); k_fa set so the slip force F_s is reached at the engage limit, so
        # a clutch dragged across its capture window enters slip. Pereverzev shape stays molecular.
        self.engage_range = float(c_adh)
        self.basal_band = self.engage_range            # node within this of z0 = basal
        self.k_fa = self.p.F_s / max(self.engage_range, 1e-12)
        # force→k_off table over [0, ~5·F_s] (analytic, but tabulate to skip the overflow guard)
        self._fs = np.linspace(0.0, 5.0 * self.p.F_s, 256)
        self._koff = pereverzev_k_off(self._fs, self._pp)
        # clutch set: node index + anchor (3,) fixed substrate site
        self.node_idx = np.zeros(0, dtype=np.int64)
        self.anchor = np.zeros((0, 3), dtype=np.float64)
        self.n_engaged = 0
        self.n_slipped = 0
        self._dev = None

    def _koff_of(self, F):
        return np.interp(F, self._fs, self._koff)

    def update(self, P: np.ndarray) -> None:
        """One FA tick: catch-slip BREAK of loaded clutches, then ENGAGE slack basal nodes."""
        P = np.asarray(P, dtype=np.float64)
        cof = self.cof
        # --- BREAK (Pereverzev catch-slip on the clutch load) ---
        if self.node_idx.size:
            ni = self.node_idx
            live = cof[ni] >= 0
            ni = ni[live]; anc = self.anchor[live]
            F = self.k_fa * np.linalg.norm(P[ni] - anc, axis=1)
            p_break = 1.0 - np.exp(-self._koff_of(F) * self.dt_batch)
            keep = self._rng.random(ni.size) >= p_break
            self.n_slipped += int((~keep).sum())
            self.node_idx = ni[keep]; self.anchor = anc[keep]
        # --- ENGAGE (basal, unbonded nodes grip the dish at their xy, z=z0) ---
        engaged = np.zeros(P.shape[0], dtype=bool)
        if self.node_idx.size:
            engaged[self.node_idx] = True
        basal = (cof >= 0) & (np.abs(P[:, 2] - self.z0) <= self.basal_band) & (~engaged)
        cand = np.flatnonzero(basal)
        if cand.size:
            # C4: engagement on-rate ∝ ligand density (more ligand sites → faster binding)
            p_on = 1.0 - np.exp(-self.p.k_on * self.p.ligand_density * self.dt_batch)
            fire = cand[self._rng.random(cand.size) < p_on]
            if fire.size:
                sites = P[fire].copy()
                sites[:, 2] = self.z0                  # the dish ligand it gripped (z = z0)
                self.node_idx = np.concatenate([self.node_idx, fire])
                self.anchor = np.concatenate([self.anchor, sites], axis=0)
                self.n_engaged += int(fire.size)

    @property
    def n_clutches(self) -> int:
        return int(self.node_idx.size)

    def upload(self, device):
        import warp as wp
        M = self.n_clutches
        if self._dev is None or self._dev["cap"] < M:
            cap = max(M, 1, (self._dev["cap"] * 2 if self._dev else 0))
            self._dev = {"cap": cap,
                         "node": wp.zeros(cap, dtype=wp.int32, device=device),
                         "anchor": wp.zeros(cap, dtype=wp.vec3d, device=device)}
        if M:
            self._dev["node"].assign(self.node_idx.astype(np.int32))
            self._dev["anchor"].assign(np.ascontiguousarray(self.anchor))
        self._dev["n"] = M
        return self._dev
