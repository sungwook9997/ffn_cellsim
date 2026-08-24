"""T-3 EDGE-ADVANCING CLUTCH TRACTION — the pancake flattening driver (PI 2026-07-13).

The native full-compartment result (2026-07-13, `native_spread_*`): T-1 (substrate adhesion) +
T-2 (active area generation) alone leave the DCM cell a rounded WETTING CAP (dome, AR 1.0–1.8),
NOT a pancake — because nothing organises the T-2 membrane-area budget into a flat footprint.
Real cells spread PAST their passive wetting angle by active lamellipodial protrusion + molecular
clutch: the leading edge advances along the substrate and traction pulls the rim outward.

THE TERM (this module): the basal-rim LEADING nodes are pulled radially OUTWARD in the substrate
plane (z0) toward clutch-anchored actin sites that sit AHEAD of the current rim and ADVANCE as the
cell spreads (the ratchet). This converts protrusion into net edge advance and HOLDS the flattened
footprint — driving the cell from the wetting cap into a thin pancake (a fried-egg lamella). It is
mechanistic, not a body-force proxy: the tether is a spring to a clutch-gripped substrate site
whose reaction is the rigid dish (`lamellipodium_tether_accum`, the committed parity kernel).

Composition: reuses the committed ``gather_lead_pos`` + ``lamellipodium_tether_accum`` device
kernels; this host re-detects the rim + advances the actin ring at low cadence (the clutch ratchet).

GROUNDING (no magic numbers):
  * traction_cap (per leading node) = per-FA traction ~5 nN (KB-2.12: per-FA 1–10 nN). Total edge
    traction = n_lead × cap ~ 10–100 nN/cell (KB-2.12 per-cell). Clutch-limited (capped), so the
    magnitude is the grounded quantity; ``k_tether`` need only be stiff enough to reach the cap.
  * advance_gap = lamellipodial reach ahead of the rim ~1–2 µm (KB-3.6 lamellipodium ~µm).
  * r_max_spread = max footprint radius, tied to the T-2 reservoir so area-budget (T-2) and
    edge-advance (T-3) are consistent: R_max = sqrt(2·reservoir)·R (thin-pancake A/A0 = 2·reservoir).
  * lead_frac / basal_band = which basal nodes are the advancing lamellipodial rim (geometry).

Sanity Gate:
  * Dimensional: traction_cap [N], k_tether [N/m], advance_gap/r_max [m]. ✓
  * Bounded advance: the actin ring never passes r_max_spread → the footprint is reservoir-capped
    (no runaway spread). ✓
  * Outward-only: anchors are placed OUTWARD of the rim (radius = rim+gap) and the tether only
    pulls toward anchors with larger radial projection → net OUTWARD edge advance. ✓
  * Off-invariance: edge_traction disabled ⇒ no leading nodes uploaded ⇒ kernel is a no-op. ✓
"""

from __future__ import annotations

import numpy as np

import warp as wp

wp.init()

#: KB-2.12 per-FA traction (~5 nN) — the per-leading-node clutch traction cap.
TRACTION_CAP_N: float = 5.0e-9


class EdgeTractionHost:
    """Host-managed advancing clutch-traction rim (low cadence) + device arrays for the tether.

    Detects the basal-rim leading nodes, places an outward-advancing actin ring in the z0 plane,
    and exposes the device arrays consumed by ``gather_lead_pos`` + ``lamellipodium_tether_accum``.
    Call :meth:`update` at the ratchet cadence (host), :meth:`upload` after each update.
    """

    def __init__(
        self,
        *,
        z0: float,
        R: float,
        reservoir: float = 1.4,
        k_tether: float = 1.0e-2,          # stiff → saturates to the clutch-limited cap
        traction_cap: float = TRACTION_CAP_N,
        total_traction: float | None = None,   # GRID-INVARIANT: per-node cap = total/n_lead (KB-2.12
                                               # per-cell total ~100 nN; Gil-Redondo MCF7 102 nN)
        lead_frac: float = 0.7,
        basal_band: float = 2.0e-6,
        advance_gap: float = 1.5e-6,       # lamellipodial reach ahead of the rim
        n_actin: int = 32,
        r_max_spread: float | None = None,
    ):
        self.z0 = float(z0)
        self.R = float(R)
        self.k_tether = float(k_tether)
        self.traction_cap = float(traction_cap)
        self.total_traction = (None if total_traction is None else float(total_traction))
        self.lead_frac = float(lead_frac)
        self.basal_band = float(basal_band)
        self.advance_gap = float(advance_gap)
        self.n_actin = int(n_actin)
        # footprint cap consistent with the T-2 membrane-area reservoir (thin-pancake A/A0=2·res)
        self.r_max_spread = float(r_max_spread if r_max_spread is not None
                                  else np.sqrt(2.0 * reservoir) * R)
        self.tether_radius = 3.0 * R
        self._dev: dict | None = None

    def update(self, pos: np.ndarray, cof: np.ndarray) -> None:
        """Re-detect the basal-rim leading nodes + advance the actin ring (host, low cadence)."""
        active = cof >= 0
        idx_all = np.where(active)[0]
        if idx_all.size == 0:
            self._geom = None
            return
        p = pos[idx_all]
        centroid = p.mean(axis=0)
        rel = pos - centroid
        rxy = np.linalg.norm(rel[:, :2], axis=1)
        basal = active & (np.abs(pos[:, 2] - self.z0) <= self.basal_band)
        rmax = float(rxy[basal].max()) if basal.any() else float(rxy[active].max())
        lead = basal & (rxy >= self.lead_frac * max(rmax, 1e-18))
        li = np.where(lead)[0]
        if li.size == 0:                                  # no rim yet → fall back to lowest basal ring
            li = np.where(basal)[0]
        if li.size == 0:
            self._geom = None
            return
        ox = np.where(rxy[li] > 0, rel[li, 0] / np.maximum(rxy[li], 1e-18), 0.0)
        oy = np.where(rxy[li] > 0, rel[li, 1] / np.maximum(rxy[li], 1e-18), 0.0)
        proj = rel[li, 0] * ox + rel[li, 1] * oy
        # advancing actin ring: OUTWARD of the current rim by advance_gap, capped at r_max_spread
        r_ring = min(rmax + self.advance_gap, self.r_max_spread)
        ang = np.linspace(0.0, 2.0 * np.pi, self.n_actin, endpoint=False)
        actin = np.stack([centroid[0] + r_ring * np.cos(ang),
                          centroid[1] + r_ring * np.sin(ang),
                          np.full(self.n_actin, self.z0)], axis=1)
        # GRID-INVARIANT traction: distribute the grounded TOTAL edge traction over the detected
        # leading nodes → per-node cap = total/n_lead (so the total is subdiv-independent; a finer
        # mesh has more, weaker leading nodes summing to the same physiological cell traction).
        if self.total_traction is not None:
            self.traction_cap = self.total_traction / max(int(li.size), 1)
        self._geom = dict(
            lead_idx=li.astype(np.int32),
            ccx=np.full(li.size, centroid[0]), ccy=np.full(li.size, centroid[1]),
            ox=ox, oy=oy, proj=proj, actin=actin, r_ring=r_ring, rmax=rmax, n_lead=int(li.size),
        )

    def upload(self, device) -> dict:
        """Push the current rim geometry to device arrays for the tether kernel."""
        g = getattr(self, "_geom", None)
        if g is None:
            self._dev = {"n_lead": 0}
            return self._dev
        L = g["n_lead"]
        self._dev = {
            "n_lead": L,
            "lead_idx": wp.array(g["lead_idx"], dtype=wp.int32, device=device),
            "lead_rp": wp.zeros(L, dtype=wp.vec3d, device=device),
            "ccx": wp.array(np.ascontiguousarray(g["ccx"]), dtype=wp.float64, device=device),
            "ccy": wp.array(np.ascontiguousarray(g["ccy"]), dtype=wp.float64, device=device),
            "ox": wp.array(np.ascontiguousarray(g["ox"]), dtype=wp.float64, device=device),
            "oy": wp.array(np.ascontiguousarray(g["oy"]), dtype=wp.float64, device=device),
            "proj": wp.array(np.ascontiguousarray(g["proj"]), dtype=wp.float64, device=device),
            "actin": wp.array(np.ascontiguousarray(g["actin"]), dtype=wp.vec3d, device=device),
            "r_ring": g["r_ring"], "rmax": g["rmax"],
        }
        return self._dev


if __name__ == "__main__":
    import sys
    print("=" * 68)
    print("T-3 edge-advancing clutch traction — self-test / Sanity Gate")
    print("=" * 68)
    R = 7.5e-6
    # a hemispherical cap resting on z0 (dome, like the native T-1+T-2 result)
    rng = np.random.default_rng(0)
    th = rng.uniform(0, np.pi / 2, 200); ph = rng.uniform(0, 2 * np.pi, 200)
    pos = np.stack([R * np.sin(th) * np.cos(ph), R * np.sin(th) * np.sin(ph),
                    R * (1 - np.cos(th))], axis=1)  # bottom at z=0
    cof = np.zeros(200, dtype=np.int64)
    h = EdgeTractionHost(z0=0.0, R=R, reservoir=2.0)
    h.update(pos, cof)
    g = h._geom
    print(f"\n[grounding] traction_cap={h.traction_cap:.1e} N (per-FA ~5nN, KB-2.12); "
          f"r_max_spread={h.r_max_spread*1e6:.1f}µm (sqrt(2·2.0)·R); advance_gap={h.advance_gap*1e6:.1f}µm")
    print(f"[A] rim detected: n_lead={g['n_lead']} of 200; rmax={g['rmax']*1e6:.2f}µm; "
          f"actin ring at r={g['r_ring']*1e6:.2f}µm (OUTWARD of rim: {g['r_ring']>g['rmax']})")
    # outward check: ring radius > rim radius, and bounded by r_max
    outward = g["r_ring"] > g["rmax"]
    bounded = g["r_ring"] <= h.r_max_spread + 1e-12
    # leading nodes are basal (low z) + peripheral
    lead_z = pos[g["lead_idx"], 2]
    basal_ok = bool(lead_z.max() <= h.basal_band + 1e-9)
    dev = h.upload("cpu")
    upload_ok = dev["n_lead"] == g["n_lead"]
    # off-invariance: empty geom → no-op device
    h2 = EdgeTractionHost(z0=0.0, R=R); h2._geom = None
    off_ok = h2.upload("cpu")["n_lead"] == 0
    ok = outward and bounded and basal_ok and upload_ok and off_ok
    print(f"[B] outward advance={outward}  reservoir-bounded={bounded}  leading-are-basal={basal_ok}")
    print(f"[C] upload={upload_ok}  off-invariance(no rim→no-op)={off_ok}")
    print("\n" + "=" * 68)
    print(f"SELF-TEST {'PASS' if ok else 'FAIL'}")
    print("=" * 68)
    sys.exit(0 if ok else 1)
