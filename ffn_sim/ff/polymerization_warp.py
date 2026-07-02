"""Actin barbed-end polymerization ratchet (FF, Warp) — piece 1/5 of active cell movement (PI 2026-07-02).

The protrusion ENGINE the FF engine lacked (FF_ACTIVE_MOVEMENT_ASSESSMENT): a load-dependent barbed-end
ELONGATION that pushes the leading edge forward. Mechanism = the Mogilner-Oster Brownian/elastic ratchet
(the standard modeling form, already the FF-corpus convention: Mueller 2017, Betorz 2023, Popov 2016):

    v(f) = v0 · exp(−f·δ / kBT)          [µm/s]   (f = compressive load on the barbed end, ≥0)
    v0   = (k_on·M − k_off) · δ_full               (Pollard 1986 barbed-end kinetics)

Each step the barbed-end SEGMENT rest length grows by v(f)·dt (kinematic elongation); the reshape then
advances the barbed node → protrusion. Load-dependence is the ratchet: a heavier load slows growth,
stalling as f → (kBT/δ)·ln(v0/v).

**Lit-anchored constants (DOIs verified; NOT invented — flagged for PI KB-registration):**
- k_on = 11.6 µM⁻¹s⁻¹, k_off = 1.4 s⁻¹ (barbed end, ATP-actin) — Pollard 1986 JCB 10.1083/jcb.103.6.2747
  (THE canonical actin kinetics; restated in Pollard-Borisy 2003).
- δ_full = 2.7 nm axial rise/subunit; ratchet δ = δ_full/2 = 1.35 nm (double-helix half-monomer) — Mogilner-
  Oster 1996 10.1016/S0006-3495(96)79496-1.
- kBT = 4.28 pN·nm at 310 K (physiological body temp, per the physiological-baseline rule).
- Force-velocity form = exponential BR (Mogilner-Oster 2003 10.1016/S0006-3495(03)74969-8, in corpus).

**Analytic ground truth (oracle-is-crosscheck rule):** the always-on validation is (1) v(f) = v0·exp(−fδ/kBT)
exactly, and (2) the e-fold force kBT/δ ≈ 3.2 pN. Footer 2007 (~1 pN free single-filament stall) / M-O (~5-7 pN
anchored) are context-dependent cross-checks, NOT hard-coded (a free filament buckles — real protrusion needs
bundling/branching, and the membrane load is SHARED across N barbed ends → per-filament f = f_load/N).

Units: FF is µm, pN, s; kBT in pN·µm, δ in µm, v0 in µm/s.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import warp as wp

wp.init()

KBT_310K_PN_UM = 4.28e-3          # k_B·T at 310 K [pN·µm] (=4.28 pN·nm)
ACTIN_RISE_UM = 2.7e-3            # δ_full monomer axial rise [µm] (kim_network ACTIN_RISE_UM)
K_ON_BARBED = 11.6               # µM⁻¹ s⁻¹ (Pollard 1986)
K_OFF_BARBED = 1.4               # s⁻¹     (Pollard 1986)


@dataclass(frozen=True, slots=True)
class ResolvedPolymerization:
    """Barbed-end ratchet params in FF units (µm, pN, s)."""
    v0_um_s: float       # unloaded elongation velocity (k_on·M − k_off)·δ_full
    delta_um: float      # ratchet step (half-monomer) 1.35 nm
    kBT_pN_um: float     # 4.28e-3 at 310 K
    efold_force_pN: float  # kBT/δ ≈ 3.2 pN — the analytic ground-truth e-fold force
    G_actin_uM: float


def resolve_polymerization(*, G_actin_uM: float = 20.0, k_on: float = K_ON_BARBED, k_off: float = K_OFF_BARBED,
                           delta_full_um: float = ACTIN_RISE_UM, kBT: float = KBT_310K_PN_UM) -> ResolvedPolymerization:
    """Resolve the ratchet from Pollard barbed-end kinetics at local free G-actin M [µM] (physiological ~10-40;
    default 20 → v0 ≈ 0.62 µm/s). δ = δ_full/2 (ratchet half-monomer). Grid/lit-anchored; no tuning."""
    if not (1.0 <= G_actin_uM <= 100.0):
        raise ValueError(f"G_actin {G_actin_uM} µM outside physiological band [1,100]")
    v0 = (k_on * G_actin_uM - k_off) * delta_full_um       # µm/s
    if v0 <= 0.0:
        raise ValueError(f"v0 ≤ 0 (G-actin {G_actin_uM} below the critical concentration)")
    delta = 0.5 * delta_full_um
    return ResolvedPolymerization(v0_um_s=v0, delta_um=delta, kBT_pN_um=kBT,
                                  efold_force_pN=kBT / delta, G_actin_uM=G_actin_uM)


@wp.func
def _ratchet_velocity(f_comp: wp.float64, v0: wp.float64, delta: wp.float64, kBT: wp.float64) -> wp.float64:
    """Mogilner-Oster barbed-end velocity v(f) = v0·exp(−f·δ/kBT); f_comp is the compressive load (≥0)."""
    fc = f_comp
    if fc < wp.float64(0.0):
        fc = wp.float64(0.0)                              # tension does not speed polymerization past v0
    return v0 * wp.exp(-fc * delta / kBT)


@wp.kernel
def polymerization_kernel(
    pos: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),          # current force on each node [pN] (load felt at the barbed end)
    barbed_node: wp.array(dtype=wp.int32),    # (B,) barbed-end node index of each growing filament
    prev_node: wp.array(dtype=wp.int32),      # (B,) the node one in from the barbed end (defines the axis)
    barbed_seg: wp.array(dtype=wp.int32),     # (B,) the seg_rest index of the barbed-end segment
    load_share: wp.array(dtype=wp.float64),   # (B,) per-filament load divisor N (barbed-end density sharing)
    seg_rest: wp.array(dtype=wp.float64),     # (S,) segment rest lengths [µm] (grown in place)
    v0: wp.float64, delta: wp.float64, kBT: wp.float64, dt: wp.float64,
    v_out: wp.array(dtype=wp.float64),        # (B,) the realized v(f) [µm/s] (diagnostic)
):
    """Grow each barbed-end segment by v(f)·dt (load-dependent ratchet). f = compressive load along the
    OUTWARD barbed axis, divided by the shared barbed-end count (f_load/N). Kinematic elongation → the
    reshape advances the barbed node → protrusion."""
    b = wp.tid()
    nb = barbed_node[b]
    npv = prev_node[b]
    ax = pos[nb] - pos[npv]                                # outward barbed axis (unnormalised)
    L = wp.length(ax)
    if L < wp.float64(1e-12):
        v_out[b] = wp.float64(0.0)
        return
    axn = ax / L
    # compressive load on the barbed end = force pushing it INWARD (−axis), per shared barbed end
    f_along = -wp.dot(force[nb], axn)                      # >0 ⇒ compressive (opposes growth)
    share = load_share[b]
    if share < wp.float64(1.0):
        share = wp.float64(1.0)
    f_comp = f_along / share
    v = _ratchet_velocity(f_comp, v0, delta, kBT)
    v_out[b] = v
    seg_rest[barbed_seg[b]] = seg_rest[barbed_seg[b]] + v * dt


def ratchet_velocity_np(f_comp, spec: ResolvedPolymerization):
    """Host reference for the analytic ground truth v(f)=v0·exp(−fδ/kBT) (validation)."""
    f = np.maximum(np.asarray(f_comp, dtype=np.float64), 0.0)
    return spec.v0_um_s * np.exp(-f * spec.delta_um / spec.kBT_pN_um)
