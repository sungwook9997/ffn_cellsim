r"""Focal-adhesion maturation (FF, Warp, µm·pN·s) — talin unfolding + vinculin reinforcement + Hill FA growth.

The live FF integrin clutch (``fa_clutch_warp``: Chan-Odde binding + Pereverzev catch-slip) is a bare bond —
the MECHANOSENSOR layer that turns a nascent adhesion into a mature, reinforced, force-growing focal adhesion
was absent. This adds it, all KB-2.x anchored, as per-clutch/per-FA state on top of the existing clutch:

  1. TALIN unfolding (KB-2.6): a rod domain under tension unfolds at ``k_unfold(F)=k_u0·exp(F·dx/kBT)`` (R3, the
     lowest: F_th≈5 pN, k_u0≈3e-4 s⁻¹, dx≈0.5 nm), exposing a vinculin-binding site. Force-gated mechanosensor.
  2. VINCULIN recruitment (KB-2.7): ``dN_vin/dt = k_rec·p_unf·(N_max−N_vin) − k_diss·N_vin`` reinforces the
     clutch stiffness ``k_int_eff = k_int·(1 + α·N_vin)`` (α≈0.01–0.05). [k_diss is DOI-unreconciled → SURFACED.]
  3. FA GROWTH (KB-2.17): the adhesion area grows Hill-wise under force ``dA/dt = (k_g(F) − k_d)·A``,
     ``k_g(F) = k_g0·Fⁿ/(Fⁿ+F_thⁿ)`` (n≈3, k_g0≈0.01, k_d≈0.005 s⁻¹, F_th≈5 pN) — below F_th the FA disassembles.

Emergent biphasic-stiffness traction + durotaxis come from this maturation running on the compliant substrate
(``ff/substrate.py``). No lumped mechanosensing law — force → unfolding → binding → stiffness/area, molecule by
molecule. No magic numbers: every constant is a KB-2.x value; the ONLY gap (vinculin k_diss) is flagged.

Verification gates (analytic): (1) talin Bell parity — kernel ``k_unfold(5 pN)`` equals the closed form
``k_u0·exp(5·dx/kBT)``; (2) Hill FIXED POINT — dA/dt = 0 at F = F_th = 5.0 pN exactly (k_g(F_th)=k_g0/2 must
equal k_d ⇒ the KB constraint k_g0=2·k_d is self-consistent), growing above and disassembling below.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

from aleph.laws.units import KBT

# KB-2.6 talin (R3, the lowest-force rod domain — the default one-state mechanosensor)
TALIN_KU0 = 3.0e-4              # s⁻¹   unfolding rate at zero force
TALIN_DX_UM = 5.0e-4           # µm     unfolding distance 0.5 nm
TALIN_FTH_PN = 5.0             # pN     characteristic unfolding force
TALIN_KREFOLD = 1.0e-2         # s⁻¹    refolding rate at zero load (slow; load suppresses folding)
# KB-2.7 vinculin
VIN_KREC = 0.03                # s⁻¹    recruitment rate to an unfolded talin
VIN_ALPHA = 0.03               # dimensionless  stiffness gain per vinculin (band 0.01–0.05)
VIN_NMAX = 5.0                 # vinculin sites per talin rod (R1–R13 carry VBS; working default)
VIN_KDISS = 0.03               # s⁻¹    ⚠ SURFACED: KB-2.7 gives no k_diss (Han2021/TapiaRojo2019 NO_DOI) —
#                                        placeholder = k_rec (⇒ steady N_vin≈N_max·p/2); NOT a validated value.
# KB-2.17 FA growth (Hill)
FA_KG0 = 0.01                  # s⁻¹    max growth rate
FA_KD = 0.005                  # s⁻¹    disassembly rate  (KB constraint: k_g0 = 2·k_d ⇒ fixed point at F_th)
FA_N = 3.0                     # Hill coefficient (band 2–4)
FA_FTH_PN = 5.0                # pN     per-clutch growth threshold


@dataclass
class MaturationParams:
    ku0: float = TALIN_KU0; dx_um: float = TALIN_DX_UM; k_refold: float = TALIN_KREFOLD
    k_rec: float = VIN_KREC; alpha: float = VIN_ALPHA; n_max: float = VIN_NMAX; k_diss: float = VIN_KDISS
    kg0: float = FA_KG0; kd: float = FA_KD; n_hill: float = FA_N; fth: float = FA_FTH_PN
    kT: float = KBT


def talin_unfold_rate(F_pn: float, p: MaturationParams | None = None) -> float:
    """Bell unfolding rate ``k_u0·exp(F·dx/kBT)`` [s⁻¹] — the talin mechanosensor closed form (KB-2.6 gate)."""
    p = p or MaturationParams()
    return float(p.ku0 * np.exp(F_pn * p.dx_um / p.kT))


def fa_growth_rate(F_pn, p: MaturationParams | None = None):
    """Net FA area growth rate ``k_g(F) − k_d`` [s⁻¹] (KB-2.17); zero (fixed point) at F = F_th."""
    p = p or MaturationParams()
    F = np.asarray(F_pn, float)
    kg = p.kg0 * F ** p.n_hill / (F ** p.n_hill + p.fth ** p.n_hill)
    return kg - p.kd


@wp.kernel
def talin_vinculin_kernel(load: wp.array(dtype=wp.float64), p_unf: wp.array(dtype=wp.float64),
                          n_vin: wp.array(dtype=wp.float64), bound: wp.array(dtype=wp.int32),
                          ku0: wp.float64, dx: wp.float64, kT: wp.float64, k_refold: wp.float64,
                          k_rec: wp.float64, k_diss: wp.float64, n_max: wp.float64, dt: wp.float64):
    """Advance per-clutch talin-unfolded fraction ``p_unf`` (Bell forward k_unfold(F)·(1−p) − k_refold·p) and
    vinculin count ``n_vin`` (k_rec·p·(N_max−N_vin) − k_diss·N_vin). ``load[k]`` = current clutch tension [pN].
    Unbound clutches refold + shed vinculin (no load)."""
    k = wp.tid()
    F = load[k]
    if bound[k] == 0:
        F = wp.float64(0.0)
    kf = ku0 * wp.exp(F * dx / kT)                         # Bell unfolding rate — force ENHANCES (del Rio 2009)
    kr = k_refold * wp.exp(-F * dx / kT)                   # refolding — force SUPPRESSES (Yao 2016; two-state Bell)
    p = p_unf[k]
    p_unf[k] = wp.min(wp.max(p + dt * (kf * (wp.float64(1.0) - p) - kr * p), wp.float64(0.0)), wp.float64(1.0))
    nv = n_vin[k]
    n_vin[k] = wp.max(nv + dt * (k_rec * p_unf[k] * (n_max - nv) - k_diss * nv), wp.float64(0.0))


@wp.kernel
def clutch_load_kernel(pos: wp.array(dtype=wp.vec3d), ac_idx: wp.array(dtype=wp.int32),
                       anchor: wp.array(dtype=wp.vec3d), bound: wp.array(dtype=wp.int32),
                       k_int: wp.float64, rest: wp.float64, load: wp.array(dtype=wp.float64)):
    """Per-clutch tension magnitude ``F = k_int·|L − rest|`` [pN] — the mechanosensor input for talin/FA growth.
    Unbound → 0."""
    k = wp.tid()
    if bound[k] == 0:
        load[k] = wp.float64(0.0)
        return
    L = wp.length(pos[ac_idx[k]] - anchor[k])
    load[k] = k_int * wp.abs(L - rest)


@wp.kernel
def fa_disassemble_kernel(area: wp.array(dtype=wp.float64), bound: wp.array(dtype=wp.int32), a_min: wp.float64):
    """Load-dependent adhesion: an FA whose area has shrunk below ``a_min`` (sustained sub-threshold force,
    KB-2.17) DISASSEMBLES — its clutch unbinds. The mechanosensor that makes adhesion force-gated."""
    k = wp.tid()
    if bound[k] == 1 and area[k] < a_min:
        bound[k] = 0
        area[k] = wp.float64(1.0)                          # reset for a future nascent re-bind


@wp.kernel
def fa_growth_kernel(load: wp.array(dtype=wp.float64), area: wp.array(dtype=wp.float64),
                     kg0: wp.float64, kd: wp.float64, n_hill: wp.float64, fth: wp.float64, dt: wp.float64):
    """Hill FA-area growth (KB-2.17): ``dA/dt = (k_g0·Fⁿ/(Fⁿ+F_thⁿ) − k_d)·A``; grows above F_th, disassembles
    below. ``load`` = per-FA force [pN]. Area floored at 0 (a disassembled FA)."""
    k = wp.tid()
    F = load[k]
    Fn = wp.pow(F, n_hill)
    kg = kg0 * Fn / (Fn + wp.pow(fth, n_hill))
    area[k] = wp.max(area[k] * (wp.float64(1.0) + dt * (kg - kd)), wp.float64(0.0))


def k_int_effective(k_int_bare: float, n_vin, p: MaturationParams | None = None):
    """Vinculin-reinforced clutch stiffness ``k_int·(1 + α·N_vin)`` [pN/µm] (KB-2.7)."""
    p = p or MaturationParams()
    return k_int_bare * (1.0 + p.alpha * np.asarray(n_vin, float))


__all__ = ["MaturationParams", "talin_unfold_rate", "fa_growth_rate", "talin_vinculin_kernel",
           "fa_growth_kernel", "clutch_load_kernel", "fa_disassemble_kernel", "k_int_effective"]
