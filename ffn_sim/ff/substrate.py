r"""Compliant elastic substrate + ligand coat (FF, Warp, µm·pN·s) — the PI experimental axis (PAA gel).

The FF cell adhered onto a RIGID plane (`substrate_plane_kernel`) with FA clutches anchored to FIXED points —
so substrate stiffness E, durotaxis, and ligand density (the PI's Bare/Pre/Lam4 PAA-gel axis) were structurally
impossible. This makes the substrate ELASTIC: each FA-clutch anchor is a movable point tied to the dish by a
Winkler spring ``k_sub`` derived from the substrate Young's modulus E; when the motor-clutch pulls, the anchor
displaces by F/k_sub, so the whole molecular-clutch traction becomes E-dependent (Bangasser-Odde biphasic) and
a stiffness gradient drives durotaxis — all EMERGENT from the fine-grained clutches, no lumped traction law.

KB-anchored (stiffness axis is port-ready): substrate E range **KB-1.5** (0.1–100 kPa biological; PAA 0.1–50;
Phase-1 default **E=5 kPa, ν=0.45**); molecular-clutch constants **KB-PIV-4** (k_on 0.3, k_off 0.1 s⁻¹, F_b 2 pN,
κ_c 0.8 pN/nm, v_u 120 nm/s; optimal substrate stiffness 2–300 kPa); ligand identities **KB-PIV-5/6**
(α6β1 for Lam4, α5β1-FN for Pre; MCF7 lacks α6β4). The absolute magnitude of ``k_sub`` needs the adhesion patch
radius ``a`` — the SAME family as the missing ρ_L ligand-density datum (KB-PIV-5: "α6β1–laminin single-molecule
force kinetics GENUINELY ABSENT"). So ``a`` is a SURFACED parameter (default = nascent-adhesion radius), and only
the RELATIVE stiffness axis (biphasic shape, durotaxis sign, ligand ordering) is claimed as a gate — a magnitude
match awaits KB registration. NOT tuned to pass anything.

E→k_sub bridge: a rigid circular adhesion of radius ``a`` on an elastic half-space (Young E, Poisson ν) has
normal stiffness ``k_sub = 2·E·a/(1−ν²)`` (Boussinesq/Sneddon flat-punch). Rigid limit E→∞ ⇒ k_sub→∞ ⇒ anchor
immovable ⇒ the fixed-anchor path is recovered byte-for-byte (the off-path invariance gate).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import warp as wp

# KB-1.5 substrate modulus (Pa → pN/µm²: 1 Pa = 1 pN/µm²); default physiological PAA gel.
E_SUB_DEFAULT_PA = 5.0e3        # 5 kPa (KB-1.5 Phase-1 default)
NU_SUB_DEFAULT = 0.45          # KB-1.5
# Adhesion patch radius setting the E→k_sub magnitude — SURFACED (same missing-datum family as ρ_L, KB-PIV-5).
A_ADHESION_UM_DEFAULT = 0.05   # nascent adhesion ~50 nm (order-of-magnitude; NOT a fitted constant — flagged)


@dataclass
class SubstrateParams:
    """Elastic substrate in FF units (pN, µm)."""
    E_pn_um2: float             # Young's modulus [pN/µm² = Pa]
    nu: float                   # Poisson ratio
    a_um: float                 # adhesion patch radius [µm] (magnitude lever — surfaced, KB-blocked)
    k_sub: float                # Winkler anchor stiffness [pN/µm] = 2·E·a/(1−ν²)


def resolve_substrate(E_pa: float = E_SUB_DEFAULT_PA, nu: float = NU_SUB_DEFAULT,
                      a_um: float = A_ADHESION_UM_DEFAULT) -> SubstrateParams:
    """Substrate spring from Young's modulus via the flat-punch half-space relation k_sub = 2·E·a/(1−ν²).

    ``E_pa`` in Pa (= pN/µm²). KB-1.5 range 0.1–100 kPa. The E→k_sub magnitude carries the surfaced ``a`` lever;
    the E-DEPENDENCE (∝E, hence the biphasic/durotaxis physics) is exact regardless of ``a``."""
    if not (E_pa > 0 and 0 <= nu < 0.5 and a_um > 0):
        raise ValueError(f"E_pa={E_pa}, nu={nu}, a_um={a_um} out of range")
    k_sub = 2.0 * E_pa * a_um / (1.0 - nu ** 2)
    return SubstrateParams(E_pn_um2=float(E_pa), nu=float(nu), a_um=float(a_um), k_sub=float(k_sub))


@wp.kernel
def substrate_anchor_relax_kernel(anchor: wp.array(dtype=wp.vec3d), rest: wp.array(dtype=wp.vec3d),
                                  bound: wp.array(dtype=wp.int32), clutch_force: wp.array(dtype=wp.vec3d),
                                  k_sub: wp.float64, dt_over_gamma: wp.float64):
    """Overdamped relaxation of each movable substrate anchor (legacy; prefer the closed-form equilibrium below):
    anchor pulled by the clutch reaction ``clutch_force[k]``, restored to ``rest[k]`` by the Winkler spring
    k_sub. Explicit — CFL-bound; can be unstable if the clutch is stiff → use the equilibrium kernel instead."""
    k = wp.tid()
    d = anchor[k] - rest[k]
    f = clutch_force[k] - k_sub * d
    if bound[k] == 0:
        f = -k_sub * d
    anchor[k] = anchor[k] + dt_over_gamma * f


@wp.kernel
def substrate_anchor_equilibrium_kernel(pos: wp.array(dtype=wp.vec3d), ac_idx: wp.array(dtype=wp.int32),
                                        anchor: wp.array(dtype=wp.vec3d), rest: wp.array(dtype=wp.vec3d),
                                        bound: wp.array(dtype=wp.int32), k_int: wp.float64, k_sub: wp.float64):
    """Set each movable anchor to the clutch↔substrate SERIES equilibrium (no explicit iteration → stable):
    ``anch = (k_int·p_actin + k_sub·p_rest)/(k_int + k_sub)``, i.e. the anchor where the clutch spring force
    balances the substrate spring force. The actin then feels the series stiffness k_int·k_sub/(k_int+k_sub) —
    the Bangasser-Odde compliant-substrate physics, quasi-static (the substrate equilibrates within a step).
    Rigid limit k_sub→∞ ⇒ anch→rest (fixed-pin path recovered). Unbound anchors relax to rest."""
    k = wp.tid()
    if bound[k] == 0:
        anchor[k] = rest[k]
        return
    anchor[k] = (k_int * pos[ac_idx[k]] + k_sub * rest[k]) / (k_int + k_sub)


@wp.kernel
def clutch_anchor_reaction_kernel(pos: wp.array(dtype=wp.vec3d), ac_idx: wp.array(dtype=wp.int32),
                                  anchor: wp.array(dtype=wp.vec3d), bound: wp.array(dtype=wp.int32),
                                  k_int: wp.float64, rest: wp.float64, out: wp.array(dtype=wp.vec3d)):
    """The clutch-spring reaction ON each substrate anchor = k_int·(L−rest)·(actin−anchor)/L (Newton pair of the
    force the clutch exerts on the actin node). Drives the movable-anchor relaxation. Unbound → 0."""
    k = wp.tid()
    if bound[k] == 0:
        out[k] = wp.vec3d(0.0, 0.0, 0.0)
        return
    d = pos[ac_idx[k]] - anchor[k]
    L = wp.length(d)
    if L < wp.float64(1e-12):
        out[k] = wp.vec3d(0.0, 0.0, 0.0)
        return
    out[k] = (k_int * (L - rest) / L) * d                 # pulls the anchor toward the loaded actin


def substrate_cfl_dt(k_sub: float, gamma: float, safety: float = 0.1) -> float:
    """Explicit CFL for the anchor Winkler spring: dt < safety·γ/k_sub. (k_sub in the implicit K is the alt.)"""
    return float(safety * gamma / max(k_sub, 1e-30))


__all__ = ["SubstrateParams", "resolve_substrate", "substrate_anchor_relax_kernel",
           "substrate_anchor_equilibrium_kernel", "clutch_anchor_reaction_kernel", "substrate_cfl_dt",
           "E_SUB_DEFAULT_PA", "NU_SUB_DEFAULT", "A_ADHESION_UM_DEFAULT"]
