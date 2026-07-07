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
    """Overdamped relaxation of each movable substrate anchor: the anchor is pulled by the clutch reaction
    ``clutch_force[k]`` and restored to its dish rest point ``rest[k]`` by the Winkler spring k_sub. Equilibrium
    displacement = F_clutch / k_sub — the substrate compliance the clutch feels. Rigid limit k_sub→∞ pins the
    anchor at rest (byte-identical to the fixed-anchor path). Bound anchors only (unbound relax to rest)."""
    k = wp.tid()
    d = anchor[k] - rest[k]
    f = clutch_force[k] - k_sub * d                       # net force on the anchor DOF
    if bound[k] == 0:
        f = -k_sub * d                                    # unbound: just relax to rest
    anchor[k] = anchor[k] + dt_over_gamma * f


def substrate_cfl_dt(k_sub: float, gamma: float, safety: float = 0.1) -> float:
    """Explicit CFL for the anchor Winkler spring: dt < safety·γ/k_sub. (k_sub in the implicit K is the alt.)"""
    return float(safety * gamma / max(k_sub, 1e-30))


__all__ = ["SubstrateParams", "resolve_substrate", "substrate_anchor_relax_kernel", "substrate_cfl_dt",
           "E_SUB_DEFAULT_PA", "NU_SUB_DEFAULT", "A_ADHESION_UM_DEFAULT"]
