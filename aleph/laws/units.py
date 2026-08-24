"""FF engine unit system — nondimensionalize to **pN · µm · s** (Stage 6c-a).

WHY (the scale finding, ENGINE.md §2). In SI a single actin filament's bending force is ~1e-12 N
and its energy ~1e-21 J. The reused implicit solver (``dcm.dcm_warp_implicit``) carries absolute
thresholds calibrated for DCM whole-cell forces (newton_tol=1e-10, cg_tol=1e-8); against
SI-scale FF forces the residual ‖G‖ ~1e-12 is already below newton_tol on entry → instant FALSE
convergence (Δx=0, no relaxation). The fix is to carry the FF layer in a unit system where the
characteristic force is O(1).

THE UNIT SYSTEM (chosen so single-filament forces are O(0.1–100), well-conditioned):

    length   L*  = 1 µm      = 1e-6  m
    time     t*  = 1 s       = 1      s
    force    F*  = 1 pN      = 1e-12 N
    energy   E*  = 1 pN·µm   = 1e-18 J          (= F*·L*)
    viscosity    1 Pa·s      = 1 pN·s·µm⁻²       (N·s·m⁻² → identity in these units, see below)

Conversions (SI → FF):  x[µm]=x[m]/1e-6 ; F[pN]=F[N]/1e-12 ; E[pN·µm]=E[J]/1e-18 ;
κ[pN·µm²]=κ[N·m²]/1e-24 (E·I = F·L²) ; η[pN·s·µm⁻²]=η[Pa·s]·(1e12 pN/N)·(1e-6 m/µm)² = η[Pa·s]·1
(so Pa·s is numerically unchanged — a convenient coincidence of this system).

GROUNDED CONSTANTS (literature inputs — none invented; see refs in each entry).

The force kernels (``forces_warp``) and the constraint projector (``constraints``) are unit-AGNOSTIC
(pure arithmetic on the arrays), so nondimensionalizing is purely a matter of *building the network
in these units* and feeding the matching γ/κ — no kernel change. This module is the single source of
the convention + the grounded physical constants + the SI↔FF converters + the NF2007 mobility law.
"""

from __future__ import annotations

import numpy as np

# ── SI base scales of the FF unit system ────────────────────────────────────────────────────────
UM = 1e-6          # m  per FF length unit
PN = 1e-12         # N  per FF force unit
S = 1.0            # s  per FF time unit
PN_UM = PN * UM    # J  per FF energy unit (= 1e-18 J)
PN_UM2 = PN * UM * UM  # N·m² per FF bending-modulus unit (= 1e-24 N·m²)

# ── Physical constants ──────────────────────────────────────────────────────────────────────────
KB_SI = 1.380649e-23           # J/K  (CODATA)
T_PHYSIO_K = 310.0             # K  (37 °C, physiological)
KBT_SI = KB_SI * T_PHYSIO_K    # J   ≈ 4.28e-21
KBT = KBT_SI / PN_UM           # pN·µm  ≈ 4.28e-3  (thermal energy in FF units)

# Actin filament (the cortex/network fiber).
LP_ACTIN_UM = 17.0             # µm   persistence length (Gittes et al. 1993, J. Cell Biol. 120:923; ~17.7±0.9)
ACTIN_DIAMETER_UM = 0.008      # µm   filament diameter ≈ 8 nm (standard; e.g. Gittes 1993)
KAPPA_ACTIN = KBT * LP_ACTIN_UM   # pN·µm²  bending modulus κ = k_B·T·L_p  (≈ 0.0728)

# Microtubule (the value NF2007 itself supplies, p9 — a cross-check anchor, not used in the cortex).
KAPPA_MT = 20.0               # pN·µm²  (Nédélec & Foethke 2007, NJP 9:427, p9)

# Cytoplasm viscosity (physiological-baseline HARD rule; same datum the DCM layer uses).
ETA_CYTOPLASM = 65.9          # pN·s·µm⁻² (= 65.9 Pa·s, MCF7 cytoplasm, Dessard 2024)

# Solvent (interstitial water/cytosol) viscosity — the SKELETON-ONLY drag medium once the Biot pore
# fluid supplies the Darcy resistance explicitly (FSI wiring, DYNAMIC_FEM_CFD_UPGRADE_ROADMAP 2026-07-16).
# The measured effective cytoplasm η=65.9 Pa·s IS the poroelastic composite (fiber friction + Darcy
# resistance of solvent percolating the ~ξ mesh); with the Biot fluid explicit, the bare fibers feel
# ONLY the solvent and the 65.9 EMERGES from the coupling — no η double-count (the retire-6πηR guard).
# Value: aqueous cytosol ≈ water; water at 37 °C = 6.9e-4 Pa·s → rounded to the roadmap-ratified 1e-3.
ETA_SOLVENT = 1.0e-3          # pN·s·µm⁻² (= 1e-3 Pa·s ≈ water @ 37 °C; FSI skeleton drag medium)


def fiber_point_drag_solvent(L_um: float, p_segments: int, *,
                             diameter_um: float = ACTIN_DIAMETER_UM,
                             l_hydro_um: float | None = None) -> float:
    """γ_solid — the SKELETON-only per-point drag = :func:`fiber_point_drag` at the SOLVENT viscosity.

    The FSI double-count resolution (roadmap 2026-07-16): once the spatial Biot fluid supplies the
    Darcy dissipation, the fiber points feel only the solvent (η≈1e-3 Pa·s), NOT the lumped cytoplasm
    η=65.9. This is the drag that REPLACES the whole-cell 6πηR time-scale when ``biot_fsi`` is ON;
    the effective 65.9 then emerges from the pore-pressure back-reaction rather than being imposed.
    """
    return fiber_point_drag(L_um, p_segments, eta=ETA_SOLVENT,
                            diameter_um=diameter_um, l_hydro_um=l_hydro_um)


# ── SI ↔ FF converters ───────────────────────────────────────────────────────────────────────────
def length_to_ff(x_m: np.ndarray | float) -> np.ndarray | float:
    """metres → FF length (µm)."""
    return np.asarray(x_m) / UM if not np.isscalar(x_m) else x_m / UM


def length_to_si(x_um: np.ndarray | float) -> np.ndarray | float:
    """FF length (µm) → metres."""
    return np.asarray(x_um) * UM if not np.isscalar(x_um) else x_um * UM


def force_to_ff(f_N):
    """newtons → FF force (pN)."""
    return np.asarray(f_N) / PN if not np.isscalar(f_N) else f_N / PN


def force_to_si(f_pN):
    """FF force (pN) → newtons."""
    return np.asarray(f_pN) * PN if not np.isscalar(f_pN) else f_pN * PN


def kappa_to_ff(kappa_Nm2):
    """bending modulus N·m² → FF (pN·µm²)."""
    return kappa_Nm2 / PN_UM2


def kappa_to_si(kappa_ff):
    """bending modulus FF (pN·µm²) → N·m²."""
    return kappa_ff * PN_UM2


def energy_to_ff(e_J):
    """energy J → FF (pN·µm)."""
    return e_J / PN_UM


# ── NF2007 §5.2 mobility / drag (in FF units) ─────────────────────────────────────────────────────
def fiber_mobility(L_um: float, *, eta: float = ETA_CYTOPLASM,
                   diameter_um: float = ACTIN_DIAMETER_UM, l_hydro_um: float | None = None) -> float:
    """Cytosim per-point mobility μ for a fiber (NF2007 §5.2, p10), in FF units (µm·pN⁻¹·s⁻¹).

        μ = log(L_h/δ) / (3π η L),   per-point μ_p = (p+1) μ

    Args:
        L_um: fiber length L [µm].
        eta: medium viscosity [pN·s·µm⁻²] (= Pa·s).
        diameter_um: filament diameter δ [µm].
        l_hydro_um: hydrodynamic cutoff L_h [µm]; the paper uses L_h = min(L, cutoff). Default
            ``L_h = L`` (no cutoff). Must satisfy L_h > δ.

    Returns:
        μ, the single-point mobility [µm·pN⁻¹·s⁻¹].
    """
    L_h = L_um if l_hydro_um is None else min(L_um, l_hydro_um)
    if L_h <= diameter_um:
        raise ValueError("hydrodynamic length L_h must exceed the filament diameter δ")
    return float(np.log(L_h / diameter_um) / (3.0 * np.pi * eta * L_um))


def fiber_point_drag(L_um: float, p_segments: int, *, eta: float = ETA_CYTOPLASM,
                     diameter_um: float = ACTIN_DIAMETER_UM, l_hydro_um: float | None = None) -> float:
    """Per-model-point overdamped drag γ = 1/μ_p for the implicit step (γẋ = F), in FF units.

    NF2007 derives a SINGLE mobility for the p+1 points of a fiber: μ_p = (p+1)·μ (§5.2). The
    overdamped law the implicit solver integrates is γ ẋ = F with γ the per-point drag, so
    γ = 1/μ_p [pN·s·µm⁻¹].
    """
    mu = fiber_mobility(L_um, eta=eta, diameter_um=diameter_um, l_hydro_um=l_hydro_um)
    return float(1.0 / ((p_segments + 1) * mu))
