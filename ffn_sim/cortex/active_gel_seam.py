"""H.7 active-gel γ-seam — M1 0-D DIAGNOSTIC postprocessor (NOT a solution).

PI 2026-06-07: implement the active-gel seam M1 as a *diagnostic*, not a fix. It
takes the fine-grained (FG) cell's measured active-stress drive and relaxes a 0-D
active-Maxwell gel over the seconds timescale the MD cannot reach, then asks the
decisive question the instantaneous MD γ cannot:

  Is the active cortical-tension floor a TIMESCALE GAP (the myosin *could* make
  band-level tension, the µs-dt MD just couldn't relax to it) or a GENERATION
  LIMIT (the myosin force-dipole density is itself too small — no amount of
  relaxation reaches the band)?

The diagnostic computes the steady tension three ways and compares to the band:
  * γ_ss_realized — relax the FG-measured active stress (γ_soft/h). What the cell
    currently generates.
  * γ_ss_capacity — the engaged heads (N_eng) pulling at the per-head stall force.
  * γ_ss_ceiling  — ALL heads (N_total) at stall: the absolute generation ceiling.

If even the CEILING is below the band, the floor is a generation limit and the
seam (which only adds time) cannot close it → go to M2 (1-D shell/flow) or
re-target the gate observable (PI contingency). This module draws no KU-3.5
conclusion; the acceptance oracle is the analytic active-Maxwell solution.

Active-gel constitutive law (0-D, ε̇=0; Jülicher 2018): τ·σ̇ + σ = ζΔμ, so
σ(t) = ζΔμ + (σ_prestress − ζΔμ)·e^(−t/τ), steady σ_ss = ζΔμ. γ = σ·h.
The active-stress drive ζΔμ ≡ σ_a is the FG export (force-dipole density), NOT a
free fit; τ is the turnover relaxation time.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Cortical-tension reference band [N/m] (literature overlay, not a fit).
BAND_N_PER_M = (0.35e-3, 0.65e-3)
# MCF7-specific interphase datum (Hosseini 2020, Adv Sci, AFM) [N/m].
MCF7_INTERPHASE_GAMMA = 0.27e-3


def active_stress_from_dipoles(
    n_heads: int, f_per_head: float, ell_dipole: float, area: float, thickness: float
) -> float:
    """Actomyosin active stress σ_a [Pa] = Σ(f·ℓ) force-dipole density.

    σ_a = n_heads · f_per_head · ell_dipole / (area · thickness). The
    corresponding cortical tension is γ = σ_a · thickness = n·f·ℓ/area [N/m].
    """
    return n_heads * f_per_head * ell_dipole / (area * thickness)


def relax_active_maxwell(
    sigma_a: float, sigma_prestress: float, tau: float, *, t_max: float | None = None,
    n_steps: int = 400,
) -> tuple[np.ndarray, np.ndarray, float]:
    """0-D active-Maxwell relaxation τσ̇ + σ = σ_a (ε̇=0).

    Returns (t, sigma(t), sigma_steady). Analytic: σ(t) = σ_a + (σ_pre−σ_a)e^(−t/τ).
    The transient illustrates the seconds-scale buildup the µs-dt MD cannot reach;
    the steady value is σ_a (the active drive), independent of the prestress.
    """
    if t_max is None:
        t_max = 6.0 * tau
    t = np.linspace(0.0, t_max, n_steps)
    sigma = sigma_a + (sigma_prestress - sigma_a) * np.exp(-t / tau)
    return t, sigma, float(sigma_a)


@dataclass(slots=True)
class SeamDiagnosis:
    gamma_soft: float        # FG instantaneous active tension [N/m]
    gamma_rigid: float       # FG structural (M-SHAKE) tension [N/m]
    n_eng: int               # engaged myosin heads (bound attach bonds)
    n_total_heads: int
    f_stall_per_head: float  # [N]
    ell_dipole: float        # [m]
    area: float              # cortex area [m²]
    thickness: float         # cortex thickness [m]
    tau: float               # turnover relaxation time [s]
    motor_density_per_um2: float
    gamma_ss_realized: float  # [N/m]
    gamma_ss_capacity: float  # [N/m]
    gamma_ss_ceiling: float   # [N/m]
    verdict: str

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.__slots__}


def diagnose_seam(
    *, gamma_soft: float, gamma_rigid: float, n_eng: int, n_total_heads: int,
    f_stall_per_head: float, ell_dipole: float, area: float, thickness: float,
    tau: float,
) -> SeamDiagnosis:
    """Run the M1 diagnostic: realized / capacity / ceiling steady tension + verdict."""
    # γ = σ_a · h = n·f·ℓ/area for the capacity/ceiling; realized = the FG soft tension.
    g_ceiling = n_total_heads * f_stall_per_head * ell_dipole / area
    g_capacity = n_eng * f_stall_per_head * ell_dipole / area
    g_realized = gamma_soft  # σ_a_realized·h = (γ_soft/h)·h
    lo, hi = BAND_N_PER_M
    motor_density = (n_total_heads / 20.0) / (area * 1e12)  # heads→minifilaments (20/mf) per µm²

    if g_ceiling < lo:
        verdict = (
            f"GENERATION-LIMITED: even the ceiling (all {n_total_heads} heads at stall) "
            f"= {g_ceiling*1e3:.4f} mN/m is BELOW band [{lo*1e3:.2f},{hi*1e3:.2f}] "
            f"({lo/g_ceiling:.0f}x short). The seam (time only) CANNOT close this — "
            "fix the myosin generation (density/force) or re-target the observable (M2)."
        )
    elif g_capacity < lo <= g_ceiling:
        verdict = (
            "ENGAGEMENT-LIMITED: the ceiling reaches band but only "
            f"{n_eng}/{n_total_heads} heads engage → realized stays low. Fix myosin binding."
        )
    elif g_realized < lo <= g_capacity:
        verdict = (
            "TIMESCALE-GAP: the engaged heads' capacity reaches band but the MD-instant "
            "γ_soft is low → the seam's seconds-scale relaxation should close it."
        )
    else:
        verdict = "γ reaches band in the FG measure already (no seam needed)."

    return SeamDiagnosis(
        gamma_soft=gamma_soft, gamma_rigid=gamma_rigid, n_eng=int(n_eng),
        n_total_heads=int(n_total_heads), f_stall_per_head=f_stall_per_head,
        ell_dipole=ell_dipole, area=area, thickness=thickness, tau=tau,
        motor_density_per_um2=motor_density,
        gamma_ss_realized=g_realized, gamma_ss_capacity=g_capacity,
        gamma_ss_ceiling=g_ceiling, verdict=verdict,
    )
