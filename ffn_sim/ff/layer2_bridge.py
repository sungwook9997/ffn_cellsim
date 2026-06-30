"""FF → Layer-2 scale bridge (Stage 6i / Decision 3a) — FF single-cell γ feeds the spheroid A/A0 law.

The Layer-2 spheroid line reproduces the PI experimental curvature law A/A0 = a + b/R + c/R² from a
single-cell cortical tension γ via Young-Laplace (the b/R term) + cell–cell adhesion β (the c/R²
term). The existing bridge (`scripts/layer2_surface_tension_bridge.py`) fed a FIXED
`resolved.cortical_tension` ≈ 0.57 mN/m — which this session reclassified as the g_rigid turgor
double-book. This module supplies γ from the FF FIRST-PRINCIPLES single cell instead:

  γ_cortical = γ_passive (state-dependent osmotic turgor, ΔP·R/2; Guo 2017, 40 Pa) + γ_active
               (actomyosin, the floored method-of-planes value)

so the spheroid FORM is now grounded in the FF cortex, and the MAGNITUDE inherits the FF single-cell
γ — including its force-bearing-density floor (FF_STAGE6D/6H). Key honest point: the A/A0 FORM
(a + b/R + c/R²) follows from the Young-Laplace + adhesion STRUCTURE regardless of the γ magnitude
(the Layer-2 line already showed r²≈0.998); FF sets the COEFFICIENTS (b ∝ γ_cortical, c ∝ β), which
carry the single-cell floor up to the spheroid — magnitude is NOT a free knob.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ffn_sim.ff.cortex_assembly import CortexParams
from ffn_sim.ff.gamma_floor import (
    TURGOR_DP0,
    build_crosslinked_cortex,
    equilibrate,
    gamma_passive_young_laplace,
    measure_gamma,
)


@dataclass(frozen=True, slots=True)
class FFCellTension:
    """FF single-cell cortical-tension components [mN/m] feeding the Layer-2 bridge."""

    gamma_active_mN_m: float     # actomyosin (method-of-planes) — the floored channel
    gamma_passive_mN_m: float    # state-dependent osmotic turgor Young-Laplace ΔP·R/2
    gamma_total_mN_m: float      # active + passive (what enters the tissue Young-Laplace)
    R_um: float


def ff_single_cell_tension(*, f_myo: float = 5.0, n_filaments: int = 100, n_xl: int = 300,
                           n_myo: int = 100, seed: int = 0, n_steps: int = 300) -> FFCellTension:
    """Compute the FF first-principles single-cell cortical tension (active + passive turgor) [mN/m].

    Builds + settles the FF crosslinked cortex, applies the lit myosin prestress, reads the
    method-of-planes actomyosin γ and the state-dependent osmotic turgor passive γ. 1 pN/µm = 1e-3 mN/m.
    """
    rng = np.random.default_rng(seed)
    cortex = build_crosslinked_cortex(CortexParams(), n_filaments=n_filaments, n_xl=n_xl,
                                      n_myo=n_myo, rng=rng)
    equilibrate(cortex, 0.0, n_steps=n_steps, turgor=False)        # resting shell
    m = measure_gamma(cortex, f_myo, turgor=True)
    ga = m["gamma_active"] * 1e-3
    gp = m["gamma_passive"] * 1e-3
    return FFCellTension(gamma_active_mN_m=ga, gamma_passive_mN_m=gp,
                         gamma_total_mN_m=ga + gp, R_um=cortex.R_um)


def aa0_law(R_um: np.ndarray, gamma_cortical_mN_m: float, beta_mN_m: float,
            a: float = 1.0, k_b: float = 1.0, k_c: float = 1.0) -> np.ndarray:
    """Spheroid A/A0 curvature law A/A0 = a + b/R + c/R² with b ∝ γ_cortical, c ∝ adhesion β.

    The Young-Laplace 1/R term (b = k_b·γ_cortical·µm) flattens larger spheroids less; the cell–cell
    interfacial-tension term (c = k_c·β·µm²) is the higher-order correction. ``k_b``/``k_c`` are the
    geometric proportionality constants from the Layer-2 surface-tension bridge (NOT tuned here —
    placeholders 1.0; the real values come from `scripts/layer2_surface_tension_bridge`). FORM is
    robust; FF sets the coefficient magnitudes.
    """
    R = np.asarray(R_um, dtype=np.float64)
    b = k_b * gamma_cortical_mN_m
    c = k_c * beta_mN_m
    return a + b / R + c / R**2


def bridge_summary(*, beta_frac: float = 0.4, **kw) -> dict:
    """FF γ → A/A0 coefficients over a spheroid radius range. ``beta_frac`` = β/γ (cell–cell vs cortex)."""
    t = ff_single_cell_tension(**kw)
    beta = beta_frac * t.gamma_total_mN_m
    R = np.linspace(20.0, 120.0, 11)        # spheroid radii [µm] (PI experimental range)
    aa0 = aa0_law(R, t.gamma_total_mN_m, beta)
    return {"tension": t, "beta_mN_m": beta, "R_um": R.tolist(), "aa0": aa0.tolist(),
            "b_coeff": t.gamma_total_mN_m, "c_coeff": beta}
