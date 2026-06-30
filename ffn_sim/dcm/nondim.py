"""Shared nondimensionalization for the DCM (and, by design, the FF) engine.

PI 2026-06-29 ("FF와 공유할 무차원화"): give both engines ONE dimensionless language so a
cortical surface tension γ — whichever engine produces it — maps onto the same physics. The
**FF engine measures γ** from its Cytosim filament network (the γ-floor work, Stage-6d/6e on
`dcm/main`); the **DCM engine consumes γ** as the cortical-tension input to faceting. This module
is the bridge: it converts a γ [N/m] into the DCM faceting group γ̃ = γ/(K·ℓ) and reports which
SimuCell3D regime it lands in — so FF's γ slots directly onto the DCM γ-sweep
(``dcm.gamma_sweep``) without re-deriving scales on each side.

Nothing here is tuned to an outcome — it is a scale registry + dimensionless-group calculators
anchored to measured values (provenance inline). γ stays a controlled/measured INPUT, never a
knob fit to make a morphology appear.

Promotion: this is prototyped in ``dcm/`` (the DCM session's owned domain). Once FF reviews the
shared scale set it should move to ``ffn_sim/common/`` (engine-agnostic) and be imported by both.

Run:
    PYTHONPATH=. python -m ffn_sim.dcm.nondim                 # print the scale table + γ ladder
    PYTHONPATH=. python -m ffn_sim.dcm.nondim --gamma 1e-3    # locate an (FF-measured) γ
    PYTHONPATH=. python -m ffn_sim.dcm.nondim --gamma 1e-3 --k-vol 2.5e3
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np

#: SimuCell3D faceting band in γ̃ = γ/(K·ℓ) (Runser, Vetter & Iber 2024, Nat. Comput. Sci.).
FACETING_BAND = (0.02, 0.10)


@dataclass(frozen=True)
class CellScales:
    """Measured physical scales of the MCF7 DCM cell (SI). Provenance inline — no tuned values.

    These are the *reference* scales; a specific run may differ (e.g. K is a controlled variable
    in the γ-sweep, ℓ is recomputed from the actual mesh). Defaults mirror the engine
    (``geometry.ResolvedDCM`` + ``dcm_warp_decohesion``).
    """
    R_cell: float = 7.5e-6        # m   MCF7 radius (Wagner 2011) — geometry.ResolvedDCM.R_cell
    K_vol: float = 7.73e5         # Pa  osmotic/bulk modulus (driver default; single-cell-spread
    #                                   tuning). SimuCell3D/Fischer-Friedrich faceting K ≈ 2.5e3.
    K_simucell3d: float = 2.5e3   # Pa  SimuCell3D/Fischer-Friedrich 2014 cytoplasm bulk modulus.
    turgor_dP0: float = 133.0     # Pa  MCF7 baseline osmotic turgor (geometry.ResolvedDCM)
    w_cs: float = 2.85e-3         # J/m² MCF7 cell–substrate adhesion energy density
    eta_cytoplasm: float = 65.9   # Pa·s MCF7 cytoplasm viscosity (Dessard 2024)

    @property
    def V0(self) -> float:
        """Rest cell volume [m³] of the reference sphere."""
        return (4.0 / 3.0) * np.pi * self.R_cell ** 3

    @property
    def ell(self) -> float:
        """Cell length scale ℓ = V₀^(1/3) [m] (the faceting-group length; ≈ 12.1 µm)."""
        return self.V0 ** (1.0 / 3.0)


# ---------------------------------------------------------------------------
# dimensionless groups (the shared language)
# ---------------------------------------------------------------------------
def gamma_tilde(gamma: float, K: float, ell: float) -> float:
    """DCM faceting group γ̃ = γ/(K·ℓ). In [0.02, 0.10] ⇒ SimuCell3D foam-like faceting."""
    denom = K * ell
    return float(gamma / denom) if denom > 0 else 0.0


def gamma_from_tilde(gt: float, K: float, ell: float) -> float:
    """Inverse: the γ [N/m] that yields a target γ̃ at (K, ℓ) — for designing a band-crossing sweep."""
    return float(gt * K * ell)


def faceting_regime(gt: float) -> str:
    """Classify a γ̃ into the SimuCell3D morphology regime (controlled read-out, not a gate)."""
    lo, hi = FACETING_BAND
    if gt < lo:
        return "rounded (turgor-dominated marbles)"
    if gt <= hi:
        return "FACETING band (foam-like polygonal junctions)"
    return "cortex-crushed (cortex overwhelms turgor)"


def elastocapillary_length(gamma: float, K: float) -> float:
    """λ_ec = γ/K [m]: below this length, cortical tension dominates bulk elasticity."""
    return float(gamma / K) if K > 0 else float("inf")


def turgor_capillary_ratio(dP0: float, gamma: float, R: float) -> float:
    """Young–Laplace balance dP₀·R/(2γ): >1 turgor wins (rounds/inflates), <1 cortex wins."""
    denom = 2.0 * gamma
    return float(dP0 * R / denom) if denom > 0 else float("inf")


def douezan_spreading(w_cs: float, gamma: float) -> float:
    """Douezan dimensionless spreading s = w_cs/(2γ) − 1 (S = w_cs − 2γ). >0 ⇒ wetting/spread."""
    return float(w_cs / (2.0 * gamma) - 1.0) if gamma > 0 else float("inf")


def viscous_capillary_time(eta: float, ell: float, gamma: float) -> float:
    """Capillary-viscous relaxation time τ = η·ℓ/γ [s]: cortical-tension shape-relaxation scale."""
    return float(eta * ell / gamma) if gamma > 0 else float("inf")


def describe_gamma(gamma: float, scales: CellScales, K: float | None = None) -> dict:
    """Full dimensionless characterization of a cortical tension γ [N/m] (e.g. FF-measured).

    Returns the faceting group + regime + the companion balances, at bulk modulus ``K``
    (defaults to the driver K_vol; pass ``scales.K_simucell3d`` for the faceting-physics K).
    """
    K = scales.K_vol if K is None else K
    ell = scales.ell
    gt = gamma_tilde(gamma, K, ell)
    return {
        "gamma_Npm": gamma,
        "K_Pa": K,
        "ell_um": ell * 1e6,
        "gamma_tilde": gt,
        "regime": faceting_regime(gt),
        "in_faceting_band": FACETING_BAND[0] <= gt <= FACETING_BAND[1],
        "elastocapillary_length_um": elastocapillary_length(gamma, K) * 1e6,
        "turgor_capillary_ratio": turgor_capillary_ratio(scales.turgor_dP0, gamma, scales.R_cell),
        "douezan_s": douezan_spreading(scales.w_cs, gamma),
        "viscous_capillary_time_s": viscous_capillary_time(scales.eta_cytoplasm, ell, gamma),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gamma", type=float, default=None,
                    help="a cortical tension γ [N/m] (e.g. FF-measured) to locate in the regime")
    ap.add_argument("--k-vol", type=float, default=None,
                    help="bulk modulus K [Pa] (default: driver K_vol 7.73e5; "
                         "use 2.5e3 for SimuCell3D faceting physics)")
    args = ap.parse_args()
    s = CellScales()
    K = args.k_vol if args.k_vol is not None else s.K_vol

    print("=== DCM/FF shared cell scales (SI, measured) ===")
    print(f"  R_cell        = {s.R_cell*1e6:.2f} µm   (Wagner 2011)")
    print(f"  ℓ = V₀^(1/3)  = {s.ell*1e6:.2f} µm")
    print(f"  K_vol (driver)= {s.K_vol:.3e} Pa   (single-cell-spread tuning)")
    print(f"  K_SimuCell3D  = {s.K_simucell3d:.3e} Pa   (Fischer-Friedrich 2014; faceting physics)")
    print(f"  turgor dP₀    = {s.turgor_dP0:.1f} Pa")
    print(f"  w_cs          = {s.w_cs:.3e} J/m²   (MCF7)")
    print(f"  η_cytoplasm   = {s.eta_cytoplasm:.1f} Pa·s   (Dessard 2024)")
    print(f"  faceting band γ̃ ∈ {FACETING_BAND}  (SimuCell3D)")

    if args.gamma is not None:
        print(f"\n=== locate γ = {args.gamma:.3e} N/m  (K = {K:.3e} Pa) ===")
        for k, v in describe_gamma(args.gamma, s, K=K).items():
            print(f"  {k:28s} = {v}")
    else:
        print(f"\n=== γ ladder → γ̃ regime  (K = {K:.3e} Pa, ℓ = {s.ell*1e6:.2f} µm) ===")
        print("  to cross the faceting band at this K, γ must span "
              f"[{gamma_from_tilde(FACETING_BAND[0], K, s.ell):.2e}, "
              f"{gamma_from_tilde(FACETING_BAND[1], K, s.ell):.2e}] N/m")
        for g in (1e-4, 5e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2):
            gt = gamma_tilde(g, K, s.ell)
            print(f"  γ={g:.1e} N/m  γ̃={gt:.4f}  {faceting_regime(gt)}")


if __name__ == "__main__":
    main()
