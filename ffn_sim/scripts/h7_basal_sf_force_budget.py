"""B5 (i) — quasi-static SF TENSION force-budget (KU-3.5, PI directive (a)/(A)).

SKETCH / DIAGNOSTIC. The full dynamic B5 hits two known walls (PLATFORM_PI_QUEUE /
BASAL_MESH_DESIGN): (1) the physical EA backbone k≈1.75 N/m makes explicit BAOAB
dt≈2e-11 s → ~1e10 steps for a 0.5 s contraction (needs the FROZEN native constrained
M-SHAKE integrator), and (2) the mesoscale myosin force budget is far below the Kumar
10-30 nN band — the SAME generation-limit as the cortical γ-floor. This script
sidesteps the integrator wall and quantifies wall (2) directly: it computes the
STEADY-STATE single-SF tension as a motor FORCE BUDGET and compares it to Kumar.

Key physics (why this is N-independent):
  At force balance a contractile bundle's tension = the aggregate MOTOR force pulling
  across a cross-section, NOT μ_SF·ε with μ_SF set by N_filaments. μ_SF only sets the
  STRAIN ε = T/μ_SF at that tension; the tension itself is motor-generation-limited.
  So sweeping N changes ε, not the tension ceiling — the Kumar comparison is a
  GENERATION question, exactly like the cortical γ-floor.

Budget model (assumptions stated explicitly; PROVISIONAL):
  * engaged fraction φ = k_on/(k_on + k_off0)  (zero-load steady; an UPPER estimate).
  * per-bipolar-minifilament axial contractile force
        f_mini = φ · n_heads_per_side · F_stall_per_head · ⟨|cos|⟩_axial.
  * along ONE mesoscale cable, minifilaments are in SERIES → same tension →
        T_cable(mesoscale) ≈ f_mini · (mean # minifilaments engaged on the cable, as a
        coherence/duty proxy ≥ 1).
  * the mesoscale cable stands for ~N_filaments native filaments in PARALLEL across the
    bundle cross-section, plus the native:effective minifilament ratio → the
    PHYSICAL bundle tension needs a mesoscale FORCE-SCALE factor (Route B). We report
    the RAW mesoscale tension AND the factor needed to reach the Kumar band.

This is option (i): no frozen integrator, decisive on whether SF hits the γ-floor.
Run:  python ffn_sim/scripts/h7_basal_sf_force_budget.py
"""

from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

_HERE = Path(__file__).resolve()
if str(_HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(_HERE.parents[2]))

import numpy as np
import yaml

from ffn_sim.cell.basal_surface import build_flat_basal_surface
from ffn_sim.cell.basal_mesh import (
    build_basal_filament_network,
    place_sf_myosin_on_apparatus,
)
from ffn_sim.cortex.myosin import resolve_cortex_myosin

_OUT = _HERE.parents[1] / "outputs" / "h7"

ELL0 = 0.5e-6
FOOT_R = 5.0e-6
Z_BASAL = -7.0e-6
N_FA = 60
N_CABLES = 12
N_INFILL = 500
N_MOTORS = 40
KUMAR_BAND_N = (10.0e-9, 30.0e-9)

# --- Literature reference minifilament (the SAME dipole the cortical γ-floor uses) ---
# One-sided NMII-IIA bipolar minifilament axial (dipole) force = heads-per-side ×
# per-head stall, the quantity in the cortical active-gel envelope γ=½·n2D·f·ℓ
# (H7_ACTIVE_GAMMA_SYNTHESIS_2026-06-09). This is the cross-line consistency anchor:
# the basal SF budget below uses the BRIEF-LITERAL minifilament (10 heads × 0.5 pN),
# which is ~11× SMALLER than this literature minifilament — exactly the per-head /
# stiffness fidelity the cortical line corrected in the §9 continuous_stroke redesign.
LIT_HEADS_PER_SIDE = 28           # Billington 2013 (NMII-IIA ~28-30 heads/side)
LIT_F_PER_HEAD_N = 2.0e-12        # 2 pN/head one-sided dipole (Chugh 2017 × Billington 2013)
LIT_F_MINIFILAMENT_N = LIT_HEADS_PER_SIDE * LIT_F_PER_HEAD_N   # ≈ 56 pN
# Steady-state REALISED engagement on the fine-grained cortex (geometric/availability
# limited, NOT the kinetic φ≈0.99): h7_engagement_saturation.json, 11.4% plateau.
CORTEX_REALISED_ENGAGEMENT = 0.114


def _engaged_fraction(p_sf) -> float:
    """Zero-load steady bound fraction φ = k_on/(k_on+k_off0) (UPPER estimate)."""
    return float(p_sf.head_actin_k_on / (p_sf.head_actin_k_on + p_sf.head_actin_k_off0))


def force_budget(app, p_sf, myo) -> dict:
    """Quasi-static single-SF (cable) tension force budget vs Kumar."""
    phi = _engaged_fraction(p_sf)
    H = int(p_sf.n_heads_per_side)
    F_stall = float(p_sf.F_stall_per_head)

    # per-bipolar-minifilament axial contractile force (one side pulls inward;
    # ⟨|cos|⟩ ≈ 1 for minifilaments laid along their host filament tangent).
    cos_axial = 1.0
    f_mini = phi * H * F_stall * cos_axial

    # raw single-mesoscale-cable tension: minifilaments along a cable are in series
    # (same tension); take f_mini as the per-cross-section motor force.
    T_cable_raw = f_mini

    lo, hi = KUMAR_BAND_N
    in_band_raw = bool(lo <= T_cable_raw <= hi)
    gap_to_lo = lo / T_cable_raw if T_cable_raw > 0 else float("inf")

    # mesoscale FORCE-SCALE factor needed to reach the Kumar band centre.
    band_centre = 0.5 * (lo + hi)
    factor_needed = band_centre / T_cable_raw if T_cable_raw > 0 else float("inf")

    return {
        "phi_engaged": phi,
        "n_heads_per_side": H,
        "F_stall_per_head_N": F_stall,
        "f_per_minifilament_N": f_mini,
        "T_cable_raw_N": T_cable_raw,
        "kumar_band_N": list(KUMAR_BAND_N),
        "in_kumar_band_raw": in_band_raw,
        "gap_to_band_lo": gap_to_lo,
        "mesoscale_force_factor_needed_to_band_centre": factor_needed,
        "n_motors": int(p_sf.n_motors_per_cell),
    }


def decompose_generation_gap(budget: dict) -> dict:
    """Factor the raw SF generation gap into LABELED physical sources (mandate ①).

    The raw budget (``force_budget``) is N_filaments-INDEPENDENT and ~2014× under the
    Kumar single-SF band. This decomposes that gap so each multiplier is attributable
    to a named physical/parameter cause rather than one opaque "force-scale factor":

      gap = (per-minifilament fidelity) × (cross-sectional NMII count) ,  evaluated at
            the UPPER engagement φ≈0.99; a realistic engagement makes it WORSE.

    Factors (all assumptions stated, none tuned to pass):
      * ``per_minifilament_fidelity`` — the brief-literal minifilament (10 heads ×
        0.5 pN ≈ 5 pN) vs the LITERATURE minifilament (28 heads × 2 pN ≈ 56 pN) the
        cortical γ-floor uses. ~11×. This is the SAME per-head/stiffness correction the
        cortical line executed in the §9 continuous_stroke redesign — i.e. ~11× of the
        SF gap is a cross-line parameter-fidelity item already being closed elsewhere.
      * ``cross_sectional_nmii_count`` — after the per-minifilament fix, the residual
        gap = how many LITERATURE minifilaments must act coherently across ONE SF
        cross-section to reach the Kumar band centre (Kumar/f_lit). This is the
        Route-B native:effective force-scale factor = the MISSING MCF7 SF NMII density
        datum (the deep-research target of the sibling "close-the-floor" session).
      * ``engagement_realism_penalty`` — the budget uses the zero-load UPPER φ≈0.99;
        the fine-grained cortex realises only ~11% engagement (geometric/availability
        limited). A realistic SF would be ~``0.99/0.114`` ≈ 9× WORSE than this budget,
        not better — reported so the bound is not read optimistically.

    Returns the factors + the implied per-SF minifilament count, with the explicit
    note that ALL of this is an UPPER bound on a COHERENT bundle: the §11 SF-2c result
    (random mixed-polarity → −61 pN slackening, not contraction) shows the budget is
    only realised with SARCOMERIC POLARITY ORGANISATION (Hotulainen-Lappalainen).
    """
    band_lo, band_hi = budget["kumar_band_N"]
    band_centre = 0.5 * (band_lo + band_hi)

    # per-minifilament fidelity compares the FULL (un-engaged) brief minifilament
    # (H·F_stall) against the FULL literature dipole — engagement φ is factored ONCE,
    # in the cross-sectional count below, so fidelity × count == the raw factor exactly.
    f_brief_full = budget["n_heads_per_side"] * budget["F_stall_per_head_N"]
    per_mini_fidelity = LIT_F_MINIFILAMENT_N / f_brief_full if f_brief_full > 0 else float("inf")
    # residual gap once the minifilament is the literature dipole (engagement still φ≈0.99)
    f_lit_engaged = budget["phi_engaged"] * LIT_F_MINIFILAMENT_N
    cross_section_count_to_centre = band_centre / f_lit_engaged if f_lit_engaged > 0 else float("inf")
    cross_section_count_to_lo = band_lo / f_lit_engaged if f_lit_engaged > 0 else float("inf")
    engagement_penalty = budget["phi_engaged"] / CORTEX_REALISED_ENGAGEMENT

    return {
        "per_minifilament_fidelity_factor": per_mini_fidelity,
        "lit_minifilament_N": LIT_F_MINIFILAMENT_N,
        "lit_heads_per_side": LIT_HEADS_PER_SIDE,
        "lit_f_per_head_N": LIT_F_PER_HEAD_N,
        "cross_sectional_nmii_count_to_band_centre": cross_section_count_to_centre,
        "cross_sectional_nmii_count_to_band_lo": cross_section_count_to_lo,
        "engagement_realism_penalty_factor": engagement_penalty,
        "cortex_realised_engagement": CORTEX_REALISED_ENGAGEMENT,
        "decomposition_check_product": per_mini_fidelity * cross_section_count_to_centre,
        "unified_with_cortical_gamma_floor": (
            "SAME ½·n·f·ℓ generation budget. Cortex: γ_active=½·n2D·f·ℓ, ~36-80× under "
            "the MCF7 active target → needs n2D≈16-47/µm² vs the HeLa-proxy 0.6/µm² "
            "(NO MCF7 datum). SF: T=N_cross·f, ~180× under Kumar after the per-minifilament "
            "fix → needs N_cross≈"
            f"{cross_section_count_to_centre:.0f} literature minifilaments across one SF "
            "cross-section (NO MCF7 SF-NMII datum). Same f≈56 pN dipole, same density-datum "
            "gap (areal for cortex, cross-sectional for SF), same architecture caveat: count "
            "is the WRONG lever — cortex tension is overlap-set (Chugh 2017 / Truong Quang "
            "2021), SF tension is sarcomeric-polarity-set (§11 SF-2c: random polarity → "
            "−61 pN slackening; Hotulainen-Lappalainen 2006)."
        ),
        "upper_bound_caveat": (
            "All factors evaluated at the UPPER engagement φ≈0.99 AND assume coherent "
            "axial summation (cos=1). Both are optimistic: realised engagement ~11% (×9 "
            "worse) and a random mixed-polarity bundle does NOT rectify (§11). The budget "
            "is the CEILING a perfectly-organised sarcomeric SF could reach, not a prediction "
            "of the current random-polarity construction."
        ),
    }


def run() -> dict:
    surf = build_flat_basal_surface(
        footprint_radius=FOOT_R, z_basal=Z_BASAL, n_rings=10, n_fa=N_FA,
        rng=np.random.default_rng(7),
    )
    app = build_basal_filament_network(
        surf, ell0=ELL0, n_cables=N_CABLES, n_infill=N_INFILL,
        rng=np.random.default_rng(8),
    )
    cfg = deepcopy(yaml.safe_load(open(_HERE.parents[1] / "configs" / "phase1_h3.yaml")))
    cfg["cortex"]["myosin"]["prefix"] = "sf_myosin_"
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = N_MOTORS
    p_sf = resolve_cortex_myosin(cfg, dt=1e-9)
    myo = place_sf_myosin_on_apparatus(app, p_sf, rng=np.random.default_rng(11))

    budget = force_budget(app, p_sf, myo)
    decomposition = decompose_generation_gap(budget)
    verdict = "REFUTE" if not budget["in_kumar_band_raw"] else "PASS"
    return {
        "gate": "B5(i) SF tension force-budget (quasi-static, generation-limit probe)",
        "verdict": verdict,
        "physics_claim": True,
        "provisional": True,
        "budget": budget,
        "decomposition": decomposition,
        "note": (
            "QUASI-STATIC FORCE BUDGET (no stiff-backbone dynamics, no frozen "
            "integrator). The single-SF tension is MOTOR-generation-limited and "
            "N_filaments-INDEPENDENT (N sets strain ε=T/μ_SF, not the tension). "
            f"Raw mesoscale tension {budget['T_cable_raw_N']:.2e} N vs Kumar "
            f"{KUMAR_BAND_N[0]:.0e}-{KUMAR_BAND_N[1]:.0e} N → "
            f"{budget['gap_to_band_lo']:.0f}x under band: the SAME generation-limit as "
            "the cortical γ-floor. Reaching Kumar needs either the mesoscale "
            f"force-scale factor (~{budget['mesoscale_force_factor_needed_to_band_centre']:.0f}x, "
            "Route B native:effective parallel-bundle scaling) or a native/density "
            "fix — NOT a backbone-stiffness (N_filaments) change. ASSUMPTIONS: φ "
            "zero-load upper estimate; per-minifilament series tension; provisional "
            "pending the cross-simulator deep-research + B5-method decision."
        ),
    }


def main() -> int:
    rep = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_basal_sf_force_budget.json"
    with open(path, "w") as fh:
        json.dump(rep, fh, indent=2, default=str)
    b = rep["budget"]
    print(f"[B5(i) SF force-budget] {rep['verdict']} (provisional)")
    print(f"  φ_engaged={b['phi_engaged']:.3f}, {b['n_heads_per_side']} heads/side, "
          f"F_stall={b['F_stall_per_head_N']:.1e} N")
    print(f"  f/minifilament={b['f_per_minifilament_N']:.2e} N → "
          f"T_cable(raw mesoscale)={b['T_cable_raw_N']:.2e} N")
    print(f"  Kumar band {b['kumar_band_N'][0]:.0e}-{b['kumar_band_N'][1]:.0e} N → "
          f"in_band={b['in_kumar_band_raw']}, {b['gap_to_band_lo']:.0f}x UNDER lo")
    print(f"  mesoscale force-scale factor to band centre: "
          f"~{b['mesoscale_force_factor_needed_to_band_centre']:.0f}x (Route B)")
    d = rep["decomposition"]
    print("  --- gap decomposition (labeled, N-independent) ---")
    print(f"    per-minifilament fidelity (brief 5pN → lit {d['lit_minifilament_N']*1e12:.0f}pN): "
          f"~{d['per_minifilament_fidelity_factor']:.0f}x "
          f"(= the §9 cortical per-head fix)")
    print(f"    cross-sectional NMII count to band centre: "
          f"~{d['cross_sectional_nmii_count_to_band_centre']:.0f} lit-minifilaments/SF "
          f"(= Route-B MCF7 density datum, session (i))")
    print(f"    engagement realism penalty (φ0.99 → realised "
          f"{d['cortex_realised_engagement']*100:.0f}%): ~{d['engagement_realism_penalty_factor']:.0f}x WORSE")
    print(f"  → SF tension is GENERATION-limited (same ½·n·f·ℓ as cortical γ-floor), "
          f"N_filaments-INDEPENDENT.")
    print(f"  ⚠ UPPER bound on a COHERENT bundle; §11 SF-2c: random polarity → "
          f"−61 pN slackening → needs SARCOMERIC organisation.")
    print(f"  json: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
