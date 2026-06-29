"""B5 (i) — single-SF TENSION force-budget + Route-B NMII force-scaling probe.

PI directive (a)/(A), SESSION (i) "CLOSE THE FLOOR" (PLATFORM_PI_QUEUE.md). This
script quantifies the SF generation-limit (the SF instance of the cortical γ-floor)
and tests whether a REAL native-NMII density datum closes the Kumar 10-30 nN single-SF
tension band via the Route-B mesoscale force-scale factor — WITHOUT a frozen integrator
and WITHOUT back-solving the factor to hit the band.

Two walls (PLATFORM_PI_QUEUE / BASAL_MESH_DESIGN): (1) the physical EA backbone makes
explicit BAOAB dt tiny — SOLVED by the bending wire (bend-before-stretch → AFINES
soft-stretch); (2) the mesoscale myosin force budget is far below Kumar — the SAME
generation-limit as the cortical γ-floor. This script sidesteps (1) and quantifies (2).

═══════════════════════════════════════════════════════════════════════════════════
DEEP-RESEARCH RESULT (2026-06-09, 104-agent harness, 3-vote adversarial verify).
The literature decomposes the single-SF tension into pieces of VERY different quality:

  ANCHOR 2a — per-minifilament MOLECULAR content (SOLID, 3-0):
    Native NMII bipolar minifilament = ~28-30 molecules = ~56-60 heads, ~30 heads/SIDE,
    300 nm long (Billington 2013 JBC EM; Hu 2017 NCB 3D-SIM U2OS; Melli 2018 eLife EM;
    Niederman & Pollard 1975 platelet). The runtime models 10 heads/side @ 0.5 pN/head.

  ANCHOR 2b — per-minifilament FORCE (ORDER / EXTRAPOLATED, self-flagged):
    The contractile stall force of a single NMII minifilament is NOT directly measured.
    The ONLY primary estimate: fs ≈ 17 pN = 10 heads/side × 1.7 pN muscle-myosin head
    (Stachowiak & O'Shaughnessy 2009 Biophys J, explicitly "using ... muscle myosin II
    since nonmuscle myosin II forces have not been directly measured"). Per-head NMII
    force spans ~0.7-10 pN; unloaded duty ratio NM2-A ~0.05 (RISES under load — Kovacs
    2007), so the engaged-head fraction in a LOADED fiber is the open variable.

  ANCHOR 1 — parallel minifilaments per CROSS-SECTION (THE MISSING DATUM):
    *** No surviving primary source gives a direct measured count. *** The one explicit
    estimate (~50/cross-section, Stachowiak-derived from ~100 actin filaments ÷ 2) was
    REFUTED 0-3 in verification. Back-solving 10-30 nN ÷ 17 pN ⇒ ~590-1760 parallel
    minifilaments — GEOMETRICALLY IMPOSSIBLE: an SF cross-section is 50-250 nm radius,
    ~10-30 actin filaments across (MBC 2021) → can host O(5-15) minifilaments, not ~600.

  REFRAME — single-fiber ACTIVE vs NETWORK total (HIGH, 3-0 / 2-1):
    Kassianidou, Brand, Schwarz & Kumar 2017 (PNAS 114:2622, U2OS, the Kumar lab's OWN
    active-Kelvin-Voigt model): single-fiber aggregate motor stall force Fs = k·Lo ≈
    6 nN (k=3 nN/µm, Lo=2 µm); a CONNECTING SF adds only ~5 nN of ACTIVE myosin force,
    while the length-defined SF reaches ~25 nN at center — the rest is NETWORK/PRESTRESS.
    ⇒ Kumar 10-30 nN is NOT a pure single-fiber active-generation target; the active
    component is ~5-6 nN. Kumar 2006 itself reports STRESS (Pa), not a clean per-SF nN.

VERDICT LOGIC (stated BEFORE the run, no gate-loosening):
  * The DERIVABLE Route-B piece is the per-minifilament MOLECULAR correction (Anchor 2a/2b):
    bring the model minifilament (5 pN) up to the literature per-minifilament estimate
    (~17 pN). factor_mini = 17/5 ≈ 3.4×, grid-invariant, lit-anchored — NOT back-solved.
  * The parallel-cross-section count (Anchor 1) — the OTHER multiplier Route B needs — has
    NO usable datum (refuted; back-solve geometrically impossible). So the FULL factor that
    would reach Kumar CANNOT be derived from a density datum.
  * Even granting the molecular correction AND the (refuted-generous) ~50 parallel count:
    50 × 17 pN ≈ 0.85 nN — still ~7× under the ~6 nN single-fiber ACTIVE target and
    ~12-35× under Kumar 10-30 nN. With the geometric max (~5-15): ~0.09-0.26 nN.
  ⇒ REFUTE: the SF generation-limit is REAL and is NOT closable from a measured density
    datum. This is the SF instance of the cortical γ-floor (same ½·n·f·ℓ budget; same
    MISSING motor-density datum; same insight that the literature band includes passive/
    prestress the active motors alone do not supply). HALT→PI on the gate-reframe.

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

from ffn_sim.archive.hoomd_legacy.cell.basal_surface import build_flat_basal_surface
from ffn_sim.archive.hoomd_legacy.cell.basal_mesh import (
    build_basal_filament_network,
    place_sf_myosin_on_apparatus,
)
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin

_OUT = _HERE.parents[1] / "outputs" / "h7"

ELL0 = 0.5e-6
FOOT_R = 5.0e-6
Z_BASAL = -7.0e-6
N_FA = 60
N_CABLES = 12
N_INFILL = 500
N_MOTORS = 40
KUMAR_BAND_N = (10.0e-9, 30.0e-9)

# ── Literature anchors (deep-research 2026-06-09; provenance + confidence labelled) ──
# Anchor 2a — native minifilament molecular content (SOLID, 3-0).
NATIVE_HEADS_PER_SIDE = 30          # ~30 heads/side (≈30 molecules, ~60 heads total)
#   Billington 2013 JBC (PMC3829186) NM2A 29/58, NM2B 30/60; Hu 2017 NCB (28114270)
#   ~30 molecules, 300±20 nm; Melli 2018 eLife (32871) 30 (A/B); Niederman&Pollard
#   1975 JCB 67:72 28/56. Confidence SOLID.
# Anchor 2b — per-minifilament STALL force (ORDER / EXTRAPOLATED, self-flagged).
LIT_PER_MINIFILAMENT_STALL_N = 17.0e-12   # fs ≈ 17 pN (Stachowiak & O'Shaughnessy 2009
#   Biophys J PMC2711311 = 10 heads/side × 1.7 pN muscle-myosin; the ONLY primary
#   per-minifilament estimate, authors flag it is NOT a direct NMII measurement).
#   Confidence ORDER. (Canonical 30 heads/side × 1.7 pN ⇒ ~51 pN optimistic stall, but
#   unloaded duty 0.05 lowers it; 17 pN is the conservative literature anchor.)
# Anchor 1 — parallel minifilaments per CROSS-SECTION: *** NO USABLE DATUM ***.
PARALLEL_PER_CROSS_SECTION = None         # MISSING — the ~50 estimate was REFUTED 0-3;
#   back-solve (590-1760) is geometrically impossible. NONE-GATED (no magic number).
PARALLEL_REFUTED_GENEROUS = 50            # the REFUTED literature value (for the bound only)
PARALLEL_GEOMETRIC_MAX = (5, 15)          # SF cross-section 50-250 nm r, ~10-30 actin →
#   O(5-15) minifilaments (2 actin/minifilament). Geometric ceiling, not a count datum.
# Reframe — single-fiber ACTIVE myosin force (HIGH, 3-0): Kassianidou/Schwarz/Kumar 2017
SF_ACTIVE_TARGET_N = 6.0e-9               # ~5-6 nN single-fiber ACTIVE (PNAS 114:2622)
SF_NETWORK_TOTAL_N = 25.0e-9             # ~25 nN length-defined SF center (network+prestress)

# --- Literature reference minifilament (the SAME dipole the cortical γ-floor uses) ---
# One-sided NMII-IIA bipolar minifilament axial (dipole) force = heads-per-side ×
# per-head stall, the quantity in the cortical active-gel envelope γ=½·n2D·f·ℓ
# (H7_ACTIVE_GAMMA_SYNTHESIS_2026-06-09). This is the cross-line consistency anchor:
# the basal SF budget below uses the BRIEF-LITERAL minifilament (10 heads × 0.5 pN),
# which is ~11× SMALLER than this literature minifilament. That ~11× factors cleanly as
# (28/10 heads-per-side) × (2/0.5 pN per-head stall) = ~2.8× × ~4× — both are PARAMETER
# choices, NOT a delivery/stiffness issue: the quasi-static budget already uses F_stall
# directly (not the soft-k cap k·r), i.e. it already assumes the §9 continuous_stroke
# delivery fix. So the ~11× is the brief-vs-literature minifilament PARAMETERISATION, of
# the same magnitude/category as (but distinct mechanism from) the §9 per-head recovery.
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
    """Quasi-static single-SF tension force budget + Route-B molecular correction.

    Reports (a) the RAW mesoscale per-minifilament tension, (b) the DERIVABLE
    per-minifilament molecular correction toward the literature ~17 pN, and (c) the
    HONEST accounting of the MISSING parallel-cross-section count — never a back-solve.
    """
    phi = _engaged_fraction(p_sf)
    H = int(p_sf.n_heads_per_side)
    F_stall = float(p_sf.F_stall_per_head)

    # per-bipolar-minifilament axial contractile force (one side pulls inward;
    # ⟨|cos|⟩ ≈ 1 for minifilaments laid along their host filament tangent).
    cos_axial = 1.0
    f_mini = phi * H * F_stall * cos_axial      # model raw per-minifilament force [N]

    # raw single-mesoscale-cable tension: minifilaments along a cable are in series
    # (same tension); take f_mini as the per-cross-section motor force.
    T_cable_raw = f_mini

    lo, hi = KUMAR_BAND_N
    band_centre = 0.5 * (lo + hi)
    in_band_raw = bool(lo <= T_cable_raw <= hi)
    gap_to_lo = lo / T_cable_raw if T_cable_raw > 0 else float("inf")
    factor_needed = band_centre / T_cable_raw if T_cable_raw > 0 else float("inf")

    # ── Route-B DERIVABLE piece: per-minifilament MOLECULAR correction (Anchor 2a/2b) ──
    # Bring the model minifilament (10 heads/side × 0.5 pN) up to the literature
    # per-minifilament STALL estimate (~17 pN, Stachowiak 2009). This factor is
    # grid-invariant and lit-anchored (NOT chosen to pass the band).
    factor_mini = LIT_PER_MINIFILAMENT_STALL_N / f_mini if f_mini > 0 else float("inf")
    f_mini_corrected = LIT_PER_MINIFILAMENT_STALL_N      # = f_mini × factor_mini

    # ── The OTHER Route-B multiplier — parallel cross-section count — is MISSING ──
    # We do NOT pick a value. We report the single-fiber active tension under the
    # bounding parallel counts (refuted-generous 50, geometric 5-15) to show the
    # residual gap is irreducible without a measured count.
    gmin, gmax = PARALLEL_GEOMETRIC_MAX
    T_active_geom_lo = gmin * f_mini_corrected
    T_active_geom_hi = gmax * f_mini_corrected
    T_active_refuted = PARALLEL_REFUTED_GENEROUS * f_mini_corrected   # 50 × 17 pN

    # gap of the MOST GENEROUS (refuted-50) single-fiber active estimate vs targets.
    gap_active_vs_kassianidou = SF_ACTIVE_TARGET_N / T_active_refuted
    gap_active_vs_kumar_lo = lo / T_active_refuted

    return {
        # raw model minifilament
        "phi_engaged": phi,
        "n_heads_per_side_model": H,
        "F_stall_per_head_model_N": F_stall,
        "f_per_minifilament_model_N": f_mini,
        "T_cable_raw_N": T_cable_raw,
        "kumar_band_N": list(KUMAR_BAND_N),
        "in_kumar_band_raw": in_band_raw,
        "gap_to_band_lo": gap_to_lo,
        "mesoscale_force_factor_needed_to_band_centre": factor_needed,
        "n_motors": int(p_sf.n_motors_per_cell),
        # DERIVABLE molecular correction (Anchor 2a/2b)
        "native_heads_per_side": NATIVE_HEADS_PER_SIDE,
        "lit_per_minifilament_stall_N": LIT_PER_MINIFILAMENT_STALL_N,
        "factor_mini_molecular": factor_mini,
        "f_per_minifilament_corrected_N": f_mini_corrected,
        # MISSING parallel count (Anchor 1) — honest accounting, no back-solve
        "parallel_per_cross_section_datum": PARALLEL_PER_CROSS_SECTION,  # None = MISSING
        "parallel_refuted_generous": PARALLEL_REFUTED_GENEROUS,
        "parallel_geometric_max": list(PARALLEL_GEOMETRIC_MAX),
        "parallel_backsolved_to_kumar": [lo / LIT_PER_MINIFILAMENT_STALL_N,
                                         hi / LIT_PER_MINIFILAMENT_STALL_N],
        # single-fiber active tension under bounding parallel counts
        "T_active_geometric_N": [T_active_geom_lo, T_active_geom_hi],
        "T_active_refuted50_N": T_active_refuted,
        "sf_active_target_N": SF_ACTIVE_TARGET_N,
        "sf_network_total_N": SF_NETWORK_TOTAL_N,
        "gap_active_refuted50_vs_kassianidou": gap_active_vs_kassianidou,
        "gap_active_refuted50_vs_kumar_lo": gap_active_vs_kumar_lo,
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
        cortical γ-floor uses. ~11× = (28/10 heads) × (2/0.5 pN per-head) = ~2.8× × ~4×,
        both PARAMETER choices (NOT delivery: the budget uses F_stall directly, already
        assuming the §9 fix). Same magnitude/category as the cortical §9 per-head
        recovery, distinct mechanism — a cross-line parameter-fidelity item.
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
    # accept either key convention: the production force_budget() emits the *_model
    # suffixed keys (session i); the decomposition fixture/contract (session ii) uses
    # the bare names. Bridge both so the merged consumer works with either producer.
    n_heads = budget.get("n_heads_per_side_model", budget.get("n_heads_per_side"))
    f_stall = budget.get("F_stall_per_head_model_N", budget.get("F_stall_per_head_N"))
    f_brief_full = n_heads * f_stall
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
    # REFUTE: the raw tension is far under band AND the full Route-B factor is not
    # derivable (parallel count missing/refuted/geometrically impossible). The
    # decomposition factors that same gap into labeled physical sources (session ii).
    decomposition = decompose_generation_gap(budget)
    verdict = "REFUTE" if not budget["in_kumar_band_raw"] else "PASS"
    return {
        "gate": "B5(i) SF tension force-budget + Route-B NMII force-scaling probe",
        "verdict": verdict,
        "physics_claim": True,
        "provisional": False,   # the MISSING-datum finding is decisive, not provisional
        "halt_to_pi": True,
        "budget": budget,
        "decomposition": decomposition,
        "note": (
            "REFUTE — the SF generation-limit is NOT closable from a measured NMII "
            "density datum. The DERIVABLE Route-B piece is the per-minifilament "
            f"MOLECULAR correction (model {budget['f_per_minifilament_model_N']:.2e} N → "
            f"literature {budget['lit_per_minifilament_stall_N']:.2e} N, "
            f"factor_mini≈{budget['factor_mini_molecular']:.1f}×, Stachowiak 2009 ORDER; "
            f"native {budget['native_heads_per_side']} heads/side SOLID). The OTHER "
            "Route-B multiplier — parallel minifilaments per cross-section — has NO "
            "usable datum: the ~50 literature estimate was REFUTED 0-3, and the "
            "back-solve (~590-1760) is geometrically impossible (cross-section hosts "
            f"O(5-15)). Even the refuted-generous 50 × {budget['lit_per_minifilament_stall_N']:.0e} N "
            f"= {budget['T_active_refuted50_N']:.2e} N is "
            f"~{budget['gap_active_refuted50_vs_kassianidou']:.0f}× under the single-fiber "
            f"ACTIVE target (~6 nN, Kassianidou/Schwarz/Kumar 2017 PNAS) and "
            f"~{budget['gap_active_refuted50_vs_kumar_lo']:.0f}× under Kumar 10 nN. REFRAME: "
            "Kumar 10-30 nN conflates single-fiber active myosin (~5-6 nN) with NETWORK/"
            "PRESTRESS (~25 nN); the active component is itself unclosed (load-dependent "
            "duty ratio is the open variable). SF instance of the cortical γ-floor — "
            "same ½·n·f·ℓ budget, same missing motor-density datum. HALT→PI on the gate-"
            "reframe (active-only target ~6 nN vs network-total 10-30 nN)."
        ),
    }


def main() -> int:
    rep = run()
    (_OUT / "production").mkdir(parents=True, exist_ok=True)
    path = _OUT / "production" / "h7_basal_sf_force_budget.json"
    with open(path, "w") as fh:
        json.dump(rep, fh, indent=2, default=str)
    b = rep["budget"]
    print(f"[B5(i) SF force-budget + Route-B probe] {rep['verdict']} (HALT→PI)")
    print(f"  model: φ={b['phi_engaged']:.3f}, {b['n_heads_per_side_model']} heads/side, "
          f"F_stall={b['F_stall_per_head_model_N']:.1e} N → "
          f"f/minifilament={b['f_per_minifilament_model_N']:.2e} N")
    print(f"  Kumar band {b['kumar_band_N'][0]:.0e}-{b['kumar_band_N'][1]:.0e} N → "
          f"raw {b['gap_to_band_lo']:.0f}× under lo (in_band={b['in_kumar_band_raw']})")
    print(f"  mesoscale force-scale factor to band centre: "
          f"~{b['mesoscale_force_factor_needed_to_band_centre']:.0f}x (Route B)")
    # ── session (i): molecular correction + missing parallel-count accounting ──
    print(f"  DERIVABLE molecular correction (Anchor 2a/2b, Stachowiak 2009 ORDER): "
          f"native {b['native_heads_per_side']} heads/side, lit "
          f"{b['lit_per_minifilament_stall_N']:.1e} N/minifilament → "
          f"factor_mini≈{b['factor_mini_molecular']:.1f}×")
    print(f"  MISSING datum (Anchor 1): parallel/cross-section = {b['parallel_per_cross_section_datum']} "
          f"(refuted ~50; geometric {b['parallel_geometric_max']}; back-solve "
          f"{b['parallel_backsolved_to_kumar'][0]:.0f}-{b['parallel_backsolved_to_kumar'][1]:.0f} impossible)")
    print(f"  refuted-generous 50 × {b['lit_per_minifilament_stall_N']:.0e} N = "
          f"{b['T_active_refuted50_N']:.2e} N → "
          f"~{b['gap_active_refuted50_vs_kassianidou']:.0f}× under ~6 nN ACTIVE target, "
          f"~{b['gap_active_refuted50_vs_kumar_lo']:.0f}× under Kumar lo")
    print(f"  REFRAME: single-fiber ACTIVE ~{b['sf_active_target_N']:.0e} N (Kassianidou 2017) "
          f"vs NETWORK total ~{b['sf_network_total_N']:.0e} N → Kumar band is not pure active")
    # ── session (ii): labeled gap decomposition (cross-line with cortical γ-floor) ──
    d = rep["decomposition"]
    print("  --- gap decomposition (labeled, N-independent) ---")
    print(f"    per-minifilament fidelity (brief 10×0.5pN → lit 28×2pN={d['lit_minifilament_N']*1e12:.0f}pN): "
          f"~{d['per_minifilament_fidelity_factor']:.0f}x "
          f"(=2.8× heads × 4× per-head stall; params, budget already assumes §9 delivery)")
    print(f"    cross-sectional NMII count to band centre: "
          f"~{d['cross_sectional_nmii_count_to_band_centre']:.0f} lit-minifilaments/SF "
          f"(= Route-B MCF7 density datum, session (i))")
    print(f"    engagement realism penalty (φ0.99 → realised "
          f"{d['cortex_realised_engagement']*100:.0f}%): ~{d['engagement_realism_penalty_factor']:.0f}x WORSE")
    print(f"  → SF tension is GENERATION-limited (same ½·n·f·ℓ as cortical γ-floor), "
          f"N_filaments-INDEPENDENT. HALT→PI.")
    print(f"  ⚠ UPPER bound on a COHERENT bundle; §11 SF-2c: random polarity → "
          f"−61 pN slackening → needs SARCOMERIC organisation.")
    print(f"  json: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
