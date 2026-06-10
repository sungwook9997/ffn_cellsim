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
    # derivable (parallel count missing/refuted/geometrically impossible).
    verdict = "REFUTE"
    return {
        "gate": "B5(i) SF tension force-budget + Route-B NMII force-scaling probe",
        "verdict": verdict,
        "physics_claim": True,
        "provisional": False,   # the MISSING-datum finding is decisive, not provisional
        "halt_to_pi": True,
        "budget": budget,
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
          f"raw {b['gap_to_band_lo']:.0f}× under lo")
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
    print(f"  → SF generation-limit REAL, not density-closable. HALT→PI. json: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
