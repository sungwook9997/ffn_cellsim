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
    verdict = "REFUTE" if not budget["in_kumar_band_raw"] else "PASS"
    return {
        "gate": "B5(i) SF tension force-budget (quasi-static, generation-limit probe)",
        "verdict": verdict,
        "physics_claim": True,
        "provisional": True,
        "budget": budget,
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
    print(f"  → SF tension is GENERATION-limited (same as cortical γ-floor), "
          f"N_filaments-INDEPENDENT.")
    print(f"  json: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
