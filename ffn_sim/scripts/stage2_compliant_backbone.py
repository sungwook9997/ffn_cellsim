"""STAGE-2 / Decision-4(a) — does a COMPLIANT (harmonic) backbone let the cortex
build ACTIVE shell tension, or does it only RELABEL the structural tension channel?

Context (verified 2026-06-04). The active cortical-tension floor lives in ``g_soft``
(harmonic-bond method-of-planes), while ``g_rigid`` (the M-SHAKE rigid-backbone Lagrange
tension) is already IN the KU-3.5 band [0.35,0.65] mN/m. In every constrained run
``r/r0 = 1.00000`` (the rigid shell never contracts). PI Decision-4(a): allow the cortex
to actually contract by using a COMPLIANT backbone.

KEY PHYSICS CAVEAT (the reason this script exists). A round cell is a PRESSURIZED sphere:
its cortical tension is the tangential hoop stress balanced by the enclosed-volume pressure
(Young-Laplace ΔP = 2γ/R) — it exists AT r/r0 = 1.0, it does NOT require radial shrink. With
a RIGID backbone that hoop stress is carried by the constraint Lagrange multipliers →
measured in ``g_rigid`` (in-band). If we simply make the backbone COMPLIANT (harmonic,
``cell.build(constrained=False)`` → backbone k = bond_k), that SAME structural hoop stress
moves into the harmonic bonds → measured in ``g_soft``. That would make ``g_soft`` jump into
band WITHOUT any new active contraction — a measurement-channel RELABEL, not a fix.

So Decision-4(a) is only a real fix if the compliant backbone enables NEW active contraction
(e.g. Murrell-Gardel buckling-mediated contractility, suppressed by the rigid constraint).
The decisive test is the motors ON vs OFF CONTROL on the SAME compliant backbone:

    g_soft(OFF) = structural hoop stress relabelled into the harmonic channel (the baseline)
    g_soft(ON)  = structural + any NEW active contraction
    Δ_active = g_soft(ON) − g_soft(OFF)   ← the genuine active cortical tension

  - Δ_active reaches band  → compliant backbone enables real active contractility (4a WORKS)
  - g_soft(OFF) jumps, Δ_active floored → pure RELABEL artifact (4a is NOT a fix; report so)

This reuses the existing full-cell builder (NO frozen-integrator change: ``constrained=False``
is an already-supported mode — it is what the warm-up runs). NOT a contract change.

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.stage2_compliant_backbone --smoke
      python -m ffn_sim.scripts.stage2_compliant_backbone --n-fil 1000 --n-samples 6
Outputs: ffn_sim/outputs/h3/production/stage2_compliant_backbone.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import hoomd

from ffn_sim.common.production_policy import (
    add_production_device_args,
    require_full_cell_physiological_baseline,
    validate_production_device_args,
)
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.scripts.h3_ku35_tension import _tension_method_of_planes
from ffn_sim.scripts.mcf7_fullcell_stage1 import (
    _build,
    _cfg_for,
    _resolve_compartments,
    _tagpos,
)

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "h3" / "production"
BAND = (0.35, 0.65)   # KU-3.5 mN/m


def _run_compliant(n_fil, n_motors, n_xl, *, n_warmup, n_sample, interval,
                   bind_scale, kon_scale, device="gpu", seed=1):
    """Build + run a COMPLIANT-backbone (constrained=False) full cell; sample g_soft.

    Returns (samples list of (r/r0, g_soft_mN_m), exploded: bool).
    """
    cfg = _cfg_for(n_fil, n_motors, n_xl, "grip_walk", force_scaling=True, backbone_nm=300)
    p0 = resolve_h3_derived(cfg)
    comp = _resolve_compartments(p0)
    require_full_cell_physiological_baseline(comp)
    # constrained=False → cortex backbone bond k = p_cortex.bond_k (harmonic = COMPLIANT);
    # equilibrate softstarts the stiff harmonic stretch mode at the small unconstrained dt.
    p, _p_myo, _p_xl, dtc, hw = _build(
        cfg, stepping_mode="grip_walk", force_scaling=True, constrained=False,
        compartments=comp, equilibrate=True, n_warmup=n_warmup, device=device,
        kon_scale=kon_scale, bind_scale=bind_scale, seed=seed,
    )
    sim = hw["sim"]
    sim.run(0)
    nca = p.n_filaments * p.beads_per_filament
    pos0 = _tagpos(sim)
    r0 = float(np.linalg.norm(pos0[:nca], axis=1).mean())
    samples = []
    exploded = False
    for k in range(n_sample):
        sim.run(interval)
        pos = _tagpos(sim)
        if not np.all(np.isfinite(pos)):
            exploded = True
            break
        rmean = float(np.linalg.norm(pos[:nca], axis=1).mean())
        g_soft = _tension_method_of_planes(sim, p.R_cell)
        samples.append((rmean / r0, g_soft * 1e3))   # mN/m
        print(f"    [compliant n_motors={n_motors}] s={k+1}/{n_sample} "
              f"r/r0={rmean/r0:.5f} g_soft={g_soft*1e3:.4e} mN/m", flush=True)
    return samples, exploded


def main() -> int:
    ap = argparse.ArgumentParser(description="Compliant-backbone motors ON/OFF control (Decision-4a).")
    ap.add_argument("--smoke", action="store_true", help="tiny: n_fil=300, 3 samples")
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--n-motors", type=int, default=100)
    ap.add_argument("--n-samples", type=int, default=6)
    ap.add_argument("--interval", type=int, default=8000)
    ap.add_argument("--n-warmup", type=int, default=8000)
    ap.add_argument("--bind-scale", type=float, default=6.0, help="mesoscale-consistent reach")
    ap.add_argument("--kon-scale", type=float, default=300.0, help="couple_accel → percolate")
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args, hoomd_module=hoomd)

    if args.smoke:
        n_fil, n_motors, n_sample, interval, n_warmup = 300, 60, 3, 3000, 3000
    else:
        n_fil, n_motors, n_sample, interval, n_warmup = (
            args.n_fil, args.n_motors, args.n_samples, args.interval, args.n_warmup)
    n_xl = int(1.5 * n_fil)   # n_xl/n_fil=1.5 → z≈3 (Ennomani optimum) with bind_scale=6

    print(f"[compliant-backbone] n_fil={n_fil} n_xl={n_xl} bind_scale={args.bind_scale} "
          f"kon_scale={args.kon_scale} | samples={n_sample}×{interval} after warmup={n_warmup}")
    print(f"[control] measuring g_soft with motors ON (n={n_motors}) vs OFF (n=0) on the "
          f"SAME compliant backbone; Δ_active = ON − OFF is the genuine active tension.\n")

    print("  --- ARM A: motors OFF (structural-relabel baseline) ---")
    off, off_exp = _run_compliant(n_fil, 0, n_xl, n_warmup=n_warmup, n_sample=n_sample,
                                  interval=interval, bind_scale=args.bind_scale,
                                  kon_scale=args.kon_scale, device=args.device)
    print("  --- ARM B: motors ON (structural + active) ---")
    on, on_exp = _run_compliant(n_fil, n_motors, n_xl, n_warmup=n_warmup, n_sample=n_sample,
                                interval=interval, bind_scale=args.bind_scale,
                                kon_scale=args.kon_scale, device=args.device)

    def _mean_gsoft(s):
        return float(np.mean([g for _, g in s])) if s else float("nan")

    def _mean_rr0(s):
        return float(np.mean([r for r, _ in s])) if s else float("nan")

    gs_off, gs_on = _mean_gsoft(off), _mean_gsoft(on)
    delta = gs_on - gs_off
    verdict = (
        "EXPLODED (compliant stretch unstable at this dt — the reason production uses M-SHAKE)"
        if (off_exp or on_exp) else
        "REAL active contractility (Δ_active in band)" if BAND[0] <= delta <= BAND[1] else
        "PARTIAL (Δ_active lifted but under band)" if delta > 0.05 else
        "RELABEL artifact (Δ_active floored; any g_soft rise is structural, not active)"
    )

    res = {
        "config": {"n_fil": n_fil, "n_motors_on": n_motors, "n_xl": n_xl,
                   "bind_scale": args.bind_scale, "kon_scale": args.kon_scale,
                   "n_sample": n_sample, "interval": interval, "n_warmup": n_warmup,
                   "backbone": "compliant (constrained=False, harmonic k=bond_k)"},
        "motors_off": {"g_soft_mean_mN_m": gs_off, "r_over_r0_mean": _mean_rr0(off),
                       "samples": off, "exploded": off_exp},
        "motors_on": {"g_soft_mean_mN_m": gs_on, "r_over_r0_mean": _mean_rr0(on),
                      "samples": on, "exploded": on_exp},
        "delta_active_mN_m": delta,
        "ku35_band_mN_m": list(BAND),
        "verdict": verdict,
        "interpretation": (
            "g_soft(OFF)=structural hoop stress relabelled into the harmonic channel; "
            "g_soft(ON)=structural+active; Δ_active=ON−OFF is the genuine active cortical "
            "tension. Compliant backbone is a real fix ONLY if Δ_active reaches band."
        ),
    }
    print(f"\n=== RESULT (compliant backbone, motors ON/OFF control) ===")
    print(f"  g_soft(OFF, structural) = {gs_off:.4e} mN/m   (r/r0={_mean_rr0(off):.5f})")
    print(f"  g_soft(ON,  total)      = {gs_on:.4e} mN/m   (r/r0={_mean_rr0(on):.5f})")
    print(f"  Δ_active = ON − OFF     = {delta:.4e} mN/m   (band {BAND})")
    print(f"  → VERDICT: {verdict}")

    _OUT.mkdir(parents=True, exist_ok=True)
    out = _OUT / "stage2_compliant_backbone.json"
    out.write_text(json.dumps(res, indent=2))
    print(f"\n[json] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
