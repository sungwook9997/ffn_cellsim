"""STAGE-2 force-budget bottleneck diagnostic — does signed contractility SCALE with motors?

The SIGNED-contractility measurement (stage2_signed_contractility.py) showed the cortex is
already NET CONTRACTILE (coherent +, no dipole cancellation) but ~100x (vs Warmt single-cell)
to ~2700x (vs the slow-rounded-cell band) too WEAK. So the wall is the FORCE BUDGET, not motor
arrangement. This sweep pins WHICH budget lever is the bottleneck by sweeping the one directly
tunable knob — motor count n_motors — and reading, per condition:

  * signed_gamma  — net contractile hoop stress (the stage2_signed_contractility observable)
  * n_bound / bound_frac — how many myosin heads are actually engaged (engaged-head budget)
  * mean F/F_stall on bound heads — realized per-head load (is it stall-limited or slack?)

Diagnosis logic:
  - signed_gamma ∝ n_motors (LINEAR, unsaturating) AND bound_frac ~const ⇒ pure ENGAGED-COUNT
    lever: more motors ⇒ proportionally more contraction; extrapolate n_motors → band/Warmt.
    The fix is recruitment/density (Bohec 2026 sustained recruitment), a clean magnitude axis.
  - signed_gamma SATURATES while n_motors rises ⇒ a downstream CEILING (transmission, realized
    load, or measurement) caps it — adding motors won't help; the lever is elsewhere.
  - bound_frac falls as n_motors rises ⇒ engagement (binding) is the limiter, not raw count.
  - F/F_stall << 1 and flat ⇒ heads run slack (sustained-load/transport lever, not count).

NOT a fix and NOT integrator-touching — a read-only-physics diagnostic that tells us exactly
which 'core part' to change. Validated by the SIGNED observable (contraction, not |magnitude|).

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.stage2_forcebudget_sweep --n-fil 300
Outputs: ffn_sim/outputs/h3/production/stage2_forcebudget_sweep.json
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
from ffn_sim.scripts.mcf7_fullcell_stage1 import _build, _cfg_for, _resolve_compartments, _tagpos
from ffn_sim.scripts.stage2_signed_contractility import signed_cortical_stress

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "h3" / "production"


def _bound_metrics(ma, sim):
    """(bound_frac, mean_F_over_stall_on_bound, n_bound) from the myosin action."""
    p = ma.p
    bound = ma._head_bound_to_actin >= 0
    n_bound = int(bound.sum())
    n_heads = int(ma._head_bound_to_actin.size)
    if n_bound == 0:
        return 0.0, 0.0, 0, n_heads
    pos = _tagpos(sim)
    bidx = np.flatnonzero(bound)
    htags = np.array([ma._head_global_tag(int(h)) for h in bidx], dtype=np.int64)
    atags = ma._head_bound_to_actin[bidx]
    r = np.linalg.norm(pos[htags] - pos[atags], axis=1)
    F = p.k_head_actin * np.clip(r, 0.0, None)
    return float(n_bound / n_heads), float(np.mean(F / p.F_stall_per_head)), n_bound, n_heads


def _run_condition(n_fil, n_motors, n_xl, *, n_warmup, n_sample, interval,
                   bind_scale, kon_scale, device="gpu", seed=1):
    cfg = _cfg_for(n_fil, n_motors, n_xl, "grip_walk", force_scaling=True, backbone_nm=300)
    p0 = resolve_h3_derived(cfg)
    comp = _resolve_compartments(p0)
    require_full_cell_physiological_baseline(comp)
    _, _, _, dtc, hw = _build(cfg, stepping_mode="grip_walk", force_scaling=True,
                              constrained=False, compartments=comp, equilibrate=True,
                              n_warmup=n_warmup, device=device, kon_scale=kon_scale,
                              bind_scale=bind_scale, seed=seed)
    hw["sim"].run(0)
    pos_warm = _tagpos(hw["sim"])
    del hw
    p, _pm, _px, dtc, hc = _build(cfg, stepping_mode="grip_walk", force_scaling=True,
                                  constrained=True, compartments=comp, dtc=dtc, device=device,
                                  kon_scale=kon_scale, bind_scale=bind_scale, seed=seed)
    sim = hc["sim"]
    ma = hc.get("myosin_action")
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm
    sim.state.set_snapshot(snap)
    sim.run(0)
    signed, boundf, fos, nb = [], [], [], []
    for _ in range(n_sample):
        sim.run(interval)
        st = signed_cortical_stress(sim, p.R_cell)
        signed.append(st["signed_mN_m"])
        if ma is not None and n_motors > 0:
            bf, f_over, n_bound, n_heads = _bound_metrics(ma, sim)
            boundf.append(bf); fos.append(f_over); nb.append(n_bound)
    return {
        "n_motors": n_motors,
        "signed_gamma_mN_m": float(np.mean(signed)) if signed else float("nan"),
        "bound_frac": float(np.mean(boundf)) if boundf else 0.0,
        "mean_F_over_stall": float(np.mean(fos)) if fos else 0.0,
        "n_bound": float(np.mean(nb)) if nb else 0.0,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Force-budget bottleneck sweep (signed contractility vs motor count).")
    ap.add_argument("--n-fil", type=int, default=300)
    ap.add_argument("--motors", type=int, nargs="+", default=[0, 50, 100, 200, 400])
    ap.add_argument("--n-samples", type=int, default=4)
    ap.add_argument("--interval", type=int, default=6000)
    ap.add_argument("--n-warmup", type=int, default=6000)
    ap.add_argument("--bind-scale", type=float, default=6.0)
    ap.add_argument("--kon-scale", type=float, default=300.0)
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args, hoomd_module=hoomd)
    n_xl = int(1.5 * args.n_fil)

    print(f"[force-budget sweep] n_fil={args.n_fil} n_xl={n_xl} motors={args.motors} | "
          f"signed γ + engaged-head + F/F_stall vs motor count\n")
    conditions = []
    for nm in args.motors:
        r = _run_condition(args.n_fil, nm, n_xl, n_warmup=args.n_warmup, n_sample=args.n_samples,
                           interval=args.interval, bind_scale=args.bind_scale,
                           kon_scale=args.kon_scale, device=args.device)
        conditions.append(r)
        print(f"  n_motors={nm:4d}: signed γ={r['signed_gamma_mN_m']:+.4e} mN/m  "
              f"bound_frac={r['bound_frac']:.3f}  n_bound={r['n_bound']:.0f}  "
              f"F/F_stall={r['mean_F_over_stall']:.3f}", flush=True)

    # scaling analysis (exclude n_motors=0 baseline)
    act = [c for c in conditions if c["n_motors"] > 0 and np.isfinite(c["signed_gamma_mN_m"])]
    verdict = "inconclusive"
    if len(act) >= 2:
        nm = np.array([c["n_motors"] for c in act], float)
        g = np.array([c["signed_gamma_mN_m"] for c in act], float)
        # per-motor contractility (slope sense): is γ/n_motors ~const (linear) or falling (saturating)?
        per_motor = g / nm
        ratio = per_motor[-1] / per_motor[0] if per_motor[0] != 0 else float("nan")
        bf = np.array([c["bound_frac"] for c in act], float)
        if ratio > 0.7:
            verdict = (f"LINEAR in motor count (per-motor γ ~const, ratio={ratio:.2f}) → "
                       f"ENGAGED-COUNT/recruitment is the lever; extrapolate n_motors→target")
        elif ratio < 0.4:
            verdict = (f"SATURATING (per-motor γ falls {ratio:.2f}×) → a downstream CEILING caps "
                       f"contraction; more motors won't help (bound_frac trend {bf[0]:.2f}→{bf[-1]:.2f})")
        else:
            verdict = f"SUB-LINEAR (per-motor γ ratio={ratio:.2f}) → partial ceiling; mixed lever"

    res = {
        "config": {"n_fil": args.n_fil, "n_xl": n_xl, "motors": args.motors,
                   "n_sample": args.n_samples, "interval": args.interval},
        "conditions": conditions,
        "verdict": verdict,
    }
    print(f"\n=== FORCE-BUDGET VERDICT ===\n  {verdict}")
    _OUT.mkdir(parents=True, exist_ok=True)
    out = _OUT / "stage2_forcebudget_sweep.json"
    out.write_text(json.dumps(res, indent=2))
    print(f"\n[json] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
