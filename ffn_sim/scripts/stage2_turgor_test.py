"""STAGE-2 turgor test — does baseline osmotic pressure PRE-TENSION the cortex to band?

PI-ratified 2026-06-04: a rounded cell is an osmotically-PRESSURISED drop, so the cortex
should carry a Young-Laplace hoop tension γ ≈ Π₀·R/2 from the resting intracellular pressure
Π₀ ALONE (Stewart 2011: ~40 Pa interphase → ~400 Pa metaphase), before any myosin. The prior
enclosed-volume model had Π₀ = 0 (force-free at construction) — the missing pressure that the
saturation diagnostic implicated as the ceiling. This adds the literature turgor (opt-in,
``turgor_dP0``) and tests, MOTORS OFF, whether cortical tension tracks Π₀·R/2.

Prediction (R_cell = 7.5 µm MCF7): γ = Π₀·R/2 ⇒
    Π₀ =  40 Pa → γ ≈ 0.15 mN/m   (interphase; just under the band [0.35,0.65])
    Π₀ = 200 Pa → γ ≈ 0.75 mN/m   (above band)
    Π₀ = 400 Pa → γ ≈ 1.50 mN/m   (metaphase, mitotic rounding)
    band [0.35,0.65] mN/m ⇔ Π₀ ∈ [93, 173] Pa.

Measures, motors OFF, vs turgor Π₀: g_soft (harmonic method-of-planes, gated), g_rigid
(rigid-backbone Lagrange tension), g_total, signed γ, and r/r0. If g_total ≈ Π₀·R/2, the
pressure-coupling hypothesis is CONFIRMED: cortical tension is dominantly turgor-set (the
structural channel), reachable to band by the literature osmotic pressure — and the myosin
'active floor' was the wrong target (myosin only modulates on top, Bohec 2026).

Run:  conda activate ffn_sim
      python -m ffn_sim.scripts.stage2_turgor_test --smoke
      python -m ffn_sim.scripts.stage2_turgor_test --n-fil 300 --turgor 0 40 200 400
Outputs: ffn_sim/outputs/h3/production/stage2_turgor_test.json
"""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.scripts.h3_ku35_tension import (
    _tension_method_of_planes,
    _tension_method_of_planes_rigid,
)
from ffn_sim.scripts.mcf7_fullcell_stage1 import (
    _build,
    _cfg_for,
    _resolve_compartments,
    _tagpos,
)
from ffn_sim.scripts.stage2_signed_contractility import signed_cortical_stress

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "h3" / "production"
BAND = (0.35, 0.65)


def _run_turgor(n_fil, n_motors, n_xl, turgor_dP0, *, n_warmup, n_sample, interval,
                bind_scale, kon_scale, device="cpu", seed=1):
    cfg = _cfg_for(n_fil, n_motors, n_xl, "grip_walk", force_scaling=True, backbone_nm=300)
    p0 = resolve_h3_derived(cfg)
    comp = _resolve_compartments(p0)
    R = comp["p_enclosed_volume"].R_cell
    # inject the baseline turgor into the enclosed-volume compartment (opt-in)
    comp["p_enclosed_volume"] = replace(comp["p_enclosed_volume"], turgor_dP0=float(turgor_dP0))
    yl_pred = turgor_dP0 * R / 2.0 * 1e3  # Young-Laplace γ = Π₀·R/2 [mN/m]

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
    act = hc["baoab_action"]
    if hasattr(act, "record_lambda"):
        act.record_lambda = True
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm
    sim.state.set_snapshot(snap)
    sim.run(0)
    nca = p.n_filaments * p.beads_per_filament
    r0 = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
    gs, gr, sg, rr = [], [], [], []
    for _ in range(n_sample):
        sim.run(interval)
        gs.append(_tension_method_of_planes(sim, p.R_cell) * 1e3)
        gr.append(_tension_method_of_planes_rigid(sim, act, p.R_cell, dtc) * 1e3)
        sg.append(signed_cortical_stress(sim, p.R_cell)["signed_mN_m"])
        rmean = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
        rr.append(rmean / r0)
    g_soft = float(np.mean(gs)); g_rigid = float(np.mean(gr))
    return {
        "turgor_dP0_Pa": float(turgor_dP0),
        "yl_pred_mN_m": yl_pred,
        "g_soft_mN_m": g_soft,
        "g_rigid_mN_m": g_rigid,
        "g_total_mN_m": g_soft + g_rigid,
        "signed_soft_mN_m": float(np.mean(sg)),
        "r_over_r0": float(np.mean(rr)),
        "n_motors": n_motors,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Turgor pre-tension test (cortical tension vs Π₀).")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--n-fil", type=int, default=300)
    ap.add_argument("--turgor", type=float, nargs="+", default=[0, 40, 200, 400])
    ap.add_argument("--n-motors", type=int, default=0, help="0 = isolate pure pressure effect")
    ap.add_argument("--n-samples", type=int, default=4)
    ap.add_argument("--interval", type=int, default=6000)
    ap.add_argument("--n-warmup", type=int, default=6000)
    ap.add_argument("--bind-scale", type=float, default=6.0)
    ap.add_argument("--kon-scale", type=float, default=300.0)
    ap.add_argument("--device", default="cpu", choices=["cpu", "gpu"])
    args = ap.parse_args()
    if args.smoke:
        n_fil, turgor, n_sample, interval, n_warmup = 200, [0, 40, 400], 2, 2500, 2500
    else:
        n_fil, turgor, n_sample, interval, n_warmup = (
            args.n_fil, args.turgor, args.n_samples, args.interval, args.n_warmup)
    n_xl = int(1.5 * n_fil)

    print(f"[turgor test] n_fil={n_fil} n_motors={args.n_motors} turgor(Pa)={turgor} | "
          f"cortical tension vs baseline osmotic pressure (Young-Laplace γ=Π₀·R/2)\n")
    conditions = []
    for dp in turgor:
        r = _run_turgor(n_fil, args.n_motors, n_xl, dp, n_warmup=n_warmup, n_sample=n_sample,
                        interval=interval, bind_scale=args.bind_scale, kon_scale=args.kon_scale,
                        device=args.device)
        conditions.append(r)
        inb = BAND[0] <= r["g_total_mN_m"] <= BAND[1]
        print(f"  Π₀={dp:6.0f} Pa: g_total={r['g_total_mN_m']:.4e} "
              f"(soft={r['g_soft_mN_m']:.3e} rigid={r['g_rigid_mN_m']:.3e}) "
              f"YL_pred={r['yl_pred_mN_m']:.4e} mN/m  r/r0={r['r_over_r0']:.5f} "
              f"{'[IN BAND]' if inb else ''}", flush=True)

    # does g_total track the Young-Laplace prediction Π₀·R/2?
    act = [c for c in conditions if c["turgor_dP0_Pa"] > 0]
    verdict = "inconclusive"
    if act:
        ratios = [c["g_total_mN_m"] / c["yl_pred_mN_m"] for c in act if c["yl_pred_mN_m"] > 0]
        mean_ratio = float(np.mean(ratios)) if ratios else float("nan")
        any_band = any(BAND[0] <= c["g_total_mN_m"] <= BAND[1] for c in conditions)
        if 0.5 <= mean_ratio <= 2.0:
            verdict = (f"CONFIRMED — cortical tension TRACKS Young-Laplace Π₀·R/2 "
                       f"(g_total/YL ≈ {mean_ratio:.2f}); turgor pre-tensions the cortex, "
                       f"band reached at Π₀≈{2*BAND[0]/ (act[0]['g_total_mN_m']/act[0]['turgor_dP0_Pa']) if act[0]['g_total_mN_m']>0 else float('nan'):.0f} Pa. "
                       f"The myosin 'active floor' was the wrong target.")
        elif mean_ratio < 0.5:
            verdict = (f"PARTIAL — turgor lifts tension but BELOW Young-Laplace "
                       f"(g_total/YL ≈ {mean_ratio:.2f}); the shell partly inflates/relaxes "
                       f"rather than holding the full hoop tension (check r/r0).")
        else:
            verdict = f"OVER prediction (g_total/YL ≈ {mean_ratio:.2f}) — check for prestress contamination."

    res = {"config": {"n_fil": n_fil, "n_motors": args.n_motors, "turgor_Pa": turgor,
                      "R_cell_um": 7.5, "band_mN_m": list(BAND)},
           "conditions": conditions, "verdict": verdict}
    print(f"\n=== TURGOR VERDICT ===\n  {verdict}")
    _OUT.mkdir(parents=True, exist_ok=True)
    out = _OUT / "stage2_turgor_test.json"
    out.write_text(json.dumps(res, indent=2))
    print(f"\n[json] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
