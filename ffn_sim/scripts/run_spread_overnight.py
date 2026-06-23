"""Overnight proxy-free spreading run — the full mechanistic stack with the cohesion-bundle FORCE
fix (cad×40 / ecm×167, KB-anchored), measured with the peeling-aware spread_eval (so a 7-cell
basal-rim PEELING artifact is NOT misread as collective spread — the retracted A/A0→1.94 lesson).

NO proxy: wetting / substrate-well OFF; spreading must come from the ECM clutch (Pereverzev catch-
slip) + lamellipodium traction. NO outcome-tuning: every value is the committed lit/KB anchor.

Usage (on gbook A5000):
  python -m ffn_sim.scripts.run_spread_overnight --n 100 --steps 20000 --tag n100_proxyfree_spread
"""
from __future__ import annotations
import argparse
import json
import time
from ffn_sim.warp_port.dcm_warp_decohesion import run_decohesion

OUT = "/home/sungwook/ffn_phase_c/ffn_sim/outputs/warp_decohesion"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--settle", type=int, default=1000)
    ap.add_argument("--frames", type=int, default=60)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--tag", default="n100_proxyfree_spread")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--ipc", action="store_true",
                    help="M1: use IPC node-face contact (barrier+CCD) instead of the capped penalty — "
                         "the decisive validation that the strong cadherin/ecm bundle no longer tunnels "
                         "(penalty gave pen 3.1; IPC target pen<0.3, V/V0 stable, A/A0 unconfounded)")
    a = ap.parse_args()
    npz = f"{a.out}/{a.tag}.npz"
    t0 = time.time()
    out = run_decohesion(
        n_cells=a.n, subdiv=2, steps=a.steps, frames=a.frames, device=a.device,
        dt=8e-6, warmup=200, settle_steps=a.settle, gap=2.05,
        rep_strength=2e8, adh_strength=5e7,
        surface_tension=True, gamma_surf=1e-4, bending=True, edge_edge=True, nucleus=True,
        cadherin=True, ecm_clutch=True, lamellipodium=True,
        cad_bundle=40.0, ecm_bundle=167.0,
        substrate_wetting=False, use_substrate_well=False,   # NO proxy
        builder="fcc", integrator="implicit", accel_dt=8e-4,
        ipc=a.ipc,
        save_frames=npz)
    dt = time.time() - t0
    rec = {"tag": a.tag, "n": a.n, "steps": a.steps, "elapsed_s": round(dt, 1),
           "aa0_final": out.get("aa0_final"), "aa0_peak": out.get("aa0_peak"),
           "vv0_final": out.get("vv0_final"), "pen_final": out.get("pen_frac_final"),
           "pen_peak": out.get("pen_frac_peak"), "drift_um": out.get("drift_final_um"),
           "npz": npz}
    with open(f"{a.out}/{a.tag}.json", "w") as f:
        json.dump(rec, f, indent=2)
    print("SPREAD_RUN_DONE", json.dumps(rec))
    # peeling-aware verdict inline
    try:
        from ffn_sim.scripts.spread_eval import evaluate, verdict
        m = evaluate(npz)
        print("SPREAD_EVAL", verdict(m),
              f"aa0={m['aa0_final']:.3f} maxZdrop={m['maxZ_drop_pct']:.1f}% "
              f"peel_idx={m['peeling_index']:.1f} Psi={m['psi_init']:.3f}->{m['psi_final']:.3f}")
    except Exception as e:
        print("SPREAD_EVAL_FAILED", repr(e))


if __name__ == "__main__":
    main()
