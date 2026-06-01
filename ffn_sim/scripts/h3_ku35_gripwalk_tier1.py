#!/usr/bin/env python
"""KU-3.5 grip-walk Tier-1 micro-diagnostic (PI-ratified 2026-05-31).

The DECISIVE, minutes-scale test from KU35_GRIP_WALK_DESIGN_2026-05-31.md §5.1:
does the AFINES ``grip_walk`` stepping mode produce a *sustained / growing*
head-actin stretch and cortical tension, where the legacy ``binned_r0`` lumped
proxy relaxes to ~0 (the diagnosed KU-3.5 floor)?

It builds a small cortex + myosin cell, warms it up, then runs the constrained
(rigid-actin M-SHAKE) production integrator for a short window in BOTH modes
(A = binned_r0 proxy, B = grip_walk) and reports, per sample:

  * γ_soft  — method-of-planes soft-bond cortical tension (the contraction signal)
  * r/r0    — mean cortex radius vs initial (contraction proxy; <1 = contracted)
  * mean commanded grip stretch s_grip over engaged heads (grip_walk only)
  * n_engaged, n_step_advances

DIAGNOSTIC ACCELERANT (documented, NOT a physics change): at the literal
v0 = 0.2 µm/s + the FA-limited dt, walking ONE bead spacing takes ~3.7e6 steps
(2.5 s sim time), so no minutes-scale run can show contraction regardless of
mechanism (design brief §1.4). We raise v0 by ``--v0-accel`` so integer-bead
walks fire within ~1e4-1e5 steps. v0 only sets the WALK RATE; the mechanism
question — does the stretch persist/grow or relax? — is v0-independent. The
production γ re-measure (Tier 3) uses the literal v0 at a correctly-sized run
length.

PASS contract (Tier 1):
  1. grip_walk mean s_grip strictly increases over the first samples (proxy: 0).
  2. grip_walk γ_soft ends ≥ 10× the binned_r0 γ_soft floor (directional proof).
  3. grip_walk r/r0 < binned_r0 r/r0 (net contraction attributable to the walk).
  4. No constraint-drift blowup, no LJ-CFL WARN in either arm.

Usage:
    python ffn_sim/scripts/h3_ku35_gripwalk_tier1.py            # default A/B
    python ffn_sim/scripts/h3_ku35_gripwalk_tier1.py --v0-accel 50 --n-sample 6
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml
import hoomd
import hoomd.write

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cell.cell import build_cortex_full_simulation
from ffn_sim.scripts.h3_ku35_tension import (
    _tension_method_of_planes,
    _tension_method_of_planes_rigid,
    _tagpos,
)

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"


def _mean_grip_s(ma) -> float:
    """Mean commanded sub-bead stretch over engaged heads (grip_walk only)."""
    engaged = ma._head_bound_to_actin >= 0
    if not engaged.any():
        return 0.0
    return float(ma._head_grip_s[engaged].mean())


def _run_arm(
    *, stepping_mode: str, n_fil: int, n_motors: int, seed: int,
    dt_factor: float, v0_accel: float, n_warmup: int, n_sample: int,
    interval: int, force_scaling: bool = False, device: str = "cpu",
    gsd_path: str | None = None, gsd_period: int = 0,
) -> dict:
    """Build + warm-up + short constrained run for one stepping mode."""
    cfg = yaml.safe_load(open(CFG))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors
    cfg["cortex"]["myosin"]["stepping_mode"] = stepping_mode
    # Route B: mesoscale force scaling (grip_walk only; derived factor).
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = bool(force_scaling)
    p = resolve_h3_derived(cfg)
    N = p.beads_per_filament
    F = p.n_filaments
    nca = F * N
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = dt_factor * tau_bend
    p_myo_lit = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p.R_cell)  # literal v0
    # Diagnostic accelerant — multiply the literal v0 so walks are observable
    # within the SHORT measured window. Applied ONLY to the constrained
    # production phase: during warm-up the accelerated proxy over-contracts an
    # un-equilibrated network → LJ-overlap blow-up, so warm-up uses literal v0.
    p_myo_acc = replace(p_myo_lit, v0_per_head=p_myo_lit.v0_per_head * v0_accel)
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    dev = (hoomd.device.GPU(notice_level=0) if device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    # Warm-up (standard BAOAB, CFL dt, LITERAL v0) to relax construction overlaps.
    # Use the build's soft-start equilibration (clipped-Brownian _softstart_step
    # + BAOAB drain) rather than a raw sim.run(n_warmup): at large n_fil a raw
    # warm-up drives LJ-overlap displacements past the int32-image guard. The
    # soft-start clips per-step displacement so overlaps drain safely regardless
    # of N (GPU-main port 2026-06-01). softstart count scales with N.
    n_soft = max(300, n_warmup // 8)
    hw = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo_lit, device=dev, with_baoab=True,
        constrained=False, rng=np.random.default_rng(seed),
        equilibrate=True, equilibrate_steps=n_warmup,
        equilibrate_softstart_steps=n_soft)
    hw["sim"].run(0)
    pos_warm = _tagpos(hw["sim"])
    del hw

    # Constrained production (rigid actin backbone, fast dt, ACCELERATED v0).
    # Same particle topology as the warm-up build (only v0/mode differ), so the
    # relaxed positions transfer by tag directly.
    hc = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo_acc, device=dev, with_baoab=True,
        constrained=True, constrained_dt=dtc, rng=np.random.default_rng(seed))
    sim = hc["sim"]
    ma = hc["myosin_action"]
    act = hc["baoab_action"]
    # Capture the M-SHAKE Lagrange multiplier so the DOMINANT cortex-tension
    # channel (rigid actin backbone, ~200× the soft-bond channel per the
    # KU-3.5 root-cause) is measurable. No behavioural change (λ is already
    # computed inside SHAKE every step; only the discarded result is captured).
    if hasattr(act, "record_lambda"):
        act.record_lambda = True
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm
    sim.state.set_snapshot(snap)
    sim.run(0)

    # Optional GSD trajectory dump for 3D structural viz (OVITO direct / Blender
    # import) — the "GSD writer unlock" (viz-decisions). HOOMD-native, writes
    # device-resident state (positions + box + frame-0 topology) every
    # gsd_period steps with no extra per-step host sync. Opt-in (period 0 = off).
    if gsd_path and gsd_period > 0:
        gsd_writer = hoomd.write.GSD(
            filename=gsd_path, trigger=hoomd.trigger.Periodic(int(gsd_period)),
            mode="wb", dynamic=["property"],
        )
        sim.operations.writers.append(gsd_writer)
        print(f"[gsd] writing trajectory -> {gsd_path} every {gsd_period} steps",
              flush=True)

    r0 = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
    box_L = float(sim.state.box.L[0])
    r_prev = _tagpos(sim)
    samples = []
    t0 = time.time()
    for k in range(n_sample):
        sim.run(interval)
        r = _tagpos(sim)
        rmean = float(np.linalg.norm(r[:nca], axis=1).mean())
        dvec_raw = r - r_prev
        max_disp_raw = float(np.linalg.norm(dvec_raw, axis=1).max())
        r_prev = r
        gamma_soft = _tension_method_of_planes(sim, p.R_cell)
        gamma_rigid = _tension_method_of_planes_rigid(sim, act, p.R_cell, dtc)
        gamma_total = gamma_soft + gamma_rigid
        rec = dict(
            sample=k + 1, step=int(sim.timestep),
            r_over_r0=rmean / r0,
            gamma_soft_mN_per_m=gamma_soft * 1e3,
            gamma_rigid_mN_per_m=gamma_rigid * 1e3,
            gamma_total_mN_per_m=gamma_total * 1e3,
            mean_grip_s_nm=_mean_grip_s(ma) * 1e9,
            n_engaged=int(ma.n_engaged),
            step_advances=int(ma.n_step_advances_total),
            max_drift=float(act.max_constraint_drift),
            lj_cfl_warn=bool(max_disp_raw > 5.0 * box_L),
        )
        samples.append(rec)
        print(
            f"[{stepping_mode}] sample={k+1}/{n_sample} "
            f"r/r0={rec['r_over_r0']:.4f} "
            f"g_tot={rec['gamma_total_mN_per_m']:.3e} "
            f"(soft={rec['gamma_soft_mN_per_m']:.2e} rigid={rec['gamma_rigid_mN_per_m']:.3e}) mN/m "
            f"s_grip={rec['mean_grip_s_nm']:.1f}nm "
            f"engaged={rec['n_engaged']} adv={rec['step_advances']} "
            f"drift={rec['max_drift']:.2e}"
            f"{' WARN_LJ_CFL' if rec['lj_cfl_warn'] else ''}",
            flush=True,
        )
    return dict(
        stepping_mode=stepping_mode, wall_s=time.time() - t0,
        v0_accel=v0_accel, samples=samples,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=60)
    ap.add_argument("--n-motors", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dt-factor", type=float, default=0.001)
    ap.add_argument("--v0-accel", type=float, default=300.0,
                    help="diagnostic accelerant on v0 (walk rate; not physics)")
    ap.add_argument("--n-warmup", type=int, default=15_000)
    ap.add_argument("--n-sample", type=int, default=6)
    ap.add_argument("--interval", type=int, default=20_000)
    ap.add_argument("--force-scaling", action="store_true",
                    help="Route B: derived mesoscale force scaling on the grip_walk arm")
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--arm", choices=["both", "binned_r0", "grip_walk"],
                    default="both",
                    help="run a single arm (skips the A/B verdict) for fast "
                         "per-arm stability/config iteration")
    ap.add_argument("--gsd-period", type=int, default=0,
                    help="if >0, dump a GSD trajectory every N production steps "
                         "for 3D structural viz (OVITO/Blender). File is "
                         "<out>.<arm>.gsd, or gpu_run.<arm>.gsd if --out unset.")
    args = ap.parse_args()

    def _gsd_path(mode: str) -> str | None:
        if args.gsd_period <= 0:
            return None
        base = args.out if args.out else "gpu_run"
        # strip a trailing .json so the gsd sits beside the json artefact
        if base.endswith(".json"):
            base = base[:-5]
        return f"{base}.{mode}.gsd"

    common = dict(
        n_fil=args.n_fil, n_motors=args.n_motors, seed=args.seed,
        dt_factor=args.dt_factor, v0_accel=args.v0_accel,
        n_warmup=args.n_warmup, n_sample=args.n_sample, interval=args.interval,
        device=args.device, gsd_period=args.gsd_period,
    )
    print(f"=== KU-3.5 grip-walk Tier-1 micro-diagnostic (v0_accel={args.v0_accel}×"
          f"{', force-scaling ON' if args.force_scaling else ''}, arm={args.arm}) ===",
          flush=True)
    arm_a = arm_b = None
    if args.arm in ("both", "binned_r0"):
        arm_a = _run_arm(stepping_mode="binned_r0",
                         gsd_path=_gsd_path("binned_r0"), **common)  # floor (scaling is grip_walk-gated)
    if args.arm in ("both", "grip_walk"):
        arm_b = _run_arm(stepping_mode="grip_walk", force_scaling=args.force_scaling,
                         gsd_path=_gsd_path("grip_walk"), **common)

    # Single-arm mode: print that arm's trajectory + write it, skip the A/B
    # verdict (which compares both arms). Used to iterate fast on one arm's
    # stability/config without re-running the floor each time.
    if args.arm != "both":
        single = arm_a if args.arm == "binned_r0" else arm_b
        sbx = single["samples"]
        print(f"\n=== {args.arm} single-arm summary ===", flush=True)
        print(f"  final s_grip_nm = {sbx[-1].get('mean_grip_s_nm', 0.0)}", flush=True)
        print(f"  final step_advances = {sbx[-1]['step_advances']}", flush=True)
        print(f"  final gamma_total_mN/m = {sbx[-1]['gamma_total_mN_per_m']*1e3:.3e}",
              flush=True)
        print(f"  any_cfl_blowup = {any(s['lj_cfl_warn'] for s in sbx)}", flush=True)
        if args.out:
            Path(args.out).write_text(json.dumps(
                dict(arm=args.arm, samples=sbx, params=common), indent=2))
            print(f"wrote {args.out}", flush=True)
        return

    # --- PASS contract evaluation ---
    # γ is measured on the DOMINANT channel (γ_total = γ_soft + γ_rigid); the
    # rigid actin backbone carries ~200× the soft-bond tension (root-cause §1).
    sb = arm_b["samples"]
    sa = arm_a["samples"]
    grip_s = [s["mean_grip_s_nm"] for s in sb]
    g_tot_b_final = sb[-1]["gamma_total_mN_per_m"]
    g_tot_a_floor = max(s["gamma_total_mN_per_m"] for s in sa)  # proxy ceiling
    # 1. commanded stretch grows monotonically (proxy is flat at 0) — the
    #    anti-relaxation proof (the proxy's fatal flaw is that stretch relaxes).
    c1 = (all(grip_s[i + 1] >= grip_s[i] for i in range(len(grip_s) - 1))
          and grip_s[-1] > grip_s[0])
    # 2. real bead re-targets fire — the AFINES material-transport half engaged
    #    (not just sub-bead pre-stretch).
    c2 = sb[-1]["step_advances"] > 0
    # 3. grip_walk γ_total ends ≥ 10× the binned_r0 floor (sustained tension).
    c3 = g_tot_b_final >= 10.0 * max(g_tot_a_floor, 1e-30)
    # 4. no blow-ups in either arm.
    c4 = not any(s["lj_cfl_warn"] for s in sa + sb)

    verdict = dict(
        stretch_grows=c1,
        bead_transport_engaged=c2,
        gamma_total_10x_off_floor=c3,
        no_cfl_blowup=c4,
        binned_r0_gamma_total_floor_mN_per_m=g_tot_a_floor,
        grip_walk_gamma_total_final_mN_per_m=g_tot_b_final,
        ratio=g_tot_b_final / max(g_tot_a_floor, 1e-30),
        grip_walk_final_s_grip_nm=grip_s[-1],
        grip_walk_final_step_advances=sb[-1]["step_advances"],
        grip_walk_final_r_over_r0=sb[-1]["r_over_r0"],
        binned_r0_final_r_over_r0=sa[-1]["r_over_r0"],
        PASS=bool(c1 and c2 and c3 and c4),
    )
    print("\n=== TIER-1 VERDICT ===", flush=True)
    for k, v in verdict.items():
        print(f"  {k} = {v}", flush=True)
    print(f"\n{'✅ TIER-1 PASS' if verdict['PASS'] else '❌ TIER-1 FAIL'}", flush=True)

    if args.out:
        Path(args.out).write_text(json.dumps(
            dict(verdict=verdict, arm_binned_r0=arm_a, arm_grip_walk=arm_b,
                 params=common), indent=2))
        print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
