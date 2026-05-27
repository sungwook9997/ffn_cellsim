#!/usr/bin/env python
"""KU-3.5 — myosin-active cortex under constrained-BD (foundational driver).

Builds a cortex + myosin cell, warms it up with the standard BAOAB at the
CFL dt to relax construction overlaps, then hands the relaxed configuration
to the constrained integrator (rigid actin backbone, M-SHAKE + Fixman) for a
fast-dt production run where myosin actively contracts the network.

Saves the cortex-actin trajectory (for visualisation) plus per-snapshot
diagnostics: myosin engagement, total binding/stepping, mean cortex radius
(a contraction proxy), constraint drift.

NOTE — cortical tension γ: the rigorous Young-Laplace / method-of-planes
measurement with RIGID constraints (the actin "tension" is a Lagrange
multiplier, not a harmonic-bond force) is a Sanity-Gate measurement-protocol
decision pending PI ratification (like the L_p band). This driver reports a
CONTRACTION PROXY (mean radius vs time) and myosin engagement, NOT a
ratified γ — KU-3.5 PASS is not declared here.

ERM is OFF by default: at the LJ-limited production dt (~0.001·τ_bend = 0.7 μs)
the ratified k_ERM=1e-4 (τ_ERM=3.9 μs) re-violates CFL; the k_ERM for this dt
is a PI item. Myosin-active contraction is demonstrable without ERM.

Usage:
    python ffn_sim/scripts/h3_ku35_tension.py --n-fil 150 --dt-factor 0.001
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.cell.cell import build_cortex_full_simulation

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"


def _tagpos(sim) -> np.ndarray:
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position); tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        return pos[inv].copy()


def run(n_fil: int, seed: int, *, dt_factor: float = 0.001,
        n_warmup: int = 40_000, n_sample: int = 80, interval: int = 5_000) -> dict:
    cfg = yaml.safe_load(open(CFG))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    N = p.beads_per_filament; F = p.n_filaments; nca = F * N
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = dt_factor * tau_bend
    p_myo = resolve_cortex_myosin(cfg, dt=dtc)
    dev = hoomd.device.CPU(notice_level=0)

    # Phase 1 — warm-up (standard BAOAB, CFL dt) to relax overlaps.
    hw = build_cortex_full_simulation(
        p, p_myosin=p_myo, device=dev, with_baoab=True, constrained=False,
        rng=np.random.default_rng(seed))
    hw["sim"].run(0); hw["sim"].run(n_warmup)
    pos_warm = _tagpos(hw["sim"])
    del hw

    # Phase 2 — constrained production (rigid actin backbone, fast dt).
    hc = build_cortex_full_simulation(
        p, p_myosin=p_myo, device=dev, with_baoab=True, constrained=True,
        constrained_dt=dtc, rng=np.random.default_rng(seed))
    sim = hc["sim"]; ma = hc["myosin_action"]; act = hc["baoab_action"]
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm
    sim.state.set_snapshot(snap)
    sim.run(0)

    frames = np.empty((n_sample, F, N, 3))
    diag = []
    r0 = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
    for k in range(n_sample):
        sim.run(interval)
        r = _tagpos(sim)
        frames[k] = r[:nca].reshape(F, N, 3)
        rmean = float(np.linalg.norm(r[:nca], axis=1).mean())
        diag.append(dict(
            step=int(sim.timestep), r_cortex_um=rmean * 1e6,
            r_over_r0=rmean / r0, myosin_engaged=int(ma.n_engaged),
            bind_total=int(ma.n_bind_total), step_advances=int(ma.n_step_advances_total),
            max_drift=float(act.max_constraint_drift),
        ))

    outdir = PKG / "outputs" / "h3" / "production" / "ku35"
    outdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(outdir / f"ku35_frames_n{n_fil}_s{seed}.npz", frames=frames)
    result = dict(
        n_fil=int(F), seed=int(seed), dt_s=float(dtc), dt_factor=float(dt_factor),
        dt_speedup_vs_cfl=float(dtc / p.dt_cfl), r0_um=r0 * 1e6,
        r_final_over_r0=diag[-1]["r_over_r0"], myosin_engaged_final=diag[-1]["myosin_engaged"],
        n_motors=int(p_myo.n_motors_per_cell), max_drift=diag[-1]["max_drift"],
        note="contraction proxy + myosin engagement; ratified gamma pending PI protocol",
        diag=diag,
    )
    (outdir / f"ku35_n{n_fil}_s{seed}.json").write_text(json.dumps(result, indent=2))
    print("RESULT " + json.dumps({k: v for k, v in result.items() if k != "diag"}), flush=True)
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dt-factor", type=float, default=0.001)
    ap.add_argument("--n-warmup", type=int, default=40_000)
    ap.add_argument("--n-sample", type=int, default=80)
    ap.add_argument("--interval", type=int, default=5_000)
    args = ap.parse_args()
    run(args.n_fil, args.seed, dt_factor=args.dt_factor,
        n_warmup=args.n_warmup, n_sample=args.n_sample, interval=args.interval)


if __name__ == "__main__":
    main()
