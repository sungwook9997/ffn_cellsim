#!/usr/bin/env python
"""Profile the constrained-BD production step to data-ground GPU porting.

PI 2026-05-31: before porting anything to GPU, quantify WHERE per-step wall time
goes — the prior finding (2026-05-30) was that SHAKE+Fixman dominate; this
confirms/quantifies it so the GPU port targets the dominant kernel, not a
sideshow (porting BAOAB-only would be wasted if SHAKE+Fixman dominate).

Method: build a production-ish constrained cortex (n_fil, n_motors), short
warm-up to relax overlaps, then cProfile a fixed number of constrained steps.
Reports: steps/sec, top Python functions by tottime, and the Python-attributed
vs HOOMD-native (C++ force/nlist) split. Also a myosin ON/OFF A/B to isolate
the binding-updater cost.

Usage: python -m ffn_sim.scripts.h3_profile_constrained
"""
from __future__ import annotations

import argparse
import cProfile
import io
import pstats
import time
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cell.cell import build_cortex_full_simulation
from ffn_sim.scripts.h3_ku35_tension import _tagpos

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"

# Python functions whose tottime we attribute to named cost buckets.
BUCKETS = {
    "shake": ("shake_project", "shake_project_chains"),
    "fixman": ("fixman_logdet_and_force",),
    "baoab": ("act",),            # BAOAB act (filtered to baoab.py below)
    "wrap": ("_wrap_into_box",),
    "binding(myosin/xlink)": ("query_ball_point", "cKDTree", "__init__"),
}


def _make_device(device: str):
    if device == "gpu":
        return hoomd.device.GPU(notice_level=0)
    return hoomd.device.CPU(notice_level=0)


def _build(n_fil, n_motors, seed, dt_factor, with_myosin, device="cpu"):
    cfg = yaml.safe_load(open(CFG))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    cfg["cortex"]["myosin"]["n_motors_per_cell"] = n_motors if with_myosin else 0
    p = resolve_h3_derived(cfg)
    nca = p.n_filaments * p.beads_per_filament
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = dt_factor * tau_bend
    # NB: profiling uses default binned_r0 + no force scaling, so R_cell is not
    # needed here (keeps this script runnable against older synced trees too).
    p_myo = resolve_cortex_myosin(cfg, dt=dtc) if with_myosin else None
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    dev = _make_device(device)
    return p, p_myo, p_xl, dtc, nca, dev


def _warm_and_constrain(p, p_myo, p_xl, dtc, dev, seed, n_warmup):
    hw = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
        constrained=False, rng=np.random.default_rng(seed))
    hw["sim"].run(0)
    hw["sim"].run(n_warmup)
    pos_warm = _tagpos(hw["sim"])
    del hw
    hc = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
        constrained=True, constrained_dt=dtc, rng=np.random.default_rng(seed))
    sim = hc["sim"]
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = pos_warm
    sim.state.set_snapshot(snap)
    sim.run(0)
    return sim


def _timed(sim, steps):
    t0 = time.time()
    sim.run(steps)
    return time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=150)
    ap.add_argument("--n-motors", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dt-factor", type=float, default=0.001)
    ap.add_argument("--n-warmup", type=int, default=8000)
    ap.add_argument("--profile-steps", type=int, default=1500)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--no-cprofile", action="store_true",
                    help="skip the cProfile pass (steps/s timing only — for GPU smoke)")
    args = ap.parse_args()

    print(f"=== constrained-BD profile: device={args.device} n_fil={args.n_fil} "
          f"n_motors={args.n_motors} profile_steps={args.profile_steps} "
          f"hoomd_gpu_build={hoomd.version.gpu_enabled} ===", flush=True)

    # --- A/B: myosin ON vs OFF (isolate binding-updater cost) ---
    for with_myo in (True, False):
        p, p_myo, p_xl, dtc, nca, dev = _build(
            args.n_fil, args.n_motors, args.seed, args.dt_factor, with_myo,
            device=args.device)
        sim = _warm_and_constrain(p, p_myo, p_xl, dtc, dev, args.seed, args.n_warmup)
        _timed(sim, 100)  # warm the step path
        wall = _timed(sim, args.profile_steps)
        sps = args.profile_steps / wall
        label = "myosin ON " if with_myo else "myosin OFF"
        print(f"[{label}] {sps:7.1f} steps/s   ({wall:.2f}s / {args.profile_steps} steps)",
              flush=True)
        if with_myo:
            sim_for_profile = (p, p_myo, p_xl, dtc, dev)

    if args.no_cprofile:
        return
    # --- cProfile the myosin-ON constrained run (Python-side breakdown) ---
    p, p_myo, p_xl, dtc, dev = sim_for_profile
    sim = _warm_and_constrain(p, p_myo, p_xl, dtc, dev, args.seed, args.n_warmup)
    _timed(sim, 100)
    pr = cProfile.Profile()
    t0 = time.time()
    pr.enable()
    sim.run(args.profile_steps)
    pr.disable()
    wall = time.time() - t0

    st = pstats.Stats(pr)
    # Total Python tottime across all funcs (≈ host-side cost; remainder = C++).
    total_tottime = sum(v[2] for v in st.stats.values())
    print(f"\n--- cProfile ({args.profile_steps} steps, wall {wall:.2f}s) ---", flush=True)
    print(f"Python-attributed tottime = {total_tottime:.2f}s "
          f"({100*total_tottime/wall:.0f}% of wall); "
          f"HOOMD-native (C++ force/nlist/integrate) ≈ {wall-total_tottime:.2f}s "
          f"({100*(wall-total_tottime)/wall:.0f}%)", flush=True)

    # Bucket the named cost functions.
    print("\n--- named-function tottime (host-side hotspots) ---", flush=True)
    rows = []
    for (fpath, line, fname), v in st.stats.items():
        tot = v[2]
        if tot < 0.02:
            continue
        rows.append((tot, fname, Path(fpath).name, line))
    rows.sort(reverse=True)
    for tot, fname, fname_file, line in rows[:20]:
        print(f"  {tot:7.2f}s  {100*tot/wall:4.0f}%  {fname}  ({fname_file}:{line})",
              flush=True)


if __name__ == "__main__":
    main()
