#!/usr/bin/env python
"""KU-3.5 cortical-tension SATURATION-TIMESCALE diagnostic (H.3).

Purpose (PI 2026-05-28): the KU-3.5 production gate needs a myosin-active
cortex run to *steady state*, and the only way to know whether that is a
"now-feasible" local run or a constrained-BD-dependent one is to MEASURE
how many BAOAB steps the cortical tension takes to saturate.

This driver builds a cortex + myosin (+ ERM) Cell, runs it in bursts, and
records a panel of tension proxies vs timestep:

* mean cortex bead radius ⟨|r|⟩ and radius of gyration  (contraction)
* total potential energy                                  (network loading)
* scalar pressure + full pressure_tensor (config. stress) (in-plane tension)

It does NOT compute the final Young-Laplace γ — that calibration is the
next step IF the trajectory plateaus inside a feasible step budget. This
is a *diagnostic*, not the gate (no-gate-loosening: ratification stays
with the pytest gate once the run protocol is settled).

The myosin tension-buildup timescale is set by motor kinetics + network
relaxation, not by filament count, so the default ``smoke`` scale (300
filaments, full 100 minifilaments) characterises the timescale cheaply.

Usage:
    python ffn_sim/scripts/h3_ku35_tension_diag.py --scale smoke --max-steps 1500000
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import yaml
import hoomd
from hoomd import md

from ffn_sim.cell import Cell, CellBuildOptions
from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.cortex.erm import resolve_erm
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"

SCALES = {
    "smoke":  dict(n_fil=300),
    "medium": dict(n_fil=500),
    "full":   dict(n_fil=1000),
}


def _attach_thermo(sim: hoomd.Simulation, every: int) -> md.compute.ThermodynamicQuantities:
    """ThermodynamicQuantities + a devnull Table writer to force virial calc."""
    import os
    tq = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(tq)
    logger = hoomd.logging.Logger(categories=["scalar"])
    logger.add(tq, quantities=["pressure", "potential_energy"])
    writer = hoomd.write.Table(
        trigger=hoomd.trigger.Periodic(every),
        logger=logger,
        output=open(os.devnull, "w"),
    )
    sim.operations.writers.append(writer)
    return tq


def _cortex_radial_stats(sim: hoomd.Simulation, n_cortex_actin: int) -> tuple[float, float]:
    """Return (mean |r|, radius_of_gyration) of the cortex actin beads."""
    with sim.state.cpu_local_snapshot as snap:
        tag = np.asarray(snap.particles.tag)
        pos = np.asarray(snap.particles.position)
        inv = np.empty_like(tag)
        inv[tag] = np.arange(tag.size, dtype=tag.dtype)
        cortex = pos[inv][:n_cortex_actin]
    r = np.linalg.norm(cortex, axis=1)
    com = cortex.mean(axis=0)
    rg = float(np.sqrt(((cortex - com) ** 2).sum(axis=1).mean()))
    return float(r.mean()), rg


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scale", choices=SCALES, default="smoke")
    add_production_device_args(ap, default="gpu")
    ap.add_argument("--max-steps", type=int, default=1_500_000)
    ap.add_argument("--burst", type=int, default=25_000)
    ap.add_argument("--with-erm", action="store_true", default=True)
    ap.add_argument("--no-erm", dest="with_erm", action="store_false")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    with open(CFG) as f:
        cfg = yaml.safe_load(f)
    cfg["cortex"]["n_filaments"] = SCALES[args.scale]["n_fil"]
    cfg["cortex"]["demo_mode"] = True

    p = resolve_h3_derived(cfg)
    p_myo = resolve_cortex_myosin(cfg, dt=p.dt_cfl)
    p_erm = resolve_erm(cfg, kT=p.kT, R_cell=p.R_cell) if args.with_erm else None

    device = hoomd.device.GPU() if args.device == "gpu" else hoomd.device.CPU()
    opts = CellBuildOptions(
        with_baoab=True, with_myosin=True, with_erm=args.with_erm,
    )
    cell = Cell.build(
        p, p_myosin=p_myo, p_erm=p_erm, options=opts, device=device,
        rng=np.random.default_rng(p.seed),
    )
    sim = cell.simulation
    n_cortex_actin = cell.n_cortex_actin
    n_total = int(sim.state.N_particles)

    every = max(1, args.burst // 2)
    tq = _attach_thermo(sim, every)

    outdir = PKG / "outputs" / "h3" / "production"
    outdir.mkdir(parents=True, exist_ok=True)

    print(
        f"[start] KU-3.5 tension diag scale={args.scale} device={args.device} "
        f"n_cortex_actin={n_cortex_actin} n_total={n_total} "
        f"with_erm={args.with_erm} dt_cfl={p.dt_cfl:.3e}s "
        f"hoomd={hoomd.version.version} gpu_build={hoomd.version.gpu_enabled}",
        flush=True,
    )

    rows: list[dict] = []
    t0 = time.time()
    done = 0
    # Record an initial (pre-dynamics) baseline.
    sim.run(0)
    while done <= args.max_steps:
        rmean, rg = _cortex_radial_stats(sim, n_cortex_actin)
        try:
            ptensor = [float(x) for x in tq.pressure_tensor]
            pscalar = float(tq.pressure)
            pe = float(tq.potential_energy)
        except Exception:
            ptensor, pscalar, pe = [float("nan")] * 6, float("nan"), float("nan")
        sim_time_s = done * p.dt_cfl
        el = time.time() - t0
        row = dict(
            step=done, sim_time_s=sim_time_s,
            r_mean_um=rmean * 1e6, rg_um=rg * 1e6,
            potential_energy=pe, pressure_pa=pscalar,
            pressure_tensor=ptensor, wall_s=el,
        )
        rows.append(row)
        tps = done / el if el > 0 else 0.0
        print(
            f"[t={done:>9d} | {sim_time_s*1e3:7.3f} ms] "
            f"r_mean={rmean*1e6:6.3f}um rg={rg*1e6:6.3f}um "
            f"PE={pe:.4e} P={pscalar:+.4e}Pa "
            f"wall={el:6.1f}s tps={tps:6.0f}",
            flush=True,
        )
        if done >= args.max_steps:
            break
        step = min(args.burst, args.max_steps - done)
        sim.run(step)
        done += step

    result = dict(
        scale=args.scale, device=args.device,
        n_cortex_actin=int(n_cortex_actin), n_total=int(n_total),
        with_erm=bool(args.with_erm), dt_cfl_s=float(p.dt_cfl),
        max_steps=int(args.max_steps), burst=int(args.burst),
        wall_s=time.time() - t0, rows=rows,
    )
    out = outdir / f"ku35_tension_diag_{args.scale}_{args.device}.json"
    out.write_text(json.dumps(result, indent=2))
    print("RESULT_FILE " + str(out), flush=True)
    print("RESULT " + json.dumps({k: v for k, v in result.items() if k != "rows"}),
          flush=True)


if __name__ == "__main__":
    main()
