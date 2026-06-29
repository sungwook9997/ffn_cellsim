"""H.7 hot-loop scaling benchmark.

This is a timing harness, not a validation gate. It measures two ceilings:

* ``baoab-fast``: the current fastest Python/CuPy Gate-A path
  (device BAOAB, xlink frozen, snapshot-free myosin).
* ``native-dt0``: HOOMD built-in force-stack scheduling with ``dt=0`` and
  ``kT=0`` so particles do not move. This is a timing ceiling for the current
  explicit topology, not a physical run.

The script records failures instead of hiding them. Large ``n_fil`` values may
fail the physical BAOAB path before speed is measurable; that is itself useful
production information.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from copy import deepcopy
from pathlib import Path
from typing import Any

import hoomd
import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin
from ffn_sim.scripts.h7_gate_a import run as run_gate_a
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"


def _device(name: str):
    if name == "gpu":
        return hoomd.device.GPU(notice_level=0)
    return hoomd.device.CPU(notice_level=0)


def _resolved(n_fil: int, batch_steps: int):
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = int(n_fil)
    cfg["cortex"]["demo_mode"] = True
    cfg.setdefault("cortex", {}).setdefault("myosin", {})
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = True
    cfg["cortex"]["myosin"]["batch_steps"] = int(batch_steps)
    p_cortex = resolve_h3_derived(cfg)
    p_myo = resolve_cortex_myosin(
        cfg, dt=p_cortex.dt_cfl, R_cell=p_cortex.R_cell
    )
    p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
    return p_cortex, p_myo, p_xl


def _count_topology(sim: hoomd.Simulation) -> dict[str, int]:
    snap = sim.state.get_snapshot()
    if snap.communicator.rank != 0:
        return {}
    return {
        "particles": int(snap.particles.N),
        "bonds": int(snap.bonds.N),
        "angles": int(snap.angles.N),
        "particle_types": len(list(snap.particles.types)),
        "bond_types": len(list(snap.bonds.types)),
    }


def _bench_native_dt0(
    *,
    n_fil: int,
    steps: int,
    batch_steps: int,
    device_name: str,
    seed: int,
) -> dict[str, Any]:
    p_cortex, p_myo, p_xl = _resolved(n_fil, batch_steps)
    build_t0 = time.time()
    hw = build_cortex_full_simulation(
        p_cortex,
        p_xlinks=p_xl,
        p_myosin=p_myo,
        device=_device(device_name),
        with_baoab=False,
        connected_mesh=True,
        rng=np.random.default_rng(seed),
    )
    sim = hw["sim"]
    while len(sim.operations.updaters) > 0:
        sim.operations.updaters.remove(sim.operations.updaters[0])
    sim.operations.integrator.dt = 0.0
    sim.operations.integrator.methods.append(
        hoomd.md.methods.Brownian(
            filter=hoomd.filter.All(),
            kT=0.0,
            default_gamma=p_cortex.gamma_b,
        )
    )
    build_wall = time.time() - build_t0
    sim.run(0)
    topo = _count_topology(sim)
    run_t0 = time.time()
    sim.run(int(steps))
    wall = time.time() - run_t0
    return {
        "mode": "native-dt0",
        "n_fil": int(n_fil),
        "steps": int(steps),
        "wall_s": wall,
        "steps_per_s": float(steps / wall) if wall > 0 else float("nan"),
        "build_wall_s": build_wall,
        "topology": topo,
        "forces": [type(f).__name__ for f in sim.operations.integrator.forces],
        "device": device_name,
    }


def _bench_baoab_fast(
    *,
    n_fil: int,
    steps: int,
    batch_steps: int,
    device_name: str,
    seed: int,
) -> dict[str, Any]:
    if device_name == "gpu":
        os.environ.setdefault("FFN_GPU_DEVICE_BAOAB", "1")
    ticks = max(1, int(round(steps / batch_steps)))
    res = run_gate_a(
        n_fil=int(n_fil),
        batch_steps=int(batch_steps),
        ticks=ticks,
        measure_every=ticks,
        seed=int(seed),
        device=_device(device_name),
        connected_mesh=True,
        measure_mode="final",
        freeze_xlinks=True,
        snapshot_free_myosin=True,
        myosin_warm_bind_ticks=0,
        snapshot_free_grip_clock=True,
    )
    return {
        "mode": "baoab-fast",
        "n_fil": int(n_fil),
        "requested_steps": int(steps),
        "batch_steps": int(batch_steps),
        "ticks": ticks,
        "wall_s": res["wall_s"],
        "steps_per_s": res["steps_per_s"],
        "simulated_baoab_steps": res["simulated_baoab_steps"],
        "physics_mode": res["physics_mode"],
        "verdict": res["verdict"],
        "last_row": res["rows"][-1],
        "device": device_name,
    }


def _run_one(mode: str, **kwargs) -> dict[str, Any]:
    try:
        if mode == "native-dt0":
            return _bench_native_dt0(**kwargs)
        if mode == "baoab-fast":
            return _bench_baoab_fast(**kwargs)
        raise ValueError(f"unknown mode {mode!r}")
    except Exception as exc:  # noqa: BLE001 - benchmark must record failures
        return {
            "mode": mode,
            "n_fil": int(kwargs["n_fil"]),
            "status": "failed",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "steps": int(kwargs["steps"]),
            "batch_steps": int(kwargs["batch_steps"]),
            "device": kwargs["device_name"],
        }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, nargs="+", default=[150, 1000, 3000])
    ap.add_argument(
        "--mode",
        choices=("baoab-fast", "native-dt0"),
        nargs="+",
        default=["baoab-fast", "native-dt0"],
    )
    ap.add_argument("--steps", type=int, default=60000)
    ap.add_argument("--batch-steps", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument(
        "--json",
        default=str(PKG / "outputs" / "h7" / "production" / "hotloop_scaling.json"),
    )
    add_production_device_args(ap, default="gpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)

    rows = []
    for n_fil in args.n_fil:
        for mode in args.mode:
            print(f"RUN mode={mode} n_fil={n_fil}", flush=True)
            row = _run_one(
                mode,
                n_fil=n_fil,
                steps=args.steps,
                batch_steps=args.batch_steps,
                device_name=args.device,
                seed=args.seed,
            )
            rows.append(row)
            if row.get("status") == "failed":
                print(
                    f"  FAIL {row['error_type']}: {row['error']}",
                    flush=True,
                )
            else:
                print(
                    f"  steps/s={row['steps_per_s']:.1f} wall={row['wall_s']:.2f}s",
                    flush=True,
                )

    out = Path(args.json)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "rows": rows,
        "steps": int(args.steps),
        "batch_steps": int(args.batch_steps),
        "device": args.device,
        "seed": int(args.seed),
    }
    json.dump(payload, open(out, "w"), indent=2)
    print(f"WROTE {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
