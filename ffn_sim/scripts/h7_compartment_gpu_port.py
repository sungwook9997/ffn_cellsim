"""H.7 compartment-force GPU port — throughput + γ-parity matrix (full cell).

The dominant full-cell step cost at the real operating point is the three
compartment ``md.force.Custom`` (turgor / membrane-tension / nucleus-confinement)
reading positions via ``cpu_local_snapshot`` every step (a GPU→host sync;
≈72% of the native-integrator step — H7_NATIVE_FULLCELL_GO_2026-06-07.md). Their
device-resident twins (``*_gpu.py``, gpu_local_snapshot + cupy, NO host-sync)
already exist; this driver wires them on (``FFN_GPU_DEVICE_COMPARTMENTS=1``) and
measures the win, with γ parity to prove identical physics.

Arms (same stable two-phase build, identical config/seed/steps):
  * ``cupy_cpucomp`` — cupy integrator + CPU compartment forces (current default).
  * ``cupy_gpucomp`` — cupy integrator + GPU compartment forces (compartment port
    alone).
  * ``native_gpucomp`` — native constrained integrator + GPU compartment forces
    (the full GPU-main stack).

Reports steps/s + 2×10⁸ ETA per arm + γ (soft+rigid) parity vs the baseline.
GPU-only. ``PYTHONPATH`` must include the native plugin build+python for the
native arm.

Usage (gbook):
    export PYTHONPATH="$PWD/native/ffn_hoomd_plugin/build:$PWD/native/ffn_hoomd_plugin/python:$PYTHONPATH"
    python -m ffn_sim.scripts.h7_compartment_gpu_port --n-fil 1000 --warmup 2000 \
        --steps 1500 --interval 200 --device gpu
"""

from __future__ import annotations

import argparse
import json
import os
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.scripts.h3_ku35_tension import (
    _tension_method_of_planes,
    _tension_method_of_planes_rigid,
)
from ffn_sim.scripts.h7_native_fullcell_go import (
    _NativeLambdaAdapter,
    _tagpos,
    _to_host,
)

_OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "production"


def _build_constrained(manifest, *, device, seed, warm_pos):
    cell = build_baseline_cell(
        manifest=deepcopy(manifest), device=device, seed=seed,
        constrained=True, equilibrate=False,
    )
    if hasattr(cell.baoab_action, "record_lambda"):
        cell.baoab_action.record_lambda = True
    snap = cell.simulation.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = warm_pos
    cell.simulation.state.set_snapshot(snap)
    cell.simulation.run(0)
    return cell


def _time_arm(cell, *, steps, interval, R_cell, dtc, rigid_action):
    sim = cell.simulation
    sim.run(interval)  # warm / JIT
    t0 = time.perf_counter()
    sim.run(int(steps))
    sps = steps / (time.perf_counter() - t0)
    samples = []
    for _ in range(5):
        sim.run(interval)
        samples.append((
            float(_tension_method_of_planes(sim, R_cell)),
            float(_tension_method_of_planes_rigid(sim, rigid_action, R_cell, dtc)),
        ))
    arr = np.asarray(samples)
    return sps, float(arr[:, 0].mean()), float(arr[:, 1].mean())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--warmup", type=int, default=2000)
    ap.add_argument("--steps", type=int, default=1500)
    ap.add_argument("--interval", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    ap.add_argument("--allow-cpu-dev", action="store_true")
    ap.add_argument("--gate-steps", type=float, default=2.0e8)
    args = ap.parse_args()

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    manifest.setdefault("cortex_overrides", {}).setdefault(
        "cortex", {})["n_filaments"] = int(args.n_fil)
    softstart = max(300, args.warmup // 8)

    # ---- warm-up once (unconstrained; compartments parity-equal) ----
    os.environ.pop("FFN_GPU_DEVICE_COMPARTMENTS", None)
    cell_w = build_baseline_cell(
        manifest=deepcopy(manifest), device=dev, seed=args.seed,
        constrained=False, equilibrate=True, equilibrate_steps=args.warmup,
        equilibrate_softstart_steps=softstart,
    )
    cell_w.simulation.run(0)
    warm_pos = _tagpos(cell_w.simulation)
    R_cell = float(cell_w.p_cortex.R_cell)
    n_part = int(cell_w.simulation.state.N_particles)
    del cell_w

    results = {}

    def _run(label, *, gpu_comp, native_int):
        if gpu_comp:
            os.environ["FFN_GPU_DEVICE_COMPARTMENTS"] = "1"
        else:
            os.environ.pop("FFN_GPU_DEVICE_COMPARTMENTS", None)
        cell = _build_constrained(manifest, device=dev, seed=args.seed, warm_pos=warm_pos)
        sim = cell.simulation
        act = cell.baoab_action
        dtc = float(act.dt)
        rigid_action = act
        if native_int:
            from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater
            chains_stacked = np.asarray(act.chains_tag_stacked)
            sim.run(0)  # ensure inv_gamma populated
            inv_gamma = [float(x) for x in np.asarray(
                _to_host(act._inv_gamma_by_tag), dtype=np.float64)]
            F, npc = chains_stacked.shape
            sim.operations.updaters.remove(cell.baoab_updater)
            nat = NativeConstrainedBaoabUpdater(
                dt=dtc, kT=float(act.kT), inv_gamma_by_tag=inv_gamma,
                chains_tag=chains_stacked.reshape(-1).astype(np.int32),
                n_chains=F, bonds_per_chain=npc - 1,
                rest_length=float(act._chain_rest_length),
                seed=int(act._seed), tol=1.0e-9, max_iter=200,
                trigger=hoomd.trigger.Periodic(1),
            )
            sim.operations.updaters.append(nat)
            sim.run(0)
            rigid_action = _NativeLambdaAdapter(
                nat, chains_stacked, float(act._chain_rest_length))
        sps, g_soft, g_rigid = _time_arm(
            cell, steps=args.steps, interval=args.interval,
            R_cell=R_cell, dtc=dtc, rigid_action=rigid_action)
        nonconv = int(nat.nonconverged_count) if native_int else None
        results[label] = {
            "steps_per_s": sps,
            "eta_hours_at_gate": args.gate_steps / sps / 3600.0 if sps > 0 else None,
            "g_soft_mN_m": g_soft * 1e3, "g_rigid_mN_m": g_rigid * 1e3,
            "nonconverged_count": nonconv,
        }
        del cell
        os.environ.pop("FFN_GPU_DEVICE_COMPARTMENTS", None)

    _run("cupy_cpucomp", gpu_comp=False, native_int=False)
    _run("cupy_gpucomp", gpu_comp=True, native_int=False)
    native_avail = True
    try:
        import ffn_hoomd_plugin  # noqa: F401
    except Exception as exc:  # noqa: BLE001
        native_avail = False
        results["native_gpucomp"] = {"error": f"native unavailable: {exc!r}"}
    if native_avail:
        _run("native_gpucomp", gpu_comp=True, native_int=True)

    base = results["cupy_cpucomp"]["steps_per_s"]
    summary = {
        "config": {"n_fil": args.n_fil, "n_part": n_part, "R_cell_m": R_cell,
                   "warmup": args.warmup, "timed_steps": args.steps,
                   "device": args.device, "gate_steps": args.gate_steps},
        "arms": results,
        "speedups_vs_cupy_cpucomp": {
            k: (v["steps_per_s"] / base)
            for k, v in results.items() if "steps_per_s" in v
        },
    }
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = _OUT_DIR / "h7_compartment_gpu_port.json"
    out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"\n{'='*64}")
    for k, v in results.items():
        if "steps_per_s" in v:
            sp = v["steps_per_s"] / base
            print(f"{k:18s} {v['steps_per_s']:7.1f} steps/s  {sp:5.2f}x  "
                  f"ETA {v['eta_hours_at_gate']:.1f}h  g_soft={v['g_soft_mN_m']:.3e} "
                  f"g_rigid={v['g_rigid_mN_m']:.3e}"
                  + (f"  nonconv={v['nonconverged_count']}" if v.get('nonconverged_count') is not None else ""))
    print(f"→ {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
