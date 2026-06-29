"""LARGE native-mesh + tent DCM spheroid on the GPU (gbook RTX A5000) — run-script.

This is the script the PI runs on the gbook A5000 to exercise the GPU-resident DCM
stack: native md.mesh turgor/edge shell (already GPU-native) + the device-dispatched
``DcmTentContactGPU`` cell-cell contact (cupy K-grid neighbour list, NO host sync) +
the optional ``DcmSubstrateForceGPU``. It builds a spheroid of up to a few hundred
cells, runs it, and measures the spreading observables (A/A0, effective radius,
necrotic fraction) reusing the validated capstone measurement helpers.

DEVICE: uses ``hoomd.device.GPU()`` when available; if HOOMD reports no GPU build
(e.g. this dev Mac — no CUDA), it FALLS BACK to ``hoomd.device.CPU()`` with a clear
log line and the bit-identical CPU dispatch path. The GPU speedup is a gbook-only
result; on CPU this is a structural/finite smoke only.

Run on the gbook A5000 (production):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python -m ffn_sim.scripts.dcm_gpu_run \
        --n-cells 200 --equil-blocks 18 --block 1000

Run on a CPU dev box (structural/finite smoke, small N):
    PYTHONPATH=. ~/miniconda3/envs/ffn_sim/bin/python -m ffn_sim.scripts.dcm_gpu_run \
        --n-cells 12 --equil-blocks 6 --block 500 --force-cpu

Outputs: outputs/h_dcm_gpu/dcm_gpu_run.json  (+ a one-line summary to stdout).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import hoomd  # noqa: E402

from ffn_sim.archive.hoomd_legacy.cell.dcm_native_shell import (  # noqa: E402
    ResolvedNativeDCM,
    build_native_dcm_simulation,
)
from ffn_sim.archive.hoomd_legacy.cell.dcm_gpu_forces import (  # noqa: E402
    DcmTentContactGPU,
    DcmSubstrateForceGPU,
    on_gpu,
)
# Reuse the VALIDATED capstone measurement helpers (does not modify them).
from ffn_sim.scripts.dcm_native_capstone import (  # noqa: E402
    _positions,
    _footprint_area,
    _cell_centroids,
    _effective_radius,
    _necrosis_zones,
)

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "h_dcm_gpu")
os.makedirs(OUT, exist_ok=True)

_UM = 1.0e6


# ---------------------------------------------------------------------------
# Device selection
# ---------------------------------------------------------------------------
def select_device(force_cpu: bool):
    """Return (device, on_gpu_bool, note). Prefer GPU; fall back to CPU + log."""
    if force_cpu:
        return hoomd.device.CPU(notice_level=0), False, "forced CPU (--force-cpu)"
    try:
        if hoomd.device.GPU.is_available():
            return hoomd.device.GPU(notice_level=0), True, "GPU (CUDA available)"
        note = "no CUDA GPU build available — FALLBACK to CPU"
    except Exception as e:  # noqa: BLE001
        note = f"GPU init failed ({e}) — FALLBACK to CPU"
    return hoomd.device.CPU(notice_level=0), False, note


# ---------------------------------------------------------------------------
# Build a large spheroid with the GPU-capable forces wired in
# ---------------------------------------------------------------------------
def build_gpu_spheroid(n_cells: int, *, device, gpu_substrate: bool,
                       dt: float = 3.0e-10, seed: int = 7):
    """Build the native shell, then attach DcmTentContactGPU (+ substrate GPU).

    We let ``build_native_dcm_simulation`` build the native mesh shell + BAOAB on
    the chosen device with ``contact=False, substrate=False``, then append our
    device-dispatched GPU-capable forces. This keeps the GPU contact/substrate as
    NEW wiring without touching the live ``DcmTentContact`` build path.
    """
    p = ResolvedNativeDCM(subdivisions=2, cluster="3d", spacing_factor=2.3,
                          adh_strength=1.0e8, rep_strength=1.0e8,
                          W_cs_Jm2=0.5e-3, dt=dt, seed=seed)
    h = build_native_dcm_simulation(p, n_cells, device=device,
                                    substrate=False, contact=False)
    sim = h["sim"]
    ig = sim.operations.integrator
    nv = h["nv"]
    patch_area = 4.0 * np.pi * p.R_cell ** 2 / nv

    tent = DcmTentContactGPU(
        cell_of_node=h["cell_of_tag"], r_contact=1.05 * h["mean_edge"],
        c_adh=p.c_adh, rep_strength=p.rep_strength,
        adh_strength=p.adh_strength, patch_area=patch_area,
        force_cap=p.force_cap)
    ig.forces.append(tent)

    sub = None
    if gpu_substrate:
        area_per_node = 4.0 * np.pi * p.R_cell ** 2 / nv
        W_cs = p.W_cs_Jm2 * p.ligand_density * area_per_node
        sub = DcmSubstrateForceGPU(z0=p.z_substrate, W_cs=W_cs,
                                   adh_range=p.R_cell, k_sub=p.k_sub_Nm)
        ig.forces.append(sub)

    sim.run(0)
    h["tent"] = tent
    h["substrate"] = sub
    return h, p


# ---------------------------------------------------------------------------
# Run + measure
# ---------------------------------------------------------------------------
def run(args) -> dict:
    device, is_gpu, dev_note = select_device(args.force_cpu)
    print(f"[device] {dev_note}  -> on_gpu={is_gpu}")
    print(f"[build]  n_cells={args.n_cells} subdiv=2 ...")

    t_build0 = time.time()
    h, p = build_gpu_spheroid(args.n_cells, device=device,
                              gpu_substrate=not args.no_gpu_substrate,
                              dt=args.dt, seed=args.seed)
    sim, ranges, vol = h["sim"], h["ranges"], h["volume"]
    R = p.R_cell
    assert on_gpu(sim) == is_gpu
    n_part = sim.state.N_particles
    print(f"[build]  done in {time.time()-t_build0:.1f}s  "
          f"({n_part} particles, {args.n_cells} cells)")

    # --- finite gate: run(0) -> run(2000) ---
    sim.run(0)
    pg0 = _positions(sim)
    A0_self = _footprint_area(pg0, p.z_substrate, R)
    t0 = time.time()
    sim.run(2000)
    pg = _positions(sim)
    finite = bool(np.all(np.isfinite(pg)) and np.all(np.isfinite(vol.volume)))
    print(f"[gate]   run(2000) finite={finite}  "
          f"({2000/(time.time()-t0):.0f} steps/s)")
    if not finite:
        print("[gate]   NON-FINITE -> aborting (drop --dt or check params)")
        return {"finite_gate": False, "n_cells": args.n_cells,
                "on_gpu": is_gpu, "device_note": dev_note}

    # --- equilibrate, track footprint ---
    areas, areas_abs = [], []
    t_run0 = time.time()
    total_steps = 0
    for _ in range(args.equil_blocks):
        sim.run(args.block)
        total_steps += args.block
        pg = _positions(sim)
        if not np.all(np.isfinite(pg)):
            print("[run]    NON-FINITE during equilibration -> abort")
            return {"finite_gate": False, "n_cells": args.n_cells,
                    "on_gpu": is_gpu, "device_note": dev_note}
        A_abs = _footprint_area(pg, p.z_substrate, R)
        areas.append(A_abs / A0_self)
        areas_abs.append(A_abs * _UM ** 2)
    wall = time.time() - t_run0
    sps = total_steps / wall if wall > 0 else float("nan")

    tail = areas[max(1, 2 * len(areas) // 3):]
    AoverA0 = float(np.mean(tail))
    AoverA0_std = float(np.std(tail))

    cents = _cell_centroids(pg, ranges)
    R_max, R_rms, _cc, _r = _effective_radius(cents)
    R_eff_um = (R_max + R) * _UM
    states, depth_um, R_cluster_um, frac = _necrosis_zones(cents, R)

    print(f"[result] R_eff={R_eff_um:.1f} µm  A/A0={AoverA0:.2f}±{AoverA0_std:.2f}  "
          f"necrotic_frac={frac['necrotic']:.2f}  "
          f"({sps:.0f} steps/s over {total_steps} steps)")

    return {
        "finite_gate": True,
        "on_gpu": is_gpu,
        "device_note": dev_note,
        "n_cells": args.n_cells,
        "n_particles": int(n_part),
        "subdivisions": 2,
        "dt": args.dt,
        "seed": args.seed,
        "steps_per_s": sps,
        "total_steps": total_steps,
        "R_eff_um": R_eff_um,
        "R_max_um": float(R_max * _UM),
        "R_rms_um": float(R_rms * _UM),
        "R_cluster_um": R_cluster_um,
        "AoverA0": AoverA0,
        "AoverA0_tailstd": AoverA0_std,
        "A_traj": [float(a) for a in areas],
        "A_abs_um2_traj": [float(a) for a in areas_abs],
        "A0_self_um2": float(A0_self * _UM ** 2),
        "necrotic_frac": frac["necrotic"],
        "quiescent_frac": frac["quiescent"],
        "prolif_frac": frac["proliferating"],
        "gpu_substrate": not args.no_gpu_substrate,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Large native+tent DCM spheroid on GPU (gbook A5000).")
    ap.add_argument("--n-cells", type=int, default=64,
                    help="number of cells in the spheroid (up to a few hundred).")
    ap.add_argument("--equil-blocks", type=int, default=18,
                    help="number of equilibration blocks after the finite gate.")
    ap.add_argument("--block", type=int, default=1000,
                    help="BAOAB steps per equilibration block.")
    ap.add_argument("--dt", type=float, default=3.0e-10, help="BAOAB timestep [s].")
    ap.add_argument("--seed", type=int, default=7, help="RNG seed.")
    ap.add_argument("--force-cpu", action="store_true",
                    help="force the CPU device (dev/parity smoke).")
    ap.add_argument("--no-gpu-substrate", action="store_true",
                    help="omit the GPU substrate force (contact-only).")
    args = ap.parse_args(argv)

    result = run(args)
    jpath = os.path.join(OUT, "dcm_gpu_run.json")
    with open(jpath, "w") as f:
        json.dump(result, f, indent=2)
    print(f"[out]    json -> {os.path.relpath(jpath)}")
    return 0 if result.get("finite_gate") else 1


if __name__ == "__main__":
    raise SystemExit(main())
