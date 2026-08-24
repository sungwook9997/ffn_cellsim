"""D0 driver — run the membrane+cytosol CUDA trajectory, dump frames + gate report.

Usage (on a CUDA host, e.g. gbook A5000)::

    python -m aleph.archive.lane_d.run_d0 --device cuda \
        --outdir aleph/outputs/lane_d --steps 60

Produces (under ``--outdir``):
    d0_trajectory.npz   per-frame membrane pos/disp/traction + midplane pressure + scalar ledger
    d0_gates.json       the D0 sanity-gate report
    d0_meta.json        run metadata (device, sizes, params, timing)
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import warp as wp

from .membrane_cytosol_rig import (
    MembraneCytosolRig, DP_HYD_REST_PA, DPI_OSM_REST_PA, DX_UM,
)
from . import gates_d0


def _require_cuda(device):
    wp.init()
    dev = wp.get_device(device)
    if not dev.is_cuda:
        raise RuntimeError(f"Lane D requires a CUDA device; got {dev} (type={dev.arch}). "
                           "Run on the A5000 (gbook) — there is no CPU simulation path.")
    return str(dev)


def run(device: str, outdir: str, steps: int = 60, dt_phys: float = 0.02,
        inv_gamma: float = 0.05, pulse_steps: int = 25, pulse_amp: float = 120.0,
        pulse_radius: float = 2.0) -> dict:
    dev = _require_cuda(device)
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()

    rig = MembraneCytosolRig(device=dev)
    resid = rig.relax_to_rest()
    print(f"pre-relax residual max force: {resid:.4f} pN")
    patch = rig.patch_centre()

    frames_pos, frames_disp, frames_trac, frames_pmid = [], [], [], []
    ledger = {k: [] for k in ("t", "w_membrane", "w_fluid", "flux", "content",
                              "p_excess_max", "disp_max", "bulge_x", "trac_plus", "trac_minus")}
    p_bar_h = rig.grid.p_bar.numpy()
    c = rig.verts0.mean(axis=0)
    cap_p = (rig.verts0[:, 0] - c[0]) > 0.5 * rig.r_mem_um
    cap_m = (rig.verts0[:, 0] - c[0]) < -0.5 * rig.r_mem_um

    # record rest frame — pmid stores the pressure EXCESS (p − p_bar) so the perturbation pops
    def record(t):
        st = rig.membrane_state()
        frames_pos.append(st["pos"].astype(np.float32))
        frames_disp.append(st["disp"].astype(np.float32))
        frames_trac.append(st["traction"].astype(np.float32))
        pex_mid = (rig.grid.pressure_to_host() - p_bar_h)[:, :, rig._nx // 2]
        frames_pmid.append(pex_mid.astype(np.float32))
        pex = rig.grid.pressure_to_host() - p_bar_h
        outward = np.sum(st["traction"] * st["normals0"], axis=1)   # outward traction per node [pN]
        radial = np.sum(st["disp"] * st["normals0"], axis=1)        # radial displacement from rest [µm]
        ledger["t"].append(float(t))
        ledger["p_excess_max"].append(float(np.max(pex)))
        ledger["disp_max"].append(float(np.max(np.linalg.norm(st["disp"], axis=1))))
        ledger["bulge_x"].append(float(np.mean(radial[cap_p])))
        ledger["trac_plus"].append(float(np.mean(outward[cap_p])))
        ledger["trac_minus"].append(float(np.mean(outward[cap_m])))

    record(0.0)
    # perturbation: a sustained localized cytosol overpressure (+pulse_amp Pa, a fraction of the 40 Pa
    # turgor) HELD in the +x patch for `pulse_steps` outer steps so the stiff membrane can bulge against a
    # real load, then RELEASED — the excess then spreads through the cytosol and the traction redistributes.
    t = 0.0
    for s in range(steps):
        hold = (patch, pulse_radius, pulse_amp) if s < pulse_steps else None
        rec = rig.outer_step(dt_phys, inv_gamma=inv_gamma, dpi_osm=DPI_OSM_REST_PA,
                             move_membrane=True, hold_patch=hold)
        t += dt_phys
        ledger["w_membrane"].append(rec["w_membrane_pn_um_s"])
        ledger["w_fluid"].append(rec["w_fluid_pn_um_s"])
        ledger["flux"].append(rec["integrated_flux_um3_s"])
        ledger["content"].append(rec["total_content"])
        record(t)

    npz_path = os.path.join(outdir, "d0_trajectory.npz")
    np.savez_compressed(
        npz_path,
        pos=np.stack(frames_pos), disp=np.stack(frames_disp), traction=np.stack(frames_trac),
        pmid=np.stack(frames_pmid), faces=rig.faces_np, verts0=rig.verts0.astype(np.float32),
        origin=rig._origin, dx=DX_UM, r_mem=rig.r_mem_um, patch=np.asarray(patch),
        pulse_steps=pulse_steps, dt_phys=dt_phys,
        t=np.asarray(ledger["t"]),
        w_membrane=np.asarray(ledger["w_membrane"]), w_fluid=np.asarray(ledger["w_fluid"]),
        flux=np.asarray(ledger["flux"]), content=np.asarray(ledger["content"]),
        p_excess_max=np.asarray(ledger["p_excess_max"]), disp_max=np.asarray(ledger["disp_max"]),
        bulge_x=np.asarray(ledger["bulge_x"]),
        trac_plus=np.asarray(ledger["trac_plus"]), trac_minus=np.asarray(ledger["trac_minus"]),
    )

    # gates on fresh rigs
    gate_report = gates_d0.run_all(dev)
    with open(os.path.join(outdir, "d0_gates.json"), "w") as fh:
        json.dump(gate_report, fh, indent=2)

    meta = {
        "device": dev, "n_nodes": rig.n_nodes, "n_faces": rig.n_faces,
        "grid_nx": rig._nx, "n_fluid_cells": rig.n_fluid, "steps": steps, "dt_phys": dt_phys,
        "dP_hyd_rest_pa": DP_HYD_REST_PA, "dPi_osm_rest_pa": DPI_OSM_REST_PA,
        "kappa_tilde": rig.kappa_tilde, "wall_time_s": time.time() - t0,
        "peak_gpu_bytes": _peak_bytes(dev), "trajectory": npz_path,
        "all_gates_passed": gate_report["all_passed"],
    }
    with open(os.path.join(outdir, "d0_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(json.dumps({"D0": meta}, indent=2))
    print("GATES:", json.dumps({g["name"]: g["passed"] for g in gate_report["gates"]}, indent=2))
    return meta


def _peak_bytes(dev):
    try:
        return int(wp.get_device(dev).total_memory) - int(wp.get_device(dev).free_memory)
    except Exception:
        return -1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default=None)
    ap.add_argument("--outdir", default="aleph/outputs/lane_d")
    ap.add_argument("--steps", type=int, default=60)
    ap.add_argument("--dt-phys", type=float, default=0.02)
    ap.add_argument("--inv-gamma", type=float, default=0.05)
    ap.add_argument("--pulse-steps", type=int, default=25)
    ap.add_argument("--pulse-amp", type=float, default=120.0)
    ap.add_argument("--pulse-radius", type=float, default=2.0)
    a = ap.parse_args()
    run(a.device, a.outdir, steps=a.steps, dt_phys=a.dt_phys, inv_gamma=a.inv_gamma,
        pulse_steps=a.pulse_steps, pulse_amp=a.pulse_amp, pulse_radius=a.pulse_radius)
