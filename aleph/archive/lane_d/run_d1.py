"""D1 driver — run the MT/IF/LINC force-path CUDA trajectory (pulses A & B), dump frames + gate report.

Usage (CUDA host, gbook A5000)::

    python -m aleph.archive.lane_d.run_d1 --device cuda --outdir aleph/outputs/lane_d

Two diagnostic perturbations, each ramped and settled quasi-statically:
    A. membrane-side patch pushed outward → reaction transmitted through capture→anchor→IF/LINC→nucleus.
    B. prescribed nuclear +x displacement pulse (EXTERNAL input, logged) → LINC/IF → anchor → membrane reaction.

Produces (under --outdir):
    d1_trajectory.npz   per-frame global positions + per-channel forces + scalar ledger for A and B
    d1_gates.json       the D1 sanity-gate report
    d1_meta.json        run metadata
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import warp as wp

from .force_path_rig import ForcePathRig, R_MEM_UM
from . import gates_d1


def _require_cuda(device):
    wp.init()
    dev = wp.get_device(device)
    if not dev.is_cuda:
        raise RuntimeError(f"Lane D requires a CUDA device; got {dev}. Run on the A5000 (gbook).")
    return str(dev)


def _channel_forces(rig):
    return {c: rig.channel_force(c) for c in ("linc", "if", "mt", "capture", "nucleus")}


def run(device, outdir, n_ramp=24, push_A=0.8, pull_B=0.6, settle_iter=60) -> dict:
    dev = _require_cuda(device)
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()

    rig = ForcePathRig(device=dev)
    resid = rig.relax_to_rest(n_iter=300)
    print(f"D1 rest residual: {resid:.4f} pN | N={rig.N} (mem {rig.Nm}, nuc {rig.Nn}, anc {rig.Na}, "
          f"MT {rig.Nmt}, IF {rig.Nif}), capture {rig.n_capture}, LINC {rig.n_linc}")
    rest_pos = rig.pos_d.numpy().copy()
    nuc = rig.nucleus_nodes()
    patch = rig.off_mem + rig.membrane_patch_nodes()

    # ---------------- Pulse A: membrane push → nucleus reaction ----------------
    rig.set_fixed(nuc, True)
    framesA, ledgerA = [], {k: [] for k in ("step", "applied_um", "nuc_reaction_pn",
                                            "linc_pn", "if_pn", "mt_pn", "load_at_patch_pn")}
    for i in range(n_ramp):
        frac = (i + 1) / n_ramp
        rig.pos_d.assign(np.ascontiguousarray(rest_pos))  # reset then displace patch to fractional push
        p = rest_pos.copy(); p[patch, 0] += frac * push_A
        rig.pos_d.assign(np.ascontiguousarray(p)); rig.set_fixed(patch, True)
        rig.settle(settle_iter)
        rig.accumulate_all(); f = rig.f_d.numpy()
        ch = _channel_forces(rig)
        framesA.append(rig.pos_d.numpy().astype(np.float32))
        ledgerA["step"].append(i)
        ledgerA["applied_um"].append(frac * push_A)
        ledgerA["nuc_reaction_pn"].append(float(np.linalg.norm(np.sum(f[nuc], axis=0))))
        ledgerA["linc_pn"].append(float(np.max(np.linalg.norm(ch["linc"], axis=1))))
        ledgerA["if_pn"].append(float(np.max(np.linalg.norm(ch["if"], axis=1))))
        ledgerA["mt_pn"].append(float(np.max(np.linalg.norm(ch["mt"], axis=1))))
        ledgerA["load_at_patch_pn"].append(float(np.linalg.norm(np.sum(f[patch], axis=0))))
    chA = _channel_forces(rig)  # final-frame channel forces for the force-path figure

    # ---------------- Pulse B: prescribed nuclear displacement → membrane-anchor reaction ----------------
    rig.restore({"pos": rest_pos, "linc_n": rig.linc_n_d.numpy(), "cap_bound": rig.cap_bound_d.numpy(),
                 "linc_bound": rig.linc_bound, "epoch": rig.topology_epoch, "fixed": np.zeros(rig.N, np.int32)})
    rig.set_fixed(patch, True)  # membrane frame (incl. patch) pinned; nucleus prescribed
    anchors = np.arange(rig.Na) + rig.off_anc
    framesB, ledgerB = [], {k: [] for k in ("step", "prescribed_um", "anchor_reaction_pn",
                                            "linc_pn", "if_pn", "membrane_reaction_pn")}
    for i in range(n_ramp):
        frac = (i + 1) / n_ramp
        p = rest_pos.copy(); p[nuc, 0] += frac * pull_B    # EXTERNAL prescribed nuclear +x displacement pulse
        rig.pos_d.assign(np.ascontiguousarray(p)); rig.set_fixed(nuc, True)
        rig.settle(settle_iter)
        rig.accumulate_all(); f = rig.f_d.numpy()
        ch = _channel_forces(rig)
        framesB.append(rig.pos_d.numpy().astype(np.float32))
        ledgerB["step"].append(i)
        ledgerB["prescribed_um"].append(frac * pull_B)
        ledgerB["anchor_reaction_pn"].append(float(np.linalg.norm(np.sum(f[anchors], axis=0))))
        ledgerB["linc_pn"].append(float(np.max(np.linalg.norm(ch["linc"], axis=1))))
        ledgerB["if_pn"].append(float(np.max(np.linalg.norm(ch["if"], axis=1))))
        ledgerB["membrane_reaction_pn"].append(float(np.linalg.norm(np.sum(f[patch], axis=0))))
    chB = _channel_forces(rig)

    npz = os.path.join(outdir, "d1_trajectory.npz")
    np.savez_compressed(
        npz,
        rest_pos=rest_pos.astype(np.float32), framesA=np.stack(framesA), framesB=np.stack(framesB),
        nuc_slice=np.array([rig.off_nuc, rig.off_nuc + rig.Nn]),
        anc_slice=np.array([rig.off_anc, rig.off_anc + rig.Na]),
        mem_n=rig.Nm, off_mt=rig.off_mt, n_mt=rig.n_mt,
        mt_fiber_off=np.asarray(rig.mt.fiber_off), patch=patch.astype(np.int64),
        linc_n=rig._linc_n, linc_a=rig._linc_a, cap_ia=rig.cap_ia_d.numpy(), cap_ib=rig.cap_ib_d.numpy(),
        nuc_faces=rig.nuc_faces, mem_faces=rig.mem_faces, r_mem=R_MEM_UM,
        chA_linc=chA["linc"], chA_if=chA["if"], chA_mt=chA["mt"], chA_capture=chA["capture"],
        chB_linc=chB["linc"], chB_if=chB["if"],
        **{f"A_{k}": np.asarray(v) for k, v in ledgerA.items()},
        **{f"B_{k}": np.asarray(v) for k, v in ledgerB.items()},
    )

    gate_report = gates_d1.run_all(dev)
    with open(os.path.join(outdir, "d1_gates.json"), "w") as fh:
        json.dump(gate_report, fh, indent=2)

    meta = {
        "device": dev, "N": rig.N, "Nm": rig.Nm, "Nn": rig.Nn, "Na": rig.Na, "Nmt": rig.Nmt, "Nif": rig.Nif,
        "n_capture": rig.n_capture, "n_linc": rig.n_linc, "rest_residual_pn": resid,
        "push_A_um": push_A, "pull_B_um": pull_B, "n_ramp": n_ramp,
        "final_nuc_reaction_A_pn": ledgerA["nuc_reaction_pn"][-1],
        "final_membrane_reaction_B_pn": ledgerB["membrane_reaction_pn"][-1],
        "wall_time_s": time.time() - t0, "trajectory": npz,
        "all_gates_passed": gate_report["all_passed"],
    }
    with open(os.path.join(outdir, "d1_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(json.dumps({"D1": meta}, indent=2))
    print("GATES:", json.dumps({g["name"]: g["passed"] for g in gate_report["gates"]}, indent=2))
    return meta


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default=None)
    ap.add_argument("--outdir", default="aleph/outputs/lane_d")
    ap.add_argument("--n-ramp", type=int, default=24)
    ap.add_argument("--push-A", type=float, default=0.8)
    ap.add_argument("--pull-B", type=float, default=0.6)
    ap.add_argument("--settle-iter", type=int, default=60)
    a = ap.parse_args()
    run(a.device, a.outdir, n_ramp=a.n_ramp, push_A=a.push_A, pull_B=a.pull_B, settle_iter=a.settle_iter)
