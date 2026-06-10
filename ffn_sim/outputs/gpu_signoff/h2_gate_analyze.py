"""H.2 single-filament L_p gate metrics from a trajectory npz (device-agnostic).

Reuses the SAME estimators the pytest gate uses (test_persistence_length.py
TestH2Production): C(1)-local L_p and fit-tail L_p via filament_math, plus 3D
equipartition. "driver, not gate" — informational for the GPU sign-off; the
authoritative ratification stays with the pytest gate. Usage:

    PYTHONPATH=. python h2_gate_analyze.py <traj_gpu.npz> <traj_cpu.npz>
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

from ffn_sim.common.filament_math import (
    tangent_correlation,
    fit_persistence_length,
    equipartition_check,
)

# Gate bands (from configs/phase1_h2.yaml, mirrored in the pytest gate).
BAND_C1_UM = (14.8, 18.0)
BAND_TAIL_UM = (8.0, 33.0)
EQUIPART_TARGET_KT = 0.9898   # 3D analytical mean
EQUIPART_TOL = 0.60           # test_equipartition_within_tol tolerance


def analyze(npz_path: str) -> dict:
    d = np.load(npz_path)
    pos = d["positions"]                       # (n_frames, 1, N, 3)
    rest_length = float(d["rest_length"])
    angle_k = float(d["angle_k"])
    kT = float(d["kT"])
    pos_all = pos.reshape(-1, pos.shape[-2], 3)  # (n_frames, N, 3)

    # C(1)-local L_p (exact gate formula).
    C = tangent_correlation(pos_all, box=None, max_separation=1)
    L_p_C1 = float(-rest_length / math.log(C[1])) if 0.0 < C[1] < 1.0 else float("nan")
    # Fit-tail L_p.
    fit = fit_persistence_length(pos_all, rest_length=rest_length, box=None)
    L_p_tail = float(fit.L_p_m)

    # 3D equipartition: per-interior-bond bending energy ½ k_angle θ² → compare via
    # equipartition_check (mean vs kT/2; report mean/kT for the 3D target).
    # bending energy per interior bead from the bond angle.
    e_samples = []
    for f in range(pos_all.shape[0]):
        p = pos_all[f]
        v1 = p[1:-1] - p[:-2]
        v2 = p[2:] - p[1:-1]
        n1 = np.linalg.norm(v1, axis=1); n2 = np.linalg.norm(v2, axis=1)
        costh = np.clip((v1 * v2).sum(1) / (n1 * n2).clip(min=1e-30), -1.0, 1.0)
        theta = np.arccos(costh)
        e_samples.append(0.5 * angle_k * theta ** 2)         # (N-2,)
    eq = equipartition_check(e_samples, kT_J=kT)
    mean_over_kT = float(eq.mean_J / kT)

    return {
        "npz": Path(npz_path).name,
        "n_frames": int(pos_all.shape[0]),
        "L_p_C1_um": L_p_C1 * 1e6,
        "L_p_tail_um": L_p_tail * 1e6,
        "equipartition_mean_over_kT": mean_over_kT,
        "C1": float(C[1]),
        "in_band_C1": bool(BAND_C1_UM[0] <= L_p_C1 * 1e6 <= BAND_C1_UM[1]),
        "in_band_tail": bool(BAND_TAIL_UM[0] <= L_p_tail * 1e6 <= BAND_TAIL_UM[1]),
        "equipart_within_tol": bool(abs(mean_over_kT - EQUIPART_TARGET_KT)
                                    <= EQUIPART_TOL * EQUIPART_TARGET_KT),
    }


def main() -> int:
    paths = sys.argv[1:]
    if not paths:
        print("usage: h2_gate_analyze.py <npz> [<npz> ...]")
        return 1
    out = {}
    for p in paths:
        dev = "gpu" if "gpu" in Path(p).name.lower() else (
            "cpu" if "cpu" in Path(p).name.lower() else Path(p).stem)
        r = analyze(p)
        out[dev] = r
        print(f"[{dev}] frames={r['n_frames']}  L_p_C1={r['L_p_C1_um']:.2f}um "
              f"(band {BAND_C1_UM}, in={r['in_band_C1']})  "
              f"L_p_tail={r['L_p_tail_um']:.2f}um (in={r['in_band_tail']})  "
              f"<E>/kT={r['equipartition_mean_over_kT']:.3f} "
              f"(tol@{EQUIPART_TARGET_KT}±{EQUIPART_TOL*100:.0f}%, in={r['equipart_within_tol']})")
    if "gpu" in out and "cpu" in out:
        dC1 = out["gpu"]["L_p_C1_um"] - out["cpu"]["L_p_C1_um"]
        print(f"\nGPU vs CPU: ΔL_p_C1 = {dC1:+.2f} um "
              f"({abs(dC1)/max(out['cpu']['L_p_C1_um'],1e-9)*100:.1f}% — cupy RNG≠numpy → statistical)")
    Path(__file__).with_name("h2_gate_result.json").write_text(json.dumps(out, indent=2))
    print(f"json -> {Path(__file__).with_name('h2_gate_result.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
