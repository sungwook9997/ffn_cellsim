"""Evaluate a DCM spread run for GENUINE collective spreading (vs peeling vs no-spread).

Reports, from a --save-frames npz:
  • A/A0 rise + maxZ DROP (the spheroid→pancake flattening signal; maxZ holding = no flatten)
  • per-cell rim-vs-bulk horizontal displacement (peeling = rim flies, bulk frozen; cohesive =
    both move together) + a peeling index = rim_max / bulk_mean
  • cell sphericity Ψ (1=sphere, <1=flattened tissue) initial vs final
  • V/V0, interpenetration pen

Run:  PYTHONPATH=. python -m ffn_sim.scripts.spread_eval --npz <file> [--json <file>]
"""

from __future__ import annotations

import argparse
import json as _json
import numpy as np
from scipy.spatial import cKDTree

R = 7.5e-6
UM = 1e6


def _cell_psi(P, faces, fcell, c):
    f = faces[fcell == c]
    v = P[f]
    a, b, cc = v[:, 0], v[:, 1], v[:, 2]
    V = abs(np.einsum("ij,ij->i", np.cross(a, b), cc).sum() / 6.0)
    A = 0.5 * np.linalg.norm(np.cross(b - a, cc - a), axis=1).sum()
    return np.pi ** (1 / 3) * (6 * V) ** (2 / 3) / A


def evaluate(npz_path: str) -> dict:
    d = np.load(npz_path, allow_pickle=True)
    F = d["frames"].astype(np.float64)
    cof = d["cof"]; faces = d["faces"]; fcell = cof[faces[:, 0]]
    cells = np.unique(cof[cof >= 0])
    nc = cells.size
    P0, PE = F[0], F[-1]
    z0 = P0[cof >= 0][:, 2].min()
    c0 = np.array([P0[cof == c].mean(0) for c in cells])
    cE = np.array([PE[cof == c].mean(0) for c in cells])
    h0 = (c0[:, 2] - z0) / R
    rim = h0 < 1.5
    dxy = np.linalg.norm(cE[:, :2] - c0[:, :2], axis=1) * UM
    rim_max = dxy[rim].max() if rim.any() else 0.0
    bulk_mean = dxy[~rim].mean() if (~rim).any() else 0.0
    peel_idx = rim_max / max(bulk_mean, 1e-6)
    psi0 = np.mean([_cell_psi(P0, faces, fcell, c) for c in cells])
    psiE = np.mean([_cell_psi(PE, faces, fcell, c) for c in cells])
    aa0 = d["aa0"]; maxZ = d["maxZ"]; vv0 = d["vv0"]; pen = d["pen_frac"]
    return {
        "n_cells": int(nc), "aa0_final": float(aa0[-1]), "aa0_peak": float(aa0.max()),
        "maxZ_init": float(maxZ[0]), "maxZ_final": float(maxZ[-1]),
        "maxZ_drop_pct": float(100 * (maxZ[0] - maxZ[-1]) / maxZ[0]),
        "rim_disp_max_um": float(rim_max), "bulk_disp_mean_um": float(bulk_mean),
        "peeling_index": float(peel_idx), "psi_init": float(psi0), "psi_final": float(psiE),
        "vv0_final": float(vv0[-1]), "pen_final": float(pen[-1]), "pen_peak": float(pen.max()),
    }


def verdict(m: dict) -> str:
    # GENUINE spreading: A/A0 up, maxZ drops, low peeling, cells flatten, volume held
    spreads = m["aa0_final"] > 1.10
    flattens = m["maxZ_drop_pct"] > 3.0 or m["psi_final"] < m["psi_init"] - 0.01
    cohesive = m["peeling_index"] < 8.0
    stable = abs(m["vv0_final"] - 1.0) < 0.05
    if not stable:
        return "UNSTABLE"
    if m["peeling_index"] > 20:
        return "PEELING (rim tears, bulk frozen)"
    if spreads and flattens and cohesive:
        return "GENUINE COLLECTIVE SPREAD"
    if m["aa0_final"] > 1.03 and cohesive:
        return "modest cohesive spread"
    return "no/low spread"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--label", default="")
    args = ap.parse_args()
    m = evaluate(args.npz)
    lab = args.label or args.npz.split("/")[-1]
    print(f"=== {lab} (N={m['n_cells']}) ===")
    print(f"  A/A0 final={m['aa0_final']:.3f} peak={m['aa0_peak']:.3f}")
    print(f"  maxZ {m['maxZ_init']:.1f}→{m['maxZ_final']:.1f}µm (drop {m['maxZ_drop_pct']:.1f}% — flattening signal)")
    print(f"  per-cell: rim_max={m['rim_disp_max_um']:.1f}µm bulk_mean={m['bulk_disp_mean_um']:.2f}µm  peeling_index={m['peeling_index']:.1f}")
    print(f"  sphericity Ψ {m['psi_init']:.3f}→{m['psi_final']:.3f}  V/V0={m['vv0_final']:.3f}  pen={m['pen_final']:.2f}(peak {m['pen_peak']:.2f})")
    print(f"  VERDICT: {verdict(m)}")


if __name__ == "__main__":
    main()
