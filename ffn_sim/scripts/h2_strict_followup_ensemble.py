"""H.2 strict-PASS follow-up — 5-seed ensemble runner.

Runs the canonical H.2 production simulation across 5 seeds in parallel
via parallel_ensemble.run_seed_pool, with the current yaml (which the
strict-PASS session sets to slab Lz = 10 μm for the slab-confinement
diagnostic).  Each seed produces an independent trajectory + summary;
results are aggregated into a single JSON for downstream metric
comparison vs first-principles.

PI 2026-05-25: strict-PASS follow-up diagnostic.  Single-seed canonical
already showed +50% equipartition deviation; ensemble quantifies the
seed-to-seed scatter on that deviation and tightens the slab-diagnostic
conclusion.
"""

from __future__ import annotations

import json
import math
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
import yaml


CONFIG_PATH = Path(__file__).resolve().parents[1] / "configs" / "phase1_h2.yaml"
OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h2" / "strict_followup"


def _seed_task(seed: int, run_label: str) -> dict:
    """One worker: build H.2 sim at this seed, sample, summarise."""
    from ffn_sim.common.filament_math import (
        bending_energy_per_bond,
        equipartition_check,
        fit_persistence_length,
        tangent_correlation,
    )
    from ffn_sim.scripts.h2_single_filament import (
        resolve_h2_derived,
        run_h2_and_sample,
    )

    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    cfg = deepcopy(cfg)
    cfg["filament"]["seed"] = int(seed)
    p = resolve_h2_derived(cfg)

    t0 = time.time()
    result = run_h2_and_sample(p)
    wall_s = time.time() - t0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(
        OUT_DIR / f"{run_label}_seed{seed:02d}.npz",
        positions=result["positions"],
        frame_steps=result["frame_steps"],
        wall_s=wall_s,
        seed=seed,
    )

    pos_all = result["positions"].reshape(-1, p.beads_per_fiber, 3)
    C = tangent_correlation(pos_all, box=None, max_separation=10)
    L_p_C1 = float(-p.rest_length / math.log(C[1])) if 0 < C[1] < 1 else float("nan")
    tail = fit_persistence_length(pos_all, rest_length=p.rest_length, box=None)
    energies = [
        bending_energy_per_bond(frame, angle_k=p.angle_k, angle_t0=p.angle_t0, box=None)
        for frame in result["positions"]
    ]
    eq = equipartition_check(energies, kT_J=p.kT)
    target_J = p.equipartition_target_3d_kT * p.kT

    return {
        "seed": int(seed),
        "wall_s": wall_s,
        "C_s": [float(c) for c in C.tolist()],
        "L_p_C1_m": L_p_C1,
        "L_p_tail_m": float(tail.L_p_m),
        "E_mean_kT": float(eq.mean_J / p.kT),
        "rel_vs_3d_kT": float((eq.mean_J - target_J) / target_J),
        "box_Lz_m": float(p.box_Lz),
    }


def main(run_label: str = "slab_Lz10um", seeds: list[int] = None) -> dict:
    from ffn_sim.scripts.parallel_ensemble import run_seed_pool

    if seeds is None:
        seeds = [42, 43, 44, 45, 46]
    print(f"[ensemble] launching {len(seeds)} seeds: {seeds}", flush=True)
    print(f"[ensemble] yaml at {CONFIG_PATH}", flush=True)

    results = run_seed_pool(
        _seed_task, seeds=seeds, n_workers=min(len(seeds), 6),
        extra_kwargs={"run_label": run_label},
    )

    aggregate = {
        "run_label": run_label,
        "seeds": seeds,
        "per_seed": results,
        "summary": {
            "L_p_C1_mean_m":   float(np.mean([r["L_p_C1_m"]  for r in results])),
            "L_p_C1_std_m":    float(np.std (([r["L_p_C1_m"] for r in results]), ddof=1)),
            "L_p_tail_mean_m": float(np.mean([r["L_p_tail_m"] for r in results])),
            "L_p_tail_std_m":  float(np.std (([r["L_p_tail_m"] for r in results]), ddof=1)),
            "E_mean_kT_mean":  float(np.mean([r["E_mean_kT"]  for r in results])),
            "E_mean_kT_std":   float(np.std (([r["E_mean_kT"] for r in results]), ddof=1)),
            "rel_vs_3d_kT_mean": float(np.mean([r["rel_vs_3d_kT"] for r in results])),
            "rel_vs_3d_kT_std":  float(np.std (([r["rel_vs_3d_kT"] for r in results]), ddof=1)),
            "wall_max_s": float(max(r["wall_s"] for r in results)),
        },
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_json = OUT_DIR / f"{run_label}_aggregate.json"
    with open(out_json, "w") as f:
        json.dump(aggregate, f, indent=2)
    print(f"[ensemble] aggregate written to {out_json}", flush=True)

    s = aggregate["summary"]
    print(
        f"[ensemble] summary ({run_label}, n={len(seeds)} seeds):\n"
        f"  L_p_C1   = {s['L_p_C1_mean_m']*1e6:.3f} ± {s['L_p_C1_std_m']*1e6:.3f} μm   "
        f"(first-principles 16.44 ± 0.13 μm)\n"
        f"  L_p_tail = {s['L_p_tail_mean_m']*1e6:.3f} ± {s['L_p_tail_std_m']*1e6:.3f} μm   "
        f"(first-principles 16.44 μm)\n"
        f"  ⟨E⟩      = {s['E_mean_kT_mean']:.4f} ± {s['E_mean_kT_std']:.4f} kT   "
        f"(first-principles 0.99 kT)\n"
        f"  rel_3d   = {s['rel_vs_3d_kT_mean']:+.3f} ± {s['rel_vs_3d_kT_std']:.3f}",
        flush=True,
    )
    return aggregate


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="slab_Lz10um")
    ap.add_argument("--seeds", type=int, nargs="+", default=[42, 43, 44, 45, 46])
    args = ap.parse_args()
    main(run_label=args.label, seeds=args.seeds)
