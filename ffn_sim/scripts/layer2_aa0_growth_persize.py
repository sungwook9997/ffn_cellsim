"""L2.4 memory-safe per-size driver for the proliferation-driven A/A0(R0) growth sweep.

Why this exists alongside ``layer2_aa0_growth_sweep`` (the all-in-one sweep): the all-in-one
sweep grows all 5 sizes x 3 seeds in ONE process, and the HOOMD epoch-rebuild loop accumulates
enough resident memory across sizes that the OS kills it around N0~250-400 on a 16 GB CPU box
(SIGKILL/137). This driver runs ONE (n_cells_init, seed) per process so memory is fully
reclaimed between sizes; a thin shell/Python loop fans out the 15 runs and aggregates the
JSONL. It uses the SAME tested ``run_growth`` + ``core_projected_area`` path, so the numbers
are canonical — it is a hardware workaround, not a different measurement.

Run one cell:
    PYTHONPATH=<repo> python -m ffn_sim.scripts.layer2_aa0_growth_persize <N0> <seed>
        -> one JSON line: {n, seed, R0, aa0_hull, aa0_core, growth, rim}

Aggregate + fit (a + b/R + c/R^2) over a JSONL of those lines:
    python -m ffn_sim.scripts.layer2_aa0_growth_persize --fit results.jsonl
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.observables import effective_radius
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.archive.hoomd_legacy.spheroid.proliferation import run_growth_pooled
from ffn_sim.validation.oracles.spheroid.aa0_law import fit_aa0

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


def run_one(n: int, seed: int) -> dict:
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    total_time = 2.0 * prolif.cycle_time_mean
    # leak-free pooled growth (single Simulation + pre-allocated pool) — flat ~240 MB/run,
    # so per-process isolation is now belt-and-suspenders rather than a necessity.
    res = run_growth_pooled(resolved, prolif, n_cells_init=n, total_time=total_time, seed=seed)
    return {
        "n": n,
        "seed": seed,
        "R0": float(effective_radius(res["a0"])),
        "aa0_hull": float(res["area_over_a0"][-1]),
        "aa0_core": float(res["area_core_over_a0"][-1]),
        "growth": float(res["growth_factor"]),
        "rim": (
            float(res["rim_fraction_mean"])
            if not np.isnan(res["rim_fraction_mean"])
            else None
        ),
    }


def fit_jsonl(path: str) -> dict:
    rows = [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
    by_n: dict[int, list] = defaultdict(list)
    for r in rows:
        by_n[r["n"]].append(r)
    R0, core_m, hull_m = [], [], []
    for n in sorted(by_n):
        g = by_n[n]
        R0.append(float(np.mean([x["R0"] for x in g])))
        cm = [x["aa0_core"] for x in g]
        hm = [x["aa0_hull"] for x in g]
        core_m.append(float(np.mean(cm)))
        hull_m.append(float(np.mean(hm)))
        print(
            f"  N0={n:4d}  R0={R0[-1]*1e6:5.1f} um  core={np.mean(cm):.2f}+/-{np.std(cm):.2f}"
            f"  hull={np.mean(hm):.2f}+/-{np.std(hm):.2f}  (nseeds={len(g)})"
        )
    R0a, core_a, hull_a = np.array(R0), np.array(core_m), np.array(hull_m)
    fc, fh = fit_aa0(R0a, core_a), fit_aa0(R0a, hull_a)
    print(
        f"\n[fit CORE] A/A0 = {fc['a']:.3f} + ({fc['b']*1e6:.3f} um)/R + "
        f"({fc['c']*1e12:.3f} um^2)/R^2   r^2={fc['r_squared']:.3f}"
    )
    print(
        f"[fit HULL] A/A0 = {fh['a']:.3f} + ({fh['b']*1e6:.3f} um)/R + "
        f"({fh['c']*1e12:.3f} um^2)/R^2   r^2={fh['r_squared']:.3f}"
    )
    g3_min = 0.95
    print(
        f"\n[G3 fit r^2 >= {g3_min}] {'PASS' if fc['r_squared'] >= g3_min else 'FAIL'} "
        f"(core r^2={fc['r_squared']:.3f}) -- gate scored on the fragmentation-robust CORE area"
    )
    return {"R0": R0, "core": core_m, "hull": hull_m, "fit_core": fc, "fit_hull": fh}


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if args and args[0] == "--fit":
        fit_jsonl(args[1])
        return 0
    n, seed = int(args[0]), int(args[1])
    print(json.dumps(run_one(n, seed)), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
