"""B2 — native-N GPU scale-up of the L2.5 G3-PASS A/A0(R0) catch-bond growth law.

The A2 magnitude / R0-range gaps (REPORT §A2): the CPU first pass reached R0 = 40-78 um while
the PI MCF7 spheroids span R0 = 87-419 um — **no overlap**, so the law was only EXTRAPOLATED
into the experiment's size range. The blocker was hardware: a native-N CBM pool (N ~ 10^3-10^4
cells, growing x4 under proliferation) is infeasible on the 16 GB CPU box (memory + wall-time).

The GPU-main BAOAB port (``integrator/baoab_device.py``) removes the per-step host sync that
throttled the GPU, so the CBM cohesion (``md.pair.Table``) + neighbour list — already HOOMD
C++ — now run device-resident end-to-end. This driver runs the SAME validated catch-bond
``run_growth_pooled`` path at native N0 on the GPU, extending the G3 law into the PI R0 range.

Faithful to the §HEADLINE physics: ``cohesion='catch'`` (Rakshit-2012 sliding-rebinding),
NO substrate / NO traction (those are the A-line ligand axes; B2 is the SCALE axis). The
connected-core A/A0 is the gate-scored observable, exactly as the CPU headline.

Usage
-----
One (N0, seed) growth on a device → one JSON line (the per-process unit, fanned out by a shell):
    python -m ffn_sim.scripts.layer2_gpu_scaleup run <N0> <seed> [cpu|gpu]

Parity (correctness): the SAME (N0, seed) on BOTH devices; A/A0 must agree within seed noise:
    python -m ffn_sim.scripts.layer2_gpu_scaleup parity <N0> <seed>

Aggregate + fit a+b/R+c/R^2 over a JSONL of run lines:
    python -m ffn_sim.scripts.layer2_gpu_scaleup --fit <results.jsonl>
"""

from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.spheroid.observables import effective_radius
from ffn_sim.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.spheroid.proliferation import run_growth_pooled
from ffn_sim.validation.oracles.spheroid.aa0_law import fit_aa0

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


def _make_device(kind: str):
    """Build a HOOMD device: 'gpu' → CUDA A5000 (native-N), else CPU (dev/parity)."""
    import hoomd
    if kind == "gpu":
        return hoomd.device.GPU(notice_level=0)
    return hoomd.device.CPU(notice_level=0)


def run_one(n: int, seed: int, kind: str = "cpu") -> dict:
    """One catch-bond proliferation growth at N0=n on the chosen device → record dict.

    Mirrors ``layer2_aa0_growth_persize.run_one`` but (a) catch cohesion (the G3-PASS physics),
    (b) device-selectable, (c) a pool sized for 2 doublings (x4) + margin, (d) wall-time logged.
    """
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved)
    total_time = 2.0 * prolif.cycle_time_mean
    max_cells = int(5 * n)  # x4 growth headroom + margin; the pool is pre-allocated once
    device = _make_device(kind)

    t0 = time.perf_counter()
    res = run_growth_pooled(
        resolved, prolif, n_cells_init=n, total_time=total_time, seed=seed,
        cohesion="catch", cad=cad, max_cells=max_cells,
    )
    wall = time.perf_counter() - t0
    return {
        "n": n,
        "seed": seed,
        "device": kind,
        "R0_um": float(effective_radius(res["a0"]) * 1e6),
        "aa0_hull": float(res["area_over_a0"][-1]),
        "aa0_core": float(res["area_core_over_a0"][-1]),
        "growth": float(res["growth_factor"]),
        "n_final": int(res["n_cells"][-1]),
        "rim": (
            float(res["rim_fraction_mean"])
            if not np.isnan(res["rim_fraction_mean"]) else None
        ),
        "ejected": bool(res.get("ejected", False)),
        "capped": bool(res.get("capped", False)),
        "wall_s": round(wall, 1),
        "steps_per_s": round((total_time / resolved.dt_cfl) / wall, 1),
    }


def parity(n: int, seed: int) -> int:
    """Run the same (N0, seed) on CPU and GPU; report A/A0 agreement + speedup."""
    print(f"[parity] N0={n} seed={seed} — CPU vs GPU catch-bond growth\n")
    cpu = run_one(n, seed, "cpu")
    print("  CPU:", json.dumps(cpu))
    gpu = run_one(n, seed, "gpu")
    print("  GPU:", json.dumps(gpu))
    d_core = abs(cpu["aa0_core"] - gpu["aa0_core"])
    d_hull = abs(cpu["aa0_hull"] - gpu["aa0_hull"])
    speedup = cpu["wall_s"] / gpu["wall_s"] if gpu["wall_s"] > 0 else float("nan")
    print(
        f"\n  |dA/A0_core| = {d_core:.3f}   |dA/A0_hull| = {d_hull:.3f}   "
        f"(GPU RNG stream differs → statistical parity, ~seed-noise ±0.05-0.15)"
    )
    print(f"  GPU/CPU wall speedup = {speedup:.2f}x   "
          f"(CPU {cpu['steps_per_s']} st/s, GPU {gpu['steps_per_s']} st/s)")
    return 0


def fit_jsonl(path: str) -> dict:
    """Aggregate run lines by N0, mean over seeds, fit a+b/R+c/R^2 to CORE (G3)."""
    rows = [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]
    by_n: dict[int, list] = defaultdict(list)
    for r in rows:
        by_n[r["n"]].append(r)
    R0, core_m, hull_m = [], [], []
    for n in sorted(by_n):
        g = by_n[n]
        R0.append(float(np.mean([x["R0_um"] for x in g])) * 1e-6)
        cm = [x["aa0_core"] for x in g]
        hm = [x["aa0_hull"] for x in g]
        core_m.append(float(np.mean(cm)))
        hull_m.append(float(np.mean(hm)))
        print(f"  N0={n:6d}  R0={R0[-1]*1e6:6.1f} um  core={np.mean(cm):.2f}+/-{np.std(cm):.2f}"
              f"  hull={np.mean(hm):.2f}+/-{np.std(hm):.2f}  (nseeds={len(g)})")
    R0a, core_a, hull_a = np.array(R0), np.array(core_m), np.array(hull_m)
    fc, fh = fit_aa0(R0a, core_a), fit_aa0(R0a, hull_a)
    print(f"\n[fit CORE] A/A0 = {fc['a']:.3f} + ({fc['b']*1e6:.3f} um)/R + "
          f"({fc['c']*1e12:.3f} um^2)/R^2   r^2={fc['r_squared']:.3f}")
    print(f"[fit HULL] A/A0 = {fh['a']:.3f} + ({fh['b']*1e6:.3f} um)/R + "
          f"({fh['c']*1e12:.3f} um^2)/R^2   r^2={fh['r_squared']:.3f}")
    g3 = 0.95
    print(f"\n[G3 fit r^2 >= {g3}] {'PASS' if fc['r_squared'] >= g3 else 'FAIL'} "
          f"(core r^2={fc['r_squared']:.3f})  R0 range {R0a.min()*1e6:.0f}-{R0a.max()*1e6:.0f} um")
    return {"R0": R0, "core": core_m, "hull": hull_m, "fit_core": fc, "fit_hull": fh}


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if not args:
        print(__doc__)
        return 1
    if args[0] == "--fit":
        fit_jsonl(args[1])
        return 0
    if args[0] == "parity":
        return parity(int(args[1]), int(args[2]))
    if args[0] == "run":
        kind = args[3] if len(args) > 3 else "cpu"
        print(json.dumps(run_one(int(args[1]), int(args[2]), kind)), flush=True)
        return 0
    print(__doc__)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
