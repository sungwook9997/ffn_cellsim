"""DECISIVE experiment (PI 2026-06-04): does the lamellipodium-ANCHORED active traction close
the ~5-9x magnitude gap, or is the residual the CBM 1-particle structural limit?

Bracket the two PI-ratified anchors (resolve_active_traction, ligand_traction.py) on the
substrate-confined catch-cohesion CBM at a matched R0:
  baseline  f_active = 0      substrate ON   (substrate-only spreading, no active traction)
  arm A     f_active = 1.6 nN substrate ON   normal dt  (whole-cell MCF7 overlay anchor; in B1 band)
  arm B     f_active = 9.4 nN substrate ON   dt/5       (protrusion Bieling anchor; > B1 ceiling,
            PI-authorized smaller-dt / sub-step for THIS Layer-2 arm only — integrator-freeze
            stays intact for the single-cell line; the smaller dt resolves boundary-cell peeling
            finer so a marginally-detaching cell can be re-captured by the tension-strengthening
            catch bond instead of being flung F/gamma per coarse epoch step).

dt-substep is implemented WITHOUT touching the frozen integrator: a copy of the resolved params
with dt_cfl/factor + epoch_steps*factor keeps the biological time identical, only finer.

READ THE FORK: if arm B's A/A0 jumps toward the PI medians (7-10) -> the gap was MISSING ACTIVE
MAGNITUDE (the heuristic <=3 nN cap suppressed it). If arm B ejects boundary cells or stays
~1-2 even held -> the gap is the 1-particle CBM ABSTRACTION (cannot crawl-while-attached /
resolve contact-line traction) -> a Layer-2-vs-fine-grained scope finding, not a parameter.

Usage:  python -m ffn_sim.scripts.layer2_decisive_traction <N0> <seed> [cpu|gpu]
        -> one JSON line per arm (baseline/A/B) with core+raw A/A0, ejected, wall.
"""

from __future__ import annotations

import dataclasses
import json
import sys
import time
from pathlib import Path

import yaml

from ffn_sim.spheroid.cadherin_bonds import resolve_cadherin
from ffn_sim.spheroid.ligand_traction import resolve_active_traction
from ffn_sim.spheroid.observables import effective_radius
from ffn_sim.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.spheroid.proliferation import run_growth_pooled
from ffn_sim.spheroid.substrate import resolve_substrate

_CFG = Path(__file__).resolve().parents[1] / "configs" / "layer2_cbm.yaml"


def _device(kind: str):
    import hoomd
    return hoomd.device.GPU(notice_level=0) if kind == "gpu" else hoomd.device.CPU(notice_level=0)


def _one_arm(resolved, prolif, cad, sub, *, n, seed, f_active, dt_factor, device, total_time,
             epoch_steps=1200, max_cells=None):
    """Run one arm; dt_factor>1 => smaller dt + proportionally more steps (same biological time)."""
    r = resolved
    es = epoch_steps
    if dt_factor != 1:
        r = dataclasses.replace(resolved, dt_cfl=resolved.dt_cfl / dt_factor)
        es = epoch_steps * dt_factor
    mc = max_cells if max_cells is not None else int(2.5 * n) + 1000
    t0 = time.perf_counter()
    res = run_growth_pooled(
        r, prolif, n_cells_init=n, total_time=total_time, seed=seed,
        cohesion="catch", cad=cad, substrate=sub, f_traction=f_active, max_cells=mc,
        epoch_steps=es,
    )
    return {
        "f_active_nN": round(f_active * 1e9, 2),
        "dt_factor": dt_factor,
        "R0_um": round(float(effective_radius(res["a0"]) * 1e6), 1),
        "aa0_core": round(float(res["area_core_over_a0"][-1]), 3),
        "aa0_raw": round(float(res["area_raw_over_a0"][-1]), 3),
        "n_final": int(res["n_cells"][-1]),
        "ejected": bool(res.get("ejected", False)),
        "capped": bool(res.get("capped_at_max_cells", False)),
        "wall_s": round(time.perf_counter() - t0, 1),
    }


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    n, seed = int(args[0]), int(args[1])
    kind = args[2] if len(args) > 2 else "cpu"
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = resolve_cadherin(resolved)
    sub = resolve_substrate(resolved)              # common substrate (adhesion_ratio=1.0)
    at = resolve_active_traction(resolved)
    total_time = 2.0 * prolif.cycle_time_mean
    device = _device(kind)

    arms = [
        ("baseline", 0.0, 1),
        ("armA_wholecell", at.f_wholecell, 1),
        ("armB_protrusion", at.f_protrusion, 5),   # PI-authorized dt/5 for the > ceiling arm
    ]
    for name, f, dtf in arms:
        rec = _one_arm(resolved, prolif, cad, sub, n=n, seed=seed, f_active=f, dt_factor=dtf,
                       device=device, total_time=total_time)
        rec = {"arm": name, "n": n, "seed": seed, "device": kind, **rec}
        print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
