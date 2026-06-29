"""Stage 3 (force stack): LJ neighbour-list choice — md.nlist.Tree vs Cell.

The build uses md.nlist.Tree(buffer=0.5 sigma) for the excluded-volume LJ (the
606 us / 45% slice of the post-integrator full-cell step), though cortex.py:955
notes "Cell is fine here". This times the LJ-only step (full-cell topology) with
Tree vs Cell on the SAME config to recommend the faster list. Config-only — a
build recommendation, edits no build file.

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_lj_nlist --n-fil 1000 --steps 800
"""

from __future__ import annotations

import argparse
import time
from copy import deepcopy
from pathlib import Path

import hoomd
import hoomd.md as md
import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin

CFG = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"


def _time(sim, steps):
    sim.run(40)
    t0 = time.perf_counter(); sim.run(steps)
    return 1e6 * (time.perf_counter() - t0) / steps


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--steps", type=int, default=800)
    args = ap.parse_args()
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = args.n_fil
    cfg["cortex"]["demo_mode"] = True
    cfg.setdefault("cortex", {}).setdefault("myosin", {})["batch_steps"] = 20000
    p = resolve_h3_derived(cfg)
    p_myo = resolve_cortex_myosin(cfg, dt=p.dt_cfl, R_cell=p.R_cell)
    p_xl = resolve_crosslinkers(cfg, dt=p.dt_cfl)
    hw = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=hoomd.device.GPU(notice_level=0),
        with_baoab=False, constrained=False, connected_mesh=True,
        rng=np.random.default_rng(1))
    sim = hw["sim"]
    integ = sim.operations.integrator
    lj_tree = [f for f in integ.forces if isinstance(f, md.pair.LJ)]
    if not lj_tree:
        print("no LJ in build"); return 2
    lj_tree = lj_tree[0]
    N = int(sim.state.get_snapshot().particles.N)
    types = list(sim.state.particle_types)

    # isolate the LJ, time with the build's Tree list
    integ.forces = [lj_tree]
    t_tree = _time(sim, args.steps)

    # rebuild the SAME LJ on a Cell list, copy every type-pair param + r_cut
    cell_nl = md.nlist.Cell(buffer=0.5 * p.lj_sigma)
    lj_cell = md.pair.LJ(nlist=cell_nl, default_r_cut=0.0)
    for i, a in enumerate(types):
        for b in types[i:]:
            lj_cell.params[(a, b)] = lj_tree.params[(a, b)]
            lj_cell.r_cut[(a, b)] = lj_tree.r_cut[(a, b)]
    integ.forces = [lj_cell]
    t_cell = _time(sim, args.steps)

    print(f"LJ nlist: N={N} types={len(types)}  (full-cell connected mesh)", flush=True)
    print(f"[LJ-only step]  Tree {t_tree:7.1f} us  |  Cell {t_cell:7.1f} us  "
          f"({t_tree/max(t_cell,1e-9):.2f}x {'Cell faster' if t_cell < t_tree else 'Tree faster'})",
          flush=True)
    rec = "Cell" if t_cell < t_tree * 0.97 else "Tree (keep)"
    print(f"[recommend] {rec}  (saves {max(t_tree-t_cell,0):.0f} us/step if switched)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
