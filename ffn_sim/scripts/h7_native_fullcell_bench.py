"""Full-cell end-to-end throughput: native vs cupy constrained BAOAB (Stage 1b).

Builds the FULL cell (build_cortex_full_simulation: cortex + myosin/xlink/integrin
updaters + ERM/turgor/membrane custom forces + connected mesh), soft-starts
(equilibrate_no_shear), then times the cupy ConstrainedLM-BAOAB vs the native
FFNConstrainedBaoabUpdater on the SAME topology. This is the real end-to-end
number the integrator micro-bench (44x) does NOT capture — the per-step custom
forces + (rare) binder host-syncs are still there. Read-only; edits no build file.

  PYTHONPATH="$HOME/ffn_cellsim/native/ffn_hoomd_plugin/build:\
$HOME/ffn_cellsim/native/ffn_hoomd_plugin/python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_fullcell_bench --n-fil 1000 --steps 1000
"""

from __future__ import annotations

import argparse
import time
from copy import deepcopy
from pathlib import Path

import hoomd
import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.cell.cell import build_cortex_full_simulation
from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.archive.hoomd_legacy.cortex.myosin import resolve_cortex_myosin
from ffn_sim.archive.hoomd_legacy.ecm.equilibrate import equilibrate_no_shear

CFG = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"


def _resolved(n_fil):
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = int(n_fil)
    cfg["cortex"]["demo_mode"] = True
    cfg.setdefault("cortex", {}).setdefault("myosin", {})
    cfg["cortex"]["myosin"]["stepping_mode"] = "grip_walk"
    cfg["cortex"]["myosin"]["mesoscale_force_scaling"] = True
    cfg["cortex"]["myosin"]["batch_steps"] = 20000
    p = resolve_h3_derived(cfg)
    p_myo = resolve_cortex_myosin(cfg, dt=p.dt_cfl, R_cell=p.R_cell)
    p_xl = resolve_crosslinkers(cfg, dt=p.dt_cfl)
    return p, p_myo, p_xl


def _build(n_fil, seed):
    p, p_myo, p_xl = _resolved(n_fil)
    hw = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=hoomd.device.GPU(notice_level=0),
        with_baoab=True, constrained=True, connected_mesh=True,
        rng=np.random.default_rng(seed))
    sim = hw["sim"]
    # soft-start: relax construction overlaps (clipped Brownian) before the hot loop.
    equilibrate_no_shear(sim, hw["baoab_action"], hw["baoab_updater"],
                         n_softstart=500, n_baoab=4000,
                         rest_length=p.rest_length, gamma_b=p.gamma_b)
    sim.run(0)
    return hw, p


def _time(sim, steps):
    sim.run(100)
    t0 = time.perf_counter(); sim.run(steps)
    return steps / (time.perf_counter() - t0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--steps", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    # cupy full-cell
    hw, p = _build(args.n_fil, args.seed)
    sim = hw["sim"]
    act = hw["baoab_action"]
    N = int(sim.state.get_snapshot().particles.N)
    chains = list(act._chains_tag)
    F = len(chains); m = int(chains[0].shape[0]) - 1
    uniform = len({int(c.shape[0]) for c in chains}) == 1
    print(f"full-cell: N={N} F={F} x {m+1} uniform={uniform} dt={act.dt:.3e}", flush=True)
    if not uniform:
        print("variable-length chains -> native (uniform) N/A here.", flush=True)
        return 2
    rest = float(act._chain_rest_length)
    inv_g = np.asarray(act._xp.asnumpy(act._inv_gamma_by_tag) if act._on_gpu
                       else act._inv_gamma_by_tag).astype(np.float64)
    sps_cupy = _time(sim, args.steps)
    print(f"[fullcell] cupy   {sps_cupy:.1f} steps/s ({1e6/sps_cupy:.1f} us)", flush=True)

    # native full-cell (swap baoab, same everything else)
    from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater
    hw2, _ = _build(args.n_fil, args.seed)
    sim2 = hw2["sim"]
    sim2.operations.updaters.remove(hw2["baoab_updater"])
    nat = NativeConstrainedBaoabUpdater(
        dt=float(act.dt), kT=float(act.kT), inv_gamma_by_tag=inv_g.tolist(),
        chains_tag=np.concatenate([np.asarray(c) for c in chains]).astype(np.int32),
        n_chains=F, bonds_per_chain=m, rest_length=rest, seed=int(act._seed),
        tol=1e-9, max_iter=200, trigger=hoomd.trigger.Periodic(1))
    sim2.operations.updaters.insert(0, nat)
    sim2.run(0)
    sps_nat = _time(sim2, args.steps)
    print(f"[fullcell] native {sps_nat:.1f} steps/s ({1e6/sps_nat:.1f} us)  "
          f"nonconv={nat.nonconverged_count}", flush=True)
    print(f"[fullcell] END-TO-END speedup native/cupy = {sps_nat/sps_cupy:.2f}x", flush=True)
    print(f"[fullcell] 2e8 steps: cupy {2e8/sps_cupy/3600:.1f} h -> native "
          f"{2e8/sps_nat/3600:.1f} h", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
