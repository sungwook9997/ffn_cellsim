"""Scope Stage 2: decompose the native full-cell step (what is the remaining cost?).

After the native integrator removed the constrained-BAOAB wall, the full-cell step
is ~1.36 ms. This times the native full-cell with force subsets pruned (timing
only — physics is irrelevant here) to attribute the remainder to: native
integrator, built-in force stack (bond/angle/LJ, already C++), and the per-step
Python custom forces (ERM/turgor/membrane = md.force.Custom — the Stage-2 target).

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_fullcell_profile --n-fil 1000 --steps 600
"""

from __future__ import annotations

import argparse
import time

import hoomd
import hoomd.md as md
import numpy as np

from ffn_sim.scripts.h7_native_fullcell_bench import _build


def _is_custom(f) -> bool:
    return isinstance(f, md.force.Custom)


def _time(sim, steps):
    sim.run(60)
    t0 = time.perf_counter(); sim.run(steps)
    return 1e6 * (time.perf_counter() - t0) / steps  # us/step


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--steps", type=int, default=600)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater
    hw, p = _build(args.n_fil, args.seed)
    sim = hw["sim"]
    act = hw["baoab_action"]
    chains = list(act._chains_tag)
    F = len(chains); m = int(chains[0].shape[0]) - 1
    rest = float(act._chain_rest_length)
    inv_g = np.asarray(act._xp.asnumpy(act._inv_gamma_by_tag) if act._on_gpu
                       else act._inv_gamma_by_tag).astype(np.float64)
    sim.operations.updaters.remove(hw["baoab_updater"])
    nat = NativeConstrainedBaoabUpdater(
        dt=float(act.dt), kT=float(act.kT), inv_gamma_by_tag=inv_g.tolist(),
        chains_tag=np.concatenate([np.asarray(c) for c in chains]).astype(np.int32),
        n_chains=F, bonds_per_chain=m, rest_length=rest, seed=int(act._seed),
        tol=1e-9, max_iter=200, trigger=hoomd.trigger.Periodic(1))
    sim.operations.updaters.insert(0, nat)

    integ = sim.operations.integrator
    all_forces = list(integ.forces)
    custom = [f for f in all_forces if _is_custom(f)]
    builtin = [f for f in all_forces if not _is_custom(f)]
    names_c = [type(f).__name__ for f in custom]
    names_b = [type(f).__name__ for f in builtin]
    print(f"full-cell native: {len(all_forces)} forces = {len(builtin)} built-in "
          f"{names_b} + {len(custom)} custom {names_c}", flush=True)

    t_all = _time(sim, args.steps)
    integ.forces = builtin            # drop Python custom forces (timing only)
    t_nocustom = _time(sim, args.steps)
    integ.forces = []                 # integrator + updaters only
    t_none = _time(sim, args.steps)
    integ.forces = all_forces

    custom_us = t_all - t_nocustom
    builtin_us = t_nocustom - t_none
    print(f"[profile] full step          {t_all:8.1f} us  (native {1e6/t_all:.0f} steps/s)", flush=True)
    print(f"[profile]  custom forces     {custom_us:8.1f} us  ({100*custom_us/t_all:4.1f}%)  <- Stage 2 target", flush=True)
    print(f"[profile]  built-in forces   {builtin_us:8.1f} us  ({100*builtin_us/t_all:4.1f}%)  (LJ/angle/bond, already C++)", flush=True)
    print(f"[profile]  integrator+rest   {t_none:8.1f} us  ({100*t_none/t_all:4.1f}%)  (native BAOAB + updater triggers)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
