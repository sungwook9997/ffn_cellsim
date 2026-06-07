"""Scope the remaining ~2008 us/step of the everything-native Gate-A step.

Builds the everything-native cell (native integrator + native compartment forces,
the verified 6x ARM C), then prunes force/updater subsets (timing only) to split
the residual into: per-step binder updaters (myosin/xlink/integrin/ligand-pin
host-sync), the LJ excluded volume, and the integrator + built-in bond/angle.
Tells us whether the next ~Nx lever is the binder host-sync or the LJ wall.

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_gate_a_profile --n-fil 1000 --steps 500
"""

from __future__ import annotations

import argparse
import time
from copy import deepcopy

import numpy as np

from ffn_sim.cell.manifest import load_manifest
from ffn_sim.scripts.h7_native_fullcell_go import _build_two_phase, _to_host
from ffn_sim.scripts.h7_native_gate_a_verify import _native_integrator, _swap_compartment_forces


def _time(sim, steps):
    sim.run(40)
    t0 = time.perf_counter(); sim.run(steps)
    return 1e6 * (time.perf_counter() - t0) / steps


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--warmup", type=int, default=4000)
    ap.add_argument("--steps", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    import hoomd
    import hoomd.md as md
    dev = hoomd.device.GPU(notice_level=0)

    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})["n_filaments"] = args.n_fil
    build_constrained = _build_two_phase(manifest, device=dev, warmup=args.warmup,
                                         softstart=max(300, args.warmup // 8), seed=args.seed)
    cell = build_constrained()
    act = cell.baoab_action
    P = dict(dt=float(act.dt), kT=float(act.kT),
             inv_gamma=[float(x) for x in np.asarray(_to_host(act._inv_gamma_by_tag), np.float64)],
             chains_stacked=np.asarray(act.chains_tag_stacked),
             F=int(np.asarray(act.chains_tag_stacked).shape[0]),
             m=int(np.asarray(act.chains_tag_stacked).shape[1]) - 1,
             r0=float(act._chain_rest_length), seed=int(act._seed))
    from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater
    nat = _native_integrator(cell, P, hoomd)
    nswap = _swap_compartment_forces(cell)
    cell.simulation.run(0)

    sim = cell.simulation
    integ = sim.operations.integrator
    all_forces = list(integ.forces)
    all_upd = list(sim.operations.updaters)
    binders = [u for u in all_upd if not isinstance(u, NativeConstrainedBaoabUpdater)]
    lj = [f for f in all_forces if isinstance(f, md.pair.LJ)]
    N = int(sim.state.N_particles)
    print(f"everything-native: N={N} comp-swapped={nswap}  {len(all_forces)} forces "
          f"({len(lj)} LJ) + {len(binders)} non-integrator updaters", flush=True)

    t_full = _time(sim, args.steps)
    # remove the per-step binder updaters (keep native integrator)
    for u in binders:
        sim.operations.updaters.remove(u)
    t_no_binder = _time(sim, args.steps)
    # also drop the LJ
    integ.forces = [f for f in all_forces if not isinstance(f, md.pair.LJ)]
    t_no_binder_no_lj = _time(sim, args.steps)
    # integrator + bare bond/angle only (drop ALL forces)
    integ.forces = []
    t_int_only = _time(sim, args.steps)

    binder_us = t_full - t_no_binder
    lj_us = t_no_binder - t_no_binder_no_lj
    rest_forces_us = t_no_binder_no_lj - t_int_only
    print(f"[profile] full everything-native   {t_full:8.1f} us  ({1e6/t_full:.0f} steps/s)", flush=True)
    print(f"  binder updaters (per-step)        {binder_us:8.1f} us  ({100*binder_us/t_full:4.1f}%)  <- host-sync lever", flush=True)
    print(f"  LJ excluded volume                {lj_us:8.1f} us  ({100*lj_us/t_full:4.1f}%)  (HOOMD wall / mesh research)", flush=True)
    print(f"  bond/angle/compartment forces     {rest_forces_us:8.1f} us  ({100*rest_forces_us/t_full:4.1f}%)", flush=True)
    print(f"  integrator + updater triggers     {t_int_only:8.1f} us  ({100*t_int_only/t_full:4.1f}%)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
