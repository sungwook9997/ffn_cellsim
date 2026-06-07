"""Stage 4 (step-count): largest STABLE dt for the native constrained full-cell.

The constrained integrator's purpose is to remove the stiff stretch-bond CFL so a
much larger dt is admissible (the cortex-only bench ran at 1607x CFL). If the
full-cell still runs at the bond-CFL dt, the constraint's dt headroom is unused —
a pure step-count win (fewer steps for the same physical time). This resets to a
soft-started config and sweeps dt upward, recording where the native run stays
stable (nonconverged==0, finite, bond drift small). Read-only; edits no build file.

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_dt_sweep --n-fil 1000
"""

from __future__ import annotations

import argparse
import time

import hoomd
import numpy as np

from ffn_sim.scripts.h7_native_fullcell_bench import _build


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--check", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mults", type=str, default="1,3,10,30,100,300,1000")
    args = ap.parse_args()
    from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater

    hw, p = _build(args.n_fil, args.seed)
    sim = hw["sim"]
    act = hw["baoab_action"]
    base_dt = float(act.dt)
    chains = list(act._chains_tag)
    F = len(chains); m = int(chains[0].shape[0]) - 1
    rest = float(act._chain_rest_length)
    inv_g = np.asarray(act._xp.asnumpy(act._inv_gamma_by_tag) if act._on_gpu
                       else act._inv_gamma_by_tag).astype(np.float64)
    chains_flat = np.concatenate([np.asarray(c) for c in chains]).astype(np.int32)
    N = int(sim.state.get_snapshot().particles.N)
    sim.operations.updaters.remove(hw["baoab_updater"])
    relaxed = sim.state.get_snapshot()  # soft-started config
    print(f"full-cell n_fil={args.n_fil} N={N}  base dt={base_dt:.3e} (build dt_cfl)", flush=True)

    prev = None
    best = None
    for mult in [int(x) for x in args.mults.split(",")]:
        dt = base_dt * mult
        sim.state.set_snapshot(relaxed)
        if prev is not None:
            sim.operations.updaters.remove(prev)
        sim.operations.integrator.dt = dt
        nat = NativeConstrainedBaoabUpdater(
            dt=dt, kT=float(act.kT), inv_gamma_by_tag=inv_g.tolist(), chains_tag=chains_flat,
            n_chains=F, bonds_per_chain=m, rest_length=rest, seed=int(act._seed),
            tol=1e-9, max_iter=200, trigger=hoomd.trigger.Periodic(1))
        sim.operations.updaters.insert(0, nat)
        prev = nat
        try:
            sim.run(0); sim.run(50)  # settle at this dt
            t0 = time.perf_counter(); sim.run(args.check); wall = time.perf_counter() - t0
            pos = sim.state.get_snapshot().particles.position
            finite = bool(np.all(np.isfinite(np.asarray(pos))))
            nc = int(nat.nonconverged_count)
            cf = chains_flat.reshape(F, m + 1)
            s = np.asarray(pos)[cf[:, :-1]] - np.asarray(pos)[cf[:, 1:]]
            drift = float(np.max(np.abs(np.sqrt(np.sum(s * s, 2)) - rest)) / rest)
            sps = args.check / wall
            ok = finite and nc == 0 and drift < 1e-6
            tag = "STABLE" if ok else "UNSTABLE"
            print(f"  dt={dt:.3e} ({mult:>4d}x cfl)  {tag}  nonconv={nc} drift={drift:.2e} "
                  f"finite={finite}  {sps:.0f} steps/s", flush=True)
            if ok:
                best = (mult, dt, sps)
            else:
                break
        except Exception as e:
            print(f"  dt={dt:.3e} ({mult:>4d}x cfl)  DIVERGED ({type(e).__name__})", flush=True)
            break

    if best:
        mult, dt, sps = best
        eq_phys = "same physical time at %dx fewer steps" % mult
        print(f"[dt-sweep] max stable dt = {dt:.3e} = {mult}x build dt_cfl  ({sps:.0f} steps/s)", flush=True)
        print(f"[dt-sweep] step-count lever: {mult}x ({eq_phys}) on top of the {('6.3' if mult else '')}x "
              f"integrator throughput", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
