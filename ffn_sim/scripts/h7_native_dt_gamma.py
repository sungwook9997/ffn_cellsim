"""Stage 4 physics check: does the 30x dt preserve gamma (bond tension)?

The 30x-dt step-count lever is only usable if the larger dt samples the SAME
equilibrium tension. T_bond = lambda*r0/dt is the physical (dt-normalised) bond
tension, so its ensemble rms must be dt-INVARIANT. Collects rms T_bond from the
native full-cell at 1x and 30x dt_cfl and compares (Welch). BAOAB-limit is
designed for correct configurational sampling at large dt, so this is expected
to pass; this makes it explicit on the full cell (LJ + binders included).

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_dt_gamma --n-fil 1000 --samples 250
"""

from __future__ import annotations

import argparse
import math

import hoomd
import numpy as np

from ffn_sim.scripts.h7_native_fullcell_bench import _build


def _collect(sim, nat, r0, dt, n_eq, samples, every):
    # n_eq is in PHYSICAL-time-matched steps (scaled 1/mult by the caller) so 1x
    # and 30x dt equilibrate to the SAME physical time before sampling.
    sim.run(n_eq)
    rms = []
    for _ in range(samples):
        sim.run(every)
        lam = np.asarray(nat.lambda_buf, dtype=np.float64).ravel()
        rms.append(np.sqrt((lam * lam).mean()) * r0 / dt)
    return np.array(rms)


def _welch(a, b):
    se = math.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size)
    return a.mean(), b.mean(), ((a.mean() - b.mean()) / se if se > 0 else 0.0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--samples", type=int, default=250)
    ap.add_argument("--every", type=int, default=300)
    ap.add_argument("--phys-eq", type=int, default=60000)
    ap.add_argument("--seed", type=int, default=1)
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
    cflat = np.concatenate([np.asarray(c) for c in chains]).astype(np.int32)
    sim.operations.updaters.remove(hw["baoab_updater"])
    relaxed = sim.state.get_snapshot()

    # Physical-time-matched equilibration: PHYS_EQ in dt_cfl(1x)-steps; the Nx-dt
    # run uses PHYS_EQ/N steps so both settle the SAME physical time. Same for
    # the sample spacing (every/N) so the two samplings span equal physical time.
    PHYS_EQ = args.phys_eq

    def run_dt(mult):
        dt = base_dt * mult
        sim.state.set_snapshot(relaxed)
        for u in list(sim.operations.updaters):
            if isinstance(u, NativeConstrainedBaoabUpdater):
                sim.operations.updaters.remove(u)
        sim.operations.integrator.dt = dt
        nat = NativeConstrainedBaoabUpdater(
            dt=dt, kT=float(act.kT), inv_gamma_by_tag=inv_g.tolist(), chains_tag=cflat,
            n_chains=F, bonds_per_chain=m, rest_length=rest, seed=int(act._seed),
            tol=1e-9, max_iter=200, trigger=hoomd.trigger.Periodic(1))
        sim.operations.updaters.insert(0, nat)
        sim.run(0)
        n_eq = max(1, PHYS_EQ // mult)
        every = max(1, args.every // mult)
        return _collect(sim, nat, rest, dt, n_eq, args.samples, every)

    r1 = run_dt(1)
    r30 = run_dt(30)
    m1, m30, z = _welch(r1, r30)
    rel = abs(m1 - m30) / (abs(m1) + 1e-30)
    ok = abs(z) < 3.0
    print(f"dt-gamma: n_fil={args.n_fil} samples={args.samples}", flush=True)
    print(f"[rms T_bond] 1x_dt={m1:.4e}  30x_dt={m30:.4e}  rel={rel:.2%}  Welch_z={z:+.2f}  "
          f"({'PASS' if ok else 'FAIL'} |z|<3)", flush=True)
    print("30x-dt gamma preservation " + ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
