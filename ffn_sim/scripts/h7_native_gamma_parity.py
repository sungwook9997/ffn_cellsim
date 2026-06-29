"""gamma (bond-tension) statistical parity: native vs cupy constrained BAOAB.

The native (splitmix64 counter-RNG) and cupy (default_rng stream) draw DIFFERENT
noise, so trajectories diverge — but both sample the same equilibrium, so the
ENSEMBLE bond-tension statistic that feeds the method-of-planes gamma must agree
within sampling error. Collects the per-step M-SHAKE Lagrange multipliers lambda
(T_bond = lambda*r0/dt) over a warmed cortex from BOTH integrators and compares
the mean / rms tension (Welch two-sample on the per-step means).

  PYTHONPATH="$HOME/ffn_cellsim/native/ffn_hoomd_plugin/build:\
$HOME/ffn_cellsim/native/ffn_hoomd_plugin/python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_gamma_parity --n-fil 300 --samples 400
"""

from __future__ import annotations

import argparse
import math
from copy import deepcopy
from pathlib import Path

import hoomd
import hoomd.md as md
import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.cortex.cortex import resolve_h3_derived
from ffn_sim.archive.hoomd_legacy.integrator.constrained_baoab import make_constrained_baoab_updater
from ffn_sim.scripts.h7_native_constrained_bench import _warm_and_state, _make_sim

CFG = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"


def _collect_cupy(p, snap, dt, seed, pairs, chains, l0, n_eq, samples, every):
    sim = _make_sim(p, snap, dt, seed)
    act, upd = make_constrained_baoab_updater(
        kT=p.kT, gamma={"actin_cortex": p.gamma_b}, dt=dt, constraint_pairs=pairs,
        constraint_lengths=np.full(pairs.shape[0], l0), chains=chains, seed=seed,
        shake_tol=1.0e-9, shake_max_iter=200, record_lambda=True)
    sim.operations.updaters.append(upd)
    sim.run(0); sim.run(n_eq)
    means, rms = [], []
    for _ in range(samples):
        sim.run(every)
        lam = np.asarray(act.lambda_buf, dtype=np.float64).ravel()
        means.append(lam.mean()); rms.append(np.sqrt((lam * lam).mean()))
    fac = l0 / dt
    return np.array(means) * fac, np.array(rms) * fac


def _collect_native(p, snap, dt, seed, chains_flat, F, N, l0, n_eq, samples, every):
    from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater
    sim = _make_sim(p, snap, dt, seed)
    nat = NativeConstrainedBaoabUpdater(
        dt=dt, kT=p.kT, inv_gamma_by_tag=[1.0 / p.gamma_b] * (F * N), chains_tag=chains_flat,
        n_chains=F, bonds_per_chain=N - 1, rest_length=l0, seed=seed, tol=1.0e-9,
        max_iter=200, trigger=hoomd.trigger.Periodic(1))
    sim.operations.updaters.append(nat)
    sim.run(0); sim.run(n_eq)
    means, rms = [], []
    for _ in range(samples):
        sim.run(every)
        lam = np.asarray(nat.lambda_buf, dtype=np.float64).ravel()
        means.append(lam.mean()); rms.append(np.sqrt((lam * lam).mean()))
    fac = l0 / dt
    return np.array(means) * fac, np.array(rms) * fac


def _welch(a, b):
    ma, mb = a.mean(), b.mean()
    va, vb = a.var(ddof=1), b.var(ddof=1)
    se = math.sqrt(va / a.size + vb / b.size)
    return ma, mb, (ma - mb) / se if se > 0 else 0.0, se


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=300)
    ap.add_argument("--warm", type=int, default=30000)
    ap.add_argument("--n-eq", type=int, default=2000)
    ap.add_argument("--samples", type=int, default=400)
    ap.add_argument("--every", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = args.n_fil
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    l0 = p.rest_length
    tau_bend = p.gamma_b * l0 ** 3 / p.bending_modulus
    dt = 0.03 * tau_bend
    snap, topology, N, F, n_part = _warm_and_state(p, args.n_fil, args.seed, args.warm)
    pairs = np.asarray(topology.bond_groups, dtype=np.int64)
    chains = [np.arange(f * N, (f + 1) * N, dtype=np.int64) for f in range(F)]
    chains_flat = np.concatenate(chains).astype(np.int32)

    cm, cr = _collect_cupy(p, snap, dt, args.seed, pairs, chains, l0,
                           args.n_eq, args.samples, args.every)
    nm, nr = _collect_native(p, snap, dt, args.seed, chains_flat, F, N, l0,
                             args.n_eq, args.samples, args.every)

    print(f"gamma-parity: n_fil={args.n_fil} samples={args.samples} (every {args.every} steps)",
          flush=True)
    for name, c, n in (("mean T_bond", cm, nm), ("rms  T_bond", cr, nr)):
        mc, mn, z, se = _welch(c, n)
        rel = abs(mc - mn) / (abs(mc) + 1e-30)
        ok = abs(z) < 3.0
        print(f"[{name}] cupy={mc:.4e} native={mn:.4e} rel={rel:.2%} "
              f"Welch_z={z:+.2f} (|z|<3 {'PASS' if ok else 'FAIL'})", flush=True)
    mc, mn, z, _ = _welch(cr, nr)
    print("gamma statistical parity " + ("PASS" if abs(z) < 3.0 else "FAIL"), flush=True)
    return 0 if abs(z) < 3.0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
