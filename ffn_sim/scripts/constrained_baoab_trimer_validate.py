#!/usr/bin/env python
"""Empirical Fixman-correctness gate for constrained-BD (trimer angle PDF).

The decisive Milestone-1 validation (constrained_baoab.py §Sanity Gate §6):
a rigid-bond + soft-angle trimer integrated with the constrained L-M Action
(+ Fixman pseudo-force) must reproduce the SAME bending-angle distribution
as the stiff-harmonic-bond + soft-angle trimer integrated with the frozen
baoab.py — because the H.2/H.3 gates are anchored on the stiff-spring
marginal. The Fixman correction is exactly what converts the rigid-chain
equilibrium into the stiff-spring one (Fixman 1974; Hinch 1994).

Runs ONE (mode, seed) per invocation so many can be fanned out as
independent concurrent processes (one per core) — see
[[feedback-cpu-jobs-gbook-concurrent]]. Aggregate the per-seed JSONs to get
ensemble mean ± stderr of the bending energy and compare modes.

Usage:
    python ffn_sim/scripts/constrained_baoab_trimer_validate.py \
        --mode fixman --seed 7 --n-step 2000000 --out OUT.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import hoomd
from hoomd import md

from ffn_sim.archive.hoomd_legacy.integrator.baoab import make_baoab_updater
from ffn_sim.archive.hoomd_legacy.integrator.constrained_baoab import (
    FIXMAN_SIGN, make_constrained_baoab_updater,
)

PKG = Path(__file__).resolve().parents[1]


def run(mode: str, seed: int, n_step: int, *,
        L0=1.0, kT=1.0, gamma=1.0, kth=2.0, kstr=2000.0, dt=1e-5,
        burn=100_000, samples=400) -> dict:
    snap = hoomd.Snapshot()
    snap.particles.N = 3
    snap.particles.types = ["A"]
    snap.particles.position[:] = [[-L0, 0, 0], [0, 0, 0], [L0, 0, 0]]
    snap.particles.typeid[:] = [0, 0, 0]
    snap.configuration.box = [100, 100, 100, 0, 0, 0]
    snap.bonds.N = 2; snap.bonds.types = ["b"]
    snap.bonds.group[:] = [[0, 1], [1, 2]]; snap.bonds.typeid[:] = [0, 0]
    snap.angles.N = 1; snap.angles.types = ["a"]
    snap.angles.group[:] = [[0, 1, 2]]; snap.angles.typeid[:] = [0]
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(snap)

    forces = []
    angle = md.angle.Harmonic(); angle.params["a"] = dict(k=kth, t0=np.pi)
    forces.append(angle)
    if mode == "harmonic":
        bond = md.bond.Harmonic(); bond.params["b"] = dict(k=kstr, r0=L0)
        forces.append(bond)
    sim.operations.integrator = md.Integrator(dt=dt, forces=forces, methods=[])

    if mode == "harmonic":
        _, upd = make_baoab_updater(kT=kT, gamma={"A": gamma}, dt=dt, seed=seed)
    elif mode in ("fixman", "nofixman"):
        chains = [np.array([0, 1, 2])] if mode == "fixman" else None
        _, upd = make_constrained_baoab_updater(
            kT=kT, gamma={"A": gamma}, dt=dt,
            constraint_pairs=np.array([[0, 1], [1, 2]]),
            constraint_lengths=np.array([L0, L0]), chains=chains, seed=seed,
        )
    else:
        raise ValueError(mode)

    sim.operations.updaters.append(upd)
    sim.run(0); sim.run(burn)
    th = np.empty(samples)
    per = n_step // samples
    for s in range(samples):
        sim.run(per)
        with sim.state.cpu_local_snapshot as snp:
            p = np.asarray(snp.particles.position); tg = np.asarray(snp.particles.tag)
            inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
            r = p[inv]
        b1 = r[0] - r[1]; b2 = r[2] - r[1]
        c = np.dot(b1, b2) / (np.linalg.norm(b1) * np.linalg.norm(b2))
        th[s] = np.arccos(np.clip(c, -1, 1))
    dev2 = (np.pi - th) ** 2
    return dict(
        mode=mode, seed=int(seed), n_step=int(n_step), samples=int(samples),
        kth=kth, fixman_sign=FIXMAN_SIGN,
        theta_mean=float(th.mean()),
        E_bend=float(0.5 * kth * dev2.mean()),
        E_bend_sem=float(0.5 * kth * dev2.std(ddof=1) / np.sqrt(samples)),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["harmonic", "fixman", "nofixman"], required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n-step", type=int, default=2_000_000)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    res = run(args.mode, args.seed, args.n_step)
    out = args.out or str(
        PKG / "outputs" / "h3" / "production"
        / f"cbd_trimer_{args.mode}_seed{args.seed}.json"
    )
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(res, indent=2))
    print("RESULT " + json.dumps(res), flush=True)


if __name__ == "__main__":
    main()
