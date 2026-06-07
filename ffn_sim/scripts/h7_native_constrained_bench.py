"""Stable real-cortex throughput: native vs cupy constrained BAOAB (Stage 1b).

Uses the validated stability pattern (constrained_baoab_cortex_validate):
  Phase 1 — warm-up with the unconstrained BAOAB (stiff bond, CFL dt) to relax
            construction-time overlaps (else the large production dt diverges);
  Phase 2 — constrained production (M-SHAKE + Fixman, NO stretch bond) at the
            large dt = dt_factor·τ_bend, with the angle force stack.
Times the cupy ConstrainedLM-BAOAB then the native FFNConstrainedBaoabUpdater on
the SAME warmed config + topology, and checks the native keeps the rigid bonds
satisfied (nonconverged=0, bond drift). Read-only; edits no build file.

  PYTHONPATH="$HOME/ffn_cellsim/native/ffn_hoomd_plugin/build:\
$HOME/ffn_cellsim/native/ffn_hoomd_plugin/python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_constrained_bench --n-fil 1000 --steps 2000
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

from ffn_sim.cortex.cortex import (
    resolve_h3_derived, build_cortex_state, build_cortex_simulation,
)
from ffn_sim.integrator.constrained_baoab import make_constrained_baoab_updater

CFG = Path(__file__).resolve().parents[1] / "configs" / "phase1_h3.yaml"


def _tag_pos(sim, n):
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position); tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        return pos[inv][:n].copy()


def _warm_and_state(p, n_fil, seed, warm_steps):
    N = p.beads_per_filament
    F = p.n_filaments
    n_part = F * N
    # Warm on CPU: soft-start (clipped Brownian relaxes severe construction
    # overlaps that the raw BAOAB cannot) + BAOAB settle. Phase 2 runs GPU.
    from ffn_sim.ecm.equilibrate import equilibrate_no_shear
    simw, updater, action, topology, _ = build_cortex_simulation(
        p, device=hoomd.device.CPU(notice_level=0), with_baoab=True,
        with_crosslinkers=False, rng=np.random.default_rng(seed))
    equilibrate_no_shear(
        simw, action, updater, n_softstart=500, n_baoab=4000,
        rest_length=p.rest_length, gamma_b=p.gamma_b)
    if warm_steps > 0:
        simw.run(warm_steps)
    pos_warm = _tag_pos(simw, n_part)
    del simw
    snap, topology, _ = build_cortex_state(p, with_crosslinkers=False,
                                           rng=np.random.default_rng(seed))
    snap.particles.position[:n_part] = pos_warm
    return snap, topology, N, F, n_part


def _make_sim(p, snap, dt, seed):
    sim = hoomd.Simulation(device=hoomd.device.GPU(notice_level=0), seed=seed)
    sim.create_state_from_snapshot(snap)
    angle = md.angle.Harmonic()
    angle.params["cortex-angle"] = dict(k=p.angle_k, t0=p.angle_t0)
    sim.operations.integrator = md.Integrator(dt=dt, forces=[angle], methods=[])
    return sim


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--warm", type=int, default=50000)
    ap.add_argument("--dt-factor", type=float, default=0.03)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = args.n_fil
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    l0 = p.rest_length
    tau_bend = p.gamma_b * l0 ** 3 / p.bending_modulus
    dt = args.dt_factor * tau_bend

    snap, topology, N, F, n_part = _warm_and_state(p, args.n_fil, args.seed, args.warm)
    pairs = np.asarray(topology.bond_groups, dtype=np.int64)
    chains = [np.arange(f * N, (f + 1) * N, dtype=np.int64) for f in range(F)]
    chains_flat = np.concatenate(chains).astype(np.int32)
    print(f"cortex: n_part={n_part} F={F} x {N} beads  dt={dt:.3e} "
          f"(={dt/p.dt_cfl:.0f}x CFL)  warmed {args.warm} steps", flush=True)

    res = {}
    # --- cupy constrained ---
    sim = _make_sim(p, snap, dt, args.seed)
    _, upd = make_constrained_baoab_updater(
        kT=p.kT, gamma={"actin_cortex": p.gamma_b}, dt=dt, constraint_pairs=pairs,
        constraint_lengths=np.full(pairs.shape[0], l0), chains=chains, seed=args.seed,
        shake_tol=1.0e-9, shake_max_iter=200)
    sim.operations.updaters.append(upd)
    sim.run(0); sim.run(200)
    t0 = time.perf_counter(); sim.run(args.steps); wall = time.perf_counter() - t0
    res["cupy"] = args.steps / wall
    print(f"[bench] cupy   constrained {res['cupy']:.1f} steps/s ({1e6/res['cupy']:.1f} us)",
          flush=True)

    # --- native constrained (same warmed config) ---
    from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater
    sim2 = _make_sim(p, snap, dt, args.seed)
    nat = NativeConstrainedBaoabUpdater(
        dt=dt, kT=p.kT, inv_gamma_by_tag=[1.0 / p.gamma_b] * n_part, chains_tag=chains_flat,
        n_chains=F, bonds_per_chain=N - 1, rest_length=l0, seed=args.seed, tol=1.0e-9,
        max_iter=200, trigger=hoomd.trigger.Periodic(1))
    sim2.operations.updaters.append(nat)
    sim2.run(0); sim2.run(200)
    t0 = time.perf_counter(); sim2.run(args.steps); wall = time.perf_counter() - t0
    res["native"] = args.steps / wall
    # constraint satisfaction + λ availability
    pos = _tag_pos(sim2, n_part)
    s = pos[chains_flat.reshape(F, N)[:, :-1]] - pos[chains_flat.reshape(F, N)[:, 1:]]
    drift = float(np.max(np.abs(np.sqrt(np.sum(s * s, 2)) - l0)) / l0)
    lam = np.asarray(nat.lambda_buf)
    print(f"[bench] native constrained {res['native']:.1f} steps/s ({1e6/res['native']:.1f} us)  "
          f"nonconv={nat.nonconverged_count} bond_drift={drift:.2e} lambda_shape={lam.shape}",
          flush=True)
    print(f"[bench] constrained throughput speedup native/cupy = "
          f"{res['native']/res['cupy']:.2f}x", flush=True)
    print(f"[bench] 2e8 steps @ dt={dt:.2e}s: cupy {2e8/res['cupy']/3600:.1f} h -> "
          f"native {2e8/res['native']/3600:.1f} h", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
