#!/usr/bin/env python
"""constrained-BD Milestone 2 — single-filament L_p re-validation (Fixman sign).

The DECISIVE empirical test of the constrained integrator + Fixman +1 sign
(the trimer could not arbitrate: its metric effect ±3-4% < seed noise). A
21-bead H.2 filament has 19 cumulative bending angles + the tight KU-1.1
L_p band [15.3, 18.7] μm, so the metric correction IS resolvable here.

Modes (one (mode, seed) per process → fan out concurrently on gbook):
  baseline   : stiff harmonic stretch bond + angle, standard BAOAB at dt_cfl.
  constrained: rigid bond (SHAKE) + angle + Fixman(+1), constrained BAOAB.
               --dt-mode cfl   → same dt as baseline (isolates Fixman).
               --dt-mode bend  → dt = safety·τ_bend (CFL-win, ~hundreds×).

Gate (no oracle/contract change — standalone validation only): L_p ∈ KU-1.1
band [15.3, 18.7] μm AND ⟨E_bend⟩ ∈ 0.9898 ±5% kT. If constrained+Fixman
reproduces baseline (and the band/target), the +1 sign is confirmed.

Usage:
  python ffn_sim/scripts/constrained_baoab_lp_validate.py \
      --mode constrained --dt-mode cfl --seed 1 --out OUT.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import yaml
import hoomd
import gsd.hoomd
from hoomd import md

from ffn_sim.scripts.h2_single_filament import resolve_h2_derived
from ffn_sim.integrator.baoab import make_baoab_updater
from ffn_sim.integrator.constrained_baoab import make_constrained_baoab_updater
from ffn_sim.common.filament_math import fit_persistence_length

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h2.yaml"


def _build_state(p) -> gsd.hoomd.Frame:
    N = p.beads_per_fiber
    pos = np.zeros((N, 3), dtype=np.float64)
    pos[:, 0] = np.linspace(-0.5 * p.L_fiber, 0.5 * p.L_fiber, N)
    snap = gsd.hoomd.Frame()
    snap.particles.N = N
    snap.particles.types = ["actin"]
    snap.particles.typeid = np.zeros(N, dtype=np.uint32)
    snap.particles.position = pos
    snap.particles.mass = np.ones(N)
    snap.bonds.N = N - 1
    snap.bonds.types = ["actin-bond"]
    snap.bonds.typeid = np.zeros(N - 1, dtype=np.uint32)
    snap.bonds.group = np.stack(
        [np.arange(N - 1, dtype=np.uint32), np.arange(1, N, dtype=np.uint32)], -1)
    snap.angles.N = N - 2
    snap.angles.types = ["actin-angle"]
    snap.angles.typeid = np.zeros(N - 2, dtype=np.uint32)
    snap.angles.group = np.stack(
        [np.arange(N - 2, dtype=np.uint32), np.arange(1, N - 1, dtype=np.uint32),
         np.arange(2, N, dtype=np.uint32)], -1)
    snap.configuration.box = [p.box_Lx, p.box_Ly, p.box_Lz, 0, 0, 0]
    return snap


def run(mode: str, dt_mode: str, seed: int, *,
        n_eq: int, n_sample: int, interval: int,
        dt_factor: float | None = None) -> dict[str, Any]:
    with open(CFG) as f:
        cfg = yaml.safe_load(f)
    p = resolve_h2_derived(cfg)
    N = p.beads_per_fiber
    l0 = p.rest_length

    dt = p.dt_cfl
    if mode == "constrained" and dt_mode == "bend":
        # rigid → CFL set by τ_bend, but the explicit-predictor + SHAKE step
        # also needs per-step displacement ≲ ℓ₀, which can cap dt below
        # cfl_safety·τ_bend. dt_factor overrides for the ceiling scan.
        # 0.03 is the validated ceiling: M-SHAKE stays machine-precision and
        # equipartition/L_p match H.2 (factor 0.1 = cfl_safety diverges — the
        # explicit-predictor displacement exceeds what SHAKE can project).
        f = dt_factor if dt_factor is not None else 0.03
        dt = f * p.tau_bend

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=seed)
    sim.create_state_from_snapshot(_build_state(p))
    angle = md.angle.Harmonic()
    angle.params["actin-angle"] = dict(k=p.angle_k, t0=p.angle_t0)

    if mode == "baseline":
        bond = md.bond.Harmonic()
        bond.params["actin-bond"] = dict(k=p.bond_k, r0=l0)
        sim.operations.integrator = md.Integrator(dt=dt, forces=[bond, angle], methods=[])
        _, upd = make_baoab_updater(kT=p.kT, gamma={"actin": p.gamma_b}, dt=dt, seed=seed)
    elif mode == "constrained":
        # NO stretch bond force — the bond is the rigid constraint.
        sim.operations.integrator = md.Integrator(dt=dt, forces=[angle], methods=[])
        pairs = np.stack([np.arange(N - 1), np.arange(1, N)], -1)
        _, upd = make_constrained_baoab_updater(
            kT=p.kT, gamma={"actin": p.gamma_b}, dt=dt,
            constraint_pairs=pairs, constraint_lengths=np.full(N - 1, l0),
            chains=[np.arange(N)], seed=seed,
            shake_tol=1.0e-9, shake_max_iter=3000,
        )
    else:
        raise ValueError(mode)
    sim.operations.updaters.append(upd)
    sim.run(0)
    sim.run(n_eq)

    # Collect ALL sample frames then fit the POOLED tangent correlation once
    # (H.2's robust protocol — per-snapshot single-filament fits are far too
    # noisy). C(1) = ⟨t̂_i·t̂_{i+1}⟩ gives the low-scatter local L_p_C1;
    # equipartition ⟨E_bend⟩ is the robust, fast-equilibrating primary gate.
    frames = np.empty((n_sample, N, 3))
    tang_dot = []   # adjacent unit-tangent dots → C(1)
    Ebend = []
    for k in range(n_sample):
        sim.run(interval)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position); tg = np.asarray(s.particles.tag)
            inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
            r = pos[inv]
        frames[k] = r
        bv = r[1:] - r[:-1]
        bn = bv / np.linalg.norm(bv, axis=1, keepdims=True)
        d = np.einsum("ij,ij->i", bn[:-1], bn[1:])     # t̂_i·t̂_{i+1} = cos(bend)
        tang_dot.append(d)
        th = np.arccos(np.clip(-d, -1, 1))             # interior angle (t0=π)
        Ebend.append(0.5 * p.angle_k * (th - np.pi) ** 2 / p.kT)
    fit = fit_persistence_length(frames, rest_length=l0)      # pooled over frames
    Eb = np.concatenate(Ebend)
    C1 = float(np.mean(np.concatenate(tang_dot)))
    L_p_C1_um = float(-l0 / math.log(C1) * 1e6) if 0.0 < C1 < 1.0 else float("nan")
    return dict(
        mode=mode, dt_mode=dt_mode, seed=int(seed), dt_s=float(dt),
        n_eq=int(n_eq), n_sample=int(n_sample), interval=int(interval),
        phys_time_s=float(n_sample * interval * dt),
        L_p_pooled_um=float(fit.L_p_m * 1e6),
        L_p_C1_um=L_p_C1_um, C1=C1,
        E_bend_kT=float(Eb.mean()), E_bend_sem=float(Eb.std(ddof=1) / math.sqrt(Eb.size)),
        KU11_band_um=[15.3, 18.7], L_p_C1_ref_um=16.56, eq_target_kT=0.9898, eq_tol=0.05,
        max_drift=float(getattr(upd.action, "max_constraint_drift", 0.0)),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["baseline", "constrained"], required=True)
    ap.add_argument("--dt-mode", choices=["cfl", "bend"], default="cfl")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n-eq", type=int, default=200_000)
    ap.add_argument("--n-sample", type=int, default=400)
    ap.add_argument("--interval", type=int, default=10_000)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    # For the bend-dt constrained run, match physical time, not step count.
    n_eq, interval = args.n_eq, args.interval
    res = run(args.mode, args.dt_mode, args.seed,
              n_eq=n_eq, n_sample=args.n_sample, interval=interval)
    tag = f"{args.mode}_{args.dt_mode}_s{args.seed}"
    out = args.out or str(PKG / "outputs" / "h3" / "production" / "cbd_lp" / f"{tag}.json")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(res, indent=2))
    print("RESULT " + json.dumps(res), flush=True)


if __name__ == "__main__":
    main()
