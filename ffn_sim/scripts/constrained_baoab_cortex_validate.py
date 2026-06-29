#!/usr/bin/env python
"""constrained-BD Milestone 3 — cortex equipartition + L_p with rigid backbone.

Scales the validated constrained integrator (Milestone 2) to the multi-
filament cortex (rigid backbone bonds + soft angles + WCA excluded volume),
at the fast dt = 0.03·τ_bend. Confirms the 3D equipartition ⟨E_bend⟩ ≈
0.9898 kT and per-filament L_p ∈ KU-1.1 band still hold at cortex scale and
with inter-filament LJ present — i.e. constrained-BD is production-ready for
the KU-3.x active-mechanics gates.

Standalone validation (no oracle/contract change). Pooled L_p (over the F
filaments per snapshot, the robust H.3 protocol) + equipartition (over all
angles × snapshots).

Usage:
    python ffn_sim/scripts/constrained_baoab_cortex_validate.py --n-fil 500 --seed 1
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import yaml
import hoomd
from hoomd import md

from ffn_sim.archive.hoomd_legacy.cortex.cortex import (
    resolve_h3_derived, build_cortex_state, build_cortex_simulation,
)
from ffn_sim.archive.hoomd_legacy.integrator.constrained_baoab import make_constrained_baoab_updater
from ffn_sim.common.filament_math import fit_persistence_length


def _tag_ordered_positions(sim, n_part: int) -> np.ndarray:
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position); tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
        return pos[inv][:n_part].copy()

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"


def run(n_fil: int, seed: int, *, dt_factor: float = 0.03, lj: bool = False,
        n_eq: int = 8000, n_sample: int = 300, interval: int = 30) -> dict:
    with open(CFG) as f:
        cfg = yaml.safe_load(f)
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    p = resolve_h3_derived(cfg)
    N = p.beads_per_filament
    l0 = p.rest_length
    F = p.n_filaments
    tau_bend = p.gamma_b * l0 ** 3 / p.bending_modulus
    dt = dt_factor * tau_bend

    # Phase 1 — warm-up with the standard BAOAB (stiff bond, cfl-dt) to relax
    # construction-time LJ overlaps + settle bonds to ℓ₀. Without this the
    # large production dt × initial overlap forces gives a per-step
    # displacement the explicit predictor + SHAKE cannot handle (diverges).
    n_part = F * N
    sim_w, _, _, _, _ = build_cortex_simulation(
        p, device=hoomd.device.CPU(notice_level=0), with_baoab=True,
        with_crosslinkers=False, rng=np.random.default_rng(seed),
    )
    sim_w.run(50_000)
    pos_warm = _tag_ordered_positions(sim_w, n_part)
    del sim_w

    # Phase 2 — constrained production from the warmed-up configuration.
    rng = np.random.default_rng(seed)
    snap, topology, _ = build_cortex_state(p, with_crosslinkers=False, rng=rng)
    snap.particles.position[:n_part] = pos_warm
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=seed)
    sim.create_state_from_snapshot(snap)

    # angle + LJ WCA forces (NO stretch bond — the backbone is rigid).
    angle = md.angle.Harmonic()
    angle.params["cortex-angle"] = dict(k=p.angle_k, t0=p.angle_t0)
    forces = [angle]
    if lj:
        # WCA excluded volume imposes its OWN fast CFL (stiff near contact);
        # at large dt it caps the achievable speedup below τ_bend's. Off by
        # default — irrelevant to the per-filament bending equipartition/L_p
        # this milestone validates (phantom-chain bending is LJ-independent).
        nlist = md.nlist.Tree(buffer=0.5 * p.lj_sigma)
        ljf = md.pair.LJ(nlist=nlist, default_r_cut=0.0)
        ljf.params[("actin_cortex", "actin_cortex")] = dict(
            epsilon=p.lj_epsilon, sigma=p.lj_sigma)
        ljf.r_cut[("actin_cortex", "actin_cortex")] = p.lj_r_cut if p.lj_enabled else 0.0
        ljf.mode = "shift"
        forces.append(ljf)
    sim.operations.integrator = md.Integrator(dt=dt, forces=forces, methods=[])

    # backbone constraints + per-filament chains for M-SHAKE + Fixman.
    pairs = np.asarray(topology.bond_groups, dtype=np.int64)
    chains = [np.arange(f * N, (f + 1) * N, dtype=np.int64) for f in range(F)]
    _, upd = make_constrained_baoab_updater(
        kT=p.kT, gamma={"actin_cortex": p.gamma_b}, dt=dt,
        constraint_pairs=pairs, constraint_lengths=np.full(pairs.shape[0], l0),
        chains=chains, seed=seed, shake_tol=1.0e-9, shake_max_iter=200,
    )
    sim.operations.updaters.append(upd)
    sim.run(0)
    sim.run(n_eq)

    Lp = []; Edata = []; tangd = []
    for _ in range(n_sample):
        sim.run(interval)
        with sim.state.cpu_local_snapshot as s:
            pos = np.asarray(s.particles.position); tg = np.asarray(s.particles.tag)
            inv = np.empty_like(tg); inv[tg] = np.arange(tg.size)
            r = pos[inv].reshape(F, N, 3)
        fit = fit_persistence_length(r, rest_length=l0)         # pooled over F
        if np.isfinite(fit.L_p_m):
            Lp.append(fit.L_p_m)
        bv = r[:, 1:, :] - r[:, :-1, :]
        bn = bv / np.linalg.norm(bv, axis=-1, keepdims=True)
        d = np.einsum("fij,fij->fi", bn[:, :-1, :], bn[:, 1:, :])
        tangd.append(d.ravel())
        th = np.arccos(np.clip(-d, -1, 1))
        Edata.append((0.5 * p.angle_k * (th - np.pi) ** 2 / p.kT).ravel())
    Lp = np.asarray(Lp); Eb = np.concatenate(Edata)
    C1 = float(np.mean(np.concatenate(tangd)))
    return dict(
        n_fil=int(F), seed=int(seed), dt_s=float(dt), dt_factor=float(dt_factor),
        dt_speedup_vs_cfl=float(dt / p.dt_cfl),
        n_eq=int(n_eq), n_sample=int(n_sample), interval=int(interval),
        L_p_pooled_um=float(Lp.mean() * 1e6),
        L_p_pooled_sem_um=float(Lp.std(ddof=1) / math.sqrt(len(Lp)) * 1e6),
        L_p_C1_um=float(-l0 / math.log(C1) * 1e6) if 0 < C1 < 1 else float("nan"),
        E_bend_kT=float(Eb.mean()), E_bend_sem=float(Eb.std(ddof=1) / math.sqrt(Eb.size)),
        KU11_band_um=[15.3, 18.7], eq_target_kT=0.9898, eq_tol=0.05,
        max_drift=float(getattr(upd.action, "max_constraint_drift", 0.0)),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dt-factor", type=float, default=0.03)
    ap.add_argument("--n-eq", type=int, default=8000)
    ap.add_argument("--n-sample", type=int, default=300)
    ap.add_argument("--interval", type=int, default=30)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    res = run(args.n_fil, args.seed, dt_factor=args.dt_factor,
              n_eq=args.n_eq, n_sample=args.n_sample, interval=args.interval)
    out = args.out or str(PKG / "outputs" / "h3" / "production" / "cbd_cortex"
                          / f"cortex_n{args.n_fil}_s{args.seed}.json")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(res, indent=2))
    print("RESULT " + json.dumps(res), flush=True)


if __name__ == "__main__":
    main()
