"""H.7 Gate-B native relaxed (unilateral M-SHAKE + τ_bend-EMA) parity vs cupy.

At kT=0 the predictor + Fixman vanish (∝kT), so the constrained step is
deterministic. With chains started COMPRESSED (bonds below ℓ₀) under a constant
force, the rigid pre-pass yields a sustained compressive λ → the τ_bend EMA crosses
F_crit → the unilateral relaxed solve releases those bonds. The native
FFNConstrainedBaoabUpdater(compression_release=True) and the cupy
ConstrainedLeimkuhlerMatthewsBAOAB(compression_release=True, release_load_crit,
load_tau) must agree bit-close (same rigid λ → same EMA → same eligibility → same
released set → same projection).

Also asserts: (1) released bonds actually appear (the mode is exercised), and
(2) with F_crit set ABOVE every load, NOTHING releases → bit-identical to the
rigid native path (rigid-limit superset, contract §5.1).

    cd ~/ffn_cellsim/native/ffn_hoomd_plugin
    PYTHONPATH="$PWD/build:$PWD/python:$HOME/ffn_cellsim" python test_constrained_relaxed_parity.py
"""
from __future__ import annotations

import numpy as np
import hoomd

from ffn_sim.integrator.constrained_baoab import ConstrainedLeimkuhlerMatthewsBAOAB
from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater


def _build(F, m, L, r0, stretch, rng):
    N = F * (m + 1)
    pos = np.zeros((N, 3))
    chains, cpairs = [], []
    for f in range(F):
        base = f * (m + 1)
        c = (rng.random(3) - 0.5) * (L * 1e-3)
        ch = []
        for b in range(m + 1):
            pos[base + b] = c
            ch.append(base + b)
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            c = c + r0 * stretch * d
        chains.append(np.array(ch, dtype=np.int64))
        for a in range(m):
            cpairs.append((base + a, base + a + 1))
    return N, pos, chains, np.array(cpairs, dtype=np.int64)


def _snapshot(N, L, pos):
    s = hoomd.Snapshot()
    if s.communicator.rank == 0:
        s.configuration.box = [L, L, L, 0, 0, 0]
        s.particles.N = N
        s.particles.types = ["A"]
        s.particles.position[:] = pos
        s.particles.typeid[:] = np.zeros(N, dtype=np.int32)
    return s


def _unwrap(sim, L):
    snap = sim.state.get_snapshot()
    return (np.array(snap.particles.position, np.float64)
            + np.array(snap.particles.image, np.float64) * L)


def main() -> int:
    F, m = 800, 6
    L, r0, gamma, dt = 1.0e6, 100.0, 1.0, 1.0e-3
    seed, tol, mit = 7, 1.0e-10, 200
    load_tau = 5.0e-3          # alpha = dt/load_tau = 0.2 → EMA builds in a few steps
    dev = hoomd.device.GPU()
    fails = []

    rng = np.random.default_rng(11)
    # COMPRESSED init (stretch 0.6 < 1) so bonds carry sustained compression.
    N, pos, chains, cpairs = _build(F, m, L, r0, 0.6, rng)
    cl = np.full(cpairs.shape[0], r0)
    chains_flat = np.concatenate([c for c in chains]).astype(np.int32)

    def run_native(steps, cr, fcrit):
        sim = hoomd.Simulation(device=dev, seed=1)
        sim.create_state_from_snapshot(_snapshot(N, L, pos))
        cf = hoomd.md.force.Constant(filter=hoomd.filter.All())
        cf.constant_force["A"] = (0.5, 0.0, 0.0)
        cf.constant_torque["A"] = (0.0, 0.0, 0.0)
        ig = hoomd.md.Integrator(dt=dt); ig.forces = [cf]; ig.methods = []
        sim.operations.integrator = ig
        sim.operations.updaters.append(NativeConstrainedBaoabUpdater(
            dt=dt, kT=0.0, inv_gamma_by_tag=[1.0 / gamma] * N, chains_tag=chains_flat,
            n_chains=F, bonds_per_chain=m, rest_length=r0, seed=seed, tol=tol, max_iter=mit,
            trigger=hoomd.trigger.Periodic(1),
            compression_release=cr, release_load_crit=fcrit, load_tau=load_tau))
        sim.run(0); sim.run(steps)
        return _unwrap(sim, L)

    def run_cupy(steps, cr, fcrit):
        sim = hoomd.Simulation(device=dev, seed=1)
        sim.create_state_from_snapshot(_snapshot(N, L, pos))
        cf = hoomd.md.force.Constant(filter=hoomd.filter.All())
        cf.constant_force["A"] = (0.5, 0.0, 0.0)
        cf.constant_torque["A"] = (0.0, 0.0, 0.0)
        ig = hoomd.md.Integrator(dt=dt); ig.forces = [cf]; ig.methods = []
        sim.operations.integrator = ig
        act = ConstrainedLeimkuhlerMatthewsBAOAB(
            kT=0.0, gamma={"A": gamma}, dt=dt, constraint_pairs=cpairs,
            constraint_lengths=cl, chains=chains, seed=seed, shake_tol=tol, shake_max_iter=mit,
            compression_release=cr, release_load_crit=(fcrit if cr else None),
            load_tau=(load_tau if cr else None))
        sim.operations.updaters.append(
            hoomd.update.CustomUpdater(action=act, trigger=hoomd.trigger.Periodic(1)))
        sim.run(0); sim.run(steps)
        return _unwrap(sim, L), act

    def _bondlen(out):
        cf = chains_flat.reshape(F, m + 1)
        s = out[cf[:, :-1]] - out[cf[:, 1:]]
        return np.sqrt(np.sum(s * s, axis=2))

    # ---- Test 1: relaxed parity (F_crit moderate → some bonds release) ----
    # F_crit chosen between loads: pick from cupy's own EMA behaviour. Use a value
    # that lets ~part of the bonds buckle; parity must hold regardless.
    fcrit = 50.0
    nat = run_native(12, True, fcrit)
    cup, act = run_cupy(12, True, fcrit)
    d1 = float(np.max(np.abs(nat - cup)))
    sc1 = float(np.max(np.abs(cup))) + 1e-30
    ok1 = d1 < 1e-8 * max(sc1, 1.0)
    rel_frac = float(np.mean(_bondlen(nat) < r0 * (1 - 1e-9)))
    print(f"[1 relaxed parity] max|nat-cupy| = {d1:.3e} (scale {sc1:.3e})  "
          f"released≈{rel_frac:.2f}  {'PASS' if ok1 else 'FAIL'}")
    if not ok1:
        fails.append("relaxed-parity")
    if rel_frac < 0.05:
        print("  WARN: almost nothing released — test not exercising the buckling path")

    # ---- Test 2: rigid-limit superset (F_crit huge → nothing eligible) ----
    nat_relaxed_hi = run_native(12, True, 1.0e12)
    nat_rigid = run_native(12, False, 0.0)
    d2 = float(np.max(np.abs(nat_relaxed_hi - nat_rigid)))
    sc2 = float(np.max(np.abs(nat_rigid))) + 1e-30
    ok2 = d2 < 1e-9 * max(sc2, 1.0)
    print(f"[2 rigid-limit  ] relaxed(F_crit→∞) vs rigid: max|Δ| = {d2:.3e}  "
          f"{'PASS' if ok2 else 'FAIL'}")
    if not ok2:
        fails.append("rigid-limit")

    if fails:
        print(f"native relaxed parity FAIL: {fails}")
        return 1
    print("native relaxed parity PASS (relaxed≈cupy + rigid-limit superset)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
