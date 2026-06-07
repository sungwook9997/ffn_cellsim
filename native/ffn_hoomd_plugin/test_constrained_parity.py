"""End-to-end native constrained L-M BAOAB (Stage 1b assembly) validation.

Test A (kT=0 + constant force): native FFNConstrainedBaoabUpdater vs the cupy
  ConstrainedLeimkuhlerMatthewsBAOAB. At kT=0 the noise prefactor and the Fixman
  force both vanish (∝ kT), so the step is deterministic — predictor + M-SHAKE +
  wrap + tag→row — and must agree bit-for-bit (no RNG-stream mismatch). Chains
  start STRETCHED so M-SHAKE is exercised. (Fixman + the splitmix64 noise are
  validated separately: fixman/shake/baoab parity tests.)
Test B (kT>0, native alone): the full step (Fixman + noise + M-SHAKE) runs stable
  and keeps the rigid bonds satisfied over many steps.

    cd ~/ffn_cellsim/native/ffn_hoomd_plugin
    PYTHONPATH="$PWD/build:$PWD/python:$HOME/ffn_cellsim" python test_constrained_parity.py
"""

from __future__ import annotations

import numpy as np
import hoomd

from ffn_sim.integrator.constrained_baoab import ConstrainedLeimkuhlerMatthewsBAOAB
from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater


def _build(F, m, L, r0, stretch, gamma, rng):
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
    F, m = 1000, 6
    L, r0, gamma, dt = 1.0e6, 100.0, 1.0, 1.0e-3
    seed, tol, mit = 7, 1.0e-10, 200
    dev = hoomd.device.GPU()
    fails = []

    # ---- Test A: kT=0, constant force, stretched init → deterministic parity ----
    rng = np.random.default_rng(11)
    N, pos, chains, cpairs = _build(F, m, L, r0, 1.08, gamma, rng)
    cl = np.full(cpairs.shape[0], r0)
    chains_flat = np.concatenate([c for c in chains]).astype(np.int32)

    def run_native(kT, steps):
        sim = hoomd.Simulation(device=dev, seed=1)
        sim.create_state_from_snapshot(_snapshot(N, L, pos))
        cf = hoomd.md.force.Constant(filter=hoomd.filter.All())
        cf.constant_force["A"] = (1.0, 0.0, 0.0)
        cf.constant_torque["A"] = (0.0, 0.0, 0.0)
        ig = hoomd.md.Integrator(dt=dt); ig.forces = [cf]; ig.methods = []
        sim.operations.integrator = ig
        sim.operations.updaters.append(NativeConstrainedBaoabUpdater(
            dt=dt, kT=kT, inv_gamma_by_tag=[1.0 / gamma] * N, chains_tag=chains_flat,
            n_chains=F, bonds_per_chain=m, rest_length=r0, seed=seed, tol=tol, max_iter=mit,
            trigger=hoomd.trigger.Periodic(1)))
        sim.run(0); sim.run(steps)
        return _unwrap(sim, L)

    def run_cupy(kT, steps):
        sim = hoomd.Simulation(device=dev, seed=1)
        sim.create_state_from_snapshot(_snapshot(N, L, pos))
        cf = hoomd.md.force.Constant(filter=hoomd.filter.All())
        cf.constant_force["A"] = (1.0, 0.0, 0.0)
        cf.constant_torque["A"] = (0.0, 0.0, 0.0)
        ig = hoomd.md.Integrator(dt=dt); ig.forces = [cf]; ig.methods = []
        sim.operations.integrator = ig
        act = ConstrainedLeimkuhlerMatthewsBAOAB(
            kT=kT, gamma={"A": gamma}, dt=dt, constraint_pairs=cpairs,
            constraint_lengths=cl, chains=chains, seed=seed, shake_tol=tol, shake_max_iter=mit)
        sim.operations.updaters.append(
            hoomd.update.CustomUpdater(action=act, trigger=hoomd.trigger.Periodic(1)))
        sim.run(0); sim.run(steps)
        return _unwrap(sim, L)

    nat = run_native(0.0, 10)
    cup = run_cupy(0.0, 10)
    dA = float(np.max(np.abs(nat - cup)))
    scA = float(np.max(np.abs(cup))) + 1e-30
    okA = dA < 1e-9 * max(scA, 1.0)
    print(f"[A kT=0 parity] max|nat-cupy| = {dA:.3e} (scale {scA:.3e})  "
          f"{'PASS' if okA else 'FAIL'}")
    if not okA:
        fails.append("A")

    # ---- Test B: kT>0 native alone → constraint stays satisfied + finite ----
    out = run_native(1.0, 500)
    s = out[chains_flat.reshape(F, m + 1)[:, :-1]] - out[chains_flat.reshape(F, m + 1)[:, 1:]]
    bond = np.sqrt(np.sum(s * s, axis=2))
    drift = float(np.max(np.abs(bond - r0)) / r0)
    finite = bool(np.all(np.isfinite(out)))
    okB = drift < 1e-7 and finite
    print(f"[B kT>0 native ] 500 steps  max bond drift = {drift:.3e}  finite={finite}  "
          f"{'PASS' if okB else 'FAIL'}")
    if not okB:
        fails.append("B")

    if fails:
        print(f"constrained assembly FAIL: {fails}")
        return 1
    print("native-constrained assembly PASS (A parity + B stability)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
