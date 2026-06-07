"""End-to-end proof: pool binder (FFNAttachmentSpringForce) vs set_snapshot binder.

Same system, same bind/unbind schedule, a real integrator running — two arms:
  A  set_snapshot binder : bind/unbind = mutate HOOMD bonds via get/set_snapshot
                           (the production pattern; the 51.9 ms/firing wall)
  B  pool binder         : bind/unbind = set_attachments toggle on the native pool
                           (no snapshot, no bond mutation)
Measures end-to-end throughput so the ~8.4x binder lever is demonstrated in a loop,
not just a microbench. Read-only; new file.

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_binder_pool_demo --n 6000 --pool 3000 --steps 2000
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import hoomd
import hoomd.md as md


def _toggle_schedule(M, n_fire, frac, seed):
    """Precompute identical bind/unbind events for both arms: per firing, a
    boolean active-mask over the M slots (≈frac active, churned each firing)."""
    rng = np.random.default_rng(seed)
    masks = []
    active = rng.random(M) < frac
    for _ in range(n_fire):
        flip = rng.random(M) < 0.2 * frac  # churn ~20% of the active fraction
        active = np.where(flip, ~active, active)
        masks.append(active.copy())
    return masks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6000)
    ap.add_argument("--pool", type=int, default=3000)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--period", type=int, default=100)
    ap.add_argument("--seed", type=int, default=4)
    args = ap.parse_args()
    N, M = args.n, args.pool
    L, k, r0, kT, gamma = 1.0e-5, 2.0e-3, 1.0e-7, 4.0e-21, 1.0e-6
    dev = hoomd.device.GPU(notice_level=0)
    rng = np.random.default_rng(args.seed)
    pos = (rng.random((N, 3)) - 0.5) * (0.4 * L)
    perm = rng.permutation(N)[: 2 * M]
    head, actin = perm[:M].astype(np.int32), perm[M:].astype(np.int32)
    n_fire = args.steps // args.period + 2
    masks = _toggle_schedule(M, n_fire, 0.5, args.seed)

    def base_snap(with_bonds):
        s = hoomd.Snapshot()
        if s.communicator.rank == 0:
            s.configuration.box = [L, L, L, 0, 0, 0]
            s.particles.N = N; s.particles.types = ["A"]
            s.particles.position[:] = pos; s.particles.typeid[:] = np.zeros(N, np.int32)
            if with_bonds:
                s.bonds.types = ["spring"]  # start empty; binder adds
        return s

    # ---- ARM A: set_snapshot binder ----
    class SnapBinder(hoomd.custom.Action):
        def __init__(self):
            self.i = 0
        def act(self, timestep):
            sim = self._state._simulation
            m = masks[min(self.i, len(masks) - 1)]; self.i += 1
            grp = np.stack([head[m], actin[m]], axis=1).astype(np.int32)
            snap = sim.state.get_snapshot()            # global read
            if snap.communicator.rank == 0:
                snap.bonds.N = int(grp.shape[0])
                snap.bonds.types = ["spring"]
                if grp.shape[0]:
                    snap.bonds.group[:] = grp
                    snap.bonds.typeid[:] = np.zeros(grp.shape[0], np.int32)
            sim.state.set_snapshot(snap)               # global rebuild (the wall)

    simA = hoomd.Simulation(device=dev, seed=1)
    simA.create_state_from_snapshot(base_snap(True))
    hb = md.bond.Harmonic(); hb.params["spring"] = dict(k=k, r0=r0)
    igA = md.Integrator(dt=1e-10, forces=[hb], methods=[
        md.methods.Langevin(filter=hoomd.filter.All(), kT=kT)])
    igA.methods[0].gamma.default = gamma
    simA.operations.integrator = igA
    simA.operations.updaters.append(hoomd.update.CustomUpdater(
        action=SnapBinder(), trigger=hoomd.trigger.Periodic(args.period)))
    simA.run(args.period)  # warm (1 firing)
    t0 = time.perf_counter(); simA.run(args.steps); spsA = args.steps / (time.perf_counter() - t0)

    # ---- ARM B: pool binder ----
    from ffn_hoomd_plugin import NativeAttachmentSpringForce
    nat = NativeAttachmentSpringForce(pool_size=M)

    class PoolBinder(hoomd.custom.Action):
        def __init__(self):
            self.i = 0; self.r0v = np.full(M, r0)
        def act(self, timestep):
            m = masks[min(self.i, len(masks) - 1)]; self.i += 1
            kv = np.where(m, k, 0.0)                    # toggle k (no snapshot)
            nat.set_attachments(head, actin, kv, self.r0v)

    simB = hoomd.Simulation(device=dev, seed=1)
    simB.create_state_from_snapshot(base_snap(False))
    igB = md.Integrator(dt=1e-10, forces=[nat], methods=[
        md.methods.Langevin(filter=hoomd.filter.All(), kT=kT)])
    igB.methods[0].gamma.default = gamma
    simB.operations.integrator = igB
    simB.operations.updaters.append(hoomd.update.CustomUpdater(
        action=PoolBinder(), trigger=hoomd.trigger.Periodic(args.period)))
    simB.run(args.period)
    t0 = time.perf_counter(); simB.run(args.steps); spsB = args.steps / (time.perf_counter() - t0)

    print(f"binder pool demo: N={N} pool={M} period={args.period} steps={args.steps}", flush=True)
    print(f"  ARM A set_snapshot binder  {spsA:8.1f} steps/s ({1e6/spsA:.0f} us/step)", flush=True)
    print(f"  ARM B pool binder          {spsB:8.1f} steps/s ({1e6/spsB:.0f} us/step)", flush=True)
    print(f"  pool/set_snapshot end-to-end speedup = {spsB/spsA:.1f}x", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
