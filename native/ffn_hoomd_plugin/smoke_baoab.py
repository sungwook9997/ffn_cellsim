"""Native BAOAB (Stage 1a) smoke + parity vs the validated cupy device path.

Three checks (GPU):
  * drift      — constant force, kT=0: Δx == (F/γ)·t exactly (deterministic).
  * diffusion  — zero force, kT>0: D == kT/γ within tolerance (Einstein).
  * parity     — zero force, kT>0: native FFNBaoabUpdater vs the cupy
                 OverdampedBAOABDevice are BIT-IDENTICAL (same splitmix64 RNG
                 keyed by (seed, timestep, tag), same two-Gaussian step).

Run on gbook:
    cd ~/ffn_cellsim/native/ffn_hoomd_plugin
    PYTHONPATH="$PWD/build:$PWD/python:$HOME/ffn_cellsim" \
        python smoke_baoab.py --test all
"""

from __future__ import annotations

import argparse

import hoomd
import numpy as np

from ffn_hoomd_plugin import NativeBaoabUpdater


def _make_snapshot(N: int, L: float, positions: np.ndarray) -> hoomd.Snapshot:
    snap = hoomd.Snapshot()
    if snap.communicator.rank == 0:
        snap.configuration.box = [L, L, L, 0.0, 0.0, 0.0]
        snap.particles.N = N
        snap.particles.types = ["A"]
        snap.particles.position[:] = positions
        snap.particles.typeid[:] = np.zeros(N, dtype=np.int32)
    return snap


def _integrator(dt: float, forces: list) -> hoomd.md.Integrator:
    integ = hoomd.md.Integrator(dt=dt)
    integ.forces = list(forces)
    integ.methods = []  # BAOAB Updater owns the position step
    return integ


def _unwrapped(sim: hoomd.Simulation, L: float) -> np.ndarray:
    snap = sim.state.get_snapshot()
    pos = np.array(snap.particles.position, dtype=np.float64)
    img = np.array(snap.particles.image, dtype=np.float64)
    return pos + img * L


def _run_native(*, device, N, L, pos0, dt, kT, gamma, seed, steps, force=None):
    sim = hoomd.Simulation(device=device, seed=1)
    sim.create_state_from_snapshot(_make_snapshot(N, L, pos0))
    forces = []
    if force is not None:
        cf = hoomd.md.force.Constant(filter=hoomd.filter.All())
        cf.constant_force["A"] = tuple(force)
        cf.constant_torque["A"] = (0.0, 0.0, 0.0)
        forces = [cf]
    sim.operations.integrator = _integrator(dt, forces)
    upd = NativeBaoabUpdater(
        dt=dt, kT=kT, gamma_by_tag=[gamma] * N, seed=seed,
        trigger=hoomd.trigger.Periodic(1),
    )
    sim.operations.updaters.append(upd)
    sim.run(0)  # seed net_force F(r0)
    sim.run(steps)
    return _unwrapped(sim, L)


def _run_cupy(*, device, N, L, pos0, dt, kT, gamma, seed, steps):
    from ffn_sim.integrator.baoab_device import OverdampedBAOABDevice

    sim = hoomd.Simulation(device=device, seed=1)
    sim.create_state_from_snapshot(_make_snapshot(N, L, pos0))
    sim.operations.integrator = _integrator(dt, [])
    action = OverdampedBAOABDevice(kT=kT, gamma={"A": gamma}, dt=dt, seed=seed)
    sim.operations.updaters.append(
        hoomd.update.CustomUpdater(action=action, trigger=hoomd.trigger.Periodic(1))
    )
    sim.run(0)
    sim.run(steps)
    return _unwrapped(sim, L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    ap.add_argument(
        "--test", choices=["drift", "diffusion", "parity", "bench", "all"], default="all"
    )
    ap.add_argument("--n", type=int, default=4096)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--dt", type=float, default=1e-3)
    ap.add_argument("--kT", type=float, default=1.0)
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    dev = hoomd.device.GPU() if args.device == "gpu" else hoomd.device.CPU()
    L = 1.0e6  # huge box → wrapping never fires in these short tests
    rng = np.random.default_rng(12345)
    pos0 = (rng.random((args.n, 3)) - 0.5) * (L * 0.1)

    fails = []

    if args.test in ("drift", "all"):
        F = 2.0
        out = _run_native(
            device=dev, N=args.n, L=L, pos0=pos0.copy(), dt=args.dt, kT=0.0,
            gamma=args.gamma, seed=args.seed, steps=args.steps, force=(F, 0.0, 0.0),
        )
        dx = float(np.mean(out[:, 0] - pos0[:, 0]))
        expect = (F / args.gamma) * args.dt * args.steps
        rel = abs(dx - expect) / abs(expect)
        ok = rel < 1e-6
        print(f"[drift]     Δx={dx:.10g} expect={expect:.10g} rel={rel:.2e} "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append("drift")

    if args.test in ("diffusion", "all"):
        out = _run_native(
            device=dev, N=args.n, L=L, pos0=pos0.copy(), dt=args.dt, kT=args.kT,
            gamma=args.gamma, seed=args.seed, steps=args.steps, force=None,
        )
        disp = out - pos0
        msd = float(np.mean(np.sum(disp * disp, axis=1)))
        t = args.dt * args.steps
        D_meas = msd / (6.0 * t)
        D_exp = args.kT / args.gamma
        rel = abs(D_meas - D_exp) / D_exp
        ok = rel < 0.05
        print(f"[diffusion] D_meas={D_meas:.6g} D_exp={D_exp:.6g} rel={rel:.2e} "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append("diffusion")

    if args.test in ("parity", "all"):
        nat = _run_native(
            device=dev, N=args.n, L=L, pos0=pos0.copy(), dt=args.dt, kT=args.kT,
            gamma=args.gamma, seed=args.seed, steps=args.steps, force=None,
        )
        cup = _run_cupy(
            device=dev, N=args.n, L=L, pos0=pos0.copy(), dt=args.dt, kT=args.kT,
            gamma=args.gamma, seed=args.seed, steps=args.steps,
        )
        max_abs = float(np.max(np.abs(nat - cup)))
        scale = float(np.max(np.abs(cup))) + 1e-30
        ok = max_abs < 1e-9 * max(scale, 1.0)
        print(f"[parity]    max|native-cupy|={max_abs:.3e} scale={scale:.3e} "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append("parity")

    if args.test == "bench":
        import time

        from ffn_sim.integrator.baoab_device import OverdampedBAOABDevice

        def _build(use_native):
            sim = hoomd.Simulation(device=dev, seed=1)
            sim.create_state_from_snapshot(_make_snapshot(args.n, L, pos0.copy()))
            sim.operations.integrator = _integrator(args.dt, [])
            if use_native:
                sim.operations.updaters.append(
                    NativeBaoabUpdater(
                        dt=args.dt, kT=args.kT, gamma_by_tag=[args.gamma] * args.n,
                        seed=args.seed, trigger=hoomd.trigger.Periodic(1),
                    )
                )
            else:
                act = OverdampedBAOABDevice(
                    kT=args.kT, gamma={"A": args.gamma}, dt=args.dt, seed=args.seed
                )
                sim.operations.updaters.append(
                    hoomd.update.CustomUpdater(action=act, trigger=hoomd.trigger.Periodic(1))
                )
            sim.run(0)
            sim.run(200)  # warm autotuners / JIT before timing
            return sim

        res = {}
        for label, native in (("native", True), ("cupy", False)):
            sim = _build(native)
            t0 = time.perf_counter()
            sim.run(args.steps)
            wall = time.perf_counter() - t0
            res[label] = args.steps / wall
            print(f"[bench] {label:6s} N={args.n} steps={args.steps} "
                  f"{res[label]:.1f} steps/s ({1e6*wall/args.steps:.1f} us/step)")
        print(f"[bench] integrator-only speedup native/cupy = "
              f"{res['native'] / res['cupy']:.2f}x")
        return 0

    if fails:
        print(f"SMOKE FAIL: {fails}")
        return 1
    print("native-baoab smoke PASS (all)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
