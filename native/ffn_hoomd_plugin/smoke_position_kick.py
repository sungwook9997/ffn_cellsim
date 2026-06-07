"""Smoke-test native position mutation on CPU and GPU HOOMD devices."""

from __future__ import annotations

import argparse
import math

import hoomd

from ffn_hoomd_plugin import NativePositionKickUpdater


def build_simulation(device: hoomd.device.Device) -> hoomd.Simulation:
    sim = hoomd.Simulation(device=device, seed=17)
    snap = hoomd.Snapshot()
    if snap.communicator.rank == 0:
        snap.configuration.box = [10.0, 10.0, 10.0, 0.0, 0.0, 0.0]
        snap.particles.N = 1
        snap.particles.types = ["A"]
        snap.particles.position[:] = [[0.0, 0.0, 0.0]]
    sim.create_state_from_snapshot(snap)
    return sim


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    parser.add_argument("--steps", type=int, default=8)
    parser.add_argument("--dx", type=float, default=0.125)
    args = parser.parse_args()

    device = hoomd.device.GPU() if args.device == "gpu" else hoomd.device.CPU()
    sim = build_simulation(device)
    updater = NativePositionKickUpdater(
        delta=(args.dx, 0.0, 0.0), trigger=hoomd.trigger.Periodic(1)
    )
    sim.operations.updaters.append(updater)
    sim.run(args.steps)

    snap = sim.state.get_snapshot()
    x_final = float(snap.particles.position[0, 0])
    expected = args.steps * args.dx
    if not math.isclose(x_final, expected, rel_tol=1e-7, abs_tol=1e-7):
        raise SystemExit(
            f"Position mismatch: x_final={x_final:.12g}, expected={expected:.12g}"
        )
    if updater.update_count != args.steps:
        raise SystemExit(
            f"Native updater count mismatch: {updater.update_count} != {args.steps}"
        )

    print(
        "native-position-kick smoke PASS "
        f"device={args.device} steps={args.steps} "
        f"x_final={x_final:.8g} count={updater.update_count}"
    )


if __name__ == "__main__":
    main()
