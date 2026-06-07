"""Smoke-test that the native updater attaches and runs inside HOOMD."""

from __future__ import annotations

import argparse

import hoomd

from ffn_hoomd_plugin import NativeNoOpUpdater


def build_simulation(device: hoomd.device.Device) -> hoomd.Simulation:
    sim = hoomd.Simulation(device=device, seed=13)
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
    parser.add_argument("--steps", type=int, default=16)
    args = parser.parse_args()

    device = hoomd.device.GPU() if args.device == "gpu" else hoomd.device.CPU()
    sim = build_simulation(device)
    updater = NativeNoOpUpdater(trigger=hoomd.trigger.Periodic(1))
    sim.operations.updaters.append(updater)
    sim.run(args.steps)

    if updater.update_count != args.steps:
        raise SystemExit(
            f"Native updater count mismatch: {updater.update_count} != {args.steps}"
        )
    print(
        "native-noop-updater smoke PASS "
        f"device={args.device} steps={args.steps} "
        f"count={updater.update_count} last={updater.last_timestep}"
    )


if __name__ == "__main__":
    main()
