"""Measure Python callback vs native updater overhead in HOOMD run loops."""

from __future__ import annotations

import argparse
import json
import time

import hoomd

from ffn_hoomd_plugin import NativeNoOpUpdater, NativePositionKickUpdater


class PythonNoOpAction(hoomd.custom.Action):
    def __init__(self) -> None:
        self.count = 0

    def act(self, timestep: int) -> None:
        self.count += 1


def build_simulation(device: hoomd.device.Device, n_particles: int) -> hoomd.Simulation:
    sim = hoomd.Simulation(device=device, seed=23)
    snap = hoomd.Snapshot()
    if snap.communicator.rank == 0:
        snap.configuration.box = [100.0, 100.0, 100.0, 0.0, 0.0, 0.0]
        snap.particles.N = n_particles
        snap.particles.types = ["A"]
        snap.particles.position[:] = [[0.0, 0.0, 0.0]] * n_particles
    sim.create_state_from_snapshot(snap)
    return sim


def attach_mode(sim: hoomd.Simulation, mode: str):
    if mode == "none":
        return None
    if mode == "python-noop":
        action = PythonNoOpAction()
        updater = hoomd.update.CustomUpdater(
            action=action, trigger=hoomd.trigger.Periodic(1)
        )
        sim.operations.updaters.append(updater)
        return action
    if mode == "native-noop":
        updater = NativeNoOpUpdater(trigger=hoomd.trigger.Periodic(1))
        sim.operations.updaters.append(updater)
        return updater
    if mode == "native-position-kick":
        updater = NativePositionKickUpdater(
            delta=(0.0, 0.0, 0.0), trigger=hoomd.trigger.Periodic(1)
        )
        sim.operations.updaters.append(updater)
        return updater
    raise ValueError(f"Unknown mode: {mode}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["none", "python-noop", "native-noop", "native-position-kick"],
        required=True,
    )
    parser.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    parser.add_argument("--steps", type=int, default=200_000)
    parser.add_argument("--n-particles", type=int, default=1)
    args = parser.parse_args()

    device = hoomd.device.GPU() if args.device == "gpu" else hoomd.device.CPU()
    sim = build_simulation(device, args.n_particles)
    op = attach_mode(sim, args.mode)
    sim.run(0)

    start = time.perf_counter()
    sim.run(args.steps)
    elapsed_s = time.perf_counter() - start
    count = getattr(op, "update_count", getattr(op, "count", None))
    result = {
        "mode": args.mode,
        "device": args.device,
        "steps": args.steps,
        "n_particles": args.n_particles,
        "elapsed_s": elapsed_s,
        "steps_per_s": args.steps / elapsed_s,
        "count": count,
    }
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
