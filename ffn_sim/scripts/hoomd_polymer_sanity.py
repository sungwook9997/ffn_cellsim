"""Phase 0.2 — HOOMD-blue install verification (Plan Unit H.0.2).

A 100-bead linear polymer chain with harmonic bonds and harmonic angles
under a Langevin thermostat. Confirms the v2 environment (HOOMD 7.0.1,
osx-arm64 CPU mode, Python 3.13) is functional end-to-end:

  * `hoomd.Simulation` constructs and steps
  * `hoomd.md.bond.Harmonic` and `hoomd.md.angle.Harmonic` evaluate
  * `hoomd.md.methods.Langevin` thermostats correctly
  * `gsd.hoomd` writes trajectories the rest of the toolchain can read

This is a *smoke test*, not a science run. Output goes to
ffn_sim/outputs/h_0_2/.

Run:
    conda activate ffn_sim
    python ffn_sim/scripts/hoomd_polymer_sanity.py --steps 5000

For Plan-spec 1M-step bench:
    python ffn_sim/scripts/hoomd_polymer_sanity.py --steps 1000000 --bench
"""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import gsd.hoomd
import hoomd
import numpy as np

N_BEADS = 100
BOND_LENGTH = 1.0
BOND_K = 100.0
ANGLE_T0 = math.pi
ANGLE_K = 5.0
KT = 1.0
GAMMA = 1.0
DT = 0.005
LOG_EVERY = 1000

OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h_0_2"


def build_initial_snapshot() -> gsd.hoomd.Frame:
    snap = gsd.hoomd.Frame()
    snap.particles.N = N_BEADS
    snap.particles.types = ["A"]
    snap.particles.typeid = np.zeros(N_BEADS, dtype=np.uint32)
    xs = np.arange(N_BEADS, dtype=np.float64) * BOND_LENGTH - (N_BEADS - 1) * BOND_LENGTH / 2.0
    pos = np.column_stack([xs, np.zeros(N_BEADS), np.zeros(N_BEADS)])
    snap.particles.position = pos
    snap.particles.mass = np.ones(N_BEADS, dtype=np.float64)

    snap.bonds.N = N_BEADS - 1
    snap.bonds.types = ["polymer"]
    snap.bonds.typeid = np.zeros(N_BEADS - 1, dtype=np.uint32)
    snap.bonds.group = np.column_stack([np.arange(N_BEADS - 1), np.arange(1, N_BEADS)]).astype(np.uint32)

    snap.angles.N = N_BEADS - 2
    snap.angles.types = ["bend"]
    snap.angles.typeid = np.zeros(N_BEADS - 2, dtype=np.uint32)
    snap.angles.group = np.column_stack([
        np.arange(N_BEADS - 2),
        np.arange(1, N_BEADS - 1),
        np.arange(2, N_BEADS),
    ]).astype(np.uint32)

    box_side = max(2.0 * N_BEADS * BOND_LENGTH, 50.0)
    snap.configuration.box = [box_side, box_side, box_side, 0.0, 0.0, 0.0]
    return snap


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--bench", action="store_true", help="print wall-time per 1k steps at the end")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    initial_path = OUT_DIR / "initial.gsd"
    trajectory_path = OUT_DIR / "trajectory.gsd"

    with gsd.hoomd.open(initial_path, mode="w") as f:
        f.append(build_initial_snapshot())

    device = hoomd.device.CPU()
    sim = hoomd.Simulation(device=device, seed=args.seed)
    sim.create_state_from_gsd(filename=str(initial_path))

    harmonic_bond = hoomd.md.bond.Harmonic()
    harmonic_bond.params["polymer"] = dict(k=BOND_K, r0=BOND_LENGTH)
    harmonic_angle = hoomd.md.angle.Harmonic()
    harmonic_angle.params["bend"] = dict(k=ANGLE_K, t0=ANGLE_T0)

    integrator = hoomd.md.Integrator(dt=DT)
    integrator.forces.append(harmonic_bond)
    integrator.forces.append(harmonic_angle)
    langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=KT)
    langevin.gamma["A"] = GAMMA
    integrator.methods.append(langevin)
    sim.operations.integrator = integrator

    writer = hoomd.write.GSD(
        trigger=hoomd.trigger.Periodic(LOG_EVERY),
        filename=str(trajectory_path),
        mode="wb",
        filter=hoomd.filter.All(),
        dynamic=["property"],
    )
    sim.operations.writers.append(writer)

    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    print(f"Phase 0.2 polymer sanity — {N_BEADS} beads, {args.steps} steps, dt={DT}, kT={KT}")
    print(f"Device: {device} (CPU)  HOOMD {hoomd.version.version}")

    sim.run(0)
    kE0 = float(thermo.kinetic_energy)
    pE0 = float(thermo.potential_energy)
    print(f"  step 0      | KE={kE0:.3f}  PE={pE0:.3f}")

    t0 = time.perf_counter()
    sim.run(args.steps)
    elapsed = time.perf_counter() - t0

    kE = float(thermo.kinetic_energy)
    pE = float(thermo.potential_energy)
    T_measured = (2.0 / 3.0) * kE / N_BEADS
    print(f"  step {args.steps:>6} | KE={kE:.3f}  PE={pE:.3f}  T_eff={T_measured:.3f} (target {KT:.3f})")
    print(f"Wall time: {elapsed:.2f} s for {args.steps} steps ({args.steps / elapsed:.0f} steps/s)")

    if args.bench:
        per_1k = elapsed / args.steps * 1000.0
        print(f"Bench: {per_1k:.3f} s / 1k steps")

    print(f"\nArtifacts:\n  - initial:    {initial_path}\n  - trajectory: {trajectory_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
