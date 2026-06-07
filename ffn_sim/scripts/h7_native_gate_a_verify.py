"""Gate-A EVERYTHING-NATIVE verification: native integrator + native compartment forces.

Confirms the ~8.3x full-cell estimate end-to-end on the REAL production path
(build_baseline_cell, the Lead's two-phase warm->constrained harness), swapping
BOTH the constrained integrator (FFNConstrainedBaoabUpdater) AND the 3 compartment
md.force.Custom (-> FFNRadialShellForce) onto a built cell at runtime. Reuses the
Lead's GO driver helpers; adds the compartment-force swap + a 3rd arm.

Arms (same warmed config, same seed):
  A  cupy integrator + cpu compartment forces      (production baseline)
  B  native integrator + cpu compartment forces     (Lead's 2.43x GO arm)
  C  native integrator + NATIVE compartment forces  (everything-native; the 8.3x claim)
Reports pure-stepping throughput + ETA + the g_soft/g_rigid parity across arms
(g_soft includes the compartment-force tension, so C must match A/B).

  PYTHONPATH=".../build:.../python:$HOME/ffn_cellsim" \
    python -m ffn_sim.scripts.h7_native_gate_a_verify --n-fil 1000 --steps 2000
"""

from __future__ import annotations

import argparse
from copy import deepcopy

import numpy as np

from ffn_sim.cell.nucleus import NucleusConfinement
from ffn_sim.cell.membrane_surface import MembraneSurfaceTension
from ffn_sim.cortex.enclosed_volume import EnclosedVolumePressure
from ffn_sim.scripts.h7_native_fullcell_go import (
    _build_two_phase, _run_arm, _to_host, _NativeLambdaAdapter,
)
from ffn_sim.cell.manifest import load_manifest


def _swap_compartment_forces(cell):
    """Replace the 3 compartment Custom forces with native FFNRadialShellForce.
    Returns the count swapped."""
    from ffn_hoomd_plugin import NativeRadialShellForce
    integ = cell.simulation.operations.integrator
    swapped, kept = 0, []
    for f in list(integ.forces):
        rng = (getattr(f, "tag_start", None), getattr(f, "tag_end", None))
        if isinstance(f, NucleusConfinement) and None not in rng:
            kept.append(NativeRadialShellForce.from_nucleus(f.p, rng)); swapped += 1
        elif isinstance(f, MembraneSurfaceTension) and None not in rng:
            kept.append(NativeRadialShellForce.from_membrane(f.p, rng)); swapped += 1
        elif isinstance(f, EnclosedVolumePressure) and None not in rng:
            kept.append(NativeRadialShellForce.from_turgor(f.p, rng)); swapped += 1
        else:
            kept.append(f)
    integ.forces = kept
    return swapped


def _native_integrator(cell, P, hoomd):
    from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater
    sim = cell.simulation
    sim.operations.updaters.remove(cell.baoab_updater)
    nat = NativeConstrainedBaoabUpdater(
        dt=P["dt"], kT=P["kT"], inv_gamma_by_tag=P["inv_gamma"],
        chains_tag=P["chains_stacked"].reshape(-1).astype(np.int32),
        n_chains=P["F"], bonds_per_chain=P["m"], rest_length=P["r0"],
        seed=P["seed"], tol=1.0e-9, max_iter=200, trigger=hoomd.trigger.Periodic(1))
    sim.operations.updaters.append(nat)
    return nat


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--warmup", type=int, default=4000)
    ap.add_argument("--steps", type=int, default=2000)
    ap.add_argument("--interval", type=int, default=200)
    ap.add_argument("--gate-steps", type=float, default=2.0e8)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    import hoomd
    dev = hoomd.device.GPU(notice_level=0)

    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})["n_filaments"] = args.n_fil
    softstart = max(300, args.warmup // 8)
    build_constrained = _build_two_phase(manifest, device=dev, warmup=args.warmup,
                                         softstart=softstart, seed=args.seed)

    # ARM A: baseline (cupy integrator + cpu compartment) — extract params.
    cellA = build_constrained()
    R_cell = float(cellA.p_cortex.R_cell)
    act = cellA.baoab_action
    P = dict(dt=float(act.dt), kT=float(act.kT),
             inv_gamma=[float(x) for x in np.asarray(_to_host(act._inv_gamma_by_tag), np.float64)],
             chains_stacked=np.asarray(act.chains_tag_stacked),
             F=int(np.asarray(act.chains_tag_stacked).shape[0]),
             m=int(np.asarray(act.chains_tag_stacked).shape[1]) - 1,
             r0=float(act._chain_rest_length), seed=int(act._seed))
    sA, spsA = _run_arm(cellA, steps=args.steps, interval=args.interval,
                        R_cell=R_cell, dtc=P["dt"], rigid_action=act)
    ncomp = sum(isinstance(f, (NucleusConfinement, MembraneSurfaceTension, EnclosedVolumePressure))
                for f in cellA.simulation.operations.integrator.forces)
    N = int(cellA.simulation.state.N_particles)
    del cellA

    # ARM B: native integrator + cpu compartment.
    cellB = build_constrained()
    natB = _native_integrator(cellB, P, hoomd)
    cellB.simulation.run(0)
    sB, spsB = _run_arm(cellB, steps=args.steps, interval=args.interval, R_cell=R_cell,
                        dtc=P["dt"], rigid_action=_NativeLambdaAdapter(natB, P["chains_stacked"], P["r0"]))
    ncB = int(natB.nonconverged_count)
    del cellB

    # ARM C: native integrator + NATIVE compartment forces (everything-native).
    cellC = build_constrained()
    natC = _native_integrator(cellC, P, hoomd)
    nswap = _swap_compartment_forces(cellC)
    cellC.simulation.run(0)
    sC, spsC = _run_arm(cellC, steps=args.steps, interval=args.interval, R_cell=R_cell,
                        dtc=P["dt"], rigid_action=_NativeLambdaAdapter(natC, P["chains_stacked"], P["r0"]))
    ncC = int(natC.nonconverged_count)
    del cellC

    def gmean(s):
        a = np.asarray(s) * 1e3  # mN/m
        return (float(a[:, 0].mean()), float(a[:, 1].mean())) if len(a) else (float("nan"),) * 2

    gA, gB, gC = gmean(sA), gmean(sB), gmean(sC)
    eta = lambda sps: args.gate_steps / sps / 3600.0
    print(f"Gate-A verify: N={N}  compartment-forces={ncomp} (swapped to native: {nswap})  "
          f"dt={P['dt']:.3e}", flush=True)
    print(f"  ARM A cupy-int  + cpu-comp     {spsA:7.1f} steps/s  ETA {eta(spsA):6.1f} h  "
          f"g_soft={gA[0]:.3f} g_rigid={gA[1]:.3f} mN/m", flush=True)
    print(f"  ARM B nat-int   + cpu-comp     {spsB:7.1f} steps/s  ETA {eta(spsB):6.1f} h  "
          f"g_soft={gB[0]:.3f} g_rigid={gB[1]:.3f}  ({spsB/spsA:.2f}x A)  nonconv={ncB}", flush=True)
    print(f"  ARM C nat-int   + NATIVE-comp  {spsC:7.1f} steps/s  ETA {eta(spsC):6.1f} h  "
          f"g_soft={gC[0]:.3f} g_rigid={gC[1]:.3f}  ({spsC/spsA:.2f}x A)  nonconv={ncC}", flush=True)
    dg = abs(gC[0] - gA[0]) / (abs(gA[0]) + 1e-30)
    print(f"  g_soft parity C-vs-A: rel {dg:.2%} ({'PASS' if dg < 0.15 else 'CHECK'})", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
