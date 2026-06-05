"""STAGE-2 diagnostics — crosslink coordination z (Ennomani contractility gate).

Literature (Ennomani 2016, nihms803347): actin-network contractility is maximal at
crosslink coordination z = (avg connectors per filament) in [2, 4], optimum ~3; this
= the rigidity-percolation threshold. Below 2 the network doesn't percolate (no global
contraction); above 4 it jams. (Belmonte 2017 corroborates: contractility maximal at
an intermediate crosslinker quantity; network must be elastically percolated.)

This measures the EMERGENT z of our cortex after the dynamic crosslinkers bind, so we
can confirm we're in [2,4] (the real contractility gate) before STAGE-2. z is computed
as 2 * (bound dimeric crosslinks) / n_filaments — each fully-bound xlink connects 2
filaments. We also report the single-head attach count.

Run:  conda run -n ffn_sim python -m ffn_sim.scripts.stage2_diagnostics --n-fil 1000
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers, xlink_attach_bin_names
from ffn_sim.cell.cell import build_cortex_full_simulation

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"


def measure_z(n_fil: int, n_xl: int | None, equil_steps: int, seed: int = 1,
              kon_scale: float = 1.0, bind_scale: float = 1.0):
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    if n_xl is not None:
        cfg["cortex"]["dynamic_crosslinkers"]["n_xl"] = int(n_xl)
    n_xl_eff = int(cfg["cortex"]["dynamic_crosslinkers"]["n_xl"])

    p = resolve_h3_derived(cfg)
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = 0.001 * tau_bend
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    if kon_scale != 1.0:  # accelerate binding (couple_accel-style; equilibrium probe)
        p_xl = replace(p_xl, k_on=p_xl.k_on * kon_scale)
    if bind_scale != 1.0:  # mesoscale-consistent reach: at ×40 areal coarse-graining the
        # inter-filament spacing grows ~sqrt(40)≈6.3x, so the physical 60nm filamin reach
        # cannot bridge the sparse mesoscale mesh. Scaling max_bind_dist by the spacing
        # ratio restores percolation at mesoscale WITHOUT stiffening the crosslinker (k
        # untouched — only the search radius for a binding partner). This is the geometric
        # dual of the ×40 filament coarse-graining already ratified (Plan v2 §3 H.3).
        p_xl = replace(p_xl, max_bind_dist=p_xl.max_bind_dist * bind_scale)
    dev = hoomd.device.CPU(notice_level=0)
    hw = build_cortex_full_simulation(
        p, p_xlinks=p_xl, device=dev, with_baoab=True, constrained=False,
        rng=np.random.default_rng(seed))
    sim = hw["sim"]
    sim.run(equil_steps)   # let the dynamic crosslinkers bind to steady state

    # Count xlink_attach bonds (one per bound head); a dimeric crosslink with both
    # heads bound = 2 attach bonds connecting 2 filaments.
    attach_names = set(xlink_attach_bin_names(p_xl.n_bins))
    with sim.state.cpu_local_snapshot as s:
        btypes = list(sim.state.bond_types)
        bt = np.asarray(s.bonds.typeid)
        attach_tids = [i for i, nm in enumerate(btypes) if nm in attach_names]
        n_attach = int(np.isin(bt, attach_tids).sum())

    # Each bound dimeric crosslink ~ 2 attach bonds -> connects 2 filaments.
    bound_dimers = n_attach / 2.0
    z = 2.0 * bound_dimers / n_fil          # avg filament-filament connections / filament
    return dict(n_fil=n_fil, n_xl=n_xl_eff, n_attach=n_attach,
                bound_frac=bound_dimers / max(n_xl_eff, 1), z=z)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--n-xl", type=int, default=None, help="override; default = config")
    ap.add_argument("--equil-steps", type=int, default=20000)
    ap.add_argument("--sweep", action="store_true",
                    help="sweep n_xl = {0.5,1,1.5,2}*n_fil to find z~3")
    ap.add_argument("--kon-scale", type=float, default=1.0,
                    help="scale crosslinker k_on (couple_accel-style binding accel)")
    ap.add_argument("--kon-sweep", action="store_true",
                    help="sweep kon_scale {1,30,100,300} at n_xl=n_fil to reach z in [2,4]")
    ap.add_argument("--bind-scale", type=float, default=1.0,
                    help="scale crosslinker max_bind_dist (mesoscale-consistent reach)")
    ap.add_argument("--bind-sweep", action="store_true",
                    help="sweep bind_scale {1,3,6,10} at fixed kon to find mesoscale percolation")
    args = ap.parse_args()

    print(f"=== STAGE-2 crosslink coordination z (gate [2,4], opt ~3; Ennomani) ===")
    print(f"HOOMD {hoomd.version.version}", flush=True)
    if args.bind_sweep:
        print(f"  (n_fil={args.n_fil}, kon_scale={args.kon_scale}x; scaling max_bind_dist)",
              flush=True)
        for bs in (1.0, 3.0, 6.0, 10.0):
            r = measure_z(args.n_fil, args.n_fil, args.equil_steps,
                          kon_scale=args.kon_scale, bind_scale=bs)
            gate = "IN-BAND" if 2 <= r["z"] <= 4 else ("LOW" if r["z"] < 2 else "JAMMED")
            print(f"  bind_scale={bs:5.1f}x (~{bs*60:.0f}nm reach)  n_attach={r['n_attach']:6d}  "
                  f"bound_frac={r['bound_frac']:.3f}  z={r['z']:.2f}  [{gate}]", flush=True)
    elif args.kon_sweep:
        for ks in (1.0, 30.0, 100.0, 300.0):
            r = measure_z(args.n_fil, args.n_fil, args.equil_steps, kon_scale=ks,
                          bind_scale=args.bind_scale)
            gate = "IN-BAND" if 2 <= r["z"] <= 4 else ("LOW" if r["z"] < 2 else "JAMMED")
            print(f"  kon_scale={ks:5.0f}x  n_attach={r['n_attach']:6d}  "
                  f"bound_frac={r['bound_frac']:.3f}  z={r['z']:.2f}  [{gate}]", flush=True)
    elif args.sweep:
        for mult in (0.5, 1.0, 1.5, 2.0):
            r = measure_z(args.n_fil, int(mult * args.n_fil), args.equil_steps,
                          kon_scale=args.kon_scale, bind_scale=args.bind_scale)
            gate = "IN-BAND" if 2 <= r["z"] <= 4 else ("LOW" if r["z"] < 2 else "JAMMED")
            print(f"  n_xl={r['n_xl']:5d} ({mult:.1f}x)  n_attach={r['n_attach']:5d}  "
                  f"bound_frac={r['bound_frac']:.2f}  z={r['z']:.2f}  [{gate}]", flush=True)
    else:
        r = measure_z(args.n_fil, args.n_xl, args.equil_steps, kon_scale=args.kon_scale,
                      bind_scale=args.bind_scale)
        gate = "IN-BAND" if 2 <= r["z"] <= 4 else ("LOW (sub-percolation)" if r["z"] < 2 else "JAMMED")
        print(f"  n_fil={r['n_fil']} n_xl={r['n_xl']} n_attach={r['n_attach']} "
              f"bound_frac={r['bound_frac']:.2f}", flush=True)
        print(f"  ==> z = {r['z']:.2f}  (gate [2,4], opt ~3)  [{gate}]", flush=True)
        if r["z"] < 2:
            need = int(np.ceil(3.0 / max(r["z"], 0.01) * r["n_xl"]))
            print(f"      to reach z~3: n_xl ~ {need} ({need/r['n_fil']:.1f}x n_fil)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
