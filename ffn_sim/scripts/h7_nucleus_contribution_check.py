"""Diagnostic: how much does the NUCLEUS STIFFNESS contribute to cortical tension?

The nucleus's stiff lamin mode sets the constrained-integrator CFL (the dt wall),
but does that stiffness actually transmit to the cortical γ we measure in Gate-B?
If γ_structural is ~flat across the KU-3.B2.1 nuclear-modulus band [1, 10] kPa, the
nucleus is "numerically loud, mechanically quiet" for cortical tension — justifying
a multi-timescale integrator (keep the nucleus for volume/confinement, sub-cycle
its fast mode without paying the CFL tax).

Suspended/rounded MCF7 (FA off, turgor on; Gate-B §6 operating point), MYOSIN OFF
(isolates the PASSIVE structural cortical tension), constrained cupy integrator,
CPU/smoke. Sweeps E_nuc over the band; reports γ_soft / γ_rigid / γ_structural /
γ_passive for each.

    python -m ffn_sim.scripts.h7_nucleus_contribution_check --device cpu --allow-cpu-dev
"""
from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.cortex.cortical_tension import measure_cortical_tension
from ffn_sim.scripts.h7_native_fullcell_go import _tagpos

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "production"
_MN = 1.0e3


def _manifest(*, n_fil, e_nuc):
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["optional_subsystems"]["fa"]["enabled"] = False
    co = m.setdefault("cortex_overrides", {}).setdefault("cortex", {})
    co["n_filaments"] = int(n_fil)
    co["demo_mode"] = True
    co.setdefault("myosin", {})["n_motors_per_cell"] = 0  # myo OFF → passive structure
    m["compartments"]["nucleus"]["E_nuc"] = float(e_nuc)  # the swept variable
    return m


def _run(*, e_nuc, n_fil, warmup, softstart, ticks, interval, measure_last, device, seed,
         dt_safety):
    manifest = _manifest(n_fil=n_fil, e_nuc=e_nuc)
    cell_w = build_baseline_cell(
        manifest=deepcopy(manifest), device=device, seed=seed, constrained=False,
        with_baoab=True, equilibrate=True, equilibrate_steps=warmup,
        equilibrate_softstart_steps=softstart, connected_mesh=True)
    cell_w.simulation.run(0)
    warm_pos = _tagpos(cell_w.simulation)
    del cell_w

    # Reduced constrained dt so the STIFF end of the band (10 kPa) is also CFL-stable
    # (the nucleus CFL ∝ 1/E_nuc); γ at equilibrium is dt-independent, so the
    # across-band comparison stays valid with all conditions at the same safe dt.
    cell = build_baseline_cell(
        manifest=deepcopy(manifest), device=device, seed=seed, constrained=True,
        equilibrate=False, connected_mesh=True, constrained_dt_safety=dt_safety)
    act = cell.baoab_action
    if hasattr(act, "record_lambda"):
        act.record_lambda = True
    sim = cell.simulation
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = warm_pos
    sim.state.set_snapshot(snap)
    sim.run(0)

    R_cell = float(cell.p_cortex.R_cell)
    dtc = float(act.dt)
    pev = getattr(cell, "p_enclosed_volume", None)
    chans = {k: [] for k in ("gamma_soft", "gamma_rigid", "gamma_structural", "gamma_passive")}
    for k in range(ticks):
        sim.run(interval)
        if k >= ticks - measure_last:
            ct = measure_cortical_tension(sim, R_cell=R_cell, p_enclosed_volume=pev,
                                          lambda_accumulator=act, dt=dtc, n_planes=12)
            for key in chans:
                chans[key].append(ct[key] * _MN)
    out = {k: float(np.mean(v)) if v else 0.0 for k, v in chans.items()}
    out["E_nuc_kPa"] = e_nuc / 1e3
    out["tau_nuc_us"] = None
    try:
        # τ_nuc the CFL sees, for context
        out["tau_nuc_us"] = float(getattr(cell, "_nuc_tau", np.nan)) * 1e6
    except Exception:
        pass
    del cell
    print(f"  E_nuc={e_nuc/1e3:.1f}kPa  γ_struct={out['gamma_structural']:.4e}  "
          f"(soft={out['gamma_soft']:.3e} rigid={out['gamma_rigid']:.3e}) "
          f"γ_passive={out['gamma_passive']:.3e} mN/m", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=160)
    ap.add_argument("--warmup", type=int, default=600)
    ap.add_argument("--softstart", type=int, default=150)
    ap.add_argument("--ticks", type=int, default=6)
    ap.add_argument("--interval", type=int, default=2000)
    ap.add_argument("--measure-last", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--e-nuc-kpa", type=float, nargs="+", default=[1.0, 4.7, 10.0],
                    help="nuclear moduli to sweep (KU-3.B2.1 band [1,10] kPa)")
    ap.add_argument("--dt-safety", type=float, default=0.4,
                    help="constrained_dt_safety so the stiff (10 kPa) end is CFL-stable")
    add_production_device_args(ap, default="cpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)
    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    print("=" * 72, flush=True)
    print("NUCLEUS-STIFFNESS → cortical-γ contribution (suspended, myo OFF, passive)",
          flush=True)
    print(f"  n_fil={args.n_fil}  band sweep E_nuc(kPa)={args.e_nuc_kpa}", flush=True)
    print("-" * 72, flush=True)
    rows = []
    for ek in args.e_nuc_kpa:
        rows.append(_run(e_nuc=ek * 1e3, n_fil=args.n_fil, warmup=args.warmup,
                         softstart=args.softstart, ticks=args.ticks, interval=args.interval,
                         measure_last=args.measure_last, device=dev, seed=args.seed,
                         dt_safety=args.dt_safety))
    print("-" * 72, flush=True)
    gs = [r["gamma_structural"] for r in rows]
    span = (max(gs) - min(gs))
    rel = span / (np.mean(gs) + 1e-30)
    print(f"  γ_structural across band: min={min(gs):.4e} max={max(gs):.4e} mN/m", flush=True)
    print(f"  → spread = {rel*100:.1f}% of mean (×{max(gs)/(min(gs)+1e-30):.2f})", flush=True)
    if rel < 0.15:
        print("  ⇒ cortical γ is ~INSENSITIVE to nucleus stiffness across the band —", flush=True)
        print("    nucleus is mechanically QUIET for cortical tension (numerically loud only).", flush=True)
        print("    → multi-timescale / nucleus-subcycle integrator is PHYSICALLY JUSTIFIED.", flush=True)
    else:
        print("  ⇒ cortical γ DEPENDS on nucleus stiffness — nucleus is mechanically", flush=True)
        print("    coupled to the cortex; sub-cycling would change the observable.", flush=True)
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "h7_nucleus_contribution.json").write_text(json.dumps(
        {"rows": rows, "gamma_struct_spread_frac": rel}, indent=2))
    print(f"  json → {_OUT / 'h7_nucleus_contribution.json'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
