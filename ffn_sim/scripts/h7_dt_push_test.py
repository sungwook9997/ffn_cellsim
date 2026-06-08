"""Empirical dt-push test (no nucleus): how far past the conservative CFL guards
can the constrained dt go while keeping cortical γ ACCURATE (not just stable)?

The compartment CFL gates use a 0.1 safety on τ = γ_b/k_eff — that's an ACCURACY
margin, well below the overdamped stability limit (~2τ). The membrane gate flags
10× as a violation, but the membrane mode may still be both stable AND accurate
there. This bypasses the guards (cfl_strict=False) and measures, at each dt
multiple, the constraint drift + whether γ_structural stays consistent with the 1×
baseline (same warm state + seed, so dt is the only variable). The usable dt is the
largest multiple where γ holds AND drift stays tiny.

Nucleus OFF (PI 2026-06-08 route), suspended/rounded, myo OFF, constrained cupy, CPU.

    python -m ffn_sim.scripts.h7_dt_push_test --device cpu --allow-cpu-dev
"""
from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import numpy as np

import ffn_sim.cell.manifest as M
import ffn_sim.cell.cell as C
from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)
from ffn_sim.cortex.cortical_tension import measure_cortical_tension
from ffn_sim.scripts.h7_native_fullcell_go import _tagpos

_OUT = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "production"
_MN = 1.0e3
_orig_resolve_nucleus = M.resolve_nucleus


def _bypass_cfl_and_nucleus():
    """No nucleus (p_nucleus=None) + bypass the conservative compartment CFL guards
    so we can probe the true (accuracy-limited) dt ceiling empirically."""
    M.resolve_nucleus = lambda *a, **k: None

    def _nostrict(fn):
        def w(*a, **k):
            k["cfl_strict"] = False
            return fn(*a, **k)
        return w
    for nm in ("attach_enclosed_volume_to_simulation", "attach_erm_to_simulation",
               "attach_membrane_surface", "attach_nucleus_confinement"):
        if hasattr(C, nm):
            setattr(C, nm, _nostrict(getattr(C, nm)))


def _manifest(n_fil):
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["optional_subsystems"]["fa"]["enabled"] = False
    co = m.setdefault("cortex_overrides", {}).setdefault("cortex", {})
    co["n_filaments"] = int(n_fil); co["demo_mode"] = True
    co.setdefault("myosin", {})["n_motors_per_cell"] = 0
    return m


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=160)
    ap.add_argument("--warmup", type=int, default=800)
    ap.add_argument("--softstart", type=int, default=200)
    ap.add_argument("--settle", type=int, default=12000, help="constrained steps to settle/measure")
    ap.add_argument("--equal-phys-time", action="store_true",
                    help="run each dt for the SAME physical time (steps=settle/mult) so the dt "
                         "discretization error is isolated from extra settling")
    ap.add_argument("--mults", type=float, nargs="+", default=[1, 5, 10, 15, 20])
    ap.add_argument("--seed", type=int, default=1)
    add_production_device_args(ap, default="cpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)
    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    _bypass_cfl_and_nucleus()
    manifest = _manifest(args.n_fil)

    # One shared warm state (no nucleus) so dt is the ONLY variable across conditions.
    cw = build_baseline_cell(manifest=deepcopy(manifest), device=dev, seed=args.seed,
                             constrained=False, with_baoab=True, equilibrate=True,
                             equilibrate_steps=args.warmup, equilibrate_softstart_steps=args.softstart,
                             connected_mesh=True)
    cw.simulation.run(0); warm = _tagpos(cw.simulation); del cw

    print("=" * 72, flush=True)
    print("dt-PUSH (no nucleus, CFL guards bypassed) — stability + γ accuracy", flush=True)
    print(f"  n_fil={args.n_fil} settle={args.settle} steps  (membrane CFL guard flags ~5.8x)",
          flush=True)
    print("-" * 72, flush=True)
    rows = []
    dt0 = None
    g_ref = None
    for mult in args.mults:
        try:
            cell = build_baseline_cell(manifest=deepcopy(manifest), device=dev, seed=args.seed,
                                       constrained=True, equilibrate=False, connected_mesh=True,
                                       constrained_dt_safety=float(mult))
        except Exception as e:
            print(f"  {mult:>4.0f}x: BUILD FAIL → {str(e).splitlines()[0][:90]}", flush=True)
            break
        act = cell.baoab_action
        if hasattr(act, "record_lambda"):
            act.record_lambda = True
        if dt0 is None:
            dt0 = float(act.dt) / mult
        sim = cell.simulation
        snap = sim.state.get_snapshot()
        if snap.communicator.rank == 0:
            snap.particles.position[:] = warm
        sim.state.set_snapshot(snap)
        R_cell = float(cell.p_cortex.R_cell); dtc = float(act.dt)
        pev = getattr(cell, "p_enclosed_volume", None)
        steps = int(round(args.settle / mult)) if args.equal_phys_time else args.settle
        try:
            sim.run(0); sim.run(steps)
            ct = measure_cortical_tension(sim, R_cell=R_cell, p_enclosed_volume=pev,
                                          lambda_accumulator=act, dt=dtc, n_planes=12)
            gstruct = ct["gamma_structural"] * _MN
            drift = float(getattr(act, "max_constraint_drift", np.nan))
            p = _tagpos(sim); finite = bool(np.all(np.isfinite(p)))
            if g_ref is None:
                g_ref = gstruct
            dg = abs(gstruct - g_ref) / (abs(g_ref) + 1e-30)
            ok = finite and drift < 1e-6 and dg < 0.15
            print(f"  {mult:>4.0f}x (dt={dtc:.2e}, {steps}st): γ_struct={gstruct:.4e} "
                  f"(Δvs1x={dg*100:5.1f}%) drift={drift:.1e} finite={finite}  "
                  f"{'OK' if ok else 'BREAK'}", flush=True)
            rows.append({"mult": mult, "dt": dtc, "gamma_struct": gstruct, "dgamma_frac": dg,
                         "drift": drift, "finite": finite, "ok": ok})
            del cell
            if not ok:
                break
        except Exception as e:
            print(f"  {mult:>4.0f}x (dt={dtc:.2e}): RUN FAIL → {str(e).splitlines()[0][:90]}",
                  flush=True)
            rows.append({"mult": mult, "dt": dtc, "ok": False, "error": str(e).splitlines()[0][:120]})
            break
    okmults = [r["mult"] for r in rows if r.get("ok")]
    best = max(okmults) if okmults else 1
    print("-" * 72, flush=True)
    print(f"  → usable dt multiple (stable AND γ-accurate, no nucleus) = {best:.0f}x", flush=True)
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "h7_dt_push_test.json").write_text(json.dumps({"dt0": dt0, "rows": rows, "best": best},
                                                          indent=2, default=str))
    print(f"  json → {_OUT / 'h7_dt_push_test.json'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
