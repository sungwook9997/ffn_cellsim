"""Diagnostic for the 'remove the nucleus for cortical-tension' route (PI 2026-06-08).

The nucleus stiffness contributes 0.0% to cortical γ across the KU-3.B2.1 band
(h7_nucleus_contribution_check), and the enclosed-volume turgor uses the CORTEX
shell volume (V0), NOT nucleus-subtracted — so the nucleus is decoupled from the
cortical observable at the suspended/rounded operating point. This script:

  (1) VERIFIES it directly — cortical γ with nucleus ON vs fully OFF must match.
  (2) Finds the dt CEILING with the nucleus removed — the constrained dt is set by
      a hardcoded 0.001·τ_bend (very conservative; CFL norm is ~0.1·τ), with the
      nucleus the first wall just above. Removing it, how far can dt go before the
      next compartment CFL (enclosed-volume / ERM) or numerical instability?

Suspended/rounded MCF7 (FA off, turgor on), MYOSIN OFF (passive structural γ),
constrained cupy, CPU. No-nucleus is realised by p_nucleus=None (the build skips
all nucleus blocks — cell.py:1061); the manifest baseline guard is bypassed ONLY
for this measurement-specific diagnostic.

    python -m ffn_sim.scripts.h7_no_nucleus_dt --device cpu --allow-cpu-dev
"""
from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

import numpy as np

import ffn_sim.cell.manifest as M
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


def _set_nucleus(enabled: bool):
    """Toggle the nucleus globally: enabled→real resolve; disabled→p_nucleus=None
    (build skips it). Manifest keeps nucleus 'enabled' so the baseline guard passes;
    only the resolved object is nulled — the sanctioned measurement-specific path."""
    M.resolve_nucleus = _orig_resolve_nucleus if enabled else (lambda *a, **k: None)


def _manifest(n_fil):
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["optional_subsystems"]["fa"]["enabled"] = False
    co = m.setdefault("cortex_overrides", {}).setdefault("cortex", {})
    co["n_filaments"] = int(n_fil)
    co["demo_mode"] = True
    co.setdefault("myosin", {})["n_motors_per_cell"] = 0
    return m


def _measure_gamma(*, nucleus, n_fil, warmup, softstart, ticks, interval, measure_last,
                   device, seed):
    _set_nucleus(nucleus)
    manifest = _manifest(n_fil)
    cw = build_baseline_cell(manifest=deepcopy(manifest), device=device, seed=seed,
                             constrained=False, with_baoab=True, equilibrate=True,
                             equilibrate_steps=warmup, equilibrate_softstart_steps=softstart,
                             connected_mesh=True)
    cw.simulation.run(0)
    warm = _tagpos(cw.simulation); del cw
    cell = build_baseline_cell(manifest=deepcopy(manifest), device=device, seed=seed,
                               constrained=True, equilibrate=False, connected_mesh=True)
    act = cell.baoab_action
    if hasattr(act, "record_lambda"):
        act.record_lambda = True
    sim = cell.simulation
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = warm
    sim.state.set_snapshot(snap); sim.run(0)
    R_cell = float(cell.p_cortex.R_cell); dtc = float(act.dt)
    pev = getattr(cell, "p_enclosed_volume", None)
    chans = {k: [] for k in ("gamma_soft", "gamma_rigid", "gamma_structural", "gamma_passive")}
    for k in range(ticks):
        sim.run(interval)
        if k >= ticks - measure_last:
            ct = measure_cortical_tension(sim, R_cell=R_cell, p_enclosed_volume=pev,
                                          lambda_accumulator=act, dt=dtc, n_planes=12)
            for key in chans:
                chans[key].append(ct[key] * _MN)
    out = {k: float(np.mean(v)) for k, v in chans.items()}
    out["n_part"] = int(sim.state.N_particles)
    del cell
    return out


def _dt_ceiling(*, n_fil, warmup, softstart, device, seed, stab_steps):
    """No-nucleus: raise constrained_dt_safety until a compartment CFL or numerical
    instability breaks. Returns the max stable multiple + what bound it."""
    _set_nucleus(False)
    manifest = _manifest(n_fil)
    cw = build_baseline_cell(manifest=deepcopy(manifest), device=device, seed=seed,
                             constrained=False, with_baoab=True, equilibrate=True,
                             equilibrate_steps=warmup, equilibrate_softstart_steps=softstart,
                             connected_mesh=True)
    cw.simulation.run(0); warm = _tagpos(cw.simulation); del cw
    dt0 = None
    results = []
    for mult in [1, 2, 5, 10, 20, 30, 50, 100]:
        try:
            cell = build_baseline_cell(manifest=deepcopy(manifest), device=device, seed=seed,
                                       constrained=True, equilibrate=False, connected_mesh=True,
                                       constrained_dt_safety=float(mult))
            if dt0 is None:
                dt0 = float(cell.baoab_action.dt) / mult
            sim = cell.simulation
            snap = sim.state.get_snapshot()
            if snap.communicator.rank == 0:
                snap.particles.position[:] = warm
            sim.state.set_snapshot(snap); sim.run(0)
            sim.run(stab_steps)
            p = _tagpos(sim)
            finite = bool(np.all(np.isfinite(p)))
            drift = float(getattr(cell.baoab_action, "max_constraint_drift", np.nan))
            ok = finite and (not np.isnan(drift)) and drift < 1e-6
            results.append((mult, "OK" if ok else f"UNSTABLE(finite={finite},drift={drift:.1e})"))
            print(f"    {mult:>3}x (dt={dt0*mult:.2e}): {'STABLE' if ok else 'UNSTABLE'} "
                  f"(drift={drift:.2e} finite={finite})", flush=True)
            del cell
            if not ok:
                break
        except Exception as e:
            msg = str(e).split(chr(10))[0][:120]
            results.append((mult, f"CFL_FAIL: {msg}"))
            print(f"    {mult:>3}x (dt={(dt0 or 0)*mult:.2e}): CFL_FAIL → {msg}", flush=True)
            break
    return dt0, results


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=160)
    ap.add_argument("--warmup", type=int, default=600)
    ap.add_argument("--softstart", type=int, default=150)
    ap.add_argument("--ticks", type=int, default=6)
    ap.add_argument("--interval", type=int, default=2000)
    ap.add_argument("--measure-last", type=int, default=3)
    ap.add_argument("--stab-steps", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1)
    add_production_device_args(ap, default="cpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)
    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    print("=" * 72, flush=True)
    print("NO-NUCLEUS ROUTE: (1) cortical-γ ON vs OFF  (2) dt ceiling without nucleus",
          flush=True)
    print("-" * 72, flush=True)
    kw = dict(n_fil=args.n_fil, warmup=args.warmup, softstart=args.softstart, ticks=args.ticks,
              interval=args.interval, measure_last=args.measure_last, device=dev, seed=args.seed)
    on = _measure_gamma(nucleus=True, **kw)
    off = _measure_gamma(nucleus=False, **kw)
    print(f"  nucleus ON : γ_struct={on['gamma_structural']:.4e}  γ_passive={on['gamma_passive']:.3e}  "
          f"n_part={on['n_part']}", flush=True)
    print(f"  nucleus OFF: γ_struct={off['gamma_structural']:.4e}  γ_passive={off['gamma_passive']:.3e}  "
          f"n_part={off['n_part']}", flush=True)
    dg = abs(on["gamma_structural"] - off["gamma_structural"]) / (abs(on["gamma_structural"]) + 1e-30)
    dp = abs(on["gamma_passive"] - off["gamma_passive"]) / (abs(on["gamma_passive"]) + 1e-30)
    print(f"  Δγ_struct = {dg*100:.2f}%   Δγ_passive(turgor) = {dp*100:.2f}%", flush=True)
    safe = dg < 0.02 and dp < 0.02
    print(f"  ⇒ removing the nucleus is {'SAFE (γ + turgor unchanged)' if safe else 'NOT safe — it shifts the observable'} "
          "for cortical tension", flush=True)
    print("-" * 72, flush=True)
    print("  dt CEILING without the nucleus (raise safety until CFL/instability):", flush=True)
    dt0, results = _dt_ceiling(n_fil=args.n_fil, warmup=args.warmup, softstart=args.softstart,
                               device=dev, seed=args.seed, stab_steps=args.stab_steps)
    max_stable = max([m for m, s in results if s in ("OK",)] + [m for m, s in results if "STABLE" in str(s)]
                     + [0]) if results else 0
    # last STABLE multiple
    stable_mults = [m for m, s in results if s == "OK"]
    max_stable = max(stable_mults) if stable_mults else 1
    print("-" * 72, flush=True)
    print(f"  → max STABLE dt multiple without nucleus = {max_stable}x  (dt {dt0*max_stable:.2e} s)",
          flush=True)
    print(f"  → step-count reduction vs current = ~{max_stable}x  (+ fewer particles: "
          f"{on['n_part']}→{off['n_part']})", flush=True)
    _OUT.mkdir(parents=True, exist_ok=True)
    (_OUT / "h7_no_nucleus_dt.json").write_text(json.dumps(
        {"gamma_on": on, "gamma_off": off, "d_gamma_struct_frac": dg, "d_turgor_frac": dp,
         "removal_safe": safe, "dt0": dt0, "dt_sweep": results, "max_stable_mult": max_stable},
        indent=2, default=str))
    print(f"  json → {_OUT / 'h7_no_nucleus_dt.json'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
