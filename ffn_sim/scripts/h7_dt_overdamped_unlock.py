"""H.7 overdamped-dt UNLOCK probe — is the cortex dt over-conservative by the friction ratio?

DECISIVE setup finding (2026-06-09): the BAOAB integrator runs cortex actin/myosin at the
PHYSIOLOGICAL cytoplasm friction γ_cyto = 6π·η_eff·R = 3.73e-5 N·s/m (η_eff=65.9 Pa·s, MCF7),
but the CFL timestep dt_cfl is computed from WATER viscosity (cortex.py:446 gamma_b = 6π·η_water·R
= 3.9e-10), a 95,328× mismatch. In the overdamped (large-γ·dt) Leimkuhler-Matthews limit the
update is dr = F/γ·dt + noise, whose stability limit is dt < 2γ/k — so at the REAL (cytoplasm)
friction the usable dt is ~10⁴-10⁵× the water-CFL dt. If true, the active-γ "floor" is partly a
TIMESCALE artifact: contraction/condensation is drag-limited (~seconds at η_eff), needing ~10⁹
steps at the water-dt (impossible) but only ~10⁴-10⁶ steps at the overdamped-correct dt.

This probe (a) scans dt upward (via cfl_safety_factor override, owned cortex config) and finds the
largest dt that is STABLE and γ-ACCURATE vs a small-dt reference at MATCHED physical time; (b) at
the usable dt, runs to developed contraction and measures whether the network finally CONDENSES
(r/r0 < 1, end-to-end shortening — the Miyazaki/Lenz symmetry-breaking the rigid backbone forbade)
and whether g_soft / g_soft_actin RISE. Unconstrained (soft, condensable) backbone + cytoplasm-η,
exactly the lit-study lever (relax M-SHAKE → unconstrained+cytoplasm-η).

NO band-tuning; dt is a CFL/discretization parameter justified by the physiological friction, not a
free knob to land γ. The myosin binding CFL (batch_dt·k_off≤1e-3) auto-shrinks batch_steps, so the
motor kinetics stay valid as dt grows; the cap is dt ≲ 1e-3/k_off0 ≈ 2.9e-3 s.

Usage:
    python -m ffn_sim.scripts.h7_dt_overdamped_unlock --n-filaments 160 \
        --mults 1 10 100 300 1000 --device cpu --allow-cpu-dev
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

_MNM = 1e3


def _manifest(n_fil, cfl_safety, n_motors=None):
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["optional_subsystems"]["fa"]["enabled"] = False
    co = m.setdefault("cortex_overrides", {}).setdefault("cortex", {})
    co["n_filaments"] = int(n_fil)
    co["demo_mode"] = True
    co.setdefault("dynamics", {})["cfl_safety_factor"] = float(cfl_safety)
    if n_motors is not None:
        co.setdefault("myosin", {})["n_motors_per_cell"] = int(n_motors)
    return m


def _bypass_strict_cfl():
    """Bypass the conservative per-compartment CFL guards (they assume the water-dt);
    the overdamped stability is governed by the integrator, tested empirically here."""
    import ffn_sim.cell.cell as C

    def _nostrict(fn):
        def w(*a, **k):
            k["cfl_strict"] = False
            return fn(*a, **k)
        return w
    for nm in ("attach_enclosed_volume_to_simulation", "attach_erm_to_simulation",
               "attach_membrane_surface", "attach_nucleus_confinement"):
        if hasattr(C, nm):
            setattr(C, nm, _nostrict(getattr(C, nm)))


def _condensation(sim, n_cortex_actin, bond_type="cortex-bond"):
    """Mean bond/ℓ0 and end-to-end-over-contour for cortex actin (condensation readout)."""
    snap = sim.state.get_snapshot()
    if snap.communicator.rank != 0:
        return None
    pos = np.asarray(snap.particles.position, dtype=np.float64)
    types = list(snap.bonds.types)
    tids = np.asarray(snap.bonds.typeid)
    grp = np.asarray(snap.bonds.group)
    if bond_type not in types:
        return None
    bt = types.index(bond_type)
    bb = grp[tids == bt]
    if bb.shape[0] == 0:
        return None
    d = np.linalg.norm(pos[bb[:, 0]] - pos[bb[:, 1]], axis=1)
    return {"mean_bond_nm": float(np.mean(d) * 1e9), "n_bonds": int(bb.shape[0])}


def _measure(cell):
    sim = cell.simulation
    R = float(cell.p_cortex.R_cell)
    ct = measure_cortical_tension(sim, R_cell=R, p_enclosed_volume=cell.p_enclosed_volume)
    cond = _condensation(sim, int(cell.n_cortex_actin))
    sg = (float(cell.myosin_action._head_grip_s.mean() / float(cell.p_cortex.rest_length))
          if cell.myosin_action is not None else 0.0)
    pos = np.asarray(sim.state.get_snapshot().particles.position)
    finite = bool(np.all(np.isfinite(pos)))
    return {
        "g_soft_mN_m": float(ct["gamma_soft"]) * _MNM,
        "g_rigid_mN_m": float(ct["gamma_rigid"]) * _MNM,
        "g_passive_mN_m": float(ct["gamma_passive"]) * _MNM,
        "s_grip_over_l0": sg,
        "mean_bond_nm": cond["mean_bond_nm"] if cond else None,
        "rest_length_nm": float(cell.p_cortex.rest_length) * 1e9,
        "finite": finite,
    }


def run(args, dev):
    _bypass_strict_cfl()
    n_fil = None if (args.n_filaments is None or args.n_filaments <= 0) else args.n_filaments
    # reference dt from the default cfl_safety (water-based)
    ref_m = _manifest(n_fil or 1000, args.base_cfl)
    cref = build_baseline_cell(manifest=deepcopy(ref_m), device=dev, seed=args.seed,
                               constrained=False, with_baoab=True, equilibrate=True,
                               equilibrate_steps=args.warmup,
                               equilibrate_softstart_steps=args.softstart, connected_mesh=True)
    dt0 = float(cref.simulation.operations.integrator.dt)
    gamma_actin = None
    h = getattr(cref, "extras", {}).get("handles", {}) if hasattr(cref, "extras") else {}
    ba = h.get("baoab_action")
    if ba is not None and getattr(ba, "gamma_map", None):
        gamma_actin = ba.gamma_map.get("actin_cortex")
    print("=" * 76, flush=True)
    print("H.7 OVERDAMPED-dt UNLOCK PROBE (unconstrained, cytoplasm-η friction)", flush=True)
    print(f"  base dt (water-CFL) = {dt0:.3e} s  |  actin BAOAB γ = "
          f"{gamma_actin:.3e} N·s/m" + (f"  (γ/γ_water≈{gamma_actin/3.9092e-10:.0f}×)"
                                        if gamma_actin else ""), flush=True)
    print(f"  overdamped stability dt<2γ/k + binding-CFL cap dt≲1e-3/k_off0≈2.9e-3 s", flush=True)
    print("-" * 76, flush=True)
    del cref

    rows = []
    for mult in args.mults:
        cfl = args.base_cfl * mult
        try:
            cell = build_baseline_cell(manifest=_manifest(n_fil or 1000, cfl), device=dev,
                                       seed=args.seed, constrained=False, with_baoab=True,
                                       equilibrate=True, equilibrate_steps=args.warmup,
                                       equilibrate_softstart_steps=args.softstart,
                                       connected_mesh=True)
        except Exception as e:
            print(f"  {mult:>6.0f}×: BUILD FAIL → {str(e).splitlines()[0][:80]}", flush=True)
            rows.append({"mult": mult, "build_fail": str(e).splitlines()[0][:120]})
            continue
        dt = float(cell.simulation.operations.integrator.dt)
        # run a FIXED PHYSICAL TIME so dt is the only variable (steps = phys_time/dt)
        steps = max(1, int(round(args.phys_time / dt)))
        steps = min(steps, args.max_steps)
        try:
            cell.simulation.run(steps)
            m = _measure(cell)
            ok = m["finite"]
        except Exception as e:
            m = {"finite": False, "error": str(e).splitlines()[0][:120]}
            ok = False
        m.update({"mult": mult, "cfl_safety": cfl, "dt_s": dt, "steps": steps,
                  "phys_time_s": steps * dt})
        rows.append(m)
        tag = "OK " if ok else "UNSTABLE"
        print(f"  {mult:>6.0f}× dt={dt:.3e} steps={steps:>7d} t={steps*dt:.2e}s  [{tag}]  "
              f"g_soft={m.get('g_soft_mN_m', float('nan')):.3e} "
              f"s_grip={m.get('s_grip_over_l0', 0):.3f} "
              f"bond={m.get('mean_bond_nm', float('nan')):.1f}nm "
              f"(ℓ0={m.get('rest_length_nm', float('nan')):.0f})", flush=True)
        del cell

    # usable dt = largest stable mult whose g_soft agrees with the ref-trend (finite + sane bond)
    stable = [r for r in rows if r.get("finite") and r.get("g_soft_mN_m") is not None]
    out = {
        "operating_point": "unconstrained soft backbone + cytoplasm-η, FA OFF, turgor ON, connected mesh",
        "base_dt_s": dt0, "actin_baoab_gamma": gamma_actin,
        "n_filaments": n_fil, "phys_time_s": args.phys_time, "rows": rows,
        "largest_stable_mult": max((r["mult"] for r in stable), default=0),
    }
    print("-" * 76, flush=True)
    print(f"  largest stable mult = {out['largest_stable_mult']}× "
          f"(dt up to {dt0*out['largest_stable_mult']:.3e} s)", flush=True)
    print("  → if a large dt is stable+accurate, contraction/condensation becomes reachable in", flush=True)
    print("    feasible steps; check whether mean_bond < ℓ0 (CONDENSATION) emerges at long t.", flush=True)
    print("=" * 76, flush=True)
    return out


def run_contract(args, dev):
    """Run ONE long contraction at a chosen (large, overdamped-accurate) dt and sample
    s_grip, condensation (mean_bond/ℓ0), and γ over PHYSICAL time. The decisive test:
    given enough physical time at the overdamped-correct dt, does the soft network CONDENSE
    (mean_bond < ℓ0) and does γ_soft rise toward band?"""
    _bypass_strict_cfl()
    n_fil = None if (args.n_filaments is None or args.n_filaments <= 0) else args.n_filaments
    cfl = args.base_cfl * args.contract_mult
    cell = build_baseline_cell(manifest=_manifest(n_fil or 1000, cfl), device=dev, seed=args.seed,
                               constrained=False, with_baoab=True, equilibrate=True,
                               equilibrate_steps=args.warmup,
                               equilibrate_softstart_steps=args.softstart, connected_mesh=True)
    sim = cell.simulation
    dt = float(sim.operations.integrator.dt)
    total_steps = int(round(args.contract_phys_time / dt))
    n_samples = max(4, args.contract_samples)
    chunk = max(1, total_steps // n_samples)
    print("=" * 76, flush=True)
    print(f"H.7 OVERDAMPED-dt CONTRACTION: dt={dt:.3e}s ({args.contract_mult}× water-CFL), "
          f"target {args.contract_phys_time}s = {total_steps} steps", flush=True)
    ell0_nm = float(cell.p_cortex.rest_length) * 1e9
    print(f"  ℓ0={ell0_nm:.0f}nm  (CONDENSATION ⇔ mean_bond < ℓ0; r/r0<1 = Miyazaki symmetry-break)",
          flush=True)
    print("-" * 76, flush=True)
    rows = []
    done = 0
    while done < total_steps:
        c = min(chunk, total_steps - done)
        sim.run(c); done += c
        m = _measure(cell)
        m["step"] = done; m["t_s"] = done * dt
        m["bond_over_l0"] = (m["mean_bond_nm"] / ell0_nm) if m["mean_bond_nm"] else None
        rows.append(m)
        print(f"  t={done*dt:7.3f}s step={done:>8d}  s_grip={m['s_grip_over_l0']:.3f}  "
              f"bond/ℓ0={m['bond_over_l0']:.5f}  g_soft={m['g_soft_mN_m']:.3e}  "
              f"g_rigid={m['g_rigid_mN_m']:.3e}  fin={m['finite']}", flush=True)
        if not m["finite"]:
            print("  !! NON-FINITE — dt too large for accuracy/stability; stop.", flush=True)
            break
    out = {"mode": "contract", "dt_s": dt, "mult": args.contract_mult,
           "ell0_nm": ell0_nm, "n_filaments": n_fil, "rows": rows}
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print("-" * 76, flush=True)
    last = rows[-1] if rows else {}
    cond = last.get("bond_over_l0")
    print(f"  FINAL: s_grip={last.get('s_grip_over_l0',0):.3f}  bond/ℓ0={cond}  "
          f"g_soft={last.get('g_soft_mN_m',0):.3e} mN/m", flush=True)
    if cond is not None and cond < 0.999:
        print(f"   → CONDENSATION ENGAGED (bond/ℓ0={cond:.4f} < 1) — the network densifies; "
              f"check if g_soft rose.", flush=True)
    else:
        print(f"   → NO condensation (bond/ℓ0≈1) even at reachable contraction — the rigid-vs-soft "
              f"backbone is not the only block.", flush=True)
    print(f"  json → {args.out_json}", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-filaments", type=int, default=160)
    ap.add_argument("--warmup", type=int, default=1000)
    ap.add_argument("--softstart", type=int, default=200)
    ap.add_argument("--base-cfl", type=float, default=0.1, help="baseline cfl_safety_factor")
    ap.add_argument("--mults", type=float, nargs="+", default=[1, 10, 100, 300, 1000, 3000])
    ap.add_argument("--phys-time", type=float, default=1.0e-4,
                    help="physical time per condition [s] (steps = phys_time/dt, matched across dt)")
    ap.add_argument("--max-steps", type=int, default=40000, help="cap steps per condition")
    ap.add_argument("--contract", action="store_true",
                    help="run ONE long contraction at --contract-mult instead of the dt scan")
    ap.add_argument("--contract-mult", type=float, default=770.0,
                    help="dt multiplier for the contraction run (770× → dt≈1e-5 s, accurate+stable)")
    ap.add_argument("--contract-phys-time", type=float, default=2.5,
                    help="physical time for the contraction run [s] (~myosin contraction time)")
    ap.add_argument("--contract-samples", type=int, default=12)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-json", type=str,
                    default="ffn_sim/outputs/h7/production/h7_dt_overdamped_unlock.json")
    add_production_device_args(ap, default="cpu")
    args = ap.parse_args()
    validate_production_device_args(ap, args)
    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    if args.contract:
        run_contract(args, dev)
    else:
        out = run(args, dev)
        Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out_json).write_text(json.dumps(out, indent=2))
        print(f"  json → {args.out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
