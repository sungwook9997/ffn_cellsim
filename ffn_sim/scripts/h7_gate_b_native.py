"""H.7 GATE-B production run on the full ×40 native GPU stack — does permitting
filament buckling (relaxed M-SHAKE, compression side) transmit Gate-A's contraction
into spanning cortical tension?

Contract: docs/v2_audit/H7_GATE_B_CONTRACT_2026-06-08.md (LOCKED). Mechanism: the
native unilateral M-SHAKE + τ_bend-EMA Euler gate (validated bit-parity vs cupy,
test_constrained_relaxed_parity.py). This is the authoritative full-scale run the
smoke driver (h7_gate_b_relaxed_smoke.py, cupy) defers to.

Observable (§4): γ_active = [γ_soft + γ_rigid](myo ON) − [...](myo OFF), suspended/
rounded cell (FA off, turgor on, §6), at the settled operating point where s_grip
has developed (as in Gate-A). One condition per invocation so the four
(rigid/relaxed × myoON/OFF) run as checkpointed jobs on the single GPU; then
``--aggregate`` computes γ_active + the verdict.

Pass/fail (§7, on γ_active mN/m): CONFIRM [0.18,0.40] / PARTIAL ≥0.03 / REFUTE <0.03.

Usage (gbook, per condition):
    export PYTHONPATH="$PWD/native/ffn_hoomd_plugin/build:$PWD/native/ffn_hoomd_plugin/python:$PYTHONPATH"
    export FFN_GPU_DEVICE_COMPARTMENTS=1
    python -u -m ffn_sim.scripts.h7_gate_b_native --mode relaxed --myo on \
        --n-fil 1000 --warmup 4000 --ticks 600 --interval 5000 --tag gateB
Aggregate once all four exist:
    python -m ffn_sim.scripts.h7_gate_b_native --aggregate --tag gateB
"""
from __future__ import annotations

import argparse
import json
import os
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.scripts.h3_ku35_tension import (
    _tension_method_of_planes,
    _tension_method_of_planes_rigid,
)
from ffn_sim.scripts.h7_native_fullcell_go import (
    _NativeLambdaAdapter,
    _tagpos,
    _to_host,
)

_OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "production"
_FIG_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "figs"
CONFIRM_BAND = (0.18, 0.40)   # Hosseini MCF7 interphase IQR (contract §7)
PARTIAL_FLOOR = 0.03
GATE_A_PLATEAU_MN_M = 3.06e-3
_CONDS = ("rigid_myoON", "rigid_myoOFF", "relaxed_myoON", "relaxed_myoOFF")


def _cond_name(mode, myo_on):
    return f"{mode}_myo{'ON' if myo_on else 'OFF'}"


def _mean_s_grip(ma, ell0):
    if ma is None:
        return 0.0, 0.0
    bound = np.asarray(_to_host(ma._head_bound_to_actin)) >= 0
    if not np.any(bound):
        return 0.0, 0.0
    s = float(np.mean(np.asarray(_to_host(ma._head_grip_s))[bound]))
    return s, s / ell0


def _make_manifest(*, n_fil, myo_on, no_nucleus=False):
    """Suspended/rounded MCF7 (FA off, turgor on; §6), myosin per condition.

    ``no_nucleus`` (PI 2026-06-08): drop the nucleus for the cortical-tension run —
    it is decoupled from cortical γ (0% across the band; removal Δγ within seed
    noise) and its stiff-lamin CFL is the dt bottleneck. Removing it lets dt rise to
    the membrane CFL (~5×). MUST NOT be used for whole-cell/confinement observables."""
    m = deepcopy(load_manifest("mcf7_baseline.yaml"))
    m["optional_subsystems"]["fa"]["enabled"] = False
    co = m.setdefault("cortex_overrides", {}).setdefault("cortex", {})
    co["n_filaments"] = int(n_fil)
    if not myo_on:
        co.setdefault("myosin", {})["n_motors_per_cell"] = 0  # §5.4 baseline
    if no_nucleus:
        m["compartments"]["nucleus"]["enabled"] = False
    return m


def _build_native(manifest, *, device, seed, warm_pos, relaxed, no_nucleus=False, dt_safety=1.0,
                  connected_mesh=True):
    """Constrained full cell + native integrator swap; relaxed adds the unilateral
    M-SHAKE + τ_bend-EMA Euler gate (F_crit, τ_bend DERIVED from κ_B, ℓ₀).
    ``connected_mesh`` (default TRUE) builds the percolated spanning cortex (the
    2026-06-04 rebuild: bridge-different-filament crosslinks, giant ≥0.9) — REQUIRED
    for a transmission measurement; False is the fragmented mesh (giant ~7%, the
    γ-floor artifact). ``no_nucleus`` + ``dt_safety`` are the optional optimizations."""
    import hoomd
    from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater

    cell = build_baseline_cell(
        manifest=deepcopy(manifest), device=device, seed=seed,
        constrained=True, equilibrate=False, connected_mesh=connected_mesh,
        allow_no_nucleus=no_nucleus, constrained_dt_safety=dt_safety)
    if hasattr(cell.baoab_action, "record_lambda"):
        cell.baoab_action.record_lambda = True
    sim = cell.simulation
    snap = sim.state.get_snapshot()
    if snap.communicator.rank == 0:
        snap.particles.position[:] = warm_pos
    sim.state.set_snapshot(snap)
    sim.run(0)

    act = cell.baoab_action
    dtc = float(act.dt)
    chains_stacked = np.asarray(act.chains_tag_stacked)
    F, npc = chains_stacked.shape
    inv_gamma = [float(x) for x in np.asarray(
        _to_host(act._inv_gamma_by_tag), dtype=np.float64)]
    r0 = float(act._chain_rest_length)
    kappa_B = float(cell.p_cortex.bending_modulus)
    F_crit = float(np.pi ** 2 * kappa_B / r0 ** 2)         # Euler buckling (§8)
    tau_bend = float(cell.p_cortex.tau_bend)               # EMA window (PI 2026-06-08)

    sim.operations.updaters.remove(cell.baoab_updater)
    nat = NativeConstrainedBaoabUpdater(
        dt=dtc, kT=float(act.kT), inv_gamma_by_tag=inv_gamma,
        chains_tag=chains_stacked.reshape(-1).astype(np.int32),
        n_chains=F, bonds_per_chain=npc - 1, rest_length=r0,
        seed=int(act._seed), tol=1.0e-9, max_iter=200,
        trigger=hoomd.trigger.Periodic(1),
        compression_release=bool(relaxed),
        release_load_crit=(F_crit if relaxed else 0.0),
        load_tau=(tau_bend if relaxed else 0.0))
    sim.operations.updaters.append(nat)
    sim.run(0)
    adapter = _NativeLambdaAdapter(nat, chains_stacked, r0)
    info = {"F_crit_pN": F_crit * 1e12, "tau_bend_us": tau_bend * 1e6,
            "ell0_nm": r0 * 1e9, "n_cortex_actin": int(F * npc),
            "no_nucleus": bool(no_nucleus), "dt_safety": float(dt_safety),
            "connected_mesh": bool(connected_mesh),
            "dt_s": float(dtc), "n_part": int(cell.simulation.state.N_particles)}
    return cell, nat, adapter, dtc, r0, info


def _condensation(sim, chains_stacked, ell0):
    pos = _tagpos(sim)
    cf = np.asarray(chains_stacked)
    seg = pos[cf[:, :-1]] - pos[cf[:, 1:]]
    blen = np.linalg.norm(seg, axis=2)
    e2e = np.linalg.norm(pos[cf[:, -1]] - pos[cf[:, 0]], axis=1)
    contour = (cf.shape[1] - 1) * ell0
    return {"mean_bond_over_l0": float(np.mean(blen) / ell0),
            "frac_bonds_compressed": float(np.mean(blen < ell0)),
            "mean_e2e_over_contour": float(np.mean(e2e) / contour)}


def _save_ckpt(path, *, tick, positions, myo_s, myo_bound, rows):
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, tick=tick, positions=positions,
             myo_s=myo_s if myo_s is not None else np.empty(0),
             myo_bound=myo_bound if myo_bound is not None else np.empty(0),
             rows=json.dumps(rows))
    os.replace(tmp, path)


def _run_condition(args, dev):
    import hoomd
    mode = args.mode
    myo_on = (args.myo == "on")
    name = _cond_name(mode, myo_on)
    relaxed = (mode == "relaxed")
    manifest = _make_manifest(n_fil=args.n_fil, myo_on=myo_on, no_nucleus=args.no_nucleus)
    softstart = max(300, args.warmup // 8)
    ckpt = _OUT_DIR / f"h7_{args.tag}_{name}.ckpt.npz"
    cond_json = _OUT_DIR / f"h7_{args.tag}_{name}.cond.json"
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    _FIG_DIR.mkdir(parents=True, exist_ok=True)

    rows, start_tick, warm_pos = [], 0, None
    ckpt_myo_s = ckpt_myo_bound = None
    if ckpt.exists():
        z = np.load(ckpt, allow_pickle=True)
        warm_pos = z["positions"]; start_tick = int(z["tick"]); rows = json.loads(str(z["rows"]))
        ckpt_myo_s = z["myo_s"] if z["myo_s"].size else None
        ckpt_myo_bound = z["myo_bound"] if z["myo_bound"].size else None
        print(f"[RESUME] {name} tick {start_tick}; {len(rows)} rows", flush=True)
    else:
        cell_w = build_baseline_cell(
            manifest=deepcopy(manifest), device=dev, seed=args.seed,
            constrained=False, equilibrate=True, equilibrate_steps=args.warmup,
            equilibrate_softstart_steps=softstart, allow_no_nucleus=args.no_nucleus,
            connected_mesh=(not args.no_connected_mesh))
        cell_w.simulation.run(0)
        warm_pos = _tagpos(cell_w.simulation)
        del cell_w

    cell, nat, adapter, dtc, ell0, info = _build_native(
        manifest, device=dev, seed=args.seed, warm_pos=warm_pos, relaxed=relaxed,
        no_nucleus=args.no_nucleus, dt_safety=args.dt_safety,
        connected_mesh=(not args.no_connected_mesh))
    sim = cell.simulation
    R_cell = float(cell.p_cortex.R_cell)
    if ckpt_myo_s is not None and cell.myosin_action is not None:
        try:
            cell.myosin_action._head_grip_s[:] = ckpt_myo_s
            cell.myosin_action._head_bound_to_actin[:] = ckpt_myo_bound
        except Exception as exc:  # noqa: BLE001
            print(f"[RESUME] WARN myosin restore failed: {exc}", flush=True)
    chains_stacked = adapter.chains_tag_stacked

    print(f"[gate-b/{name}] n_part={sim.state.N_particles} R={R_cell*1e6:.2f}µm "
          f"dt={dtc:.2e}s F_crit={info['F_crit_pN']:.2f}pN τ_bend={info['tau_bend_us']:.0f}µs "
          f"relaxed={relaxed} start={start_tick}/{args.ticks}", flush=True)

    t0 = time.perf_counter()
    for k in range(start_tick, args.ticks):
        sim.run(args.interval)
        if (k + 1) % args.measure_every == 0 or k == args.ticks - 1:
            g_soft = float(_tension_method_of_planes(sim, R_cell))
            g_rigid = float(_tension_method_of_planes_rigid(sim, adapter, R_cell, dtc))
            s_abs, s_rel = _mean_s_grip(cell.myosin_action, ell0)
            cond = _condensation(sim, chains_stacked, ell0)
            row = dict(tick=k + 1, g_soft_mN_m=g_soft * 1e3, g_rigid_mN_m=g_rigid * 1e3,
                       g_struct_mN_m=(g_soft + g_rigid) * 1e3, s_grip_over_l0=s_rel,
                       nonconv=int(nat.nonconverged_count), **cond)
            rows.append(row)
            sps = ((k + 1 - start_tick) * args.interval) / max(time.perf_counter() - t0, 1e-9)
            print(f"  [{name}] tick {k+1}/{args.ticks} γ_struct={row['g_struct_mN_m']:.4e} "
                  f"(soft={row['g_soft_mN_m']:.3e} rigid={row['g_rigid_mN_m']:.3e}) "
                  f"s_grip/ℓ0={s_rel:.3f} bond/ℓ0={cond['mean_bond_over_l0']:.4f} "
                  f"nonconv={row['nonconv']} {sps:.0f} st/s", flush=True)
        if (k + 1) % args.ckpt_every == 0 or k == args.ticks - 1:
            ma = cell.myosin_action
            _save_ckpt(ckpt, tick=k + 1, positions=_tagpos(sim),
                       myo_s=(np.asarray(_to_host(ma._head_grip_s)) if ma else None),
                       myo_bound=(np.asarray(_to_host(ma._head_bound_to_actin)) if ma else None),
                       rows=rows)
            cond_json.write_text(json.dumps({"name": name, "mode": mode, "myo_on": myo_on,
                                             "info": info, "rows": rows}, indent=1))
    cond_json.write_text(json.dumps({"name": name, "mode": mode, "myo_on": myo_on,
                                     "info": info, "rows": rows}, indent=1))
    print(f"  cond → {cond_json}", flush=True)
    return 0


def _plateau_mean(rows, last_n, key):
    vals = [r[key] for r in rows[-last_n:]] if rows else []
    return float(np.mean(vals)) if vals else 0.0


def _aggregate(args):
    last_n = args.plateau_last
    conds = {}
    for name in _CONDS:
        p = _OUT_DIR / f"h7_{args.tag}_{name}.cond.json"
        if not p.exists():
            print(f"  missing {name} ({p.name})", flush=True)
            continue
        d = json.loads(p.read_text())
        d["g_struct_plateau"] = _plateau_mean(d["rows"], last_n, "g_struct_mN_m")
        d["s_grip_plateau"] = _plateau_mean(d["rows"], last_n, "s_grip_over_l0")
        d["bond_over_l0_last"] = d["rows"][-1]["mean_bond_over_l0"] if d["rows"] else None
        conds[name] = d
    missing = [n for n in _CONDS if n not in conds]
    if missing:
        print(f"[aggregate] incomplete — missing {missing}", flush=True)
        return 1

    def _s(n):
        return conds[n]["g_struct_plateau"]
    ga_rigid = _s("rigid_myoON") - _s("rigid_myoOFF")
    ga_relaxed = _s("relaxed_myoON") - _s("relaxed_myoOFF")
    g = ga_relaxed
    if CONFIRM_BAND[0] <= g <= CONFIRM_BAND[1]:
        verdict = (f"CONFIRM: γ_active={g:.4f} ∈ {CONFIRM_BAND} — buckling IS the lever.")
    elif g >= PARTIAL_FLOOR:
        verdict = (f"PARTIAL: γ_active={g:.4f} ≥ {PARTIAL_FLOOR} but < {CONFIRM_BAND[0]} — "
                   "contributes, insufficient.")
    else:
        verdict = (f"REFUTE: γ_active={g:.4f} < {PARTIAL_FLOOR} — buckling NOT the lever.")

    summary = {
        "tag": args.tag, "confirm_band_mN_m": list(CONFIRM_BAND),
        "gamma_active_rigid_mN_m": ga_rigid, "gamma_active_relaxed_mN_m": ga_relaxed,
        "gate_a_plateau_mN_m": GATE_A_PLATEAU_MN_M, "verdict": verdict,
        "conditions": {n: {"g_struct_plateau_mN_m": conds[n]["g_struct_plateau"],
                           "s_grip_plateau": conds[n]["s_grip_plateau"],
                           "bond_over_l0_last": conds[n]["bond_over_l0_last"]}
                       for n in _CONDS},
    }
    out = _OUT_DIR / f"h7_{args.tag}_SUMMARY.json"
    out.write_text(json.dumps(summary, indent=2))
    print("=" * 72, flush=True)
    print(f"  γ_active(rigid)   = {ga_rigid:.4e} mN/m  (control ≈ Gate-A {GATE_A_PLATEAU_MN_M:.2e})",
          flush=True)
    print(f"  γ_active(relaxed) = {ga_relaxed:.4e} mN/m  ← GATE", flush=True)
    for n in _CONDS:
        c = conds[n]
        print(f"    {n:16s} γ_struct={c['g_struct_plateau']:.4e} s_grip={c['s_grip_plateau']:.3f} "
              f"bond/ℓ0={c['bond_over_l0_last']:.4f}", flush=True)
    print(f"  VERDICT: {verdict}", flush=True)
    print(f"  → {out}", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["rigid", "relaxed"], default="relaxed")
    ap.add_argument("--myo", choices=["on", "off"], default="on")
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--warmup", type=int, default=4000)
    ap.add_argument("--ticks", type=int, default=600)
    ap.add_argument("--interval", type=int, default=5000)
    ap.add_argument("--measure-every", type=int, default=2)
    ap.add_argument("--ckpt-every", type=int, default=5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    ap.add_argument("--no-connected-mesh", action="store_true",
                    help="DIAGNOSTIC ONLY — use the fragmented mesh (giant ~7%). Default is "
                         "the connected spanning mesh (2026-06-04 rebuild), REQUIRED for a "
                         "transmission measurement. Gate-A/old runs wrongly used fragmented.")
    ap.add_argument("--no-nucleus", action="store_true",
                    help="drop the nucleus (sanctioned cortical-tension exception, PI "
                         "2026-06-08) — removes its stiff-lamin CFL so dt rises to the "
                         "membrane CFL; pair with --dt-safety ~5")
    ap.add_argument("--dt-safety", type=float, default=1.0,
                    help="constrained_dt_safety multiple (no-nucleus allows ~5×, the "
                         "membrane CFL limit; every remaining compartment stays in-guard)")
    ap.add_argument("--tag", type=str, default="gateB")
    ap.add_argument("--plateau-last", type=int, default=20, help="ticks averaged for the plateau")
    ap.add_argument("--aggregate", action="store_true")
    args = ap.parse_args()

    if args.aggregate:
        return _aggregate(args)
    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    return _run_condition(args, dev)


if __name__ == "__main__":
    raise SystemExit(main())
