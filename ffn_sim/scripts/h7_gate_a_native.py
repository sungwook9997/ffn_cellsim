"""H.7 Gate-A on the full GPU-main stack — does active γ rise as s_grip develops?

Gate-A question (H7_SESSION_HANDOFF §4#3): the active cortical tension γ_soft is
floored because the myosin grip-stretch ``s_grip`` stays ≈0 (heads bind+load but
don't WALK). Run LONG enough for ``s_grip`` to grow toward ℓ₀ and watch γ vs
s_grip:
  * CONFIRM (Gate-A was the wall): γ_soft / γ_rigid climb toward the band as
    s_grip → ℓ₀.
  * REFUTE  (deeper wall): γ stays flat while s_grip → ℓ₀ → Gate-B (buckling).

This is the production Gate-A run on the **full GPU-main stack** (native
constrained integrator + native compartment ForceCompute, 5.86× / ~4.7 d per
2×10⁸ — H7_NATIVE_FULLCELL_GO_2026-06-07.md), at the real physiological operating
point (build_baseline_cell: all compartments + grip_walk myosin + xlinks ON). It
reuses the validated two-phase build (unconstrained warm-up → constrained
production seeded from warm positions) and the native integrator swap from
``h7_native_fullcell_go``. Checkpoints positions + myosin grip state + trajectory
each ``--ckpt-every`` ticks so a multi-day run survives a crash (resume reseeds
both, NOT just positions — s_grip is the whole point).

Usage (gbook):
    export PYTHONPATH="$PWD/native/ffn_hoomd_plugin/build:$PWD/native/ffn_hoomd_plugin/python:$PYTHONPATH"
    export FFN_GPU_DEVICE_COMPARTMENTS=1
    python -u -m ffn_sim.scripts.h7_gate_a_native --n-fil 1000 --warmup 4000 \
        --ticks 4000 --interval 5000 --device gpu --tag gateA_prod
"""

from __future__ import annotations

import argparse
import json
import math
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
BAND = (0.35, 0.65)  # active cortical-tension band [mN/m]


def _mean_s_grip_over_l0(ma, rest_length: float) -> tuple[float, float]:
    """Mean grip-stretch of BOUND myosin heads, absolute [m] and / ℓ₀."""
    if ma is None:
        return 0.0, 0.0
    bound = np.asarray(_to_host(ma._head_bound_to_actin)) >= 0
    if not np.any(bound):
        return 0.0, 0.0
    s = float(np.mean(np.asarray(_to_host(ma._head_grip_s))[bound]))
    return s, s / rest_length


def _build_native_constrained(manifest, *, device, seed, warm_pos):
    """Constrained full cell (compartments via env) + native integrator swap."""
    import hoomd
    from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater

    cell = build_baseline_cell(
        manifest=deepcopy(manifest), device=device, seed=seed,
        constrained=True, equilibrate=False,
    )
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
    sim.operations.updaters.remove(cell.baoab_updater)
    nat = NativeConstrainedBaoabUpdater(
        dt=dtc, kT=float(act.kT), inv_gamma_by_tag=inv_gamma,
        chains_tag=chains_stacked.reshape(-1).astype(np.int32),
        n_chains=F, bonds_per_chain=npc - 1, rest_length=r0,
        seed=int(act._seed), tol=1.0e-9, max_iter=200,
        trigger=hoomd.trigger.Periodic(1),
    )
    sim.operations.updaters.append(nat)
    sim.run(0)
    adapter = _NativeLambdaAdapter(nat, chains_stacked, r0)
    return cell, nat, adapter, dtc, r0


def _save_ckpt(path: Path, *, tick, positions, myo_s, myo_bound, rows):
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, tick=tick, positions=positions,
             myo_s=myo_s if myo_s is not None else np.empty(0),
             myo_bound=myo_bound if myo_bound is not None else np.empty(0),
             rows=json.dumps(rows))
    os.replace(tmp, path)


def _make_figure(rows, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = np.array([r["tick"] for r in rows])
    gs = np.array([r["g_soft_mN_m"] for r in rows])
    gr = np.array([r["g_rigid_mN_m"] for r in rows])
    sg = np.array([r["s_grip_over_l0"] for r in rows])
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ax.axhspan(BAND[0], BAND[1], color="0.85", label=f"band {BAND} mN/m")
    ax.plot(t, gs, "-o", ms=3, color="tab:blue", label="γ_soft (mN/m)")
    ax.plot(t, gr, "-s", ms=3, color="tab:purple", alpha=0.7, label="γ_rigid (mN/m)")
    ax.set_xlabel("tick"); ax.set_ylabel("γ (mN/m)"); ax.set_yscale("log")
    ax2 = ax.twinx()
    ax2.plot(t, sg, "-^", ms=3, color="green", alpha=0.7, label="s_grip/ℓ₀")
    ax2.set_ylabel("s_grip / ℓ₀", color="green"); ax2.set_ylim(0, 1.05)
    ax.set_title("H.7 Gate-A (GPU-main): does γ rise as s_grip → ℓ₀?")
    ax.legend(loc="upper left", fontsize=8); ax2.legend(loc="lower right", fontsize=8)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--warmup", type=int, default=4000)
    ap.add_argument("--ticks", type=int, default=4000)
    ap.add_argument("--interval", type=int, default=5000, help="steps per tick")
    ap.add_argument("--measure-every", type=int, default=1, help="ticks between γ samples")
    ap.add_argument("--ckpt-every", type=int, default=5, help="ticks between checkpoints")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    ap.add_argument("--allow-cpu-dev", action="store_true")
    ap.add_argument("--tag", type=str, default="gateA_native")
    args = ap.parse_args()

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    manifest.setdefault("cortex_overrides", {}).setdefault(
        "cortex", {})["n_filaments"] = int(args.n_fil)
    softstart = max(300, args.warmup // 8)
    ckpt_path = _OUT_DIR / f"h7_gate_a_native_{args.tag}.ckpt.npz"
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    _FIG_DIR.mkdir(parents=True, exist_ok=True)

    # ---- warm-up (unconstrained; or restore from checkpoint) ----
    rows, start_tick, warm_pos = [], 0, None
    if ckpt_path.exists():
        z = np.load(ckpt_path, allow_pickle=True)
        warm_pos = z["positions"]
        start_tick = int(z["tick"])
        rows = json.loads(str(z["rows"]))
        ckpt_myo_s = z["myo_s"] if z["myo_s"].size else None
        ckpt_myo_bound = z["myo_bound"] if z["myo_bound"].size else None
        print(f"[RESUME] checkpoint at tick {start_tick}; {len(rows)} rows", flush=True)
    else:
        cell_w = build_baseline_cell(
            manifest=deepcopy(manifest), device=dev, seed=args.seed,
            constrained=False, equilibrate=True, equilibrate_steps=args.warmup,
            equilibrate_softstart_steps=softstart,
        )
        cell_w.simulation.run(0)
        warm_pos = _tagpos(cell_w.simulation)
        del cell_w
        ckpt_myo_s = ckpt_myo_bound = None

    # ---- constrained production on the GPU-main stack ----
    cell, nat, adapter, dtc, r0 = _build_native_constrained(
        manifest, device=dev, seed=args.seed, warm_pos=warm_pos)
    sim = cell.simulation
    R_cell = float(cell.p_cortex.R_cell)
    ell0 = float(r0)
    # Restore myosin grip state on resume (s_grip is the whole point).
    if ckpt_myo_s is not None and cell.myosin_action is not None:
        try:
            cell.myosin_action._head_grip_s[:] = ckpt_myo_s
            cell.myosin_action._head_bound_to_actin[:] = ckpt_myo_bound
            print("[RESUME] myosin grip state restored", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[RESUME] WARN myosin restore failed: {exc}", flush=True)

    nca = cell.p_cortex.n_filaments * cell.p_cortex.beads_per_filament
    r0_mean = float(np.linalg.norm(_tagpos(sim)[:nca], axis=1).mean())
    print(f"[gate-a] n_part={sim.state.N_particles} R_cell={R_cell*1e6:.2f}µm "
          f"dtc={dtc:.3e}s ℓ0={ell0:.3e}m start_tick={start_tick}/{args.ticks}", flush=True)

    t0 = time.perf_counter()
    for k in range(start_tick, args.ticks):
        sim.run(args.interval)
        if (k + 1) % args.measure_every == 0 or k == args.ticks - 1:
            g_soft = float(_tension_method_of_planes(sim, R_cell))
            g_rigid = float(_tension_method_of_planes_rigid(sim, adapter, R_cell, dtc))
            s_abs, s_rel = _mean_s_grip_over_l0(cell.myosin_action, ell0)
            r = _tagpos(sim)
            r_over_r0 = float(np.linalg.norm(r[:nca], axis=1).mean()) / r0_mean
            row = dict(tick=k + 1, g_soft_mN_m=g_soft * 1e3, g_rigid_mN_m=g_rigid * 1e3,
                       s_grip_m=s_abs, s_grip_over_l0=s_rel,
                       backbone_r_over_r0=r_over_r0,
                       nonconv=int(nat.nonconverged_count))
            rows.append(row)
            sps = ((k + 1 - start_tick) * args.interval) / max(time.perf_counter() - t0, 1e-9)
            print(f"  [gateA] tick {k+1}/{args.ticks} γ_soft={row['g_soft_mN_m']:.4e} "
                  f"γ_rigid={row['g_rigid_mN_m']:.4e} s_grip/ℓ0={s_rel:.3f} "
                  f"r/r0={r_over_r0:.4f} nonconv={row['nonconv']} {sps:.0f} steps/s",
                  flush=True)
        if (k + 1) % args.ckpt_every == 0 or k == args.ticks - 1:
            ma = cell.myosin_action
            _save_ckpt(ckpt_path, tick=k + 1, positions=_tagpos(sim),
                       myo_s=(np.asarray(_to_host(ma._head_grip_s)) if ma else None),
                       myo_bound=(np.asarray(_to_host(ma._head_bound_to_actin)) if ma else None),
                       rows=rows)
            _make_figure(rows, _FIG_DIR / f"h7_gate_a_native_{args.tag}.png")
            (_OUT_DIR / f"h7_gate_a_native_{args.tag}.json").write_text(json.dumps(rows, indent=1))

    # ---- verdict ----
    last = rows[-1] if rows else {}
    s_grown = last.get("s_grip_over_l0", 0.0) >= 0.5
    g_band = last.get("g_soft_mN_m", 0.0) >= BAND[0]
    if s_grown and g_band:
        verdict = "CONFIRM: s_grip developed AND γ_soft reached the band → Gate-A was the wall."
    elif s_grown and not g_band:
        verdict = ("REFUTE: s_grip developed but γ_soft floored below band → deeper wall "
                   "(Gate-B buckling / lever).")
    else:
        verdict = (f"INCONCLUSIVE: s_grip did not develop enough "
                   f"(s_grip/ℓ0={last.get('s_grip_over_l0', 0):.2f}) — extend --ticks.")
    print(f"\n[gate-a VERDICT] {verdict}", flush=True)
    (_OUT_DIR / f"h7_gate_a_native_{args.tag}.json").write_text(json.dumps(
        {"verdict": verdict, "rows": rows}, indent=1))
    _make_figure(rows, _FIG_DIR / f"h7_gate_a_native_{args.tag}.png")
    print(f"→ figs/h7_gate_a_native_{args.tag}.png", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
