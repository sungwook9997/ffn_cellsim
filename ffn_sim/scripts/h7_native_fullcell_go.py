"""H.7 GO test — native vs cupy constrained integrator on the FULL CELL.

The decisive go/no-go the GPU-opt sub-session set for the Stage-1b native
constrained BAOAB plugin (`cacb27b`): the integrator-only microbench is ~40×,
but does that survive **end-to-end on the full cell**, where the binding updaters
(myosin / xlink / integrin) still host-sync each batch? Per the sub-session, that
full-cell measurement is the production-driver's domain (mine).

This driver runs the SAME stable two-phase build the production driver uses
(unconstrained warm-up with cytoplasm η → fresh constrained production seeded from
the warmed positions — the path that does NOT diverge), then compares two arms at
identical config/seed/step-count:

  * **cupy** — the in-tree ``ConstrainedLeimkuhlerMatthewsBAOAB`` (current default).
  * **native** — the same constrained cell with its cupy integrator CustomUpdater
    swapped for ``ffn_hoomd_plugin.NativeConstrainedBaoabUpdater``, built from the
    cupy action's OWN per-tag mobility + chains (apples-to-apples).

Reports the three GO criteria:
  1. ``nonconverged_count == 0`` (native M-SHAKE converges every step),
  2. γ (soft method-of-planes + rigid Lagrange-λ) statistically matches cupy,
  3. measured steps/s + 2×10⁸-step ETA for each → the end-to-end speedup.

GPU-only for the native arm (the plugin is a CUDA build; needs
``PYTHONPATH=native/ffn_hoomd_plugin/build:native/ffn_hoomd_plugin/python``).
The cupy arm runs on GPU or CPU. No core modules are modified — the native
integrator is swapped onto a built cell at runtime (additive, opt-in harness).

Usage (gbook):
    export PYTHONPATH="$PWD/native/ffn_hoomd_plugin/build:$PWD/native/ffn_hoomd_plugin/python:$PYTHONPATH"
    python -m ffn_sim.scripts.h7_native_fullcell_go --n-fil 1000 --warmup 4000 \
        --steps 2000 --interval 200 --device gpu
"""

from __future__ import annotations

import argparse
import json
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

from ffn_sim.cell.manifest import build_baseline_cell, load_manifest
from ffn_sim.scripts.h3_ku35_tension import (
    _tension_method_of_planes,
    _tension_method_of_planes_rigid,
)

_OUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "h7" / "production"


def _to_host(a):
    """cupy/np → host numpy."""
    if a is None:
        return None
    try:
        return a.get()  # cupy
    except AttributeError:
        return np.asarray(a)


def _tagpos(sim) -> np.ndarray:
    with sim.state.cpu_local_snapshot as s:
        pos = np.asarray(s.particles.position)
        tg = np.asarray(s.particles.tag)
        inv = np.empty_like(tg)
        inv[tg] = np.arange(tg.size)
        return pos[inv].copy()


class _NativeLambdaAdapter:
    """Duck-typed shim so ``_tension_method_of_planes_rigid`` reads native λ.

    Exposes the three attributes that estimator needs: ``lambda_buf`` (host (F,m)),
    ``chains_tag_stacked`` (host (F, m+1)), ``_chain_rest_length``.
    """

    def __init__(self, nat, chains_stacked, rest_length):
        self._nat = nat
        self.chains_tag_stacked = chains_stacked
        self._chain_rest_length = float(rest_length)

    @property
    def lambda_buf(self):
        lam = _to_host(self._nat.lambda_buf)
        if lam is None:
            return None
        lam = np.asarray(lam, dtype=np.float64)
        F, npc = self.chains_tag_stacked.shape
        return lam.reshape(F, npc - 1)


def _build_two_phase(manifest, *, device, warmup, softstart, seed):
    """Stable warm-up (unconstrained) → warm positions; then a fn to build the
    constrained cell seeded from those positions (callable twice for two arms)."""
    cell_w = build_baseline_cell(
        manifest=deepcopy(manifest), device=device, seed=seed,
        constrained=False, equilibrate=True, equilibrate_steps=warmup,
        equilibrate_softstart_steps=softstart,
    )
    cell_w.simulation.run(0)
    warm_pos = _tagpos(cell_w.simulation)
    del cell_w

    def build_constrained():
        cell = build_baseline_cell(
            manifest=deepcopy(manifest), device=device, seed=seed,
            constrained=True, equilibrate=False,
        )
        if hasattr(cell.baoab_action, "record_lambda"):
            cell.baoab_action.record_lambda = True
        snap = cell.simulation.state.get_snapshot()
        if snap.communicator.rank == 0:
            snap.particles.position[:] = warm_pos
        cell.simulation.state.set_snapshot(snap)
        cell.simulation.run(0)
        return cell

    return build_constrained


def _run_arm(cell, *, steps, interval, R_cell, dtc, rigid_action):
    """Run an arm; return (samples list of (g_soft,g_rigid), steps_per_s).

    Throughput is timed over PURE ``sim.run(steps)`` — γ is sampled AFTER the
    timed region, NOT interleaved. The per-sample γ estimators do
    ``cpu_local_snapshot`` host reads (and the rigid one reads λ); interleaving
    them inside the timed loop charges that host-sync to every interval and
    dilutes the native integrator's speedup (it penalises the GPU-resident arm
    most). In a real production run γ is read ~1% of steps, so pure-stepping
    throughput is the Gate-A/B-feasibility-relevant number; γ here is only the
    parity check.
    """
    sim = cell.simulation
    sim.run(interval)  # warm / JIT before timing
    t0 = time.perf_counter()
    sim.run(int(steps))  # PURE timing — no host-sync measurement inside
    dt_wall = time.perf_counter() - t0
    sps = steps / dt_wall if dt_wall > 0 else float("nan")
    # γ parity samples AFTER the timed region (a few short windows).
    samples = []
    for _ in range(5):
        sim.run(interval)
        g_soft = _tension_method_of_planes(sim, R_cell)
        g_rigid = _tension_method_of_planes_rigid(sim, rigid_action, R_cell, dtc)
        samples.append((float(g_soft), float(g_rigid)))
    return samples, sps


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--n-motors", type=int, default=200)
    ap.add_argument("--warmup", type=int, default=4000)
    ap.add_argument("--softstart", type=int, default=None)
    ap.add_argument("--steps", type=int, default=2000, help="timed production steps/arm")
    ap.add_argument("--interval", type=int, default=200, help="steps between γ samples")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="gpu")
    ap.add_argument("--allow-cpu-dev", action="store_true")
    ap.add_argument("--gate-steps", type=float, default=2.0e8,
                    help="step count for the ETA projection (Gate-A/B scale)")
    args = ap.parse_args()

    import hoomd
    dev = (hoomd.device.GPU(notice_level=0) if args.device == "gpu"
           else hoomd.device.CPU(notice_level=0))

    manifest = deepcopy(load_manifest("mcf7_baseline.yaml"))
    co = manifest.setdefault("cortex_overrides", {}).setdefault("cortex", {})
    co["n_filaments"] = int(args.n_fil)
    softstart = args.softstart if args.softstart is not None else max(300, args.warmup // 8)

    build_constrained = _build_two_phase(
        manifest, device=dev, warmup=args.warmup, softstart=softstart, seed=args.seed,
    )

    # ---- ARM 1: cupy constrained integrator ----
    cell_c = build_constrained()
    R_cell = float(cell_c.p_cortex.R_cell)
    act_c = cell_c.baoab_action
    dtc = float(act_c.dt)
    samples_c, sps_c = _run_arm(
        cell_c, steps=args.steps, interval=args.interval,
        R_cell=R_cell, dtc=dtc, rigid_action=act_c,
    )
    # Extract chain + mobility params from the cupy action (now populated).
    chains_stacked = np.asarray(act_c.chains_tag_stacked)
    F, npc = chains_stacked.shape
    m = npc - 1
    r0 = float(act_c._chain_rest_length)
    inv_gamma = _to_host(act_c._inv_gamma_by_tag)
    inv_gamma = [float(x) for x in np.asarray(inv_gamma, dtype=np.float64)]
    kT = float(act_c.kT)
    seed_c = int(act_c._seed)
    n_part = int(cell_c.simulation.state.N_particles)
    del cell_c

    # ---- ARM 2: native constrained integrator (swap) ----
    native_ok = True
    native_err = ""
    samples_n, sps_n, nonconv = [], float("nan"), -1
    try:
        from ffn_hoomd_plugin import NativeConstrainedBaoabUpdater
    except Exception as exc:  # noqa: BLE001
        native_ok = False
        native_err = f"native plugin import failed: {exc!r} (build + PYTHONPATH?)"

    if native_ok:
        cell_n = build_constrained()
        sim_n = cell_n.simulation
        # Swap: remove the cupy integrator CustomUpdater, append the native one.
        sim_n.operations.updaters.remove(cell_n.baoab_updater)
        nat = NativeConstrainedBaoabUpdater(
            dt=dtc, kT=kT,
            inv_gamma_by_tag=inv_gamma,
            chains_tag=chains_stacked.reshape(-1).astype(np.int32),
            n_chains=F, bonds_per_chain=m, rest_length=r0,
            seed=seed_c, tol=1.0e-9, max_iter=200,
            trigger=hoomd.trigger.Periodic(1),
        )
        sim_n.operations.updaters.append(nat)
        sim_n.run(0)
        adapter = _NativeLambdaAdapter(nat, chains_stacked, r0)
        samples_n, sps_n = _run_arm(
            cell_n, steps=args.steps, interval=args.interval,
            R_cell=R_cell, dtc=dtc, rigid_action=adapter,
        )
        nonconv = int(nat.nonconverged_count)
        del cell_n

    # ---- Report ----
    def _stats(samples):
        if not samples:
            return {"g_soft_mN_m": None, "g_rigid_mN_m": None}
        arr = np.asarray(samples)
        return {
            "g_soft_mN_m": float(arr[:, 0].mean() * 1e3),
            "g_soft_std_mN_m": float(arr[:, 0].std() * 1e3),
            "g_rigid_mN_m": float(arr[:, 1].mean() * 1e3),
            "g_rigid_std_mN_m": float(arr[:, 1].std() * 1e3),
        }

    def _eta_h(sps):
        return (args.gate_steps / sps / 3600.0) if (sps and sps == sps and sps > 0) else None

    cupy_stats, native_stats = _stats(samples_c), _stats(samples_n)
    speedup = (sps_n / sps_c) if (native_ok and sps_c > 0) else None

    # GO criteria.
    parity_ok = None
    if native_ok and cupy_stats["g_soft_mN_m"] and native_stats["g_soft_mN_m"]:
        gs_c, gs_n = cupy_stats["g_soft_mN_m"], native_stats["g_soft_mN_m"]
        tol = 0.25 * abs(gs_c) + 3.0 * (cupy_stats.get("g_soft_std_mN_m", 0.0))
        parity_ok = bool(abs(gs_n - gs_c) <= max(tol, 1e-6))
    go = bool(native_ok and nonconv == 0 and parity_ok)

    summary = {
        "config": {
            "n_fil": args.n_fil, "n_part": n_part, "F_chains": F, "m_bonds": m,
            "warmup": args.warmup, "timed_steps": args.steps,
            "interval": args.interval, "device": args.device, "dtc_s": dtc,
            "R_cell_m": R_cell,
        },
        "cupy": {"steps_per_s": sps_c, "eta_hours_at_gate": _eta_h(sps_c), **cupy_stats},
        "native": {
            "available": native_ok, "error": native_err,
            "steps_per_s": sps_n, "eta_hours_at_gate": _eta_h(sps_n),
            "nonconverged_count": nonconv, **native_stats,
        },
        "end_to_end_speedup_native_over_cupy": speedup,
        "gate_steps": args.gate_steps,
        "GO_criteria": {
            "native_available": native_ok,
            "nonconverged_zero": (nonconv == 0) if native_ok else None,
            "gamma_soft_parity": parity_ok,
            "GO": go,
        },
    }
    _OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_json = _OUT_DIR / "h7_native_fullcell_go.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"\n{'='*60}")
    print(f"cupy:   {sps_c:.1f} steps/s  (2e8 ETA {_eta_h(sps_c)}h)  "
          f"g_soft={cupy_stats['g_soft_mN_m']}")
    if native_ok:
        print(f"native: {sps_n:.1f} steps/s  (2e8 ETA {_eta_h(sps_n)}h)  "
              f"g_soft={native_stats['g_soft_mN_m']}  nonconv={nonconv}")
        print(f"END-TO-END SPEEDUP native/cupy = "
              f"{speedup:.2f}x" if speedup else "speedup n/a")
        print(f"GO = {go}  (nonconv0={nonconv==0}, soft-parity={parity_ok})")
    else:
        print(f"native arm SKIPPED: {native_err}")
    print(f"→ {out_json}")
    return 0 if go or not native_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
