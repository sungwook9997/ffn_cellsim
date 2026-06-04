#!/usr/bin/env python
"""Track-1 g_soft binding-equilibrium VERIFICATION run (KU-3.5, READ-ONLY physics).

WHY: The months-long KU-3.5 active cortical-tension floor (`g_soft` ~2e-4 mN/m,
~1000x under the band [0.35, 0.65] mN/m) was DIAGNOSED (see
`docs/CORTICAL_TENSION_TRIAGE_2026-06-03.md` §0-2, `track1_binding_diag.py`,
`outputs/h3/track1_binding_diag.json`) as a **batch_dt rate-throttle ARTIFACT,
not a missing mechanism**:

    p_bind/tick = 1 - exp(-k_on * batch_dt) = 1 - exp(-50 * 1.30e-6) = 6.5e-5

so the head->actin bound fraction needs ~1/(k_on+k_off0) / batch_dt ~ 1.3e4 ticks
(~1.7e9 BAOAB steps) to reach its k_on/(k_on+k_off0) ~ 0.83 equilibrium. Short
production runs sampled only ~0.1% bound -> the ~10^3 g_soft deficit. Geometry,
gates, off-rate, and s_grip stepping were all verified fine by the diagnostic.

THIS SCRIPT empirically tests the artifact hypothesis: it builds the SAME small
cortex+myosin cell and runs it to head->actin BINDING EQUILIBRIUM, measuring at
intervals (a) bound fraction, (b) mean F/F_stall, (c) g_soft [mN/m] via the
KU-3.5 method-of-planes estimator (the exact `_tension_method_of_planes` from
`h3_ku35_tension.py`).

It reaches equilibrium FAST by raising ONLY the NUMERICAL batching parameter
`batch_steps` (a script-level config override -- NO committed config is edited;
this decoupling IS the proposed fix). The D2 CFL constraint
``batch_dt * k_off_max <= 1e-3`` (the `resolve_cortex_myosin` §2 gate) is checked
and respected -- the resolver auto-clamps batch_steps if a chosen value would
violate it, and this script reports the realised CFL product.

Single Simulation, pre-allocated, NO rebuilds, NO parameter sweep (one cortex,
one long run) -- respects the HOOMD-rebuild memory-leak rule. CPU is fine.

Usage (small / fast, CPU):
    python ffn_sim/scripts/track1_gsoft_verify.py --n-fil 300 --batch-steps 4000 \
        --ticks 600 --measure-every 20
"""
from __future__ import annotations

import argparse
import json
import math
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
import yaml
import hoomd

from ffn_sim.cortex.cortex import resolve_h3_derived
from ffn_sim.cortex.myosin import resolve_cortex_myosin
from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
from ffn_sim.cell.cell import build_cortex_full_simulation

# Reuse the EXACT KU-3.5 g_soft estimator from the production driver (no copy).
from ffn_sim.scripts.h3_ku35_tension import _tension_method_of_planes

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"
OUT_JSON = PKG / "outputs" / "h3" / "track1_gsoft_verify.json"
OUT_PNG = PKG / "outputs" / "h3" / "track1_gsoft_verify.png"

# KU-3.5 active cortical-tension band [mN/m] and the diagnosed floor.
BAND = (0.35, 0.65)
FLOOR_MN_PER_M = 2.0e-4  # the ~1000x-under-band g_soft floor (short-run sampling)


def _tagpos(snap) -> np.ndarray:
    """Positions in tag order (single-rank get_snapshot is already tag-ordered)."""
    return np.asarray(snap.particles.position, dtype=np.float64)


def _bound_metrics(ma, snap) -> tuple[float, float, int]:
    """Return (bound_fraction, mean_F_over_stall_on_bound, n_bound)."""
    p = ma.p
    pos = _tagpos(snap)
    bound = ma._head_bound_to_actin >= 0
    n_bound = int(bound.sum())
    n_heads = ma._head_bound_to_actin.size
    if n_bound == 0:
        return 0.0, 0.0, 0
    bidx = np.flatnonzero(bound)
    htags = np.array([ma._head_global_tag(int(h)) for h in bidx], dtype=np.int64)
    atags = ma._head_bound_to_actin[bidx]
    r = np.linalg.norm(pos[htags] - pos[atags], axis=1)
    F = p.k_head_actin * np.clip(r, 0.0, None)  # r0 ~ 0 (grip_walk)
    f_over_stall = float(np.mean(F / p.F_stall_per_head))
    return float(n_bound / n_heads), f_over_stall, n_bound


def run(
    n_fil: int,
    batch_steps_override: int,
    ticks: int,
    measure_every: int,
    seed: int,
    stepping_mode: str,
    n_motors: int | None = None,
    with_xlinks: bool = True,
) -> dict:
    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["n_filaments"] = n_fil
    cfg["cortex"]["demo_mode"] = True
    cfg.setdefault("cortex", {}).setdefault("myosin", {})
    cfg["cortex"]["myosin"]["stepping_mode"] = stepping_mode
    if n_motors is not None:
        # Numerical/system-size knob (smaller cell -> faster CPU run-to-eq).
        # Does NOT change per-head binding kinetics, only the head count the
        # bound FRACTION is computed over -> the equilibrium fraction is
        # n_motors-invariant (it is k_on/(k_on+k_off0)).
        cfg["cortex"]["myosin"]["n_motors_per_cell"] = int(n_motors)
    # ---- THE FIX UNDER TEST: raise ONLY the numerical batch_steps ----------
    # Script-level override; the resolver re-reads it and the §2 D2 CFL gate
    # auto-clamps if batch_dt*k_off_max would exceed 1e-3. No committed file
    # is edited -- this exact decoupling is the proposed config change.
    cfg["cortex"]["myosin"]["batch_steps"] = int(batch_steps_override)

    p_cortex = resolve_h3_derived(cfg)
    p_myo = resolve_cortex_myosin(cfg, dt=p_cortex.dt_cfl, R_cell=p_cortex.R_cell)
    p_xl = None
    if with_xlinks:
        try:
            p_xl = resolve_crosslinkers(cfg, dt=p_cortex.dt_cfl)
        except Exception:
            p_xl = None

    # --- CFL verification (D2 batch): batch_dt * k_off_max <= 1e-3 ----------
    cfl_product = p_myo.batch_dt * p_myo.k_off_max_at_zero_load
    cfl_ok = cfl_product <= 1.0e-3 + 1e-12
    clamped = "batch_steps_shrunk_from" in p_myo.extras
    p_bind_tick = 1.0 - math.exp(-p_myo.head_actin_k_on * p_myo.batch_dt)
    p_break_tick = 1.0 - math.exp(-p_myo.head_actin_k_off0 * p_myo.batch_dt)
    naive_eq = p_myo.head_actin_k_on / (
        p_myo.head_actin_k_on + p_myo.head_actin_k_off0
    )
    # Relaxation time of the two-state binding ODE ~ 1/(k_on + k_off0).
    tau_relax_s = 1.0 / (p_myo.head_actin_k_on + p_myo.head_actin_k_off0)
    ticks_per_tau = tau_relax_s / p_myo.batch_dt
    total_steps = ticks * p_myo.batch_steps

    print(
        f"[build] n_fil={n_fil} beads/fil={p_cortex.beads_per_filament} "
        f"n_motors={p_myo.n_motors_per_cell} "
        f"heads={2*p_myo.n_heads_per_side*p_myo.n_motors_per_cell} "
        f"mode={p_myo.stepping_mode}\n"
        f"        requested batch_steps={batch_steps_override} -> "
        f"RESOLVED batch_steps={p_myo.batch_steps} "
        f"(clamped_by_CFL={clamped})\n"
        f"        dt_cfl={p_cortex.dt_cfl:.3e}s batch_dt={p_myo.batch_dt:.3e}s "
        f"k_on={p_myo.head_actin_k_on}/s k_off0={p_myo.head_actin_k_off0}/s\n"
        f"        CFL product batch_dt*k_off_max={cfl_product:.3e} "
        f"(<=1e-3 ? {cfl_ok})\n"
        f"        p_bind/tick={p_bind_tick:.4e} p_break/tick={p_break_tick:.4e} "
        f"naive_eq_bound_frac={naive_eq:.3f}\n"
        f"        tau_relax={tau_relax_s:.3e}s = {ticks_per_tau:.1f} ticks; "
        f"running {ticks} ticks = {ticks/ticks_per_tau:.1f}x tau "
        f"({total_steps:,} BAOAB steps)\n"
        f"        hoomd={hoomd.version.version}",
        flush=True,
    )
    if not cfl_ok:
        raise RuntimeError(
            f"D2 CFL violated even after resolver clamp: product={cfl_product:.3e}"
        )

    device = hoomd.device.CPU(notice_level=0)
    hw = build_cortex_full_simulation(
        p_cortex, p_xlinks=p_xl, p_myosin=p_myo,
        device=device, with_baoab=True,
        rng=np.random.default_rng(seed),
    )
    sim = hw["sim"]
    ma = hw["myosin_action"]
    if ma is None:
        raise RuntimeError("Could not locate MyosinStepUpdater on the sim.")

    sim.run(0)  # initialize updater

    # g_soft of the freshly-built (essentially unbound) cortex == the floor we
    # are trying to lift. Record it as the baseline.
    snap0 = sim.state.get_snapshot()
    g_soft0 = _tension_method_of_planes(sim, p_cortex.R_cell) * 1e3  # mN/m
    bf0, fos0, nb0 = _bound_metrics(ma, snap0)
    print(
        f"[t=0] bound_frac={bf0:.4f} (n_bound={nb0}) "
        f"F/F_stall={fos0:.3f} g_soft={g_soft0:.4e} mN/m",
        flush=True,
    )

    rows = []
    t0 = time.time()
    # One long run, sampled every `measure_every` ticks (= measure_every batches).
    n_measures = max(1, ticks // measure_every)
    for m in range(n_measures):
        sim.run(measure_every * p_myo.batch_steps)
        snap = sim.state.get_snapshot()
        bf, fos, nb = _bound_metrics(ma, snap)
        g_soft = _tension_method_of_planes(sim, p_cortex.R_cell) * 1e3  # mN/m
        tick_now = (m + 1) * measure_every
        wall = time.time() - t0
        eta = wall * (n_measures - (m + 1)) / max(m + 1, 1)
        rows.append(dict(
            tick=int(tick_now),
            step=int(sim.timestep),
            bound_frac=float(bf),
            n_bound=int(nb),
            mean_F_over_stall=float(fos),
            g_soft_mN_per_m=float(g_soft),
            tau_multiple=float(tick_now / ticks_per_tau),
        ))
        print(
            f"PROGRESS tick={tick_now}/{ticks} ({tick_now/ticks_per_tau:.2f}x tau) "
            f"bound_frac={bf:.4f} (n={nb}) F/F_stall={fos:.3f} "
            f"g_soft={g_soft:.4e} mN/m  wall={wall:.0f}s eta={eta:.0f}s",
            flush=True,
        )
        # Incremental partial dump (survives an early stop on slow CPU): the
        # final structured result is still written at loop end, but this keeps
        # the latest trajectory on disk after every measurement.
        try:
            OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
            (OUT_JSON.with_name("track1_gsoft_verify_partial.json")).write_text(
                json.dumps(dict(
                    status="running",
                    config_batch_dt_s=float(p_myo.batch_dt),
                    cfl_product=float(cfl_product),
                    naive_eq_bound_frac=float(naive_eq),
                    ticks_per_tau=float(ticks_per_tau),
                    baseline_g_soft_mN_per_m=float(g_soft0),
                    rows=rows,
                ), indent=2)
            )
        except Exception:
            pass

    # Plateau metrics: last third of the measurements (binding equilibrium).
    last_third = rows[max(0, len(rows) * 2 // 3):]
    eq_bound_frac = float(np.mean([r["bound_frac"] for r in last_third]))
    eq_g_soft = float(np.mean([r["g_soft_mN_per_m"] for r in last_third]))
    eq_f_over_stall = float(np.mean([r["mean_F_over_stall"] for r in last_third]))
    peak_g_soft = float(np.max([r["g_soft_mN_per_m"] for r in rows]))
    peak_bound_frac = float(np.max([r["bound_frac"] for r in rows]))

    # 2C (PI Decision, 2026-06-04): physical-velocity acceptance filter. A myosin head
    # cannot deliver more than its stall force; any sample with mean_F_over_stall > 1 is
    # OVER-DRIVEN (the grip-walk delivered force k_head_actin·r exceeds the per-head Hill
    # stall — the documented "CH2 over-driving" leak that bought the prior v0x3000=0.030
    # mN/m number). Report g_soft over ONLY the physically-valid (F/F_stall<=1) plateau
    # samples so no g_soft can ever be claimed via super-stall force. For a literal-v0 run
    # (F/F_stall~0.2) this is a no-op; it only bites over-driven (high-v0_accel) runs.
    physical = [r for r in rows if r["mean_F_over_stall"] <= 1.0]
    n_physical = len(physical)
    phys_plateau = physical[max(0, n_physical * 2 // 3):]
    phys_g_soft = (float(np.mean([r["g_soft_mN_per_m"] for r in phys_plateau]))
                   if phys_plateau else float("nan"))
    overdriven_fraction = (1.0 - n_physical / len(rows)) if rows else float("nan")

    factor_vs_floor = eq_g_soft / FLOOR_MN_PER_M if FLOOR_MN_PER_M > 0 else float("nan")
    in_band = BAND[0] <= eq_g_soft <= BAND[1]
    above_floor = eq_g_soft > 10.0 * FLOOR_MN_PER_M  # >>10x the floor = lifted
    toward_band = eq_g_soft >= BAND[0]  # reached/exceeded band floor

    result = dict(
        verdict_hypothesis=(
            "CONFIRMED" if (above_floor and toward_band)
            else ("PARTIAL" if above_floor else "REFUTED")
        ),
        config=dict(
            n_fil=int(n_fil),
            beads_per_filament=int(p_cortex.beads_per_filament),
            n_motors=int(p_myo.n_motors_per_cell),
            n_heads=int(2 * p_myo.n_heads_per_side * p_myo.n_motors_per_cell),
            stepping_mode=p_myo.stepping_mode,
            seed=int(seed),
            R_cell_um=float(p_cortex.R_cell * 1e6),
            k_on_per_s=float(p_myo.head_actin_k_on),
            k_off0_per_s=float(p_myo.head_actin_k_off0),
            dt_cfl_s=float(p_cortex.dt_cfl),
        ),
        batch_steps_override=dict(
            requested=int(batch_steps_override),
            resolved=int(p_myo.batch_steps),
            clamped_by_CFL=bool(clamped),
            batch_dt_s=float(p_myo.batch_dt),
            cfl_product=float(cfl_product),
            cfl_limit=1.0e-3,
            cfl_ok=bool(cfl_ok),
            p_bind_per_tick=float(p_bind_tick),
            p_break_per_tick=float(p_break_tick),
            ticks=int(ticks),
            total_baoab_steps=int(total_steps),
            tau_relax_s=float(tau_relax_s),
            ticks_per_tau=float(ticks_per_tau),
            run_tau_multiple=float(ticks / ticks_per_tau),
        ),
        baseline=dict(
            g_soft_mN_per_m=float(g_soft0),
            bound_frac=float(bf0),
            floor_reference_mN_per_m=FLOOR_MN_PER_M,
        ),
        equilibrium=dict(
            bound_frac=eq_bound_frac,
            naive_eq_bound_frac=float(naive_eq),
            mean_F_over_stall=eq_f_over_stall,
            g_soft_mN_per_m=eq_g_soft,
            peak_g_soft_mN_per_m=peak_g_soft,
            peak_bound_frac=peak_bound_frac,
        ),
        comparison=dict(
            band_mN_per_m=list(BAND),
            floor_mN_per_m=FLOOR_MN_PER_M,
            factor_vs_floor=float(factor_vs_floor),
            eq_in_band=bool(in_band),
            eq_reached_band_floor=bool(toward_band),
            eq_lifted_above_floor=bool(above_floor),
        ),
        physical_validity=dict(  # 2C: Hill-bounded (F/F_stall<=1) g_soft only
            n_samples=len(rows),
            n_physical_samples=int(n_physical),
            overdriven_fraction=float(overdriven_fraction),
            phys_g_soft_mN_per_m=float(phys_g_soft),
            phys_in_band=bool(np.isfinite(phys_g_soft) and BAND[0] <= phys_g_soft <= BAND[1]),
            note=("g_soft over ONLY physically-valid samples (per-head force <= stall). "
                  "If overdriven_fraction>0 the unfiltered equilibrium g_soft was partly "
                  "bought with super-stall (non-Hill) force and is NOT a physical tension."),
        ),
        rows=rows,
    )

    print(
        f"[2C physical-validity] {n_physical}/{len(rows)} samples Hill-bounded "
        f"(F/F_stall<=1); over-driven fraction={overdriven_fraction:.2f}; "
        f"PHYSICAL g_soft={phys_g_soft:.4e} mN/m "
        f"(unfiltered eq={eq_g_soft:.4e}) -> "
        f"{'IN band' if (np.isfinite(phys_g_soft) and BAND[0]<=phys_g_soft<=BAND[1]) else 'under band'}",
        flush=True,
    )
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=2))
    print("RESULT_FILE " + str(OUT_JSON), flush=True)

    _plot(rows, result)
    return result


def _plot(rows: list[dict], result: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = np.array([r["tick"] for r in rows], dtype=float)
    bf = np.array([r["bound_frac"] for r in rows], dtype=float)
    gs = np.array([r["g_soft_mN_per_m"] for r in rows], dtype=float)
    naive_eq = result["batch_steps_override"]
    naive_eq_bf = result["equilibrium"]["naive_eq_bound_frac"]

    fig, axes = plt.subplots(2, 1, figsize=(9, 9), sharex=True)

    ax0 = axes[0]
    ax0.plot(t, bf, "-o", color="tab:blue", ms=3, label="bound fraction (measured)")
    ax0.axhline(
        naive_eq_bf, color="tab:blue", ls="--", lw=1,
        label=f"naive eq k_on/(k_on+k_off0) = {naive_eq_bf:.3f}",
    )
    ax0.set_ylabel("head->actin bound fraction")
    ax0.set_ylim(bottom=0.0)
    ax0.set_title(
        "Track-1 g_soft binding-equilibrium verification "
        f"(batch_steps={naive_eq['resolved']}, "
        f"batch_dt={naive_eq['batch_dt_s']:.2e}s, "
        f"CFL={naive_eq['cfl_product']:.1e}<=1e-3)"
    )
    ax0.legend(loc="best", fontsize=8)
    ax0.grid(alpha=0.3)

    ax1 = axes[1]
    ax1.plot(t, gs, "-o", color="tab:red", ms=3, label="g_soft (method-of-planes)")
    ax1.axhspan(BAND[0], BAND[1], color="green", alpha=0.15,
                label=f"KU-3.5 band [{BAND[0]}, {BAND[1]}] mN/m")
    ax1.axhline(FLOOR_MN_PER_M, color="black", ls=":", lw=1,
                label=f"diagnosed floor ~{FLOOR_MN_PER_M:.0e} mN/m")
    eq_g = result["equilibrium"]["g_soft_mN_per_m"]
    ax1.axhline(eq_g, color="tab:red", ls="--", lw=1,
                label=f"eq g_soft = {eq_g:.3e} mN/m")
    ax1.set_xlabel("myosin updater ticks")
    ax1.set_ylabel("g_soft  [mN/m]")
    ax1.set_yscale("log")
    ax1.legend(loc="best", fontsize=8)
    ax1.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130)
    plt.close(fig)
    print("FIG " + str(OUT_PNG), flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-fil", type=int, default=300)
    ap.add_argument(
        "--batch-steps", type=int, default=4000,
        help="NUMERICAL batch_steps override (the proposed fix). The §2 D2 CFL "
             "gate auto-clamps to <=1e-3/(k_off0*dt) if too large.",
    )
    ap.add_argument("--ticks", type=int, default=600,
                    help="total myosin updater ticks (>= ~3x relaxation time)")
    ap.add_argument("--measure-every", type=int, default=20,
                    help="measure bound-frac + g_soft every N ticks")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--stepping-mode", choices=["grip_walk", "binned_r0"],
                    default="grip_walk")
    ap.add_argument("--n-motors", type=int, default=None,
                    help="override n_motors_per_cell (system-size/speed knob; "
                         "eq bound FRACTION is n_motors-invariant)")
    ap.add_argument("--no-xlinks", dest="with_xlinks", action="store_false",
                    default=True,
                    help="omit crosslinkers (irrelevant to head->actin binding; "
                         "removes a soft-bond confound + per-tick overhead)")
    args = ap.parse_args()

    run(
        n_fil=args.n_fil,
        batch_steps_override=args.batch_steps,
        ticks=args.ticks,
        measure_every=args.measure_every,
        seed=args.seed,
        stepping_mode=args.stepping_mode,
        n_motors=args.n_motors,
        with_xlinks=args.with_xlinks,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
