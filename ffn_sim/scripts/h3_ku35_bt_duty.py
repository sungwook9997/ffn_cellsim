#!/usr/bin/env python
"""KU-3.5-active TRACK 1 — head bound-fraction (DUTY) vs NMIIA literature.

This is a LIT-ANCHORED override assay (NOT a band-tuning sweep). It asks one
question of the binding-throughput floor:

    Is the model's STEADY-STATE head bound-fraction (occupancy / duty ratio)
    physiological for non-muscle myosin IIA (NMIIA), or is it UNDER-set so that
    raising k_on TO THE LIT-ANCHORED value would raise occupancy → bipolar
    completion → g_soft?

It overrides binding-KINETICS params (head_actin_k_off0, head_actin_k_on) IN
THIS FIXTURE ONLY — never the core modules (cortex/myosin.py, cell/cell.py are
read+reuse). Each override value is DERIVED FROM LITERATURE and cited inline; no
value is chosen to hit the KU-3.5 band (BAND = 0.35-0.65 mN/m).

--------------------------------------------------------------------------------
LITERATURE ANCHOR (the model's own AFINES / Stam-Hocky myosin lineage)
--------------------------------------------------------------------------------
The model re-implements the AFINES motor (CLAUDE.md stack line: "Active filament
model: AFINES ... re-implemented on HOOMD"). AFINES defines the myosin DUTY
RATIO explicitly as

    r_D = k_on^m / (k_on^m + k_end^m)          (Freedman 2017, Biophys J, p8)

with the canonical NMIIA kinetic anchors:

  * NMIIA actin-detachment reference off-rate  k_off(0) = 0.35 s^-1
        — Stam, Hocky et al. 2015, citing Wang 2003 (JBC, NMII kinetics) and
          Kovács 2003 (JBC 278:38132, NMIIA in-vitro motility); reproduced in
          Tam 2021 (bioRxiv 2021.02.23.432588, Table p8/p35 and text p37:
          "the reference off-rate koff(0) for non-muscle myosin is 0.35 s^-1
          (IIA) and 1.71 s^-1 (IIB) ... we adopt the value for myosin-IIA").
  * AFINES attachment rate  k_on^m = 1 s^-1  (range [0.001, 2] s^-1),
        detachment k_off^m = 0.1 s^-1  (Freedman 2017 Table 1).

NMIIA is the canonical HIGH-DUTY-RATIO (processive ensemble) cortical motor; the
ENSEMBLE duty of a minifilament is near-unity, while a SINGLE HEAD at zero load
is low-duty 0.04-0.05 (Harris & Warshaw 1993 JBC 268:14764; the cell-confinement
 refs). The model binds PER HEAD with a two-state attach/detach, so the relevant
comparison is the per-head two-state steady-state bound fraction

    p_bound = k_on / (k_on + k_off0).

--------------------------------------------------------------------------------
MODEL CURRENT VALUE vs LIT (computed in `duty_table`)
--------------------------------------------------------------------------------
  MODEL:  k_on = 50 s^-1, k_off0 = 10 s^-1  →  p_bound = 50/60 = 0.833.
  LIT  :  k_off0(NMIIA) = 0.35 s^-1 (Stam-Hocky / Kovács 2003).
          The MODEL off-rate is ~29x FASTER than the NMIIA reference.
          At the lit off-rate, p_bound is HIGHER, not lower:
            k_on=1  → 0.741 ;  k_on=10 → 0.966 ;  k_on=50 → 0.993.

So the model duty is NOT under-set — it is already 0.83 and the lit-anchored
off-rate pushes it UP toward saturation. Therefore raising k_on cannot raise
duty meaningfully (it is near the equilibrium ceiling). The TRACK-1 hypothesis
"duty is under-set, raise k_on to the lit value to lift it" is FALSIFIED at the
equilibrium level.

HOWEVER there is a second, real, lit-anchored lever the off-rate controls: the
D2-batch CFL clamp in resolve_cortex_myosin sizes batch_dt by
``batch_dt * k_off0 <= 1e-3`` so the per-tick fresh-binding probability is
``p_bind = 1 - exp(-k_on * batch_dt)``. At the model's FAST k_off0=10 the clamp
forces batch_dt <= 1e-4 s → p_bind ~ 0.5% per tick (the per-tick k_on occupancy
THROTTLE = limiter (a) in the GATE-B verdict). The lit-anchored SLOWER off-rate
0.35 s^-1 RELAXES that clamp to batch_dt <= 2.86e-3 s, raising per-tick p_bind
~28x and letting the finite recruitment window fill. This assay measures whether
that lit-anchored off-rate (and the AFINES k_on) raises bound-fraction →
frac_complete_pairs → g_soft, and reports the honest remaining gap.

ACTIVE readout only: g_soft via the method-of-planes over the myosin ATTACH
bonds (``aggregation_ledger`` -> ``tension.g_mop_attach_Npm``). NEVER g_rigid.

Run:
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_bt_duty --duty-table
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_bt_duty --probe --quick
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_bt_duty --probe --fig
"""
from __future__ import annotations

import argparse
import sys
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parents[1]
CFG = PKG / "configs" / "phase1_h3.yaml"

# KU-3.5 cortical-tension band [N/m] (Salbreux 2012; same constant as the gate).
BAND = (0.35e-3, 0.65e-3)

# --------------------------------------------------------------------------- #
# Literature anchors (cited; values DERIVED FROM LITERATURE, not band-tuned).
# --------------------------------------------------------------------------- #
# NMIIA actin-detachment reference off-rate (Stam-Hocky 2015 / Wang 2003 /
# Kovács 2003; via Tam 2021 Table p8/p35, text p37).
K_OFF0_NMIIA = 0.35          # s^-1
# AFINES motor attachment rate (Freedman 2017 Biophys J, Table 1): default 1
# s^-1, range [0.001, 2]. The DUTY-RATIO definition r_D = k_on/(k_on+k_end) is
# Freedman 2017 p8. We use the AFINES default attachment rate.
K_ON_AFINES = 1.0            # s^-1


# --------------------------------------------------------------------------- #
# (0) Pure-arithmetic duty table (HOOMD-free, unit-testable)
# --------------------------------------------------------------------------- #
def steady_state_duty(k_on: float, k_off0: float) -> float:
    """Two-state zero-load steady-state bound fraction p_bound = k_on/(k_on+k_off0).

    This is the model's per-head occupancy ceiling (continuous-time two-state
    detailed balance) and equals the AFINES duty-ratio definition r_D =
    k_on/(k_on + k_end) (Freedman 2017 p8) with k_end ↦ k_off0.
    """
    if k_on + k_off0 <= 0:
        return float("nan")
    return float(k_on / (k_on + k_off0))


def batch_dt_cap(k_off0: float, cfl_product: float = 1.0e-3) -> float:
    """D2-batch CFL cap on batch_dt: batch_dt * k_off0 <= cfl_product.

    Mirrors resolve_cortex_myosin §2 (batch_dt * k_off_max <= 1e-3). The
    per-tick fresh-binding probability p_bind = 1 - exp(-k_on * batch_dt) is
    bounded by this cap; a SLOWER (lit) off-rate RELAXES the cap → higher p_bind
    per tick (attacks the per-tick k_on occupancy throttle, limiter (a)).
    """
    return float(cfl_product / k_off0) if k_off0 > 0 else float("inf")


def duty_table(verbose: bool = True) -> list[dict]:
    """Model vs lit-anchored duty + per-tick p_bind. Pure arithmetic."""
    rows = []
    cases = [
        ("MODEL (yaml)",                       50.0, 10.0),
        ("lit off-rate, model k_on",           50.0, K_OFF0_NMIIA),
        ("lit off-rate, AFINES k_on=1",        K_ON_AFINES, K_OFF0_NMIIA),
        ("lit off-rate, k_on=10",              10.0, K_OFF0_NMIIA),
    ]
    for name, k_on, k_off0 in cases:
        bdt = batch_dt_cap(k_off0)
        p_bind = 1.0 - np.exp(-k_on * bdt)
        rows.append(dict(
            case=name, k_on=k_on, k_off0=k_off0,
            duty=steady_state_duty(k_on, k_off0),
            batch_dt_cap=bdt, p_bind_per_tick=float(p_bind),
        ))
    if verbose:
        print("=== TRACK 1: head bound-fraction (DUTY) vs NMIIA literature ===")
        print("  AFINES duty def: r_D = k_on/(k_on+k_end)  [Freedman 2017 p8]")
        print(f"  NMIIA k_off(0) = {K_OFF0_NMIIA} /s  "
              "[Stam-Hocky 2015 / Wang 2003 / Kovács 2003; via Tam 2021 p37]")
        print(f"  AFINES k_on^m  = {K_ON_AFINES} /s  (range 0.001-2)  "
              "[Freedman 2017 Table 1]")
        print()
        print(f"  {'case':32s} {'k_on':>6s} {'k_off0':>7s} {'duty':>7s} "
              f"{'bdt_cap[s]':>11s} {'p_bind/tick':>12s}")
        for r in rows:
            print(f"  {r['case']:32s} {r['k_on']:6.1f} {r['k_off0']:7.2f} "
                  f"{r['duty']:7.4f} {r['batch_dt_cap']:11.3e} "
                  f"{r['p_bind_per_tick']:12.5f}")
        print()
        print("  VERDICT: model duty 0.833 is NOT under-set; the lit off-rate "
              "(0.35/s) is\n  ~29x SLOWER than the model (10/s) so the lit duty "
              "is HIGHER (-> saturation),\n  not lower. The lit lever is the "
              "RELAXED CFL cap -> higher per-tick p_bind.")
    return rows


# --------------------------------------------------------------------------- #
# (1) Live probe: build the connected-mesh cortex with lit-anchored overrides
# --------------------------------------------------------------------------- #
def _build_cortex(*, k_on: float, k_off0: float, n_fil: int, n_motors: int,
                  seed: int, device: str = "cpu", n_ticks: int = 40,
                  run_steps: int = 0, batch_steps: int | str | None = None):
    """Build the production connected-mesh cortex+myosin with binding-kinetics
    overrides applied THROUGH THE CONFIG DICT (no core edit), then populate the
    bipolar bookkeeping by running myosin recruitment ticks (+ optional short
    integration so attach bonds carry tension for g_soft).

    ``batch_steps``: optional override of the D2-batch tick stride.
      * int  → use that literal value.
      * "cfl_max" → fill the (relaxed-by-lit-off-rate) D2 CFL headroom:
        batch_steps = floor((1e-3/k_off0)/dt). At the lit off-rate (0.35/s) this
        is ~28x the model's effective batch_dt → per-tick p_bind ~28x higher
        (the principled lever on limiter (a), the per-tick k_on throttle).
      * None → yaml default (100).

    Returns (sim, p_cortex, p_myo, topology, myosin_action, batch_dt).
    """
    import yaml
    import hoomd

    from ffn_sim.cortex.cortex import resolve_h3_derived
    from ffn_sim.cortex.myosin import resolve_cortex_myosin
    from ffn_sim.cortex.crosslinkers import resolve_crosslinkers
    from ffn_sim.cell.cell import build_cortex_full_simulation

    cfg = deepcopy(yaml.safe_load(open(CFG)))
    cfg["cortex"]["R_cell"] = 7.5e-6                       # MCF7
    cfg["cortex"]["n_filaments"] = int(n_fil)
    cfg["cortex"]["demo_mode"] = True
    myo = cfg["cortex"]["myosin"]
    myo["n_motors_per_cell"] = int(n_motors)
    myo["stepping_mode"] = "grip_walk"
    myo["mesoscale_force_scaling"] = True
    myo["backbone_length"] = 700.0e-9
    # ---- LIT-ANCHORED binding-KINETICS overrides (config dict, NOT core) ----
    myo["head_actin_k_off0"] = float(k_off0)              # NMIIA 0.35 /s (lit)
    myo["head_actin_k_on"] = float(k_on)                  # AFINES k_on (lit)

    p = resolve_h3_derived(cfg)
    tau_bend = p.gamma_b * p.rest_length ** 3 / p.bending_modulus
    dtc = 0.001 * tau_bend
    # Optional batch_steps override (config dict, NOT core). "cfl_max" fills the
    # lit-off-rate-relaxed D2 CFL headroom so the per-tick p_bind throttle lifts.
    if batch_steps == "cfl_max":
        import math as _m
        myo["batch_steps"] = max(1, int(_m.floor((1.0e-3 / float(k_off0)) / dtc)))
    elif batch_steps is not None:
        myo["batch_steps"] = int(batch_steps)
    p_myo = resolve_cortex_myosin(cfg, dt=dtc, R_cell=p.R_cell)
    p_xl = resolve_crosslinkers(cfg, dt=dtc)
    dev = (hoomd.device.GPU(notice_level=0) if device == "gpu"
           else hoomd.device.CPU(notice_level=0))
    # Soft-start equilibration in the build resolves the LJ cold-start overlap so
    # a subsequent short run is BAOAB-stable at production scale (n_fil=1000).
    # This is numerical hygiene (no physics change); recruitment ticks + the
    # g_soft read happen AFTER it.
    eq_steps = max(800, int(run_steps)) if run_steps > 0 else 0
    hw = build_cortex_full_simulation(
        p, p_xlinks=p_xl, p_myosin=p_myo, device=dev, with_baoab=True,
        rng=np.random.default_rng(seed), connected_mesh=True,
        cm_z_struct=3.7, cm_bundle_mult=2,
        equilibrate=(eq_steps > 0), equilibrate_steps=eq_steps,
        equilibrate_softstart_steps=max(300, eq_steps // 4),
    )
    sim = hw["sim"]
    myo_act = hw.get("myosin_action")
    sim.run(0)
    # Recruitment ticks (de-novo k_on binding) on the static frame — same GATE-B
    # read as h3_ku35_completion_diag (recruitment geometry, integrator not
    # advanced for the binding read).
    if myo_act is not None:
        if getattr(myo_act, "_sim_ref", None) is None:
            myo_act._sim_ref = sim
        for tick in range(int(n_ticks)):
            myo_act.act(tick)
    # Short integration so the attach bonds carry a real tension for the g_soft
    # method-of-planes (g_soft is 0 on the force-free construction frame). The
    # build already soft-start-equilibrated, so this is stable.
    if run_steps > 0:
        try:
            sim.run(int(run_steps))
        except FloatingPointError as exc:        # cold-start instability guard
            print(f"    [WARN] integration FPE ({exc}); reporting recruitment-"
                  "frame g_soft (force-free → ~0).", flush=True)
    return sim, p, p_myo, hw["topology"], myo_act, float(p_myo.batch_dt)


def _head_bound_fraction(myo_act) -> float:
    """Fraction of myosin heads currently bound to actin (occupancy)."""
    bound = getattr(myo_act, "_head_bound_to_actin", None)
    if bound is None or len(bound) == 0:
        return float("nan")
    return float(np.mean(np.asarray(bound) >= 0))


def probe(*, k_on: float, k_off0: float, n_fil: int, n_motors: int,
          seed: int = 1, device: str = "cpu", n_ticks: int = 40,
          run_steps: int = 4000, batch_steps: int | str | None = None,
          verbose: bool = True) -> dict:
    """One condition: build with the (k_on, k_off0) override, measure
    bound-fraction, frac_complete_pairs, and g_soft (mop over attach bonds)."""
    from ffn_sim.scripts.h3_ku35_stresslet import stresslet_ledger
    from ffn_sim.scripts.h3_ku35_aggregation import aggregation_ledger

    t0 = time.time()
    sim, p, p_myo, topology, myo_act, bdt = _build_cortex(
        k_on=k_on, k_off0=k_off0, n_fil=n_fil, n_motors=n_motors, seed=seed,
        device=device, n_ticks=n_ticks, run_steps=run_steps,
        batch_steps=batch_steps)

    bf = _head_bound_fraction(myo_act)
    sled = stresslet_ledger(sim, p_myo=p_myo, myosin_action=myo_act,
                            beads_per_filament=p.beads_per_filament)
    summ = sled["summary"]
    agg = aggregation_ledger(sim, R_cell=p.R_cell, p_myo=p_myo,
                             myosin_action=myo_act)
    g_soft_attach = float(agg["tension"]["g_mop_attach_Npm"])    # ACTIVE readout
    g_ceiling = float(agg["tension"]["g_ceiling_attach_Npm"])
    n_engaged = int(summ["n_engaged_motors"])
    fc = float(summ["frac_complete_pairs"])
    coh = float(summ["coherence"])
    maxF = max((x["max_Fbond_over_Fstall"] for x in sled["stresslets"]
                if np.isfinite(x["max_Fbond_over_Fstall"])), default=float("nan"))
    p_bind_tick = float(1.0 - np.exp(-k_on * bdt))
    row = dict(
        k_on=k_on, k_off0=k_off0, batch_dt=bdt, p_bind_tick=p_bind_tick,
        bound_fraction=bf, n_engaged=n_engaged,
        frac_complete=fc, coherence=coh,
        g_soft_attach_Npm=g_soft_attach,
        g_soft_over_band_lo=g_soft_attach / BAND[0],
        g_ceiling_Npm=g_ceiling,
        max_Fbond_over_Fstall=maxF,
        duty_eq=steady_state_duty(k_on, k_off0),
        wall_s=time.time() - t0,
    )
    if verbose:
        print(f"  [k_on={k_on:>5.1f} k_off0={k_off0:>5.2f}  duty_eq={row['duty_eq']:.3f}"
              f" bdt={bdt:.2e}s p_bind/tick={p_bind_tick:.4f}]  bound_frac={bf*100:5.1f}%"
              f"  engaged={n_engaged:3d}  complete={fc*100:5.1f}%  coh={coh:+.2f}"
              f"  g_soft={g_soft_attach*1e3:.3e}mN/m"
              f" ({row['g_soft_over_band_lo']:.2e}x band_lo)"
              f"  maxF/Fst={maxF:.2f}  [{row['wall_s']:.0f}s]", flush=True)
    del sim
    return row


def run_probe(*, n_fil: int, n_motors: int, seed: int, device: str,
              n_ticks: int, run_steps: int) -> list[dict]:
    """MODEL baseline vs lit-anchored conditions, head-to-head.

    Conditions (batch_steps None = yaml default 100; "cfl_max" = fill the lit
    off-rate's relaxed D2 CFL headroom so per-tick p_bind lifts ~28x):
    """
    conds = [
        ("MODEL (k_on=50, k_off0=10)",                  50.0, 10.0, None),
        ("lit off-rate + model k_on (50, 0.35)",        50.0, K_OFF0_NMIIA, None),
        ("lit off-rate + model k_on + CFL-max bdt",     50.0, K_OFF0_NMIIA, "cfl_max"),
        ("lit off-rate + AFINES k_on (1, 0.35)",        K_ON_AFINES, K_OFF0_NMIIA, None),
        ("lit off-rate + AFINES k_on + CFL-max bdt",    K_ON_AFINES, K_OFF0_NMIIA, "cfl_max"),
    ]
    rows = []
    print("=== TRACK 1 PROBE: bound-fraction / completion / g_soft "
          f"(connected-mesh, n_fil={n_fil}, n_motors={n_motors}, device={device}) ===",
          flush=True)
    for name, k_on, k_off0, bs in conds:
        print(f"\n-- {name} --", flush=True)
        rows.append(probe(k_on=k_on, k_off0=k_off0, n_fil=n_fil,
                          n_motors=n_motors, seed=seed, device=device,
                          n_ticks=n_ticks, run_steps=run_steps, batch_steps=bs))
    return rows


def make_figure(rows: list[dict], out_png: Path) -> None:
    """bound-fraction, frac_complete, g_soft (vs band) across the conditions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [f"k_on={r['k_on']:.0f}\nk_off0={r['k_off0']:.2f}" for r in rows]
    x = np.arange(len(rows))
    bf = [r["bound_fraction"] * 100 for r in rows]
    fc = [r["frac_complete"] * 100 for r in rows]
    gs = [r["g_soft_attach_Npm"] * 1e3 for r in rows]

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    ax[0].bar(x, bf, color="C0")
    ax[0].set_ylabel("head bound-fraction [%]")
    ax[0].set_title("Duty / occupancy")
    ax[0].set_ylim(0, 100)
    ax[1].bar(x, fc, color="C2")
    ax[1].set_ylabel("frac_complete_pairs [%]")
    ax[1].set_title("Bipolar completion")
    ax[1].set_ylim(bottom=0)
    ax[2].bar(x, gs, color="C3")
    ax[2].axhspan(BAND[0] * 1e3, BAND[1] * 1e3, color="0.7", alpha=0.5,
                  label="KU-3.5 band 0.35-0.65 mN/m")
    ax[2].set_ylabel("g_soft (mop over attach bonds) [mN/m]")
    ax[2].set_title("ACTIVE shell tension g_soft")
    ax[2].set_yscale("log")
    ax[2].legend(fontsize=8)
    for a in ax:
        a.set_xticks(x)
        a.set_xticklabels(labels, fontsize=8)
        a.grid(alpha=0.3, axis="y")
    fig.suptitle("KU-3.5-active TRACK 1 — head duty (lit-anchored NMIIA k_off0=0.35/s) "
                 "→ completion → g_soft", fontsize=11)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    print(f"\n[FIG] wrote {out_png}", flush=True)


# --------------------------------------------------------------------------- #
# Self-test (pure-arithmetic duty invariants; HOOMD-free)
# --------------------------------------------------------------------------- #
def self_test(verbose: bool = True) -> bool:
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        if verbose:
            print(f"    [{'PASS' if cond else 'FAIL'}] {name}")

    # Model duty is 50/60 = 0.8333...
    check("model duty == 0.8333",
          abs(steady_state_duty(50.0, 10.0) - 50.0 / 60.0) < 1e-9)
    # Lit off-rate (slower) RAISES duty at fixed k_on (duty is monotone-decreasing
    # in k_off0) → confirms duty is NOT under-set by the model's fast off-rate.
    check("lit off-rate raises duty at k_on=50 (NOT under-set)",
          steady_state_duty(50.0, K_OFF0_NMIIA) > steady_state_duty(50.0, 10.0))
    # AFINES duty-ratio identity r_D = k_on/(k_on+k_end).
    check("AFINES duty identity r_D(1,1)=0.5",
          abs(steady_state_duty(1.0, 1.0) - 0.5) < 1e-12)
    # The lit (slower) off-rate RELAXES the CFL cap → bigger batch_dt → bigger
    # per-tick p_bind (the real lit lever on limiter (a)).
    bdt_model = batch_dt_cap(10.0)
    bdt_lit = batch_dt_cap(K_OFF0_NMIIA)
    check("lit off-rate relaxes batch_dt cap (>= 10x)",
          bdt_lit >= 10.0 * bdt_model)
    p_model = 1 - np.exp(-50.0 * bdt_model)
    p_lit = 1 - np.exp(-50.0 * bdt_lit)
    check("lit off-rate raises per-tick p_bind", p_lit > p_model)
    # Bounds.
    check("duty in [0,1]", 0.0 <= steady_state_duty(1.0, 0.35) <= 1.0)
    if verbose:
        print(f"\n{'='*56}\nTRACK-1 DUTY SELF-TEST "
              f"{'PASSED' if ok else 'FAILED'}\n{'='*56}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--duty-table", action="store_true",
                    help="pure-arithmetic model-vs-lit duty table (no HOOMD)")
    ap.add_argument("--probe", action="store_true",
                    help="live connected-mesh probe: bound-frac/completion/g_soft")
    ap.add_argument("--quick", action="store_true",
                    help="smaller probe (n_fil=300, n_motors=60) for a fast pass")
    ap.add_argument("--n-fil", type=int, default=1000)
    ap.add_argument("--n-motors", type=int, default=100)
    ap.add_argument("--n-ticks", type=int, default=40)
    ap.add_argument("--run-steps", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--device", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--fig", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return 0 if self_test() else 1
    if args.duty_table:
        duty_table()
        return 0
    if args.probe:
        duty_table()
        print()
        n_fil, n_motors = args.n_fil, args.n_motors
        if args.quick:
            n_fil, n_motors = 300, 60
        rows = run_probe(n_fil=n_fil, n_motors=n_motors, seed=args.seed,
                         device=args.device, n_ticks=args.n_ticks,
                         run_steps=args.run_steps)
        print("\n=== TRACK 1 SUMMARY ===", flush=True)
        base = rows[0]
        for r in rows:
            print(f"  k_on={r['k_on']:>5.1f} k_off0={r['k_off0']:>5.2f} "
                  f"p_bind/tick={r['p_bind_tick']:.4f}: "
                  f"bound={r['bound_fraction']*100:5.1f}%  "
                  f"complete={r['frac_complete']*100:5.1f}%  "
                  f"g_soft={r['g_soft_attach_Npm']*1e3:.3e}mN/m "
                  f"({r['g_soft_over_band_lo']:.2e}x band_lo)", flush=True)
        # Best lit-anchored condition (max bound-fraction among non-model rows).
        lit_rows = rows[1:]
        lit = max(lit_rows, key=lambda r: (r["bound_fraction"], r["frac_complete"]))
        d_bf = lit["bound_fraction"] - base["bound_fraction"]
        d_fc = lit["frac_complete"] - base["frac_complete"]
        gs_ratio = (lit["g_soft_attach_Npm"] / base["g_soft_attach_Npm"]
                    if base["g_soft_attach_Npm"] not in (0.0,) else float("inf"))
        print(f"\n  BEST lit-anchored: k_on={lit['k_on']} k_off0={lit['k_off0']} "
              f"p_bind/tick={lit['p_bind_tick']:.4f}", flush=True)
        print(f"  lit-anchored Δbound_fraction = {d_bf*100:+.1f} pts", flush=True)
        print(f"  lit-anchored Δfrac_complete  = {d_fc*100:+.1f} pts", flush=True)
        print(f"  lit-anchored g_soft ratio    = {gs_ratio:.2f}x model", flush=True)
        in_band = BAND[0] <= lit["g_soft_attach_Npm"] <= BAND[1]
        print(f"  lit g_soft in KU-3.5 band? {in_band} "
              f"(still {lit['g_soft_over_band_lo']:.2e}x band_lo)", flush=True)
        if args.fig:
            make_figure(rows, PKG / "outputs" / "h3" / "figs"
                        / "ku35_bt_duty_probe.png")
        return 0

    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
