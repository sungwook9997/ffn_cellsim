#!/usr/bin/env python
"""KU-3.5-active TRACK 2 — completion vs binding-tick window (steady-state test).

The diagnosed KU-3.5-active floor is GENERATION-LIMITED via binding-throughput-
limited bipolar completion: in the connected-mesh cortex only ~18-26 % of ENGAGED
minifilaments form a COMPLETE bipolar pair (a real contractile dipole). The
GATE-B read (``h3_ku35_completion_diag``) populated that number by running the
myosin recruitment Updater for only ``n_ticks = 40`` ticks on the static
construction frame.

This module answers ONE clean question:

    Is 23.5 % completion a SHORT-WINDOW artifact (the binding has not reached
    equilibrium in 40 ticks), or is it the true STEADY STATE at physiological
    binding kinetics?

It does so WITHOUT touching any binding-kinetics parameter — it runs at the
LITERATURE ``k_on`` / ``k_off0`` from ``phase1_h3.yaml`` (no override) and only
varies the number of recruitment ticks the Updater is run for (40 → 200 → 1000 →
4000 → 8000). For each window it records:

  * ``frac_complete_pairs`` (the ACTIVE-readout, via ``h3_ku35_stresslet``),
  * the head BOUND-FRACTION (de-novo occupancy, the upstream driver of
    completion),
  * ``coherence`` over complete pairs (must stay contractile, < 0).

It overlays the ANALYTIC two-state binding equilibrium for an isolated head:

    bound_frac_eq  = k_on / (k_on + k_off0)            (single-head occupancy)
    bound_frac(t)  = bound_frac_eq · (1 − e^{−(k_on+k_off0) t})
    tau_eq         = 1 / (k_on + k_off0)               (equilibration timescale)

so the window can be read directly in units of the binding equilibration time.
A COMPLETE bipolar pair needs BOTH a + head and a − head bound simultaneously on
DIFFERENT antiparallel filaments — the Stam-Hocky / Miyazaki "effective myosin
concentration ∝ (duty ratio)²" picture (Stam 2017 PNAS; Miyazaki 2015) — so the
completion ceiling driven purely by occupancy scales ~ (bound_frac)², further
filtered by the antiparallel-partner geometry + the bipolar VETO. The point of
this script is NOT to predict that product but to show empirically whether the
PLATEAU of frac_complete vs window sits at the 23.5 % floor or materially above
it once binding has equilibrated.

WHY THIS IS NOT BAND-TUNING (Magic-Number Block):
  * NO binding-kinetics parameter is overridden. ``k_on=50/s``, ``k_off0=10/s``,
    ``x_beta``, ``batch_steps`` are read straight from ``phase1_h3.yaml`` — the
    physiological-anchor values (NMII duty cycle; Veigel 2002 / Kovács 2003).
  * The ONLY swept quantity is wall-/sim-time (number of recruitment ticks),
    which is grid-invariant and has a literature reference timescale
    ``tau_eq = 1/(k_on+k_off0)`` it is measured against.
  * Letting binding reach its OWN equilibrium and reading the plateau is the
    OPPOSITE of tuning a knob to hit the KU-3.5 band: if completion plateaus AT
    ~24 %, that is the honest steady-state floor (a NEGATIVE result for "longer
    time fixes it"); if it climbs materially, the 40-tick read was an artifact.
    Either outcome is reported as-is.

Lit anchors (binding kinetics, NOT changed here, only read):
  * ``k_off0 = 10 /s`` head-actin unloaded off-rate — NMII duty cycle anchor,
    Veigel et al. 2002 Nat Cell Biol 4:59 (single-molecule NMII detachment).
  * ``k_on  = 50 /s`` per-head bind rate when an acceptor is present.
  * ⇒ single-head steady-state duty (bound fraction) = 50/(50+10) = 0.833;
    NMII is a HIGH-duty motor (Kovács 2003 JBC 278:38132). ``(duty)² = 0.69``
    is the two-heads-bound-simultaneously upper bound before geometry/veto.
  * ``tau_eq = 1/(k_on+k_off0) = 16.7 ms`` ≈ 239 recruitment ticks at the
    resolved ``batch_dt`` — so the prior 40-tick read = 0.17 tau (NOT
    equilibrated).

ACTIVE readout only (frac_complete_pairs / coherence via stresslet ledger);
NEVER g_rigid. NO core edits, NO shared state, small + short sims.

Run:
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_bt_timescale --self-test
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_bt_timescale --run --quick
  conda run -n ffn_sim python -m ffn_sim.scripts.h3_ku35_bt_timescale --run --fig
"""
from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ffn_sim.common.production_policy import (
    add_production_device_args,
    validate_production_device_args,
)

PKG = Path(__file__).resolve().parents[1]
TRACK = "track2_window"   # track-unique output key


# --------------------------------------------------------------------------- #
# Analytic two-state binding equilibration (HOOMD-free, unit-tested)
# --------------------------------------------------------------------------- #
@dataclass(slots=True)
class BindingEquilibrium:
    """Two-state (bound/unbound) head occupancy equilibration at fixed F=0.

    A single head with on-rate ``k_on`` (acceptor present) and unloaded off-rate
    ``k_off0`` relaxes to its steady-state bound fraction exponentially with the
    equilibration timescale ``tau_eq = 1/(k_on+k_off0)``. This is the upstream
    driver of bipolar completion (completion needs both sides bound).
    """

    k_on: float       # 1/s
    k_off0: float     # 1/s

    @property
    def bound_frac_eq(self) -> float:
        return self.k_on / (self.k_on + self.k_off0)

    @property
    def tau_eq(self) -> float:
        return 1.0 / (self.k_on + self.k_off0)

    def bound_frac(self, t: np.ndarray | float) -> np.ndarray | float:
        """bound_frac(t) starting from zero occupancy."""
        t = np.asarray(t, dtype=np.float64)
        return self.bound_frac_eq * (1.0 - np.exp(-(self.k_on + self.k_off0) * t))

    def ticks_to_fraction(self, frac_of_eq: float, batch_dt: float) -> float:
        """Number of recruitment ticks to reach ``frac_of_eq`` of equilibrium."""
        if not (0.0 < frac_of_eq < 1.0):
            raise ValueError("frac_of_eq must be in (0,1)")
        t = -self.tau_eq * np.log(1.0 - frac_of_eq)
        return float(t / batch_dt)


def equilibrium_from_params(p_myo) -> BindingEquilibrium:
    """Build the analytic equilibrium model from a resolved myosin parameter set
    (reads the LITERATURE k_on/k_off0 — no override)."""
    return BindingEquilibrium(
        k_on=float(p_myo.head_actin_k_on),
        k_off0=float(p_myo.head_actin_k_off0),
    )


# --------------------------------------------------------------------------- #
# Live window sweep on the connected-mesh build (reuses completion_diag builder)
# --------------------------------------------------------------------------- #
def _measure_window(*, n_fil, n_motors, seed, device, n_ticks,
                    allow_cpu_dev: bool = False):
    """Build the connected-mesh cortex + run ``n_ticks`` recruitment ticks at the
    LITERATURE k_on, return (frac_complete, bound_frac, coherence, n_engaged).

    Reuses ``h3_ku35_completion_diag._build_cortex_for_sweep`` (same connected-
    mesh build path as ``mcf7_fullcell_stage1``), with ``reach_scale=1.0`` (no
    reach override) — the ONLY varied quantity is ``n_ticks``.
    """
    from ffn_sim.scripts.h3_ku35_completion_diag import _build_cortex_for_sweep
    from ffn_sim.scripts.h3_ku35_stresslet import stresslet_ledger

    sim, p, p_myo, _topology, myo_act = _build_cortex_for_sweep(
        n_fil=n_fil, n_motors=n_motors, reach_scale=1.0, seed=seed,
        device=device, n_ticks=n_ticks, allow_cpu_dev=allow_cpu_dev,
    )
    sled = stresslet_ledger(
        sim, p_myo=p_myo, myosin_action=myo_act,
        beads_per_filament=p.beads_per_filament)
    summ = sled["summary"]
    bound = np.asarray(myo_act._head_bound_to_actin)
    bound_frac = float((bound >= 0).mean()) if bound.size else float("nan")
    eq = equilibrium_from_params(p_myo)
    out = dict(
        n_ticks=int(n_ticks),
        sim_time_s=float(n_ticks * p_myo.batch_dt),
        n_engaged=int(summ["n_engaged_motors"]),
        frac_complete=float(summ["frac_complete_pairs"]),
        bound_frac=bound_frac,
        coherence=float(summ["coherence"]),
        n_break_total=int(getattr(myo_act, "_n_break_total", 0)),
        k_on=eq.k_on, k_off0=eq.k_off0,
        bound_frac_eq=eq.bound_frac_eq, tau_eq=eq.tau_eq,
        batch_dt=float(p_myo.batch_dt),
        window_over_tau=float(n_ticks * p_myo.batch_dt / eq.tau_eq),
    )
    del sim
    return out


def window_sweep(*, windows, n_fil, n_motors, n_seeds=2, device="gpu",
                 allow_cpu_dev: bool = False, verbose=True):
    """Sweep the recruitment-tick window at the literature k_on; ensemble over
    seeds. Returns per-window aggregated rows (mean ± sd over seeds)."""
    rows = []
    for w in windows:
        per_seed = []
        t0 = time.time()
        for s in range(n_seeds):
            per_seed.append(_measure_window(
                n_fil=n_fil, n_motors=n_motors, seed=1 + s, device=device,
                n_ticks=int(w), allow_cpu_dev=allow_cpu_dev))
        def agg(key):
            vals = np.array([r[key] for r in per_seed], dtype=np.float64)
            vals = vals[np.isfinite(vals)]
            return (float(np.mean(vals)) if vals.size else float("nan"),
                    float(np.std(vals)) if vals.size else float("nan"))
        fc_m, fc_s = agg("frac_complete")
        bf_m, bf_s = agg("bound_frac")
        co_m, co_s = agg("coherence")
        ref = per_seed[0]
        row = dict(
            n_ticks=int(w),
            sim_time_s=ref["sim_time_s"],
            window_over_tau=ref["window_over_tau"],
            n_engaged_mean=float(np.mean([r["n_engaged"] for r in per_seed])),
            frac_complete_mean=fc_m, frac_complete_sd=fc_s,
            bound_frac_mean=bf_m, bound_frac_sd=bf_s,
            coherence_mean=co_m, coherence_sd=co_s,
            n_break_total=int(np.sum([r["n_break_total"] for r in per_seed])),
            bound_frac_eq=ref["bound_frac_eq"], tau_eq=ref["tau_eq"],
            batch_dt=ref["batch_dt"], k_on=ref["k_on"], k_off0=ref["k_off0"],
            wall_s=time.time() - t0,
        )
        rows.append(row)
        if verbose:
            print(
                f"  ticks={w:>5d} (t={row['sim_time_s']*1e3:7.2f} ms ="
                f" {row['window_over_tau']:5.2f}·tau)  "
                f"bound={bf_m*100:5.1f}±{bf_s*100:.1f}%  "
                f"complete={fc_m*100:5.1f}±{fc_s*100:.1f}%  "
                f"coh={co_m:+.2f}  engaged={row['n_engaged_mean']:.0f} "
                f"breaks={row['n_break_total']} [{row['wall_s']:.1f}s]",
                flush=True,
            )
    return rows


def analyse_plateau(rows: list[dict]) -> dict:
    """Decide steady-state-limited vs window-artifact from the window sweep.

    A window-ARTIFACT means completion is still climbing materially at the
    largest window (the 40-tick read undershot). STEADY-STATE-limited means
    completion has plateaued (last two large windows agree within tolerance)
    AND that plateau is near the original 40-tick value.
    """
    if len(rows) < 3:
        return dict(verdict="INSUFFICIENT-WINDOWS")
    rows = sorted(rows, key=lambda r: r["n_ticks"])
    fc = np.array([r["frac_complete_mean"] for r in rows])
    bf = np.array([r["bound_frac_mean"] for r in rows])
    # 40-tick (or smallest) reference vs the plateau (largest two windows).
    fc_small = float(fc[0])
    fc_plateau = float(np.mean(fc[-2:]))
    bf_plateau = float(np.mean(bf[-2:]))
    bf_eq = float(rows[-1]["bound_frac_eq"])
    # plateau test: last two windows agree (relative change small)?
    last = float(fc[-1]); prev = float(fc[-2])
    rel_change_top = abs(last - prev) / max(prev, 1e-9)
    plateaued = rel_change_top < 0.10        # < 10 % change across top decade-step
    # did completion move materially from the 40-tick read to the plateau?
    rel_gain = (fc_plateau - fc_small) / max(fc_small, 1e-9)
    material_gain = rel_gain > 0.25          # > 25 % relative climb = artifact-ish
    # is bound-fraction near its single-head equilibrium at the plateau?
    bf_equilibrated = bf_plateau >= 0.8 * bf_eq

    if material_gain and not plateaued:
        verdict = ("WINDOW-ARTIFACT(still-climbing): completion has NOT plateaued "
                   f"and rose {rel_gain*100:.0f}% from the 40-tick read — longer "
                   "time moves it materially")
    elif material_gain and plateaued:
        verdict = ("WINDOW-SENSITIVE-BUT-PLATEAUED: the 40-tick read undershot "
                   f"by {rel_gain*100:.0f}%, but completion DOES plateau at the "
                   f"equilibrated value ({fc_plateau*100:.1f}%) — that plateau is "
                   "the true steady state")
    else:
        verdict = ("STEADY-STATE-LIMITED: completion plateaus at "
                   f"{fc_plateau*100:.1f}% (≈ the 40-tick read {fc_small*100:.1f}%); "
                   "longer (physiological) time alone does NOT move completion "
                   "materially — the floor is a true steady state, not a window "
                   "artifact")
    return dict(
        verdict=verdict,
        fc_40tick=fc_small,
        fc_plateau=fc_plateau,
        bound_frac_plateau=bf_plateau,
        bound_frac_eq=bf_eq,
        bound_frac_equilibrated=bf_equilibrated,
        rel_gain_small_to_plateau=rel_gain,
        rel_change_top_two=rel_change_top,
        plateaued=plateaued,
        material_gain=material_gain,
    )


# --------------------------------------------------------------------------- #
# Figure (window curve + analytic equilibration overlay) — viz-integrity rules
# --------------------------------------------------------------------------- #
def make_figure(rows: list[dict], out_png: Path, *, title: str = "") -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = sorted(rows, key=lambda r: r["n_ticks"])
    ticks = np.array([r["n_ticks"] for r in rows], dtype=float)
    t_s = np.array([r["sim_time_s"] for r in rows], dtype=float)
    fc = np.array([r["frac_complete_mean"] for r in rows]) * 100
    fc_sd = np.array([r["frac_complete_sd"] for r in rows]) * 100
    bf = np.array([r["bound_frac_mean"] for r in rows]) * 100
    bf_sd = np.array([r["bound_frac_sd"] for r in rows]) * 100

    eq = BindingEquilibrium(k_on=rows[0]["k_on"], k_off0=rows[0]["k_off0"])
    tau = eq.tau_eq
    batch_dt = rows[0]["batch_dt"]
    t_dense = np.linspace(0, max(t_s.max(), 5 * tau), 400)
    bf_analytic = eq.bound_frac(t_dense) * 100

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12.5, 4.8))

    # LEFT: bound fraction vs sim time (ms) — measured vs analytic equilibration.
    ax_l.errorbar(t_s * 1e3, bf, yerr=bf_sd, fmt="s-", color="C0",
                  capsize=3, label="head bound-fraction (measured)")
    ax_l.plot(t_dense * 1e3, bf_analytic, "--", color="C0", alpha=0.6,
              label=f"analytic 1−e^(−t/τ), τ={tau*1e3:.1f} ms")
    ax_l.axhline(eq.bound_frac_eq * 100, color="grey", ls=":",
                 label=f"single-head duty k_on/(k_on+k_off)={eq.bound_frac_eq:.2f}")
    ax_l.axvline(tau * 1e3, color="green", ls=":", alpha=0.5,
                 label=f"τ_eq = {tau*1e3:.1f} ms")
    ax_l.set_xlabel("sim time (ms)  [recruitment ticks × batch_dt]")
    ax_l.set_ylabel("head bound-fraction (%)")
    ax_l.set_ylim(0, 100)
    ax_l.set_title("Binding occupancy equilibration (literature k_on, no override)")
    ax_l.grid(alpha=0.3)
    ax_l.legend(fontsize=7, loc="best")

    # RIGHT: frac_complete_pairs vs ticks (log-x), with 40-tick + plateau marks.
    ax_r.errorbar(ticks, fc, yerr=fc_sd, fmt="o-", color="C3", capsize=3,
                  label="frac_complete_pairs (ACTIVE readout)")
    ax_r.axhline(23.5, color="purple", ls="--", alpha=0.7,
                 label="GATE-B 40-tick read (23.5%)")
    if fc.size >= 2:
        ax_r.axhline(float(np.mean(fc[-2:])), color="black", ls=":",
                     alpha=0.7, label=f"plateau ≈ {np.mean(fc[-2:]):.1f}%")
    ax_r.set_xscale("log")
    ax_r.set_xlabel("recruitment ticks (log)")
    ax_r.set_ylabel("frac_complete_pairs (%)")
    ax_r.set_ylim(0, max(60, fc.max() * 1.3 if fc.size else 60))
    ax_r.set_title("Bipolar completion vs run window")
    ax_r.grid(alpha=0.3, which="both")
    ax_r.legend(fontsize=7, loc="best")

    fig.suptitle(title or "KU-3.5-active TRACK 2: completion vs binding-tick window")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Self-test (analytic equilibration sanity — HOOMD-free)
# --------------------------------------------------------------------------- #
def self_test(verbose: bool = True) -> bool:
    ok = True

    def check(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        if verbose:
            print(f"    [{'PASS' if cond else 'FAIL'}] {name}")

    # (A) equilibrium fraction + timescale match the closed forms exactly.
    eq = BindingEquilibrium(k_on=50.0, k_off0=10.0)
    check("bound_frac_eq = k_on/(k_on+k_off)", abs(eq.bound_frac_eq - 50 / 60) < 1e-12)
    check("tau_eq = 1/(k_on+k_off)", abs(eq.tau_eq - 1.0 / 60.0) < 1e-12)

    # (B) bound_frac(0)=0, bound_frac(∞)→eq, monotone increasing.
    check("bound_frac(0) == 0", abs(float(eq.bound_frac(0.0))) < 1e-12)
    check("bound_frac(large t) → eq",
          abs(float(eq.bound_frac(100 * eq.tau_eq)) - eq.bound_frac_eq) < 1e-6)
    ts = np.linspace(0, 5 * eq.tau_eq, 50)
    bf = eq.bound_frac(ts)
    check("bound_frac monotone increasing", bool(np.all(np.diff(bf) >= -1e-15)))

    # (C) at t = tau, fraction is (1−1/e) of equilibrium (the defining property).
    check("bound_frac(τ) = (1−1/e)·eq",
          abs(float(eq.bound_frac(eq.tau_eq)) - eq.bound_frac_eq * (1 - 1 / np.e))
          < 1e-9)

    # (D) ticks_to_fraction inverts bound_frac (round-trip).
    batch_dt = 6.98e-5
    n = eq.ticks_to_fraction(0.95, batch_dt)
    t = n * batch_dt
    check("ticks_to_fraction(0.95) round-trips to 0.95·eq",
          abs(float(eq.bound_frac(t)) - 0.95 * eq.bound_frac_eq)
          < 1e-3 * eq.bound_frac_eq)
    if verbose:
        print(f"    [info] τ_eq={eq.tau_eq*1e3:.2f} ms; "
              f"40 ticks = {40*batch_dt/eq.tau_eq:.2f}·τ; "
              f"95%-of-eq needs {n:.0f} ticks")

    # (E) analyse_plateau verdict logic on synthetic curves.
    def synth(fcs):
        return [dict(n_ticks=t, frac_complete_mean=f, frac_complete_sd=0.0,
                     bound_frac_mean=b, bound_frac_sd=0.0, bound_frac_eq=0.83,
                     sim_time_s=t * batch_dt, k_on=50.0, k_off0=10.0,
                     batch_dt=batch_dt, tau_eq=eq.tau_eq, window_over_tau=0.0,
                     n_engaged_mean=30, coherence_mean=-0.8, coherence_sd=0.0,
                     n_break_total=0)
                for t, f, b in fcs]
    # flat → steady-state-limited
    flat = synth([(40, 0.235, 0.20), (200, 0.24, 0.55), (1000, 0.235, 0.69),
                  (4000, 0.236, 0.69)])
    rA = analyse_plateau(flat)
    check("flat completion → STEADY-STATE-LIMITED",
          rA["verdict"].startswith("STEADY-STATE-LIMITED"))
    # still climbing → window-artifact
    climb = synth([(40, 0.10, 0.05), (200, 0.25, 0.30), (1000, 0.45, 0.55),
                   (4000, 0.65, 0.69)])
    rB = analyse_plateau(climb)
    check("climbing completion → WINDOW-ARTIFACT",
          rB["verdict"].startswith("WINDOW-ARTIFACT"))

    if verbose:
        print(f"\n{'='*56}\nBT-TIMESCALE SELF-TEST "
              f"{'PASSED' if ok else 'FAILED'}\n{'='*56}")
    return ok


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--run", action="store_true",
                    help="live connected-mesh completion-vs-window sweep "
                         "(literature k_on, no override)")
    ap.add_argument("--quick", action="store_true",
                    help="fewer/smaller windows + smaller cortex for a fast pass")
    ap.add_argument("--n-fil", type=int, default=None)
    ap.add_argument("--n-motors", type=int, default=None)
    ap.add_argument("--n-seeds", type=int, default=2)
    add_production_device_args(ap, default="gpu")
    ap.add_argument("--fig", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return 0 if self_test() else 1

    if args.run:
        validate_production_device_args(ap, args)
        if args.quick:
            windows = [40, 200, 1000, 4000]
            n_fil = args.n_fil or 300
            n_motors = args.n_motors or 40
        else:
            windows = [40, 200, 1000, 4000, 8000]
            n_fil = args.n_fil or 1000
            n_motors = args.n_motors or 100
        print("=== KU-3.5-active TRACK 2: completion vs binding-tick window "
              f"(connected-mesh, device={args.device}) ===", flush=True)
        # report the equilibration timescale up front (literature k_on/k_off).
        from ffn_sim.scripts.h3_ku35_completion_diag import _build_cortex_for_sweep
        sim, p, p_myo, _t, _a = _build_cortex_for_sweep(
            n_fil=n_fil, n_motors=n_motors, reach_scale=1.0, seed=1,
            device=args.device, n_ticks=0,
            allow_cpu_dev=args.allow_cpu_dev)
        eq = equilibrium_from_params(p_myo)
        del sim
        print(f"  k_on={eq.k_on:.0f}/s  k_off0={eq.k_off0:.0f}/s  "
              f"single-head duty (bound_frac_eq)={eq.bound_frac_eq:.3f}\n"
              f"  τ_eq = 1/(k_on+k_off0) = {eq.tau_eq*1e3:.2f} ms = "
              f"{eq.tau_eq/p_myo.batch_dt:.0f} ticks; 40-tick read = "
              f"{40*p_myo.batch_dt/eq.tau_eq:.2f}·τ\n", flush=True)

        rows = window_sweep(
            windows=windows, n_fil=n_fil, n_motors=n_motors,
            n_seeds=args.n_seeds, device=args.device,
            allow_cpu_dev=args.allow_cpu_dev)
        res = analyse_plateau(rows)
        print(f"\n  --- PLATEAU ANALYSIS ---")
        print(f"  40-tick frac_complete = {res['fc_40tick']*100:.1f}%")
        print(f"  plateau frac_complete = {res['fc_plateau']*100:.1f}%")
        print(f"  plateau bound-frac    = {res['bound_frac_plateau']*100:.1f}% "
              f"(single-head eq = {res['bound_frac_eq']*100:.1f}%, "
              f"equilibrated={res['bound_frac_equilibrated']})")
        print(f"  rel gain 40-tick→plateau = {res['rel_gain_small_to_plateau']*100:+.0f}%")
        print(f"\n  VERDICT: {res['verdict']}", flush=True)

        if args.fig:
            out = (PKG / "outputs" / "h3" / "figs"
                   / f"ku35_bt_{TRACK}_completion_vs_window.png")
            make_figure(rows, out,
                        title="KU-3.5-active TRACK 2: completion vs binding-tick "
                              f"window  ({n_fil} fil / {n_motors} motors)")
            print(f"\n[FIG] wrote {out}", flush=True)
        return 0

    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
