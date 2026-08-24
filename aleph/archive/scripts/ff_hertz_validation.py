"""Piece #1 — Hertz-inversion validation harness (the gating diagnostic).

Extract an EFFECTIVE apparent Young's modulus E_eff from the FF parallel-plate compression at SMALL
strain by inverting the SAME parabolic (Sneddon) Hertz model an AFM operator uses, so the model is
compared modulus-to-modulus against the measured MCF7 E≈249 Pa (Zbiral 2023 IJMS 24:12208, 10µm-DIAMETER
colloidal bead; KB-6.1.1) — instead of dividing a large-strain nominal stress σ=F/πR² by a small-strain
fixed-E Hertz prediction (the large-strain-vs-small-strain category error flagged in FF_RESULTS_LOG 2026-07-07).

Parabolic Hertz:  F = (4/3)·E*·√R·δ1^1.5,  E* = E/(1−ν²),  ν=0.5 (incompressible cell, E=0.75·E*).
Whole-cell sphere-between-plates: per-contact indent δ1 = δ_total/2 = R_cell·strain; invert with R = R_cell = 7.5 µm.
Report a SLOPE fit of F vs δ1^1.5 over the small-strain sweep (not one point) → E_fit, vs the 224–279 Pa band.

Validity gates (Zbiral δ/R=0.10; we stay an order inside): strain 2–3% ⇒ δ1/R ≤ 0.03. E_pt must be flat vs strain
(a rising E_pt = undrained turgor spike leaking in → not Hertzian → the 249-Pa comparison is invalid).
Run undrained (baseline) and drained (long load_time + K_drained, or pressure_setpoint) for the before/after.
"""
import argparse
import numpy as np

from aleph.laws.gamma_floor import CortexParams, NMIIA_MINIFIL_STALL_PN, build_crosslinked_cortex
from aleph.laws.network_warp import simulate_whole_cell_compression_on_device
from aleph.laws.compartments import resolve_membrane, resolve_nucleus

NU = 0.5
E_MCF7_ZBIRAL = 249.0            # Pa, 10µm colloidal bead whole-cell (KB-6.1.1)
E_MCF7_BAND = (224.0, 279.0)     # Zbiral breast panel: MDA-MB-231 → MCF10A
# Lit-anchored physical values (FF units; see FF_RESULTS_LOG + spec):
LP_EXOSMOTIC = 1.6e-8            # µm/(s·Pa) — COS-7 efflux-rectified (default for compression)
LP_CENTRAL = 1.0e-7             # µm/(s·Pa) — COS-7 + airway Pf
K_DRAINED = 300.0               # Pa — Moeendarbary soft-epithelial drained bulk modulus


def hertz_E_from_force(F_pN: float, R_um: float, delta1_um: float, nu: float = NU) -> float:
    """Invert F = (4/3)·E*·√R·δ1^1.5 → apparent E [Pa]. F [pN], R,δ1 [µm] → E [pN/µm²] = Pa."""
    if delta1_um <= 0.0 or F_pN <= 0.0:
        return 0.0
    e_star = F_pN / ((4.0 / 3.0) * np.sqrt(R_um) * delta1_um ** 1.5)
    return e_star * (1.0 - nu ** 2)


def slope_fit_E(deltas1, forces, R_um, nu=NU):
    """Least-squares slope of F vs δ1^1.5 through the origin → E_fit [Pa]."""
    x = np.asarray(deltas1) ** 1.5
    y = np.asarray(forces)
    m = float(np.sum(x * y) / np.sum(x * x))    # slope through origin
    e_star = m / ((4.0 / 3.0) * np.sqrt(R_um))
    return e_star * (1.0 - nu ** 2)


def _fresh_cortex(n_filaments: int, seed: int = 1):
    cx = build_crosslinked_cortex(CortexParams(), n_filaments=n_filaments, n_xl=n_filaments,
                                  n_myo=n_filaments // 10, rng=np.random.default_rng(seed))
    cx.R0_mean = float(np.linalg.norm(cx.net.pos - cx.net.pos.mean(0), axis=1).mean())
    return cx


def run(strains=(0.02, 0.025, 0.03), n_filaments=2000, n_steps=1600,
        Lp=None, K_drained_Pa=None, load_time_s=None, v_load_um_s=None, xl_koff_per_s=None,
        pressure_setpoint=None, f_excess=0.0, device="cpu", quiet=False):
    """Small-strain Hertz E_fit. v_load_um_s (if set) fixes the loading RATE → per-strain load_time =
    δ_total/v = 2·R·strain/v (physical: reaching larger strain takes longer). Overridden by load_time_s."""
    R0 = _fresh_cortex(n_filaments).R0_mean
    mem = resolve_membrane(f_excess=f_excess) if f_excess > 0.0 else None
    mode = ("setpoint@%.0f" % pressure_setpoint) if pressure_setpoint is not None else (
        "v=%.3gµm/s" % v_load_um_s if v_load_um_s is not None else
        "drained(Lp=%.1e,t=%.1fs)" % (Lp, load_time_s) if Lp is not None and load_time_s is not None else "undrained")
    if not quiet:
        print(f"# Hertz | N_fil={n_filaments} R0={R0:.3f}µm mode={mode} K_drained={K_drained_Pa} "
              f"koff={xl_koff_per_s} f_excess={f_excess}", flush=True)
        print("strain δ1[µm]  F[pN]     dP[Pa] drainf turn  E_pt[Pa]  E_pt/249", flush=True)
    d1s, Fs, rows = [], [], []
    for s in strains:
        lt = load_time_s if load_time_s is not None else (2.0 * R0 * s / v_load_um_s if v_load_um_s else None)
        cx = _fresh_cortex(n_filaments)
        _, m = simulate_whole_cell_compression_on_device(
            cx, NMIIA_MINIFIL_STALL_PN, strain=s, n_steps=n_steps, membrane=mem,
            Lp_um_s_Pa=Lp, K_drained_Pa=K_drained_Pa, load_time_s=lt, xl_koff_per_s=xl_koff_per_s,
            pressure_setpoint=pressure_setpoint, device=device)
        delta1 = R0 * s
        F = m["F_plate_pN"]
        E_pt = hertz_E_from_force(F, R0, delta1)
        d1s.append(delta1); Fs.append(F)
        rows.append({"strain": s, "delta1_um": delta1, "F_pN": F, "E_pt_Pa": E_pt,
                     "dP_Pa": m["dP_turgor_Pa"], "drained_frac": m.get("drained_frac", 0.0)})
        if not quiet:
            print(f"{s*100:5.1f} {delta1:5.3f}  {F:8.1f}  {m['dP_turgor_Pa']:6.0f} {m.get('drained_frac',0):5.2f} "
                  f"{int(m.get('xl_turnover_on',0)):3d}  {E_pt:8.1f}  {E_pt/E_MCF7_ZBIRAL:7.1f}×", flush=True)
    E_fit = slope_fit_E(d1s, Fs, R0)
    flat = (max(r["E_pt_Pa"] for r in rows) / max(min(r["E_pt_Pa"] for r in rows), 1e-9)) if rows else 0
    verdict = "IN BAND" if E_MCF7_BAND[0] <= E_fit <= E_MCF7_BAND[1] else (
        "TOO STIFF" if E_fit > E_MCF7_BAND[1] else "TOO SOFT")
    if not quiet:
        print(f"# E_fit (slope) = {E_fit:.1f} Pa  ({E_fit/E_MCF7_ZBIRAL:.1f}× MCF7)  [{verdict}]  flatness={flat:.2f}", flush=True)
    return {"E_fit_Pa": E_fit, "rows": rows, "verdict": verdict, "R0": R0}


def rate_sweep(rates=(5.0, 0.5, 0.05, 0.005, 0.0005), Lp=LP_CENTRAL, K_drained_Pa=K_DRAINED,
               xl_koff_per_s=0.066, f_excess=0.25, n_filaments=2000, n_steps=1600, device="cpu"):
    """E_fit vs loading rate with FULL physics (drainage + drained solid + α-actinin turnover + membrane
    reservoir). Reveals where the residual closes: fast (Zbiral 5µm/s) → undrained+elastic; slow → drained+remodeled."""
    print(f"# RATE SWEEP (full physics: Lp={Lp:.0e}, K_drained={K_drained_Pa}, koff={xl_koff_per_s}/s, f_excess={f_excess})", flush=True)
    print(f"#   τ_osm≈{(1.0-0.30)* (4/3*3.14159*7.5**3) / (Lp*4*3.14159*7.5**2*5.15e5):.0f}s (drainage), 1/koff≈{1/xl_koff_per_s:.0f}s (turnover)", flush=True)
    print("v[µm/s]  load_t@3%[s]  E_fit[Pa]  ×MCF7   verdict", flush=True)
    out = []
    for v in rates:
        r = run(v_load_um_s=v, Lp=Lp, K_drained_Pa=K_drained_Pa, xl_koff_per_s=xl_koff_per_s,
                f_excess=f_excess, n_filaments=n_filaments, n_steps=n_steps, device=device, quiet=True)
        lt3 = 2.0 * r["R0"] * 0.03 / v
        out.append({"v": v, "E_fit": r["E_fit_Pa"], "verdict": r["verdict"]})
        print(f"{v:7.4g}  {lt3:11.2f}  {r['E_fit_Pa']:9.1f}  {r['E_fit_Pa']/E_MCF7_ZBIRAL:5.1f}×  {r['verdict']}", flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["undrained", "drained", "setpoint", "ratesweep"], default="undrained")
    ap.add_argument("--lp", type=float, default=LP_EXOSMOTIC)
    ap.add_argument("--load-time", type=float, default=300.0, help="physical load time [s] (drained mode)")
    ap.add_argument("--kdrained", type=float, default=K_DRAINED)
    ap.add_argument("--fexcess", type=float, default=0.0)
    ap.add_argument("--nfil", type=int, default=2000)
    ap.add_argument("--steps", type=int, default=1600)
    ap.add_argument("--device", default="cpu")
    a = ap.parse_args()
    kw = dict(n_filaments=a.nfil, n_steps=a.steps, f_excess=a.fexcess, device=a.device)
    if a.mode == "undrained":
        run(**kw)
    elif a.mode == "drained":
        run(Lp=a.lp, K_drained_Pa=a.kdrained, load_time_s=a.load_time, xl_koff_per_s=0.066, **kw)
    elif a.mode == "setpoint":
        run(pressure_setpoint=40.0, K_drained_Pa=a.kdrained, **kw)
    else:
        rate_sweep(Lp=a.lp, K_drained_Pa=a.kdrained, f_excess=a.fexcess, n_filaments=a.nfil,
                   n_steps=a.steps, device=a.device)
    print("# done", flush=True)
