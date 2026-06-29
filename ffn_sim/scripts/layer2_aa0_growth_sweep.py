"""L2.4: proliferation-driven A/A0(R0) spreading law, fit to A/A0 = a + b/R + c/R² (G3+G4).

The minimal CBM (cohesion + edge-traction) is COHESION-LOCKED — A/A0 ≈ 1 and the PI law does
NOT emerge (L2.3, ``layer2_aa0_sweep.py``, kept as the baseline). This sweep adds the
MECHANISM the minimal model lacked: **contact-inhibited proliferation** (``proliferation.py``).
The experiment runs over days; spreading is partly proliferation-driven, and proliferation is
RIM-localised (contact-inhibited bulk) so the proliferating fraction ∝ surface/volume ∝ 1/R.
That surface-to-volume geometry is the EMERGENT origin of the law's 1/R, 1/R² terms — small
spheroids grow/spread more (relatively) than large ones — NOT a hard-coded term.

Sweeps initial size R0 (via cell count), grows each spheroid for a fixed biological time,
measures the steady spread ratio A/A0, ensemble-averages, and fits the PI's novel form. PI
poster A/A0 is overlay-only at comparison time (we hold only the functional form — never fit
to PI data). Also checks the G4 mechanism gates (rim-localised division, sub-exponential
growth). MCF7 is slow-cycling (~30 h doubling); this is a first-pass CPU harness — production
(more seeds, longer biological time, native counts) runs on the GPU box.

Usage: python -m ffn_sim.scripts.layer2_aa0_growth_sweep
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

from ffn_sim.archive.hoomd_legacy.spheroid.observables import effective_radius
from ffn_sim.archive.hoomd_legacy.spheroid.params import resolve_layer2, resolve_proliferation
from ffn_sim.archive.hoomd_legacy.spheroid.proliferation import run_growth_pooled
from ffn_sim.validation.oracles.spheroid.aa0_law import aa0_model, fit_aa0, term_contributions

_ROOT = Path(__file__).resolve().parents[1]
_CFG = _ROOT / "configs" / "layer2_cbm.yaml"
_OCFG = _ROOT / "validation" / "oracles" / "configs" / "layer2_cbm.yaml"
_FIG_DIR = _ROOT / "outputs" / "layer2" / "figs"
# Per-size checkpoint so a killed long-running sweep resumes instead of restarting (the big
# sizes accumulate wall-time/memory and the sweep can be killed by the environment).
_CKPT = _ROOT / "outputs" / "layer2" / "growth_sweep_core.checkpoint.json"


def _load_ckpt() -> dict:
    if _CKPT.exists():
        return json.loads(_CKPT.read_text())
    return {}


def _save_ckpt(d: dict) -> None:
    _CKPT.write_text(json.dumps(d, indent=2))


def main(argv: list[str] | None = None) -> int:
    import sys as _sys
    args = _sys.argv[1:] if argv is None else argv
    cohesion = "catch" if "--catch" in args else ("morse" if "--morse" in args else "catch")
    cfg = yaml.safe_load(_CFG.read_text())
    resolved = resolve_layer2(cfg)
    prolif = resolve_proliferation(cfg, resolved)
    cad = None
    if cohesion == "catch":
        from ffn_sim.archive.hoomd_legacy.spheroid.cadherin_bonds import resolve_cadherin
        cad = resolve_cadherin(resolved)
    gates = yaml.safe_load(_OCFG.read_text())["spheroid"]["acceptance"]
    g3, g4 = gates["g3"], gates["g4"]
    # cohesion-tagged checkpoint so morse/catch runs don't collide
    global _CKPT
    _CKPT = _ROOT / "outputs" / "layer2" / f"growth_sweep_{cohesion}.checkpoint.json"
    print(f"[growth-sweep] cohesion={cohesion}")

    n_list = [60, 120, 250, 400, 600]
    n_seeds = 3
    total_time = 2.0 * prolif.cycle_time_mean  # 2 MCF7 doublings of biological time (~60 h)
    print(f"[growth-sweep] cycle={prolif.cycle_time_mean/3600:.0f} h, total={total_time/3600:.0f} h "
          f"(~{total_time/prolif.cycle_time_mean:.1f} doublings), n_seeds={n_seeds}, "
          f"min_gap={prolif.min_gap/resolved.morse_r0:.2f}·r0")

    exp_ceiling = 2.0 ** (total_time / prolif.cycle_time_mean)
    ckpt = _load_ckpt()
    for n in n_list:
        key = str(n)
        if key in ckpt:
            row = ckpt[key]
            print(f"  N0={n:4d}  [resumed]  R0={row['R0']*1e6:6.1f} µm  "
                  f"A/A0(core)={row['aac_m']:.2f}±{row['aac_e']:.2f}")
            continue
        r0r, aar, acr, rimr, grr = [], [], [], [], []
        n_ejected = 0
        for s in range(n_seeds):
            res = run_growth_pooled(
                resolved, prolif, n_cells_init=n, total_time=total_time,
                epoch_steps=1200, settle_steps=1000, seed=1000 + s, max_cells=4000,
                cohesion=cohesion, cad=cad,
            )
            if res.get("ejected"):
                n_ejected += 1  # surfaced below (no silent truncation): a runaway realisation
            r0r.append(effective_radius(res["a0"]))
            aar.append(float(res["area_over_a0"][-1]))
            acr.append(float(res["area_core_over_a0"][-1]))  # fragmentation-robust
            if not np.isnan(res["rim_fraction_mean"]):
                rimr.append(res["rim_fraction_mean"])
            grr.append(res["growth_factor"])
        row = {
            "R0": float(np.mean(r0r)),
            "aa_m": float(np.mean(aar)), "aa_e": float(np.std(aar)),
            "aac_m": float(np.mean(acr)), "aac_e": float(np.std(acr)),
            "rim": float(np.mean(rimr)) if rimr else float("nan"),
            "growth": float(np.mean(grr)),
            "sub": bool(np.mean(grr) < exp_ceiling),
            "ejected": n_ejected,
        }
        ckpt[key] = row
        _save_ckpt(ckpt)  # checkpoint after EACH size so a kill resumes here
        ej = f"  ⚠EJECTED {n_ejected}/{n_seeds}" if n_ejected else ""
        print(f"  N0={n:4d}  R0={row['R0']*1e6:6.1f} µm  A/A0(hull)={row['aa_m']:.2f}±{row['aa_e']:.2f}  "
              f"A/A0(core)={row['aac_m']:.2f}±{row['aac_e']:.2f}  growth={row['growth']:.2f}  rim={row['rim']:.2f}{ej}")

    rows = [ckpt[str(n)] for n in n_list]
    R0s = [r["R0"] for r in rows]
    AA0m = [r["aa_m"] for r in rows]; AA0e = [r["aa_e"] for r in rows]
    AAcm = [r["aac_m"] for r in rows]; AAce = [r["aac_e"] for r in rows]
    rims = [r["rim"] for r in rows]; growths = [r["growth"] for r in rows]
    subexp = [r["sub"] for r in rows]
    R0 = np.array(R0s); AA0 = np.array(AA0m); AA0err = np.array(AA0e)
    AAcore = np.array(AAcm); AAcerr = np.array(AAce)
    fit = fit_aa0(R0, AA0)
    fit_core = fit_aa0(R0, AAcore)
    a, b, c, r2 = fit["a"], fit["b"], fit["c"], fit["r_squared"]
    ac, bc, cc, r2c = fit_core["a"], fit_core["b"], fit_core["c"], fit_core["r_squared"]
    print(f"\n[fit hull] A/A0 = {a:.3f} + ({b*1e6:.3f} µm)/R + ({c*1e12:.3f} µm²)/R²   r²={r2:.3f}")
    print(f"[fit core] A/A0 = {ac:.3f} + ({bc*1e6:.3f} µm)/R + ({cc*1e12:.3f} µm²)/R²   r²={r2c:.3f}"
          f"   <-- fragmentation-robust (largest connected component)")
    # the headline fit for the gate is the robust CORE area (hull is fragmentation-inflated)
    fit, r2 = fit_core, r2c
    AA0, AA0err = AAcore, AAcerr
    # edge-vs-bulk readout at the smallest and largest R0
    for R_eval in (R0.min(), R0.max()):
        tc = term_contributions(R_eval, a, b, c)
        print(f"   R0={R_eval*1e6:5.1f} µm: dominant term = {tc['dominant']!r}  "
              f"fractions a/b/c = {tc['fractions']['a']:.2f}/{tc['fractions']['b']:.2f}/"
              f"{tc['fractions']['c']:.2f}")

    # ---- gate verdicts ----
    ok_r2 = r2 >= g3["aa0_fit_r_squared_min"]
    ok_n = len(n_list) >= g3["n_radii_min"]
    # rim gate evaluated on the LARGEST spheroid (a small spheroid is ~all rim → no real bulk)
    rim_big = rims[-1]
    ok_rim = (not np.isnan(rim_big)) and rim_big >= g4["rim_division_fraction_min"]
    ok_sub = all(subexp)
    signal = AA0.max() - AA0.min()
    print(f"\n[G3 fit r² ≥ {g3['aa0_fit_r_squared_min']}] {'PASS' if ok_r2 else 'FAIL'}  (r²={r2:.3f})")
    print(f"[G3 ≥{g3['n_radii_min']} radii        ] {'PASS' if ok_n else 'FAIL'}  ({len(n_list)})")
    print(f"[G4 rim-localised   ≥ {g4['rim_division_fraction_min']}] "
          f"{'PASS' if ok_rim else 'FAIL'}  (largest-spheroid rim frac={rim_big:.2f})")
    print(f"[G4 sub-exponential       ] {'PASS' if ok_sub else 'FAIL'}  (all sizes below 2^(t/τ))")
    print(f"[signal Δ(A/A0)] {signal:.3f}  ({'measurable' if signal > 0.05 else 'small'})")

    try:
        _make_figure(R0, AA0, AA0err, fit, growths, rims, prolif, total_time, signal)
    except Exception as exc:  # noqa: BLE001 — viz is best-effort
        print(f"[viz] skipped figure: {exc}")
    return 0 if (ok_r2 and ok_n and ok_sub) else 1


def _make_figure(R0, AA0, AA0err, fit, growths, rims, prolif, total_time, signal) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    Rg = np.linspace(R0.min() * 0.9, R0.max() * 1.05, 200)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))

    # left: emergent A/A0(R0) + fit
    ax.errorbar(R0 * 1e6, AA0, yerr=AA0err, fmt="o", ms=8, color="seagreen", capsize=4,
                zorder=3, label="CBM + proliferation (mean ± sd)")
    ax.plot(Rg * 1e6, aa0_model(Rg, fit["a"], fit["b"], fit["c"]), "-", color="crimson",
            label=f"fit a+b/R+c/R²  (r²={fit['r_squared']:.3f})")
    ax.axhline(fit["a"], color="gray", ls=":", lw=1, label=f"a (R→∞ baseline) = {fit['a']:.2f}")
    ax.axhline(1.0, color="black", ls="--", lw=0.8, label="A/A₀ = 1 (no spreading)")
    ax.set_xlabel("initial effective radius R₀ (µm)")
    ax.set_ylabel("A / A₀  (spread ratio after growth)")
    ax.set_title(
        f"L2.4 — proliferation-driven MCF7 spreading law\n"
        f"A/A0 = {fit['a']:.2f} + ({fit['b']*1e6:.2f} µm)/R + ({fit['c']*1e12:.1f} µm²)/R²"
        f"   |  Δ(A/A0)={signal:.2f}", fontsize=9,
    )
    ax.legend(fontsize=8)

    # right: the mechanism — growth factor + rim fraction vs R0 (the 1/R surface origin)
    ax2.plot(R0 * 1e6, growths, "s-", color="darkorange", label="population growth factor")
    ax2.axhline(2.0 ** (total_time / prolif.cycle_time_mean), color="darkorange", ls=":",
                lw=1, label=f"exp. ceiling 2^(t/τ) = {2.0**(total_time/prolif.cycle_time_mean):.1f}")
    ax2.plot(R0 * 1e6, rims, "^-", color="purple", label="rim-localised division fraction")
    ax2.axhline(0.70, color="purple", ls=":", lw=1, label="G4 rim ≥ 0.70")
    ax2.set_xlabel("initial effective radius R₀ (µm)")
    ax2.set_ylabel("growth factor  /  rim fraction")
    ax2.set_title("mechanism: smaller R₀ → more relative growth (∝ surface/volume ∝ 1/R)\n"
                  "+ rim-localised, sub-exponential (contact-inhibited bulk)", fontsize=9)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out = _FIG_DIR / "fig_layer2_aa0_growth_law.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
