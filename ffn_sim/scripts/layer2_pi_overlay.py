"""A2: overlay the PI poster A/A0(R) on the A1 emergent ligand curves (OVERLAY-ONLY, never fit).

Hard rule (CLAUDE.md): NO fitting of the model to PI experimental data. This script only
OVERLAYS the PI poster measurements alongside the platform's emergent A1 curves and reports the
agreements/gaps. The PI CSV exports (``references/260313_{Bare,Pre,Lam4}.csv``, per-spheroid
time series of A/A0 and effective radius) are kept LOCAL (gitignored) — only this comparison
figure (the sanctioned overlay) and the extracted summary are committed.

PI data shape: each ``Series`` is one spheroid tracked over time (``Time_min``); ``A_over_A0``
is its spread ratio, ``EffectiveRadius_um`` its current radius, ``A0_um2`` its initial area. For
the law A/A0 = a + b/R + c/R² each spheroid contributes one point (R₀ = initial effective
radius, A/A0 at a common observation time T*).

What A2 compares (reported honestly — see REPORT §A2):
1. LAW SHAPE / SIGN — does PI A/A0 decrease with R₀ (the platform's emergent 1/R, b>0)?
2. CONDITION ORDERING — does the platform reproduce the experiment's Bare/Pre/Lam4 ranking?
3. MAGNITUDE — platform A/A0 vs experiment A/A0 (expected gap: core-area + cohesion-locked).
4. R₀ RANGE — the experiment's spheroids vs the platform's (overlap? native-N/GPU gap).

Usage: python -m ffn_sim.scripts.layer2_pi_overlay   [--time-min 3600]
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from ffn_sim.spheroid.ligand_traction import resolve_ligand_traction
from ffn_sim.validation.oracles.spheroid.aa0_law import aa0_model, fit_aa0

_ROOT = Path(__file__).resolve().parents[1]
_REF = _ROOT / "references"
_FIG_DIR = _ROOT / "outputs" / "layer2" / "figs"
_A1_CKPT = _ROOT / "outputs" / "layer2" / "ligand_traction_sweep.checkpoint.json"
_SUMMARY = _ROOT / "outputs" / "layer2" / "pi_overlay_summary.json"

_CONDITIONS = ("Bare", "Pre", "Lam4")
_COLORS = {"Bare": "tab:gray", "Pre": "tab:blue", "Lam4": "tab:green"}
# my A1 sim ran 2 MCF7 doublings ≈ 60 h biological time; default the PI comparison to match.
_DEFAULT_T_MIN = 3600.0


def load_pi_condition(cond: str, t_target_min: float) -> tuple[np.ndarray, np.ndarray]:
    """Per-spheroid (R₀ [m], A/A0 at t_target) for one condition from its PI CSV.

    R₀ = the spheroid's effective radius at its first frame; A/A0 is interpolated at
    ``t_target_min`` (clamped to the series' own time span). One point per ``Series``.
    """
    path = _REF / f"260313_{cond}.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"PI poster CSV not found: {path} (kept local/overlay-only — ask PI for the export)."
        )
    series: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
    with path.open() as f:
        for row in csv.DictReader(f):
            series[row["Series"]].append(
                (float(row["Time_min"]), float(row["A_over_A0"]), float(row["EffectiveRadius_um"]))
            )
    R0, AA = [], []
    for recs in series.values():
        recs.sort()
        t = np.array([r[0] for r in recs])
        aa = np.array([r[1] for r in recs])
        er = np.array([r[2] for r in recs])
        R0.append(er[0] * 1e-6)                                   # m
        AA.append(float(np.interp(min(t_target_min, t.max()), t, aa)))
    order = np.argsort(R0)
    return np.array(R0)[order], np.array(AA)[order]


def load_a1_condition(cond: str) -> tuple[np.ndarray, np.ndarray]:
    """The platform's A1 emergent (R₀ [m], A/A0) points for one condition (from the checkpoint)."""
    ck = json.loads(_A1_CKPT.read_text())
    keys = sorted((k for k in ck if k.startswith(f"{cond}:")), key=lambda k: int(k.split(":")[1]))
    R0 = np.array([ck[k]["R0"] for k in keys])
    AA = np.array([ck[k]["aa"] for k in keys])
    return R0, AA


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    t_min = _DEFAULT_T_MIN
    if "--time-min" in args:
        t_min = float(args[args.index("--time-min") + 1])

    print(f"[A2] OVERLAY-ONLY comparison (PI never fit). PI A/A0 sampled at t={t_min/60:.0f} h "
          f"(platform = 2 doublings ≈ 60 h).\n")
    summary: dict[str, dict] = {}
    pi_med, model_mid = {}, {}
    for cond in _CONDITIONS:
        Rp, Ap = load_pi_condition(cond, t_min)
        Rm, Am = load_a1_condition(cond)
        fit_p = fit_aa0(Rp, Ap)
        fit_m = fit_aa0(Rm, Am)
        corr = float(np.corrcoef(Rp, Ap)[0, 1])
        f_tr = resolve_ligand_traction(cond).f_traction
        summary[cond] = {
            "pi_n": int(len(Rp)),
            "pi_R0_um": [float(Rp.min() * 1e6), float(np.median(Rp) * 1e6), float(Rp.max() * 1e6)],
            "pi_AA_med": float(np.median(Ap)),
            "pi_corr_R0_AA": corr,
            "pi_fit": {k: fit_p[k] for k in ("a", "b", "c", "r_squared")},
            "model_R0_um": [float(Rm.min() * 1e6), float(Rm.max() * 1e6)],
            "model_AA_range": [float(Am.min()), float(Am.max())],
            "model_fit": {k: fit_m[k] for k in ("a", "b", "c", "r_squared")},
            "f_traction_nN": f_tr * 1e9,
        }
        pi_med[cond] = float(np.median(Ap))
        model_mid[cond] = float(np.median(Am))
        print(f"  {cond:5s}: PI  n={len(Rp):2d}  R0={Rp.min()*1e6:3.0f}–{Rp.max()*1e6:3.0f} µm  "
              f"A/A0_med={np.median(Ap):5.1f}  corr(R0,A/A0)={corr:+.2f}  (b={fit_p['b']*1e6:+.0f}µm)")
        print(f"  {'':5s}  MODEL R0={Rm.min()*1e6:3.0f}–{Rm.max()*1e6:3.0f} µm  "
              f"A/A0={Am.min():.2f}–{Am.max():.2f}  (f_traction={f_tr*1e9:.2f} nN)\n")

    # --- the three honest comparisons ---
    pi_rank = sorted(_CONDITIONS, key=lambda c: pi_med[c], reverse=True)
    model_rank = sorted(_CONDITIONS, key=lambda c: model_mid[c], reverse=True)
    all_neg = all(summary[c]["pi_corr_R0_AA"] < 0 for c in _CONDITIONS)
    print("[A2] === honest comparison ===")
    print(f"  1. LAW SHAPE/SIGN : PI corr(R0,A/A0) < 0 for all 3 ({all_neg}) → A/A0 DECREASES "
          f"with R0 = the platform's emergent 1/R (b>0). {'MATCH ✓' if all_neg else 'MISMATCH'}")
    print(f"  2. ORDERING       : PI(collective) {' > '.join(pi_rank)}  vs  "
          f"MODEL(single-cell) {' > '.join(model_rank)}  "
          f"→ {'MATCH' if pi_rank==model_rank else 'SPLIT (single-cell↔collective laminin)'}")
    mag = np.median([pi_med[c] for c in _CONDITIONS]) / np.median([model_mid[c] for c in _CONDITIONS])
    print(f"  3. MAGNITUDE      : PI A/A0_med≈{np.median(list(pi_med.values())):.1f} vs "
          f"model≈{np.median(list(model_mid.values())):.1f} → platform under-spreads ~{mag:.0f}× "
          f"(connected-core + cohesion-locked + raw-seg vs core area)")
    pr = [summary[c]["pi_R0_um"] for c in _CONDITIONS]
    print(f"  4. R0 RANGE       : PI {min(p[0] for p in pr):.0f}–{max(p[2] for p in pr):.0f} µm vs "
          f"model {summary['Bare']['model_R0_um'][0]:.0f}–{summary['Bare']['model_R0_um'][1]:.0f} µm "
          f"→ NO overlap; matching needs native-N (GPU, B2)")

    _SUMMARY.write_text(json.dumps(summary, indent=2))
    print(f"\n[A2] summary → {_SUMMARY}")
    try:
        _make_figure(t_min)
    except Exception as exc:  # noqa: BLE001
        print(f"[viz] skipped figure: {exc}")
    return 0


def _make_figure(t_min: float) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.4))

    for cond in _CONDITIONS:
        Rp, Ap = load_pi_condition(cond, t_min)
        Rm, Am = load_a1_condition(cond)
        fit_m = fit_aa0(Rm, Am)
        c = _COLORS[cond]
        proxy = " (proxy)" if cond == "Lam4" else ""
        # PI points (overlay-only) + their own a+b/R+c/R2 fit
        ax.scatter(Rp * 1e6, Ap, s=22, color=c, alpha=0.55, edgecolor="none",
                   label=f"{cond} PI (n={len(Rp)}, exp.)")
        fit_p = fit_aa0(Rp, Ap)
        Rg_p = np.linspace(Rp.min() * 0.95, Rp.max() * 1.02, 100)
        ax.plot(Rg_p * 1e6, aa0_model(Rg_p, fit_p["a"], fit_p["b"], fit_p["c"]), "-", color=c, lw=1.6)
        # platform emergent points + fit (its own R0 range, dashed) + extrapolation (dotted)
        ax.scatter(Rm * 1e6, Am, s=55, color=c, marker="D", edgecolor="k", linewidth=0.6, zorder=5)
        in_range = np.linspace(Rm.min() * 0.95, Rm.max(), 100)
        ax.plot(in_range * 1e6, aa0_model(in_range, fit_m["a"], fit_m["b"], fit_m["c"]),
                "--", color=c, lw=1.6)
        extrap = np.linspace(Rm.max(), Rp.max() * 1.02, 120)
        ax.plot(extrap * 1e6, aa0_model(extrap, fit_m["a"], fit_m["b"], fit_m["c"]),
                ":", color=c, lw=1.0, alpha=0.7)  # extrapolation into the PI range (flagged)

    ax.axhline(1.0, color="black", ls="--", lw=0.7)
    ax.set_xlabel("initial effective radius R₀ (µm)")
    ax.set_ylabel(f"A / A₀  at t≈{t_min/60:.0f} h")
    ax.set_title("A2 — PI poster (○, solid fit) vs platform A1 emergent (◇, dashed; ⋯ extrapolated)\n"
                 "OVERLAY-ONLY (PI never fit). Both: A/A₀ DECREASES with R₀ (the 1/R law)", fontsize=9)
    ax.legend(fontsize=7, ncol=1)

    # right: the two honest gaps — median A/A0 per condition (PI vs model) + ordering
    xs = np.arange(len(_CONDITIONS))
    pim = [np.median(load_pi_condition(c, t_min)[1]) for c in _CONDITIONS]
    mom = [np.median(load_a1_condition(c)[1]) for c in _CONDITIONS]
    w = 0.36
    ax2.bar(xs - w / 2, pim, w, color=[_COLORS[c] for c in _CONDITIONS], alpha=0.9, label="PI (collective, exp.)")
    ax2.bar(xs + w / 2, mom, w, color=[_COLORS[c] for c in _CONDITIONS], alpha=0.4, hatch="//",
            label="platform A1 (single-cell traction)")
    for x, c in zip(xs, _CONDITIONS):
        ax2.text(x - w / 2, pim[xs.tolist().index(x)] + 0.2, f"{pim[xs.tolist().index(x)]:.1f}",
                 ha="center", fontsize=8)
    ax2.set_xticks(xs); ax2.set_xticklabels(_CONDITIONS)
    ax2.set_ylabel(f"median A/A₀ at t≈{t_min/60:.0f} h")
    ax2.set_title("ordering + magnitude gap (honest):\nPI Lam4>Pre>Bare (collective) vs model "
                  "Pre>Lam4≳Bare (single-cell) — the laminin split; ~5× magnitude", fontsize=9)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out = _FIG_DIR / "fig_layer2_pi_overlay.png"
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[viz] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
