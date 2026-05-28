#!/usr/bin/env python
"""KU-3.5 canonical multi-seed sweep analysis + figures.

Reads N ``seed{i}.json`` outputs from
``ffn_sim/outputs/h3/production/ku35_canonical/`` (produced by
:mod:`h3_ku35_tension`) and renders:

  fig_h3_ku35_tension_sweep.png   tension γ vs sim time, per-seed
                                  trajectories + ensemble mean ±1σ band
                                  + KU-3.5 target band [0.35, 0.65] mN/m.

  fig_h3_ku35_engagement_sweep.png  myosin engagement + step_advances vs
                                    time + radius contraction r/r_0.

  REPORT_ku35.md                   per-seed plateau γ + ensemble mean +
                                   gate decision (in-band / below band
                                   / above band).

Run AFTER all seeds complete:

    python ffn_sim/scripts/h3_ku35_sweep_analysis.py \
        --indir ffn_sim/outputs/h3/production/ku35_canonical

The script tolerates missing seeds (silently skips) so it can be run as
an in-progress preview as well.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PKG = Path(__file__).resolve().parents[1]
FIGS = PKG / "outputs" / "h3" / "figs"


def _load_seeds(indir: Path) -> list[dict]:
    seeds = []
    for p in sorted(indir.glob("seed*.json")):
        try:
            d = json.loads(p.read_text())
            d["_path"] = str(p)
            seeds.append(d)
        except Exception as e:
            print(f"WARN: skipping {p} ({e})", flush=True)
    return seeds


def _per_seed_arrays(seed: dict) -> tuple[np.ndarray, ...]:
    diag = seed["diag"]
    step = np.array([d["step"] for d in diag])
    t_s = step * seed["dt_s"]
    t_ms = t_s * 1e3
    gamma = np.array([d["tension_mN_per_m"] for d in diag])
    eng = np.array([d["myoss_engaged"] if "myoss_engaged" in d else d.get("myosin_engaged", 0) for d in diag])
    bind_total = np.array([d.get("bind_total", 0) for d in diag])
    step_advances = np.array([d.get("step_advances", 0) for d in diag])
    r_over_r0 = np.array([d["r_over_r0"] for d in diag])
    return t_ms, gamma, eng, bind_total, step_advances, r_over_r0


def render_tension(seeds: list[dict], out_path: Path) -> dict:
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.axhspan(0.35, 0.65, color="#a8d8a8", alpha=0.35,
               label="KU-3.5 band [0.35, 0.65] mN/m")
    ax.axhline(0.5, color="#308030", ls="--", lw=1.0, alpha=0.6,
               label="KU-3.5 target = 0.5 mN/m")
    all_t = []; all_g = []
    for sd in seeds:
        t_ms, gamma, *_ = _per_seed_arrays(sd)
        ax.plot(t_ms, gamma, lw=0.9, alpha=0.6,
                label=f"seed {sd['seed']}")
        all_t.append(t_ms); all_g.append(gamma)
    # Ensemble mean ±1σ on shared time grid (truncate to shortest).
    if seeds:
        nmin = min(len(t) for t in all_t)
        G = np.stack([g[:nmin] for g in all_g], axis=0)
        T = all_t[0][:nmin]
        gm = G.mean(axis=0); gs = G.std(axis=0, ddof=1) if G.shape[0] > 1 else np.zeros_like(gm)
        ax.plot(T, gm, color="black", lw=2.0, label=f"ensemble mean (n={G.shape[0]})")
        ax.fill_between(T, gm - gs, gm + gs, color="black", alpha=0.15)
        # Plateau (last third) statistics.
        plateau_slice = slice(2 * nmin // 3, nmin)
        plateau_mean = float(gm[plateau_slice].mean())
        plateau_seed_means = G[:, plateau_slice].mean(axis=1)
        # Annotate decision.
        in_band = 0.35 <= plateau_mean <= 0.65
        verdict = ("PASS" if in_band else
                   "FAIL (below)" if plateau_mean < 0.35 else "FAIL (above)")
        ax.text(0.02, 0.95,
                f"plateau (last 1/3) ⟨γ⟩ = {plateau_mean:.3e} mN/m\n"
                f"per-seed plateau means: {plateau_seed_means.tolist()}\n"
                f"verdict: {verdict}",
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox=dict(boxstyle="round", fc="#fffae0", ec="#888", alpha=0.92))
    ax.set_yscale("symlog", linthresh=1e-4)
    ax.set_xlabel("sim time $t$ [ms]")
    ax.set_ylabel(r"cortical tension $\gamma$ [mN/m]")
    ax.set_title("KU-3.5 — method-of-planes cortical tension (soft-bond contribution)")
    ax.grid(alpha=0.25); ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight"); plt.close(fig)
    print(f"WROTE_TENSION_SWEEP {out_path}", flush=True)
    if not seeds:
        return {}
    return dict(
        n_seeds=int(len(seeds)),
        plateau_mean=plateau_mean,
        plateau_per_seed=plateau_seed_means.tolist(),
        in_band=bool(in_band),
        verdict=verdict,
    )


def render_engagement(seeds: list[dict], out_path: Path) -> None:
    fig, axs = plt.subplots(3, 1, figsize=(8.5, 9), sharex=True)
    a_eng, a_step, a_r = axs
    for sd in seeds:
        t_ms, _, eng, bind, sa, rr = _per_seed_arrays(sd)
        a_eng.plot(t_ms, eng, lw=1.0, marker=".", ms=2.5,
                   label=f"seed {sd['seed']}")
        a_step.plot(t_ms, sa, lw=1.0, marker=".", ms=2.5,
                    label=f"seed {sd['seed']}")
        a_r.plot(t_ms, rr, lw=1.0, marker=".", ms=2.5,
                 label=f"seed {sd['seed']}")
    a_eng.set_ylabel("# engaged motor-actin\nattach bonds")
    a_eng.set_title("KU-3.5 — myosin engagement + Hill stepping + cortex radius (all seeds)")
    a_eng.grid(alpha=0.25); a_eng.legend(fontsize=8, loc="lower right")

    a_step.set_ylabel("cumulative Hill\nstep advances")
    a_step.grid(alpha=0.25); a_step.legend(fontsize=8, loc="lower right")

    a_r.axhline(1.0, color="gray", ls="--", lw=0.6, alpha=0.5,
                label="r/r₀ = 1 (no contraction)")
    a_r.set_xlabel("sim time $t$ [ms]")
    a_r.set_ylabel(r"$\langle r_{cortex}\rangle / r_0$")
    a_r.grid(alpha=0.25); a_r.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight"); plt.close(fig)
    print(f"WROTE_ENGAGEMENT_SWEEP {out_path}", flush=True)


def write_report(seeds: list[dict], summary: dict, out_path: Path) -> None:
    lines = [
        "# KU-3.5 canonical multi-seed sweep — REPORT",
        "",
        f"n_seeds: {len(seeds)}",
    ]
    if summary:
        lines += [
            f"plateau ensemble ⟨γ⟩ = {summary['plateau_mean']:.4e} mN/m",
            f"per-seed plateau γ: {summary['plateau_per_seed']}",
            f"target band: [0.35, 0.65] mN/m",
            f"verdict: **{summary['verdict']}**",
            "",
            "## Per-seed metadata",
        ]
        for sd in seeds:
            lines += [
                f"- seed {sd['seed']}: n_fil={sd['n_fil']}, dt={sd['dt_s']:.3e} s, "
                f"dt_factor={sd['dt_factor']:.5f}, n_motors={sd['n_motors']}, "
                f"final myosin_engaged={sd['myosin_engaged_final']}, "
                f"step_advances_final={sd['step_advances_final']}, "
                f"r_final/r₀={sd['r_final_over_r0']:.6f}, "
                f"max_drift={sd['max_drift']:.3e}",
            ]
    out_path.write_text("\n".join(lines) + "\n")
    print(f"WROTE_REPORT {out_path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", type=Path,
                    default=PKG / "outputs" / "h3" / "production" / "ku35_canonical")
    args = ap.parse_args()
    FIGS.mkdir(parents=True, exist_ok=True)

    seeds = _load_seeds(args.indir)
    if not seeds:
        print("NO_SEEDS — sweep not yet complete or wrong --indir.", flush=True)
        return
    print(f"LOADED {len(seeds)} seeds", flush=True)

    summary = render_tension(seeds,
                             FIGS / "fig_h3_ku35_tension_sweep.png")
    render_engagement(seeds,
                       FIGS / "fig_h3_ku35_engagement_sweep.png")
    write_report(seeds, summary, args.indir / "REPORT_ku35.md")
    print("KU35_ANALYSIS_DONE", flush=True)


if __name__ == "__main__":
    main()
