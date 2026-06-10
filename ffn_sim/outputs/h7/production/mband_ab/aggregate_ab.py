"""Aggregate the M-band-targeting A/B (loop24d) and plot the comparison.

Reads ab_s{SEED}_{on,off}.json produced by run_ab.sh and compares the COHERENT
contractile-traction differential (ON-force − OFF-force, pN) between M-band
targeting (mband=on) and random greedy placement (mband=off), per seed and as
an ensemble. Improvement = M-band raises the contractile traction.

Visualization integrity (CLAUDE.md): per-seed points + ensemble mean±sem,
SI/pN units annotated, no axis truncation, the random-placement arm overlaid as
the reference baseline.
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def _load(mode: str):
    rows = []
    for f in sorted(glob.glob(str(HERE / f"ab_s*_{mode}.json"))):
        d = json.load(open(f))
        rows.append(dict(
            seed=int(Path(f).stem.split("_")[1][1:]),
            diff=float(d["coherent_differential_pN"]),
            sem=float(d["coherent_differential_sem_pN"]),
            eng=float(d["engaged_heads_mean"]),
            off=float(d["coherent_off_traction_pN"]),
            on=float(d["coherent_on_traction_pN"]),
            verdict=d.get("verdict"),
        ))
    return rows


def main() -> int:
    on, off = _load("on"), _load("off")
    if not on or not off:
        print("no A/B jsons yet"); return 1
    seeds = sorted({r["seed"] for r in on} & {r["seed"] for r in off})
    on = {r["seed"]: r for r in on}
    off = {r["seed"]: r for r in off}

    don = np.array([on[s]["diff"] for s in seeds])
    dof = np.array([off[s]["diff"] for s in seeds])
    n = len(seeds)
    on_mean = float(don.mean()); of_mean = float(dof.mean())
    on_sem = float(don.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    of_sem = float(dof.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    # paired difference (same seed → same bundle/thermostat; only placement differs)
    paired = don - dof
    p_mean = float(paired.mean())
    p_sem = float(paired.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    improved = p_mean > 2 * p_sem if p_sem > 0 else None

    print(f"  seeds={seeds}")
    for s in seeds:
        print(f"   seed {s}: ON diff={on[s]['diff']:+8.2f}±{on[s]['sem']:.2f} pN "
              f"(eng {on[s]['eng']:.0f}) | OFF diff={off[s]['diff']:+8.2f}±{off[s]['sem']:.2f} pN "
              f"(eng {off[s]['eng']:.0f}) | paired Δ={on[s]['diff']-off[s]['diff']:+8.2f}")
    print(f"  ENSEMBLE ON  (M-band)  = {on_mean:+.2f} ± {on_sem:.2f} pN")
    print(f"  ENSEMBLE OFF (random)  = {of_mean:+.2f} ± {of_sem:.2f} pN")
    print(f"  PAIRED Δ (ON − OFF)    = {p_mean:+.2f} ± {p_sem:.2f} pN  "
          f"⇒ {'IMPROVEMENT' if improved else ('NO IMPROVEMENT' if improved is False else 'n<2')}")

    summary = dict(
        config="sarcomeric n_fil=36 n_beads=144 n_motors=16 equil=120000 contract=120000 nsamp=20",
        seeds=seeds, n_seeds=n,
        on_diff_pN=[on[s]["diff"] for s in seeds], off_diff_pN=[off[s]["diff"] for s in seeds],
        on_ensemble_mean_pN=on_mean, on_ensemble_sem_pN=on_sem,
        off_ensemble_mean_pN=of_mean, off_ensemble_sem_pN=of_sem,
        paired_delta_mean_pN=p_mean, paired_delta_sem_pN=p_sem,
        improved=bool(improved) if improved is not None else None,
        on_engaged=[on[s]["eng"] for s in seeds], off_engaged=[off[s]["eng"] for s in seeds],
    )
    (HERE / "AB_SUMMARY.json").write_text(json.dumps(summary, indent=2))
    print(f"  json → {HERE/'AB_SUMMARY.json'}")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  (matplotlib unavailable: {e}; skipped figure)"); return 0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    x = np.arange(n)
    w = 0.36
    ax1.bar(x - w / 2, dof, w, yerr=[off[s]["sem"] for s in seeds], capsize=3,
            color="#bbbbbb", label="random (mband off)")
    ax1.bar(x + w / 2, don, w, yerr=[on[s]["sem"] for s in seeds], capsize=3,
            color="#2c7fb8", label="M-band (mband on)")
    ax1.axhline(0, color="k", lw=0.8)
    ax1.set_xticks(x); ax1.set_xticklabels([f"seed {s}" for s in seeds])
    ax1.set_ylabel("coherent contractile differential  ON−OFF force  [pN]")
    ax1.set_title("Per-seed paired same-seed differential\n(+ = contractile traction at FA anchors)")
    ax1.legend(fontsize=8)

    # ensemble
    ax2.errorbar([0], [of_mean], yerr=[of_sem], fmt="o", color="#bbbbbb",
                 capsize=4, ms=9, label="random (off)")
    ax2.errorbar([1], [on_mean], yerr=[on_sem], fmt="o", color="#2c7fb8",
                 capsize=4, ms=9, label="M-band (on)")
    for s in seeds:  # thin per-seed connectors
        ax2.plot([0, 1], [off[s]["diff"], on[s]["diff"]], color="0.7", lw=0.8, alpha=0.7)
    ax2.axhline(0, color="k", lw=0.8)
    ax2.set_xlim(-0.5, 1.5); ax2.set_xticks([0, 1]); ax2.set_xticklabels(["random", "M-band"])
    ax2.set_ylabel("coherent contractile differential [pN]")
    ax2.set_title(f"Ensemble (n={n})  paired Δ = {p_mean:+.1f} ± {p_sem:.1f} pN")
    ax2.legend(fontsize=8)

    fig.suptitle("M-band targeting A/B (loop24d) — per-SF coherent traction, CPU preliminary",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_png = HERE / "mband_ab.png"
    fig.savefig(out_png, dpi=130)
    print(f"  fig → {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
