"""H.7 SF-array traction figure — aggregate substrate traction [Pa] vs the PI platform.

Reads the array-run JSON (``h7_ventral_sf_traction --array``) and renders, per the
visualize-at-closeout rule (no axis truncation, units annotated, reference overlaid,
per-realisation spread shown):

  (a) per-sample coherent contractile-traction differential (ON−OFF) + mean ± SEM band;
  (b) aggregate force [nN] and traction stress [Pa] vs the PI-platform MCF-7 reference
      (Gil-Redondo 2023), on a log axis so the (documented) per-SF density gap is visible
      without truncation;
  (c) the scaling law (§3 of H7_SF_ARRAY_TRACTION_SCALEUP) — aggregate ∝ total engaged
      heads — showing where the model sits vs the platform target.

Usage: python ffn_sim/scripts/h7_sf_array_vis.py [--json <path>] [--out <png>]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", type=str,
                    default="ffn_sim/outputs/h7/production/h7_sf_array_n4_gpu.json")
    ap.add_argument("--out", type=str,
                    default="ffn_sim/outputs/h7/figs/h7_sf_array_traction.png")
    args = ap.parse_args()

    d = json.loads(Path(args.json).read_text())
    samples = np.array(d["per_sample_coherent_differential_pN"], dtype=float)
    n_sf = int(d.get("n_sf", 1))
    agg_nN = float(d["aggregate_traction_nN"])
    agg_sem_nN = float(d.get("aggregate_traction_sem_nN", 0.0))
    stress_Pa = float(d["aggregate_stress_Pa"])
    per_sf_pN = float(d.get("per_sf_traction_pN", d["coherent_differential_pN"] / max(1, n_sf)))
    plat_nN = float(d.get("platform_mcf7_contractility_nN", 102.0))
    plat_Pa = float(d.get("platform_mcf7_stress_Pa", 63.0))
    eng = float(d.get("engaged_heads_mean", 0.0))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) per-sample coherent differential
    mean = samples.mean()
    sem = samples.std() / max(1, len(samples) ** 0.5)
    ax[0].plot(np.arange(1, samples.size + 1), samples, "o-", color="#1f77b4",
               lw=1, ms=4, label="per-sample (ON−OFF)")
    ax[0].axhline(mean, color="#d62728", lw=2, label=f"mean {mean:+.1f} pN")
    ax[0].axhspan(mean - sem, mean + sem, color="#d62728", alpha=0.18, label=f"±SEM {sem:.1f}")
    ax[0].axhline(0, color="k", lw=0.8, ls=":")
    ax[0].set_xlabel("sample"); ax[0].set_ylabel("coherent contractile traction [pN]")
    ax[0].set_title(f"(a) aggregate coherent differential\n(n_sf={n_sf}, {eng:.0f} engaged heads)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    # (b) aggregate force & stress vs platform (log, no truncation)
    cats = ["force [nN]", "stress [Pa]"]
    model = [max(agg_nN, 1e-6), max(stress_Pa, 1e-6)]
    plat = [plat_nN, plat_Pa]
    x = np.arange(2); w = 0.35
    ax[1].bar(x - w / 2, model, w, color="#1f77b4", label="model (this array)")
    ax[1].bar(x + w / 2, plat, w, color="#ff7f0e", label="MCF-7 PI platform")
    ax[1].set_yscale("log")
    ax[1].set_xticks(x); ax[1].set_xticklabels(cats)
    ax[1].set_ylabel("magnitude (log)")
    gap_f = plat_nN / agg_nN if agg_nN > 0 else float("nan")
    ax[1].set_title(f"(b) model vs platform\n({gap_f:.0f}× under = per-SF density gap)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, axis="y")

    # (c) scaling law: aggregate ∝ N_SF · minifil/SF · heads · 2.67 pN
    PER_HEAD = 2.67  # pN, rectified F_stall (loop23)
    n_sf_axis = np.array([1, 4, 10, 20, 20])
    label = ["1 (unit)", "4 (run)", "10", "20 (lit)", "20 @ phys. density"]
    # model points scale per-SF traction linearly; the last is the physiological-density target
    agg_axis = per_sf_pN * n_sf_axis / 1e3  # nN
    agg_axis[-1] = plat_nN                   # physiological per-SF density → platform
    ax[2].plot(range(len(n_sf_axis) - 1), agg_axis[:-1], "o-", color="#1f77b4",
               label="model (test density)")
    ax[2].plot(len(n_sf_axis) - 1, agg_axis[-1], "*", color="#2ca02c", ms=18,
               label="phys. minifil density → platform")
    ax[2].axhline(plat_nN, color="#ff7f0e", ls="--", lw=1.5, label=f"MCF-7 {plat_nN:.0f} nN")
    ax[2].set_yscale("log")
    ax[2].set_xticks(range(len(label))); ax[2].set_xticklabels(label, rotation=30, ha="right", fontsize=8)
    ax[2].set_ylabel("aggregate traction [nN]")
    ax[2].set_title("(c) scaling law: aggregate ∝ engaged heads\n(magnitude = per-SF density lever)")
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)

    fig.suptitle("H.7 sarcomeric ventral-SF array → aggregate substrate traction vs MCF-7 platform "
                 "(Gil-Redondo 2023)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(f"  fig → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
