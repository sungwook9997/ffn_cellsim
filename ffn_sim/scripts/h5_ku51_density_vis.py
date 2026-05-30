#!/usr/bin/env python
"""KU-5.1 dendritic density single-/multi-seed analysis + figures.

Reads N ``ku51_*.json`` outputs from
``ffn_sim/outputs/h5/production/ku51/`` (produced by
:mod:`h5_ku51_density`) and renders:

  fig_h5_ku51_density.png       barbed-end density [per μm²] vs sim time,
                                per-seed trajectories + ensemble mean ±1σ
                                band (when N ≥ 2) + KU-5.1 target band.

  fig_h5_ku51_capping.png       n_barbed_ends + n_capped vs sim time
                                (per-seed) — diagnostic for the
                                elongation / branching / capping balance.

  REPORT_ku51.md                per-seed plateau density + ensemble mean +
                                gate decision vs KU-5.1 target ≈ 100/μm².

Run AFTER any seed completes (the driver auto-vizes on each seed end
via subprocess; multi-seed sweeps converge to the full ensemble after the
last seed finishes). The script tolerates missing seeds (silently skips),
so it can be run as an in-progress preview as well.

Acceptance band (PI-ratified 2026-05-30, audit A3):
- KU-5.1 acceptance band = [141, 341] barbed ends/µm², the in-cell
  leading-edge filament-end areal density 241 ± 100/µm² from Abraham
  et al. 1999 (Biophys J 77:1721) ±1σ. This replaces the prior un-ratified
  placeholder ±50% window [50, 150]/µm². The H.5-brief target ≈ 100/µm²
  (Bieling 2016 reconstituted) is drawn as the dashed line; the gate band
  is the Abraham in-cell range. Derivation lives in configs/phase1_h5.yaml
  acceptance block (Magic-Number Block).
- No Bieling F-V or Funk abortive viz here (those are separate gates
  KU-5.2 / KU-5.3 with their own drivers + vis scripts).

Usage:
    python ffn_sim/scripts/h5_ku51_density_vis.py \
        --indir ffn_sim/outputs/h5/production/ku51
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
FIGS = PKG / "outputs" / "h5" / "figs"

# H.5 brief KU-5.1 target (≈ 100 barbed ends per μm² of WAVE plane;
# Bieling 2016 reconstituted) drawn as the dashed reference line. The
# acceptance BAND is the in-cell leading-edge filament-end areal density
# 241 ± 100/µm² (Abraham et al. 1999, Biophys J 77:1721) ±1σ → [141, 341]
# /µm². PI-ratified 2026-05-30 (audit A3), replacing the prior un-ratified
# placeholder [50, 150]. Derivation: configs/phase1_h5.yaml acceptance
# block Magic-Number Block.
KU51_TARGET = 100.0
KU51_BAND = (141.0, 341.0)   # Abraham 1999 241 ± 100 /µm² in-cell ±1σ


def _load_seeds(indir: Path) -> list[dict]:
    seeds = []
    for p in sorted(indir.glob("ku51_*.json")):
        try:
            d = json.loads(p.read_text())
            d["_path"] = str(p)
            seeds.append(d)
        except Exception as e:
            print(f"WARN: skipping {p} ({e})", flush=True)
    return seeds


def _per_seed_arrays(seed: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    diag = seed["diag"]
    step = np.array([d["step"] for d in diag])
    t_ms = step * seed["dt_s"] * 1e3
    density = np.array([d["density_per_um2"] for d in diag])
    n_barbed = np.array([d["n_barbed_ends"] for d in diag])
    n_capped = np.array([d["n_capped"] for d in diag])
    return t_ms, density, n_barbed, n_capped


def render_density(seeds: list[dict], out_path: Path) -> dict:
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.axhspan(*KU51_BAND, color="#a8c8d8", alpha=0.25,
               label=(f"KU-5.1 band {KU51_BAND}/μm² "
                      "(Abraham 1999 241±100, ±1σ)"))
    ax.axhline(KU51_TARGET, color="#306080", ls="--", lw=1.0, alpha=0.7,
               label=f"KU-5.1 target ≈ {KU51_TARGET:.0f}/μm² (H.5 brief)")
    all_t: list[np.ndarray] = []; all_d: list[np.ndarray] = []
    for sd in seeds:
        t_ms, density, *_ = _per_seed_arrays(sd)
        ax.plot(t_ms, density, lw=0.9, alpha=0.65, label=f"seed {sd['seed']}")
        all_t.append(t_ms); all_d.append(density)

    summary: dict = {}
    if seeds:
        nmin = min(len(t) for t in all_t)
        D = np.stack([d[:nmin] for d in all_d], axis=0)
        T = all_t[0][:nmin]
        dm = D.mean(axis=0)
        ds = D.std(axis=0, ddof=1) if D.shape[0] > 1 else np.zeros_like(dm)
        if D.shape[0] > 1:
            ax.plot(T, dm, color="black", lw=2.0,
                    label=f"ensemble mean (n={D.shape[0]})")
            ax.fill_between(T, dm - ds, dm + ds, color="black", alpha=0.15)
        plateau_slice = slice(2 * nmin // 3, nmin)
        plateau_mean = float(dm[plateau_slice].mean())
        plateau_seed_means = D[:, plateau_slice].mean(axis=1)
        in_band = KU51_BAND[0] <= plateau_mean <= KU51_BAND[1]
        verdict = ("PASS" if in_band else
                   "FAIL (below)" if plateau_mean < KU51_BAND[0]
                   else "FAIL (above)")
        ax.text(0.02, 0.97,
                f"plateau (last 1/3) ⟨ρ⟩ = {plateau_mean:.2f}/μm²\n"
                f"per-seed plateau: {[f'{x:.1f}' for x in plateau_seed_means.tolist()]}\n"
                f"verdict (Abraham 1999 band): {verdict}",
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox=dict(boxstyle="round", fc="#fffae0", ec="#888", alpha=0.92))
        summary = dict(
            n_seeds=int(D.shape[0]),
            plateau_mean=plateau_mean,
            plateau_per_seed=plateau_seed_means.tolist(),
            in_band=bool(in_band),
            verdict=verdict,
        )

    ax.set_xlabel("sim time $t$ [ms]")
    ax.set_ylabel(r"dendritic density [barbed ends / μm²]")
    ax.set_title("KU-5.1 — barbed-end density at WAVE plane (lamellipodium production)")
    ax.grid(alpha=0.25); ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight"); plt.close(fig)
    print(f"WROTE_DENSITY {out_path}", flush=True)
    return summary


def render_capping(seeds: list[dict], out_path: Path) -> None:
    fig, axs = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True)
    a_be, a_cap = axs
    for sd in seeds:
        t_ms, _, n_be, n_cap = _per_seed_arrays(sd)
        a_be.plot(t_ms, n_be, lw=1.0, marker=".", ms=2.5,
                  label=f"seed {sd['seed']}")
        a_cap.plot(t_ms, n_cap, lw=1.0, marker=".", ms=2.5,
                   label=f"seed {sd['seed']}")
    a_be.set_ylabel("# eligible barbed ends")
    a_be.set_title("KU-5.1 — elongation / branching / capping balance (all seeds)")
    a_be.grid(alpha=0.25); a_be.legend(fontsize=8, loc="lower right")
    a_cap.set_ylabel("# capped barbed ends")
    a_cap.set_xlabel("sim time $t$ [ms]")
    a_cap.grid(alpha=0.25); a_cap.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight"); plt.close(fig)
    print(f"WROTE_CAPPING {out_path}", flush=True)


def write_report(seeds: list[dict], summary: dict, out_path: Path) -> None:
    lines = [
        "# KU-5.1 dendritic density — REPORT",
        "",
        f"n_seeds: {len(seeds)}",
    ]
    if summary:
        lines += [
            f"plateau ensemble ⟨ρ⟩ = {summary['plateau_mean']:.2f} barbed ends / μm²",
            f"per-seed plateau ρ: {summary['plateau_per_seed']}",
            f"target (H.5 brief, Bieling 2016 reconstituted): ≈ {KU51_TARGET:.0f}/μm²",
            f"acceptance band: {list(KU51_BAND)} barbed ends / μm² "
            "(Abraham et al. 1999 Biophys J 77:1721 in-cell 241±100, ±1σ; "
            "PI-ratified 2026-05-30 audit A3)",
            f"verdict (Abraham 1999 band): **{summary['verdict']}**",
            "",
            "## Per-seed metadata",
        ]
        for sd in seeds:
            lines += [
                f"- seed {sd['seed']}: n_fil={sd['n_fil']}, dt={sd['dt_s']:.3e} s, "
                f"dt_factor={sd['dt_factor']:.5f}, "
                f"n_WAVE={sd['n_WAVE']}, wave_area_m²={sd['wave_area_m2']:.3e}, "
                f"final n_barbed={sd['n_barbed_ends_final']}, "
                f"final n_capped={sd['n_capped_final']}, "
                f"final density={sd['final_density_per_um2']:.2f}/μm²",
            ]
    out_path.write_text("\n".join(lines) + "\n")
    print(f"WROTE_REPORT {out_path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", type=Path,
                    default=PKG / "outputs" / "h5" / "production" / "ku51")
    args = ap.parse_args()
    FIGS.mkdir(parents=True, exist_ok=True)

    seeds = _load_seeds(args.indir)
    if not seeds:
        print("NO_SEEDS — KU-5.1 sweep not yet complete or wrong --indir.",
              flush=True)
        return
    print(f"LOADED {len(seeds)} seeds", flush=True)

    summary = render_density(seeds, FIGS / "fig_h5_ku51_density.png")
    render_capping(seeds, FIGS / "fig_h5_ku51_capping.png")
    write_report(seeds, summary, args.indir / "REPORT_ku51.md")
    print("KU51_ANALYSIS_DONE", flush=True)


if __name__ == "__main__":
    main()
