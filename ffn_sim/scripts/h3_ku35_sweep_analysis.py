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

  fig_h3_ku35_gamma_breakdown.png   (R1 schema only) γ_soft / γ_rigid /
                                    γ_total stacked + rigid fraction with
                                    Luo 2013 1/7 force-sharing reference.

  REPORT_ku35.md                   per-seed plateau γ + ensemble mean +
                                   gate decision (in-band / below band
                                   / above band) + Luo 2013 oracle
                                   comparison if R1 fields present.

Run AFTER all seeds complete:

    python ffn_sim/scripts/h3_ku35_sweep_analysis.py \
        --indir ffn_sim/outputs/h3/production/ku35_canonical

The script tolerates missing seeds (silently skips) so it can be run as
an in-progress preview as well. The γ-breakdown panel is rendered only
when seed JSONs carry the R1 schema fields (tension_{soft,rigid,total}_
mN_per_m). Pre-R1 sweeps populate only tension_mN_per_m (soft-bond only,
no rigid backbone contribution); the breakdown panel is skipped silently.
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
    eng = np.array([d.get("myosin_engaged", d.get("myoss_engaged", 0)) for d in diag])
    bind_total = np.array([d.get("bind_total", 0) for d in diag])
    step_advances = np.array([d.get("step_advances", 0) for d in diag])
    r_over_r0 = np.array([d["r_over_r0"] for d in diag])
    return t_ms, gamma, eng, bind_total, step_advances, r_over_r0


def _per_seed_gamma_breakdown(seed: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Extract R1 γ_{soft,rigid,total} arrays per seed.

    Returns None when the seed JSON lacks the R1 schema fields (pre-R1 sweep —
    only ``tension_mN_per_m`` populated = soft-bond only). Callers use this
    to gate-render the breakdown panel without crashing on legacy JSONs.
    """
    diag = seed["diag"]
    if not diag or "tension_soft_mN_per_m" not in diag[0]:
        return None
    step = np.array([d["step"] for d in diag])
    t_ms = step * seed["dt_s"] * 1e3
    g_soft = np.array([d["tension_soft_mN_per_m"] for d in diag])
    g_rigid = np.array([d["tension_rigid_mN_per_m"] for d in diag])
    g_total = np.array([d["tension_total_mN_per_m"] for d in diag])
    return t_ms, g_soft, g_rigid, g_total


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


def render_gamma_breakdown(seeds: list[dict], out_path: Path) -> dict | None:
    """γ_soft + γ_rigid + γ_total breakdown vs time + rigid-fraction panel.

    Pedagogical for the R1 (RIGID_LAGRANGE_TENSION_DESIGN.md) split: shows
    that for our cortex parameters γ_rigid dominates (~200× soft per the
    R1 motivation §1), which is the structural reason pre-R1 KU-3.5 under-
    reported tension by ~200×. Plots also overlay Luo et al. 2013 (Nature
    Materials 12:1064) WT Dictyostelium force-sharing fit ``ζ = 1/7``
    (myosin II carries ~14% of cortical tension, crosslinkers ~86%) as a
    soft-bond-fraction reference ``γ_soft/γ_total ≈ 6/7 = 0.857``. NB our
    γ_soft includes xlinks + myosin head-actin attach + ERM + myosin
    internal bonds together — to map cleanly onto Luo's ζ we'd need
    finer per-bond-type instrumentation (KU-3.21 candidate, see
    Notion 단계 14c). This panel reports the geometric (soft vs rigid)
    decomposition and overlays Luo's reference for orientation, NOT as
    a strict gate.
    """
    broken = [(sd, _per_seed_gamma_breakdown(sd)) for sd in seeds]
    broken = [(sd, b) for sd, b in broken if b is not None]
    if not broken:
        print("SKIP_GAMMA_BREAKDOWN — no R1 schema fields in any seed JSON.",
              flush=True)
        return None

    fig, axs = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True,
                            gridspec_kw=dict(height_ratios=[2.0, 1.0]))
    a_g, a_frac = axs

    all_t = []; all_soft = []; all_rigid = []; all_total = []
    for sd, (t_ms, g_soft, g_rigid, g_total) in broken:
        a_g.plot(t_ms, g_total, lw=1.0, alpha=0.75,
                 label=f"seed {sd['seed']} γ_total")
        a_g.plot(t_ms, g_soft, lw=0.7, alpha=0.5, ls="--",
                 label=f"seed {sd['seed']} γ_soft")
        a_g.plot(t_ms, g_rigid, lw=0.7, alpha=0.5, ls=":",
                 label=f"seed {sd['seed']} γ_rigid")
        all_t.append(t_ms); all_soft.append(g_soft)
        all_rigid.append(g_rigid); all_total.append(g_total)

    # KU-3.5 band overlay on absolute γ panel.
    a_g.axhspan(0.35, 0.65, color="#a8d8a8", alpha=0.25,
                label="KU-3.5 band [0.35, 0.65] mN/m")
    a_g.set_yscale("symlog", linthresh=1e-4)
    a_g.set_ylabel(r"γ [mN/m]")
    a_g.set_title("KU-3.5 — R1 γ breakdown (soft / rigid / total) + Luo 2013 oracle overlay")
    a_g.grid(alpha=0.25); a_g.legend(fontsize=7, loc="lower right", ncol=2)

    # Bottom panel: γ_soft / γ_total fraction per seed.
    # Luo 2013 WT Dictyostelium: ζ = F_myosin / F_internal = 1/7
    # → crosslinker share (1 - ζ) = 6/7 = 0.857. Our γ_soft includes
    # xlinks + motors + ERM together so this is a PROXY, not a strict
    # gate (motor fraction is bundled inside γ_soft).
    LUO_CROSSLINKER_SHARE = 6.0 / 7.0
    nmin = min(len(t) for t in all_t)
    T = all_t[0][:nmin]
    soft_frac_per_seed = []
    for sd, soft, total in zip([b[0] for b in broken], all_soft, all_total):
        st = soft[:nmin]; to = total[:nmin]
        # Guard against zero γ_total at the very first samples (no motor
        # engagement yet → both γ_soft and γ_rigid ~ 0).
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(np.abs(to) > 1e-20, st / to, np.nan)
        a_frac.plot(T, frac, lw=0.9, alpha=0.65, label=f"seed {sd['seed']}")
        soft_frac_per_seed.append(frac)

    a_frac.axhline(LUO_CROSSLINKER_SHARE, color="#a0522d", ls="--", lw=1.4,
                   label=f"Luo 2013 WT ζ=1/7 → soft-share ≈ {LUO_CROSSLINKER_SHARE:.3f}")
    a_frac.axhline(0.5, color="gray", ls=":", lw=0.6, alpha=0.5,
                   label="50/50 split")
    a_frac.set_ylabel(r"$\gamma_{soft} / \gamma_{total}$")
    a_frac.set_xlabel("sim time $t$ [ms]")
    a_frac.set_ylim(-0.05, 1.10)
    a_frac.grid(alpha=0.25); a_frac.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight"); plt.close(fig)
    print(f"WROTE_GAMMA_BREAKDOWN {out_path}", flush=True)

    # Plateau (last third) summary of soft fraction per seed.
    plateau_slice = slice(2 * nmin // 3, nmin)
    plateau_soft_frac = float(np.nanmean(
        np.stack([f[plateau_slice] for f in soft_frac_per_seed], axis=0)))
    plateau_total = float(np.nanmean(
        np.stack([t[plateau_slice] for t in all_total], axis=0)))
    plateau_soft = float(np.nanmean(
        np.stack([s[plateau_slice] for s in all_soft], axis=0)))
    plateau_rigid = float(np.nanmean(
        np.stack([r[plateau_slice] for r in all_rigid], axis=0)))
    return dict(
        n_seeds=len(broken),
        plateau_soft_mN_per_m=plateau_soft,
        plateau_rigid_mN_per_m=plateau_rigid,
        plateau_total_mN_per_m=plateau_total,
        plateau_soft_fraction=plateau_soft_frac,
        luo_wt_crosslinker_share=LUO_CROSSLINKER_SHARE,
        soft_vs_luo_delta=plateau_soft_frac - LUO_CROSSLINKER_SHARE,
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


def write_report(seeds: list[dict], summary: dict,
                  breakdown: dict | None, out_path: Path) -> None:
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
    if breakdown is not None:
        lines += [
            "",
            "## R1 γ breakdown (soft + rigid + total)",
            "",
            "Per RIGID_LAGRANGE_TENSION_DESIGN.md (PI 2026-05-28 verbal):",
            "the rigid actin backbone Lagrange contribution is now exposed",
            "alongside the soft-bond method-of-planes sum. Plateau (last 1/3):",
            "",
            f"- ⟨γ_soft⟩  = {breakdown['plateau_soft_mN_per_m']:.4e} mN/m",
            f"- ⟨γ_rigid⟩ = {breakdown['plateau_rigid_mN_per_m']:.4e} mN/m",
            f"- ⟨γ_total⟩ = {breakdown['plateau_total_mN_per_m']:.4e} mN/m",
            f"- soft fraction γ_soft/γ_total = {breakdown['plateau_soft_fraction']:.4f}",
            "",
            "### Luo 2013 oracle comparison (Nature Materials 12:1064–1071)",
            "",
            "Luo, Mohan, Iglesias & Robinson 2013 fit Dictyostelium WT cortex to",
            "ζ = F_myosin / F_internal = 1/7 — myosin II carries ~14% of cortical",
            "tension, crosslinkers ~86%. Our γ_soft bundles xlinks + motor-actin",
            "attach + ERM + myosin internal bonds together, so:",
            "",
            f"- Luo expected ‘soft share’ proxy (1 - ζ) = 6/7 = {breakdown['luo_wt_crosslinker_share']:.4f}",
            f"- Our γ_soft/γ_total                          = {breakdown['plateau_soft_fraction']:.4f}",
            f"- Δ (ours − Luo)                              = {breakdown['soft_vs_luo_delta']:+.4f}",
            "",
            "Note: this is a PROXY, not a strict gate. Mapping cleanly to Luo's",
            "ζ requires per-bond-type itemisation of γ_soft (KU-3.21 candidate).",
            "Even so, the magnitude is informative: γ_soft/γ_total ≪ 1 would",
            "indicate the rigid backbone dominates (consistent with the R1",
            "motivation's ~200× under-report claim for soft-only KU-3.5).",
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
    breakdown = render_gamma_breakdown(
        seeds, FIGS / "fig_h3_ku35_gamma_breakdown.png")
    write_report(seeds, summary, breakdown, args.indir / "REPORT_ku35.md")
    print("KU35_ANALYSIS_DONE", flush=True)


if __name__ == "__main__":
    main()
