#!/usr/bin/env python
r"""Apply the pre-registered reading rule to the ``n_inner`` DOWN sweep, and render it.

The rule is NOT chosen here — it was committed in
``aleph/outputs/ac/inner_budget/PREREGISTRATION.md`` before the sweep ran, and this module only
executes it:

  primary observable   ``max_pf_pn`` at the FINAL accepted step, read from the stored per-step
                       trajectory the driver already writes — not a scalar the run chose to keep.
  judgement            against the ACROSS-SEED scatter at ``--outer 40``, never a within-run ``sem``
                       (``STATE.md`` MEASUREMENT RULE: one seed cannot resolve a difference below ~5%).
  INSENSITIVE at p     ⇔ ``|mean(p) − mean(40)| <= S₄₀``, with ``S₄₀`` the sample standard deviation of
                       max\|PF\| across the three seeds at ``--outer 40``.
  reportable result    the LOWEST ``--outer`` still INSENSITIVE, and nothing beyond that.
  gamma                excluded in magnitude AND ratio form — ``STATE.md`` (c) 17.

What the answer may not be turned into is in the pre-registration too, and is repeated on the figure:
insensitivity means 40 is not buying what it costs; it does NOT mean the lower budget converges, and
nothing here converges at all (``STATE.md`` (c) 4).

Reads committed JSON on the dev machine; touches no device.
"""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REFERENCE_OUTER = 40  # the incumbent default, and the point the scatter is taken at


def _final_max_pf(record: dict) -> float:
    """Return max|PF| at the final ACCEPTED step from the stored trajectory."""
    traj = record["measurements"]["trajectory"]
    if not traj:
        raise ValueError("record has an empty trajectory; the run produced no accepted step")
    return float(traj[-1]["max_pf_pn"])


def collect(directory: Path) -> dict[int, list[dict]]:
    """Group every ``outer<N>_seed<S>.json`` record in ``directory`` by its RECORDED ``n_inner``."""
    by_outer: dict[int, list[dict]] = {}
    for path in sorted(directory.glob("outer*_seed*.json")):
        rec = json.loads(path.read_text(encoding="utf-8"))
        # Group on what the run RECORDED, never on the filename — a mislabelled file would otherwise
        # land in the wrong bin silently.
        outer = int(rec["config"]["n_inner"])
        by_outer.setdefault(outer, []).append({
            "path": path.name,
            "seed": int(rec["config"]["seed"]),
            "max_pf_pn": _final_max_pf(rec),
            "steps_run": int(rec["config"]["steps_run"]),
            "wall_s_per_step": float(rec["timing"]["wall_seconds"]) / max(int(rec["timing"]["n_steps"]), 1),
            "cortex_filaments": int(rec["census"]["cortex_filaments"]),
            "build": rec.get("build", {}).get("commit"),
            "device": rec.get("device"),
        })
    return by_outer


def judge(by_outer: dict[int, list[dict]]) -> dict:
    """Execute the pre-registered rule. No threshold is chosen here; S₄₀ comes from the run itself."""
    if REFERENCE_OUTER not in by_outer:
        raise ValueError(f"no --outer {REFERENCE_OUTER} point; the scatter has no reference to come from")
    ref = [r["max_pf_pn"] for r in by_outer[REFERENCE_OUTER]]
    if len(ref) < 2:
        raise ValueError(f"--outer {REFERENCE_OUTER} has {len(ref)} seed(s); a scatter needs at least 2")
    scatter = statistics.stdev(ref)
    ref_mean = statistics.fmean(ref)

    points = {}
    for outer, runs in sorted(by_outer.items()):
        vals = [r["max_pf_pn"] for r in runs]
        mean = statistics.fmean(vals)
        points[outer] = {
            "n_seeds": len(vals),
            "max_pf_mean_pn": mean,
            "max_pf_stdev_pn": statistics.stdev(vals) if len(vals) > 1 else None,
            "abs_shift_from_reference_pn": abs(mean - ref_mean),
            "verdict": "INSENSITIVE" if abs(mean - ref_mean) <= scatter else "SENSITIVE",
            "wall_s_per_step_mean": statistics.fmean([r["wall_s_per_step"] for r in runs]),
            "runs": runs,
        }
    insensitive = [o for o, p in points.items() if p["verdict"] == "INSENSITIVE"]
    return {
        "reference_outer": REFERENCE_OUTER,
        "reference_mean_max_pf_pn": ref_mean,
        "seed_scatter_S40_pn": scatter,
        "points": points,
        "lowest_insensitive_outer": min(insensitive) if insensitive else None,
        "up_control_outer80": points.get(80, {}).get("verdict"),
    }


def render(result: dict, out_png: Path) -> None:
    """Two linear panels from zero: the observable with its scatter band, and the cost it buys."""
    outers = sorted(result["points"])
    means = [result["points"][o]["max_pf_mean_pn"] for o in outers]
    walls = [result["points"][o]["wall_s_per_step_mean"] for o in outers]
    errs = [result["points"][o]["max_pf_stdev_pn"] or 0.0 for o in outers]
    ref, scat = result["reference_mean_max_pf_pn"], result["seed_scatter_S40_pn"]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    ax.axhspan(ref - scat, ref + scat, color="#2c7fb8", alpha=0.16,
               label=f"seed scatter at outer={REFERENCE_OUTER}  (S={scat:.4g} pN)")
    ax.axhline(ref, color="#2c7fb8", lw=1.1, ls="--", label=f"mean at outer={REFERENCE_OUTER}")
    ax.errorbar(outers, means, yerr=errs, marker="o", ms=7, lw=1.4, capsize=4,
                color="#252525", label="3-seed mean ± stdev")
    for o in outers:
        p = result["points"][o]
        ax.annotate(p["verdict"][:5], (o, p["max_pf_mean_pn"]), textcoords="offset points",
                    xytext=(0, 12), ha="center", fontsize=8,
                    color="#2c7fb8" if p["verdict"] == "INSENSITIVE" else "#d94801")
    ax.set_xlabel("n_inner  (driver flag --outer) [iterations per accepted step]")
    ax.set_ylabel("max|PF| at final accepted step [pN]")
    ax.set_title("Primary observable — projection declared before the run")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8.5, loc="best")
    ax.spines[["top", "right"]].set_visible(False)

    ax2.plot(outers, walls, marker="s", ms=7, lw=1.4, color="#d94801")
    ax2.set_xlabel("n_inner  (driver flag --outer) [iterations per accepted step]")
    ax2.set_ylabel("wall-clock per accepted step [s]")
    ax2.set_title("What the budget costs (reported, never argued from)")
    ax2.set_ylim(bottom=0)
    ax2.spines[["top", "right"]].set_visible(False)

    any_run = next(iter(result["points"][outers[0]]["runs"]))
    lowest = result["lowest_insensitive_outer"]
    fig.suptitle(
        "n_inner sensitivity, DOWNWARD — pre-registered "
        "(`outputs/ac/inner_budget/PREREGISTRATION.md`)\n"
        f"{any_run['cortex_filaments']:,} cortical filaments (full native) · {any_run['device']} · "
        f"build {any_run['build']} · 30 accepted steps, dt=0.01 s · 3 seeds per point · "
        f"up-control outer=80: {result['up_control_outer80']}",
        fontsize=10)
    fig.text(0.5, 0.015,
             "INSENSITIVE means 40 is not buying what it costs. It does NOT mean the lower budget "
             "converges — nothing here converges (STATE.md (c) 4). No γ appears, in magnitude or ratio "
             "form (STATE.md (c) 17). kind: diagnostic.",
             ha="center", fontsize=8.3, style="italic", color="#4d4d4d", wrap=True)
    fig.tight_layout(rect=(0, 0.06, 1, 0.88))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=160)
    print(f"wrote {out_png}   lowest INSENSITIVE outer = {lowest}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Apply the pre-registered n_inner reading rule.")
    ap.add_argument("--dir", default="aleph/outputs/ac/inner_budget")
    ap.add_argument("--out-json", default="aleph/outputs/ac/inner_budget/native_record.json")
    ap.add_argument("--out-png", default="aleph/outputs/ac/inner_budget/inner_budget_sensitivity.png")
    args = ap.parse_args()

    by_outer = collect(Path(args.dir))
    result = judge(by_outer)
    record = {
        "schema": "ac.engine.observe/run-record@2",
        "run_label": "inner_budget_sensitivity_down",
        "kind": "diagnostic",
        "evidence": "NATIVE",
        "quantitative_claim_status": "BLOCKED",
        "evidence_basis": (
            "max|PF| at the final accepted step vs n_inner, judged against the across-seed scatter at "
            "the incumbent default. Pre-registered in outputs/ac/inner_budget/PREREGISTRATION.md before "
            "the sweep ran. BLOCKED: nothing converged, so no magnitude here is physical and the "
            "wall-clock is not a speed number for the engine — a timing is a property of an engine at a "
            "stated convergence, and this run states none."
        ),
        "pre_registration": "aleph/outputs/ac/inner_budget/PREREGISTRATION.md",
        "result": result,
    }
    Path(args.out_json).write_text(json.dumps(record, indent=2), encoding="utf-8")
    render(result, Path(args.out_png))
    print(json.dumps({k: v for k, v in result.items() if k != "points"}, indent=2))
    for o in sorted(result["points"]):
        p = result["points"][o]
        print(f"  outer={o:3d}  n={p['n_seeds']}  max|PF| = {p['max_pf_mean_pn']:.6g} pN  "
              f"shift {p['abs_shift_from_reference_pn']:.4g}  {p['verdict']:12s}  "
              f"{p['wall_s_per_step_mean']:.3f} s/step")


if __name__ == "__main__":
    main()
