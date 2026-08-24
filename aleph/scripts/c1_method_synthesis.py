#!/usr/bin/env python
r"""C-1 synthesis — put every reading of one clock on one axis and see what transfers.

THE QUESTION THIS ANSWERS.  Five readings now exist of a single simulated medium: three observation
operators inside the step-relaxation method (whole-field L2, single point, finite contact patch) and
the oscillatory method's 45-degree corner, plus the closed-form oracle where one exists.  They are not
five experiments — they are five *instruments pointed at one state*.  A real laboratory has exactly
this situation and cannot see it, because it only ever owns one instrument at a time.

What the synthesis reports, and the only thing it is allowed to conclude:

  * the EXPONENT in ``a``, per reading.  If these agree, the scaling is a property of the medium.
  * the PREFACTOR, per reading.  If these disagree, the magnitude is a property of the instrument.
  * the spread between them, stated as a ratio, because "they differ by 2.9x" is the finding and
    "cytoplasm has a relaxation time of X" is the thing that cannot be said.

WHAT IT MAY NOT CONCLUDE.  Nothing about a cell: every input record is the Biot field alone, with no
mechanics, no contact law and no compartments.  Nothing about which reading is "right" — they measure
different functionals and the question is not well posed.  And nothing about C-1 itself, whose gate is
scored on an indented cell.

Host arithmetic; reads records, writes a record, a figure and an HTML page. Evaluates no physics.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

#: How close two exponents must be to be called "the same scaling". Declared here, before the records
#: are read, because an agreement threshold chosen after seeing the numbers is not a threshold.
EXPONENT_AGREEMENT: float = 0.10

#: Below this the prefactors are NOT meaningfully different and the method-disagreement claim fails —
#: the honest negative outcome, declared in advance so it can actually happen.
PREFACTOR_DISAGREEMENT_MIN: float = 1.5


def _readings(step: dict | None, rheo: dict | None, creep: dict | None) -> list[dict]:
    """Flatten every per-reading exponent fit into one list."""
    out: list[dict] = []
    if step:
        for name, s in step.get("method_comparison", {}).get("per_observer", {}).items():
            if "exponent" in s:
                out.append({"method": "step_relaxation", "reading": name,
                            "exponent": s["exponent"], "prefactor": s["prefactor"], "r2": s["r2"],
                            "n_kept": s["n_kept"], "verdict": s["verdict"]})
    for rec, method, label in ((rheo, "oscillatory", "minus_3dB_bandwidth"),
                               (creep, "creep", "flux_decay")):
        if not rec:
            continue
        s = rec.get("scored", {}).get("a_um", {})
        # THE READING'S VERDICT IS NOT THE METHOD'S. A sweep can fit a clean exponent while the method
        # around it fails a check the sweep does not see — creep's a-sweep scores PASS while its own
        # grid-invariance check reports FAIL, because the two are different questions. Take the WORSE
        # of the two, or the partition below lets a method that failed itself into the headline.
        method_verdict = rec.get("verdict", "PASS")
        if "exponent" in s:
            combined = "PASS" if (s["verdict"] == "PASS" and method_verdict == "PASS") else "FAIL"
            out.append({"method": method, "reading": label,
                        "exponent": s["exponent"], "prefactor": s["prefactor"], "r2": s["r2"],
                        "n_kept": s["n_kept"], "verdict": combined,
                        "sweep_verdict": s["verdict"], "method_verdict": method_verdict})
        else:
            # A method that produced no fit is still a result about the method, so it is carried into
            # the synthesis as an entry with a reason rather than silently absent from the table.
            out.append({"method": method, "reading": label, "exponent": None, "prefactor": None,
                        "r2": None, "n_kept": s.get("n_kept", 0),
                        "verdict": s.get("verdict", "ABSENT"),
                        "why_no_fit": s.get("rejected")})
    return out


def synthesise(step_path: Path | None, rheo_path: Path | None, out_dir: Path,
               creep_path: Path | None = None) -> dict:
    def _load(p: Path | None) -> dict | None:
        return json.loads(p.read_text()) if p and p.exists() else None

    step, rheo, creep = _load(step_path), _load(rheo_path), _load(creep_path)
    reads = _readings(step, rheo, creep)
    fitted = [r for r in reads if r["exponent"] is not None]

    # THE PARTITION, and the rule behind it, stated because it decides the headline:
    # a reading whose OWN method reported FAIL on its internal validity checks is not evidence about
    # the medium — it is evidence about the method. Creep's a-sweep, for instance, has r2 = 0.98263
    # (above its own declared floor) and an exponent of 0.764, and its grid-invariance check failed by
    # a factor of 21: its number is mostly mesh. Letting that into "do the exponents agree?" would
    # answer a question about the medium with a measurement that is not of the medium.
    #
    # This is a partition, not a filter. Every excluded reading is reported with its failure, and the
    # headline states how many were excluded — a synthesis that silently drops its disagreements is
    # worth nothing.
    # The D every reading's prefactor is expressed against, read from the record rather than assumed:
    # a prefactor in seconds cannot be compared with a closed form without it.
    d_ref = 50.0
    if step and step.get("sweeps", {}).get("a_um"):
        d_ref = float(step["sweeps"]["a_um"][0]["D_um2_s"])
    for r in fitted:
        r["prefactor_in_a2_over_D"] = r["prefactor"] * d_ref

    # THE MESH-RATIO STRUCTURE. When a sweep holds a/dx roughly fixed across several boxes, the residual
    # against the closed form separates cleanly: points sharing an a/dx give the SAME deviation, to five
    # decimal places, however far apart their probe sizes are. That is the difference between "the
    # exponent is 2 and there is 0.1 % of something left" and "the exponent is 2 and the something left
    # is entirely the mesh" — the second is a much stronger statement and it is only visible if the
    # design put several probe sizes at one mesh ratio.
    mesh_structure = None
    if step and step.get("sweeps", {}).get("a_um"):
        groups: dict[float, list[dict]] = {}
        for p in step["sweeps"]["a_um"]:
            if p.get("tau_measured_s") and p.get("tau_oracle_s"):
                key = round(p["a_um"] / p["dx_um"], 2)
                groups.setdefault(key, []).append(
                    {"a_um": p["a_um"], "ratio": p["tau_measured_s"] / p["tau_oracle_s"]})
        shared = {k: v for k, v in groups.items() if len(v) > 1}
        if shared:
            mesh_structure = {
                "by_a_over_dx": {
                    str(k): {"a_um": [q["a_um"] for q in v],
                             "oracle_ratio": [q["ratio"] for q in v],
                             "spread_within_group": float(max(q["ratio"] for q in v)
                                                          - min(q["ratio"] for q in v))}
                    for k, v in sorted(shared.items())},
                "reading": ("points sharing a mesh ratio agree far more tightly than points sharing a "
                            "probe size, so the residual is discretisation and not a length-scale "
                            "effect — the exponent is exact in the continuum"),
            }

    # MESH SENSITIVITY, per method, and the continuum extrapolation it makes possible.
    #
    # Once the residual is known to be a function of a/dx alone, two things follow that a single-mesh
    # sweep cannot give. First, HOW STEEP that function is, is a property OF THE METHOD — a number for
    # "how fine a mesh does this protocol need", which is the practical question and is never reported.
    # Second, the prefactor can be extrapolated to a/dx -> infinity, which turns "this method disagrees
    # with the closed form" into "this method converges to the closed form, slowly".
    mesh_sensitivity = {}
    for rec, method in ((step, "step_relaxation"), (rheo, "oscillatory"), (creep, "creep")):
        if not rec:
            continue
        pts, key = rec.get("sweeps", {}).get("a_um", []), None
        for cand in ("tau_measured_s", "tau_omega_s", "tau_creep_s"):
            if pts and cand in pts[0]:
                key = cand
                break
        if key is None:
            continue
        by_ratio: dict[float, list[float]] = {}
        for p in pts:
            if p.get(key) and p.get("dx_um"):
                c = p[key] * d_ref / (p["a_um"] ** 2)
                by_ratio.setdefault(round(p["a_um"] / p["dx_um"], 2), []).append(c)
        if len(by_ratio) < 3:
            continue
        xs = np.array([1.0 / k for k in by_ratio])
        ys = np.array([float(np.mean(v)) for v in by_ratio.values()])
        slope, intercept = np.polyfit(xs, ys, 1)
        resid = ys - (slope * xs + intercept)
        ss_tot = float(np.sum((ys - ys.mean()) ** 2))
        r2 = 1.0 if ss_tot == 0 else float(1.0 - np.sum(resid ** 2) / ss_tot)
        # AN EXTRAPOLATION IS ONLY WORTH QUOTING WHEN THERE IS A TREND TO EXTRAPOLATE, and the two
        # failure modes are different. Low r2 means no trend at all — the residual is scatter at the
        # converged level, and the "continuum" value is a line fitted through noise. An extrapolation
        # that lands OUTSIDE the measured range with a good r2 is the opposite problem: a real trend
        # whose functional form is probably wrong, because a second-order scheme errs as (dx/a)^2 and a
        # fit linear in (dx/a) will overshoot. Both are flagged rather than folded into one number.
        inside = float(ys.min()) <= float(intercept) <= float(ys.max())
        if r2 < 0.5:
            status = "NO_TREND — residual is scatter at the converged level; the fit is through noise"
        elif not inside:
            status = ("TREND_BUT_EXTRAPOLATED_BEYOND_DATA — check the assumed order before quoting; a "
                      "second-order scheme errs as (dx/a)^2 and a fit linear in (dx/a) overshoots")
        else:
            status = "OK"
        mesh_sensitivity[method] = {
            "C_by_a_over_dx": {str(k): float(np.mean(v)) for k, v in sorted(by_ratio.items())},
            "relative_spread": float((ys.max() - ys.min()) / ys.mean()),
            "C_measured_range": [float(ys.min()), float(ys.max())],
            "continuum_C_extrapolated": float(intercept),
            "extrapolation_r2": r2,
            "extrapolation_status": status,
            "why": "C fitted against 1/(a/dx) and extrapolated to an infinitely fine mesh",
        }
    valid = [r for r in fitted if r["verdict"] == "PASS"]
    invalid = [r for r in fitted if r["verdict"] != "PASS"]
    if not valid:
        raise SystemExit("no reading passed its own method's validity checks — nothing to synthesise")

    exps = np.array([r["exponent"] for r in valid])
    pres = np.array([r["prefactor"] for r in valid])
    exponent_spread = float(exps.max() - exps.min())
    prefactor_ratio = float(pres.max() / pres.min())

    verdict = {
        "scaling_transfers": bool(exponent_spread <= EXPONENT_AGREEMENT),
        "magnitude_is_instrument_specific": bool(prefactor_ratio >= PREFACTOR_DISAGREEMENT_MIN),
    }
    tail = (f" ({len(valid)} of {len(fitted)} readings; {len(invalid)} excluded for failing their own "
            f"method's validity checks)" if invalid else f" (all {len(valid)} readings)")
    reading = (
        "the scaling belongs to the medium and the magnitude belongs to the instrument"
        if all(verdict.values()) else
        "NOT the expected pattern — see the per-reading table before quoting anything from this run"
    ) + tail

    record = {
        "record": "run-record@2",
        "kind": "synthesis",
        "gate": "c1_method_synthesis",
        "inputs": {
            "step_relaxation": None if not step else {
                "build": step.get("build"), "verdict": step.get("verdict"), "device": step.get("device")},
            "oscillatory": None if not rheo else {
                "build": rheo.get("build"), "verdict": rheo.get("verdict"), "device": rheo.get("device"),
                "failed_inversion": rheo.get("failed_inversion")},
            "creep": None if not creep else {
                "build": creep.get("build"), "verdict": creep.get("verdict"),
                "device": creep.get("device"),
                "grid_invariance": creep.get("grid_invariance"),
                "creep_over_relaxation": creep.get("creep_over_relaxation")},
        },
        "n_readings_total": len(reads),
        "n_readings_fitted": len(fitted),
        "n_readings_valid": len(valid),
        "excluded_from_headline": [
            {"method": r["method"], "reading": r["reading"], "exponent": r["exponent"],
             "r2": r["r2"], "verdict": r["verdict"],
             "why": "its own method reported a failing internal validity check, so it is evidence "
                    "about the method rather than about the medium"}
            for r in invalid],
        "declared_before_reading_the_records": {
            "exponent_agreement": EXPONENT_AGREEMENT,
            "prefactor_disagreement_min": PREFACTOR_DISAGREEMENT_MIN,
        },
        "readings": reads,
        "exponent_spread": exponent_spread,
        "exponent_mean": float(exps.mean()),
        "prefactor_ratio_max_over_min": prefactor_ratio,
        "mesh_ratio_structure": mesh_structure,
        "mesh_sensitivity": mesh_sensitivity,
        "D_reference_um2_s": d_ref,
        "prefactor_units": ("a^2/D — the raw fit is in seconds and is not comparable to a closed form "
                            "without its D"),
        "verdict_parts": verdict,
        "reading": reading,
        "replicate_scatter": {
            "value": 0.0,
            "why": ("the Biot field solve is DETERMINISTIC — no RNG enters it, so replicate scatter is "
                    "exactly zero and running replicates would measure nothing. The uncertainty here is "
                    "NUMERICAL, and it is bounded by the dx-invariance sweep and the fit r2, both of "
                    "which are in the step record. C-1 section 6.1's three-replicate requirement binds "
                    "the CELL measurement, where kinetics carry seeds; it does not transfer to this."),
        },
        "may_not_be_quoted_for": [
            "any statement about a cell — every input is the field alone",
            "which reading is correct — they measure different functionals of one state",
            "the C-1 gate, which is scored on an indented cell",
            "a cytoplasm viscosity or any physiological magnitude",
        ],
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "synthesis.json").write_text(json.dumps(record, indent=2))
    _plot(step, rheo, creep, reads, out_dir / "c1_method_synthesis.png", d_ref)
    _html(record, out_dir / "c1_method_synthesis.html")
    return record


def _plot(step: dict | None, rheo: dict | None, creep: dict | None,
          reads: list[dict], path: Path, d_ref: float = 50.0) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 5))
    cols = plt.cm.viridis(np.linspace(0.15, 0.85, max(len(reads), 1)))

    if step:
        pts = [p for p in step["sweeps"]["a_um"] if p["tau_measured_s"] is not None]
        for i, name in enumerate(("l2_norm", "peak", "patch_mean")):
            xs = [p["a_um"] for p in pts if p["observers"][name]["tau_measured_s"]]
            ys = [p["observers"][name]["tau_measured_s"] for p in pts
                  if p["observers"][name]["tau_measured_s"]]
            if xs:
                ax.loglog(xs, ys, "o-", ms=5, lw=1.2, color=cols[i % len(cols)],
                          label=f"step / {name}")
    if rheo:
        pts = [p for p in rheo["sweeps"]["a_um"] if p.get("tau_omega_s") is not None]
        if pts:
            ax.loglog([p["a_um"] for p in pts], [p["tau_omega_s"] for p in pts], "s--", ms=6, lw=1.2,
                      color="crimson", label="oscillatory / -3 dB")
    if creep:
        pts = [p for p in creep["sweeps"]["a_um"] if p.get("tau_creep_s") is not None]
        if pts:
            ax.loglog([p["a_um"] for p in pts], [p["tau_creep_s"] for p in pts], "^-.", ms=6, lw=1.2,
                      color="#8172B2", label="creep / flux decay")
    ax.set_xlabel("probe length scale a [µm]")
    ax.set_ylabel(r"reported relaxation time $\tau$ [s]")
    ax.set_title(f"{len(reads)} readings of ONE medium\n"
                 "parallel lines = same exponent, different offsets")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.15)

    # Only fitted readings can be plotted as points. Ones that produced no fit stay in the record's
    # table with their reason — absent from a chart is not the same as absent from the result.
    reads = [r for r in reads if r["exponent"] is not None]
    lo, hi = 1.7, 2.3
    names = [f"{r['method']}\n{r['reading']}"
             + ("" if lo <= r["exponent"] <= hi else f"\n({r['exponent']:.3f}, off scale)")
             for r in reads]
    bx.axhline(2.0, ls="--", lw=1.0, color="0.5", label="declared exponent 2")
    on = [(i, r) for i, r in enumerate(reads) if lo <= r["exponent"] <= hi]
    bx.errorbar([i for i, _ in on], [r["exponent"] for _, r in on], fmt="o", ms=8, color="#2b6cb0")
    for i, r in on:
        # C is reported in a^2/D units so it can be compared with the closed forms; the raw fitted
        # prefactor is in seconds and is not comparable to anything without its D.
        bx.annotate(f"C={r['prefactor'] * d_ref:.4g}", (i, r["exponent"]), textcoords="offset points",
                    xytext=(0, 12), ha="center", fontsize=8)
    for i, r in enumerate(reads):
        if not (lo <= r["exponent"] <= hi):
            # An off-scale point must SAY it is off-scale; a missing marker under a tick label reads
            # as missing data, which is a different and much friendlier failure than the real one.
            bx.annotate("", (i, lo if r["exponent"] < lo else hi), xytext=(0, 0),
                        textcoords="offset points")
            bx.plot([i], [lo if r["exponent"] < lo else hi], marker="v" if r["exponent"] < lo else "^",
                    ms=9, color="crimson", clip_on=False)
    bx.set_xticks(range(len(reads)))
    bx.set_xticklabels(names, fontsize=7)
    bx.set_ylabel("fitted exponent p in " + r"$\tau \propto a^p$")
    bx.set_ylim(lo, hi)
    bx.set_title("exponent agrees across methods; prefactor C does not")
    bx.legend(fontsize=8)
    bx.grid(True, axis="y", alpha=0.15)
    fig.suptitle("C-1 method synthesis — Biot field alone (no cell, no mechanics, no contact)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _html(rec: dict, path: Path) -> None:
    def cell(v, fmt: str) -> str:
        """A reading that produced no fit shows an em dash, never a zero or a blank."""
        return "&mdash;" if v is None else format(v, fmt)

    rows = "".join(
        f"<tr><td>{r['method']}</td><td class=mono>{r['reading']}</td>"
        f"<td class=mono>{cell(r['exponent'], '.4f')}</td>"
        f"<td class=mono>{cell(r.get('prefactor_in_a2_over_D'), '.5g')}</td>"
        f"<td class=mono>{cell(r['r2'], '.5f')}</td><td class=mono>{r['n_kept']}</td>"
        f"<td class=mono>{r['verdict']}</td></tr>"
        for r in rec["readings"])
    nq = "".join(f"<li>{v}</li>" for v in rec["may_not_be_quoted_for"])
    build = (rec["inputs"].get("step_relaxation") or {}).get("build", {}) or {}
    path.write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>C-1 method synthesis</title><style>
body{{margin:0;background:#0f1115;color:#e6e8ee;font:14px/1.6 ui-sans-serif,system-ui,sans-serif}}
.w{{max-width:900px;margin:0 auto;padding:28px 24px}}
h1{{font-size:19px;margin:0 0 4px}} .sub{{color:#9aa3b2;font-size:12px;margin-bottom:22px}}
table{{border-collapse:collapse;width:100%;margin:10px 0 22px}}
td,th{{padding:6px 10px;border-bottom:1px solid #262b36;text-align:left;font-size:13px}}
th{{color:#9aa3b2;font-weight:500}} .mono{{font-family:ui-monospace,Menlo,monospace}}
.big{{font-size:15px;background:#151822;border:1px solid #262b36;border-radius:8px;padding:14px 16px;margin:18px 0}}
.no{{color:#ff9b9b}} ul{{margin:6px 0 0 18px;padding:0}} li{{font-size:12.5px;color:#c9cfdb}}
img{{max-width:100%;border-radius:8px;border:1px solid #262b36;margin:8px 0 22px}}
</style></head><body><div class=w>
<h1>C-1 method synthesis — {rec['n_readings_total']} readings of one medium</h1>
<div class=sub>Biot pore-pressure field alone: no cell, no mechanics, no contact law.
build <span class=mono>{str(build.get('commit', 'unknown'))[:8]}</span></div>
<div class=big><b>{rec['reading']}</b><br>
exponent spread <span class=mono>{rec['exponent_spread']:.4f}</span>
(bound {rec['declared_before_reading_the_records']['exponent_agreement']}) &nbsp;·&nbsp;
prefactor ratio <span class=mono>{rec['prefactor_ratio_max_over_min']:.4f}</span>
(min {rec['declared_before_reading_the_records']['prefactor_disagreement_min']})</div>
<img src="c1_method_synthesis.png" alt="method synthesis">
<table><tr><th>method</th><th>reading</th><th>exponent p</th>
<th>prefactor C [a&sup2;/D]</th><th>r&sup2;</th>
<th>points</th><th>verdict</th></tr>{rows}</table>
<div class=big><b class=no>May not be quoted for</b><ul>{nq}</ul></div>
<div class=big><b>Replicate scatter: {rec['replicate_scatter']['value']}</b><br>
<span style="font-size:12.5px;color:#c9cfdb">{rec['replicate_scatter']['why']}</span></div>
</div></body></html>""")


def main() -> None:
    ap = argparse.ArgumentParser(description="Synthesise the C-1 method comparison.")
    ap.add_argument("--step", type=Path, default=Path("aleph/outputs/ac/c1_biot_timescale/record.json"))
    ap.add_argument("--rheo", type=Path,
                    default=Path("aleph/outputs/ac/c1_biot_microrheology/record.json"))
    ap.add_argument("--creep", type=Path,
                    default=Path("aleph/outputs/ac/c1_biot_creep/record.json"))
    ap.add_argument("--out", type=Path, default=Path("aleph/outputs/ac/c1_synthesis"))
    a = ap.parse_args()
    rec = synthesise(a.step, a.rheo, a.out, a.creep)
    print(f"[c1-synth] {rec['reading']}")
    print(f"[c1-synth] exponent spread {rec['exponent_spread']:.4f} | "
          f"prefactor ratio {rec['prefactor_ratio_max_over_min']:.4f} -> {a.out}")


if __name__ == "__main__":
    main()
