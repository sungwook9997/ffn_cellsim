#!/usr/bin/env python3
"""Statistics + threat-investigation for the scaled KB benchmark.

Turns a single benchmark run (n≈100) into defensible numbers and probes three
threats-to-validity from the report:

  T1  statistical power — bootstrap 95% CIs on accuracy/hallucination per
       condition (micro = per-question; macro = per-class mean, removes the
       class-size imbalance from the generated set); McNemar paired test on
       per-question correctness between conditions.
  T5  C4 ceiling = syn quality — for C4, parse the synthesised SQL from the
       stored context and measure how often it returned ZERO rows (a syn miss,
       not a data gap); cross-tabulate with correctness.
  T8  citation metric — classify every cited identifier as resolves-in-KB vs
       unresolvable, per condition (the unresolvable bucket is the hallucination
       signal; "real-but-external" would need CrossRef and is flagged separately).

RUN:  conda activate ffn_sim && python kb_benchmark_stats.py [--in FILE]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
from collections import defaultdict

import numpy as np

HERE = pathlib.Path(__file__).parent
CONDS = ["C1", "C2", "C3", "C4"]
LAB = {"C1": "no-KB", "C2": "RAG", "C3": "+Obsidian", "C4": "+TAG"}

# fixed bootstrap seed for reproducibility (no Date/random-at-runtime concerns here)
RNG = np.random.default_rng(20260602)


def correct(r):
    """unified per-question correctness: programmatic accuracy if present,
    else judge 'correct'. For trap/negative, hallucination==0 counts as correct."""
    if r.get("cls") in ("citation_trap", "negative"):
        return int(not r.get("hallucination", 0))
    if r.get("accuracy") is not None:
        return int(r["accuracy"])
    j = r.get("judge") or {}
    return int(bool(j.get("correct")))


def boot_ci(vals, B=5000):
    vals = np.asarray(vals, float)
    if len(vals) == 0:
        return (float("nan"),) * 3
    means = [RNG.choice(vals, len(vals), replace=True).mean() for _ in range(B)]
    return 100 * vals.mean(), 100 * np.percentile(means, 2.5), 100 * np.percentile(means, 97.5)


def macro_ci(rows, cond, metric, B=5000):
    """bootstrap a per-class-mean (macro) metric: resample within each class."""
    by_cls = defaultdict(list)
    for r in rows:
        if r["cond"] == cond:
            by_cls[r["cls"]].append(metric(r))
    classes = list(by_cls)
    point = 100 * np.mean([np.mean(by_cls[c]) for c in classes])
    draws = []
    for _ in range(B):
        cms = [RNG.choice(by_cls[c], len(by_cls[c]), replace=True).mean() for c in classes]
        draws.append(100 * np.mean(cms))
    return point, np.percentile(draws, 2.5), np.percentile(draws, 97.5)


def mcnemar(rows, ca, cb):
    """paired test on per-question correctness; returns (b, c, chi2_cc, sig)."""
    a = {r["qid"]: correct(r) for r in rows if r["cond"] == ca}
    bb = {r["qid"]: correct(r) for r in rows if r["cond"] == cb}
    qs = set(a) & set(bb)
    b = sum(1 for q in qs if a[q] == 1 and bb[q] == 0)   # ca right, cb wrong
    c = sum(1 for q in qs if a[q] == 0 and bb[q] == 1)   # ca wrong, cb right
    if b + c == 0:
        return b, c, 0.0, False
    chi2 = (abs(b - c) - 1) ** 2 / (b + c)               # continuity-corrected
    return b, c, chi2, chi2 > 3.841                      # p<0.05, df=1


def extract_sql(ctx):
    m = re.search(r"query:\s*(.+?)\n\s*rows:", ctx, re.S)
    return m.group(1).strip() if m else None


def tag_returned_rows(ctx):
    """did the C4 TAG SQL return any data rows? (context shows 'rows:' then data)"""
    m = re.search(r"rows:\n(.*?)(?:\n\n|\Z)", ctx, re.S)
    if not m:
        return None
    body = [ln for ln in m.group(1).splitlines() if ln.strip()]
    # first line is the header; data rows are the rest
    return len(body) > 1


DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\)\]\}>,;\"']+", re.I)
KBID_RE = re.compile(r"KB-\d+(?:\.\d+)+", re.I)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=str(HERE / "benchmark_results_scaled.json"))
    args = ap.parse_args()
    d = json.loads(pathlib.Path(args.inp).read_text())
    rows = d["results"]
    n = d.get("n_questions", len({r["qid"] for r in rows}))
    print(f"=== scaled stats: model={d['model']}  n={n} questions × {len(CONDS)} conditions ===\n")

    # ---- T1: accuracy with CIs (micro + macro) + McNemar -------------------
    print("ACCURACY (correct %), bootstrap 95% CI:")
    print(f"{'cond':12s} {'micro [95% CI]':>26s}   {'macro (per-class) [95% CI]':>30s}")
    for cnd in CONDS:
        cs = [correct(r) for r in rows if r["cond"] == cnd]
        mi, lo, hi = boot_ci(cs)
        ma, mlo, mhi = macro_ci(rows, cnd, correct)
        print(f"{LAB[cnd]:12s} {mi:5.1f} [{lo:4.1f}, {hi:4.1f}]      "
              f"{ma:5.1f} [{mlo:4.1f}, {mhi:4.1f}]")
    print("\nMcNemar paired test (per-question correctness):")
    for ca, cb in [("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C1", "C4")]:
        b, c, chi2, sig = mcnemar(rows, ca, cb)
        print(f"  {LAB[ca]:10s} vs {LAB[cb]:10s}: {cb}>{ca} on {c}q, {ca}>{cb} on {b}q, "
              f"χ²={chi2:5.1f}  {'SIGNIFICANT (p<.05)' if sig else 'n.s.'}")

    # ---- per-class accuracy table ------------------------------------------
    classes = sorted({r["cls"] for r in rows})
    print("\nACCURACY by question class × condition:")
    print(f"{'class':16s} " + " ".join(f"{LAB[c]:>10s}" for c in CONDS) + "   n")
    for cls in classes:
        line = f"{cls:16s} "
        nq = len({r["qid"] for r in rows if r["cls"] == cls})
        for cnd in CONDS:
            cs = [correct(r) for r in rows if r["cond"] == cnd and r["cls"] == cls]
            line += f"{100*np.mean(cs):9.0f}%" if cs else f"{'-':>10s}"
        print(line + f"   {nq}")

    # ---- T5: C4 syn-quality ceiling ----------------------------------------
    print("\nT5  C4 ceiling = syn quality:")
    c4 = [r for r in rows if r["cond"] == "C4"]
    got_rows = [r for r in c4 if tag_returned_rows(r.get("context", "")) is True]
    no_rows = [r for r in c4 if tag_returned_rows(r.get("context", "")) is False]
    structured = [r for r in c4 if r["cls"] in ("aggregation", "relational", "lookup")]
    s_norows = [r for r in structured if tag_returned_rows(r.get("context", "")) is False]
    print(f"  C4 questions where TAG SQL returned ZERO rows: {len(no_rows)}/{len(c4)} "
          f"({100*len(no_rows)/max(len(c4),1):.0f}%)")
    print(f"  of structured (agg/rel/lookup) Qs, zero-row (syn miss candidates): "
          f"{len(s_norows)}/{len(structured)} ({100*len(s_norows)/max(len(structured),1):.0f}%)")
    acc_rows = np.mean([correct(r) for r in got_rows]) if got_rows else float("nan")
    acc_norows = np.mean([correct(r) for r in no_rows]) if no_rows else float("nan")
    print(f"  C4 accuracy WHEN SQL returned rows: {100*acc_rows:.0f}%  |  "
          f"when it returned none: {100*acc_norows:.0f}%")
    print("  -> if 'returned rows' accuracy ≫ 'no rows', the C4 ceiling is syn (NL→SQL), not data.")

    # ---- T8: citation classification ---------------------------------------
    print("\nT8  citation grounding (cited identifiers, per condition):")
    print(f"{'cond':12s} {'#cited':>7s} {'resolves-in-KB':>15s} {'unresolvable':>13s}")
    for cnd in CONDS:
        cited = res = unres = 0
        for r in rows:
            if r["cond"] != cnd:
                continue
            cited += r.get("n_cited", 0)
            res += r.get("n_resolved", 0)
            unres += r.get("n_unresolved", 0)
        print(f"{LAB[cnd]:12s} {cited:7d} {res:15d} {unres:13d}")
    print("  (unresolvable = hallucination signal; a real-but-external paper not in the\n"
          "   KB would also land here — disambiguating that needs CrossRef, flagged as future work.)")


if __name__ == "__main__":
    main()
