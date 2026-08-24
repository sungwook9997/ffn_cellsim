#!/usr/bin/env python3
"""Run the G-Eval judge on a STRATIFIED SUBSET of an already-answered run.

Reuses the stored answers + contexts (no re-answering), so the RAGAS/FActScore
metrics align with the same n questions the programmatic backbone scored, at a
fraction of the cost. Judges `per_class` questions per class × all 4 conditions.

RUN:  conda activate ffn_sim
      python judge_subset.py benchmark_results_scaled.json --per-class 2
"""
from __future__ import annotations

import argparse
import json
import pathlib
from collections import defaultdict

from kb_benchmark import judge

HERE = pathlib.Path(__file__).parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--questions", default=str(HERE / "questions_gen.json"))
    ap.add_argument("--per-class", type=int, default=2)
    ap.add_argument("--judge-model", default="claude-opus-4-8")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    d = json.loads(pathlib.Path(args.results).read_text())
    qmap = {q["id"]: q for q in json.loads(pathlib.Path(args.questions).read_text())}

    # pick the first `per_class` question ids per class (deterministic)
    by_cls = defaultdict(list)
    for q in qmap.values():
        by_cls[q["cls"]].append(q["id"])
    chosen = set()
    for cls, ids in by_cls.items():
        for qid in sorted(ids)[: args.per_class]:
            chosen.add(qid)
    print(f"judging {len(chosen)} questions × 4 conditions = "
          f"{sum(1 for r in d['results'] if r['qid'] in chosen)} answers")

    n = 0
    for r in d["results"]:
        if r["qid"] not in chosen:
            continue
        g = qmap[r["qid"]]
        r["judge"] = judge(g, r.get("context", ""), r["answer"], args.judge_model)
        n += 1
        jc = r["judge"].get("correct")
        print(f"  {r['qid']:5s} {r['cond']:3s} correct={jc} "
              f"faith={r['judge'].get('faithfulness')} "
              f"crecall={r['judge'].get('context_recall')}")

    out = args.out or args.results
    pathlib.Path(out).write_text(json.dumps(d, indent=2))
    print(f"\nattached judge to {n} records -> {out}")


if __name__ == "__main__":
    main()
