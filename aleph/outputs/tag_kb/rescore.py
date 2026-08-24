#!/usr/bin/env python3
"""Re-apply programmatic scoring to a stored benchmark run — NO LLM calls.

Use after fixing a scorer bug (e.g. the DOI-prefix normalisation in
score_groundedness): recomputes n_cited / n_resolved / n_unresolved and the
hallucination flag from the already-stored answers, in place.

RUN:  python rescore.py benchmark_results_scaled.json [--questions questions_gen.json]
"""
from __future__ import annotations

import argparse
import json
import pathlib

import duckdb
from kb_benchmark import score_groundedness, score_trap, DB_PATH

HERE = pathlib.Path(__file__).parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--questions", default=str(HERE / "questions_gen.json"))
    args = ap.parse_args()
    d = json.loads(pathlib.Path(args.results).read_text())
    qmap = {}
    qpath = pathlib.Path(args.questions)
    if qpath.exists():
        for q in json.loads(qpath.read_text()):
            qmap[q["id"]] = q

    con = duckdb.connect(str(DB_PATH), read_only=True)
    chg = 0
    for r in d["results"]:
        ans = r["answer"]
        nc, nr, nu = score_groundedness(con, ans)
        old = (r.get("n_cited"), r.get("n_resolved"), r.get("n_unresolved"),
               r.get("hallucination"))
        r["n_cited"], r["n_resolved"], r["n_unresolved"] = nc, nr, nu
        q = qmap.get(r["qid"], {})
        if q.get("fake"):
            hall, caught = score_trap(q, ans)
            r["hallucination"], r["caught_fabrication"] = hall, caught
        else:
            r["hallucination"] = int(nu > 0)
        if (nc, nr, nu, r["hallucination"]) != old:
            chg += 1
    d["rescored"] = True
    pathlib.Path(args.results).write_text(json.dumps(d, indent=2))
    print(f"rescored {len(d['results'])} records ({chg} changed) -> {args.results}")


if __name__ == "__main__":
    main()
