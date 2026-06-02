#!/usr/bin/env python3
"""TAG-bench-style eval for the ffn_cellsim TAG engine.

A small gold set of NL questions spanning the question classes TAG must handle —
aggregation, relational join, multi-hop, audit-table lookup, parameter lookup,
and PDF-content synthesis — each with substrings the answer MUST contain. Runs
`tag_query.py` end-to-end (syn -> exec -> gen) and checks the generated answer.

This is a regression guard: re-run after engine/schema/KB changes.

RUN:  conda activate ffn_sim && python tag_eval.py
"""
from __future__ import annotations

import re
import subprocess
import sys

HERE = __import__("pathlib").Path(__file__).parent

# (question, [substrings the gen answer must contain — case-insensitive], note)
GOLD = [
    ("How many SourceEvidence rows have no DOI? Give the number.",
     [r"\d"], "aggregation / count"),
    ("Which single SourceEvidence paper supports the most KnowledgeClaims?",
     ["broedersz", "licup"], "relational join + group-by (any top paper)"),
    ("Using the source_audit table, which citation keys have web_verdict = HALLUCINATION?",
     ["yao2011", "yapkovacs", "nanoconvergence"], "audit-table lookup (all 3)"),
    ("What alpha-actinin k_off value does KnowledgeClaim KB-3.19 give?",
     ["0.066"], "parameter lookup (post-fix value)"),
    ("List ValidationGates whose status is failing.",
     ["vg-", "fail"], "status filter"),
    ("According to the reference PDFs, what sets a spheroid or tissue's surface tension?",
     ["cortic", "tension"], "PDF-content synthesis (BM25)"),
]


def run(q):
    r = subprocess.run([sys.executable, str(HERE / "tag_query.py"), q],
                       capture_output=True, text=True, timeout=300)
    out = r.stdout
    # the answer is everything after the gen header
    ans = out.split("gen (answer)")[-1].lower()
    return out, ans


def check(ans, needles):
    """pass if ANY needle group matches (substring or regex)."""
    hits = []
    for n in needles:
        ok = bool(re.search(n, ans)) if any(c in n for c in ".\\[]+*") else (n.lower() in ans)
        hits.append((n, ok))
    return any(ok for _, ok in hits), hits


def main():
    only = None
    if "--n" in sys.argv:
        only = int(sys.argv[sys.argv.index("--n") + 1])
    npass = 0
    gold = GOLD[:only] if only else GOLD
    for i, (q, needles, note) in enumerate(gold, 1):
        try:
            _, ans = run(q)
        except Exception as e:
            print(f"[{i}] ERROR  {note}: {e}")
            continue
        ok, hits = check(ans, needles)
        npass += ok
        mark = "\033[32mPASS\033[0m" if ok else "\033[31mFAIL\033[0m"
        got = ", ".join(f"{n}{'✓' if h else '✗'}" for n, h in hits)
        print(f"[{i}] {mark}  {note}\n      Q: {q[:70]}\n      need: {got}")
    print(f"\n=== {npass}/{len(gold)} gold questions passed ===")


if __name__ == "__main__":
    main()
