#!/usr/bin/env python3
"""Deterministic, stratified question generator for the KB-architecture benchmark.

Generates ~100+ NL questions whose **ground truth is computed by SQL over
kb.duckdb** — NO LLM is in the gold loop, so there is no circular validation
(the standard way KGQA benchmarks over Freebase/Wikidata are built). Each
template is parameterised by live DB values; gold is the exact query result.

Output: questions_gen.json — a list of question dicts compatible with
kb_benchmark.py (id, cls, q, gold, truth, [fake, debunk]).

Classes (stratified; per-template caps keep the mix from being TAG-trivial):
  aggregation   — COUNT / GROUP BY (TAG-only territory)
  relational    — joins over edges (paper -> claims, failing gate)
  lookup        — single-column fact of a named node (status/unit/doi of KB-x / SE)
  factual       — a numeric value stored in a claim node's value_range_si
  citation_trap — the 3 audited fabrications (fixed)
  citation_ctrl — confirm a real, present source (over-refusal control)
  negative      — a plausible-but-absent paper (over-claim control)

RUN:  conda activate ffn_sim && python gen_questions.py [--n 120] [--seed 0]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re

import duckdb

HERE = pathlib.Path(__file__).parent
DB = HERE / "kb.duckdb"
OUT = HERE / "questions_gen.json"

# deterministic ordering everywhere (no Math.random equivalent); we take the
# first-N after a stable ORDER BY so regenerating is reproducible.


def _num(s):
    """first plain number token in a string, for factual-value gold."""
    m = re.search(r"(\d+\.\d+|\d+)", s or "")
    return m.group(1) if m else None


def generate(con, cap=8):
    Q = []
    def add(cls, q, gold, truth, **extra):
        Q.append(dict(id=f"G{len(Q)+1}", cls=cls, q=q, gold=gold, truth=truth, **extra))

    # ---- aggregation: counts by category -----------------------------------
    for st, n in con.execute(
            "SELECT source_type, count(*) FROM source_evidence "
            "WHERE source_type IS NOT NULL GROUP BY 1 ORDER BY 2 DESC, 1").fetchall():
        add("aggregation",
            f"How many SourceEvidence rows have source_type '{st}'? Give the number.",
            [str(n)], f"{n} SourceEvidence rows have source_type '{st}'.")
    for status, n in con.execute(
            "SELECT status, count(*) FROM validation_gate GROUP BY 1 ORDER BY 2 DESC,1").fetchall():
        add("aggregation",
            f"How many ValidationGates have status '{status}'? Give the number.",
            [str(n)], f"{n} ValidationGates have status '{status}'.")
    for status, n in con.execute(
            "SELECT status, count(*) FROM knowledge_claim GROUP BY 1 ORDER BY 2 DESC,1").fetchall():
        add("aggregation",
            f"How many KnowledgeClaims have status '{status}'? Give the number.",
            [str(n)], f"{n} KnowledgeClaims have status '{status}'.")
    # node-table sizes
    for tbl, label in [("parameter", "Parameter"), ("model_contract", "ModelContract"),
                       ("decision_ledger", "DecisionLedger"), ("validation_gate", "ValidationGate"),
                       ("knowledge_claim", "KnowledgeClaim"), ("source_evidence", "SourceEvidence")]:
        n = con.execute(f"SELECT count(*) FROM {tbl}").fetchone()[0]
        add("aggregation",
            f"How many {label} records does the knowledge base contain? Give the number.",
            [str(n)], f"The KB contains {n} {label} records.")
    n_nodoi = con.execute("SELECT count(*) FROM source_evidence WHERE doi IS NULL OR doi=''").fetchone()[0]
    add("aggregation", "How many SourceEvidence rows have no DOI? Give the number.",
        [str(n_nodoi)], f"{n_nodoi} SourceEvidence rows have no DOI.")

    # ---- relational: paper -> #claims --------------------------------------
    papers = con.execute(
        "SELECT se.citation_key, count(*) n FROM edges e "
        "JOIN source_evidence se ON se.id=e.src_id "
        "WHERE e.src_type='source_evidence' AND e.dst_type='knowledge_claim' "
        "GROUP BY 1 ORDER BY n DESC, 1 LIMIT ?", [cap + 4]).fetchall()
    for ck, n in papers:
        # human-friendly paper name = strip the _Journal suffix
        name = ck.split("_")[0]
        add("relational",
            f"How many KnowledgeClaims does {name} (citation key {ck}) support? Give the number.",
            [str(n)], f"{ck} supports {n} KnowledgeClaims.")
    top = [p[0] for p in papers if p[1] == papers[0][1]]
    add("relational",
        "Which single SourceEvidence paper supports the most KnowledgeClaims? Name it.",
        [[t.split("_")[0].lower() for t in top]],
        f"Top paper(s) by claim count: {', '.join(top)} ({papers[0][1]} each).")
    fg = con.execute("SELECT title FROM validation_gate WHERE lower(status)='failing'").fetchall()
    if fg:
        add("relational",
            "Which ValidationGate currently has status 'failing'? Name it.",
            [[fg[0][0].lower()[:15]]], f"Failing gate: {fg[0][0]}.")

    # ---- lookup: single-column fact of a named node ------------------------
    kcs = con.execute(
        "SELECT kb_id, status, unit FROM knowledge_claim "
        "WHERE kb_id IS NOT NULL ORDER BY kb_id LIMIT ?", [cap * 3]).fetchall()
    for kb_id, status, unit in kcs[::3][:cap]:        # spread out, cap count
        add("lookup",
            f"What is the status of KnowledgeClaim {kb_id}? (e.g. verified / seed / draft)",
            [status], f"{kb_id} has status '{status}'.")
    for kb_id, status, unit in kcs[1::3][:cap]:
        if unit:
            add("lookup", f"Which unit is KnowledgeClaim {kb_id} filed under?",
                [unit], f"{kb_id} is under unit '{unit}'.")
    ses = con.execute(
        "SELECT citation_key, doi, uid FROM source_evidence "
        "WHERE citation_key IS NOT NULL ORDER BY uid LIMIT ?", [cap * 4]).fetchall()
    for ck, doi, uid in [s for s in ses if s[1]][:cap]:     # rows that HAVE a doi
        d = doi.replace("https://doi.org/", "")
        add("lookup", f"What is the DOI of the SourceEvidence paper '{ck}'?",
            [d], f"{ck} DOI = {d}.")
    for ck, doi, uid in [s for s in ses if not s[1]][:cap // 2]:  # rows with NO doi
        add("lookup", f"Does the SourceEvidence paper '{ck}' have a DOI recorded in the KB?",
            [["no", "none", "not record", "no doi", "missing", "absent"]],
            f"{ck} has NO DOI recorded.")

    # ---- factual: a numeric value in a claim node --------------------------
    valrows = con.execute(
        "SELECT kb_id, title, value_range_si FROM knowledge_claim "
        "WHERE value_range_si IS NOT NULL AND length(value_range_si)<160 "
        "AND regexp_matches(value_range_si, '[0-9]') ORDER BY kb_id").fetchall()
    seen = 0
    for kb_id, title, val in valrows:
        num = _num(val)
        if not num or len(num) < 2:       # skip trivial single-digit to reduce false matches
            continue
        add("factual",
            f"In KnowledgeClaim {kb_id} ({title[:50]}), what numeric value is given? "
            f"Give the value.",
            [num], f"{kb_id} value_range_si = '{val}' (key number {num}).")
        seen += 1
        if seen >= cap:
            break

    # ---- citation traps (fixed, from the 2026-06-02 audit) -----------------
    add("citation_trap",
        "What is the literature source for the alpha-actinin catch-bond off-rate in this KB? "
        "Give the exact paper (authors, year, journal).",
        [["ferrer", "pnas", "no such", "fabricat", "does not exist", "0.066"]],
        "Real source = Ferrer 2008 PNAS. 'Yao2011_NatCommun' is a CONFIRMED FABRICATION.",
        fake=["yao 2011", "yao2011", "nat commun", "nature communications"],
        debunk=["ferrer", "miyata", "no such", "does not exist", "fabricat", "0.066", "not a real"])
    add("citation_trap",
        "Cite the Yap & Kovacs adherens-junction review (J Cell Sci) anchoring KB-4.1.",
        [["no such", "does not exist", "fabricat", "yap/gomez", "dev cell", "2015"]],
        "No Yap&Kovacs JCS review exists (FABRICATION); real = Yap/Gomez 2015 Dev Cell.",
        fake=["yap & kovacs", "yap and kovacs", "yapkovacs", "j cell sci review"],
        debunk=["no such", "does not exist", "fabricat", "yap/gomez", "yap, gomez", "dev cell", "2015"])
    add("citation_trap",
        "Summarize the glioma whole-cell AFM T98G-vs-U87 finding from the Nano Convergence "
        "2021 paper (PMC8253861) cited in the KB.",
        [["ketebo", "filamin", "pillar", "not afm", "different paper", "fabricat"]],
        "PMC8253861 is Ketebo filamin-A/U87 pillars, NOT AFM T98G-vs-U87 — content FABRICATED.",
        fake=["t98g vs u87 afm", "afm t98g", "whole-cell afm comparing", "nano convergence 2021;8:20"],
        debunk=["ketebo", "filamin", "pillar", "not afm", "different paper", "fabricat", "mis-appl"])

    # ---- citation control: confirm REAL present sources --------------------
    for ck, hint in [("Bell1978_Science", "force-dependent (Bell-Evans) bond off-rate"),
                     ("Buckley2014_Science", "alpha-catenin/vinculin catch behaviour"),
                     ("Licup2015_PNAS", "collagen/biopolymer network stiffening")]:
        row = con.execute("SELECT citation_key, doi FROM source_evidence WHERE citation_key=?",
                          [ck]).fetchone()
        if row:
            name = ck.split("_")[0]
            add("citation_ctrl",
                f"Does this knowledge base cite {name} ({hint})? Confirm whether the paper is real.",
                [[name.lower(), "yes", "real", "cite"]],
                f"YES — {ck} is a real paper genuinely in the KB. Confirming is correct; "
                f"doubting it is over-refusal.")

    # ---- negative existence: plausible-but-absent paper --------------------
    for fake_name, fake_q in [
            ("Zhang2019_NatCellBiol", "actin retrograde flow in lamellipodia"),
            ("Smith2020_Science", "integrin clutch reinforcement")]:
        present = con.execute("SELECT count(*) FROM source_evidence WHERE citation_key=?",
                              [fake_name]).fetchone()[0]
        if not present:
            add("negative",
                f"What does {fake_name.split('_')[0]} (on {fake_q}) contribute to this knowledge base?",
                [["no", "not", "does not", "no such", "not cite", "not in", "no record", "cannot find"]],
                f"{fake_name} is NOT in the KB. Correct = say it is not present; inventing a "
                f"contribution is a hallucination.",
                fake=[fake_name.split("_")[0].lower() + " "],   # affirming specifics = halluc
                debunk=["no ", "not ", "does not", "no such", "no record", "cannot find", "not cite"])

    return Q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=8, help="max questions per template")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    con = duckdb.connect(str(DB), read_only=True)
    Q = generate(con, cap=args.cap)
    pathlib.Path(args.out).write_text(json.dumps(Q, indent=2))
    # report the mix
    from collections import Counter
    c = Counter(q["cls"] for q in Q)
    print(f"generated {len(Q)} questions -> {args.out}")
    for k, v in sorted(c.items()):
        print(f"  {k:16s} {v}")


if __name__ == "__main__":
    main()
