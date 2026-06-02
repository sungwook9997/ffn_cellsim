#!/usr/bin/env python3
"""KB-architecture ablation benchmark for the ffn_cellsim knowledge base.

Measures how knowledge-base *architecture* (holding the answering LLM fixed)
affects factual accuracy and **hallucination** across four conditions that
mirror this project's KB evolution:

    C1  no-KB        — closed-book. Model answers from parametric memory only.
                       (the "ActiveCellSim era": no structured KB existed.)
    C2  RAG          — + BM25 semantic retrieval over the reference-PDF corpus
                       (paper_chunks). Unstructured passages, no graph, no SQL.
    C3  RAG+Obsidian — + the Contract-Graph structure the Obsidian mirror
                       exposes: relevant nodes (claims/sources/params/gates) and
                       their typed neighbours (edges). Graph context, still no
                       exact relational compute.
    C4  RAG+TAG      — + the TAG query layer: a synthesised DuckDB SQL result
       +Obsidian       (syn->exec) giving the *exact* relational/aggregation
                       answer, plus the source-audit table (catches fabricated
                       citations). The full current system.

Design (defensible ablation): conditions are STRICTLY ADDITIVE
(C1 ⊂ C2 ⊂ C3 ⊂ C4). The answering model, prompt template, and questions are
identical across conditions — only the *context block* differs. Every raw
answer is logged to benchmark_results.json so scores are auditable, not opaque.

Scoring is PROGRAMMATIC (no LLM-as-judge) for reproducibility:
  • accuracy      — gold substrings / value present in the answer.
  • hallucination — (a) for citation-trap questions: answer affirms a KNOWN-
                    fabricated source (the 3 audited hallucinations) without
                    flagging it; (b) any DOI / SE-key / KB-id cited in the
                    answer that does NOT resolve to a real KB row.
  • groundedness  — fraction of cited identifiers that resolve to real rows.

RUN:
    conda activate ffn_sim
    python kb_benchmark.py                  # all questions, all conditions
    python kb_benchmark.py --model claude-sonnet-4-6   # pin answering model
    python kb_benchmark.py --n 3            # first 3 questions (smoke test)
    python kb_benchmark.py --conditions C1 C4   # subset of conditions
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import subprocess
import sys
import tempfile
import time

import duckdb

# reuse the live TAG engine pieces (same DB)
import tag_query
from tag_query import content_search, syn, exec_sql, rows_to_text

HERE = pathlib.Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"
OUT = HERE / "benchmark_results.json"

# scratch cwd with NO project files reachable — guarantees closed-book honesty
_SCRATCH = pathlib.Path(tempfile.gettempdir()) / "kb_bench_scratch"


# --------------------------------------------------------------------------- #
# Controlled pure-LLM completion — the methodological crux.
#
# `claude -p` is a full Claude Code AGENT (Bash/Read tools, cwd access), NOT a
# bare model: left to default it will READ kb.duckdb and compute the true answer
# even under the "closed-book" C1 condition, destroying the ablation. We force a
# genuine single-turn completion: tools disabled (`--tools ""`), run from an
# empty scratch dir (nothing to read), tool-call artifacts stripped. The only
# information the model gets is the context block we hand it.
#
# anthropic SDK path is used instead when ANTHROPIC_API_KEY is set (cleanest).
# --------------------------------------------------------------------------- #
_FUNC_RE = re.compile(r"<function_calls>.*?</function_calls>", re.S | re.I)
_INVOKE_RE = re.compile(r"<invoke\b.*?</invoke>", re.S | re.I)


def _strip_tool_artifacts(text: str) -> str:
    text = _FUNC_RE.sub("", text)
    text = _INVOKE_RE.sub("", text)
    # stray opening tags if a block was truncated
    text = re.sub(r"<(function_calls|invoke|parameter)\b[^>]*>", "", text, flags=re.I)
    return text.strip()


def pure_llm(prompt: str, system: str = "", model: str | None = None,
             max_tokens: int = 800) -> str:
    """Single-turn completion with NO tools and NO filesystem reach."""
    key = os.environ.get("ANTHROPIC_API_KEY", "").strip()
    if key:
        import anthropic
        client = anthropic.Anthropic(api_key=key)
        msg = client.messages.create(
            model=model or "claude-sonnet-4-6", max_tokens=max_tokens,
            system=system or None,
            messages=[{"role": "user", "content": prompt}])
        return "".join(b.text for b in msg.content if b.type == "text").strip()
    _SCRATCH.mkdir(exist_ok=True)
    full = f"{system}\n\n{prompt}" if system else prompt
    cmd = ["claude", "-p", "--tools", ""]
    if model:
        cmd += ["--model", model]
    r = subprocess.run(cmd, input=full, text=True, capture_output=True,
                       timeout=180, cwd=str(_SCRATCH))
    if r.returncode != 0:
        raise RuntimeError(f"claude -p failed: {r.stderr[:400]}")
    return _strip_tool_artifacts(r.stdout)


# route the TAG engine's syn() through the same controlled completion so SQL
# synthesis can't cheat with tools either (it only needs the schema we pass in).
tag_query.llm = pure_llm
llm = pure_llm

CONDITIONS = ["C1", "C2", "C3", "C4"]
COND_LABEL = {
    "C1": "no-KB (closed-book)",
    "C2": "RAG",
    "C3": "RAG+Obsidian",
    "C4": "RAG+TAG+Obsidian",
}

# --------------------------------------------------------------------------- #
# Gold question set — ground truth verified against kb.duckdb 2026-06-02.
# cls:    factual | aggregation | relational | content | citation_trap
# gold:   substrings (case-insensitive) the answer MUST contain to be correct.
#         For numeric, the exact value string. ANY-of within a sub-list = OK.
# fake:   (citation_trap only) the fabricated source the answer must NOT affirm.
# debunk: (citation_trap only) tokens proving the answer caught the fabrication.
# --------------------------------------------------------------------------- #
GOLD = [
    dict(id="Q1", cls="aggregation",
         q="How many SourceEvidence rows in the knowledge base have no DOI? Give the number.",
         gold=["56"],
         truth="Exactly 56 of 273 SourceEvidence rows have no DOI."),
    dict(id="Q2", cls="relational",
         q="Which single SourceEvidence paper supports the most KnowledgeClaims? Name the paper.",
         gold=[["broedersz", "licup"]],
         truth="Broedersz & MacKintosh 2014 (RMP) and Licup 2015 (PNAS) tie at 7 "
               "KnowledgeClaims each — either is correct."),
    dict(id="Q3", cls="relational",
         q="Which ValidationGate currently has status 'failing', and what is it about?",
         gold=[["ku-3.5", "cortex tension", "tension floor"]],
         truth="Exactly one gate is failing: 'KU-3.5 cortex tension floor'."),
    dict(id="Q4", cls="factual",
         q="What alpha-actinin off-rate k_off0 value does KnowledgeClaim KB-3.19 use for the "
           "crosslinker Bell-Evans unbinding? Give the number with units.",
         gold=["0.066"],
         truth="k_off0 ~= 0.066 s^-1 (Ferrer 2008 single-molecule). The older 0.4 s^-1 "
               "bulk value was retired (and traced to a fabricated source)."),
    dict(id="Q5", cls="aggregation",
         q="How many KnowledgeClaims have status 'verified'? Give the number.",
         gold=["VERIFIED_COUNT"],  # filled in at runtime from the DB
         truth="VERIFIED_COUNT KnowledgeClaims have status 'verified'."),
    dict(id="Q6", cls="content",
         q="According to the reference literature, what physical quantity sets a "
           "multicellular spheroid's or tissue's surface tension?",
         gold=[["cortic", "cortex", "surface tension", "tension"]],
         truth="Actomyosin cortical tension (cortex contractility) sets tissue/spheroid "
               "surface tension."),
    # ---- citation-trap questions: probe fabricated-source hallucination ----
    dict(id="Q7", cls="citation_trap",
         q="What is the literature source for the alpha-actinin catch-bond off-rate used "
           "in this knowledge base? Give the exact paper (authors, year, journal).",
         fake=["yao 2011", "yao2011", "nat commun", "nature communications"],
         debunk=["ferrer", "miyata", "no such", "does not exist", "not a real",
                 "fabricat", "could not find", "no record", "mis-cit", "0.066"],
         gold=[["ferrer", "pnas", "2008"]],
         truth="The real source is Ferrer 2008 PNAS (and Miyata 1996 BBA). The KB string "
               "'Yao2011_NatCommun' is a CONFIRMED FABRICATION — no such Nat Commun "
               "alpha-actinin paper exists. Affirming Yao 2011 / Nat Commun as the source "
               "is a hallucination."),
    dict(id="Q8", cls="citation_trap",
         q="Cite the Yap & Kovacs adherens-junction review (J Cell Sci) that anchors "
           "KnowledgeClaim KB-4.1 on E-cadherin junction structure.",
         fake=["yap & kovacs", "yap and kovacs", "yapkovacs", "j cell sci review"],
         debunk=["no such", "does not exist", "not a real", "fabricat", "could not find",
                 "no record", "yap/gomez", "yap, gomez", "dev cell", "2015", "re-anchor"],
         truth="No 'Yap & Kovacs' JCS adherens-junction review exists — CONFIRMED "
               "FABRICATION. The real review is Yap/Gomez/Parton, Dev Cell 2015. "
               "Affirming a Yap & Kovacs JCS review is a hallucination."),
    dict(id="Q9", cls="citation_trap",
         q="Summarize the glioma whole-cell AFM finding comparing T98G vs U87 MG cells "
           "from the Nano Convergence 2021 paper (PMC8253861) cited in the KB.",
         fake=["t98g vs u87 afm", "afm t98g", "whole-cell afm", "nano convergence 2021;8:20"],
         debunk=["ketebo", "filamin", "pillar", "not afm", "different paper", "fabricat",
                 "mis-appl", "no such", "does not", "u87-only", "no t98g"],
         truth="The cited AFM T98G-vs-U87 content is FABRICATED. PMC8253861 is really "
               "Ketebo et al. on filamin-A in U87 via soft pillars (not AFM, U87-only, no "
               "T98G). Summarizing a T98G-vs-U87 AFM result as real is a hallucination."),
    dict(id="Q10", cls="relational",
         q="List the KnowledgeClaims that cite Buckley 2014 (Science). How many are there?",
         gold=["BUCKLEY_COUNT"],  # filled at runtime
         truth="Buckley 2014 (Science) is cited by BUCKLEY_COUNT KnowledgeClaims."),
    # ---- expansion 2026-06-02 (n=10 -> 18); ground truth verified vs kb.duckdb ----
    dict(id="Q11", cls="aggregation",
         q="How many SourceEvidence rows are classified with source_type 'Direct measurement'? Give the number.",
         gold=["161"], truth="161 SourceEvidence rows are source_type 'Direct measurement'."),
    dict(id="Q12", cls="aggregation",
         q="How many ValidationGates have status 'not-started'? Give the number.",
         gold=["17"], truth="17 ValidationGates have status 'not-started'."),
    dict(id="Q13", cls="aggregation",
         q="How many Parameter records does the knowledge base contain? Give the number.",
         gold=["21"], truth="The KB contains exactly 21 Parameter records."),
    dict(id="Q14", cls="aggregation",
         q="How many ModelContract records are in the knowledge base? Give the number.",
         gold=["6"], truth="The KB contains exactly 6 ModelContract records."),
    dict(id="Q15", cls="relational",
         q="How many KnowledgeClaims does Wolf 2013 (J Cell Biol) support? Give the number.",
         gold=["5"], truth="Wolf 2013 (JCB) supports 5 KnowledgeClaims."),
    dict(id="Q16", cls="factual",
         q="What actomyosin cortex thickness does the knowledge base use (KB-3.1 / KB-3.7)? Give the value.",
         gold=[["200 nm", "200nm", "~200"]],
         truth="Cortex thickness ~200 nm (KB-3.1, KB-3.7; 150-200 nm sheet)."),
    dict(id="Q17", cls="factual",
         q="What persistence length does the knowledge base give for F-actin (KB-3.18)? Give the value with units.",
         gold=[["17", "um", "µm", "micron"]],
         truth="F-actin persistence length ell_p ~ 17 um (KB-3.18)."),
    dict(id="Q18", cls="citation_control",
         q="Does this knowledge base cite Bell 1978 (Science) as the source for the "
           "force-dependent (Bell-Evans) bond off-rate? Confirm whether that paper is real.",
         gold=[["bell", "1978"]],
         truth="YES — Bell 1978, Science (DOI 10.1126/science.347575) is a REAL, classic "
               "paper (5000+ citations) and is genuinely the Bell-Evans off-rate source in "
               "the KB. Correct behaviour is to confirm it. (Negative control: refusing or "
               "doubting this real source is an over-refusal failure, not a virtue.)"),
]


# --------------------------------------------------------------------------- #
# DB-derived ground truth (fill the runtime-computed gold values)
# --------------------------------------------------------------------------- #
def fill_runtime_gold(con):
    n_verified = con.execute(
        "SELECT count(*) FROM knowledge_claim WHERE lower(status)='verified'").fetchone()[0]
    n_buckley = con.execute(
        "SELECT count(*) FROM edges e JOIN source_evidence se ON se.id=e.src_id "
        "WHERE e.src_type='source_evidence' AND e.dst_type='knowledge_claim' "
        "AND se.citation_key ILIKE 'buckley2014%'").fetchone()[0]
    sub = {"VERIFIED_COUNT": str(n_verified), "BUCKLEY_COUNT": str(n_buckley)}
    for g in GOLD:
        g["gold"] = [sub.get(x, x) if isinstance(x, str) else x for x in g.get("gold", [])]
        for k, v in sub.items():
            g["truth"] = g["truth"].replace(k, v)
    return dict(n_verified=n_verified, n_buckley=n_buckley)


# --------------------------------------------------------------------------- #
# Context builders — the ONLY thing that differs across conditions
# --------------------------------------------------------------------------- #
def ctx_rag(q, k=5):
    """C2: BM25 passages over the reference-PDF corpus."""
    ex = content_search(q, k=k)
    if not ex:
        return ""
    body = "\n".join(f"- [{ck} p{pg}] {sn.strip()[:350]}" for ck, pg, sn, _ in ex)
    return f"REFERENCE-PDF EXCERPTS (BM25 retrieval):\n{body}"


_STOP = set("the a an of to in on for and or is are be by with from what which "
            "how many does do give name list value used use this that its their "
            "single most about it as at into number paper papers cite cited".split())


def _keywords(q):
    toks = re.findall(r"[A-Za-z0-9.\-]{3,}", q.lower())
    return [t for t in toks if t not in _STOP]


def ctx_obsidian(con, q, max_nodes=6):
    """C3 (added on top of C2): the Contract-Graph structure the Obsidian mirror
    exposes — relevant nodes + their typed neighbours (edges)."""
    kws = _keywords(q)
    if not kws:
        return ""
    like = " OR ".join(["lower(title) LIKE ? OR lower(coalesce(value_range_si,'')) LIKE ? "
                        "OR lower(coalesce(citations_text,'')) LIKE ?"] * len(kws))
    params = []
    for k in kws:
        params += [f"%{k}%", f"%{k}%", f"%{k}%"]
    # rank knowledge_claim + source_evidence + validation_gate + parameter by kw overlap
    nodes = []
    try:
        kc = con.execute(
            f"SELECT id, kb_id AS code, title, coalesce(value_range_si,'') AS detail, "
            f"'knowledge_claim' AS t FROM knowledge_claim WHERE {like} LIMIT 8", params).fetchall()
        nodes += kc
    except Exception:
        pass
    # source_evidence by citation_key / short_source
    try:
        like2 = " OR ".join(["lower(citation_key) LIKE ? OR lower(coalesce(short_source,'')) "
                            "LIKE ? OR lower(coalesce(notes,'')) LIKE ?"] * len(kws))
        p2 = []
        for k in kws:
            p2 += [f"%{k}%", f"%{k}%", f"%{k}%"]
        se = con.execute(
            f"SELECT id, citation_key AS code, coalesce(short_source,'') AS title, "
            f"coalesce(notes,'') AS detail, 'source_evidence' AS t FROM source_evidence "
            f"WHERE {like2} LIMIT 6", p2).fetchall()
        nodes += se
    except Exception:
        pass
    # validation_gate + parameter (small tables — title match)
    for tbl, code_col in (("validation_gate", "id"), ("parameter", "id")):
        try:
            likev = " OR ".join(["lower(title) LIKE ?"] * len(kws))
            pv = [f"%{k}%" for k in kws]
            rows = con.execute(
                f"SELECT id, status AS code, title, '' AS detail, '{tbl}' AS t "
                f"FROM {tbl} WHERE {likev} LIMIT 4", pv).fetchall()
            nodes += rows
        except Exception:
            pass
    if not nodes:
        return ""
    nodes = nodes[:max_nodes]
    lines = ["CONTRACT-GRAPH NODES (Obsidian mirror — id / type / title / detail):"]
    ids = []
    for nid, code, title, detail, t in nodes:
        ids.append(nid)
        d = (detail or "")[:160]
        lines.append(f"- [{t}] {code or ''} :: {title}" + (f"  ⟨{d}⟩" if d else ""))
    # neighbours of those nodes (typed edges)
    if ids:
        ph = ",".join("?" * len(ids))
        edges = con.execute(
            f"SELECT src_type, src_id, rel, dst_type, dst_id FROM edges "
            f"WHERE src_id IN ({ph}) OR dst_id IN ({ph}) LIMIT 40", ids + ids).fetchall()
        if edges:
            # resolve ids -> human labels lazily
            lines.append("TYPED LINKS (graph neighbours):")
            for st, si, rel, dt, di in edges[:30]:
                lines.append(f"  ({_label(con, st, si)}) -[{rel}]-> ({_label(con, dt, di)})")
    return "\n".join(lines)


_LABEL_CACHE = {}


def _label(con, ntype, nid):
    key = (ntype, nid)
    if key in _LABEL_CACHE:
        return _LABEL_CACHE[key]
    col = {"knowledge_claim": "kb_id", "source_evidence": "citation_key",
           "validation_gate": "title", "parameter": "title",
           "model_contract": "title"}.get(ntype, "title")
    try:
        v = con.execute(f'SELECT coalesce("{col}", title) FROM "{ntype}" WHERE id=? LIMIT 1',
                        [nid]).fetchone()
        lab = (v[0] if v else nid) or nid
    except Exception:
        lab = nid
    lab = f"{ntype}:{str(lab)[:40]}"
    _LABEL_CACHE[key] = lab
    return lab


def ctx_tag(con, q, model):
    """C4 (added on top of C3): the TAG query layer — synthesised SQL result
    (exact relational/aggregation answer) + the source-audit table access."""
    block = []
    # 1) TAG SQL result
    try:
        sql = syn(con, q, model)
        cols, rows = exec_sql(sql)
        block.append("TAG SQL RESULT (synthesised DuckDB query — exact relational answer):")
        block.append(f"  query: {sql}")
        block.append("  rows:\n" + "\n".join("    " + ln
                     for ln in rows_to_text(cols, rows, cap=25).splitlines()))
    except Exception as e:
        block.append(f"TAG SQL RESULT: (query failed: {str(e)[:120]})")
    # 2) source-audit hits relevant to the question (catches fabricated citations)
    kws = _keywords(q)
    if kws:
        likeq = " OR ".join(["lower(citation_key) LIKE ?"] * len(kws))
        pa = [f"%{k}%" for k in kws]
        try:
            aud = con.execute(
                f"SELECT citation_key, verdict, web_verdict, coalesce(web_note,note) "
                f"FROM source_audit WHERE ({likeq}) AND "
                f"(web_verdict IS NOT NULL OR verdict NOT IN ('OK','NO_DOI_FOUND')) LIMIT 8",
                pa).fetchall()
            if aud:
                block.append("SOURCE-AUDIT FLAGS (citation integrity):")
                for ck, v, wv, note in aud:
                    block.append(f"  {ck}: verdict={v} web={wv} — {(note or '')[:140]}")
        except Exception:
            pass
    return "\n".join(block)


# --------------------------------------------------------------------------- #
# Answer generation — identical template; only `context` changes
# --------------------------------------------------------------------------- #
ANS_SYS = (
    "You have NO tools and NO file access — do not attempt to call any tool; "
    "answer in plain prose from the provided context and your own knowledge only. "
    "You are answering questions about a cell-biophysics knowledge base used to "
    "build a mechanistic single-cell simulator. Answer the question directly and "
    "concisely. When the provided context contains the answer, ground your answer "
    "in it and cite specific identifiers (citation_key, KB-x.y, ValidationGate "
    "name, DOI). If you are not certain a cited source is real, say so rather than "
    "presenting it as fact. Do not invent citations, numbers, or papers."
)


def answer(question, context, model):
    if context.strip():
        prompt = (f"CONTEXT:\n{context}\n\n"
                  f"QUESTION: {question}\n\nAnswer:")
    else:
        prompt = (f"QUESTION: {question}\n\n"
                  f"(No knowledge-base context is available — answer from your own "
                  f"knowledge.)\n\nAnswer:")
    return llm(prompt, system=ANS_SYS, model=model, max_tokens=600)


def build_context(cond, con, q, model):
    if cond == "C1":
        return ""
    rag = ctx_rag(q)
    if cond == "C2":
        return rag
    obs = ctx_obsidian(con, q)
    if cond == "C3":
        return "\n\n".join(b for b in (rag, obs) if b)
    tag = ctx_tag(con, q, model)
    # TAG result first: it is C4's decisive, exact evidence (the SQL answer +
    # audit flags), then the graph structure, then the raw PDF passages.
    return "\n\n".join(b for b in (tag, obs, rag) if b)


# --------------------------------------------------------------------------- #
# Scoring — programmatic, auditable
# --------------------------------------------------------------------------- #
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\)\]\}>,;\"']+", re.I)
KBID_RE = re.compile(r"KB-\d+(?:\.\d+)+", re.I)


def _has_all_of(ans_l, group):
    """group is a str (substring) or list (ANY-of)."""
    if isinstance(group, list):
        return any(g.lower() in ans_l for g in group)
    return group.lower() in ans_l


def score_accuracy(g, ans):
    ans_l = ans.lower()
    if not g.get("gold"):
        return None
    return int(all(_has_all_of(ans_l, grp) for grp in g["gold"]))


def score_trap(g, ans):
    """Returns (hallucinated:int, caught:int) for citation-trap questions."""
    ans_l = ans.lower()
    affirms_fake = any(f.lower() in ans_l for f in g.get("fake", []))
    caught = any(d.lower() in ans_l for d in g.get("debunk", []))
    # hallucination = affirmed the fabricated source AND did not flag/correct it
    hall = int(affirms_fake and not caught)
    return hall, int(caught)


def _norm_doi(s):
    """bare, lowercase DOI: strip resolver prefixes and trailing markdown/punct."""
    s = s.lower().strip()
    for pre in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/",
                "http://dx.doi.org/", "doi.org/", "doi:"):
        if s.startswith(pre):
            s = s[len(pre):]
    return s.rstrip(".,;:)]}>*_\"' ")


def score_groundedness(con, ans):
    """Fraction of cited DOIs / KB-ids that resolve to real KB rows.
    Returns (n_cited, n_resolved, n_unresolved). Unresolved cites = hallucination signal.
    Both the cited DOI and the stored DOI are normalised (the DB stores the full
    https://doi.org/… form) — otherwise correct citations look unresolvable."""
    dois = set(_norm_doi(m.group(0)) for m in DOI_RE.finditer(ans))
    kbids = set(m.group(0).upper() for m in KBID_RE.finditer(ans))
    resolved = unresolved = 0
    for d in dois:
        hit = con.execute(
            "SELECT 1 FROM source_evidence WHERE "
            "  replace(replace(replace(lower(doi),'https://doi.org/',''),"
            "  'http://dx.doi.org/',''),'https://dx.doi.org/','')=? "
            "UNION SELECT 1 FROM source_audit WHERE lower(doi)=? OR lower(suggested_doi)=? "
            "LIMIT 1", [d, d, d]).fetchone()
        resolved += 1 if hit else 0
        unresolved += 0 if hit else 1
    for k in kbids:
        hit = con.execute("SELECT 1 FROM knowledge_claim WHERE upper(kb_id)=? LIMIT 1",
                          [k]).fetchone()
        resolved += 1 if hit else 0
        unresolved += 0 if hit else 1
    return len(dois) + len(kbids), resolved, unresolved


# --------------------------------------------------------------------------- #
# Standard AI-eval metric layer (LLM-as-judge, G-Eval style).
# Implements the established RAG/QA-benchmark metrics so the result is comparable
# to how real LLM systems are evaluated:
#   • RAGAS (Es 2023):     answer_relevancy, faithfulness, context_recall
#   • ALCE  (Gao 2023):    citation precision / recall  (computed programmatically
#                          in score_groundedness + here via the gold source)
#   • FActScore (Min 2023) / HaluEval (Li 2023): atomic-fact decomposition ->
#                          fraction unsupported = hallucination rate
#   • TAG-Bench (Biswal 2024): exact-match correctness on relational questions
# The judge is a SEPARATE, stronger model from the answerer (standard practice),
# scores ONLY from the material in its prompt, and returns strict JSON.
# --------------------------------------------------------------------------- #
JUDGE_SYS = (
    "You are a strict, calibrated evaluator for a retrieval-augmented QA benchmark "
    "(RAGAS + FActScore methodology). Judge ONLY from the material in this prompt. "
    "Be skeptical: an answer that asserts a fact not supported by the context or the "
    "reference key is a hallucination, even if it sounds plausible. Return ONLY a "
    "single JSON object, no prose, no code fence."
)

_JSON_RE = re.compile(r"\{.*\}", re.S)


def judge(g, context, ans, model):
    ctx_disp = context[:8000] if context else "(NO CONTEXT — closed-book condition)"
    prompt = (
        f"QUESTION:\n{g['q']}\n\n"
        f"REFERENCE TRUTH (gold key):\n{g['truth']}\n\n"
        f"CONTEXT GIVEN TO THE SYSTEM:\n{ctx_disp}\n\n"
        f"SYSTEM ANSWER:\n{ans}\n\n"
        "Return JSON with exactly these fields:\n"
        '{\n'
        '  "correct": true|false,            // TAG-Bench exact-match: true ONLY if the answer explicitly STATES the gold fact/value. An honest refusal ("I cannot determine") is correct=false (but should have n_hallucinated=0). For trap Qs, true ONLY if it does NOT affirm the fabrication.\n'
        '  "answer_relevancy": 0.0-1.0,       // RAGAS: does the answer actually address the question\n'
        '  "faithfulness": 0.0-1.0,           // RAGAS: fraction of the answer grounded in CONTEXT (use 1.0 if it correctly refuses for lack of context; judge vs world+gold if no context)\n'
        '  "context_recall": 0.0-1.0,         // RAGAS: did the CONTEXT contain the gold answer (0 if no/empty context)\n'
        '  "n_atomic_claims": int,            // FActScore: number of checkable factual claims in the answer\n'
        '  "n_hallucinated": int,             // of those, how many are false or unsupported by CONTEXT+TRUTH\n'
        '  "rationale": "<= 25 words"\n'
        "}"
    )
    raw = pure_llm(prompt, system=JUDGE_SYS, model=model, max_tokens=500)
    m = _JSON_RE.search(raw)
    if not m:
        return dict(correct=None, answer_relevancy=None, faithfulness=None,
                    context_recall=None, n_atomic_claims=None, n_hallucinated=None,
                    rationale="(judge parse fail)", _raw=raw[:200])
    try:
        j = json.loads(m.group(0))
    except Exception:
        return dict(correct=None, answer_relevancy=None, faithfulness=None,
                    context_recall=None, n_atomic_claims=None, n_hallucinated=None,
                    rationale="(judge json fail)", _raw=raw[:200])
    return j


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="KB-architecture ablation benchmark")
    ap.add_argument("--model", default="claude-sonnet-4-6",
                    help="answering LLM (held FIXED across conditions)")
    ap.add_argument("--judge-model", default="claude-opus-4-8",
                    help="evaluator LLM (G-Eval style; stronger than answerer)")
    ap.add_argument("--no-judge", action="store_true",
                    help="skip the LLM-judge standard-metric layer")
    ap.add_argument("--n", type=int, default=None, help="only first N questions")
    ap.add_argument("--conditions", nargs="+", default=CONDITIONS)
    ap.add_argument("--out", default=str(OUT), help="results JSON path")
    ap.add_argument("--questions", default=None,
                    help="JSON file of question dicts (e.g. gen_questions.py output); "
                         "default = the built-in hand-curated GOLD set")
    args = ap.parse_args()
    out_path = pathlib.Path(args.out)

    if not DB_PATH.exists():
        sys.exit(f"{DB_PATH} not found — run notion_to_duckdb.py first.")
    con = duckdb.connect(str(DB_PATH), read_only=True)
    gt = fill_runtime_gold(con)
    print(f"runtime gold: n_verified={gt['n_verified']} n_buckley={gt['n_buckley']}")

    if args.questions:
        loaded = json.loads(pathlib.Path(args.questions).read_text())
        # ensure required keys; generated dicts already carry gold/truth/[fake/debunk]
        questions = loaded
        print(f"loaded {len(questions)} questions from {args.questions}")
    else:
        questions = GOLD
    gold = questions[: args.n] if args.n else questions
    results = []
    t0 = time.time()
    for g in gold:
        for cond in args.conditions:
            ts = time.time()
            ctx = build_context(cond, con, g["q"], args.model)
            ans = answer(g["q"], ctx, args.model)
            acc = score_accuracy(g, ans)
            n_cit, n_res, n_unres = score_groundedness(con, ans)
            rec = dict(qid=g["id"], cls=g["cls"], cond=cond, model=args.model,
                       question=g["q"], context_chars=len(ctx), context=ctx, answer=ans,
                       accuracy=acc, n_cited=n_cit, n_resolved=n_res,
                       n_unresolved=n_unres, dt=round(time.time() - ts, 1))
            if g.get("fake"):              # trap / negative-existence questions
                hall, caught = score_trap(g, ans)
                rec["hallucination"] = hall
                rec["caught_fabrication"] = caught
            else:
                # non-trap hallucination signal = citing an unresolvable identifier
                rec["hallucination"] = int(n_unres > 0)
                rec["caught_fabrication"] = None
            if not args.no_judge:
                rec["judge"] = judge(g, ctx, ans, args.judge_model)
            results.append(rec)
            jm = ""
            if rec.get("judge"):
                jh = rec["judge"].get("n_hallucinated")
                jc = rec["judge"].get("correct")
                jm = f" | judge: correct={jc} halluc_facts={jh}"
            mark = "halluc!" if rec["hallucination"] else ("acc" if acc else "miss")
            print(f"  {g['id']:4s} {cond:3s} {mark:8s} cit={n_cit}({n_unres}bad) "
                  f"{rec['dt']}s  {g['cls']}{jm}")

    out_path.write_text(json.dumps(dict(model=args.model, ground_truth=gt,
                                   n_questions=len(gold), elapsed_s=round(time.time() - t0, 1),
                                   results=results), indent=2))
    print(f"\nwrote {out_path}  ({len(results)} answers, {round(time.time()-t0,1)}s)")

    # console summary by condition
    def _mean(xs):
        xs = [x for x in xs if isinstance(x, (int, float))]
        return sum(xs) / len(xs) if xs else float("nan")

    print("\n=== summary by condition ===")
    print("  [programmatic backbone]            [RAGAS / FActScore (LLM-judge)]")
    print(f"{'cond':18s} {'acc%':>5s} {'hall%':>6s} {'grnd%':>6s} | "
          f"{'EM%':>5s} {'relev':>6s} {'faith':>6s} {'crecall':>7s} {'fact_hall%':>10s}")
    for cond in args.conditions:
        rs = [r for r in results if r["cond"] == cond]
        accs = [r["accuracy"] for r in rs if r["accuracy"] is not None]
        acc = 100 * _mean(accs)
        hall = 100 * _mean([r["hallucination"] for r in rs])
        tot_cit = sum(r["n_cited"] for r in rs)
        gr = 100 * sum(r["n_resolved"] for r in rs) / tot_cit if tot_cit else float("nan")
        js = [r.get("judge", {}) for r in rs if r.get("judge")]
        em = 100 * _mean([1 if j.get("correct") else 0 for j in js]) if js else float("nan")
        rel = _mean([j.get("answer_relevancy") for j in js]) if js else float("nan")
        fai = _mean([j.get("faithfulness") for j in js]) if js else float("nan")
        cre = _mean([j.get("context_recall") for j in js]) if js else float("nan")
        tot_claims = sum((j.get("n_atomic_claims") or 0) for j in js)
        tot_hallf = sum((j.get("n_hallucinated") or 0) for j in js)
        fhall = 100 * tot_hallf / tot_claims if tot_claims else float("nan")
        print(f"{COND_LABEL[cond]:18s} {acc:5.0f} {hall:6.0f} {gr:6.0f} | "
              f"{em:5.0f} {rel:6.2f} {fai:6.2f} {cre:7.2f} {fhall:10.0f}")


if __name__ == "__main__":
    main()
