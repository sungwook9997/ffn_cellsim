#!/usr/bin/env python
"""Build the KB-upload payload from the 2026-07-15 missing-inventory audit (wf_072fc794).

Reads the workflow journal (inventory + KB-check + source-verify agent results) and emits
a structured `kb_ingest_payload_2026-07-15.json` with two lists:

  se_new[]   {citation_key, doi, title, authors_year, doi_confidence, gap, compartment}
  ku_draft[] {prov_kb_id, unit_partition, title, statement, value_range_si, params_text,
              confidence, citations_text, gap, compartment, collision_flag}

Nothing is pushed here — this only assembles the reviewable payload. The pusher
`kb_ingest_audit_2026-07-15.py` CrossRef-verifies every DOI (auto-rejects the hallucinated
ones the manifest flagged for DISCARD) and de-dups against the live Notion SoT.

Provisional KB IDs are assigned per compartment→Unit-partition next-free index; genuine
collisions (>1 gap wanting the same next index) are flagged for PI renumber at ingest.
"""
from __future__ import annotations
import json, re, sys
from pathlib import Path

JOURNAL = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/Users/sw1/.claude/projects/-Users-sw1-ffn-cellsim/"
    "5a116646-948c-473f-a664-87a30b8d22af/subagents/workflows/wf_072fc794-e2e/journal.jsonl")
OUT = Path(__file__).with_name("kb_ingest_payload_2026-07-15.json")

# Compartment name -> KnowledgeClaim "Unit" select-partition (from the live KU schema)
UNIT = {
    "membrane": "Unit 3.B Compartment", "cortex": "Unit 3 Cell", "actin": "Unit 3 Cell",
    "microtub": "Unit 3 Cell", "intermediate": "Unit 3.B Compartment", "myosin": "Unit 3 Cell",
    "cytoplasm": "Unit 3.B Compartment", "cfd": "Unit 3.B Compartment", "nucleus": "Unit 3.B Compartment",
    "volume": "Unit 3.B Compartment", "osmotic": "Unit 3.B Compartment",
    "ecm-adhesion": "Unit 2 ECM-Cell", "adhesion": "Unit 2 ECM-Cell",
    "junction": "Unit 4 Junction", "cell–cell": "Unit 4 Junction",
    "ecm": "Unit 1 ECM", "dynamic": "Unit 3 Cell", "remodel": "Unit 3 Cell",
    "hazard": "Unit 3 Cell", "integration": "Unit 3 Cell",
}
def unit_for(comp: str) -> str:
    c = comp.lower()
    for k, v in UNIT.items():
        if k in c:
            return v
    return "Unit 3 Cell"

def load_results(path: Path):
    inv, sv, mani = [], [], None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            o = json.loads(line)
        except Exception:
            continue
        if o.get("type") != "result":
            continue
        r = o.get("result")
        if isinstance(r, dict) and "gaps" in r and any("kb_status" in g for g in r.get("gaps", [])):
            inv.append(r)   # KB-check dicts (augmented with kb_status)
        elif isinstance(r, dict) and "doi_confidence" in r:
            sv.append(r)    # source-verify dicts
        elif isinstance(r, str) and len(r) > 2000 and "SourceEvidence" in r:
            mani = r        # the manifest synthesis text (curated §1 tables)
    seen = {}
    for r in inv:
        seen[r["compartment"]] = r
    return list(seen.values()), sv, mani


def parse_manifest_se(mani: str):
    """Parse the curated §1 SourceEvidence tables from the manifest.

    Rows: | Paper | DOI | Status | Unblocks |. The manifest already applied the
    citation-integrity corrections (the DISCARD replacements live in §1 as the RIGHT
    paper), so parsing §1 is the safe SE source — NOT regex over raw agent prose.
    FIX rows (metadata backfill to an EXISTING SE) are collected separately, not created.
    """
    i1, i2 = mani.find("# (1)"), mani.find("# (2)")
    block = mani[i1:i2] if i1 >= 0 and i2 > i1 else ""
    comp, se, fixes = "", [], []
    for line in block.splitlines():
        s = line.strip()
        if s.startswith("### "):
            comp = s[4:].strip()
            continue
        if not s.startswith("|") or s.startswith("|---") or s.lower().startswith("| paper"):
            continue
        cols = [c.strip() for c in s.strip("|").split("|")]
        if len(cols) < 3:
            continue
        paper, doi_cell, status = cols[0], cols[1], cols[2]
        unblocks = cols[3] if len(cols) > 3 else ""
        m = DOI_TOKEN.search(doi_cell)
        if not m:
            continue
        doi = strip_doi(m.group(0))
        st = status.upper()
        if "DISCARD" in st:
            continue                                  # never upload
        paper_clean = re.sub(r"\*+|SE\d+\s*", "", paper).strip()
        aym = re.match(r"([A-Za-zÀ-ſ'’&/-]+)[^0-9]*((?:19|20)\d\d)", paper_clean)
        key = (re.sub(r"[^A-Za-z]", "", aym.group(1)) + aym.group(2)) if aym else make_key("", paper_clean)
        row = {"citation_key": key, "doi": doi, "title": paper_clean[:180],
               "authors_year": paper_clean, "doi_confidence": "HIGH" if "HIGH" in st else st.strip(" []"),
               "gap": unblocks[:120], "compartment": comp}
        if "FIX" in st:
            fixes.append(row)                          # update to an existing SE — report, don't create
        else:
            se.append(row)
    return se, fixes

def strip_doi(tok: str) -> str:
    """Clean a DOI token: keep balanced internal parens (old Elsevier DOIs like
    10.1016/S0006-3495(99)77040-2), strip trailing punctuation + any UNbalanced
    trailing ')' that came from prose like '(see 10.x/y)'."""
    d = tok.strip().rstrip(".,;]}").lower()
    while d.endswith(")") and d.count(")") > d.count("("):
        d = d[:-1].rstrip(".,;]}")
    return d

# DOIs end at whitespace/semicolon; internal () are allowed (stripped if unbalanced)
DOI_TOKEN = re.compile(r"10\.\d{4,9}/[^\s;|*]+")


def make_key(authors_year: str, gap: str) -> str:
    m = re.match(r"([A-Za-zÀ-ſ'-]+)", (authors_year or "").strip())
    a = (m.group(1) if m else re.sub(r"[^A-Za-z]", "", (gap or "X"))[:8]) or "Ref"
    y = re.search(r"(19|20)\d\d", authors_year or "")
    return f"{a}{y.group(0) if y else ''}"

def main():
    inv, sv, mani = load_results(JOURNAL)
    print(f"loaded {len(inv)} compartments, {len(sv)} source-verify rows, manifest={'yes' if mani else 'NO'}")

    # ---- SourceEvidence candidates ----
    # (1) curated manifest §1 tables (compartments 1-2, corrections applied) — authoritative.
    # (2) source-verify PRIMARY DOI per gap for the compartments §1 didn't tabulate — take
    #     ONLY the first clean DOI in the `doi` field, and SKIP it if a reject-keyword
    #     ("reject/wrong/invalid/misattribut/not resolve/substitut/corrected from") precedes
    #     it, so we never pick up a DOI the agent explicitly rejected. CrossRef gates the rest.
    se_raw, fixes = parse_manifest_se(mani or "")
    se, seen_doi = [], set()
    for r in se_raw:
        if r["doi"] in seen_doi:
            continue
        seen_doi.add(r["doi"])
        se.append(r)

    REJECT = re.compile(r"(reject|wrong|invalid|misattribut|not\s+resolve|does\s+not|substitut|corrected\s+from|instead\s+of)", re.I)
    for r in sv:
        conf = r.get("doi_confidence", "")
        if conf == "UNVERIFIED":
            continue
        field = str(r.get("doi", ""))
        m = DOI_TOKEN.search(field)
        if not m:
            continue
        # skip if a reject-keyword appears before this first DOI
        if REJECT.search(field[:m.start()]):
            continue
        d = strip_doi(m.group(0))
        if d in seen_doi:
            continue
        seen_doi.add(d)
        se.append({
            "citation_key": make_key(r.get("authors_year", ""), r.get("gap", "")),
            "doi": d,
            "title": (r.get("source_title") or "")[:180],
            "authors_year": r.get("authors_year", ""),
            "doi_confidence": conf,
            "gap": r.get("gap", ""),
            "compartment": r.get("compartment", ""),
        })

    # ---- KnowledgeClaim drafts (from inventory gaps needing a KU) ----
    NEED = {"NOT-IN-TAG", "IN-CORPUS-NOT-MATERIALIZED", "PROJECT-INTERNAL-ONLY"}
    ku, next_idx = [], {}
    for comp in inv:
        cname = comp["compartment"]
        for g in comp.get("gaps", []):
            ua = (g.get("upload_action") or "").lower()
            if g.get("kb_status") not in NEED:
                continue
            if "knowledgeclaim" not in ua and "ku" not in ua and "knowledge claim" not in ua:
                continue
            u = unit_for(cname)
            n = next_idx.get(u, 1)
            next_idx[u] = n + 1
            params = "; ".join(
                f"{p.get('param','')}={p.get('value','')}{(' '+p['unit']) if p.get('unit') else ''}".strip()
                for p in (g.get("key_params") or []))
            ku.append({
                "prov_kb_id": f"KB-DRAFT-{u.split()[1]}-{n:02d}",   # provisional, PI renumbers
                "unit_partition": u,
                "title": (g.get("name") or "")[:180],
                "statement": (g.get("correct_model") or "")[:1800],
                "value_range_si": params[:1800],
                "params_text": params,
                "confidence": "Medium",
                "citations_text": (g.get("lit_source") or "")[:400],
                "gap": g.get("name", ""),
                "compartment": cname,
                "kb_status": g.get("kb_status"),
                "collision_flag": "PROVISIONAL id — PI assigns final KB number at ingest",
            })

    payload = {"generated": "2026-07-15", "source": "wf_072fc794 missing-inventory audit",
               "se_new": se, "ku_draft": ku, "se_fixes": fixes}
    OUT.write_text(json.dumps(payload, indent=1, ensure_ascii=False))
    print(f"wrote {OUT}")
    print(f"  se_new: {len(se)} papers (curated manifest §1, DISCARD-excluded; pusher CrossRef-gates + de-dups)")
    print(f"  se_fixes: {len(fixes)} metadata backfills to EXISTING SE (reported, not auto-created)")
    print(f"  ku_draft: {len(ku)} claims (status=draft, provisional ids per Unit partition)")
    from collections import Counter
    print("  KU by Unit:", dict(Counter(k["unit_partition"] for k in ku)))
    print("  SE by confidence:", dict(Counter(s["doi_confidence"] for s in se)))

if __name__ == "__main__":
    main()
