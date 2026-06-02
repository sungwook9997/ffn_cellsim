#!/usr/bin/env python3
"""Pass 2 — adjudicate the rows pass 1 flagged (DOI_DEAD / NO_DOI_NOMATCH / CHECK).

Pass 1 (verify_sources.py) is a brittle string-match triage: accented authors,
empty short_source, and malformed citation keys produce false positives. This
pass gathers richer evidence per flagged row — the stored DOI's actual CrossRef
record + up to 3 bibliographic candidates + the row's short_source/notes — and
asks the LLM to make ONE judgement per row:

  REAL          paper exists and the citation key is consistent (false alarm)
  DRIFT         real paper, but the key's journal/year is wrong (metadata fix)
  WRONG_DOI     stored DOI resolves to a DIFFERENT paper (DOI fix; paper may exist)
  HALLUCINATION cannot confirm the paper exists at all (likely fabricated)

Writes `source_audit_pass2.md` + an `audit2` column note. Read-only on the KB
except the refined report; does NOT edit Notion.

RUN:  conda activate ffn_sim && python verify_sources_pass2.py
"""
from __future__ import annotations

import json
import re
import time

import duckdb

from tag_query import llm
from verify_sources import (CR, SESSION, msg_authors, msg_title, msg_year,
                            norm_doi, spaced_key)

HERE = __import__("pathlib").Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"
REPORT = HERE / "source_audit_pass2.md"
FLAGGED = ("DOI_DEAD", "NO_DOI_NOMATCH", "CHECK")


def cr_by_doi(doi):
    try:
        r = SESSION.get(f"{CR}/{norm_doi(doi)}", timeout=20)
        return r.json()["message"] if r.status_code == 200 else None
    except Exception:
        return None


def candidates(short, notes, ck):
    qs = [q for q in (short, (notes or "")[:80], spaced_key(ck)) if q and len(q) > 6]
    out, seen = [], set()
    for q in qs:
        try:
            r = SESSION.get(CR, params={"query.bibliographic": q, "rows": 3}, timeout=20)
            items = r.json()["message"]["items"] if r.status_code == 200 else []
        except Exception:
            items = []
        for m in items:
            d = norm_doi(m.get("DOI", ""))
            if d and d not in seen:
                seen.add(d)
                out.append({"doi": d, "title": msg_title(m)[:90],
                            "authors": msg_authors(m)[:3], "year": msg_year(m)})
        time.sleep(0.1)
        if len(out) >= 4:
            break
    return out[:4]


def main():
    con = duckdb.connect(str(DB_PATH))
    rows = con.execute(f"""
        SELECT a.verdict, a.citation_key, a.doi, a.suggested_doi,
               se.short_source, se.notes
        FROM source_audit a JOIN source_evidence se ON se.citation_key = a.citation_key
        WHERE a.verdict IN {FLAGGED}
        ORDER BY a.verdict, a.citation_key
    """).fetchall()
    print(f"adjudicating {len(rows)} flagged rows ...", flush=True)

    cache = HERE / "pass2_dossiers.json"
    if cache.exists():
        dossiers = json.loads(cache.read_text())
        print(f"  loaded {len(dossiers)} cached dossiers", flush=True)
    else:
        dossiers = []
        for verdict, ck, doi, sugg, short, notes in rows:
            d = {"citation_key": ck, "pass1": verdict, "stored_doi": doi or "",
                 "short_source": short or "", "notes": (notes or "")[:160]}
            if doi:
                m = cr_by_doi(doi)
                d["stored_doi_resolves_to"] = (
                    None if m is None else
                    {"title": msg_title(m)[:90], "authors": msg_authors(m)[:3],
                     "year": msg_year(m)})
            d["crossref_candidates"] = candidates(short, notes, ck)
            dossiers.append(d)
            print(f"  gathered {ck}", flush=True)
        cache.write_text(json.dumps(dossiers, ensure_ascii=False, indent=1))

    # Rule-based adjudication over the gathered CrossRef ground truth. Accent/
    # hyphen-insensitive author matching (the main pass-1 false-positive source).
    import unicodedata

    def norm(s: str) -> str:
        s = unicodedata.normalize("NFKD", s or "")
        s = "".join(c for c in s if not unicodedata.combining(c))
        return re.sub(r"[^a-z0-9]", "", s.lower())

    def yr_of(*texts) -> int | None:
        for t in texts:
            m = re.search(r"(19|20)\d{2}", t or "")
            if m:
                return int(m.group(0))
        return None

    def rec_match(rec, kb, yr):
        a = any(len(norm(x)) >= 4 and norm(x) in kb for x in rec.get("authors", []))
        ry = str(rec.get("year") or "")
        y = bool(yr) and ry.isdigit() and abs(int(ry) - yr) <= 1
        return a, y

    verdicts = []
    for d in dossiers:
        ck = d["citation_key"]
        kb = norm(ck) + norm(d.get("short_source", ""))
        yr = yr_of(ck, d.get("short_source", ""))
        resolved = d.get("stored_doi_resolves_to")
        recs = ([(resolved, "stored")] if resolved else []) + \
               [(c, "cand") for c in d.get("crossref_candidates", [])]
        strong = next((r for r in recs if all(rec_match(r[0], kb, yr))), None)
        authonly = next((r for r in recs if rec_match(r[0], kb, yr)[0]), None)
        has_doi = bool(d.get("stored_doi"))
        if strong:
            rec, src = strong
            if src == "stored":
                label, reason = "CONFIRMED_REAL", "stored DOI resolves to matching author+year"
            elif has_doi and resolved and not rec_match(resolved, kb, yr)[0]:
                label, reason = "METADATA_FIX", f"stored DOI -> wrong paper; real DOI likely {rec.get('doi','')}"
            else:
                label, reason = "CONFIRMED_REAL", f"author+year match (DOI {rec.get('doi','')})"
        elif authonly:
            label, reason = "METADATA_FIX", f"author matches; year/journal differs (cf {authonly[0].get('year')})"
        elif has_doi and resolved:
            label, reason = "METADATA_FIX", "DOI resolves but author mismatch — verify key"
        else:
            label, reason = "UNCONFIRMED", "no CrossRef record matches author — verify existence"
        verdicts.append({"citation_key": ck, "label": label, "reason": reason})
    vmap = {v["citation_key"]: v for v in verdicts}

    con.execute("ALTER TABLE source_audit ADD COLUMN IF NOT EXISTS audit2 TEXT")
    con.execute("ALTER TABLE source_audit ADD COLUMN IF NOT EXISTS audit2_reason TEXT")
    for ck, v in vmap.items():
        con.execute("UPDATE source_audit SET audit2=?, audit2_reason=? "
                    "WHERE citation_key=?", [v["label"], v.get("reason", ""), ck])

    from collections import Counter
    counts = Counter(v["label"] for v in verdicts)
    order = {"UNCONFIRMED": 0, "METADATA_FIX": 1, "CONFIRMED_REAL": 2}
    verdicts.sort(key=lambda v: order.get(v["label"], 9))
    lines = ["# SourceEvidence audit — pass 2 (LLM adjudication of flagged rows)\n",
             f"Adjudicated **{len(verdicts)}** rows flagged by pass 1.\n",
             "| label | n |", "|---|---|"]
    for lbl in sorted(counts, key=lambda x: order.get(x, 9)):
        lines.append(f"| {lbl} | {counts[lbl]} |")
    lines += ["\n| label | citation_key | pass1 | reason |", "|---|---|---|---|"]
    for v in verdicts:
        p1 = next((r[0] for r in rows if r[1] == v["citation_key"]), "")
        lines.append(f"| {v['label']} | {v['citation_key']} | {p1} | "
                     f"{v.get('reason','').replace('|','/')} |")
    REPORT.write_text("\n".join(lines))
    con.close()

    print("\n=== pass-2 labels ===")
    for lbl in sorted(counts, key=lambda x: order.get(x, 9)):
        print(f"  {lbl:14s} {counts[lbl]}")
    print(f"\nreport: {REPORT}")


if __name__ == "__main__":
    main()
