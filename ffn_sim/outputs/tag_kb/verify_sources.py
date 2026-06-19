#!/usr/bin/env python3
"""Hallucination audit of the SourceEvidence table against CrossRef.

The 243 SourceEvidence rows were built by prior LLM sessions (KB migration +
deep-research dossiers) WITHOUT the actual PDFs, so some may be hallucinated:
fabricated citations, wrong/dead DOIs, or a real DOI pointing at a different
paper than the citation key claims. This is the cheap first pass — it checks
*existence + metadata*, not claim support (that needs the PDF, a later pass).

For each SourceEvidence row:
  • has DOI  -> GET api.crossref.org/works/{doi}; compare the resolved paper's
      author-surname + year against the citation key / short source.
  • no DOI   -> bibliographic search; report the best candidate + a suggested DOI.

Verdicts (suspicion-ranked):
  DOI_DEAD        DOI does not resolve (404/None)            <- high suspicion
  DOI_MISMATCH    resolves but author AND year both differ   <- high suspicion
  NO_DOI_NOMATCH  no DOI and no confident CrossRef match     <- high suspicion
  CHECK           resolves, only author XOR year matches      <- review
  NO_DOI_FOUND    no DOI but a real paper matches (suggests DOI) <- low (just add DOI)
  OK              resolves, author + year both match          <- verified

Writes the `source_audit` table into kb.duckdb + `source_audit_report.md`.

RUN:
  conda activate ffn_sim
  python verify_sources.py              # audits all rows
  python verify_sources.py --limit 20   # quick sample

Limited runs write source_audit_report_limit<N>.md and source_audit_sample,
leaving the canonical source_audit_report.md/source_audit table untouched.
"""
from __future__ import annotations

import re
import sys
import time

import duckdb
import requests

HERE = __import__("pathlib").Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"
REPORT = HERE / "source_audit_report.md"
CR = "https://api.crossref.org/works"
UA = "ffn_cellsim-source-audit/1.0 (mailto:sungwook999@gmail.com)"
SESSION = requests.Session()
SESSION.headers["User-Agent"] = UA

SUSPICION = {"DOI_DEAD": 0, "DOI_MISMATCH": 1, "NO_DOI_NOMATCH": 2,
             "CHECK": 3, "NO_DOI_FOUND": 4, "OK": 5}


def norm_doi(s: str | None) -> str:
    if not s:
        return ""
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", s.strip().lower()).rstrip(".)")


def year_of(*texts: str) -> str:
    for t in texts:
        m = re.search(r"\b(19|20)\d{2}\b", t or "")
        if m:
            return m.group(0)
    return ""


def cr_get_doi(doi: str):
    try:
        r = SESSION.get(f"{CR}/{doi}", timeout=20)
        if r.status_code == 404:
            return None
        r.raise_for_status()
        return r.json()["message"]
    except Exception:
        return None


def cr_search(query: str, rows: int = 5) -> list:
    try:
        r = SESSION.get(CR, params={"query.bibliographic": query, "rows": rows},
                        timeout=20)
        r.raise_for_status()
        return r.json()["message"]["items"]
    except Exception:
        return []


def spaced_key(ck: str) -> str:
    """'BroederszMacKintosh2014_RMP' -> 'Broedersz Mac Kintosh 2014' (drop journal
    token) so a bibliographic search on concatenated author keys can match."""
    head = ck.split("_")[0]                                   # drop _journal
    head = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", head)          # camel -> spaces
    head = re.sub(r"(?<=[A-Za-z])(?=\d)", " ", head)          # author|year split
    return head.strip()


def msg_year(m: dict) -> str:
    for k in ("published", "published-print", "published-online", "issued"):
        try:
            y = m[k]["date-parts"][0][0]
            if y is not None and str(y).isdigit():
                return str(y)
        except Exception:
            pass
    return ""


def msg_authors(m: dict) -> list[str]:
    return [a.get("family", "").lower() for a in m.get("author", []) if a.get("family")]


def msg_title(m: dict) -> str:
    t = m.get("title") or [""]
    return t[0] if t else ""


def author_hit(authors: list[str], *texts: str) -> bool:
    blob = " ".join(t.lower() for t in texts if t)
    return any(a and a in blob for a in authors)


def audit_row(uid, ck, doi, short, notes):
    nd = norm_doi(doi)
    yr_ck = year_of(ck, short)
    if nd:
        m = cr_get_doi(nd)
        if m is None:
            return ("DOI_DEAD", nd, "", "DOI did not resolve on CrossRef")
        cr_y, cr_a, cr_t = msg_year(m), msg_authors(m), msg_title(m)
        a_ok = author_hit(cr_a, ck, short)
        y_ok = bool(yr_ck) and abs(int(yr_ck) - int(cr_y or 0)) <= 1 if cr_y else False
        if a_ok and y_ok:
            v = "OK"
        elif a_ok or y_ok:
            v = "CHECK"
        else:
            v = "DOI_MISMATCH"
        note = f"CrossRef: {(cr_a[0] if cr_a else '?')} {cr_y} — {cr_t[:70]}"
        return (v, nd, "", note)
    # no DOI -> bibliographic search; scan top candidates across two query forms
    queries = [q for q in (short, spaced_key(ck)) if q and len(q) > 6]
    strong = weak = None
    for q in queries:
        for m in cr_search(q):
            cr_y, cr_a = msg_year(m), msg_authors(m)
            a_ok = author_hit(cr_a, ck, short)
            y_ok = bool(yr_ck) and bool(cr_y) and abs(int(yr_ck) - int(cr_y)) <= 1
            if a_ok and y_ok:
                strong = m
                break
            if a_ok and weak is None:
                weak = m
        if strong:
            break
    hit = strong or weak
    if hit is None:
        return ("NO_DOI_NOMATCH", "", "", "no CrossRef author match in top hits")
    cr_y, cr_a, cr_t = msg_year(hit), msg_authors(hit), msg_title(hit)
    cand = norm_doi(hit.get("DOI", ""))
    note = f"candidate: {(cr_a[0] if cr_a else '?')} {cr_y} — {cr_t[:60]}"
    if strong is not None:
        return ("NO_DOI_FOUND", "", cand, note)
    return ("CHECK", "", cand, "author-only match (year differs); " + note)


def check() -> int:
    """Non-destructive audit-status read for refresh.sh — NO web calls, NO writes.

    Reports the materialised `source_audit` verdict distribution, lists suspect
    rows, and warns if the audit is partial/stale vs the current source_evidence
    count (so a citation-integrity check rides every refresh; the slow CrossRef
    re-verification stays an explicit `verify_sources.py` run). Always exits 0.
    """
    # Robust to an absent / un-materialized DB: the citation audit reads from
    # kb.duckdb (Notion -> duckdb), which is gitignored and absent on a fresh
    # clone / CI / any machine without a Notion pull. Fail OPEN (skip, exit 0) so
    # `make kb-check` and refresh.sh never hard-crash where the DB isn't built;
    # CI enforcement of the production-cited subset rides verify_params instead.
    if not DB_PATH.exists():
        print("  [audit] kb.duckdb absent — run refresh.sh (Notion->duckdb) to "
              "materialize source_evidence before the citation audit; skipped.")
        return 0
    con = duckdb.connect(str(DB_PATH), read_only=True)

    def _has_table(name: str) -> bool:
        return con.execute("SELECT count(*) FROM information_schema.tables "
                           "WHERE table_name=?", [name]).fetchone()[0] > 0

    if not _has_table("source_evidence"):
        print("  [audit] source_evidence not materialized in kb.duckdb — run "
              "refresh.sh (Notion=SoT) first; citation audit skipped.")
        con.close()
        return 0
    has = _has_table("source_audit")
    n_se = con.execute("SELECT count(*) FROM source_evidence").fetchone()[0]
    if not has:
        print("  [audit] source_audit table MISSING — run `python verify_sources.py` "
              f"to audit all {n_se} SourceEvidence rows.")
        con.close()
        return 0
    rows = con.execute("SELECT verdict, count(*) FROM source_audit "
                       "GROUP BY 1 ORDER BY 2 DESC").fetchall()
    n_aud = con.execute("SELECT count(*) FROM source_audit").fetchone()[0]
    print(f"  [audit] source_audit: {n_aud}/{n_se} SourceEvidence rows audited")
    for verdict, n in rows:
        flag = "  <-- SUSPECT" if SUSPICION.get(verdict, 9) <= 3 else ""
        print(f"           {verdict:16s} {n}{flag}")
    suspects = con.execute(
        "SELECT citation_key, verdict FROM source_audit "
        "WHERE verdict IN ('DOI_DEAD','DOI_MISMATCH','NO_DOI_NOMATCH','CHECK') "
        "ORDER BY citation_key").fetchall()
    for ck, v in suspects:
        print(f"           ⚠️  {v}: {ck}")
    if n_aud < n_se:
        print(f"  [audit] ⚠️ PARTIAL/STALE — {n_se - n_aud} rows unaudited "
              "(likely a --limit sample). Run full `python verify_sources.py`.")
    con.close()
    return 0


def main():
    if "--check" in sys.argv:
        sys.exit(check())
    limit = None
    if "--limit" in sys.argv:
        limit = int(sys.argv[sys.argv.index("--limit") + 1])
    partial = limit is not None
    report_path = REPORT if not partial else HERE / f"source_audit_report_limit{limit}.md"
    table_name = "source_audit" if not partial else "source_audit_sample"

    con = duckdb.connect(str(DB_PATH))
    rows = con.execute(
        "SELECT uid, citation_key, doi, short_source, notes FROM source_evidence "
        "ORDER BY TRY_CAST(regexp_replace(uid, '[^0-9]', '', 'g') AS INTEGER) "
        "NULLS LAST, citation_key").fetchall()
    if limit is not None:
        rows = rows[:limit]

    results = []
    for i, (uid, ck, doi, short, notes) in enumerate(rows, 1):
        verdict, used_doi, suggest, note = audit_row(uid, ck, doi, short, notes)
        results.append((uid, ck, doi or "", verdict, suggest, note))
        flag = "" if verdict in ("OK", "NO_DOI_FOUND") else "  <-- SUSPECT"
        print(f"[{i:3d}/{len(rows)}] {verdict:14s} {ck[:38]:38s}{flag}", flush=True)
        time.sleep(0.12)  # polite to CrossRef

    con.execute(f"DROP TABLE IF EXISTS {table_name}")
    con.execute(f"CREATE TABLE {table_name} (uid TEXT, citation_key TEXT, "
                "doi TEXT, verdict TEXT, suggested_doi TEXT, note TEXT)")
    con.executemany(f"INSERT INTO {table_name} VALUES (?,?,?,?,?,?)", results)

    # report
    from collections import Counter
    counts = Counter(r[3] for r in results)
    results.sort(key=lambda r: (SUSPICION.get(r[3], 9), r[1]))
    lines = ["# SourceEvidence hallucination audit (CrossRef)\n",
             f"Audited **{len(results)}** SourceEvidence rows.\n", "## Summary\n",
             "| verdict | n | meaning |", "|---|---|---|"]
    if partial:
        lines.insert(
            2,
            f"> Sample run from `--limit {limit}`. Canonical source_audit_report.md "
            "and source_audit were not overwritten.\n",
        )
    meaning = {"DOI_DEAD": "DOI does not resolve — likely fabricated/wrong",
               "DOI_MISMATCH": "DOI resolves to a DIFFERENT paper (author+year both off)",
               "NO_DOI_NOMATCH": "no DOI and no confident match — unverifiable",
               "CHECK": "partial match (author XOR year) — review",
               "NO_DOI_FOUND": "real paper found — just missing DOI (suggested)",
               "OK": "DOI resolves, author+year match — verified"}
    for v in sorted(counts, key=lambda v: SUSPICION.get(v, 9)):
        lines.append(f"| {v} | {counts[v]} | {meaning.get(v,'')} |")
    susp = [r for r in results if r[3] in ("DOI_DEAD", "DOI_MISMATCH", "NO_DOI_NOMATCH")]
    lines += [f"\n**{len(susp)} high-suspicion rows** (top of list).\n",
              "## All rows (suspicion-ranked)\n",
              "| verdict | citation_key | DOI | suggested | CrossRef note |",
              "|---|---|---|---|---|"]
    for uid, ck, doi, verdict, suggest, note in results:
        lines.append(f"| {verdict} | {ck} | {doi[:40]} | {suggest} | {note.replace('|','/')} |")
    report_path.write_text("\n".join(lines))
    con.close()

    print("\n=== verdict counts ===")
    for v in sorted(counts, key=lambda v: SUSPICION.get(v, 9)):
        print(f"  {v:14s} {counts[v]}")
    print(f"\nhigh-suspicion: {len(susp)}  ->  report: {report_path}")
    print(f"duckdb table: {table_name}")


if __name__ == "__main__":
    main()
