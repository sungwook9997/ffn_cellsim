#!/usr/bin/env python3
"""Ingest the references/ PDF corpus into the TAG content layer.

Adds two tables to `kb.duckdb` (built first by notion_to_duckdb.py):
  • references   — one row per PDF (citation_key, title, doi, path, n_pages,
      sha1, source, se_uid, se_citation_key, n_chunks). se_uid/se_citation_key
      are filled when the PDF's DOI matches an existing SourceEvidence row.
  • paper_chunks — page-level text chunks (chunk_id, citation_key, se_uid, page,
      text) + a DuckDB FTS (BM25) index, so TAG `gen` can pull evidence text.

Also writes `references/tag_corpus.json` as the local PDF↔metadata↔SE index.

Safe by construction: reads the PDFs, never deletes them; exact byte-duplicate
files (same sha1) are skipped (first path wins) and reported.

RUN (after notion_to_duckdb.py):
  conda activate ffn_sim
  python references_ingest.py
  python references_ingest.py --stats
  python references_ingest.py --check
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import sys

import duckdb
import fitz  # PyMuPDF

HERE = pathlib.Path(__file__).parent
DB_PATH = HERE / "kb.duckdb"
REF_DIR = pathlib.Path("/Users/sw1/ffn_cellsim/aleph/references")
MANIFEST = REF_DIR / "analysis" / "_manifest.json"
CORPUS_JSON = REF_DIR / "tag_corpus.json"

DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+")
CHUNK_CHARS = 1800


def sha1(p: pathlib.Path) -> str:
    h = hashlib.sha1()
    h.update(p.read_bytes())
    return h.hexdigest()


def norm_doi(s: str | None) -> str:
    if not s:
        return ""
    s = s.strip().lower()
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", s)
    return s.rstrip(".)").strip()


def slug(s: str, n: int = 60) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "-", s).strip("-").lower()
    return s[:n] or "untitled"


def pdf_paths() -> list[pathlib.Path]:
    paths = sorted(REF_DIR.glob("*.pdf"))
    for sub in ("cellpress_bundle", "downloaded", "2026_07_10"):
        d = REF_DIR / sub
        if d.exists():
            paths += sorted(d.glob("*.pdf"))
    return paths


def load_titles() -> dict[str, str]:
    if not MANIFEST.exists():
        return {}
    return {e["path"]: e.get("title", "") for e in json.loads(MANIFEST.read_text())}


def first_doi(text: str) -> str:
    m = DOI_RE.search(text)
    return norm_doi(m.group(0)) if m else ""


def page_chunks(text: str) -> list[str]:
    text = re.sub(r"[ \t]+", " ", text).strip()
    if len(text) <= CHUNK_CHARS:
        return [text] if text else []
    return [text[i:i + CHUNK_CHARS] for i in range(0, len(text), CHUNK_CHARS)]


def existing_corpus_link_count() -> int | None:
    if not CORPUS_JSON.exists():
        return None
    try:
        rows = json.loads(CORPUS_JSON.read_text())
    except Exception:
        return None
    if not isinstance(rows, list):
        return None
    return sum(1 for row in rows if row.get("linked_to_source_evidence"))


def assert_no_link_regression(corpus: list[dict], *, allow_link_drop: bool) -> None:
    previous = existing_corpus_link_count()
    if previous is None:
        return
    current = sum(1 for row in corpus if row.get("linked_to_source_evidence"))
    if current >= previous or allow_link_drop:
        return
    raise SystemExit(
        "Refusing to overwrite tag_corpus.json because linked_to_source_evidence "
        f"would decrease ({previous} -> {current}). Re-run with --allow-link-drop "
        "only if this source-link drop is intentional."
    )


def table_exists(con, name: str) -> bool:
    return bool(
        con.execute(
            "SELECT count(*) FROM information_schema.tables WHERE table_name = ?",
            [name],
        ).fetchone()[0]
    )


def count_rows(con, table: str, where: str | None = None) -> int:
    if not table_exists(con, table):
        return 0
    sql = f"SELECT count(*) FROM {table}"
    if where:
        sql += f" WHERE {where}"
    return int(con.execute(sql).fetchone()[0])


def sanity_summary() -> None:
    if not DB_PATH.exists():
        sys.exit(f"{DB_PATH} not found — run notion_to_duckdb.py first.")
    con = duckdb.connect(str(DB_PATH), read_only=True)
    se_rows = count_rows(con, "source_evidence")
    kc_rows = count_rows(con, "knowledge_claim")
    ref_rows = count_rows(con, "paper_refs")
    chunk_rows = count_rows(con, "paper_chunks")
    linked_refs = count_rows(con, "paper_refs", "se_uid IS NOT NULL")
    orphan_refs = count_rows(con, "paper_refs", "se_uid IS NULL")
    if table_exists(con, "source_evidence") and table_exists(con, "paper_refs"):
        orphan_sources = int(
            con.execute(
                """
                SELECT count(*)
                FROM source_evidence se
                LEFT JOIN (
                  SELECT DISTINCT se_uid FROM paper_refs WHERE se_uid IS NOT NULL
                ) pr ON pr.se_uid = se.uid
                WHERE pr.se_uid IS NULL
                """
            ).fetchone()[0]
        )
    else:
        orphan_sources = 0
    con.close()

    corpus_total = corpus_linked = None
    if CORPUS_JSON.exists():
        try:
            corpus = json.loads(CORPUS_JSON.read_text())
            if isinstance(corpus, list):
                corpus_total = len(corpus)
                corpus_linked = sum(1 for row in corpus if row.get("linked_to_source_evidence"))
        except Exception:
            pass

    print("TAG sanity summary")
    print(f"  SourceEvidence rows: {se_rows}")
    print(f"  KnowledgeClaim rows: {kc_rows}")
    print(f"  paper_refs rows: {ref_rows}")
    print(f"  paper_chunks rows: {chunk_rows}")
    print(f"  linked paper_refs: {linked_refs}")
    print(f"  unlinked/orphan paper_refs: {orphan_refs}")
    print(f"  SourceEvidence rows without linked PDF: {orphan_sources}")
    if corpus_total is not None:
        print(f"  tag_corpus.json linked entries: {corpus_linked}/{corpus_total}")


def ingest(*, allow_link_drop: bool = False):
    if not DB_PATH.exists():
        sys.exit(f"{DB_PATH} not found — run notion_to_duckdb.py first.")
    con = duckdb.connect(str(DB_PATH))

    # SourceEvidence DOI -> (uid, citation_key) for dedup linking
    se_by_doi: dict[str, tuple] = {}
    for uid, ck, doi in con.execute(
            "SELECT uid, citation_key, doi FROM source_evidence").fetchall():
        nd = norm_doi(doi)
        if nd:
            se_by_doi[nd] = (uid, ck)

    titles = load_titles()
    seen_sha: dict[str, str] = {}
    ref_rows, chunk_rows, corpus = [], [], []
    skipped_dups = []

    for p in pdf_paths():
        sha = sha1(p)
        if sha in seen_sha:
            skipped_dups.append((str(p.name), seen_sha[sha]))
            continue
        seen_sha[sha] = p.name
        rel = p.name if p.parent == REF_DIR else f"{p.parent.name}/{p.name}"
        source = "new" if p.parent == REF_DIR else p.parent.name
        try:
            doc = fitz.open(p)
        except Exception as e:
            print(f"  !! open failed {rel}: {e}")
            continue
        pages_text = [doc[i].get_text() for i in range(doc.page_count)]
        n_pages = doc.page_count
        doc.close()
        full = "\n".join(pages_text)
        doi = first_doi(full)
        title = (titles.get(p.name) or "").strip()
        if not title:                                  # fall back to first text line
            for ln in pages_text[0].splitlines() if pages_text else []:
                if len(ln.strip()) > 15:
                    title = ln.strip()
                    break
        se_uid, se_ck = se_by_doi.get(doi, ("", ""))
        ck = se_ck or f"{slug(title or p.stem, 40)}-{sha[:6]}"

        nc = 0
        for pi, ptext in enumerate(pages_text):
            for ci, ch in enumerate(page_chunks(ptext)):
                # include the file sha so two distinct PDFs that resolve to the
                # SAME citation_key (e.g. publisher + preprint of one DOI) do not
                # collide on chunk_id — a duplicate PK breaks the BM25 FTS macro.
                chunk_rows.append((f"{ck}:{sha[:6]}:p{pi}:{ci}", ck, se_uid or None, pi, ch))
                nc += 1
        ref_rows.append((ck, title, doi or None, rel, n_pages, sha, source,
                         se_uid or None, se_ck or None, nc))
        corpus.append({"citation_key": ck, "title": title, "doi": doi,
                       "path": rel, "n_pages": n_pages, "source": source,
                       "se_uid": se_uid, "se_citation_key": se_ck, "n_chunks": nc,
                       "linked_to_source_evidence": bool(se_uid)})

    try:
        assert_no_link_regression(corpus, allow_link_drop=allow_link_drop)
    except SystemExit:
        con.close()
        raise

    con.execute("DROP TABLE IF EXISTS paper_chunks")
    con.execute("DROP TABLE IF EXISTS paper_refs")
    con.execute(
        "CREATE TABLE paper_refs (citation_key TEXT, title TEXT, doi TEXT, "
        "path TEXT, n_pages INTEGER, sha1 TEXT, source TEXT, se_uid TEXT, "
        "se_citation_key TEXT, n_chunks INTEGER)")
    con.execute(
        "CREATE TABLE paper_chunks (chunk_id TEXT, citation_key TEXT, "
        "se_uid TEXT, page INTEGER, text TEXT)")

    con.executemany("INSERT INTO paper_refs VALUES (?,?,?,?,?,?,?,?,?,?)", ref_rows)
    con.executemany("INSERT INTO paper_chunks VALUES (?,?,?,?,?)", chunk_rows)
    # BM25 full-text index for the semantic content step
    con.execute("INSTALL fts; LOAD fts;")
    con.execute("PRAGMA create_fts_index('paper_chunks', 'chunk_id', 'text', "
                "overwrite=1)")
    con.close()

    CORPUS_JSON.write_text(json.dumps(corpus, indent=2, ensure_ascii=False))
    linked = sum(1 for c in corpus if c["linked_to_source_evidence"])
    print(f"ingested {len(ref_rows)} PDFs -> {len(chunk_rows)} chunks "
          f"({linked} linked to SourceEvidence by DOI)")
    if skipped_dups:
        print(f"skipped {len(skipped_dups)} exact byte-duplicate(s):")
        for dup, orig in skipped_dups:
            print(f"  {dup}  == {orig}")
    print(f"local index -> {CORPUS_JSON}")


def stats():
    con = duckdb.connect(str(DB_PATH), read_only=True)
    nref = con.execute("SELECT count(*) FROM paper_refs").fetchone()[0]
    nch = con.execute("SELECT count(*) FROM paper_chunks").fetchone()[0]
    nlink = con.execute("SELECT count(*) FROM paper_refs WHERE se_uid IS NOT NULL").fetchone()[0]
    print(f"references: {nref} PDFs | paper_chunks: {nch} | linked to SE: {nlink}")
    print("\nsample BM25 search 'cortical tension actin':")
    rows = con.execute("""
        SELECT citation_key, page, substr(text,1,90) AS snippet,
               fts_main_paper_chunks.match_bm25(chunk_id, 'cortical tension actin') AS score
        FROM paper_chunks
        WHERE score IS NOT NULL ORDER BY score DESC LIMIT 5
    """).fetchall()
    for ck, pg, sn, sc in rows:
        print(f"  [{sc:.2f}] {ck} p{pg}: {sn.strip()}")
    con.close()


if __name__ == "__main__":
    if "--check" in sys.argv:
        sanity_summary()
    elif "--stats" in sys.argv:
        stats()
    else:
        ingest(allow_link_drop="--allow-link-drop" in sys.argv)
        print()
        sanity_summary()
        print()
        stats()
