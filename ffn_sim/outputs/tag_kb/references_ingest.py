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
REF_DIR = pathlib.Path("/Users/sw1/ffn_cellsim/ffn_sim/references")
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
    bundle = REF_DIR / "cellpress_bundle"
    if bundle.exists():
        paths += sorted(bundle.glob("*.pdf"))
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


def ingest():
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
    con.execute("DROP TABLE IF EXISTS paper_chunks")
    con.execute("DROP TABLE IF EXISTS paper_refs")
    con.execute(
        "CREATE TABLE paper_refs (citation_key TEXT, title TEXT, doi TEXT, "
        "path TEXT, n_pages INTEGER, sha1 TEXT, source TEXT, se_uid TEXT, "
        "se_citation_key TEXT, n_chunks INTEGER)")
    con.execute(
        "CREATE TABLE paper_chunks (chunk_id TEXT, citation_key TEXT, "
        "se_uid TEXT, page INTEGER, text TEXT)")

    seen_sha: dict[str, str] = {}
    ref_rows, chunk_rows, corpus = [], [], []
    skipped_dups = []

    for p in pdf_paths():
        sha = sha1(p)
        if sha in seen_sha:
            skipped_dups.append((str(p.name), seen_sha[sha]))
            continue
        seen_sha[sha] = p.name
        rel = p.name if p.parent == REF_DIR else f"cellpress_bundle/{p.name}"
        source = "new" if p.parent == REF_DIR else "bundle"
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
                chunk_rows.append((f"{ck}:p{pi}:{ci}", ck, se_uid or None, pi, ch))
                nc += 1
        ref_rows.append((ck, title, doi or None, rel, n_pages, sha, source,
                         se_uid or None, se_ck or None, nc))
        corpus.append({"citation_key": ck, "title": title, "doi": doi,
                       "path": rel, "n_pages": n_pages, "source": source,
                       "se_uid": se_uid, "se_citation_key": se_ck, "n_chunks": nc,
                       "linked_to_source_evidence": bool(se_uid)})

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
    if "--stats" in sys.argv:
        stats()
    else:
        ingest()
        print()
        stats()
